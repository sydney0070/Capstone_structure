"""Compare two JEZ-1 residual densities in OpenSees with identical DP parameters."""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import openseespy.opensees as ops

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from capmat.materials import dp_constants, material


COORDINATES = (
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
)
NODES = tuple(range(1, 9))
CONFINEMENTS_PA = (5_000, 15_000, 25_000, 50_000, 75_000)


def start_specimen(parameters, density):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    dp = dp_constants(parameters['phi_deg'], parameters['c_Pa'], parameters['psi_deg'])
    ops.nDMaterial(
        'DruckerPrager', 1, parameters['K_Pa'], parameters['G_Pa'],
        dp['sigmaY'], dp['rho'], dp['rhoBar'], 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, density, 101000.0,
    )
    for tag, xyz in zip(NODES, COORDINATES):
        ops.node(tag, *xyz)


def static_analysis(increment):
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormUnbalance', 1e-4, 100, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', increment)
    ops.analysis('Static')


def analyze(steps=1):
    if ops.analyze(steps):
        raise RuntimeError(f'OpenSees failed at pseudo-time {ops.getTime()}')


def triaxial(parameters, density, confinement):
    start_specimen(parameters, density)
    for tag, (x, y, z) in zip(NODES, COORDINATES):
        ops.fix(tag, int(x == 0), int(y == 0), int(z == 0))
    ops.element('bbarBrick', 1, *NODES, 1, 0.0, 0.0, 0.0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for tag, (x, y, z) in zip(NODES, COORDINATES):
        ops.load(tag, -confinement * x / 4, -confinement * y / 4, -confinement * z / 4)
    static_analysis(0.1)
    analyze(10)
    ops.loadConst('-time', 0.0)

    axial_reference = ops.nodeDisp(5, 3)
    axial_targets = axial_reference - np.linspace(0, 0.03, 601)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *axial_targets.tolist(), '-useLast')
    ops.pattern('Plain', 2, 2)
    for tag in NODES[4:]:
        ops.sp(tag, 3, 1.0)
    ops.wipeAnalysis()
    static_analysis(0.0)
    analyze()
    ops.integrator('LoadControl', 1.0)

    rows = []
    for step in range(len(axial_targets)):
        if step:
            analyze()
        stresses = np.asarray([
            ops.eleResponse(1, 'material', gp, 'stress') for gp in range(1, 9)
        ])
        if stresses.shape != (8, 6) or not np.isfinite(stresses).all():
            raise RuntimeError('Incomplete stress response at the eight integration points')
        compression = -stresses[:, :3]
        if not np.allclose(compression[:, :2], confinement, rtol=1e-5, atol=1e-3):
            raise RuntimeError('Lateral confinement did not remain constant')
        if not np.allclose(compression, compression[0], rtol=1e-7, atol=1e-3):
            raise RuntimeError('Nonuniform stress in the single-element specimen')
        q = float((compression[:, 2] - compression[:, :2].mean(axis=1)).mean())
        rows.append((axial_reference - ops.nodeDisp(5, 3), q))

    curve = np.asarray(rows)
    strain, q = curve.T
    if abs(q[0]) > 1e-3 or not np.allclose(strain, np.linspace(0, 0.03, 601), atol=1e-9):
        raise RuntimeError('Invalid hydrostatic reference or imposed displacement history')
    plateau = float(q[-1])
    if not np.allclose(q[-20:], plateau, rtol=1e-5, atol=1e-3):
        raise RuntimeError('No plastic plateau reached within the applied strain range')
    return curve, {
        'initial_axial_modulus_MPa': float(np.polyfit(strain[:5], q[:5], 1)[0] / 1e6),
        'plateau_q_kPa': plateau / 1000,
    }


def check_density(parameters, density):
    # A free rigid translation of a 1 m³ brick against a known spring measures its mass.
    start_specimen(parameters, density)
    for tag in NODES:
        ops.fix(tag, 0, 1, 1)
    for tag in NODES[1:]:
        ops.equalDOF(1, tag, 1)
    ops.element('bbarBrick', 1, *NODES, 1, 0.0, 0.0, 0.0)
    ops.node(9, 0, 0, 0)
    ops.fix(9, 1, 1, 1)
    spring = 1e6
    ops.uniaxialMaterial('Elastic', 2, spring)
    ops.element('zeroLength', 2, 9, 1, '-mat', 2, '-dir', 1)
    static_analysis(0.0)
    eigenvalue = ops.eigen('-fullGenLapack', 1)[0]
    recovered_mass = spring / eigenvalue
    if not math.isclose(recovered_mass, density, rel_tol=1e-8):
        raise RuntimeError(f'Incorrect mass for density {density}: {recovered_mass}')
    return recovered_mass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    output = parser.parse_args().output_dir
    if output.exists() and any(output.iterdir()):
        parser.error('choose an empty output directory')
    output.mkdir(parents=True, exist_ok=True)

    parameters = material('mars_jez1_residual')
    reference_density = parameters['density_kg_m3']
    comparison_density = 1300
    densities = (reference_density, comparison_density)
    sine = math.sin(math.radians(parameters['phi_deg']))
    compression_ratio = (1 + sine) / (1 - sine)
    summary = {
        'analysis': 'OpenSees single-element triaxial compression, no specimen selfweight',
        'scope': 'Density-only numerical sensitivity with shared constitutive parameters',
        'opensees_version': ops.version(),
        'material_id': 'mars_jez1_residual',
        'shared_parameters': {
            key: parameters[key] for key in ('E_Pa', 'nu', 'G_Pa', 'K_Pa', 'phi_deg', 'c_Pa', 'psi_deg')
        },
        'reference_density_kg_m3': reference_density,
        'comparison_density_kg_m3': comparison_density,
        'cases': [],
    }
    try:
        summary['mass_from_eigenvalue_kg_per_m3'] = {
            density: check_density(parameters, density) for density in densities
        }
        with (output / 'curves.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['density_kg_m3', 'confinement_kPa', 'axial_compression_strain', 'q_kPa'])
            for pressure in CONFINEMENTS_PA:
                curves, metrics = {}, {}
                for density in densities:
                    curve, result = triaxial(parameters, density, pressure)
                    curves[density] = curve
                    metrics[density] = result
                    writer.writerows((density, pressure / 1000, strain, q / 1000) for strain, q in curve)
                expected_q = (compression_ratio - 1) * pressure + 2 * parameters['c_Pa'] * math.sqrt(compression_ratio)
                for result in metrics.values():
                    if not math.isclose(result['initial_axial_modulus_MPa'], parameters['E_Pa'] / 1e6, rel_tol=1e-5):
                        raise RuntimeError('Unexpected axial stiffness')
                    if not math.isclose(result['plateau_q_kPa'], expected_q / 1000, rel_tol=1e-4):
                        raise RuntimeError('Computed strength differs from the analytical DP prediction')
                gap = float(np.max(np.abs(curves[reference_density][:, 1] - curves[comparison_density][:, 1])))
                summary['cases'].append({
                    'confinement_kPa': pressure / 1000,
                    'specimens': metrics,
                    'max_curve_difference_kPa': gap / 1000,
                })
        (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    finally:
        ops.wipe()


if __name__ == '__main__':
    main()
