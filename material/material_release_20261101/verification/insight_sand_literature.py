"""Audit the initial InSight sand inputs against their PDFs and simple-shear response."""

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import openseespy.opensees as ops

from jez1_density import COORDINATES, NODES, static_analysis

PROJECT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT / 'ground' / 'ground_modeling'))
from ground_hybrid.materials import record, register
from ground_dome_footing import SANDY_REGOLITH_CANDIDATES

PRESSURES_PA = (97.0, 250.0, 500.0, 680.0)


def plot_comparison(summary, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    colors = {'mgs1': '#226A83', 'mgs1s': '#AC5127', 'jez1': '#547437'}
    figure, axes = plt.subplots(1, 3, figsize=(10.5, 3.7), sharex=True, sharey=True, layout='constrained')
    for ax, candidate in zip(axes, SANDY_REGOLITH_CANDIDATES):
        rows = [row for row in summary['cases'] if row['candidate'] == candidate]
        x = [row['normal_stress_kPa'] for row in rows]
        ax.plot(x, [row['dp_plateau_tau_kPa'] for row in rows], '-o', color=colors[candidate], label='OpenSees DP')
        ax.plot(x, [row['graph_envelope_tau_kPa'] for row in rows], '--s', color='#52616B', label='Paper envelope')
        ax.set_title({'mgs1': 'MGS-1', 'mgs1s': 'MGS-1S', 'jez1': 'JEZ-1'}[candidate])
        ax.set_xlabel('Normal stress (kPa)')
        ax.set_xlim(0, 0.75)
        ax.set_ylim(0, 1.0)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('Shear resistance (kPa)')
    axes[-1].legend(loc='lower right', frameon=False)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def elastic_constants(density, vs, vp):
    shear = density * vs**2
    bulk = density * (vp**2 - 4 * vs**2 / 3)
    return dict(
        E_Pa=9 * bulk * shear / (3 * bulk + shear),
        nu=(3 * bulk - 2 * shear) / (2 * (3 * bulk + shear)),
        G_Pa=shear, K_Pa=bulk,
    )


def simple_shear(material_id, pressure):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    parameters = register(ops, material_id, 1)
    for tag, point in zip(NODES, COORDINATES):
        ops.node(tag, *point)
        ops.fix(tag, *((1, 1, 1) if tag <= 4 else (0, 1, 0)))
    for tag in (6, 7, 8):
        ops.equalDOF(5, tag, 1, 3)
    ops.element('bbarBrick', 1, *NODES, 1, 0.0, 0.0, 0.0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(5, 0.0, 0.0, -pressure)
    static_analysis(0.1)
    ops.test('NormUnbalance', 1e-6, 80)
    if ops.analyze(10):
        raise RuntimeError('Normal load initialization failed')
    ops.loadConst('-time', 0.0)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    ops.load(5, 1.0, 0.0, 0.0)

    rows = []
    targets = np.r_[np.linspace(0, 1e-5, 11), np.linspace(2e-5, 0.002, 100)]
    for target in targets:
        if target:
            ops.integrator('DisplacementControl', 5, 1, float(target - ops.nodeDisp(5, 1)))
            if ops.analyze(1):
                raise RuntimeError('Shear displacement increment failed')
        stress = np.asarray([ops.eleResponse(1, 'material', gp, 'stress') for gp in range(1, 9)])
        if stress.shape != (8, 6) or not np.isfinite(stress).all():
            raise RuntimeError('Missing integration-point response')
        if not np.allclose(stress[:, 2], -pressure, rtol=1e-5, atol=1e-3):
            raise RuntimeError('Vertical pressure changed')
        ops.reactions()
        force = ops.getLoadFactor(2)
        reaction = sum(ops.nodeReaction(tag, 1) for tag in NODES[:4])
        if not math.isclose(force + reaction, 0, abs_tol=1e-3):
            raise RuntimeError('Shear reaction equilibrium failed')
        tau = float(stress[:, 5].mean())
        if not math.isclose(tau, force, rel_tol=1e-5, abs_tol=1e-3):
            raise RuntimeError('Shear stress and applied traction disagree')
        rows.append((float(target), tau))

    curve = np.asarray(rows)
    stiffness = float(np.polyfit(curve[:5, 0], curve[:5, 1], 1)[0])
    if not math.isclose(stiffness, parameters['G_Pa'], rel_tol=1e-5):
        raise RuntimeError('Unexpected initial shear modulus')
    if not np.allclose(curve[-10:, 1], curve[-1, 1], rtol=1e-5, atol=1e-3):
        raise RuntimeError('Shear plateau was not reached')
    return curve, stiffness


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output-dir', type=Path, required=True)
    args = cli.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        cli.error('choose an empty output directory')
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdf_dir = Path('/Users/sydney/Downloads')
    pdf_names = ('2403.11029v1.pdf', 's11214-023-00964-0.pdf', 'InSightLanderFlexibleModes_Changes.pdf')
    sources = {
        name: dict(path=str(pdf_dir / name), sha256=hashlib.sha256((pdf_dir / name).read_bytes()).hexdigest())
        for name in pdf_names
    }
    summary = dict(
        sources=sources,
        opensees_version=ops.version(),
        method='Homogeneous plane-strain simple shear at constant vertical stress on a 1 m cube',
        comparison_scope='Approximate Dotson graph-derived Mohr–Coulomb envelopes compared with the implemented DP law',
        normal_stress_range_Pa=list(PRESSURES_PA),
        dotson_graph_readings=dict(
            mgs1=dict(density_kg_m3=1320, phi_deg=35.5, c_Pa=280),
            mgs1s=dict(density_kg_m3=1310, phi_deg=34.0, c_Pa=280),
            jez1=dict(density_kg_m3=1380, phi_deg=24.0, c_Pa=230),
        ),
        input_density_kg_m3=1300,
        murdoch_table2_1m=dict(density_kg_m3=1307, Vs_m_s=94.6, Vp_m_s=157.9, confinement_Pa=5440),
        robotic_arm_section_7_7=dict(
            pdf_page=27, sand_density_range_kg_m3=[1300, 1500],
            bulk_top_40cm_density_kg_m3=1211, bulk_top_40cm_E_Pa=11.48e6,
            bulk_top_40cm_G_Pa=4.47e6, bulk_top_40cm_K_Pa=7.79e6,
            bulk_top_40cm_nu=0.29,
        ),
        cases=[],
    )
    original = elastic_constants(1307, 94.6, 157.9)
    target = elastic_constants(1300, 94.6, 157.9)
    summary['murdoch_1m_derived_moduli'] = original
    summary['murdoch_1m_moduli_at_target_density'] = target
    summary['shared_input_moduli'] = {key: record(SANDY_REGOLITH_CANDIDATES['mgs1'])[key] for key in target}
    summary['relative_rounding_errors'] = {
        key: abs(summary['shared_input_moduli'][key] / value - 1) for key, value in target.items()
    }
    try:
        with (args.output_dir / 'shear_comparison.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['candidate', 'normal_stress_kPa', 'shear_strain', 'dp_tau_kPa', 'graph_envelope_tau_kPa'])
            for candidate, material_id in SANDY_REGOLITH_CANDIDATES.items():
                params = record(material_id)
                cases = []
                for pressure in PRESSURES_PA:
                    curve, stiffness = simple_shear(material_id, pressure)
                    mc_tau = params['c_Pa'] + pressure * math.tan(math.radians(params['phi_deg']))
                    writer.writerows((candidate, pressure / 1000, strain, tau / 1000, mc_tau / 1000) for strain, tau in curve)
                    case = dict(
                        candidate=candidate, normal_stress_kPa=pressure / 1000,
                        dp_plateau_tau_kPa=float(curve[-1, 1] / 1000),
                        graph_envelope_tau_kPa=mc_tau / 1000,
                        relative_strength_difference=float(curve[-1, 1] / mc_tau - 1),
                        initial_G_Pa=stiffness,
                    )
                    cases.append(case)
                    summary['cases'].append(case)
                regression = np.polyfit(PRESSURES_PA, [row['dp_plateau_tau_kPa'] * 1000 for row in cases], 1)
                summary.setdefault('effective_simple_shear_parameters', {})[candidate] = dict(
                    phi_deg=math.degrees(math.atan(regression[0])), c_Pa=float(regression[1]),
                )
        (args.output_dir / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
        plot_comparison(summary, args.output_dir / 'shear_comparison.png')
    finally:
        ops.wipe()


if __name__ == '__main__':
    main()
