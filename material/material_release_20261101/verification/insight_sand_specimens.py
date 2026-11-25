"""Build the three InSight sandy-regolith trial specimens and run triaxial tests."""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np
import openseespy.opensees as ops

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from jez1_density import check_density, triaxial

PROJECT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT / 'ground' / 'ground_modeling'))
from ground_hybrid.materials import record


CANDIDATES = {
    'mgs1': ('MGS-1', 'mars_insight_sand_mgs1_dp_1300'),
    'mgs1s': ('MGS-1S', 'mars_insight_sand_mgs1s_dp_1300'),
    'jez1': ('JEZ-1', 'mars_insight_sand_jez1_dp_1300'),
}
# Predicted confinement at 1 m and 2 m in Murdoch et al., Table 2.
CONFINEMENTS_PA = (5440.0, 10300.0)


def run_specimen(candidate, output):
    label, material_id = CANDIDATES[candidate]
    parameters = record(material_id)
    density = parameters['density_kg_m3']
    result = {
        'label': label,
        'material_id': material_id,
        'material': parameters,
        'opensees_version': ops.version(),
        'element_type': 'bbarBrick',
        'loading': 'Constant lateral confinement, then displacement-controlled axial compression without selfweight',
        'scope': 'Numerical specimen using transferred initial parameters',
        'tests': [],
    }
    curves = {}
    with (output / f'{candidate}_curves.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['confinement_kPa', 'axial_compression_strain', 'q_kPa'])
        for pressure in CONFINEMENTS_PA:
            curve, metrics = triaxial(parameters, density, pressure)
            strain, q = curve.T
            confinement_kpa = pressure / 1000
            curves[confinement_kpa] = (100 * strain, q / 1000)
            writer.writerows((confinement_kpa, x, y / 1000) for x, y in curve)
            if not math.isclose(metrics['initial_axial_modulus_MPa'], parameters['E_Pa'] / 1e6, rel_tol=1e-5):
                raise RuntimeError(f'{material_id}: unexpected initial axial stiffness')
            result['tests'].append({
                'confinement_kPa': confinement_kpa,
                'maximum_axial_strain': float(strain[-1]),
                **metrics,
            })

    node_tags = list(ops.eleNodes(1))
    coordinates = [list(ops.nodeCoord(tag)) for tag in node_tags]
    result['node_coordinates_m'] = dict(zip(node_tags, coordinates))
    result['element_connectivity'] = node_tags
    result['dimensions_m'] = np.ptp(np.asarray(coordinates), axis=0).tolist()
    result['mass_from_eigenvalue_kg'] = check_density(parameters, density)
    (output / f'{candidate}_specimen.json').write_text(
        json.dumps(result, indent=2, allow_nan=False) + '\n'
    )
    return curves


def plot_comparison(curves_by_candidate, path):
    colors = {'mgs1': '#226A83', 'mgs1s': '#AC5127', 'jez1': '#547437'}
    figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.9), sharey=True, layout='constrained')
    for ax, pressure in zip(axes, CONFINEMENTS_PA):
        confinement_kpa = pressure / 1000
        for candidate, curves in curves_by_candidate.items():
            strain, q = curves[confinement_kpa]
            ax.plot(strain, q, label=CANDIDATES[candidate][0], color=colors[candidate], linewidth=2)
        ax.set_title(f'Confining pressure {confinement_kpa:g} kPa')
        ax.set_xlabel('Axial compression strain (%)')
        ax.set_xlim(0, 0.2)
        ax.set_xticks([0.0, 0.05, 0.10, 0.15, 0.20])
        ax.grid(alpha=0.2)
    peak = max(float(q.max()) for curves in curves_by_candidate.values() for _, q in curves.values())
    axes[0].set_ylim(0, 1.12 * peak)
    axes[0].set_ylabel('Deviator stress q (kPa)')
    axes[1].legend(loc='lower right', frameon=False)
    figure.savefig(path, dpi=185, facecolor='white')
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', choices=CANDIDATES, help='run one specimen instead of all three')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error('choose an empty output directory')
    args.output_dir.mkdir(parents=True, exist_ok=True)

    candidates = (args.candidate,) if args.candidate else CANDIDATES
    try:
        curves = {candidate: run_specimen(candidate, args.output_dir) for candidate in candidates}
        plot_comparison(curves, args.output_dir / 'triaxial_comparison.png')
    finally:
        ops.wipe()


if __name__ == '__main__':
    main()
