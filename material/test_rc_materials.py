"""Direct material tests (no structure model required).

Usage: python test_rc_materials.py
Dependencies: openseespy, numpy, matplotlib
Expected output: CSV files, PNG plots and parameter JSON in results_rc/.
"""
from pathlib import Path
import csv

import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
from rc_materials import ConcreteParams, SteelParams, register_rc_materials, write_parameters_json


OUT = Path(__file__).resolve().parent / 'results_rc'


def strain_history_concrete():
    # Compression negative; unload and tension after compression.
    knots = [0., -0.0001, -0.001, -0.002, -0.004, 0.,
             0.00005, 0.0002, 0., -0.002, -0.006]
    return np.concatenate([np.linspace(a, b, 40, endpoint=False)
                           for a, b in zip(knots[:-1], knots[1:])] + [np.array([knots[-1]])])


def strain_history_steel():
    knots = [0., 0.0005, 0.003, 0.01, 0., -0.005,
             0.005, -0.015, 0.015, 0.]
    return np.concatenate([np.linspace(a, b, 40, endpoint=False)
                           for a, b in zip(knots[:-1], knots[1:])] + [np.array([knots[-1]])])


def material_history(mat_tag, history):
    ops.testUniaxialMaterial(mat_tag)
    rows = []
    for eps in history:
        # setStrain automatically commits the state according to OpenSees documentation.
        ops.setStrain(float(eps))
        rows.append((float(eps), float(ops.getStress()), float(ops.getTangent())))
    return rows


def save_results(stem, rows):
    with (OUT / (stem + '.csv')).open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['strain', 'stress_Pa', 'tangent_Pa'])
        writer.writerows(rows)
    data = np.asarray(rows)
    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    ax.plot(data[:, 0], data[:, 1] / 1e6, lw=1.8)
    ax.set(xlabel='Strain (-)', ylabel='Stress (MPa)', title=stem)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / (stem + '.png'), dpi=170)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    concrete = ConcreteParams()
    steel = SteelParams()
    register_rc_materials(ops, concrete=concrete, steel=steel)
    write_parameters_json(OUT / 'parameters.json', concrete, steel)

    crows = material_history(1, strain_history_concrete())
    srows = material_history(2, strain_history_steel())
    save_results('Concrete02_C30_provisional', crows)
    save_results('Steel02_SD400_provisional', srows)

    # Small-strain initial secant checks use fresh material states.
    # Re-register before each check to avoid comparing an already damaged state.
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    register_rc_materials(ops, concrete=concrete, steel=steel)
    ops.testUniaxialMaterial(1)
    ops.setStrain(-1e-8)
    concrete_secant = ops.getStress() / -1e-8
    ops.testUniaxialMaterial(2)
    ops.setStrain(1e-8)
    steel_secant = ops.getStress() / 1e-8
    print('Concrete02 tangent estimate:', concrete_secant / 1e9, 'GPa; expected', concrete.initial_tangent / 1e9)
    print('Steel02 tangent estimate:', steel_secant / 1e9, 'GPa; nominal', steel.E0 / 1e9)
    print('Saved:', OUT)


if __name__ == '__main__':
    main()
