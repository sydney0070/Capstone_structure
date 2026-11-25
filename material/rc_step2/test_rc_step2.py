"""Step 2: independent monotonic and cyclic *material-level* verification.

Execute: python test_rc_step2.py
Dependencies: openseespy, numpy, matplotlib
No geometry/element is created. Units: Pa, dimensionless strain.
All non-strength parameters are provisional and NOT experimentally calibrated.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import openseespy.opensees as ops
from rc_materials import ConcreteParams, SteelParams, register_rc_materials, write_parameters_json

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_rc_step2'
CONCRETE = ConcreteParams()
STEEL = SteelParams()


def path_from_knots(knots, per_segment=100):
    """Continuous path, including each end knot once. per_segment >= 2."""
    if per_segment < 2:
        raise ValueError('per_segment must be >= 2')
    parts = [np.linspace(a, b, per_segment, endpoint=False)
             for a, b in zip(knots[:-1], knots[1:])]
    return np.concatenate(parts + [np.asarray([knots[-1]], dtype=float)])


# Tests with fresh uniaxial material states, one test per ops.wipe().
TESTS = {
    'concrete_monotonic_compression': ('concrete', [0.0, -0.0001, -0.002, -0.006]),
    'concrete_monotonic_tension': ('concrete', [0.0, 0.00005, 0.0003, 0.0025]),
    'concrete_cyclic': ('concrete', [0.0, -0.001, 0.0, -0.002, 0.0,
                                     -0.0035, 0.0, -0.005, 0.0]),
    'steel_monotonic_tension': ('steel', [0.0, 0.0005, 0.002, 0.006, 0.02]),
    'steel_monotonic_compression': ('steel', [0.0, -0.0005, -0.002, -0.006, -0.02]),
    'steel_cyclic': ('steel', [0.0, 0.005, -0.005, 0.005, -0.005,
                               0.010, -0.010, 0.010, -0.010, 0.0]),
}


def run_material(material: str, history: np.ndarray) -> np.ndarray:
    """Returns columns: strain, stress(Pa), tangent(Pa)."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    register_rc_materials(ops, concrete=CONCRETE, steel=STEEL)
    tag = {'concrete': 1, 'steel': 2}[material]
    ops.testUniaxialMaterial(tag)
    rows = []
    for eps in history:
        ops.setStrain(float(eps))
        stress, tangent = float(ops.getStress()), float(ops.getTangent())
        if not (np.isfinite(stress) and np.isfinite(tangent)):
            raise RuntimeError(f'Non-finite response for {material} at eps={eps}')
        rows.append((float(eps), stress, tangent))
    return np.asarray(rows, dtype=float)


def zero_stress_crossings(data):
    """Interpolated zero-stress strain on unloading steps (not necessarily residual strain)."""
    eps, sig = data[:, 0], data[:, 1]
    crossings = []
    for i in range(1, len(eps)):
        if i < 2 or eps[i] == eps[i - 1]:
            continue
        # Unloading: strain motion goes toward zero stress at the preceding state.
        if (sig[i] == 0 or sig[i - 1] * sig[i] < 0) and sig[i - 1] * (eps[i] - eps[i - 1]) < 0:
            if sig[i] == sig[i - 1]:
                continue
            root = eps[i - 1] - sig[i - 1] * (eps[i] - eps[i - 1]) / (sig[i] - sig[i - 1])
            # Ignore the zero initial state; keep nonzero crossings as a proxy.
            if abs(root) > 1e-10:
                crossings.append(float(root))
    return crossings


def cycle_work(data, history_knots, n_per_segment):
    """Net ∫sigma*dε on completed closed-strain cycles; NOT always exact dissipation.

    Concrete cyclic paths return to 0 after each compression pulse.
    Steel cyclic test has multiple ± alternating segments, not a truly periodic state.
    Report work only where the endpoints share the same strain.
    """
    eps, sig = data[:, 0], data[:, 1]
    records = []
    for j in range(len(history_knots) - 1):
        # Find the last return to the starting knot of a candidate multi-leg cycle.
        start = j * n_per_segment
        if start >= len(eps):
            continue
        for k in range(j + 2, len(history_knots)):
            if not np.isclose(history_knots[j], history_knots[k], atol=1e-12):
                continue
            end = k * n_per_segment
            if end >= len(eps):
                continue
            # Keep shortest adjacent closed excursions, avoid nested duplicates.
            if any(np.isclose(history_knots[j], history_knots[q], atol=1e-12)
                   for q in range(j + 1, k)):
                break
            work = float(np.trapezoid(sig[start:end + 1], eps[start:end + 1])) if hasattr(np, 'trapezoid') else float(np.trapz(sig[start:end + 1], eps[start:end + 1]))
            records.append({'from_knot': j, 'to_knot': k,
                            'net_work_J_per_m3': work,
                            'note': 'Net path work; not exact dissipated energy if internal state changes'})
            break
    return records


def summarize(data, material, name, knots, n_per_segment):
    eps, sig, tangent = data.T
    expected_E = CONCRETE.initial_tangent if material == 'concrete' else STEEL.E0
    near_idx = next((i for i in range(1, len(eps)) if abs(eps[i]) >= 1e-8), 1)
    initial_secant = float(sig[near_idx] / eps[near_idx]) if eps[near_idx] != 0 else None
    summary = {
        'test': name, 'material': material, 'n_points': int(len(eps)),
        'min_strain': float(np.min(eps)), 'max_strain': float(np.max(eps)),
        'min_stress_MPa': float(np.min(sig) / 1e6),
        'max_stress_MPa': float(np.max(sig) / 1e6),
        'initial_secant_GPa': None if initial_secant is None else initial_secant / 1e9,
        'reference_initial_tangent_GPa': float(expected_E / 1e9),
        'zero_stress_unloading_strains': zero_stress_crossings(data),
        'strain_path_work_J_per_m3': float(np.trapezoid(sig, eps)) if hasattr(np, 'trapezoid') else float(np.trapz(sig, eps)),
        'closed_strain_excursion_work': cycle_work(data, knots, n_per_segment) if 'cyclic' in name else [],
        'status': 'NUMERICAL_MATERIAL_RESPONSE_ONLY_NOT_EXPERIMENTAL_VALIDATION',
    }
    if name == 'concrete_monotonic_compression':
        peak_idx = int(np.argmin(sig))
        summary['compressive_peak_strain'] = float(eps[peak_idx])
        summary['compressive_peak_stress_MPa'] = float(sig[peak_idx] / 1e6)
        summary['expected_peak_MPa'] = CONCRETE.fpc / 1e6
    if name == 'concrete_monotonic_tension':
        summary['expected_tensile_strength_MPa'] = CONCRETE.ft / 1e6
    if name == 'steel_monotonic_tension':
        summary['nominal_yield_strain'] = STEEL.Fy / STEEL.E0
        summary['nominal_yield_stress_MPa'] = STEEL.Fy / 1e6
    return summary


def save_csv(path: Path, data: np.ndarray):
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['strain', 'stress_Pa', 'tangent_Pa'])
        writer.writerows(data)


def plot_series(path: Path, data: np.ndarray, title: str):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(data[:, 0], data[:, 1] / 1e6, lw=1.6)
    ax.set(xlabel='Strain (-)', ylabel='Stress (MPa)', title=title)
    ax.axhline(0, color='gray', lw=0.7)
    ax.axvline(0, color='gray', lw=0.7)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    write_parameters_json(OUT / 'material_parameters.json', CONCRETE, STEEL)
    report = {'metadata': {'units': 'Pa, strain (-), J/m^3',
                           'experimentally_calibrated': False,
                           'warning': 'Materials are uniaxial prototypes. Not an RC shell/section.'},
              'tests': {}}
    for name, (material, knots) in TESTS.items():
        knots = list(knots)
        n = 100
        history = path_from_knots(knots, n)
        data = run_material(material, history)
        save_csv(OUT / f'{name}.csv', data)
        plot_series(OUT / f'{name}.png', data, name)
        metrics = summarize(data, material, name, knots, n)
        # A finer strain increment is a numerical resolution check, not experimental validation.
        fine = run_material(material, path_from_knots(knots, 2*n))
        # samples are aligned at identical strain knots and segment midsteps
        coarse_at_knots = data[::n, 1]
        fine_at_knots = fine[::2*n, 1]
        scale = max(1e6, float(np.max(np.abs(fine_at_knots))))
        metrics['step_refinement_max_knot_difference_MPa'] = float(np.max(np.abs(coarse_at_knots - fine_at_knots)) / 1e6)
        metrics['step_refinement_relative_knot_difference'] = float(np.max(np.abs(coarse_at_knots - fine_at_knots)) / scale)
        report['tests'][name] = metrics
        print(f'{name}: {len(data)} states; max stress={metrics["max_stress_MPa"]:.3f} MPa; '
              f'min stress={metrics["min_stress_MPa"]:.3f} MPa; '
              f'refinement Δ={metrics["step_refinement_max_knot_difference_MPa"]:.5f} MPa')
    (OUT / 'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'Saved: {OUT}')
    print('NOTE: Work integrals are path work, not automatically verified dissipated energy.')
    print('NOTE: Results do not establish physical RC performance or experimental calibration.')


if __name__ == '__main__':
    main()
