"""Analytical-energy and verdict tests; no OpenSees dependency."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from material_evaluation import cycle_metrics, evaluate, load_targets, json_safe


def specimen(damping=0.05):
    q, n, amp, G = 200, 3, 0.001, 10000.0
    phase = np.linspace(0, 2 * np.pi * n, 4 * q * n + 1)
    gamma = amp * np.sin(phase)
    # Ellipse: integral tau dgamma = 2*pi*D*G*amp**2.
    tau = G * amp * (np.sin(phase) + 2 * damping * np.cos(phase))
    stress = np.zeros((len(gamma), 6))
    stress[:, :3] = [-45, -45, -100]
    stress[:, 5] = tau
    strain = np.zeros_like(stress)
    strain[:, 5] = gamma
    displacement = np.zeros((len(gamma), 8, 3))
    displacement[:, 4:, 0] = gamma[:, None]
    reaction = np.zeros((len(gamma), 3))
    reaction[:, 0], reaction[:, 2] = tau, -100
    return dict(gamma=gamma, tau=tau, stress=stress, strain=strain,
                displacement=displacement, reaction=reaction,
                eta=np.full(len(gamma), 0.5), target=gamma.copy(),
                steps=q, cycles=n, amplitude=amp, sigma_v=100.0,
                height=1.0, rho=1.9, completed=True, error=None,
                material={}, opensees_version="synthetic")


class EvaluationTests(unittest.TestCase):
    def test_ellipse_energy_and_reversal_modulus(self):
        cycles = cycle_metrics(specimen())
        self.assertIsNone(cycles[0]["D_percent"])
        self.assertAlmostEqual(cycles[-1]["D_percent"], 5.0, places=3)
        self.assertAlmostEqual(cycles[-1]["Gsec_kPa"], 10000.0)

    def test_elastic_loop_zero_damping(self):
        self.assertAlmostEqual(cycle_metrics(specimen(0))[-1]["D_percent"], 0, places=10)

    def test_negative_energy_not_hidden(self):
        report = evaluate(specimen(-0.05), {}, {})
        self.assertEqual(report["consistency_status"], "FAIL")
        self.assertLess(report["metrics"]["D_last_percent"], 0)

    def test_volume_drift_withholds_damping(self):
        run = specimen()
        run["strain"][:, 2] = np.linspace(0, -0.003, len(run["gamma"]))
        self.assertIsNone(cycle_metrics(run)[-1]["D_percent"])

    def test_open_stress_loop_withholds_damping(self):
        run = specimen()
        run["tau"] += np.linspace(0, 10, len(run["gamma"]))
        self.assertIsNone(cycle_metrics(run)[-1]["D_percent"])

    def test_missing_targets_not_pass(self):
        report = evaluate(specimen(), {}, {})
        self.assertEqual(report["target_status"], "NOT_EVALUATED")
        self.assertEqual(report["numerical_status"], "NOT_EVALUATED")

    def target(self, lower, upper):
        return {"context": {"conditions_confirmed": True, "soil_description": "synthetic test only",
                            "test": "drained_constant_sigma_v_simple_shear", "cycles": 3,
                            "gamma_amplitude": 0.001, "sigma_v_kPa": 100.0},
                "bounds": {"D_last_percent": {"lower": lower, "upper": upper, "source": "analytical ellipse"}}}

    def test_target_pass_and_fail(self):
        for low, high, status in [(4.9, 5.1, "PASS"), (1, 2, "FAIL")]:
            report = evaluate(specimen(), {}, self.target(low, high))
            row = next(r for r in report["target_checks"] if r["name"] == "D_last_percent")
            self.assertEqual(row["status"], status)

    def test_wrong_conditions_not_compared(self):
        target = self.target(1, 2)
        target["context"]["sigma_v_kPa"] = 200
        self.assertEqual(evaluate(specimen(), {}, target)["target_status"], "NOT_EVALUATED")

    def test_incomplete_run_fails_numerical(self):
        run = specimen()
        run["completed"] = False
        report = evaluate(run, {}, self.target(4, 6))
        self.assertEqual(report["numerical_status"], "FAIL")
        self.assertEqual(report["target_status"], "NOT_EVALUATED")

    def test_eta_violation_fails_consistency(self):
        run = specimen()
        run["eta"][-1] = 1.2
        self.assertEqual(evaluate(run, {}, {})["consistency_status"], "FAIL")

    def test_invalid_target_rejected(self):
        target = self.target(6, 4)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "target.json"
            path.write_text(json.dumps(target))
            with self.assertRaises(ValueError):
                load_targets(path)

    def test_nonfinite_json_becomes_null(self):
        self.assertEqual(json.dumps(json_safe({"bad": float("nan")}), allow_nan=False), '{"bad": null}')

    def test_truncated_cycle_not_counted(self):
        run = specimen()
        run["gamma"] = run["gamma"][:-1]
        self.assertEqual(len(cycle_metrics(run)), 2)


if __name__ == "__main__":
    unittest.main()
