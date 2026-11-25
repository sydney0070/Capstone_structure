# RC Step 3 — Final numerical checks

## Purpose
1. Apply the validated nine-stage cyclic tip displacement history to the provisional C30/SD400 3 m RC beam.
2. After returning to tip displacement zero, unload the horizontal **load factor** to zero using 120 `LoadControl` steps and record remaining tip displacement (true zero-lateral-force residual displacement under this model).
3. Independently compare discretizations: fiber patches 24x12 vs 48x24; Lobatto integration points 5 vs 7; 2x2 matrix.
4. Save each run separately in `results_rc_final_validation/`, including partial data if analysis fails.

## How to run
```bash
cd ~/Desktop/capstone/Capstone_structure/material/rc_step3_final_validation
python validate_final.py
```
Requires `openseespy` and `numpy`. Copy this entire folder under your `material` directory.

## Outputs
- `results_rc_final_validation/final_summary.json`
- `results_rc_final_validation/<case>/history.csv`
- `partial_history.csv` for failed cases

## Interpretation
- `final_tip_displacement_mm` in baseline is residual displacement **only if** `final_force_N` is near zero and support force equilibrium passes.
- Compare force histories and peak forces across the four discretizations. No universal tolerance is assumed; decide acceptance thresholds before declaring convergence.
- The beam test has no vertical/axial preload, even though the separate moment-curvature test had -360 kN.
- All dimensions and constitutive parameters remain provisional. No experimental calibration, physical validation, RC Shell verification, or cracking/objectivity guarantee.
- A failed zero-load unloading may reflect difficulty with load control or a limit point: inspect `partial_history.csv` rather than calling it a physical failure.
