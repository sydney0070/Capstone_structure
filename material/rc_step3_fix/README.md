# RC Step 3 diagnostic fix

This package retains the same provisional Concrete02 and Steel02 parameters and section geometry as Step 3.

## Changes
1. Define `LoadControl` integrator **before** creating Static analysis; removes the default-integrator warning.
2. Apply the intended axial load in moment–curvature with a load increment of 1.0 (the old `LoadControl(0.0)` applied no axial load).
3. Try `NormDispIncr` at 1e-9 with 80 iterations for numerical diagnosis. Passing this criterion does **not** establish force equilibrium; inspect reactions and sensitivity before accepting final results.
4. Reduce beam lateral displacement increment from 0.1 mm to 0.05 mm.
5. If beam analysis fails, save completed points to `beam_cyclic_partial.csv`.

## Run (macOS)
`cd ~/Desktop/capstone/Capstone_structure/material/rc_step3_fix`
`python test_rc_step3.py`

The moment–curvature data **must be regenerated**, because axial load application was corrected. Values will differ from the previous 320.94 kN·m. The beam case still has no axial preload (same as original). Numerical execution is not experimental material calibration.
