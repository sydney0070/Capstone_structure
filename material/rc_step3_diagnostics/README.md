# RC Step 3 — Load reversal diagnostic v2

## Purpose
Reproduce 24×12, 36×18, 48×24 fiber discretizations with Lobatto integration 5 and 7. For every successful substep, store tip displacement, load factor, base horizontal reaction, force imbalance, base moment, and beam local force. A case is **not** verified by convergence alone.

## Run on Mac
```bash
cd ~/Desktop/capstone/Capstone_structure/material/rc_step3_diagnostics
python diagnose_rc.py
```
For recovery experiments only:
```bash
python diagnose_rc.py --adaptive
```

Inspect `results_rc_diagnostics/diagnostic_summary.json` and each case CSV. Keep previous folders intact. If you run `--adaptive`, its outputs overwrite those of the default run; rename the results folder first if you need both.

## Caveats
- Uses the prior provisional C30/SD400 parameters and unchanged rectangle + rebar fiber section; this is deliberate to reproduce the discrepancy before changing model assumptions.
- The `Concrete02` patch includes the full section area; steel fiber area is superimposed. Fix overlapping material area as a later controlled comparison, not while diagnosing the existing discrepancy.
- This is 2D cantilever without imposed axial load, not the separate moment-curvature test and not a 3D RC dome shell.
- No residual displacement unloading here: focus is load reversal and force balance first.
- Successful convergence does not certify equilibrium. `NormDispIncr` alone is not an equilibrium check.
- If the reaction and load-factor outputs disagree, check whether force recovery, time-series/loading state, or output sign conventions are consistent before concluding loss of physical equilibrium.
