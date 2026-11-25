# RC forceBeamColumn vs dispBeamColumn diagnostic comparison

Purpose: isolate element internal compatibility issues, numerical force imbalance, and fiber subdivision effects. This is not a calibration or physical validation.

## Run on macOS
```
cd ~/Desktop/capstone/Capstone_structure/material/rc_element_comparison
python compare_elements.py
```
Optional recovery run (do not overwrite baseline output before archiving):
```
python compare_elements.py --retry
```

## What's preserved
- Original user-uploaded `rc_section.py` (unmodified, including concrete/steel overlap)
- Original provisional Concrete02/Steel02 parameters
- Same imposed displacement history and PDelta transform
- `diagnose_rc_original.py` copy for reference

## What's new
- `forceBeamColumn -iter 30 1e-10` and 50 iterations in one comparison
- `dispBeamColumn` alternative formulation
- Same 24x12, 36x18, 48x24 subdivision choices
- Equilibrium gate stops any run at relative imbalance above 1e-3 (this threshold is a diagnostic hard stop, not a passing criterion)
- PASS requires COMPLETE and max relative imbalance below 1e-4
- Per-step CSV including applied lateral load, reaction, deformation and step location

## Interpret with care
- Large errors after load reversal may be numerical. More iterations do not prove physical validity.
- Displacement-based vs force-based elements need not give identical answers at the same discretization.
- The concrete rectangular patch overlaps steel reinforcement volume. This is held unchanged for diagnostic isolation.
- Do not classify non-converged results as failure of RC material; numerical convergence and force equilibrium are distinct checks.
- The 10-MN force guard prevents extreme unphysical values from propagating during diagnostics.
- This code has been syntax checked only. Run with your OpenSeesPy environment and share summary.json + CSV files.
