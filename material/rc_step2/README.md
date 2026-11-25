# RC material verification — Step 2

## Files
- `rc_materials.py` — unchanged original Concrete02/Steel02 definitions; same provisional parameters as Step 1.
- `test_rc_step2.py` — six **independent fresh-state** material tests: concrete monotonic compression, concrete monotonic tension, concrete cyclic; steel monotonic tension, steel monotonic compression, steel cyclic.
- `results_rc_step2/` — created locally with CSV, PNG, JSON.

## Run on macOS (from directory containing both .py files)

```bash
python -m pip install numpy matplotlib openseespy
python test_rc_step2.py
```

Or copy both `.py` files to the existing `material` folder, replacing `rc_materials.py` only if it is the same original version.

## Source and status
User-selected: concrete nominal compressive strength 30 MPa and SD400 nominal yield stress 400 MPa.
All other constitutive parameters are **provisional**, not experimental calibration. Details and sources are in `material_parameters.json` generated on run.

- Concrete02: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Concrete02.html
- Steel02: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Steel02.html
- Uniaxial material testing: https://opensees.github.io/OpenSeesDocumentation/user/manual/misc/testUniaxialMaterial.html

## Interpretation
- Monotonic tension and compression are separated to avoid carrying damage history into the other loading direction.
- Cyclic histories include load reversal and unloading, for hysteretic response.
- `zero_stress_unloading_strains` finds approximate zero-stress crossings along unloading segments; it is a *proxy* only, not a verified permanent strain after stress-controlled unloading.
- `strain_path_work_J_per_m3` is numerical ∫σ dε over the specified path; **do not automatically call it energy dissipated** when the path is nonclosed in material state or strain is not periodic. Some closed-strain excursions are also reported as net path work, with the same caveat.
- Step refinement compares the stress at identical strain knots using 100 versus 200 subdivisions per segment. These are not experimental validation results.
- Direct `testUniaxialMaterial` testing verifies material response; `RC` composite action, confinement, anchorage, plate/shell multiaxial stress states and structural failure are **not tested**.
