# RC Step 3 (SI units)

This extends Step 2: C30 Concrete02 + SD400 Steel02, **all other constitutive parameters provisional**.

Files:
- `rc_materials.py`: Step 2 materials, preserved.
- `rc_section.py`: assumed 400×600 mm RC rectangular section, top/bottom 3-D20 rebar, 45 mm cover to rebar outer surface. **No confinement enhancement, bond slip, shear failure, or damage calibration.**
- `test_rc_step3.py`: two independent OpenSeesPy tests.

Run from the `rc_step3` folder:

```bash
python test_rc_step3.py
```

Outputs in `results_rc_step3/`: CSV + PNG for (1) moment–curvature with assumed axial compression 5% of fcAg, (2) single forceBeamColumn 3m cantilever with PDelta and cyclic lateral displacement. `summary.json` reports errors independently. No full-Shell validation.

The results only demonstrate a numerically operational 2D RC fiber section and beam element. They are NOT experimental calibration or structural safety certification.

Model notes: moment-curvature follows OpenSeesPy official zeroLengthSection approach. Beam uses `forceBeamColumn` and Lobatto integration. Cyclic displacement is imposed quasi-statically, not actual earthquake ground motion. Units are N, m, Pa. 

If numerical failure occurs, share full terminal log, `summary.json` and OpenSees version. Do not mark FAIL as PASS.

Reference:
- https://openseespydoc.readthedocs.io/en/latest/src/MomentCurvature.html
- https://opensees.github.io/OpenSeesDocumentation/user/manual/model/beamIntegrations/Lobatto.html
