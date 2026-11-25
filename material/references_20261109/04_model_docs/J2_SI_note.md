# J2Plasticity SI diagnostic

The tested OpenSees3.4 engine stayed elastic above yield in the Pa-unit shell coupon.
The original result is preserved in `J2Plasticity_SI_rejected.json`.
The [J2Plasticity implementation](https://github.com/OpenSees/OpenSees/blob/master/SRC/material/nD/J2Plasticity.cpp)
sets tolerance=1e-8*sigma0 and initial residual=1 before a while(abs(residual)>tolerance) loop.
Thus fy>=100MPa in Pa can skip the plastic multiplier iteration.

The release uses the zero-friction DruckerPrager limit:
rho=rhoBar=H=Kinf=Ko=0 and sigmaY=fy, yielding
f=||s||-sqrt(2/3)*fy. This is the selected perfect-plastic von Mises law,
with density supplied by that material. `metal_shell_plastic.py` checks400 strain
increments and reversals for each metal against an independent1D return map.
The very large DP tensile invariant cutoff remains a limitation outside the
ordinary low-pressure scope. This is not a calibrated metal damage model.
