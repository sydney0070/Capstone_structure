# DP cohesion conversion: equation and executable check

OpenSees' implemented perfect-plastic yield equation is
`||s|| + rho*I1 - sqrt(2/3)*sigmaY = 0` (tension positive).
For triaxial compression with positive compressive magnitudes S1,S3,
`||s||=sqrt(2/3)*(S1-S3)` and `I1=-(S1+2*S3)`.
Let `a=rho*sqrt(3/2)=2*sin(phi)/(3-sin(phi))`. Then
`S1=(1+2*a)/(1-a)*S3 + sigmaY/(1-a)`.
Matching `S1=Nphi*S3+2*c*sqrt(Nphi)` gives
`rho=2*sqrt(2/3)*sin(phi)/(3-sin(phi))`,
**`sigmaY=6*c*cos(phi)/(3-sin(phi))`**.

The [online documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html)
prints an additional sqrt(2) in sigmaY's denominator, inconsistent with that invariant normalization.
The [implementation](https://github.com/OpenSees/OpenSees/blob/master/SRC/material/nD/UWmaterials/DruckerPrager.cpp)
uses `root23*Kiso` and `Kiso(0)=sigmaY`. The adapter uses the derived conversion.

`dp_document_mapping_before_correction.json` preserves the failed check:
for JEZ-1 peak, c=3100Pa, phi40.6deg and S3=15000Pa,
the printed-doc mapping gives80435.94Pa instead of MC84384.08Pa.
Cohesive contribution is reduced by1/sqrt(2). This discrepancy vanishes when c=0;
therefore a cohesionless-only test would miss it. The corrected converter is tested
at3 confinements with both c=0 and c>0 in `verification/soil_material.py`.
This establishes the mapping, NOT direct-shear/triaxial equivalence of a soil.
