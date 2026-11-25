"""Gravity and overburden bookkeeping, not a planetary soil property database."""
GRAVITY=dict(earth=9.81,moon=1.62,mars=3.71)

def vertical_stress(planet,rho,depth,surcharge=0.):
    """Dry total/effective overburden kPa. rho is tonne/m³, NOT unit weight.
    For layered soils integrate rho*g*dz. For saturated Earth subtract pore pressure
    separately; do not use this dry formula as an effective-stress groundwater model.
    """
    if rho<=0 or depth<0 or surcharge<0: raise ValueError('invalid overburden inputs')
    return rho*GRAVITY[planet]*depth+surcharge

def nodal_weight(mass_tonne,planet):
    return -mass_tonne*GRAVITY[planet]

def shear_dashpot(rho,Vs,area):
    """First-order absorbing boundary estimate C_s=rho*Vs*A [kN s/m].
    Do not count the same radiation damping in both SSI springs and boundaries.
    """
    if min(rho,Vs,area)<=0: raise ValueError('positive rho, Vs, area required')
    return rho*Vs*area
