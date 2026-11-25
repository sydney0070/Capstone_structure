"""Geometry-dependent contact penalties and radiation impedances (SI)."""
import math

def _positive(**values):
 for k,v in values.items():
  if not math.isfinite(v) or v<=0:raise ValueError(f'{k} must be finite and positive')
def contact_parameters(area_m2,adjacent_E_Pa,adjacent_G_Pa,element_height_m,max_pressure_Pa,allowed_penetration_m,mu,*,penalty_factor=10.,adhesion_Pa=0.):
 _positive(area_m2=area_m2,adjacent_E_Pa=adjacent_E_Pa,adjacent_G_Pa=adjacent_G_Pa,element_height_m=element_height_m,max_pressure_Pa=max_pressure_Pa,allowed_penetration_m=allowed_penetration_m,penalty_factor=penalty_factor)
 if not math.isfinite(mu) or mu<0 or not math.isfinite(adhesion_Pa) or adhesion_Pa<0:raise ValueError('nonnegative mu and adhesion required')
 kn=max(penalty_factor*adjacent_E_Pa/element_height_m,max_pressure_Pa/allowed_penetration_m)
 kt=penalty_factor*adjacent_G_Pa/element_height_m
 return dict(Kn_N_m=kn*area_m2,Kt_N_m=kt*area_m2,mu=mu,cohesive_force_N=adhesion_Pa*area_m2,area_m2=area_m2,estimated_penetration_m=max_pressure_Pa/kn,
             scope='interface penalty only, NOT foundation impedance; test penalty0.5/1/2 and mesh sensitivity; no cohesive tension bond implied')
def radiation_dashpots(density_kg_m3,K_Pa,G_Pa,area_m2):
 _positive(density_kg_m3=density_kg_m3,K_Pa=K_Pa,G_Pa=G_Pa,area_m2=area_m2)
 vs=math.sqrt(G_Pa/density_kg_m3);vp=math.sqrt((K_Pa+4*G_Pa/3)/density_kg_m3)
 return dict(Vs_m_s=vs,Vp_m_s=vp,Ct_N_s_m=density_kg_m3*vs*area_m2,Cn_N_s_m=density_kg_m3*vp*area_m2,
             scope='local plane-wave absorbing boundary approximation, not static foundation spring; assign only on appropriate far boundaries; do not duplicate radiation damping')
