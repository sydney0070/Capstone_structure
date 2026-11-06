"""Elastic bonded RC continuum for foundation/skin initial-stiffness comparisons.
Steel volume fractions along global/material x,y,z are explicit geometry inputs.
Voigt strain compatibility: concrete occupies1-sum(rho_i), each bar family adds
axial Es*rho_i only. No smeared cracking, bond slip, yielding, or confinement law.
"""
import numpy as np
from .materials import material,_tag

def register_elastic_rc(ops,tag,reinforcement_volume_fractions):
 _tag(tag);ratios=np.asarray(reinforcement_volume_fractions,float)
 if ratios.shape!=(3,) or not np.all(np.isfinite(ratios)) or np.any(ratios<0) or sum(ratios)>=1:raise ValueError('three nonnegative steel volume fractions with sum<1')
 c=material('concrete_c30_prototype');s=material('rebar_sd400_prototype');fraction=1-float(sum(ratios));E=c['E_Pa'];nu=c['nu'];G=c['G_Pa'];lam=E*nu/((1+nu)*(1-2*nu))
 Cn=fraction*(np.full((3,3),lam)+np.eye(3)*2*G)+np.diag(ratios*s['params']['E']);S=np.linalg.inv(Cn)
 Ex,Ey,Ez=1/np.diag(S);nuxy=-S[0,1]/S[0,0];nuyz=-S[1,2]/S[1,1];nuzx=-S[2,0]/S[2,2];Gs=fraction*G
 density=fraction*c['density_kg_m3']+float(sum(ratios))*s['density_kg_m3']
 ops.nDMaterial('ElasticOrthotropic',tag,float(Ex),float(Ey),float(Ez),float(nuxy),float(nuyz),float(nuzx),Gs,Gs,Gs,density)
 return dict(tag=tag,reinforcement_volume_fractions=ratios.tolist(),normal_stiffness_Pa=Cn.tolist(),shear_modulus_Pa=Gs,density_kg_m3=density,mass_in_material=True,scope='initial elastic bonded RC Voigt reference, x/y/z bar families; NOT nonlinear RC solid or cracked stiffness; do not add steel mass twice')
