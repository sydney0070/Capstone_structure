"""Exact principal-axis0/90 ply registration; arbitrary-angle ABD in mechanics.
Angles refer to EACH shell element's local x axis. No hidden isotropization.
"""
import math
import numpy as np
from .materials import material,_tag
from .mechanics import laminate_abd

def register_crossply(ops,section_tag,material_tag_start,angles_deg,thicknesses_m=None,*,name='cfrp_as4_8552_ud_rtd',sublayers_per_ply=32):
 _tag(section_tag);_tag(material_tag_start)
 lam=laminate_abd(angles_deg,thicknesses_m,name=name);m=material(name)
 if type(sublayers_per_ply) is not int or sublayers_per_ply<1 or len(angles_deg)*sublayers_per_ply<3:raise ValueError('at least3 integration layers; positive integer sublayers per ply')
 if any(min(abs(a%180),abs(a%180-90),abs(a%180-180))>1e-9 for a in angles_deg):raise ValueError('OpenSees adapter supports0/90 only; arbitrary-angle ABD calculator is separate; no silent approximate rotation')
 E1,E2,G12,nu12=m['E1_Pa'],m['E2_Pa'],m['G12_Pa'],m['nu12'];E3=E2;nu23=m['nu23'];G23=E2/(2*(1+nu23));nu31=nu12*E3/E1
 # Symmetric compliance; reciprocity is essential for the minor Poisson ratio.
 S=np.diag([1/E1,1/E2,1/E3,1/G23,1/G12,1/G12]);S[0,1]=S[1,0]=-nu12/E1;S[0,2]=S[2,0]=-nu31/E3;S[1,2]=S[2,1]=-nu23/E2
 if np.min(np.linalg.eigvalsh(S))<=0:raise ValueError('orthotropic compliance not positive definite')
 layers=[]
 for i,a in enumerate(angles_deg):
  tag=material_tag_start+2*i;ptag=tag+1
  if abs((a%180)-90)<1e-9:
   args=[E2,E1,E3,nu12*E2/E1,nu12,nu23*E3/E2,G12,G12,G23]
  else:args=[E1,E2,E3,nu12,nu23,nu31,G12,G23,G12]
  ops.nDMaterial('ElasticOrthotropic',tag,*args,m['density_kg_m3']);ops.nDMaterial('PlateFiber',ptag,tag)
  for _ in range(sublayers_per_ply):layers.extend([ptag,float(lam['thicknesses_m'][i])/sublayers_per_ply])
 ops.section('LayeredShell',section_tag,len(angles_deg)*sublayers_per_ply,*layers)
 return dict(section_tag=section_tag,areal_mass_kg_m2=float(lam['areal_mass_kg_m2']),material_tags=list(range(material_tag_start,material_tag_start+2*len(angles_deg))),mass_in_material=True,angles_deg=list(angles_deg),sublayers_per_ply=sublayers_per_ply,local_axis_requirement='orient element localx consistently with declared laminate0deg; z bottom-to-top must match ply sequence')
