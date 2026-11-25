"""Geometry-dependent conversions and independent laminate calculations (SI)."""
import math
import numpy as np
from .materials import material

def _positive(**values):
 for k,v in values.items():
  if not math.isfinite(v) or v<=0:raise ValueError(f'{k} must be finite and positive')
def shielding(thickness_m,area_m2,gravity_m_s2,*,density_kg_m3=1700):
 _positive(thickness_m=thickness_m,area_m2=area_m2,gravity_m_s2=gravity_m_s2,density_kg_m3=density_kg_m3)
 mass=density_kg_m3*thickness_m*area_m2
 return dict(mass_kg=mass,weight_N=mass*gravity_m_s2,areal_mass_kg_m2=density_kg_m3*thickness_m,pressure_vertical_Pa=density_kg_m3*thickness_m*gravity_m_s2,
             convention='thickness measured normal to local actual surface; area is actual surface area; weight acts globally vertical, not surface-normal; mass added once')
from SSI.parameters import contact_parameters, radiation_dashpots

def reduced_Q(m):
 E1,E2,nu,G=m['E1_Pa'],m['E2_Pa'],m['nu12'],m['G12_Pa'];den=1-nu*nu*E2/E1
 if den<=0:raise ValueError('nonpositive orthotropic compliance')
 return np.array([[E1/den,nu*E2/den,0],[nu*E2/den,E2/den,0],[0,0,G]])
def qbar(Q,theta_deg):
 c=math.cos(math.radians(theta_deg));s=math.sin(math.radians(theta_deg));q11,q22,q12,q66=Q[0,0],Q[1,1],Q[0,1],Q[2,2]
 return np.array([[q11*c**4+2*(q12+2*q66)*s*s*c*c+q22*s**4,(q11+q22-4*q66)*s*s*c*c+q12*(s**4+c**4),(q11-q12-2*q66)*c**3*s-(q22-q12-2*q66)*c*s**3],
 [(q11+q22-4*q66)*s*s*c*c+q12*(s**4+c**4),q11*s**4+2*(q12+2*q66)*s*s*c*c+q22*c**4,(q11-q12-2*q66)*c*s**3-(q22-q12-2*q66)*c**3*s],
 [(q11-q12-2*q66)*c**3*s-(q22-q12-2*q66)*c*s**3,(q11-q12-2*q66)*c*s**3-(q22-q12-2*q66)*c**3*s,(q11+q22-2*q12-2*q66)*s*s*c*c+q66*(s**4+c**4)]])
def laminate_abd(angles_deg,thicknesses_m=None,*,name='cfrp_as4_8552_ud_rtd'):
 if not len(angles_deg):raise ValueError('nonempty stacking sequence required')
 m=material(name);t=np.asarray(thicknesses_m if thicknesses_m is not None else [m['ply_thickness_m']]*len(angles_deg),float)
 if len(t)!=len(angles_deg) or not np.all(np.isfinite(t)) or np.any(t<=0) or not np.all(np.isfinite(angles_deg)):raise ValueError('invalid plies')
 z=np.r_[-sum(t)/2,-sum(t)/2+np.cumsum(t)];A=np.zeros((3,3));B=A.copy();D=A.copy();Q=reduced_Q(m)
 for k,theta in enumerate(angles_deg):
  qb=qbar(Q,theta);A+=qb*(z[k+1]-z[k]);B+=qb*(z[k+1]**2-z[k]**2)/2;D+=qb*(z[k+1]**3-z[k]**3)/3
 return dict(A_N_m=A,B_N=B,D_N_m=D,ABD=np.block([[A,B],[B,D]]),z_m=z,areal_mass_kg_m2=m['density_kg_m3']*sum(t),angles_deg=list(angles_deg),thicknesses_m=t)
def ply_screen(angles_deg,membrane_strain,curvature_per_m,thicknesses_m=None,*,name='cfrp_as4_8552_ud_rtd'):
 m=material(name);lam=laminate_abd(angles_deg,thicknesses_m,name=name);Q=reduced_Q(m);out=[]
 e0=np.asarray(membrane_strain,float);kappa=np.asarray(curvature_per_m,float)
 if e0.shape!=(3,) or kappa.shape!=(3,) or not np.all(np.isfinite(np.r_[e0,kappa])):raise ValueError('three finite engineering strain/curvature components required')
 for i,a in enumerate(angles_deg):
  c=math.cos(math.radians(a));s=math.sin(math.radians(a));T=np.array([[c*c,s*s,c*s],[s*s,c*c,-c*s],[-2*c*s,2*c*s,c*c-s*s]])
  for side,z in [('bottom',lam['z_m'][i]),('top',lam['z_m'][i+1])]:
   stress=Q@T@(e0+z*kappa);x,y,t=stress
   ratios=[abs(x)/(m['Xt_mean_Pa'] if x>=0 else m['Xc_mean_Pa']),abs(y)/(m['Yt_mean_Pa'] if y>=0 else m['Yc_mean_Pa']),abs(t)/m['S_offset02_Pa']]
   out.append(dict(ply=i,side=side,stress_12_Pa=stress.tolist(),max_ratio=max(ratios),scope='mean max-stress/0.2%-shear-onset screen, NOT design safety factor or progressive failure'))
 return out
