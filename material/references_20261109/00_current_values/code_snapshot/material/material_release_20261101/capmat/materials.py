"""OpenSees registration in m,N,s,kg,Pa. Explicit tags, no wipe, no fallback laws."""
import copy,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def catalog():return json.loads((ROOT/'configs'/'catalog.json').read_text())
def material(name):return copy.deepcopy(catalog()['materials'][name])
def _tag(tag):
 if type(tag) is not int or tag<=0:raise ValueError('tag must be a positive integer')
def register_uniaxial(ops,name,tag,*,elastic=False):
 _tag(tag);m=material(name);law=m['law']
 if elastic:
  E=m.get('E_Pa',m.get('params',{}).get('E'))
  if E is None:raise ValueError('No single isotropic axial modulus: specify ply axis or constituent')
  ops.uniaxialMaterial('Elastic',tag,E)
 elif law=='Steel02':
  p=m['params'];args=[p[k] for k in ['Fy','E','b','R0','cR1','cR2']]
  if 'a1' in p:args += [p[k] for k in ['a1','a2','a3','a4']]
  ops.uniaxialMaterial('Steel02',tag,*args)
 elif law=='J2PerfectPlastic':ops.uniaxialMaterial('Steel01',tag,m['fy_Pa'],m['E_Pa'],0.)
 elif law=='Concrete02':
  ops.uniaxialMaterial('Concrete02',tag,*[m['params'][k] for k in ['fpc','epsc0','fpcu','epsU','lam','ft','Ets']])
 elif law=='ElasticIsotropic':ops.uniaxialMaterial('Elastic',tag,m['E_Pa'])
 else:raise ValueError(f'{name}: no verified uniaxial mapping for {law}')
 return dict(tag=tag,density_kg_m3=m['density_kg_m3'],mass_in_material=False,law='Elastic' if elastic else law)
def dp_constants(phi_deg,c_Pa,psi_deg):
 """Match implemented f=||s||+rho I1-sqrt(2/3)*sigmaY to MC compression.
 sigmaY=6c cos(phi)/(3-sin(phi)); the online doc's extra sqrt(2)
 denominator reduces the cohesive contribution to 1/sqrt(2). Derived
 from the invariant equation and verified against OpenSees3.4 triaxial
 tests. See research/dp_mapping_note.md. psi mapping is a convention.
 """
 if not 0<phi_deg<89 or not 0<=psi_deg<=phi_deg or c_Pa<0:raise ValueError('invalid friction/cohesion/dilation')
 def slope(a):s=math.sin(math.radians(a));return 2*math.sqrt(2)*s/(math.sqrt(3)*(3-s))
 p=math.radians(phi_deg)
 return dict(rho=slope(phi_deg),rhoBar=slope(psi_deg),sigmaY=6*c_Pa*math.cos(p)/(3-math.sin(p)))
def register_nd(ops,name,tag,*,elastic=False,stiffness_factor=1.):
 _tag(tag);m=material(name);law=m['law']
 if not math.isfinite(stiffness_factor) or stiffness_factor<=0:raise ValueError('positive stiffness factor required')
 rho=m['density_kg_m3'];mass=True
 if elastic or law=='ElasticIsotropic':
  if 'E_Pa' not in m:raise ValueError('No isotropic elastic modulus in record')
  ops.nDMaterial('ElasticIsotropic',tag,m['E_Pa']*stiffness_factor,m['nu'],rho)
 elif law=='J2PerfectPlastic':
  # J2Plasticity3.4 initializes local residual=1 but tolerance=1e-8*fy;
  # for fy>=100MPa in Pa its return loop is skipped. Use the exact
  # zero-friction DP limit f=||s||-sqrt(2/3)*fy, H=0 instead.
  ops.nDMaterial('DruckerPrager',tag,m['K_Pa']*stiffness_factor,m['G_Pa']*stiffness_factor,m['fy_Pa'],0.,0.,0.,0.,0.,0.,0.,1.,rho,101000.)
  mass=True
 elif law=='DruckerPrager':
  p=dp_constants(m['phi_deg'],m['c_Pa'],m['psi_deg'])
  ops.nDMaterial('DruckerPrager',tag,m['K_Pa']*stiffness_factor,m['G_Pa']*stiffness_factor,p['sigmaY'],p['rho'],p['rhoBar'],0.,0.,0.,0.,0.,1.,rho,101000.)
 elif law=='PressureDependMultiYield02':
  ops.nDMaterial('PressureDependMultiYield02',tag,3,rho,m['G_Pa']*stiffness_factor,m['K_Pa']*stiffness_factor,m['phi_deg'],m['peak_shear_strain'],m['p_ref_Pa'],m['pressure_exponent'],m['PTAngle_deg'],m['contrac1'],m['contrac3'],m['dilat1'],m['dilat3'],m['n_surfaces'],m['contrac2'],m['dilat2'],m['liquefac1'],m['liquefac2'],m['e'],m['cs1'],m['cs2'],m['cs3'],m['pa_Pa'],m['numerical_c_Pa'])
 else:raise ValueError(f'{name}: {law} has no isotropic solid law; use laminate/RC constituent API')
 return dict(tag=tag,law='ElasticIsotropic' if elastic else law,density_kg_m3=rho,mass_in_material=mass,external_mass_required=not mass,implementation='DruckerPrager_zero_friction_J2_limit' if law=='J2PerfectPlastic' and not elastic else ('ElasticIsotropic' if elastic else law))
def register_plate(ops,name,material_tag,plate_tag,section_tag,thickness_m,*,elastic=False):
 if not math.isfinite(thickness_m) or thickness_m<=0:raise ValueError('positive thickness required')
 for t in [material_tag,plate_tag,section_tag]:_tag(t)
 if material_tag==plate_tag:raise ValueError('distinct nD tags required')
 meta=register_nd(ops,name,material_tag,elastic=elastic)
 ops.nDMaterial('PlateFiber',plate_tag,material_tag)
 ops.section('PlateFiber',section_tag,plate_tag,thickness_m)
 return {**meta,'section_tag':section_tag,'areal_mass_kg_m2':meta['density_kg_m3']*thickness_m}
