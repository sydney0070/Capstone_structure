"""Rigid-translation solid mass check from an independent SDOF eigenfrequency.
All8 brick nodes translate together; added springk gives omega²=k/(rhoV).
"""
import json,sys,traceback,math
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.materials import register_nd,material
from capmat.mass import attach_material_volume_mass,add_translational_mass
XYZ=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
def case(name,added_mass=0.):
 ops.wipe();ops.model('basic','-ndm',3,'-ndf',3);meta=register_nd(ops,name,1);ledger={}
 for i,xyz in enumerate(XYZ,1):ops.node(i,*xyz);ops.fix(i,0,1,1)
 for i in range(2,9):ops.equalDOF(1,i,1)
 ops.element('bbarBrick',1,*range(1,9),1,0.,0.,0.)
 if meta['external_mass_required']:attach_material_volume_mass(ops,list(range(1,9)),1.,meta,source_key='solid1_self',ledger=ledger)
 if added_mass:add_translational_mass(ops,list(range(1,9)),added_mass,source_key='shield',ledger=ledger)
 ops.node(99,0.,0.,0.);ops.fix(99,1,1,1);ops.uniaxialMaterial('Elastic',99,1e6);ops.element('zeroLength',99,99,1,'-mat',99,'-dir',1)
 ops.constraints('Transformation');ops.numberer('RCM');ops.system('FullGeneral');ops.algorithm('Linear');ops.integrator('LoadControl',0.);ops.analysis('Static')
 lam=ops.eigen('-fullGenLapack',1)[0];mass=meta['density_kg_m3']+added_mass;expected=1e6/mass;err=abs(lam-expected)/expected
 return dict(name=name,added_mass_kg=added_mass,eigenvalue_per_s2=lam,expected_per_s2=expected,mass_from_eigen_kg=1e6/lam,expected_mass_kg=mass,relative_error=err,ledger=ledger,**{'pass':bool(math.isfinite(lam) and err<1e-7)})
if __name__=='__main__':
 r=dict(runtime=ops.version(),cases=[])
 for name in ['steel_s355_nominal','al6061_t6_extrusion','sintered_jsc1a_high_porosity','moon_jsc1a_18_peak','mars_jez1_peak','earth_sand_dr50_dry_reference']:
  for extra in [0.,1700.]:
   try:c=case(name,extra)
   except Exception as e:c=dict(name=name,error=str(e),traceback=traceback.format_exc(),**{'pass':False})
   r['cases'].append(c);print(c,flush=True)
 r['all_pass']=all(c['pass'] for c in r['cases']);(ROOT/'results'/'dynamic_mass.json').write_text(json.dumps(r,indent=2))
