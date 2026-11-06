"""Rigid-translation eigenvalue validates actual shell self mass, including CFRP."""
import json,sys,math
from pathlib import Path
import openseespy.opensees as o
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat import register_plate
from capmat.laminate import register_crossply
from capmat.mass import add_translational_mass
r=dict(runtime=o.version(),cases=[])
for name in ['steel_s355_nominal','al6061_t6_extrusion','sintered_jsc1a_high_porosity','cfrp']:
 o.wipe();o.model('basic','-ndm',3,'-ndf',6)
 m=register_crossply(o,10,100,[0,90,90,0]) if name=='cfrp' else register_plate(o,name,100,101,10,.01)
 for n,(x,y) in enumerate([(0,0),(1,0),(1,1),(0,1)],1):o.node(n,x,y,0.);o.fix(n,0,1,1,1,1,1)
 for n in [2,3,4]:o.equalDOF(1,n,1)
 o.element('ShellMITC4',1,1,2,3,4,10)
 if not m['mass_in_material']:add_translational_mass(o,[1,2,3,4],m['areal_mass_kg_m2'],source_key='shell1_self',ledger={})
 o.node(99,0.,0.,0.);o.fix(99,1,1,1,1,1,1);o.uniaxialMaterial('Elastic',99,1e6);o.element('zeroLength',99,99,1,'-mat',99,'-dir',1)
 o.constraints('Plain');o.numberer('RCM');o.system('FullGeneral');o.algorithm('Linear');o.integrator('LoadControl',0.);o.analysis('Static');eig=o.eigen('-fullGenLapack',1)[0];actual=1e6/eig;expected=m['areal_mass_kg_m2'];err=abs(actual-expected)/expected
 r['cases'].append(dict(name=name,computed_kg=actual,expected_kg=expected,relative_error=err,**{'pass':bool(math.isfinite(eig) and err<1e-7)}))
r['all_pass']=all(x['pass'] for x in r['cases']);(ROOT/'results'/'shell_mass.json').write_text(json.dumps(r,indent=2));print('shell_mass',r['all_pass'])
