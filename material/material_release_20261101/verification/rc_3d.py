"""Biaxial/torsion/axial elastic-limit checks and beam mass for3D RC API."""
import json,sys,math
from pathlib import Path
import openseespy.opensees as o
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.rc import register_rc_hinge_beam
BARS=[(y,z,.02) for y in [-.245,.245] for z in [-.145,0,.145]]
def build():
 o.wipe();o.model('basic','-ndm',3,'-ndf',6);o.node(1,0.,0.,0.);o.node(2,3.,0.,0.);o.geomTransf('Linear',1,0.,0.,1.)
 return register_rc_hinge_beam(o,1,1,2,1,1,10,20,1,2,width_m=.4,depth_m=.6,bars=BARS,hinge_i_m=.3,hinge_j_m=.3,GJ_N_m2=1e7,ny=36,nz=8)
def static():
 m=build();o.fix(1,1,1,1,1,1,1);o.timeSeries('Linear',1);o.pattern('Plain',1,1);o.load(2,100.,100.,50.,100.,0.,0.)
 o.constraints('Plain');o.numberer('RCM');o.system('BandGeneral');o.test('NormUnbalance',1e-6,100);o.algorithm('Newton');o.integrator('LoadControl',.1);o.analysis('Static');ok=o.analyze(10)
 expected=[100*3/m['EA_N'],100*27/(3*m['EIz_N_m2']),50*27/(3*m['EIy_N_m2']),100*3/1e7];u=o.nodeDisp(2);err=[abs(u[i]-e)/abs(e) for i,e in enumerate(expected)]
 return dict(name='elastic_limit_biaxial_axial_torsion',relative_errors=err,computed=u[:4],expected=expected,**{'pass':bool(ok==0 and max(err)<1e-3)})
def mass():
 m=build();o.fix(1,0,1,1,1,1,1);o.fix(2,0,1,1,1,1,1);o.equalDOF(1,2,1);o.node(99,0.,0.,0.);o.fix(99,1,1,1,1,1,1);o.uniaxialMaterial('Elastic',99,1e6);o.element('zeroLength',99,99,1,'-mat',99,'-dir',1)
 o.constraints('Plain');o.numberer('RCM');o.system('FullGeneral');o.algorithm('Linear');o.integrator('LoadControl',0.);o.analysis('Static');lam=o.eigen('-fullGenLapack',1)[0];expected=m['mass_per_length_kg_m']*3;actual=1e6/lam;err=abs(actual-expected)/expected
 return dict(name='rigid_translation_mass',expected_kg=expected,computed_kg=actual,relative_error=err,**{'pass':bool(math.isfinite(lam) and err<1e-7)})
r=dict(runtime=o.version(),cases=[static(),mass()]);r['all_pass']=all(c['pass'] for c in r['cases']);(ROOT/'results'/'rc_3d.json').write_text(json.dumps(r,indent=2));print(r)
