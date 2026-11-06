"""Pure curvature shell energy vs independent classical laminate D matrix."""
import json,sys,math
from pathlib import Path
import openseespy.opensees as o
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.laminate import register_crossply
from capmat.mechanics import laminate_abd

def case(angles,sub):
 o.wipe();o.model('basic','-ndm',3,'-ndf',6);register_crossply(o,10,100,angles,sublayers_per_ply=sub)
 for i,(x,y) in enumerate([(0,0),(1,0),(1,1),(0,1)],1):o.node(i,x,y,0.)
 o.element('ShellMITC4',1,1,2,3,4,10)
 o.node(99,0.,0.,0.);o.fix(99,0,1,1,1,1,1);o.uniaxialMaterial('Elastic',99,1.);o.element('zeroLength',99,1,99,'-mat',99,'-dir',1)
 o.timeSeries('Linear',1);o.pattern('Plain',1,1);k=.001;vals={}
 for i,(x,y) in enumerate([(0,0),(1,0),(1,1),(0,1)],1):
  disp=[0.,0.,k*x*x/2,0.,-k*x,0.];vals[i]=disp
  for dof,value in enumerate(disp,1):o.sp(i,dof,value)
 o.constraints('Transformation');o.numberer('RCM');o.system('BandGeneral');o.test('NormUnbalance',1e-8,30);o.algorithm('Newton');o.integrator('LoadControl',1.);o.analysis('Static');ok=o.analyze(1);o.reactions()
 energy=.5*sum(sum(u*f for u,f in zip(vals[n],o.nodeReaction(n))) for n in vals);expected=.5*k*k*laminate_abd(angles)['D_N_m'][0,0];err=abs(energy-expected)/expected
 return dict(angles=angles,sublayers=sub,energy_J=energy,expected_J=float(expected),relative_error=float(err),**{'pass':bool(ok==0 and err<.002)})
r=dict(runtime=o.version(),cases=[case(a,32) for a in [[0],[90],[0,90,90,0]]]);r['all_pass']=all(c['pass'] for c in r['cases']);(ROOT/'results'/'laminate_bending.json').write_text(json.dumps(r,indent=2));print(r)
