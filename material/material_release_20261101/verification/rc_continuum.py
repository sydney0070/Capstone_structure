import json,sys
from pathlib import Path
import numpy as np
import openseespy.opensees as o
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.elastic_rc import register_elastic_rc
XYZ=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
r=dict(runtime=o.version(),cases=[])
for strain in [[1e-5,0,0,0,0,0],[0,1e-5,0,0,0,0],[0,0,1e-5,0,0,0],[0,0,0,1e-5,0,0]]:
 o.wipe();o.model('basic','-ndm',3,'-ndf',3);m=register_elastic_rc(o,1,[.005,.005,0.])
 for n,xyz in enumerate(XYZ,1):o.node(n,*xyz)
 o.element('SSPbrick',1,*range(1,9),1,0.,0.,0.);o.node(99,0.,0.,0.);o.fix(99,0,1,1);o.uniaxialMaterial('Elastic',99,1.);o.element('zeroLength',99,1,99,'-mat',99,'-dir',1)
 o.timeSeries('Linear',1);o.pattern('Plain',1,1);ex,ey,ez,gxy,gyz,gzx=strain
 for n,(x,y,z) in enumerate(XYZ,1):
  for dof,u in enumerate([ex*x+gxy*y,ey*y+gyz*z,ez*z+gzx*x],1):o.sp(n,dof,u)
 o.constraints('Transformation');o.numberer('RCM');o.system('BandGeneral');o.test('NormUnbalance',1e-8,30);o.algorithm('Linear');o.integrator('LoadControl',1.);o.analysis('Static');ok=o.analyze(1)
 # Independent constituent calculation: concrete stress plus axial steel-family stress.
 Ec=30e9;nu=.2;G=Ec/(2*(1+nu));lam=Ec*nu/((1+nu)*(1-2*nu));e=np.array(strain)
 expected=np.r_[.99*(2*G*e[:3]+lam*sum(e[:3]))+200e9*np.array([.005,.005,0])*e[:3],.99*G*e[3:]];actual=np.array(o.eleResponse(1,'stress')[:6]);err=float(max(abs(actual-expected))/max(max(abs(expected)),1))
 r['cases'].append(dict(strain=strain,stress_Pa=actual.tolist(),expected_Pa=expected.tolist(),relative_error=err,**{'pass':bool(ok==0 and err<1e-8)}))
r['all_pass']=all(x['pass'] for x in r['cases']);(ROOT/'results'/'rc_continuum.json').write_text(json.dumps(r,indent=2));print('rc_continuum',r['all_pass'])
