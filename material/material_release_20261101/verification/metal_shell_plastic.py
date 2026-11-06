"""J2 perfect plastic shell coupon vs independent1D return map, including reversals."""
import json,sys
from pathlib import Path
import numpy as np
import openseespy.opensees as o
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat import register_plate,material
r=dict(runtime=o.version(),cases=[])
for name in ['steel_s355_nominal','al6061_t6_extrusion']:
 o.wipe();o.model('basic','-ndm',3,'-ndf',6);register_plate(o,name,100,101,10,.01);m=material(name)
 for n,(x,y) in enumerate([(0,0),(1,0),(1,1),(0,1)],1):o.node(n,x,y,0.);o.fix(n,int(x==0),int(y==0),1,1,1,1)
 o.element('ShellMITC4',1,1,2,3,4,10);targets=[0.,.006,-.006,.003,0.];x=np.r_[0.,np.concatenate([np.linspace(a,b,101)[1:] for a,b in zip(targets[:-1],targets[1:])])]
 o.timeSeries('Path',1,'-dt',1.,'-values',*x.tolist(),'-useLast');o.pattern('Plain',1,1)
 for n in [2,3]:o.sp(n,1,1.)
 o.constraints('Transformation');o.numberer('RCM');o.system('BandGeneral');o.test('NormUnbalance',.001,100);o.algorithm('Newton');o.integrator('LoadControl',1.);o.analysis('Static')
 ep=0.;rows=[];ok=0
 for e in x[1:]:
  ok=o.analyze(1)
  if ok:break
  s=float(np.clip(m['E_Pa']*(e-ep),-m['fy_Pa'],m['fy_Pa']));ep=e-s/m['E_Pa'];o.reactions();actual=-sum(o.nodeReaction(n,1) for n in [1,4])/.01;rows.append([e,s,actual])
 err=max([abs(a-b)/m['fy_Pa'] for e,a,b in rows] or [1.]);r['cases'].append(dict(name=name,steps=len(rows),max_error_over_yield=err,**{'pass':bool(ok==0 and len(rows)==400 and err<2e-4)}))
r['all_pass']=all(c['pass'] for c in r['cases']);(ROOT/'results'/'metal_shell_plastic.json').write_text(json.dumps(r,indent=2));print(r)
