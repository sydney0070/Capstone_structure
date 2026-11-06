"""Contact normal orientation, nodal area scaling, shear capacity, opening/reclosure."""
import csv,json,sys,traceback
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.ssi import register_contact
from capmat.mechanics import contact_parameters,radiation_dashpots
OUT=ROOT/'results'/'ssi_contact';OUT.mkdir(exist_ok=True)
def case(area,cohesion,axis=3):
 ops.wipe();ops.model('basic','-ndm',3,'-ndf',3)
 for n in [1,2,99]:ops.node(n,0.,0.,0.)
 ops.fix(1,1,1,1);ops.fix(99,0,1,1)
 tangent=1 if axis!=1 else 2
 ops.fix(2,*[int(d not in [axis,tangent]) for d in [1,2,3]])
 ops.uniaxialMaterial('Elastic',99,1.);ops.element('zeroLength',99,1,99,'-mat',99,'-dir',1)
 p=contact_parameters(area,30e6,12e6,1,1e5,.0001,.34,adhesion_Pa=cohesion)
 register_contact(ops,1,2,1,p,normal_axis=axis,auxiliary_nodes=(101,102),normal_material_tag=10,rotation_material_tag=11,friction_tag=1)
 endpoints=[(0.,0.),(0.,-.0001),(.002,-.0001),(-.002,-.0001),(-.002,.0002),(.002,.0002),(.002,-.00005),(-.002,-.00005)]
 points=[endpoints[0]];phases=[0]
 for j,(start,end) in enumerate(zip(endpoints[:-1],endpoints[1:]),1):
  for t in np.linspace(0,1,41)[1:]:points.append(tuple((1-t)*np.array(start)+t*np.array(end)));phases.append(j)
 arr=np.array(points)
 for tag,d,values in [(1,tangent,arr[:,0]),(2,axis,arr[:,1])]:
  ops.timeSeries('Path',tag,'-dt',1.,'-values',*values.tolist(),'-useLast');ops.pattern('Plain',tag,tag);ops.sp(2,d,1.)
 ops.constraints('Transformation');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',1e-7,30);ops.algorithm('Newton');ops.integrator('LoadControl',1.);ops.analysis('Static')
 rows=[]
 for i in range(1,len(points)):
  if ops.analyze(1):raise RuntimeError(f'contact step{i}')
  ops.reactions();ef=np.array(ops.eleForce(1));r=ef[6:9];g=ef[:3];N=-r[axis-1];T=r[tangent-1]
  if max(abs(np.array(ops.nodeDisp(102)[:3])-ops.nodeDisp(2)))>1e-10:raise RuntimeError('MP displacement transfer failed')
  rows.append(dict(phase=phases[i],ux=arr[i,0],un=arr[i,1],normal_N=N,tangent_N=T,expected_normal_N=p['Kn_N_m']*max(-arr[i,1],0),balance_N=float(np.linalg.norm(r+g))))
 err=max(abs(r['normal_N']-r['expected_normal_N']) for r in rows)/max(p['Kn_N_m']*.0001,1)
 bound=max(max(0,abs(r['tangent_N'])-(p['mu']*max(r['normal_N'],0)+p['cohesive_force_N'])) for r in rows)
 opening=max(max(abs(r['normal_N']),abs(r['tangent_N'])) for r in rows if r['un']>1e-9)
 slide=[r for r in rows if (r['phase']==2 and r['ux']>.0015) or (r['phase']==3 and r['ux']<-.0015)]
 slide_error=max(abs(abs(r['tangent_N'])-(p['mu']*r['normal_N']+p['cohesive_force_N'])) for r in slide)
 balance=max(r['balance_N'] for r in rows)
 name=f'area{area}_adhesion{cohesion}_axis{axis}'
 with (OUT/(name+'.csv')).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 return dict(name=name,parameters=p,normal_relative_error=err,capacity_excess_N=bound,opening_force_N=opening,sliding_capacity_error_N=slide_error,action_reaction_error_N=balance,**{'pass':bool(err<1e-8 and max(bound,opening,slide_error,balance)<1e-5)})
def viscous():
 from capmat.ssi import register_dashpot
 r=radiation_dashpots(1800,30e6,12e6,2);C=r['Ct_N_s_m'];data=[]
 for v in [.01,-.01,.2,-.2]:
  ops.wipe();ops.model('basic','-ndm',3,'-ndf',3)
  for n in [1,2,99]:ops.node(n,0.,0.,0.)
  ops.fix(1,1,1,1);ops.fix(2,0,1,1);ops.fix(99,0,1,1);ops.mass(99,1.,0.,0.)
  ops.uniaxialMaterial('Elastic',99,1.);ops.element('zeroLength',99,1,99,'-mat',99,'-dir',1)
  register_dashpot(ops,1,1,2,1,C,1)
  ops.timeSeries('Path',1,'-dt',.01,'-values',*[v*.01*k for k in range(12)],'-useLast');ops.pattern('Plain',1,1);ops.sp(2,1,1.)
  ops.setNodeVel(2,1,v,'-commit')
  ops.constraints('Transformation');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',1e-8,30);ops.algorithm('Newton');ops.integrator('Newmark',.5,.25);ops.analysis('Transient')
  for k in range(10):
   if ops.analyze(1,.01):raise RuntimeError('dashpot transient failed')
   ops.reactions('-dynamic');data.append([v,ops.nodeVel(2,1),ops.nodeReaction(2,1),C*v])
 return dict(name='linear_radiation_dashpot_transient',samples=data,**{'pass':bool(max(abs(f-expect) for v,actual,f,expect in data)<1e-6 and max(abs(v-a) for v,a,f,e in data)<1e-10)},scope='imposed constant-velocity transient force; far-boundary absorption needs domain-size/frequency check')
if __name__=='__main__':
 result=dict(runtime=ops.version(),cases=[])
 for area,c,axis in [(1,0,3),(2,0,3),(1,0,1),(1,0,2)]:
  try:r=case(area,c,axis)
  except Exception as e:r=dict(error=str(e),traceback=traceback.format_exc(),**{'pass':False})
  result['cases'].append(r);print(r,flush=True)
 result['cases'].append(viscous());result['all_pass']=all(r['pass'] for r in result['cases']);(OUT/'summary.json').write_text(json.dumps(result,indent=2))
