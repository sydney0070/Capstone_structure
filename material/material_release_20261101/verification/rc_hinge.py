"""Production RC hinge formulation: fiber/time convergence and axial preload.
Uses full path and force-free lateral unloading, not peak-only screening.
"""
import json,sys
from pathlib import Path
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'verification')]
from capmat.rc import register_rc_hinge_beam
from capmat.analysis import advance_displacement
import rc_convergence as c
OUT=ROOT/'results'/'rc_hinge';OUT.mkdir(exist_ok=True)
BARS=[(y,z,.02) for y in [-.245,.245] for z in [-.145,0,.145]]
LP=.3;P=0.;META={}
def model(ny,nip,form='force',ne=1):
 global META
 ops.wipe();ops.model('basic','-ndm',2,'-ndf',3);ops.node(1,0.,0.);ops.node(2,0.,3.);ops.fix(1,1,1,1);ops.geomTransf('PDelta',1)
 META=register_rc_hinge_beam(ops,1,1,2,1,1,10,20,1,2,width_m=.4,depth_m=.6,bars=BARS,hinge_i_m=LP,hinge_j_m=LP,ny=ny,nz=2)
 ops.constraints('Plain');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',.05,100,0);ops.algorithm('Newton');ops.integrator('LoadControl',.05);ops.analysis('Static')
 if P:
  ops.timeSeries('Linear',9);ops.pattern('Plain',9,9);ops.load(2,0.,-P,0.)
  if ops.analyze(20):raise RuntimeError('axial preload failed')
  ops.loadConst('-time',0.)
 ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(2,1.,0.,0.)
 return 2

def snapshot(top,segment,fraction,phase):
 ops.reactions();u,v,_=ops.nodeDisp(top);f=ops.getLoadFactor(1);rx,ry,m=ops.nodeReaction(1)
 return dict(phase=phase,segment=segment,fraction=fraction,u_m=u,v_m=v,force_N=f,rx_N=rx,ry_N=ry,mz_Nm=m,force_error_N=abs(rx+f),vertical_error_N=abs(ry-P),moment_error_Nm=abs(m-3*f-P*u))
c.model=model;c.snapshot=snapshot
if __name__=='__main__':
 result=dict(runtime=ops.version(),scope='3m cantilever,0.4x0.6m,6D20 unconfined nominalC30/SD400, PDelta small-rotation, fixed explicit hinge length. No population/experimental validation. Prior arbitrary-Lobatto collapse release rejected.',cases=[],convergence=[],hinge_sensitivity=[])
 runs=[]
 for ny,steps,lp,ratio in [(24,60,.3,0),(36,60,.3,0),(48,60,.3,0),(36,120,.3,0),(36,60,.15,0),(36,60,.35,0),(36,60,.3,.05),(36,60,.3,.10),(36,120,.3,.10),(36,60,.3,.20)]:
  LP=lp;P=ratio*.24*30e6;c.OUT=OUT/f'lp{lp}_axial{ratio}';c.OUT.mkdir(exist_ok=True)
  stats={}
  def balance():
   r=snapshot(2,0,0,'internal');scale=max(abs(r['force_N']),P,1000.)
   if max(r['force_error_N'],r['vertical_error_N'],r['moment_error_Nm']/3)>1e-4*scale:raise RuntimeError('internal adaptive-step equilibrium failure')
  def advance(top,desired):return advance_displacement(ops,top,1,desired,on_commit=balance,statistics=stats)
  info,rows=c.run(ny,6,steps,advance=advance);info['adaptive_steps']=stats;info['hinge_length_m']=lp;info['axial_load_N']=P;info['axial_ratio_gross_fc']=ratio;info['production_section']=META;info.pop('section_invariants',None);info.pop('section',None)
  runs.append((info,rows));result['cases'].append(info)
  (OUT/'summary.json').write_text(json.dumps(result,indent=2))
 for i,j in [(0,1),(1,2),(1,3),(7,8)]:result['convergence'].append(c.compare(runs[i],runs[j]))
 result['hinge_sensitivity']=[dict(lp_m=runs[i][0]['hinge_length_m'],peak_force_N=runs[i][0]['peak_force_N'],residual_m=runs[i][0]['residual_m']) for i in [4,1,5]]
 result['all_pass']=all(i['equilibrium_pass'] for i in result['cases']) and all(i['pass'] for i in result['convergence'])
 (OUT/'summary.json').write_text(json.dumps(result,indent=2))
