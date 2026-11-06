"""DP triaxial-compression analytical yield; PDMY02 drained cyclic simple shear.
This verifies implemented constitutive laws; it does not prove in-situ calibration.
"""
import csv,json,math,sys,traceback
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.materials import register_nd,material
XYZ=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
OUT=ROOT/'results'/'soil_material';OUT.mkdir(exist_ok=True)
def analysis():
 ops.constraints('Transformation');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',1e-4,100,0);ops.algorithm('Newton');ops.integrator('LoadControl',.1);ops.analysis('Static')
def checked(n=1):
 if ops.analyze(n):raise RuntimeError('analysis failed')
def triaxial(name,pressure):
 ops.wipe();ops.model('basic','-ndm',3,'-ndf',3);register_nd(ops,name,1)
 for tag,(x,y,z) in enumerate(XYZ,1):ops.node(tag,x,y,z);ops.fix(tag,int(x==0),int(y==0),int(z==0))
 ops.element('SSPbrick',1,*range(1,9),1,0.,0.,0.)
 ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
 for tag,(x,y,z) in enumerate(XYZ,1):ops.load(tag,-pressure*x/4,-pressure*y/4,-pressure*z/4)
 analysis();checked(10);ops.loadConst('-time',0.)
 initial=np.array(ops.eleResponse(1,'stress'));u0=ops.nodeDisp(5,3)
 values=u0-np.linspace(0,.03,601)
 ops.timeSeries('Path',2,'-dt',1.,'-values',*values.tolist(),'-useLast');ops.pattern('Plain',2,2)
 for tag in range(5,9):ops.sp(tag,3,1.)
 ops.wipeAnalysis();analysis();ops.integrator('LoadControl',0.);checked();ops.integrator('LoadControl',1.)
 rows=[]
 for k in range(1,len(values)):
  checked();s=np.array(ops.eleResponse(1,'stress'));e=np.array(ops.eleResponse(1,'strain'));ops.reactions()
  rows.append([k,-(ops.nodeDisp(5,3)-u0),*s[:6],*e[:6]])
 m=material(name);sn=math.sin(math.radians(m['phi_deg']));N=(1+sn)/(1-sn);expected=pressure*N+2*m['c_Pa']*math.sqrt(N)
 actual=-rows[-1][4];err=abs(actual-expected)/expected;lat=max(abs(rows[-1][2]+pressure),abs(rows[-1][3]+pressure))/pressure
 with (OUT/f'{name}_{int(pressure)}.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['step','axial_compression_increment','sx','sy','sz','sxy','syz','szx','ex','ey','ez','gxy','gyz','gzx']);w.writerows(rows)
 return dict(name=name,confinement_Pa=pressure,initial_hydrostatic_max_error_Pa=float(max(abs(initial[:3]+pressure))),expected_sigma1_Pa=expected,computed_sigma1_Pa=actual,yield_relative_error=err,lateral_stress_relative_error=lat,pass_=True,**{'pass':bool(err<.005 and lat<1e-5 and np.all(np.isfinite(rows)))})
def pdmy():
 name='earth_sand_dr50_dry_reference';ops.wipe();ops.model('basic','-ndm',3,'-ndf',3);register_nd(ops,name,1)
 for tag,xyz in enumerate(XYZ,1):ops.node(tag,*xyz);ops.fix(tag,*((1,1,1) if tag<=4 else (0,1,0)))
 for tag in [6,7,8]:ops.equalDOF(5,tag,1,3)
 ops.element('SSPbrick',1,*range(1,9),1,0.,0.,0.);ops.updateMaterialStage('-material',1,'-stage',0)
 ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(5,0.,0.,-50000.)
 analysis();checked(10);ops.loadConst('-time',0.);ops.updateMaterialStage('-material',1,'-stage',1);ops.integrator('LoadControl',0);checked()
 ops.timeSeries('Linear',2);ops.pattern('Plain',2,2);ops.load(5,1.,0.,0.)
 targets=[.0001,-.0001,.001,-.001,.005,-.005,.01,-.01,0];rows=[];start=0.
 for target in targets:
  for u in np.linspace(start,target,51)[1:]:
   ops.integrator('DisplacementControl',5,1,float(u-ops.nodeDisp(5,1)));checked();s=np.array(ops.eleResponse(1,'stress'));ops.reactions();f=ops.getLoadFactor(2);r=sum(ops.nodeReaction(n,1) for n in [1,2,3,4]);rows.append([u,s[5],s[2],f,r,float(np.sum(ops.eleResponse(1,'strain')[:3]))])
  start=target
 a=np.array(rows);err=float(np.max(abs(a[:,3]+a[:,4]))/max(np.max(abs(a[:,3])),1));vert=float(np.max(abs(a[:,2]+50000))/50000)
 with (OUT/(name+'_cyclic.csv')).open('w',newline='') as f:w=csv.writer(f);w.writerow(['shear_strain','shear_stress_Pa','vertical_stress_Pa','applied_shear_N','base_reaction_N','volume_strain']);w.writerows(rows)
 return dict(name=name,steps=len(rows),peak_tau_Pa=float(max(abs(a[:,1]))),equilibrium_relative_error=err,vertical_stress_relative_error=vert,**{'pass':bool(np.all(np.isfinite(a)) and err<1e-4 and vert<1e-4)},scope='constant vertical effective stress, free volume change, no pore fluid or liquefaction validation')
if __name__=='__main__':
 result=dict(runtime=ops.version(),cases=[])
 tasks=[(n,p) for n in ['moon_jsc1a_18_peak','moon_jsc1a_18_residual','mars_jez1_peak','mars_jez1_residual'] for p in [15000.,30000.,60000.]]
 for n,p in tasks:
  try:r=triaxial(n,p)
  except Exception as e:r=dict(name=n,pressure_Pa=p,error=str(e),traceback=traceback.format_exc(),**{'pass':False})
  result['cases'].append(r);(OUT/'summary.json').write_text(json.dumps(result,indent=2));print(r,flush=True)
 try:r=pdmy()
 except Exception as e:r=dict(name='pdmy',error=str(e),traceback=traceback.format_exc(),**{'pass':False})
 result['cases'].append(r);result['all_pass']=all(r['pass'] for r in result['cases']);(OUT/'summary.json').write_text(json.dumps(result,indent=2));print(r,flush=True)
