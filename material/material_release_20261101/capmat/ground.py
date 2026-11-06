"""Adapter to the selected user's genuine 3D mesh generator; original is not edited.
Build BEFORE structure: source builder wipes domain. Empty-node guard prevents
accidental destruction of existing structures. Boundaries are static only.
"""
import hashlib,importlib.util,math
from pathlib import Path
import numpy as np
from .materials import material,register_nd,catalog
SOURCE=Path(__file__).resolve().parents[3]/'ground'/'ground_modeling'/'ground_modeling.py'
DEFAULTS={'Earth':'earth_sand_dr50_dry_reference','Moon':'moon_jsc1a_18_peak','Mars':'mars_jez1_peak'}

def initialize_gravity(ops,loads,pdmy_tags,*,n_steps=10,balance_tolerance=1e-4):
 """Step to actual factor1; never replay already-committed gravity increments."""
 if type(n_steps) is not int or n_steps<1:raise ValueError('positive n_steps')
 for tag in pdmy_tags:ops.updateMaterialStage('-material',tag,'-stage',0)
 ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
 total=np.zeros(3)
 for tag,v in loads.items():ops.load(tag,*v);total+=v
 ops.wipeAnalysis();ops.constraints('Transformation');ops.numberer('RCM');ops.system('BandGeneral')
 ops.test('NormUnbalance',max(.001,float(np.linalg.norm(total))*1e-9),100,0);ops.algorithm('Newton');ops.integrator('LoadControl',0.);ops.analysis('Static')
 attempts=[];step=1/n_steps
 while ops.getLoadFactor(1)<1-1e-12:
  before=ops.getLoadFactor(1);delta=min(step,1-before);ops.integrator('LoadControl',delta)
  ok=ops.analyze(1);after=ops.getLoadFactor(1)
  attempts.append(dict(before=before,delta=delta,after=after,ok=ok))
  if ok:
   if abs(after-before)>1e-10:raise RuntimeError('failed step changed committed factor')
   step=delta/2
   if step<1e-6:raise RuntimeError('gravity failed at minimum increment')
  elif abs(after-(before+delta))>1e-9:raise RuntimeError('unexpected gravity factor')
 if abs(ops.getLoadFactor(1)-1)>1e-9:raise RuntimeError('gravity load factor not1')
 ops.loadConst('-time',0.)
 for tag in pdmy_tags:ops.updateMaterialStage('-material',tag,'-stage',1)
 ops.integrator('LoadControl',0.)
 if ops.analyze(1):raise RuntimeError('post-stage-switch re-equilibration failed')
 ops.reactions();reactions=np.sum([ops.nodeReaction(t)[:3] for t in ops.getNodeTags()],axis=0)
 disp=np.array([ops.nodeDisp(t)[:3] for t in ops.getNodeTags()]);err=float(np.linalg.norm(reactions+total)/max(np.linalg.norm(total),1))
 if not np.all(np.isfinite(disp)) or err>balance_tolerance:raise RuntimeError(f'gravity balance/finite gate failed: {err}')
 result=dict(analysis_ok=0,analysis_steps=sum(a['ok']==0 for a in attempts),load_step=1/n_steps,final_factor=ops.getLoadFactor(1),total_applied_load=total.tolist(),reaction_vector=reactions.tolist(),balance_vector=(reactions+total).tolist(),balance_norm=float(np.linalg.norm(reactions+total)),relative_balance_error=err,max_displacement_m=float(np.max(np.linalg.norm(disp,axis=1))),mean_displacement_m=float(np.mean(np.linalg.norm(disp,axis=1))),max_vertical_displacement_m=float(np.max(abs(disp[:,2]))),stage_after_gravity=1,post_stage_re_equilibrated=True,attempts=attempts)
 ops.wipeAnalysis();return result

def build_ground(X,Y,Z,*,body,depth_m,n_depth,name=None,elastic=False,stiffness_factor=1.,source=SOURCE):
 import openseespy.opensees as ops
 if ops.getNodeTags() or ops.getEleTags():raise RuntimeError('Build ground first in a fresh domain: source builder wipes domain')
 if body not in DEFAULTS:raise ValueError('body Earth/Moon/Mars')
 if not math.isfinite(depth_m) or depth_m<=0 or type(n_depth) is not int or n_depth<1:raise ValueError('positive depth/n_depth')
 name=name or DEFAULTS[body];m=material(name)
 if name=='earth_sand_dr50_pdmy02':raise ValueError('saturated-density record needs pore-pressure/effective-stress initialization; dry wrapper cannot assume it')
 if m['law'] not in ['DruckerPrager','PressureDependMultiYield02','ElasticIsotropic']:raise ValueError('not a soil continuum reference')
 source=Path(source);spec=importlib.util.spec_from_file_location('_capstone_ground_release',source);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 def registry(materials,order,analysis_nd=None):
  out=[]
  for key in order:
   meta=register_nd(ops,key,materials[key]['tag'],elastic=elastic,stiffness_factor=stiffness_factor)
   if meta['law']=='PressureDependMultiYield02':out.append(meta['tag'])
  return out
 def gravity(nodal_self_weight_loads,nonlinear_material_tags,ndm=3,n_steps=10,d_lambda=.1):
  if ndm!=3:raise ValueError('3D only')
  return initialize_gravity(ops,nodal_self_weight_loads,nonlinear_material_tags,n_steps=n_steps)
 mod.register_opensees_soil_materials=registry;mod.run_elastic_gravity_then_plastic_stage=gravity
 # Descriptive fields satisfy original mesh checks; actual law uses SI catalog.
 rec={**m,'tag':101,'model':'ElasticIsotropic' if elastic else m['law'],'rho':m['density_kg_m3'],'E':m.get('E_Pa'), 'refShearModul':m['G_Pa']*stiffness_factor,'refBulkModul':m['K_Pa']*stiffness_factor}
 if elastic:rec['E']=m['E_Pa']*stiffness_factor
 cfg=dict(body_name=body,gravity=catalog()['gravity_m_s2'][body],materials={name:rec},material_order=[name],assign_material=lambda **kwargs:name,base_depth=depth_m,n_depth=n_depth,vertical_mesh={'mode':'uniform'},boundary_condition_mode='gravity_static_box',material_confidence={name:m['status']})
 data=mod.build_opensees_soil_model_3d(X,Y,Z,cfg,run_gravity=True)
 data['release']=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),material_name=name,stiffness_factor=stiffness_factor,elastic=elastic,scope='gravity/static mesh only; DP active during gravity; PDMY stage0-to1 with re-equilibration; not dynamic SSI boundary validation')
 data['analysis_stage']='gravity_equilibrated_static_reference'
 return data
