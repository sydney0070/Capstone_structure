"""Runnable research draft: Mars raft / Moon annular shell + same scaled steel dome.
Default is a compact 4m-diameter benchmark, not the 40m design or a crater model.
No Martian accelerogram is silently assigned to the Moon. No pressure by default.
"""
from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PROJECT_ROOT), str(PROJECT_ROOT/'structure'),
               str(PROJECT_ROOT/'ground/ground_modeling'),
               str(PROJECT_ROOT/'material/material_release_20261101')]

import sys,argparse,json,time,math,hashlib,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'ground/ground_modeling')]
import openseespy.opensees as ops
from ground_hybrid import default_config
from ground_20261107 import (
    GroundHistory, build_ground_20261107, ground_state, insight_reference_site,
    static_ramp, validate_ground_20261107,
)
from SSI.interface import contact_summary
from ground_hybrid.run_journal import RunJournal,atomic_json
from planet_foundations_20261109 import FoundationConfig,add_foundation,connect_dome_port
from mars_dome_ring_ground_20261107 import axis_with_ports
import dome_in_mars_tuned as dome

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--body',choices=['Mars','Moon'],default='Mars');p.add_argument('--kind',choices=['raft','ring']);p.add_argument('--rotation',choices=['fixed','hinged','spring'],default='fixed');p.add_argument('--rotational-stiffness',type=float)
 p.add_argument('--radius',type=float,default=2.);p.add_argument('--width',type=float,default=1.2);p.add_argument('--thickness',type=float,default=.15);p.add_argument('--spacing',type=float,default=.75);p.add_argument('--soil-factor',type=float,default=1.);p.add_argument('--penalty-factor',type=float,default=10.);p.add_argument('--mass-ton',type=float,default=.2);p.add_argument('--domain-width',type=float,default=8.);p.add_argument('--depth',type=float,default=3.);p.add_argument('--vertical-spacing',type=float,default=1.);p.add_argument('--static-steps',type=int,default=40)
 p.add_argument('--worker',action='store_true',help=argparse.SUPPRESS);p.add_argument('--max-runtime',type=float,default=180.)
 p.add_argument('--mesh',choices=['cartesian','conforming'],default='cartesian');p.add_argument('--sectors',type=int,default=24);p.add_argument('--radial-subdivisions',type=int,default=1)
 p.add_argument('--plan-refinement',type=int,default=1);p.add_argument('--vertical-subdivisions',type=int,default=1)
 p.add_argument('--max-ground-elements',type=int,default=20000,help='explicit memory/work budget, not a convergence tolerance')
 p.add_argument('--ground-history',action='store_true')
 p.add_argument('--initial-load-increment',type=float)
 p.add_argument('--linear-system',choices=['UmfPack','SuperLU'],default='UmfPack')
 p.add_argument('--initial-tangent-seed',action='store_true',help='first accepted load step uses initial tangent, while all materials remain nonlinear')
 p.add_argument('--symmetric-gravity',action='store_true',help='research-only mirror-plane tangential constraints; reactions must be negligible; static only')
 order=p.add_mutually_exclusive_group();order.add_argument('--newton-first',action='store_true');order.add_argument('--krylov-first',action='store_true');order.add_argument('--initial-tangent-first',action='store_true')
 p.add_argument('--foundation-only',action='store_true');p.add_argument('--anchor-stiffness',type=float);p.add_argument('--anchor-capacity',type=float);p.add_argument('--output-dir',type=Path,required=True)
 a=p.parse_args()
 if any(not math.isfinite(v) or v<=0 for v in [a.radius,a.width,a.thickness,a.spacing,a.soil_factor,a.penalty_factor,a.mass_ton,a.domain_width,a.depth,a.vertical_spacing]):p.error('positive finite geometry/material parameters')
 if a.static_steps<1:p.error('positive static-steps required')
 if a.max_ground_elements<1:p.error('positive element budget required')
 if a.domain_width<=2*a.radius+a.width:p.error('domain must enclose foundation')
 if a.output_dir.exists() and any(a.output_dir.iterdir()):p.error('fresh output directory required')
 a.output_dir.mkdir(parents=True,exist_ok=True);start=time.perf_counter();result={'parameters':vars(a)|{'output_dir':str(a.output_dir)},'complete':False,'safe_verdict':False}
 try:
  result['runtime']={'python':sys.version,'opensees':ops.version(),'numpy':np.__version__}
  from ground_hybrid.process_control import THREAD_VARIABLES
  result['runtime']['thread_environment']={name:os.environ.get(name) for name in THREAD_VARIABLES}
  sources=[Path(__file__),ROOT/'structure/planet_foundations_20261109.py',ROOT/'structure/dome_in_mars_tuned.py',ROOT/'ground/ground_modeling/ground_20261107.py',ROOT/'ground/ground_modeling/ground_hybrid/mesh.py',ROOT/'ground/ground_modeling/ground_hybrid/interface.py',ROOT/'material/material_release_20261101/capmat/ssi.py']
  sources += [ROOT/'structure/mars_dome_ring_ground_20261107.py',ROOT/'ground/ground_modeling/ground_modeling.py']
  sources += list((ROOT/'SSI').glob('*.py'))+list((ROOT/'structure/sections').glob('*.py'))
  sources += [ROOT/'integrated/mars_ssi_integrated_20261109.py']
  sources += [ROOT/'material/materials_20261107.py']
  sources += list((ROOT/'ground/ground_modeling/ground_hybrid').rglob('*.py'))+list((ROOT/'ground/ground_modeling/ground_hybrid/configs').glob('*.json'))
  sources += list((ROOT/'material/material_release_20261101/capmat').rglob('*.py'))+[ROOT/'material/material_release_20261101/configs/catalog.json']
  result['code_sha256']={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(set(sources))}
  result['model_validation']={'pressure_closure':False,'foundation_material_calibrated_for_planet':False,'anchor_pullout_calibrated':False,'nonlinear_mesh_independence_established':False}
  ops.wipe();ops.logFile(str(a.output_dir/'opensees.log'));theta=np.arange(12)*2*np.pi/12;ports=np.c_[a.radius*np.cos(theta),a.radius*np.sin(theta)]
  axis=axis_with_ports(a.domain_width,ports.ravel(),a.spacing);X,Y=np.meshgrid(axis,axis)
  c=default_config(a.body);c.update(width_m=a.domain_width,length_m=a.domain_width,depth_m=a.depth,core_width_m=.7*a.domain_width,core_length_m=.7*a.domain_width,max_vertical_spacing_m=a.vertical_spacing,stiffness_factor=a.soil_factor,lenses=[],max_elements=a.max_ground_elements)
  kind=a.kind or ('raft' if a.body=='Mars' else 'ring')
  if a.mesh=='conforming':
   from ground_hybrid.conforming_mesh import build_conforming_ground
   g=validate_ground_20261107(build_conforming_ground(c,radius=a.radius,width=a.width,kind=kind,sectors=a.sectors,radial_subdivisions=a.radial_subdivisions,plan_refinement=a.plan_refinement,vertical_subdivisions=a.vertical_subdivisions))
  else:g=build_ground_20261107(c,(X,Y,np.zeros_like(X)))
  result['mesh']=g['mesh'];result['ground_elements']=len(g['elements']);result['ground_revision']=g['revision_validation']
  if a.body=='Mars':result['site'],samples=insight_reference_site(a.domain_width);np.savetxt(a.output_dir/'site_dem_samples.csv',samples,delimiter=',')
  else:result['site']={'kind':'synthetic_lunar_reference','coordinate':None,'crater_geometry_modeled':False,'reason':'attachment is conceptual lunar ring; no selected lunar seismic record or validated crater survey'}
  config=FoundationConfig(body=a.body,kind=a.kind or ('raft' if a.body=='Mars' else 'ring'),radius_m=a.radius,width_m=a.width,thickness_m=a.thickness,penalty_factor=a.penalty_factor,anchor_stiffness_N_m=a.anchor_stiffness,anchor_capacity_N=a.anchor_capacity)
  f=add_foundation(g,config,ports);loads=dict(g['gravity_loads']);loads.update(f['gravity_loads']);dome_mass={};connectors=[]
  symmetry=[]
  if a.symmetric_gravity:
   # Full domain, only mirror-plane tangential foundation displacements fixed.
   # No radial/vertical restraint. This is valid only for symmetric gravity;
   # post-solve reaction checks reject its use to prop an asymmetric solution.
   for n in f['nodes']:
    x,y,z=ops.nodeCoord(n);fix=[int(abs(x)<1e-9),int(abs(y)<1e-9),0,0,0,0]
    if any(fix):ops.fix(n,*fix);symmetry.extend((n,dof+1) for dof,value in enumerate(fix) if value)
   result['static_symmetry']={'constraints':symmetry,'scope':'static gravity mirror-plane tangential restraint, no seismic use; not a physical anchor'}
  if not a.foundation_only:
   scale=a.radius/20.;coords={n:np.array([v[1],v[0],v[2]])*scale for n,v in dome.dome_NODES.items()}
   ops.model('basic','-ndm',3,'-ndf',6);mapping={}
   for n,xyz in coords.items():
    t=800000+n;mapping[n]=t;ops.node(t,*xyz)
    if abs(xyz[2])<1e-9:
     matches=[fp for fp in f['dome_ports'] if np.linalg.norm(np.array(ops.nodeCoord(fp))-xyz)<1e-7]
     if len(matches)!=1:raise ValueError('exact dome/foundation port mismatch')
     connectors.append(connect_dome_port(matches[0],t,rotation=a.rotation,rotational_stiffness_Nm_rad=a.rotational_stiffness,tag=900000+n))
   section,_=dome.make_section_for_target_mass(a.mass_ton/scale);dome.register_tuned_shs(section)
   ops.geomTransf('PDelta',820000,0.,0.,1.);ops.geomTransf('PDelta',820001,1.,0.,0.)
   for e,i,j in dome.DOME_EDGES:
    tr=820001 if np.linalg.norm((coords[j]-coords[i])[:2])<1e-9 else 820000
    ops.element('forceBeamColumn',800000+e,mapping[i],mapping[j],tr,dome.STEEL_INTEGRATION_TAG,'-iter',40,1e-10)
   dome_mass=dome.compute_tributary_mass(coords,dome.DOME_EDGES,mapping,section['A'])
   steel_mass=sum(dome_mass.values())
   if not math.isclose(steel_mass,a.mass_ton*1000,rel_tol=1e-10):raise RuntimeError('scaled steel target mass mismatch')
   for n,m in dome.compute_membrane_node_mass(mapping,coords,dome.DOME_EDGES,verbose=False).items():dome_mass[n]=dome_mass.get(n,0.)+m
   for n,m in dome_mass.items():ops.mass(n,m,m,m,0.,0.,0.);loads[n]=[0.,0.,-m*c['gravity_m_s2'],0.,0.,0.]
   result['steel_mass_kg']=steel_mass
  np.savetxt(a.output_dir/'foundation_nodes.csv',[[n,*ops.nodeCoord(n),f['nodal_mass_kg'][n],f['weights_m2'][n]] for n in f['nodes']],delimiter=',',header='node,east_m,north_m,up_m,mass_kg,tributary_area_m2',comments='')
  np.savetxt(a.output_dir/'foundation_shells.csv',[[e,*ops.eleNodes(e)] for e in f['elements']],delimiter=',',header='element,n1,n2,n3,n4',comments='')
  result['ground_laws']={n:m['active_law'] for n,m in g['materials'].items()}
  result['foundation']={k:f[k] for k in ['config','area_m2','target_area_m2','area_relative_error','total_mass_kg','scope','anchors']};result['connectors']=connectors
  journal=RunJournal(a.output_dir,result);history=None
  if a.ground_history:
   history=GroundHistory(g,a.output_dir/'static_ground_history.csv');history.capture()
  algorithms=('NewtonInitial','KrylovNewton','NewtonLineSearch','Newton') if a.initial_tangent_first else ('KrylovNewton','NewtonLineSearch','Newton') if a.krylov_first else ('Newton','KrylovNewton','NewtonLineSearch') if a.newton_first else ('NewtonLineSearch','Newton','KrylovNewton')
  result['gravity']=static_ramp(g,loads,steps=a.static_steps,max_attempts=max(300,a.static_steps*3),algorithms=algorithms,attempt_observer=journal,on_commit=history.capture if history else None,initial_increment=a.initial_load_increment,linear_system=a.linear_system,first_step_algorithm='NewtonInitial' if a.initial_tangent_seed else None)
  if history:result['ground_history']=history.save(a.output_dir/'static_ground_history.csv')
  result['ground_state']=ground_state(g,gauss_points=8);result['contact']=contact_summary(f)
  if symmetry:
   ops.reactions();reactions=[float(ops.nodeReaction(n,dof)) for n,dof in symmetry]
   structural_weight=sum(-loads[n][2] for n in list(f['nodes'])+list(dome_mass))
   error=sum(abs(v) for v in reactions)/structural_weight
   result['static_symmetry'].update(reactions_N=reactions,total_absolute_reaction_to_structural_weight=error)
   if error>1e-6:raise RuntimeError('symmetry restraints carry non-negligible load; symmetry experiment rejected')
  u=np.array([ops.nodeDisp(n) for n in f['nodes']]);result['settlement_m']=dict(maximum=float(-u[:,2].min()),differential=float(np.ptp(u[:,2])))
  weights=np.array([f['weights_m2'][n] for n in f['nodes']]);port_u=np.array([ops.nodeDisp(n) for n in f['dome_ports']])
  if not np.isfinite(u).all() or not np.isfinite(port_u).all():raise RuntimeError('finite foundation response required')
  result['settlement_m'].update(area_weighted_mean=float(-weights@u[:,2]/sum(weights)),dome_ports_maximum=float(-port_u[:,2].min()),dome_ports_differential=float(np.ptp(port_u[:,2])))
  result['dome_port_displacement_m_rad']=[dict(node=n,coordinate_m=list(ops.nodeCoord(n)),displacement_m_rad=list(ops.nodeDisp(n))) for n in f['dome_ports']]
  expected=f['total_mass_kg']+sum(dome_mass.values());actual=sum(ops.nodeMass(n,1) for n in list(f['nodes'])+list(dome_mass))
  if not math.isclose(actual,expected,rel_tol=1e-10):raise RuntimeError('structural nodal mass check')
  result['mass_check_kg']=[actual,expected];result['complete']=True
 except BaseException as exc:result['error']=str(exc);raise
 finally:
  result['elapsed_s']=time.perf_counter()-start;atomic_json(a.output_dir/'summary.json',result);ops.wipe()
def supervised_cli():
 import os
 from ground_hybrid.process_control import MANAGED_ENV,run_isolated
 p=argparse.ArgumentParser(add_help=False);p.add_argument('--max-runtime',type=float,default=180.);p.add_argument('--output-dir',type=Path)
 a,_=p.parse_known_args()
 if '--worker' in sys.argv:return main()
 if not math.isfinite(a.max_runtime) or a.max_runtime<=0:raise SystemExit('positive finite max-runtime required')
 if a.output_dir and a.output_dir.exists() and any(a.output_dir.iterdir()):
  print('fresh output directory required',file=sys.stderr);return 2
 status=run_isolated([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--worker'],a.max_runtime,join_managed_group=os.environ.get(MANAGED_ENV)=='1')
 if a.output_dir:
  a.output_dir.mkdir(parents=True,exist_ok=True)
  tmp=a.output_dir/'run_status.json.tmp';tmp.write_text(json.dumps(status|{'safe_verdict':False},indent=2));tmp.replace(a.output_dir/'run_status.json')
 return status['returncode']
if __name__=='__main__':raise SystemExit(supervised_cli())
