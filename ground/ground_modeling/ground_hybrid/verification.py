"""Independent mechanics checks and assembled integration gates. Run via --verify."""
import contextlib
import importlib.util
import json
import math
import sys
import time
import traceback
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
from . import default_config,build_ground,add_raft,initialize_gravity,apply_static_loads
from .materials import record,register,RELEASE
from .interface import pressure_loads,raft_response
from .solver import initial_stress_summary
from .io import write_json,export_results

OUT=Path(__file__).resolve().parents[1]/'hybrid_verification'


def compact(body,*,factor=1.,elastic=False,spacing=1.,width=6.):
    c=default_config(body)
    c.update(width_m=width,length_m=width,core_spacing_m=spacing,max_outer_spacing_m=3.,
             max_vertical_spacing_m=2.,stiffness_factor=factor,elastic=elastic)
    return c


def analytical_column(body):
    c=compact(body,elastic=True);c['lenses']=[]
    X,Y=np.meshgrid([0.,1.],[0.,1.]);Z=np.zeros_like(X)
    ops.wipe();g=build_ground(c,(X,Y,Z));initialize_gravity(g)
    expected=0.;overburden=0.;gravity=c['gravity_m_s2']
    for e in sorted(g['elements'].values(),key=lambda e:e['depth_m']):
        m=g['materials'][e['material']];h=e['height_m'];M=m['K_Pa']+4*m['G_Pa']/3
        expected+=(overburden*h+m['density_kg_m3']*gravity*h*h/2)/M
        overburden+=m['density_kg_m3']*gravity*h
    top=g['index'][0,0,0];actual=-ops.nodeDisp(top,3)
    relative=abs(actual-expected)/expected
    return dict(body=body,analytical_settlement_m=expected,actual_settlement_m=actual,
                relative_error=relative,**{'pass':relative<1e-7})


def geometry_consistency():
    c=compact('Mars',elastic=True);c.update(depth_m=3.);c['lenses']=[]
    X,Y=np.meshgrid([-3.,-.5,1.,3.],[-3.,0.,3.]);Z=.1*X+.03*Y
    ops.wipe();g=build_ground(c,(X,Y,Z));weights=g['gravity_loads']
    total=np.sum(list(weights.values()),axis=0)
    nodalmoment=np.sum([np.cross(g['nodes'][t],v[:3]) for t,v in weights.items()],axis=0)
    # Independent geometric integration of z-directed uniform body forces: 
    # horizontal first moments for each vertical-prism hexahedron via 2x2 quad.
    expected=np.zeros(3);mass=0.;a=1/math.sqrt(3)
    signs=np.array([[-1,-1],[1,-1],[1,1],[-1,1]])
    for e in g['elements'].values():
        xyz=np.array([g['nodes'][t] for t in e['nodes']]);rho=g['materials'][e['material']]['density_kg_m3']
        dx=xyz[1,0]-xyz[0,0];dy=xyz[3,1]-xyz[0,1]
        for r in [-a,a]:
            for s in [-a,a]:
                N=np.prod(1+signs*np.array([r,s]),axis=1)/4
                xy=N@xyz[:4,:2];h=N@(xyz[4:,2]-xyz[:4,2]);dm=rho*dx*dy/4*h;mass+=dm
                expected+=np.cross([xy[0],xy[1],0],[0,0,-dm*c['gravity_m_s2']])
    moment_error=float(np.linalg.norm(expected-nodalmoment)/max(np.linalg.norm(expected),1.))
    mass_error=abs(mass-g['mass_kg'])/mass
    # Invalid terrain must be rejected, rather than silently filling voids.
    rejected=False
    try:
        bad=Z.copy();bad[0,0]=np.nan;build_ground(c,(X,Y,bad))
    except ValueError:rejected=True
    return dict(mass_relative_error=mass_error,moment_relative_error=moment_error,
                invalid_terrain_rejected=rejected,**{'pass':mass_error<1e-10 and moment_error<1e-10 and rejected})


def assembled(body,factor=1.,*,spacing=1.,width=6.,frame=False,penalty=10.,export=False):
    ops.wipe();g=build_ground(compact(body,factor=factor,spacing=spacing,width=width));raft=add_raft(g,penalty_factor=penalty)
    f=None
    if frame:
        from .cli import add_demo_frame
        f=add_demo_frame(raft,g)
    initialize_gravity(g,extra_loads=f['gravity_loads'] if f else None)
    gravity_stress=initial_stress_summary(g)
    loads=pressure_loads(raft,10000.)
    if f:
        for t in f['top_nodes']:loads[t]=[50.,0.,0.,0.,0.,0.]
    stat=apply_static_loads(g,loads);r=raft_response(g,raft);contact=r['contact']
    expected=raft['mass_kg']*g['config']['gravity_m_s2']+40000.
    if f:expected+=(f['member_mass_kg']+f['lumped_mass_kg'])*g['config']['gravity_m_s2']
    normal_error=abs(contact['total_normal_N']-expected)/expected
    friction_excess=max(max(0,row['shear_N']-row['friction_limit_N']) for row in contact['rows'])
    passed=(normal_error<1e-5 and contact['max_mpc_error_m']<1e-9 and friction_excess<.05
            and r['mean_additional_settlement_m']>0 and len(g['material_volume_fractions'])>=3)
    if export:export_results(g,OUT/f'{body}_assembled',raft=raft,extra=dict(static=stat,frame=f),html=False)
    return dict(body=body,factor=factor,spacing=spacing,width=width,frame=frame,penalty=penalty,
                soil_mass_kg=g['mass_kg'],materials=g['material_volume_fractions'],
                additional_settlement_m=r['mean_additional_settlement_m'],
                differential_settlement_m=r['differential_settlement_m'],
                force_error=stat['force_balance']['relative_error'],contact_normal_relative_error=normal_error,
                friction_capacity_excess_N=friction_excess,contact_mpc_error_m=contact['max_mpc_error_m'],
                max_contact_penetration_m=contact['max_penetration_m'],
                gravity_stress_samples=len(gravity_stress),**{'pass':passed})


def boundary_rest(body):
    from .dynamics import activate_absorbing_boundaries,run_transient
    ops.wipe();g=build_ground(compact(body));raft=add_raft(g);initialize_gravity(g)
    before={t:np.array(ops.nodeDisp(t)) for t in raft['nodes']}
    boundary=activate_absorbing_boundaries(g)
    r=run_transient(g,dt_s=.001,steps=25,record_nodes=raft['nodes'])
    drift=max(np.linalg.norm(np.array(row['nodes'][t]['u'])-before[t]) for row in r['rows'] for t in raft['nodes'])
    return dict(body=body,max_unforced_drift_m=drift,boundary=boundary,**{'pass':drift<1e-7})


def wave_column():
    from .dynamics import activate_absorbing_boundaries,add_incident_base_velocity,run_transient
    c=compact('Mars',elastic=True);c.update(depth_m=10.,max_vertical_spacing_m=.125)
    c['layers']=[dict(bottom_m=None,material='mars_jez1_peak')];c['lenses']=[]
    X,Y=np.meshgrid([0.,1.],[0.,1.]);Z=np.zeros_like(X);ops.wipe();g=build_ground(c,(X,Y,Z))
    # A horizontally periodic, elastic shear column. Y/Z translations fixed;
    # X tied across each plane. This benchmark isolates the base impedance.
    for t,fixed in g['fixities'].items():
        for d,yes in enumerate(fixed,1):
            if yes:ops.remove('sp',t,d)
        bottom=abs(g['nodes'][t][2]+10)<1e-9
        retained_base=bottom and t==g['index'][0,0,g['mesh']['nz']]
        ops.fix(t,int(retained_base),1,1);g['fixities'][t]=[int(retained_base),0,0]
        g['gravity_loads'][t]=np.zeros(3)
    for k in range(g['mesh']['nz']+1):
        retained=g['index'][0,0,k]
        for i,j in [(1,0),(1,1),(0,1)]:ops.equalDOF(retained,g['index'][i,j,k],1)
    g['faces']=[f for f in g['faces'] if f['side'] in ('bottom','top')]
    initialize_gravity(g);activate_absorbing_boundaries(g)
    dt=.00025;steps=2800;t=np.arange(steps+1)*dt;T=.08
    v=.001*np.sin(2*np.pi*t/T)*np.sin(np.pi*t/T)**2;v[t>T]=0
    add_incident_base_velocity(g,v,dt,axis=1)
    top=g['index'][0,0,0];r=run_transient(g,dt_s=dt,steps=steps,record_nodes=[top])
    m=record('mars_jez1_peak');Vs=math.sqrt(m['G_Pa']/m['density_kg_m3']);delay=10/Vs
    tt=np.array([row['time_s'] for row in r['rows']]);actual=np.array([row['nodes'][top]['v'][0] for row in r['rows']])
    expected=2*np.interp(tt-delay,t,v,left=0,right=0)
    error=np.linalg.norm(actual-expected)/np.linalg.norm(expected)
    late=float(max(abs(actual[tt>.5]))/max(abs(actual)))
    np.savetxt(OUT/'shear_wave_column.csv',np.column_stack([tt,actual,expected]),delimiter=',',header='time_s,computed_top_velocity_m_s,analytical_top_velocity_m_s',comments='')
    return dict(Vs_m_s=Vs,travel_time_s=delay,relative_L2_error=float(error),late_reflection_ratio=late,
                scope='normal-incidence homogeneous elastic S-wave; not layered oblique-wave/free-field validation',
                **{'pass':error<.08 and late<.03})


def material_checks():
    path=RELEASE/'verification'/'soil_material.py';spec=importlib.util.spec_from_file_location('_soil_tests',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.OUT=OUT/'material';mod.OUT.mkdir(exist_ok=True)
    mod.register_nd=register;mod.material=record
    rows=[]
    for name in ['earth_sand_dr30_dp_reference','earth_sand_dr50_dp_reference','earth_sand_dr75_dp_reference',
                 'moon_jsc1a_17_peak','moon_jsc1a_18_peak','moon_jsc1a_19_peak',
                 'mars_jez1_peak','mars_insight_duricrust_proxy',
                 'moon_field_proxy_17','moon_field_proxy_18','moon_field_proxy_19']:
        for p in [1000.,15000.,60000.]:rows.append(mod.triaxial(name,p))
    for dr in [30,50,75]:
        name=f'earth_sand_dr{dr}_dry_reference'
        mod.register_nd=lambda ops,ignored,tag,chosen=name:register(ops,chosen,tag)
        row=mod.pdmy();row['name']=name;rows.append(row)
        source=mod.OUT/'earth_sand_dr50_dry_reference_cyclic.csv'
        source.rename(mod.OUT/f'{name}_checked_cyclic.csv')
    write_json(OUT/'material_summary.json',rows)
    return dict(cases=rows,**{'pass':all(r['pass'] for r in rows)})


def public_api_guards():
    ops.wipe();ops.model('basic','-ndm',3,'-ndf',6);ops.node(9,100.,100.,100.);ops.fix(9,1,1,1,1,1,1)
    g=build_ground(compact('Moon',elastic=True));preserved=9 in ops.getNodeTags()
    before=len(ops.getNodeTags());collision=False
    try:build_ground(compact('Moon',elastic=True))
    except ValueError:collision=True
    intact=len(ops.getNodeTags())==before
    initialize_gravity(g);staged=False
    try:add_raft(g)
    except ValueError:staged=True
    return dict(existing_structure_preserved=preserved,collision_rejected=collision,domain_intact=intact,
                unsupported_staged_installation_rejected=staged,**{'pass':preserved and collision and intact and staged})


def structure_port_transfer():
    from .interface import connect_matching_nodes
    ops.wipe();g=build_ground(compact('Earth',elastic=True));raft=add_raft(g)
    port=raft['nodes'][0];external=490000
    ops.model('basic','-ndm',3,'-ndf',6);ops.node(external,*ops.nodeCoord(port))
    connect_matching_nodes(raft,((a,b) for a,b in [(port,external)]))
    initialize_gravity(g)
    result=apply_static_loads(g,{external:[100.,50.,-100.,0.,100.,0.]})
    difference=float(np.max(np.abs(np.asarray(ops.nodeDisp(port))-ops.nodeDisp(external))))
    rotation=float(np.linalg.norm(ops.nodeDisp(port)[3:]))
    return dict(max_six_dof_difference=difference,raft_rotation_norm_rad=rotation,
                force_balance=result['force_balance'],**{'pass':difference<1e-10 and rotation>1e-12})


def main():
    OUT.mkdir(exist_ok=True);start=time.time();results=[]
    write_json(OUT/'summary.json',dict(engine=ops.version(),complete=False,cases=[],all_pass=False))
    tasks=[('geometry_consistent_gravity',geometry_consistency),('API_guards',public_api_guards),
           ('structure_port_transfer',structure_port_transfer),('material_checks',material_checks)]
    for b in ['Earth','Moon','Mars']:
        tasks.extend([(b+'_analytical_layered_column',lambda body=b:analytical_column(body)),
                      (b+'_nonlinear_raft_frame',lambda body=b:assembled(body,frame=True,export=True)),
                      (b+'_absorbing_rest',lambda body=b:boundary_rest(body))])
    tasks.append(('incident_shear_wave',wave_column))
    for name,fn in tasks:
        t=time.time()
        with (OUT/(name+'.log')).open('w') as log,contextlib.redirect_stdout(log):
            try:result=fn()
            except Exception as e:result=dict(error=str(e),traceback=traceback.format_exc(),**{'pass':False})
        results.append(dict(name=name,seconds=time.time()-t,**result))
        write_json(OUT/'summary.json',dict(engine=ops.version(),elapsed_s=time.time()-start,cases=results,all_pass=all(r['pass'] for r in results),complete=len(results)==len(tasks)))
        print(name,'PASS' if result['pass'] else 'FAIL',round(time.time()-t,2),result.get('error',''),flush=True)
    return 0 if all(r['pass'] for r in results) else 1

if __name__=='__main__':raise SystemExit(main())
