"""20261109 integrated Mars ground–RC foundation–dome analysis.

InSight/ELYSE-centred flat reference site by default, ENU, SI. One domain, one mass/selfweight assignment and one
3-component earthquake input. Default: unpressurized, uncovered coupling check.
Not an optimal/safe design; no claim of reconstructed InSight site geology.
"""
from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PROJECT_ROOT), str(PROJECT_ROOT/'structure'),
               str(PROJECT_ROOT/'ground/ground_modeling'),
               str(PROJECT_ROOT/'material/material_release_20261101')]

import argparse,json,math,sys,time,hashlib
import numpy as np
import openseespy.opensees as ops

HERE=Path(__file__).resolve().parent
GROUND=HERE.parent/'ground'/'ground_modeling'
sys.path.insert(0,str(GROUND))
import dome_in_mars_tuned as dome
from ring_foundation_rc import RingConfig,build_ring
from ground_hybrid import default_config
from insight_motion import INSIGHT_ACC_PATH
from ground_20261107 import build_ground_20261107 as build_ground
from ground_20261107 import (
    GroundHistory, advance, ground_state, insight_reference_site, lateral_dashpots,
    prepare_motion, remove_base_rigid_motion, require_finite, resolution,
    save_response_viewer, static_ramp,
)
sys.path.insert(0,str(HERE.parent/"material"))
from materials_20261107 import audit_catalogs



def axis_with_ports(width,ports,max_spacing):
    seeds=sorted(set([round(float(v),12) for v in ports]+[-width/2,0.,width/2]))
    if seeds[0]<-width/2 or seeds[-1]>width/2:raise ValueError('Ground must contain ring')
    values=[seeds[0]]
    for a,b in zip(seeds[:-1],seeds[1:]):
        values.extend(np.linspace(a,b,max(1,math.ceil((b-a)/max_spacing))+1)[1:])
    return np.array(values)


def assemble(args, *, stage="full"):
    """Caller owns wipe. Reuse RC ring ports as dome base nodes; no MP chains."""
    if stage not in ('ground','ring','full'):raise ValueError('Unknown assembly stage')
    theta=np.arange(args.ring_segments)*2*np.pi/args.ring_segments
    axis=axis_with_ports(args.ground_width,np.r_[20*np.cos(theta),20*np.sin(theta)],args.spacing)
    X,Y=np.meshgrid(axis,axis);Z=np.zeros_like(X)
    c=default_config('Mars')
    c["stiffness_factor"]=args.ground_stiffness_factor
    c.update(width_m=args.ground_width,length_m=args.ground_width,depth_m=getattr(args,'depth',c['depth_m']),core_width_m=40.,core_length_m=40.,
             core_spacing_m=args.spacing,max_outer_spacing_m=args.spacing,max_vertical_spacing_m=4.,
             tag_start=10000,max_elements=10000,elastic=args.elastic_ground)
    g=build_ground(c,(X,Y,Z))
    if stage=='ground':return dict(ground=g,extra_gravity={})
    r=build_ring(RingConfig(segments=args.ring_segments),demo_supports=False)
    surface=[t for (i,j,k),t in g['index'].items() if k==0]
    contacts=[];pairs=[]
    if args.interface=='contact':
        from SSI.contact import register_contact
        from SSI.parameters import contact_parameters
    for i,n in enumerate(r['nodes']):
        xyz=np.array(ops.nodeCoord(n));dist=[np.linalg.norm(g['nodes'][t]-xyz) for t in surface]
        j=int(np.argmin(dist));gn=surface[j]
        if dist[j]>1e-6:raise ValueError('No coincident ground node; no nearest-node snapping allowed')
        pairs.append((gn,n))
        if args.interface=='tied':ops.equalDOF(gn,n,1,2,3)
        else:
            m=g['materials']['mars_insight_duricrust_proxy']
            area=r['config']['width_m']*r['chord_m']
            p=contact_parameters(area,m['E_Pa'],m['G_Pa'],.2,10000.,.001,.35,penalty_factor=10.)
            contacts.append(register_contact(ops,700000+i,n,gn,p,
                auxiliary_nodes=(710000+2*i,710001+2*i),normal_material_tag=720000+2*i,
                rotation_material_tag=720001+2*i,friction_tag=730000+i,restore_ndf=6))
    if stage=='ring':
        return dict(ground=g,ring=r,contacts=contacts,ground_ring_pairs=pairs,extra_gravity={t:[0.,0.,-mass*g['config']['gravity_m_s2'],0.,0.,0.] for t,mass in r['nodal_mass_kg'].items()})
    ops.model('basic','-ndm',3,'-ndf',6)
    # Original dome uses North/East/Up. Convert coordinates, vectors and bases to ENU.
    coords={n:np.array([v[1],v[0],v[2]]) for n,v in dome.dome_NODES.items()}
    node_map={};ring_xyz={n:np.array(ops.nodeCoord(n)) for n in r['dome_ports']}
    for n,xyz in coords.items():
        if abs(xyz[2])<1e-8:
            match=[t for t,v in ring_xyz.items() if np.linalg.norm(v-xyz)<1e-6]
            if len(match)!=1:raise ValueError('Dome base/RC port mismatch')
            node_map[n]=match[0]
        else:
            node_map[n]=800000+n;ops.node(node_map[n],*xyz)
    section,_=dome.make_section_for_target_mass(args.mass_ton)
    dome.register_tuned_shs(section)
    ops.geomTransf('PDelta',820000,0.,0.,1.)
    ops.geomTransf('PDelta',820001,1.,0.,0.)
    dome_eles=[]
    for e,a,b in dome.DOME_EDGES:
        tr=820001 if np.linalg.norm((coords[b]-coords[a])[:2])<1e-9 else 820000
        tag=800000+e;ops.element('forceBeamColumn',tag,node_map[a],node_map[b],tr,dome.STEEL_INTEGRATION_TAG,'-iter',40,1e-10)
        dome_eles.append(tag)
    masses=dome.compute_tributary_mass(dome.dome_NODES,dome.DOME_EDGES,node_map,section['A'])
    membrane=dome.compute_membrane_node_mass(node_map,dome.dome_NODES,dome.DOME_EDGES,verbose=False)
    cover=dome.compute_cover_node_mass(node_map,dome.dome_NODES,dome.DOME_EDGES,args.cover_thickness,verbose=False)
    for part in (membrane,cover):
        for t,m in part.items():masses[t]=masses.get(t,0.)+m
    # ops.mass replaces existing mass: include ring mass at reused base ports exactly once.
    for t,m in masses.items():
        total=m+r['nodal_mass_kg'].get(t,0.)
        ops.mass(t,total,total,total,0.,0.,0.)
    extra={t:[0.,0.,-m*g['config']['gravity_m_s2'],0.,0.,0.] for t,m in r['nodal_mass_kg'].items()}
    for t,m in masses.items():
        if t not in extra:extra[t]=[0.]*6
        extra[t][2]-=m*g['config']['gravity_m_s2']
    return dict(ground=g,ring=r,node_map=node_map,dome_elements=dome_eles,dome_coords=coords,
                ground_ring_pairs=pairs,contacts=contacts,extra_gravity=extra,
                structure_mass_kg=sum(masses.values()),cover_mass_kg=sum(cover.values()),section=section)


def write_view(m,visual_tags,frames,path,scale,*,reference,metadata):
    """Same material-panel viewer as the ground-only path; real committed data."""
    g,r=m['ground'],m['ring']
    coordinates=[ops.nodeCoord(t) for t in visual_tags]
    displacement=[reference+np.asarray(f['du']) for f in frames]
    structures=[dict(name='RC 링 · C30 / SD400',color='#437e9e',
        segments=[(r['nodes'][i],r['nodes'][(i+1)%len(r['nodes'])]) for i in range(len(r['nodes']))]),
        dict(name='강재 돔 · S355',color='#147d86',
        segments=[(m['node_map'][a],m['node_map'][b]) for _,a,b in dome.DOME_EDGES])]
    site=metadata.get('site',{})
    location=(f"ELYSE {site['latitude_deg']}°N, {site['longitude_deg']}°E" if site.get('latitude_deg') is not None else '위치 미지정 합성 지반')
    context=dict(title='화성 통합 응답 · 지반–RC 링–돔',subtitle='비선형 지반 + 기초 + 구조물 통합본',
        summary=f"지반 {len(g['elements']):,} 요소 · RC 링 · S355 돔",static_only=len(frames)==1,
        scope=location+' 기준 평탄 지반. 문헌 기반 기준 층서. 실제 DEM 요철 미적용. '+
            '관측 지표 가속도의 기저 입력 가정 및 근사 경계. 링은 점 접촉 모델. ETFE는 질량만 반영. 색상은 재료 구분이며 파괴 판정이 아닙니다.')
    return save_response_viewer(g,visual_tags,coordinates,displacement,[f['t'] for f in frames],path,
        context=context,structures=structures,scale=scale)



def material_state(m):
    """Outer-step section diagnostics; not a failure/energy criterion."""
    E=dome.TUNED_STEEL['params']['E'];smax=ep=rc_c=rc_t=rebar=0.
    for e in m['dome_elements']:
        for ip in (1,2,3):
            x=np.asarray(ops.eleResponse(e,'section',ip,'fiberData'),float)
            if not len(x) or len(x)%5:raise RuntimeError('Steel fiber output missing')
            x=require_finite(x,f'fiber element={e} section={ip} time={ops.getTime()}').reshape(-1,5);smax=max(smax,float(max(abs(x[:,3]))));ep=max(ep,float(max(abs(x[:,4]-x[:,3]/E))))
    nb=len(m['ring']['bars_y_z_d_m'])
    for e in m['ring']['elements']:
        for ip in (1,2,5,6):
            x=np.asarray(ops.eleResponse(e,'section',ip,'fiberData'),float)
            if not len(x):continue  # Elastic interior has no fibers.
            if len(x)%5:raise RuntimeError('Unexpected RC fiber output')
            x=require_finite(x,f'fiber element={e} section={ip} time={ops.getTime()}').reshape(-1,5)
            if len(x)<=nb:raise RuntimeError('RC constituent fiber output incomplete')
            rc_c=max(rc_c,float(max(-x[:-nb,3])));rc_t=max(rc_t,float(max(x[:-nb,3])));rebar=max(rebar,float(max(abs(x[-nb:,3]))))
    return dict(steel_abs_stress_Pa=smax,steel_apparent_nonelastic_strain=ep,
        concrete_compression_Pa=max(rc_c,0.),concrete_tension_Pa=max(rc_t,0.),rebar_abs_stress_Pa=rebar)


def main(argv=None):
    from datetime import datetime
    import os,platform,uuid
    p=argparse.ArgumentParser(description=__doc__)
    for name,default in [('ground-width',60.),('spacing',5.),('depth',12.),('mass-ton',200.),('cover-thickness',0.),('net-pressure',0.),('duration',2.),('warmup',5.),('motion-scale',1.),('visual-scale',1e6),('freqmax',1.),('dt',.05),('ground-stiffness-factor',1.)]:
        p.add_argument('--'+name,type=float,default=default)
    p.add_argument('--max-runtime',type=float,default=180.,help='Wall-clock deadline in seconds; worker is terminated on timeout')
    p.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--quick',action='store_true',help='1 s warmup + 1 s target for connection checks, not history-converged design results')
    p.add_argument('--ring-segments',type=int,default=12)
    p.add_argument('--max-frames',type=int,default=41)
    p.add_argument('--interface',choices=['tied','contact'],default='contact')
    p.add_argument('--boundary',choices=['fixed-box','viscous-sides'],default='viscous-sides')
    p.add_argument('--target-pga',type=float,default=None,help='Optional target ENU vector PGA [m/s2], scaled virtual scenario, not design hazard')
    p.add_argument('--motion-start',type=float,default=None,help='Target interval source start [s], not warmup start')
    p.add_argument('--ground-history',action='store_true',help='All 8 GP work/plastic-strain diagnostics at every committed static and dynamic substep')
    p.add_argument('--static-steps',type=int,default=20)
    p.add_argument('--site',choices=['insight','synthetic-flat'],default='insight')
    p.add_argument('--allow-short-warmup',action='store_true',help='Explicitly permit insufficient prehistory, recorded in metadata')
    p.add_argument('--static-only',action='store_true');p.add_argument('--elastic-ground',action='store_true');p.add_argument('--no-html',action='store_true')
    p.add_argument('--allow-unclosed-pressure',action='store_true',help='Explicit research-only open dome pressure experiment; floor/closure missing')
    p.add_argument('--output-dir',type=Path,default=None)
    a=p.parse_args(argv)
    raw=sys.argv[1:] if argv is None else argv
    if a.quick:
        if '--warmup' not in raw:a.warmup=1.
        if '--duration' not in raw:a.duration=1.
    nums=[v for v in vars(a).values() if isinstance(v,float)]
    if not np.isfinite(nums).all() or a.ground_width<=42 or min(a.depth,a.spacing,a.mass_ton,a.duration,a.dt,a.freqmax,a.visual_scale,a.ground_stiffness_factor)<=0 or min(a.motion_scale,a.warmup,a.cover_thickness)<0:p.error('Invalid finite geometry/time/scale inputs')
    if not .05<a.freqmax<5 or a.max_frames<2 or a.static_steps<1 or a.ring_segments<12 or a.ring_segments%12:p.error('Require .05<freqmax<5, max-frames>=2, ring segments multiple of 12')
    if a.dt*a.freqmax>.05+1e-12:p.error('Use dt<=1/(20*freqmax); temporal resolution rejected before assembly')
    if a.target_pga is not None and (a.target_pga<0 or a.motion_scale!=1):p.error('target-pga >=0 and mutually exclusive with non-default motion-scale')
    if a.net_pressure and not a.allow_unclosed_pressure:p.error('Pressure closure/floor not modeled: use --allow-unclosed-pressure only for an explicit research scenario')
    if a.cover_thickness:p.error('Cover distribution/retention model not verified for coupled ground: dated structural baseline requires zero cover')
    out=a.output_dir.resolve() if a.output_dir else HERE/'results'/'mars_ssi_20261109'/('run_'+datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6])
    if out.exists() and any(out.iterdir()):p.error('Output directory is not empty: choose a fresh directory to preserve run provenance')
    out.mkdir(parents=True,exist_ok=True)
    status=dict(status='running',complete=False,safe_verdict=False,run_id=out.name)
    result=dict(date='2026-11-07',parameters={**vars(a),'output_dir':str(out)},runtime=dict(python=sys.version,opensees=ops.version(),platform=platform.platform()),
        assumptions=['ELYSE-centred flat DEM reference by default; relief idealized, no reconstructed geology or excavation history',
        'DP reference stiffness/strength, transferred S355 SHS model and nominal C30/SD400; no new experimental calibration',
        'Point contact along ring centerline, not resolved footing width/bearing pressure',
        'Local viscous sides are approximate, base remains fixed; no free-field boundary or deconvolved base motion',
        'Windowed strong-motion experiment, not full-record maximum displacement or design earthquake',
        'PDelta geometric transformation is not postbuckling/local buckling validation',
        'No pressure-vessel safety, optimal mass, RC shear/bond/delamination/fracture conclusions',
        'Dynamic response/material/contact extrema sampled at every committed adaptive substep'])
    started=time.perf_counter()
    def atomic(path,data):
        tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist() if isinstance(x,np.ndarray) else str(x) if isinstance(x,Path) else (_ for _ in ()).throw(TypeError(type(x).__name__))));tmp.replace(path)
    def save():
        try:atomic(out/'summary.json',result)
        except Exception as exc:
            status.update(status='failed',complete=False,serialization_error=str(exc))
            atomic(out/'run_status.json',status);raise
        atomic(out/'run_status.json',status)
    save();ops.wipe();ops.logFile(str(out/'opensees.log'))
    try:
        result['revision']='2026-11-09 fixes and ELYSE geographic reference'
        result['material_audit']=audit_catalogs()
        sources=[Path(__file__),GROUND/'ground_20261107.py',HERE.parent/'material'/'materials_20261107.py',PROJECT_ROOT/'structure/ring_foundation_rc.py',PROJECT_ROOT/'structure/dome_in_mars_tuned.py',GROUND/'ground_hybrid'/'mesh.py']
        sources += [GROUND/'insight_dynamic_viewer.py',GROUND/'insight_motion.py',GROUND/'insight_site_dem.py']
        sources += list((PROJECT_ROOT/'SSI').glob('*.py'))
        sources += list((PROJECT_ROOT/'structure/sections').glob('*.py'))
        sources += list((GROUND/'ground_hybrid').rglob('*.py'))
        sources += list((GROUND/'ground_hybrid'/'configs').glob('*.json'))
        sources += list((HERE.parent/'material'/'material_release_20261101'/'capmat').rglob('*.py'))
        sources += [HERE.parent/'material'/'material_release_20261101'/'configs'/'catalog.json']
        result['code_sha256']={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(set(sources))}
        from importlib.metadata import version
        result['runtime']['packages']={n:version(n) for n in ('numpy','scipy','obspy','plotly')}
        result['damping_assumption']='No imposed structural Rayleigh damping; nonlinear constitutive/contact effects and viscous sides only'
        if a.site=='insight':
            result['site'],samples=insight_reference_site(a.ground_width)
            np.savetxt(out/'site_dem_samples.csv',samples,delimiter=',',header='east_m,north_m,dem_elevation_m,dem_relative_to_station_m,modeled_surface_local_z_m',comments='')
        else:
            result['site']=dict(kind='synthetic_flat',geographic_origin=None)
        save()
        m=assemble(a);g,r=m['ground'],m['ring']
        result['resolution']=resolution(g,a.freqmax,a.dt)
        active_materials={name:g['materials'][name] for name in g['material_volume_fractions']}
        result['ground_laws']={name:m['active_law'] for name,m in active_materials.items()}
        if not a.elastic_ground and not any(v=='DruckerPrager' for v in result['ground_laws'].values()):
            raise RuntimeError('Nonlinear ground required; no active DruckerPrager soil')
        # Assert no old standalone supports and no duplicated nodal structural mass.
        if r['fixed'] or r['support_elements']:raise RuntimeError('Legacy ring springs present in assembled 3D ground')
        structural=set(r['nodes'])|set(m['node_map'].values())
        actual=sum(ops.nodeMass(n,1) for n in structural);expected=r['total_mass_kg']+m['structure_mass_kg']
        if not math.isclose(actual,expected,rel_tol=1e-12):raise RuntimeError('Duplicated/missing structural mass')
        result['assembly']=dict(ground_elements=len(g['elements']),ground_nodes=len(g['nodes']),ring_elements=len(r['elements']),steel_elements=len(m['dome_elements']),
            mass_check_kg=[actual,expected],ground_revision=g['revision_validation'],mesh=g['mesh'],ground_materials=active_materials,
            ground_layers=g['config']['layers'],material_volume_fractions=g['material_volume_fractions'],layer_assumption=g['config']['assumption'],
            ring_config=r['config'],ring_section=r['section'],ring_bars=r['bars_y_z_d_m'],steel_section=m['section'])
        loads={t:np.asarray(v) for t,v in g['gravity_loads'].items()}
        for t,v in m['extra_gravity'].items():
            if t in loads:raise RuntimeError('Unexpected structural/ground load tag collision')
            loads[t]=np.array(v)
        print('[assembly]',len(g['elements']),'ground bricks; output',out,flush=True)
        static_history=None
        if a.ground_history:
            static_history=GroundHistory(g,out/'static_ground_history.csv');static_history.capture()
        from ground_hybrid.run_journal import RunJournal
        save()
        journal=RunJournal(out,dict(parameters=result['parameters'],runtime=result['runtime'],code_sha256=result['code_sha256']))
        result['gravity']=static_ramp(g,loads,steps=a.static_steps,max_attempts=max(150,a.static_steps*3),on_commit=static_history.capture if static_history else None,attempt_observer=journal)
        if static_history:result['static_ground_history']=static_history.save(out/'static_ground_history.csv')
        if a.net_pressure:
            f=dome.pressure_node_forces(dome.dome_NODES,dome.DOME_EDGES,a.net_pressure)
            loads={m['node_map'][n]:[float(v[1]),float(v[0]),float(v[2]),0.,0.,0.] for n,v in f.items()}
            result['pressure']=static_ramp(g,loads,pattern=800002)
        if a.interface=='tied' and result['gravity']['physical_support_moment_error']>1e-5:raise RuntimeError('Physical support moment equilibrium failed')
        result['static_material_state']=material_state(m)
        result['static_ground_state']=ground_state(g)
        from SSI.interface import contact_summary
        if m['contacts']:result['static_contact']=contact_summary({'contacts':m['contacts']})
        else:
            err=max(np.linalg.norm(np.array(ops.nodeDisp(x)[:3])-ops.nodeDisp(y)[:3]) for x,y in m['ground_ring_pairs'])
            if err>1e-8:raise RuntimeError('Tied interface mismatch')
        tags=sorted(set(g['nodes'])|set(r['nodes'])|set(m['node_map'].values()))
        if a.max_frames*len(tags)*3*8>200_000_000:p.error('Frame memory estimate >200 MB; reduce --max-frames')
        ix={n:i for i,n in enumerate(tags)};xyz=np.array([ops.nodeCoord(n) for n in tags]);U0=require_finite([ops.nodeDisp(n)[:3] for n in tags],'static displacement')
        base_idx=[ix[n] for n in r['dome_ports']];dome_idx=[ix[n] for n in m['node_map'].values()];ground_idx=[ix[n] for n in g['nodes']];ring_idx=[ix[n] for n in r['nodes']]
        B0=xyz[base_idx]+U0[base_idx];D0=xyz[dome_idx]+U0[dome_idx]
        frames=[dict(t=0.,du=np.zeros_like(U0).tolist())]
        result['static_displacement_m']=dict(ground=float(max(np.linalg.norm(U0[ground_idx],axis=1))),ring=float(max(np.linalg.norm(U0[ring_idx],axis=1))),dome=float(max(np.linalg.norm(U0[dome_idx],axis=1))))
        result['boundary']=dict(kind='static_box')
        if not a.static_only:
            if a.boundary=='viscous-sides':result['boundary']=lateral_dashpots(g)
            ops.loadConst('-time',0.)
            input_dt,acc,meta=prepare_motion(INSIGHT_ACC_PATH,out,duration=a.duration,warmup=a.warmup,start=a.motion_start,freqmax=a.freqmax,allow_short_warmup=a.allow_short_warmup)
            if a.site=='insight' and any(c.split('.')[1]!='ELYSE' for c in meta['channels']):
                raise ValueError('Motion station does not match ELYSE site reference')
            meta['site_anchor']=dict(station='ELYSE',site_mode=a.site,location_matches_station_reference=a.site=='insight',surface_record_is_not_verified_base_input=True)
            target=meta['target_start_model_s']
            target_index=round(target/input_dt);raw_pga=float(max(np.linalg.norm(acc[target_index:],axis=1)))
            if a.target_pga is not None and raw_pga<=0:raise ValueError('Cannot normalize zero input PGA')
            scale=a.motion_scale if a.target_pga is None else a.target_pga/raw_pga
            meta.update(effective_scale=scale,target_vector_pga_m_s2=scale*raw_pga,pga_mars_g=scale*raw_pga/g['config']['gravity_m_s2'],scenario='observed_waveform' if scale==1 else 'zero_input_control' if scale==0 else 'scaled_virtual_waveform_NOT_hazard_design_event')
            total=(len(acc)-1)*input_dt;nsteps=int(round(total/a.dt))
            if not math.isclose(nsteps*a.dt,total,abs_tol=1e-8):raise ValueError('--dt must exactly divide selected record length')
            for i in range(3):
                t=900000+i;ops.timeSeries('Path',t,'-dt',input_dt,'-values',*acc[:,i].tolist(),'-factor',scale)
                ops.pattern('UniformExcitation',t,i+1,'-accel',t)
            ops.wipeAnalysis();ops.constraints('Plain');ops.numberer('RCM');ops.system('UmfPack')
            tolerance=max(.001,sum(np.linalg.norm(v[:3]) for v in g['applied_loads'].values())*1e-10)
            ops.test('NormUnbalance',tolerance,60);ops.algorithm('Newton');ops.integrator('Newmark',.5,.25);ops.analysis('Transient')
            dynamic_history=None
            if a.ground_history:
                dynamic_history=GroundHistory(g,out/'dynamic_ground_history.csv');dynamic_history.capture()
            rows=[];attempts=[];peaks=None;steel_peak={};contact_peaks=dict(penetration_m=0.,relative_shear_motion_m=0.,mpc_error_m=0.);apex=m['node_map'][6]
            frame_steps=set(np.linspace(0,nsteps,min(nsteps+1,a.max_frames),dtype=int))
            import csv
            columns=['model_time_s','target_time_s','in_target','ground_max_increment_m','ring_max_increment_m','dome_max_increment_m','dome_after_base_rigid_fit_m','apex_absolute_acceleration_m_s2']
            committed_samples=[0];cache={}
            with (out/'coupled_response.csv').open('w',newline='') as f,(out/'committed_response.csv').open('w',newline='') as cf:
                w=csv.writer(f);w.writerow(columns);cw=csv.writer(cf);cw.writerow(columns)
                def observe_commit():
                    nonlocal peaks
                    if dynamic_history:dynamic_history.capture()
                    t=ops.getTime();U=require_finite([ops.nodeDisp(n)[:3] for n in tags],f'displacement time={t}');du=U-U0
                    residual=remove_base_rigid_motion(B0,xyz[base_idx]+U[base_idx],D0,xyz[dome_idx]+U[dome_idx])
                    input_a=np.array([np.interp(t,np.arange(len(acc))*input_dt,acc[:,i]) for i in range(3)])*scale
                    absolute=require_finite(ops.nodeAccel(apex)[:3],f'apex acceleration time={t}')+input_a
                    values=[float(max(np.linalg.norm(du[indices],axis=1))) for indices in [ground_idx,ring_idx,dome_idx]]+[float(max(np.linalg.norm(residual,axis=1))),float(np.linalg.norm(absolute))]
                    in_target=t>=target-1e-9
                    row=[t,t-target,int(in_target),*values];cw.writerow(row);cf.flush();committed_samples[0]+=1
                    if in_target:
                        peaks=np.array(values) if peaks is None else np.maximum(peaks,values)
                        for k,v in material_state(m).items():steel_peak[k]=max(steel_peak.get(k,0.),v)
                        if m['contacts']:
                            c=contact_summary({'contacts':m['contacts']});contact_peaks['penetration_m']=max(contact_peaks['penetration_m'],c['max_penetration_m']);contact_peaks['mpc_error_m']=max(contact_peaks['mpc_error_m'],c['max_mpc_error_m'])
                            contact_peaks['relative_shear_motion_m']=max(contact_peaks['relative_shear_motion_m'],max(np.linalg.norm(v['relative_displacement_m'][:2]) for v in c['rows']))
                    cache.update(t=t,du=du,values=values,row=row,in_target=in_target)
                    atomic(out/'dynamic_checkpoint.json',dict(complete=False,safe_verdict=False,committed_samples=committed_samples[0],last_time_s=t,target_peaks=peaks.tolist() if peaks is not None else None,material_state_target_max=steel_peak,contact_target_extrema=contact_peaks))
                for step in range(nsteps):
                    advance(a.dt,attempts=attempts,on_commit=observe_commit)
                    t,du,values,row,in_target=(cache[k] for k in ('t','du','values','row','in_target'))
                    w.writerow(row);f.flush();rows.append(row)
                    if not a.no_html and step+1 in frame_steps:frames.append(dict(t=t,du=du.tolist()))
                    if (step+1)%max(1,nsteps//10)==0:print(f'[dynamic] {step+1}/{nsteps}; target={in_target}, ground={values[0]:.3e}m',flush=True)
            if dynamic_history:result['dynamic_ground_history']=dynamic_history.save(out/'dynamic_ground_history.csv')
            result['final_ground_state']=ground_state(g)
            result['dynamic']=dict(completed_steps=nsteps,committed_samples=committed_samples[0],extrema_sampling='every committed adaptive substep',integration_dt_s=a.dt,input_dt_s=input_dt,substep_attempts=attempts,
                target_peaks=peaks.tolist(),peak_order=['ground_m','ring_m','dome_m','dome_base_rigid_fit_removed_m','apex_absolute_acceleration_m_s2'],
                material_state_target_max=steel_peak,contact_target_extrema=contact_peaks,motion=meta,
                input_pga_vector_m_s2=float(scale*max(np.linalg.norm(acc,axis=1))),scale=scale,
                response='UniformExcitation translations are relative to base; apex absolute acceleration adds imposed ENU input')
        if not a.no_html:result['viewer']=write_view(m,tags,frames,out/'mars_integrated_response.html',a.visual_scale,reference=U0,metadata=result)
        result['elapsed_s']=time.perf_counter()-started;status.update(status='complete',complete=True,elapsed_s=result['elapsed_s'])
        print('[complete]',out,'elapsed',result['elapsed_s'],flush=True)
    except BaseException as exc:
        status.update(status='failed',complete=False,error=str(exc),elapsed_s=time.perf_counter()-started)
        raise
    finally:save();ops.wipe()

def supervised_cli():
    """Separate process enforces a deadline even inside a long native solver call."""
    import os,uuid
    from ground_hybrid.process_control import MANAGED_ENV,run_isolated
    from datetime import datetime
    raw=sys.argv[1:]
    if '--worker' in raw or '--help' in raw or '-h' in raw:return main()
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--max-runtime',type=float,default=180.)
    parser.add_argument('--output-dir',type=Path,default=None)
    known,_=parser.parse_known_args(raw)
    if not math.isfinite(known.max_runtime) or known.max_runtime<=0:
        print('max-runtime must be positive and finite',file=sys.stderr);return 2
    out=known.output_dir.resolve() if known.output_dir else HERE/'results'/'mars_ssi_20261109'/('run_'+datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6])
    if out.exists() and any(out.iterdir()):
        print('Output directory is not empty: choose a fresh directory',file=sys.stderr);return 2
    cmd=[sys.executable,str(Path(__file__).resolve()),*raw,'--worker']
    if known.output_dir is None:cmd+=['--output-dir',str(out)]
    status=run_isolated(cmd,known.max_runtime,join_managed_group=os.environ.get(MANAGED_ENV)=='1')
    out.mkdir(parents=True,exist_ok=True)
    worker_status=json.loads((out/'run_status.json').read_text()) if (out/'run_status.json').exists() else {}
    status=dict(worker_status,**status,safe_verdict=False,max_runtime_s=known.max_runtime)
    tmp=out/'run_status.json.tmp';tmp.write_text(json.dumps(status,indent=2));tmp.replace(out/'run_status.json')
    if status['status'] in ('timeout','interrupted'):
        print('[stopped]',status['status'],out,file=sys.stderr)
    return status['returncode']

if __name__=='__main__':raise SystemExit(supervised_cli())
