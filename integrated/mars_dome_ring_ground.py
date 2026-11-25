"""Unified tuned-steel dome + RC ring + nonlinear 3D ground research model.

Flat synthetic site, ENU, SI. One domain, one mass/selfweight assignment and one
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
from ground_hybrid import default_config,build_ground
from ground_hybrid.solver import initialize_gravity,apply_static_loads
from insight_motion import INSIGHT_ACC_PATH, preprocess_acc_mseed, select_motion_window, load_acceleration_csv


def axis_with_ports(width,ports,max_spacing):
    seeds=sorted(set([round(float(v),12) for v in ports]+[-width/2,0.,width/2]))
    if seeds[0]<-width/2 or seeds[-1]>width/2:raise ValueError('Ground must contain ring')
    values=[seeds[0]]
    for a,b in zip(seeds[:-1],seeds[1:]):
        values.extend(np.linspace(a,b,max(1,math.ceil((b-a)/max_spacing))+1)[1:])
    return np.array(values)


def assemble(args):
    """Caller owns wipe. Reuse RC ring ports as dome base nodes; no MP chains."""
    theta=np.arange(args.ring_segments)*2*np.pi/args.ring_segments
    axis=axis_with_ports(args.ground_width,np.r_[20*np.cos(theta),20*np.sin(theta)],args.spacing)
    X,Y=np.meshgrid(axis,axis);Z=np.zeros_like(X)
    c=default_config('Mars')
    c.update(width_m=args.ground_width,length_m=args.ground_width,depth_m=12.,core_width_m=40.,core_length_m=40.,
             core_spacing_m=args.spacing,max_outer_spacing_m=args.spacing,max_vertical_spacing_m=4.,lenses=[],
             tag_start=10000,max_elements=10000,elastic=args.elastic_ground)
    g=build_ground(c,(X,Y,Z))
    r=build_ring(RingConfig(segments=args.ring_segments),demo_supports=False)
    surface=[t for (i,j,k),t in g['index'].items() if k==0]
    contacts=[];pairs=[]
    if args.interface=='contact':
        from capmat.ssi import register_contact
        from capmat.mechanics import contact_parameters
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


def write_view(m,visual_tags,frames,path,scale):
    import plotly.graph_objects as go
    g,r=m['ground'],m['ring'];index={t:i for i,t in enumerate(visual_tags)}
    xyz=np.array([ops.nodeCoord(t) for t in visual_tags])
    # Top surface grid edges, not arbitrary synthetic oscillation.
    edges=set()
    for face in g['faces']:
        if face['side'] in ('top','bottom','xmin','xmax','ymin','ymax'):
            q=face['nodes']
            for a,b in zip(q,q[1:]+q[:1]):edges.add(tuple(sorted((a,b))))
    re=[(r['nodes'][i],r['nodes'][(i+1)%len(r['nodes'])]) for i in range(len(r['nodes']))]
    de=[(m['node_map'][a],m['node_map'][b]) for _,a,b in dome.DOME_EDGES]
    def lines(points,links,name,color):
        p=[]
        for a,b in links:p.extend([points[index[a]].tolist(),points[index[b]].tolist(),[None]*3])
        return go.Scatter3d(x=[v[0] for v in p],y=[v[1] for v in p],z=[v[2] for v in p],mode='lines',name=name,line=dict(color=color,width=2 if name=='Ground' else 5))
    def data(f):
        points=xyz+scale*np.array(f['du'])
        return [lines(points,edges,'Ground','#a78b68'),lines(points,re,'RC ring','#4a8bca'),lines(points,de,'Steel dome','#d85f36')]
    pf=[go.Frame(name=str(i),data=data(f),layout=go.Layout(title=f"Coupled ground–RC ring–dome | display ×{scale:g} | actual max increment {np.max(np.linalg.norm(f['du'],axis=1))*1000:.6g} mm")) for i,f in enumerate(frames)]
    fig=go.Figure(data=data(frames[0]),frames=pf)
    fig.update_layout(title=f"Coupled ground–RC ring–dome | displacement increment ×{scale:g} | research model",
        scene=dict(aspectmode='data',xaxis_title='East [m]',yaxis_title='North [m]',zaxis_title='Up [m]'),
        updatemenus=[dict(type='buttons',buttons=[dict(label='Play',method='animate',args=[None,dict(frame=dict(duration=100,redraw=True),transition=dict(duration=0),fromcurrent=True)]),dict(label='Pause',method='animate',args=[[None],dict(mode='immediate',frame=dict(duration=0,redraw=False))])])],
        sliders=[dict(steps=[dict(label=f"{f['t']:.2f}s",method='animate',args=[[str(i)],dict(mode='immediate',frame=dict(duration=0,redraw=True),transition=dict(duration=0))]) for i,f in enumerate(frames)])])
    fig.write_html(path,include_plotlyjs=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ground-width',type=float,default=60.)
    p.add_argument('--spacing',type=float,default=5.)
    p.add_argument('--ring-segments',type=int,default=12)
    p.add_argument('--interface',choices=['tied','contact'],default='tied')
    p.add_argument('--mass-ton',type=float,default=200.)
    p.add_argument('--cover-thickness',type=float,default=0.)
    p.add_argument('--net-pressure',type=float,default=0.,help='Pa; default is deliberately unpressurized connection check')
    p.add_argument('--duration',type=float,default=2.)
    p.add_argument('--motion-start',type=float,default=None)
    p.add_argument('--motion-scale',type=float,default=1.)
    p.add_argument('--visual-scale',type=float,default=1000000.,help='Display only; actual displacement remains in meters in CSV')
    p.add_argument('--static-only',action='store_true')
    p.add_argument('--elastic-ground',action='store_true')
    p.add_argument('--no-html',action='store_true')
    p.add_argument('--output-dir',type=Path,default=HERE/'results'/'mars_integrated_legacy')
    a=p.parse_args(argv)
    if not np.isfinite([a.ground_width,a.spacing,a.mass_ton,a.cover_thickness,a.net_pressure,a.duration,a.motion_scale,a.visual_scale]).all():p.error('Finite numeric inputs required')
    if a.ground_width<=42 or min(a.spacing,a.mass_ton,a.duration,a.visual_scale)<=0 or a.motion_scale<0 or a.cover_thickness<0:p.error('Ground width >42 m; positive spacing/mass/duration/scales and nonnegative cover required')
    if a.ring_segments<12 or a.ring_segments%12:p.error('ring-segments must be a multiple of 12')
    out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    status=dict(status='running',complete=False,safe_verdict=False)
    result=dict(parameters={**vars(a),'output_dir':str(out)},assumptions=[
        'Synthetic flat layered Mars analogue, not geolocated DEM at ELYSE',
        'Default tied translational interface can carry tension: no uplift or slip claim',
        'Optional contact is lumped at beam centerline nodes with tributary areas; not resolved footing pressure distribution',
        'Default 12-element polygonal RC ring is coarser than separate 48-element ring; hinge/mesh convergence not established',
        'C30/SD400 nominal catalogue RC; transferred Steel02 S355 coupon used for SHS',
        'All selfweights ramped together; no staged installation or excavation',
        'Fixed bottom, lateral rollers; one 3-component UniformExcitation input for whole assembled model; reflects waves',
        'Surface observation treated as base input for experiment; not deconvolved incident wave',
        'Default uncovered/unpressurized baseline; not prior pressurized 7.5 m cover design',
        'No automatic Rayleigh damping; hysteretic constitutive behavior only',
        'Apex minus mean base removes translation only, not full rigid-body rotation',
        'No strength/safety/optimal mass conclusion'])
    started=time.perf_counter();ops.wipe()
    def save():
        (out/'run_status.json').write_text(json.dumps(status,indent=2))
        (out/'summary.json').write_text(json.dumps(result,indent=2))
    save()
    try:
        m=assemble(a);g,r=m['ground'],m['ring']
        result['assembly']=dict(ground_nodes=len(g['nodes']),ground_elements=len(g['elements']),ring_nodes=len(r['nodes']),dome_members=len(m['dome_elements']),
            ground_mass_kg=g['mass_kg'],ring_mass_kg=r['total_mass_kg'],dome_membrane_cover_mass_kg=m['structure_mass_kg'],
            ground_materials=g['materials'],ring_materials=['concrete_c30_prototype','rebar_sd400_prototype'],steel=dome.TUNED_STEEL,
            mesh=g['mesh'],node_map=m['node_map'],ground_ring_pairs=m['ground_ring_pairs'])
        print('[assembly]',len(g['elements']),'soil bricks,',len(r['elements']),'RC beams,',len(m['dome_elements']),'steel beams',flush=True)
        result['gravity']=initialize_gravity(g,extra_loads=m['extra_gravity'])
        print('[gravity] equilibrated',result['gravity']['force_balance']['relative_error'],flush=True)
        if a.net_pressure:
            f=dome.pressure_node_forces(dome.dome_NODES,dome.DOME_EDGES,a.net_pressure)
            loads={m['node_map'][n]:[float(v[1]),float(v[0]),float(v[2]),0.,0.,0.] for n,v in f.items()}
            result['pressure']=apply_static_loads(g,loads)
        # Ground/ring interface constraints and actual dome base reuse.
        result['interface_error_m']=max(float(np.linalg.norm(np.array(ops.nodeDisp(x)[:3])-ops.nodeDisp(y)[:3])) for x,y in m['ground_ring_pairs'])
        if a.interface=='tied' and result['interface_error_m']>1e-8:raise RuntimeError('Ground/ring tied kinematics failed')
        if m['contacts']:
            from ground_hybrid.interface import contact_summary
            result['static_contact']=contact_summary({'contacts':m['contacts']})
        surface=sorted({n for f in g['faces'] if f['side']=='top' for n in f['nodes']})
        boundary_nodes={n for f in g['faces'] for n in f['nodes']}
        tags=sorted(set(list(boundary_nodes)+r['nodes']+list(m['node_map'].values())))
        U0=np.array([ops.nodeDisp(n)[:3] for n in tags]);ix={n:i for i,n in enumerate(tags)}
        frames=[dict(t=0.,du=np.zeros_like(U0).tolist())]
        result['static_displacement_m']=dict(ground=max(float(np.linalg.norm(ops.nodeDisp(n)[:3])) for n in surface),
            ring=max(float(np.linalg.norm(ops.nodeDisp(n)[:3])) for n in r['nodes']),dome=max(float(np.linalg.norm(ops.nodeDisp(n)[:3])) for n in m['node_map'].values()))
        if not a.static_only:
            full=preprocess_acc_mseed(INSIGHT_ACC_PATH,out/'motion_full_enu.csv')
            window=select_motion_window(full,out/'motion_window_enu.csv',duration=a.duration,start=a.motion_start)
            dt,acc=load_acceleration_csv(window)
            ops.loadConst('-time',0.)
            for i in range(3):
                tag=900000+i;ops.timeSeries('Path',tag,'-dt',dt,'-values',*acc[:,i].tolist(),'-factor',a.motion_scale)
                ops.pattern('UniformExcitation',tag,i+1,'-accel',tag)
            ops.wipeAnalysis();ops.constraints('Plain');ops.numberer('RCM');ops.system('UmfPack')
            ops.test('NormUnbalance',max(.001,g['mass_kg']*g['config']['gravity_m_s2']*1e-10),60)
            ops.algorithm('Newton');ops.integrator('Newmark',.5,.25);ops.analysis('Transient')
            apex=m['node_map'][6];base=r['dome_ports'];rows=[]
            dynamic_started=time.perf_counter()
            for step in range(len(acc)-1):
                if ops.analyze(1,dt):raise RuntimeError(f'Nonconvergence at dynamic step {step+1}; no failure/safety inference')
                U=np.array([ops.nodeDisp(n)[:3] for n in tags]);du=U-U0
                ground_max=max(float(np.linalg.norm(du[ix[n]])) for n in surface)
                ring_max=max(float(np.linalg.norm(du[ix[n]])) for n in r['nodes'])
                dome_max=max(float(np.linalg.norm(du[ix[n]])) for n in m['node_map'].values())
                relative=du[ix[apex]]-np.mean([du[ix[n]] for n in base],axis=0)
                rows.append([ops.getTime(),ground_max,ring_max,dome_max,float(np.linalg.norm(relative))])
                frames.append(dict(t=float(ops.getTime()),du=du.tolist()))
                print(f'[dynamic] {step+1}/{len(acc)-1}, ground={ground_max:.3e} m, dome={dome_max:.3e} m',flush=True)
            np.savetxt(out/'coupled_response.csv',rows,delimiter=',',header='time_s,ground_max_increment_m,ring_max_increment_m,dome_max_increment_m,apex_minus_mean_base_increment_m',comments='')
            if m['contacts']:result['final_contact']=contact_summary({'contacts':m['contacts']})
            result['dynamic']=dict(completed_steps=len(rows),dt_s=dt,
                input_peak_acceleration_vector_m_s2=float(a.motion_scale*np.max(np.linalg.norm(acc,axis=1))),peak_response_m=np.max(np.array(rows)[:,1:],axis=0).tolist(),
                peak_response_order=['ground','ring','dome','apex_relative_to_mean_base'],elapsed_s=time.perf_counter()-dynamic_started,
                motion=json.loads(window.with_suffix('.json').read_text()),scale=a.motion_scale)
        if not a.no_html:write_view(m,tags,frames,out/'mars_integrated_response.html',a.visual_scale)
        result['elapsed_s']=time.perf_counter()-started;status.update(status='complete',complete=True,elapsed_s=result['elapsed_s'])
        print('[complete]',result['elapsed_s'],'s',flush=True)
    except Exception as exc:
        status.update(status='failed',error=str(exc),elapsed_s=time.perf_counter()-started);raise
    finally:save();ops.wipe()

if __name__=='__main__':main()
