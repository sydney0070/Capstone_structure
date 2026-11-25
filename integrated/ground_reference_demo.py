"""Runnable entry point for mixed 3D ground and an illustrative raft/frame."""
import argparse
import json
import math
from datetime import datetime,timezone
from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PROJECT_ROOT), str(PROJECT_ROOT/'structure'),
               str(PROJECT_ROOT/'ground/ground_modeling'),
               str(PROJECT_ROOT/'material/material_release_20261101')]

import numpy as np
import openseespy.opensees as ops
from ground_hybrid import default_config,build_ground,add_raft,initialize_gravity,apply_static_loads
from ground_hybrid.interface import pressure_loads,raft_response
from ground_hybrid.io import export_results,load_dem,write_json


def add_demo_frame(raft,g,*,height_m=3.,total_mass_kg=4000.,tag_start=450000):
    """Four steel cantilevers reusing corner raft ports; not a dome surrogate.
    Independent columns verify translation, rotation and mass load transfer.
    """
    from capmat.materials import material
    bounds=raft['bounds_m'];corners=[(x,y) for x in bounds[:2] for y in bounds[2:]]
    ports=[]
    for x,y in corners:
        matches=[t for t in raft['nodes'] if abs(ops.nodeCoord(t)[0]-x)<1e-9 and abs(ops.nodeCoord(t)[1]-y)<1e-9]
        if len(matches)!=1:raise ValueError('raft corner port missing')
        ports.append(matches[0])
    ops.model('basic','-ndm',3,'-ndf',6);ops.geomTransf('Linear',tag_start,1.,0.,0.)
    m=material('steel_s355_nominal');tops=[];loads={}
    A=.01;Iy=Iz=8e-5;J=1.6e-4;member_mass=m['density_kg_m3']*A*height_m
    for i,p in enumerate(ports):
        top=tag_start+i;x,y,z=ops.nodeCoord(p);ops.node(top,x,y,z+height_m);tops.append(top)
        ops.element('elasticBeamColumn',tag_start+i,p,top,A,m['E_Pa'],m['G_Pa'],J,Iy,Iz,tag_start,'-mass',m['density_kg_m3']*A)
        # Element has mass; gravity still needs explicit equivalent joint loads.
        top_mass=total_mass_kg/4;ops.mass(top,top_mass,top_mass,top_mass,0.,0.,0.)
        loads[top]=[0.,0.,-(top_mass+member_mass/2)*g['config']['gravity_m_s2'],0.,0.,0.]
        loads[p]=[0.,0.,-member_mass/2*g['config']['gravity_m_s2'],0.,0.,0.]
    return dict(ports=ports,top_nodes=tops,gravity_loads=loads,member_mass_kg=4*member_mass,
                lumped_mass_kg=total_mass_kg,scope='four independent elastic steel columns, illustrative connection test only')


def main(argv=None):
    p=argparse.ArgumentParser(description='3D mixed Earth/Moon/Mars ground + SSI reference (SI).')
    p.add_argument('--body',choices=['Earth','Moon','Mars'],default='Earth')
    p.add_argument('--earth-model',choices=['dp','pdmy02'],default='dp',help='dp: static reference; pdmy02: advanced cyclic skeleton option')
    p.add_argument('--config',type=Path,help='full config JSON; overrides body/default geometry')
    p.add_argument('--terrain',choices=['flat','slope','dem'],default='flat')
    p.add_argument('--prepare-pad',action='store_true',help='explicit initial graded pad for DEM (not staged excavation)')
    p.add_argument('--no-raft',action='store_true')
    p.add_argument('--demo-frame',action='store_true')
    p.add_argument('--elastic',action='store_true',help='explicit elastic diagnostic instead of nonlinear soil')
    p.add_argument('--stiffness-factor',type=float,default=None,help='override JSON stiffness multiplier only when supplied')
    p.add_argument('--pressure-kpa',type=float,default=10.,help='illustrative extra uniform foundation pressure')
    p.add_argument('--gravity-only',action='store_true')
    p.add_argument('--no-html',action='store_true')
    p.add_argument('--no-open',action='store_true',help='save HTML without automatically opening the browser')
    p.add_argument('--html-frames',type=int,default=151,help='maximum saved animation frames, including initial state; does not change integration steps')
    p.add_argument('--output',type=Path)
    p.add_argument('--dynamic',choices=['none','rest','pulse','incident'],default='none')
    p.add_argument('--dt',type=float,default=.001)
    p.add_argument('--steps',type=int,default=100)
    p.add_argument('--verify',action='store_true')
    args=p.parse_args(argv)
    if args.html_frames<2:p.error('--html-frames must be at least2')
    if args.verify:
        from ground_hybrid.verification import main as verify
        return verify()
    c=json.loads(args.config.read_text()) if args.config else default_config(args.body,earth_model=args.earth_model)
    if args.elastic:c['elastic']=True
    if args.stiffness_factor is not None:c['stiffness_factor']=args.stiffness_factor
    if not math.isfinite(args.pressure_kpa):p.error('--pressure-kpa must be finite')
    if args.dynamic!='none' and (not math.isfinite(args.dt) or args.dt<=0 or args.steps<1):
        p.error('dynamic analyses require positive --dt and --steps')
    output=args.output or PROJECT_ROOT/'integrated/results/reference_demo'/c['body']
    output.mkdir(parents=True,exist_ok=True)
    started=datetime.now(timezone.utc).isoformat()
    write_json(output/'RUN_STATUS.json',dict(status='running',started_utc=started,complete=False))
    terrain=None;terrain_meta=dict(kind='synthetic_flat',native_resolution_m=None)
    ops.wipe() # only the CLI owns/clears its complete domain
    try:
        if args.terrain=='dem':terrain,terrain_meta=load_dem(c,prepare_pad=args.prepare_pad)
        elif args.terrain=='slope':
            from ground_hybrid.config import graded_axis
            x=graded_axis(c['width_m'],c['core_width_m'],c['core_spacing_m'],c['max_outer_spacing_m'])
            y=graded_axis(c['length_m'],c['core_length_m'],c['core_spacing_m'],c['max_outer_spacing_m'])
            X,Y=np.meshgrid(x,y);Z=.02*X+.01*Y;terrain=(X,Y,Z);terrain_meta=dict(kind='synthetic_slope',slope_x=.02,slope_y=.01)
        g=build_ground(c,terrain);raft=None;frame=None
        if not args.no_raft:
            bounds=(-c['core_width_m']/2,c['core_width_m']/2,-c['core_length_m']/2,c['core_length_m']/2)
            raft=add_raft(g,bounds)
        if args.demo_frame:
            if raft is None:raise ValueError('--demo-frame needs raft')
            frame=add_demo_frame(raft,g)
        print(f"Built {c['body']}: {len(g['nodes'])} nodes, {len(g['elements'])} bricks; materials={list(g['material_volume_fractions'])}",flush=True)
        initialize_gravity(g,extra_loads=frame['gravity_loads'] if frame else None)
        from ground_hybrid.solver import initial_stress_summary
        write_json(output/'gravity_stress.json',initial_stress_summary(g))
        extra=dict(terrain=terrain_meta,frame=frame);visualization=None
        if not args.gravity_only and raft and args.pressure_kpa:
            extra['static']=apply_static_loads(g,pressure_loads(raft,args.pressure_kpa*1000))
        if args.dynamic!='none':
            from ground_hybrid.dynamics import activate_absorbing_boundaries,run_transient,add_incident_base_velocity,mesh_frequency_guidance
            extra['dynamic_boundary']=activate_absorbing_boundaries(g)
            nodes=(frame['top_nodes'] if frame else raft['nodes'] if raft else [max(g['nodes'])])
            if args.dynamic=='incident':
                t=np.arange(args.steps+1)*args.dt
                v=.001*np.sin(2*np.pi*5*t)*np.sin(np.pi*t/t[-1])**2
                extra['input']=add_incident_base_velocity(g,v,args.dt)
            elif args.dynamic=='pulse':
                t=np.arange(args.steps+1)*args.dt
                f=100.*np.sin(np.pi*np.minimum(t/.02,1.))**2;f[t>.02]=0.
                ops.timeSeries('Path',800030,'-dt',args.dt,'-values',*f.tolist());ops.pattern('Plain',800030,800030)
                for n in nodes:ops.load(n,*([1/len(nodes),0.,0.]+([0.,0.,0.] if len(ops.nodeDisp(n))==6 else [])))
                extra['input']=dict(kind='100N total synthetic horizontal pulse',duration_s=.02)
            extra['frequency_guidance']=mesh_frequency_guidance(g,5.)
            from ground_hybrid.viewer import visualization_node_tags
            visual_nodes=None if args.no_html else visualization_node_tags(g,frame)
            transient=run_transient(g,dt_s=args.dt,steps=args.steps,record_nodes=nodes,
                                    visualization_nodes=visual_nodes,max_visualization_frames=args.html_frames)
            visualization=transient['visualization']
            write_json(output/'transient.json',transient);extra['dynamic']=dict(steps=args.steps,dt_s=args.dt,input=args.dynamic)
        summary=export_results(g,output,raft=raft,extra=extra,html=not args.no_html,visualization=visualization)
        (output/'FAILED.json').unlink(missing_ok=True)
        write_json(output/'RUN_STATUS.json',dict(status='success',complete=True,started_utc=started,
                                                finished_utc=datetime.now(timezone.utc).isoformat()))
        if raft:
            r=summary['raft_response'];print('Additional mean settlement [mm]:',1000*r['mean_additional_settlement_m'],flush=True)
        print('Saved:',output,flush=True)
        if not args.no_html:
            path=(output/'ground_viewer.html').resolve();print('HTML viewer:',path,flush=True)
            if not args.no_open:
                import webbrowser
                try:
                    if not webbrowser.open(path.as_uri()):print('Browser could not open automatically; open the HTML file above.',flush=True)
                except Exception as e:print('Browser open failed (analysis succeeded):',e,flush=True)
        return 0
    except Exception as e:
        import traceback
        write_json(output/'FAILED.json',dict(error=str(e),traceback=traceback.format_exc(),config=c))
        write_json(output/'RUN_STATUS.json',dict(status='failed',complete=False,started_utc=started,error=str(e)))
        raise

if __name__=='__main__':raise SystemExit(main())
