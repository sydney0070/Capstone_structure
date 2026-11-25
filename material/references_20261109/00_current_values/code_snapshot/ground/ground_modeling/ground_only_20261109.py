"""Mars ground ONLY: nonlinear gravity + ELYSE motion, no foundation/structure.

Public entry point: ground_20261107.py. This runner owns its OpenSees domain.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, sys, time, uuid
import numpy as np
import openseespy.opensees as ops
from ground_hybrid import default_config
from ground_hybrid.config import graded_axis
from ground_hybrid.run_journal import RunJournal, atomic_json
from ground_hybrid.io import serial
from ground_20261107 import (build_ground_20261107, static_ramp, ground_state,
    insight_reference_site, prepare_motion, resolution, lateral_dashpots, advance, require_finite)
from response_viewer_20261109 import save_response_viewer

HERE = Path(__file__).resolve().parent


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    for name,default in [('ground-width',60.),('spacing',5.),('depth',12.),
                         ('duration',2.),('warmup',5.),('dt',.05),('freqmax',1.),
                         ('max-runtime',180.)]:
        p.add_argument('--'+name,type=float,default=default)
    p.add_argument('--static-steps',type=int,default=20)
    p.add_argument('--max-frames',type=int,default=101)
    p.add_argument('--quick',action='store_true',help='1 s warmup + 1 s target; connection check only')
    p.add_argument('--static-only',action='store_true')
    p.add_argument('--ground-history',action='store_true')
    p.add_argument('--no-html',action='store_true')
    p.add_argument('--visual-scale',type=float,default=None,help='display only; default automatic')
    p.add_argument('--boundary',choices=['viscous-sides','fixed-box'],default='viscous-sides')
    p.add_argument('--output-dir',type=Path)
    p.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    return p


def main(argv=None):
    raw = sys.argv[1:] if argv is None else argv
    p = parser(); a = p.parse_args(raw)
    if a.quick:
        if '--warmup' not in raw:a.warmup=1.
        if '--duration' not in raw:a.duration=1.
    if (not np.isfinite([a.ground_width,a.spacing,a.depth,a.duration,a.warmup,a.dt,a.freqmax,a.max_runtime]).all()
            or min(a.ground_width,a.spacing,a.depth,a.duration,a.dt,a.max_runtime)<=0
            or a.warmup<0 or not .05<a.freqmax<5 or a.dt*a.freqmax>.05+1e-12
            or a.static_steps<1 or a.max_frames<2
            or (a.visual_scale is not None and (not math.isfinite(a.visual_scale) or a.visual_scale<=0))):
        p.error('invalid finite geometry, time or display parameters')
    out = a.output_dir.resolve() if a.output_dir else new_output()
    if out.exists() and any(out.iterdir()):p.error('choose a fresh output directory')
    out.mkdir(parents=True,exist_ok=True); started=time.perf_counter()
    result=dict(model='Mars ground only',complete=False,safe_verdict=False,
                parameters={**vars(a),'output_dir':str(out)},structure_count=0,foundation_count=0,contact_count=0)
    result['parameters']['output_dir']=str(out)
    def save():
        result['elapsed_s']=time.perf_counter()-started
        atomic_json(out/'summary.json',serial(result))
    ops.wipe(); ops.logFile(str(out/'opensees.log'))
    try:
        sources=[Path(__file__),HERE/'ground_20261107.py',HERE/'response_viewer_20261109.py',HERE/'insight_dynamic_viewer.py',HERE/'insight_motion.py']
        sources += [HERE/'ground_history_20261109.py',HERE/'insight_site_dem.py',HERE/'ground_modeling.py',HERE.parents[1]/'material/materials_20261107.py']
        sources += list((HERE.parents[1]/'SSI').glob('*.py'))
        sources += list((HERE/'ground_hybrid').rglob('*.py'))+list((HERE/'ground_hybrid/configs').glob('*.json'))
        result['code_sha256']={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
        result['runtime']=dict(python=sys.version,opensees=ops.version(),numpy=np.__version__)
        result['site'],samples=insight_reference_site(a.ground_width)
        np.savetxt(out/'site_dem_samples.csv',samples,delimiter=',',header='east_m,north_m,dem_elevation_m,dem_relative_m,modeled_z_m',comments='')
        c=default_config('Mars'); c.update(width_m=a.ground_width,length_m=a.ground_width,
            depth_m=a.depth,core_width_m=.7*a.ground_width,core_length_m=.7*a.ground_width,
            core_spacing_m=a.spacing,max_outer_spacing_m=a.spacing,max_vertical_spacing_m=4.,
            lenses=[],max_elements=20000,elastic=False,
            site_scope='ELYSE-centred flat ground only; analogue layers, not measured site stratigraphy')
        x=graded_axis(c['width_m'],c['core_width_m'],a.spacing,a.spacing)
        X,Y=np.meshgrid(x,x); g=build_ground_20261107(c,(X,Y,np.zeros_like(X)))
        tags=sorted(g['nodes']); coordinates=np.asarray([g['nodes'][t] for t in tags])
        # This gate is before adding any ground radiation-boundary dashpots.
        if set(ops.getEleTags())!=set(g['elements']) or any(len(ops.nodeDisp(t))!=3 for t in tags) or g['attachments']:
            raise RuntimeError('ground-only domain contains non-ground model components')
        result.update(ground_elements=len(g['elements']),ground_nodes=len(tags),mesh=g['mesh'],
                      ground_materials=g['materials'],ground_revision=g['revision_validation'],
                      assembly_gate='only 3DOF ground nodes and bbarBrick ground elements before boundary activation')
        result['resolution']=resolution(g,a.freqmax,a.dt); save()
        journal=RunJournal(out,dict(model=result['model'],parameters=result['parameters'],code_sha256=result['code_sha256']))
        history=None
        if a.ground_history:
            from ground_history_20261109 import GroundHistory
            history=GroundHistory(g,out/'static_ground_history.csv');history.capture()
        result['gravity']=static_ramp(g,g['gravity_loads'],steps=a.static_steps,
            on_commit=history.capture if history else None,attempt_observer=journal)
        if history:result['static_history']=history.save(out/'static_ground_history.csv')
        reference=require_finite([ops.nodeDisp(t)[:3] for t in tags],'ground gravity')
        frames=[reference.copy()];times=[0.];result['boundary']={'kind':'static_box'}
        if not a.static_only:
            if a.boundary=='viscous-sides':result['boundary']=lateral_dashpots(g)
            ops.loadConst('-time',0.)
            input_dt,acc,motion=prepare_motion(HERE/'MARS_DATA/S1222a_ACC.mseed',out,
                duration=a.duration,warmup=a.warmup,freqmax=a.freqmax)
            if any(channel.split('.')[1]!='ELYSE' for channel in motion['channels']):
                raise ValueError('ground site and motion station mismatch')
            result['motion']=motion
            for d in range(3):
                tag=900000+d;ops.timeSeries('Path',tag,'-dt',input_dt,'-values',*acc[:,d].tolist())
                ops.pattern('UniformExcitation',tag,d+1,'-accel',tag)
            count=int(round((len(acc)-1)*input_dt/a.dt))
            if not math.isclose(count*a.dt,(len(acc)-1)*input_dt,abs_tol=1e-8):
                raise ValueError('dt must divide selected motion duration')
            ops.wipeAnalysis();ops.constraints('Plain');ops.numberer('RCM');ops.system('UmfPack')
            ops.test('NormUnbalance',max(.001,g['mass_kg']*c['gravity_m_s2']*1e-10),80)
            ops.algorithm('Newton');ops.integrator('Newmark',.5,.25);ops.analysis('Transient')
            dynamic_history=None
            if a.ground_history:
                dynamic_history=GroundHistory(g,out/'dynamic_ground_history.csv');dynamic_history.capture()
            attempts=[];peak=[0.];commits=[0]
            selected=set(np.linspace(0,count,min(count+1,a.max_frames),dtype=int))
            with (out/'committed_response.csv').open('w',newline='') as stream:
                writer=csv.writer(stream);writer.writerow(['time_s','target_time_s','ground_max_increment_m'])
                def observe():
                    u=require_finite([ops.nodeDisp(t)[:3] for t in tags],'ground dynamic')
                    value=float(np.linalg.norm(u-reference,axis=1).max());peak[0]=max(peak[0],value);commits[0]+=1
                    writer.writerow([ops.getTime(),ops.getTime()-motion['target_start_model_s'],value]);stream.flush()
                    if dynamic_history:dynamic_history.capture()
                for step in range(count):
                    advance(a.dt,attempts=attempts,on_commit=observe)
                    if not a.no_html and step+1 in selected:
                        frames.append(require_finite([ops.nodeDisp(t)[:3] for t in tags],'viewer frame').copy());times.append(ops.getTime())
            if dynamic_history:result['dynamic_history']=dynamic_history.save(out/'dynamic_ground_history.csv')
            result['dynamic']=dict(steps=count,committed_samples=commits[0],max_increment_m=peak[0],attempts=attempts,
                scope='every committed substep; ground response only; relative to imposed base excitation')
        result['final_ground_state']=ground_state(g,gauss_points=8)
        if not a.no_html:
            site=result['site'];context=dict(title='화성 지반 · 동적 응답' if not a.static_only else '화성 지반 · 자중 응답',
                subtitle='InSight / ELYSE · 지반 단독 비선형 해석',
                summary=f"지반 {len(g['elements']):,} 요소 · 구조물·기초·연결부 없음",static_only=a.static_only,
                scope=f"ELYSE {site['latitude_deg']}°N, {site['longitude_deg']}°E 기준 평탄 지반. 문헌·모사토 기반 가정 층서. 실제 DEM 요철 미적용. 비선형 DP 토층·탄성 기반암. "+
                    ('자중 정적 응답.' if a.static_only else 'SEIS 지표 관측 가속도를 기저 운동으로 가정. 고정 바닥·'+('근사 점성 측면.' if a.boundary=='viscous-sides' else '측면 롤러.')))
            result['viewer']=save_response_viewer(g,tags,coordinates,frames,times,out/'ground_response.html',context=context,scale=a.visual_scale)
        result['complete']=True; save();print('[complete ground only]',out,flush=True)
    except BaseException as exc:
        result['error']=str(exc);save();raise
    finally:ops.wipe()


def new_output():
    from datetime import datetime
    return HERE.parent/'results'/'ground_20261109'/('run_'+datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6])


def supervised_cli():
    from ground_hybrid.process_control import run_isolated,MANAGED_ENV
    raw=sys.argv[1:]
    if '--worker' in raw or '--help' in raw or '-h' in raw:return main()
    a=parser().parse_args(raw)
    if not math.isfinite(a.max_runtime) or a.max_runtime<=0:parser().error('max-runtime must be positive')
    out=a.output_dir.resolve() if a.output_dir else new_output()
    if out.exists() and any(out.iterdir()):parser().error('choose a fresh output directory')
    command=[sys.executable,str(Path(__file__)),*raw,'--worker']
    if a.output_dir is None:command+=['--output-dir',str(out)]
    status=run_isolated(command,a.max_runtime,join_managed_group=os.environ.get(MANAGED_ENV)=='1')
    out.mkdir(parents=True,exist_ok=True);atomic_json(out/'run_status.json',dict(status,safe_verdict=False))
    return status['returncode']


if __name__=='__main__':raise SystemExit(supervised_cli())
