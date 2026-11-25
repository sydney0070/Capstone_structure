"""Element-formulation diagnostics for the user's RC fiber cantilever.
Only for numerical verification. Original Concrete02/Steel02 parameters remain unchanged.
SI: N, m, Pa. Requires openseespy and numpy.
"""
import argparse, csv, json, math, traceback
from dataclasses import replace, asdict
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
from rc_section import SectionGeometry, create_rc_fiber_section

TARGETS = [0.003,-0.003,0.006,-0.006,0.012,-0.012,0.018,-0.018,0.0]
ROOT = Path(__file__).resolve().parent
FIELDS=['segment','substep','target_m','u_m','force_N','reaction_x_N','balance_N','base_moment_Nm','algorithm','status']


def section_invariants(geom):
    # Geometry is held invariant when only fiber subdivisions change.
    Acon = geom.width*geom.depth
    Abar = math.pi*geom.bar_diameter**2/4
    Asteel = (geom.bars_top+geom.bars_bottom)*Abar
    # Existing patch covers the full section; steel is SUPERIMPOSED.
    return {'nominal_concrete_patch_area_m2':Acon, 'steel_area_m2':Asteel,
            'steel_area_fraction_of_concrete_patch':Asteel/Acon,
            'note':'Concrete patch includes steel locations; overlap is present in all runs. No correction performed.'}


def model(ny,nz,nip,form,element_iter):
    ops.wipe()
    ops.model('basic','-ndm',2,'-ndf',3)
    geom=replace(SectionGeometry(),ny=ny,nz=nz)
    create_rc_fiber_section(ops,geom=geom)
    ops.node(1,0.0,0.0);ops.node(2,0.0,3.0)
    ops.fix(1,1,1,1)
    ops.geomTransf('PDelta',1)
    ops.beamIntegration('Lobatto',1,10,nip)
    if form=='force':
        ops.element('forceBeamColumn',1,1,2,1,1,'-iter',element_iter,1e-10)
    elif form=='disp':
        ops.element('dispBeamColumn',1,1,2,1,1)
    else:raise ValueError(form)
    ops.timeSeries('Linear',1)
    ops.pattern('Plain',1,1)
    ops.load(2,1.0,0.0,0.0)
    ops.constraints('Plain'); ops.numberer('RCM');ops.system('BandGeneral')
    ops.test('NormDispIncr',1e-9,80,0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl',0.0)
    ops.analysis('Static')
    return geom


def record(seg,sub,target,alg):
    ops.reactions()
    u=float(ops.nodeDisp(2,1)); f=float(ops.getLoadFactor(1)); r=float(ops.nodeReaction(1,1))
    return {'segment':seg,'substep':sub,'target_m':target,'u_m':u,'force_N':f,
            'reaction_x_N':r,'balance_N':f+r,'base_moment_Nm':float(ops.nodeReaction(1,3)),
            'algorithm':alg,'status':'converged'}


def step(du, retry):
    methods=[('Newton',1e-9,80)]
    if retry:methods += [('NewtonLineSearch',1e-9,120),('KrylovNewton',1e-9,120)]
    for alg,tol,iterations in methods:
        ops.test('NormDispIncr',tol,iterations,0)
        ops.algorithm(alg)
        ops.integrator('DisplacementControl',2,1,du)
        if ops.analyze(1)==0:return True,alg
    return False,'all_methods_failed'


def run(name,ny,nz,nip,form,niter,steps,retry):
    geom=model(ny,nz,nip,form,niter)
    history=[]; fail=None; previous=0.0
    for segment,target in enumerate(TARGETS,1):
        du=(target-previous)/steps
        for substep in range(1,steps+1):
            good,alg=step(du,retry)
            if not good:
                fail={'segment':segment,'substep':substep,'target_m':target,'du_m':du}
                break
            snapshot=record(segment,substep,target,alg)
            history.append(snapshot)
            # Abort even on success if reactions diverge from applied force.
            denom=max(abs(snapshot['force_N']),abs(snapshot['reaction_x_N']),1000.)
            if not np.isfinite(snapshot['force_N']) or abs(snapshot['balance_N'])/denom > 1e-3:
                fail={'segment':segment,'substep':substep,'reason':'equilibrium_gate',
                      'force_N':snapshot['force_N'],'reaction_x_N':snapshot['reaction_x_N']}
                break
            if max(abs(snapshot['force_N']),abs(snapshot['reaction_x_N'])) > 1e7:
                fail={'segment':segment,'substep':substep,'reason':'force_guard_10MN'}
                break
        if fail:break
        previous=target
    force=np.array([h['force_N'] for h in history])
    rx=np.array([h['reaction_x_N'] for h in history])
    mismatch=np.array([h['balance_N'] for h in history])
    balance_rel=float(np.max(np.abs(mismatch))/max(np.max(np.abs(force)),np.max(np.abs(rx)),1.0)) if len(history) else None
    info={'case':name,'formulation':form,'section':asdict(geom),'section_invariants':section_invariants(geom),
          'nip':nip,'element_iterations':niter if form=='force' else None,
          'status':'COMPLETED' if fail is None else 'FAILED','failure':fail,'points':len(history),
          'peak_force_kN':float(np.max(np.abs(force))/1000) if len(history) else None,
          'balance_rel_peak':balance_rel, 'equilibrium_pass': bool(fail is None and balance_rel is not None and balance_rel<1e-4)}
    return history,info


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--steps',type=int,default=60)
    p.add_argument('--retry',action='store_true',help='add NewtonLineSearch and KrylovNewton recovery')
    args=p.parse_args()
    if args.steps<5:raise ValueError('steps must >= 5')
    out=ROOT/'results_rc_element_comparison';out.mkdir(exist_ok=True)
    cases=[('force_24x12_ip5_it30',24,12,5,'force',30),
           ('force_36x18_ip5_it30',36,18,5,'force',30),
           ('force_48x24_ip5_it30',48,24,5,'force',30),
           ('force_48x24_ip5_it50',48,24,5,'force',50),
           ('force_24x12_ip7_it30',24,12,7,'force',30),
           ('force_48x24_ip7_it30',48,24,7,'force',30),
           ('disp_24x12_ip5',24,12,5,'disp',0),
           ('disp_36x18_ip5',36,18,5,'disp',0),
           ('disp_48x24_ip5',48,24,5,'disp',0)]
    summary={'status':'DIAGNOSTIC_ONLY_NOT_CALIBRATED','steps':args.steps,'retry':args.retry,'cases':{}}
    for name,ny,nz,nip,form,niter in cases:
        try:
            rows,info=run(name,ny,nz,nip,form,niter,args.steps,args.retry)
            with (out/f'{name}.csv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=FIELDS);writer.writeheader();writer.writerows(rows)
        except Exception as e:
            info={'case':name,'status':'EXCEPTION','error':str(e),'traceback':traceback.format_exc()}
        summary['cases'][name]=info
        (out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
        print(name,info['status'],'peak_kN=',info.get('peak_force_kN'), 'balance=',info.get('balance_rel_peak'), 'failure=',info.get('failure'))
    print('Saved:',out)

if __name__=='__main__':main()
