"""Diagnostic reproduction of RC cyclic beam: force balance, reversal, mesh sensitivity.
OpenSeesPy required. SI N-m-Pa. Results are NOT experimental validation.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from dataclasses import asdict, replace
from pathlib import Path
import traceback
import numpy as np
import openseespy.opensees as ops
from rc_section import SectionGeometry, create_rc_fiber_section

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results_rc_diagnostics'
TARGETS = (0.003, -0.003, 0.006, -0.006, 0.012, -0.012, 0.018, -0.018, 0.0)
FIELDS = ['case','phase','target_m','segment','substep','u_m','load_factor_N','reaction_x_N','force_imbalance_N','relative_force_imbalance','base_moment_Nm','element_local_0','element_local_1','element_local_2','algorithm','status']


def setup(ny, nz, nip):
    ops.wipe()
    ops.model('basic','-ndm',2,'-ndf',3)
    geom=replace(SectionGeometry(),ny=ny,nz=nz)
    create_rc_fiber_section(ops,geom=geom)
    ops.node(1,0.0,0.0); ops.node(2,0.0,3.0)
    ops.fix(1,1,1,1)
    ops.geomTransf('PDelta',1)
    ops.beamIntegration('Lobatto',1,10,nip)
    ops.element('forceBeamColumn',1,1,2,1,1)
    ops.timeSeries('Linear',1)
    ops.pattern('Plain',1,1)
    ops.load(2,1.0,0.0,0.0) # actual reference horizontal load = 1 N
    ops.constraints('Plain');ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr',1e-9,80)
    ops.algorithm('Newton')
    ops.integrator('LoadControl',1.0)
    ops.analysis('Static')
    return geom


def snapshot(case,phase,target,segment,substep,algorithm='Newton',status='ok'):
    ops.reactions()
    u=float(ops.nodeDisp(2,1))
    force=float(ops.getLoadFactor(1))
    rx=float(ops.nodeReaction(1,1))
    bm=float(ops.nodeReaction(1,3))
    elem=list(map(float,ops.eleResponse(1,'localForce')))
    mismatch=force+rx
    # relative-to-peak force computed after the run; this local scaling helps diagnose.
    return dict(case=case,phase=phase,target_m=target,segment=segment,substep=substep,u_m=u,
                load_factor_N=force,reaction_x_N=rx,force_imbalance_N=mismatch,
                relative_force_imbalance=abs(mismatch)/max(abs(force),abs(rx),1.0),
                base_moment_Nm=bm,
                element_local_0=elem[0] if len(elem)>0 else math.nan,
                element_local_1=elem[1] if len(elem)>1 else math.nan,
                element_local_2=elem[2] if len(elem)>2 else math.nan,
                algorithm=algorithm,status=status)


def attempt(du):
    # Same trial state is rolled back by OpenSees on unsuccessful analyze(1).
    strategies=[('Newton','NormDispIncr',1e-9,80),
                ('NewtonLineSearch','NormDispIncr',1e-8,100),
                ('KrylovNewton','NormDispIncr',1e-8,100)]
    for algo,test,tol,iters in strategies:
        ops.test(test,tol,iters,0)
        ops.algorithm(algo)
        ops.integrator('DisplacementControl',2,1,du)
        if ops.analyze(1)==0:
            return True,algo
    return False,'all_failed'


def run(case,ny,nz,nip,steps,adaptive=False):
    geom=setup(ny,nz,nip)
    rows=[]
    rows.append(snapshot(case,'start',0,0,0))
    previous=0.0
    failure=None
    for segment,target in enumerate(TARGETS,1):
        du=(target-previous)/steps
        for substep in range(1,steps+1):
            ok,algo=attempt(du)
            if not ok and adaptive:
                # Recovery only for the failing *substep*, using two half steps.
                # If either half step fails the full case is invalidated.
                ok1,a1=attempt(du/2)
                if ok1:
                    rows.append(snapshot(case,'cyclic_half_step',target,segment,substep,a1))
                    ok2,a2=attempt(du/2)
                    ok=ok2
                    algo=f'{a1}+{a2}' if ok2 else 'half_step_2_failed'
                else:
                    algo='half_step_1_failed'
            if not ok:
                failure=f'Failed at target={target}, segment={segment}, substep={substep}, du={du}'
                break
            rows.append(snapshot(case,'cyclic',target,segment,substep,algo))
        if failure:break
        previous=target
    return rows,{'case':case,'geometry':asdict(geom),'nip':nip,'steps':steps,'status':'FAILED' if failure else 'COMPLETED','error':failure}


def summarize(rows,info):
    good=[r for r in rows if r['phase'].startswith('cyclic')]
    if not good: return info
    fx=np.asarray([r['load_factor_N'] for r in good],dtype=float)
    residual=np.asarray([r['force_imbalance_N'] for r in good],dtype=float)
    rx=np.asarray([r['reaction_x_N'] for r in good],dtype=float)
    denom=max(float(np.max(np.abs(fx))),float(np.max(np.abs(rx))),1.0)
    idx=int(np.argmax(np.abs(residual)))
    info.update({'points':len(good),'peak_abs_force_kN':float(np.max(np.abs(fx))/1000),
                 'max_balance_abs_N':float(np.max(np.abs(residual))),
                 'max_balance_rel_peak':float(np.max(np.abs(residual))/denom),
                 'worst_balance_at':{k:good[idx][k] for k in ('segment','substep','u_m','load_factor_N','reaction_x_N')},
                 'endpoint_u_m':float(good[-1]['u_m']),
                 'endpoint_F_kN':float(good[-1]['load_factor_N']/1000),
                 'equilibrium_pass':bool(np.max(np.abs(residual))/denom<1e-4)})
    info['validation_pass']=bool(info['status']=='COMPLETED' and info['equilibrium_pass'])
    return info


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--adaptive',action='store_true',help='Try half steps only if a full step fails')
    parser.add_argument('--steps',type=int,default=60)
    args=parser.parse_args()
    if args.steps<2:raise ValueError('steps must be >=2')
    OUT.mkdir(exist_ok=True)
    cases=[('baseline_24x12_ip5',24,12,5),('mid_36x18_ip5',36,18,5),
           ('refined_48x24_ip5',48,24,5),('base_24x12_ip7',24,12,7),
           ('refined_48x24_ip7',48,24,7)]
    overall={'model_status':'PROVISIONAL_NOT_CALIBRATED','adaptive':args.adaptive,'steps':args.steps,'cases':{}}
    for name,ny,nz,nip in cases:
        try:
            rows,info=run(name,ny,nz,nip,args.steps,args.adaptive)
            info=summarize(rows,info)
            write_csv(OUT/f'{name}.csv',rows)
        except Exception as exc:
            info={'status':'EXCEPTION','error':str(exc),'traceback':traceback.format_exc()}
        overall['cases'][name]=info
        print(name,info.get('status'), 'balance_rel=',info.get('max_balance_rel_peak'),'peak_kN=',info.get('peak_abs_force_kN'))
        (OUT/'diagnostic_summary.json').write_text(json.dumps(overall,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Saved:',OUT)
    print('Do NOT classify COMPLETED as physically valid unless equilibrium_pass is true.')

if __name__=='__main__':main()
