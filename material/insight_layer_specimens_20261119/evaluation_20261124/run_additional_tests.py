"""Literature-pressure monotonic shear tests of the four unchanged MGS records.

Homogeneous plane-strain simple shear is an envelope comparison, not a direct
reproduction of Dotson's 10.17cm shear-box displacement curves.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import specimens as base
ops=base.ops


def test(record,pressure):
    base.cube(record)
    for n in range(1,9):
        ops.fix(n,*((1,1,1) if n<=4 else (0,1,0)))
    for n in (6,7,8):
        ops.equalDOF(5,n,1,3)
    pdmy=record['family']=='PressureDependMultiYield02'
    ops.element('SSPbrick' if pdmy else 'bbarBrick',1,*range(1,9),1,0.,0.,0.)
    if pdmy:ops.updateMaterialStage('-material',1,'-stage',0)
    ops.timeSeries('Linear',1)
    ops.pattern('Plain',1,1)
    ops.load(5,0.,0.,-pressure)
    base.configure(.1)
    ops.test('NormDispIncr',1e-11,100,0) if pdmy else ops.test('NormUnbalance',1e-6,100,0)
    base.checked(10)
    ops.loadConst('-time',0.)
    if pdmy:
        ops.updateMaterialStage('-material',1,'-stage',1)
        ops.integrator('LoadControl',0.)
        base.checked()
    ops.timeSeries('Linear',2)
    ops.pattern('Plain',2,2)
    ops.load(5,1.,0.,0.)
    refinements=0
    def advance(target,depth=0):
        nonlocal refinements
        previous=ops.nodeDisp(5,1)
        ops.integrator('DisplacementControl',5,1,float(target-previous))
        if ops.analyze(1)==0:return
        if depth>=8:raise RuntimeError('Nonconvergence after 8 bisections')
        refinements+=1
        advance((target+previous)/2,depth+1);advance(target,depth+1)
    rows=[];failure=None
    for gamma in np.r_[0.,np.linspace(1e-6,1e-5,10),np.linspace(2e-5,.002,100)]:
        try:
            if gamma:advance(gamma)
        except RuntimeError as e:
            failure=str(e);break
        s=base.gp_stress().mean(axis=0)
        ops.reactions()
        rows.append([gamma,s[5],s[2],ops.getLoadFactor(2),sum(ops.nodeReaction(n,1) for n in range(1,5)),-float(np.mean(s[:3]))])
    out=HERE/'additional_tests'/record['id'];out.mkdir(parents=True,exist_ok=True)
    with (out/f'monotonic_p{int(pressure)}.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['gamma','tau_Pa','sigma_z_Pa','shear_force_N','base_reaction_N','mean_pressure_Pa']);w.writerows(rows)
    a=np.array(rows);scale=max(np.max(abs(a[:,1])),1.)
    eq=float(np.max(abs(a[:,3]+a[:,4]))/scale)
    vert=float(np.max(abs(a[:,2]+pressure))/pressure)
    tail=a[-20:,1]
    metrics=dict(id=record['id'],pressure_Pa=pressure,initial_G_Pa=float(np.polyfit(a[:6,0],a[:6,1],1)[0]),
                 peak_tau_Pa=float(np.max(a[:,1])),final_tau_Pa=float(a[-1,1]),final_gamma=float(a[-1,0]),
                 tail_range_relative=float(np.ptp(tail)/max(abs(tail.mean()),1.)),
                 force_balance_relative_error=eq,vertical_stress_relative_error=vert,
                 bisections=refinements,completed=failure is None,error=failure,
                 numerical_pass=bool(failure is None and eq<1e-5 and vert<1e-5))
    ops.wipe()
    return metrics


def main():
    records=base.load_materials()
    names=['sand_mgs1_dp_1300','sand_mgs1s_dp_1300','sand_mgs1_pdmy02_1300','sand_mgs1s_pdmy02_1300']
    result=dict(date='2026-11-24',scope='MC strength-envelope comparison only; not experimental stress-strain curve validation',cases=[])
    for name in names:
        record=next(r for r in records if r['id']==name)
        for pressure in [97.,250.,500.,680.]:
            try:r=test(record,pressure)
            except Exception as e:r=dict(id=name,pressure_Pa=pressure,numerical_pass=False,error=str(e))
            result['cases'].append(r)
            print(json.dumps(r),flush=True)
    result['all_numerical_pass']=all(r['numerical_pass'] for r in result['cases'])
    (HERE/'additional_tests/summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['all_numerical_pass']:raise SystemExit(1)


if __name__=='__main__':main()
