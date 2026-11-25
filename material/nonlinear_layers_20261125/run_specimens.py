"""Data -> 1m cube -> multi-path tests -> response files (OpenSees SI).

Independent worker per material prevents a non-returning native integration call
from blocking all other candidates. Failure is retained, not fitted away.
"""
import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import openseespy.opensees as ops

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'insight_layer_specimens_20261119'))
import specimens as base


def records():return json.loads((HERE/'materials.json').read_text())['materials']


def save(path,data):path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def curve(path,header,rows):
    with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(header);w.writerows(rows)


def checked():
    if ops.analyze(1):raise RuntimeError('OpenSees analysis failed')


def axial(record,p,out,axis=2,tension=False):
    base.cube(record)
    for n,xyz in enumerate(base.XYZ,1):ops.fix(n,*[int(x==0) for x in xyz])
    ops.element('bbarBrick',1,*range(1,9),1,0.,0.,0.)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
    for n,xyz in enumerate(base.XYZ,1):ops.load(n,*[-p*x/4 for x in xyz])
    base.configure(.1)
    if ops.analyze(10):raise RuntimeError('Hydrostatic initialization failed')
    ops.loadConst('-time',0.)
    face=[n for n,xyz in enumerate(base.XYZ,1) if xyz[axis]==1]
    u0=ops.nodeDisp(face[0],axis+1)
    par=record['parameters'];s=math.sin(math.radians(par['phi_deg']));c=par['c_Pa'];cs=math.cos(math.radians(par['phi_deg']))
    qy=6*c*cs/(3+s) if tension else (2*p*s+2*c*cs)/(1-s)
    ey=qy/par['E_Pa'];sign=-1. if tension else 1.
    strains=sign*np.r_[np.linspace(0,.8*ey,41),np.linspace(.8*ey,3*ey,111)[1:],np.linspace(3*ey,2*ey,61)[1:]]
    ops.timeSeries('Path',2,'-dt',1.,'-values',*(u0-strains).tolist(),'-useLast');ops.pattern('Plain',2,2)
    for n in face:ops.sp(n,axis+1,1.)
    ops.wipeAnalysis();base.configure(0.);checked();ops.integrator('LoadControl',1.)
    rows=[];lat=[i for i in range(3) if i!=axis];max_lat=0.;max_balance=0.;max_spread=0.
    for i in range(len(strains)):
        if i:checked()
        stress=base.gp_stress();m=stress.mean(axis=0);q=-m[axis]+m[lat].mean()
        max_lat=max(max_lat,float(np.max(abs(stress[:,lat]+p))));max_spread=max(max_spread,float(np.max(abs(stress-m))))
        ops.reactions();max_balance=max(max_balance,abs(sum(ops.nodeReaction(n,axis+1) for n in face)+q))
        state=np.asarray([ops.eleResponse(1,'material',gp,'state') for gp in range(1,9)])
        rows.append([i,u0-ops.nodeDisp(face[0],axis+1),q,*m,float(state[:,4].max())])
    a=np.array(rows);E=float(np.polyfit(a[:20,1],a[:20,2],1)[0]);Eu=float(np.polyfit(a[151:,1],a[151:,2],1)[0]);actual=float(a[150,2]);scale=max(p,qy,1.)
    metrics=dict(test='uniaxial_tension_unload' if tension else 'triaxial_compression_unload',axis='xyz'[axis],pressure_Pa=p,
                 expected_yield_q_Pa=sign*qy,computed_yield_q_Pa=actual,yield_error_pct=100*abs(actual/(sign*qy)-1),
                 initial_E_Pa=E,initial_E_error_pct=100*abs(E/par['E_Pa']-1),unloading_E_Pa=Eu,unloading_E_error_pct=100*abs(Eu/par['E_Pa']-1),
                 final_q_Pa=float(a[-1,2]),final_axial_strain=float(a[-1,1]),max_plastic_strain_norm=float(a[:,9].max()),
                 lateral_error_relative=max_lat/max(p,1.),force_balance_relative=max_balance/scale,gauss_point_spread_relative=max_spread/scale)
    checks=dict(E=abs(E/par['E_Pa']-1)<1e-5,unload_E=abs(Eu/par['E_Pa']-1)<1e-5,yield_=abs(actual/(sign*qy)-1)<1e-4,
                confinement=max_lat/max(p,1.)<1e-5,force_balance=max_balance/scale<1e-5,homogeneous=max_spread/scale<1e-5,
                finite=bool(np.isfinite(a).all()),plastic_exercised=float(a[:,9].max())>1e-9)
    metrics.update(checks={k:bool(v) for k,v in checks.items()},passed=all(checks.values()))
    if 'constants' in record['calibration'] and not tension:
        from build_materials import hb_sigma1
        cal=record['calibration'];hb=float(hb_sigma1(p,cal['constants'],cal['sigma_ci_Pa']))
        metrics.update(HB_sigma1_Pa=hb,computed_sigma1_Pa=actual+p,HB_sigma1_difference_pct=100*((actual+p)/hb-1),HB_fit_within_10pct=abs((actual+p)/hb-1)<=.10)
    if tension and 'constants' in record['calibration']:
        cal=record['calibration'];hb=cal['s'] if 's' in cal else cal['constants']['s']
        hb_strength=hb*cal['sigma_ci_Pa']/cal['constants']['mb']
        metrics.update(HB_tensile_strength_Pa=hb_strength,DP_tension_to_HB_difference_pct=100*(abs(actual)/hb_strength-1),
                       tension_scope='Compression-fitted DP cannot independently match HB tensile cutoff. Do not claim tensile validation.')
    filename=('tension' if tension else 'triaxial')+f'_axis{axis}_p{int(p)}.csv'
    curve(out/filename,['step','axial_compression_strain','q_Pa','sx_Pa','sy_Pa','sz_Pa','sxy_Pa','syz_Pa','szx_Pa','plastic_strain_norm'],rows)
    return metrics


def shear(record,p,out,cyclic=False):
    base.cube(record)
    for n in range(1,9):ops.fix(n,*((1,1,1) if n<=4 else (0,1,0)))
    for n in (6,7,8):ops.equalDOF(5,n,1,3)
    ops.element('bbarBrick',1,*range(1,9),1,0.,0.,0.)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(5,0.,0.,-p)
    base.configure(.1)
    if ops.analyze(10):raise RuntimeError('Shear normal loading failed')
    ops.loadConst('-time',0.);ops.timeSeries('Linear',2);ops.pattern('Plain',2,2);ops.load(5,1.,0.,0.)
    par=record['parameters'];phi=math.radians(par['phi_deg'])
    tau_theory=2*math.sqrt(3)*math.cos(phi)/(3-math.sin(phi))*(par['c_Pa']+p*math.tan(phi))
    gy=tau_theory/par['G_Pa']
    if cyclic:
        amplitude=max(5*gy,.002);targets=np.r_[np.linspace(0,.001*gy,21),np.linspace(.001*gy,amplitude,81)[1:],np.linspace(amplitude,-amplitude,161)[1:],np.linspace(-amplitude,amplitude,161)[1:],np.linspace(amplitude,0,81)[1:]]
    else:targets=np.r_[np.linspace(0,.001*gy,21),np.linspace(.02*gy,.8*gy,40),np.linspace(gy,max(.002,20*gy),151)]
    rows=[];bisections=0
    def advance(target,depth=0):
        nonlocal bisections
        previous=ops.nodeDisp(5,1);ops.integrator('DisplacementControl',5,1,float(target-previous))
        if ops.analyze(1)==0:return
        if depth>=8:raise RuntimeError('Shear failed after8 bisections')
        bisections+=1;advance((previous+target)/2,depth+1);advance(target,depth+1)
    for target in targets:
        if target:advance(target)
        elif rows:advance(target)
        s=base.gp_stress().mean(axis=0);ops.reactions();F=ops.getLoadFactor(2);R=sum(ops.nodeReaction(n,1) for n in range(1,5))
        state=np.asarray([ops.eleResponse(1,'material',gp,'state') for gp in range(1,9)])
        rows.append([target,s[5],s[2],F,R,-float(s[:3].mean()),ops.nodeDisp(5,3),float(state[:,4].max())])
    a=np.asarray(rows);scale=max(float(np.max(abs(a[:,1]))),1.);eq=float(np.max(abs(a[:,3]+a[:,4]))/scale);vertical=float(np.max(abs(a[:,2]+p)))/p
    metrics=dict(test='cyclic_simple_shear' if cyclic else 'monotonic_simple_shear',pressure_Pa=p,steps=len(rows),bisections=bisections,
                 initial_G_Pa=float(np.polyfit(a[:6,0],a[:6,1],1)[0]),peak_tau_Pa=float(np.max(abs(a[:,1]))),final_tau_Pa=float(a[-1,1]),final_gamma=float(a[-1,0]),
                 DP_steady_shear_Pa=tau_theory,steady_shear_formula_difference_pct=None if cyclic else 100*(float(a[-1,1])/tau_theory-1),
                 force_balance_relative=eq,vertical_stress_relative=vertical,final_plastic_strain_norm=float(a[-1,7]),
                 shear_work_Pa=float(np.trapz(a[:,1],a[:,0])),work_scope='Shear work along this path, not a measured damping ratio or complete dissipation measure')
    cal=record['calibration'];target=cal['target_MC_c_Pa']+p*math.tan(math.radians(cal['target_MC_phi_deg']))
    if not cyclic:metrics.update(target_MC_tau_Pa=target,target_MC_difference_pct=100*(float(a[-1,1])/target-1),target_scope='Calibration/analogue envelope, not independent raw experimental stress-strain points')
    checks=dict(finite=bool(np.isfinite(a).all()),force_balance=eq<1e-5,vertical=vertical<1e-5,
                initial_G=abs(metrics['initial_G_Pa']/par['G_Pa']-1)<1e-4)
    if not cyclic:checks['DP_formula']=abs(float(a[-1,1])/tau_theory-1)<1e-4
    metrics.update(checks={k:bool(v) for k,v in checks.items()},passed=all(checks.values()))
    curve(out/(('cyclic' if cyclic else 'shear')+f'_p{int(p)}.csv'),['gamma','tau_Pa','sigma_z_Pa','force_N','reaction_N','mean_pressure_Pa','vertical_displacement_m','plastic_strain_norm'],rows)
    return metrics


def worker(name):
    r=next(x for x in records() if x['id']==name);out=HERE/'results'/name;out.mkdir(parents=True,exist_ok=True)
    save(out/'material.json',r);result=dict(id=name,tests=[],passed=False,independent_Mars_validation=False)
    save(out/'summary.json',result)
    try:
        result['mass']=base.verify_mass(r);save(out/'summary.json',result)
        layer=r['layer']
        pressures=[0.,50000.,200000.,500000.] if layer=='fractured_basalt' else [0.,10e6,20e6] if layer=='laboratory_reference_only' else [5000.,50000.,200000.] if layer=='blocky_ejecta' else [500.,5000.,50000.]
        for p in pressures:
            result['tests'].append(axial(r,p,out));save(out/'summary.json',result)
        # One matched pressure in all three axial directions verifies isotropy.
        for axis in (0,1):
            result['tests'].append(axial(r,pressures[1],out,axis=axis));save(out/'summary.json',result)
        if layer in ['fractured_basalt','laboratory_reference_only']:
            result['tests'].append(axial(r,0.,out,tension=True));save(out/'summary.json',result)
        if layer=='sand':shear_pressures=[97.,250.,500.,680.]
        elif layer=='duricrust':shear_pressures=[500.,5000.]
        elif layer=='blocky_ejecta':shear_pressures=[5000.,50000.]
        else:shear_pressures=[]
        for p in shear_pressures:
            result['tests'].append(shear(r,p,out));save(out/'summary.json',result)
        if shear_pressures:
            result['tests'].append(shear(r,shear_pressures[-1],out,cyclic=True));save(out/'summary.json',result)
        result['passed']=all(t['passed'] for t in result['tests'])
    except Exception as e:result.update(error=str(e),traceback=traceback.format_exc())
    finally:ops.wipe()
    save(out/'summary.json',result);print(json.dumps(dict(id=name,passed=result['passed'],tests=len(result['tests']),error=result.get('error'))),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case');parser.add_argument('--timeout-s',type=float,default=45.);args=parser.parse_args()
    if args.case:worker(args.case);return
    output=HERE/'results';output.mkdir(exist_ok=True);all_results=[]
    for r in records():
        out=output/r['id'];out.mkdir(exist_ok=True);start=time.monotonic()
        with (out/'execution.log').open('w') as f:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--case',r['id']],stdout=f,stderr=f,timeout=args.timeout_s)
                result=json.loads((out/'summary.json').read_text());result['returncode']=p.returncode
            except subprocess.TimeoutExpired:
                result=json.loads((out/'summary.json').read_text());result.update(passed=False,timed_out=True,error='Native material case exceeded its explicit worker time limit; completed partial tests retained')
        result['elapsed_s']=time.monotonic()-start;save(out/'summary.json',result);all_results.append(result)
        save(output/'summary.json',dict(cases=all_results,passed=sum(x['passed'] for x in all_results),count=len(all_results),complete=len(all_results)==len(records())))
        print(r['id'],result['passed'],len(result['tests']),f'{result["elapsed_s"]:.2f}s',flush=True)


if __name__=='__main__':main()
