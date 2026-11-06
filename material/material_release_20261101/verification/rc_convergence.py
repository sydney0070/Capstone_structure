"""Fixed-law RC verification; SI. No experimental calibration is claimed.
Numerical convergence includes complete prescribed path AND force-free unloading.
Run from release root: python verification/rc_convergence.py
"""
import csv, json, math, sys, time
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'rc_convergence'; OUT.mkdir(parents=True,exist_ok=True)
TARGETS=[.003,-.003,.006,-.006,.012,-.012,.018,-.018,0.]

def model(ny,nip,form='force',ne=1):
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
    ops.uniaxialMaterial('Concrete02',1,-30e6,-.002,-6e6,-.006,.1,3e6,1.5e9)
    ops.uniaxialMaterial('Steel02',2,400e6,200e9,.01,18,.925,.15)
    ops.section('Fiber',10)
    # Preserve the original gross concrete approximation to isolate numerics.
    ops.patch('rect',1,ny,4,-.3,-.2,.3,.2)
    for y in [-.245,.245]:ops.layer('straight',2,3,math.pi*.02**2/4,y,-.145,y,.145)
    for j in range(ne+1):ops.node(j+1,0,3*j/ne)
    ops.fix(1,1,1,1);ops.geomTransf('PDelta',1);ops.beamIntegration('Lobatto',1,10,nip)
    for j in range(ne):
        args=['forceBeamColumn',j+1,j+1,j+2,1,1,'-iter',30,1e-12] if form=='force' else ['dispBeamColumn',j+1,j+1,j+2,1,1]
        ops.element(*args)
    top=ne+1
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(top,1.,0.,0.)
    ops.constraints('Plain');ops.numberer('RCM');ops.system('BandGeneral')
    ops.test('NormUnbalance',.05,100,0);ops.algorithm('Newton');ops.integrator('LoadControl',0);ops.analysis('Static')
    return top

def snapshot(top,segment,fraction,phase):
    ops.reactions(); u,v,_=ops.nodeDisp(top);f=ops.getLoadFactor(1)
    rx,ry,m=ops.nodeReaction(1)
    return dict(phase=phase,segment=segment,fraction=fraction,u_m=u,v_m=v,force_N=f,rx_N=rx,ry_N=ry,mz_Nm=m,
                force_error_N=abs(rx+f),vertical_error_N=abs(ry),moment_error_Nm=abs(m-3*f),deformed_lever_diagnostic_Nm=abs(m-(3+v)*f))

def run(ny,nip,steps,form='force',ne=1,advance=None):
    name=f'{form}_ny{ny}_ip{nip}_steps{steps}_ne{ne}'
    top=model(ny,nip,form,ne); rows=[];failure=None;started=time.perf_counter()
    for seg,target in enumerate(TARGETS):
        begin=ops.nodeDisp(top,1)
        for k in range(1,steps+1):
            desired=begin+(target-begin)*k/steps
            if advance is None:
                ops.integrator('DisplacementControl',top,1,desired-ops.nodeDisp(top,1))
                accepted=ops.analyze(1)==0
            else:accepted=advance(top,desired)
            if not accepted:
                failure=f'path segment {seg} step {k}';break
            r=snapshot(top,seg,k/steps,'path'); rows.append(r)
            if not all(math.isfinite(x) for x in [r['force_N'],r['u_m'],r['mz_Nm']]):failure='nonfinite';break
            if abs(r['u_m']-desired)>1e-8:failure='displacement target';break
        if failure:break
    # The last prescribed displacement is zero, NOT the zero-force residual.
    end_force=ops.getLoadFactor(1)
    if not failure:
        for k in range(1,61):
            target_force=end_force*(1-k/60)
            ops.integrator('LoadControl',target_force-ops.getLoadFactor(1))
            if ops.analyze(1)!=0:failure=f'force-free unload {k}';break
            rows.append(snapshot(top,len(TARGETS),k/60,'unload'))
    path=[r for r in rows if r['phase']=='path']
    peak=max([abs(r['force_N']) for r in path] or [1.])
    info=dict(name=name,ny=ny,nip=nip,steps=steps,form=form,ne=ne,completed=failure is None,failure=failure,
              path_points=len(path),peak_force_N=peak,force_at_zero_displacement_N=end_force,
              residual_m=ops.nodeDisp(top,1) if not failure else None,
              final_force_N=ops.getLoadFactor(1),elapsed_s=time.perf_counter()-started)
    for key in ['force_error_N','vertical_error_N','moment_error_Nm']:
        info['max_'+key]=max([r[key] for r in rows] or [float('inf')])
    info['equilibrium_pass']=bool(not failure and info['max_force_error_N']/peak<1e-4 and info['max_vertical_error_N']/peak<1e-4 and info['max_moment_error_Nm']/(3*peak)<1e-4 and abs(info['final_force_N'])<.05)
    if rows:
        with (OUT/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print(name, json.dumps(info),flush=True)
    return info,rows

def compare(a,b):
    ia,ra=a;ib,rb=b
    if not (ia['completed'] and ib['completed']):return {'pass':False,'reason':'incomplete'}
    err=[]
    for seg in range(len(TARGETS)):
        x=[r for r in ra if r['phase']=='path' and r['segment']==seg]
        y=[r for r in rb if r['phase']=='path' and r['segment']==seg]
        grid=np.linspace(max(x[0]['fraction'],y[0]['fraction']),1,121)
        fa=np.interp(grid,[r['fraction'] for r in x],[r['force_N'] for r in x])
        fb=np.interp(grid,[r['fraction'] for r in y],[r['force_N'] for r in y])
        err.extend(fa-fb)
    scale=max(ia['peak_force_N'],ib['peak_force_N']);res=abs(ia['residual_m']-ib['residual_m'])
    d=dict(a=ia['name'],b=ib['name'],peak_relative=abs(ia['peak_force_N']-ib['peak_force_N'])/scale,
           path_rms_over_peak=float(np.sqrt(np.mean(np.square(err)))/scale),path_max_over_peak=float(np.max(np.abs(err))/scale),
           residual_difference_m=res)
    d['pass']=bool(ia['equilibrium_pass'] and ib['equilibrium_pass'] and d['peak_relative']<.02 and d['path_rms_over_peak']<.02 and d['path_max_over_peak']<.05 and res<.0002)
    return d

if __name__=='__main__':
    spec=[(24,5,60),(36,5,60),(48,5,60),(72,5,60),(48,7,60),(48,9,60),(48,10,60),(48,5,120),(48,7,120),
          (48,5,60,'disp',4),(48,5,60,'disp',8),(48,5,60,'disp',16),(48,5,60,'disp',32),(48,5,60,'disp',64)]
    result=dict(runtime=ops.version(),units='N m Pa s kg',scope='PDelta equilibrium uses initial lever3m; deformed-lever residual diagnostic only (PDelta is not finite-rotation Corotational). NO axial preload; 3 m cantilever; fixed uncalibrated C30/SD400 prototype; gross concrete overlap 0.7854%; no confinement',
                criteria=dict(peak_relative=.02,path_rms_over_peak=.02,path_max_over_peak=.05,residual_difference_m=.0002,equilibrium_relative=1e-4),cases=[],comparisons=[])
    runs=[]
    for s in spec:
        r=run(*s);runs.append(r);result['cases'].append(r[0]);(OUT/'summary.json').write_text(json.dumps(result,indent=2))
    for i,j in [(0,1),(1,2),(2,3),(2,4),(4,5),(5,6),(2,7),(4,8),(9,10),(10,11),(11,12),(12,13),(2,13)]:result['comparisons'].append(compare(runs[i],runs[j]))
    result['all_comparisons_pass']=all(d['pass'] for d in result['comparisons'])
    (OUT/'summary.json').write_text(json.dumps(result,indent=2))
