"""Owned, monitored ground-profile comparisons, without altering material values."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ENTRY=ROOT/'ground/ground_modeling/ground_dome_footing.py'


def save(path,data):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(path)


def one(profile,radial,out,timeout,native_timeout,footing_steps=32,algorithm='Newton',policy='primary_only',gravity_algorithm='KrylovNewton',static_steps=100,linear_system='UmfPack'):
    dest=out/(profile+'_r'+str(radial));dest.mkdir(parents=True,exist_ok=True)
    run=dest/'run'
    if run.exists() and any(run.iterdir()):raise ValueError('Fresh case required: '+str(run))
    cmd=[sys.executable,str(ENTRY),'--material-profile',profile,'--radial-subdivisions',str(radial),'--footing-steps',str(footing_steps),'--static-algorithm',algorithm,'--solver-policy',policy,'--html','--output-dir',str(run)]
    cmd+=['--static-steps',str(static_steps)]
    cmd+=['--linear-system',linear_system]
    if gravity_algorithm:cmd+=['--gravity-algorithm',gravity_algorithm]
    env=dict(os.environ);env.update({k:'1' for k in ['OPENBLAS_NUM_THREADS','OPENBLAS_DEFAULT_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']});env['PYTHONDONTWRITEBYTECODE']='1';env['MPLCONFIGDIR']='/private/tmp/capstone_specimen_mpl'
    result=dict(profile=profile,radial_subdivisions=radial,footing_steps=footing_steps,static_steps=static_steps,algorithm=algorithm,gravity_algorithm=gravity_algorithm or algorithm,solver_policy=policy,linear_system=linear_system,command=cmd,complete=False,entry_sha256=hashlib.sha256(ENTRY.read_bytes()).hexdigest(),started_on='2026-11-25',run_dir=str(run))
    (dest/'ground_dome_footing_snapshot.py').write_bytes(ENTRY.read_bytes())
    save(dest/'execution_status.json',result)
    with (dest/'execution.log').open('w') as f:
        started=time.monotonic();job=subprocess.Popen(cmd,stdout=f,stderr=f,env=env)
        last_change=started;signature=None;reason=None;checkpoint=None
        try:
            while job.poll() is None:
                time.sleep(2.)
                cp=run/'solver_checkpoint.json'
                if cp.exists():
                    try:checkpoint=json.loads(cp.read_text())
                    except (OSError,json.JSONDecodeError):continue
                    current=(checkpoint.get('event'),checkpoint.get('pattern'),checkpoint.get('before'),checkpoint.get('increment'),checkpoint.get('algorithm'))
                    if current!=signature:signature=current;last_change=time.monotonic()
                age=time.monotonic()-last_change
                if time.monotonic()-started>timeout:reason='case_wallclock_limit';break
                if checkpoint and checkpoint.get('event')=='attempt_started' and age>native_timeout:reason='native_attempt_no_return_limit';break
            if reason:
                job.terminate()
                try:job.wait(timeout=2.)
                except subprocess.TimeoutExpired:job.kill();job.wait()
            result.update(returncode=job.wait(),elapsed_s=time.monotonic()-started,limit_reason=reason,last_checkpoint=checkpoint)
        finally:
            if job.poll() is None:job.kill();job.wait()
    summary=run/'summary.json'
    if summary.exists():
        d=json.loads(summary.read_text());result.update(complete=bool(d.get('complete')) and result['returncode']==0,stage=d.get('stage'),error=d.get('error'))
        if result['complete']:
            result.update(max_settlement_mm=d['depressurized']['maximum_mm'],pressurized_max_mm=d['pressurized']['maximum_mm'])
    save(dest/'execution_status.json',result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profiles',nargs='+',required=True);p.add_argument('--radial',type=int,default=1);p.add_argument('--workers',type=int,default=2);p.add_argument('--timeout-s',type=float,default=1800.);p.add_argument('--native-timeout-s',type=float,default=180.);p.add_argument('--output',type=Path,required=True);p.add_argument('--footing-steps',type=int,default=32);p.add_argument('--algorithm',choices=['KrylovNewton','Newton','NewtonLineSearch'],default='Newton');p.add_argument('--policy',choices=['fallback','primary_only'],default='primary_only')
    p.add_argument('--gravity-algorithm',choices=['KrylovNewton','NewtonLineSearch','Newton'],default='KrylovNewton');p.add_argument('--static-steps',type=int,default=100)
    p.add_argument('--linear-system',choices=['UmfPack','SuperLU'],default='UmfPack')
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True);results=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(one,name,args.radial,args.output,args.timeout_s,args.native_timeout_s,args.footing_steps,args.algorithm,args.policy,args.gravity_algorithm,args.static_steps,args.linear_system):name for name in args.profiles}
        for future in as_completed(jobs):
            try:r=future.result()
            except Exception as e:r=dict(profile=jobs[future],complete=False,error=str(e))
            results.append(r);save(args.output/'batch_summary.json',dict(cases=results,finished=len(results),requested=len(args.profiles),complete=len(results)==len(args.profiles)))
            print(json.dumps(r,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
