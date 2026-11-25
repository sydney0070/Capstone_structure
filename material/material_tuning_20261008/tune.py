"""CLI for forward runs, parameter sweeps and history-preserving calibration."""
import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from copy import deepcopy
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'capstone-mpl'))
import numpy as np
import openseespy.opensees as ops
from presets import CASES,get_case,validate
from experiments import simulate

BASE=Path(__file__).resolve().parent

def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def units(c):
    if c['kind']=='viscous': return 'velocity (m/s)','force (kN)'
    if c['kind'] in ('spring','contact'): return 'displacement (m)','force (kN)'
    if c['kind'].startswith('soil'): return 'engineering shear strain','shear stress (kPa)'
    return 'axial strain','axial stress (kPa)'

def save(c,r,out,target=None):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    write_json(out/'config.json',c)
    with (out/'history.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(r.keys()); w.writerows(zip(*r.values()))
    meta=dict(case=c['name'],status=c.get('status','NOT_CALIBRATED'),
              opensees=ops.version(),python=sys.version,points=len(r['x']),
              units=units(c),y_min=float(np.min(r['y'])),y_max=float(np.max(r['y'])),
              physical_validation='NOT_EVALUATED_WITHOUT_INDEPENDENT_EXPERIMENT',
              history_sha256=hashlib.sha256((out/'history.csv').read_bytes()).hexdigest())
    write_json(out/'summary.json',meta)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4.8),layout='constrained')
    ax.plot(r['x'],r['y'],color='#156b82',lw=1.5,label='OpenSees')
    if target is not None: ax.plot(target['x'],target['y'],'k.',ms=2,label='Target')
    ax.set(xlabel=units(c)[0],ylabel=units(c)[1],title=c['name'])
    ax.grid(alpha=.2);ax.legend();fig.savefig(out/'response.png',dpi=160);plt.close(fig)

def read_csv(path):
    a=np.genfromtxt(path,delimiter=',',names=True,encoding='utf-8-sig')
    if a.ndim!=1 or len(a)<2 or not {'x','y'}.issubset(a.dtype.names or ()):
        raise ValueError('CSV needs >=2 rows with x,y headers')
    d={name:np.asarray(a[name],float) for name in a.dtype.names}
    if any(not np.isfinite(v).all() for v in d.values()):
        raise ValueError('CSV contains NaN/Inf/missing values')
    if abs(d['x'][0])>1e-12: raise ValueError('CSV must start at x=0 for a virgin specimen')
    return d

def fit(c,manifest_path,out,maxiter=30):
    """Bounded search, serial OpenSees calls; each trial creates new specimens.
    Fits continuous material params only. Density/mesh/damping cannot be inferred
    from one static stress-strain curve. Multiple tests can have different pressure.
    """
    from scipy.optimize import differential_evolution
    manifest_path=Path(manifest_path).resolve()
    m=json.loads(manifest_path.read_text())
    if m.get('units')!='kN,m,s,tonne,kPa' or m.get('case')!=c['name']:
        raise ValueError('manifest case/units mismatch')
    if not m.get('source') or not m.get('conditions_confirmed'):
        raise ValueError('Supply data source and confirm specimen/test compatibility in manifest')
    tests=[]
    for spec in m['tests']:
        cc=deepcopy(c);cc['test'].update(spec.get('test',{}))
        d=read_csv(manifest_path.parent/spec['csv'])
        channels=spec.get('channels',{'y':1.})
        if not channels or any(v<=0 for v in channels.values()):
            raise ValueError('positive channel weights required')
        scales={k:max(float(np.ptp(d[k])),float(np.max(np.abs(d[k]))),1e-12)
                for k in channels}
        role=spec.get('role','train')
        if role not in ('train','validation'): raise ValueError('unknown test role')
        tests.append((cc,d,channels,scales,role))
    if not any(v[-1]=='train' for v in tests): raise ValueError('at least one train test needed')
    keys=list(c['bounds'])
    disallowed={'rho','density','soil_type','qz_type','surfaces','c','Kn','thickness'}
    if not keys or any(k in disallowed or k not in c['params'] for k in keys):
        raise ValueError('bounds must select identifiable continuous material parameters')
    bounds=np.asarray([c['bounds'][k] for k in keys],float)
    if bounds.shape!=(len(keys),2) or np.any(bounds[:,0]>=bounds[:,1]):
        raise ValueError('bounds must be finite increasing [lower,upper] pairs')
    if not np.isfinite(bounds).all(): raise ValueError('nonfinite bounds')
    failures=[]; evaluations=0
    def evaluate(v,role='train',return_curves=False):
        nonlocal evaluations
        evaluations+=1;errors=[];curves=[]
        for cc,d,channels,scales,test_role in tests:
            if role!=test_role: continue
            cc=deepcopy(cc); cc['params'].update(zip(keys,map(float,v)))
            rr=simulate(cc,d['x'])
            errors.append(sum(w*np.mean(((rr[k]-d[k])/scales[k])**2)
                              for k,w in channels.items())/sum(channels.values()))
            curves.append((cc,rr,d))
        if not errors: return (None,[]) if return_curves else None
        score=float(np.sqrt(np.mean(errors)))
        return (score,curves) if return_curves else score
    def objective(z):
        v=bounds[:,0]+z*(bounds[:,1]-bounds[:,0])
        try: return evaluate(v)
        except (ValueError,RuntimeError,ops.OpenSeesError) as exc:
            if len(failures)<50: failures.append(str(exc))
            return 1.e6
    result=differential_evolution(objective,[(0.,1.)]*len(keys),seed=20261008,
            maxiter=maxiter,popsize=6,tol=1e-5,atol=1e-7,polish=True,workers=1)
    best=bounds[:,0]+result.x*(bounds[:,1]-bounds[:,0])
    if result.fun>=1.e6: raise RuntimeError('No successful calibration candidate')
    c=deepcopy(c);c['params'].update(zip(keys,map(float,best)))
    c['status']='FIT_TO_SYNTHETIC_ONLY' if m.get('synthetic',False) else 'FITTED_REQUIRES_INDEPENDENT_VALIDATION'
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    write_json(out/'best_config.json',c)
    write_json(out/'manifest_used.json',m)
    info=dict(parameters=dict(zip(keys,map(float,best))),train_nrmse=float(result.fun),
              optimizer_success=bool(result.success),optimizer_message=str(result.message),
              evaluations=evaluations,failed_trial_examples=failures,
              seed=20261008,opensees=ops.version(),synthetic=bool(m.get('synthetic',False)),
              validation_nrmse=None)
    for role in ('train','validation'):
        score,curves=evaluate(best,role,True)
        info[role+'_nrmse']=score
        for i,(cc,rr,d) in enumerate(curves):
            cc['status']=c['status']; save(cc,rr,out/f'{role}_{i+1}',d)
    write_json(out/'fit_report.json',info)
    return info

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['list','init','run','all','sweep','fit'])
    ap.add_argument('--case',choices=list(CASES));ap.add_argument('--config',type=Path)
    ap.add_argument('--out',type=Path,default=BASE/'results')
    ap.add_argument('--param');ap.add_argument('--values',nargs='+',type=float)
    ap.add_argument('--targets',type=Path);ap.add_argument('--maxiter',type=int,default=30)
    args=ap.parse_args(argv)
    if args.action=='list':
        for n,c in CASES.items(): print(f'{n:26} {c["material"]:30} {c["kind"]}')
        return
    if args.action=='init':
        args.out.mkdir(parents=True,exist_ok=True)
        for n in CASES: write_json(args.out/f'{n}.json',get_case(n))
        return
    if args.action=='all':
        failed=[]
        for n in CASES:
            c=get_case(n)
            try: save(c,simulate(c),args.out/n);print('OK',n,flush=True)
            except Exception as e: failed.append(dict(case=n,error=str(e)));print('FAILED',n,e)
        write_json(args.out/'failures.json',failed)
        if failed: raise SystemExit(2)
        return
    if args.config:
        c=json.loads(args.config.read_text())
    elif args.case: c=get_case(args.case)
    else: ap.error('--case or --config required')
    validate(c)
    if args.action=='run': save(c,simulate(c),args.out)
    elif args.action=='sweep':
        if args.param not in c['params'] or not args.values: ap.error('--param and --values required')
        for i,v in enumerate(args.values):
            cc=deepcopy(c);cc['params'][args.param]=v
            save(cc,simulate(cc),args.out/f'{i:02}_{args.param}_{v:g}')
    elif args.action=='fit':
        if not args.targets: ap.error('--targets manifest.json required')
        print(json.dumps(fit(c,args.targets,args.out,args.maxiter),ensure_ascii=False,indent=2))

if __name__=='__main__': main()
