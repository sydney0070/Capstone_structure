"""Model sensitivities, reported even when a recommended numerical gate fails."""
import argparse
import json
import subprocess
import sys
import time
import traceback
from pathlib import Path
from .verification import assembled
from .io import write_json

OUT=Path(__file__).resolve().parents[1]/'hybrid_verification'
VARIANTS=dict(baseline={},stiffness_half=dict(factor=.5),stiffness_double=dict(factor=2.),
              penalty_half=dict(penalty=5.),penalty_double=dict(penalty=20.),
              horizontal_refined=dict(spacing=.5),wider_domain=dict(width=12.))


def comparisons_for(cases):
    comparisons=[]
    for body in ['Earth','Moon','Mars']:
        rows={r['variant']:r for r in cases if r['body']==body}
        if 'baseline' not in rows:continue
        for name in ['penalty_half','penalty_double','horizontal_refined','wider_domain']:
            if name not in rows:continue
            limit=.05 if name.startswith('penalty') else .10
            valid=rows['baseline']['pass'] and rows[name]['pass']
            change=(abs(rows[name]['additional_settlement_m']-rows['baseline']['additional_settlement_m'])/
                    rows['baseline']['additional_settlement_m']) if valid else None
            comparisons.append(dict(body=body,variant=name,relative_settlement_change=change,
                                     recommended_limit=limit,analysis_complete=bool(valid),
                                     **{'pass':bool(valid and change<limit)}))
    return comparisons

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--timeout-s',type=float,default=120.)
    p.add_argument('--case-body',choices=['Earth','Moon','Mars'])
    p.add_argument('--case-variant',choices=list(VARIANTS))
    args=p.parse_args(argv)
    if not 0<args.timeout_s<float('inf'):p.error('finite positive timeout required')
    OUT.mkdir(exist_ok=True);cases=[];comparisons=[];start=time.time()
    if args.case_body:
        if args.case_variant is None:p.error('--case-variant required')
        try:r=assembled(args.case_body,**VARIANTS[args.case_variant])
        except Exception as e:r=dict(body=args.case_body,**{'pass':False},error=str(e),traceback=traceback.format_exc())
        r['variant']=args.case_variant
        write_json(OUT/f'sensitivity_{args.case_body}_{args.case_variant}.json',r)
        return 0 if r['pass'] else 1
    for body in ['Earth','Moon','Mars']:
        for name in VARIANTS:
            result_path=OUT/f'sensitivity_{body}_{name}.json'
            result_path.unlink(missing_ok=True)
            with (OUT/f'sensitivity_{body}_{name}.log').open('w') as log:
                try:
                    child=subprocess.run([sys.executable,'-m',__spec__.name,'--case-body',body,'--case-variant',name],
                                         stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout_s)
                    if not result_path.exists():raise RuntimeError(f'child exited {child.returncode} without result')
                    r=json.loads(result_path.read_text())
                except (subprocess.TimeoutExpired,RuntimeError) as e:
                    r=dict(body=body,variant=name,**{'pass':False},error=str(e),timeout_s=args.timeout_s)
            cases.append(r);comparisons=comparisons_for(cases)
            write_json(OUT/'sensitivity.json',dict(cases=cases,comparisons=comparisons,elapsed_s=time.time()-start,complete=False))
            print(body,name,r.get('additional_settlement_m'),r['pass'],flush=True)
    write_json(OUT/'sensitivity.json',dict(cases=cases,comparisons=comparisons,elapsed_s=time.time()-start,complete=True,
       all_analyses_pass=all(r['pass'] for r in cases),numerical_sensitivity_pass=len(comparisons)==12 and all(r['pass'] for r in comparisons),
       scope='stiffness factors are physical uncertainty, not a convergence pass target; mesh/domain tolerances are preliminary project gates, not standards'))
    return 0 if all(r['pass'] for r in cases) and all(r['pass'] for r in comparisons) else 1
if __name__=='__main__':raise SystemExit(main())
