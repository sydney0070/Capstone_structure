"""Factor-isolation probe; source module stays unchanged on disk."""
import csv, inspect, json, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'rc_element_comparison'))
import compare_elements as c
import openseespy.opensees as ops

src=inspect.getsource(c.model)
out=ROOT/'research'/'rc_transform_results';out.mkdir(exist_ok=True)
summary={'version':ops.version(),'criteria':'NormUnbalance 0.05 (mixed N and Nm raw norm; separate balance gate still active)','cases':{}}
for transf in ['Linear','PDelta','Corotational']:
    for n in [24,36,48]:
        ns=c.__dict__
        exec(src.replace("'PDelta'",repr(transf)).replace('element_iter,1e-10','element_iter,1e-12'),ns)
        def step(du,retry):
            ops.test('NormUnbalance',0.05,100,0)
            ops.algorithm('Newton')
            ops.integrator('DisplacementControl',2,1,du)
            return ops.analyze(1)==0,'Newton'
        c.step=step
        name=f'{transf}_ny{n}_ip5'
        start=time.perf_counter()
        rows,info=c.run(name,n,4,5,'force',30,60,False)
        info['elapsed_s']=time.perf_counter()-start
        summary['cases'][name]=info
        with (out/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=c.FIELDS);w.writeheader();w.writerows(rows)
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        print(name,info['status'],info['points'],info['peak_force_kN'],info['balance_rel_peak'],flush=True)
