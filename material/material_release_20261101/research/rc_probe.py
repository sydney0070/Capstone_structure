"""Read-only-source solver diagnosis; no coefficients or archived outputs modified."""
import csv, importlib.util, json, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'rc_element_comparison'
sys.path.insert(0, str(SOURCE))
import compare_elements as c
import openseespy.opensees as ops

def main():
    out=ROOT/'research'/'rc_probe_results';out.mkdir(exist_ok=True)
    tests=[('NormDispIncr',1e-9,36,18,5),('NormUnbalance',1e-3,24,12,5),
           ('NormUnbalance',1e-3,36,18,5),('NormUnbalance',1e-3,48,24,5),
           ('NormUnbalance',1e-3,24,12,7),('NormUnbalance',1e-3,48,24,7)]
    summary={'opensees_version':ops.version(),'python':sys.version,'cases':{}}
    for criterion,tol,ny,nz,nip in tests:
        name=f'{criterion}_{ny}x{nz}_ip{nip}'
        def step(du,retry):
            ops.test(criterion,tol,100,0)
            ops.algorithm('Newton')
            ops.integrator('DisplacementControl',2,1,du)
            return ops.analyze(1)==0,'Newton'
        c.step=step
        start=time.perf_counter()
        rows,info=c.run(name,ny,nz,nip,'force',30,60,False)
        info['elapsed_s']=time.perf_counter()-start
        summary['cases'][name]=info
        with (out/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=c.FIELDS);w.writeheader();w.writerows(rows)
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        print(name,info['status'],info['points'],info['balance_rel_peak'],flush=True)

if __name__=='__main__': main()
