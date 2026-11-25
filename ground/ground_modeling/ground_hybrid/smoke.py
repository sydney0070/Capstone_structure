"""End-to-end CLI checks of three planets, flat terrain and supplied DEMs.

Each example includes the raft, four frame columns, gravity and extra10kPa.
Results replace the named example folders under hybrid_results.
"""
import json
from pathlib import Path
import subprocess
import sys
import time
from .io import write_json

ROOT=Path(__file__).resolve().parents[1]


def main():
    verification=ROOT/'hybrid_verification';verification.mkdir(exist_ok=True)
    rows=[];target=verification/'smoke.json'
    write_json(target,dict(complete=False,all_pass=False,cases=[]))
    for terrain in ['flat','dem']:
        for body in ['Earth','Moon','Mars']:
            name=body if terrain=='flat' else body+'_DEM'
            output=ROOT/'hybrid_results'/name
            command=[sys.executable,'-m','ground_hybrid.bounded_run','--timeout-s','90',
                     '--body',body,'--demo-frame','--terrain',terrain,'--no-open','--output',str(output)]
            if terrain=='dem':command+=['--prepare-pad']
            start=time.monotonic()
            with (verification/(name+'_cli.log')).open('w') as log:
                result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)
            status=output/'RUN_STATUS.json'
            complete=(result.returncode==0 and status.exists() and json.loads(status.read_text()).get('complete',False))
            row=dict(body=body,terrain=terrain,returncode=result.returncode,elapsed_s=time.monotonic()-start,
                     output=str(output),**{'pass':complete})
            if complete:
                data=json.loads((output/'summary.json').read_text())
                row.update(elements=data['elements'],nodes=data['nodes'],
                           settlement_mm=1000*data['raft_response']['mean_additional_settlement_m'])
            rows.append(row)
            write_json(target,dict(complete=len(rows)==6,all_pass=all(r['pass'] for r in rows),cases=rows))
            print(name,'PASS' if complete else 'FAIL',round(row['elapsed_s'],2),flush=True)
    return 0 if all(r['pass'] for r in rows) else 1


if __name__=='__main__':raise SystemExit(main())
