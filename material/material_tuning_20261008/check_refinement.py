"""Compare prescribed-path increment sizes; tolerance is a project choice."""
from copy import deepcopy
import numpy as np
from presets import get_case
from experiments import simulate
from tune import BASE,write_json

if __name__=='__main__':
    rows=[]
    for name in ('sand_pdmy02','clay_pimy','moon_regolith_dp','mars_regolith_dp',
                 'sand_pm4','silt_pm4','shell_steel','shell_aluminum','interface_friction'):
        c=get_case(name)
        try:
            r=simulate(c);cc=deepcopy(c);cc['test']['substeps']*=2
            s=simulate(cc,r['x'])
            err=float(np.max(abs(r['y']-s['y']))/max(np.max(abs(s['y'])),1e-12))
            row=dict(case=name,substeps_base=c['test']['substeps'],
                substeps_fine=cc['test']['substeps'],relative_max_path_error=err,
                project_tolerance=.02,passed=bool(err<=.02))
            if 'eps_v' in r:
                row['max_eps_v_difference']=float(max(abs(r['eps_v']-s['eps_v'])))
            if 'ru_equivalent' in r:
                row['max_ru_equivalent_difference']=float(max(abs(r['ru_equivalent']-s['ru_equivalent'])))
        except Exception as exc:
            row=dict(case=name,passed=False,error=str(exc))
        rows.append(row)
        print(row,flush=True)
    write_json(BASE/'refinement_checks.json',rows)
    if not all(r['passed'] for r in rows):raise SystemExit(2)
