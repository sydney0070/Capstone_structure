"""Stress and implemented DP plastic-strain diagnostics at all brick Gauss points."""
import csv
import math
from collections import defaultdict

import numpy as np
import openseespy.opensees as ops


def collect_material_state(ground, path, *, reference_norms=None, threshold=1e-8):
    rows=[];by_material=defaultdict(list);norms={}
    for element,e in ground['elements'].items():
        name=e['material'];law=ground['materials'][name]['active_law']
        for gp in range(1,9):
            s=np.asarray(ops.eleResponse(element,'material',gp,'stress'),dtype=float)
            if s.shape!=(6,) or not np.isfinite(s).all():
                raise RuntimeError(f'Invalid stress response at {element}/{gp}')
            p=-float(np.mean(s[:3]));dev=s.copy();dev[:3]+=p
            q=math.sqrt(1.5*float(np.dot(dev[:3],dev[:3])+2*np.dot(dev[3:],dev[3:])))
            plastic=None;delta=None
            if law=='DruckerPrager':
                state=np.asarray(ops.eleResponse(element,'material',gp,'state'),dtype=float)
                if state.shape!=(5,) or not np.isfinite(state).all():
                    raise RuntimeError(f'Invalid DP state response at {element}/{gp}')
                # OpenSees3.4.0 mState[4] is norm_ep, not a damage/failure flag.
                plastic=float(state[4]);norms[element,gp]=plastic
                if reference_norms is not None:
                    delta=plastic-reference_norms[element,gp]
            row=[element,gp,name,law,e['depth_m'],p,q,*s,plastic,delta]
            rows.append(row);by_material[name].append((element,p,q,plastic,delta))
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['element','gauss_point','material','law','center_depth_m','mean_compression_Pa','q_Pa','sx_Pa','sy_Pa','sz_Pa','sxy_Pa','syz_Pa','szx_Pa','DP_plastic_strain_norm','DP_norm_change_since_gravity']);w.writerows(rows)
    result={}
    for name,data in by_material.items():
        plastic=[x for x in data if x[3] is not None]
        result[name]=dict(elements=len({x[0] for x in data}),gauss_points=len(data),
                          max_mean_compression_Pa=max(x[1] for x in data),max_q_Pa=max(x[2] for x in data),
                          dp_state_recorded=bool(plastic))
        if plastic:
            result[name].update(plastic_strain_norm_threshold=threshold,max_plastic_strain_norm=max(x[3] for x in plastic),
                                gauss_points_with_plastic_norm_above_threshold=sum(x[3]>threshold for x in plastic),
                                elements_with_plastic_norm_above_threshold=len({x[0] for x in plastic if x[3]>threshold}))
            if reference_norms is not None:
                result[name].update(max_abs_norm_change_since_gravity=max(abs(x[4]) for x in plastic),
                                    gauss_points_with_changed_norm_since_gravity=sum(abs(x[4])>threshold for x in plastic))
    return dict(by_material=result,scope='Stress and DP plastic-strain norm diagnostics. Norm change is not a plastic-strain tensor increment, bearing-capacity factor or collapse check.'),norms
