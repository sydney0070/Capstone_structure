"""Adapter from current FE results to the familiar offline material viewer.

No model assembly or synthetic displacement is performed here. Structure lines
are optional data supplied only by the integrated entry point.
"""
from pathlib import Path
import numpy as np
from insight_dynamic_viewer import save_dynamic_viewer

COLORS = {'mars_insight_duricrust_proxy':'#E6A23C',
          'mars_insight_sandy_regolith':'#63B5D6', 'mars_insight_blocky_ejecta':'#9B6451',
          'mars_insight_fractured_basalt':'#596579',
          'mars_jez1_peak':'#63B5D6', 'basalt_elastic_analogue':'#596579'}


def save_response_viewer(g, tags, coordinates, displacements, times, path, *,
                         context, structures=(), scale=None):
    tags = [int(t) for t in tags]
    xyz = np.asarray(coordinates, dtype=float)
    u = np.asarray(displacements, dtype=float)
    times = np.asarray(times, dtype=float)
    if (len(set(tags)) != len(tags) or xyz.shape != (len(tags),3)
            or u.shape != (len(times),len(tags),3) or not len(times)
            or not all(np.isfinite(a).all() for a in (xyz,u,times))
            or np.any(np.diff(times) <= 0)):
        raise ValueError('invalid aligned physical viewer frames')
    if not set(g['nodes']) <= set(tags):
        raise ValueError('record all physical ground nodes for the ground maximum')
    # The legacy exporter sorts tags: reorder all supplied arrays explicitly.
    order = np.argsort(tags); tags = [tags[i] for i in order]
    xyz = xyz[order]; u = u[:,order,:]
    index = {t:i for i,t in enumerate(tags)}
    lines = []
    for item in structures:
        lines.append(dict(name=item['name'], color=item['color'],
            segments=[[index[a],index[b]] for a,b in item['segments']]))
    names = list(dict.fromkeys(e['material'] for e in g['elements'].values()))
    palette = ['#c8754a','#a65132','#7c4633','#4b3028']
    materials = {name:dict(E=g['materials'][name]['E_Pa'],
                         rho=g['materials'][name]['density_kg_m3'],
                         law=g['materials'][name]['active_law'],
                         color=COLORS.get(name,palette[i % len(palette)]))
                 for i,name in enumerate(names)}
    if scale is None:
        delta = u if context.get('static_only') else u-u[0]
        peak = float(np.linalg.norm(delta,axis=2).max())
        span = max(float(np.ptp(xyz,axis=0).max()),1.)
        scale = span*.06/peak if peak > 1e-15 else 1.
    model = dict(node_coords={t:p.tolist() for t,p in zip(tags,xyz)},
        elements=g['elements'], materials=materials, material_order=names,
        body_name=g['config']['body'],
        dynamic_material_model=' · '.join(dict.fromkeys(m['law'] for m in materials.values())),
        viewer_context=context, viewer_structures=lines, save_viewer_payload=True)
    frames = [dict(time=float(t),disp=a.tolist()) for t,a in zip(times,u)]
    save_dynamic_viewer(model,frames,[],tags[0],path,visual_scale=scale)
    return dict(file=str(Path(path)),frames=len(times),nodes=len(tags),
                structures=len(lines),offline=True,reference='total' if context.get('static_only') else 'dynamic')
