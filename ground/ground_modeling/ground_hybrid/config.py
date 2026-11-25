"""Explicit, editable reference scenarios. Depths/distribution are assumptions."""
import copy
import json
import math
from pathlib import Path
import numpy as np
from .materials import catalog


def default_config(body='Earth', *, earth_model='dp', mars_profile='nonlinear_20261125'):
    body = body.capitalize()
    if body not in ('Earth', 'Moon', 'Mars'):
        raise ValueError('body must be Earth, Moon or Mars')
    if earth_model not in ('dp','pdmy02'):raise ValueError('earth_model dp or pdmy02')
    if body == 'Mars':
        if mars_profile not in ('nonlinear_20261125','reference'):
            raise ValueError('mars_profile nonlinear_20261125 or reference')
        filename='Mars_reference_20261116.json' if mars_profile=='reference' else 'Mars.json'
        config_path = Path(__file__).parent / 'configs' / filename
        return {**json.loads(config_path.read_text()), 'earth_model': earth_model}
    layers = {
        'Earth': [(2., 'earth_sand_dr30_dry_reference'), (6., 'earth_sand_dr50_dry_reference'),
                  (10., 'earth_sand_dr75_dry_reference'), (None, 'basalt_elastic_analogue')],
        'Moon': [(.3, 'moon_jsc1a_17_peak'), (1., 'moon_jsc1a_18_peak'),
                 (8., 'moon_jsc1a_19_peak'), (None, 'basalt_elastic_analogue')],
    }[body]
    if body=='Earth' and earth_model=='dp':
        layers=[(d,m.replace('_dry_reference','_dp_reference')) for d,m in layers]
    return dict(schema=1, body=body, gravity_m_s2=catalog()['gravity_m_s2'][body],
                scenario='layered_with_soil_lens_and_basal_rock', units='m N s kg Pa',
                depth_m=12., width_m=6., length_m=6., core_width_m=2., core_length_m=2.,
                core_spacing_m=1., max_outer_spacing_m=3., max_vertical_spacing_m=2.,
                mesh_scope='Compact integration example, NOT a domain/mesh-converged foundation design. The off-centre lens is clipped by the domain.',
                layers=[dict(bottom_m=d, material=m) for d,m in layers],
                # Physical coordinates stay fixed when the mesh is refined.
                lenses=[dict(center_m=[2.5,0.,4.5], radii_m=[1.5,2.,1.5],
                             material={'Earth':'earth_sand_dr30_'+('dp_reference' if earth_model=='dp' else 'dry_reference'),
                                       'Moon':'moon_jsc1a_17_peak'}[body],
                             evidence='assumed soil/crust pocket; distribution and size not measured')],
                layer_geometry='surface_following', boundary='static_box',
                elastic=False, earth_model=earth_model, stiffness_factor=1., pdmy_elastic_seating_fraction=.001, tag_start=10000,
                max_elements=12000, max_aspect_ratio=40.,
                assumption='Synthetic dry reference deposit; layer depths/lens location NOT inferred from DEM or a borehole.',
                material_scope='Finite continuum regions; not grain-scale mixing. Rock/soil internal interfaces bonded.',
                site_scope={'Earth':'dry sand over basalt analogue; not Pohang site reconstruction',
                            'Moon':'JSC-1A density states with basalt analogue; not polar ice/regolith reconstruction'}[body])


def validate_config(config):
    c = copy.deepcopy(config)
    if c['body'] not in ('Earth','Moon','Mars'):raise ValueError('unsupported body')
    for k in ('gravity_m_s2','depth_m','width_m','length_m','core_width_m','core_length_m',
              'core_spacing_m','max_outer_spacing_m','max_vertical_spacing_m','stiffness_factor'):
        if not math.isfinite(c[k]) or c[k] <= 0: raise ValueError(f'positive {k} required')
    if c['core_width_m'] >= c['width_m'] or c['core_length_m'] >= c['length_m']:
        raise ValueError('core must lie inside outer domain')
    if c['layer_geometry'] != 'surface_following':
        raise ValueError('only explicitly surface-following layers implemented')
    if c['boundary'] != 'static_box':
        raise ValueError('build with static_box; transition explicitly after gravity for dynamics')
    if type(c['tag_start']) is not int or c['tag_start'] <= 0: raise ValueError('positive tag_start')
    if type(c['max_elements']) is not int or c['max_elements']<1:raise ValueError('positive integer max_elements')
    if not math.isfinite(c['max_aspect_ratio']) or c['max_aspect_ratio']<1:raise ValueError('max_aspect_ratio >=1 required')
    names = catalog()['materials']
    bottoms=[]
    for i, layer in enumerate(c['layers']):
        if layer['material'] not in names: raise ValueError('unknown material')
        d = layer['bottom_m']
        if d is None:
            if i != len(c['layers'])-1: raise ValueError('open layer must be last')
        else:
            if not math.isfinite(d) or d <= 0 or (bottoms and d <= bottoms[-1]):
                raise ValueError('layer depths must increase')
            bottoms.append(d)
    if not c['layers'] or c['layers'][-1]['bottom_m'] is not None: raise ValueError('last layer must be open')
    for lens in c['lenses']:
        if lens['material'] not in names: raise ValueError('unknown lens material')
        if len(lens['center_m']) != 3 or len(lens['radii_m']) != 3 or not np.all(np.isfinite(lens['center_m'])) or not all(math.isfinite(r) and r>0 for r in lens['radii_m']):
            raise ValueError('finite lens center, positive radii')
    return c


def graded_axis(width, core_width, core_spacing, outer_spacing):
    n = max(2, math.ceil(core_width / core_spacing))
    if n % 2: n += 1
    positive = list(np.linspace(0., core_width/2, n//2+1))
    step = core_width/n
    while positive[-1] < width/2-1e-10:
        step = min(step*1.5, outer_spacing)
        remaining = width/2-positive[-1]
        if remaining < 1.5*step: step = remaining
        positive.append(positive[-1]+step)
    return np.array([-x for x in positive[:0:-1]]+positive)


def depth_levels(c):
    points = sorted({0., c['depth_m'], min(1.,c['depth_m']),
                     *[l['bottom_m'] for l in c['layers'] if l['bottom_m'] is not None and l['bottom_m'] < c['depth_m']]})
    out = [0.]
    for a,b in zip(points[:-1],points[1:]):
        step = min(.5,c['max_vertical_spacing_m']) if a < 1. else c['max_vertical_spacing_m']
        n = max(1, math.ceil((b-a)/step))
        out.extend(np.linspace(a,b,n+1)[1:].tolist())
    return out


def assigned_material(c,x,y,depth):
    for lens in c['lenses']:
        distance = sum(((q-q0)/r)**2 for q,q0,r in zip((x,y,depth),lens['center_m'],lens['radii_m']))
        if distance <= 1.: return lens['material']
    return next(l['material'] for l in c['layers'] if l['bottom_m'] is None or depth < l['bottom_m'])
