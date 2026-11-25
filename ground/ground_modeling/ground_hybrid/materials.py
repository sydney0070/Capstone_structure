"""Versioned, source-traceable material records; no silent replacement laws."""
import copy
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT.parents[1] / 'material' / 'material_release_20261101'
if str(RELEASE) not in sys.path:
    sys.path.insert(0, str(RELEASE))
from capmat.materials import dp_constants

def catalog():
    directory = Path(__file__).parent / 'configs'
    base = json.loads((directory / 'materials.json').read_text())
    # Audited specimen records are additive; original/saved records cannot be overwritten.
    for filename in ('materials_recommended_20261124.json','materials_nonlinear_20261125.json'):
        specimen_file = directory / filename
        if specimen_file.exists():
            extra = json.loads(specimen_file.read_text())['materials']
            if set(extra) & set(base['materials']):
                raise ValueError('Specimen catalog must not replace an existing material id')
            base['materials'].update(extra)
    return base

def record(name, stiffness_factor=1.):
    if not math.isfinite(stiffness_factor) or stiffness_factor <= 0:
        raise ValueError('stiffness_factor must be positive')
    m = copy.deepcopy(catalog()['materials'][name])
    for key in ('E_Pa', 'G_Pa', 'K_Pa'):
        m[key] *= stiffness_factor
    return m

def register(ops, name, tag, *, elastic=False, stiffness_factor=1.):
    m = record(name, stiffness_factor)
    law = m['law']
    rho = m['density_kg_m3']
    if elastic or law == 'ElasticIsotropic':
        ops.nDMaterial('ElasticIsotropic', tag, m['E_Pa'], m['nu'], rho)
        law = 'ElasticIsotropic'
    elif law == 'DruckerPrager':
        p = dp_constants(m['phi_deg'], m['c_Pa'], m['psi_deg'])
        ops.nDMaterial('DruckerPrager', tag, m['K_Pa'], m['G_Pa'],
                       p['sigmaY'], p['rho'], p['rhoBar'], 0., 0., 0., 0., 0., 1., rho, 101000.)
    elif law == 'PressureDependMultiYield02':
        ops.nDMaterial(law, tag, 3, rho, m['G_Pa'], m['K_Pa'], m['phi_deg'],
                       m['peak_shear_strain'], m['p_ref_Pa'], m['pressure_exponent'],
                       m['PTAngle_deg'], m['contrac1'], m['contrac3'], m['dilat1'], m['dilat3'],
                       m['n_surfaces'], m['contrac2'], m['dilat2'], m['liquefac1'], m['liquefac2'],
                       m['e'], m['cs1'], m['cs2'], m['cs3'], m['pa_Pa'], m['numerical_c_Pa'])
    else:
        raise ValueError(f'Unsupported material law: {law}')
    return {**m, 'name': name, 'tag': tag, 'active_law': law, 'mass_in_material': True}
