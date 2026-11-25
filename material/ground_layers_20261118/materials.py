"""UNFINISHED DRAFT: not a released or calibrated material set.

Independent OpenSees materials for the non-sandy InSight reference layers.

Units: m, N, s, kg, Pa. The two deeper strengths are explicitly trial inputs.
Importing this module does not build a model or change OpenSees state.
"""

import copy
import json
import math
from pathlib import Path

import openseespy.opensees as ops

HERE = Path(__file__).resolve().parent
MATERIAL_IDS = {
    'duricrust': 'mars_insight_duricrust_proxy',
    'blocky_ejecta': 'mars_insight_blocky_ejecta_dp_trial',
    'fractured_basalt': 'mars_insight_fractured_basalt_dp_trial',
}


def material(layer):
    """Return a copy of the frozen, source-labelled material record."""
    records = json.loads((HERE / 'materials.json').read_text())['materials']
    return copy.deepcopy(records[MATERIAL_IDS[layer]])


def register(layer, tag):
    """Register one 3D, non-dilatant, perfectly plastic Drucker–Prager material."""
    if type(tag) is not int or tag < 1:
        raise ValueError('tag must be a positive integer')
    m = material(layer)
    phi = math.radians(m['phi_deg'])
    psi = math.radians(m['psi_deg'])
    friction = 2 * math.sqrt(2 / 3) * math.sin(phi) / (3 - math.sin(phi))
    dilation = 2 * math.sqrt(2 / 3) * math.sin(psi) / (3 - math.sin(psi))
    # Match the implemented DP invariant equation to MC triaxial compression.
    # See material_release_20261101/research/dp_mapping_note.md.
    sigma_y = 6 * m['c_Pa'] * math.cos(phi) / (3 - math.sin(phi))
    ops.nDMaterial(
        'DruckerPrager', tag, m['K_Pa'], m['G_Pa'], sigma_y,
        friction, dilation, 0., 0., 0., 0., 0., 1., m['density_kg_m3'], 101000.,
    )
    return m


def hoek_brown_sigma1(sigma3_Pa, recipe):
    """Compression-positive rock-mass envelope; Hoek et al. (2002), Eqs. 2–5.

    UCS, mi, GSI and D in this project's recipe are assumed scenario inputs.
    """
    ucs = recipe['intact_ucs_Pa']
    mi, gsi, disturbance = recipe['mi'], recipe['GSI'], recipe['D']
    if not (ucs > 0 and mi > 0 and 0 <= gsi <= 100 and 0 <= disturbance <= 1):
        raise ValueError('invalid Hoek–Brown scenario')
    if sigma3_Pa < 0:
        raise ValueError('this conversion supports compression only')
    mb = mi * math.exp((gsi - 100) / (28 - 14 * disturbance))
    s = math.exp((gsi - 100) / (9 - 3 * disturbance))
    a = .5 + (math.exp(-gsi / 15) - math.exp(-20 / 3)) / 6
    return sigma3_Pa + ucs * (mb * sigma3_Pa / ucs + s) ** a


def rock_secant_mc(recipe):
    """Endpoint secant over [0, max confinement], not the paper's mean fit.

    The concave HB envelope lies above this line inside the specified interval.
    This bound only concerns the triaxial-compression path, not other Lode angles.
    """
    limit = recipe['maximum_confinement_Pa']
    if limit <= 0:
        raise ValueError('positive confinement interval required')
    intercept = hoek_brown_sigma1(0., recipe)
    slope = (hoek_brown_sigma1(limit, recipe) - intercept) / limit
    phi = math.degrees(math.asin((slope - 1) / (slope + 1)))
    cohesion = intercept / (2 * math.sqrt(slope))
    return dict(phi_deg=phi, c_Pa=cohesion)
