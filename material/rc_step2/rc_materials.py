"""RC constituent materials for OpenSeesPy (SI units: N, m, Pa).

Not experimentally calibrated. Concrete C30 and SD400 are user-selected
nominal strengths; all other numeric inputs are explicitly provisional.
"""
from dataclasses import dataclass, asdict
import json
from pathlib import Path


@dataclass(frozen=True)
class ConcreteParams:
    fpc: float = -30e6       # Pa; user-selected 30 MPa, OpenSees negative compression
    epsc0: float = -0.002   # PROVISIONAL peak compressive strain
    fpcu: float = -6e6      # PROVISIONAL residual compressive strength
    epsU: float = -0.006    # PROVISIONAL crushing strain
    lam: float = 0.10       # PROVISIONAL unloading stiffness ratio
    ft: float = 3e6        # PROVISIONAL tensile strength
    Ets: float = 1.5e9     # PROVISIONAL tension-softening slope, Pa

    @property
    def initial_tangent(self):
        return 2.0 * self.fpc / self.epsc0

    def validate(self):
        assert self.fpc < 0 and self.epsc0 < 0
        assert self.fpcu <= 0 and self.epsU < self.epsc0
        assert 0 <= self.lam <= 1 and self.ft >= 0 and self.Ets > 0


@dataclass(frozen=True)
class SteelParams:
    Fy: float = 400e6       # Pa; user-selected SD400 nominal yield strength
    E0: float = 200e9       # PROVISIONAL steel elastic modulus
    b: float = 0.01         # PROVISIONAL strain hardening ratio
    R0: float = 18.0        # official recommended range 10~20; chosen 18 is PROVISIONAL
    cR1: float = 0.925      # OpenSees recommended transition parameter
    cR2: float = 0.15       # OpenSees recommended transition parameter

    def validate(self):
        assert self.Fy > 0 and self.E0 > 0 and 0 <= self.b < 1
        assert self.R0 > 0 and self.cR1 > 0 and self.cR2 > 0


def register_rc_materials(ops, concrete_tag=1, steel_tag=2,
                          concrete=None, steel=None):
    """Call after ops.wipe()/ops.model(); defines *uniaxial materials only*."""
    concrete = concrete or ConcreteParams()
    steel = steel or SteelParams()
    concrete.validate()
    steel.validate()
    ops.uniaxialMaterial('Concrete02', concrete_tag,
                        concrete.fpc, concrete.epsc0, concrete.fpcu,
                        concrete.epsU, concrete.lam, concrete.ft, concrete.Ets)
    ops.uniaxialMaterial('Steel02', steel_tag,
                        steel.Fy, steel.E0, steel.b,
                        steel.R0, steel.cR1, steel.cR2)
    return {'concrete_tag': concrete_tag, 'steel_tag': steel_tag,
            'concrete': concrete, 'steel': steel}


def write_parameters_json(path, concrete=None, steel=None):
    concrete = concrete or ConcreteParams()
    steel = steel or SteelParams()
    data = {
        'units': {'length': 'm', 'force': 'N', 'stress': 'Pa', 'strain': 'dimensionless'},
        'status': 'PROTOTYPE_ONLY_NOT_CALIBRATED',
        'user_selected': {'concrete_fc_MPa': 30, 'rebar_grade': 'SD400', 'rebar_Fy_MPa': 400},
        'concrete_Concrete02': asdict(concrete),
        'steel_Steel02': asdict(steel),
        'concrete_initial_tangent_Pa': concrete.initial_tangent,
        'documentation': {
            'Concrete02': 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Concrete02.html',
            'Steel02': 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Steel02.html'
        }
    }
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
