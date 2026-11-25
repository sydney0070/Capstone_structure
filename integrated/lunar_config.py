from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Geometry:
    length: float = 8.0
    radius: float = 2.0
    n_bays: int = 8
    n_ring_nodes: int = 16
    leg_height: float = 1.0


@dataclass(frozen=True)
class StructuralMaterial:
    E: float = 69.0e9
    nu: float = 0.33
    rho: float = 2700.0

    @property
    def G(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))


@dataclass(frozen=True)
class Section:
    A_ring: float = 2.0e-2
    Iy_ring: float = 2.0e-4
    Iz_ring: float = 2.0e-4
    J_ring: float = 4.0e-4

    A_long: float = 1.5e-2
    Iy_long: float = 1.2e-4
    Iz_long: float = 1.2e-4
    J_long: float = 2.4e-4

    A_leg: float = 3.0e-2
    Iy_leg: float = 3.0e-4
    Iz_leg: float = 3.0e-4
    J_leg: float = 6.0e-4


@dataclass(frozen=True)
class Loading:
    g_moon: float = 1.62
    internal_pressure: float = 55.0e3
    regolith_cover_thickness: float = 1.0
    regolith_density: float = 1500.0
    equipment_mass: float = 4000.0
    structural_mass_scale: float = 1.0


@dataclass(frozen=True)
class GroundContact:
    # zeroLengthImpact3D parameters
    initial_gap: float = 1.0e-3
    friction_ratio: float = 0.80
    tangential_penalty: float = 2.0e6
    normal_penalty_1: float = 8.0e6
    normal_penalty_2: float = 1.5e7
    yield_deformation: float = 2.0e-3
    cohesion_force: float = 300.0

    # Rotational resistance at each support
    krx: float = 1.5e6
    kry: float = 1.5e6
    krz: float = 5.0e5

    # Contact plane normal direction: 3 = global Z
    direction: int = 3


@dataclass(frozen=True)
class AnalysisConfig:
    n_modes: int = 8
    output_dir: str = "lunar_habitat_results"
    deformation_scale: float = 0.0
    use_nonlinear_contact: bool = True


class TagManager:
    """Single shared tag source for the entire OpenSees domain."""

    def __init__(self) -> None:
        self.node = 0
        self.element = 0
        self.material = 0
        self.transf = 0

    def next_node(self) -> int:
        self.node += 1
        return self.node

    def next_element(self) -> int:
        self.element += 1
        return self.element

    def next_material(self) -> int:
        self.material += 1
        return self.material

    def next_transf(self) -> int:
        self.transf += 1
        return self.transf


GEO = Geometry()
MAT = StructuralMaterial()
SEC = Section()
LOAD = Loading()
GROUND = GroundContact()
ANA = AnalysisConfig()
