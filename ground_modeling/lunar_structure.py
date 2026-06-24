from __future__ import annotations

import math
from typing import Callable, Dict, List, Tuple

import openseespy.opensees as ops

from lunar_config import GEO, LOAD, MAT, SEC, TagManager

Vector3 = Tuple[float, float, float]


def vector_sub(a: Vector3, b: Vector3) -> Vector3:
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def norm(v: Vector3) -> float:
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def unit(v: Vector3) -> Vector3:
    length = norm(v)
    if length <= 1.0e-12:
        raise ValueError("Zero-length direction vector encountered.")
    return v[0] / length, v[1] / length, v[2] / length


def cross(a: Vector3, b: Vector3) -> Vector3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def choose_vecxz(point_i: Vector3, point_j: Vector3) -> Vector3:
    """Choose a stable local x-z plane vector for a 3D beam element."""
    local_x = unit(vector_sub(point_j, point_i))
    candidates = [
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
        (1.0, 0.0, 0.0),
    ]
    reference = min(
        candidates,
        key=lambda r: abs(
            local_x[0] * r[0]
            + local_x[1] * r[1]
            + local_x[2] * r[2]
        ),
    )
    local_y = unit(cross(reference, local_x))
    return unit(cross(local_x, local_y))


class LunarStructure:
    """Build the 3D habitat structure and DEM-based support nodes only.

    This class must not call ops.wipe(), ops.model(), or analysis commands.
    The OpenSees 3D/6-DOF domain must already exist before build() is called.
    """

    def __init__(
        self,
        tags: TagManager,
        terrain_z: Callable[[float, float], float] | None = None,
        x_origin: float = 0.0,
        y_origin: float = 0.0,
    ) -> None:
        self.tags = tags
        self.terrain_z = terrain_z or (lambda x, y: 0.0)
        self.x_origin = float(x_origin)
        self.y_origin = float(y_origin)

        self.coords: Dict[int, Vector3] = {}
        self.ring_nodes: List[List[int]] = []
        self.base_nodes: List[int] = []
        self.support_top_nodes: List[int] = []
        self.frame_elements: List[int] = []
        self.element_data: Dict[int, Tuple[int, int, float, float]] = {}
        self.node_mass: Dict[int, float] = {}
        self.pressure_force_by_node: Dict[int, Vector3] = {}

        # Filled by _prepare_terrain_placement().
        self.support_locations: List[Tuple[float, float]] = []
        self.support_elevations: List[float] = []
        self.terrain_reference_z = 0.0

    def build(self) -> None:
        self._prepare_terrain_placement()
        self._create_habitat_nodes()
        self._create_ring_elements()
        self._create_longitudinal_elements()
        self._create_support_legs_and_base_nodes()
        self._assemble_nodal_masses()
        self._apply_mass_to_domain()
        self._validate_built_structure()

    def _support_xy_locations(self) -> List[Tuple[float, float]]:
        """Return the exact four support XY locations."""
        dx = GEO.length / GEO.n_bays
        axial_ring_ids = [1, GEO.n_bays - 1]
        lower_angles = [
            5.0 * math.pi / 4.0,
            7.0 * math.pi / 4.0,
        ]

        locations: List[Tuple[float, float]] = []

        for ring_id in axial_ring_ids:
            x = self.x_origin + ring_id * dx

            for angle in lower_angles:
                j = round(
                    angle
                    / (2.0 * math.pi)
                    * GEO.n_ring_nodes
                ) % GEO.n_ring_nodes

                theta = (
                    2.0
                    * math.pi
                    * j
                    / GEO.n_ring_nodes
                )

                y = (
                    self.y_origin
                    + GEO.radius * math.cos(theta)
                )

                locations.append((float(x), float(y)))

        if len(locations) != 4:
            raise RuntimeError(
                f"Exactly four support locations were expected, got {len(locations)}."
            )

        return locations

    def _prepare_terrain_placement(self) -> None:
        """Sample the DEM once and store the four support elevations."""
        self.support_locations = self._support_xy_locations()
        self.support_elevations = []

        for support_index, (x, y) in enumerate(
            self.support_locations,
            start=1,
        ):
            z = float(self.terrain_z(x, y))

            if not math.isfinite(z):
                raise ValueError(
                    f"Support {support_index} received invalid DEM elevation: "
                    f"x={x}, y={y}, z={z}"
                )

            self.support_elevations.append(z)

        # The cylindrical frame is placed relative to the highest support.
        # Its lowest ring point is reference_z + GEO.leg_height.
        self.terrain_reference_z = max(self.support_elevations)

        print("[LunarStructure] DEM support placement")
        for index, ((x, y), z) in enumerate(
            zip(self.support_locations, self.support_elevations),
            start=1,
        ):
            print(
                f"  Support {index}: "
                f"x={x:.4f} m, y={y:.4f} m, z={z:.4f} m"
            )
        print(
            "  Structure terrain reference z:",
            f"{self.terrain_reference_z:.4f} m",
        )

    def add_node(self, x: float, y: float, z: float) -> int:
        tag = self.tags.next_node()
        ops.node(tag, float(x), float(y), float(z))
        self.coords[tag] = (float(x), float(y), float(z))
        self.node_mass[tag] = 0.0
        return tag

    def _create_habitat_nodes(self) -> None:
        dx = GEO.length / GEO.n_bays

        for i in range(GEO.n_bays + 1):
            x = self.x_origin + i * dx
            ring: List[int] = []

            for j in range(GEO.n_ring_nodes):
                theta = (
                    2.0
                    * math.pi
                    * j
                    / GEO.n_ring_nodes
                )

                y = (
                    self.y_origin
                    + GEO.radius * math.cos(theta)
                )

                z = (
                    self.terrain_reference_z
                    + GEO.radius * math.sin(theta)
                    + GEO.radius
                    + GEO.leg_height
                )

                ring.append(self.add_node(x, y, z))

            self.ring_nodes.append(ring)

    def _add_elastic_beam(
        self,
        node_i: int,
        node_j: int,
        area: float,
        iy: float,
        iz: float,
        torsion_j: float,
    ) -> int:
        point_i = self.coords[node_i]
        point_j = self.coords[node_j]

        element_length = norm(vector_sub(point_j, point_i))
        if element_length <= 1.0e-9:
            raise ValueError(
                f"Zero-length beam requested between nodes "
                f"{node_i} and {node_j}."
            )

        transf_tag = self.tags.next_transf()
        vecxz = choose_vecxz(point_i, point_j)
        ops.geomTransf("Linear", transf_tag, *vecxz)

        element_tag = self.tags.next_element()
        ops.element(
            "elasticBeamColumn",
            element_tag,
            node_i,
            node_j,
            area,
            MAT.E,
            MAT.G,
            torsion_j,
            iy,
            iz,
            transf_tag,
        )

        self.frame_elements.append(element_tag)
        self.element_data[element_tag] = (
            node_i,
            node_j,
            area,
            MAT.rho,
        )
        return element_tag

    def _create_ring_elements(self) -> None:
        for ring in self.ring_nodes:
            for j in range(GEO.n_ring_nodes):
                self._add_elastic_beam(
                    ring[j],
                    ring[(j + 1) % GEO.n_ring_nodes],
                    SEC.A_ring,
                    SEC.Iy_ring,
                    SEC.Iz_ring,
                    SEC.J_ring,
                )

    def _create_longitudinal_elements(self) -> None:
        for i in range(GEO.n_bays):
            for j in range(GEO.n_ring_nodes):
                self._add_elastic_beam(
                    self.ring_nodes[i][j],
                    self.ring_nodes[i + 1][j],
                    SEC.A_long,
                    SEC.Iy_long,
                    SEC.Iz_long,
                    SEC.J_long,
                )

    def _create_support_legs_and_base_nodes(self) -> None:
        axial_ring_ids = [1, GEO.n_bays - 1]
        lower_angles = [
            5.0 * math.pi / 4.0,
            7.0 * math.pi / 4.0,
        ]

        support_index = 0

        for ring_id in axial_ring_ids:
            for angle in lower_angles:
                j = round(
                    angle
                    / (2.0 * math.pi)
                    * GEO.n_ring_nodes
                ) % GEO.n_ring_nodes

                top_node = self.ring_nodes[ring_id][j]
                top_x, top_y, top_z = self.coords[top_node]

                # Use the already sampled elevation. Do not sample the DEM again.
                expected_x, expected_y = self.support_locations[support_index]
                base_z = self.support_elevations[support_index]

                if (
                    abs(top_x - expected_x) > 1.0e-8
                    or abs(top_y - expected_y) > 1.0e-8
                ):
                    raise RuntimeError(
                        "Support coordinate mismatch between the prepared DEM "
                        "location and the actual structural ring node."
                    )

                if top_z <= base_z:
                    raise ValueError(
                        f"Support {support_index + 1} has non-positive leg length: "
                        f"top_z={top_z}, base_z={base_z}"
                    )

                base_node = self.add_node(
                    top_x,
                    top_y,
                    base_z,
                )

                self.support_top_nodes.append(top_node)
                self.base_nodes.append(base_node)

                self._add_elastic_beam(
                    top_node,
                    base_node,
                    SEC.A_leg,
                    SEC.Iy_leg,
                    SEC.Iz_leg,
                    SEC.J_leg,
                )

                support_index += 1

    def _assemble_nodal_masses(self) -> None:
        for node_i, node_j, area, density in self.element_data.values():
            element_length = norm(
                vector_sub(
                    self.coords[node_j],
                    self.coords[node_i],
                )
            )
            element_mass = (
                area
                * element_length
                * density
                * LOAD.structural_mass_scale
            )

            self.node_mass[node_i] += 0.5 * element_mass
            self.node_mass[node_j] += 0.5 * element_mass

        habitat_nodes = [
            node
            for ring in self.ring_nodes
            for node in ring
        ]

        equipment_per_node = (
            LOAD.equipment_mass / len(habitat_nodes)
        )

        for node in habitat_nodes:
            self.node_mass[node] += equipment_per_node

        dx = GEO.length / GEO.n_bays
        arc = (
            2.0
            * math.pi
            * GEO.radius
            / GEO.n_ring_nodes
        )

        for i, ring in enumerate(self.ring_nodes):
            axial_tributary = (
                0.5 * dx
                if i in (0, GEO.n_bays)
                else dx
            )
            area_tributary = axial_tributary * arc
            regolith_mass = (
                LOAD.regolith_density
                * LOAD.regolith_cover_thickness
                * area_tributary
            )

            for node in ring:
                self.node_mass[node] += regolith_mass

    def _apply_mass_to_domain(self) -> None:
        for node, mass in self.node_mass.items():
            rotational_mass = max(
                mass * GEO.radius**2 * 1.0e-6,
                1.0e-9,
            )

            ops.mass(
                node,
                mass,
                mass,
                mass,
                rotational_mass,
                rotational_mass,
                rotational_mass,
            )

    def _validate_built_structure(self) -> None:
        expected_ring_nodes = (
            (GEO.n_bays + 1)
            * GEO.n_ring_nodes
        )

        actual_ring_nodes = sum(
            len(ring) for ring in self.ring_nodes
        )

        if actual_ring_nodes != expected_ring_nodes:
            raise RuntimeError(
                f"Ring-node count mismatch: expected "
                f"{expected_ring_nodes}, got {actual_ring_nodes}."
            )

        if len(self.base_nodes) != 4:
            raise RuntimeError(
                f"Expected 4 base nodes, got {len(self.base_nodes)}."
            )

        if len(self.frame_elements) == 0:
            raise RuntimeError("No frame elements were created.")

        print("[LunarStructure] Build completed")
        print(f"  Structural nodes : {len(self.coords)}")
        print(f"  Frame elements   : {len(self.frame_elements)}")
        print(f"  Base nodes       : {self.base_nodes}")

    def apply_static_loads(
        self,
        time_series_tag: int = 1,
        pattern_tag: int = 1,
    ) -> None:
        ops.timeSeries("Linear", time_series_tag)
        ops.pattern("Plain", pattern_tag, time_series_tag)

        for node, mass in self.node_mass.items():
            ops.load(
                node,
                0.0,
                0.0,
                -mass * LOAD.g_moon,
                0.0,
                0.0,
                0.0,
            )

        dx = GEO.length / GEO.n_bays
        arc = (
            2.0
            * math.pi
            * GEO.radius
            / GEO.n_ring_nodes
        )

        for i, ring in enumerate(self.ring_nodes):
            axial_tributary = (
                0.5 * dx
                if i in (0, GEO.n_bays)
                else dx
            )
            tributary_area = axial_tributary * arc
            radial_force = (
                LOAD.internal_pressure
                * tributary_area
            )

            for j, node in enumerate(ring):
                theta = (
                    2.0
                    * math.pi
                    * j
                    / GEO.n_ring_nodes
                )

                fy = radial_force * math.cos(theta)
                fz = radial_force * math.sin(theta)

                ops.load(
                    node,
                    0.0,
                    fy,
                    fz,
                    0.0,
                    0.0,
                    0.0,
                )

                self.pressure_force_by_node[node] = (
                    0.0,
                    fy,
                    fz,
                )
