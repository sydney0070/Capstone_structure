from __future__ import annotations

from typing import Dict, List, Tuple

import openseespy.opensees as ops

from lunar_config import GROUND, GroundContact, TagManager

Vector3 = Tuple[float, float, float]


class LunarGround:
    """Builds only the support-ground interface.

    It receives base nodes created by LunarStructure. It does not initialize
    the OpenSees model and does not run an analysis.
    """

    def __init__(
        self,
        tags: TagManager,
        base_nodes: List[int],
        structure_coords: Dict[int, Vector3],
        parameters: GroundContact = GROUND,
    ) -> None:
        self.tags = tags
        self.base_nodes = base_nodes
        self.structure_coords = structure_coords
        self.parameters = parameters

        self.ground_nodes: List[int] = []
        self.contact_elements: List[int] = []
        self.rotation_elements: List[int] = []
        self.coords: Dict[int, Vector3] = {}

    def build(self, nonlinear_contact: bool = True) -> None:
        for base_node in self.base_nodes:
            x, y, z = self.structure_coords[base_node]
            ground_node = self.tags.next_node()
            ops.node(ground_node, x, y, z)
            ops.fix(ground_node, 1, 1, 1, 1, 1, 1)

            self.ground_nodes.append(ground_node)
            self.coords[ground_node] = (x, y, z)

            if nonlinear_contact:
                self._add_nonlinear_contact(ground_node, base_node)
            else:
                self._add_linear_spring_contact(ground_node, base_node)

            self._add_rotational_springs(ground_node, base_node)

    def _add_nonlinear_contact(self, ground_node: int, base_node: int) -> None:
        p = self.parameters
        element_tag = self.tags.next_element()
        ops.element(
            "zeroLengthImpact3D",
            element_tag,
            ground_node,
            base_node,
            p.direction,
            p.initial_gap,
            p.friction_ratio,
            p.tangential_penalty,
            p.normal_penalty_1,
            p.normal_penalty_2,
            p.yield_deformation,
            p.cohesion_force,
        )
        self.contact_elements.append(element_tag)

    def _add_linear_spring_contact(self, ground_node: int, base_node: int) -> None:
        """Fallback/debug model using three translational elastic springs."""
        p = self.parameters
        stiffnesses = [p.tangential_penalty, p.tangential_penalty, p.normal_penalty_1]
        material_tags: List[int] = []
        for stiffness in stiffnesses:
            material_tag = self.tags.next_material()
            ops.uniaxialMaterial("Elastic", material_tag, stiffness)
            material_tags.append(material_tag)

        element_tag = self.tags.next_element()
        ops.element(
            "zeroLength",
            element_tag,
            ground_node,
            base_node,
            "-mat",
            *material_tags,
            "-dir",
            1,
            2,
            3,
        )
        self.contact_elements.append(element_tag)

    def _add_rotational_springs(self, ground_node: int, base_node: int) -> None:
        p = self.parameters
        material_tags: List[int] = []
        for stiffness in (p.krx, p.kry, p.krz):
            material_tag = self.tags.next_material()
            ops.uniaxialMaterial("Elastic", material_tag, stiffness)
            material_tags.append(material_tag)

        element_tag = self.tags.next_element()
        ops.element(
            "zeroLength",
            element_tag,
            ground_node,
            base_node,
            "-mat",
            *material_tags,
            "-dir",
            4,
            5,
            6,
        )
        self.rotation_elements.append(element_tag)
