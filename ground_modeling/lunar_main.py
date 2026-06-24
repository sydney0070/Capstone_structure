from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple

import openseespy.opensees as ops
import plotly.graph_objects as go

from lunar_config import ANA, GEO, LOAD, TagManager
from lunar_ground import LunarGround
from lunar_structure import LunarStructure

Vector3 = Tuple[float, float, float]


def get_output_dir() -> Path:
    """Return a result folder next to this script, independent of terminal cwd."""
    return Path(__file__).resolve().parent / ANA.output_dir


def initialize_domain() -> TagManager:
    """The only place where the OpenSees domain is initialized."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    return TagManager()


def configure_static_analysis(load_step: float = 0.02) -> None:
    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("NewtonLineSearch", 0.8)
    ops.integrator("LoadControl", load_step)
    ops.analysis("Static")


def run_gravity_and_pressure_analysis() -> None:
    configure_static_analysis(load_step=0.02)
    ok = ops.analyze(50)
    if ok == 0:
        return

    print("[WARN] NewtonLineSearch failed. Retrying with ModifiedNewton.")
    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-7, 200, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    ops.integrator("LoadControl", 0.005)
    ops.analysis("Static")
    ok = ops.analyze(200)
    if ok != 0:
        raise RuntimeError(f"Static analysis failed. OpenSees return code = {ok}")


def run_eigen_analysis(number_of_modes: int) -> List[float]:
    eigenvalues = ops.eigen("-fullGenLapack", number_of_modes)
    return [
        math.sqrt(value) / (2.0 * math.pi) if value > 0.0 else float("nan")
        for value in eigenvalues
    ]


def write_csv_results(
    structure: LunarStructure,
    ground: LunarGround,
    frequencies: List[float],
) -> Path:
    output_dir = get_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    with (output_dir / "nodal_displacements.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "node",
                "x_m",
                "y_m",
                "z_m",
                "ux_m",
                "uy_m",
                "uz_m",
                "rx_rad",
                "ry_rad",
                "rz_rad",
                "lumped_mass_kg",
            ]
        )
        for node in sorted(structure.coords):
            writer.writerow(
                [
                    node,
                    *structure.coords[node],
                    *ops.nodeDisp(node),
                    structure.node_mass[node],
                ]
            )

    ops.reactions()
    with (output_dir / "ground_reactions.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.writer(file)
        writer.writerow(
            ["ground_node", "Rx_N", "Ry_N", "Rz_N", "RMx_Nm", "RMy_Nm", "RMz_Nm"]
        )
        for node in ground.ground_nodes:
            writer.writerow([node, *ops.nodeReaction(node)])

    with (output_dir / "modal_frequencies.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.writer(file)
        writer.writerow(["mode", "frequency_Hz"])
        for mode, frequency in enumerate(frequencies, start=1):
            writer.writerow([mode, frequency])

    return output_dir


def write_html_visualization(
    structure: LunarStructure,
    ground: LunarGround,
    output_dir: Path,
) -> Path:
    """
    Create one interactive HTML containing both the structural model and
    the regolith/contact model.

    Important:
    - The contact springs and ground nodes are actual OpenSees model objects.
    - The translucent regolith block and its grid are visualization geometry.
      They do not represent SSPbrick elements in the current spring model.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    html_file = output_dir / "lunar_habitat_structure_and_ground.html"

    structural_nodes = sorted(structure.coords)
    displacement: Dict[int, Vector3] = {}
    max_translation = 0.0

    for node in structural_nodes:
        response = ops.nodeDisp(node)
        value = (float(response[0]), float(response[1]), float(response[2]))
        displacement[node] = value
        magnitude = math.sqrt(sum(component**2 for component in value))
        max_translation = max(max_translation, magnitude)

    if ANA.deformation_scale > 0.0:
        scale = ANA.deformation_scale
    elif max_translation > 1.0e-14:
        scale = 0.08 * max(GEO.length, 2.0 * GEO.radius) / max_translation
    else:
        scale = 1.0

    def deformed_coordinate(node: int) -> Vector3:
        x, y, z = structure.coords[node]
        ux, uy, uz = displacement[node]
        return x + scale * ux, y + scale * uy, z + scale * uz

    def frame_lines(
        use_deformed: bool,
    ) -> Tuple[List[float | None], List[float | None], List[float | None]]:
        xs: List[float | None] = []
        ys: List[float | None] = []
        zs: List[float | None] = []

        for element in structure.frame_elements:
            node_i, node_j, _, _ = structure.element_data[element]
            point_i = (
                deformed_coordinate(node_i)
                if use_deformed
                else structure.coords[node_i]
            )
            point_j = (
                deformed_coordinate(node_j)
                if use_deformed
                else structure.coords[node_j]
            )
            xs.extend([point_i[0], point_j[0], None])
            ys.extend([point_i[1], point_j[1], None])
            zs.extend([point_i[2], point_j[2], None])

        return xs, ys, zs

    ux_line, uy_line, uz_line = frame_lines(use_deformed=False)
    dx_line, dy_line, dz_line = frame_lines(use_deformed=True)

    undeformed_hover: List[str] = []
    deformed_hover: List[str] = []
    for node in structural_nodes:
        x, y, z = structure.coords[node]
        ux, uy, uz = displacement[node]
        magnitude = math.sqrt(ux**2 + uy**2 + uz**2)
        undeformed_hover.append(
            f"Structure node {node}<br>"
            f"x={x:.4f} m<br>y={y:.4f} m<br>z={z:.4f} m"
        )
        deformed_hover.append(
            f"Structure node {node}<br>"
            f"ux={ux:.4e} m<br>uy={uy:.4e} m<br>uz={uz:.4e} m<br>"
            f"|u|={magnitude:.4e} m"
        )

    # ------------------------------------------------------------
    # Actual contact model: base nodes, fixed ground nodes, springs
    # ------------------------------------------------------------
    contact_x: List[float | None] = []
    contact_y: List[float | None] = []
    contact_z: List[float | None] = []
    contact_hover: List[str | None] = []

    spring_turns = 5
    spring_amplitude = 0.08

    for index, (base_node, ground_node) in enumerate(
        zip(structure.base_nodes, ground.ground_nodes),
        start=1,
    ):
        bx, by, bz = structure.coords[base_node]
        gx, gy, gz = ground.coords[ground_node]

        # zeroLength elements have coincident coordinates. Offset the displayed
        # ground point downward and draw a coil symbol so the connection is visible.
        display_gz = gz - 0.30
        segments = 32
        for k in range(segments + 1):
            ratio = k / segments
            angle = 2.0 * math.pi * spring_turns * ratio
            x = bx + spring_amplitude * math.sin(angle)
            y = by + spring_amplitude * math.cos(angle)
            z = bz + (display_gz - bz) * ratio
            contact_x.append(x)
            contact_y.append(y)
            contact_z.append(z)
            contact_hover.append(
                f"Support {index}<br>"
                f"Contact element {ground.contact_elements[index - 1]}<br>"
                f"Base node {base_node}<br>Ground node {ground_node}"
            )

        contact_x.append(None)
        contact_y.append(None)
        contact_z.append(None)
        contact_hover.append(None)

    # ------------------------------------------------------------
    # Regolith visualization volume (not analysis brick elements)
    # ------------------------------------------------------------
    all_structure_x = [point[0] for point in structure.coords.values()]
    all_structure_y = [point[1] for point in structure.coords.values()]

    x_min = min(all_structure_x) - 1.5
    x_max = max(all_structure_x) + 1.5
    y_min = min(all_structure_y) - 1.5
    y_max = max(all_structure_y) + 1.5
    surface_z = min(point[2] for point in structure.coords.values())
    soil_depth = 2.5
    bottom_z = surface_z - soil_depth

    # Eight vertices of the regolith block.
    block_x = [x_min, x_max, x_max, x_min, x_min, x_max, x_max, x_min]
    block_y = [y_min, y_min, y_max, y_max, y_min, y_min, y_max, y_max]
    block_z = [
        surface_z, surface_z, surface_z, surface_z,
        bottom_z, bottom_z, bottom_z, bottom_z,
    ]

    # Triangular faces of the rectangular block.
    block_i = [0, 0, 4, 4, 0, 0, 1, 1, 2, 2, 3, 3]
    block_j = [1, 2, 5, 6, 1, 5, 2, 6, 3, 7, 0, 4]
    block_k = [2, 3, 6, 7, 5, 4, 6, 5, 7, 6, 4, 7]

    # 3D grid drawn through the block to make the soil volume visible.
    grid_x: List[float | None] = []
    grid_y: List[float | None] = []
    grid_z: List[float | None] = []

    nx_grid = 9
    ny_grid = 7
    nz_grid = 5

    x_values = [
        x_min + (x_max - x_min) * i / (nx_grid - 1)
        for i in range(nx_grid)
    ]
    y_values = [
        y_min + (y_max - y_min) * j / (ny_grid - 1)
        for j in range(ny_grid)
    ]
    z_values = [
        bottom_z + (surface_z - bottom_z) * k / (nz_grid - 1)
        for k in range(nz_grid)
    ]

    # Horizontal grids at several depths.
    for z_value in z_values:
        for x_value in x_values:
            grid_x.extend([x_value, x_value, None])
            grid_y.extend([y_min, y_max, None])
            grid_z.extend([z_value, z_value, None])
        for y_value in y_values:
            grid_x.extend([x_min, x_max, None])
            grid_y.extend([y_value, y_value, None])
            grid_z.extend([z_value, z_value, None])

    # Vertical grid lines.
    for x_value in x_values:
        for y_value in y_values:
            grid_x.extend([x_value, x_value, None])
            grid_y.extend([y_value, y_value, None])
            grid_z.extend([bottom_z, surface_z, None])

    # Foundation pad plates shown above each spring location.
    pad_x: List[float | None] = []
    pad_y: List[float | None] = []
    pad_z: List[float | None] = []
    half_pad = 0.30
    for base_node in structure.base_nodes:
        x, y, z = structure.coords[base_node]
        corners = [
            (x - half_pad, y - half_pad, z),
            (x + half_pad, y - half_pad, z),
            (x + half_pad, y + half_pad, z),
            (x - half_pad, y + half_pad, z),
            (x - half_pad, y - half_pad, z),
        ]
        for px, py, pz in corners:
            pad_x.append(px)
            pad_y.append(py)
            pad_z.append(pz)
        pad_x.append(None)
        pad_y.append(None)
        pad_z.append(None)

    figure = go.Figure()

    # Trace 0: undeformed structure
    figure.add_trace(
        go.Scatter3d(
            x=ux_line, y=uy_line, z=uz_line,
            mode="lines",
            name="Undeformed structure",
            line=dict(width=4),
            hoverinfo="skip",
            visible=True,
        )
    )

    # Trace 1: undeformed structural nodes
    figure.add_trace(
        go.Scatter3d(
            x=[structure.coords[node][0] for node in structural_nodes],
            y=[structure.coords[node][1] for node in structural_nodes],
            z=[structure.coords[node][2] for node in structural_nodes],
            mode="markers",
            name="Structural nodes",
            marker=dict(size=2.5),
            text=undeformed_hover,
            hovertemplate="%{text}<extra></extra>",
            visible=True,
        )
    )

    # Trace 2: deformed structure
    figure.add_trace(
        go.Scatter3d(
            x=dx_line, y=dy_line, z=dz_line,
            mode="lines",
            name=f"Deformed structure ×{scale:.3g}",
            line=dict(width=6),
            hoverinfo="skip",
            visible=True,
        )
    )

    # Trace 3: deformed nodes
    figure.add_trace(
        go.Scatter3d(
            x=[deformed_coordinate(node)[0] for node in structural_nodes],
            y=[deformed_coordinate(node)[1] for node in structural_nodes],
            z=[deformed_coordinate(node)[2] for node in structural_nodes],
            mode="markers",
            name=f"Deformed nodes ×{scale:.3g}",
            marker=dict(size=3.0),
            text=deformed_hover,
            hovertemplate="%{text}<extra></extra>",
            visible=True,
        )
    )

    # Trace 4: foundation pads
    figure.add_trace(
        go.Scatter3d(
            x=pad_x, y=pad_y, z=pad_z,
            mode="lines",
            name="Foundation pads",
            line=dict(width=7),
            hoverinfo="skip",
            visible=True,
        )
    )

    # Trace 5: base nodes
    figure.add_trace(
        go.Scatter3d(
            x=[structure.coords[node][0] for node in structure.base_nodes],
            y=[structure.coords[node][1] for node in structure.base_nodes],
            z=[structure.coords[node][2] for node in structure.base_nodes],
            mode="markers",
            name="Base nodes",
            marker=dict(size=7, symbol="diamond"),
            text=[f"Base node {node}" for node in structure.base_nodes],
            hovertemplate="%{text}<extra></extra>",
            visible=True,
        )
    )

    # Trace 6: fixed ground nodes, displayed below coincident base nodes
    figure.add_trace(
        go.Scatter3d(
            x=[ground.coords[node][0] for node in ground.ground_nodes],
            y=[ground.coords[node][1] for node in ground.ground_nodes],
            z=[ground.coords[node][2] - 0.30 for node in ground.ground_nodes],
            mode="markers",
            name="Fixed ground nodes",
            marker=dict(size=7, symbol="square"),
            text=[f"Fixed ground node {node}" for node in ground.ground_nodes],
            hovertemplate="%{text}<extra></extra>",
            visible=True,
        )
    )

    # Trace 7: contact springs
    figure.add_trace(
        go.Scatter3d(
            x=contact_x, y=contact_y, z=contact_z,
            mode="lines",
            name="Nonlinear ground-contact springs",
            line=dict(width=6),
            text=contact_hover,
            hovertemplate="%{text}<extra></extra>",
            visible=True,
        )
    )

    # Trace 8: translucent regolith volume
    figure.add_trace(
        go.Mesh3d(
            x=block_x,
            y=block_y,
            z=block_z,
            i=block_i,
            j=block_j,
            k=block_k,
            name="Regolith volume (visual aid)",
            opacity=0.28,
            flatshading=True,
            hovertemplate=(
                "Regolith visualization block<br>"
                "Current analysis model: discrete contact springs"
                "<extra></extra>"
            ),
            visible=True,
        )
    )

    # Trace 9: regolith grid
    figure.add_trace(
        go.Scatter3d(
            x=grid_x, y=grid_y, z=grid_z,
            mode="lines",
            name="Regolith visualization grid",
            line=dict(width=1),
            hoverinfo="skip",
            visible=True,
        )
    )

    overlay_title = (
        "3D Lunar Habitat + Regolith Foundation Model"
        f"<br><sup>Maximum translation: {max_translation:.4e} m · "
        f"Deformation magnification: ×{scale:.4g}</sup>"
    )

    # Visibility order corresponds to 10 traces above.
    show_all = [True] * 10
    undeformed_view = [
        True, True, False, False, True, True, True, True, True, True
    ]
    deformed_view = [
        False, False, True, True, True, True, True, True, True, True
    ]
    structure_only = [
        True, True, False, False, True, True, False, False, False, False
    ]
    ground_only = [
        False, False, False, False, True, True, True, True, True, True
    ]

    figure.update_layout(
        title=overlay_title,
        template="plotly_white",
        margin=dict(l=0, r=0, b=0, t=95),
        legend=dict(
            x=0.01,
            y=0.99,
            bgcolor="rgba(255,255,255,0.82)",
        ),
        scene=dict(
            aspectmode="data",
            xaxis_title="X — habitat axis (m)",
            yaxis_title="Y (m)",
            zaxis_title="Z (m)",
            camera=dict(eye=dict(x=1.55, y=1.55, z=1.10)),
        ),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.5,
                xanchor="center",
                y=1.08,
                yanchor="bottom",
                buttons=[
                    dict(
                        label="Structure + ground",
                        method="update",
                        args=[{"visible": show_all}, {"title": overlay_title}],
                    ),
                    dict(
                        label="Undeformed",
                        method="update",
                        args=[
                            {"visible": undeformed_view},
                            {"title": "Undeformed Structure + Regolith Model"},
                        ],
                    ),
                    dict(
                        label="Deformed",
                        method="update",
                        args=[
                            {"visible": deformed_view},
                            {
                                "title": (
                                    "Deformed Structure + Regolith Model"
                                    f"<br><sup>Magnification ×{scale:.4g}</sup>"
                                )
                            },
                        ],
                    ),
                    dict(
                        label="Structure only",
                        method="update",
                        args=[
                            {"visible": structure_only},
                            {"title": "Structural Model"},
                        ],
                    ),
                    dict(
                        label="Ground only",
                        method="update",
                        args=[
                            {"visible": ground_only},
                            {"title": "Regolith + Foundation Contact Model"},
                        ],
                    ),
                ],
            )
        ],
        annotations=[
            dict(
                text=(
                    "The soil block/grid is visualization geometry; "
                    "the current OpenSees soil model consists of contact springs."
                ),
                x=0.5,
                y=0.015,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="bottom",
                showarrow=False,
            )
        ],
    )

    figure.write_html(
        str(html_file),
        include_plotlyjs=True,
        full_html=True,
        auto_open=False,
        config={
            "responsive": True,
            "displaylogo": False,
            "scrollZoom": True,
        },
    )
    return html_file

def print_summary(
    structure: LunarStructure,
    ground: LunarGround,
    frequencies: List[float],
) -> None:
    habitat_nodes = [node for ring in structure.ring_nodes for node in ring]
    max_node = max(
        habitat_nodes,
        key=lambda node: math.sqrt(sum(value**2 for value in ops.nodeDisp(node)[:3])),
    )
    maximum = math.sqrt(sum(value**2 for value in ops.nodeDisp(max_node)[:3]))

    print("\n" + "=" * 66)
    print("MODULAR 3D LUNAR HABITAT + REGOLITH CONTACT MODEL")
    print("=" * 66)
    print(f"Structural nodes          : {len(structure.coords)}")
    print(f"Frame elements            : {len(structure.frame_elements)}")
    print(f"Ground contact elements   : {len(ground.contact_elements)}")
    print(f"Rotational spring elements: {len(ground.rotation_elements)}")
    print(f"Maximum translation       : {maximum:.6e} m at node {max_node}")
    if frequencies:
        print("Natural frequencies")
        for mode, frequency in enumerate(frequencies, start=1):
            print(f"  Mode {mode:02d}: {frequency:.6f} Hz")
    else:
        print("Natural frequencies      : not available")
    print("=" * 66)


def main() -> None:
    tags = initialize_domain()

    structure = LunarStructure(tags)
    structure.build()

    ground = LunarGround(tags, structure.base_nodes, structure.coords)
    ground.build(nonlinear_contact=ANA.use_nonlinear_contact)

    output_dir = get_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create the combined structure-ground HTML immediately.
    # Therefore, the geometry file remains available even if analysis fails.
    html_file = write_html_visualization(structure, ground, output_dir)
    print(f"\nGeometry HTML created before analysis:")
    print(f"  {html_file.resolve()}")

    frequencies: List[float] = []

    try:
        structure.apply_static_loads()
        run_gravity_and_pressure_analysis()

        try:
            frequencies = run_eigen_analysis(ANA.n_modes)
        except Exception as eigen_error:
            print(f"[WARN] Eigenvalue analysis failed: {eigen_error}")
            frequencies = []

        # Rewrite the same HTML after analysis so it contains deformed geometry.
        html_file = write_html_visualization(structure, ground, output_dir)

        if frequencies:
            write_csv_results(structure, ground, frequencies)
        else:
            # Still save displacement and reaction CSV files with no modal rows.
            write_csv_results(structure, ground, [])

        print_summary(structure, ground, frequencies)

    except Exception as analysis_error:
        print("\n[WARN] OpenSees analysis did not complete.")
        print(f"Reason: {analysis_error}")
        print("The undeformed structure-ground HTML was still created successfully.")

    print(f"\nResults folder:")
    print(f"  {output_dir.resolve()}")
    print(f"Combined HTML:")
    print(f"  {html_file.resolve()}")


if __name__ == "__main__":
    main()
