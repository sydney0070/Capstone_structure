# ============================================================
# PressureDependMultiYield02
# 3D Single Element Cyclic Simple Shear Test
#
# + 3D deformation animation
# + tau-gamma hysteresis animation
#
# Units
#   Length : m
#   Force  : kN
#   Stress : kPa = kN/m^2
#   Mass   : ton
# ============================================================

import openseespy.opensees as ops

import numpy as np

import argparse
from pathlib import Path

from material_evaluation import evaluate, save_report


# ============================================================
# 1. USER SETTINGS
# ============================================================

# ------------------------------------------------------------
# Element geometry
# ------------------------------------------------------------

Lx = 1.0
Ly = 1.0
H = 1.0

AREA = Lx * Ly


# ------------------------------------------------------------
# Initial vertical stress
# ------------------------------------------------------------

sigma_v = 100.0       # kPa

total_vertical_load = sigma_v * AREA


# ------------------------------------------------------------
# Cyclic shear
# ------------------------------------------------------------

gamma_max = 0.01      # 1 %

n_cycles = 3

steps_per_quarter_cycle = 40


# ------------------------------------------------------------
# Visualization
#
# Actual deformation is very small.
# This factor ONLY magnifies the plotted deformation.
# ------------------------------------------------------------

VIS_SCALE = 15.0

ANIMATION_INTERVAL = 30       # ms

SHOW_ANIMATION = True

SAVE_ANIMATION = False

ANIMATION_FILENAME = "PDMY02_cyclic_shear.gif"


# ============================================================
# 2. PressureDependMultiYield02 PARAMETERS
# ============================================================

MAT_TAG = 1

nd = 3


# ------------------------------------------------------------
# Density
# ------------------------------------------------------------

rho = 1.9              # ton/m^3


# ------------------------------------------------------------
# Reference stiffness
# ------------------------------------------------------------

G_ref = 1.0e5          # kPa

B_ref = 2.33e5         # kPa


# ------------------------------------------------------------
# Strength
# ------------------------------------------------------------

friction_angle = 33.5

peak_shear_strain = 0.10


# ------------------------------------------------------------
# Pressure dependency
# ------------------------------------------------------------

ref_pressure = 101.0

pressure_depend_coeff = 0.5


# ------------------------------------------------------------
# Phase transformation
# ------------------------------------------------------------

PT_angle = 25.5


# ------------------------------------------------------------
# Contraction
# ------------------------------------------------------------

contrac1 = 0.045

contrac3 = 0.15


# ------------------------------------------------------------
# Dilation
# ------------------------------------------------------------

dilat1 = 0.06

dilat3 = 0.15


# ------------------------------------------------------------
# Yield surfaces
# ------------------------------------------------------------

n_yield_surface = 20


def run_test(amplitude=gamma_max, quarters=steps_per_quarter_cycle, cycles=n_cycles):
    """Independent specimen; original material/BCs, no plotting or evaluation."""
    if amplitude <= 0 or quarters < 2 or cycles < 1:
        raise ValueError("amplitude > 0, quarters >= 2, cycles >= 1 required")
    gamma_max = amplitude
    steps_per_quarter_cycle = quarters
    n_cycles = cycles
    # ============================================================
    # 3. INITIALIZE OPENSEES
    # ============================================================

    ops.wipe()

    ops.model(
        "basic",
        "-ndm", 3,
        "-ndf", 3
    )


    # ============================================================
    # 4. ORIGINAL NODE COORDINATES
    # ============================================================
    #
    #         z
    #         ↑
    #
    #       8 -------- 7
    #      /|         /|
    #     5 -------- 6 |
    #     | |        | |
    #     | 4 -------|-3
    #     |/         |/
    #     1 -------- 2  → x
    #
    # y direction = depth
    #

    original_coordinates = {

        1: np.array([0.0, 0.0, 0.0]),

        2: np.array([Lx, 0.0, 0.0]),

        3: np.array([Lx, Ly, 0.0]),

        4: np.array([0.0, Ly, 0.0]),

        5: np.array([0.0, 0.0, H]),

        6: np.array([Lx, 0.0, H]),

        7: np.array([Lx, Ly, H]),

        8: np.array([0.0, Ly, H]),
    }


    # ============================================================
    # 5. CREATE NODES
    # ============================================================

    for node_tag, coord in original_coordinates.items():

        ops.node(
            node_tag,
            coord[0],
            coord[1],
            coord[2]
        )


    # ============================================================
    # 6. BOUNDARY CONDITIONS
    # ============================================================
    #
    # Bottom:
    #
    # ux = 0
    # uy = 0
    # uz = 0
    #
    #
    # Top:
    #
    # ux = free
    # uy = fixed
    # uz = free
    #
    # The top four nodes are tied in x and z.
    # ============================================================


    # Bottom nodes
    for node in [1, 2, 3, 4]:

        ops.fix(
            node,
            1,
            1,
            1
        )


    # Top nodes
    for node in [5, 6, 7, 8]:

        ops.fix(
            node,
            0,
            1,
            0
        )


    # ------------------------------------------------------------
    # Make top surface move together
    #
    # Node 5 = master
    # ------------------------------------------------------------

    for node in [6, 7, 8]:

        ops.equalDOF(
            5,
            node,
            1,
            3
        )


    # ============================================================
    # 7. MATERIAL
    # ============================================================

    ops.nDMaterial(

        "PressureDependMultiYield02",

        MAT_TAG,

        nd,

        rho,

        G_ref,

        B_ref,

        friction_angle,

        peak_shear_strain,

        ref_pressure,

        pressure_depend_coeff,

        PT_angle,

        contrac1,

        contrac3,

        dilat1,

        dilat3,

        n_yield_surface
    )


    # ============================================================
    # 8. SSPbrick ELEMENT
    # ============================================================

    ELE_TAG = 1


    ops.element(

        "SSPbrick",

        ELE_TAG,

        1, 2, 3, 4,

        5, 6, 7, 8,

        MAT_TAG,

        0.0,
        0.0,
        0.0
    )


    # ============================================================
    # 9. INITIAL ELASTIC STAGE
    # ============================================================

    ops.updateMaterialStage(

        "-material",
        MAT_TAG,

        "-stage",
        0
    )


    # ============================================================
    # 10. INITIAL VERTICAL CONFINEMENT
    # ============================================================

    ops.timeSeries(
        "Linear",
        1
    )

    ops.pattern(
        "Plain",
        1,
        1
    )


    # ------------------------------------------------------------
    # Since top nodes are tied,
    # apply resultant vertical force to master node.
    # ------------------------------------------------------------

    ops.load(

        5,

        0.0,

        0.0,

        -total_vertical_load
    )


    # ============================================================
    # 11. INITIAL COMPRESSION ANALYSIS
    # ============================================================

    ops.constraints(
        "Transformation"
    )

    ops.numberer(
        "RCM"
    )

    ops.system(
        "BandGeneral"
    )

    ops.test(

        "NormDispIncr",

        1.0e-10,

        50
    )

    ops.algorithm(
        "Newton"
    )

    ops.integrator(

        "LoadControl",

        0.01
    )

    ops.analysis(
        "Static"
    )


    print()
    print("========================================")
    print(" INITIAL CONFINEMENT")
    print("========================================")


    ok = ops.analyze(100)


    if ok != 0:

        raise RuntimeError(
            "Initial compression analysis failed."
        )


    print("Initial confinement completed.")


    # ============================================================
    # 12. FIX INITIAL LOAD
    # ============================================================

    ops.loadConst(
        "-time",
        0.0
    )


    # ============================================================
    # 13. SWITCH TO ELASTOPLASTIC MATERIAL
    # ============================================================

    ops.updateMaterialStage(

        "-material",
        MAT_TAG,

        "-stage",
        1
    )


    print(
        "Material stage -> elastoplastic"
    )


    # ============================================================
    # 14. DEFINE CYCLIC SHEAR HISTORY
    # ============================================================
    #
    # One cycle:
    #
    #       +gamma_max
    #          /\
    #         /  \
    # 0 -----    -----
    #              \  /
    #               \/
    #          -gamma_max
    #
    # ============================================================

    gamma_targets = [0.0]


    for cycle in range(n_cycles):

        points = [

            gamma_max,

            0.0,

            -gamma_max,

            0.0
        ]


        for target in points:

            start = gamma_targets[-1]

            segment = np.linspace(

                start,

                target,

                steps_per_quarter_cycle + 1

            )[1:]


            gamma_targets.extend(
                segment.tolist()
            )


    gamma_targets = np.array(
        gamma_targets
    )


    # ============================================================
    # 15. STORAGE
    # ============================================================

    result_gamma = []

    result_tau = []

    result_sigma_x = []

    result_sigma_y = []

    result_sigma_z = []

    result_eta = []
    result_strain = []
    result_stress = []
    result_reaction = []
    result_target = []
    analysis_error = None


    # ------------------------------------------------------------
    # Nodal displacement history
    #
    # Shape later:
    #
    # step × node × xyz
    # ------------------------------------------------------------

    node_displacement_history = []


    # ============================================================
    # 16. FUNCTION TO RECORD CURRENT NODE DISPLACEMENTS
    # ============================================================

    def get_all_node_displacements():

        displacement_step = []


        for node in range(1, 9):

            disp = ops.nodeDisp(node)


            displacement_step.append(

                [
                    disp[0],
                    disp[1],
                    disp[2]
                ]
            )


        return np.array(
            displacement_step
        )


    # ============================================================
    # 17. RECORD INITIAL STATE
    # ============================================================

    stress = ops.eleResponse(

        ELE_TAG,

        "stress"
    )


    strain = ops.eleResponse(

        ELE_TAG,

        "strain"
    )


    result_sigma_x.append(
        stress[0]
    )

    result_sigma_y.append(
        stress[1]
    )

    result_sigma_z.append(
        stress[2]
    )


    result_tau.append(
        stress[5]
    )


    result_gamma.append(
        strain[5]
    )


    if len(stress) >= 7:

        result_eta.append(
            stress[6]
        )


    result_stress.append(list(stress[:6]))
    result_strain.append(list(strain[:6]))
    ops.reactions()
    result_reaction.append([-sum(ops.nodeReaction(n, d) for n in [1, 2, 3, 4]) / AREA
                            for d in [1, 2, 3]])
    result_target.append(0.0)

    node_displacement_history.append(

        get_all_node_displacements()
    )


    # ============================================================
    # 18. CREATE REFERENCE LOAD FOR DISPLACEMENT CONTROL
    # ============================================================

    ops.timeSeries(
        "Linear",
        2
    )

    ops.pattern(
        "Plain",
        2,
        2
    )


    ops.load(

        5,

        1.0,

        0.0,

        0.0
    )


    # ============================================================
    # 19. CYCLIC SHEAR ANALYSIS
    # ============================================================

    print()
    print("========================================")
    print(" CYCLIC SIMPLE SHEAR")
    print("========================================")


    current_gamma = 0.0
    ux_at_shear_start = ops.nodeDisp(5, 1)


    for step in range(

        1,

        len(gamma_targets)

    ):


        # --------------------------------------------------------
        # Target strain
        # --------------------------------------------------------

        target_gamma = gamma_targets[step]


        delta_ux = ux_at_shear_start + target_gamma * H - ops.nodeDisp(5, 1)

        # --------------------------------------------------------
        # Reset analysis objects only
        # --------------------------------------------------------

        ops.wipeAnalysis()


        # --------------------------------------------------------
        # Analysis configuration
        # --------------------------------------------------------

        ops.constraints(
            "Transformation"
        )

        ops.numberer(
            "RCM"
        )

        ops.system(
            "BandGeneral"
        )

        ops.test(

            "NormDispIncr",

            1.0e-9,

            100
        )

        ops.algorithm(
            "Newton"
        )


        # --------------------------------------------------------
        # Displacement-controlled shear
        # --------------------------------------------------------

        ops.integrator(

            "DisplacementControl",

            5,

            1,

            delta_ux
        )


        ops.analysis(
            "Static"
        )


        # --------------------------------------------------------
        # Run one increment
        # --------------------------------------------------------

        ok = ops.analyze(1)


        if ok != 0:

            print(
                f"WARNING: analysis failed at step {step}"
            )

            analysis_error = f"cyclic step {step}: analyze returned {ok}"
            break


        # --------------------------------------------------------
        # Update target
        # --------------------------------------------------------

        current_gamma = target_gamma


        # ========================================================
        # MATERIAL RESPONSE
        # ========================================================

        stress = ops.eleResponse(

            ELE_TAG,

            "stress"
        )


        strain = ops.eleResponse(

            ELE_TAG,

            "strain"
        )


        # --------------------------------------------------------
        # PDMY02 3D stress
        #
        # [0] sigma_xx
        # [1] sigma_yy
        # [2] sigma_zz
        # [3] sigma_xy
        # [4] sigma_yz
        # [5] sigma_zx
        # [6] eta_r
        #
        #
        # strain
        #
        # [0] epsilon_xx
        # [1] epsilon_yy
        # [2] epsilon_zz
        # [3] gamma_xy
        # [4] gamma_yz
        # [5] gamma_zx
        # --------------------------------------------------------

        sigma_x = stress[0]

        sigma_y = stress[1]

        sigma_z = stress[2]

        tau_zx = stress[5]

        gamma_zx = strain[5]


        # --------------------------------------------------------
        # Save results
        # --------------------------------------------------------

        result_sigma_x.append(
            sigma_x
        )

        result_sigma_y.append(
            sigma_y
        )

        result_sigma_z.append(
            sigma_z
        )

        result_tau.append(
            tau_zx
        )

        result_gamma.append(
            gamma_zx
        )


        if len(stress) >= 7:

            result_eta.append(
                stress[6]
            )


        # --------------------------------------------------------
        # Save all node displacements for animation
        # --------------------------------------------------------

        result_stress.append(list(stress[:6]))
        result_strain.append(list(strain[:6]))
        ops.reactions()
        result_reaction.append([-sum(ops.nodeReaction(n, d) for n in [1, 2, 3, 4]) / AREA
                                for d in [1, 2, 3]])
        result_target.append(target_gamma)

        node_displacement_history.append(

            get_all_node_displacements()

        )


    print(
        "Cyclic shear analysis completed." if analysis_error is None
        else f"Cyclic shear analysis INCOMPLETE: {analysis_error}"
    )


    # ============================================================
    # 20. CONVERT RESULTS TO NUMPY
    # ============================================================

    result_gamma = np.array(
        result_gamma
    )

    result_tau = np.array(
        result_tau
    )

    result_sigma_x = np.array(
        result_sigma_x
    )

    result_sigma_y = np.array(
        result_sigma_y
    )

    result_sigma_z = np.array(
        result_sigma_z
    )

    result_eta = np.array(
        result_eta
    )


    node_displacement_history = np.array(

        node_displacement_history

    )


    # ============================================================
    # 21. RESULT SUMMARY
    # ============================================================

    print()
    print("========================================")
    print(" RESULT SUMMARY")
    print("========================================")


    print(

        "Maximum shear strain :",

        np.max(
            np.abs(result_gamma)
        )

    )


    print(

        "Maximum shear stress :",

        np.max(
            np.abs(result_tau)
        ),

        "kPa"

    )


    print(

        "Final sigma_z :",

        result_sigma_z[-1],

        "kPa"

    )



    return {
        "gamma": result_gamma, "tau": result_tau,
        "stress": np.asarray(result_stress), "strain": np.asarray(result_strain),
        "eta": result_eta, "reaction": np.asarray(result_reaction),
        "displacement": node_displacement_history, "target": np.asarray(result_target),
        "coordinates": original_coordinates,
        "completed": analysis_error is None and len(result_gamma) == len(gamma_targets),
        "error": analysis_error, "steps": steps_per_quarter_cycle,
        "cycles": n_cycles, "amplitude": gamma_max, "sigma_v": sigma_v,
        "height": H, "rho": rho, "opensees_version": ops.version(),
        "material": {"G_ref": G_ref, "B_ref": B_ref, "ref_pressure": ref_pressure,
                     "friction_angle": friction_angle, "PT_angle": PT_angle,
                     "peak_shear_strain": peak_shear_strain,
                     "pressure_depend_coeff": pressure_depend_coeff,
                     "contrac1": contrac1, "contrac3": contrac3,
                     "dilat1": dilat1, "dilat3": dilat3,
                     "n_yield_surface": n_yield_surface},
    }


def plot_results(result):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    result_gamma = result["gamma"]
    result_tau = result["tau"]
    result_sigma_x, result_sigma_y, result_sigma_z = result["stress"][:, :3].T
    node_displacement_history = result["displacement"]
    original_coordinates = result["coordinates"]
    gamma_max = result["amplitude"]
    # ============================================================
    # 22. STATIC HYSTERESIS GRAPH
    # ============================================================

    plt.figure(
        figsize=(7, 6)
    )


    plt.plot(

        result_gamma * 100.0,

        result_tau,

        linewidth=1.5

    )


    plt.xlabel(

        r"Shear strain $\gamma_{zx}$ (%)"

    )


    plt.ylabel(

        r"Shear stress $\tau_{zx}$ (kPa)"

    )


    plt.title(

        "PressureDependMultiYield02\n"
        "Cyclic Simple Shear"

    )


    plt.grid(
        True
    )


    plt.tight_layout()


    plt.show()


    # ============================================================
    # 23. NORMAL STRESS GRAPH
    # ============================================================

    plt.figure(
        figsize=(7, 6)
    )


    plt.plot(

        result_gamma * 100.0,

        result_sigma_x,

        label=r"$\sigma_{xx}$"

    )


    plt.plot(

        result_gamma * 100.0,

        result_sigma_y,

        label=r"$\sigma_{yy}$"

    )


    plt.plot(

        result_gamma * 100.0,

        result_sigma_z,

        label=r"$\sigma_{zz}$"

    )


    plt.xlabel(

        r"Shear strain $\gamma_{zx}$ (%)"

    )


    plt.ylabel(
        "Normal stress (kPa)"
    )


    plt.title(
        "Normal Stress Evolution"
    )


    plt.legend()


    plt.grid(
        True
    )


    plt.tight_layout()


    plt.show()


    # ============================================================
    # 24. 3D ELEMENT ANIMATION
    # ============================================================

    if SHOW_ANIMATION:


        # --------------------------------------------------------
        # Node order for plotting
        # --------------------------------------------------------

        node_array = np.array(

            [
                original_coordinates[i]
                for i in range(1, 9)
            ]

        )


        # --------------------------------------------------------
        # Element edges
        # --------------------------------------------------------

        edges = [

            # Bottom
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),

            # Top
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),

            # Vertical
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7)
        ]


        # ========================================================
        # CREATE FIGURE
        # ========================================================

        fig = plt.figure(

            figsize=(13, 6)

        )


        # --------------------------------------------------------
        # Left: 3D element
        # --------------------------------------------------------

        ax3d = fig.add_subplot(

            1,
            2,
            1,

            projection="3d"

        )


        # --------------------------------------------------------
        # Right: hysteresis
        # --------------------------------------------------------

        ax_hyst = fig.add_subplot(

            1,
            2,
            2

        )


        # ========================================================
        # AXIS LIMITS
        # ========================================================

        max_visual_shear = (

            gamma_max
            *
            H
            *
            VIS_SCALE

        )


        margin = 0.2


        ax3d.set_xlim(

            -margin,

            Lx
            +
            max_visual_shear
            +
            margin

        )


        ax3d.set_ylim(

            -margin,

            Ly + margin

        )


        ax3d.set_zlim(

            -margin,

            H + margin

        )


        ax3d.set_xlabel(
            "X"
        )

        ax3d.set_ylabel(
            "Y"
        )

        ax3d.set_zlabel(
            "Z"
        )


        ax3d.set_title(

            f"SSPbrick deformation\n"
            f"visual scale = {VIS_SCALE}x"

        )


        # --------------------------------------------------------
        # Try equal visual proportions
        # --------------------------------------------------------

        ax3d.set_box_aspect(

            [
                1 + max_visual_shear,
                1,
                1
            ]

        )


        # ========================================================
        # ORIGINAL ELEMENT
        # ========================================================

        for i, j in edges:

            ax3d.plot(

                [
                    node_array[i, 0],
                    node_array[j, 0]
                ],

                [
                    node_array[i, 1],
                    node_array[j, 1]
                ],

                [
                    node_array[i, 2],
                    node_array[j, 2]
                ],

                linestyle="--",

                linewidth=0.8,

                alpha=0.4

            )


        # ========================================================
        # DEFORMED ELEMENT LINES
        # ========================================================

        deformed_lines = []


        initial_displacement = (

            node_displacement_history[0]

        )


        initial_deformed = (

            node_array
            +
            VIS_SCALE
            *
            initial_displacement

        )


        for i, j in edges:

            line, = ax3d.plot(

                [
                    initial_deformed[i, 0],
                    initial_deformed[j, 0]
                ],

                [
                    initial_deformed[i, 1],
                    initial_deformed[j, 1]
                ],

                [
                    initial_deformed[i, 2],
                    initial_deformed[j, 2]
                ],

                linewidth=2.0

            )


            deformed_lines.append(
                line
            )


        # --------------------------------------------------------
        # Nodes
        # --------------------------------------------------------

        node_scatter = ax3d.scatter(

            initial_deformed[:, 0],

            initial_deformed[:, 1],

            initial_deformed[:, 2],

            s=35

        )


        # ========================================================
        # HYSTERESIS AXIS
        # ========================================================

        gamma_percent = (

            result_gamma
            *
            100.0

        )


        gamma_limit = (

            max(
                np.max(
                    np.abs(gamma_percent)
                ),
                0.1
            )

            *
            1.15

        )


        tau_limit = (

            max(
                np.max(
                    np.abs(result_tau)
                ),
                1.0
            )

            *
            1.15

        )


        ax_hyst.set_xlim(

            -gamma_limit,

            gamma_limit

        )


        ax_hyst.set_ylim(

            -tau_limit,

            tau_limit

        )


        ax_hyst.set_xlabel(

            r"Shear strain $\gamma_{zx}$ (%)"

        )


        ax_hyst.set_ylabel(

            r"Shear stress $\tau_{zx}$ (kPa)"

        )


        ax_hyst.set_title(

            r"$\tau_{zx}$ - $\gamma_{zx}$ hysteresis"

        )


        ax_hyst.grid(
            True
        )


        hysteresis_line, = ax_hyst.plot(

            [],

            [],

            linewidth=2

        )


        current_point, = ax_hyst.plot(

            [],

            [],

            marker="o",

            markersize=7

        )


        # ========================================================
        # STATUS TEXT
        # ========================================================

        status_text = fig.text(

            0.5,

            0.02,

            "",

            ha="center"

        )


        # ========================================================
        # ANIMATION UPDATE FUNCTION
        # ========================================================

        def update_animation(frame):


            # ----------------------------------------------------
            # Current displacements
            # ----------------------------------------------------

            displacement = (

                node_displacement_history[
                    frame
                ]

            )


            # ----------------------------------------------------
            # Exaggerated deformation for visualization only
            # ----------------------------------------------------

            deformed = (

                node_array

                +

                VIS_SCALE
                *
                displacement

            )


            # ----------------------------------------------------
            # Update element edges
            # ----------------------------------------------------

            for line, edge in zip(

                deformed_lines,

                edges

            ):


                i, j = edge


                line.set_data(

                    [
                        deformed[i, 0],
                        deformed[j, 0]
                    ],

                    [
                        deformed[i, 1],
                        deformed[j, 1]
                    ]

                )


                line.set_3d_properties(

                    [
                        deformed[i, 2],
                        deformed[j, 2]
                    ]

                )


            # ----------------------------------------------------
            # Update nodes
            # ----------------------------------------------------

            node_scatter._offsets3d = (

                deformed[:, 0],

                deformed[:, 1],

                deformed[:, 2]

            )


            # ----------------------------------------------------
            # Update hysteresis history
            # ----------------------------------------------------

            hysteresis_line.set_data(

                gamma_percent[
                    :frame + 1
                ],

                result_tau[
                    :frame + 1
                ]

            )


            # ----------------------------------------------------
            # Current state marker
            # ----------------------------------------------------

            current_point.set_data(

                [
                    gamma_percent[
                        frame
                    ]
                ],

                [
                    result_tau[
                        frame
                    ]
                ]

            )


            # ----------------------------------------------------
            # Update status text
            # ----------------------------------------------------

            status_text.set_text(

                f"Step {frame:4d}    |    "

                f"γzx = "
                f"{gamma_percent[frame]:+.3f} %    |    "

                f"τzx = "
                f"{result_tau[frame]:+.2f} kPa    |    "

                f"σzz = "
                f"{result_sigma_z[frame]:+.2f} kPa"

            )


            return (

                deformed_lines

                +
                [
                    node_scatter,
                    hysteresis_line,
                    current_point,
                    status_text
                ]

            )


        # ========================================================
        # CREATE ANIMATION
        # ========================================================

        animation = FuncAnimation(

            fig,

            update_animation,

            frames=len(
                node_displacement_history
            ),

            interval=ANIMATION_INTERVAL,

            blit=False,

            repeat=True

        )


        plt.tight_layout(

            rect=[
                0,
                0.05,
                1,
                1
            ]

        )


        # ========================================================
        # OPTIONAL GIF EXPORT
        # ========================================================

        if SAVE_ANIMATION:

            from matplotlib.animation import PillowWriter


            print()

            print(
                "Saving animation..."
            )


            animation.save(

                ANIMATION_FILENAME,

                writer=PillowWriter(
                    fps=25
                )

            )


            print(

                f"Saved: "
                f"{ANIMATION_FILENAME}"

            )


        plt.show()



def main():
    parser = argparse.ArgumentParser(description="PDMY02 cyclic test + six-item evaluation")
    parser.add_argument("--no-plots", action="store_true", help="evaluate without GUI")
    parser.add_argument("--quick", action="store_true", help="skip independent probe/refinement runs")
    parser.add_argument("--gamma", type=float, default=gamma_max, help="strain fraction, not percent")
    parser.add_argument("--steps", type=int, default=steps_per_quarter_cycle)
    parser.add_argument("--cycles", type=int, default=n_cycles)
    parser.add_argument("--targets", type=Path, default=Path(__file__).with_name("material_targets.json"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("material_results"))
    args = parser.parse_args()
    # Reject invalid targets before spending time on analysis.
    from material_evaluation import load_targets
    targets = load_targets(args.targets)
    runs = {}
    try:
        base = run_test(args.gamma, args.steps, args.cycles)
        if not args.quick and base["completed"]:
            import contextlib
            import io
            cases = [("probe", 1e-7, 20, 2), ("probe_half", 5e-8, 20, 2),
                     ("fine", args.gamma, args.steps * 2, args.cycles),
                     ("finer", args.gamma, args.steps * 4, args.cycles)]
            for name, amplitude, steps, cycles in cases:
                print(f"Evaluation run: {name}, gamma={amplitude:g}, quarter_steps={steps}", flush=True)
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        runs[name] = run_test(amplitude, steps, cycles)
                except RuntimeError as exc:
                    runs[name] = {"completed": False, "error": str(exc)}
        report = evaluate(base, runs, targets)
        save_report(report, base, runs, args.output_dir)
        if not args.no_plots:
            plot_results(base)
        return 2 if report["consistency_status"] == "FAIL" or report["target_status"] == "FAIL" else 0
    finally:
        ops.wipe()


if __name__ == "__main__":
    raise SystemExit(main())
