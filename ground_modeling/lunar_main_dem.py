import ground_modeling as gm
import lunar_main as lm

from lunar_config import ANA, GEO
from lunar_dem import TerrainSampler
from lunar_ground import LunarGround
from lunar_structure import LunarStructure


def main():
    # ========================================================
    # 1. 달 DEM 읽기
    # ========================================================

    print("[1] Reading lunar DEM...")

    (
        X,
        Y,
        Z_rel,
        Z_abs,
        rows,
        cols,
    ) = gm.extract_lunar_area_patch(
        gm.XML_PATH,
        gm.IMG_PATH,
    )

    print("DEM shape:", Z_rel.shape)
    print(
        "DEM X range:",
        float(X.min()),
        "~",
        float(X.max()),
    )
    print(
        "DEM Y range:",
        float(Y.min()),
        "~",
        float(Y.max()),
    )
    print(
        "DEM Z range:",
        float(Z_rel.min()),
        "~",
        float(Z_rel.max()),
    )

    # ========================================================
    # 2. DEM 보간기 생성
    # ========================================================

    terrain = TerrainSampler(
        X,
        Y,
        Z_rel,
    )

    terrain_center_x = 0.5 * (
        float(X.min()) + float(X.max())
    )
    terrain_center_y = 0.5 * (
        float(Y.min()) + float(Y.max())
    )

    # 구조물 길이 중심을 DEM 중앙에 배치
    structure_origin_x = (
        terrain_center_x - GEO.length / 2.0
    )
    structure_origin_y = terrain_center_y

    print(
        "Structure origin:",
        structure_origin_x,
        structure_origin_y,
    )

    # ========================================================
    # 3. OpenSees domain 초기화
    # ========================================================

    print("[2] Initializing 3D OpenSees domain...")

    tags = lm.initialize_domain()

    # ========================================================
    # 4. 구조물 생성
    # ========================================================

    print("[3] Building lunar structure...")

    structure = LunarStructure(
        tags=tags,
        terrain_z=terrain,
        x_origin=structure_origin_x,
        y_origin=structure_origin_y,
    )

    structure.build()

    print(
        "Structural nodes:",
        len(structure.coords),
    )
    print(
        "Frame elements:",
        len(structure.frame_elements),
    )
    print(
        "Foundation base nodes:",
        structure.base_nodes,
    )

    # ========================================================
    # 5. 기초 접촉 모델 생성
    # ========================================================

    print("[4] Building ground-contact model...")

    ground = LunarGround(
        tags,
        structure.base_nodes,
        structure.coords,
    )

    ground.build(
        nonlinear_contact=ANA.use_nonlinear_contact
    )

    # ========================================================
    # 6. 하중 적용 및 정적해석
    # ========================================================

    print("[5] Applying loads...")

    structure.apply_static_loads()

    print("[6] Running static analysis...")

    lm.run_gravity_and_pressure_analysis()

    # ========================================================
    # 7. 고유치 해석
    # ========================================================

    print("[7] Running eigenvalue analysis...")

    try:
        frequencies = lm.run_eigen_analysis(
            ANA.n_modes
        )

    except Exception as error:
        print(
            "[WARN] Eigenvalue analysis failed:",
            error,
        )
        frequencies = []

    # ========================================================
    # 8. 결과 저장
    # ========================================================

    print("[8] Saving results...")

    output_dir = lm.write_csv_results(
        structure,
        ground,
        frequencies,
    )

    # 현재 lunar_main.py의 기존 시각화 함수 호출
    html_file = lm.write_html_visualization(
        structure,
        ground,
        output_dir,
    )

    lm.print_summary(
        structure,
        ground,
        frequencies,
    )

    print("\nCompleted")
    print("Results:", output_dir.resolve())
    print("HTML:", html_file.resolve())


if __name__ == "__main__":
    main()

