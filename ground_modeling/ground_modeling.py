import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import plotly.graph_objects as go
import webbrowser

import openseespy.opensees as ops


# ============================================================
# 0. 실행 선택
# ============================================================

# "Earth", "Mars", "Moon" 중 하나
TARGET_BODY = "Moon"

# 한 번 실행할 때 하나만 생성:
# "surface" 또는 "solid"
VISUALIZATION_MODE = "solid"

# 선택한 천체의 DEM과 지반 재료를 이용해 OpenSees 2D 지반모델을 생성한다.
BUILD_OPENSEES_MODEL = True

# HTML을 생성한 뒤 브라우저에서 자동으로 열기
AUTO_OPEN_HTML = True


# ============================================================
# 1. 지구 설정
# ============================================================

# Earth_file : NASADEM HGT 파일 경로
HGT_PATH = Path("n36e129.hgt")

# HGT 타일 정보
# n36e129.hgt는 위도 36~37, 경도 129~130 영역
TILE_LAT_MIN = 36.0
TILE_LON_MIN = 129.0

# 포항 흥해읍 근처 단면 좌표
# 동서 방향 단면
AREA_CENTER_LAT = 36.11
AREA_CENTER_LON = 129.36

# 시각화 영역 크기
AREA_WIDTH_M = 3000.0
AREA_HEIGHT_M = 3000.0
N_X = 70
N_Y = 70

# OpenSees 지반 모델 깊이
BASE_DEPTH = 30.0  # m

# 깊이 방향 요소 개수₩₩
N_DEPTH = 12

# 지반 물성값: 1차 단순 모델
E_SOIL = 50.0e6       # Pa = N/m^2
NU_SOIL = 0.30
RHO_SOIL = 1800.0     # kg/m^3


# ============================================================
# 1. 달 설정
# ============================================================

# Moon_file : IMG, XML 파일 경로
IMG_PATH = Path("ldem_1024_90s_75s_000_030.img")
XML_PATH = Path("ldem_1024_90s_75s_000_030.xml")


# 달 PDS4 DEM은 매우 크므로 전체 원본을 그대로 Plotly에 넣지 않는다.
# 먼저 원본 raster에서 일부 영역을 읽고, 지구와 동일한 N_X × N_Y 격자로 변환한다.
MOON_READ_MODE = "center_crop"       # "center_crop" 또는 "full_downsample"
MOON_CROP_WIDTH_PIXELS = 3000
MOON_CROP_HEIGHT_PIXELS = 3000

# None이면 파일 중앙을 사용한다.
MOON_CROP_COL_OFFSET = None
MOON_CROP_ROW_OFFSET = None

# ============================================================
# 1-1. 화성 설정
# ============================================================

# 화성 DEM 파일 경로는 사용자가 자료를 받은 뒤 직접 입력한다.
# 예시: MARS_DEM_PATH = Path("mars_dem.tif")
MARS_LABEL_PATH = Path("megt00n000hb.lbl")
MARS_IMAGE_PATH = Path("megt00n000hb.img")

# 화성 DEM이 매우 큰 경우 일부 영역만 읽거나 전체를 다운샘플링한다.
MARS_READ_MODE = "center_crop"       # "center_crop" 또는 "full_downsample"
MARS_CROP_WIDTH_PIXELS = 3000
MARS_CROP_HEIGHT_PIXELS = 3000
MARS_CROP_COL_OFFSET = None
MARS_CROP_ROW_OFFSET = None


# ============================================================
# 1-1. 지반 재료 물성 및 재료 모델 설정
# ============================================================
# 단위계:
#   힘: N
#   길이: m
#   질량: kg
#   응력/탄성계수: Pa = N/m^2
#
# 모델링 방침:
#   - 점착성/암반계 재료: ElasticIsotropic
#   - 비점착성 재료 sand, gravel: PressureDependMultiYield02
#   - PDMY02 재료는 gravity 단계에서 stage 0, 이후 stage 1로 전환

USE_PDMY02_FOR_GRANULAR_SOIL = True

SOIL_MATERIALS = {
    "organic": {
        "tag": 1,
        "model": "ElasticIsotropic",
        "label": "Organic / topsoil",
        "E": 5.0e6,
        "nu": 0.35,
        "rho": 1300.0,
        "color": "#5b3a29",
    },
    "clay": {
        "tag": 2,
        "model": "ElasticIsotropic",
        "label": "Clay",
        "E": 20.0e6,
        "nu": 0.45,
        "rho": 1700.0,
        "color": "#8b6f47",
    },
    "silt": {
        "tag": 3,
        "model": "ElasticIsotropic",
        "label": "Silt",
        "E": 35.0e6,
        "nu": 0.38,
        "rho": 1750.0,
        "color": "#c2b280",
    },
    "sand": {
        "tag": 4,
        "model": "PressureDependMultiYield02",
        "label": "Sand",
        "rho": 1850.0,

        # PDMY02 핵심 파라미터
        # 아래 값은 초기 예시값이다. 실제 값은 지반조사 자료로 보정해야 한다.
        "nd": 2,
        "refShearModul": 90.0e6,
        "refBulkModul": 220.0e6,
        "frictionAng": 32.0,
        "peakShearStra": 0.10,
        "refPress": 101.0e3,
        "pressDependCoe": 0.50,
        "PTAng": 26.0,
        "contrac1": 0.067,
        "contrac3": 0.23,
        "dilat1": 0.06,
        "dilat3": 0.27,
        "noYieldSurf": 20,
        "contrac2": 5.0,
        "dilat2": 3.0,
        "liquefac1": 1.0,
        "liquefac2": 0.0,
        "e": 0.77,
        "cs1": 0.9,
        "cs2": 0.02,
        "cs3": 0.7,
        "pa": 101.0e3,
        "c": 0.1e3,

        "color": "#f2d16b",
    },
    "gravel": {
        "tag": 5,
        "model": "PressureDependMultiYield02",
        "label": "Gravel",
        "rho": 2000.0,

        # 조밀한 비점착성 재료 예시값
        "nd": 2,
        "refShearModul": 130.0e6,
        "refBulkModul": 260.0e6,
        "frictionAng": 36.5,
        "peakShearStra": 0.10,
        "refPress": 101.0e3,
        "pressDependCoe": 0.50,
        "PTAng": 26.0,
        "contrac1": 0.013,
        "contrac3": 0.0,
        "dilat1": 0.30,
        "dilat3": 0.0,
        "noYieldSurf": 20,
        "contrac2": 5.0,
        "dilat2": 3.0,
        "liquefac1": 1.0,
        "liquefac2": 0.0,
        "e": 0.55,
        "cs1": 0.9,
        "cs2": 0.02,
        "cs3": 0.7,
        "pa": 101.0e3,
        "c": 0.1e3,

        "color": "#9e9e9e",
    },
    "weathered_rock": {
        "tag": 6,
        "model": "ElasticIsotropic",
        "label": "Weathered rock",
        "E": 500.0e6,
        "nu": 0.25,
        "rho": 2200.0,
        "color": "#6f7f80",
    },
    "rock": {
        "tag": 7,
        "model": "ElasticIsotropic",
        "label": "Rock",
        "E": 2.0e9,
        "nu": 0.22,
        "rho": 2500.0,
        "color": "#3f4a4d",
    },
}

SOIL_ORDER = [
    "organic",
    "clay",
    "silt",
    "sand",
    "gravel",
    "weathered_rock",
    "rock",
]


# ============================================================
# 1-2. 달 지반 재료 물성 및 재료 모델 설정
# ============================================================
# 주의:
#   아래 값은 달 레골리스 해석 구조를 구현하기 위한 초기 가정값이다.
#   최종 연구에서는 대상 착륙지 및 참고 문헌에 따라 보정해야 한다.
#
# 모델링 방침:
#   - 느슨한/조밀한 레골리스: PressureDependMultiYield02
#   - 파쇄 암반/기반암: ElasticIsotropic

MOON_SOIL_MATERIALS = {
    "loose_regolith": {
        "tag": 101,
        "model": "PressureDependMultiYield02",
        "label": "Loose lunar regolith",
        "rho": 1500.0,
        "nd": 2,
        "refShearModul": 35.0e6,
        "refBulkModul": 80.0e6,
        "frictionAng": 30.0,
        "peakShearStra": 0.10,
        "refPress": 101.0e3,
        "pressDependCoe": 0.50,
        "PTAng": 24.0,
        "contrac1": 0.08,
        "contrac3": 0.20,
        "dilat1": 0.03,
        "dilat3": 0.15,
        "noYieldSurf": 20,
        "contrac2": 5.0,
        "dilat2": 3.0,
        "liquefac1": 0.0,
        "liquefac2": 0.0,
        "e": 0.90,
        "cs1": 0.9,
        "cs2": 0.02,
        "cs3": 0.7,
        "pa": 101.0e3,
        "c": 0.05e3,
        "color": "#c9b18f",
    },
    "dense_regolith": {
        "tag": 102,
        "model": "PressureDependMultiYield02",
        "label": "Dense lunar regolith",
        "rho": 1800.0,
        "nd": 2,
        "refShearModul": 80.0e6,
        "refBulkModul": 170.0e6,
        "frictionAng": 38.0,
        "peakShearStra": 0.08,
        "refPress": 101.0e3,
        "pressDependCoe": 0.50,
        "PTAng": 28.0,
        "contrac1": 0.03,
        "contrac3": 0.10,
        "dilat1": 0.12,
        "dilat3": 0.10,
        "noYieldSurf": 20,
        "contrac2": 5.0,
        "dilat2": 3.0,
        "liquefac1": 0.0,
        "liquefac2": 0.0,
        "e": 0.65,
        "cs1": 0.9,
        "cs2": 0.02,
        "cs3": 0.7,
        "pa": 101.0e3,
        "c": 0.10e3,
        "color": "#9b8265",
    },
    "fractured_rock": {
        "tag": 103,
        "model": "ElasticIsotropic",
        "label": "Fractured lunar rock",
        "E": 300.0e6,
        "nu": 0.25,
        "rho": 2200.0,
        "color": "#6f6f6f",
    },
    "bedrock": {
        "tag": 104,
        "model": "ElasticIsotropic",
        "label": "Lunar bedrock",
        "E": 2.0e9,
        "nu": 0.22,
        "rho": 2700.0,
        "color": "#3f3f3f",
    },
}

MOON_SOIL_ORDER = [
    "loose_regolith",
    "dense_regolith",
    "fractured_rock",
    "bedrock",
]




# ============================================================
# 1-3. 화성 지반 재료 물성 및 재료 모델 설정
# ============================================================
# 보고서의 화성 레골리스 범위에 맞춘 초기 모델이다.
#   밀도: 1500~1800 kg/m^3
#   탄성계수: 약 10~100 MPa
#   내부마찰각: 30~45 deg
#   포아송비: 0.25~0.35
#   점착력: 0~5 kPa
#
# 아래 PDMY02 파라미터는 해석 구조 구현용 초기값이다.
# 실제 DEM 위치와 탐사 자료가 확정되면 반드시 보정해야 한다.

MARS_SOIL_MATERIALS = {
    "loose_regolith": {
        "tag": 201,
        "model": "PressureDependMultiYield02",
        "label": "Loose Martian regolith",
        "rho": 1500.0,
        "nd": 2,
        "refShearModul": 12.0e6,
        "refBulkModul": 30.0e6,
        "frictionAng": 32.0,
        "peakShearStra": 0.10,
        "refPress": 101.0e3,
        "pressDependCoe": 0.50,
        "PTAng": 24.0,
        "contrac1": 0.08,
        "contrac3": 0.20,
        "dilat1": 0.03,
        "dilat3": 0.15,
        "noYieldSurf": 20,
        "contrac2": 5.0,
        "dilat2": 3.0,
        "liquefac1": 0.0,
        "liquefac2": 0.0,
        "e": 0.90,
        "cs1": 0.9,
        "cs2": 0.02,
        "cs3": 0.7,
        "pa": 101.0e3,
        "c": 0.5e3,
        "color": "#c96f45",
    },
    "dense_regolith": {
        "tag": 202,
        "model": "PressureDependMultiYield02",
        "label": "Dense Martian regolith",
        "rho": 1800.0,
        "nd": 2,
        "refShearModul": 40.0e6,
        "refBulkModul": 90.0e6,
        "frictionAng": 40.0,
        "peakShearStra": 0.08,
        "refPress": 101.0e3,
        "pressDependCoe": 0.50,
        "PTAng": 28.0,
        "contrac1": 0.03,
        "contrac3": 0.10,
        "dilat1": 0.12,
        "dilat3": 0.10,
        "noYieldSurf": 20,
        "contrac2": 5.0,
        "dilat2": 3.0,
        "liquefac1": 0.0,
        "liquefac2": 0.0,
        "e": 0.65,
        "cs1": 0.9,
        "cs2": 0.02,
        "cs3": 0.7,
        "pa": 101.0e3,
        "c": 2.0e3,
        "color": "#a94f32",
    },
    "cemented_regolith": {
        "tag": 203,
        "model": "ElasticIsotropic",
        "label": "Cemented Martian regolith",
        "E": 100.0e6,
        "nu": 0.30,
        "rho": 1850.0,
        "color": "#7f3f2c",
    },
    "bedrock": {
        "tag": 204,
        "model": "ElasticIsotropic",
        "label": "Martian bedrock",
        "E": 2.0e9,
        "nu": 0.25,
        "rho": 2600.0,
        "color": "#4f342d",
    },
}

MARS_SOIL_ORDER = [
    "loose_regolith",
    "dense_regolith",
    "cemented_regolith",
    "bedrock",
]


def assign_soil_material(depth_from_surface, x_norm=0.0, y_norm=0.0, elevation_norm=0.0):
    """
    깊이와 위치에 따라 지반 재료를 결정한다.

    depth_from_surface:
        지표면에서 아래 방향 깊이 [m]
    x_norm, y_norm:
        0~1 범위의 정규화 좌표
    elevation_norm:
        0~1 범위의 상대 고도
    """

    d = float(depth_from_surface)

    # 랜덤이 아니라 deterministic heterogeneity.
    # 매번 실행해도 같은 재료 분포가 나온다.
    lens = (
        0.9 * np.sin(2.0 * np.pi * x_norm)
        + 0.6 * np.cos(2.0 * np.pi * y_norm)
        + 0.4 * np.sin(4.0 * np.pi * (x_norm + y_norm))
    )

    # 0 ~ 1.5 m: 유기물/표토층
    if d < 1.5:
        return "organic"

    # 1.5 ~ 5 m: 점토-실트 혼재층
    if d < 5.0:
        if lens + 0.3 * elevation_norm > 0.35:
            return "silt"
        return "clay"

    # 5 ~ 10 m: 실트-모래 전이층
    if d < 10.0:
        if lens > 0.45:
            return "sand"
        elif lens < -0.55:
            return "clay"
        return "silt"

    # 10 ~ 16 m: 모래 우세, 일부 자갈
    if d < 16.0:
        if lens > 0.25:
            return "gravel"
        return "sand"

    # 16 ~ 24 m: 풍화암 + 자갈 전이층
    if d < 24.0:
        if lens > 0.55:
            return "gravel"
        return "weathered_rock"

    # 24 m 이하: 기반암
    return "rock"



def assign_moon_material(
    depth_from_surface,
    x_norm=0.0,
    y_norm=0.0,
    elevation_norm=0.0,
):
    """
    깊이에 따라 달 지반 재료를 결정한다.

    현재는 단순한 층상 모델을 사용한다.
    최종 연구에서는 착륙 후보지의 레골리스 두께와 물성 자료로 보정해야 한다.
    """
    d = float(depth_from_surface)

    if d < 2.0:
        return "loose_regolith"
    if d < 8.0:
        return "dense_regolith"
    if d < 15.0:
        return "fractured_rock"
    return "bedrock"




def assign_mars_material(
    depth_from_surface,
    x_norm=0.0,
    y_norm=0.0,
    elevation_norm=0.0,
):
    """깊이에 따라 화성 레골리스 및 암반 재료를 배정한다."""
    d = float(depth_from_surface)

    if d < 3.0:
        return "loose_regolith"
    if d < 12.0:
        return "dense_regolith"
    if d < 22.0:
        return "cemented_regolith"
    return "bedrock"


# 중력가속도
G_EARTH = 9.81        # m/s^2
G_MARS = 3.71         # m/s^2
G_MOON = 1.62         # m/s^2


EARTH_CONFIG = {
    "body_name": "Earth",
    "gravity": G_EARTH,
    "materials": SOIL_MATERIALS,
    "material_order": SOIL_ORDER,
    "assign_material": assign_soil_material,
    "base_depth": BASE_DEPTH,
}

MOON_CONFIG = {
    "body_name": "Moon",
    "gravity": G_MOON,
    "materials": MOON_SOIL_MATERIALS,
    "material_order": MOON_SOIL_ORDER,
    "assign_material": assign_moon_material,
    "base_depth": BASE_DEPTH,
}


MARS_CONFIG = {
    "body_name": "Mars",
    "gravity": G_MARS,
    "materials": MARS_SOIL_MATERIALS,
    "material_order": MARS_SOIL_ORDER,
    "assign_material": assign_mars_material,
    "base_depth": BASE_DEPTH,
}


def register_opensees_soil_materials(materials, material_order):
    """
    OpenSees nDMaterial 등록.

    입력된 재료 딕셔너리와 순서에 따라 OpenSees nDMaterial을 등록한다.

    - ElasticIsotropic 재료는 탄성체로 등록
    - PressureDependMultiYield02 재료는 비선형 재료 태그 목록에 추가
    """

    nonlinear_material_tags = []

    for name in material_order:
        m = materials[name]
        mat_tag = int(m["tag"])
        model = m["model"]

        if model == "ElasticIsotropic":
            ops.nDMaterial(
                "ElasticIsotropic",
                mat_tag,
                float(m["E"]),
                float(m["nu"]),
                float(m["rho"]),
            )

        elif model == "PressureDependMultiYield02":
            ops.nDMaterial(
                "PressureDependMultiYield02",
                mat_tag,
                int(m["nd"]),
                float(m["rho"]),
                float(m["refShearModul"]),
                float(m["refBulkModul"]),
                float(m["frictionAng"]),
                float(m["peakShearStra"]),
                float(m["refPress"]),
                float(m["pressDependCoe"]),
                float(m["PTAng"]),
                float(m["contrac1"]),
                float(m["contrac3"]),
                float(m["dilat1"]),
                float(m["dilat3"]),
                int(m["noYieldSurf"]),
                float(m["contrac2"]),
                float(m["dilat2"]),
                float(m["liquefac1"]),
                float(m["liquefac2"]),
                float(m["e"]),
                float(m["cs1"]),
                float(m["cs2"]),
                float(m["cs3"]),
                float(m["pa"]),
                float(m["c"]),
            )

            nonlinear_material_tags.append(mat_tag)

        else:
            raise ValueError(f"지원하지 않는 재료 모델: {model}")

    return nonlinear_material_tags


def set_pdmy_material_stage(material_tags, stage):
    """
    PDMY/PDMY02 재료의 stage를 변경한다.

    stage 0:
        선형 탄성 상태.
        중력 재하 및 초기 응력장 형성용.

    stage 1:
        소성 상태.
        이후 말뚝 압입, 인장 재하, 지진 동적 재하 등에 사용.
    """

    for mat_tag in material_tags:
        ops.updateMaterialStage(
            "-material",
            int(mat_tag),
            "-stage",
            int(stage),
        )


def run_elastic_gravity_then_plastic_stage(
    nodal_self_weight_loads,
    nonlinear_material_tags,
    n_steps=10,
    d_lambda=0.1,
):
    """
    탄성 중력 재하 단계 → loadConst → 소성 전환 단계.

    1) PDMY/PDMY02 재료를 stage 0으로 설정
       - 전 소성 거동 비활성화
       - 순수 선형 탄성체로 자중 응력 상태 형성

    2) pattern Plain으로 지반 자중 하중 적용

    3) Static analysis 수행

    4) loadConst('-time', 0.0)으로 중력 재하 결과 고정

    5) PDMY/PDMY02 재료를 stage 1로 전환
       - 이후 외력 재하, 말뚝 압입, 지진 동적 재하 가능
    """

    print("\n[Stage 0] Set PDMY/PDMY02 materials to elastic stage.")
    set_pdmy_material_stage(nonlinear_material_tags, stage=0)

    # 중력 하중 패턴
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)

    for node_tag, load_y in nodal_self_weight_loads.items():
        ops.load(int(node_tag), 0.0, float(load_y))

    # 정적 탄성 중력 해석 설정
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-8, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", float(d_lambda))
    ops.analysis("Static")

    print("[Stage 0] Running elastic gravity analysis...")
    ok = ops.analyze(int(n_steps))

    if ok != 0:
        print("[Warning] Gravity analysis failed with Newton. Trying ModifiedNewton...")

        ops.wipeAnalysis()
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1.0e-7, 100)
        ops.algorithm("ModifiedNewton")
        ops.integrator("LoadControl", float(d_lambda))
        ops.analysis("Static")

        ok = ops.analyze(int(n_steps))

    if ok != 0:
        raise RuntimeError("탄성 중력 재하 해석이 수렴하지 않았음.")

    print("[Stage 0] Gravity analysis converged.")

    # 중력 응력장과 변위 상태 고정
    ops.loadConst("-time", 0.0)
    print("[Stage 0] Gravity state fixed by loadConst('-time', 0.0).")

    # 이후 해석을 위해 analysis object 초기화
    ops.wipeAnalysis()

    # 소성 상태로 전환
    print("[Stage 1] Switch PDMY/PDMY02 materials to plastic stage.")
    set_pdmy_material_stage(nonlinear_material_tags, stage=1)

    print("[Stage 1] Plastic material stage is ready.")
    print("         이후 말뚝 압입, 인장 재하, 지진 동적 가속도 재하를 추가하면 됨.\n")

# 모델 두께
THICKNESS = 1.0       # m


# ============================================================
# 2. NASADEM HGT 읽기 함수
# ============================================================

def read_hgt_file(hgt_path: Path):
    """
    NASADEM/SRTM 계열 .hgt 파일을 읽는다.
    .hgt는 big-endian signed 16-bit integer 형식이다.
    """
    raw = np.fromfile(hgt_path, dtype=">i2")
    n = int(np.sqrt(raw.size))

    if n * n != raw.size:
        raise ValueError(
            f"HGT 파일 크기가 정사각 격자와 맞지 않음: size={raw.size}"
        )

    dem = raw.reshape((n, n)).astype(float)

    # 결측값 처리
    dem[dem == -32768] = np.nan

    return dem


# ============================================================
# 3. 위도/경도 → DEM 인덱스 변환 및 보간
# ============================================================

def latlon_to_pixel(lat, lon, n, lat_min, lon_min):
    """
    HGT 파일은 일반적으로 북쪽에서 남쪽 방향으로 row가 증가한다.
    n36e129의 경우:
      lat: 36~37
      lon: 129~130

    row = 북쪽 37도에서 아래로 내려오는 방향
    col = 서쪽 129도에서 동쪽으로 가는 방향
    """
    lat_max = lat_min + 1.0
    lon_max = lon_min + 1.0

    row = (lat_max - lat) / (lat_max - lat_min) * (n - 1)
    col = (lon - lon_min) / (lon_max - lon_min) * (n - 1)

    return row, col


def bilinear_sample(dem, lat, lon, lat_min, lon_min):
    """
    DEM에서 특정 위도/경도의 고도값을 bilinear interpolation으로 추출한다.
    """
    n = dem.shape[0]
    row, col = latlon_to_pixel(lat, lon, n, lat_min, lon_min)

    r0 = int(np.floor(row))
    c0 = int(np.floor(col))
    r1 = r0 + 1
    c1 = c0 + 1

    if r0 < 0 or c0 < 0 or r1 >= n or c1 >= n:
        return np.nan

    dr = row - r0
    dc = col - c0

    z00 = dem[r0, c0]
    z01 = dem[r0, c1]
    z10 = dem[r1, c0]
    z11 = dem[r1, c1]

    values = np.array([z00, z01, z10, z11], dtype=float)

    if np.any(np.isnan(values)):
        return np.nan

    z = (
        z00 * (1 - dr) * (1 - dc)
        + z01 * (1 - dr) * dc
        + z10 * dr * (1 - dc)
        + z11 * dr * dc
    )

    return z


# ============================================================
# 4. 지구 단면 추출
# ============================================================

def extract_profile(dem):
    """
    포항 흥해읍 근처의 동서 방향 고도 단면을 추출한다.
    """
    lons = np.linspace(PROFILE_LON_START, PROFILE_LON_END, N_PROFILE)
    lats = np.full_like(lons, PROFILE_LAT)

    elevations = np.array([
        bilinear_sample(
            dem,
            lat=lat,
            lon=lon,
            lat_min=TILE_LAT_MIN,
            lon_min=TILE_LON_MIN,
        )
        for lat, lon in zip(lats, lons)
    ])

    # 거리 계산
    # 위도/경도 차이를 근사적으로 m 단위 거리로 변환
    earth_radius = 6371000.0  # m
    lat_rad = np.deg2rad(PROFILE_LAT)

    x = earth_radius * np.cos(lat_rad) * np.deg2rad(lons - lons[0])

    # 결측 제거
    valid = ~np.isnan(elevations)
    x = x[valid]
    lons = lons[valid]
    lats = lats[valid]
    elevations = elevations[valid]

    # 너무 큰 절대 고도보다 상대 고도가 모델링에 편함
    elevations_rel = elevations - np.nanmin(elevations)

    return x, elevations_rel, elevations, lats, lons

# ============================================================
# 4-1. 면적 영역 DEM 추출
# ============================================================

def extract_area_patch(dem):
    """
    중심 위도/경도를 기준으로 일정 크기의 직사각형 DEM 영역을 추출한다.
    결과는 Plotly 3D surface/solid 시각화에 사용한다.
    """

    earth_radius = 6371000.0
    lat_rad = np.deg2rad(AREA_CENTER_LAT)

    # m → degree 변환
    dlat = (AREA_HEIGHT_M / 2.0) / earth_radius * (180.0 / np.pi)
    dlon = (AREA_WIDTH_M / 2.0) / (earth_radius * np.cos(lat_rad)) * (180.0 / np.pi)

    lat_min = AREA_CENTER_LAT - dlat
    lat_max = AREA_CENTER_LAT + dlat
    lon_min = AREA_CENTER_LON - dlon
    lon_max = AREA_CENTER_LON + dlon

    lats = np.linspace(lat_min, lat_max, N_Y)
    lons = np.linspace(lon_min, lon_max, N_X)

    # 실제 거리 좌표계 생성
    x = earth_radius * np.cos(lat_rad) * np.deg2rad(lons - lons[0])
    y = earth_radius * np.deg2rad(lats - lats[0])

    X, Y = np.meshgrid(x, y)

    Z_abs = np.zeros_like(X, dtype=float)

    for iy, lat in enumerate(lats):
        for ix, lon in enumerate(lons):
            Z_abs[iy, ix] = bilinear_sample(
                dem,
                lat=lat,
                lon=lon,
                lat_min=TILE_LAT_MIN,
                lon_min=TILE_LON_MIN,
            )

    # 결측값이 있으면 주변 평균으로 간단 보정
    if np.any(np.isnan(Z_abs)):
        mean_val = np.nanmean(Z_abs)
        Z_abs[np.isnan(Z_abs)] = mean_val

    # 상대 고도
    Z_rel = Z_abs - np.nanmin(Z_abs)

    return X, Y, Z_rel, Z_abs, lats, lons



# ============================================================
# 4-2. 달 영역 추출
# ============================================================

def extract_lunar_area_patch(xml_path, img_path):
    """
    PDS4 XML 라벨을 통해 같은 폴더의 IMG DEM을 읽고,
    지구 extract_area_patch()와 동일하게 다음을 반환한다.

        X, Y, Z_rel, Z_abs, row_indices, col_indices

    핵심:
    - XML이 raster 크기, 자료형, scale/offset, 좌표계를 제공한다.
    - IMG는 실제 DEM 값을 저장한다.
    - 최종 출력 shape는 지구와 동일한 N_Y × N_X이다.
    """

    try:
        import rasterio
        from rasterio.enums import Resampling
        from rasterio.windows import Window
        from rasterio.transform import xy as transform_xy
    except ImportError as exc:
        raise ImportError(
            "달 PDS4 DEM을 읽으려면 rasterio/GDAL이 필요합니다.\n"
            "설치 명령:\n"
            "  conda install -c conda-forge rasterio gdal"
        ) from exc

    xml_path = Path(xml_path)
    img_path = Path(img_path)

    if not xml_path.exists():
        raise FileNotFoundError(f"달 XML 파일이 없습니다: {xml_path}")
    if not img_path.exists():
        raise FileNotFoundError(f"달 IMG 파일이 없습니다: {img_path}")

    try:
        src = rasterio.open(xml_path)
    except Exception as exc:
        raise RuntimeError(
            "PDS4 XML을 열지 못했습니다.\n"
            "XML과 IMG가 같은 폴더에 있는지, GDAL PDS4 드라이버가 "
            "설치되어 있는지 확인하세요.\n"
            f"원래 오류: {exc}"
        ) from exc

    with src:
        if src.count < 1:
            raise RuntimeError("PDS4 XML에서 raster band를 찾지 못했습니다.")

        if MOON_READ_MODE == "center_crop":
            crop_width = min(MOON_CROP_WIDTH_PIXELS, src.width)
            crop_height = min(MOON_CROP_HEIGHT_PIXELS, src.height)

            if MOON_CROP_COL_OFFSET is None:
                col_off = max((src.width - crop_width) // 2, 0)
            else:
                col_off = max(
                    min(MOON_CROP_COL_OFFSET, src.width - crop_width),
                    0,
                )

            if MOON_CROP_ROW_OFFSET is None:
                row_off = max((src.height - crop_height) // 2, 0)
            else:
                row_off = max(
                    min(MOON_CROP_ROW_OFFSET, src.height - crop_height),
                    0,
                )

            window = Window(
                col_off=col_off,
                row_off=row_off,
                width=crop_width,
                height=crop_height,
            )

        elif MOON_READ_MODE == "full_downsample":
            window = Window(
                col_off=0,
                row_off=0,
                width=src.width,
                height=src.height,
            )

        else:
            raise ValueError(
                "MOON_READ_MODE는 'center_crop' 또는 "
                "'full_downsample'이어야 합니다."
            )

        # 지구와 동일한 시각화 격자 크기 N_Y × N_X로 직접 읽는다.
        raw = src.read(
            1,
            window=window,
            out_shape=(N_Y, N_X),
            resampling=Resampling.bilinear,
            masked=False,
        ).astype(np.float64)

        nodata = src.nodata
        if nodata is not None and np.isfinite(nodata):
            raw[np.isclose(raw, nodata)] = np.nan

        # PDS4/GDAL metadata scale 및 offset 적용
        scale = float(src.scales[0]) if src.scales else 1.0
        offset = float(src.offsets[0]) if src.offsets else 0.0
        Z_abs = raw * scale + offset
        Z_abs[~np.isfinite(Z_abs)] = np.nan

        if not np.any(np.isfinite(Z_abs)):
            raise RuntimeError("달 DEM에서 유효한 고도값을 찾지 못했습니다.")

        # 시각화용 결측 보정
        if np.any(np.isnan(Z_abs)):
            Z_abs = Z_abs.copy()
            Z_abs[np.isnan(Z_abs)] = np.nanmean(Z_abs)

        # 읽은 window와 downsample 결과에 맞는 transform 생성
        base_transform = src.window_transform(window)
        display_transform = base_transform * rasterio.Affine.scale(
            float(window.width) / N_X,
            float(window.height) / N_Y,
        )

        rows, cols = np.indices((N_Y, N_X))
        x_values, y_values = transform_xy(
            display_transform,
            rows,
            cols,
            offset="center",
        )

        # rasterio.transform.xy()는 환경/버전에 따라
        # 1차원 좌표 배열을 반환할 수 있으므로 DEM shape로 복원한다.
        X = np.asarray(x_values, dtype=np.float64).reshape(N_Y, N_X)
        Y = np.asarray(y_values, dtype=np.float64).reshape(N_Y, N_X)

        # 지구 시각화와 동일하게 local distance 좌표로 변환
        X = X - X[0, 0]
        Y = Y - Y[0, 0]

        # shape map이 달 중심 반경값이어도 local terrain에서는 상대고도 사용
        Z_rel = Z_abs - np.nanmin(Z_abs)

        print("Moon PDS4 raster information")
        print("  Driver:", src.driver)
        print("  Native shape:", (src.height, src.width))
        print("  Display shape:", Z_rel.shape)
        print("  CRS:", src.crs)
        print("  Scale:", scale)
        print("  Offset:", offset)
        print(
            "  Window:",
            (
                int(window.col_off),
                int(window.row_off),
                int(window.width),
                int(window.height),
            ),
        )

    return X, Y, Z_rel, Z_abs, rows, cols




# ============================================================
# 4-3. 화성 DEM 영역 추출
# ============================================================


def extract_mars_area_patch(label_path, image_path):
    """
    PDS3 MOLA MEGDR의 LBL + IMG 파일을 읽어
    X, Y, Z_rel, Z_abs, rows, cols를 반환한다.

    처리 순서
    1) 우선 GDAL/rasterio의 PDS 드라이버로 LBL을 연다.
    2) 현재 GDAL 빌드가 PDS3를 지원하지 않으면 LBL을 직접 읽고
       IMG를 NumPy memmap으로 해석한다.
    3) 최종 좌표는 OpenSees에서 사용할 수 있도록 m 단위의
       지역 좌표계로 변환한다.
    """

    import re

    label_path = Path(label_path)
    image_path = Path(image_path)

    if not label_path.exists():
        raise FileNotFoundError(
            f"화성 PDS3 라벨 파일이 없습니다: {label_path.resolve()}"
        )

    if not image_path.exists():
        raise FileNotFoundError(
            f"화성 IMG 파일이 없습니다: {image_path.resolve()}"
        )

    print("Mars label path:", label_path.resolve())
    print("Mars image path:", image_path.resolve())

    try:
        import rasterio
        from rasterio.enums import Resampling
        from rasterio.windows import Window
        from rasterio.transform import xy as transform_xy
    except ImportError as exc:
        raise ImportError(
            "화성 DEM을 읽으려면 rasterio/GDAL이 필요합니다.\n"
            "설치 명령: conda install -c conda-forge rasterio gdal"
        ) from exc

    def make_window(width, height):
        """설정값에 따라 읽을 raster window를 생성한다."""
        if MARS_READ_MODE == "center_crop":
            crop_width = min(int(MARS_CROP_WIDTH_PIXELS), int(width))
            crop_height = min(int(MARS_CROP_HEIGHT_PIXELS), int(height))

            if MARS_CROP_COL_OFFSET is None:
                col_off = max((int(width) - crop_width) // 2, 0)
            else:
                col_off = max(
                    min(int(MARS_CROP_COL_OFFSET), int(width) - crop_width),
                    0,
                )

            if MARS_CROP_ROW_OFFSET is None:
                row_off = max((int(height) - crop_height) // 2, 0)
            else:
                row_off = max(
                    min(int(MARS_CROP_ROW_OFFSET), int(height) - crop_height),
                    0,
                )

            return Window(
                col_off=col_off,
                row_off=row_off,
                width=crop_width,
                height=crop_height,
            )

        if MARS_READ_MODE == "full_downsample":
            return Window(
                col_off=0,
                row_off=0,
                width=int(width),
                height=int(height),
            )

        raise ValueError(
            "MARS_READ_MODE는 'center_crop' 또는 "
            "'full_downsample'이어야 합니다."
        )

    def fill_invalid(values):
        values = np.asarray(values, dtype=np.float64)
        values[~np.isfinite(values)] = np.nan

        if not np.any(np.isfinite(values)):
            raise RuntimeError("화성 DEM에서 유효한 고도값을 찾지 못했습니다.")

        if np.any(np.isnan(values)):
            values = values.copy()
            values[np.isnan(values)] = np.nanmean(values)

        return values

    def coordinates_to_local_m(X_raw, Y_raw, crs):
        """
        raster 좌표를 지역 m 좌표로 변환한다.

        - geographic CRS: 화성 평균반지름으로 degree → m 변환
        - projected CRS: 좌표 단위가 m라고 간주
        - CRS 없음: 좌표 범위를 검사해 degree 가능성이 높으면 변환
        """
        mars_radius = 3_396_190.0  # m

        x_range = float(np.nanmax(X_raw) - np.nanmin(X_raw))
        y_range = float(np.nanmax(Y_raw) - np.nanmin(Y_raw))

        is_geographic = bool(
            crs is not None
            and getattr(crs, "is_geographic", False)
        )

        # PDS3 라벨을 GDAL이 읽어도 CRS가 비어 있을 수 있어
        # 좌표 범위가 경위도 범위이면 geographic으로 간주한다.
        looks_like_degrees = (
            x_range <= 360.0
            and y_range <= 180.0
            and np.nanmax(np.abs(Y_raw)) <= 90.0 + 1.0e-6
        )

        if is_geographic or (crs is None and looks_like_degrees):
            lon0 = float(X_raw[0, 0])
            lat0 = float(Y_raw[0, 0])
            lat0_rad = np.deg2rad(lat0)

            X_local = (
                mars_radius
                * np.cos(lat0_rad)
                * np.deg2rad(X_raw - lon0)
            )
            Y_local = mars_radius * np.deg2rad(Y_raw - lat0)
        else:
            X_local = X_raw - X_raw[0, 0]
            Y_local = Y_raw - Y_raw[0, 0]

        return X_local, Y_local

    # --------------------------------------------------------
    # 1) rasterio/GDAL PDS 드라이버로 LBL 열기
    # --------------------------------------------------------
    try:
        with rasterio.open(label_path) as src:
            if src.count < 1:
                raise RuntimeError(
                    "화성 PDS3 라벨에서 raster band를 찾지 못했습니다."
                )

            window = make_window(src.width, src.height)

            raw = src.read(
                1,
                window=window,
                out_shape=(N_Y, N_X),
                resampling=Resampling.bilinear,
                masked=False,
            ).astype(np.float64)

            nodata = src.nodata
            if nodata is not None and np.isfinite(nodata):
                raw[np.isclose(raw, nodata)] = np.nan

            scale = float(src.scales[0]) if src.scales else 1.0
            offset = float(src.offsets[0]) if src.offsets else 0.0
            Z_abs = fill_invalid(raw * scale + offset)

            base_transform = src.window_transform(window)
            display_transform = base_transform * rasterio.Affine.scale(
                float(window.width) / N_X,
                float(window.height) / N_Y,
            )

            rows, cols = np.indices((N_Y, N_X))
            x_values, y_values = transform_xy(
                display_transform,
                rows,
                cols,
                offset="center",
            )

            X_raw = np.asarray(
                x_values,
                dtype=np.float64,
            ).reshape(N_Y, N_X)
            Y_raw = np.asarray(
                y_values,
                dtype=np.float64,
            ).reshape(N_Y, N_X)

            X, Y = coordinates_to_local_m(X_raw, Y_raw, src.crs)
            Z_rel = Z_abs - np.nanmin(Z_abs)

            print("Mars raster information")
            print("  Reader: rasterio/GDAL PDS driver")
            print("  Driver:", src.driver)
            print("  Native shape:", (src.height, src.width))
            print("  Display shape:", Z_rel.shape)
            print("  CRS:", src.crs)
            print("  Transform:", src.transform)
            print("  Scale:", scale)
            print("  Offset:", offset)
            print(
                "  Window:",
                (
                    int(window.col_off),
                    int(window.row_off),
                    int(window.width),
                    int(window.height),
                ),
            )

            return X, Y, Z_rel, Z_abs, rows, cols

    except Exception as rasterio_error:
        print(
            "[Warning] rasterio가 PDS3 LBL을 직접 열지 못했습니다."
        )
        print("          NumPy 기반 PDS3 fallback reader를 사용합니다.")
        print("          rasterio error:", rasterio_error)

    # --------------------------------------------------------
    # 2) PDS3 LBL 직접 파싱 + IMG memmap fallback
    # --------------------------------------------------------
    label_text = label_path.read_text(
        encoding="ascii",
        errors="ignore",
    )

    def get_label_value(key, required=True, default=None):
        pattern = rf"(?mi)^\s*{re.escape(key)}\s*=\s*(.+?)\s*$"
        match = re.search(pattern, label_text)

        if match is None:
            if required:
                raise RuntimeError(
                    f"화성 LBL에서 필수 항목 '{key}'를 찾지 못했습니다."
                )
            return default

        value = match.group(1).strip()
        # 줄 끝 주석 제거
        value = value.split("/*", 1)[0].strip()
        return value

    def parse_number(value):
        if value is None:
            return None
        cleaned = re.sub(r"<[^>]+>", "", str(value))
        cleaned = cleaned.strip().strip('"').strip("'")
        return float(cleaned)

    lines = int(parse_number(get_label_value("LINES")))
    samples = int(parse_number(get_label_value("LINE_SAMPLES")))
    sample_bits = int(parse_number(get_label_value("SAMPLE_BITS")))
    sample_type = get_label_value("SAMPLE_TYPE").upper()

    if sample_bits not in (8, 16, 32, 64):
        raise RuntimeError(
            f"지원하지 않는 SAMPLE_BITS입니다: {sample_bits}"
        )

    if "REAL" in sample_type or "FLOAT" in sample_type:
        kind = "f"
    elif "UNSIGNED" in sample_type:
        kind = "u"
    else:
        kind = "i"

    if "LSB" in sample_type or "PC" in sample_type:
        endian = "<"
    elif sample_bits == 8:
        endian = "|"
    else:
        endian = ">"

    dtype = np.dtype(f"{endian}{kind}{sample_bits // 8}")

    scaling_factor = parse_number(
        get_label_value(
            "SCALING_FACTOR",
            required=False,
            default="1.0",
        )
    )
    value_offset = parse_number(
        get_label_value(
            "OFFSET",
            required=False,
            default="0.0",
        )
    )

    missing_value = get_label_value(
        "MISSING_CONSTANT",
        required=False,
        default=None,
    )
    if missing_value is None:
        missing_value = get_label_value(
            "CORE_NULL",
            required=False,
            default=None,
        )
    missing_value = (
        parse_number(missing_value)
        if missing_value is not None
        else None
    )

    expected_bytes = lines * samples * dtype.itemsize
    actual_bytes = image_path.stat().st_size

    if actual_bytes < expected_bytes:
        raise RuntimeError(
            "화성 IMG 파일 크기가 LBL 정보보다 작습니다.\n"
            f"LBL 예상 크기: {expected_bytes:,} bytes\n"
            f"실제 IMG 크기: {actual_bytes:,} bytes"
        )

    raster = np.memmap(
        image_path,
        dtype=dtype,
        mode="r",
        shape=(lines, samples),
    )

    window = make_window(samples, lines)

    row_start = int(window.row_off)
    row_stop = row_start + int(window.height)
    col_start = int(window.col_off)
    col_stop = col_start + int(window.width)

    # 전체 crop을 메모리에 복사하지 않고 최종 출력 위치만 샘플링한다.
    row_indices_native = np.linspace(
        row_start,
        row_stop - 1,
        N_Y,
    ).round().astype(int)
    col_indices_native = np.linspace(
        col_start,
        col_stop - 1,
        N_X,
    ).round().astype(int)

    raw = np.asarray(
        raster[np.ix_(row_indices_native, col_indices_native)],
        dtype=np.float64,
    )

    if missing_value is not None:
        raw[np.isclose(raw, missing_value)] = np.nan

    Z_abs = fill_invalid(raw * scaling_factor + value_offset)
    Z_rel = Z_abs - np.nanmin(Z_abs)

    # MAP_SCALE이 있으면 우선 사용한다.
    map_scale_text = get_label_value(
        "MAP_SCALE",
        required=False,
        default=None,
    )
    map_resolution_text = get_label_value(
        "MAP_RESOLUTION",
        required=False,
        default=None,
    )

    pixel_size_m = None

    if map_scale_text is not None:
        map_scale = parse_number(map_scale_text)
        if "KM" in map_scale_text.upper():
            pixel_size_m = map_scale * 1000.0
        else:
            pixel_size_m = map_scale

    if pixel_size_m is None and map_resolution_text is not None:
        pixels_per_degree = parse_number(map_resolution_text)
        mars_radius = 3_396_190.0
        pixel_size_m = (
            2.0 * np.pi * mars_radius / 360.0 / pixels_per_degree
        )

    if pixel_size_m is None:
        raise RuntimeError(
            "LBL에서 MAP_SCALE 또는 MAP_RESOLUTION을 찾지 못해 "
            "화성 수평좌표를 m 단위로 생성할 수 없습니다."
        )

    x_native = (col_indices_native - col_indices_native[0]) * pixel_size_m
    y_native = (row_indices_native - row_indices_native[0]) * pixel_size_m

    X, Y = np.meshgrid(x_native, y_native)
    rows, cols = np.indices((N_Y, N_X))

    print("Mars raster information")
    print("  Reader: NumPy PDS3 fallback")
    print("  Native shape:", (lines, samples))
    print("  Display shape:", Z_rel.shape)
    print("  SAMPLE_TYPE:", sample_type)
    print("  SAMPLE_BITS:", sample_bits)
    print("  dtype:", dtype)
    print("  Scale:", scaling_factor)
    print("  Offset:", value_offset)
    print("  Pixel size:", pixel_size_m, "m")
    print(
        "  Window:",
        (
            int(window.col_off),
            int(window.row_off),
            int(window.width),
            int(window.height),
        ),
    )

    return X, Y, Z_rel, Z_abs, rows, cols

# ============================================================
# 5. OpenSeesPy 2D 지반 모델 생성
# ============================================================

def build_opensees_soil_model(x, ground_y, config, run_gravity=True):
    """
    DEM 지표면을 상부 경계로 하는 2D plane strain 지반 모델 생성.

    반영 사항:
    1. config에 따라 지구, 화성 또는 달 재료와 중력가속도를 선택
    2. ElasticIsotropic 및 PressureDependMultiYield02 재료 사용
    3. PDMY02 재료는 중력 재하 단계에서 stage 0, 이후 stage 1로 전환
    4. 자중은 pattern Plain 기반의 등가 절점하중으로 적용
    """

    ops.wipe()

    # 2D plane strain, 각 노드 자유도 ux, uy
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    body_name = config["body_name"]
    gravity = float(config["gravity"])
    materials = config["materials"]
    material_order = config["material_order"]
    assign_material = config["assign_material"]
    base_depth = float(config["base_depth"])

    # 재료 정의
    nonlinear_material_tags = register_opensees_soil_materials(
        materials,
        material_order,
    )

    nx = len(x)
    nz = N_DEPTH

    base_y = np.nanmin(ground_y) - base_depth

    node_tags = {}
    node_tag = 1

    # --------------------------------------------------------
    # 1) 노드 생성
    # --------------------------------------------------------

    for i in range(nx):
        for j in range(nz + 1):
            ratio = j / nz

            # j = 0   : 지표면
            # j = nz  : 하부 경계
            y = (1.0 - ratio) * ground_y[i] + ratio * base_y

            ops.node(node_tag, float(x[i]), float(y))
            node_tags[(i, j)] = node_tag
            node_tag += 1

    # --------------------------------------------------------
    # 2) 요소 생성 + 요소별 재료 배정
    # --------------------------------------------------------

    ele_tag = 1
    ele_materials = {}

    # 자중을 pattern Plain 안의 nodal load로 넣기 위한 저장소
    nodal_self_weight_loads = {}

    x_min = float(np.nanmin(x))
    x_max = float(np.nanmax(x))
    g_min = float(np.nanmin(ground_y))
    g_max = float(np.nanmax(ground_y))

    for i in range(nx - 1):
        for j in range(nz):
            # quad 요소 노드 순서: 반시계 방향
            n1 = node_tags[(i, j + 1)]      # bottom-left
            n2 = node_tags[(i + 1, j + 1)]  # bottom-right
            n3 = node_tags[(i + 1, j)]      # top-right
            n4 = node_tags[(i, j)]          # top-left

            # 요소 중심 깊이 계산
            x_center = 0.5 * (x[i] + x[i + 1])
            ground_center = 0.5 * (ground_y[i] + ground_y[i + 1])

            ratio_mid = (j + 0.5) / nz
            y_center = (1.0 - ratio_mid) * ground_center + ratio_mid * base_y

            depth_from_surface = ground_center - y_center

            x_norm = (x_center - x_min) / (x_max - x_min + 1.0e-12)
            elevation_norm = (ground_center - g_min) / (g_max - g_min + 1.0e-12)

            mat_name = assign_material(
                depth_from_surface=depth_from_surface,
                x_norm=x_norm,
                y_norm=0.0,
                elevation_norm=elevation_norm,
            )

            mat = materials[mat_name]
            mat_tag = int(mat["tag"])
            rho = float(mat["rho"])

            # 여기서는 quad의 body force를 0으로 둔다.
            # 이유:
            #   자중을 pattern Plain 안에서 명시적으로 적용하기 위해서.
            #   body force와 nodal load를 동시에 넣으면 자중이 중복될 수 있다.
            ops.element(
                "quad",
                ele_tag,
                n1, n2, n3, n4,
                THICKNESS,
                "PlaneStrain",
                mat_tag,
                0.0,
                0.0,
                0.0,
                0.0,
            )

            # ------------------------------------------------
            # 요소 자중을 등가 nodal load로 변환
            # ------------------------------------------------
            # 2D plane strain 단면에서 요소 면적 * 두께 * 밀도 * 중력
            # 현재 mesh는 거의 사각형이지만, DEM 지형 때문에 일반 quad로 보고
            # polygon area 공식으로 면적을 계산한다.

            coords = np.array([
                ops.nodeCoord(n1),
                ops.nodeCoord(n2),
                ops.nodeCoord(n3),
                ops.nodeCoord(n4),
            ], dtype=float)

            x_coords = coords[:, 0]
            y_coords = coords[:, 1]

            area = 0.5 * abs(
                np.dot(x_coords, np.roll(y_coords, -1))
                - np.dot(y_coords, np.roll(x_coords, -1))
            )

            element_weight = rho * gravity * area * THICKNESS

            # 네 개 노드에 균등 분배
            nodal_weight = -element_weight / 4.0

            for nd in [n1, n2, n3, n4]:
                nodal_self_weight_loads[nd] = nodal_self_weight_loads.get(nd, 0.0) + nodal_weight

            ele_materials[ele_tag] = {
                "material": mat_name,
                "material_tag": mat_tag,
                "model": mat["model"],
                "rho": rho,
                "depth_m": depth_from_surface,
                "x_center_m": x_center,
                "area_m2": area,
            }

            ele_tag += 1

    # --------------------------------------------------------
    # 3) 경계조건
    # --------------------------------------------------------

    # 바닥: x, y 고정
    for i in range(nx):
        bottom_node = node_tags[(i, nz)]
        ops.fix(bottom_node, 1, 1)

    # 좌우 측면: x 방향 고정, y 방향 자유
    for j in range(nz):
        left_node = node_tags[(0, j)]
        right_node = node_tags[(nx - 1, j)]

        ops.fix(left_node, 1, 0)
        ops.fix(right_node, 1, 0)

    # --------------------------------------------------------
    # 4) 탄성 중력 재하 stage 0 → loadConst → 소성 stage 1
    # --------------------------------------------------------

    if run_gravity:
        run_elastic_gravity_then_plastic_stage(
            nodal_self_weight_loads=nodal_self_weight_loads,
            nonlinear_material_tags=nonlinear_material_tags,
            n_steps=10,
            d_lambda=0.1,
        )

    print(f"{body_name} OpenSees soil model generated.")
    print("Number of nodes:", len(ops.getNodeTags()))
    print("Number of elements:", len(ops.getEleTags()))
    print("Nonlinear PDMY/PDMY02 material tags:", nonlinear_material_tags)

    return node_tags, ele_materials, nodal_self_weight_loads, nonlinear_material_tags


# ============================================================
# 6. 고유진동수 해석
# ============================================================

def run_modal_analysis(num_modes=3):
    """
    생성된 OpenSees 모델에 대해 고유치 해석을 수행한다.
    """
    try:
        eigen_values = ops.eigen(num_modes)

        freqs = []

        for lam in eigen_values:
            omega = np.sqrt(lam)
            freq = omega / (2 * np.pi)
            freqs.append(freq)

        return freqs

    except Exception as e:
        print("Modal analysis failed:", e)
        return []


# ============================================================
# 7. 모델 시각화
# ============================================================

def plot_profile(x, ground_y, title="Extracted DEM Profile"):
    plt.figure(figsize=(10, 4))
    plt.plot(x, ground_y, marker="o", markersize=3)
    plt.xlabel("Distance along profile x (m)")
    plt.ylabel("Relative elevation (m)")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_mesh(x, ground_y):
    base_y = np.nanmin(ground_y) - BASE_DEPTH
    nz = N_DEPTH

    plt.figure(figsize=(11, 5))

    # 수직선
    for i in range(len(x)):
        ys = []
        xs = []

        for j in range(nz + 1):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            ys.append(y)
            xs.append(x[i])

        plt.plot(xs, ys, linewidth=0.5)

    # 수평선
    for j in range(nz + 1):
        xs = []
        ys = []

        for i in range(len(x)):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            xs.append(x[i])
            ys.append(y)

        plt.plot(xs, ys, linewidth=0.5)

    plt.plot(x, ground_y, linewidth=2, label="DEM ground surface")

    plt.xlabel("Distance x (m)")
    plt.ylabel("Elevation y (m)")
    plt.title("OpenSees 2D Soil Mesh from NASADEM")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
# ============================================================
# 7-4. Plotly HTML 3D 지형 영역 저장
# ============================================================

def save_area_3d_html(X, Y, Z_rel, body_name="Earth"):
    """
    DEM 면적 영역을 3D surface로 시각화한다.
    정사각형 또는 3:4 직사각형 형태의 지형 확인용.
    """

    fig = go.Figure()

    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z_rel,
            colorscale="Earth",
            opacity=0.95,
            colorbar=dict(title="Relative<br>Elevation (m)"),
            name="DEM surface",
        )
    )

    fig.update_layout(
        title=f"{body_name} DEM 3D Terrain Area",
        template="plotly_white",
        width=1100,
        height=850,
        scene=dict(
            xaxis_title="East-West distance x (m)",
            yaxis_title="North-South distance y (m)",
            zaxis_title="Relative elevation z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.5, y=-1.8, z=1.1),
            ),
        ),
    )

    html_path = Path(f"{body_name.lower()}_area_3d_terrain.html").resolve()
    fig.write_html(str(html_path), auto_open=False)

    print(f"Saved: {html_path}")
    return html_path
    

# ============================================================
# 7-5. Plotly HTML 3D 지반 Solid Block 저장
# ============================================================

def save_area_3d_solid_html(X, Y, Z_rel, body_name="Earth"):
    """
    DEM 면적 영역을 두께가 있는 지반 block처럼 시각화한다.
    Abaqus viewport에서 solid part를 보는 느낌에 가까움.
    """

    base_z = np.nanmin(Z_rel) - BASE_DEPTH

    ny, nx = Z_rel.shape

    x_flat = []
    y_flat = []
    z_flat = []

    top_node = {}
    bottom_node = {}

    node_id = 0

    # top nodes
    for iy in range(ny):
        for ix in range(nx):
            x_flat.append(float(X[iy, ix]))
            y_flat.append(float(Y[iy, ix]))
            z_flat.append(float(Z_rel[iy, ix]))
            top_node[(iy, ix)] = node_id
            node_id += 1

    # bottom nodes
    for iy in range(ny):
        for ix in range(nx):
            x_flat.append(float(X[iy, ix]))
            y_flat.append(float(Y[iy, ix]))
            z_flat.append(float(base_z))
            bottom_node[(iy, ix)] = node_id
            node_id += 1

    I = []
    J = []
    K = []
    intensity = []

    def add_quad(n1, n2, n3, n4, value):
        I.append(n1)
        J.append(n2)
        K.append(n3)
        intensity.append(value)

        I.append(n1)
        J.append(n3)
        K.append(n4)
        intensity.append(value)

    # top surface
    for iy in range(ny - 1):
        for ix in range(nx - 1):
            n1 = top_node[(iy, ix)]
            n2 = top_node[(iy, ix + 1)]
            n3 = top_node[(iy + 1, ix + 1)]
            n4 = top_node[(iy + 1, ix)]

            z_mean = np.mean([z_flat[n1], z_flat[n2], z_flat[n3], z_flat[n4]])
            add_quad(n1, n2, n3, n4, z_mean)

    # bottom surface
    for iy in range(ny - 1):
        for ix in range(nx - 1):
            n1 = bottom_node[(iy, ix)]
            n2 = bottom_node[(iy + 1, ix)]
            n3 = bottom_node[(iy + 1, ix + 1)]
            n4 = bottom_node[(iy, ix + 1)]

            add_quad(n1, n2, n3, n4, base_z)

    # front/back/left/right side surfaces
    # y-min side
    iy = 0
    for ix in range(nx - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = bottom_node[(iy, ix + 1)]
        n3 = top_node[(iy, ix + 1)]
        n4 = top_node[(iy, ix)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n3], z_flat[n4]]))

    # y-max side
    iy = ny - 1
    for ix in range(nx - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = top_node[(iy, ix)]
        n3 = top_node[(iy, ix + 1)]
        n4 = bottom_node[(iy, ix + 1)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n2], z_flat[n3]]))

    # x-min side
    ix = 0
    for iy in range(ny - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = top_node[(iy, ix)]
        n3 = top_node[(iy + 1, ix)]
        n4 = bottom_node[(iy + 1, ix)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n2], z_flat[n3]]))

    # x-max side
    ix = nx - 1
    for iy in range(ny - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = bottom_node[(iy + 1, ix)]
        n3 = top_node[(iy + 1, ix)]
        n4 = top_node[(iy, ix)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n3], z_flat[n4]]))

    fig = go.Figure()

    fig.add_trace(
        go.Mesh3d(
            x=x_flat,
            y=y_flat,
            z=z_flat,
            i=I,
            j=J,
            k=K,
            intensity=intensity,
            colorscale="Earth",
            opacity=0.95,
            showscale=True,
            colorbar=dict(title="Elevation<br>(m)"),
            name="Soil solid block",
        )
    )

    # element edge 느낌의 선 추가
    edge_x = []
    edge_y = []
    edge_z = []

    def add_edge(n1, n2):
        edge_x.extend([x_flat[n1], x_flat[n2], None])
        edge_y.extend([y_flat[n1], y_flat[n2], None])
        edge_z.extend([z_flat[n1], z_flat[n2], None])

    # top grid edges
    for iy in range(ny):
        for ix in range(nx - 1):
            add_edge(top_node[(iy, ix)], top_node[(iy, ix + 1)])

    for ix in range(nx):
        for iy in range(ny - 1):
            add_edge(top_node[(iy, ix)], top_node[(iy + 1, ix)])

    # vertical side edges만 일부 표시
    for iy in [0, ny - 1]:
        for ix in range(nx):
            add_edge(top_node[(iy, ix)], bottom_node[(iy, ix)])

    for ix in [0, nx - 1]:
        for iy in range(ny):
            add_edge(top_node[(iy, ix)], bottom_node[(iy, ix)])

    fig.add_trace(
        go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            line=dict(color="black", width=2),
            name="Mesh edges",
            hoverinfo="skip",
        )
    )

    fig.update_layout(
        title=f"Abaqus-like 3D Solid Terrain Block — {body_name} DEM",
        template="plotly_white",
        width=1200,
        height=850,
        scene=dict(
            xaxis_title="East-West distance x (m)",
            yaxis_title="North-South distance y (m)",
            zaxis_title="Elevation z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.4, y=-1.7, z=1.0),
            ),
        ),
    )

    html_path = Path(f"{body_name.lower()}_area_3d_solid.html").resolve()
    fig.write_html(str(html_path), auto_open=False)

    print(f"Saved: {html_path}")
    return html_path
    
    
# ============================================================
# 7-6. Plotly HTML 3D 재료 분포 저장
# ============================================================

def save_area_3d_material_distribution_html(X, Y, Z_rel):
    """
    DEM 기반 지반 block 내부의 재료 분포를 색깔별로 시각화한다.

    주의:
    - 이 함수는 해석 모델이 아니라 시각화용이다.
    - OpenSees의 material tag가 자동으로 색을 만드는 게 아니므로,
      assign_soil_material() 결과를 Plotly에 따로 반영해야 한다.
    """

    base_z = np.nanmin(Z_rel) - BASE_DEPTH

    ny, nx = Z_rel.shape
    nz = N_DEPTH

    fig = go.Figure()

    material_points = {}

    for name in SOIL_ORDER:
        material_points[name] = {
            "x": [],
            "y": [],
            "z": [],
            "text": [],
        }

    z_min = float(np.nanmin(Z_rel))
    z_max = float(np.nanmax(Z_rel))

    for iy in range(ny):
        for ix in range(nx):
            x = float(X[iy, ix])
            y = float(Y[iy, ix])
            ground_z = float(Z_rel[iy, ix])

            x_norm = ix / max(nx - 1, 1)
            y_norm = iy / max(ny - 1, 1)
            elevation_norm = (ground_z - z_min) / (z_max - z_min + 1.0e-12)

            for kz in range(nz):
                ratio_mid = (kz + 0.5) / nz
                depth = ratio_mid * BASE_DEPTH

                z = (1.0 - ratio_mid) * ground_z + ratio_mid * base_z

                mat_name = assign_soil_material(
                    depth_from_surface=depth,
                    x_norm=x_norm,
                    y_norm=y_norm,
                    elevation_norm=elevation_norm,
                )

                mat = SOIL_MATERIALS[mat_name]

                material_points[mat_name]["x"].append(x)
                material_points[mat_name]["y"].append(y)
                material_points[mat_name]["z"].append(float(z))
                material_points[mat_name]["text"].append(
                    f"Material: {mat['label']}<br>"
                    f"Depth: {depth:.2f} m<br>"
                    f"Model: {mat['model']}<br>"
                    f"rho: {mat['rho']:.0f} kg/m³"
                )

    for name in SOIL_ORDER:
        pts = material_points[name]

        if len(pts["x"]) == 0:
            continue

        mat = SOIL_MATERIALS[name]

        fig.add_trace(
            go.Scatter3d(
                x=pts["x"],
                y=pts["y"],
                z=pts["z"],
                mode="markers",
                marker=dict(
                    size=3.5,
                    color=mat["color"],
                    opacity=0.75,
                ),
                name=mat["label"],
                text=pts["text"],
                hovertemplate="%{text}<extra></extra>",
            )
        )

    # 지표면 추가
    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z_rel,
            colorscale="Greys",
            opacity=0.35,
            showscale=False,
            name="DEM ground surface",
            hovertemplate="Ground surface<extra></extra>",
        )
    )

    fig.update_layout(
        title=(
            "Pohang Soil Material Distribution"
            "<br><sup>Organic / Clay / Silt / Sand / Gravel / Weathered Rock / Rock</sup>"
        ),
        template="plotly_white",
        width=1250,
        height=850,
        scene=dict(
            xaxis_title="East-West distance x (m)",
            yaxis_title="North-South distance y (m)",
            zaxis_title="Elevation / Depth z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.6, y=-1.8, z=1.1),
            ),
        ),
        legend=dict(
            title="Soil material",
            x=0.02,
            y=0.98,
        ),
    )

    html_path = Path("pohang_area_3d_soil_material_distribution.html").resolve()
    fig.write_html(str(html_path), auto_open=True)

    print(f"Saved material distribution HTML: {html_path}")
    webbrowser.open(html_path.as_uri())
    

# ============================================================
# 8. 단면 좌표 저장
# ============================================================

def save_profile_csv(x, ground_y, elevations_abs, lats, lons):
    data = np.column_stack([x, ground_y, elevations_abs, lats, lons])

    header = "x_m,relative_elevation_m,absolute_elevation_m,latitude,longitude"

    np.savetxt(
        "pohang_nasadem_profile.csv",
        data,
        delimiter=",",
        header=header,
        comments="",
    )

    print("Saved: pohang_nasadem_profile.csv")


# ============================================================
# 9. 메인 실행
# ============================================================

def main():
    body = TARGET_BODY.strip().lower()

    if body == "earth":
        print("Reading Earth HGT file...")
        dem = read_hgt_file(HGT_PATH)

        print("DEM shape:", dem.shape)
        print("DEM elevation min:", np.nanmin(dem))
        print("DEM elevation max:", np.nanmax(dem))

        print("Extracting Earth area patch...")
        (
            X_area,
            Y_area,
            Z_area_rel,
            Z_area_abs,
            area_lats,
            area_lons,
        ) = extract_area_patch(dem)

        body_name = "Earth"

    elif body == "mars":
        print("Reading Mars DEM...")
        (
            X_area,
            Y_area,
            Z_area_rel,
            Z_area_abs,
            area_rows,
            area_cols,
        ) = extract_mars_area_patch(
            MARS_LABEL_PATH,
            MARS_IMAGE_PATH,
        )

        body_name = "Mars"

    elif body == "moon":
        print("Reading Moon PDS4 XML + IMG...")
        (
            X_area,
            Y_area,
            Z_area_rel,
            Z_area_abs,
            area_rows,
            area_cols,
        ) = extract_lunar_area_patch(XML_PATH, IMG_PATH)

        body_name = "Moon"

    else:
        raise ValueError(
            "TARGET_BODY는 'Earth', 'Mars', 'Moon' 중 하나여야 합니다."
        )

    print("\nCommon terrain grid")
    print("  Body:", body_name)
    selected_config = {"Earth": EARTH_CONFIG, "Mars": MARS_CONFIG, "Moon": MOON_CONFIG}[body_name]
    print("  Gravity:", selected_config["gravity"], "m/s^2")
    print("  Shape:", Z_area_rel.shape)
    print(
        "  X length:",
        float(np.nanmax(X_area) - np.nanmin(X_area)),
        "m",
    )
    print(
        "  Y length:",
        float(np.nanmax(Y_area) - np.nanmin(Y_area)),
        "m",
    )
    print(
        "  Relative elevation range:",
        float(np.nanmin(Z_area_rel)),
        "~",
        float(np.nanmax(Z_area_rel)),
        "m",
    )

    mode = VISUALIZATION_MODE.strip().lower()

    if mode == "surface":
        print("\nSaving one 3D surface HTML...")
        html_path = save_area_3d_html(
            X_area,
            Y_area,
            Z_area_rel,
            body_name=body_name,
        )

    elif mode == "solid":
        print("\nSaving one 3D solid-block HTML...")
        html_path = save_area_3d_solid_html(
            X_area,
            Y_area,
            Z_area_rel,
            body_name=body_name,
        )

    else:
        raise ValueError(
            "VISUALIZATION_MODE는 'surface' 또는 'solid'이어야 합니다."
        )

    print("\nVisualization created:")
    print(" ", html_path)

    if AUTO_OPEN_HTML:
        webbrowser.open(html_path.as_uri())

    # 지구 또는 달 설정을 선택하여 OpenSees 지반모델 생성
    if BUILD_OPENSEES_MODEL:
        print(f"\nBuilding {body_name} OpenSees 2D soil model...")

        center_iy = Z_area_rel.shape[0] // 2
        x_profile = X_area[center_iy, :]
        ground_profile = Z_area_rel[center_iy, :]

        if body_name == "Earth":
            soil_config = EARTH_CONFIG
        elif body_name == "Mars":
            soil_config = MARS_CONFIG
        else:
            soil_config = MOON_CONFIG

        build_opensees_soil_model(
            x_profile,
            ground_profile,
            config=soil_config,
            run_gravity=True,
        )

        print(f"{body_name} OpenSees soil model is ready.")





# ============================================================
# 10. 시각화
# ============================================================

def save_profile_html(x, ground_y, elevations_abs, lats, lons):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=x,
            y=ground_y,
            mode="lines+markers",
            name="DEM ground surface",
            text=[
                f"lat={lat:.5f}<br>lon={lon:.5f}<br>abs elev={elev:.2f} m"
                for lat, lon, elev in zip(lats, lons, elevations_abs)
            ],
            hovertemplate=(
                "x = %{x:.2f} m<br>"
                "relative elevation = %{y:.2f} m<br>"
                "%{text}<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title="Pohang Heunghae-eup NASADEM 2D Ground Profile",
        xaxis_title="Distance along profile x (m)",
        yaxis_title="Relative elevation (m)",
        template="plotly_white",
        width=1100,
        height=500,
    )

    html_path = Path("pohang_profile_2d.html").resolve()

    fig.write_html(str(html_path), auto_open=False)

    print(f"Saved: {html_path}")


# ============================================================
# 7-2. Plotly HTML 2D Mesh 저장
# ============================================================

def save_mesh_2d_html(x, ground_y):
    base_y = np.nanmin(ground_y) - BASE_DEPTH
    nz = N_DEPTH

    fig = go.Figure()

    # 수직 mesh line
    for i in range(len(x)):
        xs = []
        ys = []

        for j in range(nz + 1):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            xs.append(x[i])
            ys.append(y)

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(width=1, color="gray"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    # 수평 mesh line
    for j in range(nz + 1):
        xs = []
        ys = []

        for i in range(len(x)):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            xs.append(x[i])
            ys.append(y)

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(width=1, color="gray"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    # 지표면 강조
    fig.add_trace(
        go.Scatter(
            x=x,
            y=ground_y,
            mode="lines+markers",
            line=dict(width=4, color="black"),
            marker=dict(size=4),
            name="DEM ground surface",
        )
    )

    fig.update_layout(
        title="OpenSees 2D Soil Mesh from NASADEM",
        xaxis_title="Distance x (m)",
        yaxis_title="Elevation y (m)",
        template="plotly_white",
        width=1100,
        height=600,
        yaxis=dict(scaleanchor="x", scaleratio=1),
    )

    html_path = Path("pohang_soil_mesh_2d.html").resolve()

    fig.write_html(str(html_path), auto_open=False)

    print(f"Saved: {html_path}")


# ============================================================
# 7-3. Abaqus 느낌의 3D Solid Mesh HTML 저장
# ============================================================

def save_mesh_3d_html(x, ground_y):
    """
    2D plane strain 지반 단면을 시각화용 두께 방향으로 extrusion하여
    Abaqus viewport처럼 보이는 3D solid mesh HTML을 생성한다.

    주의:
    - 해석 모델은 여전히 2D plane strain이다.
    - 여기서 y 방향 두께는 시각화용 가상 두께이다.
    """
    base_z = np.nanmin(ground_y) - BASE_DEPTH
    nx = len(x)
    nz = N_DEPTH

    thickness_vis = 80.0  # 시각화용 두께, 실제 해석 두께 아님
    y_front = -thickness_vis / 2.0
    y_back = thickness_vis / 2.0

    # --------------------------------------------------------
    # 1) 3D 시각화용 노드 생성
    # 좌표계:
    #   X = 단면 진행 방향
    #   Y = 시각화용 두께 방향
    #   Z = 고도/깊이 방향
    # --------------------------------------------------------

    X = []
    Y = []
    Z = []

    front_node = {}
    back_node = {}

    node_id = 0

    for side_name, y_side, node_map in [
        ("front", y_front, front_node),
        ("back", y_back, back_node),
    ]:
        for i in range(nx):
            for j in range(nz + 1):
                ratio = j / nz
                z = (1 - ratio) * ground_y[i] + ratio * base_z

                X.append(float(x[i]))
                Y.append(float(y_side))
                Z.append(float(z))

                node_map[(i, j)] = node_id
                node_id += 1

    # --------------------------------------------------------
    # 2) 표면 face 생성
    # Mesh3d는 삼각형 face만 받으므로 quad를 삼각형 2개로 분할
    # --------------------------------------------------------

    I = []
    J = []
    K = []
    intensity = []

    def add_quad(n1, n2, n3, n4, value):
        # quad: n1-n2-n3-n4
        # triangle 1: n1-n2-n3
        I.append(n1)
        J.append(n2)
        K.append(n3)
        intensity.append(value)

        # triangle 2: n1-n3-n4
        I.append(n1)
        J.append(n3)
        K.append(n4)
        intensity.append(value)

    # front/back face
    for i in range(nx - 1):
        for j in range(nz):
            # front face
            f1 = front_node[(i, j + 1)]
            f2 = front_node[(i + 1, j + 1)]
            f3 = front_node[(i + 1, j)]
            f4 = front_node[(i, j)]

            z_mean_front = np.mean([Z[f1], Z[f2], Z[f3], Z[f4]])
            add_quad(f1, f2, f3, f4, z_mean_front)

            # back face
            b1 = back_node[(i, j + 1)]
            b2 = back_node[(i + 1, j + 1)]
            b3 = back_node[(i + 1, j)]
            b4 = back_node[(i, j)]

            z_mean_back = np.mean([Z[b1], Z[b2], Z[b3], Z[b4]])
            add_quad(b4, b3, b2, b1, z_mean_back)

    # top surface
    for i in range(nx - 1):
        f1 = front_node[(i, 0)]
        f2 = front_node[(i + 1, 0)]
        b2 = back_node[(i + 1, 0)]
        b1 = back_node[(i, 0)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, f2, b2, b1, z_mean)

    # bottom surface
    for i in range(nx - 1):
        f1 = front_node[(i, nz)]
        f2 = front_node[(i + 1, nz)]
        b2 = back_node[(i + 1, nz)]
        b1 = back_node[(i, nz)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, b1, b2, f2, z_mean)

    # left side surface
    for j in range(nz):
        f1 = front_node[(0, j + 1)]
        f2 = front_node[(0, j)]
        b2 = back_node[(0, j)]
        b1 = back_node[(0, j + 1)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, f2, b2, b1, z_mean)

    # right side surface
    for j in range(nz):
        f1 = front_node[(nx - 1, j + 1)]
        f2 = front_node[(nx - 1, j)]
        b2 = back_node[(nx - 1, j)]
        b1 = back_node[(nx - 1, j + 1)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, b1, b2, f2, z_mean)

    fig = go.Figure()

    # --------------------------------------------------------
    # 3) Solid surface 추가
    # --------------------------------------------------------

    fig.add_trace(
        go.Mesh3d(
            x=X,
            y=Y,
            z=Z,
            i=I,
            j=J,
            k=K,
            intensity=intensity,
            colorscale="Earth",
            opacity=0.92,
            name="Soil solid",
            showscale=True,
            colorbar=dict(title="Elevation<br>(m)"),
            hovertemplate=(
                "x = %{x:.2f} m<br>"
                "y = %{y:.2f} m<br>"
                "z = %{z:.2f} m<br>"
                "<extra></extra>"
            ),
        )
    )

    # --------------------------------------------------------
    # 4) Mesh line 추가
    # Abaqus의 element edge처럼 보이게 하는 부분
    # --------------------------------------------------------

    edge_x = []
    edge_y = []
    edge_z = []

    def add_edge(n1, n2):
        edge_x.extend([X[n1], X[n2], None])
        edge_y.extend([Y[n1], Y[n2], None])
        edge_z.extend([Z[n1], Z[n2], None])

    # front/back mesh lines
    for node_map in [front_node, back_node]:
        # vertical lines
        for i in range(nx):
            for j in range(nz):
                add_edge(node_map[(i, j)], node_map[(i, j + 1)])

        # horizontal lines
        for j in range(nz + 1):
            for i in range(nx - 1):
                add_edge(node_map[(i, j)], node_map[(i + 1, j)])

    # thickness direction connection lines
    for i in range(nx):
        for j in range(nz + 1):
            add_edge(front_node[(i, j)], back_node[(i, j)])

    fig.add_trace(
        go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            line=dict(color="black", width=2),
            name="Element edges",
            hoverinfo="skip",
        )
    )

    # --------------------------------------------------------
    # 5) 지표면 선 강조
    # --------------------------------------------------------

    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=np.full_like(x, y_front),
            z=ground_y,
            mode="lines",
            line=dict(color="red", width=6),
            name="Front ground surface",
        )
    )

    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=np.full_like(x, y_back),
            z=ground_y,
            mode="lines",
            line=dict(color="red", width=6),
            name="Back ground surface",
            showlegend=False,
        )
    )

    # --------------------------------------------------------
    # 6) Layout
    # --------------------------------------------------------

    fig.update_layout(
        title="Abaqus-like 3D Visualization of NASADEM-based OpenSees Soil Mesh",
        template="plotly_white",
        width=1200,
        height=800,
        scene=dict(
            xaxis_title="Distance x (m)",
            yaxis_title="Visualization thickness (m)",
            zaxis_title="Elevation z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.8, y=-1.8, z=1.1),
                center=dict(x=0.0, y=0.0, z=0.0),
            ),
        ),
        legend=dict(
            x=0.02,
            y=0.98,
        ),
    )

    html_path = Path("pohang_soil_mesh_3d.html").resolve()

    fig.write_html(str(html_path), auto_open=True)

    print(f"Saved: {html_path}")
    webbrowser.open(html_path.as_uri())
    
if __name__ == "__main__":
    main()
    