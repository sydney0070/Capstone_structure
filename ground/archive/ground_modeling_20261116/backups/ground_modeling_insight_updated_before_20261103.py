import os
import tempfile
from pathlib import Path
import argparse
import copy
import importlib.util
import math
import time

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "matplotlib"),
)

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import webbrowser
from insight_motion import run_uniform_excitation, inspect_mseed_headers

try:
    import openseespy.opensees as ops
except ImportError:
    ops = None

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR


# ============================================================
# 0. 실행 선택
# ============================================================

# "Earth", "Mars", "Moon" 중 하나
TARGET_BODY = "Earth"

# 한 번 실행할 때 하나만 생성:
# "surface" 또는 "solid"
VISUALIZATION_MODE = "solid"

# 선택한 천체의 DEM과 지반 재료를 이용해 OpenSees 3D 지반모델을 생성한다.
BUILD_OPENSEES_MODEL = True

# 3D 모델을 만든 뒤 중력 자중 해석까지 수행한다.
RUN_GRAVITY_ANALYSIS = True

# HTML을 생성한 뒤 브라우저에서 자동으로 열기
AUTO_OPEN_HTML = True

# 전역 DEM은 지형 확인과 후보지 선택에 사용하고, OpenSees 해석은 기본적으로
# 작은 local SSI domain에서 수행한다.
USE_LOCAL_OPENSEES_DOMAIN = True
LOCAL_DOMAIN_WIDTH_M = 150.0
LOCAL_DOMAIN_HEIGHT_M = 150.0
LOCAL_N_X = 13
LOCAL_N_Y = 13

# 구조물 크기가 정해지면 foundation footprint 검토에 사용한다.
FOUNDATION_FOOTPRINT_WIDTH_M = None
FOUNDATION_FOOTPRINT_LENGTH_M = None

BOUNDARY_CONDITION_MODE = "gravity_static_box"

DEFAULT_VERTICAL_MESH = {
    "mode": "layer_aware",
    "near_surface_depth": 3.0,
    "near_surface_max_dz": 1.5,
    "max_dz": 3.0,
}

VALIDATION_LIMITS = {
    "min_volume_m3": 1.0e-9,
    "min_jacobian_det": 1.0e-12,
    "aspect_ratio_warning": 12.0,
    "aspect_ratio_fail": 40.0,
    "force_balance_relative_warning": 1.0e-4,
    "force_balance_relative_fail": 1.0e-2,
    "max_gravity_displacement_warning_m": 0.50,
}


# ============================================================
# 1-1. 지구 설정
# ============================================================

# Earth_file : NASADEM HGT 파일 경로
HGT_PATH = SCRIPT_DIR / "n36e129.hgt"

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

# 깊이 방향 요소 개수
N_DEPTH = 12

# OpenSees 해석 모델은 3D solid brick 기반으로 생성한다.
OPENSEES_MODEL_DIMENSION = "3D"
BRICK_ELEMENT_TYPE = "bbarBrick"
STATIC_ANALYSIS_SYSTEM = "UmfPack"

# ============================================================
# 1-2. 달 설정
# ============================================================

# Moon_file : IMG, XML 파일 경로
IMG_PATH = SCRIPT_DIR / "ldem_1024_90s_75s_000_030.img"
XML_PATH = SCRIPT_DIR / "ldem_1024_90s_75s_000_030.xml"


# 달 PDS4 DEM은 매우 크므로 전체 원본을 그대로 Plotly에 넣지 않는다.
# 먼저 원본 raster에서 일부 영역을 읽고, 지구와 동일한 N_X × N_Y 격자로 변환한다.
MOON_READ_MODE = "center_crop"       # "center_crop" 또는 "full_downsample"
MOON_CROP_WIDTH_PIXELS = 3000
MOON_CROP_HEIGHT_PIXELS = 3000

# None이면 파일 중앙을 사용한다.
MOON_CROP_COL_OFFSET = None
MOON_CROP_ROW_OFFSET = None

# ============================================================
# 1-3. 화성 설정
# ============================================================

# 화성 DEM 파일 경로는 사용자가 자료를 받은 뒤 직접 입력한다.
# 예시: MARS_DEM_PATH = Path("mars_dem.tif")
MARS_LABEL_PATH = SCRIPT_DIR / "megt00n000hb.lbl"
MARS_IMAGE_PATH = SCRIPT_DIR / "megt00n000hb.img"

# 화성 DEM이 매우 큰 경우 일부 영역만 읽거나 전체를 다운샘플링한다.
MARS_READ_MODE = "center_crop"       # "center_crop" 또는 "full_downsample"
MARS_CROP_WIDTH_PIXELS = 3000
MARS_CROP_HEIGHT_PIXELS = 3000
MARS_CROP_COL_OFFSET = None
MARS_CROP_ROW_OFFSET = None


# ============================================================
# 2-1. 지반 재료 물성 및 재료 모델 설정
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
        "nd": 3,
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
        "nd": 3,
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
# 2-2. 달 지반 재료 물성 및 재료 모델 설정
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
        "nd": 3,
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
        "nd": 3,
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
# 2-3. 화성 지반 재료 물성 및 재료 모델 설정
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
        "nd": 3,
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
        "nd": 3,
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


MATERIAL_CONFIDENCE = {
    "Earth": {
        "global": "reasonable_initial_assumption",
        "PDMY02_reference": "OpenSees_suggested_sand_parameters_require_site_calibration",
    },
    "Moon": {
        "global": "reasonable_initial_assumption",
        "rho": "requires_calibration",
        "refShearModul": "requires_calibration",
        "refBulkModul": "requires_calibration",
        "frictionAng": "requires_calibration",
        "peakShearStra": "requires_calibration",
        "refPress": "model_reference_pressure_not_atmospheric_pressure",
        "pa": "critical_state_reference_pressure_not_atmospheric_pressure",
        "c": "requires_calibration",
        "contraction_dilation": "requires_calibration",
        "layer_depths": "assumed_until_site_data_available",
    },
    "Mars": {
        "global": "reasonable_initial_assumption",
        "rho": "requires_calibration",
        "refShearModul": "requires_calibration",
        "refBulkModul": "requires_calibration",
        "frictionAng": "requires_calibration",
        "peakShearStra": "requires_calibration",
        "refPress": "model_reference_pressure_not_atmospheric_pressure",
        "pa": "critical_state_reference_pressure_not_atmospheric_pressure",
        "c": "requires_calibration",
        "contraction_dilation": "requires_calibration",
        "layer_depths": "assumed_until_site_data_available",
    },
}


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
    "n_depth": N_DEPTH,
    "layer_boundaries": [1.5, 5.0, 10.0, 16.0, 24.0],
    "vertical_mesh": copy.deepcopy(DEFAULT_VERTICAL_MESH),
    "boundary_condition_mode": BOUNDARY_CONDITION_MODE,
    "material_confidence": MATERIAL_CONFIDENCE["Earth"],
}

MOON_CONFIG = {
    "body_name": "Moon",
    "gravity": G_MOON,
    "materials": MOON_SOIL_MATERIALS,
    "material_order": MOON_SOIL_ORDER,
    "assign_material": assign_moon_material,
    "base_depth": BASE_DEPTH,
    "n_depth": N_DEPTH,
    "layer_boundaries": [2.0, 8.0, 15.0],
    "vertical_mesh": copy.deepcopy(DEFAULT_VERTICAL_MESH),
    "boundary_condition_mode": BOUNDARY_CONDITION_MODE,
    "material_confidence": MATERIAL_CONFIDENCE["Moon"],
}


MARS_CONFIG = {
    "body_name": "Mars",
    "gravity": G_MARS,
    "materials": MARS_SOIL_MATERIALS,
    "material_order": MARS_SOIL_ORDER,
    "assign_material": assign_mars_material,
    "base_depth": BASE_DEPTH,
    "n_depth": N_DEPTH,
    "layer_boundaries": [3.0, 12.0, 22.0],
    "vertical_mesh": copy.deepcopy(DEFAULT_VERTICAL_MESH),
    "boundary_condition_mode": BOUNDARY_CONDITION_MODE,
    "material_confidence": MATERIAL_CONFIDENCE["Mars"],
}


def require_opensees():
    if ops is None:
        raise ImportError(
            "OpenSeesPy가 현재 Python 환경에 설치되어 있지 않습니다. "
            "이 파일의 DEM 처리, geometry 검증, synthetic terrain 생성 함수는 "
            "사용할 수 있지만 OpenSees 모델 생성/해석에는 openseespy가 필요합니다."
        )
    return ops


def register_opensees_soil_materials(materials, material_order, analysis_nd=None):
    """
    OpenSees nDMaterial 등록.

    입력된 재료 딕셔너리와 순서에 따라 OpenSees nDMaterial을 등록한다.

    - ElasticIsotropic 재료는 탄성체로 등록
    - PressureDependMultiYield02 재료는 비선형 재료 태그 목록에 추가
    - analysis_nd가 주어지면 PDMY02의 nd를 해석 차원에 맞춰 덮어쓴다.
    """

    require_opensees()
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
            nd = int(analysis_nd if analysis_nd is not None else m["nd"])

            ops.nDMaterial(
                "PressureDependMultiYield02",
                mat_tag,
                nd,
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

    require_opensees()
    for mat_tag in material_tags:
        ops.updateMaterialStage(
            "-material",
            int(mat_tag),
            "-stage",
            int(stage),
        )


def configure_static_analysis(system_name, test_tol, max_iter, algorithm_name, d_lambda):
    require_opensees()
    ops.constraints("Transformation")
    ops.numberer("RCM")

    try:
        ops.system(system_name)
    except Exception as exc:
        print(
            f"[Warning] OpenSees system '{system_name}' is unavailable: {exc}. "
            "Falling back to BandGeneral."
        )
        ops.system("BandGeneral")

    ops.test("NormDispIncr", float(test_tol), int(max_iter))
    ops.algorithm(algorithm_name)
    ops.integrator("LoadControl", float(d_lambda))
    ops.analysis("Static")


def run_elastic_gravity_then_plastic_stage(
    nodal_self_weight_loads,
    nonlinear_material_tags,
    ndm=3,
    n_steps=10,
    d_lambda=0.1,
):
    """
    탄성 중력 재하 단계 → loadConst → 소성 전환 단계.

    1) PDMY/PDMY02 재료를 stage 0으로 설정
       - 전 소성 거동 비활성화
       - 순수 선형 탄성체로 자중 응력 상태 형성

    2) pattern Plain으로 지반 자중 하중 적용
       - 2D: {node_tag: Fy}
       - 3D: {node_tag: (Fx, Fy, Fz)}

    3) Static analysis 수행

    4) loadConst('-time', 0.0)으로 중력 재하 결과 고정

    5) PDMY/PDMY02 재료를 stage 1로 전환
       - 이후 외력 재하, 말뚝 압입, 지진 동적 재하 가능
    """

    require_opensees()
    ndm = int(ndm)
    total_applied_load = np.zeros(ndm, dtype=float)

    print("\n[Stage 0] Set PDMY/PDMY02 materials to elastic stage.")
    set_pdmy_material_stage(nonlinear_material_tags, stage=0)

    # 중력 하중 패턴
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)

    for node_tag, load_value in nodal_self_weight_loads.items():
        if np.isscalar(load_value):
            if int(ndm) == 2:
                load_vector = (0.0, float(load_value))
            elif int(ndm) == 3:
                load_vector = (0.0, 0.0, float(load_value))
            else:
                raise ValueError("ndm은 2 또는 3이어야 합니다.")
        else:
            load_vector = tuple(float(value) for value in load_value)
            if len(load_vector) != int(ndm):
                raise ValueError(
                    f"절점 {node_tag}의 하중 벡터 길이({len(load_vector)})가 "
                    f"ndm={ndm}와 맞지 않습니다."
                )

        ops.load(int(node_tag), *load_vector)
        total_applied_load += np.asarray(load_vector, dtype=float)

    # 정적 탄성 중력 해석 설정
    configure_static_analysis(
        system_name=STATIC_ANALYSIS_SYSTEM,
        test_tol=1.0e-8,
        max_iter=50,
        algorithm_name="Newton",
        d_lambda=d_lambda,
    )

    print("[Stage 0] Running elastic gravity analysis...")
    ok = ops.analyze(int(n_steps))

    if ok != 0:
        print("[Warning] Gravity analysis failed with Newton. Trying ModifiedNewton...")

        ops.wipeAnalysis()
        configure_static_analysis(
            system_name=STATIC_ANALYSIS_SYSTEM,
            test_tol=1.0e-7,
            max_iter=100,
            algorithm_name="ModifiedNewton",
            d_lambda=d_lambda,
        )

        ok = ops.analyze(int(n_steps))

    if ok != 0:
        raise RuntimeError("탄성 중력 재하 해석이 수렴하지 않았음.")

    print("[Stage 0] Gravity analysis converged.")

    ops.reactions()
    reaction_vector = np.zeros(ndm, dtype=float)
    for node_tag in ops.getNodeTags():
        reaction = np.asarray(ops.nodeReaction(int(node_tag)), dtype=float)
        if reaction.size >= ndm:
            reaction_vector += reaction[:ndm]

    displacement_norms = []
    vertical_displacements = []
    for node_tag in ops.getNodeTags():
        displacement = np.asarray(ops.nodeDisp(int(node_tag)), dtype=float)
        if displacement.size >= ndm:
            displacement_norms.append(float(np.linalg.norm(displacement[:ndm])))
            vertical_displacements.append(
                float(displacement[2] if ndm == 3 else displacement[1])
            )

    max_displacement = (
        float(np.nanmax(displacement_norms))
        if displacement_norms
        else 0.0
    )
    mean_displacement = (
        float(np.nanmean(displacement_norms))
        if displacement_norms
        else 0.0
    )
    max_vertical_displacement = (
        float(np.nanmax(np.abs(vertical_displacements)))
        if vertical_displacements
        else 0.0
    )

    balance_vector = reaction_vector + total_applied_load
    applied_norm = float(np.linalg.norm(total_applied_load))
    balance_norm = float(np.linalg.norm(balance_vector))
    relative_balance_error = (
        balance_norm / applied_norm
        if applied_norm > 0.0
        else 0.0
    )

    gravity_summary = {
        "analysis_ok": int(ok),
        "analysis_steps": int(n_steps),
        "load_step": float(d_lambda),
        "total_applied_load": tuple(float(v) for v in total_applied_load),
        "reaction_vector": tuple(float(v) for v in reaction_vector),
        "balance_vector": tuple(float(v) for v in balance_vector),
        "balance_norm": balance_norm,
        "relative_balance_error": relative_balance_error,
        "max_displacement_m": max_displacement,
        "mean_displacement_m": mean_displacement,
        "max_vertical_displacement_m": max_vertical_displacement,
        "stage_after_gravity": 0,
    }

    print("[Stage 0] Gravity equilibrium summary")
    print("  Total applied load:", gravity_summary["total_applied_load"])
    print("  Total reaction:", gravity_summary["reaction_vector"])
    print("  Relative balance error:", relative_balance_error)
    print("  Max displacement:", max_displacement, "m")

    # 중력 응력장과 변위 상태 고정
    ops.loadConst("-time", 0.0)
    print("[Stage 0] Gravity state fixed by loadConst('-time', 0.0).")

    # 이후 해석을 위해 analysis object 초기화
    ops.wipeAnalysis()

    # 소성 상태로 전환
    print("[Stage 1] Switch PDMY/PDMY02 materials to plastic stage.")
    set_pdmy_material_stage(nonlinear_material_tags, stage=1)
    gravity_summary["stage_after_gravity"] = 1

    print("[Stage 1] Plastic material stage is ready.")
    print("         이후 말뚝 압입, 인장 재하, 지진 동적 가속도 재하를 추가하면 됨.\n")

    return gravity_summary

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

def extract_area_patch(
    dem,
    center_lat=AREA_CENTER_LAT,
    center_lon=AREA_CENTER_LON,
    width_m=AREA_WIDTH_M,
    height_m=AREA_HEIGHT_M,
    n_x=N_X,
    n_y=N_Y,
):
    """
    중심 위도/경도를 기준으로 일정 크기의 직사각형 DEM 영역을 추출한다.
    결과는 Plotly 3D surface/solid 시각화에 사용한다.
    """

    earth_radius = 6371000.0
    lat_rad = np.deg2rad(center_lat)

    # m → degree 변환
    dlat = (float(height_m) / 2.0) / earth_radius * (180.0 / np.pi)
    dlon = (float(width_m) / 2.0) / (earth_radius * np.cos(lat_rad)) * (180.0 / np.pi)

    lat_min = center_lat - dlat
    lat_max = center_lat + dlat
    lon_min = center_lon - dlon
    lon_max = center_lon + dlon

    lats = np.linspace(lat_min, lat_max, int(n_y))
    lons = np.linspace(lon_min, lon_max, int(n_x))

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

    rasterio = None
    Resampling = None
    Window = None
    transform_xy = None
    try:
        import rasterio
        from rasterio.enums import Resampling
        from rasterio.windows import Window
        from rasterio.transform import xy as transform_xy
    except ImportError as exc:
        print(
            "[Warning] rasterio/GDAL이 없어 화성 PDS3 DEM은 "
            "NumPy fallback reader로 읽습니다."
        )
        print("          import error:", exc)

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

            return {
                "col_off": col_off,
                "row_off": row_off,
                "width": crop_width,
                "height": crop_height,
            }

        if MARS_READ_MODE == "full_downsample":
            return {
                "col_off": 0,
                "row_off": 0,
                "width": int(width),
                "height": int(height),
            }

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
        if rasterio is None:
            raise RuntimeError("rasterio/GDAL unavailable")

        with rasterio.open(label_path) as src:
            if src.count < 1:
                raise RuntimeError(
                    "화성 PDS3 라벨에서 raster band를 찾지 못했습니다."
                )

            window = make_window(src.width, src.height)
            raster_window = Window(
                col_off=window["col_off"],
                row_off=window["row_off"],
                width=window["width"],
                height=window["height"],
            )

            raw = src.read(
                1,
                window=raster_window,
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

            base_transform = src.window_transform(raster_window)
            display_transform = base_transform * rasterio.Affine.scale(
                float(window["width"]) / N_X,
                float(window["height"]) / N_Y,
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
                    int(window["col_off"]),
                    int(window["row_off"]),
                    int(window["width"]),
                    int(window["height"]),
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

    row_start = int(window["row_off"])
    row_stop = row_start + int(window["height"])
    col_start = int(window["col_off"])
    col_stop = col_start + int(window["width"])

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
            int(window["col_off"]),
            int(window["row_off"]),
            int(window["width"]),
            int(window["height"]),
        ),
    )

    return X, Y, Z_rel, Z_abs, rows, cols


def polygon_area_xy(points):
    """x-y 평면에 투영한 다각형 면적을 계산한다."""
    coords = np.asarray(points, dtype=float)
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]
    return 0.5 * abs(
        np.dot(x_coords, np.roll(y_coords, -1))
        - np.dot(y_coords, np.roll(x_coords, -1))
    )


def signed_polygon_area_xy(points):
    """x-y 평면 투영 다각형의 부호 있는 면적을 계산한다."""
    coords = np.asarray(points, dtype=float)
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]
    return 0.5 * (
        np.dot(x_coords, np.roll(y_coords, -1))
        - np.dot(y_coords, np.roll(x_coords, -1))
    )


def add_nodal_load_vector(loads, node_tag, vector):
    previous = loads.get(int(node_tag), (0.0, 0.0, 0.0))
    loads[int(node_tag)] = tuple(
        previous[i] + float(vector[i])
        for i in range(3)
    )


def ensure_finite_terrain_grid(X, Y, Z):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)

    if X.shape != Y.shape or X.shape != Z.shape:
        raise ValueError(
            "X, Y, Z는 같은 2D shape이어야 합니다. "
            f"현재 shape: X={X.shape}, Y={Y.shape}, Z={Z.shape}"
        )

    if X.ndim != 2:
        raise ValueError("terrain grid는 2D 배열이어야 합니다.")

    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(Y)):
        raise ValueError("terrain X/Y 좌표에 NaN 또는 inf가 있습니다.")

    if not np.any(np.isfinite(Z)):
        raise ValueError("terrain Z에 유효한 고도값이 없습니다.")

    if not np.all(np.isfinite(Z)):
        Z = Z.copy()
        Z[~np.isfinite(Z)] = np.nanmean(Z)

    return X, Y, Z


def orient_terrain_grid_for_opensees(X, Y, Z):
    """
    OpenSees brick ordering에 맞도록 local x/y 격자 방향을 정리한다.

    bbarBrick node 1-4는 하부면, node 5-8은 상부면이며, 각 면의
    1-2-3-4 방향이 양의 z 방향 법선을 만들도록 x-y cell을 둔다.
    """
    X, Y, Z = ensure_finite_terrain_grid(X, Y, Z)
    operations = []

    mean_dx = float(np.nanmean(np.diff(X, axis=1))) if X.shape[1] > 1 else 0.0
    mean_dy = float(np.nanmean(np.diff(Y, axis=0))) if Y.shape[0] > 1 else 0.0

    if mean_dx < 0.0:
        X = np.flip(X, axis=1)
        Y = np.flip(Y, axis=1)
        Z = np.flip(Z, axis=1)
        operations.append("flipped_x_columns")

    if mean_dy < 0.0:
        X = np.flip(X, axis=0)
        Y = np.flip(Y, axis=0)
        Z = np.flip(Z, axis=0)
        operations.append("flipped_y_rows")

    if X.shape[0] > 1 and X.shape[1] > 1:
        first_cell = [
            (X[0, 0], Y[0, 0]),
            (X[0, 1], Y[0, 1]),
            (X[1, 1], Y[1, 1]),
            (X[1, 0], Y[1, 0]),
        ]
        if signed_polygon_area_xy(first_cell) < 0.0:
            X = np.flip(X, axis=0)
            Y = np.flip(Y, axis=0)
            Z = np.flip(Z, axis=0)
            operations.append("flipped_y_rows_for_positive_xy_area")

    X = X - float(X[0, 0])
    Y = Y - float(Y[0, 0])

    metadata = {
        "operations": operations,
        "mean_dx_m": float(np.nanmean(np.diff(X, axis=1))) if X.shape[1] > 1 else 0.0,
        "mean_dy_m": float(np.nanmean(np.diff(Y, axis=0))) if X.shape[0] > 1 else 0.0,
    }
    return X, Y, Z, metadata


def resample_rectilinear_grid(X, Y, Z, target_x, target_y):
    """
    거의 rectilinear인 DEM local grid를 target x/y 축으로 bilinear resampling한다.
    """
    X, Y, Z = ensure_finite_terrain_grid(X, Y, Z)
    x_axis = np.asarray(np.nanmean(X, axis=0), dtype=float)
    y_axis = np.asarray(np.nanmean(Y, axis=1), dtype=float)

    x_order = np.argsort(x_axis)
    y_order = np.argsort(y_axis)
    x_axis = x_axis[x_order]
    y_axis = y_axis[y_order]
    Z_sorted = Z[np.ix_(y_order, x_order)]

    if np.any(np.diff(x_axis) <= 0.0) or np.any(np.diff(y_axis) <= 0.0):
        raise ValueError(
            "local site 추출은 단조 증가하는 rectilinear X/Y grid가 필요합니다."
        )

    target_x = np.asarray(target_x, dtype=float)
    target_y = np.asarray(target_y, dtype=float)

    x_pass = np.empty((len(y_axis), len(target_x)), dtype=float)
    for iy in range(len(y_axis)):
        x_pass[iy, :] = np.interp(target_x, x_axis, Z_sorted[iy, :])

    z_out = np.empty((len(target_y), len(target_x)), dtype=float)
    for ix in range(len(target_x)):
        z_out[:, ix] = np.interp(target_y, y_axis, x_pass[:, ix])

    X_out, Y_out = np.meshgrid(target_x, target_y)
    return X_out, Y_out, z_out


def extract_local_site(
    X,
    Y,
    Z,
    center_x=None,
    center_y=None,
    width_m=LOCAL_DOMAIN_WIDTH_M,
    height_m=LOCAL_DOMAIN_HEIGHT_M,
    target_nx=LOCAL_N_X,
    target_ny=LOCAL_N_Y,
    structure_footprint=None,
):
    """
    Global DEM grid에서 OpenSees SSI용 local terrain patch를 추출한다.

    아직 foundation/contact를 만들지는 않고, 나중에 footprint 아래 surface node를
    찾을 수 있도록 중심과 footprint metadata만 함께 반환한다.
    """
    X, Y, Z, orientation = orient_terrain_grid_for_opensees(X, Y, Z)

    x_min = float(np.nanmin(X))
    x_max = float(np.nanmax(X))
    y_min = float(np.nanmin(Y))
    y_max = float(np.nanmax(Y))

    if center_x is None:
        center_x = 0.5 * (x_min + x_max)
    if center_y is None:
        center_y = 0.5 * (y_min + y_max)

    half_width = 0.5 * float(width_m)
    half_height = 0.5 * float(height_m)

    local_x_min = max(float(center_x) - half_width, x_min)
    local_x_max = min(float(center_x) + half_width, x_max)
    local_y_min = max(float(center_y) - half_height, y_min)
    local_y_max = min(float(center_y) + half_height, y_max)

    if local_x_max <= local_x_min or local_y_max <= local_y_min:
        raise ValueError("요청한 local site 범위가 global terrain 밖에 있습니다.")

    target_x = np.linspace(local_x_min, local_x_max, int(target_nx))
    target_y = np.linspace(local_y_min, local_y_max, int(target_ny))
    X_local, Y_local, Z_local = resample_rectilinear_grid(
        X,
        Y,
        Z,
        target_x,
        target_y,
    )

    Z_local = Z_local - float(np.nanmin(Z_local))

    source_dx = np.abs(np.diff(np.nanmean(X, axis=0)))
    source_dy = np.abs(np.diff(np.nanmean(Y, axis=1)))
    source_dx = source_dx[np.isfinite(source_dx) & (source_dx > 0.0)]
    source_dy = source_dy[np.isfinite(source_dy) & (source_dy > 0.0)]
    target_dx = (
        float(local_x_max - local_x_min) / max(int(target_nx) - 1, 1)
    )
    target_dy = (
        float(local_y_max - local_y_min) / max(int(target_ny) - 1, 1)
    )
    source_dx_mean = float(np.mean(source_dx)) if source_dx.size else math.nan
    source_dy_mean = float(np.mean(source_dy)) if source_dy.size else math.nan

    metadata = {
        "source_shape": tuple(int(v) for v in Z.shape),
        "local_shape": tuple(int(v) for v in Z_local.shape),
        "center_x_m": float(center_x),
        "center_y_m": float(center_y),
        "x_range_m": (float(local_x_min), float(local_x_max)),
        "y_range_m": (float(local_y_min), float(local_y_max)),
        "width_m": float(local_x_max - local_x_min),
        "height_m": float(local_y_max - local_y_min),
        "target_nx": int(target_nx),
        "target_ny": int(target_ny),
        "target_dx_m": target_dx,
        "target_dy_m": target_dy,
        "source_dx_mean_m": source_dx_mean,
        "source_dy_mean_m": source_dy_mean,
        "oversampled_from_source_grid": bool(
            (
                np.isfinite(source_dx_mean)
                and target_dx < 0.75 * source_dx_mean
            )
            or (
                np.isfinite(source_dy_mean)
                and target_dy < 0.75 * source_dy_mean
            )
        ),
        "orientation": orientation,
        "structure_footprint": structure_footprint,
    }

    return X_local, Y_local, Z_local, metadata


def make_vertical_surface_depths(config):
    """
    지표면 기준 depth level을 만든다.

    마지막 바닥면은 평면 `base_z`로 따로 생성하므로 반환 리스트에는
    `base_depth` 자체를 포함하지 않는다.
    """
    base_depth = float(config["base_depth"])
    vertical_mesh = config.get("vertical_mesh", {})
    mode = vertical_mesh.get("mode", "uniform")

    if base_depth <= 0.0:
        raise ValueError("base_depth는 양수이어야 합니다.")

    if mode == "uniform":
        n_depth = int(config.get("n_depth", N_DEPTH))
        if n_depth < 1:
            raise ValueError("n_depth는 1 이상이어야 합니다.")
        depths = np.linspace(0.0, base_depth, n_depth + 1)[:-1]
        return [float(depth) for depth in depths]

    if mode != "layer_aware":
        raise ValueError(f"지원하지 않는 vertical mesh mode: {mode}")

    layer_boundaries = [
        float(depth)
        for depth in config.get("layer_boundaries", [])
        if 0.0 < float(depth) < base_depth
    ]
    breakpoints = sorted({0.0, base_depth, *layer_boundaries})

    max_dz = float(vertical_mesh.get("max_dz", 3.0))
    near_surface_depth = float(vertical_mesh.get("near_surface_depth", 3.0))
    near_surface_max_dz = float(vertical_mesh.get("near_surface_max_dz", max_dz))

    if max_dz <= 0.0 or near_surface_max_dz <= 0.0:
        raise ValueError("vertical mesh dz는 양수이어야 합니다.")

    depths = []
    for start, end in zip(breakpoints[:-1], breakpoints[1:]):
        interval = end - start
        dz_limit = near_surface_max_dz if start < near_surface_depth else max_dz
        n_seg = max(1, int(math.ceil(interval / dz_limit)))
        for i in range(n_seg):
            depth = start + interval * i / n_seg
            if not depths or abs(depth - depths[-1]) > 1.0e-9:
                depths.append(float(depth))

    return [
        float(depth)
        for depth in depths
        if depth < base_depth - 1.0e-9
    ]


BRICK_LOCAL_COORDS = np.array(
    [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ],
    dtype=float,
)

BRICK_EDGE_PAIRS = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
)


def brick_shape_derivatives(xi, eta, zeta):
    derivatives = np.empty((8, 3), dtype=float)
    for i, (xi_i, eta_i, zeta_i) in enumerate(BRICK_LOCAL_COORDS):
        derivatives[i, 0] = 0.125 * xi_i * (1.0 + eta * eta_i) * (1.0 + zeta * zeta_i)
        derivatives[i, 1] = 0.125 * eta_i * (1.0 + xi * xi_i) * (1.0 + zeta * zeta_i)
        derivatives[i, 2] = 0.125 * zeta_i * (1.0 + xi * xi_i) * (1.0 + eta * eta_i)
    return derivatives


def brick_jacobian_determinants(coords):
    coords = np.asarray(coords, dtype=float)
    if coords.shape != (8, 3):
        raise ValueError(f"brick coords shape는 (8, 3)이어야 합니다: {coords.shape}")

    gauss = (-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0))
    dets = []
    for xi in gauss:
        for eta in gauss:
            for zeta in gauss:
                derivatives = brick_shape_derivatives(xi, eta, zeta)
                jacobian = derivatives.T @ coords
                dets.append(float(np.linalg.det(jacobian)))
    return np.asarray(dets, dtype=float)


def brick_geometry_metrics(coords):
    coords = np.asarray(coords, dtype=float)
    dets = brick_jacobian_determinants(coords)
    edge_lengths = np.array(
        [
            np.linalg.norm(coords[j] - coords[i])
            for i, j in BRICK_EDGE_PAIRS
        ],
        dtype=float,
    )

    min_edge = float(np.nanmin(edge_lengths))
    max_edge = float(np.nanmax(edge_lengths))
    aspect_ratio = (
        max_edge / min_edge
        if min_edge > 0.0
        else math.inf
    )

    return {
        "signed_volume_m3": float(np.sum(dets)),
        "absolute_gauss_volume_m3": float(np.sum(np.abs(dets))),
        "min_jacobian_det": float(np.nanmin(dets)),
        "max_jacobian_det": float(np.nanmax(dets)),
        "jacobian_dets": tuple(float(v) for v in dets),
        "min_edge_m": min_edge,
        "max_edge_m": max_edge,
        "aspect_ratio": float(aspect_ratio),
    }


def build_boundary_node_groups(node_tags, nx, ny, nz):
    surface = [
        node_tags[(ix, iy, 0)]
        for iy in range(ny)
        for ix in range(nx)
    ]
    bottom = [
        node_tags[(ix, iy, nz)]
        for iy in range(ny)
        for ix in range(nx)
    ]
    x_min = [
        node_tags[(0, iy, iz)]
        for iy in range(ny)
        for iz in range(nz + 1)
    ]
    x_max = [
        node_tags[(nx - 1, iy, iz)]
        for iy in range(ny)
        for iz in range(nz + 1)
    ]
    y_min = [
        node_tags[(ix, 0, iz)]
        for ix in range(nx)
        for iz in range(nz + 1)
    ]
    y_max = [
        node_tags[(ix, ny - 1, iz)]
        for ix in range(nx)
        for iz in range(nz + 1)
    ]

    side = sorted(set(x_min + x_max + y_min + y_max))
    return {
        "surface": surface,
        "bottom": bottom,
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
        "side": side,
    }


def apply_boundary_conditions(node_coords, node_groups, mode):
    """
    현재 구현된 경계조건은 정적 중력/기초하중용 단순 box boundary이다.
    동적 SSI에는 absorbing/free-field/dashpot boundary가 별도로 필요하다.
    """
    require_opensees()
    fixities = {
        int(tag): [0, 0, 0]
        for tag in node_coords
    }

    if mode not in {"gravity_static_box", "fixed_base_lateral_rollers", "experimental_fixed_base_free_sides"}:
        raise ValueError(f"지원하지 않는 boundary_condition_mode: {mode}")

    for tag in node_groups["bottom"]:
        fixities[int(tag)] = [1, 1, 1]

    if mode != "experimental_fixed_base_free_sides":
        for tag in node_groups["x_min"] + node_groups["x_max"]:
            fixities[int(tag)][0] = 1

        for tag in node_groups["y_min"] + node_groups["y_max"]:
            fixities[int(tag)][1] = 1

    for tag, fixity in fixities.items():
        if any(fixity):
            ops.fix(int(tag), *fixity)

    return fixities


def select_surface_nodes_within_footprint(
    model_data,
    center_x,
    center_y,
    width_m,
    length_m,
):
    """향후 foundation 연결을 위한 surface node selector."""
    half_width = 0.5 * float(width_m)
    half_length = 0.5 * float(length_m)
    selected = []

    for tag in model_data.get("surface_nodes", []):
        x, y, z = model_data["node_coords"][tag]
        if (
            abs(x - center_x) <= half_width
            and abs(y - center_y) <= half_length
        ):
            selected.append(
                {
                    "tag": int(tag),
                    "coord": (float(x), float(y), float(z)),
                }
            )

    return selected


def build_opensees_soil_model_3d(X, Y, ground_z, config, run_gravity=True):
    """
    DEM 지표면을 상부 경계로 하는 3D OpenSees solid 지반 모델 생성.

    생성 방식:
    - OpenSees 도메인: -ndm 3, -ndf 3
    - 요소: 8절점 bbarBrick/stdBrick 계열 solid 요소
    - 수평 격자: DEM X, Y, Z 좌표를 그대로 사용
    - 깊이 격자: 지표면에서 base_z까지 선형 보간
    - 자중: 각 brick 체적에서 계산한 등가 절점하중으로 적용
    """

    require_opensees()
    generation_start = time.perf_counter()

    X, Y, ground_z, terrain_orientation = orient_terrain_grid_for_opensees(
        X,
        Y,
        ground_z,
    )

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    body_name = config["body_name"]
    gravity = float(config["gravity"])
    materials = config["materials"]
    material_order = config["material_order"]
    assign_material = config["assign_material"]
    base_depth = float(config["base_depth"])
    vertical_surface_depths = make_vertical_surface_depths(config)

    if not vertical_surface_depths or abs(vertical_surface_depths[0]) > 1.0e-9:
        vertical_surface_depths = [0.0, *vertical_surface_depths]

    nonlinear_material_tags = register_opensees_soil_materials(
        materials,
        material_order,
        analysis_nd=3,
    )

    ny, nx = ground_z.shape
    if nx < 2 or ny < 2:
        raise ValueError("3D 지반 모델에는 최소 2 x 2 이상의 terrain grid가 필요합니다.")

    base_z = float(np.nanmin(ground_z) - base_depth)
    nz = len(vertical_surface_depths)
    nominal_vertical_depths = [
        *[float(depth) for depth in vertical_surface_depths],
        base_depth,
    ]

    node_tags = {}
    node_coords = {}
    node_depths_from_surface = {}
    node_tag = 1

    # --------------------------------------------------------
    # 1) 3D 노드 생성
    # --------------------------------------------------------
    for iy in range(ny):
        for ix in range(nx):
            for iz in range(nz + 1):
                if iz < nz:
                    depth = float(vertical_surface_depths[iz])
                    z = float(ground_z[iy, ix] - depth)
                else:
                    z = base_z
                    depth = float(ground_z[iy, ix] - base_z)

                x = float(X[iy, ix])
                y = float(Y[iy, ix])

                ops.node(node_tag, x, y, z)
                node_tags[(ix, iy, iz)] = node_tag
                node_coords[node_tag] = (x, y, z)
                node_depths_from_surface[node_tag] = depth
                node_tag += 1

    # --------------------------------------------------------
    # 2) 3D brick 요소 생성 + 요소별 재료 배정
    # --------------------------------------------------------
    ele_tag = 1
    ele_materials = {}
    nodal_self_weight_loads = {}

    x_min = float(np.nanmin(X))
    x_max = float(np.nanmax(X))
    y_min = float(np.nanmin(Y))
    y_max = float(np.nanmax(Y))
    z_min = float(np.nanmin(ground_z))
    z_max = float(np.nanmax(ground_z))

    for iy in range(ny - 1):
        for ix in range(nx - 1):
            top_cell_xy = [
                (X[iy, ix], Y[iy, ix]),
                (X[iy, ix + 1], Y[iy, ix + 1]),
                (X[iy + 1, ix + 1], Y[iy + 1, ix + 1]),
                (X[iy + 1, ix], Y[iy + 1, ix]),
            ]
            cell_area_xy = polygon_area_xy(top_cell_xy)

            ground_corners = np.array(
                [
                    ground_z[iy, ix],
                    ground_z[iy, ix + 1],
                    ground_z[iy + 1, ix + 1],
                    ground_z[iy + 1, ix],
                ],
                dtype=float,
            )
            ground_center = float(np.nanmean(ground_corners))
            x_center = float(
                np.nanmean(
                    [
                        X[iy, ix],
                        X[iy, ix + 1],
                        X[iy + 1, ix + 1],
                        X[iy + 1, ix],
                    ]
                )
            )
            y_center = float(
                np.nanmean(
                    [
                        Y[iy, ix],
                        Y[iy, ix + 1],
                        Y[iy + 1, ix + 1],
                        Y[iy + 1, ix],
                    ]
                )
            )

            x_norm = (x_center - x_min) / (x_max - x_min + 1.0e-12)
            y_norm = (y_center - y_min) / (y_max - y_min + 1.0e-12)
            elevation_norm = (ground_center - z_min) / (z_max - z_min + 1.0e-12)

            for iz in range(nz):
                # bbarBrick/stdBrick 표준 8절점 순서:
                # 하부면 1-2-3-4, 상부면 5-6-7-8.
                # OpenSees 공식 brick figure의 local zeta=-1/+1에 맞춘 순서이다.
                element_nodes = (
                    node_tags[(ix, iy, iz + 1)],
                    node_tags[(ix + 1, iy, iz + 1)],
                    node_tags[(ix + 1, iy + 1, iz + 1)],
                    node_tags[(ix, iy + 1, iz + 1)],
                    node_tags[(ix, iy, iz)],
                    node_tags[(ix + 1, iy, iz)],
                    node_tags[(ix + 1, iy + 1, iz)],
                    node_tags[(ix, iy + 1, iz)],
                )

                coords = np.array(
                    [node_coords[tag] for tag in element_nodes],
                    dtype=float,
                )
                geometry = brick_geometry_metrics(coords)

                if (
                    geometry["signed_volume_m3"] <= VALIDATION_LIMITS["min_volume_m3"]
                    or geometry["min_jacobian_det"] <= VALIDATION_LIMITS["min_jacobian_det"]
                ):
                    raise ValueError(
                        "유효하지 않은 bbarBrick geometry가 감지되었습니다. "
                        f"grid=(ix={ix}, iy={iy}, iz={iz}), "
                        f"volume={geometry['signed_volume_m3']}, "
                        f"min_detJ={geometry['min_jacobian_det']}"
                    )

                z_center = float(np.nanmean(coords[:, 2]))
                depth_from_surface = ground_center - z_center

                mat_name = assign_material(
                    depth_from_surface=depth_from_surface,
                    x_norm=x_norm,
                    y_norm=y_norm,
                    elevation_norm=elevation_norm,
                )

                if mat_name not in materials:
                    raise ValueError(f"정의되지 않은 지반 재료명: {mat_name}")

                mat = materials[mat_name]
                mat_tag = int(mat["tag"])
                rho = float(mat["rho"])

                # 자중은 아래 등가 절점하중으로 명시 적용하므로
                # brick 요소의 body force 항은 0으로 둔다.
                ops.element(
                    BRICK_ELEMENT_TYPE,
                    ele_tag,
                    *element_nodes,
                    mat_tag,
                    0.0,
                    0.0,
                    0.0,
                )

                top_z = coords[4:, 2]
                bottom_z = coords[:4, 2]
                element_height = max(float(np.nanmean(top_z - bottom_z)), 0.0)
                projected_volume = float(cell_area_xy * element_height)
                volume = float(geometry["signed_volume_m3"])

                element_weight = rho * gravity * volume
                nodal_weight = -element_weight / 8.0

                for nd in element_nodes:
                    add_nodal_load_vector(
                        nodal_self_weight_loads,
                        nd,
                        (0.0, 0.0, nodal_weight),
                    )

                ele_materials[ele_tag] = {
                    "element_type": BRICK_ELEMENT_TYPE,
                    "nodes": element_nodes,
                    "grid_index": (ix, iy, iz),
                    "material": mat_name,
                    "material_tag": mat_tag,
                    "model": mat["model"],
                    "rho": rho,
                    "depth_m": float(depth_from_surface),
                    "x_center_m": x_center,
                    "y_center_m": y_center,
                    "z_center_m": float(z_center),
                    "volume_m3": volume,
                    "projected_volume_m3": projected_volume,
                    "volume_method": "2x2x2_gauss_integrated_hexahedron",
                    "volume_error_fraction_vs_projected": (
                        (volume - projected_volume) / volume
                        if abs(volume) > 0.0
                        else math.inf
                    ),
                    "geometry": {
                        key: value
                        for key, value in geometry.items()
                        if key != "jacobian_dets"
                    },
                }

                ele_tag += 1

    # --------------------------------------------------------
    # 3) 3D 경계조건
    # --------------------------------------------------------
    boundary_mode = config.get("boundary_condition_mode", BOUNDARY_CONDITION_MODE)
    node_groups = build_boundary_node_groups(node_tags, nx, ny, nz)
    fixities = apply_boundary_conditions(node_coords, node_groups, boundary_mode)

    used_material_tags = {
        int(element["material_tag"])
        for element in ele_materials.values()
    }
    used_nonlinear_material_tags = [
        tag
        for tag in nonlinear_material_tags
        if int(tag) in used_material_tags
    ]

    gravity_summary = None
    if run_gravity:
        gravity_summary = run_elastic_gravity_then_plastic_stage(
            nodal_self_weight_loads=nodal_self_weight_loads,
            nonlinear_material_tags=used_nonlinear_material_tags,
            ndm=3,
            n_steps=10,
            d_lambda=0.1,
        )

    generation_time = time.perf_counter() - generation_start
    total_mass = float(
        sum(
            element["rho"] * element["volume_m3"]
            for element in ele_materials.values()
        )
    )
    expected_weight = float(total_mass * gravity)
    total_nodal_gravity_load = np.zeros(3, dtype=float)
    for load_vector in nodal_self_weight_loads.values():
        total_nodal_gravity_load += np.asarray(load_vector, dtype=float)

    model_data = {
        "dimension": 3,
        "ndm": 3,
        "ndf": 3,
        "body_name": body_name,
        "element_type": BRICK_ELEMENT_TYPE,
        "node_tags": node_tags,
        "node_coords": node_coords,
        "node_depths_from_surface": node_depths_from_surface,
        "elements": ele_materials,
        "nodal_self_weight_loads": nodal_self_weight_loads,
        "nonlinear_material_tags": nonlinear_material_tags,
        "used_nonlinear_material_tags": used_nonlinear_material_tags,
        "materials": materials,
        "material_order": material_order,
        "material_confidence": config.get("material_confidence", {}),
        "foundation_footprint": config.get("foundation_footprint"),
        "local_site": config.get("local_site"),
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "base_z": base_z,
        "base_depth": base_depth,
        "layer_boundaries": tuple(
            float(depth)
            for depth in config.get("layer_boundaries", [])
        ),
        "gravity": gravity,
        "vertical_surface_depths_m": tuple(float(v) for v in vertical_surface_depths),
        "nominal_vertical_depths_m": tuple(float(v) for v in nominal_vertical_depths),
        "terrain_orientation": terrain_orientation,
        "boundary_condition_mode": boundary_mode,
        "boundary_node_groups": node_groups,
        "surface_nodes": node_groups["surface"],
        "bottom_nodes": node_groups["bottom"],
        "side_boundary_nodes": node_groups["side"],
        "fixities": fixities,
        "total_soil_mass_kg": total_mass,
        "expected_total_weight_n": expected_weight,
        "total_nodal_gravity_load": tuple(float(v) for v in total_nodal_gravity_load),
        "gravity_summary": gravity_summary,
        "analysis_stage": (
            "gravity_initialized_plastic_stage_ready"
            if run_gravity
            else "model_generated_without_gravity"
        ),
        "generation_time_sec": float(generation_time),
    }

    validation_report = validate_soil_model(model_data)
    model_data["validation_report"] = validation_report
    print_validation_report(validation_report)

    print(f"{body_name} OpenSees 3D soil model generated.")
    print("OpenSees domain: ndm=3, ndf=3")
    print("Element type:", BRICK_ELEMENT_TYPE)
    print("Number of nodes:", len(ops.getNodeTags()))
    print("Number of elements:", len(ops.getEleTags()))
    print("Nonlinear PDMY/PDMY02 material tags:", nonlinear_material_tags)
    print("Total soil mass:", total_mass, "kg")
    print("Expected total weight:", expected_weight, "N")
    print("Model generation time:", generation_time, "s")

    return model_data


def validation_item(status, category, check, message, value=None):
    return {
        "status": status,
        "category": category,
        "check": check,
        "message": message,
        "value": value,
    }


def summarize_status(items):
    if any(item["status"] == "FAIL" for item in items):
        return "FAIL"
    if any(item["status"] == "WARNING" for item in items):
        return "WARNING"
    return "PASS"


def stats_tuple(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "min": math.nan,
            "mean": math.nan,
            "max": math.nan,
        }
    return {
        "min": float(np.min(values)),
        "mean": float(np.mean(values)),
        "max": float(np.max(values)),
    }


def compute_mesh_adequacy(model_data):
    node_tags = model_data["node_tags"]
    node_coords = model_data["node_coords"]
    nx = int(model_data["nx"])
    ny = int(model_data["ny"])
    nz = int(model_data["nz"])

    dx_values = []
    dy_values = []
    dz_values = []

    for iy in range(ny):
        for ix in range(nx - 1):
            n1 = node_tags[(ix, iy, 0)]
            n2 = node_tags[(ix + 1, iy, 0)]
            dx_values.append(
                np.linalg.norm(
                    np.asarray(node_coords[n2]) - np.asarray(node_coords[n1])
                )
            )

    for iy in range(ny - 1):
        for ix in range(nx):
            n1 = node_tags[(ix, iy, 0)]
            n2 = node_tags[(ix, iy + 1, 0)]
            dy_values.append(
                np.linalg.norm(
                    np.asarray(node_coords[n2]) - np.asarray(node_coords[n1])
                )
            )

    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                n1 = node_tags[(ix, iy, iz)]
                n2 = node_tags[(ix, iy, iz + 1)]
                dz_values.append(
                    abs(float(node_coords[n1][2]) - float(node_coords[n2][2]))
                )

    aspect_values = [
        element["geometry"]["aspect_ratio"]
        for element in model_data["elements"].values()
    ]
    volume_values = [
        element["volume_m3"]
        for element in model_data["elements"].values()
    ]
    det_values = [
        element["geometry"]["min_jacobian_det"]
        for element in model_data["elements"].values()
    ]

    x_values = [coord[0] for coord in node_coords.values()]
    y_values = [coord[1] for coord in node_coords.values()]
    z_values = [coord[2] for coord in node_coords.values()]

    return {
        "domain_x_length_m": float(np.nanmax(x_values) - np.nanmin(x_values)),
        "domain_y_length_m": float(np.nanmax(y_values) - np.nanmin(y_values)),
        "domain_depth_m": float(np.nanmax(z_values) - np.nanmin(z_values)),
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "dx_m": stats_tuple(dx_values),
        "dy_m": stats_tuple(dy_values),
        "dz_m": stats_tuple(dz_values),
        "element_aspect_ratio": stats_tuple(aspect_values),
        "element_volume_m3": stats_tuple(volume_values),
        "min_jacobian_det": stats_tuple(det_values),
        "node_count": len(node_coords),
        "element_count": len(model_data["elements"]),
    }


def validate_soil_model(model_data):
    items = []
    node_coords = model_data["node_coords"]
    elements = model_data["elements"]
    materials = model_data["materials"]
    material_tags = {
        int(material["tag"])
        for material in materials.values()
    }

    coord_array = np.asarray(list(node_coords.values()), dtype=float)
    if np.all(np.isfinite(coord_array)):
        items.append(
            validation_item(
                "PASS",
                "Geometry",
                "finite_node_coordinates",
                "All node coordinates are finite.",
            )
        )
    else:
        items.append(
            validation_item(
                "FAIL",
                "Geometry",
                "finite_node_coordinates",
                "NaN or inf node coordinates were found.",
            )
        )

    min_volume = min(
        element["volume_m3"]
        for element in elements.values()
    )
    if min_volume <= VALIDATION_LIMITS["min_volume_m3"]:
        items.append(
            validation_item(
                "FAIL",
                "Geometry",
                "positive_element_volume",
                "At least one element has zero or negative volume.",
                min_volume,
            )
        )
    else:
        items.append(
            validation_item(
                "PASS",
                "Geometry",
                "positive_element_volume",
                "All elements have positive Gauss-integrated volume.",
                min_volume,
            )
        )

    min_det = min(
        element["geometry"]["min_jacobian_det"]
        for element in elements.values()
    )
    if min_det <= VALIDATION_LIMITS["min_jacobian_det"]:
        items.append(
            validation_item(
                "FAIL",
                "Geometry",
                "positive_jacobian",
                "At least one Gauss point has non-positive Jacobian determinant.",
                min_det,
            )
        )
    else:
        items.append(
            validation_item(
                "PASS",
                "Geometry",
                "positive_jacobian",
                "All Gauss point Jacobian determinants are positive.",
                min_det,
            )
        )

    max_aspect = max(
        element["geometry"]["aspect_ratio"]
        for element in elements.values()
    )
    if max_aspect > VALIDATION_LIMITS["aspect_ratio_fail"]:
        aspect_status = "FAIL"
        aspect_message = "Element aspect ratio is extreme."
    elif max_aspect > VALIDATION_LIMITS["aspect_ratio_warning"]:
        aspect_status = "WARNING"
        aspect_message = "Element aspect ratio is high; refine horizontal or vertical mesh if SSI stresses are sensitive."
    else:
        aspect_status = "PASS"
        aspect_message = "Element aspect ratios are within the configured preliminary limit."
    items.append(
        validation_item(
            aspect_status,
            "Geometry",
            "element_aspect_ratio",
            aspect_message,
            max_aspect,
        )
    )

    max_projected_error = max(
        abs(element["volume_error_fraction_vs_projected"])
        for element in elements.values()
    )
    projected_status = "PASS" if max_projected_error < 1.0e-6 else "WARNING"
    items.append(
        validation_item(
            projected_status,
            "Geometry",
            "gauss_volume_vs_projected_volume",
            "Gauss volume is compared with the previous projected-area approximation.",
            max_projected_error,
        )
    )

    all_nodes = set(node_coords)
    referenced_nodes = []
    invalid_refs = []
    element_keys = set()
    duplicate_elements = []

    for ele_tag, element in elements.items():
        nodes = tuple(int(tag) for tag in element["nodes"])
        referenced_nodes.extend(nodes)
        invalid_refs.extend(tag for tag in nodes if tag not in all_nodes)
        key = tuple(sorted(nodes))
        if key in element_keys:
            duplicate_elements.append(ele_tag)
        element_keys.add(key)

    orphan_nodes = sorted(all_nodes - set(referenced_nodes))
    if invalid_refs:
        items.append(
            validation_item(
                "FAIL",
                "Connectivity",
                "valid_node_references",
                "Elements reference nodes that do not exist.",
                invalid_refs[:10],
            )
        )
    else:
        items.append(
            validation_item(
                "PASS",
                "Connectivity",
                "valid_node_references",
                "All element node references are valid.",
            )
        )

    orphan_status = "WARNING" if orphan_nodes else "PASS"
    items.append(
        validation_item(
            orphan_status,
            "Connectivity",
            "orphan_nodes",
            "No orphan nodes were found." if not orphan_nodes else "Nodes not used by any element were found.",
            orphan_nodes[:10],
        )
    )

    duplicate_status = "FAIL" if duplicate_elements else "PASS"
    items.append(
        validation_item(
            duplicate_status,
            "Connectivity",
            "duplicate_elements",
            "No duplicate element connectivity was found." if not duplicate_elements else "Duplicate element connectivity was found.",
            duplicate_elements[:10],
        )
    )

    rounded_coords = [
        tuple(round(float(value), 9) for value in coord)
        for coord in node_coords.values()
    ]
    duplicate_coord_count = len(rounded_coords) - len(set(rounded_coords))
    coord_status = "WARNING" if duplicate_coord_count else "PASS"
    items.append(
        validation_item(
            coord_status,
            "Connectivity",
            "duplicate_node_coordinates",
            "No duplicate node coordinates were found." if not duplicate_coord_count else "Duplicate node coordinates were found.",
            duplicate_coord_count,
        )
    )

    undefined_material_tags = sorted(
        {
            int(element["material_tag"])
            for element in elements.values()
            if int(element["material_tag"]) not in material_tags
        }
    )
    material_status = "FAIL" if undefined_material_tags else "PASS"
    items.append(
        validation_item(
            material_status,
            "Materials",
            "defined_material_tags",
            "All element material tags are defined." if not undefined_material_tags else "Some element material tags are not defined.",
            undefined_material_tags,
        )
    )

    bad_materials = []
    for name, material in materials.items():
        rho = float(material.get("rho", math.nan))
        if not np.isfinite(rho) or rho <= 0.0:
            bad_materials.append((name, "rho", rho))

        if material["model"] == "ElasticIsotropic":
            E = float(material.get("E", math.nan))
            nu = float(material.get("nu", math.nan))
            if not np.isfinite(E) or E <= 0.0:
                bad_materials.append((name, "E", E))
            if not np.isfinite(nu) or not (-1.0 < nu < 0.5):
                bad_materials.append((name, "nu", nu))

        elif material["model"] == "PressureDependMultiYield02":
            nd = int(material.get("nd", -1))
            if nd != int(model_data["ndm"]):
                bad_materials.append((name, "nd", nd))
            for key in ("refShearModul", "refBulkModul", "refPress", "pa"):
                value = float(material.get(key, math.nan))
                if not np.isfinite(value) or value <= 0.0:
                    bad_materials.append((name, key, value))
            friction = float(material.get("frictionAng", math.nan))
            if not np.isfinite(friction) or not (0.0 < friction < 60.0):
                bad_materials.append((name, "frictionAng", friction))
            for key in ("contrac1", "contrac2", "contrac3", "dilat1", "dilat2", "dilat3"):
                value = float(material.get(key, math.nan))
                if not np.isfinite(value) or value < 0.0:
                    bad_materials.append((name, key, value))

    items.append(
        validation_item(
            "FAIL" if bad_materials else "PASS",
            "Materials",
            "material_parameter_ranges",
            "Material parameters are inside basic numerical ranges." if not bad_materials else "Invalid material parameter values were found.",
            bad_materials[:15],
        )
    )

    total_mass = float(model_data["total_soil_mass_kg"])
    expected_weight = float(model_data["expected_total_weight_n"])
    applied_load = np.asarray(model_data["total_nodal_gravity_load"], dtype=float)
    expected_load_z = -expected_weight
    load_error = abs(float(applied_load[2]) - expected_load_z)
    load_rel_error = load_error / max(abs(expected_load_z), 1.0)

    if load_rel_error > VALIDATION_LIMITS["force_balance_relative_fail"]:
        load_status = "FAIL"
    elif load_rel_error > VALIDATION_LIMITS["force_balance_relative_warning"]:
        load_status = "WARNING"
    else:
        load_status = "PASS"

    items.append(
        validation_item(
            load_status,
            "Loads",
            "self_weight_load_sum",
            "Applied nodal self-weight matches mass times gravity.",
            {
                "total_mass_kg": total_mass,
                "expected_weight_n": expected_weight,
                "applied_load_z_n": float(applied_load[2]),
                "relative_error": load_rel_error,
            },
        )
    )

    gravity_summary = model_data.get("gravity_summary")
    if gravity_summary is None:
        items.append(
            validation_item(
                "WARNING",
                "Loads",
                "reaction_equilibrium",
                "Gravity analysis was not run; reaction equilibrium is not verified.",
            )
        )
    else:
        reaction_error = float(gravity_summary["relative_balance_error"])
        if reaction_error > VALIDATION_LIMITS["force_balance_relative_fail"]:
            reaction_status = "FAIL"
        elif reaction_error > VALIDATION_LIMITS["force_balance_relative_warning"]:
            reaction_status = "WARNING"
        else:
            reaction_status = "PASS"
        items.append(
            validation_item(
                reaction_status,
                "Loads",
                "reaction_equilibrium",
                "Sum of reactions balances applied gravity load.",
                {
                    "relative_balance_error": reaction_error,
                    "reaction_vector": gravity_summary["reaction_vector"],
                    "applied_load": gravity_summary["total_applied_load"],
                },
            )
        )

        max_disp = float(gravity_summary["max_displacement_m"])
        disp_status = (
            "WARNING"
            if max_disp > VALIDATION_LIMITS["max_gravity_displacement_warning_m"]
            else "PASS"
        )
        items.append(
            validation_item(
                disp_status,
                "Loads",
                "gravity_displacement_size",
                "Gravity displacement magnitude is within preliminary sanity limit." if disp_status == "PASS" else "Gravity displacement is large; check stiffness, mesh, and boundary assumptions.",
                max_disp,
            )
        )

    fixities = model_data["fixities"]
    bottom_unfixed = [
        tag
        for tag in model_data["bottom_nodes"]
        if fixities.get(int(tag)) != [1, 1, 1]
    ]
    items.append(
        validation_item(
            "FAIL" if bottom_unfixed else "PASS",
            "Boundary",
            "bottom_fixity",
            "All bottom nodes are fixed in x, y, z." if not bottom_unfixed else "Some bottom nodes are not fully fixed.",
            bottom_unfixed[:10],
        )
    )

    side_nodes = model_data["side_boundary_nodes"]
    items.append(
        validation_item(
            "PASS" if side_nodes else "FAIL",
            "Boundary",
            "side_boundary_groups",
            "Side boundary node groups are available." if side_nodes else "Side boundary node groups are missing.",
            len(side_nodes),
        )
    )

    if model_data["boundary_condition_mode"] == "gravity_static_box":
        items.append(
            validation_item(
                "WARNING",
                "Boundary",
                "dynamic_ssi_boundary",
                "Current fixed-base/lateral-roller boundary is for gravity and preliminary static loading, not dynamic SSI.",
            )
        )

    mesh = compute_mesh_adequacy(model_data)
    items.append(
        validation_item(
            "PASS",
            "Mesh",
            "mesh_counts",
            "Mesh counts were computed.",
            {
                "nodes": mesh["node_count"],
                "elements": mesh["element_count"],
                "nx": mesh["nx"],
                "ny": mesh["ny"],
                "nz": mesh["nz"],
            },
        )
    )

    vertical_depths = list(model_data.get("vertical_surface_depths_m", []))
    base_depth = float(model_data.get("base_depth", math.inf))
    layer_boundaries = [
        float(depth)
        for depth in model_data.get("layer_boundaries", [])
        if 0.0 < float(depth) < base_depth
    ]
    missing_boundaries = [
        depth
        for depth in layer_boundaries
        if not any(abs(depth - mesh_depth) < 1.0e-8 for mesh_depth in vertical_depths)
    ]
    items.append(
        validation_item(
            "WARNING" if missing_boundaries else "PASS",
            "Mesh",
            "layer_boundary_alignment",
            "Layer boundaries are represented by vertical mesh levels." if not missing_boundaries else "Some layer boundaries do not coincide with vertical mesh levels.",
            missing_boundaries,
        )
    )

    local_site = model_data.get("local_site")
    if local_site is not None:
        if local_site.get("oversampled_from_source_grid", False):
            local_status = "WARNING"
            local_message = (
                "Local FEM grid is finer than the source DEM grid; terrain detail is interpolated, not measured."
            )
        else:
            local_status = "PASS"
            local_message = "Local FEM spacing is compatible with the source DEM grid spacing."

        items.append(
            validation_item(
                local_status,
                "Mesh",
                "local_dem_resolution",
                local_message,
                {
                    "source_dx_mean_m": local_site.get("source_dx_mean_m"),
                    "source_dy_mean_m": local_site.get("source_dy_mean_m"),
                    "target_dx_m": local_site.get("target_dx_m"),
                    "target_dy_m": local_site.get("target_dy_m"),
                },
            )
        )

    footprint = model_data.get("foundation_footprint")
    if footprint is not None:
        dx_mean = float(mesh["dx_m"]["mean"])
        dy_mean = float(mesh["dy_m"]["mean"])
        width_ratio = float(footprint["width_m"]) / max(dx_mean, 1.0e-12)
        length_ratio = float(footprint["length_m"]) / max(dy_mean, 1.0e-12)
        min_ratio = min(width_ratio, length_ratio)

        if min_ratio < 2.0:
            footprint_status = "FAIL"
            footprint_message = "Foundation footprint is smaller than about two soil elements across."
        elif min_ratio < 4.0:
            footprint_status = "WARNING"
            footprint_message = "Foundation footprint has fewer than four soil elements across; SSI mesh may be too coarse."
        else:
            footprint_status = "PASS"
            footprint_message = "Foundation footprint is resolved by the current surface mesh."

        items.append(
            validation_item(
                footprint_status,
                "Mesh",
                "foundation_width_over_element_size",
                footprint_message,
                {
                    "width_over_dx": width_ratio,
                    "length_over_dy": length_ratio,
                },
            )
        )

    report = {
        "overall_status": summarize_status(items),
        "items": items,
        "mesh": mesh,
        "limits": copy.deepcopy(VALIDATION_LIMITS),
    }
    return report


def print_validation_report(report):
    items = report["items"]
    counts = {
        status: sum(1 for item in items if item["status"] == status)
        for status in ("PASS", "WARNING", "FAIL")
    }

    print("\n[Validation] Soil model sanity check")
    print("  Overall:", report["overall_status"])
    print(
        "  Counts:",
        f"PASS={counts['PASS']}",
        f"WARNING={counts['WARNING']}",
        f"FAIL={counts['FAIL']}",
    )

    for item in items:
        if item["status"] in {"WARNING", "FAIL"}:
            print(
                f"  [{item['status']}] "
                f"{item['category']} / {item['check']}: "
                f"{item['message']}"
            )
            if item["value"] is not None:
                print("    value:", item["value"])


def save_opensees_3d_model_html(model_data):
    """
    실제 OpenSees 3D 모델의 노드 좌표와 brick 요소 연결을 그대로 사용해
    외곽 형상, 요소 경계, 재료 분포를 HTML로 저장한다.
    """

    node_coords = model_data["node_coords"]
    elements = model_data["elements"]
    materials = model_data["materials"]
    material_order = model_data["material_order"]
    body_name = model_data["body_name"]

    node_id_to_plot_index = {
        tag: idx
        for idx, tag in enumerate(sorted(node_coords))
    }
    x_flat = [
        node_coords[tag][0]
        for tag in sorted(node_coords)
    ]
    y_flat = [
        node_coords[tag][1]
        for tag in sorted(node_coords)
    ]
    z_flat = [
        node_coords[tag][2]
        for tag in sorted(node_coords)
    ]

    face_definitions = (
        ("bottom", (0, 3, 2, 1)),
        ("top", (4, 5, 6, 7)),
        ("y_min", (0, 1, 5, 4)),
        ("x_max", (1, 2, 6, 5)),
        ("y_max", (2, 3, 7, 6)),
        ("x_min", (3, 0, 4, 7)),
    )

    face_registry = {}
    for ele_tag, element in elements.items():
        nodes = element["nodes"]
        for face_name, face_indices in face_definitions:
            face_nodes = tuple(nodes[index] for index in face_indices)
            face_key = tuple(sorted(face_nodes))

            if face_key in face_registry:
                face_registry[face_key]["count"] += 1
            else:
                face_registry[face_key] = {
                    "count": 1,
                    "nodes": face_nodes,
                    "element_tag": ele_tag,
                    "face_name": face_name,
                    "material": element["material"],
                }

    exterior_faces = [
        face
        for face in face_registry.values()
        if face["count"] == 1
    ]

    I = []
    J = []
    K = []
    face_colors = []

    fallback_colors = [
        "#4e79a7",
        "#f28e2b",
        "#e15759",
        "#76b7b2",
        "#59a14f",
        "#edc948",
        "#b07aa1",
        "#ff9da7",
    ]
    fallback_by_material = {
        name: fallback_colors[index % len(fallback_colors)]
        for index, name in enumerate(material_order)
    }

    def material_color(material_name):
        material = materials.get(material_name, {})
        return material.get(
            "color",
            fallback_by_material.get(material_name, "#999999"),
        )

    for face in exterior_faces:
        plot_nodes = [
            node_id_to_plot_index[node]
            for node in face["nodes"]
        ]
        color = material_color(face["material"])

        I.extend([plot_nodes[0], plot_nodes[0]])
        J.extend([plot_nodes[1], plot_nodes[2]])
        K.extend([plot_nodes[2], plot_nodes[3]])
        face_colors.extend([color, color])

    edge_keys = set()
    edge_x = []
    edge_y = []
    edge_z = []

    def add_edge(node_a, node_b):
        edge_key = tuple(sorted((node_a, node_b)))
        if edge_key in edge_keys:
            return
        edge_keys.add(edge_key)
        xa, ya, za = node_coords[node_a]
        xb, yb, zb = node_coords[node_b]
        edge_x.extend([xa, xb, None])
        edge_y.extend([ya, yb, None])
        edge_z.extend([za, zb, None])

    for face in exterior_faces:
        nodes = face["nodes"]
        add_edge(nodes[0], nodes[1])
        add_edge(nodes[1], nodes[2])
        add_edge(nodes[2], nodes[3])
        add_edge(nodes[3], nodes[0])

    fig = go.Figure()

    fig.add_trace(
        go.Mesh3d(
            x=x_flat,
            y=y_flat,
            z=z_flat,
            i=I,
            j=J,
            k=K,
            facecolor=face_colors,
            opacity=0.88,
            flatshading=True,
            name="OpenSees exterior brick faces",
            hovertemplate=(
                "OpenSees 3D soil mesh<br>"
                "x=%{x:.2f} m<br>"
                "y=%{y:.2f} m<br>"
                "z=%{z:.2f} m"
                "<extra></extra>"
            ),
            showscale=False,
        )
    )

    fig.add_trace(
        go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            line=dict(color="black", width=1.5),
            name="Exterior element edges",
            hoverinfo="skip",
        )
    )

    used_materials = {
        element["material"]
        for element in elements.values()
    }
    for material_name in material_order:
        if material_name not in used_materials:
            continue

        material = materials[material_name]
        fig.add_trace(
            go.Scatter3d(
                x=[None],
                y=[None],
                z=[None],
                mode="markers",
                marker=dict(size=8, color=material_color(material_name)),
                name=material.get("label", material_name),
                visible="legendonly",
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        title=(
            f"{body_name} OpenSees 3D Soil Model "
            f"(ndm={model_data['ndm']}, ndf={model_data['ndf']}, "
            f"grid={model_data['nx']} x {model_data['ny']}, "
            f"depth elements={model_data['nz']}, "
            f"nodes={len(node_coords)}, "
            f"elements={len(elements)} {model_data['element_type']})"
        ),
        template="plotly_white",
        width=1250,
        height=900,
        scene=dict(
            xaxis_title="x distance (m)",
            yaxis_title="y distance (m)",
            zaxis_title="Elevation / depth z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.45, y=-1.65, z=1.05),
            ),
        ),
        legend=dict(
            x=0.02,
            y=0.98,
        ),
    )

    html_path = (OUTPUT_DIR / f"{body_name.lower()}_opensees_3d_soil_model.html").resolve()
    fig.write_html(str(html_path), auto_open=False)

    print(f"Saved exact OpenSees model HTML: {html_path}")
    return html_path


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

    html_path = (OUTPUT_DIR / f"{body_name.lower()}_area_3d_terrain.html").resolve()
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

    html_path = (OUTPUT_DIR / f"{body_name.lower()}_area_3d_solid.html").resolve()
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

    html_path = (OUTPUT_DIR / "pohang_area_3d_soil_material_distribution.html").resolve()
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


def make_synthetic_terrain(
    kind="flat",
    nx=5,
    ny=5,
    width_m=40.0,
    height_m=40.0,
):
    x = np.linspace(0.0, float(width_m), int(nx))
    y = np.linspace(0.0, float(height_m), int(ny))
    X, Y = np.meshgrid(x, y)

    kind = kind.lower()
    if kind == "flat":
        Z = np.zeros_like(X)
    elif kind == "slope":
        Z = 0.08 * X + 0.04 * Y
    elif kind == "irregular":
        Z = (
            0.35 * np.sin(2.0 * np.pi * X / max(float(width_m), 1.0))
            + 0.25 * np.cos(2.0 * np.pi * Y / max(float(height_m), 1.0))
            + 0.10 * np.sin(2.0 * np.pi * (X + Y) / max(float(width_m + height_m), 1.0))
        )
        Z = Z - float(np.nanmin(Z))
    else:
        raise ValueError("synthetic terrain kind는 flat, slope, irregular 중 하나여야 합니다.")

    return X, Y, Z


def run_synthetic_smoke_tests(kind="all", run_gravity=True):
    require_opensees()

    if kind == "all":
        terrain_kinds = ["flat", "slope", "irregular"]
    else:
        terrain_kinds = [kind]

    body_configs = [
        ("Earth", EARTH_CONFIG),
        ("Moon", MOON_CONFIG),
        ("Mars", MARS_CONFIG),
    ]

    results = []
    for body_name, base_config in body_configs:
        for terrain_kind in terrain_kinds:
            print(
                "\n[Synthetic Smoke]",
                f"body={body_name}",
                f"terrain={terrain_kind}",
            )
            X, Y, Z = make_synthetic_terrain(
                terrain_kind,
                nx=5,
                ny=5,
                width_m=40.0,
                height_m=40.0,
            )
            config = copy.deepcopy(base_config)
            config["base_depth"] = 18.0
            config["vertical_mesh"] = copy.deepcopy(DEFAULT_VERTICAL_MESH)

            model_data = build_opensees_soil_model_3d(
                X,
                Y,
                Z,
                config=config,
                run_gravity=run_gravity,
            )
            results.append(
                {
                    "body": body_name,
                    "terrain": terrain_kind,
                    "overall_status": model_data["validation_report"]["overall_status"],
                    "node_count": len(model_data["node_coords"]),
                    "element_count": len(model_data["elements"]),
                    "total_mass_kg": model_data["total_soil_mass_kg"],
                    "gravity_balance_error": (
                        model_data["gravity_summary"]["relative_balance_error"]
                        if model_data["gravity_summary"] is not None
                        else None
                    ),
                    "max_gravity_displacement_m": (
                        model_data["gravity_summary"]["max_displacement_m"]
                        if model_data["gravity_summary"] is not None
                        else None
                    ),
                }
            )

    print("\n[Synthetic Smoke] Summary")
    for result in results:
        print(
            " ",
            result["body"],
            result["terrain"],
            result["overall_status"],
            "nodes=",
            result["node_count"],
            "elements=",
            result["element_count"],
            "eq_error=",
            result["gravity_balance_error"],
        )

    return results


def collect_required_file_status():
    required = {
        "Earth NASADEM HGT": HGT_PATH,
        "Moon PDS4 XML": XML_PATH,
        "Moon PDS4 IMG": IMG_PATH,
        "Mars PDS3 LBL": MARS_LABEL_PATH,
        "Mars PDS3 IMG": MARS_IMAGE_PATH,
    }

    present = []
    missing = []
    for label, path in required.items():
        entry = {
            "label": label,
            "path": str(Path(path).resolve()),
        }
        if Path(path).exists():
            present.append(entry)
        else:
            missing.append(entry)

    dependencies = {
        "openseespy": ops is not None,
        "rasterio": importlib.util.find_spec("rasterio") is not None,
    }

    missing_but_non_blocking = []
    missing_required_for_final_verification = []

    if not dependencies["rasterio"]:
        missing_but_non_blocking.append(
            {
                "label": "rasterio/GDAL",
                "reason": "Mars PDS3 fallback and synthetic tests can still run.",
            }
        )
        missing_required_for_final_verification.append(
            {
                "label": "rasterio/GDAL",
                "reason": "Required to verify Moon PDS4 XML/IMG reader and the rasterio Mars reader path.",
            }
        )

    if not dependencies["openseespy"]:
        missing_required_for_final_verification.append(
            {
                "label": "openseespy",
                "reason": "Required to build the OpenSees domain and verify gravity convergence/reactions.",
            }
        )

    return {
        "present": present,
        "missing": missing,
        "missing_but_non_blocking": missing_but_non_blocking,
        "missing_and_required_for_final_verification": missing_required_for_final_verification,
        "dependencies": dependencies,
    }


def estimate_hgt_native_spacing_m(dem_shape, latitude_deg):
    n_rows, n_cols = dem_shape
    earth_radius = 6371000.0
    lat_spacing_deg = 1.0 / max(int(n_rows) - 1, 1)
    lon_spacing_deg = 1.0 / max(int(n_cols) - 1, 1)
    dx = (
        earth_radius
        * math.cos(math.radians(float(latitude_deg)))
        * math.radians(lon_spacing_deg)
    )
    dy = earth_radius * math.radians(lat_spacing_deg)
    return float(dx), float(dy)


def make_direct_local_site_metadata(
    source_shape,
    local_shape,
    center_x,
    center_y,
    width_m,
    height_m,
    target_nx,
    target_ny,
    source_dx_mean_m,
    source_dy_mean_m,
    structure_footprint=None,
):
    target_dx = float(width_m) / max(int(target_nx) - 1, 1)
    target_dy = float(height_m) / max(int(target_ny) - 1, 1)
    return {
        "source_shape": tuple(int(v) for v in source_shape),
        "local_shape": tuple(int(v) for v in local_shape),
        "center_x_m": float(center_x),
        "center_y_m": float(center_y),
        "x_range_m": (0.0, float(width_m)),
        "y_range_m": (0.0, float(height_m)),
        "width_m": float(width_m),
        "height_m": float(height_m),
        "target_nx": int(target_nx),
        "target_ny": int(target_ny),
        "target_dx_m": target_dx,
        "target_dy_m": target_dy,
        "source_dx_mean_m": float(source_dx_mean_m),
        "source_dy_mean_m": float(source_dy_mean_m),
        "oversampled_from_source_grid": bool(
            target_dx < 0.75 * float(source_dx_mean_m)
            or target_dy < 0.75 * float(source_dy_mean_m)
        ),
        "orientation": {"operations": ["direct_source_sampling"]},
        "structure_footprint": structure_footprint,
    }


def print_required_file_status(report):
    print("\n[Required files]")
    print("  Present")
    for item in report["present"]:
        print("   -", item["label"], ":", item["path"])

    print("  Missing")
    if report["missing"]:
        for item in report["missing"]:
            print("   -", item["label"], ":", item["path"])
    else:
        print("   - None")

    print("  Missing but non-blocking")
    if report["missing_but_non_blocking"]:
        for item in report["missing_but_non_blocking"]:
            print("   -", item["label"], ":", item["reason"])
    else:
        print("   - None")

    print("  Missing and required for final verification")
    if report["missing_and_required_for_final_verification"]:
        for item in report["missing_and_required_for_final_verification"]:
            print("   -", item["label"], ":", item["reason"])
    else:
        print("   - None")


def get_body_config(body_name):
    body_name = body_name.strip().lower()
    if body_name == "earth":
        return "Earth", copy.deepcopy(EARTH_CONFIG)
    if body_name == "moon":
        return "Moon", copy.deepcopy(MOON_CONFIG)
    if body_name == "mars":
        return "Mars", copy.deepcopy(MARS_CONFIG)
    raise ValueError("TARGET_BODY는 'Earth', 'Mars', 'Moon' 중 하나여야 합니다.")


def parse_cli_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Earth/Moon/Mars OpenSees 3D soil model generator."
    )
    parser.add_argument(
        "--body",
        choices=["Earth", "Moon", "Mars", "earth", "moon", "mars"],
        default=TARGET_BODY,
    )
    parser.add_argument(
        "--synthetic-test",
        choices=["none", "flat", "slope", "irregular", "all"],
        default="none",
    )
    parser.add_argument("--no-gravity", action="store_true")
    parser.add_argument("--inspect-insight-mseed", nargs=3, type=Path, metavar="MSEED", default=None,
                        help="세 관측 파일 헤더 검사 후 종료; ObsPy 없어도 됨")
    parser.add_argument("--dynamic-csv", type=Path, default=None,
                        help="StationXML로 기기보정하고 ENU 변환한 가속도 CSV (m/s²)")
    parser.add_argument("--experimental-surface-as-base", action="store_true",
                        help="SEIS 지표면 관측운동을 기반암 입력으로 가정 (검증용 전용)")
    parser.add_argument("--experimental-free-sides", action="store_true",
                        help="반사 경계인 측면 자유 경계를 사용 (연구결과 해석 금지)")
    parser.add_argument("--axes-enu-verified", action="store_true",
                        help="DEM과 FEM 축이 각각 동/북/상임을 검증한 경우만 지정")
    parser.add_argument("--dynamic-scale", type=float, default=1.0)
    parser.add_argument("--dynamic-max-steps", type=int, default=None)
    parser.add_argument("--dynamic-animation-stride", type=int, default=10,
                        help="해석 결과 프레임 저장 간격 (time steps)")
    parser.add_argument("--dynamic-visual-scale", type=float, default=100.0,
                        help="3D 표시 전용 변형 확대배율")
    parser.add_argument("--no-html", action="store_true")
    parser.add_argument("--no-open", action="store_true")
    parser.add_argument("--global-fem", action="store_true")
    parser.add_argument("--local-width", type=float, default=LOCAL_DOMAIN_WIDTH_M)
    parser.add_argument("--local-height", type=float, default=LOCAL_DOMAIN_HEIGHT_M)
    parser.add_argument("--local-nx", type=int, default=LOCAL_N_X)
    parser.add_argument("--local-ny", type=int, default=LOCAL_N_Y)
    parser.add_argument("--foundation-width", type=float, default=FOUNDATION_FOOTPRINT_WIDTH_M)
    parser.add_argument("--foundation-length", type=float, default=FOUNDATION_FOOTPRINT_LENGTH_M)
    return parser.parse_args(argv)


# ============================================================
# 9. 메인 실행
# ============================================================

def main(argv=None):
    args = parse_cli_args(argv)
    if args.inspect_insight_mseed is not None:
        import json
        print(json.dumps(inspect_mseed_headers(args.inspect_insight_mseed), indent=2, ensure_ascii=False))
        return
    if args.dynamic_csv is not None:
        if args.body.lower() != "mars":
            raise ValueError("InSight 지진파 입력은 --body Mars 모델에 한정합니다.")
        if args.no_gravity or not RUN_GRAVITY_ANALYSIS:
            raise ValueError("동적해석 전에 중력 초기화가 필요합니다.")
        if not args.experimental_free_sides:
            raise ValueError("현재 측면 롤러 경계에서 수평 가진 불가. 실험용 --experimental-free-sides가 필요합니다.")
    required_file_report = collect_required_file_status()
    print_required_file_status(required_file_report)

    if args.synthetic_test != "none":
        run_synthetic_smoke_tests(
            kind=args.synthetic_test,
            run_gravity=not args.no_gravity,
        )
        return

    body = args.body.strip().lower()

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
    _, selected_config = get_body_config(body_name)
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
    html_path = None

    if args.no_html:
        print("\nHTML visualization skipped by --no-html.")

    elif mode == "surface":
        print("\nSaving one 3D surface HTML...")
        html_path = save_area_3d_html(
            X_area,
            Y_area,
            Z_area_rel,
            body_name=body_name,
        )

    elif mode == "solid" and not BUILD_OPENSEES_MODEL:
        print("\nSaving DEM-only 3D solid-block HTML...")
        html_path = save_area_3d_solid_html(
            X_area,
            Y_area,
            Z_area_rel,
            body_name=body_name,
        )

    elif mode == "solid":
        print("\nOpenSees 3D model HTML will be saved after the model is built.")

    else:
        raise ValueError(
            "VISUALIZATION_MODE는 'surface' 또는 'solid'이어야 합니다."
        )

    # 선택한 천체 설정으로 OpenSees 3D 지반모델 생성
    if BUILD_OPENSEES_MODEL:
        print(f"\nBuilding {body_name} OpenSees 3D soil model...")

        _, soil_config = get_body_config(body_name)
        if args.dynamic_csv is not None:
            soil_config["boundary_condition_mode"] = "experimental_fixed_base_free_sides"

        local_site_metadata = None
        footprint = None
        if args.foundation_width is not None and args.foundation_length is not None:
            footprint = {
                "width_m": float(args.foundation_width),
                "length_m": float(args.foundation_length),
            }

        if USE_LOCAL_OPENSEES_DOMAIN and not args.global_fem:
            print("\nExtracting local SSI analysis domain...")
            if body_name == "Earth":
                (
                    X_model,
                    Y_model,
                    Z_model,
                    _Z_model_abs,
                    _local_lats,
                    _local_lons,
                ) = extract_area_patch(
                    dem,
                    width_m=args.local_width,
                    height_m=args.local_height,
                    n_x=args.local_nx,
                    n_y=args.local_ny,
                )
                native_dx, native_dy = estimate_hgt_native_spacing_m(
                    dem.shape,
                    AREA_CENTER_LAT,
                )
                local_site_metadata = make_direct_local_site_metadata(
                    source_shape=dem.shape,
                    local_shape=Z_model.shape,
                    center_x=0.5 * args.local_width,
                    center_y=0.5 * args.local_height,
                    width_m=args.local_width,
                    height_m=args.local_height,
                    target_nx=args.local_nx,
                    target_ny=args.local_ny,
                    source_dx_mean_m=native_dx,
                    source_dy_mean_m=native_dy,
                    structure_footprint=footprint,
                )
            else:
                (
                    X_model,
                    Y_model,
                    Z_model,
                    local_site_metadata,
                ) = extract_local_site(
                    X_area,
                    Y_area,
                    Z_area_rel,
                    width_m=args.local_width,
                    height_m=args.local_height,
                    target_nx=args.local_nx,
                    target_ny=args.local_ny,
                    structure_footprint=footprint,
                )
            print("  Local shape:", Z_model.shape)
            print("  Local width:", local_site_metadata["width_m"], "m")
            print("  Local height:", local_site_metadata["height_m"], "m")
        else:
            X_model, Y_model, Z_model, _ = orient_terrain_grid_for_opensees(
                X_area,
                Y_area,
                Z_area_rel,
            )

        if footprint is not None:
            soil_config["foundation_footprint"] = footprint
        if local_site_metadata is not None:
            soil_config["local_site"] = local_site_metadata

        model_data = build_opensees_soil_model_3d(
            X_model,
            Y_model,
            Z_model,
            config=soil_config,
            run_gravity=RUN_GRAVITY_ANALYSIS and not args.no_gravity,
        )

        if local_site_metadata is not None:
            model_data["local_site"] = local_site_metadata

        if (
            args.foundation_width is not None
            and args.foundation_length is not None
        ):
            center_x = 0.5 * (
                float(np.nanmin(X_model)) + float(np.nanmax(X_model))
            )
            center_y = 0.5 * (
                float(np.nanmin(Y_model)) + float(np.nanmax(Y_model))
            )
            footprint_nodes = select_surface_nodes_within_footprint(
                model_data,
                center_x=center_x,
                center_y=center_y,
                width_m=args.foundation_width,
                length_m=args.foundation_length,
            )
            model_data["foundation_footprint_nodes"] = footprint_nodes
            print("  Foundation footprint surface nodes:", len(footprint_nodes))

        if args.dynamic_csv is not None:
            model_data["dynamic_axes_verified_enu"] = bool(args.axes_enu_verified)
            dynamic_summary = run_uniform_excitation(
                model_data, args.dynamic_csv,
                acknowledge_surface_as_base=args.experimental_surface_as_base,
                scale=args.dynamic_scale,
                max_steps=args.dynamic_max_steps,
                output_csv=(OUTPUT_DIR / "insight_experimental_dynamic_response.csv"),
                animation_html=(OUTPUT_DIR / "insight_dynamic_animation.html"),
                animation_stride=args.dynamic_animation_stride,
                visual_scale=args.dynamic_visual_scale,
            )
            print("[EXPERIMENT ONLY]", dynamic_summary)

        if mode == "solid" and not args.no_html:
            html_path = save_opensees_3d_model_html(model_data)

        print(f"{body_name} OpenSees 3D soil model is ready.")

    if html_path is not None:
        print("\nVisualization created:")
        print(" ", html_path)

        if AUTO_OPEN_HTML and not args.no_open:
            webbrowser.open(html_path.as_uri())





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

    html_path = (OUTPUT_DIR / "pohang_profile_2d.html").resolve()

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

    html_path = (OUTPUT_DIR / "pohang_soil_mesh_2d.html").resolve()

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

    html_path = (OUTPUT_DIR / "pohang_soil_mesh_3d.html").resolve()

    fig.write_html(str(html_path), auto_open=True)

    print(f"Saved: {html_path}")
    webbrowser.open(html_path.as_uri())

if __name__ == "__main__":
    main()
