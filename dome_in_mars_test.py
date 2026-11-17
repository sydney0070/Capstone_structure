"""
화성 지오데식 돔 — 링 기초 + 레골리스 덮개 + 인장 불가 지반 + 실측 하중(S1222a, TWINS, PS)
지진: 응답스펙트럼 기반 등가정적하중
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import json
import math
import zipfile
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import openseespy.opensees as ops
import plotly.graph_objects as go


NODE_OFFSET = 1000  # 구조물 절점 태그 = 원래 nodeTag + NODE_OFFSET, 태그 중복 방지

# ----------------------------------------------------------------------
# 데이터 경로 — 환경변수 MARS_DATA_DIR이 있으면 그 경로를 우선 사용
# ----------------------------------------------------------------------
DATA_DIR = Path(os.environ.get("MARS_DATA_DIR", Path.home() / "Desktop" / "MARS_DATA"))
QUAKE_FILE = "S1222a_ACC.mseed"
TWINS_FILES = ("twins_calib_1217_01.csv", "twins_calib_1218_01.csv")
PS_FILES = ("ps_calib_1217_01.csv", "ps_calib_1218_01.csv")
OUTPUT_ROOT = Path("results_mars_dome")   # 결과는 OUTPUT_ROOT/<돔 형상>_<재료>/ 에 저장 (형상·재료별로 안 덮어씀)

# ======================================================================
# [재료 파트] 돔 부재 재료 — 팀 재료 파일(catalog.json)에서 읽는다
#   MEMBER_MATERIAL: catalog.json "materials"의 이름 (예: "al6061_t6_extrusion", "steel_s355_nominal")
#   재료 폴더: 환경변수 MATERIAL_RELEASE_DIR이 있으면 그 경로, 없으면 바탕화면/structure/material_release_20261101
#              (팀원이 새 버전을 보내면 폴더 이름만 바꾸면 됨)
#   catalog에 없는 값은 아래 MEMBER_ALPHA_BY_NAME(열팽창계수)만 코드에서 지정
# ======================================================================
MATERIAL_RELEASE_DIR = Path(os.environ.get(
    "MATERIAL_RELEASE_DIR", Path.home() / "Desktop" / "structure" / "material_release_20261101"))
MATERIAL_CATALOG_FILE = MATERIAL_RELEASE_DIR / "configs" / "catalog.json"
MEMBER_MATERIAL = "al6061_t6_extrusion"   # <- 재료는 여기만 바꾸면 됨

MEMBER_ALPHA_BY_NAME = {                  # 1/K, catalog에 없음 -> 일반 문헌값
    "al6061_t6_extrusion": 2.36e-5,
    "steel_s355_nominal": 1.2e-5,
    "steel_s355_coupon": 1.2e-5,
}
REBAR_DENSITY = 7850.0                    # kg/m^3, 링 철근 (돔 부재 재료와 무관)


def load_member_material(name):
    """catalog.json에서 돔 부재(보 요소)용 E, nu, fy, rho를 읽는다. SI 단위."""
    path = MATERIAL_CATALOG_FILE
    if not path.exists():
        raise FileNotFoundError(f"재료 파일을 찾지 못했습니다: {path}\n"
                                "MATERIAL_RELEASE_DIR(폴더 이름) 또는 같은 이름의 환경변수를 확인하세요.")
    materials = json.loads(path.read_text(encoding="utf-8"))["materials"]
    if name not in materials:
        raise KeyError(f"catalog에 '{name}'이 없습니다. 사용 가능: {', '.join(materials)}")
    rec, p = materials[name], materials[name].get("params", {})
    E, fy = rec.get("E_Pa", p.get("E")), rec.get("fy_Pa", p.get("Fy"))
    if E is None or fy is None or "nu" not in rec:
        raise ValueError(f"'{name}'은 등방성 E·nu·항복강도가 없어 돔 부재(보)로 쓸 수 없습니다 "
                         f"(law={rec.get('law')}).")
    if name not in MEMBER_ALPHA_BY_NAME:
        raise KeyError(f"'{name}'의 열팽창계수가 없습니다. MEMBER_ALPHA_BY_NAME에 추가하세요.")
    return {"name": name, "E": float(E), "nu": float(rec["nu"]), "fy": float(fy),
            "rho": float(rec["density_kg_m3"]), "alpha": MEMBER_ALPHA_BY_NAME[name],
            "status": rec.get("status", ""), "sources": rec.get("sources", []), "file": str(path)}


_MAT = load_member_material(MEMBER_MATERIAL)
MEMBER_E = _MAT["E"]                      # Pa
MEMBER_NU = _MAT["nu"]
MEMBER_FY = _MAT["fy"]                    # Pa
MEMBER_DENSITY = _MAT["rho"]              # kg/m^3
MEMBER_ALPHA = _MAT["alpha"]              # 1/K
# 결과 파일(csv, json, png, html, npz)에 같이 남기는 재료 정보
MEMBER_MATERIAL_LABEL = (f"{MEMBER_MATERIAL} (E {MEMBER_E/1e9:.1f} GPa, fy {MEMBER_FY/1e6:.0f} MPa, "
                         f"ρ {MEMBER_DENSITY:.0f} kg/m³)")
MEMBER_MATERIAL_INFO = {"name": MEMBER_MATERIAL, "E_GPa": MEMBER_E / 1e9, "nu": MEMBER_NU,
                        "fy_MPa": MEMBER_FY / 1e6, "density_kg_m3": MEMBER_DENSITY, "alpha_1_K": MEMBER_ALPHA,
                        "status": _MAT["status"], "sources": _MAT["sources"], "catalog_file": _MAT["file"]}

# ----------------------------------------------------------------------
# 화성 실측 지진 설정 (S1222a 원본 그대로 사용, 배율/PGA 스케일링 없음)
# ----------------------------------------------------------------------
MARS_EVENT = "S1222a"           # InSight 최대 화성지진 (2022-05-04)
DIR_TO_CHANNEL = {1: "BHN", 2: "BHE", 3: "BHZ"}   # X=북남, Y=동서, Z=수직
QUAKE_DT = 0.05                 # mseed 원본 20 Hz와 동일 -> 리샘플링 없음

# ----------------------------------------------------------------------
# 지진 등가정적하중 설정 (응답스펙트럼)
#   지진 기록 전체(20 Hz, 수 시간)를 프로그램 시작 시 한 번만 1자유도계 정확해로 돌려
#   "어떤 고유주기의 구조물이든 받을 수 있는 최대 가속도" Sa_max 하나로 압축한다.
#   이후 해석에서는 지진 = 절점 질량 x Sa_max 의 정적 하중 (바람·순압과 같은 단위하중 방식).
# ----------------------------------------------------------------------
SEISMIC_DAMPING = 0.05                           # 감쇠비 5% (강구조 + 흙 덮개, 내진설계 표준값)
SPECTRUM_PERIODS_S = np.geomspace(0.05, 20.0, 300)   # 스펙트럼 주기 범위 (2·dt ~ 20 s)
SEISMIC_USE_PERIOD_ENVELOPE = True               # True: 주기 무관 최대값 사용 (고유주기 몰라도 보수적)
                                                 # False: 해석마다 1차 고유주기의 Sa(T1) 사용
SEISMIC_HORIZONTAL_ENVELOPE = True               # 수평 두 성분(BHN, BHE) 중 큰 Sa를 X·Y 모두에 적용
SEISMIC_ORTHO_FACTOR = 0.3                       # 직교방향 조합 100% + 30% + 30% (KDS 41 17 00, ASCE 7). 1.0이면 세 성분 절대합
SEISMIC_MODE_STEPS = (60, 100, 125, 150, 200)    # 모드 수를 늘려가며 아래 두 조건을 만족할 때까지 고유치 해석
SEISMIC_RIGID_FACTOR = 1.5                       # 마지막 모드 진동수 >= 1.5 x 기록 나이퀴스트(10 Hz) -> 남은 모드는 강체처럼 거동
SEISMIC_MASS_TARGET = 0.99                       # 방향별 모드 유효질량 합 >= 99% (나머지는 '잔여질량' 하중으로 항상 포함)
EIGEN_MODES = 6                                  # 보고용 고유진동수 개수
AIR_GAMMA = 1.4                                  # 내부 공기 비열비 (빠른 진동 = 단열 압축, 보수적)

# 장시간 진동(S1222a 코다 4~10시간) 피로 상한 검토 — EN 1993-1-9
FATIGUE_DETAIL_CATEGORY_MPA = 36.0 * min(1.0, MEMBER_E / 2.0e11)  # 강재 최불리 상세범주 36 MPa를 E비로 환산 (강재는 36 그대로, 알루미늄 약 12.7 — EN 1999-1-3 정식 값 아님)
FATIGUE_GAMMA_MF = 1.35                          # 피로강도 부분안전계수 (손상허용·파괴 영향 큼)
FATIGUE_GAMMA_FF = 1.0                           # 피로하중 부분안전계수
FATIGUE_DESIGN_LIFE_MARS_YEARS = 20.0            # 설계수명 (화성년)
MARS_YEAR_SOLS = 668.6

# 선택: 최종 설계 비선형 시간이력 비교 (등가정적이 실제 동적 응답보다 큰지 숫자로 보이고 싶을 때 한 번만 True)
#   판정에는 필요 없다 — 등가정적 상한은 선형 응답의 수학적 상한이고, 선형성은 verify_quake_linearity가 매번 확인한다.
#   켜면 원본 20 Hz 6만 스텝을 방향마다 돌리므로 방향당 10분 이상 걸린다.
VERIFY_TIME_HISTORY = False
VERIFY_TH_DIRECTIONS = (1, 3)                    # 1=X(BHN), 2=Y(BHE), 3=Z(BHZ) (돔이 축대칭이라 Y는 X와 같은 거동)

# ----------------------------------------------------------------------
# 실측 환경하중 설정 (TWINS 온도·바람 + PS 기압, Sol 1217~1218 연속 약 23시간)
# ----------------------------------------------------------------------
SOL_SECONDS = 88775.244          # 화성 1솔 = 24시간 39분 35초 (지구 초)
ENV_DT_S = 600.0                 # 환경 시계열 집계 간격 (10분). 기온·기압은 평균, 풍압은 구간 최대
MARS_GAS_CONSTANT = 191.8        # J/(kg·K), 화성 대기 평균 분자량 43.34 g/mol 기준
TWINS_TEMP_VALID_K = (160.0, 300.0)   # 막대온도 물리적 유효범위 (센서 튐 제거)
PS_VALID_PA = (300.0, 1500.0)         # 대기압 유효범위

USE_WIND = True                  # 덮개가 있으면 자동으로 꺼짐 (effective_env)
WIND_DIR_IS_FROM = True          # TWINS 풍향 = 불어오는 방향(북 0°, 시계방향)으로 가정
WIND_SECTORS = 16                # 풍향을 22.5° 간격 16방위로 나눠 단위풍하중을 미리 계산
WIND_REF_HEIGHT_M = 1.2          # TWINS 붐 측정 높이
WIND_ROUGHNESS_M = 0.01          # 지표 거칠기 길이 z0 (가정)
USE_WIND_HEIGHT_PROFILE = True   # 로그 법칙으로 패널 높이별 풍속 보정
WIND_CP_WINDWARD = 0.8           # 풍상측 압력계수 (cosθ에 비례해 감소)
WIND_CP_SUCTION = -0.5           # 풍하측/측면/정상부 흡입계수 (단순 가정)

USE_INTERNAL_PRESSURE = True
INTERNAL_PRESSURE_PA = 101325.0  # 내부 1기압. 순압 = 내부압 - 실측 외기압(t)

USE_THERMAL = True               # 열팽창계수는 [재료 파트]의 MEMBER_ALPHA
THERMAL_REF_TEMP_K = None        # 무응력 기준온도. None이면 타임라인 평균 기온 사용

# ======================================================================
# [외피 파트] 패널 외피 = Combitherm 기밀층 + Kevlar 구속층 (NASA TransHab 그대로)
#   구속층만 하중을 받는다. 기밀층은 구속층보다 크게 재단해 공기만 막는다 (하중·강도 검토 없음)
#   외피는 부재 길이 전체에 고정 -> 내압·덮개·바람은 부재 선하중으로 전달 (5-0절, USE_LINE_LOADS)
# ======================================================================
USE_MEMBRANE = True

# (1) 기밀층: Combitherm 3겹 (폴리에틸렌·나일론·EVOH 적층, TransHab 3중 기밀층)
BLADDER_NAME = "Combitherm"
BLADDER_LAYERS = 3
BLADDER_LAYER_T = 12.75 * 25.4e-6   # m, 1겹 12.75 mil = 0.32 mm (TransHab 기밀층 후보 필름 시험 두께 중 최대, Pedley & Mayeaux 2001)
BLADDER_DENSITY = 1050.0            # kg/m^3, PE·나일론·EVOH 적층 평균 (가정)

# (2) 구속층: Kevlar 1인치 웨빙(12K 등급)을 두 방향으로 엮음 — 하중을 받는 유일한 층 (TransHab SDU II, 4기압 수압시험 통과)
RESTRAINT_NAME = "Kevlar"
RESTRAINT_STRAP_WIDTH = 0.0254      # m, 1인치 (TransHab 표준 웨빙)
RESTRAINT_STRAP_BREAK_N = 53.4e3    # N, 12,000 lbf. 12K급 1인치 Kevlar 웨빙 공칭 12,500 lbf(Jones & Doggett 2014)보다 낮게 잡음 (보수적)
RESTRAINT_STRAP_KG_PER_M = 0.052    # kg/m, 웨빙 1가닥 (Kevlar 49 섬유 3.0 GPa x 직조 효율 0.5, 1,440 kg/m^3로 역산)
RESTRAINT_STRAPS_PER_M = 38         # 가닥/m (한 방향), 고정값. 1인치 웨빙은 1 m에 최대 39가닥 (1기압 지오데식 DCR 약 0.75)
RESTRAINT_PLIES = 1                 # 웨빙 겹 수 (두 방향 1세트 = 1겹). 1 m에 다 안 들어가면 겹을 늘린다
RESTRAINT_SF = 4.0                  # TransHab 구속층 안전율

PANEL_SAG_RATIO = 0.10              # 패널이 부푼 높이 / 패널 폭(내접원 지름). 공기막 쿠션 적정 높이/폭 비 0.1 준용 (Gu et al. 2012). 0 < 값 <= 0.5
USE_LINE_LOADS = True               # 외피 하중(내압·덮개·바람)을 부재 길이 전체에 선하중으로 건다
                                    # (False: 이전처럼 패널 꼭짓점 절점에만 -> 반영 전/후 비교용)

if RESTRAINT_STRAPS_PER_M * RESTRAINT_STRAP_WIDTH >= 1.0:
    raise ValueError("RESTRAINT_STRAPS_PER_M x 웨빙 폭이 1 m 이상입니다 (1인치 웨빙은 1 m에 최대 39가닥). "
                     "가닥 수를 줄이고 RESTRAINT_PLIES를 늘리세요.")
if not 0.0 < PANEL_SAG_RATIO <= 0.5:
    raise ValueError("PANEL_SAG_RATIO는 0보다 크고 0.5(반구) 이하여야 합니다.")

BLADDER_AREAL = BLADDER_LAYERS * BLADDER_DENSITY * BLADDER_LAYER_T
RESTRAINT_AREAL = 2 * RESTRAINT_PLIES * RESTRAINT_STRAPS_PER_M * RESTRAINT_STRAP_KG_PER_M   # 두 방향
MEMBRANE_AREAL_DENSITY = BLADDER_AREAL + RESTRAINT_AREAL                                   # kg/m^2
RESTRAINT_ALLOW_N_M = (RESTRAINT_PLIES * RESTRAINT_STRAPS_PER_M * RESTRAINT_STRAP_BREAK_N
                       / RESTRAINT_SF)                                                     # N/m, 허용 장력
MEMBRANE_NAME = (f"외피 ({BLADDER_NAME} {BLADDER_LAYERS}겹 + {RESTRAINT_NAME} 웨빙 "
                 f"{RESTRAINT_STRAPS_PER_M}가닥/m x {RESTRAINT_PLIES}겹)")
ENVELOPE_INFO = {
    "bladder": {"material": BLADDER_NAME, "layers": BLADDER_LAYERS, "layer_t_mm": BLADDER_LAYER_T * 1e3,
                "density_kg_m3": BLADDER_DENSITY},
    "restraint": {"material": f"{RESTRAINT_NAME} 12K 1-inch webbing", "straps_per_m": RESTRAINT_STRAPS_PER_M,
                  "plies": RESTRAINT_PLIES, "strap_break_kN": RESTRAINT_STRAP_BREAK_N / 1e3, "SF": RESTRAINT_SF,
                  "allow_kN_m": RESTRAINT_ALLOW_N_M / 1e3},
    "areal_kg_m2": {"bladder": BLADDER_AREAL, "restraint": RESTRAINT_AREAL, "total": MEMBRANE_AREAL_DENSITY},
    "panel_sag_ratio": PANEL_SAG_RATIO,
}

# ----------------------------------------------------------------------
# 1솔 타임라인 출력 설정 (피로 해석 입력용)
#   시각은 환경 타임라인 시작(Sol 1217 LMST 07:00) 기준 지구 초
# ----------------------------------------------------------------------
SOL_QUAKE_STARTS_S = [20000.0, 45000.0, 70000.0]  # 1솔 중 S1222a급 지진 발생 시각 목록 (선형 가정으로 같은 응답 재사용)
SAVE_SOL_TIMELINE = True


# ----------------------------------------------------------------------
# 0. 돔 원본 데이터
# ----------------------------------------------------------------------
dome_NODES = {
    1: (13.85640646, -8.0, 12.0),
    2: (8.0, -13.85640646, 12.0),
    3: (4.3589, -7.54983627, 18.0),
    4: (8.7178, 0.0, 18.0),
    5: (16.0, 0.0, 12.0),
    6: (0.0, 0.0, 20.0),
    7: (-4.3589, -7.54983627, 18.0),
    8: (0.0, -16.0, 12.0),
    9: (4.9379568, -18.42870565, 6.0),
    10: (0.0, -20.0, 0.0),
    11: (10.0, -17.32050808, 0.0),
    12: (17.32050808, -10.0, 0.0),
    13: (18.42870565, -4.9379568, 6.0),
    14: (13.85640646, 8.0, 12.0),
    15: (8.0, 13.85640646, 12.0),
    16: (0.0, 16.0, 12.0),
    17: (4.3589, 7.54983627, 18.0),
    18: (13.49074886, -13.49074886, 6.0),
    19: (-8.7178, 0.0, 18.0),
    20: (-4.3589, 7.54983627, 18.0),
    21: (-4.9379568, -18.42870565, 6.0),
    22: (-8.0, -13.85640646, 12.0),
    23: (-13.49074886, -13.49074886, 6.0),
    24: (-10.0, -17.32050808, 0.0),
    25: (20.0, 0.0, 0.0),
    26: (18.42870565, 4.9379568, 6.0),
    27: (4.9379568, 18.42870565, 6.0),
    28: (10.0, 17.32050808, 0.0),
    29: (13.49074886, 13.49074886, 6.0),
    30: (17.32050808, 10.0, 0.0),
    31: (-8.0, 13.85640646, 12.0),
    32: (-13.85640646, 8.0, 12.0),
    33: (-16.0, 0.0, 12.0),
    34: (0.0, 20.0, 0.0),
    35: (-10.0, 17.32050808, 0.0),
    36: (-13.49074886, 13.49074886, 6.0),
    37: (-18.42870565, 4.9379568, 6.0),
    38: (-13.85640646, -8.0, 12.0),
    39: (-18.42870565, -4.9379568, 6.0),
    40: (-20.0, 0.0, 0.0),
    41: (-17.32050808, -10.0, 0.0),
    42: (-17.32050808, 10.0, 0.0),
    43: (-4.9379568, 18.42870565, 6.0),
}
"""(x, y, z) [m] — 40 m급 돔 (밑면 반지름 20 m, 높이 20 m)"""

DOME_EDGES = [
    (1, 1, 2), (2, 3, 2), (3, 4, 3), (4, 4, 5), (5, 1, 5), (6, 3, 1), (7, 6, 3), (8, 6, 7),
    (9, 8, 7), (10, 8, 2), (11, 9, 2), (12, 9, 10), (13, 10, 11), (14, 12, 11), (15, 12, 13),
    (16, 1, 13), (17, 4, 1), (18, 14, 4), (19, 14, 15), (20, 16, 15), (21, 16, 17), (22, 4, 17),
    (23, 6, 4), (24, 3, 7), (25, 8, 3), (26, 9, 8), (27, 5, 14), (28, 5, 13), (29, 18, 13),
    (30, 18, 12), (31, 1, 18), (32, 6, 19), (33, 6, 20), (34, 20, 16), (35, 7, 19), (36, 8, 21),
    (37, 22, 8), (38, 22, 23), (39, 23, 21), (40, 10, 21), (41, 10, 24), (42, 21, 24), (43, 21, 9),
    (44, 9, 11), (45, 11, 18), (46, 18, 9), (47, 2, 18), (48, 25, 12), (49, 13, 25), (50, 13, 26),
    (51, 27, 15), (52, 28, 27), (53, 28, 29), (54, 29, 30), (55, 26, 30), (56, 16, 27), (57, 16, 31),
    (58, 32, 31), (59, 33, 32), (60, 19, 33), (61, 15, 17), (62, 29, 15), (63, 27, 29), (64, 27, 34),
    (65, 35, 34), (66, 36, 35), (67, 37, 36), (68, 32, 37), (69, 7, 22), (70, 38, 7), (71, 39, 38),
    (72, 40, 39), (73, 37, 40), (74, 14, 29), (75, 14, 26), (76, 25, 26), (77, 26, 5), (78, 26, 29),
    (79, 6, 17), (80, 20, 32), (81, 17, 20), (82, 20, 19), (83, 21, 22), (84, 41, 23), (85, 41, 40),
    (86, 40, 42), (87, 36, 42), (88, 23, 24), (89, 23, 39), (90, 33, 39), (91, 30, 25), (92, 43, 27),
    (93, 43, 35), (94, 35, 42), (95, 42, 37), (96, 33, 37), (97, 34, 28), (98, 30, 28), (99, 16, 43),
    (100, 20, 31), (101, 38, 33), (102, 22, 38), (103, 34, 43), (104, 31, 36), (105, 43, 31),
    (106, 39, 37), (107, 38, 23), (108, 39, 41), (109, 17, 14), (110, 19, 32), (111, 19, 38),
    (112, 41, 24), (113, 43, 36), (114, 36, 32),
]
"""(element번호,노드1,노드2):노드1과 노드2를 이은 부재번호"""

# ======================================================================
# [돔 형상 파트] DOME_TYPE 한 줄로 지오데식 / 곡선 리브 돔 전환
#   geodesic   : 위 dome_NODES / DOME_EDGES (곧은 부재를 꺾인 마디로 이어 구를 근사)
#   curved_rib : 같은 반구(R = 20 m) 위에 자오선·링·대각 리브를 구면을 따라 휘게 배치.
#                OpenSees 보 요소는 직선이라 리브 하나를 CURVED_N_SEG개의 짧은 요소로 나눠 곡선을 근사한다.
#                리브 중간 절점은 접합부가 아니라 같은 휜 관의 일부 -> 실제 부재(리브) 수는 늘지 않는다.
# ======================================================================
DOME_TYPE = "geodesic"            # <- "geodesic" 또는 "curved_rib"
CURVED_N_MERIDIANS = 12           # 자오선 리브 수 (지오데식 밑단 지지점 12개와 같게 -> 링 기초 조건 동일)
CURVED_N_BANDS = 4                # 밑단 ~ 정점 층 수 (지오데식과 같은 4층)
CURVED_N_SEG = 4                  # 리브 하나를 나누는 직선 요소 수 (1~9). 늘릴수록 곡선에 가까움
CURVED_DIAGONALS = True           # 사각 패널 대각 리브 (없으면 면내 전단을 휨으로 버텨 약해짐)
CURVED_ELE_BASE = 10000           # 요소 태그 = CURVED_ELE_BASE + 리브번호*10 + 구간번호 (링·지반 요소 4000~5999와 안 겹침)


def make_curved_rib_dome(R, n_meridians, n_bands, n_seg, diagonals=True):
    """중심 원점, 반지름 R인 반구 위 곡선 리브 돔.
      접합부(joint): 층 k(고도각 90°·k/n_bands) x 자오선 i(방위각 360°·i/n_meridians) + 정점
      리브: 자오선(대원호), 링(위도선 원호), 대각(대원호). 각 리브를 n_seg개 직선 요소로 나눈다.
    반환: nodes {id: (x,y,z)}, edges [(ele, i, j)], panels [패널 경계 절점 루프],
          buckling_length {ele: 리브 전체 호 길이}, info {members, ele_member, joints, n_seg}"""
    if not 1 <= n_seg <= 9:
        raise ValueError("CURVED_N_SEG는 1~9 (요소 태그 규칙: 리브번호*10 + 구간번호)")
    if n_meridians < 3 or n_bands < 1:
        raise ValueError("자오선 3개 이상, 층 1개 이상이어야 합니다.")

    def sph(phi, th):
        return (R * math.cos(phi) * math.cos(th), R * math.cos(phi) * math.sin(th), R * math.sin(phi))

    phis = [0.5 * math.pi * k / n_bands for k in range(n_bands)]
    ths = [2.0 * math.pi * i / n_meridians for i in range(n_meridians)]
    nodes, jid = {}, {}
    for k in range(n_bands):
        for i in range(n_meridians):
            jid[k, i] = 1 + k * n_meridians + i
            nodes[jid[k, i]] = sph(phis[k], ths[i])
    apex = n_bands * n_meridians + 1
    nodes[apex] = (0.0, 0.0, float(R))
    joints = set(nodes)

    def J(k, i):
        return apex if k >= n_bands else jid[k, i % n_meridians]

    def great_arc(a, b):          # 대원호 위 중간점 (자오선, 대각)
        u0, u1 = np.asarray(nodes[a]) / R, np.asarray(nodes[b]) / R
        w = math.acos(max(-1.0, min(1.0, float(u0 @ u1))))
        return [R * (math.sin((1 - s / n_seg) * w) * u0 + math.sin(s / n_seg * w) * u1) / math.sin(w)
                for s in range(1, n_seg)]

    def lat_arc(k, i):            # 위도선 원호 위 중간점 (링)
        dth = 2.0 * math.pi / n_meridians
        return [sph(phis[k], ths[i] + dth * s / n_seg) for s in range(1, n_seg)]

    members, path = {}, {}
    next_id = [apex + 1]

    def add_member(kind, a, b, inner):
        chain = [a]
        for p in inner:
            nodes[next_id[0]] = tuple(float(v) for v in p)
            chain.append(next_id[0])
            next_id[0] += 1
        chain.append(b)
        mid = len(members) + 1
        members[mid] = (kind, chain)
        path[a, b], path[b, a] = chain, chain[::-1]

    for i in range(n_meridians):
        for k in range(n_bands):
            add_member("자오선", J(k, i), J(k + 1, i), great_arc(J(k, i), J(k + 1, i)))
    for k in range(n_bands):
        for i in range(n_meridians):
            add_member("링", J(k, i), J(k, i + 1), lat_arc(k, i))
    if diagonals:
        for k in range(n_bands - 1):
            for i in range(n_meridians):      # 섹터마다 방향을 번갈아 -> 좌우 대칭 마름모 패턴
                a, b = (J(k, i), J(k + 1, i + 1)) if i % 2 == 0 else (J(k, i + 1), J(k + 1, i))
                add_member("대각", a, b, great_arc(a, b))

    def loop(*corners):
        out = []
        for a, b in zip(corners, corners[1:] + corners[:1]):
            out += path[a, b][:-1]
        return out

    panels = []
    for k in range(n_bands):
        for i in range(n_meridians):
            p00, p01 = J(k, i), J(k, i + 1)
            if k == n_bands - 1:
                panels.append(loop(p00, p01, apex))
                continue
            p10, p11 = J(k + 1, i), J(k + 1, i + 1)
            if not diagonals:
                panels.append(loop(p00, p01, p11, p10))
            elif i % 2 == 0:
                panels += [loop(p00, p01, p11), loop(p00, p11, p10)]
            else:
                panels += [loop(p00, p01, p10), loop(p01, p11, p10)]

    edges, blen, ele_member = [], {}, {}
    for mid, (kind, chain) in members.items():
        L = sum(math.dist(nodes[p], nodes[q]) for p, q in zip(chain, chain[1:]))
        for s, (p, q) in enumerate(zip(chain, chain[1:]), start=1):
            ele = CURVED_ELE_BASE + mid * 10 + s
            edges.append((ele, p, q))
            blen[ele], ele_member[ele] = L, (mid, s)
    if max(nodes) >= 1000:
        raise ValueError("돔 절점이 태그 범위(1~999)를 넘습니다. 자오선 수, 층 수, CURVED_N_SEG를 줄이세요.")
    return nodes, edges, panels, blen, {"members": members, "ele_member": ele_member, "joints": joints, "n_seg": n_seg}


DOME_PANELS = None                # None: 부재 연결에서 삼각 패널 자동 탐색 (지오데식)
DOME_BUCKLING_LENGTH = None       # None: 요소 길이 = 좌굴 길이 (지오데식)
DOME_MEMBER_INFO = None
if DOME_TYPE == "curved_rib":
    _R_SPHERE = float(np.mean([np.linalg.norm(p) for p in dome_NODES.values()]))   # 지오데식과 같은 구 (20 m)
    dome_NODES, DOME_EDGES, DOME_PANELS, DOME_BUCKLING_LENGTH, DOME_MEMBER_INFO = make_curved_rib_dome(
        _R_SPHERE, CURVED_N_MERIDIANS, CURVED_N_BANDS, CURVED_N_SEG, CURVED_DIAGONALS)
    DOME_N_MEMBERS = len(DOME_MEMBER_INFO["members"])
    DOME_LABEL = (f"곡선 리브 돔 (자오선 {CURVED_N_MERIDIANS}, {CURVED_N_BANDS}층, 대각 {'있음' if CURVED_DIAGONALS else '없음'}, "
                  f"리브 {DOME_N_MEMBERS}개 x 요소 {CURVED_N_SEG}개)")
elif DOME_TYPE == "geodesic":
    DOME_N_MEMBERS = len(DOME_EDGES)
    DOME_LABEL = f"지오데식 돔 (곧은 부재 {DOME_N_MEMBERS}개)"
else:
    raise ValueError("DOME_TYPE은 'geodesic' 또는 'curved_rib'이어야 합니다.")
OUTPUT_DIR = OUTPUT_ROOT / f"{DOME_TYPE}_{MEMBER_MATERIAL}"
DOME_INFO = {"type": DOME_TYPE, "label": DOME_LABEL, "n_nodes": len(dome_NODES), "n_elements": len(DOME_EDGES),
             "n_members": DOME_N_MEMBERS}

ELE_TAGS = [ele for ele, _, _ in DOME_EDGES]   # 부재력 배열의 행 순서


def member_label(ele):
    """보고용 부재 이름. 곡선 리브 돔은 요소가 어느 리브의 몇 번째 구간인지 같이 적는다."""
    if DOME_MEMBER_INFO is None:
        return f"{ele}"
    mid, s = DOME_MEMBER_INFO["ele_member"][ele]
    return f"{ele} ({DOME_MEMBER_INFO['members'][mid][0]} 리브 #{mid}, {s}/{DOME_MEMBER_INFO['n_seg']}구간)"


# ----------------------------------------------------------------------
# 0. 내부 체적 계산 (참고용)
# ----------------------------------------------------------------------
def _find_triangles_from_edges(nodes, edges):
    """엣지 목록에서 삼각형(a,b,c 세 노드가 서로 다 연결된 경우)을 찾는다."""
    adj = defaultdict(set)
    for _, i, j in edges:
        adj[i].add(j)
        adj[j].add(i)

    triangles = set()
    for _, i, j in edges:
        for k in adj[i] & adj[j]:   # i, j에 공통으로 연결된 k -> i-j-k 삼각 패널
            triangles.add(tuple(sorted((i, j, k))))
    return sorted(triangles)        # 정렬해서 실행마다 순서가 같도록 고정


def _base_center_radius(nodes):
    """밑면 중심 (cx, cy, 0)과 구 반지름 (전 절점까지 거리 평균)"""
    P = np.array(list(nodes.values()), dtype=float)
    base = P[np.abs(P[:, 2]) < 1e-6]
    c = np.array([base[:, 0].mean(), base[:, 1].mean(), 0.0])
    return c, float(np.linalg.norm(P - c, axis=1).mean())


_SHELL_CACHE = {}


def _shell_triangles(nodes, edges):
    """외향 정렬된 쉘 삼각형 [(분배 절점들, 분배 비율들, (v0, v1, v2)), ...]
      - 지오데식(3절점 평면 패널): 패널 그대로, 하중은 세 절점에 1/3씩 (기존과 동일)
      - 곡선 리브 돔(경계에 리브 중간 절점이 있는 곡면 패널): 패널 중심을 구면 위로 올린 점에서 부채꼴로 나눠
        곡면을 근사하고, 각 조각의 하중은 조각이 닿은 경계 두 절점에 1/2씩 -> 패널 하중이 리브를 따라 나눠 전달된다"""
    key = (id(nodes), id(edges))
    if key in _SHELL_CACHE:
        return _SHELL_CACHE[key]
    P = {k: np.asarray(v, dtype=float) for k, v in nodes.items()}
    c, R = _base_center_radius(nodes)
    loops = DOME_PANELS if DOME_PANELS is not None else [list(t) for t in _find_triangles_from_edges(nodes, edges)]
    out = []
    for lp in loops:
        V = [P[k] for k in lp]
        g = np.mean(V, axis=0)
        if len(lp) > 3:
            g = c + R * (g - c) / np.linalg.norm(g - c)      # 패널 중심을 곡면 위로
        area2 = sum(np.cross(V[i] - g, V[(i + 1) % len(V)] - g) for i in range(len(V)))
        if np.dot(area2, g - c) < 0:                         # 외향 정렬
            lp, V = lp[::-1], V[::-1]
        if len(lp) == 3:
            out.append((tuple(lp), (1.0 / 3.0,) * 3, (V[0], V[1], V[2])))
        else:
            n = len(lp)
            out += [((lp[i], lp[(i + 1) % n]), (0.5, 0.5), (g, V[i], V[(i + 1) % n])) for i in range(n)]
    _SHELL_CACHE[key] = out
    return out


def compute_enclosed_volume(nodes, edges):
    """돔 쉘이 감싸는 내부 체적(m^3). 밑면 중심을 원점으로 둔 발산정리 V = Σ (1/6) v0·(v1 x v2).
    밑면(z = 0)은 원점을 지나는 평면이라 기여가 0 -> 외향 쉘 삼각형만 더하면 된다."""
    c, _ = _base_center_radius(nodes)
    return abs(sum(float(np.dot(v0 - c, np.cross(v1 - c, v2 - c)))
                   for _, _, (v0, v1, v2) in _shell_triangles(nodes, edges)) / 6.0)


DOME_ENCLOSED_VOLUME_M3 = compute_enclosed_volume(dome_NODES, DOME_EDGES)


# ----------------------------------------------------------------------
# 1. 행성 파라미터
# ----------------------------------------------------------------------
PLANETS = {
    "Mars":  {"g": 3.71, "label": "화성"},
}
CURRENT_PLANET = "Mars"   # 화성데이터(InSight)이므로 화성 고정


def set_planet(name: str):
    global CURRENT_PLANET
    if name not in PLANETS:
        raise ValueError(f"알 수 없는 행성: {name} (실측 환경·지진이 화성(InSight) 자료라 Mars만 가능)")
    CURRENT_PLANET = name
    p = PLANETS[name]
    print(f"[ground] 행성 = {p['label']} (g={p['g']} m/s^2)")


def planet_g():
    return PLANETS[CURRENT_PLANET]["g"]


# ----------------------------------------------------------------------
# 1-2. 실측 환경 데이터 로딩 (TWINS + PS) -> 10분 간격 환경 타임라인
# ----------------------------------------------------------------------
def _require_files(*paths):
    missing = [str(p) for p in paths if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            "입력 파일을 찾지 못했습니다. DATA_DIR 또는 MARS_DATA_DIR 환경변수를 확인하세요.\n  "
            + "\n  ".join(missing))


def _parse_pds_utc(text):
    """PDS UTC 표기 '2022-119T23:21:01.322Z' (연-연중일) -> timezone-aware datetime"""
    return datetime.strptime(str(text).rstrip("Z"), "%Y-%jT%H:%M:%S.%f").replace(tzinfo=timezone.utc)


def _lmst_to_hours(lmst, sol0):
    """'01217M07:00:45.023' -> sol0 자정 기준 연속 시간(h). 솔이 바뀌어도 24, 25...로 이어진다."""
    s = pd.Series(lmst).astype(str)
    sol = s.str.slice(0, 5).astype(int)
    hh = s.str.slice(6, 8).astype(int)
    mm = s.str.slice(9, 11).astype(int)
    ss = s.str.slice(12).astype(float)
    return ((sol - sol0) * 24 + hh + mm / 60.0 + ss / 3600.0).to_numpy(dtype=float)


def format_lmst(hours, sol0):
    """연속 LMST 시간(h) -> 'Sol 1218 03:50'"""
    total_min = int(round(hours * 60.0))
    day, minute = divmod(total_min, 24 * 60)
    return f"Sol {sol0 + day} {minute // 60:02d}:{minute % 60:02d}"


def load_twins(paths):
    """TWINS 보정 데이터 -> 1초 간격 기온·풍속·풍향.

    - 기온: 보정 데이터에는 '대기온도' 열이 없어서, 공기에 가장 가까운 막대 끝단(TIP_ROD)
      온도를 두 붐(BMY, BPY) 평균으로 사용한다. 유효범위 밖 값(센서 튐)은 제거.
    - 바람: 착륙선 뒤쪽(풍하측) 붐은 착륙선 몸체에 가려 풍속이 낮게 측정되므로,
      매 시각 두 붐 중 풍속이 큰 쪽(풍상측)의 풍속과 풍향을 사용한다.
    """
    cols = ["AOBT", "LMST", "UTC",
            "BMY_HORIZONTAL_WIND_SPEED", "BMY_WIND_DIRECTION", "BMY_TIP_ROD_TEMP",
            "BPY_HORIZONTAL_WIND_SPEED", "BPY_WIND_DIRECTION", "BPY_TIP_ROD_TEMP"]
    df = pd.concat([pd.read_csv(p, usecols=cols) for p in paths], ignore_index=True)
    df = df.sort_values("AOBT").drop_duplicates("AOBT").reset_index(drop=True)
    sol0 = int(str(df["LMST"].iloc[0])[:5])

    lo, hi = TWINS_TEMP_VALID_K
    tips_raw = df[["BMY_TIP_ROD_TEMP", "BPY_TIP_ROD_TEMP"]].astype(float)
    tips = tips_raw.where((tips_raw >= lo) & (tips_raw <= hi))
    n_rejected = int(tips.isna().sum().sum() - tips_raw.isna().sum().sum())

    v = df[["BMY_HORIZONTAL_WIND_SPEED", "BPY_HORIZONTAL_WIND_SPEED"]].to_numpy(dtype=float)
    d = df[["BMY_WIND_DIRECTION", "BPY_WIND_DIRECTION"]].to_numpy(dtype=float)
    pick = np.argmax(np.where(np.isnan(v), -np.inf, v), axis=1)   # 풍속이 큰 붐
    rows = np.arange(len(df))

    out = pd.DataFrame({
        "aobt": df["AOBT"].to_numpy(dtype=float),
        "lmst_h": _lmst_to_hours(df["LMST"], sol0),
        "T_air": tips.mean(axis=1).to_numpy(),
        "V": v[rows, pick],
        "dir_from": d[rows, pick],
    })
    meta = {
        "sol0": sol0,
        "aobt_last": float(out["aobt"].iloc[-1]),
        "lmst_last": float(out["lmst_h"].iloc[-1]),
        "utc_first": _parse_pds_utc(df["UTC"].iloc[0]),
        "utc_last": _parse_pds_utc(df["UTC"].iloc[-1]),
        "temp_rejected": n_rejected,
        "n_samples": len(out),
    }
    return out, meta


def load_ps(paths):
    """PS 보정 데이터 -> 10 Hz 대기압(Pa)."""
    df = pd.concat([pd.read_csv(p, usecols=["AOBT", "PRESSURE"]) for p in paths], ignore_index=True)
    df = df.sort_values("AOBT").drop_duplicates("AOBT").reset_index(drop=True)
    lo, hi = PS_VALID_PA
    df = df[(df["PRESSURE"] >= lo) & (df["PRESSURE"] <= hi)]
    return pd.DataFrame({"aobt": df["AOBT"].to_numpy(dtype=float),
                         "P": df["PRESSURE"].to_numpy(dtype=float)})


def _wind_sector(az_to_deg):
    step = 360.0 / WIND_SECTORS
    return (np.rint(np.asarray(az_to_deg) / step).astype(int)) % WIND_SECTORS


def build_env_timeline(twins, twins_meta, ps, dt=ENV_DT_S):
    """TWINS(1 Hz)와 PS(10 Hz)를 기내시계(AOBT)로 맞춰 dt 간격 환경 타임라인을 만든다.

    매 1초 샘플에서 ρ = P/(R·T), q = ½ρV² 를 먼저 계산한 뒤 구간별로
      - 기온, 외기압: 평균
      - 풍압 q: 구간 내 최대값(순간 최대풍) + 그 순간의 풍향
    을 취한다. 바람은 짧은 돌풍이 하중을 지배하므로 평균 대신 최대를 쓴다(보수적).
    """
    t0 = max(twins["aobt"].iloc[0], ps["aobt"].iloc[0])
    t1 = min(twins["aobt"].iloc[-1], ps["aobt"].iloc[-1])
    tw = twins[(twins["aobt"] >= t0) & (twins["aobt"] <= t1)].copy()
    pr = ps[(ps["aobt"] >= t0) & (ps["aobt"] <= t1)].copy()

    tw["P"] = np.interp(tw["aobt"], pr["aobt"], pr["P"])
    tw["rho"] = tw["P"] / (MARS_GAS_CONSTANT * tw["T_air"])
    tw["q"] = 0.5 * tw["rho"] * tw["V"] ** 2
    tw["bin"] = ((tw["aobt"] - t0) // dt).astype(int)
    pr["bin"] = ((pr["aobt"] - t0) // dt).astype(int)

    g = tw.groupby("bin")
    w = pd.DataFrame({
        "aobt": g["aobt"].mean(),
        "lmst_h": g["lmst_h"].mean(),
        "T_K": g["T_air"].mean(),
        "V_mean": g["V"].mean(),
    })
    w["P_Pa"] = pr.groupby("bin")["P"].mean()

    wind = tw.dropna(subset=["q", "dir_from"])
    gust = wind.loc[wind.groupby("bin")["q"].idxmax(), ["bin", "q", "V", "dir_from"]].set_index("bin")
    w = w.join(gust.rename(columns={"q": "q_Pa", "V": "V_gust", "dir_from": "dir_from_deg"}))

    n_missing = int(w[["T_K", "P_Pa", "q_Pa"]].isna().any(axis=1).sum())
    num_cols = ["T_K", "P_Pa", "V_mean", "q_Pa", "V_gust"]
    w[num_cols] = w[num_cols].interpolate(limit_direction="both")
    w["dir_from_deg"] = w["dir_from_deg"].ffill().bfill()   # 각도는 선형보간하지 않는다

    T = w["T_K"].to_numpy()
    P = w["P_Pa"].to_numpy()
    T_ref = float(np.mean(T)) if THERMAL_REF_TEMP_K is None else float(THERMAL_REF_TEMP_K)
    dir_from = w["dir_from_deg"].to_numpy()
    az_to = (dir_from + 180.0) % 360.0 if WIND_DIR_IS_FROM else dir_from % 360.0
    sector = _wind_sector(az_to)

    p_in = INTERNAL_PRESSURE_PA if USE_INTERNAL_PRESSURE else 0.0
    dp = p_in - P if USE_INTERNAL_PRESSURE else np.zeros_like(P)

    env = {
        "t_s": w["aobt"].to_numpy() - t0,
        "aobt": w["aobt"].to_numpy(),
        "lmst_h": w["lmst_h"].to_numpy(),
        "T_K": T,
        "P_Pa": P,
        "rho": P / (MARS_GAS_CONSTANT * T),
        "V_mean": w["V_mean"].to_numpy(),
        "V_gust": w["V_gust"].to_numpy(),
        "q_Pa": w["q_Pa"].to_numpy(),
        "dir_from_deg": dir_from,
        "az_to_deg": az_to,
        "sector": sector,
        "dT": T - T_ref,
        "T_ref": T_ref,
        "dp": dp,
        "dp_ref": float(np.mean(dp)),
        "dt": float(dt),
        "duration_s": float(t1 - t0),
        "wind_active": USE_WIND,
        "thermal_factor": 1.0,
        "meta": {**twins_meta, "t0_aobt": float(t0), "missing_bins": n_missing},
    }
    return env


def print_env_summary(env):
    m = env["meta"]
    s0 = m["sol0"]
    print(f"[env] 구간: {format_lmst(env['lmst_h'][0], s0)} ~ {format_lmst(env['lmst_h'][-1], s0)} "
          f"(지구시간 {env['duration_s']/3600:.2f} h, {len(env['t_s'])}개 x {env['dt']/60:.0f}분)")
    print(f"[env] 기온 {env['T_K'].min():.1f} ~ {env['T_K'].max():.1f} K "
          f"(기준 {env['T_ref']:.1f} K -> ΔT {env['dT'].min():+.1f} ~ {env['dT'].max():+.1f} K)")
    print(f"[env] 외기압 {env['P_Pa'].min():.1f} ~ {env['P_Pa'].max():.1f} Pa, "
          f"대기밀도 {env['rho'].min():.4f} ~ {env['rho'].max():.4f} kg/m^3")
    print(f"[env] 순간최대풍속 {np.max(env['V_gust']):.1f} m/s, 최대 풍압(측정높이) {np.max(env['q_Pa']):.2f} Pa")
    if USE_INTERNAL_PRESSURE:
        print(f"[env] 순압 {env['dp'].min():,.0f} ~ {env['dp'].max():,.0f} Pa (평균 {env['dp_ref']:,.0f} Pa)")
    print(f"[env] 센서 튐 제거 {m['temp_rejected']}개, 결측 구간 보간 {m['missing_bins']}개")


def save_env_csv(env, path):
    s0 = env["meta"]["sol0"]
    df = pd.DataFrame({
        "t_s": env["t_s"], "LMST": [format_lmst(h, s0) for h in env["lmst_h"]],
        "T_K": env["T_K"], "dT_K": env["dT"], "P_Pa": env["P_Pa"], "rho_kg_m3": env["rho"],
        "V_mean_m_s": env["V_mean"], "V_gust_m_s": env["V_gust"], "q_Pa": env["q_Pa"],
        "wind_from_deg": env["dir_from_deg"], "wind_sector": env["sector"], "dp_Pa": env["dp"],
    })
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    print(f"[env] 환경 타임라인 저장: {path.resolve()}")


def utc_to_lmst_hours(env, utc_dt):
    """지구 UTC 시각 -> 환경 데이터와 같은 기준의 연속 LMST 시간(h). 1솔 = SOL_SECONDS로 외삽."""
    m = env["meta"]
    utc_minus_aobt = m["utc_last"].timestamp() - m["aobt_last"]
    aobt = utc_dt.timestamp() - utc_minus_aobt
    return m["lmst_last"] + (aobt - m["aobt_last"]) / SOL_SECONDS * 24.0


def env_indices_at_local_time(env, lmst_start_h, duration_s):
    """하루 중 같은 지방시(hour-of-day)에 해당하는 환경 구간 인덱스.
    S1222a가 일어난 지방시의 환경 조건을 타임라인에서 찾을 때 사용한다."""
    dur_h = duration_s / SOL_SECONDS * 24.0
    half_bin_h = env["dt"] / SOL_SECONDS * 24.0 / 2.0
    hod = env["lmst_h"] % 24.0
    rel = (hod - (lmst_start_h % 24.0)) % 24.0
    return np.where((rel <= dur_h + half_bin_h) | (rel >= 24.0 - half_bin_h))[0]


# ----------------------------------------------------------------------
# 1-3. 실측 지진 로딩 (mseed -> 가속도 배열)
# ----------------------------------------------------------------------
def load_ground_motion(path, channel, dt=QUAKE_DT):
    """S1222a_ACC.mseed에서 한 성분을 읽는다. 파일은 이미 가속도(m/s^2)로 변환된 상태.
    평균만 제거하고 크기는 그대로 둔다(스케일링 없음)."""
    from obspy import read

    st = read(str(path)).select(channel=channel)
    if len(st) == 0:
        raise ValueError(f"{path}에 {channel} 채널이 없습니다.")
    st.merge(method=1, fill_value="interpolate")
    tr = st[0].copy()
    tr.detrend("demean")
    if not math.isclose(tr.stats.delta, dt, rel_tol=1e-6):
        tr.filter("lowpass", freq=0.4 / dt, zerophase=True)   # 앨리어싱 방지 후 리샘플링
        tr.resample(1.0 / dt)

    acc = tr.data.astype(float)
    t = np.arange(len(acc)) * dt
    i_pk = int(np.argmax(np.abs(acc)))
    start = tr.stats.starttime.datetime.replace(tzinfo=timezone.utc)
    info = {
        "id": tr.id,
        "start_utc": start,
        "peak_utc": start + timedelta(seconds=float(t[i_pk])),
        "pga": float(abs(acc[i_pk])),
        "duration_s": float(t[-1]),
    }
    return t, acc, info


# ----------------------------------------------------------------------
# 1-4. 응답스펙트럼 — 지진 기록 전체를 "최대 가속도 하나"로 압축 (프로그램 시작 시 1회)
# ----------------------------------------------------------------------
SPECTRUM_PEAK_SAMPLING = 0.03   # 응답을 고유주기의 3%보다 촘촘히 읽어 샘플 사이 최대값을 놓치지 않게 함 (오차 <= 0.5%)


def _sdof_peak_disp(acc, dt, period, xi):
    """질량 1인 1자유도계 u'' + 2ξω u' + ω² u = -a_g(t) 의 최대 |u|.
    지반가속도를 샘플 사이에서 직선으로 잇는다고 보고 '정확해'로 푼다
    (Nigam & Jennings 1969과 같은 구간선형 가정, 1차 홀드 이산화). dt가 고유주기에 비해 커도 적분 오차가 없다.
    짧은 주기에서는 응답이 샘플 사이에서 정점을 찍을 수 있으므로, 같은 직선 가정으로 입력을 잘게 나눠
    (고유주기의 SPECTRUM_PEAK_SAMPLING 이하 간격) 최대값을 읽는다."""
    from scipy import signal

    n_sub = max(1, int(math.ceil(dt / (SPECTRUM_PEAK_SAMPLING * period))))
    h = dt / n_sub
    p = np.concatenate([[0.0], -np.asarray(acc, dtype=float)])   # t=0에 정지 상태에서 출발 (0 -> a_0 직선)
    if n_sub > 1:
        p = np.interp(np.arange((len(p) - 1) * n_sub + 1) * h, np.arange(len(p)) * dt, p)
    w = 2.0 * math.pi / period
    A = np.array([[0.0, 1.0], [-w * w, -2.0 * xi * w]])
    B = np.array([[0.0], [1.0]])
    C = np.array([[1.0, 0.0]])
    D = np.array([[0.0]])
    Ad, Bd, Cd, Dd, _ = signal.cont2discrete((A, B, C, D), h, method="foh")
    num, den = signal.ss2tf(Ad, Bd, Cd, Dd)
    u = signal.lfilter(num[0], den, p)
    return float(np.max(np.abs(u)))


def response_spectrum(acc, dt, periods=SPECTRUM_PERIODS_S, xi=SEISMIC_DAMPING):
    """의사가속도 응답스펙트럼 PSA(T) = ω²·max|u| (m/s^2). 기록 길이와 무관하게 수 초면 끝난다."""
    periods = np.asarray(periods, dtype=float)
    sd = np.array([_sdof_peak_disp(acc, dt, T, xi) for T in periods])
    psa = (2.0 * math.pi / periods) ** 2 * sd
    return {"T": periods, "Sd": sd, "PSA": psa}


def spectral_value(spec, T):
    """주기 T의 Sa (로그 주기 보간). 범위 밖은 끝값 (짧은 주기 쪽은 PGA에 수렴)."""
    lt = np.log(spec["T"])
    return float(np.interp(math.log(max(T, 1e-6)), lt, spec["PSA"]))


def build_seismic_input(quake_path, dt=QUAKE_DT, xi=SEISMIC_DAMPING, verbose=True):
    """세 성분(BHN/BHE/BHZ)을 읽고 응답스펙트럼과 설계 가속도 Sa를 정한다.
      Sa_env(성분) = max(PGA, 모든 주기의 PSA)  -> 고유주기와 무관한 상한
      설계 Sa: X·Y = (수평 두 성분 중 큰 값) [SEISMIC_HORIZONTAL_ENVELOPE], Z = BHZ"""
    comps, records = {}, {}
    for d, ch in DIR_TO_CHANNEL.items():
        t, acc, info = load_ground_motion(quake_path, ch, dt=dt)
        spec = response_spectrum(acc, dt, xi=xi)
        k = int(np.argmax(spec["PSA"]))
        comps[d] = {"channel": ch, "info": info, "spec": spec, "pga": info["pga"],
                    "Sa_env": max(info["pga"], float(spec["PSA"][k])), "T_peak": float(spec["T"][k]),
                    "amp": max(info["pga"], float(spec["PSA"][k])) / max(info["pga"], 1e-30)}
        records[d] = (t, acc)
    Sa_own = {d: comps[d]["Sa_env"] for d in comps}
    Sa = dict(Sa_own)
    if SEISMIC_HORIZONTAL_ENVELOPE:
        Sa[1] = Sa[2] = max(Sa_own[1], Sa_own[2])
    seis = {"comps": comps, "records": records, "Sa": Sa, "Sa_own": Sa_own, "xi": xi, "dt": dt,
            "nyquist_hz": 0.5 / dt, "source": "period_envelope" if SEISMIC_USE_PERIOD_ENVELOPE else "Sa(T1)"}
    if verbose:
        print_seismic_summary(seis)
    return seis


def print_seismic_summary(seis):
    """입력 지진 요약 (성분별 PGA, Sa_max). 설계 가속도는 최종 보고(print_seismic_report)에서 한 번만 출력"""
    print(f"[seismic] {MARS_EVENT} 응답스펙트럼 (감쇠 {seis['xi']*100:.0f}%, 주기 {SPECTRUM_PERIODS_S[0]:.2f}~"
          f"{SPECTRUM_PERIODS_S[-1]:.0f} s, 구간선형 정확해, 원본 {1/seis['dt']:.0f} Hz 그대로)")
    for d, c in seis["comps"].items():
        print(f"  {c['channel']}: 길이 {c['info']['duration_s']/3600:.2f} h, PGA {c['pga']:.3e} m/s^2, "
              f"Sa_max {c['Sa_env']:.3e} m/s^2 (주기 {c['T_peak']:.2f} s, 증폭 {c['amp']:.2f}배)")


def plot_response_spectrum(seis, out_path, T1=None):
    """세 성분 응답스펙트럼 + 설계 가속도 + (있으면) 돔 1차 고유주기."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _setup_korean_font(plt)

    colors = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a"}
    names = {1: "BHN (X, 북남)", 2: "BHE (Y, 동서)", 3: "BHZ (Z, 수직)"}
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for d, c in seis["comps"].items():
        ax.loglog(c["spec"]["T"], c["spec"]["PSA"], color=colors[d], lw=2.0, label=names[d])
    sa_h, sa_z = seis["Sa"][1], seis["Sa"][3]
    x_lab = SPECTRUM_PERIODS_S[-1] * 0.95
    ax.axhline(sa_h, color="#52514e", lw=1.2, ls="--")
    ax.text(x_lab, sa_h * 1.08, f"설계 수평 Sa = {sa_h:.2e}", va="bottom", ha="right", fontsize=9, color="#52514e")
    ax.axhline(sa_z, color="#52514e", lw=1.2, ls=":")
    ax.text(x_lab, sa_z * 1.08, f"설계 수직 Sa = {sa_z:.2e}", va="bottom", ha="right", fontsize=9, color="#52514e")
    if T1 is not None:
        ax.axvline(T1, color="#0b0b0b", lw=1.0, alpha=0.6)
        ax.text(T1, ax.get_ylim()[0], f" 돔 1차 주기 {T1:.2f} s", rotation=90, va="bottom", fontsize=9)
    from matplotlib.ticker import FuncFormatter, NullFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))     # 로그축 지수 표기(유니코드 마이너스) 경고 방지
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0e}"))
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("고유주기 T (s)")
    ax.set_ylabel("의사가속도 PSA (m/s²)")
    ax.set_title(f"{MARS_EVENT} 응답스펙트럼 (감쇠 {seis['xi']*100:.0f}%)", fontsize=12)
    ax.grid(alpha=0.25, which="major")
    ax.grid(alpha=0.08, which="minor")
    ax.set_ylim(top=max(sa_h, sa_z) * 2.0)
    ax.legend(fontsize=9, frameon=False, loc="lower left")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[seismic] 응답스펙트럼 그림 저장: {Path(out_path).resolve()}")


def save_spectrum_csv(seis, path):
    T = SPECTRUM_PERIODS_S
    df = pd.DataFrame({"T_s": T, **{f"PSA_{c['channel']}_m_s2": c["spec"]["PSA"] for c in seis["comps"].values()}})
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


# ----------------------------------------------------------------------
# 2. 기초 · 덮개 · 지반 설정
# ----------------------------------------------------------------------
FOUNDATION_MODE = "ring"   # "ring"  : 링 기초 + 인장 불가 지반 (화성 현실 모델)
                           # "none"  : 기초 없음. 돔 밑단 절점이 원지반 위에 그냥 놓임 (들림은 덮개+외피+돔 무게로만 버틸 경우)
                           
if FOUNDATION_MODE not in ("ring", "none"):
    raise ValueError("FOUNDATION_MODE는 'ring' 또는 'none'이어야 합니다.")

# 2-1. 원지반 레골리스 (링이 묻히는 흙)
SITE_REGOLITH_DENSITY = 1300.0    # kg/m^3, InSight 착륙지 추정 약 1.3 Mg/m^3
SITE_PHI_DEG = 30.0               # 내부마찰각 (가정, 민감도 대상)
SITE_COHESION_PA = 0.0            # InSight 착륙지: 무점착 ~ 수 kPa 미만 -> 0 (보수적)
SITE_UPLIFT_K = 0.0               # 인발 측면 마찰 제외. TM 5-818-1 6-4 식 6-1: 무게만으로 FS 1.5 
SITE_SUBGRADE_K = 2.0e7           # N/m^3, 지반반력계수 (가정 -> 팀원 지반 코드 값으로 교체)
SITE_H_TO_V_STIFFNESS = 0.5       # 수평/수직 지반 강성비 (가정)
SITE_BASE_FRICTION_RATIO = 2.0 / 3.0   # 바닥 마찰각 δ = (2/3)φ (관례)
SITE_REGOLITH_THICKNESS_M = 3   # 모래질 레골리스 두께 약 3 m, 그 아래 약 15 m 암괴층 (Hobiger et al. 2021,
                                  # Nat. Commun. 12:6756). 링 깊이 상한 -> 링이 φ·ρ 가정한 레골리스 안에만 놓임
UPLIFT_MOBILIZE_M = 0.02          # 플랜지 위 흙 인발 저항이 다 발휘되는 변위 (가정)

# 2-2. 레골리스 덮개
COVER_DENSITY = 1450.0            # MMS-1 상대밀도 80% 건조밀도 (Fujikawa et al., 2019)
COVER_MIN_THICKNESS_M = 2.0       # 방사선 차폐 하한 (Zhang et al. 2022: 1~2 m -> 상한값)
COVER_INERTIA = True              # 덮개 질량을 지진 관성에 포함 (보수적). False면 무게만 하중으로
COVER_T_STEP_M = 0.1              # t_min 탐색 간격
COVER_T_MAX_M = 15.0              # t_min 탐색 상한
DIURNAL_SKIN_DEPTH_M = 0.05       # 일교차 열침투 깊이 δ = sqrt(κP/π), κ ≈ 3e-8 ~ 1e-7 m^2/s (가정)
                                  # 덮개 아래 부재 ΔT = 기온 ΔT x exp(-t/δ)

# 2-3. 링 기초 (황 콘크리트, 역T 단면). 링 윗면 = 지표 = 돔 밑단
RING_N_SUB = 4                    # 밑단 노드 사이 원호 등분 수 (12 x 4 = 48개 링 노드)
RING_TOTAL_NODES = 48             # 링 절점 목표 수. 밑단 절점이 많은 곡선 리브 돔도 같은 48개로 맞춘다
RING_STEM_B = 1.0                 # m, 기둥부 폭 (돔 밑단이 앉는 부분)
RING_FLANGE_T = 0.5               # m, 바닥판 최소 두께 (휨 검토에서 필요하면 자동으로 두꺼워짐)
RING_B_GRID = tuple(float(b) for b in np.round(np.arange(2.0, 15.01, 1.0), 3))   # 탐색할 바닥판 폭 B (m). 0.5 간격이면 시간 2배

# 2-3-1. 바닥판(기둥부 양옆 외팔보) 휨
FLANGE_MODE = "plain"             # "plain": 무근 레골리스 황 콘크리트만 (기초에 금속을 쓰지 않음). 안 들어가면 그 칸 탈락
                                  # "auto" : 무근으로 두께를 키우고, 깊이 안에 안 들어가면 철근으로 보강
                                  # "rebar": 두께 RING_FLANGE_T 고정 + 철근
FLANGE_T_STEP = 0.1               # m, 바닥판 두께 올림 간격
STEM_H_MIN = 0.5                  # m, 기둥부 최소 높이. D - 바닥판 두께 >= 이 값이어야 링이 성립
FS_FLANGE = 1.875                 # 바닥판 휨 안전율. 무근: EN 1992-1-1 12.3.1 (γc 1.5 / αct,pl 0.8) -> f_t / FS
FLANGE_REBAR_COVER = 0.075        # m, 철근 피복 (유효깊이 d = 두께 - 피복)

RING_MATERIAL_CASES = {
    #        f_c (Pa)  f_t (Pa)   E (Pa)   rho     근거
    "낙관": (50.0e6,  3.9e6,    5.0e9,   2200.0),  # Wan et al. 2016 (f_c, 쪼갬인장 f_t 실측), E·rho: Okutsu et al. 2024 (황 50%)
    "기준": (25.0e6,  2.5e6,    5.0e9,   2200.0),  # Snehal et al. 2024 (27 MPa), f_t = f_c/10 가정, E·rho: Okutsu et al. 2024
    "비관": (3.7e6,   0.37e6,   9.6e9,   2200.0),  # Tute & Goulas 2024 (CO2 조건, E 실측), f_t = f_c/10 가정, rho: Okutsu et al. 2024
}
RING_MATERIAL = "기준"
RING_NU = 0.2                     # EN 1992-1-1 3.1.3: 균열 없는 콘크리트
RING_REBAR_FY = 500.0e6           # Pa, 링 인장 보강 철근 (무근 모드에서는 미사용)

# 2-4. 기초 안전율
FS_UPLIFT = 1.5                   # TM 5-818-1 (1983) 6-4: 지속 들림 하중, 무게 기준
FS_BEARING = 3.0                  # AASHTO (1996)
FS_SLIDING = 1.5                  # AASHTO (1996)
# (돔 질량 600 t 가정은 없앰: 칸마다 0 t로 시작해 실제 최소 질량으로 기초를 다시 설계, 바뀌지 않을 때까지 반복)

# 2-5. 탐색 설정
RING_DEPTH_GRID = (1.0, 2.0, 3.0)            # 링 깊이 D 후보 (D <= 레골리스 두께). D x B 모든 칸을 계산한다
CELL_ITER_MAX = 6                             # 칸 하나에서 '기초 설계 <-> 돔 부재 질량' 반복 상한 (보통 2~3바퀴에 수렴)
                                              # 상한까지 수렴하지 않은 칸은 결과에서 제외한다

# 2-5-1. 결과 정리 기준 (가중치 없음)
#   두 가지 재료를 따로 본다: ① 돔 재료 = 돔 부재 + 철근 (t), ② 기초 재료 = 링 황 콘크리트 (t)
#   파레토 해 = 두 값 모두 더 작은 다른 칸이 없는 칸 (탐색 결과 표·그림에 기록)
#     material_min : 돔 재료 최소안 (기초는 그만큼 커짐)  <- 연구 목적 = 돔 부재 최소 -> 최종 설계
#     concrete_min : 기초 재료 최소안 (돔은 그만큼 무거워짐) <- 비교용. 넣으면 같은 방식으로 함께 해석
#   돔 재료가 같은 칸이 여럿이면 링 콘크리트가 적은 칸을 고른다
#   덮개·굴착 부피는 따로 기록만 한다 (재료 종류가 달라 합치지 않음)
FINAL_ENDS = ("material_min",)                # 비교까지 보려면 ("material_min", "concrete_min")
END_LABELS = {"material_min": "돔 재료 최소안", "concrete_min": "기초 재료 최소안"}

# 2-5-2. 침하 (감압 LC2, 링 지지점 연직 변위). 지반반력계수 SITE_SUBGRADE_K가 가정값이므로 결과도 그 범위에서 해석
SETTLEMENT_LIMIT_M = 0.05                     # 최대 침하 (얕은 기초 관례 25~50 mm 중 상한)
ANGULAR_DISTORTION_LIMIT = 1.0 / 500.0        # 이웃 지지점 사이 부등침하 / 거리 (Bjerrum 1963: 균열 방지 한계)
SEARCH_WITH_QUAKE = False                     # False: 탐색은 지진 제외(LC1·LC2·환경)로 빠르게, 최적안만 지진 상한으로 재확인
                                              #        (지진 여유배율 λ >> 1 이면 지진이 최소 질량을 바꾸지 않음 -> 결과 동일)
SEARCH_VERBOSE = False                        # True: 질량 후보마다 판정 줄 + 해석 상세 로그
                                              # False: 칸마다 '[cell] n바퀴'와 '[cell] 완료' 줄만
MASS_SEARCH_GUESS_STEP_TON = 20.0             # 이전 결과(2바퀴째, 옆 칸, 지진 확인) 근처에서 다시 찾을 때 브라케팅 폭 (가까우니 좁게)
FOUNDATION_RETRY_STEP_M = 0.5                 # 해석의 절점별 기초 검토에서 탈락하면 덮개를 이만큼 올려 재시도
FOUNDATION_RETRY_MAX = 6

# 2-6. 태그 범위 (돔 절점 1001~1043, 돔 부재 1~114와 겹치지 않게)
RING_NODE_BASE = 2000
GROUND_NODE_BASE = 3000
RING_ELE_BASE = 4000
GROUND_ELE_BASE = 5000
GROUND_MAT_BASE = 100000          # 절점 i의 재료 = GROUND_MAT_BASE + 10*i + (1~8)


@dataclass(frozen=True)
class Design:
    """기초·덮개 설계값. 부재 질량은 이분 탐색이 따로 정한다."""
    ring_depth: float = 0.0       # m, 지표 ~ 링 바닥 (D)
    cover_t: float = 0.0          # m, 덮개 두께 (t)
    ring_b: float = 0.0           # m, 바닥판 폭 (B)
    flange_t: float = RING_FLANGE_T   # m, 바닥판 두께 (휨 검토로 결정)
    flange_rebar_kg: float = 0.0  # 바닥판 철근. 무근이면 0

    def label(self):
        if FOUNDATION_MODE != "ring":
            return f"{FOUNDATION_MODE}_t{self.cover_t:.1f}"
        return f"D{self.ring_depth:.1f}_B{self.ring_b:.1f}_t{self.cover_t:.1f}_f{self.flange_t:.1f}"


class FoundationInfeasible(Exception):
    """부재 질량을 키워도 해결되지 않는 기초 탈락 (들림, 지지력, 활동, 링 압괴)."""

    def __init__(self, design, reason, mass_ton=None):
        super().__init__(f"{design.label()}: {reason}")
        self.design = design
        self.reason = reason
        self.mass_ton = mass_ton      # 탈락한 질량 후보 (재시도 탐색을 여기서 시작)


def get_base_node_tags(nodes, tol=1e-6):
    """z=0인 노드(돔 밑단 링)의 원래 nodeTag 목록을 자동으로 찾는다."""
    return [tag for tag, (x, y, z) in nodes.items() if abs(z) < tol]


# ----------------------------------------------------------------------
# 2-7. 링 기초 형상 · 단면 · 지반 저항
# ----------------------------------------------------------------------
def ring_layout(nodes, n_sub=RING_N_SUB):
    """돔 밑단 노드를 각도 순으로 잇고, 사이 원호를 n_sub 등분한 링 절점 배치.
    반환: points [(x, y, 돔밑단태그 or None)], L_trib(각 절점이 맡는 원호 길이), R, arc_length, chord_length"""
    c, _ = _base_center_radius(nodes)
    base = get_base_node_tags(nodes)
    ang = {t: math.atan2(nodes[t][1] - c[1], nodes[t][0] - c[0]) for t in base}
    base.sort(key=lambda t: ang[t])
    R = float(np.mean([math.hypot(nodes[t][0] - c[0], nodes[t][1] - c[1]) for t in base]))

    points, theta = [], []
    for i, t in enumerate(base):
        a0, a1 = ang[t], ang[base[(i + 1) % len(base)]]
        while a1 <= a0:
            a1 += 2.0 * math.pi
        points.append((nodes[t][0], nodes[t][1], t))
        theta.append(a0)
        for j in range(1, n_sub):
            a = a0 + (a1 - a0) * j / n_sub
            points.append((c[0] + R * math.cos(a), c[1] + R * math.sin(a), None))
            theta.append(a)

    th = np.array(theta)
    d_next = np.diff(np.append(th, th[0] + 2.0 * math.pi))
    d_prev = np.roll(d_next, 1)
    L_trib = R * (d_next + d_prev) / 2.0
    xy = np.array([(x, y) for x, y, _ in points])
    chord = float(np.sum(np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1)))
    return {"points": points, "L_trib": L_trib, "R": R, "center": c,
            "arc_length": float(L_trib.sum()), "chord_length": chord}


def _rect_torsion(a, b):
    """직사각형 비틀림 상수 (긴 변 a, 짧은 변 b 근사식)"""
    lo, hi = min(a, b), max(a, b)
    return hi * lo ** 3 / 3.0 * (1.0 - 0.63 * lo / hi)


def ring_section(D, B, stem_b=RING_STEM_B, tf=RING_FLANGE_T):
    """역T 링 단면. 기둥부 높이 h = D - tf. 성립하지 않으면 None.
    I_horiz: 수평축(연직면 휨, 부재 국부 y축), I_vert: 연직축(링 평면 내 휨, 부재 국부 z축)"""
    h = D - tf
    if h <= 0.0 or B < stem_b:
        return None
    A = stem_b * h + B * tf
    y_f, y_s = tf / 2.0, tf + h / 2.0
    yc = (B * tf * y_f + stem_b * h * y_s) / A
    I_horiz = (B * tf ** 3 / 12.0 + B * tf * (y_f - yc) ** 2
               + stem_b * h ** 3 / 12.0 + stem_b * h * (y_s - yc) ** 2)
    I_vert = h * stem_b ** 3 / 12.0 + tf * B ** 3 / 12.0
    J = _rect_torsion(B, tf) + _rect_torsion(h, stem_b)
    return {"A": A, "I_horiz": I_horiz, "I_vert": I_vert, "J": J, "h": h, "D": D, "B": B}


def ring_ground_per_m(D, B, cover_t, tf=RING_FLANGE_T):
    """링 단위길이(1 m)당 콘크리트 무게, 플랜지 위 흙 무게, 인발 저항, 활동 저항 (N/m).

    인발: 플랜지 위 흙 기둥 무게(양쪽) + 덮개 상재압(바깥쪽만) + 양끝 연직 미끄럼면 마찰
          (안쪽은 돔 내부 바닥이라 덮개 없음. 내압이 바닥을 누르는 효과는 무시 -> 보수적)
    활동: 안쪽 면 수동토압(덮개 없음, 보수적) + (링 + 플랜지 위 흙) 무게 x 바닥 마찰"""
    sec = ring_section(D, B, tf=tf)
    if sec is None:
        return None
    g = planet_g()
    _, _, _, rho_c = RING_MATERIAL_CASES[RING_MATERIAL]
    gam = SITE_REGOLITH_DENSITY * g
    q0 = COVER_DENSITY * g * cover_t
    phi = math.radians(SITE_PHI_DEG)
    h = sec["h"]
    over = (B - RING_STEM_B) / 2.0                       # 한쪽 플랜지 돌출 길이
    conc_w = sec["A"] * rho_c * g
    soil_w = gam * h * 2.0 * over + q0 * over
    friction = SITE_UPLIFT_K * math.tan(phi) * ((gam * h ** 2 / 2.0 + q0 * h) + gam * h ** 2 / 2.0)
    Kp = math.tan(math.pi / 4.0 + phi / 2.0) ** 2
    mu = math.tan(SITE_BASE_FRICTION_RATIO * phi)
    slide = 0.5 * Kp * gam * D ** 2 + mu * (conc_w + soil_w)
    return {"sec": sec, "conc_w": conc_w, "soil_w": soil_w, "uplift": soil_w + friction, "slide": slide}


def bearing_capacity_pa(B, D, phi_deg=SITE_PHI_DEG, c=SITE_COHESION_PA, rho=SITE_REGOLITH_DENSITY):
    """띠 기초 극한 지지력 q_ult = c·Nc + γ·D·Nq + ½·γ·B·Nγ (Terzaghi-Meyerhof, Vesic Nγ). 화성 g."""
    phi = math.radians(phi_deg)
    gam = rho * planet_g()
    Nq = math.exp(math.pi * math.tan(phi)) * math.tan(math.pi / 4.0 + phi / 2.0) ** 2
    Ng = 2.0 * (Nq + 1.0) * math.tan(phi)
    Nc = (Nq - 1.0) / math.tan(phi)
    return c * Nc + gam * D * Nq + 0.5 * gam * B * Ng


def cover_curvature_factor(t, R):
    """반지름 R 구 위에 두께 t로 덮을 때 단위 쉘 면적당 흙 부피 / t"""
    return 1.0 + t / R + t ** 2 / (3.0 * R ** 2)


# ----------------------------------------------------------------------
# 3. 부재 단면(각형강관, SHS(!!square hollow section!!))
# ----------------------------------------------------------------------
# 단면 두께비 t/D: 재료별 국부좌굴 한계 (AISC 360-16 Table B4.1a Case 6: b/t <= 1.40√(E/Fy), b = D - 3t), 소수 셋째 자리 올림
MEMBER_T_D_RATIO = math.ceil(1000.0 / (1.40 * math.sqrt(MEMBER_E / MEMBER_FY) + 3.0)) / 1000.0

# 자동 최소질량 탐색 설정
MASS_SEARCH_START_TON = 200.0          # 시작 질량 (SAFE면 자동으로 아래로 내려가며 FAIL을 찾는다)
MASS_SEARCH_COARSE_STEP_TON = 100.0    # 브라케팅 고정 증감폭
MASS_SEARCH_TOL_TON = 5.0              # 최종 질량 해상도 (ton)
MASS_SEARCH_MAX_TON = 5000.0           # 무한 반복 방지용 탐색 상한
MASS_SEARCH_DIRECTION = 1              # 지진 발생 지방시 계산·영상에 쓰는 기준 성분 (1=X 북남)
DISP_LIMIT_RATIO = None                # 예: 1/400 -> 지진 등가정적 최대 수평변위(기초 이동 포함, 보수적)/스팬 초과 시 FAIL (None이면 보고만)


def total_member_length(nodes, edges):
    """모든 부재의 중심선 길이 합계(m)."""
    return sum(float(np.linalg.norm(np.asarray(nodes[nj], float) - np.asarray(nodes[ni], float)))
               for _, ni, nj in edges)


def size_uniform_shs_for_target_mass(nodes, edges, target_mass_ton,
                                     t_over_D=MEMBER_T_D_RATIO, density=MEMBER_DENSITY):
    """동일 SHS를 쓰는 구조의 총 강재 질량이 목표값이 되도록 D와 t를 계산한다.
    A_target = M_target / (rho * sum(L)),  t = alpha*D 일 때 A = 4*alpha*(1-alpha)*D^2"""
    if target_mass_ton <= 0.0:
        raise ValueError("target_mass_ton은 0보다 커야 합니다.")
    if not 0.0 < t_over_D < 0.5:
        raise ValueError("t_over_D는 0과 0.5 사이여야 합니다.")
    length_sum = total_member_length(nodes, edges)
    if length_sum <= 0.0:
        raise ValueError("총 부재 길이가 0입니다.")
    area = target_mass_ton * 1000.0 / (density * length_sum)
    D = np.sqrt(area / (4.0 * t_over_D * (1.0 - t_over_D)))
    return float(D), float(t_over_D * D), float(area), float(length_sum)


DOME_TOTAL_MEMBER_LENGTH = total_member_length(dome_NODES, DOME_EDGES)


def steel_shs_section(D, t, E=None, G=None):
    """정사각형 각관(SHS) 물성 (모서리 라운드 무시). 재료는 [재료 파트]의 MEMBER_MATERIAL
        A = D^2 - (D-2t)^2,  I = (D^4 - (D-2t)^4)/12,  J = t(D-t)^3 (Bredt),  c = D/2"""
    if E is None:
        E = MEMBER_E
    if G is None:
        G = E / (2 * (1 + MEMBER_NU))
    Do = D - 2 * t
    A = D ** 2 - Do ** 2
    I = (D ** 4 - Do ** 4) / 12.0
    J = t * (D - t) ** 3
    return {"A": A, "E": E, "G": G, "J": J, "Iy": I, "Iz": I, "c": D / 2.0, "D": D, "t": t}


def make_section_for_target_mass(target_mass_ton):
    """목표 강재질량에 맞는 동일 SHS 단면과 단면 산정정보를 반환한다."""
    D, t, area, length_sum = size_uniform_shs_for_target_mass(
        dome_NODES, DOME_EDGES, target_mass_ton=target_mass_ton, t_over_D=MEMBER_T_D_RATIO)
    section = steel_shs_section(D=D, t=t)
    info = {"target_mass_ton": float(target_mass_ton), "D": D, "t": t, "A": area,
            "length_sum": length_sum}
    return section, info


def compute_tributary_mass(nodes, edges, node_tags, A, density=MEMBER_DENSITY):
    """부재 자중(단면적 x 길이 x 밀도)을 양 끝 절점에 절반씩 분배해서 절점질량을 구한다."""
    mass = defaultdict(float)
    for _, ni, nj in edges:
        L = np.linalg.norm(np.asarray(nodes[nj], float) - np.asarray(nodes[ni], float))
        m_ele = A * L * density
        mass[node_tags[ni]] += m_ele / 2
        mass[node_tags[nj]] += m_ele / 2
    return dict(mass)


# ----------------------------------------------------------------------
# 4. 구조물 + 링 기초 + 지반 생성
# ----------------------------------------------------------------------
TRANSF_VERTICAL = NODE_OFFSET + 9001    # 부재축이 global Z와 나란할 때
TRANSF_GENERAL = NODE_OFFSET + 9002     # 그 외 (링 부재 포함: 국부 z = 연직)


def build_structure(section, design, transf_type="PDelta", contact=None, verbose=True):
    """돔 + (링 기초) + 지반 스프링을 만든다.
      contact=None : 비선형 지반 (실제 해석용. ENT + 흙 인발 한계 + 수평 한계)
      contact=배열 : 선형 지반 (단위하중 해석용). 링 노드별 수직 접선강성을 LC1 기준상태로 선택
                     (True = 눌려 있음 -> 지반반력 강성, False = 들려서 플랜지 위 흙이 잡는 중 -> 인발 강성)"""
    if not isinstance(section, dict):
        raise ValueError("section은 steel_shs_section(D, t) 또는 make_section_for_target_mass()의 결과여야 합니다.")
    A, E, G, J, Iy, Iz = (section[k] for k in ("A", "E", "G", "J", "Iy", "Iz"))

    ops.wipe()
    _LINE_LOAD_STATE.clear()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.geomTransf(transf_type, TRANSF_VERTICAL, 1, 0, 0)
    ops.geomTransf(transf_type, TRANSF_GENERAL, 0, 0, 1)

    node_tags = {}
    for tag, (x, y, z) in dome_NODES.items():
        ops.node(NODE_OFFSET + tag, x, y, z)
        node_tags[tag] = NODE_OFFSET + tag
    base_tags = get_base_node_tags(dome_NODES)
    base_node_tags = {tag: node_tags[tag] for tag in base_tags}

    elements = {}
    for ele_tag, ni, nj in DOME_EDGES:
        xi, yi, _ = dome_NODES[ni]
        xj, yj, _ = dome_NODES[nj]
        is_vertical = abs(xi - xj) < 1e-9 and abs(yi - yj) < 1e-9
        ops.element("elasticBeamColumn", ele_tag, node_tags[ni], node_tags[nj],
                    A, E, G, J, Iy, Iz, TRANSF_VERTICAL if is_vertical else TRANSF_GENERAL)
        elements[ele_tag] = ele_tag

    mass_dict = compute_tributary_mass(dome_NODES, DOME_EDGES, node_tags, A)
    steel_mass = A * DOME_TOTAL_MEMBER_LENGTH * MEMBER_DENSITY
    model = {"node_tags": node_tags, "base_node_tags": base_node_tags, "elements": elements,
             "mass_dict": mass_dict, "section": (A, E, G, J, Iy, Iz), "c": section["c"],
             "steel_mass_kg": steel_mass}

    if FOUNDATION_MODE == "ring":
        model["ground"] = _build_ring_and_ground(model, design, contact)
    else:
        model["ground"] = _build_bare_ground(model, section, contact)

    if verbose:
        gr = model["ground"]
        extra = ((f", 링 {len(gr['ring_eles'])}구간 (D={design.ring_depth:.1f} m, B={design.ring_b:.1f} m, "
                  f"{gr['ring_mass_kg']/1000:,.0f} ton)") if gr["mode"] == "ring" else
                 ", 기초 없음 (밑단 절점이 지반 위에 놓임)")
        print(f"[model] 돔 절점 {len(node_tags)}개, 부재 {DOME_N_MEMBERS}개 (요소 {len(elements)}개, SHS {MEMBER_MATERIAL} {steel_mass/1000:.1f} ton)"
              f"{extra}, 지반 지지점 {len(gr['fixed'])}개 ({'선형' if contact is not None else '비선형'})")
    return model


def _build_ring_and_ground(model, design, contact):
    lay = RING_LAYOUT
    gp = ring_ground_per_m(design.ring_depth, design.ring_b, design.cover_t, design.flange_t)
    if gp is None:
        raise ValueError(f"링 단면이 성립하지 않습니다: {design.label()} (D > 바닥판 두께, B >= 기둥부 폭)")
    sec = gp["sec"]
    fc, ft, Ec, rho_c = RING_MATERIAL_CASES[RING_MATERIAL]
    Gc = Ec / (2.0 * (1.0 + RING_NU))
    B, L = design.ring_b, lay["L_trib"]

    ring_nodes = []
    for i, (x, y, bt) in enumerate(lay["points"]):
        if bt is None:
            tag = RING_NODE_BASE + i + 1
            ops.node(tag, x, y, 0.0)
        else:
            tag = model["node_tags"][bt]              # 돔 밑단 노드 = 링 노드 (강접합)
        ring_nodes.append(tag)

    n = len(ring_nodes)
    ring_eles = []
    for i in range(n):
        e = RING_ELE_BASE + i + 1
        ops.element("elasticBeamColumn", e, ring_nodes[i], ring_nodes[(i + 1) % n],
                    sec["A"], Ec, Gc, sec["J"], sec["I_horiz"], sec["I_vert"], TRANSF_GENERAL)
        ring_eles.append(e)

    ring_mass = 0.0
    for i, tag in enumerate(ring_nodes):
        m = sec["A"] * rho_c * L[i]
        model["mass_dict"][tag] = model["mass_dict"].get(tag, 0.0) + m
        ring_mass += m

    fixed, springs = [], []
    cap_up, cap_h, soil_w, kz_all = (np.zeros(n) for _ in range(4))
    for i, (tag, (x, y, _)) in enumerate(zip(ring_nodes, lay["points"])):
        kz = SITE_SUBGRADE_K * B * L[i]
        kh = SITE_H_TO_V_STIFFNESS * kz
        k_rock = kz * B ** 2 / 12.0
        cap_up[i] = gp["uplift"] * L[i]
        cap_h[i] = gp["slide"] * L[i]
        soil_w[i] = gp["soil_w"] * L[i]
        kz_all[i] = kz
        k_up = cap_up[i] / UPLIFT_MOBILIZE_M

        mb = GROUND_MAT_BASE + 10 * i
        if contact is None:
            ops.uniaxialMaterial("ElasticPP", mb + 1, kh, cap_h[i] / kh)          # 수평: 한계에서 미끄러짐
            ops.uniaxialMaterial("ENT", mb + 7, kz)                               # 수직 눌림: 링 밑 흙
            ops.uniaxialMaterial("ElasticPPGap", mb + 8, k_up, cap_up[i], 0.0)    # 수직 들림: 플랜지 위 흙
            ops.uniaxialMaterial("Parallel", mb + 3, mb + 7, mb + 8)
        else:
            ops.uniaxialMaterial("Elastic", mb + 1, kh)
            ops.uniaxialMaterial("Elastic", mb + 3, kz if contact[i] else k_up)
        ops.uniaxialMaterial("Elastic", mb + 4, k_rock)
        ops.uniaxialMaterial("Elastic", mb + 6, 0.1 * k_rock)

        f_tag = GROUND_NODE_BASE + i + 1
        ops.node(f_tag, x, y, 0.0)
        ops.fix(f_tag, 1, 1, 1, 1, 1, 1)
        s_tag = GROUND_ELE_BASE + i + 1
        ops.element("zeroLength", s_tag, f_tag, tag,
                    "-mat", mb + 1, mb + 1, mb + 3, mb + 4, mb + 4, mb + 6, "-dir", 1, 2, 3, 4, 5, 6)
        fixed.append(f_tag)
        springs.append(s_tag)

    return {"mode": "ring", "fixed": fixed, "springs": springs, "ring_nodes": ring_nodes,
            "ring_eles": ring_eles, "ring_sec": sec, "ring_mass_kg": ring_mass,
            "ring_fc": fc, "ring_ft": ft, "L_trib": L.copy(), "B": B,
            "cap_up": cap_up, "cap_h": cap_h, "soil_w": soil_w, "kz": kz_all,
            "q_allow": bearing_capacity_pa(B, design.ring_depth) / FS_BEARING}


def _build_bare_ground(model, section, contact):
    """기초 없음: 돔 밑단 절점이 원지반 위에 그냥 놓인다. 링·흙·마찰 저항 없음.
    지표에 놓인 밑단 둘레 부재가 폭 = 부재 바깥 폭 D, 길이 = 절점이 맡는 둘레 길이만큼 흙을 누른다고 보고
    링 기초와 같은 지반 값(SITE_SUBGRADE_K, SITE_H_TO_V_STIFFNESS)으로 강성을 정한다.
    수직은 눌림만 (ENT, 인장 0), 수평·회전은 탄성. contact=배열이면 선형 지반 (눌림 = kz, 들림 = 거의 0)"""
    b = section["D"]
    L_of = {bt: L for (_, _, bt), L in zip(BARE_LAYOUT["points"], BARE_LAYOUT["L_trib"])}
    fixed, springs = [], []
    for i, (tag, stag) in enumerate(model["base_node_tags"].items()):
        x, y, _ = dome_NODES[tag]
        kz = SITE_SUBGRADE_K * b * L_of[tag]
        kh = SITE_H_TO_V_STIFFNESS * kz
        k_rock = kz * b ** 2 / 12.0
        mb = GROUND_MAT_BASE + 10 * i
        ops.uniaxialMaterial("Elastic", mb + 1, kh)                                # 수평
        if contact is None:
            ops.uniaxialMaterial("ENT", mb + 3, kz)                                # 수직: 눌림만
        else:
            ops.uniaxialMaterial("Elastic", mb + 3, kz if contact[i] else 1e-6 * kz)
        ops.uniaxialMaterial("Elastic", mb + 4, k_rock)
        ops.uniaxialMaterial("Elastic", mb + 6, 0.1 * k_rock)
        f_tag = GROUND_NODE_BASE + i + 1
        ops.node(f_tag, x, y, 0.0)
        ops.fix(f_tag, 1, 1, 1, 1, 1, 1)
        s_tag = GROUND_ELE_BASE + i + 1
        ops.element("zeroLength", s_tag, f_tag, stag,
                    "-mat", mb + 1, mb + 1, mb + 3, mb + 4, mb + 4, mb + 6, "-dir", 1, 2, 3, 4, 5, 6)
        fixed.append(f_tag)
        springs.append(s_tag)
    return {"mode": "none", "fixed": fixed, "springs": springs, "ring_nodes": [], "ring_eles": [],
            "ring_mass_kg": 0.0}


# ----------------------------------------------------------------------
# 5. 하중 -> 절점력 변환 (모두 {원래 nodeTag: np.array([Fx, Fy, Fz])} 형태로 반환)
# ----------------------------------------------------------------------
def _oriented_shell_panels(nodes, edges):
    """외향 정렬된 쉘 조각 [(분배 절점들, 분배 비율들, 외향 면적벡터, 중심높이 z), ...]"""
    panels = []
    for ids, w, (v0, v1, v2) in _shell_triangles(nodes, edges):
        area_vec = 0.5 * np.cross(v1 - v0, v2 - v0)
        if np.linalg.norm(area_vec) < 1e-12:
            continue
        panels.append((ids, w, area_vec, float((v0[2] + v1[2] + v2[2]) / 3.0)))
    return panels


def _sum_forces(*force_dicts):
    total = defaultdict(lambda: np.zeros(3))
    for fd in force_dicts:
        for k, f in fd.items():
            total[k] = total[k] + np.asarray(f, dtype=float)
    return dict(total)


def compute_membrane_node_mass(node_tags, nodes, edges, areal_density=MEMBRANE_AREAL_DENSITY, verbose=True):
    """막 자중(면적 x 면밀도)을 패널 절점에 분배비율대로 (지오데식 1/3씩). 키는 구조물 절점 태그."""
    mass = defaultdict(float)
    for ids, w, area_vec, _ in _oriented_shell_panels(nodes, edges):
        m = areal_density * float(np.linalg.norm(area_vec))
        for k, wk in zip(ids, w):
            mass[node_tags[k]] += m * wk
    if verbose:
        print(f"[membrane] {MEMBRANE_NAME} (면밀도 {areal_density:.2f} kg/m^2) 총 질량 {sum(mass.values())/1000:.2f} ton")
    return dict(mass)


def compute_cover_node_mass(node_tags, nodes, edges, t, rho=COVER_DENSITY, verbose=True):
    """레골리스 덮개: 패널 면적 x t x 곡률보정 x rho 를 세 절점에 1/3씩. 하중은 apply_gravity가 수직 하향으로 건다."""
    if t <= 0.0:
        return {}
    _, R = _base_center_radius(nodes)
    f = cover_curvature_factor(t, R)
    mass = defaultdict(float)
    for ids, w, area_vec, _ in _oriented_shell_panels(nodes, edges):
        m = rho * t * f * float(np.linalg.norm(area_vec))
        for k, wk in zip(ids, w):
            mass[node_tags[k]] += m * wk
    if verbose:
        tot = sum(mass.values())
        print(f"[cover] 덮개 t={t:.2f} m, rho={rho:.0f} kg/m^3, 곡률보정 {f:.3f} -> "
              f"{tot/1000:,.0f} ton ({tot*planet_g()/1e6:,.1f} MN)")
    return dict(mass)


def pressure_node_forces(nodes, edges, p_net):
    """순압 p_net(+면 바깥으로 미는 힘)을 패널 외향 방향으로 걸고 세 절점에 1/3씩 분배."""
    force = defaultdict(lambda: np.zeros(3))
    for ids, w, area_vec, _ in _oriented_shell_panels(nodes, edges):
        for k, wk in zip(ids, w):
            force[k] = force[k] + p_net * area_vec * wk
    return dict(force)


def wind_height_factor(z):
    """로그 법칙 풍속 분포에 따른 풍압 배율 (V(z)/V(z_ref))^2. 측정 높이 아래는 측정값 그대로."""
    if not USE_WIND_HEIGHT_PROFILE:
        return 1.0
    z_eff = max(z, WIND_REF_HEIGHT_M)
    return (math.log(z_eff / WIND_ROUGHNESS_M) / math.log(WIND_REF_HEIGHT_M / WIND_ROUGHNESS_M)) ** 2


def wind_node_forces(nodes, edges, q_ref, az_to_deg):
    """측정높이 풍압 q_ref(Pa), 바람이 불어가는 방위 az_to_deg(북=0°, 동=90°)에 대한 절점력.
    모델 좌표는 지진 채널과 같게 X=북, Y=동 이다."""
    az = math.radians(az_to_deg)
    wdir = np.array([math.cos(az), math.sin(az), 0.0])
    force = defaultdict(lambda: np.zeros(3))
    for ids, w, area_vec, zc in _oriented_shell_panels(nodes, edges):
        n_hat = area_vec / np.linalg.norm(area_vec)
        cos_t = float(np.dot(n_hat, wdir))
        cp = -WIND_CP_WINDWARD * cos_t if cos_t < 0.0 else WIND_CP_SUCTION
        f_tri = -q_ref * wind_height_factor(zc) * cp * area_vec    # 압력(+)은 안쪽, 흡입(-)은 바깥쪽
        for k, wk in zip(ids, w):
            force[k] = force[k] + f_tri * wk
    return dict(force)


def thermal_node_forces(nodes, edges, axial_force):
    """균일 온도변화 ΔT의 등가 절점력. axial_force = E·A·α·ΔT.
    온도가 오르면 부재가 늘어나려 하므로 양 끝 절점을 부재축 방향 바깥으로 민다.
    실제 부재력 = 해석 부재력 - E·A·α·ΔT (correct_thermal_forces에서 보정). 링은 묻혀 있어 ΔT = 0."""
    force = defaultdict(lambda: np.zeros(3))
    for _, ni, nj in edges:
        p1, p2 = np.asarray(nodes[ni], float), np.asarray(nodes[nj], float)
        ex = (p2 - p1) / np.linalg.norm(p2 - p1)
        force[ni] = force[ni] - axial_force * ex
        force[nj] = force[nj] + axial_force * ex
    return dict(force)


# ----------------------------------------------------------------------
# 5-0. 외피 하중 -> 부재 선하중
#   외피(구속층)는 부재 길이 전체에 고정되어 있어, 패널 하중(내압, 덮개, 바람)이 부재에 등분포하중으로 전달된다.
#   패널 하중 합력 F를 패널 둘레에 고르게 나눠 (둘레 1 m당 F/둘레, 삼각형이면 p·a/2) 각 변의 부재에 건다.
#   - 막 장력의 면내 성분(패널 안쪽으로 당기는 힘)은 이웃 패널끼리 상쇄된다고 본다.
#   - 밑단 변(양 끝이 지표 z = 0)은 외피를 링 기초에 직접 고정한다고 보고, 그 몫은 밑단 절점(링)에 절점력으로 준다.
#   - 부재 중앙 단면도 검토하도록 부재력 배열은 [I단 6 | J단 6 | 중앙 6] 18열이다 (element_local_forces).
# ----------------------------------------------------------------------
MEMBER_SECTIONS = ("I", "J", "M")         # 부재 검토 단면: I단, J단, 중앙
SECTION_LABEL = {"I": "I단", "J": "J단", "M": "중앙"}
_LINE_LOAD_STATE = {}                     # 지금 모델에 걸린 선하중 {요소: 국부 (wx, wy, wz) N/m}, 중앙 단면 계산용


def _element_axes():
    """부재 국부축 (OpenSees geomTransf와 같은 규칙: y = vecxz x x, z = x x y)과 길이"""
    axes, lengths = {}, {}
    for e, ni, nj in DOME_EDGES:
        p1, p2 = np.asarray(dome_NODES[ni], float), np.asarray(dome_NODES[nj], float)
        L = float(np.linalg.norm(p2 - p1))
        x = (p2 - p1) / L
        vertical = abs(p1[0] - p2[0]) < 1e-9 and abs(p1[1] - p2[1]) < 1e-9
        vecxz = np.array([1.0, 0.0, 0.0]) if vertical else np.array([0.0, 0.0, 1.0])
        y = np.cross(vecxz, x)
        y /= np.linalg.norm(y)
        axes[e], lengths[e] = (x, y, np.cross(x, y)), L
    return axes, lengths


ELE_AXES, ELE_LENGTH = _element_axes()


def _build_shell_loops():
    """외피 패널(지오데식 삼각형, 곡선 리브 돔은 리브로 둘러싸인 곡면 패널)마다 외향 면적벡터, 면적, 둘레, 변 목록.
    변: (절점 a, 절점 b, 요소 태그, 길이, 밑단 변 여부)"""
    P = {k: np.asarray(v, dtype=float) for k, v in dome_NODES.items()}
    c, R_sph = _base_center_radius(dome_NODES)
    loops = DOME_PANELS if DOME_PANELS is not None else [list(t) for t in _find_triangles_from_edges(dome_NODES, DOME_EDGES)]
    ele_of = {frozenset((i, j)): e for e, i, j in DOME_EDGES}
    out = []
    for lp in loops:
        lp = list(lp)
        g = np.mean([P[k] for k in lp], axis=0)
        if len(lp) > 3:
            g = c + R_sph * (g - c) / np.linalg.norm(g - c)
        n = len(lp)
        fan = [0.5 * np.cross(P[lp[i]] - g, P[lp[(i + 1) % n]] - g) for i in range(n)]
        avec = np.sum(fan, axis=0)
        if np.dot(avec, g - c) < 0:                          # 외향 정렬
            lp, avec = lp[::-1], -avec
        edges = []
        for i in range(n):
            a, b = lp[i], lp[(i + 1) % n]
            edges.append((a, b, ele_of.get(frozenset((a, b))), float(np.linalg.norm(P[b] - P[a])),
                          abs(P[a][2]) < 1e-6 and abs(P[b][2]) < 1e-6))
        out.append({"nodes": lp, "area_vec": avec, "area": float(sum(np.linalg.norm(f) for f in fan)),
                    "perimeter": sum(e[3] for e in edges), "zc": float(g[2]), "edges": edges})
    missing = [(a, b) for lp in out for a, b, e, _, base in lp["edges"] if e is None and not base]
    if missing:
        raise ValueError(f"패널 변 중 부재가 없는 변이 있습니다: {missing[:5]}")
    return out


SHELL_LOOPS = _build_shell_loops()


def pressure_panel_forces(p_net):
    """순압 p_net(+ 바깥으로 미는 힘) -> 패널별 합력 [(패널, 전역 힘벡터)]"""
    return [(lp, p_net * lp["area_vec"]) for lp in SHELL_LOOPS]


def cover_panel_forces(t, rho=COVER_DENSITY):
    """레골리스 덮개 무게 (연직 하향) -> 패널별 합력. compute_cover_node_mass와 같은 질량"""
    if t <= 0.0:
        return []
    w = rho * planet_g() * t * cover_curvature_factor(t, DOME_SPHERE_RADIUS_M)
    return [(lp, np.array([0.0, 0.0, -w * lp["area"]])) for lp in SHELL_LOOPS]


def wind_panel_forces(q_ref, az_to_deg):
    """측정높이 풍압 q_ref, 바람이 불어가는 방위 az_to_deg -> 패널별 합력 (wind_node_forces와 같은 풍압계수)"""
    az = math.radians(az_to_deg)
    wdir = np.array([math.cos(az), math.sin(az), 0.0])
    out = []
    for lp in SHELL_LOOPS:
        cos_t = float(np.dot(lp["area_vec"] / np.linalg.norm(lp["area_vec"]), wdir))
        cp = -WIND_CP_WINDWARD * cos_t if cos_t < 0.0 else WIND_CP_SUCTION
        out.append((lp, -q_ref * wind_height_factor(lp["zc"]) * cp * lp["area_vec"]))
    return out


def panel_line_loads(panel_forces):
    """패널별 합력 -> 부재 선하중 {요소: 전역 (wx, wy, wz) N/m} + 밑단 변 몫의 절점력 {원래 절점: (Fx, Fy, Fz)}"""
    line, nodal = defaultdict(lambda: np.zeros(3)), defaultdict(lambda: np.zeros(3))
    for lp, F in panel_forces:
        q = np.asarray(F, dtype=float) / lp["perimeter"]
        for a, b, ele, L, base in lp["edges"]:
            if base:
                nodal[a] = nodal[a] + q * L / 2.0
                nodal[b] = nodal[b] + q * L / 2.0
            else:
                line[ele] = line[ele] + q
    return dict(line), dict(nodal)


def apply_shell_loads(model, panel_forces, legacy_nodal, pattern_tag, ts_tag=1):
    """외피 하중을 새 하중패턴으로 건다.
    USE_LINE_LOADS: 부재 선하중(eleLoad -beamUniform) + 밑단 변 절점력. 아니면 legacy_nodal(이전 꼭짓점 절점력)."""
    ops.pattern("Plain", pattern_tag, ts_tag)
    if not USE_LINE_LOADS:
        for k, f in legacy_nodal.items():
            ops.load(model["node_tags"][k], float(f[0]), float(f[1]), float(f[2]), 0.0, 0.0, 0.0)
        return
    line, nodal = panel_line_loads(panel_forces)
    for e, wg in line.items():
        x, y, z = ELE_AXES[e]
        wl = np.array([wg @ x, wg @ y, wg @ z])
        ops.eleLoad("-ele", e, "-type", "-beamUniform", float(wl[1]), float(wl[2]), float(wl[0]))
        _LINE_LOAD_STATE[e] = _LINE_LOAD_STATE.get(e, np.zeros(3)) + wl
    for k, f in nodal.items():
        ops.load(model["node_tags"][k], float(f[0]), float(f[1]), float(f[2]), 0.0, 0.0, 0.0)


def correct_thermal_forces(local_forces, axial_force):
    """localForce [FxI ... | FxJ ... | 중앙 N ...]에 열하중 고정단력 보정. N_실제 = N_해석 - EAαΔT
    (N = -FxI = FxJ = 중앙 N 이므로 FxI += EAαΔT, FxJ -= EAαΔT, 중앙 N -= EAαΔT)"""
    f = np.array(local_forces, dtype=float, copy=True)
    f[..., 0] += axial_force
    f[..., 6] -= axial_force
    if f.shape[-1] > 12:
        f[..., 12] -= axial_force
    return f


def apply_nodal_forces(node_tags, forces, pattern_tag, ts_tag=1):
    ops.pattern("Plain", pattern_tag, ts_tag)
    for tag, f in forces.items():
        ops.load(node_tags[tag], float(f[0]), float(f[1]), float(f[2]), 0.0, 0.0, 0.0)


def element_local_forces(ele_tags=ELE_TAGS):
    """(부재수, 18) = [FxI FyI FzI MxI MyI MzI | FxJ FyJ FzJ MxJ MyJ MzJ | 중앙 N Vy Vz T My Mz]
    중앙 단면은 J단과 같은 부호(N > 0 인장)로, I단 힘과 지금 걸린 선하중(국부 wx, wy, wz)에서 평형으로 구한다."""
    out = np.zeros((len(ele_tags), 18))
    zero = np.zeros(3)
    for r, e in enumerate(ele_tags):
        F = np.asarray(ops.eleResponse(e, "localForce"), dtype=float)
        wx, wy, wz = _LINE_LOAD_STATE.get(e, zero)
        h = ELE_LENGTH[e] / 2.0
        out[r, :12] = F
        out[r, 12:] = (-(F[0] + wx * h), -(F[1] + wy * h), -(F[2] + wz * h), -F[3],
                       -(F[4] + h * F[2] + wz * h * h / 2.0), -(F[5] - h * F[1] - wy * h * h / 2.0))
    return out


def ground_reactions(model):
    """지반 지지점 반력 (n, 3). +z = 지반이 링을 밀어 올림(눌림), -z = 플랜지 위 흙이 잡아당김(들림)."""
    ops.reactions()
    return np.array([ops.nodeReaction(f)[:3] for f in model["ground"]["fixed"]], dtype=float)


def capture_node_disp(node_tags):
    """돔 절점 변위 (절점수, 3)"""
    return np.array([ops.nodeDisp(s)[:3] for s in node_tags.values()], dtype=float)


def ring_axial(model):
    """링 부재 축력 (m,). +인장"""
    eles = model["ground"]["ring_eles"]
    return np.array([ops.eleResponse(e, "localForce")[6] for e in eles], dtype=float)


# 형상에서 한 번만 계산하는 상수
RING_LAYOUT = ring_layout(dome_NODES, n_sub=max(1, RING_TOTAL_NODES // len(get_base_node_tags(dome_NODES))))
BARE_LAYOUT = ring_layout(dome_NODES, n_sub=1)   # 기초 없음: 돔 밑단 절점만, 절점별 둘레 길이
DOME_SHELL_AREA_M2 = float(sum(np.linalg.norm(p[2]) for p in _oriented_shell_panels(dome_NODES, DOME_EDGES)))
DOME_UPLIFT_PER_PA = float(sum(f[2] for f in pressure_node_forces(dome_NODES, DOME_EDGES, 1.0).values()))   # m^2
DOME_SPHERE_RADIUS_M = _base_center_radius(dome_NODES)[1]


# ----------------------------------------------------------------------
# 5-1. 외피 강도 검토 — 패널마다 구속층(Kevlar) 장력
#   부푼 패널 = 반지름 R인 구면 조각 -> 막 장력 N = p·R/2 (풍선 식, 단위폭당)
#   R은 패널 내접원 반지름 a와 부푼 높이 f = 2a x PANEL_SAG_RATIO 로 정한다: R = (a² + f²) / 2f
# ----------------------------------------------------------------------
def envelope_panels(nodes=None, edges=None):
    """패널(지오데식 삼각형, 곡선 리브 돔은 리브로 둘러싸인 곡면 패널)마다 면적, 둘레, 내접원 반지름 a, 부푼 반지름 R.
    a = 2 x 면적 / 둘레 (삼각형은 정확, 다각형 패널은 근사)."""
    nodes = dome_NODES if nodes is None else nodes
    edges = DOME_EDGES if edges is None else edges
    P = {k: np.asarray(v, dtype=float) for k, v in nodes.items()}
    c, R_sph = _base_center_radius(nodes)
    loops = DOME_PANELS if DOME_PANELS is not None else [list(t) for t in _find_triangles_from_edges(nodes, edges)]
    out = []
    for lp in loops:
        V = [P[k] for k in lp]
        g = np.mean(V, axis=0)
        if len(lp) > 3:
            g = c + R_sph * (g - c) / np.linalg.norm(g - c)
        n = len(V)
        area = sum(0.5 * float(np.linalg.norm(np.cross(V[i] - g, V[(i + 1) % n] - g))) for i in range(n))
        perim = sum(float(np.linalg.norm(V[(i + 1) % n] - V[i])) for i in range(n))
        a = 2.0 * area / perim
        f = 2.0 * a * PANEL_SAG_RATIO
        out.append({"nodes": tuple(lp), "area": area, "perimeter": perim, "a": a,
                    "R": (a * a + f * f) / (2.0 * f), "z": float(g[2])})
    return out


ENVELOPE_PANELS = envelope_panels()


def check_envelope(dp_max, cover_t=0.0, verbose=True):
    """외피 강도 검토.
      가압 (LC1): N = Δp_max·R/2. 덮개가 눌러 주는 효과는 무시 (보수적)
      감압 (LC2): N = q·R/2, q = ρ·g·t·곡률보정. 흙 무게를 전부 막에 수직으로 건다 (보수적).
                  구속층은 인장만 받으므로 패널이 안쪽으로 처지며 해먹처럼 버틴다 (대변형은 범위 밖)
      허용 장력: 웨빙 가닥/m x 겹 x 정격 파단하중 / 안전율
      기밀층은 구속층보다 크게 재단해 하중을 받지 않으므로 검토하지 않는다 (TransHab)"""
    R = np.array([p["R"] for p in ENVELOPE_PANELS])
    q = (COVER_DENSITY * planet_g() * cover_t * cover_curvature_factor(cover_t, DOME_SPHERE_RADIUS_M)
         if cover_t > 0 else 0.0)
    N1, N2 = dp_max * R / 2.0, q * R / 2.0
    i1, i2 = int(np.argmax(N1)), int(np.argmax(N2))
    gap = 1.0 / RESTRAINT_STRAPS_PER_M - RESTRAINT_STRAP_WIDTH
    per_strap = RESTRAINT_STRAP_BREAK_N / RESTRAINT_SF
    need = math.ceil(max(N1[i1], N2[i2]) / per_strap)                     # 필요 가닥/m (모든 겹 합)
    out = {"dp_max_Pa": dp_max, "cover_t_m": cover_t, "q_cover_Pa": q,
           "N_lc1_max_N_m": float(N1[i1]), "N_lc2_max_N_m": float(N2[i2]), "N_allow_N_m": RESTRAINT_ALLOW_N_M,
           "dcr_lc1": float(N1[i1] / RESTRAINT_ALLOW_N_M), "dcr_lc2": float(N2[i2] / RESTRAINT_ALLOW_N_M),
           "gov_panel_lc1": ENVELOPE_PANELS[i1]["nodes"], "a_gov_m": ENVELOPE_PANELS[i1]["a"],
           "R_gov_m": float(R[i1]), "a_range_m": (float(min(p["a"] for p in ENVELOPE_PANELS)),
                                                  float(max(p["a"] for p in ENVELOPE_PANELS))),
           "straps_needed_per_m": need, "gap_mm": gap * 1e3,
           "lc2_over_lc1": float(q / dp_max) if dp_max > 0 else np.inf}
    out["dcr"] = max(out["dcr_lc1"], out["dcr_lc2"])
    out["ok"] = out["dcr"] < 1.0
    if verbose:
        print_envelope_check(out)
    return out


def print_envelope_check(ev):
    a0, a1 = ev["a_range_m"]
    print(f"[envelope] {MEMBRANE_NAME}, 면밀도 {MEMBRANE_AREAL_DENSITY:.2f} kg/m^2 "
          f"(기밀층 {BLADDER_AREAL:.2f} + 구속층 {RESTRAINT_AREAL:.2f}), "
          f"돔 전체 {MEMBRANE_AREAL_DENSITY * DOME_SHELL_AREA_M2 / 1000:.1f} t")
    print(f"[envelope] 패널 {len(ENVELOPE_PANELS)}개, 내접원 반지름 {a0:.2f}~{a1:.2f} m, 부푼 높이/폭 {PANEL_SAG_RATIO:.0%} "
          f"-> 최대 R {ev['R_gov_m']:.2f} m | 구속층 허용 장력 {ev['N_allow_N_m']/1e3:.0f} kN/m "
          f"(파단 {RESTRAINT_STRAP_BREAK_N/1e3:.1f} kN x {RESTRAINT_STRAPS_PER_M}가닥 x {RESTRAINT_PLIES}겹 "
          f"/ 안전율 {RESTRAINT_SF})")
    print(f"[envelope] 가압 Δp {ev['dp_max_Pa']/1e3:.1f} kPa: 장력 {ev['N_lc1_max_N_m']/1e3:.0f} kN/m, DCR {ev['dcr_lc1']:.3f}"
          + (f" | 감압 덮개 {ev['cover_t_m']:.1f} m ({ev['q_cover_Pa']/1e3:.1f} kPa): 장력 {ev['N_lc2_max_N_m']/1e3:.0f} kN/m, "
             f"DCR {ev['dcr_lc2']:.3f}" if ev["cover_t_m"] > 0 else ""))
    print(f"[envelope] 웨빙 틈 {ev['gap_mm']:.1f} mm | 필요 웨빙 {ev['straps_needed_per_m']}가닥/m "
          f"(현재 {RESTRAINT_STRAPS_PER_M * RESTRAINT_PLIES}) -> {'OK' if ev['ok'] else 'FAIL'}")


def cover_volume_m3(t):
    return DOME_SHELL_AREA_M2 * t * cover_curvature_factor(t, DOME_SPHERE_RADIUS_M) if t > 0 else 0.0


def effective_env(env, design):
    """덮개가 있으면 부재 ΔT를 exp(-t/δ)로 줄이고 바람을 끈다 (흙이 단열·차폐)."""
    e = dict(env)
    covered = design.cover_t > 0.0
    factor = math.exp(-design.cover_t / DIURNAL_SKIN_DEPTH_M) if covered else 1.0
    e["dT"] = np.asarray(env["dT"]) * factor
    e["thermal_factor"] = factor
    e["wind_active"] = bool(USE_WIND and not covered)
    if not e["wind_active"]:
        e["q_Pa"] = np.zeros_like(np.asarray(env["q_Pa"], dtype=float))
    return e


# ----------------------------------------------------------------------
# 5-2. 칸(D, B) 하나의 기초 설계 — 덮개 t, 바닥판 두께를 돔 질량(가정 없이 실제 값)으로 정한다
# ----------------------------------------------------------------------
def foundation_loads(D, B, t, tf, dome_mass_ton):
    """링 1개 전체에 대한 정적 평형 값 (N, Pa). LC1 = 외기압 최저 시각 가압, LC2 = 감압 사고(내압 0)."""
    g = planet_g()
    Lr = RING_LAYOUT["arc_length"]
    gp = ring_ground_per_m(D, B, t, tf)
    if gp is None:
        return None
    W = {"cover": COVER_DENSITY * cover_volume_m3(t) * g,
         "membrane": MEMBRANE_AREAL_DENSITY * DOME_SHELL_AREA_M2 * g if USE_MEMBRANE else 0.0,
         "ring": gp["conc_w"] * Lr, "soil": gp["soil_w"] * Lr, "dome": dome_mass_ton * 1000.0 * g}
    return {"gp": gp, "W": W,
            "uplift_resist": W["cover"] + W["membrane"] + W["ring"] + gp["uplift"] * Lr + W["dome"],
            "q_lc2": (W["cover"] + W["membrane"] + W["ring"] + W["soil"] + W["dome"]) / (B * Lr),
            "q_allow": bearing_capacity_pa(B, D) / FS_BEARING}


def flange_bending(D, B, t, tf, q_lc2, sig_hoop=0.0):
    """바닥판 = 기둥부 면에서 길이 L = (B - 기둥부 폭)/2 로 튀어나온 외팔보 (링 1 m당).
      LC2: 지반이 밑에서 q로 밀어 올림 - 위 흙·바닥판 자중. 안쪽 바닥판은 덮개가 없어 더 불리
      LC1: 링이 들릴 때 바닥판 위 흙 + 덮개 + 자중이 위에서 누름 (바깥쪽)
      sig_hoop: 해석에서 나온 링 축인장응력 (Pa). 같은 무근 단면에 겹치므로 휨 허용응력에서 뺀다
                -> FS x (6M/h^2 + sig_hoop) <= f_t
    반환: 설계 휨모멘트 M (N·m/m, 두 경우 중 큰 값), 무근으로 버틸 최소 두께 h_plain"""
    g = planet_g()
    _, ft, _, rho_c = RING_MATERIAL_CASES[RING_MATERIAL]
    h = D - tf
    L = max((B - RING_STEM_B) / 2.0, 0.0)
    soil, cover, self_w = SITE_REGOLITH_DENSITY * g * h, COVER_DENSITY * g * t, rho_c * g * tf
    m_lc2 = max(q_lc2 - soil - self_w, 0.0) * L ** 2 / 2.0
    m_lc1 = (soil + cover + self_w) * L ** 2 / 2.0
    M = max(m_lc2, m_lc1)
    allow = ft / FS_FLANGE - sig_hoop                  # 휨에 쓸 수 있는 인장응력 (링 인장 몫을 뺌)
    h_plain = math.sqrt(6.0 * M / allow) if allow > 0.0 else math.inf
    return {"M": M, "M_lc2": m_lc2, "M_lc1": m_lc1, "L": L, "h_plain": h_plain, "sig_hoop": sig_hoop}


def ring_hoop_stress(design, case):
    """해석 결과의 링 최대 축인장 / 링 단면적 (Pa). 바닥판 휨 인장과 같은 단면에 겹치는 인장."""
    if FOUNDATION_MODE != "ring" or not case or "ring" not in case:
        return 0.0
    sec = ring_section(design.ring_depth, design.ring_b, tf=design.flange_t)
    return float(case["ring"]["T_max"]) / sec["A"] if sec else 0.0


def flange_rebar_kg(D, B, tf, M_lc2, M_lc1):
    """철근 바닥판: 위쪽 면(LC2) + 아래쪽 면(LC1) 철근을 바닥판 전체 폭 B에 걸쳐 링 둘레만큼"""
    d = tf - FLANGE_REBAR_COVER
    if d <= 0:
        return np.inf
    As = FS_FLANGE * (M_lc2 + M_lc1) / (0.9 * d * RING_REBAR_FY)          # m^2 / (링 1 m)
    return As * B * RING_LAYOUT["arc_length"] * REBAR_DENSITY


def design_cell(D, B, env, dome_mass_ton=0.0, t_floor=None, sig_hoop=0.0):
    """칸(D, B) 하나의 기초 설계. 돔 질량 dome_mass_ton은 가정값이 아니라 이전 바퀴에서 찾은 실제 최소 질량(처음엔 0).
      ① 들림: 덮개를 방사선 하한부터 올리며 (덮개+막+링+바닥판 위 흙·마찰+돔) >= FS x 들림 인 최소 t
      ② 지지력(감압 사고): 바닥 압력 <= 극한지지력 / FS. 못 넘으면 이 칸 탈락 (B가 좁음)
      ③ 바닥판 휨: 필요한 두께를 구하고, 두께가 바뀌면 ①부터 다시 (링 무게·흙 높이가 바뀌므로)
         sig_hoop = 이전 바퀴 해석의 링 축인장응력. 바닥판 휨 인장에 더해 함께 검토 (처음엔 0)"""
    F_up = DOME_UPLIFT_PER_PA * float(np.max(env["dp"])) if USE_INTERNAL_PRESSURE else 0.0
    tf = RING_FLANGE_T
    t_start = max(COVER_MIN_THICKNESS_M, t_floor or 0.0)
    for _ in range(20):
        if D - tf < STEM_H_MIN - 1e-9:
            return {"ok": False, "reason": f"바닥판 {tf:.1f} m가 깊이 {D:.1f} m 안에 안 들어감"}
        t, fl = t_start, None
        while t <= COVER_T_MAX_M + 1e-9:
            fl = foundation_loads(D, B, round(t, 3), tf, dome_mass_ton)
            if fl is None:
                return {"ok": False, "reason": "링 단면 불성립"}
            if fl["uplift_resist"] >= FS_UPLIFT * F_up:
                break
            t += COVER_T_STEP_M
        else:
            return {"ok": False, "reason": f"덮개 {COVER_T_MAX_M} m까지 들림 불성립"}
        t = round(t, 3)
        ev = check_envelope(0.0, cover_t=t, verbose=False) if USE_MEMBRANE else {"dcr_lc2": 0.0}
        if ev["dcr_lc2"] >= 1.0:                               # 덮개를 더 올려도 나빠지기만 함 -> 이 칸 탈락
            return {"ok": False, "reason": f"외피 감압 장력 초과 (덮개 {t:.1f} m, DCR {ev['dcr_lc2']:.2f})", "t": t}
        if fl["q_lc2"] > fl["q_allow"]:
            return {"ok": False, "reason": f"지지력 부족 ({fl['q_lc2']/1e3:.0f} > {fl['q_allow']/1e3:.0f} kPa)", "t": t}
        fb = flange_bending(D, B, t, tf, fl["q_lc2"], sig_hoop)
        rebar = 0.0
        if FLANGE_MODE == "rebar":
            rebar = flange_rebar_kg(D, B, tf, fb["M_lc2"], fb["M_lc1"])
        else:
            if not math.isfinite(fb["h_plain"]):
                return {"ok": False, "reason": f"링 인장 {sig_hoop/1e6:.2f} MPa가 바닥판 허용 인장응력을 다 씀"}
            need = max(RING_FLANGE_T, math.ceil(fb["h_plain"] / FLANGE_T_STEP - 1e-9) * FLANGE_T_STEP)
            fit = math.floor((D - STEM_H_MIN) / FLANGE_T_STEP + 1e-9) * FLANGE_T_STEP
            if need > tf + 1e-9:
                if need <= fit + 1e-9:
                    tf = round(need, 3)
                    continue                                # 두꺼워진 바닥판으로 ①부터 다시
                if fit > tf + 1e-9:
                    tf = round(fit, 3)                      # 들어가는 최대 두께로 한 번 더 확인 (두꺼우면 자중이 늘어 필요 두께가 줄 수 있음)
                    continue
                if FLANGE_MODE == "plain":
                    return {"ok": False, "reason": f"무근 바닥판 {need:.1f} m가 깊이 {D:.1f} m 안에 안 들어감"}
                rebar = flange_rebar_kg(D, B, tf, fb["M_lc2"], fb["M_lc1"])
        return {"ok": True, "design": Design(D, t, B, tf, float(rebar)), "F_up": F_up,
                "uplift_ratio": fl["uplift_resist"] / max(F_up, 1e-9), "q_lc2": fl["q_lc2"], "q_allow": fl["q_allow"],
                "flange": fb, "W": fl["W"], "envelope_lc2_dcr": ev["dcr_lc2"],
                "t_governs": "방사선" if t <= COVER_MIN_THICKNESS_M + 1e-9 else "구조(들림)"}
    return {"ok": False, "reason": "바닥판 두께 반복이 수렴하지 않음"}


def bare_cover(env, dome_mass_ton=0.0, t_floor=None):
    """기초 없음: (덮개 + 외피 + 돔) 무게 >= FS_UPLIFT x 들림 인 최소 덮개 t (방사선 하한부터).
    design_cell ①에서 링·바닥판 위 흙·마찰 몫을 뺀 것과 같다."""
    g = planet_g()
    F_up = DOME_UPLIFT_PER_PA * float(np.max(env["dp"])) if USE_INTERNAL_PRESSURE else 0.0
    W_fix = ((MEMBRANE_AREAL_DENSITY * DOME_SHELL_AREA_M2 if USE_MEMBRANE else 0.0) + dome_mass_ton * 1000.0) * g
    t = max(COVER_MIN_THICKNESS_M, t_floor or 0.0)
    while t <= COVER_T_MAX_M + 1e-9:
        t = round(t, 3)
        if COVER_DENSITY * cover_volume_m3(t) * g + W_fix >= FS_UPLIFT * F_up:
            break
        t += COVER_T_STEP_M
    else:
        raise RuntimeError(f"기초 없음: 덮개 {COVER_T_MAX_M} m까지 들림을 못 막음 -> COVER_T_MAX_M을 올리세요.")
    ev = check_envelope(0.0, cover_t=t, verbose=False) if USE_MEMBRANE else {"dcr_lc2": 0.0}
    if ev["dcr_lc2"] >= 1.0:
        raise RuntimeError(f"기초 없음: 외피 감압 장력 초과 (덮개 {t:.1f} m, DCR {ev['dcr_lc2']:.2f})")
    print(f"[none] 돔 {dome_mass_ton:.0f} t -> 덮개 {t:.1f} m (누르는 무게 "
          f"{(COVER_DENSITY * cover_volume_m3(t) * g + W_fix)/1e6:,.1f} MN >= {FS_UPLIFT} x 들림 {F_up/1e6:,.1f} MN)")
    return t


# ----------------------------------------------------------------------
# 6. 해석 단계
# ----------------------------------------------------------------------
def apply_gravity(node_mass, extra_weight=None, ts_tag=1, pattern_tag=1, verbose=True, weight_exempt=None):
    """node_mass: 질량(관성 + 중력하중), extra_weight: 관성 없이 무게만 거는 질량 (덮개 COVER_INERTIA=False)
    weight_exempt: 질량(관성)은 두되 무게는 절점에 걸지 않을 몫 (덮개 -> 부재 선하중으로 따로 건다)"""
    g = planet_g()
    ops.timeSeries("Linear", ts_tag)
    ops.pattern("Plain", pattern_tag, ts_tag)
    weight = defaultdict(float)
    for node, mass in node_mass.items():
        ops.mass(node, mass, mass, mass, 0.0, 0.0, 0.0)
        weight[node] += mass
    for node, mass in (extra_weight or {}).items():
        weight[node] += mass
    for node, mass in (weight_exempt or {}).items():
        weight[node] -= mass
    for node, mass in weight.items():
        if abs(mass) > 1e-9:
            ops.load(node, 0.0, 0.0, -mass * g, 0.0, 0.0, 0.0)
    if verbose:
        W = (sum(node_mass.values()) + sum((extra_weight or {}).values())) * g
        W_line = sum((weight_exempt or {}).values()) * g
        print(f"[gravity] 중력 하중 {W/1e6:,.2f} MN (g={g}) -> {len(node_mass)}개 절점"
              + (f", 그중 덮개 {W_line/1e6:,.2f} MN은 부재 선하중" if W_line else ""))


def apply_dead_loads(model, design, verbose=True):
    """LC2 하중: 강재·링·외피 무게(절점) + 덮개 무게(USE_LINE_LOADS면 부재 선하중, 아니면 절점). 반환: assemble_masses와 같음"""
    masses, extra_weight, mass_parts, cover = assemble_masses(model, design, verbose=verbose)
    line_cover = USE_LINE_LOADS and design.cover_t > 0.0
    apply_gravity(masses, extra_weight, verbose=verbose, weight_exempt=cover if line_cover else None)
    if line_cover:
        apply_shell_loads(model, cover_panel_forces(design.cover_t), {}, pattern_tag=3, ts_tag=1)
    return masses, extra_weight, mass_parts


def run_static_robust(n_steps=10, max_halving=8, tol=1.0e-6, verbose=True, label="정적"):
    """하중계수 0 -> 1. 수렴하지 않으면 증분을 반으로 줄여 다시 시도한다(실패한 스텝은 자동으로 되돌려짐).
    끝나면 loadConst로 현재 하중을 고정하고 시간을 0으로 되돌린다."""
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", tol, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / n_steps)
    ops.analysis("Static")

    lam, dl, halvings = 0.0, 1.0 / n_steps, 0
    while lam < 1.0 - 1e-10:
        step = min(dl, 1.0 - lam)
        ops.integrator("LoadControl", step)
        ok = ops.analyze(1)
        if ok != 0:
            ops.algorithm("ModifiedNewton", "-initial")
            ok = ops.analyze(1)
            ops.algorithm("Newton")
        if ok == 0:
            lam = ops.getTime()
            dl = min(2.0 * dl, 1.0 / n_steps)     # 수렴하면 증분을 다시 키운다
            continue
        halvings += 1
        if halvings > max_halving:
            if verbose:
                print(f"[static] {label} 해석 불수렴 (하중계수 {lam:.3f}에서 중단)")
            return -1
        dl = step / 2.0
    ops.loadConst("-time", 0.0)
    if verbose:
        print(f"[static] {label} 해석 성공" + (f" (증분 {halvings}회 축소)" if halvings else ""))
    return 0


def _run_linear_static(setup=True):
    """선형 정적해석 1회. setup=False면 이미 만든 해석(강성 1회 분해, -factorOnce)을 재사용해 하중만 바꿔 푼다."""
    if setup:
        ops.wipeAnalysis()
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1.0e-12, 6)
        ops.algorithm("Linear", "-factorOnce")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("단위하중 선형 정적해석 실패")


SEISMIC_UNIT_CASES = ("Q1u", "Q2u", "Q3u")   # 방향(1=X, 2=Y, 3=Z)별 질량비례 단위가속도 (잔여질량 계산에 사용)


def apply_tag_forces(forces_by_tag, pattern_tag, ts_tag=1):
    """OpenSees 절점 태그 기준 절점력 {tag: (Fx, Fy, Fz)}"""
    ops.pattern("Plain", pattern_tag, ts_tag)
    for tag, f in forces_by_tag.items():
        ops.load(int(tag), float(f[0]), float(f[1]), float(f[2]), 0.0, 0.0, 0.0)


def seismic_unit_node_forces(node_mass, direction):
    """지반가속도 1 m/s^2의 관성력 F_i = m_i (질량비례) {OpenSees 절점 태그: (Fx, Fy, Fz)}.
    질량에는 강재·막·덮개(관성 포함 시)·링이 모두 들어간다."""
    out = {}
    for t, m in node_mass.items():
        f = np.zeros(3)
        f[direction - 1] = m
        out[t] = f
    return out


_MODE_COUNT_HINT = [None]   # 직전 해석에서 조건을 만족한 모드 수 (다음 해석은 여기서 바로 시작)


def modal_analysis(node_mass, f_rigid_hz=None, steps=SEISMIC_MODE_STEPS):
    """현재 (선형) 모델에 질량을 넣고 고유치 해석. 반환: 주기, 참여계수 Γ(모드, 방향), 유효질량비, 모드형상.
    Γ_id = Σ m φ_d / Σ m |φ|²,  유효질량 M_id = (Σ m φ_d)² / Σ m |φ|²  (병진 질량만, 회전 관성 0)
    모드 수는 (마지막 모드 진동수 >= f_rigid_hz) 그리고 (방향별 유효질량 >= SEISMIC_MASS_TARGET)이 될 때까지 늘린다.
    그러면 빠진 모드는 지진 기록에 없는 높은 진동수라 지반과 함께 움직이고(강체), 잔여질량 하중이 그 몫을 정확히 덮는다."""
    tags = list(node_mass.keys())
    m = np.array([node_mass[t] for t in tags], dtype=float)
    for t, mi in zip(tags, m):
        ops.mass(int(t), mi, mi, mi, 0.0, 0.0, 0.0)
    lam, n_modes = None, 0
    if _MODE_COUNT_HINT[0]:
        steps = (_MODE_COUNT_HINT[0],) + tuple(n for n in steps if n > _MODE_COUNT_HINT[0])
    for n_try in steps:
        n_try = int(min(n_try, 3 * len(tags) - 2))
        if n_try <= n_modes:
            break
        ops.wipeAnalysis()
        try:
            res = ops.eigen(n_try)
        except Exception:
            res = None
        if not res or len(res) < n_try:
            if n_modes:                              # 고유치 해석 한계 -> 직전 모드 수로 다시 풀어 고유벡터 복원
                ops.wipeAnalysis()
                lam = np.array(ops.eigen(n_modes), dtype=float)
            break
        lam, n_modes = np.array(res, dtype=float), n_try
        mass_ok = _modal_mass_ratio(tags, m, n_modes).min() >= SEISMIC_MASS_TARGET
        if (f_rigid_hz is None or math.sqrt(max(lam[-1], 0.0)) / (2 * math.pi) >= f_rigid_hz) and mass_ok:
            break
    if lam is None:
        raise RuntimeError("지진용 고유치 해석에 실패했습니다.")
    _MODE_COUNT_HINT[0] = n_modes
    phi = np.array([[ops.nodeEigenvector(int(t), i + 1)[:3] for t in tags] for i in range(n_modes)])  # (n, 절점, 3)
    mphi = m[None, :, None] * phi
    gen_mass = np.einsum("nkd,nkd->n", mphi, phi)
    L = mphi.sum(axis=1)                               # (n, 3) = Σ m φ_d
    gamma = L / gen_mass[:, None]
    meff = L ** 2 / gen_mass[:, None]
    omega = np.sqrt(np.maximum(lam, 0.0))
    return {"tags": tags, "m": m, "omega": omega, "T": 2.0 * math.pi / np.maximum(omega, 1e-12),
            "phi": phi, "gamma": gamma, "meff_ratio": meff / m.sum(), "M_total": float(m.sum())}


def _modal_mass_ratio(tags, m, n_modes):
    phi = np.array([[ops.nodeEigenvector(int(t), i + 1)[:3] for t in tags] for i in range(n_modes)])
    mphi = m[None, :, None] * phi
    L = mphi.sum(axis=1)
    return (L ** 2 / np.einsum("nkd,nkd->n", mphi, phi)[:, None]).sum(axis=0) / m.sum()


def compute_unit_load_responses(section, design, contact, env_e, node_mass=None, f_rigid_hz=None):
    """LC1 기준상태 접선강성의 선형 모델에서 단위하중별 응답. 실측 하중 = 단위하중 x 실측값.
      P: 순압 +1 Pa,  T: 부재 온도 +1 K (열 보정 포함, 링은 ΔT=0),  W[s]: 측정높이 풍압 1 Pa, s번째 방위
      (node_mass가 있으면 지진용으로 추가)
      Q1u/Q2u/Q3u: 지반가속도 1 m/s^2 질량비례 관성력,  M{i}: i번째 모드 관성력 형상 m·φ_i
    부재력(F), 지반 반력(R), 링 축력(N), 돔 절점변위를 모두 기록한다."""
    model = build_structure(section, design, transf_type="Linear",
                            contact=contact, verbose=False)
    EA_alpha = section["E"] * section["A"] * MEMBER_ALPHA
    step = 360.0 / WIND_SECTORS

    modal = None
    cases = []   # (이름, 하중, 종류: "shell" = 외피 하중(선하중), "node" = 원래 절점 번호, "tag" = OpenSees 태그)
    if USE_INTERNAL_PRESSURE:
        cases.append(("P", (pressure_panel_forces(1.0), pressure_node_forces(dome_NODES, DOME_EDGES, 1.0)), "shell"))
    if USE_THERMAL:
        cases.append(("T", thermal_node_forces(dome_NODES, DOME_EDGES, EA_alpha * 1.0), "node"))
    if env_e.get("wind_active", USE_WIND):
        cases += [(f"W{s}", (wind_panel_forces(1.0, s * step), wind_node_forces(dome_NODES, DOME_EDGES, 1.0, s * step)),
                   "shell") for s in range(WIND_SECTORS)]
    if node_mass:
        modal = modal_analysis(node_mass, f_rigid_hz)  # 질량은 정적해석에 영향 없음
        for name in SEISMIC_UNIT_CASES:
            cases.append((name, seismic_unit_node_forces(node_mass, int(name[1])), "tag"))
        for i in range(len(modal["T"])):
            mphi = modal["m"][:, None] * modal["phi"][i]
            cases.append((f"M{i}", dict(zip(modal["tags"], mphi)), "tag"))

    ops.timeSeries("Linear", 1)
    F, R, N, disp = {}, {}, {}, {}
    for k, (name, forces, kind) in enumerate(cases):
        tag = 100 + k
        if kind == "tag":
            apply_tag_forces(forces, pattern_tag=tag, ts_tag=1)
        elif kind == "shell":
            apply_shell_loads(model, forces[0], forces[1], pattern_tag=tag, ts_tag=1)
        else:
            apply_nodal_forces(model["node_tags"], forces, pattern_tag=tag, ts_tag=1)
        _run_linear_static(setup=(k == 0))      # 강성은 첫 케이스에서 한 번만 분해
        F[name] = element_local_forces()
        R[name] = ground_reactions(model)
        N[name] = ring_axial(model)
        disp[name] = capture_node_disp(model["node_tags"])
        ops.remove("loadPattern", tag)
        _LINE_LOAD_STATE.clear()          # 이 케이스의 선하중도 같이 지운다
        ops.reset()                       # 변위·시간을 0으로 되돌려 다음 단위하중과 섞이지 않게 한다
    ops.wipe()

    def pack(d, correct=None):
        out = {"P": d.get("P"), "T": d.get("T"), "W": None,
               "Q": {n: d[n] for n in SEISMIC_UNIT_CASES if n in d} or None, "modes": None}
        if correct is not None and out["T"] is not None:
            out["T"] = correct(out["T"])
        if "W0" in d:
            out["W"] = np.stack([d[f"W{s}"] for s in range(WIND_SECTORS)])
        if "M0" in d:
            out["modes"] = np.stack([d[f"M{i}"] for i in range(len(modal["T"]))])
        return out

    return {
        "F": pack(F, correct=lambda f: correct_thermal_forces(f, EA_alpha * 1.0)),
        "R": pack(R),
        "N": pack(N),
        "D": pack(disp),                           # 돔 절점변위 (절점수, 3)
        "modal": modal,
        "EA_alpha": EA_alpha,
        "node_order": list(model["node_tags"].keys()),
    }


# ----------------------------------------------------------------------
# 6-2. 지진 등가정적 상한 — 모드별 최대값의 절대합 + 잔여질량, 방향 100/30/30
#
#   선형 구조의 지진 응답 r(t) = Σ_i r_i · Γ_i · A_i(t)   (r_i: 모드 i 관성력 형상 m·φ_i에 대한 정적 응답,
#                                                        A_i: 모드 i 의사가속도, |A_i| <= PSA(T_i) <= Sa_max)
#   -> |r(t)| <= Σ_i |r_i · Γ_i| · Sa_max  (모든 모드가 같은 순간 최대·같은 부호라고 가정한 절대합)
#   각 항을 '생성자' g 로 두면 지진 응답은 {Σ s_g g : s_g ∈ [-1, 1]} (zonotope) 안에 있다.
#   잔여질량(모드가 못 담은 질량)은 지반과 함께 움직인다고 보고 질량비례 하중 x Sa_max로 더한다.
# ----------------------------------------------------------------------
def seismic_spectral_matrix(seis, T, sa):
    """모드 주기마다 방향별 스펙트럼 값 S(i, d). 기본은 주기 무관 최대값(모든 모드에 Sa_max)."""
    n = len(T)
    if SEISMIC_USE_PERIOD_ENVELOPE:
        return np.tile([sa[1], sa[2], sa[3]], (n, 1)).astype(float)
    S = np.array([[spectral_value(seis["comps"][d]["spec"], t) for d in (1, 2, 3)] for t in T])
    if SEISMIC_HORIZONTAL_ENVELOPE:
        S[:, 0] = S[:, 1] = S[:, :2].max(axis=1)
    return S


def build_quake_generators(unit, seis, sa, ortho=SEISMIC_ORTHO_FACTOR, directions=(1, 2, 3)):
    """지진 상한의 생성자 (주방향 p, 생성자 g, ...) — 주방향 X·Y·Z 각각 100%, 나머지 ortho배.
      모드 생성자  : (Σ_d f_pd |Γ_id| S_id) x r_i
      잔여질량 생성자: f_pd Sa_d x (r_uniform,d - Σ_i Γ_id r_i)
    directions로 성분을 제한하면 그 성분만의 상한(시간이력 비교용)이 된다."""
    md = unit["modal"]
    G = md["gamma"]
    S = seismic_spectral_matrix(seis, md["T"], sa)
    use = np.array([1.0 if d in directions else 0.0 for d in (1, 2, 3)])
    prim = [d for d in (1, 2, 3) if d in directions]
    fac = np.array([[(1.0 if d == p else ortho) * use[d - 1] for d in (1, 2, 3)] for p in prim])   # (P, 3)
    a_mode = (np.abs(G)[None, :, :] * S[None, :, :] * fac[:, None, :]).sum(axis=2)                # (P, n)
    sa_vec = np.array([sa[1], sa[2], sa[3]], dtype=float)
    res_coef = fac * sa_vec[None, :]                                                             # (P, 3)
    out = {"coef": np.concatenate([a_mode, res_coef], axis=1), "primary": prim,
           "labels": [f"mode{i + 1}" for i in range(len(md["T"]))] + ["residual_X", "residual_Y", "residual_Z"]}
    for key in ("F", "R", "N", "D"):
        modes = unit[key]["modes"]                                                               # (n, ...)
        uni = np.stack([unit[key]["Q"][f"Q{d}u"] for d in (1, 2, 3)])                            # (3, ...)
        res = uni - np.tensordot(G.T, modes, axes=(1, 0))                                        # (3, ...)
        tail = (1,) * (modes.ndim - 1)
        out[key] = np.stack([np.concatenate([a_mode[k].reshape((-1,) + tail) * modes,
                                             res_coef[k].reshape((-1,) + tail) * res])
                             for k in range(len(prim))])                                         # (P, n+3, ...)
    return out


def quake_bound_abs(gen_arr):
    """생성자 절대합 = 각 성분의 상한 (주방향 중 최대). gen_arr: (P, G, ...) -> (...)"""
    return np.abs(gen_arr).sum(axis=1).max(axis=0)


def volume_change_m3(disp, node_order):
    """돔 절점변위 (절점수, 3)에 의한 내부 체적 변화 ΔV = Σ 패널 외향면적 · 평균변위 (바닥 = 지반, 고정)."""
    idx = {tag: i for i, tag in enumerate(node_order)}
    dv = 0.0
    for ids, w, area_vec, _ in _oriented_shell_panels(dome_NODES, DOME_EDGES):
        u = sum(wk * disp[idx[k]] for k, wk in zip(ids, w))
        dv += float(np.dot(area_vec, u))
    return dv


def combine_env(base, parts, env):
    """환경 타임라인 응답 (시각수, *base.shape)
    X(t) = X_기준 + (Δp(t)-Δp평균)·P + ΔT(t)·T + q(t)·W[풍향(t)]"""
    base = np.asarray(base, dtype=float)
    nt = len(env["t_s"])
    out = np.repeat(base[None], nt, axis=0)
    if base.size == 0:
        return out
    shp = (nt,) + (1,) * base.ndim
    if USE_INTERNAL_PRESSURE and parts.get("P") is not None:
        out = out + (env["dp"] - env["dp_ref"]).reshape(shp) * parts["P"][None]
    if USE_THERMAL and parts.get("T") is not None:
        out = out + np.asarray(env["dT"]).reshape(shp) * parts["T"][None]
    if env.get("wind_active", USE_WIND) and parts.get("W") is not None:
        out = out + np.asarray(env["q_Pa"]).reshape(shp) * parts["W"][env["sector"]]
    return out


# [시간이력 검증 전용] Rayleigh 감쇠 (1, 3차 모드 기준 5%)
def setup_rayleigh_damping(xi=0.05, mode_i=1, mode_j=3):
    lam = ops.eigen(mode_j)
    wi = lam[mode_i - 1] ** 0.5
    wj = lam[mode_j - 1] ** 0.5
    a0 = 2 * xi * wi * wj / (wi + wj)
    a1 = 2 * xi / (wi + wj)
    ops.rayleigh(a0, 0.0, 0.0, a1)
    print(f"[seismic] Rayleigh 감쇠 적용: xi={xi} (a0={a0:.6f}, a1={a1:.6f})")
    return a0, a1


# [시간이력 검증 전용] 지진파 입력
def apply_earthquake(acc, dt, direction=1, ts_tag=2, pattern_tag=2):
    ops.timeSeries("Path", ts_tag, "-dt", dt, "-values", *acc, "-factor", 1.0)
    ops.pattern("UniformExcitation", pattern_tag, direction, "-accel", ts_tag)
    print(f"[seismic] 지진 하중 등록 ({len(acc)}개 값, dt={dt})")


# [시간이력 검증 전용] 동적 해석
def run_dynamic_robust(num_steps, dt, primary="KrylovNewton", fallback="ModifiedNewton", progress_every=5000):
    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-8, 50)
    ops.algorithm(primary)
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    fail_count = 0
    for step in range(num_steps):
        ok = ops.analyze(1, dt)
        if ok != 0:
            ops.algorithm(fallback, "-initial")
            ok = ops.analyze(1, dt)
            ops.algorithm(primary)
            if ok != 0:
                fail_count += 1
                print(f"[seismic] 경고: 수렴 실패 step {step} (t={ops.getTime():.2f}s)")
        if progress_every and step % progress_every == 0:
            print(f"[seismic] 진행률 {step/num_steps*100:.0f}% (t={ops.getTime():.1f}s)")

    print(f"[seismic] 동적 해석 완료: {'성공' if fail_count == 0 else f'{fail_count}스텝 수렴실패'}")
    return fail_count == 0


# ----------------------------------------------------------------------
# 7. recorder / 결과 읽기
#    지진 증분(수십 N)이 기준 축력(수 MN)에 묻히지 않도록 유효숫자 12자리로 기록
# ----------------------------------------------------------------------
RECORDER_PRECISION = 12


def setup_member_force_recorders(elements, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for ele_tag in elements.values():
        path = os.path.join(out_dir, f"force_ele{ele_tag}.txt")
        ops.recorder("Element", "-file", path, "-precision", RECORDER_PRECISION,
                     "-time", "-ele", ele_tag, "localForce")
        paths[ele_tag] = path
    print(f"[member] 부재 {len(elements)}개의 부재력(localForce) recorder 등록 완료")
    return paths


def setup_foundation_envelope_recorders(model, out_dir):
    """지진 중 지반 반력·링 부재력의 최소/최대만 기록 (Envelope: 행 = min, max, absmax)."""
    os.makedirs(out_dir, exist_ok=True)
    g = model["ground"]
    paths = {"R": os.path.join(out_dir, "env_ground_reaction.txt")}
    ops.recorder("EnvelopeNode", "-file", paths["R"], "-precision", RECORDER_PRECISION,
                 "-node", *g["fixed"], "-dof", 1, 2, 3, "reaction")
    if g["ring_eles"]:
        paths["N"] = os.path.join(out_dir, "env_ring_force.txt")
        ops.recorder("EnvelopeElement", "-file", paths["N"], "-precision", RECORDER_PRECISION,
                     "-ele", *g["ring_eles"], "localForce")
    return paths


def _read_table(path):
    """공백 구분 recorder 파일을 빠르게 읽는다 (np.loadtxt보다 수 배 빠름)."""
    return pd.read_csv(path, sep=r"\s+", header=None, dtype=float).to_numpy()


def _read_envelope(path):
    return np.atleast_2d(np.loadtxt(path))


def get_member_lengths(nodes, edges):
    """요소별 좌굴 길이. 곡선 리브 돔은 요소 하나가 아니라 접합부 사이 리브 전체(호 길이)가 좌굴 길이."""
    if DOME_BUCKLING_LENGTH is not None and edges is DOME_EDGES:
        return {ele: DOME_BUCKLING_LENGTH[ele] for ele, _, _ in edges}
    return {ele: float(np.linalg.norm(np.asarray(nodes[nj], float) - np.asarray(nodes[ni], float)))
            for ele, ni, nj in edges}


# ----------------------------------------------------------------------
# 8. 판정
#
#    부재: DCR = |N|/Pn + |My|/Mny + |Mz|/Mnz  (인장 Pn=Fy·A, 압축 Pn=min(Fy·A, Euler Pcr))
#          DCR은 (N, My, Mz)에 대한 볼록함수 -> '환경(t) + 지진증분(τ)' 조합의 최대는
#          지진증분 점들의 볼록껍질 꼭짓점에서만 나온다.
#          LC1(가압 타임라인 [+지진])과 LC2(감압) 중 큰 값이 부재의 DCR.
#    기초: 들림(플랜지 위 흙 / FS), 지지력(q_ult / FS), 활동(수평 한계 / FS), 링 압축·인장(철근)
# ----------------------------------------------------------------------
BUCKLING_K = 1.0   # 항복강도는 [재료 파트]의 MEMBER_FY


def _end_forces(F, end):
    """부재력(…, 12 또는 18) -> 단면의 (N, My, Mz). N>0 인장. end: "I", "J", "M"(중앙, 18열일 때)"""
    if end == "I":
        return np.stack([-F[..., 0], F[..., 4], F[..., 5]], axis=-1)
    if end == "J":
        return np.stack([F[..., 6], F[..., 10], F[..., 11]], axis=-1)
    return np.stack([F[..., 12], F[..., 16], F[..., 17]], axis=-1)


def _combo_extremes(E, Q, cap, max_pairs=2_000_000):
    """모든 (E[i] + Q[j]) 조합에서 최대 DCR과 최대 |응력|(부호 포함, 3D 그림 색)을 구한다."""
    Pt, Pc, Mny, Mnz, A, Iy, Iz, c = cap
    best = {"dcr": -np.inf, "i": 0, "j": 0, "NMM": (0.0, 0.0, 0.0), "sigma": 0.0}
    rows = max(1, max_pairs // max(len(Q), 1))
    for s in range(0, len(E), rows):
        S = E[s:s + rows, None, :] + Q[None, :, :]
        N, My, Mz = S[..., 0], S[..., 1], S[..., 2]

        dcr = np.where(N >= 0.0, N / Pt, -N / Pc) + np.abs(My) / Mny + np.abs(Mz) / Mnz
        k = int(np.argmax(dcr))
        if dcr.flat[k] > best["dcr"]:
            i, j = np.unravel_index(k, dcr.shape)
            best.update(dcr=float(dcr.flat[k]), i=s + int(i), j=int(j),
                        NMM=(float(N[i, j]), float(My[i, j]), float(Mz[i, j])))

        bend = np.abs(My) * c / Iy + np.abs(Mz) * c / Iz
        sp, sm = N / A + bend, N / A - bend
        sig = np.where(np.abs(sp) >= np.abs(sm), sp, sm)
        k2 = int(np.argmax(np.abs(sig)))
        if abs(sig.flat[k2]) > abs(best["sigma"]):
            best["sigma"] = float(sig.flat[k2])
    return best


def _dcr_planes(Pt, Pc, Mny, Mnz):
    """DCR = max(N/Pt, -N/Pc) + |My|/Mny + |Mz|/Mnz 를 8개 평면의 최대값으로 쓴 계수 (8, 3).
    DCR이 이 평면들의 최대값이라 '환경 + λ·지진' 조합의 한계 λ를 닫힌 식으로 구할 수 있다."""
    return np.array([[aN, sy / Mny, sz / Mnz]
                     for aN in (1.0 / Pt, -1.0 / Pc) for sy in (1.0, -1.0) for sz in (1.0, -1.0)])


def _quake_margin(E_pts, Q_pts, planes):
    """λ_crit: 하루 중 어느 시각이든 '환경(t) + λ·지진(조합 j)'의 DCR이 1에 닿는 최소 배율.
    DCR_max(λ) = max_s [max_t c_s·E_t + λ·max_j c_s·Q_j] 이므로 λ = min_s (1 - max_t c_s·E_t) / max_j c_s·Q_j."""
    cE = (E_pts @ planes.T).max(axis=0)
    cQ = (Q_pts @ planes.T).max(axis=0)
    if np.any(cE >= 1.0):
        return 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(cQ > 0.0, (1.0 - cE) / cQ, np.inf)
    return float(lam.min())


def _quake_dcr_increment(Q_pts, planes):
    """지진만으로 생기는 DCR 증가분의 상한 max_j max_s c_s·Q_j (등가정적 vs 시간이력 비교용)"""
    return float(max((Q_pts @ planes.T).max(), 0.0))


def _zonotope_vertices(gen, end, planes):
    """지진 상한 zonotope {Σ s_g g}에서 DCR 평면 8개 각각을 최대로 만드는 꼭짓점 (주방향 P개 x 8, 3).
    DCR이 평면들의 최대값이므로 이 꼭짓점만 보면 zonotope 전체의 최대 DCR이 정확히 나온다.
    gen: (P, G, 12) 한 부재의 생성자 localForce. 반환: (꼭짓점, 축력 증분 상한)"""
    V, dN = [], 0.0
    for gp in gen:
        Gp = _end_forces(gp, end)                       # (G, 3) = (N, My, Mz)
        S = Gp @ planes.T                               # (G, 8)
        sgn = np.where(S >= 0.0, 1.0, -1.0)
        V.append(sgn.T @ Gp)                            # (8, 3)
        dN = max(dN, float(np.abs(Gp[:, 0]).sum()))
    return np.vstack(V), dN


def assess_members(F_env, F_lc2, Q_gen, section, c, lengths, env,
                   actual_idx=None, fy=MEMBER_FY, K=BUCKLING_K):
    """부재별 파손 판정.
      - LC1 환경 단독 : 타임라인 모든 시각
      - LC1 + 지진    : 지진 상한 zonotope(모드 절대합 + 잔여질량, 100/30/30)을 타임라인 모든 시각에 더함
                        = 지진이 하루 중 어느 시각에 와도 견디는지 (Q_gen=None이면 생략)
      - LC2 감압      : 자중 + 덮개, 내압 0
    Q_gen: (주방향, 생성자, 부재수, 12) 지진 상한 생성자 부재력 (LC1 기준상태 대비 증분)
    반환: (부재별 결과 dict, 시각별 LC1 환경 단독 최대 DCR 배열)"""
    A, E, G, J, Iy, Iz = section
    I_min = min(Iy, Iz)
    r = (I_min / A) ** 0.5
    Pt, Mny, Mnz = fy * A, fy * Iy / c, fy * Iz / c
    zero = np.zeros((1, 3))

    env_dcr_t = np.zeros(F_env.shape[0])
    result = {}
    for e_idx, ele in enumerate(ELE_TAGS):
        L = lengths[ele]
        p_cr = math.pi ** 2 * E * I_min / (K * L) ** 2
        Pc = min(Pt, p_cr)
        cap = (Pt, Pc, Mny, Mnz, A, Iy, Iz, c)
        planes = _dcr_planes(Pt, Pc, Mny, Mnz)

        gen = None if Q_gen is None else Q_gen[:, :, e_idx, :]

        best1, best2 = None, None
        best_actual, env_max, dN_quake, sigma_any = -np.inf, 0.0, 0.0, 0.0
        lam, dcr_q = np.inf, 0.0
        for end in MEMBER_SECTIONS:
            E_pts = _end_forces(F_env[:, e_idx, :], end)
            if gen is not None:
                Q_pts, dN_b = _zonotope_vertices(gen, end, planes)
                dN_quake = max(dN_quake, dN_b)
                Q = np.vstack([zero, Q_pts])                     # 0 = 지진 없는 순간
                lam = min(lam, _quake_margin(E_pts, Q_pts, planes))
                dcr_q = max(dcr_q, _quake_dcr_increment(Q_pts, planes))
            else:
                Q = zero
            b1 = _combo_extremes(E_pts, Q, cap)
            env_only = b1 if gen is None else _combo_extremes(E_pts, zero, cap)
            env_max = max(env_max, env_only["dcr"])
            b2 = _combo_extremes(_end_forces(F_lc2[e_idx][None, :], end), zero, cap)

            N, My, Mz = E_pts[:, 0], E_pts[:, 1], E_pts[:, 2]
            env_dcr_t = np.maximum(env_dcr_t, np.where(N >= 0, N / Pt, -N / Pc)
                                   + np.abs(My) / Mny + np.abs(Mz) / Mnz)
            if gen is not None and actual_idx is not None and len(actual_idx):
                best_actual = max(best_actual, _combo_extremes(E_pts[actual_idx], Q[1:], cap)["dcr"])

            for b in (b1, b2):
                if abs(b["sigma"]) > abs(sigma_any):
                    sigma_any = b["sigma"]
            if best1 is None or b1["dcr"] > best1["dcr"]:
                best1 = {**b1, "end": end}
            if best2 is None or b2["dcr"] > best2["dcr"]:
                best2 = {**b2, "end": end}

        if best2["dcr"] > best1["dcr"]:
            gov, gov_case = best2, "LC2 감압"
        else:
            gov, gov_case = best1, ("LC1 가압+지진" if best1["j"] > 0 else "LC1 가압")
        max_dcr = gov["dcr"]
        result[ele] = {
            "failed": max_dcr >= 1.0,
            "max_dcr": max_dcr,
            "max_dcr_lc2": best2["dcr"],
            "max_dcr_env_only": env_max,
            "max_dcr_actual_time": best_actual if np.isfinite(best_actual) else None,
            "quake_dcr_increment": dcr_q,
            "quake_margin": lam,
            "gov_case": gov_case,
            "gov_env_index": best1["i"],
            "gov_lmst_h": float(env["lmst_h"][best1["i"]]),
            "gov_end": gov["end"],
            "gov_end_lc1": best1["end"],
            "critical_N": gov["NMM"][0],
            "critical_My": gov["NMM"][1],
            "critical_Mz": gov["NMM"][2],
            "signed_stress": sigma_any,
            "slenderness": K * L / r,
            "compression_capacity": Pc,
            "quake_max_dN": dN_quake,
        }
    return result, env_dcr_t


def _bisect_margin(f, lam_hi=1e12, iters=200):
    """f(λ) <= 1을 만족하는 최대 λ (f는 λ에 대해 증가). f(0) >= 1이면 0."""
    if f(0.0) >= 1.0:
        return 0.0
    hi = 1.0
    while f(hi) < 1.0:
        hi *= 10.0
        if hi > lam_hi:
            return np.inf
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if f(mid) < 1.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-6 * hi:
            break
    return lo


_H_DIRS = np.stack([np.cos(np.linspace(0.0, 2.0 * np.pi, 361)[:-1]),
                    np.sin(np.linspace(0.0, 2.0 * np.pi, 361)[:-1])], axis=1)   # 수평 합력 방향 1° 간격


def _assess_bare_ground(out, R_env, R_lc2, R_gen=None):
    """기초 없음: ① 어느 발(밑단 절점)도 들리면 안 됨 (인장 0)
      ② 활동: 링 모드와 같은 방식. 돔이 바깥으로 미는 힘은 지표에 놓인 밑단 둘레 부재가 후프처럼 스스로 평형을 이루므로
         (발 하나가 국부 마찰 한계에 닿는 것은 허용), 돔 전체가 미끄러지는지를 '수평 반력 합력' vs 'μ x 발 수직반력 합'으로 본다.
    지진 상한을 타임라인 모든 시각에 더한다. 지지력·침하는 검토하지 않는다 (덮개 두께 비교용)."""
    n = R_lc2.shape[0]
    Bz = np.zeros(n) if R_gen is None else quake_bound_abs(R_gen[..., 2])
    if R_gen is None:
        supp = np.zeros(len(_H_DIRS))
    else:
        h = np.stack([R_gen[..., 0].sum(axis=2), R_gen[..., 1].sum(axis=2)], axis=-1)   # (P, G, 2) 수평 합력 생성자
        supp = np.abs(h @ _H_DIRS.T).sum(axis=1).max(axis=0)
    Rz_env_min = R_env[..., 2].min(axis=0)
    uH = np.stack([R_env[..., 0].sum(axis=1), R_env[..., 1].sum(axis=1)], axis=1) @ _H_DIRS.T
    net_lc2 = float(np.hypot(R_lc2[:, 0].sum(), R_lc2[:, 1].sum()))
    mu = math.tan(SITE_BASE_FRICTION_RATIO * math.radians(SITE_PHI_DEG))

    def dcr(L):
        Rz = np.minimum(Rz_env_min - L * Bz, R_lc2[:, 2])          # 발마다 가장 덜 눌린 순간
        if np.any(Rz <= 0.0):
            return np.inf, np.inf
        H = max(float((uH + L * supp[None, :]).max()), net_lc2)
        return 0.0, FS_SLIDING * H / (mu * float(Rz.sum()))

    up, slide = dcr(1.0)
    H_node = np.maximum(np.hypot(R_env[..., 0], R_env[..., 1]).max(axis=0), np.hypot(R_lc2[:, 0], R_lc2[:, 1]))
    reasons = [name for name, v in (("들림", up), ("활동", slide)) if v >= 1.0]
    out.update(mode="none", ok=not reasons, reason="/".join(reasons), uplift_dcr=up, bearing_dcr=0.0,
               sliding_dcr=slide, n_flip=0, n_flip_quake=0, q_max=0.0,
               n_h_limit=int(np.sum(H_node >= mu * np.minimum(Rz_env_min, R_lc2[:, 2]))),
               n_lifted=int(np.sum(Rz_env_min - Bz <= 0.0)), H_net=max(float(uH.max()), net_lc2),
               margin={"uplift": _bisect_margin(lambda L: dcr(L)[0]),
                       "sliding": _bisect_margin(lambda L: dcr(L)[1])})
    return out


def assess_ground(model, R_env, R_lc2, R_base, R_gen=None, U_lc2=None):
    """지반 판정. R_env (시각수, n, 3), R_lc2/R_base (n, 3), R_gen = 지진 상한 생성자 반력 (주방향, 생성자, n, 3).
    +z = 눌림(지반이 밀어 올림), -z = 들림(플랜지 위 흙이 잡아당김).
    지진 상한을 타임라인 모든 시각에 더한다(지진이 언제 와도). λ = 각 판정이 1에 닿는 지진 배율."""
    g = model["ground"]
    Rz_env = R_env[..., 2]
    out = {"sum_Rz_lc2": float(R_lc2[:, 2].sum()), "sum_Rz_lc1": float(R_base[:, 2].sum()),
           "min_Rz": float(Rz_env.min()), "max_Rz": float(max(Rz_env.max(), R_lc2[:, 2].max()))}
    if g["mode"] == "none":
        return _assess_bare_ground(out, R_env, R_lc2, R_gen)

    n = R_base.shape[0]
    Bz = np.zeros(n) if R_gen is None else quake_bound_abs(R_gen[..., 2])          # 지점별 |ΔRz| 상한
    if R_gen is None:
        supp = np.zeros(len(_H_DIRS))
    else:
        h = np.stack([R_gen[..., 0].sum(axis=2), R_gen[..., 1].sum(axis=2)], axis=-1)   # (P, G, 2) 수평 합력 생성자
        supp = np.abs(h @ _H_DIRS.T).sum(axis=1).max(axis=0)                         # 방향별 지지함수 (360,)
    Rz_env_min, Rz_env_max = Rz_env.min(axis=0), Rz_env.max(axis=0)
    Rz_min = Rz_env_min - Bz
    Rz_max = np.maximum(Rz_env_max + Bz, R_lc2[:, 2])

    uplift_dcr = np.maximum(-Rz_min, 0.0) * FS_UPLIFT / g["cap_up"]
    area = g["B"] * g["L_trib"]
    q = (np.maximum(Rz_max, 0.0) + g["soil_w"]) / area
    bearing_dcr = q / g["q_allow"]
    # 활동: 돔이 바깥으로 미는 힘은 링 후프 인장이 스스로 평형을 이룬다(지반은 일부만 분담, 국부 한계 도달 허용).
    #       링 전체가 미끄러지는지는 '지반 수평 반력의 합력' vs '링 전체 수평 저항'으로 본다.
    #       지진 상한 집합 위의 최대 |합력| = max_방향 (u·H_env(t) + 지지함수(u))
    Henv = np.stack([R_env[..., 0].sum(axis=1), R_env[..., 1].sum(axis=1)], axis=1)     # (nt, 2)
    uH = Henv @ _H_DIRS.T                                                               # (nt, 360)
    net_lc2 = float(np.hypot(R_lc2[:, 0].sum(), R_lc2[:, 1].sum()))
    H_of = lambda L: float((uH + L * supp[None, :]).max())
    H_net = max(H_of(1.0), net_lc2)
    cap_h_sum = float(g["cap_h"].sum())
    sliding_dcr = H_net * FS_SLIDING / cap_h_sum
    H_node = np.maximum(np.hypot(R_env[..., 0], R_env[..., 1]).max(axis=0), np.hypot(R_lc2[:, 0], R_lc2[:, 1]))
    n_h_limit = int(np.sum(H_node >= 0.98 * g["cap_h"]))

    # 접촉 상태: LC1 기준(눌림/들림)이 하루 중 + 지진 상한으로 바뀔 수 있는 지지점 수
    contact0 = R_base[:, 2] > 0.0
    flips_env = np.any((Rz_env > 0.0) != contact0[None, :], axis=0)
    flips_all = flips_env | (contact0 & (Rz_env_min - Bz <= 0.0)) | (~contact0 & (Rz_env_max + Bz > 0.0))

    # 지진 여유배율 λ (환경은 그대로, 지진만 λ배)
    lam_up = _bisect_margin(lambda L: float((np.maximum(-(Rz_env_min - L * Bz), 0.0)
                                             * FS_UPLIFT / g["cap_up"]).max()))
    lam_bear = _bisect_margin(lambda L: float(((np.maximum(Rz_env_max + L * Bz, 0.0)
                                                + g["soil_w"]) / area / g["q_allow"]).max()))
    lam_slide = _bisect_margin(lambda L: H_of(L) * FS_SLIDING / cap_h_sum)

    out.update(uplift_dcr=float(uplift_dcr.max()), bearing_dcr=float(bearing_dcr.max()),
               sliding_dcr=float(sliding_dcr), H_net=H_net, n_h_limit=n_h_limit, n_flip=int(flips_env.sum()),
               n_flip_quake=int((flips_all & ~flips_env).sum()), dRz_bound=Bz,
               min_contact_margin=float(np.min(np.where(contact0, Rz_env_min, -Rz_env_max) / np.maximum(Bz, 1e-30))),
               n_lifted=int((~contact0).sum()), q_max=float(q.max()), mode="ring",
               margin={"uplift": lam_up, "bearing": lam_bear, "sliding": lam_slide})
    # 침하 (감압 LC2 = 가장 무거운 상태): 최대 침하, 이웃 지지점 사이 각변위
    settle = {"max_m": 0.0, "diff_m": 0.0, "angular": 0.0, "dcr": 0.0}
    if U_lc2 is not None:
        sz = -np.asarray(U_lc2, dtype=float)                                   # 아래로 + 
        xy = np.array([(x, y) for x, y, _ in RING_LAYOUT["points"]])
        chord = np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1)
        ang = np.abs(np.roll(sz, -1) - sz) / chord
        settle = {"max_m": float(sz.max()), "diff_m": float(sz.max() - sz.min()), "angular": float(ang.max())}
        settle["dcr"] = max(settle["max_m"] / SETTLEMENT_LIMIT_M if SETTLEMENT_LIMIT_M else 0.0,
                            settle["angular"] / ANGULAR_DISTORTION_LIMIT if ANGULAR_DISTORTION_LIMIT else 0.0)
    out["settle"] = settle
    reasons = [name for name, v in (("들림", out["uplift_dcr"]), ("지지력", out["bearing_dcr"]),
                                    ("활동", out["sliding_dcr"]), ("침하", settle["dcr"])) if v >= 1.0]
    out["ok"] = not reasons
    out["reason"] = "/".join(reasons)
    return out


def assess_ring(model, N_env, N_lc2, N_gen=None):
    """링 축력 판정: 압축 이용률(f_c), 인장이 f_t·A를 넘으면 필요한 철근 질량. N_gen = 지진 상한 생성자 (주방향, 생성자, m)."""
    g = model["ground"]
    if g["mode"] != "ring":
        return {"ok": True, "T_max": 0.0, "C_max": 0.0, "comp_dcr": 0.0, "rebar_kg": 0.0, "need_rebar": False,
                "margin": np.inf}
    A = g["ring_sec"]["A"]
    BN = np.zeros_like(N_lc2) if N_gen is None else quake_bound_abs(N_gen)
    N_hi = np.maximum(N_env.max(axis=0) + BN, N_lc2)
    N_lo = np.minimum(N_env.min(axis=0) - BN, N_lc2)
    T_max, C_max = float(max(N_hi.max(), 0.0)), float(max(-N_lo.min(), 0.0))
    need = T_max > g["ring_ft"] * A
    A_s = T_max / RING_REBAR_FY if need else 0.0
    rebar_kg = A_s * RING_LAYOUT["chord_length"] * REBAR_DENSITY
    comp_dcr = C_max / (g["ring_fc"] * A)
    N_env_min = N_env.min(axis=0)
    lam = _bisect_margin(lambda L: float(max(-(N_env_min - L * BN).min(), 0.0)) / (g["ring_fc"] * A))
    return {"ok": comp_dcr < 1.0, "T_max": T_max, "C_max": C_max, "comp_dcr": comp_dcr,
            "need_rebar": need, "A_s": A_s, "rebar_kg": rebar_kg, "margin": lam}


# ----------------------------------------------------------------------
# 9. Abaqus 스타일 3D solid mesh 시각화 (면=삼각형 분할 Mesh3d + 검은 wireframe)
# ----------------------------------------------------------------------
def _local_axes(p1, p2):
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    axis = p2 - p1
    L = np.linalg.norm(axis)
    if L < 1e-9:
        raise ValueError("부재 길이가 0입니다.")
    ex = axis / L
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(ex, ref)) > 0.99:
        ref = np.array([1.0, 0.0, 0.0])
    ey = np.cross(ref, ex)
    ey /= np.linalg.norm(ey)
    ez = np.cross(ex, ey)
    return ex, ey, ez


_BOX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
_BOX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
_SQ = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]


class _MeshBuf:
    """Mesh3d 정점·삼각형·면값·호버 문구와 검은 모서리선을 모으는 버퍼"""

    def __init__(self):
        self.X, self.Y, self.Z, self.I, self.J, self.K, self.val, self.txt = ([] for _ in range(8))
        self.ex, self.ey, self.ez = [], [], []

    def add(self, verts, tris, value, text, edges=()):
        base = len(self.X)
        for v in verts:
            self.X.append(float(v[0])); self.Y.append(float(v[1])); self.Z.append(float(v[2]))
            self.txt.append(text)
        for a, b, c in tris:
            self.I.append(base + a); self.J.append(base + b); self.K.append(base + c); self.val.append(value)
        for a, b in edges:
            self.ex += [verts[a][0], verts[b][0], None]
            self.ey += [verts[a][1], verts[b][1], None]
            self.ez += [verts[a][2], verts[b][2], None]

    def line(self, pts):
        for a, b in zip(pts[:-1], pts[1:]):
            self.ex += [a[0], b[0], None]; self.ey += [a[1], b[1], None]; self.ez += [a[2], b[2], None]


def _box_tris():
    return [(f[0], f[1], f[2]) for f in _BOX_FACES] + [(f[0], f[2], f[3]) for f in _BOX_FACES]


def shs_tube_mesh(p1, p2, D, t):
    """선정된 정사각 강관(SHS)의 실제 단면으로 부재를 만든다: 바깥 한 변 D, 안쪽 D-2t (속이 빈 관).
    정점 16개 = 양 끝 x (바깥 4 + 안쪽 4). 반환: (정점, 삼각형, 모서리)"""
    _, ey, ez = _local_axes(p1, p2)
    p1, p2 = np.asarray(p1, dtype=float), np.asarray(p2, dtype=float)
    ho, hi = D / 2.0, D / 2.0 - t
    V = np.array([p + sy * h * ey + sz * h * ez for p in (p1, p2) for h in (ho, hi) for sy, sz in _SQ])

    def ix(e, layer, k):
        return e * 8 + layer * 4 + (k % 4)

    tris = []
    for k in range(4):
        a, b, c, d = ix(0, 0, k), ix(0, 0, k + 1), ix(1, 0, k + 1), ix(1, 0, k)      # 바깥 벽
        tris += [(a, b, c), (a, c, d)]
        a, b, c, d = ix(0, 1, k), ix(1, 1, k), ix(1, 1, k + 1), ix(0, 1, k + 1)      # 안쪽 벽
        tris += [(a, b, c), (a, c, d)]
        for e in (0, 1):                                                             # 끝면 (사각 고리)
            a, b, c, d = ix(e, 0, k), ix(e, 0, k + 1), ix(e, 1, k + 1), ix(e, 1, k)
            tris += [(a, b, c), (a, c, d)]
    edges = ([(ix(e, layer, k), ix(e, layer, k + 1)) for e in (0, 1) for layer in (0, 1) for k in range(4)]
             + [(ix(0, 0, k), ix(1, 0, k)) for k in range(4)])
    return V, tris, edges


def _ring_prism(c, a0, a1, r_in, r_out, z_top, z_bot):
    """링 원호 한 구간 (각 a0 -> a1)의 직사각 단면 프리즘 정점 8개 (_BOX_FACES 순서)"""
    def pts(a):
        u = np.array([math.cos(a), math.sin(a), 0.0])
        return [c + r_in * u + [0, 0, z_top], c + r_out * u + [0, 0, z_top],
                c + r_out * u + [0, 0, z_bot], c + r_in * u + [0, 0, z_bot]]
    return np.array(pts(a0) + pts(a1))


def _xanchor(fig, row, col):
    """2D 칸의 x축 이름 ('x', 'x2', ...) — 가로세로 비율 1:1 고정용"""
    return fig.get_subplot(row, col).xaxis.plotly_name.replace("xaxis", "x")


def plot_dome_3d_abaqus(dome, nodes, edges, member_stress, section_info=None, design=None, assessment=None,
                        save_html_path=None, title=None, show=True):
    """Abaqus 스타일 3D 그림 — 해석에서 정해진 치수 그대로
      - 돔 부재: 최소 질량으로 선정된 SHS 단면(바깥 D, 두께 t, 속이 빈 관)으로 그리고 지배 조합 응력으로 색칠
      - 링 기초: 역T 단면(기둥부 폭 RING_STEM_B, 바닥판 폭 B, 바닥판 두께 design.flange_t)을
                 지표(z=0)부터 깊이 D까지 실제 위치에 그림. 돔 밑단 12개 절점이 링 기둥부 위에 강접
      - 원지반: 링이 보이도록 반투명. 레골리스 덮개도 반투명으로 함께 표시"""
    sec_D = section_info["D"] if section_info else 2.0 * dome["c"]
    sec_t = section_info["t"] if section_info else 0.04 * sec_D
    ring = FOUNDATION_MODE == "ring" and design is not None and design.ring_depth > 0
    lay = RING_LAYOUT
    c0, R = np.asarray(lay["center"], dtype=float), lay["R"]

    # ---- 돔 부재 (실제 SHS 단면) ----
    dome_buf = _MeshBuf()
    for ele_tag, ni, nj in edges:
        V, tris, eds = shs_tube_mesh(nodes[ni], nodes[nj], sec_D, sec_t)
        sig = member_stress.get(ele_tag, 0.0) / 1e6
        a = assessment.get(ele_tag) if assessment else None
        text = (f"부재 {member_label(ele_tag)} ({ni}-{nj})<br>SHS □{sec_D*1000:.1f} × {sec_t*1000:.1f} mm<br>"
                f"지배 응력 {sig:+.1f} MPa" + (f"<br>DCR {a['max_dcr']:.3f} ({a['gov_case']})" if a else ""))
        dome_buf.add(V, tris, sig, text, eds)

    from plotly.subplots import make_subplots
    fig = make_subplots(rows=2, cols=2, column_widths=[0.72, 0.28], horizontal_spacing=0.06, vertical_spacing=0.12,
                        specs=[[{"type": "scene", "rowspan": 2}, {"type": "xy"}], [None, {"type": "xy"}]],
                        subplot_titles=("", "선정 부재 단면 (SHS)", "링 기초 단면 (역T)" if ring else ""))
    _orig_add = fig.add_trace
    fig.add_trace = lambda tr, row=1, col=1, **kw: _orig_add(tr, row=row, col=col, **kw)   # 기본은 3D 칸
    cmax = max((abs(v) for v in dome_buf.val), default=1.0) or 1.0
    fig.add_trace(go.Mesh3d(
        x=dome_buf.X, y=dome_buf.Y, z=dome_buf.Z, i=dome_buf.I, j=dome_buf.J, k=dome_buf.K,
        intensity=dome_buf.val, intensitymode="cell", colorscale="RdBu", cmin=-cmax, cmax=cmax,
        opacity=1.0, flatshading=True, showscale=True,
        colorbar=dict(title="응력 (MPa)<br>+인장 / −압축", x=0.70, xanchor="left", y=0.5, yanchor="middle",
                      len=0.6, thickness=14),
        text=dome_buf.txt, hovertemplate="%{text}<extra></extra>",
        name=f"돔 부재 SHS □{sec_D*1000:.0f}×{sec_t*1000:.1f} mm ({MEMBER_MATERIAL})", showlegend=True))
    fig.add_trace(go.Scatter3d(x=dome_buf.ex, y=dome_buf.ey, z=dome_buf.ez, mode="lines",
                               line=dict(color="black", width=1.5), name="부재 모서리", hoverinfo="skip"))

    z_floor = -(design.ring_depth + 2.0) if ring else -0.1 * 2 * R
    extent = R + (design.ring_b / 2.0 if ring else 0.0) + 4.0
    if ring:
        D, B, sb, tf = design.ring_depth, design.ring_b, RING_STEM_B, design.flange_t
        h = D - tf
        ang = []
        for x, y, _bt in lay["points"]:
            ang.append(math.atan2(y - c0[1], x - c0[0]))
        ang = np.unwrap(np.array(ang))
        ring_buf = _MeshBuf()
        rtext = (f"링 기초 (황 콘크리트, 역T 단면)<br>깊이 D = {D:.1f} m (지표 ~ 바닥)<br>"
                 f"바닥판 폭 B = {B:.1f} m, 두께 {tf:.1f} m<br>기둥부 폭 {sb:.1f} m, 높이 {h:.1f} m")
        n = len(ang)
        for i in range(n):
            a0, a1 = ang[i], ang[(i + 1) % n] + (2 * math.pi if i == n - 1 else 0.0)
            for r_in, r_out, zt, zb in ((R - sb / 2, R + sb / 2, 0.0, -h), (R - B / 2, R + B / 2, -h, -D)):
                V = _ring_prism(c0, a0, a1, r_in, r_out, zt, zb)
                ring_buf.add(V, _box_tris(), 0.0, rtext, [(0, 4), (1, 5), (2, 6), (3, 7)])
        # 돔 밑단 12개 절점 위치에서 역T 단면 윤곽 (기초가 내려가는 12개 구역)
        base_xy = []
        for x, y, bt in lay["points"]:
            if bt is None or (DOME_MEMBER_INFO is not None and bt not in DOME_MEMBER_INFO["joints"]):
                continue
            u = np.array([x - c0[0], y - c0[1], 0.0]); u /= np.linalg.norm(u)
            prof = [(R - sb / 2, 0.0), (R + sb / 2, 0.0), (R + sb / 2, -h), (R + B / 2, -h), (R + B / 2, -D),
                    (R - B / 2, -D), (R - B / 2, -h), (R - sb / 2, -h), (R - sb / 2, 0.0)]
            ring_buf.line([c0 + r * u + [0, 0, z] for r, z in prof])
            base_xy.append((x, y, bt))
        fig.add_trace(go.Mesh3d(
            x=ring_buf.X, y=ring_buf.Y, z=ring_buf.Z, i=ring_buf.I, j=ring_buf.J, k=ring_buf.K,
            color="#9a8f84", opacity=1.0, flatshading=True, text=ring_buf.txt, hovertemplate="%{text}<extra></extra>",
            name=f"링 기초 (깊이 {D:.1f} m, 바닥판 폭 {B:.1f} m)", showlegend=True))
        fig.add_trace(go.Scatter3d(x=ring_buf.ex, y=ring_buf.ey, z=ring_buf.ez, mode="lines",
                                   line=dict(color="#3d3832", width=2), name=f"링 단면 윤곽 (돔 밑단 {len(base_xy)}곳)",
                                   hoverinfo="skip"))
        fig.add_trace(go.Scatter3d(
            x=[p[0] for p in base_xy], y=[p[1] for p in base_xy], z=[0.0] * len(base_xy), mode="markers",
            marker=dict(size=5, color="#eb6834", symbol="diamond"), name=f"돔 밑단 접합부 {len(base_xy)}개 (링 기둥부에 강접)",
            text=[f"돔 밑단 절점 {bt}" for _, _, bt in base_xy], hovertemplate="%{text}<extra></extra>"))
    # ---- 원지반 (반투명) ----
    gV = np.array([[c0[0] - extent, c0[1] - extent, 0.0], [c0[0] + extent, c0[1] - extent, 0.0],
                   [c0[0] + extent, c0[1] + extent, 0.0], [c0[0] - extent, c0[1] + extent, 0.0],
                   [c0[0] - extent, c0[1] - extent, z_floor], [c0[0] + extent, c0[1] - extent, z_floor],
                   [c0[0] + extent, c0[1] + extent, z_floor], [c0[0] - extent, c0[1] + extent, z_floor]])
    gt = _box_tris()
    fig.add_trace(go.Mesh3d(x=gV[:, 0], y=gV[:, 1], z=gV[:, 2], i=[a for a, _, _ in gt], j=[b for _, b, _ in gt],
                            k=[cc for _, _, cc in gt], color="#d9c4a0", opacity=0.10, name="원지반 레골리스 (반투명)",
                            hoverinfo="skip", showlegend=True))
    gbuf = _MeshBuf()
    gbuf.add(gV, [], 0.0, "", _BOX_EDGES)
    fig.add_trace(go.Scatter3d(x=gbuf.ex, y=gbuf.ey, z=gbuf.ez, mode="lines", line=dict(color="#b8a684", width=1.5),
                               showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter3d(x=[c0[0] - extent * 0.95], y=[c0[1] - extent * 0.95], z=[0.4], mode="text",
                               text=["지표 z = 0"], textfont=dict(size=12, color="#6b5d45"), showlegend=False,
                               hoverinfo="skip"))

    # ---- 레골리스 덮개 (반투명, 기본 표시) ----
    if design is not None and design.cover_t > 0:
        Rc = DOME_SPHERE_RADIUS_M + design.cover_t
        th = np.linspace(0.0, math.pi / 2, 19)
        ph = np.linspace(0.0, 2 * math.pi, 49)
        CX = [c0[0] + Rc * math.cos(a) * math.cos(b) for a in th for b in ph]
        CY = [c0[1] + Rc * math.cos(a) * math.sin(b) for a in th for b in ph]
        CZ = [Rc * math.sin(a) for a in th for b in ph]
        ci, cj, ck = [], [], []
        m = len(ph)
        for a in range(len(th) - 1):
            for b in range(m - 1):
                p00, p01, p10, p11 = a * m + b, a * m + b + 1, (a + 1) * m + b, (a + 1) * m + b + 1
                ci += [p00, p00]; cj += [p01, p11]; ck += [p11, p10]
        fig.add_trace(go.Mesh3d(x=CX, y=CY, z=CZ, i=ci, j=cj, k=ck, color="#b08d57", opacity=0.18,
                                name=f"레골리스 덮개 t = {design.cover_t:.1f} m",
                                hoverinfo="skip", showlegend=True))

    # ---- 2D: 선정 부재 단면 (실제 치수, mm) ----
    Dm, tm = sec_D * 1000.0, sec_t * 1000.0
    sq = lambda h: ([-h, h, h, -h, -h], [-h, -h, h, h, -h])
    ox, oy = sq(Dm / 2); ix_, iy_ = sq(Dm / 2 - tm)
    fig.add_trace(go.Scatter(x=ox, y=oy, mode="lines", fill="toself", fillcolor="#2a78d6", line=dict(color="#1a1a1a", width=1.5),
                             showlegend=False, hoverinfo="skip"), row=1, col=2)
    fig.add_trace(go.Scatter(x=ix_, y=iy_, mode="lines", fill="toself", fillcolor="white", line=dict(color="#1a1a1a", width=1.5),
                             showlegend=False, hoverinfo="skip"), row=1, col=2)
    A_mm2 = (sec_D ** 2 - (sec_D - 2 * sec_t) ** 2) * 1e6
    fig.add_trace(go.Scatter(
        x=[0.0, 0.0, Dm / 2 - tm / 2], y=[-Dm / 2 - 70, 0.0, Dm / 2 + 55], mode="text",
        text=[f"D = {Dm:.1f} mm", f"A = {A_mm2/1e3:,.1f}×10³ mm²", f"t = {tm:.1f} mm"],
        textfont=dict(size=12, color="#1a1a1a"), showlegend=False, hoverinfo="skip"), row=1, col=2)
    lim = Dm / 2 + 120
    fig.update_xaxes(range=[-lim, lim], visible=False, row=1, col=2)
    fig.update_yaxes(range=[-lim, lim], visible=False, scaleanchor=_xanchor(fig, 1, 2), scaleratio=1, row=1, col=2)

    # ---- 2D: 링 기초 단면 (실제 치수, m). 가로 = 링 중심선에서 반경 방향 거리 ----
    if ring:
        prof_x = [-sb / 2, sb / 2, sb / 2, B / 2, B / 2, -B / 2, -B / 2, -sb / 2, -sb / 2]
        prof_z = [0.0, 0.0, -h, -h, -D, -D, -h, -h, 0.0]
        gx = B / 2 + 1.5
        fig.add_trace(go.Scatter(x=[-gx, gx, gx, -gx, -gx], y=[0.0, 0.0, -D - 1.5, -D - 1.5, 0.0], mode="lines", fill="toself",
                                 fillcolor="rgba(217,196,160,0.35)", line=dict(width=0), showlegend=False,
                                 hoverinfo="skip"), row=2, col=2)
        fig.add_trace(go.Scatter(x=[-gx, gx], y=[0.0, 0.0], mode="lines", line=dict(color="#6b5d45", width=2),
                                 showlegend=False, hoverinfo="skip"), row=2, col=2)
        fig.add_trace(go.Scatter(x=prof_x, y=prof_z, mode="lines", fill="toself", fillcolor="#9a8f84",
                                 line=dict(color="#1a1a1a", width=1.5), showlegend=False, hoverinfo="skip"),
                      row=2, col=2)
        fig.add_trace(go.Scatter(x=[0.0], y=[0.0], mode="markers", marker=dict(size=9, color="#eb6834", symbol="diamond"),
                                 showlegend=False, hoverinfo="skip"), row=2, col=2)
        dim = dict(color="#1a1a1a", width=1.2)
        fig.add_trace(go.Scatter(x=[B / 2 + 0.6] * 2, y=[0.0, -D], mode="lines+markers", line=dim,
                                 marker=dict(size=5, symbol="line-ew-open"), showlegend=False, hoverinfo="skip"),
                      row=2, col=2)
        fig.add_trace(go.Scatter(x=[-B / 2, B / 2], y=[-D - 0.6] * 2, mode="lines+markers", line=dim,
                                 marker=dict(size=5, symbol="line-ns-open"), showlegend=False, hoverinfo="skip"),
                      row=2, col=2)
        fig.add_trace(go.Scatter(
            x=[B / 2 + 0.75, 0.0, sb / 2 + 0.25, -B / 2 - 0.15, 0.0], y=[-D / 2, -D - 1.05, -h / 2, -h - tf / 2, 0.45],
            mode="text", textposition=["middle right", "middle center", "middle right", "middle left", "middle center"],
            text=[f"깊이 D = {D:.1f} m", f"바닥판 폭 B = {B:.1f} m", f"기둥부 {sb:.1f} m", f"바닥판 두께 {tf:.1f} m",
                  "지표 z = 0 · 돔 밑단 절점"],
            textfont=dict(size=11, color="#1a1a1a"), showlegend=False, hoverinfo="skip"), row=2, col=2)
        fig.update_xaxes(range=[-gx - 2.6, gx + 2.6], visible=False, row=2, col=2)
        fig.update_yaxes(range=[-D - 1.5, 1.0], visible=False, scaleanchor=_xanchor(fig, 2, 2), scaleratio=1,
                         row=2, col=2)

    if title is None:
        title = (f"{DOME_LABEL} — 부재 SHS □{sec_D*1000:.1f} × {sec_t*1000:.1f} mm"
                 + (f" · 링 기초 깊이 {design.ring_depth:.1f} m, 바닥판 폭 {design.ring_b:.1f} m" if ring else "")
                 + " · 지배 조합 응력"
                 + f"<br><sup>돔 부재 재료: {MEMBER_MATERIAL_LABEL}</sup>")
    zmin = z_floor
    fig.add_trace = _orig_add
    fig.update_layout(
        title=title, template="plotly_white", width=1500, height=950,
        scene=dict(xaxis_title="X 북 (m)", yaxis_title="Y 동 (m)", zaxis_title="Z (m)", aspectmode="data",
                   zaxis=dict(range=[zmin - 1.0, DOME_SPHERE_RADIUS_M + (design.cover_t if design else 0) + 1.0]),
                   camera=dict(eye=dict(x=1.5, y=-1.5, z=0.55), center=dict(x=0, y=0, z=-0.15))),
        legend=dict(x=0.0, y=0.99, bgcolor="rgba(255,255,255,0.8)"))
    if save_html_path:
        out_path = Path(save_html_path).resolve()
        fig.write_html(str(out_path), auto_open=False)
        print(f"[viz] 저장: {out_path} (부재 SHS □{sec_D*1000:.1f}×{sec_t*1000:.1f} mm"
              + (f", 링 깊이 {design.ring_depth:.1f} m · 폭 {design.ring_b:.1f} m" if ring else "") + ")")
    if show:
        fig.show()
    return fig


# ----------------------------------------------------------------------
# 9-3. 그림: 한글 폰트 / 환경 타임라인 / 탐색 결과
# ----------------------------------------------------------------------
def _setup_korean_font(plt):
    from matplotlib import font_manager
    avail = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Malgun Gothic", "AppleGothic", "NanumGothic", "Noto Sans CJK KR"):
        if name in avail:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams["axes.unicode_minus"] = False


def plot_env_timeline(env, env_dcr_t=None, quake_lmst_h=None, out_path=None):
    """기온 / 외기압·순압 / 바람 / (LC1 환경 단독) 최대 DCR 을 LMST 축으로 그린다."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _setup_korean_font(plt)

    x = env["lmst_h"]
    n = 4 if env_dcr_t is not None else 3
    fig, axes = plt.subplots(n, 1, figsize=(11, 2.4 * n + 0.6), sharex=True)

    ax = axes[0]
    ax.plot(x, env["T_K"], color="tab:red", lw=1.4)
    ax.axhline(env["T_ref"], color="0.4", ls="--", lw=0.9)
    ax.text(x[0], env["T_ref"], f" 무응력 기준 {env['T_ref']:.0f} K", va="bottom", fontsize=8, color="0.3")
    ax.set_ylabel("기온 (K)")

    ax = axes[1]
    ax.plot(x, env["P_Pa"], color="tab:blue", lw=1.4)
    ax.set_ylabel("외기압 (Pa)")

    ax = axes[2]
    ax.plot(x, env["V_gust"], color="tab:green", lw=1.2, label="순간최대풍속")
    ax.plot(x, env["V_mean"], color="tab:green", lw=0.8, alpha=0.5, label="평균풍속")
    ax.set_ylabel("풍속 (m/s)")
    ax2 = ax.twinx()
    ax2.plot(x, env["q_Pa"], color="tab:purple", lw=0.9, ls=":", label="풍압 q")
    ax2.set_ylabel("풍압 (Pa)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    if env_dcr_t is not None:
        ax = axes[3]
        ax.plot(x, env_dcr_t, color="black", lw=1.3)
        ax.set_ylabel("LC1 환경 단독\n최대 DCR")

    if quake_lmst_h is not None:
        hod = quake_lmst_h % 24.0
        for xq in (hod + 24.0 * k for k in range(-1, 3)):
            if x[0] <= xq <= x[-1]:
                for a in axes:
                    a.axvline(xq, color="tab:orange", lw=1.0, alpha=0.8)
                axes[0].text(xq, axes[0].get_ylim()[1], " S1222a 지방시", color="tab:orange",
                             va="top", fontsize=8)

    ticks = np.arange(math.ceil(x[0] / 3) * 3, x[-1] + 0.01, 3)
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels([f"{int(t % 24):02d}:00" for t in ticks])
    axes[-1].set_xlabel(f"LMST (Sol {env['meta']['sol0']} -> {env['meta']['sol0'] + 1})")
    for a in axes:
        a.grid(alpha=0.3)
    fig.suptitle("InSight 실측 환경하중 타임라인 (TWINS + PS)", fontsize=12)
    fig.tight_layout()
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150)
        print(f"[env] 그림 저장: {Path(out_path).resolve()}")
    plt.close(fig)


def plot_search_results(df, out_path):
    """D x B 격자 결과: (1) B에 따른 돔 재료, (2) B에 따른 덮개·바닥판 두께, (3) 돔 재료 vs 기초 재료 (파레토)"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _setup_korean_font(plt)

    ok = df[df["feasible"]].sort_values(["D_m", "B_m"])
    if ok.empty:
        print("[search] 성립하는 칸이 없어 그림을 생략합니다.")
        return
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, max(ok["D_m"].nunique(), 1)))
    for c, (D, g) in zip(colors, ok.groupby("D_m")):
        axes[0].plot(g["B_m"], g["material_ton"], "o-", color=c, label=f"D = {D:.0f} m")
        axes[1].plot(g["B_m"], g["t_m"], "o-", color=c, label=f"덮개 (D = {D:.0f} m)")
        axes[1].plot(g["B_m"], g["flange_t_m"], "s--", color=c, label=f"바닥판 (D = {D:.0f} m)")
        axes[2].scatter(g["ring_concrete_ton"], g["material_ton"], color=c, s=26, alpha=0.6, label=f"D = {D:.0f} m")
    par = df[df.get("pareto", False) == True].sort_values("ring_concrete_ton") if "pareto" in df else df.iloc[0:0]
    if not par.empty:
        axes[2].plot(par["ring_concrete_ton"], par["material_ton"], "-", color="tab:red", lw=1.2, label="파레토 해")
    marks = {"material_min": ("*", "tab:red"), "concrete_min": ("D", "tab:blue")}
    if "end" in df:
        for name, (mk, col) in marks.items():
            if name not in FINAL_ENDS:
                continue
            r = df[df["end"].str.contains(name, na=False)]
            for ax, x in ((axes[0], "B_m"), (axes[2], "ring_concrete_ton")):
                ax.scatter(r[x], r["material_ton"], s=180, marker=mk, color=col, zorder=5, label=END_LABELS[name])
    axes[0].set_xlabel("바닥판 폭 B (m)"); axes[0].set_ylabel("돔 재료 (t) = 돔 부재 + 철근")
    axes[1].set_xlabel("바닥판 폭 B (m)"); axes[1].set_ylabel("두께 (m)")
    axes[2].set_xlabel("기초 재료 = 링 황 콘크리트 (t)"); axes[2].set_ylabel("돔 재료 (t)")
    axes[0].set_title("B를 넓히면 돔 재료는?", fontsize=10)
    axes[1].set_title("B를 넓히면 덮개는 얇아지고 바닥판은 두꺼워짐", fontsize=10)
    axes[2].set_title("빨간 선 = 둘 다 더 줄일 수 없는 설계 (파레토 해)", fontsize=10)
    for ax in axes:
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle(f"링 깊이 D x 바닥판 폭 B 격자 탐색 ({DOME_LABEL}, 내압 {INTERNAL_PRESSURE_PA/1000:.0f} kPa, "
                 f"링 황 콘크리트 '{RING_MATERIAL}' f_c {RING_MATERIAL_CASES[RING_MATERIAL][0]/1e6:.0f} MPa, "
                 f"바닥판 {'무근' if FLANGE_MODE == 'plain' else '무근/철근 자동'})", fontsize=11)
    fig.text(0.995, 0.005, f"돔 부재 재료: {MEMBER_MATERIAL_LABEL}", ha="right", va="bottom", fontsize=7, color="0.4")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[search] 그림 저장: {Path(out_path).resolve()}")


def _validate_quake_starts(starts, quake_len, duration):
    starts = sorted(float(s) for s in starts)
    if starts and (starts[0] < 0.0 or starts[-1] + quake_len > duration):
        raise ValueError(f"지진 구간이 타임라인(0~{duration:.0f}s)을 벗어납니다. "
                         f"(마지막 지진 {starts[-1]:.0f}~{starts[-1] + quake_len:.0f}s)")
    for a, b in zip(starts[:-1], starts[1:]):
        if a + quake_len > b:
            raise ValueError(f"지진 구간이 겹칩니다: {a:.0f}~{a + quake_len:.0f}s 와 {b:.0f}s 시작.")
    return starts


def save_sol_timeline(path, env, F_env, gens=None, th=None, quake_starts=SOL_QUAKE_STARTS_S):
    """npz(zip) 구조
        t_env, lmst_env, T_K, P_Pa, q_Pa, dT, dp : 환경 타임라인 (dT, q_Pa는 덮개 반영 후 값)
        env_ele_<tag>   (시각수, 18)              : 환경 부재력 [I단 6 | J단 6 | 중앙 단면 6]
        eqs_gen_ele_<tag> (주방향, 생성자, 18)     : 지진 상한 생성자 부재력. 지진 증분은 Σ s_g g (|s_g| <= 1) 안에 있다
        eqs_gen_labels, eqs_gen_primary          : 생성자 이름(mode1.., residual_X..), 주방향
      th(시간이력 검증을 했을 때만): {"force_paths": ..., "F_base": ...}
        t_quake, quake_starts_s                  : 지진 증분 시간축, 지진 시작 시각
        quake_ele_<tag> (지진스텝, 12) float32    : 실측 시간이력에 의한 부재력 증분
    복원은 reconstruct_member_history() 사용."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    starts = []
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        def put(name, arr):
            with zf.open(f"{name}.npy", "w", force_zip64=True) as f:
                np.lib.format.write_array(f, np.asarray(arr))

        put("t_env", env["t_s"]); put("lmst_env", env["lmst_h"])
        put("member_material", json.dumps(MEMBER_MATERIAL_INFO, ensure_ascii=False))
        put("dome", json.dumps(DOME_INFO, ensure_ascii=False))
        for k in ("T_K", "P_Pa", "q_Pa", "dT", "dp"):
            put(k, env[k])
        for e_idx, ele in enumerate(ELE_TAGS):
            put(f"env_ele_{ele}", F_env[:, e_idx, :])
            if gens is not None:
                put(f"eqs_gen_ele_{ele}", gens["F"][:, :, e_idx, :])
        if gens is not None:
            put("eqs_gen_labels", np.array(gens["labels"])); put("eqs_gen_primary", np.array(gens["primary"]))
        if th is not None:
            t_quake = _read_table(next(iter(th["force_paths"].values())))[:, 0]
            starts = _validate_quake_starts(quake_starts, float(t_quake[-1]), env["duration_s"])
            put("t_quake", t_quake); put("quake_starts_s", np.array(starts))
            for e_idx, ele in enumerate(ELE_TAGS):
                put(f"quake_ele_{ele}", (_read_table(th["force_paths"][ele])[:, 1:13]
                                         - th["F_base"][e_idx, :12]).astype(np.float32))

    s0 = env["meta"]["sol0"]
    if starts:
        spans = ", ".join(f"{format_lmst(np.interp(s, env['t_s'], env['lmst_h']), s0)}" for s in starts)
        print(f"[sol] 1솔 부재력 타임라인 저장: {path.resolve()} (시간이력 지진 {len(starts)}회: {spans})")
    else:
        print(f"[sol] 1솔 부재력 타임라인 저장: {path.resolve()} (환경 10분 + 지진 등가정적 상한 생성자)")
    return path


def reconstruct_member_history(npz_path, ele_tag):
    """저장된 타임라인에서 부재 하나의 전체 이력 (t, 부재력)을 복원한다 (시간이력 지진이 있으면 I·J단 12열).
    시간이력 지진이 저장돼 있으면 그 구간은 환경(선형보간) + 지진 증분, 없으면 환경 부재력만 반환한다
    (등가정적 지진은 시각이 없는 '상한'이므로 이력에 넣지 않는다 -> 피로 상한은 fatigue_upper_bound 참조)."""
    d = np.load(npz_path)
    t_env, F_env = d["t_env"], d[f"env_ele_{ele_tag}"]
    if "t_quake" not in d.files:
        return t_env, F_env
    t_q, dF = d["t_quake"], d[f"quake_ele_{ele_tag}"].astype(float)
    F_env = F_env[:, :dF.shape[1]]

    t_parts, f_parts, cursor = [], [], -np.inf
    for st in d["quake_starts_s"]:
        m = (t_env > cursor) & (t_env < st)
        t_parts.append(t_env[m]); f_parts.append(F_env[m])
        tq = st + t_q
        env_at_q = np.column_stack([np.interp(tq, t_env, F_env[:, k]) for k in range(F_env.shape[1])])
        t_parts.append(tq); f_parts.append(env_at_q + dF)
        cursor = tq[-1]
    m = t_env > cursor
    t_parts.append(t_env[m]); f_parts.append(F_env[m])
    return np.concatenate(t_parts), np.vstack(f_parts)


# ----------------------------------------------------------------------
# 11. 실행 파이프라인 — 질량 하나에 대해 해석하고 파손 여부까지 반환
# ----------------------------------------------------------------------
def assemble_masses(model, design, verbose=True):
    """절점 질량(강재 + 링 + 막 [+ 덮개]) 과 관성 없는 추가 무게(덮개, COVER_INERTIA=False일 때)"""
    mass = dict(model["mass_dict"])
    parts = {"steel": model["steel_mass_kg"], "ring": model["ground"]["ring_mass_kg"], "membrane": 0.0, "cover": 0.0}
    if USE_MEMBRANE:
        mem = compute_membrane_node_mass(model["node_tags"], dome_NODES, DOME_EDGES, verbose=verbose)
        for k, m in mem.items():
            mass[k] = mass.get(k, 0.0) + m
        parts["membrane"] = sum(mem.values())
    cover = compute_cover_node_mass(model["node_tags"], dome_NODES, DOME_EDGES, design.cover_t, verbose=verbose)
    parts["cover"] = sum(cover.values())
    extra = {}
    if COVER_INERTIA:
        for k, m in cover.items():
            mass[k] = mass.get(k, 0.0) + m
    else:
        extra = cover
    return mass, extra, parts, cover


def _failed_static_case(target_mass_ton, design, status, reason, section_info):
    return {"target_mass_ton": float(target_mass_ton), "design": design, "status": status, "reason": reason,
            "members_ok": status != "member_fail", "foundation_ok": status != "foundation_fail",
            "disp_ok": True, "passed": False, "failed_members": [], "failed_count": 0,
            "max_dcr": np.inf, "max_dcr_env_only": np.inf, "disp": None, "section_info": section_info}


def _classify_static_failure(model):
    """정적해석 불수렴 원인 추정: 어느 지지점이든 흙 인발 한계에 도달했으면 기초(들림), 아니면 부재/수치"""
    g = model["ground"]
    if g["mode"] != "ring":
        return "member_fail", "정적해석 불수렴"
    R = ground_reactions(model)
    if np.any(-R[:, 2] >= 0.98 * g["cap_up"]):
        return "foundation_fail", "들림: 플랜지 위 흙 인발 한계 도달 (정적해석 불수렴)"
    if np.hypot(R[:, 0].sum(), R[:, 1].sum()) >= 0.98 * g["cap_h"].sum():
        return "foundation_fail", "활동: 링 전체 수평 저항 한계 도달 (정적해석 불수렴)"
    return "member_fail", "정적해석 불수렴 (부재 불안정 또는 수치 문제)"


def build_lc1_state(section, design, env, verbose=True, label_prefix=""):
    """비선형 모델(PDelta + ENT/인발한계 지반)에서 1단계 LC2(자중+막+덮개+링), 2단계 LC1(+평균 순압).
    LC1이 '내부기압이 만든 기준상태'이고, 지반 접촉 상태(눌림/들림)가 여기서 정해진다.
    끝나면 모델은 LC1 상태로 남아 있다 (호출한 쪽이 ops.wipe() 책임). 실패하면 status/reason을 돌려준다."""
    g = planet_g()
    model = build_structure(section, design, verbose=verbose)
    masses, extra_weight, mass_parts = apply_dead_loads(model, design, verbose=verbose)
    st = {"model": model, "masses": masses, "extra_weight": extra_weight, "mass_parts": mass_parts, "ok": False}
    if run_static_robust(verbose=verbose, label=f"{label_prefix}LC2 (자중+덮개+링, 내압 0)") != 0:
        st["status"], reason = _classify_static_failure(model)
        st["reason"] = "LC2: " + reason
        return st
    st.update(F_lc2=element_local_forces(), N_lc2=ring_axial(model), R_lc2=ground_reactions(model))
    if model["ground"]["mode"] == "ring":
        st["U_lc2"] = np.array([ops.nodeDisp(int(t), 3) for t in model["ground"]["ring_nodes"]])   # 링 지지점 연직 변위
    st["W_total"] = (sum(masses.values()) + sum(extra_weight.values())) * g

    F_up = 0.0
    if USE_INTERNAL_PRESSURE:
        F_up = DOME_UPLIFT_PER_PA * env["dp_ref"]
        apply_shell_loads(model, pressure_panel_forces(env["dp_ref"]),
                          pressure_node_forces(dome_NODES, DOME_EDGES, env["dp_ref"]), pattern_tag=4)
        if verbose:
            print(f"[pressure] 기준 순압 {env['dp_ref']:,.0f} Pa -> 들림 {F_up/1e6:,.1f} MN "
                  f"(무게 {st['W_total']/1e6:,.1f} MN)")
        if run_static_robust(verbose=verbose, label=f"{label_prefix}LC1 (평균 순압 가압)") != 0:
            st["status"], reason = _classify_static_failure(model)
            st["reason"] = "LC1: " + reason
            return st
    st.update(F_base=element_local_forces(), N_base=ring_axial(model), R_base=ground_reactions(model),
              F_up=F_up, ok=True)
    st["contact"] = st["R_base"][:, 2] > 0.0
    return st


def run_case(target_mass_ton, design, seismic, env, include_quake=True, quake_actual_idx=None, verbose=True):
    """질량 하나 x 설계(D, t, B) 하나의 전체 해석.
    seismic: build_seismic_input() 결과(S1222a 원본 응답스펙트럼, 배율 없음), env: build_env_timeline() 결과(덮개 반영 전).
    지진은 동적해석 대신 등가정적하중(단위가속도 x Sa)으로 넣으므로 해석 시간이 지진 길이와 무관하다."""
    if target_mass_ton <= 0.0:
        raise ValueError("target_mass_ton은 0보다 커야 합니다.")
    section, section_info = make_section_for_target_mass(target_mass_ton)
    env_e = effective_env(env, design)

    # ---- A) 비선형 모델: LC2 감압 -> LC1 가압 기준상태 ----
    st = build_lc1_state(section, design, env, verbose=verbose)
    if not st["ok"]:
        ops.wipe()
        return _failed_static_case(target_mass_ton, design, st["status"], st["reason"], section_info)
    model, masses = st["model"], st["masses"]
    F_lc2, N_lc2, R_lc2 = st["F_lc2"], st["N_lc2"], st["R_lc2"]
    F_base, N_base, R_base = st["F_base"], st["N_base"], st["R_base"]
    W_total, F_up, contact = st["W_total"], st["F_up"], st["contact"]
    eq_err = abs(R_base[:, 2].sum() - (W_total - F_up)) / max(W_total, F_up, 1.0)
    if verbose:
        print(f"[equilibrium] 지반 반력 합 {R_base[:, 2].sum()/1e6:,.2f} MN = 무게 {W_total/1e6:,.2f} - "
              f"들림 {F_up/1e6:,.2f} MN (오차 {eq_err:.1e}) | 들린 지지점 {int((~contact).sum())}/{len(contact)}")
    ops.wipe()

    # ---- B) 단위하중 (LC1 기준상태 접선강성의 선형 모델): 환경 + 지진(질량비례, 모드별) ----
    inertia = masses if include_quake else None
    unit = compute_unit_load_responses(section, design, contact, env_e, node_mass=inertia,
                                       f_rigid_hz=SEISMIC_RIGID_FACTOR * seismic["nyquist_hz"] if include_quake else None)

    # ---- C) 환경 타임라인 (선형 중첩) ----
    F_env = combine_env(F_base, unit["F"], env_e)
    R_env = combine_env(R_base, unit["R"], env_e)
    N_env = combine_env(N_base, unit["N"], env_e)

    # ---- D) 지진 등가정적 상한 (모드 절대합 + 잔여질량, 100/30/30) ----
    sa = dict(seismic["Sa"]) if include_quake else None
    gens = build_quake_generators(unit, seismic, sa) if include_quake else None
    Q_gen = gens["F"] if gens else None
    R_gen = gens["R"] if gens else None
    N_gen = gens["N"] if gens else None

    # ---- E) 판정 ----
    assessment, env_dcr_t = assess_members(
        F_env, F_lc2, Q_gen, model["section"], model["c"],
        get_member_lengths(dome_NODES, DOME_EDGES), env_e, actual_idx=quake_actual_idx)
    ground = assess_ground(model, R_env, R_lc2, R_base, R_gen, U_lc2=st.get("U_lc2"))
    ring = assess_ring(model, N_env, N_lc2, N_gen)

    # ---- F) 자기점검: 모드 질량, 지진 변위, 내부 체적·기압 변화, 지진/내부기압 하중비 ----
    checks, freqs = {}, []
    if include_quake:
        md = unit["modal"]
        order = unit["node_order"]
        idx_long = int(np.argmax(md["T"]))
        T1 = float(md["T"][idx_long])
        freqs = sorted((1.0 / md["T"]).tolist())[:EIGEN_MODES]
        D_gen = gens["D"]
        dV_gen = np.array([[volume_change_m3(D_gen[p, k], order) for k in range(D_gen.shape[1])]
                           for p in range(D_gen.shape[0])])
        dV_bound = float(np.abs(dV_gen).sum(axis=1).max())
        Bd = quake_bound_abs(D_gen)                                           # (절점, 3) 성분별 변위 상한
        X0 = np.array([dome_NODES[k] for k in order], dtype=float)
        span = max(np.ptp(X0[:, 0]), np.ptp(X0[:, 1]))
        h_disp = float(np.hypot(Bd[:, 0], Bd[:, 1]).max())
        m_inertia = md["M_total"]
        mass_sum = md["meff_ratio"].sum(axis=0)
        checks = {
            "Sa": sa, "T1": T1, "freqs": freqs, "n_modes": len(md["T"]),
            "mass_ratio": {d: float(mass_sum[d - 1]) for d in (1, 2, 3)},
            "Sa_T1": {d: spectral_value(seismic["comps"][d]["spec"], T1) for d in (1, 2, 3)},
            "inertia_mass_kg": m_inertia,
            "V_h_N": m_inertia * sa[1], "V_v_N": m_inertia * sa[3],
            "ratio_Vh_Fup": (m_inertia * sa[1]) / F_up if F_up > 0 else np.inf,
            "ratio_Vh_W": (m_inertia * sa[1]) / W_total,
            "dV_max_m3": dV_bound,
            "dp_internal_max_Pa": AIR_GAMMA * INTERNAL_PRESSURE_PA * dV_bound / DOME_ENCLOSED_VOLUME_M3,
            "h_disp_max_m": h_disp, "h_disp_ratio": h_disp / span, "span_m": span,
            "f_max_mode_hz": float(1.0 / md["T"].min()),
            "f1_below_nyquist": (1.0 / T1) < seismic["nyquist_hz"],
        }
        if verbose and min(mass_sum) < SEISMIC_MASS_TARGET:
            print(f"[warn] 모드 유효질량 합이 {SEISMIC_MASS_TARGET:.0%} 미만 (X {mass_sum[0]:.1%}, Y {mass_sum[1]:.1%}, "
                  f"Z {mass_sum[2]:.1%}) — 나머지는 잔여질량 하중으로 포함됨. SEISMIC_MODE_STEPS를 늘리면 더 정확")
        if verbose and checks["f_max_mode_hz"] < SEISMIC_RIGID_FACTOR * seismic["nyquist_hz"]:
            print(f"[warn] 마지막 모드 {checks['f_max_mode_hz']:.1f} Hz < {SEISMIC_RIGID_FACTOR * seismic['nyquist_hz']:.0f} Hz "
                  f"— 잔여질량 하중의 강체 가정이 덜 엄밀함 (잔여질량에도 Sa_max를 써서 보수적으로 처리)")
    disp_ok = (not include_quake) or DISP_LIMIT_RATIO is None or checks["h_disp_ratio"] <= DISP_LIMIT_RATIO

    failed_members = [k for k, v in assessment.items() if v["failed"]]
    gov_ele = max(assessment, key=lambda k: assessment[k]["max_dcr"])
    gov = assessment[gov_ele]
    max_dcr = gov["max_dcr"]
    members_ok = len(failed_members) == 0 and max_dcr < 1.0
    foundation_ok = ground["ok"] and ring["ok"]
    reasons = []
    if not ground["ok"]:
        reasons.append(f"지반 {ground['reason']}")
    if not ring["ok"]:
        reasons.append(f"링 압축 이용률 {ring['comp_dcr']:.2f}")
    actual_vals = [v["max_dcr_actual_time"] for v in assessment.values() if v["max_dcr_actual_time"] is not None]
    margin = None
    if include_quake:
        lam_m = min(v["quake_margin"] for v in assessment.values())
        margin = {"members": lam_m, **ground.get("margin", {}), "ring": ring.get("margin", np.inf)}
        margin["overall"] = min(margin.values())
    s0 = env["meta"]["sol0"]
    if verbose:
        where = format_lmst(gov["gov_lmst_h"], s0) if gov["gov_case"].startswith("LC1") else "-"
        print(f"[check] 부재 최대 DCR {max_dcr:.4f} (부재 {gov_ele}, {gov['gov_case']}, {where}) | "
              f"LC1 환경 단독 {max(v['max_dcr_env_only'] for v in assessment.values()):.4f} | "
              f"LC2 감압 {max(v['max_dcr_lc2'] for v in assessment.values()):.4f}"
              + (f" | 실제 지방시 지진 {max(actual_vals):.4f}" if actual_vals else ""))
        if ground["mode"] == "ring":
            print(f"[check] 지반 DCR: 들림 {ground['uplift_dcr']:.3f}, 지지력 {ground['bearing_dcr']:.3f}, "
                  f"활동 {ground['sliding_dcr']:.3f} | 접촉상태 변화 {ground['n_flip']}곳"
                  f" (지진으로 추가 {ground['n_flip_quake']}곳), 국부 수평 한계 {ground['n_h_limit']}곳 | "
                  f"링 압축 이용률 {ring['comp_dcr']:.3f}, 최대 인장 {ring['T_max']/1e6:.2f} MN"
                  f"{' -> 철근 ' + format(ring['rebar_kg']/1000, '.1f') + ' ton' if ring['need_rebar'] else ''}")
        if include_quake:
            print(f"[check] 지진 상한 (모드 {checks['n_modes']}개 + 잔여질량) | 질량 x Sa = {checks['V_h_N']:,.3g} N = "
                  f"내부기압 들림의 {checks['ratio_Vh_Fup']:.2e}배 | 지진 여유배율 λ = {margin['overall']:.3g}")

    status = "ok" if (members_ok and foundation_ok and disp_ok) else (
        "member_fail" if not (members_ok and disp_ok) else "foundation_fail")
    return {
        "target_mass_ton": float(target_mass_ton),
        "design": design,
        "status": status,
        "reason": "; ".join(reasons),
        "actual_steel_mass_ton": model["steel_mass_kg"] / 1000.0,
        "mass_parts_kg": st["mass_parts"],
        "section": section,
        "section_info": section_info,
        "dome": model,
        "freqs": freqs,
        "unit": unit,
        "env_e": env_e,
        "F_base": F_base,
        "F_lc2": F_lc2,
        "F_env": F_env,
        "R_base": R_base,
        "R_lc2": R_lc2,
        "N_base": N_base,
        "Q_gen": Q_gen,
        "quake_gens": gens,
        "node_mass": masses,
        "seismic_checks": checks,
        "quake_margin": margin,
        "W_total_N": W_total,
        "F_up_N": F_up,
        "equilibrium_error": eq_err,
        "env_dcr_t": env_dcr_t,
        "assessment": assessment,
        "ground": ground,
        "ring": ring,
        "member_stress": {k: v["signed_stress"] for k, v in assessment.items()},
        "failed_members": failed_members,
        "failed_count": len(failed_members),
        "max_dcr": max_dcr,
        "max_dcr_env_only": max(v["max_dcr_env_only"] for v in assessment.values()),
        "max_dcr_lc2": max(v["max_dcr_lc2"] for v in assessment.values()),
        "max_dcr_actual_time": max(actual_vals) if actual_vals else None,
        "gov_ele": gov_ele,
        "members_ok": members_ok,
        "foundation_ok": foundation_ok,
        "passed": status == "ok",
        "disp": None,
        "disp_ok": disp_ok,
        "include_quake": include_quake,
    }


def _apply_env_step(model, env_e, k, dp, EA_aT):
    """검증용: 시각 k의 순압 dp·바람(외피 하중, 패턴 7)과 온도(절점, 패턴 5)를 건다. 건 하중이 있으면 True"""
    shell, legacy = [], []
    if USE_INTERNAL_PRESSURE:
        shell += pressure_panel_forces(dp)
        legacy.append(pressure_node_forces(dome_NODES, DOME_EDGES, dp))
    if env_e["wind_active"]:
        az = env_e["sector"][k] * 360.0 / WIND_SECTORS
        shell += wind_panel_forces(env_e["q_Pa"][k], az)
        legacy.append(wind_node_forces(dome_NODES, DOME_EDGES, env_e["q_Pa"][k], az))
    if shell:
        apply_shell_loads(model, shell, _sum_forces(*legacy), pattern_tag=7)
    if USE_THERMAL:
        apply_nodal_forces(model["node_tags"], thermal_node_forces(dome_NODES, DOME_EDGES, EA_aT), pattern_tag=5)
    return bool(shell) or USE_THERMAL


def verify_superposition(case, env, k=None):
    """선형 중첩 검증: LC1 지배 시각 k의 실측 하중을 비선형 지반 모델에 직접 걸어 푼 부재력과
    '기준 + 단위하중 x 실측값'으로 계산한 부재력을 비교한다. (PDelta + ENT/인발한계 지반 포함)"""
    design = case["design"]
    env_e = effective_env(env, design)
    if k is None:
        k = int(np.argmax(case["env_dcr_t"]))
    model = build_structure(case["section"], design, verbose=False)
    apply_dead_loads(model, design, verbose=False)
    if run_static_robust(verbose=False) != 0:
        raise RuntimeError("검증용 직접해석(LC2 단계)이 수렴하지 않았습니다.")

    EA_aT = case["unit"]["EA_alpha"] * env_e["dT"][k] if USE_THERMAL else 0.0
    if _apply_env_step(model, env_e, k, dp=env_e["dp"][k], EA_aT=EA_aT):
        if run_static_robust(verbose=False) != 0:
            raise RuntimeError("검증용 직접해석(LC1 단계)이 수렴하지 않았습니다.")
    F_direct = correct_thermal_forces(element_local_forces(), EA_aT)
    ops.wipe()

    F_sup = case["F_env"][k]
    err = float(np.abs(F_direct - F_sup).max() / np.abs(F_direct).max())
    print(f"[verify] 선형 중첩 vs 비선형 직접해석 ({format_lmst(env['lmst_h'][k], env['meta']['sol0'])}): "
          f"최대 상대오차 {err:.2e}")
    return err


def quake_vertex(case, p, signs):
    """지진 상한의 한 꼭짓점(주방향 p, 생성자 부호 signs)에 해당하는 절점 하중과 선형 예측 응답.
    모드 생성자 하중 = 계수 x m·φ_i,  잔여질량 생성자 하중 = 계수 x (m·r_d - Σ_i Γ_id m·φ_i)"""
    md, gens = case["unit"]["modal"], case["quake_gens"]
    coef = gens["coef"][p]
    m, phi, G = md["m"], md["phi"], md["gamma"]
    n = len(md["T"])
    mphi = m[None, :, None] * phi                                   # (n, 절점, 3)
    load = np.tensordot(signs[:n] * coef[:n], mphi, axes=(0, 0))
    for d in range(3):
        res = -np.tensordot(G[:, d], mphi, axes=(0, 0))
        res[:, d] += m
        load += signs[n + d] * coef[n + d] * res
    pred = {k: np.tensordot(signs, gens[k][p], axes=(0, 0)) for k in ("F", "R")}
    return dict(zip(md["tags"], load)), pred


def verify_quake_linearity(case, env, scale=None):
    """지진 선형 중첩 검증: 지배 부재·시각의 '환경 + 지진 상한 꼭짓점'을 비선형 지반 + PDelta 모델에 직접 걸어
    지진에 의한 부재력 증분을 구하고, 선형 단위하중 결과와 비교한다.
    지진이 너무 작아 수치오차에 묻히지 않도록 지진을 scale배로 키워 건 뒤 1/scale로 되돌린다
    (scale배에서도 선형이면 1배에서는 당연히 선형). 접촉 상태가 바뀌는지도 함께 본다."""
    if case["quake_gens"] is None:
        return None
    design, env_e = case["design"], effective_env(env, case["design"])
    gov_ele = case["gov_ele"]
    gov = case["assessment"][gov_ele]
    k = gov["gov_env_index"]
    e = ELE_TAGS.index(gov_ele)
    A, E, G, J, Iy, Iz = case["dome"]["section"]
    c = case["dome"]["c"]
    planes = _dcr_planes(MEMBER_FY * A, gov["compression_capacity"], MEMBER_FY * Iy / c, MEMBER_FY * Iz / c)
    gen = case["quake_gens"]["F"][:, :, e, :]
    S_all = [_end_forces(gp, gov["gov_end_lc1"]) @ planes.T for gp in gen]          # 주방향별 (G, 8)
    p, sidx = max(((pp, ss) for pp in range(len(S_all)) for ss in range(8)),
                  key=lambda t: np.abs(S_all[t[0]][:, t[1]]).sum())
    signs = np.where(S_all[p][:, sidx] >= 0.0, 1.0, -1.0)
    lam = case["quake_margin"]["overall"] if case["quake_margin"] else np.inf
    if scale is None:
        scale = float(np.clip(0.5 * lam if np.isfinite(lam) else 1.0e3, 1.0, 1.0e3))
    load, pred = quake_vertex(case, p, signs)

    st = build_lc1_state(case["section"], design, env, verbose=False)
    if not st["ok"]:
        ops.wipe()
        raise RuntimeError("지진 선형성 검증용 직접해석(LC1)이 수렴하지 않았습니다.")
    model = st["model"]
    EA_aT = case["unit"]["EA_alpha"] * env_e["dT"][k] if USE_THERMAL else 0.0
    if _apply_env_step(model, env_e, k, dp=env_e["dp"][k] - env["dp_ref"], EA_aT=EA_aT):
        if run_static_robust(verbose=False, tol=1e-10) != 0:
            ops.wipe()
            raise RuntimeError("지진 선형성 검증용 직접해석(환경 단계)이 수렴하지 않았습니다.")
    F1, R1 = element_local_forces(), ground_reactions(model)
    apply_tag_forces({t: scale * f for t, f in load.items()}, pattern_tag=6)
    if run_static_robust(verbose=False, n_steps=1, tol=1e-10) != 0:
        ops.wipe()
        raise RuntimeError("지진 선형성 검증용 직접해석(지진 단계)이 수렴하지 않았습니다.")
    F2, R2 = element_local_forces(), ground_reactions(model)
    ops.wipe()

    dF = (F2 - F1) / scale
    ref = pred["F"]
    err_all = float(np.abs(dF - ref).max() / max(np.abs(ref).max(), 1e-30))
    err_gov = float(np.abs(dF[e] - ref[e]).max() / max(np.abs(ref[e]).max(), 1e-30))
    flips = int(np.sum((R2[:, 2] > 0.0) != (R1[:, 2] > 0.0)))
    name = f"주{'XYZ'[case['quake_gens']['primary'][p] - 1]} 꼭짓점 (부재 {gov_ele} DCR 최대)"
    out = {"scale": scale, "combo": name, "env_index": k, "err_all": err_all, "err_gov": err_gov, "flips": flips}
    print(f"[verify] 지진 선형성 (비선형 지반 + PDelta 직접해석, {name}, "
          f"{format_lmst(env['lmst_h'][k], env['meta']['sol0'])}, 지진 x{scale:,.0f}로 시험): "
          f"부재력 증분 최대 상대오차 {err_all:.2e} (지배 부재 {member_label(gov_ele)}: {err_gov:.2e}), 접촉 상태 변화 {flips}곳")
    return out


def verify_with_time_history(case, seismic, env, directions=VERIFY_TH_DIRECTIONS):
    """최종 설계 1회 검증: 실측 S1222a 원본(20 Hz, 배율 없음)을 비선형 지반 + PDelta 모델에 시간이력으로 넣고,
    부재별 '지진만에 의한 DCR 증가분'을 같은 성분만의 등가정적 상한(그 성분의 Sa_max, 모드 절대합 + 잔여질량)과 비교한다.
      비 = 시간이력 / 등가정적 <= 1 이면 등가정적이 실제 동적 응답보다 크다(보수적).
    방향마다 동적해석 1회 (지금까지 탐색 중 반복하던 계산을 마지막에 한 번만 한다)."""
    if case["quake_gens"] is None:
        return None
    design = case["design"]
    A, E, G, J, Iy, Iz = case["dome"]["section"]
    c = case["dome"]["c"]
    results = {}
    for d in directions:
        t, acc = seismic["records"][d]
        dt = seismic["dt"]
        ch = DIR_TO_CHANNEL[d]
        print(f"\n[TH] 시간이력 검증 {ch} ({'XYZ'[d - 1]}방향): {t[-1]:.0f} s, {len(acc):,}스텝, 원본 그대로")
        st = build_lc1_state(case["section"], design, env, verbose=False)
        if not st["ok"]:
            ops.wipe()
            print(f"[TH] LC1 정적해석 실패 ({st['reason']}) -> 이 방향 검증 생략")
            continue
        model = st["model"]
        out_dir = OUTPUT_DIR / "th_verify" / design.label() / f"dir_{d}"
        force_paths = setup_member_force_recorders(model["elements"], out_dir=str(out_dir))
        env_paths = setup_foundation_envelope_recorders(model, str(out_dir))
        setup_rayleigh_damping(xi=SEISMIC_DAMPING, mode_i=1, mode_j=3)
        apply_earthquake(acc, dt=dt, direction=d)
        ok = run_dynamic_robust(num_steps=len(acc) - 1, dt=dt)
        ops.wipe()
        if not ok:
            print("[TH] 동적해석 일부 스텝 불수렴 -> 결과는 참고용")

        own = build_quake_generators(case["unit"], seismic, seismic["Sa_own"], directions=(d,))
        ratios, dN_th, dN_eqs = {}, 0.0, 0.0
        for e_idx, ele in enumerate(ELE_TAGS):
            a = case["assessment"][ele]
            planes = _dcr_planes(MEMBER_FY * A, a["compression_capacity"], MEMBER_FY * Iy / c, MEMBER_FY * Iz / c)
            dF = _read_table(force_paths[ele])[:, 1:13] - st["F_base"][e_idx, :12]   # 시간이력 기록은 I·J단만
            th_inc = max(_quake_dcr_increment(_end_forces(dF, end), planes) for end in ("I", "J"))
            eqs_inc, dNb = 0.0, 0.0
            for end in ("I", "J"):
                V, dn = _zonotope_vertices(own["F"][:, :, e_idx, :], end, planes)
                eqs_inc, dNb = max(eqs_inc, _quake_dcr_increment(V, planes)), max(dNb, dn)
            ratios[ele] = th_inc / eqs_inc if eqs_inc > 0 else np.inf
            dN_th = max(dN_th, float(max(np.abs(_end_forces(dF, end)[:, 0]).max() for end in ("I", "J"))))
            dN_eqs = max(dN_eqs, dNb)

        Renv = _read_envelope(env_paths["R"]).reshape(3, -1, 3)
        dRz_th = float(np.max(np.abs(np.stack([Renv[0, :, 2], Renv[1, :, 2]]) - st["R_base"][None, :, 2])))
        dRz_eqs = float(quake_bound_abs(own["R"][..., 2]).max())
        contact0 = st["R_base"][:, 2] > 0.0
        flips = int(np.sum(((Renv[0, :, 2] > 0.0) != contact0) | ((Renv[1, :, 2] > 0.0) != contact0)))
        r_ele = max(ratios, key=ratios.get)
        r_max = ratios[r_ele]
        results[d] = {"channel": ch, "ratio_max": r_max, "ratio_ele": r_ele, "ratios": ratios,
                      "ratio_median": float(np.median(list(ratios.values()))),
                      "dN_th": dN_th, "dN_eqs": dN_eqs, "dRz_th": dRz_th, "dRz_eqs": dRz_eqs,
                      "flips": flips, "dynamic_ok": ok, "force_paths": force_paths, "F_base": st["F_base"]}
        lam = case["quake_margin"]["overall"]
        verdict = ("등가정적이 시간이력보다 큼 (보수적)" if r_max <= 1.0 else
                   f"일부 부재에서 시간이력이 큼 -> 지진을 {r_max:.2f}배 해도 λ={lam:.3g} {'> 그 배율, 여전히 SAFE' if lam > r_max else '<= 그 배율, 설계 재검토 필요'}")
        print(f"[TH] {ch}: 지진 DCR 증가분 비(시간이력/등가정적) 최대 {r_max:.3f} (부재 {r_ele}), "
              f"중앙값 {results[d]['ratio_median']:.3f} | 축력 증분 최대 {dN_th:,.3g} / {dN_eqs:,.3g} N | "
              f"지반 수직반력 증분 최대 {dRz_th:,.3g} / {dRz_eqs:,.3g} N | 접촉 상태 변화 {flips}곳")
        print(f"[TH] 판정: {verdict}")
    return results


def fatigue_upper_bound(case, seismic):
    """장시간 진동(S1222a 코다 4~10시간) 피로 상한 — 일부러 비현실적으로 불리하게 잡는다.
      응력범위 Δσ = 2 x (지진 상한의 최대 연단응력, 생성자 절대합) : 매 주기 최대 진폭으로 양진
      주기 수   N = (기록이 담을 수 있는 최고 진동수 = 나이퀴스트) x 1솔 x 설계수명 : 매일 하루 종일 최대 진동
    EN 1993-1-9: 가장 불리한 상세범주 36, γMf = 1.35.
      (1) γFf·Δσ <= ΔσL/γMf (절단한계) 이면 피로 손상 0
      (2) 절단한계를 무시하고 기울기 m=5 S-N 곡선을 끝까지 연장한 Miner 손상도 D도 함께 계산 (더 보수적)"""
    if case["quake_gens"] is None:
        return None
    A, E, G, J, Iy, Iz = case["dome"]["section"]
    c = case["dome"]["c"]
    gen = case["quake_gens"]["F"]                         # (주방향, 생성자, 부재, 12)
    sig = 0.0
    for end in MEMBER_SECTIONS:
        NMM = _end_forces(gen, end)                       # (주방향, 생성자, 부재, 3)
        sg = np.abs(NMM[..., 0]) / A + np.abs(NMM[..., 1]) * c / Iy + np.abs(NMM[..., 2]) * c / Iz
        sig = max(sig, float(sg.sum(axis=1).max()))
    dsig = 2.0 * sig * FATIGUE_GAMMA_FF
    f_cyc = seismic["nyquist_hz"]
    n_sol = f_cyc * SOL_SECONDS
    n_life = n_sol * MARS_YEAR_SOLS * FATIGUE_DESIGN_LIFE_MARS_YEARS
    dsC = FATIGUE_DETAIL_CATEGORY_MPA * 1e6
    dsD = (2.0 / 5.0) ** (1.0 / 3.0) * dsC          # 5백만 회 일정진폭 피로한계
    dsL = (5.0 / 100.0) ** (1.0 / 5.0) * dsD        # 1억 회 절단한계
    lim = dsL / FATIGUE_GAMMA_MF
    N_R = 5.0e6 * ((dsD / FATIGUE_GAMMA_MF) / max(dsig, 1e-30)) ** 5
    return {"sigma_amp_Pa": sig, "dsig_Pa": dsig, "f_cycle_hz": f_cyc, "n_per_sol": n_sol, "n_life": n_life,
            "dsC_Pa": dsC, "dsD_Pa": dsD, "dsL_Pa": dsL, "limit_Pa": lim, "below_cutoff": dsig <= lim,
            "cutoff_margin": lim / max(dsig, 1e-30), "miner_D_no_cutoff": n_life / N_R}


def print_seismic_report(case, seismic, fat=None):
    """지진 결과 요약. 선형성 오차([verify])와 시간이력 비([TH])는 각 검증 함수가 이미 출력하므로 여기서 다시 출력하지 않는다."""
    chk = case["seismic_checks"]
    if not chk:
        print("[seismic] 지진 제외 해석")
        return
    mg = case["quake_margin"]
    a = case["assessment"]
    print("\n========== 지진 (등가정적) ==========")
    sa = chk["Sa"]
    print(f"[seismic] 설계 가속도 X {sa[1]:.3e} / Y {sa[2]:.3e} / Z {sa[3]:.3e} m/s^2 "
          f"({'주기 무관 상한' if SEISMIC_USE_PERIOD_ENVELOPE else 'Sa(T1)'}, 감쇠 {SEISMIC_DAMPING*100:.0f}%"
          + (", 수평 두 성분 중 큰 값을 X·Y에 모두 적용" if SEISMIC_HORIZONTAL_ENVELOPE else "")
          + f"), = 화성 중력의 {max(sa.values())/planet_g():.2e}배, 타임라인 {len(case['env_e']['t_s'])}시각 전부에 더함")
    mr = chk["mass_ratio"]
    print(f"[seismic] 모드 {chk['n_modes']}개 (유효질량 X {mr[1]:.1%}, Y {mr[2]:.1%}, Z {mr[3]:.1%}) + 잔여질량 3방향, "
          f"모드 절대합(모든 모드 동시 최대·같은 부호) x 주방향 3가지(100/30/30)")
    T1, sat = chk["T1"], chk["Sa_T1"]
    print(f"[seismic] 돔 고유진동수(LC1 선형) {[round(f, 3) for f in chk['freqs']]} Hz -> T1 = {T1:.3f} s, "
          f"Sa(T1) X {sat[1]:.3e} / Y {sat[2]:.3e} / Z {sat[3]:.3e} -> 설계값이 Sa(T1)의 "
          f"{sa[1]/max(sat[1], sat[2], 1e-30):.2f}배 (수평)")
    if not chk["f1_below_nyquist"]:
        print(f"[warn] 1차 고유진동수가 기록 나이퀴스트({seismic['nyquist_hz']:.0f} Hz)보다 높습니다. "
              f"이 경우 Sa는 PGA에 수렴하므로 주기 무관 상한(PGA 포함)을 쓰는 지금 설정이 맞습니다.")
    print(f"[seismic] 관성 질량 {chk['inertia_mass_kg']/1000:,.0f} t -> 질량 x Sa: 수평 {chk['V_h_N']:,.4g} N, "
          f"수직 {chk['V_v_N']:,.4g} N")
    print(f"[seismic] 비교: 내부기압 들림 {case['F_up_N']/1e6:,.1f} MN, 전체 무게 {case['W_total_N']/1e6:,.1f} MN "
          f"-> 지진 수평력 = 들림의 {chk['ratio_Vh_Fup']:.2e}배, 무게의 {chk['ratio_Vh_W']:.2e}배")
    print(f"[seismic] 지진에 의한 내부 체적 변화 최대 {chk['dV_max_m3']:.3e} m^3 -> 내부기압 변화 "
          f"{chk['dp_internal_max_Pa']:.3e} Pa (내부 {INTERNAL_PRESSURE_PA:,.0f} Pa, 단열 γ={AIR_GAMMA})")
    print(f"[seismic] 지진 최대 수평변위 상한(지반 대비, 기초 이동 포함) {chk['h_disp_max_m']*1000:.4e} mm "
          f"(스팬의 {chk['h_disp_ratio']:.2e})")
    gr = case["ground"]
    print(f"[seismic] 지진에 의한 접촉 상태 추가 변화 {gr.get('n_flip_quake', 0)}곳 "
          f"(0이면 지진이 LC1 기준상태를 바꾸지 않음 -> 선형 중첩 성립)"
          + (f", 접촉이 바뀌려면 지진이 최소 {gr['min_contact_margin']:.3g}배 커야 함" if "min_contact_margin" in gr else ""))
    dq = max(v["quake_dcr_increment"] for v in a.values())
    print(f"[seismic] 지진에 의한 부재 DCR 증가분 최대 {dq:.3e}, 축력 증분 최대 "
          f"{max(v['quake_max_dN'] for v in a.values()):,.4g} N")
    print(f"[margin] 지진 여유배율 λ (S1222a를 몇 배까지 키워도 견디나): 부재 {mg['members']:.3g}, "
          f"들림 {mg.get('uplift', np.inf):.3g}, 지지력 {mg.get('bearing', np.inf):.3g}, "
          f"활동 {mg.get('sliding', np.inf):.3g}, 링 압축 {mg.get('ring', np.inf):.3g} -> 전체 {mg['overall']:.3g}배")
    if fat:
        print(f"[fatigue] 상한 가정: 지진 최대 응력범위 {fat['dsig_Pa']/1e6:.3e} MPa가 {fat['f_cycle_hz']:.0f} Hz로 "
              f"매일 하루 종일 {FATIGUE_DESIGN_LIFE_MARS_YEARS:.0f}화성년 ({fat['n_life']:.2e}회)")
        print(f"[fatigue] EN 1993-1-9 곡선, 상세범주 {FATIGUE_DETAIL_CATEGORY_MPA:.1f} MPa ({MEMBER_MATERIAL}): 절단한계/γMf = "
              f"{fat['limit_Pa']/1e6:.2f} MPa -> {'절단한계 이하, 피로 손상 0' if fat['below_cutoff'] else '절단한계 초과'} "
              f"(여유 {fat['cutoff_margin']:.2e}배) | 절단한계 무시 Miner 손상도 {fat['miner_D_no_cutoff']:.2e}")


def seismic_summary_dict(case, seismic, fat=None, lin=None, th=None):
    """Word 보고서에 옮겨 적을 핵심 수치 (JSON 저장용)"""
    chk = case["seismic_checks"] or {}
    a = case["assessment"]
    out = {
        "event": MARS_EVENT, "design": case["design"].label(), "steel_ton": case["actual_steel_mass_ton"],
        "member_material": MEMBER_MATERIAL_INFO,
        "dome": {**DOME_INFO, "enclosed_volume_m3": DOME_ENCLOSED_VOLUME_M3, "shell_area_m2": DOME_SHELL_AREA_M2},
        "damping": SEISMIC_DAMPING, "period_envelope": SEISMIC_USE_PERIOD_ENVELOPE,
        "components": {c["channel"]: {"pga_m_s2": c["pga"], "Sa_max_m_s2": c["Sa_env"], "T_peak_s": c["T_peak"],
                                      "duration_h": c["info"]["duration_s"] / 3600.0}
                       for c in seismic["comps"].values()},
        "Sa_design_m_s2": {"X": chk.get("Sa", {}).get(1), "Y": chk.get("Sa", {}).get(2), "Z": chk.get("Sa", {}).get(3)},
        "freqs_hz": chk.get("freqs"), "T1_s": chk.get("T1"), "n_modes": chk.get("n_modes"),
        "modal_mass_ratio": chk.get("mass_ratio"),
        "inertia_mass_t": chk.get("inertia_mass_kg", 0.0) / 1000.0,
        "V_h_N": chk.get("V_h_N"), "F_up_MN": case["F_up_N"] / 1e6, "W_total_MN": case["W_total_N"] / 1e6,
        "ratio_Vh_Fup": chk.get("ratio_Vh_Fup"), "dp_internal_max_Pa": chk.get("dp_internal_max_Pa"),
        "n_flip_quake": case["ground"].get("n_flip_quake"),
        "contact_flip_margin": case["ground"].get("min_contact_margin"),
        "max_dcr": case["max_dcr"], "max_dcr_env_only": case["max_dcr_env_only"],
        "max_dcr_lc2": case["max_dcr_lc2"], "gov_ele": case["gov_ele"], "gov_case": a[case["gov_ele"]]["gov_case"],
        "quake_dcr_increment_max": max(v["quake_dcr_increment"] for v in a.values()),
        "quake_margin": case["quake_margin"],
    }
    if fat:
        out["fatigue"] = {k: fat[k] for k in ("dsig_Pa", "n_life", "limit_Pa", "below_cutoff", "cutoff_margin",
                                              "miner_D_no_cutoff")}
    if lin:
        out["linearity"] = {k: lin[k] for k in ("scale", "combo", "err_all", "err_gov", "flips")}
    if th:
        out["time_history"] = {r["channel"]: {"ratio_max": r["ratio_max"], "ratio_ele": r["ratio_ele"],
                                              "ratio_median": r["ratio_median"], "flips": r["flips"]}
                               for r in th.values()}
    return out


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def find_minimum_safe_mass(design, seismic, env,
                           include_quake=True,
                           start_mass=MASS_SEARCH_START_TON,
                           coarse_step=MASS_SEARCH_COARSE_STEP_TON,
                           tol=MASS_SEARCH_TOL_TON,
                           max_mass=MASS_SEARCH_MAX_TON,
                           quake_actual_idx=None,
                           known_failed_ton=None):
    """설계(D, t, B) 하나에서 고정폭 브라케팅 + 이분법으로 최소 SAFE 부재 질량을 찾는다.
      - 시작 질량이 SAFE면 아래로, FAIL이면 위로 브라케팅한다.
      - known_failed_ton: 이미 FAIL로 아는 질량 (지진 제외 탐색의 최대 FAIL). 지진은 하중만 더하므로
                          지진 포함에서도 FAIL이다 -> 시작 질량이 SAFE면 아래로 다시 찾지 않고 이 값을 FAIL 경계로 쓴다.
      - 부재는 통과했는데 기초가 탈락하면 FoundationInfeasible (부재 질량을 키워도 해결 안 됨).
    반환되는 min_safe_mass_ton은 최종 구간의 SAFE측 상한이다. (max_failed < 임계질량 <= min_safe)
    질량 후보마다의 판정 줄은 SEARCH_VERBOSE = True일 때만 출력한다."""
    if start_mass <= 0.0 or coarse_step <= 0.0 or tol <= 0.0:
        raise ValueError("start_mass, coarse_step, tol은 0보다 커야 합니다.")
    if max_mass <= start_mass:
        raise ValueError("max_mass는 start_mass보다 커야 합니다.")
    for name, val in (("start_mass", start_mass), ("coarse_step", coarse_step)):
        if not math.isclose(val / tol, round(val / tol)):
            raise ValueError(f"{name}는 {tol:.1f} ton 단위여야 합니다.")

    cache = {}
    tag = "지진 포함" if include_quake else "지진 제외"

    def evaluate(mass_ton, stage):
        mass_ton = round(float(mass_ton), 9)
        if mass_ton in cache:
            return cache[mass_ton]
        case = run_case(mass_ton, design, seismic, env, include_quake=include_quake,
                        quake_actual_idx=quake_actual_idx, verbose=SEARCH_VERBOSE)
        members_ok = case["members_ok"] and case["disp_ok"]
        reason = case["reason"] or ""
        # 지지력·링 압축 탈락은 강재를 늘려도 (무게가 늘어) 나빠질 뿐이므로 바로 기초를 바꾼다.
        # 들림·활동은 강재 무게가 조금 도움이 되므로 부재가 통과한 뒤에만 기초 탈락으로 본다.
        # (탈락 이유는 부르는 쪽(solve_cell 등)이 출력한다)
        if not case["foundation_ok"] and (members_ok or "지지력" in reason or "링 압축" in reason or "침하" in reason):
            raise FoundationInfeasible(design, reason or "기초 탈락", mass_ton)
        safe = bool(case["passed"])
        if SEARCH_VERBOSE:
            verdict = "SAFE" if safe else (f"FAIL ({case['reason']})" if case["reason"]
                                            else f"FAIL (파손 부재 {case['failed_count']}개)")
            print(f"[{stage} | {design.label()} | {tag}] 질량 {mass_ton:.1f} ton -> {verdict} "
                  f"(최대 DCR {case['max_dcr']:.4f})")
        if not safe:   # 탈락 케이스의 큰 배열은 바로 버려 메모리 절약
            case = {k: case[k] for k in ("max_dcr", "failed_members", "passed", "status", "reason")}
        cache[mass_ton] = (safe, case)
        return cache[mass_ton]

    if known_failed_ton is not None and known_failed_ton >= start_mass:
        known_failed_ton = None
    if evaluate(start_mass, "브라케팅")[0]:
        upper_safe, lower_failed = float(start_mass), 0.0
        if known_failed_ton is not None:
            lower_failed = float(known_failed_ton)
        else:
            m = start_mass - coarse_step
            while m >= tol - 1e-9:
                if evaluate(m, "브라케팅(하향)")[0]:
                    upper_safe = m
                    m -= coarse_step
                else:
                    lower_failed = m
                    break
    else:
        lower_failed = float(start_mass)
        while True:
            candidate = lower_failed + coarse_step
            if candidate > max_mass:
                raise RuntimeError(f"탐색 상한 {max_mass:.1f} ton까지 SAFE 질량을 찾지 못했습니다. ({design.label()})")
            if evaluate(candidate, "브라케팅")[0]:
                upper_safe = candidate
                break
            lower_failed = candidate

    low_index, high_index = int(round(lower_failed / tol)), int(round(upper_safe / tol))
    iteration = 0
    while high_index - low_index > 1:
        iteration += 1
        middle_index = (low_index + high_index) // 2
        if evaluate(middle_index * tol, f"이분법 {iteration}")[0]:
            high_index = middle_index
        else:
            low_index = middle_index

    lower_failed, upper_safe = low_index * tol, high_index * tol
    final_safe, final_case = evaluate(upper_safe, "최종 확인")
    check_lower = lower_failed > 0 and not (known_failed_ton is not None and math.isclose(lower_failed, known_failed_ton))
    lower_safe = evaluate(lower_failed, "최종 확인")[0] if check_lower else False
    if not final_safe or lower_safe:
        raise RuntimeError("최종 FAIL-SAFE 경계 검증이 일관되지 않습니다.")

    if SEARCH_VERBOSE:
        print(f"[search] 최소 SAFE 질량 ({design.label()}, {tag}): {upper_safe:.1f} ton "
              f"(최대 FAIL {lower_failed:.1f} ton, 이분법 {iteration}회)")
    return {"min_safe_mass_ton": float(upper_safe), "max_failed_mass_ton": float(lower_failed),
            "final_case": final_case, "design": design, "include_quake": include_quake}


def solve_cell(D, B, seismic, env, include_quake=SEARCH_WITH_QUAKE, quake_actual_idx=None,
               start_mass=MASS_SEARCH_START_TON, coarse_step=MASS_SEARCH_COARSE_STEP_TON):
    """칸(D, B) 하나를 끝까지 푼다 (안쪽 반복).
      1바퀴: 돔 0 t로 기초 설계 -> 그 기초 위에서 돔 부재 최소 질량 M
      2바퀴: M으로 기초 다시 설계 -> 덮개·바닥판이 그대로면 끝, 바뀌면 M 근처에서 다시 찾는다 (좁은 브라케팅)"""
    M, s, prev, notes, t_floor = 0.0, None, None, [], None
    sig_hoop = 0.0                                       # 이전 바퀴 해석의 링 축인장응력 (바닥판 휨과 합산)
    start, step = start_mass, coarse_step
    for it in range(1, CELL_ITER_MAX + 1):
        cell = design_cell(D, B, env, dome_mass_ton=M, t_floor=t_floor, sig_hoop=sig_hoop)
        if not cell["ok"]:
            reason = cell["reason"] if s is None else f"{it}바퀴 기초 재설계 불성립: {cell['reason']}"
            return {"ok": False, "reason": " | ".join(notes + [reason])}
        design = cell["design"]
        if s is not None and design == prev:
            break                                        # 실제 질량으로 다시 설계해도 기초가 같음 -> 수렴
        if it == CELL_ITER_MAX:
            continue                                     # 마지막 바퀴인데 아직 바뀜 -> 수렴 실패로 처리
        print(f"[cell D{D:.1f} B{B:.1f}] {it}바퀴: 돔 {M:.0f} t로 기초 설계 -> 덮개 {design.cover_t:.1f} m, "
              f"바닥판 {design.flange_t:.1f} m" + (f" + 철근 {design.flange_rebar_kg/1000:.1f} t" if design.flange_rebar_kg else "")
              + f" (들림 저항/들림 {cell['uplift_ratio']:.2f}, 지지 {cell['q_lc2']/1e3:.0f}/{cell['q_allow']/1e3:.0f} kPa"
              + (f", 링 인장 {sig_hoop/1e6:.2f} MPa 합산" if sig_hoop else "") + ")")
        prev = design
        try:
            s = find_minimum_safe_mass(design, seismic, env, include_quake=include_quake, start_mass=start,
                                       coarse_step=step, quake_actual_idx=quake_actual_idx)
        except FoundationInfeasible as exc:
            if "지지력" in exc.reason or "링 압축" in exc.reason or "침하" in exc.reason:
                return {"ok": False, "reason": " | ".join(notes + [f"해석 단계 기초 탈락 ({exc.reason})"])}
            t_floor = design.cover_t + FOUNDATION_RETRY_STEP_M         # 들림·활동: 덮개를 올려 다시
            notes.append(f"{exc.reason} -> 덮개 {t_floor:.1f} m 이상으로")
            print(f"[cell D{D:.1f} B{B:.1f}] 해석 단계 기초 탈락 ({exc.reason}) -> 덮개를 {t_floor:.1f} m 이상으로")
            start, step = exc.mass_ton or start, MASS_SEARCH_GUESS_STEP_TON
            s, prev = None, None
            continue
        M = s["min_safe_mass_ton"]
        sig_hoop = ring_hoop_stress(design, s["final_case"])
        start = min(max(M, MASS_SEARCH_TOL_TON), MASS_SEARCH_MAX_TON - MASS_SEARCH_COARSE_STEP_TON)
        step = MASS_SEARCH_GUESS_STEP_TON
    else:
        return {"ok": False, "reason": " | ".join(notes + [f"기초-돔 질량 반복이 {CELL_ITER_MAX}바퀴 안에 수렴하지 않음 (제외)"])}
    return {"ok": True, "search": s, "cell": cell, "iterations": it, "notes": notes}


def _grid_row(D, B, res):
    s, cell = res["search"], res["cell"]
    fc, design = s["final_case"], s["design"]
    ring, ground = fc["ring"], fc["ground"]
    gov = fc["assessment"][fc["gov_ele"]]
    Lr = RING_LAYOUT["arc_length"]
    hoop_t, flange_t_rebar = ring["rebar_kg"] / 1000.0, design.flange_rebar_kg / 1000.0
    cov = cover_volume_m3(design.cover_t)
    exc = design.ring_b * design.ring_depth * Lr
    conc = ring_section(design.ring_depth, design.ring_b, tf=design.flange_t)["A"] * Lr
    st_ = ground.get("settle", {})
    return {
        "D_m": D, "B_m": B, "feasible": True, "t_m": design.cover_t, "flange_t_m": design.flange_t,
        "t_governs": cell["t_governs"], "flange_type": "철근" if design.flange_rebar_kg > 0 else "무근",
        "steel_ton": s["min_safe_mass_ton"], "max_failed_ton": s["max_failed_mass_ton"],
        "ring_rebar_ton": hoop_t, "flange_rebar_ton": flange_t_rebar, "rebar_ton": hoop_t + flange_t_rebar,
        "material_ton": s["min_safe_mass_ton"] + hoop_t + flange_t_rebar,
        "cover_volume_m3": cov, "excavation_m3": exc, "ring_concrete_m3": conc,
        "dome_kg_per_m3": s["min_safe_mass_ton"] * 1000.0 / DOME_ENCLOSED_VOLUME_M3,
        "cover_ton": fc["mass_parts_kg"]["cover"] / 1000.0, "ring_concrete_ton": fc["mass_parts_kg"]["ring"] / 1000.0,
        "envelope_ton": fc["mass_parts_kg"]["membrane"] / 1000.0,
        "envelope_lc2_dcr": cell.get("envelope_lc2_dcr", 0.0),
        "max_dcr": fc["max_dcr"], "gov_case": gov["gov_case"], "gov_ele": fc["gov_ele"],
        "uplift_ratio_static": cell["uplift_ratio"], "bearing_q_kPa": cell["q_lc2"] / 1e3,
        "bearing_q_allow_kPa": cell["q_allow"] / 1e3, "flange_M_kNm_per_m": cell["flange"]["M"] / 1e3,
        "uplift_dcr": ground["uplift_dcr"], "bearing_dcr": ground["bearing_dcr"], "sliding_dcr": ground["sliding_dcr"],
        "n_flip": ground["n_flip"], "n_lifted": ground["n_lifted"], "ring_comp_dcr": ring["comp_dcr"],
        "settle_max_mm": st_.get("max_m", 0.0) * 1000.0, "settle_diff_mm": st_.get("diff_m", 0.0) * 1000.0,
        "angular_distortion": st_.get("angular", 0.0),
        "cell_iterations": res["iterations"], "quake_included": s["include_quake"], "note": " | ".join(res["notes"]),
    }


def search_ring_cover(seismic, env, quake_actual_idx=None, depth_grid=RING_DEPTH_GRID, b_grid=RING_B_GRID,
                      include_quake=SEARCH_WITH_QUAKE):
    """링 깊이 D x 바닥판 폭 B 모든 칸을 끝까지 풀고(바깥 반복), 선택 기준으로 최종안을 고른다."""
    rows, searches = [], {}
    print("\n" + "=" * 70)
    print(f"[search] 링 깊이 D {list(depth_grid)} m x 바닥판 폭 B {b_grid[0]:.1f}~{b_grid[-1]:.1f} m "
          f"({len(b_grid)}개) = {len(depth_grid) * len(b_grid)}칸, "
          f"{'지진 포함' if include_quake else '지진 제외 (최종안만 2단계에서 지진 확인)'}")
    print(f"[search] 링 재료 '{RING_MATERIAL}', 바닥판 {'무근' if FLANGE_MODE == 'plain' else '무근/철근 자동'}, FS 들림 {FS_UPLIFT}/지지력 {FS_BEARING}/"
          f"활동 {FS_SLIDING}/바닥판 {FS_FLANGE}, 방사선 하한 {COVER_MIN_THICKNESS_M} m, 레골리스 {SITE_REGOLITH_THICKNESS_M} m")
    print("=" * 70)
    for D in depth_grid:
        if D > SITE_REGOLITH_THICKNESS_M + 1e-9 or D < RING_FLANGE_T + STEM_H_MIN - 1e-9:
            rows += [{"D_m": D, "B_m": B, "feasible": False, "note": "링 깊이 범위 밖"} for B in b_grid]
            continue
        start, step = MASS_SEARCH_START_TON, MASS_SEARCH_COARSE_STEP_TON
        for B in b_grid:
            res = solve_cell(D, B, seismic, env, include_quake, quake_actual_idx, start_mass=start, coarse_step=step)
            if not res["ok"]:
                rows.append({"D_m": D, "B_m": B, "feasible": False, "note": res["reason"]})
                print(f"[cell D{D:.1f} B{B:.1f}] 탈락: {res['reason']}")
                continue
            row = _grid_row(D, B, res)
            rows.append(row)
            searches[res["search"]["design"]] = res["search"]
            start = min(max(round(row["steel_ton"] / MASS_SEARCH_TOL_TON) * MASS_SEARCH_TOL_TON, MASS_SEARCH_TOL_TON),
                        MASS_SEARCH_MAX_TON - MASS_SEARCH_COARSE_STEP_TON)   # 옆 칸은 비슷한 질량에서 시작 (시간 절약)
            step = MASS_SEARCH_GUESS_STEP_TON                                 # 가까운 값에서 시작하므로 좁게 브라케팅
            print(f"[cell D{D:.1f} B{B:.1f}] 완료: 덮개 {row['t_m']:.1f} m, 바닥판 {row['flange_t_m']:.1f} m "
                  f"({row['flange_type']}), 돔 부재 {row['steel_ton']:.0f} t ({row['dome_kg_per_m3']:.2f} kg/m^3), 링 콘크리트 {row['ring_concrete_ton']:,.0f} t, "
                  f"침하 {row['settle_max_mm']:.1f} mm ({row['cell_iterations']}바퀴)")

    df = pd.DataFrame(rows)
    df["member_material"] = MEMBER_MATERIAL
    df["dome_type"] = DOME_TYPE
    ok = df[df["feasible"]].copy()
    if ok.empty:
        raise RuntimeError("성립하는 (D, B) 칸이 없습니다. B 격자, 덮개 상한, 바닥판 방식(FLANGE_MODE)을 확인하세요.")

    # ---- 파레토 해: 돔 재료(돔 부재 + 철근)와 기초 재료(링 콘크리트)가 둘 다 더 작은 다른 칸이 없는 칸 ----
    X = ok[["material_ton", "ring_concrete_ton"]].to_numpy()
    is_pareto = [not np.any(np.all(X <= x, axis=1) & np.any(X < x, axis=1)) for x in X]
    df["pareto"] = False
    df.loc[ok.index[is_pareto], "pareto"] = True
    end_rows = {"material_min": ok.sort_values(["material_ton", "ring_concrete_ton"]).iloc[0],
                "concrete_min": ok.sort_values(["ring_concrete_ton", "material_ton"]).iloc[0]}
    df["end"] = ""
    for name, r in end_rows.items():
        df.loc[r.name, "end"] = (df.loc[r.name, "end"] + "+" if df.loc[r.name, "end"] else "") + name
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / "ring_cover_search.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")

    cols = ["D_m", "B_m", "t_m", "flange_t_m", "flange_type", "steel_ton", "rebar_ton", "material_ton",
            "ring_concrete_ton", "cover_ton", "settle_max_mm", "gov_case", "cell_iterations"]
    print("\n" + ok[cols].to_string(index=False, float_format=lambda v: f"{v:,.2f}"))
    print(f"\n[search] 결과 표 저장: {csv_path.resolve()}")
    par = df[df["pareto"]].sort_values("ring_concrete_ton")
    print(f"[pareto] 돔 재료 vs 기초 재료 파레토 해 {len(par)}개 / 성립 {len(ok)}개:")
    for _, r in par.iterrows():
        print(f"         D={r['D_m']:.1f} B={r['B_m']:.1f} m: 돔 재료 {r['material_ton']:.1f} t "
              f"({r['dome_kg_per_m3']:.2f} kg/m^3), 링 콘크리트 {r['ring_concrete_ton']:,.0f} t, 덮개 {r['t_m']:.1f} m"
              + (f"  <- {', '.join(END_LABELS[e] + ('' if e in FINAL_ENDS else ' (참고)') for e in r['end'].split('+'))}"
                 if r["end"] else ""))
    for name, r in end_rows.items():
        if name not in FINAL_ENDS:
            continue
        edge = [f"B가 탐색 상한 {max(b_grid):.1f} m" for _ in [0] if math.isclose(r["B_m"], max(b_grid))] + \
               [f"D가 탐색 상한 {max(depth_grid):.1f} m" for _ in [0] if math.isclose(r["D_m"], max(depth_grid))]
        if edge:
            print(f"[pareto] 주의: {END_LABELS[name]}이 탐색 범위 끝에 있음 ({', '.join(edge)}) -> 상한의 근거를 보고서에 명시")
    ends = {}
    for name, r in end_rows.items():
        d = next(x for x in searches if math.isclose(x.ring_depth, r["D_m"]) and math.isclose(x.ring_b, r["B_m"]))
        ends[name] = (d, searches[d])
    plot_search_results(df, OUTPUT_DIR / "ring_cover_search.png")
    return df, ends


def confirm_with_quake(design, search_no_quake, seismic, env, quake_actual_idx=None):
    """2단계: 지진 제외로 찾은 최종안을 지진 포함으로 다시 확인한다.
    지진은 하중만 더하므로 지진 제외에서 FAIL인 질량은 지진 포함에서도 FAIL이다.
    -> 지진 제외 최소 질량 M 하나만 지진 포함으로 풀어 SAFE면 그대로 끝 (해석 1회). FAIL이면 M부터 위로 찾는다.
    지진 중 기초가 탈락하면 같은 (D, B)에서 덮개를 올려 다시 설계한다 (설계가 바뀌므로 그때는 처음부터 탐색)."""
    tol = MASS_SEARCH_TOL_TON
    start = max(search_no_quake["min_safe_mass_ton"], tol)
    known = search_no_quake["max_failed_mass_ton"]
    for _ in range(FOUNDATION_RETRY_MAX + 1):
        try:
            s = find_minimum_safe_mass(design, seismic, env, include_quake=True, start_mass=start,
                                       coarse_step=MASS_SEARCH_GUESS_STEP_TON, quake_actual_idx=quake_actual_idx,
                                       known_failed_ton=known)
            add = s["min_safe_mass_ton"] - search_no_quake["min_safe_mass_ton"]
            print(f"[confirm] 지진 포함 최소 질량 {s['min_safe_mass_ton']:.1f} ton "
                  f"(지진 제외 대비 {add:+.1f} ton)")
            return s
        except FoundationInfeasible as exc:
            start, known = exc.mass_ton or start, None
            cell = design_cell(design.ring_depth, design.ring_b, env, dome_mass_ton=start,
                               t_floor=design.cover_t + FOUNDATION_RETRY_STEP_M,
                               sig_hoop=ring_hoop_stress(design, search_no_quake["final_case"]))
            print(f"[confirm] 지진 포함 시 기초 탈락 ({exc.reason}) -> "
                  + (f"덮개 {cell['design'].cover_t:.1f} m로 재설계" if cell["ok"] else cell["reason"]))
            if not cell["ok"]:
                raise RuntimeError(f"기초 재설계 불성립: {cell['reason']}") from exc
            design = cell["design"]
    raise RuntimeError("지진 포함 확인에서 재시도 한도를 넘었습니다.")


def solve_bare(seismic, env, quake_actual_idx=None):
    """기초 없음 (FOUNDATION_MODE = "none"): 덮개 <-> 돔 질량 반복 (solve_cell과 같은 방식) 후 지진 포함 최종 탐색."""
    M, t_floor, prev, s = 0.0, None, None, None
    start, step = MASS_SEARCH_START_TON, MASS_SEARCH_COARSE_STEP_TON
    for _ in range(CELL_ITER_MAX):
        design = Design(cover_t=bare_cover(env, M, t_floor))
        if s is not None and design == prev:
            break                                        # 실제 돔 질량으로 다시 구해도 덮개가 같음 -> 수렴
        prev = design
        try:
            s = find_minimum_safe_mass(design, seismic, env, include_quake=SEARCH_WITH_QUAKE, start_mass=start,
                                       coarse_step=step, quake_actual_idx=quake_actual_idx)
        except FoundationInfeasible as exc:
            t_floor = design.cover_t + FOUNDATION_RETRY_STEP_M
            print(f"[none] 지반 탈락 ({exc.reason}) -> 덮개를 {t_floor:.1f} m 이상으로")
            start, step, s, prev = exc.mass_ton or start, MASS_SEARCH_GUESS_STEP_TON, None, None
            continue
        M = s["min_safe_mass_ton"]
        start = min(max(M, MASS_SEARCH_TOL_TON), MASS_SEARCH_MAX_TON - MASS_SEARCH_COARSE_STEP_TON)
        step = MASS_SEARCH_GUESS_STEP_TON
    else:
        raise RuntimeError(f"기초 없음: 덮개-돔 질량 반복이 {CELL_ITER_MAX}바퀴 안에 수렴하지 않음")
    print(f"[none] 수렴: 덮개 {design.cover_t:.1f} m, 돔 부재 {M:.0f} t"
          + ("" if SEARCH_WITH_QUAKE else " (지진 제외) -> 지진 포함 확인"))
    if SEARCH_WITH_QUAKE:
        return s
    # 지진은 하중만 더하므로 지진 제외 최소 질량 하나만 지진 포함으로 확인 (confirm_with_quake와 같은 방식)
    return find_minimum_safe_mass(design, seismic, env, include_quake=True, start_mass=max(M, MASS_SEARCH_TOL_TON),
                                  coarse_step=MASS_SEARCH_GUESS_STEP_TON, quake_actual_idx=quake_actual_idx,
                                  known_failed_ton=s["max_failed_mass_ton"])


# ----------------------------------------------------------------------
# 12. 보고
# ----------------------------------------------------------------------
def print_final_report(case, env, quake_info=None, quake_lmst_h=None):
    s0 = env["meta"]["sol0"]
    a = case["assessment"]
    gov = a[case["gov_ele"]]
    si = case["section_info"]
    d = case["design"]
    mp = case["mass_parts_kg"]
    ring, ground = case["ring"], case["ground"]
    env_e = case["env_e"]
    print("\n========== 최종 연구 결과 ==========")
    if quake_info is not None:
        print(f"지진: {MARS_EVENT} ({quake_info['id']}, 최대가속도 {quake_info['pga']:.3e} m/s^2"
              + (f", 지방시 {format_lmst(quake_lmst_h, s0)}" if quake_lmst_h is not None else "") + ")")
    print(f"환경: TWINS + PS 실측 {format_lmst(env['lmst_h'][0], s0)} ~ {format_lmst(env['lmst_h'][-1], s0)} "
          f"(덮개 아래 ΔT 배율 exp(-t/δ) = {env_e['thermal_factor']:.2e}"
          + (" -> 덮개가 일교차를 사실상 차단, 부재 열응력 0" if env_e["thermal_factor"] < 1e-6 else "")
          + f", 바람 {'반영' if env_e['wind_active'] else '차폐'})")
    if FOUNDATION_MODE == "ring":
        print(f"[design] 링 깊이 D={d.ring_depth:.1f} m, 바닥판 폭 B={d.ring_b:.1f} m, 바닥판 두께 {d.flange_t:.1f} m"
              + (f" (+ 철근 {d.flange_rebar_kg/1000:.1f} t)" if d.flange_rebar_kg else " (무근)")
              + f", 덮개 t={d.cover_t:.1f} m (링 재료 '{RING_MATERIAL}', 원지반 φ={SITE_PHI_DEG:.0f}°, ρ={SITE_REGOLITH_DENSITY:.0f})")
    else:
        print(f"[design] 기초 없음: 돔 밑단 절점이 원지반 위에 놓임, 덮개 t={d.cover_t:.1f} m "
              f"(들림은 덮개+외피+돔 무게로만 버팀)")
    print(f"[section] D={si['D']:.4f} m, t={si['t']:.4f} m, A={si['A']:.6f} m^2 | 재료 {MEMBER_MATERIAL_LABEL}")
    print(f"[mass] 재료 질량: 돔 부재({MEMBER_MATERIAL}) {mp['steel']/1000:,.1f} t + 링 철근 {ring['rebar_kg']/1000:,.1f} t "
          f"+ 바닥판 철근 {d.flange_rebar_kg/1000:,.1f} t = {(mp['steel'] + ring['rebar_kg'] + d.flange_rebar_kg)/1000:,.1f} t | 기초 공사: 링 황 콘크리트 {mp['ring']/1000:,.0f} t, "
          f"덮개 {mp['cover']/1000:,.0f} t ({cover_volume_m3(d.cover_t):,.0f} m^3) | 외피 {mp['membrane']/1000:.2f} t")
    print(f"[equilibrium] 무게 {case['W_total_N']/1e6:,.2f} MN, 들림(평균 순압) {case['F_up_N']/1e6:,.2f} MN, "
          f"지반 반력 합 LC2 {ground['sum_Rz_lc2']/1e6:,.2f} / LC1 {ground['sum_Rz_lc1']/1e6:,.2f} MN "
          f"(평형 오차 {case['equilibrium_error']:.1e})")
    where = (f" at {format_lmst(gov['gov_lmst_h'], s0)}" if gov["gov_case"].startswith("LC1") else "")
    print(f"[check] 지배 부재 {member_label(case['gov_ele'])} ({SECTION_LABEL[gov['gov_end']]}): DCR {gov['max_dcr']:.4f}, "
          f"{gov['gov_case']}{where}")
    print(f"        N={gov['critical_N']/1e3:,.1f} kN, My={gov['critical_My']/1e3:,.2f} kN·m, "
          f"Mz={gov['critical_Mz']/1e3:,.2f} kN·m, 세장비 {gov['slenderness']:.1f}")
    print(f"[check] 케이스별 최대 DCR: LC1 환경 단독 {case['max_dcr_env_only']:.4f}, "
          f"LC2 감압 {case['max_dcr_lc2']:.4f}, 전체 {case['max_dcr']:.4f}")
    if case["max_dcr_actual_time"] is not None:
        print(f"[check] 실제 지진 지방시의 환경 + 지진 최대 DCR {case['max_dcr_actual_time']:.4f}")
    if ground["mode"] == "none":
        print(f"[foundation] 기초 없음: 발 들림 {'있음' if ground['n_lifted'] else '없음'}, "
              f"활동 DCR {ground['sliding_dcr']:.3f} (수평 합력 {ground['H_net']/1e3:,.1f} kN), "
              f"국부 마찰 한계 도달 {ground['n_h_limit']}곳 (밑단 둘레 부재가 분담) | 발 수직반력 "
              f"{ground['min_Rz']/1e6:,.2f} ~ {ground['max_Rz']/1e6:,.2f} MN (지지력·침하는 미검토, 덮개 두께 비교용)")
    if ground["mode"] == "ring":
        print(f"[foundation] DCR 들림 {ground['uplift_dcr']:.3f}, 지지력 {ground['bearing_dcr']:.3f} "
              f"(최대 {ground['q_max']/1e3:,.0f} kPa), 활동 {ground['sliding_dcr']:.3f} "
              f"(수평 합력 {ground['H_net']/1e3:,.1f} kN) | LC1 들린 지지점 {ground['n_lifted']}곳, "
              f"하루 중 접촉상태 변화 {ground['n_flip']}곳, 국부 수평 한계 도달 {ground['n_h_limit']}곳 "
              f"(링 후프가 분담)")
        st_ = ground.get("settle", {})
        if st_:
            print(f"[settle] 감압(LC2) 침하: 최대 {st_['max_m']*1000:.1f} mm (한계 {SETTLEMENT_LIMIT_M*1000:.0f} mm), "
                  f"부등 {st_['diff_m']*1000:.1f} mm, 각변위 1/{1/max(st_['angular'], 1e-12):,.0f} "
                  f"(한계 1/{1/ANGULAR_DISTORTION_LIMIT:.0f}) | 지반반력계수 {SITE_SUBGRADE_K:.1e} N/m^3 가정")
        print(f"[ring] 최대 인장 {ring['T_max']/1e6:,.2f} MN, 최대 압축 {ring['C_max']/1e6:,.2f} MN "
              f"(압축 이용률 {ring['comp_dcr']:.3f}) -> "
              + (f"철근 {ring['A_s']*1e4:,.0f} cm^2 ({ring['rebar_kg']/1000:,.1f} t)" if ring["need_rebar"]
                 else "인장이 황 콘크리트 인장강도 이내 (철근 불필요)"))
        if not d.flange_rebar_kg:
            ft_r = RING_MATERIAL_CASES[RING_MATERIAL][1]
            fl_r = foundation_loads(d.ring_depth, d.ring_b, d.cover_t, d.flange_t, case["target_mass_ton"])
            fb_r = flange_bending(d.ring_depth, d.ring_b, d.cover_t, d.flange_t, fl_r["q_lc2"])
            sig_b, sig_h = 6.0 * fb_r["M"] / d.flange_t ** 2, ring_hoop_stress(d, case)
            dcr_r = FS_FLANGE * (sig_b + sig_h) / ft_r
            print(f"[ring] 바닥판 휨 {sig_b/1e6:.3f} + 링 인장 {sig_h/1e6:.3f} = {(sig_b + sig_h)/1e6:.3f} MPa "
                  f"(허용 f_t/FS {ft_r/FS_FLANGE/1e6:.3f} MPa, DCR {dcr_r:.3f}, 두 하중 경우 최대값끼리 합산)"
                  + (" -> [warn] 초과" if dcr_r > 1.0 else ""))
        if ground["n_flip"]:
            print("[warn] 하루 중 접촉 상태가 바뀌는 지지점이 있어 선형 중첩 오차가 커질 수 있습니다. "
                  "verify_superposition 결과를 확인하세요.")


def print_load_breakdown(case, env):
    """지배 부재의 부재력을 하중별로 분해. 덮개가 막아 0인 하중(온도·바람)은 줄을 생략한다."""
    ele = case["gov_ele"]
    gov = case["assessment"][ele]
    e, end = ELE_TAGS.index(ele), gov["gov_end"]
    if gov["gov_case"] == "LC2 감압":
        N, My, Mz = _end_forces(case["F_lc2"][e], end)
        print(f"[breakdown] 부재 {member_label(ele)} {SECTION_LABEL[end]}, LC2 감압 (자중+막+덮개+링, 내압 0)")
        print(f"  LC2           : N={N/1e3:11.1f} kN | My={My/1e3:8.2f} | Mz={Mz/1e3:8.2f} kN·m")
        N, My, Mz = _end_forces(case["F_base"][e], end)
        print(f"  (참고) LC1 기준: N={N/1e3:11.1f} kN | My={My/1e3:8.2f} | Mz={Mz/1e3:8.2f} kN·m")
        return
    k = gov["gov_env_index"]
    env_e = case["env_e"]
    u = case["unit"]["F"]
    parts = {
        "LC2(자중+덮개+링)": case["F_lc2"][e],
        "평균 순압 가압": case["F_base"][e] - case["F_lc2"][e],
    }
    if u["P"] is not None:
        parts["순압 변동"] = (env_e["dp"][k] - env_e["dp_ref"]) * u["P"][e]
    thermal_on = u["T"] is not None and env_e["thermal_factor"] >= 1e-6
    if thermal_on:
        parts["온도"] = env_e["dT"][k] * u["T"][e]
    if u["W"] is not None and env_e["wind_active"]:
        parts["바람"] = env_e["q_Pa"][k] * u["W"][env_e["sector"][k], e]
    print(f"[breakdown] 부재 {ele} {SECTION_LABEL[end]}, {format_lmst(env['lmst_h'][k], env['meta']['sol0'])}"
          + (f" (부재 ΔT {env_e['dT'][k]:+.2e} K)" if thermal_on else ""))
    for name, F in parts.items():
        N, My, Mz = _end_forces(F, end)
        print(f"  {name:<14}: N={N/1e3:11.1f} kN | My={My/1e3:8.2f} | Mz={Mz/1e3:8.2f} kN·m")


def main():
    global OUTPUT_DIR          # 양 끝 설계 결과를 하위 폴더에 나눠 저장하기 위해 바꿈
    set_planet(CURRENT_PLANET)
    print(f"[dome] {DOME_LABEL}: 절점 {len(dome_NODES)}개, 요소 {len(DOME_EDGES)}개 | 결과 폴더 {OUTPUT_DIR.resolve()}")
    print(f"[volume] 돔 내부 감싸인 체적 = {DOME_ENCLOSED_VOLUME_M3:.1f} m^3, 쉘 면적 {DOME_SHELL_AREA_M2:.1f} m^2, "
          f"순압 1 Pa당 들림 {DOME_UPLIFT_PER_PA:.1f} N")
    print(f"[material] 돔 부재 {MEMBER_MATERIAL}: E {MEMBER_E/1e9:.2f} GPa, fy {MEMBER_FY/1e6:.0f} MPa, "
          f"ρ {MEMBER_DENSITY:.1f} kg/m^3, α {MEMBER_ALPHA:.2e} 1/K")
    print(f"[material] 출처 {_MAT['file']} | status {_MAT['status']} | sources {_MAT['sources']}")
    print(f"[data] 입력 폴더: {DATA_DIR}")

    twins_paths = [DATA_DIR / f for f in TWINS_FILES]
    ps_paths = [DATA_DIR / f for f in PS_FILES]
    quake_path = DATA_DIR / QUAKE_FILE
    _require_files(*twins_paths, *ps_paths, quake_path)

    # ---- 실측 환경 타임라인 (10분) ----
    twins, twins_meta = load_twins(twins_paths)
    env = build_env_timeline(twins, twins_meta, load_ps(ps_paths))
    print_env_summary(env)
    save_env_csv(env, OUTPUT_DIR / "env_timeline_sol1217_1218.csv")

    # ---- 외피 가압 검토 (돔 형상과 내압만으로 정해짐 -> 부재 질량 탐색 전에 1회) ----
    #      통과하면 조용히 넘어가고, 결과는 최종 보고의 외피 검토(가압 + 감압)에서 한 번에 출력
    if USE_MEMBRANE:
        ev = check_envelope(float(np.max(env["dp"])) if USE_INTERNAL_PRESSURE else 0.0, verbose=False)
        if not ev["ok"]:
            raise ValueError(f"외피가 가압을 못 버팁니다 (DCR {ev['dcr']:.3f}). 웨빙이 {ev['straps_needed_per_m']}가닥/m 필요 "
                             f"(현재 {RESTRAINT_STRAPS_PER_M}가닥 x {RESTRAINT_PLIES}겹). RESTRAINT_STRAPS_PER_M(최대 39)·"
                             f"RESTRAINT_PLIES를 늘리거나 INTERNAL_PRESSURE_PA를 낮추세요.")

    # ---- 실측 지진: 원본 그대로 -> 응답스펙트럼 (1회) ----
    seismic = build_seismic_input(quake_path)
    save_spectrum_csv(seismic, OUTPUT_DIR / f"response_spectrum_{MARS_EVENT}.csv")
    quake_info = seismic["comps"][MASS_SEARCH_DIRECTION]["info"]
    quake_lmst_h = utc_to_lmst_hours(env, quake_info["peak_utc"])
    quake_lmst_start = utc_to_lmst_hours(env, quake_info["start_utc"])
    actual_idx = env_indices_at_local_time(env, quake_lmst_start, quake_info["duration_s"])
    print(f"[mars] {quake_info['id']} (스케일링 없음) 최대가속도 시각 {quake_info['peak_utc']:%Y-%m-%d %H:%M:%S} UTC "
          f"= 지방시 {format_lmst(quake_lmst_h, env['meta']['sol0'])} "
          f"-> 타임라인에서 같은 지방시 구간 {len(actual_idx)}개를 '실제 시각' 조합에 사용")

    # ---- 탐색 (D x B 격자, 지진 제외) -> 양 끝 설계 ----
    base_dir = OUTPUT_DIR
    if FOUNDATION_MODE == "ring":
        df, ends = search_ring_cover(seismic, env, quake_actual_idx=actual_idx)
        targets, seen = [], set()
        for name in FINAL_ENDS:
            d, s = ends[name]
            if d in seen:
                print(f"[final] {END_LABELS[name]}은 다른 끝과 같은 설계 ({d.label()}) -> 한 번만 해석")
                continue
            seen.add(d)
            targets.append((name, d, s))
    else:
        targets = [("none", Design(), None)]

    finals = {}
    for name, design, s in targets:
        OUTPUT_DIR = base_dir / name if FOUNDATION_MODE == "ring" else base_dir
        print("\n" + "#" * 70 + f"\n# 최종 설계 [{END_LABELS.get(name, name)}] {design.label()}\n" + "#" * 70)
        if FOUNDATION_MODE == "none":
            final_search = solve_bare(seismic, env, quake_actual_idx=actual_idx)
        elif SEARCH_WITH_QUAKE:
            final_search = s
        else:
            final_search = confirm_with_quake(design, s, seismic, env, quake_actual_idx=actual_idx)
        finals[name] = _final_report(final_search, env, seismic, quake_info, quake_lmst_h)
    OUTPUT_DIR = base_dir

    # ---- 양 끝 설계 비교 ----
    if len(finals) > 1:
        print("\n========== 양 끝 설계 비교 (같은 해석 기준) ==========")
        for name, f in finals.items():
            print(f"[{END_LABELS.get(name, name)}] {f['design']}: 돔 부재 {f['dome_ton']:.1f} t ({f['dome_kg_per_m3']:.2f} kg/m^3), "
                  f"철근 {f['rebar_ton']:.1f} t, 링 콘크리트 {f['ring_concrete_ton']:,.0f} t, 덮개 {f['cover_ton']:,.0f} t, "
                  f"침하 {f['settle_max_mm']:.1f} mm, 지진 여유배율 {f['quake_margin']:.3g}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_DIR / "ends_summary.json", "w", encoding="utf-8") as fp:
        json.dump({"dome": DOME_INFO, "member_material": MEMBER_MATERIAL_INFO, "envelope": ENVELOPE_INFO, "ends": finals},
                  fp, ensure_ascii=False, indent=2, default=_json_default)
    print(f"[summary] 최종 설계 요약 저장: {(OUTPUT_DIR / 'ends_summary.json').resolve()}")
    return finals


def _final_report(final_search, env, seismic, quake_info, quake_lmst_h):
    """최종 설계 하나에 대한 보고 + 자기점검 + 그림 (결과는 현재 OUTPUT_DIR에 저장)"""
    final_case = final_search["final_case"]
    final_design = final_case["design"]
    print_final_report(final_case, env, quake_info, quake_lmst_h)
    env_chk = None
    if USE_MEMBRANE:
        print("\n========== 외피 (가압 + 최종 덮개 두께로 감압) ==========")
        env_chk = check_envelope(float(np.max(env["dp"])) if USE_INTERNAL_PRESSURE else 0.0, final_design.cover_t)
    verify_superposition(final_case, env)
    lin = verify_quake_linearity(final_case, env)
    fat = fatigue_upper_bound(final_case, seismic)
    th = None
    if VERIFY_TIME_HISTORY:
        th = verify_with_time_history(final_case, seismic, env)
    print_seismic_report(final_case, seismic, fat)
    print_load_breakdown(final_case, env)

    summary = seismic_summary_dict(final_case, seismic, fat, lin, th)
    summary["envelope"] = {**ENVELOPE_INFO, "check": env_chk}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_DIR / "seismic_eqs_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2, default=_json_default)
    print(f"[summary] 보고서용 핵심 수치 저장: {(OUTPUT_DIR / 'seismic_eqs_summary.json').resolve()}")

    T1 = final_case["seismic_checks"].get("T1") if final_case["seismic_checks"] else None
    plot_response_spectrum(seismic, OUTPUT_DIR / f"response_spectrum_{MARS_EVENT}.png", T1=T1)
    plot_env_timeline(env, final_case["env_dcr_t"], quake_lmst_h, OUTPUT_DIR / "env_timeline.png")
    th1 = th.get(MASS_SEARCH_DIRECTION) if th else None
    if SAVE_SOL_TIMELINE:
        save_sol_timeline(OUTPUT_DIR / f"sol_timeline_{MARS_EVENT}.npz", final_case["env_e"], final_case["F_env"],
                          gens=final_case["quake_gens"],
                          th=({"force_paths": th1["force_paths"], "F_base": th1["F_base"]} if th1 else None))
    plot_dome_3d_abaqus(final_case["dome"], dome_NODES, DOME_EDGES, final_case["member_stress"],
                        section_info=final_case["section_info"], design=final_design,
                        assessment=final_case["assessment"], save_html_path=OUTPUT_DIR / "dome_stress.html")

    mp, ring, ground = final_case["mass_parts_kg"], final_case["ring"], final_case["ground"]
    margin = final_case["quake_margin"] or {}
    return {"design": final_design.label(), "D_m": final_design.ring_depth, "B_m": final_design.ring_b,
            "cover_t_m": final_design.cover_t, "flange_t_m": final_design.flange_t,
            "dome_ton": mp["steel"] / 1000.0, "dome_kg_per_m3": mp["steel"] / DOME_ENCLOSED_VOLUME_M3,
            "rebar_ton": (ring["rebar_kg"] + final_design.flange_rebar_kg) / 1000.0,
            "ring_concrete_ton": mp["ring"] / 1000.0, "cover_ton": mp["cover"] / 1000.0,
            "settle_max_mm": ground.get("settle", {}).get("max_m", 0.0) * 1000.0,
            "angular_distortion": ground.get("settle", {}).get("angular", 0.0),
            "envelope_ton": mp["membrane"] / 1000.0, "envelope_dcr": env_chk["dcr"] if env_chk else None,
            "max_dcr": final_case["max_dcr"], "quake_margin": margin.get("overall", np.inf)}


if __name__ == "__main__":
    main()
