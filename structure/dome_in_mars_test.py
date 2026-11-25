"""
화성 지오데식 돔 — 링 기초 + 레골리스 덮개 + 인장 불가 지반 + 실측 하중(S1222a, TWINS, PS)

입력 데이터 (MARS_DATA 폴더 하나에 모두 모음)
  - S1222a_ACC.mseed                : SEIS(ELYSE) S1222a 지반가속도 (m/s^2, BHN/BHE/BHZ, 20 Hz)
  - twins_calib_1217_01.csv / 1218 : TWINS 보정 데이터 (붐 막대온도, 수평풍속, 풍향, 1 Hz)
  - ps_calib_1217_01.csv / 1218    : PS 보정 데이터 (대기압, 10 Hz)
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
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
STEEL_DENSITY = 7850.0  # kg/m^3, 일반구조용 압연강재

# ----------------------------------------------------------------------
# 데이터 경로 — 환경변수 MARS_DATA_DIR이 있으면 그 경로를 우선 사용
# ----------------------------------------------------------------------
DATA_DIR = Path(os.environ.get("MARS_DATA_DIR", Path.home() / "Desktop" / "MARS_DATA"))
QUAKE_FILE = "S1222a_ACC.mseed"
TWINS_FILES = ("twins_calib_1217_01.csv", "twins_calib_1218_01.csv")
PS_FILES = ("ps_calib_1217_01.csv", "ps_calib_1218_01.csv")
OUTPUT_DIR = Path("results_mars_dome")

# ----------------------------------------------------------------------
# 화성 실측 지진 설정 (S1222a 원본 그대로 사용, 배율/PGA 스케일링 없음)
# ----------------------------------------------------------------------
MARS_EVENT = "S1222a"           # InSight 최대 화성지진 (2022-05-04)
DIR_TO_CHANNEL = {1: "BHN", 2: "BHE", 3: "BHZ"}   # X=북남, Y=동서, Z=수직
QUAKE_DT = 0.05                 # mseed 원본 20 Hz와 동일 -> 리샘플링 없음

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

USE_THERMAL = True
STEEL_ALPHA = 1.2e-5             # 1/K, 강재 열팽창계수
THERMAL_REF_TEMP_K = None        # 무응력 기준온도. None이면 타임라인 평균 기온 사용

USE_MEMBRANE = True
MEMBRANE_NAME = "ETFE 3겹"
MEMBRANE_DENSITY = 1750.0           # kg/m^3 (ETFE)
MEMBRANE_LAYER_THICKNESS = 0.2e-3   # m, 1겹 두께 (가정)
MEMBRANE_LAYERS = 3
MEMBRANE_AREAL_DENSITY = MEMBRANE_DENSITY * MEMBRANE_LAYER_THICKNESS * MEMBRANE_LAYERS  # kg/m^2

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
    1: (1385.640646, -800.0, 1200.0),
    2: (800.0, -1385.640646, 1200.0),
    3: (435.89, -754.983627, 1800.0),
    4: (871.78, 0.0, 1800.0),
    5: (1600.0, 0.0, 1200.0),
    6: (0.0, 0.0, 2000.0),
    7: (-435.89, -754.983627, 1800.0),
    8: (0.0, -1600.0, 1200.0),
    9: (493.79568, -1842.870565, 600.0),
    10: (0.0, -2000.0, 0.0),
    11: (1000.0, -1732.050808, 0.0),
    12: (1732.050808, -1000.0, 0.0),
    13: (1842.870565, -493.79568, 600.0),
    14: (1385.640646, 800.0, 1200.0),
    15: (800.0, 1385.640646, 1200.0),
    16: (0.0, 1600.0, 1200.0),
    17: (435.89, 754.983627, 1800.0),
    18: (1349.074886, -1349.074886, 600.0),
    19: (-871.78, 0.0, 1800.0),
    20: (-435.89, 754.983627, 1800.0),
    21: (-493.79568, -1842.870565, 600.0),
    22: (-800.0, -1385.640646, 1200.0),
    23: (-1349.074886, -1349.074886, 600.0),
    24: (-1000.0, -1732.050808, 0.0),
    25: (2000.0, 0.0, 0.0),
    26: (1842.870565, 493.79568, 600.0),
    27: (493.79568, 1842.870565, 600.0),
    28: (1000.0, 1732.050808, 0.0),
    29: (1349.074886, 1349.074886, 600.0),
    30: (1732.050808, 1000.0, 0.0),
    31: (-800.0, 1385.640646, 1200.0),
    32: (-1385.640646, 800.0, 1200.0),
    33: (-1600.0, 0.0, 1200.0),
    34: (0.0, 2000.0, 0.0),
    35: (-1000.0, 1732.050808, 0.0),
    36: (-1349.074886, 1349.074886, 600.0),
    37: (-1842.870565, 493.79568, 600.0),
    38: (-1385.640646, -800.0, 1200.0),
    39: (-1842.870565, -493.79568, 600.0),
    40: (-2000.0, 0.0, 0.0),
    41: (-1732.050808, -1000.0, 0.0),
    42: (-1732.050808, 1000.0, 0.0),
    43: (-493.79568, 1842.870565, 600.0),
}
"""(x,y,z)"""

dome_SCALE = 1 / 100  # 4000m급 돔 -> 40m급 돔으로 (반지름 2000m -> 20m)
dome_NODES = {tag: (x * dome_SCALE, y * dome_SCALE, z * dome_SCALE)
              for tag, (x, y, z) in dome_NODES.items()}

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

ELE_TAGS = [ele for ele, _, _ in DOME_EDGES]   # 부재력 배열의 행 순서


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


def _orient_triangles_outward(nodes, triangles, center):
    """각 삼각형의 법선이 center에서 바깥을 향하도록 정점 순서를 보정한다."""
    center = np.asarray(center, dtype=float)
    fixed = []
    for a, b, c in triangles:
        v0, v1, v2 = (np.asarray(nodes[k], dtype=float) for k in (a, b, c))
        normal = np.cross(v1 - v0, v2 - v0)
        centroid = (v0 + v1 + v2) / 3.0
        fixed.append((a, c, b) if np.dot(normal, centroid - center) < 0 else (a, b, c))
    return fixed


def _mesh_signed_volume(nodes, triangles):
    """발산정리: V = (1/6) * sum(v0 . (v1 x v2))"""
    total = 0.0
    for a, b, c in triangles:
        v0, v1, v2 = (np.asarray(nodes[k], dtype=float) for k in (a, b, c))
        total += np.dot(v0, np.cross(v1, v2))
    return total / 6.0


def compute_enclosed_volume(nodes, edges):
    """돔 쉘(삼각망)이 감싸는 내부 체적(m^3).
    밑단 링(z=0)을 중심 가상노드로 부채꼴 삼각형화해 닫은 뒤 발산정리로 계산한다."""
    shell_tris = _find_triangles_from_edges(nodes, edges)

    base_tags = [t for t, (x, y, z) in nodes.items() if abs(z) < 1e-6]
    base_tags.sort(key=lambda t: np.arctan2(nodes[t][1], nodes[t][0]))

    CENTER_TAG = -1  # 실제 절점 태그와 겹치지 않는 가상 노드
    cx = float(np.mean([nodes[t][0] for t in base_tags]))
    cy = float(np.mean([nodes[t][1] for t in base_tags]))
    nodes_ext = dict(nodes)
    nodes_ext[CENTER_TAG] = (cx, cy, 0.0)

    n = len(base_tags)
    cap_tris = [(CENTER_TAG, base_tags[i], base_tags[(i + 1) % n]) for i in range(n)]

    zs = [z for (x, y, z) in nodes.values()]
    ref_center = (cx, cy, (min(zs) + max(zs)) / 2.0)
    oriented = _orient_triangles_outward(nodes_ext, shell_tris + cap_tris, ref_center)
    return abs(_mesh_signed_volume(nodes_ext, oriented))


DOME_ENCLOSED_VOLUME_M3 = compute_enclosed_volume(dome_NODES, DOME_EDGES)


# ----------------------------------------------------------------------
# 1. 행성 파라미터
# ----------------------------------------------------------------------
PLANETS = {
    "Earth": {"g": 9.81, "label": "지구"},
    "Mars":  {"g": 3.71, "label": "화성"},
}
CURRENT_PLANET = "Mars"   # 실측 지진·환경 데이터가 모두 화성(InSight)이므로 화성 고정


def set_planet(name: str):
    global CURRENT_PLANET
    if name not in PLANETS:
        raise ValueError(f"알 수 없는 행성: {name} (Earth / Mars 중 선택)")
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
        "pga_g_local": float(abs(acc[i_pk]) / planet_g()),
        "duration_s": float(t[-1]),
    }
    return t, acc, info


# ----------------------------------------------------------------------
# 2. 기초 · 덮개 · 지반 설정
# ----------------------------------------------------------------------
FOUNDATION_MODE = "ring"   # "ring"  : 링 기초 + 인장 불가 지반 (화성 현실 모델)
                           # "legacy": 이전 탄성 스프링 (인장 무한 저항, 이전 결과 재현·비교용)

# 2-1. 원지반 레골리스 (링이 묻히는 흙)
SITE_REGOLITH_DENSITY = 1300.0    # kg/m^3, InSight 착륙지 추정 약 1.3 Mg/m^3
SITE_PHI_DEG = 30.0               # 내부마찰각 (가정, 민감도 대상)
SITE_COHESION_PA = 0.0            # InSight 착륙지: 무점착 ~ 수 kPa 미만 -> 0 (보수적)
SITE_UPLIFT_K = 1.0               # 인발 미끄럼면 토압계수 (가정)
SITE_SUBGRADE_K = 2.0e7           # N/m^3, 지반반력계수 (가정 -> 팀원 지반 코드 값으로 교체)
SITE_H_TO_V_STIFFNESS = 0.5       # 수평/수직 지반 강성비 (가정)
SITE_BASE_FRICTION_RATIO = 2.0 / 3.0   # 바닥 마찰각 δ = (2/3)φ (관례)
SITE_REGOLITH_THICKNESS_M = 5.0   # 레골리스 두께 (InSight 3~5 m 이상 추정). 링 깊이 상한
UPLIFT_MOBILIZE_M = 0.02          # 플랜지 위 흙 인발 저항이 다 발휘되는 변위 (가정)

# 2-2. 레골리스 덮개
COVER_DENSITY = 1500.0            # kg/m^3, 쌓은 레골리스 (로드맵 범위 1200~1600)
COVER_MIN_THICKNESS_M = 2.0       # 방사선 차폐 하한 (Zhang et al. 2022: 1~2 m -> 상한값)
COVER_INERTIA = True              # 덮개 질량을 지진 관성에 포함 (보수적). False면 무게만 하중으로
COVER_T_STEP_M = 0.1              # t_min 탐색 간격
COVER_T_MAX_M = 15.0              # t_min 탐색 상한
DIURNAL_SKIN_DEPTH_M = 0.05       # 일교차 열침투 깊이 δ = sqrt(κP/π), κ ≈ 3e-8 ~ 1e-7 m^2/s (가정)
                                  # 덮개 아래 부재 ΔT = 기온 ΔT x exp(-t/δ)

# 2-3. 링 기초 (황 콘크리트, 역T 단면). 링 윗면 = 지표 = 돔 밑단
RING_N_SUB = 4                    # 밑단 노드 사이 원호 등분 수 (12 x 4 = 48개 링 노드)
RING_STEM_B = 1.0                 # m, 기둥부 폭 (돔 밑단이 앉는 부분)
RING_FLANGE_T = 0.5               # m, 바닥판 두께
RING_B_GRID = tuple(float(b) for b in np.round(np.arange(2.0, 20.01, 0.5), 3))   # 바닥판 폭 후보 (지지력으로 자동 선택)

RING_MATERIAL_CASES = {
    #        f_c (Pa)  f_t (Pa)   E (Pa)   rho     근거
    "낙관": (50.0e6,  5.0e6,    2.0e10,  2300.0),  # Wan, Wendner, Cusatis 2016 (지구 실험실)
    "기준": (25.0e6,  2.5e6,    2.0e10,  2300.0),  # Snehal et al. 2023 (27 MPa), f_t = f_c/10
    "비관": (3.7e6,   0.37e6,   9.6e9,   2300.0),  # Tute & Goulas 2024 (CO2 조건), f_t 가정
}
RING_MATERIAL = "기준"
RING_NU = 0.2                     # 가정 (일반 콘크리트 관례)
RING_REBAR_FY = 500.0e6           # Pa, 링 인장 보강 철근 (지구 반입)

# 2-4. 기초 안전율
FS_UPLIFT = 1.5
FS_BEARING = 3.0
FS_SLIDING = 1.5
BEARING_STEEL_ALLOWANCE_TON = 600.0   # 사전 지지력 검토용 돔 강재 추정치 (실제 값은 해석에서 재확인)

# 2-5. 탐색 설정
RING_DEPTH_GRID = (1.0, 2.0, 3.0, 4.0, 5.0)   # 링 깊이 D 후보 (RING_FLANGE_T < D <= 레골리스 두께)
COVER_CHECK_OFFSETS = (0.0, 1.0)              # t_min, t_min+1 m -> '덮개는 얇을수록 유리' 검증
SEARCH_WITH_QUAKE = False                     # True면 탐색 단계부터 지진 포함 (매우 느림)
FOUNDATION_RETRY_STEP_M = 0.5                 # 해석의 절점별 기초 검토에서 탈락하면 덮개를 이만큼 올려 재시도
FOUNDATION_RETRY_MAX = 6

# 2-6. 태그 범위 (돔 절점 1001~1043, 돔 부재 1~114와 겹치지 않게)
RING_NODE_BASE = 2000
GROUND_NODE_BASE = 3000
RING_ELE_BASE = 4000
GROUND_ELE_BASE = 5000
GROUND_MAT_BASE = 100000          # 절점 i의 재료 = GROUND_MAT_BASE + 10*i + (1~8)

DEFAULT_SPRING_STIFFNESS = {      # legacy 모드 전용 (이전 코드와 동일)
    1: 5.0e7, 2: 5.0e7, 3: 8.0e7,   # kx, ky, kz
    4: 1.0e6, 5: 1.0e6, 6: 5.0e5,   # rx, ry, rz
}


@dataclass(frozen=True)
class Design:
    """기초·덮개 설계값. 부재 질량은 이분 탐색이 따로 정한다."""
    ring_depth: float = 0.0   # m, 지표 ~ 링 바닥 (D)
    cover_t: float = 0.0      # m, 덮개 두께 (t)
    ring_b: float = 0.0       # m, 바닥판 폭 (B, 지지력으로 자동 산정)

    def label(self):
        if FOUNDATION_MODE != "ring":
            return f"legacy_t{self.cover_t:.1f}"
        return f"D{self.ring_depth:.1f}_t{self.cover_t:.1f}_B{self.ring_b:.1f}"


class FoundationInfeasible(Exception):
    """부재 질량을 키워도 해결되지 않는 기초 탈락 (들림, 지지력, 활동, 링 압괴)."""

    def __init__(self, design, reason):
        super().__init__(f"{design.label()}: {reason}")
        self.design = design
        self.reason = reason


def get_base_node_tags(nodes, tol=1e-6):
    """z=0인 노드(돔 밑단 링)의 원래 nodeTag 목록을 자동으로 찾는다."""
    return [tag for tag, (x, y, z) in nodes.items() if abs(z) < tol]


def dome_center_radius(nodes):
    """밑면 중심 (cx, cy, 0)과 구 반지름 (전 절점까지 거리 평균)"""
    P = np.array(list(nodes.values()), dtype=float)
    base = P[np.abs(P[:, 2]) < 1e-6]
    c = np.array([base[:, 0].mean(), base[:, 1].mean(), 0.0])
    return c, float(np.linalg.norm(P - c, axis=1).mean())


# ----------------------------------------------------------------------
# 2-7. 링 기초 형상 · 단면 · 지반 저항
# ----------------------------------------------------------------------
def ring_layout(nodes, n_sub=RING_N_SUB):
    """돔 밑단 노드를 각도 순으로 잇고, 사이 원호를 n_sub 등분한 링 절점 배치.
    반환: points [(x, y, 돔밑단태그 or None)], L_trib(각 절점이 맡는 원호 길이), R, arc_length, chord_length"""
    c, _ = dome_center_radius(nodes)
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


def ring_ground_per_m(D, B, cover_t):
    """링 단위길이(1 m)당 콘크리트 무게, 플랜지 위 흙 무게, 인발 저항, 활동 저항 (N/m).

    인발: 플랜지 위 흙 기둥 무게(양쪽) + 덮개 상재압(바깥쪽만) + 양끝 연직 미끄럼면 마찰
          (안쪽은 돔 내부 바닥이라 덮개 없음. 내압이 바닥을 누르는 효과는 무시 -> 보수적)
    활동: 안쪽 면 수동토압(덮개 없음, 보수적) + (링 + 플랜지 위 흙) 무게 x 바닥 마찰"""
    sec = ring_section(D, B)
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
MEMBER_T_D_RATIO = 0.04           # 단면 형상비 t/D (구조별로 동일하게 유지) 두께와 외경의 비

# 자동 최소질량 탐색 설정
MASS_SEARCH_START_TON = 200.0          # 시작 질량 (SAFE면 자동으로 아래로 내려가며 FAIL을 찾는다)
MASS_SEARCH_COARSE_STEP_TON = 100.0    # 브라케팅 고정 증감폭
MASS_SEARCH_TOL_TON = 5.0              # 최종 질량 해상도 (ton)
MASS_SEARCH_MAX_TON = 5000.0           # 무한 반복 방지용 탐색 상한
MASS_SEARCH_DIRECTION = 1              # 1=X(북남), 2=Y(동서), 3=Z
DISP_LIMIT_RATIO = None                # 예: 1/400 -> 구조변형 최대 수평변위/스팬 초과 시 FAIL (None이면 보고만)


def total_member_length(nodes, edges):
    """모든 부재의 중심선 길이 합계(m)."""
    return sum(float(np.linalg.norm(np.asarray(nodes[nj], float) - np.asarray(nodes[ni], float)))
               for _, ni, nj in edges)


def size_uniform_shs_for_target_mass(nodes, edges, target_mass_ton,
                                     t_over_D=MEMBER_T_D_RATIO, density=STEEL_DENSITY):
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


def steel_shs_section(D, t, E=2.0e11, G=None):
    """정사각형 강관(SHS) 물성 (모서리 라운드 무시)
        A = D^2 - (D-2t)^2,  I = (D^4 - (D-2t)^4)/12,  J = t(D-t)^3 (Bredt),  c = D/2"""
    if G is None:
        G = E / (2 * 1.3)  # 포아송비 0.3
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


def compute_tributary_mass(nodes, edges, node_tags, A, density=STEEL_DENSITY, extra_mass_per_node=0.0):
    """부재 자중(단면적 x 길이 x 밀도)을 양 끝 절점에 절반씩 분배해서 절점질량을 구한다."""
    mass = defaultdict(float)
    for _, ni, nj in edges:
        L = np.linalg.norm(np.asarray(nodes[nj], float) - np.asarray(nodes[ni], float))
        m_ele = A * L * density
        mass[node_tags[ni]] += m_ele / 2
        mass[node_tags[nj]] += m_ele / 2
    if extra_mass_per_node:
        for stag in node_tags.values():
            mass[stag] += extra_mass_per_node
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
    steel_mass = A * DOME_TOTAL_MEMBER_LENGTH * STEEL_DENSITY
    model = {"node_tags": node_tags, "base_node_tags": base_node_tags, "elements": elements,
             "mass_dict": mass_dict, "section": (A, E, G, J, Iy, Iz), "c": section["c"],
             "steel_mass_kg": steel_mass}

    if FOUNDATION_MODE == "ring":
        model["ground"] = _build_ring_and_ground(model, design, contact)
    else:
        model["ground"] = _build_legacy_ground(model)

    if verbose:
        gr = model["ground"]
        extra = (f", 링 {len(gr['ring_eles'])}구간 (D={design.ring_depth:.1f} m, B={design.ring_b:.1f} m, "
                 f"{gr['ring_mass_kg']/1000:,.0f} ton)") if gr["mode"] == "ring" else ", 기존 탄성 스프링"
        print(f"[model] 돔 절점 {len(node_tags)}개, 부재 {len(elements)}개 (SHS {steel_mass/1000:.1f} ton)"
              f"{extra}, 지반 지지점 {len(gr['fixed'])}개 ({'선형' if contact is not None else '비선형'})")
    return model


def _build_ring_and_ground(model, design, contact):
    lay = RING_LAYOUT
    gp = ring_ground_per_m(design.ring_depth, design.ring_b, design.cover_t)
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


def _build_legacy_ground(model):
    """이전 코드와 같은 탄성 스프링 (인장 무한 저항). 비교·검증 전용."""
    k = DEFAULT_SPRING_STIFFNESS
    for dof, kval in k.items():
        ops.uniaxialMaterial("Elastic", GROUND_MAT_BASE + dof, kval)
    fixed, springs = [], []
    for i, (tag, stag) in enumerate(model["base_node_tags"].items()):
        x, y, _ = dome_NODES[tag]
        f_tag = GROUND_NODE_BASE + i + 1
        ops.node(f_tag, x, y, 0.0)
        ops.fix(f_tag, 1, 1, 1, 1, 1, 1)
        s_tag = GROUND_ELE_BASE + i + 1
        ops.element("zeroLength", s_tag, f_tag, stag,
                    "-mat", *[GROUND_MAT_BASE + d for d in k], "-dir", *list(k.keys()))
        fixed.append(f_tag)
        springs.append(s_tag)
    return {"mode": "legacy", "fixed": fixed, "springs": springs, "ring_nodes": [], "ring_eles": [],
            "ring_mass_kg": 0.0}


# ----------------------------------------------------------------------
# 5. 하중 -> 절점력 변환 (모두 {원래 nodeTag: np.array([Fx, Fy, Fz])} 형태로 반환)
# ----------------------------------------------------------------------
def _oriented_shell_panels(nodes, edges):
    """외향 정렬된 삼각패널 [(a, b, c, 외향 면적벡터, 중심높이 z), ...]"""
    tris = _find_triangles_from_edges(nodes, edges)
    cx = float(np.mean([p[0] for p in nodes.values()]))
    cy = float(np.mean([p[1] for p in nodes.values()]))
    tris = _orient_triangles_outward(nodes, tris, (cx, cy, 0.0))   # 반구형: 밑면 중심 기준 외향
    panels = []
    for a, b, c in tris:
        v0, v1, v2 = (np.asarray(nodes[k], dtype=float) for k in (a, b, c))
        area_vec = 0.5 * np.cross(v1 - v0, v2 - v0)
        if np.linalg.norm(area_vec) < 1e-12:
            continue
        panels.append((a, b, c, area_vec, float((v0[2] + v1[2] + v2[2]) / 3.0)))
    return panels


def _sum_forces(*force_dicts):
    total = defaultdict(lambda: np.zeros(3))
    for fd in force_dicts:
        for k, f in fd.items():
            total[k] = total[k] + np.asarray(f, dtype=float)
    return dict(total)


def compute_membrane_node_mass(node_tags, nodes, edges, areal_density=MEMBRANE_AREAL_DENSITY, verbose=True):
    """막 자중(면적 x 면밀도)을 패널 세 절점에 1/3씩 분배. 키는 구조물 절점 태그."""
    mass = defaultdict(float)
    for a, b, c, area_vec, _ in _oriented_shell_panels(nodes, edges):
        m = areal_density * float(np.linalg.norm(area_vec))
        for k in (a, b, c):
            mass[node_tags[k]] += m / 3.0
    if verbose:
        print(f"[membrane] {MEMBRANE_NAME} (면밀도 {areal_density:.2f} kg/m^2) 총 질량 {sum(mass.values())/1000:.2f} ton")
    return dict(mass)


def compute_cover_node_mass(node_tags, nodes, edges, t, rho=COVER_DENSITY, verbose=True):
    """레골리스 덮개: 패널 면적 x t x 곡률보정 x rho 를 세 절점에 1/3씩. 하중은 apply_gravity가 수직 하향으로 건다."""
    if t <= 0.0:
        return {}
    _, R = dome_center_radius(nodes)
    f = cover_curvature_factor(t, R)
    mass = defaultdict(float)
    for a, b, c, area_vec, _ in _oriented_shell_panels(nodes, edges):
        m = rho * t * f * float(np.linalg.norm(area_vec))
        for k in (a, b, c):
            mass[node_tags[k]] += m / 3.0
    if verbose:
        tot = sum(mass.values())
        print(f"[cover] 덮개 t={t:.2f} m, rho={rho:.0f} kg/m^3, 곡률보정 {f:.3f} -> "
              f"{tot/1000:,.0f} ton ({tot*planet_g()/1e6:,.1f} MN)")
    return dict(mass)


def pressure_node_forces(nodes, edges, p_net):
    """순압 p_net(+면 바깥으로 미는 힘)을 패널 외향 방향으로 걸고 세 절점에 1/3씩 분배."""
    force = defaultdict(lambda: np.zeros(3))
    for a, b, c, area_vec, _ in _oriented_shell_panels(nodes, edges):
        for k in (a, b, c):
            force[k] = force[k] + p_net * area_vec / 3.0
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
    w = np.array([math.cos(az), math.sin(az), 0.0])
    force = defaultdict(lambda: np.zeros(3))
    for a, b, c, area_vec, zc in _oriented_shell_panels(nodes, edges):
        n_hat = area_vec / np.linalg.norm(area_vec)
        cos_t = float(np.dot(n_hat, w))
        cp = -WIND_CP_WINDWARD * cos_t if cos_t < 0.0 else WIND_CP_SUCTION
        f_tri = -q_ref * wind_height_factor(zc) * cp * area_vec    # 압력(+)은 안쪽, 흡입(-)은 바깥쪽
        for k in (a, b, c):
            force[k] = force[k] + f_tri / 3.0
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


def correct_thermal_forces(local_forces, axial_force):
    """localForce [FxI ... | FxJ ...]에 열하중 고정단력 보정. N_실제 = N_해석 - EAαΔT
    (N = -FxI = FxJ 이므로 FxI += EAαΔT, FxJ -= EAαΔT)"""
    f = np.array(local_forces, dtype=float, copy=True)
    f[..., 0] += axial_force
    f[..., 6] -= axial_force
    return f


def apply_nodal_forces(node_tags, forces, pattern_tag, ts_tag=1):
    ops.pattern("Plain", pattern_tag, ts_tag)
    for tag, f in forces.items():
        ops.load(node_tags[tag], float(f[0]), float(f[1]), float(f[2]), 0.0, 0.0, 0.0)


def element_local_forces(ele_tags=ELE_TAGS):
    """(부재수, 12) = [FxI FyI FzI MxI MyI MzI | FxJ FyJ FzJ MxJ MyJ MzJ]"""
    return np.array([ops.eleResponse(e, "localForce") for e in ele_tags], dtype=float)


def ground_reactions(model):
    """지반 지지점 반력 (n, 3). +z = 지반이 링을 밀어 올림(눌림), -z = 플랜지 위 흙이 잡아당김(들림)."""
    ops.reactions()
    return np.array([ops.nodeReaction(f)[:3] for f in model["ground"]["fixed"]], dtype=float)


def ring_axial(model):
    """링 부재 축력 (m,). +인장"""
    eles = model["ground"]["ring_eles"]
    return np.array([ops.eleResponse(e, "localForce")[6] for e in eles], dtype=float)


# 형상에서 한 번만 계산하는 상수
RING_LAYOUT = ring_layout(dome_NODES)
DOME_SHELL_AREA_M2 = float(sum(np.linalg.norm(p[3]) for p in _oriented_shell_panels(dome_NODES, DOME_EDGES)))
DOME_UPLIFT_PER_PA = float(sum(f[2] for f in pressure_node_forces(dome_NODES, DOME_EDGES, 1.0).values()))   # m^2
DOME_SPHERE_RADIUS_M = dome_center_radius(dome_NODES)[1]


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
# 5-2. 기초 사전 설계 (정적 평형) — 부재 질량과 무관하게 D, t, B의 성립 여부를 판정
# ----------------------------------------------------------------------
def predesign_foundation(D, t, env, B_min=0.0):
    """링 깊이 D, 덮개 t에서 (B_min: 해석 검토에서 지지력이 모자라 폭을 넓혀야 할 때의 하한)
      (1) 바닥판 폭 B = 감압 사고(LC2) 지지력을 만족하는 최소 폭 (폭은 지지력이 정한다)
            (덮개 + 막 + 링 + 플랜지 위 흙 + 강재 추정) / (B x 둘레) <= q_ult / FS
      (2) 그 B에서 정상 가압(LC1, 외기압 최저 시각) 들림 검토 (들림은 덮개·깊이가 맡는다)
            덮개 + 막 + 링 + 플랜지 위 흙 인발저항 >= FS x 들림
    강재 자중은 들림 저항에서 빼서 보수적으로 둔다."""
    g = planet_g()
    Lr = RING_LAYOUT["arc_length"]
    F_up = DOME_UPLIFT_PER_PA * float(np.max(env["dp"])) if USE_INTERNAL_PRESSURE else 0.0
    W_cover = COVER_DENSITY * cover_volume_m3(t) * g
    W_mem = MEMBRANE_AREAL_DENSITY * DOME_SHELL_AREA_M2 * g if USE_MEMBRANE else 0.0
    W_steel = BEARING_STEEL_ALLOWANCE_TON * 1000.0 * g
    out = {"D": D, "t": t, "F_up": F_up, "W_cover": W_cover, "ok": False,
           "reason": f"지지력 부족 (B {RING_B_GRID[-1]:.1f} m까지)"}
    for B in RING_B_GRID:
        gp = ring_ground_per_m(D, B, t) if B >= B_min - 1e-9 else None
        if gp is None:
            continue
        W_ring, W_soil, R_up = gp["conc_w"] * Lr, gp["soil_w"] * Lr, gp["uplift"] * Lr
        q_app = (W_cover + W_mem + W_ring + W_soil + W_steel) / (B * Lr)
        q_allow = bearing_capacity_pa(B, D) / FS_BEARING
        if q_app > q_allow:
            continue
        resist = W_cover + W_mem + W_ring + R_up
        ratio = resist / max(F_up, 1e-9)
        out.update(B=B, W_ring=W_ring, R_up=R_up, resist=resist, uplift_ratio=ratio, q_app=q_app, q_allow=q_allow)
        if resist >= FS_UPLIFT * F_up:
            out.update(ok=True, reason="", design=Design(D, t, B))
        else:
            out.update(reason=f"들림 저항 부족 (저항/들림 {ratio:.2f} < FS {FS_UPLIFT})")
        return out
    return out


def adjust_design_after_failure(design, reason, env):
    """해석의 절점별 기초 검토에서 탈락했을 때 다음 설계안.
      지지력 / 링 압축 -> 같은 덮개에서 바닥판 폭 B를 한 칸 넓힌다 (덮개를 올리면 지지력은 더 나빠진다)
      들림 / 활동 / 불수렴 -> 덮개를 FOUNDATION_RETRY_STEP_M 올린다 (B는 지금 값 이상 유지)"""
    D, t, B = design.ring_depth, design.cover_t, design.ring_b
    if "지지력" in reason or "링 압축" in reason:
        nxt = [b for b in RING_B_GRID if b > B + 1e-9]
        if not nxt:
            return {"ok": False, "reason": f"바닥판 폭이 격자 상한 {RING_B_GRID[-1]:.1f} m에 도달"}
        pre = predesign_foundation(D, t, env, B_min=nxt[0])
        action = f"바닥판 폭 {B:.1f} -> {pre.get('B', nxt[0]):.1f} m"
    else:
        t_new = round(t + FOUNDATION_RETRY_STEP_M, 3)
        pre = predesign_foundation(D, t_new, env, B_min=B)
        action = f"덮개 {t:.1f} -> {t_new:.1f} m"
    pre["action"] = action
    return pre


def find_min_cover(D, env, t_start=None):
    """링 깊이 D에서 기초가 성립하는 최소 덮개 두께 (방사선 하한 이상). 없으면 (None, 마지막 판정)."""
    t = COVER_MIN_THICKNESS_M if t_start is None else max(t_start, COVER_MIN_THICKNESS_M)
    last = None
    while t <= COVER_T_MAX_M + 1e-9:
        last = predesign_foundation(D, round(t, 3), env)
        if last["ok"]:
            return round(t, 3), last
        t += COVER_T_STEP_M
    return None, last


# ----------------------------------------------------------------------
# 6. 해석 단계
# ----------------------------------------------------------------------
def apply_gravity(node_mass, extra_weight=None, ts_tag=1, pattern_tag=1, verbose=True):
    """node_mass: 질량(관성 + 중력하중), extra_weight: 관성 없이 무게만 거는 질량 (덮개 COVER_INERTIA=False)"""
    g = planet_g()
    ops.timeSeries("Linear", ts_tag)
    ops.pattern("Plain", pattern_tag, ts_tag)
    for node, mass in node_mass.items():
        ops.mass(node, mass, mass, mass, 0.0, 0.0, 0.0)
        ops.load(node, 0.0, 0.0, -mass * g, 0.0, 0.0, 0.0)
    for node, mass in (extra_weight or {}).items():
        ops.load(node, 0.0, 0.0, -mass * g, 0.0, 0.0, 0.0)
    if verbose:
        W = (sum(node_mass.values()) + sum((extra_weight or {}).values())) * g
        print(f"[gravity] 중력 하중 {W/1e6:,.2f} MN (g={g}) -> {len(node_mass)}개 절점")


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


def _run_linear_static():
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-12, 6)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("단위하중 선형 정적해석 실패")


def compute_unit_load_responses(section, design, contact, env_e):
    """LC1 기준상태 접선강성의 선형 모델에서 단위하중별 응답. 실측 하중 = 단위하중 x 실측값.
      P: 순압 +1 Pa,  T: 부재 온도 +1 K (열 보정 포함, 링은 ΔT=0),  W[s]: 측정높이 풍압 1 Pa, s번째 방위
    부재력(F), 지반 반력(R), 링 축력(N), 돔 절점변위를 모두 기록한다."""
    model = build_structure(section, design, transf_type="Linear",
                            contact=contact if FOUNDATION_MODE == "ring" else None, verbose=False)
    EA_alpha = section["E"] * section["A"] * STEEL_ALPHA
    step = 360.0 / WIND_SECTORS

    cases = []
    if USE_INTERNAL_PRESSURE:
        cases.append(("P", pressure_node_forces(dome_NODES, DOME_EDGES, 1.0)))
    if USE_THERMAL:
        cases.append(("T", thermal_node_forces(dome_NODES, DOME_EDGES, EA_alpha * 1.0)))
    if env_e.get("wind_active", USE_WIND):
        cases += [(f"W{s}", wind_node_forces(dome_NODES, DOME_EDGES, 1.0, s * step)) for s in range(WIND_SECTORS)]

    ops.timeSeries("Linear", 1)
    F, R, N, disp = {}, {}, {}, {}
    for k, (name, forces) in enumerate(cases):
        tag = 100 + k
        apply_nodal_forces(model["node_tags"], forces, pattern_tag=tag, ts_tag=1)
        _run_linear_static()
        F[name] = element_local_forces()
        R[name] = ground_reactions(model)
        N[name] = ring_axial(model)
        disp[name] = capture_node_disp(model["node_tags"])
        ops.remove("loadPattern", tag)
        ops.reset()                       # 변위·시간을 0으로 되돌려 다음 단위하중과 섞이지 않게 한다
    ops.wipe()

    def pack(d, correct=None):
        out = {"P": d.get("P"), "T": d.get("T"), "W": None}
        if correct is not None and out["T"] is not None:
            out["T"] = correct(out["T"])
        if "W0" in d:
            out["W"] = np.stack([d[f"W{s}"] for s in range(WIND_SECTORS)])
        return out

    return {
        "F": pack(F, correct=lambda f: correct_thermal_forces(f, EA_alpha * 1.0)),
        "R": pack(R),
        "N": pack(N),
        "EA_alpha": EA_alpha,
        "T_disp": disp.get("T"),                   # 부재 +1 K일 때 돔 절점변위 (절점수, 3)
        "node_order": list(model["node_tags"].keys()),
    }


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


def combine_env_forces(F_base, unit, env):
    return combine_env(F_base, unit["F"], env)


# [3단계] 고유 진동수
def eigen(num_modes=3, verbose=True):
    lam = ops.eigen(num_modes)
    freqs = [(l ** 0.5) / (2 * math.pi) for l in lam]
    if verbose:
        print(f"[seismic] 고유진동수(Hz): {[round(f, 3) for f in freqs]}")
    return freqs


# [4단계] Rayleigh 감쇠 (1, 3차 모드 기준 5%)
def setup_rayleigh_damping(xi=0.05, mode_i=1, mode_j=3):
    lam = ops.eigen(mode_j)
    wi = lam[mode_i - 1] ** 0.5
    wj = lam[mode_j - 1] ** 0.5
    a0 = 2 * xi * wi * wj / (wi + wj)
    a1 = 2 * xi / (wi + wj)
    ops.rayleigh(a0, 0.0, 0.0, a1)
    print(f"[seismic] Rayleigh 감쇠 적용: xi={xi} (a0={a0:.6f}, a1={a1:.6f})")
    return a0, a1


# [5단계] 지진파 입력
def apply_earthquake(acc, dt, direction=1, ts_tag=2, pattern_tag=2):
    ops.timeSeries("Path", ts_tag, "-dt", dt, "-values", *acc, "-factor", 1.0)
    ops.pattern("UniformExcitation", pattern_tag, direction, "-accel", ts_tag)
    print(f"[seismic] 지진 하중 등록 ({len(acc)}개 값, dt={dt})")


# [6단계] 동적 해석
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
STEEL_FY = 235.0e6  # Pa
BUCKLING_K = 1.0


def _end_forces(F, end):
    """localForce(…, 12) -> 한쪽 끝단의 (N, My, Mz). N>0 인장."""
    if end == "I":
        return np.stack([-F[..., 0], F[..., 4], F[..., 5]], axis=-1)
    return np.stack([F[..., 6], F[..., 10], F[..., 11]], axis=-1)


def _hull_vertices(P):
    """점군의 볼록껍질 꼭짓점. scipy가 없거나 퇴화(평면 등)하면 전체 점을 그대로 쓴다(느리지만 정확)."""
    P = np.asarray(P, dtype=float)
    if len(P) <= 64:
        return P
    try:
        from scipy.spatial import ConvexHull
        scale = np.abs(P).max(axis=0)
        scale[scale == 0.0] = 1.0
        return P[np.unique(ConvexHull(P / scale, qhull_options="QJ").vertices)]
    except Exception:
        return P


def _combo_extremes(E, Q, cap, max_pairs=2_000_000):
    """모든 (E[i] + Q[j]) 조합에서 최대 DCR, 최대 |응력|(부호 포함), 최대 압축/좌굴 비를 구한다."""
    Pt, Pc, Pcr, Mny, Mnz, A, Iy, Iz, c = cap
    best = {"dcr": -np.inf, "i": 0, "j": 0, "NMM": (0.0, 0.0, 0.0), "sigma": 0.0, "comp": 0.0}
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
        best["comp"] = max(best["comp"], float(max(-N.min(), 0.0)) / Pcr)
    return best


def assess_members(F_env, F_lc2, force_paths, F_base, section, c, lengths, env,
                   actual_idx=None, fy=STEEL_FY, K=BUCKLING_K):
    """부재별 파손 판정.
      - LC1 환경 단독 : 타임라인 모든 시각
      - LC1 + 지진    : 지진이 하루 중 어느 시각에 와도 견디는지 (force_paths=None이면 생략)
      - LC2 감압      : 자중 + 덮개, 내압 0
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
        cap = (Pt, Pc, p_cr, Mny, Mnz, A, Iy, Iz, c)

        dF = None
        if force_paths is not None:
            dF = _read_table(force_paths[ele])[:, 1:13] - F_base[e_idx]   # 지진에 의한 부재력 증분

        best1, best2 = None, None
        best_actual, env_max, dN_quake, sigma_any, comp_any = -np.inf, 0.0, 0.0, 0.0, 0.0
        for end in ("I", "J"):
            E_pts = _end_forces(F_env[:, e_idx, :], end)
            if dF is not None:
                Q_pts = _end_forces(dF, end)
                dN_quake = max(dN_quake, float(np.abs(Q_pts[:, 0]).max()))
                Q = np.vstack([zero, _hull_vertices(Q_pts)])   # 0 = 지진 없는 순간
            else:
                Q = zero
            b1 = _combo_extremes(E_pts, Q, cap)
            env_only = b1 if dF is None else _combo_extremes(E_pts, zero, cap)
            env_max = max(env_max, env_only["dcr"])
            b2 = _combo_extremes(_end_forces(F_lc2[e_idx][None, :], end), zero, cap)

            N, My, Mz = E_pts[:, 0], E_pts[:, 1], E_pts[:, 2]
            env_dcr_t = np.maximum(env_dcr_t, np.where(N >= 0, N / Pt, -N / Pc)
                                   + np.abs(My) / Mny + np.abs(Mz) / Mnz)
            if dF is not None and actual_idx is not None and len(actual_idx):
                best_actual = max(best_actual, _combo_extremes(E_pts[actual_idx], Q[1:], cap)["dcr"])

            for b in (b1, b2):
                if abs(b["sigma"]) > abs(sigma_any):
                    sigma_any = b["sigma"]
                comp_any = max(comp_any, b["comp"])
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
            "interaction_fail": max_dcr >= 1.0,
            "yield_fail": abs(sigma_any) >= fy,
            "buckling_fail": comp_any >= 1.0,
            "max_dcr": max_dcr,
            "max_dcr_lc1": best1["dcr"],
            "max_dcr_lc2": best2["dcr"],
            "max_dcr_env_only": env_max,
            "max_dcr_actual_time": best_actual if np.isfinite(best_actual) else None,
            "gov_case": gov_case,
            "gov_env_index": best1["i"],
            "gov_lmst_h": float(env["lmst_h"][best1["i"]]),
            "gov_with_quake": gov_case == "LC1 가압+지진",
            "gov_end": gov["end"],
            "critical_N": gov["NMM"][0],
            "critical_My": gov["NMM"][1],
            "critical_Mz": gov["NMM"][2],
            "signed_stress": sigma_any,
            "max_abs_stress": abs(sigma_any),
            "max_compression_ratio": comp_any,
            "slenderness": K * L / r,
            "p_cr": p_cr,
            "compression_capacity": Pc,
            "quake_max_dN": dN_quake,
        }
    return result, env_dcr_t


def assess_ground(model, R_env, R_lc2, R_base, dR=None):
    """지반 판정. R_env (시각수, n, 3), R_lc2/R_base (n, 3), dR = (지진 중 최소 증분, 최대 증분) 각 (n, 3).
    +z = 눌림(지반이 밀어 올림), -z = 들림(플랜지 위 흙이 잡아당김)."""
    g = model["ground"]
    Rz_env = R_env[..., 2]
    out = {"sum_Rz_lc2": float(R_lc2[:, 2].sum()), "sum_Rz_lc1": float(R_base[:, 2].sum()),
           "min_Rz": float(Rz_env.min()), "max_Rz": float(max(Rz_env.max(), R_lc2[:, 2].max()))}
    if g["mode"] != "ring":
        out.update(ok=True, uplift_dcr=0.0, bearing_dcr=0.0, sliding_dcr=0.0, n_flip=0, n_lifted=0,
                   n_h_limit=0, H_net=0.0, q_max=0.0, reason="", mode="legacy")
        return out

    dmin = dR[0] if dR is not None else np.zeros_like(R_base)
    dmax = dR[1] if dR is not None else np.zeros_like(R_base)
    Rz_min = Rz_env.min(axis=0) + dmin[:, 2]
    Rz_max = np.maximum(Rz_env.max(axis=0) + dmax[:, 2], R_lc2[:, 2])

    uplift_dcr = np.maximum(-Rz_min, 0.0) * FS_UPLIFT / g["cap_up"]
    q = (np.maximum(Rz_max, 0.0) + g["soil_w"]) / (g["B"] * g["L_trib"])
    bearing_dcr = q / g["q_allow"]
    # 활동: 돔이 바깥으로 미는 힘은 링 후프 인장이 스스로 평형을 이룬다(지반은 일부만 분담, 국부 한계 도달 허용).
    #       링 전체가 미끄러지는지는 '지반 수평 반력의 합력' vs '링 전체 수평 저항'으로 본다.
    net_env = float(np.hypot(R_env[..., 0].sum(axis=1), R_env[..., 1].sum(axis=1)).max())
    net_quake = float(np.hypot(np.maximum(np.abs(dmin[:, 0]), np.abs(dmax[:, 0])).sum(),
                               np.maximum(np.abs(dmin[:, 1]), np.abs(dmax[:, 1])).sum()))
    net_lc2 = float(np.hypot(R_lc2[:, 0].sum(), R_lc2[:, 1].sum()))
    H_net = max(net_env + net_quake, net_lc2)
    sliding_dcr = H_net * FS_SLIDING / g["cap_h"].sum()
    H_node = np.maximum(np.hypot(R_env[..., 0], R_env[..., 1]).max(axis=0), np.hypot(R_lc2[:, 0], R_lc2[:, 1]))
    n_h_limit = int(np.sum(H_node >= 0.98 * g["cap_h"]))

    contact0 = R_base[:, 2] > 0.0
    flips = np.any((Rz_env > 0.0) != contact0[None, :], axis=0)
    out.update(uplift_dcr=float(uplift_dcr.max()), bearing_dcr=float(bearing_dcr.max()),
               sliding_dcr=float(sliding_dcr), H_net=H_net, n_h_limit=n_h_limit, n_flip=int(flips.sum()),
               n_lifted=int((~contact0).sum()), q_max=float(q.max()), mode="ring")
    reasons = [name for name, v in (("들림", out["uplift_dcr"]), ("지지력", out["bearing_dcr"]),
                                    ("활동", out["sliding_dcr"])) if v >= 1.0]
    out["ok"] = not reasons
    out["reason"] = "/".join(reasons)
    return out


def assess_ring(model, N_env, N_lc2, dN=None):
    """링 축력 판정: 압축 이용률(f_c), 인장이 f_t·A를 넘으면 필요한 철근 질량."""
    g = model["ground"]
    if g["mode"] != "ring":
        return {"ok": True, "T_max": 0.0, "C_max": 0.0, "comp_dcr": 0.0, "rebar_kg": 0.0, "need_rebar": False}
    A = g["ring_sec"]["A"]
    dmin = dN[0] if dN is not None else np.zeros_like(N_lc2)
    dmax = dN[1] if dN is not None else np.zeros_like(N_lc2)
    N_hi = np.maximum(N_env.max(axis=0) + dmax, N_lc2)
    N_lo = np.minimum(N_env.min(axis=0) + dmin, N_lc2)
    T_max, C_max = float(max(N_hi.max(), 0.0)), float(max(-N_lo.min(), 0.0))
    need = T_max > g["ring_ft"] * A
    A_s = T_max / RING_REBAR_FY if need else 0.0
    rebar_kg = A_s * RING_LAYOUT["chord_length"] * STEEL_DENSITY
    comp_dcr = C_max / (g["ring_fc"] * A)
    return {"ok": comp_dcr < 1.0, "T_max": T_max, "C_max": C_max, "comp_dcr": comp_dcr,
            "need_rebar": need, "A_s": A_s, "rebar_kg": rebar_kg}


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


def _section_corners(width, height):
    hw, hh = width / 2.0, height / 2.0
    return [(-hw, -hh), (hw, -hh), (hw, hh), (-hw, hh)]


def element_box_vertices(p1, p2, width, height):
    _, ey, ez = _local_axes(p1, p2)
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    return np.array([p + oy * ey + oz * ez for p in (p1, p2) for (oy, oz) in _section_corners(width, height)])


_BOX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
_BOX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]


def _add_box(verts, value, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z):
    base = len(X)
    for v in verts:
        X.append(float(v[0])); Y.append(float(v[1])); Z.append(float(v[2]))
    for (a, b, c, d) in _BOX_FACES:
        I.append(base + a); J.append(base + b); K.append(base + c); intensity.append(value)
        I.append(base + a); J.append(base + c); K.append(base + d); intensity.append(value)
    for (a, b) in _BOX_EDGES:
        edge_x.extend([verts[a][0], verts[b][0], None])
        edge_y.extend([verts[a][1], verts[b][1], None])
        edge_z.extend([verts[a][2], verts[b][2], None])


def _add_ground_block(footing_xy, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z, depth=200.0, pad=100.0):
    xs = [p[0] for p in footing_xy]
    ys = [p[1] for p in footing_xy]
    x0, x1 = min(xs) - pad, max(xs) + pad
    y0, y1 = min(ys) - pad, max(ys) + pad
    verts = np.array([[x0, y0, 0.0], [x1, y0, 0.0], [x1, y1, 0.0], [x0, y1, 0.0],
                      [x0, y0, -depth], [x1, y0, -depth], [x1, y1, -depth], [x0, y1, -depth]])
    _add_box(verts, 0.0, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z)


def plot_dome_3d_abaqus(dome, nodes, edges, member_stress, colorscale="RdBu",
                        ground_depth=None, ground_pad=None, vis_side=None, save_html_path=None,
                        design=None, title="지오데식 돔 — 지배 조합 부재 응력"):
    if vis_side is None:
        vis_side = 2.0 * dome["c"]
    base_tags = list(dome["base_node_tags"].keys())
    footing_xy = [(nodes[t][0], nodes[t][1]) for t in base_tags]
    xs = [p[0] for p in nodes.values()]
    ys = [p[1] for p in nodes.values()]
    span = max(max(xs) - min(xs), max(ys) - min(ys))
    ground_depth = ground_depth or span * 0.1
    ground_pad = ground_pad or span * 0.05

    X, Y, Z, I, J, K, intensity = [], [], [], [], [], [], []
    edge_x, edge_y, edge_z = [], [], []
    _add_ground_block(footing_xy, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z,
                      depth=ground_depth, pad=ground_pad)
    for ele_tag, ni, nj in edges:
        verts = element_box_vertices(nodes[ni], nodes[nj], vis_side, vis_side)
        _add_box(verts, member_stress.get(ele_tag, 0.0), X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z)

    cmax = max((abs(v) for v in intensity), default=1.0) or 1.0
    fig = go.Figure()
    fig.add_trace(go.Mesh3d(
        x=X, y=Y, z=Z, i=I, j=J, k=K, intensity=intensity, intensitymode="cell",
        colorscale=colorscale, cmin=-cmax, cmax=cmax, opacity=0.95, showscale=True,
        colorbar=dict(title="Stress (Pa)<br>+인장 / -압축"), name="Dome solid",
        hovertemplate="x=%{x:.2f} m<br>y=%{y:.2f} m<br>z=%{z:.2f} m<extra></extra>"))
    fig.add_trace(go.Scatter3d(x=edge_x, y=edge_y, z=edge_z, mode="lines",
                               line=dict(color="black", width=2), name="Element edges", hoverinfo="skip"))
    if FOUNDATION_MODE == "ring" and design is not None:
        pts = RING_LAYOUT["points"] + RING_LAYOUT["points"][:1]
        fig.add_trace(go.Scatter3d(
            x=[p[0] for p in pts], y=[p[1] for p in pts], z=[-0.5 * design.ring_depth] * len(pts),
            mode="lines", line=dict(color="saddlebrown", width=10),
            name=f"링 기초 (D={design.ring_depth:.1f} m, B={design.ring_b:.1f} m)"))
    fig.update_layout(
        title=title, template="plotly_white", width=1200, height=900,
        scene=dict(xaxis_title="X 북 (m)", yaxis_title="Y 동 (m)", zaxis_title="Z (m)", aspectmode="data",
                   camera=dict(eye=dict(x=1.8, y=-1.8, z=1.1), center=dict(x=0, y=0, z=0))),
        legend=dict(x=0.02, y=0.98))
    if save_html_path:
        out_path = Path(save_html_path).resolve()
        fig.write_html(str(out_path), auto_open=False)
        print(f"[viz] 저장: {out_path}")
    fig.show()
    return fig


# ----------------------------------------------------------------------
# 9-2. 최대변위: 노드 변위 recorder / 기초 강체운동 제거 / 요약 (지진 응답)
# ----------------------------------------------------------------------
def setup_node_disp_recorder(node_tags, out_dir):
    """전 절점 (ux, uy, uz). 열 순서: time | node1(ux uy uz) | node2(...) ..."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "node_disp.txt")
    ops.recorder("Node", "-file", path, "-precision", RECORDER_PRECISION, "-time",
                 "-node", *node_tags.values(), "-dof", 1, 2, 3, "disp")
    return path


def capture_node_disp(node_tags):
    return np.array([ops.nodeDisp(s)[:3] for s in node_tags.values()], dtype=float)


def load_node_disp(path, n_nodes):
    data = _read_table(path)
    return data[:, 0], data[:, 1:].reshape(len(data), n_nodes, 3)


def remove_base_rigid_motion(X0, du, base_idx):
    """밑단 링 변위에 강체운동(병진 u0 + 회전 theta x r)을 최소제곱으로 맞춰 전 절점에서 제거."""
    R = X0 - X0[base_idx].mean(axis=0)
    M = np.zeros((len(X0), 3, 3))  # M @ theta = theta x r
    M[:, 0, 1], M[:, 0, 2] = R[:, 2], -R[:, 1]
    M[:, 1, 0], M[:, 1, 2] = -R[:, 2], R[:, 0]
    M[:, 2, 0], M[:, 2, 1] = R[:, 1], -R[:, 0]
    nb = len(base_idx)
    A = np.concatenate([np.tile(np.eye(3), (nb, 1, 1)), M[base_idx]], axis=2).reshape(3 * nb, 6)
    b = du[:, base_idx, :].reshape(len(du), 3 * nb)
    params = b @ np.linalg.pinv(A).T  # (nt, 6) = [u0, theta]
    rigid = params[:, None, :3] + np.einsum("nij,tj->tni", M, params[:, 3:])
    return du - rigid


def summarize_peak_displacement(nodes, node_tags, base_tags, t, U, U0):
    """abs: 지반 대비(U-U0). rel: abs에서 기초 병진/락킹 제거(구조 자체 변형)."""
    tags = list(node_tags.keys())
    idx = {tag: i for i, tag in enumerate(tags)}
    X0 = np.array([nodes[k] for k in tags], dtype=float)
    du = U - U0[None, :, :]
    du_rel = remove_base_rigid_motion(X0, du, [idx[b] for b in base_tags])
    span = max(np.ptp(X0[:, 0]), np.ptp(X0[:, 1]))
    height = X0[:, 2].max()
    apex_tag = max(tags, key=lambda k: nodes[k][2])

    def peak(d):
        h = np.hypot(d[..., 0], d[..., 1])
        ti, ni = np.unravel_index(int(np.argmax(h)), h.shape)
        vz = np.abs(d[..., 2])
        tj, nj = np.unravel_index(int(np.argmax(vz)), vz.shape)
        return {"h": float(h[ti, ni]), "h_node": tags[ni], "h_time": float(t[ti]),
                "uz": float(vz[tj, nj]), "uz_node": tags[nj], "uz_time": float(t[tj]),
                "apex_h": float(h[:, idx[apex_tag]].max())}

    a, r = peak(du), peak(du_rel)
    return {"abs": a, "rel": r, "span": float(span), "height": float(height), "apex_tag": apex_tag,
            "gravity_max_uz": float(np.abs(U0[:, 2]).max()), "total_max_uz": float(np.abs(U[..., 2]).max()),
            "ratio_h_abs": a["h"] / span, "ratio_h_rel": r["h"] / span, "ratio_uz_rel": r["uz"] / span}


def print_disp_summary(s):
    if s is None:
        print("[disp] 지진 해석 생략 (지진 응답 변위 없음)")
        return

    def line(name, p):
        return (f"  {name}: 수평 {p['h']*1000:.4f} mm (절점 {p['h_node']}, t={p['h_time']:.2f}s) | "
                f"수직 {p['uz']*1000:.4f} mm (절점 {p['uz_node']}, t={p['uz_time']:.2f}s) | "
                f"정점(절점 {s['apex_tag']}) 수평 {p['apex_h']*1000:.4f} mm")
    print("[disp] 지진에 의한 최대변위 (기준 정적처짐 제외)")
    print(line("지반 대비   ", s["abs"]))
    print(line("기초운동 제거", s["rel"]))
    r = s["ratio_h_rel"]
    if r > 0:
        print(f"[disp] 구조 변형 최대 수평변위/스팬 = {r:.3e} (스팬/{1/r:,.0f}) "
              f"[스팬 {s['span']:.1f} m, 높이 {s['height']:.1f} m]")
    print(f"[disp] LC1 기준 정적변위(자중+덮개+링+평균순압) 수직 최대 {s['gravity_max_uz']*1000:.2f} mm, "
          f"지진 포함 수직 최대 {s['total_max_uz']*1000:.2f} mm")


# ----------------------------------------------------------------------
# 9-3. 그림: 한글 폰트 / 환경 타임라인 / 지진 응답 영상 / 탐색 결과
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
    """링 깊이 D별 최소 덮개·링 폭, 최소 부재 질량(지구 반입 질량), 현지 작업량"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _setup_korean_font(plt)

    ok = df[df["feasible"] & df["is_min_cover"]].sort_values("D_m")
    if ok.empty:
        print("[search] 성립하는 설계가 없어 그림을 생략합니다.")
        return
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    ax = axes[0]
    ax.plot(ok["D_m"], ok["t_m"], "o-", color="tab:brown", label="최소 덮개 t")
    ax.axhline(COVER_MIN_THICKNESS_M, color="0.5", ls="--", lw=0.9, label="방사선 하한")
    ax.set_xlabel("링 깊이 D (m)"); ax.set_ylabel("덮개 두께 (m)")
    ax2 = ax.twinx()
    ax2.plot(ok["D_m"], ok["B_m"], "s:", color="tab:gray", label="바닥판 폭 B")
    ax2.set_ylabel("바닥판 폭 B (m)")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8)

    ax = axes[1]
    ax.plot(ok["D_m"], ok["steel_ton"], "o-", label="돔 강재")
    ax.plot(ok["D_m"], ok["earth_ton"], "s-", label="지구 반입 (강재 + 철근)")
    chk = df[df["feasible"] & ~df["is_min_cover"]].sort_values("D_m")
    if not chk.empty:
        ax.plot(chk["D_m"], chk["steel_ton"], "x", color="tab:red", label="덮개 +Δ 검증")
    ax.set_xlabel("링 깊이 D (m)"); ax.set_ylabel("질량 (ton)"); ax.legend(fontsize=8)

    ax = axes[2]
    ax.plot(ok["D_m"], ok["cover_volume_m3"], "o-", label="덮개 부피")
    ax.plot(ok["D_m"], ok["excavation_m3"], "s-", label="링 굴착 부피")
    ax.plot(ok["D_m"], ok["local_work_m3"], "^-", label="합계")
    ax.set_xlabel("링 깊이 D (m)"); ax.set_ylabel("현지 작업량 (m³)"); ax.legend(fontsize=8)
    for a in axes:
        a.grid(alpha=0.3)
    fig.suptitle(f"링 깊이 x 덮개 두께 탐색 (내압 {INTERNAL_PRESSURE_PA/1000:.0f} kPa, 링 재료 '{RING_MATERIAL}')")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[search] 그림 저장: {Path(out_path).resolve()}")


def animate_dome_response(nodes, edges, node_tags, t, U, U0,
                          ground_t=None, ground_acc=None, member_value=None,
                          direction=1, out_path="dome_shaking.mp4",
                          fps=20, stride=5, scale=None, dpi=90,
                          title="지진 응답", value_label="최대 DCR"):
    """stride*dt 간격으로 프레임 생성. scale=None이면 최대변위가 스팬의 8%로 보이도록 자동 확대."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    _setup_korean_font(plt)

    tags = list(node_tags.keys())
    idx = {tag: i for i, tag in enumerate(tags)}
    X0 = np.array([nodes[k] for k in tags], dtype=float)
    du = U - U0[None, :, :]
    ia = np.array([idx[ni] for _, ni, _ in edges])
    ib = np.array([idx[nj] for _, _, nj in edges])

    span = max(np.ptp(X0[:, 0]), np.ptp(X0[:, 1]))
    height = X0[:, 2].max()
    peak = float(np.abs(du).max())
    if scale is None:
        scale = float(np.clip(0.08 * span / peak, 1.0, 1e9)) if peak > 0 else 1.0

    apex = max(tags, key=lambda k: nodes[k][2])
    base_idx = [i for i, k in enumerate(tags) if abs(nodes[k][2]) < 1e-6]
    apex_u = remove_base_rigid_motion(X0, du, base_idx)[:, idx[apex], direction - 1]

    if member_value:
        vals = np.array([member_value.get(e, 0.0) for e, _, _ in edges], dtype=float)
        norm = Normalize(0.0, max(1.0, float(vals.max())))
        cmap = plt.get_cmap("turbo")
        colors = cmap(norm(vals))
    else:
        colors = "tab:blue"

    frames = np.arange(0, len(t), max(1, int(stride)))
    fig = plt.figure(figsize=(12, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.0, 1.0], hspace=0.35, wspace=0.42)
    ax = fig.add_subplot(gs[:, 0], projection="3d")
    ax_u = fig.add_subplot(gs[0, 1])
    ax_a = fig.add_subplot(gs[1, 1])

    def segments(k):
        P = X0 + scale * du[k]
        return np.stack([P[ia], P[ib]], axis=1)

    lc = Line3DCollection(segments(0), colors=colors, linewidths=2.2)
    ax.add_collection3d(lc)
    R = np.hypot(X0[:, 0], X0[:, 1]).max()
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(R * np.cos(th), R * np.sin(th), 0 * th, color="0.5", lw=1.0)
    lim = 0.55 * span
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(0, height * 1.15)
    ax.set_box_aspect((2 * lim, 2 * lim, height * 1.15), zoom=1.1)
    ax.view_init(elev=18, azim=-58)
    ax.set_xlabel("X 북 (m)"); ax.set_ylabel("Y 동 (m)"); ax.set_zlabel("Z (m)")

    p0 = X0[idx[apex]] + scale * du[0, idx[apex]]
    (apex_pt,) = ax.plot([p0[0]], [p0[1]], [p0[2]], "o", color="red", ms=8)
    if member_value:
        fig.colorbar(ScalarMappable(norm, cmap), ax=ax, location="left", shrink=0.55, pad=0.0, label=value_label)

    ax_u.plot(t, apex_u * 1000.0, color="tab:blue", lw=1.0)
    (mk_u,) = ax_u.plot([t[0]], [apex_u[0] * 1000.0], "o", color="red")
    vl_u = ax_u.axvline(t[0], color="red", lw=0.8, alpha=0.6)
    ax_u.set_title(f"정점 변위 ({'XYZ'[direction - 1]}방향, 기초운동 제거)", fontsize=10)
    ax_u.set_ylabel("mm"); ax_u.grid(alpha=0.3)

    vl_a = None
    if ground_t is not None and ground_acc is not None:
        ax_a.plot(ground_t, ground_acc, color="tab:gray", lw=0.8)
        vl_a = ax_a.axvline(t[0], color="red", lw=0.8, alpha=0.6)
    ax_a.set_title("입력 지반가속도", fontsize=10)
    ax_a.set_xlabel("t (s)"); ax_a.set_ylabel("m/s²"); ax_a.grid(alpha=0.3)

    def update(n):
        k = frames[n]
        lc.set_segments(segments(k))
        p = X0[idx[apex]] + scale * du[k, idx[apex]]
        apex_pt.set_data([p[0]], [p[1]])
        apex_pt.set_3d_properties([p[2]])
        ax.set_title(f"{title}   t = {t[k]:7.1f} s   (변형 ×{scale:,.0f})", fontsize=12)
        mk_u.set_data([t[k]], [apex_u[k] * 1000.0])
        vl_u.set_xdata([t[k], t[k]])
        if vl_a is not None:
            vl_a.set_xdata([t[k], t[k]])
        return lc, apex_pt

    anim = animation.FuncAnimation(fig, update, frames=len(frames), interval=1000 / fps, blit=False)
    ext = os.path.splitext(out_path)[1].lower()
    if ext == ".mp4" and not animation.writers.is_available("ffmpeg"):
        out_path = os.path.splitext(out_path)[0] + ".gif"
        ext = ".gif"
        print("[anim] ffmpeg를 찾지 못해 GIF로 저장합니다.")
    writer = animation.FFMpegWriter(fps=fps, bitrate=3500) if ext == ".mp4" else animation.PillowWriter(fps=fps)
    anim.save(out_path, writer=writer, dpi=dpi,
              progress_callback=lambda i, n: print(f"[anim] 프레임 {i + 1}/{n}") if (i + 1) % 200 == 0 else None)
    plt.close(fig)
    print(f"[anim] 저장: {os.path.abspath(out_path)} ({len(frames)}프레임, 변형 ×{scale:,.0f})")
    return out_path


# ----------------------------------------------------------------------
# 10. 1솔 타임라인 저장 (피로 해석 입력용)
#     환경 부재력(10분 간격) + 지진 증분(0.05 s)을 따로 저장하고, 필요할 때 합쳐서 복원한다.
# ----------------------------------------------------------------------
def _validate_quake_starts(starts, quake_len, duration):
    starts = sorted(float(s) for s in starts)
    if starts and (starts[0] < 0.0 or starts[-1] + quake_len > duration):
        raise ValueError(f"지진 구간이 타임라인(0~{duration:.0f}s)을 벗어납니다. "
                         f"(마지막 지진 {starts[-1]:.0f}~{starts[-1] + quake_len:.0f}s)")
    for a, b in zip(starts[:-1], starts[1:]):
        if a + quake_len > b:
            raise ValueError(f"지진 구간이 겹칩니다: {a:.0f}~{a + quake_len:.0f}s 와 {b:.0f}s 시작.")
    return starts


def save_sol_timeline(path, env, F_env, force_paths, F_base, quake_starts=SOL_QUAKE_STARTS_S):
    """npz(zip) 구조
        t_env, lmst_env, T_K, P_Pa, q_Pa, dT, dp : 환경 타임라인 (dT, q_Pa는 덮개 반영 후 값)
        t_quake, quake_starts_s                  : 지진 증분 시간축, 지진 시작 시각
        env_ele_<tag>   (시각수, 12)              : 환경 부재력 (localForce)
        quake_ele_<tag> (지진스텝, 12) float32    : 지진에 의한 부재력 증분
    복원은 reconstruct_member_history() 사용."""
    t_quake = _read_table(next(iter(force_paths.values())))[:, 0]
    starts = _validate_quake_starts(quake_starts, float(t_quake[-1]), env["duration_s"])

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        def put(name, arr):
            with zf.open(f"{name}.npy", "w", force_zip64=True) as f:
                np.lib.format.write_array(f, np.asarray(arr))

        put("t_env", env["t_s"]); put("lmst_env", env["lmst_h"])
        for k in ("T_K", "P_Pa", "q_Pa", "dT", "dp"):
            put(k, env[k])
        put("t_quake", t_quake); put("quake_starts_s", np.array(starts))
        for e_idx, ele in enumerate(ELE_TAGS):
            put(f"env_ele_{ele}", F_env[:, e_idx, :])
            put(f"quake_ele_{ele}", (_read_table(force_paths[ele])[:, 1:13] - F_base[e_idx]).astype(np.float32))

    s0 = env["meta"]["sol0"]
    spans = ", ".join(f"{format_lmst(np.interp(s, env['t_s'], env['lmst_h']), s0)}" for s in starts)
    print(f"[sol] 1솔 부재력 타임라인 저장: {path.resolve()} (지진 {len(starts)}회: {spans})")
    return path


def reconstruct_member_history(npz_path, ele_tag):
    """저장된 타임라인에서 부재 하나의 전체 이력 (t, localForce(nt,12))을 복원한다.
    지진 구간: 환경(선형보간) + 지진 증분, 그 밖: 환경 부재력."""
    d = np.load(npz_path)
    t_env, F_env = d["t_env"], d[f"env_ele_{ele_tag}"]
    t_q, dF = d["t_quake"], d[f"quake_ele_{ele_tag}"].astype(float)

    t_parts, f_parts, cursor = [], [], -np.inf
    for st in d["quake_starts_s"]:
        m = (t_env > cursor) & (t_env < st)
        t_parts.append(t_env[m]); f_parts.append(F_env[m])
        tq = st + t_q
        env_at_q = np.column_stack([np.interp(tq, t_env, F_env[:, k]) for k in range(12)])
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
    return mass, extra, parts


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


def run_case(target_mass_ton, design, ground_motion, env, include_quake=True, dt=QUAKE_DT, direction=1,
             quake_actual_idx=None, verbose=True):
    """질량 하나 x 설계(D, t, B) 하나의 전체 해석.
    ground_motion=(t, acc): S1222a 실측 기록 그대로(배율 없음), env: build_env_timeline() 결과(덮개 반영 전)."""
    if target_mass_ton <= 0.0:
        raise ValueError("target_mass_ton은 0보다 커야 합니다.")
    if direction not in (1, 2, 3):
        raise ValueError("direction은 1(X), 2(Y), 3(Z) 중 하나여야 합니다.")
    t, acc = ground_motion
    if include_quake and not math.isclose(t[1] - t[0], dt, rel_tol=1e-6):
        raise ValueError(f"실측 지진 dt({t[1]-t[0]})와 해석 dt({dt})가 다릅니다.")

    section, section_info = make_section_for_target_mass(target_mass_ton)
    env_e = effective_env(env, design)
    g = planet_g()

    # ---- A) 비선형 모델: 1단계 자중(+막+덮개+링) = LC2 감압, 2단계 평균 순압 = LC1 기준 ----
    model = build_structure(section, design, verbose=verbose)
    masses, extra_weight, mass_parts = assemble_masses(model, design, verbose=verbose)
    apply_gravity(masses, extra_weight, verbose=verbose)
    if run_static_robust(verbose=verbose, label="LC2 (자중+덮개+링, 내압 0)") != 0:
        status, reason = _classify_static_failure(model)
        ops.wipe()
        return _failed_static_case(target_mass_ton, design, status, "LC2: " + reason, section_info)
    F_lc2, N_lc2, R_lc2 = element_local_forces(), ring_axial(model), ground_reactions(model)
    W_total = (sum(masses.values()) + sum(extra_weight.values())) * g

    F_up = 0.0
    if USE_INTERNAL_PRESSURE:
        F_up = DOME_UPLIFT_PER_PA * env["dp_ref"]
        apply_nodal_forces(model["node_tags"], pressure_node_forces(dome_NODES, DOME_EDGES, env["dp_ref"]),
                           pattern_tag=4)
        if verbose:
            print(f"[pressure] 기준 순압 {env['dp_ref']:,.0f} Pa -> 들림 {F_up/1e6:,.1f} MN "
                  f"(무게 {W_total/1e6:,.1f} MN)")
        if run_static_robust(verbose=verbose, label="LC1 (평균 순압 가압)") != 0:
            status, reason = _classify_static_failure(model)
            ops.wipe()
            return _failed_static_case(target_mass_ton, design, status, "LC1: " + reason, section_info)
    F_base, N_base, R_base = element_local_forces(), ring_axial(model), ground_reactions(model)
    contact = R_base[:, 2] > 0.0
    eq_err = abs(R_base[:, 2].sum() - (W_total - F_up)) / max(W_total, F_up, 1.0)
    if verbose:
        print(f"[equilibrium] 지반 반력 합 {R_base[:, 2].sum()/1e6:,.2f} MN = 무게 {W_total/1e6:,.2f} - "
              f"들림 {F_up/1e6:,.2f} MN (오차 {eq_err:.1e}) | 들린 지지점 {int((~contact).sum())}/{len(contact)}")
    freqs = eigen(3, verbose=verbose) if include_quake else []

    # ---- B) 지진: LC1 기준상태 위에 S1222a (비선형 지반 그대로) ----
    force_paths, disp, t_disp, U, U0, dR, dN = None, None, None, None, None, None, None
    if include_quake:
        out_dir = OUTPUT_DIR / "mass_search" / design.label() / f"mass_{target_mass_ton:g}ton".replace(".", "_") / f"dir_{direction}"
        force_paths = setup_member_force_recorders(model["elements"], out_dir=str(out_dir))
        disp_path = setup_node_disp_recorder(model["node_tags"], str(out_dir))
        env_paths = setup_foundation_envelope_recorders(model, str(out_dir))
        U0 = capture_node_disp(model["node_tags"])
        print(f"[seismic] 실측 지진 {MARS_EVENT}: {t[-1]:.1f}s, PGA={np.max(np.abs(acc)):.3e} m/s^2 (스케일링 없음)")
        setup_rayleigh_damping(xi=0.05, mode_i=1, mode_j=3)
        apply_earthquake(acc, dt=dt, direction=direction)
        dynamic_ok = run_dynamic_robust(num_steps=len(acc) - 1, dt=dt)
        ops.wipe()   # recorder 파일을 닫아 내용을 모두 디스크에 기록
        if not dynamic_ok:
            status = "foundation_fail" if FOUNDATION_MODE == "ring" else "member_fail"
            return _failed_static_case(target_mass_ton, design, status,
                                       "지진 해석 불수렴 (들림·미끄러짐 또는 수치 문제)", section_info)
        t_disp, U = load_node_disp(disp_path, len(model["node_tags"]))
        disp = summarize_peak_displacement(dome_NODES, model["node_tags"], model["base_node_tags"], t_disp, U, U0)
        if verbose:
            print_disp_summary(disp)
        Renv = _read_envelope(env_paths["R"]).reshape(3, -1, 3)
        dR = (Renv[0] - R_base, Renv[1] - R_base)
        if "N" in env_paths:
            Nenv = _read_envelope(env_paths["N"]).reshape(3, -1, 12)[:, :, 6]
            dN = (Nenv[0] - N_base, Nenv[1] - N_base)
    else:
        ops.wipe()
    disp_ok = disp is None or DISP_LIMIT_RATIO is None or disp["ratio_h_rel"] <= DISP_LIMIT_RATIO

    # ---- C) 단위하중 (LC1 기준상태 접선강성의 선형 모델) ----
    unit = compute_unit_load_responses(section, design, contact, env_e)

    # ---- D) 환경 타임라인 (선형 중첩) ----
    F_env = combine_env(F_base, unit["F"], env_e)
    R_env = combine_env(R_base, unit["R"], env_e)
    N_env = combine_env(N_base, unit["N"], env_e)

    # ---- E) 판정 ----
    assessment, env_dcr_t = assess_members(
        F_env, F_lc2, force_paths, F_base, model["section"], model["c"],
        get_member_lengths(dome_NODES, DOME_EDGES), env_e, actual_idx=quake_actual_idx)
    ground = assess_ground(model, R_env, R_lc2, R_base, dR)
    ring = assess_ring(model, N_env, N_lc2, dN)

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
    s0 = env["meta"]["sol0"]
    if verbose:
        where = format_lmst(gov["gov_lmst_h"], s0) if gov["gov_case"].startswith("LC1") else "-"
        print(f"[check] 부재 최대 DCR {max_dcr:.4f} (부재 {gov_ele}, {gov['gov_case']}, {where}) | "
              f"LC1 환경 단독 {max(v['max_dcr_env_only'] for v in assessment.values()):.4f} | "
              f"LC2 감압 {max(v['max_dcr_lc2'] for v in assessment.values()):.4f}"
              + (f" | 실제 지방시 지진 {max(actual_vals):.4f}" if actual_vals else ""))
        if ground["mode"] == "ring":
            print(f"[check] 지반 DCR: 들림 {ground['uplift_dcr']:.3f}, 지지력 {ground['bearing_dcr']:.3f}, "
                  f"활동 {ground['sliding_dcr']:.3f} | 접촉상태 변화 {ground['n_flip']}곳, "
                  f"국부 수평 한계 {ground['n_h_limit']}곳 | "
                  f"링 압축 이용률 {ring['comp_dcr']:.3f}, 최대 인장 {ring['T_max']/1e6:.2f} MN"
                  f"{' -> 철근 ' + format(ring['rebar_kg']/1000, '.1f') + ' ton' if ring['need_rebar'] else ''}")

    status = "ok" if (members_ok and foundation_ok and disp_ok) else (
        "member_fail" if not (members_ok and disp_ok) else "foundation_fail")
    return {
        "target_mass_ton": float(target_mass_ton),
        "design": design,
        "status": status,
        "reason": "; ".join(reasons),
        "actual_steel_mass_ton": model["steel_mass_kg"] / 1000.0,
        "mass_parts_kg": mass_parts,
        "direction": direction,
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
        "W_total_N": W_total,
        "F_up_N": F_up,
        "equilibrium_error": eq_err,
        "env_dcr_t": env_dcr_t,
        "force_paths": force_paths,
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
        "disp": disp,
        "disp_ok": disp_ok,
        "disp_t": t_disp,
        "disp_U": U,
        "disp_U0": U0,
        "ground_t": t,
        "ground_acc": acc,
        "include_quake": include_quake,
    }


def verify_superposition(case, env, k=None):
    """선형 중첩 검증: LC1 지배 시각 k의 실측 하중을 비선형 지반 모델에 직접 걸어 푼 부재력과
    '기준 + 단위하중 x 실측값'으로 계산한 부재력을 비교한다. (PDelta + ENT/인발한계 지반 포함)"""
    design = case["design"]
    env_e = effective_env(env, design)
    if k is None:
        k = int(np.argmax(case["env_dcr_t"]))
    model = build_structure(case["section"], design, verbose=False)
    masses, extra, _ = assemble_masses(model, design, verbose=False)
    apply_gravity(masses, extra, verbose=False)
    if run_static_robust(verbose=False) != 0:
        raise RuntimeError("검증용 직접해석(LC2 단계)이 수렴하지 않았습니다.")

    parts = []
    if USE_INTERNAL_PRESSURE:
        parts.append(pressure_node_forces(dome_NODES, DOME_EDGES, env_e["dp"][k]))
    if env_e["wind_active"]:
        parts.append(wind_node_forces(dome_NODES, DOME_EDGES, env_e["q_Pa"][k],
                                      env_e["sector"][k] * 360.0 / WIND_SECTORS))
    EA_aT = case["unit"]["EA_alpha"] * env_e["dT"][k] if USE_THERMAL else 0.0
    if USE_THERMAL:
        parts.append(thermal_node_forces(dome_NODES, DOME_EDGES, EA_aT))
    if parts:
        apply_nodal_forces(model["node_tags"], _sum_forces(*parts), pattern_tag=5)
        if run_static_robust(verbose=False) != 0:
            raise RuntimeError("검증용 직접해석(LC1 단계)이 수렴하지 않았습니다.")
    F_direct = correct_thermal_forces(element_local_forces(), EA_aT)
    ops.wipe()

    F_sup = case["F_env"][k]
    err = float(np.abs(F_direct - F_sup).max() / np.abs(F_direct).max())
    print(f"[verify] 선형 중첩 vs 비선형 직접해석 ({format_lmst(env['lmst_h'][k], env['meta']['sol0'])}): "
          f"최대 상대오차 {err:.2e}")
    return err


def find_minimum_safe_mass(design, ground_motion, env,
                           include_quake=True,
                           start_mass=MASS_SEARCH_START_TON,
                           coarse_step=MASS_SEARCH_COARSE_STEP_TON,
                           tol=MASS_SEARCH_TOL_TON,
                           direction=MASS_SEARCH_DIRECTION,
                           max_mass=MASS_SEARCH_MAX_TON,
                           dt=QUAKE_DT, quake_actual_idx=None):
    """설계(D, t, B) 하나에서 고정폭 브라케팅 + 이분법으로 최소 SAFE 부재 질량을 찾는다.
      - 시작 질량이 SAFE면 아래로, FAIL이면 위로 브라케팅한다.
      - 부재는 통과했는데 기초가 탈락하면 FoundationInfeasible (부재 질량을 키워도 해결 안 됨).
    반환되는 min_safe_mass_ton은 최종 구간의 SAFE측 상한이다. (max_failed < 임계질량 <= min_safe)"""
    if start_mass <= 0.0 or coarse_step <= 0.0 or tol <= 0.0:
        raise ValueError("start_mass, coarse_step, tol은 0보다 커야 합니다.")
    if max_mass <= start_mass:
        raise ValueError("max_mass는 start_mass보다 커야 합니다.")
    for name, val in (("start_mass", start_mass), ("coarse_step", coarse_step)):
        if not math.isclose(val / tol, round(val / tol)):
            raise ValueError(f"{name}는 {tol:.1f} ton 단위여야 합니다.")

    cache, history = {}, []
    tag = "지진 포함" if include_quake else "지진 제외"

    def evaluate(mass_ton, stage):
        mass_ton = round(float(mass_ton), 9)
        if mass_ton in cache:
            return cache[mass_ton]
        print(f"\n[{stage} | {design.label()} | {tag}] 질량 {mass_ton:.1f} ton 해석 시작")
        case = run_case(mass_ton, design, ground_motion, env, include_quake=include_quake, dt=dt,
                        direction=direction, quake_actual_idx=quake_actual_idx)
        members_ok = case["members_ok"] and case["disp_ok"]
        if members_ok and not case["foundation_ok"]:
            raise FoundationInfeasible(design, case["reason"] or "기초 탈락")
        safe = bool(case["passed"])
        history.append({"mass_ton": mass_ton, "safe": safe, "status": case["status"],
                        "failed_count": case["failed_count"], "max_dcr": case["max_dcr"],
                        "max_dcr_env_only": case["max_dcr_env_only"], "stage": stage})
        verdict = "SAFE" if safe else (f"FAIL ({case['reason']})" if case["reason"]
                                        else f"FAIL (파손 부재 {case['failed_count']}개)")
        print(f"[결과] 질량={mass_ton:.1f} ton -> {verdict} (최대 DCR={case['max_dcr']:.4f})")
        if not safe:   # 탈락 케이스의 큰 배열은 바로 버려 메모리 절약
            case = {k: case[k] for k in ("max_dcr", "failed_members", "passed", "status", "reason")}
        cache[mass_ton] = (safe, case)
        return cache[mass_ton]

    if evaluate(start_mass, "브라케팅")[0]:
        upper_safe, lower_failed = float(start_mass), 0.0
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
    print(f"\n========== 브라케팅 완료: FAIL {lower_failed:.1f} ton / SAFE {upper_safe:.1f} ton ==========")

    low_index, high_index = int(round(lower_failed / tol)), int(round(upper_safe / tol))
    iteration = 0
    while high_index - low_index > 1:
        iteration += 1
        middle_index = (low_index + high_index) // 2
        if evaluate(middle_index * tol, f"이분법 {iteration}")[0]:
            high_index = middle_index
        else:
            low_index = middle_index
        print(f"[범위 갱신] {low_index * tol:.1f} ton FAIL < M_min <= {high_index * tol:.1f} ton SAFE")

    lower_failed, upper_safe = low_index * tol, high_index * tol
    final_safe, final_case = evaluate(upper_safe, "최종 확인")
    lower_safe = evaluate(lower_failed, "최종 확인")[0] if lower_failed > 0 else False
    if not final_safe or lower_safe:
        raise RuntimeError("최종 FAIL-SAFE 경계 검증이 일관되지 않습니다.")

    print(f"\n========== 최소 SAFE 질량 ({design.label()}, {tag}): {upper_safe:.1f} ton "
          f"(최대 FAIL {lower_failed:.1f} ton, 이분법 {iteration}회) ==========")
    return {"min_safe_mass_ton": float(upper_safe), "max_failed_mass_ton": float(lower_failed),
            "interval_ton": float(upper_safe - lower_failed), "tolerance_ton": float(tol),
            "iterations": iteration, "history": history, "final_case": final_case,
            "design": design, "include_quake": include_quake}


def _search_row(D, t_min, offset, pre, s):
    fc = s["final_case"]
    design = s["design"]
    ring = fc["ring"]
    ground = fc["ground"]
    gov = fc["assessment"][fc["gov_ele"]]
    rebar_t = ring["rebar_kg"] / 1000.0
    exc = design.ring_b * design.ring_depth * RING_LAYOUT["arc_length"]
    vol = cover_volume_m3(design.cover_t)
    return {
        "D_m": D, "t_min_m": t_min, "t_offset_m": offset, "t_m": design.cover_t, "B_m": design.ring_b,
        "is_min_cover": offset == 0.0, "feasible": True,
        "t_governs": "방사선" if t_min <= COVER_MIN_THICKNESS_M + 1e-9 else "구조(들림)",
        "steel_ton": s["min_safe_mass_ton"], "max_failed_ton": s["max_failed_mass_ton"],
        "rebar_ton": rebar_t, "earth_ton": s["min_safe_mass_ton"] + rebar_t,
        "ring_concrete_ton": fc["mass_parts_kg"]["ring"] / 1000.0, "cover_ton": fc["mass_parts_kg"]["cover"] / 1000.0,
        "cover_volume_m3": vol, "excavation_m3": exc, "local_work_m3": vol + exc,
        "max_dcr": fc["max_dcr"], "gov_case": gov["gov_case"], "gov_ele": fc["gov_ele"],
        "uplift_dcr": ground["uplift_dcr"], "bearing_dcr": ground["bearing_dcr"],
        "sliding_dcr": ground["sliding_dcr"], "n_flip": ground["n_flip"], "n_lifted": ground["n_lifted"],
        "n_h_limit": ground["n_h_limit"],
        "ring_comp_dcr": ring["comp_dcr"], "ring_T_max_MN": ring["T_max"] / 1e6,
        "static_uplift_ratio": pre.get("uplift_ratio", np.nan), "quake_included": s["include_quake"],
        "note": "",
    }


def search_ring_cover(ground_motion, env, quake_actual_idx=None,
                      depth_grid=RING_DEPTH_GRID, offsets=COVER_CHECK_OFFSETS, include_quake=SEARCH_WITH_QUAKE):
    """링 깊이 D 격자 x (t_min, t_min+Δ)에서 최소 부재 질량 탐색. 결과 표와 최적 설계를 반환한다."""
    rows, searches = [], {}
    print("\n" + "=" * 70)
    print(f"[search] 링 깊이 {list(depth_grid)} m, 덮개 = 최소 덮개 + {list(offsets)} m, "
          f"{'지진 포함' if include_quake else '지진 제외 (최적 설계만 2단계에서 지진 확인)'}")
    print(f"[search] 링 재료 '{RING_MATERIAL}', FS 들림 {FS_UPLIFT}/지지력 {FS_BEARING}/활동 {FS_SLIDING}, "
          f"방사선 하한 {COVER_MIN_THICKNESS_M} m, 레골리스 두께 {SITE_REGOLITH_THICKNESS_M} m")
    print("=" * 70)

    for D in depth_grid:
        if D <= RING_FLANGE_T or D > SITE_REGOLITH_THICKNESS_M + 1e-9:
            rows.append({"D_m": D, "feasible": False, "is_min_cover": True,
                         "note": f"범위 밖 (바닥판 두께 {RING_FLANGE_T} m < D <= 레골리스 {SITE_REGOLITH_THICKNESS_M} m)"})
            continue
        t_min, pre = find_min_cover(D, env)
        if t_min is None:
            rows.append({"D_m": D, "feasible": False, "is_min_cover": True,
                         "note": f"덮개 {COVER_T_MAX_M} m까지 기초 불성립 ({pre['reason']})"})
            print(f"[search] D={D:.1f} m: 덮개 {COVER_T_MAX_M} m까지 기초 불성립 ({pre['reason']})")
            continue
        print(f"\n[search] D={D:.1f} m: 정적 평형 기준 최소 덮개 {t_min:.1f} m, 링 폭 B={pre['B']:.1f} m "
              f"(저항/들림 {pre['uplift_ratio']:.2f}, 지지 압력 {pre['q_app']/1e3:,.0f}/{pre['q_allow']/1e3:,.0f} kPa)")

        for off in offsets:
            pre_t = pre if off == 0.0 else predesign_foundation(D, round(t_min + off, 3), env)
            s, notes = None, []
            for _ in range(FOUNDATION_RETRY_MAX + 1):
                if not pre_t["ok"]:
                    notes.append(f"사전설계 불성립({pre_t['reason']})")
                    break
                try:
                    s = find_minimum_safe_mass(pre_t["design"], ground_motion, env, include_quake=include_quake,
                                               quake_actual_idx=quake_actual_idx)
                    break
                except FoundationInfeasible as exc:
                    pre_t = adjust_design_after_failure(exc.design, exc.reason, env)
                    notes.append(f"{exc.design.label()} 기초 탈락({exc.reason}) -> {pre_t['action']}")
                    print(f"[search] {notes[-1]}")
            if s is None:
                rows.append({"D_m": D, "t_min_m": t_min, "t_offset_m": off, "feasible": False,
                             "is_min_cover": off == 0.0, "note": " | ".join(notes) or "재시도 한도 초과"})
                continue
            note = " | ".join(notes)
            row = _search_row(D, t_min, off, pre, s)
            row["note"] = note
            rows.append(row)
            searches[s["design"]] = s

    df = pd.DataFrame(rows)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / "ring_cover_search.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"\n[search] 결과 표 저장: {csv_path.resolve()}")

    ok = df[df["feasible"]] if "feasible" in df else df.iloc[0:0]
    if ok.empty:
        raise RuntimeError("성립하는 (링 깊이, 덮개) 조합이 없습니다. 링 폭 격자, 덮개 상한, 내압 설정을 확인하세요.")

    # '덮개는 얇을수록 유리' 검증 (같은 D에서 t_min+Δ가 더 가벼우면 경고)
    for D, grp in ok.groupby("D_m"):
        base = grp[grp["is_min_cover"]]
        thick = grp[~grp["is_min_cover"]]
        if not base.empty and not thick.empty and thick["steel_ton"].min() < base["steel_ton"].iloc[0] - 1e-9:
            print(f"[search] 경고: D={D:.1f} m에서 덮개를 두껍게 하니 부재가 더 가벼워짐 "
                  f"({base['steel_ton'].iloc[0]:.1f} -> {thick['steel_ton'].min():.1f} ton). "
                  f"이 D에서는 덮개 두께도 격자로 탐색해야 합니다.")
        elif not base.empty and not thick.empty:
            print(f"[search] 확인: D={D:.1f} m에서 덮개 +Δ -> 부재 {base['steel_ton'].iloc[0]:.1f} -> "
                  f"{thick['steel_ton'].min():.1f} ton (최소 덮개가 유리)")

    best_row = ok.sort_values(["earth_ton", "local_work_m3"]).iloc[0]
    best = next(d for d in searches if math.isclose(d.ring_depth, best_row["D_m"])
                and math.isclose(d.cover_t, best_row["t_m"]) and math.isclose(d.ring_b, best_row["B_m"]))
    cols = ["D_m", "t_m", "B_m", "t_governs", "steel_ton", "rebar_ton", "earth_ton", "gov_case",
            "uplift_dcr", "bearing_dcr", "ring_comp_dcr", "local_work_m3", "note"]
    print("\n" + df[[c for c in cols if c in df]].to_string(index=False, float_format=lambda v: f"{v:,.2f}"))
    print(f"\n[search] 최적: {best.label()} -> 지구 반입 {best_row['earth_ton']:.1f} ton "
          f"(강재 {best_row['steel_ton']:.1f} + 철근 {best_row['rebar_ton']:.1f})")
    plot_search_results(df, OUTPUT_DIR / "ring_cover_search.png")
    return df, best, searches[best]


def confirm_with_quake(design, search_no_quake, ground_motion, env, quake_actual_idx=None):
    """2단계: 지진 제외로 찾은 최적 설계를 지진 포함으로 다시 탐색한다.
    지진은 하중만 더하므로 지진 제외 최대 FAIL 질량부터 위로 올라가면 된다.
    지진 중 기초가 탈락하면 덮개를 올려 다시 시도한다."""
    tol = MASS_SEARCH_TOL_TON
    start = max(search_no_quake["max_failed_mass_ton"], tol)
    step = max(tol, 4 * tol)
    for _ in range(FOUNDATION_RETRY_MAX + 1):
        try:
            s = find_minimum_safe_mass(design, ground_motion, env, include_quake=True, start_mass=start,
                                       coarse_step=step, quake_actual_idx=quake_actual_idx)
            add = s["min_safe_mass_ton"] - search_no_quake["min_safe_mass_ton"]
            print(f"[confirm] 지진 포함 최소 질량 {s['min_safe_mass_ton']:.1f} ton "
                  f"(지진 제외 대비 {add:+.1f} ton)")
            return s
        except FoundationInfeasible as exc:
            pre = adjust_design_after_failure(exc.design, exc.reason, env)
            print(f"[confirm] 지진 포함 시 기초 탈락 ({exc.reason}) -> {pre.get('action', '')}")
            if not pre["ok"]:
                raise RuntimeError(f"기초 사전설계 불성립: {pre['reason']}") from exc
            design = pre["design"]
    raise RuntimeError("지진 포함 확인에서 재시도 한도를 넘었습니다.")


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
          f"(덮개 아래 ΔT 배율 {env_e['thermal_factor']:.2e}, 바람 {'반영' if env_e['wind_active'] else '차폐'})")
    if FOUNDATION_MODE == "ring":
        print(f"[design] 링 깊이 D={d.ring_depth:.1f} m, 바닥판 폭 B={d.ring_b:.1f} m, 덮개 t={d.cover_t:.1f} m "
              f"(링 재료 '{RING_MATERIAL}', 원지반 φ={SITE_PHI_DEG:.0f}°, ρ={SITE_REGOLITH_DENSITY:.0f})")
    else:
        print("[design] legacy 탄성 스프링 (인장 무한 저항, 비교용)")
    print(f"[section] D={si['D']:.4f} m, t={si['t']:.4f} m, A={si['A']:.6f} m^2")
    print(f"[mass] 지구 반입: 돔 강재 {mp['steel']/1000:,.1f} t + 링 철근 {ring['rebar_kg']/1000:,.1f} t "
          f"= {(mp['steel'] + ring['rebar_kg'])/1000:,.1f} t | 현지: 링 황 콘크리트 {mp['ring']/1000:,.0f} t, "
          f"덮개 {mp['cover']/1000:,.0f} t ({cover_volume_m3(d.cover_t):,.0f} m^3) | 막 {mp['membrane']/1000:.2f} t")
    print(f"[equilibrium] 무게 {case['W_total_N']/1e6:,.2f} MN, 들림(평균 순압) {case['F_up_N']/1e6:,.2f} MN, "
          f"지반 반력 합 LC2 {ground['sum_Rz_lc2']/1e6:,.2f} / LC1 {ground['sum_Rz_lc1']/1e6:,.2f} MN "
          f"(평형 오차 {case['equilibrium_error']:.1e})")
    where = (f" at {format_lmst(gov['gov_lmst_h'], s0)}" if gov["gov_case"].startswith("LC1") else "")
    print(f"[check] 지배 부재 {case['gov_ele']} ({gov['gov_end']}단): DCR {gov['max_dcr']:.4f}, "
          f"{gov['gov_case']}{where}")
    print(f"        N={gov['critical_N']/1e3:,.1f} kN, My={gov['critical_My']/1e3:,.2f} kN·m, "
          f"Mz={gov['critical_Mz']/1e3:,.2f} kN·m, 세장비 {gov['slenderness']:.1f}")
    print(f"[check] 케이스별 최대 DCR: LC1 환경 단독 {case['max_dcr_env_only']:.4f}, "
          f"LC2 감압 {case['max_dcr_lc2']:.4f}, 전체 {case['max_dcr']:.4f}")
    if case["max_dcr_actual_time"] is not None:
        print(f"[check] 실제 지진 지방시의 환경 + 지진 최대 DCR {case['max_dcr_actual_time']:.4f}")
    if case["include_quake"]:
        dN = max(v["quake_max_dN"] for v in a.values())
        print(f"[check] 지진에 의한 최대 축력 변동 {dN:,.1f} N (기준 상태 대비)")
    if ground["mode"] == "ring":
        print(f"[foundation] DCR 들림 {ground['uplift_dcr']:.3f}, 지지력 {ground['bearing_dcr']:.3f} "
              f"(최대 {ground['q_max']/1e3:,.0f} kPa), 활동 {ground['sliding_dcr']:.3f} "
              f"(수평 합력 {ground['H_net']/1e3:,.1f} kN) | LC1 들린 지지점 {ground['n_lifted']}곳, "
              f"하루 중 접촉상태 변화 {ground['n_flip']}곳, 국부 수평 한계 도달 {ground['n_h_limit']}곳 "
              f"(링 후프가 분담)")
        print(f"[ring] 최대 인장 {ring['T_max']/1e6:,.2f} MN, 최대 압축 {ring['C_max']/1e6:,.2f} MN "
              f"(압축 이용률 {ring['comp_dcr']:.3f}) -> "
              + (f"철근 {ring['A_s']*1e4:,.0f} cm^2 ({ring['rebar_kg']/1000:,.1f} t)" if ring["need_rebar"]
                 else "인장이 황 콘크리트 인장강도 이내 (철근 불필요)"))
        if ground["n_flip"]:
            print("[warn] 하루 중 접촉 상태가 바뀌는 지지점이 있어 선형 중첩 오차가 커질 수 있습니다. "
                  "verify_superposition 결과를 확인하세요.")


def print_thermal_effect(case, env):
    """기온 변화에 따른 부재 팽창·수축과 그로 인한 부재력 요약 (덮개 반영 ΔT)"""
    if not USE_THERMAL or case["unit"]["T_disp"] is None:
        return
    env_e = case["env_e"]
    s0, lm, dT = env["meta"]["sol0"], env["lmst_h"], env_e["dT"]
    k_cold, k_hot = int(np.argmin(dT)), int(np.argmax(dT))
    unit = case["unit"]
    print("\n========== 온도 (기온에 따른 팽창·수축) ==========")
    print(f"[thermal] 기온 {env['T_K'].min():.1f} ~ {env['T_K'].max():.1f} K, 무응력 기준 {env['T_ref']:.1f} K, "
          f"덮개 {case['design'].cover_t:.1f} m 아래 부재 ΔT 배율 exp(-t/δ) = {env_e['thermal_factor']:.2e} "
          f"(δ = {DIURNAL_SKIN_DEPTH_M} m)")
    if env_e["thermal_factor"] < 1e-6:
        print("[thermal] 덮개가 일교차를 사실상 차단 -> 부재 열응력 0. 아래는 부재 +1 K당 구속 정도(참고)")
    else:
        print(f"[thermal] 부재 ΔT {dT[k_cold]:+.2f} K ({format_lmst(lm[k_cold], s0)}) ~ "
              f"{dT[k_hot]:+.2f} K ({format_lmst(lm[k_hot], s0)})")

    L = np.array(list(get_member_lengths(dome_NODES, DOME_EDGES).values()))
    print(f"[thermal] 부재 자유 신축 ΔL = α·ΔT·L (최장 {L.max():.2f} m): "
          f"{STEEL_ALPHA*dT[k_cold]*L.max()*1000:+.3f} ~ {STEEL_ALPHA*dT[k_hot]*L.max()*1000:+.3f} mm")

    U1, order = unit["T_disp"], unit["node_order"]
    apex = order.index(max(order, key=lambda n: dome_NODES[n][2]))
    base = [i for i, n in enumerate(order) if abs(dome_NODES[n][2]) < 1e-6]
    xy = np.array([dome_NODES[order[i]][:2] for i in base], dtype=float)
    R = np.linalg.norm(xy, axis=1)
    ur1 = float(np.mean(np.sum(U1[base, :2] * (xy / R[:, None]), axis=1)))
    ur_free = STEEL_ALPHA * float(R.mean())
    print(f"[thermal] 부재 +1 K당 밑단 반경방향 변위 {ur1*1000:+.4f} mm (구속 없으면 {ur_free*1000:+.4f} mm "
          f"-> 링·지반이 신축의 {(1 - ur1/ur_free)*100:.0f}%를 막음), 정점 수직 {U1[apex, 2]*1000:+.4f} mm")

    N1 = np.stack([_end_forces(unit["F"]["T"], e)[:, 0] for e in ("I", "J")])
    N1 = N1[np.argmax(np.abs(N1), axis=0), np.arange(N1.shape[1])]
    j = int(np.argmax(np.abs(N1)))
    Py = STEEL_FY * case["section"]["A"]
    for label, k in (("가장 추울 때", k_cold), ("가장 더울 때", k_hot)):
        N = N1[j] * dT[k]
        print(f"[thermal] {label}: 부재 {ELE_TAGS[j]} 축력 {N/1e3:+,.3f} kN, Fy·A 대비 {abs(N)/Py:.2e}")


def print_load_breakdown(case, env):
    """지배 부재의 부재력을 하중별로 분해"""
    ele = case["gov_ele"]
    gov = case["assessment"][ele]
    e, end = ELE_TAGS.index(ele), gov["gov_end"]
    if gov["gov_case"] == "LC2 감압":
        N, My, Mz = _end_forces(case["F_lc2"][e], end)
        print(f"[breakdown] 부재 {ele} {end}단, LC2 감압 (자중+막+덮개+링, 내압 0)")
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
        "순압 변동": (env_e["dp"][k] - env_e["dp_ref"]) * u["P"][e] if u["P"] is not None else np.zeros(12),
        "온도": env_e["dT"][k] * u["T"][e] if u["T"] is not None else np.zeros(12),
        "바람": env_e["q_Pa"][k] * u["W"][env_e["sector"][k], e] if u["W"] is not None else np.zeros(12),
    }
    print(f"[breakdown] 부재 {ele} {end}단, {format_lmst(env['lmst_h'][k], env['meta']['sol0'])} "
          f"(부재 ΔT {env_e['dT'][k]:+.2e} K)")
    for name, F in parts.items():
        N, My, Mz = _end_forces(F, end)
        print(f"  {name:<14}: N={N/1e3:11.1f} kN | My={My/1e3:8.2f} | Mz={Mz/1e3:8.2f} kN·m")


if __name__ == "__main__":
    set_planet(CURRENT_PLANET)
    if CURRENT_PLANET != "Mars":
        raise ValueError("실측 환경·지진 데이터가 화성(InSight) 자료이므로 CURRENT_PLANET은 'Mars'여야 합니다.")
    print(f"[volume] 돔 내부 감싸인 체적 = {DOME_ENCLOSED_VOLUME_M3:.1f} m^3, 쉘 면적 {DOME_SHELL_AREA_M2:.1f} m^2, "
          f"순압 1 Pa당 들림 {DOME_UPLIFT_PER_PA:.1f} N")
    print(f"[data] 입력 폴더: {DATA_DIR}")

    twins_paths = [DATA_DIR / f for f in TWINS_FILES]
    ps_paths = [DATA_DIR / f for f in PS_FILES]
    quake_path = DATA_DIR / QUAKE_FILE
    _require_files(*twins_paths, *ps_paths, quake_path)

    # ---- 실측 환경 타임라인 ----
    twins, twins_meta = load_twins(twins_paths)
    env = build_env_timeline(twins, twins_meta, load_ps(ps_paths))
    print_env_summary(env)
    save_env_csv(env, OUTPUT_DIR / "env_timeline_sol1217_1218.csv")

    # ---- 실측 지진: 원본 그대로 ----
    t_m, acc_m, quake_info = load_ground_motion(quake_path, DIR_TO_CHANNEL[MASS_SEARCH_DIRECTION], dt=QUAKE_DT)
    quake_lmst_h = utc_to_lmst_hours(env, quake_info["peak_utc"])
    quake_lmst_start = utc_to_lmst_hours(env, quake_info["start_utc"])
    actual_idx = env_indices_at_local_time(env, quake_lmst_start, quake_info["duration_s"])
    print(f"[mars] {quake_info['id']} PGA = {quake_info['pga']:.3e} m/s^2 "
          f"({quake_info['pga_g_local']:.3e} g, 화성 g 기준, 스케일링 없음)")
    print(f"[mars] 최대가속도 시각 {quake_info['peak_utc']:%Y-%m-%d %H:%M:%S} UTC "
          f"= 지방시 {format_lmst(quake_lmst_h, env['meta']['sol0'])} "
          f"-> 타임라인에서 같은 지방시 구간 {len(actual_idx)}개를 '실제 시각' 조합에 사용")

    # ---- 탐색 ----
    if FOUNDATION_MODE == "ring":
        df, best_design, best_search = search_ring_cover((t_m, acc_m), env, quake_actual_idx=actual_idx)
        if SEARCH_WITH_QUAKE:
            final_search = best_search
        else:
            final_search = confirm_with_quake(best_design, best_search, (t_m, acc_m), env, quake_actual_idx=actual_idx)
    else:
        final_search = find_minimum_safe_mass(Design(), (t_m, acc_m), env, include_quake=True,
                                              quake_actual_idx=actual_idx)
    final_case = final_search["final_case"]
    final_design = final_case["design"]

    print_final_report(final_case, env, quake_info, quake_lmst_h)
    print_disp_summary(final_case["disp"])
    verify_superposition(final_case, env)
    print_thermal_effect(final_case, env)
    print_load_breakdown(final_case, env)

    plot_env_timeline(env, final_case["env_dcr_t"], quake_lmst_h, OUTPUT_DIR / "env_timeline.png")
    if SAVE_SOL_TIMELINE and final_case["force_paths"] is not None:
        save_sol_timeline(OUTPUT_DIR / f"sol_timeline_{MARS_EVENT}.npz", final_case["env_e"], final_case["F_env"],
                          final_case["force_paths"], final_case["F_base"])

    plot_dome_3d_abaqus(final_case["dome"], dome_NODES, DOME_EDGES, final_case["member_stress"],
                        save_html_path=OUTPUT_DIR / "dome_stress.html", design=final_design)
    if final_case["disp_U"] is not None:
        animate_dome_response(
            dome_NODES, DOME_EDGES, final_case["dome"]["node_tags"],
            final_case["disp_t"], final_case["disp_U"], final_case["disp_U0"],
            ground_t=final_case["ground_t"], ground_acc=final_case["ground_acc"],
            member_value={k: v["max_dcr"] for k, v in final_case["assessment"].items()},
            direction=MASS_SEARCH_DIRECTION, out_path=str(OUTPUT_DIR / f"dome_shaking_{MARS_EVENT}.mp4"),
            title=f"{PLANETS[CURRENT_PLANET]['label']} {MARS_EVENT}",
            stride=40)   # dt=0.05 x 40 = 2 s/프레임 -> 50분 약 1,500프레임

    print(f"[dome] 링 기초 + 레골리스 덮개 + 실측 환경하중 + S1222a 기반 최소 SAFE 질량 탐색 완료 ({final_design.label()}).")
