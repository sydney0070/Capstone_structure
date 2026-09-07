"""
라멘(Rahmen) 구조물
"""

import os
import math
from pathlib import Path
from collections import defaultdict

import numpy as np
import openseespy.opensees as ops
import plotly.graph_objects as go


NODE_OFFSET = 1000  # 구조물 절점 태그 = 원래 nodeTag + NODE_OFFSET
STEEL_DENSITY = 7850.0  # kg/m^3, 일반구조용 압연강재


# ----------------------------------------------------------------------
# 0. 라멘 원본 데이터
# ----------------------------------------------------------------------
RAHMEN_NODES = {
    1: (0.0, 0.0, 200.0),
    2: (0.0, 0.0, 300.0),
    3: (0.0, 100.0, 300.0),
    4: (100.0, 100.0, 300.0),
    5: (100.0, 0.0, 300.0),
    6: (0.0, -100.0, 300.0),
    7: (0.0, -100.0, 200.0),
    8: (0.0, -100.0, 100.0),
    9: (0.0, -100.0, 0.0),
    10: (100.0, -100.0, 0.0),
    11: (100.0, 0.0, 0.0),
    12: (100.0, 0.0, 100.0),
    13: (100.0, 100.0, 100.0),
    14: (100.0, 100.0, 200.0),
    15: (0.0, 100.0, 200.0),
    16: (0.0, 100.0, 100.0),
    17: (-100.0, 100.0, 100.0),
    18: (-100.0, 100.0, 200.0),
    19: (-100.0, 0.0, 200.0),
    20: (100.0, 0.0, 200.0),
    21: (0.0, 0.0, 100.0),
    22: (-100.0, 0.0, 100.0),
    23: (-100.0, 0.0, 0.0),
    24: (0.0, 0.0, 0.0),
    25: (0.0, 100.0, 0.0),
    26: (100.0, 100.0, 0.0),
    27: (-100.0, 0.0, 300.0),
    28: (-100.0, -100.0, 300.0),
    29: (100.0, -100.0, 300.0),
    30: (-100.0, 100.0, 300.0),
    31: (100.0, -100.0, 200.0),
    32: (-100.0, -100.0, 200.0),
    33: (100.0, -100.0, 100.0),
    34: (-100.0, -100.0, 0.0),
    35: (-100.0, 100.0, 0.0),
    36: (-100.0, -100.0, 100.0),
}

RAHMEN_EDGES = [
    (1, 1, 2), (2, 3, 2), (3, 4, 3), (4, 4, 5), (5, 2, 5),
    (6, 2, 6), (7, 6, 7), (8, 7, 8), (9, 8, 9), (10, 10, 9),
    (11, 11, 10), (12, 12, 11), (13, 13, 12), (14, 13, 14), (15, 14, 15),
    (16, 16, 15), (17, 16, 17), (18, 18, 17), (19, 19, 18), (20, 19, 1),
    (21, 1, 7), (22, 1, 20), (23, 21, 1), (24, 21, 22), (25, 23, 22),
    (26, 24, 23), (27, 24, 25), (28, 25, 26), (29, 26, 13), (30, 13, 16),
    (31, 21, 16), (32, 27, 2), (33, 28, 27), (34, 6, 28), (35, 29, 6),
    (36, 5, 29), (37, 15, 3), (38, 15, 18), (39, 30, 18), (40, 27, 30),
    (41, 4, 14), (42, 5, 20), (43, 7, 31), (44, 32, 7), (45, 32, 19),
    (46, 22, 19), (47, 22, 17), (48, 8, 33), (49, 8, 21), (50, 12, 21),
    (51, 12, 33), (52, 10, 33), (53, 9, 34), (54, 9, 24), (55, 11, 24),
    (56, 11, 26), (57, 17, 35), (58, 15, 1), (59, 24, 21), (60, 34, 23),
    (61, 34, 36), (62, 36, 22), (63, 25, 16), (64, 35, 25), (65, 19, 27),
    (66, 28, 32), (67, 29, 31), (68, 3, 30), (69, 20, 14), (70, 20, 31),
    (71, 20, 12), (72, 32, 36), (73, 36, 8), (74, 33, 31), (75, 23, 35),
]

# ----------------------------------------------------------------------
# 0. 내부 체적 계산 (단위체적당 강재량 q=M/V 비교용)
# ----------------------------------------------------------------------
def compute_bounding_box_volume(nodes):
    """노드 좌표의 x/y/z 각 축 최대-최소 폭을 곱한 바운딩 박스 체적(m^3)을 반환한다.
    라멘처럼 삼각화된 감싸는 면이 없는 순수 골조 구조에 쓰는 대체 정의."""
    xs = [x for x, y, z in nodes.values()]
    ys = [y for x, y, z in nodes.values()]
    zs = [z for x, y, z in nodes.values()]
    Lx = max(xs) - min(xs)
    Ly = max(ys) - min(ys)
    Lz = max(zs) - min(zs)
    return float(Lx * Ly * Lz)


RAHMEN_BOUNDING_VOLUME_M3 = compute_bounding_box_volume(RAHMEN_NODES)

# ----------------------------------------------------------------------
# 1. 행성 파라미터
# ----------------------------------------------------------------------
PLANETS = {
    "Earth": {"g": 9.81, "atm_rho": 1.225, "label": "지구"},
    "Mars":  {"g": 3.71, "atm_rho": 0.020, "label": "화성"},
    "Moon":  {"g": 1.62, "atm_rho": 0.000, "label": "달"},
}
CURRENT_PLANET = "Earth"


def set_planet(name: str):
    global CURRENT_PLANET
    if name not in PLANETS:
        raise ValueError(f"알 수 없는 행성: {name} (Earth / Mars / Moon 중 선택)")
    CURRENT_PLANET = name
    p = PLANETS[name]
    print(f"[ground] 행성 = {p['label']} (g={p['g']} m/s^2)")


def planet_g():
    return PLANETS[CURRENT_PLANET]["g"]


# ----------------------------------------------------------------------
# 2. 지반 (스프링 기초) — 라멘 밑단(z=0) 노드 개수/좌표에 맞춰 자동 생성
# ----------------------------------------------------------------------
SOIL_MAT_TAG_BASE = 1
INTERFACE_NODES = {}

DEFAULT_SPRING_STIFFNESS = {
    1: 5.0e7, 2: 5.0e7, 3: 8.0e7,   # kx, ky, kz
    4: 1.0e6, 5: 1.0e6, 6: 5.0e5,   # rx, ry, rz
}


def get_base_node_tags(nodes, tol=1e-6):
    """z=0인 노드(라멘 밑단)의 원래 nodeTag 목록을 자동으로 찾는다."""
    return [tag for tag, (x, y, z) in nodes.items() if abs(z) < tol]


def build_foundation(footing_xy, stiffness=None):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    k = stiffness or DEFAULT_SPRING_STIFFNESS
    mat_tags = []
    for dof, kval in k.items():
        mat_tag = SOIL_MAT_TAG_BASE + dof
        ops.uniaxialMaterial("Elastic", mat_tag, kval)
        mat_tags.append(mat_tag)

    INTERFACE_NODES.clear()
    for i, (x, y) in enumerate(footing_xy):
        fixed_tag = 10 * (i + 1)
        iface_tag = 10 * (i + 1) + 1
        ops.node(fixed_tag, x, y, 0.0)
        ops.fix(fixed_tag, 1, 1, 1, 1, 1, 1)
        ops.node(iface_tag, x, y, 0.0)
        ops.element("zeroLength", 1000000 + i, fixed_tag, iface_tag,
                    "-mat", *mat_tags, "-dir", *list(k.keys()))
        INTERFACE_NODES[i] = iface_tag

    print(f"[ground] 기초 {len(footing_xy)}개소 생성 완료 (라멘 밑단 노드 수와 일치)")
    return dict(INTERFACE_NODES)


def connect_structure_base(foundation_index, structure_base_node, dofs=(1, 2, 3, 4, 5, 6)):
    iface_tag = INTERFACE_NODES[foundation_index]
    ops.equalDOF(iface_tag, structure_base_node, *dofs)


# ----------------------------------------------------------------------
# 3. 부재 단면(각형강관, SHS)
# ----------------------------------------------------------------------
MEMBER_T_D_RATIO = 0.04           # 단면 형상비 t/D (기둥·보 공통, 돔과 동일하게 유지)
COLUMN_BEAM_D_RATIO = 1.5

# 자동 최소질량 탐색 설정 (돔/하이퍼볼로이드 스크립트와 동일한 이름/기본값 체계)
MASS_SEARCH_PGA_G = 0.6
MASS_SEARCH_START_TON = 2_0000.0       # 반드시 한계 초과가 발생하는 양의 질량에서 시작
MASS_SEARCH_COARSE_STEP_TON = 1_0000.0 # 최초 SAFE 상한을 찾기 위한 고정 증가폭
MASS_SEARCH_TOL_TON = 10.0            # 질량 탐색 간격(결과를 10 ton 단위로 제한)
MASS_SEARCH_MAX_TON = 200_000.0       # 무한 반복 방지용 탐색 상한 (라멘은 총질량이 커서 여유있게)
MASS_SEARCH_SEED = 42
MASS_SEARCH_DIRECTION = 1             # 1=X, 2=Y, 3=Z


def steel_shs_section(D, t, E=2.0e11, G=None):
    """정사각형 강관(SHS, Square Hollow Section) hollow 사각관임

    각 물성치 계산식 (모서리 라운드는 무시하고 완전한 사각형으로 근사):
        A  = D^2 - (D-2t)^2                         (단면적)
        I  = (D^4 - (D-2t)^4) / 12                  (Iy = Iz, 대칭 단면)
        J  = t * (D-t)^3                             (Bredt 박벽 폐단면 근사, 중심선 기준)
        c  = D / 2                                    (굽힘응력용 최외단거리 = 외변폭의 절반)
    """
    if G is None:
        G = E / (2 * 1.3)  # 포아송비 0.3 가정
    Do = D - 2 * t
    A = D ** 2 - Do ** 2
    I = (D ** 4 - Do ** 4) / 12.0
    J = t * (D - t) ** 3
    c = D / 2.0
    return {"A": A, "E": E, "G": G, "J": J, "Iy": I, "Iz": I, "c": c, "D": D, "t": t}


def is_column(nodes, ni, nj):
    """부재가 기둥(수직부재)인지 판정 — 두 절점의 x,y가 같으면 수직(=기둥)."""
    xi, yi, _ = nodes[ni]
    xj, yj, _ = nodes[nj]
    return abs(xi - xj) < 1e-9 and abs(yi - yj) < 1e-9


def _section_for_ele(nodes, ni, nj, column_section, beam_section):
    return column_section if is_column(nodes, ni, nj) else beam_section


def total_member_length_by_type(nodes, edges):
    """기둥 부재 길이 합과 보 부재 길이 합을 각각 구한다. 돔의 total_member_length()를
    기둥/보로 나눠 집계하도록 확장한 것."""
    column_length = 0.0
    beam_length = 0.0
    for ele_tag, ni, nj in edges:
        L = float(np.linalg.norm(np.asarray(nodes[nj], dtype=float)
                                  - np.asarray(nodes[ni], dtype=float)))
        if is_column(nodes, ni, nj):
            column_length += L
        else:
            beam_length += L
    return column_length, beam_length


def size_column_beam_shs_for_target_mass(nodes, edges, target_mass_ton,
                                          t_over_D=MEMBER_T_D_RATIO,
                                          column_beam_ratio=COLUMN_BEAM_D_RATIO,
                                          density=STEEL_DENSITY):
    """기둥:보 외변폭 비율(column_beam_ratio)을 고정한 채, 총 강재 질량(기둥+보)이
    목표값이 되도록 보의 외변폭(beam_D)을 역산하고 기둥은 그 비율만큼 키운다.

    A(D) = D^2 - (D-2t)^2 = k*D^2  (단, t=t_over_D*D, k=4*t_over_D*(1-t_over_D))
    column_D = ratio * beam_D 이므로 A_column = k * ratio^2 * beam_D^2

    총질량 = density * (beam_len*A_beam + column_len*A_column)
           = density * k * beam_D^2 * (beam_len + ratio^2 * column_len)
    위 식을 beam_D에 대해 풀어서 역산
    """
    if target_mass_ton <= 0.0:
        raise ValueError("target_mass_ton은 0보다 커야 합니다.")
    if not 0.0 < t_over_D < 0.5:
        raise ValueError("t_over_D는 0과 0.5 사이여야 합니다.")

    column_len, beam_len = total_member_length_by_type(nodes, edges)
    if column_len <= 0.0 or beam_len <= 0.0:
        raise ValueError("기둥 또는 보의 총 길이가 0입니다.")

    k = 4.0 * t_over_D * (1.0 - t_over_D)
    target_mass_kg = target_mass_ton * 1000.0
    denom = density * k * (beam_len + (column_beam_ratio ** 2) * column_len)

    beam_D = float(np.sqrt(target_mass_kg / denom))
    beam_t = t_over_D * beam_D
    column_D = column_beam_ratio * beam_D
    column_t = t_over_D * column_D

    beam_A = k * beam_D ** 2
    column_A = k * column_D ** 2
    return {
        "beam_D": beam_D, "beam_t": beam_t, "beam_A": beam_A, "beam_length_sum": beam_len,
        "column_D": column_D, "column_t": column_t, "column_A": column_A,
        "column_length_sum": column_len,
    }


def make_section_for_target_mass(target_mass_ton):
    sizing = size_column_beam_shs_for_target_mass(
        RAHMEN_NODES, RAHMEN_EDGES, target_mass_ton=target_mass_ton,
        t_over_D=MEMBER_T_D_RATIO, column_beam_ratio=COLUMN_BEAM_D_RATIO)

    column_section = steel_shs_section(D=sizing["column_D"], t=sizing["column_t"])
    beam_section = steel_shs_section(D=sizing["beam_D"], t=sizing["beam_t"])

    section_info = {
        "target_mass_ton": float(target_mass_ton),
        "column_D": sizing["column_D"], "column_t": sizing["column_t"],
        "column_A": sizing["column_A"],
        "beam_D": sizing["beam_D"], "beam_t": sizing["beam_t"], "beam_A": sizing["beam_A"],
        "column_length_sum": sizing["column_length_sum"],
        "beam_length_sum": sizing["beam_length_sum"],
        "D": sizing["beam_D"], "t": sizing["beam_t"], "A": sizing["beam_A"],
    }
    return column_section, beam_section, section_info


def compute_tributary_mass(nodes, edges, node_tags, section_by_tag,
                            density=STEEL_DENSITY, extra_mass_per_node=0.0):
    mass = defaultdict(float)
    for ele_tag, ni, nj in edges:
        p1 = np.asarray(nodes[ni], dtype=float)
        p2 = np.asarray(nodes[nj], dtype=float)
        L = np.linalg.norm(p2 - p1)
        A_ele = section_by_tag[ele_tag]["A"]
        m_ele = A_ele * L * density
        mass[node_tags[ni]] += m_ele / 2
        mass[node_tags[nj]] += m_ele / 2
    if extra_mass_per_node:
        for stag in node_tags.values():
            mass[stag] += extra_mass_per_node
    return dict(mass)


# ----------------------------------------------------------------------
# 3 - 유효좌굴길이계수(K) — AISC 정렬차트(alignment chart) 근사식
# ----------------------------------------------------------------------
BASE_G = 10.0  # 기초(스프링 지지) 절점의 G값 — 상대적으로 유연하므로 핀에 가깝게 근사

def compute_alignment_chart_K(nodes, edges, section_by_tag, base_tags,
                               base_G=BASE_G, beam_K=1.0):
    """절점별 G(기둥 I/L 합 / 보 I/L 합)를 구하고, 비횡구속 골조용 근사식으로
    기둥마다 유효좌굴길이계수 K를 계산한다. 반환: {ele_tag: K}"""
    base_set = set(base_tags)
    stiff_c = defaultdict(float)
    stiff_b = defaultdict(float)
    for ele_tag, ni, nj in edges:
        sec = section_by_tag[ele_tag]
        L = float(np.linalg.norm(np.asarray(nodes[ni]) - np.asarray(nodes[nj])))
        if is_column(nodes, ni, nj):
            stiff_c[ni] += sec["Iz"] / L
            stiff_c[nj] += sec["Iz"] / L
        else:
            stiff_b[ni] += sec["Iz"] / L
            stiff_b[nj] += sec["Iz"] / L

    def G_at(node):
        if node in base_set:
            return base_G
        sb = stiff_b.get(node, 0.0)
        sc = stiff_c.get(node, 0.0)
        if sb <= 1e-12:
            return base_G  # 보가 없으면 사실상 핀에 가깝다고 근사
        return sc / sb

    K_by_tag = {}
    for ele_tag, ni, nj in edges:
        if is_column(nodes, ni, nj):
            Ga, Gb = G_at(ni), G_at(nj)
            K_by_tag[ele_tag] = float(
                ((1.6 * Ga * Gb + 4.0 * (Ga + Gb) + 7.5) / (Ga + Gb + 7.5)) ** 0.5)
        else:
            K_by_tag[ele_tag] = beam_K
    return K_by_tag


# ----------------------------------------------------------------------
# 4. 라멘 구조물 생성
# ----------------------------------------------------------------------
def build_rahmen_structure(nodes, edges, column_section=None, beam_section=None,
                            extra_mass_per_node=0.0, node_offset=NODE_OFFSET):
    if column_section is None or beam_section is None:
        raise ValueError(
            "column_section과 beam_section을 모두 지정해야 합니다. "
            "steel_shs_section(D, t)로 만들거나, "
            "make_section_for_target_mass(target_mass_ton)으로 목표질량에서 "
            "역산해서 넘기세요.")
        
    base_tags = get_base_node_tags(nodes)
    footing_xy = [(nodes[tag][0], nodes[tag][1]) for tag in base_tags]
    build_foundation(footing_xy)

    transf_vertical = node_offset + 9001    # 부재축이 global Z와 거의 나란할 때
    transf_general = node_offset + 9002     # 그 외 일반적인 경우
    ops.geomTransf("Linear", transf_vertical, 1, 0, 0)
    ops.geomTransf("Linear", transf_general, 0, 0, 1)

    node_tags = {}
    for tag, (x, y, z) in nodes.items():
        stag = node_offset + tag
        ops.node(stag, x, y, z)
        node_tags[tag] = stag

    base_node_tags = {tag: node_tags[tag] for tag in base_tags}
    for i, tag in enumerate(base_tags):
        connect_structure_base(i, node_tags[tag])

    section_by_tag = {}
    elements = {}
    for ele_tag, ni, nj in edges:
        sec = _section_for_ele(nodes, ni, nj, column_section, beam_section)
        section_by_tag[ele_tag] = sec
        A, E, G, J, Iy, Iz = sec["A"], sec["E"], sec["G"], sec["J"], sec["Iy"], sec["Iz"]

        transf_tag = transf_vertical if is_column(nodes, ni, nj) else transf_general
        ops.element("elasticBeamColumn", ele_tag, node_tags[ni], node_tags[nj],
                    A, E, G, J, Iy, Iz, transf_tag)
        elements[ele_tag] = ele_tag

    mass_dict = compute_tributary_mass(nodes, edges, node_tags, section_by_tag,
                                        extra_mass_per_node=extra_mass_per_node)
    total_mass = sum(mass_dict.values())

    K_by_tag = compute_alignment_chart_K(nodes, edges, section_by_tag, base_tags)

    n_col = sum(1 for t, ni, nj in edges if is_column(nodes, ni, nj))
    n_beam = len(edges) - n_col
    print(f"[rahmen] 절점 {len(node_tags)}개 (기초연결 {len(base_node_tags)}곳), "
          f"부재 {len(elements)}개(기둥 {n_col} + 보 {n_beam})로 라멘 구조물 생성 완료 "
          f"(구조 총질량 ≈ {total_mass/1000:.1f} ton, 부재 자중 기반 계산)")

    return {
        "node_tags": node_tags,
        "base_node_tags": base_node_tags,
        "elements": elements,
        "mass_dict": mass_dict,
        "section_by_tag": section_by_tag,
        "K_by_tag": K_by_tag,
        "column_section": column_section,
        "beam_section": beam_section,
    }


# ----------------------------------------------------------------------
# 5. 환경 하중 / 해석 실행
# ----------------------------------------------------------------------

# [1단계] 중력을 가해서 구조물에 질량과 하중을 부여
def apply_gravity(node_mass: dict, ts_tag=1, pattern_tag=1):
    g = planet_g()
    ops.timeSeries("Linear", ts_tag)
    ops.pattern("Plain", pattern_tag, ts_tag)
    for node, mass in node_mass.items():
        ops.mass(node, mass, mass, mass, 0.0, 0.0, 0.0)
        ops.load(node, 0.0, 0.0, -mass * g, 0.0, 0.0, 0.0)
    print(f"[ground] 중력 하중 적용 (g={g}) -> {len(node_mass)}개 절점")


# [2단계] 중력에 의해 건물이 안정적으로 자리잡도록 정적 해석
def run_static(steps=10, reset_time=True):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-6, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    ok = ops.analyze(steps)
    if reset_time:
        ops.loadConst("-time", 0.0)
    print(f"[ground] 정적 해석 {'성공' if ok == 0 else '실패'}")
    return ok


# [3단계] 무게가 실린 건물의 고유 진동수(주파수) 찾기
def eigen(num_modes=3):
    lam = ops.eigen(num_modes)
    freqs = [(l ** 0.5) / (2 * 3.14159265) for l in lam]
    print(f"[ground] 고유진동수(Hz): {[round(f, 3) for f in freqs]}")
    return freqs


# [4단계] 진동수를 찾았으니, 이에 맞춰 에너지를 흡수할 감쇠 장치(브레이크) 세팅
def setup_rayleigh_damping(xi=0.05, mode_i=1, mode_j=3):
    lam = ops.eigen(mode_j)
    wi = lam[mode_i - 1] ** 0.5
    wj = lam[mode_j - 1] ** 0.5
    a0 = 2 * xi * wi * wj / (wi + wj)
    a1 = 2 * xi / (wi + wj)
    ops.rayleigh(a0, 0.0, 0.0, a1)
    print(f"[seismic] Rayleigh 감쇠 적용: xi={xi} (a0={a0:.6f}, a1={a1:.6f})")
    return a0, a1


# [5단계] 가상의 지진파 데이터 생성 — 행성별 지반 주파수 특성(omega_g, xi_g)을 반영
def generate_earthquake(dt, t_total, pga_g=0.5, seed=42,
                         omega_g=None, xi_g=None, t1=2.0, t2=12.0,
                         pga_reference="earth"):
    """가상의 지진파 데이터 생성 (Kanai-Tajimi 필터 기반).
    omega_g/xi_g를 지정하지 않으면 CURRENT_PLANET에 맞는 기본값을 사용한다
    (돔 스크립트와 동일한 행성별 지반 필터 특성)."""

    planet = CURRENT_PLANET

    if omega_g is None or xi_g is None:
        if planet == "Earth":
            default_omega_g, default_xi_g = 15.0, 0.60
        elif planet == "Mars":
            default_omega_g, default_xi_g = 8.5, 0.405
        elif planet == "Moon":
            default_omega_g, default_xi_g = 5.0, 0.25
        else:
            raise ValueError(f"알 수 없는 행성: {planet}")
        if omega_g is None:
            omega_g = default_omega_g
        if xi_g is None:
            xi_g = default_xi_g

    rng = np.random.default_rng(seed)
    n_steps = int(round(t_total / dt)) + 1
    t = np.arange(n_steps) * dt

    white_noise = rng.standard_normal(n_steps)
    freqs = np.fft.fftfreq(n_steps, d=dt) * 2 * np.pi

    H = np.zeros(n_steps, dtype=complex)
    for i, w in enumerate(freqs):
        num = omega_g ** 2 + 2j * xi_g * omega_g * w
        den = omega_g ** 2 - w ** 2 + 2j * xi_g * omega_g * w
        H[i] = num / den if abs(den) > 1e-10 else 0.0

    X = np.fft.fft(white_noise)
    acc = np.real(np.fft.ifft(X * H))

    envelope = np.zeros(n_steps)
    for i, ti in enumerate(t):
        if ti < t1:
            envelope[i] = (ti / t1) ** 2
        elif ti < t2:
            envelope[i] = 1.0
        else:
            envelope[i] = np.exp(-0.3 * (ti - t2))
    acc = acc * envelope

    if pga_reference == "earth":
        target_accel = pga_g * 9.81
    elif pga_reference == "local":
        target_accel = pga_g * planet_g()
    else:
        raise ValueError(f"알 수 없는 pga_reference: {pga_reference} ('earth' 또는 'local')")

    pga_now = np.max(np.abs(acc))
    if pga_now > 0:
        acc = acc * target_accel / pga_now

    print(f"[seismic] 합성 지진파 생성: PGA={pga_g}g ({pga_reference} 기준, "
          f"target={target_accel:.3f} m/s^2), planet={planet} (omega_g={omega_g}, xi_g={xi_g}), "
          f"{t_total}s ({n_steps}step)")
    return t, acc


# [6단계] 생성한 지진파 데이터를 구조물 바닥에 입력(장착)
def apply_earthquake(acc, dt, direction=1, ts_tag=2, pattern_tag=2):
    """지진파를 파일로 저장하지 않고 배열(acc)을 그대로 OpenSees에 전달한다."""
    ops.timeSeries("Path", ts_tag, "-dt", dt, "-values", *acc, "-factor", 1.0)
    ops.pattern("UniformExcitation", pattern_tag, direction, "-accel", ts_tag)
    print(f"[ground] 지진 하중 등록 (파일 미사용, 배열 {len(acc)}개 값 직접 전달, dt={dt})")


# [7단계] 모든 준비를 마쳤으니, 시간에 따라 지진 시뮬레이션 가동! (동적 해석)
def run_dynamic_robust(num_steps, dt, primary="KrylovNewton",
                        fallback="ModifiedNewton", progress_every=200):
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
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
            print(f"[seismic] 진행률 {step/num_steps*100:.0f}% (t={ops.getTime():.2f}s)")

    status = "성공" if fail_count == 0 else f"{fail_count}스텝 수렴실패"
    print(f"[seismic] 동적 해석 완료: {status}")
    return fail_count == 0


# ----------------------------------------------------------------------
# 6. 부재 응력 recorder / 집계 — 연구 목적(응력값 자체, elasticBeamColumn 전 부재 공통)
# ----------------------------------------------------------------------
def setup_member_force_recorders(elements, out_dir="results_members"):
    """반드시 apply_earthquake() 전, run_static() 이후에 호출할 것 (시간이력 누락 방지)."""
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for ele_tag in elements.values():
        path = os.path.join(out_dir, f"force_ele{ele_tag}.txt")
        ops.recorder('Element', '-file', path, '-time', '-ele', ele_tag, 'localForce')
        paths[ele_tag] = path
    print(f"[member] 부재 {len(elements)}개의 부재력(localForce) recorder 등록 완료")
    return paths


def get_member_peak_axial(force_paths):
    result = {}
    for ele_tag, path in force_paths.items():
        data = np.loadtxt(path)
        if data.ndim == 1:
            data = data.reshape(1, -1)

        Ni = -data[:, 1]
        Nj = data[:, 7]

        combined = np.concatenate([Ni, Nj])

        # 가장 음수인 값 = 최대 압축력
        result[ele_tag] = float(np.min(combined))

    return result


def get_member_peak_stress(force_paths, section_by_tag):
    result = {}
    for ele_tag, path in force_paths.items():
        sec = section_by_tag[ele_tag]
        A, Iy, Iz, c = sec["A"], sec["Iy"], sec["Iz"], sec["c"]

        data = np.loadtxt(path)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        Ni, My_i, Mz_i = -data[:, 1], data[:, 5], data[:, 6]
        Nj, My_j, Mz_j = data[:, 7], data[:, 11], data[:, 12]
        """time | FxI FyI FzI MxI MyI MzI | FxJ FyJ FzJ MxJ MyJ MzJ"""

        N = np.concatenate([Ni, Nj])
        My = np.concatenate([My_i, My_j])
        Mz = np.concatenate([Mz_i, Mz_j])

        bending = np.abs(My) * c / Iy + np.abs(Mz) * c / Iz
        axial = N / A
        sigma_plus = axial + bending    # 한쪽 극단섬유
        sigma_minus = axial - bending   # 반대쪽 극단섬유

        pick_plus = np.abs(sigma_plus) >= np.abs(sigma_minus)
        sigma_extreme = np.where(pick_plus, sigma_plus, sigma_minus)

        idx = int(np.argmax(np.abs(sigma_extreme)))
        result[ele_tag] = float(sigma_extreme[idx])
    return result


def get_member_lengths(nodes, edges):
    """RAHMEN_NODES/RAHMEN_EDGES 원본 좌표로 부재별 실제 길이(m)를 계산. 좌굴 계산의 L에 사용."""
    lengths = {}
    for ele_tag, ni, nj in edges:
        p1 = np.asarray(nodes[ni], dtype=float)
        p2 = np.asarray(nodes[nj], dtype=float)
        lengths[ele_tag] = float(np.linalg.norm(p2 - p1))
    return lengths


STEEL_FY = 235.0e6  # Pa


def assess_member_failure(force_paths, lengths, section_by_tag, K_by_tag, fy=STEEL_FY):
    result = {}
    for ele_tag, path in force_paths.items():
        sec = section_by_tag[ele_tag]
        A, E, Iy, Iz, c = sec["A"], sec["E"], sec["Iy"], sec["Iz"], sec["c"]
        I_min = min(Iy, Iz)
        r = (I_min / A) ** 0.5
        K = K_by_tag[ele_tag]

        tensile_capacity = fy * A
        moment_capacity_y = fy * Iy / c
        moment_capacity_z = fy * Iz / c

        data = np.loadtxt(path)
        if data.ndim == 1:
            data = data.reshape(1, -1)

        # 부호 통일: N>0 인장, N<0 압축. I/J단을 같은 배열에 이어 붙인다.
        N = np.concatenate([-data[:, 1], data[:, 7]])
        My = np.concatenate([data[:, 5], data[:, 11]])
        Mz = np.concatenate([data[:, 6], data[:, 12]])

        L = lengths[ele_tag]
        slenderness = (K * L) / r
        p_cr = (np.pi ** 2 * E * I_min) / ((K * L) ** 2)
        compression_capacity = min(tensile_capacity, p_cr)

        axial_capacity = np.where(
            N >= 0.0,
            tensile_capacity,
            compression_capacity,
        )
        dcr_history = (
            np.abs(N) / axial_capacity
            + np.abs(My) / moment_capacity_y
            + np.abs(Mz) / moment_capacity_z
        )

        critical_index = int(np.argmax(dcr_history))
        max_dcr = float(dcr_history[critical_index])
        critical_N = float(N[critical_index])
        critical_My = float(My[critical_index])
        critical_Mz = float(Mz[critical_index])

        sigma_plus = N / A + np.abs(My) * c / Iy + np.abs(Mz) * c / Iz
        sigma_minus = N / A - np.abs(My) * c / Iy - np.abs(Mz) * c / Iz
        max_abs_stress = float(max(np.max(np.abs(sigma_plus)),
                                   np.max(np.abs(sigma_minus))))

        compression_ratio = np.where(N < 0.0, np.abs(N) / p_cr, 0.0)
        max_compression_ratio = float(np.max(compression_ratio))
        yield_fail = max_abs_stress >= fy
        buckling_fail = max_compression_ratio >= 1.0
        interaction_fail = max_dcr >= 1.0

        result[ele_tag] = {
            "yield_fail": yield_fail,
            "buckling_fail": buckling_fail,
            "interaction_fail": interaction_fail,
            "failed": interaction_fail,
            "max_dcr": max_dcr,
            "critical_N": critical_N,
            "critical_My": critical_My,
            "critical_Mz": critical_Mz,
            "max_abs_stress": max_abs_stress,
            "slenderness": slenderness,
            "p_cr": p_cr,
            "compression_capacity": compression_capacity,
            "max_compression_ratio": max_compression_ratio,
            "K": K,
        }
    return result


# ----------------------------------------------------------------------
# 7. Abaqus 스타일 3D solid mesh 시각화 (면=삼각형 분할 Mesh3d + 검은 wireframe)
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
    ex, ey, ez = _local_axes(p1, p2)
    corners = _section_corners(width, height)
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    verts = []
    for p in (p1, p2):
        for (oy, oz) in corners:
            verts.append(p + oy * ey + oz * ez)
    return np.array(verts)


_BOX_FACES = [
    (0, 1, 2, 3), (4, 5, 6, 7),
    (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
]
_BOX_EDGES = [
    (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
]


def _add_box(verts, value, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z):
    base = len(X)
    for v in verts:
        X.append(float(v[0])); Y.append(float(v[1])); Z.append(float(v[2]))
    for (a, b, c, d) in _BOX_FACES:
        I.append(base + a); J.append(base + b); K.append(base + c)
        intensity.append(value)
        I.append(base + a); J.append(base + c); K.append(base + d)
        intensity.append(value)
    for (a, b) in _BOX_EDGES:
        edge_x.extend([verts[a][0], verts[b][0], None])
        edge_y.extend([verts[a][1], verts[b][1], None])
        edge_z.extend([verts[a][2], verts[b][2], None])


def _add_ground_block(footing_xy, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z,
                       depth=200.0, pad=100.0):
    xs = [p[0] for p in footing_xy]
    ys = [p[1] for p in footing_xy]
    x0, x1 = min(xs) - pad, max(xs) + pad
    y0, y1 = min(ys) - pad, max(ys) + pad
    verts = np.array([
        [x0, y0, 0.0], [x1, y0, 0.0], [x1, y1, 0.0], [x0, y1, 0.0],
        [x0, y0, -depth], [x1, y0, -depth], [x1, y1, -depth], [x0, y1, -depth],
    ])
    _add_box(verts, 0.0, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z)


def plot_rahmen_3d_abaqus(rahmen, nodes, edges, member_stress,
                           colorscale="RdBu", ground_depth=None, ground_pad=None,
                           column_vis_side=None, beam_vis_side=None, save_html_path=None):
    """부재 단면을 실제 두께가 있는 solid로 압출해서 Abaqus viewport처럼 렌더링.
    색상은 판정 등급이 아니라 부재 응력(Pa) 그 자체 (압축=음수/인장=양수, RdBu).

    column_vis_side/beam_vis_side : 시각화용 단면 폭/높이(m). 기둥/보 실제 단면(D)을
               기본값으로 쓰되, 순수하게 눈으로 보기 위한 두께라 실제 응력 계산과는 무관.
    """
    column_vis_side = column_vis_side or rahmen["column_section"]["D"]
    beam_vis_side = beam_vis_side or rahmen["beam_section"]["D"]

    node_coord = {tag: ops.nodeCoord(stag) for tag, stag in rahmen["node_tags"].items()}
    base_tags = list(rahmen["base_node_tags"].keys())
    footing_xy = [(nodes[t][0], nodes[t][1]) for t in base_tags]

    span = max(
        max(x for x, y, z in nodes.values()) - min(x for x, y, z in nodes.values()),
        max(y for x, y, z in nodes.values()) - min(y for x, y, z in nodes.values()),
    )
    ground_depth = ground_depth or span * 0.1
    ground_pad = ground_pad or span * 0.05

    X, Y, Z, I, J, K, intensity = [], [], [], [], [], [], []
    edge_x, edge_y, edge_z = [], [], []

    _add_ground_block(footing_xy, X, Y, Z, I, J, K, intensity, edge_x, edge_y, edge_z,
                       depth=ground_depth, pad=ground_pad)

    for ele_tag, ni, nj in edges:
        p1 = node_coord[ni]
        p2 = node_coord[nj]
        w = h = column_vis_side if is_column(nodes, ni, nj) else beam_vis_side
        verts = element_box_vertices(p1, p2, w, h)
        _add_box(verts, member_stress.get(ele_tag, 0.0), X, Y, Z, I, J, K, intensity,
                  edge_x, edge_y, edge_z)

    cmax = max((abs(v) for v in intensity), default=1.0) or 1.0

    fig = go.Figure()
    fig.add_trace(go.Mesh3d(
        x=X, y=Y, z=Z, i=I, j=J, k=K,
        intensity=intensity, intensitymode="cell",  # [버그 수정 — 돔 스크립트와 동일]
        colorscale=colorscale, cmin=-cmax, cmax=cmax,
        opacity=0.95, showscale=True,
        colorbar=dict(title="Stress (Pa)<br>+인장 / -압축"),
        name="Rahmen solid",
        hovertemplate="x=%{x:.2f} m<br>y=%{y:.2f} m<br>z=%{z:.2f} m<extra></extra>",
    ))
    fig.add_trace(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z, mode="lines",
        line=dict(color="black", width=2), name="Element edges", hoverinfo="skip",
    ))
    fig.update_layout(
        title="OpenSeesPy 라멘(Rahmen) 구조물 — Abaqus-style 3D 응력 시각화 (연구 목적: 부재 응력 분포)",
        template="plotly_white", width=1200, height=900,
        scene=dict(xaxis_title="X (m)", yaxis_title="Y (m)", zaxis_title="Z (m)",
                   aspectmode="data",
                   camera=dict(eye=dict(x=1.8, y=-1.8, z=1.1), center=dict(x=0, y=0, z=0))),
        legend=dict(x=0.02, y=0.98),
    )

    if save_html_path:
        out_path = Path(save_html_path).resolve()
        fig.write_html(str(out_path), auto_open=False)
        print(f"[viz] 저장: {out_path}")

    fig.show()
    return fig


# ----------------------------------------------------------------------
# 8. 실행 파이프라인 — PGA 하나를 넣으면 그 강도에서의 파손 여부까지 반환
# ----------------------------------------------------------------------
def run_case(pga_g, target_mass_ton,
             dt=0.01, t_total=15.0, seed=42, direction=1,
             column_section=None, beam_section=None, extra_mass_per_node=0.0,
             pga_reference="earth"):
    if target_mass_ton <= 0.0:
        raise ValueError("target_mass_ton은 0보다 커야 합니다.")
    if direction not in (1, 2, 3):
        raise ValueError("direction은 1(X), 2(Y), 3(Z) 중 하나여야 합니다.")

    if column_section is None and beam_section is None:
        column_section, beam_section, section_info = make_section_for_target_mass(
            target_mass_ton)
    else:
        _fallback = column_section if column_section is not None else beam_section
        section_info = {
            "target_mass_ton": float(target_mass_ton),
            "D": _fallback.get("D") if isinstance(_fallback, dict) else None,
            "t": _fallback.get("t") if isinstance(_fallback, dict) else None,
            "A": _fallback.get("A") if isinstance(_fallback, dict) else None,
        }

    pga_tag = f"{pga_g:.4f}".rstrip("0").rstrip(".").replace(".", "_")
    mass_tag = f"{target_mass_ton:.3f}".rstrip("0").rstrip(".").replace(".", "_")

    rahmen = build_rahmen_structure(RAHMEN_NODES, RAHMEN_EDGES,
                                     column_section=column_section,
                                     beam_section=beam_section,
                                     extra_mass_per_node=extra_mass_per_node)

    apply_gravity(rahmen["mass_dict"])
    static_ok = run_static()
    if static_ok != 0:
        raise RuntimeError(
            f"[PGA={pga_g}g] 중력 정적해석이 수렴하지 않았습니다. "
            f"파손 판정을 진행하지 않습니다."
        )
    eigen(3)

    ops.wipeAnalysis()  # 정적해석 analysis 객체 정리 후 동적해석으로 전환

    # recorder는 지진하중 적용 '전'에 등록해야 시간이력이 전부 기록됨
    force_paths = setup_member_force_recorders(
        rahmen["elements"],
        out_dir=os.path.join(
            "results_mass_search",
            f"mass_{mass_tag}ton",
            f"pga_{pga_tag}g_seed_{seed}_dir_{direction}",
        ),
    )

    t, acc = generate_earthquake(dt, t_total, pga_g=pga_g, seed=seed,
                                  pga_reference=pga_reference)

    setup_rayleigh_damping(xi=0.05, mode_i=1, mode_j=3)
    apply_earthquake(acc, dt=dt, direction=direction)
    dynamic_ok = run_dynamic_robust(num_steps=len(acc) - 1, dt=dt)

    if not dynamic_ok:
        raise RuntimeError(
            f"[PGA={pga_g}g] 동적 해석이 일부 스텝에서 수렴하지 않았습니다. "
            f"이 상태의 결과로 파손 판정을 하지 않습니다. "
            f"(수치적 문제일 가능성이 높으며 실제 구조 붕괴로 해석하지 말 것)"
        )

    member_stress = get_member_peak_stress(force_paths, rahmen["section_by_tag"])
    member_axial = get_member_peak_axial(force_paths)
    lengths = get_member_lengths(RAHMEN_NODES, RAHMEN_EDGES)

    # 돔/하이퍼볼로이드와 동일한 interaction DCR 판정 (부재별 K/단면만 라멘에 맞춤)
    assessment = assess_member_failure(force_paths, lengths,
                                        rahmen["section_by_tag"], rahmen["K_by_tag"])
    failed_members = [tag_ for tag_, info in assessment.items() if info["failed"]]
    max_dcr = max((info["max_dcr"] for info in assessment.values()), default=0.0)

    steel_mass_kg = sum(
        rahmen["section_by_tag"][ele_tag]["A"] * lengths[ele_tag] * STEEL_DENSITY
        for ele_tag in rahmen["elements"].values()
    )

    return {
        "pga_g": pga_g,
        "target_mass_ton": float(target_mass_ton),
        "actual_steel_mass_ton": steel_mass_kg / 1000.0,
        "seed": seed,
        "direction": direction,
        "section_info": section_info,
        "rahmen": rahmen,
        "member_stress": member_stress,
        "member_axial": member_axial,
        "assessment": assessment,
        "failed_members": failed_members,
        "failed_count": len(failed_members),
        "max_dcr": max_dcr,
        "passed": len(failed_members) == 0 and max_dcr < 1.0,
    }


def find_minimum_safe_mass(pga_g=MASS_SEARCH_PGA_G,
                           start_mass=MASS_SEARCH_START_TON,
                           coarse_step=MASS_SEARCH_COARSE_STEP_TON,
                           tol=MASS_SEARCH_TOL_TON,
                           seed=MASS_SEARCH_SEED,
                           direction=MASS_SEARCH_DIRECTION,
                           max_mass=MASS_SEARCH_MAX_TON,
                           dt=None, t_total=None,
                           pga_reference="earth"):

    if start_mass <= 0.0:
        raise ValueError("start_mass는 0보다 커야 합니다.")
    if coarse_step <= 0.0:
        raise ValueError("coarse_step은 0보다 커야 합니다.")
    if tol <= 0.0:
        raise ValueError("tol은 0보다 커야 합니다.")
    if max_mass <= start_mass:
        raise ValueError("max_mass는 start_mass보다 커야 합니다.")

    run_case_time_kwargs = {}
    if dt is not None:
        run_case_time_kwargs["dt"] = dt
    if t_total is not None:
        run_case_time_kwargs["t_total"] = t_total

    cache = {}
    history = []

    def evaluate(mass_ton, stage):
        """동일 질량의 중복 해석을 피하면서 FAIL/SAFE를 판정한다."""
        mass_ton = round(float(mass_ton), 9)
        if mass_ton in cache:
            return cache[mass_ton]

        print(f"\n[{stage}] 질량 {mass_ton:.3f} ton 해석 시작")
        case = run_case(
            pga_g=pga_g,
            target_mass_ton=mass_ton,
            seed=seed,
            direction=direction,
            pga_reference=pga_reference,
            **run_case_time_kwargs,
        )
        safe = bool(case["passed"])
        record = {
            "mass_ton": mass_ton,
            "safe": safe,
            "failed_count": case["failed_count"],
            "max_dcr": case["max_dcr"],
            "stage": stage,
        }
        history.append(record)
        cache[mass_ton] = (safe, case)
        if safe:
            print(
                f"[결과] 질량={mass_ton:.3f} ton, PGA={pga_g}g "
                f"-> 파손 부재 없음 "
                f"(최대 DCR={record['max_dcr']:.4f})"
            )
        else:
            print(
                f"[결과] 질량={mass_ton:.3f} ton, PGA={pga_g}g "
                f"-> 파손 부재 {record['failed_count']}개: "
                f"{case['failed_members']} "
                f"(최대 DCR={record['max_dcr']:.4f})"
            )
        return cache[mass_ton]

    # 1단계: start_mass부터 고정폭으로 올려 최초 SAFE 상한을 찾는다.
    start_safe, _ = evaluate(start_mass, "브라케팅")
    if start_safe:
        raise ValueError(
            f"시작 질량 {start_mass:.3f} ton이 이미 SAFE입니다. "
            "더 작은 FAIL 질량을 start_mass로 입력하십시오."
        )

    lower_failed = float(start_mass)
    while True:
        candidate_mass = lower_failed + coarse_step
        if candidate_mass > max_mass:
            raise RuntimeError(
                f"탐색 상한 {max_mass:.1f} ton까지 SAFE 질량을 찾지 못했습니다."
            )
        candidate_safe, _ = evaluate(candidate_mass, "브라케팅")
        if candidate_safe:
            upper_safe = candidate_mass
            break
        lower_failed = candidate_mass

    print("\n========== 브라케팅 완료 ==========")
    print(f"FAIL 하한: {lower_failed:.3f} ton")
    print(f"SAFE 상한: {upper_safe:.3f} ton")

    # 2단계: tol(기본 10 ton) 단위의 정수 격자에서 이분탐색한다.
    if not math.isclose(start_mass / tol, round(start_mass / tol)):
        raise ValueError(f"start_mass는 {tol:.1f} ton 단위여야 합니다.")
    if not math.isclose(coarse_step / tol, round(coarse_step / tol)):
        raise ValueError(f"coarse_step은 {tol:.1f} ton 단위여야 합니다.")

    low_index = int(round(lower_failed / tol))
    high_index = int(round(upper_safe / tol))
    initial_index_width = high_index - low_index
    expected_iterations = max(0, math.ceil(math.log2(initial_index_width)))
    iteration = 0

    while high_index - low_index > 1:
        iteration += 1
        middle_index = (low_index + high_index) // 2
        middle_mass = middle_index * tol
        middle_safe, _ = evaluate(middle_mass, f"이분법 {iteration}")
        if middle_safe:
            high_index = middle_index
        else:
            low_index = middle_index

        lower_failed = low_index * tol
        upper_safe = high_index * tol
        print(
            f"[범위 갱신] {lower_failed:.3f} ton FAIL "
            f"< M_min <= {upper_safe:.3f} ton SAFE"
        )

    lower_failed = low_index * tol
    upper_safe = high_index * tol

    final_safe, final_case = evaluate(upper_safe, "최종 확인")
    lower_safe, _ = evaluate(lower_failed, "최종 확인")
    if not final_safe or lower_safe:
        raise RuntimeError("최종 FAIL-SAFE 경계 검증이 일관되지 않습니다.")

    result = {
        "min_safe_mass_ton": float(upper_safe),
        "max_failed_mass_ton": float(lower_failed),
        "interval_ton": float(upper_safe - lower_failed),
        "tolerance_ton": float(tol),
        "iterations": iteration,
        "expected_iterations": expected_iterations,
        "history": history,
        "final_case": final_case,
    }

    print("\n========== 최소 SAFE 질량 탐색 완료 ==========")
    print(f"최대 FAIL 질량 : {lower_failed:.3f} ton")
    print(f"최소 SAFE 질량 : {upper_safe:.3f} ton")
    print(f"최종 구간폭    : {upper_safe-lower_failed:.3f} ton")
    print(f"이분법 횟수    : {iteration}회")
    print(f"최대 DCR       : {final_case['max_dcr']:.4f}")
    return result


if __name__ == "__main__":
    set_planet(CURRENT_PLANET)
    print(f"[volume] 라멘 바운딩박스 체적 = {RAHMEN_BOUNDING_VOLUME_M3:.1f} m^3 ")

    # ------------------------------------------------------------------
    # 최소 안전질량 자동 탐색
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print(f"[search] 최소 안전질량 자동 탐색 시작 (PGA={MASS_SEARCH_PGA_G}g, "
          f"{MASS_SEARCH_START_TON:.1f}ton부터, {MASS_SEARCH_TOL_TON:.1f}ton 단위로 확정)")
    print("=" * 60)
    search_result = find_minimum_safe_mass(
        pga_g=MASS_SEARCH_PGA_G,
        start_mass=MASS_SEARCH_START_TON,
        coarse_step=MASS_SEARCH_COARSE_STEP_TON,
        tol=MASS_SEARCH_TOL_TON,
        seed=MASS_SEARCH_SEED,
        direction=MASS_SEARCH_DIRECTION,
        max_mass=MASS_SEARCH_MAX_TON,
    )

    minimum_mass = search_result["min_safe_mass_ton"]
    maximum_failed_mass = search_result["max_failed_mass_ton"]
    final_case = search_result["final_case"]
    section_info = final_case["section_info"]

    q = minimum_mass * 1000.0 / RAHMEN_BOUNDING_VOLUME_M3
    print("\n========== 최종 연구 지표 ==========")
    print(
        f"최종 경계: {maximum_failed_mass:.3f} ton FAIL "
        f"< M_min <= {minimum_mass:.3f} ton SAFE"
    )
    print(
        f"[section] 기둥 D={section_info['column_D']:.6f} m, t={section_info['column_t']:.6f} m, "
        f"A={section_info['column_A']:.6f} m^2 / "
        f"보 D={section_info['beam_D']:.6f} m, t={section_info['beam_t']:.6f} m, "
        f"A={section_info['beam_A']:.6f} m^2 (비율 {COLUMN_BEAM_D_RATIO:.2f}:1 고정)"
    )
    print(f"[efficiency] 단위체적당 강재량 q = M/V = {q:.4f} kg/m^3 "
          f"(V=바운딩박스)")

    # 최소 SAFE 질량 케이스의 부재 응력 분포를 시각화한다.
    plot_rahmen_3d_abaqus(
        final_case["rahmen"],
        RAHMEN_NODES,
        RAHMEN_EDGES,
        final_case["member_stress"],
    )

    print("[rahmen] 브라케팅-이분법 기반 최소 SAFE 질량 탐색 완료.")
