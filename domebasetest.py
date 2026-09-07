import math
import os
import webbrowser
from pathlib import Path
from collections import defaultdict

import numpy as np
import openseespy.opensees as ops
import plotly.graph_objects as go


NODE_OFFSET = 1000  # 구조물 절점 태그 = 원래 nodeTag + NODE_OFFSET, 태그 중복 방지
STEEL_DENSITY = 7850.0  # kg/m^3, 일반구조용 압연강재


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

dome_SCALE = 1 / 20  # 4000m급 돔 -> 200m급 돔으로 (반지름 2000m -> 100m)
dome_NODES = {tag: (x * dome_SCALE, y * dome_SCALE, z * dome_SCALE)
              for tag, (x, y, z) in dome_NODES.items()}

DOME_EDGES = [
    (1, 1, 2),
    (2, 3, 2),
    (3, 4, 3),
    (4, 4, 5),
    (5, 1, 5),
    (6, 3, 1),
    (7, 6, 3),
    (8, 6, 7),
    (9, 8, 7),
    (10, 8, 2),
    (11, 9, 2),
    (12, 9, 10),
    (13, 10, 11),
    (14, 12, 11),
    (15, 12, 13),
    (16, 1, 13),
    (17, 4, 1),
    (18, 14, 4),
    (19, 14, 15),
    (20, 16, 15),
    (21, 16, 17),
    (22, 4, 17),
    (23, 6, 4),
    (24, 3, 7),
    (25, 8, 3),
    (26, 9, 8),
    (27, 5, 14),
    (28, 5, 13),
    (29, 18, 13),
    (30, 18, 12),
    (31, 1, 18),
    (32, 6, 19),
    (33, 6, 20),
    (34, 20, 16),
    (35, 7, 19),
    (36, 8, 21),
    (37, 22, 8),
    (38, 22, 23),
    (39, 23, 21),
    (40, 10, 21),
    (41, 10, 24),
    (42, 21, 24),
    (43, 21, 9),
    (44, 9, 11),
    (45, 11, 18),
    (46, 18, 9),
    (47, 2, 18),
    (48, 25, 12),
    (49, 13, 25),
    (50, 13, 26),
    (51, 27, 15),
    (52, 28, 27),
    (53, 28, 29),
    (54, 29, 30),
    (55, 26, 30),
    (56, 16, 27),
    (57, 16, 31),
    (58, 32, 31),
    (59, 33, 32),
    (60, 19, 33),
    (61, 15, 17),
    (62, 29, 15),
    (63, 27, 29),
    (64, 27, 34),
    (65, 35, 34),
    (66, 36, 35),
    (67, 37, 36),
    (68, 32, 37),
    (69, 7, 22),
    (70, 38, 7),
    (71, 39, 38),
    (72, 40, 39),
    (73, 37, 40),
    (74, 14, 29),
    (75, 14, 26),
    (76, 25, 26),
    (77, 26, 5),
    (78, 26, 29),
    (79, 6, 17),
    (80, 20, 32),
    (81, 17, 20),
    (82, 20, 19),
    (83, 21, 22),
    (84, 41, 23),
    (85, 41, 40),
    (86, 40, 42),
    (87, 36, 42),
    (88, 23, 24),
    (89, 23, 39),
    (90, 33, 39),
    (91, 30, 25),
    (92, 43, 27),
    (93, 43, 35),
    (94, 35, 42),
    (95, 42, 37),
    (96, 33, 37),
    (97, 34, 28),
    (98, 30, 28),
    (99, 16, 43),
    (100, 20, 31),
    (101, 38, 33),
    (102, 22, 38),
    (103, 34, 43),
    (104, 31, 36),
    (105, 43, 31),
    (106, 39, 37),
    (107, 38, 23),
    (108, 39, 41),
    (109, 17, 14),
    (110, 19, 32),
    (111, 19, 38),
    (112, 41, 24),
    (113, 43, 36),
    (114, 36, 32),
]
"""(element번호,노드1,노드2):노드1과 노드2를 이은 부재번호"""


# ----------------------------------------------------------------------
# 0. 내부 체적 계산 (단위체적당 강재량 q=M/V 비교용)
# ----------------------------------------------------------------------
def _find_triangles_from_edges(nodes, edges):
    """엣지 목록에서 삼각형(a,b,c 세 노드가 서로 다 연결된 경우)을 찾는다."""
    adj = defaultdict(set)
    # 1. 모든 노드의 연결 상태(Adjacency List)를 만들기
    for _, i, j in edges:
        adj[i].add(j) # i 노드에 j가 연결
        adj[j].add(i) # j 노드에 i가 연결
        
    triangles = set()
    # 2. 다시 모든 엣지(i, j)를 순회합니다.
    for _, i, j in edges:
        # 3. i와 j에 공통으로 연결된 노드(k)를 찾습니다. (교집합)
        common = adj[i] & adj[j]
        for k in common:
            # 4. i-j, j-k, k-i가 모두 연결되어 있으므로 이는 하나의 삼각형 패널입니다.
            # 중복 방지를 위해 노드 번호를 정렬(sorted)해서 튜플로 저장합니다.
            triangles.add(tuple(sorted((i, j, k))))
    return list(triangles)


def _orient_triangles_outward(nodes, triangles, center):
    """각 삼각형의 법선이 center에서 바깥을 향하도록 정점 순서를 보정한다."""
    center = np.asarray(center, dtype=float) # 돔 내부의 기준 중심점
    fixed = []
    
    for a, b, c in triangles:
        # 1. 삼각형을 이루는 세 꼭짓점의 3D 좌표(벡터)를 가져옴.
        v0 = np.asarray(nodes[a], dtype=float)
        v1 = np.asarray(nodes[b], dtype=float)
        v2 = np.asarray(nodes[c], dtype=float)
        
        # 2. 두 벡터(v1-v0, v2-v0)를 외적(Cross Product)하여 표면에 수직인 법선 벡터를 구함.
        normal = np.cross(v1 - v0, v2 - v0)
        
        # 3. 삼각형의 무게중심(Centroid)을 구합니다.
        centroid = (v0 + v1 + v2) / 3.0
        
        # 4. 중심점에서 표면을 향하는 벡터(centroid - center)와 법선 벡터를 내적(Dot Product)
        # 내적값이 0보다 작으면(음수면) 법선이 돔 안쪽을 향하고 있다는 뜻
        if np.dot(normal, centroid - center) < 0:
            # 5. 안쪽을 향한다면, 정점의 순서를 바꿔(a, c, b) 외적 방향(오른손 법칙)을 뒤집음.
            fixed.append((a, c, b))
        else:
            # 6. 이미 바깥을 향한다면 그대로 유지.
            fixed.append((a, b, c))
    return fixed


def _mesh_signed_volume(nodes, triangles):
    """발산정리: V = (1/6) * sum(v0 . (v1 x v2))"""
    total = 0.0
    for a, b, c in triangles:
        # 1. 세 정점의 좌표(원점에서의 위치 벡터)를 가져옴.
        v0 = np.asarray(nodes[a], dtype=float)
        v1 = np.asarray(nodes[b], dtype=float)
        v2 = np.asarray(nodes[c], dtype=float)
        
        # 2. 스칼라 삼중적 [v0, v1, v2] 계산.
        # v1과 v2를 외적하여 평행사변형 면적 벡터를 구하고, 이를 v0와 내적하여 평행육면체의 체적을 구함.
        total += np.dot(v0, np.cross(v1, v2))
        
    # 3. 평행육면체 체적을 6으로 나누면 사면체의 부피. 이를 모두 합산하여 반환.
    return total / 6.0


def compute_enclosed_volume(nodes, edges):
    """돔 쉘(삼각망)이 감싸는 내부 체적(m^3)을 계산한다.
    1) DOME_EDGES에서 삼각형(쉘 패널)을 찾는다.
    2) z=0 밑단 링 노드를 각도순 정렬 후, 그 중심점을 가상 노드로 추가해서
       부채꼴로 삼각형화한다 (쉘만으로는 밑이 뚫린 사발 모양이라 반드시 닫아야 함).
    3) 전체 삼각형(쉘+바닥캡)의 법선을 바깥쪽으로 통일한 뒤 발산정리로 체적 계산.
    """
    shell_tris = _find_triangles_from_edges(nodes, edges)

    base_tags = [t for t, (x, y, z) in nodes.items() if abs(z) < 1e-6]
    base_tags.sort(key=lambda t: np.arctan2(nodes[t][1], nodes[t][0]))

    CENTER_TAG = -1  # dome_NODES/DOME_EDGES의 실제 태그(1~114)와 겹치지 않는 가상 노드
    cx = float(np.mean([nodes[t][0] for t in base_tags]))
    cy = float(np.mean([nodes[t][1] for t in base_tags]))
    nodes_ext = dict(nodes)
    nodes_ext[CENTER_TAG] = (cx, cy, 0.0)

    n = len(base_tags)
    cap_tris = [(CENTER_TAG, base_tags[i], base_tags[(i + 1) % n]) for i in range(n)]

    all_tris = shell_tris + cap_tris

    zs = [z for (x, y, z) in nodes.values()]
    ref_center = (cx, cy, (min(zs) + max(zs)) / 2.0)  # 돔 내부 대략 중심 (법선 방향 판별용)
    oriented = _orient_triangles_outward(nodes_ext, all_tris, ref_center)

    return abs(_mesh_signed_volume(nodes_ext, oriented))


DOME_ENCLOSED_VOLUME_M3 = compute_enclosed_volume(dome_NODES, DOME_EDGES)


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
# 2. 지반 (스프링 기초) — 돔 밑단 링(z=0) 노드 개수/좌표에 맞춰 자동 생성
# ----------------------------------------------------------------------
SOIL_MAT_TAG_BASE = 1
INTERFACE_NODES = {}

DEFAULT_SPRING_STIFFNESS = {
    1: 5.0e7, 2: 5.0e7, 3: 8.0e7,   # kx, ky, kz
    4: 1.0e6, 5: 1.0e6, 6: 5.0e5,   # rx, ry, rz
}


def get_base_node_tags(nodes, tol=1e-6):
    return [tag for tag, (x, y, z) in nodes.items() if abs(z) < tol]
"""z=0인 노드(돔 밑단 링)의 원래 nodeTag 목록을 자동으로 찾는다."""


def build_foundation(footing_xy, stiffness=None):
    """stiffness=none > 따로 언급이 없으면 default stiffness를 쓸 것."""
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

        # 1. 완전 고정점 생성
        ops.node(fixed_tag, x, y, 0.0)
        ops.fix(fixed_tag, 1, 1, 1, 1, 1, 1)

        # 2. 구조물과 연결될 인터페이스 노드 생성
        ops.node(iface_tag, x, y, 0.0)

        # 3. 두 노드 사이에 zeroLength 스프링 요소 연결
        ops.element("zeroLength", 1000000 + i, fixed_tag, iface_tag,
                    "-mat", *mat_tags, "-dir", *list(k.keys()))
        INTERFACE_NODES[i] = iface_tag

    print(f"[ground] 기초 {len(footing_xy)}개소 생성 완료 (돔 밑단 링 노드 수와 일치)")
    return dict(INTERFACE_NODES)


def connect_structure_base(foundation_index, structure_base_node, dofs=(1, 2, 3, 4, 5, 6)):
    iface_tag = INTERFACE_NODES[foundation_index]
    ops.equalDOF(iface_tag, structure_base_node, *dofs)


# ----------------------------------------------------------------------
# 3. 부재 단면(각형강관, SHS(!!square hollow section!!))
# ----------------------------------------------------------------------
MEMBER_T_D_RATIO = 0.04           # 단면 형상비 t/D (구조별로 동일하게 유지) 두께와 외경의 비
"""mass와 t/D에 맞게 자동으로 단면적이 적용 될 것"""

# 자동 최소질량 탐색 설정
MASS_SEARCH_PGA_G = 0.6
MASS_SEARCH_START_TON = 1_000.0       # 반드시 한계 초과가 발생하는 양의 질량에서 시작
MASS_SEARCH_COARSE_STEP_TON = 1_000.0 # 최초 SAFE 상한을 찾기 위한 고정 증가폭
MASS_SEARCH_TOL_TON = 10.0            # 질량 탐색 간격(결과를 10 ton 단위로 제한)
MASS_SEARCH_MAX_TON = 50_000.0        # 무한 반복 방지용 탐색 상한
MASS_SEARCH_SEED = 42
MASS_SEARCH_DIRECTION = 1             # 1=X, 2=Y, 3=Z

def total_member_length(nodes, edges):
    """모든 부재의 중심선 길이 합계(m)."""
    return sum(
        float(np.linalg.norm(
            np.asarray(nodes[nj], dtype=float) - np.asarray(nodes[ni], dtype=float)
        ))
        for _, ni, nj in edges
    )


def size_uniform_shs_for_target_mass(nodes, edges, target_mass_ton,
                                     t_over_D=MEMBER_T_D_RATIO,
                                     density=STEEL_DENSITY):
    """동일 SHS를 쓰는 구조의 총 강재 질량이 목표값이 되도록 D와 t를 계산한다.

    A_target = M_target / (rho * sum(L))
    t = alpha*D일 때 A = 4*alpha*(1-alpha)*D^2 관계를 역산한다.
    이 방식은 A뿐 아니라 D, t, I, J, c가 하나의 실제 SHS 형상으로 함께 바뀐다.
    """
    if target_mass_ton <= 0.0:
        raise ValueError("target_mass_ton은 0보다 커야 합니다.")
    if not 0.0 < t_over_D < 0.5:
        raise ValueError("t_over_D는 0과 0.5 사이여야 합니다.")

    length_sum = total_member_length(nodes, edges)
    if length_sum <= 0.0:
        raise ValueError("총 부재 길이가 0입니다.")

    target_mass_kg = target_mass_ton * 1000.0
    area = target_mass_kg / (density * length_sum)
    D = np.sqrt(area / (4.0 * t_over_D * (1.0 - t_over_D)))
    t = t_over_D * D
    return float(D), float(t), float(area), float(length_sum)

DOME_TOTAL_MEMBER_LENGTH = total_member_length(dome_NODES, DOME_EDGES)

def steel_shs_section(D, t, E=2.0e11, G=None):
    """정사각형 강관(SHS, Square Hollow Section) hollow 사각관임
    각 물성치 계산식 (모서리 라운드는 무시하고 완전한 사각형으로 근사):
        A  = D^2 - (D-2t)^2                         (단면적)
        I  = (D^4 - (D-2t)^4) / 12                  (Iy = Iz, 대칭 단면),(bh^3/12; b=h) 
        J  = t * (D-t)^3                             (Bredt 박벽 폐단면 근사, 중심선 기준)
        c  = D / 2                                    (굽힘응력용 최외단거리 = 외변폭의 절반)
    """
    if G is None:
        G = E / (2 * 1.3)  # 포아송비 0.3 가정 (기존 G=7.7e10과 동일한 가정)
    Do = D - 2 * t  # 내부 개구부 변폭
    A = D ** 2 - Do ** 2
    I = (D ** 4 - Do ** 4) / 12.0
    J = t * (D - t) ** 3  # 박벽 사각 폐단면 "비틀림상수" (Bredt 공식, 중심선 변폭 D-t 사용)
    c = D / 2.0            # 굽힘응력 계산용 최외단거리 (기존 sqrt(A)/2 근사 대체)
    return {"A": A, "E": E, "G": G, "J": J, "Iy": I, "Iz": I, "c": c, "D": D, "t": t}


def make_section_for_target_mass(target_mass_ton):
    """목표 강재질량에 맞는 동일 SHS 단면과 단면 산정정보를 반환한다."""
    D, t, area, length_sum = size_uniform_shs_for_target_mass(
        dome_NODES,
        DOME_EDGES,
        target_mass_ton=target_mass_ton,
        t_over_D=MEMBER_T_D_RATIO,
    )
    section = steel_shs_section(D=D, t=t)
    section_info = {
        "target_mass_ton": float(target_mass_ton),
        "D": D,
        "t": t,
        "A": area,
        "length_sum": length_sum,
    }
    return section, section_info


def compute_tributary_mass(nodes, edges, node_tags, A, density=STEEL_DENSITY,
                            extra_mass_per_node=0.0):
    """부재 자중(단면적 x 길이 x 밀도)을 양 끝 절점에 절반씩 분배해서 절점질량을 구한다."""
    mass = defaultdict(float)
    for ele_tag, ni, nj in edges:
        p1 = np.asarray(nodes[ni], dtype=float) 
        """절점1"""
        p2 = np.asarray(nodes[nj], dtype=float)
        """절점2"""
        L = np.linalg.norm(p2 - p1)
        """두절점사이의 길이 - 특정 ele의 길이"""
        m_ele = A * L * density
        """구한 특정 ele의 길이와 단면적 밀도를 곱해 ele 1개의 총 질량"""
        mass[node_tags[ni]] += m_ele / 2
        mass[node_tags[nj]] += m_ele / 2
        """계산 단순화를 위한 ele 양쪽 node에 질량을 절반씩 부가"""
    
    if extra_mass_per_node:
        for stag in node_tags.values(): 
            mass[stag] += extra_mass_per_node
            """부가 질량이 있으면 더해줄것"""
    return dict(mass)
    

# ----------------------------------------------------------------------
# 4. 돔 구조물 생성
# ----------------------------------------------------------------------
def build_dome_structure(nodes, edges, section=None, extra_mass_per_node=0.0,
                          node_offset=NODE_OFFSET):
    """dome_NODES / DOME_EDGES를 nodes, edges로 함수내부 범용 매개변수로 처리(종속되지 않게 하기 위함!!)
    section은 필수 인자다 — steel_shs_section(D, t)로 만들어서 넘길 것."""
    if section is None:
        raise ValueError(
            "section을 지정해야 합니다. steel_shs_section(D, t)로 원하는 단면을 만들거나, "
            "make_section_for_target_mass(target_mass_ton)으로 목표질량에서 역산해서 넘기세요.")
    if isinstance(section, dict):
        A, E, G, J, Iy, Iz = (section["A"], section["E"], section["G"],
                              section["J"], section["Iy"], section["Iz"])
        c_outer = section.get("c", (A ** 0.5) / 2)
    else:
        A, E, G, J, Iy, Iz = section
        c_outer = (A ** 0.5) / 2  # 단면 상세정보가 없을 때만 쓰는 근사

    # ---- 모델(ops.model)은 build_foundation() 내부에서 ops.wipe()와 함께 생성된다.
    #      따라서 geomTransf/node/element보다 반드시 먼저 호출해야 한다.
    
    base_tags = get_base_node_tags(nodes)
    footing_xy = [(nodes[tag][0], nodes[tag][1]) for tag in base_tags]
    build_foundation(footing_xy)
    """앞서 구한 z=0지점의 테그들을 불러와 build_foundation(강성포함) 만들기"""

    transf_vertical = node_offset + 9001    # 부재축이 global Z와 거의 나란할 때(10001:수직)
    transf_general = node_offset + 9002     # 그 외 일반적인 경우(10002:수직하지않은 이외의 모든)
    """node offset=1000, 부재가 수직한지 아닌지를 10001,10002로 구분"""
    # PDelta: 자중 축력에 의한 2차 모멘트 효과를 동적해석 강성/저항력에 반영한다.
    ops.geomTransf("PDelta", transf_vertical, 1, 0, 0)
    ops.geomTransf("PDelta", transf_general, 0, 0, 1)
    """회전하는 끝면의 단면적의 기준을 정해줌 - geomtransf"""

    node_tags = {}
    for tag, (x, y, z) in nodes.items():
        stag = node_offset + tag
        ops.node(stag, x, y, z)
        node_tags[tag] = stag
        """노드 반영"""

    base_node_tags = {tag: node_tags[tag] for tag in base_tags}
    for i, tag in enumerate(base_tags):
        connect_structure_base(i, node_tags[tag])
        """z=0의 노드를 찾아 지반에 연결"""

    elements = {}
    for ele_tag, ni, nj in edges:
        xi, yi, zi = nodes[ni]
        xj, yj, zj = nodes[nj]
        # z의 값만 다를시 vertical element, 1e-9 = 미세한 차이를 고려한 거의 1
        is_vertical = abs(xi - xj) < 1e-9 and abs(yi - yj) < 1e-9
        transf_tag = transf_vertical if is_vertical else transf_general
        ops.element("elasticBeamColumn", ele_tag, node_tags[ni], node_tags[nj],
                    A, E, G, J, Iy, Iz, transf_tag)
        elements[ele_tag] = ele_tag
        """ele 조립"""

    mass_dict = compute_tributary_mass(nodes, edges, node_tags, A,
                                        extra_mass_per_node=extra_mass_per_node)
    """앞서 정한 질량(extra mass per node)를 각각 elements에 반영"""

    total_mass = sum(mass_dict.values())
    steel_mass = A * total_member_length(nodes, edges) * STEEL_DENSITY
    extra_mass = total_mass - steel_mass

    """순서 = 노드 - 베이스 고정 - element - mass"""

    print(f"[dome] 절점 {len(node_tags)}개 (기초연결 {len(base_node_tags)}곳), "
          f"부재 {len(elements)}개로 돔 구조물 생성 완료 "
          f"(SHS={steel_mass/1000:.2f} ton, 추가질량={extra_mass/1000:.2f} ton, "
          f"동적 총질량={total_mass/1000:.2f} ton)")

    return {
        "node_tags": node_tags,
        "base_node_tags": base_node_tags,
        "elements": elements,
        "mass_dict": mass_dict,
        "section": (A, E, G, J, Iy, Iz),
        "c": c_outer,
        "steel_mass_kg": steel_mass,
        "extra_mass_kg": extra_mass,
        "total_mass_kg": total_mass,
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


# [2단계] 중력에 의해 건물이 안정적으로 자리잡도록 정적 해석 (준비 운동)
def run_static(steps=10, reset_time=True):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-6, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    ok = ops.analyze(steps)
    if reset_time:
        ops.loadConst("-time", 0.0)
    print(f"[ground] 정적 해석 {'성공' if ok == 0 else '실패'}")
    return ok
"""건물에 중력과 하중을 반영한 뒤에 지진을 부여하기전 상태로 준비"""

# [3단계] 무게가 실린 건물의 고유 진동수(주파수) 찾기
def eigen(num_modes=3):
    lam = ops.eigen(num_modes)
    freqs = [(l ** 0.5) / (2 * 3.14159265) for l in lam] #lam^0.5/2Pi
    """eigen을 통해 고유치를 구하고 이를 통해 고유진동수(frequency) 구함 """

    print(f"[ground] 고유진동수(Hz): {[round(f, 3) for f in freqs]}")
    return freqs


# [4단계] 진동수를 찾았으니, 이에 맞춰 에너지를 흡수할 감쇠 장치(브레이크) 세팅
#xi = 0.05는 철골, 콘크리트 실험 설계 기준 입증된 평균 방어력 값 5% , 지진이 자연스럽게 멈추기 위한 레일리 감쇠를 도입(대략 1번,3번 고유치로)
def setup_rayleigh_damping(xi=0.05, mode_i=1, mode_j=3):
    lam = ops.eigen(mode_j)
    wi = lam[mode_i - 1] ** 0.5 #각속도로 변환
    wj = lam[mode_j - 1] ** 0.5
    a0 = 2 * xi * wi * wj / (wi + wj)
    a1 = 2 * xi / (wi + wj)
    ops.rayleigh(a0, 0.0, 0.0, a1)
    print(f"[seismic] Rayleigh 감쇠 적용: xi={xi} (a0={a0:.6f}, a1={a1:.6f})")
    return a0, a1


# [5단계] 가상의 지진파 데이터 생성
def generate_earthquake(dt, t_total, pga_g=0.5, seed=42,
                         omega_g=None, xi_g=None, t1=2.0, t2=12.0,
                         pga_reference="earth"):
    """가상의 지진파 데이터 생성 (Kanai-Tajimi 필터 기반)"""

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

    rng = np.random.default_rng(seed)  # 랜덤 세팅
    n_steps = int(round(t_total / dt)) + 1
    t = np.arange(n_steps) * dt

    white_noise = rng.standard_normal(n_steps)
    """무작위 정규분포 숫자를 연속으로 뽑아내어 백색소음 생성"""
    freqs = np.fft.fftfreq(n_steps, d=dt) * 2 * np.pi
    """np의 주파수 함수(fftfreq)를 통해 얻은 지진파 주파수 눈금에 2파이를 곱해 각주파수를 구함"""

    H = np.zeros(n_steps, dtype=complex)
    for i, w in enumerate(freqs):
        num = omega_g ** 2 + 2j * xi_g * omega_g * w
        den = omega_g ** 2 - w ** 2 + 2j * xi_g * omega_g * w
        H[i] = num / den if abs(den) > 1e-10 else 0.0
    """주파수 응답 함수(Frequency Response Function, H) -> 카나이-타지미(Kanai-Tajimi) 공식"""

    X = np.fft.fft(white_noise)
    acc = np.real(np.fft.ifft(X * H))
    """백색소음의 파장을 지진화 시킴"""

    envelope = np.zeros(n_steps)
    for i, ti in enumerate(t):
        if ti < t1:
            envelope[i] = (ti / t1) ** 2
            """ti(initial)부터 t1까지 서서히 커지는 지진파 생성을 위한 제곱식"""
        elif ti < t2:
            envelope[i] = 1.0
            """t1부터 t2까지 지진파가 가장 강한 순간"""
        else:
            envelope[i] = np.exp(-0.3 * (ti - t2))
            """t2이후에 점차 줄어들어 0에 수렴하는 exp식"""
    acc = acc * envelope
    """t에 따른 지진파 변형을 반영"""

    if pga_reference == "earth":
        target_accel = pga_g * 9.81
    elif pga_reference == "local":
        target_accel = pga_g * planet_g()
        """planet g 와 current planet이랑 이어져"""
    else:
        raise ValueError(f"알 수 없는 pga_reference: {pga_reference} ('earth' 또는 'local')")

    pga_now = np.max(np.abs(acc))  # 진동크기가 가장 큰 순간의 가속도
    if pga_now > 0:
        acc = acc * target_accel / pga_now
        """최대크기로 나누어 최대값을 1로 만들어 진동을 유지하며 최대 1을 넘지 않도록"""

    print(f"[seismic] 합성 지진파 생성: PGA={pga_g}g ({pga_reference} 기준, "
          f"target={target_accel:.3f} m/s^2), {t_total}s ({n_steps}step)")
    return t, acc


# [6단계] 생성한 지진파 데이터를 구조물 바닥에 입력(장착)
def apply_earthquake(acc, dt, direction=1, ts_tag=2, pattern_tag=2):
    ops.timeSeries("Path", ts_tag, "-dt", dt, "-values", *acc, "-factor", 1.0)
    ops.pattern("UniformExcitation", pattern_tag, direction, "-accel", ts_tag)
    print(f"[ground] 지진 하중 등록 (파일 미사용, 배열 {len(acc)}개 값 직접 전달, dt={dt})")


# [7단계] 모든 준비를 마쳤으니, 시간에 따라 지진 시뮬레이션 가동! (동적 해석)
def run_dynamic_robust(num_steps, dt, primary="KrylovNewton",
                        fallback="ModifiedNewton", progress_every=200):
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
            print(f"[seismic] 진행률 {step/num_steps*100:.0f}% (t={ops.getTime():.2f}s)")

    status = "성공" if fail_count == 0 else f"{fail_count}스텝 수렴실패"
    print(f"[seismic] 동적 해석 완료: {status}")
    return fail_count == 0


# ----------------------------------------------------------------------
# 6. 부재 응력 recorder / 집계 — 연구 목적(응력값 자체, elasticBeamColumn 전 부재 공통)
# ----------------------------------------------------------------------
def setup_member_force_recorders(elements, out_dir="results_members"):
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for ele_tag in elements.values():
        path = os.path.join(out_dir, f"force_ele{ele_tag}.txt")
        """부재별 force"""
        ops.recorder('Element', '-file', path, '-time', '-ele', ele_tag, 'localForce')
        paths[ele_tag] = path
    print(f"[member] 부재 {len(elements)}개의 부재력(localForce) recorder 등록 완료")
    return paths


def get_member_peak_axial(force_paths):
    """localForce 시간이력에서 부재별 최대(절대값) 축력 N을 뽑는다.
    좌굴 판정은 굽힘이 섞인 응력이 아니라 순수 축력으로 해야 하므로 별도 함수로 분리.
    OpenSees 3D localForce는 I단/J단 순서로 FxI FyI FzI MxI MyI MzI | FxJ FyJ FzJ MxJ MyJ MzJ
    Fx로 부재파손 확인"""
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


def get_member_peak_stress(force_paths, A, Iy, Iz, cy, cz):
    """localForce로부터 sigma = N/A ± (|My|*cz/Iy + |Mz|*cy/Iz) 계산.
    A, Iy, Iz, cy, cz는 build_dome_structure()에서 쓴 section/c와 일치시킬 것.
    반환: {ele_tag: stress(Pa)} (인장=양수, 압축=음수)"""
    result = {}
    for ele_tag, path in force_paths.items():
        data = np.loadtxt(path)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        Ni, My_i, Mz_i = -data[:, 1], data[:, 5], data[:, 6]
        Nj, My_j, Mz_j = data[:, 7], data[:, 11], data[:, 12]
        """time | FxI FyI FzI MxI MyI MzI | FxJ FyJ FzJ MxJ MyJ MzJ"""

        N = np.concatenate([Ni, Nj])
        My = np.concatenate([My_i, My_j])
        Mz = np.concatenate([Mz_i, Mz_j])

        bending = np.abs(My) * cz / Iy + np.abs(Mz) * cy / Iz
        axial = N / A
        sigma_plus = axial + bending    # 한쪽 극단섬유
        sigma_minus = axial - bending   # 반대쪽 극단섬유

        pick_plus = np.abs(sigma_plus) >= np.abs(sigma_minus)
        sigma_extreme = np.where(pick_plus, sigma_plus, sigma_minus)

        idx = int(np.argmax(np.abs(sigma_extreme)))
        result[ele_tag] = float(sigma_extreme[idx])
    return result


def get_member_lengths(nodes, edges):
    """dome_NODES/DOME_EDGES 원본 좌표로 부재별 실제 길이(m)를 계산. 좌굴 계산의 L에 사용."""
    lengths = {}
    for ele_tag, ni, nj in edges:
        p1 = np.asarray(nodes[ni], dtype=float)
        p2 = np.asarray(nodes[nj], dtype=float)
        lengths[ele_tag] = float(np.linalg.norm(p2 - p1))
    return lengths


STEEL_FY = 235.0e6  # Pa
BUCKLING_K = 1.0


def assess_member_failure(force_paths, lengths, section, c,
                           fy=STEEL_FY, K=BUCKLING_K):
    """각 시간스텝·각 부재 끝단에서 축력-2축휨 조합 DCR을 계산한다.

    비교 연구용 단순 선형 상관식:
        DCR = |N|/Pn + |My|/Mny + |Mz|/Mnz

    인장일 때 Pn=Fy*A, 압축일 때 Pn=min(Fy*A, Euler Pcr)로 둔다.
    서로 다른 시간의 축력과 모멘트 최대값을 섞지 않고 같은 시각의 하중을 사용한다.
    DCR >= 1이면 한계 초과로 판정한다. 실제 파단 과정을 직접 모사하는 식은 아니다.
    """
    A, E, G, J, Iy, Iz = section
    I_min = min(Iy, Iz)
    r = (I_min / A) ** 0.5

    tensile_capacity = fy * A
    moment_capacity_y = fy * Iy / c
    moment_capacity_z = fy * Iz / c

    result = {}
    for ele_tag, path in force_paths.items():
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
        }
    return result

# ----------------------------------------------------------------------
# 7. Abaqus 스타일 3D solid mesh 시각화 (면=삼각형 분할 Mesh3d + 검은 wireframe)
#    임의 방향 부재를 압출하는 방식이라 돔처럼 부재 방향이 제각각이어도 그대로 동작.
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


def plot_dome_3d_abaqus(dome, nodes, edges, member_stress,
                         colorscale="RdBu", ground_depth=None, ground_pad=None,
                         vis_side=None, save_html_path=None):
    if vis_side is None:
        vis_side = 2.0 * dome["c"]  # c = D/2 (steel_shs_section 참고)

    node_coord = {tag: ops.nodeCoord(stag) for tag, stag in dome["node_tags"].items()}
    base_tags = list(dome["base_node_tags"].keys())
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

    w = h = vis_side
    for ele_tag, ni, nj in edges:
        p1 = node_coord[ni]
        p2 = node_coord[nj]
        verts = element_box_vertices(p1, p2, w, h)
        _add_box(verts, member_stress.get(ele_tag, 0.0), X, Y, Z, I, J, K, intensity,
                  edge_x, edge_y, edge_z)

    cmax = max((abs(v) for v in intensity), default=1.0) or 1.0

    fig = go.Figure()
    fig.add_trace(go.Mesh3d(
        x=X, y=Y, z=Z, i=I, j=J, k=K,
        intensity=intensity, intensitymode="cell",  # [버그 수정 5]
        colorscale=colorscale, cmin=-cmax, cmax=cmax,
        opacity=0.95, showscale=True,
        colorbar=dict(title="Stress (Pa)<br>+인장 / -압축"),
        name="Dome solid",
        hovertemplate="x=%{x:.2f} m<br>y=%{y:.2f} m<br>z=%{z:.2f} m<extra></extra>",
    ))

    fig.add_trace(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z, mode="lines",
        line=dict(color="black", width=2), name="Element edges", hoverinfo="skip",
    ))
    fig.update_layout(
        title="OpenSeesPy 지오데식 돔 — Abaqus-style 3D 응력 시각화 (연구 목적: 부재 응력 분포)",
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
             section=None, extra_mass_per_node=0.0, pga_reference="earth"):
    if target_mass_ton <= 0.0:
        raise ValueError("target_mass_ton은 0보다 커야 합니다.")
    if direction not in (1, 2, 3):
        raise ValueError("direction은 1(X), 2(Y), 3(Z) 중 하나여야 합니다.")

    if section is None:
        section, section_info = make_section_for_target_mass(target_mass_ton)
    else:
        section_info = {
            "target_mass_ton": float(target_mass_ton),
            "D": section.get("D") if isinstance(section, dict) else None,
            "t": section.get("t") if isinstance(section, dict) else None,
            "A": section.get("A") if isinstance(section, dict) else section[0],
            "length_sum": total_member_length(dome_NODES, DOME_EDGES),
        }

    pga_tag = f"{pga_g:.4f}".rstrip("0").rstrip(".").replace(".", "_")
    mass_tag = f"{target_mass_ton:.3f}".rstrip("0").rstrip(".").replace(".", "_")

    dome = build_dome_structure(dome_NODES, DOME_EDGES, section=section,
                                 extra_mass_per_node=extra_mass_per_node)
    """앞서 정한 범용 매개변수에 실제 노드 대입(node -> dome_Nodes)"""

    apply_gravity(dome["mass_dict"])
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
        dome["elements"],
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

    A, E, G, J, Iy, Iz = dome["section"]
    c = dome["c"]  # 기존 c=sqrt(A)/2 근사 대신, 각형강관 외변폭 D/2를 그대로 사용
    member_stress = get_member_peak_stress(force_paths, A=A, Iy=Iy, Iz=Iz, cy=c, cz=c)
    member_axial = get_member_peak_axial(force_paths)
    lengths = get_member_lengths(dome_NODES, DOME_EDGES)

    assessment = assess_member_failure(
        force_paths,
        lengths,
        dome["section"],
        c,
    )
    failed_members = [tag_ for tag_, info in assessment.items() if info["failed"]]
    max_dcr = max((info["max_dcr"] for info in assessment.values()), default=0.0)

    return {
        "pga_g": pga_g,
        "target_mass_ton": float(target_mass_ton),
        "actual_steel_mass_ton": dome["steel_mass_kg"] / 1000.0,
        "seed": seed,
        "direction": direction,
        "section_info": section_info,
        "dome": dome,
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
                           dt=0.01, t_total=15.0,
                           pga_reference="earth"):
    """고정폭 브라케팅 후 이분법으로 최소 SAFE 질량을 찾는다.

    1) start_mass는 반드시 FAIL이어야 한다.
    2) coarse_step씩 증가시켜 최초 SAFE 상한을 찾는다.
    3) FAIL 하한과 SAFE 상한의 차가 tol 이하가 될 때까지 이분한다.

    반환되는 min_safe_mass_ton은 최종 구간의 SAFE측 상한이다. 즉,
    max_failed_mass_ton < 실제 임계질량 <= min_safe_mass_ton 관계를 갖는다.
    """
    if start_mass <= 0.0:
        raise ValueError("start_mass는 0보다 커야 합니다.")
    if coarse_step <= 0.0:
        raise ValueError("coarse_step은 0보다 커야 합니다.")
    if tol <= 0.0:
        raise ValueError("tol은 0보다 커야 합니다.")
    if max_mass <= start_mass:
        raise ValueError("max_mass는 start_mass보다 커야 합니다.")

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
            dt=dt,
            t_total=t_total,
            seed=seed,
            direction=direction,
            pga_reference=pga_reference,
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
    # 따라서 모든 시험 질량과 최종 최소 SAFE 질량은 tol의 정수배가 된다.
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

    # 양 경계는 evaluate 결과로 이미 검증됐으며, 캐시에서 최종 SAFE 케이스를 얻는다.
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
    print(f"[volume] 돔 내부 감싸인 체적 = {DOME_ENCLOSED_VOLUME_M3:.1f} m^3")

    print(
        f"[mass search] PGA={MASS_SEARCH_PGA_G}g, "
        f"시작={MASS_SEARCH_START_TON:.1f} ton, "
        f"고정 증가폭={MASS_SEARCH_COARSE_STEP_TON:.1f} ton, "
        f"허용오차={MASS_SEARCH_TOL_TON:.1f} ton"
    )

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

    q = minimum_mass * 1000.0 / DOME_ENCLOSED_VOLUME_M3
    print("\n========== 최종 연구 지표 ==========")
    print(
        f"최종 경계: {maximum_failed_mass:.3f} ton FAIL "
        f"< M_min <= {minimum_mass:.3f} ton SAFE"
    )
    print(
        f"[section] D={section_info['D']:.6f} m, "
        f"t={section_info['t']:.6f} m, A={section_info['A']:.6f} m^2"
    )
    print(f"[efficiency] 단위체적당 강재량 q = M/V = {q:.4f} kg/m^3")

    # 최소 SAFE 질량 케이스의 부재 응력 분포를 시각화한다.
    plot_dome_3d_abaqus(
        final_case["dome"],
        dome_NODES,
        DOME_EDGES,
        final_case["member_stress"],
    )

    print("[dome] 브라케팅-이분법 기반 최소 SAFE 질량 탐색 완료.")
