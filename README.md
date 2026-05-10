import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os

# ─────────────────────────────────────────────────────────────────────
# 0. 전역 파라미터
# ─────────────────────────────────────────────────────────────────────
# 구조물 파라미터
NUM_STORIES   = 3          # 층수
STORY_HEIGHT  = 3.5        # 층고 [m]
BAY_WIDTH     = 6.0        # 스팬 [m]
NUM_BAYS      = 2          # 스팬 수

# RC 단면 (보/기둥)
col_b, col_h  = 0.5, 0.5   # 기둥 단면 [m]
beam_b, beam_h = 0.3, 0.5  # 보 단면 [m]
fc_MPa        = 28.0        # 콘크리트 압축강도 [MPa]
fy_MPa        = 400.0       # 철근 항복강도 [MPa]
Es_GPa        = 200.0       # 철근 탄성계수 [GPa]

# 재료 상수 (SI: N, m 단위)
MPa           = 1e6
GPa           = 1e9
Ec            = 4700 * (fc_MPa**0.5) * MPa   # 콘크리트 탄성계수 [Pa]
Es            = Es_GPa * GPa                  # 철근 탄성계수 [Pa]

# 지반 파라미터 (3개 층)
# [층번호, 두께(m), 단위중량(kN/m³), 전단파속도Vs(m/s), 감쇠비(%), 비선형여부]
SOIL_LAYERS = [
    {"id": 1, "thickness": 3.0,  "gamma": 18.0, "Vs": 150, "xi": 0.05, "nonlinear": True},   # 연약 표층
    {"id": 2, "thickness": 7.0,  "gamma": 19.5, "Vs": 280, "xi": 0.04, "nonlinear": True},   # 중간층
    {"id": 3, "thickness": 15.0, "gamma": 21.0, "Vs": 500, "xi": 0.02, "nonlinear": False},  # 기반암층
]

# 지진 해석 파라미터
dt            = 0.01        # 시간 간격 [s]
t_total       = 20.0        # 해석 총 시간 [s]
SCALE_PGA     = 0.3         # 목표 PGA [g]

# ─────────────────────────────────────────────────────────────────────
# 1. 지진파 생성 (El Centro 1940 모사 - Kanai-Tajimi 필터)
# ─────────────────────────────────────────────────────────────────────
def generate_earthquake(dt, t_total, pga_g=0.3, seed=42):
    """
    Kanai-Tajimi 스펙트럼 기반 인공 지진파 생성
    - ωg=15.0 rad/s, ξg=0.6 (연약지반 특성)
    - Envelope 함수로 비정상 특성 부여
    """
    np.random.seed(seed)
    n_steps = int(t_total / dt)
    t = np.linspace(0, t_total, n_steps)

    # 백색잡음 생성
    white_noise = np.random.randn(n_steps)

    # Kanai-Tajimi 필터 적용 (주파수 영역)
    omega_g = 15.0   # 지반 탁월 각주파수 [rad/s]
    xi_g    = 0.60   # 지반 감쇠비

    from numpy.fft import fft, ifft, fftfreq
    freqs = fftfreq(n_steps, d=dt) * 2 * np.pi
    H_kt  = np.zeros(n_steps, dtype=complex)
    for i, w in enumerate(freqs):
        num = omega_g**2 + 2j * xi_g * omega_g * w
        den = omega_g**2 - w**2 + 2j * xi_g * omega_g * w
        H_kt[i] = num / den if abs(den) > 1e-10 else 0.0

    # 필터링
    X   = fft(white_noise)
    acc = np.real(ifft(X * H_kt))

    # Envelope 함수 (상승-지속-감쇠)
    envelope = np.zeros(n_steps)
    t1, t2 = 2.0, 12.0
    for i, ti in enumerate(t):
        if ti < t1:
            envelope[i] = (ti / t1) ** 2
        elif ti < t2:
            envelope[i] = 1.0
        else:
            envelope[i] = np.exp(-0.3 * (ti - t2))

    acc = acc * envelope

    # PGA 정규화
    pga_current = np.max(np.abs(acc))
    if pga_current > 0:
        acc = acc * (pga_g * 9.81) / pga_current

    return t, acc


# ─────────────────────────────────────────────────────────────────────
# 2. OpenSeesPy 모델 구축
# ─────────────────────────────────────────────────────────────────────
def build_model():
    """
    구조물 + 지반 통합 모델 구축
    
    노드 번호 체계:
    - 1xx : 지반 노드 (101~103: 각 지반층 경계)
    - 1~99: 구조물 노드 (층×(스팬+1) 배열)
    - 기초 노드: 구조물 1층과 지반 상단 공유
    """
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)  # 2D, DOF=3 (Ux, Uy, Rz)

    # ── 재료 정의 ────────────────────────────────────────────────────
    # 콘크리트 (Concrete02: *비선형 압축 + 인장 softening*)
    #uniaxialMaterial('Concrete02', matTag, fpc, epsc0, fpcu, epsU, lamda, ft, Ets)
    fc   = -fc_MPa * MPa #콘크리트 최대 압축 강도(-)
    ec0  = -0.002 #최대 압축 강도 변형률
    fcu  = 0.2 * fc #crushing strength
    ecu  = -0.01 # ultimate strength
    #lamda = 하역 경사와 초기 경사의 비율
    ft   = 0.1 * abs(fc) #tensile strength
    Ets  = 0.1 * Ec #tension softening stiffness
    ops.uniaxialMaterial('Concrete02', 1, fc, ec0, fcu, ecu, 0.1, ft, Ets)

    # 철근 (Steel02: Giuffré-Menegotto-Pinto)
    # uniaxialMaterial('Steel02', matTag, Fy, E0, b, *params, a1=a2*Fy/E0, a2=1.0, a3=a4*Fy/E0, a4=1.0, sigInit=0.0)
    fy  = fy_MPa * MPa
    b   = 0.01   # 변형경화율
    ops.uniaxialMaterial('Steel02', 2, fy, Es, b, 18, 0.925, 0.15)

    # 지반 재료 (각 층 개별 Elastic + Hysteretic 비선형 스프링)
    for i, layer in enumerate(SOIL_LAYERS):
        mat_id = 10 + i + 1
        rho    = layer["gamma"] * 1000 / 9.81   # 밀도 [kg/m³]
        G      = rho * layer["Vs"] ** 2          # 전단탄성계수 [Pa]
        K      = G * 2 * (1 + 0.3) / (3 * (1 - 2 * 0.3))  # 체적탄성계수
        if layer["nonlinear"]:
            # Hysteretic: 비선형 이력 거동 (감쇠+강성저하)
            tau_y  = G * 0.005          # 항복 전단응력
            tau_u  = G * 0.02           # 극한 전단응력
            gam_y  = tau_y / G
            gam_u  = tau_u / G
            ops.uniaxialMaterial('Hysteretic', mat_id,
                                 tau_y, gam_y, tau_u, gam_u,
                                 -tau_y, -gam_y, -tau_u, -gam_u,
                                 0.8, 0.5, 0.0, 0.0, 0.1)
        else:
            ops.uniaxialMaterial('Elastic', mat_id, G)

    # ── 단면 정의 ────────────────────────────────────────────────────
    # 기둥 단면 (Fiber Section) - A, 2모멘트 Iz
    col_A   = col_b * col_h
    col_Iz  = col_b * col_h ** 3 / 12
    beam_A  = beam_b * beam_h
    beam_Iz = beam_b * beam_h ** 3 / 12

    # 기둥 (secTag=1): Fiber 단면
    ops.section('Fiber', 1)
    ops.patch('rect', 1, 10, 10,
              -col_h/2, -col_b/2, col_h/2, col_b/2)
    rebar_A = 0.0025 * col_A / 8   # 8개 철근
    rebar_d = col_h / 2 * 0.85
    for y_sign in [-1, 1]:
        for x_sign in [-1, 0, 1]:
            if x_sign == 0:
                ops.fiber(y_sign * rebar_d, 0, rebar_A * 2, 2)
            else:
                ops.fiber(y_sign * rebar_d, x_sign * rebar_d * 0.6, rebar_A, 2)

    # 보 (secTag=2): 탄성 단면 (보수적)
    ops.section('Elastic', 2, Ec, beam_A, beam_Iz)

    # ── 노드 생성 ────────────────────────────────────────────────────
    # 구조물 노드: 노드번호 = (층번호)*(NUM_BAYS+1) + (기둥번호+1)
    # ex 101- 1층 1번 노드 , 305- 3층 5번 노드
    # 0층(기초): floor=0
    node_id = {}
    for floor in range(NUM_STORIES + 1):
        y = floor * STORY_HEIGHT
        for col in range(NUM_BAYS + 1):
            x   = col * BAY_WIDTH
            nid = floor * 100 + col + 1
            ops.node(nid, x, y)
            node_id[(floor, col)] = nid

    # 지반 노드 (1D 등가선형 컬럼 - 좌측 기준)
    soil_nodes = []
    y_soil = 0.0
    soil_nodes.append(200)
    ops.node(200, 0.0, 0.0)   # 지표면 (구조물 기초와 공유)

    cumulative_depth = 0.0
    for i, layer in enumerate(SOIL_LAYERS):
        cumulative_depth += layer["thickness"]
        snid = 201 + i
        ops.node(snid, 0.0, -cumulative_depth)
        soil_nodes.append(snid)

    # 기반암 고정 (최하부 지반 노드)
    ops.fix(soil_nodes[-1], 1, 1, 1)

    # 기초 고정 여부: SSI 모델 → 스프링으로 연결
    # 구조물 기초 노드는 지반 상단 노드와 연결 (제약 없음)
    # 구조물 상단 자유도 유지

    # ── 경계조건 ─────────────────────────────────────────────────────
    # 기초층 (floor=0): 지반 스프링으로 연결 (직접 고정 안 함)
    # 지반 최하단: 고정 (already done)
    # 구조물 기초와 지반 표면 연결 (equalDOF)
    for col in range(NUM_BAYS + 1):
        foundation_node = node_id[(0, col)]
        ops.fix(foundation_node, 0, 1, 1)  # 수직·회전 고정, 수평 자유

    # ── 기하변환 ─────────────────────────────────────────────────────
    ops.geomTransf('PDelta', 1)   # P-Delta 효과 포함
    ops.geomTransf('Linear', 2)   # 보: 선형

    # ── 기둥 요소 ────────────────────────────────────────────────────
    ele_id = 1
    int_pts = 5   # 적분점 수
    for floor in range(NUM_STORIES):
        for col in range(NUM_BAYS + 1):
            n_i = node_id[(floor, col)]
            n_j = node_id[(floor + 1, col)]
            ops.element('nonlinearBeamColumn', ele_id, n_i, n_j,
                        int_pts, 1, 1)
            ele_id += 1

    # ── 보 요소 ──────────────────────────────────────────────────────
    for floor in range(1, NUM_STORIES + 1):
        for bay in range(NUM_BAYS):
            n_i = node_id[(floor, bay)]
            n_j = node_id[(floor, bay + 1)]
            ops.element('elasticBeamColumn', ele_id, n_i, n_j,
                        beam_A, Ec, beam_Iz, 2)
            ele_id += 1

    # ── 지반 기둥 요소 (등가 전단 빔) ────────────────────────────────
    for i in range(len(soil_nodes) - 1):
        n_i = soil_nodes[i]
        n_j = soil_nodes[i + 1]
        layer = SOIL_LAYERS[i]
        rho   = layer["gamma"] * 1000 / 9.81
        G     = rho * layer["Vs"] ** 2
        A_soil = 1.0          # 단위 면적 [m²]
        I_soil = 1e-4         # 무시할 관성모멘트
        ops.element('elasticBeamColumn', ele_id, n_i, n_j,
                    A_soil, G, I_soil, 2)
        ele_id += 1

    # ── SSI 스프링 (기초-지반 접면) ───────────────────────────────────
    # 각 기초 노드와 지반 상단 사이에 수평 스프링
    rho_top = SOIL_LAYERS[0]["gamma"] * 1000 / 9.81
    G_top   = rho_top * SOIL_LAYERS[0]["Vs"] ** 2
    k_ssi   = G_top * BAY_WIDTH / SOIL_LAYERS[0]["thickness"]

    ops.uniaxialMaterial('Elastic', 99, k_ssi)
    for col in range(NUM_BAYS + 1):
        fn = node_id[(0, col)]
        ops.element('zeroLength', ele_id, fn, 200, '-mat', 99, '-dir', 1)
        ele_id += 1

    # ── 질량 할당 ────────────────────────────────────────────────────
    rho_c   = 2400       # 콘크리트 밀도 [kg/m³]
    for floor in range(1, NUM_STORIES + 1):
        # 각 층 질량 = 보 질량 + 슬래브 질량
        floor_mass = (
            rho_c * beam_A * BAY_WIDTH * NUM_BAYS    # 보
            + rho_c * 0.2 * BAY_WIDTH * NUM_BAYS     # 슬래브 (200mm)
            + 500 * BAY_WIDTH * NUM_BAYS              # 활하중 환산 질량
        )
        mass_per_node = floor_mass / (NUM_BAYS + 1)
        for col in range(NUM_BAYS + 1):
            nid = node_id[(floor, col)]
            ops.mass(nid, mass_per_node, mass_per_node, 0.0)

    return node_id, soil_nodes


# ─────────────────────────────────────────────────────────────────────
# 3. 고유치 해석 (모달 분석)
# ─────────────────────────────────────────────────────────────────────
def modal_analysis(node_id):
    """고유진동수·모드형상 계산"""
    num_modes = 3
    lam = ops.eigen(num_modes)
    freqs_rad = [np.sqrt(abs(l)) for l in lam]
    freqs_hz  = [w / (2 * np.pi) for w in freqs_rad]
    periods   = [1.0 / f if f > 0 else 0 for f in freqs_hz]
    print("\n━━━ 고유치 해석 결과 ━━━")
    for i, (T, f) in enumerate(zip(periods, freqs_hz)):
        print(f"  모드 {i+1}: T = {T:.3f} s  (f = {f:.2f} Hz)")
    return periods, freqs_hz


# ─────────────────────────────────────────────────────────────────────
# 4. 비선형 시간이력 해석
# ─────────────────────────────────────────────────────────────────────
def time_history_analysis(node_id, soil_nodes, t_arr, acc_arr):
    """
    Newmark β법 기반 비선형 시간이력 해석
    - 감쇠: Rayleigh (1·3모드 기준, ξ=5%)
    - 하중: UniformExcitation (지진파 입력)
    """
    dt_analysis = t_arr[1] - t_arr[0]

    # Rayleigh 감쇠 설정 (1·3모드)
    xi   = 0.05
    lam  = ops.eigen(3)
    w1   = np.sqrt(abs(lam[0]))
    w3   = np.sqrt(abs(lam[2]))
    a0   = 2 * xi * w1 * w3 / (w1 + w3)
    a1   = 2 * xi / (w1 + w3)
    ops.rayleigh(a0, 0, 0, a1)

    # 중력 하중 (정적)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for floor in range(1, NUM_STORIES + 1):
        floor_wt = (
            2400 * beam_A * BAY_WIDTH * NUM_BAYS
            + 2400 * 0.2 * BAY_WIDTH * NUM_BAYS
            + 500  * BAY_WIDTH * NUM_BAYS
        ) * 9.81
        load_per_node = -floor_wt / (NUM_BAYS + 1)
        for col in range(NUM_BAYS + 1):
            ops.load(node_id[(floor, col)], 0.0, load_per_node, 0.0)

    # 정적 해석 (중력)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-6, 100)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / 10)
    ops.analysis('Static')
    ops.analyze(10)
    ops.loadConst('-time', 0.0)

    # 지진파 TimeSeries 등록
    ops.timeSeries('Path', 2, '-dt', dt_analysis,
                   '-values', *acc_arr.tolist(),
                   '-factor', 1.0)
    ops.pattern('UniformExcitation', 2, 1, '-accel', 2)

    # 동적 해석 설정
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-8, 50)
    ops.algorithm('KrylovNewton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')

    # 결과 저장
    results = {
        "time": [],
        "disp_roof": [],       # 옥상 수평 변위
        "disp_story": {i: [] for i in range(1, NUM_STORIES + 1)},
        "drift": {i: [] for i in range(1, NUM_STORIES + 1)},
        "accel_roof": [],
        "soil_disp": [],
    }

    n_steps   = len(t_arr) - 1
    roof_node = node_id[(NUM_STORIES, 1)]
    base_node = node_id[(0, 1)]

    prev_roof_disp = 0.0
    ok = 0
    print("\n━━━ 시간이력 해석 진행 ━━━")
    for step in range(n_steps):
        ok = ops.analyze(1, dt_analysis)
        if ok != 0:
            # 수렴 실패 시 알고리즘 전환
            ops.algorithm('ModifiedNewton', '-initial')
            ok = ops.analyze(1, dt_analysis)
            ops.algorithm('KrylovNewton')
            if ok != 0:
                print(f"  ⚠ 수렴 실패: step {step}, t={t_arr[step]:.2f}s")

        t_cur = ops.getTime()
        results["time"].append(t_cur)

        # 각 층 변위
        roof_d = ops.nodeDisp(roof_node, 1)
        results["disp_roof"].append(roof_d)
        for floor in range(1, NUM_STORIES + 1):
            d = ops.nodeDisp(node_id[(floor, 1)], 1)
            results["disp_story"][floor].append(d)

        # 층간 변위비 (Inter-story Drift Ratio)
        for floor in range(1, NUM_STORIES + 1):
            if floor == 1:
                d_bot = ops.nodeDisp(node_id[(0, 1)], 1)
            else:
                d_bot = ops.nodeDisp(node_id[(floor - 1, 1)], 1)
            d_top = ops.nodeDisp(node_id[(floor, 1)], 1)
            drift = (d_top - d_bot) / STORY_HEIGHT
            results["drift"][floor].append(drift)

        # 지반 변위
        results["soil_disp"].append(ops.nodeDisp(200, 1))

        if step % 200 == 0:
            progress = (step / n_steps) * 100
            print(f"  진행: {progress:.0f}% (t={t_cur:.2f}s, "
                  f"roof={roof_d*1000:.1f}mm)")

    print("  완료!")
    return results


# ─────────────────────────────────────────────────────────────────────
# 5. 안정성 평가
# ─────────────────────────────────────────────────────────────────────
def stability_assessment(results, periods):
    """
    KDS 기준 안정성 평가
    - 층간변위비 한계: 0.5% (즉시거주), 1.5% (인명안전), 2.5% (붕괴방지)
    - 최대 지붕 변위
    - 잔류 변형
    """
    print("\n━━━ 구조물 안정성 평가 ━━━")
    print(f"  기본 주기: T1 = {periods[0]:.3f} s")

    # 층간 변위비 최대값
    max_drifts = {}
    for floor in range(1, NUM_STORIES + 1):
        max_drift = max(abs(d) for d in results["drift"][floor])
        max_drifts[floor] = max_drift

    # 한계값 정의 (KDS 41 17 00)
    limits = {
        "IO": 0.005,   # Immediate Occupancy
        "LS": 0.015,   # Life Safety
        "CP": 0.025,   # Collapse Prevention
    }

    assessment = {}
    for floor, drift in max_drifts.items():
        if drift <= limits["IO"]:
            status = "✅ 즉시거주 (IO)"
            level  = "IO"
        elif drift <= limits["LS"]:
            status = "⚠️ 인명안전 (LS)"
            level  = "LS"
        elif drift <= limits["CP"]:
            status = "🔶 붕괴방지 (CP)"
            level  = "CP"
        else:
            status = "❌ 붕괴 위험"
            level  = "FAIL"
        assessment[floor] = {"drift": drift, "status": status, "level": level}
        print(f"  {floor}층 최대 층간변위비: {drift*100:.3f}% → {status}")

    # 최대 지붕 변위
    max_roof = max(abs(d) for d in results["disp_roof"])
    H_total  = NUM_STORIES * STORY_HEIGHT
    roof_ratio = max_roof / H_total
    print(f"\n  최대 지붕 변위: {max_roof*1000:.1f} mm "
          f"(H/{H_total/max_roof:.0f})")

    # 잔류 변형
    residual = abs(results["disp_roof"][-1])
    print(f"  잔류 변형: {residual*1000:.2f} mm")

    return assessment, max_drifts, limits


# ─────────────────────────────────────────────────────────────────────
# 6. 시각화
# ─────────────────────────────────────────────────────────────────────
def visualize_results(t_arr, acc_arr, results, assessment, max_drifts,
                      limits, periods):
    """결과 종합 시각화 (6 패널)"""
    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor('#0f1117')
    gs  = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

    PANEL_BG = '#1a1d27'
    C_MAIN   = '#4fc3f7'   # 청색
    C_WARN   = '#ffb74d'   # 주황
    C_DANGER = '#ef5350'   # 적색
    C_OK     = '#66bb6a'   # 녹색
    C_SOIL   = '#a5d6a7'   # 연두

    def style_ax(ax, title):
        ax.set_facecolor(PANEL_BG)
        ax.set_title(title, color='white', fontsize=11, fontweight='bold', pad=8)
        ax.tick_params(colors='#aaaaaa', labelsize=8)
        ax.spines['bottom'].set_color('#444')
        ax.spines['top'].set_color('#444')
        ax.spines['left'].set_color('#444')
        ax.spines['right'].set_color('#444')
        ax.xaxis.label.set_color('#aaaaaa')
        ax.yaxis.label.set_color('#aaaaaa')

    # ① 입력 지진파
    ax1 = fig.add_subplot(gs[0, :2])
    style_ax(ax1, '① 입력 지진파 가속도 이력')
    pga = max(abs(acc_arr)) / 9.81
    ax1.plot(t_arr, acc_arr / 9.81, color=C_MAIN, lw=0.8, alpha=0.9)
    ax1.axhline(0, color='#444', lw=0.5)
    ax1.fill_between(t_arr, acc_arr / 9.81, 0,
                     alpha=0.2, color=C_MAIN)
    ax1.set_xlabel('시간 [s]')
    ax1.set_ylabel('가속도 [g]')
    ax1.text(0.98, 0.92, f'PGA = {pga:.3f} g',
             transform=ax1.transAxes, color=C_WARN,
             ha='right', fontsize=10, fontweight='bold')

    # ② 지붕 변위 이력
    ax2 = fig.add_subplot(gs[0, 2])
    style_ax(ax2, '② 지붕 수평 변위')
    time_r = results["time"]
    disp_mm = [d * 1000 for d in results["disp_roof"]]
    ax2.plot(time_r, disp_mm, color=C_MAIN, lw=1.0)
    ax2.axhline(0, color='#444', lw=0.5)
    ax2.set_xlabel('시간 [s]')
    ax2.set_ylabel('변위 [mm]')
    max_d = max(abs(d) for d in disp_mm)
    ax2.text(0.98, 0.92, f'Max: {max_d:.1f} mm',
             transform=ax2.transAxes, color=C_WARN,
             ha='right', fontsize=9)

    # ③ 층간 변위비 이력 (3개 층 오버레이)
    ax3 = fig.add_subplot(gs[1, :2])
    style_ax(ax3, '③ 층간 변위비 이력 (Inter-story Drift)')
    colors_floor = [C_OK, C_WARN, C_DANGER]
    for floor in range(1, NUM_STORIES + 1):
        drift_pct = [d * 100 for d in results["drift"][floor]]
        ax3.plot(time_r, drift_pct,
                 color=colors_floor[floor - 1],
                 lw=0.9, alpha=0.85,
                 label=f'{floor}층')
    # 한계선
    ax3.axhline(limits["IO"] * 100, color='white',
                lw=1.0, ls='--', alpha=0.4, label='IO 한계 0.5%')
    ax3.axhline(limits["LS"] * 100, color=C_WARN,
                lw=1.0, ls='--', alpha=0.6, label='LS 한계 1.5%')
    ax3.axhline(limits["CP"] * 100, color=C_DANGER,
                lw=1.0, ls='--', alpha=0.6, label='CP 한계 2.5%')
    ax3.set_xlabel('시간 [s]')
    ax3.set_ylabel('층간변위비 [%]')
    ax3.legend(loc='upper right', fontsize=7,
               facecolor='#222', edgecolor='#444',
               labelcolor='white')

    # ④ 최대 층간 변위비 바 차트 + 안정성 판정
    ax4 = fig.add_subplot(gs[1, 2])
    style_ax(ax4, '④ 안정성 판정')
    floors     = list(range(1, NUM_STORIES + 1))
    drift_vals = [max_drifts[f] * 100 for f in floors]
    bar_colors = []
    for f in floors:
        lv = assessment[f]["level"]
        bar_colors.append(
            C_OK if lv == "IO" else
            C_WARN if lv == "LS" else
            C_DANGER
        )
    bars = ax4.barh(floors, drift_vals, color=bar_colors,
                    height=0.5, alpha=0.85)
    ax4.axvline(limits["IO"] * 100, color='white',
                lw=1, ls='--', alpha=0.5)
    ax4.axvline(limits["LS"] * 100, color=C_WARN,
                lw=1, ls='--', alpha=0.7)
    ax4.axvline(limits["CP"] * 100, color=C_DANGER,
                lw=1, ls='--', alpha=0.7)
    ax4.set_xlabel('층간변위비 [%]')
    ax4.set_ylabel('층')
    ax4.set_yticks(floors)
    ax4.set_yticklabels([f'{f}층' for f in floors], color='white')
    for bar, val in zip(bars, drift_vals):
        ax4.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                 f'{val:.3f}%', va='center', color='white', fontsize=8)

    # ⑤ 지반-구조물 상대 변위
    ax5 = fig.add_subplot(gs[2, :2])
    style_ax(ax5, '⑤ 지반-구조물 상대 변위 (SSI 효과)')
    soil_mm   = [d * 1000 for d in results["soil_disp"]]
    struct_mm = [d * 1000 for d in results["disp_story"][1]]
    rel_mm    = [s - g for s, g in zip(struct_mm, soil_mm)]
    ax5.plot(time_r, soil_mm,   color=C_SOIL, lw=0.9, alpha=0.8,
             label='지반 표면 변위')
    ax5.plot(time_r, struct_mm, color=C_MAIN, lw=0.9, alpha=0.8,
             label='구조물 1층 변위')
    ax5.plot(time_r, rel_mm,    color=C_WARN, lw=0.9, alpha=0.8,
             ls='--', label='상대 변위 (SSI)')
    ax5.axhline(0, color='#444', lw=0.5)
    ax5.set_xlabel('시간 [s]')
    ax5.set_ylabel('변위 [mm]')
    ax5.legend(loc='upper right', fontsize=7,
               facecolor='#222', edgecolor='#444', labelcolor='white')

    # ⑥ 변형 모드 개요 (구조물 다이어그램)
    ax6 = fig.add_subplot(gs[2, 2])
    style_ax(ax6, '⑥ 최대 변형 형상')
    ax6.set_xlim(-1.5, NUM_BAYS * BAY_WIDTH + 1.5)
    ax6.set_ylim(-2, (NUM_STORIES + 0.5) * STORY_HEIGHT)
    ax6.set_aspect('equal')
    ax6.axis('off')

    # 지반 블록
    soil_h_draw = 2.0
    ax6.add_patch(mpatches.FancyBboxPatch(
        (-0.5, -soil_h_draw), NUM_BAYS * BAY_WIDTH + 1.0, soil_h_draw,
        boxstyle="round,pad=0.1", fc='#2e4a2e', ec=C_SOIL, lw=1.0, alpha=0.7))
    ax6.text(NUM_BAYS * BAY_WIDTH / 2, -soil_h_draw / 2,
             '다층 지반', ha='center', va='center',
             color=C_SOIL, fontsize=8)

    # 스케일 변위 (과장 50배)
    scale = 50.0
    max_t = results["time"].index(
        max(results["time"],
            key=lambda t: abs(results["disp_roof"][results["time"].index(t)])))
    peak_idx = np.argmax(np.abs(results["disp_roof"]))

    deformed_x = {}
    for floor in range(NUM_STORIES + 1):
        if floor == 0:
            dx = results["soil_disp"][peak_idx] * scale
        else:
            dx = results["disp_story"][floor][peak_idx] * scale
        for col in range(NUM_BAYS + 1):
            x_orig = col * BAY_WIDTH
            deformed_x[(floor, col)] = x_orig + dx

    # 기둥 그리기
    for floor in range(NUM_STORIES):
        for col in range(NUM_BAYS + 1):
            x1 = deformed_x[(floor, col)]
            x2 = deformed_x[(floor + 1, col)]
            y1 = floor * STORY_HEIGHT
            y2 = (floor + 1) * STORY_HEIGHT
            ax6.plot([x1, x2], [y1, y2], color=C_MAIN, lw=2.5, alpha=0.9)

    # 보 그리기
    for floor in range(1, NUM_STORIES + 1):
        for bay in range(NUM_BAYS):
            x1 = deformed_x[(floor, bay)]
            x2 = deformed_x[(floor, bay + 1)]
            y  = floor * STORY_HEIGHT
            ax6.plot([x1, x2], [y, y], color=C_WARN, lw=2.0, alpha=0.9)

    # 층 레이블
    for floor in range(1, NUM_STORIES + 1):
        ax6.text(NUM_BAYS * BAY_WIDTH + 0.8,
                 floor * STORY_HEIGHT - STORY_HEIGHT / 2,
                 f'{floor}F', color='#aaa', fontsize=7, va='center')

    ax6.text(NUM_BAYS * BAY_WIDTH / 2, -0.3,
             f'(변형 {scale:.0f}배 과장)',
             ha='center', color='#777', fontsize=7)

    # 전체 제목
    fig.suptitle(
        f'OpenSeesPy  |  {NUM_STORIES}층 RC 프레임 + 비선형 다층 지반  |  '
        f'지진 PGA = {SCALE_PGA}g',
        color='white', fontsize=13, fontweight='bold', y=0.98
    )

    out_path = '/mnt/user-data/outputs/seismic_analysis_result.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"\n  그래프 저장 완료: {out_path}")
    return out_path


# ─────────────────────────────────────────────────────────────────────
# 7. 요약 리포트 출력
# ─────────────────────────────────────────────────────────────────────
def print_summary(results, assessment, periods, t_arr, acc_arr):
    pga = max(abs(acc_arr)) / 9.81
    print("\n" + "="*55)
    print("   OpenSeesPy 구조-지반 복합 지진 해석 요약 리포트")
    print("="*55)
    print(f"  구조물   : {NUM_STORIES}층 RC 모멘트 프레임 ({NUM_BAYS}스팬)")
    print(f"  지반     : 비선형 다층 지반 ({len(SOIL_LAYERS)}층)")
    print(f"  입력 지진: PGA = {pga:.3f} g")
    print(f"  기본 주기: T1 = {periods[0]:.3f} s")
    print("-"*55)
    overall = "통과"
    for f, info in assessment.items():
        print(f"  {f}층: {info['status']}  (IDR={info['drift']*100:.3f}%)")
        if info["level"] in ("CP", "FAIL"):
            overall = "주의 필요"
    print("-"*55)
    print(f"  종합 판정: {overall}")
    print("="*55)


# ─────────────────────────────────────────────────────────────────────
# 8. 메인 실행
# ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("━━━ 지진파 생성 ━━━")
    t_arr, acc_arr = generate_earthquake(dt, t_total, pga_g=SCALE_PGA)

    print("━━━ 모델 구축 ━━━")
    node_id, soil_nodes = build_model()

    print("━━━ 고유치 해석 ━━━")
    periods, freqs = modal_analysis(node_id)

    print("━━━ 시간이력 해석 ━━━")
    ops.wipe()
    node_id, soil_nodes = build_model()   # 모델 재구축 (eigen 후 초기화)
    results = time_history_analysis(node_id, soil_nodes, t_arr, acc_arr)

    print("━━━ 안정성 평가 ━━━")
    assessment, max_drifts, limits = stability_assessment(results, periods)

    print("━━━ 시각화 ━━━")
    img_path = visualize_results(
        t_arr, acc_arr, results, assessment, max_drifts, limits, periods)

    print_summary(results, assessment, periods, t_arr, acc_arr)
    print(f"\n✅ 결과 이미지: {img_path}")
