import openseespy.opensees as ops
import numpy as np

# [1] 클래스 정의: 지구 환경 엔진 (Earth Environment Engine)
class EarthEnvironment:
    """
    지구 환경 물리 상수를 관리하고 파생 물리량을 계산하는 클래스
    단위계: SI (m, N, kg, s, °C) 사용
    """
    def __init__(self, temperature_now=15.0, wind_speed=10.0):
        # --- (A) 고정 파라미터: 지구의 불변 상수 (Immutable Constants) ---
        self.g = 9.81                 # 중력가속도 [m/s^2]
        self.patm = 101325.0          # 표준 대기압 [Pa] (1 atm)
        self.air_density = 1.225      # 공기 밀도 [kg/m^3] (15°C 해수면 기준)
        self.thermal_alpha = 1.2e-5   # 선열팽창계수 [1/°C] (Structural Steel/Concrete 근사치)
        
        # --- (B) 범위형/변동 파라미터: 환경에 따라 재정의 가능 (Variable States) ---
        self.rho = 1800.0             # 지반 밀도 [kg/m^3] (Bulk Density)
        self.nu = 0.30                # 포아송비 [-] (Poisson's Ratio)
        self.Vs = 200.0               # 전단파 속도 [m/s] (지반 강성 결정 인자)
        self.temp_ref = 15.0          # 해석 기준 온도 [°C] (무응력 상태 온도)
        
        # 실시간 변동 변수
        self.temp_now = temperature_now # 현재 외기 온도 [°C]
        self.wind_speed = wind_speed  # 현재 풍속 [m/s]
        self.Cd = 1.2                 # 구조물 항력 계수 [-]
        
        # 초기 파생 물리량 계산 호출
        self._update_derived_properties()

    def _update_derived_properties(self):
        """환경 변화(온도, 풍속 등)에 따른 종속 공학 물리량 갱신"""
        # G (전단탄성계수) = ρ * Vs^2 [Pa]
        self.G = self.rho * (self.Vs ** 2)
        # E (탄성계수) = 2G(1+ν) [Pa] (Young's Modulus)
        self.E = 2.0 * self.G * (1.0 + self.nu)
        # dT (온도 변화량) = 현재온도 - 기준온도 [°C]
        self.dT = self.temp_now - self.temp_ref
        # q (설계 풍압) = 0.5 * ρ_air * V^2 [Pa]
        self.q_wind = 0.5 * self.air_density * (self.wind_speed ** 2)
        
# [2] OpenSees 초기 모델링 (반드시 루프 밖에서 정의)
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

# [3] 모델 기하학적 정의 (Node & Element)
# 분석(analyze)을 돌리기 전에 반드시 계산 대상이 존재해야 합니다.
ops.node(1, 0.0, 0.0, 0.0)    # 1번 노드 (좌표: 0, 0, 0)
ops.node(2, 1.0, 0.0, 0.0)    # 2번 노드 (좌표: 1, 0, 0)
ops.fix(1, 1, 1, 1, 1, 1, 1)  # 1번 노드 완전 고정

env = EarthEnvironment()
# nDMaterial 정의: (태그, 유형, E, nu, rho)
ops.nDMaterial('ElasticIsotropic', 1, env.E, env.nu, env.rho)
# 요소 정의: (유형, 태그, 노드1, 노드2, 단면적, 재료태그)
ops.element('truss', 1, 1, 2, 1.0, 1) 

# [4] 타임 시리즈 및 하중 패턴 설정
dt = 600.0                    # 분석 시간 간격 [s] (10분)
total_time = 86400.0          # 전체 시뮬레이션 시간 [s] (24시간)
time_steps = np.arange(0, total_time + dt, dt)

# Path 데이터 미리 생성 (if문 활용)
temp_values = [25.0 if 21600 <= t <= 64800 else 5.0 for t in time_steps]

# -time 옵션은 삭제하여 -dt와 충돌 방지
ops.timeSeries('Path', 1, '-dt', dt, '-values', temp_values)
ops.pattern('Plain', 1, 1)
ops.load(2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0) # 기준 하중 설정

# [5] 분석 알고리즘 설정
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-6, 10)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0)
ops.analysis('Static')

# [6] 통합 분석 루프 (Rule-based Environment Control)
print(f"{'Time(h)':>8} | {'State':>6} | {'Temp(°C)':>8} | {'Vs(m/s)':>8} | {'E(Pa)':>12}")
print("-" * 65)

for i, current_time in enumerate(time_steps[:-1]):
    # OpenSees 내부 시간 동기화
    ops.setTime(current_time)
    
    # --- 실시간 환경 규칙 적용 (낮/밤 분기) ---
    if 21600.0 <= current_time <= 64800.0:
        state, env.temp_now, env.Vs = "DAY", 25.0, 200.0
    else:
        state, env.temp_now, env.Vs = "NIGHT", 5.0, 220.0
        
    # 물리량 갱신 (E값 등 재계산)
    env._update_derived_properties()
    
    # 분석 실행 (1 Step)
    ok = ops.analyze(1)
    
    if ok == 0:
        hour = current_time / 3600.0
        print(f"{hour:>8.1f} | {state:>6} | {env.temp_now:>8.1f} | {env.Vs:>8.1f} | {env.E:>12.2e}")
    else:
        print(f"Error at {current_time} sec"); break

print("-" * 65)
print("지구 환경 주기 시뮬레이션 완료!")


