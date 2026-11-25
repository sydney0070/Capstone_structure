# 돔–RC 링–3D 지반 통합 연구 코드

실행 파일: `mars_dome_ring_ground.py`. 기존 돔·황콘크리트·독립 링·지반 파일은 변경하지 않았다.

## 어떻게 연결했는가

- 지반: 기존 ground_hybrid의 3DOF bbarBrick 연속체. 표층 0.2 m InSight duricrust proxy, 8 m까지 JEZ-1, 아래 탄성 현무암. DP 2종은 비선형이다. 평면 지표, 60 × 60 × 12 m 가정 영역이며 DEM 관측소 위치 재현이 아니다.
- 링: 새 ring_foundation_rc.py의 C30/SD400 섬유단면 HingeRadau 비선형 RC 링. 예제 배근·직사각 1 × 0.6 m 단면, R=20 m. 빠른 기본 통합은 12분할이다. 별도 48분할 링의 메시 수렴성을 그대로 승계한 것이 아니다.
- 돔: dome_in_mars_tuned.py의 S355 Steel02 섬유단면을 재사용한다. 43개 돔 점/114개 부재. 기반부 12점은 RC 링의 절점을 **직접 재사용**한다. 원본의 강접합 가정이다.
- 좌표는 동/북/상으로 통일한다. 원본 돔의 북/동 축은 좌표와 압력 벡터를 함께 변환한다.
- 지반 격자에 링 좌표를 명시적으로 포함한다. 임의 최근접 지반점으로 이동시키지 않고 일치 좌표를 확인한다.
- 기본 `--interface tied`: 지반 3DOF와 링 6DOF의 병진 1/2/3만 equalDOF 연결한다. 링 회전은 임의 고정하지 않는다. 인장 전달을 허용하는 연결 기준모델이므로 들림·미끄러짐 결과로 해석하면 안 된다.
- `--interface contact`: 검증 패키지의 flatSliderBearing 기반 수직 접촉/분리 및 Coulomb 마찰을 사용한다. 독립 원지반 스프링을 추가하는 것이 아니라 지반과 링 사이의 접촉 penalty다. 마찰계수 0.35는 가정이다. 링 중심선 점마다 폭×현 길이의 분담면적을 배정한 집중 접촉이며, 실제 유한 폭 기초판 접촉압을 면으로 해상한 모델은 아니다. 큰 미끄러짐·기초 지지력·파괴 검증은 남았다.
- 하나의 모델에 지반·기초·돔을 먼저 조립한 뒤 함께 자중을 재하한다. 돔 코드를 그대로 호출해 기존 지반을 지우지 않는다. 독립 링의 demo_supports 및 기존 돔 지반 스프링은 사용하지 않는다.
- 지반은 요소 밀도 관성 질량, RC 링과 돔은 절점 질량을 쓴다. 돔이 공유하는 링 절점에는 두 질량을 합쳐 한 번 부여한다. 자중도 각각 한 번만 부여한다.

## 빠른 실행

```bash
/Users/sydney/anaconda3/envs/cap/bin/python /Users/sydney/Desktop/capstone/Capstone_structure/structure/mars_dome_ring_ground.py
```

기본은 **무가압·무복토**의 연결 검증이다. 이전 7.5 m 복토/1기압 사례의 완성 설계가 아니다. 최대 가속도 주변 2초를 사용한다. 결과는 results_mars_integrated/에 저장되며, mars_integrated_response.html을 열면 지반 경계·링·돔이 같은 계산 시간축에서 움직인다. HTML은 동적 시작 평형상태 이후의 실제 변위 증분을 확대한다. 기본 표시배율 1,000,000은 **화면 전용**이고 CSV는 실제 m다. 실제 변위는 HTML 프레임 제목과 CSV에서 확인한다. 정적 실행의 HTML은 동적 변위 증분이 없으므로 변형 0을 표시한다; 정적 최대 변위는 summary.json에 있다.

접촉을 포함한 실행:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python /Users/sydney/Desktop/capstone/Capstone_structure/structure/mars_dome_ring_ground.py --interface contact --output-dir /Users/sydney/Desktop/capstone/Capstone_structure/structure/results_mars_integrated/contact_run
```

실측 파형을 증폭한 **가상 강진** 예제:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python /Users/sydney/Desktop/capstone/Capstone_structure/structure/mars_dome_ring_ground.py --interface contact --motion-scale 1000 --visual-scale 1000 --output-dir /Users/sydney/Desktop/capstone/Capstone_structure/structure/results_mars_integrated/strong_run
```

입력 1,000배를 실제 관측 지진이나 특정 재현주기 지진으로 표현하지 않는다. 비선형 재료가 활성화되어 있어도 소성·손상이 크게 발생했다는 뜻은 아니며, 이번 출력만으로 파괴를 판정하지 않는다.

추가 옵션: `--duration 10`, `--motion-start 초`, `--static-only`, `--no-html`, `--net-pressure Pa`, `--cover-thickness m`, `--mass-ton`, `--spacing m`, `--ground-width m`, `--ring-segments 24`. 압력/복토를 넣는다고 이전에 완료한 독립 스프링 모델의 가압 경로가 보장되지는 않는다. 새로운 하중 경로의 수렴·평형·접촉 상태를 확인해야 한다. 폭·배근·RC 기초 형상 확정, 복토의 지반 측 하중 및 시공 이력은 후속 검토 대상이다.

## 입력·경계 및 검증 범위

S1222a_ACC를 보정 완료 m/s² 제품으로 가정해 기존 0.05–2 Hz 필터 및 10 Hz 재표본화를 적용하고 최대 벡터 가속도 주변 구간을 추출한다. 세 성분의 UniformExcitation을 **통합 모델 전체에 한 번씩만** 적용한다. 같은 입력을 지반 출력에 다시 적용하지 않는다. 응답은 공통 기반 운동에 대한 상대 응답이다. 지표 관측을 기반 입력으로 사용하는 실험이며 역산 기반암/입사파 입력이 아니다.

현재는 고정 바닥·측면 롤러의 유한 경계다. 파동 반사·공간/시간 해상도·영역 및 링 메시 수렴·감쇠가 충분히 검증된 실제 지진 현장 모델이 아니다. 기본 2초 구간은 강한 부분을 빠르게 확인하는 용도이고, 선행 진동의 이력과 전체 기록 최대 응답을 포함한다는 보장은 없다. 추가 Rayleigh 감쇠는 자동으로 넣지 않는다. 지반/강재/RC 구성식의 이력효과는 유지한다.

C30/SD400은 명목 prototype, S355는 기존 HEM320 시편 계수의 SHS 전용 가정, 지반은 문헌 및 proxy 상수다. 전 재료의 실제 구조물·화성 현장 보정 완료가 아니다. RC 전단/정착/비틀림 파괴, SHS 국부좌굴·접합부, 실제 기초 폭의 접촉압·지지력·큰 미끄러짐 및 내압 밀폐 경로도 검증하지 않았다.

실행 성공은 run_status.json의 complete=true로 확인한다. 실패하면 성공 HTML로 간주하지 않으며 과거 HTML이 남아 있어도 현재 결과로 쓰지 않는다. SAFE나 최소 질량을 출력하지 않는다.

## 이번에 완료한 확인

기본 조립: 지반 solid 1,176개, RC 링 12개 보, 튜닝 강재 돔 114개 보.

- 병진 연결 정적: 약 11.9초, 전체 힘 평형 상대오차 약 3.05e-15, 지반–링 병진 차이 0.
- 접촉 정적: 약 13.3초, 힘 평형 약 2.04e-15, penalty에 의한 지반–링 상대변위 존재.
- 실측 배율 1, 병진 연결, 2초 20단계 및 HTML: 약 25.6초.
- 실측 배율 1, 접촉, 2초 20단계: 약 23.6초.
- 가상 배율 1,000, 접촉, 2초 20단계 및 HTML: 약 28.8초.

이 시간은 확인한 작은 모델과 기기 상태의 관측값이며 실행시간 보증이 아니다. 큰 입력에서 불수렴이 나면 구조 파괴로 바로 판정하지 않는다. 더 긴 구간·촘촘한 격자·큰 영역은 계산량이 늘어난다.

지반/링/돔의 max increment는 각 집합 절점의 변위 벡터 크기 최대다. apex_minus_mean_base는 돔 꼭대기에서 기초 평균 병진을 뺀 값이며, 전체 강체 회전을 제거한 순수 골조변형은 아니다.

확인 결과와 화면:

- `results_mars_integrated/tied_quake_check/`
- `results_mars_integrated/contact_quake_check/`
- `results_mars_integrated/strong_contact_check/mars_integrated_response.html`

원본 파일들은 그대로 유지했다. 통합을 위해 새 파일만 추가했다.

추가 확인: 무입력 접촉 모델의 20단계 최대 변위 증분은 4.276e-14 m였다. 원래 스프링 요소가 남지 않았고 전체 구조물 절점 질량이 링+돔+막 질량 합계와 일치함을 별도 조립 검사로 확인했다. 근거는 `results_mars_integrated/verification_manifest.json`에 있다.
