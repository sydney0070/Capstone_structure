# JEZ-1 해석 지연 원인과 후보 선정 판단

2026-11-18. 사용자 요청으로 진행 중인 두 JEZ-1 해석을 모두 중지했다. 도구 세션의 종료 코드 143과 해당 실행 경로의 활성 프로세스가 없음을 확인했다. 아래 비교는 추가 지반 해석 없이 저장된 입력·로그·전단 시험을 후처리한 결과다.

## 판단

현재 0.2–3 m 모래층의 기본 해석 후보에는 **MGS-1을 추천**, **MGS-1S도 대체 후보로 사용 가능**하다. **현재 1,300 kg/m³ JEZ-1 초기 매개변수 조합은 기본 후보 채택을 보류**하고 약한 지반 조건의 비교 후보로 남기는 편이 타당하다. JEZ-1 모사토라는 재료 자체를 부적합하다고 판정한 것은 아니다. 세 후보 모두 InSight 밀도·예측 강성과 다른 모사토의 직접전단 강도를 조합한 초기 재료이며 현장 보정을 완료하지 않았다.

## 같은 모델에서의 결과

공통: 밀도 1,300 kg/m³, E 28.4 MPa, ν 0.22, ψ 0°, Drucker–Prager 무경화 모델. 영역 90 × 90 × 30 m, 4,093절점, 3,504요소, 지반 총질량 542,431,655.75 kg가 같다. 후보별 변경은 0.2–3 m 층의 마찰각과 점착력이다.

| 후보 | φ | c | 실행 상태 | 감압 최대 침하 | 가압 최대 침하 |
|---|---:|---:|---|---:|---:|
| MGS-1 | 35.5° | 280 Pa | 전체 완료, 약 5.7분 | 1.4338 mm | 0.2827 mm |
| MGS-1S | 34.0° | 280 Pa | 전체 완료, 약 6.5분 | 1.4339 mm | 0.2828 mm |
| JEZ-1 | 24.0° | 230 Pa | 약 21.3분 뒤 사용자 중지, 자중 49.20% | 미산출 | 미산출 |

감압은 내압 효과를 빼기 전 하향 150.3 MN을 가한 상태이고, 가압은 120.72 MN만큼 하중을 줄인 상태다. 침하는 굴착 지반의 자중 평형을 기준으로 한 링 바닥의 증분값이다. 동시 실행 조건이 달라 위 시간은 엄밀한 성능 벤치마크가 아니다. 반복 횟수와 재시도 기록을 함께 사용해 지연 원인을 판단했다.

## 확인된 지연 원인

1. **비선형 평형 계산의 반복 실패.** MGS-1·MGS-1S는 자중 20증분을 모두 통과했고 실패가 0회였다. 성공 증분당 평균 반복은 각각 5.9회와 7.0회였다. JEZ-1 증분 축소 실행은 자중 절반에 못 미친 상태에서 성공 62회와 실패 13회가 발생했다. 실패한 13회에는 약 560.8초를 썼다. 완료된 증분 계산 시간의 약 45.7%가 실패한 계산이었다.
2. **필요한 증분이 계속 작아짐.** JEZ-1 초기 증분은 자중의 1.25%였지만 마지막 성공 증분은 약 0.301%였다. 30% 부근 실패 이후 증분을 절반으로 줄이면 다시 진행됐고 최종 49.20%까지 갔다. 초기 실패 하중을 물리적 한계 하중으로 볼 수 없다. 기본 해법의 실패 시 세 알고리즘을 최대 80회씩 시도하는 정책도 시간을 늘렸다. 비교용 Krylov 차원 확대 실행은 자중 34.38%에서 중지했으며 확실한 속도 개선을 입증하지 못했다.
3. **엄격한 오차 기준도 비용에 기여.** 자중의 NormUnbalance 허용값은 약 0.2013 N이다. 실패 기록에는 잔차가 0.612 N 또는 1.001 N인 경우와 10⁶ N 수준까지 커진 경우가 함께 있다. 반복 상한이나 허용오차에 가까워서 중단된 시도와 뚜렷하게 발산한 시도를 구분해야 한다. 이번에는 허용오차를 완화하지 않았다. 전역 반력의 합이 맞는 것과 모든 자유도의 불평형력 노름이 작은 것은 별개의 검사다.

## 재료·모델 차이의 해석

- **가장 유력한 입력 차이는 낮은 강도다.** 같은 밀도와 탄성강성에서 JEZ-1의 φ와 c가 모두 작다. 문헌 판독값의 Mohr–Coulomb 포락선 τ=c+σn tanφ에 σn=0.5 kPa를 대입하면 MGS-1 0.6366 kPa, MGS-1S 0.6173 kPa, JEZ-1 0.4526 kPa다. JEZ-1은 이 응력에서 MGS-1보다 약 28.9% 낮다. 이 식은 이번 DP 모델의 모든 응력 경로에 그대로 적용되는 강도식은 아니다.
- **무경화·영 팽창각 가정은 수렴을 어렵게 할 수 있다.** 세 재료 모두 경화 계수들이 0이고 ψ=0인 비연관 흐름을 사용한다. 낮은 강도와 결합해 소성 변형이 일찍 발생하거나 집중되는 것이 현재의 가장 유력한 설명이다. 다만 저장된 전체 지반 로그에는 적분점 소성변형률·항복 상태 이력이 없어 실제 소성 집중 위치와 파괴 메커니즘은 확인하지 못했다. 추가 해석 없이 이를 확정하지 않는다.
- **굴착 초기화도 검토 대상이다.** 현재 코드는 링 부분을 깊이 3 m까지 처음부터 제거한 형상에 자중을 증가시킨다. 굴착 벽을 지지하는 기초 구조물이나 접촉 요소를 넣지 않고, 뒤 단계에서 링 바닥에만 등가 하중을 가한다. 저구속 굴착 경계의 영향과 초기응력 경로가 수렴에 영향을 줄 가능성이 있다. 자연 지반의 초기 평형 → 굴착 → 기초 설치 순서를 검증한 모델은 아니다.

따라서 지연은 **약한 초기 재료 조합·굴착 초기화·하중 증분 및 수렴 정책이 결합된 수치적 문제**로 보는 것이 타당하다. 현 기록만으로 JEZ-1 지반의 실제 붕괴나 돔 지지 불능을 확정할 수 없다. 중지 때 기록된 최대 전체 변위 약 2.983 mm는 부분 자중에서의 전체 변위이며 완주 후보의 링 침하와 비교할 수 없다.

## 사용 적합성의 범위

MGS-1과 MGS-1S의 침하가 비슷한 이유에는 같은 탄성계수와 같은 하부층이 있다. 링 바닥은 깊이 3 m이며 그 아래는 E=520 MPa의 기존 탄성 조립질층이다. 이번 후보 재료는 그 위 0.2–3 m 층에만 들어간다. 두 결과가 가깝다고 해서 세 후보의 지지력이나 전체 깊이의 재료 성능이 같다는 뜻은 아니다.

세 후보의 단일 요소 전단 비교에서도 현재 삼축압축 대응 DP는 그래프 기반 직접전단 포락선보다 약 16.6–22.0% 큰 저항을 냈다. 소스의 수직응력 범위는 0.097–0.68 kPa다. 실제 모델의 모든 응력을 그 범위 안에서 검증하지 않았다. 따라서 MGS-1/MGS-1S는 **캡스톤의 초기 수치해석 후보로 사용할 수 있는 상태**이고, 실구조물의 안전성·지지력 검증을 통과한 재료라는 판단은 할 수 없다.

JEZ-1을 다시 기본 후보로 검토하려면 직접전단 응력 경로에 맞춘 DP 강도 보정, 실험 근거가 있는 팽창·경화 설정, 초기응력·굴착·기초 설치 단계와 적분점 상태 기록을 먼저 정리하는 것이 좋다. 수렴을 얻기 위해 c나 φ를 임의로 올리지는 않았다. 이번 요청에 따라 재해석은 진행하지 않았다.

![수렴 진행과 증분 축소](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/jez1_solver_progress.png)

## 확인한 자료

- [원시 진단 지표](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/jez1_runtime_diagnosis.json)
- [JEZ-1 21분 실행 요약](/Users/sydney/Desktop/capstone/Capstone_structure/ground/results/ground_dome_footing/analysis_jez1_adaptive_20261118/summary.json)
- [증분별 반복·잔차 기록](/Users/sydney/Desktop/capstone/Capstone_structure/ground/results/ground_dome_footing/analysis_jez1_adaptive_20261118/solver_attempts.jsonl)
- [지반·문헌 대조 보고서](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/report_KO.md)
- [재료 등록 코드](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_hybrid/materials.py:32)
- [굴착·경계조건 코드](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_hybrid/conforming_mesh.py:119)
- [해석 및 하중 적용 코드](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_dome_footing.py:191)
- [OpenSees DP 공식 문서](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html): 경화·소성 체적변화 계수, state 출력 및 하중 제어 예제의 수렴 제한.
- [NormUnbalance 공식 문서](https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/test/NormUnbalance.html): 반복 종료 기준.
