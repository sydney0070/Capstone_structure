# 초기 레골리스 후보의 3D 해석 및 문헌 대조

2026-11-18 실행. 대상은 90 × 90 m, 깊이 30 m의 기존 층상 지반이다. 0.2–3 m 모래층에 후보별 1,300 kg/m³ DP 재료를 적용했다. 매입 깊이 3 m, 폭 14 m 링 바닥에 하향 총하중 150.3 MN을 걸고, 이어서 내압 순들림 120.72 MN만큼 하중을 줄였다.

## 3D 계산 결과

| 후보 | 계산 상태 | 감압 최대 침하 | 가압 최대 침하 |
|---|---|---:|---:|
| MGS-1 | 완료 | 1.4338 mm | 0.2827 mm |
| MGS-1S | 완료 | 1.4339 mm | 0.2828 mm |
| JEZ-1 | 추가 계산 후 자중 49.20%에서 사용자 중지 | — | — |

MGS-1과 MGS-1S는 자중, 감압, 가압의 목표 하중배율 1.0을 모두 달성했고 지점 반력의 상대 평형 오차는 4×10⁻¹⁵ 이하였다. 침하는 굴착 지반의 자중 평형 이후 링 바닥 위치의 증분값이다. 기초 지지 깊이 3 m 아래는 기존 탄성 조립질 레골리스이므로 이번 재료 변경은 그 위의 모래층에 적용된다.

JEZ-1의 기본 20증분 해석은 자중의 15.0%를 확정한 상태에서 600초 제한에 도달했다. 80증분 Newton 재시도는 23.8%를 확정한 상태에서 같은 제한에 도달했다. 두 실행 모두 기초 재하 단계에 도달하지 못했다. 이 계산 상태만으로 지반의 실제 붕괴 여부를 판정할 수 없다.

세 논문은 이 규모 링 기초의 실측 침하를 제공하지 않는다. 위 침하값은 이번 유한요소 모델의 수치 결과이며, 현장 침하와의 독립 비교는 별도 자료가 필요하다.

## 입력값 대조

- 로봇팔 논문 §7.7, PDF 27쪽은 천부 경화층 아래 모래의 밀도를 1,300–1,500 kg/m³로 추정한다. 공통 밀도 1,300 kg/m³는 이 범위의 하한이다. 이 절의 관측 범위는 천부 약 40 cm다.
- Murdoch 표 2, PDF 16쪽의 1 m 깊이 예측값은 밀도 1,307 kg/m³, Vs 94.6 m/s, Vp 157.9 m/s, 구속압 5.44 kPa다. Vs와 Vp를 유지하고 목표 밀도 1,300 kg/m³에서 환산하면 E 28.38780 MPa, ν 0.220046, G 11.63391 MPa, K 16.90026 MPa다. 입력한 E 28.4 MPa, ν 0.22, G 11.63934 MPa, K 16.90476 MPa와의 반올림 오차는 모두 0.05% 이내다.
- 로봇팔 논문의 E 11.48 MPa, G 4.47 MPa, K 7.79 MPa는 대표 밀도 1,211 kg/m³에서 상부 40 cm 전체의 파속으로 얻은 값이다. 이번의 1 m 예측 강성 28.4 MPa는 약 2.47배이며 관측 깊이와 상태가 다르다. 28.4 MPa의 선정 근거는 Murdoch 표 2다.
- Dotson 그림 5·6은 후보별 저응력 직접전단 최대 저항 포락선의 c와 φ를 제공한다. 선택한 MGS-1, MGS-1S, JEZ-1 점은 각각 약 1,320, 1,310, 1,380 kg/m³이고, 모두 1,300 kg/m³로 이전한 초기 입력이다. c와 φ는 그래프 판독 근사값이다. 원시 수치와 회귀 공분산을 확보한 완전 보정은 아니다.
- ψ = 0과 경화계수 0은 현재 모델의 가정이다. 세 논문에서 이 값들을 모두 실험으로 식별한 것은 아니다.

## 전단 응답 대조

Dotson의 수직응력 범위인 0.097–0.68 kPa 안에서 12개 단일 요소 단순전단 시험을 수행했다. 1 m 정육면체의 하부를 고정하고 상부의 수평면 변위를 묶었으며, 수직응력을 일정하게 유지해 공학전단변형률 0.2%까지 가했다. 각 적분점의 유한 응답, 수직응력 유지와 전단력 평형을 확인했다. 아래 비교 기준은 그래프 판독 c·φ로 계산한 τ = c + σn tanφ다. 단순전단 경계조건은 분할 전단상자의 국부 변형과 차이가 있다.

수직응력 0.5 kPa에서 비교한 결과다.

| 후보 | 그래프 기반 포락선 | OpenSees DP 전단 저항 | 차이 |
|---|---:|---:|---:|
| MGS-1 | 0.6366 kPa | 0.7421 kPa | +16.57% |
| MGS-1S | 0.6173 kPa | 0.7263 kPa | +17.66% |
| JEZ-1 | 0.4526 kPa | 0.5523 kPa | +22.03% |

네 수직응력 모두에서 같은 비율의 차이가 나타났다. 현재 삼축압축 대응 DP 변환은 해당 단순전단 경로에서 저항을 16.6–22.0% 크게 계산한다. 따라서 삼축압축의 해석해 일치와 직접전단 실험의 재현을 구분해야 한다. 강성 환산과 수치 평형은 확인됐고, 직접전단을 목표로 하는 강도 보정은 남아 있다. 이번 비교에서는 승인된 초기 물성값을 그대로 사용했다.

## 재현과 결과 파일

[물성 및 시편 비교 코드](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/verification/insight_sand_literature.py)

[지반 해석 코드](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_dome_footing.py)

[전단 곡선과 비교 수치](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/shear_comparison.csv)

[입력값과 원본 해시](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/summary.json)

[지반 실행 상태와 침하](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/ground_comparison.json)

![전단 저항 비교](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/shear_comparison.png)

지반 실행 폴더에는 summary.json, console.log, opensees.log, solver_attempts.jsonl과 solver_checkpoint.json을 남겼다. 완주한 후보 폴더에는 footing_settlement.html과 물리 변위 payload도 있다.

## 참고 자료

1. [Dotson 외 원문](/Users/sydney/Downloads/2403.11029v1.pdf), 그림 5·6 및 §2.4. 최대 저항 포락선, 저응력 직접전단.
2. [InSight 로봇팔 연구 원문](/Users/sydney/Downloads/s11214-023-00964-0.pdf), §7.7 및 PDF 27쪽. 천부 밀도·파속·탄성계수.
3. [Murdoch 외 원문](/Users/sydney/Downloads/InSightLanderFlexibleModes_Changes.pdf), 식 4–6 및 표 2, PDF 15–16쪽. 깊이별 예측 물성.
4. [OpenSees DP 문서](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html), 항복면과 응력경로별 변환. 프로젝트는 기존 capmat의 삼축압축 대응식을 사용한다.

## 시간 제한 해제 후 중지 및 원인 분석

사용자 요청으로 추가 계산을 중지했다. JEZ-1 증분 축소 실행은 약 21.3분 동안 자중 49.20%까지 도달했고 기초 재하는 시작하지 않았다. 이전 600초 실행과 이번 사용자 중지를 구분한다. 모든 관련 해석 프로세스는 종료했다.

현재 모래층 기본 후보는 MGS-1, 대체 후보는 MGS-1S를 추천하며, JEZ-1 초기 조합은 기본 후보 채택을 보류한다. 실제 붕괴를 판정한 결과는 아니다.

[시간 지연 원인 및 재료 선정 판단](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/results/insight_sand_literature_20261118/jez1_diagnosis_KO.md)
