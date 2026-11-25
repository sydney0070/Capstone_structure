# 철근 Steel02 재보정 — 2026-10-12

상태: **REUSED_HOLDOUT_SCREENING_PASSED_NEW_VALIDATION_REQUIRED**

기존 원자료·4학습/2검증 분할·탐색범위·합격 기준을 유지했다. 학습 시편 A65의 극값 오차도 기준을 넘었으므로, 학습 곡선과 인장/압축 극값을 함께 고려했다.
목적함수는 시편별 0.8×곡선 NRMSE² + 0.1×인장 극값 오차² + 0.1×압축 극값 오차²의 평균이다. 5개 시작점 중 학습 목적함수로만 선택하고 선택값 저장 후 검증을 실행했다.
**기존 검증 시편의 결과는 이전 작업에서 이미 관측했다. 이번 결과를 새 독립 블라인드 검증으로 부르지 않으며, 검증 결과를 보고 재시도하지 않았다.**

| 시편 | 역할 | 이전 NRMSE | 이번 NRMSE | 이전 극값 오차 | 이번 극값 오차 | 기준 통과 |
|---|---|---:|---:|---:|---:|---|
| 01_A64.csv | train | 2.519% | 3.321% | 9.217% | 2.191% | True |
| 02_A65.csv | train | 2.067% | 2.380% | 12.246% | 6.428% | True |
| 03_A66.csv | train | 2.840% | 3.054% | 4.658% | 4.738% | True |
| 04_A67.csv | train | 2.040% | 2.312% | 7.172% | 2.773% | True |
| 05_A68.csv | validation | 2.603% | 3.207% | 8.239% | 4.076% | True |
| 06_A69.csv | validation | 2.373% | 2.530% | 11.084% | 5.320% | True |

| 계수 | 값 |
|---|---:|
| E | 218990793.2 |
| Fy | 512792.1296 |
| b | 0.01288371987 |
| cR1 | 0.9266119551 |
| cR2 | 0.1430886481 |
| a1 | 0.06331854194 |
| a3 | 0.04637773159 |

E·Fy는 kPa, 나머지는 무차원. R0·a2·a4 등 고정 계수는 reviewed_config.json에 포함된다.
경계 부근 계수: 없음. 정규화 국부 감도 조건수: 35.09 (유일해·신뢰구간을 의미하지 않음).
최대 증분 세분화 차이: 0.000000%. 최대 반력–응력 상대 차이: 0.
검토 기준은 기존과 동일: 검증 시편별 NRMSE ≤5%, 극값 오차 ≤10%, 증분 차이 ≤1%. 이번에는 반력–응력 일치도 확인했다.
적용 범위: UCSD Grade 60 A6 그룹의 제공된 반복 단축 이력. 전체 Grade 60 제품, 좌굴·파단·피로·다축 거동 검증이 아니다.
기존 fitted/ 및 reviewed_config.json은 보존했다. 새 후보는 이 폴더의 reviewed_config.json이며 추가 독립 검증 전에는 최종 채택값으로 취급하지 않는다.

재현: 패키지 루트에서 `python material_calibration/refine_rebar.py --out material_calibration/rebar_ucsd_a6/새폴더명`.
실험 출처: [NHERI SimCenter qfem-0028](https://nheri-simcenter.github.io/quoFEM-Documentation/common/user_manual/examples/desktop/qfem-0028/README.html). 원본 해시·라이선스는 상위 provenance.json과 raw/rebar_ucsd/SOURCE_LICENSE.txt를 따른다.
