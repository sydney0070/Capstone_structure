# rebar_ucsd_a6 보정 결과

상태: **FITTED_BUT_NOT_ACCEPTED**

Material: `Steel02`. Element: `truss`.

UCSD Grade 60 A6 group supplied histories; ASTM subclass mapping not asserted; not buckling, fracture or all Grade 60 products.

| 변수 | 보정값 |
|---|---:|
| E | 215963757 |
| Fy | 512903.843 |
| b | 0.0221485218 |
| cR1 | 0.912542236 |
| cR2 | 0.0752671623 |
| a1 | 0.0470995896 |
| a3 | 0.0346058577 |

응력·탄성계수 단위 kPa. 그 외 단위는 starting_config.json을 따른다.

| 시편 | 역할 | NRMSE | 극값 오차 | 증분 오차 | 기준 통과 |
|---|---|---:|---:|---:|---|
| 01_A64.csv | train | 2.519% | 9.217% | 0.00000% | True |
| 02_A65.csv | train | 2.067% | 12.246% | 0.00000% | False |
| 03_A66.csv | train | 2.840% | 4.658% | 0.00000% | True |
| 04_A67.csv | train | 2.040% | 7.172% | 0.00000% | True |
| 05_A68.csv | validation | 2.603% | 8.239% | 0.00000% | True |
| 06_A69.csv | validation | 2.373% | 11.084% | 0.00000% | False |

경계 부근 매개변수: 없음.
추가 확인사항: 선택 시편 범위 밖 거동은 검증하지 않음.

사전 기준: 검증 시편 각각 NRMSE ≤5%, 극값 오차 ≤10%, 내부 증분 세분화 오차 ≤1%. 설계 규격이 아닌 프로젝트 검토 기준이다. 평균 곡선 오차 통과만으로 파괴·피로·다축·환경 거동의 검증을 주장하지 않는다.

실험 출처·처리: targets.json, provenance.json. 곡선: fitted/train_*, fitted/validation_*의 response.png. 초기값·고정값·탐색범위: starting_config.json.
