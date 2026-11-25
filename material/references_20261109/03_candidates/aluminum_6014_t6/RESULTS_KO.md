# aluminum_6014_t6 보정 결과

상태: **FITTED_REQUIRES_PARAMETER_OR_MEASUREMENT_REVIEW**

Material: `J2Plasticity`. Element: `ShellMITC4 + PlateFiber`.

AW6014 T6, 23 C, rolling 0 degrees, 15 mm/min; engineering strain <=2%; axial monotonic response only. Elastic strain/compliance review required. nu=.33 assumed, H=0 fixed; no multiaxial/isotropy/cyclic validation.

| 변수 | 보정값 |
|---|---:|
| E | 60000000 |
| sig0 | 171541.853 |
| sigInf | 300000 |
| delta | 23.8772626 |

응력·탄성계수 단위 kPa. 그 외 단위는 starting_config.json을 따른다.

| 시편 | 역할 | NRMSE | 극값 오차 | 증분 오차 | 기준 통과 |
|---|---|---:|---:|---:|---|
| specimen_1.csv | train | 2.554% | 2.016% | 0.00000% | True |
| specimen_2.csv | train | 2.568% | 2.234% | 0.00000% | True |
| specimen_3.csv | validation | 2.449% | 1.330% | 0.00000% | True |

경계 부근 매개변수: ['E', 'sigInf'].
추가 확인사항: ELASTIC_STRAIN_AND_MACHINE_COMPLIANCE_NOT_CONFIRMED.

사전 기준: 검증 시편 각각 NRMSE ≤5%, 극값 오차 ≤10%, 내부 증분 세분화 오차 ≤1%. 설계 규격이 아닌 프로젝트 검토 기준이다. 평균 곡선 오차 통과만으로 파괴·피로·다축·환경 거동의 검증을 주장하지 않는다.

실험 출처·처리: targets.json, provenance.json. 곡선: fitted/train_*, fitted/validation_*의 response.png. 초기값·고정값·탐색범위: starting_config.json.
