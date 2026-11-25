# aluminum_6014_t4 보정 결과

상태: **FITTED_BUT_NOT_ACCEPTED**

Material: `J2Plasticity`. Element: `ShellMITC4 + PlateFiber`.

AW6014 T4, 23 C, rolling 0 degrees, 15 mm/min; engineering strain <=2%; axial monotonic response only. Elastic strain/compliance review required. nu=.33 assumed, H=0 fixed; no multiaxial/isotropy/cyclic validation.

| 변수 | 보정값 |
|---|---:|
| E | 60000000 |
| sig0 | 114684.254 |
| sigInf | 300000 |
| delta | 16.4331785 |

응력·탄성계수 단위 kPa. 그 외 단위는 starting_config.json을 따른다.

| 시편 | 역할 | NRMSE | 극값 오차 | 증분 오차 | 기준 통과 |
|---|---|---:|---:|---:|---|
| specimen_1.csv | train | 4.349% | 2.352% | 0.00000% | True |
| specimen_2.csv | train | 3.236% | 2.894% | 0.00000% | True |
| specimen_3.csv | train | 8.613% | 2.802% | 0.00000% | False |
| specimen_4.csv | validation | 4.155% | 4.148% | 0.00000% | True |
| specimen_5.csv | validation | 5.158% | 2.830% | 0.00000% | False |

경계 부근 매개변수: ['E', 'sigInf'].
추가 확인사항: ELASTIC_STRAIN_AND_MACHINE_COMPLIANCE_NOT_CONFIRMED.

사전 기준: 검증 시편 각각 NRMSE ≤5%, 극값 오차 ≤10%, 내부 증분 세분화 오차 ≤1%. 설계 규격이 아닌 프로젝트 검토 기준이다. 평균 곡선 오차 통과만으로 파괴·피로·다축·환경 거동의 검증을 주장하지 않는다.

실험 출처·처리: targets.json, provenance.json. 곡선: fitted/train_*, fitted/validation_*의 response.png. 초기값·고정값·탐색범위: starting_config.json.
