# PDMY02 재료 평가: 문헌 근거와 코드 판정

조사일: 2026-09-17. 대상: `material_proto.py`의 3D SSPbrick, 일정 수직응력, 변형률제어 반복단순전단.

## 먼저 구분할 것

문헌을 조사한 결과, 모든 흙에 통용되는 **강성·감쇠·강도·체적변형·반복 변화의 고정 합격값**은 확인되지 않았다. 서로 다른 흙/시험의 관측값을 이 재료의 합격 기준으로 옮겨 쓰지 않는다. 현재 입력값은 PDMY02 권장 예제에 가깝고, 대상 흙의 실험자료는 주어지지 않았다.

따라서 출력은 다음 세 가지를 구분한다.

- **물리/모델 일관성**: 양의 할선강성, 강도면 범위, 닫힌 루프의 비음수 에너지.
- **수치 검증**: 해석 완료, 반력 평형, 변형률 일치, 독립 미소변형 시험, 증분 세분화.
- **대상 흙 적합성**: 같은 흙·초기 상태·경계조건·진폭·사이클의 독립 비교자료와 허용구간을 지정해야 PASS/FAIL 판정. 없으면 NOT_EVALUATED.

PASS는 표시된 검사만 통과했다는 뜻이다. 전체 지반 모델 사용 승인이나 안전 판정이 아니다. `INFO`는 참고값이고 `NOT_APPLICABLE`은 현 시험으로 판정할 수 없는 항목이다. 종료코드 0도 자료 부족을 포함할 수 있으며, 실패 검사가 있으면 2이다.

## 여섯 항목의 문헌 조사 및 판정값

| 항목 | 관련 문헌에서 확인한 값/관계 | 코드에 적용하는 판정 | 한계 |
|---|---|---|---|
| 강성 | PDMY02 권장표의 상대밀도 30~75% 모래 예시는 기준 전단강성 60,000~130,000 kPa. 작은 변형 강성과 압력 의존성을 맞추고 MRD 곡선을 검증해야 함 [1–3] | 모든 완결 사이클 Gsec > 0. 별도 미소변형 시험에서 G0 측정. G0, 마지막 Gsec, Gsec/G0를 실험 lower~upper와 비교 | 권장표 범위는 합격 구간이 아님. 표의 강성 기준압력 80 kPa와 입력 refPress 101 kPa 표기 차이에 주의. 현재 p0는 약 63.6 kPa |
| 감쇠 | 닫힌 반복 루프의 D=ΔW/(4πWs), Ws=Gsec·γa²/2. 강성감소와 감쇠곡선을 함께 맞춰야 함. 일반 MRD 관계는 약 0.3~0.5% 초과에서 신뢰 제한 [1,2] | 첫 사이클 제외. 응력 폐합/체적변형 드리프트 검사 후 D 산정. 닫힌 루프 D≥0. 실험 D 허용구간과 비교 | 5%, 10%, 20% 같은 고정 합격값 없음. 현재 1% 시험에 소변형 경험곡선을 외삽하여 합격시키지 않음 |
| 강도 | PDMY02 출력 정의상 0≤eta_r≤1. 권장 마찰각 31~36.5°는 상대밀도별 예시 [3]. 큰 변형 지반응답에는 강도 제약이 필요 [1,2] | eta_r 범위 검사. 반복시험 관측 최대응력을 같은 시험의 실험값과 비교 | eta_r는 입력 재료의 자체 강도면에 대한 비이므로 독립 검증이 아님. 관측 최대 τ를 극한강도로 부르지 않음. τ/σv=tanφ로 단정하지 않음 |
| 체적변형 | Hsu & Vucetic(2004): 모래의 체적변형 문턱 전단변형률 γtv≈0.01~0.02%; PI≈30 점토는 약 0.04~0.09% [4] | 전단 시작 상태를 뺀 Δεv의 마지막 값·최솟값·최댓값 기록. 같은 시험의 Δεv 허용구간과 비교 | 문턱 **전단변형률**은 허용 **체적변형률**이 아님. 현재 1%는 이 문턱보다 크지만, 이것만으로 수축량 적합성을 판단할 수 없음 |
| 반복 변화 | Carotenuto 등(2022): 조밀한 규사에서 선행 반복 이력에 따라 이후 비배수 30회 이내 Gcy가 8~45 MPa 또는 76~82 MPa, 파괴회수는 60 또는 245로 달라짐 [5] | G_N/G_1, τa,N/τa,1, 사이클별 Δεv 기록. 두 비율을 같은 조건 시험의 허용구간과 비교 | 이 문헌의 삼축/선행재하/비배수 관측값은 현재 배수 단순전단의 합격값이 아님. 변화율 0%, 강성비 1이 보편 목표가 아님. 3회 시험으로 장기 수명 판정 불가 |
| 수치 안정성 | 문헌의 재료 보정/검증 요구와 별개로, 이 프로토타입에 보편적인 오차 합격 숫자는 제시되지 않음 [1,2] | 아래 **프로젝트 수치 허용오차**로 직접 평가 | 아래 2%, 5% 등을 논문/설계기준에서 정한 합격값으로 인용하면 안 됨 |

## 수치 허용오차: 구현자가 명시적으로 선택한 개발 기준

`material_evaluation.py`의 `TOL`에 모았다. 대상 해석의 정확도 요구에 따라 조정할 수 있다.

| 검사 | 기준 |
|---|---|
| 해석 완료 | 요청 증분 전부 수렴, NaN/Inf 없음 |
| 반력 평형 | max\|τreaction−τmaterial\| / max(σv, max\|τ\|, 1 kPa) ≤ 1e−4 |
| 일정 수직응력 | max\|σzz+σv\| / σv ≤ 1e−4 |
| 절점–재료 및 목표–응답 변형률 | 절대차 ≤ 1e−8 (변형률 fraction) |
| PDMY02 eta_r | [0,1]에 수치 여유 ±1e−6 |
| 미소변형 G0 | 독립 시편 γa=1e−7 및 5e−8, 두 G0 모두 양수, 상대차 ≤2% |
| 40→80→160 분기당 증분 | 각 단계에서 응력경로 정규화 차이, 최대응력·마지막 Gsec 상대차 ≤2% |
| 사이클 전단일 | 모든 완결 사이클에서 세분화 전후 상대차 ≤5%; 0 부근 분모 하한 1e−8 kPa |
| 체적변형 경로 | 세분화 전후 절대차 ≤1e−5 = 0.001%p |
| 감쇠 계산용 루프 판별 | 첫 사이클 제외, 끝–시작 응력 차이 / (2τa) ≤1%, 체적변형 드리프트 ≤1e−5 |

감쇠용 체적변형 드리프트 조건은 정상상태 감쇠로 오인하지 않기 위한 **보수적인 프로젝트 선택**이다. 문헌의 보편적 감쇠 측정 표준이 아니다. 조건을 못 맞춘 경우 에너지 부호를 절댓값으로 바꾸거나 열린 루프를 강제로 닫지 않는다. `shear_work_kPa`는 그대로 저장하되 소산에너지라고 단정하지 않는다. 순환하는 응력–변형률 루프만으로 내부 상태 전체의 주기성을 증명하지도 못한다.

## 산정 정의

- `G0_kPa`: 초기 압밀을 새로 수행한 미소변형 시험의 첫 양의 반전점까지의 할선강성. 두 진폭에서 수렴 확인 후 작은 진폭 값을 사용한다. `G_ref`를 분모에 그대로 대입하지 않는다.
- `Gsec_last_kPa`: 마지막 완결 사이클의 양·음 **변형률 반전점에서** `(τ+−τ−)/(γ+−γ−)`. 응력의 별도 최댓값/최솟값을 조합하지 않는다.
- `G_ratio_last`: 마지막 Gsec/G0. 같은 초기 압력의 독립 미소변형 시편을 기준으로 한다.
- `D_last_percent`: 마지막 완결 사이클의 감쇠비(%). 평가 부적합이면 JSON null.
- `tau_peak_observed_kPa`: 전체 성공 증분 중 max|τ|. 극한강도나 단조전단 최대강도가 아니다.
- `eps_v_final_percent`: 100×[tr(εfinal)−tr(εinitial)]. 음수=수축, 양수=팽창. 압밀 변형 제외.
- `G_last_over_first`, `tau_last_over_first`: 마지막/첫 완결 사이클 비율. 일정 변형률 진폭에서는 두 비율이 사실상 같으므로 독립적인 두 검증으로 과대평가하지 않는다.
- 현재 static 해석의 OpenSees time은 하중계수이며 물리적인 재하 주파수로 해석하지 않는다.
- SSPbrick에는 간극수압 자유도가 없으므로 ru, CRR, 액상화 파괴회수는 판정하지 않는다.

## 사용법

기존처럼 실행하면 원래 그래프/애니메이션과 함께 평가를 수행한다. GUI 없는 실행:

```bash
python material_proto.py --no-plots
```

자동으로 원 시험 + 독립 미소변형 2개 + 증분 2배/4배 시험을 수행한다. 재료 상수와 원 시험의 경계조건은 유지했다.

```bash
python material_proto.py --no-plots --quick
python material_proto.py --no-plots --gamma 0.001 --steps 40 --cycles 10 --output-dir results_01pct
python -m unittest -v test_material_evaluation.py
```

`--quick`에서는 G0와 증분 수렴 판정을 생략하므로 NOT_EVALUATED가 정상이다. 각 실행은 `--output-dir`의 같은 이름 결과를 갱신하므로 비교할 실행은 서로 다른 폴더를 사용한다.

생성 파일:

- `evaluation.txt`: 한국어 판정·측정값·허용값·근거.
- `evaluation.json`: 판정, 설정, 버전, 출처, 수치 허용오차, 타깃 설정.
- `cycles.csv`: 사이클별 강성·전단일·감쇠·폐합 오차·체적변형.
- `base_history.npz`, `probe_history.npz`, `probe_half_history.npz`, `fine_history.npz`, `finer_history.npz`: 응력/변형률 전체 성분, 반력, 목표변형률, 절점변위.

NPZ는 실행 시 생성된 경우만 읽어야 한다. `--quick` 실행 후 같은 폴더에 남아 있는 이전 보조시험 파일은 이번 실행 결과가 아니다. 이번 실행 보조시험 목록은 evaluation.json의 auxiliary_runs에 명시된다.

## 대상 흙 자료 입력

`material_targets.json`에서 `context.conditions_confirmed`를 true로 설정하기 전에 흙의 종류, 밀도/간극비, 초기 응력 상태, 배수 조건, 시험 이력, 진폭·회수를 확인한다. `soil_description`에 상태와 자료 식별자를 기록한다. 코드가 자동으로 확인하는 것은 시험 종류·수직응력·진폭·사이클 일치이며, 흙의 동일성은 입력자의 확인에 의존한다.

각 null을 다음 구조로 교체한다.

```json
{"lower": 4.0, "upper": 6.0, "source": "형식 설명용 예시: 실제 시험 보고서/시편/사이클/측정 산포로 교체"}
```

위 숫자는 **D_last_percent의 파일 형식 설명용 가상 예시**이며 권장 감쇠 범위가 아니다. 실제 데이터가 없으므로 기본 파일에는 넣지 않았다. 허용구간은 독립 시험 반복의 산포나 프로젝트의 명시적 정확도 요구로 정한다. 모델 계산값을 보고 합격하도록 범위를 역설정하지 않는다. 모든 지표를 제공해도 판정은 그 압력·진폭의 시험에 한정된다.

## 이번 검증 결과

임시 가상환경의 OpenSeesPy/OpenSees 3.8.0에서 기본값(100 kPa, ±1%, 3회, 분기당 40/80/160 증분)으로 실행했다.

- 수치 안정성: PASS.
- G0 ≈79,389.96 kPa, 마지막 Gsec ≈7,439.37 kPa, Gsec/G0 ≈0.09371.
- 관측 최대 |τ| ≈75.11199 kPa, max eta_r ≈0.95052.
- 최종 전단 유발 Δεv ≈−0.19578% (수축).
- 마지막/첫 사이클 G 비율 ≈1.00055.
- 정상 루프 감쇠: 체적변형 드리프트 조건 미충족으로 NOT_EVALUATED.
- 실험자료 적합성: 독립 비교자료 미입력으로 NOT_EVALUATED.

이 결과는 새 문헌의 보편적 재료 합격범위를 뜻하지 않는다. 구현 검증에는 알려진 면적을 갖는 타원 루프, 탄성 루프, 음의 에너지, 열린 루프, 중도 실패, 누락/불일치 타깃 등을 포함한 단위시험도 사용했다.

## 참고문헌과 확인 범위

1. Stewart, J. P., Afshari, K., & Hashash, Y. M. A. (2014). *Guidelines for Performing Hazard-Consistent One-Dimensional Ground Response Analysis for Ground Motion Prediction*, PEER 2014-16. 발행기관 공개 초록의 권고 1–3 확인: 실측 Vs, MRD의 약 0.3~0.5% 한계, 큰 변형에서 강도 고려. [발행기관](https://peer.berkeley.edu/publications/2014-16).
2. Stewart, J. P., Kwok, A. O.-L., Hashash, Y. M. A., Matasovic, N., Pyke, R., Wang, Z., & Yang, Z. (2008). *Benchmarking of Nonlinear Geotechnical Ground Response Analysis Procedures*, PEER 2008-04. 발행기관 공개 초록: 작은 변형 <0.3%, 약 1%의 큰 변형 강도, backbone/감쇠 동시 보정의 절충. [발행기관](https://peer.berkeley.edu/publications/2008-04).
3. OpenSees, *Pressure Dependent Multi Yield 02*, 공식 재료 문서. 출력 성분·eta_r 정의, 기준 압력, 권장 파라미터표 확인. 논문 기반 대상지반 시험값과 구분. [문서](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html).
4. Hsu, C.-C., & Vucetic, M. (2004). *Volumetric Threshold Shear Strain for Cyclic Settlement*, Journal of Geotechnical and Geoenvironmental Engineering 130(1), 58–70. 출판사 검색에 제공된 공개 초록 확인; 유료 본문의 모든 표를 검토한 것은 아님. 7종 흙, 11개 다단계 직접단순전단시험의 문턱 범위. [DOI](https://doi.org/10.1061/(ASCE)1090-0241(2004)130:1(58)).
5. Carotenuto, P., Jostad, H. P., & Teng, Y. (2022). *Effect of drained cyclic preshearing on the undrained cyclic response of dense sand in triaxial tests*, Proceedings of the 20th ICSMGE, pp.53–58. 공개 논문 초록, 시험 조건, 결과 부분 확인: 선행 이력과 반복응답의 차이. [ISSMGE 원문](https://www.issmge.org/uploads/publications/1/120/ICSMGE_2022-10.pdf).

감쇠의 면적 정의 참고: Stewart 등(2001), PEER 2001-09, 지반 동적 물성 부분. [PEER 원문](https://peer.berkeley.edu/sites/default/files/0109_j._stewart_s._chiou_j._bray_r._graves_p._somerville_n._abrahamson.pdf).

Darendeli(2001) 전체 경험식이나 Bolton(1986)의 광물별 마찰각을 자동 합격판정에 구현하지 않았다. 전자는 원문 접근과 필요한 토질 입력 확인이 부족하고, 후자는 현재 DSS 최대응력의 직접 합격값으로 사용할 수 없기 때문이다.
