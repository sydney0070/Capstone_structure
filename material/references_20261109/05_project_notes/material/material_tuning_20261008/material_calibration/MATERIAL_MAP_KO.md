# 재료별 OpenSees material·element와 보정 상태

2026-10-19 추가 자료 검토: [달·화성 지반 원문 검토](regolith_sources_20261019/REVIEW_KO.md). JSC-1A 삼축시험 요약 46행 추출, JEZ-1 접촉시험 응력 범위 및 압축·제하 조건 보완. 새로운 계수 보정·독립 검증 완료를 뜻하지 않는다.

2026-10-12 재료 전용 후속 갱신. `material`은 응력–변형률 또는 힘–변위 법칙이고, `element`는 그 법칙을 사용하는 시편/구조 요소다. 아래 element는 현재 실행 코드 기준이다. 전체 돔 형상 없이도 시편 보정이 가능하다. 최신 작업 요약은 [재료 전용 진행 결과](PROGRESS_20261012_KO.md)를 따른다.

| 실제 재료/거동 | OpenSees material | 현재 시편 element | 보정 대상·제한 |
|---|---|---|---|
| S355 구조용 강재 | uniaxialMaterial `Steel02` | `truss` | 기존 선택 시편 검증 통과. 다축 쉘 검증 아님 |
| Grade 60 철근 UCSD A6 그룹 | uniaxialMaterial `Steel02` | `truss` | 7변수 재보정으로 기존 검증 2시편의 수치 기준 통과. 검증 자료 재사용으로 새 독립 검증 필요 |
| 콘크리트 피복 | uniaxialMaterial `Concrete02` | `truss` | 압축강도·최대강도 변형률·연화. 인장·반복 자료 별도 필요 |
| 구속 콘크리트 코어 | uniaxialMaterial `Concrete02` | `truss` | 구속 조건 시험 필요. 피복 시험값으로 완료 판정 불가 |
| 강재 쉘 | nDMaterial `J2Plasticity` → `PlateFiber` | `ShellMITC4` | 다축 등방 탄소성. Steel02 보정값을 직접 복사할 수 없음 |
| 알루미늄 AW6014 T4 / T6 | nDMaterial `J2Plasticity` → `PlateFiber` | `ShellMITC4` | 이번 단조 인장 2% 이내 보정. E·포화강도 식별/측정 검토 필요 |
| CFRP | nDMaterial `ElasticOrthotropic` → `PlateFiber` | `ShellMITC4` | 정렬된 재료축 탄성. Ex·Ey·전단계수·포아송비, 방향별 시험 필요. 손상·적층 파괴 미구현 |
| 콘크리트 쉘 | nDMaterial `ElasticIsotropic` → `PlateFiber` | `ShellMITC4` | 현재 E·ν 탄성 비교만 가능. 균열·압괴 모델 아님 |
| 소결 레골리스 | nDMaterial `ElasticIsotropic` → `PlateFiber` | `ShellMITC4` | 제조조건별 탄성 출발 모델. 취성강도·균열 미구현 |
| 레골리스 지오폴리머 | nDMaterial `ElasticIsotropic` → `PlateFiber` | `ShellMITC4` | 배합·양생별 탄성 출발 모델. 연화·환경 의존성 미구현 |
| 황콘크리트 | nDMaterial `ElasticIsotropic` → `PlateFiber` | `ShellMITC4` | 배합·온도별 탄성 출발 모델. 열·파괴 거동 미구현 |
| 지구 지반 탄성 기준 | nDMaterial `ElasticIsotropic` | `SSPbrick` | E·ν·밀도. 항복/액상화 재현 불가 |
| 지구 모래 | nDMaterial `PressureDependMultiYield02` | `SSPbrick` | G·K·마찰·수축·팽창. 현재 단순전단 시편 |
| 지구 점토 | nDMaterial `PressureIndependMultiYield` | `SSPbrick` | G·K·비배수 강도 중심. 압밀·OCR의 일반 모델 아님 |
| 지구 포화 모래 대안 | nDMaterial `PM4Sand` | `SSPquad` | Dr·G0·hpo 등. 현재 일정체적 근사이며 u-p 유체 모델 아님 |
| 지구 실트 대안 | nDMaterial `PM4Silt` | `SSPquad` | Su·G0·hpo 등. 실트 자료를 일반 점토 자료로 혼용하지 않음 |
| 달 건조 레골리스 | nDMaterial `DruckerPrager` | `SSPbrick` | c·φ·K·G·팽창·경화. 특정 시뮬런트/밀도별 보정 필요 |
| 화성 건조 레골리스 | nDMaterial `DruckerPrager` | `SSPbrick` | JEZ-1 문헌 강도값 추가. 곡선 보정·독립 검증 미완료 |
| SSI 수직 접촉·간극 | uniaxialMaterial `ElasticPPGap` | `zeroLength` | 압축 접촉 스프링. 초기 간극·압축강성·항복력 |
| SSI 마찰 접촉 | 별도 material 없음 | `zeroLengthContactASDimplex` | element에 Kn·Kt·μ 직접 입력. 현재 부착력 없는 마찰 모델 |
| 말뚝 횡방향 지반 반력 | uniaxialMaterial `PySimple1` | `zeroLength` | pult·y50·gap/drag. 지반+말뚝 치수 의존 |
| 말뚝 주면 마찰 | uniaxialMaterial `TzSimple1` | `zeroLength` | tult·z50. 둘레·분담 길이 변환 필요 |
| 말뚝 선단 지지 | uniaxialMaterial `QzSimple1` | `zeroLength` | qult·z50·suction. 선단 면적 변환 필요 |
| SSI 방사감쇠 근사 | uniaxialMaterial `Viscous` | 직접 material 속도 시험 | C·α. 일반 재료강도 보정과 구분, 전체 모델에서는 zeroLength로 배치 가능 |

## 이번 실행 결과

- [철근 후속 재보정](rebar_ucsd_a6/refinement_20261012/RESULTS_KO.md): 4학습·기존 2검증 유지. 학습 곡선과 인장/압축 극값의 가중 목적함수로 재보정. 검증 NRMSE 3.207%/2.530%, 극값 오차 4.076%/5.320%로 기존 수치 기준 통과. 평균 오차는 이전보다 다소 증가했으나 극값 오차가 개선됨. 기존 결과를 이미 확인한 검증 시편이므로 `REUSED_HOLDOUT_SCREENING_PASSED_NEW_VALIDATION_REQUIRED`이며 새 블라인드 검증으로 주장하지 않는다. [이전 결과](rebar_ucsd_a6/RESULTS_KO.md)와 기존 fitted/ 설정은 보존했다.
- [AW6014 T4](aluminum_6014_t4/RESULTS_KO.md): 3학습·2검증. 검증 NRMSE 4.155%/5.158%, 2번째 곡선 기준 미달.
- [AW6014 T6](aluminum_6014_t6/RESULTS_KO.md): 2학습·1검증. NRMSE 2.449%, 극값 오차 1.330%. 수치 기준은 통과했지만 E·sigInf가 탐색 경계에 도달하여 완료 보류.
- T4/T6 모두 초기 탄성 구간과 시험기 변형/변형률 정의 확인이 필요하다. 현 자료는 초기 구간 점수가 적고 시편 간 차이도 있다. 오차를 줄이려고 탄성계수를 근거 없이 낮추거나 검증 시편에 맞춰 재보정하지 않았다. 실제 재료의 E·포화강도로 채택하지 않는다.
- [알루미늄 후속 식별성 점검](aluminum_audit_20261012/RESULTS_KO.md): 학습자료만으로 E=60/70/80 GPa 고정 재보정 6건 수행. 모두 sigInf가 기존 탐색 하한에 도달. 기존 보정점에서 sigInf–delta 감도 방향 코사인이 T4 0.999576, T6 0.999182로 유사하여 2% 범위에서 두 계수의 분리 식별이 어렵다. 10–50 MPa 구간의 학습 측정점은 시편당 1~2개. 기존 보류 판정을 유지하며 E 범위 확대·임의 변형률 보정은 하지 않았다.
- [JEZ-1 강도 입력](mars_jez1_strength/README.md): 논문의 c·φ·접촉마찰각을 기록했다. 강성·팽창·연화는 아직 식별하지 않았다.

## 나머지 재료의 다음 작업

콘크리트 micro-concrete 원자료는 헤더가 없어 단위 확인을 계속 필요로 한다. 저자 공개본 링크는 확인했지만 원문 서버 접근이 제한되었다. 별도 검증 없이 일반 구조 콘크리트로 대체하지 않는다.

지구 모래의 새 후보는 [LEAP-2020 GWU 삼축·직접단순전단 데이터](https://doi.org/10.17603/ds2-dyph-pf84)다. 실제 파일·배수조건·압밀응력·밀도 확인과 단순전단/삼축 시편 구분이 선행되어야 한다. 이 목록 추가는 다운로드/보정 완료가 아니다. 점토·실트는 별도 같은 상태의 반복시험 자료가 필요하다.

달은 기존 레골리스 데이터베이스의 물성 범위를 유지한다. 서로 다른 시뮬런트나 밀도 행을 평균하여 가상의 하나의 재료로 만들지 않는다. CFRP의 추가 자료 후보는 [Smartfan 횡방향 인장 데이터](https://doi.org/10.5281/zenodo.5876612)이나, 아직 원자료의 섬유/수지·방향 정의를 검토하지 않았고 횡방향 하나로 전체 직교이방성 계수를 식별할 수 없다.

소결체·지오폴리머·황콘크리트는 배합·제조 조건과 탄성 이후 거동을 확인한 뒤 비선형 모델을 선택한다. 기존 ElasticIsotropic 설정을 강도까지 보정한 모델로 부르지 않는다.

## 재현

패키지 루트에서 기존 requirements.txt 환경을 사용한다.

이번 후속 작업은 기존 결과를 덮어쓰지 않는 별도 폴더에서 수행한다. 아래 `새폴더명`은 실제로 존재하지 않는 경로를 지정해야 한다.

```bash
python material_calibration/refine_rebar.py --out material_calibration/rebar_ucsd_a6/새폴더명
python material_calibration/audit_aluminum.py --out material_calibration/새진단폴더명
```

아래 명령은 **이전 보정 절차**다. 이전 결과를 재현하며 이번 후속 후보를 생성하지 않는다.

```bash
python material_calibration/prepare_rebar.py
python tune.py fit --config material_calibration/rebar_ucsd_a6/starting_config.json --targets material_calibration/rebar_ucsd_a6/targets.json --maxiter 200 --out material_calibration/rebar_ucsd_a6/fitted
python material_calibration/assess_case.py rebar_ucsd_a6
```

알루미늄 원본 재추출에는 numpy·openpyxl이 있는 Python으로 `material_calibration/prepare_aluminum.py`를 실행한다. 이미 변환한 CSV를 사용하면 추가 Excel 패키지 없이 보정 가능하다. 각 aluminum_6014_t4/t6 폴더의 starting_config.json과 targets.json으로 같은 fit 명령을 실행하고 `--maxiter 120`을 사용한다. 평가에는 `assess_case.py aluminum_6014_t4` 또는 `aluminum_6014_t6`를 실행한다. 동일 경로 실행은 결과를 갱신한다.

## 출처

철근: [SimCenter qfem-0028](https://nheri-simcenter.github.io/quoFEM-Documentation/common/user_manual/examples/desktop/qfem-0028/README.html), [Carreño 등 (2020)](https://doi.org/10.1061/(ASCE)ST.1943-541X.0002505). 공식 예제와 같은 Steel02 계수 집합을 선택했으나 여기서는 결정론적 다중 시편 보정을 수행했다. 베이지안 불확실성 분석을 수행한 것은 아니다. 실제 저장소 commit과 파일 해시는 provenance.json에 기록했다.

알루미늄: Böhnke 등, [데이터 v2](https://doi.org/10.5281/zenodo.21993433). 원본 및 변환 데이터의 라이선스는 CC BY-SA 4.0. 시험 조건은 23°C, 압연 0°, 15 mm/min이다. 원본 source.xlsx를 수정하지 않았고 단위 변환·초기 오프셋·2% 절단을 provenance.json에 기록했다. 변환 데이터도 CC BY-SA 4.0으로 제공한다.

JEZ-1: Yin 등, [원 연구](https://doi.org/10.3390/app11157052), CC BY 4.0. [접촉 요소 공식 문서](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/zeroLengthContactASDimplex.html)와 비교해 현재 구현의 부착력 누락을 확인했다. 그 밖의 모델 근거와 기존 연구는 상위 SOURCES.md·RESEARCH_REPORT_KO.md에 있다.
