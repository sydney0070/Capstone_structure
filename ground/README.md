# 지반

## 현재 비선형 재료와 검증

2026-11-25부터 화성 기본 설정과 링 기초 해석의 네 층은 모두 `nDMaterial('DruckerPrager', ...)`를 사용한다. 현재 기준 비교 조합은 `nl_ss_block35`다. 사용 하중에서 경화층·현무암이 탄성 구간에 머물러도 재료는 탄소성이다. 모래와 암괴층에서는 소성변형을 확인했다.

| 깊이 | 재료 시편 | 물성 설정의 근거 |
|---|---|---|
| 0–0.2 m | `crust_tc` | Spohn의 경화층 강도 추정과 Brinkman의 천부 유효 강성 이전 |
| 0.2–3 m | `mgs1_ss` | 밀도 1,300 kg/m³, Dotson 직접전단 포락선에 맞춘 DP 변환 |
| 3–13 m | `blocky_phi35` | 무점착 암괴 시나리오, 다른 화성 지역·지구 유사재 마찰각 범위의 비교값 |
| 13 m 아래 | `basalt_hb_gsi55` | 지구 현무암과 Hoek–Brown 암반식을 압축 구간에서 DP로 근사 |

이 조합은 현장 시험으로 확정한 최적 재료가 아니다. 모래 강도는 시험 최소 밀도보다 조금 낮은 밀도와 높은 기초 하중 범위로 외삽한다. 암괴층 c·φ와 암반 GSI도 InSight 실측값이 아니다. 압축에 맞춘 현무암 DP는 인장강도를 과대평가하므로 인장·균열 검증에 사용하지 않는다.

추가한 소성 재하 접선 진단에서는 항복면 근처 모래 일부 점에서 국부화 가능성이 나타났다. 팽창각0°·경화0인 비연관 DP의 조건부 진단이며 지반 전체의 붕괴 판정은 아니다. 구성식·평형 검사 통과를 전체 지반의 안정성 검증으로 해석하지 않는다. [접선 계산·음향 텐서 검사](../material/nonlinear_layers_20261125/stability/stability_assessment.json).

[시편 20종·167개 재하 이력 평가](../material/nonlinear_layers_20261125/evaluation/report.html) · [3D 지반 조합과 메시·증분 비교](results/nonlinear_layers_20261125/comparison/report.html) · [분류한 논문·검증 방법 자료](../docs/references/index.html)

시편 검사는 구성식 항복강도·초기/제하 강성·질량·평형과 문헌 포락선 대조를 구분한다. 작은 구성식 오차가 독립 화성 실험 검증을 의미하지 않는다. 실제 기초 지지력, 인발, 접촉 분리와 동적 감쇠 검증은 별도다.

링 기초 기본 해석은 자중 100증분(`KrylovNewton`), 재하·가압 각32증분(`Newton`)을 사용한다. `--footing-steps`와 `--static-steps`로 수치 민감도를 비교하며, 실패 시 강도나 허용 잔차를 늘리지 않는다. 단일 알고리즘 안에서도 실패한 증분은 줄여 재시도한다. 결과는 기존 실행과 섞이지 않도록 빈 폴더에 저장한다.

```bash
python ground/ground_modeling/ground_dome_footing.py --html \
  --output-dir ground/results/my_nonlinear_run
```

`--material-profile nl_ss_block30`, `nl_ss_block40`, `nl_mgs1s_block35` 등으로 비교 조합을 명시할 수 있다. `--material-profile recommended_20261124`는 하부 두 층이 탄성인 직전 기준, `reference`는 초기 기준이다. 이전 탄성·JEZ-1 비교 결과는 아래에 이력으로 남겼다. `--sandy-regolith-density`는 과거 탄성 모래 후보를 명시하는 옵션이므로 새 비선형 조합 실행에는 사용하지 않는다.

미사용 초기 코드·지형 데이터와 중복 파일 정리 내용은 [해시를 포함한 정리 목록](../material/nonlinear_layers_20261125/cleanup/plan.json)에 있다. 활성 DEM·지진 입력과 두 저장 후보는 보존했다.

## PPT 최종 돔의 링 기초 침하

[ground_dome_footing.py](ground_modeling/ground_dome_footing.py)는 구조물팀 PPT 16–18쪽의 최종안을 사용한다. 기존 돔 좌표의 반경 20 m에 폭 14 m의 바닥판을 적용하면 기초 외경은 54 m, 원형 환대 면적은 약 1,759 m²다. 매입 깊이 3 m, 바닥판 두께 2.5 m, 레골리스 덮개 6.2 m, 돔 부재 질량 200 t를 입력 근거로 기록한다. 19쪽 응력 그림에는 바닥판 폭 12 m인 별도 해석이 표시돼 있다.

기준 해석 영역은 **90 × 90 m, 깊이 30 m**다. 링 외곽에서 측면까지 18 m, 기초 밑면에서 고정 바닥까지 27 m를 확보한다. 링을 따라 깊이 3 m까지 굴착한 메시에서 자중 평형을 구한 다음, PPT의 하향 중량 150.3 MN을 바닥판 아래에 면적으로 분배한다. 가압 조건은 순압 100.6 kPa와 돔 바닥 투영 면적 1,200 m²로 계산한 들림 120.72 MN만큼 하중을 줄인다. 이 절차는 수직 등가하중에 따른 침하를 계산한다. 허용 침하, 하부 레골리스 강도, 링–지반 접촉·인발 저항이 있어야 지지력과 들뜸을 판정할 수 있다.

`cap` 환경에서 프로젝트 루트를 기준으로 실행한다.

```bash
python ground/ground_modeling/ground_dome_footing.py --html \
  --output-dir ground/results/ground_dome_footing/my_run
```

`summary.json`에 재하 면적, 하중, 반력 평형과 감압·가압 상태별 침하를 저장한다. `footing_settlement.html`은 감압 상태의 침하를 표시한다. `--ground-width`, `--depth`, `--radial-subdivisions`로 영역과 메시를 비교할 수 있다.

과거 탄성 기준에서 `--material-profile reference --sandy-regolith-density 1300`은 0.2–3 m 모래층의 저밀도 탄성 재료를 선택한다. 3D 화면은 실제 크기에서 시작하며, 재료별 표시와 변형 확대 배율을 조절할 수 있다. [과거 탄성 1,300 kg/m³ 실행 결과](results/ground_dome_footing/ppt_final_regolith_1300/summary.json) · [과거 3D 지반 보기](results/ground_dome_footing/ppt_final_regolith_1300/footing_settlement.html)

아래는 이전 탄성 기준의 영역 비교 결과다. 새 비선형 조합의 영역·메시 수렴 결과로 재사용하지 않는다. 최대 침하는 굴착 지반의 자중 평형 이후 링 바닥 위치의 증분값이다.

| 지반 해석 영역 | 감압 최대 침하 | 가압 최대 침하 |
|---|---:|---:|
| 60 × 60 m, 깊이 12 m | 1.266 mm | 0.249 mm |
| 90 × 90 m, 깊이 30 m | 1.433 mm | 0.282 mm |
| 110 × 110 m, 깊이 30 m | 1.429 mm | 0.281 mm |
| 90 × 90 m, 깊이 42 m | 1.457 mm | 0.287 mm |

90 × 90 m, 깊이 30 m에서 감압 최대 침하는 폭을 늘리면 0.3%, 깊이를 늘리면 1.6% 변했다. 같은 영역에서 반경 방향 분할을 2에서 3으로 늘리면 1.433 mm에서 1.428 mm로 바뀌었다. [기준 실행 결과](results/ground_dome_footing/ppt_final/summary.json) · [침하 시각화](results/ground_dome_footing/ppt_final/footing_settlement.html)

## 화성 지반 단독 해석

[ground_20261107.py](ground_modeling/ground_20261107.py)의 `main()`은 데이터 불러오기 → 모델링 → 해석 → 시각화 순서로 실행한다. 절점·재료·`bbarBrick` 요소는 [mesh.py](ground_modeling/ground_hybrid/mesh.py)에서 만들고, 입력 전처리와 응답 이력·HTML 출력은 각각 `ground_hybrid/mars_inputs.py`와 `ground_hybrid/mars_results.py`에서 처리한다.

기준 모델은 InSight ELYSE 관측소를 중심으로 **60 × 60 m, 깊이 12 m의 평탄한 3차원 지반**을 만든다. DEM 표고는 표본으로 저장하고 계산 지표면은 `z=0`으로 둔다. 아래 깊이는 계산 지표면 기준이다.

| 깊이 | 등가 재료 | OpenSees 구성 법칙 |
|---|---|---|
| 0–0.2 m | 천부 경화층 · `crust_tc` | `DruckerPrager` |
| 0.2–3 m | 모래질 레골리스 · `mgs1_ss` | `DruckerPrager` |
| 3–13 m | 암괴가 섞인 조립질 레골리스 · `blocky_phi35` | `DruckerPrager` |
| 13 m 아래 | 파쇄 현무암 · `basalt_hb_gsi55` | `DruckerPrager` |

표면의 얇은 느슨한 모래층은 생략했다. 기본 깊이에서는 조립질 레골리스 중 3–12 m 구간까지 조립한다. 바닥을 고정하고 측면 롤러에서 자중 평형을 구한 뒤, 기본 동적 해석은 측면에 점성 경계를 적용한다. ELYSE 지표 가속도 3성분을 바닥 가진으로 가정한다.

## 과거 모래질 레골리스 초기 비교 이력

0.2–3 m 모래층의 비선형 비교용으로 [재료 카탈로그](ground_modeling/ground_hybrid/configs/materials.json)에 후보 세 개를 등록했다. 원지반 벌크밀도는 InSight 로봇팔 논문 §7.7의 경화층 하부 모래 범위인 1,300–1,500 kg/m³에서 1,300 kg/m³로 설정했다. Murdoch 외 표 2의 1 m 예측 파속 \(V_s=94.6\) m/s, \(V_p=157.9\) m/s를 이 밀도에서 탄성계수로 환산해 **E = 28.4 MPa, ν = 0.22, G = 11.64 MPa, K = 16.90 MPa**를 공통 적용했다. [InSight 로봇팔 연구](https://doi.org/10.1007/s11214-023-00964-0) · [Murdoch 외, 표 2 및 식 4–6](https://doi.org/10.1007/s11214-018-0553-y)

| 비교 후보 | Dotson 시험 밀도 (kg/m³) | φ (°) | c (kPa) | 재료 ID |
|---|---:|---:|---:|---|
| MGS-1 · 현무암질 모사토 | 약 1,320 | 약 35.5 | 약 0.28 | `mars_insight_sand_mgs1_dp_1300` |
| MGS-1S · 황산염 첨가 모사토 | 약 1,310 | 약 34 | 약 0.28 | `mars_insight_sand_mgs1s_dp_1300` |
| JEZ-1 · 기존 JEZ-1과 비교 | 약 1,380 | 약 24 | 약 0.23 | `mars_insight_sand_jez1_dp_1300` |

표의 강도는 Dotson 외 그림 5·6의 저응력 직접전단(수직응력 0.097–0.68 kPa) 최대 저항 포락선에서 읽은 값이다. OpenSees에서는 등가 Drucker–Prager로 변환하고 팽창각과 경화는 0으로 가정했다. 실제 시편 밀도에서 읽은 강도와 InSight의 예측 강성을 조합한 초기값이며, 링 기초 하중 범위와 현장 강도에 대한 보정은 별도다. [Dotson 외, 그림 5·6과 시험 방법](https://arxiv.org/pdf/2403.11029v1)

기존 선형탄성 기준은 `mars_insight_sandy_regolith_1300`이다. 비선형 후보를 비교할 때는 `Mars.json`의 `bottom_m: 3.0` 층에 후보 재료 ID를 지정한다.

세 후보 모두 OpenSees 등록과 구속압 5 kPa의 단일 요소 수치시험에서 초기 강성 및 해석해와의 항복강도 일치를 확인했다. 2026-11-17의 80초 시험 실행은 세 후보 모두 자중 단계 도중 시간 제한에 도달했다.

2026-11-18에 600초로 확장한 3D 계산에서 MGS-1과 MGS-1S는 자중과 두 압력 조건의 링 재하를 완료했다. 감압 최대 침하는 각각 1.4338 mm와 1.4339 mm, 가압 최대 침하는 0.2827 mm와 0.2828 mm였다. 논문 그림의 근사 강도 포락선과 대조한 단순전단 시험에서는 현재 삼축압축 대응 DP 변환의 저항이 16.6–22.0% 크게 나왔다. [문헌 대조 및 지반 실행 보고서](../material/material_release_20261101/results/insight_sand_literature_20261118/report_KO.md)

같은 지반 영역의 초기 형상과 재료 배정을 확인하는 3D 모델도 생성했다. 각 모델은 평면 90 × 90 m, 깊이 30 m, 절점 4,093개와 요소 3,504개를 사용하고 0.2–3 m 모래층에 후보를 적용한다. 아래 파일의 `summary.json`에는 `analysis_performed: false`로 조립 단계를 기록했다.

- [MGS-1 기반 지반](results/ground_dome_footing/candidate_mgs1_90x90x30/ground_geometry.html)
- [MGS-1S 기반 지반](results/ground_dome_footing/candidate_mgs1s_90x90x30/ground_geometry.html)
- [JEZ-1 기반 지반](results/ground_dome_footing/candidate_jez1_90x90x30/ground_geometry.html)

```bash
python ground/ground_modeling/ground_dome_footing.py \
  --sandy-regolith-candidate mgs1 --model-only --html \
  --output-dir ground/results/ground_dome_footing/my_mgs1_geometry
```

`--sandy-regolith-candidate` 값은 `mgs1`, `mgs1s`, `jez1` 중에서 지정한다. 하중 해석을 실행할 때는 `--model-only`를 생략한다.

## 1차원 지반 기둥

[ground_1d_column.py](ground_modeling/ground_1d_column.py)는 같은 층서·재료·깊이 방향 메시와 지진 입력을 사용한다. 수평 방향은 1 m² 단면에 요소 하나만 두고, 같은 깊이의 네 절점에서 X·Y·Z 변위를 `equalDOF`로 묶는다. 고정 바닥 위에서 깊이에 따라 달라지는 3성분 응답을 계산한다. 이 모델의 가정은 수직 전파와 수평으로 균일한 변위다. 지형, 수평 위치에 따른 응답 차이, 측면 점성 경계는 1차원 모델 범위에서 제외된다.

`cap` 환경에서 프로젝트 루트를 기준으로 실행한다.

```bash
python ground/ground_modeling/ground_20261107.py --quick --ground-history
python ground/ground_modeling/ground_1d_column.py --quick --ground-history
```

각 실행의 `ground_response.html`에서 변위를 보고 `summary.json`과 `committed_response.csv`에서 수치를 확인한다. 결과는 각각 `ground/results/ground_20261109/`, `ground/results/ground_1d_column/` 아래의 새 폴더에 저장된다. `--quick`은 선행 1초와 응답 1초를 계산한다.

PPT 링 기초 해석과 같은 깊이의 자유장 지진 응답은 1차원 기둥에 `--depth 30`을 지정한다.

## 입력과 파일 위치

- `ground_modeling/ground_hybrid/configs/`: 지반 재료·층서·메시 설정
- `ground_modeling/insight_motion.py`: 지진 입력
- `ground_modeling/insight_site_dem.py`: DEM와 위치
- `ground_modeling/material_ground_adapter.py`: 공통 재료를 지반 빌더에 연결
- [화성 지진·환경 입력](../structure/MARS_DATA/)과 `ground_modeling/`의 지형 파일

지반–기초/구조물 연결 구현은 `../SSI/`, 조립 실행은 `../integrated/`에 있다. 기존 `ground_hybrid/interface.py`, `cli.py`, `bounded_run.py`와 달 통합 파일은 호환 경로이며 실제 구현 위치를 파일 상단에 표시했다. 과거 지반 연구·관측 입력 결과는 원래 위치에 보존했다.

### JEZ-1 추가 실행 중지 및 판단 (2026-11-18)

시간 제한을 해제한 추가 JEZ-1 실행은 사용자 요청으로 약 21.3분, 자중 49.20%에서 중지했다. 기초 하중 단계에는 도달하지 않았다. 모든 관련 실행은 종료했다. 현재 모래층 기본 후보는 MGS-1, 대체 후보는 MGS-1S를 추천하며, JEZ-1 초기 조합은 약한 지반의 비교 후보로 보류한다. [수렴 지연 원인과 선정 근거](../material/material_release_20261101/results/insight_sand_literature_20261118/jez1_diagnosis_KO.md).

### 선택한 두 후보 저장 (2026-11-18)

사용자 선택에 따라 MGS-1·MGS-1S만 별도 저장본에 정리했다. [선택 재료, 시편 및 완료 지반 해석 결과](saved_candidates_20261118/README.md). 저장된 후보 해석의 0–3 m는 탄소성 DP이고, 3–30 m는 탄성이다.
