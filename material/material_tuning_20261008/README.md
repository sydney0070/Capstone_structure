# 캡스톤 재료 튜닝 기반 패키지

대상: **하이퍼볼로이드·돔 / 지구·달·화성 / 구조 재료 비교**. 조사일 2026-10-08.

2026-10-12 현재 우선 작업은 **재료 시편의 실험자료 보정**입니다. [실험 보정 실행법](material_calibration/README.md), [S355 실제 보정·독립 검증 결과](material_calibration/s355_hem320_web/RESULTS_KO.md), [전체 재료 완료 현황](material_calibration/MATERIAL_STATUS_KO.md)을 확인하세요. 기준 구조물은 이 작업에 필요하지 않습니다. 기존 전체 구조물·접촉 확장은 중단했으며 재료 전용 배포본에서 제외했습니다.

먼저 [조사 보고서](RESEARCH_REPORT_KO.md)를 읽고, [출처 목록](SOURCES.md)에서 논문과 공식 명령 문서를 확인하세요. 이 폴더는 기존 `Capstone_structure`를 수정하지 않는 독립 작업물입니다.

**23개 시편 설정**을 제공합니다. 재료·쉘·SSI 시편을 실제 OpenSees로 재하하고, JSON 변경 → 곡선 확인 → 매개변수 비교 → 실험 CSV 보정으로 이어집니다. 최적 재료를 선정한 결과나 전체 돔의 안전성 검증 결과는 아닙니다.

## 빠른 실행

이 컴퓨터에서 검증한 임시 Python은 `/private/tmp/capstone-materials-venv/bin/python`입니다. 임시 폴더는 나중에 삭제될 수 있으므로 계속 사용할 때는 아래처럼 프로젝트 가상환경을 만드세요.

```bash
cd /Users/sydney/Desktop/capstone/material_tuning_20261008
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python tune.py list
python tune.py run --case moon_regolith_dp --out my_results/moon
python tune.py run --case shell_aluminum --out my_results/aluminum
python tune.py all --out my_results/all
```

각 실행의 `history.csv`, `response.png`, `config.json`, `summary.json`을 확인합니다. 동일 출력 폴더로 재실행하면 해당 결과가 갱신되므로 비교 시 다른 폴더를 지정하세요. 결과의 재료값은 **DEMONSTRATION_NOT_CALIBRATED**입니다.

## 재료별 실행 파일

`cases/`의 파일은 해당 재료를 선택해 공통 시험기를 실행하는 진입점입니다. 실제 재료 정의는 `materials.py`, 시험 경계조건은 `experiments.py`에 한 번씩 구현해 잘못된 복사본이 생기지 않게 했습니다.

```bash
python cases/tune_moon_regolith_dp.py --out my_results/moon
python cases/tune_shell_steel.py --out my_results/steel
python cases/tune_interface_friction.py --out my_results/interface
```

| 목적 | 설정 이름 |
|---|---|
| 지구 토질·암반의 탄성 기준 | `elastic_soil` |
| 지구 모래 / 점토 | `sand_pdmy02`, `clay_pimy` |
| 달 / 화성의 건조 레골리스 후보 | `moon_regolith_dp`, `mars_regolith_dp` |
| 지구 포화토의 일정체적 시험 대안 | `sand_pm4`, `silt_pm4` |
| 콘크리트 피복 / 구속 코어 / 철근 | `concrete_cover`, `concrete_core`, `rebar_steel02` |
| 쉘 강재 / 알루미늄 / 탄소섬유복합재 | `shell_steel`, `shell_aluminum`, `shell_cfrp_elastic` |
| 쉘 콘크리트 / 소결체 / 지오폴리머 / 황콘크리트 | `shell_concrete`, `shell_sintered_regolith`, `shell_geopolymer`, `shell_sulfur_concrete` |
| 말뚝 횡 / 주면 / 선단 | `pile_py`, `pile_tz`, `pile_qz` |
| 압축 접촉 / 접촉 마찰 / 방사감쇠 | `normal_gap`, `interface_friction`, `radiation_dashpot` |

## 수동 보정과 민감도 비교

`configs/*.json`에서 `params`를 바꾸고 실행합니다. `test`는 시험 조건, `bounds`는 자동 보정할 변수와 탐색 범위입니다. 범위 역시 예시이며 문헌의 보편적 권장구간이 아닙니다.

`test.amplitude`는 일반 순환 시험의 진폭입니다. Concrete02 기본 경로는 인장 0.0001과 압축 −0.001/−0.002/−0.004/−0.008의 고정 단계이며, Viscous 기본 경로는 최대 속도 0.1 m/s입니다. 이 두 경로는 `experiments.path()`에서 변경하거나, 보정 때 실험 CSV의 `x` 이력을 사용합니다. `steps`는 단계당 출력 점수, `substeps`는 출력 사이 내부 증분 수입니다. 기본 설정 전체를 다시 생성하려면 `python tune.py init --out configs`를 실행합니다(기존 설정이 덮어써집니다).

```bash
python tune.py run --config configs/moon_regolith_dp.json --out my_results/moon_changed
python tune.py sweep --case shell_aluminum --param sig0 --values 180000 250000 300000 --out my_results/aluminum_yield
python tune.py sweep --case interface_friction --param mu --values 0.3 0.5 0.7 --out my_results/friction
```

달·화성 설정의 `planet`은 출처·조건 기록용입니다. 실제 시편 초기응력은 **`test.sigma_v`가 결정**합니다. 깊이·중력을 바꾸면 `planet_conditions.vertical_stress()`로 다시 계산해 이 값을 갱신하세요. 재료 쿠폰에는 자중을 중복 재하하지 않습니다. 전체 지반 모델에서는 각 요소·절점의 질량과 중력으로 초기응력을 만들어야 합니다.

## 실험 데이터로 자동 보정

CSV 열은 최소 `x,y`입니다. 한 행이 하나의 재하 이력 지점이며 **시간 순서를 유지**합니다. 응력–변형률을 변형률 크기순으로 정렬하면 이력거동이 사라집니다. 처음부터의 이력을 제공하고, 첫 `x`는 0이어야 합니다. 미리 재하된 시편은 그 이력도 포함해야 합니다.

| 시험 | x | y | 추가 보정 가능 열 |
|---|---|---|---|
| 3D 지반 | 공학전단변형률, fraction | 전단응력 kPa | `eps_v`(압밀 이후 체적변형률), `sigma_vertical` |
| PM4 일정체적 | 공학전단변형률 | 전단응력 kPa | `ru_equivalent`(유효 수직응력 감소율) |
| Concrete02 / Steel02 / 쉘 | 축변형률, fraction | 축응력 kPa | 쉘 `transverse_strain` |
| Py/Tz/Qz/Gap/마찰 | 변위 m | 힘 kN | 마찰 `normal_displacement` |
| Viscous | **속도 m/s** | 힘 kN | 없음 |

인장·팽창은 양수, 압축·수축은 음수입니다. Qz의 선단 압축은 음의 변위·힘입니다. `%` 변형률은 100으로 나누고, MPa는 1000을 곱해 kPa로 바꾸세요. 지반 실험의 압축 양수 체계도 변환해야 합니다. 지반·구조 시편의 재하 단계는 물리적인 시간이 아닙니다. PM4도 `SSPquad`의 정적 일정체적 시험이며 간극수압 자유도가 없습니다. `ru_equivalent`를 실제 수압으로 해석하면 안 됩니다.

실험 CSV와 같은 폴더에 아래 형식의 `targets.json`을 만듭니다.

```json
{
  "case": "moon_regolith_dp",
  "units": "kN,m,s,tonne,kPa",
  "source": "시험 보고서 번호, 시편명, 밀도, 온도, 데이터 파일 출처",
  "synthetic": false,
  "conditions_confirmed": true,
  "tests": [
    {"csv": "test_3kPa.csv", "role": "train",
     "test": {"sigma_v": 3.0}, "channels": {"y": 1.0, "eps_v": 1.0}},
    {"csv": "test_8kPa.csv", "role": "validation",
     "test": {"sigma_v": 8.0}, "channels": {"y": 1.0, "eps_v": 1.0}}
  ]
}
```

여기서 3·8 kPa는 파일 형식 설명용 가상 시험조건입니다. `conditions_confirmed`는 같은 재료 상태·시험 종류·단위로 비교함을 데이터 작성자가 확인하는 항목입니다. 코드가 재료 동일성을 자동 판정하지는 않습니다.

```bash
python tune.py fit --config configs/moon_regolith_dp.json --targets /absolute/path/targets.json --out my_results/moon_fit
```

`best_config.json`은 최적 후보, `fit_report.json`은 학습·독립 검증 오차와 최적화 종료 사유입니다. 오차는 각 채널의 max(범위, 최대 절대값)로 정규화한 RMSE이며, 채널 가중치를 적용한 뒤 시험들에 같은 가중치를 줍니다. 압력별 시험 전체를 독립 시편으로 다시 계산합니다. 누락 열·NaN·잘못된 부호·불수렴은 감추지 않습니다. 다른 시편·압력·진폭 검증 없이 작은 학습 오차만으로 채택하지 마세요.

이미 준비한 합성자료 예제로 보정 도구를 확인할 수 있습니다.

```bash
python make_demo.py
python tune.py fit --config synthetic_demo/starting_config.json --targets synthetic_demo/targets.json --maxiter 40 --out synthetic_demo/fitted
```

이는 정답을 알고 있는 **합성자료 회수 시험**이며, 실제 철근 시험 결과가 아닙니다. 예제 진값은 Fy=380,000 kPa, b=0.015입니다.

## 전체 구조물로 옮기기

- 쉘 돔·하이퍼볼로이드: `materials.shell_section()` → `ShellMITC4` 시편 경로. 전체 모델의 초기 결함·좌굴·큰 회전은 별도 검증한 기하비선형 쉘 경로가 필요합니다.
- 격자 돔·하이퍼볼로이드: `Steel02`/`Concrete02` → fiber section → `forceBeamColumn` 또는 `dispBeamColumn`; 축력만 받는 부재라면 `corotTruss` 후보. 재료 비선형성과 기하비선형성은 별개입니다.
- `shell_concrete` 등 취성 쉘 4종은 **탄성 강성 보정**만 합니다. 철근·균열·파괴를 포함한 RC 쉘 모델이 아닙니다. `Concrete02`를 `ShellMITC4`나 `SSPbrick`에 직접 붙일 수 없습니다.
- `shell_cfrp_elastic`은 국부 축 방향 직교이방성만 나타냅니다. 적층각 변환, 파손·박리·수지의 온도 의존성을 구현한 복합재 설계 모델이 아닙니다.
- `shell_steel`/`shell_aluminum`의 J2 모델은 등방경화입니다. 반복 하중의 바우싱거 효과가 필요하면 다른 다축 구성모델과 실험 검증을 추가해야 합니다.
- J2에는 질량밀도 입력 인자가 없으므로 **전체 쉘 모델의 절점 질량을 별도로 적분**해야 합니다. 예제 JSON의 `rho`는 비교·질량 조립용이며 J2가 자동 사용하지 않습니다.
- `rc_section()`/`rc_beam()`은 2D RC 보·기둥 조립 예시입니다. 쉘 보강재나 실제 돔을 대신하지 않습니다.
- `pile_spring()`은 동일 위치의 지반·구조 절점을 연결합니다. 실제 깊이별 용량은 `pile_capacities()`로 힘 단위로 환산하고, 얕은 링기초에 말뚝 스프링 상수를 그대로 사용하지 않습니다.
- `simulate()`는 `ops.wipe()`를 호출합니다. 기존 전체 모델을 만든 뒤 호출하지 마세요. `materials.py`의 생성 함수들은 wipe하지 않습니다.

## 확인한 범위

OpenSees 3.8.0 / Python 3.12 환경에서 23개 기본 실행 완료. `verify.py`는 탄성 해석해, 콘크리트 피크, 반력 평형, 마찰 상한, 압축 접촉, 단위 환산 등을 확인합니다. 상세 기록은 [검증 문서](VALIDATION.md)와 `verification.log`에 있습니다. 증분·메시·조건 변화 및 실험 타당성은 최종 사용 시 추가 확인 대상입니다.

현재 원시 실험자료가 없어 문헌 수치 전체를 보정값으로 확정하지 않았습니다. 특히 소결체·황콘크리트·복합재·행성 레골리스 설정은 후보 비교를 위한 초기값입니다.
