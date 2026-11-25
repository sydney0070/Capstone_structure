# 현재 1107 계열 모델 안내: 재료·위치·코드 구조

확인일: 2026-11-09  
프로젝트: `/Users/sydney/Desktop/capstone/Capstone_structure`

지반 단독 실행은 **ground_20261107.py**, 통합 실행은 **mars_ssi_integrated_20261109.py**로 분리했다. 첨부 화면 형식의 재료별 3D HTML을 두 경로에서 생성한다. 이 문서의 돔·링 설명은 통합본에만 해당한다.

이 문서는 **1109 수정사항이 반영된 통합 실행의 기본 설정**을 설명한다. 재료 카탈로그에 존재하는 모든 후보를 현재 채택한 재료로 나열하지 않는다. 실행 옵션을 바꾼 경우 실제 적용값은 해당 실행의 `summary.json`과 소스 해시로 확인해야 한다.

## 1. 현재 실행하는 모델

화성의 3차원 지반 위에 직경 40 m 강재 돔과 RC 링을 결합한 연구 모델이다. 기본 지반 영역은 60×60 m, 깊이는 12 m이다. 중력은 3.711 m/s²이며 단위는 m·N·s·kg·Pa이다. 돔 강재 목표 질량은 200 t이고, ETFE와 RC 링 질량은 별도로 더해진다. 기본 링은 반지름 20 m, 폭 1 m, 높이 0.6 m, 12개 보 요소다.

실행 진입점은 [mars_ssi_integrated_20261109.py](/Users/sydney/Desktop/capstone/Capstone_structure/structure/mars_ssi_integrated_20261109.py)이다. [ground_20261107.py](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_20261107.py)는 지반 구성·입력·평형 검사를 연결하며, 직접 실행하면 구조물·기초·연결부 없는 지반 단독 해석을 수행한다.

기본 통합 실행은 압력 0, 레골리스 덮개 두께 0이다. ETFE는 질량만 반영하며 막의 인장 강성·기밀·압력 폐합 구조를 해석하지 않는다.

## 2. 실제 채택한 지반 재료

깊이는 평탄한 계산 지표면에서 아래 방향으로 잰다. 아래 값은 기본 강성 배율 1.0 기준이다.

| 깊이 | 재료 ID / 의미 | 구성 법칙 | E | ν | 밀도 [kg/m³] | φ | c |
|---|---|---|---:|---:|---:|---:|---:|
| 0–0.2 m | `mars_insight_duricrust_proxy` / InSight 천부 경화층 대용 | DruckerPrager | 11.48 MPa | 0.28 | 1,200 | 30° | 5.8 kPa |
| 0.2–8 m | `mars_jez1_peak` / JEZ-1 화성 모사토의 최대강도 기준 | DruckerPrager | 30 MPa | 0.30 | 1,840 | 40.6° | 3.1 kPa |
| 8–12 m | `basalt_elastic_analogue` / 현무암 기반암 대용 | ElasticIsotropic | 40 GPa | 0.30 | 2,900 | 해당 없음 | 해당 없음 |

두 토층의 팽창각 ψ는 0이며, 현재 등록 경로의 경화 계수는 0이다. 토층은 비선형 소성을 유지하되 탄성계수는 구속압에 따라 자동 변화하지 않는 설정이다. 하부 기반암은 탄성이며 절리·파쇄·인장 파괴를 포함하지 않는다.

재료값은 [지반 재료 카탈로그](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_hybrid/configs/materials.json), 층 구분의 기본값은 [mars.json](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_hybrid/configs/mars.json)에 있다. 통합 실행은 `mars.json`의 소형 영역을 60×60×12 m로 바꾸고, 예시 렌즈 영역은 `lenses=[]`로 비활성화한다. 설정 파일만 보고 렌즈가 현재 통합 모델에 들어간 것으로 판단하면 안 된다.

**이 세 층은 ELYSE에서 측정한 하나의 지층 기둥이 아니다.** InSight 천부 자료·JEZ-1 모사토·현무암 대용값을 결합한 가정 모델이다. 특히 JEZ-1의 E=30 MPa와 ν=0.30은 가정값이다. 경화층 E도 정적 경화층 시험으로 직접 보정한 값이 아니다. 지층 경계 깊이 역시 현장 시추로 확인되지 않았다.

## 3. 실제 채택한 구조 재료와 접촉

| 부위 | 재료 ID / 재료 | 현재 반영 방식 | 주요 값 |
|---|---|---|---|
| 돔 부재 | `steel_s355_coupon` / S355 계열 쿠폰 보정 강재 | Steel02 섬유 단면, SHS 중공 정사각형 부재, forceBeamColumn | E≈204.17 GPa, Fy≈384.17 MPa, 밀도 7,850 kg/m³ |
| 링 콘크리트 | `concrete_c30_prototype` / C30 예제 | Concrete02 섬유, 양단 HingeRadau 영역과 탄성 내부 | 압축강도 30 MPa, E=30 GPa, 밀도 2,400 kg/m³ |
| 링 철근 | `rebar_sd400_prototype` / SD400 예제 | Steel02 철근 섬유 | E=200 GPa, Fy=400 MPa, 밀도 7,850 kg/m³ |
| 외피 | ETFE 3겹 | 등가 절점 질량만 반영 | 1겹 0.2 mm, 밀도 1,750 kg/m³, 면밀도 1.05 kg/m² |
| 지반–링 접촉 | 비접착 Coulomb 마찰·들림 | flatSliderBearing, 접촉 강성으로 하중 전달 | 기본 마찰계수 μ=0.35 |

구조 재료의 실제 수치는 [공통 카탈로그](/Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101/configs/catalog.json)에서 읽는다. S355라는 명칭과 보정된 Fy≈384.17 MPa를 구분해야 한다. 쿠폰 보정을 실제 SHS 부재의 국부 좌굴·용접·반복 파괴 검증으로 확대할 수 없다.

현재 통합 링은 **C30/SD400 예제 RC**를 사용한다. 재사용하는 큰 돔 파일의 다른 함수·출력에 등장하는 황 콘크리트 등의 재료를 이 실행이 채택한 것으로 읽으면 안 된다. C30·SD400 역시 화성 시공 재료로 선정·검증된 값은 아니다. 링의 힌지 길이와 철근 부착·전단·구속 효과도 충분히 보정되지 않았다.

링 접촉은 기본 12개 위치의 점 하중 전달이다. 명목 접촉 면적을 사용하더라도 유한 폭의 면 접촉을 직접 해상한 모델은 아니다. 별도의 전면·환형 판기초 초안과 구분한다.

## 4. 선정한 위치와 지진자료

| 항목 | 현재 기본 설정 |
|---|---|
| 행성 | 화성 |
| 위치 기준 | InSight 착륙지 주변 ELYSE 관측소 StationXML |
| 위도 / 경도 | 북위 4.502384° / 동경 135.623447° |
| 좌표축 | x: 동쪽, y: 북쪽, z: 위쪽 |
| 계산 지표면 | 해당 위치 DEM 표고를 기준으로 하는 평탄한 국부 z=0 면 |
| 지형자료 | `DTEEC_037262_1845_036761_1845_U01.IMG` |
| 위치 메타데이터 | `stationxml.xb.elyse.2019.151.xml` |
| 지진자료 | `MARS_DATA/S1222a_ACC.mseed`의 ELYSE 3성분 |

위경도는 StationXML에서 직접 확인했다. 착륙체와 실제 설치된 SEIS 사이의 측량 좌표를 새로 확정한 것은 아니다. DEM의 실제 요철은 표본으로 저장하지만 현재 통합 지반 요소 표면에는 적용하지 않는다. DEM 중심 표고와 StationXML 표고의 차이도 결과에 보존한다.

지진자료는 관측 파형을 사용하는 연구 입력이다. 지표 관측 가속도를 기저 운동으로 사용하는 가정이 남아 있으며, 역산된 기저 운동이나 설계 지진을 확보한 것으로 표현하지 않는다. 기본 실행은 파형 배율 1이고, 별도 PGA 증폭은 가상 시나리오로 기록한다.

## 5. 주요 코드 트리

아래는 주요 실행 의존성과 별도 연구 경로만 추린 트리다. 모든 파일을 열거한 것은 아니다.

```text
/Users/sydney/Desktop/capstone/Capstone_structure/
├── MODEL_OVERVIEW_20261109_KO.md       ← 이 안내서
├── structure/
│   ├── mars_ssi_integrated_20261109.py     ← 새 화성 RC 링–돔 통합 진입점
│   ├── mars_dome_ring_ground_20261107.py  ← 이전 이름의 호환 진입점
│   ├── ring_foundation_rc.py             ← C30/SD400 링·섬유 단면 구성
│   ├── dome_in_mars_tuned.py             ← 돔 형상·SHS·강재 단면·질량 함수 재사용
│   ├── verify_revision_20261107.py       ← 수정사항 회귀 검사
│   ├── planet_foundation_demo_20261109.py ← 별도 화성/달 판기초 연구 진입점
│   ├── planet_foundations_20261109.py     ← 별도 전면·환형 판기초와 연결부
│   └── results_mars_ssi_integrated_20261109/ ← 새 통합 해석 결과
├── ground/ground_modeling/
│   ├── ground_20261107.py                ← 지반 단독 실행·공용 지반 함수
│   ├── ground_only_20261109.py            ← 구조물 없는 해석·실행 관리
│   ├── response_viewer_20261109.py        ← 현재 결과를 재료별 HTML에 연결
│   ├── ground_history_20261109.py        ← 적분점·소성·일 이력
│   ├── insight_site_dem.py              ← DEM 읽기
│   ├── insight_motion.py                ← 지진자료 전처리
│   ├── ground_modeling.py               ← 재사용하는 요소 형상 검사 등
│   ├── ground_hybrid/
│   │   ├── config.py                    ← 설정 읽기·검증
│   │   ├── materials.py                 ← 카탈로그 → OpenSees 재료 등록
│   │   ├── mesh.py                      ← 기본 3차원 지반 요소·질량·자중
│   │   ├── conforming_mesh.py           ← 별도 고정 다각형 접촉 경계 메시
│   │   ├── solver.py                    ← 단계 재하·재시도·평형 검사
│   │   ├── interface.py                 ← 일반 기초 연결·응답 보조
│   │   ├── process_control.py           ← 시간 제한·프로세스 정리
│   │   ├── run_journal.py               ← 해석 시도 로그·체크포인트
│   │   └── configs/
│   │       ├── mars.json                ← 화성 기본 층서·영역 설정
│   │       ├── moon.json                ← 별도 달 기본 설정
│   │       └── materials.json           ← 지반 재료값·출처·가정
│   ├── stationxml.xb.elyse.2019.151.xml  ← 관측소 위치
│   ├── DTEEC_037262_1845_036761_1845_U01.IMG ← DEM
│   └── MARS_DATA/S1222a_ACC.mseed        ← 지진 입력
├── material/
│   ├── materials_20261107.py            ← 재료 범위·중복·출처 감사
│   └── material_release_20261101/
│       ├── configs/catalog.json         ← 구조·지반 공통 재료 카탈로그
│       └── capmat/
│           ├── materials.py             ← 공통 재료 조회·등록
│           ├── rc.py                    ← RC 섬유 단면
│           ├── ssi.py                   ← 마찰·들림 접촉 요소
│           └── mechanics.py             ← 접촉 강성 등 계산
└── research_20261109/
    ├── FOLLOWUP_REPORT_KO.md             ← 보완 작업·검증 결과
    ├── FAILURE_CAUSE_ANALYSIS_KO.md      ← 반복 실패 원인·버전 판단
    ├── followup_20261109/                ← 후속 검증·메시 비교
    └── failure_analysis_20261109/        ← 반복 잔차·재료·복원 진단
```

통합본 실행 흐름은 `통합 진입점 → 지반/링/돔 조립 → 전체 자중 재하 → 지진 입력 → 결과 저장`이다. 현재 RC 링 통합 경로가 `conforming_mesh.py`나 별도 판기초 빌더를 자동으로 사용하는 것은 아니다. 지반 자중 평형 후 기초를 설치하는 별도 재하 이력은 아직 적용 전이다.

## 6. 실행과 확인할 결과

**지반만 실행** — 구조물·기초·연결부 없음:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_20261107.py \
  --quick --ground-history
```

지반 단독 결과는 `ground/ground_modeling/results_ground_20261109/` 아래 새 실행 폴더에 저장된다. `ground_response.html`을 연다. 실제 총변위·시간·재료 정보를 `ground_response.json`에도 보존한다. 고정 바닥·점성 측면의 지반 경계 요소는 포함하지만 기초 접촉·돔·RC 요소는 만들지 않는다.

**통합본 실행** — 지반·RC 링·돔·접촉을 함께 해석:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/structure/mars_ssi_integrated_20261109.py \
  --site insight --quick --ground-history
```

`--quick`는 선행 1초＋응답 1초다. 생략하면 기본 선행 5초＋응답 2초다. `--static-only`는 지진 해석을 생략하고 전체 자중 평형만 수행한다. `--elastic-ground`는 명시적인 탄성 대조 옵션이며 기본값이 아니다. 기본 시간 제한은 180초이고 제한 내 완료되지 않으면 미완료 상태를 남긴다.

위치 점검만 실행:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_20261107.py \
  --inspect-site
```

통합본 결과는 `structure/results_mars_ssi_integrated_20261109/` 아래의 새 실행 폴더에 저장된다. 아래 표는 통합본 출력이다. 과거 결과는 원래 폴더에 보존했다.

| 결과 | 용도 |
|---|---|
| `summary.json` | 위치·지반 재료·단면·실행 옵션·소스 해시·응답 |
| `run_status.json` | 완료·실패·시간 초과 여부; 결과 판단 전에 확인 |
| `site_dem_samples.csv` | 실제 DEM 표본과 계산 평면 비교 |
| `static_ground_history.csv` | `--ground-history` 사용 시 자중 재하 적분점 이력 |
| `dynamic_ground_history.csv` | 같은 옵션의 동적 적분점 이력 |
| `coupled_response.csv` | 외부 출력 시간의 응답 |
| `committed_response.csv` | 성공한 세부 시간 단계의 응답 |
| `mars_integrated_response.html` | 지반·링·돔 응답 시각화 |

동적 파일과 HTML 등은 해당 단계가 실행·완료됐을 때 생성된다. 중간 파일이 있다는 사실만으로 전체 실행 완료를 판단하지 않는다.

## 7. 현재 판단

현재 수정된 1107 계열은 연구 기준으로 유지할 만하다. 40 m 달 기초 결합의 수렴, 비선형 메시 독립성, 현장 지층·낮은 구속압 물성 보정, 실제 부재·접촉 검증은 미완료다. 실행 경로와 HTML을 분리·수정했으며 재료 계수는 유지했다. [시각화 및 실행 분리 보고서](/Users/sydney/Desktop/capstone/Capstone_structure/research_20261109/VISUALIZATION_SPLIT_20261109_KO.md)에 변경·검증 범위를 기록했다.
