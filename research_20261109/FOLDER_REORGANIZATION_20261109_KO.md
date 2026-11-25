# 폴더 정리 보고서

작성일: 2026-11-09  
대상: `/Users/sydney/Desktop/capstone/Capstone_structure`

## 1. 현재 기준 구조

```text
Capstone_structure/
├── README.md                 ← 전체 안내·현재 실행 명령
├── ground/                   ← 지반·지형·지진 입력
│   ├── ground_modeling/      ← 지반 단독 실행·구성 코드·원본 데이터
│   └── results/              ← 현재 지반 단독 결과
├── material/                 ← 재료 카탈로그·보정·재료 시험
│   ├── material_release_20261101/
│   ├── material_tuning_20261008/
│   ├── tests/                ← 재료 평가 검사
│   ├── docs/                 ← 재료 평가 보고서
│   └── archive/              ← 원본 재료 ZIP
├── structure/                ← 돔·링·판기초·구조 단면
│   └── sections/             ← RC 섬유 단면·탄성 등가 RC
├── SSI/                      ← 지반–기초·구조물 연결 구현
│   ├── contact.py
│   ├── parameters.py
│   ├── interface.py
│   └── connectors_20261109.py
├── integrated/               ← 지반·기초·구조물 조립 해석
│   ├── mars_ssi_integrated_20261109.py
│   ├── planet_foundation_demo_20261109.py
│   ├── lunar_main.py / lunar_main_dem.py
│   ├── ground_reference_demo.py
│   ├── tests/
│   └── results/              ← 현재·과거 통합 결과
├── docs/                     ← 모델 안내·검토·참고 보고서
├── archive/                  ← 기타 기존 보관 파일
└── research_20261109/         ← 연구 기록·비교 결과·소스 스냅샷
```

각 주요 폴더에 README를 추가했다. 현재 실행 파일·재료·위치·전체 트리는 [모델 안내](/Users/sydney/Desktop/capstone/Capstone_structure/docs/MODEL_OVERVIEW_20261109_KO.md)에 반영했다.

## 2. 실제 구현의 이동

| 이전 구현 위치 | 현재 실제 구현 |
|---|---|
| `structure/mars_ssi_integrated_20261109.py` | `integrated/mars_ssi_integrated_20261109.py` |
| `structure/mars_dome_ring_ground.py` | `integrated/mars_dome_ring_ground.py` |
| `structure/planet_foundation_demo_20261109.py` | `integrated/planet_foundation_demo_20261109.py` |
| `ground/ground_modeling/lunar_main*.py` | `integrated/lunar_main*.py` |
| `ground/ground_modeling/lunar_config.py` | `integrated/lunar_config.py` |
| `ground_hybrid/cli.py`, `bounded_run.py` | `integrated/ground_reference_demo.py`, `reference_bounded.py` |
| `ground/ground_modeling/lunar_structure.py` | `structure/lunar_structure.py` |
| `capmat/rc.py`, `elastic_rc.py` | `structure/sections/rc.py`, `elastic_rc.py` |
| `capmat/ssi.py` | `SSI/contact.py` |
| `ground_hybrid/interface.py` | `SSI/interface.py` |
| `capmat/mechanics.py`의 접촉·방사 계수 함수 | `SSI/parameters.py` |
| 판기초 파일의 돔 연결·인장 타이 생성 | `SSI/connectors_20261109.py` |
| `capmat/ground.py` | `ground/ground_modeling/material_ground_adapter.py` |

파일 단위 구현 이동은 14개이며, 접촉 계수·돔 연결부 함수도 별도 추출했다. 구조물 형상과 판기초 구조 요소는 `structure`에 두고 접촉·힌지·스프링·타이는 `SSI`에서 등록한다. 지반 공용 패키지의 SSI 편의 이름은 지연 조회하도록 해 지반 기본 import가 구조물 인터페이스를 즉시 불러오지 않게 했다.

기존 연구·검증 스크립트가 옛 경로를 많이 사용하므로 옛 위치에는 짧은 **호환 모듈**만 남겼다. 실제 구현을 두 군데 복사해 관리하지 않는다. 호환 파일 상단에는 현재 구현 위치를 표시했다. 새 실행·수정은 위 표의 실제 구현을 기준으로 한다.

`ground_modeling.py`의 기본 실행도 지반 단독 경로로 연결했다. 구형 DEM 기능은 명시적인 기존 선택 경로로 보존했다. 달 공통 설정은 지반·구조물·하중을 함께 포함하므로 통합 폴더에 두었다.

## 3. 결과·문서·보관 파일

- 현재 지반 결과: `ground/results/ground_20261109/`
- 현재 화성 통합 결과: `integrated/results/mars_ssi_20261109/`
- 기존 통합 결과: `integrated/results/mars_integrated_20261107/`, `mars_integrated_legacy/`, `reference_legacy/`
- 루트의 코드 검토·개정 보고서: `docs/reviews/`
- 사용자 제공 지반·재료 검토 문서: `docs/reference/`
- 현재 모델 안내: `docs/MODEL_OVERVIEW_20261109_KO.md`
- 재료 평가 문서·검사: `material/docs/`, `material/tests/`
- 루트의 `material.zip`: `material/archive/material_original.zip`
- 기타 기존 INI 보관 파일: `archive/misc/`

문서·자료·결과 이동은 13개 파일/폴더다. 원본 지진·DEM 데이터와 재료 보정 자료는 기존 도메인 폴더 안에 유지했다. 수 GB 규모의 원본 자료를 복제하지 않았다. 기존 구조물 단독 결과와 날짜별 연구·보정 자료도 각 도메인 안에서 보존했다.

과거 보고서·실행 결과에 기록된 옛 절대 경로와 소스 해시는 당시 기록이므로 일괄 고쳐 쓰지 않았다. 현재 안내서와 새 README는 새 위치를 사용한다. 과거 경로를 찾을 때는 아래 이동 대응표를 참고한다.

- [소스 이동 대응표](/Users/sydney/Desktop/capstone/Capstone_structure/research_20261109/folder_reorganization_20261109/source_moves.json)
- [자료·결과 이동 대응표](/Users/sydney/Desktop/capstone/Capstone_structure/research_20261109/folder_reorganization_20261109/asset_moves.json)
- 정리 전 소스: `folder_reorganization_20261109/before/`

## 4. 경로 수정과 검증

모듈 검색·직접 실행·원본 데이터·소스 해시·기본 결과 폴더 경로를 함께 수정했다. 프로젝트 외부 작업 폴더에서도 절대 경로 명령이 동작하도록 했다. 초기 회귀 검사에서 호환 파일의 프로젝트 루트 계산이 한 단계 짧은 오류를 발견해 수정했고, 외부 작업 폴더 실행 검사를 추가했다. 최초 실패 로그는 그대로 보존했다.

최종 회귀 결과는 기존 38개, 시각화·분리 4개, 새 폴더 경로 검사 2개, 재료 평가 13개로 **총 57개 통과**다. 새 검사는 옛 import가 실제 새 구현으로 연결되는지와 프로젝트 밖에서 실행 경로를 찾는지를 확인한다. 구형 달·일반 예제 모듈도 import를 확인했다. 이들 구형 모델 전체의 물리 검증을 새로 완료했다는 뜻은 아니다.

실제 지반·화성 통합 실행은 `--quick --ground-history`로 프로젝트 밖에서 수행했다. 지반 단독 결과의 구조물·기초·접촉 개수는 모두 0이었다. 정리 전 완료 결과와 비교한 지반 응답 CSV와 통합 응답 CSV의 최대 수치 차이는 각각 **0**이었다. 작은 달 판기초·힌지 사례도 별도로 실행해 SSI 이동 경로를 점검했다. 개별 감독 상태와 최종 증거는 아래 검증 파일을 기준으로 한다.

검증 기록은 `folder_reorganization_20261109/verification.json`, `response_comparison.json`과 개별 로그에 남긴다. 재료 JSON과 지반 층서·중력 설정은 이번 정리로 변경하지 않았다. `.git`과 기존 사용자 파일 삭제 상태는 건드리지 않았다.

## 5. 현재 실행 명령

지반만:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_20261107.py \
  --quick --ground-history
```

통합본:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/integrated/mars_ssi_integrated_20261109.py \
  --quick --ground-history
```

이번 작업은 코드 소유 위치·실행 경로·자료 위치를 정리한 것이다. 비선형 메시 독립성, 대형 달 모델 수렴, 실측 물성 보정 등 기존 미완료 연구 문제를 해결한 것으로 표시하지 않는다.
