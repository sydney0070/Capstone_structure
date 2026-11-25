# Capstone 구조·지반 해석

## 2026-11-25 재료·지반·문헌 작업

현재 화성 링 기초 지반은 네 층 모두 탄소성 DruckerPrager를 사용한다. 연구용 기준은 경화층, MGS-1 단순전단 보정 모래, 암괴층 φ35°, 현무암 GSI55 압축 근사다. 구성식·평형 검사와 현장 타당성·소성 국부화·메시 검증은 구분한다.

- [비선형 수치시편 20세트 평가 및 검증 기준](material/nonlinear_layers_20261125/evaluation/report.html)
- [90 × 90 × 30 m 지반 조합과 수치 민감도 비교](ground/results/nonlinear_layers_20261125/comparison/report.html)
- [논문·보고서·공식 자료 분류 목록](docs/references/index.html)
- [GitHub 업로드 범위와 별도 보관한 대용량 파일](docs/GITHUB_UPLOAD_20261125.json)

프로젝트 루트의 OpenSeesPy 환경에서 `python ground/ground_modeling/ground_dome_footing.py --html --output-dir ground/results/my_run`으로 실행한다. 논문 HTML과 보고서는 내려받은 뒤 브라우저에서 열 수 있다. 팀 공동 화성 DEM·관측·지진 입력 자료와 큰 백업 ZIP은 로컬에 보존하며 업로드 목록에 제외 범위를 남겼다.

현재 코드는 **지반 / 재료 / 구조물 / SSI 연결부 / 통합 해석**으로 구분한다. 최신 실행은 아래 두 명령을 사용한다.

| 폴더 | 역할 | 주요 파일 |
|---|---|---|
| `ground/` | 지반 단독, 층서·지형·지진 입력·지반 해석 | `ground_modeling/ground_20261107.py` |
| `material/` | 재료 카탈로그·보정·재료 검증 | `materials_20261107.py`, `material_release_20261101/configs/catalog.json` |
| `structure/` | 돔·RC 링·판기초·부재/단면 생성 | `dome_in_mars_tuned.py`, `ring_foundation_rc.py`, `planet_foundations_20261109.py` |
| `SSI/` | 접촉·미끄럼·들림·기초–지반 및 돔 연결 | `contact.py`, `parameters.py`, `interface.py`, `connectors_20261109.py` |
| `integrated/` | 지반·기초·구조물을 조립해 해석 | `mars_ssi_integrated_20261109.py`, `planet_foundation_demo_20261109.py` |
| `docs/` | 현재 모델 안내·기존 검토·참고 보고서 | `MODEL_OVERVIEW_20261109_KO.md` |
| `research_20261109/` | 연구 기록·원본 스냅샷·비교 결과 | `folder_reorganization_20261109/` 등 |
| `archive/` | 분류가 애매한 기존 보관 파일 | `misc/` |

## 지반만 실행

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_20261107.py \
  --quick --ground-history
```

구조물·기초·접촉을 만들지 않는다. 결과는 `ground/results/ground_20261109/`의 새 실행 폴더에 저장한다. `ground_response.html`을 연다.

층서와 지진 입력을 공유하는 [1차원 지반 기둥 비교 모델](ground/README.md#1차원-지반-기둥)도 실행할 수 있다.
구조물팀 PPT 최종 돔의 링 바닥 하중과 침하는 [링 기초 지반 해석](ground/README.md#ppt-최종-돔의-링-기초-침하)에서 확인한다.

## 통합본 실행

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/integrated/mars_ssi_integrated_20261109.py \
  --quick --ground-history
```

결과는 `integrated/results/mars_ssi_20261109/`의 새 실행 폴더에 저장한다. `mars_integrated_response.html`을 연다. 별도 화성/달 판기초 초안은 `integrated/planet_foundation_demo_20261109.py`다.

`--quick`는 선행 1초＋응답 1초의 연결 점검이다. 현재 기본 화성 층서는 모두 탄소성 DP로 지정돼 있다. 지반 단독 12 m 모델은 암괴층 안에서 끝나며, 링 기초 모델은 30 m까지 해석한다.

## 자료 찾기

- [현재 재료·위치·코드 안내](docs/MODEL_OVERVIEW_20261109_KO.md)
- [재료 선정·채택값 참고자료 모음](material/references_20261109/README.md)
- [폴더별 설명: ground](ground/README.md), [material](material/README.md), [structure](structure/README.md), [SSI](SSI/README.md), [integrated](integrated/README.md)
- [폴더 정리 보고서](research_20261109/FOLDER_REORGANIZATION_20261109_KO.md)

이전 실행·검증 스크립트를 위해 일부 옛 경로에는 짧은 호환 모듈을 남겼다. 실제 구현은 위 폴더가 기준이다. 날짜가 붙은 과거 재료 보정·검증 패키지와 연구 소스 스냅샷은 재현성을 위해 내부 구조를 보존한다. 과거 결과의 원래 절대 경로·소스 해시는 재작성하지 않으며, 이동 대응표는 정리 보고서에서 찾을 수 있다.
