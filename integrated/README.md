# 통합 해석

| 파일 | 용도 |
|---|---|
| `mars_ssi_integrated_20261109.py` | 현재 화성 지반–RC 링–돔 통합본 |
| `planet_foundation_demo_20261109.py` | 화성 전면 / 달 환형 판기초–돔 연구 초안 |
| `mars_dome_ring_ground.py` | 이전 합성 지반 통합 모델 |
| `lunar_main.py`, `lunar_main_dem.py` | 기존 달 지반–구조물 통합 실행 |
| `lunar_config.py` | 기존 달 지반·구조물·하중의 공통 설정 |
| `ground_reference_demo.py`, `reference_bounded.py` | 일반 지반–판기초·예제 기둥 연결 시험 |

현재 화성 통합 결과는 `results/mars_ssi_20261109/`에 저장한다. 이전 통합 결과는 `results/mars_integrated_20261107/`, `results/mars_integrated_legacy/`, `results/reference_legacy/`에 보존했다. 지반만 보려면 `../ground/ground_modeling/ground_20261107.py`를 사용한다.

모든 통합본을 검증된 현재 해법으로 취급하지 않는다. 달 구형 통합 모델과 일반 예제는 보존한 연구 경로이며, 현재 기본은 날짜가 붙은 화성 통합본과 별도 판기초 초안이다.
