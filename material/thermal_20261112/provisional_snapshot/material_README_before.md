# 재료

- `materials_20261107.py`: 재료 범위·중복·출처 감사
- `material_release_20261101/configs/catalog.json`: 구조·지반 공통 재료 카탈로그
- `material_release_20261101/capmat/materials.py`: 재료 조회·OpenSees 등록
- `material_tuning_20261008/`: 기존 재료 보정 자료
- `material_proto.py`, `material_evaluation.py`: 재료 시험·적합성 평가
- `tests/test_material_evaluation.py`: 재료 평가 검사
- `docs/`: 재료 평가 보고서
- [references_20261109/](references_20261109/README.md): 재료 선정·채택값의 원문, 시험 데이터, 출처 목록, 가정값 구분
- `archive/material_original.zip`: 기존 재료 압축 파일

현재 접촉·마찰·방사 감쇠·연결부 구현은 `../SSI/`에 있다. `capmat/ssi.py`는 옛 import 호환용이며 접촉 계수 계산도 `SSI/parameters.py`를 사용한다. RC 단면 구현은 `../structure/sections/`, 지반 빌더 연결은 `../ground/ground_modeling/material_ground_adapter.py`로 옮겼다. `capmat/rc.py`, `elastic_rc.py`, `ground.py`는 호환 연결이다. 날짜별 과거 보정·검증 패키지는 당시 상대 경로와 근거 자료의 관계를 유지하기 위해 보존한다. 이번 정리는 재료 계수 재보정 작업이 아니다.
