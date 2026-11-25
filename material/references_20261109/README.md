# 재료 선정·채택값 참고자료 모음

정리일: **2026-11-09 (Asia/Seoul)**

현재 지반·구조물·SSI 코드의 재료 선정과 채택 물성값을 추적할 수 있도록 원문, 실험 데이터, 기존 보정 기록, 입력값 사본을 모았다. 현재 카탈로그의 **고유 재료 ID 29종**을 포함한다. 실행 기본값으로 사용하는 재료와 이전 후보 자료는 아래 안내로 구분한다.

먼저 [현재 채택값과 근거](ADOPTED_VALUES_KO.md)를 보고, [참고자료 목록](SOURCES_KO.md)에서 해당 원문을 연다. 직접적인 실험 근거가 부족한 값은 [가정값·자료 부족 목록](GAPS_AND_ASSUMPTIONS_KO.md)에 모았다.

| 폴더 | 내용 |
|---|---|
| `00_current_values/` | 두 재료 카탈로그, 행성 설정, 채택값을 확인할 코드 사본, 29종 전체 기록 |
| `01_ground/` | InSight HP³, JEZ-1, JSC-1A LASP, 현무암, NASA 달 시뮬런트 가이드, 달 깊이별 물성 자료와 추출값 |
| `02_structure/` | S355 선택 시편 원본·변환 CSV, 출처·라이선스·보정값·검증 결과 |
| `03_candidates/` | 소결 JSC-1A·6061·AS4 참고자료와 이전 철근·6014·CFRP·micro-concrete 조사/보정 자료 |
| `04_model_docs/` | OpenSees DP·Steel02·Concrete02·PDMY02·접촉 공식 문서, DP v3.4.0 구현, 기존 변환식 검토 메모 |
| `05_project_notes/` | 내부 재료 검토·선정 보고서, 기존 참고문헌, 사용자 제공 기초 선행연구 요약 |
| `06_site_and_input/` | ELYSE StationXML, 선택된 지진 입력의 처리 메타데이터, 완료된 지반 실행의 위치·설정 맥락 |

## 값과 문헌을 연결하는 방법

- [ADOPTED_VALUES_KO.md](ADOPTED_VALUES_KO.md): 현재 화성 지반, 돔, RC 링, 달 지반 초안의 주요 값과 채택 성격.
- [SOURCES_KO.md](SOURCES_KO.md): 출처 ID → 적용 재료 → 원문 링크 → 이 폴더의 사본.
- [source_registry.json](source_registry.json): 두 원 카탈로그의 출처 메타데이터와 모은 파일 경로를 구조화한 목록.
- [all_material_records.json](00_current_values/all_material_records.json): 카탈로그별 원 기록을 모두 보존한 29종 목록. 중복 ID도 어느 카탈로그에서 왔는지 남겼다.
- [copy_manifest.json](copy_manifest.json): 기존 자료 **148개**의 원래 절대 경로, 사본 경로, 크기, SHA-256.
- [download_log.json](download_log.json): 추가 다운로드 성공 **11개**와 미확보 **1개**의 URL·응답 상태·체크섬.
- [verification.json](verification.json): 사본 체크섬, 출처 연결, 로컬 파일 존재 여부의 최종 검사 결과.

## 보관 범위

기존 코드·카탈로그·자료는 이동하거나 수정하지 않고 사본을 만들었다. 이 폴더는 날짜가 고정된 근거 자료 모음이며, 실제 실행은 기존 `ground/`, `integrated/` 진입점을 사용한다. `00_current_values/code_snapshot/`은 값 추적용 사본으로 실행용 패키지가 아니다.

현재 기본 화성 지반은 **비선형 DP 토층 2개와 탄성 현무암 기반암 1개**다. InSight 위치와 문헌을 사용하지만 모든 층이 ELYSE 현장에서 측정된 하나의 지층 자료라는 뜻은 아니다. C30·SD400·ETFE·접촉값 중에는 문헌에서 직접 추출하지 않은 가정값이 남아 있다.

관측 MiniSEED·대형 DEM·강재 전체 데이터베이스 ZIP처럼 용량이 큰 파일은 이 재료 자료 폴더에 중복 저장하지 않았다. 실제 사용하는 강재 선택 시편 원본과 지진·위치 출처 메타데이터는 포함하며, [대형 원본 위치 목록](06_site_and_input/large_originals.json)에서 나머지 원본 경로·해시를 찾을 수 있다.

원문 PDF와 시험 데이터의 저자·출처·라이선스는 기존 기록을 보존했다. 기존 내부 보고서의 과거 실행 지시나 후속 계획은 당시 기록이며, 이번 수집 작업에서 새 작업 지시로 취급하지 않았다. HTML 문서는 외부 그림·스타일을 참조할 수 있다. NIAR AS4 보고서는 원문 링크를 제공하며, PDF 파일 미확보 상태를 유지했다.
