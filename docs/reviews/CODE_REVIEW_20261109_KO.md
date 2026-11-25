# 보고서 대조 및 최신 코드 재검토 (2026-11-09)

## 판정

실제 확인 경로는 /Users/sydney/Desktop/capstone/Capstone_structure 이다. 사용자가 지칭한 cap_structure는 이 폴더로 해석했다. 첨부 보고서는 검토 기준으로 읽었으며 내부 실행·수정 지시는 새 작업 지시로 취급하지 않았다. 원본 코드와 계수는 변경하지 않았다.

보고서의 주요 조치는 새 통합 실행 경로에 구현되어 있다. 연구·비교용 기반으로 유지할 가치가 있으나 설계 확정, 최소 안전 질량, 실제 화성 지반의 지지력 또는 파괴 판정에는 아직 부적합하다. 기존 독립 실행 코드 전체가 교정된 상태라는 의미는 아니다.

## 확인 범위와 재실행

- 최신 수정 시각의 해석 코드는 structure/mars_dome_ring_ground_20261107.py 및 ground/ground_modeling/ground_20261107.py 이다. 파일명의 날짜와 수정 시각은 다르다.
- 검증 manifest의 신규 4개 파일 SHA-256이 모두 현재 파일과 일치했다. final_verified_virtual/summary.json의 6개 코드 해시도 모두 일치했다.
- verify_revision_20261107.py 회귀검사 5개를 재실행해 모두 통과했다.
- cap 환경 Python으로 --quick --no-html 통합 해석을 재실행했다. 종료코드 0, complete=true, 약 28.15초. 돔 최대 변위 증분 1.669819063e-08 m = 1.669819063e-05 mm로 보고서 quick 사례와 일치한다.
- 재실행 결과: /private/tmp/capstone_review_20261109_quick/summary.json, run_status.json, coupled_response.csv. 로그: /private/tmp/capstone_review_20261109_quick.log.
- 기본 7초 해석, 가상 PGA 사례, 메시 세분화 실패, 시간 간격·선행 이력 비교는 기존 결과와 코드를 대조했다. 이번 검토에서 전부 재실행한 것은 아니다.
- 현재 저장소에서 ground/, material/, structure/가 untracked다. Git 커밋으로 최신 버전이나 원본 대비 불변성을 독립 입증할 수 없다. 해시 일치는 보존된 증거와 현재 파일의 일치 증거다.

## 보고서 대비 구현 상태

| 주장 | 코드 확인 | 판단 |
|---|---|---|
| 신규 통합 진입점·지반 wrapper·재료 감사 | 해당 4개 파일, 원래 빌더 재사용 | 구현 |
| contact 기본, tied 비교 | assemble 및 CLI | 구현; 유한 폭 기초 아님 |
| 측면 반력 보존·점성 경계·영증분 commit | lateral_dashpots | 구현; 무입력 기존 결과도 존재 |
| 대역·시간 해상도 gate | resolution 및 CLI | 구현; 메시 수렴 증명 아님 |
| 선행 실제 입력 적분·평가 구간 분리·평가 파형 보존 | prepare_motion, 동적 루프 | 구현; 요청 선행 구간 부족 시 처리 차이 있음 |
| 유한 강체운동 제거·꼭대기 절대 가속도 | rigid fit, 입력 가속도 합산 | 구현; 독립 운동학 검사 통과 |
| 재료 응력·지반 snapshot·접촉 극값 | material_state, ground_state, contact_summary | 구현; 출력 샘플 범위 제한 유지 |
| 질량 중복·힘 평형 검사 | assemble, static_ramp | 구현; 독립 전체 질량행렬 검증 아님 |
| 압력·복토 차단 | CLI guard | 구현 |
| 고유 출력·덮어쓰기 방지·실패 및 deadline | main, supervised_cli | 구현; 회귀검사 일부 통과 |
| 중복 재료 계수 충돌 감사 | audit_catalogs | 최상위 숫자만 비교; 중첩 계수 누락 |
| 비선형 세분 메시 검증 미완료 | 실패 상태·manifest | 보고서가 한계를 정확히 공개함 |

## 재현한 코드 문제

### 1. 중첩 재료 계수 충돌을 감사가 통과시킨다 — 우선 수정

material/materials_20261107.py:50–56은 공통 최상위 숫자와 law만 비교한다. Steel02/Concrete02의 params는 dict여서 비교 대상에서 빠진다. 임시 카탈로그 두 개에 동일 ID의 steel_s355_coupon을 두고 한쪽 params.Fy를 10% 바꿨을 때 audit_catalogs가 성공했다. 원본 카탈로그는 변경하지 않았다. 현재 실제 카탈로그의 충돌이 발견됐다는 뜻은 아니며, 충돌 방지 기능의 결함이다.

조치: 계수 dict를 재귀 비교하고, 모델별 필수 필드·허용 law를 명시한다. 중첩 값 변경 거부 검사를 추가한다. 출처도 덮어쓴 merged만이 아니라 각 원본 record에서 검증하는 편이 안전하다.

### 2. 요청한 선행 구간이 부족해도 완료된다 — 보고서 표현과 불일치

ground/ground_modeling/ground_20261107.py:72에서 first=max(0, ...)로 자른다. start=0.5초, warmup=5초, duration=1초를 넣으면 오류 없이 effective_warmup 약 0.5초가 된다. requested/effective 값은 JSON에 남으므로 완전히 숨겨진 변경은 아니지만, 보고서의 '범위를 조용히 자르는 대신 오류'는 선행 구간 전체에 적용되지 않는다.

조치: 기본은 선행 구간 부족 오류로 하고, 축소 허용은 명시 옵션으로 제공한다. duration/start/warmup의 0.1초 샘플 반올림도 요청값과 실효값을 모두 보존한다.

### 3. 비유한 수치 차단이 모든 원시 응답을 포괄하지 않는다 — 정적 검토

structure/mars_dome_ring_ground_20261107.py:269–281에서 전체 U, 절대가속도, 섬유 응답을 즉시 isfinite 검사하지 않는다. Python max는 NaN 위치에 따라 정상 숫자를 남길 수 있어 최종 JSON의 allow_nan=False만으로는 누락된 비정상 응답을 보장해 잡지 못한다. 지반 snapshot에는 별도 검사가 있으나 모든 구조 섬유·모든 시간의 검사는 아니다. 실제 실행에서 NaN이 발생했다는 증거는 없다.

조치: 매 단계 원시 U/acceleration/fiber/contact 값에 유한성 검사를 먼저 적용하고 실패 단계·절점·요소를 기록한다.

## 모델과 사용적합성의 핵심 한계

1. 비선형 spacing4m 사례는 완료되지 않았다. 탄성 세분 모델의 정적 성공은 비선형 메시 독립성을 보장하지 않는다. 지반 단독→링 결합→전체 조립 순으로 동일 입력·재료 조건을 유지하며 초기응력 경로와 국부 응력집중을 분리해야 한다.
2. 중심선 점 접촉과 폭×chord의 분담면적은 유한 폭 접촉압 분포를 해상하지 않는다. contact penalty와 마찰계수도 독립 보정·민감도 근거가 필요하다. 침투·전단 상대운동을 실제 지지력이나 영구 미끄러짐으로 해석하면 안 된다.
3. 지표 관측 입력을 기반에 적용하고 저면을 고정한다. 국부 점성 측면은 자유장·입사파 검증을 대체하지 않는다. 구조물 없는 파동 전달, 영역 크기·반사 민감도, 기반 입력 정합을 확인해야 한다.
4. 최신 통합 경로에는 구조 Rayleigh 감쇠 호출이 없다. 외곽 dashpot 및 구성식의 비탄성 효과와 구조 자체 감쇠를 구분하고, 감쇠 가정·에너지 수지·민감도를 명시해야 한다. 임의 5%를 넣는 것이 자동 해결책은 아니다.
5. S355 시편 모델의 SHS 전용, 명목 C30/SD400, 가정 힌지·배근, 지반 강성·반복 감쇠는 보정 범위가 제한된다. 내압 폐합, 복토, 접합부·좌굴·전단·정착·파괴 검증도 남았다.
6. 외부 출력 단계 극값과 첫 적분점 초기/최종 지반 snapshot은 시간·공간 전체 소성/손상 극값이 아니다. 강한 비선형 연구에는 하위 적분 단계 관측과 소성 일·잔류 변형·에너지 검사가 필요하다.

## 앞으로 가져갈 코드와 우선순위

- 기준 실행 경로: mars_dome_ring_ground_20261107.py + ground_20261107.py + materials_20261107.py. ring_foundation_rc.py와 capmat은 재사용하되 검증된 가정 범위를 유지한다.
- 이전 mars_dome_ring_ground.py 및 독립 실행기는 회귀 비교·보관용으로 구분한다. dome_in_mars_tuned.py의 과거 SAFE/minimum mass 출력은 최신 통합 모델의 검증된 설계 판정으로 승격하지 않는다.
- 첫 단계: 위 감사·선행 이력·유한성 문제 수정과 실패 경로 검사. Git 추적, 의존성 버전 고정, 상대경로 실행 안내, 검증 자동화를 마련한다.
- 다음 단계: 비선형 지반 초기응력 및 메시/시간/링 분할 수렴. 같은 입력 대역과 평가 구간으로 비교한다.
- 이후: 파동·경계 benchmark, 유한 폭 접촉, 실제 하중 폐합, 재료 및 부재 한계상태 검증. 이 근거가 쌓인 뒤 설계·최적 질량 판단으로 확장한다.

OpenSees 공식 근거: UniformExcitation은 상대 응답을 반환한다(https://opensees.github.io/OpenSeesDocumentation/user/manual/model/pattern/uniformExcitationPattern.html). Plain은 identity equalDOF와 동차 SP를 처리하며 MP 연쇄를 따라가지 않는다(https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/constraint/PlainConstraints.html). PDelta는 2차 효과를 고려하는 기하변환이다(https://opensees.github.io/OpenSeesDocumentation/user/manual/model/geomTransf/PDelta.html).
