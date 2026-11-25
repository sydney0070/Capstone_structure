# 지반 단독·통합본 분리 및 HTML 시각화 수정

작성일: 2026-11-09

## 변경 내용

앞선 실행 안내는 공용 지반 모듈을 설명하면서 통합 해석 명령을 제시해 지반 단독 실행과 혼동할 수 있었다. 실행 경로를 다음과 같이 분리했다.

| 구분 | 실행 파일 | 포함하는 모델 |
|---|---|---|
| 지반 단독 | `ground/ground_modeling/ground_20261107.py` | 비선형 토층·탄성 기반암·지반 자중·ELYSE 지진 입력·지반 경계 |
| 새 통합본 | `structure/mars_ssi_integrated_20261109.py` | 지반·C30/SD400 RC 링·S355 돔·기초 접촉 |
| 이전 통합 이름 | `structure/mars_dome_ring_ground_20261107.py` | 새 통합본으로 연결하는 호환 진입점 |

지반 단독 실행에는 돔·RC 링·기초 판·기초 접촉을 조립하지 않는다. 지반 경계의 점성 요소는 포함할 수 있다. 실행 옵션에서도 구조물 질량·링 분할·기초 접촉 선택을 제외했다. 기존 공용 수치 함수와 연구용 라이브러리는 의존성을 위해 유지했다. 과거 지반·SSI 예제 파일을 전부 삭제하거나 다른 모델의 물리 설정을 바꾸는 작업은 하지 않았다.

통합 구현은 새 파일에 실제로 옮겼으며 이전 이름은 모듈 정체성을 유지하는 호환 연결로 바꿨다. 기존 연구 스크립트의 import·patch와 이전 CLI 명령은 계속 동작한다. 새 결과는 두 경로의 서로 다른 결과 폴더에 저장된다.

## 시각화

첨부 화면에 쓰였던 `insight_dynamic_viewer.py`의 UI를 현재 결과에 연결했다.

- 왼쪽 설정 패널과 지층별 색상·재료 체크박스
- 동적 시작 이후 / 중력 포함 전체 변위 선택
- 배율 입력, 실제 크기 ×1, 자동 확대
- 재생·일시정지·처음·시간 슬라이더·재생 속도
- 최대 변위 프레임 이동과 저장 프레임별 최대 변위 그래프
- 절점에 올리면 확대 전 실제 변위를 mm로 표시
- 현재 ELYSE 위치·경계·층서 가정 표시

통합본은 같은 UI에 RC 링과 강재 돔 표시 선택을 추가했다. 지반 단독 HTML에는 이 구조물 그룹 자체가 없다. 두 경로 모두 실제 저장 총변위를 넣고 동적 기준에서는 첫 프레임 변위를 빼도록 연결했다. 배율은 화면 좌표에만 적용한다.

재료가 다른 층 사이의 면도 각 재료의 경계로 유지하도록 수정했다. 윗층을 숨겼을 때 아래층의 윗면이 사라지는 문제를 줄인다. 같은 재료끼리의 내부 면은 계속 제거한다. 표시 재료값은 현재 카탈로그를 읽으며 예전 캡처의 다른 물성값을 복사하지 않았다.

HTML은 Plotly를 포함한 오프라인 파일이다. 같은 이름의 JSON에는 실제 좌표·변위·시간·재료·표시 그룹을 저장했다. HTML 제목·설명도 지반 단독과 통합본을 구분한다. 단일 정적 결과에서는 동적 재생을 비활성화한다.

## 실제 실행 확인

| 실행 | 결과 |
|---|---|
| 지반 단독 `--quick --ground-history` | 완료. 지반 1,176개 요소, 1,575개 절점. 구조물·기초·접촉 개수 각각 0. 41개 프레임 |
| 새 통합본 `--quick --ground-history` | 완료. 41개 프레임, 지반 재료 3그룹＋RC 링·돔 2그룹 |
| 기존 통합 실행과 새 통합 실행의 응답 비교 | 동일 옵션의 `coupled_response.csv` 40×8 배열 차이 0 |

지반 단독은 조립 직후 모든 절점이 3자유도이며 요소 태그 집합이 지반 요소 집합과 같은지 확인한다. 이 검사는 점성 경계를 추가하기 전에 수행한다. 새 경로도 비선형 DP를 유지하며 재료 계수·허용오차를 수렴 목적으로 바꾸지 않았다.

검증 결과 경로:

- [지반 HTML](/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/results_ground_20261109/viewer_split_verified_v2/ground_response.html)
- [통합 HTML](/Users/sydney/Desktop/capstone/Capstone_structure/structure/results_mars_ssi_integrated_20261109/viewer_split_verified/mars_integrated_response.html)
- 검사 로그·제어 코드: `research_20261109/viewer_split_20261109/`
- 수정 전 소스: 같은 폴더의 `before/`

지반 단독 첫 개발 실행은 NumPy 논리값을 JSON으로 직접 저장해 실패했다. 공용 직렬화 함수를 적용한 뒤 새 폴더의 두 번째 실행이 완료됐다. 첫 실패 결과를 성공으로 덮어쓰지 않았다.

## 검사 범위와 한계

추가 Python 검사 4개에서 재료 경계면, 정렬된 절점에 대한 변위 대응, 총변위·증분변위 보존, 스크립트 문자열 이스케이프, 이전 통합 모듈 호환, 지반 CLI 옵션 분리를 확인했다. 기존 회귀 11개·연구 9개·후속 18개도 통과해 **총 42개 검사를 통과했다.** 최종 결과는 이 작업 폴더의 로그와 `verification.json`에 남겼다.

실제 생성된 두 HTML의 JavaScript를 DOM·Plotly 대체 객체와 실행해 변위 기준, 배율에 대한 실제 변위값 불변, 재료·구조물 숨김, 시간 이동, 최대 프레임, 재생·일시정지 동작을 확인했다. 이는 제어 로직·데이터 검사이며 실제 WebGL 렌더링 검사가 아니다.

인앱 브라우저는 `file://` 로컬 HTML 열기를 보안 정책으로 거부했다. 다른 브라우저 경로나 우회 서버로 이를 우회하지 않았다. 따라서 실제 Chrome의 최종 화면 배치·WebGL 렌더링은 직접 확인하지 못했다. 기존 화면의 템플릿을 재사용했고 생성 파일은 사용자가 Chrome에서 열 수 있도록 제공한다.

이 작업은 실행 분리와 시각화를 다룬다. 40 m 달 모델 수렴, 비선형 메시 독립성, 현장 물성 보정 등의 기존 미완료 문제를 해결한 것으로 표시하지 않는다.

## 실행 명령

지반만 실행:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling/ground_20261107.py \
  --quick --ground-history
```

통합본 실행:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python \
  /Users/sydney/Desktop/capstone/Capstone_structure/structure/mars_ssi_integrated_20261109.py \
  --quick --ground-history
```

`--quick`는 선행 1초＋응답 1초의 연결 점검이다. 지반 단독 결과는 `ground/ground_modeling/results_ground_20261109/`, 통합본 결과는 `structure/results_mars_ssi_integrated_20261109/`의 새 실행 폴더에 저장된다. 각각 `ground_response.html`, `mars_integrated_response.html`을 연다. 과거 HTML은 수정 전 결과로 보존돼 있으므로 새로 생성한 파일을 열어야 한다.

## 추가 수정: 세 지층의 색상 구분

기존 갈색 계열이 비슷하게 보인다는 피드백에 따라 경화층은 황금색 `#E6A23C`, JEZ-1은 하늘색 `#63B5D6`, 기반암은 짙은 회색 `#596579`으로 변경했다. 초기 형상의 재료별 면 범위는 각각 z=-0.2–0, -8–-0.2, -12–-8 m로 확인했다. 첫 층 두께 0.2 m는 그대로이며 보기 편하게 층 두께를 바꾸지 않았다.

앞서 제공한 두 HTML·표시 JSON에도 색상을 반영했다. 갱신 전 파일은 `viewer_split_20261109/layer_colors/`에 보존했고, 색상 이외 표시 데이터가 전부 같은지 비교했다. 좌표·변위·시간·재료값·삼각형·구조물 데이터는 변경하지 않았으며 제어 로직 검사도 두 파일에서 다시 통과했다. 해석은 재실행하지 않았다. 새 기본 팔레트는 이후 생성하는 HTML에도 적용된다.
