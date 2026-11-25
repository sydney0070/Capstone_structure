# 지반·재료 보완과 화성/달 기초 연결부 초안

기준 코드: 기존 `mars_dome_ring_ground_20261107.py`. 신규 연구와 결과는 이 폴더에 분리했다. 원래 재료 카탈로그 계수는 바꾸지 않았다. 연구 논문의 기초 개념을 참고했으며, 논문의 구조물·치수·재료·실험 결과를 재현한 모델은 아니다.

현재 통합 실행이 채택한 재료값·위치·주요 코드 트리·실행 명령은 [현재 모델 안내](/Users/sydney/Desktop/capstone/Capstone_structure/docs/MODEL_OVERVIEW_20261109_KO.md)에 정리했다. 기존 RC 링 통합 경로와 별도 판기초 연구 경로를 구분한다.

현재 실행은 **지반 단독 `ground_20261107.py` / 통합본 `mars_ssi_integrated_20261109.py`**로 분리했다. 첨부 화면 형식으로 바꾼 HTML과 검증 결과·명령은 [시각화 및 실행 분리 보고서](/Users/sydney/Desktop/capstone/Capstone_structure/research_20261109/VISUALIZATION_SPLIT_20261109_KO.md)를 참고한다. 이전 통합 이름은 호환 진입점이다.

후속 핵심 문제 보완과 최신 검증 결과는 [후속 보고서](/Users/sydney/Desktop/capstone/Capstone_structure/research_20261109/FOLLOWUP_REPORT_KO.md) 및 [재현 안내](/Users/sydney/Desktop/capstone/Capstone_structure/research_20261109/followup_20261109/README.md)를 참고한다. 고정 접촉 경계 메시, 실제 지형 적분, 세부 시간 단계 최대값, 프로세스 정리를 추가했다. 비선형 메시 독립성과 대형 달 모델은 여전히 미해결이다.

반복 내부 잔차·무점착 완전소성 기작·실패 후 복원 진단과 현재 1107 계열의 사용 판단은 [반복 실패 원인 분석](/Users/sydney/Desktop/capstone/Capstone_structure/research_20261109/FAILURE_CAUSE_ANALYSIS_KO.md)에 기록했다. 연구 기준은 **1109 수정사항이 반영된 1107 계열**이며, 대형 모델까지 검증된 최종 설계 버전이라는 의미는 아니다.

## 바로 실행할 수 있는 초안

프로젝트 루트 `/Users/sydney/Desktop/capstone/Capstone_structure`에서 cap 환경으로 실행한다. 출력은 새로운 빈 폴더를 지정한다.

```bash
# 화성: 전면 판기초–돔, 상대 회전 고정
/Users/sydney/anaconda3/envs/cap/bin/python structure/planet_foundation_demo_20261109.py --body Mars --output-dir research_20261109/artifacts/my_mars

# 달: 환형 판기초–돔, 회전 힌지
/Users/sydney/anaconda3/envs/cap/bin/python structure/planet_foundation_demo_20261109.py --body Moon --rotation hinged --output-dir research_20261109/artifacts/my_moon

# 달의 상대 회전 고정 비교
/Users/sydney/anaconda3/envs/cap/bin/python structure/planet_foundation_demo_20261109.py --body Moon --rotation fixed --output-dir research_20261109/artifacts/my_moon_fixed

# 회전 스프링: 강성은 사용자가 근거를 갖고 명시해야 한다.
# --rotation spring --rotational-stiffness <Nm/rad>
```

기본값은 돔 직경 **4m**, 강재 목표 질량 **0.2t**, 판 두께 **0.15m**, 환형 폭 **1.2m**, 최대 수평 간격 **0.75m**, 지반 깊이 **3m**의 소형 실행 예제다. 이 치수는 추천 설계치가 아니다. 기존 40m 돔과 혼동하지 않는다. 달에 화성 지진파를 자동 적용하지 않으며, 기본 실행은 자중 정적 SSI다.

직경 40m의 시도도 결과에 남겼으나 완료되지 않은 사례가 있다. `RESEARCH_REPORT_KO.md`의 상태표를 먼저 확인한다. `--max-runtime` 기본 180초이며 원시 solver도 별도 worker 종료로 제한한다.

## 모델 구성과 한계

- `structure/planet_foundations_20261109.py`: 면적을 해상한 ShellMITC4 전면/환형 판 + 각 지반 절점의 면적 가중 압축/마찰/들림 접촉.
- `connect_dome_port`: 상대 회전 고정 / 힌지 / 명시한 회전 스프링. **고정은 기초와 돔 사이의 상대 고정**이고 지반 전체를 고정한다는 뜻은 아니다.
- 기초 재료는 균질 탄성 가정이다. E·ν·밀도가 출력에 명시된다. 화성/달 현장 보정 RC 또는 논문 속 지오폴리머로 인증된 재료가 아니다.
- 판 중립면과 접촉면을 같은 높이로 둔 이상화이며 두께 편심을 생략했다.
- 원형/환형 외곽은 직교 격자의 계단 모양이다. 목표 면적 대비 오차 10%를 넘으면 기본 거부한다. 이 예비 형상 gate가 형상·메시 독립성의 증명은 아니다.
- 셸 요소에는 밀도를 0으로 주고 질량과 자중을 절점에 한 번만 배분한다. 지반은 요소 재료의 밀도로 관성을 가진다.
- 기본 토층은 Drucker–Prager 비선형, 하부 암반 아날로그는 탄성이다. 달/화성 계수는 각 카탈로그를 사용한다.
- 앵커 옵션은 두 양의 입력 `--anchor-stiffness`, `--anchor-capacity`가 함께 있어야 한다. 압축을 전달하지 않는 인장 연결부 초안이며 매입 길이·부착·인발 파괴의 현장 검증 모델은 아니다. 기본적으로 사용하지 않는다.
- 내압·바닥의 압력 폐합, 복토 보유·분포, 기초 설치·굴착 이력, 전단/정착/접합부 파괴는 구현 완료로 표시하지 않는다.
- 실제 지형에 곧바로 수평 기초를 얹지 않는다. 별도 지반 단독 `site_ground_diagnostic.py`가 실제 ELYSE DEM 요철을 검사한다.

## 검증 및 연구 도구

```bash
/Users/sydney/anaconda3/envs/cap/bin/python structure/verify_revision_20261107.py
/Users/sydney/anaconda3/envs/cap/bin/python research_20261109/test_research.py
/Users/sydney/anaconda3/envs/cap/bin/python research_20261109/material_benchmarks.py
/Users/sydney/anaconda3/envs/cap/bin/python research_20261109/cyclic_soil.py
/Users/sydney/anaconda3/envs/cap/bin/python research_20261109/wave_column.py
```

결과를 덮어쓰지 않고 재실행하려면 각 연구 스크립트의 OUT 또는 사본 경로를 바꾼다. 기초 CLI와 통합 CLI는 비어 있지 않은 출력 폴더를 거부한다. 연구 benchmark의 정해진 출력은 재실행하면 해당 연구 출력만 갱신한다.

단계 진단은 `run_bounded_cases.py`로 감싼다. 예:

```bash
/Users/sydney/anaconda3/envs/cap/bin/python research_20261109/run_bounded_cases.py --timeout 90 --out research_20261109/artifacts/my_ground4 research_20261109/diagnose_stages.py --stage ground --spacing 4 --out research_20261109/artifacts/my_ground4
```

기존 통합 해석의 `--ground-history`는 모든 8개 적분점과 성공한 적응 세분 단계의 소성변형/일을 검사한다. `static_ground_history.csv`, `dynamic_ground_history.csv`는 매 성공 단계 저장한다. 지반 일 분해는 **전체 모델의 에너지 보존 검증이 아니다**. 현재 이 도구의 일 적분은 직교 직육면체 요소만 지원한다.

전체 결과와 미해결 항목은 `RESEARCH_REPORT_KO.md`, 기계 판독 요약은 `artifacts/research_manifest.json`에서 확인한다. 신규 형상 후보를 현장 안전 설계 또는 최적 질량으로 판정하지 않는다.

추가 초기화 실험 `dp_seating_probe.py`는 검증된 요소 `materialState` parameter를 사용한다. DP에 대한 `updateMaterialStage` wrapper는 현재 런타임에서 효과가 없었으며 관련 잘못된 해석은 거부 기록에 분리했다. 이 실험은 기본값이 아니며 일반 실행은 비선형 DP 경로를 유지한다.
