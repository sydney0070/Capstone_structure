# 독립 RC 링 기초 코드

실행 파일: `ring_foundation_rc.py`. 원본 황콘크리트 돔과 튜닝 강재 돔 파일은 변경하지 않았다.

## 선택한 재료

기존 `material_release_20261101/configs/catalog.json`의 `concrete_c30_prototype` 및 `rebar_sd400_prototype`을 사용한다. 각각 Concrete02와 Steel02다. 카탈로그의 명목 기준 재료이며 실험곡선 보정 완료 재료가 아니다. 새 파일은 황콘크리트와 다른 RC 비교안이다.

## 예제 형상과 해석

- 반지름 20 m, 48개 직선 보로 근사한 폐합 링, 직사각 단면 폭 1 m × 높이 0.6 m.
- 원본의 역T형 황콘크리트 단면을 복제한 모델이 아니다.
- 상·하부 각 4개, 직경 20 mm의 원주방향 철근. 피복은 철근 표면까지 50 mm. 이는 예제 배근이며 설계·정착·이음·스터럽/구속 검토 결과가 아니다.
- 기존 패키지의 철근 면적 중복을 제거한 RC 섬유단면을 재사용한다. HingeRadau 적분, 단부 비선형 섬유단면 및 중앙 변환 탄성단면, 탄성 비틀림 GJ를 사용한다. 단부 힌지 길이 0.15 m는 명시적 가정이며 물리적 보정값이 아니다. 요소 길이에 맞지 않는 힌지 설정은 차단한다.
- 별도 예제 실행에만 압축 전용 ENT 수직 스프링과 탄성 수평 스프링을 사용한다. k=20 MN/m³ 및 수평/수직 비 0.5는 가정이다. 실제 3D 지반이 아니며 마찰 미끄러짐·인발 저항은 예제에서 모델링하지 않는다.
- 원주상의 12개 연결점에 각각 100 kN 하향 예제 하중 및 링 자중을 준다. 실제 돔 반력을 추출한 하중은 아니다.
- 질량 및 자중은 모델에 사용한 다각형의 실제 현 길이에 따라 계산한다. 보 요소 질량과 절점 질량을 중복 부여하지 않는다.
- 균열/콘크리트 연화와 강재 비선형 구성식은 연결되어 있지만, 이번 작은 하중 실행이 큰 손상 상태를 검증하지는 않는다. 전단 파괴·철근 좌굴·정착·연화 객관성·큰 변형·전체 설계 안전성은 미검증이다.

## 실행

```bash
/Users/sydney/anaconda3/envs/cap/bin/python /Users/sydney/Desktop/capstone/Capstone_structure/structure/ring_foundation_rc.py
```

기본 출력: `structure/results_ring_rc/`. `ring_rc.html`, `ring_displacements.csv`, `ring_local_forces.csv`, `summary.json`, `run_status.json`을 생성한다. HTML은 실제 정적 계산 변위의 100배 표시다. 동적 애니메이션은 아직 포함하지 않는다. 완료는 SAFE 판정이 아니다.

`--radius`, `--segments`, `--width`, `--depth`, `--bar-diameter`, `--bars-per-face`, `--hinge-length`, `--vertical-load-per-port`로 예제 설정을 바꾼다. `--elastic`은 동일 배근·물성의 탄성 비교이며 `--no-html`은 화면 출력을 생략한다. 결과 보존에는 `--output-dir`을 사용한다.

## 통합 API

```python
from ring_foundation_rc import RingConfig, build_ring
ring = build_ring(RingConfig(tag_base=600000), demo_supports=False)
# ring['nodes']: 6자유도 링 절점
# ring['dome_ports']: 30도 간격의 12개 구조물 연결점
# ring['nodal_mass_kg']: 이미 관성 질량은 부여됨; 자중은 호출자가 한 번 부여
```

생성 함수는 `ops.wipe()`를 호출하지 않는다. 지반·구조물과 합칠 때에는 전체 모델 소유자가 초기화와 자중 재하를 담당한다. 기본 API는 스프링·하중을 만들지 않는다. 링은 동/북/상 좌표의 원점 중심이며, 실제 돔의 기초 좌표와 축·위치를 맞춘 뒤 연결해야 한다. 3DOF 지반과 6DOF 링은 접촉 및 병진 하중 전달을 따로 구현해야 하며 회전을 임의로 같게 묶지 않는다. 돔과 링에 같은 물리적 질량·자중을 두 번 부여하지 않는다.

## 이번 확인 결과

- 자중 + 총 1.2 MN 예제 하중, 20단계 정적 해석 완료.
- 총 수직 하중 1.87725 MN, 최대 침하 0.974936 mm.
- 힘 평형 상대오차 약 4.96e-16, 모멘트 평형 정규화 오차 약 1.31e-16.
- 동일 탄성 비교 최대 침하 0.974847 mm, 비선형 결과와 약 0.00916% 차이. 작은 하중에서의 일관성 확인이며 비선형 손상 검증이 아니다.
- 기존 임의 절점 보존, 지반 스프링 기본 미생성, 12개 연결 포트 및 반복 생성 태그 충돌 차단 확인.
- 메시·힌지·배근 민감도, 실제 돔 하중 및 실제 지반과의 통합은 아직 하지 않았다.
