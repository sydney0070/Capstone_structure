# 구조물

- `dome_in_mars_tuned.py`: 현재 돔 형상·SHS 부재·강재 단면·질량
- `ring_foundation_rc.py`: C30/SD400 RC 링
- `planet_foundations_20261109.py`: 화성 전면 / 달 환형 판기초 형상·구조 요소
- `lunar_structure.py`: 기존 달 구조물 빌더
- `sections/rc.py`, `sections/elastic_rc.py`: RC 섬유 단면·탄성 등가 RC 구성
- 기존 `jjun*`, `dome*test`, `Rahmenbasetest`와 Abaqus 파일: 구조물 연구 원본

고정·힌지·회전 스프링 연결부는 `../SSI/connectors_20261109.py`를 사용한다. 지반과 구조물을 함께 실행하는 구현은 `../integrated/`에 있다. `mars_ssi_integrated_20261109.py`, `mars_dome_ring_ground*.py`, `planet_foundation_demo_20261109.py`의 옛 경로에는 호환 연결만 남겼다. 기존 구조물 단독 결과는 이 폴더에 보존하고, 통합 결과는 `../integrated/results/`로 옮겼다.
