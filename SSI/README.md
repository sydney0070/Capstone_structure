# SSI · 연결부

| 파일 | 역할 |
|---|---|
| `contact.py` | 마찰·들림 접촉과 점성 경계 요소 등록 |
| `parameters.py` | 면적·인접 강성 기반 접촉 penalty와 방사 감쇠 계수 |
| `interface.py` | 지반–기초 조립·반력·접촉 응답·하중 전달 |
| `connectors_20261109.py` | 기초–돔 고정·힌지·회전 스프링 및 인장 타이 연결 |

이 폴더가 현재 SSI 구현의 기준이다. 재료 카탈로그는 `material/`, 구조물·기초 형상은 `structure/`, 전체 해석은 `integrated/`를 사용한다. 이전 `capmat.ssi`, `capmat.mechanics`의 SSI 함수, `ground_hybrid.interface`는 이 구현으로 연결한다. 중복된 실제 접촉 구현을 두지 않는다.
