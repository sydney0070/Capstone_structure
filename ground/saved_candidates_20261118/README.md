# 저장 후보: MGS-1 · MGS-1S

2026-11-18 사용자 선택에 따라 두 후보의 초기 물성, 시편 시험과 완료된 3D 지반 결과를 복사해 보관했다. 이 저장본의 선택 후보는 두 개다.

| 후보 | 밀도 (kg/m³) | E (MPa) | ν | φ | c (Pa) | 감압 최대 침하 (mm) | 가압 최대 침하 (mm) |
|---|---:|---:|---:|---:|---:|---:|---:|
| MGS-1 | 1300 | 28.4 | 0.22 | 35.5° | 280 | 1.433758 | 0.282701 |
| MGS-1S | 1300 | 28.4 | 0.22 | 34.0° | 280 | 1.433890 | 0.282829 |

두 후보는 ψ=0°, 무경화 Drucker–Prager 초기 재료다. 현장 보정을 완료한 물성은 아니다. 삼축압축 및 직접전단 수치 곡선과 물성 출처를 함께 저장했다. 원본으로부터 복사한 파일의 SHA-256은 manifest.json에서 확인할 수 있다.

## 파일

- [선택 재료 카탈로그](selected_materials.json)
- [저장 목록과 원본 대조 해시](manifest.json)
- [두 후보의 단순전단 대조 수치](shear_comparison.csv)
- [MGS-1 물성](mgs1/material.json) · [시편](mgs1/specimen.json) · [완료 해석 요약](mgs1/ground_analysis/summary.json) · [3D 침하 보기](mgs1/ground_analysis/footing_settlement.html)
- [MGS-1S 물성](mgs1s/material.json) · [시편](mgs1s/specimen.json) · [완료 해석 요약](mgs1s/ground_analysis/summary.json) · [3D 침하 보기](mgs1s/ground_analysis/footing_settlement.html)

각 후보의 ground_analysis/에는 원본 summary.json, HTML/JSON 응답과 수렴 로그를 보관했다. model_checkpoint.json은 실행 시작 당시의 체크포인트이므로 완료 여부는 summary.json과 execution_status.json을 기준으로 본다. 3D 뷰어는 감압 상태의 자중 평형 이후 증분 변위를 표시한다.

## 저장된 지반의 현재 구성

영역은 90 × 90 m, 깊이 30 m이고 기초 밑면은 깊이 3 m다.

| 깊이 | 재료 | 구성 모델 |
|---|---|---|
| 0–0.2 m | 천부 경화층 | Drucker–Prager 탄소성 |
| 0.2–3 m | 선택한 MGS-1 또는 MGS-1S 후보 | Drucker–Prager 탄소성 |
| 3–13 m | 암괴가 섞인 조립질 레골리스 | ElasticIsotropic 탄성 |
| 13–30 m | 파쇄 현무암 | ElasticIsotropic 탄성 |

구성 모델이 탄소성이라는 사실과 실제 재하 중 소성 변형이 발생했다는 판단은 구분해야 한다. 현재 저장본은 기존 해석을 그대로 보관한 것으로 하부층 모델 변경이나 추가 재해석을 수행하지 않았다.

## 탄소성 확장 검토

기초 아래의 잔류 침하와 항복 가능성을 보려면 우선 3–13 m 조립질층을 탄소성으로 확장한 별도 모델을 현재 결과와 비교하는 접근을 추천한다. 이 층의 c, φ, 팽창각 및 경화 설정에는 층에 맞는 시험·문헌 근거가 필요하다. 상부 저밀도 모사토의 강도를 하부층에 그대로 복사하지 않는다. 13–30 m 파쇄 현무암은 우선 기존 탄성 기준으로 두고 하부층의 응력과 변형을 확인한 뒤 파괴 모델의 필요성을 판단한다.

[OpenSees Drucker–Prager 문서](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html)는 항복·경화·소성 체적변화 파라미터 및 소성 상태 기록을 설명한다. 탄소성 전환만으로 지지력·붕괴·접촉 문제의 검증이 완료되는 것은 아니다.
