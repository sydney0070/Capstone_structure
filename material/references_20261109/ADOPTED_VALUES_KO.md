# 현재 채택값과 근거

기준: 2026-11-09의 실제 카탈로그·코드. 단위는 별도 표시가 없으면 m·N·s·kg·Pa이다. **문헌에 있는 숫자, 그 숫자의 모델 전용, 자체 가정, 실험 보정은 서로 다른 근거 수준**이다. 이 구분은 현재 코드의 기록을 추적한 것이며 새로운 물성 재보정 결과가 아니다.

## 1. 현재 화성 지반 단독 실행

진입점은 프로젝트의 `ground/ground_modeling/ground_20261107.py`, 실행 구현은 `ground_only_20261109.py`이다. 기본 지반은 60×60×12 m, 수평 간격 5 m, 비선형 자중 20단계 후 지진 입력이다. 카탈로그와 기본 행성 설정 위에 실행 코드의 크기·메시·렌즈 제거 설정이 덮어써진다. 따라서 `Mars.json`의 작은 예제 영역을 현재 실행 영역으로 읽으면 안 된다.

| 층 / 재료 ID | 채택값 | 참고자료 | 채택 성격·한계 |
|---|---|---|---|
| 0–0.2 m / `mars_insight_duricrust_proxy` | DP; E=11.48 MPa, ν=0.28, ρ=1,200 kg/m³, φ=30°, c=5.8 kPa, ψ=0° | [Spohn 등, 2022 PDF](01_ground/spohn_2022_insight_hp3.pdf), [출판사](https://link.springer.com/article/10.1007/s11214-022-00941-z); 카탈로그의 §6.4–6.6 연결 | 문헌 기반 proxy. E·ν는 얕은 지반의 파동 기반 추정값을 표층에 전용. 정적 duricrust 단독 E의 직접 측정값이 아님. c는 가정 φ·ρ에서 추론한 값. ψ=0은 자체 가정. |
| 0.2–8 m / `mars_jez1_peak` | DP; E=30 MPa, ν=0.30, ρ=1,840 kg/m³, φ=40.6°, c=3.1 kPa, ψ=0° | [JEZ-1 논문 PDF](01_ground/jez1_2021.pdf), [추출값](01_ground/jez1_inputs/published_inputs.json), [DOI](https://doi.org/10.3390/app11157052); Table 4, Fig.10, §4.3 | 밀도·강도는 지구에서 시험한 시뮬런트 자료의 전용. **E=30 MPa·ν=0.30·ψ=0은 가정**. ELYSE 해당 깊이의 직접 조사값이 아님. |
| 8–12 m / `basalt_elastic_analogue` | 탄성; E=40 GPa, ν=0.30, ρ=2,900 kg/m³ | [Musiol 등, 2016 PDF](01_ground/musiol_2016_basalt.pdf), [DOI](https://doi.org/10.1002/2015JE004896); Table 1의 지각 모델 값 | 발표된 모델의 탄성 상수를 국부 기반암 유사재료로 전용. ELYSE 시추 코어 시험값이 아님. 균열·절리·암반 파괴 미반영. |

층 경계 0.2 m·8 m, 전체 깊이 12 m, 평탄한 지표 및 층의 결합은 이 코드의 시나리오다. 서로 다른 문헌을 결합한 층서이며 하나의 관측 지층으로 검증되지 않았다. 토층은 DP 비선형을 유지하고, 기반암만 탄성으로 둔다.

두 토층의 ψ=0, 경화 인자 0 및 일정한 탄성 강성은 구성 모델의 선택이다. 문헌 강도를 DP로 변환했다고 반복감쇠·체적변형·저구속압 응답까지 실험 보정되는 것은 아니다. [DP 공식 명령 문서 사본](04_model_docs/DruckerPrager.html), [기존 변환식 검토](04_model_docs/dp_mapping_note.md), [실행 엔진 3.4.0 구현 사본](04_model_docs/DruckerPrager_v3.4.0.cpp)을 함께 둔다.

## 2. 현재 화성 통합 모델의 구조 재료

아래 구조 재료는 `integrated/mars_ssi_integrated_20261109.py`의 통합 실행에 해당한다. **지반 단독 실행에는 구조물·RC 링·접촉이 없다.**

| 대상 | 채택값 | 참고자료 | 채택 성격·한계 |
|---|---|---|---|
| 돔 / `steel_s355_coupon` | Steel02; E=204.167883 GPa, Fy=384.169992 MPa, b=0.0052586813, ρ=7,850 kg/m³, ν=0.30 | [S355 데이터 DOI](https://doi.org/10.5281/zenodo.6965147), [원본 출처·시험 조건](02_structure/s355/provenance.json), [보정 결과](02_structure/s355/RESULTS_KO.md), [실제 채택 계수](02_structure/s355/reviewed_config.json) | S355 J2+M HEM320 웹 선택 시편의 상온·파단 전 구간으로 보정. SHS 단면에 전용한 상태. ρ·ν는 시편에서 식별한 값이 아닌 관례 가정. |
| RC 링 콘크리트 / `concrete_c30_prototype` | Concrete02; fpc=−30 MPa, epsc0=−0.002, fpcu=−6 MPa, epsU=−0.006, λ=0.1, ft=3 MPa, Ets=1.5 GPa; E=30 GPa, ν=0.20, ρ=2,400 kg/m³ | [실제 카탈로그](00_current_values/release_catalog.json), [Concrete02 공식 문서](04_model_docs/Concrete02.html) | **고정된 예제 가정값**. 공식 문서는 법칙·인자의 의미를 설명하며 이 콘크리트 계수의 시험 근거는 아님. |
| RC 링 철근 / `rebar_sd400_prototype` | Steel02; Fy=400 MPa, E=200 GPa, b=0.01, R0=18, cR1=0.925, cR2=0.15, ρ=7,850 kg/m³ | [실제 카탈로그](00_current_values/release_catalog.json), [Steel02 공식 문서](04_model_docs/Steel02.html) | **명목 예제 가정값**. 이전 UCSD A6 재보정 계수와 다르며 해당 데이터를 SD400 검증 근거로 연결하지 않는다. |
| 외피 / ETFE | ρ=1,750 kg/m³, 1겹 0.2 mm×3겹 → 면밀도 1.05 kg/m² | [값이 선언된 돔 코드 사본](00_current_values/code_snapshot/structure/dome_in_mars_tuned.py) | 기존 설계 가정. 직접 연결된 제품 시험자료는 확보하지 못했다. 통합 해석에는 질량만 반영하며 막 강성·내압 폐합은 해석하지 않는다. |
| 지반–RC 링 접촉 | μ=0.35, 부착력 0, penalty factor=10 | [현재 통합 코드 사본](00_current_values/code_snapshot/integrated/mars_ssi_integrated_20261109.py), [접촉 계수 계산](00_current_values/code_snapshot/SSI/parameters.py) | 마찰계수와 penalty factor는 자체 가정. JEZ-1 매끈한 강판 시험의 μ≈0.3424·부착력 1.5 kPa를 그대로 채택한 것이 아니다. RC–토양 접촉 보정 미완료. |

통합 모델의 돔 강재 목표 질량 200 t, 링 반지름 20 m·폭 1 m·깊이 0.6 m, 기본 12개 링 요소/접촉 위치, 피복 0.05 m·철근 지름 0.02 m·상하 각 4개, 소성힌지 길이 0.15 m도 연구용 설계 입력이다. 문헌에서 검증된 최적 치수가 아니다. 단독 RC 링 빌더의 기본 분할 48개와 통합 실행의 12개를 구분한다.

## 3. 달 기초 초안의 지반·기초값

달 예제는 현재 화성 지반 CLI와 다른 실행이다. `Moon.json`을 사용하는 기본 층별 JSC-1A peak 상태를 아래에 적었다. `moon_field_proxy_*`는 NASA 점착력을 조합한 **별도 후보**이며 기본 무점착 JSC-1A의 보정값을 덮어쓴 것이 아니다.

| 재료 ID / 기본 깊이 | E [MPa] | ρ [kg/m³] | φ [°] | c [Pa] | 근거 |
|---|---:|---:|---:|---:|---|
| `moon_jsc1a_17_peak` / 0–0.3 m | 15.10 | 1,705 | 42.9720 | 0 | LASP 삼축시험 요약의 밀도군별 강도 적합·30 kPa군 정적 접선 E 요약 |
| `moon_jsc1a_18_peak` / 0.3–1 m | 20.45 | 1,780 | 45.5722 | 0 | 동일 |
| `moon_jsc1a_19_peak` / 1–8 m | 17.25 | 1,880 | 49.4849 | 0 | 동일 |
| `basalt_elastic_analogue` / 8 m 이하 | 40,000 | 2,900 | — | — | 위 현무암 유사재료 |

모든 JSC-1A군에 ν=0.30·ψ=0을 가정한다. [LASP 원문](01_ground/jsc1a_lasp/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf), [46개 시험 요약](01_ground/jsc1a_lasp/triaxial_observations.json), [추출 조건·한계](01_ground/jsc1a_lasp/provenance.json), [실제 강도 적합 결과](01_ground/moon_strength_fits.json)를 함께 확인한다. 원 시험 구속압은 15–60 kPa이며 달 얕은 지반의 더 낮은 압력으로 전용하는 한계가 남는다.

`moon_field_proxy_*`의 점착력 520/900/3,000 Pa는 [NASA 2024 가이드](01_ground/nasa_lunar_simulant_guide_2024.pdf)의 Table 18, 보고서 p.48 평균값 연결이다. 원문 깊이 범위와 현재 두꺼운 모델 층은 같지 않아 층별 배정은 추가 가정이다. [Connolly–Carrier 깊이별 물성 초록](01_ground/connolly_carrier_2023_lunar_profile.pdf)은 물성이 깊이에 따라 달라진다는 배경 근거이며 현재 층 경계의 직접 근거는 아니다.

화성 전면판·달 환형판 초안은 기본 E=30 GPa·ν=0.20·ρ=2,400 kg/m³·두께 0.6 m인 **명목 균질 탄성 기초**다. 달 콘크리트의 확정 재료 선정을 의미하지 않는다. 고정/힌지는 연결 이상화이고, 회전스프링 및 앵커의 강성·용량은 사용자가 명시해야 하는 입력으로 남아 있다. [초안 설정 코드](00_current_values/code_snapshot/structure/planet_foundations_20261109.py), [연결부 코드](00_current_values/code_snapshot/SSI/connectors_20261109.py).

## 4. 후보 재료와 이전 보정값

소결 JSC-1A, Al6061-T6, AS4/8552, 명목 S355, 지구 모래, 차폐용 느슨한 JSC-1A, peak/residual·field proxy 등은 두 카탈로그에 등록되어 있으나 모두 현재 기본 돔이나 지반에 동시에 사용되는 것은 아니다.

각 후보의 수치·law·status·assumptions·sources는 [29종 전체 기록](00_current_values/all_material_records.json)에 정확히 보존했다. [자료 목록](SOURCES_KO.md)에서 [소결체 원문](03_candidates/sintered_jsc1a_leag2016.pdf), [Hydro 6061 원문](03_candidates/hydro_6061_2019.pdf), [NIST 6061 자료](03_candidates/nist_6061_t6.html), [AS4/8552 원문 링크](https://www.wichita.edu/industry_and_defense/NIAR/Research/hexcel-8552/AS4-Unitape-2.pdf)를 찾을 수 있다. 이전 6014 실험은 6061의 근거로, 이전 Smartfan 횡인장 자료는 AS4 전체 물성의 근거로 바꾸어 연결하지 않는다.

## 5. 위치·입력 및 수치 설정

현재 화성 위치는 [ELYSE StationXML 사본](06_site_and_input/stationxml.xb.elyse.2019.151.xml)의 **위도 4.502384°N, 경도 135.623447°E**이며 좌표는 East–North–Up이다. DEM 높이를 기준으로 한 평면 지반이고 실제 DEM 기복을 현재 형상에 적용한 것은 아니다.

지진은 `S1222a_ACC.mseed`에서 처리한 ELYSE 관측 가속도를 사용한다. [처리 메타데이터](06_site_and_input/motion_full_enu.json), [선택 구간 메타데이터](06_site_and_input/motion_window_enu.json), [해당 완료 실행 맥락](06_site_and_input/ground_run_context.json)에 단위·처리·원본 해시·위치를 보관했다. 관측 지표 기록을 기저 입력으로 사용하는 가정이 남아 있다.

중력 Mars=3.711 m/s²·Moon=1.62 m/s²는 시나리오 설정이다. 지반 실행 기본 dt=0.05 s·상한 주파수 1 Hz·수평 간격 5 m·자중 20단계, 기본 선행 5초＋평가 2초, `--quick`의 선행 1초＋평가 1초는 수치·시간 구간 설정이며 문헌 물성이나 메시 독립성 입증값이 아니다. HTML 변형 확대배율은 표시용이다.
