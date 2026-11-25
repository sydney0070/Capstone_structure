# 참고자료 목록

정리일: 2026-11-09 (Asia/Seoul). 원 카탈로그의 출처 연결을 보존했다. 이 목록은 모든 수치를 새로 검증했다는 뜻이 아니다.

| 출처 ID | 사용 재료 / 목적 | 원문 | 모은 파일 |
|---|---|---|---|
| `s355_coupon` | steel_s355_coupon | [공식 출처 / DOI](https://doi.org/10.5281/zenodo.6965147) | [provenance.json](02_structure/s355/provenance.json)<br>[reviewed_config.json](02_structure/s355/reviewed_config.json) |
| `hydro6061` | al6061_t6_extrusion | [공식 출처 / DOI](https://www.hydro.com/globalassets/01-products--services/extruded-profiles/americas/ena-resources/alloy-data-sheets/hydro_2019_data_sheet_6061.pdf) | [hydro_6061_2019.pdf](03_candidates/hydro_6061_2019.pdf) |
| `nist6061` | al6061_t6_extrusion | [공식 출처 / DOI](https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061) | [nist_6061_t6.html](03_candidates/nist_6061_t6.html) |
| `niar_as4` | cfrp_as4_8552_ud_rtd | [공식 출처 / DOI](https://www.wichita.edu/industry_and_defense/NIAR/Research/hexcel-8552/AS4-Unitape-2.pdf) | 원문 링크만 보관; 로컬 원문 미확보 |
| `sintered2016` | sintered_jsc1a_high_porosity, sintered_jsc1a_low_porosity | [공식 출처 / DOI](https://www.hou.usra.edu/meetings/leag2016/eposter/5082.pdf) | [sintered_jsc1a_leag2016.pdf](03_candidates/sintered_jsc1a_leag2016.pdf) |
| `lasp_jsc1a` | moon_field_proxy_17, moon_field_proxy_18, moon_field_proxy_19, moon_jsc1a_17_peak, moon_jsc1a_17_residual, moon_jsc1a_18_peak, moon_jsc1a_18_residual, moon_jsc1a_19_peak, moon_jsc1a_19_residual, shield_loose_jsc1a | [공식 출처 / DOI](https://ares.jsc.nasa.gov/projects/simulants/attachments/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf) | [JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf](01_ground/jsc1a_lasp/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf)<br>[triaxial_observations.json](01_ground/jsc1a_lasp/triaxial_observations.json)<br>[moon_strength_fits.json](01_ground/moon_strength_fits.json) |
| `jez1` | mars_jez1_peak, mars_jez1_residual | [공식 출처 / DOI](https://www.mdpi.com/2076-3417/11/15/7052) | [jez1_2021.pdf](01_ground/jez1_2021.pdf)<br>[published_inputs.json](01_ground/jez1_inputs/published_inputs.json) |
| `pdmy02` | earth_sand_dr50_dp_reference, earth_sand_dr50_dry_reference, earth_sand_dr50_pdmy02 | [공식 출처 / DOI](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html) | [PressureDependentMultiYield02.html](04_model_docs/PressureDependentMultiYield02.html) |
| `dp` | mars_jez1_peak, mars_jez1_residual, moon_field_proxy_17, moon_field_proxy_18, moon_field_proxy_19, moon_jsc1a_17_peak, moon_jsc1a_17_residual, moon_jsc1a_18_peak, moon_jsc1a_18_residual, moon_jsc1a_19_peak, moon_jsc1a_19_residual | [공식 출처 / DOI](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html) | [DruckerPrager.html](04_model_docs/DruckerPrager.html) |
| `pdmy02_table` | earth_sand_dr30_dp_reference, earth_sand_dr30_dry_reference, earth_sand_dr75_dp_reference, earth_sand_dr75_dry_reference | [공식 출처 / DOI](https://opensees.berkeley.edu/wiki/index.php/PressureDependMultiYield02_Material) | 원문 링크만 보관; 로컬 원문 미확보 |
| `musiol2016` | basalt_elastic_analogue | [공식 출처 / DOI](https://gfzpublic.gfz.de/rest/items/item_1802896_4/component/file_1802898/content) | [musiol_2016_basalt.pdf](01_ground/musiol_2016_basalt.pdf) |
| `spohn2022` | mars_insight_duricrust_proxy | [공식 출처 / DOI](https://link.springer.com/article/10.1007/s11214-022-00941-z) | [spohn_2022_insight_hp3.pdf](01_ground/spohn_2022_insight_hp3.pdf) |
| `lunar_profile` | 구성 법칙·접촉 구현 참고 | [공식 출처 / DOI](https://ntrs.nasa.gov/api/citations/20220018853/downloads/LPSC%202023%20Abstract_Connolly_Carrier_v2.pdf) | [connolly_carrier_2023_lunar_profile.pdf](01_ground/connolly_carrier_2023_lunar_profile.pdf) |
| `layered_example` | 구성 법칙·접촉 구현 참고 | [공식 출처 / DOI](https://opensees.ist.berkeley.edu/wiki/index.php?title=Effective_Stress_Site_Response_Analysis_of_a_Layered_Soil_Column) | 원문 링크만 보관; 로컬 원문 미확보 |
| `nasa_lunar_guide2024` | moon_field_proxy_17, moon_field_proxy_18, moon_field_proxy_19 | [공식 출처 / DOI](https://ntrs.nasa.gov/api/citations/20240011783/downloads/Lunar_Regolith_Simulant_Users_Guide_Rev_A_28OCT.pdf) | [nasa_lunar_simulant_guide_2024.pdf](01_ground/nasa_lunar_simulant_guide_2024.pdf) |
| `steel02` | 구성 법칙·접촉 구현 참고 | [공식 출처 / DOI](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Steel02.html) | [Steel02.html](04_model_docs/Steel02.html) |
| `concrete02` | 구성 법칙·접촉 구현 참고 | [공식 출처 / DOI](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Concrete02.html) | [Concrete02.html](04_model_docs/Concrete02.html) |
| `flat_slider` | 구성 법칙·접촉 구현 참고 | [공식 출처 / DOI](https://opensees.berkeley.edu/wiki/index.php/Flat_Slider_Bearing_Element) | [FlatSliderBearing.html](04_model_docs/FlatSliderBearing.html) |
| `dp_engine_340` | 구성 법칙·접촉 구현 참고 | [공식 출처 / DOI](https://raw.githubusercontent.com/OpenSees/OpenSees/v3.4.0/SRC/material/nD/UWmaterials/DruckerPrager.cpp) | [DruckerPrager_v3.4.0.cpp](04_model_docs/DruckerPrager_v3.4.0.cpp) |

## 함께 보관한 이전 조사 자료

- [이전 전체 참고문헌](05_project_notes/material/material_tuning_20261008/SOURCES.md): 재료·기초 형식의 배경 조사. 목록의 모든 논문에서 계수를 가져온 것은 아니다.
- [사용자 제공 기초 선행연구 요약](05_project_notes/mars_moon_dome_foundation_research_summary.md): 화성 전면기초·달 링기초의 개념 참고.
- UCSD 철근·6014 알루미늄·Smartfan CFRP·micro-concrete·달 시뮬런트 데이터베이스는 `03_candidates/`의 이전 조사·보정 자료다. 현재 SD400·6061·AS4 카탈로그와 혼용하지 않는다.
- 기존 보고서 사본은 당시 문맥을 보존했다. 그 내부의 옛 상대 경로는 이동으로 끊길 수 있으므로 이 목록과 `copy_manifest.json`으로 실제 사본을 찾는다.
- HTML 사본은 원문 응답 보관용이며 CSS·그림 등 외부 자산까지 저장한 완전한 오프라인 패키지는 아니다.
- NIAR AS4 원문은 웹 도구에서 조회했으나 파일 다운로드 응답이 PDF가 아니어서 저장하지 않았다. 원문을 확보한 것으로 표시하지 않는다.
- 최초 일반 실행의 DNS 제한 기록은 `download_log_attempt1.json`, 네트워크 다운로드 결과는 `download_log.json`에 남겼다.
