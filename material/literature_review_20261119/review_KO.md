# 화성 지반 층별 구성과 OpenSees 단축 재료 후보

조사일: 2026-11-19. 범위는 현재 프로젝트가 따르는 **InSight 착륙지의 천부 지반**이다. 화성 전체에 같은 층서가 있다는 가정은 하지 않는다. 문헌 조사와 명령 인식 확인만 수행했으며, 새 강도값 선정·층별 재료 보정·전체 지반 재해석은 수행하지 않았다.

## 핵심 판단

- 현재 층서는 현무암질 물질이 입도, 결합, 파쇄 상태에 따라 달라지는 지반으로 해석하는 것이 적절하다. 깊은 층의 광물별 함량을 직접 측정한 자료는 이번 조사에서 확보하지 못했다.
- **3D 지반층의 구성 법칙과 기초의 단축 스프링은 서로 다른 모델이다.** 현재 bbarBrick에는 nDMaterial이 필요하다. uniaxialMaterial 후보는 기초 접촉·지반 반력을 단축 힘–변위 관계로 줄여 표현하는 용도다. [D1](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/bbarBrick.html)
- 단축 비교 모델의 1순위는 수직 방향 **QzSimple2**, 단순 압축 접촉 비교는 **ElasticPPGap**, 실험 곡선을 직접 맞추려면 **Hysteretic**이다. 수평 거동까지 포함하면 **PySimple2·TzSimple2**를 함께 검토한다. 이 순위는 현재 링 기초의 목적을 고려한 판단이며 화성 보정 결과가 아니다.

## 1. 다른 층은 무엇으로 이루어져 있는가

다음 깊이는 **현재 코드의 층 구분**이다. 관측으로 확정된 균일한 경계 깊이라는 뜻은 아니다.

| 현재 모델 구간 | 구성과 상태 | 관측·해석의 구분 | 모델에 필요한 거동 |
|---|---|---|---|
| 0–0.2 m: 경화층, duricrust | 모래·잔입자 속에 자갈과 작은 암편이 결합된 층 | 사진에서 결합된 입자와 급경사 벽이 관찰됨. 염류에 의한 결합은 제안된 설명이며 결합제의 종류·함량이 확정된 것은 아님 | 결합에 따른 강도, 압축 항복, 인장 분리·결합 손상 가능성 |
| 0.2–3 m: 모래질 레골리스 | 충돌로 부서진 비교적 세립의 모래와 일부 자갈·암편 | 착륙지 관측과 충돌구 해석에 근거한 근지표 해석 | 마찰, 구속압 의존 강도, 다짐과 반복 전단 |
| 3–13 m: 조립질 분출물, blocky ejecta | 비교적 큰 현무암질 암편과 입상 기질이 섞인 조립질·각력질 물질 | 더 깊은 층의 구성은 충돌구·지형 해석으로 추정. 시추 코어와 층별 직접전단 자료를 확보한 상태는 아님 | 입자 간 마찰, 구속압 효과, 재배열·다짐; 결합·입자 파쇄는 추가 자료 필요 |
| 13–30 m: 파쇄 현무암 | 균열·절리로 나뉜 현무암질 용암류를 등가 암반으로 표현 | 깊이와 파쇄 정도는 모델화된 해석. 입상 모래나 무균열 현무암 코어와 구분 | 암석 자체의 강도와 절리 미끄럼·벌어짐, 암반의 방향성 |

Golombek 외(2020) 그림 9는 얇은 모래 아래의 경화층과, 약 3 m 세립 레골리스 아래 조립질 분출물·파쇄 현무암이라는 해석 단면을 제시한다. 근지표 약 10 cm의 직접 관찰과 깊은 부분의 추정을 구분한다. 해당 그림의 경화층은 약 5–10 cm이며, 현재 코드의 0.2 m는 유효층 설정이다. [G1](https://elib.dlr.de/139605/1/Golombek_et_al.InSight_Landing_Site.Nat_Comm_2020.pdf) Warner 외(2020)의 충돌구 기반 층서도 참고했다. 13 m 경계는 현재 참조 모델의 경계로 취급한다. [G2](https://ntrs.nasa.gov/api/citations/20205003548/downloads/Warner_impact_origin_InSight_2020_2019JE006333_optimized.pdf)

광물 수준에서는 **사장석(plagioclase), 휘석(pyroxene), 감람석(olivine)**과 유리질·비정질 성분, 철 산화물 등을 포함하는 현무암질 물질이 비교 기준이다. Cannon 외(2019)의 MGS-1은 Gale 분화구 Rocknest 조성을 바탕으로 설계됐다. 이 조성이나 모사토의 배합률을 InSight 깊은 층의 실측 함량으로 옮기지 않는다. 조립질층과 암반층을 구분하는 데는 광물 이름 외에도 입도·공극·결합·절리 상태가 중요하다. [G3](https://sciences.ucf.edu/class/wp-content/uploads/sites/23/2019/01/Cannon_MGS-1.pdf)

Morgan 외(2018)는 착륙 전 추정에서 낮은 점착력의 모래, 결합된 토양과 밀도·마찰각의 관계를 정리한다. 그 자료도 현재 3–13 m층이나 13–30 m층의 항복값을 직접 보정해 주지는 않는다. [G4](https://elib.dlr.de/122371/1/Revised_InSight%20Landing%20Site%20Regolith%20Physical%20Properties%20Complete%202018-6-6.pdf)

## 2. OpenSees 연구에서 실제로 사용한 재료

공개 검색에서 원문과 구현을 확인할 수 있었던 범위에서는 **InSight의 위 층들을 각각 보정한 OpenSees uniaxialMaterial 사례를 확인하지 못했다.** 아래는 지구의 모래·점토 얕은기초 연구와 공식 구현 자료다. 화성 전용 검증 자료로 분류하지 않는다.

| 문헌 | 실제 모델과 검증 | 확인한 재료·역할 | 이번 프로젝트에 가져올 수 있는 부분 |
|---|---|---|---|
| Harden 외, PEER 2005/04 | OpenSees BNWF, 지구 모래·점토 실험 비교 | §3.3에서 ElasticPPGap, 일반 Hysteretic, QzSimple1, PySimple1, TzSimple1 검토 | 들림·누적침하·수평저항을 나누는 구성. 기본 EPP-gap의 누적침하 제약도 검토함 |
| Raychowdhury & Hutchinson, 2009 | 사각·띠기초 원심모형시험과 비교 | 비선형 스프링·감쇠기·갭의 분포 모델 | 얕은기초의 침하·회전·미끄럼 검증 방법. 링 기초로의 적용은 추가 검토 |
| NIST GCR 12-917-21, 2012 | 연구와 실무를 연결하는 기술 지침 | q–z: QzSimple2, p–x: PxSimple1, t–x: TxSimple1로 설명 | 수직 지지와 들림, 측면 수동저항, 밑면 미끄럼을 구분하는 체계 |
| OpenSees ShallowFoundationGen 문서 | 공식 얕은기초 생성 예제 | 예제 명령은 QzSimple2, PySimple2, TzSimple2 | 설치본의 명령 이름과 입력 형식을 확인하는 자료 |

근거: [R1](https://peer.berkeley.edu/sites/default/files/peer_504_c._harden_t._hutchinson_g._martin_b._kutter.pdf) §3.3, 인쇄 27–37쪽; [R2](https://onlinelibrary.wiley.com/doi/10.1002/eqe.902) 공개 초록; [R3](https://nehrp-consultants.org/publications/download/nistgcr12-917-21REV1.pdf) 인쇄 2-36–2-38쪽; [R4](https://opensees.berkeley.edu/OpenSees/manuals/usermanual/documents/ShallowDocumentation.pdf) 예제 명령. R2의 정확한 명령 이름은 공개 초록만으로 판정하지 않고 R3·R4와 교차 확인했다. R3의 미국 정부 서버 PDF는 웹 파서 용량 제한이 있어 작성기관의 수정본을 읽었다.

## 3. uniaxialMaterial 후보

**다음은 재료 종류의 후보이며 아직 보정된 수치 세트가 아니다.** 토층별 화학 성분과 재료 명령을 일대일로 대응시키지 않는다.

| 후보 | 적합한 역할 | 선정에 필요한 데이터 | 주의할 점 | 추천 위치 |
|---|---|---|---|---|
| **QzSimple2** | 기초 밑 수직 지지, 항복, 들림 및 반복 변형 | 분담 극한수직력, 기준 변위, 초기강성·배면곡선, 인장저항과 감쇠 | 얕은기초 모델의 적용 근거가 있으나 지구 실험 기반. 링 형상·화성 중력에 맞춘 보정 필요 | **수직 스프링 1순위** |
| **ElasticPPGap** | 압축 접촉, 분리, 단순 탄소성 비교 | 스프링 강성, 항복력, 초기 갭, 경화비, damage/recentering 규칙 | 기본 복원 규칙과 damage 옵션을 구분. 선택한 규칙의 잔류변형을 확인해야 함 | **가장 단순한 접촉 비교** |
| **Hysteretic** | 측정한 압축·전단 반복 곡선을 직접 맞춤 | 양·음 방향 힘–변위 점, 핀칭·손상·제하강성 | 구속압이나 절리 상태가 자동 반영되지 않음. 실험 곡선이 있어야 장점이 큼 | **데이터를 확보한 뒤 보정용** |
| **QzSimple1** | 수직 q–z 비교 기준 | qzType, 극한력, Z50, 인장저항비, 감쇠 | 원래 말뚝 선단의 모래·점토 배면곡선. 얕은 링 기초에는 그대로 적용하지 않음 | 기존 연구 재현·비교용 |
| **PySimple2** | 매입 기초 옆면의 수동저항과 갭 | 극한수평력, 기준 변위, drag와 감쇠 | 수직 지지용 재료가 아님. NIST의 PxSimple1 표기와 설치 명령을 구분 | 수평 거동 추가 시 |
| **TzSimple2** | 기초 밑면 전단·미끄럼 | 극한전단력, 기준 변위, 접촉면 마찰·감쇠 | 수직 접촉력이 줄거나 들릴 때 전단저항도 줄도록 연계할 필요 | 수평 거동 추가 시 |

ElasticPPGap의 압축 부호와 복원·손상 옵션은 공식 문서를 따른다. [D5](https://opensees.berkeley.edu/OpenSees/manuals/usermanual/170.htm) Hysteretic은 양·음 배면곡선, 핀칭과 손상 및 제하강성 규칙을 지정하는 일반 모델이다. [D6](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Hysteretic.html) QzSimple1의 극한값은 응력이 아닌 분담 **힘**이고 Z50은 변위다. [D2](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/QzSimple1.html) PySimple1과 TzSimple1도 원래 말뚝용 관계와 힘 단위 변환을 확인해야 한다. [D3](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/PySimple1.html) [D4](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/TzSimple1.html)

### 현재 설치본에서의 이름 확인

`cap` 환경의 OpenSees 엔진은 **3.4.0**이다. 임의의 구문 확인용 값으로 생성 명령의 인식 여부만 확인했다. 하중 해석·매개변수 보정·거동 검증을 한 결과가 아니다.

- 인식: QzSimple1, QzSimple2, PySimple1, PySimple2, TzSimple1, TzSimple2, ElasticPPGap, Hysteretic.
- 미인식: 문헌에 표기된 이름 그대로의 PxSimple1, TxSimple1.
- 따라서 실제 코드 후보 이름은 R4에서 사용한 QzSimple2/PySimple2/TzSimple2를 기준으로 검토한다. 문헌 이름만 보고 함수명이나 입력값을 그대로 옮기지 않는다.

[명령 인식 확인 기록](uniaxial_command_availability.json)

## 4. 현재 3D 지반과 연결하는 방법

현재 bbarBrick의 재료 태그는 **nDMaterial**을 가리킨다. `uniaxialMaterial` 태그를 바로 넣어서 층을 바꾸는 방법은 사용할 수 없다. [D1](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/bbarBrick.html)

| 목적 | 우선 검토 모델 | 적용 범위 |
|---|---|---|
| 경화층·조립질층의 3D 탄소성 | DruckerPrager | 압력 의존 항복의 단순 기준. c·φ·팽창·경화 보정 필요 |
| 조립질층의 반복 전단·다짐 관련 거동 | PressureDependMultiYield02 | 3D 가능. 모래 기반 반복 거동 후보이며 저압·건조 화성 조건과 각력질 입자에 대한 보정 필요 |
| 파쇄 현무암의 3D 응답 | 탄성 기준 + 보정된 암반 구성모델/절리 인터페이스 | 모래 재료를 그대로 적용하지 않음. 절리·암석 강도와 모델화할 스케일을 먼저 결정 |
| 층상 지반 전체를 줄인 기초 반력 모델 | 분포 QzSimple2 + 필요 시 수평 스프링 | 층별 기여를 반영한 등가 힘–변위 관계를 보정 |

PressureDependMultiYield02는 3D 차원과 압력 의존 강성, 반복 수축·팽창 관련 계수를 제공하지만 공식 예시의 지구 모래값이 화성의 보정값은 아니다. [D7](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html) DP의 항복·경화·소성 체적변화 설정도 별도의 입력이 필요하다. [D8](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html)

### 단축 후보를 수치 세트로 만들기 전 필요한 것

1. **0–0.2 m 경화층:** 제하를 포함한 압축 곡선, 결합 파괴·인장 분리 특성. 현재 층 두께는 유효화 가정으로 관리.
2. **3–13 m 조립질층:** 입도·공극률·밀도, 해당 구속압에서의 삼축/직접전단 강도, 압축·제하 및 반복 곡선. 입자 재배열과 파쇄를 포함할지 결정.
3. **13–30 m 파쇄 현무암:** 암석 코어의 압축·인장 강도, 암반의 절리 간격·방향·상태와 절리 전단/수직 강성. Hoek–Brown 등을 쓸 경우 GSI·mi·D와 적용 응력 범위도 근거 필요.
4. **기초 스프링 보정:** 각 스프링의 분담 면적, 하중–침하·제하 관계, 극한력과 기준 변위. 현재 완료된 탄성 하부층 해석은 초기강성 비교 자료가 될 수 있지만 극한지지력이나 항복변위를 결정해 주지는 않는다.

스프링의 강성 단위는 보통 N/m이며 기존 E의 단위 Pa와 다르다. 단순한 수평변형 구속 1D 기둥이라면 M=K+4G/3, k=M A/h로 초기 압축강성을 근사할 수 있다. 이는 링 기초의 3D 영향 범위와 층간 전단을 자동으로 표현하는 식이 아니다. 또한 스프링 자체에 3D 지반의 밀도·질량을 옮기는 것은 아니므로 질량 배분도 따로 정한다.

3D 지반의 변형성과 등가 지반 스프링을 함께 넣으면 같은 유연성을 중복 계산할 수 있다. 실제 적용 시에는 기초 접촉만 스프링으로 표현하는지, 층상 지반 전체를 스프링으로 축약하는지를 먼저 정해야 한다. 수평 스프링의 용량과 수직 접촉력의 연계도 검토한다.

## 5. 프로젝트에 권하는 다음 선택

**3D 층별 지반 모델은 nDMaterial로 유지하고, 단축 비교 모델은 QzSimple2부터 검토**하는 접근을 추천한다. 접촉의 단순 비교용으로 ElasticPPGap을 추가하고, 층별 시험 곡선을 확보하면 Hysteretic을 보정한다. 매입 링 기초의 수평 응답까지 필요해질 때 PySimple2·TzSimple2를 더한다. 현재는 새 c·φ·항복력·변위를 지정하지 않는다.

[선별 후보와 요구 입력](uniaxial_candidates.json) · [문헌별 확인 범위](sources.json)

## 6. 문헌 및 확인 범위

- **G1** Golombek et al., 2020. [Geology of the InSight landing site on Mars](https://elib.dlr.de/139605/1/Golombek_et_al.InSight_Landing_Site.Nat_Comm_2020.pdf). Figure 9; discussion on printed pp. 6–7
- **G2** Warner et al., 2020. [An Impact Crater Origin for the InSight Landing Site at Homestead Hollow, Mars](https://ntrs.nasa.gov/api/citations/20205003548/downloads/Warner_impact_origin_InSight_2020_2019JE006333_optimized.pdf). Near-surface stratigraphy interpretation
- **G3** Cannon et al., 2019. [Mars global simulant MGS-1: A Rocknest-based open standard for basaltic martian regolith simulants](https://sciences.ucf.edu/class/wp-content/uploads/sites/23/2019/01/Cannon_MGS-1.pdf). Sections 1–2; mineralogical design
- **G4** Morgan et al., 2018. [A Pre-Landing Assessment of Regolith Properties at the InSight Landing Site](https://elib.dlr.de/122371/1/Revised_InSight%20Landing%20Site%20Regolith%20Physical%20Properties%20Complete%202018-6-6.pdf). Sections 1.3 and 2.3–2.4
- **R1** Harden, Hutchinson, Martin and Kutter, 2005. [Numerical Modeling of the Nonlinear Cyclic Response of Shallow Foundations](https://peer.berkeley.edu/sites/default/files/peer_504_c._harden_t._hutchinson_g._martin_b._kutter.pdf). PEER 2005/04, Chapter 3, Sections 3.3.1–3.3.5; printed pp. 27–37, PDF pp. 46–56
- **R2** Raychowdhury and Hutchinson, 2009. [Performance evaluation of a nonlinear Winkler-based shallow foundation model using centrifuge test results](https://onlinelibrary.wiley.com/doi/10.1002/eqe.902). Earthquake Engineering & Structural Dynamics 38, 679–698; accessible abstract
- **R3** NEHRP Consultants Joint Venture / NIST, 2012. [Soil-Structure Interaction for Building Structures](https://nehrp-consultants.org/publications/download/nistgcr12-917-21REV1.pdf). GCR 12-917-21, printed pp. 2-36–2-38, PDF pp. 68–70
- **R4** Raychowdhury / OpenSees, 2008. [ShallowFoundationGen implementation documentation](https://opensees.berkeley.edu/OpenSees/manuals/usermanual/documents/ShallowDocumentation.pdf). Example material commands for shallow foundations
- **D1** , . [OpenSees bbarBrick documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/bbarBrick.html). 
- **D2** , . [OpenSees QzSimple1 documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/QzSimple1.html). 
- **D3** , . [OpenSees PySimple1 documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/PySimple1.html). 
- **D4** , . [OpenSees TzSimple1 documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/TzSimple1.html). 
- **D5** , . [OpenSees ElasticPPGap documentation](https://opensees.berkeley.edu/OpenSees/manuals/usermanual/170.htm). 
- **D6** , . [OpenSees Hysteretic documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Hysteretic.html). 
- **D7** , . [OpenSees PressureDependMultiYield02 documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html). 
- **D8** , . [OpenSees DruckerPrager documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html). 

## 조사 기록과 제한

검색어에는 `InSight subsurface layers duricrust blocky ejecta fractured basalt`, `Mars OpenSees regolith`, `Martian OpenSees soil material`, `OpenSees shallow foundation QzSimple2`, `Harden Hutchinson 2005 BNWF` 등을 사용했다. 검색 결과의 2차 요약만으로 특정 화성 보정 사례를 주장하지 않았다. R2는 공개 초록과 서지·참고문헌을 확인했고, 재료별 상세 구현은 읽을 수 있는 R1·R3·R4로 확인했다. 이것은 모든 비공개 연구를 포함하는 체계적 문헌고찰은 아니다.

기존 MGS-1/MGS-1S 저장본과 지반 재료 카탈로그는 이 조사에서 변경하지 않았다. 전날 작성 중이던 `ground_layers_20261118/materials.py`는 입력 JSON과 보정이 없는 미완성 초안이며, 문헌 조사 완료 전 실사용 재료로 채택하지 않았다.

## 원문 중앙 보관 위치

외부 원문은 [기존 논문 모음](../../docs/references/index.html)의 지반·SSI·공식 문서 분류에 보관했다. 원문 미확보 자료는 해당 목록에 별도로 표시했다. 조사 ID와 보관 ID는 [연결 기록](../../docs/references/_catalog/research_20261119_crosswalk.json)에서 확인한다.

## 층별 적용 추가 정리

[층·부위별 추천과 근거](layer_recommendations_KO.md)에 3D 층 재료와 단축 반력·접촉 후보를 구분해 정리했다.

## 첨부 자료 추가 검토

[첨부 요약을 반영한 원문 대조·층별 후보](attached_summary_assessment_KO.md)에서 천부 수십 cm와 심부 층서를 분리하고 3D 재료 후보를 추가 정리했다.
