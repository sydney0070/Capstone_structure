# 하이퍼볼로이드·돔 재료 선정과 지구·달·화성 지반 모델 조사

조사일: 2026-10-08. 목적: 재료를 아직 정하지 않은 캡스톤의 **후보 선정·OpenSeesPy 구현·시편 보정** 기반 마련.

추천하는 출발점은 **구조 재료 7종 비교 + 지구 모래/점토 + 달·화성 건조 레골리스 + 기초 접촉·마찰**이다. 지반의 레골리스와 이를 소결·결합해 만든 구조용 레골리스는 서로 다른 재료로 다룬다. 지하수·액상화와 말뚝은 지구 조건 또는 해당 기초 형식을 선택했을 때 추가한다.

현재 문헌의 실제 채택 사실, 공식 소프트웨어 지원, 이번 작업의 구현 제안을 구분했다. 특정 소재나 OpenSees 모델이 논문에 등장한다는 사실만으로 캡스톤 조건에서 검증되었다고 보지 않는다. 문헌 전체의 체계적 전수조사나 모든 원시 실험 데이터의 재현까지 수행한 것은 아니다.

## 1. 무엇을 재료로 준비해야 하는가

### 구조물 후보

| 후보 | 돔·하이퍼볼로이드에서의 역할 | 먼저 보정할 값 | 이번 구현 |
|---|---|---|---|
| 구조용 강재 | 격자 돔 부재, 링, 지지부, 금속 쉘 | E, ν, ρ, 항복강도, 경화, 반복응답 | `shell_steel`: J2 쉘; 격자 부재는 Steel02 경로 응용 |
| 알루미늄 합금 | 경량 골조·패널·기밀 라이너 후보 | 합금/열처리별 E, 항복강도, 경화, 용접부 저하 | `shell_aluminum`: J2 쉘 |
| 일반 콘크리트·철근콘크리트 | 지구 돔, 링기초, RC 지지부; 행성 환경에서는 비교 기준 | 압축·인장·연화, 철근 항복, 배근·구속 | 콘크리트/철근 단축 시험 + 탄성 콘크리트 쉘 |
| 소결 레골리스·소결 블록 | 현지 자원 기반 압축 부재, 외피·차폐체 | 공극률·제조 조건별 E, 압축/직접인장강도, 파괴에너지 | `shell_sintered_regolith`: 탄성 쉘 강성 보정 |
| 레골리스 지오폴리머 | 결합재를 이용한 외피·블록·돔 | 배합·양생별 E, 강도, 건조/진공/열 이력 | `shell_geopolymer`: 탄성 쉘 강성 보정 |
| 황–레골리스 콘크리트 | 화성 ISRU 구조재 비교 후보 | 황 비율, 입도, 온도별 강성·강도·취성 | `shell_sulfur_concrete`: 탄성 쉘 강성 보정 |
| CFRP 등 섬유복합재 | 경량 쉘·패널, 금속 라이너 보강 | E1/E2, G12, ν12, 적층 방향, 인장/압축·층간강도 | `shell_cfrp_elastic`: 직교이방성 탄성 쉘 |

섬유막·Kevlar·폴리에틸렌·단열재도 거주구 연구에 등장한다. 다만 하중을 지지하는 주구조, 기밀층, 단열·방사선 차폐층의 기능이 다르다. 이번의 7종은 주구조의 역학 비교를 위한 첫 묶음이다. 팽창식 구조를 선택하면 프리스트레스·인장전용·주름을 다루는 막 모델을 별도로 준비해야 한다. NASA의 복합재 거주구는 탄소복합재 압력용기와 폴리에틸렌 차폐층을 구분했고, HLIS 연구는 Aluminum 7075-T6와 Kevlar 29를 선택했다. [NASA 보고서](https://ntrs.nasa.gov/citations/20050185204), [HLIS 연구](https://www.sciencedirect.com/science/article/abs/pii/S0094576520306342).

### 지반 후보

| 환경·재료 | 준비할 상태 정보 | 적합한 첫 모델 | 주의할 조건 |
|---|---|---|---|
| 지구 모래·비소성 입상토 | 상대밀도/간극비, Vs, 입도, φ, 배수, 유효응력 | PDMY02; 액상화 중심이면 PM4Sand | 자갈로 확대할 때 입자파쇄·입도 효과 재검증 |
| 지구 점토 | su, PI, OCR, 압밀·배수, 강성감소 | PIMY: 빠른 비배수 전단의 간략화 | 장기 압밀·크리프·OCR 효과의 완전한 모델이 아님 |
| 지구 저소성 실트·점토 | su/σ′v, 반복연화, 수축, 응력이력 | PM4Silt | 모든 점토에 일반화하지 않음 |
| 지구 암반·매우 단단한 층 | 암반 E/ν/ρ 또는 Vs, 절리 상태 | ElasticIsotropic 기준 모델 | 암반 파괴·절리 미끄럼은 추가 모델 필요 |
| 달 레골리스 | mare/highland 구분, 시뮬런트 종류, 건조 밀도·입도·각짐, 낮은 압력의 c/φ/ψ | DruckerPrager 시작; PDMY02는 압력 의존/반복성 대안 | 수분 없이도 입자 맞물림·응집 효과가 가능하며 낮은 구속압 실험이 중요 |
| 화성 레골리스 | 풍성 모래/먼지/굳은 표층 구분, 밀도, c/φ/ψ, 빙결 여부 | 건조 층별 DruckerPrager 시작 | 광물 조성이 비슷한 시뮬런트가 역학적으로도 같은 것은 아님 |
| 결합된 표층·기반암·빙결층 | 층 두께, 결합 강도, 온도, 시간 의존성 | 초기에는 탄성/마찰성 층으로 구분해 민감도 조사 | 얼음·시멘트화 거동을 일반 모래로 확정하지 않음 |

달과 화성을 각각 하나의 고정 물성 세트로 표현하지 않는다. 특히 **화성의 ‘점토 광물 존재’와 지구 포화 점토의 PIMY/PM4Silt 거동은 동의어가 아니다.** 물리적 시험 조건에 맞춰 분류한다.

### SSI 후보

| 현상 | 실제 준비할 시험/상수 | 구현 후보 |
|---|---|---|
| 링기초·받침 아래 압축과 들뜸 | 압축 하중–침하, 무인장 여부, 초기 틈 | 압축 스프링/접촉 요소; 필요 시 지반 연속체 |
| 기초 바닥 미끄럼 | 접촉 재료 조합·거칠기별 직접전단, δ 또는 μ | zeroLengthContact2D/3D; ContactMaterial + 적합한 접촉 요소 |
| 분포 지반 반력과 rocking | 수직·수평·회전 재하 시험, 접지면 이동 | 얕은기초 BNWF 또는 연속체+접촉 |
| 말뚝 횡저항/주면마찰/선단 | 깊이·직경·지반별 p–y/t–z/q–z | PySimple1/TzSimple1/QzSimple1 + zeroLength |
| 동적 방사감쇠 | 밀도·파속·경계 면적 또는 동적 임피던스 | Viscous + zeroLength 또는 흡수 경계 |

SSI의 마찰계수는 지반 내부마찰각 φ와 같다고 가정하지 않는다. 바닥면의 접촉각 δ를 측정하고 μ=tanδ로 정한다. 패널 사이 조인트나 철근 부착은 구조 내부 접합 문제로, 지반 SSI와 따로 관리한다.

## 2. 기존 연구에서 실제로 사용한 것

아래 표의 ‘OpenSees 적용’은 문헌의 채택 사실과 이번 제안을 구분한다. 상세 서지는 [SOURCES.md](SOURCES.md)에 있다.

| 연구·확인 범위 | 실제 소재/모델·지반 처리 | 캡스톤에 가져올 내용 |
|---|---|---|
| **Mrozek et al. (2025), 달 분화구 덮개 돔** — 공개 본문 | 레골리스 지오폴리머. ABAQUS S4R 쉘, 선형탄성. 지반–기초 상호작용 상세 모델 없음; 링 경계 고정/힌지 비교 | 돔·ISRU 소재의 직접 사례. OpenSees의 탄성 쉘로 옮기는 것은 **이번 구현 제안**. 이 논문을 달 지반의 비선형 SSI 검증 사례로 인용하면 안 됨 |
| **Ruess, Schaenzlin & Benaroya (2006)** — 대학 메타데이터·공개 본문 발췌 | 지구 제작 고강도 알루미늄 부재의 달 거주구; 아치·타이 설계 | 알루미늄을 후보군에 넣을 근거. 코드의 예시 항복강도는 이 논문의 물성 재현값이 아님 |
| **NASA/CR-2005-213845** — NASA 초록 | 탄소복합재 거주구·차폐층의 다층 설계 | CFRP를 등방 금속과 구분하고 적층·차폐 기능을 분리 |
| **Wan, Wendner & Cusatis (2016)** — 저자 공개 초록 | 화성 시뮬런트+황. 압축·쪼갬·3점휨 시험, LDPM 수치모델 | 황콘크리트의 실험 근거. **Concrete02나 OpenSees를 사용한 연구가 아님** |
| **Gupta et al. (2023)** — 저자 공개 초록 | LHS/LMS/MGS 시뮬런트 소결, 수치 격자 모델과 시험 | 느슨한 레골리스와 소결 구조재의 구분. 시험별 최대 압축강도를 모든 소결체 공통값으로 복사하지 않음 |
| **Elgamal et al. (2003)** — 출판사 공개 교정 부분 | 포화 비점착토 cyclic mobility, 다중항복면 소성; Nevada No.120 모래 약 Dr=40% 시험자료 | PDMY 계열의 실험·이론 계보. 오늘의 PDMY02 명령이 논문에서 그대로 사용됐다고 소급하지 않음 |
| **Boulanger & Ziotopoulou (2017/2018)** — 개발자/공식 문서 | PM4Sand/PM4Silt, 지진 시 모래·저소성 실트/점토의 반복응답 | 지구 포화토의 대안. 개발자 FLAC 최신 버전과 OpenSees 구현 버전은 별개 |
| **Boulanger et al. (1999)** — 출판사 초록 | 연약점토 위 조밀한 모래의 원심모형, 동적 p–y 모델 검증 | 말뚝 SSI의 실험 근거. 이후 Py/Tz/Qz 구현 계보의 근거이며 1999년에 Python 코드를 썼다는 뜻이 아님 |
| **Harden et al. (2005); Gajan et al. (2008)** — PEER 공개 보고서/초록 | 모래·점토 위 얕은기초, 분포 비선형 Winkler와 들뜸·침하·회전 | 돔 링기초 SSI에서 우선 참고. 말뚝용 스프링을 무조건 재사용하는 접근과 구분 |
| **McGann & Arduino, OpenSees 공식 말뚝 예제** — 코드/설명 | dispBeamColumn 말뚝, PySimple1/TzSimple1/QzSimple1, zeroLength, equalDOF | 명령 조합을 직접 참고할 수 있는 구현 예제. 연구 논문의 새로운 검증 결과와는 구분 |
| **Kramer, Arduino & Shin (2008)** — PEER 초록 | 액상화/비액상화 지반·말뚝·교대·교량을 OpenSees에 통합 | 지구 전체 SSI의 체계 참고. 초록만으로 특정 요소·상수까지 확정하지 않음 |
| **Luu, Mo & Hsu (2017)** — 출판사 공개 설명 | CSMM 기반 적층 RC 쉘, OpenSees를 기반으로 SCS-3D 구현; 패널·벽·탱크 검증 | RC 쉘은 1D Concrete02+Steel02만으로 완성되지 않음. 사용자 구현 연구를 기본 설치 명령과 혼동하지 않음 |
| **Karakas, Ozgan & Daloglu (2018)** — 출판사 초록 | 하이퍼볼릭 냉각탑, 3변수 Vlasov 지반; SAP2000–MATLAB 연계 | 링기초/지반 유연성 비교 사례. **OpenSees 연구가 아님** |
| **OpenSees 공식 RC·쉘 예제** — 명령 문서/예제 | Concrete02/Steel02 + fiber 보·기둥; 탄성 쉘 + ShellMITC4 | 구현 문법과 요소 연결 경로의 직접 근거 |

위 첫 행의 직접 출처: [달 돔 논문](https://www.nature.com/articles/s41598-025-07901-x). 지오폴리머 E=8.13 GPa는 논문에서 강도 기반 경험식으로 산정한 값이다. 5.5 MPa는 본문 시험 설명상 **휨인장강도**이므로 직접인장 `ft`로 바로 넣지 않았다. 밀도는 본문에서 1750/2050 및 단위 표기에 불일치가 있어 정밀 재현 전 확인이 필요하다. 이번 `shell_geopolymer`의 ρ=1.75 tonne/m³는 초기 비교 가정이고, ν=0.2·시편 두께 0.02 m도 이번 시험용이다.

검색한 범위에서 **‘동일한 돔/하이퍼볼로이드를 지구·달·화성에 두고, 각 지반을 보정된 OpenSeesPy 재료로 비교한 공개 코드’는 확인하지 못했다.** 따라서 행성 구조 논문에서 소재·하중을, 지구 OpenSees SSI 연구에서 구성모델·요소 조합을 각각 가져와 검증 단계를 거치는 것이 현실적이다. 검색되지 않았다는 것이 관련 연구가 전혀 없다는 뜻은 아니다.

## 3. 행성 레골리스 물성의 근거를 어떻게 선택할 것인가

### 달

NASA의 레골리스 기계적 성질 자료는 실제 달 토질에서 밀도·입도·입자 형태와 낮은 구속압이 중요하고, 지구 시험에서는 막압·시험기 마찰·시편 자중의 영향이 커질 수 있음을 설명한다. [NASA 기계적 성질 자료](https://ntrs.nasa.gov/citations/19940004207).

JSC-1A 연구는 상태 의존성의 좋은 예다. NASA의 ring shear 연구에서는 미리 전단한 시료가 대략 무점착·φ≈40°인 반면, 진동 다짐 후 시험에서는 약 c≈1 kPa·φ≈60°의 다른 결과가 보고됐다. 이들은 **같은 시험 조건의 보편적인 상수가 아니다**. 이를 보고 ‘달의 마찰각=60°’로 고정하면 안 된다. [NASA/TM-2009-215625](https://ntrs.nasa.gov/citations/20090023168).

국내 재료를 고려하면 KLS-1도 후보이다. 개발 연구는 Apollo/JSC-1/FJS-1과 화학·입도·비중·전단 특성을 비교했다. 국내 실험실에서 구할 수 있는 시료의 특성을 확인하고, 선택한 밀도·압력에서 직접전단 또는 삼축시험을 수행하는 방향으로 연결할 수 있다. [KLS-1 개발 연구](https://ascelibrary.org/doi/10.1061/%28ASCE%29AS.1943-5525.0000798).

### 화성

시뮬런트에는 JSC Mars-1, MMS 계열, MGS-1 등이 있지만 목적이 광물학·자원 처리인지, 지반의 강도·주행성인지 확인해야 한다. NASA의 geotechnical review는 이런 대표성 문제를 다룬다. [NASA 화성 시뮬런트 지반 검토](https://ntrs.nasa.gov/api/citations/20200003046/downloads/20200003046.pdf).

다짐한 달·화성 시뮬런트에 대한 2024년 연구도 입도·밀도가 점착력과 전단강도에 영향을 준다고 보고한다. 이를 바탕으로 ‘행성명별 상수’보다 **시뮬런트–밀도–구속압–제조/다짐 상태별 자료**를 준비해야 한다. [직접전단 연구](https://www.sciencedirect.com/science/article/pii/S0019103524000010).

### 이번 레골리스 코드의 의미

`moon_regolith_dp`, `mars_regolith_dp`는 c=0.5 kPa, φ=35°, ψ=5°를 둔 **가상 시작점**이다. 특정 달/화성 착륙지, NASA 시뮬런트 또는 위 논문의 교정값을 의미하지 않는다. 같은 기본 골격을 사용해 입력을 쉽게 바꿀 수 있게 했으며, 두 사례의 응답 차이를 실제 행성의 재료 차이라고 해석할 수 없다.

Drucker–Prager 선택은 낮은 점착력과 마찰·팽창을 3D에서 조절하기 쉽다는 구현상 제안이다. Mohr–Coulomb과의 변환은 **삼축압축 자오선 일치**이며, 모든 응력경로에서 동일한 파괴면이 아니다. `rho`라는 이름도 DP의 마찰 계수와 질량밀도에 중복 사용될 수 있어 코드에서는 `friction`과 `rho`를 구분했다. DP의 기본 모델만으로 입자파쇄, 진공 응집, 온도, 순환 압밀을 충분히 설명한다고 가정하지 않는다. [DP 공식 문서](https://opensees.berkeley.edu/wiki/index.php/Drucker_Prager).

## 4. OpenSeesPy material–section–element 연결

| 목표 | material / section | element 및 자유도 | 제공 코드 |
|---|---|---|---|
| 3D 탄성 지반 | ElasticIsotropic nD | SSPbrick, `ndm=3, ndf=3` | `elastic_soil` |
| 3D 모래 | PressureDependMultiYield02 nD | SSPbrick; 전체 모델 대안 bbarBrick | `sand_pdmy02` |
| 3D 간략 점토 | PressureIndependMultiYield nD | SSPbrick | `clay_pimy` |
| 달·화성 건조 레골리스 | DruckerPrager nD | SSPbrick | 행성 DP 2종 |
| PM4 재료 보정 경로 | PM4Sand/PM4Silt nD | 제공 쿠폰: SSPquad, `ndm=2, ndf=2`; 수압 연성 전체 모델 후보: SSPquadUP | PM4 2종, 수압 DOF 없는 일정체적 시편 |
| 콘크리트·철근 단축 | Concrete02, Steel02 uniaxial | truss 시편; fiber section → force/dispBeamColumn | 3종 단축 시편 |
| 금속 쉘 | J2Plasticity nD → PlateFiber wrapper → PlateFiber section | ShellMITC4, `ndm=3, ndf=6` | 강재·알루미늄 2종 |
| 취성재 초기 쉘 비교 | ElasticIsotropic → PlateFiber → PlateFiber section | ShellMITC4 | 콘크리트·소결·지오폴리머·황콘크리트 4종 |
| 직교이방성 쉘 | ElasticOrthotropic → PlateFiber → PlateFiber section | ShellMITC4 | CFRP 탄성 1종 |
| 말뚝 스프링 | PySimple1, TzSimple1, QzSimple1 | zeroLength, 방향별 재료 | 3종 스프링 |
| 무인장 압축 틈 | ElasticPPGap, Fy/gap 음수 | zeroLength | `normal_gap` |
| 압력에 종속되는 미끄럼 | 요소 자체에 Kn/Kt/μ | zeroLengthContactASDimplex, implicit 적분 | `interface_friction` |
| 방사감쇠 | Viscous | zeroLength 동적 연결; 쿠폰은 속도 직접 입력 | `radiation_dashpot` |

공식 경로: [SSPbrick](https://openseespydoc.readthedocs.io/en/latest/src/SSPbrick.html), [PlateFiber section](https://openseespydoc.readthedocs.io/en/latest/src/plateFiberSection.html), [J2Plasticity](https://openseespydoc.readthedocs.io/en/latest/src/J2Plasticity.html), [접촉 요소](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/zeroLengthContactASDimplex.html). 접촉 시편은 `-intType 0`을 사용한다. 이 선택은 공식 기능과 이번 해석해·증분 검증에 근거한 구현 제안이며, 달·화성 돔 논문에서 이 요소를 사용했다는 뜻은 아니다.

**1D 재료와 쉘의 중간 단계가 중요하다.** Concrete02를 정의했다고 콘크리트 돔의 다축 균열이 생기는 것은 아니다. RC 쉘을 최종 선택하면 적합한 평면응력 콘크리트와 `PlateRebar`, `LayeredShell` 등의 조합을 정하고 판/쉘 검증 시험을 추가한다. 공식 RC 벽 예제도 이러한 별도 연결을 보여준다. [공식 RC 적층 쉘 예제](https://openseespydoc.readthedocs.io/en/latest/src/RCshearwall.html). 이번에는 단축 콘크리트 보정과 탄성 쉘 비교를 제공하며, 비선형 RC 쉘을 완성했다고 표시하지 않았다.

ShellMITC4 쿠폰이 실행된 것은 재료–단면–요소의 연결 및 막응답을 확인한 것이다. 돔 좌굴, snap-through, 초기 형상 오차를 검증한 것이 아니다. 형상 비교 시 큰 회전/기하비선형을 지원하는 쉘 선택과 벤치마크를 별도로 확인한다.

PM4Sand·PM4Silt의 여기서 사용한 OpenSees 인터페이스는 2D 평면변형률이다. 개발자 웹사이트에 PM4Sand3D가 존재한다고 해서 이 명령을 그대로 SSPbrick에 넣을 수 있는 것은 아니다. [OpenSees PM4Sand](https://openseespydoc.readthedocs.io/en/latest/src/PM4Sand.html), [개발자 사이트](https://pm4sand.engr.ucdavis.edu/).

제공 PM4 코드는 초기압밀 후 높이를 고정한 `SSPquad` 시편이다. `ru_equivalent=1−σ′v/σ′v0`는 유효 수직응력 감소를 나타내는 대리지표이며 실제 간극수압이 아니다. 수압/배수/압밀/CRR를 연구하려면 `SSPquadUP` 등의 별도 동적 u–p 시험을 구축·검증해야 한다. 검증 중 u–p 경계조건의 증분 세분화 불안정을 확인했으므로 그 시험을 검증된 코드로 배포하지 않았다.

## 5. 재료별 보정 순서와 필요한 실험

| 모델 | 1차 보정 | 2차 보정·검증 | 피해야 할 보정 |
|---|---|---|---|
| ElasticIsotropic 지반 | 밀도 독립 측정, Vs와 G=ρVs²; ν 별도 | 다른 압력·방향의 초기 기울기 | 하나의 전단곡선으로 E·ν·ρ를 동시에 유일하게 결정 |
| PDMY02 모래 | 같은 초기 유효응력의 G, K, 강도면 | 여러 압력의 전단/체적 변화, 수축·팽창, 반복 강성·감쇠 | 단일 곡선으로 모든 수축·팽창·액상화 인자를 최적화 |
| PIMY 점토 | 초기 G 및 전단강도 | 진폭별 이력, 강성감소 | 장기 압밀 자료를 빠른 전단 모델로 맞춤 |
| 레골리스 DP | 다압력 강도 포락선의 c·φ, 초기 G/K | 체적변형으로 ψ, 변형경화 H, 매우 낮은 압력의 독립 시편 | 지구 100 kPa 결과만으로 수 kPa의 행성 얕은지반 판정 |
| PM4Sand/PM4Silt | 상대밀도 또는 su, 작은 변형 강성 | 독립 시험의 반복연화·응력이력; 별도 응력제어 CRR | 이번 일정체적 쿠폰을 간극수압/CRR 실험과 동일시 |
| Concrete02 | fpc, epsc0, 압축 후연화; 인장 자료가 있으면 ft/Ets | 재하–제하 기울기 λ, 반복 이력 | E=2fpc/epsc0 관계 무시, 휨인장값을 그대로 직접인장 ft로 입력 |
| Steel02 | E, Fy, b | R0/cR1/cR2, 필요 시 등방경화 인자 | 한 방향 인장시험만으로 반복 이력·좌굴·저주기피로까지 주장 |
| 금속 J2 쉘 | E·ν, sig0, H/포화경화 | 두 축 하중·반복시험, 실제 합금 열처리 | 일반 강재 계수를 알루미늄 합금의 검증값으로 사용 |
| 취성 쉘 4종 | 단축·초음파 등 독립 측정한 E·ν·ρ | 직접인장/압축, 균열·파괴에너지, 환경 시험 | 탄성 쿠폰으로 파괴강도나 기밀성을 보정했다고 주장 |
| CFRP 탄성 | 0°/90°/전단 시험의 E1/E2/G12/ν | 적층 방향·굽힘·환경·접합 | 0° 축시험 하나로 모든 이방성 상수를 결정 |
| Py/Tz/Qz | 용량과 50% 동원 변위 | 반복·틈·인발, 여러 깊이·말뚝 시험 | 재료 `soilType`를 각 명령에서 같은 의미로 간주 |
| 접촉·마찰 | 하중별 μ, 접선 초기기울기 | 다른 정상력·들뜸·재접촉 | penalty Kn을 실제 지반 강성과 무조건 동일시 |
| Viscous | 동적 속도–힘 및 파속 기반 C | 주파수·경계 변화 | 정적 응력–변형률 곡선으로 C를 식별 |

PDMY02 권장표에는 강성 기준압력 80 kPa와 `refPress` 101 kPa가 함께 표기된다. 따라서 예시의 수치를 단일 정밀 교정 세트처럼 복사하지 않았다. 코드의 PDMY02 값은 기존 사용자 시험과 공식 예시 계열의 시작점이며, G/K가 정의된 `p_ref`와 실제 시험의 평균 유효응력을 구분하고 압력 의존식에 맞춰 기준 강성을 다시 보정해야 한다. 시험마다 `p_ref`를 실제 압력으로 덮어쓰라는 뜻은 아니다. [공식 PDMY02](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html).

3D 전단 시편의 초기 압축은 측방 변형을 구속한다. 초기 수평응력은 입력 G/K 또는 E/ν에서 정해지므로, 현장 K0를 독립적으로 지정한 시험은 아니다. `sigma_v`와 평균응력은 같지 않으며 후자는 `history.csv`의 `p_mean`으로 확인한다. 현장과 다른 초기응력비가 중요하면 초기 재하 경계를 확장해야 한다.

Concrete02 구속 코어와 피복은 별도 재료 태그다. ‘코어 압축강도를 임의로 올린 것’을 구속 검증으로 보지 않고, 배근·횡구속 조건에 맞는 데이터/구속모델을 사용한다. 이번의 코어 설정은 동작 예시다. [Concrete02](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Concrete02.html), [Steel02](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Steel02.html).

## 6. 세 행성을 공정하게 비교하는 방법

먼저 동일한 재료·형상으로 중력만 바꾸는 시험과, 각 행성에 맞는 재료·환경을 함께 바꾸는 시험을 나눈다. 전자는 중력 효과를 분리하고 후자는 시나리오를 비교한다. 둘을 섞으면 재료와 중력의 기여를 알 수 없다.

1. 재료 시편 보정: 압력·밀도·온도·배수 조건별로 독립 수행.
2. 형상 고정 비교: 동일 형상·두께·하중에서 재료 교체. 질량·변형·응력 비교.
3. 동등 질량 또는 동등 성능 비교: 두께·단면을 조절. ‘같은 두께’ 결과와 분리.
4. 고정/힌지 지지 기준 모델 → 선형 지반 스프링 → 접촉/마찰 → 필요하면 3D 지반 연속체.
5. 사용 하중: 자중, 내부–외부 **압력차**, 레골리스 차폐층 무게, 장비·탑재, 온도, 필요 시 지진/월진/충격을 별도 정의.
6. 최종 판단: 질량·응력뿐 아니라 취성, 좌굴, 균열·누설, 제조·운송, 접합을 함께 평가.

단위계를 kN–m–s–tonne–kPa로 통일했다. ρ=1.5 tonne/m³, 깊이 1 m라는 **동일한 가상 건조 지반**의 자중응력은 지구 14.715 kPa, 달 2.43 kPa, 화성 5.565 kPa이다. 이는 물성 추천값이 아닌 σv=ρgz 계산 예다. 행성 중력은 근삿값 9.81/1.62/3.71 m/s²를 사용했다. [NASA 행성 비교표](https://nssdc.gsfc.nasa.gov/planetary/factsheet/), [NASA 화성 중력](https://solarsystem.nasa.gov/mars-by-the-numbers/).

**질량을 중력 비율로 줄이지 않는다.** 달에서도 같은 물체의 질량은 같고 무게가 작아진다. 지반의 강성·강도가 압력 의존적이므로 중력이 바뀐 후 초기 평형을 다시 형성한다. PDMY/DP의 `pa` 또는 `refPress`는 구성모델의 기준압력이며 달의 진공이나 화성의 대기압으로 자동 치환하지 않는다.

내압의 힘은 중력에 비례해 감소하지 않는다. 따라서 행성 비교에서 자중만 바꾸고 동일한 내부 압력차를 유지하면 들뜸·앵커 인장·링의 응력이 중요해질 수 있다. 무인장 접촉에 필요한 인장 지지력을 임의 suction으로 넣는 대신 실제 앵커·타이의 존재와 성능을 따로 모델링한다.

지반·쉘의 자유도도 다르다. 3D 지반 절점은 보통 병진 3개, 쉘·보 절점은 병진+회전 6개이다. 접촉용 3자유도 절점을 별도로 만들고 구조 절점과 필요한 병진만 연결하는 등 일관된 조립이 필요하다. 회전까지 지반과 `equalDOF`로 묶으면 실제보다 고정된 기초가 될 수 있다.

## 7. SSI 구현에서 특히 많이 틀리는 점

- `PySimple1.pult`: 흔한 설계식의 kN/m 값에 **담당 말뚝 길이**를 곱한 kN 입력.
- `TzSimple1.tult`: 주면 전단강도 kPa × 담당 주면적 m² → kN.
- `QzSimple1.qult`: 선단 지지응력 kPa × 선단 면적 m² → kN.
- Py의 type 1/2는 Matlock 점토/API 모래 관계이며, Tz와 Qz의 type 1/2는 각각 다른 출처 관계를 뜻한다. 하나의 동일한 soilType 사전으로 일괄 처리하지 않는다.
- 얕은기초 `ShallowFoundationGen` 문서에는 QzSimple2/PxSimple1/TxSimple1 등 수정 모델이 등장한다. `QzSimple1` 하나를 붙였다고 문헌의 링기초 BNWF를 재현한 것은 아니다.
- 접촉의 Kn/Kt는 penalty 및 접촉 변형을 정하는 값이다. Kn을 늘려 침투를 줄이는 수치 수렴 연구와 실제 지반 탄성/침하 보정은 구분한다.
- 접촉 면이 열렸을 때 미끄럼 저항도 달라져야 한다. 독립적인 수평 스프링만으로 이를 자동 구현했다고 가정하지 않는다.
- 스프링의 dashpot, 지반 흡수 경계, Rayleigh 감쇠를 함께 쓸 때 같은 감쇠 기여를 중복 계산하지 않는다.

근거: [Py](https://openseespydoc.readthedocs.io/en/latest/src/PySimple1.html), [Tz](https://openseespydoc.readthedocs.io/en/latest/src/TzSimple1.html), [Qz](https://openseespydoc.readthedocs.io/en/latest/src/QzSimple1.html), [얕은기초 생성기](https://opensees.berkeley.edu/wiki/index.php/ShallowFoundationGen).

## 8. 제공물과 남아 있는 실험 입력

`cases/`와 `configs/`에는 23개 개별 시작점이 있다. `materials.py`는 재료·단면 생성 함수, `experiments.py`는 시편과 재하 경로, `tune.py`는 실행·sweep·실험 CSV 최적화, `planet_conditions.py`는 중력/응력·단위 보조 함수다. [실행 방법](README.md), [검증 결과](VALIDATION.md).

23개 기본 시편은 실행 완료했고 해석해·평형 등 구현 검사를 수행했다. 그러나 **아직 어떤 소재도 해당 돔과 행성 조건의 최종 재료로 채택하거나 실제 시험에 보정 완료한 것은 아니다.** 이를 결정하려면 후보별 시편 곡선과 선정한 형상의 부재·쉘·기초 검증이 필요하다.

먼저 확보할 자료는 (1) 지반의 시료명·밀도·입도·낮은 압력의 전단/체적변형, (2) 구조재의 압축/인장·E·ρ·온도 및 제조 이력, (3) 기초 바닥 조합의 정상력별 마찰·침하, (4) 구조 형상·두께·내압·차폐층·기초 방식이다. 이 자료가 들어오면 현재 JSON/CSV 인터페이스로 후보별 보정을 진행할 수 있다.
