# 3차원 혼합 지반·SSI 기준 모델

실행 진입점은 **`ground_modeling.py`**, 실제 새 구현은 **`ground_hybrid/`**이다. 지구·화성·달을 같은 API로 전환한다. 단위는 **m · N · s · kg · Pa**다.

기존 파일의 DEM·기하 함수는 유지하고 기본 실행 경로를 새 코드로 연결했다. 이전 전체 파일은 [이전 코드](../archive/ground_modeling_20261116/backups/ground_modeling_before_hybrid_20261102.py)에 보관했다. `--legacy`는 이전 실행 방식이다. 새 builder는 `ops.wipe()`를 호출하지 않는다. CLI만 자신의 도메인을 초기화한다.

## 바로 실행

검증한 환경: `/Users/sydney/anaconda3/envs/cap/bin/python`, OpenSees **엔진 3.4.0**. 배포 패키지 버전과 엔진 버전은 다를 수 있다. 설치 변경은 하지 않았다.

```sh
cd /Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Earth --demo-frame
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Moon --demo-frame
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Mars --demo-frame
```

기본 예제는 **6×6m, 최소 깊이 12m의 작은 연결 시험 모델**에 2×2m, 두께 0.4m 탄성 RC 판기초를 둔다. 지반·기초 자중 후 **추가 10kPa**를 재하한다. `--demo-frame`은 기초 모서리에 3m 강재 기둥 4개와 총 4,000kg 집중질량을 추가한다. **돔/하이퍼볼로이드의 대체 설계나 확정 기초 치수가 아니다.** 기본 메시의 영역 크기 독립성은 확보되지 않았다.

12×12m의 더 큰 시험 설정은 `ground_hybrid/configs/*_large_candidate.json`에 보존했다. 지구·화성의 큰 기초/기둥 예제는 실행됐지만, 달의 큰 비선형 모델은 수렴 문제가 남았다. 작고 검증 가능한 기본 예제로 바꾼 것은 **연결·실행 범위를 명확히 하기 위한 선택**이며, 큰 모델의 문제가 해결됐다는 뜻은 아니다. 수렴을 위해 달의 기존 점착력·마찰각을 임의로 변경하지 않았다.

큰 모델이나 변경한 하중을 시험할 때는 별도 프로세스의 실행 시간 제한을 권장한다. 제한을 넘으면 해당 자식 계산만 종료하고 종료 코드 124를 반환한다.

```sh
/Users/sydney/anaconda3/envs/cap/bin/python -m ground_hybrid.bounded_run --timeout-s 180 --config ground_hybrid/configs/Moon_large_candidate.json --demo-frame --output hybrid_results/my_moon_trial
```

초기에는 `--no-html`로 계산만 실행할 수 있다. 세부 설정은 `ground_hybrid/configs/Earth.json`, `Moon.json`, `Mars.json`을 복사·수정하여 `--config 파일경로`로 전달한다. 출력 기본 위치는 `hybrid_results/천체/`이고 재실행 시 같은 결과를 갱신하므로 별도 보존에는 `--output 다른폴더`를 사용한다.

## HTML에서 지형·재료·변형 보기

기본 실행은 계산이 성공하면 **`ground_viewer.html`을 만들고 기본 브라우저로 연다.** 외부 서버나 인터넷 연결은 필요하지 않다. 브라우저를 자동으로 열지 않으려면 `--no-open`, HTML 자체를 생략하려면 `--no-html`을 사용한다. 기존 `ground_cutaway.html`도 동일한 뷰어 사본으로 유지한다.

- 전체 3D 지반, Y 중간면의 반쪽 단면, 지표면만 보기 전환. 단면은 요소 단위이며 임의 평면 재메시는 아니다.
- 재료별 표시/숨기기, 체적비·입력 E·밀도 확인, 요소 경계 표시.
- 재료 색상, 초기 표고, 변위 크기, 수직 변위 색상 전환.
- 미변형/전체 변위/조립체 자중 평형 이후/동적 시작 이후의 변위 기준 선택.
- 화면 변형 확대 배율, 실제 크기 ×1, 자동 확대. 실제 변위는 **mm**, 표시 좌표는 **m**다.
- 동적 해석 시 재생·일시정지·시간 슬라이더·속도 조절, 저장 프레임별 최대 변위 그래프.
- 기초 판과 `--demo-frame` 기둥도 절점 변위로 표시한다. 회전 곡률·균열이나 입자 이동을 그리는 모델은 아니다.

기본 `--terrain flat`은 평평한 기준 지형이다. 제공된 실제 DEM을 보려면 다음처럼 실행한다.

```sh
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Moon --terrain dem --prepare-pad --demo-frame
```

동적 형상 재생 예제:

```sh
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Mars --terrain dem --prepare-pad --demo-frame --dynamic pulse --dt 0.001 --steps 300 --output hybrid_results/Mars_dynamic_viewer
```

`pulse`는 총100N의 합성 수평 펄스이며 실제 지진 기록이 아니다. 기본 실행에는 동하중이 없으므로 재생 버튼도 비활성화된다. HTML의 조작은 저장 결과의 표시만 바꾸며 새로운 하중으로 재해석하지 않는다.

애니메이션은 OpenSees에서 실제로 수렴한 절점 병진 변위를 사용한다. `transient.json`의 `visualization`에 절점 태그, 실제 시간, 3축 변위와 동적 시작 기준을 저장한다. 기본은 초기·마지막 상태를 포함해 최대151프레임이며 `--html-frames 301`처럼 조절할 수 있다. 저장 간격을 바꿔도 적분 시간 간격이나 기존 관측점의 매 스텝 기록은 바뀌지 않는다. 저장 프레임 사이의 순간 최대 응답은 화면 그래프에서 놓칠 수 있다. `--no-html`은 화면용 변위 이력을 생략한다.

표시 좌표는 `초기 좌표 + 배율 × 선택한 기준 대비 변위`다. 따라서 동적 시작 기준으로 보면 초기 형상 위에서 진동 성분을 분리해 보여준다. 전체 침하와 진동을 함께 보려면 전체 변위 기준을 선택한다. 매우 작은 평형 드리프트를 큰 진동처럼 보이지 않도록 자동 확대는 변위1e-10m 미만에서 ×1을 유지한다. 색상 범위와 축 범위는 같은 기준의 저장 시간 전체에서 고정한다.

구조물 API를 직접 사용하는 경우, 동적 적분 전에 `ground_hybrid.viewer.visualization_node_tags(ground, frame)`로 화면 절점을 구하고 `run_transient(..., visualization_nodes=tags)`에 넘긴다. 반환된 `visualization`을 `export_results(..., visualization=result['visualization'])`에 전달한다. 기본 도우미가 그리는 구조물은 기초와 예제 기둥이며, 다른 구조물의 요소 형상을 표시하려면 해당 연결 정보를 별도로 확장해야 한다.

## 채택한 단순화와 혼합 방식

- 3D `bbarBrick` 연속체. 하나의 균질 등가재료로 평균내지 않고 **층과 렌즈 영역마다 별도 재료 태그**를 배정한다.
- 층 경계는 지표면 아래 깊이를 기준으로 한다. 마지막 바닥은 평면이다. 내부 재료 경계는 절점을 공유하는 완전 부착이다. 암석 접촉면의 분리·입자 굴림·미세입도 혼합은 해석하지 않는다.
- 천체별 물성은 연구용 기준값이며 특정 부지의 시추 결과가 아니다. 층 깊이·렌즈 위치·기반암 깊이는 설정 파일에 명시한 **가정**이다. DEM은 지표 형상만 제공한다.

| 천체 | 기본 층 구성 | 비균질 영역 |
|---|---|---|
| 지구 | 0–2m 느슨한 모래, 2–6m 중간 밀도 모래, 6–10m 조밀한 모래, 10m 아래 암석 | 중간층 내 느슨한 모래 렌즈 |
| 달 | 0–0.3m JSC-1A 1.7군, 0.3–1m 1.8군, 1–8m 1.9군, 8m 아래 암석 | 깊은 레골리스 내 1.7군 렌즈 |
| 화성 | 0–0.2m 표층 경화층, 0.2–3m 모래 우세 레골리스, 3–13m 조립질 레골리스, 13m 아래 파쇄 현무암 | — |

지구와 달의 렌즈 중심 `(2.5,0,4.5)`m의 세 번째 값은 **지표면 아래 깊이**이며, 반경은 `(1.5,2,1.5)`m다. 작은 기본 영역에서는 렌즈 일부가 경계에서 잘린다. 실제 체적분율은 요소별 계산 후 출력하며 메시마다 달라질 수 있다. 화성의 3 m와 13 m는 InSight 착륙지 문헌에 근거한 기준 경계이며, 기본 12 m 영역은 조립질 레골리스 안에서 끝난다. [화성 층서 및 물성 근거](../../docs/MODEL_OVERVIEW_20261109_KO.md#2-실제-채택한-지반-재료)

## 재료 근거와 추가 보정의 의미

`ground_hybrid/configs/materials.json`에 원본 재료 카탈로그 SHA256, 출처, 상수, 가정을 기록했다. 기존 `material_release_20261101`을 덮어쓰지 않았다.

- **지구:** OpenSees 공식 Dr30/50/75 제안값. G/K의 80→101kPa 기준압력을 일관되게 환산하고, 가정한 입자밀도 2650kg/m³와 공극비로 건조 체적밀도를 계산했다. 기본은 이 E/G/K·마찰각을 사용하는 **Drucker–Prager 정적 단순화**다. 일정 강성과 완전소성 마찰을 사용하며 PDMY02의 압력 의존 강성·수축·팽창·반복 이력을 재현하지 않는다. [공식 계수 표](https://opensees.berkeley.edu/wiki/index.php/PressureDependMultiYield02_Material)
- **달:** 이전 JSC-1A 삼축 요약자료 보정의 밀도별 세 재료를 재사용했다. φ≈42.97/45.57/49.48°, E≈15.10/20.45/17.25MPa다. 조밀할수록 E가 반드시 커지도록 값을 조작하지 않았다. 시험 구속압은 15–60kPa이며 얕은 달 지반으로의 외삽은 불확실하다. 밀도 증가의 일반 방향은 [NASA/Connolly·Carrier 자료](https://ntrs.nasa.gov/api/citations/20220018853/downloads/LPSC%202023%20Abstract_Connolly_Carrier_v2.pdf)와 비교하되, JSC 실험밀도를 실제 달 밀도분포로 오인하지 않는다.
- **화성 레골리스:** 모래 우세층과 조립질 층은 등가 탄성체로 모델링한다. E는 각각 28 MPa와 520 MPa다. [Kenda 외, 2020 표 1](https://doi.org/10.1029/2020JE006387)
- **추가 표층 재료:** InSight 연구의 c5.8kPa·φ30°·ρ1200kg/m³ 조건, 천부 평균 파동 기반 E11.48MPa·ν0.28을 대용재로 채택했다. 파동 평균 물성을 표층의 정적 강성에 옮긴 것은 별도 가정이다. 첫 1cm 먼지층은 정지된 표면 기준에서 생략했다. [Spohn 외, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9734249/)
- **암석:** 지구와 달은 E40GPa의 현무암 대용재를, 화성은 13 m 아래에 E13GPa의 파쇄 현무암 대용재를 사용한다. 모두 등가 탄성체다. [Musiol 외, 2016](https://gfzpublic.gfz.de/rest/items/item_1802896_4/component/file_1802898/content) · [Kenda 외, 2020 표 1](https://doi.org/10.1029/2020JE006387)

새 실험곡선을 얻어 모든 재료를 재최적화한 것이 아니다. **기존 보정 재사용 + 문헌 상수 확장 + 구성식 환산의 수치 검증**을 수행했다. DP는 삼축압축 Mohr–Coulomb 강도에 맞춘 변환이며 모든 응력 경로의 직접전단 결과를 그대로 재현하지 않는다. 강성 0.5/1/2배 민감도 결과를 함께 제공한다.

`moon_field_proxy_17/18/19`는 별도 연구 후보다. [NASA Lunar Regolith Simulant User's Guide, Table 18](https://ntrs.nasa.gov/api/citations/20240011783/downloads/Lunar_Regolith_Simulant_Users_Guide_Rev_A_28OCT.pdf)의 현장 점착력 0.52/0.90/3.0kPa를 JSC-1A 강성·밀도·마찰각과 조합했다. **다른 자료의 결합이며 JSC 실험 재보정이 아니다.** 원래 문헌의 깊이 구간과 현재 층 구간도 다르다. 큰 모델의 자중 단계는 통과했지만 추가 하중은 완료되지 않아 기본값으로 채택하지 않았다. 단위요소 강도 검증의 통과가 조립 모델의 검증을 대신하지 않는다.

지구의 PDMY02 선택형:

```sh
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Earth --earth-model pdmy02 --no-html
```

PDMY02는 단위요소 반복전단 검증을 포함한다. 복합 모델의 초기화는 기본적으로 0.1% 탄성 자리잡기 하중 후 stage1로 전환하고, 실제 중력을 100%까지 비선형 재하한다. **전체 중력을 탄성으로 넣은 뒤 단계만 바꾸는 방식과 경로가 다르다.** 강도나 허용오차를 낮춰 수렴을 맞추지 않는다. 어떤 모델이든 실패하면 `FAILED.json`을 저장하고 종료한다. 일반 기본 검증의 완전 통과를 선택형의 모든 형상·하중 보증으로 확대하지 않는다.

## 구조물 코드 연결

```python
import openseespy.opensees as ops
from ground_hybrid import default_config, build_ground, add_raft, initialize_gravity
from ground_hybrid.interface import pressure_loads, raft_response
from ground_hybrid.solver import apply_static_loads

ops.wipe()                       # 전체 모델 소유자만 초기화
config = default_config('Moon')
ground = build_ground(config)    # builder 내부에는 wipe가 없음
raft = add_raft(ground, bounds=(-1, 1, -1, 1), thickness_m=0.4)

# 여기에서 실제 구조물을 생성한다.
# raft['nodes']: 6자유도 구조물 연결 절점.
# ops.nodeCoord(tag)로 확인하고 beam/shell의 기초 절점으로 직접 재사용한다.
# 새 노드를 쓸 때는 일치 좌표를 확인하는 connect_matching_nodes를 사용한다.
# 구조물 질량과 자중은 모델 작성자가 각각 한 번만 넣는다.
structure_selfweights = {}       # {node: [Fx,Fy,Fz,Mx,My,Mz]}
initialize_gravity(ground, extra_loads=structure_selfweights)
apply_static_loads(ground, pressure_loads(raft, 10_000.))
print(raft_response(ground, raft))
```

- 지반은 3자유도, 기초 판과 구조물 연결 절점은 6자유도다. 기초 회전을 임의로 고정하지 않는다.
- SI 질량밀도는 요소에 포함한다. 중력 하중은 별도 적분한다. 같은 질량을 절점에 다시 추가하지 않는다.
- 하중 제어 접촉은 **`Plain` 제약 처리기 + identity `equalDOF`**를 사용한다. 이 경로에 `rigidLink`, `rigidDiaphragm`, 연쇄 MP 또는 비영 SP를 섞으면 별도 검증이 필요하다.
- 동일 좌표가 아닌 절점을 가까운 절점에 자동 연결하지 않는다. 기초 외곽은 지반 메시 선과 일치해야 한다.
- 초기 지반·기초·구조물을 먼저 조립한 뒤 중력을 넣는다. 이미 자중 변형한 지반에 무응력 기초를 나중에 설치하는 시공 단계는 구현하지 않았으며 API가 거부한다. 추가 침하는 **조립체 자중 평형 이후**를 기준으로 계산한다.
- 기본 판은 탄성 RC, 무균열, 접촉면과 판 기준면의 편심 0 가정이다. 실제 기초 두께·배근·편심·균열강성은 구조물 설계 시 교체한다.
- 기초 미끄러짐/들뜸은 이전 검증된 점접촉을 사용한다. 기본 μ는 지구·달 0.35 가정, 화성은 JEZ-1/매끈한 강판의 tan18.9°를 전용한 가정이다. 콘크리트 기초 접촉의 실측 계수로 표시하지 않는다. 점착은 0이다.
- `displacement_path`는 **실험용** 단일 절점 병진 제어 함수다. 조립 모델의 큰 미끄러짐·반전 경로는 실패했으며 설계용 미끄러짐 용량으로 사용하면 안 된다. 여러 지점의 강제 부등침하나 다중 지진 입력도 구현하지 않았다. 기존 점접촉 단품 검증과 새 조립체 전체의 검증 범위는 다르다.
- 예약 태그 기본: 지반 10000부터, 판/접촉 300000부터, 예제 기둥 450000부터, 경계 600000부터, 하중 800001 등. 실제 구조물과 중복되지 않게 예약한다. 반복 하중 패턴에는 서로 다른 태그가 필요하다.

## DEM / 실제 지형

```sh
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Moon --terrain dem --prepare-pad
```

HGT/PDS 파일을 NumPy memmap으로 직접 읽는다. 이 실행환경에는 rasterio가 없지만 새 경로에는 필요하지 않다. 지원 대상은 제공된 HGT, 달 SignedLSB2 PDS4, 화성 MSB_INTEGER16 PDS3 형식이다. 다른 지도 투영/형식은 조용히 추정하지 않고 거부한다.

X=동쪽, Y=북쪽이며 국소 접평면 근사를 쓴다. 달의 위도 약 −82.5°에서 동서 간격에는 cos(latitude)를 적용한다. 기본 지점은 지구 기존 포항 좌표, 달/화성 제공 타일 중앙이다. `site_latitude_deg`, `site_longitude_deg`로 타일 내 위치를 지정할 수 있다.

기본 위치의 원래 간격은 지구 약25×31m, 달 약3.9×29.6m, 화성 약429×463m다. **0.5m 해석 메시가 0.5m 실측 지형이라는 뜻은 아니다.** 원본 샘플 해시, 형식, 좌표, 해상도를 결과에 기록한다. 결측값과 타일 밖 보간은 거부한다.

`--prepare-pad`는 초기 지형을 기초 아래 평탄하게 하고 주변에 천이 구간을 만든다. 실제 굴착·성토의 응력 이력 계산이 아니다. 이 옵션 없이 경사진 지형에 평면 접촉기초를 연결하면 오류로 멈춘다. 지반만 볼 때는 `--no-raft`를 사용한다.

## 동적 단계

```sh
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Mars --dynamic rest --steps 25 --no-html
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --body Mars --dynamic pulse --steps 100 --no-html
```

자중·정적 하중 평형의 지지반력을 외력으로 보존한 뒤 경계 구속을 해제하고, 면적 적분된 ρVsA/ρVpA 점성 경계를 적용한다. 모서리에는 인접한 각 면의 기여를 합한다. 초기 평형 유지와 균질 지반의 수직 입사 S파를 검증했다. 자유장 경계나 임의 각도·층상 지진 입력의 완성 검증은 아니다.

`rest`는 추가 동하중이 없는 평형 검사, `pulse`는 총100N 합성 수평 펄스, `incident`는 합성 상향 입사 속도파다. 실제 지진 기록을 의미하지 않는다. `add_incident_base_velocity` 입력은 **입사 속도 m/s**이며 기저 힘은 2ρcAv이다. 노두 가속도를 그대로 전달하거나 해제된 바닥에 UniformExcitation을 함께 넣지 않는다.

Newmark 평균가속도법과 실패 단계 이분을 사용한다. Rayleigh 감쇠는 기본 0이며, 물리적 자료 없이 세 천체에 같은 감쇠율을 넣지 않는다. 필요 시 `run_transient(rayleigh=(감쇠비,f1,f2))`로 명시한다. 동적 경계로 전환한 모델은 강체 모드 때문에 정적 해석을 다시 실행하지 않는다. 모델을 재생성해야 한다.

## 결과·검증

- `summary.json`: 물성별 체적비, 질량, 자중 평형, 기초 추가/부등침하, 접촉력.
- `structure_ports.json`: 6자유도 연결 절점, 좌표, 분담면적, 자중 기준 변위.
- `gravity_stress.json`: 자중 직후 응력. K0를 임의로 덮어쓰지 않고 계산된 수평/수직 응력비를 확인한다.
- `ground.vtu`: ParaView용 실제 3D 지반 요소·재료 ID·변위. 판/기둥 요소는 이 지반 VTU에 포함하지 않는다.
- `ground_viewer.html`: 외부 접속 없이 볼 수 있는 지형·단면·재료 분포·변형 결과 뷰어. `ground_cutaway.html`은 호환 사본.
- `manifest.json`: 실행 엔진, 입력/코드 해시.
- `RUN_STATUS.json`: CLI 실행 완료 여부. `success`·`complete=true`인 실행만 완료 결과로 취급한다. 실패 시 `FAILED.json`이 생긴다. 강제 시간 종료 시 `running`·`complete=false`가 남을 수 있으므로 이전 출력 파일을 새 결과로 오인하지 않는다. 성공한 재실행은 이전 실패 표시를 제거한다.

```sh
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --verify
/Users/sydney/anaconda3/envs/cap/bin/python -m ground_hybrid.smoke
/Users/sydney/anaconda3/envs/cap/bin/python -m ground_hybrid.viewer_verification
/Users/sydney/anaconda3/envs/cap/bin/python -m ground_hybrid.sensitivity --timeout-s 60
```

기존 검증의 [요약](../archive/ground_modeling_20261116/hybrid_verification/summary.json)과 [민감도 결과](../archive/ground_modeling_20261116/hybrid_verification/sensitivity.json)를 보관했다. **해석 성공과 메시·영역 크기 독립성은 별도 판정**이다. 지하수/액상화, 온도·진공 변화, 절리·입자파괴, 실제 돔의 하중조합/안전성은 아직 포함하지 않는다. 실제 구조물 치수를 받으면 기초·영역·메시·접촉계수·입력파에 대해 다시 확인해야 한다.

민감도 시험은 각 조건을 별도 프로세스로 실행한다. 현재 21조건 중 20조건이 해석을 완료했고 달의 넓은 영역은 실패/시간 제한에 걸렸다. 지구 영역 확대와 화성 수평 메시 세분화에서 침하 변화가 임시 10% 기준을 초과했다. 영역 확대는 외곽 메시 및 영역 안에 포함되는 렌즈 체적도 바꾸므로 순수한 경계거리 효과만 분리한 시험은 아니다. **현재 결과는 메시 독립 설계값이 아니다.** 자세한 수치와 다음 작업은 `GROUND_HYBRID_REPORT_KO.md`를 참조한다.
