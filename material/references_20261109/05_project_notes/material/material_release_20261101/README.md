# 캡스톤 재료 적용 준비 패키지

**SI: m · N · s · kg · Pa.** 화면/보고서에서는 mm·kN·MPa로 변환한다. 원본 material/ground 코드는 수정하지 않았다. 재료 상수·근거는 `configs/catalog.json`, 검증 판정은 `results/verification_manifest.json`이 기준이다.

이 패키지는 구조물 코드가 오기 전 사용할 **문헌 기반 기준 재료와 검증된 계산 경로**를 제공한다. 모든 재료가 실험 보정되었다는 뜻은 아니다. 실제 돔/하이퍼볼로이드의 안정성·단면·기초 크기·동적 경계·충돌 파괴 검증은 포함하지 않는다.

## 실행

현재 확인한 실행환경은 Python3.10, OpenSees 엔진3.4.0이다. 설치된 `openseespy` 배포판 표시는3.7.1.2이므로 패키지 버전만 보고 엔진 버전을 판단하면 안 된다.

```sh
cd /Users/sydney/Desktop/capstone/Capstone_structure/material/material_release_20261101
/Users/sydney/anaconda3/envs/cap/bin/python run_verification.py
```

선정한 구현의 필수 검사가 모두 통과하면 종료코드0이다. 이전 RC Lobatto 방식의 실패 진단은 삭제하지 않고 별도로 기록한다. 다른 환경에서는 반드시 위 명령을 다시 실행한다. `build_catalog.py`는 원본 보정/관측자료에서 설정을 다시 생성하는 명령이며, 정상 사용에는 실행할 필요가 없다.

## 무엇을 사용할지

| 대상 | 기준 설정/구현 | 근거와 사용 범위 |
|---|---|---|
| 일반 강재 | `steel_s355_nominal`: E200GPa, fy355MPa, ρ7850 | S355 명목 기준. 실제 강종·두께별 보증값이 아님. 완전탄소성 비교 모델 |
| 보정 강재 | `steel_s355_coupon`: Steel02 | 기존 S355 HEM320 web 동일 출처 시편의 선택 이력에 한정. 현재 환경에서 독립 시편 재검증 |
| 6061-T6 | `al6061_t6_extrusion`: E70.359GPa@293K, fy240MPa, ρ2712.63 | Hydro 압출재 최소 항복강도 + NIST 탄성계수. 용접부·판재·열화·피로에 자동 적용 금지 |
| RC 프레임 | `register_rc_hinge_beam` | Concrete02 C30/Steel02 SD400 명목 가정, 비구속 콘크리트, 철근 면적 중복 제거. 힌지 길이와 실제 배근은 명시적 입력 |
| RC 기초 solid | `register_elastic_rc` | 철근 체적비 x/y/z 입력. 부착 완전·초기 탄성 Voigt 기준. 균열/부착슬립/다축 파괴 모델 아님 |
| CFRP 외피 | AS4/8552 UD RTD, `[0/90]s` 예시 | NIAR 시험 평균값. 허용응력 아님. OpenSees 연결은0/90 적층; 임의각 적층은 ABD 계산까지 제공 |
| 지구 지반 | `earth_sand_dr50_dry_reference` | 공식 PDMY02 제안 계수 + 건조 체적밀도1558.82kg/m³ 가정. 현장 보정 아님 |
| 달 지반 | `moon_jsc1a_18_peak` 및 `_residual` | JSC-1A 건조 삼축 요약자료. ρ1780, E20.45MPa, peak φ45.572°, residual φ38.1°, c0. 정적 E와 ψ0는 범위를 제한해서 사용 |
| 화성 지반 | `mars_jez1_peak` 및 `_residual` | JEZ-1 건조 직접전단: ρ1840, peak φ40.6°/c3.1kPa. E30MPa·ν0.3·ψ0는 가정; 강성0.5/1/2배 비교 |
| 차폐용 느슨한 레골리스 | `shielding` | 기본ρ1700kg/m³, 실제 면적·두께로 질량과 수직 자중 산정. 구조용 쉘 강성을 부여하지 않음 |
| 소결 구조용 레골리스 | `sintered_jsc1a_high_porosity` | 같은2016 시험군: E5.475GPa, ν.246, ρ2200, 압축강도 평균84.3MPa. 인장강도 미확인; 파괴 안전율 산정에 사용 불가 |
| SSI 접촉 | `register_contact` | 검증된 무점착 Coulomb 마찰/들뜸의 병진 접촉 연결. 면적별 Kn/Kt와 μ를 입력. 접착/인장 결합 미지원 |
| SSI 방사 감쇠 | `radiation_dashpots`, `register_dashpot` | ρVsA, ρVpA로 형상에 따라 환산. 속도-힘 관계 검증; 실제 경계 흡수 성능은 통합모델에서 확인 |

CFRP는 E1≈127.28GPa, E2≈9.24GPa, G12≈4.83GPa, ν12=.302, ρ1590kg/m³이다. 시험 평균 강도와 일부 가정된 두께방향 상수는 catalog에서 별도로 표시했다. 소결재와 결합재 레골리스는 다른 재료이므로 이번 선정은 **소결 JSC-1A**다.

## API와 연결 순서

- `capmat.materials`: `register_uniaxial`, `register_nd`, `register_plate`. 호출자의 도메인을 지우지 않는다. 태그는 호출자가 지정한다.
- `capmat.rc`: 실제 단면/철근 좌표를 받는 RC 단면과2D/3D HingeRadau 부재. 3D에는 GJ를 입력한다. 부재에 질량을 넣었으면 같은 자중 질량을 다시 절점에 추가하지 않는다.
- `capmat.laminate`: 물리적 ply의0/90 방향과 두께를 입력한다. 각 ply를 두께방향32개 계산층으로 나누어 굽힘 적분 오차를 줄였다. 요소의 local x축을 섬유0도 방향과 맞춰야 한다.
- `capmat.mechanics`: 임의 적층 ABD와 ply 응력 검토, 차폐 질량, 접촉 penalty, 감쇠계수 환산. `ply_screen`은 평균 강도 기반 초기 선별이며 설계 안전율이 아니다.
- `capmat.mass`: 별도 추가 질량을 기존 절점 질량에 누적하고 동일 source_key 중복을 거부한다. 재료가 자체 밀도를 제공하는지 반환 metadata를 확인한다.
- `capmat.ground.build_ground`: 선택한 원본 `ground/ground_modeling/ground_modeling.py`를 불러와3D bbarBrick 지반을 생성한다. **구조물을 만들기 전에 빈 도메인에서 호출**한다. 원본 builder는 wipe를 사용하므로 기존 절점/요소가 있으면 어댑터가 거부한다.
- `capmat.ssi`:3DOF 지반과3/6DOF 구조물의 병진 접촉. 보조6DOF 절점 태그·재료 태그·마찰 태그를 명시한다. 물리적 구조물 회전을 고정하지 않는다.

```python
import openseespy.opensees as ops
from capmat import register_uniaxial, register_plate
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)
steel = register_uniaxial(ops, 'steel_s355_coupon', 1001)
al_skin = register_plate(ops, 'al6061_t6_extrusion', 1101, 1102, 1201,
                         thickness_m=0.005)  # 예시 두께: 실제 구조물 값으로 교체
```

지반 예시:

```python
import numpy as np
from capmat.ground import build_ground
ops.wipe()
X, Y = np.meshgrid(np.linspace(0, 4, 3), np.linspace(0, 4, 3))
ground = build_ground(X, Y, np.zeros_like(X), body='Moon',
                      depth_m=3, n_depth=3)  # 검증용 소형 지반, 실제 기초 치수 아님
```

## 반드시 유지할 계산 조건

1. 금속 완전탄소성 nD 재료는 `DruckerPrager`의 마찰0·팽창0·경화0인 von Mises 극한으로 구현했다. 기존 엔진의 `J2Plasticity`는 Pa 단위에서 항복응력이100MPa 이상이면 내부 반복문이 건너뛰어지는 문제가 실제 시험에서 재현됐다. 현재 구현은 금속 쉘 반복 인장·압축400단계를 독립1D 해석해와 비교했다. 내부 DP의 큰 정수압 인장 cutoff는 일반 소변형 재료 비교 범위 밖이며, 충돌/고압 상태까지 확장하지 않는다.
2. DP 점착력 변환은 구현된 항복식에서 유도한 `sigmaY=6*c*cos(phi)/(3-sin(phi))`를 사용한다. 웹 문서에 인쇄된 추가√2 분모를 그대로 쓰면 점착력 기여가 작아진다. `research/dp_mapping_note.md`와 수정 전/후 검증을 참고한다.
3. RC는 임의 Lobatto 적분점 수를 바꿔 같은 소프트닝 응답이 나온다고 가정하지 않는다. 선정한 HingeRadau 모델에서 **힌지 길이를 물리적 입력으로 고정**한다. 검증 예제의0.30m를 모든 부재에 복사하지 않는다. 현재 wrapper는 `4*(lpI+lpJ)<L`인 양의 탄성 내부 구간을 요구한다.
4. 접촉은 작은 상대슬립·고정 법선·무점착 점 접촉 모델이다. 보조 절점의 bearing 회전 반력은 물리 절점으로 전달하지 않는다. 큰 미끄럼/회전/접촉면 이동에는 별도 표면 접촉이 필요하다.
5. **현재 엔진에서 SSI 힘 제어 조립은 `constraints('Plain')` + identity equalDOF + 동차 고정을 사용한다.** 보조절점 연결에 Transformation을 사용하면 물리 절점 변위 누적 오류가 재현됐다. 접촉 단독의 비영 SP 패치시험 설정을 힘 제어 모델에 복사하면 안 된다. 부등침하는 DisplacementControl 경로로 적용하거나 별도 검증된 다중지점 변위 입력을 구성한다.
6. 지반 builder의 고정저면·측면 롤러는 자중/정적 초기화용이다. 실제 지진 SSI에는 도메인 크기·Vs·주파수·입력운동에 맞는 동적 경계를 추가한다. PDMY 포화 밀도 원본 record를 건조 bbarBrick에 무심코 적용하지 않는다.
7. 차폐층은 수직 자중과 관성 질량을 모두 반영한다. 압력은 차압과 닫힌 압력 경계를 사용한다. 달의 중력을 적용한다고 관성 질량까지 줄이지 않는다.

## 실제 구조물 코드가 오면 이어서 할 일

1. 부재별 재료, 단면/두께, 배근, CFRP 적층·방향, 차폐 두께를 배정한다.
2. RC 힌지 길이/GJ, 기초 면적·깊이·접촉 분담면적, 접촉 재료쌍/μ를 결정한다. JEZ-1–매끈한 강재 μ≈.3424는 해당 실험 조합에만 근거가 있다. 나머지 조합은 미끄럼0/중간/상한 민감도와 접촉시험 근거를 구분한다.
3. 자중→내압→정적 평형→모드→기본 동적응답 순서로 질량·반력·변위·단위·경계 영향 확인 후 통합한다. 바닥/폐합부를 포함한 내압 전달을 확인한다.
4. 풍하중은 대기/장소 조건을 별도로 입력한다. 부등침하는 기초 경계 변위 시나리오다. 운석은 전역 충격응답과 국부 관통/파괴를 분리한다. 현재 재료 패키지로 관통·파편·고변형률 강도를 예측했다고 주장하지 않는다.

## 재현성/외부 의존

원본에서 읽는 자료는 S355 선정 시편 CSV 및 reviewed_config, JSC-1A triaxial_observations, 선택한 ground_modeling.py이다. 삭제/이동하지 않는다. `configs/catalog.json`의 sources와 `results/verification_manifest.json`의 코드 해시를 함께 보관한다. 권장 실행환경: numpy2.2.6, matplotlib3.10.9(그림만), plotly6.7.0(원본 지반 import), 현재 설치된 OpenSees 엔진3.4.0.

실험자료가 없는 6061-T6 이력경화, RC 실제 배근체/구속, CFRP 손상·박리, 소결 레골리스 인장/파괴, 달·화성 현장 저구속 지반은 **새 실험 보정이 끝난 것으로 표시하지 않는다.**
