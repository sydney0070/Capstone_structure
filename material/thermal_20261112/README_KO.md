# 19개 재료 열팽창 입력과 재료 단위 재검증

2026-11-12. 사용자가 구조물 코드의 기존 온도 설정을 기준으로 사용하도록 지정했다. 해당 환경 함수만 별도로 실행해 **PS 공통 구간의 10분 평균185.531–264.906K, 무응력 기준219.637899K**를 추출했다. 이 기준에서19개 재료의 열입력·재료점 수치검증을 수행했다. 무응력 기준은 기존 코드가 평균 기온으로 정한 가정이며 실측 조립온도가 아니다.

## 이번에 설정한 내용

`catalog_with_thermal.json`은 기존 19개 기계 재료에 `thermal` 객체를 추가한 버전이다. 기존 검증 패키지의 E·강도·경화 등은 바꾸지 않았다. `thermal_catalog.json`에는 열물성, 출처, 기준온도, 적용 범위와 가정을 별도 관리한다. 기계 기준 카탈로그가 바뀌면 SHA-256 검사에서 재검증을 요구한다.

| 재료 record | 선팽창 입력, 10⁻⁶/K | 근거 상태 |
|---|---:|---|
| S355 시편 보정형 | 12 | 기존 돔의 관용값, 열시험 적합 아님 |
| S355 명목형 | 12 | 구조용 강재 관용 참고값 |
| Al6061-T6 | NIST α(T), 4–300K | 발표된 길이 변화 함수로 계산 |
| C30 | 10 | 보통중량 콘크리트 관용값, 배합별 미보정 |
| SD400 | 12 | 일반 강재값 전용, 제품시험 아님 |
| CFRP AS4/8552 UD | α1≈0, α2=α3=32.61 | 방향별 경화공정 연구 모델값 전용, 저온 실측값 아님 |
| 소결 JSC-1A 고·저공극률 2개 | 각8, 비교4/8/12 | 현무암계 대용 가정, 공극률별 CTE 미식별 |
| 달 JSC-1A 밀도별 peak/residual 6개 | 각8, 비교4/8/12 | 유효 벌크 가정, 열·다짐 미보정 |
| 화성 JEZ-1 peak/residual 2개 | 각8, 비교4/8/12 | 유효 벌크 대용 가정, JEZ-1 실측 아님 |
| 지구 모래 포화밀도형·건조형 2개 | 각11, 비교5.5/11/16.5 | 골격 대용 가정, 열간극수압 제외 |
| 느슨한 JSC-1A 차폐층 | 8, 비교4/8/12 | MassOnly 장부용, 열력·강성 부여 없음 |

레골리스의 입자 팽창과 조립체의 유효 벌크 팽창은 다르다. peak/residual은 강도 경계 차이이므로 같은 재료군의 열입력을 동일하게 둔다. 공극률·밀도·수분·구속압·열주기별 실제 CTE는 추가 시험이 필요하다. 가정값에 대한 수치 검사를 통과해도 실험 보정으로 승격하지 않는다.

## 열변형 적용 방법

`thermal.py`는 `기계변형률 = 전체변형률 - 열변형률`을 계산한다. Al은 NIST의 L(T)/L(Tref)-1을 쓰며, α(T)=L'(T)/L(T)이다. CFRP는 재료축과 적층 방향을 구분하며 CLT의 열막력·열모멘트 계산도 제공한다. 차폐재는 질량을 보존하고 가상의 쉘 강성이나 열력을 만들지 않는다.

`verify_thermal.py`는 변환한 기계변형률을 기존 실제 OpenSees Steel02·Concrete02·DP·PDMY02·탄성·직교이방 재료에 전달한다. 3D 재료점은 SSPbrick의 아핀 변위 패치로 실행한다. 재료점 시험이며 전체 SSI의 요소별 온도 연성은 아니다.

현재 OpenSees3.4.0에는 확인한 3D `InitStrain` 래퍼가 없다. `thermal.require_native_3d_thermal_support`는 이를 명시적인 오류로 처리한다. 기존 통합본에 계수만 추가했다고 열연성이 활성화됐다고 표시하지 않는다. 3.4의 uniaxial InitStrainMaterial은 초기 상태를 commit하므로 비선형 열주기를 자유팽창과 함께 적용하는 전역 FE 경로도 자동 승인하지 않는다.

## 기존 구조물 설정과 재검증 근거

- `extract_structure_temperature.py`는 기존 구조물의 환경 함수만 AST로 추출해 그대로 실행한다. OpenSees 모델이나 질량 최적화는 실행하지 않는다. 소스·TWINS·PS 해시를 `structure_environment.json`에 저장했다.
- 기존 `THERMAL_REF_TEMP_K=None`은 TWINS/PS 공통 구간의600초 평균·결측 보간 후 전체 평균을 쓴다. 결과는185.5312308–264.9059608K, 기준219.6378991K,142개 환경 표본이다. 원시 센서 필터160–300K는 해석 온도 범위가 아니다.
- 돔은 ΔT에 exp(-덮개두께/0.05m)를 곱한다. 덮개0m는 원래 변동을 받고,2m와7.5m는 이 가정식에서 거의0이다. 이는 단순 감쇠 가정이며 열전도 검증 결과가 아니다. 기존 매립 링은 ΔT=0 가정, 지반의 깊이별 온도는 미정이다.
- `structure_environment_temperature.csv`가 기존 코드와 일치하는10분 평균 이력이다. 모든 재료가 이 온도를 실제로 받는다고 확정하지 않고, 공통 재료점 비교 입력으로 사용했다.
- `results_structure_v1`:19개 record,99개 수치 검사 통과. 자유팽창, 미소 구속 가열의 부호·단위, 온도주기, 온도 보상 후 기존 기계 이력 재현, CFRP 방향/대칭 적층, 질량 보존, PDMY stage0/1의 실제 전단응답 차이 확인.
- `results_v1`: 초기 PDMY 단계 설정 경고가 있어 합격 근거에서 제외했다. `REJECTED.md` 참고.
- `structure_history_free_v1`, `structure_history_restrained_v1`:19개에 기존 구조물 환경 이력을 가한 자유/구속 재료점 계산, 각142개 환경 표본.
- `structure_history_alpha05_v1`, `structure_history_alpha15_v1`:PROVISIONAL 가정 계수만0.5/1.5배로 바꾼 같은 구속 조건의 민감도 계산. 원본 계수1배와 비교한다.
- 이전 `results_v2`, `results_v3_provisional`, `mars_history_*`는 원시 센서 범위와5분 평균을 사용한 예비 결과다. 최종 구조물 설정의 결과로 혼용하지 않는다. `provisional_snapshot`에 당시 카탈로그를 보존했다.

기존 기계 계수는 새 온도시험에 맞춰 재식별하지 않았다. 온도별 E·항복/압축강도·손상·크리프·열압밀·얼음/물 상변화는 이번 검증에서 보정되지 않았다.

## 기존 구조물 설정으로 재실행

```sh
/Users/sydney/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 extract_structure_temperature.py
/Users/sydney/anaconda3/envs/cap/bin/python prepare_catalog.py --structure-environment structure_environment.json
/Users/sydney/anaconda3/envs/cap/bin/python verify_thermal.py --output-dir <새결과폴더>
```

## 별도 숫자 온도로 변경할 경우

숫자는 모두 K이다. 아래 옵션에 사용자가 지정한 실제 최저·최고·무응력 기준온도를 넣는다. 허용 범위 밖이면 임의 외삽하지 않고 새 물성 근거를 준비한다.

```sh
/Users/sydney/anaconda3/envs/cap/bin/python prepare_catalog.py --minimum-temperature-K <최저K> --maximum-temperature-K <최고K> --reference-temperature-K <기준K>
/Users/sydney/anaconda3/envs/cap/bin/python verify_thermal.py --output-dir <새결과폴더>
```

실제 지정 온도이력이 있다면 별도의 부재·토층 온도 입력으로 연결해야 한다. 기존 기온을 모든 재료 온도로 동일 적용하는 것은 공통 재료점 연구 시나리오의 가정일 뿐이다.

## 주요 출처

- JRC, [Handbook3 IV-7](https://eurocodes.jrc.ec.europa.eu/sites/default/files/2021-12/handbook3.pdf)
- NIST, [6061-T6 저온 길이 변화 함수](https://trc.nist.gov/cryogenics/materials/6061%20Aluminum/6061_T6Aluminum_rev.htm)
- FHWA, [콘크리트·골재·강재 열팽창 참고자료](https://www.fhwa.dot.gov/publications/research/infrastructure/pavements/pccp/thermal.cfm)
- [AS4/8552 경화공정 연구모델 Table1](https://www.mdpi.com/2504-477x/5/2/63)
- OpenSees, [현행 3D InitStrain 정의](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/InitStrain.html)
