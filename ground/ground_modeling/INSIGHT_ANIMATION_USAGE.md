# 빠른 비선형 확인

반복 작업에는 `--insight-local --quick --no-open`을 사용하세요. 비선형 재료와 원본 필터/시간 간격은 유지하고, 기본 40 × 40 m / 9 × 9 절점 및 최대 가속도 주변 10초를 선택합니다. 형상·시간 구간을 축소한 확인용 결과이며 전체 기록이나 영역 수렴 해석을 대신하지 않습니다. 명시한 영역/격자/관찰 길이는 우선합니다.

기본 출력은 `insight_results_quick/`이며 전체 결과와 분리합니다. `--quick`을 빼면 기존 100 m / 17격자 및 전체 기록 설정입니다. 동적 로그에 경과시간과 단계당 관측 속도 기반 잔여시간 추정치를 출력합니다.

---

# 현재 비선형 재료 설정 (2026-11-04 후속 변경)

`ground_modeling_insight_updated.py --insight-local`은 이제 카탈로그의 `mars_insight_duricrust_proxy`, `mars_jez1_peak`를 DruckerPrager 비선형으로, `basalt_elastic_analogue`를 탄성 기반암으로 사용합니다. `--dynamic-material-model` 기본은 `nonlinear`입니다. 아래 과거 문서의 ‘기본 탄성’ 및 ‘1 km² 기본 영역’ 설명은 현재 설정에 적용되지 않습니다.

새 기본 확인 영역은 100 × 100 m, 17 × 17 절점입니다. 0.2 m 표층을 포함하도록 변경했으며 메시·영역 수렴이 확보된 최종 설계 영역은 아닙니다. 중력 하중계수와 수렴 진행을 출력합니다. 실제 재료는 `insight_materials_used.json`에 기록합니다. 재료 변경 후 기존 HTML은 재실행 전까지 이전 결과입니다.

빠른 실행은 `--dynamic-max-steps 100 --no-open`을 추가하고, 선형 비교가 필요할 때만 `--dynamic-material-model elastic`을 명시합니다. 화성 재료의 문헌 강도와 강성 가정을 연결한 것으로 현장 반복 지진 보정 완료를 뜻하지 않습니다. 전체 재료 적합성 보고서는 `../../material/재료_튜닝_적합성_및_InSight_적용_보고서_20261104.md`에 있습니다.

---

# 2026-11-04 지진 입력 교체: S1222a

현재 기본 동적 입력은 [S1222a_ACC.mseed](../../structure/MARS_DATA/S1222a_ACC.mseed)입니다. 아래의 2019년 UVW 설명은 `--insight-mseed`로 직접 지정하는 이전 입력에만 적용됩니다.

```bash
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling.py --insight-local --no-open
```

`ground_modeling.py`의 정적 기본 실행은 유지하며, `--insight-local`을 지정하면 동적 실행 코드로 연결합니다. `ground_modeling_insight_updated.py --insight-local`도 동일한 새 입력을 사용합니다.

- 새 기록: XB.ELYSE.02.BHE/BHN/BHZ, 2022-05-04 23:13:07.901 UTC부터 약 3,000초, 20 Hz, 성분별 59,999샘플.
- `ACC` 제품을 기기보정 완료된 m/s² 가속도로 가정합니다. MiniSEED 자체로는 물리 단위를 검증할 수 없습니다. 응답 보정과 UVW 회전은 적용하지 않습니다.
- 선형 추세 제거, 5% taper, 0.05–2 Hz 필터, 10 Hz 재표본화 후 ENU 순서의 CSV로 저장합니다. 필터 및 표본화 옵션은 기존과 같습니다.
- `--preprocess-only`로 전처리만 실행할 수 있습니다. `--mars-acc-mseed PATH`로 다른 보정 완료 ENZ 가속도 파일을 지정할 수 있습니다.
- PS/TWINS CSV는 기압·온도·바람 자료이며 이번 지진 입력에는 합치지 않습니다. 지진과 다른 관측 시기의 자료입니다.
- 지표 관측을 기반 입력으로 사용하는 기존 실험용 해석 및 DEM 설정은 유지합니다.

---

# HiRISE IMG + InSight MiniSEED 동적 해석

`ground_modeling_insight_updated.py`에서 지정한 IMG와 같은 폴더의 BHU/BHV/BHW 기록을 직접 연결합니다.
원래 Earth/Moon/Mars 정적 실행 옵션도 유지합니다. 수정 전 코드는 [백업 폴더](../archive/ground_modeling_20261116/backups/)에 있습니다.

## 실행

이 컴퓨터에서 필요한 패키지가 설치된 환경은 `/Users/sydney/anaconda3/envs/cap/bin/python`입니다.

```bash
cd /Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling_insight_updated.py --insight-local
```

`--insight-local`은 다음을 선택합니다.

- 화성, `DTEEC_037262_1845_036761_1845_U01.IMG` 영상 중앙의 1,000 × 1,000 m 영역 (1 km²)
- 수평 41 × 41 절점, 25 m 간격의 초기 확인용 격자 (기존 150 m 모델은 12.5 m 간격)
- 같은 폴더의 `xb.elyse.03.bhu/bhv/bhw.mqs2019drgb.1.mseed` 세 파일
- `stationxml.xb.elyse.2019.151.xml`의 기기 응답과 센서 방향
- 선형 추세 제거, taper, 기기보정(가속도 m/s²), 0.05–2 Hz 필터, 10 Hz 재표본화, UVW → ENU 회전
- 기존 기준 G/K를 E/ν로 환산한 선형 탄성 재료 (`--dynamic-material-model elastic`)
- 중력 초기화 후 고정 바닥/자유 측면의 실험용 UniformExcitation 해석
- 기본 원기록 전체 약 1시간(10 Hz에서 약 36,014단계), HTML은 약 180개 구간으로 자동 저장; 초기·최종 프레임 포함

HTML이 생성되면 자동으로 브라우저에서 열립니다. 자동 열기를 끄려면 `--no-open`을 붙이세요.
기본 저장 위치: `insight_results/` (`--output-dir`로 변경 가능).

- `insight_dynamic_animation.html`: Plotly를 포함하는 오프라인 HTML
- `insight_acceleration_enu.csv`, `.json`: 보정 가속도와 처리 메타데이터
- `insight_experimental_dynamic_response.csv`: 중앙 지표 절점의 변위 및 상대 가속도

## 다른 MiniSEED 세트 / 해석 시간

동일 관측소·이벤트의 U/V/W 세 성분과 그 관측 시각에 유효한 StationXML을 지정합니다.
임의의 단일 성분 또는 ENZ 형식은 이 전처리 경로의 대상이 아닙니다.

```bash
/Users/sydney/anaconda3/envs/cap/bin/python ground_modeling_insight_updated.py \
  --insight-local \
  --insight-mseed /absolute/path/event.bhu.mseed /absolute/path/event.bhv.mseed /absolute/path/event.bhw.mseed \
  --insight-stationxml /absolute/path/event.station.xml \
  --no-open
```

기본은 원기록 전체입니다. `--dynamic-max-steps`를 지정하면 해석 단계 수를 제한합니다.
`--dynamic-scale`은 물리 입력을, `--dynamic-visual-scale`은 표시 배율만 바꿉니다.
필터/표본화는 `--motion-freqmin`, `--motion-freqmax`, `--motion-rate`로 변경합니다.

빠른 연결 검증은 `--local-nx 5 --local-ny 5 --dynamic-max-steps 10 --dynamic-animation-stride 1`을 사용하세요.
전처리만 실행하려면 `--preprocess-only`를 추가합니다.

## HTML 조작

- 시간 슬라이더, 재생/일시정지, 처음, 0.02×–10× 재생 속도
- 동적 해석 시작 이후 변위 / 중력 포함 전체 변위 선택
- 변형 표시 배율 직접 입력, 실제 크기 ×1, 자동 확대
- 재료 표시 켜기/끄기, 3D 회전/확대, 절점별 실제 변위 hover
- 저장 프레임의 전체 절점 최대 변위 그래프와 현재 시간선, 최대 프레임 이동

재생 시간 간격은 저장 프레임의 실제 시간 차이를 따릅니다. 브라우저 처리 속도에 따라 느려질 수 있습니다.
프레임 사이에 발생한 응답 극값은 HTML 그래프에서 누락될 수 있습니다.
표시 변위는 OpenSees의 상대 절점 변위이며 지반 입력 이동을 더한 절대 이동이 아닙니다.

## 지형 위치와 모델 범위

PDS3 헤더의 자료형, raster offset, 배율, 결측, 해상도를 읽고 메모리 매핑과 bilinear 보간으로 필요한 영역만 추출합니다.
X/Y/Z는 동/북/상이며, 기준 고도는 추출 영역 평균입니다.
`--dem-center-row`, `--dem-center-col`은 0부터 시작하는 영상 행/열입니다. `--local-width/height`는 m입니다.
기본 영역은 **영상 중앙**입니다. 실제 SEIS 관측소 위치와의 일치는 보장하지 않습니다.

지층/재료는 기존 코드의 화성 가정 물성입니다. IMG는 표면 고도만 제공하며 지층 정보를 제공하지 않습니다.
`--insight-local`은 SEIS 지표 관측을 기반 입력으로 사용하는 실험을 선택합니다. 현재 경계는 파를 반사하고 별도 감쇠 보정도 없으므로 정량적인 현장 지반 응답 검증용으로 사용하려면 경계·물성·감쇠·입력 위치를 별도로 보정해야 합니다.

## 비선형 재료

`--dynamic-material-model nonlinear`로 기존 PDMY02를 유지할 수 있습니다. 현재 자료의 자유 측면 모델에서는 소성 단계 전환 직후 첫 동적 단계의 수렴 실패가 확인되어 기본값은 선형 탄성입니다. 비선형 옵션을 사용할 때는 중력 후 소성 평형·시간 간격·재료 보정 검토가 필요합니다. 실패한 해석은 HTML로 저장하지 않습니다.

영역 확대 후 기본 `--insight-local`은 1 km²를 사용합니다. 기존 150 m 결과 HTML은 재실행 전까지 바뀌지 않습니다. 수평 해상도 12.5 m를 유지하려면 `--local-nx 81 --local-ny 81`을 지정할 수 있지만 계산량과 메모리가 크게 증가합니다. 면적 확대만으로 해의 수렴성이나 경계 독립성이 확보되는 것은 아닙니다.

## 관찰 시간 선택

기본 `--insight-local`은 원기록 시작부터 끝까지 전체 약 1시간을 해석합니다. 현재 기록은 0–3601.4초, 36,015샘플입니다. HTML 시간 0도 원기록 시작입니다. 구간 양끝에는 5% cosine taper가 적용됩니다.

- 전체 약 1시간: `--insight-local`
- 최대 가속도 전후 총 300초만: `--insight-local --motion-duration 300`
- 원기록 2600초부터 300초만: `--insight-local --motion-start 2600 --motion-duration 300`

`--dynamic-max-steps`는 단계 수를 줄이는 옵션이므로 전체 기록을 보려면 생략하세요.
해석은 매 시간 단계마다 수행하며 HTML 저장 프레임은 기본 약 180개 구간으로 제한합니다. 프레임 사이의 응답 극값은 생략될 수 있습니다. 저장 간격은 `--dynamic-animation-stride`로 직접 설정할 수 있습니다.
원자료는 `insight_acceleration_enu.csv`, 선택 입력은 `insight_acceleration_window_enu.csv`에 저장합니다.
1 km² 모델의 약 36,014단계 해석은 계산량이 큽니다. 기존 HTML은 전체 해석이 성공적으로 완료된 뒤 갱신됩니다.
