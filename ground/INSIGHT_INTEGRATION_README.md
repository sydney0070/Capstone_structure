# InSight SEIS → OpenSeesPy 연결

## 제공 파일
- `ground_modeling_insight_updated.py`: 기존 `ground_modeling(5).py`를 기반으로 구성한 단독 실행 버전
- `insight_motion.py`: MiniSEED 헤더 검사, StationXML 기기 응답 보정, ENU 회전, 동적 시험용 입력

원본 스크립트의 기본 진입점은 별도의 `ground_hybrid` 패키지를 요구했으나 첨부되지 않았습니다. 수정본은 기존 `main()`을 직접 사용합니다. 이로 인해 프로젝트의 `ground_hybrid` CLI 기능까지 보존되는 것은 아닙니다.

## 의존성
```bash
python -m pip install numpy matplotlib plotly scipy obspy openseespy rasterio
```

## 1. 다운로드한 세 MiniSEED 파일 확인
```bash
python ground_modeling_insight_updated.py --inspect-insight-mseed \
  'xb.elyse.03.bhu.mqs2019drgb.1(1).mseed' \
  'xb.elyse.03.bhv.mqs2019drgb.1(1).mseed' \
  'xb.elyse.03.bhw.mqs2019drgb.1(1).mseed'
```
이 기능은 ObsPy 없이 MiniSEED v2의 첫 레코드 헤더만 검증합니다. 파형 샘플을 복호화하는 것은 아닙니다.

## 2. StationXML을 확보한 후 가속도 변환
StationXML은 해당 XB.ELYSE.03.BH? 채널의 관측 시각을 포함하고 응답계수 및 방위각·경사각을 제공해야 합니다.
```bash
python insight_motion.py --mseed \
  'xb.elyse.03.bhu.mqs2019drgb.1(1).mseed' \
  'xb.elyse.03.bhv.mqs2019drgb.1(1).mseed' \
  'xb.elyse.03.bhw.mqs2019drgb.1(1).mseed' \
  --stationxml insight_station.xml \
  --output insight_acceleration_enu.csv
```
출력: 시간, east/north/up 가속도(m/s²) CSV 및 보정 메타데이터 JSON. 필터(기본 0.05–2 Hz)·회전·기기응답 검증은 후속 분석 과제입니다.

## 3. 실험용 지반 가진 (최종 공학해석으로 사용하지 마세요)
현재 모델의 화성 지형은 MOLA 파일 중앙 영역을 추출하므로 실제 InSight 착륙지 지형이 아닙니다. ENU 정렬은 DEM을 교정한 다음에만 인정할 수 있습니다.
```bash
python ground_modeling_insight_updated.py --body Mars \
  --dynamic-csv insight_acceleration_enu.csv \
  --experimental-free-sides \
  --experimental-surface-as-base \
  --axes-enu-verified \
  --dynamic-max-steps 20 \
  --no-html --no-open
```
마지막 명령은 다음 세 조건을 **사용자가 별도로 확인**한 경우에만 사용하십시오.
- `--axes-enu-verified`: FEM x/y/z 방향이 동/북/상에 실제로 맞음
- `--experimental-surface-as-base`: 표면 관측파형을 기반암 운동으로 가정하며 결과를 검증된 지반 응답으로 해석하지 않음
- `--experimental-free-sides`: 3D free-side boundary의 인공 반사 문제 수용

**주의**
1. 관측된 SEIS는 지표면 파형이므로 기반암 입력운동으로 직접 사용하면 파동전파/증폭을 중복 반영할 수 있습니다.
2. 기존 `bbarBrick`에는 구성재료의 `rho`가 이미 들어가며 질량행렬을 형성할 수 있으므로 같은 질량을 다시 `ops.mass`로 주지 않았습니다.
3. 현재 환경에는 `obspy`, `openseespy`가 없어 응답 보정과 동적 수치해석을 재현하지 못했습니다.
4. OpenSees에 넣는 지진파의 주파수 콘텐츠 대비 요소 크기와 동적 경계조건을 검증해야 합니다.
5. `DTEEC` 지형과 Caltech 레골리스 분포를 연동하는 코드는 아직 이 버전에 포함되지 않았습니다.
