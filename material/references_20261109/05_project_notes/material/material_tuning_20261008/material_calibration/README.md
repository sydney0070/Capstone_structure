# 실제 실험자료를 이용한 재료 보정

후속 작업: [재료별 material·element와 최신 보정 결과](MATERIAL_MAP_KO.md). [2026-10-12 재료 전용 진행 결과](PROGRESS_20261012_KO.md)에 철근 재보정과 알루미늄 식별성 점검을 추가했다. 철근은 기존 검증 자료의 수치 기준을 통과했으나 새 독립 검증이 필요하다. 기존 결과와 원자료는 보존했다. JEZ-1의 문헌 강도 입력은 독립 곡선 보정과 구분한다.

현재 완료한 항목은 **S355 J2+M HEM320 웹의 Steel02 단축 시편 보정**이다. 공개 실험 2개로 학습하고 별도 실험 2개로 검증했다. 다른 재료들은 [현황표](MATERIAL_STATUS_KO.md)에 자료 확보와 검증 상태를 구분했다. 기존 23개 예시 설정 전체가 보정된 것은 아니다.

## 실행

상위 패키지의 requirements.txt를 설치한 Python에서, 패키지 루트를 작업 디렉터리로 사용한다.

```bash
python material_calibration/prepare_steel.py
python tune.py fit --config material_calibration/s355_hem320_web/starting_config.json --targets material_calibration/s355_hem320_web/targets.json --maxiter 100 --out material_calibration/s355_hem320_web/fitted
python material_calibration/assess_steel.py
```

동일 경로의 결과가 갱신된다. 결과만 검토하려면 다시 최적화할 필요가 없다. OpenSees 3.8.0에서 실행했으며 난수 시드는 20261008이다. 준비 코드는 전체 원자료 ZIP 또는 함께 제공하는 선택 시편 원본 4개를 읽는다. 원본 SHA256을 검증하고 변환한다.

최종 매개변수는 `s355_hem320_web/reviewed_config.json`, 오차·감도·증분 검토는 `assessment.json`, 사람이 읽는 결과는 [RESULTS_KO.md](s355_hem320_web/RESULTS_KO.md)에 있다. 이 설정을 `tune.py run --config`에 전달해 추가 단축 재하 경로를 계산할 수 있으나, 검증 범위가 자동 확대되는 것은 아니다.

## 실험과 모델의 대응

원자료의 로그변형률을 `expm1(e_true)`로, 진응력 MPa를 `Sigma_true / exp(e_true) * 1000`으로 변환했다. 균일 변형 가정의 공칭 응력 kPa·공칭 변형률이며, 처음 절대 변형률 2%를 넘기기 직전까지만 사용한다. 최초 작은 오프셋을 제거하고 취득 순서는 유지했다. 저자가 처리·다운샘플링한 실험자료이며 합성자료가 아니다.

학습 시편별 응력 범위로 잔차를 정규화하고 각 시편 평균제곱오차에 동일 가중치를 준다. 측정점 수가 많은 시편에 자동으로 더 큰 가중치를 주지 않는다. 검증 시편은 매개변수 최적화에 사용하지 않았다. 사전 기준은 각 검증 시편 NRMSE 5% 이하, 인장·압축 극값 오차 10% 이하이다. 설계 규격의 합격 기준은 아니다.

E, Fy, b, R0, a1, a3 여섯 값을 보정했다. 나머지 Steel02 계수는 설정 파일의 고정값이며 식별했다고 주장하지 않는다. 단축 강재 모델을 다축 쉘 J2 모델이나 철근으로 그대로 옮겨서 보정 완료라고 할 수 없다. 검증 곡선에는 항복 평탄부 등의 국부 차이가 남아 있다.

## 출처와 배포

선택 원본은 Hartloper, Ozden, de Castro e Sousa, Lignos (2022), [Database of Uniaxial Cyclic and Tensile Coupon Tests for Structural Metallic Materials](https://doi.org/10.5281/zenodo.6965147), v1.0.0, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)이다. 실험 원 연구는 [El Jisr (2022), EPFL thesis](https://doi.org/10.5075/epfl-thesis-8700)이다. 원본 `*_source.csv`와 변환 결과 `*.csv`를 구분했고 변경 내용은 `provenance.json`에 기록했다. 저자들의 승인·보증을 의미하지 않는다.

63 MB 전체 원자료 ZIP은 경량 배포본에 포함하지 않았다. 선택 원본 4개, 메타데이터, 참고문헌은 포함하므로 위 실행을 재현할 수 있다. `raw/stressStrain.zip`과 달 자료는 아직 보정에 사용하지 않은 조사 자료이며 각각의 출처와 한계는 현황표를 따른다.
