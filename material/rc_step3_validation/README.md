# RC Step 3 — 물리적 검증 1차 확장

- 기존 `test_rc_step3.py`, `rc_materials.py`, `rc_section.py`를 동봉. 기존 결과 폴더를 덮어쓰지 않고 `results_rc_validation/`에 별도 저장.
- `python validate_step3.py`
- 포함: 축하중 −360 kN 재하 후 지점반력 평형, `zeroLengthSection` 및 `forceBeamColumn`의 상하부 콘크리트·철근 Fiber `stressStrain`, 60 vs 120 분할 횡변위 해석 및 요약 비교.
- `axial_equilibrium`에서 압축축력의 지점반력과 상대 평형 오차 검사.
- 반복가력에는 축하중을 포함하지 않음. 제로횡변위의 최종 힘은 잔류변위가 아님. 정확한 잔류변위와 완전한 소산에너지는 다음 검증 대상.
- Fiber recorder가 빈 파일을 내는 경우 OpenSees 버전과 단면 호환성을 확인해야 함.
- 출력: `results_rc_validation/validation_summary.json`, 두 분할가력의 CSV, Fiber 텍스트.
- OpenSees 문서: https://opensees.github.io/OpenSeesDocumentation/user/manual/output/ElementRecorder.html
- **설계 안전성 검증이나 실험 데이터 Calibration이 아닌 수치 검증용 프로토타입.**
