추천 시편 적용 지반 — 2026-11-24

comparison.html: 기존 완료 MGS-1 결과와 물성·침하·평형·소성 상태 비교
run/footing_settlement.html: 감압 상태3D 결과
run/footing_settlement_pressurized.html: 내압 적용 후3D 결과
run/summary.json: 완료 해석 원본
model/ground_geometry.html: 새 재료 배치
source_snapshot: 실행 당시 코드·재료 프로파일 사본

현재 ground_dome_footing.py 기본값은 recommended_20261124.
경화층DP(ρ1025), 모래층MGS-1 DP(ρ1300), 조립질층·암반층탄성.
기존 완료 MGS-1과 실제 물성 차이는 경화층 밀도1200→1025뿐.
MGS-1의 문헌 전단강도 보정은 아직 수행하지 않음.
침하는 각 해석의 굴착 지반 자중 평형 이후 증분.
내압 적용 후에도 순압축력29.58MN이 남으며 무하중 잔류침하가 아님.
수렴·힘 평형 확인은 지지력/접촉/들림 검증을 뜻하지 않음.
