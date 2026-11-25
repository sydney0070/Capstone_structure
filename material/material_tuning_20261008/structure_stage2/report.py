"""Render Korean summary from actual saved results, without invented values."""
import json
from model import BASE

if __name__=='__main__':
    results=json.loads((BASE/'comparison.json').read_text())
    verification=json.loads((BASE/'verification.json').read_text())
    lines=['# 2단계 구조물 비교 결과','',
        '임시 치수·보정 전 재료·동일 두께의 **탄성 검증용 결과**다. 행성별 지반 특성은 동일하게 두고 중력만 바꿨다. 실제 재료 순위·설계 안전성을 뜻하지 않는다.','',
        f'자중 시나리오 {len(results)}개 실행. 검증 {len(verification["checks"])}개 중 통과 {sum(x["passed"] for x in verification["checks"])}개.',
        f'전체 자중 시나리오의 최대 상대 힘 평형 잔차: {max(r["relative_force_residual"] for r in results):.3e}.','',
        '## 알루미늄 기준: 최대 변위와 지지 유연성','',
        '| 형상 | 지지 | 지구 (mm) | 달 (mm) | 화성 (mm) | 질량 (tonne) |',
        '|---|---|---:|---:|---:|---:|']
    for shape in ('dome','hyperboloid'):
        for support in ('pinned','linear_springs'):
            rows=[next(r for r in results if (r['shape'],r['support'],r['planet'],r['material'])==
                   (shape,support,p,'shell_aluminum')) for p in ('earth','moon','mars')]
            lines.append(f'| {shape} | {support} | '+ ' | '.join(f'{r["max_displacement_m"]*1000:.6f}' for r in rows)+f' | {rows[0]["mass_tonne"]:.4f} |')
    lines+=['','질량은 행성이 바뀌어도 같다. 이 선형·자중 전용 모델의 변위는 중력에 비례한다. 이는 코드가 만족해야 할 관계이며, 실제 행성 SSI도 선형 비례한다는 결론은 아니다.','',
        '## 돔 압력차 1 kPa 보조 예제','',
        '알루미늄 돔·선형 스프링·자중과 외향 압력의 조합이다. 아래 인장 지지력은 앵커가 구현되었다는 뜻이 아니며, 무인장 접촉 모델로 옮기기 전에 확인할 요구량이다.','',
        '| 행성 | 자중 (kN) | 요구 인장 지지력 합 (kN) |','|---|---:|---:|']
    for p in ('earth','moon','mars'):
        r=json.loads((BASE/'results'/f'pressure_demo_{p}'/'summary.json').read_text())
        lines.append(f'| {p} | {r["weight_kN"]:.3f} | {r["required_hold_down_kN"]:.3f} |')
    lines+=['','중력이 작아져도 압력 하중은 줄지 않는다. 이 예제에서는 가벼워진 자중만큼 더 큰 인장 지지가 요구된다. 바닥판이 없는 돔 캡 모델이므로 완전한 밀폐 압력용기의 지반 순하중으로 일반화하지 않는다.','',
        '## 메시 검증','',
        '알루미늄 대표 조건에서 n=12/24/48/72를 비교했다. 스프링 지지는 n=24→48 차이가 약 4%여서 비교용 메시를 n=48로 올렸다. 최종 수용 기준은 n=48→72 최대 변위 차이 2% 미만으로 정했으며, 이는 이번 기준 모델의 검토 기준이다.','',
        '| 형상 | 지지 | n=48→72 상대 차이 |','|---|---|---:|']
    for row in verification['mesh']:
        a,b=row['rows'][-2:];err=abs(a['max_displacement_m']/b['max_displacement_m']-1)
        lines.append(f'| {row["shape"]} | {row["support"]} | {err*100:.4f}% |')
    lines+=['','메시 검증은 알루미늄·현재 하중·현재 기초 강성의 대표 조건이다. 재료, 형상, 집중하중, 접촉, 기초 강성 또는 비선형 모델이 달라지면 다시 검증해야 한다. 응력 집중·좌굴 모드의 수렴은 이번 변위 검증에 포함되지 않는다.','',
        '## 남은 입력','',
        '실제 형상·기초·압력차를 먼저 확정하고, 보정할 실험자료를 `experiment_inventory.json`에 정리한다. 그다음 기존 튜닝 코드로 물성을 보정하고 비선형 접촉·좌굴 모델로 확장한다.']
    (BASE/'RESULTS_KO.md').write_text('\n'.join(lines)+'\n')
