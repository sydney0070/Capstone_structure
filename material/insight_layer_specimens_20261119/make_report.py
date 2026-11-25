"""Create an offline specimen/source/result index, without rerunning analysis."""
import csv
import hashlib
import html
import json
import os
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
esc=lambda x:html.escape(str(x),quote=True)


def main():
    records=json.loads((HERE/'materials.json').read_text())['materials']
    summary=json.loads((HERE/'results/summary.json').read_text())
    results={r['id']:r for r in summary['cases']}
    sources=json.loads((HERE/'source_manifest.json').read_text())
    # These are source observations/reference values, separate from the chosen inputs.
    raw=[
        ['Spohn2022','Figure 20(d), Section 6.7','local shallow boundaries','cm','0–1 dust/sand; 1–19 duricrust; 19–31 sand; gravel below31 within penetrated interval','local interpretation; not a borehole to metres'],
        ['Spohn2022','Section 6.7','duricrust density','kg/m3','950–1100','inferred range; specimen midpoint1025 is a choice'],
        ['Spohn2022','Section 6.7','sand density','kg/m3','1300–1500','inferred; 1300 selected'],
        ['Spohn2022','Section 6.7','gravel mixture density','kg/m3','about1600','15 vol-% rocks estimate'],
        ['Spohn2022','Section 6.4','duricrust cohesion','Pa','about5800 at assumed phi30 deg, rho1200 kg/m3','conditional slope inference, transferred to selected density'],
        ['Spohn2022','Section 6.4','duricrust cohesion range','Pa','2000–15000','penetration-model dependent and conditional on phi30–40 deg'],
        ['Spohn2022','Sections 6.2–6.4','penetration resistance','Pa','crust350000–700000; deeper4900000 ±400000','NEVER used as E or UCS'],
        ['Brinkman2022','Table 4, PDF p.13','effective shallow E','MPa','11.48 (+5.91 / -2.23)','wave-based effective value, not static isolated crust'],
        ['Brinkman2022','Table 4, PDF p.13','effective shallow nu','1','.28 (+.12 / -.051)','selected with E, derive G/K consistently'],
        ['Brinkman2022','Table 4, PDF p.13','reported modal G and K','MPa','G4.47 (+2.00/-.83); K7.79 (+1.60/-1.55)','independently summarized modes; do not combine all modes as one consistent elastic tuple'],
        ['Brinkman2022','Results and Table 4','Vp; Vs','m/s','Vp119 (+45/-21); Vs63 (+11/-7)','bulk shallow volume, not separate layer velocities'],
        ['Murdoch2018','Table 2, supplied PDF pp.15–16','predicted 1m rho,Vp,Vs,confinement','SI','rho1307; Vp157.9; Vs94.6; p5440','existing rounded E28.4MPa nu.22 rho1300 candidate retained'],
        ['Dotson2024','Figures5–6','MGS1/MGS1S peak strength at low density','deg,Pa','phi35.5/34; c≈280','previous graph reading / reported friction, approximate, transferred to1300'],
        ['Kenda2020','Table1, archived PDF p.4','expected E and nu','MPa,1','surface7.5,.22; 1m28,.22; coarse520,.24; fractured13000,.28; basalt65000,.25','reference order of magnitude, not local core calibration'],
        ['Morgan2018','PDF p.39, terrestrial analogue','fractured basalt density','kg/m3','2760','terrestrial reference used in prelanding model'],
        ['Hobiger2021','Table2, PDF p.7','sediment modified-prior Vp and Vs','m/s','750 and417','prior, not final inversion stiffness or measured30–75m point'],
        ['PDMY02','suggested parameters table','Dr30/40/50/60/75 complete reference sets','SI','see materials.json and per-parameter CSV','Earth reference; G/K80kPa scaled to101kPa using d=.5; not a Mars relative-density measurement'],
        ['Shallow2008','PDF p.17, first example','QzSimple2 spring101','N,m','Qult38020; z50.00274347; type1; suction.1; damping.05','static fixture uses zero damping'],
        ['Shallow2008','PDF p.17, first example','PySimple2 spring105','N,m','Pult102000; y50.0145067; type1; drag.1; damping.05','static fixture uses zero damping'],
        ['Shallow2008','PDF p.17, first example','TzSimple2 spring106','N,m','Tult100000; z50.0142222; type1','only documented one optional damping argument; fixture uses zero'],
    ]
    with (HERE/'literature_extraction.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['source','locator','quantity','unit','reported_value','interpretation_limit']);w.writerows(raw)
    with (HERE/'results_summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['specimen','family','group','numerical_pass','mars_calibrated','density_kg_m3','E_MPa','phi_deg','c_Pa','test_count','step_bisections'])
        for r in records:
            p=r['parameters'];s=results[r['id']]
            w.writerow([r['id'],r['family'],r['group'],s['passed'],False,p.get('density_kg_m3',''),p.get('E_Pa',0)/1e6 or '',p.get('phi_deg',''),p.get('c_Pa',''),len(s['tests']),sum(t.get('step_bisections',0) for t in s['tests'])])
    labels={'first_comparison':'첫 비교군 · 탄성 / DruckerPrager 13개',
            'followup_mars_trial':'후속 후보 · 화성 적용 PDMY02 초안 4개',
            'followup_reference':'후속 후보 · 공식 지구 모래 PDMY02 기준 5개',
            'uniaxial_comparison':'단축 후보 · 접촉·반력 스프링 6개'}
    layers={
        'local_surface_0_1cm':'국부 표면 먼지·느슨한 모래 (약0–1cm)',
        'local_duricrust_1_19cm_approx':'국부 경화층 (약1–19/21cm)',
        'local_sand_19_31cm_and_metre_scale_sandy_reference':'국부 미고결 모래 / m 규모 모래질 참조',
        'sandy_regolith_1300':'모래질 레골리스 · 기존 MGS 후보',
        'local_gravel_below_31cm_penetrated_interval_only':'국부 자갈·암편 혼합층 (약31cm 아래 관입 범위)',
        'metre_scale_blocky_ejecta_3_13m_reference':'조립질 분출물 (3–13m 참조 시나리오)',
        'fractured_basalt_13_30m_reference':'파쇄 현무암 (13–30m 참조 시나리오)',
        'deep_low_velocity_sediment_hypothesis':'심부 저속도 퇴적층 후보 · 사전모델 비교',
        'earth_sand_reference_only':'공식 지구 모래 기준 · 화성층 배정 아님',
        'footing_vertical_equivalent':'기초 밑 수직 등가 반력',
        'embedded_footing_side':'매입 기초 옆면 수평 저항',
        'footing_base_tangential':'기초 밑면 접선 미끄럼',
        'conditional_normal_contact_only':'압축 접촉/법선 갭 기능 비교',
        'conditional_duricrust_curve_fitting_only':'경화층 곡선 맞춤 기능 비교 · 합성 배면곡선',
    }
    body=['''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>InSight 층별 재료 시편</title>
<style>body{margin:36px auto;max-width:1300px;padding:0 24px;font:16px/1.65 -apple-system,BlinkMacSystemFont,"Apple SD Gothic Neo",sans-serif;color:#202c3b}h1{font-size:30px}h2{margin-top:40px}a{color:#126394}p{max-width:1060px}.note{padding:17px 20px;background:#fff2df;border-left:4px solid #b36d19}.stats{display:flex;gap:16px;flex-wrap:wrap}.stats div{padding:15px 22px;background:#edf4f7;border-radius:8px}.stats b{font-size:27px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:10px 12px;text-align:left;vertical-align:top;border-bottom:1px solid #dfe5e9}th{background:#edf2f5}code,pre{font-family:ui-monospace,monospace}pre{white-space:pre-wrap;overflow-wrap:anywhere;padding:16px;background:#f3f6f8;font-size:13px}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:20px}.card{border:1px solid #dce5eb;border-radius:8px;padding:16px;min-width:0}.card h3{font-size:16px;overflow-wrap:anywhere;margin:0 0 8px}.card img{width:100%;height:auto}.muted{color:#596a78;font-size:13px}.ok{color:#18704a}summary{cursor:pointer}details{margin:12px 0}small{font-size:12px}nav a{display:inline-block;margin:5px 18px 5px 0}</style>
<h1>InSight 층별 재료 시편</h1><p>문헌 대조 시작 2026-11-19 · 시편 작성·검증 완료 2026-11-23 · N–m–kg–s–Pa · OpenSees 3.4.0</p>
<p class="note"><strong>수치 검사 통과는 화성 실험 검증 완료를 뜻하지 않습니다.</strong> 논문 참조값, 추정·이전값, 해석 가정을 매개변수별로 분리했습니다. 모든 시편의 <code>mars_calibrated</code>는 false입니다. JEZ-1은 제외했고 기존 MGS-1/MGS-1S 저장본과 전체 지반은 교체하지 않았습니다.</p>''']
    tests=sum(len(r['tests']) for r in results.values());mass=sum('mass_test' in r for r in results.values())
    body.append(f'<div class="stats"><div><b>{len(records)}</b><br>독립 시편</div><div><b>{summary["passed_count"]}/{len(records)}</b><br>수치 검사 통과</div><div><b>{tests}</b><br>하중 이력 시험</div><div><b>{mass}</b><br>1m³ 질량 검사</div></div>')
    body.append('<nav><a href="#inputs">초기 물성</a><a href="#cases">전체 후보와 곡선</a><a href="#sources">출처</a><a href="#run">실행 방법</a></nav>')
    body.append('''<p><a href="materials.json">전체 재료 JSON</a> · <a href="parameters_and_sources.csv">매개변수별 값·출처 CSV</a> · <a href="literature_extraction.csv">문헌 원값 추출표</a> · <a href="results_summary.csv">검증 요약 CSV</a> · <a href="results/summary.json">검사 상세 JSON</a></p>
<h2>시편과 시험 범위</h2><p>3D 시편은 모두 <strong>1×1×1m 정육면체</strong>입니다. 이는 층 두께나 기존 90×90m 지반을 축소한 모형이 아니라 재료 법칙을 확인하는 단일 요소 시험입니다. 탄성·DP에는 bbarBrick, PDMY02 반복 전단에는 SSPbrick을 사용했습니다. 모든 연속체 재료의 밀도는 bbarBrick 1m³의 고유치 기반 질량으로 별도 확인했습니다.</p>
<p>탄성·DP는 구속압 5/50kPa에서 압축–제하를 수행하고 E, DP 삼축압축 항복식, 측압 유지, 균일 응력, 축력 평형을 확인했습니다. PDMY02는 수직 유효응력 5/80kPa에서 ±0.001%부터 ±0.5%까지 반복 전단했습니다. 수평 변형은 구속하고 체적 변화는 허용하므로 평균압은 수직압과 다릅니다. 간극수·액상화 시험은 아닙니다.</p>
<p>단축 후보는 두 노드 zeroLength 배치를 정의하고 동일 재료를 <code>testUniaxialMaterial/setStrain</code>으로 누적·커밋해 힘–상대변위를 기록했습니다. 전역 평형 해석이나 링 분담면적 환산을 수행한 결과는 아닙니다. Qz는 수직, Py/Tz는 수평 방향입니다.</p>
<h2 id="inputs">첫 비교군의 초기 물성</h2><div class="scroll"><table><thead><tr><th>시편</th><th>ρ kg/m³</th><th>E MPa</th><th>ν</th><th>G/K MPa</th><th>φ ° / c kPa</th></tr></thead><tbody>''')
    for r in records:
        if r['group']!='first_comparison':continue
        p=r['parameters'];strength=f'{p["phi_deg"]:g} / {p["c_Pa"]/1000:g}' if 'phi_deg' in p else '—'
        body.append(f'<tr><td><a href="#case_{r["id"]}">{r["id"]}</a></td><td>{p["density_kg_m3"]:g}</td><td>{p["E_Pa"]/1e6:.4g}</td><td>{p["nu"]:.4g}</td><td>{p["G_Pa"]/1e6:.4g} / {p["K_Pa"]/1e6:.4g}</td><td>{strength}</td></tr>')
    body.append('''</tbody></table></div><p class="muted">경화층 E·ν는 천부 유효 파속 자료를 이전한 값이고 ρ1025는 추정 범위의 중간 선택값입니다. c5.8kPa는 원 논문에서 φ30°, ρ1200을 가정한 추정값입니다. 국부 자갈층 E는 모래 기질의 대용값이며, 표면·자갈·분출물 DP의 c0/φ30–35°는 시험용 가정입니다. 암반·퇴적층의 소성 강도는 만들지 않았습니다.</p>
<h2>출처값과 입력값의 구분</h2><ul><li>Brinkman 표4에서 E11.48MPa와 ν.28을 선택해 G4.484MPa, K8.697MPa를 일관되게 계산했습니다. 표의 독립적인 최빈 G/K 값 전체를 동시에 조합한 것이 아닙니다.</li><li>MGS-1/MGS-1S의 E28.4MPa·ν.22·ρ1300과 기존 c/φ를 보존했습니다. 모사토 저압 직접전단값을 이 밀도·강성에 이전한 초기값이며 고압 강도 보정은 아닙니다.</li><li>PDMY02 공식 기준 5개는 지구 모래 표입니다. G/K는 표의 80kPa 값에서 (101/80)^0.5로 환산해 refPress101kPa에 맞췄습니다. 표의 서로 다른 기준압 표기를 처리한 선택입니다.</li><li>화성 PDMY02 4개는 각 DP의 ρ/G/K/φ를 사용하고, refPress5kPa를 임의 기준으로 선택했습니다. 수축·팽창·임계상태·e·압력지수는 Dr50 기준의 미보정 이전값입니다. 물리적 Dr50을 뜻하지 않습니다. numerical_c300Pa는 DP의 c280Pa와 같은 보정값이 아닙니다.</li><li>Qz/Py/Tz 힘과 변위는 지구 얕은기초 문헌 예제입니다. QzSimple1은 QzSimple2와 같은 입력을 쓴 법칙 비교입니다. ElasticPPGap/Hysteretic은 문헌 E·c·φ로 규모를 정한 합성 기능 비교이며 실제 경화층·암반 절리의 곡선이 아닙니다.</li></ul>
<h2>현재 사용 판단</h2><p><strong>첫 탄성·DP 비교군</strong>은 층별 초기 강성·정적 항복 민감도 비교에 사용할 초안입니다. <strong>PDMY02</strong>는 반복 곡선·압력 의존 강성의 추가 보정이 필요합니다. 특히 화성 모래·자갈 초안은 자동 증분 분할이 필요했으므로 전체 지반에 바로 적용할 재료로 선정하지 않았습니다. <strong>단축 스프링</strong>은 링 기초의 실제 하중–침하·분담면적·들림 조건으로 보정해야 합니다. Py/Tz의 법선력 연계는 포함하지 않았습니다.</p>
<p>천부 수십 cm의 국부 층서와 3/13m 참조 경계는 서로 다른 근거입니다. 저속도 퇴적층의 파속은 Hobiger 표2의 사전모델을 사용했으며 실제 역산 곡선의 특정 깊이 물성을 읽은 것으로 간주하지 않습니다.</p>
<h2 id="cases">모든 후보와 응답 곡선</h2>''')
    for group,title in labels.items():
        body.append(f'<h3>{title}</h3><div class="cards">')
        for r in records:
            if r['group']!=group:continue
            name=r['id'];res=results[name]
            refinements=sum(t.get('step_bisections',0) for t in res['tests'])
            body.append(f'<article class="card" id="case_{name}"><h3>{name}</h3><div>{esc(layers.get(r["layer"],r["layer"]))}</div><div class="muted">{r["family"]} · <span class="ok">수치 검사 통과</span> · Mars 보정 미완료</div><a href="results/{name}/response.png"><img loading="lazy" src="results/{name}/response.png" alt="{name} 응답 곡선"></a><div class="muted">하중 시험 {len(res["tests"])}개 · 자동 증분 분할 {refinements}회</div><p><a href="cases/{name}.py">개별 실행 코드</a> · <a href="cases/{name}.json">물성·출처</a> · <a href="results/{name}/summary.json">검사 상세</a></p>')
            body.append('<details><summary>선정한 매개변수</summary><pre>'+esc(json.dumps(r['parameters'],ensure_ascii=False,indent=2))+'</pre></details></article>')
        body.append('</div>')
    body.append('<h2 id="sources">보관한 출처</h2><p><a href="../../docs/references/index.html">기존 논문 모음</a>에 분류했습니다. Kenda(2020)는 HAL 공개 PDF를 추가 확보했습니다. 아래 파일과 SHA-256은 <a href="source_manifest.json">출처 매니페스트</a>에 기록했습니다.</p><ul>')
    for src in sources:
        rel=os.path.relpath(src['file'],HERE)
        body.append(f'<li><a href="{esc(rel)}">{src["id"]} 보관 원문</a> · <a href="{esc(src["url"])}">원 출처</a></li>')
    body.append('''</ul><h2 id="run">실행 방법</h2><pre>cd /Users/sydney/Desktop/capstone/Capstone_structure/material/insight_layer_specimens_20261119
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \\
  /Users/sydney/anaconda3/envs/cap/bin/python specimens.py

# 하나만 실행
/Users/sydney/anaconda3/envs/cap/bin/python cases/duricrust_dp.py

# 물성 JSON·CSV 다시 생성 / 결과 색인만 다시 생성
/Users/sydney/anaconda3/envs/cap/bin/python build_catalog.py
/Users/sydney/anaconda3/envs/cap/bin/python make_report.py</pre>
<p>개별 실행은 해당 시편 결과와 summary_시편명.json을 갱신합니다. 전체 summary.json과 색인은 전체 실행 후 다시 생성합니다. 원본·저장 후보 자동 교체는 하지 않습니다. <a href="specimens.py">주 해석 코드</a>는 데이터 → 모델링 → 해석 → 시각화 순입니다.</p>
<p class="muted">모든 결과 CSV는 SI 단위입니다. 그래프에서만 kPa, kN, mm, %로 환산했습니다. 질량은 1m³ 고유치 시험, 탄성/DP 힘·구속압 오차는 10⁻⁵ 상대오차 이내, DP 항복은 10⁻⁴ 이내를 통과 기준으로 사용했습니다. PDMY02는 NormDispIncr 10⁻¹¹m와 최대 8단계 이분할을 사용하고 결과 힘 평형·수직응력 오차를 별도 검사했습니다. 실패한 재시도 로그는 results/run.log에 남아 있고 모든 최종 이력은 완료됐습니다.</p></html>''')
    (HERE/'index.html').write_text(''.join(body))
    (HERE/'README.txt').write_text('InSight 층별 재료 시편 — 2026-11-23\n\nindex.html을 열면 28개 시편·물성·출처·응답 곡선을 볼 수 있습니다.\n주 실행: /Users/sydney/anaconda3/envs/cap/bin/python specimens.py\n개별 실행: /Users/sydney/anaconda3/envs/cap/bin/python cases/duricrust_dp.py\n\n첫 비교군 13개(탄성7/DP6), PDMY02 9개(화성 초안4/공식 지구 기준5), 단축6개.\n1m 정육면체 / 별도 zeroLength 재료 시험. 전체 지반 해석이 아닙니다.\n수치검사 28/28 통과, 하중 이력50개, 질량검사22개. 실험·현장 보정 완료를 뜻하지 않습니다.\n매개변수별 문헌값/이전값/가정: parameters_and_sources.csv, materials.json\n문헌 원값: literature_extraction.csv\n검사 상세: results/summary.json\n기존 MGS 저장본과 전체 지반은 변경하지 않았고 JEZ-1은 제외했습니다.\n')
    manifest=[]
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or p.name=='file_manifest.json' or '__pycache__' in p.parts:continue
        manifest.append(dict(file=str(p.relative_to(HERE)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    (HERE/'file_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('Report, result CSV, raw extraction table and file hashes written.')


if __name__=='__main__':main()
