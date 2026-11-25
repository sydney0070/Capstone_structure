"""Compare complete and failed ground runs without confusing timeout and collapse.

All-Gauss-point stress CSVs are independently reprocessed to check the actual
range where a compression-fitted rock-mass surrogate is used.
"""
import csv
import hashlib
import html
import json
import math
from pathlib import Path

import numpy as np
from evaluate import table
from build_materials import hb_sigma1
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RUNS=ROOT/'ground/results/nonlinear_layers_20261125'
OUT=RUNS/'comparison'


def save(path,data):path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def csvout(name,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def audit(summary,run,case):
    rows=[]
    for stage in ['gravity','depressurized','pressurized']:
        path=run/(stage+'_material_state.csv')
        if not path.exists():continue
        a=np.genfromtxt(path,delimiter=',',names=True,dtype=None,encoding='utf8')
        for material in np.unique(a['material']):
            d=a[a['material']==material];m=summary['ground_materials'][material]
            stress=np.zeros((len(d),3,3));stress[:,0,0]=d['sx_Pa'];stress[:,1,1]=d['sy_Pa'];stress[:,2,2]=d['sz_Pa'];stress[:,0,1]=stress[:,1,0]=d['sxy_Pa'];stress[:,1,2]=stress[:,2,1]=d['syz_Pa'];stress[:,0,2]=stress[:,2,0]=d['szx_Pa']
            # compression-positive principal stresses, ascending sigma3..sigma1
            principal=-np.linalg.eigvalsh(stress)[:,::-1]
            s3=principal[:,0];s1=principal[:,2]
            p=d['mean_compression_Pa'];q=d['q_Pa'];phi=math.radians(m['phi_deg']);c=m['c_Pa']
            rho=2*math.sqrt(2)*math.sin(phi)/(math.sqrt(3)*(3-math.sin(phi)));sy=6*c*math.cos(phi)/(3-math.sin(phi))
            yieldf=math.sqrt(2/3)*(q-sy)-3*rho*p
            row=dict(case=case,stage=stage,material=material,gauss_points=len(d),min_sigma3_compression_Pa=float(s3.min()),max_sigma3_compression_Pa=float(s3.max()),max_sigma1_compression_Pa=float(s1.max()),tensile_gauss_points_over_1Pa=int(sum(s3 < -1.)),max_DP_yield_residual_Pa=float(yieldf.max()),DP_positive_yield_residual_relative=float(np.maximum(yieldf,0).max()/max(float(q.max()),abs(sy),1.)),max_DP_plastic_norm=float(d['DP_plastic_strain_norm'].max()),plastic_GP_over_1e_8=int(sum(d['DP_plastic_strain_norm']>1e-8)))
            row['min_vertical_compression_Pa']=float((-d['sz_Pa']).min())
            row['max_vertical_compression_Pa']=float((-d['sz_Pa']).max())
            cal=m.get('calibration',{})
            if m.get('role')=='sand' and 'A_Pa' in cal:
                row.update(vertical_stress_above_literature_680Pa_GP=int(sum(-d['sz_Pa']>680.)),vertical_stress_below_literature_97Pa_GP=int(sum(-d['sz_Pa']<97.)),stress_scope='Vertical-normal-stress range screen only; actual3D stress path differs from the shear-box test')
            if 'constants' in cal:
                low,high=cal['fit_sigma3_range_Pa'];inside=(s3>=low-1.)&(s3<=high+1.)
                hb=hb_sigma1(np.maximum(s3,0),cal['constants'],cal['sigma_ci_Pa'])
                row.update(HB_compression_range_GP=int(inside.sum()),HB_out_of_compression_fit_range_GP=int((~inside).sum()),HB_range_pass=bool(inside.all()),max_HB_compression_utilization=float(np.max(s1/hb)),interpretation='Compression-range/model-scope check only; utilization is not a bearing-capacity safety factor')
            rows.append(row)
    return rows


def main():
    OUT.mkdir(exist_ok=True);cases=[];states=[];performance=[]
    for status in sorted(RUNS.rglob('execution_status.json')):
        d=json.loads(status.read_text());run=status.parent/'run';p=run/'summary.json'
        if not p.exists():continue
        s=json.loads(p.read_text());case=str(status.parent.relative_to(RUNS));ok=bool(d.get('complete')) and s.get('complete')
        row=dict(case=case,profile=d['profile'],radial_subdivisions=d['radial_subdivisions'],complete=bool(ok),elapsed_s=d.get('elapsed_s'),stop_reason=d.get('external_stop_reason') or d.get('limit_reason') or s.get('error'),stage=s.get('stage'),elements=s.get('ground_elements'),nodes=s.get('ground_nodes'),footing_steps=s.get('numerical_settings',{}).get('footing_steps',8),gravity_steps=s.get('numerical_settings',{}).get('static_steps',20),linear_system=s.get('numerical_settings',{}).get('linear_system','UmfPack'),downforce_N=s['loads']['downforce_N'],pressurized_net_force_N=s['loads']['pressurized_net_compression_N'],footing_area_m2=s['footing'].get('modeled_footprint_m2'),independent_field_validation=False)
        if ok:
            row.update(depressurized_max_mm=s['depressurized']['maximum_mm'],pressurized_max_mm=s['pressurized']['maximum_mm'],depressurized_average_mm=s['depressurized']['average_mm'],pressurized_average_mm=s['pressurized']['average_mm'],max_force_balance_relative=max(s[st]['force_balance']['relative_error'] for st in ['gravity','depressurized','pressurized']),force_balance_pass=all(s[st]['force_balance']['relative_error']<1e-5 for st in ['gravity','depressurized','pressurized']))
            states+=audit(s,run,case)
        else:
            cp=d.get('last_checkpoint') or s.get('solver_progress') or {};row.update(last_attempt_before=cp.get('before'),last_attempt_increment=cp.get('increment'),last_algorithm=cp.get('algorithm'),interpretation='Incomplete numerical run; no bearing-capacity or collapse conclusion')
        cases.append(row)
        journal=run/'solver_attempts.jsonl'
        if journal.exists():
            events=[json.loads(line) for line in journal.read_text().splitlines() if line.strip()]
            pending=None
            for event in events:
                if event['event']=='attempt_started':pending=event
                elif event['event']=='attempt_finished' and pending:
                    performance.append(dict(case=case,stage={800001:'gravity',800002:'depressurized',800003:'pressurized'}.get(event['pattern'],str(event['pattern'])),load_factor_before=event['before'],increment=event['increment'],algorithm=event['algorithm'],iterations=event['iterations'],returncode=event['code'],attempt_elapsed_s=event['elapsed_s']-pending['elapsed_s'],unbalance_norm_N=event.get('unbalance_norm_N'),returned=True));pending=None
            if pending:
                performance.append(dict(case=case,stage={800001:'gravity',800002:'depressurized',800003:'pressurized'}.get(pending['pattern']),load_factor_before=pending['before'],increment=pending['increment'],algorithm=pending['algorithm'],iterations=None,returncode=None,attempt_elapsed_s=None,returned=False,interpretation='Native call had no completed journal event before stop/read; exact internal function not profiled'))
    csvout('ground_cases.csv',cases);csvout('all_GP_scope_audit.csv',states)
    csvout('solver_performance.csv',performance)
    complete=[r for r in cases if r['complete']]
    # Pair only cases with matching profile; display changed solver/loading too.
    differences=[]
    for i,a in enumerate(complete):
        for b in complete[i+1:]:
            if a['profile']!=b['profile']:continue
            differences.append(dict(case_a=a['case'],case_b=b['case'],radial_a=a['radial_subdivisions'],radial_b=b['radial_subdivisions'],footing_steps_a=a['footing_steps'],footing_steps_b=b['footing_steps'],depressurized_max_change_pct=100*(b['depressurized_max_mm']/a['depressurized_max_mm']-1),pressurized_max_change_pct=100*(b['pressurized_max_mm']/a['pressurized_max_mm']-1),depressurized_average_change_pct=100*(b['depressurized_average_mm']/a['depressurized_average_mm']-1),pressurized_average_change_pct=100*(b['pressurized_average_mm']/a['pressurized_average_mm']-1),interpretation='Combined discretization/algorithm/step comparison unless recorded settings match; not automatically mesh convergence'))
    csvout('discretization_comparisons.csv',differences)
    oldpath=ROOT/'ground/results/ground_dome_footing/recommended_20261124/run/summary.json';old=json.loads(oldpath.read_text())
    oldcomparison=[]
    for r in complete:
        oldcomparison.append(dict(case=r['case'],old_depressurized_max_mm=old['depressurized']['maximum_mm'],new_depressurized_max_mm=r['depressurized_max_mm'],change_pct=100*(r['depressurized_max_mm']/old['depressurized']['maximum_mm']-1),old_pressurized_max_mm=old['pressurized']['maximum_mm'],new_pressurized_max_mm=r['pressurized_max_mm'],pressurized_change_pct=100*(r['pressurized_max_mm']/old['pressurized']['maximum_mm']-1),same_radial_mesh=r['radial_subdivisions']==old['ground_mesh']['radial_subdivisions'],scope='Layer laws and sand matching/cohesion changed jointly; coarse/fine differences cannot be attributed only to material'))
    csvout('old_profile_comparison.csv',oldcomparison)
    profilepath=ROOT/'ground/ground_modeling/ground_hybrid/configs/Mars_nonlinear_20261125_profiles.json';profiles=json.loads(profilepath.read_text())
    result=dict(date='2026-11-25',complete_cases=len(complete),total_cases=len(cases),cases=cases,scope_audit=states,discretization_comparisons=differences,old_comparison=oldcomparison,current_default_profile=profiles['selected_profile'],selection_status=profiles['selection_status'],all_materials_nonlinear=True,independent_Mars_validation=False,conditional_localization_concern=True,selection_basis='Evidence/method/scope and numerical robustness, not smallest settlement; no unique measured optimum with unknown blocky/GSI parameters')
    save(OUT/'comparison.json',result)
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for r in complete:
        if r['profile']!='nl_ss_block35':continue
        a=np.load(RUNS/r['case']/'run/stage_displacements.npz');xyz=a['coordinates'];radius=np.linalg.norm(xyz[:,:2],axis=1)
        mask=np.isclose(xyz[:,2],-3.,atol=1e-6)&(radius>=13.-1e-6)&(radius<=27.+1e-6)
        labels=np.round(radius[mask],6)
        for ax,stage in zip(axes,['depressurized','pressurized']):
            values=-(a[stage+'_m'][mask,2]-a['gravity_m'][mask,2])*1000
            x=np.unique(labels);means=[float(np.mean(values[labels==v])) for v in x]
            ax.plot(x,means,'o-',label=f'r{r["radial_subdivisions"]} / load{r["footing_steps"]} / gravity{r["gravity_steps"]}')
    for ax,stage in zip(axes,['Depressurized','Pressurized']):
        ax.set(title=stage+' | same nominal material profile',xlabel='Footing radius (m)',ylabel='Angular mean settlement (mm)');ax.grid(alpha=.2)
        if ax.lines:ax.legend(fontsize=8)
    fig.savefig(OUT/'settlement_radial_profiles.png',dpi=150);plt.close(fig)
    baseline=[r for r in complete if r['radial_subdivisions']==2 and r['gravity_steps']==100 and r['footing_steps']==32 and r['profile'] in ('nl_ss_block30','nl_ss_block35','nl_ss_block40')]
    smallest=min(baseline,key=lambda r:r['depressurized_max_mm']) if baseline else None
    style='body{max-width:1200px;margin:40px auto;padding:0 24px;font:16px/1.6 -apple-system,BlinkMacSystemFont,"Apple SD Gothic Neo",sans-serif;background:#f5f7f8;color:#20333b}section{padding:24px;margin:22px 0;background:white;border-radius:12px}.scroll{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:10px;border-bottom:1px solid #d9e3e5;text-align:left;vertical-align:top}th{background:#e7f0f1}img{width:100%;height:auto}a{color:#126d80}.note{background:#fff5e6;padding:14px;border-left:4px solid #d69036}'
    page='<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>화성 3D 지반 비선형 층별 조합 비교</title><style>'+style+'</style><h1>화성 3D 지반 비선형 층별 조합 비교</h1>'
    page+=f'<p>2026년11월25일 · 완료 {len(complete)} / 실행 {len(cases)} · 기본 프로필 {html.escape(profiles["selected_profile"])} · {html.escape(profiles["selection_status"])}</p>'
    page+='<section><h2>동일한 형상과 하중</h2><p>모든 조합은90×90×30 m, 링 중심 반경20 m·바닥 폭14 m·굴착3 m,24개 원주 분할을 사용했다. 다각형 재하 면적1739.264 m², 감압 하향하중150.3 MN, 가압 후 순압축29.58 MN이다. 층 경계는0.2·3·13 m이며 모든 층이3D DruckerPrager 탄소성이다. 재하 면적과 하중은 메시 비교에서도 유지한다. 링 밑면은3 m 깊이의 암괴층 경계에 있고 접촉응력86.4·17.0 kPa는 이 밑면에 작용한다.</p><p>침하는 굴착된 지반의 자중 평형 이후 증분이다. 가압 단계는 내압에 의한 들림만큼 하중을 줄인 상태이며 하중이0인 잔류침하 시험은 아니다. 수직 등가 재하 해석으로, 링–지반 접촉/분리·인발과 구조물 설계 판정은 포함하지 않는다.</p></section>'
    page+='<section><h2>완료한 조합</h2>'+table(['실행','요소','재하 증분','감압 최대(mm)','가압 최대(mm)','힘 평형 상대오차'],[[r['case'],r['elements'],r['footing_steps'],f'{r["depressurized_max_mm"]:.6f}',f'{r["pressurized_max_mm"]:.6f}',f'{r["max_force_balance_relative"]:.3g}'] for r in complete])+'</section>'
    page+='<section><h2>어떤 조합이 좋은가</h2><p>자료와 대응하는 시험 경로를 우선하면 경화층의 기울기 추정 DP, 모래의 MGS1 또는 MGS1S 단순전단 대응 DP, 암괴층의 마찰 범위 비교군, 현무암의 GSI55 압축 대응 DP가 기준 후보군이다. 암괴 마찰각·점착력과 암반 GSI는 현장 보정값이 없어 하나의 실제 최적 조합을 확정할 수 없다.</p>'
    if smallest:page+=f'<p>같은 기준 메시(r2)·자중100증분·재하32증분으로 완료한 비교 조합 중 감압 침하가 가장 작았던 것은 <strong>{html.escape(smallest["profile"])}</strong>의{smallest["depressurized_max_mm"]:.4f} mm다. 이 값은 해당 가정에서의 응답 순위이며, 강도를 높여 얻은 작은 침하를 재료의 타당성 순위로 사용하지 않았다.</p>'
    selection=RUNS/'selection.json'
    if selection.exists():
        selected=json.loads(selection.read_text());page+='<p class="note">'+html.escape(selected.get('recommendation_ko',''))+'</p>'
    page+='<p class="note">추가 접선 안정성 진단에서 기준 지반의 항복면 근처 모래 일부 점에 국부화 가능성이 확인되었다. 이는 소성 재하 접선을 가정한 조건부 진단으로 전역 붕괴 판정은 아니다. 구성식·평형 검증 통과와 최종 재료 채택을 구분해야 한다. <a href="../../../../material/nonlinear_layers_20261125/evaluation/report.html">접선 검증과 음향 텐서 진단</a></p>'
    page+='<p>기존 하부 두 층이 탄성인 조합은 감압1.433760 mm·가압0.282703 mm였다. <a href="old_profile_comparison.csv">전체 이전값 비교</a>에는 메시가 같은지도 기록했다. 다른 메시 결과의 차이를 재료만의 효과로 해석하지 않는다.</p></section>'
    page+='<section><h2>해의 신뢰성과 재하 간격</h2><img src="settlement_radial_profiles.png">'+table(['비교 A','비교 B','감압 최대 변화(%)','가압 최대 변화(%)'],[[r['case_a'],r['case_b'],f'{r["depressurized_max_change_pct"]:.4f}',f'{r["pressurized_max_change_pct"]:.4f}'] for r in differences])+'<p>연구용 목표는 메시·증분 변경에 따른 침하 차이5% 이하다. 위 그래프는 같은 반경의 각도 평균이며 최대 절점 침하와 다르다. 알고리즘이나 자중 증분도 달라졌다면 결합된 수치 민감도이며 메시 수렴만으로 부르지 않는다. 비선형 경로 의존성과 마지막 작은 증분의 수렴을 별도로 확인한다.</p></section>'
    page+='<section><h2>모든 적분점의 재료 사용 범위</h2>'+table(['실행·단계','재료','소성 적분점/전체','최소 주압축(kPa)','최대 주압축(kPa)','HB 범위 밖'],[[r['case']+' / '+r['stage'],r['material'],str(r['plastic_GP_over_1e_8'])+'/'+str(r['gauss_points']),f'{r["min_sigma3_compression_Pa"]/1000:.3f}',f'{r["max_sigma1_compression_Pa"]/1000:.3f}',r.get('HB_out_of_compression_fit_range_GP','—')] for r in states if r['stage']=='depressurized'])+'<p>DP 상태변수의 소성변형률 노름1e-8 초과를 소성 발생 진단으로 사용했다. 이는 손상률이나 파괴 판정이 아니다. 현무암은 인장 한계를 과대평가하므로 최소 주압축이 양수이고 보정된0–500 kPa의 σ3 범위에 있는지 검사했다. HB 이용률도 지반 전체의 지지력 안전율로 해석하지 않는다.</p></section>'
    sandstates=[r for r in states if r['stage']=='depressurized' and 'vertical_stress_above_literature_680Pa_GP' in r]
    page+='<section><h2>모래 시험의 응력 범위를 넘는가</h2>'+table(['실행','모래 최대 수직압축(kPa)','680 Pa 초과 적분점/전체'],[[r['case'],f'{r["max_vertical_compression_Pa"]/1000:.4f}',str(r['vertical_stress_above_literature_680Pa_GP'])+'/'+str(r['gauss_points'])] for r in sandstates])+'<p>모래 포락선의 시험 수직응력은97–680 Pa다. 기준 조합의 모래층 수직압축은 최대 약10 kPa로 이 범위를 초과한다. 기초 밑면의86.4 kPa를 모래층 전체가 받는다고 비교하지 않고 실제 모래 적분점의 수직응력을 따로 검사했다. 이 검사는3D 응력 경로가 직접전단 시험과 같다는 증명은 아니다.</p></section>'
    failed=[r for r in cases if not r['complete']]
    page+='<section><h2>계산 시간과 미완료 실행</h2>'+table(['실행','단계','마지막 재하 이전 계수','중지 이유'],[[r['case'],r['stage'],r.get('last_attempt_before'),r['stop_reason'] or '실행 중'] for r in failed])+'<p>수렴하지 않은 실행의 중간 변위를 최종 침하 순위에 넣지 않았다. native attempt timeout은 Python으로 반환되지 않는 OpenSees 호출에 대한 자원 한계이며 물리적인 붕괴 하중이 아니다. 일부 실행의 내부 함수는 추가로1초 샘플링했다. 자중200증분·r3 실행의 정체 시점은 UmfpackGenLinSolver::solve와 umfpack_di_numeric의 행렬 분해에 집중돼 있었다. 프로세스 메모리는약2.0GB, 최대2.5GB였고 CPU사용률은98.5%였다. 이 관측은 해당 시점 한 건에 한정되며 모든 실패 원인이나 재료 적분의 무한루프를 증명하지 않는다. <a href="../performance/profiling_observation.json">프로파일 관측 기록</a> 일부 실행은80회 반복 후 큰 불평형과 함께 실패했고, 작은 자중·재하 증분의 별도 실행은 완료되었다.</p><p>현재 기본 메시의3,504개 요소·4,093개 절점에서100회 자중 증분과 각32회 이상의 재하·감압을 수행한다. 탄소성 활성점의 변화 때문에 반복 평형과 희소 비대칭 행렬의 풀이를 거듭한다. <a href="solver_performance.csv">시도별 반복 횟수와 실제 소요시간</a>을 저장했다. 재료 강도·강성·허용 잔차를 바꾸지 않고 수치 방법을 비교했다.</p></section>'
    page+='<section><h2>자료·재현·남은 검증</h2><p>독립 원시 강도·변형 곡선과 InSight 암괴층의 c·φ, 암반 GSI가 필요하다. 모래층의 실제 수직응력도97–680 Pa 시험 범위를 넘어 외삽이 남는다. 또한 기초 밑면3 m가 가정한 모래층 끝3 m와 일치하므로 층 경계의 불확실성도 별도 검증해야 한다. 단일요소와 FE 평형 검증을 통과해도 실제 지지력·인발 안전이 검증되었다고 결론내리지 않는다.</p><p><a href="../../../../material/nonlinear_layers_20261125/evaluation/report.html">시편·논문값 평가 보고서</a> · <a href="../../../../docs/references/index.html">분류한 논문 라이브러리</a> · <a href="ground_cases.csv">실행 결과 CSV</a> · <a href="all_GP_scope_audit.csv">모든 적분점 범위 감사 CSV</a> · <a href="comparison.json">전체 판정 JSON</a></p></section></html>'
    (OUT/'report.html').write_text(page)
    print(json.dumps(dict(cases=len(cases),complete=len(complete),state_groups=len(states),comparisons=len(differences)),ensure_ascii=False))


if __name__=='__main__':main()
