"""Evaluate archived FE responses; distinguish verification from data fitting.

No material parameter is fitted by this evaluator. Experimental curve RMSE is
left unavailable because independent raw Mars stress-strain curves were not supplied.
"""
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/private/tmp/capstone_specimen_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from build_materials import hb_sigma1

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'evaluation'


def save(name,data):
    (OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def write_csv(name,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def table(head,rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def main():
    OUT.mkdir(exist_ok=True)
    records=json.loads((HERE/'materials.json').read_text())['materials'];byid={r['id']:r for r in records}
    run=json.loads((HERE/'results/summary.json').read_text());byrun={r['id']:r for r in run['cases']}
    stability_path=HERE/'stability/stability_assessment.json'
    stability=json.loads(stability_path.read_text()) if stability_path.exists() else None
    criteria=[
        dict(id='N_E',metric='Initial and unloading Young modulus',formula='abs(E_FE/E_input-1)',relative_limit=1e-5,kind='numerical verification',basis='Project tolerance; load/unload path method R30/R29'),
        dict(id='N_G',metric='Initial shear modulus',formula='abs(G_FE/G_input-1)',relative_limit=1e-4,kind='numerical verification',basis='Project tolerance; elastic probe before yielding'),
        dict(id='N_Y',metric='DP compression/tension yield',formula='abs(q_FE/q_DP_analytic-1)',relative_limit=1e-4,kind='numerical verification',basis='Implemented DP surface, independently calculated scalar analytic formula'),
        dict(id='N_SS',metric='Fully developed DP simple-shear resistance',formula='abs(tau_FE/tau_DP_analytic-1)',relative_limit=1e-4,kind='numerical verification',basis='Project tolerance; equal normal stresses and zero dilation at steady shear'),
        dict(id='N_M',metric='1m3 consistent material mass',formula='abs((k/eigenvalue)/(rho*V)-1)',relative_limit=1e-7,kind='numerical verification',basis='Rigid translation eigenproblem; not a measured density test'),
        dict(id='N_F',metric='Force balance / confinement / GP spread',formula='max(abs(residual))/reference',relative_limit=1e-5,kind='numerical verification',basis='Project tolerance; FE reactions and controlled pressure'),
        dict(id='C_SS',metric='Literature MC shear envelope',formula='abs(tau_FE/(c_target+p*tan(phi_target))-1)',relative_limit=.05,kind='calibration screening',basis='5% is project screening target, not a tolerance prescribed by the papers'),
        dict(id='C_HB',metric='Compression-fitted DP against HB equation',formula='abs(sigma1_DP/sigma1_HB-1)',relative_limit=.10,kind='analogue surrogate screening',basis='10% is project target within sigma3=0–500kPa; not independent data validation'),
        dict(id='G_MESH',metric='Settlement between mesh/step comparisons',formula='abs(s_new/s_reference-1)',relative_limit=.05,kind='solution verification',basis='Project screening target; failure to meet means unresolved discretization sensitivity'),
        dict(id='S_Q',metric='Conditional plastic-loading stability screen',formula='min_sampled_normal det(Q)/G^3',relative_limit=None,threshold='Positive sampled determinants; a negative value flags possible localization, not global collapse',kind='conditional stability diagnostic',basis='R31 acoustic-tensor rationale; DP tangent derived and checked against native stress differences; limited directions/states'),
    ]
    for criterion in criteria:criterion['limit_pct']=100*criterion['relative_limit'] if criterion['relative_limit'] is not None else None
    write_csv('criteria.csv',criteria)
    rows=[];paths=[];cross=[];uncertainty=[];worst={}
    metrickeys=['yield_error_pct','initial_E_error_pct','unloading_E_error_pct','force_balance_relative','lateral_error_relative','gauss_point_spread_relative','steady_shear_formula_difference_pct']
    for r in records:
        name=r['id'];p=r['parameters'];cal=r['calibration'];result=byrun[name];tests=result['tests']
        ss=[t for t in tests if t['test']=='monotonic_simple_shear']
        errors=np.array([t['target_MC_difference_pct'] for t in ss]);diff=np.array([t['final_tau_Pa']-t['target_MC_tau_Pa'] for t in ss])
        shear_mape=float(np.mean(abs(errors))) if len(ss) else None
        shear_rmse=float(np.sqrt(np.mean(diff**2))) if len(ss) else None
        hb=[]
        if 'constants' in cal:
            # Off-grid pressures excluded from the101-point fit. This is an
            # equation interpolation check, not new measured experimental data.
            pressure=np.linspace(0,500000,1002)[1:-1]
            target=hb_sigma1(pressure,cal['constants'],cal['sigma_ci_Pa'])
            phi=math.radians(p['phi_deg']);s=math.sin(phi)
            fitted=pressure+(2*pressure*s+2*p['c_Pa']*math.cos(phi))/(1-s)
            hb=100*(fitted/target-1)
            write_csv(name+'_offgrid_HB.csv',[dict(sigma3_Pa=float(x),HB_sigma1_Pa=float(y),DP_sigma1_Pa=float(z),error_pct=float(e),kind='off-grid analytical surrogate check, not independent experiment') for x,y,z,e in zip(pressure,target,fitted,hb)])
        tc=[t for t in tests if t['test']=='triaxial_compression_unload']
        for t in tc:
            ph=math.radians(cal['target_MC_phi_deg']);s=math.sin(ph)
            targetq=(2*t['pressure_Pa']*s+2*cal['target_MC_c_Pa']*math.cos(ph))/(1-s)
            cross.append(dict(specimen=name,axis=t['axis'],pressure_Pa=t['pressure_Pa'],FE_q_Pa=t['computed_yield_q_Pa'],target_MC_q_Pa=targetq,cross_path_difference_pct=100*(t['computed_yield_q_Pa']/targetq-1),interpretation='SS-fitted DP does not also fit physical-MC triaxial compression; model-form tradeoff' if cal['target_mode']=='simple_shear' else 'TC mapping target, not independent experiment'))
        for t in tests:
            paths.append(dict(specimen=name,**{k:v for k,v in t.items() if not isinstance(v,dict)}))
            for k in metrickeys:
                value=t.get(k)
                if value is not None and (k not in worst or abs(value)>worst[k]['value']):worst[k]=dict(value=abs(value),specimen=name,path=t['test'])
        rank='reference/sensitivity'
        decision='구성식·평형 구현 확인; 현장 물성 검증은 남아 있음'
        if name=='crust_tc':rank='duricrust: baseline';decision='경사 해석 추정 강도·천부 유효 강성 사용; 밀도1200→1025 이전 가정'
        elif name=='mgs1_ss':rank='sand: basaltic-composition baseline';decision='직접전단 대응 개선; 일반 현무암질 조성을 가정한 기준 비교용'
        elif name=='mgs1s_ss':rank='sand: sulfate-rich alternative';decision='직접전단 대응 개선; 황산염 풍부 조성의 대안, 현장 조성 확인 필요'
        elif name.endswith('_ss_nearest'):rank='sand: density-transfer sensitivity';decision='가까운 시험 밀도1310/1320의 강도를1300에 옮긴 민감도 비교'
        elif name in ['mgs1_tc','mgs1s_tc']:rank='sand: TC comparison';decision='삼축압축 대응 비교군; 문헌 직접전단 포락선을 과대평가'
        elif name.startswith('blocky_'):rank='blocky: uncertainty family';decision='타지역·지구 유사재 범위 시나리오; InSight c·φ 실측값 없음'
        elif name=='basalt_hb_gsi55':rank='basalt: compression baseline';decision='평균 화산암반 GSI55 가정의 압축 비교용; 인장강도 과대평가'
        elif name.startswith('basalt_hb'):rank='basalt: GSI sensitivity';decision='GSI 불확실성 비교용; 더 큰 강도·더 작은 보정오차로 적합성을 확정할 수 없음'
        elif name=='basalt_lab_reference':rank='exclude from ground';decision='지구의 온전한 현무암 강도 참고용; 파쇄 지반에는 직접 채택하지 않음'
        tension=next((t.get('DP_tension_to_HB_difference_pct') for t in tests if 'DP_tension_to_HB_difference_pct' in t),None)
        rows.append(dict(specimen=name,layer=r['layer'],law=r['family'],numerical_pass=result['passed'],loading_paths=len(tests),mass_error_relative=result['mass']['relative_error'],shear_target_MAPE_pct=shear_mape,shear_target_RMSE_Pa=shear_rmse,HB_offgrid_max_error_pct=float(max(abs(hb))) if len(hb) else None,DP_tension_HB_difference_pct=tension,group_recommendation=rank,decision=decision,independent_Mars_validation=False,independent_experimental_curve_RMSE=None))
        if stability:
            screen=next(x for x in stability['candidate_path_screen'] if x['specimen']==name)
            rows[-1].update(conditional_pure_shear_min_det_Q_over_G3=screen['pure_shear_min_det_Q_over_G3'],conditional_stability_scope=screen['scope'],field_adoption='Research comparison only; calibration, stability and global solution verification remain separate')
        if 'A_Pa' in cal:
            A=cal['A_Pa'];k=cal['k'];dA=.34 if name.startswith('mgs1_') else 1.14;dk=.22
            uncertainty.append(dict(specimen=name,nominal_density_g_cm3=1.3,nominal_c_Pa=A*math.exp(k*1.3),lower_corner_c_Pa=(A-dA)*math.exp((k-dk)*1.3),upper_corner_c_Pa=(A+dA)*math.exp((k+dk)*1.3),scope='Deterministic coefficient-corner sensitivity; covariance absent, not a confidence interval; density extrapolation retained'))
    write_csv('material_assessment.csv',rows);write_csv('loading_path_metrics.csv',paths);write_csv('cross_path_predictions.csv',cross);write_csv('density_fit_uncertainty.csv',uncertainty)
    summary=dict(date='2026-11-25',materials=len(rows),numerically_passed=sum(r['numerical_pass'] for r in rows),loading_paths=len(paths),mass_tests=len(rows),independent_Mars_validated=0,worst_numerical_metrics=worst,materials_assessed=rows,interpretation='Numerical passes, fitted literature-envelope agreement, analogue/HB agreement and independent field validation are different claims. No independent field validation is available.')
    summary['max_mass_error_relative']=max(r['mass_error_relative'] for r in rows)
    if stability:summary['conditional_pure_shear_negative_determinant_count']=sum(x['pure_shear_min_det_Q_over_G3']<0 for x in stability['candidate_path_screen'])
    save('evaluation.json',summary)
    # Export shareable scientific figures without browser dependencies.
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':145})
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for ax,prefix in zip(axes,['mgs1','mgs1s']):
        r=byid[prefix+'_ss'];c=r['calibration'];x=np.linspace(0,700,101)
        ax.plot(x,c['target_MC_c_Pa']+x*np.tan(np.deg2rad(c['target_MC_phi_deg'])),color='black',label='Literature fit (density extrapolated)')
        for suffix,color in [('tc','#c76b29'),('ss','#167d87')]:
            t=[t for t in byrun[prefix+'_'+suffix]['tests'] if t['test']=='monotonic_simple_shear']
            ax.plot([t['pressure_Pa'] for t in t],[t['final_tau_Pa'] for t in t],'o--',label=suffix.upper()+'-matched DP',color=color)
        ax.set(title=prefix.upper()+' | density1300 kg/m3',xlabel='Normal stress (Pa)',ylabel='Steady shear stress (Pa)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(OUT/'sand_envelopes.png');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,3.8),constrained_layout=True)
    for ax,name in zip(axes,['mgs1_ss','blocky_phi35','basalt_hb_gsi55']):
        a=np.genfromtxt(HERE/'results'/name/'triaxial_axis2_p50000.csv',delimiter=',',names=True)
        ax.plot(a['axial_compression_strain']*100,a['q_Pa']/1000,color='#167d87');ax.set(title=name,xlabel='Axial compression strain (%)',ylabel='Deviator stress q (kPa)');ax.grid(alpha=.2)
    fig.savefig(OUT/'loading_unloading.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for name in ['basalt_hb_gsi40','basalt_hb_gsi55','basalt_hb_gsi70']:
        r=byid[name];cal=r['calibration'];p=r['parameters'];x=np.linspace(0,500000,301);hb=hb_sigma1(x,cal['constants'],cal['sigma_ci_Pa']);s=np.sin(np.deg2rad(p['phi_deg']));dp=x+(2*x*s+2*p['c_Pa']*np.cos(np.deg2rad(p['phi_deg'])))/(1-s)
        line=axes[0].plot(x/1000,hb/1e6,label='HB GSI'+str(cal['GSI']))[0];axes[0].plot(x/1000,dp/1e6,'--',color=line.get_color());axes[1].plot(x/1000,100*(dp/hb-1),label='GSI'+str(cal['GSI']))
    axes[0].set(xlabel='Minor principal compression (kPa)',ylabel='Major principal compression (MPa)',title='HB analogue (solid) / DP fit (dashed)');axes[1].set(xlabel='Minor principal compression (kPa)',ylabel='DP / HB difference (%)',title='Compression surrogate error');axes[1].axhline(10,color='grey',ls=':');axes[1].axhline(-10,color='grey',ls=':')
    for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(OUT/'rock_compression_fit.png');plt.close(fig)
    sources=json.loads((ROOT/'docs/references/_catalog/sources.json').read_text())
    ids=['R13','G01','G07','G13','G14','G15','G16','G17','G18','G19','G21','R15','R28','R29','R30','R31'];refs=[s for s in sources if s['id'] in ids]
    refhtml='<ul>'
    for s in refs:
        files=s.get('files',[]);links=' · '.join('<a href="../../../docs/references/'+html.escape(f)+'">'+html.escape(Path(f).name)+'</a>' for f in files if str(f).lower().endswith('.pdf'))
        refhtml+='<li><strong>'+html.escape(s['id']+' '+s['title'])+'</strong> '+links+' · <a href="'+html.escape(s['url'])+'">원 출처</a><br>'+html.escape(s.get('used_section',''))+'</li>'
    refhtml+='</ul>'
    report='''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>화성 층별 비선형 재료 시편 평가</title><style>body{max-width:1150px;margin:40px auto;padding:0 24px;font:16px/1.65 -apple-system,BlinkMacSystemFont,"Apple SD Gothic Neo",sans-serif;color:#20333b;background:#f5f7f8}h1,h2{line-height:1.3}section{background:white;padding:24px;margin:22px 0;border-radius:12px}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #dae1e4;padding:10px;text-align:left;vertical-align:top}th{background:#e7f0f1}img{width:100%;height:auto}a{color:#126d80}code{font-size:13px}p.note{border-left:4px solid #d69036;padding:12px;background:#fff5e6}.kpi{display:flex;gap:20px;flex-wrap:wrap}.kpi b{font-size:30px;color:#167d87}li{margin-bottom:14px}</style><h1>화성 층별 비선형 재료 시편 평가</h1><p>2026년11월25일 · SI 단위 · 3D nD DruckerPrager · 독립 현장 검증 전 연구용</p>'''
    report+=f'<section><div class="kpi"><span><b>{summary["numerically_passed"]}/20</b> 구성식·평형 통과</span><span><b>{len(paths)}</b> 재하 이력</span><span><b>20</b> 질량 검사</span><span><b>0</b> 독립 화성 현장 검증</span></div><p>모든 새 시편은 탄소성 재료다. 경화·연화와 팽창은 측정 자료가 없어0으로 두었다. 탄소성 재료라도 사용 하중에서 탄성 구간에 머무를 수 있다.</p><p class="note">이 보고서의 작은 수치 오차는 OpenSees가 지정한 재료식을 재현하는지에 대한 결과다. 같은 문헌값으로 보정한 강도 포락선을 다시 맞춘 것은 독립 실험 검증이 아니다. 독립 원시 응력–변형률 자료가 없어 실험 곡선 RMSE·동적 감쇠·체적변형 오차는 산출하지 않았다. 별도 접선 안정성 진단에서는 국부화 가능성이 나와 최종 지반 재료로 확정하지 않았다.</p></section>'
    report+='<section><h2>시편과 재하 방법</h2><p>1×1×1 m의 bbarBrick8절점 요소를 사용했다. 각 시편에서 재료 밀도의 질량을 고유치로 확인하고, 측압을 먼저 적용한 다음 변위 제어로 삼축압축·제하를 수행했다. 같은 압력을 X·Y·Z축으로 반복해 등방성을 확인했다. 모래·경화층·암괴층은 일정 수직응력 단순전단과 반전 재하를 추가했고 현무암은 일축 인장도 수행했다. 모래 전단의 수직응력은97·250·500·680 Pa다.</p><p>문헌의 직접전단 시험을 단일 요소의 균질 단순전단으로 근사했다. 전단상자 경계와 국부 전단대, 입자 파쇄는 재현하지 않는다. PEER2012/05 §4.5의 여러 응력 경로를 이용한 수치 검증을 참고했지만, 해당 보고서의 정수압 초기화·고정 체적 단순전단 경계조건을 그대로 복제한 시험은 아니다. NASA-STD-7009B는 구현 검증·실험 타당성·불확실성을 구분하는 근거이며 오차 허용치의 출처는 아니다.</p></section>'
    report+='<section><h2>평가 기준</h2>'+table(['항목','비교식','기준','판정 범위'],[[c['metric'],c['formula'],str(c['limit_pct'])+'%' if c['limit_pct'] is not None else c['threshold'],c['kind']] for c in criteria])+'<p>오차 허용치는 이 프로젝트의 연구용 판정값이다. 문헌의 공인 합격기준으로 표현하지 않는다. 음향 텐서는 상대오차 허용치가 아닌 별도 안정성 진단이다. 보정 기준을 넘은 비교군도 삭제하지 않고 오차를 남겼다.</p></section>'
    report+='<section><h2>구성식 검증에서 실제 나온 오차</h2>'+table(['지표','최대 오차','단위','발생 시편'],[[k,f'{v["value"]:.6g}','%' if k.endswith('_pct') else '상대값',v['specimen']] for k,v in worst.items()]+[['mass',f'{summary["max_mass_error_relative"]:.6g}','상대값','전체20종']])+'<p>위 값은 입력한 구성식·경계조건·평형에 대한 오차다. 재료의 실제 강도·강성 추정 오차나 실험 재현 오차가 아니다.</p></section>'
    chosen=[byid[name] for name in ['crust_tc','mgs1_ss','blocky_phi35','basalt_hb_gsi55']]
    report+='<section><h2>기준 조합의 실제 입력값</h2>'+table(['시편','E(MPa)','ν','ρ(kg/m³)','DP c(kPa)','DP φ(°)','목표 MC c(kPa)·φ(°)'],[[r['id'],f'{r["parameters"]["E_Pa"]/1e6:g}',r['parameters']['nu'],r['parameters']['density_kg_m3'],f'{r["parameters"]["c_Pa"]/1000:.6g}',f'{r["parameters"]["phi_deg"]:.6g}',f'{r["calibration"]["target_MC_c_Pa"]/1000:.6g} · {r["calibration"]["target_MC_phi_deg"]:.6g}'] for r in chosen])+'<p>공통 ψ=0°, 경화=0이다. 현무암의 큰 유효 마찰각은 낮은 구속압 범위의 HB 곡선을 선형으로 근사한 DP 입력이며, 현장에서 직접 측정한 마찰각이 아니다.</p><p>실행한 OpenSees3.4.0의 항복식은 ||s||+ρ<sub>DP</sub>I₁−√(2/3)σY=0이다. 삼축 대응 변환에서 σY=6c cosφ/(3−sinφ), ρ<sub>DP</sub>=2√2 sinφ/[√3(3−sinφ)]를 사용했고 별도 삼축 해석식과 일치하는지 검사했다. 온라인 문서의 σY 변환식에는 추가 √2 분모가 있어 현재 구현·이 해석해와 불일치한다. <a href="../../../docs/references/05_official_model_docs/DruckerPrager_v3.4.0.cpp">실행 버전 소스</a>와 수치 시편 결과를 함께 보존했다.</p></section>'
    report+='<section><h2>모래의 직접전단 포락선과 비교</h2><img src="sand_envelopes.png"><p>Dotson2024 Eq2/Table4의 밀도–점착력 관계로1300 kg/m³의 점착력을 산정했다. MGS1은235.664 Pa, MGS1S는238.156 Pa다. 각각 시험 최소밀도1320·1310 kg/m³보다 낮으므로 외삽이다. 마찰각35.5°·34°는 인접 밀도 값이다. 기존의 삼축압축 대응 DP는 이 직접전단 포락선을 약16.57%·17.66% 크게 계산했다. 새 단순전단 대응 변환은 이 오차를 줄인다.</p><p>단순전단 대응 DP 입력은 물리적인 MC c·φ 자체와 다르다. 두 값을 별도 열에 보존했다. 직접전단 대응을 개선하면 삼축압축 대응은 달라진다. <a href="cross_path_predictions.csv">다른 응력 경로 예측 오차</a>를 함께 확인해야 한다. 같은 데이터의 보정오차가 작다는 이유로 MGS1을 MGS1S보다 우수한 현장 재료로 단정하지 않는다.</p>'+table(['시편','포락선 MAPE(%)','포락선 RMSE(Pa)','권고'],[[r['specimen'],f'{r["shear_target_MAPE_pct"]:.8g}',f'{r["shear_target_RMSE_Pa"]:.8g}',r['group_recommendation']] for r in rows if r['layer']=='sand'])+'</section>'
    report+='<section><h2>비선형 거동과 현무암 압축 근사</h2><img src="loading_unloading.png"><p>압축에서 항복 이후 영구변형과 탄성 제하를 확인했다. 반복 전단의 전단 일은 경로에 따른 적분값이며, 실험 감쇠비나 전체 소산에너지로 해석하지 않는다.</p><img src="rock_compression_fit.png"><p>Villeneuve·Heap의 지구 현무암 강도·암반 스케일 관계와 Hoek–Brown 식을 이용해0–500 kPa의 최소 주압축응력 구간을 DP로 근사했다. GSI40·55·70의 보정 오차는 각각 약7.27%·2.58%·0.60%다. GSI55는 평균 암반 가정이고40·70은 불확실성 비교다. GSI70이 강하거나 오차가 작다고 실제 화성 암반으로 더 적합한 것은 아니다.</p><p class="note">압축에 맞춘 DP의 인장강도는 HB 절단값보다 GSI40·55·70에서 각각 약441%·301%·225% 크다. 따라서 인장·박리·균열 검증용으로 채택하지 않는다. 지반 적용에서는 암반 적분점의 주응력이 압축 상태이고 보정 구간 안에 있는지 별도로 확인한다. Villeneuve 식8의 부호는 Hoek–Brown2019의 정식 식6과 대조해 후자의 음수 부호를 적용했다.</p></section>'
    report+='<section><h2>층별 후보의 사용 판단</h2>'+table(['시편','층','수치검증','강도 대조 오차(%)','권고·한계'],[[r['specimen'],r['layer'],'PASS' if r['numerical_pass'] else 'FAIL',f'{r["HB_offgrid_max_error_pct"]:.4g}' if r['HB_offgrid_max_error_pct'] is not None else f'{r["shear_target_MAPE_pct"]:.4g}' if r['shear_target_MAPE_pct'] is not None else '—',r['decision']] for r in rows])+'<p>3–13 m 암괴층은 Hughes2022의 다른 화성 지역·지구 유사재 범위인 마찰각20–40°를 이용했다. 무점착 파쇄암 시나리오와 점착력1–5 kPa의 혼합 입상체 시나리오는 서로 다른 지반 가정이다. 수렴을 위해 실제 강도를 임의 조정한 값으로 숨기지 않는다.</p></section>'
    if stability:
        native_error=max(x['native_versus_continuum_tangent_max_error_over_G'] for x in stability['native_checks'])
        selected_states=[x for x in stability['ground_snapshot_screen'] if x['stage']=='depressurized']
        report+='<section><h2>구성식이 맞아도 전체 지반이 불안정할 수 있는가</h2><p class="note">팽창각0°·경화0의 비연관 DP는 일부 소성 응력 경로에서 국부화와 메시 의존성이 나타날 수 있다. 따라서 위의20/20 PASS는 재료식과 평형 구현에 대한 결과이며 전체 지반의 수학적 안정성까지 통과했다는 뜻이 아니다.</p><img src="../stability/acoustic_determinant.png"><p>순수전단의 소성 재하 접선으로 음향 텐서 Q를 계산했다. 그림의0 미만 영역은 국부화 가능성이 있는 방향을 뜻한다. 각도는 주응력 좌표의 법선 방향이며 지반 경사각이 아니다. 팽창각을 마찰각과 같게 둔 선은 수학적 대조군으로, 지반 재료를 이 값으로 바꾸지 않았다.</p>'
        report+=f'<p>접선 계산식은 실제 OpenSees 응력의6방향 유한차분과 대조했고 최대 차이는 G 대비{native_error:.6g}였다. 별도의 직접 텐서 연산으로 행렬식 계산도 대조했다. 기준 지반 감압 상태에서 항복면 근처 모래4,476개 점 중736개는 조사한 법선에서 음수 행렬식이 나왔다. 이는 해당 응력에서 앞으로 소성 재하가 일어난다고 가정한 조건부 진단이다.</p>'
        report+=table(['재료','항복면 근처 점','음수 행렬식 점','최소 det(Q)/G³'],[[x['material'],x['near_yield_GP'],x['negative_sampled_det_GP'],f'{x["min_det_Q_over_G3"]:.6g}' if x['min_det_Q_over_G3'] is not None else '비해당'] for x in selected_states])
        report+='<p>세 주응력 평면을0.5° 간격으로 조사했다. 축 끝점 중복을 포함한543개 표본이며 모든 방향·모든 중간 상태의 안정성 증명은 아니다. 항복면 근처라는 사실만으로 마지막 증분이 소성 재하였다고 확정할 수도 없다. 이 진단은 전역 지지력이나 native 호출 정체의 원인을 확정하지 않는다.</p><p>Sabet·de Borst2019의 안정성·정규화 논의를 참고했다. 논문의 Mohr–Coulomb 임계 경화식을 DP에 그대로 대입하지 않고3D DP 접선을 따로 유도했다. 정규화 모델·체적변형·경화 자료를 보정하기 전까지 후속 소성 변형의 최종 예측값으로 채택하는 데는 한계가 있다. <a href="../stability/stability_assessment.json">진단 전체값</a> · <a href="../stability/candidate_path_stability.csv">20종 경로별 진단</a> · <a href="../check_stability.py">재현 코드</a></p></section>'
    strain_path=HERE/'strain_audit/strain_ranges.json'
    if strain_path.exists():
        strain_data=json.loads(strain_path.read_text());strain_rows=[r for r in strain_data['stage_material_strain_ranges'] if r['stage']=='depressurized']
        report+='<section><h2>실제 지반에서 나타난 변형률 범위</h2>'+table(['층 재료','최대 전단변형률(%)','자중 이후 증분 최대 전단변형률(%)'],[[r['material'],f'{r["max_engineering_shear_strain_pct"]:.6g}',f'{r["max_postgravity_increment_shear_strain_pct"]:.6g}'] for r in strain_rows])+'<p>기준3,504개 요소·자중100증분 실행의 저장 변위에서 bbarBrick 소변형 운동학으로 복원했다. 전단변형률은 주변형률의 최댓값과 최솟값의 차이다. 원소 한 개의 비정형 형상에서 OpenSees 변형률과 비교하고, 지반 전체의 체적변형률–평균응력 관계도 대조했다. 전체 이력의 극값이 아닌 세 단계의 저장 상태다.</p><p>모래의 최대 전단변형률은 약0.189%, 자중 이후 증분은 약0.0717%다. 파속으로 환산한 초기 강성이 이 변형률 범위의 실제 강성까지 검증한다는 뜻은 아니다. 실험 G/Gmax·체적변형 곡선과 비교할 자료가 남아 있다. <a href="../strain_audit/strain_ranges.csv">층별 변형률 CSV</a> · <a href="../strain_audit/strain_ranges.json">운동학 검증과 전체 결과</a></p></section>'
    report+='<section><h2>MGS1과 MGS1S를 고르는 추가 근거</h2><p>두 후보의 작은 보정 오차만으로 실제 지반 적합성을 순위 매길 수는 없다. Long-Fox·Britt2023 Table4에서 MGS1은 현무암질 기본 조성, MGS1S는 이 기본 조성을 바꾸고 석고를40 wt% 포함한 황산염 풍부 모사토다. 따라서 일반적인 현무암질 지반이라는 가정에서는 MGS1을 기준으로, 황산염이 풍부하다는 별도 가정에서는 MGS1S를 대안으로 둔다. 이는 조성 가정에 따른 추론이며, InSight 지반이 이 석고 함량이라는 측정 근거는 아니다. 해당 논문은 광물·원소·입도 자료로, 밀도1300의 독립 응력–변형률 곡선 검증 자료로 사용하지 않았다.</p><p>현무암의 공극률4%는 Table2의 지구 유사암 연결 공극률이다. Table3의 일축압축강도 경험식 자체 산포가 크다. 전 공극률 범위의 추정 표준오차는59.5 MPa, 공극률0.1 미만 구간에서는91.8 MPa다. DP–HB 근사 오차2.58%가 이 입력 강도 불확실성을 해소하는 것은 아니다.</p></section>'
    report+='<section><h2>불확실성과 지반 적용 범위</h2><p>모래 강도 자료의 수직응력은97–680 Pa다. 기준 조합의3,504개 요소 해석에서는 모래층 수직압축이 최대 약10.0 kPa여서 시험 범위를 넘는다. 링 밑면 접촉응력86.4·17.0 kPa는3 m 깊이 암괴층에 작용하므로 이를 모래층 응력과 같다고 비교하지 않았다. 저응력 강도식의 외삽이 남는다. 경화층의 점착력은 밀도1200·φ30°를 가정한 기울기 해석 추정을1025 kg/m³에 옮겼다. 모래·현무암의 초기 강성, 심부 층서의 공간 평균, 암반 GSI와 공극률도 각각 대용·추정 가정이 남아 있다. 기초 밑면3 m가 가정한 모래층 끝3 m와 같아서, 층 경계가 달라져 밑면 아래에 모래가 남는 경우의 민감도는 별도 확인해야 한다.</p><p><a href="density_fit_uncertainty.csv">점착력 계수 민감도</a>는 계수±범위의 모서리 조합이며 신뢰구간이 아니다. 독립 압축·전단·체적변형 자료로 보정하기 전까지 지지력·인발 안전이나 최적 설계의 확정 근거로 쓰지 않는다. <a href="../../../ground/results/nonlinear_layers_20261125/comparison/report.html">3D 지반 조합 비교 보고서</a>에서 수렴과 메시·증분 의존성, 각 층의 소성 상태를 확인한다.</p><p>이전 PDMY02 시편은 보존했다. 기존 보고서에서 저응력 포락선 차이와 반복·수축·팽창 매개변수의 추가 보정 필요를 확인했다. 이번 정적 비교는 자료로 설명 가능한 최소 매개변수의 DP를 사용하며, 이 결정이 모든 동적 문제에서 DP가 PDMY02보다 우수하다는 의미는 아니다. <a href="../../insight_layer_specimens_20261119/evaluation_20261124/report.html">이전28종 평가</a></p></section>'
    report+='<section><h2>문헌과 재현 파일</h2>'+refhtml+'<p><a href="criteria.csv">평가 기준 CSV</a> · <a href="loading_path_metrics.csv">167개 재하 결과 CSV</a> · <a href="material_assessment.csv">20종 평가 CSV</a> · <a href="../parameters.csv">재료 입력값 CSV</a> · <a href="evaluation.json">기계 판독 결과</a> · <a href="../cleanup/completed.json">지반 폴더 정리 결과</a> · <a href="../verification_corrections.json">초기 검사식·측정 구간 점검 기록</a></p><p>초기 G 검사에서는 항복 이전의 작은 변형 구간을 따로 측정하도록 수정했고, DP 인장 해석식과 MC/HB 인장 대조를 분리했다. 재료 매개변수와 허용오차는 유지했다.</p><p><code>python material/nonlinear_layers_20261125/run_specimens.py</code><br><code>python material/nonlinear_layers_20261125/check_stability.py</code><br><code>python material/nonlinear_layers_20261125/strain_audit.py</code><br><code>python material/nonlinear_layers_20261125/evaluate.py</code></p></section></html>'
    (OUT/'report.html').write_text(report)
    manifest=[]
    for p in [HERE/'materials.json',HERE/'run_specimens.py',HERE/'build_materials.py',HERE/'check_stability.py',HERE/'strain_audit.py',Path(__file__),HERE/'results/summary.json',ROOT/'docs/references/_catalog/sources.json']:
        manifest.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    save('file_manifest.json',manifest)
    print(json.dumps({k:summary[k] for k in ['materials','numerically_passed','loading_paths','mass_tests','independent_Mars_validated']},ensure_ascii=False))


if __name__=='__main__':main()
