"""Compare the unchanged-geometry reference with the completed specimen-profile run."""
import argparse
import csv
import hashlib
import html
import json
import os
from pathlib import Path


ROOT=Path(__file__).resolve().parents[2]


def dump_csv(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,default=ROOT/'ground/results/ground_dome_footing/recommended_20261124/run')
    parser.add_argument('--baseline',type=Path,default=ROOT/'ground/results/ground_dome_footing/analysis_mgs1_20261118/run/summary.json')
    args=parser.parse_args();run=args.run_dir.resolve();out=run.parent
    a=json.loads(args.baseline.read_text());b=json.loads((run/'summary.json').read_text())
    if not a['complete'] or not b['complete']:raise ValueError('Both runs must be complete')
    for key in ('ground_domain_m','ground_nodes','ground_elements','ground_mesh'):
        if a[key]!=b[key]:raise ValueError(f'Unmatched geometry/mesh {key}')
    for key in ('downforce_N','dome_uplift_N'):
        if a['loads'][key]!=b['loads'][key]:raise ValueError('Unmatched footing loads')
    if a['footing']['modeled_footprint_m2']!=b['footing']['modeled_footprint_m2']:raise ValueError('Unmatched footprint area')
    params=[]
    for old_layer,new_layer in zip(a['ground_layers'],b['ground_layers']):
        if old_layer['bottom_m']!=new_layer['bottom_m']:raise ValueError('Unmatched layer boundaries')
        old=a['ground_materials'][old_layer['material']];new=b['ground_materials'][new_layer['material']]
        for key in ('law','E_Pa','nu','G_Pa','K_Pa','density_kg_m3','phi_deg','c_Pa','psi_deg'):
            if key not in old and key not in new:continue
            x,y=old.get(key),new.get(key)
            numeric=isinstance(x,(int,float)) and isinstance(y,(int,float))
            params.append(dict(layer_bottom_m=old_layer['bottom_m'],old_material=old_layer['material'],new_material=new_layer['material'],
                               parameter=key,old_value=x,new_value=y,changed=x!=y,
                               difference=y-x if numeric else None,difference_pct=100*(y/x-1) if numeric and x else None))
    settlements=[]
    for stage in ('depressurized','pressurized'):
        for key in ('average_mm','maximum_mm'):
            x,y=a[stage][key],b[stage][key]
            settlements.append(dict(stage=stage,metric=key,baseline_mm=x,new_mm=y,difference_mm=y-x,difference_pct=100*(y/x-1)))
    balance=[];states=[]
    for stage in ('gravity','depressurized','pressurized'):
        s=b[stage];balance.append(dict(stage=stage,load_factor=s['final_factor'],attempts=s['attempts'],force_balance_error=s['physical_support_force_error'],force_balance_error_pct=100*s['physical_support_force_error'],force_balance_residual_z_N=s['force_balance']['residual_N'][2]))
        for name,state in s['material_state']['by_material'].items():
            states.append(dict(stage=stage,material=name,elements=state['elements'],gauss_points=state['gauss_points'],
                               max_mean_compression_kPa=state['max_mean_compression_Pa']/1000,max_q_kPa=state['max_q_Pa']/1000,
                               DP_state_recorded=state['dp_state_recorded'],plastic_gp_count=state.get('gauss_points_with_plastic_norm_above_threshold'),
                               max_plastic_strain_norm=state.get('max_plastic_strain_norm'),
                               max_abs_norm_change_since_gravity=state.get('max_abs_norm_change_since_gravity'),
                               changed_gp_count_since_gravity=state.get('gauss_points_with_changed_norm_since_gravity')))
    dump_csv(out/'material_changes.csv',params);dump_csv(out/'settlement_comparison.csv',settlements);dump_csv(out/'force_balance.csv',balance);dump_csv(out/'material_state_summary.csv',states)
    comparison=dict(complete=True,baseline_summary=str(args.baseline.resolve()),new_summary=str((run/'summary.json').resolve()),
                    geometry_and_loads_matched=True,changed_parameters=[p for p in params if p['changed']],settlement_comparison=settlements,
                    mass_change_kg=b['ground_mass_kg']-a['ground_mass_kg'],mass_change_pct=100*(b['ground_mass_kg']/a['ground_mass_kg']-1),
                    baseline_elapsed_s=a['elapsed_s'],new_elapsed_s=b['elapsed_s'],force_balance=balance,material_states=states,
                    caveats=['Only duricrust density differs physically from the completed MGS1 baseline.',
                             'All source specimen values retained. MGS1 shear calibration has not been performed.',
                             'Same 0.2/3/13m reference layers; not a continuous observed borehole.',
                             'Displacements are increments after gravity on the excavated geometry.',
                             'Positive net compression remains in pressurized stage; it is not zero-load residual settlement.',
                             'DP strain-norm change is not a plastic-strain tensor increment or safety factor.',
                             'Equivalent distributed pressure, bonded continuum; not contact/uplift/ultimate-capacity validation.'])
    (out/'comparison.json').write_text(json.dumps(comparison,ensure_ascii=False,indent=2)+'\n')
    plot_results(out,settlements,run/'stage_displacements.npz')
    render_report(out,comparison,b)
    files=[args.baseline,run/'summary.json',ROOT/'material/insight_layer_specimens_20261119/materials.json',Path(__file__).resolve()]
    (out/'comparison_sources.json').write_text(json.dumps([dict(file=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files],indent=2)+'\n')
    print(json.dumps(dict(changes=comparison['changed_parameters'],settlement=settlements,balance=balance,elapsed_s=b['elapsed_s']),ensure_ascii=False,indent=2))


def plot_results(out,settlements,stages):
    import numpy as np
    os.environ.setdefault('MPLCONFIGDIR','/private/tmp/capstone_specimen_mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=[r for r in settlements if r['metric']=='maximum_mm']
    fig,ax=plt.subplots(figsize=(7,4.5))
    x=np.arange(2)
    ax.bar(x-.18,[r['baseline_mm'] for r in rows],.36,label='Previous MGS-1 run')
    ax.bar(x+.18,[r['new_mm'] for r in rows],.36,label='Specimen profile')
    ax.set(xticks=x,xticklabels=['Depressurized','After pressurization'],ylabel='Maximum incremental settlement (mm)')
    ax.grid(axis='y',alpha=.2);ax.legend();fig.tight_layout();fig.savefig(out/'settlement_comparison.png',dpi=150);plt.close(fig)
    with np.load(stages) as d:
        coords=d['coordinates'];reference=d['gravity_m'];values=[-(d[k]-reference)[:,2]*1000 for k in ('depressurized_m','pressurized_m')]
    assert np.isfinite(coords).all() and all(np.isfinite(v).all() for v in values)
    select=np.isclose(coords[:,1],0.,atol=1e-7)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
    for ax,v,title in zip(axes,values,['Depressurized','After pressurization']):
        points=ax.scatter(coords[select,0],coords[select,2],c=v[select],cmap='viridis',vmin=0,vmax=max(values[0]),s=17)
        ax.set(title=title,xlabel='x (m)',ylim=(-31,1),xlim=(-46,46));ax.grid(alpha=.15)
    axes[0].set_ylabel('z (m)');fig.subplots_adjust(left=.07,right=.87,bottom=.19,wspace=.1)
    fig.colorbar(points,cax=fig.add_axes([.9,.23,.015,.6]),label='Downward increment after gravity (mm)')
    fig.text(.5,.04,'Nodal values on y=0 section; excavated gravity reference; no displacement magnification.',ha='center',fontsize=8)
    fig.savefig(out/'settlement_section.png',dpi=150);plt.close(fig)


def render_report(out,result,summary):
    esc=lambda x:html.escape(str(x))
    def table(rows,fields):
        text='<div class="scroll"><table><tr>'+''.join('<th>'+label+'</th>' for key,label in fields)+'</tr>'
        for row in rows:
            text+='<tr>'
            for key,label in fields:
                v=row.get(key);value='해당 없음' if v is None else f'{v:.9g}' if isinstance(v,float) else esc(v)
                text+='<td>'+value+'</td>'
            text+='</tr>'
        return text+'</table></div>'
    parts=['''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>추천 시편 적용 지반 결과</title><style>body{max-width:1280px;margin:38px auto;padding:0 26px;font:16px/1.7 -apple-system,BlinkMacSystemFont,"Apple SD Gothic Neo",sans-serif;color:#223043}h1{font-size:29px}h2{margin-top:36px}a{color:#15659a}.note{background:#fff2dc;padding:18px;border-left:4px solid #ac7021}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}td,th{text-align:left;vertical-align:top;border-bottom:1px solid #d9e1e8;padding:10px}th{background:#edf3f7}pre{white-space:pre-wrap;padding:16px;background:#f1f5f7}.muted{color:#607080;font-size:13px}</style>
<h1>추천 시편을 적용한 3D 지반 해석</h1><p>2026-11-24 · 90×90×30m · 절점4093개 / bbarBrick3504개 · 3m 매입 링 · 바닥판 폭14m</p>
<p class="note"><strong>완료된 기존 MGS-1 해석도 같은 재료 법칙을 사용했습니다.</strong> 새 시편 레코드를 직접 연결했고 실제 변경 물성은 경화층 밀도1200→1025kg/m³입니다. 모래층의 MGS-1 DP 강도는 보정 전 값이며 하부 조립질층·암반층은 탄성입니다. 큰 재료모델 변경으로 설명하지 않습니다.</p>
<p><a href="run/footing_settlement.html">감압 상태 3D 침하 보기</a> · <a href="run/footing_settlement_pressurized.html">내압 적용 후 3D 침하 보기</a> · <a href="model/ground_geometry.html">새 재료 배치 보기</a></p>
<h2>기존 결과와 침하 비교</h2><p>감압은 내압 없이150.3MN을 재하한 상태입니다. 이후 내압 등가 들림120.72MN을 적용하면 순압축력29.58MN이 남습니다. 두 값 모두 각 모델의 자중 평형 이후 증분 침하입니다. 내압 적용 후 값은 무하중 잔류침하가 아닙니다.</p>''']
    parts.append('<p><img style="max-width:780px;width:100%" src="settlement_comparison.png" alt="기존·추천 시편 지반 최대침하 비교"></p>')
    parts.append(table(result['settlement_comparison'],[('stage','상태'),('metric','항목'),('baseline_mm','기존 mm'),('new_mm','추천 시편 mm'),('difference_mm','차이 mm'),('difference_pct','변화 %')]))
    parts.append('<p><img style="max-width:1100px;width:100%" src="settlement_section.png" alt="자중 평형 이후 y0 단면 수직변위"></p>')
    parts.append('<h2>변경한 물성</h2>'+table(result['changed_parameters'],[('parameter','변경 항목'),('old_value','기존'),('new_value','변경'),('difference_pct','변화 %')]))
    parts.append(f'<p>전체 지반 질량 변화 {result["mass_change_kg"]:,.3f}kg ({result["mass_change_pct"]:+.6f}%). 경화층의 강성·강도는 같고 자중이 달라진 비교입니다. 재료 이름만 새 시편으로 바뀐 항목은 material_changes.csv의 changed=false로 구분했습니다.</p>')
    parts.append('<h2>완료 여부와 힘 평형</h2>'+table(result['force_balance'],[('stage','단계'),('load_factor','최종 하중계수'),('attempts','시도 횟수'),('force_balance_error_pct','힘 평형 상대오차 %'),('force_balance_residual_z_N','수직 잔차 N')]))
    parts.append(f'<p>해석·결과 저장 시간 {result["new_elapsed_s"]:.1f}초. 단계별하중계수1.0과 힘 평형을 확인했습니다. 기존 결과의 실행시간과 차이를 재료 자체의 계산비용으로 단정하지 않습니다.</p>')
    parts.append('<h2>층별 응력과 소성 변형 진단</h2><p>DP 상태벡터의 소성변형률 텐서 노름을 모든 적분점에서 기록했습니다. 10⁻⁸ 초과를 진단용으로 표시했으며 손상·붕괴 판정기준은 아닙니다. 자중 단계의 소성과 추가 재하에 따른 노름 변화를 구분했습니다. 탄성층에서 DP 상태가 비어 있는 것은 재료 법칙에 따른 것입니다.</p>')
    parts.append(table(result['material_states'],[('stage','단계'),('material','재료'),('max_mean_compression_kPa','최대 평균압축 kPa'),('max_q_kPa','최대 q kPa'),('plastic_gp_count','소성 노름>10⁻⁸ 적분점'),('max_plastic_strain_norm','최대 소성 노름'),('changed_gp_count_since_gravity','자중 후 노름 변화 적분점'),('max_abs_norm_change_since_gravity','자중 후 최대 노름 변화')]))
    parts.append('''<h2>결과의 적용 범위</h2><p>현재 해석은 굴착된 지반에 등가 수직 분포하중을 가하는 침하 비교입니다. 기초–지반 접촉·들림·미끄럼, 구조물의 강성에 따른 반력 재분배와 극한지지력은 검증하지 않았습니다. 이전 시편에서 확인한 문헌 강도포락선 차이도 남아 있습니다.</p>
<p><a href="settlement_comparison.csv">침하 비교 CSV</a> · <a href="material_changes.csv">전체 물성 대조 CSV</a> · <a href="material_state_summary.csv">층별 상태 CSV</a> · <a href="force_balance.csv">평형 CSV</a> · <a href="comparison.json">종합 JSON</a> · <a href="run/summary.json">해석 원본 JSON</a></p>
<h2>재실행</h2><pre>''')
    parts.append(esc('cd '+str(ROOT)+'\nenv MPLCONFIGDIR=/private/tmp/capstone_specimen_mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \\\n  /Users/sydney/anaconda3/envs/cap/bin/python ground/ground_modeling/ground_dome_footing.py \\\n  --material-profile recommended_20261124 --html --output-dir /private/tmp/ground_recommended_fresh\n\n# 기존 재료 설정은 --material-profile reference로 선택 가능\n# 예전 MGS-1 비교 설정: --material-profile reference --sandy-regolith-candidate mgs1'))
    parts.append('</pre><p class="muted">기존 완료 해석·저장된 MGS-1/MGS-1S 재료 파일은 변경하지 않았습니다. 현재 주 실행 코드의 기본 프로파일은 recommended_20261124입니다.</p></html>')
    (out/'comparison.html').write_text(''.join(parts))


if __name__=='__main__':main()
