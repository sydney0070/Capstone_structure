"""Training-only fixed-E profiles and identifiability audit; no material acceptance.

Uses the existing chronological CSVs, J2/PlateFiber/ShellMITC4 coupon and bounds.
Never reads validation response files or changes their preprocessing.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import sys

import numpy as np
from scipy.optimize import least_squares

B=Path(__file__).resolve().parent
sys.path.insert(0,str(B.parent))
from experiments import simulate
from tune import read_csv,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=B/'aluminum_audit_20261012')
    args=p.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    plan=dict(date='2026-10-12',E_profiles_GPa=[60,70,80],
              purpose='Diagnose training-data support for E and saturation parameters; not accepted material values.',
              fixed='Original nu=0.33, H=0, shell thickness and strain preprocessing retained.',
              selection='No holdout reading or new final material selection.',max_nfev=150,
              initial_window='Strictly positive strain with corrected stress 10 to 50 MPa. Descriptive apparent slope only.',
              limitations='A sparse low-stress window cannot establish true E or quantify machine compliance.')
    write_json(out/'frozen_plan.json',plan)
    records={}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for row,state in enumerate(('t4','t6')):
        root=B/f'aluminum_6014_{state}'
        c=json.loads((root/'fitted/best_config.json').read_text())
        m=json.loads((root/'targets.json').read_text())
        specs=[s for s in m['tests'] if s['role']=='train']
        data=[(s,read_csv(root/s['csv'])) for s in specs]
        hashes={s['csv']:hashlib.sha256((root/s['csv']).read_bytes()).hexdigest() for s in specs}
        early=[]
        for s,d in data:
            mask=(d['x']>0)&(d['y']>=10000)&(d['y']<=50000)
            x,y=d['x'][mask],d['y'][mask]
            slope=float(np.polyfit(x,y,1)[0]/1e6) if len(x)>=2 and np.ptp(x)>0 else None
            early.append(dict(csv=s['csv'],points_10_to_50_MPa=int(mask.sum()),
                              apparent_slope_GPa=slope,total_points=len(d['x']),
                              negative_strain_points=int(np.sum(d['x']<0))))
            axes[row,0].plot(d['x']*100,d['y']/1000,'.-',label=s['csv'].removesuffix('.csv'))
        keys=['sig0','sigInf','delta'];bounds=np.asarray([c['bounds'][k] for k in keys])
        lo=bounds[:,0];width=np.diff(bounds,axis=1).ravel()
        z0=(np.array([c['params'][k] for k in keys])-lo)/width
        def residual(z,E):
            cc=copy.deepcopy(c);cc['params'].update(zip(keys,map(float,lo+z*width)));cc['params']['E']=E
            parts=[]
            for s,d in data:
                rr=simulate(cc,d['x'])
                scale=max(float(np.ptp(d['y'])),float(np.max(np.abs(d['y']))),1e-12)
                parts.append((rr['y']-d['y'])/scale/np.sqrt(len(data)*len(d['y'])))
            return np.concatenate(parts)
        profiles=[]
        for E_GPa in plan['E_profiles_GPa']:
            E=E_GPa*1e6
            result=least_squares(lambda z:residual(z,E),np.clip(z0,1e-7,1-1e-7),
                                 bounds=(np.zeros(3),np.ones(3)),diff_step=1e-4,
                                 max_nfev=plan['max_nfev'],ftol=1e-8,xtol=1e-8,gtol=1e-8)
            profiles.append(dict(E_GPa=E_GPa,train_nrmse=float(np.sqrt(2*result.cost)),
                                 parameters=dict(zip(keys,map(float,lo+result.x*width))),
                                 boundary_parameters=[k for k,v in zip(keys,result.x) if min(v,1-v)<.02],
                                 optimizer_success=bool(result.success),nfev=int(result.nfev)))
            print(state,profiles[-1],flush=True)
        # Finite-difference sensitivity of all 4 fitted parameters at the OLD fit.
        full_keys=list(c['bounds']);bb=np.asarray([c['bounds'][k] for k in full_keys])
        sensitivities=[]
        for j,k in enumerate(full_keys):
            h=1e-4*(bb[j,1]-bb[j,0]);low=copy.deepcopy(c);high=copy.deepcopy(c)
            low['params'][k]=max(bb[j,0],c['params'][k]-h)
            high['params'][k]=min(bb[j,1],c['params'][k]+h)
            cols=[]
            for s,d in data:
                scale=max(float(np.ptp(d['y'])),float(np.max(np.abs(d['y']))),1e-12)
                diff=simulate(high,d['x'])['y']-simulate(low,d['x'])['y']
                cols.append(diff/(high['params'][k]-low['params'][k])*(bb[j,1]-bb[j,0])/scale/np.sqrt(len(data)*len(d['y'])))
            sensitivities.append(np.concatenate(cols))
        jac=np.column_stack(sensitivities);sv=np.linalg.svd(jac,compute_uv=False)
        a=jac[:,full_keys.index('sigInf')];b=jac[:,full_keys.index('delta')]
        cosine=float(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))
        records[state]=dict(status='MEASUREMENT_AND_IDENTIFIABILITY_REVIEW_REQUIRED',training_only=True,
                            initial_window=early,E_profiles=profiles,
                            normalized_local_sensitivity_condition=float(sv[0]/sv[-1]) if sv[-1]>0 else None,
                            saturation_delta_sensitivity_cosine=cosine,
                            source_sha256=hashes,bounds_unchanged=c['bounds'])
        axes[row,0].set(xlabel='Engineering strain (%)',ylabel='Stress (MPa)',title=f'AW6014 {state.upper()}: training curves')
        axes[row,0].legend(fontsize=8);axes[row,0].grid(alpha=.2)
        axes[row,1].plot([v['E_GPa'] for v in profiles],[v['train_nrmse']*100 for v in profiles],'o-')
        axes[row,1].set(xlabel='Fixed E (GPa)',ylabel='Training NRMSE (%)',title='Local refit within original bounds')
        axes[row,1].grid(alpha=.2)
    fig.savefig(out/'diagnostics.png',dpi=160);plt.close(fig)
    write_json(out/'assessment.json',records)
    lines=['# 알루미늄 측정·식별성 점검 — 2026-10-12','',
           '**최종 물성 채택 보류 유지.** 기존 학습 시편만 사용해 E를 60/70/80 GPa로 고정하고 나머지 3계수를 기존 범위 안에서 국부 재최적화했다. 검증 시편은 사용하지 않았다.',
           'E를 가정했을 때 다른 계수가 보상하는 정도를 보는 진단이며, 전역 최적해나 물리적 E 측정이 아니다. sigInf 탐색 하한을 낮추거나 변형률 오프셋을 임의로 추가하지 않았다.','',
           '| 재료 | 고정 E (GPa) | 학습 NRMSE | 경계 부근 계수 | 최적화 수렴 |','|---|---:|---:|---|---|']
    for state,v in records.items():
        for r in v['E_profiles']:
            lines.append(f'| {state.upper()} | {r["E_GPa"]} | {r["train_nrmse"]:.3%} | {r["boundary_parameters"]} | {r["optimizer_success"]} |')
    lines += ['', '| 학습 시편 | 10–50 MPa 구간 점 수 | 해당 구간 겉보기 기울기 (GPa) |','|---|---:|---:|']
    for state,v in records.items():
        for r in v['initial_window']:
            slope='산정 안 함' if r['apparent_slope_GPa'] is None else f'{r["apparent_slope_GPa"]:.2f}'
            lines.append(f'| {state.upper()} {r["csv"]} | {r["points_10_to_50_MPa"]} | {slope} |')
    lines += ['', '위 기울기는 점 수가 적은 구간의 단순 회귀값이며 탄성계수 채택값이 아니다. 구간 선택·초기 정착·측정 변위 정의에 민감할 수 있다.']
    for state,v in records.items():
        lines.append(f'{state.upper()}: 기존 보정점의 정규화 국부 감도 조건수 {v["normalized_local_sensitivity_condition"]:.3g}, sigInf–delta 감도 방향 코사인 {v["saturation_delta_sensitivity_cosine"]:.6f}. 1에 가까우면 두 계수가 곡선에 유사한 영향을 주므로 분리 식별이 어렵다. 신뢰구간 또는 유일해 증명이 아니다.')
    lines += ['', '추가 확인: 원본 변형률이 어느 측정 길이·센서에서 계산되었는지, 초기 정착/시험기 컴플라이언스 보정 여부, 조밀한 탄성 구간. 포화강도에는 균일 변형 구간의 더 넓은 소성 자료와 응력·변형률 척도 일치가 필요하다.',
              '출처: [Böhnke 등 데이터 v2](https://zenodo.org/records/21993433), CC BY-SA 4.0. 원본 workbook 및 변환 이력은 각 재료 폴더의 source.xlsx·provenance.json. 본 진단 그림/파생 데이터도 CC BY-SA 4.0.',
              '재현: `python material_calibration/audit_aluminum.py --out material_calibration/새진단폴더명`.']
    (out/'RESULTS_KO.md').write_text('\n'.join(lines)+'\n',encoding='utf8')


if __name__=='__main__':main()
