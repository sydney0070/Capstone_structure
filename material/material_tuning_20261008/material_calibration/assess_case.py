"""Assess an experimental case without using holdout responses in optimization."""
from pathlib import Path
import sys,json,argparse
import numpy as np
B=Path(__file__).resolve().parent;sys.path.insert(0,str(B.parent))
from experiments import simulate
from tune import read_csv,write_json
p=argparse.ArgumentParser();p.add_argument('case');args=p.parse_args();root=B/args.case
c=json.loads((root/'fitted/best_config.json').read_text());m=json.loads((root/'targets.json').read_text());plan=json.loads((root/'acceptance_plan.json').read_text());fit=json.loads((root/'fitted/fit_report.json').read_text())
train=[s['csv'] for s in m['tests'] if s['role']=='train'];val=[s['csv'] for s in m['tests'] if s['role']=='validation'];assert train and val and not set(train)&set(val)
records=[]
for spec in m['tests']:
 d=read_csv(root/spec['csv']);r=simulate(c,d['x']);scale=max(np.ptp(d['y']),max(abs(d['y'])),1e-12)
 nr=float(np.sqrt(np.mean((r['y']-d['y'])**2))/scale)
 peak=float(max(abs(max(r['y'])-max(d['y'])),abs(min(r['y'])-min(d['y'])))/max(abs(d['y'])))
 cc=json.loads(json.dumps(c));cc['test']['substeps']*=2;rr=simulate(cc,d['x'])
 inc=float(max(abs(rr['y']-r['y']))/max(abs(rr['y'])))
 records.append(dict(csv=spec['csv'],role=spec['role'],nrmse=nr,peak_error=peak,increment_error=inc,passes=bool(nr<=plan['nrmse_max'] and peak<=plan['peak_error_max'] and inc<=plan['increment_error_max'])))
boundary=[k for k,(lo,hi) in c['bounds'].items() if min((c['params'][k]-lo)/(hi-lo),(hi-c['params'][k])/(hi-lo))<.02]
passed=fit['optimizer_success'] and all(r['passes'] for r in records if r['role']=='validation')
status='VALIDATED_WITHIN_SELECTED_COUPON_SCOPE' if passed else 'FITTED_BUT_NOT_ACCEPTED'
if passed and (plan.get('mandatory_review') or boundary):status='FITTED_REQUIRES_PARAMETER_OR_MEASUREMENT_REVIEW'
a=dict(status=status,numerical_screening_passed=passed,tests=records,boundary_parameters=boundary,mandatory_review=plan.get('mandatory_review'),scope=plan['scope'],optimizer_success=fit['optimizer_success'])
write_json(root/'assessment.json',a);c['status']=status;c['calibration_scope']=plan['scope'];write_json(root/'reviewed_config.json',c)
lines=[f'# {args.case} 보정 결과','',f'상태: **{status}**','',f'Material: `{c["material"]}`. Element: `'+('truss' if c['kind']=='uniaxial' else 'ShellMITC4 + PlateFiber')+'`.','',plan['scope'],'','| 변수 | 보정값 |','|---|---:|']
for k in c['bounds']:lines.append(f'| {k} | {c["params"][k]:.9g} |')
lines+=['','응력·탄성계수 단위 kPa. 그 외 단위는 starting_config.json을 따른다.','', '| 시편 | 역할 | NRMSE | 극값 오차 | 증분 오차 | 기준 통과 |','|---|---|---:|---:|---:|---|']
for r in records:lines.append(f'| {r["csv"]} | {r["role"]} | {100*r["nrmse"]:.3f}% | {100*r["peak_error"]:.3f}% | {100*r["increment_error"]:.5f}% | {r["passes"]} |')
lines+=['',f'경계 부근 매개변수: {boundary or "없음"}.',f'추가 확인사항: {plan.get("mandatory_review") or "선택 시편 범위 밖 거동은 검증하지 않음"}.','', '사전 기준: 검증 시편 각각 NRMSE ≤5%, 극값 오차 ≤10%, 내부 증분 세분화 오차 ≤1%. 설계 규격이 아닌 프로젝트 검토 기준이다. 평균 곡선 오차 통과만으로 파괴·피로·다축·환경 거동의 검증을 주장하지 않는다.','', '실험 출처·처리: targets.json, provenance.json. 곡선: fitted/train_*, fitted/validation_*의 response.png. 초기값·고정값·탐색범위: starting_config.json.']
(root/'RESULTS_KO.md').write_text('\n'.join(lines)+'\n',encoding='utf8');print(json.dumps(a,indent=2))
