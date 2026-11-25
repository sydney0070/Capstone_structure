"""Independent specimen validation and local sensitivity of the fitted coupon."""
from pathlib import Path
import sys,json
import numpy as np
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE.parent))
from experiments import simulate
from tune import read_csv,write_json

if __name__=='__main__':
    root=BASE/'s355_hem320_web'
    c=json.loads((root/'fitted/best_config.json').read_text())
    m=json.loads((root/'targets.json').read_text());plan=json.loads((root/'acceptance_plan.json').read_text())
    fit=json.loads((root/'fitted/fit_report.json').read_text());records=[]
    for spec in m['tests']:
        d=read_csv(root/spec['csv']);r=simulate(c,d['x'])
        scale=max(np.ptp(d['y']),max(abs(d['y'])),1e-12)
        err=float(np.sqrt(np.mean((r['y']-d['y'])**2))/scale)
        peak=float(max(abs(max(r['y'])-max(d['y'])),abs(min(r['y'])-min(d['y'])))/max(abs(d['y'])))
        dense=json.loads(json.dumps(c));dense['test']['substeps']=2
        rd=simulate(dense,d['x']);inc=float(max(abs(rd['y']-r['y']))/max(abs(rd['y'])))
        records.append(dict(csv=spec['csv'],role=spec['role'],range_nrmse=err,
            peak_stress_relative_error=peak,increment_refinement_relative_error=inc,
            passes_screening=bool(err<=plan['validation_range_nrmse_max'] and peak<=plan['validation_peak_stress_relative_error_max'])))
    boundary=[];cols=[]
    training=[read_csv(root/s['csv']) for s in m['tests'] if s['role']=='train']
    for key,(lo,hi) in c['bounds'].items():
        v=c['params'][key];ratio=(v-lo)/(hi-lo)
        if min(ratio,1-ratio)<.02:boundary.append(key)
        step=(hi-lo)*1e-4
        ca=json.loads(json.dumps(c));cb=json.loads(json.dumps(c))
        ca['params'][key]=min(hi,v+step);cb['params'][key]=max(lo,v-step)
        diff=[]
        for d in training:
            ra=simulate(ca,d['x']);rb=simulate(cb,d['x'])
            scale=max(np.ptp(d['y']),max(abs(d['y'])),1e-12)
            diff.extend((ra['y']-rb['y'])/(ca['params'][key]-cb['params'][key])*(hi-lo)/scale/np.sqrt(len(d['x'])))
        cols.append(diff)
    singular=np.linalg.svd(np.array(cols).T,compute_uv=False)
    passed=all(r['passes_screening'] for r in records if r['role']=='validation')
    status='VALIDATED_WITHIN_SELECTED_COUPON_SCOPE' if passed and fit['optimizer_success'] else 'FITTED_BUT_NOT_ACCEPTED'
    assessment=dict(status=status,optimizer_success=fit['optimizer_success'],tests=records,
        boundary_parameters=boundary,normalized_sensitivity_singular_values=singular.tolist(),
        sensitivity_condition_number=float(singular[0]/max(singular[-1],1e-15)),
        warning='Local sensitivity is not parameter uniqueness or statistical confidence; selected scope is not general grade certification.',
        applicability='Same S355 J2+M HEM320 web material source, selected room-temperature pre-fracture history prefixes only.')
    write_json(root/'assessment.json',assessment)
    c['status']=status;c['calibration_scope']=assessment['applicability']
    c['source_doi']='10.5281/zenodo.6965147'
    write_json(root/'reviewed_config.json',c)
    lines=['# S355 구조용 강재 시편 보정 결과','',f'상태: **{status}**','',
        '대상은 철근이 아닌 S355 J2+M HEM320 웹에서 채취한 구조용 강재다. 공개 실험 시편 2개로 학습하고 별도 시편 2개로 검증했다. Steel02 단축 모델이며, 다축 금속 쉘이나 파괴·피로까지 보정한 결과는 아니다.','',
        '| 매개변수 | 보정값 |','|---|---:|']
    for key in c['bounds']:lines.append(f'| {key} | {c["params"][key]:.8g} |')
    lines+=['','E·Fy 단위는 kPa이고 나머지는 무차원이다. 고정값도 `reviewed_config.json`에 함께 저장했다.','',
        '| 자료 | 역할 | 범위 정규화 RMSE | 피크 응력 오차 | 검토 기준 |','|---|---|---:|---:|---|']
    for r in records:lines.append(f'| {r["csv"]} | {r["role"]} | {r["range_nrmse"]*100:.3f}% | {r["peak_stress_relative_error"]*100:.3f}% | {"통과" if r["passes_screening"] else "미달"} |')
    lines+=['','사전 검토 기준은 검증 시편 각각 범위 정규화 RMSE ≤5%, 피크 응력 오차 ≤10%이다. 설계 규격이 아닌 이번 시편 보정의 작업 기준이다.',
        f'최적화 종료 성공: {fit["optimizer_success"]}. 경계 근처(탐색 폭의 2% 이내) 변수: {boundary or "없음"}.',
        f'국부 정규화 감도 행렬의 조건수: {assessment["sensitivity_condition_number"]:.3g}. 이를 유일해 또는 신뢰구간으로 해석하지 않는다.','',
        '공개 저자 처리자료의 순서를 유지했다. 진응력·로그변형률을 균일 변형 가정으로 공칭 값으로 변환하고 초기의 작은 오프셋을 제거했다. 원자료에서 처음으로 절대 공칭변형률 2%를 초과하는 순간 이후는 사용하지 않았다. 이는 파단까지 포함한 전체 실험 재현이 아니다. 추가 다운샘플링이나 변형률 정렬은 하지 않았다.','',
        '[데이터베이스·CC BY 4.0](https://doi.org/10.5281/zenodo.6965147), [해당 실험의 원 연구](https://doi.org/10.5075/epfl-thesis-8700). 원본 시편 파일, 처리 이력, 메타데이터와 해시는 `provenance.json` 및 `*_source.csv`에 보관했다.','',
        '세부 곡선은 `fitted/train_1`, `train_2`, `validation_1`, `validation_2`의 `response.png`에서 확인한다. 후보 재료를 채택하기 전 실제 사용할 재료 등급·제품·온도·제조 조건과의 일치를 확인해야 한다.']
    (root/'RESULTS_KO.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(assessment,indent=2))
