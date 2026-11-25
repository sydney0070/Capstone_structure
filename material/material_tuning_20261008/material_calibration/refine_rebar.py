"""Training-only Steel02 refinement; preserve the original fit and acceptance plan.

Run from the package root with the existing requirements environment:
  python material_calibration/refine_rebar.py
Validation histories are opened only after the training-selected parameters are saved.
The holdouts were examined in an earlier revision; this is a reused-holdout check,
not a fresh blind validation. No source histories, bounds or acceptance limits change.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import sys
import time

import numpy as np
from scipy.optimize import least_squares

B = Path(__file__).resolve().parent
sys.path.insert(0, str(B.parent))
from experiments import simulate
from tune import read_csv, save, write_json
import openseespy.opensees as ops


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metrics(y, target):
    scale = max(float(np.ptp(target)), float(np.max(np.abs(target))), 1e-12)
    peak_scale = max(float(np.max(np.abs(target))), 1e-12)
    return dict(nrmse=float(np.sqrt(np.mean((y-target)**2))/scale),
                peak_error=float(max(abs(np.max(y)-np.max(target)),
                                     abs(np.min(y)-np.min(target)))/peak_scale))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=B/'rebar_ucsd_a6/refinement_20261012')
    args = parser.parse_args()
    out = args.out.resolve()
    # Refuse to overwrite a previous experiment or adapt its frozen plan.
    out.mkdir(parents=True, exist_ok=False)
    root = B/'rebar_ucsd_a6'
    base = json.loads((root/'fitted/best_config.json').read_text())
    manifest = json.loads((root/'targets.json').read_text())
    acceptance = json.loads((root/'acceptance_plan.json').read_text())
    train_specs = [s for s in manifest['tests'] if s['role']=='train']
    val_specs = [s for s in manifest['tests'] if s['role']=='validation']
    assert train_specs and val_specs
    assert not {s['csv'] for s in train_specs} & {s['csv'] for s in val_specs}
    files = ['fitted/best_config.json', 'starting_config.json', 'targets.json',
             'acceptance_plan.json', 'provenance.json'] + [s['csv'] for s in manifest['tests']]
    hashes = {f: digest(root/f) for f in files}
    keys = list(base['bounds'])
    bounds = np.asarray([base['bounds'][k] for k in keys], float)
    lo, width = bounds[:, 0], np.diff(bounds, axis=1).ravel()
    z0 = (np.array([base['params'][k] for k in keys])-lo)/width
    plan = dict(date='2026-10-12', material='Steel02', case=manifest['case'],
                reason='Training specimen A65 also exceeded the existing 10% extremum criterion.',
                objective='mean over train specimens of 0.8*curve_NRMSE^2 + 0.1*positive_peak_error^2 + 0.1*negative_peak_error^2',
                curve_scale='max(stress range, max absolute stress)',
                peak_scale='max absolute experimental stress',
                seed=20261012, starts=5, jitter_normalized_sd=.06, max_nfev=250,
                parameter_keys=keys, bounds=base['bounds'], acceptance=acceptance,
                train=[s['csv'] for s in train_specs], validation=[s['csv'] for s in val_specs],
                holdout_policy='Previously observed holdouts, not fresh blind validation. No holdout-based selection or retry.',
                source_sha256=hashes, units=base['units'], opensees=ops.version())
    write_json(out/'frozen_plan.json', plan)
    write_json(out/'baseline_config.json', base)
    write_json(out/'acceptance_plan.json', acceptance)
    data = [(s, read_csv(root/s['csv'])) for s in train_specs]

    def configuration(z):
        c = copy.deepcopy(base)
        c['params'].update(zip(keys, map(float, lo+z*width)))
        return c

    def residual(z):
        c = configuration(z)
        pieces = []
        for spec, d in data:
            cc = copy.deepcopy(c)
            cc['test'].update(spec.get('test', {}))
            y = simulate(cc, d['x'])['y']
            scale = max(float(np.ptp(d['y'])), float(np.max(np.abs(d['y']))), 1e-12)
            peak_scale = max(float(np.max(np.abs(d['y']))), 1e-12)
            pieces.append(np.sqrt(.8/(len(data)*len(y)))*(y-d['y'])/scale)
            pieces.append(np.sqrt(.1/len(data))*np.array([
                np.max(y)-np.max(d['y']), np.min(y)-np.min(d['y'])])/peak_scale)
        return np.concatenate(pieces)

    rng = np.random.default_rng(plan['seed'])
    starts = [z0] + [np.clip(z0+rng.normal(0,.06,len(z0)),1e-6,1-1e-6) for _ in range(4)]
    baseline_objective = float(np.dot(residual(z0), residual(z0)))
    candidates = []
    begin = time.monotonic()
    for i, z in enumerate(starts, 1):
        fit = least_squares(residual, z, bounds=(np.zeros(len(z)),np.ones(len(z))),
                            max_nfev=plan['max_nfev'], diff_step=1e-4,
                            ftol=1e-8, xtol=1e-8, gtol=1e-8)
        candidates.append(dict(start=i, z=fit.x.tolist(), objective=float(2*fit.cost),
                               success=bool(fit.success), message=str(fit.message),
                               nfev=int(fit.nfev), njev=int(fit.njev)))
        write_json(out/'training_candidates.json', candidates)
        print(json.dumps(candidates[-1]), flush=True)
    successful = [r for r in candidates if r['success']]
    if not successful:
        raise RuntimeError('No optimizer convergence; holdouts have not been assessed.')
    selected = min(successful, key=lambda r:r['objective'])
    c = configuration(np.asarray(selected['z']))
    c['status'] = 'TRAINING_SELECTED_PENDING_REUSED_HOLDOUT_CHECK'
    write_json(out/'training_selected_config.json', c)
    selected_hash = digest(out/'training_selected_config.json')

    # Numerical sensitivity uses training data only and is not a confidence interval.
    z = np.asarray(selected['z'])
    jac = []
    for j in range(len(z)):
        a, b = z.copy(), z.copy()
        a[j] = max(0., z[j]-1e-4); b[j] = min(1., z[j]+1e-4)
        jac.append((residual(b)-residual(a))/(b[j]-a[j]))
    sv = np.linalg.svd(np.column_stack(jac), compute_uv=False)
    condition = float(sv[0]/sv[-1]) if sv[-1]>0 else None
    records = []
    for spec in train_specs+val_specs:
        d = read_csv(root/spec['csv'])
        cc = copy.deepcopy(c); cc['test'].update(spec.get('test', {}))
        before = copy.deepcopy(base); before['test'].update(spec.get('test', {}))
        r = simulate(cc, d['x'])
        old = metrics(simulate(before, d['x'])['y'], d['y'])
        value = metrics(r['y'], d['y'])
        fine = copy.deepcopy(cc); fine['test']['substeps'] *= 2
        finer = copy.deepcopy(cc); finer['test']['substeps'] *= 4
        r2, r4 = simulate(fine, d['x']), simulate(finer, d['x'])
        err2 = float(np.max(np.abs(r2['y']-r['y']))/max(np.max(np.abs(r2['y'])),1e-12))
        err4 = float(np.max(np.abs(r4['y']-r2['y']))/max(np.max(np.abs(r4['y'])),1e-12))
        balance = float(np.max(np.abs(r['reaction']-r['y']))/max(np.max(np.abs(r['y'])),1e-12))
        passes = (value['nrmse']<=acceptance['nrmse_max'] and
                  value['peak_error']<=acceptance['peak_error_max'] and
                  max(err2,err4)<=acceptance['increment_error_max'] and balance<=1e-8)
        rec = dict(csv=spec['csv'],role=spec['role'],**value,baseline=old,
                   increment_1_to_2=err2,increment_2_to_4=err4,reaction_error=balance,passes=bool(passes))
        records.append(rec)
        save(cc, r, out/spec['csv'].removesuffix('.csv'), d)
    boundary = [k for k,v in zip(keys,z) if min(v,1-v)<.02]
    screening = all(r['passes'] for r in records if r['role']=='validation')
    status = 'REUSED_HOLDOUT_SCREENING_PASSED_NEW_VALIDATION_REQUIRED' if screening and not boundary else 'FITTED_BUT_NOT_ACCEPTED'
    c['status'] = status
    c['calibration_scope'] = acceptance['scope']
    c['validation_note'] = plan['holdout_policy']
    write_json(out/'reviewed_config.json', c)
    assessment = dict(status=status, numerical_screening_passed=screening,
                      baseline_objective=baseline_objective, selected_training_objective=selected['objective'],
                      selected_start=selected['start'], selected_config_sha256=selected_hash,
                      tests=records, boundary_parameters=boundary,
                      normalized_local_sensitivity_condition=condition,
                      sensitivity_singular_values=sv.tolist(), elapsed_seconds=time.monotonic()-begin,
                      source_files_unchanged=all(digest(root/f)==h for f,h in hashes.items()),
                      holdout_policy=plan['holdout_policy'], parameter_values={k:c['params'][k] for k in keys})
    assert assessment['source_files_unchanged']
    write_json(out/'assessment.json', assessment)
    lines = ['# 철근 Steel02 재보정 — 2026-10-12', '', f'상태: **{status}**', '',
             '기존 원자료·4학습/2검증 분할·탐색범위·합격 기준을 유지했다. 학습 시편 A65의 극값 오차도 기준을 넘었으므로, 학습 곡선과 인장/압축 극값을 함께 고려했다.',
             '목적함수는 시편별 0.8×곡선 NRMSE² + 0.1×인장 극값 오차² + 0.1×압축 극값 오차²의 평균이다. 5개 시작점 중 학습 목적함수로만 선택하고 선택값 저장 후 검증을 실행했다.',
             '**기존 검증 시편의 결과는 이전 작업에서 이미 관측했다. 이번 결과를 새 독립 블라인드 검증으로 부르지 않으며, 검증 결과를 보고 재시도하지 않았다.**', '',
             '| 시편 | 역할 | 이전 NRMSE | 이번 NRMSE | 이전 극값 오차 | 이번 극값 오차 | 기준 통과 |',
             '|---|---|---:|---:|---:|---:|---|']
    for r in records:
        lines.append(f'| {r["csv"]} | {r["role"]} | {r["baseline"]["nrmse"]:.3%} | {r["nrmse"]:.3%} | {r["baseline"]["peak_error"]:.3%} | {r["peak_error"]:.3%} | {r["passes"]} |')
    lines += ['', '| 계수 | 값 |', '|---|---:|']
    for k in keys: lines.append(f'| {k} | {c["params"][k]:.10g} |')
    lines += ['', 'E·Fy는 kPa, 나머지는 무차원. R0·a2·a4 등 고정 계수는 reviewed_config.json에 포함된다.',
              f'경계 부근 계수: {boundary or "없음"}. 정규화 국부 감도 조건수: {condition:.4g} (유일해·신뢰구간을 의미하지 않음).',
              f'최대 증분 세분화 차이: {max(max(r["increment_1_to_2"],r["increment_2_to_4"]) for r in records):.6%}. 최대 반력–응력 상대 차이: {max(r["reaction_error"] for r in records):.3g}.',
              '검토 기준은 기존과 동일: 검증 시편별 NRMSE ≤5%, 극값 오차 ≤10%, 증분 차이 ≤1%. 이번에는 반력–응력 일치도 확인했다.',
              '적용 범위: UCSD Grade 60 A6 그룹의 제공된 반복 단축 이력. 전체 Grade 60 제품, 좌굴·파단·피로·다축 거동 검증이 아니다.',
              '기존 fitted/ 및 reviewed_config.json은 보존했다. 새 후보는 이 폴더의 reviewed_config.json이며 추가 독립 검증 전에는 최종 채택값으로 취급하지 않는다.', '',
              '재현: 패키지 루트에서 `python material_calibration/refine_rebar.py --out material_calibration/rebar_ucsd_a6/새폴더명`.',
              '실험 출처: [NHERI SimCenter qfem-0028](https://nheri-simcenter.github.io/quoFEM-Documentation/common/user_manual/examples/desktop/qfem-0028/README.html). 원본 해시·라이선스는 상위 provenance.json과 raw/rebar_ucsd/SOURCE_LICENSE.txt를 따른다.']
    (out/'RESULTS_KO.md').write_text('\n'.join(lines)+'\n', encoding='utf8')
    print(json.dumps(assessment,ensure_ascii=False,indent=2), flush=True)


if __name__ == '__main__':
    main()
