"""Research-informed diagnostics, not a universal soil qualification standard.

Strain: fraction; stress/modulus: kPa; damping and volume outputs: percent.
See MATERIAL_EVALUATION.md for provenance, scope and all engineering tolerances.
"""
import csv
import json
from pathlib import Path

import numpy as np

SOURCES = {
    "PEER2014": "https://peer.berkeley.edu/publications/2014-16",
    "PEER2008": "https://peer.berkeley.edu/publications/2008-04",
    "PEER2001": "https://peer.berkeley.edu/sites/default/files/0109_j._stewart_s._chiou_j._bray_r._graves_p._somerville_n._abrahamson.pdf",
    "PDMY02": "https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html",
    "Hsu2004": "https://doi.org/10.1061/(ASCE)1090-0241(2004)130:1(58)",
    "Carotenuto2022": "https://www.issmge.org/uploads/publications/1/120/ICSMGE_2022-10.pdf",
}
# These are explicitly project numerical tolerances, NOT literature acceptance limits.
TOL = {
    "probe_relative": 0.02, "refinement_relative": 0.02,
    "work_refinement_relative": 0.05, "stress_balance_relative": 1e-4,
    "strain_absolute": 1e-8, "eta_absolute": 1e-6,
    "loop_stress_closure_relative": 0.01,
    "loop_volume_closure_absolute": 1e-5,
    "volume_refinement_absolute": 1e-5,
}
GROUPS = {
    "강성": ["G0_kPa", "Gsec_last_kPa", "G_ratio_last"],
    "감쇠": ["D_last_percent"],
    "강도": ["tau_peak_observed_kPa"],
    "체적변형": ["eps_v_final_percent"],
    "반복에 따른 변화": ["G_last_over_first", "tau_last_over_first"],
}


def load_targets(path):
    """Each non-null target is an inclusive experimental interval with provenance."""
    with Path(path).open(encoding="utf-8") as f:
        data = json.load(f)
    allowed = {key for keys in GROUPS.values() for key in keys}
    if set(data.get("bounds", {})) - allowed:
        raise ValueError("Unknown target metric; see material_targets.json")
    for key, spec in data.get("bounds", {}).items():
        if spec is None:
            continue
        bounds = np.asarray([spec["lower"], spec["upper"]], dtype=float)
        if not np.all(np.isfinite(bounds)) or bounds[0] > bounds[1] or not spec.get("source"):
            raise ValueError(f"{key}: finite lower <= upper and source are required")
    return data


def integrate_work(tau, gamma):
    # Keep the sign. abs() would hide reversed integration or energy generation.
    return float(np.sum(0.5 * (tau[1:] + tau[:-1]) * np.diff(gamma)))


def cycle_metrics(run):
    """Zero-to-zero imposed cycles; first loading is excluded from damping.

    Gsec uses stresses AT strain reversal points, not unrelated stress extrema.
    D is withheld if stress or volume endpoints drift (a nonstationary loop).
    Work is still reported, but is not called dissipated energy on an open loop.
    """
    out = []
    q = run["steps"]
    g, t = run["gamma"], run["tau"]
    ev = np.sum(run["strain"][:, :3], axis=1)
    for i in range(run["cycles"]):
        a, b = 4 * q * i, 4 * q * (i + 1)
        if b >= len(g):
            break
        pos, neg = a + q, a + 3 * q
        ga = float((g[pos] - g[neg]) / 2)
        ta = float((t[pos] - t[neg]) / 2)
        G = ta / ga if ga > 0 else None
        work = integrate_work(t[a:b+1], g[a:b+1])
        closure = abs(float(t[b] - t[a])) / max(abs(2 * ta), 1e-12)
        volume_closure = abs(float(ev[b] - ev[a]))
        closed = (i > 0 and ga > 0 and ta > 0
                  and closure <= TOL["loop_stress_closure_relative"]
                  and volume_closure <= TOL["loop_volume_closure_absolute"])
        D = 100 * work / (4 * np.pi * 0.5 * ta * ga) if closed else None
        out.append({"cycle": i + 1, "Gsec_kPa": G, "tau_amplitude_kPa": ta,
                    "gamma_amplitude": ga, "shear_work_kPa": work,
                    "stress_closure_ratio": closure,
                    "volume_closure_strain": volume_closure,
                    "D_percent": D,
                    "eps_v_end_percent": float(100 * (ev[b] - ev[0]))})
    return out


def summarize(run):
    cycles = cycle_metrics(run)
    ev = np.sum(run["strain"][:, :3], axis=1)
    m = {key: None for keys in GROUPS.values() for key in keys}
    m.update({"tau_peak_observed_kPa": float(np.max(np.abs(run["tau"]))),
              "eps_v_final_percent": float(100 * (ev[-1] - ev[0])),
              "eps_v_min_percent": float(100 * np.min(ev - ev[0])),
              "eps_v_max_percent": float(100 * np.max(ev - ev[0])),
              "p_initial_kPa": float(-np.mean(run["stress"][0, :3])),
              "K0_x_initial": float(run["stress"][0, 0] / run["stress"][0, 2]),
              "K0_y_initial": float(run["stress"][0, 1] / run["stress"][0, 2]),
              "eta_max": float(np.max(run["eta"])) if len(run["eta"]) else None})
    if cycles:
        first, last = cycles[0], cycles[-1]
        m["Gsec_last_kPa"] = last["Gsec_kPa"]
        m["D_last_percent"] = last["D_percent"]
        if len(cycles) >= 2:
            for name, field in [("G_last_over_first", "Gsec_kPa"),
                                ("tau_last_over_first", "tau_amplitude_kPa")]:
                if first[field] is not None and first[field] > 0:
                    m[name] = last[field] / first[field]
    return m, cycles


def status_of(rows):
    states = [r["status"] for r in rows]
    if "FAIL" in states:
        return "FAIL"
    if not states or any(s != "PASS" for s in states):
        return "NOT_EVALUATED"
    return "PASS"


def evaluate(base, runs, targets):
    metrics, cycles = summarize(base)
    rows = []

    def add(group, name, status, value, criterion, basis):
        rows.append(dict(group=group, name=name, status=status, value=value,
                         criterion=criterion, basis=basis))

    def check(group, name, value, passed, criterion, basis="프로젝트 수치검증 허용오차"):
        good = value is not None and np.all(np.isfinite(value))
        add(group, name, "PASS" if good and passed else "FAIL", value, criterion, basis)

    check("수치 안정성", "전체 재하 완료", int(base["completed"]), base["completed"],
          "요청한 모든 증분 수렴; 중도 실패는 FAIL")
    finite = all(np.all(np.isfinite(base[k])) for k in
                 ["stress", "strain", "reaction", "displacement", "gamma", "tau"])
    check("수치 안정성", "유한 응답", int(finite), finite, "NaN/Inf 없음", "계산 유효성")
    stress_scale = max(base["sigma_v"], float(np.max(np.abs(base["tau"]))), 1.0)
    error = float(np.max(np.abs(base["reaction"][:, 0] - base["tau"])) / stress_scale)
    check("수치 안정성", "전단 반력-재료응력 평형", error,
          error <= TOL["stress_balance_relative"], "max|Δτ| / max(σv, |τ|, 1 kPa) <= 1e-4")
    error = float(np.max(np.abs(base["stress"][:, 2] + base["sigma_v"])) / base["sigma_v"])
    check("수치 안정성", "일정 수직응력 유지", error,
          error <= TOL["stress_balance_relative"], "max|σzz + σv| / σv <= 1e-4")
    g_nodes = (base["displacement"][:, 4, 0] - base["displacement"][:, 0, 0]) / base["height"]
    error = float(np.max(np.abs(base["gamma"] - g_nodes)))
    check("수치 안정성", "절점-재료 전단변형률", error, error <= TOL["strain_absolute"], "절대차 <= 1e-8")
    error = float(np.max(np.abs(base["gamma"] - base["gamma"][0] - base["target"])))
    check("수치 안정성", "목표 변형률 도달", error, error <= TOL["strain_absolute"], "절대차 <= 1e-8")
    positive_G = bool(cycles) and all(c["Gsec_kPa"] is not None and c["Gsec_kPa"] > 0 for c in cycles)
    check("강성", "양의 반복 할선강성", metrics["Gsec_last_kPa"], positive_G,
          "이 시험의 모든 완결 사이클 Gsec > 0 kPa", "기초 물리 일관성; 실험 적합성 아님")
    if len(base["eta"]) == len(base["gamma"]):
        eta = base["eta"]
        check("강도", "PDMY02 강도면 일관성", metrics["eta_max"],
              np.min(eta) >= -TOL["eta_absolute"] and np.max(eta) <= 1 + TOL["eta_absolute"],
              "0 <= eta_r <= 1 (수치 여유 1e-6)", "PDMY02 출력 정의; 독립 강도검증 아님")
    else:
        add("강도", "PDMY02 강도면 일관성", "NOT_EVALUATED", None,
            "7번째 stress 성분 필요", "PDMY02")

    probe_values = []
    for name in ["probe", "probe_half"]:
        run = runs.get(name)
        if run is not None and run.get("completed"):
            q = run["steps"]
            # Independent virgin tiny-strain secant; no contamination by 1% history.
            probe_values.append(float((run["tau"][q] - run["tau"][0]) /
                                      (run["gamma"][q] - run["gamma"][0])))
    if len(probe_values) == 2:
        error = abs(probe_values[0] - probe_values[1]) / max(abs(probe_values[1]), 1e-12)
        good_probe = min(probe_values) > 0 and error <= TOL["probe_relative"]
        check("수치 안정성", "미소변형 G0 진폭 수렴", error, good_probe,
              "gamma=1e-7, 5e-8 독립 시험 G0 차이 <= 2%, 두 G0 > 0")
        if good_probe:
            metrics["G0_kPa"] = probe_values[1]
            metrics["Vs_equivalent_m_s"] = float(np.sqrt(probe_values[1] / base["rho"]))
            if metrics["Gsec_last_kPa"] is not None:
                metrics["G_ratio_last"] = metrics["Gsec_last_kPa"] / probe_values[1]
    else:
        failed = any(name in runs and not runs[name].get("completed") for name in ["probe", "probe_half"])
        add("수치 안정성", "미소변형 G0 진폭 수렴", "FAIL" if failed else "NOT_EVALUATED", None,
            "독립적인 두 미소변형 시험 필요 (--quick이면 생략)", "프로젝트 수치검증")

    for left_name, right_name, left, right in [
        ("base", "fine", base, runs.get("fine")),
        ("fine", "finer", runs.get("fine"), runs.get("finer")),
    ]:
        name = f"증분 수렴 {left_name}->{right_name}"
        if left is None or right is None:
            add("수치 안정성", name, "NOT_EVALUATED", None, "2배/4배 세분화 필요", "프로젝트 수치검증")
            continue
        if not left.get("completed") or not right.get("completed"):
            add("수치 안정성", name, "FAIL", None, "두 해석 모두 완료 필요", "프로젝트 수치검증")
            continue
        lm, lc = summarize(left)
        rm, rc = summarize(right)
        # Full stress path AND per-cycle quantities; identical target history at nested points.
        error = float(np.max(np.abs(left["tau"] - right["tau"][::2])) /
                      max(np.max(np.abs(right["tau"])), 1e-12))
        check("수치 안정성", name + " 응력경로", error,
              error <= TOL["refinement_relative"], "max|Δτ| / max|τ_fine| <= 2%")
        for metric in ["tau_peak_observed_kPa", "Gsec_last_kPa"]:
            error = abs(lm[metric] - rm[metric]) / max(abs(rm[metric]), 1e-12)
            check("수치 안정성", name + " " + metric, error,
                  error <= TOL["refinement_relative"], "세분화 결과 대비 상대차 <= 2%")
        work_errors = [abs(a["shear_work_kPa"] - b["shear_work_kPa"]) /
                       max(abs(b["shear_work_kPa"]), 1e-8) for a, b in zip(lc, rc)]
        error = max(work_errors)
        check("수치 안정성", name + " 사이클 전단일", error,
              error <= TOL["work_refinement_relative"], "모든 사이클 차이 <= 5%; 분모 하한 1e-8 kPa")
        lev = np.sum(left["strain"][:, :3], axis=1)
        rev = np.sum(right["strain"][::2, :3], axis=1)
        error = float(np.max(np.abs((lev - lev[0]) - (rev - rev[0]))))
        check("수치 안정성", name + " 체적변형 경로", error,
              error <= TOL["volume_refinement_absolute"], "전단 시작점 기준 절대차 <= 1e-5 (=0.001%p)")

    last_D = metrics["D_last_percent"]
    if last_D is None:
        add("감쇠", "닫힌 반복 루프", "NOT_EVALUATED", None,
            "첫 사이클 제외; 응력 폐합 오차 <=1%, 체적변형 드리프트 <=1e-5 필요",
            "PEER2001 에너지 정의 + 프로젝트 폐합 허용오차")
    else:
        check("감쇠", "비음수 루프 에너지", last_D, last_D >= -1e-8,
              "닫힌 루프 D >= 0%; 임의 상한 없음", "에너지 일관성; 실험 적합성 아님")
    add("강도", "최대강도 도달 여부", "NOT_EVALUATED", metrics["tau_peak_observed_kPa"],
        "현재 값은 반복시험 중 관측 최대응력; 단조전단 강도시험 필요", "PEER2008/2014")
    gamma_percent = 100 * base["amplitude"]
    add("강성/감쇠", "일반 경험곡선 적용범위", "NOT_APPLICABLE" if gamma_percent > 0.5 else "INFO",
        gamma_percent, "일반 MRD 관계는 약 0.3~0.5% 초과에서 신뢰 제한; 이 범위 이하도 자동 합격 아님", "PEER2014")
    add("체적변형", "문헌의 모래 체적변형 문턱", "INFO", gamma_percent,
        "Hsu & Vucetic: gamma_tv 약 0.01~0.02%; 전단진폭 문턱이며 허용 체적변형률 아님",
        "Hsu2004; 대상 흙의 동일성 미확인, 합격판정에 전용하지 않음")
    add("반복에 따른 변화", "장기 반복/잔류전단변형", "NOT_APPLICABLE", base["cycles"],
        "단기 변형률제어 시험으로 자유 잔류전단변형 및 장기 수명 판정 불가", "시험 경계조건 / Carotenuto2022")

    # Target comparisons require explicit confirmation of matching specimen/test conditions.
    context = targets.get("context", {})
    matching = (context.get("conditions_confirmed") is True
                and bool(context.get("soil_description"))
                and context.get("test") == "drained_constant_sigma_v_simple_shear"
                and context.get("cycles") == base["cycles"]
                and isinstance(context.get("gamma_amplitude"), (int, float))
                and np.isclose(context["gamma_amplitude"], base["amplitude"], rtol=1e-6, atol=0)
                and isinstance(context.get("sigma_v_kPa"), (int, float))
                and np.isclose(context["sigma_v_kPa"], base["sigma_v"], rtol=1e-6, atol=0))
    validation = []
    for group, keys in GROUPS.items():
        group_rows = []
        for key in keys:
            spec = targets.get("bounds", {}).get(key)
            value = metrics[key]
            if spec is None:
                status, criterion, basis = "NOT_EVALUATED", "동일 조건 실험/검증자료의 lower, upper 필요", "비교자료 미입력"
            elif not matching:
                status, criterion, basis = "NOT_EVALUATED", "targets.context 시험 조건 확인/일치 필요", spec["source"]
            elif (not base["completed"] or value is None or not np.isfinite(value)
                  or any(r["status"] == "FAIL" for r in rows)):
                status, criterion, basis = "NOT_EVALUATED", "완결되고 유효한 측정값 필요", spec["source"]
            else:
                status = "PASS" if spec["lower"] <= value <= spec["upper"] else "FAIL"
                criterion, basis = f'{spec["lower"]} <= {key} <= {spec["upper"]}', spec["source"]
            item = dict(group=group, name=key, status=status, value=value, criterion=criterion, basis=basis)
            group_rows.append(item)
            validation.append(item)
        add(group, "실험자료 적합성", status_of(group_rows), None,
            "해당 항목의 모든 지정 지표를 비교; 상세 판정은 target_checks 참조", "독립 비교자료")
    numerical = status_of([r for r in rows if r["group"] == "수치 안정성"])
    target_status = status_of(validation)
    consistency = "FAIL" if any(r["status"] == "FAIL" for r in rows) else "NO_FAILURE_DETECTED"
    report = {"scope": "현재 단일 압력/진폭 시험에 한정; 지반 재료 전체 사용 승인 아님",
              "numerical_status": numerical, "target_status": target_status,
              "consistency_status": consistency,
              "metrics": metrics, "cycles": cycles, "checks": rows,
              "target_checks": validation, "targets": targets,
              "numerical_tolerances": TOL, "sources": SOURCES,
              "run": {k: base[k] for k in ["completed", "error", "steps", "cycles", "amplitude", "sigma_v", "rho", "material", "opensees_version"]},
              "auxiliary_runs": {name: {"completed": r.get("completed"), "error": r.get("error")} for name, r in runs.items()}}
    return report


def json_safe(value):
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    return value


def save_report(report, base, runs, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    report = json_safe(report)
    (directory / "evaluation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    for name, run in {"base": base, **runs}.items():
        if "gamma" in run:
            np.savez_compressed(directory / f"{name}_history.npz", **{k: run[k] for k in
                                ["gamma", "tau", "stress", "strain", "eta", "reaction", "target", "displacement"]})
    with (directory / "cycles.csv").open("w", newline="", encoding="utf-8-sig") as f:
        if report["cycles"]:
            writer = csv.DictWriter(f, fieldnames=list(report["cycles"][0]))
            writer.writeheader()
            writer.writerows(report["cycles"])
    lines = ["PDMY02 재료 평가", "PASS=해당 검사 충족 / FAIL=미충족 / NOT_EVALUATED=자료·시험 부족 / NOT_APPLICABLE=적용 불가",
             "수치검사 PASS는 대상 지반 적합성을 뜻하지 않습니다.",
             f'수치 안정성: {report["numerical_status"]}', f'실험자료 적합성: {report["target_status"]}', "", "계산 지표:"]
    for key, value in report["metrics"].items():
        lines.append(f"  {key}: {value if value is not None else '평가 불가'}")
    lines.append("\n항목별 검사:")
    for row in report["checks"] + report["target_checks"]:
        lines.append(f'[{row["status"]}] {row["group"]} / {row["name"]}: {row["value"]}\n'
                     f'  기준: {row["criterion"]}\n  근거: {row["basis"]}')
    lines.extend(["", "문헌:"] + [f"{k}: {v}" for k, v in SOURCES.items()])
    output = "\n".join(lines)
    (directory / "evaluation.txt").write_text(output, encoding="utf-8")
    print("\n" + output)
    print(f"\n평가/원시 결과 저장: {directory.resolve()}")
