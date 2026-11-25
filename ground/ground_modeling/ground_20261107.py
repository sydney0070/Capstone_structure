"""InSight 착륙지 화성 지반 단독 해석.

main()은 입력을 읽고 지반을 조립한 뒤 자중·지진 해석과 시각화를 수행한다.
지반 생성과 해석 함수는 SSI 모델에서도 쓴다.
"""

import argparse
import csv
import json
import math
import os
import sys
import tempfile
import time
import uuid
from collections import defaultdict
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import numpy as np
import openseespy.opensees as ops

from ground_hybrid import default_config
from ground_hybrid.io import serial
# SSI와 연구 스크립트가 사용하는 지반 함수의 기존 import 경로를 유지한다.
from ground_hybrid.mars_inputs import (
    file_sha256, insight_reference_site, inspect_site, prepare_motion, station_georeference,
)
from ground_hybrid.mars_results import (
    GroundHistory, ground_state, remove_base_rigid_motion, require_finite,
    save_response_viewer, write_ground_viewer,
)
from ground_hybrid.run_journal import RunJournal, atomic_json
from insight_motion import INSIGHT_ACC_PATH

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]


# 지반 단독 해석

def main(argv=None, *, column=False):
    """지반을 조립하고 자중·지진 응답을 계산해 실행 폴더에 저장한다."""
    raw_args = sys.argv[1:] if argv is None else argv
    cli_parser = parser(column=column)
    options = cli_parser.parse_args(raw_args)
    if options.quick:
        if '--warmup' not in raw_args:
            options.warmup = 1.0
        if '--duration' not in raw_args:
            options.duration = 1.0

    positive_inputs = [options.depth, options.duration, options.dt, options.max_runtime]
    if not column:
        positive_inputs.extend((options.ground_width, options.spacing))
    if (
        not np.isfinite((*positive_inputs, options.warmup, options.freqmax)).all()
        or min(positive_inputs) <= 0
        or options.warmup < 0
        or not 0.05 < options.freqmax < 5
        or options.dt * options.freqmax > 0.05 + 1e-12
        or options.static_steps < 1
        or options.max_frames < 2
        or (
            options.visual_scale is not None
            and (not math.isfinite(options.visual_scale) or options.visual_scale <= 0)
        )
    ):
        cli_parser.error('invalid finite geometry, time or display parameters')

    output_dir = options.output_dir.resolve() if options.output_dir else new_output(column=column)
    if output_dir.exists() and any(output_dir.iterdir()):
        cli_parser.error('choose a fresh output directory')
    output_dir.mkdir(parents=True, exist_ok=True)
    started_at = time.perf_counter()
    result = dict(
        model='Mars ground only (1D column)' if column else 'Mars ground only',
        complete=False,
        safe_verdict=False,
        parameters={**vars(options), 'output_dir': str(output_dir)},
        structure_count=0,
        foundation_count=0,
        contact_count=0,
    )

    def save_summary():
        result['elapsed_s'] = time.perf_counter() - started_at
        atomic_json(output_dir / 'summary.json', serial(result))

    ops.wipe()
    ops.logFile(str(output_dir / 'opensees.log'))
    try:
        # 실행에 사용한 소스와 물성 설정을 나중에 대조할 수 있도록 해시를 남긴다.
        sources = [
            Path(__file__).resolve(),
            HERE / 'insight_dynamic_viewer.py',
            HERE / 'insight_motion.py',
            HERE / 'insight_site_dem.py',
            HERE / 'ground_modeling.py',
            PROJECT_ROOT / 'material/materials_20261107.py',
            *(HERE / 'ground_hybrid').rglob('*.py'),
            *(HERE / 'ground_hybrid/configs').glob('*.json'),
        ]
        if column:
            sources.append(HERE / 'ground_1d_column.py')
        result['code_sha256'] = {str(path): file_sha256(path) for path in sources}
        result['runtime'] = dict(python=sys.version, opensees=ops.version(), numpy=np.__version__)

        # 1. 데이터 불러오기
        result['site'], dem_samples = insight_reference_site(options.ground_width)
        np.savetxt(
            output_dir / 'site_dem_samples.csv',
            dem_samples,
            delimiter=',',
            header='east_m,north_m,dem_elevation_m,dem_relative_m,modeled_z_m',
            comments='',
        )
        motion_input = None
        if not options.static_only:
            motion_input = prepare_motion(
                INSIGHT_ACC_PATH, output_dir, duration=options.duration,
                warmup=options.warmup, freqmax=options.freqmax,
            )
            result['motion'] = motion_input[2]
            if any(channel.split('.')[1] != 'ELYSE' for channel in result['motion']['channels']):
                raise ValueError('ground site and motion station mismatch')
        config = default_config('Mars',mars_profile=options.material_profile)
        result['material_profile']=dict(name=options.material_profile,scenario=config['scenario'],selection_status=config.get('selection_status','historical_reference'))
        config.update(
            depth_m=options.depth, max_vertical_spacing_m=4.0, max_elements=20000, elastic=False,
        )

        # 2. 모델링
        if column:
            from ground_1d_column import build_1d_column

            ground = build_1d_column(config)
        else:
            config.update(
                width_m=options.ground_width, length_m=options.ground_width,
                core_width_m=0.7 * options.ground_width, core_length_m=0.7 * options.ground_width,
                core_spacing_m=options.spacing, max_outer_spacing_m=options.spacing,
            )
            ground = build_ground_20261107(config)
        node_tags = sorted(ground['nodes'])

        # 측면 감쇠 요소를 추가하기 전에는 3자유도 지반 절점과 bbarBrick만 확인한다.
        if (
            set(ops.getEleTags()) != set(ground['elements'])
            or any(len(ops.nodeDisp(tag)) != 3 for tag in node_tags)
            or ground['attachments']
        ):
            raise RuntimeError('ground-only domain contains non-ground model components')
        result.update(
            ground_elements=len(ground['elements']),
            ground_nodes=len(node_tags),
            mesh=ground['mesh'],
            ground_materials={
                name: ground['materials'][name] for name in ground['material_volume_fractions']
            },
            ground_layers=config['layers'],
            material_volume_fractions=ground['material_volume_fractions'],
            ground_revision=ground['revision_validation'],
            assembly_gate='only 3DOF ground nodes and bbarBrick ground elements before boundary activation',
        )
        result['resolution'] = resolution(ground, options.freqmax, options.dt)
        save_summary()
        journal = RunJournal(
            output_dir,
            dict(
                model=result['model'],
                parameters=result['parameters'],
                code_sha256=result['code_sha256'],
            ),
        )

        # 3. 해석: 자중 평형을 구하고, 그 상태에서 지진 하중을 적용한다.
        history = None
        if options.ground_history:
            history = GroundHistory(ground, output_dir / 'static_ground_history.csv')
            history.capture()
        result['gravity'] = static_ramp(
            ground,
            ground['gravity_loads'],
            steps=options.static_steps,
            on_commit=history.capture if history else None,
            attempt_observer=journal,
        )
        if history:
            result['static_history'] = history.save(output_dir / 'static_ground_history.csv')
        reference_displacements = require_finite(
            [ops.nodeDisp(tag)[:3] for tag in node_tags], 'ground gravity'
        )
        frames = [reference_displacements.copy()]
        frame_times = [0.0]
        result['boundary'] = {'kind': ground['boundary']}

        if motion_input is not None:
            dynamic_frames, dynamic_times = _run_ground_motion(
                ground, options, output_dir, node_tags, reference_displacements, motion_input, result
            )
            frames.extend(dynamic_frames)
            frame_times.extend(dynamic_times)

        # 4. 시각화 및 결과 저장
        result['final_ground_state'] = ground_state(ground, gauss_points=8)
        if not options.no_html:
            coordinates = np.asarray([ground['nodes'][tag] for tag in node_tags])
            result['viewer'] = write_ground_viewer(
                ground, options, output_dir, node_tags, coordinates,
                frames, frame_times, result['site'],
            )
        result['complete'] = True
        save_summary()
        print('[complete ground only]', output_dir, flush=True)
    except BaseException as exc:
        result['error'] = str(exc)
        save_summary()
        raise
    finally:
        ops.wipe()


def resolution(g, freqmax, dt):
    """최고 입력 주파수의 파장과 주기를 메시·시간 간격과 비교한다."""
    from ground_hybrid.dynamics import mesh_frequency_guidance

    guidance = mesh_frequency_guidance(g, freqmax)
    guidance.update(
        integration_dt_s=dt,
        samples_per_input_cycle=1 / (dt * freqmax),
        temporal_bandwidth_pass=dt <= 1 / (20 * freqmax) + 1e-12,
        scope='Frequency/mesh heuristic using reference G,K; not a time/mesh convergence proof',
    )
    if not guidance['mesh_target_pass']:
        raise ValueError(
            'Input bandwidth exceeds 10-elements/wavelength mesh target; '
            'reduce --freqmax or refine --spacing'
        )
    if not guidance['temporal_bandwidth_pass']:
        raise ValueError('At least 20 integration samples per highest input cycle required')
    return guidance


# 지반 요소 조립과 자중 평형

def build_ground_20261107(config, terrain=None, *, roller_sides=True):
    """호출자가 관리하는 OpenSees 도메인에 지반 요소를 추가한다."""
    from ground_hybrid.mesh import build_ground

    return validate_ground_20261107(build_ground(config, terrain, roller_sides=roller_sides))


def validate_ground_20261107(g):
    """직교·기초 형상 일치 메시 모두 같은 재료 검증을 적용한다."""
    sys.path.insert(0, str(PROJECT_ROOT / 'material'))
    from materials_20261107 import validate_record

    for name, record in g['materials'].items():
        validate_record(name, record)
    g['revision_validation'] = dict(
        date='20261107',
        materials_checked=list(g['materials']),
        max_element_aspect_ratio=g['mesh']['max_aspect_ratio'],
        mass_kg=g['mass_kg'],
        site_claim=g['config'].get(
            'site_scope',
            'Synthetic reference unless caller supplies independently registered terrain/geology',
        ),
        constitutive_scope=(
            'DP pressure-independent elastic stiffness + frictional plasticity; '
            'no G/Gmax/damping curve refit'
        ),
    )
    return g


def static_ramp(
    g,
    loads,
    *,
    pattern=800001,
    steps=20,
    max_attempts=150,
    algorithms=("KrylovNewton", "NewtonLineSearch", "Newton"),
    on_commit=None,
    attempt_observer=None,
    initial_increment=None,
    linear_system='UmfPack',
    first_step_algorithm=None,
):
    """자중 또는 후속 정적 하중을 증분 재하하고 지반 지점의 평형을 확인한다."""
    from ground_hybrid.solver import continue_load, force_balance

    if linear_system not in ('UmfPack', 'SuperLU'):
        raise ValueError('validated sparse unsymmetric system required')
    if any(
        material['active_law'] == 'PressureDependMultiYield02'
        for material in g['materials'].values()
    ):
        raise ValueError('Use separately validated PDMY stage initialization')

    ops.timeSeries('Linear', pattern)
    ops.pattern('Plain', pattern, pattern)
    for node_tag, load in loads.items():
        dof_count = len(ops.nodeDisp(node_tag))
        vector = list(load)
        if len(vector) > dof_count or len(vector) < 3 or not np.isfinite(vector).all():
            raise ValueError('Invalid load dimensions/values')
        # 구조물 절점에는 필요에 따라 회전 자유도의 하중을 0으로 채운다.
        ops.load(node_tag, *(vector + [0.0] * (dof_count - len(vector))))

    ops.wipeAnalysis()
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system(linear_system)
    force_scale = sum(float(np.linalg.norm(load[:3])) for load in loads.values())
    ops.test('NormUnbalance', max(0.001, force_scale * 1e-10), 80)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.0)
    ops.analysis('Static')
    attempt_count = 0

    def report_progress(attempt):
        nonlocal attempt_count
        attempt_count += 1
        print(
            f"[static {pattern}] factor={attempt['after']:.4f}, "
            f"increment={attempt['increment']:.3g}",
            flush=True,
        )
        if attempt['code'] == 0 and on_commit is not None:
            on_commit()
        if attempt_count > max_attempts:
            raise RuntimeError('Static load attempt limit exceeded; no strength or tolerance changed')

    attempts = continue_load(
        pattern,
        steps=steps,
        min_increment=1e-5,
        progress=report_progress,
        algorithms=algorithms,
        attempt_observer=attempt_observer,
        initial_increment=initial_increment,
        max_attempts=max_attempts,
        first_step_algorithm=first_step_algorithm,
    )
    if abs(ops.getLoadFactor(pattern) - 1) > 1e-8:
        raise RuntimeError('Incomplete load factor')

    # 재하가 끝난 뒤 하중을 고정하고, 증분 0에서 최종 평형을 한 번 확인한다.
    ops.loadConst('-time', 0.0)
    ops.integrator('LoadControl', 0.0)
    if ops.analyze(1):
        raise RuntimeError('Zero-step final equilibrium failed')

    applied_loads = g.setdefault('applied_loads', {})
    for node_tag, load in loads.items():
        previous = applied_loads.get(node_tag, [])
        accumulated = np.zeros(max(len(previous), len(load)))
        accumulated[:len(previous)] += previous
        accumulated[:len(load)] += load
        applied_loads[node_tag] = accumulated

    balance = force_balance(g)
    applied_moment = sum(
        (
            np.cross(np.array(ops.nodeCoord(node_tag)), np.asarray(load[:3]))
            + np.pad(np.array(load[3:6]), (0, 3 - len(load[3:6])))
            for node_tag, load in applied_loads.items()
        ),
        np.zeros(3),
    )
    support_force = np.zeros(3)
    support_moment = np.zeros(3)
    for node_tag, fixity in g['fixities'].items():
        reaction = np.asarray(ops.nodeReaction(node_tag)[:3]) * np.array(fixity)
        support_force += reaction
        support_moment += np.cross(np.array(ops.nodeCoord(node_tag)), reaction)

    # 지반 지점의 병진 반력을 합산한다. 접촉 보조 절점의 회전 반력은 생략한다.
    force_scale = max(sum(np.linalg.norm(load[:3]) for load in applied_loads.values()), 1.0)
    applied_force = np.sum([load[:3] for load in applied_loads.values()], axis=0)
    force_error = float(np.linalg.norm(support_force + applied_force) / force_scale)
    moment_scale = force_scale * max(g['config']['width_m'], 1.0)
    moment_error = float(np.linalg.norm(support_moment + applied_moment) / moment_scale)
    if force_error > 1e-5:
        raise RuntimeError('Physical ground support force equilibrium failed')

    g['stage'] = 'gravity_equilibrated' if pattern == 800001 else 'static_loaded'
    ops.wipeAnalysis()
    return dict(
        final_factor=float(ops.getLoadFactor(pattern)),
        attempts=len(attempts),
        linear_system=linear_system,
        force_balance=balance,
        physical_support_force_error=force_error,
        physical_support_moment_error=moment_error,
        moment_caveat=(
            'Contact bearing auxiliary rotational reactions excluded; '
            'moment residual is diagnostic for contact'
        ),
    )


# 동적 경계 전환과 시간이력 적분

def lateral_dashpots(g, tag=1000000):
    """고정 바닥을 유지하며 측면 구속을 점성 경계로 전환한다."""
    if g['boundary'] != 'static_box':
        raise ValueError('Boundary transition allowed once')

    base_nodes = {
        node for face in g['faces'] if face['side'] == 'bottom' for node in face['nodes']
    }
    coefficients = defaultdict(lambda: np.zeros(3))
    for face in g['faces']:
        if face['side'] not in ('xmin', 'xmax', 'ymin', 'ymax'):
            continue
        material = g['materials'][face['material']]
        density = material['density_kg_m3']
        shear_speed = math.sqrt(material['G_Pa'] / density)
        pressure_speed = math.sqrt(
            (material['K_Pa'] + 4 * material['G_Pa'] / 3) / density
        )
        normal_axis = 0 if face['side'][0] == 'x' else 1

        # 절점 분담면적을 곱한 파동 임피던스를 경계면별로 합산한다.
        for node_tag, area in zip(face['nodes'], face['weights_m2']):
            if node_tag in base_nodes:
                continue
            damping = np.full(3, density * shear_speed * area)
            damping[normal_axis] = density * pressure_speed * area
            coefficients[node_tag] += damping

    # 해제할 자유도의 정적 반력을 외력으로 옮겨 기존 자중 평형을 유지한다.
    ops.loadConst('-time', 0.0)
    ops.reactions()
    pattern = tag + 9000
    ops.timeSeries('Constant', pattern)
    ops.pattern('Plain', pattern, pattern)
    released = []
    for node_tag, fixed in g['fixities'].items():
        if node_tag in base_nodes:
            continue
        reactions = ops.nodeReaction(node_tag)
        load = np.zeros(3)
        for axis, is_fixed in enumerate(fixed):
            if is_fixed:
                load[axis] = reactions[axis]
                released.append((node_tag, axis + 1))
        if any(fixed):
            ops.load(node_tag, *load)
    for node_tag, dof in released:
        ops.remove('sp', node_tag, dof)

    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for index, (node_tag, damping) in enumerate(sorted(coefficients.items())):
        anchor = tag + index
        if anchor in ops.getNodeTags():
            raise ValueError('Dashpot tag collision')
        ops.node(anchor, *ops.nodeCoord(node_tag))
        ops.fix(anchor, 1, 1, 1)
        for dof, coefficient in enumerate(damping, 1):
            element = tag + 10000 + 3 * index + dof
            ops.uniaxialMaterial('Viscous', element, float(coefficient), 1.0)
            ops.element('zeroLength', element, node_tag, anchor, '-mat', element, '-dir', dof)

    # Newmark 적분 전에 재평형을 확정해 이전 불평형력이 초기 가속도에 남지 않게 한다.
    physical_nodes = list(g['nodes'])
    before = np.array([ops.nodeDisp(node_tag)[:3] for node_tag in physical_nodes])
    ops.wipeAnalysis()
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('UmfPack')
    tolerance = max(0.001, g['mass_kg'] * g['config']['gravity_m_s2'] * 1e-10)
    ops.test('NormUnbalance', tolerance, 80)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.0)
    ops.analysis('Static')
    if ops.analyze(1):
        raise RuntimeError('Boundary transition equilibrium failed')

    after = np.array([ops.nodeDisp(node_tag)[:3] for node_tag in physical_nodes])
    movement = float(np.max(np.linalg.norm(after - before, axis=1)))
    if movement > 1e-8:
        raise RuntimeError(f'Boundary transition changed static state: {movement}')
    ops.wipeAnalysis()
    g['boundary'] = 'fixed_base_viscous_sides'
    return dict(
        kind=g['boundary'],
        released_dofs=len(released),
        side_nodes=len(coefficients),
        base_nodes_retained=len(base_nodes),
        transition_movement_m=movement,
        caveat=(
            'Reference impedances, fixed base, no free-field boundary; '
            'static reactions conserved; rest/attenuation checks required'
        ),
    )


def _run_ground_motion(
    ground, options, output_dir, node_tags, reference_displacements, motion_input, result,
):
    """자중 평형 이후 지진을 가하고 성공한 적분 단계의 응답을 기록한다."""
    if options.boundary == 'viscous-sides':
        result['boundary'] = lateral_dashpots(ground)
    ops.loadConst('-time', 0.0)

    input_dt, acceleration, motion = motion_input

    # 동·북·상 성분을 X·Y·Z 방향 바닥 가진으로 적용한다.
    for axis in range(3):
        tag = 900000 + axis
        ops.timeSeries('Path', tag, '-dt', input_dt, '-values', *acceleration[:, axis].tolist())
        ops.pattern('UniformExcitation', tag, axis + 1, '-accel', tag)
    total_time = (len(acceleration) - 1) * input_dt
    step_count = int(round(total_time / options.dt))
    if not math.isclose(step_count * options.dt, total_time, abs_tol=1e-8):
        raise ValueError('dt must divide selected motion duration')

    ops.wipeAnalysis()
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('UmfPack')
    tolerance = max(0.001, ground['mass_kg'] * ground['config']['gravity_m_s2'] * 1e-10)
    ops.test('NormUnbalance', tolerance, 80)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')

    history = None
    if options.ground_history:
        history = GroundHistory(ground, output_dir / 'dynamic_ground_history.csv')
        history.capture()

    attempts = []
    peak_increment = 0.0
    committed_samples = 0
    frames = []
    frame_times = []
    frame_steps = set(
        np.linspace(0, step_count, min(step_count + 1, options.max_frames), dtype=int)
    )
    with (output_dir / 'committed_response.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_s', 'target_time_s', 'ground_max_increment_m'])

        def record_committed_step():
            nonlocal peak_increment, committed_samples
            current_displacements = require_finite(
                [ops.nodeDisp(tag)[:3] for tag in node_tags], 'ground dynamic'
            )
            # 변위 증분은 자중 평형 기준이며, 성공한 세분 단계도 모두 기록한다.
            increment = float(np.linalg.norm(
                current_displacements - reference_displacements, axis=1
            ).max())
            peak_increment = max(peak_increment, increment)
            committed_samples += 1
            current_time = ops.getTime()
            writer.writerow([current_time, current_time - motion['target_start_model_s'], increment])
            stream.flush()
            if history:
                history.capture()

        for step in range(step_count):
            advance(options.dt, attempts=attempts, on_commit=record_committed_step)
            if not options.no_html and step + 1 in frame_steps:
                displacement = require_finite(
                    [ops.nodeDisp(tag)[:3] for tag in node_tags], 'viewer frame'
                )
                frames.append(displacement.copy())
                frame_times.append(ops.getTime())

    if history:
        result['dynamic_history'] = history.save(output_dir / 'dynamic_ground_history.csv')
    result['dynamic'] = dict(
        steps=step_count,
        committed_samples=committed_samples,
        max_increment_m=peak_increment,
        attempts=attempts,
        scope='every committed substep; ground response only; relative to imposed base excitation',
    )
    return frames, frame_times


def advance(dt, *, max_halvings=6, level=0, attempts=None, on_commit=None):
    """한 단계가 수렴하지 않으면 같은 시간 구간을 반으로 나눠 다시 적분한다."""
    if attempts is None:
        attempts = []
    start_time = ops.getTime()
    code = ops.analyze(1, dt)
    attempts.append(dict(t=start_time, dt=dt, code=code))
    if code:
        if abs(ops.getTime() - start_time) > 1e-9:
            raise RuntimeError('Failed step changed committed time')
        if level >= max_halvings:
            raise RuntimeError(
                'Dynamic interval failed after bounded substeps; no physical failure inference'
            )
        advance(
            dt / 2, max_halvings=max_halvings, level=level + 1,
            attempts=attempts, on_commit=on_commit,
        )
        advance(
            dt / 2, max_halvings=max_halvings, level=level + 1,
            attempts=attempts, on_commit=on_commit,
        )
    elif abs(ops.getTime() - start_time - dt) > 1e-8:
        raise RuntimeError('Unexpected committed step time')
    elif on_commit is not None:
        on_commit()
    return attempts


# 실행 옵션과 프로세스 관리

def parser(*, column=False):
    cli_parser = argparse.ArgumentParser(
        description=(
            "Mars 1D soil column: nonlinear gravity and ELYSE motion."
            if column else "Mars ground only: nonlinear gravity and ELYSE motion."
        )
    )
    numeric_options = (
        ('depth', 12.0),
        ('duration', 2.0),
        ('warmup', 5.0),
        ('dt', 0.05),
        ('freqmax', 1.0),
        ('max-runtime', 180.0),
    )
    for name, default in numeric_options:
        cli_parser.add_argument(f'--{name}', type=float, default=default)
    if column:
        cli_parser.set_defaults(ground_width=60.0, boundary='tied-column')
    else:
        cli_parser.add_argument('--ground-width', type=float, default=60.0)
        cli_parser.add_argument('--spacing', type=float, default=5.0)
        cli_parser.add_argument('--boundary', choices=['viscous-sides', 'fixed-box'], default='viscous-sides')
    cli_parser.add_argument('--static-steps', type=int, default=100)
    cli_parser.add_argument('--material-profile',choices=['nonlinear_20261125','reference'],default='nonlinear_20261125')
    cli_parser.add_argument('--max-frames', type=int, default=101)
    cli_parser.add_argument(
        '--quick', action='store_true', help='1 s warmup + 1 s target; connection check only'
    )
    cli_parser.add_argument('--static-only', action='store_true')
    cli_parser.add_argument('--ground-history', action='store_true')
    cli_parser.add_argument('--no-html', action='store_true')
    cli_parser.add_argument('--visual-scale', type=float, help='display only; default automatic')
    cli_parser.add_argument('--output-dir', type=Path)
    cli_parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    return cli_parser


def new_output(*, column=False):
    name = f"run_{datetime.now():%Y%m%d_%H%M%S}_{uuid.uuid4().hex[:6]}"
    model_folder = 'ground_1d_column' if column else 'ground_20261109'
    return HERE.parent / 'results' / model_folder / name


def supervised_cli(*, column=False):
    """같은 파일을 해석 프로세스로 실행하고 제한 시간과 종료 상태를 관리한다."""
    from ground_hybrid.process_control import MANAGED_ENV, run_isolated

    args = sys.argv[1:]
    if any(flag in args for flag in ('--worker', '--help', '-h')):
        return main(args, column=column)

    cli_parser = parser(column=column)
    options = cli_parser.parse_args(args)
    if not math.isfinite(options.max_runtime) or options.max_runtime <= 0:
        cli_parser.error('max-runtime must be positive')
    output_dir = options.output_dir.resolve() if options.output_dir else new_output(column=column)
    if output_dir.exists() and any(output_dir.iterdir()):
        cli_parser.error('choose a fresh output directory')

    script = HERE / 'ground_1d_column.py' if column else Path(__file__).resolve()
    command = [sys.executable, str(script), *args, '--worker']
    if options.output_dir is None:
        command.extend(['--output-dir', str(output_dir)])
    status = run_isolated(
        command,
        options.max_runtime,
        join_managed_group=os.environ.get(MANAGED_ENV) == '1',
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_json(output_dir / 'run_status.json', dict(status, safe_verdict=False))
    return status['returncode']


def cli(*, column=False):
    if '--inspect-site' not in sys.argv[1:]:
        return supervised_cli(column=column)

    site_parser = argparse.ArgumentParser(description='Ground-only site audit')
    site_parser.add_argument('--inspect-site', action='store_true')
    site_parser.add_argument('--output', type=Path)
    options = site_parser.parse_args()
    report = json.dumps(inspect_site(), indent=2)
    if options.output:
        options.output.parent.mkdir(parents=True, exist_ok=True)
        options.output.write_text(report)
    print(report)
    return 0


if __name__ == '__main__':
    raise SystemExit(cli())
