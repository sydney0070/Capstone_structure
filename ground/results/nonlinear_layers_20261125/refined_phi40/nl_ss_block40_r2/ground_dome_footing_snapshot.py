"""구조물팀 PPT 최종 역T형 링 기초의 등가 수직하중에 대한 화성 지반 해석."""

import argparse
import hashlib
import json
import math
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import openseespy.opensees as ops

from ground_20261107 import require_finite, save_response_viewer, static_ramp
from ground_hybrid.config import default_config
from ground_hybrid.conforming_mesh import build_conforming_ground
from ground_hybrid.run_journal import RunJournal, atomic_json
from ground_hybrid.material_audit import collect_material_state

# 구조물팀 PPT 슬라이드 16–18의 최종안. 링 중심 반경은 기존 돔 좌표 기준이다.
RING_RADIUS_M = 20.0
RING_WIDTH_M = 14.0
EMBEDMENT_M = 3.0
FLANGE_THICKNESS_M = 2.5
COVER_THICKNESS_M = 6.2
DOME_MASS_T = 200.0
RING_MASS_T = 10_260.0
PLAN_AREA_M2 = 1200.0
NET_PRESSURE_PA = 100_600.0
DOWNWARD_FORCE_N = 150.3e6
SANDY_REGOLITH_MATERIALS = {
    1400: 'mars_insight_sandy_regolith',
    1300: 'mars_insight_sandy_regolith_1300',
}
SANDY_REGOLITH_CANDIDATES = {
    'mgs1': 'mars_insight_sand_mgs1_dp_1300',
    'mgs1s': 'mars_insight_sand_mgs1s_dp_1300',
    'jez1': 'mars_insight_sand_jez1_dp_1300',
}
NONLINEAR_PROFILES = Path(__file__).parent / 'ground_hybrid/configs/Mars_nonlinear_20261125_profiles.json'


def footing_areas(ground):
    areas = defaultdict(float)
    for face in ground['faces']:
        if face.get('foundation_selected'):
            for tag, area in zip(face['nodes'], face['weights_m2']):
                areas[tag] += area
    return dict(areas)


def settlement(tags, displacements, reference, areas):
    values = -(displacements - reference)[:, 2]
    return dict(
        average_mm=float(1000 * np.average(values, weights=[areas.get(tag, 0.0) for tag in tags])),
        maximum_mm=float(1000 * max(value for tag, value in zip(tags, values) if tag in areas)),
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--ground-width', type=float, default=90.0)
    cli.add_argument('--depth', type=float, default=30.0)
    cli.add_argument('--sectors', type=int, default=24)
    cli.add_argument('--radial-subdivisions', type=int, default=2)
    cli.add_argument('--static-steps', type=int, default=20)
    cli.add_argument('--footing-steps',type=int,default=8,help='load increments for loading and pressurization; material properties unchanged')
    cli.add_argument('--solver-policy',choices=['fallback','primary_only'],default='fallback',help='primary_only disables algorithm switching; failed increments are still adaptively reduced')
    cli.add_argument('--gravity-algorithm',choices=['KrylovNewton','NewtonLineSearch','Newton'],help='optional algorithm specifically for selfweight; footing uses static-algorithm')
    nonlinear_profiles = json.loads(NONLINEAR_PROFILES.read_text())
    cli.add_argument('--material-profile', choices=['nonlinear_20261125','recommended_20261124','reference',*nonlinear_profiles['profiles']], default='nonlinear_20261125',
                     help='all-elastoplastic selected profile (default), explicit candidate, or a previous reference')
    cli.add_argument(
        '--static-algorithm', choices=['KrylovNewton', 'NewtonLineSearch', 'Newton'],
        default='KrylovNewton',
    )
    sand_options = cli.add_mutually_exclusive_group()
    sand_options.add_argument(
        '--sandy-regolith-density', type=int, choices=SANDY_REGOLITH_MATERIALS,
        help='bulk density of the 0.2–3 m sandy regolith layer in kg/m³',
    )
    sand_options.add_argument(
        '--sandy-regolith-candidate', choices=SANDY_REGOLITH_CANDIDATES,
        help='Drucker–Prager material candidate for the 0.2–3 m layer',
    )
    cli.add_argument('--model-only', action='store_true', help='assemble the mesh without static analysis')
    cli.add_argument('--html', action='store_true')
    cli.add_argument('--output-dir', type=Path, required=True)
    args = cli.parse_args()
    outer = RING_RADIUS_M + RING_WIDTH_M / 2
    inner = RING_RADIUS_M - RING_WIDTH_M / 2
    if not math.isfinite(args.ground_width) or args.ground_width <= 2 * outer:
        cli.error('ground-width must exceed the 54 m footing diameter')
    if not math.isfinite(args.depth) or args.depth <= EMBEDMENT_M:
        cli.error('depth must extend below the footing base')
    if args.static_steps < 1 or args.footing_steps < 1:
        cli.error('static-steps and footing-steps must be positive')
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        cli.error('choose a fresh output directory')

    # 1. 데이터 불러오기
    config = default_config('Mars',mars_profile='reference')
    material_profile = dict(name=args.material_profile)
    if args.material_profile == 'nonlinear_20261125' or args.material_profile in nonlinear_profiles['profiles']:
        selected=nonlinear_profiles['selected_profile'] if args.material_profile=='nonlinear_20261125' else args.material_profile
        profile=nonlinear_profiles['profiles'][selected]
        config['layers']=profile['layers'];config['scenario']=selected
        material_profile.update(name=selected,requested_name=args.material_profile,selection_status=nonlinear_profiles['selection_status'],
                                profile_path=str(NONLINEAR_PROFILES.resolve()),profile_sha256=hashlib.sha256(NONLINEAR_PROFILES.read_bytes()).hexdigest(),
                                specimen_source=profile['specimen_source'],specimen_source_sha256=profile['specimen_source_sha256'],note=profile['note'])
    elif args.material_profile == 'recommended_20261124':
        profile_path = Path(__file__).parent / 'ground_hybrid/configs/Mars_recommended_20261124.json'
        profile = json.loads(profile_path.read_text())
        config['layers'] = profile['layers']
        config['scenario'] = args.material_profile
        material_profile.update(profile_path=str(profile_path.resolve()),profile_sha256=hashlib.sha256(profile_path.read_bytes()).hexdigest(),
                                specimen_source=profile['specimen_source'],specimen_source_sha256=profile['specimen_source_sha256'],
                                note=profile['note'])
    material = SANDY_REGOLITH_CANDIDATES.get(args.sandy_regolith_candidate)
    if args.sandy_regolith_density is not None:
        material = SANDY_REGOLITH_MATERIALS[args.sandy_regolith_density]
    if material is not None:
        material_profile['explicit_sand_override']=material
        sandy_materials = {*SANDY_REGOLITH_MATERIALS.values(), *SANDY_REGOLITH_CANDIDATES.values()}
        for layer in config['layers']:
            if layer['material'] in sandy_materials or layer.get('role') == 'sand':
                layer['material'] = material
                break
        else:
            cli.error('no InSight sandy regolith layer in Mars.json')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    config.update(
        width_m=args.ground_width, length_m=args.ground_width, depth_m=args.depth,
        max_vertical_spacing_m=4.0, max_aspect_ratio=200.0,
    )
    uplift = NET_PRESSURE_PA * PLAN_AREA_M2
    result = dict(
        complete=False,
        material_profile=material_profile,
        numerical_settings=dict(static_steps=args.static_steps,footing_steps=args.footing_steps,static_algorithm=args.static_algorithm,gravity_algorithm=args.gravity_algorithm or args.static_algorithm,solver_policy=args.solver_policy),
        ground_domain_m=[args.ground_width, args.ground_width, args.depth],
        footing=dict(
            centerline_radius_m=RING_RADIUS_M, outer_diameter_m=2 * outer,
            width_m=RING_WIDTH_M, embedment_m=EMBEDMENT_M, flange_thickness_m=FLANGE_THICKNESS_M,
            cover_thickness_m=COVER_THICKNESS_M, dome_mass_t=DOME_MASS_T, ring_mass_t=RING_MASS_T,
            circular_footprint_m2=math.pi * (outer**2 - inner**2),
        ),
        loads=dict(
            source='구조물팀 PPT slides 16–18, selected 14 m footing',
            downforce_N=DOWNWARD_FORCE_N, pressure_plan_area_m2=PLAN_AREA_M2,
            net_internal_pressure_Pa=NET_PRESSURE_PA, dome_uplift_N=uplift,
            pressurized_net_compression_N=DOWNWARD_FORCE_N - uplift,
        ),
        scope=(
            'Vertical settlement under area-distributed equivalent thrust only. '
            'Initial gravity is equilibrated on the excavated ground. '
            'Foundation/soil interaction, uplift resistance, structural design and soil failure are outside this analysis.'
        ),
    )
    summary = args.output_dir / 'summary.json'
    started = time.perf_counter()
    ops.wipe()
    ops.logFile(str(args.output_dir / 'opensees.log'))
    try:
        # 2. 모델링
        ground = build_conforming_ground(
            config, radius=RING_RADIUS_M, width=RING_WIDTH_M, kind='ring',
            sectors=args.sectors, radial_subdivisions=args.radial_subdivisions,
            excavation_depth_m=EMBEDMENT_M,
        )
        areas = footing_areas(ground)
        area = sum(areas.values())
        result['footing']['modeled_footprint_m2'] = area
        result['loads']['depressurized_contact_pressure_Pa'] = DOWNWARD_FORCE_N / area
        result['loads']['pressurized_contact_pressure_Pa'] = (DOWNWARD_FORCE_N - uplift) / area
        result['ground_mesh'] = ground['mesh']
        result['ground_elements'] = len(ground['elements'])
        result['ground_nodes'] = len(ground['nodes'])
        result['ground_mass_kg'] = ground['mass_kg']
        result['ground_materials'] = {
            name: ground['materials'][name] for name in ground['material_volume_fractions']
        }
        result['ground_layers'] = [
            {**layer, 'density_kg_m3': ground['materials'][layer['material']]['density_kg_m3']}
            for layer in config['layers']
            if layer['material'] in ground['material_volume_fractions']
        ]
        nodes = sorted(ground['nodes'])
        if args.model_only:
            result['scope'] = 'Ground mesh and material assignment at assembly'
            result['analysis_performed'] = False
            if args.html:
                coordinates = np.asarray([ground['nodes'][tag] for tag in nodes])
                result['viewer'] = save_response_viewer(
                    ground, nodes, coordinates, [np.zeros_like(coordinates)], [0.0],
                    args.output_dir / 'ground_geometry.html',
                    context=dict(
                        title=f'화성 링 기초 지반 · {(args.sandy_regolith_candidate or material_profile["name"]).upper()}',
                        subtitle=f'{args.ground_width:g} × {args.ground_width:g} m · 깊이 {args.depth:g} m',
                        summary=f'절점 {len(nodes):,}개 · bbarBrick 요소 {len(ground["elements"]):,}개',
                        static_only=True,
                        geometry_only=True,
                        scope='링 굴착 깊이 3 m · 바닥판 폭 14 m. 회전하거나 층을 선택해 메시를 확인할 수 있습니다.',
                    ),
                    scale=1.0,
                )
            result['complete'] = True
            atomic_json(summary, result)
            return

        atomic_json(summary, result)

        # 3. 해석
        journal = RunJournal(args.output_dir, dict(
            candidate=args.sandy_regolith_candidate, static_steps=args.static_steps,
            primary_algorithm=args.static_algorithm,
        ))
        result['analysis_performed'] = True
        algorithms = tuple(dict.fromkeys((args.static_algorithm, 'KrylovNewton', 'NewtonLineSearch', 'Newton')))
        if args.solver_policy=='primary_only':
            algorithms=(args.static_algorithm,)

        def track_attempt(event):
            journal(event)
            result['solver_progress'] = event
            result['elapsed_s'] = time.perf_counter() - started
            atomic_json(summary, result)

        def track_commit():
            displacement = require_finite([ops.nodeDisp(tag) for tag in nodes], 'committed ground state')
            result['last_committed_ground_displacement_m'] = float(np.linalg.norm(displacement, axis=1).max())
            atomic_json(summary, result)

        result['stage'] = 'gravity'
        gravity_algorithms=algorithms
        if args.gravity_algorithm:
            gravity_algorithms=(args.gravity_algorithm,) if args.solver_policy=='primary_only' else tuple(dict.fromkeys((args.gravity_algorithm,*algorithms)))
        result['gravity'] = static_ramp(
            ground, ground['gravity_loads'], steps=args.static_steps,
            algorithms=gravity_algorithms, attempt_observer=track_attempt, on_commit=track_commit,
        )
        gravity = require_finite([ops.nodeDisp(tag) for tag in nodes], 'excavated selfweight')
        result['gravity']['material_state'], gravity_plastic_norms = collect_material_state(ground,args.output_dir/'gravity_material_state.csv')
        atomic_json(summary,result)

        # 감압 상태의 전체 중량을 링 바닥에 전달한 다음 내압 순들림만큼 접촉 하중을 줄인다.
        compression = {tag: [0.0, 0.0, -DOWNWARD_FORCE_N * weight / area] for tag, weight in areas.items()}
        result['stage'] = 'depressurized'
        result['depressurized'] = static_ramp(
            ground, compression, pattern=800002, steps=args.footing_steps,
            algorithms=algorithms, attempt_observer=track_attempt, on_commit=track_commit,
        )
        loaded = require_finite([ops.nodeDisp(tag) for tag in nodes], 'depressurized settlement')
        result['depressurized'].update(settlement(nodes, loaded, gravity, areas))
        result['depressurized']['material_state'], _ = collect_material_state(ground,args.output_dir/'depressurized_material_state.csv',reference_norms=gravity_plastic_norms)
        atomic_json(summary,result)
        unloading = {tag: [0.0, 0.0, uplift * weight / area] for tag, weight in areas.items()}
        result['stage'] = 'pressurized'
        result['pressurized'] = static_ramp(
            ground, unloading, pattern=800003, steps=args.footing_steps,
            algorithms=algorithms, attempt_observer=track_attempt, on_commit=track_commit,
        )
        pressurized = require_finite([ops.nodeDisp(tag) for tag in nodes], 'pressurized settlement')
        result['pressurized'].update(settlement(nodes, pressurized, gravity, areas))
        result['pressurized']['material_state'], _ = collect_material_state(ground,args.output_dir/'pressurized_material_state.csv',reference_norms=gravity_plastic_norms)
        np.savez_compressed(args.output_dir/'stage_displacements.npz',node_tags=nodes,coordinates=np.asarray([ground['nodes'][tag] for tag in nodes]),
                            gravity_m=gravity,depressurized_m=loaded,pressurized_m=pressurized)

        # 4. 시각화 및 결과 저장
        if args.html:
            coordinates = [ground['nodes'][tag] for tag in nodes]
            result['viewer'] = save_response_viewer(
                ground, nodes, coordinates, [loaded - gravity], [0.0],
                args.output_dir / 'footing_settlement.html',
                context=dict(
                    title='화성 돔 링 기초 · 감압 상태의 침하',
                    subtitle=f'3 m 매입 · 바닥판 폭 14 m · {args.material_profile}',
                    reference_label='자중 평형 이후 증분 변위',
                    static_only=True,
                    scope='굴착 지반의 자중 평형을 기준으로 한 등가 링 하중의 증분 변위.',
                ),
                scale=1.0,
            )
            result['pressurized_viewer'] = save_response_viewer(
                ground,nodes,coordinates,[pressurized-gravity],[0.0],args.output_dir/'footing_settlement_pressurized.html',
                context=dict(title='화성 돔 링 기초 · 가압 후 침하',subtitle=f'3 m 매입 · 바닥판 폭 14 m · {args.material_profile}',
                             reference_label='자중 평형 이후 증분 변위',static_only=True,scope='감압 상태로 재하한 뒤 내압 등가 들림 하중을 적용한 순서의 결과.'),scale=1.0)
        result['complete'] = True
        result['stage'] = 'complete'
        result['elapsed_s'] = time.perf_counter() - started
        atomic_json(summary, result)
    except Exception as exc:
        result['error'] = str(exc)
        atomic_json(summary, result)
        raise
    finally:
        ops.wipe()


if __name__ == '__main__':
    main()
