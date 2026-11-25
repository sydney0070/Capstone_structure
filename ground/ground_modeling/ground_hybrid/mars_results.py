"""Ground response histories, diagnostic snapshots and HTML visualization."""

import csv
import math
from pathlib import Path

import numpy as np
import openseespy.opensees as ops


class GroundHistory:
    """성공한 단계마다 bbarBrick의 8개 적분점에서 소성 변형률과 일을 기록한다."""

    def __init__(self, g, path=None):
        from ground_modeling import brick_jacobian_determinants

        self.path = Path(path) if path else None
        if self.path and self.path.exists():
            raise ValueError('fresh history output required')
        self.previous = None
        self.reference_elastic = None
        self.reference_plastic = None
        self.max_incremental_plastic = 0.0
        self.internal_work = 0.0
        self.plastic_work = 0.0
        self.rows = []
        self.max_plastic = 0.0
        self.max_abs_strain = 0.0

        self.layout = []
        compliance = []
        weights = []
        for element, record in g['elements'].items():
            coordinates = np.array([g['nodes'][tag] for tag in record['nodes']])
            jacobians = brick_jacobian_determinants(coordinates)
            if not np.isfinite(jacobians).all() or np.min(jacobians) <= 0:
                raise ValueError('positive finite Jacobians required for history')
            if not math.isclose(sum(jacobians), record['volume_m3'], rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError('history/mesh quadrature volume mismatch')

            material = g['materials'][record['material']]
            shear = material['G_Pa']
            bulk = material['K_Pa']
            lame = bulk - 2 * shear / 3
            stiffness = np.zeros((6, 6))
            stiffness[:3, :3] = lame
            stiffness[np.arange(3), np.arange(3)] += 2 * shear
            stiffness[3:, 3:] = np.eye(3) * shear
            elastic_compliance = np.linalg.inv(stiffness)

            # OpenSees의 적분점 순서와 detJ 배열 순서를 맞춰 체적 가중치를 배정한다.
            for ip, weight in enumerate(jacobians, start=1):
                self.layout.append((element, ip))
                compliance.append(elastic_compliance)
                weights.append(weight)
        self.compliance = np.array(compliance)
        self.weights = np.array(weights)

    def capture(self):
        stresses = []
        strains = []
        for element, ip in self.layout:
            stress = ops.eleResponse(element, 'material', ip, 'stress')
            strain = ops.eleResponse(element, 'material', ip, 'strain')
            if len(stress) != 6 or len(strain) != 6:
                raise RuntimeError(f'Missing all-GP response e={element} ip={ip}')
            stresses.append(stress)
            strains.append(strain)
        stress = np.asarray(stresses)
        strain = np.asarray(strains)
        if not np.isfinite(stress).all() or not np.isfinite(strain).all():
            raise RuntimeError('Nonfinite all-GP history')

        # 일정 탄성강성으로 탄성 변형률을 복원한다. 전단 성분은 공학전단 변형률이다.
        elastic_strain = np.einsum('nij,nj->ni', self.compliance, stress)
        plastic_strain = strain - elastic_strain
        elastic_energy = float(
            np.sum(0.5 * np.sum(stress * elastic_strain, axis=1) * self.weights)
        )

        # 연속 상태 사이의 평균 응력으로 일 증분을 적분한다.
        if self.previous is not None:
            previous_stress, previous_strain, previous_plastic = self.previous
            average_stress = (stress + previous_stress) / 2
            self.internal_work += float(np.sum(
                np.sum(average_stress * (strain - previous_strain), axis=1) * self.weights
            ))
            self.plastic_work += float(np.sum(
                np.sum(average_stress * (plastic_strain - previous_plastic), axis=1) * self.weights
            ))
        else:
            self.reference_elastic = elastic_energy
            self.reference_plastic = plastic_strain.copy()
        self.previous = (stress, strain, plastic_strain)

        plastic_norm = np.sqrt(
            np.sum(plastic_strain[:, :3] ** 2, axis=1)
            + 0.5 * np.sum(plastic_strain[:, 3:] ** 2, axis=1)
        )
        increment = plastic_strain - self.reference_plastic
        increment_norm = np.sqrt(
            np.sum(increment[:, :3] ** 2, axis=1) + 0.5 * np.sum(increment[:, 3:] ** 2, axis=1)
        )
        peak = int(np.argmax(plastic_norm))
        increment_peak = int(np.argmax(increment_norm))
        self.max_incremental_plastic = max(
            self.max_incremental_plastic, float(increment_norm[increment_peak])
        )
        self.max_plastic = max(self.max_plastic, float(plastic_norm[peak]))
        self.max_abs_strain = max(self.max_abs_strain, float(np.max(abs(strain))))

        plastic_above_threshold = plastic_norm > 1e-8
        increment_above_threshold = increment_norm > 1e-8
        elastic_change = elastic_energy - self.reference_elastic
        row = dict(
            max_plastic_strain_increment_norm=float(increment_norm[increment_peak]),
            incremental_plastic_GP_count_above_1e_8=int(sum(increment_above_threshold)),
            element_at_max_increment=self.layout[increment_peak][0],
            ip_at_max_increment=self.layout[increment_peak][1],
            time=float(ops.getTime()),
            gauss_points=len(self.layout),
            max_plastic_strain_norm=float(plastic_norm[peak]),
            plastic_GP_count_above_1e_8=int(sum(plastic_above_threshold)),
            element_at_max=self.layout[peak][0],
            ip_at_max=self.layout[peak][1],
            abs_strain_max=float(np.max(abs(strain))),
            internal_work_incremental_J=self.internal_work,
            plastic_work_incremental_J=self.plastic_work,
            elastic_energy_change_J=elastic_change,
            work_decomposition_residual_J=self.internal_work - self.plastic_work - elastic_change,
            ground_volume_m3=float(sum(self.weights)),
            plastic_volume_above_1e_8_m3=float(self.weights[plastic_above_threshold].sum()),
            incremental_plastic_volume_above_1e_8_m3=float(
                self.weights[increment_above_threshold].sum()
            ),
            integrated_plastic_strain_norm_m3=float(plastic_norm @ self.weights),
            integrated_incremental_plastic_strain_norm_m3=float(increment_norm @ self.weights),
            volume_mean_incremental_plastic_strain_norm=float(
                increment_norm @ self.weights / sum(self.weights)
            ),
        )
        if self.path:
            with self.path.open('a', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(row))
                if not self.rows:
                    writer.writeheader()
                writer.writerow(row)
        self.rows.append(row)
        return row

    def save(self, path):
        if self.rows:
            with Path(path).open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]))
                writer.writeheader()
                writer.writerows(self.rows)
        return dict(
            committed_samples=len(self.rows),
            gauss_points=len(self.layout),
            max_plastic_strain_norm=self.max_plastic,
            max_incremental_plastic_strain_norm=self.max_incremental_plastic,
            max_abs_strain=self.max_abs_strain,
            final=self.rows[-1] if self.rows else None,
            scope=(
                'all 8 GP of positive-Jacobian trilinear bbarBrick; detJ weights; '
                'e_plastic=e-C^-1 sigma for constant elastic moduli; incremental trapezoidal work, '
                'not full system energy or damage acceptance'
            ),
        )


def ground_state(g, *, gauss_points=1):
    """현재 시각의 응력·변형률과 Drucker–Prager 전단 항복함수를 기록한다."""
    from capmat.materials import dp_constants

    if gauss_points not in (1, 8):
        raise ValueError('gauss_points must be 1 or 8')

    stress_max = 0.0
    strain_max = 0.0
    yield_max = -float('inf')
    yield_element = None
    yield_ip = None
    low_pressure = {}
    constants = {
        name: dp_constants(material['phi_deg'], material['c_Pa'], material['psi_deg'])
        for name, material in g['materials'].items()
        if material['active_law'] == 'DruckerPrager'
    }
    for element, record in g['elements'].items():
        name = record['material']
        material = g['materials'][name]
        for ip in range(1, gauss_points + 1):
            stress = np.asarray(ops.eleResponse(element, 'material', ip, 'stress'), float)
            strain = np.asarray(ops.eleResponse(element, 'material', ip, 'strain'), float)
            if len(stress) != 6 or len(strain) != 6 or not np.isfinite(np.r_[stress, strain]).all():
                raise RuntimeError('Ground stress/strain response unavailable')

            stress_max = max(stress_max, float(max(abs(stress))))
            strain_max = max(strain_max, float(max(abs(strain))))
            if material['active_law'] == 'DruckerPrager':
                params = constants[name]
                trace = float(sum(stress[:3]))
                deviatoric_normal = stress[:3] - trace / 3
                shear = stress[3:6]
                deviatoric_norm = math.sqrt(float(
                    deviatoric_normal @ deviatoric_normal + 2 * (shear @ shear)
                ))
                yield_value = deviatoric_norm + params['rho'] * trace - math.sqrt(2 / 3) * params['sigmaY']
                if yield_value > yield_max:
                    yield_max = yield_value
                    yield_element = element
                    yield_ip = ip

            if name.startswith('moon_jsc1a_'):
                # 삼축시험의 최소 구속압 15 kPa를 평균 압축응력으로 간접 점검한다.
                row = low_pressure.setdefault(
                    name,
                    dict(
                        sampled_GP=0,
                        mean_compression_below_15kPa_GP=0,
                        minimum_mean_compression_Pa=float('inf'),
                        maximum_mean_compression_Pa=-float('inf'),
                    ),
                )
                mean_compression = -float(sum(stress[:3])) / 3
                row['sampled_GP'] += 1
                row['mean_compression_below_15kPa_GP'] += int(mean_compression < 15000)
                row['minimum_mean_compression_Pa'] = min(
                    row['minimum_mean_compression_Pa'], mean_compression
                )
                row['maximum_mean_compression_Pa'] = max(
                    row['maximum_mean_compression_Pa'], mean_compression
                )
    return dict(
        sampled_elements=len(g['elements']),
        integration_point=1 if gauss_points == 1 else '1..8',
        sampled_gauss_points=len(g['elements']) * gauss_points,
        abs_stress_component_max_Pa=stress_max,
        abs_strain_component_max=strain_max,
        DP_max_yield_function_Pa=yield_max if yield_element is not None else None,
        DP_element_at_max=yield_element,
        DP_ip_at_max=yield_ip,
        lunar_low_pressure_screen=low_pressure,
        scope=(
            f'{gauss_points} Gauss point(s) per element at snapshot time. '
            'DP shear surface diagnostic only; tensile surface not checked. '
            'Mean pressure screen is not triaxial confinement validation. '
            'Not time extrema, plastic work, bearing capacity or failure verdict.'
        ),
    )


def remove_base_rigid_motion(reference_base, current_base, reference_points, current_points):
    """SSI 모델에서 기초의 강체 병진과 회전을 빼고 잔여 변위를 구한다."""
    reference = np.asarray(reference_base)
    current = np.asarray(current_base)
    reference_center = reference.mean(axis=0)
    current_center = current.mean(axis=0)
    reference_centered = reference - reference_center
    current_centered = current - current_center
    u, _, vt = np.linalg.svd(reference_centered.T @ current_centered)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        # 반사가 포함되지 않도록 마지막 축을 뒤집어 순수 회전만 사용한다.
        u[:, -1] *= -1
        rotation = u @ vt
    rigid_points = (np.asarray(reference_points) - reference_center) @ rotation + current_center
    return np.asarray(current_points) - rigid_points


def require_finite(values, context):
    """최대값 등을 계산하기 전에 개별 응답의 NaN과 무한대를 확인한다."""
    response = np.asarray(values, dtype=float)
    if not np.isfinite(response).all():
        raise RuntimeError(f'Nonfinite response: {context}')
    return response


GROUND_VIEWER_COLORS = {
    'mars_insight_duricrust_proxy': '#E6A23C',
    'mars_insight_sandy_regolith': '#63B5D6',
    'mars_insight_sandy_regolith_1300': '#63B5D6',
    'mars_insight_sand_mgs1_dp_1300': '#63B5D6',
    'mars_insight_sand_mgs1s_dp_1300': '#63B5D6',
    'mars_insight_sand_jez1_dp_1300': '#63B5D6',
    'mars_insight_blocky_ejecta': '#9B6451',
    'mars_insight_fractured_basalt': '#596579',
    'mars_jez1_peak': '#63B5D6',
    'basalt_elastic_analogue': '#596579',
}


def save_response_viewer(
    g, tags, coordinates, displacements, times, path, *, context, structures=(), scale=None,
):
    """지반 및 선택적인 구조물 선분과 실제 계산 변위를 HTML로 저장한다."""
    from insight_dynamic_viewer import save_dynamic_viewer

    tags = [int(tag) for tag in tags]
    xyz = np.asarray(coordinates, dtype=float)
    displacement = np.asarray(displacements, dtype=float)
    times = np.asarray(times, dtype=float)
    if (
        len(set(tags)) != len(tags)
        or xyz.shape != (len(tags), 3)
        or displacement.shape != (len(times), len(tags), 3)
        or not len(times)
        or not all(np.isfinite(values).all() for values in (xyz, displacement, times))
        or np.any(np.diff(times) <= 0)
    ):
        raise ValueError('invalid aligned physical viewer frames')
    if not set(g['nodes']) <= set(tags):
        raise ValueError('record all physical ground nodes for the ground maximum')

    # HTML 출력의 절점 정렬 순서와 좌표·변위·구조물 선분의 인덱스를 맞춘다.
    order = np.argsort(tags)
    tags = [tags[index] for index in order]
    xyz = xyz[order]
    displacement = displacement[:, order, :]
    index_by_tag = {tag: index for index, tag in enumerate(tags)}
    lines = [
        dict(
            name=item['name'],
            color=item['color'],
            segments=[[index_by_tag[a], index_by_tag[b]] for a, b in item['segments']],
        )
        for item in structures
    ]

    names = list(dict.fromkeys(element['material'] for element in g['elements'].values()))
    palette = ['#c8754a', '#a65132', '#7c4633', '#4b3028']
    materials = {}
    for index, name in enumerate(names):
        material = g['materials'][name]
        materials[name] = dict(
            E=material['E_Pa'],
            rho=material['density_kg_m3'],
            law=material['active_law'],
            color=GROUND_VIEWER_COLORS.get(name, palette[index % len(palette)]),
        )

    if scale is None:
        # 배율은 표시용이다. 아래 프레임에는 원래 변위(m)를 전달한다.
        delta = displacement if context.get('static_only') else displacement - displacement[0]
        peak = float(np.linalg.norm(delta, axis=2).max())
        span = max(float(np.ptp(xyz, axis=0).max()), 1.0)
        scale = span * 0.06 / peak if peak > 1e-15 else 1.0
    model = dict(
        node_coords={tag: point.tolist() for tag, point in zip(tags, xyz)},
        elements=g['elements'],
        materials=materials,
        material_order=names,
        body_name=g['config']['body'],
        dynamic_material_model=' · '.join(dict.fromkeys(m['law'] for m in materials.values())),
        viewer_context=context,
        viewer_structures=lines,
        save_viewer_payload=True,
    )
    frames = [dict(time=float(time), disp=values.tolist()) for time, values in zip(times, displacement)]
    save_dynamic_viewer(model, frames, [], tags[0], path, visual_scale=scale)
    return dict(
        file=str(Path(path)),
        frames=len(times),
        nodes=len(tags),
        structures=len(lines),
        offline=True,
        reference='total' if context.get('static_only') else 'dynamic',
    )


def write_ground_viewer(ground, options, output_dir, node_tags, coordinates, frames, times, site):
    """계산이 끝난 절점 변위와 프레임 시각을 HTML 뷰어에 전달한다."""
    if options.static_only:
        title = '화성 지반 · 자중 응답'
        analysis_note = '자중 정적 응답.'
    else:
        title = '화성 지반 · 동적 응답'
        boundary_note = {
            'fixed_base_viscous_sides': '근사 점성 측면.',
            'static_box': '측면 롤러.',
            'tied_column': '깊이별 동일 수평면 변위.',
        }[ground['boundary']]
        analysis_note = f'SEIS 지표 관측 가속도를 기저 운동으로 가정. 고정 바닥·{boundary_note}'
    model_name = '1차원 지반 기둥' if ground['boundary'] == 'tied_column' else '3차원 지반'
    context = dict(
        title=title,
        subtitle=f'InSight / ELYSE · {model_name}',
        summary=f"지반 {len(ground['elements']):,} 요소 · 구조물·기초·연결부 없음",
        static_only=options.static_only,
        scope=(
            f"ELYSE {site['latitude_deg']}°N, {site['longitude_deg']}°E 기준 평탄 지반. "
            f"문헌 기반 기준 층서. 실제 DEM 요철 미적용. {analysis_note}"
        ),
    )
    return save_response_viewer(
        ground, node_tags, coordinates, frames, times,
        output_dir / 'ground_response.html', context=context, scale=options.visual_scale,
    )
