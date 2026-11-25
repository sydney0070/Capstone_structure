"""InSight 층서를 깊이 방향으로만 이산화한 1차원 지반 응답 모델."""

import numpy as np
import openseespy.opensees as ops

from ground_20261107 import build_ground_20261107, cli


def build_1d_column(config):
    """1 m² 단면에 요소 하나를 배치하고 같은 깊이의 절점 변위를 묶는다."""
    column_config = {
        **config,
        'width_m': 1.0,
        'length_m': 1.0,
        'core_width_m': 0.7,
        'core_length_m': 0.7,
    }
    east, north = np.meshgrid([-0.5, 0.5], [-0.5, 0.5])
    ground = build_ground_20261107(
        column_config, (east, north, np.zeros_like(east)), roller_sides=False,
    )

    # 바닥의 네 절점은 이미 고정되어 있다. 그 위 단면마다 세 병진 자유도를 공유한다.
    for level in range(ground['mesh']['nz']):
        master = ground['index'][0, 0, level]
        for i, j in ((1, 0), (0, 1), (1, 1)):
            ops.equalDOF(master, ground['index'][i, j, level], 1, 2, 3)
    ground['boundary'] = 'tied_column'
    return ground


if __name__ == '__main__':
    raise SystemExit(cli(column=True))
