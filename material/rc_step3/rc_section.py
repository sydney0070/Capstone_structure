"""Step 3 RC fiber section, 2D prototype (SI: N, m, Pa).
All dimensions/reinforcement/confinement assumptions are PROVISIONAL.
"""
from dataclasses import dataclass, asdict
from rc_materials import register_rc_materials

@dataclass(frozen=True)
class SectionGeometry:
    width: float = 0.40              # m, assumed
    depth: float = 0.60              # m, assumed
    cover: float = 0.045             # m, assumed distance to outer rebar centroid is cover+bar dia/2
    bar_diameter: float = 0.020      # m, assumed
    bars_top: int = 3                # assumed
    bars_bottom: int = 3             # assumed
    ny: int = 24
    nz: int = 12

    def validate(self):
        assert self.width > 0 and self.depth > 0
        assert 0 < self.cover + self.bar_diameter/2 < min(self.width,self.depth)/2
        assert self.bars_top >= 2 and self.bars_bottom >= 2
        assert self.ny >= 2 and self.nz >= 2


def create_rc_fiber_section(ops, tag=10, geom=None):
    """Concrete02 throughout: NO assumed confinement enhancement. Steel02 layers.
    2D Fiber section y = cross-section depth axis, rebar at top/bottom.
    """
    import math
    geom = geom or SectionGeometry()
    geom.validate()
    register_rc_materials(ops, concrete_tag=1, steel_tag=2)
    y0, y1 = -geom.depth/2, geom.depth/2
    z0, z1 = -geom.width/2, geom.width/2
    ops.section('Fiber', tag)
    ops.patch('rect', 1, geom.ny, geom.nz, y0, z0, y1, z1)
    yc = geom.depth/2 - geom.cover - geom.bar_diameter/2
    zc = geom.width/2 - geom.cover - geom.bar_diameter/2
    as_bar = math.pi*geom.bar_diameter**2/4
    ops.layer('straight', 2, geom.bars_top, as_bar, yc, -zc, yc, zc)
    ops.layer('straight', 2, geom.bars_bottom, as_bar, -yc, -zc, -yc, zc)
    return geom
