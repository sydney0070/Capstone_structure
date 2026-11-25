"""Ground builders and solvers; legacy SSI names are resolved lazily."""
from pathlib import Path
import sys
_project_root = str(Path(__file__).resolve().parents[3])
if _project_root not in sys.path: sys.path.insert(0, _project_root)
from .config import default_config
from .mesh import build_ground
from .solver import initialize_gravity, apply_static_loads

def __getattr__(name):
    if name in ('add_raft','connect_matching_nodes'):
        from SSI import interface
        return getattr(interface,name)
    raise AttributeError(name)
