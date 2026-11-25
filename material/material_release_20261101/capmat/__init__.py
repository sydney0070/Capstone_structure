"""Capstone SI material interface. Registration never wipes the caller's domain."""
from pathlib import Path
import sys
_project_root = str(Path(__file__).resolve().parents[3])
if _project_root not in sys.path: sys.path.insert(0, _project_root)
from .materials import catalog,material,register_uniaxial,register_nd,register_plate
from .mechanics import shielding,contact_parameters,radiation_dashpots,laminate_abd,ply_screen
