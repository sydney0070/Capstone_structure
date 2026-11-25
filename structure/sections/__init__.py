"""Structural section builders using the shared material catalog."""
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
for path in (ROOT,ROOT/'material/material_release_20261101'):
    if str(path) not in sys.path:sys.path.insert(0,str(path))
