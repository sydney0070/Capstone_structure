"""Canonical soil–foundation/structure interface implementations (SI)."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
for directory in (ROOT, ROOT/'material/material_release_20261101', ROOT/'ground/ground_modeling'):
    value = str(directory)
    if value not in sys.path:
        sys.path.insert(0, value)
