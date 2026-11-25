"""Compatibility only. Canonical implementation: integrated.mars_dome_ring_ground."""
from pathlib import Path
import importlib, sys
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
implementation = importlib.import_module('integrated.mars_dome_ring_ground')
if __name__ == '__main__': raise SystemExit(implementation.main())
sys.modules[__name__] = implementation
