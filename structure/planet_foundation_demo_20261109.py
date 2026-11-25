"""Compatibility only. Canonical implementation: integrated.planet_foundation_demo_20261109."""
from pathlib import Path
import importlib, sys
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
implementation = importlib.import_module('integrated.planet_foundation_demo_20261109')
if __name__ == '__main__': raise SystemExit(implementation.supervised_cli())
sys.modules[__name__] = implementation
