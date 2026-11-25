"""Compatibility only. Canonical implementation: integrated.ground_reference_demo."""
from pathlib import Path
import importlib, sys
ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
implementation = importlib.import_module('integrated.ground_reference_demo')
if __name__ == '__main__': raise SystemExit(implementation.main())
sys.modules[__name__] = implementation
