"""Compatibility only. Canonical implementation: integrated.reference_bounded."""
from pathlib import Path
import importlib, sys
ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
implementation = importlib.import_module('integrated.reference_bounded')
if __name__ == '__main__': raise SystemExit(implementation.main())
sys.modules[__name__] = implementation
