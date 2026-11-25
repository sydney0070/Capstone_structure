"""Compatibility only. Canonical implementation: SSI.contact."""
from pathlib import Path
import importlib, sys
ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
implementation = importlib.import_module('SSI.contact')
sys.modules[__name__] = implementation
