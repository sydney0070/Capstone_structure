"""Run this independent material specimen; SI units."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from specimens import main
if __name__ == "__main__":
    main(["--case", 'hysteretic_crust_scale_benchmark'])
