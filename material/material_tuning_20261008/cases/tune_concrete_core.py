"""Individual material entry point; shared specimen definitions in experiments.py."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tune import main
if __name__ == '__main__':
    main(['run', '--case', 'concrete_core', '--out', str(Path(__file__).resolve().parents[1]/'results'/'concrete_core'), *sys.argv[1:]])
