"""Individual material entry point; shared specimen definitions in experiments.py."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tune import main
if __name__ == '__main__':
    main(['run', '--case', 'pile_py', '--out', str(Path(__file__).resolve().parents[1]/'results'/'pile_py'), *sys.argv[1:]])
