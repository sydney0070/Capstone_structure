from pathlib import Path
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PROJECT_ROOT), str(PROJECT_ROOT/'structure'),
               str(PROJECT_ROOT/'ground/ground_modeling'),
               str(PROJECT_ROOT/'material/material_release_20261101')]

"""Run an isolated CLI analysis with a hard wall-clock limit.

Example: python -m ground_hybrid.bounded_run --timeout-s 180 --body Moon
Uses a child process because an OpenSees native call may not handle a Python
timer until after returning. Never used inside an existing structural domain.
"""
import argparse
import math
import subprocess
import sys
from ground_hybrid.process_control import run_isolated


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__,add_help=False)
    p.add_argument('--timeout-s',type=float,default=180.)
    args,remaining=p.parse_known_args(argv)
    if not math.isfinite(args.timeout_s) or args.timeout_s<=0:
        p.error('positive finite --timeout-s required')
    result=run_isolated([sys.executable,str(PROJECT_ROOT/'integrated/ground_reference_demo.py'),*remaining],args.timeout_s)
    if result['status']=='timeout':
        print(f'FAILED: analysis exceeded {args.timeout_s:g}s; child stopped. '
              'Partial/older output is not a completed result.',file=sys.stderr)
    return result['returncode']


if __name__=='__main__':raise SystemExit(main())
