"""Own and stop an entire analysis process group, including native workers.

The environment marker is scoped to the launched process tree. CLI entry points
use it to avoid spawning a second, detached supervisor inside an owned group.
"""
import math
import os
import signal
import subprocess
import time

MANAGED_ENV = 'CAPSTONE_ANALYSIS_PROCESS_GROUP'
THREAD_VARIABLES=('OPENBLAS_NUM_THREADS','OPENBLAS_DEFAULT_NUM_THREADS','OMP_NUM_THREADS',
                  'MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')


def _signal_group(job, signum):
    try:
        os.killpg(job.pid, signum)
    except ProcessLookupError:
        return False
    except PermissionError:
        # macOS can report EPERM when only adopted zombie descendants remain.
        # Do not hide a real permission failure involving a live process.
        snapshot=subprocess.run(['/bin/ps','-axo','pgid=,stat='],capture_output=True,text=True,timeout=2)
        if snapshot.returncode:raise
        members=[line.split()[1] for line in snapshot.stdout.splitlines()
                 if len(line.split())==2 and line.split()[0]==str(job.pid)]
        if any(not status.startswith('Z') for status in members):raise
        return False
    return True


def stop_group(job, grace_s=.5):
    """Signal descendants even if their direct parent has already exited."""
    _signal_group(job, signal.SIGTERM)
    deadline = time.monotonic() + grace_s
    while time.monotonic() < deadline:
        job.poll()  # reap the group leader if possible
        if not _signal_group(job, signal.SIGTERM):break
        time.sleep(.02)
    _signal_group(job, signal.SIGKILL)
    job.wait()


def run_isolated(command, timeout, *, stdout=None, stderr=None, join_managed_group=False):
    """Stop owned workers and classify timeout, interrupt or unverified cleanup."""
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('positive finite timeout required')
    if join_managed_group and os.environ.get(MANAGED_ENV)!='1':raise ValueError('joining requires an outer process-group owner')
    env = dict(os.environ, **{MANAGED_ENV: '1'})
    threads=int(env.get('CAPSTONE_SOLVER_THREADS','1'))
    if threads<1:raise ValueError('positive integer CAPSTONE_SOLVER_THREADS required')
    env.update({name:str(threads) for name in THREAD_VARIABLES})
    started = time.perf_counter()
    job = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                           start_new_session=not join_managed_group, env=env)
    cleanup_error=None
    try:
        code = job.wait(timeout=timeout)
        state = {0: 'complete', 124: 'timeout', 130: 'interrupted'}.get(code, 'failed')
    except (subprocess.TimeoutExpired, KeyboardInterrupt) as exc:
        state, code = ('timeout', 124) if isinstance(exc, subprocess.TimeoutExpired) else ('interrupted', 130)
    finally:
        # Also stop descendants of a normally exiting parent. A successful
        # launcher is not authorization to leave background solver processes.
        try:
            if join_managed_group:
                # Only the two known CLI --worker entry points use this mode.
                # They do not spawn worker processes. Their outer supervisor
                # owns the shared group and covers a simultaneous outer timeout.
                if job.poll() is None:
                    job.terminate()
                    try:job.wait(timeout=.5)
                    except subprocess.TimeoutExpired:job.kill();job.wait()
            else:stop_group(job)
        except OSError as exc:
            # Preserve timeout evidence even when a host sandbox denies process
            # inspection. Never label an unverified cleanup as successful.
            cleanup_error=str(exc)
            if code==0:state,code='cleanup_unverified',125
    return dict(status=state, returncode=code, complete=code == 0,
                elapsed_s=time.perf_counter()-started, process_group_cleanup=cleanup_error is None,
                cleanup_error=cleanup_error,solver_threads=threads,
                cleanup_scope='worker_in_outer_owned_group' if join_managed_group else 'owned_process_group')
