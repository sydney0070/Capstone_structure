"""Durable attempt journal, including entry into a possibly blocking native call."""
from pathlib import Path
import json
import time


def atomic_json(path, value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False));tmp.replace(path)


class RunJournal:
    def __init__(self, directory, metadata):
        self.directory=Path(directory);self.directory.mkdir(parents=True,exist_ok=True)
        self.path=self.directory/'solver_attempts.jsonl'
        if self.path.exists():raise ValueError('fresh journal required')
        self.start=time.perf_counter();self.metadata=metadata
        self.native_log=self.directory/'opensees.log';self.native_offset=0
        atomic_json(self.directory/'model_checkpoint.json',dict(metadata,complete=False,safe_verdict=False))

    def __call__(self, event):
        row=dict(elapsed_s=time.perf_counter()-self.start,**event)
        memory_failure=False
        if event['event']=='attempt_started' and self.native_log.exists():self.native_offset=self.native_log.stat().st_size
        if event['event']=='attempt_finished' and event.get('code') and self.native_log.exists():
            with self.native_log.open() as stream:
                stream.seek(self.native_offset);message=stream.read(1024*1024)
            # UMFPACK 5.x status -1 is allocation failure. Reducing load cannot
            # be presented as convergence recovery from this resource error.
            memory_failure='numeric analysis returns -1 -- Umfpack' in message
            if memory_failure:row['failure_class']='linear_factorization_out_of_memory'
        with self.path.open('a') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
        atomic_json(self.directory/'solver_checkpoint.json',row)
        if memory_failure:raise MemoryError('UMFPACK numeric factorization reports out-of-memory (-1); inspect opensees.log; material strength and tolerances unchanged')
