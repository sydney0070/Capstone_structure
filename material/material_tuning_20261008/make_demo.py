"""Generate clearly labeled SYNTHETIC data to verify calibration plumbing."""
from pathlib import Path
from copy import deepcopy
from presets import get_case
from experiments import simulate
from tune import save,write_json,BASE

if __name__=='__main__':
    root=BASE/'synthetic_demo';root.mkdir(exist_ok=True)
    truth=get_case('rebar_steel02');truth['params'].update(Fy=380000.,b=.015)
    truth['status']='SYNTHETIC_TRUTH_NOT_EXPERIMENT'
    save(truth,simulate(truth),root/'train')
    validation=deepcopy(truth);validation['test']['amplitude']=.008
    save(validation,simulate(validation),root/'validation')
    guess=get_case('rebar_steel02');guess['bounds']=dict(Fy=[350000.,450000.],b=[.005,.025])
    write_json(root/'starting_config.json',guess)
    write_json(root/'targets.json',dict(case=truth['name'],units=truth['units'],synthetic=True,
        source='OpenSees generated synthetic truth: Fy=380000 kPa, b=0.015; NOT experimental data',
        conditions_confirmed=True,tests=[
            dict(csv='train/history.csv',role='train',test=truth['test'],channels={'y':1.}),
            dict(csv='validation/history.csv',role='validation',test=validation['test'],channels={'y':1.})]))
