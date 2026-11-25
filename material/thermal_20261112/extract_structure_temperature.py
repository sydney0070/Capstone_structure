"""Execute only existing environment functions, without running model/search.

Requires pandas/numpy in the bundled artifact Python runtime. Does not import
the structure module, change its globals, or execute OpenSees analysis.
"""
import ast,csv,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]
SOURCE=PROJECT/'structure/dome_in_mars_tuned.py'
DATA=PROJECT/'ground/ground_modeling/MARS_DATA'
NAMES={'_parse_pds_utc','_lmst_to_hours','load_twins','load_ps','_wind_sector','build_env_timeline','effective_env'}
CONSTANTS={'TWINS_TEMP_VALID_K','PS_VALID_PA','ENV_DT_S','MARS_GAS_CONSTANT','THERMAL_REF_TEMP_K','WIND_DIR_IS_FROM','WIND_SECTORS','INTERNAL_PRESSURE_PA','USE_INTERNAL_PRESSURE','USE_WIND','DIURNAL_SKIN_DEPTH_M','TWINS_FILES','PS_FILES'}
tree=ast.parse(SOURCE.read_text())
scope={'np':np,'pd':pd,'math':math,'datetime':datetime,'timezone':timezone}
for node in tree.body:
    if isinstance(node,ast.Assign):
        names=[n.id for n in node.targets if isinstance(n,ast.Name)]
        if len(names)==1 and names[0] in CONSTANTS:scope[names[0]]=ast.literal_eval(node.value)
functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in NAMES]
if len(functions)!=len(NAMES):raise RuntimeError('Source environment API changed; re-inspect before extraction')
exec(compile(ast.Module(body=functions,type_ignores=[]),str(SOURCE),'exec'),scope)
twins_files=[DATA/f for f in scope['TWINS_FILES']];ps_files=[DATA/f for f in scope['PS_FILES']]
twins,meta=scope['load_twins'](twins_files)
env=scope['build_env_timeline'](twins,meta,scope['load_ps'](ps_files))
profile=HERE/'structure_environment_temperature.csv'
with profile.open('w',newline='') as h:
    w=csv.writer(h);w.writerow(['time_s','temperature_K','delta_temperature_K'])
    w.writerows(zip(env['t_s'],env['T_K'],env['dT']))
covers=[]
for thickness in (0.,2.,7.5):
    eff=scope['effective_env'](env,SimpleNamespace(cover_t=thickness))
    effectiveT=env['T_ref']+eff['dT']
    covers.append(dict(cover_thickness_m=thickness,amplitude_factor=eff['thermal_factor'],
        minimum_structural_temperature_K=float(min(effectiveT)),maximum_structural_temperature_K=float(max(effectiveT)),
        min_delta_K=float(min(eff['dT'])),max_delta_K=float(max(eff['dT']))))
doc=dict(reference_temperature_K=float(env['T_ref']),minimum_K=float(min(env['T_K'])),maximum_K=float(max(env['T_K'])),
    temperature_assignment_status='EXTRACTED_EXISTING_STRUCTURE_ENVIRONMENT',
    reference_basis='existing THERMAL_REF_TEMP_K=None: mean of PS-overlap,600-second-binned,interpolated TWINS tip proxy temperatures',
    source_sample_count=meta['n_samples'],timeline_sample_count=len(env['T_K']),timeline_dt_s=float(env['dt']),
    valid_sensor_range_K=list(scope['TWINS_TEMP_VALID_K']),profile_csv=str(profile),
    source_structure_path=str(SOURCE),source_structure_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    extraction_method='AST-selected exact existing load_twins/load_ps/build_env_timeline/effective_env; no OpenSees or optimization execution',
    constants={name:scope[name] for name in sorted(CONSTANTS)},cover_cases=covers,
    component_scope={'dome':'temperature change is environmental dT times exp(-cover_thickness/0.05m)',
        'buried_ring':'legacy thermal_node_forces only acts on DOME_EDGES; docstring assumes ring deltaT=0',
        'soil':'no depth-dependent temperature assignment or conduction model in current integrated path',
        'current_tuned_main':'does not call the retained thermal environment helper; full SSI has no thermal stage'},
    thermal_transfer_scope='existing structure environment proxy; covered dome attenuation is a simplified assumption, not measured component or soil temperatures',
    sources=[dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in [*twins_files,*ps_files]])
(HERE/'structure_environment.json').write_text(json.dumps(doc,ensure_ascii=False,indent=2))
print(json.dumps({k:doc[k] for k in ['minimum_K','maximum_K','reference_temperature_K','timeline_sample_count','cover_cases','component_scope']},ensure_ascii=False,indent=2))
