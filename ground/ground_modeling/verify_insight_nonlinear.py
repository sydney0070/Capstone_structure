"""Six independent analytical DP yield checks through the InSight material adapter."""
import importlib.util,json,sys
from pathlib import Path
b=Path(__file__).resolve().parent;sys.path.insert(0,str(b))
import ground_modeling_insight_updated as model
p=b.parents[1]/'material/material_release_20261101/verification/soil_material.py'
spec=importlib.util.spec_from_file_location('soil_check',p);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
v.OUT=b/'insight_verification'/'material';v.OUT.mkdir(parents=True,exist_ok=True)
from ground_hybrid.materials import record
v.material=record
def register(ops,name,tag):
 m=dict(model.MARS_SOIL_MATERIALS[name],tag=tag)
 model.register_opensees_soil_materials({name:m},[name],3)
v.register_nd=register
rows=[v.triaxial(n,p) for n in model.MARS_SOIL_ORDER[:2] for p in [1000.,15000.,60000.]]
(v.OUT/'summary.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2));assert all(r['pass'] for r in rows)
assert model.parse_cli_args(['--insight-local']).dynamic_material_model=='nonlinear'
