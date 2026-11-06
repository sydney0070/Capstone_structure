import json,sys,traceback
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.ground import build_ground
X,Y=np.meshgrid(np.linspace(0,4,3),np.linspace(0,4,3));Z=np.zeros_like(X)
result=dict(runtime=ops.version(),cases=[])
for body in ['Earth','Moon','Mars']:
 for factor in [.5,1.,2.]:
  ops.wipe()
  try:
   d=build_ground(X,Y,Z,body=body,depth_m=3.,n_depth=3,stiffness_factor=factor)
   r=dict(body=body,stiffness_factor=factor,release=d['release'],gravity=d['gravity_summary'],mass_kg=d['total_soil_mass_kg'],weight_N=d['expected_total_weight_n'],nodes=len(d['node_coords']),elements=len(d['elements']),validation=d['validation_report'],pass_gravity=d['gravity_summary']['relative_balance_error']<1e-4)
  except Exception as e:r=dict(body=body,stiffness_factor=factor,pass_gravity=False,error=str(e),traceback=traceback.format_exc())
  result['cases'].append(r);(ROOT/'results'/'ground_adapter.json').write_text(json.dumps(result,indent=2))
  print('RELEASE_CASE',body,factor,r['pass_gravity'],r.get('error'),flush=True)
