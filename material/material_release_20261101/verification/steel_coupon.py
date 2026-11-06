"""Replay existing selected S355 coupon holdouts with SI registry and current runtime."""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.materials import register_uniaxial
SOURCE=ROOT.parent/'material_tuning_20261008'/'material_calibration'/'s355_hem320_web'
OUT=ROOT/'results'/'steel_coupon';OUT.mkdir(exist_ok=True)
def response(strains,substeps):
 ops.wipe();register_uniaxial(ops,'steel_s355_coupon',1);ops.testUniaxialMaterial(1);y=[];last=0.
 for x in strains:
  for v in np.linspace(last,x,substeps+1)[1:]:ops.setStrain(float(v))
  y.append(ops.getStress());last=x
 return np.array(y)
result=dict(runtime=ops.version(),units='Pa, engineering strain',scope='selected same-source coupon histories ONLY, independent specimen IDs preserved, no refit',tests=[])
for f in ['train_LP4_id0.csv','train_LP6_id8.csv','validation_LP2_id3.csv','validation_LP9_id6.csv']:
 a=np.loadtxt(SOURCE/f,delimiter=',',skiprows=1);x=a[:,0];target=a[:,1]*1000;y=response(x,1);refined=response(x,4)
 nrmse=float(np.sqrt(np.mean((y-target)**2))/np.ptp(target));peak=float(max(abs(np.max(y)-np.max(target))/abs(np.max(target)),abs(np.min(y)-np.min(target))/abs(np.min(target))));ref=float(np.max(abs(y-refined))/max(np.max(abs(y)),1))
 r=dict(file=f,role='validation' if f.startswith('validation') else 'train',source_sha256=hashlib.sha256((SOURCE/f).read_bytes()).hexdigest(),points=len(x),range_nrmse=nrmse,peak_relative_error=peak,increment_refinement_relative_error=ref,**{'pass':bool(nrmse<.05 and peak<.10 and ref<.005)})
 result['tests'].append(r)
 with (OUT/f).open('w',newline='') as fp:w=csv.writer(fp);w.writerow(['strain','target_Pa','model_Pa','refined_Pa']);w.writerows(zip(x,target,y,refined))
 print(r)
result['all_pass']=all(r['pass'] for r in result['tests']);(OUT/'summary.json').write_text(json.dumps(result,indent=2))
