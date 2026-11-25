"""UCSD Grade 60, A6 dataset group. Source order and full supplied histories preserved."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
B=Path(__file__).resolve().parent
sys.path.insert(0,str(B.parent))
from presets import get_case
from tune import write_json
root=B/'rebar_ucsd_a6';root.mkdir(exist_ok=True)
c=get_case('rebar_steel02');c['name']='rebar_ucsd_grade60_a6'
c['params'].update(E=208e6,Fy=480e3,b=.018,R0=20.,cR1=.87,cR2=.06,a1=.02,a3=.025)
c['bounds']=dict(E=[150e6,250e6],Fy=[300e3,700e3],b=[.0001,.1],cR1=[.7,.98],cR2=[.01,.2],a1=[0.,.1],a3=[0.,.1])
c['test']['substeps']=1;c['status']='EXPERIMENTAL_CALIBRATION_START'
c['calibration_scope']='UCSD Grade 60 A6 group supplied histories; ASTM subclass mapping not asserted; not buckling, fracture or all Grade 60 products.'
tests=[];provenance=[]
for i in range(1,7):
 name=f'{i:02d}_A6{i+3}';role='train' if i<=4 else 'validation'
 xp=B/'raw/rebar_ucsd'/f'{name}_input_data.txt';yp=B/'raw/rebar_ucsd'/f'{name}_output_data.txt'
 x=np.loadtxt(xp).ravel();y=np.loadtxt(yp).ravel()*1000
 assert len(x)==len(y) and abs(x[0])<1e-4 and abs(y[0])<3000
 np.savetxt(root/f'{name}.csv',np.c_[x-x[0],y-y[0]],delimiter=',',header='x,y',comments='')
 tests.append(dict(csv=f'{name}.csv',role=role,channels={'y':1.}))
 provenance.append(dict(id=name,role=role,points=len(x),strain_offset=float(x[0]),stress_offset_kPa=float(y[0]),input_sha256=hashlib.sha256(xp.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(yp.read_bytes()).hexdigest()))
write_json(root/'starting_config.json',c)
write_json(root/'targets.json',dict(case=c['name'],units=c['units'],synthetic=False,conditions_confirmed=True,source='NHERI SimCenter quoFEM qfem-0028; UCSD Grade60 A6 group; Carreño et al 2020 DOI 10.1061/(ASCE)ST.1943-541X.0002505. Stress MPa converted to kPa, strain fraction as supplied.',tests=tests))
write_json(root/'provenance.json',dict(source='https://github.com/NHERI-SimCenter/quoFEM/tree/c30196dda348d0a33b2aeb508825c5f340a7521b/Examples/qfem-0028',commit='c30196dda348d0a33b2aeb508825c5f340a7521b',specimens=provenance,processing='Source every-20th-point files; no additional downsampling; initial small offsets removed; full chronological histories retained; 4 train / 2 holdout selected by ID before optimization.'))
write_json(root/'acceptance_plan.json',dict(criterion_defined_before_fitting=True,nrmse_max=.05,peak_error_max=.1,increment_error_max=.01,require_optimizer_success=True,scope=c['calibration_scope'],not_claimed=['grade population certification','buckling','fracture','fatigue'],notes='Project screening criteria. No closed-loop energy acceptance claim for these irregular open-ended paths.'))
print(provenance)
