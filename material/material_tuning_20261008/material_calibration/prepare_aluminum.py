"""Read-only XLSX extraction with bundled Python; writes separate analysis inputs."""
from pathlib import Path
import sys,io,json,hashlib,zipfile
import numpy as np
import openpyxl
B=Path(__file__).resolve().parent
sys.path.insert(0,str(B.parent))
from presets import get_case

def write(p,o):p.write_text(json.dumps(o,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
p=B/'raw/aluminum_TRR285-processed_data.zip'
assert hashlib.md5(p.read_bytes()).hexdigest()=='bf99e30aef562ef39ab807423e1f3d32'
z=zipfile.ZipFile(p)
for state in ['T4','T6']:
 root=B/f'aluminum_6014_{state.lower()}';root.mkdir(exist_ok=True)
 member=next(n for n in z.namelist() if f'EN-AW6014_{state}/base/' in n and n.endswith('.xlsx'))
 content=z.read(member);(root/'source.xlsx').write_bytes(content)
 w=openpyxl.load_workbook(io.BytesIO(content),read_only=True,data_only=True);s=w['Tabelle1'];rows=list(s.values)
 count=5 if state=='T4' else 3;train=3 if state=='T4' else 2
 tests=[];records=[]
 for j in range(count):
  k=j*6;assert rows[2][k+3].lower()=='mpa' and rows[2][k+4]=='%'
  a=np.array([[r[k+4]/100,r[k+3]*1000] for r in rows[3:] if isinstance(r[k+4],(int,float)) and isinstance(r[k+3],(int,float))])
  end=np.flatnonzero(a[:,0]>.02);end=int(end[0]) if len(end) else len(a)
  offset=a[0].copy();a=a[:end]-offset
  assert abs(offset[0])<1e-5 and len(a)>10
  name=f'specimen_{j+1}';role='train' if j<train else 'validation'
  np.savetxt(root/f'{name}.csv',a,delimiter=',',header='x,y',comments='')
  tests.append(dict(csv=f'{name}.csv',role=role,channels={'y':1.}))
  records.append(dict(id=rows[0][k+1],role=role,used_points=len(a),removed_initial_stress_kPa=float(offset[1]),max_strain=float(max(a[:,0])),source_columns=f'{openpyxl.utils.get_column_letter(k+4)}:{openpyxl.utils.get_column_letter(k+5)}'))
 c=get_case('shell_aluminum');c['name']=f'aluminum_aw6014_{state.lower()}_axial';c['params'].update(E=70e6,nu=.33,sig0=120e3 if state=='T4' else 200e3,sigInf=300e3,delta=80.,H=0.,thickness=.002)
 c['bounds']=dict(E=[60e6,80e6],sig0=[50e3,300e3],sigInf=[300e3,500e3],delta=[1.,400.])
 c['test']['substeps']=2;c['status']='EXPERIMENTAL_FIT_START_REQUIRES_ELASTIC_STRAIN_REVIEW'
 c['calibration_scope']='AW6014 '+state+', 23 C, rolling 0 degrees, 15 mm/min; engineering strain <=2%; axial monotonic response only. Elastic strain/compliance review required. nu=.33 assumed, H=0 fixed; no multiaxial/isotropy/cyclic validation.'
 write(root/'starting_config.json',c)
 write(root/'targets.json',dict(case=c['name'],units=c['units'],synthetic=False,conditions_confirmed=True,source='Boehnke et al DOI 10.5281/zenodo.21993433 v2; base series; supplied engineering stress MPa and strain percent, see source.xlsx and provenance.json.',tests=tests))
 write(root/'provenance.json',dict(source_doi='10.5281/zenodo.21993433',license='CC-BY-SA-4.0',source_sha256=hashlib.sha256(content).hexdigest(),zip_member=member,sheet='Tabelle1',specimens=records,processing='Stress MPa to kPa; strain percent to fraction; initial stress removed; stop before FIRST strain>0.02; preserve all earlier points including reversals/noise; no true-stress conversion for this small-strain phenomenological fit.'))
 write(root/'acceptance_plan.json',dict(criterion_defined_before_fitting=True,nrmse_max=.05,peak_error_max=.1,increment_error_max=.01,require_optimizer_success=True,scope=c['calibration_scope'],mandatory_review='ELASTIC_STRAIN_AND_MACHINE_COMPLIANCE_NOT_CONFIRMED',notes='Project screening; even numerical pass is not final material acceptance until strain provenance and elastic stiffness are reviewed.'))
 print(state,records)
