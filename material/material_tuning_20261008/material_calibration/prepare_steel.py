"""Reproducible experimental S355 coupon preparation; no structural model."""
from pathlib import Path
import sys,csv,zipfile,hashlib,io,json
import numpy as np
from contextlib import nullcontext
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE.parent))
from presets import get_case
from tune import write_json

if __name__=='__main__':
    archive=BASE/'raw/Clean_Data_v1-0-0.zip'
    if archive.exists():
        assert hashlib.md5(archive.read_bytes()).hexdigest()=='fece95a8a88c54f39ef9ffc518922b31'
    meta={r['hidden_index']:r for r in csv.DictReader((BASE/'raw/Overall_Summary.csv').open())}
    root=BASE/'s355_hem320_web';root.mkdir(exist_ok=True)
    records=[];tests=[]
    # Independent coupons from the same documented grade/source/composition.
    # Test-role selection occurs BEFORE optimization and is not tuned to its scores.
    prior=json.loads((root/'provenance.json').read_text()) if not archive.exists() else None
    with (zipfile.ZipFile(archive) if archive.exists() else nullcontext()) as z:
        for idx,role in [('0','train'),('8','train'),('3','validation'),('6','validation')]:
            m=meta[idx]
            name=f'{role}_{m["lp"]}_id{idx}'
            content=z.read(m['file']) if z else (root/f'{name}_source.csv').read_bytes()
            if not z:
                expected=next(r['source_sha256'] for r in prior['specimens'] if r['db_tag']==idx)
                assert hashlib.sha256(content).hexdigest()==expected
            assert m['grade']=='S355' and m['spec']=='J2+M' and m['source']=='HEM320 web'
            a=np.loadtxt(io.BytesIO(content),delimiter=',',skiprows=1)
            # Authors supply true stress (MPa), log strain. Invert uniform-area
            # conversion to nominal measures for the small-strain OpenSees coupon.
            strain=np.expm1(a[:,2]);stress=a[:,3]/np.exp(a[:,2])*1000.
            cutoff=np.flatnonzero(abs(strain)>.02)
            end=int(cutoff[0]) if len(cutoff) else len(a)
            x=strain[:end]-strain[0];y=stress[:end]-stress[0]
            assert abs(strain[0])<1e-4 and abs(stress[0])<1000.
            assert np.isfinite(x).all() and np.isfinite(y).all() and end>=10
            name=f'{role}_{m["lp"]}_id{idx}'
            (root/f'{name}_source.csv').write_bytes(content)
            np.savetxt(root/f'{name}.csv',np.c_[x,y],delimiter=',',header='x,y',comments='')
            record=dict(db_tag=idx,role=role,metadata=m,source_sha256=hashlib.sha256(content).hexdigest(),
                raw_points=len(a),used_points=end,engineering_strain_limit=.02,
                removed_initial_strain=float(strain[0]),removed_initial_stress_kPa=float(stress[0]),
                processing='Preserved acquisition order; uniform-area true-to-engineering conversion; initial offset removed; stop at FIRST abs(strain)>0.02. No sorting or later-cycle recovery.',
                max_abs_used_strain=float(max(abs(x))))
            records.append(record);tests.append(dict(csv=f'{name}.csv',role=role,channels={'y':1.}))
    c=get_case('rebar_steel02');c['name']='structural_steel_s355_hem320_web'
    c['status']='EXPERIMENTAL_CALIBRATION_START_NOT_REBAR_GRADE'
    c['params'].update(E=205e6,Fy=400000.,b=.01,R0=18.,a1=.02,a3=.02)
    c['bounds']=dict(E=[170e6,230e6],Fy=[330000.,460000.],b=[.001,.03],R0=[8.,30.],a1=[0.,.12],a3=[0.,.12])
    c['test']['substeps']=1
    write_json(root/'starting_config.json',c)
    write_json(root/'provenance.json',dict(dataset_doi='10.5281/zenodo.6965147',version='1.0.0',license='CC-BY-4.0',
        authors='Hartloper, Ozden, de Castro e Sousa, Lignos (2022)',
        experiment_reference_doi='10.5075/epfl-thesis-8700',specimens=records,
        applicability='Room-temperature S355 J2+M HEM320 web coupons; not rebar, aluminum, welded joints, buckling, fracture or extraterrestrial environment.',
        experimental_data_origin='Author-processed/downsampled experimental records; NOT synthetic, NOT digitized here.'))
    write_json(root/'targets.json',dict(case=c['name'],units=c['units'],synthetic=False,conditions_confirmed=True,
        source='Hartloper et al. 2022, DOI 10.5281/zenodo.6965147; El Jisr 2022 DOI 10.5075/epfl-thesis-8700; S355 J2+M HEM320 web; see provenance.json for conversion and trial IDs.',tests=tests))
    write_json(root/'acceptance_plan.json',dict(status='PROJECT_SCREENING_CRITERIA_NOT_DESIGN_STANDARD',
        validation_range_nrmse_max=.05,validation_peak_stress_relative_error_max=.10,
        require_optimizer_success=True,require_independent_specimen=True,
        criterion_defined_before_fitting=True,
        not_claimed=['material grade population certification','fracture','fatigue','temperature/vacuum','multiaxial J2 calibration']))
    print([(r['db_tag'],r['role'],r['used_points']) for r in records])
