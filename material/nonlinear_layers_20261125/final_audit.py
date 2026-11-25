"""Read-only scientific provenance and deliverable checks; save the audit result."""
import ast
import csv
import hashlib
import json
import math
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote,urlparse

import numpy as np
import openseespy.opensees as ops

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RUNS=ROOT/'ground/results/nonlinear_layers_20261125'
sys.path.insert(0,str(ROOT/'ground/ground_modeling'))
from ground_hybrid.config import default_config,validate_config
from ground_hybrid.materials import catalog
from sync_ground_catalog import expected


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


class Links(HTMLParser):
    def __init__(self):super().__init__();self.urls=[]
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ('href','src') and v:self.urls.append(v)


def main():
    checks={};materials=json.loads((HERE/'materials.json').read_text())['materials'];byid={m['id']:m for m in materials}
    run=json.loads((HERE/'results/summary.json').read_text());assert run['complete'] and run['count']==len(materials)
    assert all(c['passed'] and c['mass']['passed'] and c.get('returncode')==0 for c in run['cases'])
    assert sum(len(c['tests']) for c in run['cases'])==167
    for c in run['cases']:
        actual=json.loads((HERE/'results'/c['id']/'material.json').read_text())
        assert actual==byid[c['id']],('Specimen data drift',c['id'])
    checks['specimens']=dict(count=len(materials),loading_paths=167,mass_checks=len(materials),all_passed=True,exact_material_snapshots=True)
    source=expected();active=json.loads((ROOT/'ground/ground_modeling/ground_hybrid/configs/materials_nonlinear_20261125.json').read_text());assert source==active
    for body in ['Mars','Earth','Moon']:validate_config(default_config(body))
    config=default_config('Mars');cat=catalog()['materials'];assert all(cat[l['material']]['law']=='DruckerPrager' for l in config['layers'])
    assert sum(cat[l['material']]['law']=='ElasticIsotropic' for l in default_config('Mars',mars_profile='reference')['layers'])==3
    checks['active_catalog']=dict(specimen_parameters_identical=True,Mars_layers_all_DruckerPrager=True,selected=config['scenario'],reference_config_preserved=True)
    syntax=[]
    for p in [*HERE.glob('*.py'),ROOT/'ground/ground_modeling/ground_dome_footing.py',ROOT/'ground/ground_modeling/ground_20261107.py',ROOT/'ground/ground_modeling/ground_hybrid/materials.py',ROOT/'ground/ground_modeling/ground_hybrid/config.py']:
        ast.parse(p.read_text());syntax.append(str(p.relative_to(ROOT)))
    checks['syntax_files']=syntax
    cleanup=json.loads((HERE/'cleanup/plan.json').read_text());assert all(not (ROOT/x['path']).exists() for x in cleanup['deletions'])
    assert all((ROOT/x['path']).is_symlink() and sha(ROOT/x['path'])==x['sha256'] for x in cleanup['duplicate_aliases'])
    saved=ROOT/'ground/saved_candidates_20261118'
    assert all(sha(saved/p)==digest for p,digest in cleanup['saved_candidate_hashes'].items())
    original=json.loads((saved/'manifest.json').read_text())['copied_files'];assert all(sha(saved/x['path'])==x['sha256'] for x in original)
    assert all((ROOT/'ground'/p).is_file() for p in cleanup['preserved_active_inputs'])
    checks['cleanup']=dict(deleted_files=len(cleanup['deletions']),duplicate_paths_preserved=len(cleanup['duplicate_aliases']),saved_before_after_hash_checks=len(cleanup['saved_candidate_hashes']),original_saved_manifest_checks=len(original),active_inputs_retained=True)
    archive=json.loads((ROOT/'docs/references/_catalog/file_manifest.json').read_text());newids={'G17','G18','G19','G20','G21','R28','R29','R30','R31'};new=[r for r in archive if r.get('source_id') in newids]
    for r in new:
        p=ROOT/'docs/references'/r['file'];assert p.is_file() and sha(p)==r['sha256']
        if p.suffix=='.pdf':assert p.read_bytes().startswith(b'%PDF-')
    checks['literature']=dict(new_file_hash_checks=len(new),new_PDFs=sum(Path(r['file']).suffix=='.pdf' for r in new),unavailable_full_text='G20 Coleman2024 is a link only, not a downloaded PDF')
    stability=json.loads((HERE/'stability/stability_assessment.json').read_text())
    assert all(r['passed'] for r in stability['native_checks']) and all(r['passed'] for r in stability['determinant_formula_checks'])
    checks['additional_stability_diagnostic']=dict(native_tangent_checks=len(stability['native_checks']),direct_tensor_formula_checks=len(stability['determinant_formula_checks']),conditional_localization_concern=True,scope='Loading-tangent snapshot screen; not a global collapse or failure-cause proof')
    strain=json.loads((HERE/'strain_audit/strain_ranges.json').read_text())
    assert strain['mesh_coordinates_and_tags_match'] and strain['kinematic_verification']['native_patch_max_strain_error']<1e-12
    assert all(r['bbar_elastic_bulk_pressure_strain_error_relative']<1e-8 for r in strain['stage_material_strain_ranges'])
    checks['strain_postprocessing']=dict(native_warped_patch_verified=True,affine_completeness_verified=True,all_ground_GP_bulk_strain_pressure_consistent=True,stages=3,material_groups=len(strain['stage_material_strain_ranges']))
    pending=[];case_checks=[]
    for status in sorted(RUNS.rglob('execution_status.json')):
        execution=json.loads(status.read_text());case=str(status.parent.relative_to(RUNS))
        if execution.get('returncode') is None:pending.append(case);continue
        path=status.parent/'run/summary.json'
        if not path.exists():continue
        d=json.loads(path.read_text());assert d['ground_domain_m']==[90.,90.,30.]
        assert math.isclose(d['footing']['modeled_footprint_m2'],1739.26398308894,rel_tol=1e-12)
        assert d['loads']['downforce_N']==150300000. and d['loads']['pressurized_net_compression_N']==29580000.
        for name,record in d['ground_materials'].items():
            assert record['active_law']=='DruckerPrager'
            reference=byid[record['specimen_id']]['parameters']
            assert all(record[k]==reference[k] for k in ['E_Pa','G_Pa','K_Pa','nu','density_kg_m3','c_Pa','phi_deg','psi_deg'])
        snapshot=status.parent/'ground_dome_footing_snapshot.py'
        if not snapshot.exists():snapshot=RUNS/'source_snapshots/ground_dome_footing_initial_screen.py'
        assert sha(snapshot)==execution['entry_sha256'],('Entry snapshot mismatch',case)
        if execution['complete']:
            assert d['complete'] and execution['returncode']==0
            for stage in ['gravity','depressurized','pressurized']:
                assert math.isclose(d[stage]['final_factor'],1.,abs_tol=1e-8) and d[stage]['force_balance']['relative_error']<1e-5
            a=np.load(path.parent/'stage_displacements.npz');xyz=a['coordinates'];radius=np.linalg.norm(xyz[:,:2],axis=1)
            mask=np.isclose(xyz[:,2],-3.,atol=1e-6)&(radius>=13.-1e-6)&(radius<=27.+1e-6)
            for stage in ['depressurized','pressurized']:
                disp=a[stage+'_m'];assert np.isfinite(disp).all()
                maximum=float(np.max(-(disp-a['gravity_m'])[mask,2])*1000)
                assert math.isclose(maximum,d[stage]['maximum_mm'],rel_tol=1e-8),('Settlement recomputation differs',case,stage,maximum)
        case_checks.append(dict(case=case,complete=execution['complete'],parameters_and_entry_hash_match=True))
    checks['ground_cases']=case_checks;checks['pending_ground_cases']=pending
    checks['entrypoints']={}
    for name in ['column_static30','ground_static60','dome_default_geometry']:
        d=json.loads((RUNS/'entrypoint_checks'/name/'summary.json').read_text());assert d['complete']
        assert all(m['active_law']=='DruckerPrager' for m in d['ground_materials'].values())
        checks['entrypoints'][name]=dict(complete=True,elements=d['ground_elements'],analysis_performed=d.get('analysis_performed',name!='dome_default_geometry'))
    broken=[];linkcount=0
    for p in [HERE/'evaluation/report.html',RUNS/'comparison/report.html']:
        parser=Links();parser.feed(p.read_text())
        for url in parser.urls:
            parsed=urlparse(url)
            if parsed.scheme or not parsed.path:continue
            target=(p.parent/unquote(parsed.path)).resolve();linkcount+=1
            if not target.exists():broken.append(dict(report=str(p),url=url,target=str(target)))
    assert not broken,broken
    checks['report_local_links']=dict(checked=linkcount,broken=broken)
    checks['runtime_observed_at_audit']=dict(python=sys.version,python_executable=sys.executable,OpenSees=ops.version(),numpy=np.__version__)
    checks['complete']=not pending
    checks['scope']='File integrity, parameter transfer and implemented-formula checks; not field validation or global-stability certification'
    checks['independent_Mars_validation']=False
    (HERE/'final_audit.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(complete=not pending,specimens=len(materials),ground_runs_checked=len(case_checks),pending=pending,local_links=linkcount),ensure_ascii=False))


if __name__=='__main__':main()
