"""Keep the tested specimen parameters and active ground catalog identical.

Default is read-only checking. After intentionally changing materials.json,
--write transfers the exact values and invalidates the previous selection status.
It never silently recalibrates parameters or reuses old results as new validation.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONFIG=ROOT/'ground/ground_modeling/ground_hybrid/configs'


def expected():
    source=HERE/'materials.json';data=json.loads(source.read_text());out={}
    for r in data['materials']:
        if r['family']!='DruckerPrager':raise ValueError('This catalog is all-DP; extend the registry explicitly for another law')
        p=copy.deepcopy(r['parameters'])
        if p['hardening']!=0:raise ValueError('Current ground constructor is perfectly plastic; nonzero hardening must not be ignored')
        p.update(law=r['family'],specimen_id=r['id'],specimen_evidence=r['evidence'],calibration=r['calibration'],status=r['status'],role=r['layer'],assumptions=r['notes'],mars_calibrated=False)
        out['mars_nl25_'+r['id']]=p
    return dict(date=data['date'],specimen_source=str(source),specimen_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),materials=out)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true');args=parser.parse_args()
    payload=expected();catalog=CONFIG/'materials_nonlinear_20261125.json';old=json.loads(catalog.read_text())
    profiles_path=CONFIG/'Mars_nonlinear_20261125_profiles.json';profiles=json.loads(profiles_path.read_text())
    differences=[name for name,record in payload['materials'].items() if old['materials'].get(name)!=record]
    extra=sorted(set(old['materials'])-set(payload['materials']))
    source_changed=old.get('specimen_sha256')!=payload['specimen_sha256']
    for profile in profiles['profiles'].values():
        if any(layer['material'] not in payload['materials'] for layer in profile['layers']):raise ValueError('Profile references unregistered specimen')
        if not args.write and profile['specimen_source_sha256']!=payload['specimen_sha256']:raise ValueError('Profile/specimen provenance mismatch')
    if args.write:
        if source_changed or differences or extra:
            profiles['selection_status']='pending_specimen_and_ground_verification_after_catalog_change'
        for profile in profiles['profiles'].values():profile.update(specimen_source=payload['specimen_source'],specimen_source_sha256=payload['specimen_sha256'])
        mars_path=CONFIG/'Mars.json';mars=json.loads(mars_path.read_text());selected=profiles['selected_profile']
        mars.update(layers=profiles['profiles'][selected]['layers'],scenario=selected,selection_status=profiles['selection_status'])
        for path,data in [(catalog,payload),(profiles_path,profiles),(mars_path,mars)]:path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    elif differences or extra or source_changed:
        raise RuntimeError('Specimen/ground catalog mismatch; do not reuse old validation: '+str(dict(changed=differences,extra=extra,source_changed=source_changed)))
    print(json.dumps(dict(materials=len(payload['materials']),all_DruckerPrager=True,parameters_identical=not differences and not extra,source_hash_matches=not source_changed,wrote=args.write,selection_status=profiles['selection_status']),ensure_ascii=False))


if __name__=='__main__':main()
