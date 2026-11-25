"""Explicit, hash-audited removal of unused ground prototypes and duplicate data.

No recursive project clean: active DEMs, saved candidate files, analysis evidence
and imported legacy geometry helpers are deliberately outside the deletion list.
"""
import hashlib
import json
import os
import unicodedata
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
GROUND=ROOT/'ground'
AUDIT=HERE/'cleanup'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    AUDIT.mkdir(exist_ok=True)
    if (AUDIT/'completed.json').exists():raise RuntimeError('Cleanup already executed; read its audit instead of repeating')
    deletions=[];aliases=[]
    def remove(p,reason):
        if not p.is_file() or p.is_symlink():return
        assert p.is_relative_to(GROUND) and 'saved_candidates_20261118' not in p.parts
        deletions.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p),reason=reason))
    for name in ['earth','earth.py','soildata_fromnasa','pohang_nasadem.py']:
        remove(GROUND/name,'Standalone obsolete prototype; no active Python importer/caller; active multi-planet entrypoints and geometry helpers retained')
    for p in GROUND.iterdir():
        if p.is_file() and (p.name.startswith('pohang_') and p.suffix in ('.html','.csv')):
            remove(p,'Old Pohang prototype output, not an input to the current ground pipeline')
        if unicodedata.normalize('NFC',p.name) in ('NASADEM_HGT_n37e126','포항 지반 데이터') and p.is_dir():
            for f in p.rglob('*'):
                remove(f,'Unused n37e126 terrain / unread .num and .swb ancillary masks; active n36e129.hgt retained')
    remove(GROUND/'n37e126.hgt','Unused terrain tile; no active n37e126 reference; current Earth terrain uses n36e129.hgt')
    for p in GROUND.rglob('*'):
        if p.is_file() and (p.name=='.DS_Store' or p.suffix=='.pyc'):
            if not any(d['path']==str(p.relative_to(ROOT)) for d in deletions):
                remove(p,'Regenerable filesystem metadata or Python bytecode cache')
    def alias(duplicate,canonical):
        if not duplicate.is_file() or not canonical.is_file() or duplicate.is_symlink():return
        assert duplicate.is_relative_to(GROUND) and 'saved_candidates_20261118' not in duplicate.parts
        a=sha(duplicate);b=sha(canonical)
        if a!=b:raise RuntimeError('Duplicate content differs: '+str(duplicate))
        aliases.append(dict(path=str(duplicate.relative_to(ROOT)),canonical=str(canonical.relative_to(ROOT)),bytes=duplicate.stat().st_size,sha256=a,reason='Byte-identical duplicate; relative symlink retains old report/import paths'))
    archive=GROUND/'archive/ground_modeling_20261116'
    alias(archive/'stationxml.xb.elyse.2019.151.sxml',GROUND/'ground_modeling/stationxml.xb.elyse.2019.151.xml')
    for candidate in ['mgs1','mgs1s']:
        for ext in ['html','json']:
            alias(GROUND/f'results/ground_dome_footing/analysis_{candidate}_20261118/run/footing_settlement.{ext}',GROUND/f'saved_candidates_20261118/{candidate}/ground_analysis/footing_settlement.{ext}')
    for p in (archive/'insight_verification/material').glob('*.csv'):
        other=archive/'hybrid_verification/material'/p.name
        if other.is_file() and sha(p)==sha(other):alias(p,other)
    alias(archive/'insight_results_quick/insight_acceleration_enu.csv',archive/'insight_results/insight_acceleration_enu.csv')
    history=GROUND/'results/ground_20261109'
    canonical=history/'run_20261109_174344_f30e3d'
    for name in ['layout_verified','viewer_split_verified','viewer_split_verified_v2']:
        for filename in ['ground_response.html','ground_response.json','motion_full_enu.csv','site_dem_samples.csv']:
            p=history/name/filename;other=canonical/filename
            if p.is_file() and other.is_file() and sha(p)==sha(other):alias(p,other)
    protected=GROUND/'saved_candidates_20261118'
    saved_before={str(p.relative_to(protected)):sha(p) for p in protected.rglob('*') if p.is_file()}
    plan=dict(date='2026-11-25',authorization='User requested removal of unused ground code/data',deletions=deletions,duplicate_aliases=aliases,saved_candidate_hashes=saved_before,
              preserved_active_inputs=['ground_modeling/ldem_1024_90s_75s_000_030.img','ground_modeling/DTEEC_037262_1845_036761_1845_U01.IMG','ground_modeling/megt00n000hb.img','ground_modeling/n36e129.hgt','ground_modeling/stationxml.xb.elyse.2019.151.xml'],
              preserved_legacy_dependencies=['ground_modeling/ground_modeling.py','ground_modeling/ground_modeling_insight_updated.py'])
    (AUDIT/'plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')
    for item in deletions:
        p=ROOT/item['path']
        if sha(p)!=item['sha256']:raise RuntimeError('File changed after audit: '+str(p))
        p.unlink()
    for item in aliases:
        p=ROOT/item['path'];target=ROOT/item['canonical']
        if sha(p)!=item['sha256'] or sha(target)!=item['sha256']:raise RuntimeError('Duplicate changed')
        p.unlink();p.symlink_to(os.path.relpath(target,p.parent))
        if sha(p)!=item['sha256']:raise RuntimeError('Alias content mismatch')
    for p in sorted([p for p in GROUND.rglob('*') if p.is_dir()],key=lambda x:len(x.parts),reverse=True):
        if not any(p.iterdir()) and (p.name=='__pycache__' or p.name=='insight_opensees_integration' or unicodedata.normalize('NFC',p.name) in ('NASADEM_HGT_n37e126','포항 지반 데이터')):p.rmdir()
    saved_after={str(p.relative_to(protected)):sha(p) for p in protected.rglob('*') if p.is_file()}
    if saved_before!=saved_after:raise RuntimeError('Saved candidates changed')
    result=dict(deleted_files=len(deletions),deduplicated_files=len(aliases),removed_bytes=sum(x['bytes'] for x in deletions)+sum(x['bytes'] for x in aliases),saved_candidate_files_verified=len(saved_before),saved_candidates_unchanged=True,duplicate_aliases_verified=len(aliases),plan=str(AUDIT/'plan.json'))
    (AUDIT/'completed.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
