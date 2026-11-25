"""Verify source copies, material-source links and the new evidence guide."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
import json
import re
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


copies = json.loads((BASE / 'copy_manifest.json').read_text())
for row in copies:
    assert sha(Path(row['source'])) == row['sha256'] == sha(BASE / row['destination']), row['destination']
registry = json.loads((BASE / 'source_registry.json').read_text())
for record in registry.values():
    assert all((BASE / path).is_file() for path in record['local_files'])
records = json.loads((BASE / '00_current_values/all_material_records.json').read_text())
for variants in records.values():
    for variant in variants:
        assert set(variant['record'].get('sources', [])).issubset(registry)
pdfs = []
for path in sorted(BASE.rglob('*.pdf')):
    content = path.read_bytes()
    assert content.startswith(b'%PDF-') and b'%%EOF' in content[-4096:], path
    pdfs.append(dict(path=str(path.relative_to(BASE)), bytes=len(content), sha256=sha(path)))
downloads = json.loads((BASE / 'download_log.json').read_text())
for row in downloads:
    if row['destination']:
        assert sha(BASE / row['destination']) == row['sha256']
for name in ['README.md', 'ADOPTED_VALUES_KO.md', 'SOURCES_KO.md', 'GAPS_AND_ASSUMPTIONS_KO.md']:
    for target in re.findall(r'\]\(([^)]+)\)', (BASE / name).read_text()):
        if target.startswith(('https://', 'http://', '#')) or target == 'verification.json':
            continue
        assert (BASE / target.split('#')[0]).exists(), (name, target)
sys.path.insert(0, str(ROOT / 'material'))
from materials_20261107 import audit_catalogs
catalog_audit = audit_catalogs()
assert catalog_audit['checked_records'] == len(records) == 29
source_hashes = []
for name in ['release', 'ground']:
    catalog = json.loads((BASE / f'00_current_values/{name}_catalog.json').read_text())
    for key, source in catalog['sources'].items():
        if source.get('local') and source.get('sha256'):
            matches = sha(Path(source['local'])) == source['sha256']
            assert matches, (name, key)
            source_hashes.append(dict(catalog=name, id=key, matches_recorded_hash=matches))
for item in json.loads((BASE / '06_site_and_input/large_originals.json').read_text()):
    assert sha(Path(item['path'])) == item['sha256']
verification = dict(
    completed_on=datetime.now(ZoneInfo('Asia/Seoul')).isoformat(),
    copied_files=len(copies), all_copies_equal_original=True,
    downloaded_files=sum(bool(item['destination']) for item in downloads),
    unavailable_downloads=[item['id'] for item in downloads if not item['destination']],
    pdf_files=pdfs, source_ids=len(registry), unique_material_ids=len(records),
    all_material_source_ids_resolved=True, all_registry_local_files_exist=True,
    new_navigation_links_valid=True, large_original_hashes_verified=True,
    original_recorded_source_hashes=source_hashes, catalogue_audit=catalog_audit,
    scope='Reference organization and integrity only; no model parameters or analysis implementation changed.',
)
(BASE / 'verification.json').write_text(json.dumps(verification, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(dict(copied_files=len(copies), downloaded_files=verification['downloaded_files'],
                     pdf_files=len(pdfs), material_ids=len(records), checks='passed'), ensure_ascii=False))
