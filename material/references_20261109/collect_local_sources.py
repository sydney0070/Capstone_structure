"""Collect a dated evidence snapshot without changing model inputs or original evidence."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
import json
import shutil

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
CAL = ROOT / 'material/material_tuning_20261008/material_calibration'
RELEASE = ROOT / 'material/material_release_20261101'
copied = []


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save_json(relative, data):
    path = BASE / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')


def copy(source, relative, role):
    source = Path(source)
    if not source.is_file():
        raise FileNotFoundError(source)
    destination = BASE / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    original_hash = sha(source)
    if destination.exists() and sha(destination) != original_hash:
        raise ValueError(f'Existing snapshot differs: {destination}')
    if not destination.exists():
        shutil.copy2(source, destination)
    assert sha(destination) == original_hash
    copied.append(dict(source=str(source), destination=relative, bytes=source.stat().st_size,
                       sha256=original_hash, role=role))


def copy_top(source, destination, role):
    for path in sorted(source.iterdir()):
        if path.is_file() and path.suffix in {'.json', '.md', '.csv', '.xlsx', '.bib', '.txt'}:
            copy(path, str(Path(destination) / path.name), role)


catalog_files = {
    'release': ROOT / 'material/material_release_20261101/configs/catalog.json',
    'ground': ROOT / 'ground/ground_modeling/ground_hybrid/configs/materials.json',
}
catalogs = {key: json.loads(path.read_text()) for key, path in catalog_files.items()}
for key, path in catalog_files.items():
    copy(path, f'00_current_values/{key}_catalog.json', 'current_catalog_snapshot')
for body in ['Mars', 'Moon', 'Earth']:
    copy(ROOT / f'ground/ground_modeling/ground_hybrid/configs/{body}.json',
         f'00_current_values/{body}.json', 'base_scenario_not_runtime_override')

for relative in [
    'ground/ground_modeling/ground_20261107.py',
    'ground/ground_modeling/ground_hybrid/materials.py',
    'material/material_release_20261101/build_catalog.py',
    'material/material_release_20261101/capmat/materials.py',
    'material/materials_20261107.py',
    'structure/dome_in_mars_tuned.py',
    'structure/ring_foundation_rc.py',
    'structure/planet_foundations_20261109.py',
    'SSI/parameters.py', 'SSI/contact.py', 'SSI/connectors_20261109.py',
    'integrated/mars_ssi_integrated_20261109.py',
    'integrated/planet_foundation_demo_20261109.py',
]:
    copy(ROOT / relative, f'00_current_values/code_snapshot/{relative}', 'parameter_trace_snapshot_not_entrypoint')

for folder, destination in [
    ('jsc1a_lasp_review_20261019', '01_ground/jsc1a_lasp'),
    ('mars_jez1_strength', '01_ground/jez1_inputs'),
    ('regolith_sources_20261019', '01_ground/regolith_review'),
]:
    copy_top(CAL / folder, destination, 'ground_extraction_or_review')
for source, destination in [
    (CAL / 'regolith_sources_20261019/jez1_2021.pdf', '01_ground/jez1_2021.pdf'),
    (CAL / 'regolith_sources_20261019/nasa_guide_2024.pdf', '01_ground/nasa_lunar_simulant_guide_2024.pdf'),
    (CAL / 'jsc1a_lasp_review_20261019/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf', '01_ground/jsc1a_lasp/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf'),
    (RELEASE / 'research/sintered_jsc1a_leag2016.pdf', '03_candidates/sintered_jsc1a_leag2016.pdf'),
    (RELEASE / 'research/moon_strength_fits.json', '01_ground/moon_strength_fits.json'),
    (RELEASE / 'research/dp_mapping_note.md', '04_model_docs/dp_mapping_note.md'),
    (RELEASE / 'research/J2_SI_note.md', '04_model_docs/J2_SI_note.md'),
    (CAL / 'raw/Database_References.bib', '02_structure/s355/Database_References.bib'),
    (CAL / 'raw/Overall_Summary.csv', '02_structure/s355/Overall_Summary.csv'),
]:
    copy(source, destination, 'existing_reference_copy')
copy_top(CAL / 's355_hem320_web', '02_structure/s355', 'selected_s355_coupon_and_calibration_evidence')
for name in ['best_config.json', 'fit_report.json', 'manifest_used.json']:
    copy(CAL / f's355_hem320_web/fitted/{name}', f'02_structure/s355/fitted/{name}', 's355_fit_evidence')

for folder in ['rebar_ucsd_a6', 'aluminum_6014_t4', 'aluminum_6014_t6']:
    copy_top(CAL / folder, f'03_candidates/{folder}', 'historical_calibration_not_default_material')
copy_top(CAL / 'rebar_ucsd_a6/refinement_20261012', '03_candidates/rebar_ucsd_a6/refinement_20261012',
         'historical_rebar_refinement_not_sd400_default')
copy_top(CAL / 'raw/rebar_ucsd', '03_candidates/rebar_ucsd_raw', 'historical_rebar_source_data')
for name in ['stressStrain.zip', 'lunar_database_metadata.json', 'lunar_simulants.csv',
             'cfrp_smartfan_transversal.xlsx', 'aluminum_README.md',
             'aluminum_mapping-EN_AW6014_T4.yml', 'aluminum_mapping-EN_AW6014_T6.yml',
             'aluminum_TRR285-processed_data.zip']:
    copy(CAL / 'raw' / name, f'03_candidates/earlier_raw/{name}', 'candidate_data_not_default_calibration')

notes = [
    'docs/MODEL_OVERVIEW_20261109_KO.md',
    'docs/reference/GROUND_MATERIAL_REVIEW_20261107_KO.md',
    'material/docs/MATERIAL_EVALUATION.md',
    'material/material_tuning_20261008/SOURCES.md',
    'material/material_tuning_20261008/RESEARCH_REPORT_KO.md',
    'material/material_tuning_20261008/material_calibration/README.md',
    'material/material_tuning_20261008/material_calibration/MATERIAL_MAP_KO.md',
    'material/material_tuning_20261008/material_calibration/MATERIAL_STATUS_KO.md',
    'material/material_release_20261101/README.md',
    'research_20261109/artifacts/material_evidence_inventory.json',
]
for relative in notes:
    copy(ROOT / relative, f'05_project_notes/{relative}', 'internal_report_not_primary_experiment')
attachment = Path('/Users/sydney/Downloads/mars_moon_dome_foundation_research_summary.md')
if attachment.is_file():
    copy(attachment, '05_project_notes/mars_moon_dome_foundation_research_summary.md',
         'user_supplied_background_not_verified_material_coefficients')

copy(ROOT / 'ground/ground_modeling/stationxml.xb.elyse.2019.151.xml',
     '06_site_and_input/stationxml.xb.elyse.2019.151.xml', 'actual_station_metadata')
run = ROOT / 'ground/results/ground_20261109/run_20261109_174344_f30e3d'
if not run.is_dir():
    run = ROOT / 'ground/results/ground_20261109/layout_verified'
for name in ['motion_full_enu.json', 'motion_window_enu.json']:
    copy(run / name, f'06_site_and_input/{name}', 'processed_motion_metadata_from_named_completed_run')
summary = json.loads((run / 'summary.json').read_text())
save_json('06_site_and_input/ground_run_context.json',
          dict(source_run=str(run), site=summary['site'], parameters=summary['parameters'],
               resolution=summary['resolution'], note='A saved run context; not the universal default configuration.'))

all_materials = {}
for key, catalog in catalogs.items():
    for name, record in catalog['materials'].items():
        all_materials.setdefault(name, []).append(dict(catalog=key, record=record))
save_json('00_current_values/all_material_records.json', all_materials)

sources = {}
for key, catalog in catalogs.items():
    for name, record in catalog.get('sources', {}).items():
        existing = sources.setdefault(name, dict(id=name, catalog_records=[], local_files=[]))
        existing['catalog_records'].append(dict(catalog=key, **record))

sources['lasp_jsc1a']['url'] = 'https://ares.jsc.nasa.gov/projects/simulants/attachments/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf'
for key, local in {
    'jez1': ['01_ground/jez1_2021.pdf', '01_ground/jez1_inputs/published_inputs.json'],
    'lasp_jsc1a': ['01_ground/jsc1a_lasp/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf',
                  '01_ground/jsc1a_lasp/triaxial_observations.json', '01_ground/moon_strength_fits.json'],
    's355_coupon': ['02_structure/s355/provenance.json', '02_structure/s355/reviewed_config.json'],
    'sintered2016': ['03_candidates/sintered_jsc1a_leag2016.pdf'],
    'nasa_lunar_guide2024': ['01_ground/nasa_lunar_simulant_guide_2024.pdf'],
}.items():
    sources[key]['local_files'].extend(local)

downloads = json.loads((BASE / 'download_log.json').read_text())
for download in downloads:
    key = download['id']
    source = sources.setdefault(key, dict(id=key, catalog_records=[], local_files=[]))
    source['download'] = download
    if download['destination']:
        source['local_files'].append(download['destination'])
    source.setdefault('url', download['url'])

sources['spohn2022']['publisher_url'] = 'https://link.springer.com/article/10.1007/s11214-022-00941-z'
sources['s355_coupon']['title'] = 'Database of Uniaxial Cyclic and Tensile Coupon Tests for Structural Metallic Materials, v1.0.0'
sources['spohn2022']['title'] = 'The InSight HP3 Penetrator (Mole) on Mars: Soil Properties Derived from the Penetration Attempts and Related Activities'
sources['musiol2016']['title'] = 'Lithospheric flexure and gravity spreading of Olympus Mons volcano, Mars'
for key, source in sources.items():
    source['used_by_materials'] = sorted(name for name, variants in all_materials.items()
        if any(key in item['record'].get('sources', []) for item in variants))
save_json('source_registry.json', sources)
save_json('copy_manifest.json', copied)

lines = ['# 참고자료 목록', '', '정리일: 2026-11-09 (Asia/Seoul). 원 카탈로그의 출처 연결을 보존했다. '
         '이 목록은 모든 수치를 새로 검증했다는 뜻이 아니다.', '',
         '| 출처 ID | 사용 재료 / 목적 | 원문 | 모은 파일 |', '|---|---|---|---|']
for key, source in sources.items():
    original = source.get('publisher_url') or source.get('url') or next(
        (item.get('url') for item in source['catalog_records'] if item.get('url')), None)
    materials = ', '.join(source['used_by_materials']) or '구성 법칙·접촉 구현 참고'
    local = '<br>'.join(f'[{Path(file).name}]({file})' for file in source['local_files']) or '원문 링크만 보관; 로컬 원문 미확보'
    lines.append(f'| `{key}` | {materials} | ' + (f'[공식 출처 / DOI]({original})' if original else '아래 원 카탈로그 참조') + f' | {local} |')
lines += ['', '## 함께 보관한 이전 조사 자료', '',
          '- [이전 전체 참고문헌](05_project_notes/material/material_tuning_20261008/SOURCES.md): 재료·기초 형식의 배경 조사. 목록의 모든 논문에서 계수를 가져온 것은 아니다.',
          '- [사용자 제공 기초 선행연구 요약](05_project_notes/mars_moon_dome_foundation_research_summary.md): 화성 전면기초·달 링기초의 개념 참고.',
          '- UCSD 철근·6014 알루미늄·Smartfan CFRP·micro-concrete·달 시뮬런트 데이터베이스는 `03_candidates/`의 이전 조사·보정 자료다. 현재 SD400·6061·AS4 카탈로그와 혼용하지 않는다.',
          '- 기존 보고서 사본은 당시 문맥을 보존했다. 그 내부의 옛 상대 경로는 이동으로 끊길 수 있으므로 이 목록과 `copy_manifest.json`으로 실제 사본을 찾는다.',
          '- HTML 사본은 원문 응답 보관용이며 CSS·그림 등 외부 자산까지 저장한 완전한 오프라인 패키지는 아니다.',
          '- NIAR AS4 원문은 웹 도구에서 조회했으나 파일 다운로드 응답이 PDF가 아니어서 저장하지 않았다. 원문을 확보한 것으로 표시하지 않는다.',
          '- 최초 일반 실행의 DNS 제한 기록은 `download_log_attempt1.json`, 네트워크 다운로드 결과는 `download_log.json`에 남겼다.', '']
(BASE / 'SOURCES_KO.md').write_text('\n'.join(lines))
print(json.dumps(dict(copied_files=len(copied), unique_material_ids=len(all_materials), source_ids=len(sources),
                      download_saved=sum(bool(item['destination']) for item in downloads)), ensure_ascii=False))
