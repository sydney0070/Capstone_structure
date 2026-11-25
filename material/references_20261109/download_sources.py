"""Archive explicitly listed public references; retain failures, never save error pages as PDFs."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import urllib.request

BASE = Path(__file__).resolve().parent
JOBS = [
    ('musiol2016', '01_ground/musiol_2016_basalt.pdf', 'https://gfzpublic.gfz.de/rest/items/item_1802896_4/component/file_1802898/content', 'pdf'),
    ('spohn2022', '01_ground/spohn_2022_insight_hp3.pdf', 'https://link.springer.com/content/pdf/10.1007/s11214-022-00941-z.pdf', 'pdf'),
    ('lunar_profile', '01_ground/connolly_carrier_2023_lunar_profile.pdf', 'https://ntrs.nasa.gov/api/citations/20220018853/downloads/LPSC%202023%20Abstract_Connolly_Carrier_v2.pdf', 'pdf'),
    ('hydro6061', '03_candidates/hydro_6061_2019.pdf', 'https://www.hydro.com/globalassets/01-products--services/extruded-profiles/americas/ena-resources/alloy-data-sheets/hydro_2019_data_sheet_6061.pdf', 'pdf'),
    ('niar_as4', '03_candidates/niar_as4_8552_2011.pdf', 'https://www.wichita.edu/industry_and_defense/NIAR/Research/hexcel-8552/AS4-Unitape-2.pdf', 'pdf'),
    ('nist6061', '03_candidates/nist_6061_t6.html', 'https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061', 'html'),
    ('dp', '04_model_docs/DruckerPrager.html', 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html', 'html'),
    ('pdmy02', '04_model_docs/PressureDependentMultiYield02.html', 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html', 'html'),
    ('steel02', '04_model_docs/Steel02.html', 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Steel02.html', 'html'),
    ('concrete02', '04_model_docs/Concrete02.html', 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Concrete02.html', 'html'),
    ('flat_slider', '04_model_docs/FlatSliderBearing.html', 'https://opensees.berkeley.edu/wiki/index.php/Flat_Slider_Bearing_Element', 'html'),
    ('dp_engine_340', '04_model_docs/DruckerPrager_v3.4.0.cpp', 'https://raw.githubusercontent.com/OpenSees/OpenSees/v3.4.0/SRC/material/nD/UWmaterials/DruckerPrager.cpp', 'text'),
]


def fetch(job):
    key, rel, url, kind = job
    result = dict(id=key, url=url, destination=rel, attempted_on='2026-11-09 Asia/Seoul')
    dest = BASE / rel
    try:
        if dest.exists():
            data = dest.read_bytes()
            result['status'] = 'already_saved'
        else:
            request = urllib.request.Request(url, headers={'User-Agent': 'CapstoneReferenceArchive/1.0'})
            with urllib.request.urlopen(request, timeout=35) as response:
                data = response.read(40_000_001)
                result.update(content_type=response.headers.get('Content-Type'), final_url=response.geturl())
            if len(data) > 40_000_000:
                raise ValueError('response exceeds 40MB')
            if kind == 'pdf' and not data.startswith(b'%PDF-'):
                raise ValueError('response is not a PDF')
            if kind == 'html' and (len(data) < 1000 or any(s in data[:10000].lower() for s in
                    [b'just a moment', b'checking your browser', b'access denied', b'captcha'])):
                raise ValueError('response is not a usable HTML document')
            if kind == 'text' and b'DruckerPrager::' not in data:
                raise ValueError('unexpected implementation content')
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            result['status'] = 'saved'
        result.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    except Exception as error:
        result.update(status='unavailable', error=str(error), destination=None)
    return result


if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(fetch, JOBS))
    path = BASE / 'download_log.json'
    if path.exists():
        attempt = 1
        while path.with_name(f'download_log_attempt{attempt}.json').exists():
            attempt += 1
        path.rename(path.with_name(f'download_log_attempt{attempt}.json'))
    path.write_text(json.dumps(results, ensure_ascii=False, indent=2) + '\n')
    for result in results:
        print(result['id'], result['status'], result.get('bytes', result.get('error')))
