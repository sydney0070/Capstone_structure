"""Extract published summary observations; never treat them as fitted parameters.

Requires pypdf. Preserves source values, missing entries, and ambiguous headings.
Output directories must be new. No model configuration is changed.
"""
import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

from pypdf import PdfReader


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    reader = PdfReader(args.source)
    if len(reader.pages) != 74:
        raise ValueError('Expected 74-page LASP report; review source before extraction')
    observations = []
    for page_number, table, expected in [(21, 5, 18), (22, 6, 13), (23, 7, 15)]:
        rows = []
        for line in reader.pages[page_number - 1].extract_text().splitlines():
            if not re.match(r'^J1A_\d+_\d+[A-Z]\s', line):
                continue
            cells = line.split()
            if len(cells) != 12:
                raise ValueError(f'Unexpected row on page {page_number}: {line}')
            name, date, ton, *values = cells
            v = [None if x in ('n/a', '--') else float(x) for x in values]
            fields = ['density_g_cm3', 'phi_peak_c0_deg',
                      'column_after_phi_peak_as_printed', 'sigma1_peak_kPa',
                      'sigma3_peak_kPa', 'Etan_MPa', 'psi_deg',
                      'phi_residual_c0_deg', 'axial_strain_at_residual']
            row = dict(zip(fields, v))
            row.update(test_id=name, date_as_printed=date, batch_ton=int(ton),
                       pdf_page=page_number, printed_page=page_number-2,
                       table=table, raw_extracted_line=line,
                       nominal_confinement_kPa=int(name.split('_')[1]),
                       nominal_density_g_cm3=int(name.split('_')[2][:2])/10,
                       status='PUBLISHED_SUMMARY_OBSERVATION_NOT_CALIBRATION')
            row['column_after_phi_peak_heading'] = (
                'Etan, MPa (source heading inconsistent with neighboring tables; do not silently repair)'
                if table == 6 else 'Axial Strain, epsilon_ax')
            row['Etan_kPa_unit_conversion'] = None if row['Etan_MPa'] is None else row['Etan_MPa']*1000
            a, b = row['sigma1_peak_kPa'], row['sigma3_peak_kPa']
            # Diagnostic only: friction angle implied by printed principal stresses at c=0.
            angle = math.degrees(math.asin((a-b)/(a+b)))
            row['phi_c0_from_printed_stresses_deg_diagnostic'] = angle
            row['phi_difference_deg_diagnostic'] = angle-row['phi_peak_c0_deg']
            row['stress_angle_review_flag'] = abs(row['phi_difference_deg_diagnostic']) > 0.3
            rows.append(row)
        if len(rows) != expected:
            raise ValueError(f'Table {table}: got {len(rows)}, expected {expected}')
        observations.extend(rows)
    assert len({r['test_id'] for r in observations}) == 46
    digest = hashlib.sha256(args.source.read_bytes()).hexdigest()
    metadata = {
        'status': 'SOURCE_REVIEW_ONLY_NOT_FITTED_OR_VALIDATED',
        'source_file': args.source.name, 'source_sha256': digest,
        'source_url': 'https://ares.jsc.nasa.gov/projects/simulants/attachments/JSC-1A_LASP_Geotech_Properties_Tons1-3.pdf',
        'report': 'LASP Document 105525, Rev A, 2007-10-04',
        'observations': len(observations),
        'conditions': {'material': 'JSC-1A, Earth laboratory simulant',
                       'nominal_density_g_cm3': [1.7, 1.8, 1.9],
                       'nominal_confinement_kPa': [15, 30, 60],
                       'loading': 'monotonic triaxial compression, dry tamped specimens, pore space open to atmosphere',
                       'specimen_diameter_mm': 74.6, 'specimen_height_mm': 150,
                       'displacement_rate_in_per_hour_as_reported': 0.0197,
                       'condition_source_pdf_pages': [7, 8, 12]},
        'limitations': [
            '46 tabulated test summaries, not 46 raw stress-strain histories.',
            'Italic formatting and source-selected averaging are not encoded; do not auto-select training rows.',
            'Table 6 has an inconsistent heading; original column retained without relabeling as confirmed strain.',
            'Text on PDF page 20 reports phi=42.9 deg at 1.7 g/cm3; summary on page 28 reports 43.9 deg. Neither selected automatically.',
            'Missing values remain null. No zero filling and no manufactured curves.',
            'Source reports specimen working/retesting may reduce strength; test IDs do not establish independent virgin specimens.',
            'Etan is not established as small-strain dynamic modulus Gmax or a universal pressure-independent elastic modulus.',
            'Cyclic response, damping, uplift contact and lunar in-situ conditions are not validated.',
        ],
        'stress_angle_flags': [r['test_id'] for r in observations if r['stress_angle_review_flag']],
        'diagnostic_formula': 'phi(c=0)=asin((sigma1-sigma3)/(sigma1+sigma3)); flag >0.3 deg; no source value overwritten',
    }
    args.out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(args.source, args.out/args.source.name)
    for name, data in [('triaxial_observations.json', observations), ('provenance.json', metadata)]:
        (args.out/name).write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n', encoding='utf8')
    assert hashlib.sha256((args.out/args.source.name).read_bytes()).hexdigest() == digest
    print(json.dumps({'out': str(args.out), 'rows': len(observations),
                      'stress_angle_flags': metadata['stress_angle_flags'],
                      'source_hash_verified': True}, indent=2))


if __name__ == '__main__':
    main()
