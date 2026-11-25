# Uniaxial tensile test data for EN-AW6014 (T4, T6) and HCT590X sheet metal

## Summary

This dataset contains uniaxial tensile test data for two sheet metal materials commonly used in the automotive sector: the aluminum alloy EN-AW6014 in T4 and T6 state, and the dual-phase steel HCT590X. The data was generated within the framework of the Collaborative Research Centre TRR 285.

It comprises **43 test series covering 189 individual specimens**. Engineering stress-strain curves are provided for every specimen, alongside the underlying force-displacement records. Depending on the test regime, two kinds of full-field/optical raw data accompany the mechanical curves:

- **Digital image correlation** (DIC, via GOM Aramis) for six quasi-static and elevated-temperature series, supplied as raw `.aramis` files.
- **High-speed camera image sequences** (Photron) for the 21 high-rate series, supplied as `.cihx` sequence headers with their `.tif` frames and a companion `.csv` per sequence.

Test conditions were systematically varied around a reference condition (SEP 1230 specimen geometry, 0° to rolling direction, 15 mm/min crosshead speed, 23 °C) to study the influence of:

- **Specimen orientation:** 0°, 45°, 90° relative to rolling direction
- **Crosshead speed:** 1.5 mm/min, 15 mm/min (reference), 0.25 m/s, 2.5 m/s, 6.0 m/s; corresponding to **nominal strain rates** of 0.001, 0.01, 10, 100 and 240 s⁻¹
- **Temperature:** 23 °C (reference), 100 °C, 150 °C, 200 °C; HCT590X and EN-AW6014 T4 only
- **Sheet thickness** (HCT590X only): 0.8 mm, 1.0 mm, 1.5 mm; EN-AW6014 was tested at 2.0 mm
- **Specimen geometry** (HCT590X only): SEP 1230, DIN 50125

Room-temperature tests follow ISO 6892-1, elevated-temperature tests follow ISO 6892-2.

> **Note on the DIN 50125 series:** these were run at the same crosshead speed of 15 mm/min as the SEP 1230 reference, but the differing parallel length yields a nominal strain rate of 0.003 s⁻¹ instead of 0.01 s⁻¹.

Each test series is identified by a unique ID (e.g. `TS05`) and documented in a per-material YAML mapping file. The mapping files link test conditions to ontology terms from the PMD Core Ontology (PMDco) and the Tensile Test Ontology (TTO) to support FAIR-compliant reuse and automated processing.

`TRR285_metadata_master_v04.xlsx` is the central description of the dataset. It has two sheets, joined by the series ID:

1. `test_series`, with one row per test series, and
2. `files` with one row per file, linked to its series via `series_id`

Each row of `test_series` carries the identifiers (`material`, `material_condition`, `specimen_geometry`, `test_standard`, `test_type`), the test conditions, the number of specimens (`n_specimens_per_orientation`, either 3 or 5). Every physical quantity is stored as a value/unit pair: `temperature_value` + `temperature_unit`, `test_speed_value` + `test_speed_unit`, and likewise for sheet thickness, and orientation. As geometry independent measure of deformation rate nominal strain rate is calculated from test speed $v_0$ and parallel length $l_p$ as $\dot{\varepsilon}  = v_0 / l_p$, and stored in `nominal_strain_rate_value`. `material_condition` is populated only for the aluminum (`EN-AW6014 T4` / `T6`).

Each row of `files` is one file, keyed by `file_id` (e.g. `TS05_AR01`), with `series_id` linking back to the conditions. Three columns classify the content — `file_type`, `data_content` and `processing_level` (`raw`, `processed`, or `raw and processed`) — and `software` records the application that produced or reads it. Exact location in the archive is given by further columns; see [Locating a file via the metadata](#locating-a-file-via-the-metadata) below.

## Contents

| File | Size | Contents |
| :--- | ---: | :--- |
| `TRR285_metadata_master_v04.xlsx` | 39 KB | Master metadata: one row per test series, one row per file |
| `mapping-HCT590X.yml` | 4 KB | Ontology mapping, HCT590X |
| `mapping-EN_AW6014_T4.yml` | 3 KB | Ontology mapping, EN-AW6014 T4 |
| `mapping-EN_AW6014_T6.yml` | 3 KB | Ontology mapping, EN-AW6014 T6 |
| `TRR285-processed_data.zip` | 7 MB | **All 43 stress-strain workbooks.** Start here |
| `HCT590X-base.zip` | 5.9 GB | HCT590X reference condition, incl. DIC |
| `HCT590X-ori.zip` | 11.1 GB | HCT590X orientation variation, incl. DIC |
| `HCT590X-temp.zip` | 9.7 GB | HCT590X elevated temperature, incl. DIC |
| `HCT590X-velo.zip` | 1.2 GB | HCT590X speed variation, incl. high-speed camera |
| `HCT590X-geom.zip` | < 1 MB | HCT590X DIN 50125 geometry |
| `HCT590X-thick.zip` | < 1 MB | HCT590X sheet thickness variation |
| `EN_AW6014_T4.zip` | 2.0 GB | Complete EN-AW6014 T4 data |
| `EN_AW6014_T6.zip` | 1.7 GB | Complete EN-AW6014 T6 data |

**If you only need the mechanical curves, download `TRR285-processed_data.zip` (7 MB).** It contains the stress-strain workbook for all 43 series and requires none of the large archives.

The HCT590X data is split by variation group rather than shipped as one archive, because a single whole-material file exceeded practical upload limits.

## Data organisation

Within each archive the layout is:

```
<material>/<variation group>/<level>/
    TSxx_EX01-<condition>.xlsx          # mechanical data, stored uncompressed
    TSxx-supplementary_data.zip         # optical raw data for that series (nested archive)
```

The variation groups are `base/` (reference condition), `ori/` (orientation), `velo/` (crosshead speed), `temp/` (temperature), `thick/` (sheet thickness) and `geom/` (specimen geometry). Note that `velo/` is keyed by **crosshead speed**, not by strain rate, because the speed is the quantity actually set on the testing machine: the same 15 mm/min yields a nominal strain rate of 0.01 s⁻¹ on SEP 1230 specimens but 0.003 s⁻¹ on DIN 50125 ones, so strain rate does not partition the series cleanly. The nominal strain rate of every series is given in `test_series.nominal_strain_rate_value`; where the manuscript refers to strain rate, the corresponding speeds are those listed in the summary above.

Supplementary (optical) data is kept in a **nested archive per test series**, so the outer archive unpacks quickly and only the series you actually need has to be expanded. Inside a supplementary archive, data is grouped per specimen:

```
TSxx-supplementary_data/
    TSxx_ORnn-<condition>/              # one folder per specimen
        TSxx_CI0n-<condition>.cihx      # image sequence header
        TSxx_CS0n-<condition>.csv       # companion data for the sequence
        *.tif                           # individual frames
```

## Data types and processing levels

| Type | Suffix | Count | Level | Description |
| :--- | :--- | ---: | :--- | :--- |
| Mechanical | `.xlsx` | 43 | raw + processed | Force-displacement (raw) and stress-strain (processed), one workbook per series covering all its specimens |
| DIC | `.aramis` | 26 | raw | GOM Aramis full-field measurement, one file per specimen |
| Image sequence | `.cihx` | 105 | raw | Photron high-speed camera sequence header |
| Sequence companion | `.csv` | 105 | processed | Accompanying data for the corresponding image sequence |

The 279 files above are the units declared in the metadata. The **64,753 individual `.tif` frames** that make up the image sequences are not listed individually; each sequence is indexed by its `.cihx` header. The frames account for the bulk of the dataset volume.

Units in the mechanical workbooks: force in kN, displacement in mm, stress in MPa, strain in %.

## Locating a file via the metadata

Because most files sit inside a nested archive, three columns of the `files` sheet together give the location:

| Column | Meaning |
| :--- | :--- |
| `file_path` | Folder holding the file, or holding the archive that contains it |
| `container_archive` | Supplementary archive containing the file — **blank if the file is loose** |
| `path_in_archive` | Path inside that archive — **blank if the file is loose** |

```
loose file :  <file_path>/<file_name><file_suffix>
in archive :  <file_path>/<container_archive>  ::  <path_in_archive>
```

Also see the examples below:

```
# loose - the 43 mechanical workbooks are the only loose files
file_path         = TRR285_Materials_Testing\HCT590X\thick
file_name         = TS01_EX01-HCT590X-SEP1230-0,8mm-0°-15mm_min-23°C
file_suffix       = .xlsx
container_archive = (blank)
-> HCT590X/thick/TS01_EX01-HCT590X-SEP1230-0,8mm-0°-15mm_min-23°C.xlsx

# archived - the remaining 236 files
file_path         = TRR285_Materials_Testing\HCT590X\base
container_archive = TS05-supplementary_data.zip
path_in_archive   = TS05-supplementary_data/TS05_AR01-HCT590X-SEP1230-1,5mm-0°-15mm_min-23°C.aramis
-> open HCT590X/base/TS05-supplementary_data.zip, then the path above
```


### YAML mapping files

The three `mapping-*.yml` files express the same test matrix in a machine-readable, ontology-linked form, one file per material. Each declares a `@context` binding the metadata field names to IRIs in the PMD Core Ontology (PMDco) and the Tensile Test Ontology (TTO) — for example `orientations_value` to `pmdco:PMD_0001035` ("angle") and `test_speed_value` to `tto:TTO_0000014` ("crosshead separation rate").

The body reflects the experimental design directly: a `defaults` block holds the material's reference condition, and each entry under `folders` names a variation group with the `deviations` from that reference plus the series it contains. A series therefore inherits every parameter it does not explicitly override.

## Naming convention

```
TS05_AR01-HCT590X-SEP1230-1,5mm-0°-15mm_min-23°C.aramis
│    │    │       │        │      │  │         └── temperature
│    │    │       │        │      │  └──────────── crosshead speed
│    │    │       │        │      └─────────────── orientation to rolling direction
│    │    │       │        └────────────────────── sheet thickness
│    │    │       └─────────────────────────────── specimen geometry / standard
│    │    └─────────────────────────────────────── material
│    └──────────────────────────────────────────── file role and index within the series
└───────────────────────────────────────────────── test series ID
```

File role codes: `EX` mechanical workbook, `AR` Aramis DIC, `CI` image sequence header, `CS` sequence companion data, `OR` per-specimen folder.

Decimal separators in file names are commas (`1,5mm`), following the original acquisition convention.

## Coverage notes

Not every condition was tested for every material. In particular:

- Sheet thickness and specimen geometry were varied for **HCT590X only**; EN-AW6014 was tested exclusively at 2.0 mm with SEP 1230.
- Elevated-temperature tests exist for **HCT590X and EN-AW6014 T4 only**; EN-AW6014 T6 was tested at room temperature only.
- DIC was recorded for six HCT590X series (`TS05`–`TS07`, `TS11`–`TS13`); high-speed camera data accompanies the 21 high-rate series across all three material states.
- The low-rate point (1.5 mm/min, 0.001 s⁻¹) exists as a single HCT590X series (`TS08`).

## License

This dataset is published under the Creative Commons Attribution-ShareAlike 4.0 International license (CC BY-SA 4.0).

## Citation

**DOI:** <https://doi.org/10.5281/zenodo.21390267>

Please cite the dataset via this DOI when using the data. The Zenodo record page provides ready-made citation strings (BibTeX, DataCite and others) with the full author list.

## Funding

This research was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) — TRR 285/2 — 418701707.
