"""Load the ELYSE site reference and preprocess the InSight accelerogram."""

import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from insight_motion import load_acceleration_csv, preprocess_acc_mseed
from insight_site_dem import read_insight_hirise_patch

HERE = Path(__file__).resolve().parents[1]
HIRISE_DEM_FILENAME = 'DTEEC_037262_1845_036761_1845_U01.IMG'
STATION_XML_FILENAME = 'stationxml.xb.elyse.2019.151.xml'
STATION_XML_NS = {'s': 'http://www.fdsn.org/xml/station/1'}


def _read_elyse_station(stationxml):
    root = ET.parse(stationxml).getroot()
    return next(
        item for item in root.findall('.//s:Station', STATION_XML_NS)
        if item.get('code') == 'ELYSE'
    )


def station_georeference(img_path, stationxml):
    """StationXML 좌표를 HiRISE 영상의 0부터 시작하는 행·열 번호로 변환한다."""
    with Path(img_path).open('rb') as stream:
        label = stream.read(65536).decode('ascii', 'replace')

    def number(key):
        match = re.search(r'^\s*' + re.escape(key) + r'\s*=\s*([-+0-9.eE]+)', label, re.M)
        if not match:
            raise ValueError(f'Missing projection metadata: {key}')
        return float(match.group(1))

    if (
        not re.search(r'MAP_PROJECTION_TYPE\s*=\s*"EQUIRECTANGULAR"', label)
        or number('MAP_PROJECTION_ROTATION') != 0
        or number('CENTER_LATITUDE') != 0
    ):
        raise ValueError('Unsupported DEM orientation')

    station = _read_elyse_station(stationxml)
    latitude = float(station.find('s:Latitude', STATION_XML_NS).text)
    longitude = float(station.find('s:Longitude', STATION_XML_NS).text)
    pixels_per_degree = number('MAP_RESOLUTION')

    # 북쪽으로 갈수록 위도는 증가하고 영상의 행 번호는 감소한다.
    row = number('LINE_PROJECTION_OFFSET') - latitude * pixels_per_degree
    column = number('SAMPLE_PROJECTION_OFFSET') + (
        longitude - number('CENTER_LONGITUDE')
    ) * pixels_per_degree
    rows = int(number('LINES'))
    columns = int(number('LINE_SAMPLES'))
    if not 0 <= row < rows - 1 or not 0 <= column < columns - 1:
        raise ValueError('Station outside DEM')

    return dict(
        station='ELYSE',
        latitude_deg=latitude,
        longitude_deg=longitude,
        row_0based=row,
        col_0based=column,
        projection='north-up planetocentric equirectangular',
        pixel_scale_m=number('MAP_SCALE'),
        caveat='Coordinate conversion, not independent DEM/station datum co-registration verification',
    )


def insight_reference_site(width_m, nx=31, ny=31, root=None):
    """ELYSE를 평탄 모델의 원점으로 삼고 실제 DEM 표본을 함께 보존한다."""
    root = Path(root) if root is not None else HERE
    img = root / HIRISE_DEM_FILENAME
    xml = root / STATION_XML_FILENAME
    site = station_georeference(img, xml)
    center_row = site['row_0based']
    center_col = site['col_0based']
    east_m, north_m, _, dem_height = read_insight_hirise_patch(
        img, width_m=width_m, height_m=width_m, nx=nx, ny=ny,
        center_row=center_row, center_col=center_col,
    )
    _, _, _, center_patch = read_insight_hirise_patch(
        img, width_m=2.0, height_m=2.0, nx=3, ny=3,
        center_row=center_row, center_col=center_col,
    )
    elevation = float(center_patch[1, 1])
    station = _read_elyse_station(xml)
    station_elevation = float(station.find('s:Elevation', STATION_XML_NS).text)
    site.update(
        kind='insight_station_centered_flat_reference',
        width_m=width_m,
        local_origin='ELYSE StationXML coordinates; x East, y North, z Up',
        local_z0_dem_elevation_m=elevation,
        station_metadata_elevation_m=station_elevation,
        dem_minus_station_elevation_m=elevation - station_elevation,
        dem_elevation_range_m=[float(dem_height.min()), float(dem_height.max())],
        max_flattening_difference_m=float(np.max(abs(dem_height - elevation))),
        terrain_mode='horizontal reference at station DEM elevation; actual relief sampled but not applied',
        geology_scope='InSight landing-site reference stratigraphy inferred from geological literature',
        source_files={str(path): file_sha256(path) for path in (img, xml)},
        coordinate_caveat=(
            'Station metadata anchor near InSight lander, not an independent lander-to-SEIS survey; '
            'DEM vertical disagreement retained'
        ),
    )
    samples = np.column_stack((
        east_m.ravel(), north_m.ravel(), dem_height.ravel(),
        (dem_height - elevation).ravel(), np.zeros(dem_height.size),
    ))
    return site, samples


def inspect_site():
    """위치 점검 모드에서 StationXML과 DEM의 표고 차이를 보고한다."""
    xml = HERE / STATION_XML_FILENAME
    img = HERE / HIRISE_DEM_FILENAME
    site = station_georeference(img, xml)
    _, _, _, dem_height = read_insight_hirise_patch(
        img, width_m=60, height_m=60, nx=5, ny=5,
        center_row=site['row_0based'], center_col=site['col_0based'],
    )
    station = _read_elyse_station(xml)
    elevation = float(station.find('s:Elevation', STATION_XML_NS).text)
    site.update(
        station_metadata_elevation_m=elevation,
        dem_at_station_m=float(dem_height[2, 2]),
        dem_minus_station_elevation_m=float(dem_height[2, 2] - elevation),
        dem_patch_range_m=[float(dem_height.min()), float(dem_height.max())],
        caveat=(
            'Station-centred pixel conversion works; horizontal/vertical datum or DEM accuracy '
            'still needs external verification. Synthetic coupled model does not use this DEM.'
        ),
    )
    return site


def prepare_motion(
    source, out, *, duration=2.0, warmup=5.0, start=None,
    freqmax=1.0, rate=10.0, allow_short_warmup=False,
):
    """전체 관측파를 먼저 전처리한 뒤 선행 구간과 해석 대상 구간을 자른다."""
    out = Path(out)
    full = preprocess_acc_mseed(
        source, out / 'motion_full_enu.csv', freqmax=freqmax, target_rate=rate
    )
    dt, acceleration = load_acceleration_csv(full)
    if not np.isfinite([duration, warmup]).all() or duration < 2 * dt or warmup < 0:
        raise ValueError('Invalid motion duration/warmup')

    peak = int(np.argmax(np.linalg.norm(acceleration, axis=1)))
    count = int(round(duration / dt))
    if start is None:
        # 시작 시각이 없으면 가속도 벡터 최대값을 중심으로 구간을 잡는다.
        target = max(0, min(peak - count // 2, len(acceleration) - count - 1))
    else:
        target = int(round(start / dt))
        if not math.isfinite(start) or start < 0:
            raise ValueError('Invalid source target start')
    end = target + count
    if target < 0 or end >= len(acceleration):
        raise ValueError('Requested target interval outside record; no silent truncation')

    requested_first = target - int(round(warmup / dt))
    if requested_first < 0 and not allow_short_warmup:
        raise ValueError(
            'Requested warmup precedes record start; explicitly allow_short_warmup to truncate'
        )
    first = max(0, requested_first)
    selected = acceleration[first:end + 1].copy()
    warmup_samples = target - first

    # 입력 시작 충격을 줄이는 테이퍼는 선행 구간에만 적용한다.
    taper_samples = min(warmup_samples, int(round(min(1.0, warmup) / dt)))
    if taper_samples >= 2:
        selected[:taper_samples] *= (
            0.5 * (1 - np.cos(np.linspace(0, np.pi, taper_samples)))
        )[:, None]

    metadata = json.loads(full.with_suffix('.json').read_text())
    metadata.update(
        source_sha256=file_sha256(source),
        source_peak_s=peak * dt,
        source_start_s=first * dt,
        target_start_source_s=target * dt,
        target_end_source_s=end * dt,
        target_start_model_s=warmup_samples * dt,
        effective_warmup_s=warmup_samples * dt,
        requested_warmup_s=warmup,
        requested_duration_s=duration,
        effective_duration_s=count * dt,
        requested_target_start_s=start,
        warmup_truncated=requested_first < 0,
        allow_short_warmup=allow_short_warmup,
        selection='peak_acceleration_window_NOT_global_max_displacement',
        target_retapered=False,
        zero_warmup_warning=warmup_samples == 0,
        units_verification=(
            'Supplied ACC product assumed m/s2; original product calibration history still external'
        ),
        history_caveat='Truncated preceding history; warmup sensitivity required for nonlinear assessment',
    )
    window_path = out / 'motion_window_enu.csv'
    np.savetxt(
        window_path,
        np.column_stack((np.arange(len(selected)) * dt, selected)),
        delimiter=',',
        header='time_s,east_m_s2,north_m_s2,up_m_s2',
        comments='',
    )
    window_path.with_suffix('.json').write_text(json.dumps(metadata, indent=2))
    return dt, selected, metadata


def file_sha256(path):
    """대형 DEM에도 같은 방식으로 해시를 매길 수 있도록 나눠 읽는다."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()
