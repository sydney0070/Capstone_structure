"""Windowed reader for attached-label, north-up HiRISE PDS3 float32 DTMs.

The default patch is the image centre, NOT an asserted SEIS station location.
Coordinates are local East/North/Up metres; original DEM heights are retained.
"""
from pathlib import Path
import re
import numpy as np


def read_insight_hirise_patch(path, width_m=150., height_m=150., nx=13, ny=13,
                             center_row=None, center_col=None):
    path = Path(path)
    with path.open('rb') as f:
        label = f.read(65536).decode('ascii', errors='replace')
    label = label.split('\nEND\r')[0]
    def value(key):
        match = re.search(r'^\s*' + re.escape(key) + r'\s*=\s*([^\r\n]+)', label, re.M)
        if not match:
            raise ValueError(f'PDS3 label field missing: {key}')
        return match.group(1).strip().strip('"')
    def number(key):
        return float(value(key).split()[0])
    if value('SAMPLE_TYPE') != 'PC_REAL' or number('SAMPLE_BITS') != 32 or number('BANDS') != 1:
        raise ValueError('Only single-band PC_REAL float32 HiRISE DTM is supported')
    if value('MAP_PROJECTION_TYPE') != 'EQUIRECTANGULAR' or number('MAP_PROJECTION_ROTATION') != 0:
        raise ValueError('Expected north-up equirectangular DEM')
    rows, cols = int(number('LINES')), int(number('LINE_SAMPLES'))
    spacing = number('MAP_SCALE')
    offset = (int(number('^IMAGE')) - 1) * int(number('RECORD_BYTES'))
    if path.stat().st_size < offset + rows * cols * 4:
        raise ValueError('Truncated IMG raster')
    if nx < 2 or ny < 2 or not np.isfinite([width_m, height_m]).all() or min(width_m, height_m) <= 0:
        raise ValueError('Positive patch dimensions and nx/ny >= 2 required')
    row = (rows-1)/2 if center_row is None else float(center_row)
    col = (cols-1)/2 if center_col is None else float(center_col)
    x = np.linspace(-width_m/2, width_m/2, nx)
    y = np.linspace(-height_m/2, height_m/2, ny)
    cc, rr = np.meshgrid(col+x/spacing, row-y/spacing)
    if not np.isfinite([row, col]).all() or rr.min()<0 or cc.min()<0 or rr.max()>=rows-1 or cc.max()>=cols-1:
        raise ValueError('Requested DEM patch exceeds raster boundary; change centre/size')
    raster = np.memmap(path, dtype='<f4', mode='r', offset=offset, shape=(rows, cols))
    r0, c0 = np.floor(rr).astype(int), np.floor(cc).astype(int)
    dr, dc = rr-r0, cc-c0
    corners = [np.asarray(raster[r0+di, c0+dj], dtype=float) for di,dj in [(0,0),(0,1),(1,0),(1,1)]]
    lower, upper = number('VALID_MINIMUM'), number('VALID_MAXIMUM')
    if any(np.any(~np.isfinite(a) | (a < lower-0.01) | (a > upper+0.01)) for a in corners):
        raise ValueError('DEM patch contains missing/invalid elevations; choose another centre')
    z = ((1-dr)*((1-dc)*corners[0]+dc*corners[1]) + dr*((1-dc)*corners[2]+dc*corners[3]))
    z = z * number('SCALING_FACTOR') + number('OFFSET')
    X, Y = np.meshgrid(x, y)
    print(f'HiRISE: {rows}x{cols}, native spacing {spacing:.6f} m; patch centre row={row:g}, col={col:g} (0-based)')
    return X, Y, z-float(z.mean()), z
