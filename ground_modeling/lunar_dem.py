from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class DemSettings:
    xml_filename: str = "ldem_1024_90s_75s_000_030.xml"
    img_filename: str = "ldem_1024_90s_75s_000_030.img"

    crop_width_pixels: int = 3000
    crop_height_pixels: int = 3000
    crop_col_offset: int | None = None
    crop_row_offset: int | None = None

    display_nx: int = 70
    display_ny: int = 70


DEM = DemSettings()


def load_lunar_dem(
    base_dir: Path,
    settings: DemSettings = DEM,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """
    Read the PDS4 XML + IMG pair and return local X, Y, relative Z arrays.

    XML and IMG must be in base_dir. The XML label is opened; it references
    the IMG raster.
    """
    try:
        import rasterio
        from rasterio.enums import Resampling
        from rasterio.windows import Window
        from rasterio.transform import xy as transform_xy
    except ImportError as exc:
        raise ImportError(
            "rasterio/GDAL is required. Install with:\\n"
            "  conda install -c conda-forge rasterio gdal"
        ) from exc

    xml_path = base_dir / settings.xml_filename
    img_path = base_dir / settings.img_filename

    if not xml_path.exists():
        raise FileNotFoundError(f"Missing XML: {xml_path}")
    if not img_path.exists():
        raise FileNotFoundError(f"Missing IMG: {img_path}")

    with rasterio.open(xml_path) as src:
        width = min(settings.crop_width_pixels, src.width)
        height = min(settings.crop_height_pixels, src.height)

        if settings.crop_col_offset is None:
            col_off = max((src.width - width) // 2, 0)
        else:
            col_off = max(min(settings.crop_col_offset, src.width - width), 0)

        if settings.crop_row_offset is None:
            row_off = max((src.height - height) // 2, 0)
        else:
            row_off = max(min(settings.crop_row_offset, src.height - height), 0)

        window = Window(col_off, row_off, width, height)

        raw = src.read(
            1,
            window=window,
            out_shape=(settings.display_ny, settings.display_nx),
            resampling=Resampling.bilinear,
            masked=False,
        ).astype(np.float64)

        if src.nodata is not None and np.isfinite(src.nodata):
            raw[np.isclose(raw, src.nodata)] = np.nan

        scale = float(src.scales[0]) if src.scales else 1.0
        offset = float(src.offsets[0]) if src.offsets else 0.0
        z_abs = raw * scale + offset
        z_abs[~np.isfinite(z_abs)] = np.nan

        if not np.any(np.isfinite(z_abs)):
            raise RuntimeError("No valid DEM elevations were found.")

        z_abs[np.isnan(z_abs)] = np.nanmean(z_abs)

        base_transform = src.window_transform(window)
        display_transform = base_transform * rasterio.Affine.scale(
            float(window.width) / settings.display_nx,
            float(window.height) / settings.display_ny,
        )

        rows, cols = np.indices(
            (settings.display_ny, settings.display_nx)
        )
        xs, ys = transform_xy(
            display_transform,
            rows,
            cols,
            offset="center",
        )

        X = np.asarray(xs, dtype=float).reshape(
            settings.display_ny,
            settings.display_nx,
        )
        Y = np.asarray(ys, dtype=float).reshape(
            settings.display_ny,
            settings.display_nx,
        )

        X -= X[0, 0]
        Y -= Y[0, 0]
        Z = z_abs - np.nanmin(z_abs)

        info = {
            "xml": str(xml_path),
            "img": str(img_path),
            "native_shape": (src.height, src.width),
            "display_shape": Z.shape,
            "crs": str(src.crs),
            "window": (
                int(window.col_off),
                int(window.row_off),
                int(window.width),
                int(window.height),
            ),
            "scale": scale,
            "offset": offset,
            "x_range_m": float(np.nanmax(X) - np.nanmin(X)),
            "y_range_m": float(np.nanmax(Y) - np.nanmin(Y)),
            "z_range_m": float(np.nanmax(Z) - np.nanmin(Z)),
        }

    return X, Y, Z, info


class TerrainSampler:
    """Bilinear elevation sampler for a nearly rectilinear DEM grid."""

    def __init__(self, X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> None:
        if X.shape != Y.shape or X.shape != Z.shape:
            raise ValueError("X, Y and Z must have identical shapes.")

        self.X = X
        self.Y = Y
        self.Z = Z

        self.x_axis = np.asarray(X[0, :], dtype=float)
        self.y_axis = np.asarray(Y[:, 0], dtype=float)

        if self.x_axis[0] > self.x_axis[-1]:
            self.x_axis = self.x_axis[::-1]
            self.Z = self.Z[:, ::-1]

        if self.y_axis[0] > self.y_axis[-1]:
            self.y_axis = self.y_axis[::-1]
            self.Z = self.Z[::-1, :]

    @property
    def center(self) -> tuple[float, float]:
        return (
            0.5 * (float(self.x_axis[0]) + float(self.x_axis[-1])),
            0.5 * (float(self.y_axis[0]) + float(self.y_axis[-1])),
        )

    def __call__(self, x: float, y: float) -> float:
        x = float(np.clip(x, self.x_axis[0], self.x_axis[-1]))
        y = float(np.clip(y, self.y_axis[0], self.y_axis[-1]))

        ix = int(np.searchsorted(self.x_axis, x, side="right") - 1)
        iy = int(np.searchsorted(self.y_axis, y, side="right") - 1)

        ix = max(0, min(ix, len(self.x_axis) - 2))
        iy = max(0, min(iy, len(self.y_axis) - 2))

        x0, x1 = self.x_axis[ix], self.x_axis[ix + 1]
        y0, y1 = self.y_axis[iy], self.y_axis[iy + 1]

        tx = 0.0 if abs(x1 - x0) < 1e-14 else (x - x0) / (x1 - x0)
        ty = 0.0 if abs(y1 - y0) < 1e-14 else (y - y0) / (y1 - y0)

        z00 = self.Z[iy, ix]
        z10 = self.Z[iy, ix + 1]
        z01 = self.Z[iy + 1, ix]
        z11 = self.Z[iy + 1, ix + 1]

        return float(
            (1.0 - tx) * (1.0 - ty) * z00
            + tx * (1.0 - ty) * z10
            + (1.0 - tx) * ty * z01
            + tx * ty * z11
        )
