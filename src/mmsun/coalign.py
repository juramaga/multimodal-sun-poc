"""Phase 4: coalignment by sampling one map at another's pixel positions.

All instruments are compared in helioprojective arcsec. A candidate shift
(dx, dy) is scored by the Pearson correlation of log-intensities; the search is
coarse-to-fine (4" -> 1" -> 0.25"). Conventions used throughout:

  * a pointing correction d=(dx, dy) means   true position = header position + d
  * a map taken at t_m is sampled for reference time t_ref at
        header coord = X + r(t_ref -> t_m) - d
    where r is the solar-rotation displacement of a feature at the FOV centre.

AIA is the reference (d = 0). EIS is matched to an AIA 193 *pseudo-raster*
(one AIA frame per EIS slit step), XRT to the best of AIA 335/94/211.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates

RSUN_ARCSEC_APPROX = 965.0


# ------------------------------------------------------------------ sampling helpers
def _world_unit(m):
    import astropy.units as u

    return u.Unit(m.wcs.world_axis_units[0])


def world_arcsec(m):
    """(Tx, Ty) in arcsec of every pixel centre of map m, shape = m.data.shape."""
    import astropy.units as u

    ny, nx = m.data.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    wx, wy = m.wcs.pixel_to_world_values(xx, yy)
    unit = _world_unit(m)
    wx, wy = (wx * unit).to_value(u.deg), (wy * unit).to_value(u.deg)
    wx = (wx + 180.0) % 360.0 - 180.0          # HPLN can come back wrapped to (-360, 360)
    return wx * 3600.0, wy * 3600.0


class Sampler:
    """Bilinear sampling of a map at arbitrary arcsec positions (NaN-safe, optional blur)."""

    def __init__(self, m, blur_arcsec: float = 0.0, log: bool = False):
        import astropy.units as u

        self.map = m
        d = np.array(m.data, dtype=float)
        self.valid = np.isfinite(d)
        fill = np.nanmedian(d) if self.valid.any() else 0.0
        d = np.where(self.valid, d, fill)
        if blur_arcsec > 0:
            scale = float(np.mean(u.Quantity(m.scale).to_value(u.arcsec / u.pix)))
            sig = blur_arcsec / 2.355 / scale
            if sig > 0.3:
                d = gaussian_filter(d, sig)
        if log:
            floor = max(float(np.percentile(d, 1)), 1e-6)
            d = np.log10(np.clip(d, floor, None))
        self.data = d
        self.unit = _world_unit(m)

    def __call__(self, tx, ty):
        import astropy.units as u

        px, py = self.map.wcs.world_to_pixel_values((np.asarray(tx) * u.arcsec).to_value(self.unit),
                                                    (np.asarray(ty) * u.arcsec).to_value(self.unit))
        ny, nx = self.data.shape
        px, py = np.ravel(px), np.ravel(py)
        # positions within half a pixel of the outermost centres still belong to the image
        px = np.where((px >= -0.5) & (px < 0), 0, np.where((px > nx - 1) & (px <= nx - 0.5), nx - 1, px))
        py = np.where((py >= -0.5) & (py < 0), 0, np.where((py > ny - 1) & (py <= ny - 0.5), ny - 1, py))
        vals = map_coordinates(self.data, [np.ravel(py), np.ravel(px)], order=1, cval=np.nan)
        ok = map_coordinates(self.valid.astype(float), [np.ravel(py), np.ravel(px)], order=0, cval=0)
        vals[ok < 0.5] = np.nan
        return vals.reshape(np.shape(tx))


def rotation_offset(tx: float, ty: float, t_from: pd.Timestamp, t_to: pd.Timestamp) -> tuple[float, float]:
    """Displacement (arcsec) of an on-disk feature at (tx, ty) between two times."""
    if abs((t_to - t_from).total_seconds()) < 1:
        return 0.0, 0.0
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from sunpy.coordinates import Helioprojective
    from sunpy.physics.differential_rotation import solar_rotate_coordinate

    c = SkyCoord(tx * u.arcsec, ty * u.arcsec,
                 frame=Helioprojective(observer="earth", obstime=Time(t_from.to_pydatetime())))
    r = solar_rotate_coordinate(c, time=Time(t_to.to_pydatetime()))
    return float(r.Tx.to_value(u.arcsec) - tx), float(r.Ty.to_value(u.arcsec) - ty)


def _corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 50:
        return -1.0
    a, b = a[m] - a[m].mean(), b[m] - b[m].mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else -1.0


def search_shift(target: np.ndarray, sample_fn, max_shift: float = 40.0) -> dict:
    """Maximise corr(target, sample_fn(dx, dy)) over shifts; coarse-to-fine."""
    best = (0.0, 0.0, -2.0)
    surface = []
    for step, half in ((4.0, max_shift), (1.0, 4.0), (0.25, 1.0)):
        cx, cy = (0.0, 0.0) if step == 4.0 else best[:2]
        grid = np.arange(-half, half + 1e-9, step)
        for dx in cx + grid:
            for dy in cy + grid:
                r = _corr(target, sample_fn(dx, dy))
                if step == 4.0:
                    surface.append((dx, dy, r))
                if r > best[2]:
                    best = (float(dx), float(dy), r)
    return dict(dx=best[0], dy=best[1], r=best[2], surface=np.array(surface))


# ------------------------------------------------------------------ EIS vs AIA
def eis_vs_aia_pseudoraster(eis_map, step_times: pd.DataFrame, aia_frames: list,
                            max_shift: float = 40.0, blur_arcsec: float = 3.0) -> dict:
    """eis_map columns == EIS steps (step index = column index).

    aia_frames: list of (t_obs, sunpy map) for the pseudo-raster channel.
    The returned shift (dx, dy) is the EIS pointing correction.
    """
    tx, ty = world_arcsec(eis_map)
    target = np.log10(np.clip(np.array(eis_map.data, float), 1e-3, None))
    target[~np.isfinite(np.array(eis_map.data, float))] = np.nan
    times = pd.DatetimeIndex([pd.Timestamp(t) for t, _ in aia_frames])
    samplers = [Sampler(m, blur_arcsec=blur_arcsec, log=True) for _, m in aia_frames]
    st = step_times.set_index("step")["t_step"]
    col_frame = np.array([int(np.argmin(np.abs((times - st.loc[i]).total_seconds())))
                          for i in range(eis_map.data.shape[1])])

    def pseudo(dx, dy):
        out = np.full(target.shape, np.nan)
        for k in np.unique(col_frame):
            cols = np.where(col_frame == k)[0]
            # EIS pixel with header coord h shows the Sun at h + d -> sample AIA there
            out[:, cols] = samplers[k](tx[:, cols] + dx, ty[:, cols] + dy)
        return out

    res = search_shift(target, pseudo, max_shift)
    res["pseudoraster"] = pseudo(res["dx"], res["dy"])
    res["pseudoraster_unshifted"] = pseudo(0.0, 0.0)
    res["r_before"] = _corr(target, res["pseudoraster_unshifted"])
    res["column_frame"] = col_frame
    return res


# ------------------------------------------------------------------ generic map vs AIA
def map_vs_reference(img_map, ref_map, box, t_ref: pd.Timestamp, t_img: pd.Timestamp,
                     pix: float = 2.0, max_shift: float = 40.0, blur_arcsec: float = 3.0) -> dict:
    """Shift of img_map (e.g. XRT) relative to ref_map (AIA) inside box=(x0,x1,y0,y1)."""
    xs = np.arange(box[0], box[1], pix)
    ys = np.arange(box[2], box[3], pix)
    gx, gy = np.meshgrid(xs, ys)
    target = Sampler(ref_map, blur_arcsec=blur_arcsec, log=True)(gx, gy)
    rx, ry = rotation_offset(float(xs.mean()), float(ys.mean()), t_ref, t_img)
    smp = Sampler(img_map, blur_arcsec=blur_arcsec, log=True)
    res = search_shift(target, lambda dx, dy: smp(gx + rx - dx, gy + ry - dy), max_shift)
    res["rotation_offset"] = (rx, ry)
    res["r_before"] = _corr(target, smp(gx + rx, gy + ry))
    return res
