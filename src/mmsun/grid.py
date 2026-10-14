"""Common grid + tables: the hand-off point to Phase 5 (HATS).

Products per epoch (data/aligned/<epoch>/):
  stack.fits        one image HDU per channel on a common helioprojective grid
                    (EIS FOV, default 2"/px, reference time = EIS raster mid-time)
  pixels.parquet    one row per common-grid pixel: Tx, Ty, Carrington lon/lat, mu,
                    every channel value, and dt per channel  -> HATS-ready "object" table
  native_<mod>.parquet  each modality at native sampling inside the EIS FOV, with
                    lon/lat -> enables real heliographic cross-matching between modalities
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .coalign import RSUN_ARCSEC_APPROX, Sampler, rotation_offset, world_arcsec


def make_grid(box, pix: float):
    xs = np.arange(box[0] + pix / 2, box[1], pix)
    ys = np.arange(box[2] + pix / 2, box[3], pix)
    gx, gy = np.meshgrid(xs, ys)
    return gx, gy


def grid_header(gx, gy, pix: float, t_ref: pd.Timestamp) -> dict:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.time import Time
    from sunpy.coordinates import Helioprojective
    import sunpy.map

    ref = SkyCoord(gx[0, 0] * u.arcsec, gy[0, 0] * u.arcsec,
                   frame=Helioprojective(observer="earth", obstime=Time(t_ref.to_pydatetime())))
    hdr = sunpy.map.make_fitswcs_header(np.zeros(gx.shape), ref, reference_pixel=[0, 0] * u.pix,
                                        scale=[pix, pix] * u.arcsec / u.pix)
    return dict(hdr)


def sample_on_grid(m, gx, gy, t_ref: pd.Timestamp, t_map: pd.Timestamp, d=(0.0, 0.0),
                   pix: float = 2.0) -> np.ndarray:
    """Sample map m (pointing correction d, taken at t_map) at grid positions for t_ref."""
    rx, ry = rotation_offset(float(gx.mean()), float(gy.mean()), t_ref, t_map)
    import astropy.units as u

    scale = float(np.mean(u.Quantity(m.scale).to_value(u.arcsec / u.pix)))
    blur = pix if scale < pix else 0.0          # anti-alias only when downsampling
    return Sampler(m, blur_arcsec=blur)(gx + rx - d[0], gy + ry - d[1])


def carrington(tx, ty, t_ref: pd.Timestamp):
    """Carrington lon/lat (deg) and mu for arcsec positions; NaN off-disk."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from sunpy.coordinates import HeliographicCarrington, Helioprojective

    t = Time(t_ref.to_pydatetime())
    c = SkyCoord(np.ravel(tx) * u.arcsec, np.ravel(ty) * u.arcsec,
                 frame=Helioprojective(observer="earth", obstime=t))
    hg = c.transform_to(HeliographicCarrington(observer="earth", obstime=t))
    lon = hg.lon.to_value(u.deg).reshape(np.shape(tx))
    lat = hg.lat.to_value(u.deg).reshape(np.shape(tx))
    rho = np.hypot(tx, ty) / RSUN_ARCSEC_APPROX
    mu = np.sqrt(np.clip(1 - rho**2, 0, None))
    return lon, lat, mu


def native_table(m, box, t_ref, t_map, d=(0.0, 0.0), name="value") -> pd.DataFrame:
    """Native pixels of m that fall inside box (after correction d), with lon/lat."""
    tx, ty = world_arcsec(m)
    tx, ty = tx + d[0], ty + d[1]                 # true positions at t_map
    rx, ry = rotation_offset(float(np.mean(box[:2])), float(np.mean(box[2:])), t_map, t_ref)
    tx, ty = tx + rx, ty + ry                     # rotated to t_ref
    sel = (tx >= box[0]) & (tx <= box[1]) & (ty >= box[2]) & (ty <= box[3])
    if not sel.any():
        return pd.DataFrame()
    lon, lat, mu = carrington(tx[sel], ty[sel], t_ref)
    return pd.DataFrame({"tx": tx[sel], "ty": ty[sel], "lon_carr": lon, "lat_carr": lat,
                         "mu": mu, name: np.asarray(m.data, float)[sel],
                         "t_obs": pd.Timestamp(t_map)})
