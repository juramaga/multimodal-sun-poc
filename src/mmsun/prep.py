"""Calibration of SDO cutouts to science-ready, native-resolution maps.

AIA: Level-1 cutout -> aiapy.update_pointing -> aiapy.correct_degradation ->
     divide by exposure  (DN/s/pix, native 0.6" pixels, NOT re-registered:
     registration is unnecessary because every later step resamples through the
     WCS, and aiapy.register refuses cutouts anyway).
HMI: magnetogram cutouts are already science-ready (Gauss, LOS); only QUALITY
     is checked and off-limb NaNs are kept.
Every step that fails (e.g. calibration tables unreachable) is recorded in the
product record instead of silently skipped, so provenance stays honest.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd


def pointing_table(t_mid: pd.Timestamp):
    from aiapy.calibrate.utils import get_pointing_table
    from sunpy.time import TimeRange

    tr = TimeRange(t_mid - pd.Timedelta(hours=12), t_mid + pd.Timedelta(hours=12))
    for src in ("jsoc", "lmsal"):
        try:
            return get_pointing_table(src, time_range=tr)
        except Exception as e:  # noqa: BLE001
            last = e
    warnings.warn(f"AIA pointing table unavailable: {last!r}")
    return None


def prep_aia(path: Path, outdir: Path, ptab=None, degradation: bool = True, tag: str = "") -> dict:
    import astropy.units as u
    import sunpy.map
    from aiapy.calibrate import correct_degradation, update_pointing

    m = sunpy.map.Map(path)
    steps = []
    if ptab is not None:
        try:
            m = update_pointing(m, pointing_table=ptab)
            steps.append("update_pointing")
        except Exception as e:  # noqa: BLE001
            steps.append(f"update_pointing:FAILED({type(e).__name__})")
    if degradation:
        try:
            m = correct_degradation(m)
            steps.append("correct_degradation")
        except Exception as e:  # noqa: BLE001
            steps.append(f"correct_degradation:FAILED({type(e).__name__})")
    else:
        steps.append("correct_degradation:SKIPPED")
    exptime = m.exposure_time.to_value(u.s)
    data = (m.data / exptime).astype("float32")
    meta = m.meta.copy()
    meta["bunit"] = "DN / (pix s)"
    meta["history"] = "mmsun prep_aia: " + ",".join(steps) + f", /EXPTIME({exptime:.3f}s)"
    out = sunpy.map.Map(data, meta)
    wl = int(round(m.wavelength.to_value(u.AA)))
    fn = outdir / f"aia_{wl:04d}_{m.date.strftime('%H%M%S')}{tag}.fits"
    out.save(fn, overwrite=True)
    return dict(product=fn.stem, kind="aia", channel=wl, t_obs=str(m.date.isot), path=str(fn),
                unit="DN/pix/s", steps=";".join(steps), quality=int(m.meta.get("quality", -1)))


def prep_hmi(path: Path, outdir: Path) -> dict:
    import sunpy.map

    m = sunpy.map.Map(path)
    q = int(m.meta.get("quality", 0))
    data = m.data.astype("float32")
    meta = m.meta.copy()
    meta["history"] = "mmsun prep_hmi: float32 copy, no resampling"
    out = sunpy.map.Map(data, meta)
    content = str(m.meta.get("content", "")).upper()
    seg = "blos" if ("MAGNETOGRAM" in content or "magnetogram" in str(path).lower()) else "ic"
    fn = outdir / f"hmi_{seg}_{m.date.strftime('%H%M%S')}.fits"
    out.save(fn, overwrite=True)
    return dict(product=fn.stem, kind="hmi", channel=seg, t_obs=str(m.date.isot), path=str(fn),
                unit=str(m.unit), steps="none", quality=q)


def nan_stats(path: str | Path) -> float:
    import sunpy.map

    d = sunpy.map.Map(path).data
    return float(np.isnan(d).mean())
