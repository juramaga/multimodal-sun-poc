"""Phase 2: choose the minimal set of files that are co-temporal and co-spatial.

"Simultaneous" is defined per instrument, relative to the EIS raster:
  * AIA pseudo-raster channel: nearest frame to *each EIS slit step* (|dt| <= 1 cadence)
  * AIA context channels: one frame per channel nearest the raster mid-time
  * HMI: one frame nearest the raster mid-time
  * XRT: frames whose FOV contains >= min_fov_coverage of the EIS FOV, with dt
    measured to the raster *interval* (0 if taken during the raster); best N
    filters by preference; ties (e.g. several frames inside the raster) go to the
    frame nearest the raster mid-time; full-disk synoptic as a fallback.
dt is always stored, so later users can tighten the tolerance without re-running:
dt_s is relative to each product's own reference (step, mid-time or interval),
dt_mid_s is always t_obs - raster mid-time.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .xrt import coverage

MATCH_COLUMNS = ["epoch", "modality", "role", "channel", "t_obs", "t_rec", "dt_s", "dt_mid_s", "steps",
                 "ds_or_url", "coverage", "est_bytes", "series"]


def _nearest(df: pd.DataFrame, t: pd.Timestamp) -> tuple[int, float]:
    dt = (df["t_obs"] - t).dt.total_seconds()
    i = int(np.argmin(np.abs(dt.values)))
    return df.index[i], float(dt.iloc[i])


def aia_pseudoraster(aia: pd.DataFrame, steps: pd.DataFrame, channel: int, tol_s: float):
    """Nearest frame of `channel` to each EIS step. Returns (step_map, frames)."""
    sub = aia[(aia.channel == channel)].dropna(subset=["t_obs"])
    if sub.empty:
        return pd.DataFrame(), pd.DataFrame()
    rows = []
    for s in steps.itertuples():
        i, dt = _nearest(sub, s.t_step)
        rows.append(dict(step=s.step, t_step=s.t_step, t_rec=sub.loc[i, "t_rec"],
                         t_obs=sub.loc[i, "t_obs"], dt_s=dt, ok=abs(dt) <= tol_s))
    step_map = pd.DataFrame(rows)
    frames = (step_map[step_map.ok].groupby("t_rec")
              .agg(t_obs=("t_obs", "first"), steps=("step", lambda x: ",".join(map(str, sorted(x)))),
                   dt_s=("dt_s", lambda x: float(np.max(np.abs(x)))))
              .reset_index())
    return step_map, frames


def nearest_per_channel(df: pd.DataFrame, t: pd.Timestamp, channels, tol_s: float) -> pd.DataFrame:
    rows = []
    for ch in channels:
        sub = df[df.channel == ch].dropna(subset=["t_obs"])
        if sub.empty:
            continue
        i, dt = _nearest(sub, t)
        if abs(dt) <= tol_s:
            r = sub.loc[i].to_dict()
            r["dt_s"] = dt
            rows.append(r)
    return pd.DataFrame(rows)


def dt_to_interval(t: pd.Series, t0: pd.Timestamp, t1: pd.Timestamp) -> pd.Series:
    before = (t0 - t).dt.total_seconds().clip(lower=0)
    after = (t - t1).dt.total_seconds().clip(lower=0)
    return np.where(t < t0, -before, after)


def select_xrt(hdr: pd.DataFrame, box, t0, t1, cfg: dict, tmid=None) -> tuple[pd.DataFrame, str]:
    """Return (chosen frames, mode) where mode is 'ar', 'synoptic' or 'none'.

    Per filter: smallest dt to the raster interval, then smallest |dt| to the mid-time
    (every frame inside the raster has interval dt = 0), then longest exposure.
    """
    tmid = t0 + (t1 - t0) / 2 if tmid is None else tmid
    if hdr is None or hdr.empty or "t_obs" not in hdr:
        return pd.DataFrame(), "none"
    tol = cfg["tolerance"]
    pref = cfg["xrt"]["preferred_filters"]
    h = hdr[(hdr.error == "") & (hdr.imtype.str.lower() == "normal")].copy()
    h = h[~h["filter"].isin(["Open", "Gband"])]
    if h.empty:
        return pd.DataFrame(), "none"
    h["coverage"] = [coverage(r, box) for r in h.itertuples()]
    h["dt_s"] = dt_to_interval(h.t_obs, t0, t1)
    h["dt_mid_s"] = (h.t_obs - tmid).dt.total_seconds()
    h["rank"] = [pref.index(f) if f in pref else len(pref) for f in h["filter"]]
    h = h[h.coverage >= tol["min_fov_coverage"]]

    def best(sub, tol_s):
        sub = sub[np.abs(sub.dt_s) <= tol_s]
        if sub.empty:
            return sub
        sub = sub.assign(adt=np.abs(sub.dt_s), adt_mid=np.abs(sub.dt_mid_s))
        sub = sub.sort_values(["rank", "adt", "adt_mid", "exptime"], ascending=[True, True, True, False])
        per_filter = sub.groupby("filter", sort=False).head(1)
        return per_filter.sort_values("rank").head(cfg["xrt"]["max_filters"])

    ar = best(h[h.fovx < 1500], tol["xrt_ar_s"])
    if not ar.empty:
        return ar, "ar"
    if cfg["xrt"].get("allow_synoptic_fallback", True):
        syn = best(h[h.fovx >= 1500], tol["xrt_synoptic_s"])
        if not syn.empty:
            return syn, "synoptic"
    return pd.DataFrame(), "none"


def est_cutout_bytes(box, pix_arcsec: float, bytes_per_px: int = 2) -> int:
    w, hgt = box[1] - box[0], box[3] - box[2]
    return int((w / pix_arcsec) * (hgt / pix_arcsec) * bytes_per_px)
