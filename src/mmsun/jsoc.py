"""SDO (AIA, HMI) via JSOC: keyword-only discovery and server-side cutout export.

Discovery uses ``drms.Client().query`` (keywords only, no image data, no email
needed). Retrieval uses a JSOC *export* with the ``im_patch`` processing step so
that only a small box around the EIS field of view leaves Stanford, never the
4k x 4k full disk. Exports need an email registered at
http://jsoc.stanford.edu/ajax/register_email.html .
"""

from __future__ import annotations

from pathlib import Path

import drms
import pandas as pd


def _iso(t: pd.Timestamp) -> str:
    return t.strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_jsoc_time(values) -> pd.Series:
    """JSOC times come as '2021-03-06T06:47:09.57Z' (UTC) or '2021.03.06_06:47:15_TAI'.

    TAI is 37 s ahead of UTC in 2021, which matters for 12 s cadences, so TAI
    strings are converted properly instead of being parsed naively.
    """
    from astropy.time import Time

    out = []
    for v in pd.Series(values).astype(str):
        v = v.strip()
        if v.endswith("_TAI"):
            d, t = v[:-4].split("_", 1)
            out.append(Time(d.replace(".", "-") + "T" + t, scale="tai").utc.datetime)
        elif v in ("", "nan", "None", "MISSING"):
            out.append(pd.NaT)
        else:
            out.append(pd.Timestamp(v.replace("Z", "")).to_pydatetime())
    return pd.to_datetime(pd.Series(out))


def query_aia(series: str, t0: pd.Timestamp, t1: pd.Timestamp, wavelengths: list[int]) -> pd.DataFrame:
    """One row per AIA image in [t0, t1] for the listed channels (keywords only)."""
    dur = max(int((t1 - t0).total_seconds()), 12)
    wl = ",".join(str(w) for w in wavelengths)
    ds = f"{series}[{_iso(t0)}/{dur}s][{wl}]"
    df = drms.Client().query(ds, key="T_REC,T_OBS,WAVELNTH,QUALITY,EXPTIME,CDELT1")
    if df is None or len(df) == 0:
        return pd.DataFrame(columns=["series", "t_rec", "t_obs", "channel", "quality", "exptime"])
    df = df.rename(columns={"T_REC": "t_rec", "T_OBS": "t_obs", "WAVELNTH": "channel",
                            "QUALITY": "quality", "EXPTIME": "exptime", "CDELT1": "cdelt"})
    df["t_obs"] = parse_jsoc_time(df["t_obs"]).values
    df["series"] = series
    return df


def query_hmi(series: str, t0: pd.Timestamp, t1: pd.Timestamp) -> pd.DataFrame:
    dur = max(int((t1 - t0).total_seconds()), 45)
    ds = f"{series}[{_iso(t0)}/{dur}s]"
    df = drms.Client().query(ds, key="T_REC,T_OBS,QUALITY,CDELT1")
    if df is None or len(df) == 0:
        return pd.DataFrame(columns=["series", "t_rec", "t_obs", "quality"])
    df = df.rename(columns={"T_REC": "t_rec", "T_OBS": "t_obs", "QUALITY": "quality", "CDELT1": "cdelt"})
    df["t_obs"] = parse_jsoc_time(df["t_obs"]).values
    df["series"] = series
    df["channel"] = series.split(".")[1]
    return df


def export_cutout(ds: str, box: tuple[float, float, float, float], t_ref: pd.Timestamp,
                  email: str, outdir: Path, server_side_lev15: bool = False,
                  max_wait_s: int = 900) -> list[Path]:
    """Export a record set with im_patch and download it. box=(xmin,xmax,ymin,ymax) arcsec."""
    xmin, xmax, ymin, ymax = box
    process = {}
    if server_side_lev15 and ds.startswith("aia."):
        process["aia_scale_aialev1"] = {}
    process["im_patch"] = {
        "t_ref": t_ref.strftime("%Y-%m-%dT%H:%M:%S"),
        "t": 0, "r": 0, "c": 0,              # no tracking, no re-registration, no crop
        "locunits": "arcsec", "boxunits": "arcsec",
        "x": round((xmin + xmax) / 2, 2), "y": round((ymin + ymax) / 2, 2),
        "width": round(xmax - xmin, 2), "height": round(ymax - ymin, 2),
    }
    client = drms.Client(email=email)
    req = client.export(ds, method="url", protocol="fits", process=process)
    if not req.wait(timeout=max_wait_s):
        raise TimeoutError(f"JSOC export did not finish in {max_wait_s}s: {ds}")
    if not req.has_succeeded():
        raise RuntimeError(f"JSOC export failed: {ds} status={req.status}")
    outdir.mkdir(parents=True, exist_ok=True)
    res = req.download(str(outdir))
    return [Path(p) for p in res["download"].dropna()]
