"""Hinode/XRT: header-only discovery at the SDAC Level-1 archive, selection, prep.

The SDAC archive (one directory per UT hour) does not expose filter or pointing
in file names, and VSO cannot search on them. So discovery reads only the FITS
primary header of each candidate with an HTTP ``Range`` request (tens of kB per
file instead of MB), caches it on disk, and decides from the header which files
are worth downloading. This is the approach used by the ``hinode`` package
(sun-data/hinode); it is reimplemented here to keep the dependency surface small.
"""

from __future__ import annotations

import hashlib
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from .netutil import get

FILE_RE = re.compile(r'href="(L1_XRT(\d{8})_(\d{6}\.\d)\.fits)"')
BLOCK, CARD = 2880, 80


def hour_dirs(t0: pd.Timestamp, t1: pd.Timestamp, base: str) -> list[str]:
    h = t0.floor("h")
    out = []
    while h <= t1:
        out.append(f"{base.rstrip('/')}/{h:%Y/%m/%d}/H{h:%H}00/")
        h += pd.Timedelta(hours=1)
    return out


def list_files(t0: pd.Timestamp, t1: pd.Timestamp, base: str) -> pd.DataFrame:
    rows = []
    for d in hour_dirs(t0, t1, base):
        try:
            txt = get(d, timeout=60).text
        except Exception as e:  # missing hour (bakeout, eclipse) is normal
            print(f"  [xrt] no listing for {d}: {type(e).__name__}")
            continue
        for name, date, clock in dict.fromkeys(FILE_RE.findall(txt)):
            t = pd.Timestamp(f"{date[:4]}-{date[4:6]}-{date[6:]}T{clock[:2]}:{clock[2:4]}:{clock[4:]}")
            rows.append(dict(url=d + name, file=name, t_name=t))
    df = pd.DataFrame(rows, columns=["url", "file", "t_name"])
    return df[(df.t_name >= t0) & (df.t_name <= t1)].reset_index(drop=True)


def _read_header(url: str, cache: Path) -> fits.Header:
    key = cache / (hashlib.sha1(url.encode()).hexdigest() + ".hdr")
    if key.exists():
        return fits.Header.fromstring(key.read_text())
    nblk = 8
    while True:
        r = get(url, headers={"Range": f"bytes=0-{nblk * BLOCK - 1}"}, timeout=60)
        content = r.content
        for s in range(0, len(content) - CARD + 1, CARD):
            if content[s:s + CARD].rstrip() == b"END":
                text = content[: s + CARD].decode("ascii", errors="replace")
                key.write_text(text)
                return fits.Header.fromstring(text)
        if len(content) < nblk * BLOCK:
            raise ValueError(f"no END card in header of {url}")
        nblk *= 2


def _filter_name(h: fits.Header) -> str:
    names = [str(h.get(k, "Open")).strip() for k in ("EC_FW1_", "EC_FW2_")]
    names = [n for n in names if n and n != "Open"]
    return "/".join(names) if names else "Open"


def header_table(files: pd.DataFrame, cache: Path, max_requests: int = 400,
                 workers: int = 8) -> pd.DataFrame:
    cache.mkdir(parents=True, exist_ok=True)
    files = files.head(max_requests)

    def one(row):
        try:
            h = _read_header(row.url, cache)
        except Exception as e:
            return dict(url=row.url, error=repr(e))
        nx, ny = int(h.get("NAXIS1", 0)), int(h.get("NAXIS2", 0))
        dx, dy = float(h.get("CDELT1", np.nan)), float(h.get("CDELT2", np.nan))
        return dict(
            url=row.url, file=row.file, t_obs=pd.Timestamp(h.get("DATE_OBS", row.t_name)),
            filter=_filter_name(h), imtype=str(h.get("EC_IMTY_", "")).strip(),
            exptime=float(h.get("EXPTIME", np.nan)), naxis1=nx, naxis2=ny,
            cdelt1=dx, cdelt2=dy, xcen=float(h.get("XCEN", np.nan)), ycen=float(h.get("YCEN", np.nan)),
            fovx=nx * dx, fovy=ny * dy, chip_sum=int(h.get("CHIP_SUM", 1) or 1),
            bytes_est=nx * ny * abs(int(h.get("BITPIX", 16))) // 8, error="",
        )

    with ThreadPoolExecutor(workers) as ex:
        rows = list(ex.map(one, files.itertuples()))
    return pd.DataFrame(rows)


def coverage(row, box) -> float:
    """Fraction of box=(xmin,xmax,ymin,ymax) [arcsec] inside the XRT frame."""
    fx0, fx1 = row.xcen - row.fovx / 2, row.xcen + row.fovx / 2
    fy0, fy1 = row.ycen - row.fovy / 2, row.ycen + row.fovy / 2
    ox = max(0.0, min(fx1, box[1]) - max(fx0, box[0]))
    oy = max(0.0, min(fy1, box[3]) - max(fy0, box[2]))
    area = (box[1] - box[0]) * (box[3] - box[2])
    return float(ox * oy / area) if area > 0 else 0.0


def prepare(path: Path, outdir: Path) -> dict:
    """Level-1 -> DN/s map with a saturation mask; keeps native pixels."""
    import sunpy.map

    m = sunpy.map.Map(path)
    exptime = float(m.meta.get("exptime", np.nan))
    data = m.data.astype("float32")
    sat = data >= 2500  # XRT CCD saturates near 2500 DN (L1 values are DN)
    data = np.where(sat, np.nan, data) / exptime
    meta = m.meta.copy()
    meta["bunit"] = "DN / s"
    meta["history"] = f"mmsun: divided by EXPTIME={exptime}; saturated px (>=2500 DN) set NaN"
    out = sunpy.map.Map(data, meta)
    filt = _filter_name(fits.Header(dict((k.upper(), v) for k, v in m.meta.items()
                                         if k.upper() in ("EC_FW1_", "EC_FW2_"))))
    fn = outdir / f"xrt_{filt.replace('/', '+')}_{m.date.strftime('%H%M%S')}.fits"
    out.save(fn, overwrite=True)
    return dict(product=fn.stem, kind="xrt", filter=filt, channel=filt, t_obs=str(m.date.isot),
                path=str(fn), unit="DN / s",
                n_saturated=int(sat.sum()))
