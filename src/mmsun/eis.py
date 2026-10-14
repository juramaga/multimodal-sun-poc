"""Hinode/EIS: discovery from the small *.head.h5 file, retrieval, and preparation.

EIS is the anchor modality. Only Level-1 HDF5 pairs are used: the ``.head.h5``
file (~0.1-0.5 MB: pointing, per-step times, windows) is enough for discovery;
the ``.data.h5`` file (radiometrically calibratable count rates) is fetched only
for epochs that survive matching.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests

from .netutil import download


def _s(x):
    """Decode HDF5 byte strings / 1-element arrays into python scalars."""
    if isinstance(x, np.ndarray):
        x = x.ravel()[0] if x.size == 1 else x
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode().strip().rstrip("&").strip()
    if isinstance(x, np.generic):
        return x.item()
    return x


def bundled_path(name: str, kind: str = "data") -> Path | None:
    """Return the copy of an EIS file shipped inside eispac, if there is one."""
    try:
        import eispac
    except ImportError:
        return None
    p = Path(eispac.__file__).parent / "data" / "test" / f"{name}.{kind}.h5"
    return p if p.exists() else None


def remote_url(base: str, name: str, kind: str) -> str:
    # name = eis_YYYYMMDD_HHMMSS
    d = name.split("_")[1]
    return f"{base.rstrip('/')}/{d[:4]}/{d[4:6]}/{d[6:8]}/{name}.{kind}.h5"


def fetch(name: str, dest: Path, archives: list[str], kinds=("head", "data"),
          use_bundled: bool = True) -> dict[str, Path]:
    """Get the requested EIS files into ``dest``; bundled copy first, then mirrors."""
    out = {}
    for kind in kinds:
        target = dest / f"{name}.{kind}.h5"
        if target.exists() and target.stat().st_size > 0:
            out[kind] = target
            continue
        b = bundled_path(name, kind) if use_bundled else None
        if b is not None:
            shutil.copy2(b, target)
            out[kind] = target
            continue
        last_err = None
        for base in archives:
            try:
                download(remote_url(base, name, kind), target)
                out[kind] = target
                break
            except (requests.RequestException, OSError) as e:  # try next mirror
                last_err = e
        if kind not in out:
            raise RuntimeError(f"Could not fetch {name}.{kind}.h5 from {archives}: {last_err}")
    return out


def head_info(head_path: str | Path) -> dict:
    """Everything discovery and matching need, read straight from head.h5."""
    with h5py.File(head_path, "r") as f:
        idx = f["index"]
        times = pd.to_datetime([_s(t) for t in f["times/date_obs"][:]])
        solar_x = f["pointing/solar_x"][:].astype(float)
        solar_y = f["pointing/solar_y"][:].astype(float)
        wins = []
        for k in sorted(f["wininfo"].keys()):
            if not k.startswith("win"):
                continue
            g = f["wininfo"][k]
            wins.append(dict(iwin=int(_s(g["iwin"][:])), line_id=_s(g["line_id"][:]),
                             wvl_min=float(_s(g["wvl_min"][:])), wvl_max=float(_s(g["wvl_max"][:]))))
        # Pointing convention: eispac maps use index/xcen,ycen + pointing/offset_x,y
        # (NRL's EIS-to-AIA offsets). The pointing/solar_x,y arrays are in another
        # reference, so they are re-centred on the corrected centre below. Residual
        # errors of a few arcsec are what Phase 4 coalignment measures.
        off_x = float(_s(f["pointing/offset_x"][:]))
        off_y = float(_s(f["pointing/offset_y"][:]))
        xc = float(_s(idx["xcen"][:])) + off_x
        yc = float(_s(idx["ycen"][:])) + off_y
        info = dict(
            date_obs=pd.Timestamp(_s(idx["date_obs"][:])),
            date_end=pd.Timestamp(_s(idx["date_end"][:])),
            xcen=xc,
            ycen=yc,
            offset_x=off_x,
            offset_y=off_y,
            fovx=float(_s(f["pointing/fovx"][:])),
            fovy=float(_s(f["pointing/fovy"][:])),
            x_scale=float(_s(f["pointing/x_scale"][:])),
            y_scale=float(_s(f["pointing/y_scale"][:])),
            study_id=int(_s(idx["study_id"][:])),
            study=_s(idx["stud_acr"][:]),
            obstitle=_s(idx["obstitle"][:]),
            target=_s(idx["target"][:]),
            noaa=int(_s(idx["noaa_num"][:])),
            slit=_s(idx["slit_id"][:]),
            nsteps=len(times),
            windows=wins,
        )
    if not info["noaa"]:
        import re
        m = re.search(r"AR\s?(\d{5})", info["obstitle"])
        info["noaa"] = int(m.group(1)) if m else 0
    solar_x = solar_x - solar_x.mean() + info["xcen"]
    solar_y = solar_y - solar_y.mean() + info["ycen"]
    steps = pd.DataFrame({"step": np.arange(len(times)), "t_step": times, "x_step": solar_x})
    info["steps"] = steps.sort_values("t_step").reset_index(drop=True)
    info["solar_y"] = solar_y
    info["t_mid"] = info["date_obs"] + (info["date_end"] - info["date_obs"]) / 2
    return info


def fov_box(info: dict, margin: float = 0.0) -> tuple[float, float, float, float]:
    """(xmin, xmax, ymin, ymax) of the raster in arcsec, from the per-step positions."""
    xs, ys = info["steps"]["x_step"].values, info["solar_y"]
    hx = info["x_scale"] / 2
    return (xs.min() - hx - margin, xs.max() + hx + margin,
            ys.min() - 0.5 - margin, ys.max() + 0.5 + margin)


# --------------------------------------------------------------------- preparation
def prepare(data_path: Path, outdir: Path, fit_templates: list[str],
            integrate_all: bool = True, ncpu: int = 1) -> list[dict]:
    """Window-integrated intensity maps for every window + Gaussian-fit maps.

    Writes FITS maps into ``outdir`` and returns a list of product records.
    """
    import warnings

    import eispac

    warnings.filterwarnings("ignore")
    products = []
    wininfo = eispac.read_wininfo(str(data_path).replace(".data.h5", ".head.h5"))

    if integrate_all:
        for row in wininfo:
            iwin, line = int(row[0]), str(row[1])
            cube = eispac.read_cube(str(data_path), window=iwin)
            m = cube.sum_spectra(return_map=True)
            # sum_spectra ignores the mask: missing/bad CCD pixels (stored as -100 counts)
            # would enter the sum as large negatives. Re-sum with the mask applied and
            # rescale for the masked bins; drop pixels with >20% of the window masked.
            d = np.asarray(cube.data, float)
            bad = np.asarray(cube.mask, bool) if cube.mask is not None else np.zeros(d.shape, bool)
            bad |= d <= -50
            nbin = d.shape[-1]
            ngood = (~bad).sum(axis=-1)
            clean = np.where(bad, np.nan, d)
            summed = np.nansum(clean, axis=-1) * nbin / np.maximum(ngood, 1)
            summed[ngood < 0.8 * nbin] = np.nan
            m = type(m)(summed.astype("float32"), m.meta)
            tag = "eis_win%02d_%s" % (iwin, line.replace(" ", "_").replace(".", "p"))
            fn = outdir / f"{tag}.fits"
            m.save(fn, overwrite=True)
            products.append(dict(product=tag, kind="window_integrated", line=line, path=str(fn),
                                 unit=str(m.unit)))

    tdir = Path(eispac.__file__).parent / "data" / "templates"
    for t in fit_templates:
        tfile = tdir / f"{t}.template.h5"
        if not tfile.exists():
            raise FileNotFoundError(f"EIS template {t} not found in {tdir}")
        tmpl = eispac.read_template(str(tfile))
        cube = eispac.read_cube(str(data_path), tmpl.central_wave)
        fit = eispac.fit_spectra(cube, tmpl, ncpu=ncpu)
        for meas, short in [("intensity", "int"), ("velocity", "vel"), ("width", "wid")]:
            m = fit.get_map(component=0, measurement=meas)
            tag = f"eis_fit_{t.split('.')[0]}_{short}"
            fn = outdir / f"{tag}.fits"
            m.save(fn, overwrite=True)
            products.append(dict(product=tag, kind=f"fit_{meas}", line=t, path=str(fn),
                                 unit=str(m.unit)))
    return products
