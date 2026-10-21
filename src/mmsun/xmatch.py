"""Small EIS spectra / XRT image-patch catalogs, using public LSDB/HATS APIs.

Carrington angles index the solar sphere, NOT the celestial sky. CatalogStream
can consume the returned lazy crossmatch without materializing the full result.
"""

# LSDB 0.11's documented public re-exports lack explicit typing exports.
# pyright: reportPrivateImportUsage=false

from __future__ import annotations

from importlib.metadata import version
from pathlib import Path

import lsdb
import numpy as np
import pandas as pd
import pyarrow as pa

from . import eis
from .cli import load_json, save_json
from .coalign import rotation_offset, world_arcsec
from .config import ROOT, paths_for
from .grid import carrington


def _shift(shifts, key, limit):
    s = shifts[key]
    d = np.array([s["dx"], s["dy"]], dtype=float)
    if s.get("flagged", True) or not np.isfinite(d).all() or (np.abs(d) >= limit).any():
        raise ValueError(f"{key}: flagged, non-finite or boundary shift; inspect shifts.json")
    return d


def _new_output_path(path):
    path = Path(path).absolute()
    if path.is_relative_to((ROOT / "data" / "raw").resolve()) or path.resolve() != path:
        raise ValueError(f"Unsafe output path (raw data or symlink): {path}")
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Catalog output already exists: {path}; move it aside explicitly before rebuilding")
    return path


def _write_catalog(frame, metadata, path, options, margin):
    path = _new_output_path(path)
    if frame.empty:
        raise ValueError(f"No usable {metadata['modality']} samples")
    coords = frame[["lon_carr", "lat_carr"]].to_numpy()
    if not np.isfinite(coords).all():
        raise ValueError("Off-limb/non-finite coordinates are not supported by this prototype")
    types = ({"wavelength": pa.float64(), "spectrum": pa.float32(),
              "error": pa.float32(), "mask": pa.bool_()} if metadata["modality"] == "EIS"
             else {"image": pa.float32(), "image_mask": pa.bool_()})
    for column, element_type in types.items():
        frame[column] = pd.array(frame[column].tolist(), dtype=pd.ArrowDtype(pa.list_(element_type)))
    catalog = lsdb.from_dataframe(
        frame, ra_column="lon_carr", dec_column="lat_carr",
        catalog_name=metadata["modality"].lower(),
        partition_rows=options["partition_rows"], highest_order=options["highest_order"],
        margin_threshold=margin,
    )
    catalog.write_catalog(path, as_collection=True, overwrite=False, progress_bar=False)
    save_json(metadata, path / "mmsun.json")


def build_catalogs(cfg: dict, epoch: dict, *, fixture: bool = False) -> tuple[Path, Path]:
    """Export one prepared/coaligned epoch into a NEW data_dir/hats/<epoch>.

    EIS uses one spectral window and actual slit-step acquisition times. XRT
    uses one selected exposure: time eligibility can then be applied BEFORE
    spatial nearest-neighbor selection without a custom matching algorithm.
    """
    import astropy.units as u
    import eispac
    import sunpy.map

    p = paths_for(cfg, epoch["id"])
    options = cfg["xmatch"]
    size, stride = options["patch_size"], options["patch_stride"]
    radius = float(options["radius_heliographic_deg"])
    if not isinstance(size, int) or not isinstance(stride, int) or size < 1 or stride < 1:
        raise ValueError("patch_size and patch_stride must be positive integers")
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius_heliographic_deg must be finite and positive")
    shifts = load_json(p.aligned / "shifts.json")
    limit = cfg["coalign"]["max_shift_arcsec"]
    de = _shift(shifts, "EIS", limit)
    xp = pd.read_parquet(p.manifest / "products.parquet").query("modality == 'XRT'").copy()
    if xp.empty:
        raise ValueError("No prepared XRT exposure in products.parquet")
    preferences = cfg["xrt"]["preferred_filters"]
    xp["preference"] = xp["filter"].map(lambda f: preferences.index(f) if f in preferences else len(preferences))
    data = p.root / "raw" / "eis" / epoch["id"] / f"{epoch['eis_file']}.data.h5"
    info = eis.head_info(data.with_name(data.name.replace(".data.h5", ".head.h5")))
    tref = info["t_mid"]
    xp["distance"] = (pd.to_datetime(xp.t_obs) - tref).abs()
    selected = xp.sort_values(by=["preference", "distance"]).iloc[0]
    dx = _shift(shifts, f"XRT_{selected['filter']}", limit)
    xm = sunpy.map.Map(selected.path)
    if not isinstance(xm, sunpy.map.GenericMap) or xm.unit is None or not xm.unit.is_equivalent(u.DN / u.s):
        raise ValueError("Expected a single prepared XRT map in DN/s")
    txrt = pd.Timestamp(xm.date.isot)
    if not isinstance(txrt, pd.Timestamp):
        raise ValueError("XRT observation time is missing")
    if abs((txrt - pd.Timestamp(selected.t_obs)).total_seconds()) > 1:
        raise ValueError("XRT observation time disagrees with products.parquet")
    exptime = float((xm.meta or {}).get("exptime", np.nan))
    if not np.isfinite(exptime) or exptime <= 0:
        raise ValueError("XRT exposure time must be finite and positive")

    out = _new_output_path(p.root / "hats" / epoch["id"])
    cube = eispac.read_cube(data, window=cfg["eis"]["coalign_window"])
    ex, ey = world_arcsec(cube.sum_spectra(return_map=True))
    ex, ey = ex + de[0], ey + de[1]
    times = info["steps"].set_index("step")["t_step"]
    dt = np.empty(ex.shape)
    # shortcut: centre-displacement transport is only for small near-centre FOVs; use full reprojection for larger/near-limb rasters.
    for col in range(ex.shape[1]):
        t = pd.Timestamp(times.loc[col])
        if not isinstance(t, pd.Timestamp):
            raise ValueError(f"EIS column {col} has no acquisition time")
        rx, ry = rotation_offset(float(ex[:, col].mean()), float(ey[:, col].mean()), t, tref)
        ex[:, col] += rx
        ey[:, col] += ry
        dt[:, col] = (t - tref).total_seconds()
    lon, lat, _ = carrington(ex, ey, tref)
    flux = np.asarray(cube.data, dtype="float32")
    error = np.asarray(cube.uncertainty.array, dtype="float32")
    wave = np.asarray(cube.wavelength, dtype="float64")
    mask = np.zeros(flux.shape, dtype=bool) if cube.mask is None else np.asarray(cube.mask, dtype=bool).copy()
    mask |= ~np.isfinite(flux) | ~np.isfinite(error) | ~np.isfinite(wave)
    rows = []
    for y, x in np.ndindex(ex.shape):
        if mask[y, x].all():
            continue
        rows.append({"sample_id": f"{epoch['id']}:eis:{y}:{x}", "epoch": epoch["id"],
                     "ix": x, "iy": y, "lon_carr": lon[y, x], "lat_carr": lat[y, x], "dt_s": dt[y, x],
                     "wavelength": wave[y, x], "spectrum": flux[y, x],
                     "error": error[y, x], "mask": mask[y, x]})
    left = pd.DataFrame(rows)

    xx, xy = world_arcsec(xm)
    rx, ry = rotation_offset(float(xx.mean() + dx[0]), float(xy.mean() + dx[1]), txrt, tref)
    xx, xy = xx + dx[0] + rx, xy + dx[1] + ry
    image = np.asarray(xm.data, dtype="float32")
    rows = []
    for y in range(0, image.shape[0] - size + 1, stride):
        for x in range(0, image.shape[1] - size + 1, stride):
            sy, sx = slice(y, y + size), slice(x, x + size)
            px, py = xx[sy, sx], xy[sy, sx]
            if px.max() < ex.min() or px.min() > ex.max() or py.max() < ey.min() or py.min() > ey.max():
                continue
            patch = image[sy, sx]
            if not np.isfinite(patch).any():
                continue
            plon, plat, _ = carrington(np.array(px.mean()), np.array(py.mean()), tref)
            rows.append({"sample_id": f"{epoch['id']}:xrt:{selected['filter']}:{txrt.isoformat()}:{y}:{x}",
                         "epoch": epoch["id"], "ix": x, "iy": y,
                         "lon_carr": float(plon), "lat_carr": float(plat),
                         "dt_s": (txrt - tref).total_seconds(), "image": patch.ravel(),
                         "image_mask": (~np.isfinite(patch)).ravel()})
    right = pd.DataFrame(rows)
    if left.empty or right.empty:
        raise ValueError("No usable EIS spectra or XRT patches; no catalogs written")
    common = {"schema_version": 1, "epoch": epoch["id"], "t_ref": tref.isoformat(), "time_scale": "UTC",
              "coordinate_frame": "heliographic_carrington", "observer": "earth", "surface": "photosphere",
              "coordinate_unit": "deg", "source_kind": "synthetic-xrt-fixture" if fixture else "user-provided",
              "rotation_model": "FOV-centre differential-rotation displacement",
              "mask_convention": "True means invalid; no normalization or imputation",
              "radius_heliographic_deg": radius, "max_time_delta_s": cfg["tolerance"]["xrt_ar_s"],
              "versions": {name: version(name) for name in ("lsdb", "hats", "eispac", "sunpy")}}
    out.mkdir(parents=True, exist_ok=False)
    left_path, right_path = out / "eis", out / "xrt"
    _write_catalog(left, dict(common, modality="EIS", source_file=data.name,
                             window_angstrom=cfg["eis"]["coalign_window"], spectrum_shape=[flux.shape[-1]],
                             flux_unit=str(cube.unit), wavelength_unit="Angstrom", error_kind="standard_deviation",
                             pointing_correction_arcsec=de.tolist(),
                             calibration="eispac bundled pre-flight radiometric calibration"),
                   left_path, options, margin=None)
    _write_catalog(right, dict(common, modality="XRT", source_file=Path(selected.path).name,
                              filter=selected["filter"], observation_dt_s=(txrt - tref).total_seconds(),
                              image_shape=[size, size], image_order="C", flux_unit=str(xm.unit),
                              pointing_correction_arcsec=dx.tolist(), wcs_header=dict(xm.wcs.to_header())),
                   right_path, options, margin=radius * 3600 * 1.01)
    return left_path, right_path


def crossmatch_catalogs(eis_path, xrt_path, *, radius_deg=None, max_time_delta_s=None):
    """Reopen two local collections and return a lazy, time-eligible crossmatch.

    Contract: same epoch/frame/reference time, one XRT exposure. Radius is on the
    solar sphere. Arbitrary multi-exposure temporal matching is NOT implemented.
    """
    eis_path, xrt_path = Path(eis_path), Path(xrt_path)
    left_meta = load_json(eis_path / "mmsun.json")
    right_meta = load_json(xrt_path / "mmsun.json")
    for key in ("schema_version", "epoch", "t_ref", "time_scale", "coordinate_frame",
                "coordinate_unit", "observer", "surface"):
        if left_meta[key] != right_meta[key]:
            raise ValueError(f"Incompatible catalogs: {key} differs")
    if (left_meta["schema_version"] != 1 or left_meta["modality"] != "EIS"
            or right_meta["modality"] != "XRT"
            or left_meta["coordinate_frame"] != "heliographic_carrington"
            or left_meta["coordinate_unit"] != "deg" or left_meta["time_scale"] != "UTC"):
        raise ValueError("Expected version-1 EIS/XRT Carrington catalogs with UTC times")
    radius = float(radius_deg if radius_deg is not None else right_meta["radius_heliographic_deg"])
    tolerance = float(max_time_delta_s if max_time_delta_s is not None else right_meta["max_time_delta_s"])
    tx = float(right_meta["observation_dt_s"])
    if not np.isfinite([radius, tolerance, tx]).all() or radius <= 0 or tolerance < 0:
        raise ValueError("Radius must be positive and time tolerance nonnegative (both finite)")
    left, right = lsdb.open_catalog(eis_path), lsdb.open_catalog(xrt_path)
    epoch = left_meta["epoch"]
    left = left.query(f"epoch == {epoch!r} and dt_s >= {tx - tolerance!r} and dt_s <= {tx + tolerance!r}")
    right = right.query(f"epoch == {epoch!r} and dt_s == {tx!r}")
    return left.crossmatch(right, radius_arcsec=radius * 3600, n_neighbors=1,
                           require_right_margin=True, suffixes=("_eis", "_xrt"),
                           suffix_method="all_columns")
