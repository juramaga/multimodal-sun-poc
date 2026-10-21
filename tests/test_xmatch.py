"""Local persisted HATS, temporal eligibility, margins and real array payloads."""

# LSDB 0.11's documented public re-exports lack explicit typing exports.
# pyright: reportPrivateImportUsage=false

import subprocess
import sys
from typing import Any

import lsdb
import numpy as np
import pandas as pd
import pytest
from lsdb.streams import CatalogStream

from mmsun import eis
from mmsun.cli import load_json, save_json
from mmsun.config import ROOT, active_epochs, load_config, paths_for
from mmsun.xmatch import _write_catalog, build_catalogs, crossmatch_catalogs


def tiny_pair(tmp_path, *, margin: float | None = 40):
    # RA-like 45 degrees at the equator is a HEALPix partition boundary. These
    # are synthetic Carrington coordinates; no celestial interpretation is used.
    metadata = {"schema_version": 1, "epoch": "anchor", "t_ref": "2021-03-06T06:47:09",
                "time_scale": "UTC", "coordinate_frame": "heliographic_carrington",
                "coordinate_unit": "deg", "observer": "earth", "surface": "photosphere",
                "radius_heliographic_deg": 0.01, "max_time_delta_s": 2, "observation_dt_s": 0,
                "source_kind": "synthetic-xrt-fixture"}
    left = pd.DataFrame({"sample_id": ["a", "b", "late", "other-epoch"],
                         "epoch": ["anchor"] * 3 + ["other"],
                         "lon_carr": [44.999, 45.001, 44.999, 45.001], "lat_carr": [0.] * 4,
                         "dt_s": [0., 1., 1000., 0.], "wavelength": [[192., 193.]] * 4,
                         "spectrum": [[1., 2.]] * 4, "error": [[0.1, np.nan]] * 4,
                         "mask": [[False, True]] * 4})
    right = pd.DataFrame({"sample_id": ["valid", "closer-but-late", "closer-other-epoch"],
                          "epoch": ["anchor", "anchor", "other"],
                          "lon_carr": [45.0001, 44.999, 44.999], "lat_carr": [0.] * 3,
                          "dt_s": [0., 1000., 0.], "image": [[1., 2., 3., 4.]] * 3,
                          "image_mask": [[False, False, False, True]] * 3})
    options = {"partition_rows": 2, "highest_order": 3}
    lp, rp = tmp_path / "eis", tmp_path / "xrt"
    _write_catalog(left, dict(metadata, modality="EIS"), lp, options, None)
    _write_catalog(right, dict(metadata, modality="XRT", image_shape=[2, 2]), rp, options, margin)
    return lp, rp


def test_roundtrip_boundary_time_and_inline_payloads(tmp_path):
    left, right = tiny_pair(tmp_path)
    xr = lsdb.open_catalog(right)
    assert xr.npartitions >= 2
    assert xr.margin is not None
    assert "valid" in xr.margin.compute().sample_id.values
    matches = crossmatch_catalogs(left, right)
    chunks = list(CatalogStream(matches, partitions_per_chunk=1, shuffle=False))
    result = pd.concat(chunks).sort_values("sample_id_eis")
    assert result.sample_id_eis.tolist() == ["a", "b"]
    assert result.sample_id_xrt.tolist() == ["valid", "valid"]
    # Exact great-circle separations on this equator; the nearest eligible
    # sample is across the partition edge for 'a', not the closer late sample.
    np.testing.assert_allclose(result._dist_arcsec, [0.0011 * 3600, 0.0009 * 3600], atol=1e-6)
    assert result.iloc[0].spectrum_eis == [1., 2.]
    assert result.iloc[0].mask_eis == [False, True]
    assert np.isnan(result.iloc[0].error_eis[1])  # preserve NaN, not an Arrow null
    assert np.asarray(result.iloc[0].image_xrt).reshape(2, 2).tolist() == [[1., 2.], [3., 4.]]
    assert result.iloc[0].image_mask_xrt == [False, False, False, True]
    tighter = crossmatch_catalogs(left, right, max_time_delta_s=0).compute()
    assert tighter.sample_id_eis.tolist() == ["a"]


def test_minimal_reader_from_outside_the_repo(tmp_path):
    left, right = tiny_pair(tmp_path)
    command = [sys.executable, str(ROOT / "scripts/06_xmatch.py")]
    result = subprocess.run([*command, str(left), str(right)], cwd=tmp_path,
                            capture_output=True, text=True, check=True)
    assert "STREAM PASSED: 2 pairs" in result.stdout
    assert "spectra (1, 2), images (1, 2, 2)" in result.stdout
    assert "synthetic-xrt-fixture" in result.stdout
    missing_args = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert missing_args.returncode != 0 and "Usage:" in missing_args.stderr


def test_missing_or_insufficient_margin_is_rejected(tmp_path):
    left, right = tiny_pair(tmp_path / "missing", margin=None)
    with pytest.raises(ValueError, match="[Mm]argin"):
        crossmatch_catalogs(left, right)
    left, right = tiny_pair(tmp_path / "insufficient", margin=10)
    with pytest.raises(ValueError, match="[Mm]argin"):
        crossmatch_catalogs(left, right)


@pytest.mark.parametrize("key,value", [("epoch", "another"), ("t_ref", "2021-03-07T00:00:00"),
                                       ("coordinate_frame", "icrs"), ("coordinate_unit", "rad")])
def test_incompatible_metadata_is_rejected(tmp_path, key, value):
    left, right = tiny_pair(tmp_path)
    metadata = load_json(right / "mmsun.json")
    metadata[key] = value
    save_json(metadata, right / "mmsun.json")
    with pytest.raises(ValueError, match="Incompatible"):
        crossmatch_catalogs(left, right)


@pytest.mark.parametrize("kwargs", [{"radius_deg": 0}, {"radius_deg": np.nan},
                                    {"max_time_delta_s": -1}, {"max_time_delta_s": np.inf}])
def test_invalid_tolerances_are_rejected(tmp_path, kwargs):
    left, right = tiny_pair(tmp_path)
    with pytest.raises(ValueError, match="finite"):
        crossmatch_catalogs(left, right, **kwargs)


def test_existing_or_symlinked_outputs_are_never_overwritten(tmp_path):
    protected = tmp_path / "raw"
    protected.mkdir()
    sentinel = protected / "keep.bin"
    sentinel.write_bytes(b"original provenance")
    link = tmp_path / "linked"
    link.symlink_to(protected, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        tiny_pair(link)
    assert sorted(p.name for p in protected.iterdir()) == ["keep.bin"]
    left, _ = tiny_pair(tmp_path / "existing")
    (left / "nested-link").symlink_to(protected, target_is_directory=True)
    with pytest.raises(FileExistsError, match="already exists"):
        tiny_pair(tmp_path / "existing")
    assert sentinel.read_bytes() == b"original provenance"


def test_real_eis_spectra_and_synthetic_xrt_patch_export(tmp_path):
    import astropy.units as u
    import eispac
    import sunpy.map
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from make_fixture import header
    from sunpy.coordinates import HeliographicCarrington, Helioprojective
    from sunpy.physics.differential_rotation import solar_rotate_coordinate

    cfg = load_config()
    cfg["data_dir"] = str(tmp_path)
    epoch = active_epochs(cfg)[0]
    p = paths_for(cfg, epoch["id"])
    raw = p.raw("eis")
    eis.fetch(epoch["eis_file"], raw, [], use_bundled=True)
    data = raw / f"{epoch['eis_file']}.data.h5"
    info = eis.head_info(raw / f"{epoch['eis_file']}.head.h5")
    t = info["t_mid"] - pd.Timedelta(seconds=60)
    image = np.arange(160 * 160, dtype="float32").reshape(160, 160) / 100
    image[80, 80] = np.nan
    h = header(info["xcen"], info["ycen"], 160, 160, 1., t, instrume="XRT", telescop="HINODE",
               ec_fw1_="Al_poly", ec_fw2_="Open", exptime=1., bunit="DN / s")
    path = p.prep("xrt") / "synthetic_xrt.fits"
    native = sunpy.map.Map(image, h)
    assert isinstance(native, sunpy.map.GenericMap)
    native.save(path)
    pd.DataFrame([{"modality": "XRT", "filter": "Al_poly", "t_obs": t.isoformat(), "path": str(path)}]).to_parquet(
        p.manifest / "products.parquet")
    shifts = {"EIS": {"dx": 2., "dy": -1., "flagged": False},
              "XRT_Al_poly": {"dx": -3., "dy": 2., "flagged": False}}
    save_json(shifts, p.aligned / "shifts.json")
    left, right = build_catalogs(cfg, epoch, fixture=True)
    lm, rm = load_json(left / "mmsun.json"), load_json(right / "mmsun.json")
    assert lm["source_kind"] == rm["source_kind"] == "synthetic-xrt-fixture"
    assert rm["image_shape"] == [16, 16]
    assert lm["spectrum_shape"] == [24]
    # These reads go through persisted Parquet payloads, not original FITS/HDF5.
    cube = eispac.read_cube(data, window=cfg["eis"]["coalign_window"])
    ec = lsdb.open_catalog(left).compute()
    assert ec.dt_s.nunique() == 25
    row = ec.iloc[0]
    y, x = int(row.iy), int(row.ix)
    np.testing.assert_allclose(row.spectrum, cube.data[y, x], rtol=1e-6)
    np.testing.assert_allclose(row.wavelength, cube.wavelength[y, x])
    np.testing.assert_allclose(row.error, cube.uncertainty.array[y, x], rtol=1e-6)
    times = info["steps"].set_index("step").t_step
    assert row.dt_s == (times.loc[x] - info["t_mid"]).total_seconds()
    xr = lsdb.open_catalog(right).compute()
    patch = xr.iloc[0]
    sy, sx = slice(int(patch.iy), int(patch.iy) + 16), slice(int(patch.ix), int(patch.ix) + 16)
    np.testing.assert_allclose(np.asarray(patch.image).reshape(16, 16), image[sy, sx], equal_nan=True)
    assert any(any(mask) for mask in xr.image_mask)

    def expected_position(point, centre, correction, observed):
        # Independent SunPy transport, not the builder's rotation_offset helper.
        cx = centre.spherical.lon.to_value(u.arcsec) + correction[0]
        cy = centre.spherical.lat.to_value(u.arcsec) + correction[1]
        start = SkyCoord(cx * u.arcsec, cy * u.arcsec,
                         frame=Helioprojective(observer="earth", obstime=Time(observed)))
        # Astropy frame attributes are dynamic and absent from its type stubs.
        end: Any = solar_rotate_coordinate(start, time=Time(info["t_mid"]))
        px = point.spherical.lon.to_value(u.arcsec) + correction[0] + end.spherical.lon.to_value(u.arcsec) - cx
        py = point.spherical.lat.to_value(u.arcsec) + correction[1] + end.spherical.lat.to_value(u.arcsec) - cy
        at_ref = SkyCoord(px * u.arcsec, py * u.arcsec,
                          frame=Helioprojective(observer="earth", obstime=Time(info["t_mid"])))
        hg: Any = at_ref.transform_to(HeliographicCarrington(observer="earth", obstime=Time(info["t_mid"])))
        return hg.spherical.lon.to_value(u.deg), hg.spherical.lat.to_value(u.deg)

    base = cube.sum_spectra(return_map=True)
    for column in (0, cube.data.shape[1] - 1):
        row = ec[(ec.ix == column) & (ec.iy == 60)].iloc[0]
        expected = expected_position(base.pixel_to_world(column * u.pix, 60 * u.pix),
                                     base.pixel_to_world(column * u.pix, 59.5 * u.pix),
                                     (2., -1.), times.loc[column])
        np.testing.assert_allclose([row.lon_carr, row.lat_carr], expected, atol=1e-7, rtol=0)
    expected = expected_position(native.pixel_to_world((patch.ix + 7.5) * u.pix, (patch.iy + 7.5) * u.pix),
                                 native.pixel_to_world(79.5 * u.pix, 79.5 * u.pix), (-3., 2.), t)
    np.testing.assert_allclose([patch.lon_carr, patch.lat_carr], expected, atol=1e-7, rtol=0)
    chunks = list(CatalogStream(crossmatch_catalogs(left, right), shuffle=False))
    assert sum(map(len, chunks)) > 2000

    for key, correction in [("EIS", {"dx": 0., "dy": 0., "flagged": True}),
                            ("XRT_Al_poly", {"dx": 40., "dy": 0., "flagged": False})]:
        original = shifts[key]
        shifts[key] = correction
        save_json(shifts, p.aligned / "shifts.json")
        with pytest.raises(ValueError, match="shift"):
            build_catalogs(cfg, epoch, fixture=True)
        shifts[key] = original
    save_json(shifts, p.aligned / "shifts.json")
    existing = {path: path.read_bytes() for path in left.parent.rglob("*") if path.is_file()}
    shifts["EIS"]["dx"] += 1
    cfg["xmatch"]["patch_size"] = 9999
    save_json(shifts, p.aligned / "shifts.json")
    with pytest.raises(FileExistsError, match="already exists"):
        build_catalogs(cfg, epoch, fixture=True)
    assert all(path.read_bytes() == content for path, content in existing.items())
    alias = tmp_path / "linked-data"
    alias.symlink_to(tmp_path, target_is_directory=True)
    cfg["data_dir"] = str(alias)
    with pytest.raises(ValueError, match="symlink"):
        build_catalogs(cfg, epoch, fixture=True)
