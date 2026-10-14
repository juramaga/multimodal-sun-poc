#!/usr/bin/env python
"""Build an offline fixture: the real bundled EIS anchor raster + synthetic AIA/HMI/XRT.

The synthetic SDO/XRT images contain the EIS Fe XII intensity at its *true*
position (header + EIS_SHIFT) embedded in a smooth random corona, and the XRT
header pointing is wrong by XRT_SHIFT. The offline test then checks that
scripts 04-05 recover both shifts. No network access is needed.

Usage: python tests/make_fixture.py <data_dir>
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

EIS_SHIFT = (7.5, -4.0)    # true = header + d   (arcsec)
XRT_SHIFT = (-3.0, 5.0)
EPOCH = "eis_20210306_064444"


def truth_field(box, pix=0.5, seed=1):
    """Smooth random corona over box with the EIS Fe XII image pasted at its true position."""
    import eispac
    from scipy.ndimage import gaussian_filter, map_coordinates

    from mmsun import eis

    rng = np.random.default_rng(seed)
    xs = np.arange(box[0], box[1], pix)
    ys = np.arange(box[2], box[3], pix)
    gx, gy = np.meshgrid(xs, ys)
    field = 50 + 300 * gaussian_filter(rng.random(gx.shape), 30 / pix)
    for _ in range(25):
        x0, y0 = rng.uniform(box[0], box[1]), rng.uniform(box[2], box[3])
        field += rng.uniform(100, 800) * np.exp(-((gx - x0) ** 2 + (gy - y0) ** 2) / (2 * rng.uniform(4, 15) ** 2))
    field = gaussian_filter(field, 1.0 / pix)

    head = eis.bundled_path(EPOCH, "head")
    info = eis.head_info(head)
    cube = eispac.read_cube(str(eis.bundled_path(EPOCH, "data")), 192.394)
    img = cube.sum_spectra(return_map=True)
    from mmsun.coalign import world_arcsec
    tx, ty = world_arcsec(img)
    tx, ty = tx + EIS_SHIFT[0], ty + EIS_SHIFT[1]
    data = np.nan_to_num(np.asarray(img.data, float), nan=np.nanmedian(img.data))
    # interpolate the (4"x1") EIS image onto the truth grid inside its true footprint
    x_axis, y_axis = tx[0], ty[:, 0]
    inside = (gx >= x_axis.min()) & (gx <= x_axis.max()) & (gy >= y_axis.min()) & (gy <= y_axis.max())
    fx = np.interp(gx[inside], x_axis, np.arange(len(x_axis)))
    fy = np.interp(gy[inside], y_axis, np.arange(len(y_axis)))
    vals = map_coordinates(data, [fy, fx], order=1)
    scale = np.median(field[inside]) / np.median(vals)
    field[inside] = 0.3 * field[inside] + vals * scale * 2.0
    return xs, ys, gaussian_filter(field, 1.0 / pix), info


def sample(xs, ys, field, cx, cy, nx, ny, scale, d=(0.0, 0.0)):
    """Image with header centre (cx,cy) whose pixels see true position = header + d."""
    from scipy.ndimage import map_coordinates

    hx = cx + (np.arange(nx) - (nx - 1) / 2) * scale
    hy = cy + (np.arange(ny) - (ny - 1) / 2) * scale
    HX, HY = np.meshgrid(hx, hy)
    TX, TY = HX + d[0], HY + d[1]
    fx = (TX - xs[0]) / (xs[1] - xs[0])
    fy = (TY - ys[0]) / (ys[1] - ys[0])
    return map_coordinates(field, [fy, fx], order=1, cval=np.median(field))


def header(cx, cy, nx, ny, scale, t, **extra):
    import astropy.units as u
    import sunpy.map
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from sunpy.coordinates import Helioprojective

    ref = SkyCoord(cx * u.arcsec, cy * u.arcsec,
                   frame=Helioprojective(observer="earth", obstime=Time(t.to_pydatetime())))
    h = sunpy.map.make_fitswcs_header(np.zeros((ny, nx)), ref,
                                      reference_pixel=[(nx - 1) / 2, (ny - 1) / 2] * u.pix,
                                      scale=[scale, scale] * u.arcsec / u.pix)
    h.update(extra)
    return h


def save(data, hdr, path):
    import sunpy.map

    path.parent.mkdir(parents=True, exist_ok=True)
    sunpy.map.Map(data.astype("float32"), hdr).save(path, overwrite=True)


def main(data_dir: Path):
    from mmsun import eis

    rng = np.random.default_rng(2)
    info = eis.head_info(eis.bundled_path(EPOCH, "head"))
    b = eis.fov_box(info)
    cbox = (b[0] - 150, b[1] + 150, b[2] - 150, b[3] + 150)
    xs, ys, field, _ = truth_field((cbox[0] - 60, cbox[1] + 60, cbox[2] - 60, cbox[3] + 60))
    cx, cy = (cbox[0] + cbox[1]) / 2, (cbox[2] + cbox[3]) / 2
    t0, t1 = info["date_obs"], info["date_end"]
    tmid = info["t_mid"]
    raw = data_dir / "raw"

    # EIS: real files
    eis.fetch(EPOCH, raw / "eis" / EPOCH, [], use_bundled=True)

    # AIA pseudo-raster (193) at 12 s cadence over the raster, small box
    pbox = (b[0] - 40, b[1] + 40, b[2] - 40, b[3] + 40)
    pcx, pcy = (pbox[0] + pbox[1]) / 2, (pbox[2] + pbox[3]) / 2
    pnx, pny = int((pbox[1] - pbox[0]) / 0.6), int((pbox[3] - pbox[2]) / 0.6)
    aia_kw = dict(telescop="SDO/AIA", instrume="AIA_3", waveunit="angstrom", quality=0, exptime=2.0,
                  rsun_ref=696000000.0)
    t = t0.floor("12s") - pd.Timedelta(seconds=12)
    while t <= t1 + pd.Timedelta(seconds=12):
        img = sample(xs, ys, field, pcx, pcy, pnx, pny, 0.6) * 2.0 + rng.normal(0, 5, (pny, pnx))
        h = header(pcx, pcy, pnx, pny, 0.6, t, wavelnth=193, t_obs=t.isoformat() + "Z", **aia_kw)
        save(img, h, raw / "aia" / EPOCH / "pseudoraster" / f"aia.lev1_euv_12s.{t:%Y%m%dT%H%M%S}.193.image.fits")
        t += pd.Timedelta(seconds=12)

    # AIA context, 7 channels at one slot
    cnx, cny = int((cbox[1] - cbox[0]) / 0.6), int((cbox[3] - cbox[2]) / 0.6)
    base = sample(xs, ys, field, cx, cy, cnx, cny, 0.6)
    tc = tmid.floor("12s")
    for wl, gain in [(94, 0.05), (131, 0.1), (171, 1.5), (193, 2.0), (211, 0.8), (335, 0.08), (304, 0.6)]:
        img = base * gain * 2.0 + rng.normal(0, 1, base.shape)
        h = header(cx, cy, cnx, cny, 0.6, tc, wavelnth=wl, t_obs=tc.isoformat() + "Z", **aia_kw)
        save(img, h, raw / "aia" / EPOCH / "context" / f"aia.lev1_euv_12s.{tc:%Y%m%dT%H%M%S}.{wl}.image.fits")

    # HMI magnetogram: a signed field correlated with the corona
    hnx, hny = int((cbox[1] - cbox[0]) / 0.504), int((cbox[3] - cbox[2]) / 0.504)
    bmap = sample(xs, ys, field, cx, cy, hnx, hny, 0.504)
    blos = (bmap - np.median(bmap)) * np.sign(np.sin(np.linspace(0, 3, hnx)))[None, :] * 2
    th = tmid.floor("45s")
    h = header(cx, cy, hnx, hny, 0.504, th, telescop="SDO/HMI", instrume="HMI_FRONT2",
               content="MAGNETOGRAM", bunit="Gauss", quality=0, t_obs=th.isoformat() + "Z")
    save(blos, h, raw / "hmi" / EPOCH / f"hmi.m_45s.{th:%Y%m%d_%H%M%S}_TAI.magnetogram.fits")

    # XRT: two filters, 384"x384" partial frames at 1"/px, header pointing wrong by XRT_SHIFT
    for filt, tx_off, gain in [("Al_poly", -60, 1.0), ("Be_thin", 75, 0.3)]:
        tx = tmid + pd.Timedelta(seconds=tx_off)
        n = 384
        img = sample(xs, ys, field, cx, cy, n, n, 1.0, d=XRT_SHIFT) * gain + rng.normal(0, 2, (n, n))
        h = header(cx, cy, n, n, 1.0, tx, telescop="HINODE", instrume="XRT", ec_fw1_=filt,  # both live in filter wheel 1
                   ec_fw2_="Open", ec_imty_="normal", exptime=1.0)
        save(img, h, raw / "xrt" / EPOCH / f"L1_XRT{tx:%Y%m%d_%H%M%S}.0.fits")
    # decoy: an in-raster Be_thin frame that 02_match must NOT choose (farther from the
    # mid-time, though longer exposure) and that 04_prepare must therefore ignore
    tx = tmid + pd.Timedelta(seconds=140)
    save(rng.normal(0, 2, (384, 384)), header(cx, cy, 384, 384, 1.0, tx, telescop="HINODE", instrume="XRT",
                                              ec_fw1_="Be_thin", ec_fw2_="Open", ec_imty_="normal", exptime=4.0),
         raw / "xrt" / EPOCH / f"L1_XRT{tx:%Y%m%d_%H%M%S}.0.fits")
    write_candidates(data_dir, info, b, tmid)
    print(f"fixture written to {data_dir}")


def write_candidates(data_dir: Path, info, box, tmid):
    """Synthetic Phase-1 candidate tables, shaped like jsoc.query_* / xrt.header_table output,
    so that 02_match.py and 03_download.py --dry-run are exercised offline too."""
    man = data_dir / "manifest" / EPOCH
    man.mkdir(parents=True, exist_ok=True)
    t0, t1 = info["date_obs"], info["date_end"]
    rows = []
    t = t0.floor("12s") - pd.Timedelta(seconds=24)
    while t <= t1 + pd.Timedelta(seconds=24):
        for wl in (94, 131, 171, 193, 211, 335, 304):
            q = 0 if not (wl == 193 and t == t0.floor("12s") + pd.Timedelta(seconds=60)) else 1024
            rows.append(dict(series="aia.lev1_euv_12s", t_rec=t.strftime("%Y-%m-%dT%H:%M:%SZ"),
                             t_obs=t + pd.Timedelta(seconds={94: 2, 131: 3, 171: 0, 193: 5, 211: 7,
                                                             335: 9, 304: 4}[wl]),
                             channel=wl, quality=q, exptime=2.0, cdelt=0.6))
        t += pd.Timedelta(seconds=12)
    aia = pd.DataFrame(rows)
    aia[aia.quality == 0].to_parquet(man / "cand_aia.parquet", index=False)
    th = [tmid.floor("45s") + pd.Timedelta(seconds=45 * k) for k in (-1, 0, 1)]
    pd.DataFrame(dict(series="hmi.M_45s", t_rec=[x.strftime("%Y.%m.%d_%H:%M:%S_TAI") for x in th],
                      t_obs=th, quality=0, cdelt=0.504, channel="M_45s")).to_parquet(
        man / "cand_hmi.parquet", index=False)
    cx, cy = (box[0] + box[1]) / 2, (box[2] + box[3]) / 2
    hdr = []
    for k, (filt, dt, fov, imtype) in enumerate([
        ("Al_poly", -60, 384, "normal"), ("Be_thin", 75, 384, "normal"), ("Al_poly", 20, 384, "dark"),
        ("Ti_poly", 400, 384, "normal"), ("Gband", 10, 384, "normal"), ("Al_mesh", -1500, 2048, "normal"),
        ("Al_poly", 30, 64, "normal"),               # too small to cover the EIS FOV
        ("Be_thin", 140, 384, "normal")]):           # in-raster tie with the +75 s Be_thin: must lose
        tt = tmid + pd.Timedelta(seconds=dt)
        n = fov
        hdr.append(dict(url=f"https://example.invalid/L1_XRT{tt:%Y%m%d_%H%M%S}.0.fits",
                        file=f"L1_XRT{tt:%Y%m%d_%H%M%S}.0.fits", t_obs=tt, filter=filt, imtype=imtype,
                        exptime=4.0 if dt == 140 else 1.0, naxis1=n, naxis2=n, cdelt1=1.0, cdelt2=1.0, xcen=cx, ycen=cy,
                        fovx=float(n), fovy=float(n), chip_sum=1, bytes_est=n * n * 4, error=""))
    pd.DataFrame(hdr).to_parquet(man / "cand_xrt.parquet", index=False)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
