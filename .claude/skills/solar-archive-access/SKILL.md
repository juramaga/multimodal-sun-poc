---
name: solar-archive-access
description: Reference for querying and downloading SDO (AIA, HMI) from JSOC and Hinode (EIS, XRT) from NRL/MSSL/SDAC - metadata-only discovery, cutouts, time systems, and known pitfalls. Use when writing or debugging any code that searches or fetches solar data, or when choosing an access layer.
---

# Solar archive access layers (as used by this project)

Principle: **discover with metadata only, download only what was matched, and never
pull full-disk SDO frames.** Each layer below lists the cheapest metadata path first.

## SDO / JSOC (AIA, HMI) — `drms`

* Keywords only (no email, no images):
  `drms.Client().query("aia.lev1_euv_12s[2021-03-06T06:44:30Z/300s][193]", key="T_REC,T_OBS,WAVELNTH,QUALITY,EXPTIME")`
  Keep `QUALITY == 0`.
* Cutout export (email registered at http://jsoc.stanford.edu/ajax/register_email.html):
  `drms.Client(email=...).export(ds, method="url", protocol="fits", process={"im_patch": {...}})`
  with `t_ref`, `t=0` (no tracking), `r=0`, `c=0`, `locunits/boxunits="arcsec"`,
  `x, y, width, height`. See `src/mmsun/jsoc.py`.
* **Pitfalls**
  - `im_patch` can fail *silently* server-side (wrong parameter names → full frames or
    nothing). Always check the downloaded NAXIS against the requested box.
  - HMI `T_OBS` is in **TAI** (`2021.03.06_06:47:15_TAI`), 37 s ahead of UTC in 2021.
    Parse with `mmsun.jsoc.parse_jsoc_time`, never with naive string replacement.
  - `aiapy.calibrate.register` refuses cutouts; `update_pointing` and
    `correct_degradation` work on cutouts. Both fetch tables from JSOC/LMSAL.
  - Record sets: one T_REC slot with many wavelengths (`[t][94,131,...]`) keeps the
    context set simultaneous; a T_REC range with `@12s` gives the pseudo-raster.
* Alternative: `sunpy.net.Fido` with `a.jsoc.Series`, `a.jsoc.Notify(email)`,
  `a.jsoc.Cutout(bottom_left, top_right)`. Equivalent, less control over processing.

## Hinode/EIS — NRL Level-1 HDF5 (eispac)

* Files: `https://eis.nrl.navy.mil/level1/hdf5/YYYY/MM/DD/eis_YYYYMMDD_HHMMSS.{head,data}.h5`
  (mirror: `https://vsolar.mssl.ucl.ac.uk/eispac/hdf5/...`).
* Discovery: the `.head.h5` (0.1–0.5 MB) has per-step times (`times/date_obs`), pointing
  (`pointing/*`, `index/xcen,ycen`), windows (`wininfo/*`), study and title. Read it with
  `mmsun.eis.head_info` — no need for the data file.
* Search by time: `import eispac.net; Fido.search(a.Time(...), a.Instrument("EIS"), a.Level("1"))`
  (time-only). Search by target/study/pointing: the `eis_catalog` GUI or the as-run DB.
* **Pitfalls**
  - eispac does not build on Python 3.13 (dependency `pims`). Use 3.12.
  - `EISCube.sum_spectra()` ignores the mask; missing CCD pixels (−100 counts) become
    large negatives. `mmsun.eis.prepare` re-sums with the mask.
  - Pointing: maps use `index/xcen,ycen + pointing/offset_x,y`; the `pointing/solar_x,y`
    arrays use another reference (off by ~20″ here). Different windows also differ in y
    by the CCD offset (up to ~17″ between the short and long wavelength bands).
  - Raster steps are stored by x position; times decrease with column index for
    right-to-left rasters. Always pair columns with `times/date_obs` by index.
  - The anchor `eis_20210306_064444` is bundled in `eispac/data/test/`.

## Hinode/XRT — SDAC Level-1 + header range requests

* Archive: `https://umbra.nascom.nasa.gov/hinode/xrt/level1/YYYY/MM/DD/HHH00/L1_XRTYYYYMMDD_HHMMSS.s.fits`.
  File names carry time only. Filter (`EC_FW1_`, `EC_FW2_`), image type (`EC_IMTY_`),
  pointing (`XCEN`, `YCEN`), size and exposure are in the header.
* Discovery: HTTP `Range: bytes=0-(8*2880-1)` on each file until the `END` card →
  ~20–60 kB per file instead of MB (`mmsun.xrt.header_table`, cached on disk).
* VSO (`Fido.search(a.Instrument("XRT"))`) can only search by time — not by filter
  or FOV — so it is a fallback, not the discovery path.
* Synoptic program: full-disk images ~06 and ~18 UT daily (Al_mesh/Al_poly/Be_thin…),
  except during CCD bakeouts and eclipse season. Good spatial fallback, poor timing.
* Saturation: ~2500 DN; L1 values are DN, normalise by `EXPTIME`.

## Coverage / coordination registries (for choosing epochs)

* **HCR** (LMSAL HEK Coverage Registry): `https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&startTime=...&stopTime=...`
  returns IRIS and Hinode observations with pointing and FOV. The `instrument=` filter is
  not reliable — filter the JSON client-side.
* **DARTS Hinode search** (JAXA): cross-instrument plots of time and FOV.
* Hinode Operation Plans (HOPs) at ISAS list which instruments co-observed.

## Reachability

`python scripts/00_check_env.py` probes each endpoint. Some (SDAC, NRL) disallow
automated *crawling* (robots.txt) but serve direct file requests; keep request rates low
(8 parallel header requests, cached).
