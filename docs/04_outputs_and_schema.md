# 04 · Outputs and table schemas

```
data/
├── manifest/<epoch>/
│   ├── epoch.json            EIS raster summary (times, pointing, FOV box, windows, study, NOAA)
│   ├── eis_steps.parquet     step, t_step, x_step
│   ├── cand_aia.parquet      series, t_rec, t_obs (UTC), channel, quality, exptime, cdelt
│   ├── cand_hmi.parquet      series, t_rec, t_obs (UTC), channel, quality, cdelt
│   ├── cand_xrt.parquet      url, file, t_obs, filter, imtype, exptime, naxis*, cdelt*, xcen, ycen, fovx, fovy, coverage, error
│   ├── match.parquet         see below
│   ├── eis_step_map.parquet  step, t_step, t_rec, t_obs, dt_s, ok  (EIS step → AIA frame)
│   ├── match_report.json     per-modality verdict, n_files, est_mb
│   ├── download_log.json     files fetched, raw bytes
│   └── products.parquet      one row per calibrated map: modality, role, channel/line, t_obs, path, unit, steps (calibration provenance)
├── raw/{eis,aia,hmi,xrt}/<epoch>/   untouched downloads (provenance; never edited)
├── prep/{eis,aia,hmi,xrt}/<epoch>/  calibrated native-resolution FITS maps
├── aligned/<epoch>/
│   ├── shifts.json           pointing corrections and correlation diagnostics
│   ├── coalign_diag.npz      EIS target, AIA pseudo-raster before/after, correlation surfaces
│   ├── stack.fits            one ImageHDU per channel on the common grid (EXTNAME, BUNIT, DT_S)
│   ├── pixels.parquet        common-grid table (one row per pixel)
│   └── native.parquet        native-resolution long table (one row per native pixel and channel)
└── cache/                    XRT header cache
```

## `match.parquet`
`epoch, modality, role (anchor|pseudoraster|context|ar|synoptic), channel, t_obs, t_rec,
dt_s, dt_mid_s, steps, ds_or_url, coverage, est_bytes, series`

`dt_s` is relative to each product's own reference (EIS step for the pseudo-raster,
mid-time for context/HMI, raster interval for XRT); `dt_mid_s` = `t_obs` − raster mid-time.

## `pixels.parquet` — the HATS-ready "object" table
| column | meaning |
|---|---|
| `epoch, t_ref` | epoch id, reference time (EIS mid-time) |
| `ix, iy, tx, ty` | grid indices and helioprojective position (arcsec, coaligned) |
| `lon_carr, lat_carr` | Carrington longitude/latitude (deg) → HATS `ra/dec` columns |
| `mu` | cos θ (centre-to-limb) |
| `aia_<λ>` | AIA context channels, DN/pix/s |
| `aia_193_pseudoraster` | AIA 193 sampled at each column's EIS step time |
| `hmi_blos` | line-of-sight field, G |
| `xrt_<filter>` | DN/s |
| `eis_win<NN>_<line>` | window-summed EIS intensity (erg cm⁻² s⁻¹ sr⁻¹, not continuum-subtracted) |
| `eis_fit_<template>_{int,vel,wid}` | Gaussian-fit intensity, Doppler velocity (km/s), width (Å) |
| `dt_<channel>` | t_obs − t_ref (s); max |Δt| for the pseudo-raster |

## `native.parquet` — for cross-matching between modalities
`epoch, modality, channel, tx, ty, lon_carr, lat_carr, mu, value, t_obs` — every native
pixel of every product that falls inside the coaligned EIS FOV (rotated to t_ref).
The notebook cross-matches EIS (4″×1″) against XRT (1″) with a KD-tree and with LSDB.
