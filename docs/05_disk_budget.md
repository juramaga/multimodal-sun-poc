# 05 · Disk budget

Estimates from `python scripts/estimate_disk.py` (anchor geometry, default config), and
the measured size of the offline test run (synthetic SDO/XRT saved as float32, i.e.
somewhat larger than real int16 JSOC cutouts).

## Per epoch

| Area | Typical | Worst case | Notes |
|---|---|---|---|
| raw — EIS | 3.6 MB | 3.6 MB | HDF5 data + head |
| raw — AIA pseudo-raster | 5.2 MB | 5.2 MB | ~26 × 193 Å cutouts, 180″ × 200″ |
| raw — AIA context | 6.5 MB | 6.5 MB | 7 × cutouts, 400″ × 420″ |
| raw — HMI | 2.6 MB | 2.6 MB | 1 magnetogram cutout |
| raw — XRT | 2.1 MB | 33.6 MB | 2 filters; partial frame vs full-frame 2048² |
| XRT header cache | ~3 MB | ~3 MB | range requests, reused across epochs |
| prep (float32 copies + EIS fits) | 29 MB | 61 MB | |
| aligned (stack, tables, diagnostics) | 7 MB | 7 MB | `native.parquet` dominates (~5.5 MB) |
| **Total** | **≈ 60 MB** | **≈ 120 MB** | |

Measured (offline test, real EIS + synthetic others): raw 31.5 MB, prep 28.1 MB,
aligned 6.5 MB → **66 MB**.

## One-off

* conda environment `mmsun`: ~1.5–3 GB (the pip venv used to build this was 1.1 GB).
* sunpy/aiapy calibration tables cache: < 20 MB.

## Growth

| Scope | Data |
|---|---|
| Anchor only | 60–120 MB |
| Full 06:25–07:04 sequence (9 rasters, `--epochs seq` + anchor) | 0.5–1.1 GB (XRT frames are often shared between neighbouring rasters, so usually toward the low end) |
| Adding a large 150 MB raster (e.g. 18:12 UT) | +~0.4 GB (bigger FOV → bigger cutouts) |

Recommendation: ≥ 5 GB free for the anchor (environment included), ~6 GB for the sequence.
`make clean-products` frees prep/aligned (rebuildable) without touching raw downloads.
