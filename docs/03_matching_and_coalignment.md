# 03 · Matching (Phase 2) and coalignment (Phase 4)

## "Simultaneous" per instrument

An EIS raster is not a snapshot: each of its 25 columns has its own time (~11.6 s
apart). Matching therefore uses three different notions of time:

| Product | Reference time | Tolerance | Why |
|---|---|---|---|
| AIA 193 pseudo-raster | each EIS slit step | 12 s (one cadence) | needed to coalign EIS against the same moment in every column |
| AIA context (7 channels), HMI | raster mid-time | 12 s / 45 s | one simultaneous multi-channel snapshot; one T_REC slot for all AIA channels |
| XRT | the raster *interval* (Δt = 0 inside it); ties → nearest the mid-time | 300 s, synoptic 1800 s | XRT cadence depends on the observing program |

Spatial match: AIA/HMI are cut out around the EIS FOV, so they always contain it. XRT
frames must contain ≥ 90% of the EIS FOV (`tolerance.min_fov_coverage`).
Every selected product keeps its Δt (`dt_s`, to its own reference time) and its offset
from the raster mid-time (`dt_mid_s`), so stricter subsets can be made later without
re-downloading. For XRT, `dt_s` is 0 for every frame inside the raster; use `dt_mid_s`
to tighten. Among in-raster frames of one filter, the one nearest the mid-time wins
(then the longest exposure), so long rasters do not pick a frame from their far end.
`04_prepare.py` calibrates only the XRT frames listed in `match.parquet`; `raw/` may
hold frames from earlier matches.

## Coalignment

AIA is the reference. Search: maximise the Pearson correlation of log intensities over
shifts (4″ → 1″ → 0.25″ grids, ±40″).

* **EIS vs AIA 193 pseudo-raster**: for each EIS column, sample the AIA frame nearest
  that column's time at the EIS pixel positions (AIA blurred to ~3″). Reference
  window: Fe XII 192.39 Å. One shift is applied to all EIS windows (eispac already
  handles the relative CCD offsets).
* **XRT vs AIA**: on a 2″ grid over the EIS FOV + 30″, try AIA 335, 94, 211; keep the
  best. Solar rotation between the XRT and AIA times is applied before comparing.
* **HMI**: assumed coaligned with AIA (both JSOC level-1 with aiapy pointing).

Convention: `true position = header position + d`.

Then every channel is resampled to one helioprojective grid (default 2″/px) covering
the coaligned EIS FOV, at the EIS mid-time, with Carrington lon/lat and μ per pixel.

## Acceptance

`shifts.json` per instrument: `dx, dy, r_before, r_after, reference, flagged`.
Good: r_after ≳ 0.8, clearly above r_before, |d| well inside ±40″. Typical sizes:
EIS 3–15″, XRT 1–5″. The offline test recovers injected shifts to < 1″.
