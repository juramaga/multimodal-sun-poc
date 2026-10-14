---
name: epoch-scout
description: Find and vet new anchor epochs (EIS rasters) that maximise the chance of co-temporal, co-spatial AIA/HMI/XRT coverage, and add them to config/epochs.yaml. Use when the user wants more data, a different target type (flare, quiet Sun, limb), or when an epoch fails matching.
---

# Choosing anchor epochs

EIS is the scarcest modality, so epochs are EIS rasters. SDO is effectively always
available (check QUALITY), so the questions are: is XRT pointed at the same place at
the same time, and is the raster cheap and geometrically simple?

## Criteria, in priority order

1. **XRT co-pointing.** Hinode instruments share spacecraft pointing, so an EIS raster
   taken during a HOP or AR campaign in which XRT takes partial-frame images of the same
   target is almost always covered. Prefer observing titles mentioning a HOP, an NOAA AR,
   or "flare". Fallback: rasters bracketing the ~06/~18 UT XRT synoptic (full disk).
2. **Near disk centre** (|x|, |y| < ~600″): little foreshortening, all pixels on disk,
   clean Carrington coordinates, AIA/HMI comparisons straightforward.
3. **Short rasters** (≤ 10 min): small `dt` across modalities, small files
   (sparse rasters like the anchor: 25 steps, 3.4 MB).
4. **Windows overlapping AIA/XRT temperature range**: Fe XII 192/195 for coalignment;
   hot lines (Ca XV, Fe XVII, Fe XXIII) for XRT science.
5. Avoid XRT bakeout weeks and Hinode eclipse season (≈ May–Aug) gaps.

## Procedure

1. List rasters for a day: `https://eis.nrl.navy.mil/level1/hdf5/YYYY/MM/DD/` (file
   sizes hint at raster size), or `Fido.search(a.Time(...), a.Instrument("EIS"))`.
2. For candidates, fetch only `*.head.h5` and inspect with
   `mmsun.eis.head_info()` (duration, pointing, FOV, windows, title, NOAA).
3. Optionally check the HCR JSON for co-observations (see `solar-archive-access`).
4. Add rows to `config/epochs.yaml` with tags, then
   `python scripts/01_discover.py --epochs <id>` and `python scripts/02_match.py --epochs <id>`.
   Keep an epoch only if `match_report.json` has no MISSING modality.

## Current registry

* Anchor `eis_20210306_064444` — HOP 396 (Chandrayaan-2 XSM collaboration), AR 12807,
  flare-watch study `cam_flare_hot1_study`, near disk centre, bundled in eispac.
* Tag `seq`: the 8 sibling rasters 06:25–07:04 UT of the same sequence.
* Commented out: `eis_20210306_181244` (large 150 MB raster) for later.
