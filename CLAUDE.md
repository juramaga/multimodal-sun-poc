# CLAUDE.md — Multimodal Sun proof of concept (Phases 0–4)

You are working in a research prototype for the **Multimodal Sun (MMS)** dataset: an
ML-ready, multimodal solar dataset modelled on the Multimodal Universe. This repository
builds the **smallest possible sample that still contains every modality — SDO/AIA,
SDO/HMI, Hinode/XRT and Hinode/EIS — co-temporal and co-spatial, coaligned, and in a
table that supports a cross-match query**. It must stay easy to grow.

Scope: Phases 0–4 only. Phase 5 (writing a HATS catalogue) is *not* implemented; the
Phase 4 tables are designed so that it is a thin step later (see `docs/06_scaling_up.md`).

---

## 1. The plan at a glance

| Phase | What | Script | Output (under `data/`) | Done when |
|---|---|---|---|---|
| 0 | Pick the anchor epoch (an EIS raster) | `config/epochs.yaml` | — | anchor listed and `active_epochs` set |
| 1 | Metadata-only discovery for all four modalities | `scripts/01_discover.py` | `manifest/<epoch>/epoch.json, eis_steps.parquet, cand_{aia,hmi,xrt}.parquet` | every modality has ≥1 candidate |
| 2 | Match in time and space; choose the minimal file set | `scripts/02_match.py` | `manifest/<epoch>/match.parquet, eis_step_map.parquet, match_report.json` | exit code 0 (no MISSING modality) |
| 3a | Retrieve exactly the matched files | `scripts/03_download.py` | `raw/{eis,aia,hmi,xrt}/<epoch>/` | dry-run size ≈ estimate; files present |
| 3b | Calibrate to science-ready native maps | `scripts/04_prepare.py` | `prep/{eis,aia,hmi,xrt}/<epoch>/, manifest/<epoch>/products.parquet` | every product has a record and no FAILED step that matters |
| 4 | Coalign, common grid, tables | `scripts/05_coalign.py` | `aligned/<epoch>/{shifts.json, stack.fits, pixels.parquet, native.parquet}` | no `flagged` shift; notebook runs |
| — | Summary + cross-match demo | `notebooks/01_summary.ipynb` (`make notebook`) | executed notebook | all cells run |

One command for everything: `make all` (stops at the first failing phase).
Offline self-test (no network): `make test`.

## 2. Anchor epoch (Phase 0)

**`eis_20210306_064444`** — Hinode/EIS raster 2021-03-06 06:44:44–06:49:34 UT.

* Study `cam_flare_hot1_study` (id 598): 25 steps × 4″ with the 2″ slit, 100″×120″ FOV,
  10 s exposures, 9 windows from Fe XII 186/192 Å to Fe XXIII 263.8 Å.
* Target NOAA **AR 12807**, pointing ≈ (3″, −197″) → near disk centre.
* Observing title "HOP396 Collaboration w/ Chandrayaan-2 XSM AR12807": a Hinode flare
  watch. Hinode's instruments share pointing, so XRT is expected to be imaging the same
  AR; the XRT ~06 UT synoptic (full disk) is the fallback.
* One of **9 identical rasters every ~4 m 53 s from 06:25 to 07:04 UT** (tag `seq` in
  `config/epochs.yaml`) — the first growth step needs no new decisions.
* **Bundled inside `eispac`** (`eispac/data/test/`), so EIS never depends on the network.

Why this maximises matching: short raster (small Δt for every modality), disk-centre AR
(no limb geometry, clean heliographic coordinates), hot EIS lines that overlap XRT's
response, and a coordinated campaign. Details: `docs/01_anchor_epoch.md`.

## 3. What is (and is not) downloaded

Only science-ready products, only what matching selected:

| Modality | Downloaded | Size (typ.) | Not downloaded |
|---|---|---|---|
| EIS | Level-1 HDF5 pair (`.head.h5` + `.data.h5`) for the anchor raster | 3.6 MB | Level-0 FITS, other rasters |
| AIA | `aia.lev1_euv_12s` **cutouts** (JSOC `im_patch`): 193 Å, one frame per EIS step (~26 frames, EIS FOV + 40″); 7 EUV channels at the raster mid-time (EIS FOV + 150″) | ~12 MB | full-disk frames, UV/visible, spikes files |
| HMI | `hmi.M_45s` magnetogram cutout nearest the mid-time | ~3 MB | vector/SHARP, Dopplergrams |
| XRT | Level-1 FITS for ≤ 2 filters nearest the raster that contain ≥ 90% of the EIS FOV (headers inspected first via HTTP range requests) | 2–34 MB | Level-0, grade maps, frames not covering the FOV |

Total for the anchor: **~60 MB typical, ~120 MB worst case** (full-frame XRT), plus a
~1.5–3 GB conda environment. `python scripts/estimate_disk.py` prints the breakdown;
`make du` adds what is actually on disk. See `docs/05_disk_budget.md`.

## 4. Matching rules (Phase 2)

* AIA 193 pseudo-raster: nearest frame to **each EIS slit step**, |Δt| ≤ 12 s.
* AIA context (94, 131, 171, 193, 211, 335, 304): one T_REC slot nearest the raster
  mid-time, |Δt| ≤ 12 s, `QUALITY == 0`.
* HMI: nearest `hmi.M_45s` to the mid-time, |Δt| ≤ 45 s, `QUALITY == 0`.
* XRT: `EC_IMTY_ == normal`, FOV covers ≥ 90% of the EIS FOV, Δt to the raster
  *interval* ≤ 300 s (ties inside the raster → nearest the mid-time); best 2 filters in
  the order Al_poly, Be_thin, Al_mesh, Ti_poly…; else synoptic full-disk within 1800 s.
* Δt is stored for every product so tolerances can be tightened later without
  re-downloading. All tolerances live in `config/pipeline.yaml`.

## 5. Coalignment (Phase 4)

AIA is the reference. EIS is matched to an AIA 193 *pseudo-raster*; XRT to the
best-correlating of AIA 335/94/211; solar rotation is applied from each map's time to
the EIS mid-time. Convention: **true = header + d**. All channels are then resampled to
one 2″ helioprojective grid covering the coaligned EIS FOV, with Carrington lon/lat per
pixel. Details and acceptance criteria: `.claude/skills/solar-coalignment/SKILL.md`.

## 6. Rules for working in this repo

1. **Run phases in order** and read each phase's printed verdict before continuing.
   Never run `03_download.py` without a successful `02_match.py` for that epoch.
2. **Stop and ask the user** when: `02_match.py` exits with code 3 (missing modality);
   the download dry-run exceeds ~300 MB per epoch; a shift is `flagged` or sits at the
   search boundary; JSOC returns no files or full-disk files. Report the evidence
   (`match_report.json`, `shifts.json`) rather than working around it.
3. **Configuration over code.** Tolerances, channels, margins, filters and epochs are in
   `config/`. Don't hard-code them in scripts.
4. **Never modify or delete `data/raw/`** — it is the provenance record. Derived
   products can be rebuilt (`make clean-products`).
5. **Credentials**: the JSOC email comes only from `$JSOC_EMAIL`. Never write it into
   files, configs, notebooks or commit messages.
6. **Be polite to archives**: keep the XRT header cache, ≤ 8 parallel requests, no
   directory crawling beyond the hour folders needed.
7. **Version-sensitive APIs** (sunpy, eispac, drms, aiapy, lsdb change between minor
   versions): follow `.claude/skills/python-library-docs-first` — check the installed
   version's signature before changing a call, and say which version you checked.
8. **After any change under `src/` or `scripts/`**: run `make test` (offline, ~1 min);
   it must print `OFFLINE TEST PASSED`. If you change the notebook, edit
   `tools/build_notebook.py` and regenerate (`python tools/build_notebook.py`) — do not
   hand-edit the `.ipynb`.
9. Use the project skills: `mmsun-pipeline` (operating procedure),
   `solar-archive-access` (access layers and pitfalls), `solar-coalignment`,
   `epoch-scout` (adding epochs).

## 7. Layout

```
config/            pipeline.yaml (all knobs), epochs.yaml (epoch registry)
src/mmsun/         library: eis, xrt, jsoc, match, prep, coalign, grid, config, cli, netutil
scripts/           00_check_env, 01_discover, 02_match, 03_download, 04_prepare, 05_coalign, estimate_disk
notebooks/         01_summary.ipynb (generated by tools/build_notebook.py)
tests/             offline end-to-end test with synthetic SDO/XRT around the real EIS raster
docs/              anchor epoch, access layers, matching/coalignment, outputs/schema, disk budget, scaling up
.claude/skills/    project skills (+ one vendored from Skill Commons, see PROVENANCE.md)
data/              everything downloaded or produced (git-ignored)
```

## 8. Growth path (not now, but don't block it)

1. `--epochs seq`: the other 8 rasters of the sequence (~0.5 GB total).
2. More epochs via the `epoch-scout` skill; more channels via config.
3. Phase 5: `pixels.parquet`/`native.parquet` → `lsdb.from_dataframe(...,
   ra_column="lon_carr", dec_column="lat_carr")` → `to_hats()`; then `hats-import` at scale.

## 9. Known risks

* XRT coverage for the anchor is *expected* (co-pointed HOP) but only verified when
  `01_discover.py` runs on a networked machine. If missing → `--epochs seq`.
* JSOC `im_patch` can fail silently; check cutout sizes.
* EIS long-wavelength windows are offset in y by ~17″ (CCD offset): blank rows at the
  grid edges for those layers are expected.
* EIS pointing conventions differ between `index/` and `pointing/` groups in the HDF5
  header (~20″); the code uses eispac's convention and coalignment absorbs the residual.
