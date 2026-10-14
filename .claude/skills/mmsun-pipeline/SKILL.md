---
name: mmsun-pipeline
description: Run, debug or extend the Multimodal Sun proof-of-concept pipeline (Phases 0-4 - EIS-anchored discovery, matching, download, calibration, coalignment of AIA/HMI/XRT/EIS). Use for any request to run the pipeline, add an epoch, change tolerances, add a channel, or interpret its outputs.
---

# Multimodal Sun pipeline

The project plan and hard rules live in `CLAUDE.md`; read it first. This skill is the
operating procedure.

## Run order (each step is idempotent and reads the previous step's files)

| Step | Command | Reads | Writes |
|---|---|---|---|
| preflight | `python scripts/00_check_env.py` | config | – (exit 2 = fix something) |
| 0-1 discover | `python scripts/01_discover.py` | EIS head.h5, JSOC keywords, XRT headers | `data/manifest/<epoch>/epoch.json, eis_steps.parquet, cand_*.parquet` |
| 2 match | `python scripts/02_match.py` | `cand_*` | `match.parquet, eis_step_map.parquet, match_report.json` |
| 3a download | `python scripts/03_download.py [--dry-run]` | `match.parquet` | `data/raw/<mod>/<epoch>/` |
| 3b prepare | `python scripts/04_prepare.py` | raw | `data/prep/<mod>/<epoch>/, products.parquet` |
| 4 coalign | `python scripts/05_coalign.py` | prep | `data/aligned/<epoch>/{shifts.json, stack.fits, pixels.parquet, native.parquet}` |
| summary | `make notebook` | everything | `notebooks/01_summary.ipynb` (executed) |

All scripts take `--epochs <id|tag|anchor|all>`; default comes from
`config/pipeline.yaml:active_epochs`.

## Decision points (stop and report to the user, do not improvise)

1. `02_match.py` exits with code 3 → a modality is MISSING. Show `match_report.json`.
   The standard remedy is the next raster of the sequence (`--epochs seq`), not
   loosening tolerances. Loosen a tolerance only if the user agrees, and only in
   `config/pipeline.yaml` (never hard-code).
2. `03_download.py --dry-run` estimate > 300 MB for one epoch → something is wrong
   (probably full-disk XRT or a widened cutout). Report before downloading.
3. `shifts.json` has `flagged: true` (correlation < `coalign.min_correlation`) or a shift
   at the search boundary (|dx| or |dy| ≈ `max_shift_arcsec`) → coalignment failed;
   report r_before/r_after and do not proceed to Phase 5 work with that epoch.
4. A JSOC export returns zero files or files much larger than the cutout box → the
   `im_patch` step failed silently (a known JSOC behaviour). Re-check the box and `t_ref`.

## Extending

* **More epochs**: add rows to `config/epochs.yaml` (see the `epoch-scout` skill), then
  `--epochs <new id>`. No code change.
* **More channels**: `jsoc.aia_context_channels`, `jsoc.hmi_series` (segment names are
  mapped in `scripts/02_match.py:HMI_SEGMENT` and `scripts/03_download.py:SEGMENT`),
  `xrt.max_filters`, `eis.fit_templates` (names of files in `eispac/data/templates`).
* **New instrument** (e.g. IRIS, SOT): follow the pattern `src/mmsun/<inst>.py` with
  discovery → selection in `match.py` → retrieval → `prepare()` → a branch in
  `05_coalign.py` (reference choice + `map_vs_reference`). Keep the
  "true = header + d" convention.

## Verification after any code change

`make test` (offline: bundled EIS + synthetic SDO/XRT with injected pointing errors;
must print `OFFLINE TEST PASSED` and recover the injected shifts to < 1.5").
