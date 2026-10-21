# Original four-instrument pipeline (Phases 0–4)

The smallest sample that contains all four modalities of the planned Multimodal Sun
dataset — **SDO/AIA, SDO/HMI, Hinode/XRT, Hinode/EIS** — matched in time and space,
calibrated, coaligned, on a common grid, and in tables that support a heliographic
cross-match query. Anchor: EIS raster `eis_20210306_064444` (AR 12807, near disk
centre, 06:44:44–06:49:34 UT). Disk: **~60 MB of data per epoch** (≤ 120 MB worst case)
plus a ~1.5–3 GB environment.

* `CLAUDE.md` — the plan, phase contracts and working rules (read by Claude Code).
* `docs/` — anchor rationale, access layers, matching/coalignment, schemas, disk budget, scaling.
* `notebooks/01_summary.ipynb` — summary and visualisation, including the cross-match demo.

## Prerequisites

* macOS or Linux, `conda` (or mamba), `make`, internet access.
* An email registered with JSOC for SDO exports (one-off, free):
  http://jsoc.stanford.edu/ajax/register_email.html
* Claude Code is optional; the terminal-only commands below do not require it.

## Run it with Claude Code

```bash
conda env create -f environment.yml
conda activate mmsun
export JSOC_EMAIL="you@institution.edu"  # use your registered email
claude
```

Then a good sequence of prompts is:

1. `Read CLAUDE.md, run make test, and tell me if the offline pipeline passes.`
   (no network; real EIS + synthetic SDO/XRT; checks injected pointing errors)
2. `Run the preflight check and Phases 0-2 for the anchor epoch, then show me match_report.json.`
   (`make check discover match` — metadata only, a few MB)
3. `Show the download dry run, and if it is under 150 MB, download, prepare and coalign.`
   (`make download-dry download prepare coalign`)
4. `Execute the summary notebook and summarise the coalignment shifts and the cross-match result.`
   (`make notebook`, then open it with `jupyter lab notebooks/01_summary.ipynb`)
5. Optional growth: `Run the same pipeline for the rest of the sequence (EPOCHS=seq).`

Useful slash commands: `/mmsun-pipeline`, `/solar-archive-access`,
`/solar-coalignment`, `/epoch-scout`, `/python-library-docs-first` (project skills).

## Run it without Claude Code

```bash
conda activate mmsun
export JSOC_EMAIL=...  # use your registered email
make test             # offline self-test
make all              # check -> discover -> match -> download -> prepare -> coalign -> notebook -> du
```

If `make match` stops with "MISSING" for a modality, follow the pipeline's stop
conditions before trying sibling rasters (`make discover match EPOCHS=seq`).
Do not weaken matching or calibration requirements merely to get an export.

## Outputs

```text
data/manifest/<epoch>/   discovery + match tables, reports, provenance
data/raw/                exactly the matched files (never edited)
data/prep/               calibrated native-resolution maps
data/aligned/<epoch>/    shifts.json, stack.fits, pixels.parquet, native.parquet
```

The newer [two-modality HATS prototype](08_eis_xrt_prototype.md) consumes prepared
products but does not change this acquisition/calibration workflow.
