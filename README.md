# Multimodal Sun — proof of concept (Phases 0–4)

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
* Claude Code (needs a Pro, Max, Team, Enterprise or Console account):
  `curl -fsSL https://claude.ai/install.sh | bash` then `claude --version`.

## Run it in the terminal with Claude Code

```bash
# 1. unpack and enter the project
unzip multimodal-sun-poc.zip && cd multimodal-sun-poc

# 2. create the environment (Python 3.12; ~5-10 min the first time)
conda env create -f environment.yml
conda activate mmsun

# 3. JSOC email for the SDO cutout exports (not stored anywhere in the repo)
export JSOC_EMAIL="you@institution.edu"

# 4. optional: put it under git so you can see what Claude changes
git init && git add -A && git commit -m "Multimodal Sun PoC scaffold"

# 5. start Claude Code in the project root (it reads CLAUDE.md and .claude/skills automatically)
claude
```

Then, inside Claude Code, a good sequence of prompts is:

1. `Read CLAUDE.md, run make test, and tell me if the offline pipeline passes.`
   (no network; real EIS + synthetic SDO/XRT; checks that injected pointing errors are recovered)
2. `Run the preflight check and Phases 0-2 for the anchor epoch, then show me match_report.json.`
   (`make check discover match` — metadata only, a few MB)
3. `Show the download dry run, and if it is under 150 MB, download, prepare and coalign.`
   (`make download-dry download prepare coalign`)
4. `Execute the summary notebook and summarise the coalignment shifts and the cross-match result.`
   (`make notebook`, then open it with `jupyter lab notebooks/01_summary.ipynb`)
5. Optional growth: `Run the same pipeline for the rest of the sequence (EPOCHS=seq).`

Useful slash commands in the session: `/mmsun-pipeline`, `/solar-archive-access`,
`/solar-coalignment`, `/epoch-scout`, `/python-library-docs-first` (project skills).

## Run it without Claude Code

```bash
conda activate mmsun && export JSOC_EMAIL=...
make test          # offline self-test
make all           # check -> discover -> match -> download -> prepare -> coalign -> notebook -> du
```

If `make match` stops with "MISSING" for a modality, try the sibling rasters:
`make discover match EPOCHS=seq`.

## What gets produced

```
data/manifest/<epoch>/   discovery + match tables, reports, provenance
data/raw/                exactly the matched files (EIS HDF5, AIA/HMI cutouts, XRT L1)
data/prep/               calibrated native-resolution maps
data/aligned/<epoch>/    shifts.json, stack.fits, pixels.parquet (HATS-ready), native.parquet
```
