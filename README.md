# Multimodal Sun — minimal EIS/XRT crossmatch

A small version of the [MMU streaming workflow](https://huggingface.co/blog/hugging-science/multimodal-universe-hats):
**two local HATS catalogs → LSDB crossmatch → masked spectrum/image batches**.

**Start with [`scripts/06_xmatch.py`](scripts/06_xmatch.py).** It is the complete
reader example. Catalog preparation is separate in `scripts/build_hats.py`;
LSDB/HATS does the partitioning, margins, matching and streaming.

## Install (Python 3.12)

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.txt
```

A normal venv + pip works too. No Claude installation is required.

## Run the offline example

```bash
fixture=$(mktemp -d)
MMSUN_TEST_DATA="$fixture" make test PY=.venv/bin/python
.venv/bin/python scripts/build_hats.py --config "$fixture/pipeline_test.yaml" --fixture
.venv/bin/python scripts/06_xmatch.py "$fixture"/hats/*/eis "$fixture"/hats/*/xrt
```

Run these in order; stop if a command fails. The fresh directory keeps existing
raw data untouched. Tests print `OFFLINE TEST PASSED`; the reader prints batch
shapes and `STREAM PASSED`. **EIS is real bundled data; XRT is synthetic.**
This proves the software, not live-data science validity. No archive access or
JSOC account is needed after installation.

## Use existing data

```bash
make hats PY=.venv/bin/python EPOCHS=anchor  # requires a prepared, coaligned epoch
.venv/bin/python scripts/06_xmatch.py /path/to/eis /path/to/xrt
```

Exports refuse existing/symlinked destinations. Coordinates are **Carrington
solar-sphere angles, not celestial RA/Dec**. Scope: one EIS window and one XRT
exposure; per-column EIS times and invalid-data masks are retained.

- [Schema, safeguards, limitations and upstream references](docs/08_eis_xrt_prototype.md)
- [Original four-instrument acquisition/preparation workflow](docs/09_original_pipeline.md)
