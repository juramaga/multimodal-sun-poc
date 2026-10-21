# Multimodal Sun PoC — Phases 0-4.  `make help` lists targets.
# EPOCHS can be overridden:  make all EPOCHS="anchor"   |  EPOCHS="seq"  |  EPOCHS="all"
PY      ?= python
EPOCHS  ?= anchor
E        = --epochs $(EPOCHS)

.PHONY: help env check test notebook-test discover match download-dry download prepare coalign notebook hats all du clean-products clean-all

help:
	@echo "env            create the conda env 'mmsun' (one-off)"
	@echo "check          preflight: versions, JSOC email, archive reachability, disk"
	@echo "test           offline end-to-end test (bundled EIS + synthetic SDO/XRT; no network)"
	@echo "discover       Phase 0-1: metadata-only discovery          -> data/manifest/<epoch>/cand_*.parquet"
	@echo "match          Phase 2: choose the minimal file set        -> data/manifest/<epoch>/match.parquet"
	@echo "download-dry   show what would be downloaded and its size"
	@echo "download       Phase 3a: fetch selected files              -> data/raw/"
	@echo "prepare        Phase 3b: calibrate                         -> data/prep/"
	@echo "coalign        Phase 4: coalign, common grid, tables       -> data/aligned/<epoch>/"
	@echo "notebook       execute notebooks/01_summary.ipynb in place"
	@echo "notebook-test  execute the notebook on the offline fixture (synthetic SDO/XRT)"
	@echo "hats           export prepared EIS spectra + XRT patches as HATS collections"
	@echo "all            check -> discover -> match -> download -> prepare -> coalign -> notebook"
	@echo "du             disk estimate + measured usage"

env:
	conda env create -f environment.yml || conda env update -f environment.yml --prune

check:
	$(PY) scripts/00_check_env.py $(E)

test:
	$(PY) -m pytest -q -s tests/test_offline_pipeline.py tests/test_xmatch.py

discover:
	$(PY) scripts/01_discover.py $(E)

match:
	$(PY) scripts/02_match.py $(E)

download-dry:
	$(PY) scripts/03_download.py $(E) --dry-run

download:
	$(PY) scripts/03_download.py $(E)

prepare:
	$(PY) scripts/04_prepare.py $(E)

coalign:
	$(PY) scripts/05_coalign.py $(E)

notebook:
	MMSUN_EPOCH=$(EPOCHS) jupyter nbconvert --to notebook --execute --inplace notebooks/01_summary.ipynb

notebook-test:    ## run the notebook on the offline fixture (after `make test`); output in data_test/
	MMSUN_CONFIG=$(CURDIR)/data_test/pipeline_test.yaml jupyter nbconvert --to notebook --execute \
	  notebooks/01_summary.ipynb --output-dir data_test --output 01_summary_fixture.ipynb

hats:
	$(PY) scripts/build_hats.py $(E)

all: check discover match download prepare coalign notebook du

du:
	$(PY) scripts/estimate_disk.py $(E) --measure

clean-products:   ## keep raw downloads, drop derived products
	rm -rf data/prep data/aligned

clean-all:
	rm -rf data data_test
