#!/usr/bin/env python
"""Two local HATS catalogs -> matched, masked spectrum/image batches."""

# LSDB 0.11's documented public re-exports lack explicit typing exports.
# pyright: reportPrivateImportUsage=false

import sys
from pathlib import Path

import numpy as np
from lsdb.streams import CatalogStream

from mmsun.cli import die, load_json
from mmsun.xmatch import crossmatch_catalogs

if len(sys.argv) != 3:
    die("Usage: python scripts/06_xmatch.py EIS_HATS XRT_HATS")

left, right = map(Path, sys.argv[1:])
try:
    metadata = load_json(right / "mmsun.json")
    print(f"Source: {metadata['source_kind']}; positions are Carrington, not celestial RA/Dec")
    matched = crossmatch_catalogs(left, right)  # validates frame/time and requires the right margin
    rows = 0
    for batch in CatalogStream(matched, shuffle=False):
        if batch.empty:
            continue
        spectra = np.ma.array(batch.spectrum_eis.tolist(), mask=batch.mask_eis.tolist(), dtype="float32")
        images = np.ma.array(batch.image_xrt.tolist(), mask=batch.image_mask_xrt.tolist(),
                             dtype="float32").reshape(-1, *metadata["image_shape"])
        # Feed these masked arrays to an analysis here; do not accumulate batches.
        print(f"spectra {spectra.shape}, images {images.shape}")
        rows += len(batch)
    if not rows:
        die("No spatially and temporally eligible pairs")
    print(f"STREAM PASSED: {rows} pairs")
except (ValueError, FileNotFoundError, KeyError) as exc:
    die(str(exc))
