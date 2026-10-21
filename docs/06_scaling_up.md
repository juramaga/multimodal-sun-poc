# 06 · Scaling up (and the hand-off to Phase 5, HATS)

## Growing the sample without code changes

1. **Same sequence**: `make all EPOCHS=seq` (8 more rasters, same study/target).
   Expect EIS shifts to drift smoothly between rasters; XRT frames may be shared.
2. **More epochs**: add rows to `config/epochs.yaml` (procedure in the `epoch-scout`
   skill). Each row is just an EIS file stem plus tags.
3. **More channels**: AIA UV (1600/1700 via `aia.lev1_uv_24s`), `hmi.Ic_45s`,
   `hmi.B_720s`, more XRT filters (`xrt.max_filters`), more EIS fit templates.
4. **Tighter or looser simultaneity**: change `tolerance.*`; Δt is already stored, so
   a stricter subset can be selected from existing products.

## What changes at scale (when the sample is thousands of epochs)

* Discovery becomes a catalogue job: harvest all EIS `head.h5` files and XRT headers
  once into Parquet catalogues (they are small), then match in bulk.
* JSOC exports should be batched (one export per T_REC range and channel set per day).
* `native.parquet` should be partitioned by epoch/day; the common grid stays per epoch.

## General Phase 5 sketch (not implemented)

A narrower **EIS spectra + XRT patches** HATS writer and streaming crossmatch is
now implemented; see [the two-modality prototype](08_eis_xrt_prototype.md).
It does not replace the general pixel/native-table catalog design below.

Both tables already carry `lon_carr`/`lat_carr`, the natural "sky" for the Sun:

```python
import lsdb, pandas as pd
px = pd.concat(pd.read_parquet(p) for p in sorted(Path("data/aligned").glob("*/pixels.parquet")))
cat = lsdb.from_dataframe(px, ra_column="lon_carr", dec_column="lat_carr", catalog_name="mms_pixels")
cat.to_hats("data/hats/mms_pixels")
```

Open design questions for Phase 5 (documented in the earlier planning notes):
* Row granularity — grid pixel (as here) vs patch (e.g. 32×32 cutout with array columns).
* Time — HATS has no time axis; partition by epoch/day in directories and keep `t_ref`.
* Off-limb data — no Carrington coordinate; needs a separate table (position angle, height).
* Carrington vs Stonyhurst — Carrington lets the same feature be cross-matched across
  epochs (rotation), Stonyhurst keeps Earth-facing geometry fixed.
