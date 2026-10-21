# EIS/XRT prototype: what happens

The [reader example](../scripts/06_xmatch.py) is deliberately small. It opens two
local HATS collections, checks their shared frame/time contract, matches with
LSDB, and yields masked NumPy arrays through upstream `CatalogStream`.
See the [README](../README.md) for installation and the four demo commands.

## Conversion and schema

`scripts/build_hats.py` calls `build_catalogs` in `src/mmsun/xmatch.py`:

1. Read one raw EIS HDF5 window using EISPAC and one prepared XRT DN/s FITS map.
2. Check pointing flags/bounds, units and exposure time. Apply the saved shifts
   (`true = header + d`) and transport coordinates to the raster midpoint.
3. Make one row per EIS spatial pixel and one per overlapping XRT context patch.
4. Give arrays explicit Arrow-list types, then let `lsdb.from_dataframe` and
   `write_catalog(as_collection=True)` create spatial Parquet partitions and margins.

| Catalog | Inline arrays | Default shape / units |
| --- | --- | --- |
| EIS | `wavelength`, `spectrum`, `error`, `mask` | 24 spectral bins in the selected Fe XII window |
| XRT | `image`, `image_mask` | 16×16 native pixels, flattened in C order; DN/s |

EIS wavelengths are float64 Angstroms. Spectral values, standard-deviation errors
and images are float32; masks are boolean, **True = invalid**. EIS uses EISPAC's
bundled pre-flight radiometric calibration; this is not a new calibration-error
model. No XRT uncertainty model, normalization or imputation is supplied.
The reader uses masked arrays so invalid entries are not treated as measurements.

Both tables include sample ID, epoch, native pixel indices, Carrington longitude/
latitude and `dt_s` (acquisition time minus `t_ref`, in UTC seconds). The HEALPix
index is a storage index, not a unique sample ID. `mmsun.json` records shapes,
units, reference time, coordinate convention, source names, pointing corrections
and software versions; XRT also records its original WCS.

```text
<data_dir>/hats/<epoch>/
  eis/    HATS collection + mmsun.json
  xrt/    HATS collection + persisted right margin + mmsun.json
```

Open/copy the entire collection roots, not just their nested main catalogs.
Payloads live inside Parquet: consuming them needs no original FITS/HDF5 files.

## Checks that are not optional scaffolding

- **Coordinates:** Carrington degrees on a photospheric sphere, Earth observer.
  This is a solar-specific convention, not celestial RA/Dec or native solar-frame
  support in HATS. Inputs must share frame, observer, surface, epoch and reference
  time. Off-limb/non-finite positions are rejected.
- **Time:** one XRT exposure. EIS keeps each slit column's acquisition time.
  Eligible rows are filtered **before** spatial nearest-neighbor selection,
  including the right margin. This is not arbitrary multi-exposure time matching.
- **Margins:** the right collection includes a margin slightly wider than the
  matching radius. LSDB rejects a missing or insufficient margin.
- **Output safety:** only new destinations are allowed; existing or symlinked
  paths are refused. No raw/prepared products are rewritten. To retry an export,
  explicitly move aside the old or incomplete epoch's HATS directory first.

## Deliberate limits

- The offline fixture has real EIS and **synthetic XRT** (plus synthetic AIA/HMI
  used by the original alignment test). `--fixture` records that label. Other
  input products are labelled `user-provided`, not automatically science-validated.
- The exporter consumes the original pipeline's EIS pair, `products.parquet` and
  `shifts.json`. It does not download or repair inadequately calibrated inputs.
- The builder loads one small epoch in memory. `CatalogStream` yields partitions,
  not fixed-size minibatches or a strict RAM cap. Larger ingestion belongs in
  upstream `hats-import`, not a new framework here.
- Geometric transport uses the existing FOV/column-centre displacement model.
  It does not correct plasma evolution or recover coronal heights; large or
  near-limb rasters need full reprojection.
- Matching compares **patch centres**, not footprints. Default radius is 0.5
  heliographic degrees (~6 Mm), not 0.5 sky degrees. Patches overlap at 8-pixel
  strides and are not independent training examples.

Settings are under `xmatch` in `config/pipeline.yaml`; time tolerance uses
`tolerance.xrt_ar_s`. Tests cover persisted arrays/NaNs/masks, boundary neighbors,
wrong-time/epoch exclusion, coordinate transport, metadata/margin rejection,
CLI consumption and output safety. Live-XRT science validation and remote
publication remain out of scope.

## Upstream, not copied engines

LSDB/HATS are pinned to **0.11.0**. We use their public
[dataframe loader](https://github.com/astronomy-commons/lsdb/blob/v0.11.0/src/lsdb/loaders/dataframe/from_dataframe.py),
[collection writer](https://github.com/astronomy-commons/lsdb/blob/v0.11.0/src/lsdb/io/to_collection.py),
[catalog matching](https://github.com/astronomy-commons/lsdb/blob/v0.11.0/src/lsdb/catalog/catalog.py)
and [CatalogStream](https://github.com/astronomy-commons/lsdb/blob/v0.11.0/src/lsdb/streams/catalog_streams.py).
Both libraries are BSD-3-Clause. No AGPL benchmark implementation or unlicensed
tutorial code was copied. The original solar preparation/rotation helpers remain
unchanged.
