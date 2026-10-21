# MMU / HATS / LSDB: minimum solar prototype

## Bottom line and scope

**An in-memory crossmatch is not yet a HATS-backed dataset.** The remaining proof is persisted/reopened catalogs with margins, time-safe matching, and actual scientific payload access. This is a bounded packaging-and-semantics step, conditional on valid calibration/alignment—not production infrastructure. Local code/data were not inspected; “current in-memory demo” is supplied context, not independently verified.

## What the primary sources establish

| Concern | Evidence-backed fact | Minimum solar recommendation |
|---|---|---|
| Row grain | MMU parent datasets are survey-based, with object examples containing observations/arrays. HATS also allows individual-observation rows; `_healpix_29` is **not unique**. [1,3] | Use **epoch × spatial sample** rows, separately per modality for the demonstration. Common-grid pixels/spaxels suffice; do not invent persistent solar “objects.” Preserve sample IDs and provenance. |
| Payloads | Published MMU HATS catalogs contain real payloads: Galaxy10 stores image bytes, while SDSS stores nested lists of flux, wavelength, inverse variance and masks. HATS is not restricted to scalar metadata or URLs. [2,9] | Pixel-level scalar channels plus an EIS spectral list are sufficient if those are the intended training inputs. If models need images/cubes, retain accessible arrays with row-to-array indexing; integrated intensities alone do not preserve spectra or image context. |
| Spatial storage | HATS uses adaptive NESTED HEALPix partitions with Parquet data and catalog metadata; it is not simply an arbitrary Parquet file with longitude/latitude columns. [3] | Let the writer generate the layout and metadata; a tiny catalog needs no elaborate partition tuning. HEALPix partitions are storage units, not the science image grid. |
| Matching | LSDB matches coordinate pairs geometrically, with an explicit radius and neighbor count. Its KD-tree source converts numeric longitude/latitude to unit-sphere vectors. It does not perform solar alignment or time selection. [5,6] | Match only scientifically eligible samples; retain separation and original timestamps. Spatial proximity is not proof of the same physical plasma. |

## Carrington longitude/latitude: mechanically usable, not celestial RA/Dec

**Fact:** LSDB's numeric matching kernel can operate on Carrington longitude/latitude passed through its RA/Dec column selectors. That preserves spherical angular geometry, but it is **not a coordinate transformation**. The HATS v1.0 IVOA **Note** lists `hats_frame` values `equatorial` (ICRS), `galactic`, and `ecliptic`, not Carrington. The Note explicitly is not an adopted IVOA standard. [3,6]

**Recommendation:** a private Carrington-indexed prototype is a defensible shortcut only with explicit disclosure and same-frame inputs. Keep meaningful `lon_carr`/`lat_carr` names; never falsely label these ICRS. Record degrees, longitude wrapping, frame, reference time/time scale, observer, assumed solar radius/emission height, and rotation/alignment model in the dataset contract. `hats_frame=carrington` would be a project extension, not documented interoperability. Generic sky-catalog use requires genuine coordinates in a supported celestial frame, with the appropriate observer/time transformation—not renaming columns.

**Fact:** SunPy's Carrington frame is Sun-centered and rotating; its apparent longitude depends on observer/light travel time. Ordinary coordinate transformation across times does not itself follow differential rotation; SunPy supplies a separate surface-propagation model. [7,8]

Consequently, all compared coordinates must share the declared convention/reference epoch. A Carrington `_dist_arcsec` is **heliographic angular separation**, not image-plane arcseconds: do not reuse a 2″ projected-pixel tolerance numerically without conversion. A photospheric surface assumption for coronal emission is a modeling limitation, not something HATS solves.

## Smallest missing end-to-end proof

These are recommendations, not claims about existing implementation:

1. **Persist and reopen.** Build small modality catalogs with `lsdb.from_dataframe`, generate margins, then `Catalog.write_catalog(..., as_collection=True)` and reopen from their directories in a fresh session. The current documented writer is `write_catalog`. Include science values or resolvable payload references. [4]
2. **Make boundaries safe.** Attach the right-hand margin and choose its threshold at least as large as the matching radius. Set `require_right_margin=True` explicitly: the current crossmatch API defaults to `False`, despite older tutorial wording. Compare reopened results against a tiny brute-force angular reference, including a partition-edge case. [3,5,6]
3. **Make time safe.** Restrict the tiny demonstration to one selected epoch and already time-eligible observations, retaining exposure intervals and EIS step times. HATS has catalog-level temporal metadata, but default LSDB crossmatch has no time-window argument. Filtering a spatial nearest neighbor *afterward* can discard it while missing a different time-valid neighbor. Across repeated epochs, select time-eligible candidates before nearest-neighbor selection. [3,5]
4. **Prove ML consumption.** Load one matched four-modality example into numerical arrays with documented shapes, axes, units, masks and uncertainties. Inline Parquet arrays are sufficient; alternatively publish the referenced FITS/HDF5 payloads and deterministic slices. Do not duplicate a whole cube per pixel. HATS can carry the payload, but does not define tensor semantics, normalization or training splits. `CatalogStream` streams catalog chunks, not a complete training dataset adapter. [1,9,10]

## Access, acceptance and unknowns

For this prototype, local HATS collections plus payload files and a short schema/data card suffice. Publish the intact collection, margins and referenced payloads together only when sharing is needed; the MMU blog demonstrates Hugging Face access, but hosting is not required to establish backend correctness. [2,4]

**Acceptance:** reopen → time-safe, margin-complete crossmatch → recover one usable multimodal array example. No cluster, cloud service or new tensor storage system is necessary. Unknown here: actual persisted products, alignment quality, payload completeness and installed API versions. Carrington support beyond numeric reuse was not established; it must not be advertised as native solar-coordinate support.

## Primary sources

[1] MMU [README](https://github.com/MultimodalUniverse/MultimodalUniverse) and [design/schema](https://github.com/MultimodalUniverse/MultimodalUniverse/blob/main/DESIGN.md).

[2] MMU HATS contributors' [conversion and access report](https://huggingface.co/blog/hugging-science/multimodal-universe-hats), corroborated by implementation/docs below; laptop benchmark claims are not prototype requirements.

[3] [HATS v1.0 IVOA Note](https://www.ivoa.net/documents/Notes/HATS/20250822/NOTE-hats-ivoa-1.0-20250822.html), §§3.1–3.3 and properties table.

[4] LSDB [`from_dataframe`](https://docs.lsdb.io/en/stable/reference/api/lsdb.from_dataframe.html) and [`Catalog.write_catalog`](https://docs.lsdb.io/en/stable/reference/api/lsdb.catalog.Catalog.write_catalog.html).

[5] LSDB [`Catalog.crossmatch`](https://docs.lsdb.io/en/stable/reference/api/lsdb.catalog.Catalog.crossmatch.html).

[6] LSDB [KD-tree implementation](https://docs.lsdb.io/en/stable/_modules/lsdb/core/crossmatch/kdtree_match.html), especially validation and `_get_point_coordinates`.

[7] SunPy [`HeliographicCarrington`](https://docs.sunpy.org/en/stable/api/sunpy.coordinates.HeliographicCarrington.html).

[8] SunPy [`propagate_with_solar_surface`](https://docs.sunpy.org/en/stable/api/sunpy.coordinates.propagate_with_solar_surface.html).

[9] Published HATS dataset schemas: [Galaxy10](https://huggingface.co/datasets/UniverseTBD/mmu_gz10/blob/main/README.md), [SDSS](https://huggingface.co/datasets/UniverseTBD/mmu_sdss_sdss/blob/main/README.md).

[10] LSDB [`CatalogStream`](https://docs.lsdb.io/en/stable/reference/api/lsdb.streams.CatalogStream.html).
