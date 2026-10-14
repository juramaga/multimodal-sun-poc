# 02 · Access layers queried

| Modality | Discovery (metadata only) | Retrieval | Product level | Code |
|---|---|---|---|---|
| EIS | NRL `level1/hdf5/YYYY/MM/DD/` listing; `*.head.h5` for times, pointing, windows | `*.data.h5` from NRL (MSSL mirror fallback) or the copy bundled in eispac | Level-1 HDF5 (calibratable count rates + radcal + pointing) | `mmsun/eis.py` |
| AIA | JSOC keywords via `drms.Client().query` (T_REC, T_OBS, WAVELNTH, QUALITY, EXPTIME) | JSOC export with `im_patch` cutouts (`drms.Client(email).export`) | `aia.lev1_euv_12s` level 1 → aiapy pointing + degradation + exposure | `mmsun/jsoc.py`, `mmsun/prep.py` |
| HMI | JSOC keywords (`hmi.M_45s`) | JSOC export with `im_patch` | `hmi.M_45s` magnetogram (science-ready) | same |
| XRT | SDAC hour-folder listing + FITS header via HTTP Range requests (filter, FOV, exposure, image type) | single Level-1 FITS files from SDAC | Level 1 (xrt_prep'd) → DN/s | `mmsun/xrt.py` |
| Coverage (optional) | LMSAL HEK Coverage Registry (HCR) JSON; DARTS Hinode search | — | — | `epoch-scout` skill |

Why not one federated search (VSO/Fido) for everything? Fido works, but for Hinode it
can only search by **time**: neither EIS pointing/study nor XRT filter/FOV are
searchable, so the decisive selection would happen only after downloading. Reading
headers first (EIS `head.h5`, XRT range requests, JSOC keywords) keeps discovery to a
few MB. Fido remains a fallback (`import eispac.net` registers the EIS client).

Authentication: only JSOC exports need a registered email (`$JSOC_EMAIL`). Keyword
queries, NRL, MSSL and SDAC need nothing.

Pitfalls found while building this (details in `.claude/skills/solar-archive-access`):
TAI timestamps in HMI keywords; silent `im_patch` failures; eispac mask ignored by
`sum_spectra`; EIS pointing conventions; eispac not installable on Python 3.13.
