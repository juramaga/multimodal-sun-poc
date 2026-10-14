# 01 · Anchor epoch

## Choice

`eis_20210306_064444` — Hinode/EIS Level-1 raster, 2021-03-06 06:44:44 → 06:49:34 UT.

| Property | Value (from the file's own header) |
|---|---|
| Study | `cam_flare_hot1_study` (id 598) |
| Observing title | HOP396 Collaboration w/ Chandrayaan-2 XSM AR12807 |
| Target / science objective | Flare / AR |
| Raster | 25 slit steps × 3.99″, 2″ slit, 120 pixels × 1″ along the slit |
| FOV | 99.8″ × 120″ |
| Pointing (eispac convention) | xcen ≈ +3″, ycen ≈ −197″ (≈ 0.21 R☉ from disk centre, μ ≈ 0.98) |
| Exposure | 10 s per step, ~11.6 s cadence |
| Windows | Ca XV 181.9, Fe XII 186.75, **Fe XII 192.41**, Ar XIV 194.2, Ca XV 201.0, Fe XVII 254.95, S XIII 256.95, **Fe XXIII 263.3**, Fe XIV 270.51 |
| Files | `eis_20210306_064444.data.h5` (3.4 MB) + `.head.h5` (0.13 MB), NRL archive; also bundled in `eispac/data/test/` |

The same study ran as a sequence of 9 rasters (06:25:12, 06:30:05, 06:34:58,
06:39:51, **06:44:44**, 06:49:37, 06:54:30, 06:59:23, 07:04:16 UT), all 3.4 MB.
Later that day four large rasters (18:12, 19:40, 21:17, 22:54 UT; ~150 MB each) were taken.

## Why it maximises the chance of a full match

1. **XRT co-pointing.** It is a flare-watch HOP coordinated with an X-ray spectrometer
   (Chandrayaan-2/XSM). Hinode's instruments share the spacecraft pointing, so XRT
   imaging during the campaign sees the same AR (a published XSM–Hinode study of
   AR 12759 used XRT Al_poly and Be_thin, the two filters preferred here). The
   nominal ~06 UT XRT synoptic (full disk) is a second chance.
2. **Short raster.** 4 m 50 s end to end → every modality can be matched within one
   AIA cadence per step, and XRT within a few minutes.
3. **Disk-centre active region.** Little projection, every pixel on disk, unambiguous
   Carrington coordinates, AR 12807 also tracked by IRIS (HOP 307) the same day.
4. **Temperature overlap.** Fe XII 192 ↔ AIA 193 for coalignment; Ca XV, Fe XVII and
   Fe XXIII cover the 3–15 MK range where XRT is most sensitive.
5. **Zero-risk EIS access.** The raster is bundled with eispac, so the anchor modality
   is always available and the offline test uses real EIS data.

## What is verified vs. expected

* Verified from the archive listing and the file itself: the raster, its sequence,
  pointing, timing, windows, title (checked 2026-10-14).
* Expected but verified only on a networked run: XRT frames within 5 min that cover the
  FOV, and AIA/HMI `QUALITY == 0` frames. `scripts/01_discover.py` prints this; if XRT
  is missing, run the sequence (`--epochs seq`); its first raster (06:25 UT) lies
  within ~30 min of the nominal 06 UT synoptic, which the synoptic fallback accepts.
