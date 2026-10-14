---
name: solar-coalignment
description: Conventions and procedure for coaligning Hinode EIS/XRT with SDO AIA/HMI in this project - pseudo-rasters, shift sign conventions, solar-rotation handling, reference-channel choice and how to judge success. Use when touching src/mmsun/coalign.py, src/mmsun/grid.py or scripts/05_coalign.py, or interpreting shifts.json.
---

# Coalignment in the Multimodal Sun PoC

## Conventions (do not change without updating every caller)

* All positions are helioprojective arcsec as seen from Earth/SDO. Hinode's
  low-Earth-orbit parallax (< 0.01″) is ignored.
* A pointing correction `d = (dx, dy)` means **true position = header position + d**.
  AIA is the reference (`d = 0`; its pointing is fixed by `aiapy.update_pointing`).
  HMI is assumed co-registered with AIA at the JSOC level-1 product level.
* To sample map *m* (correction `d`, taken at `t_m`) for reference time `t_ref` at true
  position `X`: header coordinate = `X + r(t_ref → t_m) − d`, where `r` is the
  solar-rotation displacement of the FOV centre (`coalign.rotation_offset`). At disk
  centre this is ~0.15″/min, so it matters for XRT (minutes) but not within an EIS step.
* `t_ref` = EIS raster mid-time. Within the raster, each EIS column keeps its own time.

## Procedure

1. **EIS → AIA pseudo-raster.** For each EIS column (slit step), take the AIA frame
   nearest that step's time, blur it to ~3″ FWHM, and sample it at the EIS pixel
   positions. Compare with the EIS window that contains `eis.coalign_window`
   (Fe XII 192.39 Å ↔ AIA 193). Score = Pearson r of log intensities.
2. **XRT → AIA.** Compare on a 2″ grid over the EIS FOV + 30″; try AIA 335, 94, 211 and
   keep the best r (XRT filters respond to hotter plasma than 193).
3. Search: 4″ steps over ±`max_shift_arcsec`, then 1″ over ±4″, then 0.25″ over ±1″.
4. Apply: every EIS window inherits the same `d` (eispac already applies the relative
   CCD offsets between windows); XRT `d` is per frame.

## Judging the result

* Good: r_after ≳ 0.8 and r_after − r_before > 0 by a clear margin; |d| well inside the
  search box. Typical magnitudes: EIS 3–15″, XRT 1–5″.
* Suspicious: shift at the search boundary, r_after < `coalign.min_correlation`
  (flagged), or EIS shifts that differ by > 5″ between consecutive rasters of the same
  sequence (pointing should drift slowly).
* Remember that the long-wavelength EIS windows (Fe XVII, S XIII, Fe XXIII, Fe XIV) cover
  a y-range offset by the CCD offset (~17″) from the short-wavelength ones: blank rows
  at the top/bottom of those layers on the common grid are expected, not a bug.

## Tests

`make test` injects EIS d = (+7.5, −4.0)″ and XRT d = (−3.0, +5.0)″ into synthetic
AIA/XRT images built around the real EIS raster and requires recovery to < 1.5″.
