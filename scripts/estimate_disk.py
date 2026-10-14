#!/usr/bin/env python
"""Disk budget: a-priori estimate from the config, and measured usage if data exist.

The estimate is per epoch and assumes the anchor geometry (EIS 100"x120" FOV):
  python scripts/estimate_disk.py                 # estimate for active epochs
  python scripts/estimate_disk.py --epochs all    # the whole 9-raster sequence
  python scripts/estimate_disk.py --measure       # also report what is on disk now
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mmsun.cli import parser, paths_for, setup  # noqa: E402

MB = 1e6


def per_epoch(cfg, fov=(100.0, 120.0), raster_s=290.0, xrt_px=(512, 512), xrt_full=False):
    j = cfg["jsoc"]
    cm, pm = 2 * j["context_margin_arcsec"], 2 * j["pseudoraster_margin_arcsec"]
    aia_px = lambda w, h: (w / 0.6) * (h / 0.6)  # noqa: E731
    n_pr = int(raster_s / 12) + 2
    ctx_px = aia_px(fov[0] + cm, fov[1] + cm)
    pr_px = aia_px(fov[0] + pm, fov[1] + pm)
    hmi_px = ((fov[0] + cm) / 0.504) * ((fov[1] + cm) / 0.504) * len(j["hmi_series"])
    nx = cfg["xrt"]["max_filters"]
    xrt_bytes = nx * (2048 * 2048 if xrt_full else xrt_px[0] * xrt_px[1]) * 4
    raw = {
        "EIS (L1 HDF5 pair)": 3.6 * MB,
        f"AIA pseudo-raster ({n_pr} x 193 cutouts, int16)": n_pr * pr_px * 2,
        f"AIA context ({len(j['aia_context_channels'])} cutouts, int16)": len(j["aia_context_channels"]) * ctx_px * 2,
        "HMI magnetogram cutout": hmi_px * 4,
        f"XRT Level-1 ({nx} filters, {'full-frame' if xrt_full else 'partial frame'})": xrt_bytes,
        "XRT header cache (range requests)": 3 * MB,
    }
    prep = {  # float32 copies at native resolution + EIS fits
        "EIS maps (9 windows + 3 fit maps)": 1.0 * MB,
        "AIA prepped (float32)": (n_pr * pr_px + len(j["aia_context_channels"]) * ctx_px) * 4,
        "HMI prepped": hmi_px * 4,
        "XRT prepped": xrt_bytes,
    }
    aligned = {"stack.fits + pixels.parquet + native.parquet + diag": 7 * MB, "manifests": 0.2 * MB}
    return raw, prep, aligned


def show(title, d):
    tot = sum(d.values())
    print(f"  {title:9s} {tot / MB:8.1f} MB")
    for k, v in d.items():
        print(f"      {k:55s} {v / MB:7.1f}")
    return tot


def main():
    p = parser(__doc__)
    p.add_argument("--measure", action="store_true")
    args = p.parse_args()
    cfg, epochs = setup(args)

    for label, kw in (("typical (XRT partial frames)", {}), ("worst case (XRT full-frame)", {"xrt_full": True})):
        print(f"\nPer-epoch estimate, {label}:")
        raw, prep, aligned = per_epoch(cfg, **kw)
        t = show("raw", raw) + show("prep", prep) + show("aligned", aligned)
        print(f"  {'TOTAL':9s} {t / MB:8.1f} MB per epoch   x {len(epochs)} active epoch(s) = "
              f"{t * len(epochs) / MB:.0f} MB")
    print("\nSoftware (one-off): conda env 'mmsun' ~1.5-3 GB; sunpy/aiapy calibration caches < 20 MB.")
    print("Recommendation: keep >= 5 GB free for the anchor run, ~6 GB for the full 9-raster sequence.")

    if args.measure:
        print("\nMeasured on disk:")
        root = paths_for(cfg, epochs[0]["id"]).root
        for sub in ("raw", "prep", "aligned", "manifest", "cache"):
            d = root / sub
            if d.exists():
                size = sum(f.stat().st_size for f in d.rglob("*") if f.is_file())
                print(f"  {sub:9s} {size / MB:8.1f} MB")


if __name__ == "__main__":
    main()
