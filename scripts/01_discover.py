#!/usr/bin/env python
"""Phases 0-1: metadata-only discovery for each active epoch.

Downloads nothing but: the EIS *.head.h5 file (~0.1-0.5 MB, skipped if bundled),
JSOC keywords (no images), and XRT FITS headers via HTTP range requests.

Outputs (data/manifest/<epoch>/):
  epoch.json              EIS raster summary (times, pointing, FOV, windows, study)
  eis_steps.parquet       per-slit-step time and x position
  cand_aia.parquet        AIA frames (keywords) around the raster
  cand_hmi.parquet        HMI frames (keywords)
  cand_xrt.parquet        XRT Level-1 headers (filter, FOV, exposure) around the raster
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from mmsun import eis, jsoc, xrt  # noqa: E402
from mmsun.cli import parser, paths_for, save_json, setup, write_parquet  # noqa: E402


def main():
    p = parser(__doc__)
    p.add_argument("--skip", nargs="*", default=[], choices=["aia", "hmi", "xrt"],
                   help="skip a modality's discovery (e.g. when an archive is down)")
    args = p.parse_args()
    cfg, epochs = setup(args)
    tol = cfg["tolerance"]

    for ep in epochs:
        P = paths_for(cfg, ep["id"])
        print(f"\n=== {ep['id']} ===")

        # ---- EIS: head file only
        files = eis.fetch(ep["eis_file"], P.raw("eis"), cfg["eis"]["archives"], kinds=("head",),
                          use_bundled=cfg["eis"]["use_bundled_if_available"])
        info = eis.head_info(files["head"])
        box = eis.fov_box(info)
        steps = info.pop("steps")
        info.pop("solar_y")
        info["box_arcsec"] = box
        write_parquet(steps, P.manifest / "eis_steps.parquet")
        save_json(info, P.manifest / "epoch.json")
        print(f"EIS  {info['date_obs']} -> {info['date_end']}  {info['nsteps']} steps  "
              f"FOV x[{box[0]:.0f},{box[1]:.0f}] y[{box[2]:.0f},{box[3]:.0f}]\"  "
              f"study={info['study']}  NOAA {info['noaa']}")

        t0, t1, tmid = info["date_obs"], info["date_end"], info["t_mid"]

        # ---- AIA keywords: pseudo-raster channel over the raster, context around mid
        if "aia" not in args.skip:
            j = cfg["jsoc"]
            pad = pd.Timedelta(seconds=tol["aia_step_s"] + 1)
            a1 = jsoc.query_aia(j["aia_euv_series"], t0 - pad, t1 + pad, [j["aia_pseudoraster_channel"]])
            cpad = pd.Timedelta(seconds=tol["aia_context_s"] + 1)
            a2 = jsoc.query_aia(j["aia_euv_series"], tmid - cpad, tmid + cpad, j["aia_context_channels"])
            aia = pd.concat([a1, a2]).drop_duplicates(subset=["t_rec", "channel"])
            n_all = len(aia)
            if j["require_quality_zero"]:
                aia = aia[aia["quality"].astype(int) == 0]
            write_parquet(aia, P.manifest / "cand_aia.parquet")
            print(f"AIA  {len(aia)} frames with QUALITY==0 (of {n_all})  channels={sorted(aia.channel.unique())}")

        # ---- HMI keywords
        if "hmi" not in args.skip:
            hpad = pd.Timedelta(seconds=tol["hmi_s"] + 1)
            hmi = pd.concat([jsoc.query_hmi(s, tmid - hpad, tmid + hpad) for s in cfg["jsoc"]["hmi_series"]])
            if cfg["jsoc"]["require_quality_zero"]:
                hmi = hmi[hmi["quality"].astype(int) == 0]
            write_parquet(hmi, P.manifest / "cand_hmi.parquet")
            print(f"HMI  {len(hmi)} frames with QUALITY==0")

        # ---- XRT headers via range requests (AR window, then synoptic window if needed)
        if "xrt" not in args.skip:
            xc = cfg["xrt"]
            w = pd.Timedelta(seconds=tol["xrt_ar_s"])
            files = xrt.list_files(t0 - w, t1 + w, xc["level1_archive"])
            if xc.get("allow_synoptic_fallback", True):
                ws = pd.Timedelta(seconds=tol["xrt_synoptic_s"])
                more = xrt.list_files(t0 - ws, t1 + ws, xc["level1_archive"])
                files = pd.concat([files, more]).drop_duplicates("url")
            # nearest first, so the request cap never drops the most useful headers
            files = files.assign(adt=(files.t_name - tmid).abs()).sort_values("adt")
            cache = Path(cfg["_root"]) / xc["header_cache"]
            hdr = xrt.header_table(files, cache, max_requests=xc["max_header_requests"])
            if not hdr.empty and "coverage" not in hdr:
                ok = hdr[hdr.error == ""] if "error" in hdr else hdr
                hdr.loc[ok.index, "coverage"] = [xrt.coverage(r, box) for r in ok.itertuples()]
            write_parquet(hdr, P.manifest / "cand_xrt.parquet")
            if hdr.empty or "filter" not in hdr:
                print("XRT  no Level-1 files found in the search window")
            else:
                good = hdr[(hdr.error == "") & (hdr.imtype.str.lower() == "normal")]
                cov = good[good.coverage >= tol["min_fov_coverage"]]
                print(f"XRT  {len(files)} files listed, {len(good)} normal, "
                      f"{len(cov)} cover >= {tol['min_fov_coverage']:.0%} of the EIS FOV; "
                      f"filters: {cov['filter'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
