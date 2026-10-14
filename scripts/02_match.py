#!/usr/bin/env python
"""Phase 2: select the minimal co-temporal, co-spatial file set for each epoch.

Reads data/manifest/<epoch>/cand_*.parquet and writes
  match.parquet     one row per file/export to retrieve (role, dt, coverage, size estimate)
  eis_step_map.parquet  EIS step -> AIA pseudo-raster frame, with dt per step
  match_report.json     per-modality verdict; exit code 3 if any modality is missing
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from mmsun import match  # noqa: E402
from mmsun.cli import load_json, parser, paths_for, save_json, setup, write_parquet  # noqa: E402

HMI_SEGMENT = {"hmi.M_45s": "magnetogram", "hmi.Ic_45s": "continuum", "hmi.M_720s": "magnetogram"}


def read(path):
    return pd.read_parquet(path) if Path(path).exists() else pd.DataFrame()


def main():
    args = parser(__doc__).parse_args()
    cfg, epochs = setup(args)
    tol, j = cfg["tolerance"], cfg["jsoc"]
    any_missing = False

    for ep in epochs:
        P = paths_for(cfg, ep["id"])
        info = load_json(P.manifest / "epoch.json")
        t0, t1, tmid = (pd.Timestamp(info[k]) for k in ("date_obs", "date_end", "t_mid"))
        box = tuple(info["box_arcsec"])
        steps = pd.read_parquet(P.manifest / "eis_steps.parquet")
        steps["t_step"] = pd.to_datetime(steps["t_step"], format="ISO8601")
        rows, verdict = [], {}

        # EIS itself
        rows.append(dict(epoch=ep["id"], modality="EIS", role="anchor", channel="all_windows",
                         t_obs=t0, t_rec="", dt_s=0.0, steps=f"0-{len(steps) - 1}",
                         ds_or_url=ep["eis_file"], coverage=1.0, est_bytes=3_600_000, series="eis_l1_hdf5"))
        verdict["EIS"] = "ok"

        # AIA
        aia = read(P.manifest / "cand_aia.parquet")
        if not aia.empty:
            aia["t_obs"] = pd.to_datetime(aia["t_obs"], format="ISO8601")
            aia["channel"] = aia["channel"].astype(float).astype(int)
            step_map, frames = match.aia_pseudoraster(aia, steps, j["aia_pseudoraster_channel"], tol["aia_step_s"])
            write_parquet(step_map, P.manifest / "eis_step_map.parquet")
            pr_box = tuple(b + s * cfg["jsoc"]["pseudoraster_margin_arcsec"]
                           for b, s in zip(box, (-1, 1, -1, 1)))
            for f in frames.itertuples():
                rows.append(dict(epoch=ep["id"], modality="AIA", role="pseudoraster",
                                 channel=str(j["aia_pseudoraster_channel"]), t_obs=f.t_obs, t_rec=f.t_rec,
                                 dt_s=f.dt_s, steps=f.steps, ds_or_url=j["aia_euv_series"], coverage=1.0,
                                 est_bytes=match.est_cutout_bytes(pr_box, 0.6), series=j["aia_euv_series"]))
            ctx = match.nearest_per_channel(aia, tmid, j["aia_context_channels"], tol["aia_context_s"])
            cbox = tuple(b + s * cfg["jsoc"]["context_margin_arcsec"] for b, s in zip(box, (-1, 1, -1, 1)))
            for r in ctx.itertuples():
                rows.append(dict(epoch=ep["id"], modality="AIA", role="context", channel=str(r.channel),
                                 t_obs=r.t_obs, t_rec=r.t_rec, dt_s=r.dt_s, steps="",
                                 ds_or_url=j["aia_euv_series"], coverage=1.0,
                                 est_bytes=match.est_cutout_bytes(cbox, 0.6), series=j["aia_euv_series"]))
            n_ok = int(step_map["ok"].sum()) if not step_map.empty else 0
            verdict["AIA"] = (f"ok: {n_ok}/{len(steps)} steps matched, {len(frames)} frames; "
                              f"{len(ctx)}/{len(j['aia_context_channels'])} context channels")
            if n_ok == 0 or ctx.empty:
                verdict["AIA"] = "MISSING: " + verdict["AIA"]
        else:
            verdict["AIA"] = "MISSING: no candidates (run 01_discover.py)"

        # HMI
        hmi = read(P.manifest / "cand_hmi.parquet")
        if not hmi.empty:
            hmi["t_obs"] = pd.to_datetime(hmi["t_obs"], format="ISO8601")
            for series in j["hmi_series"]:
                sub = match.nearest_per_channel(hmi.assign(channel=hmi["series"]), tmid, [series], tol["hmi_s"])
                for r in sub.itertuples():
                    cbox = tuple(b + s * cfg["jsoc"]["context_margin_arcsec"] for b, s in zip(box, (-1, 1, -1, 1)))
                    rows.append(dict(epoch=ep["id"], modality="HMI", role="context",
                                     channel=HMI_SEGMENT.get(series, "image"), t_obs=r.t_obs, t_rec=r.t_rec,
                                     dt_s=r.dt_s, steps="", ds_or_url=series, coverage=1.0,
                                     est_bytes=match.est_cutout_bytes(cbox, 0.504, 4), series=series))
            verdict["HMI"] = "ok" if any(r["modality"] == "HMI" for r in rows) else "MISSING: none within tolerance"
        else:
            verdict["HMI"] = "MISSING: no candidates"

        # XRT
        hdr = read(P.manifest / "cand_xrt.parquet")
        if not hdr.empty and "t_obs" in hdr:
            hdr["t_obs"] = pd.to_datetime(hdr["t_obs"], format="ISO8601", errors="coerce")
            for c in ("exptime", "fovx", "fovy", "xcen", "ycen"):
                hdr[c] = pd.to_numeric(hdr[c], errors="coerce")
            hdr["error"] = hdr["error"].replace("nan", "")
        chosen, mode = match.select_xrt(hdr, box, t0, t1, cfg, tmid)
        for r in chosen.itertuples():
            rows.append(dict(epoch=ep["id"], modality="XRT", role=mode, channel=r.filter, t_obs=r.t_obs,
                             t_rec="", dt_s=float(r.dt_s), steps="", ds_or_url=r.url,
                             coverage=float(r.coverage), est_bytes=int(r.bytes_est), series="xrt_l1_sdac"))
        verdict["XRT"] = (f"ok ({mode}): " + ", ".join(f"{r.filter} dt={r.dt_s:+.0f}s (mid {r.dt_mid_s:+.0f}s)"
                                                         for r in chosen.itertuples())) if mode != "none" \
            else "MISSING: no frame covers the EIS FOV within tolerance"

        m = pd.DataFrame(rows, columns=match.MATCH_COLUMNS)
        m["dt_mid_s"] = (pd.to_datetime(m.t_obs) - tmid).dt.total_seconds()
        write_parquet(m, P.manifest / "match.parquet")
        missing = [k for k, v in verdict.items() if v.startswith("MISSING")]
        any_missing |= bool(missing)
        save_json(dict(epoch=ep["id"], verdict=verdict, n_files=len(m),
                       est_mb=round(m.est_bytes.sum() / 1e6, 1)), P.manifest / "match_report.json")
        print(f"\n=== {ep['id']} === ({len(m)} files, ~{m.est_bytes.sum() / 1e6:.0f} MB raw)")
        for k, v in verdict.items():
            print(f"  {k:4s} {v}")
        if missing:
            print("  -> incomplete epoch. Try the next raster in the sequence: "
                  "python scripts/01_discover.py --epochs seq")

    sys.exit(3 if any_missing else 0)


if __name__ == "__main__":
    main()
