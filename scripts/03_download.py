#!/usr/bin/env python
"""Phase 3a: retrieve exactly the files listed in match.parquet (idempotent).

  EIS  Level-1 HDF5 data+head pair (bundled copy used when available)
  AIA  JSOC export with server-side im_patch cutouts:
         - pseudo-raster channel: one export covering the raster's T_REC range
         - context channels: one export of one T_REC slot, all channels
  HMI  JSOC export, same box as the AIA context set
  XRT  Level-1 FITS from SDAC (only the selected filters/frames)

Use --dry-run to print what would be fetched and the estimated size.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from mmsun import eis, jsoc  # noqa: E402
from mmsun.cli import die, load_json, parser, paths_for, save_json, setup  # noqa: E402
from mmsun.config import jsoc_email  # noqa: E402
from mmsun.netutil import download  # noqa: E402

SEGMENT = {"hmi.M_45s": "magnetogram", "hmi.Ic_45s": "continuum", "hmi.M_720s": "magnetogram"}


def jsoc_time(t_rec: str) -> str:
    return str(t_rec).strip()


def main():
    p = parser(__doc__)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--only", nargs="*", default=["EIS", "AIA", "HMI", "XRT"])
    args = p.parse_args()
    cfg, epochs = setup(args)
    j = cfg["jsoc"]
    email = jsoc_email(cfg)
    if not args.dry_run and ({"AIA", "HMI"} & set(args.only)) and not email:
        die(f"set ${j['email_env']} to an email registered at "
            "http://jsoc.stanford.edu/ajax/register_email.html (needed for JSOC exports)")

    for ep in epochs:
        P = paths_for(cfg, ep["id"])
        m = pd.read_parquet(P.manifest / "match.parquet")
        info = load_json(P.manifest / "epoch.json")
        box = tuple(info["box_arcsec"])
        tmid = pd.Timestamp(info["t_mid"])
        cbox = tuple(b + s * j["context_margin_arcsec"] for b, s in zip(box, (-1, 1, -1, 1)))
        pbox = tuple(b + s * j["pseudoraster_margin_arcsec"] for b, s in zip(box, (-1, 1, -1, 1)))
        print(f"\n=== {ep['id']} === {len(m)} rows, est. {m.est_bytes.sum() / 1e6:.1f} MB")
        log = []

        if args.dry_run:
            print(m.groupby(["modality", "role"]).agg(n=("channel", "size"),
                                                       mb=("est_bytes", lambda x: x.sum() / 1e6)).round(2))
            continue

        if "EIS" in args.only:
            got = eis.fetch(ep["eis_file"], P.raw("eis"), cfg["eis"]["archives"],
                            use_bundled=cfg["eis"]["use_bundled_if_available"])
            log += [dict(modality="EIS", path=str(v)) for v in got.values()]
            print(f"EIS  {', '.join(v.name for v in got.values())}")

        if "AIA" in args.only:
            a = m[m.modality == "AIA"]
            pr = a[a.role == "pseudoraster"].sort_values("t_obs")
            if not pr.empty and not any(P.raw("aia").glob("pseudoraster/*.fits")):
                first, last = jsoc_time(pr.t_rec.iloc[0]), pr.t_obs.iloc[-1]
                dur = int((pd.Timestamp(last) - pd.Timestamp(pr.t_obs.iloc[0])).total_seconds()) + 13
                ds = f"{j['aia_euv_series']}[{first}/{dur}s@12s][{j['aia_pseudoraster_channel']}]{{image}}"
                print(f"AIA  export {ds}")
                files = jsoc.export_cutout(ds, pbox, tmid, email, P.raw("aia") / "pseudoraster",
                                           j["server_side_lev15"], j["max_wait_s"])
                log += [dict(modality="AIA", role="pseudoraster", path=str(f)) for f in files]
            ctx = a[a.role == "context"]
            if not ctx.empty and not any(P.raw("aia").glob("context/*.fits")):
                # one T_REC slot for all channels keeps the context set simultaneous
                t_rec = ctx.t_rec.mode().iloc[0]
                chans = ",".join(sorted(set(ctx.channel), key=int))
                ds = f"{j['aia_euv_series']}[{jsoc_time(t_rec)}][{chans}]{{image}}"
                print(f"AIA  export {ds}")
                files = jsoc.export_cutout(ds, cbox, tmid, email, P.raw("aia") / "context",
                                           j["server_side_lev15"], j["max_wait_s"])
                log += [dict(modality="AIA", role="context", path=str(f)) for f in files]

        if "HMI" in args.only:
            for r in m[m.modality == "HMI"].itertuples():
                dest = P.raw("hmi")
                if any(dest.glob(f"*{SEGMENT.get(r.series, '')}*.fits")):
                    continue
                ds = f"{r.series}[{jsoc_time(r.t_rec)}]{{{SEGMENT.get(r.series, 'magnetogram')}}}"
                print(f"HMI  export {ds}")
                files = jsoc.export_cutout(ds, cbox, tmid, email, dest, False, j["max_wait_s"])
                log += [dict(modality="HMI", path=str(f)) for f in files]

        if "XRT" in args.only:
            for r in m[m.modality == "XRT"].itertuples():
                f = download(r.ds_or_url, P.raw("xrt") / Path(r.ds_or_url).name)
                log.append(dict(modality="XRT", channel=r.channel, path=str(f)))
                print(f"XRT  {f.name}  ({r.channel}, dt={r.dt_s:+.0f}s)")

        size = sum(f.stat().st_size for mod in ("eis", "aia", "hmi", "xrt")
                   for f in P.raw(mod).rglob("*") if f.is_file())
        save_json(dict(epoch=ep["id"], files=log, raw_bytes=size), P.manifest / "download_log.json")
        print(f"raw data on disk for this epoch: {size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
