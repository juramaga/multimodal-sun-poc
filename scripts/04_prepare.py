#!/usr/bin/env python
"""Phase 3b: calibrate raw files into science-ready native-resolution maps.

  EIS  window-integrated intensity for every window + Gaussian fit (int, vel, width)
  AIA  update_pointing + correct_degradation + /exptime   (DN/pix/s)
  HMI  QUALITY check, float32 copy (Gauss)
  XRT  /exptime, saturated pixels -> NaN                  (DN/s)

Writes data/prep/<modality>/<epoch>/*.fits and data/manifest/<epoch>/products.parquet
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from mmsun import eis, prep, xrt  # noqa: E402
from mmsun.cli import load_json, parser, paths_for, setup, write_parquet  # noqa: E402


def main():
    p = parser(__doc__)
    p.add_argument("--ncpu", type=int, default=1, help="processes for EIS Gaussian fits")
    p.add_argument("--no-aiapy-tables", action="store_true",
                   help="skip pointing/degradation tables (offline runs, synthetic tests)")
    args = p.parse_args()
    cfg, epochs = setup(args)

    for ep in epochs:
        P = paths_for(cfg, ep["id"])
        info = load_json(P.manifest / "epoch.json")
        recs = []
        print(f"\n=== {ep['id']} ===")

        data = P.raw("eis") / f"{ep['eis_file']}.data.h5"
        if data.exists():
            out = eis.prepare(data, P.prep("eis"), cfg["eis"]["fit_templates"],
                              cfg["eis"]["integrate_all_windows"], ncpu=args.ncpu)
            recs += [dict(r, modality="EIS", t_obs=info["date_obs"]) for r in out]
            print(f"EIS  {len(out)} maps")

        ptab = None if args.no_aiapy_tables else prep.pointing_table(pd.Timestamp(info["t_mid"]))
        for role in ("pseudoraster", "context"):
            files = sorted((P.raw("aia") / role).glob("*.fits"))
            outdir = P.prep("aia") / role
            outdir.mkdir(exist_ok=True)
            for f in files:
                r = prep.prep_aia(f, outdir, ptab, degradation=not args.no_aiapy_tables)
                recs.append(dict(r, modality="AIA", role=role))
            if files:
                print(f"AIA  {role}: {len(files)} maps   steps: {recs[-1]['steps']}")

        for f in sorted(P.raw("hmi").glob("*.fits")):
            r = prep.prep_hmi(f, P.prep("hmi"))
            recs.append(dict(r, modality="HMI", role="context"))
            print(f"HMI  {Path(r['path']).name}  QUALITY={r['quality']}")

        # only the frames chosen by 02_match: raw/ keeps every frame ever downloaded
        match_file = P.manifest / "match.parquet"
        if match_file.exists():
            m = pd.read_parquet(match_file)
            names = [Path(u).name for u in m.loc[m.modality == "XRT", "ds_or_url"]]
            xrt_files = [P.raw("xrt") / n for n in names if (P.raw("xrt") / n).exists()]
            for n in names:
                if not (P.raw("xrt") / n).exists():
                    print(f"XRT  {n}  MISSING in raw/ (run 03_download.py)")
        else:
            print("XRT  no match.parquet: preparing every file in raw/xrt")
            xrt_files = sorted(P.raw("xrt").glob("*.fits"))
        for f in xrt_files:
            r = xrt.prepare(f, P.prep("xrt"))
            recs.append(dict(r, modality="XRT"))
            print(f"XRT  {Path(r['path']).name}  saturated px: {r['n_saturated']}")

        df = pd.DataFrame(recs)
        write_parquet(df, P.manifest / "products.parquet")
        print(f"products: {df.groupby('modality').size().to_dict()}")


if __name__ == "__main__":
    main()
