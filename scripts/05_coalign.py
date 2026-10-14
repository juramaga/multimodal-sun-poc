#!/usr/bin/env python
"""Phase 4: coalign EIS and XRT to AIA, resample everything onto one grid, write tables.

Outputs (data/aligned/<epoch>/):
  shifts.json         pointing corrections (arcsec), correlation before/after, references used
  coalign_diag.npz    EIS target, AIA pseudo-raster (before/after), correlation surfaces
  stack.fits          one HDU per channel on the common grid (EIS FOV, t_ref = raster mid-time)
  pixels.parquet      one row per grid pixel: Tx, Ty, lon_carr, lat_carr, mu, all channels, dt_*
  native.parquet      long table of native pixels inside the EIS FOV, all modalities, with lon/lat
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mmsun import coalign, grid  # noqa: E402
from mmsun.cli import load_json, parser, paths_for, save_json, setup, write_parquet  # noqa: E402


def load_map(path):
    import sunpy.map

    return sunpy.map.Map(path)


def tstamp(m) -> pd.Timestamp:
    return pd.Timestamp(m.date.datetime)


def main():
    args = parser(__doc__).parse_args()
    cfg, epochs = setup(args)
    cc, gc = cfg["coalign"], cfg["grid"]

    for ep in epochs:
        P = paths_for(cfg, ep["id"])
        info = load_json(P.manifest / "epoch.json")
        t_ref = pd.Timestamp(info["t_mid"])
        box = tuple(info["box_arcsec"])
        steps = pd.read_parquet(P.manifest / "eis_steps.parquet")
        steps["t_step"] = pd.to_datetime(steps["t_step"], format="ISO8601")
        prods = pd.read_parquet(P.manifest / "products.parquet")
        print(f"\n=== {ep['id']} ===  t_ref={t_ref}")
        shifts, diag = {}, {}

        # ------------------------------------------------ EIS vs AIA pseudo-raster
        eis_p = prods[prods.modality == "EIS"]
        cw = float(cfg["eis"]["coalign_window"])
        win = next(w for w in info["windows"] if w["wvl_min"] <= cw <= w["wvl_max"])
        eis_ref_row = eis_p[eis_p["product"].str.startswith(f"eis_win{win['iwin']:02d}_")].iloc[0]
        eis_ref = load_map(eis_ref_row.path)
        pr = prods[(prods.modality == "AIA") & (prods.get("role") == "pseudoraster")]
        aia_frames = sorted(((pd.Timestamp(r.t_obs), load_map(r.path)) for r in pr.itertuples()),
                            key=lambda x: x[0])
        if aia_frames:
            res = coalign.eis_vs_aia_pseudoraster(eis_ref, steps, aia_frames, cc["max_shift_arcsec"])
            d_eis = (res["dx"], res["dy"])
            shifts["EIS"] = dict(dx=res["dx"], dy=res["dy"], r_after=res["r"], r_before=res["r_before"],
                                 reference=f"AIA {cfg['jsoc']['aia_pseudoraster_channel']} pseudo-raster",
                                 eis_window=win["line_id"], n_aia_frames=len(aia_frames),
                                 flagged=res["r"] < cc["min_correlation"])
            diag.update(eis_target=np.asarray(eis_ref.data, float), pr_after=res["pseudoraster"],
                        pr_before=res["pseudoraster_unshifted"], eis_surface=res["surface"],
                        column_frame=res["column_frame"])
            print(f"EIS  shift dx={res['dx']:+.2f}\" dy={res['dy']:+.2f}\"  "
                  f"r {res['r_before']:.3f} -> {res['r']:.3f}  ({len(aia_frames)} AIA frames)")
        else:
            d_eis = (0.0, 0.0)
            shifts["EIS"] = dict(dx=0.0, dy=0.0, note="no AIA pseudo-raster; header pointing kept", flagged=True)
            print("EIS  no pseudo-raster frames: keeping header pointing")

        # common-grid footprint = the coalignment window's own pixel footprint, corrected
        ex, ey = coalign.world_arcsec(eis_ref)
        sx, sy = info["x_scale"], info["y_scale"]
        true_box = (ex.min() - sx / 2 + d_eis[0], ex.max() + sx / 2 + d_eis[0],
                    ey.min() - sy / 2 + d_eis[1], ey.max() + sy / 2 + d_eis[1])

        # ------------------------------------------------ XRT vs AIA context
        ctx = prods[(prods.modality == "AIA") & (prods.get("role") == "context")].copy()
        ctx["channel"] = ctx["channel"].astype(float).astype(int)
        aia_ctx = {int(r.channel): load_map(r.path) for r in ctx.itertuples()}
        xrt_maps = {}
        mbox = (true_box[0] - 30, true_box[1] + 30, true_box[2] - 30, true_box[3] + 30)
        for r in prods[prods.modality == "XRT"].itertuples():
            xm = load_map(r.path)
            best = None
            for ch in cc["xrt_reference_candidates"]:
                if ch not in aia_ctx:
                    continue
                res = coalign.map_vs_reference(xm, aia_ctx[ch], mbox, tstamp(aia_ctx[ch]), tstamp(xm),
                                               max_shift=cc["max_shift_arcsec"])
                if best is None or res["r"] > best[1]["r"]:
                    best = (ch, res)
            key = f"XRT_{r.filter}"
            if best is None:
                shifts[key] = dict(dx=0.0, dy=0.0, note="no AIA reference", flagged=True)
                xrt_maps[r.filter] = (xm, (0.0, 0.0))
                continue
            ch, res = best
            shifts[key] = dict(dx=res["dx"], dy=res["dy"], r_after=res["r"], r_before=res["r_before"],
                               reference=f"AIA {ch}", rotation_offset=res["rotation_offset"],
                               flagged=res["r"] < cc["min_correlation"])
            diag[f"xrt_surface_{r.filter}"] = res["surface"]
            xrt_maps[r.filter] = (xm, (res["dx"], res["dy"]))
            print(f"XRT  {r.filter:8s} vs AIA {ch}: dx={res['dx']:+.2f}\" dy={res['dy']:+.2f}\"  "
                  f"r {res['r_before']:.3f} -> {res['r']:.3f}")
        save_json(shifts, P.aligned / "shifts.json")
        np.savez_compressed(P.aligned / "coalign_diag.npz", **diag)

        # ------------------------------------------------ common grid
        pix = gc["pixel_arcsec"]
        gbox = tuple(b + s * gc.get("margin_arcsec", 0) for b, s in zip(true_box, (-1, 1, -1, 1)))
        gx, gy = grid.make_grid(gbox, pix)
        layers, dts, units = {}, {}, {}

        for r in eis_p.itertuples():
            m = load_map(r.path)
            name = r.product.replace("eis_", "eis_")
            layers[name] = grid.sample_on_grid(m, gx, gy, t_ref, t_ref, d_eis, pix)
            dts[name], units[name] = 0.0, r.unit

        for ch, m in aia_ctx.items():
            name = f"aia_{ch}"
            layers[name] = grid.sample_on_grid(m, gx, gy, t_ref, tstamp(m), (0, 0), pix)
            dts[name], units[name] = (tstamp(m) - t_ref).total_seconds(), "DN/pix/s"

        if aia_frames:
            # pseudo-raster on the grid: each grid column uses the AIA frame of the nearest EIS step
            # column x positions in the EIS map's own WCS (step index == column index)
            x_true = coalign.world_arcsec(eis_ref)[0][0] + d_eis[0]
            frame_of_step = diag["column_frame"]
            pr_grid = np.full(gx.shape, np.nan)
            dt_grid = np.full(gx.shape, np.nan)
            st = steps.set_index("step")["t_step"]
            for jcol in range(gx.shape[1]):
                k = int(np.argmin(np.abs(x_true - gx[0, jcol])))
                t_aia, m = aia_frames[frame_of_step[k]]
                pr_grid[:, jcol] = grid.sample_on_grid(m, gx[:, jcol:jcol + 1], gy[:, jcol:jcol + 1],
                                                       t_aia, t_aia, (0, 0), pix)[:, 0]
                dt_grid[:, jcol] = (t_aia - st.loc[k]).total_seconds()
            ch = cfg["jsoc"]["aia_pseudoraster_channel"]
            layers[f"aia_{ch}_pseudoraster"] = pr_grid
            units[f"aia_{ch}_pseudoraster"] = "DN/pix/s"
            dts[f"aia_{ch}_pseudoraster"] = float(np.nanmax(np.abs(dt_grid)))

        for r in prods[prods.modality == "HMI"].itertuples():
            m = load_map(r.path)
            name = f"hmi_{r.channel}"
            layers[name] = grid.sample_on_grid(m, gx, gy, t_ref, tstamp(m), (0, 0), pix)
            dts[name], units[name] = (tstamp(m) - t_ref).total_seconds(), r.unit

        for filt, (m, d) in xrt_maps.items():
            name = f"xrt_{filt.replace('/', '+')}"
            layers[name] = grid.sample_on_grid(m, gx, gy, t_ref, tstamp(m), d, pix)
            dts[name], units[name] = (tstamp(m) - t_ref).total_seconds(), "DN/s"

        # ------------------------------------------------ write stack.fits
        from astropy.io import fits

        hdr = grid.grid_header(gx, gy, pix, t_ref)
        hdul = [fits.PrimaryHDU(header=fits.Header({k: v for k, v in hdr.items() if k != "keycomments"}))]
        hdul[0].header["EPOCH"] = ep["id"]
        hdul[0].header["T_REF"] = t_ref.isoformat()
        for name, arr in layers.items():
            h = fits.Header({k: v for k, v in hdr.items() if k != "keycomments"})
            h["EXTNAME"] = name[:68]
            h["BUNIT"] = str(units.get(name, ""))[:68]
            h["DT_S"] = (float(dts[name]), "t_obs - t_ref [s] (max |dt| for pseudo-raster)")
            hdul.append(fits.ImageHDU(arr.astype("float32"), header=h))
        fits.HDUList(hdul).writeto(P.aligned / "stack.fits", overwrite=True)

        # ------------------------------------------------ pixels.parquet (HATS-ready)
        lon, lat, mu = grid.carrington(gx, gy, t_ref)
        iy, ix = np.mgrid[0:gx.shape[0], 0:gx.shape[1]]
        tab = pd.DataFrame({"epoch": ep["id"], "t_ref": t_ref, "ix": ix.ravel(), "iy": iy.ravel(),
                            "tx": gx.ravel(), "ty": gy.ravel(), "lon_carr": lon.ravel(),
                            "lat_carr": lat.ravel(), "mu": mu.ravel()})
        for name, arr in layers.items():
            tab[name] = arr.ravel().astype("float32")
            tab[f"dt_{name}"] = np.float32(dts[name])
        write_parquet(tab, P.aligned / "pixels.parquet")

        # ------------------------------------------------ native.parquet (for cross-matching)
        parts = []
        for r in eis_p.itertuples():
            t = grid.native_table(load_map(r.path), true_box, t_ref, t_ref, d_eis, "value")
            parts.append(t.assign(modality="EIS", channel=r.product))
        for ch, m in aia_ctx.items():
            parts.append(grid.native_table(m, true_box, t_ref, tstamp(m), (0, 0), "value")
                         .assign(modality="AIA", channel=f"aia_{ch}"))
        for r in prods[prods.modality == "HMI"].itertuples():
            m = load_map(r.path)
            parts.append(grid.native_table(m, true_box, t_ref, tstamp(m), (0, 0), "value")
                         .assign(modality="HMI", channel=f"hmi_{r.channel}"))
        for filt, (m, d) in xrt_maps.items():
            parts.append(grid.native_table(m, true_box, t_ref, tstamp(m), d, "value")
                         .assign(modality="XRT", channel=f"xrt_{filt}"))
        nat = pd.concat([p for p in parts if not p.empty], ignore_index=True)
        for c in ("tx", "ty", "lon_carr", "lat_carr", "mu", "value"):
            nat[c] = nat[c].astype("float32")
        nat["epoch"] = ep["id"]
        write_parquet(nat, P.aligned / "native.parquet")

        print(f"grid {gx.shape[1]}x{gx.shape[0]} @ {pix}\"  channels: {len(layers)}  "
              f"native rows: {len(nat)}  -> {P.aligned}")


if __name__ == "__main__":
    main()
