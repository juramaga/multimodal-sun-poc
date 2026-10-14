#!/usr/bin/env python
"""Generate notebooks/01_summary.ipynb (kept as code so it diffs cleanly and can be regenerated)."""

from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
nb = nbf.v4.new_notebook()
C = []


def md(s):
    C.append(nbf.v4.new_markdown_cell(s.strip()))


def code(s):
    C.append(nbf.v4.new_code_cell(s.strip()))


md("""
# Multimodal Sun — proof-of-concept summary (Phases 0–4)

One EIS raster is the anchor; AIA, HMI and XRT are matched to it in time and space,
calibrated, coaligned and resampled onto a common grid. This notebook reads only the
pipeline outputs (`data/manifest`, `data/prep`, `data/aligned`) and runs no downloads.

Environment variables: `MMSUN_EPOCH` (default `anchor`) picks the epoch;
`MMSUN_CONFIG` points at an alternative `pipeline.yaml` (the offline test uses this).
""")

code("""
import os, sys, json, warnings
from pathlib import Path
ROOT = Path.cwd() if (Path.cwd() / "config").exists() else Path.cwd().parent
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd, matplotlib.pyplot as plt
import sunpy, sunpy.map, sunpy.visualization.colormaps  # registers sdoaia*/hinodexrt colormaps
sunpy.log.setLevel("WARNING")
from astropy.io import fits
from mmsun.config import load_config, active_epochs, paths_for
from mmsun.coalign import world_arcsec

cfg = load_config(os.environ.get("MMSUN_CONFIG"))
ep = active_epochs(cfg, [os.environ.get("MMSUN_EPOCH", "anchor")])[0]
P = paths_for(cfg, ep["id"])
info = json.loads((P.manifest / "epoch.json").read_text())
prods = pd.read_parquet(P.manifest / "products.parquet")
shifts = json.loads((P.aligned / "shifts.json").read_text())
t_ref = pd.Timestamp(info["t_mid"])
plt.rcParams.update({"figure.dpi": 110, "image.origin": "lower", "axes.titlesize": 9,
                     "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7})

def extent(m, d=(0, 0)):
    tx, ty = world_arcsec(m)
    return [tx.min() + d[0], tx.max() + d[0], ty.min() + d[1], ty.max() + d[1]]

def box_xy(b):
    return [b[0], b[1], b[1], b[0], b[0]], [b[2], b[2], b[3], b[3], b[2]]

def robust(a, lo=1, hi=99.5):
    a = np.asarray(a, float); f = a[np.isfinite(a)]
    return (np.percentile(f, lo), np.percentile(f, hi)) if f.size else (0, 1)

print(f"epoch {ep['id']}   {info['date_obs']} -> {info['date_end']}   t_ref = {t_ref}")
print(f"EIS study {info['study']} ({info['study_id']}), {info['nsteps']} steps, slit {info['slit']}, "
      f"FOV {info['fovx']:.0f}\\"x{info['fovy']:.0f}\\" at ({info['xcen']:.0f}\\", {info['ycen']:.0f}\\"), NOAA {info['noaa']}")
print(f"title: {info['obstitle']}")
rep = P.manifest / "match_report.json"
if rep.exists():
    for k, v in json.loads(rep.read_text())["verdict"].items():
        print(f"  {k:4s} {v}")
""")

md("""
## 1. What is in the sample

Every product with its observation time and offset from the EIS raster mid-time.
The pseudo-raster channel has one frame per EIS slit step, so its offsets are
measured per step (shown in the timeline), not to the mid-time.
""")

code("""
p = prods.copy()
p["t_obs"] = pd.to_datetime(p["t_obs"], format="ISO8601")
p["dt_s"] = (p["t_obs"] - t_ref).dt.total_seconds().round(1)
p["MB"] = [Path(x).stat().st_size / 1e6 if Path(x).exists() else np.nan for x in p["path"]]
p["role"] = p.get("role", pd.Series(index=p.index, dtype=object)).replace({"nan": ""}).fillna("")
summary = (p.groupby(["modality", "role"], dropna=False)
             .agg(files=("path", "size"), dt_min_s=("dt_s", "min"), dt_max_s=("dt_s", "max"), MB=("MB", "sum"))
             .round(1))
summary
""")

code("""
steps = pd.read_parquet(P.manifest / "eis_steps.parquet"); steps["t_step"] = pd.to_datetime(steps["t_step"], format="ISO8601")
fig, ax = plt.subplots(figsize=(8, 2.6))
rows = {"EIS steps": 0, "AIA 193 pseudo-raster": 1, "AIA context": 2, "HMI": 3, "XRT": 4}
ax.vlines(steps.t_step, -0.3, 0.3, color="k", lw=1)
ax.axvspan(pd.Timestamp(info["date_obs"]), pd.Timestamp(info["date_end"]), color="0.9", zorder=0)
sel = {"AIA 193 pseudo-raster": (p.modality == "AIA") & (p.role == "pseudoraster"),
       "AIA context": (p.modality == "AIA") & (p.role == "context"),
       "HMI": p.modality == "HMI", "XRT": p.modality == "XRT"}
for name, s in sel.items():
    ax.plot(p.loc[s, "t_obs"], np.full(s.sum(), rows[name]), "o", ms=4, label=name)
for r in p[p.modality == "XRT"].itertuples():
    ax.annotate(r.filter, (r.t_obs, 4), xytext=(0, 6), textcoords="offset points", fontsize=7, ha="center")
ax.axvline(t_ref, color="C3", ls="--", lw=1)
ax.set_yticks(list(rows.values()), list(rows.keys())); ax.set_ylim(-0.6, 4.8)
ax.set_title("Observation times around the EIS raster (grey) — dashed: t_ref (raster mid-time)")
fig.autofmt_xdate(); plt.tight_layout()
""")

md("""
## 2. Spatial context

AIA 193 context cutout with the EIS field of view from the header (dashed) and after
coalignment (solid), and the XRT fields of view. HMI line-of-sight field alongside.
""")

code("""
def first(mod, role=None, chan=None):
    s = p[p.modality == mod]
    if role: s = s[s.role == role]
    if chan is not None: s = s[s.channel.astype(str) == str(chan)]
    return sunpy.map.Map(s.iloc[0].path) if len(s) else None

a193, hmi = first("AIA", "context", 193), first("HMI")
b = info["box_arcsec"]; d = (shifts["EIS"]["dx"], shifts["EIS"]["dy"])
bt = [b[0] + d[0], b[1] + d[0], b[2] + d[1], b[3] + d[1]]
fig, axs = plt.subplots(1, 2, figsize=(9, 4.2))
for ax, m, cmap, title in [(axs[0], a193, "sdoaia193", "AIA 193 (DN/pix/s, log)"),
                           (axs[1], hmi, "gray", "HMI B_los (G)")]:
    if m is None: ax.set_visible(False); continue
    dat = m.data if cmap == "gray" else np.log10(np.clip(m.data, 1, None))
    vmin, vmax = (-500, 500) if cmap == "gray" else robust(dat)
    ax.imshow(dat, extent=extent(m), cmap=cmap, vmin=vmin, vmax=vmax)
    ax.plot(*box_xy(b), "w--", lw=1, label="EIS header")
    ax.plot(*box_xy(bt), "c-", lw=1.2, label="EIS coaligned")
    for r in p[p.modality == "XRT"].itertuples():
        xm = sunpy.map.Map(r.path); key = f"XRT_{r.filter}"
        e = extent(xm, (shifts.get(key, {}).get("dx", 0), shifts.get(key, {}).get("dy", 0)))
        ax.plot(*box_xy(e), "-", lw=0.8, color="orange", label=f"XRT {r.filter}")
    ax.set_title(title); ax.set_xlabel("Solar X (arcsec)"); ax.set_ylabel("Solar Y (arcsec)")
axs[0].legend(fontsize=6, loc="lower left")
plt.tight_layout()
""")

md("""
## 3. EIS spectral products

Window-integrated intensities for every spectral window (cool → hot), and the
Fe XII 192.39 Å Gaussian-fit Doppler velocity and width.
""")

code("""
e = p[p.modality == "EIS"].copy()
wins = e[e["product"].str.startswith("eis_win")].sort_values("product")
fits_ = e[e["product"].str.startswith("eis_fit")]
items = list(wins.itertuples()) + [r for r in fits_.itertuples() if r.product.endswith(("_vel", "_wid"))]
n = len(items); ncol = 6; nrow = int(np.ceil(n / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(1.6 * ncol, 3.3 * nrow), squeeze=False)
for ax, r in zip(axs.ravel(), items):
    m = sunpy.map.Map(r.path); dat = np.asarray(m.data, float)
    if r.product.endswith("_vel"):
        ax.imshow(dat, extent=extent(m), cmap="RdBu_r", vmin=-15, vmax=15, aspect="auto"); lab = "Fe XII v (km/s)"
    elif r.product.endswith("_wid"):
        ax.imshow(dat, extent=extent(m), cmap="viridis", vmin=robust(dat)[0], vmax=robust(dat)[1], aspect="auto"); lab = "Fe XII width (Å)"
    else:
        ld = np.log10(np.clip(dat, 1e-3, None)); ld[~np.isfinite(dat)] = np.nan
        ax.imshow(ld, extent=extent(m), cmap="inferno", vmin=robust(ld, 2, 99.5)[0],
                  vmax=robust(ld, 2, 99.5)[1], aspect="auto")
        lab = r.line.replace(" ", "\\u2009")
    ax.set_title(lab, fontsize=7); ax.set_xticks([]); ax.set_yticks([])
for ax in axs.ravel()[n:]: ax.set_visible(False)
plt.suptitle("EIS raster products (native 4\\"x1\\" sampling)", fontsize=9); plt.tight_layout()
""")

md("""
## 4. Coalignment diagnostics

EIS is matched to an AIA 193 **pseudo-raster**: each EIS slit position is compared with
the AIA frame nearest in time to that slit step. XRT is matched to the best-correlating
of AIA 335/94/211. The shift table is the pointing correction (true = header + shift).
""")

code("""
rows_ = []
for k, v in shifts.items():
    rows_.append(dict(instrument=k, dx_arcsec=v.get("dx"), dy_arcsec=v.get("dy"),
                      r_before=v.get("r_before"), r_after=v.get("r_after"),
                      reference=v.get("reference", v.get("note")), flagged=v.get("flagged")))
pd.DataFrame(rows_).round(3)
""")

code("""
dg = np.load(P.aligned / "coalign_diag.npz")
if "eis_target" in dg:
    eis_t = np.log10(np.clip(dg["eis_target"], 1e-3, None))
    fig, axs = plt.subplots(1, 4, figsize=(9, 3.6), gridspec_kw={"width_ratios": [1, 1, 1, 1.6]})
    for ax, img, t in [(axs[0], eis_t, "EIS Fe XII 192 (log)"),
                       (axs[1], dg["pr_before"], "AIA 193 pseudo-raster\\nheader pointing"),
                       (axs[2], dg["pr_after"], "AIA 193 pseudo-raster\\ncoaligned")]:
        ax.imshow(img, cmap="sdoaia193", aspect="auto", vmin=robust(img)[0], vmax=robust(img)[1])
        ax.set_title(t, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
    s = dg["eis_surface"]
    sc = axs[3].scatter(s[:, 0], s[:, 1], c=s[:, 2], cmap="viridis", s=18, marker="s")
    axs[3].plot(shifts["EIS"]["dx"], shifts["EIS"]["dy"], "r+", ms=12)
    axs[3].set_title("correlation vs shift (coarse 4\\" grid)", fontsize=8)
    axs[3].set_xlabel("dx (arcsec)"); axs[3].set_ylabel("dy (arcsec)"); plt.colorbar(sc, ax=axs[3], shrink=0.8)
    plt.tight_layout()
""")

md("""
## 5. All modalities on the common grid

`stack.fits`: one layer per channel on a common helioprojective grid covering the
coaligned EIS field of view, sampled for the raster mid-time (solar rotation applied
to each map from its own observation time). Titles show Δt = t_obs − t_ref.
""")

code("""
hd = fits.open(P.aligned / "stack.fits")
keep_eis = ("192P410", "263P300", "270P510", "256P950")   # Fe XII, Fe XXIII, Fe XIV, S XIII
layers = [h for h in hd[1:] if not h.name.upper().startswith("EIS_WIN") or h.name.upper().endswith(keep_eis)]
pref = ["AIA_94", "AIA_131", "AIA_171", "AIA_193", "AIA_211", "AIA_335", "AIA_304", "AIA_193_PSEUDORASTER",
        "HMI_BLOS"]
layers.sort(key=lambda h: (pref.index(h.name.upper()) if h.name.upper() in pref else 50 + ("XRT" not in h.name.upper()), h.name))
g0 = hd[1].header
ext = [g0["CRVAL1"] - g0["CDELT1"] / 2, g0["CRVAL1"] + g0["CDELT1"] * (g0["NAXIS1"] - 0.5),
       g0["CRVAL2"] - g0["CDELT2"] / 2, g0["CRVAL2"] + g0["CDELT2"] * (g0["NAXIS2"] - 0.5)]
cmaps = {"AIA_94": "sdoaia94", "AIA_131": "sdoaia131", "AIA_171": "sdoaia171", "AIA_193": "sdoaia193",
         "AIA_211": "sdoaia211", "AIA_335": "sdoaia335", "AIA_304": "sdoaia304", "AIA_193_PSEUDORASTER": "sdoaia193"}
ncol = 7; nrow = int(np.ceil(len(layers) / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(1.45 * ncol, 2.25 * nrow), squeeze=False)
for ax, h in zip(axs.ravel(), layers):
    a = h.data.astype(float); name = h.name.upper()
    if name.startswith("HMI"):
        ax.imshow(a, extent=ext, cmap="gray", vmin=-800, vmax=800)
    elif name.endswith("_VEL"):
        ax.imshow(a, extent=ext, cmap="RdBu_r", vmin=-15, vmax=15)
    elif name.endswith("_WID"):
        ax.imshow(a, extent=ext, cmap="viridis", vmin=robust(a)[0], vmax=robust(a)[1])
    else:
        la = np.log10(np.clip(a, 1e-3, None)); la[~np.isfinite(a)] = np.nan
        cm = cmaps.get(name, "hinodexrt" if name.startswith("XRT") else "inferno")
        ax.imshow(la, extent=ext, cmap=cm, vmin=robust(la)[0], vmax=robust(la)[1])
    ax.set_title(f"{name.lower()}\\nΔt={h.header['DT_S']:+.0f}s", fontsize=6)
    ax.set_xticks([]); ax.set_yticks([])
for ax in axs.ravel()[len(layers):]: ax.set_visible(False)
plt.suptitle(f"{len(hd) - 1} co-registered channels, {g0['NAXIS1']}x{g0['NAXIS2']} px at {g0['CDELT1']:.1f}\\"", fontsize=9)
plt.tight_layout()
""")

md("""
## 6. Cross-match queries

Two kinds of query, both in **Carrington heliographic coordinates** so that they carry
over unchanged to a HATS catalogue (Phase 5), where `lon_carr`/`lat_carr` play the
role of RA/Dec.

**6a — cone query on the common-grid table** (`pixels.parquet`): everything within a
heliographic radius of the brightest XRT pixel, all modalities at once.
""")

code("""
px = pd.read_parquet(P.aligned / "pixels.parquet")
xrt_cols = [c for c in px.columns if c.startswith("xrt_") and not c.startswith("dt_")]

def hav_deg(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))

key = xrt_cols[0] if xrt_cols else "aia_94"
peak = px.loc[px[key].idxmax()]
radius_deg = 1.0     # ~12 Mm on the solar surface
cone = px[hav_deg(px.lon_carr, px.lat_carr, peak.lon_carr, peak.lat_carr) <= radius_deg]
show = ["lon_carr", "lat_carr", "mu", key, "aia_94", "aia_335", "aia_193_pseudoraster", "hmi_blos"] + \\
       [c for c in px.columns if c.startswith("eis_win07") or c.endswith("_vel")]
show = [c for c in show if c in cone and not c.startswith("dt_")]
print(f"cone of {radius_deg}° around brightest {key} pixel at (lon, lat) = "
      f"({peak.lon_carr:.2f}°, {peak.lat_carr:.2f}°): {len(cone)} grid pixels")
cone[show].describe().loc[["mean", "std", "min", "max"]].round(2)
""")

md("""
**6b — native-resolution cross-match between two modalities** (`native.parquet`).
Each EIS pixel (4″×1″) is paired with every XRT pixel whose centre lies within
2″ on the sky, done as a nearest-neighbour search on the unit sphere of Carrington
coordinates. This is the operation LSDB/HATS performs at scale.
""")

code("""
from scipy.spatial import cKDTree
nat = pd.read_parquet(P.aligned / "native.parquet")

def unit(lon, lat):
    lon, lat = np.radians(lon), np.radians(lat)
    return np.c_[np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]

eis_ch = sorted(c for c in nat.channel.unique() if c.startswith("eis_win07")) or \\
         sorted(c for c in nat.channel.unique() if c.startswith("eis_win"))
L = nat[nat.channel == eis_ch[0]].dropna(subset=["lon_carr"])
xch = sorted(c for c in nat.channel.unique() if c.startswith("xrt_"))
Rch = xch[0] if xch else "aia_94"
R = nat[nat.channel == Rch].dropna(subset=["lon_carr"])
# 2 arcsec on the sky ~ 1450 km on the Sun ~ 0.12 heliographic degrees near disk centre
rad = 2 * 725.0 / 696000.0   # radians on the unit sphere (chord ~ angle for small angles)
tree = cKDTree(unit(R.lon_carr, R.lat_carr))
pairs = tree.query_ball_point(unit(L.lon_carr, L.lat_carr), r=rad)
xm = pd.DataFrame({"eis": L["value"].values,
                   "xrt_mean": [R["value"].values[i].mean() if len(i) else np.nan for i in pairs],
                   "n_matched": [len(i) for i in pairs]}).dropna()
print(f"{eis_ch[0]} x {Rch}: {len(xm)} of {len(L)} EIS pixels matched, median {xm.n_matched.median():.0f} XRT px each")
fig, ax = plt.subplots(figsize=(4, 3.3))
ax.loglog(xm.xrt_mean.clip(1e-3), xm.eis.clip(1e-3), ".", ms=3, alpha=0.6)
ax.set_xlabel(f"{Rch} (DN/s, mean of matches)"); ax.set_ylabel(f"{eis_ch[0]}")
ax.set_title("EIS hot line vs XRT, native-resolution cross-match", fontsize=8); plt.tight_layout()
""")

md("""
**6c — the same cross-match with LSDB** (optional; the HATS engine planned for Phase 5).
LSDB expects sky coordinates, so Carrington longitude/latitude are passed as RA/Dec; the
matching radius is then an angle on the *solar* sphere (0.12° ≈ 1.45 Mm ≈ 2″ near disk centre).
""")

code("""
try:
    import lsdb
    radius = 0.12 * 3600   # arcsec of *heliographic* angle
    left = lsdb.from_dataframe(L[["lon_carr", "lat_carr", "value"]].rename(columns={"value": "intensity"}),
                               ra_column="lon_carr", dec_column="lat_carr")
    right = lsdb.from_dataframe(R[["lon_carr", "lat_carr", "value"]].rename(columns={"value": "intensity"}),
                                ra_column="lon_carr", dec_column="lat_carr",
                                margin_threshold=2 * radius)   # margin must exceed the match radius
    res = left.crossmatch(right, radius_arcsec=radius, n_neighbors=1, suffixes=("_eis", "_xrt")).compute()
    print(f"LSDB crossmatch: {len(res)} EIS pixels with a nearest XRT pixel within 0.12°")
    display(res.head())
except Exception as e:  # lsdb missing or API change: the KD-tree version above is the reference
    print(f"LSDB step skipped: {type(e).__name__}: {e}")
""")

md("## 7. Disk usage of this run")

code("""
rows_ = []
for sub in ("raw", "prep", "aligned", "manifest", "cache"):
    d = P.root / sub
    if d.exists():
        rows_.append(dict(area=sub, MB=round(sum(f.stat().st_size for f in d.rglob("*") if f.is_file()) / 1e6, 1)))
du = pd.DataFrame(rows_); du.loc[len(du)] = ["TOTAL", du.MB.sum()]; du
""")

nb["cells"] = C
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = ROOT / "notebooks" / "01_summary.ipynb"
out.parent.mkdir(exist_ok=True)
nbf.write(nb, out)
print(f"wrote {out}")
