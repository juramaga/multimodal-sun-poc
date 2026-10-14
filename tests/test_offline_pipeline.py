"""Offline end-to-end test (no network): bundled EIS raster + synthetic AIA/HMI/XRT.

Runs 01_discover (EIS part only) -> fixture -> 04_prepare -> 05_coalign and checks
that the injected pointing errors are recovered. Run with:
    pytest -q tests/test_offline_pipeline.py      (or)    python tests/test_offline_pipeline.py
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from make_fixture import EIS_SHIFT, XRT_SHIFT  # noqa: E402

DATA = Path(os.environ.get("MMSUN_TEST_DATA", ROOT / "data_test"))


def run(*cmd):
    print("$", " ".join(map(str, cmd)))
    r = subprocess.run([sys.executable, *map(str, cmd)], cwd=ROOT, capture_output=True, text=True)
    print(r.stdout[-3000:])
    if r.returncode != 0:
        print(r.stderr[-5000:])
        raise RuntimeError(f"{cmd[0]} failed")
    return r.stdout


def test_offline_pipeline():
    cfg = yaml.safe_load((ROOT / "config/pipeline.yaml").read_text())
    cfg["data_dir"] = str(DATA)
    cfg["xrt"]["header_cache"] = str(DATA / "cache/xrt_headers")
    DATA.mkdir(parents=True, exist_ok=True)
    cfgp = DATA / "pipeline_test.yaml"
    cfgp.write_text(yaml.safe_dump(cfg))
    c = ["--config", cfgp, "--epochs", "anchor"]

    run("scripts/01_discover.py", *c, "--skip", "aia", "hmi", "xrt")
    run("tests/make_fixture.py", DATA)

    # Phase 2 + dry-run on synthetic candidate tables
    run("scripts/02_match.py", *c)
    import pandas as pd
    m = pd.read_parquet(DATA / "manifest/eis_20210306_064444/match.parquet")
    pr = m[(m.modality == "AIA") & (m.role == "pseudoraster")]
    ctx = m[(m.modality == "AIA") & (m.role == "context")]
    xr = m[m.modality == "XRT"]
    assert len(ctx) == 7 and 20 <= len(pr) <= 27, (len(ctx), len(pr))
    assert (pr.dt_s.astype(float).abs() <= 12).all()
    assert list(xr.channel) == ["Al_poly", "Be_thin"] and set(xr.role) == {"ar"}, xr[["channel", "role"]]
    assert (m.modality == "HMI").sum() == 1
    run("scripts/03_download.py", *c, "--dry-run")

    run("scripts/04_prepare.py", *c, "--no-aiapy-tables")
    run("scripts/05_coalign.py", *c)

    sh = json.loads((DATA / "aligned/eis_20210306_064444/shifts.json").read_text())
    assert abs(sh["EIS"]["dx"] - EIS_SHIFT[0]) < 1.5 and abs(sh["EIS"]["dy"] - EIS_SHIFT[1]) < 1.5, sh["EIS"]
    for k in ("XRT_Al_poly", "XRT_Be_thin"):
        assert abs(sh[k]["dx"] - XRT_SHIFT[0]) < 1.5 and abs(sh[k]["dy"] - XRT_SHIFT[1]) < 1.5, sh[k]
    import pandas as pd
    px = pd.read_parquet(DATA / "aligned/eis_20210306_064444/pixels.parquet")
    for col in ("aia_193", "aia_193_pseudoraster", "hmi_blos", "xrt_Al_poly", "lon_carr"):
        assert col in px and px[col].notna().mean() > 0.9, col
    assert any(col.startswith("eis_fit_fe_12_192_394_vel") for col in px)
    print("OFFLINE TEST PASSED", json.dumps({k: (round(v.get("dx", 0), 2), round(v.get("dy", 0), 2))
                                             for k, v in sh.items()}))


if __name__ == "__main__":
    test_offline_pipeline()
