#!/usr/bin/env python
"""Preflight: package versions, JSOC email, archive reachability, free disk space.

Exit code 0 = ready, 2 = something needs attention (details printed).
"""

import importlib
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mmsun.cli import parser, setup  # noqa: E402
from mmsun.config import jsoc_email  # noqa: E402
from mmsun.netutil import reachable  # noqa: E402

PKGS = ["numpy", "scipy", "astropy", "sunpy", "aiapy", "eispac", "xrtpy", "drms", "reproject",
        "h5py", "pandas", "pyarrow", "yaml", "matplotlib", "lsdb", "hats"]

ENDPOINTS = {
    # one record by its T_REC: an open-ended record set makes JSOC scan the series (>15 s, false alarm)
    "JSOC (AIA/HMI)": "http://jsoc.stanford.edu/cgi-bin/ajax/show_info?ds=hmi.M_45s[2021.03.06_06:48:00_TAI]&key=T_REC",
    "NRL EIS archive": "https://eis.nrl.navy.mil/level1/hdf5/2021/03/06/",
    "MSSL EIS mirror": "https://vsolar.mssl.ucl.ac.uk/eispac/hdf5/2021/03/06/",
    "SDAC XRT Level-1": "https://umbra.nascom.nasa.gov/hinode/xrt/level1/2021/03/06/H0600/",
    "AIA pointing table (JSOC)": "http://jsoc.stanford.edu/",
}


def main():
    args = parser(__doc__).parse_args()
    cfg, epochs = setup(args)
    problems = []

    print(f"python {sys.version.split()[0]}")
    if sys.version_info[:2] >= (3, 13):
        problems.append("Python >= 3.13: eispac's dependency 'pims' does not build; use the 3.12 env")
    for name in PKGS:
        try:
            mod = importlib.import_module(name)
            print(f"  {name:10s} {getattr(mod, '__version__', '?')}")
        except Exception as e:  # noqa: BLE001
            optional = name in ("lsdb", "hats")
            print(f"  {name:10s} MISSING{' (optional)' if optional else ''}: {type(e).__name__}")
            if not optional:
                problems.append(f"missing package {name}")

    email = jsoc_email(cfg)
    print(f"\nJSOC email (${cfg['jsoc']['email_env']}): {'set' if email else 'NOT SET'}")
    if not email:
        problems.append(f"export {cfg['jsoc']['email_env']}=<email registered at "
                        "http://jsoc.stanford.edu/ajax/register_email.html>")

    print("\nArchive reachability:")
    for k, url in ENDPOINTS.items():
        ok, code = reachable(url)
        print(f"  {'ok ' if ok else 'XX '} {k:28s} {code}")
        if not ok and "MSSL" not in k:
            problems.append(f"cannot reach {k} ({url})")

    free = shutil.disk_usage(Path(cfg["_root"])).free / 1e9
    print(f"\nFree disk at project root: {free:.1f} GB")
    if free < 2:
        problems.append("less than 2 GB free")
    print(f"Active epochs: {[e['id'] for e in epochs]}")

    if problems:
        print("\nNEEDS ATTENTION:\n  - " + "\n  - ".join(problems))
        sys.exit(2)
    print("\nREADY")


if __name__ == "__main__":
    os.environ.setdefault("PYTHONWARNINGS", "ignore")
    main()
