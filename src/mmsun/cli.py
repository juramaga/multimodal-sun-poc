"""Shared command-line plumbing for scripts/*.py."""

from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from pathlib import Path

import pandas as pd

from .config import active_epochs, load_config, paths_for


def parser(desc: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=desc)
    p.add_argument("--config", default=None, help="path to pipeline.yaml")
    p.add_argument("--epochs", nargs="*", default=None,
                   help="epoch ids/tags (default: pipeline.yaml:active_epochs)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def setup(args):
    # sunpy installs its own Logger subclass; it must be imported before anything
    # touches logging.getLogger("sunpy"), otherwise sunpy fails to initialise.
    import astropy  # noqa: F401
    import sunpy  # noqa: F401

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    if not args.verbose:
        warnings.filterwarnings("ignore")
        for name in ("sunpy", "astropy", "drms", "parfive", "aiapy", "eispac"):
            logging.getLogger(name).setLevel(logging.WARNING)
    cfg = load_config(args.config)
    return cfg, active_epochs(cfg, args.epochs)


def save_json(obj, path: Path):
    def conv(o):
        if isinstance(o, (pd.Timestamp,)):
            return o.isoformat()
        if hasattr(o, "item"):
            return o.item()
        if isinstance(o, Path):
            return str(o)
        return str(o)

    path.write_text(json.dumps(obj, indent=2, default=conv))


def load_json(path: Path):
    return json.loads(Path(path).read_text())


def write_parquet(df: pd.DataFrame, path: Path):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str)
    df.to_parquet(path, index=False)
    return path


def die(msg: str, code: int = 1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


__all__ = ["parser", "setup", "save_json", "load_json", "write_parquet", "die", "paths_for"]
