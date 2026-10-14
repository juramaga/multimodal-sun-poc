"""Configuration loading and the on-disk layout of the project."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path else ROOT / "config" / "pipeline.yaml"
    if not path.is_absolute() and not path.exists():
        path = ROOT / path          # relative to the project root (e.g. from notebooks/)
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_root"] = str(ROOT)
    return cfg


def load_epochs(cfg: dict) -> list[dict]:
    with open(ROOT / cfg["epochs_file"]) as f:
        epochs = yaml.safe_load(f)["epochs"]
    return epochs


def active_epochs(cfg: dict, override: list[str] | None = None) -> list[dict]:
    """Resolve pipeline.yaml:active_epochs (or a CLI override) to epoch dicts."""
    epochs = load_epochs(cfg)
    wanted = override or cfg.get("active_epochs", ["anchor"])
    if wanted == ["all"] or wanted == "all":
        return epochs
    out = []
    for w in wanted:
        if w == "anchor":
            out += [e for e in epochs if "anchor" in e.get("tags", [])]
        else:
            out += [e for e in epochs if e["id"] == w or w in e.get("tags", [])]
    seen, uniq = set(), []
    for e in out:
        if e["id"] not in seen:
            uniq.append(e)
            seen.add(e["id"])
    if not uniq:
        raise SystemExit(f"No epochs match {wanted}; check config/epochs.yaml")
    return uniq


@dataclass
class Paths:
    """All locations for one epoch. Nothing outside data_dir is ever written."""

    root: Path
    epoch: str

    @property
    def data(self) -> Path:
        return self.root

    def _d(self, *parts) -> Path:
        p = self.root.joinpath(*parts)
        p.mkdir(parents=True, exist_ok=True)
        return p

    # raw downloads, one folder per modality
    def raw(self, modality: str) -> Path:
        return self._d("raw", modality.lower(), self.epoch)

    # metadata products
    @property
    def manifest(self) -> Path:
        return self._d("manifest", self.epoch)

    # calibrated, native-resolution, per-modality maps
    def prep(self, modality: str) -> Path:
        return self._d("prep", modality.lower(), self.epoch)

    # coaligned common-grid products + tables
    @property
    def aligned(self) -> Path:
        return self._d("aligned", self.epoch)

    @property
    def cache(self) -> Path:
        return self._d("cache")


def paths_for(cfg: dict, epoch_id: str) -> Paths:
    root = Path(cfg["data_dir"])
    if not root.is_absolute():
        root = Path(cfg["_root"]) / root
    return Paths(root=root, epoch=epoch_id)


def jsoc_email(cfg: dict) -> str | None:
    return os.environ.get(cfg["jsoc"].get("email_env", "JSOC_EMAIL"))
