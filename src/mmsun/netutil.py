"""Small, dependency-light HTTP helpers with retries and partial-file safety."""

from __future__ import annotations

import time
from pathlib import Path

import requests

UA = {"User-Agent": "multimodal-sun-poc/0.1 (+research prototype)"}


def get(url: str, *, headers: dict | None = None, retries: int = 4, timeout: int = 60,
        stream: bool = False) -> requests.Response:
    h = dict(UA, **(headers or {}))
    delay = 1.0
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=h, timeout=timeout, stream=stream)
            if r.status_code == 404:
                r.raise_for_status()
            if r.status_code in (408, 429) or r.status_code >= 500:
                raise requests.HTTPError(f"{r.status_code} for {url}")
            r.raise_for_status()
            return r
        except requests.HTTPError as e:
            if "404" in str(e) or attempt == retries - 1:
                raise
        except requests.RequestException:
            if attempt == retries - 1:
                raise
        time.sleep(delay)
        delay *= 2
    raise RuntimeError("unreachable")


def download(url: str, target: Path, *, chunk: int = 1 << 20) -> Path:
    """Stream to ``target.part`` then rename, so interrupted runs never leave junk."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        return target
    tmp = target.with_suffix(target.suffix + ".part")
    with get(url, stream=True, timeout=120) as r, open(tmp, "wb") as f:
        for block in r.iter_content(chunk):
            f.write(block)
    tmp.rename(target)
    return target


def reachable(url: str, timeout: int = 15) -> tuple[bool, str]:
    try:
        r = requests.get(url, headers=UA, timeout=timeout, stream=True)
        r.close()
        return r.status_code < 400, str(r.status_code)
    except requests.RequestException as e:
        return False, type(e).__name__
