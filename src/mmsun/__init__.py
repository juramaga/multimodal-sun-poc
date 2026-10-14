"""Multimodal Sun proof-of-concept: EIS-anchored co-temporal SDO + Hinode samples.

Phases implemented here (see CLAUDE.md):
  0  anchor epoch selection          (config/epochs.yaml, scripts/01_discover.py)
  1  metadata-only discovery         (scripts/01_discover.py)
  2  time/space matching             (scripts/02_match.py)
  3  retrieval + calibration         (scripts/03_download.py, scripts/04_prepare.py)
  4  coalignment + common grid       (scripts/05_coalign.py)
"""

__version__ = "0.1.0"
