#!/usr/bin/env python
"""Export a prepared EIS/XRT epoch as two local HATS collections."""

from mmsun.cli import die, parser, setup
from mmsun.xmatch import build_catalogs

p = parser(__doc__ or "Export EIS/XRT HATS collections")
p.add_argument("--fixture", action="store_true", help="label synthetic XRT explicitly")
args = p.parse_args()
try:
    cfg, epochs = setup(args)
    for epoch in epochs:
        left, right = build_catalogs(cfg, epoch, fixture=args.fixture)
        print(f"EIS HATS: {left}\nXRT HATS: {right}")
except (ValueError, FileNotFoundError, FileExistsError, KeyError) as exc:
    die(str(exc))
