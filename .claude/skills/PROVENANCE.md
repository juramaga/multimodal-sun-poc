# Skill provenance

| Skill | Origin | License |
|---|---|---|
| `python-library-docs-first/` | Vendored unchanged from [skill-commons/curated-research-skills](https://github.com/skill-commons/curated-research-skills/tree/425cb21b31ee8c648e978ca345e2c19675f3a642/skills/python-library-docs-first) (v2.0.1, commit `425cb21`), discovered via the [Skill Commons](https://github.com/skill-commons/skill-commons) catalog. Used here because sunpy, eispac, drms, aiapy and lsdb change APIs between minor versions. | MIT (see its `LICENSE`) |
| `mmsun-pipeline/` | Written for this project | BSD-3-Clause |
| `solar-archive-access/` | Written for this project | BSD-3-Clause |
| `solar-coalignment/` | Written for this project | BSD-3-Clause |
| `epoch-scout/` | Written for this project | BSD-3-Clause |

Other Skill Commons entries were reviewed and not vendored: the astronomy skills target
stellar/galactic catalogues (Gaia, RAVE, StarHorse, MUSE), and `coseecat` covers
solar energetic electron events rather than imaging/spectroscopy archives.
`tap-pyvo-adql-access` becomes relevant only if a TAP service for the solar archives is
added later. To refresh the vendored skill:

```bash
git clone --depth 1 https://github.com/skill-commons/curated-research-skills /tmp/crs
rsync -a --delete /tmp/crs/skills/python-library-docs-first/ .claude/skills/python-library-docs-first/
```
