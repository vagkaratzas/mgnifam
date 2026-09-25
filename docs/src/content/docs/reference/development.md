---
title: Development
description: Set up a checkout, run the checks, and build this site.
---

```bash
uv lock --check && uv sync --frozen
uv run pre-commit install
uv run pre-commit run --all-files
uv run pytest
```

Read [`AGENTS.md`](https://github.com/vagkaratzas/mgnifam/blob/main/AGENTS.md) before
changing `generate_families.py` or `update_families.py`: it lists the invariants each
module depends on and the test that guards each one.

## This site

The documentation lives in `docs/`, built with [Starlight](https://starlight.astro.build/)
and [starlight-pydocs](https://ewels.github.io/starlight-pydocs/). The
[Python API](/mgnifam/reference/python-api/) pages are generated from the docstrings at
build time by [Griffe](https://mkdocstrings.github.io/griffe/), which runs through `uvx`
and reads the source without importing it.

```bash
cd docs
npm ci
npm run dev      # http://localhost:4321/mgnifam/
npm run build    # static site in docs/dist/
```

Requires Node >= 22.12 and `uv` on `PATH`. Every push to `main` deploys the site to
GitHub Pages.
