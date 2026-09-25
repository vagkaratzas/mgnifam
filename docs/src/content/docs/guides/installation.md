---
title: Installation
description: Install mgnifam from PyPI, bioconda, or a locked checkout.
---

```bash
pip install mgnifam        # or: uv tool install mgnifam
```

Or from [bioconda](https://bioconda.github.io/recipes/mgnifam/README.html), into its own
environment:

```bash
conda create -n mgnifam -c conda-forge -c bioconda mgnifam
conda activate mgnifam
```

Channel order matters — put `conda-forge` before `bioconda`, per the
[bioconda setup](https://bioconda.github.io/#usage). The package is `noarch`, and conda
pulls in a Python 3.13 interpreter itself, so it does not have to be the one already on
your `PATH`. `mamba`/`micromamba` work the same way with the same flags.

Requires Python >= 3.13. Verify with `mgnifam --version`.

To work on the package itself, or to reproduce published results byte-for-byte, install
from the repository against the committed lockfile instead — see
[Reproducibility](/mgnifam/reference/reproducibility/), which is scoped to that resolved
dependency set:

```bash
git clone https://github.com/vagkaratzas/mgnifam && cd mgnifam
uv sync --frozen
```

Commands in these docs are written `uv run mgnifam ...` for the cloned checkout. On a
`pip` or `conda` install, drop the `uv run` prefix.

`mgnifam --help` lists the subcommands, and `python -m mgnifam` is equivalent to the
console script.
