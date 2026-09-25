<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/vagkaratzas/mgnifam/main/assets/logo-dark.png">
    <img src="https://raw.githubusercontent.com/vagkaratzas/mgnifam/main/assets/logo.png" width="600" alt="mgnifam">
  </picture>
</p>

# mgnifam

[![PyPI](https://img.shields.io/pypi/v/mgnifam)](https://pypi.org/project/mgnifam/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/mgnifam)](https://bioconda.github.io/recipes/mgnifam/README.html#package-package%20&#x27;mgnifam&#x27;)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22879938.svg)](https://doi.org/10.5281/zenodo.22879938)

Iterative HMM-based protein family generation over very large sequence databases.

Given a chunk of MMseqs2 clusters and a protein FASTA, `mgnifam generate_families` builds an
HMM from each cluster, recruits new members from the whole database, re-aligns, and
either converges on a family or discards the cluster. It is the core algorithm of the
[`mgnifams`](https://github.com/vagkaratzas/mgnifams) Nextflow pipeline, extracted into
a standalone, tested package.

**Documentation: [vagkaratzas.github.io/mgnifam](https://vagkaratzas.github.io/mgnifam/)**

## Install

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
[Reproducibility](https://vagkaratzas.github.io/mgnifam/reference/reproducibility/), which is scoped to that resolved dependency set:

```bash
git clone https://github.com/vagkaratzas/mgnifam && cd mgnifam
uv sync --frozen
```

## Quick start

Build families from a chunk of clusters:

```bash
mgnifam generate_families \
    --clusters_chunk clusters.tsv \
    --fasta_file mgnifams_input.fa
```

Refresh existing families against a new release:

```bash
mgnifam update_families \
    --hmm_input previous_output/hmm \
    --fasta_file new_release.fa
```

`--fasta_file` must be an **uncompressed** FASTA with unique sequence names. Pass every
threshold explicitly on a production run.

The documentation covers the rest:

- [Generating families](https://vagkaratzas.github.io/mgnifam/guides/generate-families/):
  input formats, sequence naming, every flag, and sharing one SSI index across chunks
- [Updating families](https://vagkaratzas.github.io/mgnifam/guides/update-families/):
  `--skip_refine`, family identity, and retrying a chunk in place
- [Outputs](https://vagkaratzas.github.io/mgnifam/reference/outputs/): every file, exit
  status, and the MultiQC summary
- [Reproducibility](https://vagkaratzas.github.io/mgnifam/reference/reproducibility/):
  what is byte-identical, and how outputs differ from the legacy script
- [Python API](https://vagkaratzas.github.io/mgnifam/reference/python-api/)

## Development

```bash
uv lock --check && uv sync --frozen
uv run pre-commit install
uv run pre-commit run --all-files
uv run pytest
```

The documentation site lives in `docs/`; see
[Development](https://vagkaratzas.github.io/mgnifam/reference/development/) to build it.
