<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/vagkaratzas/mgnifam/main/assets/logo-dark.png">
    <img src="https://raw.githubusercontent.com/vagkaratzas/mgnifam/main/assets/logo.png" width="600" alt="mgnifam">
  </picture>
</p>

# mgnifam

[![PyPI](https://img.shields.io/pypi/v/mgnifam)](https://pypi.org/project/mgnifam/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/mgnifam)](https://bioconda.github.io/recipes/mgnifam/README.html#package-package%20&#x27;mgnifam&#x27;)

Iterative HMM-based protein family generation over very large sequence databases.

Given a chunk of MMseqs2 clusters and a protein FASTA, `mgnifam generate_families` builds an
HMM from each cluster, recruits new members from the whole database, re-aligns, and
either converges on a family or discards the cluster. It is the core algorithm of the
[`mgnifams`](https://github.com/vagkaratzas/mgnifams) Nextflow pipeline, extracted into
a standalone, tested package.

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
[Reproducibility](#reproducibility), which is scoped to that resolved dependency set:

```bash
git clone https://github.com/vagkaratzas/mgnifam && cd mgnifam
uv sync --frozen
```

## Usage

Commands below are written `uv run mgnifam ...` for the cloned checkout. On a `pip` or
`conda` install, drop the `uv run` prefix.

```bash
uv run mgnifam generate_families \
    --clusters_chunk clusters.tsv \
    --fasta_file mgnifams_input.fa
```

Only those two are required. Every other flag defaults to the value below, so the run
above is equivalent to spelling all of them out:

```bash
uv run mgnifam generate_families \
    --clusters_chunk clusters.tsv \
    --fasta_file mgnifams_input.fa \
    --output_dir output \
    --cpus 8 \
    --chunk_num 1 \
    --discard_min_rep_length 75 \
    --discard_max_rep_length 2000 \
    --discard_min_starting_membership 0.9 \
    --max_seq_identity 0.8 \
    --max_seed_seqs 2000 \
    --max_gap_occupancy 0.5 \
    --recruit_evalue_cutoff 0.001 \
    --recruit_hit_length_percentage 0.9
```

`--clusters_chunk` is a headerless TSV of `representative<TAB>member`.
`--fasta_file` must be an **uncompressed** FASTA — Easel cannot seek within a gzip
stream — and its sequence names must be unique.

### Optional flags

Pass every threshold explicitly on a production run. The defaults exist for ad-hoc use;
relying on them means a forgotten flag produces a plausible-looking family set instead
of an error.

| flag | default | meaning |
|---|---|---|
| `--cpus` | `8` | Threads for FAMSA, `hmmsearch` and `hmmalign`. |
| `--chunk_num` | `1` | Prefix for every output file and directory. Must match `[A-Za-z0-9._-]+`. |
| `--discard_min_rep_length` | `75` | Discard a cluster whose representative is shorter than this. |
| `--discard_max_rep_length` | `2000` | Discard a cluster whose representative is longer than this. |
| `--discard_min_starting_membership` | `0.9` | Discard a family if fewer than this fraction of the original cluster members are still recruited by the final model. |
| `--max_seq_identity` | `0.8` | Redundancy cutoff when trimming a full MSA down to the next seed. |
| `--max_seed_seqs` | `2000` | Cap on sequences kept in a seed MSA. |
| `--max_gap_occupancy` | `0.5` | Trim columns off both **ends** of the seed MSA until one clears this occupancy. Interior columns are kept. |
| `--recruit_evalue_cutoff` | `0.001` | `hmmsearch` E-value threshold for recruiting new members. |
| `--recruit_hit_length_percentage` | `0.9` | Minimum hit length as a fraction of the model length. |
| `--fasta_index` | `<output_dir>/<fasta basename>.ssi` | Path to an Easel SSI index. Used exactly as given and never rebuilt; only the default path is built automatically. |
| `--output_dir` | `output` | Root directory for every generated file and folder. |
| `--batch_size` | `2 * cpus` | How many families are searched per `hmmsearch` wave. Keep it `>= cpus`. |
| `--prefetch_targets` | off | Load the database into RAM once instead of streaming it per query. Faster, `O(database)` memory, **identical results**. |

Streaming re-reads and re-parses the database once per query. `--prefetch_targets`
parses it once and keeps it in RAM; the results are byte-identical either way, so the
flag is purely a memory-vs-time dial. Leave it off unless the database fits comfortably
in RAM.

**On a production run, build the index once and share it.** Every chunk task would
otherwise re-index the whole database:

```bash
# once, upstream -- either of these
uv run python -c "from mgnifam.generate_families import build_ssi_index; \
                  build_ssi_index('db.fa', 'db.fa.ssi')"
esl-sfetch --index db.fa                     # HMMER/Easel, e.g. the nf-core module
# then, per chunk
uv run mgnifam generate_families --fasta_index db.fa.ssi ...
```

A supplied index is used as given and never rebuilt, so parallel chunk tasks can share
one read-only index safely — including one staged as a symlink by a workflow manager.
It is an error for it to be missing rather than a request to build one there, and one
that does not match the FASTA fails at the first fetch instead of being silently
replaced. The FASTA's filename need not match the one it was indexed under.

An index from `esl-sfetch --index` is interchangeable with one from `build_ssi_index`
for whole-record fetches, which is all `generate_families` performs. The two are not
byte-identical: `esl-sfetch` also records each record's `data_offset` and
`record_length`, which enables `esl-sfetch -c <from>..<to>` subsequence fetches against
its own index but not against ours, and it sizes the index's filename field from the
path you typed, so its output is not reproducible across directories. Ours is.

## Updating existing families

`mgnifam update_families` refreshes families that already exist as HMMs against a new
database. It is the answer to "a new release came out" — you do not re-derive the
families from their original clusters, you search the models you already have.

```bash
uv run mgnifam update_families \
    --hmm_input previous_output/hmm \
    --fasta_file new_release.fa
```

`--hmm_input` is either a directory of `.hmm`/`.hmm.gz` files or a single multi-model
library (`hmm.lib.gz` works). Both forms produce identical output for the same models.
`--fasta_file` is uncompressed, for the same Easel reason as above.

Every threshold flag from `generate_families` carries over with the same name and default.

| flag | meaning |
|---|---|
| `--skip_refine` | Recruit once and align. The model, seed MSA and RF line are unchanged, so only `hmm/` and `full_msa/` are written and `--max_seq_identity`, `--max_seed_seqs` and `--max_gap_occupancy` are inert. Without it, the full three-round refine loop runs and writes the complete artifact set. |
| `--chunk_num` | Labels the per-chunk aggregate files **only**. Family names come from the models, so nothing is renumbered. |

### What identity means here

A family keeps the name its model carries in its `NAME` field — `1_7` stays `1_7` across
releases, in the filenames *and* in every identity-bearing field inside the outputs. Two
consequences:

- `<chunk>_updated_metadata.csv`'s `family_id` column holds `1_7`, not a bare integer. That
  differs from `generate_families`, whose ids are a rank.
- Chunks sharing one output root **must own disjoint family names**. Nothing enforces it,
  because the names come from the input models rather than from `--chunk_num`.

A `NAME` must match `[A-Za-z0-9._-]+` and be neither `.` nor `..`. It is interpolated into
artifact paths and into CSV fields, and it arrives from a file this tool did not write.

### Outputs

Per-family artifacts land in the same `hmm/`, `full_msa/`, `seed_msa/` and `rf/`
directories, named by family. Aggregates are `<chunk>_updated_*`: `families.tsv`,
`metadata.csv`, `discarded.csv`, `successful.txt`, `converged.txt`, `reps.fasta.gz`,
`delta.csv`, and `<chunk>_updated.log`.

`<chunk>_updated_delta.csv` is what an update run is *for* — one row per family, whether it
survived or not:

```
family_id,model_length_before,model_length_after,round1_recruits,full_msa_size,retention,rounds_run,converged,outcome
```

Every field but `family_id`, `model_length_before` and `outcome` may be empty, because a
family discarded early never reached the stage that would produce one. `model_length_after`
is the length of the model that recruited the final membership. `retention` is the fraction
of round 1's own recruits still present at the end — under `--skip_refine` that is 1.0 by
construction, since there are no later rounds to drift.

`outcome` is `successful` or the discard reason. `no hits in the new database` means the
model found nothing at all in the new release; `low complexity model - confounding cluster`
means it found hits and none cleared the envelope-length filter. For an update run that
distinction is the point.

Give each run its own `--output_dir`. Re-running the same models into the same directory is
allowed, so a failed chunk can be retried in place. Running a *smaller* set of models over a
directory that still holds a larger one is refused rather than silently cleaned up:
`generate_families` can clear its own past output because it derives names as
`<chunk>_<rank>`, but an updated family keeps its model's name and `--chunk_num` never
appears in a per-family filename, so nothing on disk says which run wrote `hmm/1_7.hmm.gz`.

Input files must not overlap output paths, including through symlinks or hard links.
This is checked before writing, for both model directories and single-file libraries.

### Cost

`hmmsearch` is `O(n_families x database)` and this command does not change that.
`--skip_refine` is one database pass per family; refining is up to three. Against a
billion-sequence release that term, not the alignment, is what to budget.

`mgnifam --help` lists the subcommands, and `python -m mgnifam` is equivalent to the
console script.

## Outputs

Written under `--output_dir` (default: `output`), keyed by `--chunk_num`:

One file per family, so one directory each:

| path | contents |
|---|---|
| `seed_msa/<chunk>_<id>.sto.gz` | seed alignment |
| `full_msa/<chunk>_<id>.sto.gz` | full alignment |
| `hmm/<chunk>_<id>.hmm.gz` | the family model |
| `rf/<chunk>_<id>.txt` | reference-annotation line |

One file per chunk, so flat in the output root:

| path | contents |
|---|---|
| `<chunk>_reps.fasta.gz` | one representative per family |
| `<chunk>_families.tsv` | `family_id<TAB>sequence` |
| `<chunk>_metadata.csv` | one row per family |
| `<chunk>_successful.txt` | representatives that produced a family |
| `<chunk>_discarded.csv` | one row per discarded cluster |
| `<chunk>_converged.txt` | ids of successful families that converged naturally |
| `<chunk>.log` | run log |

Family ids are a 1-based rank among *successful* families, in cluster-file order.

Both CSVs carry a header row, so they load with `pandas.read_csv` as they are:

| file | columns |
|---|---|
| `<chunk>_metadata.csv` | `family_id,full_msa_size,protein,region,length,sequence,consensus,converged` |
| `<chunk>_discarded.csv` | `representative,reason,value` |

`protein` is quoted; `region` is `<start>-<end>` on the parent protein, or `-` when the
representative spans a whole unsliced record. The representative is the highest-scoring
reported domain of HMMER's top-ranked hit. The header is written before the run starts, so
a chunk that produces no families still yields a parseable file.

### Exit status

The status describes whether the output is safe to consume, not only whether the process
stopped:

| Code | Meaning |
|---|---|
| `0` | Chunk completed. Every family landed on exactly one side of the split (discarded or successful). Output is complete and safe to consume. |
| `1` | Fatal: the run died before finishing. **Output is incomplete and must not be consumed** — re-run the chunk. This is what a dead output sink (ENOSPC, EIO) produces, because the discard re-emit cannot record its own failure. |
| `2` | Usage error from `argparse`. Nothing ran. |
| `3` | Chunk completed, but one or more families died of an internal error and were recorded as discards. Output is complete and self-consistent, but those clusters produced no family — re-run the chunk once the cause is fixed, or accept the loss. |

## Why this is fast now

The previous implementation took roughly eight months to process the full database.
Three defects accounted for most of it:

1. **`run_initial_msa` was O(database × members), per family.** A `map()` iterator was
   rebuilt inside a comprehension's condition, turning a membership test into a full
   linear scan of the cluster for every one of the billions of database sequences.
   It is now a constant-time SSI lookup per member.
2. **Cluster selection was O(N²)** — the cluster table was boolean-masked and
   re-filtered once per family. It is now a single grouping pass.
3. **The exit-branch `hmmsearch` re-ran a search that had just been performed** with
   the identical HMM, differing only in a post-filter. Its hits are now cached and
   re-filtered, saving a full database pass per family.

On top of that, the entire FASTA was held in RAM twice — once as a
`DigitalSequenceBlock` and once as a Python `dict` of `DigitalSequence` objects. In the
default (streaming) mode both are gone: targets stream from disk, and random access goes
through an Easel SSI index. Passing `--prefetch_targets` deliberately restores the first
copy, trading that memory back for speed. Families are searched in batched waves, so `hmmsearch` uses up to `--cpus` workers
whenever enough families remain in the wave.

## Reproducibility

For the dependency set resolved in the committed `uv.lock` (install with
`uv sync --frozen`), scientific outputs are byte-identical across repeated runs, across
`PYTHONHASHSEED` values, across `--batch_size`, across `--prefetch_targets`, and — unlike
the previous implementation — across `--cpus`. The contract is scoped to that lockfile:
pyhmmer, pyfamsa and pytrimal decide hit retention, alignment and serialised bytes.
(`<chunk>.log` carries timestamps and is excluded from that contract. HMM files omit the
`DATE` and `COM` lines, which are otherwise a wall-clock and an `argv` dump.)

> **The old pipeline's recruitment depended on how many CPUs it was given.** pyhmmer
> selects `parallel="targets"` whenever the query count is below the CPU count, which was
> every call in the old family-at-a-time loop. Each worker runs its own `Pipeline` over a
> slice of the database, and the merge concatenates each slice's *stored* hits while
> re-thresholding only the reporting flags. `Z` and `domZ` come out identical, but the
> stored list grows — and the old code iterated that raw list rather than `.reported`.
> Measured on the 50 000-sequence fixture under the old pinned `pyhmmer==0.11.1`:
> `len(TopHits)` goes `26/19/55` at `--cpus 1` to `27/19/56` at `--cpus 4`, while
> `.reported` stays `26/19/54` throughout.
>
> Forcing `parallel="queries"` fixes this: the answer is the same at any core count.
> Those extra stored hits were exactly the ones failing the reporting threshold, so
> reading `.reported` closes both halves of the problem at once.

### Two bugs fixed, and what they change

**Recruitment ignored `--recruit_evalue_cutoff`.** The old code iterated the raw
`TopHits`, which retains hits pyhmmer stored but did not report. Extraction now reads
`top_hits.reported`. On the small fixture, family `4497037939_1_144` used to recruit
sequence `6320430079`, which is stored but below the reporting threshold. Families are
correspondingly smaller: on that fixture, 32/19/65 members become 31/19/61. Same
families, same representatives, fewer spurious members.

**Recruitment depended on the CPU count**, as described above. Both halves are fixed, so
`--recruit_evalue_cutoff` now means what it says, on any machine.

Outputs are therefore **not** byte-compatible with the legacy script. Every difference is
enumerated in
[CHANGELOG.md](https://github.com/vagkaratzas/mgnifam/blob/main/CHANGELOG.md).

### Indexing a very large database

Easel buffers up to 2 GB of keys in RAM before spilling to an external sort, which then
needs scratch space in `TMPDIR` plus room for the final index (roughly
`n_sequences x (name_length + 16)` bytes). Size `TMPDIR` accordingly before indexing a
billion-record FASTA.

## Development

```bash
uv lock --check && uv sync --frozen
uv run pre-commit install
uv run pre-commit run --all-files
uv run pytest
```
