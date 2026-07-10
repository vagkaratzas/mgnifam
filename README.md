# mgnifam

Iterative HMM-based protein family generation over very large sequence databases.

Given a chunk of MMseqs2 clusters and a protein FASTA, `generate-families` builds an
HMM from each cluster, recruits new members from the whole database, re-aligns, and
either converges on a family or discards the cluster. It is the core algorithm of the
[`mgnifams`](https://github.com/vagkaratzas/mgnifams) Nextflow pipeline, extracted into
a standalone, tested package.

## Install

```bash
uv sync
```

## Usage

```bash
uv run generate-families \
    --clusters_chunk clusters.tsv \
    --fasta_file mgnifams_input.fa \
    --cpus 8 \
    --chunk_num 1 \
    --discard_min_rep_length 100 \
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

| flag | default | meaning |
|---|---|---|
| `--fasta_index` | `./<fasta basename>.ssi` | Path to an Easel SSI index. Built automatically if absent. |
| `--batch_size` | `max(2 * cpus, 16)` | How many families are searched per `hmmsearch` wave. |

**On a production run, build the index once and share it.** Every chunk task would
otherwise re-index the whole database:

```bash
# once, upstream
uv run python -c "from mgnifam.generate_families import build_ssi_index; \
                  build_ssi_index('db.fa', 'db.fa.ssi')"
# then, per chunk
uv run generate-families --fasta_index db.fa.ssi ...
```

## Outputs

Written under the current working directory, keyed by `--chunk_num`:

| directory | contents |
|---|---|
| `seed_msa_sto/` | `<chunk>_<id>.sto.gz` — seed alignment |
| `full_msa_sto/` | `<chunk>_<id>.sto.gz` — full alignment |
| `hmm/` | `<chunk>_<id>.hmm.gz` — the family model |
| `rf/` | `<chunk>_<id>.txt` — reference-annotation line |
| `family_reps/` | `<chunk>.fasta.gz` — one representative per family |
| `refined_families/` | `<chunk>.tsv` — `family_id<TAB>sequence` |
| `family_metadata/` | `<chunk>.csv` |
| `successful_clusters/` | `<chunk>.txt` |
| `discarded_clusters/` | `<chunk>.csv` — `representative,reason,value` |
| `converged_families/` | `<chunk>.txt` |
| `logs/` | `<chunk>.txt` |

Family ids are a 1-based rank among *successful* families, in cluster-file order.

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
   re-filtered.

On top of that, the entire FASTA was held in RAM twice — once as a
`DigitalSequenceBlock` and once as a Python `dict` of `DigitalSequence` objects. Both
are gone: targets stream from disk, and random access goes through an Easel SSI index.
Families are searched in batched waves so `hmmsearch` still saturates every core.

## Reproducibility

Scientific outputs are byte-identical across repeated runs, across `PYTHONHASHSEED`
values, across `--batch_size`, and — unlike the previous implementation — across
`--cpus`. (`logs/` carries timestamps and is excluded from that contract.)

> **The old pipeline's results depended on how many CPUs it was given.** pyhmmer selects
> `parallel="targets"` whenever the query count is below the CPU count, which was every
> call in the old family-at-a-time loop. Each worker then ran its own `Pipeline` over a
> slice of the database and applied a *per-slice* `domZ` at the inclusion threshold, so
> extra domains were reported. Streaming targets forces `parallel="queries"`, which
> reproduces the single-threaded answer at any core count. Output from this package
> matches the old script run with `--cpus 1`, not its production output.

## Development

```bash
uv run pre-commit install
uv run pre-commit run --all-files
uv run pytest
```
