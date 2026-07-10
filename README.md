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
| `--batch_size` | `2 * cpus` | How many families are searched per `hmmsearch` wave. Keep it `>= cpus`. |
| `--prefetch_targets` | off | Load the database into RAM once instead of streaming it per query. Faster, `O(database)` memory, **identical results**. |

Streaming re-reads and re-parses the database once per query. `--prefetch_targets`
parses it once and keeps it in RAM; the results are byte-identical either way, so the
flag is purely a memory-vs-time dial. Leave it off unless the database fits comfortably
in RAM.

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
pyhmmer, pyfamsa and pytrimal decide hit retention, alignment and serialised bytes. (`logs/` carries timestamps and is excluded from that
contract. HMM files omit the `DATE` and `COM` lines, which are otherwise a wall-clock
and an `argv` dump.)

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
> Forcing `parallel="queries"` reproduces the single-threaded answer at any core count.
> On that fixture the extra hit is discarded downstream by the envelope-length filter, so
> the *final* families were unaffected — but the divergence is real at the recruitment step.

**What "matches legacy" does and does not mean.** The *search and recruitment semantics*
target the old script at `--cpus 1`. The final artifacts deliberately differ: envelope-end
clipping of seed MSAs, the new `<= 2`-sequence discard rule, gzip framing, the omitted HMM
`DATE`/`COM` lines, and upgraded dependencies all change bytes, and clipping can change
family membership and counts. This is not a byte-for-byte drop-in for the old outputs.

### A known bug, preserved

The old code recruits hits that **failed `--recruit_evalue_cutoff`**: it iterates the raw
`TopHits` list, which retains stored-but-unreported entries. On the small fixture, family
`4497037939_1_144` recruits sequence `6320430079`, which is stored but not reported, even
at `--cpus 1`.

This package **reproduces that behaviour deliberately**: it preserves legacy-at-`--cpus 1`
raw recruitment semantics rather than silently changing the science. Switching extraction to
`top_hits.reported` would make `--recruit_evalue_cutoff` mean what it says. That is a
one-line change, and a decision for the maintainers.

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
