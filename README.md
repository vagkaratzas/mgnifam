<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/logo-dark.png">
    <img src="assets/logo.png" width="600" alt="mgnifam">
  </picture>
</p>

# mgnifam

Iterative HMM-based protein family generation over very large sequence databases.

Given a chunk of MMseqs2 clusters and a protein FASTA, `mgnifam generate_families` builds an
HMM from each cluster, recruits new members from the whole database, re-aligns, and
either converges on a family or discards the cluster. It is the core algorithm of the
[`mgnifams`](https://github.com/vagkaratzas/mgnifams) Nextflow pipeline, extracted into
a standalone, tested package.

## Install

```bash
uv sync
```

Requires Python >= 3.13. Verify with `mgnifam --version`.

## Usage

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
| `--fasta_index` | `<output_dir>/<fasta basename>.ssi` | Path to an Easel SSI index. Built automatically if absent. |
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
# once, upstream
uv run python -c "from mgnifam.generate_families import build_ssi_index; \
                  build_ssi_index('db.fa', 'db.fa.ssi')"
# then, per chunk
uv run mgnifam generate_families --fasta_index db.fa.ssi ...
```

`generate_families` is the only subcommand today. `mgnifam --help` lists them, and
`python -m mgnifam` is equivalent to the console script.

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
| `<chunk>_discarded.csv` | `representative,reason,value` |
| `<chunk>_converged.txt` | ids of successful families that converged naturally |
| `<chunk>.log` | run log |

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
enumerated in [CHANGELOG.md](CHANGELOG.md).

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
