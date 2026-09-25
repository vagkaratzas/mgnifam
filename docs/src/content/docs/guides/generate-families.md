---
title: Generating families
description: Build HMM families from a chunk of MMseqs2 clusters.
---

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
    --chunk_id 1 \
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

## Sequence names

A record that is a slice of a larger protein may say so in either of two spellings, and
both are read identically:

| spelling | example | parent protein | region |
|---|---|---|---|
| `<protein>_<start>_<end>` | `3387826881_356_472` | `3387826881` | 356–472 |
| `<base>/<start>-<end>` | `3387826881/356-472` | `3387826881` | 356–472 |

The second is the form this tool *emits*, so `<chunk>_reps.fasta` from one release can be
used directly as the database for the next without its coordinates being lost. The base
keeps any slashes it carries: `3387826881/v1/356-472` is region 356–472 of the protein
`3387826881/v1`.

Bounds are read as coordinates only if they span the record exactly. `scaffold_12_34`
holding 15 residues is a whole protein named `scaffold_12_34`, not residues 12–34 of
`scaffold`. Anything else is identity and is kept whole — `3387826881/356_472`,
`3387826881/356`, `3387826881/v1` and `3387826881/356-472-243` are four distinct protein
names, none of them carrying a region.

Any other character is allowed in a name, including further slashes. Names are never
split on their first slash, so two records sharing a prefix stay distinct.
`X` and a literal `X/1_10` remain distinct even when the same family recruits residues
1–10 of `X` alongside the complete `X/1_10` record. Literal percent sequences such as
`%2F` are preserved too, in both update modes and in all emitted identities.

No name is *reserved*, but the slice spelling is not inert either. Whether a record is
independent of `3387826881` depends on which spelling it uses and on its own length:

| record, alongside `3387826881` | length | read as |
|---|---|---|
| `3387826881/356_472` | any | an unrelated protein — underscore is not the slice separator |
| `3387826881/356-472` | 117 | region 356–472 **of** `3387826881`, by its own declaration |
| `3387826881/356-472` | anything else | an unrelated protein — the bounds do not span it |

The middle row is the round-trip working as intended: a record that says it is a region
of `3387826881` is reported at those parent coordinates, exactly as the corresponding
residues of `3387826881` itself would be. If a database contains both, the same residues
are the same protein region and get the same name — they are not two things. Include the
parent and its own slices in one database only if that is what you mean.

## Optional flags

Pass every threshold explicitly on a production run. The defaults exist for ad-hoc use;
relying on them means a forgotten flag produces a plausible-looking family set instead
of an error.

| flag | default | meaning |
|---|---|---|
| `--cpus` | `8` | Threads for FAMSA, `hmmsearch` and `hmmalign`. |
| `--chunk_id` | `1` | Namespace for this chunk: it prefixes every output file and directory, and every family is named `<chunk_id>_<rank>`. Any string matching `[A-Za-z0-9._-]+` — it need not be numeric. |
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

## Indexing a very large database

Easel buffers up to 2 GB of keys in RAM before spilling to an external sort, which then
needs scratch space in `TMPDIR` plus room for the final index (roughly
`n_sequences x (name_length + 16)` bytes). Size `TMPDIR` accordingly before indexing a
billion-record FASTA.
