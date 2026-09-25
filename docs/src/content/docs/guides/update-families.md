---
title: Updating families
description: Refresh existing HMM families against a new sequence database.
---

database. It is the answer to "a new release came out" — you do not re-derive the
families from their original clusters, you search the models you already have.

```bash
uv run mgnifam update_families \
    --hmm_input previous_output/hmm \
    --fasta_file new_release.fa
```

`--hmm_input` is either a directory of `.hmm`/`.hmm.gz` files or a single multi-model
library (`hmm.lib.gz` works). Both forms produce identical output for the same models.
`--fasta_file` is uncompressed, for the same Easel reason as in [Generating families](/mgnifam/guides/generate-families/).

Every threshold flag from `generate_families` carries over with the same name and default.

| flag | meaning |
|---|---|
| `--skip_refine` | Recruit once and align. The model, seed MSA and RF line are unchanged, so only `hmm/` and `full_msa/` are written and `--max_seq_identity`, `--max_seed_seqs` and `--max_gap_occupancy` are inert. Without it, the full three-round refine loop runs and writes the complete artifact set. |
| `--chunk_id` | Labels the per-chunk aggregate files **only**. Family names come from the models, so nothing is renumbered. |

## What identity means here

A family keeps the name its model carries in its `NAME` field — `1_7` stays `1_7` across
releases, in the filenames *and* in every identity-bearing field inside the outputs. Two
consequences:

- `<chunk>_updated_metadata.csv`'s `family_id` column holds `1_7`, not a bare integer. That
  differs from `generate_families`, whose ids are a rank.
- Chunks sharing one output root **must own disjoint family names**. Nothing enforces it,
  because the names come from the input models rather than from `--chunk_id`.

A `NAME` must match `[A-Za-z0-9._-]+` and be neither `.` nor `..`. It is interpolated into
artifact paths and into CSV fields, and it arrives from a file this tool did not write.

## Update outputs

Per-family artifacts land in the same `hmm/`, `full_msa/`, `seed_msa/` and `rf/`
directories, named by family. Aggregates are `<chunk>_updated_*`: `families.tsv`,
`metadata.csv`, `discarded.csv`, `successful.txt`, `converged.txt`, `reps.fasta.gz`,
`delta.csv`, `stats.json` (see [MultiQC](/mgnifam/reference/outputs/#multiqc)), and `<chunk>_updated.log`.

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
`<chunk>_<rank>`, but an updated family keeps its model's name and `--chunk_id` never
appears in a per-family filename, so nothing on disk says which run wrote `hmm/1_7.hmm.gz`.

Input files must not overlap output paths, including through symlinks or hard links.
This is checked before writing, for both model directories and single-file libraries.
On an accepted retry, previous artifacts for the input family names are removed before
processing. Discarded families therefore leave no old models, and `--skip_refine` leaves
no seed/RF files from a previous refine run. A cleanup failure aborts the run.

## Cost

`hmmsearch` is `O(n_families x database)` and this command does not change that.
`--skip_refine` is one database pass per family; refining is up to three. Against a
billion-sequence release that term, not the alignment, is what to budget.
