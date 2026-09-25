---
title: Outputs
description: Every file mgnifam writes, its exit status, and the MultiQC summary.
---

These are the outputs of `generate_families`. `update_families` writes its own set, described
in [Updating families](/mgnifam/guides/update-families/#update-outputs).

Written under `--output_dir` (default: `output`), keyed by `--chunk_id`:

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
| `<chunk>_stats.json` | run summary for [MultiQC](#multiqc); present only after a completed run |
| `<chunk>.log` | run log |

Family ids are a 1-based rank among *successful* families, in cluster-file order.

Both CSVs carry a header row, so they load with `pandas.read_csv` as they are:

| file | columns |
|---|---|
| `<chunk>_metadata.csv` | `family_id,full_msa_size,protein,region,length,sequence,consensus,converged` |
| `<chunk>_discarded.csv` | `representative,reason,value` |

`protein` is quoted, with embedded quotes doubled; a `protein` or `representative`
containing a comma or a quote is escaped, so both files parse with a standard CSV reader.
Literal slashes stay in `protein`, including punctuation after a slash: `protein/v1,variant`
is one protein field. Only a trailing coordinate range spanning the emitted sequence is
separated into `region`.
`region` is `<start>-<end>` on the parent protein, or `-` when the
representative spans a whole unsliced record. Those two columns together are the
`<base>/<start>-<end>` spelling from [Sequence names](/mgnifam/guides/generate-families/#sequence-names), which is also how `<chunk>_reps.fasta` names its
records. The representative is the highest-scoring
reported domain of HMMER's top-ranked hit. The header is written before the run starts, so
a chunk that produces no families still yields a parseable file.

## Exit status

The status describes whether the output is safe to consume, not only whether the process
stopped:

| Code | Meaning |
|---|---|
| `0` | Chunk completed. Every family landed on exactly one side of the split (discarded or successful). Output is complete and safe to consume. |
| `1` | Fatal: the run died before finishing. **Output is incomplete and must not be consumed** — re-run the chunk. This is what a dead output sink (ENOSPC, EIO) produces, because the discard re-emit cannot record its own failure. |
| `2` | Usage error from `argparse`. Nothing ran. |
| `3` | Chunk completed, but one or more families died of an internal error and were recorded as discards. Output is complete and self-consistent, but those clusters produced no family — re-run the chunk once the cause is fixed, or accept the loss. |

## MultiQC

Every completed run (exit `0` or `3`) writes one `<chunk>_stats.json`: counts of families
in, successful, discarded, converged and crashed; discard reasons; and `{value: count}`
histograms of full-MSA size, model length and representative length. It is read back from
the chunk's own CSVs, so it never disagrees with them. One chunk is one MultiQC sample, and
MultiQC merges a pipeline's chunks into one report.

The file is written last and atomically, and a rerun removes the previous one before it
changes anything else. **Its presence therefore means the directory holds a completed,
consumable run**; exit `1` leaves none. It records only the flags that change results, so
it is byte-identical across `--cpus`, `--batch_size` and `--prefetch_targets` like the
other outputs. `"tool": "mgnifam"` is its first key, and `schema_version` changes only
when its shape does.

`update_families` writes the same summary as `<chunk>_updated_stats.json`, read back from
its delta and metadata CSVs. It adds `skip_refine` to the recorded flags and three
histograms: `model_length_change` (after minus before), `rounds_run` and `retention`.
Retention keys are the exact `delta.csv` values, and a family with no retention is left
out. In both files `converged` counts successful families only.
