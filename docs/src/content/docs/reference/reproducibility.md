---
title: Reproducibility
description: What is byte-identical, under which dependency set, and why the port differs from the legacy script.
---

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
