---
title: Computing clans
description: Family similarity and clans from residue-level overlap of family members.
---

`compute_clans` scores how much two families' members cover the same residues of the same
proteins, and groups families into clans. It reads only the member regions, so it works on
the output of either command and on any release, including one refreshed with
`update_families --skip_refine`, which never compares families.

```bash
uv run mgnifam compute_clans \
    --members output/1_families.tsv output/2_families.tsv \
    --fasta_file mgnifams_input.fa
```

## Inputs

`--members` takes one or more `family<TAB>member` files: `<chunk>_families.tsv`,
`<chunk>_updated_families.tsv`, or any file in that shape, plain or gzip. They are read as
one input, so a release made of updated chunks plus older families is just a list of
files. Every row is one member:

- `protein/start-end` is a region, 1-based and inclusive, split at the **last** slash, the
  same spelling the other commands write.
- A bare name with no trailing range is a whole protein.

`--fasta_file` is needed only when some member is a whole protein, for its length. The
FASTA may be gzipped, and its records are read only when bare members exist. The metadata
CSV is not enough: it has lengths for whole-protein *representatives* only, and other
members can be whole proteins too.

Family names must match `[A-Za-z0-9._-]+` and be neither `.` nor `..`, as in
[Updating families](/mgnifam/guides/update-families/#what-identity-means-here).

## What is measured

For family F, R(F) is the set of (protein, position) pairs its members cover. Each
protein's overlapping or adjacent regions are merged first, so repeat domains and
duplicate rows count once. For every pair of families sharing at least one protein:

| column | value |
|---|---|
| `shared_proteins` | proteins carrying both families |
| `protein_jaccard` | shared proteins over the union of the two families' proteins |
| `residue_jaccard` | \|R(A) ∩ R(B)\| / \|R(A) ∪ R(B)\| |
| `residue_containment` | \|R(A) ∩ R(B)\| / min(\|R(A)\|, \|R(B)\|) |

Positions are coordinates, not residue letters: `MGYP0000133100/10-98` covers positions 10
to 98 of that protein. Every pair that shares a protein is scored exactly. There is no HMM
search, no representative and no candidate filter.

The two Jaccards tell apart cases the protein-level measure cannot. Two families on
different domains of the same multi-domain proteins have a `protein_jaccard` of 1 and a
`residue_jaccard` of 0. A small family nested inside a large one has a low Jaccard but a
containment of 1.

## Clans

A clan is a connected component of the graph whose edges are the pairs with
`residue_jaccard >= --clan_jaccard` (default `0.5`, the release 1.0 threshold). Only
components of two or more families are written. The representative is the family with the
largest sum of `residue_jaccard` over its clan edges, and ties go to the earlier name.
Clans are numbered by size, largest first, then by representative name. Clan ids are
`<chunk>_clan_<n>` and are not stable across releases.

| flag | default | meaning |
|---|---|---|
| `--min_edge_jaccard` | `0.1` | write a pair whose `residue_jaccard` reaches this |
| `--min_edge_containment` | `0.5` | or whose `residue_containment` reaches this |
| `--clan_jaccard` | `0.5` | clan edge threshold, applied to every pair, not only written ones |
| `-n`, `--chunk_id` | `1` | prefix for the output files |
| `-p`, `--cpus` | `8` | worker processes |
| `--worker_memory_mb` | `512` | memory budget per worker |
| `--tmp_dir` | `$TMPDIR` | where the run creates its scratch directory |

Thresholds take at most six decimal places and are compared exactly, so a pair at exactly
`0.5` is a clan edge.

## Outputs

Written under `--output_dir` (default `output`):

| path | contents |
|---|---|
| `<chunk>_clan_edges.csv` | `family_1,family_2,shared_proteins,protein_jaccard,residue_jaccard,residue_containment`, one row per pair passing either edge threshold, `family_1 < family_2`, sorted |
| `<chunk>_clans.csv` | `clan_id,family_id,is_representative` |
| `<chunk>_clans_summary.json` | parameters and counts: rows, families, proteins, pair emissions, pairs sharing a protein, edges, near-redundant pairs (`residue_jaccard >= 0.95`), clans, singletons, clan size histogram |
| `<chunk>_clans.log` | run log |

Ratios are printed to six decimal places. Families sort by name in byte order, so `14`
sorts before `3`. Near-redundant pairs are counted, not removed: filter the edge list on
`residue_jaccard` to see them.

The summary is written last and atomically, and only after the run has finished reading
and checking every input. A run rejected for bad input (exit `2`) leaves the output
directory exactly as it was, including a previous run's summary. A run that fails (exit
`1`) leaves no summary, not even a previous run's, and keeps its log in the output
directory. So **a summary file means the edges and clans next to it are complete**. Edges, clans and summary are byte-identical across `--cpus`,
`--worker_memory_mb`, `--tmp_dir`, input order, how the input is split into files,
compression and `PYTHONHASHSEED`. The log has timestamps and is not.

The summary deliberately has no `"tool"` key, so the MultiQC module for the family
commands does not pick it up.

## Scale

The work is the number of family pairs on each protein, Σ k(k−1)/2, plus the number of
overlapping region pairs. It is not the number of family pairs. MGnifams release 1.1 has
2.37 billion member rows and about 60 billion family pairs over its proteins.

- **Memory** is one dense array of per-pair sums, 12 bytes per family pair (about 6 MB at
  1,000 families, 7.1 GB at 34,419), plus `--worker_memory_mb` per worker. The array is a
  memory-mapped file in scratch, not shared memory, so it works in containers with a small
  `/dev/shm`.
- **Scratch disk** is 28 bytes per member row, plus the pair array, plus 24 bytes per FASTA
  record when whole-protein members exist, plus about 10% headroom. That comes to about
  81 GB for release 1.1. Point `--tmp_dir` at a disk with that much room. The run creates its
  own directory inside it and removes only that directory, whether it succeeds or fails.
- **Time**: on one core, parsing runs at about 0.5 million rows a second and the pair phase
  at about 13 million pair emissions a second. Release 1.1 extrapolates to about four
  CPU-hours.

No protein is skipped, however many families it carries. The one input the budget cannot
handle is a single protein whose own rows exceed `--worker_memory_mb`. That run exits `1`
and asks for a larger budget.

## Limitations

- A whole unsliced protein literally named like a region, such as `X/1-5`, is read as
  positions 1–5 of `X`. The members file holds no length to tell them apart.
- Proteins are grouped by a 128-bit hash of their names. A collision would merge two
  proteins, which is negligible (about 10⁻²² at 10⁸ proteins).

## Compared with release 1.0 clans

Release 1.0 built clans from `hmmsearch` hits of each family model against the other
families' representatives, using a Jaccard over member *proteins*. This command differs in
four ways:

- It measures regions, not whole proteins, so neighbouring domains no longer join a clan.
- It scores every pair that shares a protein. In 1.0 a pair was scored only if one family's
  model hit the other's representative.
- It writes each pair once. 1.0 wrote a pair once per domain hit.
- Whole-protein representatives take part. In 1.0 they were never matched.
