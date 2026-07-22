# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

*Disclaimer*

Byte-level reproducibility is only claimed for the dependency set resolved in the
committed `uv.lock`. pyhmmer, pyfamsa and pytrimal determine hit retention, alignment
and serialised bytes, so a change in any of them may change outputs.

## [1.0.0] - unreleased

First release of `mgnifam` as a standalone package. The algorithm is a port of
`generate_families.py` from the [mgnifams](https://github.com/vagkaratzas/mgnifams)
Nextflow pipeline, preserved at `reference/legacy_generate_families.py`.

Outputs are **not** byte-compatible with the legacy script. The differences are
enumerated under *Changed* and *Fixed*, and every one of them is intentional.

### Added

- `mgnifam` console script with subcommand dispatch — `mgnifam generate_families ...`.
  Room reserved for `remove_redundant` and `merge_families`. `python -m mgnifam` works
  identically.
- `clip_env_ends()`: trims the leading and trailing `#=GC RF` `.` columns — the N- and
  C-terminal envelope overhangs — from every seed-path alignment. Interior `.` columns
  are inserts between match states and are kept. Not applied to the final full MSA.
- Discard rule: a family with `<= 2` sequences after redundancy trimming is dropped. A
  separate `< 2` guard before the initial build converts `hmmbuild`'s `eslEMEM (status
  code 5)` crash into a clean discard. The final hand build needs no such guard: round 1
  can never converge, so every seed reaching it has already passed the `<= 2` rule.
- `--fasta_index`: reuse a pre-built Easel SSI index. Production runs should build one
  index upstream and share it, or every chunk task re-indexes the whole database.
- `--output_dir`: place every generated file and directory under one root. Defaults to
  `output`.
- `--prefetch_targets`: load the database into RAM instead of streaming it per query.
  Faster, `O(database)` memory, and measured byte-identical to streaming.
- `--batch_size`: how many families are searched per `hmmsearch` wave. Defaults to
  `2 * cpus` and must stay `>= cpus`, since `hmmsearch` uses `min(cpus, n_queries)`
  workers.
- Gzip compression for the bulk outputs: `seed_msa/*.sto.gz`, `full_msa/*.sto.gz`,
  `hmm/*.hmm.gz`, `<chunk>_reps.fasta.gz`.
- Input validation before any index or output directory is created: `chunk_num` must
  match `[A-Za-z0-9._-]+` (it is interpolated into output paths), the FASTA must be
  uncompressed, percentages must lie in `[0, 1]`, length bounds must be ordered, and
  every cluster TSV row must hold exactly two non-empty fields.
- A time-throttled heartbeat, logged every 60 s during a database pass.
- `pytest` suite (32 tests) over real 50,000- and 26,949-sequence fixtures, and
  `pre-commit` with `ruff`.

### Changed

- **Flatter output layout.** Only the per-family artifacts keep a directory
  (`seed_msa/`, `full_msa/`, `hmm/`, `rf/`; the two MSA directories lose their `_sto`
  suffix, the extension already says it). The seven one-file-per-chunk outputs sit flat
  in the output root as `<chunk>_families.tsv`, `<chunk>_discarded.csv`,
  `<chunk>_successful.txt`, `<chunk>_converged.txt`, `<chunk>_metadata.csv`,
  `<chunk>_reps.fasta.gz` and `<chunk>.log` — a directory holding a single file was
  never carrying information.
- **Recruitment no longer depends on `--cpus`.** `hmmsearch` is always called with
  `parallel="queries"`. See *Fixed*.
- The FASTA is never loaded into memory. Targets stream from a `SequenceFile`, and
  individual sequences are fetched through an Easel SSI index built atomically inside a
  private temporary directory. Both the `DigitalSequenceBlock` and the Python
  `dict` of `DigitalSequence` objects are gone. `--prefetch_targets` restores the
  former on request.
- Families are processed in round-major waves rather than one at a time, so that
  streaming targets still use every core.
- HMM files omit the `DATE` and `COM` lines, which held a wall-clock timestamp
  (locale-formatted) and a dump of `sys.argv`. The family name survives as `NAME`.
  Gzip members carry neither an mtime nor a filename. Together these make the outputs
  byte-reproducible across runs, `PYTHONHASHSEED` values, `--cpus`, `--batch_size` and
  `--prefetch_targets`.
- `hmmalign` and `pyfamsa` now receive `--cpus`. Both previously used every physical
  core on the node, ignoring the task allocation. FAMSA's output is thread-count
  independent, so this changes resource use, not results.
- Output file handles are opened once on a `contextlib.ExitStack`; the shared
  `tmp/seed_msa.sto`, `tmp/full_msa.sto` and `execution.log` paths, which collided when
  chunks shared a working directory, are replaced by a per-invocation temporary
  directory and a logger writing straight to `logs/<chunk>.txt`.
- Stale `<chunk>_<n>` artifacts are cleared at start-up, so a rerun producing fewer
  families no longer leaves the surplus behind.
- Python `>= 3.12`. Dependencies updated to `pyhmmer>=0.12.1,<0.13`,
  `pyfamsa>=0.7.0,<0.8`, `pytrimal>=0.8.5,<0.9`, `numpy>=2.5.1,<3`. Upper bounds are
  deliberate; see the note at the top of this file.
- Dropped `pandas` (the cluster table is now an insertion-ordered `dict`; note that
  `itertools.groupby` would have been wrong, as it only groups *contiguous* rows while
  the legacy code gathered every row for a representative) and `biopython` (never
  imported).
- pyhmmer 0.12 exposes names, `MSA.reference` and `HMM.name` as `str` rather than
  `bytes`; every `encode`/`decode` around them is gone.

### Fixed

- **Sequences that failed `--recruit_evalue_cutoff` were being recruited.** The legacy
  code iterated the raw `TopHits`, which retains hits pyhmmer stored but did not report.
  Extraction now reads `top_hits.reported` and `hit.domains.reported`. On the small
  fixture, query `4497037939_1_144` stores 55 hits and reports 54; the 55th,
  `6320430079`, was previously recruited. Families shrink accordingly — on the fixture,
  from 32/19/65 members to 31/19/61.
- **Results depended on the CPU count.** pyhmmer selects `parallel="targets"` whenever
  the query count is below `--cpus`, which was every call in the old family-at-a-time
  loop. Its merge concatenates each worker's stored hits while re-thresholding only the
  reporting flags, so `len(TopHits)` grew with `--cpus` (measured `26/19/55` at
  `--cpus 1` versus `27/19/56` at `--cpus 4`, under the old pinned `pyhmmer==0.11.1`)
  while `.reported` stayed constant. Combined with the raw iteration above, recruitment
  varied with the machine. Both halves are fixed.
- **`run_initial_msa` was `O(database x members)` per family.** `map(str.encode, members)`
  was re-created inside a comprehension's condition, so a membership test became a full
  linear scan of the cluster for every sequence in the database — repeated for every
  cluster. It is now one indexed lookup per member. This was the dominant cost of the
  ~8-month full-database run.
- **Cluster selection was `O(N^2)`.** The cluster table was boolean-masked and dropped
  from once per family. It is now a single grouping pass.
- **The exit-branch `hmmsearch` was redundant.** It re-ran the search the preceding
  round had just performed with the identical HMM, differing only in a post-filter. Its
  hit records are cached and re-filtered, saving one full database pass per family.
- **The exported HMM did not match the model that built the family.** In the
  `family_iteration > 3` exit path the legacy loop ran round 3's redundancy trim and
  kept its result as the seed MSA, then left the loop without ever building a model
  from it. The family's hits and full MSA came from round 3's model (built from round
  2's seed), while the HMM, RF line and seed MSA written to disk came from that trim's
  output — an alignment nothing had ever searched with. The trim now runs only when a
  further `hmmbuild` will consume it, so on both exit paths the seed MSA, the HMM and
  the full MSA describe one round. This also removes round 3's `<= 2` and representative
  length checks, which no longer gate anything; `finish` re-checks the length on the full
  MSA. On the `mgnifams_v2` fixture: same 12 families, same 8 converged, same members and
  same full MSAs; the seed MSA, RF and HMM of the 4 families that exhausted three rounds
  change, and with them their consensus in `family_metadata`.
- **Discarded families were recorded as converged.** `converged_families` now contains
  only successful family ids, written after membership and length checks have passed.
- **`clip_ends()` discarded a column of every model.** It trims low-occupancy columns
  from both ends of an alignment, which it did correctly, but built an end-exclusive
  `range(start, end)` over *inclusive* bounds — so the last column that passed the
  occupancy threshold was thrown away with the failing ones. Every seed alignment, and
  therefore every HMM, was one match state short. Separately, when no column passed the
  threshold, `np.argmax` over an all-`False` array returned 0 for both scans, so the
  function reported the full span and silently trimmed nothing but the final column;
  `calculate_trim_positions` now returns `None` and the alignment is returned unchanged.

  On the small fixture the recovered column lengthens each family's model
  (representative lengths 104 → 105 and 112 → 113), and family `4497037939_1_144` now
  converges, having recruited nothing new once its model stopped losing a column. Family
  count, representatives and membership are unaffected.
- Output file handles were left unclosed, and `discard_value = 0.0` was re-initialised
  at four call sites.
- **Renumbering truncated coordinates and could drop rows.** `renumber_sto_msa` re-read
  the Stockholm file it had just written and rewrote each name by string replacement, so
  the new name had to fit the old name's column width; a name that grew was cut to
  length. On the `mgnifams_v2` fixture that shipped `4454641265/11-11` for a region that
  ends at 110, and `909822872` for a region that is `909822872/1-112` — the coordinates
  were eaten outright. The padding also leaked into `refined_families`, whose every row
  carried trailing spaces.

  Renaming now happens on the MSA object before it is written, so pyhmmer sizes the name
  column and no name is ever trimmed. `parse_protein_name` takes a row and its alignment
  instead of six positional integers, and `renumber_sto_msa` is gone along with the two
  temporary files per family.
- **Repeat domains collided and one was silently dropped.** A row's position in its
  record was recovered with `str.find` over the whole record, which returns the first
  match, so two identical repeat domains of one protein resolved to the same name and the
  second was skipped as a duplicate. The search is now bounded by the row's own envelope,
  which its name already carries, so the domains stay distinct and the duplicate skip has
  been removed. Residues that are not in their envelope raise instead of silently
  renumbering from `find`'s `-1`.
- **Stockholm output was anonymous and carried hmmalign's posteriors.** Every `#=GF`,
  `#=GS` and `#=GR` line was dropped, including the `#=GF ID`. Building the output from
  the MSA object restores `#=GF ID <family>` in both the seed and full alignments, and
  drops `#=GR PP`/`#=GC PP_cons` by construction rather than by line filtering.

  Across the `mgnifams_v2` fixture these change no membership: the same 12 families, the
  same 405 rows, and byte-identical `hmm`, `rf`, `family_metadata`, `family_reps`,
  `converged_families`, `successful_clusters` and `discarded_clusters`. The alignments
  gain their `#=GF ID`, lose the padded name column, and correct the two truncated names
  above.

## [0.1.0] - legacy

The original `mgnifams/bin/generate_families.py`, vendored unmodified at
`reference/legacy_generate_families.py` and kept as the behavioural reference for this
port.

- Iteratively builds an HMM per cluster, recruits members from the whole database with
  `hmmsearch`, re-aligns with `hmmalign`, removes redundancy with `pytrimal`, and either
  converges or is forced out after three rounds.
- Discards clusters on representative length, seed membership, or a low-complexity model.
- Emits seed and full Stockholm alignments, an HMM, an RF line, family metadata and a
  representative FASTA per chunk.
- Ran on `pyhmmer=0.11.1`, `numpy=2.3.2`, `pandas=2.3.2`, `python=3.13.5`, plus
  `biopython==1.85`, `pyfamsa==0.6.0`, `pytrimal==0.8.2`.
- Held the entire target FASTA in memory twice, as a `DigitalSequenceBlock` and as a
  Python `dict` of `DigitalSequence` objects.
