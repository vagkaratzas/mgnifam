# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

*Reproducibility*

Repeated runs are byte-identical: same inputs, same lockfile, same outputs —
regardless of `--cpus`, `--batch_size`, `--prefetch_targets` or `PYTHONHASHSEED`.

That guarantee is scoped to the dependency set resolved in the committed `uv.lock`
(`uv sync --frozen`). pyhmmer, pyfamsa and pytrimal decide hit retention, alignment
and serialised bytes, so installing from PyPI — which resolves within the declared
version ranges instead — may produce different results on a different resolution.

## [2.1.0.dev0] - unreleased

### Added

- **`mgnifam update_families`**, a subcommand that refreshes families which already exist as
  HMMs against a new database. It searches each stored model rather than
  re-deriving the family from its original clusters, so a family keeps its identity across
  releases and picks up every fix in the shared code path.

  `--hmm_input` takes a directory of `.hmm`/`.hmm.gz` files or a single multi-model library;
  both forms produce identical output. `--skip_refine` recruits once and aligns, leaving the
  model, seed MSA and RF line untouched; without it the full three-round loop runs.

  A family's name comes from its model's `NAME` field and is preserved verbatim, so
  `<chunk>_updated_metadata.csv`'s `family_id` column holds a name like `1_7` rather than a
  rank, and `--chunk_num` labels only the per-chunk aggregate files. Chunks sharing an output
  root must own disjoint family names.

  `<chunk>_updated_delta.csv` reports one row per family — the model length before and after,
  round 1's recruitment, the final size, retention, rounds run, and the outcome. `no hits in
  the new database` is distinguished from `low complexity model - confounding cluster`: the
  first means the model found nothing in the new release, the second that it found hits and
  none were long enough.

  Give each run its own `--output_dir`. Re-running the same models into the same directory
  is fine, which is what a retried chunk does. Running a *smaller* set of models over a
  directory that still holds a larger one is refused, because the dropped families'
  artifacts would be left beside aggregates that no longer list them, and an updated
  family's name carries nothing that says which run wrote it.

### Changed

- **Breaking:** a chunk that completes after containing one or more internal family
  errors now exits 3 instead of 0. Every crashed family is still recorded in
  `<chunk>_discarded.csv`, and the chunk continues through all remaining families before
  exiting, so its output is complete and self-consistent. Callers that previously ignored
  the status see unchanged artifacts; callers that check it must now re-run the degraded
  chunk or explicitly accept those lost clusters. The final `DONE.` log line includes the
  number of crashed families. Exit 1 remains fatal and marks incomplete output that must
  not be consumed; exit 2 remains an `argparse` usage error.

### Fixed

- Accepted `update_families` retries clear previous artifacts for the input family
  names before processing. A later discard no longer leaves a stale HMM/alignment,
  and recruit-only retries remove obsolete seed/RF artifacts. Foreign family names
  are still refused; cleanup failures abort before aggregate files are truncated.
- `update_families` rejects input files that overlap an output destination, including
  single-file HMM libraries and symlink/hard-link aliases, before creating any output.
  Previously a successful update could overwrite the library needed to retry it.

- **Breaking, narrowly:** if writing a row to `<chunk>_discarded.csv` fails, the run now
  exits 1 instead of trying again.

  Retrying was unsafe. A write can fail *after* putting part of the row on disk, and the
  retry appended the whole row again — leaving one and a half rows for one family. The
  entry below fixed exactly this for successful families in 2.1.0.dev0 but missed the
  discard branch; both are now handled the same way. A run that writes its output
  successfully is unaffected, and `generate_families`' output is unchanged byte-for-byte.
- A failure after a successful family first appended to shared chunk output could be
  contained as a discard, putting the same family in generated and discarded outputs.
  Shared-output commit failures now bypass per-family containment and exit 1, so exit 3
  structurally means the degraded output is coherent.
- A failed per-family artifact write could leave earlier artifacts on disk beside its
  discard row. Artifact writes now roll back every path already attempted; a rollback
  failure names the paths it could not remove and exits 1 without masking the original
  write error.
- Discarded families now release their seed and full alignments, HMM, hit records, and
  checked-sequence set immediately. This prevents a batch from pinning its dominant memory
  objects and makes containment of an oversized family's `MemoryError` capable of
  reclaiming that memory before processing continues.
- **Breaking:** representative selection now uses the highest-scoring reported domain of
  the top-ranked hit. HMMER reports a hit's domains in positional order, so a short
  leftmost fragment could previously become row 0 of the full MSA, producing incorrect
  representative metadata or discarding a healthy family. `extract_records` now orders
  each hit's reported domains by score before alignment. This intentionally changes
  scientific output for affected multi-domain hits and diverges from the legacy script;
  the regression is covered by `test_extract_records_puts_top_scoring_domain_first`.

## [2.0.0] - 2026/07/29

A major version because two documented behaviours change: the two per-chunk CSVs gain a
header row, and an explicitly supplied `--fasta_index` is no longer built when it is
missing. Both are listed under *Changed* below. Scientific outputs are unaffected for
databases whose accessions carry no underscores beyond their slice bounds, which is every
MGnifams accession — see *Fixed*.

### Changed

- **Breaking:** `<chunk>_metadata.csv` and `<chunk>_discarded.csv` now start with a header
  row, so both load with `pandas.read_csv` without `header=None` and a hand-maintained
  `names=`. Any reader that does not skip it will treat the header as data.

  ```
  family_id,full_msa_size,protein,region,length,sequence,consensus,converged
  representative,reason,value
  ```

  Column order is unchanged — only the row is new. Headers are written by `main` before
  the first batch rather than by `emit_family`, so a chunk that produces no families and
  a chunk that discards nothing both still yield a parseable file instead of an empty one.

- **Breaking:** an explicit `--fasta_index` is now used exactly as given and never
  rebuilt, which makes a single index safely shareable across parallel chunk tasks — the
  case the flag exists for. Pointing the flag at a path that does not yet exist used to
  build an index there, as the README documented; it is now rejected as a typo.
  `resolve_index` previously rebuilt any index whose mtime predated its FASTA's,
  treating a caller-supplied path as a cache it owned. An mtime comparison is not a
  staleness signal for a path this process did not create, and two consequences followed:

  - Copying, restoring from a backup or archive, or rebuilding the FASTA from identical
    bytes all reorder the two timestamps without invalidating anything. `Path.stat()`
    also follows symlinks, so a linked index reports its *target's* timestamp rather than
    the link's, and the FASTA and the index may be linked from unrelated places whose
    relative order says nothing about whether one describes the other. Every concurrent
    chunk sharing the index then re-indexed the whole database at once — precisely the
    cost `--fasta_index` is meant to avoid. The shared index itself was never corrupted
    (`os.replace` hits the link, not its target), only the work wasted.
  - An index kept somewhere the process cannot write — a shared reference directory, a
    read-only mount — could not be rebuilt at all: `build_ssi_index` creates its scratch
    directory beside the destination, so the run died with a `PermissionError` from
    `mkdtemp` after the chunk had already started.

  A missing `--fasta_index` is now rejected by `validate_inputs` as a typo instead of
  being absorbed as a build at the misspelled path, and an index that does not match its
  FASTA still surfaces at the first fetch as `IndexMismatchError`. The default path under
  `--output_dir` is unchanged: nothing else owns it, so it is still built and refreshed
  automatically. Guard: `test_supplied_fasta_index_is_used_as_given_and_never_rebuilt`.

  Verified interoperable with `esl-sfetch --index` in both directions. The two indexes
  are not byte-identical — `esl-sfetch` also fills `data_offset` and `record_length`, and
  sizes the filename field from the path it was given — but `generate_families` performs
  only whole-record fetches, for which they are interchangeable.

### Fixed

- Protein names containing underscores are no longer misread as slice bounds. A database
  record named `<protein>_<start>_<end>` is a slice of a parent protein, and
  `parse_protein_name` recovered the parent by requiring `split("_")` to yield exactly
  three fields. Three name shapes broke on that test, all reported by users running the
  tool on databases whose accessions are not bare MGnifams integers:
  - `contig_1_gene_2_88_140` — a genuine slice of a protein whose own name contains
    underscores. Four fields failed the length test, so the row was treated as a whole
    protein and renamed to `contig`, silently discarding everything after the first field.
  - `contig_1_gene_x` — three fields, but not numeric ones. `int()` raised, `family_guard`
    caught it, and the entire family was recorded as an internal-error discard.
  - `scaffold_12_34` — three numeric fields that are part of the name, not bounds. The row
    was reported at invented parent coordinates.

  Name splitting now happens from the right, via the new `split_slice_name()`, and the
  trailing two fields are accepted as bounds only when they are integers *and* span
  exactly as many residues as the record holds. That span test is the disambiguator: it
  holds for every real slice by construction, and rejects a coincidental `_12_34`. A name
  that fails it is treated as a whole protein, which is the safe reading — no crash, no
  truncation. Cost is one `len()` on a string already fetched: measured at +194 ns against
  the 11 µs `parse_protein_name` spends per row, 82% of which is its SSI fetch.

  This diverges from `reference/legacy_generate_families.py`, which has the same defect.
  Guard: `test_underscores_in_protein_names_are_not_mistaken_for_slice_bounds`.

## [1.0.0] - 2026/07/22

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
- A time-throttled heartbeat, logged every 300 s during a database pass.

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
  directory and a logger writing straight to `<chunk>.log` in the output root.
- Stale `<chunk>_<n>` artifacts are cleared at start-up, so a rerun producing fewer
  families no longer leaves the surplus behind.
- Python `>= 3.13`. Dependencies updated to `pyhmmer>=0.12.1,<0.13`,
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
