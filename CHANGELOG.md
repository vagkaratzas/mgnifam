# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

*Reproducibility*

Repeated runs are byte-identical — same inputs, same outputs — regardless of `--cpus`,
`--batch_size`, `--prefetch_targets` or `PYTHONHASHSEED`. The guarantee holds for the
dependency versions pinned in `uv.lock` (`uv sync --frozen`). Installing from PyPI resolves
pyhmmer, pyfamsa and pytrimal within their declared ranges instead, and a different
resolution may change results.

## [2.1.0.dev0] - unreleased

### Added

- **`mgnifam update_families`** refreshes existing family HMMs against a new database by
  searching the stored models, so each family keeps its name across releases.
  - `--hmm_input` takes a directory of `.hmm`/`.hmm.gz` files or one multi-model library;
    both give identical output.
  - `--skip_refine` recruits once and aligns, leaving model, seed MSA and RF line unchanged.
    Without it, the full three-round refine loop runs.
  - Family names come from each model's `NAME` field, so `family_id` in
    `<chunk>_updated_metadata.csv` holds a name like `1_7`, not a rank. Chunks sharing an
    output directory must own disjoint family names.
  - `<chunk>_updated_delta.csv` reports, per family, model length before and after, round-1
    recruits, final size, retention, rounds run and outcome.
  - Use a separate `--output_dir` per run. Retrying the same models in place is allowed and
    clears that run's previous artifacts first; running a smaller model set over a directory
    holding a larger one is refused. Inputs that overlap an output path, including through
    symlinks or hard links, are refused before anything is written.

### Changed

- **Breaking:** a chunk in which one or more families failed with an internal error now
  exits `3` instead of `0`. The chunk still finishes and its output is complete and
  consistent; the failed clusters are listed in `<chunk>_discarded.csv`, and the final
  `DONE.` log line reports how many there were. Re-run the chunk or accept the loss. See
  *Exit status* in the README.
- **Breaking:** the representative is now the highest-scoring domain of the top hit, not
  its leftmost domain. Previously a short leading fragment could become the representative,
  giving wrong metadata or discarding a healthy family. Output changes for families whose
  top hit has several domains, and diverges from `reference/legacy_generate_families.py`.

### Deprecated

- `-n, --chunk_num` is renamed `-n, --chunk_id` in both subcommands, since the value is any
  string matching `[A-Za-z0-9._-]+`. `--chunk_num` still works and will be removed in 3.0.0.

### Fixed

- A failed output write could leave a family in both the generated and discarded outputs,
  leave a partial plus a duplicate `<chunk>_discarded.csv` row, or leave its HMM and
  alignments on disk beside a discard row. A failed family's partial files are now removed;
  if a shared per-chunk file cannot be written, or cleanup itself fails, the run exits `1`.
  Runs that write successfully are unaffected.
- Discarded families no longer hold their alignments in memory for the rest of the batch,
  lowering peak memory at large `--batch_size` and letting a run recover from one oversized
  family's `MemoryError`.

## [2.0.0] - 2026-07-29

### Changed

- **Breaking:** `<chunk>_metadata.csv` and `<chunk>_discarded.csv` now start with a header
  row, so they load with `pandas.read_csv` directly. Column order is unchanged. Readers that
  do not skip the header will treat it as data.

  ```
  family_id,full_msa_size,protein,region,length,sequence,consensus,converged
  representative,reason,value
  ```

- **Breaking:** an explicit `--fasta_index` is used exactly as given and never rebuilt, so
  one index can be shared safely by parallel chunks, including from read-only locations.
  Previously a copied, restored or symlinked index could be judged stale and rebuilt by
  every chunk at once, or fail with `PermissionError` where it could not be written. A path
  that does not exist is now rejected instead of being built. The default index under
  `--output_dir` is still built and refreshed automatically. Indexes built with
  `esl-sfetch --index` are interchangeable.

### Fixed

- Protein names containing underscores are no longer misread as slice bounds. Names such as
  `contig_1_gene_2_88_140` were truncated to `contig`, `contig_1_gene_x` discarded the whole
  family as an internal error, and `scaffold_12_34` was reported at invented coordinates.
  A trailing `_<start>_<end>` is now treated as bounds only when it matches the record's
  length. MGnifams accessions are unaffected. The legacy script has the same defect, so
  output for such names diverges from `reference/legacy_generate_families.py`.

## [1.0.0] - 2026-07-22

First release of `mgnifam` as a standalone package, ported from `generate_families.py` in
the [mgnifams](https://github.com/vagkaratzas/mgnifams) Nextflow pipeline. Outputs are
intentionally **not** byte-compatible with that script; the differences are listed below.

### Added

- `mgnifam` command with subcommands (`mgnifam generate_families ...`); `python -m mgnifam`
  is equivalent.
- `--output_dir` (default `output`) places every output under one root.
- `--fasta_index` reuses a pre-built Easel SSI index. Build one upstream and share it across
  chunks, or every chunk re-indexes the whole database.
- `--prefetch_targets` loads the database into RAM: faster, database-sized memory, identical
  output.
- `--batch_size` sets how many families are searched per `hmmsearch` pass. Defaults to
  `2 * cpus`; must be `>= cpus`.
- Leading and trailing envelope-overhang columns are trimmed from seed alignments.
- Families with `<= 2` sequences after redundancy trimming are discarded, instead of
  crashing `hmmbuild`.
- Inputs are validated before anything is written: chunk name characters, uncompressed
  FASTA, percentages in `[0, 1]`, ordered length bounds, well-formed cluster TSV rows.
- Bulk outputs are gzipped: `seed_msa/*.sto.gz`, `full_msa/*.sto.gz`, `hmm/*.hmm.gz`,
  `<chunk>_reps.fasta.gz`.
- A progress heartbeat is logged every 300 s during each database pass.

### Changed

- **Flatter output layout.** Per-family artifacts live in `seed_msa/`, `full_msa/`, `hmm/`
  and `rf/`; per-chunk files sit in the output root as `<chunk>_families.tsv`,
  `<chunk>_discarded.csv`, `<chunk>_successful.txt`, `<chunk>_converged.txt`,
  `<chunk>_metadata.csv`, `<chunk>_reps.fasta.gz` and `<chunk>.log`.
- **Much lower memory.** The FASTA is streamed and fetched by index rather than held in
  memory twice.
- **Byte-reproducible outputs.** HMMs no longer carry `DATE` or `COM` lines, and gzip
  members carry no timestamp or filename.
- `hmmalign` and FAMSA now respect `--cpus` instead of using every core on the node. Results
  are unchanged.
- Chunks sharing a working directory no longer collide on temporary files or the log.
- Re-running a chunk that now produces fewer families removes the surplus artifacts.
- Requires Python `>= 3.13`, `pyhmmer>=0.12.1,<0.13`, `pyfamsa>=0.7.0,<0.8`,
  `pytrimal>=0.8.5,<0.9` and `numpy>=2.5.1,<3`. `pandas` and `biopython` are no longer
  dependencies.

### Fixed

- **Hits failing `--recruit_evalue_cutoff` were recruited.** Families shrink accordingly —
  on the small fixture from 32/19/65 members to 31/19/61.
- **Results depended on `--cpus`.** Recruitment is now identical at any CPU count.
- **Runs were far slower than necessary.** Per-family work that scaled with database size ×
  cluster size, quadratic cluster selection, and one redundant database search per family
  are gone. These accounted for most of the legacy full-database run's ~8 months.
- **The exported HMM, seed MSA and RF line could come from a different round than the
  family's hits and full MSA** when a family ran out of rounds. All outputs now describe one
  round. On the `mgnifams_v2` fixture, families, convergence, members and full MSAs are
  unchanged; seed MSA, RF, HMM and consensus change for the 4 families that exhausted three
  rounds.
- **Discarded families were listed in `converged_families`.**
- **Every model was one match state short.** End-trimming dropped the last passing column,
  and trimmed a column even when no column passed. On the small fixture representatives
  lengthen by one residue and family `4497037939_1_144` now converges; membership is
  unchanged.
- **Sequence names in alignments could be truncated**, losing coordinates (e.g.
  `4454641265/11-11` for a region ending at 110), and `refined_families` rows carried
  trailing spaces.
- **Identical repeat domains in one protein collided**, silently dropping the second.
- **Stockholm alignments lacked `#=GF ID`** and carried hmmalign posterior-probability lines.
  Family IDs are restored and posteriors dropped; on the `mgnifams_v2` fixture no membership
  changes.

## [0.1.0] - legacy

The original `mgnifams/bin/generate_families.py`, kept unmodified at
`reference/legacy_generate_families.py` as the behavioural reference.

- Builds an HMM per cluster, recruits from the database with `hmmsearch`, re-aligns with
  `hmmalign`, removes redundancy with `pytrimal`, and stops on convergence or after three
  rounds.
- Discards clusters on representative length, seed membership, or a low-complexity model.
- Emits seed and full Stockholm alignments, an HMM, an RF line, family metadata and a
  representative FASTA per chunk.
- Ran on `python=3.13.5`, `pyhmmer=0.11.1`, `numpy=2.3.2`, `pandas=2.3.2`,
  `biopython==1.85`, `pyfamsa==0.6.0`, `pytrimal==0.8.2`, holding the whole FASTA in memory
  twice.
