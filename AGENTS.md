# AGENTS.md

Notes for coding agents working in this repository. The README explains what the tool
does and how to run it; this file covers what will bite you. Read both.

## Ground rules

- **`reference/` is immutable.** `reference/legacy_generate_families.py` is the
  behavioural oracle for the whole port, vendored byte-identically from the `mgnifams`
  pipeline. It is excluded from every `pre-commit` hook. Never format, lint, or "fix" it.
  When you need to know what the code *should* do, read it.
- **Never add a dependency.** The four scientific packages are upper-bounded on purpose:
  they decide hit retention, alignment, and serialised bytes, so the reproducibility
  contract holds only for the resolved set in `uv.lock`. Use `uv run ...` for everything.
- Use `uv lock --check && uv sync --frozen` before trusting a test result.
- Commits and pushes are the human's call.

## Before you change anything in `generate_families.py`

Five invariants are load-bearing. Each was a real bug at some point, each is guarded by
exactly one test, and each survives a passing end-to-end run if you break it.

1. **`hmmsearch` must be called with `parallel="queries"`, explicitly.** Left to choose,
   pyhmmer picks `parallel="targets"` whenever the query count is below `cpus` — which
   the last wave of a batch usually is. Target-parallelism merges each worker's stored
   hits while re-thresholding only the reporting flags, so results start depending on the
   machine's core count. All searching goes through `search()`; keep it that way.

2. **Extraction reads `top_hits.reported`, never `top_hits` directly.** The raw list
   holds hits that failed `--recruit_evalue_cutoff`. Iterating it is what the legacy code
   did, and it is the bug 1.0.0 fixed. Guard:
   `test_unreported_hit_is_not_recruited_and_prefetch_matches` — query `4497037939_1_144`
   stores 55 hits, reports 54, and the extra one is sequence `6320430079`.

3. **Hit records are plain tuples, copied out immediately.** A pyhmmer `Domain` holds a
   reference to its `Hit`, which holds the whole `TopHits`. Caching `Domain` objects pins
   every result graph for the batch. Guard:
   `test_extracted_records_do_not_retain_pyhmmer_results`.

4. **Reported domains are sorted by score within each hit.** HMMER reports domains in
   positional order, but row 0 supplies representative metadata and length. The
   top-ranked hit's highest-scoring domain must lead. Guard:
   `test_extract_records_puts_top_scoring_domain_first`.

5. **The SSI index guard is an `if ... raise`, not an `assert`.** `python -O` strips
   asserts, and the failure it prevents is a silent wrong-sequence read from a stale
   index. Guard: `test_index_mismatch_guard_and_optimized_python`.

Also: `Family.advance()` and `Family.finish()` must stay side-effect free. Family ids are
the rank among *successful* families in cluster-file order, so nothing may be written
until a family's fate is known. All writing happens in `emit_family()`, called in cluster
order. A family that converges in round 1 must not record itself ahead of an earlier
family still running in round 3. Discarded families must never appear in
`converged_families`. Every family must leave `emit_family()` on exactly one side of the
generated/discarded split; an internal crash is recorded as a discard so the chunk can
finish, then also fails the completed run with exit 3. Exit 1 means containment could not
leave coherent output, so that output must not be consumed.

## Deliberate differences from the legacy script

The legacy oracle records the starting behaviour, not every current guarantee.
`CHANGELOG.md` documents the intentional fixes; the historical PLAN files describe
earlier decisions and must not be used to reverse them.

- On convergence and at `MAX_ROUNDS`, `Family.advance()` skips the re-align/trim tail.
  The seed remains the one that built the final search model. Legacy ran an unused
  final trim and exported a model built from a seed it had never searched with.
  Guard: `test_final_round_leaves_the_searched_seed_in_place`.
- `renumber_msa()` names both output alignments, so seed and full Stockholm files
  contain `#=GF ID <family>`. HMMER posterior annotations are omitted, but the RF line
  remains. Legacy's `renumber_sto_msa` stripped the family ID along with those annotations.
  Guard: `test_declared_outputs_parse_and_long_fixture_runs`.

There are two distinct clipping functions. `clip_env_ends()` reads the `#=GC RF` line and
trims envelope overhangs from seed alignments. `clip_ends()` reads gap occupancy and runs
after redundancy trimming. They are not interchangeable.

Legacy `clip_ends()` dropped the last column that passed the
occupancy threshold, and reported a full span when no column passed. Both were fixed in
1.0.0, so it now diverges from `reference/legacy_generate_families.py` on purpose. If you
diff against the legacy baseline, expect every model to be one match state wider.

`split_slice_name()` is another deliberate divergence. Legacy recovered a slice's
parent protein with `split("_")` on exactly three fields, which truncated any protein
name that itself contains underscores, raised on non-numeric trailing fields, and
invented coordinates for a name like `scaffold_12_34`. It now splits from the right and
accepts the trailing two fields as bounds only if they span the record exactly. On a
database of bare MGnifams integer accessions the two agree on every record — the span
test holds for every real slice — so this changes no output for the reference data.

## Testing

`generate_families` writes every generated artifact under `--output_dir` (default:
`output/`), including an automatically built SSI index.

`uv run pytest` runs the suite against real 50,000- and
26,949-sequence fixtures rather than toy data, because the marginal hits that several
tests depend on only exist at that scale.

`tests/fixtures/*.fa.gz` are decompressed into `tmp_path` by `conftest.py`. Easel cannot
seek inside a gzip stream, so an SSI index cannot be built over a compressed FASTA and the
CLI rejects one.

### Trying `update_families` by hand

Two fixtures exist for this and are not used by the automated suite:

- `mgnifams_v2.hmm.lib.gz` — the 14 families `generate_families` builds from
  `mgnifams_v2.tsv`, concatenated into one multi-model library. **Derived**, so regenerate
  it if `generate_families`' output ever changes:

  ```bash
  uv run mgnifam generate_families -c tests/fixtures/mgnifams_v2.tsv -f v2.fa -n v2 --output_dir gen
  for f in gen/hmm/*.hmm.gz; do zcat "$f"; done | gzip -9n > tests/fixtures/mgnifams_v2.hmm.lib.gz
  ```

- `mgnifams_v3_additions.fa` — five sequences to append to `mgnifams_v2.fa`, standing in for
  a new release. Four are real family members carrying 12–20% conservative substitutions,
  two aimed at `v2_10` and two at `v2_6`; the fifth is `v2_10`'s representative shuffled,
  which nothing should recruit. Each header says which. Kept as a separate file rather than
  a second full FASTA: it is 1 KB against 3.4 MB, and you can read it.

```bash
zcat tests/fixtures/mgnifams_v2.fa.gz > v3.fa
cat tests/fixtures/mgnifams_v3_additions.fa >> v3.fa
uv run mgnifam update_families -i tests/fixtures/mgnifams_v2.hmm.lib.gz -f v3.fa \
    -n demo --skip_refine --output_dir out
```

What to expect, and why these numbers are the point of the fixture:

| | `--skip_refine` | refine |
|---|---|---|
| `v2_10` | 4 → 6 members | 4 → 6, model 180 → 179, converges in 2 rounds |
| `v2_6` | 14 → 16 members | 14 → 16, model 125 → 124, converges in 2 rounds |
| the other 12 | unchanged | `v2_4` is **discarded** |
| decoy `9000000005` | recruited by nothing | recruited by nothing |

`v2_4` discarding under refine has nothing to do with the added sequences — none of them
are recruited by it, and it discards identically when refined against the **unchanged**
`mgnifams_v2.fa`. It is useful anyway, as the only easy way to see a discard row and the
nullable delta columns (`retention` and `full_msa_size` are empty, because it never reached
the stage that computes them).

**What it actually demonstrates: refine mode is not idempotent.** `generate_families` stops
at convergence or `MAX_ROUNDS` and deliberately skips `advance`'s re-align/trim tail on the
way out, so a finished family's seed is the one that built its final model. Refining that
model runs the tail again — it continues the iteration rather than repeating it, and the
family settles somewhere `generate_families` never took it.

That settling is a contraction, not a decay. Refining the 14 v2 families against their own
unchanged database repeatedly, every survivor reaches an exact fixed point within about ten
iterations and then does not move for at least another twenty-four. Models are not
monotonically shrinking either: `hmmalign` adds insert columns *before* the trims run, and
`hmmbuild`'s `architecture="fast"` promotes a column to a match state on occupancy, so a
model **grows** when the new recruits genuinely support it — measured at 180 → 193 for a
14-residue insertion carried by a diverse subfamily that outnumbers the original members.
`--max_seq_identity` is what keeps that honest: near-identical recruits collapse to one
representative in the seed, so a burst of duplicates cannot inflate a model, only real
diversity can.

The cost is paid once, while a family settles, and it falls on the marginal ones. Of the 14,
one died in the first two iterations:

- `v2_4` eroded — 78 → 64, representative 66, under the default `--discard_min_rep_length 75`.
  It had three residues of headroom to begin with.

The previous refine run also exposed issue #8 in `generate_families`: `v2_9`'s top-ranked
hit carried a 17-residue leftmost domain and a 99-residue higher-scoring domain. HMMER
reports those domains positionally, so the short fragment became row 0 and the healthy
family was discarded. `extract_records` now sorts each hit's reported domains by score,
keeping the 99-residue domain as the representative.

So `--skip_refine` is the right default for a release refresh: it recruits from the new
database and leaves the models exactly as they were. Reach for refine when you want the
models re-derived, and expect a one-off cull of marginal families as they settle.

When you change scientific behaviour, the honest check is a diff against the legacy
script, not a green test suite. Reproduce the baseline with a pinned environment:

```bash
uv venv -p 3.13 /tmp/legacy && VIRTUAL_ENV=/tmp/legacy uv pip install \
  "pyhmmer==0.11.1" "numpy==2.3.2" "pandas==2.3.2" "pyfamsa==0.6.0" "pytrimal==0.8.2"
/tmp/legacy/bin/python reference/legacy_generate_families.py -c ... -f ... -p 1 -n test ...
```

Run the legacy script at `-p 1`. At higher core counts its own output is not stable.

## Adding a subcommand

`mgnifam` dispatches through `src/mgnifam/cli.py`. Add a module exposing `main(argv)`, add
it to `COMMANDS`, and set its parser's `prog` to `"mgnifam <name>"`. Subcommands own their
parsers rather than registering with `add_subparsers`, so each stays directly callable and
testable without the dispatcher. `remove_redundant` and `merge_families` are expected here.

`update_families` imports the algorithm out of `generate_families` and owns only its own
`main()` loop. **When a third command lands, extract the shared machinery into a
`_pipeline.py` and make all three thin.** Not before: the shape three commands share is
worth observing rather than guessing, and `README.md` documents
`from mgnifam.generate_families import build_ssi_index` as a public entry point that a move
would break.

## Before you change anything in `update_families.py`

- **Reject input/output aliases before writing or retry cleanup.** Single-file libraries,
  model directories, FASTA and supplied SSI files must survive rejection unchanged.
  Check resolved paths and existing file identities, including output-folder symlinks.
  Guards: `test_input_models_cannot_be_overwritten`,
  `test_input_library_cannot_alias_an_aggregate`.

- **A family's identity is its model's `NAME`, and it must reach every output.** Filenames
  are the easy half. `<chunk>_updated_metadata.csv`, `<chunk>_updated_families.tsv`,
  `converged.txt` and the `reps.fasta.gz` annotation all carry it too, and the annotation is
  where a `f"{chunk}_{id}"` reconstruction hid until a review caught it: with `--chunk_id 9`
  and a preserved `1_7` it wrote `9_1_7`. Guard:
  `test_identity_is_preserved_in_every_field_not_only_in_filenames`, which deliberately runs
  a chunk number unrelated to the names.
- **Family names are untrusted input.** They come from a third-party HMM and are
  interpolated into artifact paths and unquoted CSV fields. `[A-Za-z0-9._-]+` plus rejecting
  `.`/`..` closes traversal and CSV-column corruption together. Guard:
  `test_unsafe_family_names_are_rejected_before_any_output_exists`.
- **Models are sorted by name, for a directory and a library alike.** Filename order and
  library order are different orders over the same models, and the aggregates would differ
  byte-for-byte. Guard:
  `test_library_and_directory_inputs_agree_despite_adversarial_ordering`.
- **Refuse foreign families, then clear accepted names before a retry.** An updated
  family keeps its model's name, so `--chunk_id` cannot establish artifact ownership.
  A smaller input set over a larger output set is refused before deleting anything.
  After collision validation, clear all four exact artifact paths for each input name:
  discarded families must leave no old models, and recruit-only retries must leave no
  old seed/RF files. Cleanup errors abort before aggregate files are opened. Guards:
  `test_a_directory_holding_another_runs_families_is_refused`,
  `test_a_partial_run_can_be_rerun_in_place`, `test_retry_artifacts_match_the_new_outcome`,
  `test_failed_retry_cleanup_aborts_before_aggregate_truncation`.
- **Delta metrics are captured where they are computed.** `discard()` clears the records,
  model and alignments, and `finish()` holds the membership fraction in a local before it may
  discard on representative length. Reading them back off the `Family` afterwards silently
  yields empty columns.
- Round 1 cannot converge, which is why refine mode needs no stored seed MSA. If that ever
  changes, `emit_family`'s hand build has nothing to build from on a converged-at-round-1
  family.

## Docs that must move together

`README.md`, `CHANGELOG.md`, `AGENTS.md` and the docstrings all describe the same
guarantees. `PLAN.md` and `PLAN-REVIEW-LOG.md` are a historical record of how the port was
designed and reviewed — read them for *why*, do not update them to reflect new work.
