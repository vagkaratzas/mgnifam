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

Four invariants are load-bearing. Each was a real bug at some point, each is guarded by
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

4. **The SSI index guard is an `if ... raise`, not an `assert`.** `python -O` strips
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

## Behaviour that looks like a bug and is not

Do not "fix" these. They are reproduced from the legacy script on purpose, are documented
in `CHANGELOG.md` under *Preserved deliberately*, and each carries a comment at the site:

- The `family_iteration > 3` path writes an HMM that is not the model used for the final
  search and alignment.
- `renumber_sto_msa` strips every `#=GF`/`#=GS`/`#=GR` line, so the seed Stockholm has no
  `#=GF ID`. The family name reaches the output as the HMM's `NAME` field. Asserting
  `#=GF ID` in a seed `.sto` is wrong; a plan once did, and only running the code caught it.

There are two distinct clipping functions. `clip_env_ends()` reads the `#=GC RF` line and
trims envelope overhangs from seed alignments. `clip_ends()` reads gap occupancy and runs
after redundancy trimming. They are not interchangeable.

`clip_ends()` used to be on the list above: it dropped the last column that passed the
occupancy threshold, and reported a full span when no column passed. Both were fixed in
1.0.0, so it now diverges from `reference/legacy_generate_families.py` on purpose. If you
diff against the legacy baseline, expect every model to be one match state wider.

`split_slice_name()` is the second deliberate divergence. Legacy recovered a slice's
parent protein with `split("_")` on exactly three fields, which truncated any protein
name that itself contains underscores, raised on non-numeric trailing fields, and
invented coordinates for a name like `scaffold_12_34`. It now splits from the right and
accepts the trailing two fields as bounds only if they span the record exactly. On a
database of bare MGnifams integer accessions the two agree on every record — the span
test holds for every real slice — so this changes no output for the reference data.

## Testing

`generate_families` writes every generated artifact under `--output_dir` (default:
`output/`), including an automatically built SSI index.

`uv run pytest` — 37 tests, roughly half a minute. They run against real 50,000- and
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
model runs the tail again — it continues the iteration rather than repeating it. Both trims
in that tail (`run_pytrimal_reps`, then `clip_ends`) only ever remove columns, so models
shrink monotonically. Refining the 14 v2 families against their own unchanged database
shrinks 10 of them and grows none; `v2_4` had three residues of headroom over the 75 floor,
lost 14, and died.

So `--skip_refine` is the right default for a release refresh: it recruits from the new
database and leaves the models exactly as they were. Reach for refine when you actually
want the models re-derived, and expect marginal families to be discarded by the erosion
rather than by anything in the new data.

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

- **A family's identity is its model's `NAME`, and it must reach every output.** Filenames
  are the easy half. `<chunk>_updated_metadata.csv`, `<chunk>_updated_families.tsv`,
  `converged.txt` and the `reps.fasta.gz` annotation all carry it too, and the annotation is
  where a `f"{chunk}_{id}"` reconstruction hid until a review caught it: with `--chunk_num 9`
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
- **An output directory is not cleared, it is refused.** `generate_families` clears its
  own past output from a `<chunk>_<rank>` regex, which works because it derives names. An
  updated family keeps its model's name and `--chunk_num` never enters a per-family
  filename, so the owned set cannot be derived from the directory — only recorded, and a
  record that must survive between runs and be replaced atomically is a lot of machinery
  for one case. Re-running the same models in place is allowed, so a failed chunk retries;
  a smaller set over a larger one is refused. Guards:
  `test_a_directory_holding_another_runs_families_is_refused`,
  `test_a_partial_run_can_be_rerun_in_place`.
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
