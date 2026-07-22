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
`converged_families`.

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

## Testing

`generate_families` writes every generated artifact under `--output_dir` (default:
`output/`), including an automatically built SSI index.

`uv run pytest` — 37 tests, roughly half a minute. They run against real 50,000- and
26,949-sequence fixtures rather than toy data, because the marginal hits that several
tests depend on only exist at that scale.

`tests/fixtures/*.fa.gz` are decompressed into `tmp_path` by `conftest.py`. Easel cannot
seek inside a gzip stream, so an SSI index cannot be built over a compressed FASTA and the
CLI rejects one.

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

## Docs that must move together

`README.md`, `CHANGELOG.md`, `AGENTS.md` and the docstrings all describe the same
guarantees. `PLAN.md` and `PLAN-REVIEW-LOG.md` are a historical record of how the port was
designed and reviewed — read them for *why*, do not update them to reflect new work.
