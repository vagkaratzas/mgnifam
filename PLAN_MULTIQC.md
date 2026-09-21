# MultiQC support for mgnifam

## Context

The goal is to make `mgnifam generate_families` and `mgnifam update_families` reportable in
MultiQC, and to get mgnifam into MultiQC's supported-tools list
(https://docs.seqera.io/multiqc/modules/).

**`*_mqc.*` files alone won't get us listed.** Files named `*_mqc.(tsv|csv|yaml|json|…)` are
MultiQC *custom content*. They render in any report, but the docs say explicitly that custom
content never appears in the modules list. A listing needs a **native module** merged into
`MultiQC/MultiQC`. Its requirements are:
- `multiqc/modules/<tool>/` containing the module code, a docstring, a license and a **DOI**
  (enforced by the CI linter)
- a `search_patterns.yaml` entry, a `pyproject.toml` entry point and a pytest suite
- example files in the separate `MultiQC/test-data` repo, which must be in place "before the PR will go any further"
- module requests are triaged by "tool popularity, community need, and request quality". A PR
  that arrives complete, with example files, is the strongest route for a niche tool.

MultiQC's design rule is that modules *parse* tool output and don't compute metrics
("summarise, don't replicate"). So the right mgnifam-side contribution is one small, stable,
machine-readable summary per run for the module to parse. `_mqc` files are the wrong target
(see non-goals).

**User decisions (settled):**
- **Output:** a stats JSON only, with no `_mqc` files.
- **DOI:** minted through the Zenodo–GitHub integration on the next release. The user does this.
- **Scope:** implement the mgnifam side in this repo, and draft the MultiQC module and
  test-data files in separate fork checkouts, ready to submit. **The user opens every
  external PR and issue.**

**One file per invocation = one file per chunk.** One `mgnifam` invocation processes exactly
one chunk. The mgnifams pipeline runs many chunks, and MultiQC merges their files into one
report with one row (sample) per chunk. Aggregating "per pipeline run" is MultiQC's job, and
the `<chunk>` prefix keeps chunks sharing an output root, or staged together for MultiQC,
from colliding.

## Goal and acceptance criteria

1. Every **completed** run (exit 0 or 3) writes `<output_dir>/<chunk>_stats.json`
   (`generate_families`) or `<output_dir>/<chunk>_updated_stats.json` (`update_families`).
2. A run that exits 1 (corrupted output), or dies before completing, leaves **no** stats
   file. That includes a stale one from an earlier run in the same directory, so a present
   stats file always means consumable output.
3. The file is byte-deterministic under the existing contract: identical across `--cpus`,
   `--batch_size`, `--prefetch_targets` and `PYTHONHASHSEED`. The existing reproducibility
   tests compare every file except `.log`/`.ssi`, so they cover it automatically.
4. Every count in the file equals what the aggregate CSVs of the same run say.
5. A draft MultiQC module parses both files and renders General Stats plus per-command
   sections. `multiqc --strict` passes against the drafted test data, and the module's
   pytest passes in the fork.
6. The five invariants in `generate_families`, and the rule that `emit_family` does all
   the writing, are untouched.

## Approach: the mgnifam side (this repo)

### Stats file contents (schema v1)

```json
{
  "tool": "mgnifam",
  "schema_version": 1,
  "version": "3.1.0",
  "command": "generate_families",
  "chunk_id": "1",
  "exit_status": 0,
  "parameters": { "discard_min_rep_length": 75, "...": "every scientific threshold; plus skip_refine for update" },
  "families": { "input": 14, "successful": 12, "discarded": 2, "converged": 9, "crashed": 0 },
  "discard_reasons": { "few seed sequences remained": 1, "...": 1 },
  "histograms": {
    "full_msa_size": { "4": 3, "16": 1 },
    "model_length": { "124": 1 },
    "representative_length": { "117": 2 }
  }
}
```

`update_families` adds `histograms.model_length_change` (after − before, integer keys),
`histograms.rounds_run`, and `histograms.retention`. Retention keys are the **exact**
strings `delta.csv` already holds (`str(float)`, the shortest round-trip repr). That makes
them lossless and deterministic, and the module bins them only for display. Families with
an empty retention (for example no hits) are left out of that map.
- `"tool": "mgnifam"` comes first so the MultiQC search pattern can match on contents
  within 3 lines.
- `schema_version` lets the out-of-tree module accept v1 and skip anything else with a
  warning.
- `parameters` excludes paths, `cpus`, `batch_size` and `prefetch_targets`. They don't
  change results, and including them would break acceptance criterion 3.
- `version` is `mgnifam.__version__`, never a literal. The `3.1.0` above is only
  illustrative.
- A zero-cluster chunk is valid, since an empty TSV is accepted. It yields all-zero counts
  and empty maps, and the module must render 0/0 as an empty (`null`) percentage rather
  than dividing by zero or dropping the sample.
- Histogram maps are `{value: count}`, sorted by numeric key. That is lossless and compact.
  The module does the binning at display time.
- Serialised with `json.dumps(..., indent=2)` plus `"\n"`. Key order comes from
  construction, not from `sort_keys`, so that `tool` stays first.

### Where the numbers come from: read the run's own aggregates back

**Don't touch `emit_family`, `Family` or the round loop.** Once the aggregate writers'
`ExitStack` has closed, a small function re-reads the chunk's own files:
- `generate_families`: `<chunk>_metadata.csv` gives `full_msa_size`, `length`, `converged`,
  and `len(consensus)` as the model length. `<chunk>_discarded.csv` gives the reasons.
- `update_families`: `<chunk>_updated_delta.csv` holds every family's outcome, its lengths
  before and after, `retention` and `rounds_run`. Its metadata CSV gives `full_msa_size` and
  representative length.

Why read back instead of adding counters to the loop:
- The summary cannot disagree with the CSVs (acceptance criterion 4).
- The final model length for `generate_families` only exists inside `emit_family`
  (it's a hand build), and the consensus column already carries it.
- It adds zero new state across the containment boundaries that AGENTS.md protects.

The cost is one extra read of chunk-sized CSVs, which is negligible next to `hmmsearch`.
`input` = successful + discarded rows, which equals the cluster/model count.
`converged` means **successful families that converged**, for both commands: generate's
metadata rows with `converged=True`, and update's delta rows with `outcome == successful`
and `converged == True`. A delta row can carry `converged=True` on a discard, because
`discard()` preserves `ever_converged`, so update's `converged` must never count raw delta
`converged` values. The
`crashed` count reuses the main loop's existing `crashed` variable, and it sets
`exit_status`.

### Lifecycle

**Invariant:** a stats file is present exactly when the directory holds the coherent,
completed output of the run that wrote it.

- **Remove any stale stats file as the first mutation of the output tree.** That is the
  first *mutation*, not the first statement, so every read-only refusal runs before it.
  - `generate_families`: at the top of `prepare_output_directories()`, which has no refusal
    checks.
  - `update_families`: inside `prepare_output_directories()`, *after* the foreign-family
    `strays` refusal and before the `mkdir`/artifact unlinks.
  - In both cases it comes before any per-family artifact is cleared or any aggregate is
    truncated.
  - A run rejected earlier (usage error, validation failure, refused retry) has touched
    nothing. The old output and old stats file are then still coherent, and the non-zero
    exit already tells the caller that this invocation produced nothing.
  - Removing the file earlier would destroy a valid report for no gain.
- **Write it last:** after the `DONE.` log line and before `raise
  SystemExit(EXIT_CRASHED_FAMILIES)`, inside the `try`, so a `ChunkCorrupted` path never
  reaches it.
  - **Write atomically:** `tempfile.NamedTemporaryFile(dir=root, delete=False)` (created
    exclusively under a random name, so it can never alias an existing input), then write,
    `flush`, `os.fsync`, close and `os.replace` it onto the final path.
    - On any failure, unlink the temp file and raise `ChunkCorrupted`, which gives exit 1
      with no stats file. That is consistent with the invariant.
    - A SIGKILL can leave only an orphan random-named `tmp*` file, never a partial
      `<chunk>_stats.json`.
- **Protect inputs from the stats path.**
  - `update_families.validate_output_paths()` adds `f"{prefix}_stats.json"` to its
    destination list.
  - `generate_families` has no output-alias check today. Add a minimal one just before the
    stale unlink: raise `ValueError` if the stats path resolves to, or shares an inode with,
    `--clusters_chunk`, `--fasta_file` or `--fasta_index`. This is needed because
    `os.replace` would overwrite an input sitting at the final path.
  - The same gap already exists for generate's other aggregates. It is pre-existing and out
    of scope, noted under risks.

### Code placement (the smallest diff)

- `src/mgnifam/generate_families.py`: one `write_stats(path, payload)` helper
  (`json.dumps` + `write_text`), a small `guard_stats_path(path, inputs)` alias check, and one `generate_stats(root, chunk, options, crashed)` that builds
  the dict from the CSVs. `stats_parameters(options, names)` picks the threshold flags.
  Then the stale unlink and the call site in `main()`.
- `src/mgnifam/update_families.py`: `update_stats(...)` reuses `write_stats` and
  `stats_parameters` (it already imports its machinery from `generate_families`, per
  AGENTS.md). It also adds the stale unlink, the alias-list entry and the call site.
- No new module and no new dependency. Only `json`, `csv`, `collections.Counter` and `os`
  from the stdlib. The `_pipeline.py` extraction is still deferred, because AGENTS.md
  triggers it on a third command, not on this change.

### Tests (added to the existing files, using the fixtures already there)

- `test_generate_families.py::test_stats_file_matches_aggregates`: after a clean run on
  the small fixture, assert the **whole payload**.
  - Recompute `families`, `discard_reasons` and every histogram in the test with `Counter`
    over the CSV columns, and compare them to the JSON.
  - Assert the exact `parameters` key set (no `cpus`/`batch_size`/`prefetch_targets`/paths),
    `version == mgnifam.__version__`, `chunk_id`, `exit_status == 0`, and that `tool` is
    the first key.
- Empty cluster TSV: the run completes and the stats file has zero counts and empty maps.
- Generate alias guard: `--fasta_file` placed at `<output>/<chunk>_stats.json` raises and
  is byte-unchanged.
- Extend `test_shared_append_failure_is_fatal_not_a_contained_discard`: seed a stale
  `chunk_stats.json` first, then assert it is **absent** after exit 1.
- Extend `test_console_script_returns_three_after_a_contained_family_crash`: the stats file
  exists with `exit_status == 3` and `crashed == 1`.
- `test_update_families.py`: one test for both modes (parametrised on `skip_refine`) on the
  v2 library. It compares **every** top-level count and histogram to `Counter` projections:
  `families`, `discard_reasons`, `model_length_change`, `rounds_run` and `retention` from
  `delta.csv`, and `full_msa_size`, `model_length` and `representative_length` from the
  metadata CSV. It also checks that `parameters.skip_refine` is correct.
  - Refine mode on v2 covers the discard path (`v2_4`), and `--skip_refine` covers the
    all-success path.
- An update where no model has any hits (the v2 library against a FASTA that holds only the
  shuffled decoy record from `mgnifams_v3_additions.fa`): every outcome is
  `no hits in the new database`.
  - The `delta.csv`-derived maps are **not** empty. For N families, `model_length_change`
    is `{0: N}` and `rounds_run` is `{1: N}`, because `qlen` and `rounds_run` are set
    before the no-hit discard.
  - `retention` and the metadata-derived maps (`full_msa_size`, `model_length`,
    `representative_length`) are empty.
- Stats write failure: run with a prebuilt explicit `--fasta_index` (the existing
  `shared_index` pattern), so the SSI build's own `os.replace` is never reached. Then
  monkeypatch `os.replace` in `mgnifam.generate_families` to raise. Expect exit 1, no stats
  file, no leftover `tmp*` file in the root, and a log that says the failure happened in
  the stats commit.
- Converged-then-discarded regression: call the update stats builder directly on a
  hand-written `delta.csv` containing a row with `converged=True` and a discard outcome.
  That row must not count toward `families.converged`.
- Update exit paths: extend the existing update contained-crash (exit 3) test to assert
  `<chunk>_updated_stats.json` has `exit_status == 3` and `crashed == 1`. Extend the update
  fatal-output (exit 1) test to seed a stale stats file and assert it is absent afterwards.
- Update foreign-family refusal (extend `test_a_directory_holding_another_runs_families_is_refused`):
  a pre-existing stats file survives the refusal byte-unchanged. Extend `test_input_library_cannot_alias_an_aggregate`
  (or add a case) for the stats path.
- The existing cross-`--cpus`/`--batch_size`/prefetch byte-comparison tests cover
  determinism unchanged.

### Docs (these must move together, per AGENTS.md)

- **README:** add the file to both Outputs tables, a short "MultiQC" section (schema,
  presence = complete run, one per chunk), and note that exit 1 leaves no stats file.
- **CHANGELOG:** add an `## [Unreleased]` → *Added* entry for the stats file.
- **AGENTS.md:** add a short rule. The stats file is derived by reading back the aggregates,
  is written last and atomically, is removed at start, and excludes non-scientific flags.
  Bump `schema_version` on any breaking change, because a module in another repo depends
  on it.
- Leave `PLAN.md` / `PLAN-REVIEW-LOG.md` alone; they're historical.

## Approach: the MultiQC side (drafted locally, submitted by the user)

Work in sibling checkouts outside this repo: `../MultiQC` (a fork of `MultiQC/MultiQC`) and
`../test-data` (a fork of `MultiQC/test-data`). Following `docs/markdown/development/modules.md`
of the checked-out version:

- `multiqc/modules/mgnifam/__init__.py` and `mgnifam.py` hold
  `MultiqcModule(BaseMultiqcModule)` with name, anchor, href, info, DOI (the Zenodo concept
  DOI, left as a clearly marked TODO until it is minted) and a docstring.
  - One search-pattern key, `mgnifam`: `fn: "*_stats.json"`, `contents: '"tool": "mgnifam"'`,
    `num_lines: 3`.
  - Sections are selected by the JSON `command` field, **never by filename shape**, and each
    command gets its own data dict.
    - The sample name is the filename minus `_stats.json` (`1`, `1_updated`), passed through
      `clean_s_name`.
    - A name can still collide: generate with `--chunk_id 1_updated` and update with
      `--chunk_id 1` produce the same filename. That ambiguity already exists for every
      aggregate (`1_updated_metadata.csv`), so it isn't solved here. The README note
      "one output root per run" covers it.
  - A zero-input sample renders with an empty percentage.
  - Mean retention is computed only over the retention map. When the map is empty (all
    no-hit or all early discards), the General Stats cell is null/omitted, never 0.
  - Files whose `schema_version != 1` are skipped with a log warning. Call
    `add_software_version(version, s_name)` and `add_data_source`, and honour
    `ignore_samples`. Raise `ModuleNoSamplesFound` when nothing parses.
  - **General Stats** columns: families in, % successful, converged, and mean retention for
    update runs only.
  - **Sections:**
    1. An outcome bar graph per chunk: successful stacked with each discard reason.
    2. Line graphs of `full_msa_size` and `model_length`.
    3. Update runs only: line graphs of `model_length_change` and `retention`.
- Register it in `pyproject.toml` entry points and `search_patterns.yaml`, and add it to
  `config_defaults.yaml` `module_order` if the checked-out docs still require that.
- Write `multiqc/modules/mgnifam/tests/test_mgnifam.py` in the style of an existing JSON
  module's tests.
- For test-data, add `data/modules/mgnifam/` with `v2_stats.json` (from `generate_families`
  on `mgnifams_v2.tsv`) and `demo_updated_stats.json` (from `update_families --skip_refine`
  on the v3 demo in AGENTS.md). Both are real outputs, produced by this repo's build.
- Also draft a module-request issue text (tool, repo, one-line description, key General
  Stats values), which the user can file or skip.

**External steps (the user does these, in order):**
1. Merge and release mgnifam 3.1.0 with Zenodo enabled.
2. Put the DOI into the module.
3. Open the test-data PR.
4. Open the MultiQC PR with the `module: new` label.

## Non-goals

- **`_mqc` custom-content files.** They would duplicate sections once the native module is
  merged, and they clutter every chunk's output. For the gap before the module lands, the
  mgnifams Nextflow pipeline can convert the JSON if needed.
- **`_mqc_versions.yml`.** That is a pipeline-level convention. The module reports the
  version from the JSON instead.
- **A MultiQC plugin shipped inside `mgnifam`.** It would add a `multiqc` dependency, which
  AGENTS.md forbids.
- **Parsing `.log` files**, which carry timestamps, or the per-family CSVs directly
  (large, several files per sample, no completion or version signal).
- **`CITATION.cff` / `.zenodo.json`.** Zenodo works without them. Add them only if the
  user wants rich citation metadata.
- **A version bump or release.** Those are the user's call. The CHANGELOG entry goes under
  Unreleased.

## Assumptions and risks

| Assumption / risk | Source | Mitigation |
|---|---|---|
| Custom content is never listed as a module | docs.seqera.io/multiqc/custom_content | none needed |
| DOI and license are required, and the linter enforces it | docs.seqera.io/multiqc/development/modules | Zenodo DOI before the PR |
| Test data must land first | MultiQC `.github/CONTRIBUTING.md` | draft the test-data PR alongside the module |
| Requests are triaged by popularity, so a niche tool may be deprioritised or declined | CONTRIBUTING.md, module-request.yml | submit a complete PR with example files and a DOI; if declined, the JSON still works through the pipeline or a plugin |
| The MultiQC module API drifts between versions | not verified against HEAD | follow the fork's checked-out `modules.md` and existing modules at build time |
| Generate's other aggregates can alias an input (pre-existing, no check) | `generate_families.main` | out of scope for MultiQC; only the new stats path is guarded. **Follow-up:** open a separate issue for a full generate output-alias check mirroring `update_families.validate_output_paths` |
| Read-back can't see a partially written CSV | CSVs closed by the `ExitStack` before the read | call only after the `with` block exits, on the completed path |

## Verification

In this repo:
```bash
uv lock --check && uv sync --frozen
uv run pre-commit run --all-files
uv run pytest
# manual: stats present, plausible, deterministic across cpus
uv run mgnifam generate_families -c tests/fixtures/mgnifams_v2.tsv -f v2.fa -n v2 --output_dir gen -p 1
uv run mgnifam generate_families -c tests/fixtures/mgnifams_v2.tsv -f v2.fa -n v2 --output_dir gen4 -p 4
cmp gen/v2_stats.json gen4/v2_stats.json && python -m json.tool gen/v2_stats.json
```
Expected: all tests green. For v2, `families` is `input 14, successful 12, discarded 2,
converged 8`, per `tests/test_generate_families.py` (the v2 end-to-end test). `cmp` is
silent.

In the MultiQC fork (in its own venv, not this repo's):
```bash
pip install -e . && pytest multiqc/modules/mgnifam
multiqc --strict ../test-data/data/modules/mgnifam -o /tmp/mqc-out
```
Expected: no strict-mode errors; the report shows the mgnifam General Stats columns and sections
for both `v2` and `demo_updated`. Manual visual check of the HTML report.

## Loop setup

- Host: Claude (Opus 5) plans and builds.
- Plan reviewer: Codex `gpt-5.6-luna` at `xhigh`, as the user requested.
- Final inspector: a fresh Codex session.
- Rounds: at most 5 for plan review, 2 for fixes and 2 for inspection.
- Authorization: **planning only**. The build needs explicit approval.
- The review log is kept in this file's appendix, because the repo's `PLAN-REVIEW-LOG.md`
  is historical and must not be edited.
- **First step after plan approval:** copy this plan to `<repo>/PLAN_MULTIQC.md`, and the
  full review transcript to `<repo>/PLAN_MULTIQC-REVIEW-LOG.md`. Both go in the repo root,
  as the user asked. Add both to the sdist `exclude` list in `pyproject.toml`, next to
  `PLAN.md`.

## Appendix: review log

_Rounds are appended below. The full raw responses live in the runner artifact directories
and will be copied into `PLAN_MULTIQC-REVIEW-LOG.md`._

### Round 1: Codex gpt-5.6-luna xhigh, verdict REVISE

- **Runner status:** failed ("Plan changed during the run"). The host edited the plan
  mid-run to add the copy step. The response is otherwise complete and was arbitrated.
- **Session:** `01a0c442-ef2b-7e50-b753-8fbad36f50f5`. Plan sha `4676254e…`.
- **User-accepted risk:** the Codex config had write-capable MCP servers enabled (`serena`,
  `node_repl`), and the user chose to run as-is.

**Dispositions**

| ID | Severity | Disposition |
|---|---|---|
| MQC-001 | high | **Partly accepted.** The stale-stats unlink stays in `prepare_output_directories`, which is the first mutation of the output tree. An earlier-rejected run touches nothing, so the old stats still describe coherent old output, and its non-zero exit already says this invocation produced nothing. The invariant is now stated explicitly. |
| MQC-002 | high | **Accepted.** The temp file is dropped (in-place write; exit≠0 output is never consumable anyway), so there is no extra path to alias. Generate gets a minimal alias guard for the stats path, and update adds it to `validate_output_paths`. |
| MQC-003 | high | **Partly accepted.** The module selects sections by the JSON `command` field and keeps a data dict per command. The filename ambiguity with `--chunk_id 1_updated` already exists for every aggregate, so it isn't reserved here. |
| MQC-004 | medium | **Accepted.** `version` comes from `__version__`, and the example is marked illustrative. |
| MQC-005 | medium | **Accepted.** Zero-input semantics are defined, and there is an empty-chunk test plus module handling. |
| MQC-006 | medium | **Accepted**, and verified: v2 is 14 in, 12 successful, 2 discarded, 8 converged. |
| MQC-007 | low | **Accepted.** The tests assert the whole payload against `Counter` projections of the CSVs, plus parameter keys, version and chunk_id. |

Round 1 used up one of the five plan-review rounds. Round 2 is a fresh session, because
`--resume` needs a successful prior result.

### Round 2: Codex gpt-5.6-luna xhigh, verdict REVISE

- **Runner status:** completed.
- **Session:** `01a0c44a-dcc1-7df2-ade7-b3178e80bec8`. Plan sha `71a90c4d…`.

| ID | Severity | Disposition |
|---|---|---|
| MQC-008 | high | **Accepted.** This reverses the round-1 in-place write. The file is now written through an exclusive random temp in the root, then fsync and `os.replace`. A failure unlinks the temp and exits 1. |
| MQC-009 | medium | **Accepted.** The stale unlink moves after update's foreign-family refusal, so it is the first mutation rather than the first statement. There is a test that stats survive the refusal. |
| MQC-010 | high | **Rejected for this plan (scope).** The truncation of generate's existing aggregates when an input is aliased predates this work and is unrelated to MultiQC. The new stats path is guarded, and a full check is recorded as a follow-up issue. |
| MQC-011 | medium | **Accepted.** An empty retention set gives null/omitted mean retention, and there is an all-no-hit test. |
| MQC-012 | low | **Accepted.** The update test projects every count and histogram, and all-discard, all-success and empty cases are covered. |
| MQC-013 | low | **Accepted.** Retention keys are the exact `delta.csv` strings, with no rounding. |

### Round 3: Codex gpt-5.6-luna xhigh, verdict REVISE (resumed from round 2)

- **Runner status:** completed.
- **Session:** `01a0c44a-dcc1-7df2-ade7-b3178e80bec8`. Plan sha `9bc81dd0…`.

| ID | Severity | Disposition |
|---|---|---|
| MQC-014 | medium | **Accepted**, and verified in `update_families.main`: `qlen` and `rounds_run` are set before the no-hit discard. The expected maps are now `model_length_change {0: N}` and `rounds_run {1: N}`, with retention and the metadata maps empty. |
| MQC-015 | low | **Accepted.** The failure test uses an explicit prebuilt index, so the patched `os.replace` can only hit the stats commit. |

### Round 4: Codex gpt-5.6-luna xhigh, verdict REVISE (resumed)

- **Runner status:** completed.
- **Session:** `01a0c44a-dcc1-7df2-ade7-b3178e80bec8`. Plan sha `342f31bb…`.

| ID | Severity | Disposition |
|---|---|---|
| MQC-016 | medium | **Accepted.** `converged` is defined as successful-and-converged for both commands, and there is a direct-builder regression test for a converged-then-discarded delta row. |
| MQC-017 | low | **Accepted.** The update exit-3 and exit-1 tests assert the stats file is present with `exit_status 3`, and absent after exit 1. |
