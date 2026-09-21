# PLAN_MULTIQC review log

Append-only transcript of the claudex-loop plan review of `PLAN_MULTIQC.md`.

- Roles: host/planner/builder = Claude (Opus 5); plan reviewer and final inspector = Codex.
- Requested reviewer model: `gpt-5.6-luna`, effort `xhigh` (user request). Codex CLI `codex-cli 0.155.1`. Observed model keys: not returned by the CLI.
- Scope: MultiQC support for `generate_families` and `update_families` (stats JSON + drafted native MultiQC module).
- Authorization: planning only; build requires explicit user approval. Round limit: 5.
- Risk accepted by user: Codex config had write-capable MCP servers (`serena`, `node_repl`) enabled during review; the read-only sandbox does not restrict MCP side effects.
- Host dispositions for rounds 1-4 are in the appendix of `PLAN_MULTIQC.md`.
- Known wording nit (not a reviewer finding): the "Code placement" bullet still says `json.dumps` + `write_text`; the Lifecycle section's atomic temp-file + `os.replace` write is authoritative.

## Round 1: verdict REVISE

- Runner status: `failed` — error: Plan changed during the run; result cannot approve the current plan.
- Session: `01a0c442-ef2b-7e50-b753-8fbad36f50f5`
- Plan SHA256 reviewed: `4676254e3ad8647a3f9d034f1d16bb1bdbf8c3ee416d14aef38856851050f784`
- Resumed: no
- Elapsed: 445.37 s

### Summary

The plan has material lifecycle, aliasing, naming, and edge-case defects that can leave stale or misleading MultiQC data or destroy inputs.

### Findings

#### MQC-001 (high) — `src/mgnifam/generate_families.py:1264-1272; src/mgnifam/update_families.py:156-221,387-395; plan:115-118`

**Evidence:** Stale stats are removed only inside prepare_output_directories(), but validation runs first. A rerun with a missing/malformed input, unsafe model name, or rejected output collision exits before unlinking the old stats file; MultiQC then sees a prior run as consumable output.

**Fix:** After safely validating the stats path cannot alias an input, remove the prior stats file before any remaining preflight that can fail. Add tests for validation and rejected-retry failures, not only contained runtime crashes.

#### MQC-002 (high) — `src/mgnifam/update_families.py:225-264; src/mgnifam/generate_families.py:1264-1272; plan:119-124`

**Evidence:** The plan protects only the final update stats path, but writes through <name>.tmp. A single HMM library named 1_updated_stats.json.tmp is accepted as input, then truncated by the stats writer. Generate has no output-alias validation at all; a valid FASTA or cluster file named 1_stats.json can likewise be replaced.

**Fix:** Validate final and temporary stats paths against every input using resolved paths and inode identity for both commands, or place temporary files in a private validated directory.

#### MQC-003 (high) — `src/mgnifam/generate_families.py:63; src/mgnifam/update_families.py:314-316; plan:32-36,40-41,177-178`

**Evidence:** Safe chunk IDs allow underscores. Generate with --chunk_id 1_updated and update with --chunk_id 1 both write 1_updated_stats.json; the later atomic replace destroys the earlier report, and staged files have the same MultiQC sample name.

**Fix:** Reserve the _updated suffix for update output, or introduce an unambiguous command namespace in filenames. The module should also use the JSON command field, not filename shape, to select update-only sections.

#### MQC-004 (medium) — `pyproject.toml:3; src/mgnifam/__init__.py:3-10; plan:63,218`

**Evidence:** The schema example specifies version 3.1.0 while the repository is version 3.0.0 and the plan explicitly makes a version bump a non-goal. A literal implementation would report a false software version in current builds and could drift on later releases.

**Fix:** Require the payload to use the package __version__ value; mark 3.1.0 as illustrative and regenerate external test data after the release version is actually set.

#### MQC-005 (medium) — `src/mgnifam/generate_families.py:1107-1115,1154-1190; plan:183-188`

**Evidence:** An empty cluster TSV is accepted: load_clusters() returns an empty dict and main completes with zero families. The planned General Stats percentage requires successful/input, which is 0/0 unless the module explicitly handles this valid output.

**Fix:** Define the zero-input representation (for example null or 0%) and require the MultiQC module to render it without division-by-zero or treating the sample as absent. Add an empty-chunk test.

#### MQC-006 (medium) — `plan:232-244; tests/test_generate_families.py:1552-1579`

**Evidence:** The manual verification expects families.successful == 14 for mgnifams_v2, but the repository’s end-to-end contract has 14 input clusters, 12 successful, 2 discarded, and 8 converged. The plan’s own schema example says 12/2, so the verification instructions are internally inconsistent with the fixture.

**Fix:** Change the expected values to input 14, successful 12, discarded 2, and converged 8, and use those values when validating the generated test data.

#### MQC-007 (low) — `plan:139-153`

**Evidence:** The proposed tests check outcome counts, discard reasons, and skip_refine, but do not assert histogram contents, all parameter fields/exclusions, version, chunk_id, or update model-length/round/retention bins. An implementation can therefore pass the planned tests while violating the schema’s metric requirements.

**Fix:** Assert the complete payload against projections computed from every relevant CSV column, including all histogram maps and deterministic parameter serialization; add malformed/preflight and empty-output cases.

### Coverage

- Verified the supplied plan SHA256.
- Read the generate_families and update_families lifecycle, validation, aggregate writers, emit paths, CLI options, package version, README/AGENTS contracts, and relevant existing tests.
- Traced stats placement against validation order, output alias checks, aggregate naming, and current fixture expectations.

### Limitations

- No tests or builds were run, per instruction.
- No sibling MultiQC or test-data checkout was present under the inspected project parent, so the external module API, strict-mode behavior, DOI linting, and drafted files could not be independently inspected.

## Round 2: verdict REVISE

- Runner status: `completed`
- Session: `01a0c44a-dcc1-7df2-ade7-b3178e80bec8`
- Plan SHA256 reviewed: `71a90c4dc02cdc17f8501980afa772d0db24960cf45d4f3338919cdf4c8a89b4`
- Resumed: no
- Elapsed: 655.41 s

### Summary

The plan has material lifecycle, atomicity, input-protection, and MultiQC edge-case defects.

### Findings

#### MQC-008 (high) — `PLAN_MULTIQC.md:123-135,187-188`

**Evidence:** The plan explicitly uses in-place Path.write_text(). ENOSPC, EIO, or SIGKILL after truncation can leave a partial <chunk>_stats.json, even though exit 1 or incomplete execution is required to leave no stats file. The docs also call the write atomic.

**Fix:** Write to a temporary file, flush/fsync and close it, then atomically replace the final path; remove the temporary file on failure.

#### MQC-009 (medium) — `src/mgnifam/update_families.py:267-291; PLAN_MULTIQC.md:123-129`

**Evidence:** Adding stale-stats removal as the first statement of prepare_output_directories() runs before the existing foreign-family preflight. A retry with models A only over output containing A and B will delete the old stats file, then raise the foreign-family ValueError, leaving coherent old artifacts without their completion marker.

**Fix:** Perform all read-only refusal checks, including foreign-family detection, before removing stale stats. Make stale removal the first mutation, not the first statement.

#### MQC-010 (high) — `src/mgnifam/generate_families.py:1302-1318; PLAN_MULTIQC.md:136-143`

**Evidence:** The planned generate guard protects only the new stats path. Existing aggregate writers open paths such as <chunk>_families.tsv with mode w. A valid FASTA supplied at that path is accepted, then truncated during writer setup, destroying an input. The plan acknowledges this destructive alias as out of scope but leaves it unresolved.

**Fix:** Validate every generate output destination against clusters, FASTA, explicit SSI, and hard/symlink aliases before any output mutation, or explicitly remove the unsafe invocation from the supported contract.

#### MQC-011 (medium) — `PLAN_MULTIQC.md:216-221; src/mgnifam/update_families.py:103-128,492-498`

**Evidence:** An update against a database where every model has no reported hits produces delta rows with empty retention because Family.membership remains None. The plan requires mean retention in General Stats but defines no behavior for an empty retention set; a direct mean calculation divides by zero or reports a misleading zero.

**Fix:** Define empty retention as null/omitted for that sample, filter missing values, and add an all-no-hit fixture test.

#### MQC-012 (low) — `PLAN_MULTIQC.md:174-178`

**Evidence:** Acceptance criterion 4 covers every count, but the planned update test only explicitly checks outcomes and model_length_change, rounds_run, and retention. It does not verify full_msa_size, model_length, representative_length, discard-reasons, or converged counts against metadata/delta CSVs.

**Fix:** Project and compare every top-level count and histogram from both update aggregate CSVs, including all-success, all-discard, and empty-histogram cases.

#### MQC-013 (low) — `PLAN_MULTIQC.md:78-93,216-221`

**Evidence:** Retention values are rounded to two decimal keys before serialization, while the plan calls histogram maps lossless and uses them to compute mean retention. A true value such as 1/3 becomes 0.33 and cannot recover the exact mean.

**Fix:** Either store exact retention values and round only in the MultiQC display, or add an exact retention sum/count (and document the two-decimal histogram as intentionally lossy).

### Coverage

- Verified the supplied plan SHA256.
- Inspected generate_families and update_families lifecycle, writers, aliases, delta generation, main-loop exit paths, tests, README, AGENTS.md, CHANGELOG.md, pyproject.toml, and version/CLI wiring.
- Traced shared output state through prepare_output_directories(), validate_output_paths(), emit_family(), delta_row(), and the aggregate writers.
- Consulted the current MultiQC module-development documentation: https://github.com/MultiQC/MultiQC/blob/main/docs/markdown/development/modules.md

### Limitations

- No implementation was available for inspection; no tests, builds, or linters were run, as requested.
- The planned ../MultiQC and ../test-data sibling checkouts do not exist, so the drafted external module and test-data files could not be reviewed directly.

## Round 3: verdict REVISE

- Runner status: `completed`
- Session: `01a0c44a-dcc1-7df2-ade7-b3178e80bec8`
- Plan SHA256 reviewed: `9bc81dd05b1bbd7e3b7bbbdb94426f19e379d72230b366dd4685d825c99010f7`
- Resumed: yes
- Elapsed: 389.92 s

### Summary

Two unresolved issues remain: an incorrect no-hit expectation and a stats-write failure test that can intercept SSI creation instead.

### Findings

#### MQC-014 (medium) — `PLAN:172-200; src/mgnifam/update_families.py:447-498`

**Evidence:** The plan requires update histograms to project delta.csv, but says the all-no-hit case has empty length histograms. The source sets model_length_before, family.qlen/model_length_after, and rounds_run=1 before discarding a no-hit family; delta_row then writes those fields. For unchanged loaded models, model_length_change is therefore 0 and rounds_run is 1, even though retention and metadata-derived histograms are empty.

**Fix:** Clarify and test the expected maps: for N no-hit families, model_length_change should contain 0:N and rounds_run 1:N; retention and metadata-derived histograms should be empty. Alternatively explicitly redefine the schema to exclude discarded families and update the projection requirement.

#### MQC-015 (low) — `PLAN:196-197; src/mgnifam/generate_families.py:143-198,1255-1260; tests/test_generate_families.py:39-78`

**Evidence:** The planned stats-write test monkeypatches os.replace globally. With the existing test helper's default fasta_index=None, generate_families first calls resolve_index(), whose default SSI build also calls os.replace before stats writing. The test can therefore pass by failing SSI installation and never exercising atomic stats replacement or temp cleanup.

**Fix:** Run the failure test with a valid prebuilt explicit --fasta_index, or inject/patch only the stats replacement after index setup so the failure occurs in the stats commit path.

### Coverage

- Verified the supplied plan SHA256 and reviewed the revised lifecycle, schema, tests, generate/update callers, aggregate writers, SSI creation, alias validation, and exception paths.
- Checked the update no-hit path and delta serialization beyond the plan's listed files.

### Limitations

- No implementation changes, tests, builds, or linters were run.
- The external ../MultiQC and ../test-data checkouts were not inspected; their drafted module/API compatibility remains unverified.

## Round 4: verdict REVISE

- Runner status: `completed`
- Session: `01a0c44a-dcc1-7df2-ade7-b3178e80bec8`
- Plan SHA256 reviewed: `342f31bbf6b0e2dbba9014f7345af681af3179d5483436b36d3ce6be78f8d7b9`
- Resumed: yes
- Elapsed: 507.19 s

### Summary

The revised lifecycle and no-hit cases address prior findings, but update convergence semantics remain inconsistent and update exit-status coverage is incomplete.

### Findings

#### MQC-016 (medium) — `PLAN:103-106,186-190; src/mgnifam/generate_families.py:709-717,775-777,822-835; src/mgnifam/update_families.py:115-125`

**Evidence:** The planned update projection takes families.converged from delta.csv. A family can set ever_converged=True in advance(), then be discarded by finish() for membership or representative length; discard() preserves ever_converged, and delta_row writes it even for discarded outcomes. Generate metadata, however, only contains successful families, so the same field excludes such discarded families. The report can therefore count a discarded update family as converged and make the General Stats meaning differ between commands.

**Fix:** Define converged explicitly as successful final families and count only delta rows with outcome=successful (or use converged.txt) for both commands. Add a regression case where a family converges before finish() discards it. If ever-converged semantics are intended, rename/document the field and make generate report the same semantics.

#### MQC-017 (low) — `PLAN:183-207; tests/test_update_families.py:422-470`

**Evidence:** The plan explicitly tests stats on a contained generate crash, but the existing update tests that exercise exit 3 and exit 1 only inspect delta/discarded outputs. They do not assert that update writes exit_status=3 stats or removes a stale stats file on corrupted output, so a misplaced update stats call can pass the planned update tests while violating acceptance criteria 1 and 2.

**Fix:** Extend the update contained-crash test to assert `<chunk>_updated_stats.json` exists with exit_status 3 and crashed=1, and extend the update fatal-output test to seed a stale stats file and assert it is absent after exit 1.

### Coverage

- Verified the supplied plan SHA256.
- Inspected generate_families, update_families, CLI dispatch, aggregate writers, Family state transitions, existing generate/update tests, README, AGENTS.md, and pyproject.toml.
- Checked related callers and shared emit/delta-writing paths beyond the plan's listed files.

### Limitations

- No implementation changes, tests, builds, or linters were run.
- The external MultiQC and test-data sibling checkouts were not inspected, so their drafted module/API compatibility remains unverified.

## Round 5: verdict APPROVED

- Runner status: `completed`
- Session: `01a0c44a-dcc1-7df2-ade7-b3178e80bec8`
- Plan SHA256 reviewed: `215f5d3f242ac44bf1272352fd9d9511d00e72e91ca466c44ba75c7c0a7869b3`
- Resumed: yes
- Elapsed: 84.78 s

### Summary

No material unresolved defects found in the revised plan. The previously identified convergence semantics, no-hit handling, atomic-write, aliasing, and update exit-path issues are explicitly addressed.

### Findings

None.

### Coverage

- Verified the supplied plan SHA256.
- Inspected generate_families and update_families callers, shared Family/emit/delta writers, output preparation and alias validation, CLI dispatch, relevant tests, README, AGENTS.md, and pyproject.toml.
- Traced both successful and discarded/crashed update paths and the stats lifecycle.

### Limitations

- No implementation changes, tests, builds, or linters were run.
- The external MultiQC and test-data sibling checkouts were not inspected; their drafted module/API compatibility remains unverified.

## Approval check

`runner.py check` against the round-5 result: "Approval matches the current plan." (SHA256 `215f5d3f242ac44bf1272352fd9d9511d00e72e91ca466c44ba75c7c0a7869b3`). `PLAN_MULTIQC.md` is a byte-identical copy of the approved plan.

## Post-approval plan change (user-directed, not reviewed by Codex)

During the build the user enabled Zenodo and asked for `CITATION.cff` plus a Zenodo README
badge. `PLAN_MULTIQC.md` gained a "Citation and Zenodo" section, and `CITATION.cff` left the
non-goals. This change is outside the round-5 approval (SHA `215f5d3f...`) and will be covered
by the final Codex inspection of the built diff.

## Final inspections: Codex gpt-5.6-luna xhigh, fresh sessions, one per repository

- All three ran against the pre-build commits: mgnifam `9f31bba`, MultiQC `deca1b79`, test-data `8aa7630`.
- All three returned REVISE with runner status `failed` ("Code changed during inspection"). The host committed the user-requested `v3.1.0` bump while they ran. The findings were arbitrated anyway.
- No re-inspection, at the user's instruction. The fixes below are therefore unreviewed by Codex.

### mgnifam

| ID | Severity | Disposition |
|---|---|---|
| F001 | medium | **Rejected (user decision).** The user removed the PLAN_MULTIQC exclusions from the sdist on purpose, and these files will be deleted after the work completes. If 3.1.0 is released first, they ship in its sdist. |
| F002 | medium | **Accepted.** `mkstemp` forced mode 0600 on the stats file. It now uses an exclusive `"x"` open under a random name, so the umask applies like every other output. Test: the stats file mode equals the metadata CSV mode. |
| F003 | low | **Accepted.** Temp cleanup is best-effort under `contextlib.suppress(OSError)`, only for a file this call created, and `ChunkCorrupted` is always raised. |
| F004 | low | **Accepted.** `--recruit_evalue_cutoff` must be finite in both commands (`inf`/`nan` rejected in validation), and the JSON is dumped with `allow_nan=False`. |

### MultiQC module

| ID | Severity | Disposition |
|---|---|---|
| MQC-001 | high | **Deferred (planned).** `doi=None` stays until Zenodo mints the DOI at release. |
| MQC-002 | medium | **Rejected.** Update-only sections are selected by the JSON `command` field, as planned. The shared outcome and size sections intentionally compare distinctly named generate and update samples in one plot. |
| MQC-003 | medium | **Accepted.** Retention is binned into twentieths for display; the raw data stays exact. Test added. |
| MQC-004 | medium | **Accepted.** Invalid JSON and a missing or unsupported `schema_version` are warned about and skipped. Test added. |
| MQC-005 | low | **Accepted.** Added to `module_order` after `mgikit`, although recent modules (riker, seqkit) were merged without an entry. |
| MQC-006 | medium | **Fixed**, as for test-data MQC-002. |

### test-data

| ID | Severity | Disposition |
|---|---|---|
| MQC-001 | high | **Partly accepted.** The default-flag v2 run was intentional, but it is regenerated with `--discard_min_rep_length 100` (14/12/2/8) so the generate sample shows discards. |
| MQC-002 | medium | **Accepted.** All three files are regenerated from the 3.1.0 build. |
