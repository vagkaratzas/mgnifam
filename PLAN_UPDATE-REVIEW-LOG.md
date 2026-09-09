# PLAN_UPDATE review log

Append-only transcript of the claudex-loop plan review for `PLAN_UPDATE.md`.

## Run configuration

| field | value |
|---|---|
| Host | Claude Code (this session), model Opus 5 |
| Role: requirements + plan + coordination | host |
| Role: plan reviewer | Codex CLI |
| Role: builder | none — design only, by the user's instruction |
| Role: final inspector | N/A — no build in scope |
| `PLAN_FILE` | `PLAN_UPDATE.md` |
| `LOG_FILE` | `PLAN_UPDATE-REVIEW-LOG.md` |
| `MAX_ROUNDS` | 5 |
| `mode` | `review` (a plan already existed) |
| `inspect` | N/A (no build) |
| Requested reviewer model | `gpt-6-astra` |
| Requested reviewer effort | `low` — matching the CLI default `model_reasoning_effort = "low"` in `~/.codex/config.toml`, at the user's instruction. A round was launched at `high` first and killed before completion; no output from it is used. |
| Observed reviewer CLI | `codex-cli 0.153.0` (skill's live-tested baseline is 0.153.4) |
| Runner artifacts | outside the checkout, under the session scratchpad |

## Authorization

The user asked to plan and to refine the plan with Codex. **Implementation is not
authorized by this request** and is explicitly deferred. No edits to `src/`, no commits.

## Pre-flight: MCP audit

`references/runtime.md` requires auditing write-capable MCP integrations before a review,
because the Codex read-only shell sandbox does not restrict MCP side effects.

`~/.codex/config.toml` configures:

- `serena` — exposes `replace_content`, `replace_in_files`, `replace_symbol_body`,
  `insert_after_symbol`, `insert_before_symbol`, `rename_symbol`, `safe_delete_symbol`.
  Only `write_memory` carries `approval_mode = "approve"`.
- `node_repl` — arbitrary Node execution.

These are the user's own configured servers, not unknown tooling. Rather than mutating the
user's global Codex configuration, the mitigation is verification: a SHA256 manifest of
`src/`, `tests/`, `reference/`, `pyproject.toml`, `AGENTS.md`, `README.md` and
`CHANGELOG.md` is taken before each round and re-checked after it. Any drift is reported.

## Rounds

### Round 1 — Codex, requested `gpt-6-astra`, effort `low`

- Runner result: `/tmp/claude-10825/-home-vangelis-Desktop-Projects-mgnifam/7f6048ef-b784-463c-a59f-0df2d8d96cbd/scratchpad/claudex/claudex-s_9lvldz/result.json`
- Observed CLI: `codex-cli 0.153.4` (the runner resolved a newer binary than `codex --version` on PATH reported as 0.153.0)
- Observed model: the provider returned no `observed_models` key — **model identity is unresolved**; `gpt-6-astra` is what was requested, not what is confirmed to have run.
- Plan SHA256 reviewed: `4d9dd89117e831b54a7d999ed0326993cfe5fc875e67a35bbf4b134e0fb2dddf`
- Verdict: **REVISE** — 2 high, 3 medium.
- Integrity check after the round: no drift in `src/`, `tests/`, `reference/`, `pyproject.toml`, `AGENTS.md`, `README.md`, `CHANGELOG.md`.

Full response:

```json
{
  "verdict": "REVISE",
  "summary": "The plan has concrete defects in HMM loading, path safety, deterministic ordering, identity propagation, and rerun cleanup.",
  "findings": [
    {
      "id": "UF-01",
      "severity": "high",
      "path": "PLAN_UPDATE.md:109",
      "evidence": "load_hmms specifies hmm.name.decode(), but the locked pyhmmer version is 0.12.1 and its installed plan7.pyi declares HMM.name as str (lines 366\u2013368). A normal named HMM therefore raises AttributeError before updating any family.",
      "fix": "Use hmm.name directly under the pinned API. Handle absent or empty names separately and verify whether HMMFile can actually parse a nameless model before promising filename fallback."
    },
    {
      "id": "UF-02",
      "severity": "high",
      "path": "PLAN_UPDATE.md:107-118",
      "evidence": "Names are preserved verbatim without validation, while emit_family interpolates them into filesystem paths and unquoted CSV fields. A NAME such as ../../victim escapes the artifact directory; writes and rollback unlink can affect unrelated files. A name containing a comma corrupts metadata and discard CSV columns. Checking only whether output_dir equals the input location prevents neither case.",
      "fix": "Validate every resolved family name before output creation. Reject path separators, traversal components, absolute paths and control characters; either restrict names to a documented safe identifier alphabet or quote CSV fields correctly. Verify every artifact destination remains inside its intended directory."
    },
    {
      "id": "UF-03",
      "severity": "medium",
      "path": "PLAN_UPDATE.md:107-108",
      "evidence": "Directory models are ordered by filenames, whereas library models retain library order. Files a.hmm containing NAME B and z.hmm containing NAME A produce B,A, while a library containing the same models in A,B order produces different aggregate bytes. Renaming files also changes output order despite identity being defined by NAME. This violates acceptance criterion 2.",
      "fix": "After loading and validating names, sort all models by the same canonical family-name key for both input forms. Test reversed library order and filenames whose order differs from NAME order."
    },
    {
      "id": "UF-04",
      "severity": "medium",
      "path": "src/mgnifam/generate_families.py:957",
      "evidence": "The proposed family_id override does not address the representative FASTA writer, which formats its annotation as {chunk}_{provisional_id}. With chunk 9 and preserved family_id 1_7, it emits 9_1_7 rather than 1_7, contradicting identity preservation across outputs.",
      "fix": "Write the resolved family_name into representative FASTA annotations. Its default value already preserves generate_families output. Include representative headers in the preserved-identity test."
    },
    {
      "id": "UF-05",
      "severity": "medium",
      "path": "PLAN_UPDATE.md:208-210",
      "evidence": "Clearing artifacts only by names in the current input cannot clean a shrinking rerun: a chunk first updates A and B, then reruns with only A. B's artifacts remain while the chunk aggregates are overwritten without B. Clearing all names in the shared root instead would delete other chunks' output. The plan provides no persisted ownership rule to resolve this.",
      "fix": "Read the prior chunk's successful-name list before replacing aggregates and clear artifacts owned by the union of previous and current names. Define that chunks sharing a root must own disjoint family names, and test a shrinking rerun plus preservation of another chunk's files."
    }
  ],
  "coverage": [
    "Verified PLAN_UPDATE.md SHA256 matches the supplied digest.",
    "Read the supplied plan, README, repository guidance, and relevant changelog sections.",
    "Traced Family.advance, finish, discard, membership filtering, search, emit_family, Writers, rollback, shared-output commits, and main-loop re-emission.",
    "Inspected CLI dispatch, output cleanup, index resolution, fixture definitions, relevant test references, and legacy Stockholm serialization.",
    "Checked the dependency declaration, locked pyhmmer version, and installed HMM API type declarations."
  ],
  "limitations": [
    "Static review only; no files edited, tests run, or work delegated.",
    "The proposed update_families implementation does not yet exist in the inspected source tree.",
    "External downstream consumers and the 1.6B-sequence release were not inspected.",
    "Runtime reproducibility, nameless-HMM parsing, and scientific output equivalence remain unverified."
  ]
}
```

#### Host dispositions

All five accepted. Each was independently verified against the repository before acceptance;
none was taken on the reviewer's assertion alone.

| id | disposition | verification |
|---|---|---|
| UF-01 `hmm.name.decode()` | **Accepted** | Confirmed twice: `pyhmmer/plan7.pyi:364-368` declares `name -> str`, and executing against `output/hmm/1_1.hmm.gz` under `uv run` returns `'1_1'`, a `str`. `.decode()` would have raised `AttributeError` on the first model of every run. Went further than the fix suggested: stripping `NAME` from a real model and reopening it raises `ValueError: Invalid format in file: No NAME found for HMM`, so a nameless model cannot reach the loader and the filename-stem fallback was unreachable code. It is deleted rather than corrected. |
| UF-02 unvalidated names | **Accepted** | Correct and the most serious finding. `generate_families` never needed name validation because both halves of a family name were trusted — `--chunk_num` is checked against `CHUNK_PATTERN` (`:1013`) and the id is an integer. Taking the name from an arbitrary third-party HMM removes that guarantee while `emit_family` still interpolates it into paths (`:906-923`) and into unquoted CSV fields (`:950,869`). Fixed by requiring `CHUNK_PATTERN`'s own alphabet plus rejection of `.`/`..`, which closes traversal and CSV corruption in one rule, plus a path-containment assertion as defence in depth. |
| UF-03 load ordering | **Accepted** | Sound. Filename order and library order are genuinely different orders over the same models, so acceptance criterion 2 as written was unsatisfiable. Models are now sorted by family name for both input forms, and the adversarial case (filenames sorting opposite to their `NAME`s) is a test. |
| UF-04 reps annotation | **Accepted** | Verified at `generate_families.py:956`: `f">{sequence_name}\t{chunk}_{provisional_id}\n"` rebuilds the name from `chunk` and the id. With `--chunk_num 9` and a preserved `1_7` this emits `9_1_7`. Identity would have held in filenames and broken in the one file that maps a representative back to its family. Acceptance criterion 3 was too narrow — it said "artifacts are named by NAME" and is now explicit about identity-bearing fields *inside* the artifacts. |
| UF-05 shrinking rerun | **Accepted**, fix narrowed | The defect is real and is the same mixed-two-runs state `prepare_output_directories` exists to prevent (`:1051-1056`). Adopted the reviewer's fix: clear the union of the previous run's `<chunk>_updated_successful.txt` and the current input's names. Narrowed the justification — only *successful* families write per-family artifacts, so that file is already the exact and complete ownership record; no new bookkeeping file is needed. The disjoint-names constraint is documented as a constraint, not claimed as an enforced guarantee, because nothing can enforce it when names come from the input models. |

One defect found by the host, not the reviewer, and fixed in the same revision: the round-1
plan said `emit_family` "must return `success_count + 1` unconditionally". `emit_family` has
two returns — the discard path's `return success_count` (`:872`) is correct as it stands.
Only the success-path return (`:962`) changes. The revision says so explicitly.

### Round 2 — Codex, requested `gpt-6-astra`, effort `low` (resumed session)

- Runner result: `/tmp/claude-10825/-home-vangelis-Desktop-Projects-mgnifam/7f6048ef-b784-463c-a59f-0df2d8d96cbd/scratchpad/claudex/claudex-_pmu9s99/result.json`
- Plan SHA256 reviewed: `2970002b39d5297a5246f14ab51c41f40c133228124ed13996e4e17e497c9c01`
- Verdict: **REVISE** — UF-01 through UF-05 confirmed resolved; 2 new medium findings.
- Observed model: still no `observed_models` returned — identity remains unresolved.
- Integrity check after the round: no drift.

Full response:

```json
{
  "verdict": "REVISE",
  "summary": "UF-01 through UF-04 are resolved. UF-05 handles completed reruns, but new evidence exposes an interrupted-run gap. The delta report also needs explicit handling of unavailable metrics.",
  "findings": [
    {
      "id": "UF-06",
      "severity": "medium",
      "path": "PLAN_UPDATE.md:251-254",
      "evidence": "The previous successful.txt is not a complete ownership record after a failed run. emit_family writes all per-family artifacts before appending to successful_clusters (generate_families.py:903-937). If the first shared append fails, ChunkCorrupted leaves those artifacts on disk without a successful-name entry. Rerunning with that family omitted leaves orphaned artifacts that the proposed cleanup cannot discover. Process termination before the buffered successful list is flushed has the same effect.",
      "fix": "Persist the chunk's complete intended ownership set before writing artifacts, retaining previous ownership until cleanup succeeds. Alternatively, explicitly reject reuse after an incomplete run and require a fresh output root. Add a failed-commit followed by shrinking-rerun test."
    },
    {
      "id": "UF-07",
      "severity": "medium",
      "path": "PLAN_UPDATE.md:194-203",
      "evidence": "The delta schema requires check_seed_membership's value for every discard, but a no-hit family has no round-1 members and that helper divides by zero for an empty original set (generate_families.py:529-531). Other discards lose the evidence needed to reconstruct retention: finish stores membership only in a local variable, and discard clears records, HMM and alignments (lines 648-656, 755-770). For example, a family rejected for representative length has already computed retention, then destroys the recruited records without preserving that value. The listed changes provide no metric capture or unavailable-value rule.",
      "fix": "Define nullable delta fields for stages that were never reached, including retention when round-1 membership is empty. Capture scalar metrics before discard clears scientific objects, retaining the existing memory-release behavior. Specify that delta is committed only after the final emission outcome is known, and that a delta-write failure after shared output commits is fatal rather than contained as another discard. Test no-hit, length-discard, and artifact-failure rows."
    }
  ],
  "coverage": [
    "Verified the revised PLAN_UPDATE.md SHA256 matches the supplied digest.",
    "Checked all five prior findings against the revision; confirmed the success-only return-value correction.",
    "Re-inspected Family.discard, finish, check_seed_membership, emit_family, shared-output commit ordering, rollback and main-loop re-emission.",
    "Inspected cleanup behavior and relevant artifact-failure test source; carried forward repository, CLI, fixture and dependency evidence from the prior review."
  ],
  "limitations": [
    "Static review only; no edits, tests, execution of scientific code, or delegation.",
    "update_families remains a proposed implementation.",
    "External downstream consumers, runtime reproducibility and the new protein database were not verified."
  ]
}
```

#### Host dispositions

Both accepted. UF-06 corrects a claim the host made in round 1, which is the more useful
outcome than either party being right the first time.

| id | disposition | verification |
|---|---|---|
| UF-06 ownership record | **Accepted**; the round-1 fix was wrong and is replaced. The host asserted `<chunk>_updated_successful.txt` was "the exact and complete ownership record". It is not. `emit_family` writes every per-family artifact (`:904-934`) *before* its first shared append (`:936`), so a `ChunkCorrupted` on that append — or a `SIGKILL` before the buffered handle flushes — leaves artifacts on disk with no name in that file. Those are exactly the runs that exit 1 and mandate a rerun, so the cleanup would have been blind precisely when needed. Replaced with an intended-ownership manifest, `<chunk>_updated_manifest.txt`, written and `fsync`ed at startup from the validated input name list, before any artifact exists. Chose that over the reviewer's alternative of refusing to reuse an output root after an incomplete run: the manifest costs one line per family and keeps rerun-in-place working, which is the normal recovery path for an exit-1 chunk. Ordering detail added beyond the fix as stated — the manifest is written as the union *before* clearing starts, so an interruption mid-clear still leaves a record covering both what was removed and what was not. |
| UF-07 delta metrics | **Accepted** in full; three distinct defects, all confirmed. (i) `check_seed_membership` divides by `len(original_first_parts)` (`:531`), so a no-hit family with no round-1 members yields `ZeroDivisionError`, not a zero score — `retention` must be nullable, and the schema wrongly required it for every row. (ii) `discard()` clears `records`, `hmm`, `seed_msa`, `full_msa` and `total_checked_sequences` (`:652-656`), and `finish()` holds `membership` in a local before possibly discarding on representative length two statements later (`:753-765`) — a length-discarded family has computed its retention and destroyed the evidence, so metrics must be captured at the point of computation rather than read back off the `Family`. The existing memory-release behaviour is left alone. (iii) The delta row is now committed inside `emit_family`'s shared-append `try` (`:936-961`), which makes a delta-write failure fatal for free: a delta missing for a family already in `successful.txt` is the run's summary contradicting its own output, which is the definition of incoherent output under the exit-1 contract. |

### Round 3 — Codex, requested `gpt-6-astra`, effort `low` (resumed session)

- Runner result: `/tmp/claude-10825/-home-vangelis-Desktop-Projects-mgnifam/7f6048ef-b784-463c-a59f-0df2d8d96cbd/scratchpad/claudex/claudex-jmgeqfho/result.json`
- Plan SHA256 reviewed: `31d6343154759d26f07e64c7b7305a8f09d25f32fb671fc8cddec70c53a34a62`
- Verdict: **REVISE** — earlier fixes accepted; 2 new medium findings.
- Observed model: no `observed_models` returned; identity remains unresolved.
- Integrity check after the round: no drift.

Full response:

```json
{
  "verdict": "REVISE",
  "summary": "Earlier fixes remain accepted, but the revised recovery and delta-write mechanisms leave two concrete gaps.",
  "findings": [
    {
      "id": "UF-08",
      "severity": "medium",
      "path": "PLAN_UPDATE.md:280-289",
      "evidence": "Flushing and fsyncing the replacement manifest does not protect the previous manifest during replacement. If a rerun opens the existing manifest for writing and is interrupted after truncation or a partial write, previous ownership names disappear before cleanup begins. A subsequent shrinking rerun can no longer discover those families' artifacts. The plan specifies durability after writing, but no atomic replacement.",
      "fix": "Write the sorted union to a temporary file beside the manifest, flush and fsync it, then atomically replace the manifest before cleanup. Fsync the parent directory if power-loss durability is required. Test interruption during manifest replacement, ensuring the prior ownership record remains usable."
    },
    {
      "id": "UF-09",
      "severity": "medium",
      "path": "PLAN_UPDATE.md:222-226",
      "evidence": "The specified delta-write location is unreachable for discarded families: emit_family writes discarded.csv and returns at generate_families.py:868-872, before the shared-append try at :936. Placing delta emission only in that try therefore omits every discard, including the newly specified no-hit and length-discard rows. Simply adding a delta write before the discard return also needs protection: if it fails after discarded.csv was appended, family_guard can contain the exception and re-emit a duplicate discard.",
      "fix": "Explicitly cover both success and discard emission paths. Commit the discard row and its delta row under a ChunkCorrupted boundary, preserving the unchanged success count. Test that each discard produces exactly one delta row and that a delta failure after a discard append exits 1 without re-emission."
    }
  ],
  "coverage": [
    "Verified PLAN_UPDATE.md SHA256 matches the supplied digest.",
    "Checked the revised dispositions against the prior review findings.",
    "Re-inspected emit_family's success and discard branches, family_guard, shared-write exception handling, and the proposed manifest lifecycle.",
    "Carried forward earlier inspection of recruitment, metric lifetime, cleanup, callers, fixtures, CLI and dependency declarations."
  ],
  "limitations": [
    "Static plan review only; no files edited, tests run, or work delegated.",
    "The proposed implementation and its runtime behavior remain unverified.",
    "External downstream consumers and the new protein database were not inspected."
  ]
}
```

#### Host dispositions

| id | disposition | verification |
|---|---|---|
| UF-08 manifest atomicity | **Accepted.** Durability after writing is not the same as durability *of the thing being replaced*; truncating the previous manifest and dying leaves a state strictly worse than having no manifest, because a later shrinking rerun can then never discover the stranded artifacts. The repository already contains the correct idiom for exactly this reason — `build_ssi_index` writes into a private temporary directory on the destination filesystem and `os.replace`s it into place (`:160-196`) — so the manifest follows it rather than inventing a scheme: temp file beside the destination, `flush`, `fsync`, `os.replace`, and `fsync` of the parent directory where power-loss durability is wanted. |
| UF-09 delta unreachable for discards | **Accepted**, and it surfaced a pre-existing defect in `generate_families` that this plan now fixes. Verified at `:868-872`: the discard branch appends to `discarded.csv` and returns *before* the shared-append `try` at `:936`, so a delta write placed only in that `try` would omit every discard — including the no-hit and length-discard rows the report exists to carry. Verified the second half too: that append is not under any `ChunkCorrupted` boundary, so a failure there is caught by `family_guard` in `main` (`:1256-1266`), which re-emits, which calls `emit_family` again and appends the discard row a second time. This is the same class of defect the 2.1.0.dev0 release already claims to have fixed — *"Shared-output commit failures now bypass per-family containment and exit 1"* — in the one branch that hardening missed. The discard branch therefore gains the same boundary, covering its `discarded.csv` append and its delta row together. |

**Scope note.** The UF-09 fix is the first change in this plan that is *not*
default-preserving: a failed `discarded.csv` append now exits 1 instead of being contained and
retried. It is declared as a deliberate widening rather than absorbed silently — the decisions
table, the non-goals section, a `CHANGELOG.md` *Fixed* entry and a dedicated test all name it.
Output bytes on a clean run are unaffected, so acceptance criterion 5 still holds.

### Round 4 — Codex, requested `gpt-6-astra`, effort `low` (resumed session)

- Runner result: `/tmp/claude-10825/-home-vangelis-Desktop-Projects-mgnifam/7f6048ef-b784-463c-a59f-0df2d8d96cbd/scratchpad/claudex/claudex-kss6y6go/result.json`
- Plan SHA256 reviewed: `46987fff3b01f438b8b5740ae002907fefc9fdff00d05987d39e6401faea95de`
- Verdict: **APPROVED**, zero findings.
- Approval check: `runner.py check` returned *"Approval matches the current plan."* (exit 0), binding the approval to that exact hash.
- Integrity check after the round: no drift.

Full response:

```json
{
  "verdict": "APPROVED",
  "summary": "The revision resolves UF-08 and UF-09 while retaining the earlier fixes. No material unresolved plan defects identified.",
  "findings": [],
  "coverage": [
    "Verified PLAN_UPDATE.md SHA256 matches the supplied digest.",
    "Checked atomic manifest replacement, ownership preservation and interruption recovery.",
    "Checked delta emission on both outcome paths against emit_family, family_guard and main-loop re-emission.",
    "Confirmed the discard-write behavior change is explicitly scoped and includes planned regression coverage.",
    "Carried forward prior review of loading, identity propagation, recruitment, metric lifetime, rollback, cleanup, CLI and dependency evidence."
  ],
  "limitations": [
    "Approval covers the plan, not an implementation.",
    "Read-only review; no files edited, tests run or work delegated.",
    "Runtime reproducibility, scientific equivalence, external consumers and the new protein database remain unverified."
  ]
}
```

## Outcome

**APPROVED** after 4 of 5 rounds, bound to `PLAN_UPDATE.md` at SHA256
`46987fff3b01f438b8b5740ae002907fefc9fdff00d05987d39e6401faea95de`. Nine findings were
raised across three review rounds and all nine were accepted; two of them corrected claims
the host had made, and one exposed a defect in already-shipped code.

Zero findings in the final round is a valid result and is not proof of exhaustive
correctness. The reviewer's own stated limitations carry forward.

### What the approval does not cover

- **The implementation.** Approval is of the plan only. No code was written and none was
  reviewed.
- **The three open questions.** OQ1 (additive parameters vs. a duplicate emit), OQ2 (the
  membership yardstick read as *measured* rather than *pinned*) and OQ3 (no `seed_msa/` or
  `rf/` under `--skip_refine`) were put to the reviewer explicitly. It raised no finding
  against any of them across four rounds, which is silent assent, not endorsement. **OQ2 in
  particular is a scientific choice the user must settle** — it decides whether a refined
  family may shed a sequence its own round-1 model recruited. A reviewer cannot answer that
  from the source.
- **Runtime behaviour.** Reproducibility across `--cpus` / `--batch_size` /
  `--prefetch_targets` / `PYTHONHASHSEED`, scientific equivalence, and the byte-identity
  regression baseline are all specified in the plan and all unverified until the code exists.
- **The 1.6B release.** Whether it keeps the `<protein>_<start>_<end>` fragment-naming
  convention was not verifiable here. `split_slice_name`'s span test degrades safely
  (`generate_families.py:557`), so this is a silent-precision risk, not a crash.

### Model identity

`gpt-6-astra` at effort `low` was **requested** on every round. The provider returned an
empty `observed_models` on all four, so the model that actually ran is **unresolved** and is
reported as requested-only. Observed CLI was `codex-cli 0.153.4` — note that `codex --version`
on `PATH` reports `0.153.0`; the runner resolved a different, newer binary.

### Authorization state

Implementation remains **unauthorized**. The user asked to plan and to refine the plan with
Codex; neither request authorizes building. No files under `src/` or `tests/` were modified,
and nothing was committed. Per the skill, the builder would be this host and a fresh Codex
session would inspect the result — neither has been engaged.
