# Plan Review Log: crashed families must fail the chunk loudly (2.1.0.dev0)

Act 1 (grill) complete — plan locked with the user. MAX_ROUNDS=5.

Plan under review: `PLAN-2.1.0.md`. The repo's `PLAN.md` / `PLAN-REVIEW-LOG.md` are the
historical record of the original port and are not touched (`AGENTS.md:127`).

## Act 1 — decisions the grill settled

1. **Abort point.** Finish the chunk and exit non-zero at the end, rather than aborting at
   the crash. `family_guard` exists so one family's failure does not destroy the batch's
   completed search.
2. **Detection.** Reuse `family_guard`'s existing `discard_reason` prefix via a shared
   constant — no new `Family` field, no signature change, no third `FamilyState`.
3. **Exit code.** `3`, distinct from `1` (uncaught) and `2` (argparse).
4. **Hole-1 length check.** Keep. Both paths exit 3; the check decides whether the exiting
   chunk is self-consistent.
5. **Version.** `pyproject.toml` / `uv.lock` already at `2.1.0.dev0`; CHANGELOG heading
   already present and empty.

## Round 1 — Codex

Reviewer: `gpt-5.6-sol` (config-pinned, reasoning effort medium), codex-cli 0.146.0.
Thread: `019fd1c7-0841-7c62-bee5-8fa850990aab`.

1. **Critical — the claimed output invariant is still false.** Shared writes at
   `emit_family()` lines 908–927 can raise or partially write; e.g.
   `converged_families.write()` may succeed before `successful_clusters.write()` fails,
   leaving artifacts plus a converged ID before `family_guard` reclassifies the family as
   discarded. "Uninterrupted" Python code is not atomic I/O.
   Fix: Treat all shared-output writes as a chunk-fatal commit phase, or write the
   complete chunk into staging files/directories and publish atomically only after every
   write succeeds.
2. **Critical — ENOSPC will not reliably produce exit code 3.** After an artifact write
   fails, the discard re-emission at lines 1226–1230 occurs outside `family_guard`;
   persistent ENOSPC makes that write raise immediately, so execution never increments
   `crashed`, logs `DONE`, or raises the planned `SystemExit(3)`.
   Fix: Stop classifying sink failures such as ENOSPC as family failures; propagate them
   as immediate run-level failures and reserve exit 3 for recoverable per-family
   computation exceptions.
3. **High — continuing after `MemoryError` is an unsafe assumption.** `MemoryError` is
   caught by `family_guard` because it subclasses `Exception`, but logging a traceback,
   retaining the batch, and continuing further alignments/searches all require additional
   memory and can fail unpredictably.
   Fix: Explicitly re-raise `MemoryError` (and other process/resource exhaustion
   exceptions) instead of attempting chunk completion.
4. **High — crash state is inferred from presentation text.**
   `discard_reason.startswith(INTERNAL_ERROR_PREFIX)` makes control flow depend on CSV
   wording, allowing future copy edits or legitimate reasons with the same prefix to
   silently change exit behavior; pinning two strings does not protect every call site.
   Fix: Record an explicit orthogonal `crashed: bool` on `Family`, set only by
   `family_guard`, and count that field while keeping the two-way generated/discarded
   outcome unchanged.
5. **High — rollback is neither complete nor failure-safe.** The cleanup loop stops if any
   `unlink()` raises, potentially leaving other artifacts behind, and it cannot undo shared
   writes that occurred after the artifact block.
   Fix: Stage each family under unique temporary paths, clean every staged path with
   best-effort error aggregation, and rename into final paths only during a successful
   commit.
6. **Medium — the proposed tests miss the failures motivating the release.** The updated
   integration test injects a `renumber_msa` exception before any output I/O; it does not
   exercise artifact ENOSPC, shared-stream partial writes, failed discard re-emission,
   cleanup failure, or the actual console-script exit status.
   Fix: Add fault-injection tests at every artifact/shared-output boundary plus a
   subprocess test asserting exit 3, closed/readable outputs, and `crashed=N` in the log.
7. **Medium — the README exit-status contract is incomplete.** Documenting only 0, 2, and
   3 leaves uncaught runtime/output failures unspecified even though the plan explicitly
   distinguishes them from degraded completion.
   Fix: Document exit 1 as an incomplete/fatal run whose outputs must not be consumed, and
   clarify that exit 3 applies only when chunk processing and output commit completed.

VERDICT: REVISE

### Claude's response

**Accepted — #5, and it is a live bug in code already in the working tree.** The rollback
loop calls `path.unlink(missing_ok=True)` bare; an `OSError` from one unlink replaces the
original exception and skips the remaining unlinks, so the cleanup both fails to clean and
destroys the diagnostic. `contextlib.suppress(OSError)` per unlink. One line, and the best
catch of the round.

**Accepted — #2's finding, rejected its fix.** Verified at `generate_families.py:1226-1230`:
the re-emit is outside `family_guard`, so persistent ENOSPC propagates and the process
exits 1 without ever reaching the counter, the `DONE.` line or `SystemExit(3)`. Codex is
right that the plan's contract was wrong. Its fix — classifying sink failures separately —
adds an exception-taxonomy branch to satisfy a case where the correct outcome is already
happening: a run that cannot write one CSV row is dead, and exit 1 says so. Keeping the
re-emit unguarded is the right behaviour; the plan was wrong to imply exit 3 covers it.
Fixed by documenting exit 1 and stating the boundary explicitly.

**Accepted — #7**, which is the same correction from the documentation side.

**Accepted in part — #6.** Adding a subprocess test that asserts the real console-script
exit status is cheap and tests something no in-process test can. Rejected the full
fault-injection matrix across every I/O boundary: those paths are exercised by the
existing `deterministic_gzip_binary` trap, and a test per boundary buys repetition, not
coverage.

**Rejected — #1's remedy.** The finding is technically correct: two buffered `write()`
calls are not one atomic operation, so a mid-append failure can leave an id in
`converged_families` for a family that ends up discarded. But staging the whole chunk and
publishing atomically is a rewrite of the output layer to protect a window that only opens
when the filesystem is already failing — and when it opens, the re-emit fails too and the
run exits non-zero, so the chunk is never consumed by an operator who checks the status.
The honest fix is the contract, not the architecture: exit 1 means the chunk is incomplete
and must not be read. Documented as a known limitation rather than papered over.

**Rejected — #3.** `MemoryError` from one oversized family is genuinely recoverable: the
family is discarded, its `full_msa` and records are released, and the remaining thousands
complete normally. Re-raising would destroy a chunk that would otherwise finish — the exact
failure mode `family_guard` was built to prevent — in exchange for guarding against a
thrash scenario the exit-3 signal already reports. The asymmetry favours continuing.

**Kept — #4, with a guard added.** Detection via `discard_reason` was locked with the user
in Act 1 against the explicit alternative of a `crashed: bool` field, with the string
coupling stated as the known cost. Codex's concern is fair but its remedy creates two
sources of truth for one event. Mitigation instead: a test asserting that no legitimate
discard reason begins with `INTERNAL_ERROR_PREFIX`, so a future copy edit that would
collide fails loudly.

## Round 2 — Codex

The revision addresses the missing exit-1 contract and console-entry-point coverage, and
the string-marker rejection is defensible given the shared constant plus collision test.
Three material problems remain:

1. **Exit 3 still does not guarantee coherent output.** A shared append can fail
   transiently after an earlier append succeeds, while the subsequent discard re-emit
   succeeds; that reaches exit 3 with contradictory output. The rejection assumes the
   re-emit must also fail, which is not guaranteed.
   Fix: Make any exception after the first shared append a fatal exit-1 condition that
   bypasses `family_guard` and discard re-emission.
2. **Suppressing rollback failures breaks the documented exit-3 contract.** If `unlink()`
   fails, the plan silently ignores it, records the family as discarded, and may exit 3
   while orphan artifacts remain — the exact inconsistency this release claims to fix.
   Fix: Attempt every unlink, collect cleanup failures, then raise a fatal rollback error
   that bypasses family containment and produces exit 1.
3. **The `MemoryError` rejection relies on false lifecycle reasoning.** `family.discard()`
   does not release `seed_msa`, `full_msa`, `hmm`, or `records`; the `Family` remains
   referenced by `active` through the rest of the batch, so continuing does not necessarily
   reclaim the oversized family's memory.
   Fix: Explicitly clear the discarded family's heavyweight fields before continuing after
   `MemoryError`, or re-raise it as fatal.

The broader atomic-publish rewrite can reasonably remain out of scope, but exit 3 must be
structurally prevented whenever output consistency or rollback is uncertain.

VERDICT: REVISE

### Claude's response

**#3 accepted — my Round-1 reasoning was factually wrong.** Verified at
`generate_families.py:635`: `discard()` sets `state`, `discard_reason` and `discard_value`
and nothing else. `seed_msa`, `full_msa`, `hmm`, `records` and `total_checked_sequences`
stay referenced by `active` for the rest of the batch, so the claim that a discarded
family's memory is released was false. The conclusion survives on the corrected premise
only once the premise is made true: `discard()` now clears those fields. That is safe —
`emit_family` returns at the `DISCARDED` branch before reading any of them, `finish()`
no-ops on a discarded family, and the round loop filters on `RUNNING`. It also pays for
itself independently: discarded families currently pin their MSAs for the whole batch,
which is the dominant memory term at large `--batch_size`. `MemoryError` stays catchable,
now for a reason that holds.

**#1 and #2 accepted, unified into one fix rather than two.** Codex is right that the
Round-1 rejection leaned on "the re-emit must also fail too", which a transient ENOSPC or a
per-file quota does not guarantee. Both findings are the same underlying concept: past a
certain point a failure has left output that containment *cannot* make coherent, so
containment must not apply. That is one idea, so it gets one mechanism — a `ChunkCorrupted`
exception that `family_guard` re-raises instead of converting to a discard. Raised from the
shared-append block and from a failed rollback unlink; `main` catches it, logs, and exits 1.
Roughly ten lines, and it makes exit 3 mean *coherent* structurally rather than
probabilistically, which is the whole point of the release.

Still out of scope: staging the chunk and publishing atomically. Codex agrees.

## Round 3 — Codex

All Round 2 findings are adequately addressed:

- Shared-append failures now bypass containment through `ChunkCorrupted`, preserving the
  exit-3 coherence guarantee.
- Rollback attempts every unlink and escalates any cleanup failure to exit 1.
- `discard()` releases heavyweight family state, making the chosen recoverable-`MemoryError`
  policy defensible.
- Tests cover both corruption boundaries and the real console exit status.
- Exit 1 explicitly marks potentially inconsistent output as unsafe.

No new material flaw blocks implementation. When implementing rollback diagnostics, include
failed paths in the `ChunkCorrupted` message, but that is a detail consistent with the plan
rather than a required revision.

VERDICT: APPROVED

### Claude's response

Accepted the one detail: the `ChunkCorrupted` message from a failed rollback names the
paths it could not remove, so the log says which orphans are on disk rather than only that
some are. Folded into the plan.

Converged in 3 rounds. Awaiting the user's sign-off before any implementation code.

## Act 3 — Build

Builder: Codex (`gpt-5.6-sol`), codex-cli 0.146.0, thread `019fd1e6-afef-7201-be67-2bc0e867dac5`.
Spec: `PLAN-2.1.0.md`, frozen. Baseline: `ca5b526`, clean tree. MAX_FIX_ROUNDS=2, used 1.

### Round 1 — Codex build

Implemented all six Approach steps: the `INTERNAL_ERROR_PREFIX` / `EXIT_CRASHED_FAMILIES`
constants, `ChunkCorrupted` with `family_guard` bypass, the shared-append fence, rollback
aggregation naming unremovable paths, `discard()` releasing the heavyweight fields, the
per-batch `crashed` counter, `crashed=` in the `DONE.` line, `SystemExit(3)`, five new
tests, and the CHANGELOG / README / AGENTS.md updates. Reported "deviations: none".

### Claude's verdict — round 1

Proof re-run independently: 47 passed, ruff check and format clean. Logic correct on the
diff read — the fence, the rollback aggregation and the `discard()` release all match the
spec, and `SystemExit(3)` is raised inside the `try` so the `ChunkCorrupted` handler cannot
swallow it. Two problems sent back:

1. **A deviation reported as none.** `test_legitimate_discard_reasons_do_not_collide_with_internal_errors`
   hardcoded six reason strings and asserted those literals did not start with the prefix —
   a fact about its own input. It could never fail for the reason it exists, which the spec
   stated as "a future copy edit that would silently disarm the exit code fails the suite".
   The hardcoded list was also only accidentally complete: the two rep-length reasons are
   returned by `check_rep_length` and reach `discard()` via `*reason`, never as literal
   arguments, so even an AST walk restricted to `discard()` call arguments would have
   missed them.
2. **~200 lines of avoidable diff.** The batch loop had been wrapped in a new nested `try:`
   purely to catch `ChunkCorrupted`, reindenting roughly 100 untouched lines so the whole
   loop read as changed. `main` already had an outer `try: ... finally:` that could carry
   the handler with no reindentation.

### Round 2 — Codex fixes

Both applied. The collision test now parses the module with `ast` and asserts no string
constant anywhere in it starts with the prefix, excluding the prefix's own definition by
node identity — exact, because `family_guard` builds the reason as an f-string
(`ast.JoinedStr`), which never appears as a plain `Constant`. The nested `try` is gone and
the handler moved to the existing outer `try`.

### Claude's verdict — round 2, accepted

Proof re-run: 47 passed, ruff clean. Source diff fell from 279 to 111 lines and the batch
loop is back to showing only the `crashed` accumulation and the reworded re-emit comment.

Mutation-checked the replacement test rather than trusting that it now works: replacing
`"few seed sequences remained"` with `"internal error during seed check"` in the module
source flips the assertion from pass to fail, so the test detects exactly the copy edit it
was written to catch.

Converged: Act 1 grill (4 decisions) -> Act 2 review (3 rounds, APPROVED) -> Act 3 build
(2 rounds). Awaiting the user's diff sign-off; nothing committed.
