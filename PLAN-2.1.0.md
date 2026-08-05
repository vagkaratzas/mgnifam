# Plan: crashed families must fail the chunk loudly (2.1.0.dev0)
_Locked via grill — by Claude + vagkaratzas_

> Scope note for the reviewer: `PLAN.md` and `PLAN-REVIEW-LOG.md` in this repo are the
> historical record of the original port and are read-only for new work (`AGENTS.md:127`).
> This file is the plan under review.

## Goal

A family must leave `emit_family` as **exactly one** of *discarded* or *generated*
(generated being converged or not-converged) — never both, never neither. Auditing that
invariant against `src/mgnifam/generate_families.py` found two holes, both on the failure
path, both already fixed in the working tree (they are context here, not new work):

1. **`zip(strict=True)` inside the uninterrupted region.** The row loop's strict zip ran
   *after* `converged_families` and `successful_clusters` were appended to. A raise there
   put the representative in `successful.txt` and then, through `family_guard`'s re-emit,
   in `discarded.csv` too — both outcomes at once. Fixed by hoisting a length check ahead
   of the first shared append.
2. **Orphan artifacts beside a discard row.** A failure part-way through the four
   per-family writes (`rf/`, `hmm/`, `seed_msa/`, `full_msa/`) left the earlier files on
   disk while the family was recorded as discarded — discarded *and* generated. Usually
   overwritten, because a discard does not consume `provisional_id` and the next success
   reuses the name, but nothing overwrites them when no later family in the chunk
   succeeds — and the likeliest trigger (ENOSPC) guarantees exactly that. Fixed by
   unlinking what was written when the artifact block raises.

The remaining gap is observability. `family_guard` converts any per-family exception into
a discard, so a chunk where families died of ENOSPC or `MemoryError` **exits 0** and is
indistinguishable from a chunk with legitimate discards. The operator is never told the
chunk needs re-running. This plan adds that signal: the chunk still runs to completion,
then exits with a distinct non-zero status.

## Approach

1. **Name the crash marker.** Extract the inline string at `generate_families.py:809` to
   module constants beside `MAX_ROUNDS` (`:59`):

   ```python
   INTERNAL_ERROR_PREFIX = "internal error during "
   EXIT_CRASHED_FAMILIES = 3
   ```

   `family_guard` (`:809`) uses the prefix instead of the literal. The emitted text is
   unchanged — `tests/test_generate_families.py:804` and `:826` pin that exact text and
   must keep passing untouched.

2. **Count crashes in `main`.** Reuse state already recorded; no new `Family` field, no
   `family_guard` signature change, no third `FamilyState`.
   - `:1149-1150` — add `crashed = 0` beside `success_count = 0` / `processed = 0`.
   - `:1231` — beside the existing `processed += len(active)`:

     ```python
     crashed += sum(
         1 for f in active if f.discard_reason.startswith(INTERNAL_ERROR_PREFIX)
     )
     ```

     Accumulated per batch, because `active` is rebound each batch and families are not
     retained across batches. This catches every `family_guard` call site: initialisation,
     round-N model build, round N, the exit branch, and artifact writing.

3. **Report and exit.**
   - `:1243` — add `crashed=%d` to the `DONE.` line so the log states it before exiting.
   - Immediately after that log call — still inside the `try`, so the `finally` at `:1250`
     closes the log handlers, and after the `contextlib.ExitStack` has closed every
     writer:

     ```python
     if crashed:
         raise SystemExit(EXIT_CRASHED_FAMILIES)
     ```

   `SystemExit` is a `BaseException`, and `family_guard` catches `Exception` only
   (deliberately, `:800`), so it cannot be swallowed.

4. **Make exit 3 structurally mean "coherent" — a point of no return.** Past two points, a
   failure has left output that containment cannot repair, so containment must not apply.
   Both get one mechanism, not two:

   ```python
   class ChunkCorrupted(Exception):
       """A failure that has left output no discard row can reconcile.

       `family_guard` contains a per-family failure by recording a discard, which is only
       honest while the family has written nothing a discard contradicts. Past the first
       shared append, or after a rollback that could not remove what it wrote, that is no
       longer true and the chunk is not salvageable -- so this bypasses containment and
       the run exits 1 rather than claiming a coherent exit 3.
       """
   ```

   - `family_guard` re-raises it ahead of the broad catch:

     ```python
     except ChunkCorrupted:
         raise
     except Exception:
         ...record the discard as today...
     ```

   - `emit_family` raises it from the shared-append block (`converged_families` through the
     row loop) — `except Exception as error: raise ChunkCorrupted(...) from error`. This is
     the window where `converged_families` can hold an id for a family that a re-emit then
     records as discarded.
   - `emit_family` raises it when a rollback unlink fails. Attempt **every** unlink first,
     collecting failures, so one bad path does not strand the rest; then raise
     `ChunkCorrupted` if any failed, naming the paths it could not remove so the log says
     which orphans are on disk, and chaining the original error. This also fixes a live bug
     in the working tree: the current loop calls `path.unlink(missing_ok=True)` bare, so an
     `OSError` replaces the original exception and skips the remaining unlinks.
   - `main` catches `ChunkCorrupted` around the batch loop, logs it, and exits 1 — the
     `finally` at `:1250` still closes the log handlers.

   Rejected as disproportionate: staging the chunk and publishing atomically. Codex agreed
   it can stay out of scope once exit 3 is structurally gated.

4b. **Do not guard the discard re-emit at `:1226-1230`.** It stays outside `family_guard`:
   a run that cannot append one row to `discarded.csv` has a dead output sink and nothing
   left to record the failure *with*, so it propagates and the process exits 1. Add a
   comment saying so — the current one explains only why the re-emit exists.

4c. **`discard()` must release the family's heavyweight fields.** Currently it sets
   `state`, `discard_reason` and `discard_value` only (`:635`), so a discarded family keeps
   `seed_msa`, `full_msa`, `hmm`, `records` and `total_checked_sequences` referenced by
   `active` for the rest of the batch. Clear them:

   ```python
   self.seed_msa = None
   self.full_msa = None
   self.hmm = None
   self.records = []
   self.total_checked_sequences = set()
   ```

   Safe: `emit_family` returns at the `DISCARDED` branch before reading any of them,
   `finish()` no-ops on a discarded family, and the round loop filters on `RUNNING`. Two
   payoffs — it is what makes recovering from a `MemoryError` actually reclaim the
   oversized family's memory, and independently it stops discarded families pinning their
   MSAs for the whole batch, which is the dominant memory term at large `--batch_size`.

5. **Tests** (`tests/test_generate_families.py`).
   - `test_one_failing_family_is_discarded_and_the_chunk_survives` (`:761`) — its premise
     now includes the exit code. Wrap the `run_pipeline` call in
     `pytest.raises(SystemExit)`, assert `.value.code == gf.EXIT_CRASHED_FAMILIES`, and
     keep every existing assertion: the chunk must still write complete, contiguous
     output. This is the test proving "finish the chunk" beat "abort on first crash".
   - New: a clean run exits 0 — `gf.main(...)` returns without `SystemExit`, and
     `chunk_discarded.csv` holds only legitimate reasons.
   - New: **the real console-script exit status.** A `subprocess.run(["mgnifam",
     "generate_families", ...])` on a chunk rigged to crash one family, asserting
     `returncode == 3`, `crashed=1` in `chunk.log`, and complete output on disk. In-process
     tests assert the `SystemExit` value; only a subprocess proves it survives `cli.main`
     and the console-script wrapper. The repo already runs `mgnifam` under `subprocess` at
     `:1195`, so the pattern exists.
   - New: **no legitimate discard reason collides with `INTERNAL_ERROR_PREFIX`.** Assert it
     over every literal reason `advance` and `finish` pass to `discard()`. This is what
     keeps detection-by-string honest: a future copy edit that would silently disarm the
     exit code fails the suite instead.
   - New: **a failure in the shared-append phase exits 1, not 3.** Trap a write on
     `successful_clusters` after `converged_families` has been appended; assert
     `ChunkCorrupted` is not converted to a discard and the run exits 1. This is the test
     that proves exit 3 means coherent.
   - New: **a failing rollback unlink exits 1 and does not mask the original error.**
     Covers both halves of the live bug — every unlink is attempted, and the original
     exception is chained rather than replaced.
   - The two tests added alongside the `emit_family` fixes stay as-is.

6. **Docs — `AGENTS.md:112` requires these move together.**
   - `CHANGELOG.md` — fill the existing empty `## [2.1.0.dev0] - unreleased` heading
     (`:18`). `### Fixed` for the two `emit_family` holes; `### Changed` for the new exit
     status, flagged as a behaviour change for any caller that ignored the exit code.
   - `README.md` — a short "Exit status" block under *Outputs* (`:142`). Exit codes are
     currently undocumented, so this is new documented behaviour. The contract must state
     what each code says about the *output*, not just about the run:

     | Code | Meaning |
     |------|---------|
     | `0` | Chunk completed. Every family landed on exactly one side of the split. Output is complete and safe to consume. |
     | `1` | Fatal: the run died before finishing. **Output is incomplete and must not be consumed** — re-run the chunk. This is what a dead output sink (ENOSPC, EIO) produces, because the discard re-emit cannot record its own failure. |
     | `2` | Usage error from `argparse`. Nothing ran. |
     | `3` | Chunk completed, but one or more families died of an internal error and were recorded as discards. Output is complete and self-consistent, but those clusters produced no family — re-run the chunk once the cause is fixed, or accept the loss. |

     The `1`-vs-`3` line is the one that matters operationally: `3` means *degraded but
     coherent*, `1` means *do not read this*.
   - `AGENTS.md` — extend the invariant paragraph at `:44-49` (already states "Discarded
     families must never appear in `converged_families`") with the two-way split and the
     fact that a crashed family is a discard that also fails the run.

## Key decisions & tradeoffs

- **Finish the chunk, exit at the end — not abort at the crash.** Rejected re-raising in
  `family_guard`. That guard exists precisely so one family's failure does not destroy the
  hours of search already completed by every other family in the batch
  (`family_guard` docstring `:794`; `test_one_failing_family_is_discarded_and_the_chunk_survives`).
  The requirement is *observability* — the operator must know to re-run — and an exit code
  delivers that without discarding work. Cost accepted: a chunk you will re-run burns its
  remaining CPU. Gain: nothing in-flight is lost, and salvaging N-1 good families stays an
  option the operator can take. Also rejected: stopping between batches, which bounds the
  waste but yields a partially-processed chunk that downstream could mistake for a
  complete one.
- **Detect via `discard_reason` prefix, not new state.** The information is already
  recorded by `family_guard`. A `crashed: bool` on `Family` or a counter threaded through
  `family_guard`'s five call sites would both be second sources of truth for one event.
  Accepted coupling: the exit signal now depends on an output-visible CSV string. The
  shared constant is what keeps that honest.
- **Exit code 3, not 1.** Avoids collision with `1` (uncaught exception) and `2` (argparse
  usage error), so "completed but degraded, output salvageable" stays distinguishable from
  "died". Cost: a documented contract that must be kept.
- **A crashed family stays `DISCARDED`; no `CRASHED` state.** A third `FamilyState` would
  make the per-family outcome three-way, which is the opposite of the goal. "Crashed" is a
  run-level annotation on top of a discard, not a family outcome.
- **Keep the Hole-1 length check even though it is unreachable via the real
  `renumber_msa`** (which builds the `TextMSA` from a list of `TextSequence`, so pyhmmer
  guarantees one name per row). Both paths exit 3; the check is what decides whether the
  exiting chunk is self-consistent or has the same representative in `successful.txt` and
  `discarded.csv`.

## Risks / open questions

- **The shared-append phase is still not atomic — but it can no longer produce exit 3.**
  Five buffered `write()` calls are not one operation, so an I/O failure between them can
  still leave an id in `converged_families` for a family that ends up discarded. Making
  that impossible means staging the chunk and publishing atomically, which stays out of
  scope. What changed is that the window is now *fenced*: any failure there raises
  `ChunkCorrupted`, bypasses containment, and exits 1. So incoherent output is still
  reachable, but it is never labelled coherent.
- **`MemoryError` is treated as recoverable, and 4c is what makes that true.** `family_guard`
  catches it and the family is discarded; only because `discard()` now clears the heavy
  fields does the oversized family's memory actually get reclaimed, letting the rest of the
  batch finish. Re-raising instead would destroy a chunk that would otherwise complete —
  the failure mode `family_guard` exists to prevent. Revisit if a real run shows the batch
  failing to recover.
- The crash counter keys on a string that is also public CSV output. Changing the discard
  reason text silently disarms the exit signal. Mitigated by the shared constant, by
  `:804`/`:826` pinning the text, and by the new no-collision test — but it is a real
  coupling, chosen over a second source of truth on `Family`.
- Callers that currently ignore `mgnifam`'s exit status see no change in behaviour; callers
  that check it will start failing on chunks that previously "passed". That is intended,
  but it is a breaking change for a pipeline with `errorStrategy 'terminate'`.
- Under ENOSPC the crash count inflates — every subsequent family fails for the same root
  cause — and the run most likely ends at exit 1 rather than 3, because the re-emit cannot
  write either. Harmless for the signal; both codes mean re-run.
- Not addressed: whether a chunk that crashed should leave a marker file so a downstream
  step can skip it without consulting the exit code.

## Out of scope

- Abort-on-first-crash; stopping between batches.
- Retrying a crashed family.
- A `CRASHED` `FamilyState`.
- Staging + atomic publish of the chunk's output. Reviewed in Round 1; rejected as
  disproportionate — see *Risks*.
- An exception taxonomy separating sink failures (ENOSPC/EIO) from computation failures.
  The correct outcome for a dead sink already happens: the run exits 1.
- `pyproject.toml` / `uv.lock` — already at `2.1.0.dev0`.
- Any change to scientific output. The counter is read-only over existing state and the
  exit runs after all writers close, so byte-identical output is preserved for any run
  that completes cleanly.
