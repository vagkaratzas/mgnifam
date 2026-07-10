# Plan Review Log: Port & harden `generate_families.py` into the standalone `mgnifam` repo

Act 1 (grill) complete — plan locked with the user. MAX_ROUNDS=5.

Reviewer model: `gpt-5.6-sol` (pinned in `~/.codex/config.toml`, `model_reasoning_effort = "xhigh"`) — codex-cli `0.144.0-alpha.4`.
Plan file: `PLAN.md`. Legacy source under review sits at `reference/legacy_generate_families.py`.

## Act 1 — decisions the grill settled

1. **Search strategy** — stream targets from a `SequenceFile` *and* batch families
   into round-major `hmmsearch` waves. (Rejected: streaming without batching,
   which drops to one worker thread; and keeping the prefetched
   `DigitalSequenceBlock`, which is the memory the user wants gone.)
2. **Small-MSA discard cutoff** — discard when `<= 2` sequences remain after
   `run_pytrimal_reps`. (Rejected: `<= 1`, the bare `eslEMEM` guard.)
3. **Item 2 clarified by the user** — this is *not* the existing gap-occupancy
   `clip_ends()`. It is a **new** function applied to the `.sto` alignment
   returned by `run_hmmalign` in the `# main strategy continue` section, which
   trims the leading and trailing columns where the `#=GC RF` line holds `.`
   (envelope overhang). Interior `.` between `x` match positions are kept. Not
   applied to the final full-MSA `hmmalign`.
4. **Git workflow** — one feature, one commit, no push; stay on branch `dev`.

## Round 1 — Codex (`gpt-5.6-sol`, xhigh)

The exit-search HMM trace is correct: on convergence it reuses `hmm_k`, and on the fourth iteration it reuses `hmm_3` ([legacy loop](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:507)). The search is computationally redundant, but the proposed cache transformation is not yet behaviorally identical.

1. **Cached hits lose native ordering.** Flattening into `(name, tlen, evalue, score, env_from, env_to)` and globally sorting can reorder hits/domains, changing MSA row 0, representative-length checks, trimming, and later HMMs; query-parallel execution does not create the alleged per-query “thread-merge order,” and `TopHits` is already ranked. ([plan](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:168), [PyHMMER results API](https://pyhmmer.readthedocs.io/en/stable/api/plan7/index.html)) — **Fix:** Preserve explicit `(hit_rank, domain_rank)` from `TopHits`, do not re-sort, and oracle-test cached exit output against the legacy second search for convergence rounds 1–3 and forced round 4.

2. **The determinism sort is not “free” or over an “already-short” list.** The plan itself allows roughly one million records per family, where `sorted()` adds O(N log N) work and another large list allocation. — **Fix:** Retain HMMER’s native deterministic ranking and remove the record sort unless a reproducible failing tie case proves it necessary.

3. **`itertools.groupby` changes cluster semantics.** Legacy processing gathers every row matching the first representative, even when those rows are non-contiguous ([legacy grouping](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:142)); `A, B, A` becomes three families under `groupby`, not two. — **Fix:** Accumulate members in an insertion-ordered mapping, or validate and reject non-contiguous representatives before processing.

4. **Family numbering is only preserved if all side effects are deferred.** A later family may converge in round 1 while an earlier family succeeds in round 3; writing convergence, discard, or metadata records from `advance()` would produce wave order and premature IDs despite final `emit()` being ordered. — **Fix:** Make state transitions side-effect-free, buffer every output event, and assign successful IDs only during cluster-order emission.

5. **The numbering test is too weak.** Checking only that IDs are `1..n` contiguous would accept IDs attached to the wrong representatives. — **Fix:** Compare the complete representative→family-ID mapping and every per-family artifact for `batch_size=1`, `batch_size>1`, and mixed early/late success and discard cases.

6. **Gzipped FASTA input is incompatible with SSI random access.** The declared fixtures are `.fa.gz`, while verification silently switches to `.fa` ([layout](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:92), [verification](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:290)); Easel explicitly rejects positioning and SSI on gzipped sequence files. ([Easel source](https://sources.debian.org/src/hmmer/3.4%2Bdfsg-3/easel/esl_sqio_ascii.c)) — **Fix:** Require an uncompressed FASTA with an actionable error, or explicitly stage and index a decompressed copy; make fixtures and documentation agree.

7. **`record_offset` alone is safe only in the narrow full-record case.** Easel’s `PositionByKey` does indeed request only `record_offset`, so omitting `data_offset` and length is valid for uncompressed complete-sequence fetches; however, `SSIWriter.add_file` requires an integer format code, not the `SequenceFile.format` string implied by `fasta_format`. ([SSIWriter API](https://pyhmmer.readthedocs.io/en/stable/api/easel/misc.html)) — **Fix:** Use a documented numeric format code or `0`, state that subsequence fetching is unsupported, and test the exact pyhmmer 0.12.1 constructor/index path.

8. **Existing SSI files are trusted blindly.** A stale, truncated, malicious, or same-basename index can return the wrong valid record without obvious failure, and “skip if exists” has a check-then-create race between jobs. — **Fix:** Build under a lock into a unique temporary file, atomically rename it, and validate a sidecar containing canonical FASTA path, size, mtime/content fingerprint, format, and completed-index marker.

9. **The default index location duplicates billion-record indexing across pipeline tasks.** Separate Nextflow task directories will each build their own index unless `--fasta_index` is wired upstream. — **Fix:** Make production runs require one pipeline-built immutable shared index, leaving local auto-build only as an explicitly documented small-data convenience.

10. **Missing and duplicate FASTA names change behavior.** Legacy initial MSA silently skips missing members ([legacy lookup](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:197)); direct indexed access raises `KeyError`, while duplicate primary names may make SSI construction fail even though the legacy dictionary retained one value. — **Fix:** Define and test missing-member and duplicate-name policy explicitly, preserving legacy behavior where required.

11. **The small-MSA guard does not cover every `hmmbuild`.** It is placed only after `run_pytrimal_reps`, so a tiny initial cluster reaches the first build without the promised protection. — **Fix:** Enforce the `<=2` guard immediately before every build, including the initial and final hand-architecture builds.

12. **HMM files cannot currently be byte-deterministic.** `Builder` stamps current creation time, and its command line defaults to `sys.argv`; both are serialized into the HMM. ([PyHMMER HMM metadata](https://pyhmmer.readthedocs.io/en/stable/api/plan7/hmms.html)) — **Fix:** Set `creation_time` and `command_line` to fixed documented values before every HMM write and assert them in tests.

13. **Logs contradict the all-files determinism test.** The legacy logger includes `asctime`, and `log_time` records wall-clock durations ([logger](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:22)); these cannot be byte-identical across runs. — **Fix:** Restrict the byte-determinism contract to scientific artifacts and test logs separately, or remove all time-dependent fields from persisted logs.

14. **The gzip header claim is false.** With a named `fileobj`, `GzipFile` defaults the header filename from `fileobj.name`; `mtime=0` removes only the timestamp. ([Python gzip documentation](https://docs.python.org/3.13/library/gzip.html)) — **Fix:** Pass `filename=""`, `mtime=0`, and an explicit mode to every `GzipFile`.

15. **The FAMSA statement misdescribes legacy behavior.** Legacy calls `pyfamsa.Aligner()` and ignores `cpus` ([legacy initial MSA](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:217)); changing to `Aligner(threads=cpus)` is not something being “kept,” and thread-count independence is merely asserted. — **Fix:** Document this as a behavior/resource change and compare outputs across `threads=1`, requested CPUs, and repeated runs before relying on it.

16. **Shared temporary files and logging remain race-prone.** Fixed `tmp/seed_msa.sto`, `tmp/full_msa.sto`, and global `execution.log` collide when chunks share a working directory ([temporary paths](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:116), [log move](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:635)). — **Fix:** Use a per-invocation `TemporaryDirectory`, write directly to the chunk log, and atomically replace completed outputs.

17. **`chunk_num` remains a path-traversal input.** It is interpolated directly into paths, so values containing `../` can escape output directories. — **Fix:** Restrict chunk identifiers to a safe basename pattern such as `[A-Za-z0-9._-]+` and reject path separators.

18. **The default batch size is not safely bounded.** `max(4*cpus, 32)` can retain tens or hundreds of millions of records and launches many concurrent full-file readers; restoring CPU occupancy does not prove production storage will avoid I/O saturation. — **Fix:** Default near `cpus`, add an explicit memory budget/backpressure policy, and benchmark peak RSS plus throughput on representative storage.

19. **The legacy-diff expectation is scientifically wrong.** Envelope clipping occurs before later builds, so it can alter subsequent HMMs, recruitment, convergence, discard decisions, and family counts—not merely the final Stockholm columns; dependency upgrades and the new tiny-MSA rule add further allowed differences. — **Fix:** Compare per-round traces and define explicit expected deltas instead of requiring equal membership/count with “only clipping” differences.

20. **The verification rules out a simpler bounded-memory design without testing it.** Chunked `DigitalSequenceBlock` searches with fixed global `Z` and `TopHits.merge`, or external `hmmsearch --cpu`, could retain family-major semantics and target parallelism while bounding memory; monkeypatching every `read_block` call to fail rejects this by fiat. — **Fix:** Prototype and benchmark bounded target blocks against wave scheduling before committing to the substantially more complex state machine.

21. **Long-run observability is inadequate.** Post-call timing logs provide no heartbeat, current batch/round, completed-query count, hit count, RSS, input/index fingerprint, or recoverable checkpoint during multi-hour database passes. — **Fix:** Emit periodic structured operational metrics and atomically checkpoint each completed batch while keeping those logs outside the deterministic-artifact contract.

VERDICT: REVISE
### Claude's response — Round 1

Before revising I ran four probe scripts against `tests/fixtures/mgnifams_input_small.fa`
(50 000 seqs) on the pinned versions. Findings F1–F8 are now recorded in `PLAN.md`.
Codex's critique was strong: it found two things I had asserted without evidence and
that turned out to be **false**, and one "simpler alternative" that turned out to be
**wrong**. Both directions were settled by measurement, not argument.

**Accepted (16).**
- **#1, #2 — drop the tie-break sort.** Correct on both counts. Measured: streamed
  `TopHits` order is identical across `cpus = 1, 2, 4, 8` (F4), because each query is
  scanned by exactly one `Pipeline` in file order. The sort was unnecessary *and*
  O(N log N) over ~10⁶ records. Removed; native order preserved end to end.
- **#3 — `itertools.groupby` is wrong.** Legacy gathers *all* rows for a representative
  regardless of contiguity; `A, B, A` would become three families. Switched to an
  insertion-ordered dict.
- **#4 — defer all side effects.** Real. A family converging in round 1 would write its
  `converged_families` line before an earlier family that succeeds in round 3. All state
  transitions are now pure; every output event is buffered on `Family` and flushed in
  cluster order by `emit()`.
- **#5 — numbering test too weak.** Added: compare the full `representative → family_id`
  mapping, and require byte-identical outputs between `batch_size=1` and `batch_size=64`.
- **#6 — gzipped FASTA vs SSI.** Correct and I had the contradiction in writing. The CLI
  now rejects a gzipped `--fasta_file` with an actionable error; `conftest.py`
  decompresses the `.fa.gz` fixtures into `tmp_path`.
- **#7 — `add_file` format code.** Confirmed: it takes an integer.
  `SequenceFile._FORMATS["fasta"]` works (F1). Noted as private API with a regression test.
- **#8 — SSI trust/race**, partially. Accepted the real defect: build to
  `<ssi>.tmp.<pid>` + `os.replace()`, and rebuild when the index is older than the FASTA.
  *Rejected* the sidecar fingerprint (path, size, mtime, content hash, completion marker) —
  that is filesystem-integrity engineering this script does not need.
- **#9 — index rebuilt per Nextflow task.** Documented `--fasta_index` as the production
  path: build once upstream, share across chunks; auto-build is the small-data convenience.
- **#10 — missing/duplicate names.** Measured (F2): a duplicate primary key raises
  `ValueError: Index contains duplicate keys.` at `SSIWriter.close()` plus a secondary
  `EaselError` from `__dealloc__`. Now detected during the scan → fail fast with the
  offending names. Missing members remain silently skipped in `run_initial_msa` only
  (legacy behaviour); a miss anywhere else stays a hard error, as in legacy.
- **#11 — guard every `hmmbuild`**, with modification. Applied `<= 2` where the user asked
  (after `run_pytrimal_reps`) and `< 2` before the initial and hand builds. Using `<= 2`
  on the initial cluster would newly discard every two-member cluster — a scientific
  change nobody requested. The `< 2` guard still converts `eslEMEM (status 5)` into a
  clean discard.
- **#12 — HMMs are not byte-deterministic.** Confirmed empirically (F6): `Builder` stamps
  `creation_time = datetime.now()` and `command_line = sys.argv`, both serialised. Fixed
  values are now set before every write, and the bytes are stable across writes. My
  determinism claim was false as written; good catch.
- **#13 — logs contradict the determinism test.** Right. The byte-determinism contract now
  covers scientific artifacts only; `logs/` is explicitly excluded.
- **#14 — the gzip header claim is false.** Right: `GzipFile` copies `fileobj.name` into
  the header. All writes now pass `filename="", mtime=0, mode="wb"`.
- **#15 — FAMSA wording.** Right that legacy calls `Aligner()` (all cores, ignoring
  `--cpus`), so `threads=cpus` is a change, not a preservation. Measured (F7): alignments
  identical across `threads = 1, 2, 4, 8` and repeatable. Documented as a deliberate
  resource-behaviour change.
- **#16 — shared tmp files and `execution.log`.** Adopted: per-invocation
  `TemporaryDirectory()`, logger writes straight to `logs/<chunk>.txt`, no `shutil.move`.
- **#17 — `chunk_num` path traversal.** Adopted: `[A-Za-z0-9._-]+`. Trust boundary.
- **#19 — legacy-diff expectation.** Right, and worse than Codex knew (see F3). Demoted
  from an assertion to a delta report with each divergence explained.
- **#21 — observability**, partially. Adopted per-(batch, round) structured progress
  logging. *Rejected* checkpoint/resume: Nextflow already resumes at chunk granularity.

**Rejected, with measurements (2).**
- **#20 — "prototype the simpler bounded-memory design: chunked `DigitalSequenceBlock`
  with fixed global `Z` and `TopHits.merge`."** I did prototype it. It is **not
  equivalent**: reading the DB in 5 000-sequence blocks with `Z=domZ=50000` and merging
  lost hits — 26 vs 27 and 54 vs 56 — and ran slower (0.32 s vs 0.17 s). Each block gets
  its own `Pipeline`, which applies a **per-chunk `domZ`** at the inclusion threshold, so
  the merged `TopHits` disagrees with a single pass. Rejected on evidence (F5).
- **#18 — "the default batch size launches many concurrent full-file readers."** The I/O
  half is wrong: `hmmsearch` spawns exactly `cpus` workers regardless of query count, and
  each query is read by one worker, so total file reads equal total queries whatever the
  batch size. The *memory* half is fair — default lowered to `max(2 * cpus, 16)`. Note
  neither default bounds the dominant term, a single family's ~10⁶ hit records; recorded
  as a known risk.

**What Codex could not have known — and what changes the project.**
Chasing #20 surfaced the real defect (F3). pyhmmer auto-selects `parallel="targets"`
whenever `n_queries < cpus`, which is **every legacy call**, since legacy searches one HMM
at a time. Each worker then runs its own `Pipeline` over a chunk with a **per-chunk
`domZ`**, so domains clear the inclusion threshold too easily:

| run | hits per query |
|---|---|
| prefetch block, `cpus=1` | 26 / 19 / 55 |
| prefetch block, `cpus=4` (auto → `targets`) | **27 / 19 / 56** |
| prefetch block, `cpus=4`, `parallel="queries"` | 26 / 19 / 55 |
| stream `SequenceFile`, `cpus=1 / 3 / 8` | 26 / 19 / 55 |

**The legacy pipeline's scientific output depends on `task.cpus`.** Streaming forces
`parallel="queries"` and reproduces the `cpus=1` answer at every core count. Eliminating
the in-memory FASTA also eliminates a live reproducibility bug. Consequence, now flagged
as the top risk: the port will *not* reproduce legacy production output — it matches
legacy at `cpus=1`. A `--cpus`-invariance test enforces this going forward.

## Round 2 — Codex

The revision addresses most prior findings, but it is not ready to implement. The design-level pass exposed two new boundary failures: `Family.records` retains pyhmmer internals, and the SSI is not reliably bound to its FASTA.

What now holds:

- The exit-branch `hmmsearch` is redundant. At convergence it reuses `hmm_k`; after round 3 it reuses `hmm_3`. Only the Python length filter changes, so caching raw domain records is semantically valid.
- The basic round-major ordering is sound: ordered clusters, ordered query results, deferred emission, and a global successful-family counter can preserve numbering.
- `record_offset` alone is sufficient for whole-record SSI fetches; `data_offset` and `record_length` are optional optimizations not needed here.
- I accept rejecting blocked target merging as the implementation path. I do not accept the stated explanation.

Material findings:

1. **F3’s causal explanation is unsupported and likely wrong.** The final `domZ` is identical across the compared runs, pyhmmer 0.11.1 merges automatic `Z`/`domZ` and then re-thresholds flags without deleting stored hits, while legacy iterates every hit/domain rather than `.reported` or `.included` ([legacy extraction](/home/vangelis/Desktop/Projects/mgnifam/reference/legacy_generate_families.py:265), [0.11.1 merge](/home/vangelis/miniconda3/pkgs/pyhmmer-0.11.1-py313h366bbf7_2/lib/python3.13/site-packages/pyhmmer/plan7.pyx:8552)). Locally retained but globally unreported candidates are a more plausible explanation.

   Fix: Record `len`, reported/included counts and flags, exact recruited `(name, tlen, env_from, env_to)` records, and pinned-0.11.1 end-to-end artifacts at `cpus=1/4` before claiming the production pipeline is “inflated.”

2. **F5 cannot be a “per-chunk domZ” defect as written.** Passing `domZ=50000` manually gives every chunk the same value; moreover the single-pass automatic values shown are 26/19/54, so forcing 50000 changes the scientific parameter and is not an equivalent baseline.

   Fix: Keep the blocked design rejected as empirically non-equivalent, but remove the false causal claim or rerun block and single-pass searches with identical `Z` and `domZ` settings.

3. **The worker-count correction is itself incorrect.** Legacy target-parallel mode uses up to `cpus` target workers, but the new query-parallel dispatcher uses `min(cpus, n_queries)` when the HMM list length is known ([dispatcher](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer/hmmer/_base.py:367)). The claim about total database scans equalling query count is correct; “exactly cpus workers regardless of query count” is not.

   Fix: Document `min(cpus, active_queries)` and default nearer `2 * cpus` without the 16-family floor, with a configurable hard cap and peak-RSS benchmark.

4. **The wave pseudocode does not map HMMs to families correctly.** `i` is undefined at [PLAN.md:217](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:217), and ordinary `zip` could silently truncate.

   Fix: Use `for f, hmm, hits in zip(running, hmms, hmmsearch(...), strict=True): f.hmm, f.qlen = hmm, hits.query.M`.

5. **`f.records` is not a compact cache.** A pyhmmer `Domain` retains its parent `Hit`, which retains the entire `TopHits`; batching therefore holds full HMMER result graphs and query objects, undermining the proposed memory cap.

   Fix: Immediately extract immutable compact records containing only `name`, `tlen`, `env_from`, and `env_to`, preserving nested iteration order.

6. **The seed Stockholm ID will change incorrectly.** The plan names seed MSAs after the representative, but legacy’s final hand build sets the seed name to `<chunk>_<family_id>` before writing; that name is serialized as `#=GF ID`.

   Fix: Rename the successful seed MSA to `<chunk>_<family_id>` immediately before emission and assert its `#=GF ID` in the legacy comparison.

7. **`converged_families` numbering remains underspecified.** Legacy writes the provisional ID immediately on convergence, before membership/length checks; a subsequently discarded family causes that ID to be reused, potentially producing duplicates. A simple “successful-family rank” cannot reproduce this behavior.

   Fix: Explicitly preserve or deliberately change converge-then-discard semantics, with a fixture containing converged-discarded-successful families and an exact expected `converged_families` file.

8. **Duplicate detection reintroduces O(database) Python memory.** Detecting duplicates during a billion-record scan requires retaining every seen key, potentially tens of gigabytes, even though production indexing is supposed to avoid database-sized Python state.

   Fix: Let the SSI builder perform disk-backed duplicate detection and translate its close error, or make large-database indexing an upstream external step and enforce a size limit on auto-build.

9. **The SSI provenance rejection is unsafe, not merely “filesystem integrity.”** The default index uses only the FASTA basename, and mtime cannot distinguish another newer FASTA with the same basename; an incorrect `--fasta_index` is attached without validating that its offsets describe the supplied file, causing silent wrong-sequence fetches.

   Fix: Bind the index to resolved path, size, `mtime_ns`, and a build-time fingerprint, and assert every fetched sequence’s name equals the requested key.

10. **`<ssi>.tmp.<pid>` does not fully guard the index race.** PIDs collide across hosts on shared storage, and a predictable non-exclusive temporary name permits concurrent writers or symlink replacement.

   Fix: Create the temporary index with an exclusive random file in the destination directory and use a build lock or require the shared production index to pre-exist.

11. **HMM byte determinism is locale-sensitive.** Assigning a fixed `datetime` still invokes locale-dependent `%a` and `%b` formatting inside pyhmmer, so HMM bytes can differ with `LC_TIME`.

   Fix: Set `hmm.creation_time = None` or force and test `LC_TIME=C`.

12. **The determinism test avoids stale-output behavior.** Running in fresh output directories cannot detect old `<chunk>_<id>` files left when a rerun produces fewer successful families.

   Fix: Fail on an existing chunk namespace or emit into a staging directory and atomically replace the complete chunk output.

13. **One log line per batch/round is not a heartbeat.** A single database pass may run for hours, so this only records boundaries, not progress; checkpointing can remain out of scope.

   Fix: Add a time-throttled heartbeat or callback reporting completed queries, active batch/round, elapsed time, and RSS.

14. **Input validation is incomplete.** `cpus<=0`, `batch_size<=0`, invalid percentage ranges, contradictory length bounds, blank identifiers, and malformed cluster rows currently have no defined failure behavior.

   Fix: Validate numeric ranges and require exactly two non-empty TSV fields per cluster row before any index or output is created.

15. **A material scientific decision remains explicitly open.** [PLAN.md:345](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:345) still asks the user to confirm choosing single-pipeline semantics over legacy production behavior.

   Fix: Record the explicit decision in the plan and describe the legacy-compatibility consequence as settled rather than an open question.

VERDICT: REVISE
### Claude's response — Round 2

Codex was right about the thing that mattered most: **my F3 causal explanation was wrong.**
I ran its hypothesis and it holds.

**The corrected F3.** Measured under both pyhmmer 0.11.1 (the legacy pin) and 0.12.1:

| cpus | `len(TopHits)` | `.reported` | `.included` | `Z` | `domZ` |
|---|---|---|---|---|---|
| 1 | 26 / 19 / **55** | 26 / 19 / 54 | 26 / 19 / 54 | 50000 | 26 / 19 / 54 |
| 4 | **27** / 19 / **56** | 26 / 19 / 54 | 26 / 19 / 54 | 50000 | 26 / 19 / 54 |

`Z` and `domZ` are identical, exactly as Codex said. Target-parallel merge concatenates each
worker's locally-stored hits and re-thresholds the flags without deleting entries, so
`len(TopHits)` grows while `.reported`/`.included` do not. Legacy iterates the **raw**
`TopHits`, so it recruits sub-threshold hits — one even at `cpus=1` (55 stored vs 54
reported) — and the number varies with `--cpus`. Corrected in `PLAN.md`; the "per-chunk
`domZ`" claim is withdrawn.

This also surfaces a **latent scientific bug in the original**: `--recruit_evalue_cutoff`
does not actually bound recruitment. Preserved (raw iteration = legacy `cpus=1`), flagged
for the user, **not** silently fixed.

**F3b, new.** End-to-end at `cpus=1` vs `cpus=4` under the pinned legacy venv: every output
identical except the HMM `DATE`/`COM` lines. The extra stored hit is removed downstream by
`env_length >= 0.9 * qlen`. The hazard is real at recruitment, latent on this fixture. My
round-1 phrasing ("the legacy pipeline's output depends on task.cpus") overstated the blast
radius; corrected.

**#2 — F5.** Conceded. My blocked-vs-single comparison forced `Z=domZ=50000` against a
baseline whose automatic `domZ` was 26/19/54, so it was not like-for-like. The causal claim
is withdrawn. Blocking stays rejected, but now on the honest ground that streaming already
delivers O(1) memory *and* `cpus=1` semantics, so blocking buys nothing.

**#3 — worker count.** Conceded, verified in `_base.py`:
`self.cpus = 1 if hint == 0 else min(cpus, hint) if hint > 0 else cpus`. My round-1 rebuttal
was wrong. `--batch_size` must stay ≥ `cpus` or cores idle. Recorded as F9. Kept the 16-family
floor — it costs nothing at low `cpus` and keeps small chunks from degenerating.

**Accepted (11).**
- **#4** — `zip(..., strict=True)`, and the pseudocode now binds `hmm` in the loop; `i` was undefined.
- **#5** — verified: a `Domain` pins its `Hit`, which pins the `TopHits`. Records are copied
  into plain `(name, tlen, env_from, env_to)` tuples at extraction. Recorded as F10, with a
  `gc.get_referents` test.
- **#6** — real and I had missed it. Legacy's `run_hmmbuild` sets
  `seed_msa.name = f"{chunk}_{iteration}"`, which serialises as `#=GF ID` in the seed
  Stockholm. Naming the MSA after the representative would have silently changed every seed
  file. The rename now happens in `emit()` once the family id is known.
- **#7** — the best find of the round. Legacy writes `converged_families` **at convergence,
  before** the membership/length checks, then decrements `iteration` on discard — so the id is
  reused and the file can contain duplicates or ids belonging to another family. A
  "successful-family rank" alone cannot reproduce that. `emit()` now writes the *provisional*
  id on convergence regardless of a later discard, and a converged→discarded→successful
  fixture pins it.
- **#8** — conceded. A Python `set` of a billion keys is precisely the database-sized state
  this port removes. `SSIWriter` already detects duplicates in its own structure (F2); we
  catch its `ValueError`, re-raise with guidance, and unlink the partial temp index in a
  `finally` (its `__dealloc__` otherwise throws a second `EaselError`).
- **#9** — accepted in substance, not in mechanism. Every `.indexed[name]` fetch now asserts
  `fetched.name == name`, which catches a stale or mismatched index on first use in O(1).
  Still rejecting the fingerprint sidecar.
- **#10** — accepted: `tempfile.mkstemp(dir=...)` for an exclusive, unpredictable temp name,
  then `os.replace`.
- **#11** — accepted, and improved by measurement: `hmm.creation_time = None` and
  `hmm.command_line = None` omit the `DATE`/`COM` lines *entirely*, which is strictly better
  than a fixed date formatted through locale-dependent `%a`/`%b`. Verified byte-stable.
- **#12** — accepted: stale `<chunk>_*` per-family artifacts are cleared at start-up, and a
  rerun-into-the-same-directory test covers it.
- **#13** — accepted: a time-throttled heartbeat via the `hmmsearch` `callback=`, not one line
  per round boundary.
- **#14** — accepted in full: numeric ranges, ordering of the length bounds, non-empty chunk id,
  and exactly-two-non-empty-fields per TSV row, all validated before any index or output
  directory is created.
- **#15** — accepted: the single-pipeline decision is recorded as settled under a new
  "Settled decisions" heading, with the legacy-compatibility consequence stated.

Nothing rejected this round beyond the two causal claims already withdrawn above.

## Round 3 — Codex

The corrected F3/F5 analysis is sound, and preserving raw `TopHits` iteration is the right compatibility decision. Findings #4–#7, #11, #14, and #15 are now adequately addressed. The design review still finds material failures in several accepted mechanisms.

1. **The SSI mismatch guard must not be a Python `assert`.** The plan repeatedly says `assert fetched.name == name`; `python -O` removes that check entirely, restoring silent wrong-sequence retrieval ([PLAN.md:218](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:218)).

   Fix: Use an explicit `if fetched.name != name: raise IndexMismatchError(...)` and test under `python -O`.

2. **The duplicate-SSI cleanup described cannot work through the public API.** On `eslEDUP`, `SSIWriter.close()` raises before calling `esl_newssi_Close`, leaving the native writer and external-sort files open; unlinking the main temporary index in `finally` does not close them and does not prevent the known second `__dealloc__` error ([easel.pyx:10026](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer/easel.pyx:10026)). F2 also still says “detect duplicates during the scan,” contradicting §3.

   Fix: Use a disk-backed uniqueness precheck, an external indexer, or a pyhmmer API with abort/close-on-error; update F2 and test that native external-sort scratch files are removed.

3. **The proposed heartbeat still cannot heartbeat during a long query.** Pyhmmer invokes `callback` only after `self.query(query)` completes ([base worker](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer/hmmer/_base.py:320)); if one database scan takes hours, no callback occurs for hours.

   Fix: Run an independent timer thread around `hmmsearch`; let callbacks update a thread-safe completed-query counter only.

4. **Startup deletion is destructive and namespace-unsafe.** Deleting `<chunk>_*` before successful completion destroys the prior good run if the replacement fails, races another invocation, and a broad glob for chunk `foo` can delete `foo_bar_1` belonging to chunk `foo_bar` ([PLAN.md:335](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:335)).

   Fix: Stage new outputs, lock the chunk, promote only after success, and identify old artifacts with an escaped exact numeric-suffix match or manifest—not `<chunk>_*`.

5. **The 16-family floor absolutely has a cost.** At `cpus=1`, it retains 16 potentially million-record lists instead of two while gaining no parallelism; that directly contradicts the log’s claim that it “costs nothing.” The plan’s own purpose for `batch_size` is limiting live records ([PLAN.md:302](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:302)).

   Fix: Default to `2 * cpus` or a measured memory-aware cap, and benchmark peak RSS before choosing a larger floor.

6. **“Blocking buys nothing” remains false.** Query-parallel streaming rewinds and reparses the entire uncompressed database once per query; an outer target-block loop can amortize physical reads across HMMs. It remains scientifically unsuitable because raw-`TopHits` equivalence is unproven, but it has a real I/O advantage ignored by F5 and the decision table.

   Fix: Reject blocking on semantic grounds, not “nothing to buy,” and benchmark the selected design with a database larger than page cache on representative shared storage.

7. **The billion-sequence resource claim is still unverified.** Easel holds up to 2,048 MB of keys before switching to external sort, then needs substantial scratch and final-index space ([esl_ssi.h:21](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer.libs/include/libeasel/esl_ssi.h:21)); the 50,000-record test says little about that transition.

   Fix: Document required RAM/TMPDIR/disk space and test external-sort index creation before claiming billion-record scalability.

8. **The central raw-iteration decision lacks a regression test that can detect its accidental removal.** F3b is explicitly latent—the extra hit is filtered later—so the end-to-end fixture could pass even if implementation incorrectly switched to `.reported`; the verification list tests ordering but not raw-versus-reported extraction.

   Fix: Add a fake or pinned `TopHits` fixture containing an unreported domain that survives filtering, and assert extraction includes it without accessing `.reported`.

9. **`mkstemp` needs more precise handling.** It returns an open descriptor that must be closed before `SSIWriter` reopens the path, and reopening by pathname means the absolute “cannot be symlink-swapped” claim is too strong in an untrusted shared directory.

   Fix: Close the descriptor explicitly and build inside a private `0700` temporary directory in the destination filesystem, or narrow the security claim to trusted directories.

10. **The repository documentation still publishes the withdrawn explanation.** [README.md:103](/home/vangelis/Desktop/Projects/mgnifam/README.md:103) still claims per-slice `domZ`, extra reported domains, and CPU-dependent final results—the exact claims now disproven or narrowed by F3/F3b.

   Fix: Update README alongside PLAN.md to describe retained globally unreported hits and distinguish recruitment-level dependence from the fixture’s identical final artifacts.

VERDICT: REVISE
### Claude's response — Round 3

Ten findings, all substantive. **Accepted nine, narrowed one.** Four were verified by probe
before acting.

- **#1 `assert` under `-O`.** Verified: `python -O -c "assert 'WRONG'=='want'"` passes
  silently. Changed to `if fetched.name != name: raise IndexMismatchError(...)`, with a
  test that runs under `-O`. A real hole — the guard would have evaporated in any
  optimised deployment.
- **#2 SSI cleanup through the public API.** Conceded. `close()` raises before
  `esl_newssi_Close`, so unlinking one path cannot clean up the native writer or its
  external-sort scratch. Replaced with: build inside a private `mkdtemp` directory on the
  destination filesystem, `os.replace` the finished index out, `rmtree` the directory on
  any error. Verified the duplicate-key failure leaves that directory clean. F2's
  contradictory "detect duplicates during the scan" wording is removed — `SSIWriter`
  detects them itself, and a Python `set` of a billion keys was never affordable (its own
  #8 from round 2).
- **#3 heartbeat.** Conceded — `callback` fires only after a query *completes*, so a
  multi-hour scan is silent. Now a daemon timer thread around each `hmmsearch`; the
  callback only bumps a lock-guarded counter.
- **#4 startup deletion.** The namespace bug is real and embarrassing: a `<chunk>_*` glob
  for chunk `foo` matches `foo_bar_1`. Now an exact `re.fullmatch(rf"{re.escape(chunk)}_\d+\..*")`.
  Rejected the staging/locking/promotion machinery: legacy already truncates its aggregate
  files at start-up, so the run was never atomic, and under Nextflow every task gets a
  fresh work directory.
- **#5 the 16-family floor.** Conceded, and my "costs nothing" was wrong: at `cpus=1` it
  pins 16 million-record lists for zero parallelism, defeating the flag. Default is now
  `2 * cpus`.
- **#6 "blocking buys nothing" is false.** Conceded — this one caught me contradicting my
  own earlier reasoning. Query-parallel streaming re-reads and re-parses the whole database
  **once per query**; an outer target-block loop would amortise those physical reads across
  a round's HMMs. Blocking is now rejected on **semantic** grounds only (per-block pipelines
  change the *stored* hit set, which is what raw iteration consumes). The I/O win is
  recovered honestly instead: **`--prefetch_targets`**, an opt-in flag. Measured (F11):
  `hmmsearch(hmms, block, parallel="queries")` is byte-identical to
  `hmmsearch(hmms, SequenceFile)` across the full hit/score/envelope signature, so the flag
  is a pure memory-vs-time dial with zero effect on science. The streaming tax is now stated
  plainly in `PLAN.md` and the README rather than papered over.
- **#7 billion-record SSI is unverified.** Accepted as documentation + risk: Easel buffers
  up to `eslSSI_MAXRAM` (2048 MB) of keys before spilling to an external sort. `TMPDIR` and
  index-size budgeting are now in the README; a 50 000-record probe proves nothing about
  that transition and the plan now says so.
- **#8 no regression test for raw-vs-reported.** The sharpest find of the round, precisely
  because F3b makes the end-to-end test blind to it. Found the pinned example (F12): query
  `4497037939_1_144` stores 55 hits, reports 54, and the extra one is sequence
  **`6320430079`**. The test now asserts that sequence is recruited, and fails the instant
  anyone "fixes" extraction to `.reported`.
- **#9 `mkstemp` handling.** Subsumed by #2's private-directory approach: no descriptor is
  handed to `SSIWriter` at all. The security claim is narrowed as requested — safe against
  cooperating writers and PID collisions across hosts, *not* a defence against a hostile
  actor with write access to the output directory.
- **#10 README publishes the withdrawn claim.** Correct, and the exact failure mode of
  revising a plan without revising its documentation. README now describes retained
  globally-unreported hits, separates recruitment-level `cpus` dependence from the fixture's
  identical final artifacts, and documents the preserved `--recruit_evalue_cutoff` bug with
  the `6320430079` example.

Two new decisions this round: `--prefetch_targets` (opt-in, results-identical) and
`--batch_size` defaulting to `2 * cpus`.

## Round 4 — Codex

The core reasoning is now sound: the exit search is redundant, round-major scheduling can preserve numbering, `record_offset` alone is sufficient for whole-record FASTA fetches, and I accept F11 under pinned PyHMMER when `parallel="queries"` is explicit.

Material problems remain:

1. **The pseudocode still reintroduces CPU-dependent prefetch results.** The central call omits `parallel="queries"` and only shows `SequenceFile` targets ([PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:292)). With a prefetched block and a tail wave containing fewer queries than CPUs, PyHMMER automatically selects target parallelism ([dispatcher](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer/hmmer/_hmmsearch.py:410)). The proposed end-to-end prefetch test may still pass because F3b already proves final artifacts can hide the extra raw hit.

   Fix: Route both target types through one helper that always passes `parallel="queries"`, and compare exact raw signatures for one HMM at `cpus=4` in both modes, including the `55/54` counts and `6320430079`.

2. **The converge-then-discard test bypasses the state machine.** Constructing synthetic `converged & discarded` families tests only `emit()` arithmetic ([PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:558)); it will pass even if `finish()` clears the convergence marker when changing the state to discarded.

   Fix: Drive a real `Family` through `advance()` convergence, `finish()` discard, and `emit()`, asserting that an `ever_converged` flag survives.

3. **`hmmalign` can oversubscribe the task allocation.** The plan explicitly caps FAMSA but never says to pass `cpus` to PyHMMER’s `hmmalign`; its default is `cpus=0`, which expands to all detected physical CPUs ([implementation](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer/hmmer/_hmmalign.py:27), [CPU selection](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer/hmmer/_hmmalign.py:76)).

   Fix: Pass `cpus=cpus` to every seed and final `hmmalign` call and test the forwarded value.

4. **The heartbeat lifecycle is underspecified around a lazy iterator.** `hmmsearch()` returns a generator; the scan happens while `searched` is consumed, not when it is constructed ([dispatcher](/home/vangelis/Desktop/Projects/mgnifam/.venv/lib/python3.13/site-packages/pyhmmer/hmmer/_base.py:491)). A timer “around the call” can therefore stop before any search, while an unjoined daemon can log after cleanup.

   Fix: Enclose complete iterator exhaustion in a heartbeat context that sets a stop event and joins the timer in `finally`.

5. **The Stockholm verification contradicts the preserved output semantics.** The plan correctly renames only `seed_msa` ([PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:356)), but verification requires `#=GF ID` in both `seed_msa_sto` and `full_msa_sto` ([PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:588)). A PyHMMER 0.12.1 probe confirms `hmmalign` returns the full MSA with `name=None`, matching the legacy path.

   Fix: Require `#=GF ID` only in seed MSAs; merely parse and validate full MSAs unless adding their ID is an explicit new behavior.

6. **README still overclaims legacy compatibility.** “Output … matches the old script at `cpus 1`” and “faithful `cpus 1` replacement” ([README.md](/home/vangelis/Desktop/Projects/mgnifam/README.md:121)) contradict the plan’s expected clipping, tiny-MSA, dependency, gzip, and metadata differences ([PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:593)). It also says both in-memory copies are gone although prefetch deliberately restores one.

   Fix: Say only that raw recruitment/search semantics target legacy `cpus=1`, then enumerate the documented output-changing behaviors and qualify the memory statement as the default mode.

7. **The determinism evidence is version-scoped, but the contract is not.** Measurements use locked versions while package metadata permits every future PyHMMER/PyFAMSA/PyTrimAl release ([pyproject.toml](/home/vangelis/Desktop/Projects/mgnifam/pyproject.toml:8)); those libraries determine hit retention, alignment, and serialized bytes.

   Fix: Scope byte reproducibility to the committed `uv.lock`, test with `uv sync --frozen`, and avoid promising it for unconstrained wheel installations.

VERDICT: REVISE
### Claude's response — Round 4

Seven findings, **all seven accepted**. Two were live bugs in the plan, not stylistic notes.
Three were verified by probe first.

- **#1 the pseudocode reintroduces `cpus`-dependent prefetch results.** Correct, and the
  worst remaining defect. Verified in `_hmmsearch.py`: `if _few_queries and isinstance(targets,
  DigitalSequenceBlock): parallel = "targets"`. A **tail wave** with fewer families than
  `cpus`, under `--prefetch_targets`, would silently flip to target parallelism and restore
  the F3 raw-hit divergence — and per F3b the end-to-end assertions would not notice. Both
  target types now go through a single `search()` helper that always passes
  `parallel="queries"` explicitly, and the prefetch-equivalence test compares raw signatures
  at `cpus=4` including the `55/54` counts and `6320430079`.
- **#2 the converge-then-discard test bypasses the state machine.** Correct — testing `emit()`
  arithmetic against synthetic families would pass even if `finish()` cleared the convergence
  marker. `converged` becomes a sticky `ever_converged` flag, never derived from the state
  enum, and the test drives a real `Family` through `advance()` → converge, `finish()` →
  discard, `emit()`.
- **#3 `hmmalign` oversubscribes.** Verified: `cpus = cpus if cpus > 0 else
  psutil.cpu_count(logical=False)`. Legacy never passed it because `hmmalign` only gained the
  argument in pyhmmer 0.11.4. Every call now forwards `cpus=cpus`, asserted by monkeypatch.
- **#4 heartbeat lifecycle around a lazy iterator.** Correct: `hmmsearch()` returns a
  generator and the scan happens during consumption. A timer wrapped around the *call* could
  stop before any work; an unjoined daemon could log after teardown. `search()` is now a
  context manager that starts the timer, yields the iterator, and in `finally` sets a stop
  `Event` and joins the thread, enclosing full exhaustion.
- **#5 the Stockholm verification contradicts the preserved semantics.** Verified: `hmmalign`
  returns an MSA with `name=None`, which writes no `#=GF ID`; a named seed MSA writes
  `#=GF ID chunk_1`. Requiring the ID in `full_msa_sto` would have forced a new behaviour I
  had explicitly promised not to add. Now: seed MSAs must carry it, full MSAs must merely parse.
- **#6 README overclaims legacy compatibility.** Correct on both counts. It now says the
  *search and recruitment semantics* target legacy `cpus=1`, then enumerates the deliberate
  output-changing behaviours (envelope clipping, the `<= 2` discard, gzip framing, omitted
  `DATE`/`COM`, dependency upgrades), and qualifies the "both in-memory copies are gone"
  claim as true of the default mode only — `--prefetch_targets` restores one on purpose.
- **#7 the determinism contract outran its evidence.** Correct: `>=` constraints admitted
  every future pyhmmer/pyfamsa/pytrimal, and those libraries decide hit retention, alignment
  and serialised bytes. Added compatible-release upper bounds
  (`pyhmmer>=0.12.1,<0.13`, `pyfamsa>=0.7.0,<0.8`, `pytrimal>=0.8.5,<0.9`, `numpy>=2.5.1,<3`),
  re-locked, scoped the contract to the committed `uv.lock`, and moved verification to
  `uv sync --frozen`.

Nothing rejected this round.

## Round 5 — Codex

The plan is sound enough to implement. The search semantics, family numbering, exit-hit reuse, SSI strategy, and determinism scope now hold.

Non-blocking corrections:

- [PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:588) still describes synthetic converge/discard families, contradicting the real `advance → finish → emit` test specified earlier.
  Fix: Rewrite that verification bullet to describe the real lifecycle test.

- [PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:606) still lists only end-to-end prefetch equivalence, omitting the promised raw-signature comparison.
  Fix: State the `cpus=4`, one-query raw-signature assertions, including `55/54` and `6320430079`.

- The dependency table in [PLAN.md](/home/vangelis/Desktop/Projects/mgnifam/PLAN.md:216) omits the upper bounds present in [pyproject.toml](/home/vangelis/Desktop/Projects/mgnifam/pyproject.toml:12).
  Fix: Synchronize the plan table with `pyproject.toml`.

- The `search()` implementation should explicitly close the lazy PyHMMER iterator on exceptional or partial consumption; stopping the timer alone does not clean up dispatcher workers.
  Fix: Use `contextlib.closing()` around the iterator and `stop_event.wait(60)` for promptly joinable heartbeats.

- README retains minor overstatements: “saturates every core” and “faithful `cpus=1` replacement” despite its otherwise-correct qualifications.
  Fix: Say “uses up to `cpus` workers when enough families remain” and “preserves legacy-at-`cpus=1` raw recruitment semantics.”

- For CI lock consistency, `uv sync --locked` or `uv lock --check` plus `uv sync --frozen` is stronger than `--frozen` alone, which does not verify that the lock matches project metadata.

VERDICT: APPROVED
### Claude's response — Round 5 (converged)

`VERDICT: APPROVED` after 5 rounds. All six non-blocking corrections applied anyway:

1. The converge-then-discard verification bullet still described synthetic families,
   contradicting the real `advance → finish → emit` lifecycle test agreed in round 4. Rewritten.
2. The `--prefetch_targets` verification now specifies the unit-level raw-signature comparison
   at `cpus=4` (the `55`/`54` split and `6320430079`), not just end-to-end equality — which F3b
   proves is blind to it.
3. The plan's dependency table now carries the same upper bounds as `pyproject.toml`.
4. **The best of the six**: `search()` wraps the lazy iterator in `contextlib.closing()`, so an
   exception or partial consumption tears down the dispatcher's worker threads — stopping the
   timer alone leaks them. The heartbeat now blocks on `stop_event.wait(60)` instead of
   sleeping, so it joins promptly.
5. README: "saturates every core" → "uses up to `--cpus` workers when enough families remain";
   "faithful `cpus=1` replacement" → "preserves legacy-at-`cpus=1` raw recruitment semantics".
6. Verification runs `uv lock --check && uv sync --frozen`; `--frozen` alone does not verify the
   lock against project metadata.

**Act 2 complete.** 5 rounds, 45 findings raised, 41 accepted, 2 rejected with measurements,
2 causal claims of mine withdrawn under evidence. Proceeding to Act 3 (`/codex-build`).
