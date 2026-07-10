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
