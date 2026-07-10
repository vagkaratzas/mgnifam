# Plan: Port & harden `generate_families.py` into the standalone `mgnifam` repo
_Locked via grill — by Claude + vagkaratzas. Revised after Codex review round 1._

## Goal

Port the core family-generation algorithm of the `mgnifams` Nextflow pipeline
(`reference/legacy_generate_families.py`) into this standalone repo as a
`uv`-managed, `pre-commit`-clean, `pytest`-covered package. The port must
(a) stop holding the target FASTA in memory, so it scales to the
billion-sequence database; (b) remove the quadratic hot spots that made the last
full run take ~8 months across all chunks; (c) produce byte-identical scientific
output across repeated runs *and across `--cpus` values*, without spending extra
compute to get there; and (d) land four requested behaviour changes — a new
envelope-end clipping step, a small-MSA discard guard, file-handle hygiene, and
gzipped bulk outputs.

The `generate-families` console script keeps every legacy CLI flag name and
output directory name.

---

## Empirical findings (measured, not assumed)

All of these were verified against `tests/fixtures/mgnifams_input_small.fa`
(50 000 sequences) with the exact pinned versions, before this plan was frozen.
Probe scripts are disposable; the assertions they establish become tests.

**F1 — SSI random access works and is fast.** Building an SSI by byte-scanning
for `>` and writing `SSIWriter.add_key(name, fd, record_offset)` produced a
2.35 MB index in 0.05 s; 50 000 fetches through
`SequenceFile(fa, digital=False, index=SSIReader(ssi)).indexed` took 0.40 s with
**zero mismatches** against a naive parse. `add_file` takes an *integer* format
code — `SequenceFile._FORMATS["fasta"]`. Missing key → `KeyError`.

**F2 — Duplicate FASTA names are fatal to SSI.** A duplicated primary key raises
`ValueError: Index contains duplicate keys.` at `SSIWriter.close()`, plus a
secondary `EaselError ... eslEINVAL` from `__dealloc__`. Legacy's dict silently
kept the last. → Detect duplicates during the scan and fail fast.

**F3 — `parallel="targets"` changes the science.** pyhmmer auto-selects
`parallel="targets"` whenever `n_queries < cpus` — which is *every* legacy call,
since legacy passes one HMM at a time. Each worker runs its own `Pipeline` over
a chunk and therefore computes a **per-chunk `domZ`**, so domains clear the
inclusion threshold too easily and the merged `TopHits` reports extra hits:

| run | hits per query | `domZ` |
|---|---|---|
| prefetch block, `cpus=1` | 26 / 19 / 55 | 26 / 19 / 54 |
| prefetch block, `cpus=4` (auto → `targets`) | **27 / 19 / 56** | 26 / 19 / 54 |
| prefetch block, `cpus=4`, `parallel="queries"` | 26 / 19 / 55 | 26 / 19 / 54 |
| **stream `SequenceFile`, `cpus=1/3/8`** | **26 / 19 / 55** | 26 / 19 / 54 |

**The legacy pipeline's output therefore depends on `task.cpus`.** Streaming
forces `parallel="queries"`, and streaming at `cpus=1/3/8` is byte-identical to
`prefetch, cpus=1` — the single-pipeline, `cpus`-independent answer. Fixing the
memory also fixes a live reproducibility bug. Outputs will *not* match legacy
production runs (which used `cpus > 1`); they match legacy `cpus=1`.

**F4 — Hit order is already deterministic; no sort needed.** Streamed `TopHits`
iteration order was identical across `cpus = 1, 2, 4, 8`. Each query is scanned
by exactly one `Pipeline` in file order, then ranked by HMMER's own sortkey.
Codex round 1 was right: the planned tie-break `sorted()` was both unnecessary
and O(N log N) on up to ~10⁶ records. **Removed.**

**F5 — Blocked targets + `TopHits.merge` is NOT equivalent** (Codex's proposed
"simpler alternative", rejected). Reading the DB in 5 000-sequence blocks,
searching each with `Z=domZ=50000`, and merging lost hits: 26 vs 27 and 54 vs 56.
Same per-chunk-`domZ` defect as F3. It is also *slower* here (0.32 s vs 0.17 s).

**F6 — HMM files are not byte-deterministic by default.** `Builder` stamps
`hmm.creation_time` (`datetime.now()`) and `hmm.command_line` (`sys.argv`), both
serialised into the HMM. Setting them to fixed values before writing makes the
bytes stable across writes. Confirmed.

**F7 — FAMSA is thread-count independent.** Alignments were identical across
`threads = 1, 2, 4, 8` and repeatable at `threads=4`. `Aligner(threads=cpus)` is
safe. (Legacy called `Aligner()` = all cores, ignoring `--cpus`; this is a
deliberate resource-behaviour change.)

**F8 — RF line shape confirms the `clip_env_ends` contract.** A real `hmmalign`
result gave `............xxx…xxx.xxxx.x.xxx…xxx....................` — 12 leading
dots, 33 trailing dots, 3 interior dots. Trim the leading/trailing runs, keep the
interior.

---

## Background

### The upstream libraries fixed the bug that forced the FASTA into memory

`build_sequence_dict()` carries the comment *"Safely build a name → sequence
mapping without crashing like `.indexed`"*. Resolved:

- **pyhmmer 0.12.0** added `SequenceFile.indexed` (SSI-backed `Mapping[str, Sequence]`);
  **0.12.1** fixed `SequenceFile.__init__` crashing when opening SSI indices.
- **pyhmmer 0.12.0 is breaking**: `Sequence.name`, `MSA.reference`, `MSA.names`,
  `HMM.name`, SSI keys are now `str`, not `bytes`. Every `.encode()`/`.decode()`
  in the legacy script must go. (`pyfamsa.Sequence(id=...)` still takes `bytes` —
  do not "fix" that.)
- `hmmsearch` accepts a digital-mode `SequenceFile` and streams targets from disk.

### Three real performance bugs in the legacy script

1. **`run_initial_msa` is O(database × members), per family** (`reference/legacy_generate_families.py:206`):
   ```python
   seq_dict = {seq.name: seq for seq in seqs if seq.name in map(str.encode, members)}
   ```
   `map(...)` is a fresh one-shot iterator re-created for *every* candidate
   sequence, so this linearly scans `members` for each of the billions of database
   sequences — and re-runs for every cluster. Almost certainly the bulk of the 8 months.

2. **`get_next_family` is O(N²)** over the cluster table: boolean-masks the whole
   DataFrame and calls `.drop()` once per family.

3. **The exit-branch `hmmsearch` is redundant** (Codex independently confirmed the trace):
   - convergence at round *k*: the exit branch searches with `hmm` = `hmm_k`, the
     model round *k* just used;
   - `family_iteration > 3`: the main branch is skipped, so `hmm` is still `hmm_3`.

   The only difference is `exit_flag`, which merely relaxes the
   `env_length >= recruit_hit_length_percentage * qlen` filter. Cache the
   `TopHits` records and re-filter → one fewer full-database pass per family,
   provably identical results.

### The streaming/parallelism trap

`hmmsearch` raises `RuntimeError: cannot use ``targets`` parallel mode with a
sequence file`. So streaming with **one** HMM per call occupies one worker
thread — a ~`cpus`× regression. Batching a round's HMMs into a single
`hmmsearch` call restores `parallel="queries"` scaling. `hmmsearch` yields
`TopHits` *"in the same order the queries were passed in"* (verified). Batching
is therefore both necessary and order-safe.

---

## Approach

### 1. Repo layout

```
pyproject.toml            # uv-managed, hatchling backend, console script
uv.lock
.python-version           # 3.13
.pre-commit-config.yaml
README.md
src/mgnifam/__init__.py
src/mgnifam/generate_families.py     # single implementation module
tests/conftest.py                    # decompresses .fa.gz fixtures to tmp_path
tests/fixtures/{clustering.tsv,cluster_long.tsv}
tests/fixtures/{mgnifams_input_small.fa.gz,mgnifams_extra.fa.gz}
tests/test_generate_families.py
reference/legacy_generate_families.py   # the port's source of truth
```

### 2. Dependencies

`requires-python = ">=3.12"` (numpy 2.5.1 floor); `.python-version` = `3.13`.

| package | version |
|---|---|
| `pyhmmer` | `>=0.12.1` |
| `pyfamsa` | `>=0.7.0` |
| `pytrimal` | `>=0.8.5` |
| `numpy` | `>=2.5.1` |

**Dropped**: `pandas` and `biopython` (the latter imported nowhere).
Dev group: `pytest`, `pre-commit`. `uv sync` verified clean on 3.13.

### 3. SSI-backed random access — no FASTA in memory

`build_ssi_index(fasta, ssi_path)` scans the FASTA once for `>` line offsets and
writes them via `SSIWriter`: `add_file(basename, SequenceFile._FORMATS["fasta"])`
then `add_key(name, fd, record_offset)` (F1).

- **Duplicate names** are collected during the scan; if any exist, raise with the
  first few offending names before touching `SSIWriter` (F2).
- **Atomic build**: write to `<ssi>.tmp.<pid>` then `os.replace()`. Rebuild when
  the existing index is older than the FASTA. No sidecar-fingerprint system —
  that is a filesystem-integrity problem, not this script's job.
- **`--fasta_index`** (new, optional) points at a pre-built index. Documented as
  the production path: build the index **once**, upstream, and share it across
  chunk tasks. Auto-build is the small-data convenience.
- **Gzipped FASTA is rejected** with an actionable error — Easel cannot position
  within a gzip stream. Test fixtures ship as `.fa.gz` and `conftest.py`
  decompresses them into `tmp_path`.

Two handles:
- `SequenceFile(fasta, digital=True, alphabet=ALPHABET)` → streaming `hmmsearch` targets.
- `SequenceFile(fasta, digital=False, index=SSIReader(ssi))` → `.indexed[name]`
  replaces `seq_dict` everywhere. Text mode yields `str` sequences directly, as
  `renumber_sto_msa`'s `original_seq.find(seq)` wants.

`read_pyhmmer_seqs` and `build_sequence_dict` are **deleted**.

**Missing-member policy**: `run_initial_msa` silently skips members absent from
the FASTA (legacy behaviour — `get_fasta_sequences` filtered `None`). Everywhere
else a missing name is a hard error, as it is in legacy (`[0]` → `IndexError`).

### 4. Round-major ("wave") scheduling, with all side effects deferred

Clusters are grouped into an **insertion-ordered dict** (`setdefault(rep,
[]).append(member)`), not `itertools.groupby` — legacy gathers *all* rows for a
representative even when they are non-contiguous, and `groupby` would split
`A, B, A` into three families instead of two.

```python
for batch in batched(clusters.items(), batch_size):   # cluster-file order
    active = [Family(rep, members) for rep, members in batch]
    for f in active: f.seed_msa = initial_msa(f)      # may discard (empty/1 seq)
    for round_ in (1, 2, 3):
        running = [f for f in active if f.state is RUNNING]
        if not running: break
        hmms = [hmmbuild(f.seed_msa, hand=False) for f in running]
        with SequenceFile(fasta, digital=True, alphabet=A) as targets:
            for f, hits in zip(running, hmmsearch(hmms, targets, cpus=cpus,
                                                  E=evalue, seed=42)):
                f.hmm, f.qlen = hmms[i], hits.query.M
                f.records = [rec for hit in hits for rec in hit.domains]  # native order
                f.advance(round_)          # pure: mutates f, writes nothing
    for f in active:
        if f.state in (RUNNING, CONVERGED): f.finish()   # re-filters f.records, NO new search
    for f in active:                       # cluster order → stable ids, ordered side effects
        f.emit(writers)
```

Every state transition is **side-effect free**. Convergence records, discard
records, metadata rows and family ids are buffered on the `Family` dataclass and
flushed in `emit()`, in cluster order — otherwise a family converging in round 1
would write its `converged_families` line before an earlier family that succeeds
in round 3 (Codex #4). Family ids remain a 1-based rank among *successful*
families in cluster-file order.

`--batch_size` (default `max(2 * cpus, 16)`) caps how many families' hit-record
lists are live at once. It does **not** change I/O: `hmmsearch` spawns exactly
`cpus` workers regardless of query count, and each query is read by one worker,
so total file reads == total queries either way.

`Builder.build_msa` requires a named MSA; the HMM name does not influence the
model, so the seed MSA is named after its cluster representative and
`final_hmm.name` is set to `f"{chunk}_{id}"` immediately before writing.

### 5. Exit branch reuses cached hits

`filter_hits(records, qlen, exit_flag, recruit_hit_length_percentage)` runs twice
over the same records — `exit_flag=False` during the round, `exit_flag=True` in
the exit branch. **Native `TopHits` order is preserved throughout; nothing is
re-sorted** (F4).

### 6. Determinism — by construction, at no extra compute cost

- **`parallel="queries"` is forced** (streaming mandates it), so results no longer
  depend on `--cpus` (F3). Test: `cpus=1` vs `cpus=4` → byte-identical outputs.
- **No hit sort.** Native ranking is already stable (F4).
- Cluster iteration order comes from the file, never from a `set`.
- `seed=42` on `hmmsearch` pipeline options and on `Builder`.
- **`hmm.creation_time = datetime(1970,1,1)` and `hmm.command_line =
  "generate-families"` before every write** (F6).
- **`gzip.GzipFile(fileobj=fh, mode="wb", filename="", mtime=0)`** — `filename=""`
  is required, otherwise `GzipFile` copies `fileobj.name` into the header (Codex #14).
- `total_checked_sequences` is a `set` used only for membership.
- `pyfamsa.Aligner(threads=cpus)` (F7).
- **The byte-determinism contract covers scientific artifacts only** —
  `logs/` carries timestamps and durations by design and is excluded (Codex #13).

### 7. The four requested behaviour changes

**7.1 — `clip_env_ends(msa)` (new).** Applied to the result of *every* seed-path
`hmmalign` (the `# main strategy continue` section); **not** applied to the final
full-MSA `hmmalign` in the exit branch. Reads the `#=GC RF` line
(`msa.reference`), finds the first and last non-`.` column, returns
`msa.select(columns=range(first, last + 1))`. Interior `.` columns (inserts
between match states) are untouched — only the leading and trailing runs, i.e.
the N- and C-terminal envelope overhangs (F8).

`run_hmmalign` is reduced to returning the MSA; a `msa_stats(msa) -> (num_seqs,
non_gap_rep_length)` helper is called *after* clipping on the seed path and
*directly* on the unclipped full MSA in the exit branch.

The existing gap-occupancy `clip_ends()` is a **different** function, stays after
`run_pytrimal_reps`, unchanged.

**7.2 — Discard tiny MSAs.** After `run_pytrimal_reps`: `<= 2` sequences →
discard, `discard_reason = "too few sequences after redundancy filtering"`,
`discard_value = <count>`. Additionally, a `< 2` (empty/single) guard immediately
before the *initial* and *hand* builds converts the `eslEMEM (status 5)` crash
into a clean discard. The initial guard is deliberately `< 2`, not `<= 2`, so
two-member clusters keep their legacy semantics (Codex #11, accepted with
modification).

**7.3 — File-handle and dead-state cleanup.** All output handles open once under a
single `contextlib.ExitStack` in `main()` and are passed to `emit()`. Legacy's
fixed `tmp/seed_msa.sto` / `tmp/full_msa.sto` and cwd `execution.log` (which
collide when chunks share a working directory) are replaced by a per-invocation
`tempfile.TemporaryDirectory()` and a logger writing straight to
`logs/<chunk>.txt` — no `shutil.move` at exit (Codex #16). `discard_value` is a
single `Family` field, not four re-assignments.

**7.4 — Gzip the bulk outputs.** `seed_msa_sto/*.sto.gz`, `full_msa_sto/*.sto.gz`,
`hmm/*.hmm.gz`, `family_reps/<chunk>.fasta.gz`. The small per-chunk text/CSV/TSV
outputs (`rf/`, `refined_families/`, `discarded_clusters/`,
`successful_clusters/`, `converged_families/`, `family_metadata/`, `logs/`) stay
plain — they are what a human greps.

### 8. Input validation & observability

- `--chunk_num` must match `[A-Za-z0-9._-]+`; it is interpolated into output
  paths, so `../` would escape the output directories (Codex #17).
- `--fasta_file` must not be gzipped (§3).
- Progress logging emits one structured line per (batch, round): families active,
  queries dispatched, hits recruited, elapsed. Enough to see a multi-hour database
  pass making progress. **Checkpoint/resume is out of scope** — Nextflow already
  resumes at chunk granularity (Codex #21, partially accepted).

---

## Key decisions & tradeoffs

| Decision | Alternative rejected | Why |
|---|---|---|
| Streaming `SequenceFile` + batched round-major waves | Streaming with the per-family loop | pyhmmer refuses `parallel="targets"` on a `SequenceFile`; 1 query ⇒ 1 thread. |
| Streaming targets | Keep the prefetched `DigitalSequenceBlock` | Prefetch is the O(database) memory the user asked to eliminate — and its `parallel="targets"` default makes results `cpus`-dependent (F3). |
| Streaming targets | Blocked `read_block` + `TopHits.merge` + fixed `Z` | **Measured**: loses hits (26 vs 27, 54 vs 56) via per-chunk `domZ`, and was slower (F5). |
| Cache exit-branch hits, re-filter | Re-run the exit `hmmsearch` | Same HMM as the preceding round (Codex agreed). Saves one full-DB pass per family. |
| **No** tie-break sort | Sort records by `(evalue, -score, name, env_from)` | **Measured**: native order already stable across `cpus` (F4); the sort was O(N log N) on ~10⁶ records for nothing. |
| Insertion-ordered dict for clusters | `itertools.groupby` | `groupby` only groups *contiguous* rows; legacy gathers all rows per representative. |
| Discard `<= 2` after trimal, `< 2` before initial/hand builds | `<= 2` everywhere | `<= 2` on the initial cluster would newly discard 2-member clusters — not requested. |
| Fail fast on duplicate FASTA names | Keep-last like legacy's dict | **Measured**: `SSIWriter` raises anyway (F2); a protein DB with duplicate accessions is a data bug. |
| Atomic rename + mtime staleness check | Sidecar fingerprint file + lock | Guards the real race; the rest is filesystem integrity, out of scope. |
| One implementation module | Package split | It is one script. |

## Preserved as-is (deliberately — do NOT "fix" these)

- The gap-occupancy `clip_ends()` off-by-one: `range(start_position,
  end_position)` drops the last column that passed the threshold. Out of scope;
  mark with a comment.
- In the `family_iteration > 3` exit path the written `final_hmm` (hand, from
  round 3's seed MSA) is *not* the model used for the final `hmmsearch`/`hmmalign`
  (`hmm_3`, from round 2's seed MSA). Odd, faithfully reproduced.
- All legacy CLI flag names and output directory names.
- `renumber_sto_msa`'s duplicate-name skipping and `parse_protein_name`'s
  pad/trim-to-old-length behaviour.

## Risks / open questions

- **Results will not reproduce legacy production output.** Legacy ran with
  `cpus > 1` ⇒ `parallel="targets"` ⇒ inflated hit sets (F3). The new output
  matches legacy at `cpus=1`. On top of that, `clip_env_ends` changes seed MSAs,
  hence downstream HMMs, recruitment, convergence and family counts. This is a
  correction, and the user should confirm they want it.
- **`SequenceFile._FORMATS` is private API.** Working (F1), but pinned to
  pyhmmer's internals. A regression test covers it.
- **Per-family memory** is driven by a single family's hit count (~10⁶ records),
  which `--batch_size` cannot bound. Recorded, not solved.

## Out of scope

- The `mgnifams` Nextflow module (`main.nf`, `meta.yml`, nf-test snapshots).
- Changing E-value / length-percentage defaults or any scientific threshold.
- Checkpoint/resume; multi-node distribution.
- The two "preserved as-is" quirks.

---

## Verification

1. `uv sync && uv run pre-commit run --all-files` — clean.
2. `uv run pytest`:
   - **`build_ssi_index` round-trip** — every FASTA name fetchable, sequences equal
     to a naive parse; missing key → `KeyError`; duplicate names → clear error (F1, F2).
   - **`clip_env_ends`** — RF `"...xx.xxx..."` → leading/trailing dot runs removed,
     interior dot kept (F8).
   - **`filter_hits`** — `exit_flag=True` admits short envelopes that
     `exit_flag=False` rejects; record order is the input `TopHits` order.
   - **small-MSA discard** — 2-sequence alignment out of `run_pytrimal_reps` →
     discarded, `hmmbuild` never called, `eslEMEM` never raised.
   - **`--cpus` invariance** — end-to-end at `cpus=1` and `cpus=4` → byte-identical
     scientific outputs (the bug F3 exposes).
   - **batch invariance** — end-to-end at `batch_size=1` and `batch_size=64` →
     byte-identical outputs, *and* an identical `representative → family_id`
     mapping (Codex #5: contiguity alone would accept ids on the wrong reps).
   - **determinism** — two runs, different `PYTHONHASHSEED`, different output dirs
     → every scientific artifact byte-identical; `logs/` excluded (Codex #13).
     Asserts the HMM `creation_time`/`command_line` neutralisation (F6) and the
     `filename=""` gzip header (Codex #14).
   - **no in-memory FASTA** — monkeypatch `SequenceFile.read_block` to raise; the
     end-to-end run must still pass.
   - **`chunk_num` validation** — `../evil` rejected.
   - end-to-end on `clustering.tsv` + `mgnifams_input_small.fa`, and on
     `cluster_long.tsv` + `mgnifams_extra.fa` (>200 aa): exit 0, every declared
     output dir exists, `*_msa_sto/*.sto.gz` parse as Stockholm, `hmm/*.hmm.gz`
     reloads via `HMMFile`, `family_metadata/*.csv` has one row per successful family.
3. **Legacy comparison is informational, not an assertion** (Codex #19). Run the
   legacy script at `cpus=1` on `clustering.tsv` and record the deltas in
   `refined_families/`, `successful_clusters/`, `family_metadata/`. Expected
   sources of divergence, to be explained one by one: `clip_env_ends`, the `<= 2`
   discard rule, and dependency upgrades. Any *other* delta is a port bug.

---

## Progress

Branch `dev`. **One feature, one commit, no push.**

### Act 2 — Codex adversarial review
- [x] Freeze `PLAN.md` + init `PLAN-REVIEW-LOG.md`
- [x] Round 1 → `VERDICT: REVISE` (21 findings; 16 accepted, 3 accepted-with-modification, 2 rejected with measurements)
- [ ] Round 2 → re-review until `VERDICT: APPROVED` (MAX_ROUNDS=5)

### Act 3 — Build
- [ ] 1. `chore: scaffold uv project, pre-commit, README`
- [ ] 2. `feat: port generate_families to src layout on pyhmmer 0.12 str API` (drops pandas + biopython)
- [ ] 3. `perf: SSI-indexed random access, streaming hmmsearch targets`
- [ ] 4. `perf: batch families into round-major hmmsearch waves, reuse exit-branch hits`
- [ ] 5. `feat: deterministic, cpus-invariant outputs`
- [ ] 6. `feat: clip envelope ends from seed MSAs via the RF line`
- [ ] 7. `feat: discard families with <= 2 sequences after redundancy trimming`
- [ ] 8. `feat: gzip alignment, HMM and FASTA outputs`
- [ ] 9. `fix: close output handles via context managers, single discard_value`
- [ ] 10. `test: pytest suite and fixtures`

### Act 4 — Verification
- [ ] `uv run pre-commit run --all-files` clean
- [ ] `uv run pytest` green (determinism, cpus-invariance, batch-invariance, no-`read_block`)
- [ ] Legacy-vs-new delta report
