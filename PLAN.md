# Plan: Port & harden `generate_families.py` into the standalone `mgnifam` repo
_Locked via grill — by Claude + vagkaratzas. Hardened over 5 rounds of adversarial review by Codex (`gpt-5.6-sol`); `VERDICT: APPROVED`. Full transcript: `PLAN-REVIEW-LOG.md`._

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

**F2 — Duplicate FASTA names are fatal to SSI, and `SSIWriter` detects them
itself.** A duplicated primary key raises `ValueError: Index contains duplicate
keys.` from `SSIWriter.close()` (plus a non-fatal *"Exception ignored in
`__dealloc__`"* `EaselError` on stderr). Legacy's dict silently kept the last.
No Python-side key `set` is needed — and none is affordable at a billion records
(Codex #8). Since `close()` raises *before* `esl_newssi_Close`, the native writer
and any external-sort scratch files are left open, so cleanup must not rely on
unlinking one known path (Codex #2): the build happens inside a **private
`mkdtemp` directory on the destination filesystem**, which is `rmtree`d
wholesale on any error. Verified: after the duplicate error, that directory is
empty of stray artifacts.

**F3 — `parallel="targets"` makes legacy recruitment depend on `--cpus`, and the
cause is unreported hits, not `domZ`.** pyhmmer auto-selects `parallel="targets"`
whenever `n_queries < cpus` — *every* legacy call, since legacy passes one HMM at
a time. Measured identically under **pyhmmer 0.11.1 (the legacy pin)** and 0.12.1:

| cpus | `len(TopHits)` | `.reported` | `.included` | `Z` | `domZ` |
|---|---|---|---|---|---|
| 1 | 26 / 19 / **55** | 26 / 19 / 54 | 26 / 19 / 54 | 50000 | 26 / 19 / 54 |
| 2 | 26 / 19 / 55 | — | — | 50000 | — |
| 4 | **27** / 19 / **56** | 26 / 19 / 54 | 26 / 19 / 54 | 50000 | 26 / 19 / 54 |
| 8 | **27** / 19 / **56** | — | — | 50000 | — |

`Z` and `domZ` are **identical** across runs. My round-1 "per-chunk `domZ`"
explanation was wrong; Codex's hypothesis is right. Target-parallel merge
concatenates each worker's locally-stored hits and re-thresholds the *flags*
without deleting the entries, so `len(TopHits)` grows while `.reported` /
`.included` stay fixed. Legacy iterates the **raw** `TopHits` (`for hit in
top_hits` at `reference/legacy_generate_families.py:267`), not `.reported`, so it
recruits sub-threshold hits — 1 of them even at `cpus=1` (55 stored vs 54
reported) — and the count varies with `--cpus`.

Streaming forces `parallel="queries"`; streaming at `cpus = 1, 3, 8` is
byte-identical to `prefetch, cpus=1`. So the port reproduces **legacy at
`cpus=1`** exactly, including its raw-iteration semantics, and is `cpus`-invariant.

> **Reported to the user, not silently fixed:** iterating `top_hits.reported`
> would drop the sub-threshold recruits and is arguably what
> `--recruit_evalue_cutoff` is supposed to mean. That is a scientific change
> nobody asked for. This port keeps raw iteration. See "Risks".

**F3b — end-to-end, the divergence does not reach the fixture's families.**
Running the legacy script at `cpus=1` and `cpus=4` on `clustering.tsv`, every
output is identical except the HMM `DATE`/`COM` lines: the extra stored hit is
removed downstream by the `env_length >= 0.9 * qlen` filter. The hazard is real
at the recruitment step but latent on this data.

**F4 — Hit order is already deterministic; no sort needed.** Streamed `TopHits`
iteration order was identical across `cpus = 1, 2, 4, 8`. Each query is scanned
by exactly one `Pipeline` in file order, then ranked by HMMER's own sortkey.
Codex round 1 was right: the planned tie-break `sorted()` was both unnecessary
and O(N log N) on up to ~10⁶ records. **Removed.**

**F5 — Blocked targets + `TopHits.merge` is rejected on semantic grounds.** The
round-1 measurement (blocks of 5 000 with `Z=domZ=50000` losing hits vs a single
pass) compared against a baseline whose *automatic* `domZ` was 26/19/54, so it
changed a scientific parameter and was not like-for-like — Codex is right to
reject that causal claim. Blocking is dropped because a per-block `Pipeline`
changes the *stored* hit set, and stored hits are exactly what legacy's raw
`TopHits` iteration consumes (F3); equivalence is therefore unproven and the
science could move. It is **not** dropped for "buying nothing": blocking really
would amortise physical database reads across the HMMs of a round, which
query-parallel streaming cannot. That I/O win is recovered instead by F11.

**F11 — `--prefetch_targets` is a free knob.** `hmmsearch(hmms, block, cpus=4,
parallel="queries")` is **byte-identical** to `hmmsearch(hmms, SequenceFile, cpus=4)`
across the full hit/score/envelope signature. So prefetching the database into a
`DigitalSequenceBlock` is a pure memory-vs-time tradeoff with *zero* effect on
results, as long as `parallel="queries"` is forced in both. Streaming stays the
default (the O(1)-memory requirement); the flag exists for when the database fits
in RAM. `parallel="targets"` is never used in either mode.

**F12 — A pinned example of a recruited-but-unreported hit.** For query
`4497037939_1_144` on the small fixture, `len(TopHits) == 55` and
`len(.reported) == 54`; the extra entry is sequence **`6320430079`**. Legacy's raw
iteration recruits it. This is the regression fixture that detects an accidental
switch to `.reported` (Codex #8/round 3) — the end-to-end test cannot, because
F3b shows the difference is filtered downstream.

**F6 — HMM files are not byte-deterministic by default.** `Builder` stamps
`hmm.creation_time = datetime.now()` (serialised as `DATE`, via locale-dependent
`%a`/`%b`) and `hmm.command_line = sys.argv` (`COM`). Measured: setting **both to
`None` omits the `DATE` and `COM` lines entirely** and makes the bytes stable
across writes. That is cleaner than a fixed date and sidesteps the `LC_TIME`
sensitivity Codex raised.

**F7 — FAMSA is thread-count independent.** Alignments were identical across
`threads = 1, 2, 4, 8` and repeatable at `threads=4`. `Aligner(threads=cpus)` is
safe. (Legacy called `Aligner()` = all cores, ignoring `--cpus`; this is a
deliberate resource-behaviour change.)

**F13 — Legacy's seed Stockholm has no `#=GF` lines.** `renumber_sto_msa` copies
only `# STOCKHOLM`, `#=GC RF`, `//` and sequence rows; every `#=GF`/`#=GS`/`#=GR`
line is dropped. Confirmed against a real legacy run: `grep -c '#=GF'
seed_msa_sto/test_1.sto` → `0`, while `hmm/test_1.hmm` carries `NAME  test_1`.
Setting `seed_msa.name` therefore matters only because `Builder.build_msa` requires
a name and propagates it to `HMM.name`.

**F8 — RF line shape confirms the `clip_env_ends` contract.** A real `hmmalign`
result gave `............xxx…xxx.xxxx.x.xxx…xxx....................` — 12 leading
dots, 33 trailing dots, 3 interior dots. Trim the leading/trailing runs, keep the
interior.

**F9 — `hmmsearch` uses `min(cpus, n_queries)` workers**, not `cpus`
(`_base.py`: `self.cpus = 1 if hint == 0 else min(cpus, hint) if hint > 0 else cpus`).
Codex is right and my round-1 rebuttal of #18 was wrong on this point. Consequence:
`--batch_size` must stay **≥ `cpus`** or cores idle. Total database reads still
equal total queries, independent of batch size.

**F10 — A `Domain` keeps its parent `Hit`, which keeps the whole `TopHits`.**
Caching `hit.domains` objects would pin every HMMER result graph for the batch.
Records must be copied into plain immutable tuples at extraction time.

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

| package | constraint |
|---|---|
| `pyhmmer` | `>=0.12.1,<0.13` |
| `pyfamsa` | `>=0.7.0,<0.8` |
| `pytrimal` | `>=0.8.5,<0.9` |
| `numpy` | `>=2.5.1,<3` |

Upper bounds are deliberate: these libraries decide hit retention, alignment and
serialised bytes, so the reproducibility contract is scoped to the resolved set in
`uv.lock` (Codex #7/round 4).

**Dropped**: `pandas` and `biopython` (the latter imported nowhere).
Dev group: `pytest`, `pre-commit`. `uv sync` verified clean on 3.13.

### 3. SSI-backed random access — no FASTA in memory

`build_ssi_index(fasta, ssi_path)` scans the FASTA once for `>` line offsets and
writes them via `SSIWriter`: `add_file(basename, SequenceFile._FORMATS["fasta"])`
then `add_key(name, fd, record_offset)` (F1).

- **Duplicate names**: do **not** build a Python `set` of every key — on a
  billion-record database that is database-sized Python state, exactly what this
  port exists to remove (Codex #8). `SSIWriter` detects duplicates itself and
  raises `ValueError: Index contains duplicate keys.` from `close()` (F2). Catch
  it and re-raise as `DuplicateSequenceName` with actionable guidance.
- **Atomic build, in a private directory.** `close()` raises *before*
  `esl_newssi_Close`, so the native writer and any external-sort scratch files
  stay open; unlinking one known path cannot clean up after it (Codex #2). So:
  `tempfile.mkdtemp(dir=<ssi parent>)` → build the index inside that private
  `0700` directory → `os.replace()` the finished file into place → `rmtree` the
  directory. On **any** exception, `rmtree` it regardless. This also subsumes
  Codex #9/#10: the temp name is exclusive and unpredictable, and no descriptor
  is handed to `SSIWriter` (`mkstemp` would have returned an open fd that must be
  closed before `SSIWriter` reopens the path). The narrower claim: this is safe
  against concurrent *cooperating* writers and PID collisions; it is not a
  defence against a hostile actor with write access to the output directory.
  Rebuild when the existing index is older than the FASTA.
- **Index/FASTA binding**: every `.indexed[name]` fetch checks
  `if fetched.name != name: raise IndexMismatchError(...)`. **Not an `assert`** —
  `python -O` strips those, restoring the silent wrong-sequence read the check
  exists to prevent (Codex #1, verified). A stale or mismatched index then fails
  loudly on first use. This is O(1) and catches the real hazard; a fingerprint
  sidecar (resolved path, size, `mtime_ns`, content hash, completion marker) is
  rejected as filesystem-integrity engineering this script does not need.
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
Record = tuple[str, int, int, int]     # (name, tlen, env_from, env_to) — plain, immutable

for batch in batched(clusters.items(), batch_size):   # cluster-file order
    active = [Family(rep, members) for rep, members in batch]
    for f in active: f.seed_msa = initial_msa(f)      # may discard (empty/1 seq)
    for round_ in (1, 2, 3):
        running = [f for f in active if f.state is RUNNING]
        if not running: break
        hmms = [hmmbuild(f.seed_msa, hand=False) for f in running]
        with search(hmms, targets, cpus, evalue) as searched:   # see below
            for f, hmm, hits in zip(running, hmms, searched, strict=True):
                f.hmm, f.qlen = hmm, hits.query.M
                # copy out now: a Domain pins its Hit, which pins the whole TopHits (F10)
                f.records = [(hit.name, hit.length, d.env_from, d.env_to)
                             for hit in hits for d in hit.domains]   # native order
                f.advance(round_)          # pure: mutates f, writes nothing
    for f in active:
        if f.state in (RUNNING, CONVERGED): f.finish()   # re-filters f.records, NO new search
    for f in active:                       # cluster order → stable ids, ordered side effects
        f.emit(writers)
```

`zip(..., strict=True)` — plain `zip` would silently truncate if the dispatcher
ever yielded fewer `TopHits` than queries (Codex #4).

**Both target types go through one `search()` helper that always passes
`parallel="queries"` explicitly.** Leaving `parallel=None` is a live bug: a *tail
wave* with fewer families than `cpus`, against a prefetched `DigitalSequenceBlock`,
trips `_few_queries and isinstance(targets, DigitalSequenceBlock)` in
`_hmmsearch.py` and silently reverts to target parallelism — reintroducing exactly
the `cpus`-dependent raw hit set of F3, in the one mode where it is possible
(Codex #1/round 4). And per F3b the end-to-end assertions would not notice.

`search()` is also the **heartbeat context manager**. `hmmsearch()` returns a lazy
generator: the database scan happens while the result is *consumed*, not when the
call returns (`_base.py`). A timer merely wrapped around the call could therefore
stop before any work began, and an un-joined daemon could log after teardown
(Codex #4/round 4). So `search()`:

- wraps the iterator in `contextlib.closing()`, so an exception or partial
  consumption tears down the dispatcher's worker threads instead of leaking them —
  stopping the timer alone does not do that (Codex, round 5);
- starts a daemon timer that blocks on `stop_event.wait(60)` rather than sleeping,
  so it joins promptly instead of up to a minute late;
- in `finally`, sets the stop `Event` and **joins** the timer — enclosing full
  iterator exhaustion.

Every state transition is **side-effect free**. Convergence records, discard
records, metadata rows and family ids are buffered on the `Family` dataclass and
flushed in `emit()`, in cluster order — otherwise a family converging in round 1
would write its `converged_families` line before an earlier family that succeeds
in round 3 (Codex round 1 #4).

**Family-id semantics, reproduced exactly (Codex #7).** Legacy keeps a running
`iteration` counter: `+1` per family, `-1` on discard, so a successful family's id
is its 1-based rank among *successful* families in cluster-file order. But legacy
writes the `converged_families` line **at convergence, before** the membership and
rep-length checks that may still discard the family. A converged-then-discarded
family therefore emits its provisional id, and the next family *reuses* it — so
`converged_families` can hold ids that end up belonging to a different family, or
duplicates. `emit()` reproduces this precisely:

```python
provisional = success_count + 1
if f.ever_converged: converged_file.write(f"{provisional}\n")  # even if discarded below
if f.discarded:      discarded_file.write(...)                 # success_count unchanged
else:                success_count += 1; family_id = provisional; ...
```

`ever_converged` is a **sticky flag**, never derived from the state enum: `finish()`
may subsequently move the family to `DISCARDED`, and legacy had already written the
`converged_families` line by then. A test that only exercises `emit()` arithmetic
would pass even if `finish()` cleared the marker, so the test drives a real `Family`
through `advance()` → convergence, `finish()` → discard, `emit()` (Codex #2/round 4).

A fixture with converged → discarded → successful families pins this behaviour.

`--batch_size` (default `2 * cpus`) caps how many families' record lists are live
at once. It must stay **≥ `cpus`**: `hmmsearch` uses `min(cpus, n_queries)` workers
(F9), so a smaller batch idles cores. No 16-family floor — Codex is right that at
`cpus=1` it would pin 16 potentially million-record lists for zero parallelism,
contradicting the flag's whole purpose (Codex #5/round 3).

**The streaming tax, and the escape hatch.** Query-parallel streaming rewinds and
*re-parses the entire database once per query*: `n_queries` full passes per batch,
against prefetch's single parse. Codex is right that my "blocking buys nothing"
was false — an outer target-block loop really would amortise physical reads across
HMMs. Blocking is rejected on **semantic** grounds instead: per-block `Pipeline`s
produce different *stored* hit sets, and stored hits are exactly what legacy's raw
`TopHits` iteration consumes (F3), so equivalence is unproven and the science would
move.

Instead, `--prefetch_targets` (opt-in, default off) loads the database into a
`DigitalSequenceBlock` once and searches it with `parallel="queries"`. **Measured:
prefetch + `parallel="queries"` is byte-identical to streaming.** So the flag is a
pure memory-vs-time knob with *zero* effect on results — use it when the database
fits in RAM, leave it off to honour the O(1)-memory requirement. Both paths force
`parallel="queries"`; `parallel="targets"` is never used.

`Builder.build_msa` requires a named MSA. Legacy's `run_hmmbuild` sets
`seed_msa.name = f"{chunk}_{iteration}"`, and that name is serialised into the
seed Stockholm file as `#=GF ID` (Codex #6). So the seed MSA is renamed to
`f"{chunk}_{family_id}"` in `emit()`, immediately before the hand build and the
write — *not* named after the representative. `final_hmm.name` takes the same
value. The HMM name does not influence the model, so intermediate rounds may use
any placeholder name.

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
- **`hmm.creation_time = None` and `hmm.command_line = None` before every write**,
  which omits the `DATE` and `COM` lines rather than writing a fixed date through
  locale-dependent `%a`/`%b` formatting (F6, Codex #11).
- **`gzip.GzipFile(fileobj=fh, mode="wb", filename="", mtime=0)`** — `filename=""`
  is required, otherwise `GzipFile` copies `fileobj.name` into the header (Codex #14).
- **Stale-output cleanup**: at start-up, per-family outputs in `seed_msa_sto/`,
  `full_msa_sto/`, `hmm/` and `rf/` whose names match
  `re.fullmatch(rf"{re.escape(chunk)}_\d+\..*")` are removed — an **exact**
  numeric-suffix match, because a `<chunk>_*` glob for chunk `foo` would also
  delete `foo_bar_1` belonging to chunk `foo_bar` (Codex #4/round 3). This mirrors
  legacy's start-up truncation of its per-chunk aggregate files, so a failed run
  is no more destructive than it already was. Staging + chunk locking + atomic
  promotion is rejected: under Nextflow every task already runs in a fresh work
  directory, so the only exposure is a local rerun.
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

**Every `hmmalign` call passes `cpus=cpus`.** Its default is `cpus=0`, which expands
to `psutil.cpu_count(logical=False)` — every physical core on the node, oversubscribing
the Nextflow task allocation (Codex #3/round 4). Legacy never passed it because
`hmmalign` only gained a `cpus` argument in pyhmmer 0.11.4.

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

All validated **before** any index or output directory is created (Codex #14):

- `--chunk_num` matches `[A-Za-z0-9._-]+` and is non-empty; it is interpolated
  into output paths, so `../` would escape the output directories (Codex #17).
- `--fasta_file` exists and is not gzipped (§3).
- `--cpus >= 1`; `--batch_size >= 0` (`0` = auto).
- `--discard_min_starting_membership`, `--max_seq_identity`, `--max_gap_occupancy`,
  `--recruit_hit_length_percentage` in `[0, 1]`; `--recruit_evalue_cutoff > 0`;
  `--max_seed_seqs >= 1`.
- `--discard_min_rep_length <= --discard_max_rep_length`, both `>= 1`.
- Each cluster TSV row has exactly two non-empty tab-separated fields.

Observability: pyhmmer calls `callback` only *after* a query completes
(`_base.py`: the worker invokes it after `self.query(query)` returns), so on a
database large enough for one scan to take hours the callback is silent for hours
(Codex #3/round 3). So the heartbeat is a **daemon timer thread** started around
each `hmmsearch` call, logging every 60 s: batch, round, families active, queries
completed / total, elapsed. The `hmmsearch` `callback=` only bumps a
`threading.Lock`-guarded counter that the timer reads. Heartbeat output goes to
`logs/`, outside the determinism contract.
**Checkpoint/resume is out of scope** — Nextflow already resumes at chunk
granularity.

---

## Key decisions & tradeoffs

| Decision | Alternative rejected | Why |
|---|---|---|
| Streaming `SequenceFile` + batched round-major waves | Streaming with the per-family loop | pyhmmer refuses `parallel="targets"` on a `SequenceFile`; 1 query ⇒ 1 thread. |
| Streaming targets | Keep the prefetched `DigitalSequenceBlock` | Prefetch is the O(database) memory the user asked to eliminate — and its `parallel="targets"` default makes recruitment `cpus`-dependent (F3). |
| Streaming targets, `--prefetch_targets` opt-in | Blocked `read_block` + `TopHits.merge` | Rejected on **semantic** grounds: per-block pipelines change the *stored* hit set, which is what legacy's raw iteration consumes. Blocking *would* amortise physical reads (Codex is right), so the I/O win is instead offered via `--prefetch_targets`, measured byte-identical to streaming. |
| Cache exit-branch hits, re-filter | Re-run the exit `hmmsearch` | Same HMM as the preceding round (Codex agreed). Saves one full-DB pass per family. |
| **No** tie-break sort | Sort records by `(evalue, -score, name, env_from)` | **Measured**: native order already stable across `cpus` (F4); the sort was O(N log N) on ~10⁶ records for nothing. |
| Copy hits into plain tuples at extraction | Cache `Domain` objects | A `Domain` pins its `Hit`, which pins the entire `TopHits` (F10). |
| **Keep legacy raw-`TopHits` iteration** | Iterate `.reported` / `.included` | Dropping sub-threshold recruits is a scientific change nobody requested. Flagged for the user (F3). |
| Insertion-ordered dict for clusters | `itertools.groupby` | `groupby` only groups *contiguous* rows; legacy gathers all rows per representative. |
| Discard `<= 2` after trimal, `< 2` before initial/hand builds | `<= 2` everywhere | `<= 2` on the initial cluster would newly discard 2-member clusters — not requested. |
| Let `SSIWriter` detect duplicate names | Pre-scan into a Python `set` | A `set` of a billion keys is the database-sized Python state this port exists to delete (Codex #8). |
| `if fetched.name != name: raise` on every SSI lookup | `assert`, or a fingerprint sidecar | `python -O` strips `assert` (verified). O(1), catches a mismatched index on first use; the sidecar is filesystem-integrity engineering (Codex #9). |
| Build the SSI inside a private `mkdtemp` dir, `rmtree` on error | `<ssi>.tmp.<pid>`, or `mkstemp` + unlink-on-error | `close()` raises before `esl_newssi_Close`, leaving native scratch files open; only a whole-directory `rmtree` cleans up (Codex #2/#10). |
| `--batch_size` default `2 * cpus` | `max(2 * cpus, 16)` | At `cpus=1` the floor pins 16 million-record lists for no parallelism (Codex #5/round 3). |
| Daemon timer-thread heartbeat | `hmmsearch(callback=...)` alone | The callback fires only after a query *completes* — silent for hours on a big database (Codex #3/round 3). |
| `creation_time = None` | Fixed `datetime(1970,1,1)` | Omits `DATE` entirely instead of formatting it through locale-dependent `%a`/`%b` (Codex #11). |
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

## Settled decisions (no longer open — Codex #15)

- **Single-pipeline (`cpus=1`) semantics are the reference.** The port targets
  legacy-at-`cpus=1` and is `cpus`-invariant. It will **not** reproduce legacy
  production output byte-for-byte at `cpus > 1`. Accepted: the alternative is to
  keep a result that depends on the machine's core count. On the 50 000-sequence
  fixture the divergence does not reach the final families (F3b), so the blast
  radius is expected to be small.
- **`clip_env_ends` changes seed MSAs**, hence downstream HMMs, recruitment,
  convergence and family counts. Requested by the user; accepted.
- **Raw `TopHits` iteration is preserved**, sub-threshold recruits and all.

## Risks / open questions

- **Legacy recruits hits that fail `--recruit_evalue_cutoff`** (F3): 55 stored vs
  54 reported for one fixture family even at `cpus=1`. Preserved deliberately, but
  the user may want `.reported` iteration in a follow-up. **This is a scientific
  bug in the original, surfaced here, not fixed here.**
- **`SequenceFile._FORMATS` is private API.** Working (F1), but pinned to
  pyhmmer's internals. A regression test covers it.
- **Billion-record SSI construction is unverified at scale** (Codex #7). Easel
  buffers up to `eslSSI_MAXRAM` (2048 MB) of keys before spilling to an external
  sort, which then needs scratch space in `TMPDIR` plus room for the final index.
  The 50 000-record measurement (F1) says nothing about that transition. Before a
  production run: size `TMPDIR`, and budget ~`n_keys × (keylen + 16)` bytes for
  the index. Documented in the README, not solved here.
- **The streaming tax is `n_queries` full database parses per batch.** Mitigated by
  `--prefetch_targets` when RAM allows. A peak-RSS and wall-clock benchmark against
  a database larger than page cache, on the real shared storage, is a prerequisite
  for the production rollout and is deliberately not attempted on a 10 MB fixture.
- **Per-family memory** is driven by a single family's hit count (~10⁶ records),
  which `--batch_size` cannot bound. Recorded, not solved. A peak-RSS benchmark on
  representative data is deferred to the pipeline integration.

## Out of scope

- The `mgnifams` Nextflow module (`main.nf`, `meta.yml`, nf-test snapshots).
- Changing E-value / length-percentage defaults or any scientific threshold.
- Checkpoint/resume; multi-node distribution.
- The two "preserved as-is" quirks.

---

## Verification

0. **The byte-reproducibility contract is scoped to the committed `uv.lock`.**
   pyhmmer/pyfamsa/pytrimal determine hit retention, alignment and serialised bytes,
   so a promise that spans every future release is unbackable (Codex #7/round 4).
   `pyproject.toml` carries compatible-release upper bounds, CI runs
   `uv sync --frozen`, and the README states the contract in those terms.
1. `uv lock --check && uv sync --frozen && uv run pre-commit run --all-files` — clean.
   `--frozen` alone does not verify that the lock matches project metadata (Codex, round 5).
2. `uv run pytest`:
   - **`build_ssi_index` round-trip** — every FASTA name fetchable, sequences equal
     to a naive parse; missing key → `KeyError`; duplicate names → clear error with
     no leaked temp index (F1, F2); a mismatched index → name-assertion failure
     (Codex #9); `SequenceFile._FORMATS["fasta"]` is an int (Codex #7).
   - **`clip_env_ends`** — RF `"...xx.xxx..."` → leading/trailing dot runs removed,
     interior dot kept (F8).
   - **`filter_hits`** — `exit_flag=True` admits short envelopes that
     `exit_flag=False` rejects; record order is the input `TopHits` order.
   - **small-MSA discard** — 2-sequence alignment out of `run_pytrimal_reps` →
     discarded, `hmmbuild` never called, `eslEMEM` never raised.
   - **`--cpus` invariance** — end-to-end at `cpus=1` and `cpus=4` → byte-identical
     scientific outputs (the hazard F3 exposes).
   - **batch invariance** — end-to-end at `batch_size=1` and `batch_size=64` →
     byte-identical outputs, *and* an identical `representative → family_id`
     mapping (Codex #5: contiguity alone would accept ids on the wrong reps).
   - **converge-then-discard** — drives a **real `Family`** through the actual
     lifecycle: `advance()` to convergence, `finish()` into a membership discard,
     then `emit()`. Asserts `ever_converged` survives the transition to `DISCARDED`,
     that `converged_families` receives the discarded family's *provisional* id, and
     that the next successful family **reuses** that id. Testing `emit()` arithmetic
     on synthetic families would pass even if `finish()` cleared the marker
     (Codex #2/round 4). Pins the legacy interleaving at
     `reference/legacy_generate_families.py:541` against the exit-branch discards at
     `:556-577` and the `iteration -= 1` at `:596` (Codex #7/round 2).
   - **determinism** — two runs, different `PYTHONHASHSEED`, different output dirs
     → every scientific artifact byte-identical; `logs/` excluded (Codex #13).
     Asserts the HMM has no `DATE`/`COM` line (F6) and that the gzip header carries
     no filename (Codex #14).
   - **rerun idempotence** — a second run into the *same* output directory that
     yields fewer families leaves no orphaned `<chunk>_<n>` artifacts (Codex #12).
   - **raw-vs-reported extraction** — on the small fixture, query
     `4497037939_1_144` must recruit sequence `6320430079`, which is stored but
     *not* `.reported` (F12). Fails the moment someone "fixes" extraction to
     iterate `.reported`. The end-to-end test cannot catch this (F3b).
   - **`--prefetch_targets` equivalence** — two levels. (a) Unit: one HMM at
     `cpus=4`, streaming vs prefetched, compare the **raw** `TopHits` signature
     (name, score, envelope tuples), asserting the `55`/`54` stored-vs-reported split
     and the presence of `6320430079`. (b) End-to-end with and without the flag →
     byte-identical scientific outputs (F11). The end-to-end test alone is blind here,
     because F3b shows the extra raw hit is filtered downstream.
   - **SSI guard survives `-O`** — the index-mismatch check must raise under
     `python -O`, where `assert` is stripped (Codex #1/round 3).
   - **no in-memory FASTA** — monkeypatch `SequenceFile.read_block` to raise; the
     end-to-end run must still pass with streaming (and must *not* be patched for
     the `--prefetch_targets` case).
   - **`Domain` objects are not retained** — after a batch, `f.records` holds only
     tuples (`gc.get_referents` finds no `TopHits`) (F10).
   - **input validation** — `../evil` chunk, gzipped FASTA, `cpus=0`,
     out-of-range percentages, inverted length bounds, malformed TSV rows all
     rejected before any output directory is created (Codex #14).
   - **`hmmalign` receives `cpus`** — asserted via monkeypatch; its default of `0`
     expands to every physical core (Codex #3/round 4).
   - end-to-end on `clustering.tsv` + `mgnifams_input_small.fa`, and on
     `cluster_long.tsv` + `mgnifams_extra.fa` (>200 aa): exit 0, every declared
     output dir exists, `hmm/*.hmm.gz` reloads via `HMMFile`,
     `family_metadata/*.csv` has one row per successful family.
     **Correction, found during implementation:** neither MSA carries `#=GF ID`.
     Legacy's `renumber_sto_msa` drops *every* `#=GF`/`#=GS`/`#=GR` line, including
     the `#=GF ID` that naming the seed MSA emits — verified against a real legacy
     run. The family name survives only as the HMM's `NAME` field. So: assert
     `hmm/*.hmm.gz` has `NAME <chunk>_<id>`, and that both `.sto.gz` files parse and
     contain no `#=GF`. Codex's round-2 finding #6 was right that `seed_msa.name`
     must be set before the hand build (it becomes `HMM.name`), but wrong that it
     reaches the Stockholm output; I accepted it too readily, and only building it
     exposed the error.
3. **Legacy comparison is informational, not an assertion** (Codex #19). The
   baseline exists: `reference/legacy_generate_families.py` runs under a pinned
   `pyhmmer==0.11.1 / pandas==2.3.2 / pyfamsa==0.6.0 / pytrimal==0.8.2` venv.
   Run it at `cpus=1` on `clustering.tsv` and record the deltas in
   `refined_families/`, `successful_clusters/`, `family_metadata/` and the seed
   `#=GF ID`. Expected sources of divergence, each to be explained: `clip_env_ends`,
   the `<= 2` discard rule, gzip framing, the omitted `DATE`/`COM` lines, and
   dependency upgrades. Any *other* delta is a port bug.

---

## Progress

Branch `dev`. **One feature, one commit, no push.**

### Act 2 — Codex adversarial review — **complete**
- [x] Freeze `PLAN.md` + init `PLAN-REVIEW-LOG.md`
- [x] Round 1 → `REVISE` (21 findings; 16 accepted, 3 with modification, 2 rejected with measurements)
- [x] Round 2 → `REVISE` (15 findings; all accepted — my `domZ` causal claim refuted by measurement)
- [x] Round 3 → `REVISE` (10 findings; 9 accepted, 1 narrowed — `--prefetch_targets` born here)
- [x] Round 4 → `REVISE` (7 findings; all accepted — 2 were live bugs in the pseudocode)
- [x] Round 5 → **`VERDICT: APPROVED`** (6 non-blocking corrections, all applied)

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
- [x] `uv run pre-commit run --all-files` clean
- [x] `uv run pytest` green — 21 passed
- [x] SSI mismatch guard raises under `python -O`
- [x] End-to-end invariance on the 50 000-sequence fixture, scientific artifacts
      byte-identical across: `cpus=1` vs `cpus=4` (different `PYTHONHASHSEED`),
      streaming vs `--prefetch_targets`, default vs `--batch_size 4`
- [x] Legacy-vs-new delta report (below)

### Legacy-vs-new deltas (legacy at `cpus=1`, `clustering.tsv`)

| output | result |
|---|---|
| `successful_clusters/` | **identical** |
| `converged_families/` | **identical** |
| `family_metadata/` | same 3 families, same ids, same representatives (`782510898`, `5761513631`, `1446399400`), same `converged` flags; rep regions shift |
| `refined_families/` | 115 → 116 rows; same set of member proteins |
| `rf/` | RF lengths 116/109/121 → 115/104/119 |
| `seed_msa_sto/`, `full_msa_sto/`, `hmm/` | gzipped; HMMs additionally omit `DATE`/`COM` |

Every delta traces to `clip_env_ends`: family 2's representative region moves
`863-967` → `864-967` (one N-terminal envelope column), family 3's `2-120` → `4-115`.
Family 1 gains one recruited region because a clipped seed changes its round-2 HMM
and therefore its recruitment — the documented cascade. No unexplained delta.
