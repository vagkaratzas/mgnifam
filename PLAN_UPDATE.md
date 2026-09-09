# Plan: `mgnifam update_families`

## Goal

Add a third subcommand, `mgnifam update_families`, that refreshes families which already
exist as HMMs against a newer, larger protein database — without re-deriving them from
their original MMseqs2 clusters.

The concrete driver: MGnifams were produced against the 700M-sequence MGnify Proteins
release by the legacy implementation. A 1.6B-sequence release now exists. Those families
need to recruit from the bigger database and pick up the two recruitment defects fixed in
1.0.0 (`extract_records` reading `top_hits.reported`; `clip_ends` keeping its last passing
column), both of which live in the shared code path this subcommand reuses.

### Why not re-run `generate_families` from a reconstructed cluster TSV

That route was evaluated and rejected. It fails *quietly*, which is worse than failing:

- `<chunk>_families.tsv` records members as `<protein>/<start>-<end>` in parent-protein
  coordinates (`generate_families.py:945`, via `parse_protein_name`). The FASTA holds
  fragments named `<protein>_<start>_<end>`. The two namespaces do not intersect.
- `run_initial_msa` skips members absent from the FASTA silently (`missing_ok=True`,
  `generate_families.py:254`). A family whose members all fail to resolve becomes an empty
  seed and is discarded as "too few sequences before initial hmmbuild" — a run that
  completes, reports discards, and explains nothing.
- `extract_first_part` splits on `/` (`generate_families.py:513`), so `A/1-100` collapses
  to `A`. Two fragments of one protein count once in `check_seed_membership`'s numerator
  and twice in its denominator.
- The denominator is `len(self.members)`, so members that vanished between releases count
  against the family regardless of how well the survivors align.

The stored HMMs are exactly the state that route spends its first round reconstructing.

## Acceptance criteria (observable)

1. `mgnifam update_families --hmm_input <dir|file> --fasta_file <new>.fa` completes and
   writes the artifact set in "Outputs" below.
2. `--hmm_input` accepts a directory of `.hmm`/`.hmm.gz` **and** a single multi-model
   library (`hmm.lib`, `hmm.lib.gz`). Both forms produce byte-identical outputs for the
   same set of models.
3. Family identity is preserved: every successful family's artifacts **and every
   identity-bearing field inside them** — `metadata.csv`'s `family_id`, `families.tsv`'s
   first column, `converged.txt`, and the `reps.fasta.gz` annotation — carry the HMM's
   `NAME`, not a rank. No renumbering, no mapping file.
3b. A family name that is not a safe identifier is rejected before any output exists, and
   every artifact path written resolves inside its intended directory.
4. `--skip_refine` performs exactly one `hmmsearch` pass per family; the refine path
   performs at most three. Assertable by counting `search()` invocations under monkeypatch.
5. `generate_families` output is byte-identical before and after the edits this plan makes
   to `generate_families.py`, verified against a committed SHA256 manifest.
6. The existing reproducibility contract extends to the new subcommand: identical output
   across `--cpus`, `--batch_size`, `--prefetch_targets` and `PYTHONHASHSEED`.
7. The exit-status contract (0 complete / 1 corrupted, do not consume / 2 usage / 3
   completed with contained internal family errors) holds unchanged.
8. `uv run pytest` green; `uv run pre-commit run --all-files` clean.

## Confirmed assumptions (verified in this repo, not recalled)

| Claim | Source |
|---|---|
| `pyhmmer.plan7.HMMFile` and `pyhmmer.easel.MSAFile` open `.gz` paths directly | executed against `output/hmm/1_1.hmm.gz` and `output/seed_msa/1_1.sto.gz` under `uv run`; both parsed, `reference` present on the MSA |
| A stored seed `.sto.gz` carries `#=GC RF` | `zcat output/seed_msa/1_1.sto.gz \| tail -5` |
| Legacy `renumber_sto_msa` also preserved the RF line, so legacy-produced seeds are equally usable | `reference/legacy_generate_families.py:437` |
| Legacy wrote uncompressed `.hmm`/`.sto` with `DATE`/`COM` lines | `reference/legacy_generate_families.py:602-611` |
| **Round 1 can never converge.** `Family.total_checked_sequences` starts empty, so the first `advance()` always adds names | `generate_families.py:641,703-709`; `emit_family`'s "round 1 can never converge" comment at `:877` already relies on it |
| Therefore refine mode does not need the stored seed MSA — the loaded HMM is the round-1 query and `advance()` builds the next seed. **Input contract is HMM-only.** | follows from the above |
| **`HMM.name` is `str`, not `bytes`**, under the locked pyhmmer 0.12.1 | `.venv/lib/python3.13/site-packages/pyhmmer/plan7.pyi:364-368` declares `name -> str`; executed against `output/hmm/1_1.hmm.gz` returns `'1_1'`, a `str`. `HMM.consensus` is likewise `str`. |
| **A nameless HMM cannot be parsed at all**, so `NAME` is always present | stripping `NAME` from a real model and reopening raises `ValueError: Invalid format in file: No NAME found for HMM` |
| `PLAN.md` / `PLAN-REVIEW-LOG.md` are historical for the original port and must not be edited | `AGENTS.md:117` |
| Subcommands own their parsers, expose `main(argv)`, register in `COMMANDS`, and set `prog="mgnifam <name>"` | `AGENTS.md:106-111`, `cli.py:19-21` |

Remaining risk, not verifiable here: the 1.6B release must keep the
`<protein>_<start>_<end>` fragment-naming convention for `parse_protein_name` to report
true parent coordinates. `split_slice_name`'s span test (`generate_families.py:557`)
degrades safely — a non-conforming name is treated as a whole protein rather than given
invented coordinates — so this is a silent-precision risk, not a crash.

## Key decisions

| Decision | Resolution | Trade-off accepted |
|---|---|---|
| Mode selection | `--skip_refine` (store_true). Set: one search, one `hmmalign`. Unset: the full three-round loop. | A boolean rather than `--mode`; adding a third behaviour later needs a flag change. |
| Membership yardstick | Round 1's length-filtered recruits become `family.members`. `finish()` then applies `--discard_min_starting_membership` to them. | Under `--skip_refine` retention is 1.0 by construction (see OQ2), so the check is a structural no-op there. No old-membership file is needed anywhere, and both sides of the comparison are fragment names from the *new* database, so `extract_first_part` behaves correctly. |
| Family identity | Preserved verbatim from the HMM's `NAME`. | `metadata.csv`'s `family_id` column carries `1_7`, not a bare integer — a documented format difference from `generate_families`. |
| Name source | HMM `NAME` for both input forms. No filename fallback — a nameless model cannot be parsed (see assumptions), so the fallback would be unreachable. | One rule for both forms. A renamed file does not change a family's identity. |
| Name validation | A family name must match `[A-Za-z0-9._-]+` and be neither `.` nor `..`. | That alphabet is `CHUNK_PATTERN`'s (`generate_families.py:62`), already the trusted alphabet for the other half of every artifact filename. It excludes `/`, `,`, whitespace and control characters in one rule, which closes both the path-traversal and the CSV-corruption vector below. |
| Load order | After loading and validating, models are sorted by family name, for **both** input forms. | Filename order and library order are different orders over the same models; sorting by the identity key is the only rule that makes criterion 2 true. |
| Aggregate output names | `<chunk>_updated_*`; `--chunk_num` labels aggregates only, since per-family names come from the input. | Parallel update jobs can share an output root without colliding. |
| Code reuse | Additive, behaviour-preserving parameters on `generate_families.emit_family` plus one opt-in `Family` flag; `update_families.py` imports the helpers and owns its own `main()` loop. One change is **not** default-preserving and is declared as such: the discard branch gains the `ChunkCorrupted` boundary the success branch already has. | See OQ1. Extract a shared `_pipeline.py` only when a third command (`remove_redundant`, `merge_families`) lands and the shared shape is observed rather than guessed. |
| gz inputs | HMM: yes. FASTA: **no**. | Easel cannot seek within a gzip stream, so an SSI index cannot be built over a compressed FASTA. Identical to the constraint `validate_inputs` already enforces (`generate_families.py:1018`). |

## Approach

### CLI (`src/mgnifam/update_families.py`, `parse_args`)

```
mgnifam update_families --hmm_input <dir|file> --fasta_file <new release .fa> [--skip_refine]
```

`prog="mgnifam update_families"`. Every threshold flag from `generate_families` carries
over with identical names and defaults: `--cpus`, `--chunk_num`, `--discard_min_rep_length`,
`--discard_max_rep_length`, `--discard_min_starting_membership`, `--max_seq_identity`,
`--max_seed_seqs`, `--max_gap_occupancy`, `--recruit_evalue_cutoff`,
`--recruit_hit_length_percentage`, `--fasta_index`, `--output_dir`, `--batch_size`,
`--prefetch_targets`.

`--max_seq_identity`, `--max_seed_seqs` and `--max_gap_occupancy` are validated normally
but inert under `--skip_refine` (they only act inside `advance()`'s trim tail). Document
that in `--help` and the README rather than silently ignoring them.

Register in `cli.py`'s `COMMANDS`.

### Input loading (`load_hmms`)

- Directory → `*.hmm` + `*.hmm.gz`, one `HMMFile` each.
- File → one `HMMFile` iterated to exhaustion, so a multi-model library works unchanged.
- Family name = `hmm.name`, used directly. **It is already a `str`** under the locked
  pyhmmer 0.12.1 — `.decode()` would raise `AttributeError` on every model. There is no
  filename fallback: a model with no `NAME` cannot be parsed at all, so the fallback
  would be unreachable code.
- **Names are then sorted**, for both input forms alike. Filename order and library order
  are different orders over the same models — `a.hmm` holding `NAME B` and `z.hmm` holding
  `NAME A` load as `B, A`, while a library holding them in `A, B` order loads as `A, B`,
  and the aggregate files would differ byte-for-byte. Sorting by the identity key is what
  makes acceptance criterion 2 true, and it also stops a file rename from reordering
  output when identity is supposed to come from `NAME`.
- Fatal before any output directory is created, matching `validate_inputs`' existing
  ordering (`generate_families.py:1012`, guarded by `test_validation_precedes_output_creation`):

  - an empty input;
  - a duplicate family name — an id collision would silently overwrite artifacts;
  - **a family name that is not a safe identifier.** Names arrive from an arbitrary
    third-party file and are interpolated straight into filesystem paths and into unquoted
    CSV fields. A `NAME` of `../../victim` escapes the artifact directory — and the
    rollback `unlink` would then delete outside it too; a `NAME` containing a comma
    silently shifts every column of `metadata.csv` and `discarded.csv`. Requiring
    `CHUNK_PATTERN`'s alphabet, `[A-Za-z0-9._-]+`, and rejecting the two traversal names
    `.` and `..`, closes both at once: that alphabet admits no `/`, `,`, whitespace or
    control character. As defence in depth, every artifact path is additionally asserted
    to resolve inside its intended directory before being written.
  - a model whose `alphabet` is not `ALPHABET` — a nucleotide model otherwise dies inside
    `hmmsearch` with a raw Easel message;
  - an `--output_dir` that resolves to the input HMM location.
- Models are held in memory for the whole run. An HMM is kilobytes; tens of thousands of
  them is tens of megabytes, which is negligible beside the database. This is also what
  makes the output-dir guard sufficient rather than merely advisory.

### Per-family flow

**`--skip_refine`** — one database pass:

```
search([hmm], targets)  ->  extract_records  ->  filter_hits(exit_flag=False)
   ->  family.members = unmask_sequence_names(filtered)  ->  family.finish(...)
```

`finish()` is reused verbatim. It re-filters the same cached records with
`exit_flag=True`, so the final full MSA admits partial-envelope hits exactly as
`generate_families`' exit branch does — the length filter never restricted the final full
MSA there either.

Two distinct empty-recruitment discards, because for an update run the difference between
them is the whole point of the delta report:

- `records == []` — the old model finds nothing at all in the new release → new reason
  `"no hits in the new database"`.
- records exist but none clear the envelope-length filter → the existing
  `"low complexity model - confounding cluster"`.

**Refine** — up to three database passes. Round 1 searches with the loaded HMM instead of
building one from a seed; `advance()` adopts its length-filtered recruits as `members` and
builds the next seed; rounds 2–3 and `finish()` run unchanged. Because round 1 cannot
converge (see assumptions), a seed always exists by the time `emit_family` runs on this
path.

### Outputs under `--output_dir`

| path | `--skip_refine` | refine |
|---|---|---|
| `hmm/<name>.hmm.gz` | the loaded model, re-serialised with `creation_time`/`command_line` cleared | hand-built from the final seed, as today |
| `full_msa/<name>.sto.gz` | yes | yes |
| `seed_msa/<name>.sto.gz`, `rf/<name>.txt` | **not written** (see OQ3) | yes |

Re-serialising rather than byte-copying the input HMM is deliberate: a legacy `.hmm`
carries `DATE` and `COM` lines, which are a wall clock and an `argv` dump. Clearing them is
what makes the output reproducible, and the parsed model is the same model that was
searched with.

Per-chunk aggregates, all `<chunk>_updated_*`: `families.tsv`, `metadata.csv`,
`discarded.csv`, `successful.txt`, `converged.txt`, `reps.fasta.gz`, `delta.csv`, and
`<chunk>_updated.log`. Both CSVs keep their header-before-the-run behaviour so a chunk
producing no families still parses.

`<chunk>_updated_delta.csv` is the run's actual scientific deliverable — one row per
family, successful or discarded:

```
family_id,model_length_before,model_length_after,round1_recruits,full_msa_size,retention,rounds_run,converged,outcome
```

`outcome` is `successful` or the discard reason. Under `--skip_refine`,
`model_length_before == model_length_after` and `rounds_run == 1`.

**Every field except `family_id`, `model_length_before` and `outcome` is nullable**, written as
an empty CSV field, because a family discarded early never reaches the stage that would produce
one. Two specific cases, both of which the naive schema got wrong:

- `retention` is empty when round-1 recruitment was empty. There is no yardstick in that case,
  and `check_seed_membership` divides by `len(original_first_parts)` (`:531`) — calling it with
  no members raises `ZeroDivisionError`, not a zero score.
- `model_length_after` and `full_msa_size` are empty for any discard that never built a final
  model or full MSA.

**Metrics are captured as they are computed, never read back off the `Family`.** `discard()`
deliberately clears `records`, `hmm`, `seed_msa`, `full_msa` and `total_checked_sequences`
(`:652-656`) to release memory, and `finish()` holds `membership` only in a local before
possibly discarding on representative length two statements later (`:753-765`). A family
rejected for length has computed its retention and then destroyed the evidence. So
`update_families` threads a small per-family delta record through the loop and fills each scalar
at the point of computation.

**The delta row is committed on both of `emit_family`'s exit paths**, once the family's outcome
is final, and under a `ChunkCorrupted` boundary on each. A delta write that failed after the
family was already recorded would leave the run's own summary contradicting its output, which is
the definition of incoherent output under the exit-1 contract. It must be fatal, not contained
as another discard.

The success path already has such a boundary (`:936-961`) and the delta write joins it. The
discard path does **not**: it appends to `discarded.csv` and returns at `:868-872`, before that
`try` exists. Placing the delta write only in the success block would silently omit every
discard — including the no-hit and length-discard rows this report exists to carry.

This exposes a pre-existing hole in `generate_families`, and closing it is a deliberate,
declared widening of this plan's scope rather than an additive change:

> The 2.1.0.dev0 hardening — *"Shared-output commit failures now bypass per-family containment
> and exit 1, so exit 3 structurally means the degraded output is coherent"* (`CHANGELOG.md`) —
> protected the success path only. The discard branch's single shared append at `:869` is still
> unprotected, so a failure there is caught by `family_guard` in `main` (`:1256-1266`), which
> re-emits, which calls `emit_family` again, which appends the discard row a second time. A
> partial first write followed by a successful retry duplicates the row. That is the same class
> of defect the release already claims to have fixed, in the one branch it missed.

So the discard branch gains the same `ChunkCorrupted` boundary as the success branch, covering
its `discarded.csv` append and its delta row together, and returning `success_count` unchanged.
This **does** change `generate_families`' behaviour in one respect — a failed `discarded.csv`
append now exits 1 instead of being contained and retried — so it is not covered by the
default-preserving claim below, it carries its own `CHANGELOG.md` *Fixed* entry, and its own
test. It does not affect output bytes on any clean run, so acceptance criterion 5 still holds.

Discarded rows carry the old family name in `discarded.csv`'s `representative` column;
there is no MMseqs2 representative in an update run.

### Edits to `generate_families.py`

Under `--skip_refine`, `family.seed_msa` stays `None` because `advance()` never runs.
Today's `emit_family` would then fail four ways: the hand `run_hmmbuild` (`:880`), the
`seed_msa.reference is None` check (`:888`), `renumber_msa(seed_msa.textize())` (`:891`),
and the `rf/` and `seed_msa/` writes (`:906,915`). The reuse therefore needs three
keyword-only parameters, each defaulting to exactly today's behaviour:

1. `family_name=None, family_id=None` → default to `f"{chunk}_{success_count + 1}"` and
   `success_count + 1`.
2. `final_hmm=None` → when supplied, use that model instead of the hand `hmmbuild`, still
   clearing `creation_time` and `command_line`.
3. `family.seed_msa is None` → skip the RF check and the `rf/`+`seed_msa/` writes. The
   rollback list simply holds fewer paths; the ordering guarantee (everything that can
   raise happens before the first shared append) is unchanged.

Two further sites the preserved names force, both inside `emit_family`:

- **The return value.** The success path ends `return provisional_id` (`:962`), which with
  `family_id="1_7"` would put a string into `success_count` arithmetic. That *one* return
  becomes `return success_count + 1` — identical for the existing caller, where
  `provisional_id == success_count + 1` by construction (`:867`). The discard path's early
  `return success_count` (`:872`) is correct as it stands and is not touched.
- **The representative FASTA annotation.** `:956` writes
  `f">{sequence_name}\t{chunk}_{provisional_id}\n"`, rebuilding the family name from
  `chunk` and the id rather than using the resolved name. With `--chunk_num 9` and a
  preserved `1_7` it would emit `9_1_7`, so identity would be preserved in the filenames
  and broken in the one file downstream reads to map a representative back to its family.
  It must write the resolved `family_name`, whose default already reproduces
  `f"{chunk}_{success_count + 1}"` exactly.

4. `Family` gains an explicit opt-in flag (a dataclass field, not an implicit
   `if not self.members`) causing `advance()` to adopt round 1's length-filtered recruits
   as `members`. `generate_families` never sets it.

`update_families` owns its own `main()` loop — the control flow genuinely differs (HMM-driven
round 1, optional single-pass path) — and its own stale-artifact clearing.
`prepare_output_directories`' `<chunk>_\d+\..*` regex (`:1061`) does not describe preserved
names, and clearing only the names in the *current* input is not enough either: a chunk that
updates `A` and `B`, then reruns with only `A`, would leave `B`'s artifacts on disk while the
aggregates that reference them are overwritten without `B` — the exact mixed-two-runs state
`prepare_output_directories` exists to prevent (`:1051-1056`).

`<chunk>_updated_successful.txt` is **not** a usable ownership record, despite listing exactly
the families that wrote artifacts on a *clean* run. `emit_family` writes every per-family
artifact before its first shared append (`:904-934` then `:936`), so a `ChunkCorrupted` on that
append — or a `SIGKILL` before the buffered handle flushes — leaves artifacts on disk with no
name in that file. Those are precisely the runs that exit 1 and mandate a rerun, so the cleanup
would be blind exactly when it is needed.

The record must therefore be the *intended* ownership set, written before anything else:
`<chunk>_updated_manifest.txt`, one family name per line, flushed and `fsync`ed at startup. It
is known up front — it is the validated input name list — and every artifact the run can ever
write bears a name from it, whether or not the run completes.

The clearing set is the **union of the previous manifest and the current input's names**, sorted.
The manifest is written as that union *before* clearing begins, so an interruption part-way
through the clear still leaves a record covering everything already removed and everything not
yet removed. It therefore only grows for a given output root; that is one line per family name
and is documented rather than optimised. Clearing the whole root instead is not an option — it
would delete a concurrent chunk's output.

**The replacement must be atomic.** Opening the existing manifest for writing and being
interrupted after truncation destroys the previous ownership record before cleanup has used it,
which is worse than having no manifest at all — a later shrinking rerun could then never
discover those families. The repository already has the idiom: `build_ssi_index` writes into a
private temporary directory on the destination filesystem and `os.replace`s it into place
(`generate_families.py:160-196`), for the same reason. The manifest follows it — write the
sorted union to a temporary file beside the destination, `flush`, `fsync`, then `os.replace`;
`fsync` the parent directory as well where power-loss durability is wanted.

This carries a documented constraint rather than a guarantee: **chunks sharing one output root
must own disjoint family names.** Nothing enforces it, because family names come from the input
models, not from `--chunk_num`; violating it makes two chunks fight over the same artifact
paths. It is stated in the README beside the flag.

`family_guard`, `ChunkCorrupted`, `search`, `extract_records`, `filter_hits`,
`run_hmmalign`, `clip_env_ends`, `run_pytrimal_reps`, `clip_ends`, `renumber_msa`,
`msa_stats`, `check_rep_length`, `IndexedSequences`, `resolve_index`, `build_ssi_index`,
`deterministic_gzip_*` and `Writers` are imported unchanged.

## Open questions for the reviewer

**OQ1 — additive parameters vs. a duplicate emit.** The alternative to editing
`emit_family` is `update_families` owning its own, leaving `generate_families.py`
untouched. That costs a second copy of the artifact-rollback and `ChunkCorrupted`
contract — the repo's most invariant-heavy logic (`generate_families.py:832-962`,
`AGENTS.md:44-52`), and the one place where two copies drifting apart produces output that
is neither generated nor discarded. The additive route is recommended. Challenge it.

**OQ2 — the membership yardstick reading.** The user's words were that round-1 recruits
"become the cluster members, and they should be retained in the following rounds." This
plan reads *retained* as **measured**: round-1 recruits are the yardstick
`check_seed_membership` scores the final family against, the same role MMseqs2 members
play in `generate_families`. The alternative reading is **pinned**: those sequences are
forced into every later seed and full MSA whether or not the refined model still recruits
them — a different algorithm, and one that would make a family unable to shed a bad
round-1 recruit. The plain reading is taken. Flag it if the reviewer sees the other as more
likely from the wording.

**OQ3 — no `seed_msa/` or `rf/` under `--skip_refine`.** A recruit-only update does not
change the seed alignment, and a seed cannot be reconstructed from an HMM, so those two
artifacts are simply absent and the user keeps the originals. That leaves an output
directory that is not a drop-in replacement for a `generate_families` output tree, which
downstream tooling may assume. The alternatives are a `--seed_msa_dir` copy-through flag
(more surface for one use) or emitting a one-row seed from the model consensus (a
fabricated alignment). Absence is recommended as the honest option; the reviewer should
weigh the downstream-compatibility cost.

## Non-goals

- Re-deriving a family whose seed drifted. That still needs `generate_families` over the
  original MMseqs2 clusters against the new FASTA.
- A resume / skip-existing flag. Sharding belongs to the workflow manager, and the
  HMM-library input already lets a Nextflow task carry its own subset.
- `remove_redundant` and `merge_families`, and the shared `_pipeline.py` extraction their
  arrival should trigger.
- Any change to `generate_families`' scientific behaviour. Every edit above is additive and
  default-preserving, and criterion 5 exists to prove it — with one declared exception that is
  not scientific: a failed `discarded.csv` append now exits 1 rather than being contained and
  retried, closing the branch the 2.1.0.dev0 shared-output hardening missed.
- Editing `PLAN.md` or `PLAN-REVIEW-LOG.md`.

## Toolchain

Python >= 3.13, `uv sync --frozen` against the committed lockfile. **No new dependency**
(`AGENTS.md:12`): the four scientific packages are upper-bounded because they decide hit
retention, alignment and serialised bytes. Everything this plan needs — multi-model
`HMMFile` iteration, gz-transparent opening — is already available in the pinned pyhmmer.

## Verification

### Proof commands

```bash
uv lock --check && uv sync --frozen
uv run pre-commit run --all-files
uv run pytest
```

Expected: lock unchanged, hooks clean, and the suite green at 37 existing tests plus the
new ones below.

### The byte-identity regression check for criterion 5

Stated as a procedure so it is reproducible rather than aspirational:

```bash
git worktree add /tmp/mgnifam-baseline <pre-edit-commit>
cd /tmp/mgnifam-baseline && uv sync --frozen
uv run mgnifam generate_families -c tests/fixtures/mgnifams_v2.tsv \
    -f <decompressed v2 fasta> -n baseline --output_dir /tmp/baseline-out
find /tmp/baseline-out -type f ! -name '*.log' ! -name '*.ssi' -printf '%P\n' \
    | sort | xargs -I{} sha256sum /tmp/baseline-out/{} > tests/fixtures/generate_families_manifest.txt
```

The manifest is committed as a fixture. A new test re-runs the same invocation on the
edited tree and compares. `<chunk>.log` is excluded (timestamps) and the `.ssi` index is
excluded (it embeds the FASTA's filename).

### New tests — `tests/test_update_families.py`

Reusing the existing session fixtures from `conftest.py`: generate families on
`small_fasta`, then update them against `extra_fasta` / `v2_fasta` as the "new release".

- `--skip_refine` writes an HMM whose parsed model matches its input (`M`, `consensus`,
  `NAME`) and writes no `seed_msa/` or `rf/` files.
- `--skip_refine` invokes `search()` exactly once per family; refine at most three times.
- Family names are preserved everywhere, not only in filenames: the assertion covers
  `metadata.csv`'s `family_id`, `families.tsv`'s first column, `converged.txt` and the
  `reps.fasta.gz` annotation. Run with a `--chunk_num` that differs from the names' own
  prefix, so a `{chunk}_{id}` reconstruction would be visible as `9_1_7`.
- A multi-model `hmm.lib.gz` and a directory of the same models produce byte-identical
  outputs — including the adversarial ordering case: filenames whose sort order is the
  reverse of their `NAME` order, and a library whose model order is reversed.
- A family name containing `/`, `..`, a comma, or a control character is rejected before
  any output exists. The traversal case asserts that no file appears outside the output
  root; the comma case asserts `metadata.csv` still parses with the declared column count.
- Determinism across `--cpus`, `--batch_size` and `PYTHONHASHSEED`; `--prefetch_targets`
  equivalence — mirroring `test_batch_size_invariance`,
  `test_determinism_across_hash_seeds` and `test_prefetch_end_to_end_equivalence`.
- Duplicate `NAME`, an empty input, a nucleotide model, a gzipped FASTA, and an
  `--output_dir` equal to the input location each fail before any output exists.
- A shrinking rerun (chunk updates `A` and `B`, then reruns with only `A`) leaves no `B`
  artifacts behind, and a second chunk's artifacts in the same output root survive it.
- The same shrinking rerun after an *interrupted* run: force a shared-append failure so `B`'s
  artifacts land without a `successful.txt` entry, then rerun with only `A` and assert `B`'s
  artifacts are gone — the case the manifest exists for.
- Delta rows for a no-hit family, a representative-length discard, and an artifact-write
  failure each carry the fields reachable at that stage and leave the rest empty; the file
  parses with a constant column count throughout.
- A delta-write failure after the shared appends exits 1, not 3, on **both** the success and
  the discard path, and the discard is not re-emitted — exactly one row per family per file.
- Every discarded family produces exactly one delta row.
- An interruption during manifest replacement leaves the previous manifest intact and usable:
  a shrinking rerun after it still finds and clears the stranded artifacts.
- `"no hits in the new database"` and
  `"low complexity model - confounding cluster"` are reachable and distinct, and neither
  collides with `INTERNAL_ERROR_PREFIX` — mirroring
  `test_legitimate_discard_reasons_do_not_collide_with_internal_errors`.
- A family failing mid-update is contained as a discard and the chunk exits 3 —
  mirroring `test_one_failing_family_is_discarded_and_the_chunk_survives`.

### Docs that must move with the code

`README.md` (usage, flag table, outputs, the gz-FASTA rationale), `CHANGELOG.md` (a minor
release: new subcommand, plus the `metadata.csv` `family_id` format difference),
`AGENTS.md` (the "Adding a subcommand" list, and the new command's own invariants).
`AGENTS.md:114` requires these to move together.
