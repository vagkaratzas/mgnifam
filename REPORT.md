# Repository and PR #7 review

Reviewed 2026-09-10. The original review requested changes for two output-safety
defects. Follow-up fixes are tracked below; findings describe the reviewed commit.
Scope of follow-up: F1, F2, F5 and F6. F3 and F4 remain open for separate work.

## Scope and evidence

- [PR #7: Update families](https://github.com/vagkaratzas/mgnifam/pull/7) was open at
  review time. Its head and the clean local checkout both identified commit
  `45352269cf8786590fdbbd7cb3a8cec358c4cf90`; the PR base was
  `97655b5ed34ae0d882512888857d1d56f4974f59`.
- Reviewed all production Python modules, both test modules and fixtures setup,
  packaging and dependency declarations, CI/release workflows, README, CHANGELOG,
  AGENTS instructions, the complete production-code PR diff, and relevant legacy
  algorithm and historical review sections. Read the PR's inline discussions and
  [linked issue #8](https://github.com/vagkaratzas/mgnifam/issues/8).
- `uv lock --check && uv sync --frozen` passed using `UV_CACHE_DIR=/tmp/uv-cache`.
  `uv run pytest -q` completed with exit 0; collection confirmed **70 tests**.
- Additional checks used real v2 fixtures in temporary directories: repeated updates,
  a contained transient failure, input/output collision, renamed FASTA records,
  regenerated HMMs, and the documented v3 additions in both update modes. CSV checks
  exercised the real `emit_family()` with small alignments and in-memory writers.
- No production code, tests, dependencies, reference files, or GitHub state were
  changed. This report is the only repository addition. Full pre-commit hooks,
  distribution builds, and a separate legacy-environment replay were not run.
  GitHub returned no status-check results for the reviewed head; this report relies
  on the local verification, not a claim that remote CI passed.

P1 means fix before merging; P2 means a functional defect affecting supported input;
P3 means a maintenance correction. Findings distinguish PR defects from existing issues.

## F1 — P1: Allowed reruns leave artifacts for discarded families

**DONE.** Accepted retries clear all four exact artifact paths for each input family,
after collision validation and refusal of foreign families. Cleanup failures abort
before aggregate truncation. Regression tests cover changed discard outcomes,
transient failures, refine-to-recruit-only retries, and cleanup failure.
All 35 update-command tests and all pre-commit hooks passed after F1/F2.

**Introduced by PR #7.** Location:
[`update_families.py:241–255`](src/mgnifam/update_families.py#L241), interacting with
the discard return in [`generate_families.py:926–941`](src/mgnifam/generate_families.py#L926).

`prepare_output_directories()` refuses names outside the input set, but retains all
artifacts for names inside it. The aggregate files are then truncated for the new
run. If a previously successful family now discards, `emit_family()` appends its
discard row and leaves its old HMM and alignments on disk. The artifact rollback
also covers only files attempted during the current emission, not older files.

**Reproduced with the supplied 14-model library:**

1. Update the uncompressed v2 FASTA with `--skip_refine` into a fresh directory:
   14 HMMs, 14 successes.
2. Repeat with the same models and output directory, adding
   `--discard_min_rep_length 2000`: the call returns normally, there are 14 discard
   rows and zero successes, but **all 14 old HMMs and full alignments remain**.
3. Independently, repeat a successful update with unchanged scientific parameters
   and inject one `OSError` from the first `renumber_msa()` call. The run exits **3**,
   reports `v2_1` discarded and 13 successes, but still contains **14 HMMs**.

Thus this affects the explicitly supported retry/containment path as well as a
parameter change. Consumers collecting `hmm/*.hmm.gz` silently carry a discarded
family forward, contradicting the exit-0/exit-3 coherent-output guarantee.
Changing from refine to recruit-only can likewise leave obsolete seed/RF artifacts.

**Suggested correction:** after refusing foreign families and checking input-path
collisions, remove the exact previous artifact paths belonging to accepted input
names, or reconcile the complete artifact set for each final outcome. Cleanup
failure must prevent claiming coherent output. This does not require reintroducing
an ownership manifest or allowing a smaller input set over unrelated output.

**Missing regression:** successful run followed by a legitimate discard and by a
contained failure; assert that disk artifacts correspond exactly to successful
names. Existing retry tests only repeat successful outcomes or recover to success.

## F2 — P1: A single-file input library can be overwritten by its own output

**DONE.** Input/output validation now checks artifact and aggregate destinations against
model, FASTA and supplied-index inputs, including symlink and hard-link aliases.
Regression coverage includes single-model and multi-model inputs and aggregate collisions.
The focused validation run passed all 10 selected tests.

**Introduced by PR #7.** Location:
[`update_families.py:197–204`](src/mgnifam/update_families.py#L197).

The input/output overlap guard runs only when `hmm_input.is_dir()`. A library passed
as a file may occupy an artifact path. All models are loaded before writing, so the
run can complete successfully while destroying its input.

**Reproduction:** copy `tests/fixtures/mgnifams_v2.hmm.lib.gz` to a temporary
`out/hmm/v2_1.hmm.gz`, then run:

```bash
uv run mgnifam update_families \
  -i out/hmm/v2_1.hmm.gz -f v2.fa --skip_refine --output_dir out
```

Here `v2.fa` is the decompressed v2 fixture, and `out` must be an expendable temporary
directory. The supplied filename passes the existing-artifact check because `v2_1`
is an input family. Verified result: the input's bytes change, and reopening it now
finds **only `v2_1` instead of the original 14 models**. Retrying with the same input
can no longer reproduce the run. Refine mode can also overwrite a single model
with a scientifically different model.

**Suggested correction:** before creating or truncating anything, compare input
files against actual destination paths for the loaded family names and aggregates.
Account for resolved symlink aliases, and existing-file identity where applicable.
Apply the same protection to single-file and directory inputs.

**Missing regression:** an overlapping single-model file and multi-model library,
including an alias; rejection must leave the input byte-identical.

## F3 — P2: Accepted sequence names can corrupt CSV output

**Pre-existing; shared by both commands where applicable.** Locations:
[`generate_families.py:933–935`](src/mgnifam/generate_families.py#L933) and
[`generate_families.py:1022–1025`](src/mgnifam/generate_families.py#L1022).

Discard rows interpolate the cluster representative without CSV escaping. Metadata
wraps the protein name in quotes but does not double embedded quotes. FASTA names
and cluster representatives are not restricted to the safe alphabet used for HMM
family names, so the update command's family-name validation does not solve the
protein-metadata case.

**Verified through `emit_family()`:**

- A discarded cluster named `protein,version` produces
  `protein,version,too few sequences before initial hmmbuild,1`.
  `csv.reader` sees **four fields under a three-column header**.
- A successful representative named `protein"quote` is serialized as
  `"protein"quote"`. `csv.reader` recovers **`proteinquote"`**, changing its identity.
  The HMM and alignment writes succeed; this is silent output corruption.

**Suggested correction:** use the standard-library CSV writer for these rows,
with an explicit newline policy preserving deterministic bytes. If exact existing
quoting is a contract, preserve that convention while escaping embedded quotes
correctly. No dependency is needed.

**Missing regression:** round-trip comma and quote-bearing identifiers through
`csv.DictReader`, checking both the column count and exact recovered identity.

## F4 — P2: Slash-bearing FASTA identifiers fail after recruitment

**Pre-existing, now also reachable through `update_families`.** Locations:
[`generate_families.py:578–583`](src/mgnifam/generate_families.py#L578) and the related
`extract_first_part()` / `mask_sequence()` name encoding.

The pipeline uses `/` to append envelope coordinates, but does not distinguish that
suffix from a literal slash in a FASTA identifier. `parse_protein_name()` splits at
the first slash and fetches the truncated record name. Both validators accept these
inputs; the documented FASTA restrictions do not reserve `/`.

**Reproduced:** append `/v1` to each identifier in the uncompressed v2 FASTA, preserving
all residues and descriptions, and update the original 14-model library with
`--skip_refine`. Searching succeeds, but **all 14 families fail at artifact writing
and the run exits 3**. The first traceback ends with
`KeyError: '3387826881_356_473'`; the actual record is `3387826881_356_473/v1`.

This can waste every database search before reporting an input-format limitation.
It also means unrelated identifiers sharing a prefix before `/` collapse in the
membership/convergence sets.

**Suggested correction:** retain raw database identity separately from envelope
coordinates through the shared helpers. If slash-bearing identifiers are deliberately
unsupported, state and validate that restriction before expensive processing;
do not classify ordinary input names as internal family crashes.

**Missing regression:** a real indexed FASTA containing a slash-bearing name,
including two names with the same prefix, through recruitment and emission.

## F5 — P3: AGENTS instructions contradict the implemented scientific contract

**DONE.** AGENTS now describes final-round seed preservation and Stockholm IDs as
intentional fixes, names their existing regression tests, and removes the stale
hard-coded test count. Historical PLAN files and scientific code are unchanged.

**Pre-existing documentation drift, still present in the PR.** Location:
[`AGENTS.md:59–68`](AGENTS.md#L59).

The mandatory “Do not fix” section says the final-round model mismatch is preserved
and seed Stockholm files must not contain `#=GF ID`. Current code and CHANGELOG say
both were fixed: `Family.advance()` skips the unused final trim, `renumber_msa()`
sets the alignment name, and tests explicitly require the Stockholm ID. The
referenced “Preserved deliberately” changelog section no longer exists.

A contributor following these instructions could undo fixes or remove valid tests.
Update AGENTS to match the current guarantees, keeping the historical PLAN files
unchanged. Its “37 tests” count is also stale; current collection is 70.

## F6 — P3: The derived HMM library was not fully regenerated after domain reordering

**PR fixture drift.** Location:
[`tests/fixtures/mgnifams_v2.hmm.lib.gz`](tests/fixtures/mgnifams_v2.hmm.lib.gz), with
the regeneration contract in [`AGENTS.md`](AGENTS.md#L101).

Regenerated all 14 models from `mgnifams_v2.tsv` with default thresholds and chunk
`v2`, then compared each model's serialized bytes against the committed library.
Two differ:

| Family | Library CKSUM | Freshly generated CKSUM |
|---|---:|---:|
| `v2_3` | 4088430282 | 439396970 |
| `v2_14` | 2548862710 | 854247646 |

**Only the CKSUM lines differ**; the model parameters and calibration lines match.
This is a provenance/reproducibility correction, not evidence of changed recruitment.
Regenerate the derived library using the documented command after the score-ordering
fix. The automated tests currently do not consume this library.

## Checks that did not reveal additional defects

- The shared search explicitly uses `parallel="queries"`; extraction reads reported
  hits/domains and copies plain tuples. Per-hit score sorting implements issue #8's
  requested representative rule. The SSI mismatch check survives optimized Python.
- Cluster-order emission, rank assignment, discard containment, rollback-failure
  escalation, and shared-append fatal errors are structurally sound on fresh output
  directories. F1 concerns surviving files from an earlier run.
- HMM family names are validated and sorted consistently; round-one adoption is
  cleared where it is consumed. Recruit-only output uses the loaded model.
- The documented v3 fixture produced six members for `v2_10`, sixteen for `v2_6`,
  and no recruitment of decoy `9000000005`. Refinement changed their model lengths
  from 180 to 179 and 125 to 124 in two rounds. `v2_4` discarded on representative
  length, with empty retention/full-MSA-size fields as documented.
- Refine-mode scientific artifacts were byte-identical between CPU 1 with default
  batching/streaming and CPU 4 with batch size 16 and prefetch. Final HMM lengths
  agreed with successful delta rows on the v2 probe.

The shared algorithm is a useful design boundary: search semantics and family
transitions remain centralized while the commands own their input identities.
The principal design gap is output lifecycle ownership across retries. Addressing
that at the existing preparation/emission boundary is sufficient; a broad pipeline
extraction is not needed for these findings.
