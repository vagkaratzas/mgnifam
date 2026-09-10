#!/usr/bin/env python3
"""Refresh existing protein families against a new sequence database.

`generate_families` derives a family from an MMseqs2 cluster. This command starts from
families that already exist as HMMs and searches them against a new release, so a family
recruits from the larger database and picks up every fix in the shared code path.

Re-deriving the same families by reconstructing a cluster TSV from `<chunk>_families.tsv`
does not work, and fails quietly rather than loudly: those members are recorded in
parent-protein coordinates (`<protein>/<start>-<end>`) while the FASTA holds fragments
(`<protein>_<start>_<end>`), `run_initial_msa` skips unresolvable members silently, and
`check_seed_membership` divides by members that no longer exist in the new release. The
stored HMMs are the state that route spends its first round reconstructing.

Three things shape this module, and each was a decision rather than an accident:

1. **A family's identity is its model's `NAME`, preserved verbatim.** Nothing is
   renumbered, so a family keeps the same name across releases. `--chunk_num` therefore
   labels only the per-chunk aggregate files.

2. **Round 1 searches with the loaded model, not one built from a seed.** A seed cannot be
   recovered from an HMM, and refine mode does not need one: round 1 can never converge
   (`Family.total_checked_sequences` starts empty), so `advance` always produces the seed
   the later rounds and the final hand build require.

3. **The starting membership is round 1's own recruits.** There is no cluster to score a
   refined family against, so the trusted old model's first-round recruitment becomes the
   yardstick. Both sides of the comparison are then fragment names from the new database,
   which is what makes `check_seed_membership` correct here. Under `--skip_refine` there
   are no later rounds, so retention is 1.0 by construction and the check is a no-op.
"""

import argparse
import contextlib
import itertools
import time
from collections.abc import Sequence as SequenceCollection
from dataclasses import dataclass
from pathlib import Path

import pyhmmer

from mgnifam.generate_families import (
    ALPHABET,
    CHUNK_PATTERN,
    DISCARDED_HEADER,
    EXIT_CRASHED_FAMILIES,
    FAMILY_DIRECTORIES,
    INTERNAL_ERROR_PREFIX,
    MAX_ROUNDS,
    METADATA_HEADER,
    ChunkCorrupted,
    Family,
    FamilyState,
    IndexedSequences,
    Writers,
    configure_logger,
    deterministic_gzip_text,
    emit_family,
    extract_records,
    family_guard,
    filter_hits,
    is_gzipped,
    resolve_index,
    run_hmmbuild,
    search,
    unmask_sequence_names,
)

DELTA_HEADER = (
    "family_id,model_length_before,model_length_after,round1_recruits,"
    "full_msa_size,retention,rounds_run,converged,outcome\n"
)
NO_HITS = "no hits in the new database"
# A family name is interpolated into artifact paths and into unquoted CSV fields, and it
# arrives from a third-party file rather than from this program. `CHUNK_PATTERN`'s alphabet
# is already the trusted one for the other half of every artifact filename, and it admits
# no "/", ",", whitespace or control character -- so one rule closes both path traversal
# and CSV-column corruption. Only "." and ".." get through it, and they are rejected by name.
FAMILY_NAME_PATTERN = CHUNK_PATTERN
TRAVERSAL_NAMES = frozenset({".", ".."})
# Every artifact this command can write, by directory. Used to clear a previous run's
# output exactly, without a glob that could match a name it does not own.
ARTIFACT_SUFFIXES = {"rf": ".txt", "hmm": ".hmm.gz", "seed_msa": ".sto.gz", "full_msa": ".sto.gz"}


@dataclass
class FamilyDelta:
    """The scalars of one family's update, captured as they are computed.

    Not read back off the `Family` afterwards: `discard` clears the records, model and
    alignments to release memory, and `finish` computes the membership fraction into a
    local before it may discard on representative length two statements later. A family
    rejected for length has computed its retention and then destroyed the evidence.
    """

    family_id: str
    model_length_before: int
    round1_recruits: int | None = None
    rounds_run: int = 0


def delta_row(delta: FamilyDelta, family: Family) -> str:
    """Format one delta row, once the family's outcome is final.

    Every field but the id, the starting model length and the outcome is nullable and
    written as an empty CSV field: a family discarded early never reached the stage that
    would have produced one. `retention` in particular is empty when round 1 recruited
    nothing, because there is then no yardstick at all -- `check_seed_membership` would
    divide by an empty set rather than score zero.

    `model_length_after` is the length of the model that recruited the final membership,
    i.e. the last round's query. Under `--skip_refine` that is the loaded model itself.
    """
    outcome = family.discard_reason if family.state is FamilyState.DISCARDED else "successful"
    fields: tuple[object | None, ...] = (
        delta.family_id,
        delta.model_length_before,
        family.qlen or None,
        delta.round1_recruits,
        family.full_msa_num_seqs or None,
        family.membership,
        delta.rounds_run or None,
        family.ever_converged,
        outcome,
    )
    return ",".join("" if field is None else str(field) for field in fields) + "\n"


def load_hmms(path: Path) -> list[tuple[str, pyhmmer.plan7.HMM]]:
    """Load every model under `path`, a directory of HMM files or a single HMM library.

    Returned sorted by family name, for both input forms alike. Filename order and library
    order are different orders over the same models -- `a.hmm` holding `NAME B` beside
    `z.hmm` holding `NAME A` loads as B, A, while a library holding them in A, B order
    loads as A, B, and the aggregate files would then differ byte-for-byte for the same
    input. Sorting on the identity key is what makes the two forms interchangeable, and it
    also stops renaming a file from reordering output when identity comes from `NAME`.

    `hmm.name` is used directly. It is a `str` under the pinned pyhmmer, and there is no
    filename fallback because a model without a `NAME` cannot be parsed at all.
    """
    if path.is_dir():
        # Sorted only for a deterministic parse order; the result is re-sorted by name.
        files = sorted(p for p in path.iterdir() if p.name.endswith((".hmm", ".hmm.gz")))
    else:
        files = [path]

    models: list[tuple[str, pyhmmer.plan7.HMM]] = []
    for file in files:
        with pyhmmer.plan7.HMMFile(file) as handle:
            models.extend((hmm.name, hmm) for hmm in handle)
    return sorted(models, key=lambda pair: pair[0])


def validate_inputs(options: argparse.Namespace) -> list[tuple[str, pyhmmer.plan7.HMM]]:
    """Check every input and return the loaded models, before any output exists.

    Ordering matters as much as the checks: nothing here may create a directory, so a
    rejected run leaves the output root exactly as it found it.
    """
    if CHUNK_PATTERN.fullmatch(options.chunk_num) is None:
        raise ValueError("chunk_num must match [A-Za-z0-9._-]+")
    fasta = Path(options.fasta_file)
    if not fasta.is_file():
        raise ValueError("fasta_file must be an existing file")
    if is_gzipped(fasta):
        raise ValueError("fasta_file must be uncompressed; SSI cannot seek in gzip streams")
    hmm_input = Path(options.hmm_input)
    if not hmm_input.exists():
        raise ValueError("hmm_input must be an existing file or directory")
    if options.fasta_index and not Path(options.fasta_index).is_file():
        raise ValueError("fasta_index must be an existing file; a supplied index is never built")
    if options.cpus < 1:
        raise ValueError("cpus must be at least 1")
    if options.batch_size < 0:
        raise ValueError("batch_size must be non-negative (0 selects 2 * cpus)")
    if options.batch_size and options.batch_size < options.cpus:
        raise ValueError("batch_size must be at least cpus")
    for name in (
        "discard_min_starting_membership",
        "max_seq_identity",
        "max_gap_occupancy",
        "recruit_hit_length_percentage",
    ):
        if not 0 <= getattr(options, name) <= 1:
            raise ValueError(f"{name} must be in [0, 1]")
    if options.recruit_evalue_cutoff <= 0:
        raise ValueError("recruit_evalue_cutoff must be positive")
    if options.max_seed_seqs < 1:
        raise ValueError("max_seed_seqs must be at least 1")
    if options.discard_min_rep_length < 1 or options.discard_max_rep_length < 1:
        raise ValueError("representative length bounds must be at least 1")
    if options.discard_min_rep_length > options.discard_max_rep_length:
        raise ValueError("minimum representative length cannot exceed maximum")

    # A directory of models is also where this command would write its own `hmm/` output.
    # The models are all in memory by then, so overwriting them would not break the run --
    # it would silently destroy the input the user would need to re-run it.
    if hmm_input.is_dir():
        root = options.output_dir.resolve()
        forbidden = {root} | {(root / directory).resolve() for directory in FAMILY_DIRECTORIES}
        if hmm_input.resolve() in forbidden:
            raise ValueError("hmm_input must not be the output_dir or one of its artifact folders")

    models = load_hmms(hmm_input)
    if not models:
        raise ValueError("hmm_input contains no HMMs")
    seen: set[str] = set()
    for name, hmm in models:
        if FAMILY_NAME_PATTERN.fullmatch(name) is None or name in TRAVERSAL_NAMES:
            raise ValueError(
                f"HMM NAME {name!r} is not a usable family name; it must match "
                "[A-Za-z0-9._-]+ and be neither '.' nor '..'"
            )
        if name in seen:
            raise ValueError(f"duplicate HMM NAME {name!r}; family names must be unique")
        seen.add(name)
        if hmm.alphabet != ALPHABET:
            raise ValueError(f"HMM {name!r} is not an amino-acid model")
    validate_output_paths(options, seen)
    return models


def validate_output_paths(options: argparse.Namespace, names: set[str]) -> None:
    """Reject destinations that alias an input, before output creation or retry cleanup."""
    source = Path(options.hmm_input)
    inputs = (
        [p for p in source.iterdir() if p.name.endswith((".hmm", ".hmm.gz"))]
        if source.is_dir()
        else [source]
    )
    inputs.append(Path(options.fasta_file))
    if options.fasta_index:
        inputs.append(Path(options.fasta_index))
    resolved_inputs = {p.resolve() for p in inputs}
    input_inodes = {(p.stat().st_dev, p.stat().st_ino) for p in inputs}
    root = options.output_dir
    destinations = [
        root / directory / f"{name}{suffix}"
        for directory, suffix in ARTIFACT_SUFFIXES.items()
        for name in sorted(names)
    ]
    prefix = f"{options.chunk_num}_updated"
    destinations.extend(
        root / f"{prefix}_{suffix}"
        for suffix in (
            "families.tsv",
            "metadata.csv",
            "discarded.csv",
            "successful.txt",
            "converged.txt",
            "reps.fasta.gz",
            "delta.csv",
        )
    )
    destinations.append(root / f"{prefix}.log")
    if not options.fasta_index:
        destinations.append(root / f"{Path(options.fasta_file).name}.ssi")
    for path in destinations:
        if path.resolve() in resolved_inputs or (
            path.exists() and (path.stat().st_dev, path.stat().st_ino) in input_inodes
        ):
            raise ValueError(f"output path {path} overlaps an input file")


def prepare_output_directories(root: Path, names: SequenceCollection[str]) -> None:
    """Refuse foreign families, then clear the input families' previous artifacts.

    Call only after input/output collision validation. The input names identify exactly
    which files a retry may replace. All four artifact types must be cleared even when
    this run skips refinement: a discarded family must leave no old model behind, and
    recruit-only output must not retain a previous refinement's seed or RF annotation.

    A smaller input set over a larger output set remains an error: nothing identifies
    the omitted families as ours to remove. Check the entire tree before deleting any
    accepted paths. A cleanup error propagates and aborts before aggregates are opened.
    """
    owned = set(names)
    strays = sorted(
        str(path.relative_to(root))
        for directory, suffix in ARTIFACT_SUFFIXES.items()
        for path in (root / directory).glob(f"*{suffix}")
        if path.name.removesuffix(suffix) not in owned
    )
    if strays:
        listed = ", ".join(strays[:5]) + (", ..." if len(strays) > 5 else "")
        raise ValueError(
            f"{root} already holds artifacts for families this run does not update "
            f"({listed}). Use a fresh output_dir, or remove them first."
        )
    for directory in FAMILY_DIRECTORIES:
        (root / directory).mkdir(parents=True, exist_ok=True)
    for directory, suffix in ARTIFACT_SUFFIXES.items():
        for name in names:
            (root / directory / f"{name}{suffix}").unlink(missing_ok=True)


def parse_args(args: SequenceCollection[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="mgnifam update_families",
        description="Refresh existing HMM families against a newer sequence database.",
    )
    parser.add_argument("-i", "--hmm_input", required=True, help="HMM directory or HMM library")
    parser.add_argument("-f", "--fasta_file", required=True)
    parser.add_argument(
        "--skip_refine",
        action="store_true",
        help="recruit once and align; do not rebuild the seed MSA or the model. "
        "Makes --max_seq_identity, --max_seed_seqs and --max_gap_occupancy inert.",
    )
    parser.add_argument("-p", "--cpus", type=int, default=8)
    parser.add_argument(
        "-n", "--chunk_num", default="1", help="prefix for the per-chunk aggregate files only"
    )
    parser.add_argument("--discard_min_rep_length", type=int, default=75)
    parser.add_argument("--discard_max_rep_length", type=int, default=2000)
    parser.add_argument("--discard_min_starting_membership", type=float, default=0.9)
    parser.add_argument("--max_seq_identity", type=float, default=0.8)
    parser.add_argument("--max_seed_seqs", type=int, default=2000)
    parser.add_argument("--max_gap_occupancy", type=float, default=0.5)
    parser.add_argument("--recruit_evalue_cutoff", type=float, default=0.001)
    parser.add_argument("--recruit_hit_length_percentage", type=float, default=0.9)
    parser.add_argument("--fasta_index")
    parser.add_argument("--output_dir", type=Path, default=Path("output"))
    parser.add_argument("--batch_size", type=int, default=0)
    parser.add_argument("--prefetch_targets", action="store_true")
    return parser.parse_args(args)


def open_writers(
    stack: contextlib.ExitStack, root: Path, chunk: str, indexed: IndexedSequences
) -> Writers:
    prefix = f"{chunk}_updated"

    def text(suffix: str) -> object:
        return stack.enter_context((root / f"{prefix}_{suffix}").open("w"))

    writers = Writers(
        root=root,
        indexed=indexed,
        refined_families=text("families.tsv"),  # type: ignore[arg-type]
        discarded_clusters=text("discarded.csv"),  # type: ignore[arg-type]
        successful_clusters=text("successful.txt"),  # type: ignore[arg-type]
        converged_families=text("converged.txt"),  # type: ignore[arg-type]
        family_metadata=text("metadata.csv"),  # type: ignore[arg-type]
        family_representatives=stack.enter_context(
            deterministic_gzip_text(root / f"{prefix}_reps.fasta.gz")
        ),
        family_delta=text("delta.csv"),  # type: ignore[arg-type]
    )
    # Written before any result, so a chunk that produces no families still parses.
    writers.family_metadata.write(METADATA_HEADER)
    writers.discarded_clusters.write(DISCARDED_HEADER)
    assert writers.family_delta is not None
    writers.family_delta.write(DELTA_HEADER)
    return writers


def recruit_only(
    family: Family, delta: FamilyDelta, options: argparse.Namespace, indexed: IndexedSequences
) -> None:
    """Consume round 1's hits without refining, then hand over to the exit branch.

    This is `advance`'s opening -- the same filter, the same empty-recruitment discard --
    stopping where `advance` would start rebuilding a seed for a round that will not
    happen. `finish` is then reused unchanged: it re-filters the cached records with the
    envelope requirement waived, exactly as it does for a converged family.
    """
    filtered_sequences = filter_hits(
        family.records,
        family.qlen,
        False,
        options.recruit_hit_length_percentage,
        indexed,
    )
    if not filtered_sequences:
        family.discard("low complexity model - confounding cluster", 0.0)
        return
    family.members = unmask_sequence_names(filtered_sequences)
    delta.round1_recruits = len(filtered_sequences)


def main(args: SequenceCollection[str] | None = None) -> None:
    options = parse_args(args)
    models = validate_inputs(options)
    if options.batch_size == 0:
        options.batch_size = 2 * options.cpus

    root = options.output_dir
    prepare_output_directories(root, [name for name, _ in models])
    index_path = resolve_index(options, root)
    logger = configure_logger(root / f"{options.chunk_num}_updated.log")

    started = time.monotonic()
    total_batches = -(-len(models) // options.batch_size)
    max_round = 1 if options.skip_refine else MAX_ROUNDS
    logger.info(
        "start families=%d batches=%d batch_size=%d cpus=%d skip_refine=%s prefetch=%s index=%s",
        len(models),
        total_batches,
        options.batch_size,
        options.cpus,
        options.skip_refine,
        options.prefetch_targets,
        index_path,
    )

    try:
        with contextlib.ExitStack() as stack:
            index_reader = stack.enter_context(pyhmmer.easel.SSIReader(index_path))
            indexed_file = stack.enter_context(
                pyhmmer.easel.SequenceFile(options.fasta_file, digital=False, index=index_reader)
            )
            indexed_sequences = IndexedSequences(indexed_file)
            target_file = stack.enter_context(
                pyhmmer.easel.SequenceFile(options.fasta_file, digital=True, alphabet=ALPHABET)
            )
            targets = target_file.read_block() if options.prefetch_targets else target_file
            writers = open_writers(stack, root, options.chunk_num, indexed_sequences)

            success_count = 0
            processed = 0
            crashed = 0
            for batch_number, batch in enumerate(
                itertools.batched(models, options.batch_size, strict=False), 1
            ):
                # `members` starts empty: an updated family has no cluster behind it, so
                # round 1's own recruits become the yardstick `finish` scores it against.
                # Only round 1 -- `Family.advance` clears the flag in the same block that
                # reads it (`generate_families.py`, `if self.adopt_recruits_as_members:`),
                # so rounds 2 and 3 leave `members` alone. Under --skip_refine `advance` is
                # never called at all and `recruit_only` sets `members` directly, which is
                # why the flag is off there.
                active = [
                    Family(
                        representative=name,
                        members=[],
                        adopt_recruits_as_members=not options.skip_refine,
                    )
                    for name, _ in batch
                ]
                deltas = {
                    name: FamilyDelta(family_id=name, model_length_before=hmm.M)
                    for name, hmm in batch
                }

                for round_number in range(1, max_round + 1):
                    if round_number == 1:
                        # The loaded models are the round-1 queries. No seed exists yet, and
                        # for a recruit-only update none ever will.
                        pending = list(zip(active, [hmm for _, hmm in batch], strict=True))
                    else:
                        running = [f for f in active if f.state is FamilyState.RUNNING]
                        if not running:
                            break
                        pending = []
                        for index, family in enumerate(running):
                            with family_guard(family, logger, f"round {round_number} model build"):
                                pending.append(
                                    (
                                        family,
                                        run_hmmbuild(
                                            family.seed_msa,  # type: ignore[arg-type]
                                            f"pending_{batch_number}_{round_number}_{index}",
                                        ),
                                    )
                                )
                    if not pending:
                        continue
                    searching = [family for family, _ in pending]
                    hmms = [hmm for _, hmm in pending]
                    with search(
                        hmms,
                        targets,
                        cpus=options.cpus,
                        evalue=options.recruit_evalue_cutoff,
                        logger=logger,
                        batch_number=batch_number,
                        round_number=round_number,
                    ) as searched:
                        for family, hmm, hits in zip(searching, hmms, searched, strict=True):
                            with family_guard(family, logger, f"round {round_number}"):
                                delta = deltas[family.representative]
                                family.hmm = hmm
                                family.qlen = hits.query.M
                                family.records = extract_records(hits)
                                delta.rounds_run = round_number
                                if not family.records:
                                    # Distinct from "low complexity": that one means the
                                    # model found hits and none were long enough. This one
                                    # means the model found nothing at all, which for an
                                    # update run is the answer the report exists to carry.
                                    family.discard(NO_HITS, 0.0)
                                    continue
                                if options.skip_refine:
                                    recruit_only(family, delta, options, indexed_sequences)
                                else:
                                    family.advance(options, indexed_sequences, round_number)
                                    if round_number == 1:
                                        delta.round1_recruits = len(family.members) or None
                    logger.info(
                        "batch=%d/%d round=%d searched=%d still_running=%d converged=%d "
                        "discarded=%d elapsed=%.1fs",
                        batch_number,
                        total_batches,
                        round_number,
                        len(searching),
                        sum(1 for f in active if f.state is FamilyState.RUNNING),
                        sum(1 for f in active if f.state is FamilyState.CONVERGED),
                        sum(1 for f in active if f.state is FamilyState.DISCARDED),
                        time.monotonic() - started,
                    )

                for family in active:
                    with family_guard(family, logger, "the exit branch"):
                        family.finish(options, indexed_sequences)
                for family, (_, loaded_hmm) in zip(active, batch, strict=True):
                    name = family.representative
                    # Under --skip_refine the model is unchanged, so the one that was
                    # searched with is written out; re-serialising it rather than copying
                    # the input file is what drops a legacy DATE/COM stamp.
                    final_hmm = loaded_hmm if options.skip_refine else None
                    emitted = False
                    with family_guard(family, logger, "artifact writing"):
                        success_count = emit_family(
                            family,
                            success_count,
                            options.chunk_num,
                            writers,
                            family_name=name,
                            family_id=name,
                            final_hmm=final_hmm,
                            delta_row=delta_row(deltas[name], family),
                        )
                        emitted = True
                    if not emitted:
                        # The guard turned the failure into a discard, and a discard is a
                        # result. The row is rebuilt because the outcome has changed.
                        emit_family(
                            family,
                            success_count,
                            options.chunk_num,
                            writers,
                            family_name=name,
                            family_id=name,
                            delta_row=delta_row(deltas[name], family),
                        )
                processed += len(active)
                crashed += sum(
                    1
                    for family in active
                    if family.discard_reason.startswith(INTERNAL_ERROR_PREFIX)
                )
                logger.info(
                    "batch=%d/%d written families=%d/%d successful=%d discarded=%d elapsed=%.1fs",
                    batch_number,
                    total_batches,
                    processed,
                    len(models),
                    success_count,
                    processed - success_count,
                    time.monotonic() - started,
                )
        logger.info(
            "DONE. families=%d successful=%d discarded=%d crashed=%d elapsed=%.1fs",
            processed,
            success_count,
            processed - success_count,
            crashed,
            time.monotonic() - started,
        )
        if crashed:
            raise SystemExit(EXIT_CRASHED_FAMILIES)
    except ChunkCorrupted:
        logger.exception("chunk output is corrupted and must not be consumed")
        raise SystemExit(1) from None
    finally:
        for handler in logger.handlers:
            handler.close()
        logger.handlers.clear()


if __name__ == "__main__":
    main()
