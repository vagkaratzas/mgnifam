#!/usr/bin/env python3
"""Iterative HMM-based protein family generation over very large sequence databases.

Each cluster becomes a candidate family: build an HMM from its seed alignment, search
the whole database for new members, re-align, repeat. A family converges when a round
recruits nothing new, or is forced out after three rounds; either way it exits through
a final hand-architecture build. Clusters that fail a length, membership or
sequence-count check are discarded.

Three invariants shape the design, and breaking any of them reintroduces a bug:

1. The database never enters memory. Targets stream from a `SequenceFile`; individual
   sequences are fetched through an Easel SSI index. `--prefetch_targets` opts back
   into an in-RAM block, which is a pure speed/memory trade with identical results.

2. `hmmsearch` is always called with `parallel="queries"`. Left to choose, pyhmmer
   selects `parallel="targets"` whenever the query count is below `--cpus`, and its
   merge step retains hits that no single pipeline would have reported. Results would
   then depend on the core count of the machine.

3. Families are searched in round-major waves rather than one at a time. This is what
   makes (1) affordable: `hmmsearch` uses `min(cpus, n_queries)` workers, so a single
   query per call would leave every core but one idle.

4. A family's seed MSA, its HMM and its full MSA always describe the same round. The
   redundancy trim that produces the next seed exists only to feed the next
   `hmmbuild`, so on either exit -- convergence or `MAX_ROUNDS` -- it does not run.
   See `Family.advance`.

Family ids are the 1-based rank among *successful* families in cluster-file order, so
nothing may be written until a family's fate is known. See `emit_family`.
"""

import argparse
import contextlib
import gzip
import io
import itertools
import logging
import os
import re
import shutil
import tempfile
import threading
import time
from collections.abc import Iterable, Iterator, Mapping
from collections.abc import Sequence as SequenceCollection
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import IO, Any, NamedTuple, cast

import numpy as np
import pyfamsa
import pyhmmer
import pytrimal

ALPHABET = pyhmmer.easel.Alphabet.amino()
MAX_ROUNDS = 3
CHUNK_PATTERN = re.compile(r"[A-Za-z0-9._-]+")
# Per-family artifacts: one file per family, so they get a directory each. Everything
# else is a single file per chunk and lives flat in the output root as `<chunk>_*`.
FAMILY_DIRECTORIES = ("seed_msa", "full_msa", "hmm", "rf")


class DuplicateSequenceName(ValueError):
    """Raised when an SSI index cannot be built because FASTA names repeat."""


class IndexMismatchError(RuntimeError):
    """Raised when an SSI index returns a record under the wrong key."""


class Sequence(NamedTuple):
    id: str
    seq: str


Record = tuple[str, int, int, int]


class SizedIterator:
    def __init__(self, iterator: Iterable[str], length: int) -> None:
        self.iterator = iter(iterator)
        self.length = length

    def __len__(self) -> int:
        return self.length

    def __iter__(self) -> "SizedIterator":
        return self

    def __next__(self) -> str:
        item = next(self.iterator)
        self.length -= 1
        return item


def pyfamsa_to_pyhmmer(alignment: pyfamsa.Alignment) -> pyhmmer.easel.DigitalMSA:
    return pyhmmer.easel.TextMSA(
        sequences=[
            pyhmmer.easel.TextSequence(name=seq.id.decode(), sequence=seq.sequence.decode())
            for seq in alignment
        ]
    ).digitize(ALPHABET)


def pytrimal_to_pyhmmer(alignment: pytrimal.Alignment) -> pyhmmer.easel.TextMSA:
    return pyhmmer.easel.TextMSA(
        sequences=[
            pyhmmer.easel.TextSequence(
                name=name.decode() if isinstance(name, bytes) else name,
                sequence=sequence,
            )
            for name, sequence in zip(alignment.names, alignment.sequences, strict=True)
        ]
    )


def build_ssi_index(fasta: str | os.PathLike[str], ssi_path: str | os.PathLike[str]) -> None:
    """Build an Easel SSI index for `fasta`, replacing `ssi_path` atomically.

    `fasta` must be uncompressed: Easel cannot seek within a gzip stream. Its sequence
    names must be unique, or `DuplicateSequenceName` is raised.

    The index is built inside a private directory on the destination filesystem and
    moved into place only on success, so a concurrent reader never observes a partial
    index and a crashed build leaves nothing behind. That directory is not merely
    tidiness: `SSIWriter.close()` raises *before* Easel closes its own writer, so on a
    duplicate key the native handle and any external-sort scratch files are still open
    and cannot be cleaned up by unlinking a known path.

    Safe against PID collisions and cooperating concurrent writers. Not a defence
    against an actor with write access to the destination directory.
    """
    fasta_path = Path(fasta)
    destination = Path(ssi_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_directory = Path(tempfile.mkdtemp(dir=destination.parent))
    temporary_index = temporary_directory / "index.ssi"

    try:
        writer = pyhmmer.easel.SSIWriter(temporary_index)
        # add_file wants an integer format code, not the format name.
        file_number = writer.add_file(fasta_path.name, pyhmmer.easel.SequenceFile._FORMATS["fasta"])
        # Easel positions by record_offset alone, so data_offset and record_length are
        # left at 0. That is sufficient for whole-record fetches and rules out
        # subsequence fetches, which this module never performs.
        #
        # Duplicate detection is left to SSIWriter: a Python set of every key would be
        # database-sized state, which is precisely what this module exists to avoid.
        with fasta_path.open("rb") as fasta_handle:
            while True:
                record_offset = fasta_handle.tell()
                line = fasta_handle.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    name = line[1:].split(maxsplit=1)[0].decode()
                    writer.add_key(name, file_number, record_offset)
        try:
            writer.close()
        except ValueError as error:
            if "duplicate" in str(error).lower():
                raise DuplicateSequenceName(
                    "FASTA sequence names must be unique to build an SSI index"
                ) from error
            raise
        os.replace(temporary_index, destination)
    finally:
        shutil.rmtree(temporary_directory, ignore_errors=True)


def fetch_indexed_sequence(
    indexed: Mapping[str, pyhmmer.easel.TextSequence], name: str
) -> pyhmmer.easel.TextSequence:
    """Fetch `name` from an SSI-indexed file, proving the index belongs to the file.

    Raises KeyError if absent, IndexMismatchError if the index returns a record under
    the wrong key -- which is how a stale index built against a different FASTA of the
    same basename shows up. Deliberately not an `assert`: `python -O` strips those, and
    the failure it guards against is a silent wrong-sequence read.
    """
    fetched = indexed[name]
    if fetched.name != name:
        raise IndexMismatchError(
            f"SSI index mismatch: requested {name!r}, fetched {fetched.name!r}"
        )
    return fetched


class IndexedSequences:
    def __init__(self, sequence_file: pyhmmer.easel.SequenceFile) -> None:
        self.indexed = sequence_file.indexed

    def get(self, name: str, *, missing_ok: bool = False) -> Sequence | None:
        try:
            fetched = fetch_indexed_sequence(self.indexed, name)
        except KeyError:
            if missing_ok:
                return None
            raise
        return Sequence(name, fetched.sequence)


def run_initial_msa(
    members: Iterable[str], indexed_sequences: IndexedSequences, cpus: int
) -> pyhmmer.easel.DigitalMSA | None:
    """Align a cluster's members into its first seed MSA, or None if none were found.

    Members absent from the FASTA are skipped rather than raising, matching legacy.
    Everywhere else a missing name is an error, because it can only come from a hit
    the database itself produced.
    """
    sequences = []
    for member in members:
        sequence = indexed_sequences.get(member, missing_ok=True)
        if sequence is not None:
            sequences.append(pyfamsa.Sequence(id=member.encode(), sequence=sequence.seq.encode()))
    if not sequences:
        return None
    # FAMSA's output is independent of the thread count, so honouring `cpus` costs
    # nothing in reproducibility. Legacy used every core regardless of the allocation.
    alignment = pyfamsa.Aligner(threads=cpus).align(sequences)
    return pyfamsa_to_pyhmmer(alignment)


def run_hmmbuild(
    seed_msa: pyhmmer.easel.DigitalMSA, name: str, *, hand: bool = False
) -> pyhmmer.plan7.HMM:
    """Build an HMM from a seed alignment, naming it `name`.

    `hand=True` takes the match-state architecture from the alignment's RF line instead
    of inferring it, and is used only for the final model of a successful family.

    `Builder` requires a named MSA. The name does not affect the model, so intermediate
    rounds may pass a placeholder; only the name of the final model is written out.
    """
    seed_msa.name = name
    architecture = "hand" if hand else "fast"
    builder = pyhmmer.plan7.Builder(ALPHABET, architecture=architecture, seed=42)
    background = pyhmmer.plan7.Background(ALPHABET)
    hmm, _, _ = builder.build_msa(seed_msa, background)
    return hmm


@contextlib.contextmanager
def search(
    hmms: SequenceCollection[pyhmmer.plan7.HMM],
    targets: Any,
    *,
    cpus: int,
    evalue: float,
    logger: logging.Logger,
    batch_number: int,
    round_number: int,
) -> Iterator[Iterator[Any]]:
    """Search `hmms` against `targets`, yielding TopHits in query order.

    `targets` may be a streaming SequenceFile or a prefetched DigitalSequenceBlock; the
    results are identical either way. Must be used as a context manager, and the
    iterator must be consumed inside the `with` block.

    `parallel="queries"` is passed explicitly and must never be removed: pyhmmer would
    otherwise switch to target-parallelism whenever a wave holds fewer families than
    `cpus` -- which the final wave of a batch usually does -- and target-parallelism
    reports hits that a single pipeline would not.

    A heartbeat thread logs progress every 300s. It cannot be driven by the pyhmmer
    callback alone, which fires only once a whole query has finished; on a large
    database that is silence for hours. The iterator is closed on the way out because
    abandoning it part-way would otherwise leak the dispatcher's worker threads.
    """
    completed = 0
    lock = threading.Lock()
    stop_event = threading.Event()
    started = time.monotonic()

    def callback(_query: pyhmmer.plan7.HMM, _index: int) -> None:
        nonlocal completed
        with lock:
            completed += 1

    def heartbeat() -> None:
        while not stop_event.wait(300):
            with lock:
                current = completed
            logger.info(
                "heartbeat batch=%d round=%d families=%d completed=%d/%d elapsed=%.1fs",
                batch_number,
                round_number,
                len(hmms),
                current,
                len(hmms),
                time.monotonic() - started,
            )

    timer = threading.Thread(target=heartbeat, daemon=True)
    timer.start()
    try:
        iterator = pyhmmer.hmmer.hmmsearch(
            hmms,
            targets,
            cpus=cpus,
            callback=callback,
            parallel="queries",
            E=evalue,
            seed=42,
        )
        with contextlib.closing(iterator) as searched:
            yield searched
    finally:
        stop_event.set()
        timer.join()


def extract_records(top_hits: Any) -> list[Record]:
    """Flatten a TopHits into (name, target_length, env_from, env_to) tuples.

    Only hits and domains that cleared the pipeline's reporting thresholds are
    returned, in HMMER's own ranking order.

    Two constraints drive this:

    Iterating `top_hits` directly would also yield entries HMMER stored but did
    not report, i.e. sequences that failed `--recruit_evalue_cutoff`. The legacy
    script did exactly that and recruited them (see CHANGELOG 1.0.0). Worse, the
    number of such entries depends on the parallelisation strategy, which made
    recruitment vary with `--cpus`.

    The tuples must be plain values, not pyhmmer objects: a Domain holds a
    reference to its Hit, which holds the entire TopHits. Caching Domain objects
    across a batch would pin every result graph in memory.
    """
    return [
        (hit.name, hit.length, domain.env_from, domain.env_to)
        for hit in top_hits.reported
        for domain in hit.domains.reported
    ]


def mask_sequence(sequence: Sequence, env_from: int, env_to: int) -> Sequence:
    return Sequence(f"{sequence.id}/{env_from}_{env_to}", sequence.seq[env_from - 1 : env_to])


def filter_hits(
    records: Iterable[Record],
    qlen: int,
    exit_flag: bool,
    recruit_hit_length_percentage: float,
    indexed_sequences: IndexedSequences,
) -> list[Sequence]:
    """Resolve hit records to sequences, keeping those whose envelope is long enough.

    Sequences whose envelope covers only part of the target are masked down to it, and
    renamed `<name>/<env_from>_<env_to>`.

    `exit_flag` waives the length requirement. The exit branch calls this a second time
    over the same records, which is why the records are cached rather than re-searched.
    """
    filtered_sequences = []
    for name, target_length, env_from, env_to in records:
        envelope_length = env_to - env_from + 1
        if exit_flag or envelope_length >= recruit_hit_length_percentage * qlen:
            sequence = cast(Sequence, indexed_sequences.get(name))
            if envelope_length < target_length:
                sequence = mask_sequence(sequence, env_from, env_to)
            filtered_sequences.append(sequence)
    return filtered_sequences


def run_hmmalign(
    hmm: pyhmmer.plan7.HMM, family_sequences: Iterable[Sequence], cpus: int
) -> pyhmmer.easel.TextMSA:
    sequences = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=sequence.id, sequence=sequence.seq)
        for sequence in family_sequences
    ).digitize(ALPHABET)
    return cast(
        pyhmmer.easel.TextMSA,
        pyhmmer.hmmer.hmmalign(hmm, sequences, cpus=cpus, trim=False),
    )


def msa_stats(msa: pyhmmer.easel.TextMSA) -> tuple[int, int]:
    """Return (sequence count, ungapped length of the first row).

    Row 0 is the family representative: hits arrive in HMMER's ranking order, so the
    best-scoring sequence leads the alignment.
    """
    number_of_sequences = len(msa.names)
    non_gap_representative_length = (
        len(re.sub(r"[.\-~]", "", msa.alignment[0])) if number_of_sequences else 0
    )
    return number_of_sequences, non_gap_representative_length


def clip_env_ends(msa: pyhmmer.easel.TextMSA) -> pyhmmer.easel.TextMSA:
    """Trim the envelope overhangs from a seed alignment.

    Drops the leading and trailing runs of columns where the `#=GC RF` line holds "."
    -- the N- and C-terminal residues that fell outside the model's match states.
    Interior "." columns are inserts between match states and are kept.

    Applied to every seed-path alignment, never to the final full MSA: the full MSA is
    the record of what each member sequence contributed, overhangs included.
    """
    if msa.reference is None:
        raise ValueError("hmmalign result has no RF reference annotation")
    matching_columns = [index for index, value in enumerate(msa.reference) if value != "."]
    if not matching_columns:
        return msa
    return msa.select(columns=range(matching_columns[0], matching_columns[-1] + 1))


def run_pytrimal_reps(
    full_msa: pyhmmer.easel.TextMSA, threshold: float, max_seed_seqs: int
) -> pyhmmer.easel.TextMSA:
    if full_msa.reference is None:
        raise ValueError("alignment has no RF reference annotation")
    reference = np.array(list(full_msa.reference))
    sequences = SizedIterator(
        (sequence.upper().replace(".", "-").replace("~", "-") for sequence in full_msa.alignment),
        len(full_msa.alignment),
    )
    alignment = pytrimal.Alignment([name.encode() for name in full_msa.names], sequences)
    alignment = pytrimal.RepresentativeTrimmer(identity_threshold=threshold).trim(alignment)
    reference = reference[alignment.residues_mask]

    if len(list(alignment.names)) > max_seed_seqs:
        alignment = pytrimal.RepresentativeTrimmer(clusters=max_seed_seqs).trim(alignment)
        reference = reference[alignment.residues_mask]

    result = pytrimal_to_pyhmmer(alignment)
    result.reference = "".join(reference)
    return result


def calculate_trim_positions(
    sequence_matrix: np.ndarray, occupancy_threshold: float
) -> tuple[int, int] | None:
    """Return the first and last column indices whose non-gap occupancy exceeds the threshold.

    Both bounds are inclusive. Returns None when no column qualifies, which callers must
    distinguish from a span: an `np.argmax` over an all-False array silently yields 0, and
    the legacy script reported the full alignment in that case.
    """
    numeric_matrix = np.where(sequence_matrix == "-", 0, 1)
    column_percentages = np.sum(numeric_matrix, axis=0) / numeric_matrix.shape[0]
    passing_columns = np.flatnonzero(column_percentages > occupancy_threshold)
    if passing_columns.size == 0:
        return None
    return int(passing_columns[0]), int(passing_columns[-1])


def clip_ends(msa: pyhmmer.easel.TextMSA, occupancy_threshold: float) -> pyhmmer.easel.TextMSA:
    """Trim the low-occupancy columns from both ends of an alignment.

    Every column that clears the threshold is kept, including the outermost ones. When no
    column clears it the alignment is returned unchanged, because there is no meaningful
    span to keep and an empty alignment would crash the next `hmmbuild`.

    Distinct from `clip_env_ends`, which reads the RF line rather than gap counts, and
    runs earlier in the round.
    """
    sequence_matrix = np.array([list(row) for row in msa.alignment])
    positions = calculate_trim_positions(sequence_matrix, occupancy_threshold)
    if positions is None:
        return msa
    start_position, end_position = positions
    # end_position is inclusive; range() excludes its stop value.
    return msa.select(columns=range(start_position, end_position + 1))


def extract_first_part(sequence_name: str) -> str:
    return sequence_name.split("/")[0]


def unmask_sequence_names(sequences: Iterable[Sequence]) -> list[str]:
    return [extract_first_part(name) for name, _ in sequences]


def check_seed_membership(
    original_sequence_names: Iterable[str], filtered_sequence_names: Iterable[str]
) -> float:
    """Return the fraction of the cluster's distinct proteins still recruited.

    Both sides are counted after `extract_first_part` and as sets, so the ratio cannot
    exceed 1. Dividing by the raw row count instead would let a cluster TSV that repeats
    a member report less than full membership for a family that kept every one of them.
    """
    original_first_parts = set(map(extract_first_part, original_sequence_names))
    filtered_first_parts = set(map(extract_first_part, filtered_sequence_names))
    return len(original_first_parts & filtered_first_parts) / len(original_first_parts)


def parse_protein_name(row_name: str, aligned_row: str, indexed_sequences: IndexedSequences) -> str:
    """Rename one alignment row to `<protein>/<start>-<end>` on the parent protein.

    Database records are themselves slices of a protein, named `<protein>_<start>_<end>`
    with 1-based inclusive bounds. `mask_sequence` clips a record further to a hit
    envelope and appends `/<env_from>_<env_to>`, relative to the record. Column trimming
    (`clip_env_ends`, `clip_ends`, `run_pytrimal_reps`) then drops residues from either
    end of the row, by a different amount per row, so the row's offset within its record
    can only be recovered by locating its residues.

    The envelope bounds the search. Locating the residues in the whole record -- what the
    legacy script did -- returns the first match, so two identical repeat domains of one
    protein resolved to the same name and one of them was dropped as a duplicate.
    """
    record_name, _, envelope = row_name.partition("/")
    record = cast(Sequence, indexed_sequences.get(record_name)).seq
    residues = re.sub(r"[.\-~]", "", aligned_row).upper()

    if envelope:
        env_from, env_to = (int(bound) for bound in envelope.split("_"))
    else:
        env_from, env_to = 1, len(record)
    offset = record.find(residues, env_from - 1, env_to)
    if offset < 0:
        raise ValueError(f"{row_name}: aligned residues are not in its envelope of {record_name}")

    splits = record_name.split("_")
    if len(splits) != 3:
        # A name without slice bounds is a whole protein. When the row spans all of it the
        # legacy script emitted the bare accession, which `family_metadata` records as
        # region "-"; that convention is downstream-visible, so it is kept.
        if len(residues) == len(record):
            return record_name
        return f"{splits[0]}/{offset + 1}-{offset + len(residues)}"
    start = offset + int(splits[1])
    return f"{splits[0]}/{start}-{start + len(residues) - 1}"


def renumber_msa(
    msa: pyhmmer.easel.TextMSA, family_name: str, indexed_sequences: IndexedSequences
) -> pyhmmer.easel.TextMSA:
    """Return a copy of `msa` whose rows are named in parent-protein coordinates.

    Rebuilt rather than renamed in place: `MSA.names` is read-only and `.sequences` hands
    back copies, so assigning to a row's name is silently discarded. Rebuilding also drops
    hmmalign's `#=GR PP` and `#=GC PP_cons` annotation, which would otherwise double the
    size of every stored full MSA.
    """
    renumbered = pyhmmer.easel.TextMSA(
        name=family_name.encode(),
        sequences=[
            pyhmmer.easel.TextSequence(
                name=parse_protein_name(name, row, indexed_sequences), sequence=row
            )
            for name, row in zip(msa.names, msa.alignment, strict=True)
        ],
    )
    renumbered.reference = msa.reference
    return renumbered


class FamilyState(Enum):
    RUNNING = auto()
    CONVERGED = auto()
    DISCARDED = auto()
    SUCCESSFUL = auto()


@dataclass
class Family:
    representative: str
    members: list[str]
    state: FamilyState = FamilyState.RUNNING
    seed_msa: pyhmmer.easel.DigitalMSA | None = None
    hmm: pyhmmer.plan7.HMM | None = None
    records: list[Record] = field(default_factory=list)
    qlen: int = 0
    total_checked_sequences: set[str] = field(default_factory=set)
    full_msa: pyhmmer.easel.TextMSA | None = None
    full_msa_num_seqs: int = 0
    discard_reason: str = ""
    discard_value: int | float = 0.0
    ever_converged: bool = False
    family_id: int | None = None

    def discard(self, reason: str, value: int | float) -> None:
        self.state = FamilyState.DISCARDED
        self.discard_reason = reason
        self.discard_value = value

    def initialise(self, indexed_sequences: IndexedSequences, cpus: int) -> None:
        """Build the first seed MSA from the cluster's own members.

        Discards on fewer than two sequences, which `hmmbuild` cannot model and which
        surfaces as `eslEMEM (status code 5)`. The bound is `< 2`, not `<= 2`: raising
        it here would newly discard every two-member cluster, a change to the science
        that the `<= 2` rule after redundancy trimming does not imply.
        """
        self.seed_msa = run_initial_msa(self.members, indexed_sequences, cpus)
        sequence_count = len(self.seed_msa.names) if self.seed_msa is not None else 0
        if sequence_count < 2:
            self.discard("too few sequences before initial hmmbuild", sequence_count)

    def advance(
        self, options: argparse.Namespace, indexed_sequences: IndexedSequences, round_number: int
    ) -> None:
        """Consume this round's hits and move the family to its next state.

        Writes nothing. Every side effect is deferred to `emit_family`, because a family
        that converges in round 1 must not record itself before an earlier family that
        is still running in round 3 -- family ids depend on cluster order, not on the
        order in which families finish.

        The seed MSA the family carries out of this method is always the one that built
        the model `finish` and `emit_family` will use. Two branches keep that true:

        * On convergence the seed MSA is left untouched, so the hand build sees the
          alignment that produced the converged model.
        * On `MAX_ROUNDS` the whole re-align/trim tail is skipped. Its only consumer is
          the next round's `hmmbuild`, and there is no next round. Running it anyway --
          as the legacy script did -- left the family carrying a seed one generation
          ahead of its own model, so the exported HMM described an alignment that was
          never searched with, while the full MSA beside it came from the older model.
        """
        filtered_sequences = filter_hits(
            self.records,
            self.qlen,
            False,
            options.recruit_hit_length_percentage,
            indexed_sequences,
        )
        if not filtered_sequences:
            self.discard("low complexity model - confounding cluster", 0.0)
            return

        recruited_names = set(unmask_sequence_names(filtered_sequences))
        new_recruited_sequences = recruited_names - self.total_checked_sequences
        self.total_checked_sequences.update(new_recruited_sequences)
        if not new_recruited_sequences:
            self.state = FamilyState.CONVERGED
            self.ever_converged = True
            return
        if round_number == MAX_ROUNDS:
            return

        aligned = clip_env_ends(
            run_hmmalign(cast(pyhmmer.plan7.HMM, self.hmm), filtered_sequences, options.cpus)
        )
        _, representative_length = msa_stats(aligned)
        reason = check_rep_length(representative_length, options)
        if reason is not None:
            self.discard(*reason)
            return

        trimmed = run_pytrimal_reps(aligned, options.max_seq_identity, options.max_seed_seqs)
        sequence_count = len(trimmed.names)
        if sequence_count <= 2:
            self.discard("too few sequences after redundancy filtering", sequence_count)
            return
        trimmed = clip_ends(trimmed, options.max_gap_occupancy)
        self.seed_msa = trimmed.digitize(ALPHABET)

    def finish(self, options: argparse.Namespace, indexed_sequences: IndexedSequences) -> None:
        """Run the exit branch: relax the envelope filter, then build the full MSA.

        No new search is performed. The exit branch's model is always the one the last
        round already searched with -- on convergence because that is the round that
        converged, and after `MAX_ROUNDS` because the loop stops before another build --
        so this re-filters the cached records with `exit_flag` instead, at the cost of a
        whole database pass saved per family.
        """
        if self.state not in (FamilyState.RUNNING, FamilyState.CONVERGED):
            return

        filtered_sequences = filter_hits(
            self.records,
            self.qlen,
            True,
            options.recruit_hit_length_percentage,
            indexed_sequences,
        )
        if not filtered_sequences:
            self.discard("low complexity model - confounding cluster", 0.0)
            return

        membership = check_seed_membership(self.members, unmask_sequence_names(filtered_sequences))
        if membership < options.discard_min_starting_membership:
            self.discard("few seed sequences remained", membership)
            return

        self.full_msa = run_hmmalign(
            cast(pyhmmer.plan7.HMM, self.hmm), filtered_sequences, options.cpus
        )
        self.full_msa_num_seqs, representative_length = msa_stats(self.full_msa)
        reason = check_rep_length(representative_length, options)
        if reason is not None:
            self.discard(*reason)
            return
        self.state = FamilyState.SUCCESSFUL


def check_rep_length(
    representative_length: int, options: argparse.Namespace
) -> tuple[str, int] | None:
    if representative_length < options.discard_min_rep_length:
        return "family representative length too small", representative_length
    if representative_length > options.discard_max_rep_length:
        return "family representative length too large", representative_length
    return None


@dataclass
class Writers:
    root: Path
    indexed: IndexedSequences
    refined_families: IO[str]
    discarded_clusters: IO[str]
    successful_clusters: IO[str]
    converged_families: IO[str]
    family_metadata: IO[str]
    family_representatives: IO[str]


@contextlib.contextmanager
def deterministic_gzip_text(path: Path) -> Iterator[IO[str]]:
    with (
        path.open("wb") as raw,
        gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed,
        io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text,
    ):
        yield text


@contextlib.contextmanager
def deterministic_gzip_binary(path: Path) -> Iterator[IO[bytes]]:
    with (
        path.open("wb") as raw,
        gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed,
    ):
        yield compressed


def emit_family(
    family: Family,
    success_count: int,
    chunk: str,
    writers: Writers,
) -> int:
    """Write a finished family's artifacts and return the new successful-family count.

    Must be called once per family, in cluster-file order: `provisional_id` is derived
    from how many families succeeded before this one.

    Only successful families are written to `converged_families`; discarded families
    never receive an id or appear in that file.
    """
    provisional_id = success_count + 1
    if family.state is FamilyState.DISCARDED:
        writers.discarded_clusters.write(
            f"{family.representative},{family.discard_reason},{family.discard_value}\n"
        )
        return success_count

    family.family_id = provisional_id
    if family.ever_converged:
        writers.converged_families.write(f"{provisional_id}\n")
    family_name = f"{chunk}_{provisional_id}"
    seed_msa = cast(pyhmmer.easel.DigitalMSA, family.seed_msa)
    full_msa = cast(pyhmmer.easel.TextMSA, family.full_msa)
    # No size guard: round 1 can never converge, so every seed reaching here came from a
    # `run_pytrimal_reps` that passed `advance`'s `<= 2` check and `hmmbuild` cannot hit
    # eslEMEM. A guard would have to record a discard, and discards are results.
    final_hmm = run_hmmbuild(seed_msa, family_name, hand=True)
    final_hmm.name = family_name
    # Builder stamps DATE from the wall clock and COM from sys.argv. Both are serialised
    # into the HMM, and DATE is formatted through the locale. Dropping them is what makes
    # the file byte-reproducible; the family name survives as NAME.
    final_hmm.creation_time = None
    final_hmm.command_line = None

    writers.successful_clusters.write(f"{family.representative}\n")
    if seed_msa.reference is None:
        raise ValueError("successful seed MSA has no RF reference annotation")
    (writers.root / "rf" / f"{family_name}.txt").write_text(seed_msa.reference, encoding="utf-8")
    with deterministic_gzip_binary(writers.root / "hmm" / f"{family_name}.hmm.gz") as handle:
        final_hmm.write(handle)

    indexed = cast(IndexedSequences, writers.indexed)
    renumbered_seed = renumber_msa(seed_msa.textize(), family_name, indexed)
    renumbered_full = renumber_msa(full_msa, family_name, indexed)
    for directory, msa in (
        ("seed_msa", renumbered_seed),
        ("full_msa", renumbered_full),
    ):
        path = writers.root / directory / f"{family_name}.sto.gz"
        with deterministic_gzip_binary(path) as handle:
            msa.write(handle, format="pfam")

    for row_number, (name, row) in enumerate(
        zip(renumbered_full.names, renumbered_full.alignment, strict=True)
    ):
        sequence_name = name.decode() if isinstance(name, bytes) else name
        writers.refined_families.write(f"{provisional_id}\t{sequence_name}\n")
        if row_number == 0:
            # Row 0 is the representative: hits arrive in HMMER's ranking order.
            residues = re.sub(r"[.\-~]", "", row).upper()
            protein, _, region = sequence_name.partition("/")
            writers.family_metadata.write(
                f'{provisional_id},{family.full_msa_num_seqs},"{protein}",{region or "-"},'
                f"{len(residues)},{residues},{final_hmm.consensus},{family.ever_converged}\n"
            )
            writers.family_representatives.write(
                f">{sequence_name}\t{chunk}_{provisional_id}\n{residues}\n"
            )
    return provisional_id


def load_clusters(path: str | os.PathLike[str]) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not all(fields):
                raise ValueError(
                    f"cluster TSV row {line_number} must contain exactly two non-empty fields"
                )
            representative, member = fields
            clusters.setdefault(representative, []).append(member)
    return clusters


def parse_args(args: SequenceCollection[str] | None = None) -> argparse.Namespace:
    # prog is spelled out because this parser is reached through the `mgnifam`
    # dispatcher, which would otherwise leave argv[0] in the usage line.
    parser = argparse.ArgumentParser(
        prog="mgnifam generate_families",
        description="Build protein families from a chunk of sequence clusters.",
    )
    # Only the two inputs are required; every threshold defaults to the value the README
    # documents, so `--help` doubles as the reference for what a plain run does.
    parser.add_argument("-c", "--clusters_chunk", required=True)
    parser.add_argument("-f", "--fasta_file", required=True)
    parser.add_argument("-p", "--cpus", type=int, default=8)
    parser.add_argument("-n", "--chunk_num", default="1")
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


def is_gzipped(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def validate_inputs(options: argparse.Namespace) -> dict[str, list[str]]:
    if CHUNK_PATTERN.fullmatch(options.chunk_num) is None:
        raise ValueError("chunk_num must match [A-Za-z0-9._-]+")
    fasta = Path(options.fasta_file)
    if not fasta.is_file():
        raise ValueError("fasta_file must be an existing file")
    if is_gzipped(fasta):
        raise ValueError("fasta_file must be uncompressed; SSI cannot seek in gzip streams")
    if not Path(options.clusters_chunk).is_file():
        raise ValueError("clusters_chunk must be an existing file")
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
    return load_clusters(options.clusters_chunk)


def prepare_output_directories(root: Path, chunk: str) -> None:
    """Create the output tree and clear this chunk's per-family artifacts.

    Without the clearing step, a rerun that produces fewer families leaves the surplus
    behind and the directory mixes two runs.
    """
    for directory in FAMILY_DIRECTORIES:
        (root / directory).mkdir(parents=True, exist_ok=True)
    # An exact numeric suffix, not a `<chunk>_*` glob: chunk "foo" would otherwise
    # delete "foo_bar_1", which belongs to chunk "foo_bar".
    artifact_pattern = re.compile(rf"{re.escape(chunk)}_\d+\..*")
    for directory in FAMILY_DIRECTORIES:
        for path in (root / directory).iterdir():
            if path.is_file() and artifact_pattern.fullmatch(path.name):
                path.unlink()


def configure_logger(path: Path) -> logging.Logger:
    logger = logging.getLogger(f"mgnifam.generate_families.{path.stem}.{id(path)}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def resolve_index(options: argparse.Namespace, root: Path) -> Path:
    """Return the SSI index path, building it if absent or older than the FASTA.

    Production runs should pass `--fasta_index` to share one index across chunk tasks;
    otherwise every task re-indexes the whole database under its output directory.
    """
    fasta = Path(options.fasta_file)
    index = Path(options.fasta_index) if options.fasta_index else root / f"{fasta.name}.ssi"
    if not index.exists() or index.stat().st_mtime_ns < fasta.stat().st_mtime_ns:
        build_ssi_index(fasta, index)
    return index


def main(args: SequenceCollection[str] | None = None) -> None:
    options = parse_args(args)
    clusters = validate_inputs(options)
    if options.batch_size == 0:
        options.batch_size = 2 * options.cpus

    root = options.output_dir
    prepare_output_directories(root, options.chunk_num)
    index_path = resolve_index(options, root)
    logger = configure_logger(root / f"{options.chunk_num}.log")

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

            writers = Writers(
                root=root,
                indexed=indexed_sequences,
                refined_families=stack.enter_context(
                    (root / f"{options.chunk_num}_families.tsv").open("w")
                ),
                discarded_clusters=stack.enter_context(
                    (root / f"{options.chunk_num}_discarded.csv").open("w")
                ),
                successful_clusters=stack.enter_context(
                    (root / f"{options.chunk_num}_successful.txt").open("w")
                ),
                converged_families=stack.enter_context(
                    (root / f"{options.chunk_num}_converged.txt").open("w")
                ),
                family_metadata=stack.enter_context(
                    (root / f"{options.chunk_num}_metadata.csv").open("w")
                ),
                family_representatives=stack.enter_context(
                    deterministic_gzip_text(root / f"{options.chunk_num}_reps.fasta.gz")
                ),
            )
            success_count = 0
            for batch_number, batch in enumerate(
                itertools.batched(clusters.items(), options.batch_size), 1
            ):
                active = [Family(representative, members) for representative, members in batch]
                for family in active:
                    family.initialise(indexed_sequences, options.cpus)

                for round_number in range(1, MAX_ROUNDS + 1):
                    running = [family for family in active if family.state is FamilyState.RUNNING]
                    if not running:
                        break
                    hmms = [
                        run_hmmbuild(
                            cast(pyhmmer.easel.DigitalMSA, family.seed_msa),
                            f"pending_{batch_number}_{round_number}_{index}",
                        )
                        for index, family in enumerate(running)
                    ]
                    with search(
                        hmms,
                        targets,
                        cpus=options.cpus,
                        evalue=options.recruit_evalue_cutoff,
                        logger=logger,
                        batch_number=batch_number,
                        round_number=round_number,
                    ) as searched:
                        for family, hmm, hits in zip(running, hmms, searched, strict=True):
                            family.hmm = hmm
                            family.qlen = hits.query.M
                            family.records = extract_records(hits)
                            family.advance(options, indexed_sequences, round_number)

                for family in active:
                    family.finish(options, indexed_sequences)
                for family in active:
                    success_count = emit_family(family, success_count, options.chunk_num, writers)
        logger.info("DONE.")
    finally:
        for handler in logger.handlers:
            handler.close()
        logger.handlers.clear()


if __name__ == "__main__":
    main()
