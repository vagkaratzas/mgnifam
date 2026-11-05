#!/usr/bin/env python3
"""Family similarity and clans from residue-level overlap of family members.

Two families are similar when their member regions cover the same residues of the same
proteins. For family F, `R(F)` is the set of (protein, position) pairs its member regions
cover, after merging each protein's overlapping or adjacent regions, so repeat domains and
duplicate rows count once. For every pair of families that share a protein this command
computes, exactly:

- `shared_proteins` S, and `protein_jaccard = S / (n_A + n_B - S)` over distinct proteins;
- `residue_jaccard = I / (|R_A| + |R_B| - I)` with `I = |R_A ∩ R_B|`;
- `residue_containment = I / min(|R_A|, |R_B|)`, which catches a small family nested in a
  large one.

Clans are the connected components of the pairs at `residue_jaccard >= --clan_jaccard`.

No HMM search and no candidate filter is involved: every pair that shares a protein is
scored. The cost is the number of family pairs per protein, `sum k(k-1)/2`, not the number
of family pairs, and memory is bounded by `--worker_memory_mb` per worker plus a dense
upper triangle of per-pair sums (12 bytes a pair).

The work runs in three phases, all integer arithmetic, so outputs do not depend on
`--cpus`, `--worker_memory_mb`, input order or `PYTHONHASHSEED`:

1. Members are parsed in parallel into sorted runs on scratch disk, keyed by a 128-bit
   blake2b hash of the protein name. If any member is a whole protein, the FASTA is parsed
   the same way into protein lengths.
2. Hash ranges sized to the worker budget are processed in parallel: intervals merged per
   (protein, family), then shared-protein and residue-overlap sums added to the triangle.
3. The triangle is scanned for edges and clans, and the outputs are written. The output
   tree is not touched before this phase, so a rejected run leaves a previous result as it
   was.
"""

import argparse
import array
import bisect
import contextlib
import gzip
import hashlib
import logging
import multiprocessing
import os
import re
import shutil
import sys
import tempfile
import time
from collections import Counter
from collections.abc import Iterable, Iterator
from collections.abc import Sequence as SequenceCollection
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

from mgnifam import __version__
from mgnifam.generate_families import (
    CHUNK_PATTERN,
    ChunkCorrupted,
    configure_logger,
    is_gzipped,
    write_stats,
)
from mgnifam.update_families import FAMILY_NAME_PATTERN, TRAVERSAL_NAMES

EDGES_HEADER = (
    "family_1,family_2,shared_proteins,protein_jaccard,residue_jaccard,residue_containment\n"
)
CLANS_HEADER = "clan_id,family_id,is_representative\n"
# Coordinates are stored as uint32 and widened to int64 before any arithmetic, so the
# largest accepted value leaves room for `end + 1` and differences without wrapping.
MAX_COORDINATE = 2**31 - 1
# Thresholds take at most this many decimal places, which bounds the denominator in the
# exact `num * q >= p * den` comparisons.
THRESHOLD_SCALE = 10**6
INT64_LIMIT = 2**63
NEAR_REDUNDANT = Fraction(95, 100)
SCIENTIFIC_PARAMETERS = ("min_edge_jaccard", "min_edge_containment", "clan_jaccard")
# 2^12 buckets by the top bits of the protein hash. Units group consecutive buckets by
# their actual counts, and a bucket over budget is split further by hash sub-range.
BUCKET_BITS = 12
LOCK_STRIPES = 64

RECORD = np.dtype([("h1", "<u8"), ("h2", "<u8"), ("fam", "<u4"), ("start", "<u4"), ("end", "<u4")])
LENGTH_RECORD = np.dtype([("h1", "<u8"), ("h2", "<u8"), ("length", "<u8")])
# Rows are parsed as bytes; the family rule is update_families' own, not a copy of it.
FAMILY_BYTES = re.compile(FAMILY_NAME_PATTERN.pattern.encode())
TRAVERSAL_BYTES = frozenset(name.encode() for name in TRAVERSAL_NAMES)


class InputError(ValueError):
    """Invalid input. Exit status 2; the output tree is left untouched."""


class ResourceError(RuntimeError):
    """The input cannot be processed within the worker budget. Exit status 1."""


class MissingLengths(InputError):
    """Whole-protein members whose protein has no FASTA record."""

    def __init__(self, count: int, keys: list[tuple[int, int]]) -> None:
        super().__init__(f"{count} whole-protein member rows have no FASTA record")
        self.count = count
        self.keys = keys

    def __reduce__(self) -> tuple[type, tuple[int, list[tuple[int, int]]]]:
        # Raised in a pool worker and re-raised in the parent, so it must pickle.
        return MissingLengths, (self.count, self.keys)


class Sizes(NamedTuple):
    """Work sizes derived from `--worker_memory_mb`, computed once in the parent."""

    range_bytes: int  # plain-text bytes per parse task
    run_rows: int  # records buffered before a sorted run is written
    unit_cost: int  # bytes of rows one Phase B unit may load
    emission_batch: int  # pair emissions per reduction batch


MEMBER_ROW_COST = 192  # a loaded member row with its sorted copies and int64 promotions
LENGTH_ROW_COST = 96
EMISSION_COST = 96
# A buffered member row (16 key bytes, three 8-byte values) plus the 28-byte record it
# becomes at flush peaks near 68 bytes; the rest is allocator over-allocation and slack.
RUN_ROW_COST = 96


def sizes_for(worker_memory_mb: int) -> Sizes:
    budget = worker_memory_mb * 2**20
    return Sizes(
        range_bytes=max(2**16, budget // 4),
        run_rows=max(1, budget // (2 * RUN_ROW_COST)),
        unit_cost=budget // 2,
        emission_batch=max(1, budget // (2 * EMISSION_COST)),
    )


# --- parsing -----------------------------------------------------------------------------


def parse_threshold(text: str) -> Fraction:
    """Parse a threshold in [0, 1] with at most six decimal places, exactly."""
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        raise argparse.ArgumentTypeError(f"{text!r} is not a number") from None
    if not 0 <= value <= 1:
        raise argparse.ArgumentTypeError(f"{text} must be in [0, 1]")
    if (value * THRESHOLD_SCALE).denominator != 1:
        raise argparse.ArgumentTypeError(f"{text} has more than six decimal places")
    return value


def split_member(member: bytes) -> tuple[bytes, int, int] | None:
    """Split `<protein>/<start>-<end>` at its last slash, or return None for a bare name.

    The same spelling `generate_families` emits and `SLICE_SUFFIX` matches: the protein
    keeps any earlier slashes, and only a trailing `<digits>-<digits>` is a range. A bare
    name is a whole-protein member.
    """
    protein, slash, bounds = member.rpartition(b"/")
    if not slash or not protein:
        return None
    start, dash, end = bounds.partition(b"-")
    if not dash or not start.isdigit() or not end.isdigit():
        return None
    return protein, int(start), int(end)


def protein_key(protein: bytes) -> bytes:
    """128-bit protein key, independent of `PYTHONHASHSEED` and of the process."""
    return hashlib.blake2b(protein, digest_size=16).digest()


@dataclass
class RunWriter:
    """Buffers fixed-width records and writes them as runs sorted by `h1`."""

    directory: Path
    prefix: str
    dtype: np.dtype
    run_rows: int
    keys: bytearray = field(default_factory=bytearray)
    # 8 bytes a value; a list of Python ints costs ~36.
    columns: dict[str, array.array] = field(default_factory=dict)
    runs: list[tuple[str, int]] = field(default_factory=list)
    counts: np.ndarray = field(default_factory=lambda: np.zeros(2**BUCKET_BITS, np.int64))

    def __post_init__(self) -> None:
        self.columns = {name: array.array("q") for name in self.dtype.names[2:]}  # type: ignore[index]

    def add(self, key: bytes, *values: int) -> None:
        self.keys += key
        for column, value in zip(self.columns.values(), values, strict=True):
            column.append(value)
        if len(self.keys) >= 16 * self.run_rows:
            self.flush()

    def flush(self) -> None:
        rows = len(self.keys) // 16
        if not rows:
            return
        records = np.empty(rows, self.dtype)
        hashes = np.frombuffer(self.keys, dtype="<u8").reshape(rows, 2)
        records["h1"], records["h2"] = hashes[:, 0], hashes[:, 1]
        del hashes  # release the buffer export so the bytearray can be cleared
        self.keys.clear()
        for name, column in self.columns.items():
            records[name] = column
            del column[:]
        # The buffers are freed before the sort, so its copies do not stack on them.
        records = records[np.argsort(records["h1"], kind="stable")]
        path = self.directory / f"{self.prefix}_{len(self.runs):05d}.bin"
        records.tofile(path)
        self.runs.append((str(path), rows))
        self.counts += np.bincount(
            (records["h1"] >> np.uint64(64 - BUCKET_BITS)).astype(np.int64),
            minlength=2**BUCKET_BITS,
        )


class Task(NamedTuple):
    number: int
    path: str
    start: int  # byte range of a plain file; (0, -1) streams a gzip file
    end: int


def spread(range_bytes: int, paths: Iterable[str], cpus: int) -> int:
    """Shrink plain-text ranges so a mid-sized input still yields ~4 tasks per worker."""
    total = sum(os.path.getsize(path) for path in paths)
    return max(2**16, min(range_bytes, total // (4 * cpus)))


def plan_tasks(paths: Iterable[str], range_bytes: int) -> list[Task]:
    tasks: list[Task] = []
    for path in paths:
        if is_gzipped(Path(path)):
            tasks.append(Task(len(tasks), path, 0, -1))
            continue
        size = os.path.getsize(path)
        for start in range(0, max(size, 1), range_bytes):
            tasks.append(Task(len(tasks), path, start, min(start + range_bytes, size)))
    return tasks


def task_lines(task: Task, bounded: bool = True) -> Iterator[tuple[int, bytes]]:
    """Yield `(position, line)` from the first line that starts inside the task's range.

    Position is a byte offset in a plain file and a 1-based line number in a gzip stream.
    `bounded` stops before the first line starting at or past the range end; the FASTA
    reader passes False so a record that starts in range is read to its end.
    """
    if task.end < 0:
        with gzip.open(task.path, "rb") as handle:
            yield from enumerate(handle, start=1)
        return
    with open(task.path, "rb") as handle:
        position = task.start
        if position:
            handle.seek(position - 1)
            position += len(handle.readline()) - 1
        while not bounded or position < task.end:
            line = handle.readline()
            if not line:
                return
            yield position, line
            position += len(line)


def location(task: Task, position: int) -> str:
    unit = "line" if task.end < 0 else "byte offset"
    return f"{task.path} ({unit} {position})"


class MembersResult(NamedTuple):
    number: int
    names_path: str
    runs: list[tuple[str, int]]
    counts: np.ndarray
    rows: int
    bare_rows: int


def parse_members_task(task: Task, scratch: str, run_rows: int) -> MembersResult:
    """Parse one task's member rows into sorted runs. Families get task-local indices."""
    writer = RunWriter(Path(scratch) / "members", f"m{task.number:06d}", RECORD, run_rows)
    families: dict[bytes, int] = {}
    rows = bare_rows = 0
    for position, raw in task_lines(task):
        line = raw.rstrip(b"\r\n")
        family, tab, member = line.partition(b"\t")
        if not tab or not member or b"\t" in member:
            raise InputError(f"{location(task, position)}: expected 'family<TAB>member'")
        index = families.get(family)
        if index is None:
            if FAMILY_BYTES.fullmatch(family) is None or family in TRAVERSAL_BYTES:
                raise InputError(
                    f"{location(task, position)}: family {family!r} must match "
                    "[A-Za-z0-9._-]+ and be neither '.' nor '..'"
                )
            index = families[family] = len(families)
        sliced = split_member(member)
        if sliced is None:
            writer.add(protein_key(member), index, 1, 0)  # end 0: length from the FASTA
            bare_rows += 1
        else:
            protein, start, end = sliced
            if not 1 <= start <= end <= MAX_COORDINATE:
                raise InputError(
                    f"{location(task, position)}: region {start}-{end} must satisfy "
                    f"1 <= start <= end <= {MAX_COORDINATE}"
                )
            writer.add(protein_key(protein), index, start, end)
        rows += 1
    writer.flush()
    names_path = Path(scratch) / f"names_{task.number:06d}.txt"
    names_path.write_bytes(b"\n".join(families))
    return MembersResult(task.number, str(names_path), writer.runs, writer.counts, rows, bare_rows)


class LengthsResult(NamedTuple):
    runs: list[tuple[str, int]]
    counts: np.ndarray
    records: int
    max_length: int


def parse_fasta_task(task: Task, scratch: str, run_rows: int) -> LengthsResult:
    """Record `(key, length)` for every FASTA record whose header starts in the task's range.

    Lines before the task's first header belong to the previous task's record. Lengths
    are kept for every record: which ones are needed is decided per hash range later, so
    no set of names is ever held in memory.
    """
    writer = RunWriter(Path(scratch) / "fasta", f"f{task.number:06d}", LENGTH_RECORD, run_rows)
    records = max_length = length = 0
    name: bytes | None = None
    for position, line in task_lines(task, bounded=False):
        if line.startswith(b">"):
            if task.end >= 0 and position >= task.end:
                break
            if name is not None:
                writer.add(protein_key(name), length)
                records += 1
                max_length = max(max_length, length)
            header = line[1:].split(None, 1)
            if not header:
                raise InputError(f"{location(task, position)}: FASTA header has no name")
            name, length = header[0], 0
        elif name is not None:
            length += len(line.strip())
        elif task.start == 0 and line.strip():
            raise InputError(f"{location(task, position)}: FASTA must start with '>'")
    if name is not None:
        writer.add(protein_key(name), length)
        records += 1
        max_length = max(max_length, length)
    writer.flush()
    return LengthsResult(writer.runs, writer.counts, records, max_length)


# --- per-unit work -----------------------------------------------------------------------


def group_ids(h1: np.ndarray, h2: np.ndarray) -> tuple[np.ndarray, int]:
    """Dense ids for distinct `(h1, h2)` keys, numbered in key order."""
    if not len(h1):
        return np.zeros(0, np.int64), 0
    order = np.lexsort((h2, h1))
    sorted_h1, sorted_h2 = h1[order], h2[order]
    change = (sorted_h1[1:] != sorted_h1[:-1]) | (sorted_h2[1:] != sorted_h2[:-1])
    ids_sorted = np.concatenate(([0], np.cumsum(change)))
    ids = np.empty(len(h1), np.int64)
    ids[order] = ids_sorted
    return ids, int(ids_sorted[-1]) + 1


def merge_intervals(
    protein: np.ndarray, family: np.ndarray, start: np.ndarray, end: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Merge overlapping or adjacent intervals per (protein, family).

    Inputs are int64 arrays; returns `(protein, family, start, end)` of the merged
    intervals sorted by protein, family and start. One family's merged intervals on a
    protein are disjoint and non-adjacent, which the overlap sweep relies on.
    """
    if not len(start):
        empty = np.zeros(0, np.int64)
        return empty, empty, empty, empty
    order = np.lexsort((start, family, protein))
    protein, family, start, end = protein[order], family[order], start[order], end[order]
    group_start = np.concatenate(
        ([True], (protein[1:] != protein[:-1]) | (family[1:] != family[:-1]))
    )
    group = np.cumsum(group_start) - 1
    # The group in the high bits keeps the running maximum from crossing groups.
    running = np.maximum.accumulate((group << 32) | end) & 0xFFFFFFFF
    previous = np.concatenate(([0], running[:-1]))
    new = group_start | (start > previous + 1)
    first = np.flatnonzero(new)
    return protein[first], family[first], start[first], np.maximum.reduceat(end, first)


def virtual_slice(
    counts: np.ndarray, cumulative: np.ndarray, low: int, high: int
) -> tuple[np.ndarray, np.ndarray]:
    """Elements `[low, high)` of the concatenation of `range(counts[r])` for every row r.

    Returns `(row, offset)` per element, so a batch may begin or end inside a row.
    """
    first = int(np.searchsorted(cumulative, low, side="right"))
    last = int(np.searchsorted(cumulative, high - 1, side="right"))
    rows = np.arange(first, last + 1)
    row_start = cumulative[rows] - counts[rows]
    taken = np.minimum(cumulative[rows], high) - np.maximum(row_start, low)
    row = np.repeat(rows, taken)
    offset = np.arange(low, high) - np.repeat(row_start, taken)
    return row, offset


def pair_index(a: np.ndarray, b: np.ndarray, families: int) -> np.ndarray:
    """Index of the pair `a < b` in the row-major upper triangle of `families` families."""
    return a * families - a * (a + 1) // 2 + (b - a - 1)


@dataclass
class UnitState:
    """What every Phase B worker needs, set once per process by `init_unit_worker`."""

    member_runs: list[tuple[str, int, int]]  # path, rows, remap offset
    length_runs: list[tuple[str, int]]
    remap_path: str
    families: int
    shared_path: str | None
    shared_dtype: str
    overlap_path: str | None
    emission_batch: int
    locks: list[Any]


_UNIT: dict[str, Any] = {}


def init_unit_worker(state: UnitState) -> None:
    _UNIT.clear()
    _UNIT["state"] = state
    _UNIT["maps"] = {}
    _UNIT["remap"] = np.memmap(state.remap_path, dtype="<u4", mode="r") if state.families else None
    pairs = state.families * (state.families - 1) // 2
    if state.shared_path and pairs:
        _UNIT["shared"] = np.memmap(state.shared_path, dtype=state.shared_dtype, mode="r+")
        _UNIT["overlap"] = np.memmap(state.overlap_path, dtype="<u8", mode="r+")
        _UNIT["stripes"] = (np.arange(1, LOCK_STRIPES) * pairs) // LOCK_STRIPES


def mapped(path: str, dtype: np.dtype) -> np.ndarray:
    maps = _UNIT["maps"]
    if path not in maps:
        maps[path] = np.memmap(path, dtype=dtype, mode="r")
    return maps[path]


def key_range(keys: np.ndarray, low: int, high: int) -> tuple[int, int]:
    """Rows of a run sorted by `h1` whose `h1` lies in `[low, high)`; high may be 2^64.

    `keys` is a strided field view of a memmapped run. `np.searchsorted` would copy the
    whole column to a contiguous array on every call; `bisect` reads O(log n) elements.
    """
    first = bisect.bisect_left(keys, np.uint64(low))
    last = len(keys) if high >= 2**64 else bisect.bisect_left(keys, np.uint64(high))
    return first, last


def accumulate(target: np.ndarray, index: np.ndarray, values: np.ndarray) -> None:
    """Add a batch with unique, sorted indices into a shared triangle, stripe by stripe."""
    bounds = np.searchsorted(index, _UNIT["stripes"])
    edges = np.concatenate(([0], bounds, [len(index)]))
    for stripe in np.flatnonzero(edges[1:] > edges[:-1]):
        low, high = int(edges[stripe]), int(edges[stripe + 1])
        with _UNIT["state"].locks[stripe]:
            # Unique indices, so read-add-write per element is exact.
            target[index[low:high]] += values[low:high].astype(target.dtype)


class UnitResult(NamedTuple):
    proteins_per_family: np.ndarray
    residues_per_family: np.ndarray
    proteins: int
    pair_emissions: int
    overlap_emissions: int
    max_families_per_protein: int


def load_unit(low: int, high: int) -> tuple[np.ndarray, ...]:
    state: UnitState = _UNIT["state"]
    parts = []
    for path, _, remap_offset in state.member_runs:
        run = mapped(path, RECORD)
        first, last = key_range(run["h1"], low, high)
        if last > first:
            chunk = np.array(run[first:last])
            chunk["fam"] = _UNIT["remap"][remap_offset + chunk["fam"].astype(np.int64)]
            parts.append(chunk)
    rows = np.concatenate(parts) if parts else np.zeros(0, RECORD)
    h1, h2 = rows["h1"], rows["h2"]
    family = rows["fam"].astype(np.int64)
    start = rows["start"].astype(np.int64)
    end = rows["end"].astype(np.int64)
    bare = end == 0
    if bare.any():
        end[bare] = bare_lengths(h1[bare], h2[bare], low, high)
    protein, _ = group_ids(h1, h2)
    return protein, family, start, end


def bare_lengths(h1: np.ndarray, h2: np.ndarray, low: int, high: int) -> np.ndarray:
    """FASTA lengths for whole-protein rows, joined on the full 128-bit key."""
    state: UnitState = _UNIT["state"]
    parts = []
    for path, _ in state.length_runs:
        run = mapped(path, LENGTH_RECORD)
        first, last = key_range(run["h1"], low, high)
        if last > first:
            parts.append(np.array(run[first:last]))
    lengths = np.concatenate(parts) if parts else np.zeros(0, LENGTH_RECORD)
    ids, distinct = group_ids(
        np.concatenate((h1, lengths["h1"])), np.concatenate((h2, lengths["h2"]))
    )
    needed, records = ids[: len(h1)], ids[len(h1) :]
    copies = np.bincount(records, minlength=distinct)
    duplicated = np.flatnonzero(copies[np.unique(needed)] > 1)
    if len(duplicated):
        raise InputError("the FASTA holds a whole-protein member's record more than once")
    by_id = np.zeros(distinct, np.int64)
    by_id[records] = lengths["length"].astype(np.int64)
    found = copies[needed] > 0
    if not found.all():
        missing = np.flatnonzero(~found)
        keys = sorted({(int(h1[i]), int(h2[i])) for i in missing[:1000]})
        raise MissingLengths(len(missing), keys)
    result = by_id[needed]
    if ((result < 1) | (result > MAX_COORDINATE)).any():
        raise InputError(
            f"a whole-protein member's FASTA record is empty or longer than {MAX_COORDINATE}"
        )
    return result


def process_unit(bounds: tuple[int, int]) -> UnitResult:
    """Merge one hash range's intervals and add its pair sums into the triangle."""
    state: UnitState = _UNIT["state"]
    families = state.families
    protein, family, start, end = load_unit(*bounds)
    protein, family, start, end = merge_intervals(protein, family, start, end)

    proteins_per_family = np.zeros(families, np.int64)
    residues_per_family = np.zeros(families, np.int64)
    if not len(start):
        return UnitResult(proteins_per_family, residues_per_family, 0, 0, 0, 0)
    np.add.at(residues_per_family, family, end - start + 1)
    # (protein, family) groups, still in protein-then-family order.
    group_first = np.flatnonzero(
        np.concatenate(([True], (protein[1:] != protein[:-1]) | (family[1:] != family[:-1])))
    )
    group_protein, group_family = protein[group_first], family[group_first]
    np.add.at(proteins_per_family, group_family, 1)
    protein_first = np.flatnonzero(
        np.concatenate(([True], group_protein[1:] != group_protein[:-1]))
    )
    per_protein = np.diff(np.append(protein_first, len(group_protein)))
    proteins, max_k = len(protein_first), int(per_protein.max())

    pair_emissions = overlap_emissions = 0
    if "shared" in _UNIT:
        # Stream 1: every family pair on a protein shares it once. Families are sorted
        # within a protein, so the pair list is the upper triangle of each protein's list.
        rank = np.arange(len(group_family)) - np.repeat(protein_first, per_protein)
        counts = np.repeat(per_protein, per_protein) - 1 - rank
        pair_emissions = stream(
            counts, lambda row, offset: (row, row + 1 + offset, None), group_family, "shared"
        )
        # Stream 2: overlapping intervals, swept by start within each protein. A partner
        # of i starts no later than i ends; one family's merged intervals never overlap.
        order = np.lexsort((start, protein))
        s_protein, s_family = protein[order], family[order]
        s_start, s_end = start[order], end[order]
        reach = np.searchsorted((s_protein << 32) | s_start, (s_protein << 32) | s_end, "right")
        counts = reach - np.arange(len(s_start)) - 1

        def partners(row: np.ndarray, offset: np.ndarray) -> tuple[np.ndarray, ...]:
            other = row + 1 + offset
            return row, other, np.minimum(s_end[row], s_end[other]) - s_start[other] + 1

        overlap_emissions = stream(counts, partners, s_family, "overlap")
    return UnitResult(
        proteins_per_family,
        residues_per_family,
        proteins,
        pair_emissions,
        overlap_emissions,
        max_k,
    )


def stream(counts: np.ndarray, expand: Any, family: np.ndarray, target: str) -> int:
    """Emit a virtual array of pairs in bounded batches and add them into `target`."""
    state: UnitState = _UNIT["state"]
    cumulative = np.cumsum(counts)
    total = int(cumulative[-1]) if len(cumulative) else 0
    for low in range(0, total, state.emission_batch):
        row, offset = virtual_slice(counts, cumulative, low, min(low + state.emission_batch, total))
        first, second, values = expand(row, offset)
        a, b = family[first], family[second]
        index = pair_index(np.minimum(a, b), np.maximum(a, b), state.families)
        order = np.argsort(index, kind="stable")
        index = index[order]
        heads = np.flatnonzero(np.concatenate(([True], index[1:] != index[:-1])))
        if values is None:
            sums = np.diff(np.append(heads, len(index)))
        else:
            sums = np.add.reduceat(values[order], heads)
        accumulate(_UNIT[target], index[heads].astype(np.int64), sums)
    return total


# --- planning ----------------------------------------------------------------------------


def plan_units(
    member_runs: list[tuple[str, int, int]],
    length_runs: list[tuple[str, int]],
    member_counts: np.ndarray,
    length_counts: np.ndarray,
    unit_cost: int,
    target_cost: int | None = None,
) -> list[tuple[int, int]]:
    """Cut the hash space into ranges whose rows fit `unit_cost` bytes.

    `target_cost` (at most `unit_cost`) is a softer size for spreading work over workers:
    ranges are cut to it where possible, but only `unit_cost` can fail a run.

    Consecutive buckets are grouped by their actual counts. A bucket over budget is halved
    by hash range, counted with `searchsorted` on the runs, down to a single hash value.
    """
    shift = 64 - BUCKET_BITS
    maps: dict[str, np.ndarray] = {}

    def keys(path: str, dtype: np.dtype) -> np.ndarray:
        if path not in maps:
            maps[path] = np.memmap(path, dtype=dtype, mode="r")["h1"]
        return maps[path]

    def cost(low: int, high: int) -> int:
        total = 0
        for path, _, _ in member_runs:
            first, last = key_range(keys(path, RECORD), low, high)
            total += (last - first) * MEMBER_ROW_COST
        for path, _ in length_runs:
            first, last = key_range(keys(path, LENGTH_RECORD), low, high)
            total += (last - first) * LENGTH_ROW_COST
        return total

    target = min(unit_cost, target_cost or unit_cost)

    def split(low: int, high: int) -> Iterator[tuple[int, int]]:
        total = cost(low, high)
        if not total:
            return  # an empty half would still cost every run a key lookup per unit
        if total <= target or (high - low == 1 and total <= unit_cost):
            yield low, high
        elif high - low == 1:
            raise ResourceError(
                "one protein's member rows exceed the worker budget; raise --worker_memory_mb"
            )
        else:
            middle = (low + high) // 2
            yield from split(low, middle)
            yield from split(middle, high)

    weights = member_counts * MEMBER_ROW_COST + length_counts * LENGTH_ROW_COST
    units: list[tuple[int, int]] = []
    open_low, open_cost = None, 0
    for bucket, weight in enumerate(weights.tolist()):
        low, high = bucket << shift, (bucket + 1) << shift
        if not weight:
            continue
        if weight > target:
            if open_low is not None:
                units.append((open_low, low))
                open_low, open_cost = None, 0
            units.extend(split(low, high))
            continue
        if open_low is not None and open_cost + weight > target:
            units.append((open_low, low))
            open_low, open_cost = None, 0
        if open_low is None:
            open_low = low
        open_cost += weight
    if open_low is not None:
        units.append((open_low, 2**64))
    maps.clear()
    return units


# --- clans -------------------------------------------------------------------------------


def build_clans(
    families: int, first: np.ndarray, second: np.ndarray, weight: np.ndarray
) -> list[tuple[int, list[int]]]:
    """Connected components of at least two families, with their representatives.

    `first < second` index edges in (first, second) order; families are indexed in name
    order. The representative has the largest weighted degree (sum of edge weights),
    ties going to the earlier name. Clans are ordered by size, largest first, then by
    representative name. Returns `(representative, members)` with members in name order.
    """
    parent = list(range(families))

    def root(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for a, b in zip(first.tolist(), second.tolist(), strict=True):
        ra, rb = root(a), root(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    degree = np.zeros(families, np.float64)
    # np.add.at applies in index order, so the sums are reproducible.
    np.add.at(degree, first, weight)
    np.add.at(degree, second, weight)
    members: dict[int, list[int]] = {}
    for node in sorted(set(first.tolist()) | set(second.tolist())):
        members.setdefault(root(node), []).append(node)
    clans = []
    for nodes in members.values():
        best = max(nodes, key=lambda node: (degree[node], -node))
        clans.append((best, nodes))
    clans.sort(key=lambda clan: (-len(clan[1]), clan[0]))
    return clans


# --- command -----------------------------------------------------------------------------


def parse_args(args: SequenceCollection[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="mgnifam compute_clans",
        description="Family similarity and clans from residue-level member overlap.",
    )
    parser.add_argument(
        "-i",
        "--members",
        nargs="+",
        required=True,
        help="family<TAB>member TSVs, such as <chunk>_families.tsv; plain or gzip",
    )
    parser.add_argument(
        "-f",
        "--fasta_file",
        help="database FASTA; required when a member is a whole protein (a bare name)",
    )
    parser.add_argument("-p", "--cpus", type=int, default=8)
    parser.add_argument("-n", "--chunk_id", default="1", help="prefix for the output files")
    parser.add_argument("--output_dir", type=Path, default=Path("output"))
    parser.add_argument(
        "--tmp_dir",
        type=Path,
        help="parent for this run's scratch directory (default: $TMPDIR)",
    )
    parser.add_argument("--worker_memory_mb", type=int, default=512)
    parser.add_argument("--min_edge_jaccard", type=parse_threshold, default=Fraction(1, 10))
    parser.add_argument("--min_edge_containment", type=parse_threshold, default=Fraction(1, 2))
    parser.add_argument("--clan_jaccard", type=parse_threshold, default=Fraction(1, 2))
    return parser.parse_args(args)


def output_paths(options: argparse.Namespace) -> dict[str, Path]:
    root, chunk = options.output_dir, options.chunk_id
    return {
        "edges": root / f"{chunk}_clan_edges.csv",
        "clans": root / f"{chunk}_clans.csv",
        "summary": root / f"{chunk}_clans_summary.json",
        "log": root / f"{chunk}_clans.log",
    }


def validate_options(options: argparse.Namespace) -> Path:
    """Check flags and paths before anything is created. Returns the scratch parent."""
    if CHUNK_PATTERN.fullmatch(options.chunk_id) is None:
        raise InputError("chunk_id must match [A-Za-z0-9._-]+")
    if options.cpus < 1:
        raise InputError("cpus must be at least 1")
    if options.worker_memory_mb < 16:
        raise InputError("worker_memory_mb must be at least 16")
    inputs = [Path(path) for path in options.members]
    if options.fasta_file:
        inputs.append(Path(options.fasta_file))
    for path in inputs:
        if not path.is_file():
            raise InputError(f"{path} must be an existing file")
    outputs = list(output_paths(options).values())
    resolved = [path.resolve() for path in outputs]
    if len(set(resolved)) != len(resolved):
        raise InputError("output paths alias one another")
    identities = {(p.stat().st_dev, p.stat().st_ino) for p in inputs}
    resolved_inputs = {path.resolve() for path in inputs}
    seen: set[tuple[int, int]] = set()
    for path, target in zip(outputs, resolved, strict=True):
        if target in resolved_inputs:
            raise InputError(f"output path {path} overlaps an input file")
        if path.exists():
            identity = (path.stat().st_dev, path.stat().st_ino)
            if identity in identities:
                raise InputError(f"output path {path} overlaps an input file")
            if identity in seen:
                raise InputError(f"output path {path} aliases another output")
            seen.add(identity)
    parent = Path(options.tmp_dir or tempfile.gettempdir())
    if not parent.is_dir():
        raise InputError(f"tmp_dir {parent} must be an existing directory")
    root = options.output_dir.resolve()
    if parent.resolve() == root or parent.resolve().is_relative_to(root):
        raise InputError("tmp_dir must not be the output_dir or inside it")
    return parent


@contextlib.contextmanager
def worker_pool(cpus: int, initializer: Any = None, initargs: tuple = ()) -> Iterator[Any]:
    """Unordered map over a process pool, or in-process for one cpu.

    `ProcessPoolExecutor`, not `multiprocessing.Pool`: a worker killed by a signal (OOM,
    SIGBUS) raises `BrokenProcessPool` instead of losing its task and hanging the parent.
    """
    if cpus == 1:
        if initializer:
            initializer(*initargs)
        yield lambda function, items: map(function, items)
        return
    context = multiprocessing.get_context("forkserver")
    pool = ProcessPoolExecutor(cpus, context, initializer=initializer, initargs=initargs)

    def imap(function: Any, items: Iterable[Any]) -> Iterator[Any]:
        pending = {pool.submit(function, item) for item in items}
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()

    try:
        yield imap
    finally:
        pool.shutdown(cancel_futures=True)


def star(arguments: tuple) -> Any:
    function, *rest = arguments
    return function(*rest)


def find_missing_name(paths: Iterable[str], keys: set[tuple[int, int]]) -> str:
    """Error path only: the first whole-protein member whose key has no FASTA record."""
    for path in paths:
        whole_file = Task(0, path, 0, -1 if is_gzipped(Path(path)) else 2**63)
        for _, raw in task_lines(whole_file):
            member = raw.rstrip(b"\r\n").partition(b"\t")[2]
            if split_member(member) is None:
                digest = protein_key(member)
                key = (int.from_bytes(digest[:8], "little"), int.from_bytes(digest[8:], "little"))
                if key in keys:
                    return member.decode(errors="replace")
    return "?"


@dataclass
class Totals:
    rows: int = 0
    bare_rows: int = 0
    member_runs: list[tuple[str, int, int]] = field(default_factory=list)
    member_counts: np.ndarray = field(default_factory=lambda: np.zeros(2**BUCKET_BITS, np.int64))
    length_runs: list[tuple[str, int]] = field(default_factory=list)
    length_counts: np.ndarray = field(default_factory=lambda: np.zeros(2**BUCKET_BITS, np.int64))
    fasta_records: int = 0
    max_length: int = 0


def run_phases(
    options: argparse.Namespace, scratch: Path, logger: logging.Logger, sizes: Sizes
) -> dict[str, Any]:
    """Phases A and B: parse, validate and accumulate. Writes only under `scratch`."""
    (scratch / "members").mkdir()
    (scratch / "fasta").mkdir()
    started = time.monotonic()
    totals = Totals()
    tasks = plan_tasks(options.members, spread(sizes.range_bytes, options.members, options.cpus))
    results: list[MembersResult] = []
    with worker_pool(options.cpus) as imap:
        work = [(parse_members_task, task, str(scratch), sizes.run_rows) for task in tasks]
        results.extend(imap(star, work))
    results.sort(key=lambda result: result.number)
    names: set[bytes] = set()
    for result in results:
        totals.rows += result.rows
        totals.bare_rows += result.bare_rows
        totals.member_counts += result.counts
        names.update(filter(None, Path(result.names_path).read_bytes().split(b"\n")))
    families = sorted(names)
    index = {name: number for number, name in enumerate(families)}
    del names
    remap_path = scratch / "remap.bin"
    offset = 0
    with remap_path.open("wb") as remap:
        for result in results:
            local = list(filter(None, Path(result.names_path).read_bytes().split(b"\n")))
            np.array([index[name] for name in local], dtype="<u4").tofile(remap)
            totals.member_runs.extend((path, rows, offset) for path, rows in result.runs)
            offset += len(local)
    logger.info(
        "members rows=%d whole_protein=%d families=%d tasks=%d runs=%d elapsed=%.1fs",
        totals.rows,
        totals.bare_rows,
        len(families),
        len(tasks),
        len(totals.member_runs),
        time.monotonic() - started,
    )

    if totals.bare_rows:
        if not options.fasta_file:
            raise InputError(
                f"{totals.bare_rows} members are whole proteins (bare names); "
                "pass the database FASTA with --fasta_file for their lengths"
            )
        with worker_pool(options.cpus) as imap:
            work = [
                (parse_fasta_task, task, str(scratch), sizes.run_rows)
                for task in plan_tasks(
                    [options.fasta_file],
                    spread(sizes.range_bytes, [options.fasta_file], options.cpus),
                )
            ]
            for lengths in imap(star, work):
                totals.length_runs.extend(lengths.runs)
                totals.length_counts += lengths.counts
                totals.fasta_records += lengths.records
                totals.max_length = max(totals.max_length, lengths.max_length)
        logger.info(
            "fasta records=%d max_length=%d elapsed=%.1fs",
            totals.fasta_records,
            totals.max_length,
            time.monotonic() - started,
        )

    # Every n(F) and S(A, B) is at most `rows`. Residues are bounded after Phase B, exactly.
    check_exact_bound("member rows", totals.rows)

    count = len(families)
    pairs = count * (count - 1) // 2
    shared_dtype = "<u4" if totals.rows < 2**32 else "<u8"
    shared_path = overlap_path = None
    if pairs:
        triangle_bytes = pairs * (np.dtype(shared_dtype).itemsize + 8)
        # Compared with what is still to be written, not the total: the runs already exist.
        free = shutil.disk_usage(scratch).free
        logger.info("triangle pairs=%d bytes=%d scratch_free=%d", pairs, triangle_bytes, free)
        reserve = getattr(os, "posix_fallocate", None)
        if reserve is None and free < triangle_bytes:
            raise ResourceError(f"scratch: {free} bytes free, triangle needs {triangle_bytes}")
        shared_path, overlap_path = str(scratch / "shared.bin"), str(scratch / "overlap.bin")
        for path, dtype in ((shared_path, shared_dtype), (overlap_path, "<u8")):
            size = pairs * np.dtype(dtype).itemsize
            with open(path, "wb") as handle:
                # A sparse file gets blocks only as workers write through the map, and a full
                # disk then kills them with SIGBUS. Reserving fails here with ENOSPC instead.
                if reserve:
                    reserve(handle.fileno(), 0, size)
                else:
                    # shortcut: sparse where posix_fallocate is missing (macOS); checked above.
                    handle.truncate(size)

    units = plan_units(
        totals.member_runs,
        totals.length_runs,
        totals.member_counts,
        totals.length_counts,
        sizes.unit_cost,
        # About four units per worker, so one large unit does not serialise the phase.
        int(
            totals.member_counts.sum() * MEMBER_ROW_COST
            + totals.length_counts.sum() * LENGTH_ROW_COST
        )
        // (4 * options.cpus),
    )
    context = multiprocessing.get_context("forkserver")
    locks = (
        [context.Lock() for _ in range(LOCK_STRIPES)]
        if options.cpus > 1
        else [contextlib.nullcontext() for _ in range(LOCK_STRIPES)]
    )
    state = UnitState(
        totals.member_runs,
        totals.length_runs,
        str(remap_path),
        count,
        shared_path,
        shared_dtype,
        overlap_path,
        sizes.emission_batch,
        locks,
    )
    proteins_per_family = np.zeros(count, np.int64)
    residues_per_family = np.zeros(count, np.int64)
    proteins = pair_emissions = overlap_emissions = max_k = 0
    try:
        with worker_pool(options.cpus, init_unit_worker, (state,)) as imap:
            for result in imap(process_unit, units):
                proteins_per_family += result.proteins_per_family
                residues_per_family += result.residues_per_family
                proteins += result.proteins
                pair_emissions += result.pair_emissions
                overlap_emissions += result.overlap_emissions
                max_k = max(max_k, result.max_families_per_protein)
    except MissingLengths as missing:
        name = find_missing_name(options.members, set(missing.keys))
        raise InputError(
            f"whole-protein members have no record in the FASTA ({missing.count} rows in "
            f"one hash range), for example {name!r}"
        ) from None
    finally:
        _UNIT.clear()
    logger.info(
        "pairs units=%d proteins=%d pair_emissions=%d overlap_emissions=%d "
        "max_families_per_protein=%d elapsed=%.1fs",
        len(units),
        proteins,
        pair_emissions,
        overlap_emissions,
        max_k,
        time.monotonic() - started,
    )
    # I(A, B) <= min |R|, so this bounds every residue sum and threshold cross-product. It
    # runs before Phase C, so a rejection still leaves the output tree alone.
    check_exact_bound("family residues", int(residues_per_family.max()) if count else 0)
    return {
        "families": [name.decode() for name in families],
        "proteins_per_family": proteins_per_family,
        "residues_per_family": residues_per_family,
        "shared": np.memmap(shared_path, dtype=shared_dtype, mode="r") if pairs else None,
        "overlap": np.memmap(overlap_path, dtype="<u8", mode="r") if pairs else None,
        "counts": {
            "member_rows": totals.rows,
            "whole_protein_member_rows": totals.bare_rows,
            "families": count,
            "proteins": proteins,
            "pair_emissions": pair_emissions,
            "max_families_per_protein": max_k,
        },
    }


def check_exact_bound(name: str, bound: int) -> None:
    """Reject a count whose threshold cross-products could overflow int64."""
    if bound * 2 * THRESHOLD_SCALE >= INT64_LIMIT:
        raise InputError(f"{name} ({bound}) exceed the supported bound for exact arithmetic")


def exact_at_least(
    numerator: np.ndarray, denominator: np.ndarray, threshold: Fraction
) -> np.ndarray:
    """`numerator / denominator >= threshold`, exactly, in int64 (bounds checked earlier)."""
    return numerator * threshold.denominator >= threshold.numerator * denominator


def write_outputs(
    options: argparse.Namespace, paths: dict[str, Path], data: dict[str, Any]
) -> dict[str, Any]:
    """Phase C: scan the triangle, write edges and clans. Returns the summary payload."""
    names = data["families"]
    count = len(names)
    n, residues = data["proteins_per_family"], data["residues_per_family"]
    shared, overlap = data["shared"], data["overlap"]
    row_start = np.array([a * count - a * (a + 1) // 2 for a in range(count + 1)], np.int64)
    clan_first: list[np.ndarray] = []
    clan_second: list[np.ndarray] = []
    clan_weight: list[np.ndarray] = []
    sharing = edges = near_redundant = 0
    block = 1 << 22
    with paths["edges"].open("w", encoding="utf-8") as handle:
        handle.write(EDGES_HEADER)
        a0 = 0
        while shared is not None and a0 < count - 1:
            a1 = a0 + 1
            while a1 < count - 1 and row_start[a1 + 1] - row_start[a0] <= block:
                a1 += 1
            low, high = int(row_start[a0]), int(row_start[a1])
            present = np.flatnonzero(np.asarray(shared[low:high]))
            a0 = a1
            if not len(present):
                continue
            index = present + low
            a = np.searchsorted(row_start, index, side="right") - 1
            b = index - row_start[a] + a + 1
            s = np.asarray(shared[index]).astype(np.int64)
            i = np.asarray(overlap[index]).astype(np.int64)
            union = residues[a] + residues[b] - i
            smaller = np.minimum(residues[a], residues[b])
            protein_union = n[a] + n[b] - s
            sharing += len(index)
            near_redundant += int(exact_at_least(i, union, NEAR_REDUNDANT).sum())
            keep = exact_at_least(i, union, options.min_edge_jaccard) | exact_at_least(
                i, smaller, options.min_edge_containment
            )
            clan = exact_at_least(i, union, options.clan_jaccard)
            clan_first.append(a[clan])
            clan_second.append(b[clan])
            clan_weight.append(i[clan] / union[clan])
            for row in np.flatnonzero(keep).tolist():
                handle.write(
                    f"{names[a[row]]},{names[b[row]]},{s[row]},{s[row] / protein_union[row]:.6f},"
                    f"{i[row] / union[row]:.6f},{i[row] / smaller[row]:.6f}\n"
                )
            edges += int(keep.sum())

    empty_int, empty_float = np.zeros(0, np.int64), np.zeros(0, np.float64)
    clans = build_clans(
        count,
        np.concatenate(clan_first) if clan_first else empty_int,
        np.concatenate(clan_second) if clan_second else empty_int,
        np.concatenate(clan_weight) if clan_weight else empty_float,
    )
    with paths["clans"].open("w", encoding="utf-8") as handle:
        handle.write(CLANS_HEADER)
        for number, (representative, members) in enumerate(clans, start=1):
            clan_id = f"{options.chunk_id}_clan_{number}"
            for member in members:
                handle.write(f"{clan_id},{names[member]},{member == representative}\n")

    in_clans = sum(len(members) for _, members in clans)
    sizes = Counter(len(members) for _, members in clans)
    counts = data["counts"]
    counts.update(
        {
            "pairs_sharing_proteins": sharing,
            "edges": edges,
            "near_redundant_pairs": near_redundant,
            "clans": len(clans),
            "families_in_clans": in_clans,
            "singleton_families": count - in_clans,
        }
    )
    return {
        "command": "compute_clans",
        "version": __version__,
        "chunk_id": options.chunk_id,
        "parameters": {name: float(getattr(options, name)) for name in SCIENTIFIC_PARAMETERS},
        "counts": counts,
        "clan_size_histogram": {str(size): sizes[size] for size in sorted(sizes)},
    }


def main(args: SequenceCollection[str] | None = None) -> None:
    options = parse_args(args)
    try:
        parent = validate_options(options)
    except InputError as error:
        print(f"mgnifam compute_clans: error: {error}", file=sys.stderr)
        raise SystemExit(2) from None
    paths = output_paths(options)
    scratch = Path(tempfile.mkdtemp(prefix="mgnifam_clans_", dir=parent))
    logger = configure_logger(scratch / "clans.log")
    started = time.monotonic()
    log_moved = False
    try:
        logger.info(
            "start members=%d cpus=%d worker_memory_mb=%d scratch=%s",
            len(options.members),
            options.cpus,
            options.worker_memory_mb,
            scratch,
        )
        data = run_phases(options, scratch, logger, sizes_for(options.worker_memory_mb))
        # Phase C. Nothing before this point touched the output tree.
        options.output_dir.mkdir(parents=True, exist_ok=True)
        paths["summary"].unlink(missing_ok=True)
        move_log(logger, scratch / "clans.log", paths["log"])
        log_moved = True
        payload = write_outputs(options, paths, data)
        logger.info("DONE. elapsed=%.1fs", time.monotonic() - started)
        write_stats(paths["summary"], payload)
    except InputError as error:
        logger.error("%s", error)
        print(f"mgnifam compute_clans: error: {error}", file=sys.stderr)
        raise SystemExit(2) from None
    except (ChunkCorrupted, ResourceError, BrokenProcessPool, OSError) as error:
        logger.exception("compute_clans failed; outputs must not be consumed")
        # A failed run, unlike a rejected input, must not leave a previous summary vouching
        # for it, and its traceback must outlive scratch.
        with contextlib.suppress(OSError):
            options.output_dir.mkdir(parents=True, exist_ok=True)
            paths["summary"].unlink(missing_ok=True)
            if not log_moved:
                move_log(logger, scratch / "clans.log", paths["log"])
        print(f"mgnifam compute_clans: {error}", file=sys.stderr)
        raise SystemExit(1) from None
    finally:
        for handler in logger.handlers:
            handler.close()
        logger.handlers.clear()
        shutil.rmtree(scratch, ignore_errors=True)


def move_log(logger: logging.Logger, scratch_log: Path, destination: Path) -> None:
    """Carry the run log from scratch into the output tree once Phase C begins."""
    formatter = logger.handlers[0].formatter
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    shutil.copyfile(scratch_log, destination)
    handler = logging.FileHandler(destination, mode="a", encoding="utf-8")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


if __name__ == "__main__":
    main()
