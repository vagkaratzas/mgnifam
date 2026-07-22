"""Verification suite for the family-generation port.

The end-to-end tests run against the real 50,000-sequence fixtures rather than toy
data, because several of the properties under test -- the reporting threshold, hit
ordering, envelope clipping -- only appear on a database large enough to produce
marginal hits.

Three tests exist precisely because an end-to-end run cannot observe what they check:
`test_unreported_hit_is_not_recruited_and_prefetch_matches`,
`test_index_mismatch_guard_and_optimized_python`, and
`test_extracted_records_do_not_retain_pyhmmer_results`. Each guards an invariant whose
violation leaves the final artifacts looking correct. Read their docstrings before
deleting them as redundant.
"""

import argparse
import contextlib
import gc
import gzip
import io
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyhmmer
import pytest

from mgnifam import __version__, cli
from mgnifam import generate_families as gf


def cli_args(
    clusters: Path,
    fasta: Path,
    *,
    cpus: int = 1,
    chunk: str = "chunk",
    batch_size: int = 0,
    fasta_index: Path | None = None,
    prefetch: bool = False,
) -> list[str]:
    arguments = [
        "--clusters_chunk",
        str(clusters),
        "--fasta_file",
        str(fasta),
        "--cpus",
        str(cpus),
        "--chunk_num",
        chunk,
        "--discard_min_rep_length",
        "100",
        "--discard_max_rep_length",
        "2000",
        "--discard_min_starting_membership",
        "0.9",
        "--max_seq_identity",
        "0.8",
        "--max_seed_seqs",
        "2000",
        "--max_gap_occupancy",
        "0.5",
        "--recruit_evalue_cutoff",
        "0.001",
        "--recruit_hit_length_percentage",
        "0.9",
        "--batch_size",
        str(batch_size),
    ]
    if fasta_index is not None:
        arguments.extend(["--fasta_index", str(fasta_index)])
    if prefetch:
        arguments.append("--prefetch_targets")
    return arguments


def run_pipeline(directory: Path, arguments: list[str]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    with contextlib.chdir(directory):
        gf.main(arguments)
    return directory / "output"


def scientific_artifacts(directory: Path) -> dict[str, bytes]:
    artifacts = {}
    for path in sorted(directory.glob("**/*")):
        # Logs carry timestamps and the SSI index is an input cache: neither is a result.
        if not path.is_file() or path.suffix in (".log", ".ssi"):
            continue
        artifacts[str(path.relative_to(directory))] = path.read_bytes()
    return artifacts


@pytest.fixture(scope="session")
def shared_index(tmp_path_factory: pytest.TempPathFactory, small_fasta: Path) -> Path:
    index = tmp_path_factory.mktemp("ssi") / "small.ssi"
    gf.build_ssi_index(small_fasta, index)
    return index


@pytest.fixture(scope="session")
def baseline_output(
    tmp_path_factory: pytest.TempPathFactory,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> Path:
    output = tmp_path_factory.mktemp("baseline")
    return run_pipeline(
        output,
        cli_args(
            fixture_directory / "clustering.tsv",
            small_fasta,
            fasta_index=shared_index,
        ),
    )


@pytest.fixture(scope="session")
def v2_output(
    tmp_path_factory: pytest.TempPathFactory,
    fixture_directory: Path,
    v2_fasta: Path,
) -> Path:
    return run_pipeline(
        tmp_path_factory.mktemp("v2-full-tsv"),
        cli_args(
            fixture_directory / "mgnifams_v2.tsv",
            v2_fasta,
            cpus=2,
            chunk="v2",
        ),
    )


@pytest.mark.filterwarnings(
    "ignore:Exception ignored in.*SSIWriter:pytest.PytestUnraisableExceptionWarning"
)
def test_build_ssi_index_round_trip_and_errors(tmp_path: Path, small_fasta: Path) -> None:
    assert isinstance(pyhmmer.easel.SequenceFile._FORMATS["fasta"], int)
    index = tmp_path / "small.ssi"
    gf.build_ssi_index(small_fasta, index)

    with pyhmmer.easel.SequenceFile(small_fasta, digital=False) as source:
        expected = {sequence.name: sequence.sequence for sequence in source}
    with (
        pyhmmer.easel.SSIReader(index) as reader,
        pyhmmer.easel.SequenceFile(small_fasta, digital=False, index=reader) as indexed,
    ):
        observed = {
            name: gf.fetch_indexed_sequence(indexed.indexed, name).sequence for name in expected
        }
        assert observed == expected
        with pytest.raises(KeyError):
            gf.fetch_indexed_sequence(indexed.indexed, "missing")

    duplicate_fasta = tmp_path / "duplicate.fa"
    duplicate_fasta.write_text(">same\nAAAA\n>same\nCCCC\n")
    before = set(tmp_path.iterdir())
    with pytest.raises(gf.DuplicateSequenceName, match="unique"):
        gf.build_ssi_index(duplicate_fasta, tmp_path / "duplicate.ssi")
    assert set(tmp_path.iterdir()) == before


def test_index_mismatch_guard_and_optimized_python(tmp_path: Path) -> None:
    wrong = pyhmmer.easel.TextSequence(name="wrong", sequence="AAAA")
    with pytest.raises(gf.IndexMismatchError):
        gf.fetch_indexed_sequence({"wanted": wrong}, "wanted")

    code = """
from pyhmmer.easel import TextSequence
from mgnifam.generate_families import IndexMismatchError, fetch_indexed_sequence
try:
    fetch_indexed_sequence({'wanted': TextSequence(name='wrong', sequence='AAAA')}, 'wanted')
except IndexMismatchError:
    raise SystemExit(0)
raise SystemExit(1)
"""
    result = subprocess.run(
        [sys.executable, "-O", "-c", code],
        cwd=tmp_path,
        check=False,
    )
    assert result.returncode == 0


def test_soft_masked_fasta_is_normalised_at_the_fetch_boundary(tmp_path: Path) -> None:
    """A lower-case (soft-masked) database must not break residue location.

    Lower case in a HMMER alignment marks an insert-state residue, so every row is
    upper-cased before its residues are located in the parent record. A soft-masked
    FASTA therefore used to fail `parse_protein_name`'s `str.find` and take the whole
    chunk down at emit time -- after the searches had already been paid for.
    """
    masked = tmp_path / "masked.fa"
    masked.write_text(">prot_101_140\nmktaylaagivgqqqqq\n")
    index = tmp_path / "masked.ssi"
    gf.build_ssi_index(masked, index)

    with (
        pyhmmer.easel.SSIReader(index) as reader,
        pyhmmer.easel.SequenceFile(masked, digital=False, index=reader) as handle,
    ):
        sequences = gf.IndexedSequences(handle)
        assert sequences.get("prot_101_140") == gf.Sequence("prot_101_140", "MKTAYLAAGIVGQQQQQ")
        # The alignment row is upper case with an insert column; it must still resolve.
        assert gf.parse_protein_name("prot_101_140", "MKTAY-laa", sequences) == "prot/101-108"


def text_msa(names: list[str], sequences: list[str], reference: str) -> pyhmmer.easel.TextMSA:
    msa = pyhmmer.easel.TextMSA(
        sequences=[
            pyhmmer.easel.TextSequence(name=name, sequence=sequence)
            for name, sequence in zip(names, sequences, strict=True)
        ]
    )
    msa.reference = reference
    return msa


def test_clip_env_ends_preserves_interior_insert_columns() -> None:
    msa = text_msa(["a", "b"], ["ABCDEFGHIJKL", "ABCDEFGHIJKL"], "...xx.xxx...")
    clipped = gf.clip_env_ends(msa)
    assert clipped.reference == "xx.xxx"
    assert list(clipped.alignment) == ["DEFGHI", "DEFGHI"]


class FakeSequences:
    def __init__(self, sequences: dict[str, str]) -> None:
        self.sequences = sequences

    def get(self, name: str, *, missing_ok: bool = False) -> gf.Sequence | None:
        if name not in self.sequences:
            if missing_ok:
                return None
            raise KeyError(name)
        return gf.Sequence(name, self.sequences[name])


def test_parse_protein_name_resolves_repeats_within_their_envelope() -> None:
    """Two identical repeat domains of one protein must renumber to different regions.

    The legacy script located a row's residues in the whole record, so `str.find` returned
    the first copy for both domains, they collided on one name and one was dropped as a
    duplicate. Searching within the row's own envelope keeps them apart.
    """
    repeat = "MKVLAAGIVG"
    store = FakeSequences({"prot_101_140": f"{repeat}QQQQQ{repeat}QQQQQ"})

    first = gf.parse_protein_name("prot_101_140/1_10", repeat, store)
    second = gf.parse_protein_name("prot_101_140/16_25", repeat, store)

    assert (first, second) == ("prot/101-110", "prot/116-125")

    # Gap characters are stripped, and a row shorter than its record still gets coordinates
    # -- legacy truncated this name to the bare accession to fit the old name column.
    assert gf.parse_protein_name("prot_101_140", f"-{repeat[1:]}...QQQQQ", store) == "prot/102-115"

    # A row spanning the whole of an unsliced record keeps the bare accession, which
    # `family_metadata` records as region "-".
    whole = FakeSequences({"prot": repeat})
    assert gf.parse_protein_name("prot", repeat, whole) == "prot"

    with pytest.raises(ValueError, match="not in its envelope"):
        gf.parse_protein_name("prot_101_140/1_10", "WWWWWWWWWW", store)


def test_seed_membership_counts_distinct_proteins_on_both_sides() -> None:
    """Membership is a ratio of distinct proteins, so it can never exceed 1.

    Dividing by the raw row count instead let a cluster TSV that repeats a member report
    less than full membership for a family that had in fact kept every one of them, and
    the family was discarded as "few seed sequences remained" at the 0.9 default.
    """
    assert gf.check_seed_membership(["a", "b", "c"], ["a", "b", "c"]) == 1.0
    assert gf.check_seed_membership(["a", "a", "b", "c"], ["a", "b", "c"]) == 1.0

    # Envelope suffixes collapse to the parent protein on both sides.
    assert gf.check_seed_membership(["p/1_9", "p/20_30"], ["p/1_9"]) == 1.0

    # Genuine loss still measures as loss.
    assert gf.check_seed_membership(["a", "b", "c", "d"], ["a", "b"]) == 0.5


def test_filter_hits_exit_filter_and_native_order() -> None:
    records = [("b", 10, 2, 4), ("a", 10, 1, 10), ("b", 10, 5, 10)]
    store = FakeSequences({"a": "A" * 10, "b": "B" * 10})
    regular = gf.filter_hits(records, 10, False, 0.9, store)
    exiting = gf.filter_hits(records, 10, True, 0.9, store)
    assert [sequence.id for sequence in regular] == ["a"]
    assert [sequence.id for sequence in exiting] == ["b/2_4", "a", "b/5_10"]


def family_options(**changes: object) -> argparse.Namespace:
    options = dict(
        cpus=1,
        recruit_hit_length_percentage=0.9,
        discard_min_rep_length=1,
        discard_max_rep_length=1000,
        max_seq_identity=0.8,
        max_seed_seqs=2000,
        max_gap_occupancy=0.5,
        discard_min_starting_membership=0.9,
    )
    options.update(changes)
    return argparse.Namespace(**options)


def test_small_msa_is_discarded_before_another_build(monkeypatch: pytest.MonkeyPatch) -> None:
    aligned = text_msa(["a", "b", "c"], ["AAAA", "AAAT", "AATT"], "xxxx")
    tiny = text_msa(["a", "b"], ["AAAA", "AAAT"], "xxxx")
    monkeypatch.setattr(gf, "run_hmmalign", lambda *_args, **_kwargs: aligned)
    monkeypatch.setattr(gf, "run_pytrimal_reps", lambda *_args, **_kwargs: tiny)
    monkeypatch.setattr(
        gf,
        "run_hmmbuild",
        lambda *_args, **_kwargs: pytest.fail("hmmbuild must not be called"),
    )
    family = gf.Family("a", ["a", "b", "c"], hmm=object(), records=[("a", 4, 1, 4)])
    family.advance(family_options(), FakeSequences({"a": "AAAA"}), 1)
    assert family.state is gf.FamilyState.DISCARDED
    assert family.discard_reason == "too few sequences after redundancy filtering"
    assert family.discard_value == 2


def test_final_round_leaves_the_searched_seed_in_place(monkeypatch: pytest.MonkeyPatch) -> None:
    """The last round must not re-align and trim: its output would feed no build.

    Doing it anyway is what left the legacy exit path exporting a hand-built HMM one
    generation ahead of the model that produced the family's hits and full MSA.
    """
    seed = text_msa(["a", "b", "c"], ["AAAA", "AAAT", "AATT"], "xxxx").digitize(gf.ALPHABET)
    monkeypatch.setattr(
        gf,
        "run_hmmalign",
        lambda *_args, **_kwargs: pytest.fail("the final round must not re-align"),
    )
    family = gf.Family(
        "a", ["a", "b", "c"], seed_msa=seed, hmm=object(), records=[("a", 4, 1, 4)], qlen=4
    )
    family.advance(family_options(), FakeSequences({"a": "AAAA"}), gf.MAX_ROUNDS)
    assert family.state is gf.FamilyState.RUNNING
    assert family.seed_msa is seed


def test_converged_discard_is_not_recorded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    store = FakeSequences({"a": "AAAA", "b": "AAAT"})
    seed = text_msa(["a", "b"], ["AAAA", "AAAT"], "xxxx").digitize(gf.ALPHABET)
    family = gf.Family(
        "a",
        ["a", "b"],
        seed_msa=seed,
        hmm=object(),
        records=[("a", 4, 1, 4)],
        qlen=4,
        total_checked_sequences={"a"},
    )
    options = family_options(discard_min_starting_membership=1.0)
    family.advance(options, store, 1)
    assert family.state is gf.FamilyState.CONVERGED
    family.finish(options, store)
    assert family.state is gf.FamilyState.DISCARDED
    assert family.ever_converged

    for directory in gf.FAMILY_DIRECTORIES:
        (tmp_path / directory).mkdir()
    writers = gf.Writers(
        root=tmp_path,
        indexed=store,
        refined_families=io.StringIO(),
        discarded_clusters=io.StringIO(),
        successful_clusters=io.StringIO(),
        converged_families=io.StringIO(),
        family_metadata=io.StringIO(),
        family_representatives=io.StringIO(),
    )
    success_count = gf.emit_family(family, 0, "chunk", writers)
    assert success_count == 0
    assert writers.converged_families.getvalue() == ""

    successful = gf.Family(
        "b",
        ["a", "b"],
        state=gf.FamilyState.SUCCESSFUL,
        seed_msa=seed,
        full_msa=text_msa(["a", "b"], ["AAAA", "AAAT"], "xxxx"),
        full_msa_num_seqs=2,
        ever_converged=True,
    )
    success_count = gf.emit_family(successful, success_count, "chunk", writers)
    assert success_count == 1
    assert successful.family_id == 1
    assert writers.converged_families.getvalue() == "1\n"


def test_cpus_and_sanity_anchors(
    tmp_path: Path,
    baseline_output: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    cpus_four = run_pipeline(
        tmp_path / "cpus-four",
        cli_args(
            fixture_directory / "clustering.tsv",
            small_fasta,
            cpus=4,
            fasta_index=shared_index,
        ),
    )
    assert scientific_artifacts(cpus_four) == scientific_artifacts(baseline_output)
    assert (baseline_output / "chunk_successful.txt").read_text().splitlines() == [
        "4706047775",
        "1622851798_832_939",
        "4497037939_1_144",
    ]
    metadata = (baseline_output / "chunk_metadata.csv").read_text().splitlines()
    assert [line.split(",")[2].strip('"') for line in metadata] == [
        "782510898",
        "5761513631",
        "1446399400",
    ]
    # Families 2 and 3 recruit nothing new on their second pass. Family 3 only started
    # converging once clip_ends stopped discarding a column of its model (CHANGELOG 1.0.0).
    assert (baseline_output / "chunk_converged.txt").read_text() == "2\n3\n"
    # One representative per successful family, and ids are contiguous from 1.
    assert [line.split(",")[0] for line in metadata] == ["1", "2", "3"]


def test_batch_size_invariance(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    outputs = []
    for batch_size in (1, 64):
        output = run_pipeline(
            tmp_path / f"batch-{batch_size}",
            cli_args(
                fixture_directory / "clustering.tsv",
                small_fasta,
                batch_size=batch_size,
                fasta_index=shared_index,
            ),
        )
        outputs.append(output)
    assert scientific_artifacts(outputs[0]) == scientific_artifacts(outputs[1])
    for output in outputs:
        mapping = [
            (line.split(",")[2], line.split(",", 1)[0])
            for line in (output / "chunk_metadata.csv").read_text().splitlines()
        ]
        assert mapping == [('"782510898"', "1"), ('"5761513631"', "2"), ('"1446399400"', "3")]


def test_determinism_across_hash_seeds(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    outputs = []
    arguments = cli_args(
        fixture_directory / "clustering.tsv", small_fasta, fasta_index=shared_index
    )
    for hash_seed in ("1", "987654"):
        run_directory = tmp_path / f"hash-{hash_seed}"
        run_directory.mkdir()
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = hash_seed
        subprocess.run(
            [sys.executable, "-m", "mgnifam.generate_families", *arguments],
            cwd=run_directory,
            env=environment,
            check=True,
        )
        outputs.append(run_directory / "output")
    assert scientific_artifacts(outputs[0]) == scientific_artifacts(outputs[1])
    for hmm_path in (outputs[0] / "hmm").glob("*.hmm.gz"):
        contents = gzip.decompress(hmm_path.read_bytes())
        assert not any(line.startswith(b"DATE ") for line in contents.splitlines())
        assert not any(line.startswith(b"COM ") for line in contents.splitlines())
    gzip_header = next((outputs[0] / "hmm").glob("*.hmm.gz")).read_bytes()
    assert not gzip_header[3] & 0x08


def test_rerun_removes_stale_family_artifacts(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    run_directory = tmp_path / "rerun"
    arguments = cli_args(
        fixture_directory / "clustering.tsv", small_fasta, fasta_index=shared_index
    )
    output = run_pipeline(run_directory, arguments)
    reduced_clusters = tmp_path / "reduced.tsv"
    reduced_clusters.write_text(
        "\n".join((fixture_directory / "clustering.tsv").read_text().splitlines()[:9]) + "\n"
    )
    run_pipeline(run_directory, cli_args(reduced_clusters, small_fasta, fasta_index=shared_index))
    for directory in gf.FAMILY_DIRECTORIES:
        assert all(
            "chunk_2." not in path.name and "chunk_3." not in path.name
            for path in (output / directory).iterdir()
        )


def test_one_failing_family_is_discarded_and_the_chunk_survives(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A family that raises must cost its own result, not the whole chunk's.

    A chunk is hours of searching and nothing is written until a wave ends, so an
    exception from one family used to destroy every family beside it. `emit_family` is
    the stage under test because it is the one that writes: the failure has to leave no
    half-written trail, so the representative must reach `discarded` and never
    `successful`, and the ids of the families that follow must stay contiguous.
    """
    clean = run_pipeline(
        tmp_path / "clean",
        cli_args(fixture_directory / "clustering.tsv", small_fasta, fasta_index=shared_index),
    )
    survivors = (clean / "chunk_successful.txt").read_text().splitlines()
    assert len(survivors) > 1, "fixture must produce several families for this to mean anything"
    doomed = survivors[0]

    original = gf.renumber_msa
    failed_once = False

    def explode(msa, family_name, indexed):  # type: ignore[no-untyped-def]
        # Keyed on the first call, not on the name: a failed family releases its id, so
        # the next family to succeed is renamed into it and a name-keyed trap cascades.
        nonlocal failed_once
        if not failed_once:
            failed_once = True
            raise RuntimeError("synthetic failure, comma included")
        return original(msa, family_name, indexed)

    monkeypatch.setattr(gf, "renumber_msa", explode)
    output = run_pipeline(
        tmp_path / "guarded",
        cli_args(fixture_directory / "clustering.tsv", small_fasta, fasta_index=shared_index),
    )

    discarded = (output / "chunk_discarded.csv").read_text().splitlines()
    failed_rows = [row for row in discarded if "internal error" in row]
    assert failed_rows == [f"{doomed},internal error during artifact writing,0.0"]
    # A comma in the stage text would have split the CSV.
    assert all(len(row.split(",")) == 3 for row in discarded)

    remaining = (output / "chunk_successful.txt").read_text().splitlines()
    assert doomed not in remaining
    assert remaining == survivors[1:]

    # The failure wrote nothing, so no artifact carries the id it would have taken, and
    # the families after it close the gap rather than inheriting it.
    ids = sorted(int(path.name.split("_")[1].split(".")[0]) for path in (output / "hmm").iterdir())
    assert ids == list(range(1, len(remaining) + 1))
    assert "synthetic failure" in (output / "chunk.log").read_text()


def test_family_guard_contains_only_exceptions() -> None:
    """The guard converts a failure into a discard, but must not swallow an interrupt."""
    logger = logging.getLogger("mgnifam.test.guard")
    family = gf.Family("rep", ["a"])
    with gf.family_guard(family, logger, "round 2"):
        raise ValueError("boom")
    assert family.state is gf.FamilyState.DISCARDED
    assert family.discard_reason == "internal error during round 2"

    unaffected = gf.Family("rep2", ["a"])
    with pytest.raises(KeyboardInterrupt), gf.family_guard(unaffected, logger, "round 2"):
        raise KeyboardInterrupt
    assert unaffected.state is gf.FamilyState.RUNNING


def _raw_signature(top_hits: object) -> list[tuple[str, float, tuple[tuple[int, int], ...]]]:
    """Signature over the *stored* hit list, including entries below the report cutoff.

    Deliberately not `.reported`: comparing streaming against prefetched targets must
    prove they agree on every hit pyhmmer retained, not merely on the ones that survive
    filtering.
    """
    return [
        (
            hit.name,
            hit.score,
            tuple((domain.env_from, domain.env_to) for domain in hit.domains),
        )
        for hit in top_hits
    ]


def test_unreported_hit_is_not_recruited_and_prefetch_matches(
    small_fasta: Path, shared_index: Path
) -> None:
    """Sequences that fail --recruit_evalue_cutoff must never reach a family.

    Query 4497037939_1_144 is the pinned example: pyhmmer *stores* 55 hits for it but
    *reports* 54. The 55th, 6320430079, is below the reporting threshold. The legacy
    script iterated the stored list and recruited it; extract_records must not.

    An end-to-end test cannot see this. The stray hit is removed later anyway, by the
    envelope-length filter, so the final families look the same either way. Only a
    direct assertion on extraction can catch a regression here.

    The same TopHits also proves streaming and prefetched targets agree exactly, which
    is what makes --prefetch_targets a pure memory/speed trade.
    """
    logger = logging.getLogger("test-search")
    with (
        pyhmmer.easel.SSIReader(shared_index) as reader,
        pyhmmer.easel.SequenceFile(small_fasta, digital=False, index=reader) as indexed_file,
    ):
        indexed = gf.IndexedSequences(indexed_file)
        seed = gf.run_initial_msa(
            [
                "4497037939_1_144",
                "2536035451_1_130",
                "1500182320",
                "211586000_1_131",
                "1692326086_1_136",
            ],
            indexed,
            4,
        )
        hmm = gf.run_hmmbuild(seed, "raw-query")

    results = []
    for prefetch in (False, True):
        with pyhmmer.easel.SequenceFile(
            small_fasta, digital=True, alphabet=gf.ALPHABET
        ) as targets_file:
            targets = targets_file.read_block() if prefetch else targets_file
            with gf.search(
                [hmm],
                targets,
                cpus=4,
                evalue=0.001,
                logger=logger,
                batch_number=1,
                round_number=1,
            ) as searched:
                results.append(next(searched))
    streaming, prefetched = results

    # The stored list still holds the sub-threshold hit; the reported view does not.
    assert len(streaming) == len(prefetched) == 55
    assert len(streaming.reported) == len(prefetched.reported) == 54
    assert _raw_signature(streaming) == _raw_signature(prefetched)

    stored_names = [hit.name for hit in streaming]
    extracted_names = [record[0] for record in gf.extract_records(streaming)]
    assert "6320430079" in stored_names
    assert "6320430079" not in extracted_names
    assert len(extracted_names) == len(gf.extract_records(prefetched))

    with (
        pyhmmer.easel.SSIReader(shared_index) as reader,
        pyhmmer.easel.SequenceFile(small_fasta, digital=False, index=reader) as indexed_file,
    ):
        # exit_flag=True waives the envelope-length filter, so nothing but the reporting
        # threshold can be keeping 6320430079 out.
        recruited = gf.filter_hits(
            gf.extract_records(streaming),
            streaming.query.M,
            True,
            0.9,
            gf.IndexedSequences(indexed_file),
        )
    assert "6320430079" not in gf.unmask_sequence_names(recruited)


def test_prefetch_end_to_end_equivalence(
    tmp_path: Path,
    baseline_output: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    output = run_pipeline(
        tmp_path / "prefetch",
        cli_args(
            fixture_directory / "clustering.tsv",
            small_fasta,
            fasta_index=shared_index,
            prefetch=True,
        ),
    )
    assert scientific_artifacts(output) == scientific_artifacts(baseline_output)


def test_streaming_does_not_read_block(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original = pyhmmer.easel.SequenceFile.read_block

    def fail_read_block(self: object, *args: object, **kwargs: object) -> object:
        pytest.fail("streaming mode must not call SequenceFile.read_block")

    monkeypatch.setattr(pyhmmer.easel.SequenceFile, "read_block", fail_read_block)
    run_pipeline(
        tmp_path / "streaming",
        cli_args(
            fixture_directory / "clustering.tsv",
            small_fasta,
            fasta_index=shared_index,
        ),
    )
    monkeypatch.setattr(pyhmmer.easel.SequenceFile, "read_block", original)


def test_extracted_records_do_not_retain_pyhmmer_results(
    small_fasta: Path, shared_index: Path
) -> None:
    logger = logging.getLogger("test-records")
    with (
        pyhmmer.easel.SSIReader(shared_index) as reader,
        pyhmmer.easel.SequenceFile(small_fasta, digital=False, index=reader) as indexed_file,
        pyhmmer.easel.SequenceFile(small_fasta, digital=True, alphabet=gf.ALPHABET) as targets,
    ):
        seed = gf.run_initial_msa(
            ["4706047775", "3573189919_79_196", "2632373804_177_299"],
            gf.IndexedSequences(indexed_file),
            1,
        )
        hmm = gf.run_hmmbuild(seed, "records")
        with gf.search(
            [hmm],
            targets,
            cpus=1,
            evalue=0.001,
            logger=logger,
            batch_number=1,
            round_number=1,
        ) as searched:
            records = gf.extract_records(next(searched))
    assert records and all(type(record) is tuple for record in records)
    assert not any(type(item).__name__ == "TopHits" for item in gc.get_referents(records))


@pytest.mark.parametrize(
    ("change", "expected"),
    [
        ({"chunk_num": "../evil"}, "chunk_num"),
        ({"cpus": 0}, "cpus"),
        ({"max_gap_occupancy": 1.1}, "max_gap_occupancy"),
        ({"discard_min_rep_length": 3000}, "minimum"),
    ],
)
def test_validation_precedes_output_creation(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    change: dict[str, object],
    expected: str,
) -> None:
    arguments = gf.parse_args(cli_args(fixture_directory / "clustering.tsv", small_fasta))
    for name, value in change.items():
        setattr(arguments, name, value)
    with contextlib.chdir(tmp_path), pytest.raises(ValueError, match=expected):
        gf.validate_inputs(arguments)
    assert not (tmp_path / "output").exists()


def test_gzip_and_malformed_tsv_rejected_before_outputs(
    tmp_path: Path, fixture_directory: Path, small_fasta: Path
) -> None:
    gzip_options = gf.parse_args(
        cli_args(
            fixture_directory / "clustering.tsv",
            fixture_directory / "mgnifams_input_small.fa.gz",
        )
    )
    with contextlib.chdir(tmp_path), pytest.raises(ValueError, match="uncompressed"):
        gf.validate_inputs(gzip_options)

    malformed = tmp_path / "malformed.tsv"
    malformed.write_text("a\tb\textra\n")
    malformed_options = gf.parse_args(cli_args(malformed, small_fasta))
    with contextlib.chdir(tmp_path), pytest.raises(ValueError, match="exactly two"):
        gf.validate_inputs(malformed_options)
    assert not (tmp_path / "output").exists()


def test_hmmalign_forwards_cpus(monkeypatch: pytest.MonkeyPatch) -> None:
    observed = {}

    def fake_hmmalign(
        _hmm: object, _sequences: object, *, cpus: int, trim: bool
    ) -> pyhmmer.easel.TextMSA:
        observed.update(cpus=cpus, trim=trim)
        return text_msa(["a"], ["AAAA"], "xxxx")

    monkeypatch.setattr(pyhmmer.hmmer, "hmmalign", fake_hmmalign)
    gf.run_hmmalign(object(), [gf.Sequence("a", "AAAA")], 3)
    assert observed == {"cpus": 3, "trim": False}


def test_declared_outputs_parse_and_long_fixture_runs(
    tmp_path: Path,
    baseline_output: Path,
    fixture_directory: Path,
    extra_fasta: Path,
) -> None:
    for directory in gf.FAMILY_DIRECTORIES:
        assert (baseline_output / directory).is_dir()
    for name in (
        "chunk.log",
        "chunk_families.tsv",
        "chunk_discarded.csv",
        "chunk_successful.txt",
        "chunk_converged.txt",
        "chunk_metadata.csv",
        "chunk_reps.fasta.gz",
    ):
        assert (baseline_output / name).is_file()
    for path in (baseline_output / "hmm").glob("*.hmm.gz"):
        with pyhmmer.plan7.HMMFile(path) as hmm_file:
            hmm = hmm_file.read()
        assert hmm is not None
        # The family name lives on the HMM, not in the Stockholm files.
        assert hmm.name == path.name.removesuffix(".hmm.gz")
    for directory in ("seed_msa", "full_msa"):
        for path in (baseline_output / directory).glob("*.sto.gz"):
            family_name = path.name.removesuffix(".sto.gz")
            contents = gzip.decompress(path.read_bytes())
            assert contents.startswith(b"# STOCKHOLM 1.0")
            # The family names itself, and carries no hmmalign posterior annotation.
            assert f"#=GF ID {family_name}".encode() in contents
            assert b"#=GR" not in contents
            assert b"PP_cons" not in contents
            with pyhmmer.easel.MSAFile(path, digital=False) as msa_file:
                msa = msa_file.read()
            assert msa is not None
            # Every row is renumbered onto its parent protein, with no padding left over.
            for name in msa.names:
                assert name == name.strip()
    assert len((baseline_output / "chunk_metadata.csv").read_text().splitlines()) == 3

    long_output = run_pipeline(
        tmp_path / "long",
        cli_args(fixture_directory / "cluster_long.tsv", extra_fasta, chunk="long"),
    )
    assert (long_output / "long_successful.txt").read_text().strip()
    assert (long_output / f"{extra_fasta.name}.ssi").is_file()


def test_v2_full_tsv_fixture_contains_only_clusters_with_four_or_more_members(
    fixture_directory: Path, v2_fasta: Path
) -> None:
    clusters = gf.load_clusters(fixture_directory / "mgnifams_v2.tsv")

    assert len(clusters) == 14
    assert sum(map(len, clusters.values())) == 84
    assert min(map(len, clusters.values())) >= 4
    with v2_fasta.open() as fasta:
        assert sum(line.startswith(">") for line in fasta) == 26_949


def test_v2_full_tsv_end_to_end(
    fixture_directory: Path,
    v2_fasta: Path,
    v2_output: Path,
) -> None:
    representatives = list(gf.load_clusters(fixture_directory / "mgnifams_v2.tsv"))
    discarded = (v2_output / "v2_discarded.csv").read_text().splitlines()
    discarded_representatives = {line.split(",", 1)[0] for line in discarded}
    successful = (v2_output / "v2_successful.txt").read_text().splitlines()

    assert discarded == [
        "3466270235_168_276,family representative length too small,86",
        "825086527_135_242,family representative length too small,99",
    ]
    assert successful == [
        representative
        for representative in representatives
        if representative not in discarded_representatives
    ]
    assert (v2_output / "v2_converged.txt").read_text().splitlines() == [
        str(family_id) for family_id in range(5, 13)
    ]

    metadata = (v2_output / "v2_metadata.csv").read_text().splitlines()
    assert [line.split(",", 1)[0] for line in metadata] == [str(i) for i in range(1, 13)]
    assert len((v2_output / "v2_families.tsv").read_text().splitlines()) == 405
    for directory in gf.FAMILY_DIRECTORIES:
        assert len(list((v2_output / directory).iterdir())) == 12
    assert (v2_output / f"{v2_fasta.name}.ssi").is_file()


def test_output_dir_override(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    run_directory = tmp_path / "run"
    run_directory.mkdir()
    output = tmp_path / "results"
    arguments = cli_args(
        fixture_directory / "clustering.tsv", small_fasta, fasta_index=shared_index
    )
    with contextlib.chdir(run_directory):
        gf.main([*arguments, "--output_dir", str(output)])

    assert all((output / directory).is_dir() for directory in gf.FAMILY_DIRECTORIES)
    assert not (run_directory / "output").exists()


def test_cli_dispatches_to_generate_families(monkeypatch: pytest.MonkeyPatch) -> None:
    """`mgnifam generate_families ...` forwards its arguments untouched."""
    received: list[list[str]] = []
    monkeypatch.setitem(cli.COMMANDS, "generate_families", lambda argv: received.append(list(argv)))
    cli.main(["generate_families", "--cpus", "4", "--prefetch_targets"])
    assert received == [["--cpus", "4", "--prefetch_targets"]]


def test_cli_rejects_unknown_command() -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["merge_families"])
    assert excinfo.value.code == 2


def test_cli_reports_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--version"])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip() == f"mgnifam {__version__}"


def test_console_script_and_module_entry_points_agree() -> None:
    """The `mgnifam` script and `python -m mgnifam` are the same program."""
    script = subprocess.run(["mgnifam", "--version"], capture_output=True, text=True, check=True)
    module = subprocess.run(
        [sys.executable, "-m", "mgnifam", "--version"], capture_output=True, text=True, check=True
    )
    assert script.stdout == module.stdout == f"mgnifam {__version__}\n"


def test_clip_ends_keeps_every_column_above_the_occupancy_threshold() -> None:
    """Both gappy ends are trimmed, and the outermost passing columns survive.

    The legacy script built an end-exclusive `range(start, end)` over inclusive bounds and
    so discarded the last column that passed. Fixed in 1.0.0; this pins the fix.
    """
    # Columns 0-1 and 8-9 are 25% occupied, columns 2-7 are full.
    rows = ["--ABCDEF--", "--ABCDEF--", "--ABCDEF--", "XXABCDEFXX"]
    matrix = np.array([list(row) for row in rows])

    assert gf.calculate_trim_positions(matrix, 0.5) == (2, 7)

    msa = pyhmmer.easel.TextMSA(
        name="t",
        sequences=[
            pyhmmer.easel.TextSequence(name=f"s{index}", sequence=row)
            for index, row in enumerate(rows)
        ],
    )
    clipped = list(gf.clip_ends(msa, 0.5).alignment)

    # Column 7 ("F") passed the threshold and is retained; the gappy flanks are gone.
    assert clipped == ["ABCDEF"] * 4


def test_clip_ends_returns_the_alignment_unchanged_when_no_column_passes() -> None:
    """No qualifying column means nothing to trim, not "trim everything but the last".

    `np.argmax` over an all-False array returns 0, which is how the legacy script silently
    reported the full span here and dropped the final column. `calculate_trim_positions`
    now returns None so the caller can tell the two cases apart.
    """
    rows = ["-A-", "---", "---", "---"]
    matrix = np.array([list(row) for row in rows])
    assert gf.calculate_trim_positions(matrix, 0.5) is None

    msa = pyhmmer.easel.TextMSA(
        name="t",
        sequences=[
            pyhmmer.easel.TextSequence(name=f"s{index}", sequence=row)
            for index, row in enumerate(rows)
        ],
    )
    assert list(gf.clip_ends(msa, 0.5).alignment) == rows


def test_clip_ends_is_identity_when_every_column_passes() -> None:
    rows = ["ABC"] * 4
    matrix = np.array([list(row) for row in rows])
    assert gf.calculate_trim_positions(matrix, 0.5) == (0, 2)

    msa = pyhmmer.easel.TextMSA(
        name="t",
        sequences=[
            pyhmmer.easel.TextSequence(name=f"s{index}", sequence=row)
            for index, row in enumerate(rows)
        ],
    )
    assert list(gf.clip_ends(msa, 0.5).alignment) == rows
