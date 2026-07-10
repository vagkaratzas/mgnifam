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

import pyhmmer
import pytest

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


def run_pipeline(directory: Path, arguments: list[str]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    with contextlib.chdir(directory):
        gf.main(arguments)


def scientific_artifacts(directory: Path) -> dict[str, bytes]:
    artifacts = {}
    for output_directory in gf.OUTPUT_DIRECTORIES:
        if output_directory == "logs":
            continue
        for path in sorted((directory / output_directory).glob("**/*")):
            if path.is_file():
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
    run_pipeline(
        output,
        cli_args(
            fixture_directory / "clustering.tsv",
            small_fasta,
            fasta_index=shared_index,
        ),
    )
    return output


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
    family.advance(family_options(), FakeSequences({"a": "AAAA"}))
    assert family.state is gf.FamilyState.DISCARDED
    assert family.discard_reason == "too few sequences after redundancy filtering"
    assert family.discard_value == 2


def test_converged_discard_keeps_provisional_id_for_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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
    family.advance(options, store)
    assert family.state is gf.FamilyState.CONVERGED
    family.finish(options, store)
    assert family.state is gf.FamilyState.DISCARDED
    assert family.ever_converged

    for directory in ("rf", "hmm", "seed_msa_sto", "full_msa_sto"):
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
    success_count = gf.emit_family(family, 0, "chunk", writers, tmp_path)
    assert success_count == 0
    assert writers.converged_families.getvalue() == "1\n"

    successful = gf.Family(
        "b",
        ["a", "b"],
        state=gf.FamilyState.SUCCESSFUL,
        seed_msa=seed,
        full_msa=text_msa(["a", "b"], ["AAAA", "AAAT"], "xxxx"),
        full_msa_num_seqs=2,
    )
    success_count = gf.emit_family(successful, success_count, "chunk", writers, tmp_path)
    assert success_count == 1
    assert successful.family_id == 1


def test_cpus_and_sanity_anchors(
    tmp_path: Path,
    baseline_output: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    cpus_four = tmp_path / "cpus-four"
    run_pipeline(
        cpus_four,
        cli_args(
            fixture_directory / "clustering.tsv",
            small_fasta,
            cpus=4,
            fasta_index=shared_index,
        ),
    )
    assert scientific_artifacts(cpus_four) == scientific_artifacts(baseline_output)
    assert (baseline_output / "successful_clusters" / "chunk.txt").read_text().splitlines() == [
        "4706047775",
        "1622851798_832_939",
        "4497037939_1_144",
    ]
    metadata = (baseline_output / "family_metadata" / "chunk.csv").read_text().splitlines()
    assert [line.split(",")[2].strip('"') for line in metadata] == [
        "782510898",
        "5761513631",
        "1446399400",
    ]
    assert (baseline_output / "converged_families" / "chunk.txt").read_text() == "2\n"


def test_batch_size_invariance(
    tmp_path: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    outputs = []
    for batch_size in (1, 64):
        output = tmp_path / f"batch-{batch_size}"
        run_pipeline(
            output,
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
            for line in (output / "family_metadata" / "chunk.csv").read_text().splitlines()
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
        output = tmp_path / f"hash-{hash_seed}"
        output.mkdir()
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = hash_seed
        subprocess.run(
            [sys.executable, "-m", "mgnifam.generate_families", *arguments],
            cwd=output,
            env=environment,
            check=True,
        )
        outputs.append(output)
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
    output = tmp_path / "rerun"
    arguments = cli_args(
        fixture_directory / "clustering.tsv", small_fasta, fasta_index=shared_index
    )
    run_pipeline(output, arguments)
    reduced_clusters = tmp_path / "reduced.tsv"
    reduced_clusters.write_text(
        "\n".join((fixture_directory / "clustering.tsv").read_text().splitlines()[:9]) + "\n"
    )
    run_pipeline(output, cli_args(reduced_clusters, small_fasta, fasta_index=shared_index))
    for directory in ("seed_msa_sto", "full_msa_sto", "hmm", "rf"):
        assert all(
            "chunk_2." not in path.name and "chunk_3." not in path.name
            for path in (output / directory).iterdir()
        )


def _raw_signature(top_hits: object) -> list[tuple[str, float, tuple[tuple[int, int], ...]]]:
    return [
        (
            hit.name,
            hit.score,
            tuple((domain.env_from, domain.env_to) for domain in hit.domains),
        )
        for hit in top_hits
    ]


def test_raw_unreported_hit_and_prefetch_equivalence(small_fasta: Path, shared_index: Path) -> None:
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
    assert len(streaming) == len(prefetched) == 55
    assert len(streaming.reported) == len(prefetched.reported) == 54
    assert _raw_signature(streaming) == _raw_signature(prefetched)
    assert "6320430079" in [record[0] for record in gf.extract_records(streaming)]
    with (
        pyhmmer.easel.SSIReader(shared_index) as reader,
        pyhmmer.easel.SequenceFile(small_fasta, digital=False, index=reader) as indexed_file,
    ):
        recruited = gf.filter_hits(
            gf.extract_records(streaming),
            streaming.query.M,
            True,
            0.9,
            gf.IndexedSequences(indexed_file),
        )
    assert "6320430079" in gf.unmask_sequence_names(recruited)


def test_prefetch_end_to_end_equivalence(
    tmp_path: Path,
    baseline_output: Path,
    fixture_directory: Path,
    small_fasta: Path,
    shared_index: Path,
) -> None:
    output = tmp_path / "prefetch"
    run_pipeline(
        output,
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
    assert not any((tmp_path / directory).exists() for directory in gf.OUTPUT_DIRECTORIES)


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
    assert not any((tmp_path / directory).exists() for directory in gf.OUTPUT_DIRECTORIES)


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
    for directory in gf.OUTPUT_DIRECTORIES:
        assert (baseline_output / directory).is_dir()
    for path in (baseline_output / "hmm").glob("*.hmm.gz"):
        with pyhmmer.plan7.HMMFile(path) as hmm_file:
            hmm = hmm_file.read()
        assert hmm is not None
        # The family name lives on the HMM, not in the Stockholm files.
        assert hmm.name == path.name.removesuffix(".hmm.gz")
    for path in (baseline_output / "seed_msa_sto").glob("*.sto.gz"):
        contents = gzip.decompress(path.read_bytes())
        # Legacy's renumber_sto_msa strips every #=GF line; so must ours.
        assert b"#=GF" not in contents
        assert contents.startswith(b"# STOCKHOLM 1.0")
        with pyhmmer.easel.MSAFile(path, digital=False) as msa_file:
            assert msa_file.read() is not None
    for path in (baseline_output / "full_msa_sto").glob("*.sto.gz"):
        with pyhmmer.easel.MSAFile(path, digital=False) as msa_file:
            assert msa_file.read() is not None
    assert len((baseline_output / "family_metadata" / "chunk.csv").read_text().splitlines()) == 3

    long_output = tmp_path / "long"
    run_pipeline(
        long_output,
        cli_args(fixture_directory / "cluster_long.tsv", extra_fasta, chunk="long"),
    )
    assert (long_output / "successful_clusters" / "long.txt").read_text().strip()
