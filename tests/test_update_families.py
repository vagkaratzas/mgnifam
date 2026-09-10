"""Tests for `mgnifam update_families`.

The fixtures do double duty: `small_fasta` is the "old release" a family set is generated
from, and `extra_fasta` is the "new release" that set is then updated against.
"""

import csv
import gzip
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pyhmmer
import pytest

from mgnifam import generate_families, update_families

ARTIFACT_DIRECTORIES = ("seed_msa", "full_msa", "hmm", "rf")


def generate(clusters: Path, fasta: Path, output: Path, chunk: str = "1") -> Path:
    generate_families.main(
        [
            "-c",
            str(clusters),
            "-f",
            str(fasta),
            "-n",
            chunk,
            "-p",
            "2",
            "--output_dir",
            str(output),
        ]
    )
    return output


def update(hmm_input: Path, fasta: Path, output: Path, *extra: str) -> Path:
    update_families.main(
        [
            "-i",
            str(hmm_input),
            "-f",
            str(fasta),
            "-n",
            "9",
            "-p",
            "2",
            "--output_dir",
            str(output),
            *extra,
        ]
    )
    return output


@pytest.fixture(scope="module")
def generated(tmp_path_factory: pytest.TempPathFactory, small_fasta: Path) -> Path:
    """A family set built the ordinary way, to be updated by the tests below."""
    root = tmp_path_factory.mktemp("generated")
    clusters = Path(__file__).parent / "fixtures" / "clustering.tsv"
    return generate(clusters, small_fasta, root)


def family_names(root: Path) -> list[str]:
    return sorted(path.name.removesuffix(".hmm.gz") for path in (root / "hmm").iterdir())


def read_hmm(path: Path) -> pyhmmer.plan7.HMM:
    with pyhmmer.plan7.HMMFile(path) as handle:
        return next(iter(handle))


def write_library(models: list[pyhmmer.plan7.HMM], destination: Path) -> Path:
    with gzip.open(destination, "wb") as handle:
        for hmm in models:
            hmm.write(handle)
    return destination


def test_skip_refine_preserves_the_model_and_writes_no_seed(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    """A recruit-only update changes the membership, never the model.

    The seed alignment and its RF line are unchanged by definition and cannot be rebuilt
    from an HMM, so they are absent rather than fabricated.
    """
    output = update(generated / "hmm", extra_fasta, tmp_path / "out", "--skip_refine")

    for name in family_names(output):
        before = read_hmm(generated / "hmm" / f"{name}.hmm.gz")
        after = read_hmm(output / "hmm" / f"{name}.hmm.gz")
        assert (after.name, after.M, after.consensus) == (before.name, before.M, before.consensus)

    assert list((output / "seed_msa").iterdir()) == []
    assert list((output / "rf").iterdir()) == []
    assert list((output / "full_msa").iterdir()) != []


def test_refine_writes_the_full_artifact_set(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    output = update(generated / "hmm", extra_fasta, tmp_path / "out")
    names = family_names(output)
    assert names
    suffixes = {"rf": ".txt", "hmm": ".hmm.gz", "seed_msa": ".sto.gz", "full_msa": ".sto.gz"}
    for directory, suffix in suffixes.items():
        assert sorted(p.name for p in (output / directory).iterdir()) == [
            f"{name}{suffix}" for name in names
        ]


def test_skip_refine_searches_once_and_refine_at_most_three_times(
    tmp_path: Path, generated: Path, extra_fasta: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One database pass per family recruited, three at most when refining.

    The cost of an update run is `n_families x database x rounds`, so the round count is
    the number worth asserting rather than trusting.
    """
    rounds: list[int] = []
    original = update_families.search

    def counting(*args: object, **kwargs: object) -> object:
        rounds.append(int(kwargs["round_number"]))  # type: ignore[call-overload]
        return original(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(update_families, "search", counting)

    update(generated / "hmm", extra_fasta, tmp_path / "skip", "--skip_refine")
    assert rounds == [1]

    rounds.clear()
    update(generated / "hmm", extra_fasta, tmp_path / "refine")
    assert rounds and max(rounds) <= generate_families.MAX_ROUNDS


def test_identity_is_preserved_in_every_field_not_only_in_filenames(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    """`--chunk_num` is 9 while the names carry a `1_` prefix.

    Any field that rebuilt the name as `f"{chunk}_{id}"` would read `9_1_7` here, which is
    exactly how the representative FASTA annotation was found to be wrong.
    """
    output = update(generated / "hmm", extra_fasta, tmp_path / "out", "--skip_refine")
    names = set(family_names(output))
    assert names and all(name.startswith("1_") for name in names)

    successful = (output / "9_updated_successful.txt").read_text().split()
    assert set(successful) == names

    with (output / "9_updated_metadata.csv").open() as handle:
        metadata = list(csv.DictReader(handle))
    assert {row["family_id"] for row in metadata} == names

    tsv_ids = {
        line.split("\t")[0] for line in (output / "9_updated_families.tsv").read_text().splitlines()
    }
    assert tsv_ids == names

    with gzip.open(output / "9_updated_reps.fasta.gz", "rt") as handle:
        annotations = {line.split("\t")[1].strip() for line in handle if line.startswith(">")}
    assert annotations == names


def test_library_and_directory_inputs_agree_despite_adversarial_ordering(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    """Identity comes from NAME, so neither filename order nor library order may leak in.

    The directory is rebuilt with filenames sorting opposite to their models' names, and
    the library is written in reverse, so any residual dependence on either order shows up
    as a byte difference.
    """
    names = family_names(generated)
    models = [read_hmm(generated / "hmm" / f"{name}.hmm.gz") for name in names]

    shuffled = tmp_path / "shuffled"
    shuffled.mkdir()
    for position, (name, hmm) in enumerate(zip(names, models, strict=True)):
        # Filename order is the reverse of NAME order.
        with gzip.open(shuffled / f"{len(names) - position:04d}.hmm.gz", "wb") as handle:
            hmm.write(handle)
        assert hmm.name == name

    library = write_library(list(reversed(models)), tmp_path / "library.hmm.gz")

    from_directory = update(shuffled, extra_fasta, tmp_path / "dir", "--skip_refine")
    from_library = update(library, extra_fasta, tmp_path / "lib", "--skip_refine")

    for path in sorted(from_directory.rglob("*")):
        if not path.is_file() or path.suffix == ".log" or path.suffix == ".ssi":
            continue
        other = from_library / path.relative_to(from_directory)
        assert other.read_bytes() == path.read_bytes(), path.name


def test_outputs_are_invariant_to_cpus_batch_size_prefetch_and_hash_seed(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    baseline = update(generated / "hmm", extra_fasta, tmp_path / "a", "--skip_refine")
    variants = [
        update(generated / "hmm", extra_fasta, tmp_path / "b", "--skip_refine", "-p", "1"),
        update(
            generated / "hmm", extra_fasta, tmp_path / "c", "--skip_refine", "--batch_size", "16"
        ),
        update(
            generated / "hmm", extra_fasta, tmp_path / "d", "--skip_refine", "--prefetch_targets"
        ),
    ]
    for variant in variants:
        for path in sorted(baseline.rglob("*")):
            if not path.is_file() or path.suffix in (".log", ".ssi"):
                continue
            assert (variant / path.relative_to(baseline)).read_bytes() == path.read_bytes()

    # A separate process, because PYTHONHASHSEED is read at interpreter start.
    seeded = tmp_path / "e"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "mgnifam",
            "update_families",
            "-i",
            str(generated / "hmm"),
            "-f",
            str(extra_fasta),
            "-n",
            "9",
            "-p",
            "2",
            "--skip_refine",
            "--output_dir",
            str(seeded),
        ],
        check=True,
        env={**os.environ, "PYTHONHASHSEED": "12345"},
    )
    for path in sorted(baseline.rglob("*")):
        if not path.is_file() or path.suffix in (".log", ".ssi"):
            continue
        assert (seeded / path.relative_to(baseline)).read_bytes() == path.read_bytes()


# No whitespace case: HMMER truncates NAME at the first space, so `"with space"` comes back
# from a round trip as `"with"` and can never reach the loader in that form.
@pytest.mark.parametrize(
    "name",
    ["../escape", "with/slash", "with,comma", "..", "."],
    ids=["traversal", "slash", "comma", "dotdot", "dot"],
)
def test_unsafe_family_names_are_rejected_before_any_output_exists(
    tmp_path: Path, generated: Path, small_fasta: Path, name: str
) -> None:
    """A NAME reaches artifact paths and unquoted CSV fields, and comes from a foreign file.

    Traversal would place -- and then roll back, i.e. unlink -- files outside the output
    tree; a comma would silently shift every column of two CSVs.
    """
    hmm = read_hmm(generated / "hmm" / f"{family_names(generated)[0]}.hmm.gz")
    hmm.name = name
    library = write_library([hmm], tmp_path / "bad.hmm.gz")
    output = tmp_path / "out"

    with pytest.raises(ValueError, match="not a usable family name"):
        update(library, small_fasta, output)
    assert not output.exists()


def test_rejections_happen_before_the_output_tree_is_created(
    tmp_path: Path, generated: Path, small_fasta: Path, fixture_directory: Path
) -> None:
    names = family_names(generated)
    models = [read_hmm(generated / "hmm" / f"{name}.hmm.gz") for name in names]

    duplicated = write_library([models[0], models[0]], tmp_path / "dup.hmm.gz")
    empty = tmp_path / "empty"
    empty.mkdir()

    cases: list[tuple[Path, Path, str]] = [
        (duplicated, small_fasta, "duplicate HMM NAME"),
        (empty, small_fasta, "contains no HMMs"),
        (
            generated / "hmm",
            fixture_directory / "mgnifams_input_small.fa.gz",
            "must be uncompressed",
        ),
    ]
    for index, (hmm_input, fasta, message) in enumerate(cases):
        output = tmp_path / f"out{index}"
        with pytest.raises(ValueError, match=message):
            update(hmm_input, fasta, output)
        assert not output.exists()

    # An output root whose `hmm/` is the input directory would destroy the models the user
    # needs in order to re-run.
    with pytest.raises(ValueError, match="must not be the output_dir"):
        update(generated / "hmm", small_fasta, generated)


def test_nucleotide_models_are_rejected_rather_than_dying_inside_hmmsearch(
    tmp_path: Path, small_fasta: Path
) -> None:
    alphabet = pyhmmer.easel.Alphabet.dna()
    msa = pyhmmer.easel.TextMSA(
        name=b"dna_family",
        sequences=[
            pyhmmer.easel.TextSequence(name="a", sequence="ACGTACGTACGT"),
            pyhmmer.easel.TextSequence(name="b", sequence="ACGTACGTACGA"),
        ],
    ).digitize(alphabet)
    builder = pyhmmer.plan7.Builder(alphabet, seed=42)
    hmm, _, _ = builder.build_msa(msa, pyhmmer.plan7.Background(alphabet))

    library = tmp_path / "dna.hmm.gz"
    with gzip.open(library, "wb") as handle:
        hmm.write(handle)

    output = tmp_path / "out"
    with pytest.raises(ValueError, match="not an amino-acid model"):
        update(library, small_fasta, output)
    assert not output.exists()


def test_delta_reports_every_family_with_a_constant_column_count(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    output = update(generated / "hmm", extra_fasta, tmp_path / "out", "--skip_refine")
    with (output / "9_updated_delta.csv").open() as handle:
        rows = list(csv.DictReader(handle))

    assert {row["family_id"] for row in rows} == set(family_names(generated))
    assert all(len(row) == len(update_families.DELTA_HEADER.split(",")) for row in rows)
    assert all(None not in row.values() for row in rows), "a row had fewer fields than the header"

    successful = set((output / "9_updated_successful.txt").read_text().split())
    for row in rows:
        if row["family_id"] in successful:
            assert row["outcome"] == "successful"
            # Retention is 1.0 by construction: the yardstick is round 1's own recruits,
            # and the exit branch re-filters the same records with the length requirement
            # waived, so it can only be a superset.
            assert float(row["retention"]) == 1.0
            assert row["rounds_run"] == "1"
            assert row["model_length_before"] == row["model_length_after"]


def test_a_family_with_no_hits_is_reported_distinctly_from_a_low_complexity_one(
    tmp_path: Path, generated: Path, small_fasta: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Two different answers that the shared code path would have collapsed into one.

    "no hits" means the model found nothing in the new release. "low complexity" means it
    found hits and none were long enough. For an update run that difference is the report.
    """
    monkeypatch.setattr(update_families, "extract_records", lambda _hits: [])
    output = update(generated / "hmm", small_fasta, tmp_path / "out", "--skip_refine")

    with (output / "9_updated_delta.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows and all(row["outcome"] == update_families.NO_HITS for row in rows)
    # No yardstick ever existed, so retention is empty rather than zero -- scoring it would
    # divide by an empty set.
    assert all(row["retention"] == "" for row in rows)
    assert all(row["round1_recruits"] == "" for row in rows)
    assert all(row["full_msa_size"] == "" for row in rows)

    discarded = (output / "9_updated_discarded.csv").read_text().splitlines()[1:]
    assert len(discarded) == len(rows)
    assert (output / "9_updated_successful.txt").read_text() == ""


def test_every_discard_produces_exactly_one_delta_row(
    tmp_path: Path, generated: Path, extra_fasta: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A contained per-family crash is re-emitted; the delta must not be written twice."""
    calls = {"n": 0}
    original = generate_families.Family.finish

    def exploding(self: generate_families.Family, *args: object, **kwargs: object) -> None:
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("boom")
        original(self, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(generate_families.Family, "finish", exploding)

    with pytest.raises(SystemExit) as exit_info:
        update(generated / "hmm", extra_fasta, tmp_path / "out", "--skip_refine")
    assert exit_info.value.code == generate_families.EXIT_CRASHED_FAMILIES

    output = tmp_path / "out"
    with (output / "9_updated_delta.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["family_id"] for row in rows]
    assert len(ids) == len(set(ids)) == len(family_names(generated))
    assert sum(1 for row in rows if row["outcome"].startswith("internal error")) == 1


def test_a_delta_failure_after_a_discard_append_is_fatal_and_not_re_emitted(
    tmp_path: Path, generated: Path, small_fasta: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The discard branch's shared appends are one commit, like the success branch's.

    Without a boundary here the failure is contained by `family_guard`, `main` re-emits,
    and the discard row is appended a second time.
    """
    monkeypatch.setattr(update_families, "extract_records", lambda _hits: [])

    def refuse(_writers: object, _row: str | None) -> None:
        raise OSError("no space left on device")

    monkeypatch.setattr(generate_families, "write_delta", refuse)

    with pytest.raises(SystemExit) as exit_info:
        update(generated / "hmm", small_fasta, tmp_path / "out", "--skip_refine")
    assert exit_info.value.code == 1

    discarded = (tmp_path / "out" / "9_updated_discarded.csv").read_text().splitlines()[1:]
    assert len(discarded) == 1, "the discard row was written twice by the re-emit"


def test_rerunning_the_same_models_into_the_same_directory_is_allowed(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    """The Nextflow retry path: a chunk that failed is re-run unchanged."""
    output = tmp_path / "out"
    update(generated / "hmm", extra_fasta, output, "--skip_refine")
    first = sorted(p.name for p in (output / "hmm").iterdir())

    update(generated / "hmm", extra_fasta, output, "--skip_refine")
    assert sorted(p.name for p in (output / "hmm").iterdir()) == first


@pytest.mark.parametrize("outcome", ["discard", "crash", "skip_refine"])
def test_retry_artifacts_match_the_new_outcome(
    tmp_path: Path,
    generated: Path,
    extra_fasta: Path,
    monkeypatch: pytest.MonkeyPatch,
    outcome: str,
) -> None:
    output = tmp_path / "out"
    update(generated / "hmm", extra_fasta, output)
    assert list((output / "seed_msa").iterdir())
    args = [] if outcome == "crash" else ["--skip_refine"]
    if outcome == "discard":
        args += ["--discard_min_rep_length", "2000"]
    if outcome == "crash":
        original = generate_families.renumber_msa
        failed = False

        def fail_once(*args: object, **kwargs: object) -> object:
            nonlocal failed
            if not failed:
                failed = True
                raise OSError("transient read failure")
            return original(*args, **kwargs)  # type: ignore[arg-type]

        monkeypatch.setattr(generate_families, "renumber_msa", fail_once)
        with pytest.raises(SystemExit) as error:
            update(generated / "hmm", extra_fasta, output, *args)
        assert error.value.code == 3
    else:
        update(generated / "hmm", extra_fasta, output, *args)
    successful = set((output / "9_updated_successful.txt").read_text().splitlines())
    assert bool(successful) == (outcome != "discard")
    for directory, suffix in update_families.ARTIFACT_SUFFIXES.items():
        expected = successful if directory in ("hmm", "full_msa") or outcome == "crash" else set()
        assert {p.name.removesuffix(suffix) for p in (output / directory).iterdir()} == expected
    with (output / "9_updated_discarded.csv").open() as handle:
        discarded = {row["representative"] for row in csv.DictReader(handle)}
    assert successful.isdisjoint(discarded)
    assert successful | discarded == set(family_names(generated))


def test_failed_retry_cleanup_aborts_before_aggregate_truncation(
    tmp_path: Path, generated: Path, extra_fasta: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "out"
    update(generated / "hmm", extra_fasta, output, "--skip_refine")
    before = {p: p.read_bytes() for p in output.iterdir() if p.is_file()}

    def refuse(*args: object, **kwargs: object) -> None:
        raise OSError("cleanup failed")

    monkeypatch.setattr(Path, "unlink", refuse)
    with pytest.raises(OSError, match="cleanup failed"):
        update(generated / "hmm", extra_fasta, output, "--skip_refine")
    assert {p: p.read_bytes() for p in output.iterdir() if p.is_file()} == before


def test_a_directory_holding_another_runs_families_is_refused(
    tmp_path: Path, generated: Path, extra_fasta: Path
) -> None:
    """Refused rather than cleared, because an updated family's name cannot be derived.

    `generate_families` clears its own past output from a `<chunk>_<rank>` pattern. An
    updated family keeps its model's name and `--chunk_num` never enters a per-family
    filename, so nothing on disk says which run wrote `hmm/1_7.hmm.gz`. Writing a smaller
    set over a larger one would leave the dropped families beside aggregates that no longer
    list them.
    """
    output = tmp_path / "out"
    names = family_names(generated)
    assert len(names) >= 2
    update(generated / "hmm", extra_fasta, output, "--skip_refine")

    subset = tmp_path / "subset"
    subset.mkdir()
    kept = names[0]
    with gzip.open(subset / f"{kept}.hmm.gz", "wb") as handle:
        read_hmm(generated / "hmm" / f"{kept}.hmm.gz").write(handle)

    with pytest.raises(ValueError, match="does not update"):
        update(subset, extra_fasta, output, "--skip_refine")
    # Refusing must not have half-cleared the directory it refused to write into.
    assert family_names(output) == names


def test_a_partial_run_can_be_rerun_in_place(
    tmp_path: Path, generated: Path, extra_fasta: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A run that exits 1 leaves artifacts behind; re-running it unchanged must work.

    Artifacts are written before the first shared append, so a commit failure strands them
    with no name in `successful.txt`. They all belong to this chunk's own model set, so the
    rerun overwrites them rather than being refused.
    """
    output = tmp_path / "out"

    class DeadSink:
        def write(self, _text: str) -> int:
            raise OSError("sink is dead")

    original_open_writers = update_families.open_writers

    def with_a_dead_sink(*args: object, **kwargs: object) -> object:
        writers = original_open_writers(*args, **kwargs)  # type: ignore[arg-type]
        writers.successful_clusters = DeadSink()  # type: ignore[assignment]
        return writers

    monkeypatch.setattr(update_families, "open_writers", with_a_dead_sink)
    with pytest.raises(SystemExit) as exit_info:
        update(generated / "hmm", extra_fasta, output, "--skip_refine")
    assert exit_info.value.code == 1
    assert family_names(output), "expected artifacts stranded by the failed run"
    assert (output / "9_updated_successful.txt").read_text() == ""

    monkeypatch.undo()
    update(generated / "hmm", extra_fasta, output, "--skip_refine")
    assert set((output / "9_updated_successful.txt").read_text().split()) == set(
        family_names(generated)
    )


@pytest.mark.parametrize("alias", ["direct", "symlink", "hardlink", "directory_symlink"])
@pytest.mark.parametrize("model_count", [1, 2])
def test_input_models_cannot_be_overwritten(
    tmp_path: Path, generated: Path, small_fasta: Path, alias: str, model_count: int
) -> None:
    models = [read_hmm(p) for p in sorted((generated / "hmm").iterdir())][:model_count]
    output = tmp_path / "out"
    (output / "hmm").mkdir(parents=True)
    destination = output / "hmm" / f"{models[0].name}.hmm.gz"
    source = write_library(models, destination if alias == "direct" else tmp_path / "input.hmm.gz")
    if alias == "symlink":
        destination.symlink_to(source)
    elif alias == "hardlink":
        destination.hardlink_to(source)
    elif alias == "directory_symlink":
        source.unlink()
        (output / "hmm").rmdir()
        source = tmp_path / "models"
        source.mkdir()
        write_library(models, source / destination.name)
        (output / "hmm").symlink_to(source, target_is_directory=True)
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    with pytest.raises(ValueError, match=r"overlaps an input|must not be the output_dir"):
        update(source, small_fasta, output, "--skip_refine")
    assert {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()} == before
    assert not (output / "9_updated.log").exists()


def test_input_library_cannot_alias_an_aggregate(
    tmp_path: Path, generated: Path, small_fasta: Path
) -> None:
    output = tmp_path / "out"
    output.mkdir()
    source = write_library(
        [read_hmm(next((generated / "hmm").iterdir()))], output / "9_updated_reps.fasta.gz"
    )
    before = source.read_bytes()
    with pytest.raises(ValueError, match="overlaps an input"):
        update(source, small_fasta, output, "--skip_refine")
    assert source.read_bytes() == before
    assert list(output.iterdir()) == [source]


def test_round_one_recruits_are_adopted_as_members_exactly_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Only the first round sets the yardstick; later rounds are scored against it.

    If a later round re-adopted its own recruits, `finish` would compare the final model to
    itself and the membership check could never fail -- which is the whole point of it.
    `advance` clears the flag in the same block that reads it.
    """
    rounds = iter(
        [
            [generate_families.Sequence("A_1_9/1_9", "AAAA")],
            [generate_families.Sequence("B_1_9/1_9", "BBBB")],
        ]
    )
    monkeypatch.setattr(generate_families, "filter_hits", lambda *_a, **_k: next(rounds))
    options = SimpleNamespace(recruit_hit_length_percentage=0.9)
    family = generate_families.Family(
        representative="1_7", members=[], adopt_recruits_as_members=True
    )

    family.advance(options, None, generate_families.MAX_ROUNDS)  # type: ignore[arg-type]
    assert family.members == ["A_1_9"]
    assert family.adopt_recruits_as_members is False

    family.advance(options, None, generate_families.MAX_ROUNDS)  # type: ignore[arg-type]
    assert family.members == ["A_1_9"], "a later round overwrote round 1's yardstick"


def test_cli_dispatches_to_update_families(monkeypatch: pytest.MonkeyPatch) -> None:
    from mgnifam import cli

    seen: list[list[str]] = []
    monkeypatch.setitem(cli.COMMANDS, "update_families", lambda argv: seen.append(list(argv)))
    cli.main(["update_families", "-i", "models", "-f", "db.fa"])
    assert seen == [["-i", "models", "-f", "db.fa"]]
