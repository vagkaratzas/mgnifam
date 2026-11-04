import gzip
import os
import random
import subprocess
import sys
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from mgnifam import cli
from mgnifam import compute_clans as cc

FAMILY_ROWS = list[tuple[str, str]]


def write_members(path: Path, rows: FAMILY_ROWS) -> Path:
    path.write_text("".join(f"{family}\t{member}\n" for family, member in rows))
    return path


def write_fasta(path: Path, lengths: dict[str, int]) -> Path:
    # Wrapped at 7 residues so lengths span several lines.
    with path.open("w") as handle:
        for name, length in lengths.items():
            sequence = "A" * length
            lines = [sequence[i : i + 7] for i in range(0, length, 7)]
            handle.write(f">{name} description\n" + "".join(f"{line}\n" for line in lines))
    return path


def run(tmp_path: Path, members: list[Path], *extra: str, out: str = "out", chunk: str = "t"):
    argv = ["-i", *map(str, members), "-n", chunk, "--output_dir", str(tmp_path / out), "-p", "1"]
    try:
        cc.main([*argv, *extra])
    except SystemExit as stop:
        return stop.code, tmp_path / out
    return 0, tmp_path / out


def edges(root: Path, chunk: str = "t") -> dict[tuple[str, str], list[str]]:
    lines = (root / f"{chunk}_clan_edges.csv").read_text().splitlines()
    assert lines[0] == cc.EDGES_HEADER.strip()
    return {tuple(line.split(",")[:2]): line.split(",")[2:] for line in lines[1:]}


def clans(root: Path, chunk: str = "t") -> list[tuple[str, str, str]]:
    lines = (root / f"{chunk}_clans.csv").read_text().splitlines()
    assert lines[0] == cc.CLANS_HEADER.strip()
    return [tuple(line.split(",")) for line in lines[1:]]


def outputs(root: Path, chunk: str = "t") -> list[bytes]:
    names = ("clan_edges.csv", "clans.csv", "clans_summary.json")
    return [(root / f"{chunk}_{name}").read_bytes() for name in names]


def tiny_sizes(monkeypatch: pytest.MonkeyPatch, **overrides: int) -> None:
    sizes = cc.Sizes(
        range_bytes=64, run_rows=5, unit_cost=cc.MEMBER_ROW_COST * 80, emission_batch=3
    )
    monkeypatch.setattr(cc, "sizes_for", lambda _: sizes._replace(**overrides))


# --- oracle ------------------------------------------------------------------------------


def oracle(rows: FAMILY_ROWS, lengths: dict[str, int], chunk: str = "t", **flags: Fraction):
    """Edges and clans from explicit position sets, with the command's selection rules."""
    jaccard = flags.get("min_edge_jaccard", Fraction(1, 10))
    containment = flags.get("min_edge_containment", Fraction(1, 2))
    clan_jaccard = flags.get("clan_jaccard", Fraction(1, 2))
    positions: dict[str, set[tuple[str, int]]] = defaultdict(set)
    proteins: dict[str, set[str]] = defaultdict(set)
    for family, member in rows:
        sliced = cc.split_member(member.encode())
        if sliced is None:
            protein, start, end = member, 1, lengths[member]
        else:
            protein, start, end = sliced[0].decode(), sliced[1], sliced[2]
        positions[family].update((protein, p) for p in range(start, end + 1))
        proteins[family].add(protein)
    names = sorted(positions)
    edge_lines, clan_edges = [], []
    for x, a in enumerate(names):
        for b in names[x + 1 :]:
            shared = len(proteins[a] & proteins[b])
            if not shared:
                continue
            i = len(positions[a] & positions[b])
            union = len(positions[a]) + len(positions[b]) - i
            smaller = min(len(positions[a]), len(positions[b]))
            n_union = len(proteins[a]) + len(proteins[b]) - shared
            if Fraction(i, union) >= jaccard or Fraction(i, smaller) >= containment:
                edge_lines.append(
                    f"{a},{b},{shared},{shared / n_union:.6f},{i / union:.6f},{i / smaller:.6f}"
                )
            if Fraction(i, union) >= clan_jaccard:
                clan_edges.append((a, b, i / union))
    degree: dict[str, float] = defaultdict(float)
    for a, _, weight in clan_edges:
        degree[a] += weight
    for _, b, weight in clan_edges:
        degree[b] += weight
    graph: dict[str, set[str]] = defaultdict(set)
    for a, b, _ in clan_edges:
        graph[a].add(b)
        graph[b].add(a)
    seen: set[str] = set()
    components = []
    for start in sorted(graph):
        if start in seen:
            continue
        queue, component = deque([start]), []
        seen.add(start)
        while queue:
            node = queue.popleft()
            component.append(node)
            for other in graph[node] - seen:
                seen.add(other)
                queue.append(other)
        component.sort()
        best = max(degree[f] for f in component)
        representative = min(f for f in component if degree[f] == best)
        components.append((representative, component))
    components.sort(key=lambda c: (-len(c[1]), c[0]))
    clan_rows = [
        (f"{chunk}_clan_{number}", family, str(family == representative))
        for number, (representative, component) in enumerate(components, start=1)
        for family in component
    ]
    return edge_lines, clan_rows


def random_case(seed: int) -> tuple[FAMILY_ROWS, dict[str, int]]:
    generator = random.Random(seed)
    lengths = {
        f"p{n}" + ("/v1" if n % 7 == 0 else ""): generator.randint(20, 120) for n in range(300)
    }
    proteins = list(lengths)
    rows: FAMILY_ROWS = []
    for family in range(30):
        # Families cluster on a few domains so that many pairs overlap.
        domain = generator.randint(1, 60)
        for _ in range(generator.randint(5, 60)):
            protein = generator.choice(proteins[family * 8 : family * 8 + 80])
            if generator.random() < 0.05:
                rows.append((f"f{family}", protein))
                continue
            start = max(1, min(lengths[protein], domain + generator.randint(-10, 10)))
            end = min(lengths[protein], start + generator.randint(0, 40))
            rows.append((f"f{family}", f"{protein}/{start}-{end}"))
            if generator.random() < 0.2:  # a repeat or a duplicate row
                rows.append((f"f{family}", f"{protein}/{max(1, start - 3)}-{end}"))
    # Near-copies of five families, so that clans form.
    rows += [
        (f"g{family[1:]}", member)
        for family, member in rows
        if family in {"f0", "f1", "f2", "f3", "f4"} and generator.random() < 0.8
    ]
    generator.shuffle(rows)
    return rows, lengths


# --- the issue's cases -------------------------------------------------------------------


def test_whole_protein_member_takes_its_length_from_the_fasta(tmp_path: Path, capsys) -> None:
    members = write_members(tmp_path / "m.tsv", [("A", "P"), ("B", "P/1-10")])
    code, _ = run(tmp_path, [members])
    assert code == 2 and "--fasta_file" in capsys.readouterr().err
    fasta = write_fasta(tmp_path / "db.fa", {"P": 40, "Q": 5})
    code, root = run(tmp_path, [members], "-f", str(fasta))
    assert code == 0
    # A covers 1-40, B 1-10: Jaccard 10/40, containment 1.
    assert edges(root) == {("A", "B"): ["1", "1.000000", "0.250000", "1.000000"]}

    missing = write_fasta(tmp_path / "other.fa", {"Q": 5})
    code, _ = run(tmp_path, [members], "-f", str(missing), out="missing")
    error = capsys.readouterr().err
    assert code == 2 and "'P'" in error and not (tmp_path / "missing").exists()


def test_repeat_domains_and_duplicate_rows_count_once(tmp_path: Path) -> None:
    rows = [("A", "P/1-10"), ("A", "P/5-20"), ("A", "P/5-20"), ("A", "P/21-30"), ("B", "P/1-30")]
    code, root = run(tmp_path, [write_members(tmp_path / "m.tsv", rows)])
    assert code == 0
    # A merges to 1-30 (overlapping and adjacent), identical to B.
    assert edges(root) == {("A", "B"): ["1", "1.000000", "1.000000", "1.000000"]}
    assert clans(root) == [("t_clan_1", "A", "True"), ("t_clan_1", "B", "False")]


def test_different_domains_of_the_same_proteins_are_not_similar(tmp_path: Path) -> None:
    rows = [
        (f, f"P{n}/{start}-{start + 49}") for n in range(5) for f, start in (("A", 1), ("B", 60))
    ]
    code, root = run(tmp_path, [write_members(tmp_path / "m.tsv", rows)])
    assert code == 0
    assert edges(root) == {} and clans(root) == []
    summary = (root / "t_clans_summary.json").read_text()
    assert '"pairs_sharing_proteins": 1' in summary
    # The old protein-level measure would have called them identical.
    _, protein_level = run(tmp_path, [tmp_path / "m.tsv"], "--min_edge_jaccard", "0", out="all")
    assert edges(protein_level)[("A", "B")][:3] == ["5", "1.000000", "0.000000"]


def test_nested_family_is_reported_by_containment_not_clustered(tmp_path: Path) -> None:
    rows = [("big", f"P{n}/1-100") for n in range(20)] + [("small", "P0/10-30")]
    code, root = run(tmp_path, [write_members(tmp_path / "m.tsv", rows)])
    assert code == 0
    assert edges(root) == {("big", "small"): ["1", "0.050000", "0.010500", "1.000000"]}
    assert clans(root) == []


def test_clan_representative_ordering_and_exact_thresholds(tmp_path: Path) -> None:
    rows = [
        ("A", "P/1-60"),  # A-B 50/70, B-C 45/75, A-C 35/85: a chain through B
        ("B", "P/11-70"),
        ("C", "P/26-85"),
        ("D", "Q/1-10"),  # D and E identical: equal degree, earlier name wins
        ("E", "Q/1-10"),
        ("F", "R/1-20"),  # F-G exactly 1/2
        ("G", "R/1-10"),
    ]
    code, root = run(tmp_path, [write_members(tmp_path / "m.tsv", rows)])
    assert code == 0
    assert clans(root) == [
        ("t_clan_1", "A", "False"),
        ("t_clan_1", "B", "True"),
        ("t_clan_1", "C", "False"),
        ("t_clan_2", "D", "True"),
        ("t_clan_2", "E", "False"),
        ("t_clan_3", "F", "True"),
        ("t_clan_3", "G", "False"),
    ]
    assert edges(root)[("F", "G")][2] == "0.500000"
    _, strict = run(tmp_path, [tmp_path / "m.tsv"], "--clan_jaccard", "0.500001", out="strict")
    assert ("t_clan_3", "F", "True") not in clans(strict)


def test_pair_orientation_does_not_follow_sweep_order(tmp_path: Path) -> None:
    # B starts first but A sorts first: I(A, B) must land in the (A, B) cell.
    rows = [("B", "P/1-10"), ("A", "P/5-10"), ("C", "Q/1-5")]
    code, root = run(tmp_path, [write_members(tmp_path / "m.tsv", rows)])
    assert code == 0
    assert edges(root) == {("A", "B"): ["1", "1.000000", "0.600000", "1.000000"]}


@pytest.mark.parametrize("rows", [[], [("A", "P/1-10"), ("A", "Q/1-10")]])
def test_fewer_than_two_families(tmp_path: Path, rows: FAMILY_ROWS) -> None:
    code, root = run(tmp_path, [write_members(tmp_path / "m.tsv", rows)])
    assert code == 0
    assert edges(root) == {} and clans(root) == []
    summary = (root / "t_clans_summary.json").read_text()
    assert f'"singleton_families": {len({f for f, _ in rows})}' in summary


# --- exactness and determinism -----------------------------------------------------------


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_random_inputs_match_a_position_set_oracle(tmp_path: Path, seed: int) -> None:
    rows, lengths = random_case(seed)
    members = write_members(tmp_path / "m.tsv", rows)
    fasta = write_fasta(tmp_path / "db.fa", lengths)
    flags = {"min_edge_jaccard": Fraction(1, 20), "clan_jaccard": Fraction(3, 10)}
    code, root = run(
        tmp_path, [members], "-f", str(fasta), "--min_edge_jaccard", "0.05", "--clan_jaccard", "0.3"
    )
    assert code == 0
    expected_edges, expected_clans = oracle(rows, lengths, **flags)
    assert (root / "t_clan_edges.csv").read_text().splitlines()[1:] == expected_edges
    assert clans(root) == expected_clans
    assert len(expected_edges) > 20 and expected_clans  # the case is not vacuous


def test_outputs_do_not_depend_on_execution_settings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rows, lengths = random_case(7)
    fasta = str(write_fasta(tmp_path / "db.fa", lengths))
    members = write_members(tmp_path / "m.tsv", rows)
    _, reference = run(tmp_path, [members], "-f", fasta, out="reference")
    expected = outputs(reference)

    shuffled = rows[::-1]
    half = len(shuffled) // 2
    first = write_members(tmp_path / "a.tsv", shuffled[:half])
    second = tmp_path / "b.tsv.gz"
    with gzip.open(second, "wt") as handle:
        handle.writelines(f"{f}\t{m}\n" for f, m in shuffled[half:])
    gz_fasta = tmp_path / "db.fa.gz"
    gz_fasta.write_bytes(gzip.compress(Path(fasta).read_bytes()))

    tiny_sizes(monkeypatch)
    runs = [
        ([members], fasta, "1"),
        ([first, second], str(gz_fasta), "1"),
        ([first, second], fasta, "4"),
    ]
    for number, (inputs, database, cpus) in enumerate(runs):
        code, root = run(tmp_path, inputs, "-f", database, "-p", cpus, out=f"run{number}")
        assert code == 0
        assert outputs(root) == expected


def test_outputs_do_not_depend_on_the_hash_seed(tmp_path: Path) -> None:
    rows, lengths = random_case(11)
    members = write_members(tmp_path / "m.tsv", rows)
    fasta = write_fasta(tmp_path / "db.fa", lengths)
    results = []
    for seed in ("0", "1"):
        environment = {**os.environ, "PYTHONHASHSEED": seed}
        root = tmp_path / f"seed{seed}"
        command = [
            sys.executable,
            "-m",
            "mgnifam",
            "compute_clans",
            "-i",
            str(members),
            "-f",
            str(fasta),
            "-n",
            "t",
            "--output_dir",
            str(root),
            "-p",
            "2",
        ]
        subprocess.run(command, check=True, env=environment)
        results.append(outputs(root))
    assert results[0] == results[1]


def test_work_is_bounded_and_split_anywhere(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    hub = [(f"F{f:02d}", f"H/{1 + f}-{60 + f}") for f in range(12)]
    hub += [(f"F{f:02d}", f"H/{100 + f}-{120 + f}") for f in range(12)]
    interleaved = [("X", f"L/{10 * k + 1}-{10 * k + 5}") for k in range(300)]
    interleaved += [("Y", f"L/{10 * k + 4}-{10 * k + 8}") for k in range(300)]
    spread = [(f"S{n % 3}", f"Z{n}/1-9") for n in range(200)]
    members = write_members(tmp_path / "m.tsv", hub + interleaved + spread)
    fasta = write_fasta(tmp_path / "db.fa", {f"Z{n}": 10 for n in range(500)})
    _, reference = run(tmp_path, [members], out="reference")
    log = (reference / "t_clans.log").read_text()
    # X and Y overlap 2 residues per period: 300 overlapping pairs, not 300 * 300.
    assert "overlap_emissions=" in log
    overlap = int(log.split("overlap_emissions=")[1].split()[0])
    assert overlap < 1000
    assert edges(reference)[("X", "Y")][2] == f"{600 / 2400:.6f}"

    # Two buckets only, so every bucket is over budget and is split by hash range, and
    # pair batches of 3 split inside one protein's pair rows.
    monkeypatch.setattr(cc, "BUCKET_BITS", 1)
    tiny_sizes(monkeypatch, unit_cost=cc.MEMBER_ROW_COST * 610)
    code, root = run(tmp_path, [members], out="tiny")
    assert code == 0 and outputs(root) == outputs(reference)
    assert int((root / "t_clans.log").read_text().split("units=")[1].split()[0]) > 2

    # A FASTA much larger than the members is split by its own row counts.
    bare = write_members(tmp_path / "bare.tsv", [("A", "Z1"), ("B", "Z1/2-5")])
    tiny_sizes(monkeypatch, unit_cost=cc.LENGTH_ROW_COST * 40)
    code, root = run(tmp_path, [bare], "-f", str(fasta), out="fasta")
    assert code == 0 and edges(root)[("A", "B")][2] == "0.400000"

    # One protein over budget cannot be split.
    tiny_sizes(monkeypatch, unit_cost=cc.MEMBER_ROW_COST * 10)
    code, _ = run(tmp_path, [members], out="over")
    assert code == 1


def test_spreading_work_over_workers_never_fails_a_hub(tmp_path: Path) -> None:
    # One protein holds nearly every row; ~4 units per worker would be far smaller than it.
    rows = [(f"F{f:02d}", f"H/{1 + f}-{50 + f}") for f in range(60)] + [("F00", "Q/1-9")]
    members = write_members(tmp_path / "m.tsv", rows)
    _, one = run(tmp_path, [members], out="one")
    code, many = run(tmp_path, [members], "-p", "8", out="many")
    assert code == 0 and outputs(many) == outputs(one)


def test_missing_lengths_survive_a_worker_boundary() -> None:
    import pickle

    error = pickle.loads(pickle.dumps(cc.MissingLengths(3, [(1, 2)])))
    assert (error.count, error.keys) == (3, [(1, 2)])


# --- validation and the output contract --------------------------------------------------


@pytest.mark.parametrize(
    "line",
    [
        "A\tP/1-10\tx",
        "A",
        "",
        "../x\tP/1-10",
        "A,B\tP/1-10",
        "A\tP/10-1",
        "A\tP/0-5",
        f"A\tP/1-{2**31}",
    ],
)
def test_invalid_input_leaves_the_output_tree_alone(tmp_path: Path, line: str) -> None:
    good = write_members(tmp_path / "good.tsv", [("A", "P/1-10"), ("B", "P/1-10")])
    _, root = run(tmp_path, [good])
    before = outputs(root)
    bad = tmp_path / "bad.tsv"
    # The bad row comes last, so it is found after runs were written to scratch.
    bad.write_text("A\tP/1-10\nB\tP/2-9\n" + line + "\n")
    parent = tmp_path / "scratch"
    parent.mkdir()
    (parent / "keep").write_text("mine")
    code, _ = run(tmp_path, [bad], "--tmp_dir", str(parent))
    assert code == 2
    assert outputs(root) == before
    assert [p.name for p in parent.iterdir()] == ["keep"]


def test_largest_coordinate_is_accepted(tmp_path: Path) -> None:
    rows = [("A", f"P/1-{2**31 - 1}"), ("B", f"P/{2**31 - 10}-{2**31 - 1}")]
    code, root = run(tmp_path, [write_members(tmp_path / "m.tsv", rows)])
    assert code == 0
    assert edges(root)[("A", "B")][3] == "1.000000"


def test_disjoint_and_adjacent_intervals() -> None:
    i64 = lambda *v: np.array(v, np.int64)  # noqa: E731
    _, _, start, end = cc.merge_intervals(
        i64(0, 0, 0, 0), i64(0, 0, 0, 1), i64(1, 6, 10, 1), i64(5, 8, 15, 3)
    )
    assert start.tolist() == [1, 10, 1] and end.tolist() == [8, 15, 3]
    row, offset = cc.virtual_slice(i64(3, 0, 2), np.cumsum(i64(3, 0, 2)), 2, 5)
    assert row.tolist() == [0, 2, 2] and offset.tolist() == [2, 0, 1]


@pytest.mark.parametrize("optimized", [False, True])
def test_overflow_guard_is_explicit(tmp_path: Path, optimized: bool) -> None:
    # 2,200 maximal regions exceed the bound that keeps threshold products in int64.
    rows = [("A", f"P{n}/1-{2**31 - 1}") for n in range(2200)] + [("B", "P0/1-10")]
    members = write_members(tmp_path / "m.tsv", rows)
    flags = ["-O"] if optimized else []
    command = [
        sys.executable,
        *flags,
        "-m",
        "mgnifam",
        "compute_clans",
        "-i",
        str(members),
        "--output_dir",
        str(tmp_path / "out"),
        "-p",
        "1",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 2 and "exact arithmetic" in result.stderr
    assert not (tmp_path / "out").exists()


def test_output_aliases_are_rejected(tmp_path: Path) -> None:
    members = write_members(tmp_path / "m.tsv", [("A", "P/1-10"), ("B", "P/1-10")])
    _, root = run(tmp_path, [members])
    before = outputs(root)
    edges_path, log_path = root / "t_clan_edges.csv", root / "t_clans.log"
    edges_path.unlink()
    edges_path.symlink_to(root / "t_clans_summary.json")
    assert run(tmp_path, [members])[0] == 2
    edges_path.unlink()
    os.link(members, edges_path)
    assert run(tmp_path, [members])[0] == 2
    edges_path.unlink()
    edges_path.write_bytes(before[0])
    log_path.unlink()
    os.link(root / "t_clans_summary.json", log_path)
    assert run(tmp_path, [members])[0] == 2
    assert outputs(root) == before
    assert run(tmp_path, [members], "--tmp_dir", str(root))[0] == 2


def test_empty_fasta_record_is_rejected(tmp_path: Path) -> None:
    members = write_members(tmp_path / "m.tsv", [("A", "P"), ("B", "P/1-1")])
    fasta = tmp_path / "db.fa"
    fasta.write_text(">P\n>Q\nAAA\n")
    assert run(tmp_path, [members], "-f", str(fasta))[0] == 2


def test_failure_while_writing_leaves_no_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    members = write_members(tmp_path / "m.tsv", [("A", "P/1-10"), ("B", "P/1-10")])
    _, root = run(tmp_path, [members])
    assert (root / "t_clans_summary.json").exists()

    def fail(*_: object) -> None:
        raise OSError("disk full")

    monkeypatch.setattr(cc, "build_clans", fail)
    assert run(tmp_path, [members])[0] == 1
    assert not (root / "t_clans_summary.json").exists()


def test_dispatcher_reaches_the_subcommand(capsys) -> None:
    with pytest.raises(SystemExit) as stop:
        cli.main(["compute_clans", "--help"])
    assert stop.value.code == 0
    assert "usage: mgnifam compute_clans" in capsys.readouterr().out
