import gzip
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from mgnifam import generate_families as gf


@pytest.fixture(scope="session")
def fixture_directory() -> Path:
    return Path(__file__).parent / "fixtures"


def _decompress(source: Path, destination: Path) -> Path:
    with gzip.open(source, "rb") as compressed, destination.open("wb") as uncompressed:
        shutil.copyfileobj(compressed, uncompressed)
    return destination


@pytest.fixture(scope="session")
def small_fasta(tmp_path_factory: pytest.TempPathFactory, fixture_directory: Path) -> Path:
    temporary = tmp_path_factory.mktemp("small-fasta")
    return _decompress(
        fixture_directory / "mgnifams_input_small.fa.gz",
        temporary / "mgnifams_input_small.fa",
    )


@pytest.fixture(scope="session")
def extra_fasta(tmp_path_factory: pytest.TempPathFactory, fixture_directory: Path) -> Path:
    temporary = tmp_path_factory.mktemp("extra-fasta")
    return _decompress(
        fixture_directory / "mgnifams_extra.fa.gz",
        temporary / "mgnifams_extra.fa",
    )


@pytest.fixture(scope="session")
def v2_fasta(tmp_path_factory: pytest.TempPathFactory, fixture_directory: Path) -> Path:
    temporary = tmp_path_factory.mktemp("v2-fasta")
    return _decompress(
        fixture_directory / "mgnifams_v2.fa.gz",
        temporary / "mgnifams_v2.fa",
    )


@pytest.fixture
def index_reads(monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Watch the SSI index across a run: a name read twice within one batch fails at once.

    `clears` records the memo size at each `IndexedSequences.clear`, so a caller can check
    it ran once per batch; `instances` lets the caller check the memo ended empty; `gets`
    and `reads` together show the memo was actually hit, so a pass is not vacuous.
    """
    state = SimpleNamespace(instances=[], clears=[], gets=0, reads=0)
    batch_reads: set[str] = set()
    original_fetch = gf.fetch_indexed_sequence
    original_init = gf.IndexedSequences.__init__
    original_get = gf.IndexedSequences.get
    original_clear = gf.IndexedSequences.clear

    def fetch(indexed: object, name: str) -> object:
        assert name not in batch_reads, f"{name} read from the index twice in one batch"
        fetched = original_fetch(indexed, name)  # type: ignore[arg-type]
        batch_reads.add(name)
        state.reads += 1
        return fetched

    def init(self: gf.IndexedSequences, sequence_file: object) -> None:
        original_init(self, sequence_file)  # type: ignore[arg-type]
        state.instances.append(self)

    def get(self: gf.IndexedSequences, name: str, **kwargs: bool) -> object:
        state.gets += 1
        return original_get(self, name, **kwargs)

    def clear(self: gf.IndexedSequences) -> None:
        state.clears.append(len(self.cache))
        batch_reads.clear()
        original_clear(self)

    monkeypatch.setattr(gf, "fetch_indexed_sequence", fetch)
    monkeypatch.setattr(gf.IndexedSequences, "__init__", init)
    monkeypatch.setattr(gf.IndexedSequences, "get", get)
    monkeypatch.setattr(gf.IndexedSequences, "clear", clear)
    return state
