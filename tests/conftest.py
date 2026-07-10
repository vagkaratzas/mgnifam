import gzip
import shutil
from pathlib import Path

import pytest


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
