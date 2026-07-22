"""mgnifam -- protein family generation over very large sequence databases."""

from importlib.metadata import version

# Read from the installed distribution rather than restated here, so `pyproject.toml` is
# the one place a release number is written. A literal would be a second source of truth
# that nothing checks, and the two drift silently towards whichever release forgot one of
# them. Safe under a src layout: the package cannot be imported without being installed,
# so its metadata is always present.
__version__ = version("mgnifam")

__all__ = ["__version__"]
