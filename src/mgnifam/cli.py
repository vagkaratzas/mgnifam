"""Command dispatch for the `mgnifam` executable.

Every capability is a subcommand -- `mgnifam generate_families ...` -- so that later
additions such as `remove_redundant` or `merge_families` need only a new entry in
`COMMANDS` and a module exposing a `main(argv)`.

Subcommands own their own argument parsers rather than registering with an
`add_subparsers` group. Each remains directly callable as `module.main(argv)`, which
keeps the tests independent of this dispatcher, and keeps a command's flags from
colliding with any option this layer might grow.
"""

import argparse
import sys
from collections.abc import Callable, Sequence

from mgnifam import __version__, generate_families

COMMANDS: dict[str, Callable[[Sequence[str]], None]] = {
    "generate_families": generate_families.main,
}


def main(argv: Sequence[str] | None = None) -> None:
    """Dispatch to a subcommand. Exits non-zero on an unknown or missing command."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="mgnifam",
        description="Protein family generation over very large sequence databases.",
    )
    parser.add_argument("--version", action="version", version=f"mgnifam {__version__}")
    parser.add_argument("command", choices=sorted(COMMANDS), help="subcommand to run")
    # REMAINDER hands every following token to the subcommand untouched, including
    # --help, so `mgnifam generate_families --help` reaches the subcommand's parser.
    parser.add_argument(
        "arguments", nargs=argparse.REMAINDER, help="arguments for the subcommand"
    )
    namespace = parser.parse_args(arguments)
    COMMANDS[namespace.command](namespace.arguments)


if __name__ == "__main__":
    main()
