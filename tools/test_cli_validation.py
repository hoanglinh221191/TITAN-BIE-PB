#!/usr/bin/env python3
"""End-to-end checks for TITAN terminal and arg.in syntax validation."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile


def run(binary: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(binary.resolve()), *arguments],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def output(completed: subprocess.CompletedProcess[str]) -> str:
    return completed.stdout + completed.stderr


def expect_rejected(
    binary: Path, arguments: list[str], *required_fragments: str
) -> None:
    completed = run(binary, *arguments)
    combined = output(completed)
    if completed.returncode != 3:
        raise AssertionError(
            f"command returned {completed.returncode}, expected syntax status 3: "
            f"{arguments!r}\n{combined}"
        )
    missing = [fragment for fragment in required_fragments if fragment not in combined]
    if missing:
        raise AssertionError(
            f"command {arguments!r} did not report {missing!r}\n{combined}"
        )


def expect_accepted(
    binary: Path, arguments: list[str], *required_fragments: str
) -> None:
    completed = run(binary, *arguments)
    combined = output(completed)
    if completed.returncode != 0:
        raise AssertionError(
            f"command failed with {completed.returncode}: {arguments!r}\n{combined}"
        )
    missing = [fragment for fragment in required_fragments if fragment not in combined]
    if missing:
        raise AssertionError(
            f"command {arguments!r} did not report {missing!r}\n{combined}"
        )


def write_arg_file(directory: Path, name: str, content: str) -> Path:
    path = directory / name
    path.write_text(content, encoding="utf-8", newline="\n")
    return path


def verify(binary: Path) -> None:
    expect_accepted(binary, ["-1pdb", "--version"], "TITAN ")
    expect_accepted(binary, ["-gbsa", "--version"], "TITAN ")
    expect_accepted(
        binary,
        ["-gbsa-sasa-points", "384", "-pbsa-sasa-points", "1536", "--version"],
        "TITAN ",
    )
    expect_rejected(binary, ["-1pdb", "tru", "--version"], "CLI010", "tru")
    expect_rejected(binary, ["-gbsa", "maybe", "--version"], "CLI010", "maybe")
    expect_rejected(
        binary,
        ["-gbsa-sasa-points", "384", "orphan", "--version"],
        "CLI010",
        "orphan",
    )
    expect_rejected(
        binary,
        ["-pbe-bie-opening-theta=0.6"],
        "CLI001",
        "reserved for arg.in",
    )
    expect_rejected(
        binary,
        ["-sg=", "1", "2-3", "4", "--version"],
        "CLI001",
        "reserved for arg.in",
    )
    expect_rejected(
        binary,
        ["-pbe-bie-opning-theta", "0.6"],
        "CLI002",
        "-pbe-bie-opening-theta",
    )
    expect_rejected(binary, ["-emap", "on"], "CLI006", "-edecomp")
    expect_rejected(binary, ["-pbctopo", "-pbc"], "CLI007")
    expect_rejected(
        binary,
        ["-pbe-binding-cycle", "pqr-components"],
        "CLI008",
        "complex, receptor, and ligand",
    )
    expect_accepted(
        binary,
        [
            "-sg",
            "1",
            "2-3",
            "4",
            "-pbsa-components",
            "ligand-receptor-complex",
            "--version",
        ],
        "TITAN ",
    )
    expect_rejected(
        binary,
        ["-sg", "=1", "2=3", "4", "--version"],
        "Invalid -sg expression",
    )
    expect_rejected(
        binary,
        ["-pbsa-components", "ligand-recpetor-complex", "--version"],
        "did you mean 'receptor'",
    )

    with tempfile.TemporaryDirectory(prefix="titan-cli-validation-") as temporary:
        directory = Path(temporary)
        missing_equals = write_arg_file(
            directory, "missing_equals.in", "-pbe-bie-opening-theta 0.6\n"
        )
        expect_rejected(
            binary, ["-a", str(missing_equals)], "CLI005", "'=' sign is mandatory"
        )

        implicit_switch = write_arg_file(
            directory, "implicit_switch.in", "-1pdb\n--version\n"
        )
        expect_accepted(binary, ["-a", str(implicit_switch)], "TITAN ")

        explicit_switch = write_arg_file(
            directory, "explicit_switch.in", "-1pdb=true\n--version=true\n"
        )
        expect_accepted(binary, ["-a", str(explicit_switch)], "TITAN ")

        implicit_boolean = write_arg_file(
            directory, "implicit_boolean.in", "-gbsa\n--version\n"
        )
        expect_accepted(binary, ["-a", str(implicit_boolean)], "TITAN ")

        sasa_points = write_arg_file(
            directory,
            "sasa_points.in",
            "-gbsa-sasa-points=384\n-pbsa-sasa-points=1536\n--version\n",
        )
        expect_accepted(binary, ["-a", str(sasa_points)], "TITAN ")

        bad_spacing = write_arg_file(
            directory, "bad_spacing.in", "-pbe-bie-opening-theta =0.6\n"
        )
        expect_rejected(
            binary,
            ["-a", str(bad_spacing)],
            "CLI005",
            "without whitespace before '='",
        )

        typo = write_arg_file(
            directory, "typo.in", "-pbe-bie-opning-theta=0.6\n"
        )
        expect_rejected(
            binary, ["-a", str(typo)], "CLI002", "-pbe-bie-opening-theta"
        )

        conflict = write_arg_file(
            directory,
            "conflict.in",
            "-pbctopo=true\n-pbc=true\n--version=true\n",
        )
        expect_rejected(binary, ["-a", str(conflict)], "CLI007", "-pbctopo")

        special = write_arg_file(
            directory,
            "special.in",
            "-sg= 1 2-3 4\n"
            "-pbsa-components=complex,ligand,receptor\n"
            "--version=true\n",
        )
        expect_accepted(binary, ["-a", str(special)], "TITAN ")

        complete_pbe = write_arg_file(
            directory,
            "complete_pbe.in",
            "-pbe-binding-cycle=pqr-components\n"
            "-pbe-test-pqr=complex.pqr\n"
            "-pbe-binding-receptor-pqr=receptor.pqr\n"
            "-pbe-binding-ligand-pqr=ligand.pqr\n"
            "--version=true\n",
        )
        expect_accepted(binary, ["-a", str(complete_pbe)], "TITAN ")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    verify(args.binary)
    print("PASS test_cli_validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
