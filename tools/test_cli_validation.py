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


def expect_runtime_rejected(
    binary: Path, arguments: list[str], *required_fragments: str
) -> None:
    completed = run(binary, *arguments)
    combined = output(completed)
    if completed.returncode == 0:
        raise AssertionError(
            f"command unexpectedly succeeded: {arguments!r}\n{combined}"
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
    expect_accepted(binary, ["-dham-segment-map", "metadata.txt", "--version"], "TITAN ")
    expect_accepted(binary, ["--dham-segment-map", "metadata.txt", "--version"], "TITAN ")
    expect_runtime_rejected(
        binary,
        ["-dham", "-dham-segment-map", "metadata.txt", "-dham-ener", "energy.txt"],
        "-dham-segment-map and -dham-ener are mutually exclusive",
    )
    expect_accepted(binary, ["-1pdb", "--version"], "TITAN ")
    expect_accepted(binary, ["-gbsa", "--version"], "TITAN ")
    expect_accepted(
        binary,
        ["-gbsa-sasa-points", "384", "-pbsa-sasa-points", "1536", "--version"],
        "TITAN ",
    )
    expect_runtime_rejected(
        binary,
        [
            "-pathway-file",
            "known-route.tpath",
            "-pathway-probe-lengths",
            "1.0,0.5,0.25",
            "-pathway-clearance-mode",
            "swept-linear-center",
            "-pull-direct",
            "-pull-opt",
            "corridor-grid",
        ],
        "swept-linear-center clearance requires -pbcn or -pbctopo",
    )
    expect_accepted(
        binary,
        [
            "-pathway-file",
            "known-route.tpath",
            "-pathway-route-id",
            "route-A",
            "-pathway-corridor-spacing",
            "0.5",
            "-pathway-corridor-radial-samples",
            "16",
            "-pathway-recovery-enter",
            "4.0",
            "-pathway-recovery-exit",
            "3.0",
            "-pathway-recovery-max-turn-angle",
            "60.0",
            "-pathway-energy-key-epsilon",
            "0.05",
            "-pathway-completion-min-separation",
            "4.0",
            "-pathway-min-clearance",
            "1.5",
            "-pathway-probe-lengths",
            "1.0,0.5,0.25",
            "-pathway-navigation",
            "-pathway-navigation-tube-weight",
            "1.0",
            "-pathway-navigation-backtrack-weight",
            "2.0",
            "-pathway-navigation-cost-epsilon",
            "0.1",
            "-pathway-navigation-seed-radial-samples",
            "8",
            "-pathway-navigation-seed-amplitude",
            "0.35",
            "-pathway-clearance-mode",
            "swept-linear-center",
            "-pbcn",
            "-pull-opt",
            "corridor-grid",
            "--version",
        ],
        "TITAN ",
    )
    expect_runtime_rejected(
        binary,
        [
            "-pathway-file",
            "known-route.tpath",
            "-pathway-probe-lengths",
            "1.0,not-a-number,0.25",
            "-pull-opt",
            "corridor-grid",
            "--version",
        ],
        "-pathway-probe-lengths expects positive finite distances",
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

        pathway = write_arg_file(
            directory,
            "pathway.in",
            "-pathway-file=routes/known.tpath\n"
            "-pathway-route-id=route-A\n"
            "-pathway-previous-status=states/previous.txt\n"
            "-pathway-continuity-window=1.0\n"
            "-pathway-corridor-spacing=0.5\n"
            "-pathway-corridor-radial-samples=16\n"
            "-pathway-recovery-enter=4.0\n"
            "-pathway-recovery-exit=3.0\n"
            "-pathway-recovery-max-turn-angle=60.0\n"
            "-pathway-energy-key-epsilon=0.05\n"
            "-pathway-completion-min-separation=4.0\n"
            "-pathway-min-clearance=1.5\n"
            "-pathway-probe-lengths=1.0,0.5,0.25\n"
            "-pathway-navigation\n"
            "-pathway-navigation-tube-weight=1.0\n"
            "-pathway-navigation-backtrack-weight=2.0\n"
            "-pathway-navigation-cost-epsilon=0.1\n"
            "-pathway-navigation-seed-radial-samples=8\n"
            "-pathway-navigation-seed-amplitude=0.35\n"
            "-pathway-clearance-mode=swept-linear-center\n"
            "-pbcn\n"
            "-pull-opt=corridor-grid\n"
            "--version=true\n",
        )
        expect_accepted(binary, ["-a", str(pathway)], "TITAN ")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    verify(args.binary)
    print("PASS test_cli_validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
