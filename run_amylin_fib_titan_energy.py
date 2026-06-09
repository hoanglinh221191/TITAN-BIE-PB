#!/usr/bin/env python
r"""Run TITAN protein-ligand elec/vdW energy for amylin fib systems.

This workflow targets:
  D:\Data_MD\amylin\4497_fib
  D:\Data_MD\amylin\4631_fib

Each system has a protein-containing PDB and a ligand-only ACPYPE topology
(`UNL_GMX.top`), not a shared GROMACS topol.top for the full system. TITAN is
therefore run with the protein parsed from the PDB/AMBER14SB force field and
the ligand parsed from `-topo=UNL_GMX.top`.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


OUTPUT_DIR_NAME = "titan_energy_md2_9_multi_xdr_amber14sb_ligtopo"


@dataclass(frozen=True)
class SystemConfig:
    label: str
    path: Path
    structure: str = "conf1.pdb"
    index: str = "index.ndx"
    ligand_topology: str = "UNL_GMX.top"
    protein_group: str = "Protein"
    ligand_group: str = "UNL"


def parse_index_groups(path: Path) -> dict[str, tuple[int, int]]:
    groups: dict[str, tuple[int, int]] = {}
    current: str | None = None
    counts: dict[str, int] = {}
    order: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            header = re.match(r"^\s*\[\s*(.*?)\s*\]\s*$", raw)
            if header:
                current = header.group(1)
                counts[current] = 0
                order.append(current)
                continue
            if current is None:
                continue
            counts[current] += sum(1 for token in raw.split() if token.isdigit())
    for idx, name in enumerate(order):
        groups[name] = (idx, counts[name])
    return groups


def count_topology_atoms(path: Path) -> int:
    in_atoms = False
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(";"):
                continue
            section = re.match(r"^\[\s*(.*?)\s*\]$", line)
            if section:
                in_atoms = section.group(1).strip().lower() == "atoms"
                continue
            if in_atoms and re.match(r"^\d+\s+", line):
                count += 1
    return count


def discover_trajectories(system_dir: Path) -> list[Path]:
    trajectories = []
    for path in system_dir.glob("md*_1.xtc"):
        match = re.fullmatch(r"md(\d+)_1\.xtc", path.name)
        if match and 2 <= int(match.group(1)) <= 9:
            trajectories.append(path)
    trajectories.sort(key=lambda p: int(re.fullmatch(r"md(\d+)_1\.xtc", p.name).group(1)))
    if len(trajectories) != 8:
        names = ", ".join(path.name for path in trajectories) or "none"
        raise FileNotFoundError(f"{system_dir}: expected md2_1.xtc..md9_1.xtc, found {names}")
    return trajectories


def write_input_list(system_dir: Path, trajectories: list[Path]) -> Path:
    input_path = system_dir / "input_md2_9.txt"
    input_path.write_text(
        "".join(f"{path.name}\n" for path in trajectories),
        encoding="utf-8",
    )
    return input_path


def write_arg_in(
    system: SystemConfig,
    protein_group_index: int,
    ligand_group_index: int,
    input_list: Path,
    begin_ps: float | None,
    end_ps: float | None,
    threads: int,
) -> Path:
    lines = [
        f"# TITAN native multiple-XDR protein-ligand interaction energy for {system.label}",
        "# Protein parsed from conf1.pdb with AMBER14SB; ligand parsed from UNL_GMX.top.",
        f"-s={system.structure}",
        f"-i={input_list.name}",
        f"-n={system.index}",
        f"-topo={system.ligand_topology}",
        "-ff=AMBER14SB",
        f"-g={protein_group_index}-{ligand_group_index}",
        "-mol-1-type=protein",
        "-mol-2-type=ligand",
    ]
    if begin_ps is not None:
        lines.append(f"-b={begin_ps:g}")
    if end_ps is not None:
        lines.append(f"-e={end_ps:g}")
    lines.extend(
        [
            "-pbcn",
            "-no-dssp",
            "-output-manifest=true",
            "-perf-report=summary",
            f"-nt={threads}",
            "",
        ]
    )
    arg_path = system.path / "arg.in"
    arg_path.write_text("\n".join(lines), encoding="utf-8")
    return arg_path


def rotate_existing_output(target: Path) -> None:
    if not target.exists():
        return
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = target.with_name(f"{target.name}_previous_{stamp}")
    target.rename(backup)


def move_titan_outputs(system_dir: Path, target: Path, arg_path: Path, input_path: Path) -> None:
    data_md = system_dir / "data_MD"
    target.mkdir(parents=True, exist_ok=True)
    for path in list(data_md.iterdir()):
        if path.is_file():
            shutil.move(str(path), str(target / path.name))
    root_outputs = [
        "cal.txt",
        "cal.xvg",
        "electrostatic.dat",
        "vanderwaals.dat",
        "total_energy.dat",
        "log.txt",
        "titan_outputs_manifest.csv",
        "titan_performance_modules.csv",
        "titan_performance_summary.dat",
        "titan_run_config_summary.dat",
        "fflib_assignment_report.tsv",
        "debug_helper.txt",
        "debug_topo_lg.txt",
    ]
    for name in root_outputs:
        path = system_dir / name
        if path.exists():
            shutil.move(str(path), str(target / path.name))
    shutil.copy2(arg_path, target / "arg.in")
    shutil.copy2(input_path, target / input_path.name)


def run_system(system: SystemConfig, titan_source: Path, begin_ps: float | None, end_ps: float | None, threads: int) -> None:
    required = [system.path / system.structure, system.path / system.index, system.path / system.ligand_topology]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

    groups = parse_index_groups(system.path / system.index)
    if system.protein_group not in groups:
        raise KeyError(f"{system.path}: missing index group {system.protein_group!r}")
    if system.ligand_group not in groups:
        raise KeyError(f"{system.path}: missing index group {system.ligand_group!r}")

    protein_idx, protein_atoms = groups[system.protein_group]
    ligand_idx, ligand_atoms = groups[system.ligand_group]
    topology_atoms = count_topology_atoms(system.path / system.ligand_topology)
    if topology_atoms != ligand_atoms:
        raise ValueError(
            f"{system.path}: ligand topology has {topology_atoms} atoms, "
            f"but index group {system.ligand_group} has {ligand_atoms}"
        )

    trajectories = discover_trajectories(system.path)
    titan_dest = system.path / "titan.exe"
    shutil.copy2(titan_source, titan_dest)
    ff_source = titan_source.parent / "lib_ff"
    if not ff_source.exists():
        raise FileNotFoundError(f"Missing TITAN force-field runtime directory: {ff_source}")
    shutil.copytree(ff_source, system.path / "lib_ff", dirs_exist_ok=True)
    input_path = write_input_list(system.path, trajectories)
    arg_path = write_arg_in(system, protein_idx, ligand_idx, input_path, begin_ps, end_ps, threads)

    data_md = system.path / "data_MD"
    data_md.mkdir(exist_ok=True)
    target = data_md / OUTPUT_DIR_NAME
    rotate_existing_output(target)

    stdout_path = target / "titan_stdout.log"
    stderr_path = target / "titan_stderr.log"
    target.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8", errors="replace") as stdout, stderr_path.open(
        "w", encoding="utf-8", errors="replace"
    ) as stderr:
        result = subprocess.run(
            [str(titan_dest), "-a", "arg.in"],
            cwd=system.path,
            stdout=stdout,
            stderr=stderr,
            text=True,
            check=False,
        )
    move_titan_outputs(system.path, target, arg_path, input_path)
    if stdout_path.exists() and not (target / stdout_path.name).exists():
        shutil.move(str(stdout_path), str(target / stdout_path.name))
    if stderr_path.exists() and not (target / stderr_path.name).exists():
        shutil.move(str(stderr_path), str(target / stderr_path.name))
    if result.returncode != 0:
        raise RuntimeError(f"{system.label}: TITAN exited with code {result.returncode}; see {target}")

    print(
        f"{system.label}: Protein group {protein_idx} ({protein_atoms} atoms), "
        f"UNL group {ligand_idx} ({ligand_atoms} atoms), output {target}"
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run TITAN elec/vdW calculations for amylin fib trajectories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", type=Path, default=Path(r"D:\Data_MD\amylin"))
    parser.add_argument(
        "--systems",
        nargs=2,
        default=["4497_fib", "4631_fib"],
        help="System directories under --root.",
    )
    parser.add_argument(
        "--titan",
        type=Path,
        default=Path(r"C:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw\titan.exe"),
        help="Source titan.exe to copy into each system directory.",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--begin-ps", type=float, default=None, help="Optional trajectory begin time in ps.")
    parser.add_argument("--end-ps", type=float, default=None, help="Optional trajectory end time in ps.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    systems = [SystemConfig(label=name, path=args.root / name) for name in args.systems]
    for system in systems:
        run_system(system, args.titan, args.begin_ps, args.end_ps, args.threads)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
