#!/usr/bin/env python3
"""Validate, assemble, and checksum deterministic TITAN release artifacts."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
from typing import Mapping
import zipfile

import generate_build_info as build_info


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_release_tag(repo_root: Path, tag: str) -> str:
    version, _ = build_info.read_product_version(repo_root / "VERSION", {})
    expected = f"titan-v{version}"
    if tag != expected:
        raise ValueError(f"release tag {tag!r} must exactly match {expected!r}")
    head = build_info.run_git(repo_root, "rev-parse", "HEAD")
    tag_commit = build_info.run_git(repo_root, "rev-list", "-n", "1", tag)
    if head is None or tag_commit is None or head != tag_commit:
        raise ValueError(f"release tag {tag!r} does not point at HEAD")
    status = build_info.run_git(
        repo_root, "status", "--porcelain=v1", "--untracked-files=normal"
    )
    if status is None:
        raise ValueError("release validation requires a Git worktree")
    if status:
        raise ValueError("release validation requires a clean Git worktree")
    return version


def normalized_release_files(
    repo_root: Path,
    binary: Path,
    manifest: Mapping[str, object],
) -> dict[str, tuple[bytes, int]]:
    version = str(manifest["version"])
    repository_version, _ = build_info.read_product_version(
        repo_root / "VERSION", {}
    )
    if version != repository_version:
        raise ValueError(
            "release manifest version does not match the repository VERSION file"
        )
    if manifest.get("full_version") != version or not manifest.get("release_build"):
        raise ValueError("release artifact requires a clean, exactly tagged build")
    if manifest.get("git_dirty") is not False:
        raise ValueError("release artifact manifest must record a clean Git state")
    expected_tag = f"titan-v{version}"
    if manifest.get("git_tag") != expected_tag:
        raise ValueError(
            f"release artifact manifest must record exact tag {expected_tag!r}"
        )
    commit = str(manifest.get("git_commit", ""))
    if not build_info.GIT_SHA_RE.fullmatch(commit):
        raise ValueError("release artifact manifest must contain a valid Git commit")
    if manifest.get("profile") != "portable":
        raise ValueError("release artifact must use the portable build profile")
    if not binary.is_file():
        raise ValueError(f"release binary not found: {binary}")

    manifest_bytes = (
        json.dumps(dict(manifest), indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    files: dict[str, tuple[bytes, int]] = {
        binary.name: (binary.read_bytes(), 0o755),
        "LICENSE": ((repo_root / "LICENSE").read_bytes(), 0o644),
        "VERSION": ((repo_root / "VERSION").read_bytes(), 0o644),
        "titan_build_manifest.json": (manifest_bytes, 0o644),
    }
    checksum_lines = [
        f"{sha256_bytes(content)}  {name}"
        for name, (content, _) in sorted(files.items())
    ]
    files["SHA256SUMS"] = (("\n".join(checksum_lines) + "\n").encode("utf-8"), 0o644)
    return files


def write_zip(
    output: Path,
    root_name: str,
    files: Mapping[str, tuple[bytes, int]],
    epoch: int,
) -> None:
    timestamp = dt.datetime.fromtimestamp(
        max(epoch, 315532800), tz=dt.timezone.utc
    )
    date_time = (timestamp.year, timestamp.month, timestamp.day,
                 timestamp.hour, timestamp.minute, timestamp.second)
    with zipfile.ZipFile(
        output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for name, (content, mode) in sorted(files.items()):
            entry = zipfile.ZipInfo(f"{root_name}/{name}", date_time=date_time)
            entry.create_system = 3
            entry.external_attr = (mode & 0xFFFF) << 16
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, content, compress_type=zipfile.ZIP_DEFLATED,
                             compresslevel=9)


def write_tar_gz(
    output: Path,
    root_name: str,
    files: Mapping[str, tuple[bytes, int]],
    epoch: int,
) -> None:
    with output.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=epoch) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for name, (content, mode) in sorted(files.items()):
                    entry = tarfile.TarInfo(f"{root_name}/{name}")
                    entry.size = len(content)
                    entry.mode = mode
                    entry.mtime = epoch
                    entry.uid = 0
                    entry.gid = 0
                    entry.uname = ""
                    entry.gname = ""
                    archive.addfile(entry, io.BytesIO(content))


def create_release_archive(
    repo_root: Path,
    binary: Path,
    manifest_path: Path,
    platform: str,
    output_dir: Path,
    archive_format: str,
) -> tuple[Path, Path]:
    if not build_info.SAFE_TOKEN_RE.fullmatch(platform):
        raise ValueError(f"invalid platform token {platform!r}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    version = str(manifest["version"])
    files = normalized_release_files(repo_root, binary, manifest)
    epoch = int(manifest.get("source_date_epoch", 0))
    if epoch <= 0:
        raise ValueError("release manifest must contain a positive source_date_epoch")

    root_name = f"titan-{version}-{platform}"
    suffix = ".zip" if archive_format == "zip" else ".tar.gz"
    output_dir.mkdir(parents=True, exist_ok=True)
    archive_path = output_dir / f"{root_name}{suffix}"
    if archive_format == "zip":
        write_zip(archive_path, root_name, files, epoch)
    elif archive_format == "tar.gz":
        write_tar_gz(archive_path, root_name, files, epoch)
    else:
        raise ValueError(f"unsupported archive format {archive_format!r}")

    checksum_path = archive_path.with_name(archive_path.name + ".sha256")
    checksum_path.write_text(
        f"{sha256_bytes(archive_path.read_bytes())}  {archive_path.name}\n",
        encoding="utf-8",
        newline="\n",
    )
    return archive_path, checksum_path


def verify_binary(binary: Path, manifest_path: Path) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    completed = subprocess.run(
        [str(binary.resolve()), "--version"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    first_line = completed.stdout.splitlines()[0] if completed.stdout else ""
    if first_line != f"TITAN {manifest['full_version']}":
        raise ValueError("release binary version does not match its build manifest")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--validate-tag")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--platform")
    parser.add_argument("--output-dir", type=Path, default=Path("dist"))
    parser.add_argument("--format", choices=("zip", "tar.gz"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    repo_root = args.repo_root.resolve()
    if args.validate_tag:
        version = validate_release_tag(repo_root, args.validate_tag)
        print(f"release tag: PASS ({version})")

    packaging_values = (args.binary, args.manifest, args.platform, args.format)
    if any(value is not None for value in packaging_values):
        if not all(value is not None for value in packaging_values):
            raise ValueError(
                "--binary, --manifest, --platform, and --format are required together"
            )
        verify_binary(args.binary, args.manifest)
        archive, checksum = create_release_archive(
            repo_root,
            args.binary,
            args.manifest,
            args.platform,
            args.output_dir,
            args.format,
        )
        print(archive)
        print(checksum)
    if args.validate_tag is None and not any(value is not None for value in packaging_values):
        raise ValueError("no release operation requested")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"package_release: {error}") from error
