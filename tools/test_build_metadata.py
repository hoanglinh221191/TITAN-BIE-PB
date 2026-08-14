#!/usr/bin/env python3
"""Unit and built-artifact checks for TITAN build metadata."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import generate_build_info as build_info
import package_release


class BuildMetadataTests(unittest.TestCase):
    def test_clean_release_uses_plain_semver(self) -> None:
        self.assertEqual(
            build_info.compose_full_version("1.2.3", 42, "abcdef012345", False, True),
            "1.2.3",
        )

    def test_development_version_contains_revision_and_dirty_state(self) -> None:
        self.assertEqual(
            build_info.compose_full_version("1.2.3", 42, "abcdef012345", True, False),
            "1.2.3-dev.42+gabcdef012345.dirty",
        )

    def test_environment_override_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            version_file = root / "VERSION"
            version_file.write_text("2.4.6\n", encoding="utf-8")
            env = {
                "TITAN_GIT_COMMIT": "0123456789abcdef0123456789abcdef01234567",
                "TITAN_BUILD_NUMBER": "77",
                "TITAN_GIT_DIRTY": "0",
                "TITAN_GIT_TAG": "titan-v2.4.6",
                "TITAN_RELEASE_BUILD": "1",
                "TITAN_CI_BUILD_NUMBER": "91",
                "TITAN_CI_RUN_ID": "1234",
                "SOURCE_DATE_EPOCH": "1700000000",
                "GITHUB_ACTIONS": "true",
            }
            info = build_info.collect_build_info(root, version_file, "portable", env)
            self.assertEqual(info["full_version"], "2.4.6")
            self.assertEqual(info["build_number"], 77)
            self.assertEqual(info["ci_build_number"], "91")
            self.assertEqual(info["source_date_utc"], "2023-11-14T22:13:20Z")
            self.assertTrue(info["release_build"])
            self.assertFalse(info["ci_metadata_embedded"])
            header = build_info.render_header(info)
            self.assertIn('#define TITAN_CI_PROVIDER "not-embedded"', header)
            self.assertIn('#define TITAN_CI_BUILD_NUMBER ""', header)

    def test_release_tag_must_match_version(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            version_file = root / "VERSION"
            version_file.write_text("1.0.0\n", encoding="utf-8")
            env = {
                "TITAN_GIT_COMMIT": "abcdef0123456789abcdef0123456789abcdef01",
                "TITAN_BUILD_NUMBER": "1",
                "TITAN_GIT_DIRTY": "0",
                "TITAN_GIT_TAG": "titan-v2.0.0",
                "TITAN_RELEASE_BUILD": "1",
                "SOURCE_DATE_EPOCH": "1",
            }
            with self.assertRaisesRegex(ValueError, "requires exact tag"):
                build_info.collect_build_info(root, version_file, "portable", env)

    def test_release_requires_an_exact_tag(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            version_file = root / "VERSION"
            version_file.write_text("1.0.0\n", encoding="utf-8")
            env = {
                "TITAN_GIT_COMMIT": "abcdef0123456789abcdef0123456789abcdef01",
                "TITAN_BUILD_NUMBER": "1",
                "TITAN_GIT_DIRTY": "0",
                "TITAN_RELEASE_BUILD": "1",
                "SOURCE_DATE_EPOCH": "1",
            }
            with self.assertRaisesRegex(ValueError, "requires exact tag"):
                build_info.collect_build_info(root, version_file, "portable", env)

    def test_write_if_changed_preserves_timestamp(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "generated.h"
            self.assertTrue(build_info.write_if_changed(output, "value\n"))
            first_mtime = output.stat().st_mtime_ns
            time.sleep(0.01)
            self.assertFalse(build_info.write_if_changed(output, "value\n"))
            self.assertEqual(output.stat().st_mtime_ns, first_mtime)

    def test_windows_resource_uses_bounded_numeric_build(self) -> None:
        info = {
            "version_major": 1,
            "version_minor": 2,
            "version_patch": 3,
            "build_number": 999999,
            "release_build": False,
            "git_dirty": True,
            "full_version": "1.2.3-dev.999999+gabcdef.dirty",
            "version": "1.2.3",
            "git_commit_short": "abcdef",
            "ci_build_number": "5",
            "ci_metadata_embedded": True,
        }
        resource = build_info.render_windows_resource(info, "titan.exe")
        self.assertIn("FILEVERSION 1,2,3,65535", resource)
        self.assertIn("FILEFLAGS 0x0aL", resource)
        self.assertIn('VALUE "SpecialBuild", "5\\0"', resource)

    def test_release_archive_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "LICENSE").write_text("license\n", encoding="utf-8")
            (root / "VERSION").write_text("1.2.3\n", encoding="utf-8")
            binary = root / "titan"
            binary.write_bytes(b"binary")
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "version": "1.2.3",
                        "full_version": "1.2.3",
                        "release_build": True,
                        "git_commit": "abcdef0123456789abcdef0123456789abcdef01",
                        "git_dirty": False,
                        "git_tag": "titan-v1.2.3",
                        "profile": "portable",
                        "source_date_epoch": 1700000000,
                    }
                ),
                encoding="utf-8",
            )
            archive, _ = package_release.create_release_archive(
                root, binary, manifest, "linux-x86_64", root / "dist", "tar.gz"
            )
            first = archive.read_bytes()
            archive, _ = package_release.create_release_archive(
                root, binary, manifest, "linux-x86_64", root / "dist", "tar.gz"
            )
            self.assertEqual(first, archive.read_bytes())


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
    output = completed.stdout
    expected_lines = [
        f"TITAN {manifest['full_version']}",
        f"Version: {manifest['version']}",
        f"Build number: {manifest['build_number']}",
        f"Commit: {manifest['git_commit_short']}",
        f"Profile: {manifest['profile']}",
        f"Source date: {manifest['source_date_utc']}",
    ]
    missing = [line for line in expected_lines if line not in output]
    if missing:
        raise AssertionError(
            "binary version output is missing: " + ", ".join(repr(line) for line in missing)
        )


def verify_schema_driver(driver: Path, manifest_path: Path) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="titan-build-schema-") as temporary:
        output = Path(temporary) / "titan_schema.json"
        subprocess.run(
            [str(driver.resolve()), str(output)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        schema = json.loads(output.read_text(encoding="utf-8"))
    if schema.get("schema_version") != 2:
        raise AssertionError("config schema version was not upgraded to 2")
    embedded = schema.get("build_info", {})
    expected = {
        "version": manifest["version"],
        "full_version": manifest["full_version"],
        "build_number": manifest["build_number"],
        "git_commit": manifest["git_commit"],
        "git_dirty": manifest["git_dirty"],
        "profile": manifest["profile"],
        "generated_metadata": True,
    }
    mismatches = {
        key: (embedded.get(key), value)
        for key, value in expected.items()
        if embedded.get(key) != value
    }
    if mismatches:
        raise AssertionError(f"config schema build_info mismatch: {mismatches}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--schema-driver", type=Path)
    args = parser.parse_args()
    if (args.binary is None) != (args.manifest is None):
        parser.error("--binary and --manifest must be supplied together")

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(BuildMetadataTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        return 1
    if args.binary is not None and args.manifest is not None:
        verify_binary(args.binary, args.manifest)
        print("built artifact metadata: PASS")
    if args.schema_driver is not None:
        if args.manifest is None:
            parser.error("--schema-driver requires --manifest")
        verify_schema_driver(args.schema_driver, args.manifest)
        print("config schema metadata: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
