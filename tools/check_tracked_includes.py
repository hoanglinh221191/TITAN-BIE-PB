#!/usr/bin/env python3
"""Fail if a tracked source includes a header that git does not track.

Three times in September 2026 a new engine source sat untracked while
tracked files included it, so the working tree built and a clean clone did
not.  CI catches that, but only after a push.  This runs in well under a
second and catches it at commit time.

    python tools/check_tracked_includes.py [--roots src ...]

Exit status is 1 when a tracked file includes an untracked file that exists
on disk.  Includes that resolve nowhere are ignored: those are system or
third-party headers found through -I paths, not this check's business.
"""
import argparse
import os
import re
import subprocess
import sys

INCLUDE = re.compile(r'^\s*#\s*include\s*"([^"]+)"', re.MULTILINE)
SOURCE_SUFFIXES = (".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx")


def tracked_files(repo):
    out = subprocess.run(["git", "-C", repo, "ls-files", "-z"],
                         capture_output=True, check=True).stdout
    return {p.replace("\\", "/") for p in out.decode("utf-8").split("\0") if p}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=".")
    parser.add_argument("--roots", nargs="*", default=["src"])
    args = parser.parse_args()
    repo = os.path.abspath(args.repo)

    tracked = tracked_files(repo)
    roots = tuple(r.rstrip("/") + "/" for r in args.roots)
    sources = [p for p in tracked
               if p.startswith(roots) and p.endswith(SOURCE_SUFFIXES)]

    problems = []
    for rel in sorted(sources):
        path = os.path.join(repo, rel)
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                text = handle.read()
        except OSError:
            continue
        here = os.path.dirname(rel)
        for target in INCLUDE.findall(text):
            # Resolve the way the compiler does: beside the including file
            # first, then from each root on the include path.
            for base in (here, *(r.rstrip("/") for r in roots)):
                candidate = os.path.normpath(os.path.join(base, target))
                candidate = candidate.replace("\\", "/")
                if candidate in tracked:
                    break
                if os.path.exists(os.path.join(repo, candidate)):
                    problems.append((rel, target, candidate))
                    break

    if not problems:
        print("check_tracked_includes: %d tracked sources, all includes tracked"
              % len(sources))
        return 0

    print("check_tracked_includes: %d include(s) resolve to untracked files.\n"
          "A clean clone of this commit will not compile.  Run git add on:\n"
          % len(problems), file=sys.stderr)
    for rel, target, candidate in problems:
        print("  %s includes %r -> %s (untracked)" % (rel, target, candidate),
              file=sys.stderr)
    print("\n  git add " + " ".join(sorted({c for _, _, c in problems})),
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
