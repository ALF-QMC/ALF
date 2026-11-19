#!/usr/bin/env python3
"""Check that baseline files differ from master only by the hooks listed in HOOKS.md.

Part of AI_development/ (tooling for master_ai, not part of the ALF code path).
Written with AI assistance; reviewed by: <name, date>.

For every file that exists in the base revision (default: master), the diff
between the base and the working tree must consist of added blocks only, and
each added block must be one of the blocks listed for that file in HOOKS.md.
Every listed block must occur exactly once. A Fortran hook (*.F90) must have
the form '#ifdef ALF_<FEAT>' / one line / '#endif'.

Usage (from the repository root of a master_ai checkout):
    AI_development/tools/hook_lint.py [--base master] [--hooks HOOKS.md]

Exit status 0 if the diff matches HOOKS.md, 1 otherwise.
"""
import argparse
import re
import subprocess
import sys
from collections import Counter

HEADING = re.compile(r"^### (\S+)\s*$")
HUNK = re.compile(r"^@@ -\d+(?:,\d+)? \+\d+(?:,\d+)? @@")


def listed_hooks(path):
    """Return Counter of (file, block) from HOOKS.md."""
    hooks = Counter()
    current = None
    block = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            m = HEADING.match(line)
            if m and block is None:
                current = m.group(1)
                continue
            if line.startswith("```"):
                if block is None:
                    block = []
                else:
                    if current is not None:
                        hooks[(current, "\n".join(block))] += 1
                    block = None
                    current = None
                continue
            if block is not None:
                block.append(line)
    return hooks


def diff_blocks(base):
    """Return (Counter of (file, block) added, list of problems)."""
    out = subprocess.run(["git", "diff", "-U0", "--no-color", "--no-ext-diff", base, "--"],
                         check=True, capture_output=True, text=True).stdout
    blocks = Counter()
    problems = []
    current_file = None
    new_file = False
    added = None
    for line in out.split("\n") + ["diff --git end"]:
        if line.startswith("diff --git"):
            if added and current_file:
                blocks[(current_file, "\n".join(added))] += 1
            current_file = None
            new_file = False
            added = None
            continue
        if added is None:
            # file header lines
            if line == "--- /dev/null":
                new_file = True  # added file, not a baseline file
            elif line.startswith("+++ b/") and not new_file:
                current_file = line[6:]
            elif line.startswith("+++ /dev/null"):
                problems.append("baseline file deleted (see the '--- a/' line above in git diff)")
            if not line.startswith("@@"):
                continue
        if line.startswith("@@"):
            if added and current_file:
                blocks[(current_file, "\n".join(added))] += 1
            added = []
            continue
        if added is None:
            continue
        if line.startswith("-") and current_file:
            problems.append("%s: baseline line removed or changed: %s" % (current_file, line[1:].strip()))
        elif line.startswith("+") and current_file:
            added.append(line[1:])
    return blocks, problems


def check_fortran_form(file, block):
    lines = block.split("\n")
    if not file.endswith(".F90"):
        return None
    if len(lines) != 3 or not re.match(r"^#ifdef ALF_[A-Z0-9_]+$", lines[0]) or lines[2] != "#endif":
        return "%s: hook is not of the form '#ifdef ALF_<FEAT>' / one line / '#endif':\n    %s" % (
            file, block.replace("\n", "\n    "))
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="master")
    ap.add_argument("--hooks", default="HOOKS.md")
    args = ap.parse_args()

    listed = listed_hooks(args.hooks)
    found, problems = diff_blocks(args.base)

    # New files appear in the diff with '--- /dev/null'; drop them (current_file None).
    for key, n in found.items():
        file, block = key
        form = check_fortran_form(file, block)
        if form:
            problems.append(form)
        if listed[key] != n:
            problems.append("%s: added block not listed in %s (found %d, listed %d):\n    %s" % (
                file, args.hooks, n, listed[key], block.replace("\n", "\n    ")))
    for key, n in listed.items():
        if key not in found:
            problems.append("%s: listed hook not found in the diff:\n    %s" % (
                key[0], key[1].replace("\n", "\n    ")))

    files = sorted({f for f, _ in found})
    print("hook_lint: %d hook blocks in %d baseline files (%s)" % (sum(found.values()), len(files), ", ".join(files)))
    for p in problems:
        print("ERROR " + p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
