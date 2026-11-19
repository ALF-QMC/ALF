#!/usr/bin/env python3
"""Check that a build of master_ai without AI switches is bit-identical to master.

Part of AI_development/ (tooling for master_ai, not part of the ALF code path).
Written with AI assistance; reviewed by: <name, date>.

Runs every serial scenario of testsuite/test_branch_parameters.json with two
ALF executables (reference and test) in fresh directories, with the same
parameters and the same seeds file, and compares every output file byte by
byte. The files 'info' and 'info_git' are ignored, because they record build
and run metadata. testsuite/compare_dirs.py is not used here: it compares
analysed results within a tolerance (np.allclose), not exactly.

Usage:
    baseline_guard.py --ref <master>/Prog/ALF.out --test <master_ai>/Prog/ALF.out
                      --start <ALF>/Scripts_and_Parameters_files/Start
                      [--work DIR] [--only NAME ...] [--max-bins N]

Exit status 0 if every scenario matches, 1 otherwise.
"""
import argparse
import filecmp
import json
import os
import re
import shutil
import subprocess
import sys

IGNORED = {"info", "info_git"}


def set_parameter(text, key, value):
    """Set every assignment of `key` in a namelist file; return (text, count)."""
    if isinstance(value, bool):
        val = ".T." if value else ".F."
    elif isinstance(value, str):
        val = '"%s"' % value
    else:
        val = repr(value)
    pattern = re.compile(r"^(\s*)%s(\s*)=\s*[^!\n]*" % re.escape(key),
                         re.IGNORECASE | re.MULTILINE)
    return pattern.subn(lambda m: "%s%s%s= %s " % (m.group(1), key, m.group(2), val), text)


def make_parameters(template, pars, max_bins):
    text = template
    missing = []
    for key, value in pars.items():
        if key == "Model":
            key = "ham_name"
        if key == "NBin" and max_bins is not None:
            value = min(value, max_bins)
        text, n = set_parameter(text, key, value)
        if n == 0:
            missing.append(key)
    return text, missing


def run(exe, directory, parameters, seeds):
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    with open(os.path.join(directory, "parameters"), "w") as f:
        f.write(parameters)
    shutil.copy(seeds, os.path.join(directory, "seeds"))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with open(os.path.join(directory, "stdout"), "w") as out:
        res = subprocess.run([exe], cwd=directory, stdout=out,
                             stderr=subprocess.STDOUT, env=env)
    return res.returncode


def compare(dir_ref, dir_test):
    """Return a list of differences (empty if bit-identical)."""
    diffs = []
    files_ref = {f for f in os.listdir(dir_ref) if f not in IGNORED}
    files_test = {f for f in os.listdir(dir_test) if f not in IGNORED}
    for f in sorted(files_ref ^ files_test):
        diffs.append("only in one run: %s" % f)
    for f in sorted(files_ref & files_test):
        a, b = os.path.join(dir_ref, f), os.path.join(dir_test, f)
        if os.path.isfile(a) and not filecmp.cmp(a, b, shallow=False):
            diffs.append("differs: %s" % f)
    return diffs


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="reference ALF.out (master)")
    ap.add_argument("--test", required=True, help="ALF.out under test (master_ai, no AI switch)")
    ap.add_argument("--start", required=True, help="directory with template 'parameters' and 'seeds'")
    ap.add_argument("--scenarios", default=os.path.join(here, "..", "..", "testsuite", "test_branch_parameters.json"))
    ap.add_argument("--work", default="baseline_guard_runs")
    ap.add_argument("--only", nargs="*", help="run only these scenarios")
    ap.add_argument("--max-bins", type=int, default=None, help="cap NBin to save time")
    args = ap.parse_args()

    with open(os.path.join(args.start, "parameters")) as f:
        template = f.read()
    seeds = os.path.join(args.start, "seeds")
    with open(args.scenarios) as f:
        scenarios = json.load(f)

    failed, skipped, passed = [], [], []
    for name, (_, pars) in scenarios.items():
        if args.only and name not in args.only:
            continue
        if isinstance(pars, list):
            skipped.append((name, "needs MPI (parameter list)"))
            continue
        parameters, missing = make_parameters(template, pars, args.max_bins)
        if missing:
            skipped.append((name, "keys not in template: " + ", ".join(missing)))
            continue
        dirs = {}
        rcs = {}
        for tag, exe in (("ref", args.ref), ("test", args.test)):
            dirs[tag] = os.path.join(args.work, name, tag)
            rcs[tag] = run(os.path.abspath(exe), dirs[tag], parameters, seeds)
        if rcs["ref"] != 0 and rcs["test"] == rcs["ref"]:
            # The template does not suit this scenario (for example a missing
            # lattice setting); both builds refuse it in the same way.
            skipped.append((name, "both runs stop with exit code %d" % rcs["ref"]))
            continue
        if rcs["ref"] != 0 or rcs["test"] != 0:
            failed.append((name, ["exit codes ref=%d test=%d" % (rcs["ref"], rcs["test"])]))
            continue
        diffs = compare(dirs["ref"], dirs["test"])
        (failed if diffs else passed).append((name, diffs))

    for name, _ in passed:
        print("PASS  %s" % name)
    for name, why in skipped:
        print("SKIP  %s (%s)" % (name, why))
    for name, diffs in failed:
        print("FAIL  %s" % name)
        for d in diffs:
            print("        " + d)
    print("%d passed, %d failed, %d skipped" % (len(passed), len(failed), len(skipped)))
    return 1 if failed or not passed else 0


if __name__ == "__main__":
    sys.exit(main())
