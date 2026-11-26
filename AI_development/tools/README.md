# Checks for `master_ai`

Both scripts need only Python 3 and git; no extra packages.

## `hook_lint.py`: baseline files contain only the listed hooks

Run it from the root of a `master_ai` checkout:

```bash
python3 AI_development/tools/hook_lint.py --base master --hooks HOOKS.md
```

It looks at every file that exists in `--base` and requires:
- the diff of the working tree against `--base` adds blocks only, and removes or changes no line;
- every added block is listed in `HOOKS.md`, under that file;
- every listed block is present, exactly once;
- a Fortran hook has the form `#ifdef ALF_<FEAT>` / one line / `#endif`.

New files are not checked. The script exits with 1 and prints each problem.

## `baseline_guard.py`: bit-identical output without the switch

1. Build two executables with the **same compiler and flags**:
   - one from `master` (the reference);
   - one from your branch without any AI switch, or with the switch built in and the feature switched off in `parameters`.
2. Run the guard:

```bash
python3 AI_development/tools/baseline_guard.py \
    --ref  <master>/Prog/ALF.out \
    --test <branch>/Prog/ALF.out \
    --start <master>/Scripts_and_Parameters_files/Start \
    --work /tmp/guard --max-bins 2
```

**What it does:**
- It runs every serial scenario of `testsuite/test_branch_parameters.json` with both executables, using the same parameters and the same `seeds`.
- It compares every output file byte by byte, ignoring `info` and `info_git`.

**What gets skipped:**
- Scenarios that need MPI.
- Scenarios that both builds refuse at setup.

**Notes:**
- A `Devel` build is not bit-identical to an optimised build, so compare equal flags only.
- Before you rely on the guard, show once that it fails on a planted difference, for example a run with the feature switched on.
