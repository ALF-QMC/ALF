#!/usr/bin/env python3
"""NBin and CPU_MAX are both bounds; the run stops at whichever comes first.

Runs the Start parameter set with several (NBin, CPU_MAX) pairs and checks the
bins actually written to data.h5. Needs a binary built with HDF5, and runs it
directly, so noMPI:

    . configure.sh GNU noMPI HDF5 && make program
    ./testsuite/test_nbin_cpu_max.py

Exits non-zero on the first failed check.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py

ALF_DIR = Path(__file__).resolve().parent.parent
ALF_BIN = ALF_DIR / "Prog" / "ALF.out"
START_DIR = ALF_DIR / "Scripts_and_Parameters_files" / "Start"

SMALL = {
    "VAR_lattice": {"L1": 4, "L2": 4},
    "VAR_QMC": {"NSweep": 20, "Ltau": 0},
}


def rewrite_parameters(path, overrides):
    """Rewrite ``key = value`` lines of path, per namelist."""
    namelist = None
    out = []

    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("&"):
            namelist = stripped[1:].split()[0]
        elif stripped.startswith("/"):
            namelist = None
        elif "=" in line and namelist in overrides:
            key = line.split("=", 1)[0].strip()
            for want, value in overrides[namelist].items():
                if key.lower() == want.lower():
                    line = f"{key:<20} = {value}"
                    break
        out.append(line)
    path.write_text("\n".join(out) + "\n")


def run_alf(run_dir, nbin, cpu_max):
    """Run ALF once in run_dir and return the bins in data.h5."""
    if not run_dir.exists():
        shutil.copytree(START_DIR, run_dir)
    rewrite_parameters(
        run_dir / "parameters",
        {
            "VAR_lattice": dict(SMALL["VAR_lattice"]),
            "VAR_QMC": {**SMALL["VAR_QMC"], "NBin": nbin, "CPU_MAX": f"{cpu_max}d0"},
        },
    )
    subprocess.run([str(ALF_BIN)], cwd=run_dir, check=True, timeout=900)
    with h5py.File(run_dir / "data.h5", "r") as f:
        return int(f["Ener_scal/obser"].shape[0])  # just use energy for bin count


def check(name, ok, detail):
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")
    return ok


def main():
    # check binary
    if not ALF_BIN.exists():
        sys.exit(f"{ALF_BIN} not built; run `make program` first")

    results = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # A reachable NBin must be honoured exactly. The time budget here is far
        # larger than the work, so an implementation that lets CPU_MAX discard
        # NBin runs the full budget and overshoots.
        bins = run_alf(tmp / "bin_bound", nbin=20, cpu_max=0.05)
        results.append(check("bin bound honoured", bins == 20, f"{bins} bins, want 20"))

        # A tiny budget must truncate cleanly, well short of NBin. This is the
        # graceful stop a checkpoint restart depends on, for example.
        bins = run_alf(tmp / "time_bound", nbin=1_000_000, cpu_max=0.003)
        results.append(
            check("time bound truncates", 0 < bins < 1_000_000, f"{bins} bins")
        )

        # NBin <= 0 leaves the run bounded by time alone.
        bins = run_alf(tmp / "no_bin_bound", nbin=0, cpu_max=0.003)
        results.append(check("NBin <= 0 is time-bounded", bins > 0, f"{bins} bins"))

        # NBin bounds this run alone; data.h5 accumulates across restarts.
        restart = tmp / "restart"
        bins = run_alf(restart, nbin=10, cpu_max=0.05)
        results.append(check("restart, first leg", bins == 10, f"{bins} bins, want 10"))

        # confout -> confin rename
        promoted = list(restart.glob("confout_*"))
        if not promoted:
            sys.exit("no checkpoint to promote")
        for conf in promoted:
            conf.rename(restart / f"confin_{conf.name[len('confout_') :]}")
        bins = run_alf(restart, nbin=5, cpu_max=0.05)

        # Check that we get the right total after the restart run
        results.append(
            check("restart, second leg", bins == 15, f"{bins} bins, want 15")
        )

    if not all(results):
        sys.exit(1)
    print(f"{len(results)} tests passed")


if __name__ == "__main__":
    main()
