#!/usr/bin/env python3
"""Short, deterministic sampling regressions using the normal ALF executable."""

import argparse
import itertools
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def parameters(*, ltau=0, projector=False, hmc=False, mala_sweeps=1,
               sweeps=8, bins=2):
    def logical(value):
        return ".true." if value else ".false."

    return f"""&VAR_ham_name
ham_name = 'Hubbard'
/
&VAR_lattice
Lattice_type = 'N_leg_ladder', Model = 'Hubbard', L1 = 4, L2 = 1
/
&VAR_Model_Generic
Checkerboard = .false., Symm = .true., Dtau = 0.1, Beta = 1.0
Projector = {logical(projector)}, Theta = 2.0
/
&VAR_Hubbard
ham_T = 1.0, ham_Tperp = 0.0, ham_U = 4.0, Mz = .true., Continuous = .true.
/
&VAR_QMC
Nwrap = 5, NSweep = {sweeps}, NBin = {bins}, Ltau = {ltau}
Sequential = .false., Global_MALA_moves = .true.
N_Global_MALA_sweeps = {mala_sweeps}, Delta_t_MALA_global = 0.1
Max_Force_MALA_global = 1.0
HMC = {logical(hmc)}, N_HMC_sweeps = 1, Leapfrog_steps = 1
Delta_t_Langevin_HMC = 0.1
/
"""


def run(binary, directory, *, previous=None, **options):
    directory.mkdir(parents=True)
    (directory / "parameters").write_text(parameters(**options))
    (directory / "seeds").write_text("8752143\n")
    if previous is not None:
        shutil.copyfile(previous / "confout_0", directory / "confin_0")
    with (directory / "stdout").open("w") as log:
        result = subprocess.run([str(binary)], cwd=directory, stdout=log,
                                stderr=subprocess.STDOUT, timeout=30, check=False)
    if result.returncode:
        raise AssertionError(f"ALF failed in {directory}:\n"
                             + (directory / "stdout").read_text())
    return directory


def numbers(path):
    assert path.is_file() and path.stat().st_size, f"Missing/empty observable: {path}"
    values = [float(token.replace("D", "E").replace("d", "e"))
              for token in re.split(r"[\s(),]+", path.read_text().strip()) if token]
    assert all(math.isfinite(value) for value in values), f"Non-finite output: {path}"
    return values


def close(a, b, description):
    assert len(a) == len(b), f"Different lengths: {description}"
    error = max((abs(x-y) for x, y in zip(a, b)), default=0.0)
    assert error < 1e-10, f"{description}: max absolute difference {error}"


def compare_configs(a, b):
    lines_a = (a / "confout_0").read_text().splitlines()
    lines_b = (b / "confout_0").read_text().splitlines()
    assert lines_a[0].split() == lines_b[0].split(), "Different final RNG states"
    close([float(x) for x in lines_a[1:]], [float(x) for x in lines_b[1:]],
          "Final configurations")


def measurements(binary, root):
    for projector, hmc, moves in itertools.product([False, True], [False, True], [1, 3]):
        runs = []
        for ltau in [0, 1]:
            name = f"projector-{projector}_hmc-{hmc}_moves-{moves}_ltau-{ltau}"
            directory = run(binary, root / name, projector=projector, hmc=hmc,
                            mala_sweeps=moves, ltau=ltau)
            for observable in ["Ener_scal", "Kin_scal", "Pot_scal", "Part_scal", "Green_eq"]:
                numbers(directory / observable)
            if ltau:
                numbers(directory / "Green_tau")
            runs.append(directory)
        # Time-displaced measurement must not alter the Markov chain. For finite
        # temperature it also visits the same equal-time measurement slices;
        # comparison detects extra counting of the post-MALA state in mixed runs.
        compare_configs(*runs)
        if not projector:
            for observable in ["Ener_scal", "Kin_scal", "Pot_scal", "Part_scal"]:
                close(numbers(runs[0] / observable), numbers(runs[1] / observable), observable)


def restart(binary, root):
    for ltau, moves in itertools.product([0, 1], [1, 3]):
        directory = root / f"ltau-{ltau}_moves-{moves}"
        whole = run(binary, directory / "continuous", hmc=True, ltau=ltau,
                    mala_sweeps=moves, sweeps=10, bins=1)
        previous = None
        for step in range(10):
            previous = run(binary, directory / f"restart-{step}", previous=previous,
                           hmc=True, ltau=ltau, mala_sweeps=moves, sweeps=1, bins=1)
        compare_configs(whole, previous)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("scenario", choices=["measurements", "restart"])
    args = parser.parse_args()
    # Keep failed runs for inspection; successful test runs clean up their files.
    root = Path(tempfile.mkdtemp(prefix=f"alf-mala-{args.scenario}-"))
    print(f"Test runs: {root}", flush=True)
    {"measurements": measurements, "restart": restart}[args.scenario](args.binary.resolve(), root)
    shutil.rmtree(root)
    print(f"MALA {args.scenario} checks passed")
