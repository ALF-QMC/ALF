#!/usr/bin/env python3
import sys
import shutil

import yaml

from py_alf import Simulation, ALF_source


if __name__ == "__main__":
    machine = sys.argv[1]
    with open("spec.yaml", 'r', encoding='UTF-8') as f:
        spec = yaml.safe_load(f)

    alf_src = ALF_source(alf_dir='../../..')
    sims = []
    for i, sim_dict in enumerate(spec["sim_dicts"]):
        sim = Simulation(
            alf_src,
            spec["ham_name"],
            sim_dict,
            sim_root='.',
            sim_dir=f'{i+1}',
            machine=machine,
            mpi=True,
        )
        sims.append(sim)

    sims[0].compile()
    for sim in sims:
        sim.run(only_prep=True)
    shutil.copy('../../../Prog/ALF.out', '.')
