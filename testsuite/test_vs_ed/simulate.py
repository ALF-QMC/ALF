#!/usr/bin/env python3
import yaml

from py_alf import Simulation, ALF_source


if __name__ == "__main__":
    with open("spec.yaml", 'r', encoding='UTF-8') as f:
        spec = yaml.safe_load(f)

    sim = Simulation(
        ALF_source(),
        spec["ham_name"],
        spec["sim_dict"],
        sim_root='.',
        sim_dir='.',
        machine=spec['machine'],
        mpi=True,
        n_mpi=spec.get('n_mpi', 2),
    )
    sim.compile()
    sim.run()
