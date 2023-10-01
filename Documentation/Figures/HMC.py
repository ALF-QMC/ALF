#!/usr/bin/env python3

import numpy as np                         # Numerical library
from py_alf import ALF_source, Simulation  # Interface with ALF

sim_dict = {"Model": "Hubbard", 
            "Lattice_type": "N_leg_ladder", 
            "L1": 1 , "L2": 6, 
            "Beta"      : 4.0,
            "Projector" : False, 
            "Ham_U"     : 4.0,  
            "ham_Tperp" : 1.0,  
            "Nsweep"    : 1000, 
            "NBin"      : 20,
            "Ltau"      : 0,
            "Mz"        : True,
            "HMC"       : True,
            "Sequential"       : False,
            "Delta_t_Langevin_HMC" : 0.1, 
            "Leapfrog_steps" : 50, 
            "Continuous" : True,
            }
            
sim = Simulation(ALF_source(alf_dir='/home/debian/ALF',branch='232-hybrid-monte-carlo-updates'), 'Hubbard', sim_dict,
                 machine= 'Intel',
                 mpi    = True,
                 n_mpi  = 12)

sim.compile()
sim.run()
sim.analysis() 

# Load all analysis results in a single Pandas dataframe
res = sim.get_obs()
print (res['Ener_scal0'], res['Ener_scal0_err'])

# Save all results in a single file
res.to_pickle('HMC.pkl')
