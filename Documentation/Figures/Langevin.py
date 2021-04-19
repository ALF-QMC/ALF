#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 05:20:44 2020

@author: fassaad
"""

import os
from py_alf import Simulation            # Interface with ALF
import numpy as np                       # Numerical library
sims = []      
for Time_step in [0.0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,  0.07,  0.08,  0.09,  0.1]:            # Vector of Simulation instances
    if Time_step == 0.0:
        sim_dict = {"Model": "Hubbard", 
                    "Lattice_type": "N_leg_ladder", 
                    "L1": 1 , "L2": 6, 
                    "Beta": 4.0, 
                    "Nsweep": 1000, 
                    "NBin": 50,
                    "Ltau": 0,
                    "Dtau": 0.1,
                    "Delta_t_Langevin_HMC" :  0.0,
                    }
    else:   
        sim_dict = {"Model": "Hubbard", 
                    "Lattice_type": "N_leg_ladder", 
                    "L1": 1 , "L2": 6, 
                    "Beta": 4.0, 
                    "Nsweep": 10000, 
                    "NBin": 400,
                    "Ltau": 0,
                    "Dtau": 0.1,
                    "Continuous" : True, 
                    "Langevin" : True,
                    "Delta_t_Langevin_HMC" :  Time_step, 
                    "Max_Force" : 100.0,
                    }    
    sim = Simulation('Hubbard', sim_dict,
                     #alf_dir= os.getenv('ALF_DIR', './ALF'),
                     #alf_dir='/Users/fassaad/Programs/ALF/Work',
                     branch = '165-introduce_langevin_updating_in_alf_2-0',
                     #machine= 'gnu',
                     alf_dir='/home/debian/Work',
                     machine= 'Intel',
                     mpi    = True,
                     n_mpi  = 12,
                     )
    #sim.sim_dir += 'Time_step=' + str(dict["Time_step"])
    sim.sim_dir += "_TL={}".format(sim_dict["Delta_t_Langevin_HMC"])        
    sims.append(sim)

sims[0].compile(target = "Hubbard")

En = np.empty((len(sims), 2))
TS = np.empty((len(sims), 1))
        
for i, sim in enumerate(sims):
    print (sim.sim_dir)
    sim.run()
    sim.analysis() 
    En[i] = sim.get_obs(['Ener_scalJ'])['Ener_scalJ']['obs'] 
    TS[i] = sim.sim_dict["Delta_t_Langevin_HMC"]  

with open('Langevin_ener.dat', 'w') as file:
    for i in range(len(sims) ):
        file.write('%6.6f\t' % (TS[i])) 
        file.write('%6.6f\t' % (En[i,0]))   
        file.write('%6.6f\t' % (En[i,1]))
        file.write('\n')
