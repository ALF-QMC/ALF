#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 05:20:44 2020

@author: fassaad
"""

from py_alf import Simulation            # Interface with ALF
import numpy as np                       # Numerical library
sims = []                                # Vector of Simulation instances
for Ham_Uf in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
        1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]:           # Values of Uf
    print(Ham_Uf)
    sim_dict = {"Model": "Kondo", 
                "Lattice_type": "Bilayer_square", 
                "L1": 4 , "L2": 4, 
                "Ham_Uf": Ham_Uf,  
                "Beta": 5.0, 
                "Nsweep": 500, 
                "NBin": 400,
                "Ltau": 0}
    sim = Simulation('Kondo', sim_dict,
                     alf_dir = '/home/debian/Work/',
                     #alf_dir = '/Users/fassaad/Programs/ALF/Work',
                     branch = '151-introduce-sun-kondo-hamiltonian-for-bilayer-lattices',
                     machine= 'Intel',
                     mpi    = True,
                     n_mpi  = 12)
    sims.append(sim)
#sims[0].compile(target = "Kondo")
Con = np.empty((len(sims), 2))         # Matrix for storing energy values
Uf  = np.empty((len(sims),))
Spin= np.empty((len(sims),16,2,2,4))
K   = np.empty((len(sims),16,2))           
for i, sim in enumerate(sims):
    print (sim.sim_dir)
    #sim.run()
    print (sim.sim_dir)
    sim.analysis() 
    Uf[i] = sim.sim_dict['Ham_Uf']                             # Store Uf value
    Con[i] = sim.get_obs(['Constraint_scalJ'])['Constraint_scalJ']['obs']  # Store constraint
    Spin[i]= sim.get_obs(['SpinZ_eqJK'])['SpinZ_eqJK']['dat']
    K[i]   = sim.get_obs(['SpinZ_eqJK'])['SpinZ_eqJK']['k']


with open('Constraint.dat', 'w') as file:
    for i in range(len(sims) ):
        file.write('%6.6f\t' % (Uf[i])) 
        file.write('%6.6f\t' % (Con[i,0])) 
        file.write('%6.6f\t' % (Con[i,1])) 
        file.write('%6.6f\t' % (Spin[i,15,1,1,0])) 
        file.write('%6.6f\t' % (Spin[i,15,1,1,1]))
        file.write('\n') 
