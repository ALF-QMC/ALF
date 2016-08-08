#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=06:00:00
cd /homec/hwb03/hwb034/General_QMCT_svn/Kondo_Honey_new/Dir_r
pwd 
touch RUNNING
srun   /homec/hwb03/hwb034/General_QMCT_svn/Prog_7/Kondo_Honey.out
bash out_to_in.c
rm RUNNING
