#!/bin/bash
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1 if the job is supposed to run on a single node.
# Most variables are selfexplaning, one exceptions might be Nthreads, which is refering to the number of OpenMP threads per MPI task.
# Useful configurations on SuperMUC-NG (48 cores) for Nthreads are 1,2,4,6,12,24 and NtaskPnode = 48/Nthreads
# In general we found that ALF does not profit from hyperthreading such that we suggest to only use physical cores.
#
#SBATCH --job-name  DIR_R
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assaad@physik.uni-wuerzburg.de
# Wall clock limit (HH:MM:SS):
#SBATCH --time=48:00:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pr53ju


#available partitions: test, micro, general, large, fat
#SBATCH --partition=micro
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2

#Switch off Energy Aware Runtime for profiling or benchmarking
# #SBATCH --ear=off

#Important
module load slurm_setup

module load mpi.intel
module load intel 
module load mkl 

cd /hppfs/work/pr53ju/lu57gek5/AF_F_Run/DIR_R
touch RUNNING

# the follwing eviroment variables generate an optimal pinning (to the best of our knowledge)
# This DOES NOT have to be addepted to the choice of Ntasks
# FIRST EXCEPTION: If you chose to use hyperthreading (not recommended) you should set I_MPI_PIN_CELL=cpu
# SECOND EXCEPTION: The following enviroment variables are Intel specific.
#export KMP_AFFINITY=verbose,granularity=fine,compact
export KMP_AFFINITY=granularity=fine,compact
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto:cache
export I_MPI_PIN_ORDER=scatter

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $SLURM_NTASKS $SLURM_CPUS_PER_TASK

bash ./out_to_in.sh
mpiexec -n $SLURM_NTASKS  /dss/dsshome1/05/lu57gek5/AF_F_Code/Prog/AF_F.out

rm RUNNING
