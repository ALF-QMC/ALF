#!/bin/bash
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1 if the job is supposed to run on a single node.
# Most variables are self-explanatory, one exceptions might be Nthreads, which is referring to the number of OpenMP threads per MPI task.
# Useful configurations on SuperMUC-NG (48 cores) for Nthreads are 1,2,4,6,12,24 and NtaskPnode = 48/Nthreads
# In general we found that ALF does not profit from hyper-threading such that we suggest to only use physical cores.
#
#SBATCH --job-name ##NAME##
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=##EMAIL##
# Wall clock limit (HH:MM:SS):
#SBATCH --time=##TIME##
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=##projectID##

#available partitions: test, micro, general, large, fat
#SBATCH --partition=##PARTITION##
#SBATCH --nodes=##Nnodes##
#SBATCH --ntasks-per-node=##NtaskPnode##
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=##Nthreads##

#Switch off Energy Aware Runtime for profiling or benchmarking
# #SBATCH --ear=off

#Important
module load slurm_setup

module switch mpi.intel  mpi.intel/2018
module switch intel intel/18.0
module switch mkl mkl/2018

# the follwing environment variables generate an optimal pinning (to the best of our knowledge)
# This DOES NOT have to be adapted to the choice of Ntasks
# FIRST EXCEPTION: If you chose to use hyper-threading (not recommended) you should set I_MPI_PIN_CELL=cpu
# SECOND EXCEPTION: The following environment variables are Intel specific.
#export KMP_AFFINITY=verbose,granularity=fine,compact
export KMP_AFFINITY=granularity=fine,compact
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto:cache
export I_MPI_PIN_ORDER=scatter

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


bash ./out_to_in.sh
mpiexec -n $SLURM_NTASKS ##EXECUTABLE##

