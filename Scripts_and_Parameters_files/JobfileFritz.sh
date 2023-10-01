#!/bin/bash -l
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1 if the job is supposed to run on a single node.
# Most variables are self-explanatory, one exceptions might be Nthreads, which is referring to the number of OpenMP threads per MPI task.
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

#available partitions: singlenode, multinode
#SBATCH --partition=##PARTITION##
#SBATCH --nodes=##Nnodes##
#SBATCH --ntasks-per-node=##NtaskPnode##
#SBATCH --cpus-per-task=##Nthreads##

unset SLURM_EXPORT_ENV
module load intel
module load intelmpi
module load mkl

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
srun ##EXECUTABLE##

