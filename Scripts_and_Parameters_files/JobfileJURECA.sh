#!/bin/bash -x
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1 if the job is supposed to run on a single node.
# The usual policies apply, e.g. Nnodes*NtaskPnode = Ntasks.
# Most variables are self-explanatory, one exceptions might be Nthreads, which is referring to the number of OpenMP threads per MPI task.
# Useful configurations on JURECA (24 cores) for Nthreads are 1,2,4,6,12
# In general we found the ALF does not profit from hyper-threading such that we suggest to only use physical cores.
#
#SBATCH --nodes=##Nnodes##
#SBATCH --ntasks=##Ntasks##
#SBATCH --ntasks-per-node=##NtaskPnode##
#SBATCH --cpus-per-task=##Nthreads##
#SBATCH --hint=nomultithread
#SBATCH --time=##CPUMAX##:0:00
#SBATCH --partition=batch
#SBATCH --job-name=##NAME##
#SBATCH --error=error.txt
#SBATCH --output=out.txt
### start of jobscript

module load Intel
module load IntelMPI
module load imkl

export OMP_NUM_THREADS=##Nthreads##

# the following environment variables generate an optimal pinning (to the best of our knowledge)
# This DOES NOT have to be adapted to the choice of Ntasks
# FIRST EXCEPTION: If you chose to use hyper-threading (not recommended) you should set I_MPI_PIN_CELL=cpu
# SECOND EXCEPTION: The following environment variables are Intel specific.
export KMP_AFFINITY=verbose,granularity=fine,compact
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto:cache3
export I_MPI_PIN_ORDER=scatter

./out_to_in.sh >/dev/null 2>&1
srun ##EXECUTABLE##
