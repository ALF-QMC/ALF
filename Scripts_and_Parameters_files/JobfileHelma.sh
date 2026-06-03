#!/bin/bash -l
#
# Sample jobscript for Helma CPU of NHR@FAU
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1 if the job is supposed to run on a single node.
# Most variables are self-explanatory, one exception might be Nthreads, which is referring to the number of OpenMP threads per MPI task.
# In general we found that ALF does not profit from hyper-threading such that we suggest to only use physical cores.
#
# On Helma, a single node has 384 cores.
# There are two partitions, "cpu" and "preempt_cpu".
# -"cpu": Single-node and multi-node jobs are possible. 
#         Multi-node jobs always get all 384 cores of each node,
#         while single-node jobs can request cores in multiples of 48 (one NUMA domain).
# -"preempt_cpu":
#     Only single node, but an arbitrary number of cores is allowed.
#     The walltime is up to 48 hours, but jobs might get cancelled after 2 hours
#     to make room for jobs in the "cpu" partition (with a grace time as in the GPU preempt partition).
#
# For single node jobs, the number of cores is ##NtaskPnode## * ##Nthreads##.
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

#available partitions: cpu, preempt_cpu
#SBATCH --partition=##PARTITION##
#SBATCH --nodes=##Nnodes##
#SBATCH --ntasks-per-node=##NtaskPnode##
#SBATCH --cpus-per-task=##Nthreads##

unset SLURM_EXPORT_ENV
module --force switch gpu-env/2025 cpu-env/2026
module load intel/2025.3.1
module load intelmpi/2021.17.0
module load mkl/2024.2.2

# the following environment variables generate an optimal pinning (to the best of our knowledge)
# This DOES NOT have to be adapted to the choice of Ntasks
# FIRST EXCEPTION: If you chose to use hyper-threading (not recommended) you should set I_MPI_PIN_CELL=cpu
# SECOND EXCEPTION: The following environment variables are Intel specific.
#export KMP_AFFINITY=verbose,granularity=fine,compact
export KMP_AFFINITY=granularity=fine,compact
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto:cache
export I_MPI_PIN_ORDER=scatter

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

 # Uncomment the following line to speed up MKL by making it believe it runs on an Intel CPU.
 # After `make`, `libfakeintel.so` is expected under `Libraries/Modules/`; use its absolute path here, e.g.:
 # export LD_PRELOAD=/absolute/path/to/ALF/Libraries/Modules/libfakeintel.so

bash ./out_to_in.sh
mpiexec -n $SLURM_NTASKS ##EXECUTABLE##
