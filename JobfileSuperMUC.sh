#!/bin/bash
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1 it the job is supposed to run on a single node.
# The usual policies apply, e.g. Nnodes*NtaskPnode = Ntasks.
# Most variables are selfexplaning, one exceptions might be Nthreads, which is refering to the number of OpenMP threads per MPI task.
# Useful configurations on Phase2 (28 cores) for Nthreads are 1,2,4,7,14 (4 is less suitable as one of the 7 tasked is ditributed across two sockets)
# In general we found the ALF does not profit from hyperthreading such that we suggest to only use physical cores.
#
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = micro
#@ node = ##Nnodes##
##@ total_tasks=##Ntasks##
##@ other version
#@ tasks_per_node = ##NtaskPnode##

#@ wall_clock_limit = ##CPUMAX##:00:00
#@ job_name = ##NAME##
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = ##workdir##
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=##EMAIL##
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh

#setup of environment
module switch mpi.ibm mpi.intel
module switch intel intel/17.0
module switch mkl mkl/2017

export OMP_NUM_THREADS=##Nthreads##

# the follwing eviroment variables generate an optimal pinning (to the best of our knowledge)
# This DOES NOT have to be addepted to the choice of Ntasks
# FIRST EXCEPTION: If you chose to use hyperthreading (not recommended) you should set I_MPI_PIN_CELL=cpu
# SECOND EXCEPTION: The following enviroment variables are Intel specific.
export KMP_AFFINITY=verbose,granularity=fine,compact
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto:cache3
export I_MPI_PIN_ORDER=scatter

./out_to_in.sh >/dev/null 2>&1
mpiexec -n ##Ntasks## ##EXECUTABLE##

