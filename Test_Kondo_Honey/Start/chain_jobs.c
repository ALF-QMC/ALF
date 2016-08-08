#!/bin/bash -x
# submit a chain of jobs with dependency
# number of jobs to submit
jobname=Run_name
NO_OF_JOBS=Nj_r
# define jobscript
JOB_SCRIPT=job.q
I=0
export jobname1=$jobname$I
echo "sbatch --job-name=$jobname1 ${JOB_SCRIPT}"
JOBID=$(sbatch --job-name=$jobname1 ${JOB_SCRIPT} 2>&1 | awk '{print $(NF)}')
echo $JOBID

I=1
while [ ${I} -le ${NO_OF_JOBS} ]; do
export jobname1=$jobname$I
echo "sbatch --job-name=$jobname1 --dependency=afterok:${JOBID} ${JOB_SCRIPT}"
JOBID=$(sbatch --job-name=$jobname1 --dependency=afterok:${JOBID} ${JOB_SCRIPT} 2>&1 | awk '{print $(NF)}')
let I=${I}+1
done
