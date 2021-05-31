
## KNOWN ISSUES ##

1.  For the Kondo hamiltonian  the  conduction electron green function measures the  full green function matrix.  That is, the cc,  ff,  cf  and fc  green functions.   The user should hence make sure to adapt  the analysis so as  to produce the desired  result. 

2. We have detected a bug in both 2017 versions of Intel's MPI implementation (mpi.intel/2017(default) 
and mpi.intel/2017.2) if used in combination with the parallel (threaded) MKL library. The advise is to 
use either the 2016 suite (intel/16.0 (compiler), mpi.intel/5.1 and mkl/11.3 ) or the new 2018 suite 
(intel/18.0 (compiler), mpi.intel/2018 and mkl/2018). We did not detect this issue in both enviroments. 
You should also be aware the by default, dynamic linking is used. Hence if you use the 2016 or 2018 modules 
at compilations, the bug can reenter if you still load the 2017 versions at runtime. So please adapt your
configureHPC.sh as well as your Jobfiles for the loadleveler accordingly.
Additional note: In the serial version, the bug also seems to be absent. 
If you want to use the 2017 suite, you have to use the serial version of MKL (mkl/2017_s), which means you 
cannot profit from openMP multi-threading. This library is linked statically, hence taking care of this at 
compile time is sufficient and there is no need to adapt the Jobfiles.
WARNING: Even if you do not use parallel tempering actively, we still strongly suggest to take care of 
the above bug as it is extremely hard to estimate hidden influences and correlations of this memory 
allocation bug in the rest of the program. It is possible the other parts of the algorithm might be 
effected apart from the tempering exchange step even so we have absolutely no hint of additionally 
effected sections in ALF.

3. Intel suite 2017: Apparently, there seems to be a bug in the Intel MPI threaded memory allocator if both Intel MPI 
implementation and the PARALLEL version of MKL is used. This bug corrupts the data transmission during tempering moves such 
that the Monte Carlo is broken. The current advise is to either use an earlier version (2016 and before) or a newer one
(2018 and later). It should also be OK to use the sequential MKL library if openMP multi-threading is not required. 
Setting the environment variable OMP_NUM_THREADS=1 however does not fix the bug as this still uses the threaded MKL 
library (Date: 1. June 2018)


    

