# ALF #
[![pipeline status](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/badges/master/pipeline.svg)](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/commits/master)
[![coverage report](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/badges/master/coverage.svg)](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/commits/master)
## General information ##
This version of the **A**lgorithms for **L**attice **F**ermions package provides a general code for the finite temperature auxiliary field Quantum Monte Carlo algorithm.       The code  is engineered to  be able simulate any model that can be written in terms of  sums of single body operators, of squares of single body operators and single body operators coupled to an Ising field with  given dynamics. We  provide predefined types that allow  the user to specify the model, the  Bravais lattice  as well as equal time and time displaced observables.     The code supports an MPI implementation.   Examples such as the Hubbard model on the Honeycomb lattice  as well as the Hubbard model  on the square lattice coupled to a transverse Ising field are  provided and discussed in the [documentation](https://git.physik.uni-wuerzburg.de/ALF/ALF/-/jobs/artifacts/master/raw/Documentation/doc.pdf?job=create_doc).

The Hamiltonians we can consider reads:
![The Hamiltonian0](Images/Hamiltonian0.png "The Hamiltonian")
where
![The Hamiltonian1](Images/Hamiltonian1.png "Parts explanation")

Here Z denotes an Ising spin variable with predefined dynamics. If your model can be written in this form then it will be amenable to the ALF. 

## Doxygen ##

You can find here [Doxygen](https://pawn.physik.uni-wuerzburg.de/~assaad/Doxygen_Docu/ALF/html/index.html)  formated documentation. (Work in progress)  

## PREREQUISITES ##

Libraries: Lapack and Blas

Compiler: gfortran  or ifort 


## CONFIGURATION FOR COMPILATION ##
**setenv.sh**   sets the default set of envorinment variables.  Do not change  this since this default set of  environment variables is required for the tests to run adequaltely.

**configureHPC.sh**  It is recommended to use this script to set the environment variables. Type ./configureHPC.sh to  browse through a list of options.

Once you have run the configuration script, change directory to Libraries, and to Analysis  and run the Makefiles there. In the Prog directory then type make examples.   The other programs are being updated to comply with the new version of the code.  

## FILES AND DIRECTORIES ##

**Libraries**    Libraries. Once that the environment is set in the file set_env.sh  the Libraries can be compiled with the **make** command. 

**Prog**   Main program and subroutines.  

**Analysis** Analysis programs. 

**Start**   This directory contain the files required to start a run. In particular it contains the parameter file   that specifies the model the lattice and various   parameters for the Monte Carlo run and  error analysis. 

**Examples** This directory provides a set of short example runs.  

**Documentation**  We have included in the file  [doc.pdf](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/-/jobs/artifacts/master/raw/Documentation/doc.pdf?job=create_doc)   an extensive documentation. The development of the doucmentation will take place **only** in the documentation_new branch. In this branch you all have permission to push things so that we can keep on improving  things. 
 

## TESTING ##

We have about 30 tests that test various parts of the program in the folder testsuite.
As testing framework we employ CTest.
From the subfolder testsuite the tests can be run as follows
- mkdir tests
- cd tests
- cmake ..
- make
- make test


## LICENSE ##
The various works that make up the ALF project are placed under licenses that put
a strong emphasis on the attribution of the original authors and the sharing of the contained knowledge.
To that end we have placed the ALF source code under the GPL version 3 license (see license.GPL and license.additional)
and took the liberty as per GPLv3 section 7 to include additional terms that deal with the attribution
of the original authors(see license.additional).
The Documentation of the ALF project by the ALF contributors is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License (see Documentation/license.CCBYSA)
We mention that we link against parts of lapack which licensed under a BSD license(see license.BSD).

## KNOWN ISSUES ##

We have detected a bug in both 2017 versions of Intel's MPI implementation (mpi.intel/2017(default) 
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

Intel suite 2017: Apparently, there seems to be a bug in the Intel MPI threaded memory allocator if both Intel MPI 
implementation and the PARALLEL version of MKL is used. This bug corrupts the data transmission during tempering moves such 
that the Monte Carlo is broken. The current advise is to either use an earlier version (2016 and before) or a newer one
(2018 and later). It should also be OK to use the sequential MKL library if openMP multi-threading is not required. 
Setting the environment variable OMP_NUM_THREADS=1 however does not fix the bug as this still uses the threaded MKL 
library (Date: 1. June 2018)


    

