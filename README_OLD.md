# General finite temperature auxiliary field code - ALF#

## General information ##
Last version of the code is in Prog_8

Last version of the analysis files is in Analysis_8

1)    Implemented  optimizations made by Johannes.  The code is thereby not necessarily faster !  You  may want to  try the code with the  original routines upgrade_FFA.f90  Operator_FFA.f90. You will find them in the directory Prog_8/FFA_Originals

2)    In the file parameter you can read in a  variable  CPU_MAX in the namespace &VAR_QMC. If present the code will carry out as many bins as allowed in the given amount of time.   If not present the code will just react as usual and carry out the requested number of bins. (Implemented by Martin  Bercx). 

3)    In this version,  NWRAP does not have to be a multiple of Ltrot.  This greatly facilitates temperature scans. 

4)    I have included the Kondo model on the honeycomb lattice with z-frustrating interactions.

5)    I have included  MaxEnt wrappers in the code.  You can find some minimal notes in the documentation.   I have also included some scripts in the test-run directory of the Hubbard model.  Note that you will have to recompile  the libraries. 

## PREREQUISITES ##

Libraries: Lapack and Blas

Compiler: gfortran  or ifort 


## CONFIGURATION FOR COMPILATION ##
See set_env.sh

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
To that end we have placed the ALF source code under the GPL version 3 license (see license.GPL)
and took the liberty as per GPLv3 section 7 to include additional terms that deal with the attribution
of the original authors(see license.additional).
The Documentation of the ALF project by the ALF contributors is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.
We mention that we link against parts of lapack which licensed under a BSD license(see license.BSD).

## USAGE ##

### INPUT  ###

### OUTPUT ###

### Questions ###

How do we organize the documentation? 

## TODO DOC ##


Johannes Hofmann: documentation for the  SPT model. 

Toshihiro Sato, Martin Hohenadler:  documentation for the Ising model. 

Fakher Assaad:   documentation for the Hubbard model. 

Martin Bercx and Fakher Assaad: general description for the package. 

## TODO CODE ##
cmake

