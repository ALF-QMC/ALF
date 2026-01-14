# Log of backward compatibility changes and critical bugs

## ALF 2.6
ALF 2.6 released on 2025-11-05

### 2025-09-29 Renaming Delta_S0_global to Get_Delta_S0 global
Author: J. Hofmann
Merge request !213

The function `Delta_S0_global` has been renamed to `Get_Delta_S0_global` as the new function now returns $\Delta S_0$ instead of $\exp(\Delta S_0)$. In case your model is using global moves, we advice to adapt your code accordingly. A fallback to `Delta_S0_global` is used for the time being such that backward compatibility is maintained

#### Optional changes
1) Modify and rename your implementation of `Delta_S0_global`

### 2024-07-10 Fix: Lattice in data.h5: Mixup of Norb and N_coord
Author: J. Schwab
Merge request !204

The lattice quantities `Norb` and `N_coord` were mistakenly mixed up in the HDF5 results file `data.h5`.
This has been fixed and the script `Analysis/fix-latt.py` was added to repair existing result files.

### 2024-02-01  Implement new  function  $(F,A) = \int d \omega F(\omega) A(\omega)$ in the stochastic maxent.
Author:  F. Assaad
Merge  request !196

1. Function F  is included in  parameter  list  of  the  stochastic maxent  routine


### 2024-01-25    Both  Classic MaxEnt  and  Stochastic  analytical continuation are  available <br>

Author:  J.Schwab and F. Assaad   <br>
Merge requests: !190 and !194

1. Both  approaches  allow  to  specify a  default model

2. To toggle  between the  Stochastic  and Classic approaches  a  new  logical variable,   Stochastic, has  been  introduced  in  VAR_Max_Stoch  name  space.  The  default  value  is Stochastic=.True. such that the  code  functions as  in the previous  releases.

3. There  is a  new  channel index  available P_PH for  single  particle  Green  functions  that are   particle-hole  symmetric. 


### 2023-08-14 Hubbard Stratonovich fields have been updated to be complex

Author:  F. Assaad
Merge request !176 
#### Breaking changes 

1. `Fields\_mod.F90`
    1.  `sigma%f`       is now complex    rather than real 
    2.  `nsigma%phi`    is now complex    rather than real 
2. `Hamilton\_main\_mod.F90` 
      1. `Hamiltonian\_set\_nsigma\_base(Initial\_field)`: `Initial\_field`  is  now complex, not  real 
      2. `Delta\_S0\_global\_base(Nsigma\_old)`: `Nsigma\_old`  is  complex  rather  than real 
      3. `S0\_base(n,nt,Hs\_new)` : `Hs\_new` is  now  complex rather  than real
      4. `Global\_move\_tau\_base` : The array `Flip\_value` is  now  a  complex one-dimensional array    



---
## ALF 2.5
ALF 2.5 released on 2023-06-05

### 2023-05-09 Use a RUNNING file to avoid multiple instances running in the same directory 

Author : J. Hofmann <br>
Merge request !155

### 2023-03-13  Improved  support  for  automatic HDF5 cmplilation 

Author : J. Schwab <br>
Merge request !151



---
## ALF 2.4
ALF 2.4 released on 2022-11-29



---
## ALF 2.3
ALF 2.3 released on 2022-06-24



### 2022-06-13 Work-around for (likely) preprocessor bug

Author : J.Schwab <br>
Merge request !139

### 2022-03-21 Reset fields when no update is proposed in Global_Updates

Author : A. Goetz <br>
Merge request !136

### 2022-01-31 Write parameters to HDF5 file

Author : J.Schwab <br>
Merge request !117

#### Breaking changes
1) Parameters to be formulated in format for parsing as described in Sec. 5.6 of documentation.
   Strictly speaking, it's not necessary to do that, but it simplifies the Hamiltonian,
   since the subroutine for reading parameters and writing parameters to HDF5 will be written automatically.
2) With HDF5: Add typebound procedure `write_parameters_hdf5` to Hamiltonian.

### 2021-12-08 Solves projector code runtime error

Author :  F. Parisen Toldin <br>
Merge request !129


---
## ALF 2.2
ALF 2.2 released on 2021-12-07


### 2021-11-21 Implementing HDF5

Author : J.Schwab <br>
Merge request !120

#### Breaking changes
1) In script configure.sh: The argument DEVEL/DEVELOPMENT is no longer a MACHINE name, but an optional switch

#### Optional changes
1) Added option for compiling with HDF5 by handing argument HDF5 to configure.sh


### 2021-11-21  Automatic computation of Hopping_Matrix_Type%Multiplicity

Author : F. Parisen Toldin <br>
Merge request !116

#### Breaking changes
1) Hopping_Matrix_Type%Multiplicity is now a private member, automatically initialized


### 2021-11-21  Test the checkerboard decomposition

Author : F. Parisen Toldin <br>
Merge request !124



---
## ALF 2.1
ALF 2.1 released on 2021-06-03



### 2021-03-22  Implementing Submodule Hamiltonians / All hamiltonians in one binary

Author : J. Schwab <br>
Merge request !107

#### Breaking changes
1) **In Hamiltonians** You will have to adapt your Hamiltonians to the Submodule structure
2) You will have to add your Hamiltonian name to the **Hamiltonians.list** in the Prog directory



---
## ALF 2.0
ALF 2.0 released on 2020-12-22



### 2020-11-16   Implementing  Langevin 

Author : F.F. Assaad <br>
Merge request !91 

#### Breaking changes
1) **In Hamiltonians** 

a) Mc\_step\_weight  parameter in ObserT and Obser routines <br>
b) Add 
`Subroutine Ham_Langevin_HMC_S0(Forces_0)`  <br>
Returns Bosonic forces

#### Optional changes
1) **Parameters    VAR_Hubbard**

Continuous = .F.  ! Uses (T: continuous; F: discrete) HS transformation

2) **Parameters  VAR_QMC**

a) Langevin = .F.    ! Langevin update <br>
b) Delta\_t\_Langevin\_HMC = 0.01 ! Default time step for Langevin and HMC updates <br>
c) Max\_Force            = 1.5  ! Max Force for  Langevin <br>
d) HMC     = .F.   ! HMC update <br>
e) Leapfrog_steps = 0 !  Number of leapfrog steps


### 2020-09-25   Embedding lattice information in observables 

Author :  J. Schwab <br>
Merge request !66 

#### Breaking changes
**In Hamiltonians** 

Calls to `Obser_Latt_make` should be adjusted to the subroutine's new interface
