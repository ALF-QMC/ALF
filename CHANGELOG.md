# Log of backward compatibility changes and critical bugs


## 2021-03-22  Implementing Submodule Hamiltonians / All hamiltonians in one binary

Author : J. Schwab <br>
Merge request !107

### Breaking changes
1) **In Hamiltonians** You will have to adapt your Hamiltonians to the Submodule structure
2) You will have to add your Hamiltonian name to the **Hamiltonians.list** in the Prog directory


## 2020-11-16   Implementing  Langevin 

Author : F.F. Assaad <br>
Merge request !91 

### Breaking changes
1) **In Hamiltonians** 

a) Mc\_step\_weight  parameter in ObserT and Obser routines <br>
b) Add 
`Subroutine Ham_Langevin_HMC_S0(Forces_0)`  <br>
Returns Bosonic forces

### Optional changes
1) **Parameters    VAR_Hubbard**

Continuous = .F.  ! Uses (T: continuous; F: discrete) HS transformation

2) **Parameters  VAR_QMC**

a) Langevin = .F.    ! Langevin update <br>
b) Delta\_t\_Langevin\_HMC = 0.01 ! Default time step for Langevin and HMC updates <br>
c) Max\_Force            = 1.5  ! Max Force for  Langevin <br>
d) HMC     = .F.   ! HMC update <br>
e) Leapfrog_steps = 0 !  Number of leapfrog steps



## 2020-09-25   Embedding lattice information in observables 

Author :  J. Schwab <br>
Merge request !66 

### Breaking changes
**In Hamiltonians** 

Calls to `Obser_Latt_make` should be adjusted to the subroutine's new interface
