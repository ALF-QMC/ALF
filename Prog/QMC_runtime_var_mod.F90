!  Copyright (C) 2016 - 2022 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

Module QMC_runtime_var

        Use UDV_State_mod

        Implicit none

        Integer :: Nwrap, NSweep, NBin, NBin_eff,Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
        Integer :: NTAU, NTAU1
        Real(Kind=Kind(0.d0)) :: CPU_MAX
        Character (len=64) :: file_seeds, file_para, file_dat, file_info, ham_name
        Integer :: Seed_in
        Complex (Kind=Kind(0.d0)) , allocatable, dimension(:,:) :: Initial_field

        ! Space for choosing sampling scheme
        Logical :: Propose_S0, Tempering_calc_det
        Logical :: Global_moves, Global_tau_moves
        Integer :: N_Global
        Integer :: Nt_sequential_start, Nt_sequential_end, mpi_per_parameter_set
        Integer :: N_Global_tau
        Logical :: Sequential
        real (Kind=Kind(0.d0)) ::  Amplitude  !    Needed for  update of  type  3  and  4  fields.

!  Space for reading in Langevin & HMC  parameters
        Logical                      :: Langevin,  HMC
        Integer                      :: Leapfrog_Steps, N_HMC_sweeps
        Real  (Kind=Kind(0.d0))      :: Delta_t_Langevin_HMC, Max_Force


        NAMELIST /VAR_QMC/   Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX, &
              &               Propose_S0,Global_moves,  N_Global, Global_tau_moves, &
              &               Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
              &               sequential, Langevin, HMC, Delta_t_Langevin_HMC, &
              &               Max_Force, Leapfrog_steps, N_HMC_sweeps, Amplitude

        NAMELIST /VAR_HAM_NAME/ ham_name

        !  General
        Integer :: Ierr, I,nf, nf_eff, nst, n, n1, N_op
        Logical :: Toggle,  Toggle1
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Phase, Z, Z1
        Real    (Kind=Kind(0.d0)) :: ZERO = 10D-8, X, X1
        Real    (Kind=Kind(0.d0)) :: Mc_step_weight

        ! Storage for  stabilization steps
        Integer, dimension(:), allocatable :: Stab_nt 

        ! Space for storage.
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE :: udvst

        ! For tests
        Real (Kind=Kind(0.d0)) :: Weight, Weight_tot

        ! For the truncation of the program:
        logical                   :: prog_truncation, run_file_exists
        integer (kind=kind(0.d0)) :: count_bin_start, count_bin_end
        
        ! For MPI shared memory
        character(64), parameter :: name="ALF_SHM_CHUNK_SIZE_GB"
        character(64) :: chunk_size_str
        Real    (Kind=Kind(0.d0)) :: chunk_size_gb

end Module QMC_runtime_var