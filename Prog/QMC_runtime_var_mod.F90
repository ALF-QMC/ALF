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

#ifdef MPI
        Use mpi
#endif

        Implicit none

        Integer :: Nwrap, NSweep, NBin, NBin_eff,Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
        Real(Kind=Kind(0.d0)) :: CPU_MAX
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
        
        

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Initialization of the QMC runtime variables
!>
!
!--------------------------------------------------------------------
        contains
        subroutine set_QMC_runtime_default_var()

          implicit none   
          
          ! This is a set of variables that  identical for each simulation.
          Nwrap=0;  NSweep=0; NBin=0; Ltau=0; LOBS_EN = 0;  LOBS_ST = 0;  CPU_MAX = 0.d0
          Propose_S0 = .false. ;  Global_moves = .false. ; N_Global = 0
          Global_tau_moves = .false.; sequential = .true.; Langevin = .false. ; HMC =.false.
          Delta_t_Langevin_HMC = 0.d0;  Max_Force = 0.d0 ; Leapfrog_steps = 0; N_HMC_sweeps = 1
          Nt_sequential_start = 1 ;  Nt_sequential_end  = 0;  N_Global_tau  = 0;  Amplitude = 1.d0

        end subroutine set_QMC_runtime_default_var

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Bradcastinf of the QMC runtime variables
!>
!
!--------------------------------------------------------------------
#ifdef MPI

        subroutine broadcast_QMC_runtime_var(MPI_COMM_i)

          implicit none

          Integer :: ierr, MPI_COMM_i

          CALL MPI_BCAST(Nwrap                ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(NSweep               ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(NBin                 ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Ltau                 ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(LOBS_EN              ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(LOBS_ST              ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(CPU_MAX              ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Propose_S0           ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Global_moves         ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(N_Global             ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Global_tau_moves     ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Nt_sequential_start  ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Nt_sequential_end    ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(N_Global_tau         ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(sequential           ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Langevin             ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(HMC                  ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
          CALL MPI_BCAST(Leapfrog_steps       ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)

        end subroutine broadcast_QMC_runtime_var
#endif           
           
end Module QMC_runtime_var