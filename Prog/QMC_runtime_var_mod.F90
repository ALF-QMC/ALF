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
        
        Use runtime_error_mod

#ifdef MPI
        Use mpi
#endif

        Implicit none

        Integer :: Nwrap, NSweep, NBin, NBin_eff,Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST
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


        subroutine set_default_values_measuring_interval(LOBS_ST, LOBS_EN, Thtrot, Ltrot, Projector)

            implicit none

            Integer :: LOBS_ST, LOBS_EN
            Integer, intent(in) :: Thtrot, Ltrot 
            Logical, intent(in) :: Projector

            if (Projector)  then
                if ( LOBS_ST == 0  ) then
                    LOBS_ST = Thtrot+1
                else
                    If (LOBS_ST < Thtrot+1 ) then
                        Write(error_unit,*) 'Measuring out of dedicating interval, LOBS_ST too small.'
                        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
                    endif
                endif
                if ( LOBS_EN == 0) then
                    LOBS_EN = Ltrot-Thtrot
                else
                    If (LOBS_EN > Ltrot-Thtrot ) then
                        Write(error_unit,*) 'Measuring out of dedicating interval, LOBS_EN too big.'
                        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
                    endif
                endif
            else
                if ( LOBS_ST == 0  ) then
                    LOBS_ST = 1
                endif
                if ( LOBS_EN == 0) then
                    LOBS_EN =  Ltrot
                endif
            endif

        end subroutine set_default_values_measuring_interval

        
        subroutine check_langevin_schemes_and_variables()
            if (sequential) then 
                write(output_unit,*) "Langevin mode does not allow sequential updates."
                write(output_unit,*) "Overriding Sequential=.True. from parameter files."
            endif
            if (HMC) then 
                write(output_unit,*) "Langevin mode does not allow HMC updates."
                write(output_unit,*) "Overriding HMC=.True. from parameter files."
            endif
            if (Global_moves) then 
                write(output_unit,*) "Langevin mode does not allow global updates."
                write(output_unit,*) "Overriding Global_moves=.True. from parameter files."
            endif
            if (Global_tau_moves) then 
                write(output_unit,*) "Langevin mode does not allow global tau updates."
                write(output_unit,*) "Overriding Global_tau_moves=.True. from parameter files."
            endif
#if defined(TEMPERING)
                    if ( N_exchange_steps > 0 ) then
                        write(output_unit,*) "Langevin mode does not allow tempering updates."
                        write(output_unit,*) "Overwriting N_exchange_steps to 0."
                    end if
#endif
        end subroutine check_langevin_schemes_and_variables
        
       
        ! Raise warnings for update schemes
        subroutine check_update_schemes_compatibility()

            implicit none 

            if ( .not. Sequential .and. Global_tau_moves) then
                write(output_unit,*) "Warning: Sequential = .False. and Global_tau_moves = .True."
                write(output_unit,*) "in the parameter file. Global tau updates will not occur if"
                write(output_unit,*) "Sequential is set to .False. ."
            endif

            if ( .not. Sequential .and. .not. HMC .and. .not. Langevin .and. .not. Global_moves) then
                write(output_unit,*) "Warning: no updates will occur as Sequential, HMC, Langevin, and"
                write(output_unit,*) "Global_moves are all .False. in the parameter file."
            endif

            if ( Sequential .and. Nt_sequential_end < Nt_sequential_start ) then
                write(output_unit,*) "Warning: Nt_sequential_end is smaller than Nt_sequential_start"
            endif
            
        end subroutine 
     

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Bradcast of the QMC runtime variables
!>
!
!--------------------------------------------------------------------
#ifdef MPI

        subroutine broadcast_QMC_runtime_var(MPI_COMM_i)

            implicit none

            Integer, intent(in) :: MPI_COMM_i

            Integer :: ierr

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
            CALL MPI_BCAST(N_HMC_sweeps         ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
            CALL MPI_BCAST(Max_Force            ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)
            CALL MPI_BCAST(Delta_t_Langevin_HMC ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)
            CALL MPI_BCAST(Amplitude            ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)

        end subroutine broadcast_QMC_runtime_var
#endif           
           

        
end Module QMC_runtime_var