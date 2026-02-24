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
        Use iso_fortran_env, only: output_unit, error_unit

#ifdef MPI
        Use mpi
#endif

        Implicit none

    ! Make everything private by default; expose explicit APIs below
    private

        Integer :: Nwrap, NSweep, NBin,Ltau, LOBS_EN, LOBS_ST
        Real(Kind=Kind(0.d0)) :: CPU_MAX
        ! Space for choosing sampling scheme
        Logical :: Propose_S0, Tempering_calc_det
        Logical :: Global_moves, Global_tau_moves
        Integer :: N_Global
        Integer :: Nt_sequential_start, Nt_sequential_end, mpi_per_parameter_set
        Integer :: N_Global_tau
        Logical :: Sequential
        real (Kind=Kind(0.d0)) ::  Amplitude  !    Needed for  update of  type  3  and  4  fields.

        Character (len=64) :: ham_name

        Logical :: settings_locked = .false.


        !  Space for reading in Langevin & HMC  parameters
        Logical                      :: Langevin,  HMC
        Integer                      :: Leapfrog_Steps, N_HMC_sweeps
        Real  (Kind=Kind(0.d0))      :: Delta_t_Langevin_HMC, Max_Force
 
          
#if defined(TEMPERING)
        Integer :: N_exchange_steps, N_Tempering_frequency
        NAMELIST /VAR_TEMP/  N_exchange_steps, N_Tempering_frequency, mpi_per_parameter_set, Tempering_calc_det
#endif

        NAMELIST /VAR_QMC/   Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX, &
              &               Propose_S0,Global_moves,  N_Global, Global_tau_moves, &
              &               Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
              &               sequential, Langevin, HMC, Delta_t_Langevin_HMC, &
              &               Max_Force, Leapfrog_steps, N_HMC_sweeps, Amplitude
       
        NAMELIST /VAR_HAM_NAME/ ham_name

    !=====================
    ! Public API exports
    !=====================
    public :: set_QMC_runtime_default_var
    public :: set_default_values_measuring_interval
    public :: check_langevin_schemes_and_variables
    public :: check_update_schemes_compatibility
#ifdef MPI
    public :: broadcast_QMC_runtime_var
#endif
#ifdef TEMPERING
    public :: read_and_broadcast_TEMPERING_var
#endif
    public :: read_and_broadcast_QMC_var_and_ham_name

    public :: lock_QMC_runtime_settings

    ! Accessors (getters/setters) for runtime variables
    public :: get_Nwrap,             set_Nwrap
    public :: get_NSweep,            set_NSweep
    public :: get_NBin,              set_NBin
    public :: get_Ltau,              set_Ltau
    public :: get_LOBS_EN,           set_LOBS_EN
    public :: get_LOBS_ST,           set_LOBS_ST
    public :: get_CPU_MAX,           set_CPU_MAX
    public :: get_Propose_S0,        set_Propose_S0
    public :: get_Global_moves,      set_Global_moves
    public :: get_N_Global,          set_N_Global
    public :: get_Global_tau_moves,  set_Global_tau_moves
    public :: get_Nt_sequential_start, set_Nt_sequential_start
    public :: get_Nt_sequential_end,   set_Nt_sequential_end
    public :: get_N_Global_tau,      set_N_Global_tau
    public :: get_sequential,        set_sequential
    public :: get_Langevin,          set_Langevin
    public :: get_HMC,               set_HMC
    public :: get_Leapfrog_steps,    set_Leapfrog_steps
    public :: get_N_HMC_sweeps,      set_N_HMC_sweeps
    public :: get_Max_Force,         set_Max_Force
    public :: get_Delta_t_Langevin_HMC, set_Delta_t_Langevin_HMC
    public :: get_Amplitude,         set_Amplitude
    public :: get_ham_name,          set_ham_name
    public :: get_mpi_per_parameter_set, set_mpi_per_parameter_set
    public :: get_Tempering_calc_det, set_Tempering_calc_det
#if defined(TEMPERING)
    public :: get_N_exchange_steps,  set_N_exchange_steps
    public :: get_N_Tempering_frequency, set_N_Tempering_frequency
#endif
        
        

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
            settings_locked = .false.

        end subroutine set_QMC_runtime_default_var


        subroutine set_default_values_measuring_interval(Thtrot, Ltrot, Projector)

            implicit none

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

            implicit none 
            
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

            call lock_QMC_runtime_settings()
            
        end subroutine 

        subroutine lock_QMC_runtime_settings()
            settings_locked = .true.
        end subroutine lock_QMC_runtime_settings

        subroutine ensure_settings_unlocked(caller)
            character(len=*), intent(in) :: caller

            if (settings_locked) then
                write(error_unit,*) trim(caller), ": settings are locked and cannot be modified."
                call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            endif
        end subroutine ensure_settings_unlocked

        !=========================
        ! Accessors implementation
        !=========================

        ! Integer getters/setters
        integer function get_Nwrap() result(val)
            val = Nwrap
        end function get_Nwrap

        subroutine set_Nwrap(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_Nwrap")
            Nwrap = val
        end subroutine set_Nwrap

        integer function get_NSweep() result(val)
            val = NSweep
        end function get_NSweep

        subroutine set_NSweep(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_NSweep")
            NSweep = val
        end subroutine set_NSweep

        integer function get_NBin() result(val)
            val = NBin
        end function get_NBin

        subroutine set_NBin(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_NBin")
            NBin = val
        end subroutine set_NBin

        integer function get_Ltau() result(val)
            val = Ltau
        end function get_Ltau

        subroutine set_Ltau(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_Ltau")
            Ltau = val
        end subroutine set_Ltau

        integer function get_LOBS_EN() result(val)
            val = LOBS_EN
        end function get_LOBS_EN

        subroutine set_LOBS_EN(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_LOBS_EN")
            LOBS_EN = val
        end subroutine set_LOBS_EN

        integer function get_LOBS_ST() result(val)
            val = LOBS_ST
        end function get_LOBS_ST

        subroutine set_LOBS_ST(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_LOBS_ST")
            LOBS_ST = val
        end subroutine set_LOBS_ST

        integer function get_N_Global() result(val)
            val = N_Global
        end function get_N_Global

        subroutine set_N_Global(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_N_Global")
            N_Global = val
        end subroutine set_N_Global

        integer function get_Nt_sequential_start() result(val)
            val = Nt_sequential_start
        end function get_Nt_sequential_start

        subroutine set_Nt_sequential_start(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_Nt_sequential_start")
            Nt_sequential_start = val
        end subroutine set_Nt_sequential_start

        integer function get_Nt_sequential_end() result(val)
            val = Nt_sequential_end
        end function get_Nt_sequential_end

        subroutine set_Nt_sequential_end(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_Nt_sequential_end")
            Nt_sequential_end = val
        end subroutine set_Nt_sequential_end

        integer function get_N_Global_tau() result(val)
            val = N_Global_tau
        end function get_N_Global_tau

        subroutine set_N_Global_tau(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_N_Global_tau")
            N_Global_tau = val
        end subroutine set_N_Global_tau

        integer function get_Leapfrog_steps() result(val)
            val = Leapfrog_Steps
        end function get_Leapfrog_steps

        subroutine set_Leapfrog_steps(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_Leapfrog_steps")
            Leapfrog_Steps = val
        end subroutine set_Leapfrog_steps

        integer function get_N_HMC_sweeps() result(val)
            val = N_HMC_sweeps
        end function get_N_HMC_sweeps

        subroutine set_N_HMC_sweeps(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_N_HMC_sweeps")
            N_HMC_sweeps = val
        end subroutine set_N_HMC_sweeps

        integer function get_mpi_per_parameter_set() result(val)
            val = mpi_per_parameter_set
        end function get_mpi_per_parameter_set

        subroutine set_mpi_per_parameter_set(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_mpi_per_parameter_set")
            mpi_per_parameter_set = val
        end subroutine set_mpi_per_parameter_set

#if defined(TEMPERING)
        integer function get_N_exchange_steps() result(val)
            val = N_exchange_steps
        end function get_N_exchange_steps

        subroutine set_N_exchange_steps(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_N_exchange_steps")
            N_exchange_steps = val
        end subroutine set_N_exchange_steps

        integer function get_N_Tempering_frequency() result(val)
            val = N_Tempering_frequency
        end function get_N_Tempering_frequency

        subroutine set_N_Tempering_frequency(val)
            integer, intent(in) :: val
            call ensure_settings_unlocked("set_N_Tempering_frequency")
            N_Tempering_frequency = val
        end subroutine set_N_Tempering_frequency
#endif

        ! Real getters/setters
        real(kind=kind(0.d0)) function get_CPU_MAX() result(val)
            val = CPU_MAX
        end function get_CPU_MAX

        subroutine set_CPU_MAX(val)
            real(kind=kind(0.d0)), intent(in) :: val
            call ensure_settings_unlocked("set_CPU_MAX")
            CPU_MAX = val
        end subroutine set_CPU_MAX

        real(kind=kind(0.d0)) function get_Max_Force() result(val)
            val = Max_Force
        end function get_Max_Force

        subroutine set_Max_Force(val)
            real(kind=kind(0.d0)), intent(in) :: val
            call ensure_settings_unlocked("set_Max_Force")
            Max_Force = val
        end subroutine set_Max_Force

        real(kind=kind(0.d0)) function get_Delta_t_Langevin_HMC() result(val)
            val = Delta_t_Langevin_HMC
        end function get_Delta_t_Langevin_HMC

        subroutine set_Delta_t_Langevin_HMC(val)
            real(kind=kind(0.d0)), intent(in) :: val
            call ensure_settings_unlocked("set_Delta_t_Langevin_HMC")
            Delta_t_Langevin_HMC = val
        end subroutine set_Delta_t_Langevin_HMC

        real(kind=kind(0.d0)) function get_Amplitude() result(val)
            val = Amplitude
        end function get_Amplitude

        subroutine set_Amplitude(val)
            real(kind=kind(0.d0)), intent(in) :: val
            call ensure_settings_unlocked("set_Amplitude")
            Amplitude = val
        end subroutine set_Amplitude

        ! Logical getters/setters
        logical function get_Propose_S0() result(val)
            val = Propose_S0
        end function get_Propose_S0

        subroutine set_Propose_S0(val)
            logical, intent(in) :: val
            call ensure_settings_unlocked("set_Propose_S0")
            Propose_S0 = val
        end subroutine set_Propose_S0

        logical function get_Global_moves() result(val)
            val = Global_moves
        end function get_Global_moves

        subroutine set_Global_moves(val)
            logical, intent(in) :: val
            call ensure_settings_unlocked("set_Global_moves")
            Global_moves = val
        end subroutine set_Global_moves

        logical function get_Global_tau_moves() result(val)
            val = Global_tau_moves
        end function get_Global_tau_moves

        subroutine set_Global_tau_moves(val)
            logical, intent(in) :: val
            call ensure_settings_unlocked("set_Global_tau_moves")
            Global_tau_moves = val
        end subroutine set_Global_tau_moves

        logical function get_sequential() result(val)
            val = Sequential
        end function get_sequential

        subroutine set_sequential(val)
            logical, intent(in) :: val
            call ensure_settings_unlocked("set_sequential")
            Sequential = val
        end subroutine set_sequential

        logical function get_Langevin() result(val)
            val = Langevin
        end function get_Langevin

        subroutine set_Langevin(val)
            logical, intent(in) :: val
            call ensure_settings_unlocked("set_Langevin")
            Langevin = val
        end subroutine set_Langevin

        logical function get_HMC() result(val)
            val = HMC
        end function get_HMC

        subroutine set_HMC(val)
            logical, intent(in) :: val
            call ensure_settings_unlocked("set_HMC")
            HMC = val
        end subroutine set_HMC

        logical function get_Tempering_calc_det() result(val)
            val = Tempering_calc_det
        end function get_Tempering_calc_det

        subroutine set_Tempering_calc_det(val)
            logical, intent(in) :: val
            call ensure_settings_unlocked("set_Tempering_calc_det")
            Tempering_calc_det = val
        end subroutine set_Tempering_calc_det

        ! Character getters/setters
        character(len=64) function get_ham_name() result(val)
            val = ham_name
        end function get_ham_name

        subroutine set_ham_name(val)
            character(len=*), intent(in) :: val
            call ensure_settings_unlocked("set_ham_name")
            ham_name = val
        end subroutine set_ham_name
     

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Broadcast of the QMC runtime variables
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

#ifdef TEMPERING
        subroutine read_and_broadcast_TEMPERING_var()

          use iso_fortran_env, only: error_unit
          use runtime_error_mod, only: Terminate_on_error, ERROR_FILE_NOT_FOUND

          implicit none

          Integer :: ierr, irank

          mpi_per_parameter_set = 1   ! Default value
          Tempering_calc_det = .true. ! Default value
          
          CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
          IF (irank == 0) THEN
            OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
            IF (ierr /= 0) THEN
              WRITE(error_unit,*) 'main: unable to open <parameters>',ierr
              CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
            END IF
            READ(5,NML=VAR_TEMP)
            CLOSE(5)
          END IF
          CALL MPI_BCAST(N_exchange_steps        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(N_Tempering_frequency   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(mpi_per_parameter_set   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Tempering_calc_det      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

        end subroutine read_and_broadcast_TEMPERING_var
#endif
           

        subroutine read_and_broadcast_QMC_var_and_ham_name(Group_Comm)

          use iso_fortran_env, only: error_unit
          use runtime_error_mod, only: Terminate_on_error, ERROR_FILE_NOT_FOUND

          implicit none

          Integer, intent(in) :: Group_Comm

          integer :: ierr
          Character (len=64) :: file_para
#ifdef MPI
          Integer :: MPI_COMM_i, irank, irank_g, Isize_g, igroup

          CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, Isize_g, ierr)
          igroup           = irank/isize_g

#ifdef PARALLEL_PARAMS
          MPI_COMM_i = Group_Comm
          If ( irank_g == 0 ) then
             write(file_para,'(A,I0,A)') "Temp_", igroup, "/parameters"
#else
          MPI_COMM_i = MPI_COMM_WORLD
          If ( Irank == 0 ) then
             file_para = "parameters"
#endif
#else
             file_para = "parameters"
#endif
      
             Call set_QMC_runtime_default_var()
             OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
               WRITE(error_unit,*) 'main: unable to open <parameters>', file_para, ierr
               CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF
             READ(5,NML=VAR_QMC)
             REWIND(5)
             READ(5,NML=VAR_HAM_NAME)
             CLOSE(5)
#ifdef MPI
           Endif
           call broadcast_QMC_runtime_var(MPI_COMM_i)
           CALL MPI_BCAST(ham_name,64,MPI_CHARACTER,0,MPI_COMM_i,ierr)
#endif

        end subroutine read_and_broadcast_QMC_var_and_ham_name
        
end Module QMC_runtime_var