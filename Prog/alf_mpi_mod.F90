!  Copyright (C) 2023 The ALF project
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


module alf_mpi_mod
    use iso_fortran_env, only: output_unit, error_unit
    use mpi_shared_memory
#ifdef MPI
    Use mpi
#endif
        
    Implicit none
    private
    public :: alf_mpi_type, alf_mpi, Group_comm

    integer, save :: mpi_per_parameter_set
    Integer, save :: isize, irank
    Integer, save :: Group_Comm
#ifdef MPI
    Integer, save :: irank_g, color, key, igroup, MPI_COMM_i
#endif

    type alf_mpi_type
    contains
        procedure, nopass :: init1
        procedure, nopass :: init2
        procedure, nopass :: is_main_process
        procedure, nopass :: get_mpi_per_parameter_set
    end type alf_mpi_type

    type(alf_mpi_type) :: alf_mpi
  contains

    logical function is_main_process()
        is_main_process = (Irank == 0)
    end function is_main_process

    integer function get_mpi_per_parameter_set()
        get_mpi_per_parameter_set = mpi_per_parameter_set
    end function get_mpi_per_parameter_set



    subroutine init1()
        Integer :: ierr
#ifdef MPI
        CALL MPI_INIT(ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#else
        isize = 1
        irank = 0
#endif
    end subroutine init1



    subroutine init2(N_exchange_steps, N_Tempering_frequency, Tempering_calc_det)
        integer, intent(out) :: N_exchange_steps, N_Tempering_frequency
        logical, intent(out) :: Tempering_calc_det

        integer :: ierr

#if defined(TEMPERING) && !defined(MPI)
#error Mpi has to be defined for tempering runs!
#endif
#if defined(PARALLEL_PARAMS) && !defined(TEMPERING)
#error TEMPERING has to be defined for PARALLEL_PARAMS!
#endif

#if defined(TEMPERING)
        NAMELIST /VAR_TEMP/  N_exchange_steps, N_Tempering_frequency, mpi_per_parameter_set, Tempering_calc_det
#endif

#if defined(TEMPERING)
        mpi_per_parameter_set = 1   ! Default value
        Tempering_calc_det = .true. ! Default value
        OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(error_unit,*) 'main: unable to open <parameters>',ierr
           CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
        END IF
        READ(5,NML=VAR_TEMP)
        CLOSE(5)
        CALL MPI_BCAST(N_exchange_steps        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(N_Tempering_frequency   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(mpi_per_parameter_set   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Tempering_calc_det      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if ( mod(ISIZE,mpi_per_parameter_set) .ne. 0 ) then
           Write (error_unit,*) "mpi_per_parameter_set is not a multiple of total mpi processes"
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        Call Global_Tempering_setup
#elif !defined(TEMPERING) && defined(MPI)
        mpi_per_parameter_set = Isize
#elif !defined(TEMPERING) && !defined(MPI)
        mpi_per_parameter_set = 1
#elif defined(TEMPERING)  && !defined(MPI)
        Write(error_unit,*) 'Mpi has to be defined for tempering runs'
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
#if defined(PARALLEL_PARAMS) && !defined(TEMPERING)
        Write(error_unit,*) 'TEMPERING has to be defined for PARALLEL_PARAMS'
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif

#if defined(MPI)
        color = irank/mpi_per_parameter_set
        key   =  0
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Group_comm, ierr)
        call MPI_Comm_rank(Group_Comm, Irank_g, ierr)
        igroup           = irank/mpi_per_parameter_set
        MPI_COMM_QMC = MPI_COMM_WORLD
#endif

#if defined(PARALLEL_PARAMS)
        MPI_COMM_QMC = Group_Comm
#endif
    end subroutine init2

end module alf_mpi_mod
