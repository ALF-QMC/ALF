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
!     use mpi_shared_memory
    use runtime_error_mod
#ifdef MPI
    Use mpi_f08
#endif
        
    Implicit none
    private
    public :: is_main_process, is_group_main_process, get_mpi_per_parameter_set
    public :: get_igroup, get_irank_g, get_irank, get_isize, get_isize_g
    public :: alf_mpi_init1, alf_mpi_init2
#ifdef MPI
    public :: get_group_comm, get_MPI_COMM_QMC
#endif

    integer, save :: mpi_per_parameter_set
    Integer, save :: isize, irank
#ifdef MPI
    TYPE(MPI_Comm), save :: Group_Comm, MPI_COMM_QMC
#endif
    Integer, save :: irank_g, color, key, igroup

#ifdef TEMPERING
    integer, save :: N_exchange_steps, N_Tempering_frequency
    logical, save :: Tempering_calc_det
#endif

  contains

    logical function is_main_process()
        is_main_process = (Irank == 0)
    end function is_main_process

    logical function is_group_main_process()
        is_group_main_process = (irank_g == 0)
    end function is_group_main_process

    integer function get_mpi_per_parameter_set()
        get_mpi_per_parameter_set = mpi_per_parameter_set
    end function get_mpi_per_parameter_set

#ifdef MPI
    TYPE(MPI_Comm) function get_group_comm()
        get_group_comm = Group_Comm
    end function get_group_comm

    TYPE(MPI_Comm) function get_MPI_COMM_QMC()
        get_MPI_COMM_QMC = MPI_COMM_QMC
    end function get_MPI_COMM_QMC
#endif

    integer function get_igroup()
        get_igroup = igroup
    end function get_igroup

    integer function get_irank_g()
        get_irank_g = irank_g
    end function get_irank_g

    integer function get_irank()
        get_irank = irank
    end function get_irank

    integer function get_isize()
        get_isize = isize
    end function get_isize

    integer function get_isize_g()
        get_isize_g = mpi_per_parameter_set
    end function get_isize_g



    subroutine alf_mpi_init1()
#ifdef MPI
        CALL MPI_INIT()
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK)
#else
        isize = 1
        irank = 0
#endif
    end subroutine alf_mpi_init1



    subroutine alf_mpi_init2()

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
           WRITE(error_unit,*) 'Unable to open <parameters>',ierr
           CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
        END IF
        READ(5,NML=VAR_TEMP)
        CLOSE(5)
        CALL MPI_BCAST(N_exchange_steps        ,1,MPI_INTEGER,0,MPI_COMM_WORLD)
        CALL MPI_BCAST(N_Tempering_frequency   ,1,MPI_INTEGER,0,MPI_COMM_WORLD)
        CALL MPI_BCAST(mpi_per_parameter_set   ,1,MPI_INTEGER,0,MPI_COMM_WORLD)
        CALL MPI_BCAST(Tempering_calc_det      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD)
        if ( mod(ISIZE,mpi_per_parameter_set) .ne. 0 ) then
           Write (error_unit,*) "mpi_per_parameter_set is not a multiple of total mpi processes"
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
#elif !defined(TEMPERING) && defined(MPI)
        mpi_per_parameter_set = Isize
#elif !defined(TEMPERING) && !defined(MPI)
        mpi_per_parameter_set = 1
#endif

#if defined(MPI)
        color = irank/mpi_per_parameter_set
        key   =  0
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Group_comm)
        call MPI_Comm_rank(Group_Comm, Irank_g)
        igroup           = irank/mpi_per_parameter_set
        MPI_COMM_QMC = MPI_COMM_WORLD
#else
        irank = 0
        irank_g = 0
        igroup = 0
#endif

#if defined(PARALLEL_PARAMS)
        MPI_COMM_QMC = Group_Comm
#endif
    end subroutine alf_mpi_init2

end module alf_mpi_mod
