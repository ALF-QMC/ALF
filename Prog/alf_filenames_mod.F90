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


module alf_filenames_mod
    use alf_mpi_mod
    Implicit none
    private
    public :: alf_filenames_type, alf_filenames

    Character (len=64), save :: file_seeds, file_dat, file_info
    Character (len=64), save :: file_parameters_main, file_parameters_qmc, file_parameters_hamilton

    type alf_filenames_type
    contains
        procedure, nopass :: init
        procedure, nopass :: get_file_info
    end type alf_filenames_type

    type(alf_filenames_type) :: alf_filenames
  contains


    Character(len=64) function get_file_info()
        get_file_info = file_info
    end function get_file_info

    Character(len=64) function get_file_info()
        get_file_info = file_info
    end function get_file_info

    Character(len=64) function get_file_info()
        get_file_info = file_info
    end function get_file_info


    subroutine init()
        file_parameters_main = "parameters"
        file_parameters_qmc = "parameters"
        file_parameters_hamilton = "parameters"
        file_info = "info"
        file_seeds = "seeds"

#if defined(TEMPERING)
        write(file_parameters_hamilton,'(A,I0,A)') "Temp_", igroup, "/parameters", 
        write(file_info,'(A,I0,A)') "Temp_", igroup, "/info"
#endif

#if defined(PARALLEL_PARAMS)
        write(file_parameters_qmc,'(A,I0,A)') "Temp_", igroup, "/parameters"
#endif
    end subroutine init


end module alf_filenames_mod