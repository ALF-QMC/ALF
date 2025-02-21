!  Copyright (C) 2025 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

#ifdef _OPENMP
module check_omp_num_threads_mod
    use omp_lib
    implicit none
    private
    public :: check_omp_num_threads

    contains
    subroutine check_omp_num_threads()
        character(len=64) :: value
        integer :: length
            
        call get_environment_variable("OMP_NUM_THREADS", value, length)
        if (length > 0) then
            print*, "OMP_NUM_THREADS set to ", trim(value)
        else
            print*, "OMP_NUM_THREADS unset, setting num_threads to 1"
            call omp_set_num_threads(1)
        endif
    end subroutine check_omp_num_threads
end module check_omp_num_threads_mod
#endif
