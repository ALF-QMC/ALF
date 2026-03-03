!  Copyright (C) 2025-2026 The ALF project
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
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief OpenMP thread configuration validation and safety module.
!
!> @details
!> This module provides functionality to check and configure OpenMP thread
!> settings at runtime. It ensures that parallel sections have explicit thread
!> counts, preventing unexpected behavior from unset environment variables.
!>
!> The main purpose is to provide a safety mechanism: if OMP_NUM_THREADS is not
!> explicitly set by the user, the module defaults to single-threaded execution
!> rather than using all available cores (which may not be desired for all
!> simulation scenarios).
!
!> @note
!> This module is only available when ALF is compiled with OpenMP support
!> (_OPENMP preprocessor flag). It should be called early in program initialization
!> before any parallel regions are executed.
!
!> @see
!> OpenMP specification: https://www.openmp.org/specifications/
!--------------------------------------------------------------------
    use omp_lib
    implicit none
    private
    public :: check_omp_num_threads

    contains
    
    subroutine check_omp_num_threads()
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Validates OMP_NUM_THREADS environment variable and sets safe defaults.
!
!> @details
!> This subroutine checks whether the OMP_NUM_THREADS environment variable
!> has been set by the user:
!> - If set: Reports the value to stdout (no modification)
!> - If unset: Sets the number of OpenMP threads to 1 (single-threaded mode)
!>   and reports this action to stdout
!>
!> This ensures predictable behavior and prevents accidental use of all
!> available CPU cores, which could interfere with other running jobs on
!> shared systems or exceed requested resources in HPC environments.
!
!> @note
!> Should be called before any parallel regions are encountered in the code.
!> The environment variable OMP_NUM_THREADS takes precedence if set.
!
!> @warning
!> This subroutine modifies the global OpenMP thread count via
!> omp_set_num_threads() if OMP_NUM_THREADS is not set.
!--------------------------------------------------------------------
        character(len=64) :: value   ! Buffer to store environment variable value
        integer :: length            ! Length of retrieved environment variable
            
        ! Query the OMP_NUM_THREADS environment variable
        call get_environment_variable("OMP_NUM_THREADS", value, length)
        
        if (length > 0) then
            ! Environment variable is set - respect user's choice
            print*, "OMP_NUM_THREADS set to ", trim(value)
        else
            ! Environment variable not set - use safe default of 1 thread
            print*, "OMP_NUM_THREADS unset, setting num_threads to 1"
            call omp_set_num_threads(1)
        endif
    end subroutine check_omp_num_threads
end module check_omp_num_threads_mod
#endif
