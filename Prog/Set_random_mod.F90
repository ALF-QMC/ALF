!  Copyright (C) 2016 - 2018 The ALF project
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

module set_random
   use runtime_error_mod
   implicit none
contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This  routine reads the seeds from the file File_seeds,  distributes them to mpi  processes and
!> sets the random number generator.
!
!--------------------------------------------------------------------

   subroutine Set_Random_number_Generator(File_seeds, Seed_in)

#ifdef MPI
      use mpi
#endif
      use Random_Wrap

      use iso_fortran_env, only: output_unit, error_unit

      implicit none

      character(LEN=64), intent(IN) :: File_seeds
      integer, intent(out) :: SEED_IN
      integer :: I, IERR
      integer, allocatable :: SEED_VEC(:)

#ifdef MPI
      integer        :: STATUS(MPI_STATUS_SIZE), irank_g, isize_g, igroup, ISIZE, IRANK
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
#endif

#if defined(MPI)
      if (IRANK == 0) then
         open (UNIT=5, FILE=File_seeds, STATUS='OLD', ACTION='READ', IOSTAT=IERR)
         if (IERR /= 0) then
            write (error_unit, *) 'Fields_in: unable to open <seeds>', IERR
            call Terminate_on_error(ERROR_FILE_NOT_FOUND, __FILE__, __LINE__)
         end if
         do I = ISIZE - 1, 1, -1
            read (5, *) SEED_IN
            call MPI_SEND(SEED_IN, 1, MPI_INTEGER, I, I + 1024, MPI_COMM_WORLD, IERR)
         end do
         read (5, *) SEED_IN
         close (5)
      else
         call MPI_RECV(SEED_IN, 1, MPI_INTEGER, 0, IRANK + 1024, MPI_COMM_WORLD, STATUS, IERR)
      end if
      allocate (SEED_VEC(1))
      SEED_VEC(1) = SEED_IN
      call RANSET(SEED_VEC)
      deallocate (SEED_VEC)
#else
      open (UNIT=5, FILE=FILE_seeds, STATUS='OLD', ACTION='READ', IOSTAT=IERR)
      if (IERR /= 0) then
         write (error_unit, *) 'Fields_in: unable to open <seeds>', IERR
         call Terminate_on_error(ERROR_FILE_NOT_FOUND, __FILE__, __LINE__)
      end if
      read (5, *) SEED_IN
      close (5)
      allocate (SEED_VEC(1))
      SEED_VEC(1) = SEED_IN
      call RANSET(SEED_VEC)
      deallocate (SEED_VEC)
#endif

   end subroutine Set_Random_number_Generator

end module
