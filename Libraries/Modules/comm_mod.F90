!  Copyright (C) 2026 The ALF project
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

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Thin communication layer. The callers pass the communicator explicitly and
!> do not need any #ifdef MPI: without MPI the routines behave like a single
!> rank (rank 0, size 1, root).
!--------------------------------------------------------------------

    Module comm_mod
#ifdef MPI
      Use mpi
#endif
      Implicit none
      private
      public :: comm_rank, comm_size, is_root

    contains

      !> Rank of the calling process in comm (0 without MPI).
      Integer function comm_rank(comm)
        Integer, intent(in) :: comm
#ifdef MPI
        Integer :: ierr
        call MPI_Comm_rank(comm, comm_rank, ierr)
#else
        comm_rank = 0 * comm
#endif
      end function comm_rank

      !> Number of processes in comm (1 without MPI).
      Integer function comm_size(comm)
        Integer, intent(in) :: comm
#ifdef MPI
        Integer :: ierr
        call MPI_Comm_size(comm, comm_size, ierr)
#else
        comm_size = 1 + 0 * comm
#endif
      end function comm_size

      !> True on rank 0 of comm (always true without MPI).
      Logical function is_root(comm)
        Integer, intent(in) :: comm
        is_root = (comm_rank(comm) == 0)
      end function is_root

    end Module comm_mod
