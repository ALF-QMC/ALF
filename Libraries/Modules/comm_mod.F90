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
      public :: comm_world_rank, comm_world_size
      public :: comm_reduce_sum, comm_reduce_max, comm_bcast

      !> MPI_SUM to rank 0; the result is only meaningful on rank 0.
      interface comm_reduce_sum
         module procedure reduce_sum_r8, reduce_sum_i8
      end interface comm_reduce_sum

      !> MPI_MAX to rank 0; the result is only meaningful on rank 0.
      interface comm_reduce_max
         module procedure reduce_max_r8
      end interface comm_reduce_max

      !> Broadcast from rank 0.
      interface comm_bcast
         module procedure bcast_r8
      end interface comm_bcast

      public :: comm_reduce_mean, comm_reduce_mean_array

      !> MPI_SUM to rank 0, then division by the number of processes in comm.
      !> The result is only meaningful on rank 0.
      interface comm_reduce_mean
         module procedure reduce_mean_r8
      end interface comm_reduce_mean

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


      !> Rank in MPI_COMM_WORLD (0 without MPI).
      Integer function comm_world_rank()
#ifdef MPI
        Integer :: ierr
        call MPI_Comm_rank(MPI_COMM_WORLD, comm_world_rank, ierr)
#else
        comm_world_rank = 0
#endif
      end function comm_world_rank

      !> Number of processes in MPI_COMM_WORLD (1 without MPI).
      Integer function comm_world_size()
#ifdef MPI
        Integer :: ierr
        call MPI_Comm_size(MPI_COMM_WORLD, comm_world_size, ierr)
#else
        comm_world_size = 1
#endif
      end function comm_world_size

      Subroutine reduce_sum_r8(x, comm)
        Real (Kind=Kind(0.d0)), intent(inout) :: x
        Integer, intent(in) :: comm
#ifdef MPI
        Real (Kind=Kind(0.d0)) :: tmp
        Integer :: ierr
        tmp = 0.d0
        call MPI_Reduce(x, tmp, 1, MPI_REAL8, MPI_SUM, 0, comm, ierr)
        x = tmp
#endif
      end Subroutine reduce_sum_r8

      Subroutine reduce_sum_i8(n, comm)
        Integer (Kind=8), intent(inout) :: n
        Integer, intent(in) :: comm
#ifdef MPI
        Integer (Kind=8) :: tmp
        Integer :: ierr
        tmp = 0
        call MPI_Reduce(n, tmp, 1, MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
        n = tmp
#endif
      end Subroutine reduce_sum_i8

      Subroutine reduce_max_r8(x, comm)
        Real (Kind=Kind(0.d0)), intent(inout) :: x
        Integer, intent(in) :: comm
#ifdef MPI
        Real (Kind=Kind(0.d0)) :: tmp
        Integer :: ierr
        tmp = x
        call MPI_Reduce(x, tmp, 1, MPI_REAL8, MPI_MAX, 0, comm, ierr)
        x = tmp
#endif
      end Subroutine reduce_max_r8

      Subroutine bcast_r8(x, comm)
        Real (Kind=Kind(0.d0)), intent(inout) :: x
        Integer, intent(in) :: comm
#ifdef MPI
        Integer :: ierr
        call MPI_Bcast(x, 1, MPI_REAL8, 0, comm, ierr)
#endif
      end Subroutine bcast_r8


      Subroutine reduce_mean_r8(x, comm)
        Real (Kind=Kind(0.d0)), intent(inout) :: x
        Integer, intent(in) :: comm
#ifdef MPI
        Real (Kind=Kind(0.d0)) :: tmp
        Integer :: ierr
        tmp = 0.d0
        call MPI_Reduce(x, tmp, 1, MPI_REAL8, MPI_SUM, 0, comm, ierr)
        x = tmp/dble(comm_size(comm))
#endif
      end Subroutine reduce_mean_r8

      !> Same as comm_reduce_mean for a complex array of any rank, passed with its
      !> number of elements as in MPI (not a generic: a generic name would require
      !> a fixed array rank).
      Subroutine comm_reduce_mean_array(buf, n, comm)
        Complex (Kind=Kind(0.d0)), intent(inout) :: buf(*)
        Integer, intent(in) :: n, comm
#ifdef MPI
        Complex (Kind=Kind(0.d0)), allocatable :: tmp(:)
        Integer :: ierr
        allocate(tmp(n))
        tmp = cmplx(0.d0, 0.d0, kind(0.d0))
        call MPI_Reduce(buf, tmp, n, MPI_COMPLEX16, MPI_SUM, 0, comm, ierr)
        buf(1:n) = tmp/dble(comm_size(comm))
        deallocate(tmp)
#endif
      end Subroutine comm_reduce_mean_array

    end Module comm_mod
