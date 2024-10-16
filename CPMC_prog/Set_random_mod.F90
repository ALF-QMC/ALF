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
