! Unit test for comm_mod. Works for any number of MPI ranks; without MPI the
! routines must behave like a single rank.
Program Test_comm_mod

  Use comm_mod, only: comm_rank, comm_size, is_root
#ifdef MPI
  Use mpi
#endif
  Implicit none

  Integer :: n_fail, comm, rank_ref, size_ref
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  Call MPI_Comm_rank(comm, rank_ref, ierr)
  Call MPI_Comm_size(comm, size_ref, ierr)
#else
  comm = 0
  rank_ref = 0
  size_ref = 1
#endif
  n_fail = 0

  Call check_i("comm_rank", comm_rank(comm), rank_ref)
  Call check_i("comm_size", comm_size(comm), size_ref)
  Call check_l("is_root",   is_root(comm),   rank_ref == 0)

#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  If (n_fail > 0) then
     Write(*,*) "comm_mod test: ", n_fail, " failure(s)"
     Error stop 1
  endif
  If (rank_ref == 0) Write(*,*) "comm_mod test passed"

contains

  Subroutine check_i(label, got, expected)
    Character(len=*), intent(in) :: label
    Integer, intent(in) :: got, expected
    If (got /= expected) then
       Write(*,*) "FAIL ", label, ": got ", got, " expected ", expected
       n_fail = n_fail + 1
    endif
  end Subroutine check_i

  Subroutine check_l(label, got, expected)
    Character(len=*), intent(in) :: label
    Logical, intent(in) :: got, expected
    If (got .neqv. expected) then
       Write(*,*) "FAIL ", label, ": got ", got, " expected ", expected
       n_fail = n_fail + 1
    endif
  end Subroutine check_l

end Program Test_comm_mod
