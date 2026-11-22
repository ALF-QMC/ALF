! Unit test for comm_mod. Works for any number of MPI ranks; without MPI the
! routines must behave like a single rank.
Program Test_comm_mod

  Use comm_mod, only: comm_rank, comm_size, is_root, comm_world, comm_world_rank, comm_world_size, &
       &              comm_reduce_sum, comm_reduce_max, comm_bcast, comm_reduce_mean, &
       &              comm_reduce_mean_array
#ifdef MPI
  Use mpi
#endif
  Implicit none

  Integer :: n_fail, comm, rank_ref, size_ref
  Real (Kind=Kind(0.d0)) :: x
  Integer (Kind=8) :: n8
  Complex (Kind=Kind(0.d0)) :: c(2,3)
  Integer :: i, j
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
  Call check_i("comm_world_rank", comm_world_rank(), rank_ref)
  Call check_i("comm_world_size", comm_world_size(), size_ref)
  Call check_i("comm_size(comm_world())", comm_size(comm_world()), size_ref)

  ! Rank r contributes r+1: sum = n(n+1)/2, max = n on rank 0.
  x = dble(rank_ref + 1)
  Call comm_reduce_sum(x, comm)
  If (rank_ref == 0) Call check_r("comm_reduce_sum r8", x, dble(size_ref*(size_ref+1)/2))
  n8 = int(rank_ref + 1, 8)
  Call comm_reduce_sum(n8, comm)
  If (rank_ref == 0) Call check_i("comm_reduce_sum i8", int(n8), size_ref*(size_ref+1)/2)
  x = dble(rank_ref + 1)
  Call comm_reduce_max(x, comm)
  If (rank_ref == 0) Call check_r("comm_reduce_max r8", x, dble(size_ref))

  ! Means over ranks: scalar and a rank-2 complex array (passed with its size).
  x = dble(rank_ref + 1)
  Call comm_reduce_mean(x, comm)
  If (rank_ref == 0) Call check_r("comm_reduce_mean r8", x, dble(size_ref+1)/2.d0)
  Do j = 1, 3
     Do i = 1, 2
        c(i,j) = cmplx(dble(rank_ref + 1), dble(i + 10*j), kind(0.d0))
     enddo
  enddo
  Call comm_reduce_mean_array(c, size(c), comm)
  If (rank_ref == 0) then
     Do j = 1, 3
        Do i = 1, 2
           Call check_r("comm_reduce_mean_array re", dble(c(i,j)), dble(size_ref+1)/2.d0)
           Call check_r("comm_reduce_mean_array im", aimag(c(i,j)), dble(i + 10*j))
        enddo
     enddo
  endif

  ! Broadcast from rank 0 reaches every rank.
  x = 0.d0
  If (rank_ref == 0) x = 42.d0
  Call comm_bcast(x, comm)
  Call check_r("comm_bcast r8", x, 42.d0)

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

  Subroutine check_r(label, got, expected)
    Character(len=*), intent(in) :: label
    Real (Kind=Kind(0.d0)), intent(in) :: got, expected
    If (got /= expected) then
       Write(*,*) "FAIL ", label, ": got ", got, " expected ", expected
       n_fail = n_fail + 1
    endif
  end Subroutine check_r

  Subroutine check_l(label, got, expected)
    Character(len=*), intent(in) :: label
    Logical, intent(in) :: got, expected
    If (got .neqv. expected) then
       Write(*,*) "FAIL ", label, ": got ", got, " expected ", expected
       n_fail = n_fail + 1
    endif
  end Subroutine check_l

end Program Test_comm_mod
