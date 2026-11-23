! Unit test for qmc_walker_mod: Walker_layout, Walker_alloc and Walker_dealloc
! on a small model set with Set_model_for_testing (Ndim = 4, N_FL = 2,
! Ltrot = 10, Nwrap = 3; finite temperature, no projector).
! This program should exit successfully (return code 0).
Program Test_walker
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main, only: Set_model_for_testing, Setup_flavor_map, Clear_Ham_state, N_FL_eff
  Use Operator_mod, only: Operator
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use iso_fortran_env, only: error_unit
  Implicit none

  Type (qmc_walker_t) :: walker
  Type (Operator), dimension(:,:), allocatable :: Op_T, Op_V
  Integer :: n_fail, nf_eff, n
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
#endif
  n_fail = 0

  allocate(Op_T(0,2), Op_V(0,2))
  Call Set_model_for_testing(4, 2, 1, 10, 0, .false., .false., Op_T, Op_V)
  Call Setup_flavor_map()

  ! Layout: 10 time slices in intervals of 3 -> NSTM = 4, Stab_nt = 0 3 6 9 10
  Call Walker_layout(walker, 10, 3)
  Call check(walker%NSTM == 4, "NSTM")
  Call check(all(walker%Stab_nt == [0, 3, 6, 9, 10]), "Stab_nt")
  Call check(lbound(walker%Stab_nt,1) == 0, "lbound Stab_nt")

  Call Walker_alloc(walker)
  Call check(all(shape(walker%GR) == [4, 4, 2]), "shape GR")
  Call check(all(shape(walker%GR_Tilde) == [4, 4, 2]), "shape GR_Tilde")
  Call check(all(shape(walker%Test) == [4, 4]), "shape Test")
  Call check(size(walker%Phase_array) == 2, "size Phase_array")
  Call check(N_FL_eff == 2, "N_FL_eff")
  Call check(size(walker%udvl) == N_FL_eff .and. size(walker%udvr) == N_FL_eff, "size udvl, udvr")
  Call check(all(shape(walker%udvst) == [4, N_FL_eff]), "shape udvst")
  do nf_eff = 1, N_FL_eff
     Call check(walker%udvl(nf_eff)%side == 'L' .or. walker%udvl(nf_eff)%side == 'l', "side udvl")
     Call check(walker%udvr(nf_eff)%side == 'R' .or. walker%udvr(nf_eff)%side == 'r', "side udvr")
     do n = 1, walker%NSTM
        Call check(allocated(walker%udvst(n,nf_eff)%U), "udvst allocated")
     enddo
  enddo

  Call Walker_dealloc(walker)
  Call check(.not. allocated(walker%GR) .and. .not. allocated(walker%udvst) .and. &
       &     .not. allocated(walker%Stab_nt) .and. .not. allocated(walker%Phase_array), "deallocated")
  Call Clear_Ham_state()

  If (n_fail > 0) then
     Write(error_unit,*) "walker test: ", n_fail, " failure(s)"
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  endif
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  Write(*,*) "walker test passed"

contains

  Subroutine check(ok, label)
    Logical, intent(in) :: ok
    Character(len=*), intent(in) :: label
    If (.not. ok) then
       Write(error_unit,*) "FAIL ", label
       n_fail = n_fail + 1
    endif
  end Subroutine check

End Program Test_walker
