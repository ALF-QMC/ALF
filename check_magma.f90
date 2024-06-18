! Small check if BLAS and LAPACK libraries are present and working.
program check_libs
  
  use iso_fortran_env
  use magma
  implicit none

  integer, parameter :: N=2, ix=1, iy=1, LDA=2
  integer :: INFO, LWORK
  real    (kind=kind(0.d0)) :: C(2), D(2), dot, W(N), RWORK(7*N)
  complex (kind=kind(0.d0)) :: A(2,2)
  complex (kind=kind(0.d0)), allocatable :: WORK(:)
  real    (kind=kind(0.d0)), external :: ddot
  data    A /(0,0),(0,1),(0,-1),(0,0)/,  C /3,1415/, D /42,0/
  !external magmaf_ZHEEVx
  ! external ZHEEVx

  real(kind=kind(0.d0)), parameter     :: vl=0.d0, vu=0.d0  ! Not referenced
  real(kind=kind(0.d0))      :: abstol
  integer :: m(1), ldz=N, iwork(10), ifail(N)
  integer, parameter :: il=0, iu=0  ! Not referenced
  complex (kind=kind(0.d0)) :: z(N,N)

  real :: foo
  complex (kind=kind(0.d0)) :: lwork_out(1)

  !abstol=2*lapackf77_dlamch('S')
  abstol=0.000001d0
  abstol=0.d0
  
  ! dot = magmaf_ddot(n,C,ix,D,iy)
  ! if (.not. abs(dot-126d0) < 1d-10) then
  !   write(error_unit, *) "Problem with BLAS"
  !   error stop
  ! endif
  print*, "foo"

  LWORK = -1
  
  call magmaf_ZHEEVx(MagmaVec, MagmaRangeAll, MagmaUpper, N, A, LDA, vl, vu, il, iu, abstol, m, W, Z, ldz, &
  lwork_out, LWORK, RWORK, iwork, ifail, INFO )
  
  LWORK = lwork_out(1)
  print*, "bar", lwork
  allocate(WORK(LWORK))
  call magmaf_ZHEEVx(MagmaVec, MagmaRangeAll, MagmaUpper, N, A, LDA, vl, vu, il, iu, abstol, m, W, Z, ldz, &
      WORK, LWORK, RWORK, iwork, ifail, INFO )
  if ( .not. (abs(W(1)+1d0) < 1d-10 .and. abs(W(1)+1d0) < 1d-10) ) then
    write(error_unit, *) "Problem with LAPACK"
    error stop
  endif
  
end program check_libs
