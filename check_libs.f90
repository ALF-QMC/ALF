! Small check if BLAS and LAPACK libraries are present and working.
program check_libs

  use iso_fortran_env
  integer, parameter :: N=2, ix=1, iy=1, LDA=2, LWORK=100
  integer :: INFO
  real    (kind=kind(0.d0)) :: C(2), D(2), dot, W(2), RWORK(4)
  complex (kind=kind(0.d0)) :: A(2,2), WORK(LWORK)
  real    (kind=kind(0.d0)), external :: ddot
  data    A /(0,0),(0,1),(0,-1),(0,0)/,  C /3,1415/, D /42,0/
  external ZHEEV
  
  dot = ddot(n,C,ix,D,iy)
  if (.not. abs(dot-126d0) < 1d-10) then
    write(error_unit, *) "Problem with BLAS"
    error stop
  endif

  call ZHEEV( 'V', 'U', N, A, LDA, W, WORK, LWORK, RWORK, INFO )
  if ( .not. (abs(W(1)+1d0) < 1d-10 .and. abs(W(1)+1d0) < 1d-10) ) then
    write(error_unit, *) "Problem with LAPACK"
    error stop
  endif
  
end program check_libs
