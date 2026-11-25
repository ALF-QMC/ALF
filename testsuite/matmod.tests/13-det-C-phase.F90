! Test of DET_C_LOGPHASE and DET_C_PHASE: phase and ln|det| must stay correct when the determinant itself under- or overflows.
! compile with
! gfortran -I ../../Libraries/Modules/ main.F90 ../../Libraries/Modules/modules_90.a -llapack -lblas

Program Main

  Use MyMats

  Implicit none

  Integer, parameter :: N = 400
  Complex (Kind=Kind(0.d0)) :: A(N,N), A2(2,2), A0(3,3), B(8,8), C(8,8), Phase, Phase_ref, Z
  Real (Kind=Kind(0.d0)) :: LogAbs, theta, thsum, pi
  Integer :: i, j

  pi = acos(-1.d0)

  ! 1. Upper triangular matrix with diagonal 1e-3 exp(i theta_i): |det| = 1e-1200 underflows double precision.
  !    The raw product (DET_C) gives 0, hence a phase of 0/0 = NaN; DET_C_LOGPHASE has to return the exact phase.
  A = 0.d0
  thsum = 0.d0
  do i = 1, N
     theta = 0.37d0*i
     thsum = thsum + theta
     A(i,i) = 1.d-3*cmplx(cos(theta), sin(theta), kind(0.d0))
     do j = i+1, N
        A(i,j) = 1.d-3*cmplx(1.d0, 0.5d0, kind(0.d0))   ! the strictly upper part does not change the determinant
     enddo
  enddo
  Call DET_C_LOGPHASE(A, N, LogAbs, Phase)
  Phase_ref = cmplx(cos(thsum), sin(thsum), kind(0.d0))
  if (abs(Phase - Phase_ref) > 1.d-8) then
     write (*,*) "ERROR: phase of the underflowing determinant", Phase, Phase_ref
     stop 2
  endif
  if (abs(LogAbs - N*log(1.d-3)) > 1.d-9*abs(N*log(1.d-3))) then
     write (*,*) "ERROR: ln|det| of the underflowing determinant", LogAbs, N*log(1.d-3)
     stop 2
  endif

  ! 2. A row permutation: det([[0,2],[3i,0]]) = -6i, phase -i, ln|det| = ln 6.
  A2 = 0.d0
  A2(1,2) = cmplx(2.d0, 0.d0, kind(0.d0))
  A2(2,1) = cmplx(0.d0, 3.d0, kind(0.d0))
  Call DET_C_LOGPHASE(A2, 2, LogAbs, Phase)
  if (abs(Phase - cmplx(0.d0, -1.d0, kind(0.d0))) > 1.d-13 .or. abs(LogAbs - log(6.d0)) > 1.d-13) then
     write (*,*) "ERROR: permutation sign", Phase, LogAbs
     stop 2
  endif

  ! 3. Exactly singular matrix: phase (0,0), ln|det| = -huge.
  A0 = 0.d0
  Call DET_C_LOGPHASE(A0, 3, LogAbs, Phase)
  if (abs(Phase) /= 0.d0 .or. LogAbs /= -huge(1.d0)) then
     write (*,*) "ERROR: singular matrix", Phase, LogAbs
     stop 2
  endif

  ! 4. Agreement with DET_C where DET_C is safe, and DET_C_PHASE equals the phase of DET_C_LOGPHASE.
  do j = 1, 8
     do i = 1, 8
        B(i,j) = cmplx(sin(1.3d0*i + 0.7d0*j*j), cos(0.9d0*i*j + 0.2d0), kind(0.d0))
     enddo
  enddo
  C = B
  Z = DET_C(C, 8)
  C = B
  Call DET_C_LOGPHASE(C, 8, LogAbs, Phase)
  if (abs(Z - Phase*exp(LogAbs)) > 1.d-10*abs(Z)) then
     write (*,*) "ERROR: DET_C_LOGPHASE disagrees with DET_C", Z, Phase*exp(LogAbs)
     stop 2
  endif
  C = B
  if (abs(DET_C_PHASE(C, 8) - Phase) > 1.d-13) then
     write (*,*) "ERROR: DET_C_PHASE", DET_C_PHASE(C, 8), Phase
     stop 2
  endif

  write (*, *) "success"
end Program Main
