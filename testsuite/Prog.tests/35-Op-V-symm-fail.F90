! Test that a non-commuting, non-palindromic operator sequence is detected.
! Creates two 2x2 operators sharing a site with non-commuting O matrices,
! in a non-palindromic arrangement.
! Op_V_is_symmetric should return .false., then we STOP 1.
! In CMake, mark with WILL_FAIL TRUE.
!
Program Test_Op_V_symm_fail
  Use Operator_mod
  Implicit None

  Type (Operator), allocatable :: Op_V(:,:)
  Integer :: Ndim

  Ndim = 4

  ! Two non-commuting operators sharing site 2, in non-palindromic order
  allocate(Op_V(2, 1))

  ! Op1: acts on sites {1, 2} with real off-diagonal (sigma_x)
  Call Op_make(Op_V(1, 1), 2)
  Op_V(1, 1)%P(1) = 1
  Op_V(1, 1)%P(2) = 2
  Op_V(1, 1)%O(1,2) = cmplx(1.d0, 0.d0, kind(0.d0))
  Op_V(1, 1)%O(2,1) = cmplx(1.d0, 0.d0, kind(0.d0))
  Op_V(1, 1)%g      = cmplx(0.5d0, 0.d0, kind(0.d0))
  Op_V(1, 1)%alpha   = cmplx(0.d0, 0.d0, kind(0.d0))
  Op_V(1, 1)%type    = 2

  ! Op2: acts on sites {2, 3} with imaginary off-diagonal (sigma_y)
  Call Op_make(Op_V(2, 1), 2)
  Op_V(2, 1)%P(1) = 2
  Op_V(2, 1)%P(2) = 3
  Op_V(2, 1)%O(1,2) = cmplx(0.d0, 1.d0, kind(0.d0))
  Op_V(2, 1)%O(2,1) = cmplx(0.d0, -1.d0, kind(0.d0))
  Op_V(2, 1)%g      = cmplx(0.3d0, 0.d0, kind(0.d0))
  Op_V(2, 1)%alpha   = cmplx(0.d0, 0.d0, kind(0.d0))
  Op_V(2, 1)%type    = 2

  ! These don't commute: they share site 2 and have non-commuting O matrices.
  ! Op_V_is_symmetric should return .false.
  if (Op_V_is_symmetric(Op_V, 1, Ndim)) then
     ! Incorrectly reported as symmetric — normal exit triggers WILL_FAIL failure
     write(*,*) 'FAIL: non-commuting non-palindromic NOT detected'
  else
     ! Correctly detected asymmetry — exit with non-zero code
     write(*,*) 'Asymmetry correctly detected (expected)'
     STOP 1
  endif

end Program Test_Op_V_symm_fail
