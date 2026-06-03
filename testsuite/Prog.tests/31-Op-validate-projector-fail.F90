! Test that Op_set rejects a projector with an out-of-bounds index.
! This program is EXPECTED TO FAIL (error stop) — use WILL_FAIL in CTest.
!
Program Test_Op_validate_projector_fail
  Use Operator_mod
  Implicit None

  Type (Operator) :: Op

  Call Op_make(Op, 2)
  Op%O(1,1) = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%O(2,2) = cmplx(2.d0, 0.d0, kind(0.d0))
  Op%P(1) = 1
  Op%P(2) = 0    ! Invalid: projector index must be >= 1
  Op%g     = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = 1

  ! This call should trigger Terminate_on_error due to invalid projector
  Call Op_set(Op)

  ! Should never reach here
  Write(*,*) 'ERROR: Op_set did not reject invalid projector!'

End Program Test_Op_validate_projector_fail
