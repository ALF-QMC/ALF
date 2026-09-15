! Test that Op_set rejects a 1x1 operator with non-real (complex) diagonal.
! This program is EXPECTED TO FAIL (error stop) — use WILL_FAIL in CTest.
!
Program Test_Op_validate_1x1_fail
  Use Operator_mod
  Implicit None

  Type (Operator) :: Op

  Call Op_make(Op, 1)
  Op%O(1,1) = cmplx(1.d0, 0.5d0, kind(0.d0))  ! Complex diagonal — not Hermitian
  Op%P(1) = 1
  Op%g     = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = 1

  ! This call should trigger Terminate_on_error due to complex diagonal
  Call Op_set(Op)

  ! Should never reach here
  Write(*,*) 'ERROR: Op_set did not reject non-Hermitian 1x1 operator!'

End Program Test_Op_validate_1x1_fail
