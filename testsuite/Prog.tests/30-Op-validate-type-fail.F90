! Test that Op_set rejects an invalid operator type (type=-1).
! This program is EXPECTED TO FAIL (error stop) — use WILL_FAIL in CTest.
!
Program Test_Op_validate_type_fail
  Use Operator_mod
  Implicit None

  Type (Operator) :: Op

  Call Op_make(Op, 2)
  Op%O(1,1) = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%O(2,2) = cmplx(2.d0, 0.d0, kind(0.d0))
  Op%P(1) = 1
  Op%P(2) = 2
  Op%g     = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = -1  ! Invalid type

  ! This call should trigger Terminate_on_error due to invalid type
  Call Op_set(Op)

  ! Should never reach here
  Write(*,*) 'ERROR: Op_set did not reject invalid operator type!'

End Program Test_Op_validate_type_fail
