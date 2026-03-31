! Test that Op_set rejects a non-Hermitian operator matrix.
! This program is EXPECTED TO FAIL (error stop) — use WILL_FAIL in CTest.
!
Program Test_Op_validate_hermitian_fail
  Use Operator_mod
  Implicit None

  Type (Operator) :: Op
  Integer :: I, N

  N = 3
  Call Op_make(Op, N)
  ! Construct a non-Hermitian matrix: O(1,2) != conjg(O(2,1))
  Op%O(1,1) = cmplx( 1.d0, 0.d0, kind(0.d0))
  Op%O(1,2) = cmplx( 2.d0, 1.d0, kind(0.d0))
  Op%O(2,1) = cmplx( 2.d0, 1.d0, kind(0.d0))  ! Should be (2, -1) for Hermitian
  Op%O(2,2) = cmplx( 3.d0, 0.d0, kind(0.d0))
  Op%O(3,3) = cmplx( 1.d0, 0.d0, kind(0.d0))
  do I = 1, N
     Op%P(I) = I
  enddo
  Op%g     = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = 1

  ! This call should trigger Terminate_on_error due to non-Hermitian Op%O
  Call Op_set(Op)

  ! Should never reach here
  Write(*,*) 'ERROR: Op_set did not reject non-Hermitian operator!'

End Program Test_Op_validate_hermitian_fail
