! A phase of modulus 0 (a determinant that underflowed) must stop the run.
! This program is EXPECTED TO FAIL (error stop) — use WILL_FAIL in CTest.
Program Test_Control_phase_modulus_fail
  Use Control
  Implicit None

  Call Control_Check_Phase(cmplx(0.d0, 0.d0, kind(0.d0)))

  ! Should never reach here
  Write(*,*) 'ERROR: a phase of modulus 0 was not rejected!'
End Program Test_Control_phase_modulus_fail
