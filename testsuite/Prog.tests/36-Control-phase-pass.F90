! The phase guard of the Control module must let finite phases of unit modulus pass.
Program Test_Control_phase_pass
  Use Control
  Implicit None

  Call Control_Check_Phase(cmplx(0.6d0, 0.8d0, kind(0.d0)))
  Call Control_Check_Phase(cmplx(-1.d0, 0.d0, kind(0.d0)))
  Call Control_PrecisionP(cmplx(0.6d0, 0.8d0, kind(0.d0)), cmplx(0.6d0, 0.8d0, kind(0.d0)))

  Write(*,*) 'success'
End Program Test_Control_phase_pass
