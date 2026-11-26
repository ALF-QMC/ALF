! Control_PrecisionP must stop on a NaN phase. This program is EXPECTED TO FAIL (error stop) — use WILL_FAIL in CTest.
Program Test_Control_phase_precision_nan_fail
  Use Control
  Use iso_fortran_env, only: int64
  Implicit None

  Integer (Kind=int64) :: nanbits
  Complex (Kind=Kind(0.d0)) :: Z

  nanbits = int(z'7FF8000000000000', int64)
  Z = cmplx(transfer(nanbits, 1.d0), 0.d0, kind(0.d0))

  Call Control_PrecisionP(cmplx(1.d0, 0.d0, kind(0.d0)), Z)

  ! Should never reach here
  Write(*,*) 'ERROR: a NaN phase was not rejected by Control_PrecisionP!'
End Program Test_Control_phase_precision_nan_fail
