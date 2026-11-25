! A NaN phase must stop the run. This program is EXPECTED TO FAIL (error stop) — use WILL_FAIL in CTest.
! The NaN is built from its bit pattern; the test is compiled with the flags of the program (-ffast-math), under which a guard
! written as x /= x or x > y can be folded away by the compiler.
Program Test_Control_phase_nan_fail
  Use Control
  Use iso_fortran_env, only: int64
  Implicit None

  Integer (Kind=int64) :: nanbits
  Complex (Kind=Kind(0.d0)) :: Z

  nanbits = int(z'7FF8000000000000', int64)
  Z = cmplx(transfer(nanbits, 1.d0), 0.d0, kind(0.d0))

  Call Control_Check_Phase(Z)

  ! Should never reach here
  Write(*,*) 'ERROR: a NaN phase was not rejected!'
End Program Test_Control_phase_nan_fail
