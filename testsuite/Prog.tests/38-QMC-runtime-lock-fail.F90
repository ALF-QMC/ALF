! Unit test for QMC_runtime_var: a setter called after
! lock_QMC_runtime_settings must terminate the program (WILL_FAIL test).
Program Test_QMC_runtime_lock_fail

  Use QMC_runtime_var, only: set_QMC_runtime_default_var, lock_QMC_runtime_settings, set_Nwrap
#ifdef MPI
  Use mpi
#endif
  Implicit none
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
#endif

  Call set_QMC_runtime_default_var()
  Call set_Nwrap(10)
  Call lock_QMC_runtime_settings()
  Call set_Nwrap(20)   ! must terminate

  write(*,*) "ERROR: setter after the lock did not terminate"
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
end Program Test_QMC_runtime_lock_fail
