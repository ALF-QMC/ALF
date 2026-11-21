! Unit test for QMC_runtime_var: default values, set/get round trip and
! the settings lock (getters keep working after the lock).
Program Test_QMC_runtime_var

  Use QMC_runtime_var, only: set_QMC_runtime_default_var, lock_QMC_runtime_settings, &
       &                     get_Nwrap, set_Nwrap, get_NSweep, set_NSweep, get_NBin, set_NBin, &
       &                     get_Ltau, set_Ltau, get_CPU_MAX, set_CPU_MAX, get_sequential, set_sequential, &
       &                     get_Langevin, set_Langevin, get_HMC, set_HMC, get_N_HMC_sweeps, &
       &                     get_Amplitude, set_Amplitude, get_Nt_sequential_start, get_Nt_sequential_end, &
       &                     get_N_Global_tau, get_N_Global_tau_MALA, get_Delta_t_MALA_global, &
       &                     get_Max_Force_MALA_global, get_Global_moves, set_Global_moves, &
       &                     get_ham_name, set_ham_name, apply_CPU_MAX_policy, &
       &                     apply_Langevin_overrides, get_Global_tau_moves, set_Global_tau_moves, &
       &                     get_Global_MALA_moves, set_Global_MALA_moves, get_Global_tau_MALA_moves, &
       &                     set_Global_tau_MALA_moves
#ifdef MPI
  Use mpi
#endif
  Implicit none

  Integer :: n_fail
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
#endif
  n_fail = 0

  ! Defaults
  Call set_QMC_runtime_default_var()
  Call check_i("default Nwrap", get_Nwrap(), 0)
  Call check_i("default NSweep", get_NSweep(), 0)
  Call check_i("default NBin", get_NBin(), 0)
  Call check_i("default Ltau", get_Ltau(), 0)
  Call check_i("default N_HMC_sweeps", get_N_HMC_sweeps(), 1)
  Call check_i("default Nt_sequential_start", get_Nt_sequential_start(), 1)
  Call check_i("default Nt_sequential_end", get_Nt_sequential_end(), 0)
  Call check_i("default N_Global_tau", get_N_Global_tau(), 0)
  Call check_i("default N_Global_tau_MALA", get_N_Global_tau_MALA(), 1)
  Call check_r("default CPU_MAX", get_CPU_MAX(), 0.d0)
  Call check_r("default Amplitude", get_Amplitude(), 1.d0)
  Call check_r("default Delta_t_MALA_global", get_Delta_t_MALA_global(), 0.1d0)
  Call check_r("default Max_Force_MALA_global", get_Max_Force_MALA_global(), 1.d0)
  Call check_l("default sequential", get_sequential(), .true.)
  Call check_l("default Langevin", get_Langevin(), .false.)
  Call check_l("default HMC", get_HMC(), .false.)
  Call check_l("default Global_moves", get_Global_moves(), .false.)

  ! Set/get round trip
  Call set_Nwrap(10);          Call check_i("Nwrap", get_Nwrap(), 10)
  Call set_NSweep(20);         Call check_i("NSweep", get_NSweep(), 20)
  Call set_NBin(5);            Call check_i("NBin", get_NBin(), 5)
  Call set_Ltau(1);            Call check_i("Ltau", get_Ltau(), 1)
  Call set_CPU_MAX(2.5d0);     Call check_r("CPU_MAX", get_CPU_MAX(), 2.5d0)
  Call set_Amplitude(0.25d0);  Call check_r("Amplitude", get_Amplitude(), 0.25d0)
  Call set_sequential(.false.); Call check_l("sequential", get_sequential(), .false.)
  Call set_Langevin(.true.);   Call check_l("Langevin", get_Langevin(), .true.)
  Call set_HMC(.true.);        Call check_l("HMC", get_HMC(), .true.)
  Call set_Global_moves(.true.); Call check_l("Global_moves", get_Global_moves(), .true.)
  Call set_ham_name("Hubbard")
  If (trim(get_ham_name()) /= "Hubbard") Call fail("ham_name")

  ! CPU_MAX policy: without a CPU-time limit NBin is kept, with one it is large
  Call set_NBin(5); Call set_CPU_MAX(0.d0); Call apply_CPU_MAX_policy()
  Call check_i("NBin without CPU_MAX", get_NBin(), 5)
  Call set_CPU_MAX(2.5d0); Call apply_CPU_MAX_policy()
  Call check_i("NBin with CPU_MAX", get_NBin(), 10000000)

  ! Langevin overrides switch off all other update schemes
  Call set_sequential(.true.); Call set_HMC(.true.); Call set_Global_moves(.true.)
  Call set_Global_tau_moves(.true.); Call set_Global_MALA_moves(.true.)
  Call set_Global_tau_MALA_moves(.true.)
  Call apply_Langevin_overrides()
  Call check_l("sequential after Langevin overrides", get_sequential(), .false.)
  Call check_l("HMC after Langevin overrides", get_HMC(), .false.)
  Call check_l("Global_moves after Langevin overrides", get_Global_moves(), .false.)
  Call check_l("Global_tau_moves after Langevin overrides", get_Global_tau_moves(), .false.)
  Call check_l("Global_MALA_moves after Langevin overrides", get_Global_MALA_moves(), .false.)
  Call check_l("Global_tau_MALA_moves after Langevin overrides", get_Global_tau_MALA_moves(), .false.)
  Call set_HMC(.true.)   ! restore for the lock check below

  ! After the lock the getters still return the values
  Call lock_QMC_runtime_settings()
  Call check_i("Nwrap after lock", get_Nwrap(), 10)
  Call check_l("HMC after lock", get_HMC(), .true.)

#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  If (n_fail > 0) then
     write(*,*) n_fail, " check(s) failed"
     error stop 1
  endif
  write(*,*) "success"

contains

  Subroutine fail(what)
    Character(len=*), intent(in) :: what
    write(*,*) "FAILED: ", what
    n_fail = n_fail + 1
  end Subroutine fail

  Subroutine check_i(what, got, expected)
    Character(len=*), intent(in) :: what
    Integer, intent(in) :: got, expected
    If (got /= expected) Call fail(what)
  end Subroutine check_i

  Subroutine check_r(what, got, expected)
    Character(len=*), intent(in) :: what
    Real(Kind=Kind(0.d0)), intent(in) :: got, expected
    If (got /= expected) Call fail(what)
  end Subroutine check_r

  Subroutine check_l(what, got, expected)
    Character(len=*), intent(in) :: what
    Logical, intent(in) :: got, expected
    If (got .neqv. expected) Call fail(what)
  end Subroutine check_l

end Program Test_QMC_runtime_var
