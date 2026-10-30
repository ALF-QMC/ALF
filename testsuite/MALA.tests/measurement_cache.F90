! A one-orbital model with Gaussian bosonic action and a known Green function.
! Test callbacks reject all observations inside a trial update. After the decision,
! they check the retained fields, Green function and phase, including on rejection.
module measurement_fixture
  use Hamiltonian_main
  use Fields_mod
  implicit none
  integer, parameter :: dp = kind(0.d0)
  logical :: measurements_allowed = .false.
  integer :: measurements = 0
  integer :: tau_measurements = 0
  complex(dp) :: retained(1,2)
  type, extends(ham_base) :: fixture_hamiltonian
  contains
    procedure, nopass :: Ham_Langevin_HMC_S0 => bosonic_force
    procedure, nopass :: Get_Delta_S0_global => delta_action
    procedure, nopass :: Global_MALA_move => select_fields
    procedure, nopass :: Obser => observe
    procedure, nopass :: ObserT => observe_tau
  end type
contains
  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      write(*,*) 'FAIL: ', message
      stop 1
    endif
  end subroutine

  subroutine bosonic_force(Forces_0)
    real(dp), allocatable, intent(inout) :: Forces_0(:,:)
    Forces_0 = real(nsigma%f,dp)
  end subroutine

  real(dp) function delta_action(nsigma_old)
    type(Fields), intent(in) :: nsigma_old
    delta_action = 0.5d0 * sum(real(nsigma_old%f,dp)**2 - real(nsigma%f,dp)**2)
  end function

  subroutine select_fields(Flip_list)
    integer, intent(out) :: Flip_list(:,:)
    Flip_list = 1
  end subroutine

  subroutine observe(GR, Phase, Ntau, Mc_step_weight)
    complex(dp), intent(in) :: GR(Ndim,Ndim,N_FL), Phase
    integer, intent(in) :: Ntau
    real(dp), intent(in) :: Mc_step_weight
    real(dp) :: expected_green
    call require(measurements_allowed, 'observable called before acceptance/rejection')
    call require(maxval(abs(nsigma%f-retained)) < 1.d-12, 'measuring trial fields')
    expected_green = 1.d0 / (1.d0 + exp(real(Op_V(1,1)%g * sum(retained),dp)))
    call require(abs(GR(1,1,1)-expected_green) < 1.d-10, 'Green function is not for retained fields')
    call require(abs(Phase-1.d0) < 1.d-10, 'incorrect retained phase')
    call require(abs(Mc_step_weight-1.d0) < 1.d-12, 'Metropolis measurement must have unit weight')
    call require(Ntau >= 1 .and. Ntau <= Ltrot, 'unexpected measurement time')
    measurements = measurements + 1
  end subroutine

  subroutine observe_tau(NT, GT0, G0T, G00, GTT, PHASE, Mc_step_weight)
    integer, intent(in) :: NT
    complex(dp), intent(in) :: GT0(Ndim,Ndim,N_FL), G0T(Ndim,Ndim,N_FL)
    complex(dp), intent(in) :: G00(Ndim,Ndim,N_FL), GTT(Ndim,Ndim,N_FL), PHASE
    real(dp), intent(in) :: Mc_step_weight
    real(dp) :: expected_green, propagation
    call require(measurements_allowed, 'time-displaced observable called during trial update')
    call require(maxval(abs(nsigma%f-retained)) < 1.d-12, 'time-displaced observation of trial fields')
    call require(NT >= 0 .and. NT <= Ltrot, 'unexpected time displacement')
    expected_green = 1.d0 / (1.d0 + exp(real(Op_V(1,1)%g * sum(retained),dp)))
    propagation = exp(real(Op_V(1,1)%g * sum(retained(1,1:NT)),dp))
    call require(abs(G00(1,1,1)-expected_green) < 1.d-10, 'incorrect retained G00')
    call require(abs(GTT(1,1,1)-expected_green) < 1.d-10, 'incorrect retained GTT')
    call require(abs(GT0(1,1,1)-expected_green*propagation) < 1.d-10, 'incorrect retained GT0')
    call require(abs(G0T(1,1,1)+(1.d0-expected_green)/propagation) < 1.d-10, 'incorrect retained G0T')
    call require(abs(Phase-1.d0) < 1.d-10, 'incorrect time-displaced phase')
    call require(abs(Mc_step_weight-1.d0) < 1.d-12, 'incorrect time-displaced weight')
    tau_measurements = tau_measurements + 1
  end subroutine
end module

program test_measurement_cache
  use measurement_fixture
  use Langevin_HMC_mod
  use Hop_mod
  use Control
  use Random_wrap
  use UDV_State_mod
  use Tau_m_mod
  implicit none
  type(Langevin_HMC_type) :: cached, fresh
  class(UDV_State), allocatable :: left(:), right(:), storage(:,:)
  complex(dp), allocatable :: gr(:,:,:), gr_tilde(:,:,:), work(:,:)
  complex(dp) :: phase, cached_phase, before(1,2), after(1,2)
  integer, allocatable :: stab(:), seed_before(:), seed_after(:), seed_fresh(:)
  integer :: i, scheme, k, accepted, rejected, count_before
  integer :: seed(1) = [8752143]

  Ndim=1; N_FL=1; N_FL_eff=1; N_SUN=1; Ltrot=2; Thtrot=0; Group_Comm=0
  Projector=.false.; Symm=.false.; reconstruction_needed=.false.; leap_frog_bulk=.false.
  allocate(fixture_hamiltonian :: ham)
  allocate(Calc_Fl(1),Calc_Fl_map(1),Op_T(1,1),Op_V(1,1))
  Calc_Fl=.true.; Calc_Fl_map=1
  call Op_make(Op_T(1,1),1)
  Op_T(1,1)%P=1; Op_T(1,1)%O=0.d0; Op_T(1,1)%g=1.d0; Op_T(1,1)%type=0
  call Op_set(Op_T(1,1))
  call Op_make(Op_V(1,1),1)
  Op_V(1,1)%P=1; Op_V(1,1)%O=1.d0; Op_V(1,1)%g=0.7d0
  Op_V(1,1)%alpha=-0.5d0; Op_V(1,1)%type=3
  call Op_set(Op_V(1,1))
  call nsigma%make(1,Ltrot)
  nsigma%t=3
  call Hop_mod_init
  call Control_init(Group_Comm)
  allocate(left(1),right(1),storage(2,1),stab(0:2),gr(1,1,1),gr_tilde(1,1,1),work(1,1))
  stab=[0,1,2]
  call left(1)%init(Ndim,'l')
  call right(1)%init(Ndim,'r')
  do i=1,2
    call storage(i,1)%alloc(Ndim)
  enddo
  call Get_seed_Len(k)
  allocate(seed_before(k),seed_after(k),seed_fresh(k))
  ! These tests use local samplers. Tau_m must only propagate and measure,
  ! without using or modifying the unrelated module-level Langevin/HMC object.
  call Langevin_HMC%set_Update_scheme(.false.,.false.,.false.)

  do scheme=1,2
    ! MALA uses its cache across decisions; HMC is recomputed here because its
    ! production caller refreshes the force after each completed trajectory.
    call cached%make(.false.,scheme==2,scheme==1,0.8d0,5.d0,2)
    call fresh%make(.false.,scheme==2,scheme==1,0.8d0,5.d0,2)
    call Ranset(seed)
    nsigma%f=0.2d0
    accepted=0; rejected=0
    call reset_storage()
    do i=1,200
      before=nsigma%f
      call Ranget(seed_before)
      measurements_allowed=.false.
      count_before=measurements
      call cached%update(phase,gr,gr_tilde,work,right,left,stab,storage,1,Ltrot,1)
      call require(measurements==count_before, 'update must not measure')
      after=nsigma%f; cached_phase=phase
      call Ranget(seed_after)
      if (maxval(abs(after-before)) > 1.d-12) then
        accepted=accepted+1
      else
        rejected=rejected+1
      endif

      ! Replay the identical proposal, but explicitly recompute the old forces
      ! and determinants. No observation pass is allowed to repair cached first.
      nsigma%f=before
      call Ranset(seed_before)
      call reset_storage()
      call fresh%set_L_Forces(.false.)
      call fresh%update(phase,gr,gr_tilde,work,right,left,stab,storage,1,Ltrot,1)
      call require(maxval(abs(nsigma%f-after)) < 1.d-10, 'cached and fresh proposals differ')
      call require(abs(phase-cached_phase) < 1.d-10, 'cached and fresh phases differ')
      call Ranget(seed_fresh)
      call require(all(seed_fresh==seed_after), 'cached and fresh RNG states differ')

      ! Measure only the retained state. Use fresh, leaving the MALA cache under
      ! test untouched, so a measurement pass cannot hide rejection-cache bugs.
      retained=after
      measurements_allowed=.true.
      call fresh%calc_Forces(phase,gr,gr_tilde,work,right,left,stab,storage,1,Ltrot,.true.)
      call require(measurements-count_before==Ltrot, 'must measure after rejection as well as acceptance')
      call reset_storage()
      count_before=tau_measurements
      call Tau_m(storage,gr,phase,2,1,stab,1,Ltrot)
      call require(tau_measurements-count_before==Ltrot+1, 'missing time-displaced observations')
      measurements_allowed=.false.
      call cached%set_L_Forces(scheme==1)
    enddo
    call require(accepted>0 .and. rejected>0, 'fixture must exercise acceptance and rejection')
    write(*,*) 'Scheme ',scheme,': accepted ',accepted,', rejected ',rejected
    call cached%clean()
    call fresh%clean()
  enddo
  write(*,*) 'Measurement boundary and MALA cache checks passed'
contains
  subroutine reset_storage()
    call Langevin_HMC_Reset_storage(phase,gr,right,left,stab,storage)
  end subroutine
end program
