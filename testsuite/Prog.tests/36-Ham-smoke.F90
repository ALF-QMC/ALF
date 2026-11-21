! Generic smoke test for the Hamiltonian interface.
! Usage: 36-Ham-smoke <Hamiltonian name>   (one ctest per entry of
! Prog/Hamiltonians.list, see CMakeLists.txt).
!
! Writes a parameters file with empty namelist groups, so the model runs
! with its built-in defaults (except for required nonzero Spin_Peierls
! phonon parameters), then calls Alloc_Ham, ham%Ham_set and
! Validate_Ham_Variables exactly as main.F90 does. Validate_Ham_Variables
! already checks Ndim, N_FL, N_SUN, Ltrot, Thtrot, the shapes of Op_V and
! Op_T and the operator site indices; this program adds a few generic
! consistency checks. Model-specific values are deliberately not tested:
! the bit-identical comparison against master covers them.
! This program should exit successfully (return code 0).
!
Program Test_Ham_smoke
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main
  Implicit None

  Character (len=64) :: name
  Integer :: unit_para, n, nf
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Group_Comm = MPI_COMM_WORLD
#endif

  Call get_command_argument(1, name)
  Call require(len_trim(name) > 0, 'usage: 36-Ham-smoke <Hamiltonian name>')

  Open(newunit=unit_para, file='parameters', status='replace', action='write')
  Write(unit_para,'(A)') '&VAR_lattice /'
  Write(unit_para,'(A)') '&VAR_Model_Generic /'
  if (trim(name) == 'Spin_Peierls') then
     Write(unit_para,'(A)') '&VAR_Spin_Peierls'
     Write(unit_para,'(A)') ' Ham_Lambda=1.0, Ham_Omega0=1.0'
     Write(unit_para,'(A)') '/'
  else
     Write(unit_para,'(A)') '&VAR_' // trim(name) // ' /'
  endif
  Close(unit_para)

  Call Alloc_Ham(name)
  Call ham%Ham_set()
  Call Validate_Ham_Variables()
  Call Setup_flavor_map()

  Call require(N_FL_eff >= 1 .and. N_FL_eff <= N_FL, 'N_FL_eff out of range')
  Call require(size(Calc_Fl_map) == N_FL_eff, 'Calc_Fl_map has wrong size')
  Call require(count(Calc_Fl) == N_FL_eff, 'Calc_Fl and N_FL_eff disagree')
  Call require(reconstruction_needed .eqv. (N_FL_eff /= N_FL), 'reconstruction_needed inconsistent')
  do nf = 1, N_FL_eff
     Call require(Calc_Fl(Calc_Fl_map(nf)), 'Calc_Fl_map points to a flavour that is not computed')
  enddo

  ! Every operator must act on at least one site; every interaction
  ! operator must have a field type (Op_make leaves type = 0)
  Call require(size(Op_V, 1) >= 1, 'no interaction operators')
  do nf = 1, N_FL
     do n = 1, size(Op_T, 1)
        Call require(Op_T(n,nf)%N >= 1, 'Op_T operator without sites')
     enddo
     do n = 1, size(Op_V, 1)
        Call require(Op_V(n,nf)%N >= 1, 'Op_V operator without sites')
        Call require(Op_V(n,nf)%type >= 1 .and. Op_V(n,nf)%type <= 4, &
             &       'Op_V operator without valid field type (1..4)')
     enddo
  enddo

  Write(*,'(A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'PASS: ', trim(name), &
       & '  Ndim=', Ndim, ' N_FL=', N_FL, ' N_SUN=', N_SUN, ' Ltrot=', Ltrot, &
       & ' Op_T=[', size(Op_T,1), '] Op_V=[', size(Op_V,1), ']'
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif

contains

  Subroutine require(condition, message)
    Logical, intent(in) :: condition
    Character (len=*), intent(in) :: message
    if (.not. condition) then
       Write(*,*) 'FAIL: ', trim(name), ': ', message
       error stop 1
    endif
  end Subroutine require

end Program Test_Ham_smoke
