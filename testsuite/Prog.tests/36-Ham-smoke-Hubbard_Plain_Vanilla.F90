! Smoke test for the Hamiltonian interface: Hubbard_Plain_Vanilla.
! Writes a small parameters file, calls Alloc_Ham, ham%Ham_set and
! Validate_Ham_Variables exactly as main.F90 does, and checks the
! resulting model variables against hand-computed values.
! 4x4 square lattice, Beta = 1, Dtau = 0.1, U = 4, half filling.
! This program should exit successfully (return code 0).
!
Program Test_Ham_smoke_Hubbard_Plain_Vanilla
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main
  Implicit None

  Integer, parameter :: dp = kind(0.d0)
  Integer, parameter :: L = 4, N_sites = L*L
  Real (dp), parameter :: Dtau_ref = 0.1d0, U_ref = 4.d0
  Character (len=64) :: name
  Integer :: unit_para, i, nf
  Real (dp) :: sign_ref
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Group_Comm = MPI_COMM_WORLD
#endif

  Open(newunit=unit_para, file='parameters', status='replace', action='write')
  Write(unit_para,'(A)') '&VAR_lattice'
  Write(unit_para,'(A)') 'Lattice_type = "Square"'
  Write(unit_para,'(A,I0)') 'L1 = ', L
  Write(unit_para,'(A,I0)') 'L2 = ', L
  Write(unit_para,'(A)') '/'
  Write(unit_para,'(A)') '&VAR_Hubbard_Plain_Vanilla'
  Write(unit_para,'(A)') 'ham_T = 1.d0'
  Write(unit_para,'(A)') 'Ham_chem = 0.d0'
  Write(unit_para,'(A)') 'Ham_U = 4.d0'
  Write(unit_para,'(A)') 'Dtau = 0.1d0'
  Write(unit_para,'(A)') 'Beta = 1.d0'
  Write(unit_para,'(A)') 'Projector = .false.'
  Write(unit_para,'(A)') '/'
  Close(unit_para)

  name = 'Hubbard_Plain_Vanilla'
  Call Alloc_Ham(name)
  Call ham%Ham_set()
  Call Validate_Ham_Variables()

  Call require(Ndim == N_sites, 'Ndim /= L*L')
  Call require(N_FL == 2, 'N_FL /= 2')
  Call require(N_SUN == 1, 'N_SUN /= 1')
  Call require(Ltrot == 10, 'Ltrot /= Beta/Dtau')
  Call require(Thtrot == 0, 'Thtrot /= 0 for finite temperature')
  Call require(.not. Projector, 'Projector should be .false.')

  ! Hopping: one dense operator per flavour acting on all sites
  Call require(all(shape(Op_T) == [1, 2]), 'shape(Op_T) /= [1, N_FL]')
  do nf = 1, N_FL
     Call require(Op_T(1,nf)%N == N_sites, 'Op_T%N /= Ndim')
     Call require(abs(Op_T(1,nf)%g + Dtau_ref) < 1.d-12, 'Op_T%g /= -Dtau')
  enddo

  ! Interaction: one 1x1 operator per site and flavour, g = +-sqrt(Dtau*U/2)
  Call require(all(shape(Op_V) == [N_sites, 2]), 'shape(Op_V) /= [Ndim, N_FL]')
  do nf = 1, N_FL
     sign_ref = 1.d0
     if (nf == 2) sign_ref = -1.d0
     do i = 1, N_sites
        Call require(Op_V(i,nf)%N == 1, 'Op_V%N /= 1')
        Call require(Op_V(i,nf)%P(1) == i, 'Op_V%P(1) /= site index')
        Call require(Op_V(i,nf)%type == 2, 'Op_V%type /= 2')
        Call require(abs(Op_V(i,nf)%g - sign_ref*sqrt(Dtau_ref*U_ref/2.d0)) < 1.d-12, &
             &       'Op_V%g /= +-sqrt(Dtau*U/2)')
     enddo
  enddo

  ! Half filling: flavour 1 is reconstructed from flavour 2
  Call require(allocated(Calc_Fl), 'Calc_Fl not allocated at half filling')
  Call require((.not. Calc_Fl(1)) .and. Calc_Fl(2), 'Calc_Fl /= [F, T]')

  Write(*,*) 'PASS: Hubbard_Plain_Vanilla Ham_set smoke test'
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif

contains

  Subroutine require(condition, message)
    Logical, intent(in) :: condition
    Character (len=*), intent(in) :: message
    if (.not. condition) then
       Write(*,*) 'FAIL: ', message
       error stop 1
    endif
  end Subroutine require

end Program Test_Ham_smoke_Hubbard_Plain_Vanilla
