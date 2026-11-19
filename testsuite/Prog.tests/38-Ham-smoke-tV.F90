! Smoke test for the Hamiltonian interface: tV (spinless fermions with
! nearest-neighbour interaction V, decoupled on bonds).
! Writes a small parameters file, calls Alloc_Ham, ham%Ham_set and
! Validate_Ham_Variables exactly as main.F90 does, and checks the
! resulting model variables against hand-computed values.
! 4x4 square lattice, Beta = 1, Dtau = 0.1, V = 2, N_SUN = 2, Symm = .true.,
! no checkerboard decomposition for the hopping.
! This program should exit successfully (return code 0).
!
! Op_V: the 2*L*L nearest-neighbour bonds come in 4 families. Symm = .true.
! orders them palindromically (1 2 3 4 3 2 1), so each bond of the outer
! families appears twice with weight 1/2. Instead of fixing the operator
! order, the test checks the invariant: every nearest-neighbour bond carries
! total weight sum(g**2) = Dtau*V/N_SUN, and every site has 4 bonds.
!
Program Test_Ham_smoke_tV
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main
  Implicit None

  Integer, parameter :: dp = kind(0.d0)
  Integer, parameter :: L = 4, N_sites = L*L, N_bonds = 2*N_sites
  Integer, parameter :: N_SUN_ref = 2
  Real (dp), parameter :: Dtau_ref = 0.1d0, V_ref = 2.d0
  Real (dp), parameter :: w_bond = Dtau_ref*V_ref/N_SUN_ref
  Character (len=64) :: name
  Integer :: unit_para, n, i, j
  Real (dp) :: W(N_sites, N_sites), g2
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
  Write(unit_para,'(A)') '&VAR_Model_Generic'
  Write(unit_para,'(A,I0)') 'N_SUN = ', N_SUN_ref
  Write(unit_para,'(A)') 'N_FL = 1'
  Write(unit_para,'(A)') 'Dtau = 0.1d0'
  Write(unit_para,'(A)') 'Beta = 1.d0'
  Write(unit_para,'(A)') 'Checkerboard = .false.'
  Write(unit_para,'(A)') 'Symm = .true.'
  Write(unit_para,'(A)') 'Projector = .false.'
  Write(unit_para,'(A)') '/'
  Write(unit_para,'(A)') '&VAR_tV'
  Write(unit_para,'(A)') 'ham_T = 1.d0'
  Write(unit_para,'(A)') 'Ham_chem = 0.d0'
  Write(unit_para,'(A)') 'ham_V = 2.d0'
  Write(unit_para,'(A)') '/'
  Close(unit_para)

  name = 'tV'
  Call Alloc_Ham(name)
  Call ham%Ham_set()
  Call Validate_Ham_Variables()

  Call require(Ndim == N_sites, 'Ndim /= L*L')
  Call require(N_FL == 1, 'N_FL /= 1')
  Call require(N_SUN == N_SUN_ref, 'N_SUN /= input value')
  Call require(Ltrot == 10, 'Ltrot /= Beta/Dtau')
  Call require(Thtrot == 0, 'Thtrot /= 0 for finite temperature')
  Call require(.not. Projector, 'Projector should be .false.')
  Call require(Symm, 'Symm should be .true.')

  ! Hopping without checkerboard: one dense operator on all sites
  Call require(all(shape(Op_T) == [1, 1]), 'shape(Op_T) /= [1, 1]')
  Call require(Op_T(1,1)%N == N_sites, 'Op_T%N /= Ndim')
  Call require(abs(Op_T(1,1)%g + Dtau_ref) < 1.d-12, 'Op_T%g /= -Dtau')

  ! Interaction: 4 families of L*L/2 bonds, symmetrized to 7 families
  Call require(all(shape(Op_V) == [7*N_sites/2, 1]), 'shape(Op_V) /= [56, 1]')
  W = 0.d0
  do n = 1, size(Op_V, 1)
     Call require(Op_V(n,1)%N == 2, 'Op_V%N /= 2 (bond operator)')
     Call require(Op_V(n,1)%type == 2, 'Op_V%type /= 2')
     i = Op_V(n,1)%P(1)
     j = Op_V(n,1)%P(2)
     Call require(i /= j, 'Op_V acts twice on the same site')
     Call require(abs(Op_V(n,1)%O(1,2) - 1.d0) < 1.d-12 .and. &
          &       abs(Op_V(n,1)%O(2,1) - 1.d0) < 1.d-12, 'Op_V%O off-diagonal /= 1')
     Call require(abs(aimag(Op_V(n,1)%g)) < 1.d-12, 'Op_V%g not real')
     g2 = real(Op_V(n,1)%g, dp)**2
     W(min(i,j), max(i,j)) = W(min(i,j), max(i,j)) + g2
  enddo
  Call require(count(W > 1.d-12) == N_bonds, 'number of distinct bonds /= 2*L*L')
  Call require(all(abs(pack(W, W > 1.d-12) - w_bond) < 1.d-12), &
       &       'bond weight sum(g**2) /= Dtau*V/N_SUN')
  do i = 1, N_sites
     Call require(count(W(i,:) > 1.d-12) + count(W(:,i) > 1.d-12) == 4, &
          &       'site does not have 4 nearest-neighbour bonds')
  enddo

  ! The tV model does not use flavour reconstruction
  Call require(.not. allocated(Calc_Fl), 'Calc_Fl should not be allocated')

  Write(*,*) 'PASS: tV Ham_set smoke test'
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

end Program Test_Ham_smoke_tV
