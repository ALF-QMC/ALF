! Smoke test for the Hamiltonian interface: Kondo lattice model
! (conduction c-layer plus localized f-layer, bilayer square lattice).
! Writes a small parameters file, calls Alloc_Ham, ham%Ham_set and
! Validate_Ham_Variables exactly as main.F90 does, and checks the
! resulting model variables against hand-computed values.
! 4x4 bilayer square, Beta = 1, Dtau = 0.1, N_SUN = 2, Uc = 0, Uf = 2,
! J_K = 2, no checkerboard decomposition.
! This program should exit successfully (return code 0).
!
! Expected Op_V (N_SUN = 2 branch of Ham_V):
!   1..16  : Hubbard Uf on the f-site of unit cell I, 1x1,
!            g = sqrt(-Dtau*Uf/N_SUN)  (purely imaginary)
!   17..32 : Kondo term on (c-site, f-site) of unit cell I, 2x2,
!            O = [[0,1],[1,0]], g = sqrt(Dtau*(J_K/2)/N_SUN)
!
Program Test_Ham_smoke_Kondo
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main
  Implicit None

  Integer, parameter :: dp = kind(0.d0)
  Integer, parameter :: L = 4, N_cells = L*L, N_sites = 2*N_cells
  Integer, parameter :: N_SUN_ref = 2
  Real (dp), parameter :: Dtau_ref = 0.1d0, Uf_ref = 2.d0, JK_ref = 2.d0
  Character (len=64) :: name
  Integer :: unit_para, n, a, b
  Logical :: is_f(N_sites), is_c(N_sites)
  Complex (dp) :: g_Uf, g_JK
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Group_Comm = MPI_COMM_WORLD
#endif

  Open(newunit=unit_para, file='parameters', status='replace', action='write')
  Write(unit_para,'(A)') '&VAR_lattice'
  Write(unit_para,'(A)') 'Lattice_type = "Bilayer_square"'
  Write(unit_para,'(A,I0)') 'L1 = ', L
  Write(unit_para,'(A,I0)') 'L2 = ', L
  Write(unit_para,'(A)') '/'
  Write(unit_para,'(A)') '&VAR_Model_Generic'
  Write(unit_para,'(A,I0)') 'N_SUN = ', N_SUN_ref
  Write(unit_para,'(A)') 'N_FL = 1'
  Write(unit_para,'(A)') 'Dtau = 0.1d0'
  Write(unit_para,'(A)') 'Beta = 1.d0'
  Write(unit_para,'(A)') 'Checkerboard = .false.'
  Write(unit_para,'(A)') 'Projector = .false.'
  Write(unit_para,'(A)') '/'
  Write(unit_para,'(A)') '&VAR_Kondo'
  Write(unit_para,'(A)') 'ham_T = 1.d0'
  Write(unit_para,'(A)') 'Ham_chem = 0.d0'
  Write(unit_para,'(A)') 'Ham_Uc = 0.d0'
  Write(unit_para,'(A)') 'ham_Uf = 2.d0'
  Write(unit_para,'(A)') 'Ham_JK = 2.d0'
  Write(unit_para,'(A)') '/'
  Close(unit_para)

  name = 'Kondo'
  Call Alloc_Ham(name)
  Call ham%Ham_set()
  Call Validate_Ham_Variables()

  Call require(Ndim == N_sites, 'Ndim /= 2*L*L (c- plus f-layer)')
  Call require(N_FL == 1, 'N_FL /= 1')
  Call require(N_SUN == N_SUN_ref, 'N_SUN /= input value')
  Call require(Ltrot == 10, 'Ltrot /= Beta/Dtau')
  Call require(Thtrot == 0, 'Thtrot /= 0 for finite temperature')
  Call require(.not. Projector, 'Projector should be .false.')

  ! Interaction: Uf block followed by J_K block
  g_Uf = sqrt(cmplx(-Dtau_ref*Uf_ref/N_SUN_ref, 0.d0, dp))
  g_JK = sqrt(cmplx( Dtau_ref*(JK_ref/2.d0)/N_SUN_ref, 0.d0, dp))
  Call require(all(shape(Op_V) == [2*N_cells, 1]), 'shape(Op_V) /= [2*L*L, 1]')
  is_f = .false.
  is_c = .false.
  do n = 1, N_cells
     Call require(Op_V(n,1)%N == 1, 'Uf operator is not 1x1')
     Call require(Op_V(n,1)%type == 2, 'Uf operator type /= 2')
     Call require(abs(Op_V(n,1)%g - g_Uf) < 1.d-12, 'Uf coupling /= sqrt(-Dtau*Uf/N_SUN)')
     Call require(.not. is_f(Op_V(n,1)%P(1)), 'Uf acts twice on the same f-site')
     is_f(Op_V(n,1)%P(1)) = .true.
  enddo
  do n = 1, N_cells
     associate (Op => Op_V(N_cells + n, 1))
       Call require(Op%N == 2, 'Kondo operator is not 2x2')
       Call require(Op%type == 2, 'Kondo operator type /= 2')
       Call require(abs(Op%g - g_JK) < 1.d-12, 'Kondo coupling /= sqrt(Dtau*(J_K/2)/N_SUN)')
       Call require(abs(Op%O(1,2) - 1.d0) < 1.d-12 .and. abs(Op%O(2,1) - 1.d0) < 1.d-12, &
            &       'Kondo operator off-diagonal /= 1')
       ! Same unit cell: its f-site is the one carrying Uf in the same position
       Call require(Op%P(2) == Op_V(n,1)%P(1), 'Kondo term not on the f-site of its unit cell')
       Call require(.not. is_f(Op%P(1)) .and. .not. is_c(Op%P(1)), 'Kondo c-site invalid or reused')
       is_c(Op%P(1)) = .true.
     end associate
  enddo
  Call require(all(is_c .neqv. is_f), 'sites are not split into 16 c- and 16 f-sites')

  ! Hopping without checkerboard: one dense operator; f-electrons do not hop
  Call require(all(shape(Op_T) == [1, 1]), 'shape(Op_T) /= [1, 1]')
  Call require(Op_T(1,1)%N == N_sites, 'Op_T%N /= Ndim')
  Call require(abs(Op_T(1,1)%g + Dtau_ref) < 1.d-12, 'Op_T%g /= -Dtau')
  do a = 1, Op_T(1,1)%N
     do b = 1, Op_T(1,1)%N
        if (is_f(Op_T(1,1)%P(a)) .or. is_f(Op_T(1,1)%P(b))) &
             & Call require(abs(Op_T(1,1)%O(a,b)) < 1.d-12, 'hopping matrix touches an f-site')
     enddo
  enddo
  Call require(any(abs(Op_T(1,1)%O) > 1.d-12), 'no c-electron hopping at all')

  ! The Kondo model does not use flavour reconstruction
  Call require(.not. allocated(Calc_Fl), 'Calc_Fl should not be allocated')

  Write(*,*) 'PASS: Kondo Ham_set smoke test'
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

end Program Test_Ham_smoke_Kondo
