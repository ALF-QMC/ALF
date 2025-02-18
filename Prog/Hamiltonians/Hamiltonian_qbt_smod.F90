    submodule(Hamiltonian_main) ham_qbt_smod

       use Operator_mod
       use WaveFunction_mod
       use Lattices_v3
       use MyMats
       use Random_Wrap
       use Files_mod
       use Matrix
       use Observables
       use Fields_mod
       use Predefined_Hoppings
       use LRC_Mod

       implicit none

       type, extends(ham_base) :: ham_qbt
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: Alloc_obs
          procedure, nopass :: Obser
          procedure, nopass :: ObserT
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_qbt

       !#PARAMETERS START# VAR_lattice
       character(len=64) :: Model = ''  ! Value not relevant
       character(len=64) :: Lattice_type = 'Pi_Flux'
       integer            :: L1 = 6   ! Length in direction a_1
       integer            :: L2 = 6   ! Length in direction a_2
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_Model_Generic
       !Integer              :: N_SUN        = 1        ! Number of colors
       !Integer              :: N_FL         = 1        ! Number of flavors
       real(Kind=kind(0.d0)) :: Phi_X = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
       real(Kind=kind(0.d0)) :: Phi_Y = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
       logical               :: Bulk = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
       integer               :: N_Phi = 0        ! Total number of flux quanta traversing the lattice
       real(Kind=kind(0.d0)) :: Dtau = 0.1d0    ! Thereby Ltrot=Beta/dtau
       real(Kind=kind(0.d0)) :: Beta = 5.d0     ! Inverse temperature
       logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
       !logical              :: Symm         = .true.   ! Whether symmetrization takes place
       !logical              :: Projector    = .false.  ! Whether the projective algorithm is used
       real(Kind=kind(0.d0)) :: Theta = 10.d0    ! Projection parameter
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_qbt
       real(Kind=kind(0.d0)) :: ham_T = 1.d0  ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_T2 = 1.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_T3 = 1.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_V = 4.d0  ! Hubbard interaction
       real(Kind=kind(0.d0)) :: ham_V2 = 4.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_V3 = 4.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: Ham_chem = 0.d0  ! Chemical potential
       real(Kind=kind(0.d0)) :: ham_Tperp = 1.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_Vperp = 1.d0  ! For bilayer systems
       integer               :: N_dope = 0        ! Number of doping electrons
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable   :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell
       integer, allocatable   :: bond_list(:, :, :), l_bond(:)  ! bond list
       integer, allocatable   :: bond_map_v1(:), bond_map_v2(:)  ! bond list

       type(Unit_cell), target  :: Latt_unit_qah  ! ab sublattice

    contains

       module subroutine Ham_Alloc_qbt
          allocate (ham_qbt::ham)
       end subroutine Ham_Alloc_qbt

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_qbt_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
       subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          use mpi
#endif
          implicit none

          integer                :: ierr, nf, unit_info
          character(len=64)     :: file_info

          ! L1, L2, Lattice_type, List(:,:), Invlist(:,:) -->  Lattice information
          ! Ham_T, Chem, Phi_X, XB_B, Checkerboard, Symm   -->  Hopping
          ! Interaction                              -->  Model
          ! Simulation type                          -->  Finite  T or Projection  Symmetrize Trotter.

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif

          ! From dynamically generated file "Hamiltonian_tV_read_write_parameters.F90"
          call read_parameters()

          Ltrot = nint(beta/dtau)
          Thtrot = 0
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot + 2*Thtrot
          N_FL = 1

          ! Setup the Bravais lattice
          call Ham_Latt

          ! Setup the hopping / single-particle part
          call Ham_Hop

          ! Setup the interaction.
          call Ham_Vint

          ! Setup the trival wave function, in case of a projector approach
          if (Projector) call Ham_Trial()

#ifdef MPI
          if (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write (File_info, '(A,I0,A)') "Temp_", igroup, "/info"
#endif
             open (newunit=unit_info, file=file_info, status="unknown", position="append")
             write (unit_info, *) '====================================='
             write (unit_info, *) 'Model is      : ', Model
             write (unit_info, *) 'Lattice is    : ', Lattice_type
             write (unit_info, *) '# of orbitals : ', Ndim
             write (unit_info, *) 'Flux_1        : ', Phi_X
             write (unit_info, *) 'Flux_2        : ', Phi_Y
             if (Bulk) then
                write (unit_info, *) 'Twist as phase factor in bulk'
             else
                write (unit_info, *) 'Twist as boundary condition'
             end if
             write (unit_info, *) 'Checkerboard  : ', Checkerboard
             write (unit_info, *) 'Symm. decomp  : ', Symm
             if (Projector) then
                write (unit_info, *) 'Projective version'
                write (unit_info, *) 'Theta         : ', Theta
                write (unit_info, *) 'Tau_max       : ', beta
             else
                write (unit_info, *) 'Finite temperture version'
                write (unit_info, *) 'Beta          : ', Beta
             end if
             write (unit_info, *) 'dtau,Ltrot_eff: ', dtau, Ltrot
             write (unit_info, *) 'N_SUN         : ', N_SUN
             write (unit_info, *) 'N_FL          : ', N_FL
             write (unit_info, *) 't             : ', Ham_T
             write (unit_info, *) 't2            : ', Ham_T2
             write (unit_info, *) 't3            : ', Ham_T3
             write (unit_info, *) 'Ham_V         : ', Ham_V
             write (unit_info, *) 'Ham_V2        : ', Ham_V2
             write (unit_info, *) 'Ham_V3        : ', Ham_V3
             write (unit_info, *) 'Ham_chem      : ', Ham_chem
             write (unit_info, *) 'N_dope        : ', N_dope
             if (Projector) then
                do nf = 1, N_FL
                   write (unit_info, *) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   write (unit_info, *) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                end do
             end if
             close (unit_info)
#ifdef MPI
          end if
#endif
       end subroutine Ham_Set

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
       subroutine Ham_Latt

          use Predefined_Lattices

          implicit none
          integer :: n

          ! Use predefined stuctures or set your own lattice.
          call Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)

          Latt_Unit_qah%Norb = 4
          Latt_Unit_qah%N_coord = 2
          allocate (Latt_Unit_qah%Orb_pos_p(4, 3))
          Latt_unit_qah%Orb_pos_p = 0.d0
          do n = 1, Latt_Unit_qah%Norb
             Latt_Unit_qah%Orb_pos_p(n, 3) = real(n - 1, kind(0.d0))
          end do

       end subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the Hopping
!--------------------------------------------------------------------
       subroutine Ham_Hop

          implicit none

          real(Kind=kind(0.d0)) ::  Ham_Lambda = 0.d0

          real(Kind=kind(0.d0)), allocatable :: Ham_T_vec(:), Ham_T2_vec(:), Ham_T3_vec(:), &
              & Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
              & Ham_Lambda_vec(:)
          integer, allocatable ::   N_Phi_vec(:)

          ! Use predefined stuctures or set your own hopping
          integer :: n, nth

          allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_T3_vec(N_FL), Ham_Tperp_vec(N_FL), &
              & Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL), N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL))

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_T_vec = Ham_T
          Ham_T2_vec = Ham_T2
          Ham_T3_vec = Ham_T3
          Ham_Tperp_vec = Ham_Tperp
          Ham_Chem_vec = Ham_Chem
          Phi_X_vec = Phi_X
          Phi_Y_vec = Phi_Y
          Ham_Lambda_vec = Ham_Lambda
          N_Phi_vec = N_Phi

          select case (Lattice_type)
          case ("Pi_Flux")
             call set_hopping_parameters_pi_flux_qbt(Hopping_Matrix, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, &
                 & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          end select

          call Predefined_Hoppings_set_OPT(Hopping_Matrix, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_T)

          deallocate (Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                   N_Phi_vec, Ham_Lambda_vec)

       end subroutine Ham_Hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the trial wave function
!--------------------------------------------------------------------
       subroutine Ham_Trial()

          use Predefined_Trial

          implicit none

          integer :: N_part, nf
          character(len=64) :: model_type = 'qbt'

          ! Use predefined stuctures or set your own Trial  wave function
          N_part = Ndim/2 - N_dope
          call Predefined_TrialWaveFunction(model_type, Ndim, List, Invlist, Latt, Latt_unit, &
               &                            N_part, 1.d0, N_FL, WF_L, WF_R)

       end subroutine Ham_Trial

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
       subroutine Ham_Vint

          use Predefined_Int
          implicit none

          integer                           :: I, J, I1, J1, no_I, no_J, nf, n_b, n_ops, amx, amy, npxy
          integer                           :: dx, dy, dnb, lly, n_b_t, k1, nc0
          integer                           :: n_1, n_2, Nb, n_f, l_f, n_l, N, nc, I0, J0, n_cb, nu_bond
          complex(Kind=kind(0.d0))          :: Z
          real(Kind=kind(0.d0))             :: Zero = 1.0e-6

          amx = mod(l1, 2)
          amy = mod(l2, 2)

          nu_bond = 12 + amx + 3*amy
          allocate (bond_list(nu_bond, Latt%N, 3))
          allocate (l_bond(nu_bond))

          l_bond = 0

          ! set bond list
          do n_b = 1, 12
             nc = 0
             select case (n_b)
             case (1)
                do I = 1, Latt%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = I; 
                   bond_list(n_b, nc, 2) = Invlist(I, 1); 
                   bond_list(n_b, nc, 3) = invlist(I, 2)
                end do
                l_bond(n_b) = nc
             case (2)
                do I = 1, Latt%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = I; 
                   bond_list(n_b, nc, 2) = Invlist(I, 2); 
                   bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I, 1, 0), 1)
                end do
                l_bond(n_b) = nc
             case (3)
                do I = 1, Latt%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = I; 
                   bond_list(n_b, nc, 2) = Invlist(I, 2); 
                   bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I, 1, 1), 1)
                end do
                l_bond(n_b) = nc
             case (4)
                do I = 1, Latt%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = I; 
                   bond_list(n_b, nc, 2) = Invlist(I, 2); 
                   bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I, 0, 1), 1)
                end do
                l_bond(n_b) = nc
             case (5)
                i0 = 1
                do dx = 1, L1 - amx
                do dy = 1, L2
                   if (mod(latt%list(i0, 1), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 0), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 0), 2)
                   end if

                   i0 = latt%nnlist(I0, 0, 1)
                end do
                i0 = latt%nnlist(I0, 1, 0)
                end do
                l_bond(n_b) = nc
             case (6)
                i0 = 1
                do dx = 1, L1 - amx
                do dy = 1, L2
                   if (mod(latt%list(i0, 1), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 0), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 0), 2)
                   end if

                   i0 = latt%nnlist(I0, 0, 1)
                end do
                i0 = latt%nnlist(I0, 1, 0)
                end do
                l_bond(n_b) = nc
             case (7)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt%list(i0, 2), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 0, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 0, 1), 2)
                   end if
                   i0 = latt%nnlist(I0, 1, 0)
                end do
                i0 = latt%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (8)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt%list(i0, 2), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 0, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 0, 1), 2)
                   end if
                   i0 = latt%nnlist(I0, 1, 0)
                end do
                i0 = latt%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (9)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt%list(i0, 2), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 1), 2)
                   end if
                   i0 = latt%nnlist(I0, 1, 0)
                end do
                i0 = latt%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (10)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt%list(i0, 2), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 1), 2)
                   end if
                   i0 = latt%nnlist(I0, 1, 0)
                end do
                i0 = latt%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (11)
                i0 = 1
                i0 = latt%nnlist(I0, 0, 1)
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt%list(i0, 2), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, -1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, -1), 2)
                   end if
                   i0 = latt%nnlist(I0, 1, 0)
                end do
                i0 = latt%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (12)
                i0 = 1
                i0 = latt%nnlist(I0, 0, 1)
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt%list(i0, 2), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, -1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = I0; 
                      bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, -1), 2)
                   end if
                   i0 = latt%nnlist(I0, 1, 0)
                end do
                i0 = latt%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             end select
          end do

          n_b = 12
          do npxy = 1, amx
             n_b = n_b + 1
             i0 = 1
             do dx = 1, L1 - amx
                i0 = latt%nnlist(I0, 1, 0)
             end do
             nc = 0
             do dy = 1, L2
                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 0), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 0), 2)

                i0 = latt%nnlist(I0, 0, 1)
             end do
             l_bond(n_b) = nc
          end do

          do npxy = 1, amy
             n_b = n_b + 1
             i0 = 1
             do dy = 1, L2 - amy
                i0 = latt%nnlist(I0, 0, 1)
             end do
             nc = 0
             do dx = 1, L1
                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 0, 1), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 0, 1), 2)

                i0 = latt%nnlist(I0, 1, 0)
             end do
             l_bond(n_b) = nc

             n_b = n_b + 1
             i0 = 1
             do dy = 1, L2 - amy
                i0 = latt%nnlist(I0, 0, 1)
             end do
             nc = 0
             do dx = 1, L1
                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 1), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, 1), 2)

                i0 = latt%nnlist(I0, 1, 0)
             end do
             l_bond(n_b) = nc

             n_b = n_b + 1
             i0 = 1
             nc = 0
             do dx = 1, L1
                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 1); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, -1), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = I0; 
                bond_list(n_b, nc, 2) = Invlist(I0, 2); 
                bond_list(n_b, nc, 3) = invlist(Latt%nnlist(I0, 1, -1), 2)

                i0 = latt%nnlist(I0, 1, 0)
             end do
             l_bond(n_b) = nc

          end do

          N_ops = 0
          if (abs(Ham_V) > Zero) then
             do J1 = 1, 4
                N_ops = N_ops + 2*l_bond(J1)
             end do
          end if
          if (abs(Ham_V2) > Zero) then
             do J1 = 5, 8
                N_ops = N_ops + 2*l_bond(J1)
             end do
             if (amx .eq. 1) N_ops = N_ops + 2*l_bond(13)
             if (amy .eq. 1) N_ops = N_ops + 2*l_bond(12 + 1 + amx)
          end if
          if (abs(Ham_V3) > Zero) then
             do J1 = 9, 12
                N_ops = N_ops + 2*l_bond(J1)
             end do
             if (amy .eq. 1) N_ops = N_ops + 2*l_bond(12 + 2 + amx) + 2*l_bond(12 + 3 + amx)
          end if
          allocate (Op_V(N_ops, N_FL))

          lly = 4
          allocate (bond_map_v1(lly)) ! for V_1
          lly = 4 + amx + amy
          allocate (bond_map_v2(lly)) ! for V_2

          nc = 1
          do J1 = 1, 4
             bond_map_v1(J1) = J1
             bond_map_v2(J1) = J1 + 4
          end do

          do J1 = 1, amx
             bond_map_v2(4 + amx) = 12 + amx
          end do
          do J1 = 1, amy
             bond_map_v2(4 + amx + amy) = 12 + amx + amy
          end do

          do nf = 1, N_FL
             do nc = 1, size(Op_V, 1)
                call Op_make(Op_V(nc, nf), 2)
             end do
          end do
          nc = 0

          if (abs(Ham_V) > Zero) then
             nf = 1
             lly = size(bond_map_v1)
             do n_b_t = 1, lly
                n_cb = bond_map_v1(n_b_t)
                do J0 = 1, l_bond(n_cb)
                   I = bond_list(n_cb, J0, 2)
                   J = bond_list(n_cb, J0, 3)
                   nc = nc + 1
                   Op_V(nc, nf)%P(1) = I; Op_V(nc, nf)%P(2) = J
                   Op_V(nc, nf)%O(1, 2) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%g = sqrt(cmplx(0.5d0*0.5d0*dtau*Ham_V, 0.d0, kind(0.d0)))
                   Op_V(nc, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%type = 2
                   call Op_set(Op_V(nc, nf))
                end do
             end do
          end if

          if (abs(Ham_V2) > Zero) then
             nf = 1
             lly = size(bond_map_v2)
             do n_b_t = 1, lly
                n_cb = bond_map_v2(n_b_t)
                do J0 = 1, l_bond(n_cb)
                   I = bond_list(n_cb, J0, 2)
                   J = bond_list(n_cb, J0, 3)
                   nc = nc + 1
                   Op_V(nc, nf)%P(1) = I; Op_V(nc, nf)%P(2) = J
                   Op_V(nc, nf)%O(1, 2) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%g = sqrt(cmplx(0.5d0*0.5d0*dtau*Ham_V2, 0.d0, kind(0.d0)))
                   Op_V(nc, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%type = 2
                   call Op_set(Op_V(nc, nf))
                end do
             end do

             do n_b_t = lly, 1, -1
                n_cb = bond_map_v2(n_b_t)
                do J0 = l_bond(n_cb), 1, -1
                   I = bond_list(n_cb, J0, 2)
                   J = bond_list(n_cb, J0, 3)
                   nc = nc + 1
                   Op_V(nc, nf)%P(1) = I; Op_V(nc, nf)%P(2) = J
                   Op_V(nc, nf)%O(1, 2) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%g = sqrt(cmplx(0.5d0*0.5d0*dtau*Ham_V2, 0.d0, kind(0.d0)))
                   Op_V(nc, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%type = 2
                   call Op_set(Op_V(nc, nf))
                end do
             end do
          end if

          if (abs(Ham_V) > Zero) then
             nf = 1
             lly = size(bond_map_v1)
             do n_b_t = lly, 1, -1
                n_cb = bond_map_v1(n_b_t)
                do J0 = l_bond(n_cb), 1, -1
                   I = bond_list(n_cb, J0, 2)
                   J = bond_list(n_cb, J0, 3)
                   nc = nc + 1
                   Op_V(nc, nf)%P(1) = I; Op_V(nc, nf)%P(2) = J
                   Op_V(nc, nf)%O(1, 2) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%g = sqrt(cmplx(0.5d0*0.5d0*dtau*Ham_V, 0.d0, kind(0.d0)))
                   Op_V(nc, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%type = 2
                   call Op_set(Op_V(nc, nf))
                end do
             end do
          end if

       end subroutine Ham_Vint

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
       subroutine Alloc_obs(Ltau)

          implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          integer, intent(In) :: Ltau
          integer    ::  i, N, Nt
          character(len=64) ::  Filename
          character(len=2)  ::  Channel

          ! Scalar observables
          allocate (Obs_scal(5))
          do I = 1, size(Obs_scal, 1)
             select case (I)
             case (1)
                N = 1; Filename = "Kin"
             case (2)
                N = 1; Filename = "Pot"
             case (3)
                N = 1; Filename = "Part"
             case (4)
                N = 1; Filename = "Ener"
             case (5)
                N = 1; Filename = "QAH"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             call Obser_Vec_make(Obs_scal(I), N, Filename)
          end do

          ! Equal time correlators
          allocate (Obs_eq(4))
          do I = 1, size(Obs_eq, 1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "Den"
             case (4)
                Filename = "QAH"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             !Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             if (I == 4) then
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_qah, Channel, dtau)
             else
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             end if
          end do

          if (Ltau == 1) then
             ! Time-displaced correlators
             allocate (Obs_tau(3))
             do I = 1, size(Obs_tau, 1)
                select case (I)
                case (1)
                   Channel = 'P'; Filename = "Green"
                case (2)
                   Channel = 'PH'; Filename = "SpinZ"
                case (3)
                   Channel = 'PH'; Filename = "Den"
                case default
                   write (6, *) ' Error in Alloc_obs '
                end select
                Nt = Ltrot + 1 - 2*Thtrot
                if (Projector) Channel = 'T0'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             end do
          end if

       end subroutine Alloc_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes equal time observables
!> @details
!> @param [IN] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
       subroutine Obser(GR, Phase, Ntau, Mc_step_weight)

          use Predefined_Obs

          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: PHASE
          integer, intent(IN)          :: Ntau
          real(Kind=kind(0.d0)), intent(IN) :: Mc_step_weight

          !Local
          complex(Kind=kind(0.d0)) :: GRC(Ndim, Ndim, N_FL), ZK, Zn, weight, delta, zqah, ztmp1
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY, tmp, Zv1, Zv2
          integer :: I, J, imj, nf, dec, I1, J1, no_I, no_J, n, nf2, k, k1, l, l3, I0, J0
          integer :: nb, nc, amy, lly, nb_r, is, k0, l0, js, I_nn1, I_nn2, J_nn1, J_nn2
          real(Kind=kind(0.d0)) :: X, Zero = 1.0e-6

          ZP = PHASE/real(Phase, kind(0.d0))
          ZS = real(Phase, kind(0.d0))/abs(real(Phase, kind(0.d0)))

          ZS = ZS*Mc_step_weight

          do nf = 1, N_FL
             do I = 1, Ndim
                do J = 1, Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                end do
                GRC(I, I, nf) = 1.d0 + GRC(I, I, nf)
             end do
          end do
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          ! Compute scalar observables.
          do I = 1, size(Obs_scal, 1)
             Obs_scal(I)%N = Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign = Obs_scal(I)%Ave_sign + real(ZS, kind(0.d0))
          end do

          Zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call Predefined_Hoppings_Compute_Kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin*dble(N_SUN)
          Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*ZP*ZS

          Zn = cmplx(dble(N_sun), 0.d0, kind(0.d0))
          ZPot = cmplx(0.d0, 0.d0, kind(0.d0))
          Zv1 = cmplx(0.d0, 0.d0, kind(0.d0))
          Zv2 = cmplx(0.d0, 0.d0, kind(0.d0))

          do nf = 1, N_FL
             lly = size(bond_map_v1, 1)
             do nb_r = 1, lly
                nb = bond_map_v1(nb_r)
                do nc = 1, l_bond(nb)
                   I1 = bond_list(nb, nc, 2)
                   J1 = bond_list(nb, nc, 3)
                   Zv1 = Zv1 + (GRC(I1, I1, nf)*GRC(J1, J1, nf) + GRC(I1, J1, nf)*GR(I1, J1, nf)) - &
                              & 0.5d0*(GRC(I1, I1, nf) + GRC(J1, J1, nf)) + 0.25d0
                end do
             end do

             lly = size(bond_map_v2, 1)
             do nb_r = 1, lly
                nb = bond_map_v2(nb_r)
                do nc = 1, l_bond(nb)
                   I1 = bond_list(nb, nc, 2)
                   J1 = bond_list(nb, nc, 3)
                   Zv2 = Zv2 + (GRC(I1, I1, nf)*GRC(J1, J1, nf) + GRC(I1, J1, nf)*GR(I1, J1, nf)) - &
                              & 0.5d0*(GRC(I1, I1, nf) + GRC(J1, J1, nf)) + 0.25d0
                end do
             end do
          end do
          Zpot = Zn*(Ham_V*Zv1 + Ham_V2*Zv2)
          Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Zpot*ZP*ZS

          Zrho = cmplx(0.d0, 0.d0, kind(0.d0))
          do nf = 1, N_FL
             do I = 1, Ndim
                Zrho = Zrho + Grc(i, i, nf)
             end do
          end do
          Zrho = Zrho*dble(N_SUN)
          Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Zrho*ZP*ZS

          Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

          Zqah = cmplx(0.d0, 0.d0, kind(0.d0))
          is = 1
          i0 = invlist(is, 1)
          j0 = invlist(is, 2)
          nb = 0
          do nf = 1, N_FL
             do js = 1, Latt%N
                imj = latt%imj(is, js)
                do no_J = 1, Latt_Unit_qah%Norb
                   J_nn1 = latt%nnlist(js, 1, 0)
                   J_nn2 = latt%nnlist(js, 0, -1)
                   select case (no_J)
                   case (1)
                      k0 = invlist(js, 1)
                      l0 = invlist(js, 2)
                   case (2)
                      k0 = invlist(J_nn2, 2)
                      l0 = invlist(js, 1)
                   case (3)
                      k0 = invlist(J_nn1, 1)
                      l0 = invlist(J_nn2, 2)
                   case (4)
                      k0 = invlist(Js, 2)
                      l0 = invlist(J_nn1, 1)
                   end select

                   ! real
                   ztmp1 = cmplx(0.d0, 0.d0, kind(0.d0))
                   ztmp1 = ztmp1 + (grc(i0, j0, nf) - grc(j0, i0, nf))*(grc(k0, l0, nf) - grc(l0, k0, nf))
                   ztmp1 = ztmp1 + grc(i0, l0, nf)*gr(j0, k0, nf) + grc(j0, k0, nf)*gr(i0, l0, nf) - &
                       & grc(j0, l0, nf)*gr(i0, k0, nf) - grc(i0, k0, nf)*gr(j0, l0, nf)
                   Zqah = Zqah - ztmp1

                   nb = nb + 1
                end do
             end do
          end do
          zqah = zqah/dble(nb)
          Obs_scal(5)%Obs_vec(1) = Obs_scal(5)%Obs_vec(1) + zqah*ZP*ZS

          ! Standard two-point correlations
          call Predefined_Obs_eq_Green_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(1))
          call Predefined_Obs_eq_SpinSUN_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(2))
          call Predefined_Obs_eq_Den_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(3))

          ! QAH
          obs_eq(4)%N = obs_eq(4)%N + 1
          obs_eq(4)%Ave_sign = obs_eq(4)%Ave_sign + 1.d0

          nf = 1
          do is = 1, Latt%N
          do no_I = 1, Latt_Unit_qah%Norb

             I_nn1 = latt%nnlist(is, 1, 0)
             I_nn2 = latt%nnlist(is, 0, -1)
             select case (no_I)
             case (1)
                i0 = invlist(is, 1)
                j0 = invlist(is, 2)
             case (2)
                i0 = invlist(I_nn2, 2)
                j0 = invlist(is, 1)
             case (3)
                i0 = invlist(I_nn1, 1)
                j0 = invlist(I_nn2, 2)
             case (4)
                i0 = invlist(Is, 2)
                j0 = invlist(I_nn1, 1)
             end select

             do js = 1, Latt%N
             do no_J = 1, Latt_Unit_qah%Norb
                imj = latt%imj(is, js)

                J_nn1 = latt%nnlist(js, 1, 0)
                J_nn2 = latt%nnlist(js, 0, -1)
                select case (no_J)
                case (1)
                   k0 = invlist(js, 1)
                   l0 = invlist(js, 2)
                case (2)
                   k0 = invlist(J_nn2, 2)
                   l0 = invlist(js, 1)
                case (3)
                   k0 = invlist(J_nn1, 1)
                   l0 = invlist(J_nn2, 2)
                case (4)
                   k0 = invlist(Js, 2)
                   l0 = invlist(J_nn1, 1)
                end select

                ! real
                ztmp1 = cmplx(0.d0, 0.d0, kind(0.d0))
                ztmp1 = ztmp1 + (grc(i0, j0, nf) - grc(j0, i0, nf))*(grc(k0, l0, nf) - grc(l0, k0, nf))
                ztmp1 = ztmp1 + grc(i0, l0, nf)*gr(j0, k0, nf) + grc(j0, k0, nf)*gr(i0, l0, nf) - &
                    & grc(j0, l0, nf)*gr(i0, k0, nf) - grc(i0, k0, nf)*gr(j0, l0, nf)

                obs_eq(4)%Obs_Latt(imj, 1, no_I, no_J) = Obs_eq(4)%Obs_Latt(imj, 1, no_I, no_J) - ztmp1*zp*zs
             end do
             end do

          end do
          end do

       end subroutine Obser
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!-------------------------------------------------------------------
       subroutine ObserT(NT, GT0, G0T, G00, GTT, PHASE, Mc_step_weight)

          use Predefined_Obs

          implicit none

          integer, intent(IN) :: NT
  complex(Kind=kind(0.d0)), intent(IN) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL), G00(Ndim, Ndim, N_FL), GTT(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: Phase
          real(Kind=kind(0.d0)), intent(IN) :: Mc_step_weight

          !Locals
          complex(Kind=kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY
          real(Kind=kind(0.d0)) :: X
          integer :: IMJ, I, J, I1, J1, no_I, no_J, lly

          ZP = PHASE/real(Phase, kind(0.d0))
          ZS = real(Phase, kind(0.d0))/abs(real(Phase, kind(0.d0)))

          ZS = ZS*Mc_step_weight

          ! Standard two-point correlations

          call Predefined_Obs_tau_Green_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(1))
          call Predefined_Obs_tau_SpinSUN_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(2))
          call Predefined_Obs_tau_Den_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(3))

       end subroutine OBSERT

    end submodule ham_qbt_smod
