    submodule(Hamiltonian_main) ham_cbqbt_smod

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

       implicit none

       type, extends(ham_base) :: ham_cbqbt_ob
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: Alloc_obs
          procedure, nopass :: Obser
          procedure, nopass :: ObserT
          procedure, nopass :: E0_local
          procedure, nopass :: sum_weight
          procedure, nopass :: update_fac_norm
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_cbqbt_ob

       !#PARAMETERS START# VAR_lattice
       character(len=64) :: Model = ""  ! Value not relevant
       character(len=64) :: Lattice_type = "Pi_Flux_ob"
       integer            :: L1 = 6   ! Length in direction a_1
       integer            :: L2 = 6   ! Length in direction a_2
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_Model_Generic
       !Integer              :: N_SUN        = 1        ! Number of colors
       !Integer              :: N_FL         = 1        ! Number of flavors
       !Integer              :: N_slat       = 1        ! Number of slater on trial wave function
       !Integer              :: N_wlk        = 1        ! Number of walker
       !Integer              :: ltrot        = 10       ! length of imaginary time for dynamical measure
       real(Kind=kind(0.d0)) :: Phi_X = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
       real(Kind=kind(0.d0)) :: Phi_Y = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
       logical               :: Bulk = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
       integer               :: N_Phi = 0        ! Total number of flux quanta traversing the lattice
       real(Kind=kind(0.d0)) :: Dtau = 0.1d0    ! Thereby Ltrot=Beta/dtau
       logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
       !logical              :: Symm         = .true.   ! Whether symmetrization takes place
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_cbqbt_ob
       real(Kind=kind(0.d0)) :: ham_t = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_t2 = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_t3 = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_V = 1.d0     ! V interaction
       real(Kind=kind(0.d0)) :: ham_V2 = 0.d0     ! V interaction
       real(Kind=kind(0.d0)) :: ham_V3 = 0.d0     ! V interaction
       real(Kind=kind(0.d0)) :: ham_chem = 0.d0     ! Chemical potential
       integer               :: N_dope = 0
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell
       integer, allocatable   :: bond_list(:, :, :), l_bond(:)    ! bond list
       integer, allocatable   :: bond_map_v1(:), bond_map_v2(:) ! bond list

       integer, allocatable      :: site_map(:), inv_site_map(:)

       ! 2d reference lattice
       type(Lattice), target  :: Latt_p
       type(Unit_cell), target  :: Latt_p_unit
       integer, allocatable          :: List_p(:, :), Invlist_p(:, :)
      
       ! obs
       type(Unit_cell), target  :: Latt_unit_qah, Latt_unit_bnds  ! ab sublattice 
       integer, allocatable     :: List_qah  (:,:), Invlist_qah  (:,:,:)  
       integer, allocatable     :: List_bnds(:,:), Invlist_bnds(:,:,:)

    contains

       module subroutine Ham_Alloc_cbqbt_ob
          allocate (ham_cbqbt_ob::ham)
       end subroutine Ham_Alloc_cbqbt_ob

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_cbqbt_ob_read_write_parameters.F90"

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

          call read_parameters()

          allocate (weight_k(N_wlk))
          N_grc = N_slat*N_wlk
          allocate (overlap(N_grc))
          N_wlk_mpi = N_wlk*isize_g
          N_grc_mpi = N_grc*isize_g

          ! Setup the Bravais lattice
          call Ham_Latt

          ! Setup the hopping / single-particle part
          call Ham_Hop

          ! Setup the interaction.
          call Ham_Vint

          ! Setup the trival wave function, in case of a projector approach
          call Ham_Trial

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
             write (unit_info, *) '# unit cells  : ', Latt%N
             write (unit_info, *) '# of orbitals : ', Latt_unit%Norb
             write (unit_info, *) 'Flux_1        : ', Phi_X
             write (unit_info, *) 'Flux_2        : ', Phi_Y
             if (Bulk) then
                write (unit_info, *) 'Twist as phase factor in bulk'
             else
                write (unit_info, *) 'Twist as boundary condition'
             end if
             write (unit_info, *) 'Checkerboard  : ', Checkerboard
             write (unit_info, *) 'Symm. decomp  : ', Symm
             write (unit_info, *) 'dtau,Ltrot_eff: ', dtau, Ltrot
             write (unit_info, *) 'N_SUN         : ', N_SUN
             write (unit_info, *) 'N_FL          : ', N_FL
             write (unit_info, *) 'N_wlk         : ', N_wlk
             write (unit_info, *) 'N_wlk_mpi     : ', N_wlk_mpi
             write (unit_info, *) 'N_slat        : ', N_slat
             write (unit_info, *) 'N_grc         : ', N_grc
             write (unit_info, *) 'N_grc_mpi     : ', N_grc_mpi
             write (unit_info, *) 't             : ', ham_t
             write (unit_info, *) 't1            : ', ham_t2
             write (unit_info, *) 't3            : ', ham_t3
             write (unit_info, *) 'Ham_V         : ', ham_v
             write (unit_info, *) 'Ham_V2        : ', ham_v2
             write (unit_info, *) 'Ham_V3        : ', ham_v3
             write (unit_info, *) 'Ham_chem      : ', Ham_chem
             write (unit_info, *) 'N_dope        : ', N_dope
             do nf = 1, N_FL
                write (unit_info, *) 'Degen of right trial wave function: ', wf_r(nf, 1)%degen
                write (unit_info, *) 'Degen of left  trial wave function: ', wf_l(nf, 1)%degen
             end do
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

          integer :: i, j, i1, no, a0, a1, a2, b0, b1, b2, k, k1, k2, ntype, j1, ndim_p
          integer :: ix, iy, nx, ny, nc, nc1, n, ly, ry
          integer :: no_tmp, no1_tmp, I_nn1
          real(kind=kind(0.d0)) :: pi = acos(-1.d0)
          character(len=64) :: Lattice_type_2d

          ! Use predefined stuctures or set your own lattice.
          call Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)

          Lattice_type_2d = 'Pi_Flux'
          call Predefined_Latt(Lattice_type_2d, L1, L2, ndim_p, List_p, Invlist_p, Latt_p, Latt_p_Unit)

          if (ndim_p .ne. ndim) then
             write (6, *) ' Error in set lattice mapping ', ndim, ndim_p
             error stop 1
          end if

          ! map open boundary and periodic boundary
          allocate (site_map(ndim), inv_site_map(ndim))
          Ix = 1
          do I = 1, L1
          do J = 1, L2
             do K = 1, Latt_p_Unit%Norb
                nc1 = K + (J - 1)*Latt_p_Unit%Norb
                site_map(Invlist_p(Ix, K)) = invlist(I, nc1)
                inv_site_map(invlist(I, nc1)) = invlist_p(Ix, K)
             end do
             Ix = latt_p%nnlist(Ix, 0, 1)
          end do
          Ix = latt_p%nnlist(Ix, 1, 0)
          end do
          
          !!=========================================================!!
          !!  QAH correlation
          !!=========================================================!!
          Latt_Unit_qah%norb    = (L2-1)*4+2
          Latt_Unit_qah%n_coord = 2
          allocate (Latt_unit_qah%orb_pos_p(Latt_Unit_qah%norb,3))
          Latt_unit_qah%orb_pos_p = 0.d0
          do n = 1,Latt_Unit_qah%Norb 
             Latt_Unit_qah%Orb_pos_p(n,3) = real(n-1,kind(0.d0))
          enddo

          allocate (List_qah (Latt%N*Latt_Unit_qah%Norb,2))
          allocate (Invlist_qah(Latt%N,Latt_Unit_qah%Norb,3))
          
          Ly = L2
          
          nc=0
          do I = 1, Latt%N
             
             no=0
             do ry = 2, Ly
                 !! bond 1 in plaqutte
                 no = no + 1
                 I_nn1 = latt%nnlist(I,1,0)
                 no_tmp  = (ry-1)*2+1
                 no1_tmp = (ry-1)*2+2
                 I1 = invlist(I,no_tmp )
                 J1 = invlist(I,no1_tmp)
                 
                 nc = nc + 1
                 list_qah(nc, 1) = I
                 list_qah(nc, 2) = no
                 invlist_qah(I ,no, 1) = I1
                 invlist_qah(I ,no, 2) = J1
                 invlist_qah(I ,no, 3) = nc
             
                 !! bond 2 in plaqutte
                 no = no + 1
                 I_nn1 = latt%nnlist(I,1,0)
                 no_tmp  = (ry-1)*2
                 no1_tmp = (ry-1)*2+1
                 I1 = invlist(I,no_tmp )
                 J1 = invlist(I,no1_tmp)
                 
                 nc = nc + 1
                 list_qah(nc, 1) = I
                 list_qah(nc, 2) = no
                 invlist_qah(I ,no, 1) = I1
                 invlist_qah(I ,no, 2) = J1
                 invlist_qah(I ,no, 3) = nc
             
                 !! bond 3 in plaqutte
                 no = no + 1
                 I_nn1 = latt%nnlist(I,1,0)
                 no_tmp  = (ry-1)*2+1
                 no1_tmp = (ry-1)*2
                 I1 = invlist(I_nn1,no_tmp )
                 J1 = invlist(I    ,no1_tmp)
                 
                 nc = nc + 1
                 list_qah(nc, 1) = I
                 list_qah(nc, 2) = no
                 invlist_qah(I ,no, 1) = I1
                 invlist_qah(I ,no, 2) = J1
                 invlist_qah(I ,no, 3) = nc
             
                 !! bond 4 in plaqutte
                 no = no + 1
                 I_nn1 = latt%nnlist(I,1,0)
                 no_tmp  = (ry-1)*2+2
                 no1_tmp = (ry-1)*2+1
                 I1 = invlist(I    ,no_tmp )
                 J1 = invlist(I_nn1,no1_tmp)
                 
                 nc = nc + 1
                 list_qah(nc, 1) = I
                 list_qah(nc, 2) = no
                 invlist_qah(I ,no, 1) = I1
                 invlist_qah(I ,no, 2) = J1
                 invlist_qah(I ,no, 3) = nc
             enddo

             !!=========================!!
             !! open boundary condition
             !!=========================!!
             ry = 1
             
             no = no + 1
             I_nn1 = latt%nnlist(I,1,0)
             no_tmp  = (ry-1)*2+1
             no1_tmp = (ry-1)*2+2
             I1 = invlist(I,no_tmp )
             J1 = invlist(I,no1_tmp)
             
             nc = nc + 1
             list_qah(nc, 1) = I
             list_qah(nc, 2) = no
             invlist_qah(I ,no, 1) = I1
             invlist_qah(I ,no, 2) = J1
             invlist_qah(I ,no, 3) = nc
          
             no = no + 1
             I_nn1 = latt%nnlist(I,1,0)
             no_tmp  = (ry-1)*2+2
             no1_tmp = (ry-1)*2+1
             I1 = invlist(I    ,no_tmp )
             J1 = invlist(I_nn1,no1_tmp)
             
             nc = nc + 1
             list_qah(nc, 1) = I
             list_qah(nc, 2) = no
             invlist_qah(I ,no, 1) = I1
             invlist_qah(I ,no, 2) = J1
             invlist_qah(I ,no, 3) = nc

          enddo
          
          !!=========================================================!!
          !! Bond semetic correlation
          !!=========================================================!!
          Latt_Unit_bnds%Norb    = (L2-1)*4 + 2
          Latt_Unit_bnds%N_coord = 2
          allocate (Latt_Unit_bnds%orb_pos_p(Latt_Unit_bnds%Norb,3))
          Latt_unit_bnds%orb_pos_p = 0.d0
          do n = 1,Latt_Unit_bnds%norb 
             Latt_Unit_bnds%orb_pos_p(n,3) = real(n-1,kind(0.d0))
          enddo

          allocate (List_bnds   (Latt%N*Latt_Unit_bnds%Norb,2))
          allocate (Invlist_bnds(Latt%N,Latt_Unit_bnds%Norb,3))
          
          Ly = L2
          
          nc=0
          do I = 1, Latt%N
             
             no=0
             do ry = 1, Ly-1
                 !! sub lattice a, direction x
                 no = no + 1
                 I_nn1 = latt%nnlist(I,1,0)
                 no_tmp  = (ry-1)*2+1
                 no1_tmp = (ry-1)*2+1
                 I1 = invlist(I    ,no_tmp )
                 J1 = invlist(I_nn1,no1_tmp)
                 
                 nc = nc + 1
                 list_bnds(nc, 1) = I
                 list_bnds(nc, 2) = no
                 invlist_bnds(I ,no, 1) = I1
                 invlist_bnds(I ,no, 2) = J1
                 invlist_bnds(I ,no, 3) = nc
             
                 !! sub lattice b, direction x
                 no = no + 1
                 I_nn1 = latt%nnlist(I,1,0)
                 no_tmp  = (ry-1)*2+2
                 no1_tmp = (ry-1)*2+2
                 I1 = invlist(I    ,no_tmp )
                 J1 = invlist(I_nn1,no1_tmp)
                 
                 nc = nc + 1
                 list_bnds(nc, 1) = I
                 list_bnds(nc, 2) = no
                 invlist_bnds(I ,no, 1) = I1
                 invlist_bnds(I ,no, 2) = J1
                 invlist_bnds(I ,no, 3) = nc
             
                 !! sub lattice a, direction y
                 no = no + 1
                 no_tmp  = (ry-1)*2+1
                 no1_tmp = ry*2+1
                 I1 = invlist(I, no_tmp )
                 J1 = invlist(I, no1_tmp)
                 
                 nc = nc + 1
                 list_bnds(nc, 1) = I
                 list_bnds(nc, 2) = no
                 invlist_bnds(I ,no, 1) = I1
                 invlist_bnds(I ,no, 2) = J1
                 invlist_bnds(I ,no, 3) = nc
             
                 !! sub lattice b, direction y
                 no = no + 1
                 no_tmp  = (ry-1)*2+2
                 no1_tmp = ry*2+2
                 I1 = invlist(I,no_tmp )
                 J1 = invlist(I,no1_tmp)
                 
                 nc = nc + 1
                 list_bnds(nc, 1) = I
                 list_bnds(nc, 2) = no
                 invlist_bnds(I ,no, 1) = I1
                 invlist_bnds(I ,no, 2) = J1
                 invlist_bnds(I ,no, 3) = nc
             enddo

             !!=========================!!
             !! open boundary condition
             !!=========================!!
             ry = ly
             
             !! sub lattice a, direction x
             no = no + 1
             I_nn1 = latt%nnlist(I,1,0)
             no_tmp  = (ry-1)*2+1
             no1_tmp = (ry-1)*2+1
             I1 = invlist(I    ,no_tmp )
             J1 = invlist(I_nn1,no1_tmp)
             
             nc = nc + 1
             list_bnds(nc, 1) = I
             list_bnds(nc, 2) = no
             invlist_bnds(I ,no, 1) = I1
             invlist_bnds(I ,no, 2) = J1
             invlist_bnds(I ,no, 3) = nc
          
             !! sub lattice b, direction x
             no = no + 1
             I_nn1 = latt%nnlist(I,1,0)
             no_tmp  = (ry-1)*2+2
             no1_tmp = (ry-1)*2+2
             I1 = invlist(I    ,no_tmp )
             J1 = invlist(I_nn1,no1_tmp)
             
             nc = nc + 1
             list_bnds(nc, 1) = I
             list_bnds(nc, 2) = no
             invlist_bnds(I ,no, 1) = I1
             invlist_bnds(I ,no, 2) = J1
             invlist_bnds(I ,no, 3) = nc

          enddo

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
          ham_t_vec = ham_t
          ham_t2_vec = ham_t2
          ham_t3_vec = ham_t3
          ham_tperp_vec = 0.d0
          ham_Chem_vec = ham_Chem
          Phi_X_vec = Phi_X
          Phi_Y_vec = Phi_Y
          Ham_Lambda_vec = Ham_Lambda
          N_Phi_vec = N_Phi

          select case (Lattice_type)
          case ("Pi_Flux_ob")
             call set_hopping_parameters_pi_flux_qbt_ob(Hopping_Matrix, ham_t_vec, ham_t2_vec, ham_chem_vec, &
                 & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          end select

          call Predefined_Hoppings_set_OPT(Hopping_Matrix, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_T)

          deallocate (Ham_T_vec, Ham_T2_vec, Ham_T3_vec, Ham_Tperp_vec, Ham_Chem_vec, &
              &    Phi_X_vec, Phi_Y_vec, N_Phi_vec, Ham_Lambda_vec)

       end subroutine Ham_Hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the trial wave function
!--------------------------------------------------------------------
       subroutine Ham_Trial

          use Predefined_Trial

          implicit none

          integer :: N_part, nf
          ! Use predefined stuctures or set your own Trial  wave function
          N_part = ndim/2 - n_dope
          call Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
               &                            N_part, N_FL, N_slat, WF_L, WF_R)

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
          integer                           :: dx, dy, dnb, lly, n_b_t, ly_1, ly_n
          integer                           :: n_1, n_2, Nb, n_f, l_f, n_l, N, nc, I0, J0, n_cb, nu_bond
          complex(Kind=kind(0.d0))          :: Z
          real(Kind=kind(0.d0))             :: Zero = 1.0e-6

          amx = mod(l1, 2)
          amy = mod(l2, 2)

          nu_bond = 12 + amx + 3*amy
          allocate (bond_list(nu_bond, Latt_p%N, 3))
          allocate (l_bond(nu_bond))

          l_bond = 0

          i0 = 1
          ly_1 = latt_p%list(i0, 2)
          i0 = latt_p%nnlist(i0, 0, -1)
          ly_n = latt_p%list(i0, 2)

          ! set bond list
          do n_b = 1, 12
             nc = 0
             select case (n_b)
             case (1)
                do I = 1, Latt_p%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = 0; 
                   bond_list(n_b, nc, 2) = Invlist_p(I, 1); 
                   bond_list(n_b, nc, 3) = invlist_p(I, 2)
                end do
                l_bond(n_b) = nc
             case (2)
                do I = 1, Latt_p%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = 0; 
                   bond_list(n_b, nc, 2) = Invlist_p(I, 2); 
                   bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I, 1, 0), 1)
                end do
                l_bond(n_b) = nc
             case (3)
                do I = 1, Latt_p%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = 0; 
                   if (latt_p%list(I, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                   bond_list(n_b, nc, 2) = Invlist_p(I, 2); 
                   bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I, 1, 1), 1)
                end do
                l_bond(n_b) = nc
             case (4)
                do I = 1, Latt_p%N
                   nc = nc + 1
                   bond_list(n_b, nc, 1) = 0; 
                   if (latt_p%list(I, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                   bond_list(n_b, nc, 2) = Invlist_p(I, 2); 
                   bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I, 0, 1), 1)
                end do
                l_bond(n_b) = nc
             case (5)
                i0 = 1
                do dx = 1, L1 - amx
                do dy = 1, L2
                   if (mod(latt_p%list(i0, 1), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 0), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 0), 2)
                   end if

                   i0 = latt_p%nnlist(I0, 0, 1)
                end do
                i0 = latt_p%nnlist(I0, 1, 0)
                end do
                l_bond(n_b) = nc
             case (6)
                i0 = 1
                do dx = 1, L1 - amx
                do dy = 1, L2
                   if (mod(latt_p%list(i0, 1), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 0), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 0), 2)
                   end if

                   i0 = latt_p%nnlist(I0, 0, 1)
                end do
                i0 = latt_p%nnlist(I0, 1, 0)
                end do
                l_bond(n_b) = nc
             case (7)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt_p%list(i0, 2), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 0, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 0, 1), 2)
                   end if
                   i0 = latt_p%nnlist(I0, 1, 0)
                end do
                i0 = latt_p%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (8)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt_p%list(i0, 2), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 0, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 0, 1), 2)
                   end if
                   i0 = latt_p%nnlist(I0, 1, 0)
                end do
                i0 = latt_p%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (9)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt_p%list(i0, 2), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 1), 2)
                   end if
                   i0 = latt_p%nnlist(I0, 1, 0)
                end do
                i0 = latt_p%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (10)
                i0 = 1
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt_p%list(i0, 2), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 1), 2)
                   end if
                   i0 = latt_p%nnlist(I0, 1, 0)
                end do
                i0 = latt_p%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (11)
                i0 = 1
                i0 = latt_p%nnlist(I0, 0, 1)
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt_p%list(i0, 2), 2) == 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_1) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, -1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_1) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, -1), 2)
                   end if
                   i0 = latt_p%nnlist(I0, 1, 0)
                end do
                i0 = latt_p%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             case (12)
                i0 = 1
                i0 = latt_p%nnlist(I0, 0, 1)
                do dy = 1, L2 - amy
                do dx = 1, L1
                   if (mod(latt_p%list(i0, 2), 2) .ne. 0) then
                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_1) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, -1), 1)

                      nc = nc + 1
                      bond_list(n_b, nc, 1) = 0; 
                      if (latt_p%list(I0, 2) .eq. ly_1) bond_list(n_b, nc, 1) = 1; 
                      bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                      bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, -1), 2)
                   end if
                   i0 = latt_p%nnlist(I0, 1, 0)
                end do
                i0 = latt_p%nnlist(I0, 0, 1)
                end do
                l_bond(n_b) = nc
             end select
          end do

          n_b = 12
          do npxy = 1, amx
             n_b = n_b + 1
             i0 = 1
             do dx = 1, L1 - amx
                i0 = latt_p%nnlist(I0, 1, 0)
             end do
             nc = 0
             do dy = 1, L2
                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 0), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 0), 2)

                i0 = latt_p%nnlist(I0, 0, 1)
             end do
             l_bond(n_b) = nc
          end do

          do npxy = 1, amy
             n_b = n_b + 1
             i0 = 1
             do dy = 1, L2 - amy
                i0 = latt_p%nnlist(I0, 0, 1)
             end do
             nc = 0
             do dx = 1, L1
                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 0, 1), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 0, 1), 2)

                i0 = latt_p%nnlist(I0, 1, 0)
             end do
             l_bond(n_b) = nc

             n_b = n_b + 1
             i0 = 1
             do dy = 1, L2 - amy
                i0 = latt_p%nnlist(I0, 0, 1)
             end do
             nc = 0
             do dx = 1, L1
                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 1), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                if (latt_p%list(I0, 2) .eq. ly_n) bond_list(n_b, nc, 1) = 1; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, 1), 2)

                i0 = latt_p%nnlist(I0, 1, 0)
             end do
             l_bond(n_b) = nc

             n_b = n_b + 1
             i0 = 1
             nc = 0
             do dx = 1, L1
                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                if (latt_p%list(I0, 2) .eq. ly_1) bond_list(n_b, nc, 1) = 1; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 1); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, -1), 1)

                nc = nc + 1
                bond_list(n_b, nc, 1) = 0; 
                if (latt_p%list(I0, 2) .eq. ly_1) bond_list(n_b, nc, 1) = 1; 
                bond_list(n_b, nc, 2) = Invlist_p(I0, 2); 
                bond_list(n_b, nc, 3) = invlist_p(Latt_p%nnlist(I0, 1, -1), 2)

                i0 = latt_p%nnlist(I0, 1, 0)
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
                   Op_V(nc, nf)%P(1) = site_map(I)
                   Op_V(nc, nf)%P(2) = site_map(J)
                   Op_V(nc, nf)%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 2) = cmplx(-1.d0, 0.d0, kind(0.d0))
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
                   Op_V(nc, nf)%P(1) = site_map(I)
                   Op_V(nc, nf)%P(2) = site_map(J)
                   Op_V(nc, nf)%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 2) = cmplx(-1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%g = sqrt(cmplx(0.5d0*0.5d0*dtau*Ham_V2, 0.d0, kind(0.d0)))
                   if (bond_list(n_cb, J0, 1) .eq. 1) Op_V(nc, nf)%g = cmplx(0.d0, 0.d0, kind(0.d0))
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
                   Op_V(nc, nf)%P(1) = site_map(I)
                   Op_V(nc, nf)%P(2) = site_map(J)
                   Op_V(nc, nf)%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 2) = cmplx(-1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%g = sqrt(cmplx(0.5d0*0.5d0*dtau*Ham_V2, 0.d0, kind(0.d0)))
                   if (bond_list(n_cb, J0, 1) .eq. 1) Op_V(nc, nf)%g = cmplx(0.d0, 0.d0, kind(0.d0))
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
                   Op_V(nc, nf)%P(1) = site_map(I)
                   Op_V(nc, nf)%P(2) = site_map(J)
                   Op_V(nc, nf)%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%O(2, 2) = cmplx(-1.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%g = sqrt(cmplx(0.5d0*0.5d0*dtau*Ham_V, 0.d0, kind(0.d0)))
                   if (bond_list(n_cb, J0, 1) .eq. 1) Op_V(nc, nf)%g = cmplx(0.d0, 0.d0, kind(0.d0))
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
       subroutine alloc_obs(ltau)

          implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          integer, intent(In) :: ltau
          integer    ::  i, N, Nt
          character(len=64) ::  Filename
          character(len=:), allocatable ::  Channel

          ! Scalar observables
          allocate (obs_scal(9))
          do I = 1, size(obs_scal, 1)
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
                N = 1; Filename = "qah"
             case (6)
                N = 1; Filename = "bnds"
             case (7)
                N = 1; Filename = "sni"
             case (8)
                N = ndim*ndim*n_fl; Filename = "grc"
             case (9)
                N = ndim*ndim*n_fl; Filename = "mixgrc"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             call obser_Vec_make(obs_scal(I), N, Filename)
          end do

          ! Equal time correlators
          allocate (obs_eq(4))
          do I = 1, size(obs_eq, 1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "Den"
             case (3)
                Filename = "QAH"
             case (4)
                Filename = "BNDS"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             if ( I .le. 2 ) then
                call obser_Latt_make(obs_eq(I), nt, Filename, Latt, Latt_unit, Channel, dtau)
             elseif ( I == 3 ) then
                call obser_Latt_make(obs_eq(I), nt, Filename, Latt, Latt_unit_qah, Channel, dtau)
             elseif ( I == 4 ) then
                call obser_Latt_make(obs_eq(I), nt, Filename, Latt, Latt_unit_bnds, Channel, dtau)
             endif
          end do

          if (Ltau == 1) then
             ! Time-displaced correlators
             allocate (obs_tau(1))
             do I = 1, size(obs_tau, 1)
                select case (I)
                case (1)
                   Channel = 'P'; Filename = "Green"
                case default
                   write (6, *) ' Error in Alloc_obs '
                end select
                nt = ltrot + 1
                channel = 'T0'
                call obser_Latt_make(obs_tau(I), nt, Filename, Latt, Latt_unit, Channel, dtau)
             end do
          end if

       end subroutine alloc_obs

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
       subroutine obser(gr, gr_mix, i_grc, re_w, sum_w, sum_o)

          implicit none

          complex(kind=kind(0.d0)), intent(in) :: gr(Ndim, Ndim, N_FL)
          complex(kind=kind(0.d0)), intent(in) :: gr_mix(Ndim, Ndim, N_FL)
          complex(kind=kind(0.d0)), intent(in) :: sum_w, sum_o
          real(kind=kind(0.d0)), intent(in) :: re_w
          integer, intent(in) :: i_grc

          !Local
          complex(Kind=kind(0.d0)) :: grc(Ndim, Ndim, N_FL), ZK, zone, ztmp, z_ol, zero, ztmp1, ztmp2, ztmp3, ztmp4
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, zqah, zbnds, zsni, zback, zw, z_fac, z1j, zn, zv1, zv2, ztmp5
          integer :: I, J, k, l, m, n, imj, nf, dec, i1, j1, no_I, no_J, nc, i3, j3, m2, n2, m3, n3
          integer :: i2, j2, lly, nb_r, nb, i0, j0, m1, n1, nb_qah, ipx, jpx, ipy, jpy, nc1, nc2
          real(Kind=kind(0.d0)) :: X, rsgn

          Z_ol = exp(overlap(i_grc))/sum_o
          ZW = cmplx(re_w, 0.d0, kind(0.d0))/sum_w
          Z_fac = Z_ol*ZW

          do nf = 1, N_FL
             do I = 1, Ndim
                do J = 1, Ndim
                   grc(I, J, nf) = -gr(J, I, nf)
                end do
                grc(I, I, nf) = 1.d0 + grc(I, I, nf)
             end do
          end do
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          Zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call Predefined_Hoppings_Compute_Kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin*dble(N_SUN)
          Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*Z_fac

          zn = cmplx(dble(N_sun), 0.d0, kind(0.d0))
          zpot = cmplx(0.d0, 0.d0, kind(0.d0))
          zv1 = cmplx(0.d0, 0.d0, kind(0.d0))
          zv2 = cmplx(0.d0, 0.d0, kind(0.d0))

          do nf = 1, N_FL
             lly = size(bond_map_v1, 1)
             do nb_r = 1, lly
                nb = bond_map_v1(nb_r)
                do nc = 1, l_bond(nb)
                   I0 = bond_list(nb, nc, 2)
                   J0 = bond_list(nb, nc, 3)
                   I1 = site_map(i0)
                   J1 = site_map(j0)
                   if (bond_list(nb, nc, 1) .ne. 1) then
                      zv1 = zv1 + (grc(I1, I1, nf)*grc(J1, J1, nf) + grc(I1, J1, nf)*gr(I1, J1, nf)) - &
                                & 0.5d0*(grc(I1, I1, nf) + grc(J1, J1, nf)) + 0.25d0
                   end if
                end do
             end do

             lly = size(bond_map_v2, 1)
             do nb_r = 1, lly
                nb = bond_map_v2(nb_r)
                do nc = 1, l_bond(nb)
                   I0 = bond_list(nb, nc, 2)
                   J0 = bond_list(nb, nc, 3)
                   I1 = site_map(i0)
                   J1 = site_map(j0)
                   if (bond_list(nb, nc, 1) .ne. 1) then
                      zv2 = zv2 + (grc(I1, I1, nf)*grc(J1, J1, nf) + grc(I1, J1, nf)*gr(I1, J1, nf)) - &
                                 & 0.5d0*(grc(I1, I1, nf) + grc(J1, J1, nf)) + 0.25d0
                   end if
                end do
             end do
          end do
          zpot = zn*(ham_v*zv1 + ham_v2*zv2)
          obs_scal(2)%obs_vec(1) = obs_scal(2)%obs_vec(1) + zpot*z_fac

          zrho = cmplx(0.d0, 0.d0, kind(0.d0))
          do nf = 1, n_fl
             do I = 1, ndim
                zrho = zrho + grc(i, i, nf)
             end do
          end do
          zrho = zrho*dble(N_SUN)
          obs_scal(3)%obs_vec(1) = obs_scal(3)%obs_vec(1) + zrho*z_fac

          obs_scal(4)%obs_vec(1) = obs_scal(4)%obs_vec(1) + (zkin + zpot)*z_fac

          !! set reference site
          zqah  = cmplx(0.d0,0.d0, kind(0.D0))
          i = 1; 
          nb_qah = latt_unit_qah%norb/2

          i1 = invlist_qah(i, nb_qah, 1)
          j1 = invlist_qah(i, nb_qah, 2)
          
          do j = 1,latt%N
             do no_j = 1, Latt_Unit_qah%Norb
                m1 = invlist_qah(j, no_j, 1)
                n1 = invlist_qah(j, no_j, 2)
                
                ztmp =   grc(i1,j1,1)*grc(m1,n1,1) + grc(i1,n1,1)*gr(j1,m1,1) &
                    &  + grc(j1,i1,1)*grc(n1,m1,1) + grc(j1,m1,1)*gr(i1,n1,1) &
                    &  - grc(i1,j1,1)*grc(n1,m1,1) - grc(i1,m1,1)*gr(j1,n1,1) & 
                    &  - grc(j1,i1,1)*grc(m1,n1,1) - grc(j1,n1,1)*gr(i1,m1,1)  
                zqah =  zqah - ztmp
             enddo
          enddo
          zqah = zqah/dble(latt_unit_qah%norb*latt%n) 
          obs_scal(5)%obs_vec(1) = obs_scal(5)%obs_vec(1) + zqah*z_fac
          
          zbnds = cmplx(0.d0,0.d0, kind(0.D0))
          zsni  = cmplx(0.d0,0.d0, kind(0.D0))
          nc1 = 0; nc2 = 0
          do i = 1,latt%N
          do j = 1,latt%N
             do k = 1, latt_unit%norb*latt_unit%norb
                no_i = (k - 1)/latt_unit%norb
                no_j = k - no_i*latt_unit%norb
                no_i = no_i + 1
                   
                i0 = invlist(i, no_i)   
                j0 = invlist(j, no_j)

                ipx = invlist(latt%nnlist(i,1,0),no_i)
                jpx = invlist(latt%nnlist(j,1,0),no_j)
                
                if ( mod(no_i+no_j,2) .eq. 0 ) then

                i1 = i0; j1 = ipx
                m1 = j0; n1 = jpx
                ztmp1 =   grc(i1,j1,1)*grc(m1,n1,1) + grc(i1,n1,1)*gr(j1,m1,1) &
                    &   + grc(j1,i1,1)*grc(n1,m1,1) + grc(j1,m1,1)*gr(i1,n1,1) &
                    &   + grc(i1,j1,1)*grc(n1,m1,1) + grc(i1,m1,1)*gr(j1,n1,1) & 
                    &   + grc(j1,i1,1)*grc(m1,n1,1) + grc(j1,n1,1)*gr(i1,m1,1)  
                nc1 = nc1 + 1
                
                ztmp4 = 0.d0
                if ( (no_i+2) .le. latt_unit%norb ) then
                    ipy = invlist(i, no_i+2)
                    i1 = i0; j1 = ipy
                    m1 = j0; n1 = jpx
                    ztmp4 =   grc(i1,j1,1)*grc(m1,n1,1) + grc(i1,n1,1)*gr(j1,m1,1) &
                        &   + grc(j1,i1,1)*grc(n1,m1,1) + grc(j1,m1,1)*gr(i1,n1,1) &
                        &   + grc(i1,j1,1)*grc(n1,m1,1) + grc(i1,m1,1)*gr(j1,n1,1) & 
                        &   + grc(j1,i1,1)*grc(m1,n1,1) + grc(j1,n1,1)*gr(i1,m1,1)  
                    nc1 = nc1 + 1
                endif
                
                ztmp3 = 0.d0
                if ( (no_j+2) .le. latt_unit%norb ) then
                    jpy = invlist(j, no_j+2)
                    i1 = i0; j1 = ipx
                    m1 = j0; n1 = jpy
                    ztmp3 =   grc(i1,j1,1)*grc(m1,n1,1) + grc(i1,n1,1)*gr(j1,m1,1) &
                        &   + grc(j1,i1,1)*grc(n1,m1,1) + grc(j1,m1,1)*gr(i1,n1,1) &
                        &   + grc(i1,j1,1)*grc(n1,m1,1) + grc(i1,m1,1)*gr(j1,n1,1) & 
                        &   + grc(j1,i1,1)*grc(m1,n1,1) + grc(j1,n1,1)*gr(i1,m1,1)  
                    nc1 = nc1 + 1
                endif
                
                ztmp2 = 0.d0
                if ( ( (no_i+2) .le. latt_unit%norb ) .and. ( (no_j+2) .le. latt_unit%norb ) ) then
                    i1 = i0; j1 = ipy
                    m1 = j0; n1 = jpy
                    ztmp2 =   grc(i1,j1,1)*grc(m1,n1,1) + grc(i1,n1,1)*gr(j1,m1,1) &
                        &   + grc(j1,i1,1)*grc(n1,m1,1) + grc(j1,m1,1)*gr(i1,n1,1) &
                        &   + grc(i1,j1,1)*grc(n1,m1,1) + grc(i1,m1,1)*gr(j1,n1,1) & 
                        &   + grc(j1,i1,1)*grc(m1,n1,1) + grc(j1,n1,1)*gr(i1,m1,1)  
                    nc1 = nc1 + 1
                endif
                
                zbnds =  zbnds + ztmp1 + ztmp2 - ztmp3 - ztmp4

                endif

                ztmp5 = grc(i0,i0,1)*grc(j0,j0,1) + grc(i0,j0,1)*gr(i0,j0,1)
                rsgn=1.d0-dble(mod(no_i+no_j,2))*2.d0
                zsni = zsni + rsgn*ztmp5
                nc2 = nc2 + 1
             enddo
          enddo
          enddo
          zbnds = zbnds/dble(nc1)
          obs_scal(6)%obs_vec(1) = obs_scal(6)%obs_vec(1) + zbnds*z_fac
          zsni = zsni/dble(nc2)
          obs_scal(7)%obs_vec(1) = obs_scal(7)%obs_vec(1) + zsni*z_fac

          nc = 0
          do nf = 1, n_fl
             do I = 1, ndim
             do J = 1, ndim
                nc = nc + 1
                ztmp = grc(i, j, nf)
                obs_scal(8)%obs_vec(nc) = obs_scal(5)%obs_vec(nc) + ztmp*z_fac
             end do
             end do
          end do

          nc = 0
          do nf = 1, n_fl
             do I = 1, ndim
             do J = 1, ndim
                nc = nc + 1
                zone = cmplx(0.d0, 0.d0, kind(0.d0))
                if (I .eq. J) zone = cmplx(1.d0, 0.d0, kind(0.d0))
                ztmp = zone - gr_mix(J, I, nf)
                obs_scal(9)%obs_vec(nc) = obs_scal(6)%obs_vec(nc) + ztmp*z_fac
             end do
             end do
          end do

          ! Standard two-point correlations
          z1j = cmplx(0.d0, 1.d0, kind(0.d0))
          do i = 1, latt%n

             do j = 1, latt%n

                imj = latt%imj(i, j)

                do k = 1, latt_unit%norb*latt_unit%norb
                   no_i = (k - 1)/latt_unit%norb
                   no_j = k - no_i*latt_unit%norb
                   no_i = no_i + 1

                   i1 = invlist(i, no_i)
                   j1 = invlist(j, no_j)

                     !! Green
                   ztmp = grc(i1, j1, 1)
                   obs_eq(1)%obs_Latt(imj, 1, no_i, no_j) = obs_eq(1)%obs_latt(imj, 1, no_i, no_j) + ztmp*z_fac

                     !! Den
                   ztmp = grc(i1, j1, 1)*gr(i1, j1, 1) + grc(i1, i1, 1)*grc(j1, j1, 1)
                   obs_eq(2)%obs_Latt(imj, 1, no_i, no_j) = obs_eq(2)%obs_latt(imj, 1, no_i, no_j) + ztmp*z_fac
                end do
                
                do k = 1, latt_unit_qah%norb*latt_unit_qah%norb
                   no_i = (k - 1)/latt_unit_qah%norb
                   no_j = k - no_i*latt_unit_qah%norb
                   no_i = no_i + 1

                   !! QAH
                   i1 = invlist_qah(i, no_i, 1)
                   j1 = invlist_qah(i, no_i, 2)
                   m1 = invlist_qah(j, no_j, 1)
                   n1 = invlist_qah(j, no_j, 2)
                   ztmp =   grc(i1,j1,1)*grc(m1,n1,1) + grc(i1,n1,1)*gr(j1,m1,1) &
                       &  + grc(j1,i1,1)*grc(n1,m1,1) + grc(j1,m1,1)*gr(i1,n1,1) &
                       &  - grc(i1,j1,1)*grc(n1,m1,1) - grc(i1,m1,1)*gr(j1,n1,1) & 
                       &  - grc(j1,i1,1)*grc(m1,n1,1) - grc(j1,n1,1)*gr(i1,m1,1)  
                   
                   obs_eq(3)%obs_Latt(imj, 1, no_i, no_j) = obs_eq(3)%obs_latt(imj, 1, no_i, no_j) - ztmp*z_fac
                   
                   !! BNDS
                   i1 = invlist_bnds(i, no_i, 1)
                   j1 = invlist_bnds(i, no_i, 2)
                   m1 = invlist_bnds(j, no_j, 1)
                   n1 = invlist_bnds(j, no_j, 2)
                   ztmp =   grc(i1,j1,1)*grc(m1,n1,1) + grc(i1,n1,1)*gr(j1,m1,1) &
                       &  + grc(j1,i1,1)*grc(n1,m1,1) + grc(j1,m1,1)*gr(i1,n1,1) &
                       &  + grc(i1,j1,1)*grc(n1,m1,1) + grc(i1,m1,1)*gr(j1,n1,1) & 
                       &  + grc(j1,i1,1)*grc(m1,n1,1) + grc(j1,n1,1)*gr(i1,m1,1)  
                   
                   obs_eq(4)%obs_Latt(imj, 1, no_i, no_j) = obs_eq(4)%obs_latt(imj, 1, no_i, no_j) + ztmp*z_fac
                enddo

             end do
             do no_i = 1, latt_unit%norb
                i1 = invlist(i, no_i)
                zback = grc(i1, i1, 1)
                obs_eq(2)%obs_latt0(no_i) = obs_eq(2)%obs_Latt0(no_i) + zback*z_fac
             end do
          end do

       end subroutine obser
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
       subroutine obsert(nt, gt0, g0t, g00, gtt, i_grc, re_w, sum_w, sum_o)

          implicit none

          integer, intent(in) :: nt
          complex(kind=kind(0.d0)), intent(in) :: gt0(ndim, ndim, n_fl), g0t(ndim, ndim, n_fl)
          complex(kind=kind(0.d0)), intent(in) :: g00(ndim, ndim, n_fl), gtt(ndim, ndim, n_fl)
          complex(kind=kind(0.d0)), intent(in) :: sum_w, sum_o
          real(kind=kind(0.d0)), intent(in) :: re_w
          integer, intent(in) :: i_grc

          !Locals
          complex(Kind=kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, zone, zback, z_ol, zw, z_fac, zero
          real(Kind=kind(0.d0)) :: X
          integer :: IMJ, I, J, k, l, m, n, i1, ip1, ipx, ipy, imx, imy, j1, jpx, jpy, jmx, jmy, no_I, no_J, nc

          Z_ol = exp(overlap(i_grc))/sum_o
          ZW = cmplx(re_w, 0.d0, kind(0.d0))/sum_w
          Z_fac = Z_ol*ZW

          zone = cmplx(1.d0, 0.d0, kind(0.d0))

          ! Standard two-point correlations
          do i = 1, latt%n

             do j = 1, latt%n

                imj = latt%imj(i, j)

                do k = 1, latt_unit%norb*latt_unit%norb
                   no_i = (k - 1)/latt_unit%norb
                   no_j = k - no_i*latt_unit%norb
                   no_i = no_i + 1

                   i1 = invlist(i, no_i)
                   j1 = invlist(j, no_j)

                     !! Green
                   z = gt0(i1, j1, 1)
                   obs_tau(1)%obs_Latt(imj, nt + 1, no_i, no_j) = obs_tau(1)%obs_latt(imj, nt + 1, no_i, no_j) + z*z_fac

                end do

             end do
          end do

       end subroutine obsert

       !!========================================================================!!
       !!     compute local energy of a given walker
       !!========================================================================!!
       complex(Kind=kind(0.d0)) function E0_local(gr)
          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: gr(ndim, ndim, n_fl)

          !Local
          complex(Kind=kind(0.d0)) :: grc(ndim, ndim, n_fl), ZK, zn, zv1, zv2
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY
          integer :: I, J, imj, nf, dec, I1, J1, i2, no_I, no_J, n, i0, j0
          integer :: nb, nb_r, lly, nc
          real(Kind=kind(0.d0)) :: X

          do nf = 1, N_FL
             do I = 1, Ndim
                do J = 1, Ndim
                   grc(I, J, nf) = -gr(J, I, nf)
                end do
                grc(I, I, nf) = 1.d0 + grc(I, I, nf)
             end do
          end do

          zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call predefined_hoppings_compute_kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, GRC, ZKin)
          zkin = zkin*dble(n_sun)

          zn = cmplx(dble(N_sun), 0.d0, kind(0.d0))
          zpot = cmplx(0.d0, 0.d0, kind(0.d0))
          zv1 = cmplx(0.d0, 0.d0, kind(0.d0))
          zv2 = cmplx(0.d0, 0.d0, kind(0.d0))

          do nf = 1, N_FL
             lly = size(bond_map_v1, 1)
             do nb_r = 1, lly
                nb = bond_map_v1(nb_r)
                do nc = 1, l_bond(nb)
                   I0 = bond_list(nb, nc, 2)
                   J0 = bond_list(nb, nc, 3)
                   I1 = site_map(i0)
                   J1 = site_map(j0)
                   if (bond_list(nb, nc, 1) .ne. 1) then
                      zv1 = zv1 + (grc(I1, I1, nf)*grc(J1, J1, nf) + grc(I1, J1, nf)*gr(I1, J1, nf)) - &
                                & 0.5d0*(grc(I1, I1, nf) + grc(J1, J1, nf))
                   end if
                end do
             end do

             lly = size(bond_map_v2, 1)
             do nb_r = 1, lly
                nb = bond_map_v2(nb_r)
                do nc = 1, l_bond(nb)
                   I0 = bond_list(nb, nc, 2)
                   J0 = bond_list(nb, nc, 3)
                   I1 = site_map(i0)
                   J1 = site_map(j0)
                   if (bond_list(nb, nc, 1) .ne. 1) then
                      zv2 = zv2 + (grc(I1, I1, nf)*grc(J1, J1, nf) + grc(I1, J1, nf)*gr(I1, J1, nf)) - &
                                 & 0.5d0*(grc(I1, I1, nf) + grc(J1, J1, nf))
                   end if
                end do
             end do
          end do
          zpot = zn*(ham_v*zv1 + ham_v2*zv2)

          E0_local = zpot + zkin

       end function E0_local

       !!===================================================================!!
       !!     compute the sum of the weight and rescale weight
       !!===================================================================!!
       subroutine sum_weight(z_sum_weight)
          implicit none

          complex(Kind=kind(0.d0)), intent(out) :: z_sum_weight

          !local
          integer :: i_wlk, ii, i_st, i_ed
          real(Kind=kind(0.d0)) :: X1, tot_re_weight, max_re_w, ang_w, re_lw
          complex(Kind=kind(0.d0)) :: Z1, Z2, wtmp

          complex(Kind=kind(0.d0)), allocatable :: weight_mpi(:), w_arr(:)

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif

          if (irank_g .eq. 0) then
             allocate (weight_mpi(N_wlk_mpi))
             allocate (w_arr(N_wlk))
             i_st = 1
             i_ed = N_wlk
             weight_mpi(i_st:i_ed) = weight_k(:)
          end if

          if (irank_g .eq. 0) then
             do ii = 1, isize_g - 1
                i_st = ii*N_wlk + 1
                i_ed = (ii + 1)*N_wlk
                call mpi_recv(w_arr, n_wlk, MPI_COMPLEX16, ii, 1, group_comm, status, ierr)
                weight_mpi(i_st:i_ed) = w_arr
             end do
          else
             call mpi_send(weight_k, n_wlk, MPI_COMPLEX16, 0, 1, group_comm, ierr)
          end if

          if (irank_g .eq. 0) then
             max_re_w = maxval(dble(weight_mpi))
             deallocate (weight_mpi, w_arr)
          end if
          call MPI_BCAST(max_re_w, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

          Z1 = cmplx(0.d0, 0.d0, kind(0.d0))
          Z2 = cmplx(0.d0, 0.d0, kind(0.d0))
          do i_wlk = 1, N_wlk
             ang_w = aimag(weight_k(i_wlk))
             !! Real part of the weight
             if (cos(ang_w) .gt. 0.d0) then
                weight_k(i_wlk) = weight_k(i_wlk) - max_re_w
                re_lw = dble(weight_k(i_wlk))
                z1 = z1 + exp(re_lw)*cos(ang_w)
             end if
          end do
          call MPI_REDUCE(Z1, Z2, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, IERR)
          call MPI_BCAST(Z2, 1, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr)

          z_sum_weight = Z2

       end subroutine sum_weight

       !!===================================================================!!
       !!     update the normalization factor
       !!===================================================================!!
       subroutine update_fac_norm(gr, ntw)
          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: gr(ndim, ndim, n_fl, n_grc)
          integer, intent(IN) :: ntw

          !local
          integer :: i_wlk, ii, i_st, i_ed, ns, i_grc
          real(Kind=kind(0.d0)) :: X1, re_w_tmp, max_re_w, ang_w, re_lw
          complex(Kind=kind(0.d0)) :: Z1, Z2, wtmp, el_tmp, Z, tot_ene, zr1, zr2
          character(LEN=64)  :: filename
          complex(Kind=kind(0.d0)), allocatable :: weight_mpi(:), w_arr(:)

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif
          filename = "E_local.dat"

          if (irank_g .eq. 0) then
             allocate (weight_mpi(N_wlk_mpi))
             allocate (w_arr(N_wlk))
             i_st = 1
             i_ed = N_wlk
             weight_mpi(i_st:i_ed) = weight_k(:)
          end if

          if (irank_g .eq. 0) then
             do ii = 1, isize_g - 1
                i_st = ii*N_wlk + 1
                i_ed = (ii + 1)*N_wlk
                call mpi_recv(w_arr, n_wlk, MPI_COMPLEX16, ii, 1, group_comm, status, ierr)
                weight_mpi(i_st:i_ed) = w_arr
             end do
          else
             call mpi_send(weight_k, n_wlk, MPI_COMPLEX16, 0, 1, group_comm, ierr)
          end if

          if (irank_g .eq. 0) then
             max_re_w = maxval(dble(weight_mpi))
             deallocate (weight_mpi, w_arr)
          end if
          call MPI_BCAST(max_re_w, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

          Z1 = cmplx(0.d0, 0.d0, kind(0.d0))
          Z2 = cmplx(0.d0, 0.d0, kind(0.d0))
          do i_wlk = 1, N_wlk
             ang_w = aimag(weight_k(i_wlk))
             !! Real part of the weight
             if (cos(ang_w) .gt. 0.d0) then
                weight_k(i_wlk) = weight_k(i_wlk) - max_re_w
                re_lw = dble(weight_k(i_wlk))

                z = 0.d0
                do ns = 1, N_slat
                   i_grc = ns + (i_wlk - 1)*N_slat
                   z = z + exp(overlap(i_grc))
                end do

                 !! real part of mix estimated energy
                tot_ene = cmplx(0.d0, 0.d0, kind(0.d0))
                do ns = 1, N_slat
                   i_grc = ns + (i_wlk - 1)*N_slat
                   el_tmp = dble(ham%E0_local(GR(:, :, :, i_grc)))
                   tot_ene = tot_ene + el_tmp*exp(overlap(i_grc))/Z
                end do
                re_w_tmp = exp(re_lw)*cos(ang_w)
                z1 = z1 + re_w_tmp
                z2 = z2 + re_w_tmp*tot_ene
             end if
          end do
          call mpi_reduce(z1, zr1, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, ierr)
          call mpi_reduce(z2, zr2, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, ierr)

          if (irank_g == 0) then
             zr2 = zr2/zr1
             fac_norm = dble(zr2)*dtau
             open (UNIT=77, FILE=filename, STATUS='UNKNOWN', position="append")
             write (77, *) ntw*dtau, dble(zr2)
             close (77)
          end if
          call MPI_BCAST(fac_norm, 1, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr)

       end subroutine update_fac_norm

    end submodule ham_cbqbt_smod
