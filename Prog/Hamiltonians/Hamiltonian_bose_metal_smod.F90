    submodule(Hamiltonian_main) ham_bose_metal_smod

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

       type, extends(ham_base) :: ham_bose_metal
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: Alloc_obs
          procedure, nopass :: Obser
          procedure, nopass :: ObserT
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_bose_metal

       !#PARAMETERS START# VAR_lattice
       character(len=64) :: Model = "bose_metal"  ! Value not relevant
       character(len=64) :: Lattice_type = "bilayer_square"
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

       !#PARAMETERS START# VAR_bose_metal
       real(Kind=kind(0.d0)) :: ham_t     = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_alpha = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_chem  = 0.d0     ! Chemical potential
       real(Kind=kind(0.d0)) :: ham_U     = 4.d0     ! attractive Hubbard interaction
       integer               :: N_dope    = 0
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell
       
       Integer, allocatable :: rot_del_list(:,:), rot_list(:), del_list(:,:)
       complex (Kind=Kind(0.d0)), allocatable :: ff_s(:,:)
      
       Type (Unit_cell), Target  :: Latt_unit_p  ! Unit cell for f  correlation functions
       Integer, allocatable      :: List_p(:,:), Invlist_p(:,:)

    contains

       module subroutine Ham_Alloc_bose_metal
          allocate (ham_bose_metal::ham)
       end subroutine Ham_Alloc_bose_metal

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_bose_metal_read_write_parameters.F90"

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
          
          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot + 2*Thtrot

          ! Setup the Bravais lattice
          call Ham_Latt

          ! Setup the hopping / single-particle part
          call Ham_Hop

          ! Setup the interaction.
          call Ham_V

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
             write (unit_info, *) 't             : ', ham_t
             write (unit_info, *) 'alpha         : ', ham_alpha
             write (unit_info, *) 'Ham_U         : ', Ham_U
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
          
          integer :: i, i1, no, a0, a1, a2, b0, b1, b2, k, k1, k2, ntype, j1
          integer :: ix, iy, rix, riy, rdx, rdy, ndx, ndy, nsi, nx, ny, nc
          real    (kind=kind(0.d0)) :: pi = acos(-1.d0)
          
          ! Use predefined stuctures or set your own lattice.
          call Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)
          
          !! pairing lattice
          latt_unit_p%norb    = 1
          latt_unit_p%n_coord = 2
          allocate (Latt_Unit_p%Orb_pos_p(1,3))
          latt_unit_p%Orb_pos_p(1,:) =  0.d0
          latt_unit_p%Orb_pos_p(1,3) = -1.d0
          
          allocate (list_p(latt%N*latt_unit_p%norb,2), invlist_p(latt%n,latt_unit_p%norb)) ! For measuring only on p-lattice
          nc = 0
          do I = 1,Latt%N
             do no = 1,latt_unit_p%norb
                 nc = nc + 1
                 list_p(nc,1) = i
                 list_p(nc,2) = no
                 invlist_p(i,no) = invlist(i,no)
             enddo
          enddo
            
          allocate(rot_del_list(latt%n,4))
          allocate(del_list    (latt%n,4))
          allocate(rot_list    (latt%n  ))

          del_list=0
          rot_list=0
          rot_del_list=0

          !! rotation mapping
          do i1 = 1, latt%n
             ix  = latt%list(i1,1); iy = latt%list(i1,2)
             rix = iy; riy = -ix
             rdx = rix - ix; rdy = riy - iy
             ndx = abs(rdx); ndy = abs(rdy)
             nsi = i1
             a1 = sign(1,rdx); a2 = sign(1,rdy)
             do nx = 1, ndx;
                nsi = latt%nnlist(nsi,a1,0)
             enddo
             do ny = 1, ndy;
                nsi = latt%nnlist(nsi,0,a2)
             enddo
             rot_list(i1) = nsi
          enddo
          
          !! rotation + delta mapping
          do ntype = 1, 4
             do j1 = 1, latt%n
                select case(ntype)
                case(1)
                    i1 = latt%nnlist(j1, 1, 0)
                case(2)
                    i1 = latt%nnlist(j1,-1, 0)
                case(3)
                    i1 = latt%nnlist(j1, 0, 1)
                case(4)
                    i1 = latt%nnlist(j1, 0,-1)
                end select
                ix  = latt%list(i1,1); iy = latt%list(i1,2)
                rix = iy; riy = -ix
                rdx = rix - ix; rdy = riy - iy
                ndx = abs(rdx); ndy = abs(rdy)
                nsi = i1
                a1 = sign(1,rdx); a2 = sign(1,rdy)
                do nx = 1, ndx;
                   nsi = latt%nnlist(nsi,a1,0)
                enddo
                do ny = 1, ndy;
                   nsi = latt%nnlist(nsi,0,a2)
                enddo
                rot_del_list(j1,ntype) = nsi
             enddo
          enddo
          
          !! delta mapping
          do ntype = 1, 4
             do j1 = 1, latt%n
                select case(ntype)
                case(1)
                    i1 = latt%nnlist(j1, 1, 0)
                case(2)
                    i1 = latt%nnlist(j1,-1, 0)
                case(3)
                    i1 = latt%nnlist(j1, 0, 1)
                case(4)
                    i1 = latt%nnlist(j1, 0,-1)
                end select
                del_list(j1,ntype) = i1
             enddo
          enddo
            
          !! d wave 1, px wave 2, py wave 3
          allocate(ff_s(4,3))
          do ntype = 1, 3
              select case(ntype)
              case (1)
                  ff_s(1,ntype) =  1.d0
                  ff_s(2,ntype) =  1.d0
                  ff_s(3,ntype) = -1.d0
                  ff_s(4,ntype) = -1.d0
              case (2)
                  ff_s(1,ntype) =  1.d0
                  ff_s(2,ntype) = -1.d0
                  ff_s(3,ntype) =  0.d0
                  ff_s(4,ntype) =  0.d0
              case (3)
                  ff_s(1,ntype) =  0.d0
                  ff_s(2,ntype) =  0.d0
                  ff_s(3,ntype) =  1.d0
                  ff_s(4,ntype) = -1.d0
              end select
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

          real(Kind=kind(0.d0)), allocatable :: Ham_tx_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                Ham_ty_vec(:), ham_tp_vec(:), Ham_Lambda_vec(:)
          integer, allocatable ::   N_Phi_vec(:)

          ! Use predefined stuctures or set your own hopping
          integer :: n, nth

          allocate (ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &    N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL), ham_tp_vec(n_fl))

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_tx_vec = ham_t
          Ham_ty_vec = ham_t
          Ham_Chem_vec = ham_chem
          Phi_X_vec = phi_X
          Phi_Y_vec = phi_Y
          N_Phi_vec = n_phi

          ham_tx_vec(1) = ham_t;
          ham_ty_vec(1) = ham_t*ham_alpha;
          ham_tp_vec(1) = 0.d0

          select case (Lattice_type)
          case ("bilayer_square")
            call set_hopping_parameters_bilayer_square(Hopping_Matrix, ham_tx_vec, ham_ty_vec, ham_tp_vec, Ham_Chem_vec, &
                   & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          end select

          call Predefined_Hoppings_set_OPT(Hopping_Matrix, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_T)

          deallocate (ham_tx_vec, ham_ty_vec, ham_tp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec)

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
          N_part = ndim/2-n_dope
          call Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
               &                            N_part, ham_alpha, N_FL, WF_L, WF_R)

       end subroutine Ham_Trial

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
       subroutine Ham_V

          use Predefined_Int
          implicit none

          integer :: nf, I, I1, I2, nc, J, no, N_ops
          real(Kind=kind(0.d0)) :: X, zero = 1.d-10

          N_ops = 0
          if (abs(ham_u) > zero) N_ops = N_ops + Latt%N

          allocate (op_v(n_ops, n_fl))
          if (abs(ham_u) > zero) then
             nc = 0
             do nf = 1, n_fl
             do i1 = 1, latt%N
                 nc = nc + 1
                 i = invlist(i1, 1)
                 j = invlist(i1, 2)
                 call op_make(op_v(nc,nf),2)
                 op_v(nc,nf)%p(1) = I
                 op_v(nc,nf)%p(2) = J
                 op_v(nc,nf)%o(1,2) = cmplx(0.d0, -1.d0, kind(0.D0)) 
                 op_v(nc,nf)%o(2,1) = cmplx(0.d0,  1.d0, kind(0.D0)) 
                 op_V(nc,nf)%g      = cmplx(0.d0, sqrt(dtau*ham_u/2.d0), kind(0.D0))
                 op_v(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                 op_v(nc,nf)%type   = 2
                 call op_set(op_v(nc,nf))
             end do
             end do
          end if

       end subroutine Ham_V

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
          allocate (obs_scal(4))
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
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             call obser_Vec_make(obs_scal(I), N, Filename)
          end do

          ! Equal time correlators
          allocate (obs_eq(8))
          do I = 1, size(obs_eq, 1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "Den"
             case (3)
                Filename = "SpinZ"
             case (4)
                Filename = "swave"
             case (5)
                Filename = "sfflo"
             case (6)
                Filename = "dwave"
             case (7)
                Filename = "pxwave"
             case (8)
                Filename = "pywave"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             if ( I .ge. 3 ) then
                call obser_Latt_make(obs_eq(I), nt, Filename, Latt, Latt_unit_p, Channel, dtau)
             else
                call obser_Latt_make(obs_eq(I), nt, Filename, Latt, Latt_unit  , Channel, dtau)
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
                nt = Ltrot + 1 - 2*Thtrot
                if (Projector) Channel = 'T0'
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
       subroutine obser(gr, Phase, Ntau, Mc_step_weight)

          implicit none

          complex(kind=kind(0.d0)), intent(in) :: gr(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: PHASE
          integer, intent(IN)          :: Ntau
          real(Kind=kind(0.d0)), intent(IN) :: Mc_step_weight

          !Local
          complex(Kind=kind(0.d0)) :: grc(Ndim, Ndim, N_FL), ZK, zone, ztmp, z_ol, zero, ztmp1, ztmp2, ztmp3, ztmp4
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY, zback, zw, z_fac, z1j, cpair(4)
          integer :: I, J, k, l, m, n, imj, nf, dec, i1, ipx, ipy, imx, imy, j1, jpx, jpy, jmx, jmy, no_I, no_J, nc
          integer :: i2, j2, idelta_1, jdelta_1, idelta_2, jdelta_2
          integer :: idl(4), jdl(4), irdl(4), jrdl(4), rsi, rsj, k1, k2, rsi1, rsi2, rsj1, rsj2
          real(Kind=kind(0.d0)) :: X

          ZP = PHASE/real(Phase, kind(0.d0))
          ZS = real(Phase, kind(0.d0))/abs(real(Phase, kind(0.d0)))

          ZS = ZS*Mc_step_weight

          do nf = 1, N_FL
             do I = 1, Ndim
                do J = 1, Ndim
                   grc(I, J, nf) = -gr(J, I, nf)
                end do
                grc(I, I, nf) = 1.d0 + grc(I, I, nf)
             end do
          end do
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >
          
          ! Compute scalar observables.
          do i = 1, size(obs_scal, 1)
             obs_scal(i)%n = obs_scal(i)%n + 1
             obs_scal(i)%ave_sign = obs_scal(i)%ave_sign + real(zs, kind(0.d0))
          end do
          
          ! Compute the standard two-point correlations
          Do i = 1, size(obs_eq,1)
             obs_eq(i)%n        = obs_eq(i)%n + 1
             obs_eq(i)%ave_sign = obs_eq(i)%ave_sign + real(zs,kind(0.d0))
          enddo

          Zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call Predefined_Hoppings_Compute_Kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin*dble(N_SUN)
          Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*ZP*ZS

          zpot = cmplx(0.d0, 0.d0, kind(0.d0))
          do I = 1, Latt%N
             i1 = Invlist(i,1)
             i2 = Invlist(i,2)
             zpot = zpot + 2.d0*( grc(i1,i1,1)*grc(i2,i2,1) + grc(i1,i2,1)*gr(i1,i2,1) ) - &
                 & grc(i1,i1,1) - grc(i2,i2,1) + 1.d0
          end do
          zpot = zpot*(-ham_u/2.d0)
          obs_scal(2)%Obs_vec(1) = obs_scal(2)%obs_vec(1) + zpot*ZP*ZS

          zrho = cmplx(0.d0, 0.d0, kind(0.d0))
          do nf = 1, n_fl
             do I = 1, ndim
                zrho = zrho + grc(i, i, nf)
             end do
          end do
          zrho = zrho*dble(N_SUN)
          obs_scal(3)%obs_vec(1) = obs_scal(3)%obs_vec(1) + zrho*ZP*ZS

          obs_scal(4)%obs_vec(1) = obs_scal(4)%obs_vec(1) + (zkin + zpot)*ZP*ZS

          ! Standard two-point correlations
          z1j = cmplx(0.d0,1.d0,kind(0.d0))
          do i = 1, latt%n

              rsi = rot_list(i)
              idl(:) = del_list(i,:) 
              irdl(:) = rot_del_list(i,:) 

              do j = 1, latt%n
              
                  rsj = rot_list(j)
                  jdl(:) = del_list(j,:) 
                  jrdl(:) = rot_del_list(j,:) 

                  imj  = latt%imj(i, j)
                  
                  do k = 1, 4
                     no_i = (k-1)/2
                     no_j = k-no_i*2
                     no_i = no_i + 1
              
                     i1 = invlist(i,no_i)
                     j1 = invlist(j,no_j)
                  
                     !! Green
                     ztmp = grc(i1,j1,1)
                     obs_eq(1)%obs_Latt(imj,1,no_i,no_j) = obs_eq(1)%obs_latt(imj,1,no_i,no_j) + ztmp*ZP*ZS
                     
                     !! Den
                     ztmp = grc(i1,j1,1)*gr(i1,j1,1) + grc(i1,i1,1)*grc(j1,j1,1)
                     obs_eq(2)%obs_Latt(imj,1,no_i,no_j) = obs_eq(2)%obs_latt(imj,1,no_i,no_j) + ztmp*ZP*ZS
                  enddo

                  !! SpinZ
                  i1 = invlist(i,1); j1 = invlist(j,1) !! spin up
                  i2 = invlist(i,2); j2 = invlist(j,2) !! spin dn

                  ztmp =   grc(i1,j1,1)*gr(i1,j1,1) + grc(i1,i1,1)*grc(j1,j1,1) &
                      &  + grc(i2,j2,1)*gr(i2,j2,1) + grc(i2,i2,1)*grc(j2,j2,1) &
                      &  - grc(i1,j2,1)*gr(i1,j2,1) - grc(i1,i1,1)*grc(j2,j2,1) &
                      &  - grc(i2,j1,1)*gr(i2,j1,1) - grc(i2,i2,1)*grc(j1,j1,1)
                  obs_eq(3)%obs_Latt(imj,1,1,1) = obs_eq(3)%obs_latt(imj,1,1,1) + ztmp*ZP*ZS
                  
                  !! s wave
                  ztmp = grc(i1,j1,1)*grc(i2,j2,1) - grc(i1,j2,1)*grc(i2,j1,1)
                  obs_eq(4)%obs_Latt(imj,1,1,1) = obs_eq(4)%obs_latt(imj,1,1,1) + ztmp*ZP*ZS
                  
                  cpair(:) = cmplx(0.d0,0.d0,kind(0.d0))

                  !! sfflo
                  rsi1 = invlist(rsi,1); rsj1 = invlist(rsj,1) !! spin up
                  rsi2 = invlist(rsi,2); rsj2 = invlist(rsj,2) !! spin dn
                  cpair(1) = cpair(1) + grc(i1,j1,1)*grc(rsi2,rsj2,1) - grc(i1,rsj2,1)*grc(rsi2,j1,1) &
                      &               + grc(i2,j2,1)*grc(rsi1,rsj1,1) - grc(i2,rsj1,1)*grc(rsi1,j2,1) &
                      &               - grc(i1,j2,1)*grc(rsi2,rsj1,1) + grc(i1,rsj1,1)*grc(rsi2,j2,1) &
                      &               - grc(i2,j1,1)*grc(rsi1,rsj2,1) + grc(i2,rsj2,1)*grc(rsi1,j1,1)

                  !! 4x4 nn bonds
                  do k = 1, 16

                     k2 = (k-1)/4
                     k1 = k-k2*4
                     k2 = k2 + 1

                     !! d fflo
                     idelta_1 = invlist(irdl(k1),1); jdelta_1 = invlist(jrdl(k2),1) !! spin up
                     idelta_2 = invlist(irdl(k1),2); jdelta_2 = invlist(jrdl(k2),2) !! spin dn

                     cpair(2) = cpair(2) + 0.25d0*ff_s(k1,1)*conjg(ff_s(k2,1))* &
                         & (   grc(i1,j1,1)*grc(idelta_2,jdelta_2,1) + grc(i2,j2,1)*grc(idelta_1,jdelta_1,1) &
                         &   - grc(i1,jdelta_2,1)*grc(idelta_2,j1,1) - grc(i2,jdelta_1,1)*grc(idelta_1,j2,1) &
                         &   - grc(i1,j2,1)*grc(idelta_2,jdelta_1,1) - grc(i2,j1,1)*grc(idelta_1,jdelta_2,1) &
                         &   + grc(i1,jdelta_1,1)*grc(idelta_2,j2,1) + grc(i2,jdelta_2,1)*grc(idelta_1,j1,1) )

                     !! px wave
                     idelta_1 = invlist(idl(k1),1); jdelta_1 = invlist(jdl(k2),1) !! spin up
                     idelta_2 = invlist(idl(k1),2); jdelta_2 = invlist(jdl(k2),2) !! spin dn
                     
                     cpair(3) = cpair(3) + 0.25d0*ff_s(k1,2)*conjg(ff_s(k2,2))* &
                         & (   grc(i1,j2,1)*grc(idelta_1,jdelta_2,1) - grc(i1,jdelta_2,1)*grc(idelta_1,j2,1) &
                         &   + grc(i2,j1,1)*grc(idelta_2,jdelta_1,1) - grc(i2,jdelta_1,1)*grc(idelta_2,j1,1) &
                         &   + grc(i1,j1,1)*grc(idelta_1,jdelta_1,1) + grc(i2,j2,1)*grc(idelta_2,jdelta_2,1) &
                         &   - grc(i1,jdelta_1,1)*grc(idelta_1,j1,1) - grc(i2,jdelta_2,1)*grc(idelta_2,j2,1) )
                     
                     !! py wave
                     cpair(4) = cpair(4) + 0.25d0*ff_s(k1,3)*conjg(ff_s(k2,3))* &
                         & (   grc(i1,j2,1)*grc(idelta_1,jdelta_2,1) - grc(i1,jdelta_2,1)*grc(idelta_1,j2,1) &
                         &   + grc(i2,j1,1)*grc(idelta_2,jdelta_1,1) - grc(i2,jdelta_1,1)*grc(idelta_2,j1,1) &
                         &   + grc(i1,j1,1)*grc(idelta_1,jdelta_1,1) + grc(i2,j2,1)*grc(idelta_2,jdelta_2,1) &
                         &   - grc(i1,jdelta_1,1)*grc(idelta_1,j1,1) - grc(i2,jdelta_2,1)*grc(idelta_2,j2,1) )

                  enddo
                  obs_eq(5)%obs_Latt(imj,1,1,1) = obs_eq(5)%obs_latt(imj,1,1,1) + cpair(1)*ZP*ZS
                  obs_eq(6)%obs_Latt(imj,1,1,1) = obs_eq(6)%obs_latt(imj,1,1,1) + cpair(2)*ZP*ZS
                  obs_eq(7)%obs_Latt(imj,1,1,1) = obs_eq(7)%obs_latt(imj,1,1,1) + cpair(3)*ZP*ZS
                  obs_eq(8)%obs_Latt(imj,1,1,1) = obs_eq(8)%obs_latt(imj,1,1,1) + cpair(4)*ZP*ZS
                  
              end do
              do no_i = 1, latt_unit%norb
                  i1 = invlist(i,no_i)
                  zback = grc(i1,i1,1)
                  obs_eq(2)%obs_latt0(no_i) = obs_eq(2)%obs_Latt0(no_i) + zback*zp*zs
              enddo
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
       subroutine obsert(nt, gt0, g0t, g00, gtt, PHASE, Mc_step_weight)

          implicit none

          integer, intent(IN) :: NT
          complex(Kind=kind(0.d0)), intent(IN) :: gt0(ndim, ndim, n_fl), g0t(ndim, ndim, n_fl)
          complex(Kind=kind(0.d0)), intent(IN) :: g00(ndim, ndim, n_fl), gtt(ndim, ndim, n_fl)
          complex(Kind=kind(0.d0)), intent(IN) :: Phase
          real(Kind=kind(0.d0)), intent(IN) :: Mc_step_weight

          !Locals
          complex(Kind=kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, zone, zback, z_ol, zw, z_fac, zero
          real(Kind=kind(0.d0)) :: X
          integer :: IMJ, I, J, k, l, m, n, i1, ip1, ipx, ipy, imx, imy, j1, jpx, jpy, jmx, jmy, no_I, no_J, nc

          ZP = PHASE/real(Phase, kind(0.d0))
          ZS = real(Phase, kind(0.d0))/abs(real(Phase, kind(0.d0)))
          ZS = ZS*Mc_step_weight

          zone = cmplx(1.d0, 0.d0, kind(0.d0))

          ! Standard two-point correlations
          if (nt == 0 ) then
          do i = 1, size(obs_tau,1)
             obs_tau(i)%n        = obs_tau(i)%n + 1
             obs_tau(i)%ave_sign = obs_tau(i)%ave_sign + real(zs,kind(0.d0))
          enddo
          endif

          do i = 1, latt%n

              do j = 1, latt%n

                  imj  = latt%imj(i, j)
                  
                  do k = 1, 4
                     no_i = (k-1)/2
                     no_j = k-no_i*2
                     no_i = no_i + 1
              
                     i1 = invlist(i,no_i)
                     j1 = invlist(j,no_j)
                  
                     !! Green
                     z = gt0(i1,j1,1)
                     obs_tau(1)%obs_Latt(imj,nt+1,no_i,no_j) = obs_tau(1)%obs_latt(imj,nt+1,no_i,no_j) + z*zp*zs

                  enddo

              end do
          end do

       end subroutine obsert

    end submodule ham_bose_metal_smod
