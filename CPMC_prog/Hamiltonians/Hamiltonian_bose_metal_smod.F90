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
          procedure, nopass :: E0_local
          procedure, nopass :: sum_weight
          procedure, nopass :: update_fac_norm
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_bose_metal

       !#PARAMETERS START# VAR_lattice
       character(len=64) :: Model = "bose_metal"  ! Value not relevant
       character(len=64) :: Lattice_type = "square_anisotropic"
       integer            :: L1 = 6   ! Length in direction a_1
       integer            :: L2 = 6   ! Length in direction a_2
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_Model_Generic
       !Integer              :: N_SUN        = 1        ! Number of colors
       !Integer              :: N_FL         = 2        ! Number of flavors
       !Integer              :: N_slat       = 1        ! Number of slater on trial wave function
       !Integer              :: N_wlk        = 1        ! Number of walker
       !Integer              :: ltrot        = 10       ! length of imaginary time for dynamical measure
       real(Kind=kind(0.d0)) :: Phi_X = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
       real(Kind=kind(0.d0)) :: Phi_Y = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
       logical               :: Bulk = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
       integer               :: N_Phi = 0        ! Total number of flux quanta traversing the lattice
       real(Kind=kind(0.d0)) :: Dtau = 0.1d0    ! Thereby Ltrot=Beta/dtau
       real(Kind=kind(0.d0)) :: Beta = 5.d0     ! Inverse temperature
       logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
       !logical              :: Symm         = .true.   ! Whether symmetrization takes place
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
          call Ham_V

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
             write (unit_info, *) 'alpha         : ', ham_alpha
             write (unit_info, *) 'Ham_U         : ', Ham_U
             write (unit_info, *) 'Ham_chem      : ', Ham_chem
             write (unit_info, *) 'N_dope        : ', N_dope
             do nf = 1, N_FL
                write (unit_info, *) 'Flavor of:  ', nf
                write (unit_info, *) 'Degen of right trial wave function: ', wf_r(nf, 1)%degen
                write (unit_info, *) 'Degen of left  trial wave function: ', wf_l(nf, 1)%degen
                write (unit_info, *) 'number of particle for flavor:  ', size(wf_l(nf,1)%P,2)
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
          ! Use predefined stuctures or set your own lattice.
          call Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)

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
               &                                Ham_ty_vec(:), Ham_Lambda_vec(:)
          integer, allocatable ::   N_Phi_vec(:)

          ! Use predefined stuctures or set your own hopping
          integer :: n, nth

          allocate (ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &    N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL))

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_tx_vec = ham_t
          Ham_ty_vec = ham_t
          Ham_Chem_vec = ham_chem
          Phi_X_vec = phi_X
          Phi_Y_vec = phi_Y
          N_Phi_vec = n_phi

          ham_tx_vec(1) = ham_t;
          ham_tx_vec(2) = ham_t*ham_alpha;
          ham_ty_vec(1) = ham_t*ham_alpha;
          ham_ty_vec(2) = ham_t;

          select case (Lattice_type)
          case ("square_anisotropic")
            call set_hopping_parameters_square_anisotropic(Hopping_Matrix, ham_tx_vec, ham_ty_vec, Ham_Chem_vec, &
                   & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          end select

          call Predefined_Hoppings_set_OPT(Hopping_Matrix, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_T)

          deallocate (ham_tx_vec, ham_ty_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec)

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
          ! Use predefined stuctures or set your own Trial  wave function
          N_part = ndim/2-n_dope
          call Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
               &                            N_part, ham_alpha, N_FL, N_slat, WF_L, WF_R)

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
          if (abs(ham_u) > zero) N_ops = N_ops + Latt%N*Latt_unit%Norb

          allocate (op_v(n_ops, n_fl))
          nc = 0
          do i1 = 1, latt%N
             do no = 1, latt_unit%norb
                i = invlist(i1, no)
                if (abs(ham_u) > zero) then
                   nc = nc + 1
                   do nf = 1, n_fl
                       call op_make( op_v(nc,nf), 1 )
                       op_v(nc,nf)%p(1) = I
                       op_v(nc,nf)%o(1,1) = cmplx(1.d0, 0.d0, kind(0.D0)) 
                       op_V(nc,nf)%g      = sqrt(cmplx(dtau*ham_u/2.d0, 0.D0, kind(0.D0)))
                       op_v(nc,nf)%alpha  = cmplx(-0.5d0, 0.d0, kind(0.D0))
                       op_v(nc,nf)%type   = 2
                       call op_set(op_v(nc,nf))
                   enddo
                end if
             end do
          end do

       end subroutine Ham_V

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
       subroutine alloc_obs(ltau, lmetropolis)

          implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          integer, intent(In) :: ltau, lmetropolis
          integer    ::  i, N, Nt
          character(len=64) ::  Filename
          character(len=:), allocatable ::  Channel

          ! Scalar observables
          allocate (Obs_scal(6))
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
                N = ndim*ndim*n_fl; Filename = "grc"
             case (6)
                N = ndim*ndim*n_fl; Filename = "mixgrc"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             call Obser_Vec_make(Obs_scal(I), N, Filename)
          end do

          ! Equal time correlators
          allocate (Obs_eq(9))
          do I = 1, size(Obs_eq, 1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "Den"
             case (4)
                Filename = "swave"
             case (5)
                Filename = "dxywave"
             case (6)
                Filename = "dwave"
             case (7)
                Filename = "dswave"
             case (8)
                Filename = "pxpywave"
             case (9)
                Filename = "pxipywave"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          end do

          if (Ltau == 1) then
             ! Time-displaced correlators
             allocate (Obs_tau(4))
             do I = 1, size(Obs_tau, 1)
                select case (I)
                case (1)
                   Channel = 'P'; Filename = "Green"
                case (2)
                   Channel = 'PH'; Filename = "SpinZ"
                case (3)
                   Channel = 'PH'; Filename = "Den"
                case (4)
                   Channel = 'PH'; Filename = "swave"
                case default
                   write (6, *) ' Error in Alloc_obs '
                end select
                Nt = Ltrot + 1
                Channel = 'T0'
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
       subroutine obser(GR, GR_mix, i_wlk, i_grc, sum_w, sum_o, act_mea)

          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: GR_mix(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
          integer, intent(IN) :: i_wlk, i_grc, act_mea

          !Local
          complex(Kind=kind(0.d0)) :: grc(Ndim, Ndim, N_FL), ZK, zone, ztmp, z_ol, zero, ztmp1, ztmp2, ztmp3, ztmp4
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY, zback, zw, z_fac, z1j
          integer :: I, J, k, l, m, n, imj, nf, dec, i1, ipx, ipy, imx, imy, j1, jpx, jpy, jmx, jmy, no_I, no_J, nc
          integer :: ipxpy, ipxmy, imxpy, imxmy, jpxpy, jpxmy, jmxpy, jmxmy
          real(Kind=kind(0.d0)) :: X

          Z_ol = exp(overlap(i_grc))/sum_o
          ZW = cmplx(weight_k(i_wlk), 0.d0, kind(0.d0))/sum_w
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

          if (act_mea .eq. 0) then
             ! Compute scalar observables.
             do i = 1, size(obs_scal, 1)
                obs_scal(i)%n = obs_scal(i)%n + 1
                obs_scal(i)%ave_sign = obs_scal(i)%ave_sign + 1.d0
             end do

             ! Compute the standard two-point correlations
             do i = 1, size(obs_eq, 1)
                obs_eq(i)%n = obs_eq(i)%n + 1
                obs_eq(i)%ave_sign = obs_eq(i)%ave_sign + 1.d0
             end do
          end if

          Zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call Predefined_Hoppings_Compute_Kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin*dble(N_SUN)
          Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*Z_fac

          ZPot = cmplx(0.d0, 0.d0, kind(0.d0))
          do I = 1, Latt%N
             do no_I = 1, Latt_unit%Norb
                I1 = Invlist(I, no_I)
                ZPot = ZPot - Grc(i1, i1, 1)*Grc(i1, i1, 2)*ham_U
             end do
          end do
          Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Zpot*Z_fac

          Zrho = cmplx(0.d0, 0.d0, kind(0.d0))
          do nf = 1, N_FL
             do I = 1, Ndim
                Zrho = Zrho + Grc(i, i, nf)
             end do
          end do
          Zrho = Zrho*dble(N_SUN)
          Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Zrho*Z_fac

          Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*Z_fac

          nc = 0
          do nf = 1, N_FL
             do I = 1, Ndim
             do J = 1, Ndim
                nc = nc + 1
                ztmp = grc(i, j, nf)
                Obs_scal(5)%Obs_vec(nc) = Obs_scal(5)%Obs_vec(nc) + ztmp*Z_fac
             end do
             end do
          end do

          nc = 0
          do nf = 1, N_FL
             do I = 1, Ndim
             do J = 1, Ndim
                nc = nc + 1
                zone = cmplx(0.d0, 0.d0, kind(0.d0))
                if (I .eq. J) zone = cmplx(1.d0, 0.d0, kind(0.d0))
                Ztmp = zone - GR_mix(J, I, nf)
                Obs_scal(6)%Obs_vec(nc) = Obs_scal(6)%Obs_vec(nc) + ztmp*Z_fac
             end do
             end do
          end do

          ! Standard two-point correlations
          z1j = cmplx(0.d0,1.d0,kind(0.d0))
          do i1 = 1, ndim
              i    = list(i1,1)
              no_i = list(i1,2)

              do j1 = 1, ndim
                  j    = list(j1,1)
                  no_j = list(j1,2)

                  imj  = latt%imj(i, j)
                  ztmp = grc(i1,j1,1) + grc(i1,j1,2)
                  obs_eq(1)%obs_Latt(imj,1,no_i,no_j) = obs_eq(1)%obs_latt(imj,1,no_i,no_j) + ztmp*z_fac

                  ztmp = grc(i1,j1,1)*gr(i1,j1,1) + grc(i1,j1,2)*gr(i1,j1,2) + &
                      & (grc(i1,i1,1) - grc(i1,i1,2))*(grc(j1,j1,1) - grc(j1,j1,2))
                  obs_eq(2)%obs_Latt(imj,1,no_i,no_j) = obs_eq(2)%obs_latt(imj,1,no_i,no_j) + ztmp*z_fac
                  
                  ztmp = grc(i1,j1,1)*gr(i1,j1,1) + grc(i1,j1,2)*gr(i1,j1,2) + &
                      & (grc(i1,i1,2) + grc(i1,i1,1))*(grc(j1,j1,2) + grc(j1,j1,1))
                  obs_eq(3)%obs_Latt(imj,1,no_i,no_j) = obs_eq(3)%obs_latt(imj,1,no_i,no_j) + ztmp*z_fac
                  
                  !! s wave
                  ztmp = grc(i1,j1,1)*grc(i1,j1,2)! + gr(i1,j1,1)*gr(i1,j1,2)
                  obs_eq(4)%obs_Latt(imj,1,no_i,no_j) = obs_eq(4)%obs_latt(imj,1,no_i,no_j) + ztmp*z_fac
                  
                  !! d_{xy}
                  !! \Delta_{xy} = c^{i,\uparrow}c^{i+x+y,\dnarrow}+c^{i,\uparrow}c^{i-x-y,\dnarrow}-
                  !!               c^{i,\uparrow}c^{i+x-y,\dnarrow}-c^{i,\uparrow}c^{i-x+y,\dnarrow}
                  ipxpy = invlist(latt%nnlist(i,1, 1),no_i); imxmy = invlist(latt%nnlist(i,-1,-1),no_i);
                  ipxmy = invlist(latt%nnlist(i,1,-1),no_i); imxpy = invlist(latt%nnlist(i,-1, 1),no_i);
                  jpxpy = invlist(latt%nnlist(j,1, 1),no_j); jmxmy = invlist(latt%nnlist(j,-1,-1),no_j);
                  jpxmy = invlist(latt%nnlist(j,1,-1),no_j); jmxpy = invlist(latt%nnlist(j,-1, 1),no_j);
                  
                  ztmp1 =  grc(ipxpy,jpxpy,2)*grc(i1,j1,1) + grc(imxmy,jpxpy,2)*grc(i1,j1,1) &
                      &  + grc(ipxpy,jmxmy,2)*grc(i1,j1,1) + grc(imxmy,jmxmy,2)*grc(i1,j1,1) &
                      &  - grc(ipxpy,jpxmy,2)*grc(i1,j1,1) - grc(imxmy,jpxmy,2)*grc(i1,j1,1) &
                      &  - grc(ipxpy,jmxpy,2)*grc(i1,j1,1) - grc(imxmy,jmxpy,2)*grc(i1,j1,1) &
                      &  - grc(ipxmy,jpxpy,2)*grc(i1,j1,1) - grc(imxpy,jpxpy,2)*grc(i1,j1,1) &
                      &  - grc(ipxmy,jmxmy,2)*grc(i1,j1,1) - grc(imxpy,jmxmy,2)*grc(i1,j1,1) &
                      &  + grc(ipxmy,jpxmy,2)*grc(i1,j1,1) + grc(imxpy,jpxmy,2)*grc(i1,j1,1) &
                      &  + grc(ipxmy,jmxpy,2)*grc(i1,j1,1) + grc(imxpy,jmxpy,2)*grc(i1,j1,1)
                  
                  ztmp2 =  grc(ipxpy,jpxpy,1)*grc(i1,j1,2) + grc(imxmy,jpxpy,1)*grc(i1,j1,2) &
                      &  + grc(ipxpy,jmxmy,1)*grc(i1,j1,2) + grc(imxmy,jmxmy,1)*grc(i1,j1,2) &
                      &  - grc(ipxpy,jpxmy,1)*grc(i1,j1,2) - grc(imxmy,jpxmy,1)*grc(i1,j1,2) &
                      &  - grc(ipxpy,jmxpy,1)*grc(i1,j1,2) - grc(imxmy,jmxpy,1)*grc(i1,j1,2) &
                      &  - grc(ipxmy,jpxpy,1)*grc(i1,j1,2) - grc(imxpy,jpxpy,1)*grc(i1,j1,2) &
                      &  - grc(ipxmy,jmxmy,1)*grc(i1,j1,2) - grc(imxpy,jmxmy,1)*grc(i1,j1,2) &
                      &  + grc(ipxmy,jpxmy,1)*grc(i1,j1,2) + grc(imxpy,jpxmy,1)*grc(i1,j1,2) &
                      &  + grc(ipxmy,jmxpy,1)*grc(i1,j1,2) + grc(imxpy,jmxpy,1)*grc(i1,j1,2)
                  
                  ztmp3 =  grc(ipxpy,j1,2)*grc(i1,jpxpy,1) + grc(imxmy,j1,2)*grc(i1,jpxpy,1) &
                      &  + grc(ipxpy,j1,2)*grc(i1,jmxmy,1) + grc(imxmy,j1,2)*grc(i1,jmxmy,1) &
                      &  - grc(ipxpy,j1,2)*grc(i1,jpxmy,1) - grc(imxmy,j1,2)*grc(i1,jpxmy,1) &
                      &  - grc(ipxpy,j1,2)*grc(i1,jmxpy,1) - grc(imxmy,j1,2)*grc(i1,jmxpy,1) &
                      &  - grc(ipxmy,j1,2)*grc(i1,jpxpy,1) - grc(imxpy,j1,2)*grc(i1,jpxpy,1) &
                      &  - grc(ipxmy,j1,2)*grc(i1,jmxmy,1) - grc(imxpy,j1,2)*grc(i1,jmxmy,1) &
                      &  + grc(ipxmy,j1,2)*grc(i1,jpxmy,1) + grc(imxpy,j1,2)*grc(i1,jpxmy,1) &
                      &  + grc(ipxmy,j1,2)*grc(i1,jmxpy,1) + grc(imxpy,j1,2)*grc(i1,jmxpy,1)
                  
                  ztmp4 =  grc(ipxpy,j1,1)*grc(i1,jpxpy,2) + grc(imxmy,j1,1)*grc(i1,jpxpy,2) &
                      &  + grc(ipxpy,j1,1)*grc(i1,jmxmy,2) + grc(imxmy,j1,1)*grc(i1,jmxmy,2) &
                      &  - grc(ipxpy,j1,1)*grc(i1,jpxmy,2) - grc(imxmy,j1,1)*grc(i1,jpxmy,2) &
                      &  - grc(ipxpy,j1,1)*grc(i1,jmxpy,2) - grc(imxmy,j1,1)*grc(i1,jmxpy,2) &
                      &  - grc(ipxmy,j1,1)*grc(i1,jpxpy,2) - grc(imxpy,j1,1)*grc(i1,jpxpy,2) &
                      &  - grc(ipxmy,j1,1)*grc(i1,jmxmy,2) - grc(imxpy,j1,1)*grc(i1,jmxmy,2) &
                      &  + grc(ipxmy,j1,1)*grc(i1,jpxmy,2) + grc(imxpy,j1,1)*grc(i1,jpxmy,2) &
                      &  + grc(ipxmy,j1,1)*grc(i1,jmxpy,2) + grc(imxpy,j1,1)*grc(i1,jmxpy,2)

                  obs_eq(5)%obs_Latt(imj,1,no_i,no_j) = obs_eq(5)%obs_latt(imj,1,no_i,no_j) &
                      & + (ztmp1+ztmp2+ztmp3+ztmp4)*(1.d0/16.d0)*z_fac
                  
                  !! d_{x^2-y^2}
                  !! \Delta_{x^2-y^2} = c^{i,\uparrow}c^{i\pm x,\dnarrow}-c^{i,\uparrow}c^{i\pm y,\dnarrow}
                  ipx = invlist(latt%nnlist(i, 1,0),no_i); ipy = invlist(latt%nnlist(i,0, 1),no_i)
                  imx = invlist(latt%nnlist(i,-1,0),no_i); imy = invlist(latt%nnlist(i,0,-1),no_i)
                  jpx = invlist(latt%nnlist(j, 1,0),no_j); jpy = invlist(latt%nnlist(j,0, 1),no_j)
                  jmx = invlist(latt%nnlist(j,-1,0),no_j); jmy = invlist(latt%nnlist(j,0,-1),no_j)

                  ztmp1 =  grc(ipx,jpx,2)*grc(i1,j1,1) + grc(imx,jpx,2)*grc(i1,j1,1) &
                      &  + grc(ipx,jmx,2)*grc(i1,j1,1) + grc(imx,jmx,2)*grc(i1,j1,1) &
                      &  - grc(ipx,jpy,2)*grc(i1,j1,1) - grc(imx,jpy,2)*grc(i1,j1,1) &
                      &  - grc(ipx,jmy,2)*grc(i1,j1,1) - grc(imx,jmy,2)*grc(i1,j1,1) &
                      &  - grc(ipy,jpx,2)*grc(i1,j1,1) - grc(imy,jpx,2)*grc(i1,j1,1) &
                      &  - grc(ipy,jmx,2)*grc(i1,j1,1) - grc(imy,jmx,2)*grc(i1,j1,1) &
                      &  + grc(ipy,jpy,2)*grc(i1,j1,1) + grc(imy,jpy,2)*grc(i1,j1,1) &
                      &  + grc(ipy,jmy,2)*grc(i1,j1,1) + grc(imy,jmy,2)*grc(i1,j1,1)
                  
                  ztmp2 =  grc(ipx,jpx,1)*grc(i1,j1,2) + grc(imx,jpx,1)*grc(i1,j1,2) &
                      &  + grc(ipx,jmx,1)*grc(i1,j1,2) + grc(imx,jmx,1)*grc(i1,j1,2) &
                      &  - grc(ipx,jpy,1)*grc(i1,j1,2) - grc(imx,jpy,1)*grc(i1,j1,2) &
                      &  - grc(ipx,jmy,1)*grc(i1,j1,2) - grc(imx,jmy,1)*grc(i1,j1,2) &
                      &  - grc(ipy,jpx,1)*grc(i1,j1,2) - grc(imy,jpx,1)*grc(i1,j1,2) &
                      &  - grc(ipy,jmx,1)*grc(i1,j1,2) - grc(imy,jmx,1)*grc(i1,j1,2) &
                      &  + grc(ipy,jpy,1)*grc(i1,j1,2) + grc(imy,jpy,1)*grc(i1,j1,2) &
                      &  + grc(ipy,jmy,1)*grc(i1,j1,2) + grc(imy,jmy,1)*grc(i1,j1,2)
                  
                  ztmp3 =  grc(ipx,j1,2)*grc(i1,jpx,1) + grc(imx,j1,2)*grc(i1,jpx,1) &
                      &  + grc(ipx,j1,2)*grc(i1,jmx,1) + grc(imx,j1,2)*grc(i1,jmx,1) &
                      &  - grc(ipx,j1,2)*grc(i1,jpy,1) - grc(imx,j1,2)*grc(i1,jpy,1) &
                      &  - grc(ipx,j1,2)*grc(i1,jmy,1) - grc(imx,j1,2)*grc(i1,jmy,1) &
                      &  - grc(ipy,j1,2)*grc(i1,jpx,1) - grc(imy,j1,2)*grc(i1,jpx,1) &
                      &  - grc(ipy,j1,2)*grc(i1,jmx,1) - grc(imy,j1,2)*grc(i1,jmx,1) &
                      &  + grc(ipy,j1,2)*grc(i1,jpy,1) + grc(imy,j1,2)*grc(i1,jpy,1) &
                      &  + grc(ipy,j1,2)*grc(i1,jmy,1) + grc(imy,j1,2)*grc(i1,jmy,1)
                  
                  ztmp4 =  grc(ipx,j1,1)*grc(i1,jpx,2) + grc(imx,j1,1)*grc(i1,jpx,2) &
                      &  + grc(ipx,j1,1)*grc(i1,jmx,2) + grc(imx,j1,1)*grc(i1,jmx,2) &
                      &  - grc(ipx,j1,1)*grc(i1,jpy,2) - grc(imx,j1,1)*grc(i1,jpy,2) &
                      &  - grc(ipx,j1,1)*grc(i1,jmy,2) - grc(imx,j1,1)*grc(i1,jmy,2) &
                      &  - grc(ipy,j1,1)*grc(i1,jpx,2) - grc(imy,j1,1)*grc(i1,jpx,2) &
                      &  - grc(ipy,j1,1)*grc(i1,jmx,2) - grc(imy,j1,1)*grc(i1,jmx,2) &
                      &  + grc(ipy,j1,1)*grc(i1,jpy,2) + grc(imy,j1,1)*grc(i1,jpy,2) &
                      &  + grc(ipy,j1,1)*grc(i1,jmy,2) + grc(imy,j1,1)*grc(i1,jmy,2)

                  obs_eq(6)%obs_Latt(imj,1,no_i,no_j) = obs_eq(6)%obs_latt(imj,1,no_i,no_j) & 
                      & + (ztmp1+ztmp2+ztmp3+ztmp4)*(1.d0/16.d0)*z_fac
                  
                  !! d_{x^2+y^2}
                  !! \Delta_{x^2+y^2} = c^{i,\uparrow}c^{i\pm x,\dnarrow}+c^{i,\uparrow}c^{i\pm y,\dnarrow}
                  ipx = invlist(latt%nnlist(i, 1,0),no_i); ipy = invlist(latt%nnlist(i,0, 1),no_i)
                  imx = invlist(latt%nnlist(i,-1,0),no_i); imy = invlist(latt%nnlist(i,0,-1),no_i)
                  jpx = invlist(latt%nnlist(j, 1,0),no_j); jpy = invlist(latt%nnlist(j,0, 1),no_j)
                  jmx = invlist(latt%nnlist(j,-1,0),no_j); jmy = invlist(latt%nnlist(j,0,-1),no_j)

                  ztmp1 =  grc(ipx,jpx,2)*grc(i1,j1,1) + grc(imx,jpx,2)*grc(i1,j1,1) &
                      &  + grc(ipx,jmx,2)*grc(i1,j1,1) + grc(imx,jmx,2)*grc(i1,j1,1) &
                      &  + grc(ipx,jpy,2)*grc(i1,j1,1) + grc(imx,jpy,2)*grc(i1,j1,1) &
                      &  + grc(ipx,jmy,2)*grc(i1,j1,1) + grc(imx,jmy,2)*grc(i1,j1,1) &
                      &  + grc(ipy,jpx,2)*grc(i1,j1,1) + grc(imy,jpx,2)*grc(i1,j1,1) &
                      &  + grc(ipy,jmx,2)*grc(i1,j1,1) + grc(imy,jmx,2)*grc(i1,j1,1) &
                      &  + grc(ipy,jpy,2)*grc(i1,j1,1) + grc(imy,jpy,2)*grc(i1,j1,1) &
                      &  + grc(ipy,jmy,2)*grc(i1,j1,1) + grc(imy,jmy,2)*grc(i1,j1,1)
                  
                  ztmp2 =  grc(ipx,jpx,1)*grc(i1,j1,2) + grc(imx,jpx,1)*grc(i1,j1,2) &
                      &  + grc(ipx,jmx,1)*grc(i1,j1,2) + grc(imx,jmx,1)*grc(i1,j1,2) &
                      &  + grc(ipx,jpy,1)*grc(i1,j1,2) + grc(imx,jpy,1)*grc(i1,j1,2) &
                      &  + grc(ipx,jmy,1)*grc(i1,j1,2) + grc(imx,jmy,1)*grc(i1,j1,2) &
                      &  + grc(ipy,jpx,1)*grc(i1,j1,2) + grc(imy,jpx,1)*grc(i1,j1,2) &
                      &  + grc(ipy,jmx,1)*grc(i1,j1,2) + grc(imy,jmx,1)*grc(i1,j1,2) &
                      &  + grc(ipy,jpy,1)*grc(i1,j1,2) + grc(imy,jpy,1)*grc(i1,j1,2) &
                      &  + grc(ipy,jmy,1)*grc(i1,j1,2) + grc(imy,jmy,1)*grc(i1,j1,2)
                  
                  ztmp3 = - grc(ipx,j1,2)*grc(i1,jpx,1) - grc(imx,j1,2)*grc(i1,jpx,1) &
                      &   - grc(ipx,j1,2)*grc(i1,jmx,1) - grc(imx,j1,2)*grc(i1,jmx,1) &
                      &   - grc(ipx,j1,2)*grc(i1,jpy,1) - grc(imx,j1,2)*grc(i1,jpy,1) &
                      &   - grc(ipx,j1,2)*grc(i1,jmy,1) - grc(imx,j1,2)*grc(i1,jmy,1) &
                      &   - grc(ipy,j1,2)*grc(i1,jpx,1) - grc(imy,j1,2)*grc(i1,jpx,1) &
                      &   - grc(ipy,j1,2)*grc(i1,jmx,1) - grc(imy,j1,2)*grc(i1,jmx,1) &
                      &   - grc(ipy,j1,2)*grc(i1,jpy,1) - grc(imy,j1,2)*grc(i1,jpy,1) &
                      &   - grc(ipy,j1,2)*grc(i1,jmy,1) - grc(imy,j1,2)*grc(i1,jmy,1)
                  
                  ztmp4 = - grc(ipx,j1,1)*grc(i1,jpx,2) - grc(imx,j1,1)*grc(i1,jpx,2) &
                      &   - grc(ipx,j1,1)*grc(i1,jmx,2) - grc(imx,j1,1)*grc(i1,jmx,2) &
                      &   - grc(ipx,j1,1)*grc(i1,jpy,2) - grc(imx,j1,1)*grc(i1,jpy,2) &
                      &   - grc(ipx,j1,1)*grc(i1,jmy,2) - grc(imx,j1,1)*grc(i1,jmy,2) &
                      &   - grc(ipy,j1,1)*grc(i1,jpx,2) - grc(imy,j1,1)*grc(i1,jpx,2) &
                      &   - grc(ipy,j1,1)*grc(i1,jmx,2) - grc(imy,j1,1)*grc(i1,jmx,2) &
                      &   - grc(ipy,j1,1)*grc(i1,jpy,2) - grc(imy,j1,1)*grc(i1,jpy,2) &
                      &   - grc(ipy,j1,1)*grc(i1,jmy,2) - grc(imy,j1,1)*grc(i1,jmy,2)

                  obs_eq(7)%obs_Latt(imj,1,no_i,no_j) = & 
                      & obs_eq(7)%obs_latt(imj,1,no_i,no_j) + (ztmp1+ztmp2+ztmp3+ztmp4)*(1.d0/16.d0)*z_fac
                  
                  !! p_x
                  !! \Delta_{x} = c^{i,\uparrow}c^{i+x,\uparrow}-c^{i,\uparrow}c^{i-x,\uparrow}+
                  !!              c^{i,\dnarrow}c^{i+x,\dnarrow}-c^{i,\dnarrow}c^{i-x,\dnarrow}
                  ztmp1 =   grc(ipx,jpx,1)*grc(i1,j1,1) - grc(ipx,jmx,1)*grc(i1,j1,1) &
                      &   + grc(ipx,jpx,2)*grc(i1,j1,2) - grc(ipx,jmx,2)*grc(i1,j1,2) &
                      &   - grc(imx,jpx,1)*grc(i1,j1,1) + grc(imx,jmx,1)*grc(i1,j1,1) &
                      &   - grc(imx,jpx,2)*grc(i1,j1,2) + grc(imx,jmx,2)*grc(i1,j1,2)
                  
                  !! p_y
                  !! \Delta_{y} = c^{i,\uparrow}c^{i+y,\uparrow}-c^{i,\uparrow}c^{i-y,\uparrow}+
                  !!              c^{i,\dnarrow}c^{i+y,\dnarrow}-c^{i,\dnarrow}c^{i-y,\dnarrow}
                  ztmp2 =   grc(ipy,jpy,1)*grc(i1,j1,1) - grc(ipy,jmy,1)*grc(i1,j1,1) &
                      &   + grc(ipy,jpy,2)*grc(i1,j1,2) - grc(ipy,jmy,2)*grc(i1,j1,2) &
                      &   - grc(imy,jpy,1)*grc(i1,j1,1) + grc(imy,jmy,1)*grc(i1,j1,1) &
                      &   - grc(imy,jpy,2)*grc(i1,j1,2) + grc(imy,jmy,2)*grc(i1,j1,2)

                  obs_eq(8)%obs_Latt(imj,1,no_i,no_j) = &
                      & obs_eq(8)%obs_latt(imj,1,no_i,no_j) + (ztmp1+ztmp2)*(1.d0/16.d0)*z_fac
                  
                  !! p_x+ip_y
                  !! \Delta_{x+iy} = i\Delta_x+\Delta_y
                  ztmp  = - z1j*grc(ipx,jpy,1)*grc(i1,j1,1) + z1j*grc(ipx,jmy,1)*grc(i1,j1,1) &
                      &   - z1j*grc(ipx,jpy,2)*grc(i1,j1,2) + z1j*grc(ipx,jmy,2)*grc(i1,j1,2) &
                      &   + z1j*grc(imx,jpy,1)*grc(i1,j1,1) - z1j*grc(imx,jmy,1)*grc(i1,j1,1) &
                      &   + z1j*grc(imx,jpy,2)*grc(i1,j1,2) - z1j*grc(imx,jmy,2)*grc(i1,j1,2) &
                      &   + z1j*grc(ipy,jpx,1)*grc(i1,j1,1) - z1j*grc(ipy,jmx,1)*grc(i1,j1,1) &
                      &   + z1j*grc(ipy,jpx,2)*grc(i1,j1,2) - z1j*grc(ipy,jmx,2)*grc(i1,j1,2) &
                      &   - z1j*grc(imy,jpx,1)*grc(i1,j1,1) + z1j*grc(imy,jmx,1)*grc(i1,j1,1) &
                      &   - z1j*grc(imy,jpx,2)*grc(i1,j1,2) + z1j*grc(imy,jmx,2)*grc(i1,j1,2)
                  
                  obs_eq(9)%obs_Latt(imj,1,no_i,no_j) = & 
                      & obs_eq(9)%obs_latt(imj,1,no_i,no_j) + (ztmp1+ztmp2+ztmp)*(1.d0/16.d0)*z_fac
                  
              end do
              zback = grc(i1, i1, 2) - grc(i1, i1, 1)
              obs_eq(2)%obs_latt0(no_i) = obs_eq(2)%obs_Latt0(no_i) + zback*z_fac
              zback = grc(i1,i1,1) + grc(i1,i1,2) 
              obs_eq(3)%obs_latt0(no_i) = obs_eq(3)%obs_Latt0(no_i) + zback*z_fac
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
       subroutine obsert(nt, gt0, g0t, g00, gtt, i_wlk, i_grc, sum_w, sum_o, act_mea)

          implicit none

          integer, intent(IN) :: nt
          complex(Kind=kind(0.d0)), intent(IN) :: gt0(ndim, ndim, n_fl), g0t(ndim, ndim, n_fl)
          complex(Kind=kind(0.d0)), intent(IN) :: g00(ndim, ndim, n_fl), gtt(ndim, ndim, n_fl)
          complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
          integer, intent(IN) :: i_wlk, i_grc, act_mea

          !Locals
          complex(Kind=kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, zone, zback, z_ol, zw, z_fac, zero
          real(Kind=kind(0.d0)) :: X
          integer :: IMJ, I, J, k, l, m, n, i1, ip1, ipx, ipy, imx, imy, j1, jpx, jpy, jmx, jmy, no_I, no_J, nc

          Z_ol = exp(overlap(i_grc))/sum_o
          ZW = cmplx(weight_k(i_wlk), 0.d0, kind(0.d0))/sum_w
          Z_fac = Z_ol*ZW

          zone = cmplx(1.d0, 0.d0, kind(0.d0))

          ! Standard two-point correlations
          if ((act_mea .eq. 0) .and. (nt .eq. 0)) then
             do i = 1, size(obs_tau, 1)
                obs_tau(i)%n = obs_tau(i)%n + 1
                obs_tau(i)%ave_sign = obs_tau(i)%ave_sign + 1.d0
             end do
          end if
          
          do i1 = 1, ndim
              i    = list(i1,1)
              no_i = list(i1,2)

              do j1 = 1, ndim
                  j    = list(j1,1)
                  no_j = list(j1,2)

                  imj  = latt%imj(i, j)
                  z = gt0(i1,j1,1) + gt0(i1,j1,2)
                  obs_tau(1)%obs_Latt(imj,nt+1,no_i,no_j) = obs_tau(1)%obs_latt(imj,nt+1,no_i,no_j) + z*z_fac

                  z = -g0t(j1,i1,1)*gt0(i1,j1,1) - g0t(j1,i1,2)*gt0(i1,j1,2) + &
                      & (gtt(i1,i1,1) - gtt(i1,i1,2))*(g00(j1,j1,1) - g00(j1,j1,2))
                  obs_tau(2)%obs_Latt(imj,nt+1,no_i,no_j) = obs_tau(2)%obs_latt(imj,nt+1,no_i,no_j) + z*z_fac
                  
                  z = -g0t(j1,i1,1)*gt0(i1,j1,1) - g0t(j1,i1,2)*gt0(i1,j1,2) + &
                      & (zone-gtt(i1,i1,1)+zone-gtt(i1,i1,2))*(zone-gtt(j1,j1,2)+zone-gtt(j1,j1,1))
                  obs_tau(3)%obs_Latt(imj,nt+1,no_i,no_j) = obs_tau(3)%obs_latt(imj,nt+1,no_i,no_j) + z*z_fac
                  
                  z = g0t(j1,i1,1)*g0t(j1,i1,2)! + gt0(i1,j1,1)*gt0(i1,j1,2)
                  obs_tau(4)%obs_Latt(imj,nt+1,no_i,no_j) = obs_tau(4)%obs_latt(imj,nt+1,no_i,no_j) + z*z_fac
              end do
              zback = gtt(i1, i1, 2) - gtt(i1, i1, 1)
              obs_tau(2)%obs_latt0(no_i) = obs_tau(2)%obs_Latt0(no_i) + zback*z_fac
              zback = zone - gtt(i1,i1,1) + zone - gtt(i1,i1,2) 
              obs_tau(3)%obs_latt0(no_i) = obs_tau(3)%obs_Latt0(no_i) + zback*z_fac
          end do

       end subroutine obsert
       
       complex(Kind=kind(0.d0)) function E0_local(gr)
          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: gr(ndim, ndim, n_fl)

          !Local
          complex(Kind=kind(0.d0)) :: grc(ndim, ndim, n_fl), ZK
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY
          integer :: I, J, imj, nf, dec, I1, J1, no_I, no_J, n
          real(Kind=kind(0.d0)) :: X

          do nf = 1, N_FL
             do I = 1, Ndim
                do J = 1, Ndim
                   grc(I, J, nf) = -gr(J, I, nf)
                end do
                grc(I, I, nf) = 1.d0 + grc(I, I, nf)
             end do
          end do

          Zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call Predefined_Hoppings_Compute_Kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin*dble(N_SUN)

          ZPot = cmplx(0.d0, 0.d0, kind(0.d0))
          do I = 1, Latt%N
             do no_I = 1, Latt_unit%Norb
                I1 = Invlist(I, no_I)
                ZPot = ZPot + 2.d0*Grc(i1, i1, 1)*Grc(i1, i1, 2) - &
                    & Grc(i1, i1, 1) - Grc(i1, i1, 2) + 1.d0
             end do
          end do
          zpot = zpot*(-ham_u/2.d0)

          E0_local = (Zpot + ZKin)*dtau

       end function E0_local

       subroutine sum_weight(z_sum_weight)
          implicit none

          complex(Kind=kind(0.d0)), intent(out) :: z_sum_weight

          !local
          integer                   :: i_wlk
          real(Kind=kind(0.d0)) :: X1, tot_re_weight
          complex(Kind=kind(0.d0)) :: Z1, Z2, ZP, wtmp

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif

          Z1 = cmplx(0.d0, 0.d0, kind(0.d0))
          Z2 = cmplx(0.d0, 0.d0, kind(0.d0))
          do i_wlk = 1, N_wlk
             ZP = 1.d0
             wtmp = cmplx(weight_k(i_wlk), 0.d0, kind(0.d0))
             Z1 = Z1 + wtmp*ZP
          end do
          call MPI_REDUCE(Z1, Z2, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, IERR)
          call MPI_BCAST(Z2, 1, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr)

          z_sum_weight = Z2

       end subroutine sum_weight

       subroutine update_fac_norm(GR, ntw)
          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL, N_grc)
          integer, intent(IN) :: ntw

          !local
          integer                   :: i_wlk, ns, i_grc
          real(Kind=kind(0.d0)) :: X1, tot_re_weight
          complex(Kind=kind(0.d0)) :: tot_ene, Z1, Z2, tot_c_weight, ZP, wtmp, el_tmp, Z
          character(LEN=64)  :: filename

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

         !! initial energy
          tot_ene = 0.d0
          tot_c_weight = 0.d0
          do i_wlk = 1, N_wlk

             if (weight_k(i_wlk) .ge. 0.d0) then

                Z = 0.d0
                do ns = 1, N_slat
                   i_grc = ns + (i_wlk - 1)*N_slat
                   Z = Z + exp(overlap(i_grc))
                end do

                do ns = 1, N_slat
                   i_grc = ns + (i_wlk - 1)*N_slat
                   el_tmp = dble(ham%E0_local(GR(:, :, :, i_grc)))
                   tot_ene = tot_ene + el_tmp*weight_k(i_wlk)*exp(overlap(i_grc))/Z
                end do

                tot_c_weight = tot_c_weight + weight_k(i_wlk)

             end if

          end do
          tot_re_weight = dble(tot_c_weight)

          call MPI_REDUCE(tot_ene, Z1, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, IERR)
          call MPI_REDUCE(tot_re_weight, X1, 1, MPI_REAL8, MPI_SUM, 0, Group_comm, IERR)

          if (Irank_g == 0) then
             Z1 = Z1/X1
             fac_norm = dble(Z1)
             open (UNIT=77, FILE=filename, STATUS='UNKNOWN', position="append")
             write (77, *) ntw*dtau, dble(Z1)/dtau
             close (77)
          end if
          call MPI_BCAST(fac_norm, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

       end subroutine update_fac_norm

    end submodule ham_bose_metal_smod
