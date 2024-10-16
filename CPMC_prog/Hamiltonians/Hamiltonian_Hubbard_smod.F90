    submodule(Hamiltonian_main) ham_Hubbard_smod

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

       type, extends(ham_base) :: ham_Hubbard
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: Alloc_obs
          procedure, nopass :: Obser
          procedure, nopass :: ObserT
          procedure, nopass :: S0
          procedure, nopass :: E0_local
          procedure, nopass :: sum_weight
          procedure, nopass :: update_fac_norm
          procedure, nopass :: bp_obsert
          procedure, nopass :: obsert_mc
          procedure, nopass :: init_obs_mc
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_Hubbard

       !#PARAMETERS START# VAR_lattice
       character(len=64) :: Model = 'Hubbard'  ! Value not relevant
       character(len=64) :: Lattice_type = 'Square'
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
       logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
       !logical              :: Symm         = .true.   ! Whether symmetrization takes place
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_Hubbard
       real(Kind=kind(0.d0)) :: ham_T = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: Ham_chem = 0.d0     ! Chemical potential
       real(Kind=kind(0.d0)) :: Ham_U = 4.d0     ! Hubbard interaction
       real(Kind=kind(0.d0)) :: ham_T2 = 1.d0     ! For bilayer systems
       real(Kind=kind(0.d0)) :: Ham_U2 = 4.d0     ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_Tperp = 1.d0     ! For bilayer systems
       integer               :: N_dope = 0        ! Number of doping electrons
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell

      !! mc obser
       type(obser_Latt), dimension(:), allocatable :: obs_grc
       type(obser_Latt), dimension(:), allocatable :: obs_spinz
       type(obser_Latt), dimension(:), allocatable :: obs_spinxy
       type(obser_Latt), dimension(:), allocatable :: obs_spint

    contains

       module subroutine Ham_Alloc_hubbard
          allocate (ham_Hubbard::ham)
       end subroutine Ham_Alloc_hubbard

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Hubbard_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
       subroutine Ham_Set

#if defined (MPI)
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

          ! From dynamically generated file "Hamiltonian_Hubbard_read_write_parameters.F90"
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
          call Ham_Trial()

#ifdef MPI
          if (Irank_g == 0) then
#endif
             File_info = "info"
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
             write (unit_info, *) 'HS  couples to z-component of spin'
             write (unit_info, *) 'Checkerboard  : ', Checkerboard
             write (unit_info, *) 'Symm. decomp  : ', Symm
             write (unit_info, *) 'dtau: ', dtau
             write (unit_info, *) 'ltrot         : ', ltrot
             write (unit_info, *) 'N_SUN         : ', N_SUN
             write (unit_info, *) 'N_FL          : ', N_FL
             write (unit_info, *) 'N_wlk         : ', N_wlk
             write (unit_info, *) 'N_wlk_mpi     : ', N_wlk_mpi
             write (unit_info, *) 'N_slat        : ', N_slat
             write (unit_info, *) 'N_grc         : ', N_grc
             write (unit_info, *) 'N_grc_mpi     : ', N_grc_mpi
             write (unit_info, *) 'N_dope        : ', N_dope
             write (unit_info, *) 't             : ', Ham_T
             write (unit_info, *) 'Ham_U         : ', Ham_U
             write (unit_info, *) 't2            : ', Ham_T2
             write (unit_info, *) 'Ham_U2        : ', Ham_U2
             write (unit_info, *) 'Ham_tperp     : ', Ham_tperp
             write (unit_info, *) 'Ham_chem      : ', Ham_chem
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

          real(Kind=kind(0.d0)), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                  Ham_T2_vec(:), Ham_Lambda_vec(:)
          integer, allocatable ::   N_Phi_vec(:)

          ! Use predefined stuctures or set your own hopping
          integer :: n, nth

          allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &                                   N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL))

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_T_vec = Ham_T
          Ham_Tperp_vec = Ham_Tperp
          Ham_Chem_vec = Ham_Chem
          Phi_X_vec = Phi_X
          Phi_Y_vec = Phi_Y
          Ham_T2_vec = Ham_T2
          Ham_Lambda_vec = Ham_Lambda
          N_Phi_vec = N_Phi

          select case (Lattice_type)
          case ("Square")
             call Set_Default_hopping_parameters_square(Hopping_Matrix, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          case ("Honeycomb")
             Ham_Lambda = 0.d0
      call Set_Default_hopping_parameters_honeycomb(Hopping_Matrix, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                                &                                         Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
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
          ! Use predefined stuctures or set your own Trial  wave function
          N_part = Ndim/2 - N_dope

          call Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
               &                            N_part, N_FL, N_slat, wf_l, wf_r)

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
          real(Kind=kind(0.d0)) :: X, Zero = 1.d-10
          real(Kind=kind(0.d0)), allocatable :: Ham_U_vec(:)

          allocate (Ham_U_vec(Latt_unit%Norb))

          N_ops = 0
          Ham_U_vec(:) = Ham_U
          if (abs(Ham_U) > Zero) N_ops = N_ops + Latt%N*Latt_unit%Norb

          allocate (Op_V(N_ops, N_FL))
          Ham_U_vec = Ham_U_vec/real(N_SUN, kind(0.d0))
          nc = 0
          do I1 = 1, Latt%N
             do no = 1, Latt_unit%Norb
                I = invlist(I1, no)
                if (abs(Ham_U_vec(no)) > Zero) then
                   nc = nc + 1
                   call Predefined_Int_U_MZ(OP_V(nc, 1), OP_V(nc, 2), I, DTAU, Ham_U_vec(no))
                end if
             end do
          end do

          deallocate (Ham_U_vec)

       end subroutine Ham_V

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
       subroutine Alloc_obs(ltau, lmetropolis)

          implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          integer, intent(In) :: ltau, lmetropolis
          integer    ::  i, N, Nt
          character(len=64) ::  Filename
          character(len=2)  ::  Channel

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

          allocate (Obs_eq(5))
          do I = 1, size(Obs_eq, 1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "SpinXY"
             case (4)
                Filename = "SpinT"
             case (5)
                Filename = "Den"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          end do

          if (Ltau == 1) then
             ! Time-displaced correlators
             allocate (Obs_tau(5))
             do I = 1, size(Obs_tau, 1)
                select case (I)
                case (1)
                   Channel = 'P'; Filename = "Green"
                case (2)
                   Channel = 'PH'; Filename = "SpinZ"
                case (3)
                   Channel = 'PH'; Filename = "SpinXY"
                case (4)
                   Channel = 'PH'; Filename = "SpinT"
                case (5)
                   Channel = 'PH'; Filename = "Den"
                case default
                   write (6, *) ' Error in Alloc_obs '
                end select
                Nt = Ltrot + 1
                Channel = 'T0'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             end do

             if (lmetropolis) then
                Nt = Ltrot + 1
                Channel = 'T0'
                Filename = "mc_obs"
                allocate (obs_grc(n_wlk))
                allocate (obs_spinz(n_wlk))
                allocate (obs_spinxy(n_wlk))
                allocate (obs_spint(n_wlk))
                do i = 1, n_wlk
                   call obser_latt_make(obs_grc(i), nt, Filename, Latt, Latt_unit, channel, dtau)
                   call obser_latt_make(obs_spinz(i), nt, Filename, Latt, Latt_unit, channel, dtau)
                   call obser_latt_make(obs_spinxy(i), nt, Filename, Latt, Latt_unit, channel, dtau)
                   call obser_latt_make(obs_spint(i), nt, Filename, Latt, Latt_unit, channel, dtau)
                end do
             end if

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
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
       subroutine Obser(GR, GR_mix, i_wlk, i_grc, sum_w, sum_o, act_mea)

          use Predefined_Obs

          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: GR_mix(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
          integer, intent(IN) :: i_wlk, i_grc, act_mea

          !Local
          complex(Kind=kind(0.d0)) :: GRC(Ndim, Ndim, N_FL), ZK, invsumw, zone, ztmp, z_ol
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY, ZW, Re_ZW, Z_fac
          integer :: I, J, imj, nf, dec, I1, J1, no_I, no_J, n, nc
          real(Kind=kind(0.d0)) :: X

          Z_ol = exp(overlap(i_grc))/sum_o
          ZW = cmplx(weight_k(i_wlk), 0.d0, kind(0.d0))/sum_w
          Z_fac = Z_ol*ZW

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
          if (act_mea .eq. 0) then
             do I = 1, size(Obs_scal, 1)
                Obs_scal(I)%N = Obs_scal(I)%N + 1
                Obs_scal(I)%Ave_sign = Obs_scal(I)%Ave_sign + 1.d0
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
                ZPot = ZPot + Grc(i1, i1, 1)*Grc(i1, i1, 2)*ham_U
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

          ! Compute the standard two-point correlations
          if (act_mea .eq. 0) then
             do I = 1, size(Obs_eq, 1)
                obs_eq(I)%N = obs_eq(I)%N + 1
                obs_eq(I)%Ave_sign = obs_eq(I)%Ave_sign + 1.d0
             end do
          end if

          ! Standard two-point correlations
          call Predefined_Obs_eq_Green_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, Z_fac, Obs_eq(1))
          call Predefined_Obs_eq_SpinMz_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, Z_fac, Obs_eq(2), Obs_eq(3), Obs_eq(4))
          call Predefined_Obs_eq_Den_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, Z_fac, Obs_eq(5))

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
!-------------------------------------------------------------------
       subroutine ObserT(NT, GT0, G0T, G00, GTT, i_wlk, i_grc, sum_w, sum_o, act_mea)

          use Predefined_Obs

          implicit none

          integer, intent(IN) :: NT
  complex(Kind=kind(0.d0)), intent(IN) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL), G00(Ndim, Ndim, N_FL), GTT(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
          integer, intent(IN) :: i_wlk, i_grc, act_mea

          !Locals
          complex(Kind=kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, ZW, Re_ZW, Z_fac, Z_ol
          real(Kind=kind(0.d0)) :: X
          integer :: IMJ, I, J, I1, J1, no_I, no_J

          Z_ol = exp(overlap(i_grc))/sum_o
          ZW = cmplx(weight_k(i_wlk), 0.d0, kind(0.d0))/sum_w
          Z_fac = Z_ol*ZW

          ! Compute the dynamic two-point correlations
          if ((act_mea .eq. 0) .and. (nt .eq. 0)) then
             do I = 1, size(Obs_tau, 1)
                obs_tau(I)%N = obs_tau(I)%N + 1
                obs_tau(I)%Ave_sign = obs_tau(I)%Ave_sign + 1.d0
             end do
          end if

          ! Standard two-point correlations

          call Predefined_Obs_tau_Green_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, Z_fac, Obs_tau(1))
          call Predefined_Obs_tau_SpinMz_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, Z_fac, Obs_tau(2),&
               &                                   Obs_tau(3), Obs_tau(4))
          call Predefined_Obs_tau_Den_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, Z_fac, Obs_tau(5))

       end subroutine OBSERT

       subroutine bp_obsert(i_wlk, i_grc, sum_w, act_mea)
          implicit none

          integer, intent(in) :: i_wlk, i_grc, act_mea
          complex(Kind=kind(0.d0)), intent(IN) :: sum_w

          !Locals
          complex(Kind=kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, ZW, Re_ZW, Z_fac, Z_ol
          real(Kind=kind(0.d0)) :: X
          integer :: IMJ, I, J, I1, J1, no_I, no_J, nt

          ZW = cmplx(weight_k(i_wlk), 0.d0, kind(0.d0))/sum_w
          Z_fac = ZW

          ! Compute the dynamic two-point correlations
          if (act_mea .eq. 0) then
             do I = 1, size(Obs_tau, 1)
                obs_tau(I)%N = obs_tau(I)%N + 1
                obs_tau(I)%Ave_sign = obs_tau(I)%Ave_sign + 1.d0
             end do
          end if

          obs_grc(i_wlk)%ave_sign = obs_grc(i_wlk)%ave_sign/dble(obs_grc(i_wlk)%N)
          obs_grc(i_wlk)%obs_latt(:, :, :, :) = &
              & obs_grc(i_wlk)%obs_latt(:, :, :, :)/dble(obs_grc(i_wlk)%N)/obs_grc(i_wlk)%ave_sign

          obs_tau(1)%obs_latt(:, :, :, :) = obs_tau(1)%obs_latt(:, :, :, :) + &
              & obs_grc(i_wlk)%obs_latt(:, :, :, :)*z_fac

          obs_spinz(i_wlk)%ave_sign = obs_spinz(i_wlk)%ave_sign/dble(obs_spinz(i_wlk)%N)
          obs_spinz(i_wlk)%obs_latt(:, :, :, :) = &
              & obs_spinz(i_wlk)%obs_latt(:, :, :, :)/dble(obs_spinz(i_wlk)%N)/obs_spinz(i_wlk)%ave_sign

          obs_tau(2)%obs_latt(:, :, :, :) = obs_tau(2)%obs_latt(:, :, :, :) + &
              & obs_spinz(i_wlk)%obs_latt(:, :, :, :)*z_fac

          obs_spinxy(i_wlk)%ave_sign = obs_spinxy(i_wlk)%ave_sign/dble(obs_spinxy(i_wlk)%N)
          obs_spinxy(i_wlk)%obs_latt(:, :, :, :) = &
              & obs_spinxy(i_wlk)%obs_latt(:, :, :, :)/dble(obs_spinxy(i_wlk)%N)/obs_spinxy(i_wlk)%ave_sign

          obs_tau(3)%obs_latt(:, :, :, :) = obs_tau(3)%obs_latt(:, :, :, :) + &
              & obs_spinxy(i_wlk)%obs_latt(:, :, :, :)*z_fac

          obs_spint(i_wlk)%ave_sign = obs_spint(i_wlk)%ave_sign/dble(obs_spint(i_wlk)%N)
          obs_spint(i_wlk)%obs_latt(:, :, :, :) = &
              & obs_spint(i_wlk)%obs_latt(:, :, :, :)/dble(obs_spint(i_wlk)%N)/obs_spint(i_wlk)%ave_sign

          obs_tau(4)%obs_latt(:, :, :, :) = obs_tau(4)%obs_latt(:, :, :, :) + &
              & obs_spint(i_wlk)%obs_latt(:, :, :, :)*z_fac

       end subroutine bp_obsert

       subroutine obsert_mc(nt, gt0, g0t, g00, gtt, overlap_mc)

          use Predefined_Obs

          implicit none

          integer, intent(IN) :: nt
          complex(Kind=kind(0.d0)), intent(IN) :: gt0(ndim, ndim, n_fl, n_grc), g0t(ndim, ndim, n_fl, n_grc)
          complex(Kind=kind(0.d0)), intent(IN) :: g00(ndim, ndim, n_fl, n_grc), gtt(ndim, ndim, n_fl, n_grc)
          complex(Kind=kind(0.d0)), intent(IN) :: overlap_mc(n_grc)

          !Locals
          complex(Kind=kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, ZW, Re_ZW, Z_fac, Z_ol, phase, exp_overlap(N_slat)
          complex(Kind=kind(0.d0)) :: sum_o
          real(Kind=kind(0.d0)) :: X
          integer :: imj, i, j, i1, j1, no_i, no_j, i_grc, i_wlk, ns, i_st, i_ed

          do i_wlk = 1, n_wlk

             i_st = 1 + (i_wlk - 1)*N_slat; i_ed = i_wlk*N_slat
             exp_overlap(:) = exp(overlap_mc(i_st:i_ed))
             sum_o = sum(exp_overlap(:))

             phase = sum_o/abs(sum_o)
             zp = phase/dble(phase)
             zs = dble(phase)/abs(dble(Phase))

             if (nt .eq. 0) then
                obs_grc(i_wlk)%n = obs_grc(i_wlk)%n + 1
                obs_spinz(i_wlk)%n = obs_spinz(i_wlk)%n + 1
                obs_spinxy(i_wlk)%n = obs_spinxy(i_wlk)%n + 1
                obs_spint(i_wlk)%n = obs_spint(i_wlk)%n + 1

                obs_grc(i_wlk)%ave_sign = obs_grc(i_wlk)%ave_sign + dble(zs)
                obs_spinz(i_wlk)%ave_sign = obs_spinz(i_wlk)%ave_sign + dble(zs)
                obs_spinxy(i_wlk)%ave_sign = obs_spinxy(i_wlk)%ave_sign + dble(zs)
                obs_spint(i_wlk)%ave_sign = obs_spint(i_wlk)%ave_sign + dble(zs)
             end if

             do ns = 1, n_slat
                i_grc = ns + (i_wlk - 1)*N_slat
                z_fac = zp*zs*exp(overlap_mc(i_grc))/sum_o

                call predefined_obs_tau_green_measure(Latt, Latt_unit, list, nt, gt0(:, :, :, i_grc),  &
                    &  g0t(:, :, :, i_grc), g00(:, :, :, i_grc), gtt(:, :, :, i_grc), n_sun, z_fac, obs_grc(i_wlk))
                call predefined_obs_tau_spinmz_measure(Latt, Latt_unit, list, nt, gt0(:, :, :, i_grc),  &
                    &  g0t(:, :, :, i_grc), g00(:, :, :, i_grc), gtt(:, :, :, i_grc), n_sun, z_fac, obs_spinz(i_wlk), &
                    &  obs_spinxy(i_wlk), obs_spint(i_wlk))
             end do
          end do

       end subroutine obsert_mc

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Single spin flip S0 ratio
!> @details
!> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
!> a spin flip of Operator n on time slice nt
!> @details
!--------------------------------------------------------------------
       real(Kind=kind(0.d0)) function S0(n, nt, Hs_new)
          implicit none
          !> Operator index
          integer, intent(IN) :: n
          !> Time slice
          integer, intent(IN) :: nt
          !> New local field on time slice nt and operator index n
          real(Kind=kind(0.d0)), intent(In) :: Hs_new

          integer :: nt1, I

          S0 = 1.d0

       end function S0

       complex(Kind=kind(0.d0)) function E0_local(GR)
          implicit none

          complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)

          !Local
          complex(Kind=kind(0.d0)) :: GRC(Ndim, Ndim, N_FL), ZK
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY
          integer :: I, J, imj, nf, dec, I1, J1, no_I, no_J, n
          real(Kind=kind(0.d0)) :: X

          do nf = 1, N_FL
             do I = 1, Ndim
                do J = 1, Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                end do
                GRC(I, I, nf) = 1.d0 + GRC(I, I, nf)
             end do
          end do

          Zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call Predefined_Hoppings_Compute_Kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin*dble(N_SUN)

          ZPot = cmplx(0.d0, 0.d0, kind(0.d0))
          do I = 1, Latt%N
             do no_I = 1, Latt_unit%Norb
                I1 = Invlist(I, no_I)
                ZPot = ZPot + Grc(i1, i1, 1)*Grc(i1, i1, 2)*ham_U
             end do
          end do

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

       subroutine init_obs_mc

          implicit none

         !!local
          integer :: i_wlk

         !! init observables
          do i_wlk = 1, n_wlk
             call obs_grc(i_wlk)%init
             call obs_spinz(i_wlk)%init
             call obs_spinxy(i_wlk)%init
             call obs_spint(i_wlk)%init
          end do

       end subroutine init_obs_mc

    end submodule ham_Hubbard_smod
