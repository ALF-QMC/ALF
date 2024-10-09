!  Copyright (C) 2016 - 2023 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module defines the  Hamiltonian and observables.  Here, we have included a
!> set of predefined Hamiltonians. They include the Hubbard and SU(N) tV models
!> on honeycomb, pi-flux and square lattices.

!> @details
!> The public variables of this module are the following
!>
!>
!> @param [public] OP_V
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> List of operators of type=1,2 and 3 describing the sequence of interactions on a time slice.
!> The first index runs over this sequence. The second corresponds to the flavor index.  \endverbatim
!>
!> @param [public] OP_T
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> Sequence of  operators  accounting for the  hopping on a  time slice. This can include  various
!> checkerboard decompositions. The first index runs over this sequence. The second corresponds to
!> the flavor index. \endverbatim
!> *  The progagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n}  \f$.  That is
!> first the hopping and then the potential energy.
!>
!>@param [public] WF_L
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma
!> \verbatim Type(Fields)
!> Contains all auxiliary fields in the variable f(:,:). The first index runs through the operator
!> sequence. The second through the time slices.   \endverbatim
!
!> @param [public]  Ndim
!> \verbatim Integer
!> Total number of orbitals. e.g. # unit cells * # orbitals per unit cell.  \endverbatim
!
!> @param [public]  N_FL
!> \verbatim Integer
!> # of flavors.  Propagation is block diagonal in flavors.  \endverbatim
!
!> @param [public]  N_SUN
!> \verbatim Integer
!> # of colors.  Propagation is color independent.  \endverbatim
!>
!> @param [public] Ltrot
!> \verbatim Integer
!> Available measurment interval in units of Delta Tau. \endverbatim
!>
!> @param [public] Thtrot
!>  \verbatim Integer
!> Effective projection parameter in units of Delta Tau.  (Only relevant if projective option is turned on) \endverbatim
!>
!> @param [public] Projector
!> \verbatim Logical
!> Flag for projector. If true then the total number of time slices will correspond to Ltrot + 2*Thtrot \endverbatim
!>
!> @param [public] Group_Comm
!> \verbatim Integer
!> Defines MPI communicator  \endverbatim
!
!> @param [public] Symm
!> \verbatim Logical  \endverbatim
!> If set to true then the green functions will be symmetrized
!> before being  sent to the Obser, ObserT subroutines.
!> In particular, the transformation,  \f$ \tilde{G} =  e^{-\Delta \tau T /2 } G e^{\Delta \tau T /2 } \f$
!> will be carried out  and \f$ \tilde{G} \f$  will be sent to the Obser and ObserT subroutines.  Note that
!> if you want to use this  feature, then you have to be sure the hopping and interaction terms are decomposed
!> symmetrically. If Symm is true, the propagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=N_T}^{1}e^{T_n/2} \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n/2}  \f$
!>
!>
!> You still have to add some docu for the other private variables in this module.
!>
!--------------------------------------------------------------------

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
       use LRC_Mod

       implicit none

       type, extends(ham_base) :: ham_Hubbard
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: Alloc_obs
          procedure, nopass :: Obser
          procedure, nopass :: ObserT
          procedure, nopass :: S0
          procedure, nopass :: Ham_Langevin_HMC_S0
          procedure, nopass :: Delta_S0_global
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
       !Integer              :: N_SUN        = 2        ! Number of colors
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

       !#PARAMETERS START# VAR_Hubbard
       real(Kind=kind(0.d0)) :: ham_T = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_chem = 0.d0     ! Chemical potential
       real(Kind=kind(0.d0)) :: ham_U = 4.d0     ! Hubbard interaction
       real(Kind=kind(0.d0)) :: ham_T2 = 1.d0     ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_U2 = 4.d0     ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_Tperp = 1.d0     ! For bilayer systems
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell

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

          ! From dynamically generated file "Hamiltonian_Hubbard_read_write_parameters.F90"
          call read_parameters()

          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot + 2*Thtrot
          if (Mz) then
             N_FL = 2
             if (mod(N_SUN, 2) .ne. 0) then
                write (error_unit, *) 'Ham_Set: N_SUN has to be even if Mz = True'
                call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
             end if
             N_SUN = N_SUN/2
          else
             N_FL = 1
          end if

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
             if (Mz) then
                write (unit_info, *) 'HS  couples to z-component of spin'
             else
                write (unit_info, *) 'HS  couples to density'
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
             if (Mz) then
                write (unit_info, *) 'N_SUN         : ', 2*N_SUN
             else
                write (unit_info, *) 'N_SUN         : ', N_SUN
             end if
             write (unit_info, *) 'N_FL          : ', N_FL
             write (unit_info, *) 't             : ', Ham_T
             write (unit_info, *) 'Ham_U         : ', Ham_U
             write (unit_info, *) 't2            : ', Ham_T2
             write (unit_info, *) 'Ham_U2        : ', Ham_U2
             write (unit_info, *) 'Ham_tperp     : ', Ham_tperp
             write (unit_info, *) 'Ham_chem      : ', Ham_chem
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

          select case (str_to_upper(Lattice_type))
          case ("SQUARE")
             call Set_Default_hopping_parameters_square(Hopping_Matrix, Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          case ("N_LEG_LADDER")
             call Set_Default_hopping_parameters_n_leg_ladder(Hopping_Matrix, Ham_T_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
                  &                                            Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          case ("HONEYCOMB")
             Ham_Lambda = 0.d0
      call Set_Default_hopping_parameters_honeycomb(Hopping_Matrix, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                         &                                         Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          case ("BILAYER_SQUARE")
            call Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix, Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, &
                   &                                              Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
                   &                                              List, Invlist, Latt, Latt_unit)

          case ("BILAYER_HONEYCOMB")
         call Set_Default_hopping_parameters_Bilayer_honeycomb(Hopping_Matrix, Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, &
                      &                                                 Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
                      &                                                 List, Invlist, Latt, Latt_unit)

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
          N_part = Ndim/2
          call Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
               &                            N_part, N_FL, WF_L, WF_R)

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
          if (str_to_upper(Lattice_type) == "BILAYER_SQUARE" .or. str_to_upper(Lattice_type) == "BILAYER_HONEYCOMB") then
             do no = 1, Latt_unit%Norb/2
                Ham_U_vec(no) = Ham_U
                Ham_U_vec(no + Latt_unit%Norb/2) = Ham_U2
             end do
             if (abs(Ham_U) > Zero) N_ops = N_ops + Latt%N*Latt_unit%Norb/2
             if (abs(Ham_U2) > Zero) N_ops = N_ops + Latt%N*Latt_unit%Norb/2
          else
             Ham_U_vec(:) = Ham_U
             if (abs(Ham_U) > Zero) N_ops = N_ops + Latt%N*Latt_unit%Norb
          end if
          if (Mz) then
             allocate (Op_V(N_ops, N_FL))
             Ham_U_vec = Ham_U_vec/real(N_SUN, kind(0.d0))
             nc = 0
             do I1 = 1, Latt%N
                do no = 1, Latt_unit%Norb
                   I = invlist(I1, no)
                   if (abs(Ham_U_vec(no)) > Zero) then
                      nc = nc + 1
                      if (Continuous) then
                         call Predefined_Int_U_MZ_continuous_HS(OP_V(nc, 1), OP_V(nc, 2), I, DTAU, Ham_U_vec(no))
                      else
                         call Predefined_Int_U_MZ(OP_V(nc, 1), OP_V(nc, 2), I, DTAU, Ham_U_vec(no))
                      end if
                   end if
                end do
             end do
          else
             allocate (Op_V(N_ops, N_FL))
             nc = 0
             do I1 = 1, Latt%N
                do no = 1, Latt_unit%Norb
                   I = invlist(I1, no)
                   if (abs(Ham_U_vec(no)) > Zero) then
                      nc = nc + 1
                      if (Continuous) then
                         call Predefined_Int_U_SUN_continuous_HS(OP_V(nc, 1), I, N_SUN, DTAU, Ham_U_vec(no))
                      else
                         call Predefined_Int_U_SUN(OP_V(nc, 1), I, N_SUN, DTAU, Ham_U_vec(no))
                      end if
                   end if
                end do
             end do
          end if

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
       subroutine Alloc_obs(Ltau)

          implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          integer, intent(In) :: Ltau
          integer    ::  i, N, Nt
          character(len=64) ::  Filename
          character(len=:), allocatable ::  Channel

          ! Scalar observables
          allocate (Obs_scal(4))
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
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             call Obser_Vec_make(Obs_scal(I), N, Filename)
          end do

          ! Equal time correlators
          if (Mz) then
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
                   Nt = Ltrot + 1 - 2*Thtrot
                   if (Projector) Channel = 'T0'
                   call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
                end do
             end if
          else
             ! Equal time correlators
             allocate (Obs_eq(3))
             do I = 1, size(Obs_eq, 1)
                select case (I)
                case (1)
                   Filename = "Green"
                case (2)
                   Filename = "SpinZ"
                case (3)
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
          complex(Kind=kind(0.d0)) :: GRC(Ndim, Ndim, N_FL), ZK
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY
          integer :: I, J, imj, nf, dec, I1, J1, no_I, no_J, n
          real(Kind=kind(0.d0)) :: X

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

          ZPot = cmplx(0.d0, 0.d0, kind(0.d0))
          dec = 1
          if (Mz) dec = 2
          if (str_to_upper(Lattice_type) == "BILAYER_SQUARE" .or. str_to_upper(Lattice_type) == "BILAYER_HONEYCOMB") then
             do I = 1, Latt%N
                do no_I = 1, Latt_unit%Norb
                   I1 = Invlist(I, no_I)
                   if (no_I == 1) ZPot = ZPot + Grc(i1, i1, 1)*Grc(i1, i1, dec)*ham_U
                   if (no_I == 2) ZPot = ZPot + Grc(i1, i1, 1)*Grc(i1, i1, dec)*ham_U2
                end do
             end do
          else
             do I = 1, Latt%N
                do no_I = 1, Latt_unit%Norb
                   I1 = Invlist(I, no_I)
                   ZPot = ZPot + Grc(i1, i1, 1)*Grc(i1, i1, dec)*ham_U
                end do
             end do
          end if
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

          ! Standard two-point correlations
          if (Mz) then
             call Predefined_Obs_eq_Green_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(1))
             call Predefined_Obs_eq_SpinMz_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(2), Obs_eq(3), Obs_eq(4))
             call Predefined_Obs_eq_Den_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(5))
          else
             call Predefined_Obs_eq_Green_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(1))
             call Predefined_Obs_eq_SpinSUN_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(2))
             call Predefined_Obs_eq_Den_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(3))
          end if

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
          integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/real(Phase, kind(0.d0))
          ZS = real(Phase, kind(0.d0))/abs(real(Phase, kind(0.d0)))
          ZS = ZS*Mc_step_weight

          ! Standard two-point correlations

          if (Mz) then
             call Predefined_Obs_tau_Green_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(1))
             call Predefined_Obs_tau_SpinMz_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(2),&
                  &                                   Obs_tau(3), Obs_tau(4))
             call Predefined_Obs_tau_Den_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(5))
          else
             call Predefined_Obs_tau_Green_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(1))
             call Predefined_Obs_tau_SpinSUN_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(2))
             call Predefined_Obs_tau_Den_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(3))
          end if

       end subroutine OBSERT

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
          complex(Kind=kind(0.d0)), intent(In) :: Hs_new

          integer :: nt1, I

          if (Continuous) then
             S0 = exp((-real(Hs_new, kind(0.d0))**2 + real(nsigma%f(n, nt), kind(0.d0))**2)/2.d0)
          else
             S0 = 1.d0
          end if
       end function S0

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>   Forces_0  = \partial S_0 / \partial s  are calculated and returned to  main program.
!>
!-------------------------------------------------------------------
       subroutine Ham_Langevin_HMC_S0(Forces_0)

          implicit none

          real(Kind=kind(0.d0)), intent(inout), allocatable :: Forces_0(:, :)
          !Local
          integer :: N, N_op, nt

          ! Compute \partial S_0 / \partial s
          N_op = size(nsigma%f, 1)
          Forces_0 = 0.d0
          do n = 1, N_op
             if (OP_V(n, 1)%type == 3) then
                do nt = 1, Ltrot
                   Forces_0(n, nt) = real(nsigma%f(n, nt))
                end do
             end if
          end do

       end subroutine Ham_Langevin_HMC_S0

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes the ratio exp(S0(new))/exp(S0(old))
!>
!> @details
!> This function computes the ratio \verbatim  e^{-S0(nsigma)}/e^{-S0(nsigma_old)} \endverbatim
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!-------------------------------------------------------------------
       real(Kind=kind(0.d0)) function Delta_S0_global(Nsigma_old)

          !  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          implicit none

          ! Arguments
          type(Fields), intent(IN) :: nsigma_old
          real(kind=kind(0.0d0))     :: S0_old, S0_new
          integer                    :: f, t, nfield, ntau

          Delta_S0_global = 1.d0
          nfield = size(nsigma%f, 1)
          ntau = size(nsigma%f, 2)
          S0_old = 0.0d0
          S0_new = 0.0d0
          do t = 1, ntau
             do f = 1, nfield
                S0_old = S0_old + real(nsigma_old%f(f, t), kind(0.d0))**2
                S0_new = S0_new + real(nsigma%f(f, t), kind(0.d0))**2
             end do
          end do
          S0_old = 0.5d0*S0_old
          S0_new = 0.5d0*S0_new
          Delta_S0_global = exp(-S0_new + S0_old)
          !   write(*,*) "S0 old:", S0_old, "S0 new:", S0_new
          ! S0 = exp( (-Hs_new**2  + nsigma%f(n,nt)**2 ) /2.d0 )

       end function Delta_S0_global

    end submodule ham_Hubbard_smod
