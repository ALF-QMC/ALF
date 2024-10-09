!  Copyright (C) 2016 - 2020 The ALF project
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

    submodule(Hamiltonian_main) ham_tV_smod

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

       type, extends(ham_base) :: ham_tV
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: Alloc_obs
          procedure, nopass :: Obser
          procedure, nopass :: ObserT
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_tV

       !#PARAMETERS START# VAR_lattice
       character(len=64) :: Model = ''  ! Value not relevant
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

       !#PARAMETERS START# VAR_tV
       real(Kind=kind(0.d0)) :: ham_T = 1.d0  ! Hopping parameter
       real(Kind=kind(0.d0)) :: Ham_chem = 0.d0  ! Chemical potential
       real(Kind=kind(0.d0)) :: ham_V = 4.d0  ! Hubbard interaction
       real(Kind=kind(0.d0)) :: ham_T2 = 1.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_V2 = 4.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_Tperp = 1.d0  ! For bilayer systems
       real(Kind=kind(0.d0)) :: ham_Vperp = 1.d0  ! For bilayer systems
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable   :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell

    contains

       module subroutine Ham_Alloc_tV
          allocate (ham_tV::ham)
       end subroutine Ham_Alloc_tV

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_tV_read_write_parameters.F90"

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
             write (unit_info, *) 'Ham_V         : ', Ham_V
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
       subroutine Ham_Vint

          use Predefined_Int
          implicit none

          real(Kind=kind(0.d0)), allocatable :: Ham_V_vec(:), Ham_Vperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                  Ham_V2_vec(:), Ham_Lambda_vec(:)
          integer, allocatable ::   N_Phi_vec(:)
          type(Hopping_Matrix_type), allocatable :: Bond_Matrix(:)

          integer                           :: I, J, I1, J1, no_I, no_J, nf
          integer                           :: n_1, n_2, Nb, n_f, l_f, n_l, N, nc
          complex(Kind=kind(0.d0))          :: Z
          real(Kind=kind(0.d0))             :: Zero = 1.0e-6

          allocate (Ham_V_vec(N_FL), Ham_V2_vec(N_FL), Ham_Vperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &                                   N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL))

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_V_vec = Ham_V/dble(N_SUN)
          Ham_V2_vec = Ham_V2/dble(N_SUN)
          Ham_Vperp_vec = Ham_Vperp/dble(N_SUN)
          Ham_Chem_vec = 0.0d0
          Phi_X_vec = Phi_X
          Phi_Y_vec = Phi_Y
          Ham_Lambda_vec = 0.0d0
          N_Phi_vec = N_Phi

          !Use predefined hoppings to manage the bonds since the interaction of the tV model is exactly on the hopping bonds
          select case (str_to_upper(Lattice_type))
          case ("SQUARE")
             call Set_Default_hopping_parameters_square(Bond_Matrix, Ham_V_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          case ("N_LEG_LADDER")
             call Set_Default_hopping_parameters_n_leg_ladder(Bond_Matrix, Ham_V_vec, Ham_Vperp_vec, Ham_Chem_vec, Phi_X_vec, &
                  &                                            Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          case ("HONEYCOMB")
         call Set_Default_hopping_parameters_honeycomb(Bond_Matrix, Ham_V_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                      &                                         Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          case ("BILAYER_SQUARE")
             call Set_Default_hopping_parameters_Bilayer_square(Bond_Matrix, Ham_V_vec, Ham_V2_vec, Ham_Vperp_vec, Ham_Chem_vec, &
                  &                                              Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
                  &                                              List, Invlist, Latt, Latt_unit)

          case ("BILAYER_HONEYCOMB")
            call Set_Default_hopping_parameters_Bilayer_honeycomb(Bond_Matrix, Ham_V_vec, Ham_V2_vec, Ham_Vperp_vec, Ham_Chem_vec, &
                   &                                                 Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
                   &                                                 List, Invlist, Latt, Latt_unit)

          end select

          N = 0
          do n_f = 1, Bond_Matrix(1)%N_FAM
             N = N + Bond_Matrix(1)%L_Fam(n_f)
          end do
          if (N == 0) then
             ! It is a noninteracting model and then there is no need to setup the interaction apart from one vertex per flavor
             ! for internal memory consistency (the code will always access the first vertex OP_V(1,:) hence it has to be allocated)
             ! any vertex with zero coupling strength will do, ie, the hubbard interaction with U=0
             allocate (Op_V(1, N_FL))
             do nf = 1, N_FL
                ! Fake hubbard interaction of weight 0.0 (last argument in the following call)
                call Predefined_Int_U_SUN(OP_V(1, nf), 1, N_SUN, DTAU, 0.0d0)
             end do
             write (*, *) "No interaction present"
             return
          end if

          if (Symm) call Symmetrize_families(Bond_Matrix)
          N = 0
          do n_f = 1, Bond_Matrix(1)%N_FAM
             N = N + Bond_Matrix(1)%L_Fam(n_f)
          end do

          allocate (Op_V(N, N_FL))
          do nf = 1, N_FL
             !               (not needed since we can directly access the Hamiltonian member,
             !               otherwise we risk overwriting stuff)
             !               N_Phi     = Bond_Matrix(nf)%N_Phi
             !               Phi_X     = Bond_Matrix(nf)%Phi_X
             !               Phi_Y     = Bond_Matrix(nf)%Phi_Y
             !               Bulk      = Bond_Matrix(nf)%Bulk
             do nc = 1, size(Op_V, 1)
                call Op_make(Op_V(nc, nf), 2)
             end do
             nc = 0
             do n_f = 1, Bond_Matrix(1)%N_FAM
                do l_f = 1, Bond_Matrix(1)%L_Fam(n_f)
                   I = Bond_Matrix(1)%List_Fam(n_f, l_f, 1)
                   nb = Bond_Matrix(1)%List_Fam(n_f, l_f, 2)
                   no_I = Bond_Matrix(nf)%list(Nb, 1)
                   no_J = Bond_Matrix(nf)%list(Nb, 2)
                   n_1 = Bond_Matrix(nf)%list(Nb, 3)
                   n_2 = Bond_Matrix(nf)%list(Nb, 4)
                   J = Latt%nnlist(I, n_1, n_2)
                   Z = Generic_hopping(I, no_I, n_1, n_2, no_J, N_Phi, Phi_x, Phi_y, Bulk, Latt, Latt_Unit)
                   I1 = Invlist(I, no_I)
                   J1 = Invlist(J, no_J)
                   nc = nc + 1
                   Op_V(nc, nf)%P(1) = I1
                   Op_V(nc, nf)%P(2) = J1
                   Op_V(nc, nf)%O(1, 2) = Z
                   Op_V(nc, nf)%O(2, 1) = conjg(Z)
                   Op_V(nc, nf)%g = sqrt(-Dtau*Bond_Matrix(nf)%T(Nb)*Bond_Matrix(1)%Prop_Fam(n_f))
                   Op_V(nc, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
                   Op_V(nc, nf)%type = 2
                   call Op_set(Op_V(nc, nf))
                end do
             end do
          end do

          deallocate (Ham_V_vec, Ham_V2_vec, Ham_Vperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                   N_Phi_vec, Ham_Lambda_vec)

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
          complex(Kind=kind(0.d0)) :: GRC(Ndim, Ndim, N_FL), ZK, Zn, weight, delta
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY, tmp
          integer :: I, J, imj, nf, dec, I1, J1, no_I, no_J, n, nf2, k, k1, l, l1
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

          Zn = cmplx(dble(N_sun), 0.d0, kind(0.d0))
          ZPot = cmplx(0.d0, 0.d0, kind(0.d0))
          do nf = 1, N_FL
             do n = 1, size(OP_V, 1)
                weight = -Op_V(n, nf)%g**2/dtau
                do J = 1, Op_V(n, nf)%N
                   J1 = Op_V(n, nf)%P(J)
                   do I = 1, Op_V(n, nf)%N
                      if (abs(Op_V(n, nf)%O(i, j)) >= 0.00001) then
                         I1 = Op_V(n, nf)%P(I)
                         ZPot = ZPot + N_FL*2.d0*Zn*Op_V(n, nf)%alpha*weight*Op_V(n, nf)%O(i, j)*GRC(I1, J1, nf)

                         do nf2 = 1, N_FL
                            Delta = 0.d0
                            if (nf == nf2) Delta = 1.d0
                            do K = 1, Op_V(n, nf)%N
                               K1 = Op_V(n, nf)%P(K)
                               do L = 1, Op_V(n, nf)%N
                                  if (abs(Op_V(n, nf)%O(k, l)) >= 0.00001) then
                                     L1 = Op_V(n, nf)%P(L)
                                     tmp = (delta*GRC(I1, L1, nf)*GR(J1, K1, nf) +  &
                                           &    Zn*GRC(I1, J1, nf)*GRC(K1, L1, nf2))
                                     ZPot = ZPot + weight*Op_V(n, nf)%O(i, j)*Op_V(n, nf)%O(k, l)*tmp
                                  end if
                               end do
                            end do
                         end do
                      end if
                   end do
                end do
!                 if ( .not. Projector)
                ZPot = ZPot + N_FL*weight*(Op_V(n, nf)%alpha**2)*Zn
             end do
          end do
          Zpot = Zn*Zpot
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
          call Predefined_Obs_eq_Green_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(1))
          call Predefined_Obs_eq_SpinSUN_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(2))
          call Predefined_Obs_eq_Den_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(3))

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

          call Predefined_Obs_tau_Green_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(1))
          call Predefined_Obs_tau_SpinSUN_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(2))
          call Predefined_Obs_tau_Den_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs_tau(3))

       end subroutine OBSERT

    end submodule ham_tV_smod
