!  Copyright (C) 2022 The ALF project
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
!> This File is a template for defining new models.
!> One can define a new model class by copying this file, replacing alle occurences
!> of ##NAME## by the Hamiltonian name, populating the subroutines below as needed
!> adding the Hamiltonian name to the file Prog/Hamiltonians.list.

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

    submodule(Hamiltonian_main) ham_##NAME##_smod

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

       type, extends(ham_base) :: ham_##NAME##
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: Alloc_obs
          procedure, nopass :: Obser
          procedure, nopass :: ObserT
          ! procedure, nopass :: ##PROCEDURE_NAME##  ! Some other procedure defined in ham_base
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_##NAME##

       !#PARAMETERS START# VAR_##XYZ##
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell

    contains

       module subroutine Ham_Alloc_##NAME##
          allocate (ham_##NAME##::ham)
       end subroutine Ham_Alloc_##NAME##

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_##NAME##_read_write_parameters.F90"

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

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif

!           ! From dynamically generated file "Hamiltonian_##NAME##_read_write_parameters.F90"
!           call read_parameters()
!
!           Ltrot = nint(beta/dtau)
!           if (Projector) Thtrot = nint(theta/dtau)
!           Ltrot = Ltrot+2*Thtrot
!
!           ! Setup the Bravais lattice
!           call Ham_Latt
!
!           ! Setup the hopping / single-particle part
!           call Ham_Hop
!
!           ! Setup the interaction.
!           call Ham_V
!
!           ! Setup the trival wave function, in case of a projector approach
!           if (Projector) Call Ham_Trial()

#ifdef MPI
          if (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write (File_info, '(A,I0,A)') "Temp_", igroup, "/info"
#endif
             open (newunit=unit_info, file=file_info, status="unknown", position="append")
!              Write(unit_info,*) '====================================='
!              Write(unit_info,*) 'Model is      :  ##NAME##'
!              Write(unit_info,*) 'Lattice is    : ', Lattice_type
!              Write(unit_info,*) '# unit cells  : ', Latt%N
!              Write(unit_info,*) '# of orbitals : ', Latt_unit%Norb
!              Write(unit_info,*) 'Checkerboard  : ', Checkerboard
!              Write(unit_info,*) 'Symm. decomp  : ', Symm
!              if (Projector) then
!                 Write(unit_info,*) 'Projective version'
!                 Write(unit_info,*) 'Theta         : ', Theta
!                 Write(unit_info,*) 'Tau_max       : ', beta
!              else
!                 Write(unit_info,*) 'Finite temperture version'
!                 Write(unit_info,*) 'Beta          : ', Beta
!              endif
!              Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
!              Write(unit_info,*) 'N_SUN         : ',   N_SUN
!              Write(unit_info,*) 'N_FL          : ', N_FL
!              if (Projector) then
!                 Do nf = 1,N_FL
!                    Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
!                    Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
!                 enddo
!              endif
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

!           ! Scalar observables
!           Allocate ( Obs_scal(4) )
!           Do I = 1,Size(Obs_scal,1)
!             select case (I)
!             case (1)
!               N = 1;   Filename = "Kin"
!             case (2)
!               N = 1;   Filename = "Pot"
!             case (3)
!               N = 1;   Filename = "Part"
!             case (4)
!               N = 1;   Filename = "Ener"
!             case default
!               Write(6,*) ' Error in Alloc_obs '
!             end select
!             Call Obser_Vec_make(Obs_scal(I),N,Filename)
!           enddo
!
!           ! Equal time correlators
!           Allocate ( Obs_eq(3) )
!           Do I = 1,Size(Obs_eq,1)
!             select case (I)
!             case (1)
!               Filename = "Green"
!             case (2)
!               Filename = "SpinZ"
!             case (3)
!               Filename = "Den"
!             case default
!               Write(6,*) ' Error in Alloc_obs '
!             end select
!             Nt = 1
!             Channel = '--'
!             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
!           enddo
!
!           If (Ltau == 1) then
!             ! Time-displaced correlators
!             Allocate ( Obs_tau(3) )
!             Do I = 1,Size(Obs_tau,1)
!               select case (I)
!               case (1)
!                 Channel = 'P' ; Filename = "Green"
!               case (2)
!                 Channel = 'PH'; Filename = "SpinZ"
!               case (3)
!                 Channel = 'PH'; Filename = "Den"
!               case default
!                 Write(6,*) ' Error in Alloc_obs '
!               end select
!               Nt = Ltrot+1-2*Thtrot
!               If(Projector) Channel = 'T0'
!               Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
!             enddo
!           endif

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
          integer, intent(IN) :: Ntau
          real(Kind=kind(0.d0)), intent(IN) :: Mc_step_weight

          !Local
          complex(Kind=kind(0.d0)) :: GRC(Ndim, Ndim, N_FL)
          complex(Kind=kind(0.d0)) :: ZP, ZS,
          integer :: I, J, nf
          ! Add local variables as needed

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

          ! Compute equal-time correlations

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
          complex(Kind=kind(0.d0)) :: ZP, ZS
          ! Add local variables as needed

          ZP = PHASE/real(Phase, kind(0.d0))
          ZS = real(Phase, kind(0.d0))/abs(real(Phase, kind(0.d0)))
          ZS = ZS*Mc_step_weight

          ! Compute observables

       end subroutine OBSERT

    end submodule ham_##NAME##_smod
