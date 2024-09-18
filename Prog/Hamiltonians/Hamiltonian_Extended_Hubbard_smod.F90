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

    submodule (Hamiltonian_main) ham_extended_Hubbard_smod

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3
      Use MyMats
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Fields_mod
      Use Predefined_Hoppings
      Use LRC_Mod

      Implicit none
      
      type, extends(ham_base) :: ham_extended_Hubbard
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: S0
        procedure, nopass :: Ham_Langevin_HMC_S0
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_extended_Hubbard

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = 'Extended_Hubbard'  ! Value not relevant
      Character (len=64) :: Lattice_type = 'Square'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Model_Generic
      !Integer              :: N_SUN        = 2        ! Number of colors
      !Integer              :: N_FL         = 1        ! Number of flavors
      real(Kind=Kind(0.d0)) :: Phi_X        = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
      real(Kind=Kind(0.d0)) :: Phi_Y        = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
      logical               :: Bulk         = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
      Integer               :: N_Phi        = 0        ! Total number of flux quanta traversing the lattice
      real(Kind=Kind(0.d0)) :: Dtau         = 0.1d0    ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta         = 5.d0     ! Inverse temperature
      logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
      !logical              :: Symm         = .true.   ! Whether symmetrization takes place
      !logical              :: Projector    = .false.  ! Whether the projective algorithm is used
      real(Kind=Kind(0.d0)) :: Theta        = 10.d0    ! Projection parameter
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_extended_Hubbard
      real(Kind=Kind(0.d0)) :: ham_T      = 1.d0     ! Hopping parameter
      real(Kind=Kind(0.d0)) :: Ham_chem   = 0.d0     ! Chemical potential
      real(Kind=Kind(0.d0)) :: Ham_U      = 4.d0     ! Hubbard interaction
      real(Kind=Kind(0.d0)) :: Ham_V1     = 1.d0     ! nearest neighbor interaction
      real(Kind=Kind(0.d0)) :: ham_T2     = 1.d0     ! For bilayer systems
      real(Kind=Kind(0.d0)) :: Ham_U2     = 4.d0     ! For bilayer systems
      real(Kind=Kind(0.d0)) :: Ham_V2     = 1.d0     ! For bilayer systems
      real(Kind=Kind(0.d0)) :: ham_Tperp  = 1.d0     ! For bilayer systems
      real(Kind=Kind(0.d0)) :: Ham_Vperp  = 1.d0     ! For bilayer systems
      !! ATTENTION !!
      !! Check if MZ can work!!
      logical               :: Mz         = .false.   ! When true, sets the M_z-Hubbard model: Nf=2, demands that N_sun is even, HS field couples to the z-component of magnetization; otherwise, HS field couples to the density
      logical               :: Continuous = .false.   ! Uses (T: continuous; F: discrete) HS transformation
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:) ! used for hopping terms
      Type (Hopping_Matrix_type), Allocatable :: Bond_Matrix(:)    ! used for bond interactions
      Integer, allocatable :: List(:,:), Invlist(:,:)              ! For orbital structure of Unit cell

    contains
      
      module Subroutine Ham_Alloc_extended_Hubbard
        allocate(ham_extended_Hubbard::ham)
      end Subroutine Ham_Alloc_extended_Hubbard

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Extended_Hubbard_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Implicit none

          integer                :: ierr, nf, unit_info
          Character (len=64)     :: file_info


          ! L1, L2, Lattice_type, List(:,:), Invlist(:,:) -->  Lattice information
          ! Ham_T, Chem, Phi_X, XB_B, Checkerboard, Symm   -->  Hopping
          ! Interaction                              -->  Model
          ! Simulation type                          -->  Finite  T or Projection  Symmetrize Trotter.
                  
#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          ! From dynamically generated file "Hamiltonian_extended_Hubbard_read_write_parameters.F90"
          call read_parameters()

          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot
          If ( Mz ) then
            N_SUN = 1
            N_FL  = 2
            !Write(error_unit,*) 'Ham_Set: Mz = True is not  yet  tested. '
            !CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          else
            N_FL  = 1
          endif

          ! Setup the Bravais lattice
          Call  Ham_Latt

          ! Setup the hopping / single-particle part
          Call  Ham_Hop

          ! Setup the interaction.
          call Ham_V

          ! Setup the trival wave function, in case of a projector approach
          if (Projector) Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      :  Extended Hubbard Model' 
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) '# unit cells  : ', Latt%N 
             Write(unit_info,*) '# of orbitals : ', Latt_unit%Norb
             Write(unit_info,*) 'Flux_1        : ', Phi_X
             Write(unit_info,*) 'Flux_2        : ', Phi_Y
             If (Bulk) then
                Write(unit_info,*) 'Twist as phase factor in bulk'
             Else
                Write(unit_info,*) 'Twist as boundary condition'
             endif
             If ( Mz )  then
                Write(unit_info,*) 'HS couples to z-component of spin'
             else
                Write(unit_info,*) 'HS  ccouples to density'
             endif
             Write(unit_info,*) 'Checkerboard  : ', Checkerboard
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
             else
                Write(unit_info,*) 'Finite temperture version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             if ( Mz )  then
                Write(unit_info,*) 'N_SUN         : ', 2*N_SUN
             else
                Write(unit_info,*) 'N_SUN         : ',   N_SUN
             endif
             Write(unit_info,*) 'N_FL          : ', N_FL
             Write(unit_info,*) 't             : ', Ham_T
             Write(unit_info,*) 'Ham_U         : ', Ham_U
             Write(unit_info,*) 'Ham_V1        : ', Ham_V1
             if (Index(str_to_upper(Lattice_type),'BILAYER') > 0 )  then
               Write(unit_info,*) 't2            : ', Ham_T2
               Write(unit_info,*) 'Ham_U2        : ', Ham_U2
               Write(unit_info,*) 'Ham_V2        : ', Ham_V2
               Write(unit_info,*) 'Ham_tperp     : ', Ham_tperp
               Write(unit_info,*) 'Ham_Vperp     : ', Ham_Vperp
             endif   
             Write(unit_info,*) 'Ham_chem      : ', Ham_chem
             if (Projector) then
                Do nf = 1,N_FL
                   Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
             endif
             Close(unit_info)
#ifdef MPI
          Endif
#endif
        end Subroutine Ham_Set

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
        Subroutine Ham_Latt

          Use Predefined_Lattices

          Implicit none
          ! Use predefined stuctures or set your own lattice.
          Call Predefined_Latt(Lattice_type, L1,L2,Ndim, List,Invlist,Latt,Latt_Unit)

        end Subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the Hopping
!--------------------------------------------------------------------
        Subroutine Ham_Hop

          Implicit none

          Real (Kind=Kind(0.d0) ) ::  Ham_Lambda = 0.d0

          Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                  Ham_T2_vec(:),  Ham_Lambda_vec(:)
          Integer, allocatable ::   N_Phi_vec(:)

          ! Use predefined stuctures or set your own hopping
          Integer :: n,nth

          Allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &                                   N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL) )

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_T_vec      = Ham_T
          Ham_Tperp_vec  = Ham_Tperp
          Ham_Chem_vec   = Ham_Chem
          Phi_X_vec      = Phi_X
          Phi_Y_vec      = Phi_Y
          Ham_T2_vec     = Ham_T2
          Ham_Lambda_vec = Ham_Lambda
          N_Phi_vec      = N_Phi

          Select case (str_to_upper(Lattice_type))
          Case ("SQUARE")
             Call  Set_Default_hopping_parameters_square(Hopping_Matrix,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("N_LEG_LADDER")
             Call  Set_Default_hopping_parameters_n_leg_ladder(Hopping_Matrix, Ham_T_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
                  &                                            Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("HONEYCOMB")
             Ham_Lambda = 0.d0
             Call  Set_Default_hopping_parameters_honeycomb(Hopping_Matrix, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                         Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("BILAYER_SQUARE")
             Call  Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
                  &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                  &                                              List, Invlist, Latt, Latt_unit )

          Case ("BILAYER_HONEYCOMB")
             Call  Set_Default_hopping_parameters_Bilayer_honeycomb(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
                  &                                                 Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                  &                                                 List, Invlist, Latt, Latt_unit )

          end Select

          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_T )

          Deallocate (Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                   N_Phi_vec,  Ham_Lambda_vec )

        end Subroutine Ham_Hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the trial wave function
!--------------------------------------------------------------------
        Subroutine Ham_Trial()

          Use Predefined_Trial

          Implicit none

          Integer :: N_part, nf
          ! Use predefined stuctures or set your own Trial  wave function
          N_part = Ndim/2
          Call Predefined_TrialWaveFunction(Lattice_type ,Ndim,  List,Invlist,Latt, Latt_unit, &
               &                            N_part, N_FL,  WF_L, WF_R)

        end Subroutine Ham_Trial

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
        Subroutine Ham_V

          Use Predefined_Int
          Implicit none

          Integer :: nf, I, I1, I2,  nc,  J, no,  N_ops, N_hubbard, N_Phi_loc
          Real (Kind=Kind(0.d0)) :: X,  Zero = 1.D-10,  V_eff
          Real (Kind=Kind(0.d0)), allocatable :: Ham_U_vec(:)

          Real (Kind=Kind(0.d0) ), allocatable :: Ham_V_vec(:), Ham_Vperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                  Ham_V2_vec(:),  Ham_Lambda_vec(:)
          Integer, allocatable ::   N_Phi_vec(:)

          Integer                           :: J1, no_I, no_J
          Integer                           :: n_1, n_2, Nb, n_f,l_f, n_l, N
          Real (Kind=Kind(0.d0))            :: Phi_X_loc, Phi_Y_loc
          Complex(Kind=Kind(0.d0))          :: Z


          Allocate (Ham_U_vec(Latt_unit%Norb))


          Allocate (Ham_V_vec(N_FL), Ham_V2_vec(N_FL), Ham_Vperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &                                   N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL) )

          ! Here we consider no N_FL  dependence of the hopping parameters.

          if (Mz) then
            Ham_V_vec      = -Ham_V1   /dble(2*N_FL)    !  Hopping  is  defined as -t 
            Ham_V2_vec     = -Ham_V2   /dble(2*N_FL)    !  The SU(2) spin symmetry is built into the flavor index
            Ham_Vperp_vec  = -Ham_Vperp/dble(2*N_FL)
          else
            Ham_V_vec      = -Ham_V1   /dble(2*N_SUN)    !  Hopping  is  defined as -t 
            Ham_V2_vec     = -Ham_V2   /dble(2*N_SUN)
            Ham_Vperp_vec  = -Ham_Vperp/dble(2*N_SUN)
          endif
          Ham_Chem_vec   = 0.0d0
          Phi_X_vec      = 0.d0
          Phi_Y_vec      = 0.d0
          Ham_Lambda_vec = 0.0d0
          N_Phi_vec      = 0


          !Use predefined hoppings to manage the bonds since the interaction of the tV model is exactly on the hopping bonds
          Select case (str_to_upper(Lattice_type))
          Case ("SQUARE")
             Call  Set_Default_hopping_parameters_square(Bond_Matrix,Ham_V_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("N_LEG_LADDER")
             Call  Set_Default_hopping_parameters_n_leg_ladder(Bond_Matrix, Ham_V_vec, Ham_Vperp_vec, Ham_Chem_vec, Phi_X_vec, &
                  &                                            Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("HONEYCOMB")
             Call  Set_Default_hopping_parameters_honeycomb(Bond_Matrix, Ham_V_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                         Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("BILAYER_SQUARE")
             Call  Set_Default_hopping_parameters_Bilayer_square(Bond_Matrix,Ham_V_vec,Ham_V2_vec,Ham_Vperp_vec, Ham_Chem_vec, &
                  &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                  &                                              List, Invlist, Latt, Latt_unit )

          Case ("BILAYER_HONEYCOMB")
             Call  Set_Default_hopping_parameters_Bilayer_honeycomb(Bond_Matrix,Ham_V_vec,Ham_V2_vec,Ham_Vperp_vec, Ham_Chem_vec, &
                  &                                                 Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                  &                                                 List, Invlist, Latt, Latt_unit )

          end Select

          N_ops = 0
          if ( str_to_upper(Lattice_type) == "BILAYER_SQUARE" .or. str_to_upper(Lattice_type) =="BILAYER_HONEYCOMB" ) then
             Do no = 1,  Latt_unit%Norb/2
                Ham_U_vec(no                    ) = Ham_U  
                Ham_U_vec(no + Latt_unit%Norb/2 ) = Ham_U2 
             enddo
             If (abs(Ham_U ) > Zero ) N_ops = N_ops + Latt%N*Latt_unit%Norb/2
             If (abs(Ham_U2) > Zero ) N_ops = N_ops + Latt%N*Latt_unit%Norb/2
          else
             Ham_U_vec(:)  = Ham_U
             If (abs(Ham_U ) > Zero ) N_ops = N_ops + Latt%N*Latt_unit%Norb
          endif

          do n_f = 1,Bond_Matrix(1)%N_FAM
            if (Mz) then
              ! We use real field for each kind of interaction: n_{i, sigma} n_{j, sigma'}
              N_ops = N_ops + 4*Bond_Matrix(1)%L_Fam(n_f)
            else
              N_ops = N_ops + Bond_Matrix(1)%L_Fam(n_f)
            endif
          enddo
          Allocate(Op_V(N_ops,N_FL))

          If ( Mz )  Then
             nc = 0
             Do I1 = 1,Latt%N
                do no = 1, Latt_unit%Norb
                   I = invlist(I1,no)
                   if (abs(Ham_U_vec(no)) > Zero ) then
                      nc = nc + 1
                      if (Continuous) then
                         Call Predefined_Int_U_MZ_continuous_HS( OP_V(nc,1), OP_V(nc,2), I,  DTAU, Ham_U_vec(no) )
                      else
                         Call Predefined_Int_U_MZ              ( OP_V(nc,1), OP_V(nc,2), I,  DTAU, Ham_U_vec(no) )
                      endif
                   endif
                enddo
             enddo
          else
             nc = 0
             Do I1 = 1,Latt%N
                do no = 1, Latt_unit%Norb
                   I = invlist(I1,no)
                   if (abs(Ham_U_vec(no)) > Zero ) then
                      nc = nc + 1
                      if (Continuous) then
                         Call Predefined_Int_U_SUN_continuous_HS(  OP_V(nc,1), I, N_SUN, DTAU, Ham_U_vec(no)  )
                      else
                         Call Predefined_Int_U_SUN              (  OP_V(nc,1), I, N_SUN, DTAU, Ham_U_vec(no)  )
                      endif
                   endif
                Enddo
             Enddo
          Endif

          N_hubbard = nc
          ! nf is the flavor index. n_f is the familly index.
          if (Mz) then

            do n_f = 1, Bond_Matrix(1)%N_FAM
               do l_f = 1, Bond_Matrix(1)%L_Fam(n_f)
                  I    = Bond_Matrix(1)%List_Fam(n_f,l_f,1)
                  nb   = Bond_Matrix(1)%List_Fam(n_f,l_f,2)
                  no_I = Bond_Matrix(1)%list(Nb,1)
                  no_J = Bond_Matrix(1)%list(Nb,2)
                  n_1  = Bond_Matrix(1)%list(Nb,3)
                  n_2  = Bond_Matrix(1)%list(Nb,4)
                  J    = Latt%nnlist(I,n_1,n_2)
                  I1   = Invlist(I,no_I)
                  J1   = Invlist(J,no_J)
                  ! Up-Up 
                  nc = nc + 1
                  do nf = 1,N_FL
                     call Op_make(Op_V(nc, nf), 2)
                     Op_V(nc, nf)%P(1) = I1
                     Op_V(nc, nf)%P(2) = J1
                     Op_V(nc, nf)%type = 2
                  enddo
                  Op_V(nc,1)%O(1,1) = cmplx( 1.d0, 0.d0, kind(0.D0))
                  Op_V(nc,1)%O(2,2) = cmplx(-1.d0, 0.d0, kind(0.D0))
                  V_eff = -real(Bond_Matrix(1)%T(nb),Kind(0.d0)) *Bond_Matrix(1)%Prop_Fam(n_f) 
                  Op_V(nc, 1)%g = sqrt(cmplx(-Dtau*V_eff,0.d0,Kind(0.d0)))
                  Call Op_set(Op_V(nc,1))
                  Call Op_set(Op_V(nc,2))
                  ! Do-Do  
                  nc = nc + 1
                  do nf = 1,N_FL
                     call Op_make(Op_V(nc, nf), 2)
                     Op_V(nc, nf)%P(1) = I1
                     Op_V(nc, nf)%P(2) = J1
                     Op_V(nc, nf)%type = 2
                  enddo
                  Op_V(nc,2)%O(1,1) = cmplx( 1.d0, 0.d0, kind(0.D0))
                  Op_V(nc,2)%O(2,2) = cmplx(-1.d0, 0.d0, kind(0.D0))
                  V_eff = -real(Bond_Matrix(1)%T(nb),Kind(0.d0)) *Bond_Matrix(1)%Prop_Fam(n_f) 
                  Op_V(nc, 2)%g = sqrt(cmplx(-Dtau*V_eff,0.d0,Kind(0.d0)))
                  Call Op_set(Op_V(nc,1))
                  Call Op_set(Op_V(nc,2))
                  ! Up-do
                  nc = nc + 1
                  do nf = 1,N_FL
                     call Op_make(Op_V(nc, nf), 1)
                     Op_V(nc, nf)%type = 2
                  enddo
                  Op_V(nc,1)%P(1) = I1
                  Op_V(nc,1)%O(1,1) = cmplx( 1.d0, 0.d0, kind(0.D0))
                  V_eff = -real(Bond_Matrix(1)%T(nb),Kind(0.d0)) *Bond_Matrix(1)%Prop_Fam(n_f) 
                  Op_V(nc, 1)%g     = sqrt(cmplx(-Dtau*V_eff,0.d0,Kind(0.d0)))
                  Op_V(nc, 1)%alpha = cmplx(-0.5D0,0.d0,Kind(0.d0))
                  Op_V(nc,2)%P(1) = J1
                  Op_V(nc,2)%O(1,1) = cmplx(-1.d0, 0.d0, kind(0.D0))
                  V_eff = -real(Bond_Matrix(1)%T(nb),Kind(0.d0)) *Bond_Matrix(1)%Prop_Fam(n_f) 
                  Op_V(nc, 2)%g     = sqrt(cmplx(-Dtau*V_eff,0.d0,Kind(0.d0)))
                  Op_V(nc, 2)%alpha = cmplx( 0.5D0,0.d0,Kind(0.d0))
                  Call Op_set(Op_V(nc,1))
                  Call Op_set(Op_V(nc,2))
                  ! Do-up
                  nc = nc + 1
                  do nf = 1,N_FL
                     call Op_make(Op_V(nc, nf), 1)
                     Op_V(nc, nf)%type = 2
                  enddo
                  Op_V(nc,2)%P(1) = I1
                  Op_V(nc,2)%O(1,1) = cmplx( 1.d0, 0.d0, kind(0.D0))
                  V_eff = -real(Bond_Matrix(1)%T(nb),Kind(0.d0)) *Bond_Matrix(1)%Prop_Fam(n_f) 
                  Op_V(nc, 2)%g     = sqrt(cmplx(-Dtau*V_eff,0.d0,Kind(0.d0)))
                  Op_V(nc, 2)%alpha = cmplx(-0.5D0,0.d0,Kind(0.d0))
                  Op_V(nc,1)%P(1) = J1
                  Op_V(nc,1)%O(1,1) = cmplx(-1.d0, 0.d0, kind(0.D0))
                  V_eff = -real(Bond_Matrix(1)%T(nb),Kind(0.d0)) *Bond_Matrix(1)%Prop_Fam(n_f) 
                  Op_V(nc, 1)%g     = sqrt(cmplx(-Dtau*V_eff,0.d0,Kind(0.d0)))
                  Op_V(nc, 1)%alpha = cmplx( 0.5D0,0.d0,Kind(0.d0))
                  Call Op_set(Op_V(nc,1))
                  Call Op_set(Op_V(nc,2))
               enddo
            enddo
          else
            do nf = 1,N_FL
              do nc = N_hubbard + 1, Size(Op_V,1)
                  Call Op_make(Op_V(nc,nf), 2)
              enddo

              nc = N_hubbard
              Do n_f = 1, Bond_Matrix(1)%N_FAM
                Do l_f = 1, Bond_Matrix(1)%L_Fam(n_f)
                  I    = Bond_Matrix(1)%List_Fam(n_f,l_f,1)
                  nb   = Bond_Matrix(1)%List_Fam(n_f,l_f,2)
                  no_I = Bond_Matrix(1)%list(Nb,1)
                  no_J = Bond_Matrix(1)%list(Nb,2)
                  n_1  = Bond_Matrix(1)%list(Nb,3)
                  n_2  = Bond_Matrix(1)%list(Nb,4)
                  J    = Latt%nnlist(I,n_1,n_2)
                  I1   = Invlist(I,no_I)
                  J1   = Invlist(J,no_J)
                  nc = nc + 1
                  Op_V(nc,nf)%P(1) = I1
                  Op_V(nc,nf)%P(2) = J1
                  Op_V(nc,nf)%O(1,1) = cmplx(  1.d0, 0.d0, kind(0.D0))
                  Op_V(nc,nf)%O(2,2) = cmplx(  1.d0, 0.d0, kind(0.D0))
                  V_eff = real(Bond_Matrix(1)%T(nb),Kind(0.d0)) *Bond_Matrix(1)%Prop_Fam(n_f)
                  Op_V(nc,nf)%g = sqrt(cmplx(-Dtau*V_eff,0.d0,Kind(0.d0)))
                  Op_V(nc,nf)%alpha = cmplx(-1.d0,0.d0,Kind(0.d0))
                  Op_V(nc,nf)%type=2
                  !Write(6,*)  Op_V(nc,nf)%g
                  Call Op_set(Op_V(nc,nf))
                enddo
              enddo
            enddo
          endif
          
          Deallocate (Ham_V_vec, Ham_V2_vec, Ham_Vperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                   N_Phi_vec,  Ham_Lambda_vec )

          Deallocate (Ham_U_vec)

        end Subroutine Ham_V


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
        Subroutine  Alloc_obs(Ltau)

          Implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Nt
          Character (len=64) ::  Filename
          Character (len=2)  ::  Channel


          ! Scalar observables
          Allocate ( Obs_scal(4) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename = "Kin"
             case (2)
                N = 1;   Filename = "Pot"
             case (3)
                N = 1;   Filename = "Part"
             case (4)
                N = 1;   Filename = "Ener"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo

          ! Equal time correlators
          If ( Mz ) Then
             Allocate ( Obs_eq(5) )
             Do I = 1,Size(Obs_eq,1)
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
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo

             If (Ltau == 1) then
                ! Time-displaced correlators
                Allocate ( Obs_tau(5) )
                Do I = 1,Size(Obs_tau,1)
                   select case (I)
                   case (1)
                      Channel = 'P' ; Filename = "Green"
                   case (2)
                      Channel = 'PH'; Filename = "SpinZ"
                   case (3)
                      Channel = 'PH'; Filename = "SpinXY"
                   case (4)
                      Channel = 'PH'; Filename = "SpinT"
                   case (5)
                      Channel = 'PH'; Filename = "Den"
                   case default
                      Write(6,*) ' Error in Alloc_obs '
                   end select
                   Nt = Ltrot+1-2*Thtrot
                   If(Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
                enddo
             endif
          else
             ! Equal time correlators
             Allocate ( Obs_eq(3) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Filename = "Green"
                case (2)
                   Filename = "SpinZ"
                case (3)
                   Filename = "Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo

             If (Ltau == 1) then
                ! Time-displaced correlators
                Allocate ( Obs_tau(3) )
                Do I = 1,Size(Obs_tau,1)
                   select case (I)
                   case (1)
                      Channel = 'P' ; Filename = "Green"
                   case (2)
                      Channel = 'PH'; Filename = "SpinZ"
                   case (3)
                      Channel = 'PH'; Filename = "Den"
                   case default
                      Write(6,*) ' Error in Alloc_obs '
                   end select
                   Nt = Ltrot+1-2*Thtrot
                   If(Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
                enddo
             endif
          endif

        End Subroutine Alloc_obs

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
        subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, ZV, Z, ZP,ZS, ZZ, ZXY
          Integer :: I,J, imj, nf, dec, I1, I2, J1, no_I, no_J,n, l_f, n_1, n_2, n_f, nb 
          Real    (Kind=Kind(0.d0)) :: X

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight
          
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          ! Compute scalar observables.
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo


          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin* dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          dec = 1
          If ( Mz  ) dec = 2
          if ( Lattice_type == "Bilayer_square" .or. Lattice_type =="Bilayer_honeycomb" ) then
             Do I = 1,Latt%N
                do no_I = 1,Latt_unit%Norb/2
                   I1 = Invlist(I,no_I)
                   I2 = Invlist(I,no_I + Latt_unit%Norb/2) 
                   ZPot = ZPot + Ham_U * ( dble(N_SUN)*  (Grc(i1,i1,1) - cmplx(0.5,0.d0,kind(0.d0)) )**2 +  &
                           &               Grc(i1,i1,1) *  Gr(i1,i1,1)   )                               +  &      
                                 Ham_U2* ( dble(N_SUN)*  (Grc(i2,i2,1) - cmplx(0.5,0.d0,kind(0.d0)) )**2 +  &
                           &               Grc(i2,i2,1) *  Gr(i2,i2,1)   )                               
                enddo
             Enddo
          else
             Do I = 1,Latt%N
                do no_I = 1,Latt_unit%Norb
                   I1 = Invlist(I,no_I)
                   ZPot = ZPot +  Ham_U * ( dble(N_SUN)*  (Grc(i1,i1,1) - cmplx(0.5,0.d0,kind(0.d0)) )**2 +  &
                   &               Grc(i1,i1,1) *  Gr(i1,i1,1)   )   
                enddo
             Enddo
          Endif
          If (Mz) then 
            ! To Do
            Do n_f = 1, Bond_Matrix(1)%N_FAM
               Do l_f = 1, Bond_Matrix(1)%L_Fam(n_f)
                  I    = Bond_Matrix(1)%List_Fam(n_f,l_f,1)
                  nb     = Bond_Matrix(1)%List_Fam(n_f,l_f,2)
                  no_I = Bond_Matrix(1)%list(Nb,1)
                  no_J = Bond_Matrix(1)%list(Nb,2)
                  n_1  = Bond_Matrix(1)%list(Nb,3)
                  n_2  = Bond_Matrix(1)%list(Nb,4)
                  J    = Latt%nnlist(I,n_1,n_2)
                  ! sites of the bond
                  I1   = Invlist(I,no_I)
                  J1   = Invlist(J,no_J)
                  ! interation strength on bond
                  ZV  =  2.d0*Bond_Matrix(1)%T(Nb)
                  Z  =       ( GRC(I1,I1,1) + GRC(I1,I1,2) -  cmplx(1.d0,0.d0,Kind(0.d0)) ) *  &
                       &     ( GRC(J1,J1,1) + GRC(J1,J1,2) -  cmplx(1.d0,0.d0,Kind(0.d0)) ) +  &
                       &       GRC(I1,J1,1)*GR(I1,J1,1) +  GRC(I1,J1,2)*GR(I1,J1,2)             
                  Zpot  =  Zpot + ZV*Z
               enddo
            enddo
          else
            Do n_f = 1, Bond_Matrix(1)%N_FAM
               Do l_f = 1, Bond_Matrix(1)%L_Fam(n_f)
                  I    = Bond_Matrix(1)%List_Fam(n_f,l_f,1)
                  nb     = Bond_Matrix(1)%List_Fam(n_f,l_f,2)
                  no_I = Bond_Matrix(1)%list(Nb,1)
                  no_J = Bond_Matrix(1)%list(Nb,2)
                  n_1  = Bond_Matrix(1)%list(Nb,3)
                  n_2  = Bond_Matrix(1)%list(Nb,4)
                  J    = Latt%nnlist(I,n_1,n_2)
                  ! sites of the bond
                  I1   = Invlist(I,no_I)
                  J1   = Invlist(J,no_J)
                  ! interation strength on bond
                  ZV  =  Bond_Matrix(1)%T(Nb)
                  Z   =  ( dble(N_SUN) *  ( GRC(I1,I1,1) - cmplx(0.5,0.d0,Kind(0.d0) ))        + & 
                      &                   ( GRC(J1,J1,1) - cmplx(0.5,0.d0,Kind(0.d0) )) ) **2  + &
                      &   dble(N_SUN) * GRC(I1,J1,1)*GR(I1,J1,1) 
                  Zpot  =  Zpot  +  ZV*Z
               enddo
            enddo
          endif
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf)
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS

          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

          ! Standard two-point correlations
          If ( Mz ) then
             Call Predefined_Obs_eq_Green_measure  ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
             Call Predefined_Obs_eq_SpinMz_measure ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(2),Obs_eq(3),Obs_eq(4) )
             Call Predefined_Obs_eq_Den_measure    ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(5) )
          else
             Call Predefined_Obs_eq_Green_measure  ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
             Call Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(2) )
             Call Predefined_Obs_eq_Den_measure    ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(3) )
          endif


        end Subroutine Obser
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
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,  Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY
          Real    (Kind=Kind(0.d0)) :: X
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight

          ! Standard two-point correlations

          If ( Mz ) then
             Call Predefined_Obs_tau_Green_measure  ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(1) )
             Call Predefined_Obs_tau_SpinMz_measure ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(2),&
                  &                                   Obs_tau(3), Obs_tau(4) )
             Call Predefined_Obs_tau_Den_measure    ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(5) )
          Else
             Call Predefined_Obs_tau_Green_measure  ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(1) )
             Call Predefined_Obs_tau_SpinSUN_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(2) )
             Call Predefined_Obs_tau_Den_measure    ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(3) )
          endif

        end Subroutine OBSERT

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
      Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)
        Implicit none
        !> Operator index
        Integer, Intent(IN) :: n
        !> Time slice
        Integer, Intent(IN) :: nt
        !> New local field on time slice nt and operator index n
        Complex (Kind=Kind(0.d0)), Intent(In) :: Hs_new

        Integer :: nt1,I

        if (Continuous) then
           S0 = exp( (-Real(Hs_new,Kind(0.d0))**2  +  & 
                  &    Real(nsigma%f(n,nt),kind(0.d0))**2 ) /2.d0 ) 
        else
           S0 = 1.d0
        endif
      end function S0


!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!>   Forces_0  = \partial S_0 / \partial s  are calculated and returned to  main program.
!> 
!-------------------------------------------------------------------
        Subroutine Ham_Langevin_HMC_S0(Forces_0)

          Implicit none

          Real (Kind=Kind(0.d0)), Intent(inout), allocatable :: Forces_0(:,:)
          !Local
          Integer :: N, N_op,nt
          
          ! Compute \partial S_0 / \partial s
          N_op = size(nsigma%f,1)
          Forces_0  = 0.d0
          do n = 1,N_op
             if (OP_V(n,1)%type == 3 ) then
                do nt = 1,Ltrot
                   Forces_0(n,nt) = real(nsigma%f(n,nt), kind(0.d0))
                enddo
             endif
          enddo
          
        end Subroutine Ham_Langevin_HMC_S0
        
    end submodule ham_extended_Hubbard_smod
