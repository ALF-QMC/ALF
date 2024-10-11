!  Copyright (C) 2016 - 2019 The ALF project
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

    Module Hamiltonian

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
      use iso_fortran_env, only: output_unit, error_unit


      Implicit none


      Type (Operator),     dimension(:,:), allocatable :: Op_V
      Type (Operator),     dimension(:,:), allocatable :: Op_T
      Type (WaveFunction), dimension(:),   allocatable :: WF_L
      Type (WaveFunction), dimension(:),   allocatable :: WF_R
      Type (Fields)        :: nsigma
      Integer              :: Ndim
      Integer              :: N_FL
      Integer              :: N_SUN
      Integer              :: Ltrot
      Integer              :: Thtrot
      Logical              :: Projector
      Integer              :: Group_Comm
      Logical              :: Symm


      Type (Lattice),       private :: Latt
      Type (Unit_cell),     private :: Latt_unit
      Integer,              private :: L1, L2
      Type (Hopping_Matrix_type), Allocatable, private :: Hopping_Matrix(:)
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem, Ham_h, Ham_J, Ham_xi,  Ham_tV
      real (Kind=Kind(0.d0)),        private :: ham_T2, ham_U2, ham_Tperp !  For Bilayers
      real (Kind=Kind(0.d0)),        private :: ham_alpha, Percent_change
      Real (Kind=Kind(0.d0)),        private :: Phi_Y, Phi_X
      Integer               ,        private :: N_Phi
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: Checkerboard,  Bulk,  Mz
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau

!>    Storage for the Ising action
      Real (Kind=Kind(0.d0)),  private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
      Integer,  allocatable ,  private :: L_bond(:,:), L_bond_inv(:,:), Ising_nnlist(:,:)


    contains

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

          integer                :: ierr, N_part, nf
          Character (len=64)     :: file_info, file_para


          ! L1, L2, Lattice_type, List(:,:), Invlist(:,:) -->  Lattice information
          ! Ham_T, Chem, Phi_X, XB_B, Checkerboard, Symm   -->  Hopping
          ! Interaction                              -->  Model
          ! Simulation type                          -->  Finite  T or Projection  Symmetrize Trotter.

          NAMELIST /VAR_Lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_Model_Generic/  Checkerboard, N_SUN, N_FL, Phi_X, Phi_Y, Symm, Bulk, N_Phi, Dtau, Beta, Theta, Projector

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U, ham_T2, ham_U2, ham_Tperp,  Mz

          Namelist /VAR_LRC/      ham_T, ham_T2, ham_Tperp, ham_chem, ham_U, ham_alpha, Percent_change

          Namelist /VAR_Ising/    ham_T, ham_chem, ham_U, Ham_h, Ham_J, Ham_xi

          NAMELIST /VAR_t_V/      Ham_T, ham_chem, ham_tV

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
#endif
          ! Global "Default" values.
          N_SUN        = 1
          Checkerboard = .false.
          Symm         = .false.
          Projector    = .false.
          Bulk         = .true.
          Phi_X        = 0.d0
          Phi_Y        = 0.d0
          N_Phi        = 0
          Ham_T2       = 0.d0
          Ham_Tperp    = 0.d0
          Ham_U2       = 0.d0


#ifdef MPI
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
          If (Irank_g == 0 ) then
#endif
             File_para = "parameters"
             File_info = "info"
#if defined(TEMPERING)
             write(File_para,'(A,I0,A)') "Temp_",igroup,"/parameters"
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif

             OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
             OPEN(Unit = 50,file=file_info,status="unknown",position="append")
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <parameters>',ierr
                error stop 1
             END IF
             READ(5,NML=VAR_lattice)
             READ(5,NML=VAR_Model_Generic)
             Ltrot = nint(beta/dtau)
             Thtrot = 0
             if (Projector) Thtrot = nint(theta/dtau)
             Ltrot = Ltrot+2*Thtrot


#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_SUN       ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_FL        ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_Phi       ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Phi_X       ,1  ,MPI_REAL8  ,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Phi_Y       ,1  ,MPI_REAL8  ,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
          CALL MPI_BCAST(Checkerboard,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Symm        ,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Bulk        ,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
          CALL MPI_BCAST(Ltrot       ,1,  MPI_INTEGER  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Thtrot      ,1,  MPI_INTEGER  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Projector   ,1,  MPI_LOGICAL  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau        ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta        ,1,  MPI_REAL8    , 0,Group_Comm,ierr)


#endif

          Call  Ham_Latt

#ifdef MPI
          If (Irank_g == 0) then
#endif
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) '# of orbitals : ', Ndim
             Write(50,*) 'Flux_1        : ', Phi_X
             Write(50,*) 'Flux_2        : ', Phi_Y
             If (Bulk) then
                Write(50,*) 'Twist as phase factor in bulk'
             Else
                Write(50,*) 'Twist as boundary condition'
             endif
             Write(50,*) 'Checkerboard  : ', Checkerboard
             Write(50,*) 'Symm. decomp  : ', Symm
             if (Projector) then
                Write(50,*) 'Projective version'
                Write(50,*) 'Theta         : ', Theta
                Write(50,*) 'Tau_max       : ', beta
             else
                Write(50,*) 'Finite temperture version'
                Write(50,*) 'Beta          : ', Beta
             endif
             Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
#ifdef MPI
          Endif
#endif



          ! Default is finite temperature.
          Select Case (Model)
          Case ("LRC")
             N_SUN = 2
             N_FL  = 1
#ifdef MPI
             If (Irank_g == 0 ) then
#endif
                READ(5,NML=VAR_LRC)
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_alpha     : ', Ham_alpha
                Write(50,*) 'Percent_change: ', Percent_change
                Write(50,*) 'Ham_chem      : ', Ham_chem
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(ham_T         ,1,MPI_REAL8  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_chem      ,1,MPI_REAL8  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_U         ,1,MPI_REAL8  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_alpha     ,1,MPI_REAL8  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(Percent_change,1,MPI_REAL8  ,0,Group_Comm,ierr)
#endif

          Case ("Hubbard")
#ifdef MPI
             If (Irank_g == 0 ) then
#endif
                READ(5,NML=VAR_Hubbard)
                If ( Mz ) then
                   N_FL  = 2
                   N_SUN = 1
                Else
                   N_FL=1
                Endif
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
                Write(50,*) 'Mz            : ', Mz
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(N_SUN    ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
             CALL MPI_BCAST(N_FL     ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_T    ,1 ,MPI_REAL8   ,   0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_chem ,1 ,MPI_REAL8   ,   0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_U    ,1 ,MPI_REAL8   ,   0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_T2   ,1 ,MPI_REAL8   ,   0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_U2   ,1 ,MPI_REAL8   ,   0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_Tperp,1 ,MPI_REAL8   ,   0,Group_Comm,ierr)
             CALL MPI_BCAST(Mz       ,1 ,MPI_LOGICAL ,   0,Group_Comm,IERR)

#endif
          Case ("Hubbard_SU2_Ising")
             N_FL = 1
             N_SUN = 2
#ifdef MPI
             If (Irank_g == 0 ) then
#endif
                READ(5,NML=VAR_Ising)
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
                Write(50,*) 'Ham_xi        : ', Ham_xi
                Write(50,*) 'Ham_J         : ', Ham_J
                Write(50,*) 'Ham_h         : ', Ham_h
#ifdef MPI
             Endif
#endif
             If ( Lattice_type == "Honeycomb" .or. Lattice_type == "Pi_Flux" ) then
                Write(error_unit,*) "Hubbard_SU2_Ising is only implemented for a square lattice"
                error stop 1
             Endif
#ifdef MPI
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(Ham_xi   ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,Group_Comm,ierr)
#endif
             Call Setup_Ising_action
          Case ("t_V")
             N_FL  = 1
#ifdef MPI
             If (Irank_g == 0 ) then
#endif
                READ(5,NML=VAR_t_V )
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_chem      : ', Ham_chem
                Write(50,*) 'Ham_tV         : ', Ham_tV
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_tV   ,1,MPI_REAL8,0,Group_Comm,ierr)
#endif
          Case default
             Write(error_unit,*) "Model not yet implemented!"
             error stop 1
          end Select

          If (Model == "Hubbard" )  then
             If ( Mz )  then
                Model = "Hubbard_Mz"
             else
                Model = "Hubbard_SU2"
             endif
          Endif
          Call  Ham_Hop


          if (Projector)  Then
             Call Ham_Trial
#ifdef MPI
             If (Irank_g == 0) then
#endif
                Do nf = 1,N_FL
                   Write(50,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(50,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
#ifdef MPI
             Endif
#endif
          Endif

#ifdef MPI
          If (Irank_g == 0 )  then
#endif
             close(50)
             Close(5)
#ifdef MPI
          endif
#endif

! #ifdef MPI
!           If (Irank == 0 )  then
! #endif
! #if defined(STAB1)
!              Write(50,*) 'STAB1 is defined '
! #endif
!
! #if defined(QRREF)
!              Write(50,*) 'QRREF is defined '
! #endif
!              close(50)
!              Close(5)
! #ifdef MPI
!           endif
! #endif

          If (Model == "LRC" ) Then
             Call LRC_Set_VIJ(Latt, Latt_unit, Ham_U, Ham_alpha, list, invlist)
#ifdef MPI
             If (Irank_g == 0 )  then
#endif
                Call LRC_Print(Latt, Latt_unit, list, invlist)
#ifdef MPI
             Endif
#endif

          Endif

          call Ham_V


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

          Select case (Lattice_type)
          Case ("Square")
             Call  Set_Default_hopping_parameters_square(Hopping_Matrix,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("N_leg_ladder")
             Call  Set_Default_hopping_parameters_n_leg_ladder(Hopping_Matrix, Ham_T_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
                  &                                            Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("Honeycomb")
             Ham_Lambda = 0.d0
             Call  Set_Default_hopping_parameters_honeycomb(Hopping_Matrix, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                         Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("Bilayer_square")
             Call  Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
                  &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                  &                                              List, Invlist, Latt, Latt_unit )

          Case ("Bilayer_honeycomb")
             Call  Set_Default_hopping_parameters_Bilayer_honeycomb(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
                  &                                                 Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                  &                                                 List, Invlist, Latt, Latt_unit )

          end Select

          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_T )

          !Call Predefined_Hopping(Lattice_type, Ndim, List,Invlist, Latt, Latt_unit, &
          !     &                      Dtau, Ham_T, Ham_Chem,  Phi_X, Phi_Y, Bulk, N_Phi,   &
          !      &                      N_FL,  Checkerboard, Symm, OP_T )
          !Do Nth = 1,1
          !   Do n = 1, size(OP_T(1,1)%E)
          !      Write(31,"(I4,2x,F14.7)") n_Phi, OP_T(1,1)%E(n)
          !   enddo
          !
          !   Call Op_clear (OP_T(1,1),N)
          !   Deallocate (OP_T)
          !Enddo
          !Stop

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
        Subroutine Ham_Trial


          Use Predefined_Trial

          Implicit none
          Integer :: N_part

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

          Integer :: nf, I, I1, I2,  nc, nc1,  J, no
          Real (Kind=Kind(0.d0)) :: X
          Real (Kind=Kind(0.d0)), allocatable :: Ham_U_vec(:)


          Allocate (Ham_U_vec(Latt_unit%Norb))
          Ham_U_vec(:)  = Ham_U
          if ( Lattice_type == "Bilayer_square" .or. Lattice_type =="Bilayer_honeycomb" ) then
             Ham_U_vec(1) = Ham_U
             Ham_U_vec(2) = Ham_U2
          endif


          Select case (Model)
          Case ("LRC")
             Allocate(Op_V(Ndim,N_FL))
             nc = 0
             Do i = 1,Ndim
                nc = nc + 1
                Call Predefined_Int_LRC( OP_V(nc,1), I, DTAU  )
             Enddo
          Case ("Hubbard_SU2")
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Ndim,N_FL))
             Do I1 = 1,Latt%N
                do no = 1, Latt_unit%Norb
                   I = invlist(I1,no)
                   Call Predefined_Int_U_SUN(  OP_V(I,1), I, N_SUN, DTAU, Ham_U_vec(no)  )
                enddo
             enddo
          Case ("Hubbard_Mz")
             Allocate(Op_V(Ndim,N_FL))
             Do I1 = 1,Latt%N
                do no = 1, Latt_unit%Norb
                    I = invlist(I1,no)
                    Call Predefined_Int_U_MZ ( OP_V(I,1), OP_V(I,2), I,  DTAU, Ham_U_vec(no) )
                 enddo
              enddo
           Case ("Hubbard_SU2_Ising")
              Allocate(Op_V(3*Ndim,N_FL))

             Do nc = 1,Ndim*Latt_unit%N_coord   ! Runs over bonds.  Coordination number = 2.
                ! For the square lattice Ndim = Latt%N
                I1 = L_bond_inv(nc,1)
                I2 = I1
                if ( L_bond_inv(nc,2)  == 1 ) I2 = Latt%nnlist(I1,1,0)
                if ( L_bond_inv(nc,2)  == 2 ) I2 = Latt%nnlist(I1,0,1)
                Call Predefined_Int_Ising_SUN( OP_V(nc,1), I1, I2, DTAU, Ham_xi )
             Enddo

             Do i = 1,Ndim
                nc1 = Latt_unit%N_coord*Ndim + i
                Call Predefined_Int_U_SUN(  OP_V(nc1,1), I, N_SUN, DTAU, Ham_U  )
             Enddo
          case ("t_V")
             Allocate(Op_V(Latt_unit%N_coord*Latt%N,1))
             select case (Lattice_type)
             case ("Square")
                !Write(6,*) "N_coord, Latt%N", N_coord, Latt%N, Dtau, ham_tV
                nc = 0
                do I = 1,Latt%N
                   I1 = I
                   do nc1 = 1,Latt_unit%N_coord
                      nc = nc + 1
                      if (nc1 == 1 )  I2 = latt%nnlist(I,1,0)
                      if (nc1 == 2 )  I2 = latt%nnlist(I,0,1)
                      Call Predefined_Int_V_SUN( OP_V(nc,1), I1, I2, N_SUN, DTAU, Ham_tV )
                   Enddo
                Enddo
             case ("Pi_Flux")
                nc = 0
                do I = 1,Latt%N
                   I1 = Invlist(I,1)
                   do nc1 = 1,Latt_unit%N_coord
                      nc = nc + 1
                      If (nc1 == 1 )  I2 = invlist(I,2)
                      If (nc1 == 2 )  I2 = invlist(Latt%nnlist(I,0, 1),2)
                      If (nc1 == 3 )  I2 = invlist(Latt%nnlist(I,-1,1),2)
                      If (nc1 == 4 )  I2 = invlist(Latt%nnlist(I,-1,0),2)
                      Call Predefined_Int_V_SUN( OP_V(nc,1), I1, I2, N_SUN, DTAU, Ham_tV )
                   Enddo
                Enddo
             case default
                Write(error_unit,*) "Lattice for t-V  not implemented"
                error stop 1
             End select
          case default
             Write(error_unit,*) "Model not yet implemented!"
             error stop 1
          end Select

          Deallocate (Ham_U_vec)

        end Subroutine Ham_V

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
          Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new

          Integer :: nt1,I
          !Write(6,*) "Hi1"

          S0 = 1.d0
          If ( Op_V(n,1)%type == 1 ) then
             do i = 1,4
                S0 = S0*DW_Ising_space(nsigma%i(n,nt)*nsigma%i(Ising_nnlist(n,i),nt))
             enddo
             nt1 = nt +1
             if (nt1 > Ltrot) nt1 = 1
             S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
             nt1 = nt - 1
             if (nt1 < 1  ) nt1 = Ltrot
             S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
             If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
          endif

          If (model == "LRC" ) then
             !Write(6,*) "Hi2"
             S0 = LRC_S0(n,dtau,nsigma%f(:,nt),Hs_new,N_SUN)
          Endif

        end function S0
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Setup  bonds and median lattice.
!> @details
!> Note that the median lattice on which the Ising bond fields are defined could be defined generally in
!> the predefined lattice subroutine.
!--------------------------------------------------------------------

        Subroutine Setup_Ising_action

          ! This subroutine sets up lists and arrays so as to enable an
          ! an efficient calculation of  S0(n,nt)

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation
          Real (Kind=Kind(0.d0)) :: X_p(2)

          ! Setup list of bonds for the square lattice.
          Allocate (L_Bond(Latt%N,2),  L_bond_inv(Latt%N*Latt_unit%N_coord,2) )

          nc = 0
          do nth = 1,2*Latt_unit%N_coord
             Do n1= 1, L1/2
                Do n2 = 1,L2
                   nc = nc + 1
                   n_orientation = 1
                   I1 = 1
                   If (nth == 1 ) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 2) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p  + latt%a1_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 3) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   elseif (nth == 4) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p  + latt%a2_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   endif
                   L_bond(I1,n_orientation) = nc
                   L_bond_inv(nc,1) = I1
                   L_bond_inv(nc,2) = n_orientation
                   ! The bond is given by  I1, I1 + a_(n_orientation).
                Enddo
             Enddo
          Enddo
          ! Setup the nearest neigbour lists for the Ising spins.
          allocate(Ising_nnlist(2*Latt%N,4))
          do I  = 1,Latt%N
             n  = L_bond(I,1)
             n1 = L_bond(Latt%nnlist(I, 1, 0),2)
             n2 = L_bond(Latt%nnlist(I, 0, 0),2)
             n3 = L_bond(Latt%nnlist(I, 0,-1),2)
             n4 = L_bond(Latt%nnlist(I, 1,-1),2)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
             n  = L_bond(I,2)
             n1 = L_bond(Latt%nnlist(I, 0, 1),1)
             n2 = L_bond(Latt%nnlist(I,-1, 1),1)
             n3 = L_bond(Latt%nnlist(I,-1, 0),1)
             n4 = L_bond(Latt%nnlist(I, 0, 0),1)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
          enddo
          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J)
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J)

        End Subroutine Setup_Ising_action
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
          Character (len=:), allocatable ::  Channel


          ! Scalar observables
          Allocate ( Obs_scal(4) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename ="Kin"
             case (2)
                N = 1;   Filename ="Pot"
             case (3)
                N = 1;   Filename ="Part"
             case (4)
                N = 1;   Filename ="Ener"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo

          ! Equal time correlators
          If (Model ==  "Hubbard_Mz") Then
             Allocate ( Obs_eq(5) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Filename ="Green"
                case (2)
                   Filename ="SpinZ"
                case (3)
                   Filename ="SpinXY"
                case (4)
                   Filename ="SpinT"
                case (5)
                   Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo

             If (Ltau == 1) then
                ! Equal time correlators
                Allocate ( Obs_tau(5) )
                Do I = 1,Size(Obs_tau,1)
                   select case (I)
                   case (1)
                      Channel = 'P' ; Filename ="Green"
                   case (2)
                      Channel = 'PH'; Filename ="SpinZ"
                   case (3)
                      Channel = 'PH'; Filename ="SpinXY"
                   case (4)
                      Channel = 'PH'; Filename ="SpinT"
                   case (5)
                      Channel = 'PH'; Filename ="Den"
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
                   Filename ="Green"
                case (2)
                   Filename ="SpinZ"
                case (3)
                   Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo

             If (Ltau == 1) then
                ! Equal time correlators
                Allocate ( Obs_tau(3) )
                Do I = 1,Size(Obs_tau,1)
                   select case (I)
                   case (1)
                      Channel = 'P' ; Filename ="Green"
                   case (2)
                      Channel = 'PH'; Filename ="SpinZ"
                   case (3)
                      Channel = 'PH'; Filename ="Den"
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
!> Global moves
!>
!> @details
!>  This routine generates a
!>  global update  and returns the propability T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!> @param [OUT]  T0_Proposal_ratio Real
!> \verbatimam
!>  T0_Proposal_ratio  =  T0( sigma_new -> sigma_old ) /  T0( sigma_old -> sigma_new)
!> \endverbatim
!> @param [OUT]  Size_clust Real
!> \verbatim
!>  Size of cluster that will be flipped.
!> \endverbatim
!-------------------------------------------------------------------
        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)

          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Type (Fields),  Intent(IN)  :: nsigma_old

          ! Local
          Integer :: N_op, N_tau, n1,n2, n

          T0_Proposal_ratio  = 1.d0
          size_clust         = 3.d0

          nsigma%f = nsigma_old%f
          nsigma%t = nsigma_old%t
          N_op  = size(nsigma_old%f,1)
          N_tau = size(nsigma_old%f,2)
          Do n  = 1, Nint(size_clust)
             n1 = nranf(N_op)
             n2 = nranf(N_tau)
             nsigma%f(n1,n2) = nsigma_old%flip(n1,n2)
          enddo

        End Subroutine Global_move
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
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none

          ! Arguments
          Type (Fields),  INTENT(IN) :: nsigma_old

          ! Local
          Integer :: I,n,n1,n2,n3,n4,nt,nt1, nc_F, nc_J, nc_h_p, nc_h_m


          Delta_S0_global = 1.d0
          If ( Model == "Hubbard_SU2_Ising" ) then
             nc_F = 0
             nc_J = 0
             nc_h_p = 0
             nc_h_m = 0
             Do I = 1,Latt%N
                n1  = L_bond(I,1)
                n2  = L_bond(Latt%nnlist(I,1,0),2)
                n3  = L_bond(Latt%nnlist(I,0,1),1)
                n4  = L_bond(I,2)
                do nt = 1,Ltrot
                   nt1 = nt +1
                   if (nt == Ltrot) nt1 = 1
                   if (nsigma%i(n1,nt) == nsigma%i(n1,nt1) ) then
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n1,nt) == nsigma_old%i(n1,nt1) ) then
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   if (nsigma%i(n4,nt) == nsigma%i(n4,nt1) ) then
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n4,nt) == nsigma_old%i(n4,nt1) ) then
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   nc_F = nc_F + nsigma%i    (n1,nt)*nsigma%i    (n2,nt)*nsigma%i    (n3,nt)*nsigma%i  (n4,nt)  &
                        &      - nsigma_old%i(n1,nt)*nsigma_old%i(n2,nt)*nsigma_old%i(n3,nt)*nsigma_old%i(n4,nt)

                   nc_J = nc_J + nsigma%i(n1,nt)*nsigma%i(n2,nt) + &
                        &        nsigma%i(n2,nt)*nsigma%i(n3,nt) + &
                        &        nsigma%i(n3,nt)*nsigma%i(n4,nt) + &
                        &        nsigma%i(n4,nt)*nsigma%i(n1,nt) - &
                        &        nsigma_old%i(n1,nt)*nsigma_old%i(n2,nt) - &
                        &        nsigma_old%i(n2,nt)*nsigma_old%i(n3,nt) - &
                        &        nsigma_old%i(n3,nt)*nsigma_old%i(n4,nt) - &
                        &        nsigma_old%i(n4,nt)*nsigma_old%i(n1,nt)

                enddo
             enddo
             !             Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) * &
             !                  &            exp( -Dtau*(Ham_F*real(nc_F,kind(0.d0)) -  Ham_J*real(nc_J,kind(0.d0))))
             ! No flux in example code. May want to include it.
             Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) * &
                  &            exp( Dtau* Ham_J*real(nc_J,kind(0.d0)))
             !Write(6,*) Delta_S0_global
          endif


        end Function Delta_S0_global

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
        subroutine Obser(GR,Phase,Ntau)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J,n
          Real    (Kind=Kind(0.d0)) :: X

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))


          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >

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
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" ) then
             dec = 1
             Do I = 1,Ndim
                ZPot = ZPot + Grc(i,i,1) * Grc(i,i, dec)
             Enddo
             Zpot = Zpot*ham_U
          elseif (Model == "Hubbard_Mz") then
             dec = 2
             Do I = 1,Ndim
                ZPot = ZPot + Grc(i,i,1) * Grc(i,i, dec)
             Enddo
             Zpot = Zpot*ham_U
          elseif (Model == "LRC") then
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             Do I = 1,Ndim
                ZPot = ZPot +    LRC_V_int(I,I)* Grc(i,i, 1)* Grc(i,i,1)
                Do J = I+1,Ndim
                   ZPot = ZPot + Z*LRC_V_int(I,J)*(Z *  Grc(i,i, 1)* Grc(j,j,1)  + Grc(i,j,1)*Gr(i,j,1) )
                Enddo
             Enddo
          Endif
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
          If (Model == "Hubbard_Mz" ) then
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
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase

          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY
          Real    (Kind=Kind(0.d0)) :: X
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ! Standard two-point correlations

          If (Model  == "Hubbard_Mz" ) then
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
!> Prints out the bins.  No need to change this routine.
!-------------------------------------------------------------------
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau

          !Local
          Integer :: I


          Do I = 1,Size(Obs_scal,1)
             Call  Print_bin_Vec(Obs_scal(I),Group_Comm)
          enddo
          Do I = 1,Size(Obs_eq,1)
             Call  Print_bin_Latt(Obs_eq(I),Group_Comm)
          enddo
          If (Ltau  == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call  Print_bin_Latt(Obs_tau(I),Group_Comm)
             enddo
          endif

        end Subroutine Pr_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Initializes observables to zero before each bins.  No need to change
!> this routine.
!-------------------------------------------------------------------
        Subroutine  Init_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau

          ! Local
          Integer :: I

          Do I = 1,Size(Obs_scal,1)
             Call Obser_vec_Init(Obs_scal(I))
          Enddo

          Do I = 1,Size(Obs_eq,1)
             Call Obser_Latt_Init(Obs_eq(I))
          Enddo

          If (Ltau == 1) then
             Do I = 1,Size(Obs_tau,1)
                Call Obser_Latt_Init(Obs_tau(I))
             Enddo
          Endif

        end Subroutine Init_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specify a global move on a given time slice tau.
!>
!> @details
!> @param[in] ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!> @param[out] T0_Proposal_ratio, Real
!> \verbatim
!>  T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)
!> \endverbatim
!> @param[out] S0_ratio, Real
!> \verbatim
!>  S0_ratio = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
!> \endverbatim
!> @param[out] Flip_length  Integer
!> \verbatim
!>  Number of flips stored in the first  Flip_length entries of the array Flip_values.
!>  Has to be smaller than NDIM
!> \endverbatim
!> @param[out] Flip_list  Integer(Ndim)
!> \verbatim
!>  List of spins to be flipped: nsigma%f(Flip_list(1),ntau) ... nsigma%f(Flip_list(Flip_Length),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!> @param[out] Flip_value  Real(Ndim)
!> \verbatim
!>  Flip_value(:)= nsigma%flip(Flip_list(:),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!--------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)


          Implicit none
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: T0_Proposal_ratio,  S0_ratio
          Integer                , INTENT(OUT) :: Flip_list(:)
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: Flip_value(:)
          Integer, INTENT(OUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          ! Local
          Integer :: n_op, n, ns
          Real (Kind=Kind(0.d0)) :: T0_proposal

          If (Model == "LRC" ) then
             Call LRC_draw_field(Percent_change, Dtau, nsigma%f(:,ntau), Flip_value,N_SUN)
             Do n = 1,Ndim
                Flip_list(n) = n
                !Write(6,*) Flip_value(n), nsigma%f(n,ntau)
             Enddo
             !Write(6,*)
             Flip_length    = Ndim
             ! T0_Proposal_ration exactly cancels S0_ratio
             T0_Proposal_ratio = 1.d0
             S0_ratio          = 1.d0
          else
             Flip_length = nranf(4)
             do n = 1,flip_length
                n_op = nranf(size(OP_V,1))
                Flip_list(n)  = n_op
                Flip_value(n) = nsigma%flip(n_op,ntau)
                If ( OP_V(n_op,1)%type == 1 ) then
                   S0_ratio          =   S0(n_op,ntau,Flip_value(n))
                   T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0_ratio) ! No move prob
                   If ( T0_Proposal > Ranf_wrap() ) then
                      T0_Proposal_ratio =  1.d0 / S0_ratio
                   else
                      T0_Proposal_ratio = 0.d0
                   endif
                else
                   T0_Proposal_ratio = 1.d0
                   S0_ratio          = 1.d0
                endif
             Enddo
          endif

        end Subroutine Global_move_tau

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> The user can set the initial field.
!>
!> @details
!> @param[OUT] Initial_field Real(:,:)
!> \verbatim
!>  Upon entry Initial_field is not allocated. If alloacted then it will contain the
!>  the initial field
!> \endverbatim
!--------------------------------------------------------------------
     Subroutine  Hamiltonian_set_nsigma(Initial_field)
        Implicit none

        Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(OUT) :: Initial_field


      end Subroutine Hamiltonian_set_nsigma

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This routine allows to user to  determine the global_tau sampling parameters at run time
!> It is especially usefull if these parameters are dependent on other parameters.
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Subroutine Overide_global_tau_sampling_parameters(Nt_sequential_start,Nt_sequential_end,N_Global_tau)

        Implicit none
        Integer, Intent(INOUT) :: Nt_sequential_start,Nt_sequential_end, N_Global_tau

        If ( Model  == "LRC" )  then
           Nt_sequential_start = 1
           Nt_sequential_end   = 0
           N_Global_tau   = Nint(1.d0/Percent_change)
        endif

      end Subroutine Overide_global_tau_sampling_parameters

      end Module Hamiltonian
