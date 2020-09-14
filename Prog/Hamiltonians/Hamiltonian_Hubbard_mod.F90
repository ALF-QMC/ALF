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
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem
      real (Kind=Kind(0.d0)),        private :: ham_T2, ham_U2, ham_Tperp !  For Bilayers
      real (Kind=Kind(0.d0)),        private :: Phi_Y, Phi_X
      Integer               ,        private :: N_Phi
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: Checkerboard,  Bulk, Mz
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau


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
#endif
             File_Para = "parameters"
             File_info = "info"
#if defined(TEMPERING)
             write(File_para,'(A,I0,A)') "Temp_",igroup,"/parameters"
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif

#ifdef MPI
          If (Irank_g == 0 ) then
#endif
             OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'Ham_Set: unable to open <parameters>',ierr
                error stop 1
             END IF
             READ(5,NML=VAR_lattice)
             READ(5,NML=VAR_Model_Generic)
             READ(5,NML=VAR_Hubbard)
             CLOSE(5)

             Ltrot = nint(beta/dtau)
             if (Projector) Thtrot = nint(theta/dtau)
             Ltrot = Ltrot+2*Thtrot
             If ( Mz ) then
                N_FL  = 2
                if (mod(N_SUN,2) .ne. 0 ) then
                   Write(error_unit,*) 'Ham_Set: N_SUN has to be even if Mz = True'
                   error stop 1
                endif
                N_SUN = N_SUN / 2
             else
                N_FL  = 1
             endif

#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_SUN       ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_FL        ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_Phi       ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Phi_X       ,1  ,MPI_REAL8  ,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Phi_Y       ,1  ,MPI_REAL8  ,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Bulk        ,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
          CALL MPI_BCAST(Checkerboard,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Symm        ,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
          CALL MPI_BCAST(Ltrot       ,1,  MPI_INTEGER  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Thtrot      ,1,  MPI_INTEGER  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Projector   ,1,  MPI_LOGICAL  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau        ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta        ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_T       ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_chem    ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U       ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_T2      ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U2      ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_Tperp   ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Mz          ,1,  MPI_LOGICAL  , 0,Group_Comm,IERR)
#endif

          ! Setup the Bravais lattice
          Call  Ham_Latt

          ! Setup the hopping / single-particle part
          Call  Ham_Hop


          ! Setup the interaction.
          call Ham_V

#ifdef MPI
          If (Irank_g == 0) then
#endif
             OPEN(Unit = 50,file=file_info,status="unknown",position="append")
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
             If ( Mz )  then
                Write(50,*) 'HS  couples to z-component of spin'
             else
                Write(50,*) 'HS  couples to density'
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
             if ( Mz )  then
                Write(50,*) 'N_SUN         : ', 2*N_SUN
             else
                Write(50,*) 'N_SUN         : ',   N_SUN
             endif
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_U         : ', Ham_U
             Write(50,*) 'Ham_chem      : ', Ham_chem
             Close(50)
#ifdef MPI
          Endif
#endif
          ! Setup the trival wave function, in case of a projector approach
          if (Projector)   Call Ham_Trial(File_info)


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
        Subroutine Ham_Trial(file_info)


#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Use Predefined_Trial

          Implicit none
          Character (len=64), intent(in)  :: file_info


          Integer :: N_part, nf
#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)

          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif
          ! Use predefined stuctures or set your own Trial  wave function
          N_part = Ndim/2
          Call Predefined_TrialWaveFunction(Lattice_type ,Ndim,  List,Invlist,Latt, Latt_unit, &
               &                            N_part, N_FL,  WF_L, WF_R)


#ifdef MPI
          If (Irank_g == 0) then
#endif
             OPEN(Unit = 50,file=file_info,status="unknown",position="append")
             Do nf = 1,N_FL
                Write(50,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                Write(50,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
             enddo
             close(50)
#ifdef MPI
          endif
#endif

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
          
          Integer :: nf, I, I1, I2,  nc,  J, no,  N_ops 
          Real (Kind=Kind(0.d0)) :: X,  Zero = 1.D-10
          Real (Kind=Kind(0.d0)), allocatable :: Ham_U_vec(:)


          Allocate (Ham_U_vec(Latt_unit%Norb))
          
          N_ops = 0
          if ( Lattice_type == "Bilayer_square" .or. Lattice_type =="Bilayer_honeycomb" ) then
             Ham_U_vec(1) = Ham_U
             Ham_U_vec(2) = Ham_U2
             If (abs(Ham_U ) > Zero ) N_ops = N_ops + Latt%N*Latt_unit%Norb/2
             If (abs(Ham_U2) > Zero ) N_ops = N_ops + Latt%N*Latt_unit%Norb/2
          else
             Ham_U_vec(:)  = Ham_U
             If (abs(Ham_U ) > Zero ) N_ops = N_ops + Latt%N*Latt_unit%Norb
          endif
          If ( Mz )  Then
             Allocate(Op_V(N_ops,N_FL))
             Ham_U_vec = Ham_U_vec/real(N_SUN,kind(0.d0))
             nc = 0
             Do I1 = 1,Latt%N
                do no = 1, Latt_unit%Norb
                   I = invlist(I1,no)
                   if (abs(Ham_U_vec(no)) > Zero ) then
                      nc = nc + 1
                      Call Predefined_Int_U_MZ ( OP_V(nc,1), OP_V(nc,2), I,  DTAU, Ham_U_vec(no) )
                   endif
                enddo
             enddo
          else
             Allocate(Op_V(Ndim,N_FL))
             nc = 0
             Do I1 = 1,Latt%N
                do no = 1, Latt_unit%Norb
                   I = invlist(I1,no)
                   if (abs(Ham_U_vec(no)) > Zero ) then
                      nc = nc + 1
                      Call Predefined_Int_U_SUN(  OP_V(nc,1), I, N_SUN, DTAU, Ham_U_vec(no)  )
                   endif
                Enddo
             Enddo
          Endif

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
          Integer    ::  i, N, Ns,Nt,No, Norb
          Character (len=64) ::  Filename


          Norb = Latt_unit%Norb
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
          If ( Mz ) Then
             Allocate ( Obs_eq(5) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Ns = Latt%N;  No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinT"
                case (5)
                   Ns = Latt%N;  No = Norb;  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
             enddo

             If (Ltau == 1) then
                ! Equal time correlators
                Allocate ( Obs_tau(5) )
                Do I = 1,Size(Obs_tau,1)
                   select case (I)
                   case (1)
                      Ns = Latt%N; No = Norb;  Filename ="Green"
                   case (2)
                      Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                   case (3)
                      Ns = Latt%N; No = Norb;  Filename ="SpinXY"
                   case (4)
                      Ns = Latt%N; No = Norb;  Filename ="SpinT"
                   case (5)
                      Ns = Latt%N; No = Norb;  Filename ="Den"
                   case default
                      Write(6,*) ' Error in Alloc_obs '
                   end select
                   Nt = Ltrot+1-2*Thtrot
                   Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
                enddo
             endif
          else
             ! Equal time correlators
             Allocate ( Obs_eq(3) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Ns = Latt%N;  No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N;  No = Norb;  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
             enddo

             If (Ltau == 1) then
                ! Equal time correlators
                Allocate ( Obs_tau(3) )
                Do I = 1,Size(Obs_tau,1)
                   select case (I)
                   case (1)
                      Ns = Latt%N; No = Norb;  Filename ="Green"
                   case (2)
                      Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                   case (3)
                      Ns = Latt%N; No = Norb;  Filename ="Den"
                   case default
                      Write(6,*) ' Error in Alloc_obs '
                   end select
                   Nt = Ltrot+1-2*Thtrot
                   Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
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
          dec = 1
          If ( Mz  ) dec = 2
          Do I = 1,Latt%N
             do no_I = 1,Latt_unit%Norb
                I1 = Invlist(I,no_I)
                if (no_I == 1)  ZPot = ZPot + Grc(i1,i1,1) * Grc(i1,i1, dec)* ham_U
                if (no_I == 2)  ZPot = ZPot + Grc(i1,i1,1) * Grc(i1,i1, dec)* ham_U2
             enddo
          Enddo
          
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

#include "Hamiltonian_Hubbard_include.h"


    end Module Hamiltonian
