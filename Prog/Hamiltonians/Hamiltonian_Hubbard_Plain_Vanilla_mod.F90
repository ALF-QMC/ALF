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
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Integer               ,        private :: N_part
      Character (len=64),   private :: Model, Lattice_type
     

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

          integer                :: ierr, nf
          Character (len=64)     :: file_info, file_para
          
          
          
          NAMELIST /VAR_Lattice/  L1, L2, Lattice_type, Model


          NAMELIST /VAR_Hubbard_Plain_Vanilla/  ham_T, ham_chem, ham_U, Dtau, Beta, Projector, Theta, Symm, N_part
          
          

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
#endif
          ! Global "Default" values.

          
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
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
             N_part =  L1*L2/2
             If (L1 == 1) then
                Write(6,*) 'For  one-dimensional lattices set L2=1'
                stop
             endif
             READ(5,NML=VAR_Hubbard_Plain_Vanilla)
             CLOSE(5)

             Ltrot = nint(beta/dtau)
             if (Projector) Thtrot = nint(theta/dtau)
             Ltrot = Ltrot+2*Thtrot
             N_SUN        = 1
             N_FL         = 2
          
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_SUN       ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_FL        ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
          CALL MPI_BCAST(Symm        ,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
          CALL MPI_BCAST(Ltrot       ,1,  MPI_INTEGER  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(N_part      ,1,  MPI_INTEGER  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Thtrot      ,1,  MPI_INTEGER  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Projector   ,1,  MPI_LOGICAL  , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau        ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta        ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_T       ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_chem    ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U       ,1,  MPI_REAL8    , 0,Group_Comm,ierr)
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
             Write(50,*) 'L1            : ', L1
             Write(50,*) 'L2            : ', L2
             Write(50,*) '# of orbitals : ', Ndim
             Write(50,*) 'Symm. decomp  : ', Symm
             if (Projector) then
                Write(50,*) 'Projective version'
                Write(50,*) 'Theta         : ', Theta
                Write(50,*) 'Tau_max       : ', beta
                Write(50,*) '# of particles: ', N_part
             else
                Write(50,*) 'Finite temperture version'
                Write(50,*) 'Beta          : ', Beta
             endif
             Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
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

          
          Implicit none
          
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)

          If (Lattice_Type /=  "Square")  then
             Write(6,*) 'The plain vanilla Hubbard model is only defined for the square lattice'
             stop
          Endif
          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
          Ndim = Latt%N
          
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

          Integer :: nf , I, Ix, Iy
          allocate(Op_T(1,N_FL))
          do nf = 1,N_FL
             Call Op_make(Op_T(1,nf),Ndim)
             Do I = 1,Latt%N
                Ix = Latt%nnlist(I,1,0)
                Op_T(1,nf)%O(I,  Ix) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                Op_T(1,nf)%O(Ix, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                If ( L2 > 1 ) then
                   Iy = Latt%nnlist(I,0,1)
                   Op_T(1,nf)%O(I,  Iy) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                   Op_T(1,nf)%O(Iy, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                endif
                Op_T(1,nf)%O(I,  I ) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                Op_T(1,nf)%P(i) = i 
             Enddo
             Op_T(1,nf)%g      = -Dtau
             Op_T(1,nf)%alpha  =  cmplx(0.d0,0.d0, kind(0.D0))
             Call Op_set(Op_T(1,nf))
          enddo
          

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
          
          Integer                              :: nf, Ix, Iy, I, n
          Real (Kind=Kind(0.d0)), allocatable  :: H0(:,:),  U0(:,:), E0(:)
          Real (Kind=Kind(0.d0))               :: Pi = acos(-1.d0), Delta = 0.01d0
#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)

          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif
          
          Allocate(WF_L(N_FL),WF_R(N_FL))
          do nf=1,N_FL
             Call WF_alloc(WF_L(nf),Ndim,N_part)
             Call WF_alloc(WF_R(nf),Ndim,N_part)
          enddo

          
          Allocate(H0(Ndim,Ndim),  U0(Ndim, Ndim),  E0(Ndim) )
          H0 = 0.d0; U0 = 0.d0;  E0=0.d0
          Do I = 1,Latt%N
             Ix = Latt%nnlist(I,1,0)
             H0(I,  Ix) = -Ham_T*(1.d0   +   Delta*cos(Pi*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0))))
             H0(Ix, I ) = -Ham_T*(1.d0   +   Delta*cos(Pi*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0))))
             If (L2  > 1 ) Then
                Iy = Latt%nnlist(I,0,1)
                H0(I,  Iy) = -Ham_T *(1.d0  -   Delta)
                H0(Iy, I ) = -Ham_T *(1.d0  -   Delta)
             Endif
          Enddo
          Call  Diag(H0,U0,E0)
!!$          Do I = 1,Ndim
!!$             Write(6,*) I,E0(I)
!!$          Enddo
          Do nf = 1,N_FL
             do n=1,N_part
                do I=1,Ndim
                   WF_L(nf)%P(I,n)=U0(I,n)
                   WF_R(nf)%P(I,n)=U0(I,n)
                enddo
             enddo
             WF_L(nf)%Degen = E0(N_part+1) - E0(N_part)
             WF_R(nf)%Degen = E0(N_part+1) - E0(N_part)
          enddo
          
          
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

          Deallocate(H0,  U0,  E0 )

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
          
          Integer :: nf, I
          Real (Kind=Kind(0.d0)) :: X
          

          Allocate(Op_V(Ndim,N_FL))

          do nf = 1,N_FL
             do i  = 1, Ndim
                Call Op_make(Op_V(i,nf), 1)
             enddo
          enddo
          
          Do nf = 1,N_FL
             X = 1.d0
             if (nf == 2)  X = -1.d0
             Do i = 1,Ndim
                Op_V(i,nf)%P(1)   = I
                Op_V(i,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                Op_V(i,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0))) 
                Op_V(i,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                Op_V(i,nf)%type   = 2
                Call Op_set( Op_V(i,nf) )
             Enddo
          Enddo
             
          
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
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename


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
          Allocate ( Obs_eq(5) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Ns = Latt%N;  No = 1;  Filename ="Green"
             case (2)
                Ns = Latt%N;  No = 1;  Filename ="SpinZ"
             case (3)
                Ns = Latt%N;  No = 1;  Filename ="SpinXY"
             case (4)
                Ns = Latt%N;  No = 1;  Filename ="SpinT"
             case (5)
                Ns = Latt%N;  No = 1;  Filename ="Den"
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
                   Ns = Latt%N; No = 1;  Filename ="Green"
                case (2)
                   Ns = Latt%N; No = 1;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N; No = 1;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N; No = 1;  Filename ="SpinT"
                case (5)
                   Ns = Latt%N; No = 1;  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1-2*Thtrot
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
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
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY, ZDen
          Integer :: I,J, imj, nf,  Ix, Iy
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
          Zkin = Zkin* dble(N_SUN)
          Do I = 1,Latt%N
             Ix = Latt%nnlist(I,1,0)
             Zkin = Zkin  + GRC(I,Ix,1)  + GRC(Ix,I,1)  &
                  &       + GRC(I,Ix,2)  + GRC(Ix,I,2)
          Enddo
          If (L2 > 1) then
             Do I = 1,Latt%N
                Iy = Latt%nnlist(I,0,1)
                Zkin = Zkin + GRC(I,Iy,2)  + GRC(Iy,I,2)   &
                     &      + GRC(I,Iy,1)  + GRC(Iy,I,1)  
             Enddo
          Endif
          Zkin = Zkin*cmplx(-Ham_T,0.d0,Kind(0.d0)) 
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, 2)
          Enddo
          Zpot = Zpot*ham_U
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do I = 1,Ndim
             Zrho = Zrho + Grc(i,i,1) +  Grc(i,i,2)
          enddo
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

          Do I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N         =  Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign  =  Obs_eq(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          Do I = 1,Latt%N
             Do J = 1,Latt%N
                imj  = latt%imj(I,J)
                ZXY  = GRC(I,J,1) * GR(I,J,2) +  GRC(I,J,2) * GR(I,J,1) 
                ZZ   = GRC(I,J,1) * GR(I,J,1) +  GRC(I,J,2) * GR(I,J,2)    + &
                       (GRC(I,I,2) - GRC(I,I,1))*(GRC(J,J,2) - GRC(J,J,1))  

                ZDen = (GRC(I,I,1) + GRC(I,I,2)) * (GRC(I,I,1) + GRC(I,I,2)) + &
                     &  GRC(I,J,1) * GR(I,J,1)   +  GRC(I,J,2) * GR(I,J,2)  
                Obs_eq(1)%Obs_Latt(imj,1,1,1) =  Obs_eq(1)%Obs_Latt(imj,1,1,1) + (GRC(I,J,1) + GRC(I,J,2))*ZP*ZS
                Obs_eq(2)%Obs_Latt(imj,1,1,1) =  Obs_eq(2)%Obs_Latt(imj,1,1,1) +  ZZ  *ZP*ZS
                Obs_eq(3)%Obs_Latt(imj,1,1,1) =  Obs_eq(3)%Obs_Latt(imj,1,1,1) +  ZXY *ZP*ZS
                Obs_eq(4)%Obs_Latt(imj,1,1,1) =  Obs_eq(4)%Obs_Latt(imj,1,1,1) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                Obs_eq(5)%Obs_Latt(imj,1,1,1) =  Obs_eq(5)%Obs_Latt(imj,1,1,1) +  ZDen * ZP * ZS 

                
             enddo
             Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + (GRC(I,I,1) + GRC(I,I,2)) *  ZP*ZS
          enddo
          
          
          

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
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, ZDEN
          Real    (Kind=Kind(0.d0)) :: X
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          If (NT == 0 ) then
             Do I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N         =  Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign  =  Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             Enddo
          Endif
          Do I = 1,Latt%N
             Do J = 1,Latt%N
                imj  = latt%imj(I,J)

                ZZ   =      (GTT(I,I,1) -  GTT(I,I,2) ) * ( G00(J,J,1)  -  G00(J,J,2) )   &
                     &    -  G0T(J,I,1) * GT0(I,J,1)  -  G0T(J,I,2) * GT0(I,J,2) 
                ZXY  =    -  G0T(J,I,1) * GT0(I,J,2)  -  G0T(J,I,2) * GT0(I,J,1) 
                

                ZDen =   (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I,I,1) - GTT(I,I,2) ) * &
                     &   (cmplx(2.d0,0.d0,kind(0.d0)) -  G00(J,J,1) - G00(J,J,2) )   &
                     &   -G0T(J,I,1) * GT0(I,J,1)  -  G0T(J,I,2) * GT0(I,J,2) 

                Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) + (GT0(I,J,1) + GT0(I,J,2))*ZP*ZS
                Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) +  ZZ  *ZP*ZS
                Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) +  ZXY *ZP*ZS
                Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                Obs_tau(5)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(5)%Obs_Latt(imj,NT+1,1,1) +  ZDen * ZP * ZS 

                
             enddo
             Obs_tau(5)%Obs_Latt0(1) = Obs_tau(5)%Obs_Latt0(1) + &
                  &                   (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I,I,1) - GTT(I,I,2))  *  ZP*ZS
          enddo
          

          
        end Subroutine OBSERT

#include "Hamiltonian_Hubbard_include.h"        

      
    end Module Hamiltonian
