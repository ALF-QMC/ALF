!  Copyright (C) 2016-2020 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.
!  This is for the Kondo project with z-frustration.  Kondo on the Honeycomb lattice.
    Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Observables
      Use Matrix
      Use Fields_mod
      

      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Type (WaveFunction), dimension(:),   allocatable  :: WF_L
      Type (WaveFunction), dimension(:),   allocatable  :: WF_R
      Type  (Fields)       :: nsigma
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot, Thtrot
      Logical              :: Projector
!>    Defines MPI communicator 
      Integer              :: Group_Comm
      Logical              :: Symm =.false.
      
      ! What is below is  private 
      
      Type (Lattice),       private :: Latt
      Type (Unit_cell),     private :: Latt_Unit
      !Integer, parameter,   private :: Norb=4
      Integer, allocatable, private :: List(:,:), Invlist(:,:)
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T, Ham_U,  Ham_J, Ham_Jz, del_p(2)
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Integer,              private :: Checkerboard
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord 
      Real (Kind=Kind(0.d0)),        private :: Bound

!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau

      contains 

        Subroutine Ham_Set
#ifdef MPI
          Use mpi
#endif
          Implicit none


          integer :: ierr

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_Kondo/  ham_T, Ham_U,  Ham_J, Ham_Jz,  Dtau, Beta, Checkerboard


#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
          !          NAMELIST /VAR_Model/  N_FL,  N_SUN,  ham_T , ham_xi, ham_h, ham_J,  ham_U, Ham_Vint, &
          !               &         Dtau, Beta



#ifdef MPI
          If (Irank == 0 ) Then 
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
             CLOSE(5)
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          N_FL = 2
          N_SUN = 1

          
#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_Kondo)
             CLOSE(5)
#ifdef MPI
          endif
          CALL MPI_BCAST(ham_T       ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_U       ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_J       ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_Jz      ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau        ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta        ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Checkerboard,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

          Call Ham_hop
          Ltrot = nint(beta/dtau)
          Projector = .false.
          Theta = 0.d0
          Thtrot = 0
#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) 'Lattice size  : ', L1,L2
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'Checkerboard  : ', Checkerboard
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'U             : ', Ham_U
             Write(50,*) 'J             : ', Ham_J
             Write(50,*) 'Jz            : ', Ham_Jz
#if defined(STAB1) 
             Write(50,*) 'STAB1 flag is on'
#endif
#if defined(QRREF) 
             Write(50,*) 'QRREF flag is on'
#endif
             close(50)
#ifdef MPI
          endif
#endif
          call Ham_V
        end Subroutine Ham_Set
!=============================================================================

        Subroutine Ham_Latt
          Implicit none
          !Set the lattice
          Integer :: no, I, nc
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          If ( Lattice_type =="Square" ) then
             a1_p(1) =  1.d0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.d0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             N_coord   = 2
             If ( L1 == 1 .or. L2 == 1 ) then 
                N_coord   = 1
                If (L1 == 1 ) then 
                   Write(6,*) ' For one dimensional systems set  L2 = 1 ' 
                   Stop
                endif
             endif
             Latt_Unit%Norb      = 2
             Allocate (Latt_unit%Orb_pos_p(2,2))
             Latt_Unit%Orb_pos_p(1,:) = 0.d0
             Latt_Unit%Orb_pos_p(2,:) = 0.d0
             
          elseif ( Lattice_type=="Honeycomb" ) then
             a1_p(1) =  1.d0   ; a1_p(2) =  0.d0
             a2_p(1) =  0.5d0  ; a2_p(2) =  sqrt(3.d0)/2.d0
             del_p   =  (a2_p - 0.5d0*a1_p ) * 2.d0/3.d0
             
             L1_p    =  dble(L1) * a1_p
             L2_p    =  dble(L2) * a2_p
             
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             Latt_Unit%Norb      = 4
             Allocate (Latt_unit%Orb_pos_p(4,2))
             Latt_Unit%Orb_pos_p(1,:) = 0.d0
             Latt_Unit%Orb_pos_p(2,:) = del_p(:)
             Latt_Unit%Orb_pos_p(3,:) = 0.d0
             Latt_Unit%Orb_pos_p(4,:) = del_p(:)
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif

          Ndim = Latt%N*Latt_Unit%Norb
          Allocate (List(Ndim,Latt_Unit%Norb), Invlist(Latt%N,Latt_Unit%Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Latt_unit%Norb
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
                ! no = 1 c-electron_A
                ! no = 2 c-electron_B
                ! no = 3 f-electrons_A
                ! no = 4 f-electrons_B
             Enddo
          Enddo

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !Setup the hopping
          !Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}
          Integer :: I, I1, I2,I3,no, n, Ncheck, nc, ncoord
          Real (Kind=Kind(0.d0)) :: X

          !allocate(Exp_T   (Ndim,Ndim,N_FL) )
          !allocate(Exp_T_M1(Ndim,Ndim,N_FL) )
          IF (Checkerboard == 1 ) then 
             Ncheck = 3*Latt%N
             allocate(Op_T(Ncheck,N_FL))
             do n = 1,N_FL
                nc = 0
                Do ncoord = 1,3
                   Do I = 1,Latt%N
                      nc = nc + 1
                      if      ( ncoord == 1 ) then 
                         I1 = invlist(I,1)
                         I2 = invlist(I,2)
                      elseif  ( ncoord == 2 ) then
                         I1 = invlist(I,1) 
                         I2 = Invlist( Latt%nnlist(I,1,-1),2 )
                      elseif  ( ncoord == 3 ) then
                         I1 = invlist(I,1) 
                         I2 = invlist( Latt%nnlist(I,0,-1),2 )
                      endif
                      Call Op_make(Op_T(nc,n),2)
                      Op_T(nc,n)%P(1) = I1
                      Op_T(nc,n)%P(2) = I2
                      Op_T(nc,n)%O( 1 , 2 ) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                      Op_T(nc,n)%O( 2 , 1 ) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                      Op_T(nc,n)%g=cmplx(-Dtau,0.d0,Kind(0.d0))
                      !Write(6,*) 'In Ham_hop', Ham_T
                      Call Op_set(Op_T(nc,n)) 
                      !Write(6,*) 'In Ham_hop 1'
                      !Do I = 1,2*Latt%N
                      !   Write(6,*) Op_T(nc,n)%E(i)
                      !enddo
                      !Call Op_exp( cmplx(-Dtau,0.d0), Op_T(n), Exp_T   (:,:,n) )
                      !Call Op_exp( cmplx( Dtau,0.d0), Op_T(n), Exp_T_M1(:,:,n) )
                   enddo
                enddo
             enddo
          else
             Ncheck = 1
             allocate(Op_T(Ncheck,N_FL))
             do n = 1,N_FL
                Do nc = 1,NCheck
                   Call Op_make(Op_T(nc,n),2*Latt%N)
                   DO I = 1, Latt%N
                      Do no = 0,1
                         Op_T(nc,n)%P(I + no*Latt%N ) = invlist(I,no + 1)
                      enddo
                   enddo
                   Do I = 1,Latt%N
                      I1 = I
                      Op_T(nc,n)%O( I ,I1 + Latt%N) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                      Op_T(nc,n)%O( I1+Latt%N, I  ) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                      
                      I2   = Latt%nnlist(I,1,-1)
                      Op_T(nc,n)%O( I ,I2 + Latt%N) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                      Op_T(nc,n)%O( I2+Latt%N, I  ) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                      
                      I3   = Latt%nnlist(i,0,-1)
                      Op_T(nc,n)%O( I ,I3 + Latt%N) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                      Op_T(nc,n)%O( I3+Latt%N, I  ) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                   ENDDO
                   Op_T(nc,n)%g=cmplx(-Dtau,0.d0,Kind(0.d0))
                   !Write(6,*) 'In Ham_hop', Ham_T
                   Call Op_set(Op_T(nc,n)) 
                   !Write(6,*) 'In Ham_hop 1'
                   Do I = 1,2*Latt%N
                      Write(6,*) Op_T(nc,n)%E(i)
                   enddo
                   !Call Op_exp( cmplx(-Dtau,0.d0), Op_T(n), Exp_T   (:,:,n) )
                   !Call Op_exp( cmplx( Dtau,0.d0), Op_T(n), Exp_T_M1(:,:,n) )
                enddo
             enddo
          endif

        end Subroutine Ham_hop
        
!===================================================================================           
        
        Subroutine Ham_V
          
          Use Predefined_Int

          Implicit none 
          
          Integer :: nf, nth, n, n1, n2, n3, n4, I, I1, I2, J,  Ix, Iy, nc, no
          Real (Kind=Kind(0.d0)) :: X_p(2), X1_p(2), X2_p(2), X, XJ


          ! Number of opertors. 
          ! Latt%N * 2 for the Hubbard
          ! Latt%N * 2 For the Kondo  
          ! Latt%N * 6 for the Frustrating Ising   
          ! Total of Latt%N * 10
          
          Allocate( Op_V(10*Latt%N,N_FL) )

          do nf = 1,N_fl
             do i = 4*Latt%N +1, 10*Latt%N
                Call Op_make(Op_V(i,nf),2)  !  For the Kondo and  Ising. 
             enddo
          enddo
          
          nc = 0
          !  Hubbard
          Do i = 1,Latt%N
             Do no = 3,4
                nc = nc + 1
                Call Predefined_Int_U_SUN( OP_V(nc,1), Invlist(I,no), 2 , DTAU, Ham_U  )
                Call Predefined_Int_U_SUN( OP_V(nc,2), Invlist(I,no), 2 , DTAU, Ham_U  )
             Enddo
          Enddo
          !  Kondo
          DO i = 1,Latt%N
             Do no = 1,2
                nc = nc + 1
                Call Predefined_Int_V_SUN( OP_V(nc,1), Invlist(I,no  ), Invlist(I,no+2), 4, DTAU, Ham_J  ) 
                Call Predefined_Int_V_SUN( OP_V(nc,2), Invlist(I,no  ), Invlist(I,no+2), 4, DTAU, Ham_J  ) 
             Enddo
          Enddo
          ! Sz-Sz Coupling
          Do i = 1,Latt%N
             do no = 3,4
                do nth = 1,3
                   If ( nth == 1 ) J = latt%nnlist(I, 1,0)
                   If ( nth == 2 ) J = latt%nnlist(I, 0,1)
                   If ( nth == 3 ) J = latt%nnlist(I,-1,1)
                   nc = nc + 1
                   Call Predefined_Int_Jz (  Op_V(nc,1),Op_V(nc,2), Invlist(I,no), Invlist(J,no),  DTAU, Ham_Jz  ) 
                Enddo
             Enddo
          Enddo
          
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)
          Implicit none
          Integer, Intent(IN) :: n,nt
          Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new

          S0 = 1.d0
        end function S0
!===================================================================================
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
                Ns = Latt%N;  No = Norb;  Filename ="Den"
             case (5)
                Ns = Latt%N;  No = Norb;  Filename ="SpinT"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau == 1) then 
             ! Equal time correlators
             Allocate ( Obs_tau(4) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Ns = Latt%N; No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N; No = Norb;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N; No = Norb;  Filename ="Den"
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
             Call  Print_bin_Latt(Obs_eq(I),Latt,dtau,Group_Comm)
          enddo
          If (Ltau  == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call  Print_bin_Latt(Obs_tau(I),Latt,dtau,Group_Comm)
             enddo
          endif

        end Subroutine Pr_obs

!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, G(4,4,N_FL)
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZXY, ZZ
          Integer :: I,J, no,no1, n, n1, imj, nf, I1, I2,I3, J1, J2, I0, J0, ns, nc, NC_tot
          Integer :: no_I, no_J
          
          Real (Kind=Kind(0.d0)) ::  X
          
          ZP = PHASE/cmplx(Real(Phase,Kind=Kind(0.d0)),0.d0,Kind(0.d0))
          ZS = cmplx(Real(Phase,Kind=Kind(0.d0))/Abs(Real(Phase,Kind=Kind(0.d0))), 0.d0,Kind(0.d0))
          
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   ZK = cmplx(0.d0,0.d0,Kind(0.d0))
                   If ( I == J ) ZK = cmplx(1.d0,0.d0,Kind(0.d0))
                   GRC(I,J,nf)  = ZK - GR(J,I,nf)
                Enddo
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >
          ! Compute scalar observables.

          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          
          Zkin = cmplx(0.d0,0.d0,Kind(0.d0))
          Do nf = 1,N_FL
             Do I = 1,Latt%N
                I0 = invlist(I,1)
                I1 = invlist(I,2)
                I2 = Invlist( Latt%nnlist(I,1,-1),2 )
                I3 = invlist( Latt%nnlist(I,0,-1),2 )  
                Zkin = Zkin + Grc(I0,I1,nf) +  Grc(I1,I0,nf)  + &
                     &        Grc(I0,I2,nf) +  Grc(I2,I0,nf)  + &
                     &        Grc(I0,I3,nf) +  Grc(I3,I0,nf)
             Enddo
          Enddo
          Zkin = - Zkin*cmplx( Ham_T*dble(N_SUN), 0.d0,Kind(0.d0) )
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          Zrho = cmplx(0.d0,0.d0,Kind(0.d0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho*cmplx( dble(N_SUN), 0.d0,Kind(0.d0) )
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS

          
          ZPot = cmplx(0.d0,0.d0,Kind(0.d0))
          Do no = 3,4
             Do I = 1,Latt%N
                I1 = Invlist(I,no)
                ZPot = ZPot + Grc(I1,I1,1) * Grc(I1,I1,2)
             Enddo
          Enddo
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO

          Do I1 = 1,Ndim
             I    = List(I1,1)
             no_I = List(I1,2)
             Do J1 = 1,Ndim
                J    = List(J1,1)
                no_J = List(J1,2)
                imj = latt%imj(I,J)
                ! Green
                Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                     &                          ( GRC(I1,J1,1) + GRC(I1,J1,2)) *  ZP*ZS 
                ! SpinZ
                ZZ =       GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                     &    (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))  
                Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + ZZ* ZP*ZS
                ! SpinXY
                ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                ZXY =  GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1) 
                Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + ZXY* ZP*ZS
                ! SpinT
                Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                !Den
                Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + &
                     & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                     &   (GRC(I1,I1,2) + GRC(I1,I1,1))*(GRC(J1,J1,2) + GRC(J1,J1,1))     ) * ZP*ZS
             enddo
             Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  (GRC(I1,I1,1) + GRC(I1,I1,2)) * ZP*ZS
          enddo
          
        end Subroutine Obser
!==========================================================
        
        Subroutine OBSERT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J
          

          ZP = PHASE/cmplx(Real(Phase,Kind=Kind(0.d0)),0.d0,Kind(0.d0))
          ZS = cmplx(Real(Phase,Kind=Kind(0.d0))/Abs(Real(Phase,Kind=Kind(0.d0))), 0.d0,Kind(0.d0))
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          Do I1 = 1,Ndim
             I    = List(I1,1)
             no_I = List(I1,2)
             Do J1 = 1,Ndim
                J    = List(J1,1)
                no_J = List(J1,2)
                imj = latt%imj(I,J)
                !Green
                Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                     &   +   cmplx(0.5d0,0.d0,Kind(0.d0))*( GT0(I1,J1,1) + GT0(I1,J1,2) ) * ZP* ZS 
                
                !SpinZ
                Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                     & +  ( &
                     &    (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                     &  - (G0T(J1,I1,1) * GT0(I1,J1,1)  +  G0T(J1,I1,2) * GT0(I1,J1,2) )    )*ZP*ZS
                
                !SpinXY
                Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                     &  - &
                     &   (G0T(J1,I1,1) * GT0(I1,J1,2)  +  G0T(J1,I1,2) * GT0(I1,J1,1))*ZP*ZS
                !Den
                Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                     & +  (                                        &  
                     &    (cmplx(2.D0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2) ) * &
                     &    (cmplx(2.D0,0.d0,kind(0.d0)) - G00(J1,J1,1) - G00(J1,J1,2) )   &
                     & -  ( G0T(J1,I1,1) * GT0(I1,J1,1) + G0T(J1,I1,2) * GT0(I1,J1,2) )  )*ZP*ZS     
                
             enddo
             
             Obs_tau(4)%Obs_Latt0(no_I) =  Obs_tau(4)%Obs_Latt0(no_I) + &
                  &       (cmplx(2.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2)) * ZP * ZS
          Enddo
          
        end Subroutine OBSERT
!========================================================================
        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
          
          !>  The input is the field nsigma declared in this module. This routine generates a 
          !>  global update with  and returns the propability  
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)  
          !>   
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          type (Fields),  Intent(IN)  :: nsigma_old

        End Subroutine Global_move
!========================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          type (Fields),  Intent(IN)  :: nsigma_old

        end Function Delta_S0_global
!---------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!> On input: 
!> GR(tau,m) as defined in  Global_tau_mod_PlaceGR and the direction of updating scheme
!> direction=u --> You are visiting the time slices from tau = 1  to tau =Ltrot
!> direction=d --> You are visiting the time slices from tau = Ltrot to tau = 1
!> 
!> On input the field configuration is in the array nsigma.
!> On output: 
!> Flip_list   ::  A list of spins that are to be fliped. Refers to the entires  in OP_V
!> Flip_values ::  The values of the fliped spins
!> Flip_length ::  The number of flips. The first Flip_length entries of Flip_list and Flip_values are relevant
!> S0_ratio          = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
!> T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)  
!> T0_proposal       = T0 ( sigma -> sigma_new )
!--------------------------------------------------------------------
          
          Implicit none 
          Real (Kind= kind(0.d0)), INTENT(INOUT) :: T0_Proposal_ratio,  S0_ratio
          Integer,    allocatable, INTENT(INOUT) :: Flip_list(:)
          Real (Kind= Kind(0.d0)), INTENT(INOUT) :: Flip_value(:)
          Integer, INTENT(INOUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau

        end Subroutine Global_move_tau

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
      end Subroutine Overide_global_tau_sampling_parameters


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

        ! The user can set the initial configuration
        
        Implicit none

        Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(OUT) :: Initial_field

      end Subroutine Hamiltonian_set_nsigma
      
    end Module Hamiltonian
