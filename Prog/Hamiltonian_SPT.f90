    !Model Hamiltonian for interaction-induced topological reduction
    Module Hamiltonian

      Use Operator_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      

      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot
!>    Defines MPI communicator 
      Integer              :: Group_Comm


      
      ! What is below is  private 
      
      Type (Lattice),       private :: Latt
      Integer, parameter,   private :: Norb=16
      Integer, allocatable, private :: List(:,:), Invlist(:,:)
      Integer,              private :: L1, L2, FlagSym
      real (Kind=Kind(0.d0)),        private :: Ham_T, Ham_Vint,  Ham_Lam
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta
      Character (len=64),   private :: Model, Lattice_type
      Complex (Kind=Kind(0.d0)),     private :: Gamma_M(4,4,5), Sigma_M(2,2,0:3)
      Complex (Kind=Kind(0.d0)),     private :: Gamma_13(4,4), Gamma_23(4,4), Gamma_45(4,4)


      ! Observables
      Integer,                       private :: Nobs
      Complex (Kind=Kind(0.d0)), allocatable, private :: obs_scal(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Den_eq(:,:,:), Den_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  R_eq(:,:,:), R_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1_eq(:,:,:), U1_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  L_eq(:,:,:), L_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1xy_eq(:,:,:), U1xy_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Spinz_eq(:,:,:), Spinz_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Spinxy_eq(:,:,:), Spinxy_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  TRS_eq(:,:,:), TRS_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  PHS_eq(:,:,:), PHS_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  RS_eq(:,:,:), RS_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  C4S_eq(:,:,:), C4S_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  PxS_eq(:,:,:), PxS_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  SxS_eq(:,:,:), SxS_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  SzS_eq(:,:,:), SzS_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U11S_eq(:,:,:), U11S_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U12S_eq(:,:,:), U12S_eq0(:)

      ! For time displaced
      Integer,                       private :: NobsT
      Complex (Kind=Kind(0.d0)),              private :: Phase_tau
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green_tau(:,:,:,:), Den_tau(:,:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: U1_tau(:,:,:,:), U1xy_tau(:,:,:,:), U1xyG_tau(:,:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Spinz_tau(:,:,:,:), Spinxy_tau(:,:,:,:)
      
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Den_sus(:,:,:), Den_sus0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1_sus(:,:,:), U1_sus0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1xy_sus(:,:,:), U1xy_sus0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1xyG_sus(:,:,:), U1xyG_sus0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Spinz_sus(:,:,:), Spinz_sus0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Spinxy_sus(:,:,:), Spinxy_sus0(:)

      contains 

        Subroutine Ham_Set
#ifdef MPI
          Use mpi
#endif
          Implicit none

          integer :: ierr

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_SPT/  ham_T, Ham_Vint,  Ham_Lam,  Dtau, Beta, FlagSym


#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
             CLOSE(5)
#ifdef MPI
          endif

          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt


          N_FL  = 1
          N_SUN = 1
          
#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_SPT)
             CLOSE(5)
#ifdef MPI
          endif

          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_Vint ,1,MPI_REAL8,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_Lam  ,1,MPI_REAL8,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(FlagSym  ,1,MPI_INTEGER, 0,MPI_COMM_WORLD,ierr)
#endif

          Call Ham_hop
          Ltrot = nint(beta/dtau)
#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'V             : ', Ham_Vint
             Write(50,*) 'Lambda        : ', Ham_Lam
             Write(50,*) 'FlagSym       : ', FlagSym
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
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             !Write(6,*)  'Lattice: ', Ndim
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif

          Ndim = Latt%N*Norb
          Allocate (List(Ndim,Norb), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
                ! no = 1..4   Xi_1
                ! no = 5..8   Xi_2
                ! no = 9..12  Xi_3
                ! no = 13..16 Xi_4
             Enddo
          Enddo

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !  Setup the hopping
          !  Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, I2 ,no, no1, n, Ncheck, nc
          Integer, allocatable :: Invlist_1(:,:)
          Complex (Kind=Kind(0.d0)) :: Z


          ! Setup Gamma matrices
          Gamma_M = cmplx(0.d0, 0.d0, kind(0.D0))
          Sigma_M = cmplx(0.d0, 0.d0, kind(0.D0))
          Sigma_M(1,1,0) = cmplx( 1.d0, 0.d0, kind(0.D0))
          Sigma_M(2,2,0) = cmplx( 1.d0, 0.d0, kind(0.D0))
          Sigma_M(1,2,1) = cmplx( 1.d0, 0.d0, kind(0.D0))
          Sigma_M(2,1,1) = cmplx( 1.d0, 0.d0, kind(0.D0))
          Sigma_M(1,2,2) = cmplx( 0.d0,-1.d0, kind(0.D0))
          Sigma_M(2,1,2) = cmplx( 0.d0, 1.d0, kind(0.D0))
          Sigma_M(1,1,3) = cmplx( 1.d0, 0.d0, kind(0.D0))
          Sigma_M(2,2,3) = cmplx(-1.d0, 0.d0, kind(0.D0))
          Do no = 1,2
             Do no1 = 1,2
                Gamma_M(no+2,no1  ,1)  =  Sigma_M(no,no1,0) 
                Gamma_M(no  ,no1+2,1)  =  Sigma_M(no,no1,0) 
                Gamma_M(no+2,no1  ,2)  =  cmplx( 0.d0, 1.d0, kind(0.D0))*Sigma_M(no,no1,0)
                Gamma_M(no  ,no1+2,2)  =  cmplx( 0.d0,-1.d0, kind(0.D0))*Sigma_M(no,no1,0)
                Gamma_M(no  ,no1  ,3)  =  cmplx( 1.d0, 0.d0, kind(0.D0))*Sigma_M(no,no1,1)
                Gamma_M(no+2,no1+2,3)  =  cmplx(-1.d0, 0.d0, kind(0.D0))*Sigma_M(no,no1,1)
                Gamma_M(no  ,no1  ,4)  =  cmplx( 1.d0, 0.d0, kind(0.D0))*Sigma_M(no,no1,2)
                Gamma_M(no+2,no1+2,4)  =  cmplx(-1.d0, 0.d0, kind(0.D0))*Sigma_M(no,no1,2)
                Gamma_M(no  ,no1  ,5)  =  cmplx( 1.d0, 0.d0, kind(0.D0))*Sigma_M(no,no1,3)
                Gamma_M(no+2,no1+2,5)  =  cmplx(-1.d0, 0.d0, kind(0.D0))*Sigma_M(no,no1,3)
             Enddo
          Enddo
          
          Call mmult (Gamma_13,Gamma_M(:,:,1), Gamma_M(:,:,3) )
          Call mmult (Gamma_23,Gamma_M(:,:,2), Gamma_M(:,:,3) )
          Call mmult (Gamma_45,Gamma_M(:,:,4), Gamma_M(:,:,5) )
          

          Ncheck = 4
          Allocate ( Invlist_1(Latt%N,4) )
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,NCheck
                Call Op_make(Op_T(nc,n),Ndim/4)
                I1 = 0
                Do no = 1,4
                   DO I = 1, Latt%N
                      I1 = I1 + 1
                      Invlist_1(I,no) = I1
                      Op_T(nc,n)%P(I1) = Invlist(I, no + 4*(nc -1) )
                   enddo
                enddo
                Do I = 1,Latt%N
                   do no = 1,4
                      do no1 = 1,4
                         Z =  cmplx(2.d0 + Ham_Lam, 0.d0, kind(0.D0))*Gamma_M(no,no1,3)
                         Op_T(nc,n)%O( Invlist_1(I,no) ,Invlist_1(I,no1) )  =  Z
                      enddo
                   enddo
                   I1 =  Latt%nnlist(I,1,0)
                   do no = 1,4
                      do no1 = 1,4
                         Z = (cmplx(0.d0,1.d0, kind(0.D0))*Gamma_M(no,no1,1) &
                    &    + Gamma_M(no,no1,3))*cmplx(0.5d0,0.d0, kind(0.D0))
                         Op_T(nc,n)%O( invlist_1(I ,no  ), invlist_1(I1,no1  ) )  = &
                    & Op_T(nc,n)%O( invlist_1(I ,no  ), invlist_1(I1,no1  ) ) +  Z
                         Op_T(nc,n)%O( invlist_1(I1,no1 ), invlist_1(I ,no   ) )  = &
                    & Op_T(nc,n)%O( invlist_1(I1,no1 ), invlist_1(I ,no   ) ) + conjg(Z)
                      enddo
                   enddo
                   I2   = Latt%nnlist(I,0,1)
                   do no = 1,4
                      do no1 = 1,4
                         Z =   ( cmplx(0.d0,1.d0, kind(0.D0)) * Gamma_M(no,no1,2) &
                    &    + Gamma_M(no,no1,3) )*cmplx(0.5d0,0.d0, kind(0.D0))
                         Op_T(nc,n)%O( invlist_1(I ,no ), invlist_1(I2,no1  ) )  = &
                    & Op_T(nc,n)%O( invlist_1(I ,no ), invlist_1(I2,no1  ) ) + Z
                         Op_T(nc,n)%O( invlist_1(I2,no1), invlist_1(I ,no   ) )  = &
                    & Op_T(nc,n)%O( invlist_1(I2,no1), invlist_1(I ,no   ) ) + conjg(Z)
                      enddo
                   enddo
                enddo
                Op_T(nc,n)%g=cmplx(-Dtau*Ham_T,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n)) 
                ! Just for tests
                Do I = 1, Ndim/4
                   Write(6,*) i,Op_T(nc,n)%E(i)
                enddo
             enddo
          enddo

          deallocate (Invlist_1) 

        end Subroutine Ham_hop
!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, nth, n, n1, n2, n3, n4, I, I1, I2, J,  Ix, Iy, nc, no,no1, ns, npm 
          Integer :: nxy, noff
          Real (Kind=Kind(0.d0)) :: X_p(2), X1_p(2), X2_p(2), X, XJ, Xpm

          Complex (Kind=Kind(0.d0)) :: Ps(4,4,2), Ps_G5(4,4,2), Tmp(4,4), Z
          Complex (Kind=Kind(0.d0)) :: Sx(16,16,2,2), Sy(16,16,2,2)


          Ps = cmplx(0.d0, 0.d0, kind(0.D0))
          Call mmult (Tmp,Gamma_M(:,:,3), Gamma_M(:,:,4) )
          
! !           Do n=1,5
!          write(*,*) "Gamma34"
!          do no=1,4
!          do no1=1,4
!            write(*,*) Tmp(no,no1)
!          enddo
!          write(*,*)
!          enddo
!          write(*,*)
! !           enddo
          
          do ns = 1,2
             if (ns == 1) X =  1.d0/2d0
             if (ns == 2) X = -1.d0/2.d0
             Do I = 1,4
                Do J = 1,4
                   Z = cmplx(0.d0, 0.d0, kind(0.D0))
                   if ( I == J )  Z = cmplx(1.d0/2.d0, 0.d0, kind(0.D0))
                   Ps(I,J,ns) =   Z  + cmplx(0.d0, X, kind(0.D0)) * tmp(I,J)
                Enddo
             Enddo
          Enddo
          
          Do ns = 1,2
             Call mmult ( Ps_G5(:,:,ns), Ps(:,:,ns), Gamma_M(:,:,5) )
          enddo
          
!           Do ns=1,2
!          write(*,*) "PsG5 ", ns
!          do no=1,4
!          do no1=1,4
!            write(*,*) Ps_G5(no,no1,ns)
!          enddo
!          write(*,*)
!          enddo
!          write(*,*)
!           enddo
      
          Sx = cmplx(0.d0,0.d0, kind(0.D0))
          Sy = cmplx(0.d0,0.d0, kind(0.D0))
          Do ns = 1,2
             Do npm = 1,2
                if (npm == 1) Xpm =  1.0
                if (npm == 2) Xpm = -1.0
                Do no = 1,4
                   do no1 = 1,4
                      Sx(no    , no1 + 4 ,ns,npm) =  cmplx(1.d0, 0.d0, kind(0.D0))*Ps_G5(no,no1,ns)
                      Sx(no +4 , no1     ,ns,npm) =  cmplx(1.d0, 0.d0, kind(0.D0))*Ps_G5(no,no1,ns)
                      Sx(no +8 , no1 + 12,ns,npm) =  cmplx(xpm,  0.d0, kind(0.D0))*Ps_G5(no,no1,ns)
                      Sx(no+12 , no1 + 8 ,ns,npm) =  cmplx(xpm,  0.d0, kind(0.D0))*Ps_G5(no,no1,ns)
                      
                      Sy(no    , no1 + 4 ,ns,npm) =  cmplx(0.d0, -1.d0    , kind(0.D0))*Ps_G5(no,no1,ns)
                      Sy(no +4 , no1     ,ns,npm) =  cmplx(0.d0,  1.d0    , kind(0.D0))*Ps_G5(no,no1,ns)
                      Sy(no +8 , no1 + 12,ns,npm) =  cmplx(0.d0, -1.d0*xpm, kind(0.D0))*Ps_G5(no,no1,ns)
                      Sy(no+12 , no1 + 8 ,ns,npm) =  cmplx(0.d0,  1.d0*xpm, kind(0.D0))*Ps_G5(no,no1,ns)
                   enddo
                enddo
             enddo
          enddo


          ! Number of opertors 8 per unit cell
          Allocate( Op_V(8*Latt%N,N_FL) )
          do nf = 1,N_FL
             do i  = 1, 8*Latt%N
                Call Op_make(Op_V(i,nf),Norb/2) 
             enddo
          enddo
          nc = 0
          Do nf = 1,N_FL
             do nxy = 1,2
                do ns = 1,2
                   noff=2
                   if (ns==2) noff=1
                   do npm = 1,2 
                      Xpm = 1.d0
                      if (npm == 2) Xpm = -1.d0
                      Do i = 1,Latt%N
                         nc = nc + 1 
                         Do no = 1,Norb/2
                            Op_V(nc,nf)%P(no)   = Invlist(I,2*(no-1)+noff)  
                         enddo
                         Do no = 1,Norb/2
                            Do no1 = 1,Norb/2
                               If (nxy == 1)  Op_V(nc,nf)%O(no,no1) = Sx(2*(no-1)+noff,2*(no1-1)+noff,ns,npm)
                               If (nxy == 2)  Op_V(nc,nf)%O(no,no1) = Sy(2*(no-1)+noff,2*(no1-1)+noff,ns,npm)
                            Enddo
                         Enddo
                         Op_V(nc,nf)%g = SQRT(CMPLX(-Xpm*DTAU*Ham_Vint/8.d0, 0.D0, kind(0.D0))) 
                         Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                         Op_V(nc,nf)%type   = 2
                         Call Op_set( Op_V(nc,nf) )
                         ! The operator reads: 
                         ! g*s*( c^{dagger} O c   - alpha ))
                         ! with s the HS field.
                      Enddo
                   Enddo
                Enddo
             Enddo
          Enddo
          
        end Subroutine Ham_V
        
!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt 

          S0 = 1.d0
        end function S0

!===================================================================================           
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer :: I
          Allocate ( Obs_scal(8) )
          Allocate ( Den_eq(Latt%N,1,1), Den_eq0(1) ) 
          Allocate ( U1_eq(Latt%N,1,1), U1_eq0(1) )
          Allocate ( Spinz_eq(Latt%N,1,1), spinz_eq0(1) )
          Allocate ( L_eq(Latt%N,1,1), L_eq0(1) )
          
          if (FlagSym ==1) then
           Allocate ( R_eq(Latt%N,1,1), R_eq0(1) ) 
           Allocate ( U1xy_eq(Latt%N,1,1), U1xy_eq0(1) )
           Allocate ( Spinxy_eq(Latt%N,1,1), spinxy_eq0(1) ) 
           Allocate ( TRS_eq(Latt%N,1,1), TRS_eq0(1) )
           Allocate ( PHS_eq(Latt%N,1,1), PHS_eq0(1) )
           Allocate ( RS_eq(Latt%N,1,1), RS_eq0(1) )
           Allocate ( C4S_eq(Latt%N,1,1), C4S_eq0(1) )
           Allocate ( PxS_eq(Latt%N,1,1), PxS_eq0(1) )
           Allocate ( SxS_eq(Latt%N,1,1), SxS_eq0(1) )
           Allocate ( SzS_eq(Latt%N,1,1), SzS_eq0(1) )
           Allocate ( U11S_eq(Latt%N,1,1), U11S_eq0(1) )
           Allocate ( U12S_eq(Latt%N,1,1), U12S_eq0(1) )
       endif
          
          If (Ltau == 1) then 
             Allocate ( Green_tau(Latt%N,Ltrot+1,Norb,Norb), Den_tau(Latt%N,Ltrot+1,1,1) )
             Allocate ( U1_tau(Latt%N,Ltrot+1,1,1), U1xy_tau(Latt%N,Ltrot+1,1,1), U1xyG_tau(Latt%N,Ltrot+1,1,1) )
             Allocate ( Spinz_tau(Latt%N,Ltrot+1,1,1), Spinxy_tau(Latt%N,Ltrot+1,1,1) )
          Allocate ( Den_sus(Latt%N,1,1), Den_sus0(1) ) 
          Allocate ( U1_sus(Latt%N,1,1), U1_sus0(1) )
          Allocate ( U1xy_sus(Latt%N,1,1), U1xy_sus0(1) )
           Allocate ( U1xyG_sus(Latt%N,1,1), U1xyG_sus0(1) )
          Allocate ( Spinz_sus(Latt%N,1,1), spinz_sus0(1) )
          Allocate ( Spinxy_sus(Latt%N,1,1), spinxy_sus0(1) ) 
          endif
          
        end Subroutine Alloc_obs

!===================================================================================           
        
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          Nobs = 0
          Obs_scal  = cmplx(0.d0,0.d0, kind(0.D0))
          
          Den_eq    = cmplx(0.d0,0.d0, kind(0.D0))
          Den_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
          Spinz_eq    = cmplx(0.d0,0.d0, kind(0.D0))
          Spinz_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
          U1_eq    = cmplx(0.d0,0.d0, kind(0.D0))
          U1_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
          L_eq    = cmplx(0.d0,0.d0, kind(0.D0))
          L_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
          
          if (FlagSym ==1 ) then
           R_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           R_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           U1xy_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           U1xy_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           Spinxy_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           Spinxy_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           
           TRS_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           TRS_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           PHS_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           PHS_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           RS_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           RS_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           C4S_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           C4S_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           PxS_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           PxS_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           SxS_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           SxS_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           SzS_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           SzS_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           U11S_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           U11S_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
           U12S_eq    = cmplx(0.d0,0.d0, kind(0.D0))
           U12S_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
       endif

          If (Ltau == 1) then
             NobsT = 0
             Phase_tau = cmplx(0.d0,0.d0, kind(0.D0))
             Green_tau = cmplx(0.d0,0.d0, kind(0.D0))
             Den_tau = cmplx(0.d0,0.d0, kind(0.D0))
             U1_tau = cmplx(0.d0,0.d0, kind(0.D0))
             U1xy_tau = cmplx(0.d0,0.d0, kind(0.D0))
             U1xyG_tau = cmplx(0.d0,0.d0, kind(0.D0))
             Spinz_tau = cmplx(0.d0,0.d0, kind(0.D0))
             Spinxy_tau = cmplx(0.d0,0.d0, kind(0.D0))
          
          Den_eq    = cmplx(0.d0,0.d0, kind(0.D0))
          Den_eq0   = cmplx(0.d0,0.d0, kind(0.D0))
          U1_sus    = cmplx(0.d0,0.d0, kind(0.D0))
          U1_sus0   = cmplx(0.d0,0.d0, kind(0.D0))
          U1xy_sus    = cmplx(0.d0,0.d0, kind(0.D0))
          U1xy_sus0   = cmplx(0.d0,0.d0, kind(0.D0))
           U1xyG_sus    = cmplx(0.d0,0.d0, kind(0.D0))
           U1xyG_sus0   = cmplx(0.d0,0.d0, kind(0.D0))
          Spinz_sus    = cmplx(0.d0,0.d0, kind(0.D0))
          Spinz_sus0   = cmplx(0.d0,0.d0, kind(0.D0))
          Spinxy_sus    = cmplx(0.d0,0.d0, kind(0.D0))
          Spinxy_sus0   = cmplx(0.d0,0.d0, kind(0.D0))
          endif

        end Subroutine Init_obs
        
!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, ZL, Z, ZP,ZS, weight, tmp
          Integer :: I,J, no,no1, n, n1, imj, nf, I1, I2, J1, J2, Nc, Ix, Iy, Jx, Jy, Imx, Imy, Jmx, Jmy
          Integer :: a, b, c, d, signum, K, K1, L ,L1, nf1
          
          Real (Kind=Kind(0.d0)) :: G(4,4), X, FI, FJ
          
          Nobs = Nobs + 1
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          

          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   ZK = cmplx(0.d0, 0.d0, kind(0.D0))
                   If ( I == J ) ZK = cmplx(1.d0, 0.d0, kind(0.D0))
                   GRC(I,J,nf)  = ZK - GR(J,I,nf)
                Enddo
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >
          ! Compute scalar observables. 

          Zkin = cmplx(0.d0,0.d0, kind(0.D0))

          Nc = Size( Op_T,1)
          Do nf = 1,N_FL
             Do n = 1,Nc
                Do J = 1,Op_T(n,nf)%N
                   J1 = Op_T(n,nf)%P(J)
                   DO I = 1,Op_T(n,nf)%N
                      I1 = Op_T(n,nf)%P(I)
                      Zkin  = Zkin  + Op_T(n,nf)%O(i,j)*Grc(i1,j1,nf) 
                   Enddo
                ENddo
             Enddo
          Enddo
          Zkin = Zkin*cmplx( dble(N_SUN), 0.d0 , kind(0.D0))
          
          ZL = cmplx(0.d0,0.d0, kind(0.D0))
          
!           Nc = Size( Op_T,1)
          Do nf = 1,N_FL
!              Do n = 1,Nc
                Do I = 1,Latt%N
                   DO no = 1,4
                   DO a = 1,4
                   DO b = 1,4
                      I1 = InvList(I,4*(no-1)+a)
                J1 = InvList(I,4*(no-1)+b)
                      ZL  = ZL  + Gamma_M(a,b,3)*Grc(i1,j1,nf) 
                   Enddo
                   Enddo
                   Enddo
                ENddo
!              Enddo
          Enddo
          ZL = ZL*cmplx( dble(N_SUN), 0.d0 , kind(0.D0))

          Zrho = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho*cmplx( dble(N_SUN), 0.d0 , kind(0.D0))
          ZPot = cmplx(0.d0,0.d0, kind(0.D0))

          Nc = Size( Op_V,1)
          Do nf = 1,N_FL
             Do n = 1,Nc
          weight=(-1)**((n-1)/Latt%N)/8.d0
                Do J = 1,Op_V(n,nf)%N
                   J1 = Op_V(n,nf)%P(J)
                   DO I = 1,Op_V(n,nf)%N
                      if (abs(Op_V(n,nf)%O(i,j)) >= 0.00001) then
                 I1 = Op_V(n,nf)%P(I)
                 Do K = 1,Op_V(n,nf)%N
                   K1 = Op_V(n,nf)%P(K)
                   DO L = 1,Op_V(n,nf)%N
                     if (abs(Op_V(n,nf)%O(k,l)) >= 0.00001) then
                    L1 = Op_V(n,nf)%P(L)
                    tmp =  (   GRC(I1,L1,1) * GR (J1,K1,1)      +  &
                          &     GRC(I1,J1,1) * GRC(K1,L1,1)         )
                    ZPot  = ZPot  + weight*Op_V(n,nf)%O(i,j)*Op_V(n,nf)%O(k,l)*tmp
                     endif
                   Enddo
                 ENddo
                endif
                   Enddo
                ENddo
!           write(*,*) Zpot
             Enddo
          Enddo
          ZPot = ZPot*cmplx( dble(N_SUN), 0.d0 , kind(0.D0))

          Obs_scal(1) = Obs_scal(1) + zrho * ZP*ZS
          Obs_scal(2) = Obs_scal(2) + zkin*Ham_T * ZP*ZS
          Obs_scal(3) = Obs_scal(3) + Zpot*Ham_Vint * ZP*ZS
          Obs_scal(4) = Obs_scal(4) + (zkin*Ham_T +  Zpot*Ham_Vint)*ZP*ZS
          Obs_scal(5) = Obs_scal(5) + (zkin -  Zpot)*ZP*ZS
          Obs_scal(6) = Obs_scal(6) + (Zpot)*ZP*ZS
          Obs_scal(7) = Obs_scal(7) + (ZL)*ZP*ZS
          Obs_scal(8) = Obs_scal(8) + ZS
          ! You will have to allocate more space if you want to include more  scalar observables.
          
          
          DO I1 = 1,Ndim
             I  = List(I1,1)
             no = List(I1,2)
             DO J1 = 1, Ndim
                J = List(J1,1)
                no1 = list(J1,2)
                imj = latt%imj(I,J)
          
          tmp =  (   GRC(I1,J1,1) * GR (I1,J1,1)      +  &
                &     GRC(I1,I1,1) * GRC(J1,J1,1)         ) *ZP*ZS

                DEN_Eq (imj,1,1) = DEN_Eq (imj,1,1)   +  tmp
                     
          weight=cmplx(1.d0,0.d0, kind(0.D0))
          if ( (no>=9 .and. no1<=8) .or. (no<=8 .and. no1>=9) ) weight=-weight
          Spinz_eq (imj,1,1) = Spinz_eq (imj,1,1)   +   weight * 0.25 * tmp
          
          signum = 1
          if (((no-1)/4+1==2) .or. ((no-1)/4+1==4)) signum=-1
          if (((no1-1)/4+1==2) .or. ((no1-1)/4+1 ==4)) signum=-signum
          weight = cmplx(dble(signum),0.d0, kind(0.D0))
          U1_eq (imj,1,1) = U1_eq (imj,1,1)   +  weight*tmp*0.25

             enddo
             Den_eq0(1) = Den_eq0(1) +   GRC(I1,I1,1)*ZP*ZS 
          enddo
         
!        do I=1,Latt%N
!          do J=1,Latt%N
!            imj = latt%imj(I,J)
!            do a=1,4
!            do b=1,4
!           if ( abs(gamma_M(a,b,3))>0.01 ) then
!             do c=1,4
!             do d=1,4
!               weight = gamma_M(a,b,3)*Gamma_M(c,d,3)
!               if ( abs(weight) > 0.01 ) then
!                 do no=1,4
!                do no1=1,4
!                  I1 = Invlist(I,4*(no-1)+a)
!                  I2 = Invlist(I,4*(no-1)+b)
!                  J1 = Invlist(J,4*(no1-1)+c)
!                  J2 = Invlist(J,4*(no1-1)+d)
!                  tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
!                     &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!                  L_eq (imj,1,1) = L_eq (imj,1,1)   +  weight*tmp 
!                enddo
!                 enddo
!               endif
!             enddo
!             enddo
!           endif
!            enddo
!            enddo
!          enddo
! !        enddo
!        
!        do I=1,Latt%N
!          do a=1,4
!          do b=1,4
!            if ( abs(gamma_M(a,b,3))>0.01 ) then
!           do no=1,4
!             I1 = InvList(I,4*(no-1)+a)
!             J1 = InvList(I,4*(no-1)+b)
!             L_eq0(1)  = L_eq0(1)  + Gamma_M(a,b,3)*Grc(i1,j1,1)
!           enddo
!            endif
!          enddo
!          enddo
!        enddo
          
          if (FlagSym ==1) then
         do I=1,Latt%N
           do no=1,8
            do J=1,Latt%N
              imj = latt%imj(I,J)
              do no1=1,8
               I1 = Invlist(I,no)
               I2 = Invlist(I,no+8)
               J1 = Invlist(J,no1+8)
               J2 = Invlist(J,no1)
               
               tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                   &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
               Spinxy_eq (imj,1,1) = Spinxy_eq (imj,1,1)   +  tmp
               
               if (no<=4) then
                 I1 = Invlist(I,no)
                 I2 = Invlist(I,no+4)
               else
                 I1= Invlist(I,no+4)
                 I2 = Invlist(I,no+8)
               endif
               if (no1<=4) then
                 J1 = Invlist(J,no1+4)
                 J2 = Invlist(J,no1)
               else
                 J1 = Invlist(J,no1+8)
                 J2 = Invlist(J,no1+4)
               endif
               
               tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                   &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
               U1xy_eq (imj,1,1) = U1xy_eq (imj,1,1)   +  tmp
              enddo
            enddo
            
           enddo
         enddo
         
         do I=1,Latt%N
           Ix = Latt%nnlist(I,1,0)
           Iy = Latt%nnlist(I,0,1)
           Imx = Latt%nnlist(I,-1,0)
           Imy = Latt%nnlist(I,0,-1)
           do J=1,Latt%N
            Jx = Latt%nnlist(J,1,0)
            Jy = Latt%nnlist(J,0,1)
            Jmx = Latt%nnlist(J,-1,0)
            Jmy = Latt%nnlist(J,0,-1)
            imj = latt%imj(I,J)
            do no=1,4
              do a=1,4
              do b=1,4
                do no1=1,4
               do c=1,4
               do d=1,4
                 
      !               R correlation
                 weight = -gamma_45(a,b)*Gamma_45(c,d)
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(I,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(J,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   R_eq (imj,1,1) = R_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               TR symmetry check
                 weight = gamma_13(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = gamma_13(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_13(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_13(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = gamma_23(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = gamma_23(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_23(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_23(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_13(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_13(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = gamma_13(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = gamma_13(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_23(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = -gamma_23(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = gamma_23(a,b)*Gamma_13(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
                 weight = gamma_23(a,b)*Gamma_23(c,d)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               R symmetry check X°\dag gamma_4 X as correlation
                 weight = gamma_M(a,b,4)*Gamma_M(c,d,4)
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(I,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(J,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   RS_eq (imj,1,1) = RS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               C4 symmetry check
      !               11
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               12
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               13
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               14
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               21
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               22
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               23
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               24
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               31
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               32
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               33
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               34
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               41
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               42
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               43
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               44
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               Px symmetry check
      !               11
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               12
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               13
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               14
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Iy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               21
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               22
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               23
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               24
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Ix,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               31
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               32
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               33
                 weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               34
                 weight = gamma_M(a,b,1)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imy,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               41
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               42
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               43
                 weight = gamma_M(a,b,2)*Gamma_M(c,d,1)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmy,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               44
                 weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)/cmplx(4.d0,0.d0, kind(0.D0))
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(Imx,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(Jmx,4*(no1-1)+d)
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               SU(2)_x symmetry check
                 weight = gamma_M(a,b,3)*gamma_M(c,d,3)
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(I,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(J,4*(no1-1)+d)
                   if ((no>2 .and. no1<3) .or. (no<3 .and. no1>2)) weight = -weight
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   SxS_eq (imj,1,1) = SxS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               SU(2)_z symmetry check
                 weight = gamma_M(a,b,3)*gamma_M(c,d,3)
                 if ( abs(weight) > 0.01 ) then
                   I1 = Invlist(I,4*(no-1)+a)
                   if (no<3) then
                     I2 = Invlist(I,4*(no+1)+b)
                   else
                     I2 = Invlist(I,4*(no-3)+b)
                   endif
                   J1 = Invlist(J,4*(no1-1)+c)
                   if (no1 < 3) then
                     J2 = Invlist(J,4*(no1+1)+d)
                   else
                     J2 = Invlist(J,4*(no1-3)+d)
                   endif
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   SzS_eq (imj,1,1) = SzS_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               U(1)_1 symmetry check
                 weight = gamma_M(a,b,3)*gamma_M(c,d,3)
                 if ( abs(weight) > 0.01 ) then
                   signum = (-1)**no
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(I,4*(no-signum-1)+b)
                   signum = (-1)**no1
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(J,4*(no1-signum-1)+d)
      !                 if ((no>2 .and. no1<3) .or. (no<3 .and. n1>2)) weight = -weight
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   U11S_eq (imj,1,1) = U11S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
      !               U(1)_2 symmetry check
                 if (a.eq.b .and. c.eq.d) then
                   signum = (-1)**no
                   weight = cmplx(0.d0,dble(signum), kind(0.D0))
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(I,4*(no-signum-1)+b)
                   signum = (-1)**no1
                   weight = cmplx(0.d0,dble(signum), kind(0.D0))*weight
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(J,4*(no1-signum-1)+d)
      !                 if ((no>2 .and. no1<3) .or. (no<3 .and. n1>2)) weight = -weight
                   tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                      &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
                   U12S_eq (imj,1,1) = U12S_eq (imj,1,1)   +  weight*tmp
                 endif
                 
               enddo
               enddo
                enddo
                
                
  !               if (a==b) then
  !                 I1 = Invlist(I,4*(no-1)+a)
  !                 I2 = Invlist(I,4*(no-1)+a)
  !                 signum = 1
  !                 if ((no==2) .or. (no==4)) signum=-1
  !                 weight = cmplx(dble(signum),0.d0)
  !                 tmp =   GRC(I1,I2,1)* ZP*ZS
  !                 U1_eq0 (1) = U1_eq0 (1)   +  weight*tmp*0.5
  !               endif
                
              enddo
              enddo
            enddo
           enddo
         enddo
          endif
          
!           write(*,*) U1_eq0(1)

        end Subroutine Obser
!==========================================================        

        Subroutine  Pr_obs(LTAU)

          Use Print_bin_mod
#ifdef MPI
          Use mpi
#endif
          Implicit none

          Integer,  Intent(In) ::  Ltau

          Character (len=64) :: File_pr
          Complex   (Kind=Kind(0.d0)) :: Phase_bin
#ifdef MPI
          Integer        :: Isize, Irank, Ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
!!$#ifdef MPI
!!$          Write(6,*)  Irank, 'In Pr_obs', LTAU
!!$#else
!!$          Write(6,*)  'In Pr_obs', LTAU
!!$#endif

          File_pr ="ener"
          Call Print_scal(Obs_scal, Nobs, file_pr)

          Phase_bin = Obs_scal(8)/dble(Nobs)
          File_pr ="Den_eq"
          Call Print_bin(Den_eq, Den_eq0, Latt, Nobs, Phase_bin, file_pr)
          File_pr ="U1_eq"
          Call Print_bin(U1_eq, U1_eq0, Latt, Nobs, Phase_bin, file_pr)
          File_pr ="Spinz_eq"
          Call Print_bin(Spinz_eq, Spinz_eq0, Latt, Nobs, Phase_bin, file_pr)

          if (FlagSym == 1) then
         File_pr ="R_eq"
         Call Print_bin(R_eq, R_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="U1xy_eq"
         Call Print_bin(U1xy_eq, U1xy_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="Spinxy_eq"
         Call Print_bin(Spinxy_eq, Spinxy_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheckTR_eq"
         Call Print_bin(TRS_eq, TRS_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheckR_eq"
         Call Print_bin(RS_eq, RS_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheckC4_eq"
         Call Print_bin(C4S_eq, C4S_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheckPx_eq"
         Call Print_bin(PxS_eq, PxS_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheckSx_eq"
         Call Print_bin(SxS_eq, SxS_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheckSz_eq"
         Call Print_bin(SzS_eq, SzS_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheck1U1_eq"
         Call Print_bin(U11S_eq, U11S_eq0, Latt, Nobs, Phase_bin, file_pr)
         File_pr ="SymCheck2U1_eq"
         Call Print_bin(U12S_eq, U12S_eq0, Latt, Nobs, Phase_bin, file_pr)
       endif
          
          If (Ltau == 1) then
             Phase_tau = Phase_tau/dble(NobsT)
             File_pr = "Green_tau"
             Call Print_bin_tau(Green_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
             File_pr = "Den_tau"
             Call Print_bin_tau(Den_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
             File_pr = "U1_tau"
             Call Print_bin_tau(U1_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
             File_pr = "U1xy_tau"
             Call Print_bin_tau(U1xy_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
              File_pr = "U1xyG_tau"
              Call Print_bin_tau(U1xyG_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
             File_pr = "Spinz_tau"
             Call Print_bin_tau(Spinz_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
             File_pr = "Spinxy_tau"
             Call Print_bin_tau(Spinxy_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
          File_pr ="Den_sus"
          Call Print_bin(Den_sus, Den_sus0, Latt, NobsT, Phase_tau, file_pr)
          File_pr ="U1_sus"
          Call Print_bin(U1_sus, U1_sus0, Latt, NobsT, Phase_tau, file_pr)
          File_pr ="U1xy_sus"
          Call Print_bin(U1xy_sus, U1xy_sus0, Latt, NobsT, Phase_tau, file_pr)
           File_pr ="U1xyG_sus"
           Call Print_bin(U1xyG_sus, U1xyG_sus0, Latt, NobsT, Phase_tau, file_pr)
          File_pr ="Spinz_sus"
          Call Print_bin(Spinz_sus, Spinz_sus0, Latt, NobsT, Phase_tau, file_pr)
          File_pr ="Spinxy_sus"
          Call Print_bin(Spinxy_sus, Spinxy_sus0, Latt, NobsT, Phase_tau, file_pr)
          L_eq0=L_eq0*sqrt(beta)
          File_pr ="L_eq"
          Call Print_bin(L_eq, L_eq0, Latt, NobsT, Phase_tau, file_pr)
          endif
!!$#ifdef MPI
!!$          Write(6,*)  Irank, 'out Pr_obs', LTAU
!!$#else
!!$          Write(6,*)  'out Pr_obs', LTAU
!!$#endif
        end Subroutine Pr_obs
!==========================================================        

        Subroutine OBSERT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, tmp, DeltaI, DeltaJ, weight, weightbeta
          Integer :: IMJ, I1, I, no, J1, J, no1, I2, J2, signum,a,b,c,d

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then 
             Phase_tau = Phase_tau + ZS
             NobsT     = NobsT + 1
       
          do I=1,Latt%N
            do a=1,4
            do b=1,4
           if ( abs(gamma_M(a,b,3))>0.01 ) then
             do no=1,4
               I1 = InvList(I,4*(no-1)+a)
               J1 = InvList(I,4*(no-1)+b)
               DeltaI=cmplx(0.d0,0.d0,kind(0.d0))
               if (I1==J1) DeltaI=cmplx(1.d0,0.d0,kind(0.d0))
               L_eq0(1)  = L_eq0(1)  + Gamma_M(a,b,3)*(DeltaI - G00(j1,i1,1)) * ZP* ZS
             enddo
           endif
            enddo
            enddo
          enddo
          endif
          
          weightbeta = cmplx(dtau,0.d0,kind(0.d0))
          If (NT == 0 .or. NT==Ltrot) weightbeta=0.5*weightbeta 
          
          If ( N_FL == 1 ) then 
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I1 = 1,Ndim
          I  = List(I1,1)
          no = List(I1,2)
                Do J1 = 1,Ndim
             J  = List(J1,1)
             no1 = List(J1,2)
                   imj = latt%imj(I,J)
                   Green_tau(imj,nt+1,no,no1) = green_tau(imj,nt+1,no,no1)  +  Z * GT0(I1,J1,1) * ZP* ZS
                   
                   I2=I1
                   J2=J1
                   DeltaI=cmplx(0.d0,0.d0, kind(0.D0))
                   if (I1==I2 .and. nt==0) DeltaI=cmplx(1.d0,0.d0, kind(0.D0))
                   DeltaJ=cmplx(0.d0,0.d0, kind(0.D0))
                   if (J1==J2 .and. nt==0) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0))
                   tmp =  Z * ( (DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS          !
                   
                   Den_tau  (imj,nt+1,1,1) = Den_tau  (imj,nt+1,1,1)  +  tmp  - cmplx(0.25d0,0.d0, kind(0.D0))
                   Den_sus  (imj,1,1) = Den_sus  (imj,1,1)  +  weightbeta*(tmp  - cmplx(0.25d0,0.d0, kind(0.D0)))
                     
             weight=cmplx(1.d0,0.d0,kind(0.d0))
             if ( (no>=9 .and. no1<=8) .or. (no<=8 .and. no1>=9) ) weight=-weight
             Spinz_tau (imj,nt+1,1,1) = Spinz_tau (imj,nt+1,1,1)   +   weight * 0.25 * tmp
             Spinz_sus (imj,1,1) = Spinz_sus (imj,1,1)   +  weightbeta * weight * 0.25 * tmp
              
             signum = 1
             if (((no-1)/4+1==2) .or. ((no-1)/4+1==4)) signum=-1
             if (((no1-1)/4+1==2) .or. ((no1-1)/4+1 ==4)) signum=-signum
             weight = cmplx(dble(signum),0.d0, kind(0.D0))
             U1_tau (imj,nt+1,1,1) = U1_tau (imj,nt+1,1,1)   +  weight*tmp*0.25
             U1_sus (imj,1,1) = U1_sus (imj,1,1)   + weightbeta* weight*tmp*0.25

          enddo
           enddo
           
           do I=1,Latt%N
          do no=1,8
              do J=1,Latt%N
                imj = latt%imj(I,J)
                do no1=1,8
                 I1 = Invlist(I,no)
                 I2 = Invlist(I,no+8)
                 J1 = Invlist(J,no1+8)
                 J2 = Invlist(J,no1)
                 
                 DeltaI=cmplx(0.d0,0.d0, kind(0.D0))
                 if (I1==I2 .and. nt==0) DeltaI=cmplx(1.d0,0.d0, kind(0.D0))
                 DeltaJ=cmplx(0.d0,0.d0, kind(0.D0))
                 if (J1==J2 .and. nt==0) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0))
                 tmp =  Z * ((DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS
                 Spinxy_tau (imj,nt+1,1,1) = Spinxy_tau (imj,nt+1,1,1)   +  tmp
                 Spinxy_sus (imj,1,1) = Spinxy_sus (imj,1,1)   +   weightbeta*tmp
                 
                 if (no<=4) then
                   I1 = Invlist(I,no)
                   I2 = Invlist(I,no+4)
                 else
                   I1= Invlist(I,no+4)
                   I2 = Invlist(I,no+8)
                 endif
                 if (no1<=4) then
                   J1 = Invlist(J,no1+4)
                   J2 = Invlist(J,no1)
                 else
                   J1 = Invlist(J,no1+8)
                   J2 = Invlist(J,no1+4)
                 endif
                 
                 DeltaI=cmplx(0.d0,0.d0, kind(0.D0))
                 if (I1==I2 .and. nt==0) DeltaI=cmplx(1.d0,0.d0, kind(0.D0))
                 DeltaJ=cmplx(0.d0,0.d0, kind(0.D0))
                 if (J1==J2 .and. nt==0) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0))
                 tmp =  Z * ((DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS
                 U1xy_tau (imj,nt+1,1,1) = U1xy_tau (imj,nt+1,1,1)   +  tmp
                 U1xy_sus (imj,1,1) = U1xy_sus (imj,1,1)   +  weightbeta*tmp
                  U1xyG_tau (imj,nt+1,1,1) = U1xyG_tau (imj,nt+1,1,1)   +  (-1)**(no/2+no1/2+(no-1)/4+(no1-1)/4)*tmp
                  U1xyG_sus (imj,1,1) = U1xyG_sus (imj,1,1)   +  (-1)**(no/2+no1/2+(no-1)/4+(no1-1)/4)*weightbeta*tmp
                enddo
              enddo
              
                Enddo
             Enddo
             
             
         
         do I=1,Latt%N
           do J=1,Latt%N
          imj = latt%imj(I,J)
          do a=1,4
          do b=1,4
            if ( abs(gamma_M(a,b,3))>0.01 ) then
              do c=1,4
              do d=1,4
                weight = gamma_M(a,b,3)*Gamma_M(c,d,3)
                if ( abs(weight) > 0.01 ) then
               do no=1,4
                 do no1=1,4
                   I1 = Invlist(I,4*(no-1)+a)
                   I2 = Invlist(I,4*(no-1)+b)
                   J1 = Invlist(J,4*(no1-1)+c)
                   J2 = Invlist(J,4*(no1-1)+d)
                   DeltaI=cmplx(0.d0,0.d0, kind(0.D0))
                   if (I1==I2 .and. nt==0) DeltaI=cmplx(1.d0,0.d0, kind(0.D0))
                   DeltaJ=cmplx(0.d0,0.d0, kind(0.D0))
                   if (J1==J2 .and. nt==0) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0))
                   tmp =  Z * ((DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS
                   L_eq (imj,1,1) = L_eq (imj,1,1)   +  weight*tmp*weightbeta 
                 enddo
               enddo
                endif
              enddo
              enddo
            endif
         enddo
          enddo
           enddo
         enddo
          Endif
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
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
        End Subroutine Global_move
!========================================================================

!---------------------------------------------------------------------
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
        end Function Delta_S0_global

!---------------------------------------------------------------------
        Subroutine  Hamiltonian_set_random_nsigma
          
          ! The user can set the initial configuration
          
          Implicit none
          
          Integer :: I, nt
          
          Do nt = 1,Ltrot
             Do I = 1,Size(OP_V,1)
                nsigma(I,nt)  = 1
                if ( ranf_wrap()  > 0.5D0 ) nsigma(I,nt)  = -1
             enddo
          enddo
          
        end Subroutine Hamiltonian_set_random_nsigma

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
          Integer,    allocatable, INTENT(INOUT) :: Flip_list(:), Flip_value(:)
          Integer, INTENT(INOUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau

        end Subroutine Global_move_tau


    end Module Hamiltonian
