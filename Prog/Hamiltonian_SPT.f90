    !Model Hamiltonian for interaction-induced topological reduction
    Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      

      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Type (WaveFunction), dimension(:),   allocatable  :: WF_L
      Type (WaveFunction), dimension(:),   allocatable  :: WF_R
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot, Thtrot
      Logical              :: Projector
!>    Defines MPI communicator 
      Integer              :: Group_Comm


      
      ! What is below is  private 
      
      Type (Lattice),       private :: Latt
      Integer, parameter,   private :: Norb=16, Nobs_scal=12
      Integer, allocatable, private :: List(:,:), Invlist(:,:)
      Integer,              private :: L1, L2, FlagSym
      real (Kind=Kind(0.d0)),        private :: Ham_T, Ham_Vint, Ham_V2int, Ham_Lam
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Character (len=64),   private :: Model, Lattice_type, File1
      Complex (Kind=Kind(0.d0)),     private :: Gamma_M(4,4,5), Sigma_M(2,2,0:3)
      Complex (Kind=Kind(0.d0)),     private :: Gamma_13(4,4), Gamma_23(4,4), Gamma_45(4,4)


      ! Observables
      Integer,                       private :: Nobs
      Complex (Kind=Kind(0.d0)), allocatable, private :: obs_scal(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local_eq(:,:,:), Local_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local1_eq(:,:,:,:,:), Local1_eq0(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local2_eq(:,:,:,:,:), Local2_eq0(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local3_eq(:,:,:,:,:), Local3_eq0(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local4_eq(:,:,:,:,:), Local4_eq0(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local5_eq(:,:,:,:,:), Local5_eq0(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local6_eq(:,:,:,:,:), Local6_eq0(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Local7_eq(:,:,:,:,:), Local7_eq0(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Den_eq(:,:,:), Den_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  R_eq(:,:,:), R_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1_eq(:,:,:), U1_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1G_eq(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  L_eq(:,:,:), L_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1xy_eq(:,:,:), U1xy_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  U1xyG_eq(:,:,:), U1xyG_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  Spinz_eq(:,:,:), Spinz_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private ::  SpinzG_eq(:,:,:)
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
      
      Integer, private :: Idx(4,4), OpKind(16)
      Complex (Kind=Kind(0.d0)), private :: OpWeight(4,16)

      contains 

        Subroutine Ham_Set
#ifdef MPI
          Use mpi
#endif
          Implicit none

          integer :: ierr

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_SPT/  ham_T, Ham_Vint,  Ham_V2int, Ham_Lam,  Dtau, Beta, FlagSym, projector, Theta


#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
          If (Irank == 0 ) then
#endif
             File1 = "parameters"
#if defined(TEMPERING) 
             write(File1,'(A,I0,A)') "Temp_",igroup,"/parameters"
#endif
             OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
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
          FlagSym = 0
          Ltrot = nint(beta/dtau)
          Projector = .false.
          Theta = 0.d0
          Thtrot = 0
          ham_V2int=0.d0
          
! #if defined(MPI) && !defined(TEMPERING)
!           If (Irank == 0 ) then
! #endif
! #if defined(TEMPERING)
!              write(File1,'(A,I0,A)') "Temp_",Irank,"/parameters"
!              OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
! !              Global_moves =.false.
! #else
!              OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
! #endif
#if defined(MPI) 
          If (Irank_g == 0 ) then
#endif
             File1 = "parameters"
#if defined(TEMPERING) 
             write(File1,'(A,I0,A)') "Temp_",igroup,"/parameters"
#endif
             OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_SPT)
             CLOSE(5)
#if defined(MPI)
          endif

          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_Vint ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_V2int ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_Lam  ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Theta    ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Projector,1,MPI_LOGICAL,   0,Group_Comm,ierr)
          CALL MPI_BCAST(FlagSym  ,1,MPI_INTEGER, 0,Group_Comm,ierr)
#endif

          Call Ham_hop
          
          if (Projector) Call Ham_TrialWaveFunction
          
          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot
          
#if defined(TEMPERING)
           write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
           write(File1,'(A,I0)') "info"
#endif
           
#if defined(MPI) 
           If (Irank_g == 0)  then
#endif
             Open (Unit = 50,file=File1,status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'V             : ', Ham_Vint
             Write(50,*) 'V2            : ', Ham_V2int
             Write(50,*) 'Lambda        : ', Ham_Lam
             Write(50,*) 'FlagSym       : ', FlagSym
             close(50)
#if defined(MPI)
          endif
#endif
          call Ham_V
          
          Idx(1,1) = 1
          Idx(2,1) = 2
          Idx(3,1) = 3
          Idx(4,1) = 4
          
          Idx(1,2) = 4
          Idx(2,2) = 3
          Idx(3,2) = 2
          Idx(4,2) = 1
          
          Idx(1,3) = 2
          Idx(2,3) = 1
          Idx(3,3) = 4
          Idx(4,3) = 3
          
          Idx(1,4) = 3
          Idx(2,4) = 4
          Idx(3,4) = 1
          Idx(4,4) = 2
          
          OpWeight = cmplx(1.d0,0.d0,kind(0.d0))
          
          OpWeight(1,3) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(2,3) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(3,4) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(4,4) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(3,5) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(4,5) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(1,7) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(2,7) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(1,9) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(3,9) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(2,10) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(4,10) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(1,11) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(4,11) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(2,12) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(3,12) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(2,13) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(4,13) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(1,14) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(3,14) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(1,15) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(4,15) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(1,16) = cmplx(-1.d0,0.d0,kind(0.d0))
          OpWeight(4,16) = cmplx(-1.d0,0.d0,kind(0.d0))
          
          OpWeight(:,3) = OpWeight(:,3)*cmplx(0.d0,1.d0,kind(0.d0))
          OpWeight(:,7) = OpWeight(:,7)*cmplx(0.d0,1.d0,kind(0.d0))
          OpWeight(:,9) = OpWeight(:,9)*cmplx(0.d0,1.d0,kind(0.d0))
          OpWeight(:,11) = OpWeight(:,11)*cmplx(0.d0,1.d0,kind(0.d0))
          OpWeight(:,14) = OpWeight(:,14)*cmplx(0.d0,1.d0,kind(0.d0))
          OpWeight(:,16) = OpWeight(:,16)*cmplx(0.d0,1.d0,kind(0.d0))
          
          OpKind(1)=1
          OpKind(2)=2
          OpKind(3)=2
          OpKind(4)=1
          OpKind(5)=3
          OpKind(6)=4
          OpKind(7)=4
          OpKind(8)=3
          OpKind(9)=3
          OpKind(10)=1
          OpKind(11)=3
          OpKind(12)=1
          OpKind(13)=4
          OpKind(14)=2
          OpKind(15)=2
          OpKind(16)=4
          
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
!                 Do I = 1, Ndim/4
!                    Write(6,*) i,Op_T(nc,n)%E(i)
!                 enddo
             enddo
          enddo

          deallocate (Invlist_1) 

        end Subroutine Ham_hop
        
!===================================================================================           
        Subroutine Ham_TrialWaveFunction
        
          Implicit none
          COMPLEX(Kind=Kind(0.d0)), allocatable :: H0(:,:), U(:,:), Test(:,:)
          Real(Kind=Kind(0.d0)), allocatable :: En(:)
          
          Real(Kind=Kind(0.d0)), parameter :: phi=0.0
          COMPLEX(Kind=Kind(0.d0)) :: Z, alpha, beta, kek=cmplx(-0.01d0,0.0d0,kind(0.d0))
          Integer :: n, n_part, i, I1, I2, J1, J2, no, nc1, DI1, DI2, nc, Ihex(6)
          
          N_part=Ndim/2
          
          Allocate(WF_L(N_FL),WF_R(N_FL))
          do n=1,N_FL
            Call WF_alloc(WF_L(n),Ndim,N_part)
            Call WF_alloc(WF_R(n),Ndim,N_part)
          enddo
          
          Allocate(H0(Ndim/4,Ndim/4),U(Ndim/4,Ndim/4),En(Ndim/4))
          H0=Op_T(1,1)%O
          Call Diag(H0,U,En)
          
!           write(*,*) N_part, Ndim/2, (norb/2)*Latt%N
          write(*,*) 'Gap is', En(N_part/4+1)-En(N_part/4)
          do nc=1,4
            do I2=1,Ndim/8
            do I1=1,Ndim/4
              WF_L(1)%P(Op_T(nc,1)%P(I1),(nc-1)*Ndim/8+I2)=U(I1,I2)
            enddo
            enddo
          enddo
          WF_R(1)%P=WF_L(1)%P
          
          Allocate(Test(N_part,N_part))
          alpha=1.d0
          beta=0.d0
          CALL ZGEMM('C','N',N_part,N_part,ndim,alpha,WF_L(1)%P(1,1),Ndim,WF_R(1)%P(1,1),Ndim,beta,Test(1,1),N_part)
          Z = Det_C(Test,n_part)
          write(*,*) Z
          
!           do I2=1,N_part
!             write(*,*) sum(abs(WF_L(1)%P(:,I2))), sum(abs(WF_R(1)%P(:,I2)))
!           enddo
          
          Deallocate(H0, U, En)
          
        end Subroutine Ham_TrialWaveFunction
        
!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, nth, n, n1, n2, n3, n4, I, I1, I2, J,  Ix, Iy, nc, no,no1, ns, npm 
          Integer :: nxy, noff
          Real (Kind=Kind(0.d0)) :: X_p(2), X1_p(2), X2_p(2), X, XJ, Xpm

          Complex (Kind=Kind(0.d0)) :: Ps(4,4,2), Ps_G5(4,4,2), Tmp(4,4), Z
          Complex (Kind=Kind(0.d0)) :: Sx(16,16,2,2), Sy(16,16,2,2), Sz(16,16,2)


          Ps = cmplx(0.d0, 0.d0, kind(0.D0))
          Call mmult (Tmp,Gamma_M(:,:,3), Gamma_M(:,:,4) )
          
! !           Do n=1,5
!             write(*,*) "Gamma34"
!             do no=1,4
!             do no1=1,4
!               write(*,*) Tmp(no,no1)
!             enddo
!             write(*,*)
!             enddo
!             write(*,*)
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
          
          ! In this basis we have Ps_G5(:,:,1)=diag(0,-1,0,1) and  Ps_G5(:,:,2)=diag(1,0,-1,0)
          Do ns = 1,2
             Call mmult ( Ps_G5(:,:,ns), Ps(:,:,ns), Gamma_M(:,:,5) )
          enddo
          
!           Do ns=1,2
!             write(*,*) "PsG5 ", ns
!             do no=1,4
!             do no1=1,4
!               write(*,*) Ps_G5(no,no1,ns)
!             enddo
!             write(*,*)
!             enddo
!             write(*,*)
!           enddo
      
          Sx = 0.d0
          Sy = 0.d0
          Sz = 0.d0
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
             Do no = 1,4
                do no1 = 1,4
                  Sz(no    , no1     ,ns) =  Ps(no,no1,ns)
                  Sz(no +4 , no1 + 4 ,ns) =  -Ps(no,no1,ns)
                  Sz(no +8 , no1 + 8 ,ns) =  Ps(no,no1,ns)
                  Sz(no+12 , no1 + 12,ns) =  -Ps(no,no1,ns)
                enddo
             enddo
          enddo


          ! Number of opertors 8 per unit cell
          Allocate( Op_V(10*Latt%N,N_FL) )
          do nf = 1,N_FL
             do i  = 1, 10*Latt%N
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
!                          if (nxy == 1) write(*,*) nc, ' is a x-type vertex.'
!                          if (nxy == 2) write(*,*) nc, ' is a y-type vertex.' 
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
             do ns = 1,2
                noff=2
                if (ns==2) noff=1
                  Do i = 1,Latt%N
                      nc = nc + 1 
                      Do no = 1,Norb/2
                        Op_V(nc,nf)%P(no)   = Invlist(I,2*(no-1)+noff)  
                      enddo
                      Do no = 1,Norb/2
                        Do no1 = 1,Norb/2
                            Op_V(nc,nf)%O(no,no1) = Sz(2*(no-1)+noff,2*(no1-1)+noff,ns)
                        Enddo
                      Enddo
                      Op_V(nc,nf)%g = SQRT(CMPLX(-DTAU*Ham_V2int, 0.D0, kind(0.D0))) 
                      Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                      Op_V(nc,nf)%type   = 2
                      Call Op_set( Op_V(nc,nf) )
                      ! The operator reads: 
                      ! g*s*( c^{dagger} O c   - alpha ))
                      ! with s the HS field.
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
          Allocate ( Obs_scal(Nobs_scal) )
          Allocate ( Local_eq(Latt%N,Norb,Norb), Local_eq0(Norb) )
          
          Allocate ( Local1_eq(Latt%N,3,3,3,3), Local1_eq0(3,3,3) )
          Allocate ( Local2_eq(Latt%N,1,1,3,3), Local2_eq0(1,3,3) )
          Allocate ( Local3_eq(Latt%N,1,1,3,3), Local3_eq0(1,3,3) )
          Allocate ( Local4_eq(Latt%N,3,3,3,3), Local4_eq0(3,3,3) )
          Allocate ( Local5_eq(Latt%N,1,1,3,3), Local5_eq0(1,3,3) )
          Allocate ( Local6_eq(Latt%N,1,1,3,3), Local6_eq0(1,3,3) )
          Allocate ( Local7_eq(Latt%N,2,2,3,3), Local7_eq0(2,3,3) )
          
          Allocate ( Den_eq(Latt%N,1,1), Den_eq0(1) ) 
          Allocate ( U1_eq(Latt%N,1,1), U1_eq0(1) )
          Allocate ( U1G_eq(Latt%N,1,1) )
          Allocate ( U1xy_eq(Latt%N,1,1), U1xy_eq0(1) )
          Allocate ( U1xyG_eq(Latt%N,1,1), U1xyG_eq0(1) )
          Allocate ( L_eq(Latt%N,1,1), L_eq0(1) )
          Allocate ( Spinz_eq(Latt%N,1,1), spinz_eq0(1) )
          Allocate ( SpinzG_eq(Latt%N,1,1) )
          Allocate ( Spinxy_eq(Latt%N,1,1), spinxy_eq0(1) ) 
          
          if (FlagSym ==1) then
              Allocate ( R_eq(Latt%N,1,1), R_eq0(1) ) 
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
             Allocate ( Green_tau(Latt%N,Ltrot+1-2*Thtrot,Norb,Norb), Den_tau(Latt%N,Ltrot+1-2*Thtrot,1,1) )
             Allocate ( U1_tau(Latt%N,Ltrot+1-2*Thtrot,1,1))
             Allocate ( U1xy_tau(Latt%N,Ltrot+1-2*Thtrot,1,1), U1xyG_tau(Latt%N,Ltrot+1-2*Thtrot,1,1) )
             Allocate ( Spinz_tau(Latt%N,Ltrot+1-2*Thtrot,1,1), Spinxy_tau(Latt%N,Ltrot+1-2*Thtrot,1,1) )
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
          Obs_scal  = 0.d0
          
          Local_eq    = 0.d0
          Local_eq0   = 0.d0
          
          Local1_eq    = 0.d0
          Local1_eq0   = 0.d0
          Local2_eq    = 0.d0
          Local2_eq0   = 0.d0
          Local3_eq    = 0.d0
          Local3_eq0   = 0.d0
          Local4_eq    = 0.d0
          Local4_eq0   = 0.d0
          Local5_eq    = 0.d0
          Local5_eq0   = 0.d0
          Local6_eq    = 0.d0
          Local6_eq0   = 0.d0
          Local7_eq    = 0.d0
          Local7_eq0   = 0.d0
          
          Den_eq    = 0.d0
          Den_eq0   = 0.d0
          Spinz_eq    = 0.d0
          SpinzG_eq    = 0.d0
          Spinz_eq0   = 0.d0
          U1_eq    = 0.d0
          U1G_eq    = 0.d0
          U1_eq0   = 0.d0
          U1xy_eq    = 0.d0
          U1xy_eq0   = 0.d0
          U1xyG_eq    = 0.d0
          U1xyG_eq0   = 0.d0
          Spinxy_eq    = 0.d0
          Spinxy_eq0   = 0.d0
          L_eq    = 0.d0
          L_eq0   = 0.d0
          
          if (FlagSym ==1 ) then
              R_eq    = 0.d0
              R_eq0   = 0.d0
              
              TRS_eq    = 0.d0
              TRS_eq0   = 0.d0
              PHS_eq    = 0.d0
              PHS_eq0   = 0.d0
              RS_eq    = 0.d0
              RS_eq0   = 0.d0
              C4S_eq    = 0.d0
              C4S_eq0   = 0.d0
              PxS_eq    = 0.d0
              PxS_eq0   = 0.d0
              SxS_eq    = 0.d0
              SxS_eq0   = 0.d0
              SzS_eq    = 0.d0
              SzS_eq0   = 0.d0
              U11S_eq    = 0.d0
              U11S_eq0   = 0.d0
              U12S_eq    = 0.d0
              U12S_eq0   = 0.d0
          endif

          If (Ltau == 1) then
             NobsT = 0
             Phase_tau = 0.d0
             Green_tau = 0.d0
             Den_tau = 0.d0
             U1_tau = 0.d0
             U1xy_tau = 0.d0
             U1xyG_tau = 0.d0
             Spinz_tau = 0.d0
             Spinxy_tau = 0.d0
          
             Den_eq    = 0.d0
             Den_eq0   = 0.d0
             U1_sus    = 0.d0
             U1_sus0   = 0.d0
             U1xy_sus    = 0.d0
             U1xy_sus0   = 0.d0
             U1xyG_sus    = 0.d0
             U1xyG_sus0   = 0.d0
             Spinz_sus    = 0.d0
             Spinz_sus0   = 0.d0
             Spinxy_sus    = 0.d0
             Spinxy_sus0   = 0.d0
          endif

        end Subroutine Init_obs
        
        pure Subroutine U1Structure(I,J,a,b,c,d,res,GR,GRC)
          Implicit None
          Complex (Kind=Kind(0.d0)), Intent(INOUT) :: res(3)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: GR(:,:,:),GRC(:,:,:)
          Integer, INTENT(IN)          :: I,J,a,b,c,d
          
          Integer :: no1, no2, I1, I2, J1, J2
          Complex (Kind=Kind(0.d0)) :: tmp
          
          res=0.d0
          
          do no1=0,4,4
            do no2=0,4,4
              I1 = InvList(I,no1+a)
              I2 = InvList(I,no1+b)
              J1 = InvList(J,no2+c)
              J2 = InvList(J,no2+d)
              
              tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                   &     GRC(I1,I2,1) * GRC(J1,J2,1)         )
              
              res(1)=res(1)+tmp
              if (no1==no2) then
                res(2)=res(2)+tmp
              else
                res(2)=res(2)-tmp
              endif
            enddo
          enddo
          
          do no1=0,4,4
            no2=4-no1
            I1 = InvList(I,no1+a)
            I2 = InvList(I,no2+b)
            J1 = InvList(J,no2+c)
            J2 = InvList(J,no1+d)
            
            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         )
            
            res(3)=res(3)+tmp
          enddo
          
        end Subroutine U1Structure
        
        pure Subroutine SzStructure(I,J,a,b,c,d,res,GR,GRC)
          Implicit None
          Complex (Kind=Kind(0.d0)), Intent(INOUT) :: res(3,3)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: GR(:,:,:),GRC(:,:,:)
          Integer, INTENT(IN)          :: I,J,a,b,c,d
          
          Integer :: no1, no2, no
          Complex (Kind=Kind(0.d0)) :: tmp(3)
          
          res=0.d0
          
          do no1=0,8,8
            do no2=0,8,8
              
              call U1Structure(i,j,no1+a,no1+b,no2+c,no2+d,tmp,GR,GRC)
                
              do no=1,3
                res(1,no)=res(1,no)+tmp(no)
                if (no1==no2) then
                  res(2,no)=res(2,no)+tmp(no)
                else
                  res(2,no)=res(2,no)-tmp(no)
                endif
              enddo
            enddo
          enddo
          
          do no1=0,8,8
            no2=8-no1
            
            call U1Structure(I,J,no1+a,no2+b,no2+c,no1+d,tmp,GR,GRC)
            
            do no=1,3
              res(3,no)=res(3,no)+tmp(no)
            enddo
          enddo
          
        end Subroutine SzStructure
        
!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, Res1(3,3,4,4)
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, ZL, Z, ZP,ZS, weight, tmp, ZU1, ZU1xG, ZU1yG, tmp1, tmp2, FdU
          Integer :: I,J, no,no1, n, n1, imj, nf, I1, I2, J1, J2, Nc, Ix, Iy, Jx, Jy, Imx, Imy, Jmx, Jmy
          Integer :: a, b, c, d, signum, K, K1, L ,L1, nf1,id1 ,id2
          
          Real (Kind=Kind(0.d0)) :: G(4,4), X, FI, FJ
          
          Nobs = Nobs + 1
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          

          Do nf = 1,N_FL
!$OMP parallel do default(shared) private(I,J,ZK)
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I,J,nf)  = - GR(J,I,nf)
                Enddo
                GRC(I,I,nf)  = GRC(I,I,nf) + 1
             Enddo
!$OMP end parallel do
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >
          ! Compute scalar observables. 

          Zkin = 0.d0

          Nc = Size( Op_T,1)
          Do nf = 1,N_FL
!$OMP parallel do default(shared) private(n,J,J1,I,I1) reduction(+:Zkin)
             Do n = 1,Nc
                Do J = 1,Op_T(n,nf)%N
                   J1 = Op_T(n,nf)%P(J)
                   DO I = 1,Op_T(n,nf)%N
                      I1 = Op_T(n,nf)%P(I)
                      Zkin  = Zkin  + Op_T(n,nf)%O(i,j)*Grc(i1,j1,nf) 
                   Enddo
                ENddo
             Enddo
!$OMP end parallel do
          Enddo
          Zkin = Zkin*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))
          
          ZL = 0.d0
          
!           Nc = Size( Op_T,1)
          Do nf = 1,N_FL
!              Do n = 1,Nc
!$OMP parallel do default(shared) private(I,no,a,b,I1,J1) reduction(+:ZL)
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
!$OMP end parallel do
!              Enddo
          Enddo
          ZL = ZL*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))

          Zrho = 0.d0!cmplx(0.d0, 0.d0, kind(0.D0))
          ZU1  = 0.d0
          Do nf = 1,N_FL
!$OMP parallel do default(shared) private(I,no,signum,tmp,weight) reduction(+:Zrho,ZU1)
             Do I = 1,Ndim
                tmp=GRC(I,I,nf)
                Zrho = Zrho + tmp
                no = List(I,2)
                signum = 1
                if (((no-1)/4+1==2) .or. ((no-1)/4+1==4)) signum=-1
                weight = cmplx(dble(signum),0.d0, kind(0.D0))
                ZU1 = ZU1  +  weight*tmp*0.5d0
             enddo
!$OMP end parallel do
          enddo
          Zrho = Zrho*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))
          ZU1 = ZU1*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))

          ZU1xG = 0.d0
          ZU1yG = 0.d0
!$OMP parallel do default(shared) private(I,no,I1,I2,J1,J2,tmp,weight) reduction(+:ZU1xG,ZU1yG)
          do I=1,Latt%N
            do no=1,8
                if (no<=4) then
                  I1 = Invlist(I,no)
                  I2 = Invlist(I,no+4)
                  J1 = Invlist(I,no+4)
                  J2 = Invlist(I,no)
                else
                  I1= Invlist(I,no+4)
                  I2 = Invlist(I,no+8)
                  J1 = Invlist(I,no+8)
                  J2 = Invlist(I,no+4)
                endif
!                 tmp = GRC(I1,I2,nf)
                tmp1 = GRC(I1,I2,1)
                tmp2 = GRC(J1,J2,1)
                tmp = tmp1+tmp2
                weight = (-1)**(no/2+(no-1)/4)
                ZU1xG = ZU1xG   +  weight*tmp
                tmp = cmplx(0.d0,1.d0,kind(1.d0))*(tmp1-tmp2)
                ZU1yG = ZU1yG   +  weight*tmp
                
            enddo
          enddo
!$OMP end parallel do
          ZU1xG = ZU1xG*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))
          ZU1yG = ZU1yG*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))
!           ZU1xyG = cmplx( 0.5d0*(real(ZU1xyG) + aimag(ZU1xyG)), 0.d0 , kind(0.D0))

          ZPot = 0.d0
          FdU=0.d0
          Nc = Size( Op_V,1)
          Do nf = 1,N_FL
!$OMP parallel do default(shared) private(n,weight,J,J1,I,I1,K,K1,L,L1,tmp) reduction(+:ZPot)
             Do n = 1,8*Latt%N!Nc
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
!                 write(*,*) Zpot
             Enddo
!$OMP end parallel do
!$OMP parallel do default(shared) private(n,weight,J,J1,I,I1,K,K1,L,L1,tmp) reduction(+:FdU)
             Do n = 8*Latt%N+1,Nc
                weight=1.d0
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
                                FdU  = FdU  + weight*Op_V(n,nf)%O(i,j)*Op_V(n,nf)%O(k,l)*tmp
                              endif
                            Enddo
                          ENddo
                      endif
                   Enddo
                ENddo
!                 write(*,*) Zpot
             Enddo
!$OMP end parallel do
          Enddo
          FdU = FdU*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))
          ZPot = ZPot*N_SUN!cmplx( dble(N_SUN), 0.d0 , kind(0.D0))

          Obs_scal(1) = Obs_scal(1) + zrho * ZP*ZS
          Obs_scal(2) = Obs_scal(2) + zkin*Ham_T * ZP*ZS
          Obs_scal(3) = Obs_scal(3) + (Zpot*Ham_Vint + FdU*Ham_V2int) * ZP*ZS
          Obs_scal(4) = Obs_scal(4) + (zkin*Ham_T +  Zpot*Ham_Vint + FdU*Ham_V2int)*ZP*ZS
          Obs_scal(5) = Obs_scal(5) + (zkin -  Zpot)*ZP*ZS
          Obs_scal(6) = Obs_scal(6) + (Zpot)*ZP*ZS
          Obs_scal(7) = Obs_scal(7) + (ZL)*ZP*ZS
          Obs_scal(8) = Obs_scal(8) + (ZU1)*ZP*ZS
          Obs_scal(9) = Obs_scal(9) + (ZU1xG)*ZP*ZS
          Obs_scal(10) = Obs_scal(10) + ZU1yG*ZP*ZS
          Obs_scal(11) = Obs_scal(11) + FdU*ZP*ZS
          Obs_scal(Nobs_scal) = Obs_scal(Nobs_scal) + ZS
          ! You will have to allocate more space if you want to include more  scalar observables.
          ! the last one has to be the phase!!!
          
          res1=0.d0
          Do I = 1, Latt%N
            do J = 1, Latt%N
                imj = latt%imj(I,J)
                
                do no = 1, 4
                do no1 = 1, 4
                  Call SzStructure(I,J,no,Idx(no,1),no1,Idx(no1,1),res1(:,:,1,1),GR,GRC)
                  Call SzStructure(I,J,no,Idx(no,1),no1,Idx(no1,2),res1(:,:,1,2),GR,GRC)
                  Call SzStructure(I,J,no,Idx(no,2),no1,Idx(no1,1),res1(:,:,2,1),GR,GRC)
                  Call SzStructure(I,J,no,Idx(no,2),no1,Idx(no1,2),res1(:,:,2,2),GR,GRC)
                  
                  Call SzStructure(I,J,no,Idx(no,3),no1,Idx(no1,3),res1(:,:,3,3),GR,GRC)
                  Call SzStructure(I,J,no,Idx(no,3),no1,Idx(no1,4),res1(:,:,3,4),GR,GRC)
                  Call SzStructure(I,J,no,Idx(no,4),no1,Idx(no1,3),res1(:,:,4,3),GR,GRC)
                  Call SzStructure(I,J,no,Idx(no,4),no1,Idx(no1,4),res1(:,:,4,4),GR,GRC)
                  
                  res1=res1*ZP*ZS
                  
                  do id1=1,3
                  do id2=1,3
                    do a=1,3
                    do b=1,3
                      Local1_eq(imj,a,b,id1,id2)=Local1_eq(imj,a,b,id1,id2)+OpWeight(no,a)*OpWeight(no1,b)&
                          &*res1(id1,id2,OpKind(a),OpKind(b))
                      Local4_eq(imj,a,b,id1,id2)=Local4_eq(imj,a,b,id1,id2)+OpWeight(no,5+a)*OpWeight(no1,5+b)&
                          &*res1(id1,id2,OpKind(5+a),OpKind(5+b))
                    enddo
                    enddo
                    
                    Local2_eq(imj,1,1,id1,id2)=Local2_eq(imj,1,1,id1,id2)+OpWeight(no,4)*OpWeight(no1,4)&
                          &*res1(id1,id2,OpKind(4),OpKind(4))
                    
                    Local3_eq(imj,1,1,id1,id2)=Local3_eq(imj,1,1,id1,id2)+OpWeight(no,5)*OpWeight(no1,5)&
                          &*res1(id1,id2,OpKind(5),OpKind(5))
                    
                    Local5_eq(imj,1,1,id1,id2)=Local5_eq(imj,1,1,id1,id2)+OpWeight(no,9 )*OpWeight(no1,9 )&
                          &*res1(id1,id2,OpKind(9 ),OpKind(9 ))
                    Local5_eq(imj,1,1,id1,id2)=Local5_eq(imj,1,1,id1,id2)+OpWeight(no,10)*OpWeight(no1,10)&
                          &*res1(id1,id2,OpKind(10),OpKind(10))
                    
                    Local6_eq(imj,1,1,id1,id2)=Local6_eq(imj,1,1,id1,id2)+OpWeight(no,11)*OpWeight(no1,11)&
                          &*res1(id1,id2,OpKind(11),OpKind(11))
                    Local6_eq(imj,1,1,id1,id2)=Local6_eq(imj,1,1,id1,id2)+OpWeight(no,12)*OpWeight(no1,12)&
                          &*res1(id1,id2,OpKind(12),OpKind(12))
                    
                    do a=1,2
                    do b=1,2
                      Local7_eq(imj,a,b,id1,id2)=Local7_eq(imj,a,b,id1,id2)+OpWeight(no,13+a)*OpWeight(no1,13+b)&
                          &*res1(id1,id2,OpKind(13+a),OpKind(13+b))
                      Local7_eq(imj,a,b,id1,id2)=Local7_eq(imj,a,b,id1,id2)+OpWeight(no,13+3*(a-1))*OpWeight(no1,13+3*(b-1))&
                          &*res1(id1,id2,OpKind(13+3*(a-1)),OpKind(13+3*(b-1)))
                    enddo
                    enddo
                  enddo
                  enddo
                  
                  
                enddo
                enddo
                
                do no = 1, Norb
                do no1 = 1, Norb
                  I1=Invlist(I,no)
                  j1=Invlist(I,no)
                  k1=Invlist(j,no1)
                  l1=Invlist(j,no1)
                  tmp =  (   GRC(I1,L1,1) * GR (J1,K1,1)      +  &
                      &     GRC(I1,J1,1) * GRC(K1,L1,1)         )
                  Local_eq (imj,no,no1) = Local_eq (imj,no,no1)   +   tmp*ZP*ZS
                enddo
                enddo
            enddo
            Local1_eq0(1,1,1)=Local1_eq0(1,1,1)+8.d0
            do no1=0,12,4
              do no=1,4
                I1=Invlist(I,no+no1)
                I2=Invlist(I,idx(no,OpKind(5))+no1)
                tmp = GRC(I1,I2,1)*ZP*ZS
                Local3_eq0(1,1,1)=Local3_eq0(1,1,1)+tmp*OpWeight(no,5)
              enddo
            enddo
            do no = 1, Norb
              I1=Invlist(I,no)
              tmp = GRC(I1,I1,1)*ZP*ZS
              Local_eq0(no)=Local_eq0(no)+tmp
            enddo
          enddo
          
!$OMP parallel do default(shared) private(I1,I,no,J1,J,no1,imj,tmp,weight,signum)
          DO I1 = 1,Ndim
             I  = List(I1,1)
             no = List(I1,2)
             DO J1 = 1, Ndim
                J = List(J1,1)
                no1 = list(J1,2)
                imj = latt%imj(I,J)
                
                tmp =  (   GRC(I1,J1,1) * GR (I1,J1,1)      +  &
                      &     GRC(I1,I1,1) * GRC(J1,J1,1)         ) *ZP*ZS
                      
!$OMP CRITICAL
                DEN_Eq (imj,1,1) = DEN_Eq (imj,1,1)   +  tmp
!$OMP END CRITICAL
                     
                weight=cmplx(1.d0,0.d0, kind(0.D0))
                if ( (no>=9 .and. no1<=8) .or. (no<=8 .and. no1>=9) ) weight=-weight
!$OMP CRITICAL
                Spinz_eq (imj,1,1) = Spinz_eq (imj,1,1)   +   weight * 0.25 * tmp
                SpinzG_eq (imj,1,1) = SpinzG_eq (imj,1,1)   +   weight * 0.25 * tmp*(-1)**(no/2+no1/2)
!$OMP END CRITICAL
                
                signum = 1
                if (((no-1)/4+1==2) .or. ((no-1)/4+1==4)) signum=-1
                if (((no1-1)/4+1==2) .or. ((no1-1)/4+1 ==4)) signum=-signum
                weight = cmplx(dble(signum),0.d0, kind(0.D0))
!$OMP CRITICAL
                U1_eq (imj,1,1) = U1_eq (imj,1,1)   +  weight*tmp*0.25
                U1G_eq (imj,1,1) = U1G_eq (imj,1,1)   +  weight*tmp*0.25*(-1)**(no/2+no1/2)
!$OMP END CRITICAL

             enddo
             tmp=GRC(I1,I1,1)
!$OMP CRITICAL
             Den_eq0(1) = Den_eq0(1) +   tmp*ZP*ZS 
!$OMP END CRITICAL
                     
!              weight=cmplx(1.d0,0.d0, kind(0.D0))
!              if ( (no>=9) ) weight=-weight
!              Spinz_eq0 (1) = Spinz_eq0 (1)   +   weight * 0.5 * tmp
!             
!              signum = 1
!              if (((no-1)/4+1==2) .or. ((no-1)/4+1==4)) signum=-1
!              weight = cmplx(dble(signum),0.d0, kind(0.D0))
!              U1_eq0 (1) = U1_eq0 (1)   +  weight*tmp*0.5
          enddo
!$OMP end parallel do
            
!           do I=1,Latt%N
!             do J=1,Latt%N
!               imj = latt%imj(I,J)
!               do a=1,4
!               do b=1,4
!                 if ( abs(gamma_M(a,b,3))>0.01 ) then
!                   do c=1,4
!                   do d=1,4
!                     weight = gamma_M(a,b,3)*Gamma_M(c,d,3)
!                     if ( abs(weight) > 0.01 ) then
!                       do no=1,4
!                         do no1=1,4
!                           I1 = Invlist(I,4*(no-1)+a)
!                           I2 = Invlist(I,4*(no-1)+b)
!                           J1 = Invlist(J,4*(no1-1)+c)
!                           J2 = Invlist(J,4*(no1-1)+d)
!                           tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
!                                 &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!                           L_eq (imj,1,1) = L_eq (imj,1,1)   +  weight*tmp 
!                         enddo
!                       enddo
!                     endif
!                   enddo
!                   enddo
!                 endif
!               enddo
!               enddo
!             enddo
! !           enddo
!           
!           do I=1,Latt%N
!             do a=1,4
!             do b=1,4
!               if ( abs(gamma_M(a,b,3))>0.01 ) then
!                 do no=1,4
!                   I1 = InvList(I,4*(no-1)+a)
!                   J1 = InvList(I,4*(no-1)+b)
!                   L_eq0(1)  = L_eq0(1)  + Gamma_M(a,b,3)*Grc(i1,j1,1)
!                 enddo
!               endif
!             enddo
!             enddo
!           enddo
          
!$OMP parallel do default(shared) private(I,I1,I2,no,J,J1,J2,no1,imj,tmp,weight,signum)
          do I=1,Latt%N
            do no=1,8
                I1 = Invlist(I,no)
                I2 = Invlist(I,no+8)
                do J=1,Latt%N
                  imj = latt%imj(I,J)
                  do no1=1,8
                    J1 = Invlist(J,no1+8)
                    J2 = Invlist(J,no1)
                    
                    tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                        &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                    Spinxy_eq (imj,1,1) = Spinxy_eq (imj,1,1)   +  tmp
!$OMP END CRITICAL
                  enddo
                enddo
                
            enddo
          enddo
!$OMP end parallel do
          
!$OMP parallel do default(shared) private(I,I1,I2,no,J,J1,J2,no1,imj,tmp,weight,signum)
          do I=1,Latt%N
            do no=1,8
                if (no<=4) then
                  I1 = Invlist(I,no)
                  I2 = Invlist(I,no+4)
                else
                  I1= Invlist(I,no+4)
                  I2 = Invlist(I,no+8)
                endif
                do J=1,Latt%N
                  imj = latt%imj(I,J)
                  do no1=1,8
                    if (no1<=4) then
                      J1 = Invlist(J,no1+4)
                      J2 = Invlist(J,no1)
                    else
                      J1 = Invlist(J,no1+8)
                      J2 = Invlist(J,no1+4)
                    endif
                    
                    tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                        &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                    U1xy_eq (imj,1,1) = U1xy_eq (imj,1,1)   +  tmp
                    U1xyG_eq (imj,1,1) = U1xyG_eq (imj,1,1)   +  (-1)**(no/2+no1/2+(no-1)/4+(no1-1)/4)*tmp
!$OMP END CRITICAL
                  enddo
                enddo
                
            enddo
          enddo
!$OMP end parallel do
            
          if (FlagSym ==1) then
!$OMP parallel do default(shared) private(I,I1,I2,Ix,Iy,Imx,Imy,no,J,J1,J2,jx,jy,jmx,jmy,no1,imj,a,b,c,d,tmp,weight,signum)
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
                          
      !                     R correlation
                          weight = -gamma_45(a,b)*Gamma_45(c,d)
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(I,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(J,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            R_eq (imj,1,1) = R_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     TR symmetry check
                          weight = gamma_13(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = gamma_13(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_13(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_13(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = gamma_23(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = gamma_23(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_23(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_23(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_13(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_13(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = gamma_13(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = gamma_13(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_23(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = -gamma_23(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = gamma_23(a,b)*Gamma_13(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                          weight = gamma_23(a,b)*Gamma_23(c,d)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            TRS_eq (imj,1,1) = TRS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     R symmetry check X°\dag gamma_4 X as correlation
                          weight = gamma_M(a,b,4)*Gamma_M(c,d,4)
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(I,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(J,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            RS_eq (imj,1,1) = RS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     C4 symmetry check
      !                     11
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     12
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     13
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     14
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     21
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     22
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     23
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     24
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     31
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     32
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     33
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     34
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     41
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     42
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     43
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     44
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            C4S_eq (imj,1,1) = C4S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     Px symmetry check
      !                     11
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     12
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     13
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     14
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Iy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     21
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     22
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     23
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     24
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Ix,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     31
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     32
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     33
                          weight = -gamma_M(a,b,1)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     34
                          weight = gamma_M(a,b,1)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imy,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     41
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     42
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     43
                          weight = gamma_M(a,b,2)*Gamma_M(c,d,1)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmy,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     44
                          weight = -gamma_M(a,b,2)*Gamma_M(c,d,2)*0.25d0
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(Imx,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(Jmx,4*(no1-1)+d)
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            PxS_eq (imj,1,1) = PxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     SU(2)_x symmetry check
                          weight = gamma_M(a,b,3)*gamma_M(c,d,3)
                          if ( abs(weight) > 0.01 ) then
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(I,4*(no-1)+b)
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(J,4*(no1-1)+d)
                            if ((no>2 .and. no1<3) .or. (no<3 .and. no1>2)) weight = -weight
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            SxS_eq (imj,1,1) = SxS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     SU(2)_z symmetry check
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
!$OMP CRITICAL
                            SzS_eq (imj,1,1) = SzS_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     U(1)_1 symmetry check
                          weight = gamma_M(a,b,3)*gamma_M(c,d,3)
                          if ( abs(weight) > 0.01 ) then
                            signum = (-1)**no
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(I,4*(no-signum-1)+b)
                            signum = (-1)**no1
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(J,4*(no1-signum-1)+d)
      !                       if ((no>2 .and. no1<3) .or. (no<3 .and. n1>2)) weight = -weight
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            U11S_eq (imj,1,1) = U11S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
      !                     U(1)_2 symmetry check
                          if (a.eq.b .and. c.eq.d) then
                            signum = (-1)**no
                            weight = cmplx(0.d0,dble(signum), kind(0.D0))
                            I1 = Invlist(I,4*(no-1)+a)
                            I2 = Invlist(I,4*(no-signum-1)+b)
                            signum = (-1)**no1
                            weight = cmplx(0.d0,dble(signum), kind(0.D0))*weight
                            J1 = Invlist(J,4*(no1-1)+c)
                            J2 = Invlist(J,4*(no1-signum-1)+d)
      !                       if ((no>2 .and. no1<3) .or. (no<3 .and. n1>2)) weight = -weight
                            tmp =  (   GRC(I1,J2,1) * GR (I2,J1,1)      +  &
                                  &     GRC(I1,I2,1) * GRC(J1,J2,1)         ) * ZP*ZS
!$OMP CRITICAL
                            U12S_eq (imj,1,1) = U12S_eq (imj,1,1)   +  weight*tmp
!$OMP END CRITICAL
                          endif
                          
                        enddo
                        enddo
                      enddo
                      
                      
  !                     if (a==b) then
  !                       I1 = Invlist(I,4*(no-1)+a)
  !                       I2 = Invlist(I,4*(no-1)+a)
  !                       signum = 1
  !                       if ((no==2) .or. (no==4)) signum=-1
  !                       weight = cmplx(dble(signum),0.d0)
  !                       tmp =   GRC(I1,I2,1)* ZP*ZS
  !                       U1_eq0 (1) = U1_eq0 (1)   +  weight*tmp*0.5
  !                     endif
                      
                    enddo
                    enddo
                  enddo
              enddo
            enddo
!$OMP end parallel do
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
          Integer ::  no, no1, a,b,c
          Character (len=4) :: prefix1, prefix2
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

          Phase_bin = Obs_scal(Nobs_scal)/dble(Nobs)
          
          File_pr ="ener"
          Call Print_scal(Obs_scal, Nobs, file_pr, Group_Comm)

!           do no=1,Norb
!           do no1=no,Norb
!             if(no>1 .or. no1>1) then
!               if (no==no1) then
!                 Local_eq0(1)=Local_eq0(1)+Local_eq0((no-1)*Norb+no1)**2
!               else
!                 Local_eq0(1)=Local_eq0(1)+0.5d0*(Local_eq0((no-1)*Norb+no1)+Local_eq0((no1-1)*Norb+no))**2
!                 Local_eq0(1)=Local_eq0(1)-0.5d0*(Local_eq0((no-1)*Norb+no1)-Local_eq0((no1-1)*Norb+no))**2
!               endif
!             else
!               Local_eq0(1)=Local_eq0(1)**2
!             endif
!           enddo
!           enddo
!           Local_eq0(1)=sqrt(Local_eq0(1))
          File_pr ="Local_eq"
          Call Print_bin(Local_eq, Local_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="Den_eq"
          Call Print_bin(Den_eq, Den_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="U1_eq"
          Call Print_bin(U1_eq, U1_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="U1G_eq"
          Call Print_bin(U1G_eq, U1_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="U1xy_eq"
          Call Print_bin(U1xy_eq, U1xy_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="U1xyG_eq"
          Call Print_bin(U1xyG_eq, U1xyG_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="Spinz_eq"
          Call Print_bin(Spinz_eq, Spinz_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="SpinzG_eq"
          Call Print_bin(SpinzG_eq, Spinz_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          File_pr ="Spinxy_eq"
          Call Print_bin(Spinxy_eq, Spinxy_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          
!           write(*,*) "Hi"
          do a=1,3
            if (a==1) then
              prefix1="T1_"
            elseif (a==2) then
              prefix1="Tz_"
            elseif (a==3) then
              prefix1="Txy_"
            endif
            do b=1,3
              if (b==1) then
                prefix2="S1_"
              elseif (b==2) then
                prefix2="Sz_"
              elseif (b==3) then
                prefix2="Sxy_"
              endif
              write(File_pr,'(A,A,A)') trim(ADJUSTL(prefix1)),trim(ADJUSTL(prefix2)),"1_eq"
              Call Print_bin(Local1_eq(:,:,:,a,b), Local1_eq0(:,a,b), Latt, Nobs, Phase_bin, file_pr, Group_Comm)
              write(File_pr,'(A,A,A)') trim(ADJUSTL(prefix1)),trim(ADJUSTL(prefix2)),"2_eq"
              Call Print_bin(Local2_eq(:,:,:,a,b), Local2_eq0(:,a,b), Latt, Nobs, Phase_bin, file_pr, Group_Comm)
              write(File_pr,'(A,A,A)') trim(ADJUSTL(prefix1)),trim(ADJUSTL(prefix2)),"3_eq"
              Call Print_bin(Local3_eq(:,:,:,a,b), Local3_eq0(:,a,b), Latt, Nobs, Phase_bin, file_pr, Group_Comm)
              write(File_pr,'(A,A,A)') trim(ADJUSTL(prefix1)),trim(ADJUSTL(prefix2)),"4_eq"
              Call Print_bin(Local4_eq(:,:,:,a,b), Local4_eq0(:,a,b), Latt, Nobs, Phase_bin, file_pr, Group_Comm)
              write(File_pr,'(A,A,A)') trim(ADJUSTL(prefix1)),trim(ADJUSTL(prefix2)),"5_eq"
              Call Print_bin(Local5_eq(:,:,:,a,b), Local5_eq0(:,a,b), Latt, Nobs, Phase_bin, file_pr, Group_Comm)
              write(File_pr,'(A,A,A)') trim(ADJUSTL(prefix1)),trim(ADJUSTL(prefix2)),"6_eq"
              Call Print_bin(Local6_eq(:,:,:,a,b), Local6_eq0(:,a,b), Latt, Nobs, Phase_bin, file_pr, Group_Comm)
              write(File_pr,'(A,A,A)') trim(ADJUSTL(prefix1)),trim(ADJUSTL(prefix2)),"7_eq"
              Call Print_bin(Local7_eq(:,:,:,a,b), Local7_eq0(:,a,b), Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            enddo
          enddo
!           write(*,*) "Done"

          if (FlagSym == 1) then
            File_pr ="R_eq"
            Call Print_bin(R_eq, R_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheckTR_eq"
            Call Print_bin(TRS_eq, TRS_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheckR_eq"
            Call Print_bin(RS_eq, RS_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheckC4_eq"
            Call Print_bin(C4S_eq, C4S_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheckPx_eq"
            Call Print_bin(PxS_eq, PxS_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheckSx_eq"
            Call Print_bin(SxS_eq, SxS_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheckSz_eq"
            Call Print_bin(SzS_eq, SzS_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheck1U1_eq"
            Call Print_bin(U11S_eq, U11S_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
            File_pr ="SymCheck2U1_eq"
            Call Print_bin(U12S_eq, U12S_eq0, Latt, Nobs, Phase_bin, file_pr, Group_Comm)
          endif
          
          If (Ltau == 1) then
             Phase_tau = Phase_tau/dble(NobsT)
             File_pr = "Green_tau"
             Call Print_bin_tau(Green_tau,Latt,NobsT,Phase_tau, file_pr,dtau, Group_Comm)
             File_pr = "Den_tau"
             Call Print_bin_tau(Den_tau,Latt,NobsT,Phase_tau, file_pr,dtau, Group_Comm)
             File_pr = "U1_tau"
             Call Print_bin_tau(U1_tau,Latt,NobsT,Phase_tau, file_pr,dtau, Group_Comm)
             File_pr = "U1xy_tau"
             Call Print_bin_tau(U1xy_tau,Latt,NobsT,Phase_tau, file_pr,dtau, Group_Comm)
             File_pr = "U1xyG_tau"
             Call Print_bin_tau(U1xyG_tau,Latt,NobsT,Phase_tau, file_pr,dtau, Group_Comm)
             File_pr = "Spinz_tau"
             Call Print_bin_tau(Spinz_tau,Latt,NobsT,Phase_tau, file_pr,dtau, Group_Comm)
             File_pr = "Spinxy_tau"
             Call Print_bin_tau(Spinxy_tau,Latt,NobsT,Phase_tau, file_pr,dtau, Group_Comm)
             File_pr ="Den_sus"
             Call Print_bin(Den_sus, Den_sus0, Latt, NobsT, Phase_tau, file_pr, Group_Comm)
             File_pr ="U1_sus"
             Call Print_bin(U1_sus, U1_sus0, Latt, NobsT, Phase_tau, file_pr, Group_Comm)
             File_pr ="U1xy_sus"
             Call Print_bin(U1xy_sus, U1xy_sus0, Latt, NobsT, Phase_tau, file_pr, Group_Comm)
             File_pr ="U1xyG_sus"
             Call Print_bin(U1xyG_sus, U1xyG_sus0, Latt, NobsT, Phase_tau, file_pr, Group_Comm)
             File_pr ="Spinz_sus"
             Call Print_bin(Spinz_sus, Spinz_sus0, Latt, NobsT, Phase_tau, file_pr, Group_Comm)
             File_pr ="Spinxy_sus"
             Call Print_bin(Spinxy_sus, Spinxy_sus0, Latt, NobsT, Phase_tau, file_pr, Group_Comm)
             L_eq0=L_eq0*sqrt(beta)
             File_pr ="L_eq"
             Call Print_bin(L_eq, L_eq0, Latt, NobsT, Phase_tau, file_pr, Group_Comm)
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
          
!$OMP parallel do default(shared) private(I,I1,I2,no,J,J1,J2,no1,imj,a,b,c,d,tmp,weight,weightbeta,signum,DeltaI,DeltaJ)
             do I=1,Latt%N
               do a=1,4
               do b=1,4
                 if ( abs(gamma_M(a,b,3))>0.01 ) then
                   do no=1,4
                     I1 = InvList(I,4*(no-1)+a)
                     J1 = InvList(I,4*(no-1)+b)
                     DeltaI=0.d0
                     if (I1==J1) DeltaI=cmplx(1.d0,0.d0,kind(0.D0))
!$OMP CRITICAL
                     L_eq0(1)  = L_eq0(1)  + Gamma_M(a,b,3)*(DeltaI - G00(j1,i1,1)) * ZP* ZS
!$OMP END CRITICAL
                   enddo
                 endif
               enddo
               enddo
             enddo
!$OMP end parallel do
          endif
          
          weightbeta = cmplx(dtau,0.d0,kind(0.D0))
          If (NT == 0 .or. NT==Ltrot) weightbeta=0.5*weightbeta 
          
          If ( N_FL == 1 ) then 
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
!$OMP parallel do default(shared) private(I,I1,I2,no,J,J1,J2,no1,imj,a,b,c,d,tmp,weight,weightbeta,signum,DeltaI,DeltaJ)
             Do I1 = 1,Ndim
                I  = List(I1,1)
                no = List(I1,2)
                Do J1 = 1,Ndim
                   J  = List(J1,1)
                   no1 = List(J1,2)
                   imj = latt%imj(I,J)
!$OMP CRITICAL
                   Green_tau(imj,nt+1,no,no1) = green_tau(imj,nt+1,no,no1)  +  Z * GT0(I1,J1,1) * ZP* ZS
!$OMP END CRITICAL
                   
                   I2=I1
                   J2=J1
                   DeltaI=0.d0
                   if (I1==I2) DeltaI=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                   DeltaJ=0.d0
                   if (J1==J2) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                   tmp =  Z * ( (DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS  - 0.25d0
                   
!$OMP CRITICAL
                   Den_tau  (imj,nt+1,1,1) = Den_tau  (imj,nt+1,1,1)  +  tmp
                   Den_sus  (imj,1,1) = Den_sus  (imj,1,1)  +  weightbeta*tmp
!$OMP END CRITICAL
                     
                   weight=cmplx(1.d0,0.d0,kind(0.D0))
                   if ( (no>=9 .and. no1<=8) .or. (no<=8 .and. no1>=9) ) weight=-weight
!$OMP CRITICAL
                   Spinz_tau (imj,nt+1,1,1) = Spinz_tau (imj,nt+1,1,1)   +   weight * 0.25 * tmp
                   Spinz_sus (imj,1,1) = Spinz_sus (imj,1,1)   +  weightbeta * weight * 0.25 * tmp
!$OMP END CRITICAL
                    
                   signum = 1
                   if (((no-1)/4+1==2) .or. ((no-1)/4+1==4)) signum=-1
                   if (((no1-1)/4+1==2) .or. ((no1-1)/4+1 ==4)) signum=-signum
                   weight = cmplx(dble(signum),0.d0, kind(0.D0))
!$OMP CRITICAL
                   U1_tau (imj,nt+1,1,1) = U1_tau (imj,nt+1,1,1)   +  weight*tmp*0.25
                   U1_sus (imj,1,1) = U1_sus (imj,1,1)   + weightbeta* weight*tmp*0.25
!$OMP END CRITICAL

                enddo
              enddo
!$OMP end parallel do
              
!$OMP parallel do default(shared) private(I,I1,I2,no,J,J1,J2,no1,imj,a,b,c,d,tmp,weight,weightbeta,signum,DeltaI,DeltaJ)
              do I=1,Latt%N
                do no=1,8
                    do J=1,Latt%N
                      imj = latt%imj(I,J)
                      do no1=1,8
                          I1 = Invlist(I,no)
                          I2 = Invlist(I,no+8)
                          J1 = Invlist(J,no1+8)
                          J2 = Invlist(J,no1)
                          
                          DeltaI=0.d0
                          if (I1==I2) DeltaI=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                          DeltaJ=0.d0
                          if (J1==J2) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                          tmp =  Z * ((DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS
!$OMP CRITICAL
                          Spinxy_tau (imj,nt+1,1,1) = Spinxy_tau (imj,nt+1,1,1)   +  tmp
                          Spinxy_sus (imj,1,1) = Spinxy_sus (imj,1,1)   +   weightbeta*tmp
!$OMP END CRITICAL
                          
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
                          
                          DeltaI=0.d0
                          if (I1==I2) DeltaI=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                          DeltaJ=0.d0
                          if (J1==J2) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                          tmp =  Z * ((DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS
!$OMP CRITICAL
                          U1xy_tau (imj,nt+1,1,1) = U1xy_tau (imj,nt+1,1,1)   +  tmp
                          U1xy_sus (imj,1,1) = U1xy_sus (imj,1,1)   +  weightbeta*tmp
                          U1xyG_tau (imj,nt+1,1,1) = U1xyG_tau (imj,nt+1,1,1)   +  (-1)**(no/2+no1/2+(no-1)/4+(no1-1)/4)*tmp
                          U1xyG_sus (imj,1,1) = U1xyG_sus (imj,1,1)   +  (-1)**(no/2+no1/2+(no-1)/4+(no1-1)/4)*weightbeta*tmp
!$OMP END CRITICAL
                      enddo
                    enddo
                    
                Enddo
             Enddo
!$OMP end parallel do
             
             
            
!$OMP parallel do default(shared) private(I,I1,I2,no,J,J1,J2,no1,imj,a,b,c,d,tmp,weight,weightbeta,signum,DeltaI,DeltaJ)
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
                            DeltaI=0.d0
                            if (I1==I2) DeltaI=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                            DeltaJ=0.d0
                            if (J1==J2) DeltaJ=cmplx(1.d0,0.d0, kind(0.D0)) ! .and. nt==0
                            tmp =  Z * ((DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS
!$OMP CRITICAL
                            L_eq (imj,1,1) = L_eq (imj,1,1)   +  weight*tmp*weightbeta 
!$OMP END CRITICAL
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
!$OMP end parallel do
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
          Integer :: v, t, M_v, L_trot, kind_m
          
!           nsigma = -nsigma_old

          kind_m = nranf(2)
          if (kind_m == 1) then
!             write(*,*) 'Flipping x-component'
            nsigma(1:M_v/2,:)     = -nsigma_old(1:M_v/2,:)
            nsigma(M_v/2+1:M_v,:) =  nsigma_old(M_v/2+1:M_v,:)
          else
!             write(*,*) 'Flipping y-component'
            nsigma(1:M_v/2,:)     =  nsigma_old(1:M_v/2,:)
            nsigma(M_v/2+1:M_v,:) = -nsigma_old(M_v/2+1:M_v,:)
          endif
          
!           do t=1,L_trot
!             do v=1,M_v
! !               nsigma(v,t) = -nsigma_old(v,t)
!               if (abs(aimag(Op_V(v,1)%g))>0.00001d0 ) then
!                 nsigma(v,t) =  nsigma_old(v,t)
!               else
!                 nsigma(v,t) = -nsigma_old(v,t)
!               endif
!             enddo
!           enddo
            
            
!           L_trot = size(nsigma,2)
!           M_v    = size(nsigma,1) 
!           nsigma(1:M_v/2,:) = -nsigma_old(M_v/2+1:M_v,:)
!           nsigma(M_v/2+1:M_v,:) =  nsigma_old(1:M_v/2,:)
          
          
!           kind_m = nranf(3)
!           
!           if (kind_m == 1 .or. kind_m == 3) then
!             write(*,*) 'Flipping x-component'
!             nsigma(1:M_v/2,:) = -nsigma_old(1:M_v/2,:)
!           else
!             nsigma(1:M_v/2,:) =  nsigma_old(1:M_v/2,:)
!           endif
!           if (kind_m == 2 .or. kind_m == 3) then
!             write(*,*) 'Flipping y-component'
!             nsigma(M_v/2+1:M_v,:) = -nsigma_old(M_v/2+1:M_v,:)
!           else
!             nsigma(M_v/2+1:M_v,:) =  nsigma_old(M_v/2+1:M_v,:)
!           endif
          
!           do t=1,L_trot
!             do v=1,M_v
! !               nsigma(v,t) = -nsigma_old(v,t)
!               if (v <= M_v/2 .and. kind_m .ne. 2 ) then
!                 nsigma(v,t) = -nsigma_old(v,t)
!               else
!                 nsigma(v,t) = nsigma_old(v,t)
!               endif
!               if (v > M_v/2 .and. kind_m .ne. 1 ) then
!                 nsigma(v,t) = -nsigma_old(v,t)
!               else
!                 nsigma(v,t) = nsigma_old(v,t)
!               endif
!             enddo
!           enddo
          
          T0_Proposal_ratio=1.d0
          
        End Subroutine Global_move
!========================================================================

!---------------------------------------------------------------------
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
          
          Delta_S0_global=1.0
          
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
