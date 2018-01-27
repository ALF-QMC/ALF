    !Model Hamiltonian for interaction-induced topological reduction
    Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      

      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Type (WaveFunction), dimension(:),   allocatable  :: WF_L
      Type (WaveFunction), dimension(:),   allocatable  :: WF_R
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  Ltrot, Thtrot
      Integer, parameter   :: N_FL=1,  N_SUN=2
      Logical              :: Projector
!>    Defines MPI communicator 
      Integer              :: Group_Comm


      
      ! What is below is  private 
      
      Type (Lattice),         private :: Latt
      Integer,                private :: Norb=4
      Integer, allocatable,   private :: List(:,:), Invlist(:,:)
      Integer,                private :: L1, L2
      real (Kind=Kind(0.d0)), private :: Ham_T, Ham_T2, Ham_Vint,  Ham_U, Ham_J
      real (Kind=Kind(0.d0)), private :: Dtau, Beta, Theta
      Character (len=64),     private :: Model, Lattice_type, File1
      logical,                private :: checkerboard


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau

      contains 

        Subroutine Ham_Set

#if defined(MPI)
          USE mpi
#endif   
          implicit none

          integer :: ierr

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_Z2_FLstar/  ham_T, ham_T2, Ham_Vint,  Ham_U,  Dtau, Beta, Ham_J, checkerboard, Theta, Projector

#if defined(MPI)
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
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
             CLOSE(5)
#if defined(MPI)
          endif

          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif

          Projector = .false.
          Theta = 0.d0
          Thtrot = 0
          Ham_J=0.d0
          Ham_T2=-1000000.d0
          checkerboard=.false.
          
#if defined(MPI) 
          If (Irank_g == 0 ) then
#endif
             File1 = "parameters"
#if defined(TEMPERING) 
             write(File1,'(A,I0,A)') "Temp_",igroup,"/parameters"
#endif
             OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_Z2_FLstar)
             CLOSE(5)
             if (Ham_T2==-1000000.d0) Ham_T2=Ham_T
#if defined(MPI)
          endif

          CALL MPI_BCAST(ham_T        ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_T2       ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_Vint     ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U        ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_J        ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau         ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta         ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(Theta        ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(Thtrot       ,1,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Projector    ,1,MPI_LOGICAL,   0,Group_Comm,ierr)
          CALL MPI_BCAST(checkerboard ,1,MPI_LOGICAL,   0,Group_Comm,ierr)
#endif
          if(ham_j==0.d0) Norb=2

          Call Ham_latt

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
             Write(50,*) 'L1            : ', L1
             Write(50,*) 'L2            : ', L2
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 't`            : ', Ham_T2
             Write(50,*) 'V             : ', Ham_Vint
             Write(50,*) 'U             : ', Ham_U
             Write(50,*) 'JK            : ', Ham_J
             If(checkerboard) Write(50,*) 'Using checkerboard decomposition for Ham_T'
             close(50)
#if defined(MPI)
          endif
#endif
          call Ham_V
        end Subroutine Ham_Set
!=============================================================================

        Subroutine Ham_Latt
          Implicit none
          !Set the lattice
          Integer :: no, I, nc, ndimloc, Ihex(6), I1, I2
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer, allocatable :: Counter(:)
          
          
          If  ( Lattice_type == "Kagome" ) then
             a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
             a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0

             !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
             L1_p    =  dble(L1) * a1_p
             L2_p    =  dble(L2) * a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
!              do I=1,L1*L2
!                 write(*,*) I, Latt%list(I,1), Latt%list(I,2), Latt%list(I,1)*a1_p+ Latt%list(I,2)*a2_p
!              enddo
!              write(*,*)
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif
          
!           Do I = -1,1
!              Do no = -1,1
!                 write(*,*) I,no,Latt%nnlist(1,I,no)
!              Enddo
!           Enddo

          Ndim = Latt%N*Norb
          Allocate (List(Ndim,Norb), Invlist(Latt%N,Norb))
          
          nc = 0
          ! "compact" f-electron Lattice supporting the Z2 spinliquid phase
          Do I = 1,Latt%N
             Do no = 1,2
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
             Enddo
          Enddo
          ! compact c-electrons
          if (norb==4) then
            Do I = 1,Latt%N
              Do no = 3,4
                  nc = nc + 1
                  List(nc,1) = I
                  List(nc,2) = no
                  Invlist(I,no) = nc 
              Enddo
            Enddo
          endif
          
          call Print_latt(Latt)

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !  Setup the hopping
          !  Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: n, Ncheck, nc, ncoord, I, I1, I2,o
          Complex(kind=kind(0.d0)) :: T
          
          If (checkerboard) then
            Ncheck = 3*(Norb/2)*Latt%N
            allocate(Op_T(Ncheck,N_FL))
            Do o=0,(Norb/2)-1
              T=cmplx(ham_t,0.d0,kind(0.d0))
              if (o==1) T=cmplx(ham_t2,0.d0,kind(0.d0))
              do n = 1,N_FL
                nc = 2*Latt%N*o
                Do ncoord = 1,3
                  Do I = 1,Latt%N
                    nc = nc + 1
                    I1 = invlist(I,1+2*o)
                    if      ( ncoord == 1 ) then 
                        I2 = invlist(I,2+2*o)
                    elseif  ( ncoord == 2 ) then
                        I2 = Invlist( Latt%nnlist(I,1,-1),2+2*o )
                    elseif  ( ncoord == 3 ) then
                        I2 = invlist( Latt%nnlist(I,0,-1),2+2*o )
                    endif
                    Call Op_make(Op_T(nc,n),2)
                    Op_T(nc,n)%P(1) = I1
                    Op_T(nc,n)%P(2) = I2
                    Op_T(nc,n)%O( 1 , 2 ) = -T
                    Op_T(nc,n)%O( 2 , 1 ) = -T
                    Op_T(nc,n)%g=cmplx(-Dtau,0.d0,Kind(0.d0))
                    !Write(6,*) 'In Ham_hop', Ham_T
                    Call Op_set(Op_T(nc,n)) 
                    !Write(6,*) 'In Ham_hop 1'
                    !Do I = 1,2*Latt%N
                    !   Write(6,*) Op_T(nc,n)%E(i)
                    !enddo
                  enddo
                enddo
              enddo
            enddo
          else
            Ncheck=(Norb/2)
            allocate(Op_T(Ncheck,N_FL))
            do n = 1,N_FL
              Call Op_make(Op_T(1,n),2*Latt%N)
              if(norb==4) Call Op_make(Op_T(2,n),2*Latt%N)
              Do I = 1,2*Latt%N
                Op_T(1,n)%P(I)=I
                if(norb==4) Op_T(2,n)%P(I)=2*Latt%N+I
              Enddo
              Do ncoord = 1,3
                Do I = 1,Latt%N
                  I1 = invlist(I,1)
                  if      ( ncoord == 1 ) then 
                      I2 = invlist(I,2)
                  elseif  ( ncoord == 2 ) then
                      I2 = Invlist( Latt%nnlist(I,1,-1),2 )
                  elseif  ( ncoord == 3 ) then
                      I2 = invlist( Latt%nnlist(I,0,-1),2 )
                  endif
                  Op_T(1,n)%O( I1 , I2 ) = -cmplx(ham_t,0.d0,kind(0.d0))
                  Op_T(1,n)%O( I2 , I1 ) = -cmplx(ham_t,0.d0,kind(0.d0))
                  if(norb==4) Op_T(2,n)%O( I1 , I2 ) = -cmplx(ham_t2,0.d0,kind(0.d0))
                  if(norb==4) Op_T(2,n)%O( I2 , I1 ) = -cmplx(ham_t2,0.d0,kind(0.d0))
                enddo
              enddo
              Op_T(1,n)%g=cmplx(-Dtau,0.d0,Kind(0.d0))
              Call Op_set(Op_T(1,n)) 
              if(norb==4) Op_T(2,n)%g=cmplx(-Dtau,0.d0,Kind(0.d0))
              if(norb==4) Call Op_set(Op_T(2,n)) 
            enddo
          endif

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
          
          Allocate(H0(2*Latt%N,2*Latt%N),U(2*Latt%N,2*Latt%N),En(2*Latt%N))
          H0=0.d0
          DO I = 1, Latt%N
            I1 = Invlist(I,1)
            J1 = I1
            Do nc1 = 1,3
              select case (nc1)
              case (1)
                J1 = invlist(I,2)
              case (2)
                J1 = invlist(Latt%nnlist(I,1,-1),2)
              case (3)
                J1 = invlist(Latt%nnlist(I,0,-1),2)
              case default
                Write(6,*) ' Error in  Ham_Hop '  
                Stop
              end select
              If ( Latt%list(J1,1) == 0 .and. nc1==2 ) then
                H0(I1,J1) = cmplx( -cos(phi),    -sin(phi), kind(0.D0))
                H0(J1,I1) = cmplx( -cos(phi),     sin(phi), kind(0.D0))
              else
                H0(I1,J1) = cmplx(-1.d0,    0.d0, kind(0.D0))
                H0(J1,I1) = cmplx(-1.d0,    0.d0, kind(0.D0))
              endif
            Enddo
          enddo
          
          Do I1=1,L1,3
            do J1=1,L2,3
              I=Latt%invlist(I1,J1)
              do nc=1,3
                I2=invlist(I,1)
                J2=invlist(I,2)
                H0(I2,J2) = H0(I2,J2) + kek
                H0(J2,I2) = H0(J2,I2) + kek
                
                I2=invlist(Latt%nnlist(I,1,-1),2)
                J2=invlist(Latt%nnlist(I,1,0),1)
                H0(I2,J2) = H0(I2,J2) + kek
                H0(J2,I2) = H0(J2,I2) + kek
                
                I2=invlist(Latt%nnlist(I,1,0),2)
                J2=invlist(Latt%nnlist(I,0,1),1)
                H0(I2,J2) = H0(I2,J2) + kek
                H0(J2,I2) = H0(J2,I2) + kek
                
                I=Latt%nnlist(I,1,1)
              enddo
            enddo
          enddo
          
!           do I2=1,size(H0,1)
!             write(*,*) nint(dble(H0(I2,:)))
!           enddo
          
          Call Diag(H0,U,En)
          
          write(*,*) N_part, Ndim/2, (norb/2)*Latt%N
          write(*,*) 'Gap is', En(Latt%N+1)-En(Latt%N)
          do I2=1,Latt%N
          do I1=1,2*Latt%N
            WF_L(1)%P(I1,I2)=U(I1,I2)
            WF_R(1)%P(I1,I2)=U(I1,I2)
            if(norb==4) WF_L(1)%P(I1+2*Latt%N,I2+Latt%N)=U(I1,I2)
            if(norb==4) WF_R(1)%P(I1+2*Latt%N,I2+Latt%N)=U(I1,I2)
          enddo
          enddo
          
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
          
          Integer :: nf, I, I1, I2, nc, no, ncoord
          Integer :: Ibond(2,6), Ihex(6)

          Complex (Kind=Kind(0.d0)) :: Zone


          ! Number of opertors 8 per unit cell
          write(*,*) Ham_Vint, Ham_Vint/dble(2*N_SUN)
          Allocate( Op_V((3+(norb/4)*2)*Latt%N,N_FL) )
          nc = 0
          nf = 1
          Zone=cmplx(1.d0  ,0.d0, kind(0.D0))
          Do I = 1,Latt%N
            Do ncoord = 1,3
              nc = nc + 1
              I1 = invlist(I,1)
              if      ( ncoord == 1 ) then 
                  I2 = invlist(I,2)
              elseif  ( ncoord == 2 ) then
                  I2 = Invlist( Latt%nnlist(I,1,-1),2 )
              elseif  ( ncoord == 3 ) then
                  I2 = invlist( Latt%nnlist(I,0,-1),2 )
              endif
              Call Op_make(Op_V(nc,nf),2)
              Op_V(nc,nf)%P(1) = I1
              Op_V(nc,nf)%P(2) = I2
              Op_V(nc,nf)%O( 1 , 2 ) = Zone
              Op_V(nc,nf)%O( 2 , 1 ) = Zone
              Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_Vint/dble(2*N_SUN),0.d0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
            enddo
            
            ! auxiliary interaction to project out spin fluctuations reducing spinful fermions to effective hardcore bosons
!             do no=1,2
!               nc = nc+1
! !                 write(*,*) Invlist(I,no)
!               Call Op_make(Op_V(nc,nf),1) 
!               Op_V(nc,nf)%P(1) = Invlist(I,no)
!               Op_V(nc,nf)%O( 1 , 1 )  =  Zone
!               Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_U/(beta+2*theta),0.d0, kind(0.D0)))
!               Op_V(nc,nf)%alpha  = cmplx(-0.5d0, 0.d0, kind(0.D0))
!               Op_V(nc,nf)%type   = 2
!               Call Op_set( Op_V(nc,nf) )
!             enddo
            
            ! Kondo coupling of f- and c-electrons
            if(norb==4) then
              do no=1,2
                nc = nc+1
                Call Op_make(Op_V(nc,nf),2) 
                Op_V(nc,nf)%P(1) = Invlist(I,no)
                Op_V(nc,nf)%P(2) = Invlist(I,2+no)
                Op_V(nc,nf)%O( 1 , 2 )  =  Zone
                Op_V(nc,nf)%O( 2 , 1 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_J*0.25d0,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(0.0d0, 0.d0, kind(0.D0))
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
              enddo
            endif
            
          enddo
          
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
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename

          ! Scalar observables
          Allocate ( Obs_scal(12) )
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
             case (5)
                N = 1;   Filename ="DoubleOcc"
             case (6)
                N = 1;   Filename ="SC"
             case (7)
                N = 1;   Filename ="EnerFermion"
             case (8)
                N = 1;   Filename ="EnerBoson"
             case (9)
                N = 1;   Filename ="ConductionKin"
             case (10)
                N = 1;   Filename ="ConductionPart"
             case (11)
                N = 1;   Filename ="EnerKondo"
             case (12)
                N = 1;   Filename ="ConductionDoubleOcc"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(3+Norb/4) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Ns = Latt%N;  No = Norb;  Filename ="Green"
             case (8)
                Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
             case (6)
                Ns = Latt%N;  No = Norb;  Filename ="SpinXY"
             case (2)
                Ns = Latt%N;  No = Norb;  Filename ="Den"
             case (5)
                Ns = Latt%N;  No = Norb;  Filename ="SC"
             case (3)
                Ns = Latt%N;  No = 3*(Norb/2);  Filename ="Hop"
             case (4)
                Ns = Latt%N;  No = 2;  Filename ="Hyb"
             case (7)
                Ns = Latt%N;  No = Norb;  Filename ="SuperSpin"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau > 0) then
             Allocate(Obs_tau(2))
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Ns = Latt%N; No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N; No = Norb;  Filename ="Den"
                case (3)
                   Ns = Latt%N; No = Norb;  Filename ="SC"
                case (4)
                   Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                case (5)
                   Ns = Latt%N; No = Norb;  Filename ="SpinXY"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1-2*Thtrot
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
          endif
          
        end Subroutine Alloc_obs
!===================================================================================           
        
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
!========================================================================

        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, Ztot, Phii, Phij
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZkinF, ZPot, Zocc, Z, ZP,ZS, weight, tmp,Zn
          Integer :: I,J, n, imj, nf, I1, J1, Nc, I0, I2, I3
          Integer :: K, K1, L ,L1, no_I, no_J, Ibond(2,6), Ihex(6)
          real(kind=kind(0.d0)), parameter :: pi = 4 * atan(1.d0)
          
          Zn=cmplx( dble(N_SUN), 0.d0 , kind(0.D0))

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
          
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ! Compute scalar observables. 
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
             

          Zkin = cmplx(0.d0,0.d0, kind(0.D0))
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
          
          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Latt%N
                ! interaction on the bonds representing hopping of hard core bosons
                I0 = invlist(I,1)
                I1 = invlist(I,2)
                Zpot = ZPot + Corr(Gr,Grc,I0,I1,I0,I1) + Corr(Gr,Grc,I0,I1,I1,I0) &
                       &  + Corr(Gr,Grc,I1,I0,I0,I1) + Corr(Gr,Grc,I1,I0,I1,I0)
                I1 = Invlist( Latt%nnlist(I,1,-1),2 )
                Zpot = ZPot + Corr(Gr,Grc,I0,I1,I0,I1) + Corr(Gr,Grc,I0,I1,I1,I0) &
                       &  + Corr(Gr,Grc,I1,I0,I0,I1) + Corr(Gr,Grc,I1,I0,I1,I0)
                I1 = invlist( Latt%nnlist(I,0,-1),2 ) 
                Zpot = ZPot + Corr(Gr,Grc,I0,I1,I0,I1) + Corr(Gr,Grc,I0,I1,I1,I0) &
                       &  + Corr(Gr,Grc,I1,I0,I0,I1) + Corr(Gr,Grc,I1,I0,I1,I0)
             Enddo
          Enddo
          ZPot = -Ham_Vint/dble(2*N_SUN)*ZPot
          
          ! Total energy (the kinetic part of con electrons is missing
          Nc = Size( Op_V,1)
          ZTot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do n = 1,Nc
                weight=-Op_V(n,nf)%g**2.d0 /dtau !(-1)**((n-1)/Latt%N)/8.d0
!                 write(0,*) Op_V(n,nf)%g, weight
                Do J = 1,Op_V(n,nf)%N
                   J1 = Op_V(n,nf)%P(J)
                   DO I = 1,Op_V(n,nf)%N
                      if (abs(Op_V(n,nf)%O(i,j)) >= 0.00001) then
                      I1 = Op_V(n,nf)%P(I)
                      ZTot  = ZTot  + 2.d0*Zn*Op_V(n,nf)%alpha*weight*Op_V(n,nf)%O(i,j)*GRC(I1,J1,1)
                      Do K = 1,Op_V(n,nf)%N
                        K1 = Op_V(n,nf)%P(K)
                        DO L = 1,Op_V(n,nf)%N
                          if (abs(Op_V(n,nf)%O(k,l)) >= 0.00001) then
                          L1 = Op_V(n,nf)%P(L)
                          tmp =  (   GRC(I1,L1,1) * GR (J1,K1,1)      +  &
                                &    Zn * GRC(I1,J1,1) * GRC(K1,L1,1)         )
                          ZTot  = ZTot  + weight*Op_V(n,nf)%O(i,j)*Op_V(n,nf)%O(k,l)*tmp
                          endif
                        Enddo
                      ENddo
                      endif
                   Enddo
                ENddo
                ZTot  = ZTot  + weight*(Op_V(n,nf)%alpha**2)*Zn
!           write(*,*) Zpot
             Enddo
          Enddo
          ZTot = ZTot*Zn
          Obs_scal(1)%Obs_vec(1)  =  Obs_scal(1)%Obs_vec(1) + Zkin * ZP*ZS
          Obs_scal(8)%Obs_vec(1)  =  Obs_scal(8)%Obs_vec(1) + (Zkin + Zpot) * ZP*ZS
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =  Obs_scal(4)%Obs_vec(1) + ZTot*ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Zocc = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,2*Latt%N
                Zrho = Zrho + Grc(i,i,nf) 
                Zocc = Zocc + Grc(i,i,nf)**2
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(5)%Obs_vec(1)  =    Obs_scal(5)%Obs_vec(1) - 2.d0*(Zocc - 0.5d0*Zrho) * ZP*ZS
          
          ! conduction fermions
          if (norb==4) then
            Zrho = cmplx(0.d0,0.d0, kind(0.D0))
            Zocc = cmplx(0.d0,0.d0, kind(0.D0))
            Do nf = 1,N_FL
              Do I = 2*Latt%N+1,4*Latt%N
                  Zrho = Zrho + Grc(i,i,nf) 
                  Zocc = Zocc + Grc(i,i,nf)**2
              enddo
            enddo
            Zrho = Zrho* dble(N_SUN)
            Obs_scal(10)%Obs_vec(1)  =    Obs_scal(10)%Obs_vec(1) + Zrho * ZP*ZS
            Obs_scal(12)%Obs_vec(1)  =    Obs_scal(12)%Obs_vec(1) - 2.d0*(Zocc - 0.5d0*Zrho) * ZP*ZS
            
            Zkin = cmplx(0.d0,0.d0,Kind(0.d0))
            Do nf = 1,N_FL
              Do I = 1,Latt%N
                  I0 = invlist(I,3)
                  I1 = invlist(I,4)
                  I2 = Invlist( Latt%nnlist(I,1,-1),4 )
                  I3 = invlist( Latt%nnlist(I,0,-1),4 )  
                  Zkin = Zkin + Grc(I0,I1,nf) +  Grc(I1,I0,nf)  + &
                      &        Grc(I0,I2,nf) +  Grc(I2,I0,nf)  + &
                      &        Grc(I0,I3,nf) +  Grc(I3,I0,nf)
              Enddo
            Enddo
            Zkin = - Zkin*cmplx( Ham_T2*dble(N_SUN), 0.d0,Kind(0.d0) )
            Obs_scal(9)%Obs_vec(1)  =    Obs_scal(9)%Obs_vec(1) + Zkin * ZP*ZS
            !add kin part of conduction electrons to total energy
            Obs_scal(4)%Obs_vec(1)  =  Obs_scal(4)%Obs_vec(1) + Zkin*ZP*ZS
            
            ZPot = cmplx(0.d0,0.d0,Kind(0.d0))
            Do nf = 1,N_FL
              Do I = 1,Latt%N
                Do J = 1,2
                  I0 = invlist(I,2+J)
                  I1 = invlist(I,J)  
                  Zpot = ZPot + Corr(Gr,Grc,I0,I1,I0,I1) + Corr(Gr,Grc,I0,I1,I1,I0) &
                        &  + Corr(Gr,Grc,I1,I0,I0,I1) + Corr(Gr,Grc,I1,I0,I1,I0)
                enddo
              Enddo
            Enddo
            Zpot = - ZPot*cmplx( Ham_J*0.25d0, 0.d0,Kind(0.d0) )
            Obs_scal(11)%Obs_vec(1)  =    Obs_scal(11)%Obs_vec(1) + ZPot * ZP*ZS
          endif
          
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          
          Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
          Do J1 = 1,Ndim
            J    = List(J1,1)
            no_J = List(J1,2)
            Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                imj = latt%imj(I,J)
                ! Green
                Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               Z * GRC(I1,J1,1) *  ZP*ZS 
                ! Den
                tmp=Corr(GR,GRC,I1,I1,J1,J1)* ZP*ZS
                Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + tmp
!                 Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) + tmp
!                 ! SpinZ
!                 tmp=Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
!                 Obs_eq(8)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(8)%Obs_Latt(imj,1,no_I,no_J) + tmp
!                 Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) + tmp
!                 ! SpinXY
!                 Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) + tmp
!                 Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) + tmp
!                 ! SC
!                 tmp=2.d0*GRC(I1,J1,1)**2 *  ZP*ZS
!                 Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + tmp
!                 Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) + tmp
                    
                tmp=2.d0*GRC(I1,J1,1)**2 *  ZP*ZS
                Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1) + GRC(I1,J1,1)**2/ 3.d0 /dble(Latt%N) * ZP*ZS
            ENDDO
            Obs_eq(2)%Obs_Latt0(no_J) =  Obs_eq(2)%Obs_Latt0(no_J) +  Z * cmplx(0.5d0,0.d0,kind(0.d0)) * ZP * ZS
!             Obs_eq(7)%Obs_Latt0(no_J) =  Obs_eq(7)%Obs_Latt0(no_J) +  Z * 0.5d0 * ZP * ZS
          ENDDO
          Do J = 1,Latt%N
            Do no_J=1,3*(norb/2)
              select case (no_J)
                case (1)
                  K1=invlist(J,1)
                  L1=invlist(J,2)
                  phij=cmplx(1.d0,0.d0, kind(0.d0))
                case (2)
                  K1=invlist(J,1)
                  L1=Invlist( Latt%nnlist(J,1,-1),2 )
                  phij=cmplx(cos(2.d0*Pi/3.d0),sin(2.d0*Pi/3.d0), kind(0.d0))
                case (3)
                  K1=invlist(J,1)
                  L1=Invlist( Latt%nnlist(J,0,-1),2 )  
                  phij=cmplx(cos(4.d0*Pi/3.d0),sin(4.d0*Pi/3.d0), kind(0.d0))
                case (4)
                  K1=invlist(J,3)
                  L1=invlist(J,4)
                  phij=cmplx(1.d0,0.d0, kind(0.d0))
                case (5)
                  K1=invlist(J,3)
                  L1=Invlist( Latt%nnlist(J,1,-1),4 )
                  phij=cmplx(cos(2.d0*Pi/3.d0),sin(2.d0*Pi/3.d0), kind(0.d0))
                case (6)
                  K1=invlist(J,3)
                  L1=Invlist( Latt%nnlist(J,0,-1),4 )  
                  phij=cmplx(cos(4.d0*Pi/3.d0),sin(4.d0*Pi/3.d0), kind(0.d0))
              end select
                  phij=cmplx(1.d0,0.d0, kind(0.d0))
              Do I = 1,Latt%N
                imj = latt%imj(I,J)
                Do no_i=1,3*(norb/2)
                  select case (no_I)
                    case (1)
                      I1=invlist(I,1)
                      J1=invlist(I,2)
                      phii=cmplx(1.d0,0.d0, kind(0.d0))
                    case (2)
                      I1=invlist(I,1)
                      J1=Invlist( Latt%nnlist(I,1,-1),2 )
                      phii=cmplx(cos(2.d0*Pi/3.d0),sin(2.d0*Pi/3.d0), kind(0.d0))
                    case (3)
                      I1=invlist(I,1)
                      J1=Invlist( Latt%nnlist(I,0,-1),2 )  
                      phii=cmplx(cos(4.d0*Pi/3.d0),sin(4.d0*Pi/3.d0), kind(0.d0))
                    case (4)
                      I1=invlist(I,3)
                      J1=invlist(I,4)
                      phii=cmplx(1.d0,0.d0, kind(0.d0))
                    case (5)
                      I1=invlist(I,3)
                      J1=Invlist( Latt%nnlist(I,1,-1),4 )
                      phii=cmplx(cos(2.d0*Pi/3.d0),sin(2.d0*Pi/3.d0), kind(0.d0))
                    case (6)
                      I1=invlist(I,3)
                      J1=Invlist( Latt%nnlist(I,0,-1),4 )  
                      phii=cmplx(cos(4.d0*Pi/3.d0),sin(4.d0*Pi/3.d0), kind(0.d0))
                  end select
                      phii=cmplx(1.d0,0.d0, kind(0.d0))
                  ! Hop
                  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J)  +  ZP*ZS* &
                      & (Phii*Phij       *Corr(GR,GRC,I1,J1,K1,L1)&
                      & +Phii*conjg(Phij)*Corr(GR,GRC,I1,J1,L1,K1)&
                      & +conjg(Phii)*Phij*Corr(GR,GRC,J1,I1,K1,L1)&
                      & +conjg(Phii*Phij)*Corr(GR,GRC,J1,I1,L1,K1))
                enddo
              ENDDO
              Obs_eq(3)%Obs_Latt0(no_J) = Obs_eq(3)%Obs_Latt0(no_J)+Z*(Phij*GRC(k1,L1,1)+conjg(phij)*GRC(L1,K1,1))*ZP*ZS
            enddo
          ENDDO
          if (abs(ham_J) > 0.d0) then
          Do J = 1,Latt%N
            Do no_J=1,2
              K1=invlist(J,no_J)
              L1=invlist(J,no_J+2)
              Do I = 1,Latt%N
                imj = latt%imj(I,J)
                Do no_i=1,2
                  I1=invlist(I,no_i)
                  J1=invlist(I,no_i+2)
                  ! Hyb
                  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  ZP*ZS* &
                      & (Corr(GR,GRC,I1,J1,K1,L1)&
                      & +Corr(GR,GRC,I1,J1,L1,K1)&
                      & +Corr(GR,GRC,J1,I1,K1,L1)&
                      & +Corr(GR,GRC,J1,I1,L1,K1))
                enddo
              ENDDO
              Obs_eq(4)%Obs_Latt0(no_J) = Obs_eq(4)%Obs_Latt0(no_J)+Z*(GRC(k1,L1,1)+GRC(L1,K1,1))*ZP*ZS
            enddo
          ENDDO
          endif

        end Subroutine Obser
!==========================================================   

        Complex (Kind=kind(0.d0)) Function   Corr(GR, GRC, I1, I2, I3, I4)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GRC(Ndim,Ndim,N_FL)
          Integer, INTENT(IN) :: I1,I2,I3,I4
          
          Corr=(N_SUN*Grc(I1,I2,1)*Grc(I3,I4,1)+Grc(I1,I4,1)*Gr(I2,I3,1))*N_SUN

        end function Corr
!==========================================================   

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
!==========================================================        

        Subroutine OBSERT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I1, I, J1, J, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          
          Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
          Do I1 = 1,Ndim
            I    = List(I1,1)
            no_I = List(I1,2)
            Do J1 = 1,Ndim
                J    = List(J1,1)
                no_J = List(J1,2)
                imj = latt%imj(I,J)
                ! Green
                Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & +  Z * GT0(I1,J1,1) * ZP* ZS
                
!                 ! SpinZ
!                 Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                     &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
!                 
!                 ! SpinXY
!                 Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                     &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                
                ! Den
                Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
                    &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
                    &     Z * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                
!                 ! SC
!                 Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                     & + (    G0T(J1,I1,1)**2.d0                                ) * ZP * ZS
            Enddo
            Obs_tau(2)%Obs_Latt0(no_I) = Obs_tau(2)%Obs_Latt0(no_I) + &
                  &         Z*(cmplx(0.5d0,0.d0,kind(0.d0))) * ZP * ZS
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
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
        
          T0_Proposal_ratio=0.d0
          size_clust=0.d0
          
        End Subroutine Global_move
!========================================================================
!---------------------------------------------------------------------
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
          Delta_S0_global=1.d0
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
          
          T0_Proposal_ratio=1.d0
          Flip_length=0
          S0_ratio=1.d0

        end Subroutine Global_move_tau


    end Module Hamiltonian
