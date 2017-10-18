    !Model Hamiltonian for interaction-induced topological reduction
    Module Hamiltonian

      Use Operator_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      

      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot
!>    Defines MPI communicator 
      Integer              :: Group_Comm


      
      ! What is below is  private 
      
      Type (Lattice),         private :: Latt
      Integer, parameter,     private :: Norb=5
      Integer, allocatable,   private :: List(:,:), Invlist(:,:)
      Integer,                private :: L1, L2
      real (Kind=Kind(0.d0)), private :: Ham_T, Ham_Vint,  Ham_U, Ham_JKA, Ham_JKB
      real (Kind=Kind(0.d0)), private :: Dtau, Beta
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

          NAMELIST /VAR_Z2_FLstar/  ham_T, Ham_Vint,  Ham_U,  Dtau, Beta, Ham_JKA, Ham_JKB, checkerboard


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
          Call Ham_latt


          N_FL  = 1
          N_SUN = 2
          
#if defined(MPI) 
          If (Irank_g == 0 ) then
#endif
             File1 = "parameters"
#if defined(TEMPERING) 
             write(File1,'(A,I0,A)') "Temp_",igroup,"/parameters"
#endif
             Ham_JKA=0.d0
             Ham_JKB=0.d0
             checkerboard=.true.
             OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_Z2_FLstar)
             CLOSE(5)
#if defined(MPI)
          endif

          CALL MPI_BCAST(ham_T        ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_Vint     ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U        ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_JKA      ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_JKB      ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau         ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta         ,1,MPI_REAL8,     0,Group_Comm,ierr)
          CALL MPI_BCAST(checkerboard ,1,MPI_LOGICAL,   0,Group_Comm,ierr)
#endif

          Call Ham_hop
          Ltrot = nint(beta/dtau)
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
             Write(50,*) 'V             : ', Ham_Vint
             Write(50,*) 'U             : ', Ham_U
             Write(50,*) 'JK A sublatt  : ', Ham_JKA
             Write(50,*) 'JK B sublatt  : ', Ham_JKB
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
          Integer :: no, I, nc
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          
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
             Do no = 1,3
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
             Enddo
          Enddo
          ! compact c-electrons
          Do I = 1,Latt%N
             Do no = 4,5
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
             Enddo
          Enddo

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !  Setup the hopping
          !  Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: n, Ncheck, nc, ncoord, I, I1, I2

          If (checkerboard) then
            Ncheck = 3*Latt%N
            allocate(Op_T(Ncheck,N_FL))
            do n = 1,N_FL
              nc = 0
              Do ncoord = 1,3
                Do I = 1,Latt%N
                  nc = nc + 1
                  I1 = invlist(I,4)
                  if      ( ncoord == 1 ) then 
                      I2 = invlist(I,5)
                  elseif  ( ncoord == 2 ) then
                      I2 = Invlist( Latt%nnlist(I,1,-1),5 )
                  elseif  ( ncoord == 3 ) then
                      I2 = invlist( Latt%nnlist(I,0,-1),5 )
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
                enddo
              enddo
            enddo
          else
            Ncheck=1
            allocate(Op_T(Ncheck,N_FL))
            nc = 1
            do n = 1,N_FL
              Call Op_make(Op_T(nc,n),2*Latt%N)
              Do I = 1,2*Latt%N
                 Op_T(nc,n)%P(I)=3*Latt%N+I
              Enddo
              Do ncoord = 1,3
                Do I = 1,Latt%N
                  I1 = invlist(I,4)
                  if      ( ncoord == 1 ) then 
                      I2 = invlist(I,5)
                  elseif  ( ncoord == 2 ) then
                      I2 = Invlist( Latt%nnlist(I,1,-1),5 )
                  elseif  ( ncoord == 3 ) then
                      I2 = invlist( Latt%nnlist(I,0,-1),5 )
                  endif
                  I1 = I1 - 3*Latt%N
                  I2 = I2 - 3*Latt%N
                  Op_T(nc,n)%O( I1 , I2 ) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                  Op_T(nc,n)%O( I2 , I1 ) = cmplx(-Ham_T,   0.d0,Kind(0.d0))
                enddo
              enddo
              Op_T(nc,n)%g=cmplx(-Dtau,0.d0,Kind(0.d0))
              Call Op_set(Op_T(nc,n)) 
            enddo
          endif

        end Subroutine Ham_hop
!===================================================================================   

        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, I1, I2, nc, no
          Integer :: Ibond(2,6), Ihex(6)

          Complex (Kind=Kind(0.d0)) :: Zone


          ! Number of opertors 8 per unit cell
          Allocate( Op_V((6+6+3+15+3+3)*Latt%N,N_FL) )
          nc = 0
          Zone=cmplx(1.d0  ,0.d0, kind(0.D0))
          Do nf = 1,N_FL
            Do I = 1,Latt%N
              ! interaction on the bonds representing hopping of hard core bosons
              ! t1
              Ibond(1,1)=Invlist(I,1)
              Ibond(2,1)=Invlist(I,2)
              ! t2
              Ibond(1,2)=Invlist(I,2)
              Ibond(2,2)=Invlist(I,3)
              ! t3
              Ibond(1,3)=Invlist(I,3)
              Ibond(2,3)=Invlist(Latt%nnlist(I,0,1),1)
              ! t4
              Ibond(1,4)=Invlist(I,1)
              Ibond(2,4)=Invlist(Latt%nnlist(I,-1,0),2)
              ! t5
              Ibond(1,5)=Invlist(I,3)
              Ibond(2,5)=Invlist(Latt%nnlist(I,-1,1),2)
              ! t6
              Ibond(1,6)=Invlist(I,3)
              Ibond(2,6)=Invlist(I,1)
              Do no = 1,6
!                 write(*,*) Ibond(1,no), Ibond(2,no)
                nc = Latt%N*(no-1)+I
                Call Op_make(Op_V(nc,nf),2) 
                Op_V(nc,nf)%P(1) = Ibond(1,no)
                Op_V(nc,nf)%P(2) = Ibond(2,no)
                Op_V(nc,nf)%O( 1 , 2 )  =  Zone
                Op_V(nc,nf)%O( 2 , 1 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_T*0.5d0,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
                
                nc = Latt%N*(no-1+6)+I
                Call Op_make(Op_V(nc,nf),2) 
                Op_V(nc,nf)%P(1) = Ibond(1,no)
                Op_V(nc,nf)%P(2) = Ibond(2,no)
                Op_V(nc,nf)%O( 1 , 1 )  =  Zone
                Op_V(nc,nf)%O( 2 , 2 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_T*0.25d0,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(-1.d0, 0.d0, kind(0.D0))
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
              enddo
              
              ! "true" interaction of the total density per hexagon of the kagome lattice
              Ihex(1) = Invlist(I,3)
              Ihex(2) = Invlist(I,2)
              Ihex(3) = Invlist(Latt%nnlist(I,1,0),1)
              Ihex(4) = Invlist(Latt%nnlist(I,1,0),3)
              Ihex(5) = Invlist(Latt%nnlist(I,0,1),2)
              Ihex(6) = Invlist(Latt%nnlist(I,0,1),1)
!               write(*,*) Ihex(1), Ihex(2), Ihex(3), Ihex(4), Ihex(5), Ihex(6)
              no=12
              do I1=1,5
                do I2=I1+1,6
                  nc = Latt%N*no+I
                  Call Op_make(Op_V(nc,nf),2) 
                  Op_V(nc,nf)%P(1) = Ihex(I1)
                  Op_V(nc,nf)%P(2) = Ihex(I2)
                  Op_V(nc,nf)%O( 1 , 1 )  =  Zone
                  Op_V(nc,nf)%O( 2 , 2 )  =  -Zone
                  Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_Vint*0.25d0,0.d0, kind(0.D0)))
                  Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                  Op_V(nc,nf)%type   = 2
                  Call Op_set( Op_V(nc,nf) )
                  no = no+1
                enddo
              enddo
              
              ! auxiliary interaction to project out spin fluctuations reducing spinful fermions to effective hardcore bosons
              do no=1,3
                nc = Latt%N*(no+26)+I
!                 write(*,*) Invlist(I,no)
                Call Op_make(Op_V(nc,nf),1) 
                Op_V(nc,nf)%P(1) = Invlist(I,no)
                Op_V(nc,nf)%O( 1 , 1 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_U/beta,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(-0.5d0, 0.d0, kind(0.D0))
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
              enddo
              
              ! Kondo coupling of f- and c-electrons
              do no=1,3
                nc = Latt%N*(no+29)+I
                Call Op_make(Op_V(nc,nf),2) 
                Op_V(nc,nf)%P(1) = Invlist(I,no)
                Op_V(nc,nf)%P(2) = Invlist(I,4)
                Op_V(nc,nf)%O( 1 , 2 )  =  Zone
                Op_V(nc,nf)%O( 2 , 1 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_JKA/4.d0,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(0.0d0, 0.d0, kind(0.D0))
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
              enddo
              
              do no=1,3
                nc = Latt%N*(no+32)+I
                Call Op_make(Op_V(nc,nf),2) 
                if (no == 1) then
                  Op_V(nc,nf)%P(1) = Invlist(Latt%nnlist(I,0,1),1)
                elseif (no == 2 ) then
                  Op_V(nc,nf)%P(1) = Invlist(Latt%nnlist(I,-1,1),2)
                else
                  Op_V(nc,nf)%P(1) = Invlist(I,3)
                endif
                Op_V(nc,nf)%P(2) = Invlist(I,5)
                Op_V(nc,nf)%O( 1 , 2 )  =  Zone
                Op_V(nc,nf)%O( 2 , 1 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_JKB/4.d0,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(0.0d0, 0.d0, kind(0.D0))
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
              enddo
              
            enddo
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
          Allocate ( Obs_eq(10) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Ns = Latt%N;  No = 3;  Filename ="Green"
             case (2)
                Ns = Latt%N;  No = 3;  Filename ="SpinZ"
             case (3)
                Ns = Latt%N;  No = 3;  Filename ="SpinXY"
             case (4)
                Ns = Latt%N;  No = 3;  Filename ="Den"
             case (5)
                Ns = Latt%N;  No = 3;  Filename ="SC"
             case (6)
                Ns = Latt%N;  No = 2;  Filename ="ConductionGreen"
             case (7)
                Ns = Latt%N;  No = 2;  Filename ="ConductionSpinZ"
             case (8)
                Ns = Latt%N;  No = 2;  Filename ="ConductionSpinXY"
             case (9)
                Ns = Latt%N;  No = 2;  Filename ="ConductionDen"
             case (10)
                Ns = Latt%N;  No = 2;  Filename ="ConductionSC"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau == 1) then 
             ! Equal time correlators
             Allocate ( Obs_tau(10) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Ns = Latt%N; No = 3;  Filename ="Green"
                case (2)
                   Ns = Latt%N; No = 3;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N; No = 3;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N; No = 3;  Filename ="Den"
                case (5)
                   Ns = Latt%N; No = 3;  Filename ="SC"
                case (6)
                   Ns = Latt%N; No = 2;  Filename ="ConductionGreen"
                case (7)
                   Ns = Latt%N; No = 2;  Filename ="ConductionSpinZ"
                case (8)
                   Ns = Latt%N; No = 2;  Filename ="ConductionSpinXY"
                case (9)
                   Ns = Latt%N; No = 2;  Filename ="ConductionDen"
                case (10)
                   Ns = Latt%N; No = 2;  Filename ="ConductionSC"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1
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
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, Ztot
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZkinF, ZPot, Zocc, Z, ZP,ZS, weight, tmp,Zn
          Integer :: I,J, n, imj, nf, I1, J1, Nc, I0, I2, I3
          Integer :: K, K1, L ,L1, no_I, no_J, Ibond(2,6), Ihex(6)
          
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
          ZkinF = cmplx(0.d0,0.d0, kind(0.D0))
          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          weight=Ham_Vint*0.25d0
          Do nf = 1,N_FL
             Do I = 1,Latt%N
                ! interaction on the bonds representing hopping of hard core bosons
                ! t1
                Ibond(1,1)=Invlist(I,1)
                Ibond(2,1)=Invlist(I,2)
                ! t2
                Ibond(1,2)=Invlist(I,2)
                Ibond(2,2)=Invlist(I,3)
                ! t3
                Ibond(1,3)=Invlist(I,3)
                Ibond(2,3)=Invlist(Latt%nnlist(I,0,1),1)
                ! t4
                Ibond(1,4)=Invlist(I,1)
                Ibond(2,4)=Invlist(Latt%nnlist(I,-1,0),2)
                ! t5
                Ibond(1,5)=Invlist(I,3)
                Ibond(2,5)=Invlist(Latt%nnlist(I,-1,1),2)
                ! t6
                Ibond(1,6)=Invlist(I,3)
                Ibond(2,6)=Invlist(I,1)
                
                Do J = 1,6
                  I1=Ibond(1,J)
                  J1=Ibond(2,J)
                  Zkin = Zkin - Ham_T*(GRC(I1,J1,1)**2.d0+GRC(J1,I1,1)**2.d0)
                  ZkinF= ZkinF+ Ham_T*( 3.d0*Grc(I1,J1,1)*Gr(I1,J1,1) &
                      & -0.5d0*Grc(I1,I1,1)*(Grc(I1,I1,1)-1.d0) &
                      & -0.5d0*Grc(J1,J1,1)*(Grc(J1,J1,1)-1.d0) - 1.d0 )
                enddo
                do J=1,3
                  I1=Invlist(I,J)
                  ZkinF=ZkinF-Ham_U/beta*(2.d0*Grc(I1,I1,1)*(Grc(I1,I1,1)-1.d0)+1.d0)
                enddo
                
                
                Ihex(1) = Invlist(I,3)
                Ihex(2) = Invlist(I,2)
                Ihex(3) = Invlist(Latt%nnlist(I,1,0),1)
                Ihex(4) = Invlist(Latt%nnlist(I,1,0),3)
                Ihex(5) = Invlist(Latt%nnlist(I,0,1),2)
                Ihex(6) = Invlist(Latt%nnlist(I,0,1),1)
                tmp=0.d0
                Do J = 1,6
                  J1 = Ihex(J)
                  tmp = tmp + GRC(J1,J1,1)
                enddo
                Zpot = Zpot + (Zn*tmp)**2.d0
                Zpot = Zpot - 6.d0*(Zn**2.d0)*tmp
                Zpot = Zpot + 9.d0*(Zn**2.d0)
                Do J = 1,6
                  J1 = Ihex(J)
                  I1 = J1
                  Do K = 1,6
                    K1 = Ihex(K)
                    L1 = K1
                    tmp =  GRC(I1,L1,1) * GR (J1,K1,1)
                    ZPot  = ZPot  + tmp*Zn
                  Enddo
                ENddo
             Enddo
          Enddo
          ZPot = Ham_Vint*0.25d0*ZPot
          ZkinF = ZkinF + 3.d0*Latt%N*(2.d0*Ham_T+Ham_U/beta)
          
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
                ZTot  = ZTot  + weight*(Op_V(n,nf)%alpha**2.d0)*Zn
!           write(*,*) Zpot
             Enddo
          Enddo
          ZTot = ZTot*Zn + 9.d0*Latt%N*Ham_Vint
          Obs_scal(1)%Obs_vec(1)  =  Obs_scal(1)%Obs_vec(1) + Zkin * ZP*ZS
          Obs_scal(7)%Obs_vec(1)  =  Obs_scal(7)%Obs_vec(1) + ZkinF * ZP*ZS
          Obs_scal(8)%Obs_vec(1)  =  Obs_scal(8)%Obs_vec(1) + (Zkin + Zpot) * ZP*ZS
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =  Obs_scal(4)%Obs_vec(1) + ZTot*ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Zocc = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,3*Latt%N
                Zrho = Zrho + Grc(i,i,nf) 
                Zocc = Zocc + Grc(i,i,nf)**2
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(5)%Obs_vec(1)  =    Obs_scal(5)%Obs_vec(1) - 2.d0*(Zocc - 0.5d0*Zrho) * ZP*ZS
          
          ! conduction fermions
          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Zocc = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 3*Latt%N+1,5*Latt%N
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
                I0 = invlist(I,4)
                I1 = invlist(I,5)
                I2 = Invlist( Latt%nnlist(I,1,-1),5 )
                I3 = invlist( Latt%nnlist(I,0,-1),5 )  
                Zkin = Zkin + Grc(I0,I1,nf) +  Grc(I1,I0,nf)  + &
                     &        Grc(I0,I2,nf) +  Grc(I2,I0,nf)  + &
                     &        Grc(I0,I3,nf) +  Grc(I3,I0,nf)
             Enddo
          Enddo
          Zkin = - Zkin*cmplx( Ham_T*dble(N_SUN), 0.d0,Kind(0.d0) )
          Obs_scal(9)%Obs_vec(1)  =    Obs_scal(9)%Obs_vec(1) + Zkin * ZP*ZS
          !add kin part of conduction electrons to total energy
          Obs_scal(4)%Obs_vec(1)  =  Obs_scal(4)%Obs_vec(1) + Zkin*ZP*ZS
          
          ZPot = cmplx(0.d0,0.d0,Kind(0.d0))
          Do nf = 1,N_FL
             Do I = 1,Latt%N
               Do J = 1,3
                I0 = invlist(I,4)
                I1 = invlist(I,J)  
                Zpot = ZPot + Grc(I0,I1,nf)**2 +  Grc(I1,I0,nf)**2  + &
                     &        Grc(I0,I0,nf) +  Grc(I1,I1,nf)   &
                     &        -2.d0*Grc(I0,I0,nf)*Grc(I1,I1,nf) +  zn*Grc(I0,I1,nf)*Gr(I0,I1,nf)
               enddo
             Enddo
          Enddo
          Zpot = - ZPot*cmplx( Ham_JKA*dble(N_SUN), 0.d0,Kind(0.d0) )
          Obs_scal(11)%Obs_vec(1)  =    Obs_scal(11)%Obs_vec(1) + ZPot * ZP*ZS
          ZPot = cmplx(0.d0,0.d0,Kind(0.d0))
          Do nf = 1,N_FL
             Do I = 1,Latt%N
               Do J = 1,3
                I0 = invlist(I,4)
                if (J == 1) then
                  I1 = Invlist(Latt%nnlist(I,0,1),1)
                elseif (J == 2 ) then
                  I1 = Invlist(Latt%nnlist(I,-1,1),2)
                else
                  I1 = Invlist(I,3)
                endif  
                Zpot = ZPot + Grc(I0,I1,nf)**2 +  Grc(I1,I0,nf)**2  + &
                     &        Grc(I0,I0,nf) +  Grc(I1,I1,nf)   &
                     &        -2.d0*Grc(I0,I0,nf)*Grc(I1,I1,nf) +  zn*Grc(I0,I1,nf)*Gr(I0,I1,nf)
               enddo
             Enddo
          Enddo
          Zpot = - ZPot*cmplx( Ham_JKB*dble(N_SUN), 0.d0,Kind(0.d0) )
          Obs_scal(11)%Obs_vec(1)  =    Obs_scal(11)%Obs_vec(1) + ZPot * ZP*ZS
          
          
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          
          Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
          Do J1 = 1,3*Latt%N
            J    = List(J1,1)
            no_J = List(J1,2)
            Do I1 = 1,3*Latt%N
                I    = List(I1,1)
                no_I = List(I1,2)
                imj = latt%imj(I,J)
                ! Green
                Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               Z * GRC(I1,J1,1) *  ZP*ZS 
                ! SpinZ
                Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! SpinXY
                Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! Den
                Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                    &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                    &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                    &                                   ) * Z* ZP*ZS
                ! SC
                Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               GRC(I1,J1,1)**2 *  ZP*ZS 
                Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1) + GRC(I1,J1,1)**2/ 3.d0 /dble(Latt%N) * ZP*ZS
            ENDDO
            Obs_eq(4)%Obs_Latt0(no_J) =  Obs_eq(4)%Obs_Latt0(no_J) +  Z * GRC(J1,J1,1) * ZP * ZS
          ENDDO
          Do J1 = 3*Latt%N+1,5*Latt%N
            J    = List(J1,1)
            no_J = List(J1,2)-3
            Do I1 = 3*Latt%N+1,5*Latt%N
                I    = List(I1,1)
                no_I = List(I1,2)-3
                imj = latt%imj(I,J)
                ! Green
                Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               Z * GRC(I1,J1,1) *  ZP*ZS 
                ! SpinZ
                Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! SpinXY
                Obs_eq(8)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(8)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! Den
                Obs_eq(9)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(9)%Obs_Latt(imj,1,no_I,no_J)  +  &
                    &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                    &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                    &                                   ) * Z* ZP*ZS
                ! SC
                Obs_eq(10)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(10)%Obs_Latt(imj,1,no_I,no_J) + &
                    &               GRC(I1,J1,1)**2 *  ZP*ZS 
!                 Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1) + GRC(I1,J1,1)**2/ 3.d0 /dble(Latt%N) * ZP*ZS
            ENDDO
            Obs_eq(9)%Obs_Latt0(no_J) =  Obs_eq(9)%Obs_Latt0(no_J) +  Z * GRC(J1,J1,1) * ZP * ZS
          ENDDO

        end Subroutine Obser
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
          Do I1 = 1,3*Latt%N
            I    = List(I1,1)
            no_I = List(I1,2)
            Do J1 = 1,3*Latt%N
                J    = List(J1,1)
                no_J = List(J1,2)
                imj = latt%imj(I,J)
                ! Green
                Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & +  Z * GT0(I1,J1,1) * ZP* ZS
                
                ! SpinZ
                Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                
                ! SpinXY
                Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                
                ! Den
                Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
                    &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
                    &     Z * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                
                ! SC
                Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & + (    G0T(J1,I1,1)**2.d0                                ) * ZP * ZS
            Enddo
            Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
                  &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
          Enddo
          Do I1 = 3*Latt%N+1,5*Latt%N
            I    = List(I1,1)
            no_I = List(I1,2)-3
            Do J1 = 3*Latt%N+1,5*Latt%N
                J    = List(J1,1)
                no_J = List(J1,2)-3
                imj = latt%imj(I,J)
                ! Green
                Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & +  Z * GT0(I1,J1,1) * ZP* ZS
                
                ! SpinZ
                Obs_tau(7)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(7)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                
                ! SpinXY
                Obs_tau(8)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(8)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                
                ! Den
                Obs_tau(9)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(9)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
                    &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
                    &     Z * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                
                ! SC
                Obs_tau(10)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(10)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & + (    G0T(J1,I1,1)**2.d0                                ) * ZP * ZS
            Enddo
            Obs_tau(9)%Obs_Latt0(no_I) = Obs_tau(9)%Obs_Latt0(no_I) + &
                  &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
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
          Integer :: I, length, start, t, n
          T0_Proposal_ratio=1.d0
          
          I=nranf(Latt%N)
          start=nranf(Ltrot)
          length=nranf(Ltrot)
!           write(*,*) I, start, length
          
          nsigma=nsigma_old
          do t=start,min(Ltrot,start+length-1)
!           write(*,*) t
          do n=1,15
            nsigma(Latt%N*(12+n)+I,t)=-nsigma_old(Latt%N*(12+n)+I,t)
          enddo
          enddo
          do t=1,length-(Ltrot-start+1)
!           write(*,*) t
          do n=1,15
            nsigma(Latt%N*(12+n)+I,t)=-nsigma_old(Latt%N*(12+n)+I,t)
          enddo
          enddo
          
          size_clust=dble(length)/dble(Ltrot)
          
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
          Integer                :: I
          
          Flip_length=15*Latt%N
          Do I = 1,Flip_length
            Flip_list(I)=I+12*Latt%N
            Flip_value(I)=-nsigma(I+12*Latt%N,ntau)
          Enddo
          T0_Proposal_ratio=1.d0
          S0_ratio=1.d0

        end Subroutine Global_move_tau


    end Module Hamiltonian
