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
      
      Type (Lattice),       private :: Latt
      Integer, parameter,   private :: Norb=3
      Integer, allocatable, private :: List(:,:), Invlist(:,:)
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: Ham_T, Ham_Vint,  Ham_U
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta
      Character (len=64),   private :: Model, Lattice_type, File1


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau

      contains 

        Subroutine Ham_Set

          Implicit none
#if defined(MPI)
          include 'mpif.h'
#endif   

          integer :: ierr

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_Z2_FLstar/  ham_T, Ham_Vint,  Ham_U,  Dtau, Beta


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
             OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_Z2_FLstar)
             CLOSE(5)
#if defined(MPI)
          endif

          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_Vint ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,   0,Group_Comm,ierr)
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
             Open (Unit = 50,file="info",status="unknown",position="append")
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
          Do I = 1,Latt%N
             Do no = 1,Norb
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

          Integer :: n, Ncheck, nc


          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,NCheck
                Call Op_make(Op_T(nc,n),1)
                Op_T(nc,n)%P(1) = 1
                Op_T(nc,n)%O( 1 , 1 )  =  1.d0
                Op_T(nc,n)%g=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n)) 
             enddo
          enddo

        end Subroutine Ham_hop
!===================================================================================   

        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, I1, I2, nc, no
          Integer :: Ibond(2,6), Ihex(6)

          Complex (Kind=Kind(0.d0)) :: Zone


          ! Number of opertors 8 per unit cell
          Allocate( Op_V((6+6+3+15)*Latt%N,N_FL) )
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
                nc = nc + 1
                Call Op_make(Op_V(nc,nf),2) 
                Op_V(nc,nf)%P(1) = Ibond(1,no)
                Op_V(nc,nf)%P(2) = Ibond(2,no)
                Op_V(nc,nf)%O( 1 , 2 )  =  Zone
                Op_V(nc,nf)%O( 2 , 1 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_T*0.5d0,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
                
                nc = nc + 1
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
              do I1=1,5
                do I2=I1+1,6
                  nc = nc + 1
                  Call Op_make(Op_V(nc,nf),2) 
                  Op_V(nc,nf)%P(1) = Ihex(I1)
                  Op_V(nc,nf)%P(2) = Ihex(I2)
                  Op_V(nc,nf)%O( 1 , 1 )  =  Zone
                  Op_V(nc,nf)%O( 2 , 2 )  =  -Zone
                  Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_Vint*0.25d0,0.d0, kind(0.D0)))
                  Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                  Op_V(nc,nf)%type   = 2
                  Call Op_set( Op_V(nc,nf) )
                enddo
              enddo
              
              ! auxiliary interaction to project out spin fluctuations reducing spinful fermions to effective hardcore bosons
              do no=1,3
                nc = nc + 1
!                 write(*,*) Invlist(I,no)
                Call Op_make(Op_V(nc,nf),1) 
                Op_V(nc,nf)%P(1) = Invlist(I,no)
                Op_V(nc,nf)%O( 1 , 1 )  =  Zone
                Op_V(nc,nf)%g=sqrt(cmplx(Dtau*Ham_U/beta,0.d0, kind(0.D0)))
                Op_V(nc,nf)%alpha  = cmplx(-0.5d0, 0.d0, kind(0.D0))
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
          Allocate ( Obs_scal(6) )
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
                Ns = Latt%N;  No = Norb;  Filename ="SC"
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
! <<<<<<< HEAD
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZkinF, ZPot, Zocc, Z, ZP,ZS, weight, tmp,Zn
          Integer :: I,J, n, imj, nf, I1, J1, Nc
          Integer :: K, K1, L ,L1, no_I, no_J, Ibond(2,6), Ihex(6)
! =======
!           Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, ZkinF
!           Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Zocc, Z, ZP,ZS, weight, tmp, alpha1, alpha2, Zn
!           Integer :: I,J, no,no1, n, n1, imj, nf, I1, I2, J1, J2, Nc, Ix, Iy, Jx, Jy, Imx, Imy, Jmx, Jmy
!           Integer :: a, b, c, d, signum, K, K1, L ,L1, nf1, no_I, no_J, Ibond(2,6), Ihex(6)
!           
!           Real (Kind=Kind(0.d0)) :: G(4,4), X, FI, FJ
! >>>>>>> 4104e55fe99e7a828dc61225a2dee9772cf62be9
          
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
          Nc = Size( Op_V,1)
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
!                   ZkinF= ZkinF+ Ham_T*( 3.d0*Grc(I1,J1,1)*Gr(I1,J1,1) &
!                       & -0.5d0*Grc(I1,I1,1)*(Grc(I1,I1,1)-1.d0) &
!                       & -0.5d0*Grc(J1,J1,1)*(Grc(J1,J1,1)-1.d0) - 1.d0 )
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
! !           Zkin=Zn*Zkin
!           ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
!           Do nf = 1,N_FL
!              Do n = 1,Nc
!                 weight=-Op_V(n,nf)%g**2.d0 /dtau !(-1)**((n-1)/Latt%N)/8.d0
! !                 write(0,*) Op_V(n,nf)%g, weight
!                 Do J = 1,Op_V(n,nf)%N
!                    J1 = Op_V(n,nf)%P(J)
!                    DO I = 1,Op_V(n,nf)%N
!                       if (abs(Op_V(n,nf)%O(i,j)) >= 0.00001) then
!                       I1 = Op_V(n,nf)%P(I)
!                       ZPot  = ZPot  + 2.d0*Zn*Op_V(n,nf)%alpha*weight*Op_V(n,nf)%O(i,j)*GRC(I1,J1,1)
!                       Do K = 1,Op_V(n,nf)%N
!                         K1 = Op_V(n,nf)%P(K)
!                         DO L = 1,Op_V(n,nf)%N
!                           if (abs(Op_V(n,nf)%O(k,l)) >= 0.00001) then
!                           L1 = Op_V(n,nf)%P(L)
!                           tmp =  (   GRC(I1,L1,1) * GR (J1,K1,1)      +  &
!                                 &    Zn * GRC(I1,J1,1) * GRC(K1,L1,1)         )
!                           ZPot  = ZPot  + weight*Op_V(n,nf)%O(i,j)*Op_V(n,nf)%O(k,l)*tmp
!                           endif
!                         Enddo
!                       ENddo
!                       endif
!                    Enddo
!                 ENddo
!                 ZPot  = ZPot  + weight*(Op_V(n,nf)%alpha**2.d0)*Zn
! !           write(*,*) Zpot
!              Enddo
!           Enddo
!           ZPot = ZPot*Zn + 9.d0*Latt%N*Ham_Vint
          Obs_scal(1)%Obs_vec(1)  =  Obs_scal(1)%Obs_vec(1) + Zkin * ZP*ZS
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Zocc = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
                Zocc = Zocc + Grc(i,i,nf)**2
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(5)%Obs_vec(1)  =    Obs_scal(5)%Obs_vec(1) - 2.d0*(Zocc - 0.5d0*Zrho) * ZP*ZS
          
          
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          
          Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
          Do I1 = 1,Ndim
            I    = List(I1,1)
            no_I = List(I1,2)
            Do J1 = 1,Ndim
                J    = List(J1,1)
                no_J = List(J1,2)
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
            Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
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
            Enddo
            Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
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
