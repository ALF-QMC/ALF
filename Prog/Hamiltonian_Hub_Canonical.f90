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
      Integer,              private :: L1, L2, n_Lambda
      real (Kind=8),        private :: ham_T , ham_U,  Ham_chem, Lambda, RhoD
      real (Kind=8),        private :: Dtau, Beta, Phi_x
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord, Norb 
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
    contains 


      Subroutine Ham_Set

          Implicit none

#ifdef MPI
          include 'mpif.h'
#endif   

          integer :: ierr

          
          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model, Phi_x

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U, Dtau, Beta, Lambda, n_Lambda, RhoD


#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
          
          


#ifdef MPI
          If (Irank == 0 ) then
#endif
             Phi_x = 0.d0
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
          CALL MPI_BCAST(Phi_x       ,1  ,MPI_REAL8  ,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          If ( Model == "Hubbard_Mz") then
             N_FL = 2
             N_SUN = 1
          elseif  ( Model == "Hubbard_SU2" ) then
             N_FL = 1
             N_SUN = 2
          else
             Write(6,*) "Model not yet implemented!"
             Stop
          endif
#ifdef MPI
          If (Irank == 0 ) then
#endif
             RhoD = 1.d0
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_Hubbard)
             CLOSE(5)
#ifdef MPI
          endif
          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_chem ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_U    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Lambda   ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(n_Lambda ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(RhoD     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
#endif
          Call Ham_hop

          Ltrot = nint(beta/dtau)
#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '=============Canonical========='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'L1,L2         : ', L1,L2
             Write(50,*) 'Phi_x         : ', Phi_x
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'N_SUN         : ', N_SUN
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_U         : ', Ham_U
             Write(50,*) 'Ham_chem      : ', Ham_chem
             Write(50,*) 'Desired Den   : ', RhoD
             Write(50,*) 'Lambda        : ', Lambda
             Write(50,*) 'n_Lambda      : ', n_Lambda
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
          Real (Kind=8)  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no


          If ( Lattice_type =="Square" ) then
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             !Write(6,*)  'Lattice: ', Ndim
             One_dimensional = .false.
             N_coord   = 2
             Norb = 1
             If ( L1 == 1 .or. L2 == 1 ) then 
                One_dimensional = .true.
                N_coord   = 1
                If (L1 == 1 ) then 
                   Write(6,*) ' For one dimensional systems set  L2 = 1 ' 
                   Stop
                endif
             endif
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif
          
          ! This is for the orbital structure.
          Ndim = Latt%N*Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                ! For the Honeycomb lattice no = 1,2 corresponds to the A,and B sublattices.
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

          !Setup the hopping
          !Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, I2, n, Ncheck,nc
          Complex  (Kind=8) :: Z
          Real     (Kind=8) :: Pi

          Ncheck = 1
          Pi = acos(-1.d0)
          Z = exp( cmplx(0.d0,Phi_x*2.d0*pi/real(L1,kind(0.d0)),kind(0.d0)))
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,Ncheck
                Call Op_make(Op_T(nc,n),Ndim)
                If (One_dimensional ) then 
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I, 1, 0)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Z
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Conjg(Z)
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                   ENDDO
                else
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I,1,0)
                      I2 = Latt%nnlist(I,0,1)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))*Z
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))*Conjg(Z)
                      Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                   ENDDO
                endif
                
                Do I = 1,Latt%N
                   Op_T(nc,n)%P(i) = i 
                Enddo
                if ( abs(Ham_T) < 1.E-6 .and.  abs(Ham_chem) < 1.E-6 ) then 
                   Op_T(nc,n)%g = 0.d0
                else
                   Op_T(nc,n)%g = -Dtau
                endif
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                !Write(6,*) 'In Ham_hop', Ham_T
                Call Op_set(Op_T(nc,n)) 
                !Write(6,*) 'In Ham_hop 1'
                !Do I = 1,Latt%N
                !   Write(6,*) Op_T(n)%E(i)
                !enddo
                !Call Op_exp( cmplx(-Dtau,0.d0), Op_T(n), Exp_T   (:,:,n) )
                !Call Op_exp( cmplx( Dtau,0.d0), Op_T(n), Exp_T_M1(:,:,n) )
             enddo
          enddo
        end Subroutine Ham_hop

!===================================================================================           
        
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, nc, nth
          Real (Kind=8) :: X

          
          If (Model == "Hubbard_SU2")  then
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Latt%N+n_Lambda,N_FL))
             do nf = 1,N_FL
                do i  = 1, Latt%N
                   Call Op_make(Op_V(i,nf),1)
                enddo
                do nth = 1,n_lambda
                   Call Op_make(Op_V(Latt%N+nth,nf),Latt%N)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                Do i = 1,Latt%N
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                   ! The operator reads:  
                   !  g*s*( c^{dagger} O c  + alpha ))
                   ! with s the HS field.
                Enddo
                Do nth = 1,n_Lambda
                   nc = nc + 1
                   Do  I = 1,latt%N
                      Op_V(nc,nf)%P(I) = I
                      Op_V(nc,nf)%O(I,I) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   enddo
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*Lambda/dble(n_lambda), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-DBLE(Latt%N)*0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Elseif (Model == "Hubbard_Mz")  then
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Latt%N + n_lambda ,N_FL))
             do nf = 1,N_FL
                do i  = 1, Latt%N
                   Call Op_make(Op_V(i,nf),1)
                enddo
                do nth = 1,n_lambda
                   Call Op_make(Op_V(Latt%N+nth,nf),Latt%N)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                X = 1.d0
                if (nf == 2) X = -1.d0
                Do i = 1,Latt%N
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                   ! The operator reads:  
                   ! g*s*( c^{dagger} O c   - alpha ))
                   ! with s the HS field.
                   ! Write(6,*) nc,nf, Op_V(nc,nf)%g 
                Enddo
                Do nth = 1,n_Lambda
                   nc = nc + 1
                   Do  I = 1,latt%N
                      Op_V(nc,nf)%P(I) = I
                      Op_V(nc,nf)%O(I,I) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   enddo
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*Lambda/dble(n_lambda), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-DBLE(Latt%N)*RhoD*0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                enddo
             Enddo
          Endif
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=8) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt 
          S0 = 1.d0
          
        end function S0

        Subroutine  Alloc_obs(Ltau) 

          Implicit none
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
          Allocate ( Obs_eq(4) )
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

        
!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J
          
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
          Do nf = 1,N_FL
             Do J = 1,Ndim
                Zkin = Zkin + sum(Op_T(1,nf)%O(:, j)*Grc(:, j, nf))
             ENddo
          Enddo
          Zkin = Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" ) then
            dec = 1
          else
            dec = 2
          endif
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, dec)
          Enddo
          Zpot = Zpot*ham_U
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


          
          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then 
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
                ENDDO
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
             ENDDO
          elseif (Model == "Hubard_Mz" ) Then
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
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &    (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))    ) * ZP*ZS
                   ! SpinXY
                   ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                     & (   GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1)    ) * ZP*ZS
                   !Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &   (GRC(I1,I1,2) + GRC(I1,I1,1))*(GRC(J1,J1,2) + GRC(J1,J1,1))     ) * ZP*ZS
                enddo
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  (GRC(I1,I1,1) + GRC(I1,I1,2)) * ZP*ZS
             enddo
          Endif
                

        end Subroutine Obser
!=====================================================
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then 
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
          Elseif ( Model == "Hubbard_Mz"  ) then 
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   !Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &   +   ( GT0(I1,J1,1) + GT0(I1,J1,2) ) * ZP* ZS

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
          Endif
          
        end Subroutine OBSERT

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
