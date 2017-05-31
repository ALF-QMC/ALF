!  This is for the Kondo project with z-frustration.  Kondo on the Honeycomb lattice.
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

      
      ! What is below is  private 
      
      Type (Lattice),       private :: Latt
      Integer, parameter,   private :: Norb=4
      Integer, allocatable, private :: List(:,:), Invlist(:,:)
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T, Ham_U,  Ham_J, Ham_Jz, del_p(2)
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta
      Integer,              private :: Checkerboard
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord 
      Real (Kind=Kind(0.d0)),        private :: Bound


      ! Observables
      Integer,                       private :: Nobs
      Complex (Kind=Kind(0.d0)), allocatable, private :: obs_scal(:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Spinz_eq(:,:,:), Spinz_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Spinxy_eq(:,:,:),Spinxy_eq0(:) 
      Complex (Kind=Kind(0.d0)), allocatable, private :: Den_eq(:,:,:), Den_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Dimer_eq(:,:,:), Dimer_eq0(:)

!-------------add-----------------------------------------------------------------
      Complex (Kind=Kind(0.d0)), allocatable, private :: Greenu_eq(:,:,:), Greenu_eq0(:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Greend_eq(:,:,:), Greend_eq0(:)
!---------------------------------------------------------------------------------


      ! For time displaced
      Integer,                       private :: NobsT
      Complex (Kind=Kind(0.d0)),              private :: Phase_tau
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green_tau(:,:,:,:), Den_tau(:,:,:,:)

      contains 

        Subroutine Ham_Set

          Implicit none
#ifdef MPI
          include 'mpif.h'
#endif   

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
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             !Write(6,*)  'Lattice: ', Ndim
             One_dimensional = .false.
             N_coord   = 2
             If ( L1 == 1 .or. L2 == 1 ) then 
                One_dimensional = .true.
                N_coord   = 1
                If (L1 == 1 ) then 
                   Write(6,*) ' For one dimensional systems set  L2 = 1 ' 
                   Stop
                endif
             endif
          elseif ( Lattice_type=="Honeycomb" ) then
             
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.5  ; a2_p(2) =  sqrt(3.0)/2.0
             del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
             
             L1_p    =  dble(L1) * a1_p
             L2_p    =  dble(L2) * a2_p
             
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )

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
          
          Implicit none 
          
          Integer :: nf, nth, n, n1, n2, n3, n4, I, I1, I2, J,  Ix, Iy, nc, no
          Real (Kind=Kind(0.d0)) :: X_p(2), X1_p(2), X2_p(2), X, XJ


          ! Number of opertors. 
          ! Latt%N * 2 for the Hubbard
          ! Latt%N * 2 For the Kondo  
          ! Latt%N * 6 for the Frustrating Ising   
          ! Total of Latt%N * 10
          
          Allocate( Op_V(10*Latt%N,N_FL) )
          do nf = 1,N_FL
             do i  = 1, 2*Latt%N
                Call Op_make(Op_V(i,nf),1)  ! For the Hubbard. 
             enddo
          enddo
          do nf = 1,N_fl
             do i = 2*Latt%N +1, 10*Latt%N
                Call Op_make(Op_V(i,nf),2)  !  For the Kondo and  Ising. 
             enddo
          enddo
          
          Do nf = 1,N_FL
             X = 1.d0
             if (nf == 2) X = -1.d0
             nc = 0
             !  Hubbard
             Do i = 1,Latt%N
                Do no = 3,4
                   nc = nc + 1 
                   Op_V(nc,nf)%P(1)   = Invlist(I,no)  ! f-site
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0,Kind(0.d0))
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/2.d0,0.D0,Kind(0.d0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0,Kind(0.d0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                   ! The operator reads:  
                   !  g*s*( c^{dagger} O c   - alpha ))
                   ! with s the HS field.
                Enddo
             Enddo
             !  Kondo
             DO i = 1,Latt%N
                Do no = 1,2
                   nc = nc + 1 
                   Op_V(nc,nf)%P(1) = Invlist(I,no  )
                   Op_V(nc,nf)%P(2) = Invlist(I,no+2)
                   Op_V(nc,nf)%O(1,2) = cmplx(1.d0  ,0.d0,Kind(0.d0))
                   Op_V(nc,nf)%O(2,1) = cmplx(1.d0  ,0.d0,Kind(0.d0))
                   Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_J/4.d0,0.D0,Kind(0.d0))) 
                   Op_V(nc,nf)%alpha  = cmplx(0.d0,0.d0,Kind(0.d0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
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
                      Op_V(nc,nf)%P(1) = Invlist(I,no)
                      Op_V(nc,nf)%P(2) = Invlist(J,no)
                      Op_V(nc,nf)%O(1,1) = cmplx( 1.d0  ,0.d0,Kind(0.d0))
                      Op_V(nc,nf)%O(2,2) = cmplx(-1.d0  ,0.d0,Kind(0.d0))
                      Op_V(nc,nf)%g      = X*SQRT(CMPLX(DTAU*ham_Jz/2.d0,0.D0,Kind(0.d0))) 
                      Op_V(nc,nf)%alpha  = cmplx(0.d0,0.d0,Kind(0.d0))
                      Op_V(nc,nf)%type   = 2
                      Call Op_set( Op_V(nc,nf) )
                      ! The operator reads:  
                      !  g*s*( c^{dagger} O c   - alpha ))
                      ! with s the HS field.
                   Enddo
                Enddo
             Enddo
          Enddo
          
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt 
          Integer :: i, nt1 
          S0 = 1.d0
        end function S0
!===================================================================================           
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer :: I
          Allocate ( Obs_scal(5) )
          Allocate ( Spinz_eq(Latt%N,Norb,Norb) , Spinz_eq0(Norb)   )
          Allocate ( Spinxy_eq(Latt%N,Norb,Norb), Spinxy_eq0(Norb) ) 
          Allocate ( Den_eq(Latt%N,Norb,Norb),    Den_eq0(Norb)       ) 

!-------------add-----------------------------------------------------------------
          Allocate ( Greenu_eq(Latt%N,Norb,Norb),  Greenu_eq0(Norb)     ) 
          Allocate ( Greend_eq(Latt%N,Norb,Norb),  Greend_eq0(Norb)     )
!---------------------------------------------------------------------------------

          If (Ltau == 1) then 
             Allocate ( Green_tau(Latt%N,Ltrot+1,Norb,Norb), Den_tau(Latt%N,Ltrot+1,Norb,Norb) )
          endif
          
        end Subroutine Alloc_obs

!===================================================================================           
        
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          Integer :: I,n
          
          Nobs = 0
          Obs_scal  = cmplx(0.d0,0.d0,Kind(0.d0))
          SpinZ_eq  = cmplx(0.d0,0.d0,Kind(0.d0)) 
          SpinZ_eq0 = cmplx(0.d0,0.d0,Kind(0.d0)) 
          Spinxy_eq = cmplx(0.d0,0.d0,Kind(0.d0)) 
          Spinxy_eq0= cmplx(0.d0,0.d0,Kind(0.d0)) 
          Den_eq    = cmplx(0.d0,0.d0,Kind(0.d0))
          Den_eq0   = cmplx(0.d0,0.d0,Kind(0.d0))

!-------------add-----------------------------------------------------------------
          Greenu_eq  = cmplx(0.d0,0.d0,Kind(0.d0))
          Greenu_eq0 = cmplx(0.d0,0.d0,Kind(0.d0)) 
          Greend_eq  = cmplx(0.d0,0.d0,Kind(0.d0))
          Greend_eq0 = cmplx(0.d0,0.d0,Kind(0.d0)) 
!---------------------------------------------------------------------------------

          If (Ltau == 1) then
             NobsT = 0
             Phase_tau = cmplx(0.d0,0.d0,Kind(0.d0))
             Green_tau = cmplx(0.d0,0.d0,Kind(0.d0))
             Den_tau = cmplx(0.d0,0.d0,Kind(0.d0))
          endif

        end Subroutine Init_obs
        
!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, G(4,4,N_FL)
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS
          Integer :: I,J, no,no1, n, n1, imj, nf, I1, I2,I3, J1, J2, I0, J0, ns, nc, NC_tot
          
          Real (Kind=Kind(0.d0)) ::  X
          
          Nobs = Nobs + 1
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
          !Nc_tot = Size(OP_T,1)
          !Write(6,*) 'Obser', Nc_tot
          !Do nf = 1,N_FL
          !   do nc = 1, Nc_tot
          !      Do J = 1,Op_T(nc,nf)%N
          !         J1 = Op_T(nc,nf)%P(J)
          !         DO I = 1,Op_T(nc,nf)%N
          !            I1 = Op_T(nc,nf)%P(I)
          !            Zkin  = Zkin  + Op_T(nc,nf)%O(i,j)*Grc(i1,j1,nf) 
          !         Enddo
          !      ENddo
          !   enddo
          !enddo
          !Write(6,*) 'End Compute Kin: ', Size(OP_T,1), Size(OP_T,2)
          

          Zrho = cmplx(0.d0,0.d0,Kind(0.d0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho*cmplx( dble(N_SUN), 0.d0,Kind(0.d0) )

          ZPot = cmplx(0.d0,0.d0,Kind(0.d0))
          Do no = 3,4
             Do I = 1,Latt%N
                I1 = Invlist(I,no)
                ZPot = ZPot + Grc(I1,I1,1) * Grc(I1,I1,2)
             Enddo
          Enddo
          Zpot = Zpot!*cmplx(ham_U,0.d0,Kind(0.d0))

          Obs_scal(1) = Obs_scal(1) + zrho * ZP*ZS
          Obs_scal(2) = Obs_scal(2) + zkin * ZP*ZS
          Obs_scal(3) = Obs_scal(3) + Zpot * ZP*ZS
          Obs_scal(4) = Obs_scal(4) + (zkin +  Zpot)*ZP*ZS
          Obs_scal(5) = Obs_scal(5) + ZS
          ! You will have to allocate more space if you want to include more  scalar observables.



          DO I1 = 1,Ndim
             I = List(I1,1)
             no = List(I1,2)
             DO J1 = 1, Ndim
                J = List(J1,1)
                no1 = list(J1,2)
                imj = latt%imj(I,J)
                SPINZ_Eq (imj,no,no1) = SPINZ_Eq (imj,no,no1)  +  &
                     & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                     &   (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))    ) * ZP*ZS
                ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                SPINXY_Eq (imj,no,no1) = SPINXY_Eq (imj,no,no1)  +  &
                     & (   GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1)    ) * &
                     &   cmplx(2.d0,0.d0,Kind(0.d0))* ZP*ZS

                DEN_Eq (imj,no,no1) = DEN_Eq (imj,no,no1)  +  &
                     & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                     &   (GRC(I1,I1,2) + GRC(I1,I1,1))*(GRC(J1,J1,2) + GRC(J1,J1,1))    ) * ZP*ZS

!-------------add-----------------------------------------------------------------
                Greenu_eq(imj,no,no1)=Greenu_eq(imj,no,no1)+GRC(I1,J1,1)*ZP*ZS
                Greend_eq(imj,no,no1)=Greend_eq(imj,no,no1)+GRC(I1,J1,2)*ZP*ZS
!---------------------------------------------------------------------------------

             enddo
             Den_eq0(no) = Den_eq0(no) + (GRC(I1,I1,2) + GRC(I1,I1,1)) * ZP*ZS
          enddo

        end Subroutine Obser
!==========================================================        
        Subroutine  Pr_obs(LTAU)


          Use Print_bin_mod
          Implicit none
#ifdef MPI
          include 'mpif.h'
#endif   


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
    
          Phase_bin = Obs_scal(5)/cmplx(dble(Nobs),0.d0,Kind(0.d0))

          File_pr ="SpinZ_eq"
          Call Print_bin(SpinZ_eq ,SpinZ_eq0, Latt, Nobs, Phase_bin, file_pr)

          File_pr ="SpinXY_eq"
          Call Print_bin(Spinxy_eq, Spinxy_eq0,Latt, Nobs, Phase_bin, file_pr)

          File_pr ="Den_eq"
          Call Print_bin(Den_eq   , Den_eq0, Latt, Nobs, Phase_bin, file_pr)

!-------------add-----------------------------------------------------------------
          File_pr ="Greenu_eq"
          Call Print_bin(Greenu_eq   , Greenu_eq0, Latt, Nobs, Phase_bin, file_pr)
          File_pr ="Greend_eq"
          Call Print_bin(Greend_eq   , Greend_eq0, Latt, Nobs, Phase_bin, file_pr)
!---------------------------------------------------------------------------------

          File_pr ="ener"
          Call Print_scal(Obs_scal, Nobs, file_pr)
          If (Ltau == 1) then
             Phase_tau = Phase_tau/cmplx(dble(NobsT),0.d0,Kind(0.d0))
             File_pr = "Green_tau"
             Call Print_bin_tau(Green_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
             File_pr = "Den_tau"
             Call Print_bin_tau(Den_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
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
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J

          ZP = PHASE/cmplx(Real(Phase,Kind=Kind(0.d0)),0.d0,Kind(0.d0))
          ZS = cmplx(Real(Phase,Kind=Kind(0.d0))/Abs(Real(Phase,Kind=Kind(0.d0))), 0.d0,Kind(0.d0))
          If (NT == 0 ) then 
             Phase_tau = Phase_tau + ZS
             NobsT     = NobsT + 1
          endif
          If ( N_FL == 1 ) then 
          
             Z =  ZP * ZS * N_SUN
             Do I = 1,Latt%N
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   Z = Z * GT0(I, J, 1)
                   Green_tau(imj,nt+1,1,1) = green_tau(imj,nt+1,1,1)  +  Z
                   Den_tau  (imj,nt+1,1,1) = Den_tau  (imj,nt+1,1,1)  -  Z * G0T(J,I,1)
                Enddo
             Enddo
          Endif
        end Subroutine OBSERT

!========================================================================
        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old)
          
          !>  The input is the field nsigma declared in this module. This routine generates a 
          !>  global update with  and returns the propability  
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)  
          !>   
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
        End Subroutine Global_move
!========================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
        end Function Delta_S0_global
!========================================================================
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
