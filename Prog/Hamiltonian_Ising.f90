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
      !Complex (Kind=Kind(0.d0)), dimension(:,:,:), allocatable :: Exp_T(:,:,:), Exp_T_M1(:,:,:)


      
      ! What is below is  private 
      
      Type (Lattice),       private :: Latt
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_xi, ham_h, ham_J, ham_U, Ham_Vint, Ham_chem
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta
      Character (len=64),   private :: Model, Lattice_type
      Integer, allocatable, private :: L_bond(:,:), Ising_nnlist(:,:)
      Real (Kind=Kind(0.d0)),        private :: DW_Ising_tau  (-1:1), DW_Ising_Space(-1:1)
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord 
      Real (Kind=Kind(0.d0)),        private :: Bound


      ! Observables
      Integer,                       private :: Nobs, Norb
      Complex (Kind=Kind(0.d0)), allocatable, private :: obs_scal(:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Ising_cor(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green_eq (:,:,:), Spin_eq(:,:,:), Den_eq   (:,:,:) 
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green_eq0 (:), Spin_eq0(:), Pair_eq0(:), Den_eq0(:) 
      Complex (Kind=Kind(0.d0)), allocatable, private :: Ising_cor0(:)

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


          NAMELIST /VAR_Ising/    ham_T, ham_chem, ham_xi, ham_h, ham_J, Beta, dtau, N_SUN


#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
          
          
#ifdef MPI 
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
          ENDIF
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          If ( Model == "Ising" )  then
             N_FL = 1
#ifdef MPI
             If (Irank == 0 ) then
#endif
                OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
                READ(5,NML=VAR_Ising)
                CLOSE(5)
#ifdef MPI
             endif
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_xi   ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_h    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_J    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(N_SUN    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
          else 
             Write(6,*) ' Model not yet programmed : '
             Stop
          endif
          Call Ham_hop
          Ltrot = nint(beta/dtau)
#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'N_SUN         : ', N_SUN
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'xi            : ', Ham_xi
             Write(50,*) 'h             : ', Ham_h
             Write(50,*) 'Ham_J         : ', Ham_J
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
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          If ( Lattice_type =="Square" ) then
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             Ndim = Latt%N
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
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif
        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !Setup the hopping
          
          Integer :: I, I1, I2, n, Ncheck,nc

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,Ncheck
                Call Op_make(Op_T(nc,n),Ndim)
                If (One_dimensional ) then 
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I,1,0)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem,0.d0, kind(0.D0))
                   ENDDO
                else
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I,1,0)
                      I2 = Latt%nnlist(I,0,1)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,   0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,   0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,   0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,   0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem,0.d0, kind(0.D0))
                   ENDDO
                endif
                
                Do I = 1,Latt%N
                   Op_T(nc,n)%P(i) = i 
                Enddo
                if ( abs(Ham_T) < 1.E-6 .and.  abs(Ham_chem) < 1.E-6 ) then 
                   Op_T(nc,n)%g=cmplx(0.d0 ,0.d0, kind(0.D0))
                else
                   Op_T(nc,n)%g=cmplx(-Dtau,0.d0, kind(0.D0)) 
                endif
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
          
          Integer :: nf, nth, n, n1, n2, n3, n4, I, I1, I2, nc
          Real (Kind=Kind(0.d0)) :: X_p(2)


          If (Model == "Ising"  ) then
             Allocate( Op_V(N_coord*Latt%N,N_FL) )
             do nf = 1,N_FL
                do i  = 1, N_coord*Latt%N
                   Call Op_make(Op_V(i,nf),2)
                enddo
             enddo
             Allocate (L_Bond(Latt%N,2))
             Do nf = 1,N_FL
                nc = 0
                do nth = 1,2*N_coord
                   Do n1= 1, L1/2
                      Do n2 = 1,L2
                         nc = nc + 1
                         If (nth == 1 ) then
                            X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p 
                            I1 = Inv_R(X_p,Latt)
                            I2 = Latt%nnlist(I1,1,0)
                            L_bond(I1,1) = nc
                         elseif (nth == 2) then
                            X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p  + latt%a1_p
                            I1 = Inv_R(X_p,Latt)
                            I2 = Latt%nnlist(I1,1,0)
                            L_bond(I1,1) = nc
                         elseif (nth == 3) then
                            X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p 
                            I1 = Inv_R(X_p,Latt)
                            I2 = Latt%nnlist(I1,0,1)
                            L_bond(I1,2) = nc
                         elseif (nth == 4) then
                            X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p  + latt%a2_p
                            I1 = Inv_R(X_p,Latt)
                            I2 = Latt%nnlist(I1,0,1)
                            L_bond(I1,2) = nc
                         endif
                         Op_V(nc,nf)%P(1) = I1; Op_V(nc,nf)%P(2) = I2   
                         Op_V(nc,nf)%O(1,2) = cmplx(1.d0  ,0.d0, kind(0.D0))
                         Op_V(nc,nf)%O(2,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                         Op_V(nc,nf)%g = cmplx(-dtau*Ham_xi,0.d0, kind(0.D0))
                         Op_V(nc,nf)%alpha  = cmplx(0.d0  ,0.d0, kind(0.D0))
                         Op_V(nc,nf)%type   = 1
                         Call Op_set( Op_V(nc,nf) )
                         ! For a single flavour, the operator reads:  
                         ! g*s*( c^{dagger} O c   - alpha ))
                         ! with s the HS field. 
                      Enddo
                   Enddo
                Enddo
             Enddo
!!$             Open (Unit=10,File="Latt",status="unknown")
!!$             Do I = 1,Latt%N
!!$                X_p = dble(latt%list(I,1))*latt%a1_p + dble(latt%list(I,2))*latt%a2_p
!!$                Write(10,*) X_p(1), X_p(2) 
!!$                Write(10,*) X_p(1)+ latt%a1_p(1), X_p(2) + latt%a1_p(2)
!!$                Write(10,*)
!!$                Write(10,*) X_p(1), X_p(2) 
!!$                Write(10,*) X_p(1)+ latt%a2_p(1), X_p(2) + latt%a2_p(2)
!!$                Write(10,*)
!!$             Enddo
!!$             Close(10)
             allocate(Ising_nnlist(2*Latt%N,4)) 
             do I  = 1,Latt%N
                n  = L_bond(I,1)
                n1 = L_bond(Latt%nnlist(I, 1, 0),2)
                n2 = L_bond(Latt%nnlist(I, 0, 0),2)
                n3 = L_bond(Latt%nnlist(I, 0,-1),2)
                n4 = L_bond(Latt%nnlist(I, 1,-1),2)
                Ising_nnlist(n,1) = n1
                Ising_nnlist(n,2) = n2
                Ising_nnlist(n,3) = n3
                Ising_nnlist(n,4) = n4
                n  = L_bond(I,2)
                n1 = L_bond(Latt%nnlist(I, 0, 1),1)
                n2 = L_bond(Latt%nnlist(I,-1, 1),1)
                n3 = L_bond(Latt%nnlist(I,-1, 0),1)
                n4 = L_bond(Latt%nnlist(I, 0, 0),1)
                Ising_nnlist(n,1) = n1
                Ising_nnlist(n,2) = n2
                Ising_nnlist(n,3) = n3
                Ising_nnlist(n,4) = n4
             enddo
             DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
             DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
             DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J) 
             DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J) 
!!$             Open (Unit=10,File="Ising_latt",status="unknown")
!!$             nf = 1
!!$             Do I = 1,Latt%N
!!$                I1 = Op_V(L_bond(I,1),nf)%P(2)
!!$                I2 = Op_V(L_bond(I,2),nf)%P(2)
!!$                X_p = dble(latt%list(I,1))*latt%a1_p + dble(latt%list(I,2))*latt%a2_p
!!$                X1_p = dble(latt%list(I1,1))*latt%a1_p + dble(latt%list(I1,2))*latt%a2_p
!!$                X2_p = dble(latt%list(I2,1))*latt%a1_p + dble(latt%list(I2,2))*latt%a2_p
!!$                Write(10,*) X_p (1), X_p (2) 
!!$                Write(10,*) X1_p(1), X1_p(2)
!!$                Write(10,*)
!!$                Write(10,*) X_p (1), X_p (2) 
!!$                Write(10,*) X2_p(1), X2_p(2)
!!$                Write(10,*)
!!$             Enddo
!!$             Close(10)
          endif
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt 
          Integer :: i, nt1 
          S0 = 1.d0
          If (Model == "Ising" ) then 
             do i = 1,4
                S0 = S0*DW_Ising_space(nsigma(n,nt)*nsigma(Ising_nnlist(n,i),nt))
             enddo
             nt1 = nt +1 
             if (nt1 > Ltrot) nt1 = 1
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             nt1 = nt - 1 
             if (nt1 < 1  ) nt1 = Ltrot
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
          endif
        end function S0

!===================================================================================           

        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau

          Norb = 2
          Allocate ( Obs_scal(5) )
          Allocate ( Ising_cor (Latt%N,Norb,Norb) )
          Allocate ( Green_eq(Latt%N,1,1), Spin_eq(Latt%N,1,1), Den_eq(Latt%N,1,1)  )
          Allocate ( Ising_cor0(Norb), Green_eq0(1), Spin_eq0(1), Den_eq0(1) )
          If (Ltau == 1) then 
             Allocate ( Green_tau(Latt%N,Ltrot+1,1,1), Den_tau(Latt%N,Ltrot+1,1,1) )
          endif


          
        end Subroutine Alloc_obs

!===================================================================================           
        
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          Nobs = 0
          Obs_scal  = cmplx(0.d0, 0.d0, kind(0.D0))
          Ising_cor = cmplx(0.d0, 0.d0, kind(0.D0))
          Green_eq  = cmplx(0.d0, 0.d0, kind(0.D0))
          Spin_eq   = cmplx(0.d0, 0.d0, kind(0.D0))
          Den_eq    = cmplx(0.d0, 0.d0, kind(0.D0))
          Ising_cor0= cmplx(0.d0, 0.d0, kind(0.D0))
          Green_eq0 = cmplx(0.d0, 0.d0, kind(0.D0))
          Spin_eq0  = cmplx(0.d0, 0.d0, kind(0.D0))
          Den_eq0   = cmplx(0.d0, 0.d0, kind(0.D0))

          If (Ltau == 1) then
             NobsT = 0
             Phase_tau = cmplx(0.d0, 0.d0, kind(0.D0))
             Green_tau = cmplx(0.d0, 0.d0, kind(0.D0))
             Den_tau = cmplx(0.d0, 0.d0, kind(0.D0))
          endif

        end Subroutine Init_obs
        
!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.D0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.D0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.D0)) :: GRC(Ndim,Ndim,N_FL), ZK, myZ
          Complex (Kind=Kind(0.D0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZV(4)
          Integer :: I,J, no,no1, n, n1, imj, nf, myi, imjx(4), temp
          
          Nobs = Nobs + 1
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          

          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I,J,nf)  = - GR(J,I,nf)
                Enddo
                GRC(I, I, nf) = 1.d0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >
          ! Compute scalar observables. 
          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do J = 1,Ndim
                Zkin = Zkin + sum(Op_T(1,nf)%O(:, j)*Grc(:, j, nf))
             ENddo
          Enddo
          Zkin = Zkin * dble(N_SUN)

          Zrho = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)

          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))

          myZ = ZP * ZS
          Obs_scal(1) = Obs_scal(1) + zrho * myZ
          Obs_scal(2) = Obs_scal(2) + zkin * myZ
          Obs_scal(3) = Obs_scal(3) + Zpot * myZ
          Obs_scal(4) = Obs_scal(4) + (zkin +  Zpot)*myZ
          Obs_scal(5) = Obs_scal(5) + ZS
          ! You will have to allocate more space if you want to include more  scalar observables.
          
          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
! Latt%N is according to Fakher divisible by 4 for the Ising model
             Do I = 1, Latt%N, 4
                Do J = 1, Latt%N
!                   imj = latt%imj(I,J)
                   imjx = latt%imj(I: I + 3,J)
                   GREEN_EQ(imjx,1,1) = GREEN_EQ(imjx,1,1) + Z * GRC(I:I+3,J,1) *  myZ
                   SPIN_Eq (imjx,1,1) = SPIN_Eq (imjx,1,1) + Z * GRC(I:I+3,J,1) * GR(I:I+3,J,1) * myZ
                   do myi = 1,4
                    DEN_Eq  (imjx(myi),1,1) = DEN_Eq  (imjx(myi),1,1) + (GRC(I-1+myi,I-1+myi,1) * GRC(J,J,1) *Z + &
                        &         GRC(I-1+myi,J,1) *  GR(I-1+myi,J,1)) * Z* myZ
                   enddo
                Den_eq0(1) = Den_eq0(1) + Z* myZ *(GRC(I,I,1)+GRC(I+1,I+1,1)+GRC(I+2,I+2,1)+GRC(I+3,I+3,1))
                do no = 1,Norb
                    ZV = dble(nsigma(L_bond(I : I + 3,no), ntau))*myZ
!                    n = L_bond(I,no)
                      do no1 = 1,Norb
                         myi = nsigma(L_bond(J,no1),ntau)
                         Ising_cor(imjx,no,no1) = Ising_cor(imjx,no,no1) + CMPLX(myi*dble(ZV), myi*aimag(ZV),kind(0.D0))
                      enddo
                   enddo
                enddo
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
    
          Phase_bin = Obs_scal(5)/dble(Nobs)
          File_pr ="Ising_eq"
          Call Print_bin(Ising_cor,Ising_cor0,Latt, Nobs, Phase_bin, file_pr)
          File_pr ="Green_eq"
          Call Print_bin(Green_eq, Green_eq0, Latt, Nobs, Phase_bin, file_pr)
          File_pr ="Spin_eq"
          Call Print_bin(Spin_eq,  Spin_eq0,  Latt, Nobs, Phase_bin, file_pr)
          File_pr ="Den_eq"
          Call Print_bin(Den_eq,   Den_eq0,   Latt, Nobs, Phase_bin, file_pr)


          File_pr ="ener"
          Call Print_scal(Obs_scal, Nobs, file_pr)
          If (Ltau == 1) then
             Phase_tau = Phase_tau/dble(NobsT)
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

          ZP = PHASE/Real(Phase,Kind=Kind(0.d0))
          ZS = cmplx(Real(Phase,Kind=Kind(0.d0))/Abs(Real(Phase,Kind=Kind(0.d0))), 0.d0, kind(0.D0))
          If (NT == 0 ) then 
             Phase_tau = Phase_tau + ZS
             NobsT     = NobsT + 1
          endif
          If ( N_FL == 1 ) then 
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             Do I = 1,Latt%N
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   Green_tau(imj,nt+1,1,1) = green_tau(imj,nt+1,1,1)  +  Z * GT0(I,J,1) * ZP* ZS
                   Den_tau  (imj,nt+1,1,1) = Den_tau  (imj,nt+1,1,1)  -  Z * GT0(I,J,1)*G0T(J,I,1) * ZP* ZS
                Enddo
             Enddo
          Endif
        end Subroutine OBSERT


    end Module Hamiltonian
