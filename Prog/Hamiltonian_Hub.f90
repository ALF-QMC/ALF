    Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      !, only: Op_make, Op_set, Operator
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Print_bin_mod
      
      Implicit none

      !Private 
      !Public  Op_V, Op_T, nsigma, Ndim,  N_FL,  N_SUN,  Ltrot
 
      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Type (WaveFunction), dimension(:),   allocatable  :: WF_L
      Type (WaveFunction), dimension(:),   allocatable  :: WF_R
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot, Thtrot
      Logical              :: Projector
!>    Variables for updating scheme
      Logical              :: Propose_S0, Global_moves
      Integer              :: N_Global

      !Complex (Kind=Kind(0.d0)), dimension(:,:,:), allocatable :: Exp_T(:,:,:), Exp_T_M1(:,:,:)
      
      ! ToDo.  Public and private subroutines. 

      ! What is below is  private 
      Type (Lattice),       private :: Latt
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord 


      ! Observables
      Integer,                       private :: Nobs, Norb
      Complex (Kind=Kind(0.d0)), allocatable, private :: obs_scal(:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green_eq  (:,:,:), SpinZ_eq (:,:,:), SpinXY_eq (:,:,:), &
           &                                    Den_eq(:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green_eq0 (:), SpinZ_eq0(:), SpinXY_eq0(:), &
           &                                    Den_eq0(:)
      
      ! For time displaced
      Integer,                       private :: NobsT
      Complex (Kind=Kind(0.d0)),              private :: Phase_tau
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green_tau(:,:,:,:), Den_tau(:,:,:,:), SpinZ_tau(:,:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: SpinXY_tau(:,:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Green0_tau(:), Den0_tau(:), SpinZ0_tau(:), SpinXY0_tau(:)
      
    contains 


      Subroutine Ham_Set
#ifdef MPI
          Use mpi
#endif
          Implicit none


          integer :: ierr

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U, Dtau, Beta


#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
          !          NAMELIST /VAR_Model/  N_FL,  N_SUN,  ham_T , ham_xi, ham_h, ham_J,  ham_U, Ham_Vint, &
          !               &         Dtau, Beta


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
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          Propose_S0 = .false.
          Global_moves =.false.
          N_Global = 1

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
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_Hubbard)
             CLOSE(5)
#ifdef MPI
          endif
          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
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
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'N_SUN         : ', N_SUN
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_U         : ', Ham_U
             Write(50,*) 'Ham_chem      : ', Ham_chem
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
          !Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, I2, n, Ncheck,nc

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,Ncheck
                Call Op_make(Op_T(nc,n),Ndim)
                If (One_dimensional ) then 
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I, 1, 0)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                   ENDDO
                else
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I,1,0)
                      I2 = Latt%nnlist(I,0,1)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
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
          
          Integer :: nf, I, nc
          Real (Kind=Kind(0.d0)) :: X

          
          If (Model == "Hubbard_SU2")  then
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Latt%N,N_FL))
             do nf = 1,N_FL
                do i  = 1, Latt%N
                   Call Op_make(Op_V(i,nf),1)
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
             Enddo
          Elseif (Model == "Hubbard_Mz")  then
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Latt%N,N_FL))
             do nf = 1,N_FL
                do i  = 1, Latt%N
                   Call Op_make(Op_V(i,nf),1)
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
             Enddo
          Endif
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
         
          Allocate ( Obs_scal(5) )
          Allocate ( Green_eq(Latt%N,1,1), SpinZ_eq(Latt%N,1,1), SpinXY_eq(Latt%N,1,1), &
               &     Den_eq(Latt%N,1,1) )
          Allocate ( Green_eq0(1), SpinZ_eq0(1), SpinXY_eq0(1), Den_eq0(1) )

               
          If (Ltau == 1) then 
             Allocate ( Green_tau(Latt%N,Ltrot+1,1,1), Den_tau   (Latt%N,Ltrot+1,1,1) )
             Allocate ( SpinZ_tau(Latt%N,Ltrot+1,1,1), SpinXY_tau(Latt%N,Ltrot+1,1,1) )
             Allocate ( Green0_tau(1), Den0_tau(1), SpinZ0_tau(1), SpinXY0_tau(1) )
          endif
          
        end Subroutine Alloc_obs

!===================================================================================           
        
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          Nobs = 0
          Obs_scal  = cmplx(0.d0,0.d0, kind(0.D0))
          Green_eq  = cmplx(0.d0,0.d0, kind(0.D0))
          SpinZ_eq  = cmplx(0.d0,0.d0, kind(0.D0))
          SpinXY_eq = cmplx(0.d0,0.d0, kind(0.D0))
          Den_eq    = cmplx(0.d0,0.d0, kind(0.D0))
          Green_eq0 = cmplx(0.d0,0.d0, kind(0.D0))
          SpinZ_eq0 = cmplx(0.d0,0.d0, kind(0.D0))
          SpinXY_eq0= cmplx(0.d0,0.d0, kind(0.D0))
          Den_eq0   = cmplx(0.d0,0.d0, kind(0.D0))

          If (Ltau == 1) then
             NobsT = 0
             Phase_tau  = cmplx(0.d0,0.d0, kind(0.D0))
             Green_tau  = cmplx(0.d0,0.d0, kind(0.D0))
             Den_tau    = cmplx(0.d0,0.d0, kind(0.D0))
             SpinZ_tau  = cmplx(0.d0,0.d0, kind(0.D0))
             SpinXY_tau = cmplx(0.d0,0.d0, kind(0.D0))
             Green0_tau = cmplx(0.d0,0.d0, kind(0.D0))
             Den0_tau   = cmplx(0.d0,0.d0, kind(0.D0))
             SpinZ0_tau  = cmplx(0.d0,0.d0, kind(0.D0))
             SpinXY0_tau = cmplx(0.d0,0.d0, kind(0.D0))
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
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS
          Integer :: I,J, imj, nf, dec
          
          Nobs = Nobs + 1
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
          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do J = 1,Ndim
                Zkin = Zkin + sum(Op_T(1,nf)%O(:, j)*Grc(:, j, nf))
             ENddo
          Enddo
          Zkin = Zkin * dble(N_SUN)

          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)

          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          If ( Model == "Hubbard_SU2" ) then
            dec = 1
          else
            dec = 2
          endif
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, dec)
          Enddo
          Zpot = Zpot*ham_U

          Obs_scal(1) = Obs_scal(1) + zrho * ZP*ZS
          Obs_scal(2) = Obs_scal(2) + zkin * ZP*ZS
          Obs_scal(3) = Obs_scal(3) + Zpot * ZP*ZS
          Obs_scal(4) = Obs_scal(4) + (zkin +  Zpot)*ZP*ZS
          Obs_scal(5) = Obs_scal(5) + ZS
          ! You will have to allocate more space if you want to include more  scalar observables.
          
          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          If ( Model == "Hubbard_SU2"  ) then 
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             Do I = 1,Latt%N
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   GREEN_EQ  (imj,1,1) = GREEN_EQ  (imj,1,1)  +  Z * GRC(I,J,1) *  ZP*ZS 
                   SPINXY_Eq (imj,1,1) = SPINXY_Eq (imj,1,1)  +  Z * GRC(I,J,1) * GR(I,J,1) * ZP*ZS
                   SPINZ_Eq  (imj,1,1) = SPINZ_Eq  (imj,1,1)  +  Z * GRC(I,J,1) * GR(I,J,1) * ZP*ZS
                   DEN_Eq    (imj,1,1) = DEN_Eq    (imj,1,1)  +         (      &
                        &         GRC(I,I,1) * GRC(J,J,1) *Z     + &
                        &         GRC(I,J,1) *  GR(I,J,1)          &
                        &                                   ) * Z* ZP*ZS
                ENDDO
                Den_eq0(1) = Den_eq0(1) +  Z * GRC(I,I,1) * ZP * ZS
             ENDDO
          elseif (Model == "Hubbard_Mz" ) Then
             DO I = 1,Latt%N
                DO J = 1, Latt%N
                   imj = latt%imj(I,J)
                   SPINZ_Eq (imj,1,1) = SPINZ_Eq (imj,1,1)  +  &
                        & (   GRC(I,J,1) * GR(I,J,1) +  GRC(I,J,2) * GR(I,J,2)    + &
                        &   (GRC(I,I,2) - GRC(I,I,1))*(GRC(J,J,2) - GRC(J,J,1))    ) * ZP*ZS
                   ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                   SPINXY_Eq (imj,1,1) = SPINXY_Eq (imj,1,1)  +  &
                     & (   GRC(I,J,1) * GR(I,J,2) +  GRC(I,J,2) * GR(I,J,1)    ) * ZP*ZS

                   DEN_Eq (imj,1,1) = DEN_Eq (imj,1,1)  +  &
                        & (   GRC(I,J,1) * GR(I,J,1) +  GRC(I,J,2) * GR(I,J,2)    + &
                        &   (GRC(I,I,2) + GRC(I,I,1))*(GRC(J,J,2) + GRC(J,J,1))    ) * ZP*ZS
                enddo
                Den_eq0(1) = Den_eq0(1) + (GRC(I,I,2) + GRC(I,I,1)) * ZP*ZS
             enddo
          Endif
                

        end Subroutine Obser
!==========================================================        
        Subroutine  Pr_obs(LTAU)
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
    
          Phase_bin = Obs_scal(5)/dble(Nobs)
          File_pr ="SpinZ_eq"
          Call Print_bin(SpinZ_eq ,SpinZ_eq0, Latt, Nobs, Phase_bin, file_pr)
          File_pr ="SpinXY_eq"
          Call Print_bin(Spinxy_eq, Spinxy_eq0,Latt, Nobs, Phase_bin, file_pr)
          File_pr ="Den_eq"
          Call Print_bin(Den_eq, Den_eq0, Latt, Nobs, Phase_bin, file_pr)
          File_pr ="Green_eq"
          Call Print_bin(Green_eq, Green_eq0, Latt, Nobs, Phase_bin, file_pr)

          File_pr ="ener"
          Call Print_scal(Obs_scal, Nobs, file_pr)

          If (Ltau == 1) then
             Phase_tau = Phase_tau/dble(NobsT)
             File_pr = "Green_tau"
             Call Print_bin_tau(Green_tau,Latt,NobsT,Phase_tau, file_pr,dtau       )
             File_pr = "Den_tau"
             Call Print_bin_tau(Den_tau,Latt,NobsT,Phase_tau, file_pr,dtau,Den0_tau)
             File_pr = "SpinZ_tau"
             Call Print_bin_tau(SpinZ_tau,Latt,NobsT,Phase_tau, file_pr,dtau,SpinZ0_tau)
             File_pr = "SpinXY_tau"
             Call Print_bin_tau(SpinXY_tau,Latt,NobsT,Phase_tau, file_pr,dtau,SpinXY0_tau)

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

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then 
             Phase_tau = Phase_tau + ZS
             NobsT     = NobsT + 1
          endif
          If ( Model == "Hubbard_SU2"  ) then 

             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I = 1,Latt%N
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   Green_tau(imj,nt+1,1,1) = green_tau(imj,nt+1,1,1)  +  Z * GT0(I,J,1) * ZP* ZS

                   
                   Den_tau  (imj,nt+1,1,1) = Den_tau  (imj,nt+1,1,1)  +                           ( &
                        &    Z*Z*(1.d0 - GTT(I,I,1))*(1.d0 - G00(J,J,1))  -     &
                        &    Z * GT0(I,J,1)*G0T(J,I,1)                                            ) * ZP * ZS


                   SpinZ_tau(imj,nt+1,1,1) = SpinZ_tau(imj,nt+1,1,1)  &
                       &      - Z*G0T(J,I,1) * GT0(I,J,1) *ZP*ZS

                   SpinXY_tau(imj,nt+1,1,1) = SpinXY_tau(imj,nt+1,1,1)  &
                       &      - Z*G0T(J,I,1) * GT0(I,J,1) *ZP*ZS
                Enddo
             Enddo
             
             Do I = 1,Latt%N
                Den0_tau(1) = Den0_tau(1) + Z*(1.d0 - GTT(I,I,1)) * ZP * ZS
             Enddo
          Elseif ( Model == "Hubbard_Mz"  ) then 
             Do I = 1,Latt%N
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   Green_tau(imj,nt+1,1,1) = green_tau(imj,nt+1,1,1)  +   ( GT0(I,J,1) + GT0(I,J,2) ) * ZP* ZS

                   Den_tau(imj,nt+1,1,1)  = Den_tau   (imj,nt+1,1,1)  +  (  &
                        &    (2.D0 - GTT(I,I,1) - GTT(I,I,2) ) * ( 2.D0 - G00(J,J,1) -  G00(J,J,2) ) &
                        & -  ( G0T(J,I,1) * GT0(I,J,1) + G0T(J,I,2) * GT0(I,J,2) )   )*ZP*ZS     

                   SpinZ_tau(imj,nt+1,1,1) = SpinZ_tau(imj,nt+1,1,1) + ( &
                       &    (GTT(I,I,1) -  GTT(I,I,2) ) * ( G00(J,J,1)  -  G00(J,J,2) )   &
                       &  - (G0T(J,I,1) * GT0(I,J,1)  +  G0T(J,I,2) * GT0(I,J,2) )    )*ZP*ZS

                   SpinXY_tau(imj,nt+1,1,1) = SpinXY_tau(imj,nt+1,1,1)  - &
                        &   (G0T(J,I,1) * GT0(I,J,2)  +  G0T(J,I,2) * GT0(I,J,1))*ZP*ZS

                enddo
             enddo
             
             Do I = 1,Latt%N
                Den0_tau(1) = Den0_tau(1) + (2.d0 - GTT(I,I,1) - GTT(I,I,2)) * ZP * ZS
             Enddo
          endif

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
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
        end Function Delta_S0_global
!========================================================================


      end Module Hamiltonian
