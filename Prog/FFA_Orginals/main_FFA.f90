Program Main

  Use Operator_mod
  Use Lattices_v3 
  Use MyMats 
  Use Hamiltonian
  Use Control
  Use Tau_m_mod
  Use Hop_mod
  Implicit none
#include "machine"
#ifdef MPI
  include 'mpif.h'
#endif   


  Interface
     SUBROUTINE WRAPUL(NTAU1, NTAU, UL ,DL, VL)
       Use Hamiltonian
       Implicit none
       COMPLEX (KIND=8) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
       COMPLEX (KIND=8) :: DL(Ndim,N_FL)
       Integer :: NTAU1, NTAU
     END SUBROUTINE WRAPUL
     SUBROUTINE CGR(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)
       Use UDV_Wrap_mod
       Implicit None
       COMPLEX(Kind=8), Dimension(:,:), Intent(In)    :: URUP, VRUP, ULUP, VLUP
       COMPLEX(Kind=8), Dimension(:), Intent(In)      :: DLUP, DRUP
       COMPLEX(Kind=8), Dimension(:,:), Intent(Inout) :: GRUP
       COMPLEX(Kind=8) :: PHASE
       INTEGER         :: NVAR
     END SUBROUTINE CGR
     SUBROUTINE WRAPGRUP(GR,NTAU,PHASE)
       Use Hamiltonian
       Implicit none
       COMPLEX (Kind=8), INTENT(INOUT) ::  GR(Ndim,Ndim,N_FL)
       COMPLEX (Kind=8), INTENT(INOUT) ::  PHASE
       INTEGER, INTENT(IN) :: NTAU
     END SUBROUTINE WRAPGRUP
     SUBROUTINE WRAPGRDO(GR,NTAU,PHASE)
       Use Hamiltonian 
       Implicit None
       COMPLEX (Kind=8), INTENT(INOUT) :: GR(NDIM,NDIM,N_FL)
       COMPLEX (Kind=8), INTENT(INOUT) :: PHASE
       Integer :: NTAU
     end SUBROUTINE WRAPGRDO
     SUBROUTINE WRAPUR(NTAU, NTAU1, UR, DR, VR)
       Use Hamiltonian
       Use UDV_Wrap_mod
       Implicit None
       COMPLEX (KIND=8) :: UR(Ndim,Ndim,N_FL), VR(Ndim,Ndim,N_FL)
       COMPLEX (KIND=8) :: DR(Ndim,N_FL)
       Integer :: NTAU1, NTAU
     END SUBROUTINE WRAPUR
     
  end Interface

  COMPLEX (Kind=8), Dimension(:)    , Allocatable   ::  D
  COMPLEX (KIND=8), Dimension(:,:)  , Allocatable   ::  TEST, A, U, V

  COMPLEX (Kind=8), Dimension(:,:)  , Allocatable    :: DL, DR
  COMPLEX (Kind=8), Dimension(:,:,:), Allocatable    :: UL, VL, UR, VR
  COMPLEX (Kind=8), Dimension(:,:,:), Allocatable    :: GR
  

  Integer :: Nwrap, NSweep, NBin, Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
  Integer :: NTAU, NTAU1

  NAMELIST /VAR_QMC/   Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST

  Integer :: Ierr, I,J,nf, nst, n
  Complex (Kind=8) :: Z_ONE = cmplx(1.d0,0.d0), Phase, Z, Z1
  
  ! Space for storage.
  COMPLEX (Kind=8), Dimension(:,:,:)  , Allocatable :: DST
  COMPLEX (Kind=8), Dimension(:,:,:,:), Allocatable :: UST,  VST

  ! For tests
  Integer, external :: nranf
  Real (kind=8) :: Weight
  Integer :: nr,nth
  Logical :: Log
#ifdef MPI
  Integer        ::  Isize, Irank
  INTEGER        :: STATUS(MPI_STATUS_SIZE)
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
  
  ! Write(6,*) 'Call Ham_set'
  Call Ham_set
  ! Write(6,*) 'End Call Ham_set'
  Call confin 
  Call Hop_mod_init
  !Call Hop_mod_test 
  !stop

#ifdef MPI
  If ( Irank == 0 ) then 
#endif
     OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(*,*) 'unable to open <parameters>',ierr
        STOP
     END IF
     READ(5,NML=VAR_QMC)
     CLOSE(5)
#ifdef MPI
  Endif
  CALL MPI_BCAST(Nwrap   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NSweep  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NBin    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Ltau    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(LOBS_EN ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(LOBS_ST ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
  
  
  Call control_init
  Call Alloc_obs(Ltau)
  Call Op_SetHS


!!$#ifdef Ising_test
!!$  ! Test Ising
!!$  DO  NBC = 1, NBIN
!!$     Call Init_obs
!!$     DO NSW = 1, NSWEEP
!!$        do nth = 1,Ltrot*2*Latt%N
!!$           Nt = nranf(Ltrot)
!!$           Nr = nranf(2*Latt%N)
!!$           Weight = S0(nr,nt)
!!$           log =.false.
!!$           if (Weight > ranf()) then 
!!$              nsigma(nr,nt) = - nsigma(nr,nt)
!!$              log =.true.
!!$           endif
!!$           Call Control_upgrade(log)
!!$        enddo
!!$        Call Obser
!!$     Enddo
!!$     Call Preq
!!$  Enddo
!!$  Call Ham_confout
!!$  Call control_Print
!!$  Stop
!!$  ! End Test Ising
!!$#endif

  Allocate( DL(NDIM,N_FL), DR(NDIM,N_FL) )
  Allocate( UL(NDIM,NDIM,N_FL), VL(NDIM,NDIM,N_FL), &
       &    UR(NDIM,NDIM,N_FL), VR(NDIM,NDIM,N_FL), GR(NDIM,NDIM,N_FL ) )
  NSTM = LTROT/NWRAP
#ifdef MPI
  if ( Irank == 0 ) then
#endif
     Open (Unit = 50,file="info",status="unknown",position="append")
     Write(50,*) 'Sweeps             : ', Nsweep
     Write(50,*) 'Bin                : ', NBin
     Write(50,*) 'Measure Int.       : ', LOBS_ST, LOBS_EN
     Write(50,*) 'Stabilization,Wrap : ', Nwrap
     Write(50,*) 'Nstm               : ', NSTM
     Write(50,*) 'Ltau               : ', Ltau 
     close(50)
#ifdef MPI
  endif
#endif
  
  Allocate ( UST(NDIM,NDIM,NSTM,N_FL), VST(NDIM,NDIM,NSTM,N_FL), DST(NDIM,NSTM,N_FL) )
  Allocate ( Test(Ndim,Ndim) )

  NST = NINT( DBLE(LTROT)/DBLE(NWRAP) )
  !Write(6,*) "Write UL ", NST
  Do nf = 1,N_FL
     CALL INITD(UL(:,:,Nf),Z_ONE)
     do I = 1,Ndim
        DL(I,Nf) = Z_ONE
     enddo
     CALL INITD(VL(:,:,nf),Z_ONE)
     DO I = 1,NDim
        DO J = 1,NDim
           UST(I,J,NST,nf) = UL(I,J,nf)
           VST(I,J,NST,nf) = VL(I,J,nf)
        ENDDO
     ENDDO
     DO I = 1,NDim
        DST(I,NST,nf) = DL(I,nf)
     ENDDO

     CALL INITD(UR(:,:,nf),Z_ONE)
     CALL INITD(VR(:,:,nf),Z_ONE)
     Do I = 1,Ndim
        DR(I,nf) = Z_ONE
     Enddo
  Enddo
  
  DO NT = LTROT-NWRAP,NWRAP,-1
     IF ( MOD(NT,NWRAP)  == 0 ) THEN
        NT1 = NT + NWRAP
        !Write(6,*) 'Calling Wrapul:', NT1,NT 
        CALL WRAPUL(NT1,NT,UL,DL, VL)
        NST = NINT( DBLE(NT)/DBLE(NWRAP) )
        !Write(6,*) "Write UL ", NST
        Do nf = 1,N_FL
           DO I = 1,Ndim
              DO J = 1,Ndim
                 UST(I,J,NST,nf) = UL(I,J,nf)
                 VST(I,J,NST,nf) = VL(I,J,nf)
              ENDDO
           ENDDO
           DO I = 1,Ndim
              DST(I,NST,nf) = DL(I,nf)
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  CALL WRAPUL(NWRAP,0, UL ,DL, VL)

  !WRITE(6,*) 'Filling up storage'
  !Write(6,*) 'Done wrapping'
  NVAR = 1
  Phase = cmplx(1.d0,0.d0)
  do nf = 1,N_Fl
     CALL CGR(Z, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
     Phase = Phase*Z
  Enddo
  call Op_phase(Phase,OP_V,Nsigma,N_SUN) 
#ifdef MPI 
  WRITE(6,*) 'Phase is: ', Irank, PHASE, GR(1,1,1)
#else
  WRITE(6,*) 'Phase is: ',  PHASE
!!$  if (N_FL == 1) then
!!$     Do n = 1,Ndim
!!$        Write(6,*) GR(1,n,1)
!!$     enddo
!!$  else
!!$     Do n = 1,Ndim
!!$        Write(6,*) GR(1,n,1), GR(1,n,2)
!!$     enddo
!!$  endif
#endif

  Call Control_init

  DO  NBC = 1, NBIN
     ! Here, you have the green functions on time slice 1.
     ! Set bin observables to zero.

     Call Init_obs(Ltau)
     DO NSW = 1, NSWEEP
 
        !Propagation from 1 to Ltrot
        !Set the right storage to 1
        
        do nf = 1,N_FL
           CALL INITD(UR(:,:,nf),Z_ONE)
           CALL INITD(VR(:,:,nf),Z_ONE)
           do n = 1,Ndim
              DR(n,nf)= Z_ONE
           Enddo
        Enddo
        
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           !Write(6,*) "Hi"
           CALL WRAPGRUP(GR,NTAU,PHASE) 
           !Write(6,*) "Hi1"
           IF ( MOD(NTAU1,NWRAP ) .EQ. 0 ) THEN
              NST = NINT( DBLE(NTAU1)/DBLE(NWRAP) )
              NT1 = NTAU1 - NWRAP
              CALL WRAPUR(NT1, NTAU1,UR, DR, VR)
              Z = cmplx(1.d0,0.d0)
              Do nf = 1, N_FL
                 DO J = 1,Ndim
                    DO I = 1,Ndim
                       UL(I,J,nf) = UST(I,J,NST,nf)
                       VL(I,J,nf) = VST(I,J,NST,nf)
                    ENDDO
                 ENDDO
                 DO I = 1,Ndim
                    DL(I,nf) = DST(I,NST,nf)
                 ENDDO
                 ! Write in store Right prop from 1 to LTROT/NWRAP
                 !Write(6,*)  'Write UR, read UL ', NTAU1, NST
                 DO J = 1,Ndim
                    DO I = 1,Ndim
                       UST(I,J,NST,nf) = UR(I,J,nf)
                       VST(I,J,NST,nf) = VR(I,J,nf)
                    ENDDO
                 ENDDO
                 DO I = 1,Ndim
                    DST(I,NST,nf) = DR(I,nf)
                 ENDDO
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 !Write(6,*) ' Call Cgr'
                 do J = 1,Ndim
                    do I = 1,Ndim
                       TEST(I,J) = GR(I,J,nf)
                    enddo
                 enddo
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
                 Z = Z*Z1
                 !Write(6,*) 'Calling control ',NTAU1, Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
           ENDIF

           IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
              !Write(6,*) 'Call obser ', Ntau1
              CALL Obser( GR, PHASE, Ntau1 )
              !Write(6,*) 'Return obser'
           ENDIF
           !Write(6,*) NTAU1
        ENDDO
        
        Do nf = 1,N_FL
           CALL INITD(UL(:,:,nf),Z_ONE)
           CALL INITD(VL(:,:,nf),Z_ONE)
           Do n = 1,Ndim
              DL(n,nf)  = Z_ONE
           Enddo
        ENDDO

        DO NTAU = LTROT,1,-1
           NTAU1 = NTAU - 1
           CALL WRAPGRDO(GR,NTAU, PHASE)
           IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
              CALL Obser( GR, PHASE, Ntau1 )
           ENDIF
           IF ( MOD(NTAU1,NWRAP).EQ.0 .AND. NTAU1.NE.0 ) THEN
              ! WRITE(50,*) 'Recalc at  :', NTAU1
              NST = NINT( DBLE(NTAU1)/DBLE(NWRAP) )
              NT1 = NTAU1 + NWRAP
              !Write(6,*) 'Wrapul : ', NT1, NTAU1
              CALL WRAPUL(NT1,NTAU1, UL, DL, VL )
              !Write(6,*)  'Write UL, read UR ', NTAU1, NST
              Z = cmplx(1.d0,0.d0)
              do nf = 1,N_FL
                 DO J = 1,Ndim
                    DO I = 1,Ndim
                       UR(I,J,nf) = UST(I,J,NST,nf)
                       VR(I,J,nf) = VST(I,J,NST,nf)
                    ENDDO
                 ENDDO
                 DO I = 1,Ndim
                    DR(I,nf) = DST(I,NST,nf)
                 ENDDO
                 ! WRITE in store the left prop. from LTROT/NWRAP-1 to 1
                 DO J = 1,Ndim
                    DO I = 1,Ndim
                       UST(I,J,NST,nf) =  UL(I,J,nf)
                       VST(I,J,NST,nf) =  VL(I,J,nf)
                    ENDDO
                 ENDDO
                 DO I = 1,Ndim
                    DST(I,NST,nf) =  DL(I,nf)
                 ENDDO
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 !Write(6,*) ' Call Cgr'
                 do J = 1,Ndim
                    do I = 1,Ndim
                       TEST(I,J) = GR(I,J,nf)
                    enddo
                 enddo
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
                 Z = Z*Z1
                 !Write(6,*) 'Calling control: ', NTAU1, Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
           ENDIF
        ENDDO

        !Calculate and compare green functions on time slice 0.
        NT1 = 0
        CALL WRAPUL(NWRAP,NT1, UL, DL, VL )
        
        do nf = 1,N_FL
           CALL INITD(UR(:,:,nf),Z_ONE)
           CALL INITD(VR(:,:,nf),Z_ONE)
           DO I = 1,Ndim
              DR(I,nf) = Z_ONE
           ENDDO
        ENDDO
        Z = cmplx(1.d0,0.d0)
        do nf = 1,N_FL
           do J = 1,Ndim
              do I = 1,Ndim
                 TEST(I,J) = GR(I,J,nf)
              enddo
           enddo
           NVAR = 1
           CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
           Z = Z*Z1
           !Write(6,*) 'Calling control  0', Z1
           Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
        ENDDO
        call Op_phase(Z,OP_V,Nsigma,N_SUN) 
        Call Control_PrecisionP(Z,Phase)
        Phase = Z
        NST =  NINT( DBLE(LTROT)/DBLE(NWRAP) )
        Do nf = 1,N_FL
           DO I = 1,Ndim
              DO J = 1,Ndim
                 UST(I,J,NST,nf) = CMPLX(0.D0,0.D0)
                 VST(I,J,NST,nf) = CMPLX(0.D0,0.D0)
              ENDDO
           ENDDO
           DO I = 1,Ndim
              DST(I  ,NST,nf)  = CMPLX(1.D0,0.D0)
              UST(I,I,NST,nf)  = CMPLX(1.D0,0.D0)
              VST(I,I,NST,nf)  = CMPLX(1.D0,0.D0)
           ENDDO
        enddo
        IF ( LTAU == 1 ) then
!!$#ifdef MPI
!!$           Write(6,*) Irank, 'Calling Tau_m', NWRAP, NSTM
!!$#else
!!$           Write(6,*)  'Calling Tau_m', NWRAP, NSTM
!!$#endif

           Call TAU_M( UST,DST,VST, GR, PHASE, NSTM, NWRAP ) 
!!$#ifdef MPI
!!$           Write(6,*) Irank, 'Back Calling Tau_m'
!!$#else
!!$           Write(6,*)  'Back Calling Tau_m'
!!$#endif
        endif
           
     ENDDO
     Call Pr_obs(Ltau)
     Call confout
  Enddo
  Call Control_Print

#ifdef MPI
   CALL MPI_FINALIZE(ierr)
#endif

end Program Main
