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

  COMPLEX (KIND=8), Dimension(:,:)  , Allocatable   ::  TEST

  COMPLEX (Kind=8), Dimension(:,:)  , Allocatable    :: DL, DR
  COMPLEX (Kind=8), Dimension(:,:,:), Allocatable    :: UL, VL, UR, VR
  COMPLEX (Kind=8), Dimension(:,:,:), Allocatable    :: GR
  

  Integer :: Nwrap, NSweep, NBin, NBin_eff,Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
  Integer :: NTAU, NTAU1
  Real(Kind=8) :: CPU_MAX 


  NAMELIST /VAR_QMC/   Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX

  Integer :: Ierr, I,nf, nst, n
  Complex (Kind=8) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Phase, Z, Z1
  Real    (Kind=8) :: ZERO = 10D-8
  Integer, dimension(:), allocatable :: Stab_nt
  
  ! Space for storage.
  COMPLEX (Kind=8), Dimension(:,:,:)  , Allocatable :: DST
  COMPLEX (Kind=8), Dimension(:,:,:,:), Allocatable :: UST,  VST

  ! For tests
  Integer, external :: nranf
  Real (kind=8) :: Weight, Weight_tot
  Integer :: nr,nth, nth1
  Logical :: Log
  
  ! For the truncation of the program:
  Logical        :: prog_truncation
  Real(kind=8)   :: time_bin_start,time_bin_end
    
  
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
     Nwrap=0;  NSweep=0; NBin=0; Ltau=0; LOBS_EN = 0;  LOBS_ST = 0;  CPU_MAX = 0.d0
     OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(*,*) 'unable to open <parameters>',ierr
        STOP
     END IF
     READ(5,NML=VAR_QMC)
     CLOSE(5)
     NBin_eff = NBin
#ifdef MPI
  Endif
  CALL MPI_BCAST(Nwrap   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NSweep  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NBin    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Ltau    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(LOBS_EN ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(LOBS_ST ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(CPU_MAX ,1,MPI_REAL8,  0,MPI_COMM_WORLD,ierr)
#endif
 
  IF (ABS(CPU_MAX) > Zero ) NBIN = 1000000
 
  Call control_init
  Call Alloc_obs(Ltau)
  Call Op_SetHS

  Allocate( DL(NDIM,N_FL), DR(NDIM,N_FL) )
  Allocate( UL(NDIM,NDIM,N_FL), VL(NDIM,NDIM,N_FL), &
       &    UR(NDIM,NDIM,N_FL), VR(NDIM,NDIM,N_FL), GR(NDIM,NDIM,N_FL ) )


  If ( mod(Ltrot,nwrap) == 0  ) then 
     Nstm = Ltrot/nwrap
  else
     nstm = Ltrot/nwrap + 1
  endif
  allocate ( Stab_nt(0:Nstm) )
  Stab_nt(0) = 0
  do n = 1,Nstm -1
     Stab_nt(n) = nwrap*n
  enddo
  Stab_nt(Nstm) = Ltrot
  !do n = 1,Nstm
  !   Write(6,*) n,stab_nt(n)
  !enddo

#ifdef MPI
  if ( Irank == 0 ) then
#endif
     Open (Unit = 50,file="info",status="unknown",position="append")
     Write(50,*) 'Sweeps             : ', Nsweep
     Write(50,*) 'Measure Int.       : ', LOBS_ST, LOBS_EN
     Write(50,*) 'Stabilization,Wrap : ', Nwrap
     Write(50,*) 'Nstm               : ', NSTM
     Write(50,*) 'Ltau               : ', Ltau 
#ifdef MPI
     Write(50,*) 'Number of  threads : ', ISIZE
#endif   
     If ( abs(CPU_MAX) < ZERO ) then
        Write(50,*) 'Bin                : ', NBin
        Write(50,*) 'No CPU-time limitation '
     else
        Write(50,'("Prog will stop after hours:",2x,F8.4)') CPU_MAX
     endif
     close(50)
#ifdef MPI
  endif
#endif

  
  Allocate ( UST(NDIM,NDIM,NSTM,N_FL), VST(NDIM,NDIM,NSTM,N_FL), DST(NDIM,NSTM,N_FL) )
  Allocate ( Test(Ndim,Ndim) )

  NST = NSTM 
  !Write(6,*) "Write UL ", NST
  Do nf = 1,N_FL
     CALL INITD(UL(:,:,Nf),Z_ONE)
     DL(:,nf) = Z_ONE
     CALL INITD(VL(:,:,nf),Z_ONE)
     UST(:,:,NST,nf) = UL(:,:,nf)
     VST(:,:,NST,nf) = VL(:,:,nf)
     DST(:,NST,nf) = DL(:,nf)
     CALL INITD(UR(:,:,nf),Z_ONE)
     CALL INITD(VR(:,:,nf),Z_ONE)
     DR(:,nf) = Z_ONE
  Enddo


  DO NST = NSTM-1,1,-1
     NT1 = Stab_nt(NST+1)
     NT  = Stab_nt(NST  )
     CALL WRAPUL(NT1,NT,UL,DL, VL)
     Do nf = 1,N_FL
        UST(:,:,NST,nf) = UL(:,:,nf)
        VST(:,:,NST,nf) = VL(:,:,nf)
        DST(:  ,NST,nf) = DL(:  ,nf)
     ENDDO
  ENDDO
  NT1 = stab_nt(1)
  CALL WRAPUL(NT1,0, UL ,DL, VL)

!!$  DO NT = LTROT-NWRAP,NWRAP,-1
!!$     IF ( MOD(NT,NWRAP)  == 0 ) THEN
!!$        NT1 = NT + NWRAP
!!$        !Write(6,*) 'Calling Wrapul:', NT1,NT 
!!$        CALL WRAPUL(NT1,NT,UL,DL, VL)
!!$        NST = NINT( DBLE(NT)/DBLE(NWRAP) )
!!$        !Write(6,*) "Write UL ", NST
!!$        Do nf = 1,N_FL
!!$           DO I = 1,Ndim
!!$              DO J = 1,Ndim
!!$                 UST(I,J,NST,nf) = UL(I,J,nf)
!!$                 VST(I,J,NST,nf) = VL(I,J,nf)
!!$              ENDDO
!!$           ENDDO
!!$           DO I = 1,Ndim
!!$              DST(I,NST,nf) = DL(I,nf)
!!$           ENDDO
!!$        ENDDO
!!$     ENDIF
!!$  ENDDO
!!$  CALL WRAPUL(NWRAP,0, UL ,DL, VL)
  !WRITE(6,*) 'Filling up storage'
  !Write(6,*) 'Done wrapping'
  

  NVAR = 1
  Phase = cmplx(1.d0, 0.d0, kind(0.D0))
  do nf = 1,N_Fl
     CALL CGR(Z, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
     Phase = Phase*Z
  Enddo
  call Op_phase(Phase,OP_V,Nsigma,N_SUN) 
#ifdef MPI 
  !WRITE(6,*) 'Phase is: ', Irank, PHASE, GR(1,1,1)
#else
  !WRITE(6,*) 'Phase is: ',  PHASE
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

     Call cpu_time(time_bin_start) 
     
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
        
        NST = 1
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           !Write(6,*) "Hi"
           CALL WRAPGRUP(GR,NTAU,PHASE) 
           !Write(6,*) "Hi1"
!!$           IF ( MOD(NTAU1,NWRAP ) .EQ. 0 ) THEN
!!$              NST = NINT( DBLE(NTAU1)/DBLE(NWRAP) )
!!$              NT1 = NTAU1 - NWRAP
           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              !Write(6,*) 'Wrapur : ', NT1, NTAU1
              CALL WRAPUR(NT1, NTAU1,UR, DR, VR)
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf = 1, N_FL
                 UL(:,:,nf) = UST(:,:,NST,nf)
                 VL(:,:,nf) = VST(:,:,NST,nf)
                 DL(:,  nf) = DST(:  ,NST,nf)
                 ! Write in store Right prop from 1 to LTROT/NWRAP
                 ! Write(6,*)  'Write UR, read UL ', NTAU1, NST
                 UST(:,:,NST,nf) = UR(:,:,nf)
                 VST(:,:,NST,nf) = VR(:,:,nf)
                 DST(:  ,NST,nf) = DR(:  ,nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 !Write(6,*) ' Call Cgr'
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
                 Z = Z*Z1
                 !Write(6,*) 'Calling control ',NTAU1, Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST + 1
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
           DL(:,nf)  = Z_ONE
        ENDDO
        
        NST = NSTM-1
        DO NTAU = LTROT,1,-1
           NTAU1 = NTAU - 1
           CALL WRAPGRDO(GR,NTAU, PHASE)
           IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
              CALL Obser( GR, PHASE, Ntau1 )
           ENDIF
!!$           IF ( MOD(NTAU1,NWRAP).EQ.0 .AND. NTAU1.NE.0 ) THEN
!!$              ! WRITE(50,*) 'Recalc at  :', NTAU1
!!$              NST = NINT( DBLE(NTAU1)/DBLE(NWRAP) )
!!$              NT1 = NTAU1 + NWRAP
           IF ( Stab_nt(NST) == NTAU1 .AND. NTAU1.NE.0 ) THEN
              NT1 = Stab_nt(NST+1)
              !Write(6,*) 'Wrapul : ', NT1, NTAU1
              CALL WRAPUL(NT1,NTAU1, UL, DL, VL )
              !Write(6,*)  'Write UL, read UR ', NTAU1, NST
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              do nf = 1,N_FL
                 UR(:,:,nf) = UST(:,:,NST,nf)
                 VR(:,:,nf) = VST(:,:,NST,nf)
                 DR(:  ,nf) = DST(:  ,NST,nf)
                 ! WRITE in store the left prop. from LTROT/NWRAP-1 to 1
                 UST(:,:,NST,nf) =  UL(:,:,nf)
                 VST(:,:,NST,nf) =  VL(:,:,nf)
                 DST(:  ,NST,nf) =  DL(:  ,nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 !Write(6,*) ' Call Cgr'
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
                 Z = Z*Z1
                 !Write(6,*) 'Calling control: ', NTAU1, Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST -1
           ENDIF
        ENDDO

        !Calculate and compare green functions on time slice 0.
        NT1 = Stab_nt(0)
        NT  = Stab_nt(1)
        !Write(6,*) 'Wrapul : ', NT, NT1
        CALL WRAPUL(NT,NT1, UL, DL, VL )
        
        do nf = 1,N_FL
           CALL INITD(UR(:,:,nf),Z_ONE)
           CALL INITD(VR(:,:,nf),Z_ONE)
           DR(:,nf) = Z_ONE
        ENDDO
        Z = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf = 1,N_FL
           TEST(:,:) = GR(:,:,nf)
           NVAR = 1
           CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
           Z = Z*Z1
           !Write(6,*) 'Calling control  0', Z1
           Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
        ENDDO
        call Op_phase(Z,OP_V,Nsigma,N_SUN) 
        Call Control_PrecisionP(Z,Phase)
        Phase = Z
        NST =  NSTM
        Do nf = 1,N_FL
           UST(:,:,NST,nf) = CMPLX(0.D0, 0.D0, kind(0.D0))
           VST(:,:,NST,nf) = CMPLX(0.D0, 0.D0, kind(0.D0))
           DO I = 1,Ndim
              DST(I  ,NST,nf)  = CMPLX(1.D0, 0.D0, kind(0.D0))
              UST(I,I,NST,nf)  = CMPLX(1.D0, 0.D0, kind(0.D0))
              VST(I,I,NST,nf)  = CMPLX(1.D0, 0.D0, kind(0.D0))
           ENDDO
        enddo
        IF ( LTAU == 1 ) then
!!$#ifdef MPI
!!$           Write(6,*) Irank, 'Calling Tau_m', NWRAP, NSTM
!!$#else
!!$           Write(6,*)  'Calling Tau_m', NWRAP, NSTM
!!$#endif

           Call TAU_M( UST,DST,VST, GR, PHASE, NSTM, NWRAP, STAB_NT ) 
!!$#ifdef MPI
!!$           Write(6,*) Irank, 'Back Calling Tau_m'
!!$#else
!!$           Write(6,*)  'Back Calling Tau_m'
!!$#endif
        endif
           
     ENDDO
     Call Pr_obs(Ltau)
     Call confout
          
     Call cpu_time(time_bin_end)
     prog_truncation = .false.
     if ( abs(CPU_MAX) > Zero ) then
        Call make_truncation(prog_truncation,CPU_MAX,time_bin_start,time_bin_end)
     endif
     If (prog_truncation) then 
        Nbin_eff = nbc
        exit !exit the loop over the bin index, labelled NBC.
     Endif
  Enddo
  Call Control_Print

#ifdef MPI
  If (Irank == 0 ) then
#endif
     if ( abs(CPU_MAX) > Zero ) then
        Open (Unit=50,file="info", status="unknown", position="append")
        Write(50,*)' Effective number of bins   : ', Nbin_eff
        Close(50)
     endif
#ifdef MPI
  endif
#endif
 
#ifdef MPI
   CALL MPI_FINALIZE(ierr)
#endif

end Program Main
