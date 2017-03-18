!  Copyright (C) 2016, 2017 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


Program Main

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Main program. Reads in VAR_QMC  namelist.  Calls Ham_set. Carries 
!> out the sweeps. 
!
!--------------------------------------------------------------------

  Use Operator_mod
  Use Lattices_v3 
  Use MyMats 
  Use Hamiltonian
  Use Control
  Use Tau_m_mod
  Use Hop_mod
  Use Global_mod
 
  Implicit none
#ifdef MPI
  include 'mpif.h'
#endif   
#include "git.h"


  Interface
     SUBROUTINE WRAPUL(NTAU1, NTAU, UL ,DL, VL)
       Use Hamiltonian
       Implicit none
       COMPLEX (Kind=Kind(0.d0)) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
       COMPLEX (Kind=Kind(0.d0)) :: DL(Ndim,N_FL)
       Integer :: NTAU1, NTAU
     END SUBROUTINE WRAPUL
     SUBROUTINE CGR(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)
       Use UDV_Wrap_mod
       Implicit None
       COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(In)    :: URUP, VRUP, ULUP, VLUP
       COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(In)      :: DLUP, DRUP
       COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
       COMPLEX(Kind=Kind(0.d0)) :: PHASE
       INTEGER         :: NVAR
     END SUBROUTINE CGR
     SUBROUTINE WRAPGRUP(GR,NTAU,PHASE)
       Use Hamiltonian
       Implicit none
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  GR(Ndim,Ndim,N_FL)
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  PHASE
       INTEGER, INTENT(IN) :: NTAU
     END SUBROUTINE WRAPGRUP
     SUBROUTINE WRAPGRDO(GR,NTAU,PHASE)
       Use Hamiltonian 
       Implicit None
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(NDIM,NDIM,N_FL)
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
       Integer :: NTAU
     end SUBROUTINE WRAPGRDO
     SUBROUTINE WRAPUR(NTAU, NTAU1, UR, DR, VR)
       Use Hamiltonian
       Use UDV_Wrap_mod
       Implicit None
       COMPLEX (Kind=Kind(0.d0)) :: UR(Ndim,Ndim,N_FL), VR(Ndim,Ndim,N_FL)
       COMPLEX (Kind=Kind(0.d0)) :: DR(Ndim,N_FL)
       Integer :: NTAU1, NTAU
     END SUBROUTINE WRAPUR

  end Interface

  COMPLEX (Kind=Kind(0.d0)), Dimension(:,:)  , Allocatable   ::  TEST

  COMPLEX (Kind=Kind(0.d0)), Dimension(:,:)  , Allocatable    :: DL, DR
  COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable    :: UL, VL, UR, VR
  COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable    :: GR
  

  Integer :: Nwrap, NSweep, NBin, NBin_eff,Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
  Integer :: NTAU, NTAU1
  Real(Kind=Kind(0.d0)) :: CPU_MAX 
  Character (len=64) :: file1

  
#if defined(TEMPERING)
  Integer :: N_exchange_steps, N_Tempering_frequency
  NAMELIST /VAR_TEMP/  N_exchange_steps, N_Tempering_frequency
#endif

  NAMELIST /VAR_QMC/   Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX 


  Integer :: Ierr, I,nf, nst, n
  Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Phase, Z, Z1
  Real    (Kind=Kind(0.d0)) :: ZERO = 10D-8
  Integer, dimension(:), allocatable :: Stab_nt
  
  ! Space for storage.
  COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:)  , Allocatable :: DST
  COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable :: UST,  VST

  ! For tests
  Real (Kind=Kind(0.d0)) :: Weight, Weight_tot
  Integer :: nr,nth, nth1
  Logical :: Log
  
  ! For the truncation of the program:
  logical                   :: prog_truncation
  integer (kind=kind(0.d0)) :: count_bin_start, count_bin_end

#ifdef MPI
  Integer        ::  Isize, Irank
  INTEGER        :: STATUS(MPI_STATUS_SIZE)
#endif


#ifdef MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

#ifdef MPI
  If (  Irank == 0 ) then
#endif
     write (*,*) "ALF Copyright (C) 2016,2017 The ALF project contributors"
     write (*,*) "This Program comes with ABSOLUTELY NO WARRANTY; for details see license.GPL"
     write (*,*) "This is free software, and you are welcome to redistribute it under certain conditions."
#ifdef MPI
  endif
#endif
 


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
#if defined(TEMPERING) 
     OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(*,*) 'unable to open <parameters>',ierr
        STOP
     END IF
     READ(5,NML=VAR_TEMP)
     CLOSE(5)
#endif     
#ifdef MPI
  Endif
  CALL MPI_BCAST(Nwrap          ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NSweep         ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NBin           ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Ltau           ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(LOBS_EN        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(LOBS_ST        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(CPU_MAX        ,1,MPI_REAL8,  0,MPI_COMM_WORLD,ierr)
#if defined(TEMPERING) 
   CALL MPI_BCAST(N_exchange_steps        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(N_Tempering_frequency   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif  
#endif

 
  Call Ham_set
  Call confin 
  Call Hop_mod_init

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

#if defined(TEMPERING)
           write(File1,'(A,I0,A)') "Temp_",Irank,"/info"
#else
           File1 = "info"
#endif
           
#if defined(MPI) && !defined(TEMPERING)
  if ( Irank == 0 ) then
#endif
     Open (Unit = 50,file=file1,status="unknown",position="append")
     Write(50,*) 'Sweeps             : ', Nsweep
     Write(50,*) 'Measure Int.       : ', LOBS_ST, LOBS_EN
     Write(50,*) 'Stabilization,Wrap : ', Nwrap
     Write(50,*) 'Nstm               : ', NSTM
     Write(50,*) 'Ltau               : ', Ltau     
#if defined(MPI) && !defined(TEMPERING)
     Write(50,*) 'Number of  threads : ', ISIZE
#endif   
#if defined(GIT)
     Write(50,*) 'This executable represents commit '&
&      , GIT_COMMIT_HASH , ' of branch ' , GIT_BRANCH , '.'
#endif
     If ( abs(CPU_MAX) < ZERO ) then
        Write(50,*) 'Bin                : ', NBin
        Write(50,*) 'No CPU-time limitation '
     else
        Write(50,'(" Prog will stop after hours:",2x,F8.4)') CPU_MAX
     endif
#if defined(STAB1) 
     Write(50,*) 'STAB1 is defined '
#endif
#if defined(STAB2) 
     Write(50,*) 'STAB2 is defined '
#endif
#if defined(QRREF) 
     Write(50,*) 'QRREF is defined '
#endif
#if defined(TEMPERING) 
     Write(50,*) '# of exchange steps  ',N_exchange_steps 
     Write(50,*) 'Tempering frequency  ',N_Tempering_frequency
#endif
     close(50)
#if defined(MPI) && !defined(TEMPERING)
  endif
#endif


  !Call Test_Hamiltonian

  
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
     !Write(6,*)'Hi', NT1,NT, NST
     CALL WRAPUL(NT1,NT,UL,DL, VL)
     Do nf = 1,N_FL
        UST(:,:,NST,nf) = UL(:,:,nf)
        VST(:,:,NST,nf) = VL(:,:,nf)
        DST(:  ,NST,nf) = DL(:  ,nf)
     ENDDO
  ENDDO
  NT1 = stab_nt(1)
  CALL WRAPUL(NT1,0, UL ,DL, VL)

  

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
#endif



  Call Control_init

  DO  NBC = 1, NBIN
     ! Here, you have the green functions on time slice 1.
     ! Set bin observables to zero.

     call system_clock(count_bin_start)
     
     Call Init_obs(Ltau)

     DO NSW = 1, NSWEEP

#if defined(TEMPERING)
        IF (MOD(NSW,N_Tempering_frequency) == 0) then
           !Write(6,*) "Irank, Call tempering", Irank, NSW
           CALL Exchange_Step(Phase,GR,UR,DR,VR, UL,DL,VL,Stab_nt, UST, VST, DST,N_exchange_steps)
        endif
#endif
        ! Global updates
        If (Global_moves) Call Global_Updates(Phase,GR,UR,DR,VR, UL,DL,VL,Stab_nt, UST, VST, DST)

        ! Propagation from 1 to Ltrot
        ! Set the right storage to 1
        
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
           CALL WRAPGRUP(GR,NTAU,PHASE) 
           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              CALL WRAPUR(NT1, NTAU1,UR, DR, VR)
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf = 1, N_FL
                 ! Read from storage left propagation from LTROT to  NTAU1
                 UL(:,:,nf) = UST(:,:,NST,nf)
                 VL(:,:,nf) = VST(:,:,NST,nf)
                 DL(:,  nf) = DST(:  ,NST,nf)
                 ! Write in storage right prop from 1 to NTAU1
                 UST(:,:,NST,nf) = UR(:,:,nf)
                 VST(:,:,NST,nf) = VR(:,:,nf)
                 DST(:  ,NST,nf) = DR(:  ,nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
                 Z = Z*Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST + 1
           ENDIF

           IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
              CALL Obser( GR, PHASE, Ntau1 )
           ENDIF
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
           IF ( Stab_nt(NST) == NTAU1 .AND. NTAU1.NE.0 ) THEN
              NT1 = Stab_nt(NST+1)
              !Write(6,*) 'Wrapul : ', NT1, NTAU1
              CALL WRAPUL(NT1,NTAU1, UL, DL, VL )
              !Write(6,*)  'Write UL, read UR ', NTAU1, NST
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              do nf = 1,N_FL
                 ! Read from store the right prop. from 1 to LTROT/NWRAP-1
                 UR(:,:,nf) = UST(:,:,NST,nf)
                 VR(:,:,nf) = VST(:,:,NST,nf)
                 DR(:  ,nf) = DST(:  ,NST,nf)
                 ! WRITE in store the left prop. from LTROT/NWRAP-1 to 1
                 UST(:,:,NST,nf) =  UL(:,:,nf)
                 VST(:,:,NST,nf) =  VL(:,:,nf)
                 DST(:  ,NST,nf) =  DL(:  ,nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
                 Z = Z*Z1
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
           ! Call for imaginary time displaced  correlation fuctions. 
           Call TAU_M( UST,DST,VST, GR, PHASE, NSTM, NWRAP, STAB_NT ) 
        endif
           
     ENDDO
     Call Pr_obs(Ltau)
     Call confout

     call system_clock(count_bin_end)
     prog_truncation = .false.
     if ( abs(CPU_MAX) > Zero ) then
        Call make_truncation(prog_truncation,CPU_MAX,count_bin_start,count_bin_end)        
     endif
     If (prog_truncation) then 
        Nbin_eff = nbc
        exit !exit the loop over the bin index, labelled NBC.
     Endif
  Enddo
  DEALLOCATE(DL, DR, UL, UR, VR, VL, GR, UST, VST, DST, TEST, Stab_nt)
  Call Control_Print

#if defined(MPI) && !defined(TEMPERING)
  If (Irank == 0 ) then
#endif
     if ( abs(CPU_MAX) > Zero ) then
#if defined(TEMPERING)
        write(File1,'(A,I0,A)') "Temp_",Irank,"/info"
#else
        File1 = "info"
#endif
        Open (Unit=50,file=File1, status="unknown", position="append")
        Write(50,*)' Effective number of bins   : ', Nbin_eff
        Close(50)
     endif
#if defined(MPI) && !defined(TEMPERING)
  endif
#endif
 
#ifdef MPI
   CALL MPI_FINALIZE(ierr)
#endif

end Program Main
