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
        Use UDV_State_mod
        Use Wrapgr_mod
        
        Implicit none
#ifdef MPI
        include 'mpif.h'
#endif   
#include "git.h"
        
        
        Interface
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout) :: UDVL(N_FL)
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: UDVL, UDVR
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout) :: UDVR(N_FL)
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           
        end Interface
        
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:)  , Allocatable   ::  TEST
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable    :: GR
        CLASS(UDV_State), DIMENSION(:), ALLOCATABLE :: udvl, udvr
  
        
        Integer :: Nwrap, NSweep, NBin, NBin_eff,Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
        Integer :: NTAU, NTAU1
        Real(Kind=Kind(0.d0)) :: CPU_MAX 
        Character (len=64) :: file1
  
        ! Space for choosing sampling scheme
        Logical :: Propose_S0
        Logical :: Global_moves, Global_tau_moves
        Integer :: N_Global 
        Integer :: Nt_sequential_start, Nt_sequential_end, mpi_per_parameter_set
        Integer :: N_Global_tau
        
        
#if defined(TEMPERING)
        Integer :: N_exchange_steps, N_Tempering_frequency
        NAMELIST /VAR_TEMP/  N_exchange_steps, N_Tempering_frequency, mpi_per_parameter_set
#endif
        
        NAMELIST /VAR_QMC/   Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX, &
             &               Propose_S0,Global_moves,  N_Global, Global_tau_moves, &
             &               Nt_sequential_start, Nt_sequential_end, N_Global_tau
        
        
        Integer :: Ierr, I,nf, nst, n
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Phase, Z, Z1
        Real    (Kind=Kind(0.d0)) :: ZERO = 10D-8
        Integer, dimension(:), allocatable :: Stab_nt
  
        ! Space for storage.
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE :: udvst
        
        ! For tests
        Real (Kind=Kind(0.d0)) :: Weight, Weight_tot
        Integer :: nr,nth, nth1
        Logical :: Log
        
        ! For the truncation of the program:
        logical                   :: prog_truncation
        integer (kind=kind(0.d0)) :: count_bin_start, count_bin_end
        
#ifdef MPI
        Integer        :: Isize, Irank, Irank_g, Isize_g, color, key, igroup
        INTEGER        :: STATUS(MPI_STATUS_SIZE)
#endif
        
        
#ifdef MPI
        CALL MPI_INIT(ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
        
#if defined(TEMPERING) && defined(MPI)
        mpi_per_parameter_set = 1  ! Default value
        OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(*,*) 'unable to open <parameters>',ierr
           STOP
        END IF
        READ(5,NML=VAR_TEMP)
        CLOSE(5)
        CALL MPI_BCAST(N_exchange_steps        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(N_Tempering_frequency   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(mpi_per_parameter_set   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if ( mod(ISIZE,mpi_per_parameter_set) .ne. 0 ) then
           Write (6,*) "mpi_per_parameter_set is not a multiple of total mpi processes"
           stop
        endif
        Call Global_Tempering_setup
#elif !defined(TEMPERING)  && defined(MPI)
        mpi_per_parameter_set = Isize
#elif defined(TEMPERING)  && !defined(MPI)
        Write(6,*) 'Mpi has to be defined for tempering runs'
#endif

#ifdef MPI
        color = irank/mpi_per_parameter_set
        key   =  0
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Group_comm, ierr)
        call MPI_Comm_rank(Group_Comm, Irank_g, ierr)
        call MPI_Comm_size(Group_Comm, Isize_g, ierr)
        igroup           = irank/isize_g
        !Write(6,*) 'irank, Irank_g, Isize_g', irank, irank_g, isize_g
#endif

#ifdef MPI
        If (  Irank == 0 ) then
#endif
           write (*,*) "ALF Copyright (C) 2016, 2017 The ALF project contributors"
           write (*,*) "This Program comes with ABSOLUTELY NO WARRANTY; for details see license.GPL"
           write (*,*) "This is free software, and you are welcome to redistribute it under certain conditions."
#ifdef MPI
        endif
#endif
 
#ifdef MPI
        If ( Irank == 0 ) then 
#endif
           ! This is a set of variables that  identical for each simulation.
           Nwrap=0;  NSweep=0; NBin=0; Ltau=0; LOBS_EN = 0;  LOBS_ST = 0;  CPU_MAX = 0.d0
           Propose_S0 = .false. ;  Global_moves = .false. ; N_Global = 0
           Global_tau_moves = .false. 
           Nt_sequential_start = 0 ;  Nt_sequential_end  = 0;  N_Global_tau  = 0
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
        CALL MPI_BCAST(Nwrap              ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(NSweep             ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(NBin               ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Ltau               ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(LOBS_EN            ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(LOBS_ST            ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(CPU_MAX            ,1,MPI_REAL8,  0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Propose_S0         ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Global_moves       ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(N_Global           ,1,MPI_Integer,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Global_tau_moves   ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Nt_sequential_start,1,MPI_Integer,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Nt_sequential_end  ,1,MPI_Integer,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(N_Global_tau       ,1,MPI_Integer,0,MPI_COMM_WORLD,ierr)
#endif
        
 
        Call Ham_set
        If ( .not. Global_tau_moves )  then
           ! This  corresponds to the default updating scheme
           Nt_sequential_start = 1 
           Nt_sequential_end   = Size(OP_V,1) 
           N_Global_tau        = 0
        endif

        
        Call confin 
        Call Hop_mod_init

        IF (ABS(CPU_MAX) > Zero ) NBIN = 10000000
        If (N_Global_tau > 0) then
           Call Wrapgr_alloc
        endif
 
        Call control_init
        Call Alloc_obs(Ltau)
        Call Op_SetHS

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
        write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
        File1 = "info"
#endif
           
#if defined(MPI) 
        if ( Irank_g == 0 ) then
#endif
           Open (Unit = 50,file=file1,status="unknown",position="append")
           Write(50,*) 'Sweeps                              : ', Nsweep
           If ( abs(CPU_MAX) < ZERO ) then
              Write(50,*) 'Bins                                : ', NBin
              Write(50,*) 'No CPU-time limitation '
           else
              Write(50,'(" Prog will stop after hours:",2x,F8.4)') CPU_MAX
           endif
           Write(50,*) 'Measure Int.                        : ', LOBS_ST, LOBS_EN
           Write(50,*) 'Stabilization,Wrap                  : ', Nwrap
           Write(50,*) 'Nstm                                : ', NSTM
           Write(50,*) 'Ltau                                : ', Ltau     
           Write(50,*) '# of interacting Ops per time slice : ', Size(OP_V,1)
           If ( Propose_S0 ) &
                &  Write(50,*) 'Propose Ising moves according to  bare Ising action'
           If ( Global_moves ) Then
              Write(50,*) 'Global moves are enabled   '
              Write(50,*) '# of global moves / sweep :', N_Global
           Endif
           If ( Global_tau_moves ) Then
              Write(50,*) 'Nt_sequential_start: ', Nt_sequential_start
              Write(50,*) 'Nt_sequential_end  : ', Nt_sequential_end
              Write(50,*) 'N_Global_tau       : ', N_Global_tau
           else
              Write(50,*) 'Default sequential updating '
           endif
           
#if defined(MPI) 
           Write(50,*) 'Number of mpi-processes : ', isize_g
#endif   
#if defined(GIT)
           Write(50,*) 'This executable represents commit '&
                &      , GIT_COMMIT_HASH , ' of branch ' , GIT_BRANCH , '.'
#endif
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
#if defined(MPI) 
        endif
#endif


        !Call Test_Hamiltonian
        Allocate ( Test(Ndim,Ndim), GR(NDIM,NDIM,N_FL ) )
        ALLOCATE(udvl(N_FL), udvr(N_FL), udvst(NSTM, N_FL))
        do nf = 1, N_FL
           CALL udvl(nf)%init(ndim)
           CALL udvr(nf)%init(ndim)
           do n = 1, NSTM
              CALL udvst(n,nf)%alloc(ndim)
           ENDDO
           CALL udvst(NSTM, nf)%reset
        enddo
        
        DO NST = NSTM-1,1,-1
           NT1 = Stab_nt(NST+1)
           NT  = Stab_nt(NST  )
           !Write(6,*)'Hi', NT1,NT, NST
           CALL WRAPUL(NT1, NT, UDVL)
           Do nf = 1,N_FL
              UDVST(NST, nf) = UDVL(nf)
           ENDDO
        ENDDO
        NT1 = stab_nt(1)
        CALL WRAPUL(NT1, 0, UDVL)
        
        
        
        NVAR = 1
        Phase = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf = 1,N_Fl
           CALL CGR(Z, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
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
#if defined(TEMPERING)
           Call Global_Tempering_init_obs
#endif           
           
           DO NSW = 1, NSWEEP
              
#if defined(TEMPERING)
              IF (MOD(NSW,N_Tempering_frequency) == 0) then
                 !Write(6,*) "Irank, Call tempering", Irank, NSW
                 CALL Exchange_Step(Phase,GR,udvr, udvl,Stab_nt, udvst, N_exchange_steps)
              endif
#endif
              ! Global updates
              If (Global_moves) Call Global_Updates(Phase, GR, udvr, udvl, Stab_nt, udvst,N_Global)
              
              ! Propagation from 1 to Ltrot
              ! Set the right storage to 1
              
              do nf = 1,N_FL
                 CALL udvr(nf)%reset
              Enddo
              
              NST = 1
              DO NTAU = 0, LTROT-1
                 NTAU1 = NTAU + 1
                 CALL WRAPGRUP(GR,NTAU,PHASE,Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)
                 
                 If (NTAU1 == Stab_nt(NST) ) then 
                    NT1 = Stab_nt(NST-1)
                    CALL WRAPUR(NT1, NTAU1, udvr)
                    Z = cmplx(1.d0, 0.d0, kind(0.D0))
                    Do nf = 1, N_FL
                       ! Read from storage left propagation from LTROT to  NTAU1
                       udvl(nf) = udvst(NST, nf)
                       ! Write in storage right prop from 1 to NTAU1
                       udvst(NST, nf) = udvr(nf)
                       NVAR = 1
                       IF (NTAU1 .GT. LTROT/2) NVAR = 2
                       TEST(:,:) = GR(:,:,nf)
                       CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
                       Z = Z*Z1
                       Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
                    ENDDO
                    call Op_phase(Z,OP_V,Nsigma,N_SUN) 
                    Call Control_PrecisionP(Z,Phase)
                    Phase = Z
                    NST = NST + 1
                 ENDIF
                 
                 IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
                    !Call  Global_tau_mod_Test(Gr,ntau1)
                    !Stop
                    CALL Obser( GR, PHASE, Ntau1 )
                 ENDIF
              ENDDO
              
              Do nf = 1,N_FL
                 CALL udvl(nf)%reset
              ENDDO
              
              NST = NSTM-1
              DO NTAU = LTROT,1,-1
                 NTAU1 = NTAU - 1
                 CALL WRAPGRDO(GR,NTAU, PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
                 IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
                    CALL Obser( GR, PHASE, Ntau1 )
                 ENDIF
                 IF ( Stab_nt(NST) == NTAU1 .AND. NTAU1.NE.0 ) THEN
                    NT1 = Stab_nt(NST+1)
                    !Write(6,*) 'Wrapul : ', NT1, NTAU1
                    CALL WRAPUL(NT1, NTAU1, udvl)
                    !Write(6,*)  'Write UL, read UR ', NTAU1, NST
                    Z = cmplx(1.d0, 0.d0, kind(0.D0))
                    do nf = 1,N_FL
                       ! Read from store the right prop. from 1 to LTROT/NWRAP-1
                       udvr(nf) = udvst(NST, nf)
                       ! WRITE in store the left prop. from LTROT/NWRAP-1 to 1
                       udvst(NST, nf) = udvl(nf)
                       NVAR = 1
                       IF (NTAU1 .GT. LTROT/2) NVAR = 2
                       TEST(:,:) = GR(:,:,nf)
                       CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
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
              CALL WRAPUL(NT, NT1, udvl)
        
              do nf = 1,N_FL
                 CALL udvr(nf)%reset
              ENDDO
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              do nf = 1,N_FL
                 TEST(:,:) = GR(:,:,nf)
                 NVAR = 1
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
                 Z = Z*Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST =  NSTM
              Do nf = 1,N_FL
                 CALL udvst(NST, nf)%reset
              enddo
              IF ( LTAU == 1 ) then
                 ! Call for imaginary time displaced  correlation fuctions. 
                 Call TAU_M( udvst, GR, PHASE, NSTM, NWRAP, STAB_NT ) 
              endif
              
           ENDDO
           Call Pr_obs(Ltau)
#if defined(TEMPERING)
           Call Global_Tempering_Pr
#endif           

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
        
        ! Deallocate things
        DO nf = 1, N_FL
           CALL udvl(nf)%dealloc
           CALL udvr(nf)%dealloc
           do n = 1, NSTM
              CALL udvst(n,nf)%dealloc
           ENDDO
        ENDDO
        DEALLOCATE(udvl, udvr, udvst)
        DEALLOCATE(GR, TEST, Stab_nt)
        If (N_Global_tau > 0) then
           Call Wrapgr_dealloc
        endif
        
        Call Control_Print(Group_Comm)
        
#if defined(MPI) 
        If (Irank_g == 0 ) then
#endif
           if ( abs(CPU_MAX) > Zero ) then
#if defined(TEMPERING)
              write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
              File1 = "info"
#endif
              Open (Unit=50,file=File1, status="unknown", position="append")
              Write(50,*)' Effective number of bins   : ', Nbin_eff
              Close(50)
           endif
#if defined(MPI) 
        endif
#endif
        
#ifdef MPI
        CALL MPI_FINALIZE(ierr)
#endif
        
      end Program Main
