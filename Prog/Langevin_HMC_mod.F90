!  Copyright (C) 2016 - 2020 The ALF project
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
!     along with ALF.  If not, usee http://www.gnu.org/licenses/.
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


      Module Langevin_HMC_mod
        
        Use Hamiltonian
        Use UDV_State_mod
        Use Control
        Use Hop_mod
        Use iso_fortran_env, only: output_unit, error_unit

        
        Implicit none

        Private
        
        Public :: Langevin_HMC, Langevin_HMC_type
        
        Type Langevin_HMC_type
           private
           Character (Len=64)                      :: Update_scheme
           Logical                                 :: L_Forces
           Real    (Kind=Kind(0.d0))               :: Delta_t_running, Delta_t_Langevin_HMC, Max_Force
           Complex (Kind=Kind(0.d0)),  pointer :: Forces  (:,:)
           
           Real    (Kind=Kind(0.d0)),  pointer :: Forces_0(:,:)
         CONTAINS
           procedure  ::    make        => Langevin_HMC_setup
           procedure  ::    clean       => Langevin_HMC_clear
           procedure  ::    Wrap_Forces => Wrapgrup_Forces
           procedure  ::    Update      => Langevin_HMC_update
           procedure  ::    set_L_Forces         => Langevin_HMC_set_L_Forces
           procedure  ::    get_Update_scheme    => Langevin_HMC_get_Update_scheme
           procedure  ::    set_Update_scheme    => Langevin_HMC_set_Update_scheme
           procedure  ::    get_Delta_t_running  => Langevin_HMC_get_Delta_t_running
        end type Langevin_HMC_type

        Type (Langevin_HMC_type) :: Langevin_HMC

           
      Contains

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>   Computes the  forces as well as, on demand,  observables.  
!>   On input:  a)  GR is on the first time slice and  the storage is full with left propagations.
!>              b)  Udvl is on time slice 0.
!>   On output.
!>              a)  Forces (only for field type 3 (i.e. continuous fieds)  are computed   and stored in Forces(:,:)
!>              
!> 
!--------------------------------------------------------------------
        
      SUBROUTINE  Langevin_HMC_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN,Calc_Obser_eq)
        Implicit none
        
        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
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
        end Interface
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(in), allocatable, dimension(:,:)    :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout)                     :: Phase
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:)   :: Test
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR, GR_Tilde
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
        Integer, intent(in) :: LOBS_ST, LOBS_EN
        Logical, intent(in) :: Calc_Obser_eq 
        

        !Local
        Integer :: NSTM, n, nf, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Real    (Kind=Kind(0.d0)) :: spin
        
        NSTM = Size(Stab_nt,1) - 1 
        !Do  n = 0,NSTM
        !   Write(6,*)  n, Stab_nt(n)
        !Enddo
        
        Langevin_HMC%Forces = cmplx(0.d0,0.d0,Kind(0.d0))
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        Enddo
        NST = 1
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           
           Call  Langevin_HMC%Wrap_Forces(Gr, ntau1)
           
           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              CALL WRAPUR(NT1, NTAU1, udvr)
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf = 1, N_FL
                 ! Read from storage left propagation from LTROT to  NTAU1
                 udvl(nf) = udvst(NST, nf)
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
           
           IF (NTAU1 .GE. LOBS_ST .AND. NTAU1 .LE. LOBS_EN .and. Calc_Obser_eq ) THEN
              If (Symm) then
                 Call Hop_mod_Symm(GR_Tilde,GR)
                 CALL Obser( GR_Tilde, PHASE, Ntau1,Langevin_HMC%Delta_t_running )
              else
                 CALL Obser( GR, PHASE, Ntau1, Langevin_HMC%Delta_t_running )
              endif
           endif
        enddo

      end SUBROUTINE Langevin_HMC_Forces

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>   In :  Gr is on time slice NT
!>   Out:  Gr is on time slice NT1=NT+1  and the forces on time slice NT1 are computed.
!>         and stored in Forces(:,NT1)
!>   
!--------------------------------------------------------------------

      Subroutine  Wrapgrup_Forces(this,Gr, NT1)
        
        Implicit none
        
        class (Langevin_HMC_type) :: this
        Complex (Kind=Kind(0.d0)), intent(inout), dimension(:,:,:) :: Gr
        Integer, intent(in)                                        :: nt1

        
        !Local
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Integer ::  nf, I, J, n, N_type
        Real(Kind=Kind(0.d0)) :: spin


        Do nf = 1,N_FL
           CALL HOP_MOD_mmthr   (GR(:,:,nf), nf )
           CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf )
        Enddo
        Do n = 1, size(OP_V,1) 
           this%Forces(n,nt1)  = cmplx(0.d0,0.d0,Kind(0.d0))
           Do nf = 1, N_FL
              spin = nsigma%f(n,nt1) ! Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
              N_type = 1
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
              N_type =  2
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
           enddo
           if (OP_V(n,1)%type == 3 ) then
              Do nf = 1, N_Fl
                 Z = cmplx(0.d0,0.d0,Kind(0.d0))
                 do I = 1,size(OP_V(n,nf)%P,1)
                    do J = 1,size(OP_V(n,nf)%P,1)
                       Z1 =  cmplx(0.d0,0.d0,Kind(0.d0))
                       if ( I == J ) Z1 = cmplx(1.d0,0.d0,Kind(0.d0))
                       Z  = Z +    Op_V(n,nf)%O(I,J) * ( Z1 - Gr(Op_V(n,nf)%P(J),Op_V(n,nf)%P(I), nf) )
                    Enddo
                 Enddo
                 this%Forces(n,nt1) =  this%Forces(n,nt1)  - &
                      &    Op_V(n,nf)%g * Z *  cmplx(real(N_SUN,Kind(0.d0)), 0.d0, Kind(0.d0)) 
              Enddo
           endif
        enddo
        
      end Subroutine Wrapgrup_Forces

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>   This routine is called after a  Langevin or HMC step.  On exit, the storage is full  with 
!>   ledt propagationsm,  the Green function is on time slice 0, and  both  
!>   udvl, udvr are on time slice 0. 
!--------------------------------------------------------------------
      Subroutine Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)

        Implicit none
        
        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
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
        end Interface
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR 
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
                  
        ! Local
        Integer :: NSTM, nf,  nt, nt1,  NST, NVAR
        Complex (Kind=Kind(0.d0)) :: Z

        
        NSTM = Size(Stab_nt,1) - 1 
        Do nf = 1,N_FL
           if (Projector) then
              CALL udvl(nf)%reset('l',WF_L(nf)%P)
              CALL udvst(NSTM, nf)%reset('l',WF_L(nf)%P)
           else
              CALL udvl(nf)%reset('l')
              CALL udvst(NSTM, nf)%reset('l')
           endif
        ENDDO


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
        
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        ENDDO
        
        NVAR = 1
        Phase = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf = 1,N_Fl
           CALL CGR(Z, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
           Phase = Phase*Z
        Enddo
        call Op_phase(Phase,OP_V,Nsigma,N_SUN)

      end Subroutine Langevin_HMC_Reset_storage
      
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles a  Langevin sweep.
!>   On input: a) GR is on the first time slice and  the storage is full with
!>                left propagations.   Udvr  and Udvl are on time slice 1.
!>             b) If L_Forces = .T. (.F.) Fermion_Forces  are (not)  provided      
!>                If L_Forces = .F. (.T.) equal time measurements are (not)  carried  out. 
!>   On output: The  field configuration is  updated.  GR, Udvr,  Udvl and Udvst are as on input but with the
!>              updated configuration.  
!> 
!--------------------------------------------------------------------

      SUBROUTINE  Langevin_HMC_update(this,Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN, LTAU)
        
        Implicit none
        
        class (Langevin_HMC_type) :: this
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:)   :: Test
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR, GR_Tilde
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
        Integer, intent(in) :: LOBS_ST, LOBS_EN, LTAU

        !Local
        Integer                   :: N_op, n, nt
        Real    (Kind=Kind(0.d0)) :: X, Xmax
        Logical                   :: Calc_Obser_eq
        

        select case (trim(this%Update_scheme))
        case("Langevin")
           Calc_Obser_eq = .True.
           If (LTAU == 1)   Calc_Obser_eq = .false.
           If ( .not. this%L_Forces) &
                &  Call Langevin_HMC_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                &  LOBS_ST, LOBS_EN, Calc_Obser_eq )
           
           Call Control_Langevin   ( this%Forces,Group_Comm )
           
           Call Ham_Langevin_HMC_S0( this%Forces_0)
           
           N_op = size(nsigma%f,1)
           !  Determine running time step
           Xmax = 0.d0
           do n = 1,N_op
              do nt = 1,Ltrot
                 X = abs(Real(this%Forces  (n,nt), Kind(0.d0)))
                 if (X > Xmax) Xmax = X
                 X = abs(Real(this%Forces_0(n,nt), Kind(0.d0)))
                 if (X > Xmax) Xmax = X
              enddo
           enddo
           this%Delta_t_running = this%Delta_t_Langevin_HMC 
           If ( Xmax >  this%Max_Force ) this%Delta_t_running = this%Max_Force &
                &                              * this%Delta_t_Langevin_HMC / Xmax
           
           
           do n = 1,N_op
              if (OP_V(n,1)%type == 3 ) then
                 do nt = 1,Ltrot
                    nsigma%f(n,nt)   = nsigma%f(n,nt)  -  ( this%Forces_0(n,nt) +  &
                         &  real( Phase*this%Forces(n,nt),kind(0.d0)) / Real(Phase,kind(0.d0)) ) * this%Delta_t_running + &
                         &  sqrt( 2.d0 * this%Delta_t_running) * rang_wrap()
                 enddo
              endif
           enddo
           Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
           this%L_Forces = .False. 
        case("HMC")
           WRITE(error_unit,*) 'HMC  step is not yet implemented'
           error stop 1
        case default
           WRITE(error_unit,*) 'Unknown Global_update_scheme ', trim(this%Update_scheme) 
           WRITE(error_unit,*) 'Global_update_scheme is Langevin or HMC'
           error stop 1
        end select
           
      end SUBROUTINE Langevin_HMC_update
     
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Allocates space for Forces 
!>       Checks that all fields are of tpye 3
!>       Sets default running time step 
!--------------------------------------------------------------------

      
      SUBROUTINE  Langevin_HMC_setup(this,Langevin,HMC, Delta_t_Langevin_HMC, Max_Force, Leapfrog_steps )
#ifdef MPI
        Use mpi
#endif

        Implicit none

        
        Integer :: Nr,Nt, I

        class (Langevin_HMC_type) :: this

        Logical                  , Intent(in)   :: Langevin, HMC
        Integer                  , Intent(in)   :: Leapfrog_steps
        Real    (Kind=Kind(0.d0)), Intent(in)   :: Delta_t_Langevin_HMC, Max_Force

        !Local
        Integer ::  IERR
        Logical ::  lexist
#ifdef MPI
        Real (Kind=Kind(0.d0)) :: X
        INTEGER                :: STATUS(MPI_STATUS_SIZE), ISIZE, IRANK
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
        

        If (Langevin) then 
           !  Check that all  fields are of type 3
           Nr = size(nsigma%f,1)
           Nt = size(nsigma%f,2)
           Do i = 1, Nr
              if ( nsigma%t(i) /= 3 ) then
                 WRITE(error_unit,*) 'For the Langevin runs, all fields have to be of type 3'
                 error stop 1
              endif
           enddo
           Allocate ( this%Forces(Nr,Nt),  this%Forces_0(Nr,Nt) )
           this%Update_scheme        =  "Langevin"
           this%Delta_t_Langevin_HMC =  Delta_t_Langevin_HMC
           this%Max_Force            =  Max_Force
           this%L_Forces             = .False.

           inquire (file="Langevin_time_steps",exist=lexist)
           if (lexist) then
#if defined(MPI)       
              IF (IRANK == 0) THEN
                 OPEN(UNIT=10,FILE="Langevin_time_steps",STATUS='OLD',ACTION='READ',IOSTAT=IERR)
                 Read(10,*) this%Delta_t_running
                 DO I = 1,ISIZE-1
                    Read (10,*) X
                    CALL MPI_SEND(X,1,MPI_REAL8, I, I+1024, MPI_COMM_WORLD,IERR)
                 ENDDO
                 Close(10)
              ELSE
                 CALL MPI_RECV(X, 1, MPI_REAL8,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
                 this%Delta_t_running = X
              ENDIF
#else
              OPEN(UNIT=10,FILE="Langevin_time_steps",STATUS='OLD',ACTION='READ',IOSTAT=IERR)
              Read(10,*) this%Delta_t_running
              Close(10)
#endif
           else
              this%Delta_t_running      =  Delta_t_Langevin_HMC
           endif
        elseif (HMC) then
           WRITE(error_unit,*) 'HMC  step is not yet implemented'
           error stop 1
        else
           this%Update_scheme        =  "None"
        endif
        
      end SUBROUTINE Langevin_HMC_setup

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Deallocates space for forces  and prints out running time step
!--------------------------------------------------------------------

      SUBROUTINE  Langevin_HMC_clear(this) 

#ifdef MPI
        Use mpi
#endif

        Implicit none

        class (Langevin_HMC_type) :: this

        !Local
        Integer :: IERR

#ifdef MPI
        Real (Kind=Kind(0.d0)) :: X
        INTEGER                :: I
        INTEGER                :: STATUS(MPI_STATUS_SIZE), ISIZE, IRANK
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

        select case (trim(this%Update_scheme))
        case("Langevin")
           
#if defined(MPI)       
           IF (IRANK .ne. 0) THEN
              CALL MPI_SEND(this%Delta_t_running,1,MPI_REAL8, 0, Irank + 1024, MPI_COMM_WORLD,IERR)
           ELSE
              OPEN(UNIT=10,FILE="Langevin_time_steps",STATUS='Unknown',IOSTAT=IERR)
              Write(10,*) this%Delta_t_running
              Do I = 1, Isize-1
                 CALL MPI_RECV(X, 1, MPI_REAL8,I,  I + 1024,  MPI_COMM_WORLD,STATUS,IERR)
                 Write(10,*) X
              enddo
              Close(10)
           ENDIF
#else
           OPEN(UNIT=10,FILE="Langevin_time_steps",STATUS='Unknown',IOSTAT=IERR)
           Write(10,*) this%Delta_t_running
           Close(10)
#endif
           Deallocate ( Langevin_HMC%Forces, Langevin_HMC%Forces_0 )
        case ("HMC")
           WRITE(error_unit,*) 'HMC  step is not yet implemented'
           error stop 1
        case default
        end select
      end SUBROUTINE Langevin_HMC_clear

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       sets L_Forces
!--------------------------------------------------------------------
      Subroutine Langevin_HMC_set_L_Forces(this, L_Forces)
        Implicit none
        
        class (Langevin_HMC_type) :: this
        Logical, intent(in) :: L_Forces
        Langevin_HMC%L_Forces = L_Forces
      end Subroutine Langevin_HMC_set_L_Forces

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Returns Update_scheme
!--------------------------------------------------------------------
      function Langevin_HMC_get_Update_scheme(this)
        Implicit none
        
        class (Langevin_HMC_type) :: this
        
        Character (Len=64) :: Langevin_HMC_get_Update_scheme
        
        Langevin_HMC_get_Update_scheme =  this%Update_scheme
        
      end function Langevin_HMC_get_Update_scheme

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Sets the update_scheme
!--------------------------------------------------------------------
      subroutine Langevin_HMC_set_Update_scheme(this, Langevin, HMC )
        Implicit none
        
        class (Langevin_HMC_type) :: this
        
        Logical, intent(in) :: Langevin, HMC

        If (Langevin) then
           this%Update_scheme        =  "Langevin"
        elseif (HMC)  then
           this%Update_scheme        =  "HMC"
        else
           this%Update_scheme        =  "None"
        endif
           

      end subroutine Langevin_HMC_set_Update_scheme
      
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Returns Delta_t_running
!--------------------------------------------------------------------
      function Langevin_HMC_get_Delta_t_running(this)
        Implicit none
        
        class (Langevin_HMC_type) :: this
        Real(Kind=Kind(0.d0)) :: Langevin_HMC_get_Delta_t_running
        Langevin_HMC_get_Delta_t_running = Langevin_HMC%Delta_t_running
      end function Langevin_HMC_get_Delta_t_running

    end Module Langevin_HMC_mod
