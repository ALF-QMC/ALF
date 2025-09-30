!  Copyright (C) 2016 - 2022 The ALF project
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


! TODO ATTENTION: this module still has to be updated for flavor symmetries!!!


      Module Langevin_HMC_mod
        
        Use runtime_error_mod 
        Use Hamiltonian_main
        Use UDV_State_mod
        Use Control
        Use Hop_mod
        use wrapur_mod
        use wrapul_mod
        use cgr1_mod
        Use iso_fortran_env, only: output_unit, error_unit
#ifdef MPI
        Use mpi
#endif

        
        Implicit none

        Private
        
        Public :: Langevin_HMC, Langevin_HMC_type, Langevin_HMC_Reset_storage
        
        enum, bind(c)
           enumerator :: Scheme_none = 0
           enumerator :: Scheme_Langevin = 1 
           enumerator :: Scheme_HMC = 2
        end enum
        Type Langevin_HMC_type
           private
           Integer                                 :: scheme
           Character (Len=64)                      :: Update_scheme
           Logical                                 :: L_Forces
           Real    (Kind=Kind(0.d0))               :: Delta_t_running, Delta_t_Langevin_HMC, Max_Force
           Integer                                 :: Leapfrog_Steps
           Real    (Kind=Kind(0.d0)), allocatable  :: Det_vec_old(:,:)
           Complex (Kind=Kind(0.d0)), allocatable  :: Phase_Det_old(:)
           Complex (Kind=Kind(0.d0)), allocatable  :: Forces  (:,:)
           
           Real    (Kind=Kind(0.d0)), allocatable  :: Forces_0(:,:)
         CONTAINS
           procedure  ::    make        => Langevin_HMC_setup
           procedure  ::    clean       => Langevin_HMC_clear
           procedure  ::    Wrap_Forces => Wrapgrup_Forces
           procedure  ::    Update      => Langevin_HMC_update
           procedure  ::    set_L_Forces         => Langevin_HMC_set_L_Forces
           procedure  ::    set_det_OLD          => Langevin_HMC_set_det_old
           procedure  ::    get_Update_scheme    => Langevin_HMC_get_Update_scheme
           procedure  ::    set_Update_scheme    => Langevin_HMC_set_Update_scheme
           procedure  ::    get_Delta_t_running  => Langevin_HMC_get_Delta_t_running
           procedure  ::    calc_Forces          => Langevin_HMC_Forces
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
        
      SUBROUTINE  Langevin_HMC_Forces(this, Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN,Calc_Obser_eq)
        Implicit none
        
        class (Langevin_HMC_type) :: this
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(in), allocatable, dimension(:,:)    :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout)                     :: Phase
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:)   :: Test
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR, GR_Tilde
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
        Integer, intent(in) :: LOBS_ST, LOBS_EN
        Logical, intent(in) :: Calc_Obser_eq 
        

        !Local
        Integer :: NSTM, n, nf, nf_eff, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J
        Complex (Kind=Kind(0.d0)) :: Z, Z1, Phase_array(N_FL)
        Real    (Kind=Kind(0.d0)) :: spin
        
        NSTM = Size(Stab_nt,1) - 1 
        !Do  n = 0,NSTM
        !   Write(6,*)  n, Stab_nt(n)
        !Enddo
        
        this%Forces = cmplx(0.d0,0.d0,Kind(0.d0))
        do nf_eff = 1,N_FL_eff
           if (Projector) then
              CALL udvr(nf_eff)%reset('r',WF_R(nf_eff)%P)
           else
              CALL udvr(nf_eff)%reset('r')
           endif
        Enddo
        NST = 1
        ! Check: how does wrap_forces work for forces that are zero?
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           
           Call  this%Wrap_Forces(Gr, ntau1)
           
           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              CALL WRAPUR(NT1, NTAU1, udvr)
              Phase_array = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf_eff = 1, N_FL_eff
                 nf=Calc_FL_map(nf_eff)
                 ! Read from storage left propagation from LTROT to  NTAU1
                 udvl(nf_eff) = udvst(NST, nf_eff)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf_eff), UDVL(nf_eff))
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
                 call Op_phase(Z1,OP_V,Nsigma,nf) 
                 Phase_array(nf)=Z1
              ENDDO
              if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
              Z=product(Phase_array)
              Z=Z**N_SUN 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST + 1
           ENDIF
           
           IF (NTAU1 .GE. LOBS_ST .AND. NTAU1 .LE. LOBS_EN .and. Calc_Obser_eq ) THEN
              If (Symm) then
                 Call Hop_mod_Symm(GR_Tilde,GR,ntau1)
                 If (reconstruction_needed) Call ham%GR_reconstruction( GR_Tilde )
                 CALL ham%Obser( GR_Tilde, PHASE, Ntau1,this%Delta_t_running )
              else
                 If (reconstruction_needed) Call ham%GR_reconstruction( GR )
                 CALL ham%Obser( GR, PHASE, Ntau1, this%Delta_t_running )
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
        Complex (Kind=Kind(0.d0)) :: Z(N_FL), Z1, g_loc
        Integer ::  nf, I, J, n, N_type, nf_eff
        Complex  (Kind=Kind(0.d0)) :: spin

        Do nf_eff = 1,N_FL_eff
           nf=Calc_FL_map(nf_eff)
           CALL HOP_MOD_mmthr   (GR(:,:,nf), nf, nt1 )
           CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf, nt1 )
        Enddo
        Do n = 1, size(OP_V,1) 
           this%Forces(n,nt1)  = cmplx(0.d0,0.d0,Kind(0.d0))
           Do nf_eff = 1, N_FL_eff
              nf=Calc_FL_map(nf_eff)
              spin = nsigma%f(n,nt1) ! Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
              N_type = 1
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type,nt1)
              N_type =  2
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type,nt1)
           enddo
           !CHECK: Are we sure that only Forces representing continuous fields are finite?
           !CHECK: Is it good enough if the forces on discrete fields are zero such as they do not get updated during HMC?
           if (OP_V(n,1)%type == 3 ) then
              Z = cmplx(0.d0,0.d0,Kind(0.d0))
              Do nf_eff = 1, N_Fl_eff
                 nf=Calc_FL_map(nf_eff)
                 do I = 1,size(OP_V(n,nf)%P,1)
                    do J = 1,size(OP_V(n,nf)%P,1)
                       Z1 =  cmplx(0.d0,0.d0,Kind(0.d0))
                       if ( I == J ) Z1 = cmplx(1.d0,0.d0,Kind(0.d0))
                       Z(nf)  = Z(nf) +    Op_V(n,nf)%O(I,J) * ( Z1 - Gr(Op_V(n,nf)%P(J),Op_V(n,nf)%P(I), nf) )
                    Enddo
                 Enddo
                 ! We now include alpha here, so far, alpha was not included explicitly and supposed to be part of Forces0
                 Z(nf)  = Z(nf) + Op_V(n,nf)%alpha
              Enddo
              if (reconstruction_needed) call ham%weight_reconstruction(Z)
              Do nf = 1, N_Fl
                 g_loc = Op_V(n,nf)%g
                 if (Op_V(n,nf)%get_g_t_alloc() ) g_loc = Op_V(n,nf)%g_t(nt1)
                 this%Forces(n,nt1) =  this%Forces(n,nt1)  - &
                      &    g_loc * Z(nf) *  cmplx(real(N_SUN,Kind(0.d0)), 0.d0, Kind(0.d0))
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
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR 
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
                  
        ! Local
        Integer :: NSTM, nf,  nt, nt1,  NST, NVAR, nf_eff
        Complex (Kind=Kind(0.d0)) :: Z, Phase_array(N_FL)

        
        NSTM = Size(Stab_nt,1) - 1 
        Do nf_eff = 1,N_FL_eff
           nf=Calc_FL_map(nf_eff)
           if (Projector) then
              CALL udvl(nf_eff)%reset('l',WF_L(nf)%P)
              CALL udvst(NSTM, nf_eff)%reset('l',WF_L(nf)%P)
           else
              CALL udvl(nf_eff)%reset('l')
              CALL udvst(NSTM, nf_eff)%reset('l')
           endif
        ENDDO


        DO NST = NSTM-1,1,-1
           NT1 = Stab_nt(NST+1)
           NT  = Stab_nt(NST  )
           !Write(6,*)'Hi', NT1,NT, NST
           CALL WRAPUL(NT1, NT, UDVL)
           Do nf_eff = 1,N_FL_eff
              UDVST(NST, nf_eff) = UDVL(nf_eff)
           ENDDO
        ENDDO
        NT1 = stab_nt(1)
        CALL WRAPUL(NT1, 0, UDVL)
        
        do nf_eff = 1,N_FL_eff
           nf=Calc_FL_map(nf_eff)
           if (Projector) then
              CALL udvr(nf_eff)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf_eff)%reset('r')
           endif
        ENDDO
        
        NVAR = 1
        Phase_array = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf_eff = 1,N_Fl_eff
           nf=Calc_FL_map(nf_eff)
           CALL CGR(Z, NVAR, GR(:,:,nf), UDVR(nf_eff), UDVL(nf_eff))
           call Op_phase(Z,OP_V,Nsigma,nf)
           Phase_array(nf)=Z
        Enddo
        if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
        Phase=product(Phase_array)
        Phase=Phase**N_SUN

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
        use Global_mod
        use Random_Wrap

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
        Integer                   :: N_op, n, nt, n1, n2, i, j, t_leap, nf, nf_eff
        Real    (Kind=Kind(0.d0)) :: X, Xmax,E_kin_old, E_kin_new,T0_Proposal_ratio, weight, cluster_size
        Logical                   :: Calc_Obser_eq, toggle
        Real    (Kind=Kind(0.d0)), allocatable :: Det_vec_old(:,:), Det_vec_new(:,:)
        Complex (Kind=Kind(0.d0)), allocatable :: Phase_Det_new(:), Phase_Det_old(:)
        Real    (Kind=Kind(0.d0)), allocatable :: p_tilde(:,:)
        Type    (Fields)           :: nsigma_old
        Character (Len=64)         :: storage
        Complex (Kind=Kind(0.d0))  :: Ratio(2), Phase_old, Ratiotot,Phase_new, Z
        Complex (Kind=Kind(0.d0))  :: Phase_array(N_FL)

        select case (this%scheme) !(trim(this%Update_scheme))
        case(Scheme_Langevin) !("Langevin")
           Calc_Obser_eq = .True.
           If (LTAU == 1)   Calc_Obser_eq = .false.
           If ( .not. this%L_Forces) &
                &  Call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                &  LOBS_ST, LOBS_EN, Calc_Obser_eq )
           
           Call Control_Langevin   ( this%Forces,Group_Comm )
           
           Call ham%Ham_Langevin_HMC_S0( this%Forces_0)
           ! check: should this be going over all the operators?
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
           
           ! check: I think this already ensures only continuous fields are updated
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
        case(Scheme_HMC) !("HMC")
           !Calc det[phi] if not stored (used after Leapfrog for acceptance/rejectance step)
           storage = "Full"
           Allocate ( Phase_Det_new(N_FL), Det_vec_new(NDIM,N_FL) )
           Allocate ( Phase_Det_old(N_FL), Det_vec_old(NDIM,N_FL) )
           If ( .not. this%L_Forces) then
             Call Compute_Fermion_Det(Phase_det_old,Det_Vec_old, udvl, udvst, Stab_nt, storage)
           else
             Det_vec_old = this%Det_vec_old
             Phase_det_old = this%Phase_det_old
           endif

           !store old phi
           n1 = size(nsigma%f,1)
           n2 = size(nsigma%f,2)
           call nsigma_old%make(n1, n2)
           nsigma_old%f = nsigma%f
           nsigma_old%t = nsigma%t
           Phase_old = Phase

           !Calc Forces phi (del H / del phi) (Calc_Obser_eq always false everywhere)
           Calc_Obser_eq=.False.
           If ( .not. this%L_Forces) then
              Call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
              &  LOBS_ST, LOBS_EN, Calc_Obser_eq )
           endif
           !Draw p and store E_kin_old
           allocate(p_tilde(n1,n2))
           E_kin_old=0.0d0
           Do j=1,n2
              do i=1,n1
                 p_tilde(i,j)=rang_wrap()
                 E_kin_old=E_kin_old + 0.5*p_tilde(i,j)**2
              enddo
           enddo
           
           !Apply B to Forces from phi
           Call ham%Ham_Langevin_HMC_S0( this%Forces_0)
           this%Forces_0 = this%Forces_0 + real( Phase*this%Forces,kind(0.d0)) / Real(Phase,kind(0.d0))
           call ham%Apply_B_HMC(this%Forces_0 ,.false.)
           
           !Do half step update of p
           p_tilde=p_tilde - 0.5*this%Delta_t_Langevin_HMC*this%Forces_0
! #define HMC_invertibility
#ifdef HMC_invertibility
           write(1000,*) nsigma%f
#endif

           leap_frog_bulk=.true.
           !Start Leapfrog loop (Leapfrog_steps)
           Do t_leap=1,this%Leapfrog_Steps
               ! update phi by delta t
               this%Forces_0 = p_tilde
               call ham%Apply_B_HMC(this%Forces_0 ,.True.)
               nsigma%f=nsigma%f + this%Delta_t_Langevin_HMC*this%Forces_0
#ifdef HMC_invertibility
               write(1000+t_leap,*) nsigma%f
#endif

               ! reset storage
               if (t_leap == this%Leapfrog_Steps) leap_frog_bulk=.false.
               Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
               ! if last step
               X=1.0d0
               if (t_leap == this%Leapfrog_Steps) then
                  ! calc det[phi']
                  Call Compute_Fermion_Det(Phase_det_new,Det_Vec_new, udvl, udvst, Stab_nt, storage)
                  ! LATER (optimization idea): save Phase, GR, udvr, udvl (move calc det out of the loop)
                  ! reduce delta t by 1/2
                  X=0.5d0
               endif
               ! calc forces (Calc_Obser_eq false)
               Call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                                  &  LOBS_ST, LOBS_EN, Calc_Obser_eq )
               !Apply B to Forces from phi
               Call ham%Ham_Langevin_HMC_S0( this%Forces_0)
               this%Forces_0 = this%Forces_0 + real( Phase*this%Forces,kind(0.d0)) / Real(Phase,kind(0.d0))
               call ham%Apply_B_HMC(this%Forces_0 ,.false.)

               ! update p by delta t UNLESS last step, then delta t / 2
               p_tilde=p_tilde - X*this%Delta_t_Langevin_HMC*this%Forces_0
           enddo

#ifdef HMC_invertibility
           p_tilde=-p_tilde
           !Do half step update of p
           p_tilde=p_tilde - 0.5*this%Delta_t_Langevin_HMC*this%Forces_0
           write(2000+this%Leapfrog_Steps,*) nsigma%f

           !Start Leapfrog loop (Leapfrog_steps)
           Do t_leap=1,this%Leapfrog_Steps
               ! update phi by delta t
               this%Forces_0 = p_tilde
               call ham%Apply_B_HMC(this%Forces_0 ,.True.)
               nsigma%f=nsigma%f + this%Delta_t_Langevin_HMC*this%Forces_0
               write(2000+this%Leapfrog_Steps-t_leap,*) nsigma%f

               ! reset storage
               Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
               ! if last step
               X=1.0d0
               if (t_leap == this%Leapfrog_Steps) then
                  ! calc det[phi']
                  Call Compute_Fermion_Det(Phase_det_new,Det_Vec_new, udvl, udvst, Stab_nt, storage)
                  ! LATER (optimization idea): save Phase, GR, udvr, udvl (move calc det out of the loop)
                  ! reduce delta t by 1/2
                  X=0.5d0
               endif
               ! calc forces (Calc_Obser_eq false)
               Call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                                  &  LOBS_ST, LOBS_EN, Calc_Obser_eq )
               !Apply B to Forces from phi
               Call ham%Ham_Langevin_HMC_S0( this%Forces_0)
               this%Forces_0 = this%Forces_0 + real( Phase*this%Forces,kind(0.d0)) / Real(Phase,kind(0.d0))
               call ham%Apply_B_HMC(this%Forces_0 ,.false.)

               ! update p by delta t UNLESS last step, then delta t / 2
               p_tilde=p_tilde - X*this%Delta_t_Langevin_HMC*this%Forces_0
           enddo
#endif
           !(LATER (optimization idea) restore Phase, GR, udvr, udvl if calc det moved here)
           !calc ratio
           E_kin_new=0.0d0
           Do j=1,n2
              do i=1,n1
                 E_kin_new=E_kin_new + 0.5*p_tilde(i,j)**2
              enddo
           enddo
           T0_Proposal_ratio=1.0 !exp(-E_kin_new + E_kin_old) ! this could be a exponentially large or small number
           !Compute_Ratio_Global returns Ratio(1)*exp(Ratio(2)) where Ratio(2) contains log(T0_proposal_ratio)
           Ratiotot = Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
                &                          Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio, Ratio)
           Ratiotot = Ratio(1)*exp(Ratio(2) - E_kin_new + E_kin_old)
           Weight = abs(  real( Phase_old * Ratiotot, kind=Kind(0.d0))/real(Phase_old,kind=Kind(0.d0)) )

           Phase_array=1.0d0
           Do nf_eff = 1,N_Fl_eff
              nf=Calc_Fl_map(nf_eff)
              Phase_array(nf) = Phase_det_new(nf)
              Call Op_phase(Phase_array(nf),OP_V,Nsigma,nf)
           Enddo
           if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
           Phase_new=product(Phase_array)
           Phase_new=Phase_new**N_SUN

           TOGGLE = .false.
           if ( Weight > ranf_wrap() )  Then
              TOGGLE = .true.
              ! store det[phi']
              this%Det_vec_old   = Det_vec_new
              this%Phase_Det_old = Phase_det_new
              Phase = Phase_new
           else
              ! store det[phi]
              this%Det_vec_old   = Det_vec_old
              this%Phase_Det_old = Phase_det_old
              ! restore phi
              nsigma%t = nsigma_old%t
              nsigma%f = nsigma_old%f
              Phase = Phase_old
           endif

           ! Keep track of the HMC acceptance rate and precision
           Z = Phase_old * Ratiotot/ABS(Ratiotot)

           Call Control_PrecisionP_HMC(Z,Phase_new)
           cluster_size=dble(size(nsigma%f))
           Call Control_upgrade_HMC(TOGGLE)


           ! reset storage
           Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
           !if accepted 
           ! LATER (optimization idea) restore Phase, GR, udvr, udvl and don't reset storage
        case default
           WRITE(error_unit,*) 'Unknown Global_update_scheme ', trim(this%Update_scheme) 
           WRITE(error_unit,*) 'Global_update_scheme is Langevin or HMC'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
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
                 CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
              endif
           enddo
           Allocate ( this%Forces(Nr,Nt),  this%Forces_0(Nr,Nt) )
           this%Update_scheme        =  "Langevin"
           this%scheme               =  Scheme_Langevin
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
           !  Check that all  fields are of type 3
           !  Check: Shouldn't have to be all type 3 --should be able to update discrete sequentially
           Nr = size(nsigma%f,1)
           Nt = size(nsigma%f,2)
           Do i = 1, Nr
              if ( nsigma%t(i) /= 3 ) then
                 WRITE(error_unit,*) 'For the current HMC runs, all fields have to be of type 3'
                 CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
              endif
           enddo
           Allocate ( this%Forces(Nr,Nt),  this%Forces_0(Nr,Nt) )
           Allocate ( this%Det_vec_old(NDIM,N_FL), this%Phase_Det_old(N_FL) )
           this%Update_scheme        =  "HMC"
           this%scheme               =  Scheme_HMC
           this%Delta_t_Langevin_HMC =  Delta_t_Langevin_HMC
           this%Max_Force            =  Max_Force
           this%L_Forces             = .False.
           this%Leapfrog_Steps       =  Leapfrog_steps
           this%Delta_t_running      =  1.0d0
         !   WRITE(error_unit,*) 'HMC  step is not yet implemented'
         !   CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
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

        select case (this%scheme) !(trim(this%Update_scheme))
        case(Scheme_Langevin) !("Langevin")
           
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
           Deallocate ( this%Forces, this%Forces_0 )
        case(Scheme_HMC) !("HMC")
           Deallocate ( this%Forces, this%Forces_0 )
           Deallocate ( this%Det_vec_old, this%Phase_Det_old )
         !   WRITE(error_unit,*) 'HMC  step is not yet implemented'
         !   CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
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
        this%L_Forces = L_Forces
      end Subroutine Langevin_HMC_set_L_Forces

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       sets old determinat
!--------------------------------------------------------------------
      Subroutine Langevin_HMC_set_det_old(this, Det_vec, Phase_Det)
        Implicit none
        
        class (Langevin_HMC_type) :: this
        Real    (Kind=Kind(0.d0)), allocatable  :: Det_vec(:,:)
        Complex (Kind=Kind(0.d0)), allocatable  :: Phase_Det(:)

        this%Det_vec_old = Det_vec
        this%Phase_Det_old = Phase_Det
      end Subroutine Langevin_HMC_set_det_old

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
           this%scheme               =  Scheme_Langevin
        elseif (HMC)  then
           this%Update_scheme        =  "HMC"
           this%scheme               =  Scheme_HMC
        else
           this%Update_scheme        =  "None"
           this%scheme               =  Scheme_none
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
        Langevin_HMC_get_Delta_t_running = this%Delta_t_running
      end function Langevin_HMC_get_Delta_t_running

    end Module Langevin_HMC_mod
