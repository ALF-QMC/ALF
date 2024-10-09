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

module Langevin_HMC_mod

   use runtime_error_mod
   use Hamiltonian_main
   use UDV_State_mod
   use Control
   use Hop_mod
   use wrapur_mod
   use wrapul_mod
   use cgr1_mod
   use iso_fortran_env, only: output_unit, error_unit
#ifdef MPI
   use mpi
#endif

   implicit none

   private

   public :: Langevin_HMC, Langevin_HMC_type, Langevin_HMC_Reset_storage

   enum, bind(c)
      enumerator :: Scheme_none = 0
      enumerator :: Scheme_Langevin = 1
      enumerator :: Scheme_HMC = 2
   end enum
   type Langevin_HMC_type
      private
      integer                                 :: scheme
      character(Len=64)                      :: Update_scheme
      logical                                 :: L_Forces
      real(Kind=kind(0.d0))               :: Delta_t_running, Delta_t_Langevin_HMC, Max_Force
      integer                                 :: Leapfrog_Steps
      real(Kind=kind(0.d0)), allocatable  :: Det_vec_old(:, :)
      complex(Kind=kind(0.d0)), allocatable  :: Phase_Det_old(:)
      complex(Kind=kind(0.d0)), allocatable  :: Forces(:, :)

      real(Kind=kind(0.d0)), allocatable  :: Forces_0(:, :)
   contains
      procedure  ::    make => Langevin_HMC_setup
      procedure  ::    clean => Langevin_HMC_clear
      procedure  ::    Wrap_Forces => Wrapgrup_Forces
      procedure  ::    Update => Langevin_HMC_update
      procedure  ::    set_L_Forces => Langevin_HMC_set_L_Forces
      procedure  ::    set_det_OLD => Langevin_HMC_set_det_old
      procedure  ::    get_Update_scheme => Langevin_HMC_get_Update_scheme
      procedure  ::    set_Update_scheme => Langevin_HMC_set_Update_scheme
      procedure  ::    get_Delta_t_running => Langevin_HMC_get_Delta_t_running
      procedure  ::    calc_Forces => Langevin_HMC_Forces
   end type Langevin_HMC_type

   type(Langevin_HMC_type) :: Langevin_HMC

contains

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

   subroutine Langevin_HMC_Forces(this, Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN, Calc_Obser_eq)
      implicit none

      class(Langevin_HMC_type) :: this
      class(UDV_State), intent(inout), allocatable, dimension(:) :: udvl, udvr
      class(UDV_State), intent(in), allocatable, dimension(:, :)    :: udvst
      complex(Kind=kind(0.d0)), intent(inout)                     :: Phase
      complex(Kind=kind(0.d0)), intent(inout), allocatable, dimension(:, :)   :: Test
      complex(Kind=kind(0.d0)), intent(inout), allocatable, dimension(:, :, :) :: GR, GR_Tilde
      integer, intent(in), dimension(:), allocatable :: Stab_nt
      integer, intent(in) :: LOBS_ST, LOBS_EN
      logical, intent(in) :: Calc_Obser_eq

      !Local
      integer :: NSTM, n, nf, nf_eff, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J
      complex(Kind=kind(0.d0)) :: Z, Z1, Phase_array(N_FL)
      real(Kind=kind(0.d0)) :: spin

      NSTM = size(Stab_nt, 1) - 1
      !Do  n = 0,NSTM
      !   Write(6,*)  n, Stab_nt(n)
      !Enddo

      this%Forces = cmplx(0.d0, 0.d0, kind(0.d0))
      do nf_eff = 1, N_FL_eff
         if (Projector) then
            call udvr(nf_eff)%reset('r', WF_R(nf_eff)%P)
         else
            call udvr(nf_eff)%reset('r')
         end if
      end do
      NST = 1
      ! Check: how does wrap_forces work for forces that are zero?
      do NTAU = 0, LTROT - 1
         NTAU1 = NTAU + 1

         call this%Wrap_Forces(Gr, ntau1)

         if (NTAU1 == Stab_nt(NST)) then
            NT1 = Stab_nt(NST - 1)
            call WRAPUR(NT1, NTAU1, udvr)
            Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
            do nf_eff = 1, N_FL_eff
               nf = Calc_FL_map(nf_eff)
               ! Read from storage left propagation from LTROT to  NTAU1
               udvl(nf_eff) = udvst(NST, nf_eff)
               NVAR = 1
               if (NTAU1 .gt. LTROT/2) NVAR = 2
               TEST(:, :) = GR(:, :, nf)
               call CGR(Z1, NVAR, GR(:, :, nf), UDVR(nf_eff), UDVL(nf_eff))
               call Control_PrecisionG(GR(:, :, nf), Test, Ndim)
               call Op_phase(Z1, OP_V, Nsigma, nf)
               Phase_array(nf) = Z1
            end do
            if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
            Z = product(Phase_array)
            Z = Z**N_SUN
            call Control_PrecisionP(Z, Phase)
            Phase = Z
            NST = NST + 1
         end if

         if (NTAU1 .ge. LOBS_ST .and. NTAU1 .le. LOBS_EN .and. Calc_Obser_eq) then
            if (Symm) then
               call Hop_mod_Symm(GR_Tilde, GR, ntau1)
               if (reconstruction_needed) call ham%GR_reconstruction(GR_Tilde)
               call ham%Obser(GR_Tilde, PHASE, Ntau1, this%Delta_t_running)
            else
               if (reconstruction_needed) call ham%GR_reconstruction(GR)
               call ham%Obser(GR, PHASE, Ntau1, this%Delta_t_running)
            end if
         end if
      end do

   end subroutine Langevin_HMC_Forces

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

   subroutine Wrapgrup_Forces(this, Gr, NT1)

      implicit none

      class(Langevin_HMC_type) :: this
      complex(Kind=kind(0.d0)), intent(inout), dimension(:, :, :) :: Gr
      integer, intent(in)                                        :: nt1

      !Local
      complex(Kind=kind(0.d0)) :: Z(N_FL), Z1, g_loc
      integer ::  nf, I, J, n, N_type, nf_eff
      complex(Kind=kind(0.d0)) :: spin

      do nf_eff = 1, N_FL_eff
         nf = Calc_FL_map(nf_eff)
         call HOP_MOD_mmthr(GR(:, :, nf), nf, nt1)
         call HOP_MOD_mmthl_m1(GR(:, :, nf), nf, nt1)
      end do
      do n = 1, size(OP_V, 1)
         this%Forces(n, nt1) = cmplx(0.d0, 0.d0, kind(0.d0))
         do nf_eff = 1, N_FL_eff
            nf = Calc_FL_map(nf_eff)
            spin = nsigma%f(n, nt1) ! Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
            N_type = 1
            call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), spin, Ndim, N_Type, nt1)
            N_type = 2
            call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), spin, Ndim, N_Type, nt1)
         end do
         !CHECK: Are we sure that only Forces representing continuous fields are finite?
         !CHECK: Is it good enough if the forces on discrete fields are zero such as they do not get updated during HMC?
         if (OP_V(n, 1)%type == 3) then
            Z = cmplx(0.d0, 0.d0, kind(0.d0))
            do nf_eff = 1, N_Fl_eff
               nf = Calc_FL_map(nf_eff)
               do I = 1, size(OP_V(n, nf)%P, 1)
                  do J = 1, size(OP_V(n, nf)%P, 1)
                     Z1 = cmplx(0.d0, 0.d0, kind(0.d0))
                     if (I == J) Z1 = cmplx(1.d0, 0.d0, kind(0.d0))
                     Z(nf) = Z(nf) + Op_V(n, nf)%O(I, J)*(Z1 - Gr(Op_V(n, nf)%P(J), Op_V(n, nf)%P(I), nf))
                  end do
               end do
               ! We now include alpha here, so far, alpha was not included explicitly and supposed to be part of Forces0
               Z(nf) = Z(nf) + Op_V(n, nf)%alpha
            end do
            if (reconstruction_needed) call ham%weight_reconstruction(Z)
            do nf = 1, N_Fl
               g_loc = Op_V(n, nf)%g
               if (Op_V(n, nf)%get_g_t_alloc()) g_loc = Op_V(n, nf)%g_t(nt1)
               this%Forces(n, nt1) = this%Forces(n, nt1) - &
                    &    g_loc*Z(nf)*cmplx(real(N_SUN, kind(0.d0)), 0.d0, kind(0.d0))
            end do
         end if
      end do

   end subroutine Wrapgrup_Forces

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!>   This routine is called after a  Langevin or HMC step.  On exit, the storage is full  with
!>   ledt propagationsm,  the Green function is on time slice 0, and  both
!>   udvl, udvr are on time slice 0.
!--------------------------------------------------------------------
   subroutine Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)

      implicit none

      class(UDV_State), intent(inout), allocatable, dimension(:) :: udvl, udvr
      class(UDV_State), intent(inout), allocatable, dimension(:, :) :: udvst
      complex(Kind=kind(0.d0)), intent(inout) :: Phase
      complex(Kind=kind(0.d0)), intent(inout), allocatable, dimension(:, :, :) :: GR
      integer, intent(in), dimension(:), allocatable :: Stab_nt

      ! Local
      integer :: NSTM, nf, nt, nt1, NST, NVAR, nf_eff
      complex(Kind=kind(0.d0)) :: Z, Phase_array(N_FL)

      NSTM = size(Stab_nt, 1) - 1
      do nf_eff = 1, N_FL_eff
         nf = Calc_FL_map(nf_eff)
         if (Projector) then
            call udvl(nf_eff)%reset('l', WF_L(nf)%P)
            call udvst(NSTM, nf_eff)%reset('l', WF_L(nf)%P)
         else
            call udvl(nf_eff)%reset('l')
            call udvst(NSTM, nf_eff)%reset('l')
         end if
      end do

      do NST = NSTM - 1, 1, -1
         NT1 = Stab_nt(NST + 1)
         NT = Stab_nt(NST)
         !Write(6,*)'Hi', NT1,NT, NST
         call WRAPUL(NT1, NT, UDVL)
         do nf_eff = 1, N_FL_eff
            UDVST(NST, nf_eff) = UDVL(nf_eff)
         end do
      end do
      NT1 = stab_nt(1)
      call WRAPUL(NT1, 0, UDVL)

      do nf_eff = 1, N_FL_eff
         nf = Calc_FL_map(nf_eff)
         if (Projector) then
            call udvr(nf_eff)%reset('r', WF_R(nf)%P)
         else
            call udvr(nf_eff)%reset('r')
         end if
      end do

      NVAR = 1
      Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
      do nf_eff = 1, N_Fl_eff
         nf = Calc_FL_map(nf_eff)
         call CGR(Z, NVAR, GR(:, :, nf), UDVR(nf_eff), UDVL(nf_eff))
         call Op_phase(Z, OP_V, Nsigma, nf)
         Phase_array(nf) = Z
      end do
      if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
      Phase = product(Phase_array)
      Phase = Phase**N_SUN

   end subroutine Langevin_HMC_Reset_storage

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

   subroutine Langevin_HMC_update(this, Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN, LTAU)
      use Global_mod
      use Random_Wrap

      implicit none

      class(Langevin_HMC_type) :: this
      class(UDV_State), intent(inout), allocatable, dimension(:) :: udvl, udvr
      class(UDV_State), intent(inout), allocatable, dimension(:, :) :: udvst
      complex(Kind=kind(0.d0)), intent(inout) :: Phase
      complex(Kind=kind(0.d0)), intent(inout), allocatable, dimension(:, :)   :: Test
      complex(Kind=kind(0.d0)), intent(inout), allocatable, dimension(:, :, :) :: GR, GR_Tilde
      integer, intent(in), dimension(:), allocatable :: Stab_nt
      integer, intent(in) :: LOBS_ST, LOBS_EN, LTAU

      !Local
      integer                   :: N_op, n, nt, n1, n2, i, j, t_leap, nf
      real(Kind=kind(0.d0)) :: X, Xmax, E_kin_old, E_kin_new, T0_Proposal_ratio, weight, cluster_size
      logical                   :: Calc_Obser_eq, toggle
      real(Kind=kind(0.d0)), allocatable :: Det_vec_old(:, :), Det_vec_new(:, :)
      complex(Kind=kind(0.d0)), allocatable :: Phase_Det_new(:), Phase_Det_old(:)
      real(Kind=kind(0.d0)), allocatable :: p_tilde(:, :)
      type(Fields)           :: nsigma_old
      character(Len=64)         :: storage
      complex(Kind=kind(0.d0))  :: Ratio(2), Phase_old, Ratiotot, Phase_new, Z

      select case (this%scheme) !(trim(this%Update_scheme))
      case (Scheme_Langevin) !("Langevin")
         Calc_Obser_eq = .true.
         if (LTAU == 1) Calc_Obser_eq = .false.
         if (.not. this%L_Forces) &
              &  call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
              &  LOBS_ST, LOBS_EN, Calc_Obser_eq)

         call Control_Langevin(this%Forces, Group_Comm)

         call ham%Ham_Langevin_HMC_S0(this%Forces_0)
         ! check: should this be going over all the operators?
         N_op = size(nsigma%f, 1)
         !  Determine running time step
         Xmax = 0.d0
         do n = 1, N_op
            do nt = 1, Ltrot
               X = abs(real(this%Forces(n, nt), kind(0.d0)))
               if (X > Xmax) Xmax = X
               X = abs(real(this%Forces_0(n, nt), kind(0.d0)))
               if (X > Xmax) Xmax = X
            end do
         end do
         this%Delta_t_running = this%Delta_t_Langevin_HMC
         if (Xmax > this%Max_Force) this%Delta_t_running = this%Max_Force &
                                                          &                              *this%Delta_t_Langevin_HMC/Xmax

         ! check: I think this already ensures only continuous fields are updated
         do n = 1, N_op
            if (OP_V(n, 1)%type == 3) then
               do nt = 1, Ltrot
                  nsigma%f(n, nt) = nsigma%f(n, nt) - (this%Forces_0(n, nt) +  &
                       &  real(Phase*this%Forces(n, nt), kind(0.d0))/real(Phase, kind(0.d0)))*this%Delta_t_running + &
                       &  sqrt(2.d0*this%Delta_t_running)*rang_wrap()
               end do
            end if
         end do
         call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
         this%L_Forces = .false.
      case (Scheme_HMC) !("HMC")
         !Calc det[phi] if not stored (used after Leapfrog for acceptance/rejectance step)
         storage = "Full"
         allocate (Phase_Det_new(N_FL), Det_vec_new(NDIM, N_FL))
         allocate (Phase_Det_old(N_FL), Det_vec_old(NDIM, N_FL))
         if (.not. this%L_Forces) then
            call Compute_Fermion_Det(Phase_det_old, Det_Vec_old, udvl, udvst, Stab_nt, storage)
         else
            Det_vec_old = this%Det_vec_old
            Phase_det_old = this%Phase_det_old
         end if

         !store old phi
         n1 = size(nsigma%f, 1)
         n2 = size(nsigma%f, 2)
         call nsigma_old%make(n1, n2)
         nsigma_old%f = nsigma%f
         nsigma_old%t = nsigma%t
         Phase_old = Phase

         !Calc Forces phi (del H / del phi) (Calc_Obser_eq always false everywhere)
         Calc_Obser_eq = .false.
         if (.not. this%L_Forces) then
            call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
            &  LOBS_ST, LOBS_EN, Calc_Obser_eq)
         end if
         !Draw p and store E_kin_old
         allocate (p_tilde(n1, n2))
         E_kin_old = 0.0d0
         do j = 1, n2
            do i = 1, n1
               p_tilde(i, j) = rang_wrap()
               E_kin_old = E_kin_old + 0.5*p_tilde(i, j)**2
            end do
         end do

         !Apply B to Forces from phi
         call ham%Ham_Langevin_HMC_S0(this%Forces_0)
         this%Forces_0 = this%Forces_0 + real(Phase*this%Forces, kind(0.d0))/real(Phase, kind(0.d0))
         call ham%Apply_B_HMC(this%Forces_0, .false.)

         !Do half step update of p
         p_tilde = p_tilde - 0.5*this%Delta_t_Langevin_HMC*this%Forces_0
! #define HMC_invertibility
#ifdef HMC_invertibility
         write (1000, *) nsigma%f
#endif

         leap_frog_bulk = .true.
         !Start Leapfrog loop (Leapfrog_steps)
         do t_leap = 1, this%Leapfrog_Steps
            ! update phi by delta t
            this%Forces_0 = p_tilde
            call ham%Apply_B_HMC(this%Forces_0, .true.)
            nsigma%f = nsigma%f + this%Delta_t_Langevin_HMC*this%Forces_0
#ifdef HMC_invertibility
            write (1000 + t_leap, *) nsigma%f
#endif

            ! reset storage
            if (t_leap == this%Leapfrog_Steps) leap_frog_bulk = .false.
            call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
            ! if last step
            X = 1.0d0
            if (t_leap == this%Leapfrog_Steps) then
               ! calc det[phi']
               call Compute_Fermion_Det(Phase_det_new, Det_Vec_new, udvl, udvst, Stab_nt, storage)
               ! LATER (optimization idea): save Phase, GR, udvr, udvl (move calc det out of the loop)
               ! reduce delta t by 1/2
               X = 0.5d0
            end if
            ! calc forces (Calc_Obser_eq false)
            call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                               &  LOBS_ST, LOBS_EN, Calc_Obser_eq)
            !Apply B to Forces from phi
            call ham%Ham_Langevin_HMC_S0(this%Forces_0)
            this%Forces_0 = this%Forces_0 + real(Phase*this%Forces, kind(0.d0))/real(Phase, kind(0.d0))
            call ham%Apply_B_HMC(this%Forces_0, .false.)

            ! update p by delta t UNLESS last step, then delta t / 2
            p_tilde = p_tilde - X*this%Delta_t_Langevin_HMC*this%Forces_0
         end do

#ifdef HMC_invertibility
         p_tilde = -p_tilde
         !Do half step update of p
         p_tilde = p_tilde - 0.5*this%Delta_t_Langevin_HMC*this%Forces_0
         write (2000 + this%Leapfrog_Steps, *) nsigma%f

         !Start Leapfrog loop (Leapfrog_steps)
         do t_leap = 1, this%Leapfrog_Steps
            ! update phi by delta t
            this%Forces_0 = p_tilde
            call ham%Apply_B_HMC(this%Forces_0, .true.)
            nsigma%f = nsigma%f + this%Delta_t_Langevin_HMC*this%Forces_0
            write (2000 + this%Leapfrog_Steps - t_leap, *) nsigma%f

            ! reset storage
            call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
            ! if last step
            X = 1.0d0
            if (t_leap == this%Leapfrog_Steps) then
               ! calc det[phi']
               call Compute_Fermion_Det(Phase_det_new, Det_Vec_new, udvl, udvst, Stab_nt, storage)
               ! LATER (optimization idea): save Phase, GR, udvr, udvl (move calc det out of the loop)
               ! reduce delta t by 1/2
               X = 0.5d0
            end if
            ! calc forces (Calc_Obser_eq false)
            call this%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                               &  LOBS_ST, LOBS_EN, Calc_Obser_eq)
            !Apply B to Forces from phi
            call ham%Ham_Langevin_HMC_S0(this%Forces_0)
            this%Forces_0 = this%Forces_0 + real(Phase*this%Forces, kind(0.d0))/real(Phase, kind(0.d0))
            call ham%Apply_B_HMC(this%Forces_0, .false.)

            ! update p by delta t UNLESS last step, then delta t / 2
            p_tilde = p_tilde - X*this%Delta_t_Langevin_HMC*this%Forces_0
         end do
#endif
         !(LATER (optimization idea) restore Phase, GR, udvr, udvl if calc det moved here)
         !calc ratio
         E_kin_new = 0.0d0
         do j = 1, n2
            do i = 1, n1
               E_kin_new = E_kin_new + 0.5*p_tilde(i, j)**2
            end do
         end do
         T0_Proposal_ratio = exp(-E_kin_new + E_kin_old)
         Ratiotot = Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
              &                          Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio, Ratio)
         Weight = abs(real(Phase_old*Ratiotot, kind=kind(0.d0))/real(Phase_old, kind=kind(0.d0)))

         Phase_new = cmplx(1.d0, 0.d0, kind=kind(0.d0))
         do nf = 1, N_Fl
            Phase_new = Phase_new*Phase_det_new(nf)
         end do
         call Op_phase(Phase_new, OP_V, Nsigma, N_SUN)

         TOGGLE = .false.
         if (Weight > ranf_wrap()) then
            TOGGLE = .true.
            ! store det[phi']
            this%Det_vec_old = Det_vec_new
            this%Phase_Det_old = Phase_det_new
            Phase = Phase_new
         else
            ! store det[phi]
            this%Det_vec_old = Det_vec_old
            this%Phase_Det_old = Phase_det_old
            ! restore phi
            nsigma%t = nsigma_old%t
            nsigma%f = nsigma_old%f
            Phase = Phase_old
         end if

         ! Keep track of the HMC acceptance rate and precision
         Z = Phase_old*Ratiotot/abs(Ratiotot)

         call Control_PrecisionP_HMC(Z, Phase_new)
         cluster_size = dble(size(nsigma%f))
         call Control_upgrade_HMC(TOGGLE)

         ! reset storage
         call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
         !if accepted
         ! LATER (optimization idea) restore Phase, GR, udvr, udvl and don't reset storage
      case default
         write (error_unit, *) 'Unknown Global_update_scheme ', trim(this%Update_scheme)
         write (error_unit, *) 'Global_update_scheme is Langevin or HMC'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end select

   end subroutine Langevin_HMC_update

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Allocates space for Forces
!>       Checks that all fields are of tpye 3
!>       Sets default running time step
!--------------------------------------------------------------------

   subroutine Langevin_HMC_setup(this, Langevin, HMC, Delta_t_Langevin_HMC, Max_Force, Leapfrog_steps)

      implicit none

      integer :: Nr, Nt, I

      class(Langevin_HMC_type) :: this

      logical, intent(in)   :: Langevin, HMC
      integer, intent(in)   :: Leapfrog_steps
      real(Kind=kind(0.d0)), intent(in)   :: Delta_t_Langevin_HMC, Max_Force

      !Local
      integer ::  IERR
      logical ::  lexist
#ifdef MPI
      real(Kind=kind(0.d0)) :: X
      integer                :: STATUS(MPI_STATUS_SIZE), ISIZE, IRANK
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
#endif

      if (Langevin) then
         !  Check that all  fields are of type 3
         Nr = size(nsigma%f, 1)
         Nt = size(nsigma%f, 2)
         do i = 1, Nr
            if (nsigma%t(i) /= 3) then
               write (error_unit, *) 'For the Langevin runs, all fields have to be of type 3'
               call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
            end if
         end do
         allocate (this%Forces(Nr, Nt), this%Forces_0(Nr, Nt))
         this%Update_scheme = "Langevin"
         this%scheme = Scheme_Langevin
         this%Delta_t_Langevin_HMC = Delta_t_Langevin_HMC
         this%Max_Force = Max_Force
         this%L_Forces = .false.

         inquire (file="Langevin_time_steps", exist=lexist)
         if (lexist) then
#if defined(MPI)
            if (IRANK == 0) then
               open (UNIT=10, FILE="Langevin_time_steps", STATUS='OLD', ACTION='READ', IOSTAT=IERR)
               read (10, *) this%Delta_t_running
               do I = 1, ISIZE - 1
                  read (10, *) X
                  call MPI_SEND(X, 1, MPI_REAL8, I, I + 1024, MPI_COMM_WORLD, IERR)
               end do
               close (10)
            else
               call MPI_RECV(X, 1, MPI_REAL8, 0, IRANK + 1024, MPI_COMM_WORLD, STATUS, IERR)
               this%Delta_t_running = X
            end if
#else
            open (UNIT=10, FILE="Langevin_time_steps", STATUS='OLD', ACTION='READ', IOSTAT=IERR)
            read (10, *) this%Delta_t_running
            close (10)
#endif
         else
            this%Delta_t_running = Delta_t_Langevin_HMC
         end if
      elseif (HMC) then
         !  Check that all  fields are of type 3
         !  Check: Shouldn't have to be all type 3 --should be able to update discrete sequentially
         Nr = size(nsigma%f, 1)
         Nt = size(nsigma%f, 2)
         do i = 1, Nr
            if (nsigma%t(i) /= 3) then
               write (error_unit, *) 'For the current HMC runs, all fields have to be of type 3'
               call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
            end if
         end do
         allocate (this%Forces(Nr, Nt), this%Forces_0(Nr, Nt))
         allocate (this%Det_vec_old(NDIM, N_FL), this%Phase_Det_old(N_FL))
         this%Update_scheme = "HMC"
         this%scheme = Scheme_HMC
         this%Delta_t_Langevin_HMC = Delta_t_Langevin_HMC
         this%Max_Force = Max_Force
         this%L_Forces = .false.
         this%Leapfrog_Steps = Leapfrog_steps
         this%Delta_t_running = 1.0d0
         !   WRITE(error_unit,*) 'HMC  step is not yet implemented'
         !   CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      else
         this%Update_scheme = "None"
      end if

   end subroutine Langevin_HMC_setup

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Deallocates space for forces  and prints out running time step
!--------------------------------------------------------------------

   subroutine Langevin_HMC_clear(this)

      implicit none

      class(Langevin_HMC_type) :: this

      !Local
      integer :: IERR

#ifdef MPI
      real(Kind=kind(0.d0)) :: X
      integer                :: I
      integer                :: STATUS(MPI_STATUS_SIZE), ISIZE, IRANK
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
#endif

      select case (this%scheme) !(trim(this%Update_scheme))
      case (Scheme_Langevin) !("Langevin")

#if defined(MPI)
         if (IRANK .ne. 0) then
            call MPI_SEND(this%Delta_t_running, 1, MPI_REAL8, 0, Irank + 1024, MPI_COMM_WORLD, IERR)
         else
            open (UNIT=10, FILE="Langevin_time_steps", STATUS='Unknown', IOSTAT=IERR)
            write (10, *) this%Delta_t_running
            do I = 1, Isize - 1
               call MPI_RECV(X, 1, MPI_REAL8, I, I + 1024, MPI_COMM_WORLD, STATUS, IERR)
               write (10, *) X
            end do
            close (10)
         end if
#else
         open (UNIT=10, FILE="Langevin_time_steps", STATUS='Unknown', IOSTAT=IERR)
         write (10, *) this%Delta_t_running
         close (10)
#endif
         deallocate (this%Forces, this%Forces_0)
      case (Scheme_HMC) !("HMC")
         deallocate (this%Forces, this%Forces_0)
         deallocate (this%Det_vec_old, this%Phase_Det_old)
         !   WRITE(error_unit,*) 'HMC  step is not yet implemented'
         !   CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      case default
      end select
   end subroutine Langevin_HMC_clear

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       sets L_Forces
!--------------------------------------------------------------------
   subroutine Langevin_HMC_set_L_Forces(this, L_Forces)
      implicit none

      class(Langevin_HMC_type) :: this
      logical, intent(in) :: L_Forces
      this%L_Forces = L_Forces
   end subroutine Langevin_HMC_set_L_Forces

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       sets old determinat
!--------------------------------------------------------------------
   subroutine Langevin_HMC_set_det_old(this, Det_vec, Phase_Det)
      implicit none

      class(Langevin_HMC_type) :: this
      real(Kind=kind(0.d0)), allocatable  :: Det_vec(:, :)
      complex(Kind=kind(0.d0)), allocatable  :: Phase_Det(:)

      this%Det_vec_old = Det_vec
      this%Phase_Det_old = Phase_Det
   end subroutine Langevin_HMC_set_det_old

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Returns Update_scheme
!--------------------------------------------------------------------
   function Langevin_HMC_get_Update_scheme(this)
      implicit none

      class(Langevin_HMC_type) :: this

      character(Len=64) :: Langevin_HMC_get_Update_scheme

      Langevin_HMC_get_Update_scheme = this%Update_scheme

   end function Langevin_HMC_get_Update_scheme

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Sets the update_scheme
!--------------------------------------------------------------------
   subroutine Langevin_HMC_set_Update_scheme(this, Langevin, HMC)
      implicit none

      class(Langevin_HMC_type) :: this

      logical, intent(in) :: Langevin, HMC

      if (Langevin) then
         this%Update_scheme = "Langevin"
         this%scheme = Scheme_Langevin
      elseif (HMC) then
         this%Update_scheme = "HMC"
         this%scheme = Scheme_HMC
      else
         this%Update_scheme = "None"
         this%scheme = Scheme_none
      end if

   end subroutine Langevin_HMC_set_Update_scheme

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Returns Delta_t_running
!--------------------------------------------------------------------
   function Langevin_HMC_get_Delta_t_running(this)
      implicit none

      class(Langevin_HMC_type) :: this
      real(Kind=kind(0.d0)) :: Langevin_HMC_get_Delta_t_running
      Langevin_HMC_get_Delta_t_running = this%Delta_t_running
   end function Langevin_HMC_get_Delta_t_running

end module Langevin_HMC_mod
