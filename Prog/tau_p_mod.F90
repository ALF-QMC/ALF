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

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This module handles calculation of imaginary-time-displaced Green functions and
!> calls the routine ObserT.F90 in the Hamiltonian module, so as to compute  user
!> defined time-displaced correlations functions. This module is for the projector code.
!--------------------------------------------------------------------

module Tau_p_mod
   use Hamiltonian_main
   use Operator_mod
   use Control
   use Hop_mod
   use UDV_State_mod
   use Langevin_HMC_mod
   use tau_m_mod  !, only propr, proprm1
   use cgr1_mod, only: cgrp
   use wrapur_mod

contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief      This routine computes the time displaced  zero termperature Green functions and call  obserT.
!> On input:   a) GR,  the equal time Green function,  as  well as udvl, udvr are on time slice
!>                nt_in= stab_nt(nst)   with   stab_nt(NST) <= THTROT+1  and  stab_nt( NST +1 )  > THTROT+1.
!>             b) The storage, udvst, is full with left propagations from  Ltrot to    stab_nt( NST +1 ).
!> On_input    a) GR,  the equal time Green function,  as  well as udvl, udvr are on time slice
!>                nt_in = = Stab_nt(NST_IN)  with nt_in <= THTROT+1.
!>             b) The storage is full with left propagations for all n's with stab_nt(n) > nt_in.
!>             c) udvl and udvr are on time slice nt_in  such that  a call to CGR with udvl and udvr will
!>                produce Gr.
!>
!> On_output   a) The time displaced Green functions have been computed and measurements carried out.
!>             b) If Langevin then 1)  nt_in = 0, 2) forces are computed  3)  time displaced and equal time
!>                observables are  measured.
!--------------------------------------------------------------------

   subroutine Tau_p(udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST_IN, LOBS_ST, LOBS_EN)

      implicit none

      ! Storage is full with U^{<} (left)  propagations.

      integer, intent(In) :: NSTM, NST_IN
      class(UDV_State), dimension(:), allocatable, intent(IN) :: udvl, udvr
      class(UDV_State), dimension(:, :), allocatable, intent(IN) :: udvst
      complex(Kind=kind(0.d0)), intent(in) :: GR(NDIM, NDIM, N_FL), Phase
      integer, intent(In) :: STAB_NT(0:NSTM)
      integer, intent(In) :: LOBS_ST, LOBS_EN

!       Local.
      class(UDV_State), dimension(:), allocatable :: udvr_local
      complex(Kind=kind(0.d0)) :: DETZ, ZK, DET1(2)
      complex(Kind=kind(0.d0)), dimension(:, :, :), allocatable  ::  GRUPB, GRUP
      complex(Kind=kind(0.d0)), dimension(:, :, :), allocatable  ::  G00UP, G0TUP, GT0UP, GTTUP
      complex(Kind=kind(0.d0)), dimension(:, :, :), allocatable  ::  G00UP_T, G0TUP_T, GT0UP_T, GTTUP_T
      complex(Kind=kind(0.d0)), allocatable  :: TEMP(:, :), TMPUP(:, :)

      real(Kind=kind(0.d0))  :: XMEAN_DYN, XMAX_DYN

      integer :: NTAUIN, NTDM, LFAM, NFAM, N_Part, LQ, I, NCON, NF, nf_eff, NFLAG, NL, NT1, NT_ST, NT, NTAU, NTAU1, n

      real(Kind=kind(0.d0)) :: XMEAN, XMAX
      real(Kind=kind(0.d0)) :: Mc_step_weight

      LQ = ndim

      Mc_step_weight = 1.d0
      if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN") Mc_step_weight = Langevin_HMC%get_Delta_t_running()

      allocate (GRUPB(LQ, LQ, N_FL), GRUP(LQ, LQ, N_FL), G00UP(LQ, LQ, N_FL), G0TUP(LQ, LQ, N_FL), &
           &      GT0UP(LQ, LQ, N_FL), GTTUP(LQ, LQ, N_FL), TEMP(LQ, LQ), udvr_local(N_FL_eff))

      if (Symm) then
         allocate (G00UP_T(LQ, LQ, N_FL), G0TUP_T(LQ, LQ, N_FL), GT0UP_T(LQ, LQ, N_FL), GTTUP_T(LQ, LQ, N_FL))
      end if

      do nf_eff = 1, N_FL_eff
         call udvr_local(nf_eff)%alloc(ndim, udvr(nf_eff)%N_part)
         udvr_local(nf_eff) = udvr(nf_eff)
      end do

      GTTUP = GR ! On time slice Stab_nt(NST_IN)
      NT_ST = NST_IN
      do NT = Stab_nt(NT_ST) + 1, Thtrot + 1
         if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
            & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") then
            call Langevin_HMC%Wrap_Forces(GTTUP, NT)
         else
            call PROPRM1(GTTUP, NT)
            call PROPR(GTTUP, NT)
         end if
         if ((str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN"  &
            & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") &
            .and. NT .ge. LOBS_ST .and. NT .le. LOBS_EN) then
            if (Symm) then
               call Hop_mod_Symm(GTTUP_T, GTTUP, nt)
               !call reconstruction of non-calculated flavor blocks
               if (reconstruction_needed) call ham%GR_reconstruction(GTTUP_T)
               call ham%Obser(GTTUP_T, PHASE, NT, Mc_step_weight)
            else
               !call reconstruction of non-calculated flavor blocks
               if (reconstruction_needed) call ham%GR_reconstruction(GTTUP)
               call ham%Obser(GTTUP, PHASE, NT, Mc_step_weight)
            end if
         end if
         if (NT .eq. STAB_NT(NT_ST + 1)) then
            call Wrapur(STAB_NT(NT_ST), STAB_NT(NT_ST + 1), UDVR_local)
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call CGRP(DetZ, GRUP(:, :, nf), udvr_local(nf_eff), udvst(nt_st + 1, nf_eff))
            end do
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call Control_Precision_tau(GTTUP(:, :, nf), GRUP(:, :, nf), Ndim)
            end do
            GTTUP = GRUP
            NT_ST = NT_ST + 1
         end if
      end do

      GRUPB = GTTUP
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do I = 1, Ndim
            GRUPB(I, I, nf) = GRUPB(I, I, nf) - 1.d0
         end do
      end do

      G00UP = GTTUP
      !GTTUP = GTTUP
      GT0UP = GTTUP
      G0TUP = GRUPB
      NTAU = 0
      if (Symm) then
         call Hop_mod_Symm(G00UP_T, G00UP, Thtrot + 1)
         call Hop_mod_Symm(GTTUP_T, GTTUP, Thtrot + 1)
         call Hop_mod_Symm(G0TUP_T, G0TUP, Thtrot + 1)
         call Hop_mod_Symm(GT0UP_T, GT0UP, Thtrot + 1)
         !call reconstruction of non-calculated flavor blocks
         if (reconstruction_needed) then
            call ham%GR_reconstruction(G00UP_T)
            call ham%GR_reconstruction(GTTUP_T)
            call ham%GRT_reconstruction(GT0UP_T, G0TUP_T)
         end if
         call ham%OBSERT(NTAU, GT0UP_T, G0TUP_T, G00UP_T, GTTUP_T, PHASE, Mc_step_Weight)
      else
         !call reconstruction of non-calculated flavor blocks
         if (reconstruction_needed) then
            call ham%GR_reconstruction(G00UP)
            call ham%GR_reconstruction(GTTUP)
            call ham%GRT_reconstruction(GT0UP, G0TUP)
         end if
         call ham%OBSERT(NTAU, GT0UP, G0TUP, G00UP, GTTUP, PHASE, Mc_step_Weight)
      end if
      do NT = THTROT + 1, Ltrot - THTROT
         ! UR is on time slice NT
         NTAU = NT - THTROT - 1
         if (NT .eq. STAB_NT(NT_ST + 1) .and. NTAU /= 0) then
            call Wrapur(STAB_NT(NT_ST), STAB_NT(NT_ST + 1), UDVR_local)
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call CGRP(DetZ, GRUP(:, :, nf), udvr_local(nf_eff), udvst(nt_st + 1, nf_eff))
            end do
            NT_ST = NT_ST + 1
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call Control_Precision_tau(GTTUP(:, :, nf), GRUP(:, :, nf), Ndim)
            end do
            GTTUP = GRUP

            GRUPB = -GRUP
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               do I = 1, Ndim
                  GRUPB(I, I, nf) = GRUPB(I, I, nf) + 1.d0
               end do
            end do
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call MMULT(TEMP, GRUP(:, :, nf), GT0UP(:, :, nf))
               GT0UP(:, :, nf) = TEMP
               call MMULT(TEMP, G0TUP(:, :, nf), GRUPB(:, :, nf))
               G0TUP(:, :, nf) = TEMP
            end do
         end if
         NT1 = NT + 1
         call PROPR(GT0UP, NT1)
         call PROPRM1(G0TUP, NT1)
         if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
          & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") then
            call Langevin_HMC%Wrap_Forces(GTTUP, NT1)
         else
            call PROPRM1(GTTUP, NT1)
            call PROPR(GTTUP, NT1)
         end if

         NTAU1 = NTAU + 1
         if (Symm) then
            call Hop_mod_Symm(G00UP_T, G00UP, Thtrot + 1)
            call Hop_mod_Symm(GTTUP_T, GTTUP, nt1)
            call Hop_mod_Symm(G0TUP_T, G0TUP, Thtrot + 1, nt1)
            call Hop_mod_Symm(GT0UP_T, GT0UP, nt1, Thtrot + 1)
            !call reconstruction of non-calculated flavor blocks
            if (reconstruction_needed) then
               call ham%GR_reconstruction(G00UP_T)
               call ham%GR_reconstruction(GTTUP_T)
               call ham%GRT_reconstruction(GT0UP_T, G0TUP_T)
            end if
            call ham%OBSERT(NTAU1, GT0UP_T, G0TUP_T, G00UP_T, GTTUP_T, PHASE, Mc_step_weight)
            if ((str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
                 & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC")&
                 &.and. NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN) call ham%Obser(GTTUP_T, PHASE, NT1, Mc_step_weight)
         else
            !call reconstruction of non-calculated flavor blocks
            if (reconstruction_needed) then
               call ham%GR_reconstruction(G00UP)
               call ham%GR_reconstruction(GTTUP)
               call ham%GRT_reconstruction(GT0UP, G0TUP)
            end if
            call ham%OBSERT(NTAU1, GT0UP, G0TUP, G00UP, GTTUP, PHASE, Mc_step_weight)
            if ((str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
                 & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") &
                 & .and. NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN) call ham%Obser(GTTUP, PHASE, NT1, Mc_step_weight)
         end if

      end do

      if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
         & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") then   ! Finish calculating the forces
         do NT = Ltrot - THTROT + 1, Ltrot - 1
            ! UR is on time slice NT
            if (NT .eq. STAB_NT(NT_ST + 1)) then
               call Wrapur(STAB_NT(NT_ST), STAB_NT(NT_ST + 1), UDVR_local)
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  call CGRP(DetZ, GRUP(:, :, nf), udvr_local(nf_eff), udvst(nt_st + 1, nf_eff))
               end do
               NT_ST = NT_ST + 1
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  call Control_Precision_tau(GTTUP(:, :, nf), GRUP(:, :, nf), Ndim)
               end do
               GTTUP = GRUP
            end if
            NT1 = NT + 1
            !TODO reconstruction required?
            call Langevin_HMC%Wrap_Forces(GTTUP, NT1)
            if (NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN) then
               if (Symm) then
                  call Hop_mod_Symm(GTTUP_T, GTTUP, nt1)
                  !call reconstruction of non-calculated flavor blocks
                  if (reconstruction_needed) call ham%GR_reconstruction(GTTUP_T)
                  call ham%Obser(GTTUP_T, PHASE, NT1, Mc_step_weight)
               else
                  !call reconstruction of non-calculated flavor blocks
                  if (reconstruction_needed) call ham%GR_reconstruction(GTTUP)
                  call ham%Obser(GTTUP, PHASE, NT1, Mc_step_weight)
               end if
            end if
         end do
      end if

      do nf_eff = 1, N_FL_eff
         call udvr_local(nf_eff)%dealloc
      end do
      deallocate (GRUPB, GRUP, G00UP, G0TUP, GT0UP, GTTUP, TEMP, udvr_local)
      if (Symm) then
         deallocate (G00UP_T, G0TUP_T, GT0UP_T, GTTUP_T)
      end if

   end subroutine Tau_p

end module Tau_p_mod
