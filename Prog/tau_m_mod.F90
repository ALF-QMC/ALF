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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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
!> This module  handles calculation of imagimary time displaced Green functions and
!> calls the routine ObserT.F90 in the Hamiltonian module,  so as to compute the
!> defined  time dispalced correlations functions. This modules is for the finite temperature code.
!>
!
!--------------------------------------------------------------------
module Tau_m_mod
   use Hamiltonian_main
   use Operator_mod
   use Control
   use Hop_mod
   use UDV_State_mod
   use Langevin_HMC_mod
   use wrapur_mod
   use cgr2_2_mod

contains

   subroutine TAU_M(udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN)

      implicit none

      integer, intent(In) :: NSTM, NWRAP
      class(UDV_State), dimension(:, :), allocatable, intent(IN) :: udvst
      complex(Kind=kind(0.d0)), intent(in) :: GR(NDIM, NDIM, N_FL), Phase
      integer, intent(In) :: STAB_NT(0:NSTM)
      integer, intent(In) :: LOBS_ST, LOBS_EN

      ! Local
      ! This could be placed as  private for the module
      complex(Kind=kind(0.d0))  :: GT0(NDIM, NDIM, N_FL), G00(NDIM, NDIM, N_FL), GTT(NDIM, NDIM, N_FL), G0T(NDIM, NDIM, N_FL)
      complex(Kind=kind(0.d0)), dimension(:, :, :), allocatable  :: GT0_T, G00_T, GTT_T, G0T_T
      class(UDV_State), dimension(:), allocatable :: udvr
      complex(Kind=kind(0.d0))  :: HLP4(Ndim, Ndim), HLP5(Ndim, Ndim), HLP6(Ndim, Ndim)

      complex(Kind=kind(0.d0))  ::  Z
      integer  ::  I, J, nf, nf_eff, NT, NT1, NTST, NST, N, N_type
      real(Kind=kind(0.d0))  ::  spin, Mc_step_Weight

      if (Symm) then
         allocate (G00_T(Ndim, Ndim, N_FL), G0T_T(Ndim, Ndim, N_FL), GT0_T(Ndim, Ndim, N_FL), GTT_T(Ndim, Ndim, N_FL))
      end if

      Mc_step_Weight = 1.d0
      !needed for integration weight in Langevin, dt is taken care of by acceptance ratio for HMC
      if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN") Mc_step_weight = Langevin_HMC%get_Delta_t_running()

      !Tau = 0
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do J = 1, Ndim
            do I = 1, Ndim
               Z = cmplx(0.d0, 0.d0, kind(0.d0))
               if (I == J) Z = cmplx(1.d0, 0.d0, kind(0.d0))
               G00(I, J, nf) = GR(I, J, nf)
               GT0(I, J, nf) = GR(I, J, nf)
               GTT(I, J, nf) = GR(I, J, nf)
               G0T(I, J, nf) = -(Z - GR(I, J, nf))
            end do
         end do
      end do
      NT = 0
      ! In Module Hamiltonian
      if (Symm) then
         call Hop_mod_Symm(G00_T, G00, nt)
         call Hop_mod_Symm(GTT_T, GTT, nt)
         call Hop_mod_Symm(G0T_T, G0T, nt)
         call Hop_mod_Symm(GT0_T, GT0, nt)
         !call reconstruction of non-calculated flavor blocks
         if (reconstruction_needed) then
            call ham%GR_reconstruction(G00_T)
            call ham%GR_reconstruction(GTT_T)
            call ham%GRT_reconstruction(GT0_T, G0T_T)
         end if
         call ham%OBSERT(NT, GT0_T, G0T_T, G00_T, GTT_T, PHASE, Mc_step_Weight)
      else
         !call reconstruction of non-calculated flavor blocks
         if (reconstruction_needed) then
            call ham%GR_reconstruction(G00)
            call ham%GR_reconstruction(GTT)
            call ham%GRT_reconstruction(GT0, G0T)
         end if
         call ham%OBSERT(NT, GT0, G0T, G00, GTT, PHASE, Mc_step_Weight)
      end if

      allocate (udvr(N_FL_eff))
      Z = cmplx(1.d0, 0.d0, kind(0.d0))
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         if (Projector) then
            call udvr(nf_eff)%init(ndim, 'r', WF_R(nf)%P)
         else
            call udvr(nf_eff)%init(ndim, 'r')
         end if
      end do

      NST = 1
      do NT = 0, LTROT - 1
         ! Now wrapup:
         NT1 = NT + 1
         call PROPR(GT0, NT1)
         call PROPRM1(G0T, NT1)
         if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
         & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") then
            call Langevin_HMC%Wrap_Forces(GTT, NT1)
         else
            call PROPRM1(GTT, NT1)
            call PROPR(GTT, NT1)
         end if
         ! In Module Hamiltonian
         if (Symm) then
            call Hop_mod_Symm(G00_T, G00, 0)
            call Hop_mod_Symm(GTT_T, GTT, nt1)
            call Hop_mod_Symm(G0T_T, G0T, Ltrot, nt1)
            call Hop_mod_Symm(GT0_T, GT0, nt1, Ltrot)
            !call reconstruction of non-calculated flavor blocks
            if (reconstruction_needed) then
               call ham%GR_reconstruction(G00_T)
               call ham%GR_reconstruction(GTT_T)
               call ham%GRT_reconstruction(GT0_T, G0T_T)
            end if
            call ham%OBSERT(NT1, GT0_T, G0T_T, G00_T, GTT_T, PHASE, Mc_step_weight)
            if ((str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN"  &
                 & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") &
                 &  .and. NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN) call ham%Obser(GTT_T, PHASE, NT1, Mc_step_weight)
         else
            !call reconstruction of non-calculated flavor blocks
            if (reconstruction_needed) then
               call ham%GR_reconstruction(G00)
               call ham%GR_reconstruction(GTT)
               call ham%GRT_reconstruction(GT0, G0T)
            end if
            call ham%OBSERT(NT1, GT0, G0T, G00, GTT, PHASE, Mc_step_weight)
            if ((str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
                 & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") &
                 & .and. NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN) call ham%Obser(GTT, PHASE, NT1, Mc_step_weight)
         end if

         if (Stab_nt(NST) == NT1 .and. NT1 .ne. LTROT) then
            !NTST = NT1 - NWRAP
            !NST  = NT1/(NWRAP)
            NTST = Stab_nt(NST - 1)
            ! WRITE(6,*) 'NT1, NST: ', NT1,NST
            call WRAPUR(NTST, NT1, udvr)
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               HLP4(:, :) = GTT(:, :, nf)
               HLP5(:, :) = GT0(:, :, nf)
               HLP6(:, :) = G0T(:, :, nf)
               call CGR2_2(GT0(:, :, nf), G00(:, :, nf), GTT(:, :, nf), G0T(:, :, nf), &
                    & udvr(nf_eff), udvst(NST, nf_eff), NDIM)
               call Control_Precision_tau(GR(:, :, nf), G00(:, :, nf), Ndim)
               call Control_Precision_tau(HLP4, GTT(:, :, nf), Ndim)
               call Control_Precision_tau(HLP5, GT0(:, :, nf), Ndim)
               call Control_Precision_tau(HLP6, G0T(:, :, nf), Ndim)
            end do
            NST = NST + 1
         end if
      end do

      do nf_eff = 1, N_Fl_eff
         call udvr(nf_eff)%dealloc
      end do
      deallocate (udvr)
      if (Symm) then
         deallocate (G00_T, G0T_T, GT0_T, GTT_T)
      end if

   end subroutine TAU_M

!--------------------------------------------------------------------

   subroutine PROPR(AIN, NT)

      ! Ain =       B(NT-1, NT1)
      ! Aout= Ain = B(NT  , NT1)

      implicit none
      complex(Kind=kind(0.d0)), intent(INOUT) :: Ain(Ndim, Ndim, N_FL)
      integer, intent(IN) :: NT

      !Locals
      integer :: nf, nf_eff, n

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call Hop_mod_mmthr(Ain(:, :, nf), nf, nt)
         do n = 1, size(Op_V, 1)
!                  X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
            call Op_mmultR(Ain(:, :, nf), Op_V(n, nf), nsigma%f(n, nt), 'n', nt)
         end do
      end do

   end subroutine PROPR

!--------------------------------------------------------------------

   subroutine PROPRM1(AIN, NT)

      !Ain = B^{-1}(NT-1, NT1)
      !Aout= B^{-1}(NT  , NT1)

      implicit none

      !Arguments
      complex(Kind=kind(0.d0)), intent(Inout) ::  AIN(Ndim, Ndim, N_FL)
      integer :: NT, sign = -1

      ! Locals
      integer :: nf, nf_eff, n

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         !Call MMULT(HLP4,Ain(:,:,nf),Exp_T_M1(:,:,nf) )
         call Hop_mod_mmthl_m1(Ain(:, :, nf), nf, nt)
         do n = 1, size(Op_V, 1)
            call Op_mmultL(Ain(:, :, nf), Op_V(n, nf), nsigma%f(n, nt), 'n', nt, sign)
         end do
      end do

   end subroutine PROPRM1

end module Tau_m_mod
