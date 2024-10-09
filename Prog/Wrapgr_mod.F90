!  Copyright (C) 2016 - 2022 The ALF project
!
!  This file is part of the ALF project.
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
!> Handles global updates on a single time slice
!
!--------------------------------------------------------------------

module Wrapgr_mod

   use Hamiltonian_main
   use MyMats
   use Operator_mod
   use Control
   use Random_Wrap
   use Fields_mod
   use Hamiltonian_main
   use Hop_mod
   use upgrade_mod

   implicit none

   !> Privat
   complex(Kind=kind(0.d0)), private, allocatable ::  GR_ST(:, :, :)

contains
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Allocate, Deallocate space
!--------------------------------------------------------------------
   subroutine Wrapgr_alloc
      implicit none
      allocate (GR_ST(Ndim, Ndim, N_FL))
   end subroutine Wrapgr_alloc

   subroutine Wrapgr_dealloc
      implicit none
      deallocate (GR_ST)
   end subroutine Wrapgr_dealloc

!--------------------------------------------------------------------
   subroutine WRAPGRUP(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Given the green function matrix GR at time NTAU  the routine   propagates
!> it to time  NTAU + 1 and carries  out an update of the fields at time NTAU+1
!> NTAU: [0:LTROT-1]
!
!--------------------------------------------------------------------
      implicit none

      ! Arguments
      complex(Kind=kind(0.d0)), intent(INOUT), allocatable ::  GR(:, :, :)
      complex(Kind=kind(0.d0)), intent(INOUT) ::  PHASE
      integer, intent(IN) :: NTAU
      logical, intent(IN) :: Propose_S0
      integer, intent(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

      !Local
      integer :: nf, nf_eff, N_Type, NTAU1, n, m
      complex(Kind=kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
      real(Kind=kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
      character(Len=64)        :: Mode
      logical                   :: Acc, toggle1

      ! Wrap up, upgrade ntau1.  with B^{1}(tau1)
      NTAU1 = NTAU + 1
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call HOP_MOD_mmthr(GR(:, :, nf), nf, ntau1)
         call HOP_MOD_mmthl_m1(GR(:, :, nf), nf, ntau1)
      end do
      do n = Nt_sequential_start, Nt_sequential_end
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            HS_Field = nsigma%f(n, ntau1)
            N_type = 1
            call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau1)
         end do
         nf = 1
         T0_proposal = 1.5d0
         T0_Proposal_ratio = 1.d0
         Hs_New = nsigma%flip(n, ntau1)
         S0_ratio = ham%S0(n, ntau1, Hs_New)
         if (Propose_S0) then
            if (Op_V(n, nf)%type == 1) then
               T0_proposal = 1.d0 - 1.d0/(1.d0 + S0_ratio)
               T0_Proposal_ratio = 1.d0/S0_ratio
            end if
         end if
         if (T0_proposal > ranf_wrap()) then
            !Write(6,*) 'Hi', n, Op_V(n,nf)%type, T0_Proposal_ratio, S0_ratio
            mode = "Final"
            Prev_Ratiotot = cmplx(1.d0, 0.d0, kind(0.d0))
            call Upgrade2(GR, n, ntau1, PHASE, HS_new, Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode)
         else
            toggle1 = .false.
            call Control_upgrade_eff(toggle1)
         end if
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            N_type = 2
            call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau1)
         end do
      end do

      if (N_Global_tau > 0) then
         m = Nt_sequential_end
         !if ( Nt_sequential_start >  Nt_sequential_end ) m = Nt_sequential_start
         call Wrapgr_Random_update(GR, m, ntau1, PHASE, N_Global_tau)
         call Wrapgr_PlaceGR(GR, m, size(OP_V, 1), ntau1)
      end if

   end subroutine WRAPGRUP

!--------------------------------------------------------------------
   subroutine WRAPGRDO(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Given the green function matrix GR at time NTAU  the routine   carries out an
!> update of the fields at time NTAU and  propagates  it to time NTAU-1
!> NTAU: [LTROT:1]
!
!--------------------------------------------------------------------
      implicit none

      ! Given GREEN at time NTAU => GREEN at time NTAU - 1,
      ! Upgrade NTAU  [LTROT:1]

      complex(Kind=kind(0.d0)), intent(INOUT), allocatable :: GR(:, :, :)
      complex(Kind=kind(0.d0)), intent(INOUT) :: PHASE
      integer, intent(IN) :: NTAU
      logical, intent(IN) :: Propose_S0
      integer, intent(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

      ! Local
      integer :: nf, nf_eff, N_Type, n, m
      complex(Kind=kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
      real(Kind=kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
      character(Len=64)        :: Mode
      logical                   :: Acc, toggle1

      if (N_Global_tau > 0) then
         m = size(OP_V, 1)
         !Write(6,*) 'Call Ran_up ', m,ntau
         call Wrapgr_Random_update(GR, m, ntau, PHASE, N_Global_tau)
         call Wrapgr_PlaceGR(GR, m, Nt_sequential_end, ntau)
      end if

      do n = Nt_sequential_end, Nt_sequential_start, -1
         N_type = 2
         nf = 1
         HS_Field = nsigma%f(n, ntau)
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call Op_Wrapdo(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
         end do
         !Write(6,*) 'Upgrade : ', ntau,n
         nf = 1
         T0_proposal = 1.5d0
         T0_Proposal_ratio = 1.d0
         HS_new = nsigma%flip(n, ntau)
         S0_ratio = ham%S0(n, ntau, HS_new)
         if (Propose_S0) then
            if (Op_V(n, nf)%type == 1) then
               T0_proposal = 1.d0 - 1.d0/(1.d0 + S0_ratio)
               T0_Proposal_ratio = 1.d0/S0_ratio
            end if
         end if
         if (T0_proposal > ranf_wrap()) then
            mode = "Final"
            Prev_Ratiotot = cmplx(1.d0, 0.d0, kind(0.d0))
            call Upgrade2(GR, n, ntau, PHASE, HS_new, Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode)
         else
            toggle1 = .false.
            call Control_upgrade_eff(toggle1)
         end if

         !Call Upgrade(GR,n,ntau,PHASE,Op_V(n,1)%N_non_zero)
         ! The spin has changed after the upgrade!
         nf = 1
         HS_Field = nsigma%f(n, ntau)
         N_type = 1
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call Op_Wrapdo(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
         end do
      end do
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call Hop_mod_mmthl(GR(:, :, nf), nf, ntau)
         call Hop_mod_mmthr_m1(GR(:, :, nf), nf, ntau)
      end do

   end subroutine WRAPGRDO

!--------------------------------------------------------------------
   subroutine Wrapgr_PlaceGR(GR, m, m1, ntau)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> The Green function on a given time slice reads
!> G(tau) = (  1 + B(tau) B(tau-1,0)  B(beta,tau))^(-1) with B(tau) =   U_n e^(d_n) U_n^(dag) .... U_1 e^(V_1) U_1^(dag)  e^(-dtau H_t)
!> On input you have
!> G(tau,m)  = [ 1 + U_m e^(d_m) U_m^(dag) U_m^(dag) ... U_1 e^(V_1) U_1^(dag) e^(-dtau H_t) B(tau-1,0)
!>                   B(Beta,tau)  U_n e^(d_n) U_n^(dag) ...U_(m+1) e^(d_(m+1)) U_(m+1)^(dag) U_(m+1) ]
!> On output you have
!> G(tau,  m1 )
!>
!> Note that m,m1 in [0,n]
!--------------------------------------------------------------------

      implicit none

      !Arguments
      complex(Kind=kind(0.d0)), dimension(:, :, :), intent(INOUT), allocatable :: GR
      integer, intent(IN) :: m, m1, ntau

      !Local
      integer :: n, nf, nf_eff, N_Type
      complex(Kind=kind(0.d0)) :: HS_Field

      if (m == m1) then
         return
      elseif (m1 > m) then
         !Write(6,*) "Wrapup from ",  m + 1, "to",  m1, " on tau=",  ntau
         do n = m + 1, m1
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               HS_Field = nsigma%f(n, ntau)
               N_type = 1
               call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
            end do
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               N_type = 2
               call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
            end do
         end do
      elseif (m1 < m) then
         !Write(6,*) "Wrapdo from ",  m, "to",  m1 + 1
         do n = m, m1 + 1, -1
            N_type = 2
            nf = 1
            HS_Field = nsigma%f(n, ntau)
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call Op_Wrapdo(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
            end do
            !Write(6,*) 'Upgrade : ', ntau,n
            nf = 1
            HS_Field = nsigma%f(n, ntau)
            N_type = 1
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call Op_Wrapdo(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
            end do
         end do
      end if

   end subroutine Wrapgr_PlaceGR

!--------------------------------------------------------------------
   subroutine Wrapgr_Random_update(GR, m, ntau, PHASE, N_Global_tau)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> On input:
!> GR(tau,m) as defined in  Global_tau_mod_PlaceGR and the direction of updating scheme
!> direction=u --> You are visiting the time slices from tau = 1  to tau =Ltrot
!> direction=d --> You are visiting the time slices from tau = Ltrot   to tau = 1
!>
!> The routine calls
!> Global_move_tau(T0_Proposal_ratio, Flip_list, Flip_length,ntau,m,direction)
!> in the Hamiltonian module and then carries out the update
!>
!> On output
!>
!> Flip_length==1
!>        Green function is on  GR(tau,Flip_list(1) +1 )  if direction = u
!>        Green function is on  GR(tau,Flip_list(1) -1 )  if direction = d
!>        This is valid if the move has or has not been accepted.
!>
!> Flip_length > 1
!>        Let m_min = min(Flip_list), m_max = max(Flip_list)
!>        direction = u -->  On output Green on m_max is accepted. Green is on m_min if not accepted.
!>        direction = d -->  On output Green on m_min if accepted. Green is on m_max if not accepted.
!--------------------------------------------------------------------

      implicit none

      ! Arguments
      complex(Kind=kind(0.d0)), dimension(:, :, :), intent(INOUT), allocatable :: GR
      integer, intent(INOUT) :: m
      integer, intent(IN)    :: ntau, N_Global_tau
      complex(Kind=kind(0.d0)), intent(INOUT) :: PHASE

      ! Space for local variables
      integer                   :: n, Flip_length, nf, nf_eff, N_Type, ng_c, Flip_count
      real(Kind=kind(0.d0)) :: T0_Proposal_ratio, T0_proposal, S0_ratio
      complex(Kind=kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
      logical                   :: Acc
      character(Len=64)        :: Mode
      integer, allocatable :: Flip_list(:)
      complex(Kind=kind(0.d0)), allocatable :: Flip_value(:), Flip_value_st(:)
      real(Kind=kind(0.d0)) :: Zero = 10d-8

      allocate (Flip_list(size(Op_V, 1)), Flip_value(size(Op_V, 1)), Flip_value_st(size(Op_V, 1)))

      do ng_c = 1, N_Global_tau
         ! New configuration
         call ham%Global_move_tau(T0_Proposal_ratio, S0_ratio, Flip_list, Flip_length, Flip_value, ntau)
         !Write(6,*)  "Calling global move",  m, Flip_list(1), nsigma(Flip_list(1),ntau),Flip_value(1)
         if (T0_Proposal_ratio > Zero) then
            ! Order the list
            call wrapgr_sort(Flip_length, Flip_list, Flip_value)
            if (Flip_length > 1) then
               do Flip_count = 1, Flip_length - 1
                  Flip_value_st(Flip_count) = nsigma%f(Flip_list(Flip_count), ntau)
               end do
            end if
            Prev_Ratiotot = cmplx(1.d0, 0.d0, kind(0.d0))
            !Write(6,*) "-----", Flip_length
            do Flip_count = 1, Flip_length
               n = Flip_list(Flip_count)
               !Write(6,*)  "PlaceGR",  m, n-1,ntau
               call Wrapgr_PlaceGR(GR, m, n - 1, ntau)
               !Write(6,*)  "Back from PlaceGR",  m, n-1,ntau
               if (Flip_count == 1 .and. Flip_length > 1) GR_st = Gr
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  HS_Field = nsigma%f(n, ntau)
                  N_type = 1
                  call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
               end do
               nf = 1
               if (Flip_count < Flip_length) then
                  mode = "Intermediate"
                  HS_new = Flip_value(Flip_count)
                  call Upgrade2(GR, n, ntau, PHASE, HS_new, &
                       &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode)
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     N_type = 2
                     call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
                  end do
               else
                  !Write(6,*)  "Call Up mode final", n,ntau
                  mode = "Final"
                  HS_new = Flip_value(Flip_count)
                  call Upgrade2(GR, n, ntau, PHASE, HS_new, &
                       &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode)
                  !Write(6,*)  "Back from up mode final", n,ntau
                  !Write(6,*)  "Acceptance", Acc
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     N_type = 2
                     call Op_Wrapup(Gr(:, :, nf), Op_V(n, nf), HS_Field, Ndim, N_Type, ntau)
                  end do
               end if
               m = n
            end do
            if (.not. Acc .and. Flip_length > 1) then
               Gr = Gr_st
               m = Flip_list(1) - 1
               do Flip_count = 1, Flip_length - 1
                  nsigma%f(Flip_list(Flip_count), ntau) = Flip_value_st(Flip_count)
               end do
            end if
            !If (Acc) Call Hamiltonian_Print(Ntau)
         end if
      end do

      deallocate (Flip_list, Flip_value, Flip_value_st)

   end subroutine Wrapgr_Random_update

!----------------------------------------------------------------------------
   subroutine Wrapgr_sort(Flip_length, Flip_list, Flip_value)

      ! Arguments
      integer, intent(IN) :: Flip_length
      integer, intent(INOUT), allocatable :: Flip_list(:)
      complex(Kind=kind(0.d0)), intent(INOUT), allocatable :: Flip_value(:)

      ! Local
      integer :: swaps            ! number of swaps made in one pass
      integer :: nc               ! loop variable
      integer :: temp, n          ! temporary holder for making swap
      complex(Kind=kind(0.d0))      :: X

      if (Flip_length == 1) return

      !Write(6,*) 'Before sort'
      !DO nc = 1,Flip_length
      !   Write(6,*) Flip_list(nc),  Flip_value(nc)
      !Enddo

      do
         swaps = 0           ! Initially, we've made no swaps
         do nc = 1, (Flip_length - 1)
            if (Flip_list(nc) > Flip_list(nc + 1)) then
               temp = Flip_list(nc)
               Flip_list(nc) = Flip_list(nc + 1)
               Flip_list(nc + 1) = temp
               X = Flip_value(nc)
               Flip_value(nc) = Flip_value(nc + 1)
               Flip_value(nc + 1) = X
               swaps = swaps + 1
            end if
         end do
         if (swaps == 0) exit ! do count swaps
      end do

      !Write(6,*) 'After sort'
      !DO nc = 1,Flip_length
      !   Write(6,*) Flip_list(nc),  Flip_value(nc)
      !Enddo

   end subroutine Wrapgr_sort

!----------------------------------------------------------------------------

   subroutine Wrapgr_Test(Gr, ntau)

      ! Arguments
      complex(Kind=kind(0.d0)), dimension(:, :, :), intent(INOUT), allocatable :: GR
      integer :: ntau

      ! Local
      integer :: m, m1, N_size, nf, nf_eff, nth
      real(Kind=kind(0.d0)) :: Xmean, Xmax

      !Input is the Green function on time slice tau
      N_size = size(OP_V, 1)
      GR_ST = GR
      m = N_Size
      m1 = 0
      do nth = 1, 10
         call Wrapgr_PlaceGR(GR, m, m1, ntau)
         write (6, *) m, m1
         m = m1
         m1 = nranf(N_Size) - 1
      end do
      call Wrapgr_PlaceGR(GR, m, N_Size, ntau)
      write (6, *) m, N_size

      Xmax = 0.d0
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call COMPARE(GR_st(:, :, nf), GR(:, :, nf), XMAX, XMEAN)
      end do
      write (6, *) 'Compare Global_tau_mod_Test ', Xmax
      deallocate (GR_ST)

   end subroutine Wrapgr_Test

end module Wrapgr_mod
