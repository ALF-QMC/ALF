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

module wrapul_mod
   implicit none
contains

   subroutine WRAPUL(NTAU1, NTAU, UDVL)

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Imaginary time propagation and udv decompostion from Ntau1 to Ntau, Ntau1 > Ntau
!>
!> @details
!> On input   ( B(Beta,NTAU1) )^{+}  =   UL * DL * VL \n
!> On output  ( B(Beta,NTAU ) )^{+}  =   UL * DL * VL
!>
!> @param[in] Ntau, Ntau1  Integer
!> @param[inout] UDVR Class(UDV_state)
!> \verbatim
!>  The UDV_state class carries the information of left or right propoagation.
!> \endverbatim
!--------------------------------------------------------------------

      use UDV_State_mod
#if defined(STAB2) ||  defined(STAB1)
      use Hamiltonian_main
      use Hop_mod
      use UDV_Wrap_mod

      implicit none

      ! Arguments
      class(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
      integer :: NTAU1, NTAU

      ! Working space.
      complex(Kind=kind(0.d0)) ::  U1(Ndim, Ndim), V1(Ndim, Ndim), TMP(Ndim, Ndim), TMP1(Ndim, Ndim)
      complex(Kind=kind(0.d0)) ::  Z_ONE, beta
      integer :: NT, NCON, n, nf, nf_eff, sign = 1
      real(Kind=kind(0.d0)) ::  X

      NCON = 0  ! Test for UDV ::::  0: Off,  1: On.

      Z_ONE = cmplx(1.d0, 0.d0, kind(0.d0))
      beta = 0.d0
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call INITD(TMP, Z_ONE)
         do NT = NTAU1, NTAU + 1, -1
            do n = size(Op_V, 1), 1, -1
               call Op_mmultL(Tmp, Op_V(n, nf), nsigma%f(n, nt), 'n', nt, sign)
            end do
            !CALL MMULT( TMP1,Tmp,Exp_T(:,:,nf) )
            call Hop_mod_mmthl(Tmp, nf, nt)
            ! Tmp = Tmp1
         end do

         !Carry out U,D,V decomposition.
         call ZGEMM('C', 'N', Ndim, UDVL(nf_eff)%N_part, Ndim, Z_ONE, TMP, Ndim, udvl(nf_eff)%U(1, 1), Ndim, beta, TMP1, Ndim)
         if (allocated(UDVL(nf_eff)%V)) then
            do n = 1, UDVL(nf_eff)%N_part
               TMP1(:, n) = TMP1(:, n)*udvl(nf_eff)%D(n)
            end do
            call UDV_WRAP_Pivot(TMP1, udvl(nf_eff)%U, udvl(nf_eff)%D, V1, NCON, Ndim, Ndim)
            call ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_ONE, udvl(nf_eff)%V(1, 1), Ndim, V1, Ndim, beta, TMP1, Ndim)
            udvl(nf_eff)%V = TMP1
         else
            call UDV_WRAP_Pivot(TMP1(:, 1:UDVL(nf_eff)%N_part), udvl(nf_eff)%U, udvl(nf_eff)%D, V1, NCON, Ndim, UDVL(nf_eff)%N_part)
         end if
         ! Test if scales in D are appraoching the limit of double precision.
         call UDVL(nf_eff)%testscale
      end do

#else
      use Hop_mod
      implicit none

      ! Arguments
      class(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
      integer, intent(in) :: NTAU1, NTAU

      integer :: NT, n, nf, nf_eff

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do NT = NTAU1, NTAU + 1, -1
            do n = size(Op_V, 1), 1, -1
               call Op_mmultR(udvl(nf_eff)%U, Op_V(n, nf), nsigma%f(n, nt), 'c', nt)
            end do
            call Hop_mod_mmthlc(udvl(nf_eff)%U, nf, nt)
         end do

         !Carry out U,D,V decomposition.
         call UDVL(nf_eff)%decompose
         ! Test if scales in D are appraoching the limit of double precision.
         call UDVL(nf_eff)%testscale
      end do
#endif
   end subroutine WRAPUL

end module wrapul_mod
