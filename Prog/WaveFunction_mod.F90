!  Copyright (C) 2016, 2017 The ALF project
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
!> Defines the wavefunction type
!>
!
!--------------------------------------------------------------------

module WaveFunction_mod

   use MyMats

   implicit none

   type WaveFunction
      !> P is an Ndim x N_part matrix. N_part is the number of particles
      complex(Kind=kind(0.d0)), allocatable :: P(:, :)
      !> Degeneracy of trial wave function
      real(Kind=kind(0.d0)) :: Degen
   end type WaveFunction

contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> If the left and right trial wave functions are different, then one of the trial wave functions has to
!> be multiple by a phase such that  <Psi_L | Psi_R > = 1  otherwise, the program will produce a "fake"
!> negative sign problem.
!>
!
!--------------------------------------------------------------------

   subroutine WF_overlap(WF_L, WF_R, Z_norm)
      implicit none
      type(WaveFunction), intent(IN)     :: WF_L
      type(WaveFunction), intent(INOUT)     :: WF_R
      complex(Kind=kind(0.d0)), intent(OUT) :: Z_norm

      ! Local
      integer :: N_Part, Ndim, n, ne
      complex(Kind=kind(0.d0)), allocatable ::  mat(:, :)
      complex(Kind=kind(0.d0)) :: alpha, beta

      N_part = size(WF_R%P, 2)
      Ndim = size(WF_R%P, 1)
      allocate (Mat(N_part, N_Part))

      alpha = 1.d0
      beta = 0.d0
      call ZGEMM('C', 'N', N_part, N_part, Ndim, alpha, WF_L%P(1, 1), Ndim, WF_R%P(1, 1), Ndim, beta, Mat(1, 1), N_part)
      ! Mat = (WL_L%P)^{dagger} WL_L%R

      Z_norm = Det(Mat, N_part)

      Z_norm = (cmplx(1.d0, 0.d0, kind(0.d0))/Z_norm)**(1.d0/real(N_part, kind(0.d0)))

      WF_R%P = Z_norm*WF_R%P

      deallocate (Mat)

   end subroutine WF_overlap

!--------------------------------------------------------------------

   pure subroutine WF_alloc(WF, Ndim, N_part)
      implicit none
      type(WaveFunction), intent(INOUT) :: WF
      integer, intent(IN) :: Ndim, N_part
      allocate (WF%P(Ndim, N_part))
      WF%P = cmplx(0.d0, 0.d0, kind(0.d0))
   end subroutine WF_alloc

!--------------------------------------------------------------------

   pure subroutine WF_clear(WF)
      implicit none
      type(WaveFunction), intent(INOUT) :: WF
      deallocate (WF%P)
   end subroutine WF_clear

end module WaveFunction_mod
