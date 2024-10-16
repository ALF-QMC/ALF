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
