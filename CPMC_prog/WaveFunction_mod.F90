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

   subroutine wf_overlap(wf_l, wf_r, z_norm)
      implicit none
      type(WaveFunction), dimension(:), allocatable, intent(in)    :: wf_l
      type(WaveFunction), dimension(:), allocatable, intent(inout) :: wf_r
      complex(Kind=kind(0.d0)), intent(out) :: z_norm

      ! Local
      integer :: N_Part, Ndim, n, ne, n_fl, nf
      complex(Kind=kind(0.d0)), allocatable :: tmp_mat1(:,:), tmp_mat2(:,:)
      complex(Kind=kind(0.d0)) :: alpha, beta

      N_fl   = size(wf_r, 1)
      N_part = size(wf_r(1)%p, 2)
      Ndim   = size(wf_r(1)%p, 1)
      allocate (tmp_mat1(ndim, n_part), tmp_mat2(n_part,n_part))

      alpha = 1.d0
      beta = 0.d0
      !! Z^{\dagger}*\phi_{\uparrow}
      call zgemm('c', 'n', ndim, n_part, ndim, alpha, wf_l(1)%p(1,1), ndim, wf_r(1)%p(1,1), ndim, beta, tmp_mat1(1, 1), ndim)
      !! \phi^{T}_{\downarrow}Z^{\dagger}*\phi_{\uparrow}
      call zgemm('t', 'n', n_part, n_part, ndim, alpha, wf_r(2)%p(1,1), ndim, tmp_mat1(1,1), ndim, beta, tmp_mat2(1, 1), n_part)

      z_norm = det(tmp_mat2, N_part)

      z_norm = cmplx(-1.d0,0.d0,kind(0.d0))**(dble(n_part)*dble(n_part-1)*0.5d0) * z_norm

      z_norm = (cmplx(1.d0, 0.d0, kind(0.d0))/z_norm)**(1.d0/dble(2*n_part))

      !wf_l(1)%p = z_norm*wf_l(1)%p
      wf_r(1)%p = z_norm*wf_r(1)%p(:,:)
      wf_r(2)%p = z_norm*wf_r(2)%p(:,:)

      deallocate (tmp_mat1, tmp_mat2)

   end subroutine wf_overlap

!--------------------------------------------------------------------

   pure subroutine wf_alloc(WF, Ndim, N_part)
      implicit none
      type(WaveFunction), intent(INOUT) :: WF
      integer, intent(IN) :: Ndim, N_part
      allocate (WF%P(Ndim, N_part))
      WF%P = cmplx(0.d0, 0.d0, kind(0.d0))
   end subroutine wf_alloc

!--------------------------------------------------------------------

   pure subroutine wf_clear(WF)
      implicit none
      type(WaveFunction), intent(INOUT) :: WF
      deallocate (WF%P)
   end subroutine wf_clear

end module WaveFunction_mod
