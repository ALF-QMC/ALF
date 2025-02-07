module gfun_mod
   implicit none
contains

   subroutine cgrp(zdet, gr, kappa, kappa_bar, phir_up, phir_dn, phil)
      use udv_state_mod
      class(udv_state), intent(in) :: phil, phir_up, phir_dn
      complex(kind=kind(0.d0)), dimension(:, :, :), intent(out) :: gr, kappa, kappa_bar
      complex(kind=kind(0.d0)), intent(out)  :: zdet

      complex(Kind=kind(0.d0)), allocatable, dimension(:, :) :: qmat, zphi1, zphi2
      complex(Kind=kind(0.d0)), allocatable, dimension(:, :) :: tmp_mat1, tmp_mat0
      real(kind=kind(0.d0)) :: pow
      integer, allocatable :: ipiv(:)
      complex(Kind=kind(0.d0)) :: alpha, beta, phase, pfa_phase
      integer :: Ndim, N_part, info, n

      ndim   = phir_up%ndim
      n_part = phir_up%n_part
      allocate (qmat(n_part, n_part), tmp_mat1(n_part, ndim))
      allocate (tmp_mat0(ndim,ndim))
      !! zpih1 = \phi_{\dn}^{T}*Z^{\dagger}
      !! zphi2 = \phi_{\up}^{T}*Z^{*}
      allocate (zphi1(n_part, ndim), zphi2(n_part, ndim))
      allocate (ipiv(n_part))

      ! Gr = Ur (Ul Ur)^-1 Ul
      ! Phase = 1 + Ur (Ul Ur)^-1 Ul
      ! Ul = udvl%U ^dag
      alpha = 1.d0
      beta = 0.d0
      
      !! zpih1 = \phi_{\dn}^{T}*Z^{\dagger}
      call zgemm('t', 'c', n_part, ndim, ndim, alpha, phir_dn%u(1,1), ndim, & 
          &    phil%u(1,1), ndim, beta, zphi1(1, 1), n_part)
      
      !! zphi2 = \phi_{\up}^{T}*Z^{*}
      tmp_mat0 = conjg(phil%u(:,:))
      call zgemm('t', 'n', n_part, ndim, ndim, alpha, phir_up%u(1,1), ndim, & 
          &    tmp_mat0(1,1), ndim, beta, zphi2(1, 1), n_part)
      
      !! Qmat = \phi^{T}_{\dn}Z^{\dagger}*\phi_{\up}
      call zgemm('n', 'n', n_part, n_part, ndim, alpha, zphi1(1,1), n_part, & 
          &    phir_up%u(1,1), ndim, beta, qmat(1, 1), n_part)
      
      !! inverse of Qmat
      ! ZGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call zgetrf(n_part, n_part, qmat, n_part, ipiv, info)
      ! obtain log of det
      zdet = 0.d0
      phase = 1.d0
      do n = 1, N_part
         if (ipiv(n) .ne. n) then
            phase = -phase
         end if
         zdet = zdet + log(qmat(n, n))
      end do
      zdet = zdet + log(phase)
      pfa_phase = cmplx(-1.d0,0.d0,kind(0.d0))**(dble(n_part)*dble(n_part-1)*0.5d0) 
      zdet = zdet + log(pfa_phase)

      !!! gr_{\up\up}
      tmp_mat1 = transpose(phir_up%u)
      alpha = 1.d0; beta=0.d0
      call zgetrs('t', n_part, ndim, qmat(1, 1), n_part, ipiv, tmp_mat1(1, 1), n_part, info)
      call zgemm('t', 'n', ndim, ndim, n_part, alpha, zphi1(1, 1), n_part, tmp_mat1(1, 1), n_part, beta, gr(:,:,1), ndim)
      
      !!! kappa_bar_{\up\dn}
      alpha = 1.d0
      tmp_mat1 = zphi2
      call zgetrs('t', n_part, ndim, qmat(1, 1), n_part, ipiv, tmp_mat1(1, 1), n_part, info)
      kappa_bar(:,:,1) = tmp_mat0
      alpha = -1.d0; beta = 1.d0
      call zgemm('t', 'n', ndim, ndim, n_part, alpha, zphi1(1,1), n_part, tmp_mat1(1,1), n_part, beta, kappa_bar(:,:,1), ndim)
      !!! kappa_bar_{\dn\up}
      kappa_bar(:,:,2) = -transpose(kappa_bar(:,:,1))
      
      !!! gr_{\dn\dn}
      tmp_mat1 = transpose(phir_dn%u)
      alpha = 1.d0; beta = 0.d0
      call zgetrs('n', n_part, ndim, qmat(1, 1), n_part, ipiv, tmp_mat1(1, 1), n_part, info)
      call zgemm('t', 'n', ndim, ndim, n_part, alpha, zphi2(1, 1), n_part, tmp_mat1(1, 1), n_part, beta, gr(:,:,2), ndim)
      
      !!! kappa_{\up\dn}
      alpha = -1.d0; beta=0.d0
      call zgemm('n', 'n', ndim, ndim, n_part, alpha, phir_up%u(1, 1), ndim, tmp_mat1(1, 1), n_part, beta, kappa(:,:,1), ndim)
      !!! kappa_{\dn\up}
      kappa(:,:,2) = -transpose(kappa(:,:,1))

      deallocate (qmat, zphi1, tmp_mat0, tmp_mat1)
      deallocate (ipiv)

   end subroutine cgrp

end module gfun_mod
