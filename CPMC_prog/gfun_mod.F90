module gfun_mod
   implicit none
contains

   subroutine cgrp(zdet, gr, kappa, kappa_bar, phir_up, phir_dn, phil)
      use udv_state_mod
      class(udv_state), intent(in) :: phil, phir_up, phir_dn
      complex(kind=kind(0.d0)), dimension(:, :, :), intent(out) :: gr, kappa, kappa_bar
      complex(kind=kind(0.d0)), intent(out)  :: zdet

      complex(Kind=kind(0.d0)), allocatable, dimension(:, :) :: qmat, zphi0, zphi1, zphi2
      complex(Kind=kind(0.d0)), allocatable, dimension(:, :) :: tmp_mat1, tmp_mat2
      real(kind=kind(0.d0)) :: pow
      integer, allocatable :: ipiv(:)
      complex(Kind=kind(0.d0)) :: alpha, beta, phase, pfa_phase
      integer :: Ndim, N_part, info, n

      ndim   = phir_up%ndim
      n_part = phir_up%n_part
      allocate (qmat(n_part, n_part), tmp_mat1(ndim, n_part), tmp_mat2(n_part, ndim))
      !! zphi1 = Z^{\dagger}\phi_{\up}
      !! zphi2 = Z^{*}\phi_{\dn}
      allocate (zphi0(ndim, ndim), zphi1(ndim, n_part), zphi2(ndim, n_part))

      ! Gr = Ur (Ul Ur)^-1 Ul
      ! Phase = 1 + Ur (Ul Ur)^-1 Ul
      ! Ul = udvl%U ^dag
      alpha = 1.d0
      beta = 0.d0
      
      !! zpih1 = Z^{\dagger}*\phi_{\up}
      call zgemm('c', 'n', ndim, n_part, ndim, alpha, phil%u(1,1), ndim, & 
          &    phir_up%u(1,1), ndim, beta, zphi1(1, 1), n_part)
      
      !! Qmat = \phi^{T}_{\dn}Z^{\dagger}*\phi_{\up}
      call zgemm('c', 'n', n_part, ndim, n_part, alpha, phir_dn%u(1,1), ndim, & 
          &    tmp_mat1(1,1), n_part, beta, qmat(1, 1), n_part)
      
      !! zpih2 = Z^{*}*\phi_{\dn}
      zphi0 = conjg(phil%u(1,1))
      call zgemm('n', 'n', ndim, n_part, ndim, alpha, zphi0, ndim, & 
          &    phir_dn%u(1,1), ndim, beta, zphi2(1, 1), n_part)

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
      tmp_mat2 = transpose(phir_up%u)
      alpha = 1.d0
      call zgetrs('t', n_part, ndim, qmat(1, 1), n_part, ipiv, tmp_mat2(1, 1), n_part, info)
      call zgemm('n', 'n', ndim, ndim, n_part, alpha, zphi2(1, 1), ndim, tmp_mat2(1, 1), n_part, beta, gr(:,:,1), ndim)
      
      !!! gr_{\dn\dn}
      tmp_mat2 = transpose(phir_dn%u)
      alpha = 1.d0
      call zgetrs('n', n_part, ndim, qmat(1, 1), n_part, ipiv, tmp_mat2(1, 1), n_part, info)
      call zgemm('n', 'n', ndim, ndim, n_part, alpha, zphi1(1, 1), ndim, tmp_mat2(1, 1), n_part, beta, gr(:,:,2), ndim)
      
      !!! kappa_{\up\dn}
      tmp_mat2 = transpose(phir_dn%u)
      alpha = -1.d0
      call zgetrs('n', n_part, ndim, qmat(1, 1), n_part, ipiv, tmp_mat2(1, 1), n_part, info)
      call zgemm('n', 'n', ndim, ndim, n_part, alpha, phir_up%u(1, 1), ndim, tmp_mat2(1, 1), n_part, beta, kappa(:,:,1), ndim)

      !!! kappa_{\dn\up}
      kappa(:,:,2) = transpose(kappa(:,:,1))
      
      !!! kappa_bar_{\up\dn}
      tmp_mat2 = transpose(zphi1)
      alpha = -1.d0
      call zgetrs('t', n_part, ndim, qmat(1, 1), n_part, ipiv, tmp_mat2(1, 1), n_part, info)
      call zgemm('n', 'n', ndim, ndim, n_part, alpha, zphi2(1, 1), ndim, tmp_mat2(1, 1), n_part, beta, kappa_bar(:,:,1), ndim)
      kappa_bar(:,:,1) = zphi0(:,:) + kappa_bar(:,:,1)

      !!! kappa_bar_{\dn\up}
      kappa_bar(:,:,2) = transpose(kappa_bar(:,:,1))


      deallocate (qmat, zphi0, zphi1, zphi2, tmp_mat1, tmp_mat2)
      deallocate (ipiv)

   end subroutine cgrp

end module gfun_mod
