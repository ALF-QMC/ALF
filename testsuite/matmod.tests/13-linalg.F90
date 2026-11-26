! Unit test for linalg_mod: every la_<name> must give bitwise the same result
! as the direct call of <name>.
program test_linalg
  use linalg_mod
  implicit none
  integer, parameter :: dp = kind(0.d0), n = 7, m = 5
  complex(dp), parameter :: alpha = (0.7d0, -0.3d0), beta = (0.2d0, 0.5d0)
  complex(dp) :: A(n,n), B(n,n), C1(n,n), C2(n,n), H(n,n), x(n), y1(n), y2(n)
  complex(dp) :: tau1(n), tau2(n), work1(4*n*n), work2(4*n*n)
  real(dp)    :: rwork1(5*n), rwork2(5*n), s1(n), s2(n)
  complex(dp) :: U1(n,n), U2(n,n), V1(n,n), V2(n,n), w1(n), w2(n)
  integer     :: ipiv1(n), ipiv2(n), jpvt1(n), jpvt2(n), info1, info2, n_fail
  external :: zgemm, zhemm, ztrsm, ztrmm, zgemv, zgeru, zscal, zcopy, zlacpy, &
       &      zgetrf, zgetrs, zgetri, zgeqrf, zgeqp3, zungqr, zunmqr, &
       &      zlaset, zlapmr, zlapmt, zgesvd, zheev, zgeev

  n_fail = 0
  call fill(A); call fill(B); call fill(H); H = H + conjg(transpose(H)); call fillv(x)

  C1 = B; C2 = B
  call zgemm('N', 'C', n, m, n, alpha, A, n, H, n, beta, C1, n)
  call la_zgemm('N', 'C', n, m, n, alpha, A, n, H, n, beta, C2, n)
  call check('zgemm', all(C1 == C2))

  C1 = B; C2 = B
  call zhemm('L', 'U', n, n, alpha, H, n, A, n, beta, C1, n)
  call la_zhemm('L', 'U', n, n, alpha, H, n, A, n, beta, C2, n)
  call check('zhemm', all(C1 == C2))

  C1 = B; C2 = B
  call ztrsm('L', 'U', 'N', 'N', n, n, alpha, A + 5.d0, n, C1, n)
  call la_ztrsm('L', 'U', 'N', 'N', n, n, alpha, A + 5.d0, n, C2, n)
  call check('ztrsm', all(C1 == C2))

  C1 = B; C2 = B
  call ztrmm('R', 'U', 'C', 'N', n, n, alpha, A, n, C1, n)
  call la_ztrmm('R', 'U', 'C', 'N', n, n, alpha, A, n, C2, n)
  call check('ztrmm', all(C1 == C2))

  y1 = x(n:1:-1); y2 = y1
  call zgemv('T', n, m, alpha, A, n, x, 1, beta, y1, 1)
  call la_zgemv('T', n, m, alpha, A, n, x, 1, beta, y2, 1)
  call check('zgemv', all(y1 == y2))

  C1 = B; C2 = B
  call zgeru(n, n, alpha, x, 1, y1, 1, C1, n)
  call la_zgeru(n, n, alpha, x, 1, y1, 1, C2, n)
  call check('zgeru', all(C1 == C2))

  C1 = B; C2 = B          ! a row, stride n
  call zscal(n, alpha, C1(2,1), n)
  call la_zscal(n, alpha, C2(2,1), n)
  call check('zscal', all(C1 == C2))

  y1 = 0; y2 = 0
  call zcopy(n, A(1,3), 1, y1, 1)
  call la_zcopy(n, A(1,3), 1, y2, 1)
  call check('zcopy', all(y1 == y2))

  C1 = 0; C2 = 0
  call zlacpy('U', n, n, A, n, C1, n)
  call la_zlacpy('U', n, n, A, n, C2, n)
  call check('zlacpy', all(C1 == C2))

  C1 = A; C2 = A
  call zgetrf(n, n, C1, n, ipiv1, info1)
  call la_zgetrf(n, n, C2, n, ipiv2, info2)
  call check('zgetrf', all(C1 == C2) .and. all(ipiv1 == ipiv2) .and. info1 == info2)

  y1 = B(:,1); y2 = y1    ! C1 = C2 = LU factors from here on
  call zgetrs('C', n, 1, C1, n, ipiv1, y1, n, info1)
  call la_zgetrs('C', n, 1, C2, n, ipiv2, y2, n, info2)
  call check('zgetrs', all(y1 == y2) .and. info1 == info2)

  call zgetri(n, C1, n, ipiv1, work1, size(work1), info1)
  call la_zgetri(n, C2, n, ipiv2, work2, size(work2), info2)
  call check('zgetri', all(C1 == C2) .and. info1 == info2)

  C1 = A; C2 = A
  call zgeqrf(n, n, C1, n, tau1, work1, size(work1), info1)
  call la_zgeqrf(n, n, C2, n, tau2, work2, size(work2), info2)
  call check('zgeqrf', all(C1 == C2) .and. all(tau1 == tau2) .and. info1 == info2)

  H = C1                  ! QR factors (and tau1) of zgeqrf, used by zunmqr and zungqr
  call fill(B)
  C1 = B; C2 = B
  call zunmqr('R', 'N', n, n, n, H, n, tau1, C1, n, work1, size(work1), info1)
  call la_zunmqr('R', 'N', n, n, n, H, n, tau1, C2, n, work2, size(work2), info2)
  call check('zunmqr', all(C1 == C2) .and. info1 == info2)

  C1 = H; C2 = H
  call zungqr(n, n, n, C1, n, tau1, work1, size(work1), info1)
  call la_zungqr(n, n, n, C2, n, tau1, work2, size(work2), info2)
  call check('zungqr', all(C1 == C2) .and. info1 == info2)

  C1 = A; C2 = A; jpvt1 = 0; jpvt2 = 0
  call zgeqp3(n, n, C1, n, jpvt1, tau1, work1, size(work1), rwork1, info1)
  call la_zgeqp3(n, n, C2, n, jpvt2, tau2, work2, size(work2), rwork2, info2)
  call check('zgeqp3', all(C1 == C2) .and. all(tau1 == tau2) .and. all(jpvt1 == jpvt2) .and. info1 == info2)

  C1 = B; C2 = B
  call zlaset('L', n, m, alpha, beta, C1, n)
  call la_zlaset('L', n, m, alpha, beta, C2, n)
  call check('zlaset', all(C1 == C2))

  ipiv1 = [(n + 1 - info1, info1 = 1, n)]; ipiv2 = ipiv1      ! a permutation
  C1 = A; C2 = A
  call zlapmr(.true., n, m, C1, n, ipiv1)
  call la_zlapmr(.true., n, m, C2, n, ipiv2)
  call check('zlapmr', all(C1 == C2) .and. all(ipiv1 == ipiv2))

  C1 = A; C2 = A
  call zlapmt(.false., m, n, C1, n, ipiv1)
  call la_zlapmt(.false., m, n, C2, n, ipiv2)
  call check('zlapmt', all(C1 == C2) .and. all(ipiv1 == ipiv2))

  C1 = A; C2 = A
  call zgesvd('A', 'A', n, n, C1, n, s1, U1, n, V1, n, work1, size(work1), rwork1, info1)
  call la_zgesvd('A', 'A', n, n, C2, n, s2, U2, n, V2, n, work2, size(work2), rwork2, info2)
  call check('zgesvd', all(s1 == s2) .and. all(U1 == U2) .and. all(V1 == V2) .and. info1 == info2)

  C1 = H; C2 = H
  call zheev('V', 'U', n, C1, n, s1, work1, size(work1), rwork1, info1)
  call la_zheev('V', 'U', n, C2, n, s2, work2, size(work2), rwork2, info2)
  call check('zheev', all(C1 == C2) .and. all(s1 == s2) .and. info1 == info2)

  C1 = A; C2 = A
  call zgeev('V', 'V', n, C1, n, w1, U1, n, V1, n, work1, size(work1), rwork1, info1)
  call la_zgeev('V', 'V', n, C2, n, w2, U2, n, V2, n, work2, size(work2), rwork2, info2)
  call check('zgeev', all(w1 == w2) .and. all(U1 == U2) .and. all(V1 == V2) .and. info1 == info2)

  if (n_fail > 0) then
     write(*,*) 'linalg_mod test: ', n_fail, ' failure(s)'
     error stop 1
  endif
  write(*,*) 'linalg_mod test: all wrappers bitwise equal to the direct calls'

contains

  subroutine check(name, ok)
    character(len=*), intent(in) :: name
    logical, intent(in) :: ok
    if (.not. ok) then
       write(*,*) 'FAIL: ', name
       n_fail = n_fail + 1
    endif
  end subroutine check

  subroutine fill(M)
    complex(dp), intent(out) :: M(:,:)
    real(dp) :: re(size(M,1),size(M,2)), im(size(M,1),size(M,2))
    call random_number(re); call random_number(im)
    M = cmplx(re - 0.5d0, im - 0.5d0, kind=dp)
  end subroutine fill

  subroutine fillv(v)
    complex(dp), intent(out) :: v(:)
    real(dp) :: re(size(v)), im(size(v))
    call random_number(re); call random_number(im)
    v = cmplx(re - 0.5d0, im - 0.5d0, kind=dp)
  end subroutine fillv

end program test_linalg
