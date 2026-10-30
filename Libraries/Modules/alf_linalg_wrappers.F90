!  Copyright (C) 2026 The ALF project
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

! Wrappers around all BLAS and LAPACK routines called by ALF, see alf_linalg_mod.F90.
! They are external procedures with implicit interfaces, just like BLAS/LAPACK.

#ifdef ALF_GPU
include "mkl_blas_omp_offload_lp64.f90"
include "mkl_lapack_omp_offload_lp64.f90"
#endif

!====================================================================
! BLAS level 3
!====================================================================

subroutine ALF_ZGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldc
  complex(dp), intent(in) :: alpha, beta, a(*), b(*)
  complex(dp), intent(inout) :: c(*)
#ifdef ALF_GPU
  integer(int64) :: na, nb, nc
  if (use_gpu(ID_ZGEMM, int(m, int64)*n*k)) then
    if (is_char(transa, 'N')) then; na = mspan(lda, m, k); else; na = mspan(lda, k, m); endif
    if (is_char(transb, 'N')) then; nb = mspan(ldb, k, n); else; nb = mspan(ldb, n, k); endif
    nc = mspan(ldc, m, n)
    !$omp target data map(to: a(1:na), b(1:nb)) map(tofrom: c(1:nc))
    !$omp dispatch
    call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    !$omp end target data
    return
  endif
#endif
  call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine ALF_ZGEMM

subroutine ALF_DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldc
  real(dp), intent(in) :: alpha, beta, a(*), b(*)
  real(dp), intent(inout) :: c(*)
#ifdef ALF_GPU
  integer(int64) :: na, nb, nc
  if (use_gpu(ID_DGEMM, int(m, int64)*n*k)) then
    if (is_char(transa, 'N')) then; na = mspan(lda, m, k); else; na = mspan(lda, k, m); endif
    if (is_char(transb, 'N')) then; nb = mspan(ldb, k, n); else; nb = mspan(ldb, n, k); endif
    nc = mspan(ldc, m, n)
    !$omp target data map(to: a(1:na), b(1:nb)) map(tofrom: c(1:nc))
    !$omp dispatch
    call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    !$omp end target data
    return
  endif
#endif
  call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine ALF_DGEMM

subroutine ALF_ZHEMM(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: side, uplo
  integer, intent(in) :: m, n, lda, ldb, ldc
  complex(dp), intent(in) :: alpha, beta, a(*), b(*)
  complex(dp), intent(inout) :: c(*)
#ifdef ALF_GPU
  integer(int64) :: na, nb, nc
  integer :: ka
  if (is_char(side, 'L')) then; ka = m; else; ka = n; endif
  if (use_gpu(ID_ZHEMM, int(m, int64)*n*ka)) then
    na = mspan(lda, ka, ka)
    nb = mspan(ldb, m, n)
    nc = mspan(ldc, m, n)
    !$omp target data map(to: a(1:na), b(1:nb)) map(tofrom: c(1:nc))
    !$omp dispatch
    call zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    !$omp end target data
    return
  endif
#endif
  call zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine ALF_ZHEMM

subroutine ALF_ZTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: side, uplo, transa, diag
  integer, intent(in) :: m, n, lda, ldb
  complex(dp), intent(in) :: alpha, a(*)
  complex(dp), intent(inout) :: b(*)
#ifdef ALF_GPU
  integer(int64) :: na, nb
  integer :: ka
  if (is_char(side, 'L')) then; ka = m; else; ka = n; endif
  if (use_gpu(ID_ZTRSM, int(m, int64)*n*ka)) then
    na = mspan(lda, ka, ka)
    nb = mspan(ldb, m, n)
    !$omp target data map(to: a(1:na)) map(tofrom: b(1:nb))
    !$omp dispatch
    call ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    !$omp end target data
    return
  endif
#endif
  call ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
end subroutine ALF_ZTRSM

subroutine ALF_ZTRMM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: side, uplo, transa, diag
  integer, intent(in) :: m, n, lda, ldb
  complex(dp), intent(in) :: alpha, a(*)
  complex(dp), intent(inout) :: b(*)
#ifdef ALF_GPU
  integer(int64) :: na, nb
  integer :: ka
  if (is_char(side, 'L')) then; ka = m; else; ka = n; endif
  if (use_gpu(ID_ZTRMM, int(m, int64)*n*ka)) then
    na = mspan(lda, ka, ka)
    nb = mspan(ldb, m, n)
    !$omp target data map(to: a(1:na)) map(tofrom: b(1:nb))
    !$omp dispatch
    call ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    !$omp end target data
    return
  endif
#endif
  call ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
end subroutine ALF_ZTRMM

!====================================================================
! BLAS level 2
!====================================================================

subroutine ALF_ZGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: trans
  integer, intent(in) :: m, n, lda, incx, incy
  complex(dp), intent(in) :: alpha, beta, a(*), x(*)
  complex(dp), intent(inout) :: y(*)
#ifdef ALF_GPU
  integer(int64) :: na, nx, ny
  if (use_gpu(ID_ZGEMV, int(m, int64)*n)) then
    na = mspan(lda, m, n)
    if (is_char(trans, 'N')) then
      nx = vspan(n, incx); ny = vspan(m, incy)
    else
      nx = vspan(m, incx); ny = vspan(n, incy)
    endif
    !$omp target data map(to: a(1:na), x(1:nx)) map(tofrom: y(1:ny))
    !$omp dispatch
    call zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    !$omp end target data
    return
  endif
#endif
  call zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine ALF_ZGEMV

subroutine ALF_DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: trans
  integer, intent(in) :: m, n, lda, incx, incy
  real(dp), intent(in) :: alpha, beta, a(*), x(*)
  real(dp), intent(inout) :: y(*)
#ifdef ALF_GPU
  integer(int64) :: na, nx, ny
  if (use_gpu(ID_DGEMV, int(m, int64)*n)) then
    na = mspan(lda, m, n)
    if (is_char(trans, 'N')) then
      nx = vspan(n, incx); ny = vspan(m, incy)
    else
      nx = vspan(m, incx); ny = vspan(n, incy)
    endif
    !$omp target data map(to: a(1:na), x(1:nx)) map(tofrom: y(1:ny))
    !$omp dispatch
    call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    !$omp end target data
    return
  endif
#endif
  call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine ALF_DGEMV

subroutine ALF_ZGERU(m, n, alpha, x, incx, y, incy, a, lda)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: m, n, incx, incy, lda
  complex(dp), intent(in) :: alpha, x(*), y(*)
  complex(dp), intent(inout) :: a(*)
#ifdef ALF_GPU
  integer(int64) :: na, nx, ny
  if (use_gpu(ID_ZGERU, int(m, int64)*n)) then
    na = mspan(lda, m, n)
    nx = vspan(m, incx)
    ny = vspan(n, incy)
    !$omp target data map(to: x(1:nx), y(1:ny)) map(tofrom: a(1:na))
    !$omp dispatch
    call zgeru(m, n, alpha, x, incx, y, incy, a, lda)
    !$omp end target data
    return
  endif
#endif
  call zgeru(m, n, alpha, x, incx, y, incy, a, lda)
end subroutine ALF_ZGERU

!====================================================================
! BLAS level 1
!====================================================================

subroutine ALF_ZSCAL(n, za, zx, incx)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: n, incx
  complex(dp), intent(in) :: za
  complex(dp), intent(inout) :: zx(*)
#ifdef ALF_GPU
  integer(int64) :: nx
  if (use_gpu(ID_ZSCAL, int(n, int64))) then
    nx = vspan(n, incx)
    !$omp target data map(tofrom: zx(1:nx))
    !$omp dispatch
    call zscal(n, za, zx, incx)
    !$omp end target data
    return
  endif
#endif
  call zscal(n, za, zx, incx)
end subroutine ALF_ZSCAL

subroutine ALF_ZCOPY(n, zx, incx, zy, incy)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: n, incx, incy
  complex(dp), intent(in) :: zx(*)
  complex(dp), intent(inout) :: zy(*)
#ifdef ALF_GPU
  integer(int64) :: nx, ny
  if (use_gpu(ID_ZCOPY, int(n, int64))) then
    nx = vspan(n, incx)
    ny = vspan(n, incy)
    ! tofrom, since a strided zy has gaps which must not be overwritten
    !$omp target data map(to: zx(1:nx)) map(tofrom: zy(1:ny))
    !$omp dispatch
    call zcopy(n, zx, incx, zy, incy)
    !$omp end target data
    return
  endif
#endif
  call zcopy(n, zx, incx, zy, incy)
end subroutine ALF_ZCOPY

subroutine ALF_ZSWAP(n, zx, incx, zy, incy)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: n, incx, incy
  complex(dp), intent(inout) :: zx(*), zy(*)
#ifdef ALF_GPU
  integer(int64) :: nx, ny
  if (use_gpu(ID_ZSWAP, int(n, int64))) then
    nx = vspan(n, incx)
    ny = vspan(n, incy)
    !$omp target data map(tofrom: zx(1:nx), zy(1:ny))
    !$omp dispatch
    call zswap(n, zx, incx, zy, incy)
    !$omp end target data
    return
  endif
#endif
  call zswap(n, zx, incx, zy, incy)
end subroutine ALF_ZSWAP

function ALF_DZNRM2(n, x, incx)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_blas_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  real(dp) :: ALF_DZNRM2
  integer, intent(in) :: n, incx
  complex(dp), intent(in) :: x(*)
#ifdef ALF_GPU
  integer(int64) :: nx
  if (use_gpu(ID_DZNRM2, int(n, int64))) then
    nx = vspan(n, incx)
    !$omp target data map(to: x(1:nx))
    !$omp dispatch
    ALF_DZNRM2 = dznrm2(n, x, incx)
    !$omp end target data
    return
  endif
#else
  real(dp), external :: dznrm2
#endif
  ALF_DZNRM2 = dznrm2(n, x, incx)
end function ALF_DZNRM2

!====================================================================
! LAPACK
! The offloaded calls query and allocate their own workspace, since the
! workspace required on the GPU may differ from the one on the CPU.
! Workspace queries (lwork = -1) always run on the CPU.
!====================================================================

subroutine ALF_ZGETRF(m, n, a, lda, ipiv, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: m, n, lda
  complex(dp), intent(inout) :: a(*)
  integer, intent(out) :: ipiv(*), info
#ifdef ALF_GPU
  integer(int64) :: na
  if (use_gpu(ID_ZGETRF, int(m, int64)*n*min(m, n))) then
    na = mspan(lda, m, n)
    !$omp target data map(tofrom: a(1:na)) map(from: ipiv(1:min(m, n)), info)
    !$omp dispatch
    call zgetrf(m, n, a, lda, ipiv, info)
    !$omp end target data
    return
  endif
#endif
  call zgetrf(m, n, a, lda, ipiv, info)
end subroutine ALF_ZGETRF

subroutine ALF_DGETRF(m, n, a, lda, ipiv, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: m, n, lda
  real(dp), intent(inout) :: a(*)
  integer, intent(out) :: ipiv(*), info
#ifdef ALF_GPU
  integer(int64) :: na
  if (use_gpu(ID_DGETRF, int(m, int64)*n*min(m, n))) then
    na = mspan(lda, m, n)
    !$omp target data map(tofrom: a(1:na)) map(from: ipiv(1:min(m, n)), info)
    !$omp dispatch
    call dgetrf(m, n, a, lda, ipiv, info)
    !$omp end target data
    return
  endif
#endif
  call dgetrf(m, n, a, lda, ipiv, info)
end subroutine ALF_DGETRF

subroutine ALF_ZGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
  complex(dp), intent(in) :: a(*)
  complex(dp), intent(inout) :: b(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na, nb
  if (use_gpu(ID_ZGETRS, int(n, int64)*n*nrhs)) then
    na = mspan(lda, n, n)
    nb = mspan(ldb, n, nrhs)
    !$omp target data map(to: a(1:na), ipiv(1:n)) map(tofrom: b(1:nb)) map(from: info)
    !$omp dispatch
    call zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    !$omp end target data
    return
  endif
#endif
  call zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
end subroutine ALF_ZGETRS

subroutine ALF_ZGETRI(n, a, lda, ipiv, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: n, lda, lwork, ipiv(*)
  complex(dp), intent(inout) :: a(*), work(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  complex(dp) :: wq(1)
  complex(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_ZGETRI, int(n, int64)**3)) then
    na = mspan(lda, n, n)
    !$omp target data map(tofrom: a(1:na)) map(to: ipiv(1:n)) map(from: info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call zgetri(n, a, lda, ipiv, wq, -1, info)
    !$omp end target data
    lw = max(1, int(real(wq(1))))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call zgetri(n, a, lda, ipiv, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call zgetri(n, a, lda, ipiv, work, lwork, info)
end subroutine ALF_ZGETRI

subroutine ALF_DGETRI(n, a, lda, ipiv, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: n, lda, lwork, ipiv(*)
  real(dp), intent(inout) :: a(*), work(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  real(dp) :: wq(1)
  real(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_DGETRI, int(n, int64)**3)) then
    na = mspan(lda, n, n)
    !$omp target data map(tofrom: a(1:na)) map(to: ipiv(1:n)) map(from: info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call dgetri(n, a, lda, ipiv, wq, -1, info)
    !$omp end target data
    lw = max(1, int(wq(1)))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call dgetri(n, a, lda, ipiv, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call dgetri(n, a, lda, ipiv, work, lwork, info)
end subroutine ALF_DGETRI

subroutine ALF_ZGEQRF(m, n, a, lda, tau, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: m, n, lda, lwork
  complex(dp), intent(inout) :: a(*), tau(*), work(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  complex(dp) :: wq(1)
  complex(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_ZGEQRF, int(m, int64)*n*min(m, n))) then
    na = mspan(lda, m, n)
    !$omp target data map(tofrom: a(1:na)) map(from: tau(1:max(1, min(m, n))), info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call zgeqrf(m, n, a, lda, tau, wq, -1, info)
    !$omp end target data
    lw = max(1, int(real(wq(1))))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call zgeqrf(m, n, a, lda, tau, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call zgeqrf(m, n, a, lda, tau, work, lwork, info)
end subroutine ALF_ZGEQRF

subroutine ALF_DGEQRF(m, n, a, lda, tau, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: m, n, lda, lwork
  real(dp), intent(inout) :: a(*), tau(*), work(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  real(dp) :: wq(1)
  real(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_DGEQRF, int(m, int64)*n*min(m, n))) then
    na = mspan(lda, m, n)
    !$omp target data map(tofrom: a(1:na)) map(from: tau(1:max(1, min(m, n))), info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call dgeqrf(m, n, a, lda, tau, wq, -1, info)
    !$omp end target data
    lw = max(1, int(wq(1)))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call dgeqrf(m, n, a, lda, tau, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call dgeqrf(m, n, a, lda, tau, work, lwork, info)
end subroutine ALF_DGEQRF

subroutine ALF_ZUNGQR(m, n, k, a, lda, tau, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: m, n, k, lda, lwork
  complex(dp), intent(inout) :: a(*), work(*)
  complex(dp), intent(in) :: tau(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  complex(dp) :: wq(1)
  complex(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_ZUNGQR, int(m, int64)*n*k)) then
    na = mspan(lda, m, n)
    !$omp target data map(tofrom: a(1:na)) map(to: tau(1:max(1, k))) map(from: info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call zungqr(m, n, k, a, lda, tau, wq, -1, info)
    !$omp end target data
    lw = max(1, int(real(wq(1))))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call zungqr(m, n, k, a, lda, tau, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call zungqr(m, n, k, a, lda, tau, work, lwork, info)
end subroutine ALF_ZUNGQR

subroutine ALF_DORGQR(m, n, k, a, lda, tau, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  integer, intent(in) :: m, n, k, lda, lwork
  real(dp), intent(inout) :: a(*), work(*)
  real(dp), intent(in) :: tau(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  real(dp) :: wq(1)
  real(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_DORGQR, int(m, int64)*n*k)) then
    na = mspan(lda, m, n)
    !$omp target data map(tofrom: a(1:na)) map(to: tau(1:max(1, k))) map(from: info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call dorgqr(m, n, k, a, lda, tau, wq, -1, info)
    !$omp end target data
    lw = max(1, int(wq(1)))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call dorgqr(m, n, k, a, lda, tau, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call dorgqr(m, n, k, a, lda, tau, work, lwork, info)
end subroutine ALF_DORGQR

subroutine ALF_ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: side, trans
  integer, intent(in) :: m, n, k, lda, ldc, lwork
  complex(dp), intent(in) :: a(*), tau(*)
  complex(dp), intent(inout) :: c(*), work(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na, nc
  complex(dp) :: wq(1)
  complex(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_ZUNMQR, int(m, int64)*n*k)) then
    if (is_char(side, 'L')) then; na = mspan(lda, m, k); else; na = mspan(lda, n, k); endif
    nc = mspan(ldc, m, n)
    !$omp target data map(to: a(1:na), tau(1:max(1, k))) map(tofrom: c(1:nc)) map(from: info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, wq, -1, info)
    !$omp end target data
    lw = max(1, int(real(wq(1))))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
end subroutine ALF_ZUNMQR

subroutine ALF_ZHEEV(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: jobz, uplo
  integer, intent(in) :: n, lda, lwork
  complex(dp), intent(inout) :: a(*), work(*)
  real(dp), intent(out) :: w(*)
  real(dp), intent(inout) :: rwork(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  complex(dp) :: wq(1)
  complex(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_ZHEEV, int(n, int64)**3)) then
    na = mspan(lda, n, n)
    !$omp target data map(tofrom: a(1:na)) map(from: w(1:max(1, n)), info) &
    !$omp&            map(alloc: rwork(1:max(1, 3*n - 2)))
    !$omp target data map(from: wq)
    !$omp dispatch
    call zheev(jobz, uplo, n, a, lda, w, wq, -1, rwork, info)
    !$omp end target data
    lw = max(1, int(real(wq(1))))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call zheev(jobz, uplo, n, a, lda, w, gwork, lw, rwork, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
end subroutine ALF_ZHEEV

subroutine ALF_DSYEV(jobz, uplo, n, a, lda, w, work, lwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: jobz, uplo
  integer, intent(in) :: n, lda, lwork
  real(dp), intent(inout) :: a(*), work(*)
  real(dp), intent(out) :: w(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na
  real(dp) :: wq(1)
  real(dp), allocatable :: gwork(:)
  integer :: lw
  if (lwork /= -1 .and. use_gpu(ID_DSYEV, int(n, int64)**3)) then
    na = mspan(lda, n, n)
    !$omp target data map(tofrom: a(1:na)) map(from: w(1:max(1, n)), info)
    !$omp target data map(from: wq)
    !$omp dispatch
    call dsyev(jobz, uplo, n, a, lda, w, wq, -1, info)
    !$omp end target data
    lw = max(1, int(wq(1)))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call dsyev(jobz, uplo, n, a, lda, w, gwork, lw, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
end subroutine ALF_DSYEV

subroutine ALF_ZGESVD(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
  use alf_linalg_mod
#ifdef ALF_GPU
  use onemkl_lapack_omp_offload_lp64
  use iso_fortran_env, only: int64
#endif
  implicit none
  character, intent(in) :: jobu, jobvt
  integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
  complex(dp), intent(inout) :: a(*), u(*), vt(*), work(*)
  real(dp), intent(out) :: s(*)
  real(dp), intent(inout) :: rwork(*)
  integer, intent(out) :: info
#ifdef ALF_GPU
  integer(int64) :: na, nu, nvt
  complex(dp) :: wq(1)
  complex(dp), allocatable :: gwork(:)
  integer :: lw, mn
  mn = min(m, n)
  if (lwork /= -1 .and. use_gpu(ID_ZGESVD, int(m, int64)*n*mn)) then
    na = mspan(lda, m, n)
    nu = 1
    if (is_char(jobu, 'A')) nu = mspan(ldu, m, m)
    if (is_char(jobu, 'S')) nu = mspan(ldu, m, mn)
    nvt = 1
    if (is_char(jobvt, 'A')) nvt = mspan(ldvt, n, n)
    if (is_char(jobvt, 'S')) nvt = mspan(ldvt, mn, n)
    ! u and vt tofrom, since they may have gaps which must not be overwritten
    !$omp target data map(tofrom: a(1:na), u(1:nu), vt(1:nvt)) map(from: s(1:max(1, mn)), info) &
    !$omp&            map(alloc: rwork(1:max(1, 5*mn)))
    !$omp target data map(from: wq)
    !$omp dispatch
    call zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, wq, -1, rwork, info)
    !$omp end target data
    lw = max(1, int(real(wq(1))))
    allocate(gwork(lw))
    !$omp target data map(alloc: gwork)
    !$omp dispatch
    call zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, gwork, lw, rwork, info)
    !$omp end target data
    !$omp end target data
    deallocate(gwork)
    return
  endif
#endif
  call zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
end subroutine ALF_ZGESVD

!====================================================================
! Pass-through wrappers, no OpenMP offload support in oneMKL
!====================================================================

subroutine ALF_ZGEQP3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
  use alf_linalg_mod, only: dp
  implicit none
  integer, intent(in) :: m, n, lda, lwork
  complex(dp), intent(inout) :: a(*), tau(*), work(*)
  integer, intent(inout) :: jpvt(*)
  real(dp), intent(inout) :: rwork(*)
  integer, intent(out) :: info
  call zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
end subroutine ALF_ZGEQP3

subroutine ALF_ZGEEV(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
  use alf_linalg_mod, only: dp
  implicit none
  character, intent(in) :: jobvl, jobvr
  integer, intent(in) :: n, lda, ldvl, ldvr, lwork
  complex(dp), intent(inout) :: a(*), w(*), vl(*), vr(*), work(*)
  real(dp), intent(inout) :: rwork(*)
  integer, intent(out) :: info
  call zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
end subroutine ALF_ZGEEV

subroutine ALF_ZLACPY(uplo, m, n, a, lda, b, ldb)
  use alf_linalg_mod, only: dp
  implicit none
  character, intent(in) :: uplo
  integer, intent(in) :: m, n, lda, ldb
  complex(dp), intent(in) :: a(*)
  complex(dp), intent(inout) :: b(*)
  call zlacpy(uplo, m, n, a, lda, b, ldb)
end subroutine ALF_ZLACPY

subroutine ALF_DLACPY(uplo, m, n, a, lda, b, ldb)
  use alf_linalg_mod, only: dp
  implicit none
  character, intent(in) :: uplo
  integer, intent(in) :: m, n, lda, ldb
  real(dp), intent(in) :: a(*)
  real(dp), intent(inout) :: b(*)
  call dlacpy(uplo, m, n, a, lda, b, ldb)
end subroutine ALF_DLACPY

subroutine ALF_ZLASET(uplo, m, n, alpha, beta, a, lda)
  use alf_linalg_mod, only: dp
  implicit none
  character, intent(in) :: uplo
  integer, intent(in) :: m, n, lda
  complex(dp), intent(in) :: alpha, beta
  complex(dp), intent(inout) :: a(*)
  call zlaset(uplo, m, n, alpha, beta, a, lda)
end subroutine ALF_ZLASET

subroutine ALF_ZLAPMR(forwrd, m, n, x, ldx, k)
  use alf_linalg_mod, only: dp
  implicit none
  logical, intent(in) :: forwrd
  integer, intent(in) :: m, n, ldx
  complex(dp), intent(inout) :: x(*)
  integer, intent(inout) :: k(*)
  call zlapmr(forwrd, m, n, x, ldx, k)
end subroutine ALF_ZLAPMR

subroutine ALF_ZLAPMT(forwrd, m, n, x, ldx, k)
  use alf_linalg_mod, only: dp
  implicit none
  logical, intent(in) :: forwrd
  integer, intent(in) :: m, n, ldx
  complex(dp), intent(inout) :: x(*)
  integer, intent(inout) :: k(*)
  call zlapmt(forwrd, m, n, x, ldx, k)
end subroutine ALF_ZLAPMT

subroutine ALF_ZLACRM(m, n, a, lda, b, ldb, c, ldc, rwork)
  use alf_linalg_mod, only: dp
  implicit none
  integer, intent(in) :: m, n, lda, ldb, ldc
  complex(dp), intent(in) :: a(*)
  real(dp), intent(in) :: b(*)
  complex(dp), intent(inout) :: c(*)
  real(dp), intent(inout) :: rwork(*)
  call zlacrm(m, n, a, lda, b, ldb, c, ldc, rwork)
end subroutine ALF_ZLACRM

subroutine ALF_ZLARCM(m, n, a, lda, b, ldb, c, ldc, rwork)
  use alf_linalg_mod, only: dp
  implicit none
  integer, intent(in) :: m, n, lda, ldb, ldc
  real(dp), intent(in) :: a(*)
  complex(dp), intent(in) :: b(*)
  complex(dp), intent(inout) :: c(*)
  real(dp), intent(inout) :: rwork(*)
  call zlarcm(m, n, a, lda, b, ldb, c, ldc, rwork)
end subroutine ALF_ZLARCM
