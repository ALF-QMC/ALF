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


!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Thin wrappers around the BLAS and LAPACK routines used by the QMC core.
!>
!> @details
!> Each la_<name> has exactly the argument list of the reference routine
!> <name> and does nothing but call it, so results are identical to a direct
!> call. The QMC code calls these wrappers instead of the vendor library, so
!> that a GPU backend (e.g. OpenMP `dispatch` to oneMKL) changes only this
!> module. Expected placement in such a backend:
!>  - device: the BLAS routines, zgetrf/zgetrs/zgetri;
!>  - host only (no GPU version in oneMKL): zgeqp3, zungqr, zgesvd, zheev, zgeev;
!>  - auxiliary (zlaset, zlapmr, zlapmt): simple loops, to be offloaded as kernels.
!> Unlike the implicit interface of an external routine, the interfaces
!> here are checked by the compiler (types, ranks, intents).
!--------------------------------------------------------------------
module linalg_mod
  implicit none
  private
  public :: la_zgemm, la_zhemm, la_ztrsm, la_ztrmm, la_zgemv, la_zgeru, la_zscal, la_zcopy, &
       &    la_zlacpy, la_zgetrf, la_zgetrs, la_zgetri, la_zgeqrf, la_zgeqp3, la_zungqr, la_zunmqr, &
       &    la_zlaset, la_zlapmr, la_zlapmt, la_zgesvd, la_zheev, la_zgeev
contains
  !> C = alpha op(A) op(B) + beta C (BLAS 3)
  subroutine la_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    character(len=1), intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(kind=kind(0.d0)), intent(in) :: alpha, beta
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*), b(ldb,*)
    complex(kind=kind(0.d0)), intent(inout) :: c(ldc,*)
    external :: zgemm
    call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  end subroutine la_zgemm

  !> C = alpha A B + beta C or alpha B A + beta C, A Hermitian (BLAS 3)
  subroutine la_zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    character(len=1), intent(in) :: side, uplo
    integer, intent(in) :: m, n, lda, ldb, ldc
    complex(kind=kind(0.d0)), intent(in) :: alpha, beta
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*), b(ldb,*)
    complex(kind=kind(0.d0)), intent(inout) :: c(ldc,*)
    external :: zhemm
    call zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
  end subroutine la_zhemm

  !> Triangular solve with several right-hand sides (BLAS 3)
  subroutine la_ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    character(len=1), intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    complex(kind=kind(0.d0)), intent(in) :: alpha
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*)
    complex(kind=kind(0.d0)), intent(inout) :: b(ldb,*)
    external :: ztrsm
    call ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
  end subroutine la_ztrsm

  !> Triangular matrix product (BLAS 3)
  subroutine la_ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    character(len=1), intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    complex(kind=kind(0.d0)), intent(in) :: alpha
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*)
    complex(kind=kind(0.d0)), intent(inout) :: b(ldb,*)
    external :: ztrmm
    call ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
  end subroutine la_ztrmm

  !> y = alpha op(A) x + beta y (BLAS 2)
  subroutine la_zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    character(len=1), intent(in) :: trans
    integer, intent(in) :: m, n, lda, incx, incy
    complex(kind=kind(0.d0)), intent(in) :: alpha, beta
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*), x(*)
    complex(kind=kind(0.d0)), intent(inout) :: y(*)
    external :: zgemv
    call zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
  end subroutine la_zgemv

  !> A = alpha x y^T + A (BLAS 2)
  subroutine la_zgeru(m, n, alpha, x, incx, y, incy, a, lda)
    integer, intent(in) :: m, n, incx, incy, lda
    complex(kind=kind(0.d0)), intent(in) :: alpha
    complex(kind=kind(0.d0)), intent(in) :: x(*), y(*)
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*)
    external :: zgeru
    call zgeru(m, n, alpha, x, incx, y, incy, a, lda)
  end subroutine la_zgeru

  !> x = alpha x (BLAS 1)
  subroutine la_zscal(n, alpha, x, incx)
    integer, intent(in) :: n, incx
    complex(kind=kind(0.d0)), intent(in) :: alpha
    complex(kind=kind(0.d0)), intent(inout) :: x(*)
    external :: zscal
    call zscal(n, alpha, x, incx)
  end subroutine la_zscal

  !> y = x (BLAS 1)
  subroutine la_zcopy(n, x, incx, y, incy)
    integer, intent(in) :: n, incx, incy
    complex(kind=kind(0.d0)), intent(in) :: x(*)
    complex(kind=kind(0.d0)), intent(inout) :: y(*)
    external :: zcopy
    call zcopy(n, x, incx, y, incy)
  end subroutine la_zcopy

  !> B = A, all or a triangle (LAPACK auxiliary)
  subroutine la_zlacpy(uplo, m, n, a, lda, b, ldb)
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: m, n, lda, ldb
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*)
    complex(kind=kind(0.d0)), intent(inout) :: b(ldb,*)
    external :: zlacpy
    call zlacpy(uplo, m, n, a, lda, b, ldb)
  end subroutine la_zlacpy

  !> LU factorisation with partial pivoting
  subroutine la_zgetrf(m, n, a, lda, ipiv, info)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*), info
    external :: zgetrf
    call zgetrf(m, n, a, lda, ipiv, info)
  end subroutine la_zgetrf

  !> Solve with the LU factors of zgetrf
  subroutine la_zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    character(len=1), intent(in) :: trans
    integer, intent(in) :: n, nrhs, lda, ldb
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*)
    integer, intent(in) :: ipiv(*)
    complex(kind=kind(0.d0)), intent(inout) :: b(ldb,*)
    integer, intent(out) :: info
    external :: zgetrs
    call zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
  end subroutine la_zgetrs

  !> Inverse from the LU factors of zgetrf
  subroutine la_zgetri(n, a, lda, ipiv, work, lwork, info)
    integer, intent(in) :: n, lda, lwork
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*)
    integer, intent(in) :: ipiv(*)
    complex(kind=kind(0.d0)), intent(inout) :: work(*)
    integer, intent(out) :: info
    external :: zgetri
    call zgetri(n, a, lda, ipiv, work, lwork, info)
  end subroutine la_zgetri

  !> QR factorisation
  subroutine la_zgeqrf(m, n, a, lda, tau, work, lwork, info)
    integer, intent(in) :: m, n, lda, lwork
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*)
    complex(kind=kind(0.d0)), intent(inout) :: tau(*), work(*)
    integer, intent(out) :: info
    external :: zgeqrf
    call zgeqrf(m, n, a, lda, tau, work, lwork, info)
  end subroutine la_zgeqrf

  !> QR factorisation with column pivoting (host only: no GPU version in oneMKL)
  subroutine la_zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
    integer, intent(in) :: m, n, lda, lwork
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*)
    integer, intent(inout) :: jpvt(*)
    complex(kind=kind(0.d0)), intent(inout) :: tau(*), work(*)
    real(kind=kind(0.d0)), intent(inout) :: rwork(*)
    integer, intent(out) :: info
    external :: zgeqp3
    call zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
  end subroutine la_zgeqp3

  !> Q of a QR factorisation (host only in oneMKL)
  subroutine la_zungqr(m, n, k, a, lda, tau, work, lwork, info)
    integer, intent(in) :: m, n, k, lda, lwork
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*)
    complex(kind=kind(0.d0)), intent(in) :: tau(*)
    complex(kind=kind(0.d0)), intent(inout) :: work(*)
    integer, intent(out) :: info
    external :: zungqr
    call zungqr(m, n, k, a, lda, tau, work, lwork, info)
  end subroutine la_zungqr

  !> C = op(Q) C or C op(Q), Q from a QR factorisation
  subroutine la_zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    character(len=1), intent(in) :: side, trans
    integer, intent(in) :: m, n, k, lda, ldc, lwork
    complex(kind=kind(0.d0)), intent(in) :: a(lda,*), tau(*)
    complex(kind=kind(0.d0)), intent(inout) :: c(ldc,*), work(*)
    integer, intent(out) :: info
    external :: zunmqr
    call zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
  end subroutine la_zunmqr

  !> A = beta on the diagonal and alpha off the diagonal (LAPACK auxiliary)
  subroutine la_zlaset(uplo, m, n, alpha, beta, a, lda)
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), intent(in) :: alpha, beta
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*)
    external :: zlaset
    call zlaset(uplo, m, n, alpha, beta, a, lda)
  end subroutine la_zlaset

  !> Permute the rows of X by K (LAPACK auxiliary; K is restored on exit)
  subroutine la_zlapmr(forwrd, m, n, x, ldx, k)
    logical, intent(in) :: forwrd
    integer, intent(in) :: m, n, ldx
    complex(kind=kind(0.d0)), intent(inout) :: x(ldx,*)
    integer, intent(inout) :: k(*)
    external :: zlapmr
    call zlapmr(forwrd, m, n, x, ldx, k)
  end subroutine la_zlapmr

  !> Permute the columns of X by K (LAPACK auxiliary; K is restored on exit)
  subroutine la_zlapmt(forwrd, m, n, x, ldx, k)
    logical, intent(in) :: forwrd
    integer, intent(in) :: m, n, ldx
    complex(kind=kind(0.d0)), intent(inout) :: x(ldx,*)
    integer, intent(inout) :: k(*)
    external :: zlapmt
    call zlapmt(forwrd, m, n, x, ldx, k)
  end subroutine la_zlapmt

  !> Singular value decomposition A = U S V^H (host only in oneMKL)
  subroutine la_zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
    character(len=1), intent(in) :: jobu, jobvt
    integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*), u(ldu,*), vt(ldvt,*), work(*)
    real(kind=kind(0.d0)), intent(inout) :: s(*), rwork(*)
    integer, intent(out) :: info
    external :: zgesvd
    call zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
  end subroutine la_zgesvd

  !> Eigenvalues and eigenvectors of a Hermitian matrix (host only in oneMKL)
  subroutine la_zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
    character(len=1), intent(in) :: jobz, uplo
    integer, intent(in) :: n, lda, lwork
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*), work(*)
    real(kind=kind(0.d0)), intent(inout) :: w(*), rwork(*)
    integer, intent(out) :: info
    external :: zheev
    call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
  end subroutine la_zheev

  !> Eigenvalues and eigenvectors of a general matrix (host only in oneMKL)
  subroutine la_zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    character(len=1), intent(in) :: jobvl, jobvr
    integer, intent(in) :: n, lda, ldvl, ldvr, lwork
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,*), w(*), vl(ldvl,*), vr(ldvr,*), work(*)
    real(kind=kind(0.d0)), intent(inout) :: rwork(*)
    integer, intent(out) :: info
    external :: zgeev
    call zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
  end subroutine la_zgeev

end module linalg_mod
