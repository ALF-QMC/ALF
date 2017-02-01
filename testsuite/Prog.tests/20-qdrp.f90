! compile with
! gfortran -Wall -std=f2003 -I ../../Prog_8/  -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ 18-ul-update-matrices.f90 ../../Prog_8/wrap_helpers.o ../../Prog_8/UDV_WRAP.o ../../Libraries/Modules/modules_90.a ../../../../lapack-3.6.1/liblapack.a -lblas

Program TESTQDRP
implicit none
interface
        
SUBROUTINE QDRP_decompose(Ndim, Mat, D, IPVT, TAU, WORK, LWORK)
        Use QDRP_mod
        Implicit None
        INTEGER, intent(in) :: Ndim, IPVT(Ndim)
        INTEGER, intent(inout) :: LWORK
        COMPLEX (Kind=Kind(0.d0)) :: Mat(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim), TAU(Ndim)
        COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:) WORK
        end Subroutine
end interface
        COMPLEX(Kind=Kind(0.D0)), Dimension(:,:), allocatable :: A, TEST, TMP
        COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:) :: D, TAU, WORK
        Integer, allocatable, dimension(:) :: IPVT
        Integer :: Ndim, LWORK, info, i, j
        Complex(Kind=kind(0.D0)) :: alpha
        Logical :: FWD
        
        do ndim=2, 10, 5
        Allocate(A(ndim, ndim), D(ndim), Tau(ndim), IPVT(ndim), TEST(Ndim, ndim), TMP(ndim, ndim))
        TEST = A
        call QDRP_decompose(ndim, A, D, IPVT, TAU, WORK, LWORK)
        TMP = A
        CALL ZUNGQR(ndim, ndim, ndim, A, ndim, tau, work, lwork, info)
        do i = 1, ndim
        do j = 1, ndim
        A(i,j) = A(i,j) * D(j)
        enddo
        enddo
        alpha = 1.D0
        call ztrmm('R', 'U', 'N', 'N', ndim, ndim, alpha, TMP, ndim, A, ndim)
        FWD = .true.
        call ZLAPMT(FWD, ndim, ndim, A, ndim, IPVT)
        
        deallocate(A, D, TAU IPVT, test, tmp)
        enddo
end Program TESTQDRP
