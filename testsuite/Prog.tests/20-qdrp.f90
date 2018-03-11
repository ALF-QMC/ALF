! compile with
! gfortran -Wall -std=f2003 -I ../../Prog_8/  -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ 18-ul-update-matrices.f90 ../../Prog_8/wrap_helpers.o ../../Prog_8/UDV_WRAP.o ../../Libraries/Modules/modules_90.a ../../../../lapack-3.6.1/liblapack.a -lblas

Program TESTQDRP
        Use QDRP_mod
implicit none

!interface
!        
!SUBROUTINE QDRP_decompose(Ndim, Mat, D, IPVT, TAU, WORK, LWORK)
!        Use QDRP_mod
!        Implicit None
!        INTEGER, intent(in) :: Ndim, IPVT(Ndim)
!        INTEGER, intent(inout) :: LWORK
!        COMPLEX (Kind=Kind(0.d0)) :: Mat(Ndim,Ndim)
!        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim), TAU(Ndim)
!        COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:) WORK
!        end Subroutine
!end interface
        COMPLEX(Kind=Kind(0.D0)), Dimension(:,:), allocatable :: A, TEST, TMP
        COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:) :: D, TAU, WORK
        Integer, allocatable, dimension(:) :: IPVT
        Integer :: Ndim, LWORK, info, i, j
        Complex(Kind=kind(0.D0)) :: alpha
        Logical :: FWD
        REAL(Kind=Kind(0.D0)) :: diff
        
        do ndim=2, 50, 4
        Allocate(A(ndim, ndim), D(ndim), Tau(ndim), IPVT(ndim), TEST(Ndim, ndim), TMP(ndim, ndim))
        do i = 1, ndim
        do j = 1, ndim
        A(i,j) = 1.D0/(i+j) ! Hilbert's matrix
        enddo
        enddo
        TEST = A
!        do i = 1, ndim
!        write (*,*) DBLE(TEST(i, :))
!        enddo
        call QDRP_decompose(ndim, ndim, A, D, IPVT, TAU, WORK, LWORK)
        TMP = A
        CALL ZUNGQR(ndim, ndim, ndim, A, ndim, tau, work, lwork, info)
        do i = 1, ndim
        do j = 1, ndim
        A(i,j) = A(i,j) * D(j)
        enddo
        enddo
        alpha = 1.D0
        call ztrmm('R', 'U', 'N', 'N', ndim, ndim, alpha, TMP, ndim, A, ndim)
        FWD = .false.
        call ZLAPMT(FWD, ndim, ndim, A, ndim, IPVT)
!        write(*,*) "------------------"
!        do i = 1, ndim
!        write (*,*) DBLE(A(i, :))
!        enddo
!        write(*,*) "----------------------------------------"
do i = 1, ndim
do j = 1, ndim
diff = abs(dble(test(i,j) - a(i,j)))
if (diff > max(abs(dble(test(i,j))), abs(dble(a(i,j)))) * 1D-14) then
write(*,*) "Error!", ndim, i, j, test(i,j), a(i,j)
endif
enddo
enddo
        deallocate(A, D, TAU, IPVT, test, tmp, WORK)
        enddo
end Program TESTQDRP
