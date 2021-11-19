! compile with
! gfortran -I ../../Libraries/libmscbdecomp/ -L ../../Libraries/libmscbdecomp/ 1-HomogeneousDiag.F90 -lmscbdecomp

Program HomogeneousExpTest

  Use HomogeneousSingleColExp_mod

        COMPLEX (KIND=8) :: myx
        
        type(HomogeneousSingleColExp) :: test
        integer, parameter :: nredges = 6
        integer, parameter :: ndim = 24
        Type(Node) :: nodes(nredges)
        real(kind=kind(0.D0)) :: weight, sumdiag, sumoff
        real(kind=kind(0.D0)), allocatable :: mys(:)
        COMPLEX (KIND=8), DIMENSION(ndim, ndim) :: mat
        integer :: i
        
        allocate(mys(ndim))
        nodes(1)%x=1
        nodes(1)%y=3
        nodes(1)%axy=0.1

        nodes(2)%x=5
        nodes(2)%y=7
        nodes(2)%axy=0.1
        
        nodes(3)%x=9
        nodes(3)%y=11
        nodes(3)%axy=0.1
        
        nodes(4)%x=13
        nodes(4)%y=15
        nodes(4)%axy=0.1

        nodes(5)%x=17
        nodes(5)%y=19
        nodes(5)%axy=0.1

        nodes(6)%x=21
        nodes(6)%y=23
        nodes(6)%axy=0.1
        
        weight = 1.0
        mys = 0.0 ! initialize chemical potential to zero
        mys(1) = 0.2
        mys(3) = 0.2
        mys(5) = 0.4
        mys(7) = 0.4
        mys(9) = 0.6
        mys(11) = 0.6
        mys(13) = 0.8
        mys(15) = 0.8
        mys(17) = 1
        mys(19) = 1
        mys(21) = 1.2
        mys(23) = 1.2

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo

        call test%init(nodes, nredges, mys, weight)
        call test%lmult(mat)
        call test%lmultinv(mat)
        
        ! test for Trace(mat) = ndim
        sumdiag = 0
        sumoff = 0
        do i = 1, ndim
            sumdiag = sumdiag + DBLE(mat(i,i))
        enddo
        do i = 1, ndim-3
            sumoff = sumoff + DBLE(mat(i,i+2))
        enddo
        write (*,*) sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-15) then
        ERROR STOP 2
        endif
        if (abs(sumoff) > 1E-15) then !FIXME: this limit is a bit scale less...
        ERROR STOP 4
        endif

        call test%rmult(mat)
        call test%rmultinv(mat)
        ! test for Trace(mat) = ndim
        sumdiag = 0
        sumoff = 0
        do i = 1, ndim
            sumdiag = sumdiag + DBLE(mat(i,i))
        enddo
        do i = 1, ndim-3
            sumoff = sumoff + DBLE(mat(i,i+2))
        enddo
        write (*,*) sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-15) then
        ERROR STOP 3
        endif
        if (abs(sumoff) > 1E-15) then !FIXME: this limit is a bit scale less...
        ERROR STOP 6
        endif

        call test%adjoint_over_two(mat)
        ! test for Trace(mat) = ndim
        sumdiag = 0
        sumoff = 0
        do i = 1, ndim
            sumdiag = sumdiag + DBLE(mat(i,i))
        enddo
        do i = 1, ndim-3
            sumoff = sumoff + DBLE(mat(i,i+2))
        enddo
        write (*,*) sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-15) then
        ERROR STOP 7
        endif
        if (abs(sumoff) > 1E-15) then !FIXME: this limit is a bit scale less...
        ERROR STOP 14
        endif

        write (*,*) "success"
end Program HomogeneousExpTest
