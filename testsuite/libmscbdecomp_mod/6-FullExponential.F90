! compile with
! gfortran -I ../../Libraries/libmscbdecomp/ -L ../../Libraries/libmscbdecomp/ 1-ZeroDiag-lmult.F90 -lmscbdecomp

Program FullExpTest

  Use Exponentials_mod

        COMPLEX (KIND=8) :: myx
        
        type(FullExp) :: fe
        integer, parameter :: nredges = 6
        integer, parameter :: ndim = 24
        integer, parameter :: usedcolors = 1
        integer, parameter :: method = 3
        Type(Node) :: nodes(nredges)
        real(kind=kind(0.D0)) :: weight, sumdiag, sumoff
        real(kind=kind(0.D0)), allocatable :: mys(:)
        COMPLEX (KIND=kind(0.D0)), DIMENSION(ndim, ndim) :: mat
        integer :: i
        
        allocate(mys(ndim))
        nodes(1)%x=1
        nodes(1)%y=3
        nodes(1)%col=1
        nodes(1)%axy=0.1

        nodes(2)%x=5
        nodes(2)%y=7
        nodes(2)%col=1
        nodes(2)%axy=0.1
        
        nodes(3)%x=9
        nodes(3)%y=11
        nodes(3)%col=1
        nodes(3)%axy=0.1
        
        nodes(4)%x=13
        nodes(4)%y=15
        nodes(4)%col=1
        nodes(4)%axy=0.1

        nodes(5)%x=17
        nodes(5)%y=19
        nodes(5)%col=1
        nodes(5)%axy=0.1

        nodes(6)%x=21
        nodes(6)%y=23
        nodes(6)%col=1
        nodes(6)%axy=0.1

        weight = 1.0
        mys = 0.0 ! initialize chemical potential to zero

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo

        call fe%init(nodes, usedcolors, mys, method, weight)
        
        call fe%lmult(mat)
        call fe%lmultinv(mat)

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
        


        call fe%rmult(mat)
        call fe%rmultinv(mat)
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

        call fe%adjoint_over_two(mat)
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
end Program FullExpTest
