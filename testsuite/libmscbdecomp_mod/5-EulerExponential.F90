! compile with
! gfortran -I ../../Libraries/libmscbdecomp/ -L ../../Libraries/libmscbdecomp/ 1-ZeroDiag-lmult.F90 -lmscbdecomp

subroutine exectest(nodes, nredges, ndim, mys)
  Use Exponentials_mod
  implicit none
        type(EulerExp) :: ee
        integer :: nredges
        integer :: ndim
        integer, parameter :: usedcolors = 1
        Type(Node) :: nodes(nredges)
        real(kind=kind(0.D0)) :: weight, sumdiag, sumoff
        real(kind=kind(0.D0)), allocatable :: mys(:)
        COMPLEX (KIND=kind(0.D0)), DIMENSION(ndim, ndim) :: mat
        integer :: i

        weight = 1.0

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo

        call ee%init(nodes, usedcolors, mys, weight)
        
        call ee%lmult(mat)
        call ee%lmultinv(mat)

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
        


        call ee%rmult(mat)
        call ee%rmultinv(mat)
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

        call ee%adjoint_over_two(mat)
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
        
        call ee%dealloc()
        
end subroutine 

Program EulerExpTest

  Use Exponentials_mod
  
  interface
  subroutine exectest(nodes, nredges, ndim, mys)
  Use Exponentials_mod
        integer :: nredges
        integer :: ndim
        Type(Node) :: nodes(nredges)
        real(kind=kind(0.D0)), allocatable :: mys(:)
    end subroutine
  end interface
        type(EulerExp) :: ee
        integer, parameter :: nredges = 6
        integer, parameter :: ndim = 24
        integer, parameter :: usedcolors = 1
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

        ! Start with zero diagonals
        mys = 0.0 ! initialize chemical potential to zero
        call exectest(nodes, nredges, ndim, mys)
        
        ! Now test homogeneous exponentials
        mys(1) = 0.5
        mys(3) = 0.5
        mys(5) = 0.5
        mys(7) = 0.5
        mys(9) = 0.5
        mys(11) = 0.5
        mys(13) = 0.5
        mys(15) = 0.5
        mys(17) = 0.5
        mys(19) = 0.5
        mys(21) = 0.5
        mys(23) = 0.5
        call exectest(nodes, nredges, ndim, mys)
        
        ! Now test traceless exponentials
        mys(1) = 0.5
        mys(3) = -0.5
        mys(5) = 0.5
        mys(7) = -0.5
        mys(9) = 0.5
        mys(11) = -0.5
        mys(13) = 0.5
        mys(15) = -0.5
        mys(17) = 0.5
        mys(19) = -0.5
        mys(21) = 0.5
        mys(23) = -0.5
        call exectest(nodes, nredges, ndim, mys)        
        
        ! Now test general exponentials
        mys(1) = 0.5
        mys(3) = -0.6
        mys(5) = 0.5
        mys(7) = -0.7
        mys(9) = 0.5
        mys(11) = -0.9
        mys(13) = 0.5
        mys(15) = -0.5
        mys(17) = 0.2
        mys(19) = -0.5
        mys(21) = 0.1
        mys(23) = -0.5
        call exectest(nodes, nredges, ndim, mys)        
        write (*,*) "success"
end Program EulerExpTest
