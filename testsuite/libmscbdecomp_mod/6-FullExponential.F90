! compile with
! gfortran -I ../../Libraries/libmscbdecomp/ -L ../../Libraries/libmscbdecomp/ 6-FullExponential.F90 -lmscbdecomp

subroutine exectest(gd, ndim, mys, method)
  Use Exponentials_mod
  Use colorvertex_mod
  Use graphdata_mod
  implicit none
        integer :: method
        type(FullExp) :: fe
        integer :: ndim
        integer, parameter :: usedcolors = 1
        Type(GraphData) :: gd
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
        
        fe = createFullExponentialfromGraphData(gd, mys, method)
        
        call fe%lmult(mat)
        call fe%lmult(mat)
        call fe%lmult(mat)
        call fe%lmultinv(mat)
        call fe%lmultinv(mat)
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
        if (abs(sumdiag - ndim) > ndim*1E-14) then
            ERROR STOP 2
        endif
        if (abs(sumoff) >  maxval(exp(mys))*1E-14) then
            ERROR STOP 4
        endif
        

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo
        

        call fe%rmult(mat)
        call fe%rmult(mat)
        call fe%rmult(mat)
        call fe%rmultinv(mat)
        call fe%rmultinv(mat)
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
        write(*,*) "rmult  ", sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-14) then
        ERROR STOP 3
        endif
        if (abs(sumoff) > maxval(exp(mys))*1E-14) then
        ERROR STOP 6
        endif

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo
        
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
        if (abs(sumdiag - ndim) > ndim*1E-14) then
        ERROR STOP 7
        endif
        if (abs(sumoff) > maxval(exp(mys))*1E-14) then
        ERROR STOP 14
        endif
        
        call fe%dealloc()
        
end subroutine 

Program FullExpTest

  Use Exponentials_mod
  Use colorvertex_mod
  Use MvG_mod
  Use graphdata_mod
  
  interface
    subroutine exectest(gd, ndim, mys, method)
        Use Exponentials_mod
        Use colorvertex_mod
        Use graphdata_mod
        integer :: ndim, method
        Type(GraphData) :: gd
        real(kind=kind(0.D0)), allocatable :: mys(:)
    end subroutine
  end interface
        integer, parameter :: ndim = 23
        type(GraphData) :: gd
        real(kind=kind(0.D0)) :: dt, sumdiag, sumoff
        real(kind=kind(0.D0)), allocatable :: mys(:)
        COMPLEX (KIND=kind(0.D0)), DIMENSION(:,:), allocatable :: input
        integer :: i, j
        
        allocate(mys(ndim), input(ndim, ndim) )
        input = 0
        do i = 1, ndim-1
             input(i, i + 1) = 1
             input(i + 1, i) = 1
        enddo
        
        do i = 1, ndim-2
             input(i, i + 2) = 1
             input(i + 2, i) = 1
        enddo
        dt = 0.1
        input = input * dt
        
        gd = mat2verts(input) ! convert to graphdata structure
        call MvG_decomp(gd%verts) ! perform the decomposition
        
        call determine_used_colors_of_graph(gd)
        write (*,*) "Nr edges: ", gd%nredges
        if (gd%usedcolors == gd%deltag) then
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families -> optimal decomposition"
        else
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families"
        endif

        ! Start with zero diagonals
        mys = 0.0 ! initialize chemical potential to zero
        do i = 2, 5
            call exectest(gd, ndim, mys, i)
        enddo
        ! Now test homogeneous exponentials
        mys = 0.5*dt
        do i = 2, 5
            call exectest(gd, ndim, mys, i)
        enddo

        ! Now test general exponentials
        do i = 1, ndim
            mys(i) = 0.05*i*dt
        enddo
        do i = 2, 5
            call exectest(gd, ndim, mys, i)
        enddo

        call dealloc_graphdata(gd)
        deallocate(input, mys)
        write (*,*) "success"
end Program FullExpTest
