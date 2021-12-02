! compile with
! gfortran -I ../../Libraries/libmscbdecomp/ -L ../../Libraries/libmscbdecomp/ 5-EulerExponential.F90 -lmscbdecomp

subroutine exectest(gd, ndim, mys)
  Use Exponentials_mod
  Use colorvertex_mod
  Use graphdata_mod
  implicit none
        type(EulerExp) :: ee
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
        
        ee = createEulerExponentialfromGraphData(gd, mys)
        
        call ee%lmult(mat)
        call ee%lmult(mat)
        call ee%lmult(mat)
        call ee%lmultinv(mat)
        call ee%lmultinv(mat)
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
        

        call ee%rmult(mat)
        call ee%rmult(mat)
        call ee%rmult(mat)
        call ee%rmultinv(mat)
        call ee%rmultinv(mat)
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
        
        call ee%lmult_T(mat)
        call ee%lmult_T(mat)
        call ee%lmult_T(mat)
        call ee%lmultinv_T(mat)
        call ee%lmultinv_T(mat)
        call ee%lmultinv_T(mat)

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
        ERROR STOP 102
        endif
        if (abs(sumoff) > maxval(exp(mys))*1E-14) then 
        ERROR STOP 104
        endif
        

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo
        

        call ee%rmult_T(mat)
        call ee%rmult_T(mat)
!         call ee%rmult_T(mat)
!         call ee%rmult_T(mat)
!         call ee%rmult_T(mat)
!         call ee%rmult_T(mat)
!         call ee%rmultinv_T(mat)
!         call ee%rmultinv_T(mat)
!         call ee%rmultinv_T(mat)
!         call ee%rmultinv_T(mat)
        call ee%rmultinv_T(mat)
        call ee%rmultinv_T(mat)
        ! test for Trace(mat) = ndim
        sumdiag = 0
        sumoff = 0
        do i = 1, ndim
            sumdiag = sumdiag + DBLE(mat(i,i))
        enddo
        do i = 1, ndim-3
            sumoff = sumoff + DBLE(mat(i,i+2))
        enddo
        write (*,*) "rmult_T", sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-14) then
        ERROR STOP 103
        endif
        if (abs(sumoff) > maxval(exp(mys))* 1E-14) then
        ERROR STOP 106
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
        if (abs(sumdiag - ndim) > ndim*1E-14) then
        ERROR STOP 7
        endif
        if (abs(sumoff) > maxval(exp(mys))*1E-14) then
        ERROR STOP 14
        endif
        
        call ee%dealloc()
        
end subroutine 

Program EulerExpTest

  Use Exponentials_mod
  Use colorvertex_mod
  Use graphdata_mod
  Use MvG_mod
  
  interface
  subroutine exectest(gd, ndim, mys)
  Use Exponentials_mod
  Use graphdata_mod
        integer :: ndim
        Type(GraphData) :: gd
        real(kind=kind(0.D0)), allocatable :: mys(:)
    end subroutine
  end interface
        integer, parameter :: ndim = 27
        type(GraphData) :: gd
        real(kind=kind(0.D0)) :: weight, sumdiag, sumoff
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
        input = input * 0.1
        
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
        call exectest(gd, ndim, mys)
        
        ! Now test homogeneous exponentials
        mys = 0.5
       call exectest(gd, ndim, mys)

        ! Now test general exponentials
      do i = 1, ndim
        mys(i) = 0.1*i
      enddo
      call exectest(gd, ndim, mys)

        call dealloc_graphdata(gd)
        deallocate(input, mys)
        write (*,*) "success"
end Program EulerExpTest
