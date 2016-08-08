! compile with
! gfortran -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ main.f90 ../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyLin/liblin.a ../../Libraries/MyNag/libnag.a

Program Main

  Use MyMats 

        REAL (KIND=8), DIMENSION(3,3) :: A
        REAL (KIND=8), DIMENSION(3,3) :: U
        REAL (KIND=8), DIMENSION(3,3) :: D
        REAL (KIND=8), DIMENSION(3,3) :: B
        REAL (KIND=8), DIMENSION(3) :: W
        
        A= 0d0;
        A(1,1) = 10
        A(2,1) = 1
        A(1,2) = 1
        A(2,2) = 20
        A(3,3) = 100
        D=0
        
        call DIAG_R(A, U, W)
        do i=1,3
        D(i,i) = W(i)
        enddo

        B = MATMUL(U, MATMUL(D,TRANSPOSE(U)))
        D=ABS(B-A)
        if (MAXVAL(D) > 1E-11) then
        stop 2
        endif
end Program Main
