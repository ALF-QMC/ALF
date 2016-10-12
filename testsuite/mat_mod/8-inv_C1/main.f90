! compile with
! gfortran -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ main.f90 ../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyLin/liblin.a ../../Libraries/MyNag/libnag.a

Program Test8

  Use MyMats 

        COMPLEX (KIND=8), DIMENSION(3, 3) :: A
        COMPLEX (KIND=8), DIMENSION(3, 3) :: AI
        COMPLEX (KIND=8), DIMENSION(3, 3) :: B
        COMPLEX (KIND=8), DIMENSION(3) :: W
        COMPLEX (KIND=8), DIMENSION(2) :: myDET
        
        A = (0.d0, 0.d0);
        A(1, 1) = 10
        A(2, 1) = (1,0)
        A(1, 2) = (1,0)
        A(2, 2) = 20
        A(3, 3) = 100
        call INV_C1(A, AI, myDET)
! Yes 1E-11 is really the precission that is achievable here using the linpack routines
        if (ABS(ABS(mydet(1)*(10.0**myDET(2))) - 19900) > 19900*10*EPSILON(1.d0)) then
        write (*,*) ABS(mydet(1)*(10.0**myDET(2))) - 20001
        STOP 2
        endif
        B = MATMUL(A, AI)
        if (ABS(B(2,2) - 1.0) > 10*EPSILON(1.d0) ) then
        write (*,*) B
        STOP 3
        endif
        if (ABS(B(2,1)) > EPSILON(1.d0) ) then
        write (*,*) B
        STOP 4
        endif
end Program Test8
