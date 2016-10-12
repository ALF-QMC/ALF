! compile with
! gfortran -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ main.f90 ../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyLin/liblin.a ../../Libraries/MyNag/libnag.a

Program Test2

  Use MyMats 

        REAL (KIND=8), DIMENSION(3,3) :: A
        REAL (KIND=8), DIMENSION(3,3) :: AI
        REAL (KIND=8), DIMENSION(3,3) :: B
        REAL (KIND=8), DIMENSION(3) :: W
        REAL (KIND=8), DIMENSION(2) :: myDET
        
        A= 0d0;
        A(1,1) = 10
        A(2,1) = 1
        A(1,2) = 1
        A(2,2) = 20
        A(3,3) = 100
        call INV_R1(A, AI, myDET)
! Yes 1E-11 is really the precission that is achievable here using the linpack routines
        if (ABS(mydet(1)*(10.0**myDET(2)) - 19900) > 19900*10*EPSILON(myDet(1))) then
        write (*,*) ABS(mydet(1)*(10.0**myDET(2)) - 19900)
        STOP 2
        endif
        B = MATMUL(A, AI)
        if (ABS(MAXVAL(B) - 1) > 1E-11 ) then
        write (*,*) B
        STOP 3
        endif
        if (ABS(MINVAL(B)) > 1E-11 ) then
        write (*,*) B
        STOP 4
        endif
end Program Test2
