! compile with
! gfortran -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ main.f90 ../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyLin/liblin.a ../../Libraries/MyNag/libnag.a

Program Test4

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
        call INV_R_VARIABLE_1(A, AI, myDET, 2)
! Yes 1E-11 is really the precission that is achievable here using the linpack routines
        if (ABS(myDET(1)*10.0**myDET(2) - 199) > 1E-11) then
        write (*,*) ABS(myDET(1)*10.0**myDET(2) - 199)
        STOP 2
        endif
        B = MATMUL(A, AI)
        !test negative determinants
        A(1,1) = 10
        A(2,1) = 1
        A(1,2) = 1
        A(2,2) = -20
        call INV_R_VARIABLE_1(A, AI, myDET, 2)
! Yes 1E-11 is really the precission that is achievable here using the linpack routines
        if (ABS(myDET(1)*10.0**myDET(2) + 201) > 1E-11) then
        write (*,*) "ERROR", ABS(myDET(1)*10.0**myDET(2) + 201)
        STOP 2
        endif

!        if (ABS(MAXVAL(B) - 1) > 1E-11 ) then
!        write (*,*) B
!        STOP 3
!        endif
!        if (ABS(MINVAL(B)) > 1E-11 ) then
!        write (*,*) B
!        STOP 4
!        endif
write (*,*) "success"
end Program Test4
