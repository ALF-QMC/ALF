! compile with
! gfortran -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ main.f90 ../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyLin/liblin.a ../../Libraries/MyNag/libnag.a

Program Test10

  Use MyMats

        COMPLEX (KIND=8), DIMENSION(3,3) :: A
        COMPLEX (KIND=8) :: myx
        
        myx = (1.d0,1.d0)
        CALL INITD_C(A,myx)
        IF(ABS(ABS(A(1,1)) - SQRT(2.D0)) > 10*EPSILON(1.D0)) THEN
        write (*,*) "ERROR", A
        STOP 2
        ENDIF
        IF(ABS(A(2,1)) > 10*EPSILON(1.D0)) THEN
        write (*,*) "ERROR", A
        STOP 3
        ENDIF
        write (*,*) "success"
end Program Test10
