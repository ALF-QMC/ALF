! compile with
! gfortran -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program CT11

  Use MyMats

        COMPLEX (KIND=8), DIMENSION(3,3) :: A, B
        COMPLEX (KIND=8) :: myx
        
        DO i = 1,3
	    DO j=1,3
		A(i,j) = CMPLX(DBLE(I), DBLE(j))
	    ENDDO
	ENDDO
	
        B = CT(A)
        
        DO i = 1,3
	    DO j=1,3
		if (B(j,i) .ne. CONJG(A(i,j))) then
		    STOP 2
		ENDIF
	    ENDDO
	ENDDO
end Program CT11
