! compile with
! gfortran -Wall -std=f2003  -I ../../../Prog_8/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ test2.f90 ../../../Prog_8/cgr2_1.o ../../../Libraries/Modules/modules_90.a ../../../Prog_8/UDV_WRAP.o ../../../Prog_8/cgr1.o -llapack -lblas ../../../Libraries/MyNag/libnag.a


Program TESTSCALEMATRIX

        COMPLEX (KIND=KIND(0.D0)), DIMENSION(5) :: Z
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(5,5) :: matnew, matold
        Integer :: i, j, Ndim
        COMPLEX (KIND=KIND(0.D0)) :: tmp

        Ndim = 5
        do i=1,Ndim
            do j=1,Ndim
                matnew(i, j) = CMPLX(i, j, kind(0.D0))
                matold(i, j) = CMPLX(i, j, kind(0.D0))
            enddo
            Z(i) = i
        enddo
        
        call scalematrix(matnew, Z, .TRUE., 5)
        
! check against old version
           DO J = 1,5
              DO I = 1,5
                 matold(I,J) = matold(I,J)/CONJG(Z(J))
              ENDDO
           ENDDO
    
    do i = 1,Ndim
        do j = 1,Ndim
        tmp = matold(i,j) - matnew(i,j)
        if (Aimag(tmp) > Abs(Aimag(matnew(i,j)))*1.D-15  ) then
        write (*,*) "ERROR", matold(i,j), matnew(i,j)
        STOP 2
        endif
        if (Real(tmp) > Abs(Real(matnew(i,j)))*1.D-15  ) then
        write (*,*) "ERROR", matold(i,j), matnew(i,j)
        STOP 2
        endif

        enddo
    enddo
write (*,*) "success"
end Program TESTSCALEMATRIX
