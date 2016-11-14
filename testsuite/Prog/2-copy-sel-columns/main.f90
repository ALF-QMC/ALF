! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -I ../../../Prog_8/ -L ../../../Libraries/Modules/ -L ../../../Prog_8/  main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a 

Program CSC

Use Operator_mod


        COMPLEX (KIND=KIND(0.D0)), DIMENSION(3,3) :: Vnew, Vold
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(5,5) :: mat
        Integer, DIMENSION(3) :: P
        Integer :: i,j,n
        
        do i=1,3
            do j=1,5
                mat(i,j) = CMPLX(i, j, kind(0.D0))
            enddo
        P(i) = i
        enddo
        
        call copy_select_columns(Vnew, mat, P, 3, 5)
        
! check old version
    Do n = 1,3
       Do I = 1,3
          Vold(n,I) = Mat(P(n),I)
       Enddo
    Enddo

    do i = 1,3
        do j = 1,3
        if (Vold(i,j) .ne. Vnew(i,j)) then
        STOP 2
        endif
        enddo
    enddo

end Program CSC
