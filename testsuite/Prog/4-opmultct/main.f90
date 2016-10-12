! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program OPMULTTEST

Use Operator_mod


        COMPLEX (KIND=KIND(0.D0)), DIMENSION(3,3) :: U
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(3,5) :: V
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(5,5) :: matnew, matold
        Complex (KIND = KIND(0.D0)) :: tmp
        Integer, DIMENSION(3) :: P
        Integer :: i,j,n, opn, Ndim
        
        
        opn = 3
        Ndim = 5
        do i=1,Ndim
            do j=1,Ndim
                matnew(i,j) = CMPLX(i, j, kind(0.D0))
                matold(i,j) = CMPLX(i, j, kind(0.D0))
            enddo
        enddo
        
        do i = 1, opn
        do j = 1, Ndim
        V(i,j) = CMPLX(i, j, kind(0.D0))
        enddo
        P(i) = i
        enddo
        
        do i = 1, opn
        do j = 1, opn
        U(i,j) = CMPLX(i, j, kind(0.D0))
        enddo
        enddo
        
        call opmultct(V, U, P, matnew, opn, Ndim)
        
! check against old version
    do n = 1, opn
        DO I = 1, Ndim
           tmp=cmplx(0.d0,0.d0, kind(0.D0))
       Do m = 1, opn
          tmp = tmp + conjg(U(n,m))*V(m,I)
       Enddo
       Matold(I, P(n))  = tmp 
         enddo
       enddo
    
    do i = 1,Ndim
        do j = 1,Ndim
        if (matold(i, j) .ne. matnew(i, j)) then
        STOP 2
        endif
        enddo
    enddo

end Program OPMULTTEST
