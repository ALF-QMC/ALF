! compile with
!  gfortran -std=f2003  -I ../../../Prog_8/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ test2.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a


Program OPEXPMULTTEST

Use Operator_mod
implicit none

        COMPLEX (KIND = KIND(0.D0)), DIMENSION(3,3) :: U
        COMPLEX (KIND = KIND(0.D0)), DIMENSION(3,5) :: V
        COMPLEX (KIND = KIND(0.D0)), DIMENSION(3) :: Z
        COMPLEX (KIND = KIND(0.D0)), DIMENSION(5,5) :: matnew, matold
        Complex (KIND = KIND(0.D0)) :: tmp, Z1,g
        Real (Kind = KIND(0.D0)), DIMENSION(3) :: E
        Real :: spin
        Integer, DIMENSION(3) :: P
        Integer :: i, j, n, m, opn, Ndim
        
        
        opn = 3
        Ndim = 5
        spin = -1.D0
        g = 2.D0
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
        E(i) = 3.D5*i
        enddo
        
        do i = 1, opn
        do j = 1, opn
        U(i,j) = CMPLX(i, j, kind(0.D0))
        enddo
        enddo
        
        Z = exp(g*spin*E)
        call opexpmult(V, U, P, matnew, Z, opn, Ndim)
        
! check against old version
    do n = 1,opn
       Z1 = exp(g*cmplx(E(n)*spin,0.d0, kind(0.D0)))
       Do I = 1,Ndim
          tmp = cmplx(0.d0,0.d0, kind(0.D0))
          DO m = 1,opn
             tmp  = tmp + V(m,I) * U(m,n)
          Enddo
          Matold(I,P(n))  = tmp * Z1
       enddo
    Enddo
    
    do i = 1,Ndim
        do j = 1,Ndim
        if (matold(i, j) .ne. matnew(i, j)) then
        STOP 2
        endif
        enddo
    enddo

end Program OPEXPMULTTEST
