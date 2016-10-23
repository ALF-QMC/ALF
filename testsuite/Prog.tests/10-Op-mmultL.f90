! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program OPMULTTEST

Use Operator_mod
implicit none

        Complex (Kind=Kind(0.D0)) :: Matnew(3,3), matold(3,3), VH(3,3), Z, Z1, Zre, Zim
        Real (KIND = KIND(0.D0)) :: spin
        Integer :: i, n, m, j, ndim
        Type(Operator) :: Op
        
! setup some test data
        Ndim = 3
        call op_seths()
        Call Op_make(Op, 3)
        
        do i = 1, Op%N
            Op%E(i) = 2*i-3
            Op%P(i) = i
            do n = 1,Op%N
            Op%U(i,n) = CMPLX(n, i, kind(0.D0))
            enddo
        enddo
        Op%N_non_zero = 2
        Op%g = 2.D0
        spin =-1.0
        
        do i = 1,Ndim
        do n = 1,Ndim
        matnew(i,n) = CMPLX(i,n, kind(0.D0))
        matold(i,n) = CMPLX(i,n, kind(0.D0))
        enddo
        enddo
        
        Call Op_mmultL(matnew, Op, spin, Ndim)

! check against old version from Operator_FFA.f90

VH = 0.D0
    do n = 1,Op%N
       Z = exp(Op%g*Op%E(n)*spin)! fixed in contrast to Fakhers old version
       Do m = 1,Op%N
          Z1 = Op%U(m,n)* Z
          DO I = 1,Ndim
             VH(I,n)  = VH(I,n) + Matold(I,Op%P(m)) * Z1
          Enddo
       enddo
    Enddo
    Do n = 1,Op%N
       Do I = 1,Ndim
          Matold(I,Op%P(n)) =   VH(I,n) 
       Enddo
    Enddo


    VH = 0.D0
    do n = 1,Op%N
       Do m = 1,Op%N
          Z1 = conjg(Op%U(n,m))
          DO I = 1,Ndim
             VH(I,n)  = VH(I,n) + Matold(I,Op%P(m)) * Z1
          Enddo
       enddo
    Enddo
    Do n = 1,Op%N
       Do I = 1,Ndim
          Matold(I,Op%P(n)) =   VH(I,n) 
       Enddo
    Enddo

    do i=1,3
    do j=1,3
    Zre = real(matnew(i,j)-matold(i,j))
    Zim = aimag(matnew(i,j)-matold(i,j))
    if (Abs(Zre) > MAX(ABS(real(matnew(i,j))), ABS(real(matold(i,j))) )*1D-15) then
    write (*,*) "ERROR in real part", real(matnew(i,j)), real(matold(i,j))
    STOP 2
    endif
    if (Abs(Zim) > MAX(ABS(aimag(matnew(i,j))), ABS(aimag(matold(i,j))) )*1D-15) then
    write (*,*) "ERROR in imag part", aimag(matnew(i,j)), aimag(matold(i,j))
    STOP 3
    endif
    enddo
    enddo
write (*,*) "success"
end Program OPMULTTEST
