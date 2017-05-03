! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program OPMULTTEST

Use Operator_mod
implicit none

        Complex (Kind=Kind(0.D0)) :: Matnew(3,3), matold(3,3), VH(3,3), Z, Z1, Zre, Zim
        Complex (Kind=Kind(0.D0)), allocatable, dimension(:,:) :: Uold
        Complex (Kind=Kind(0.D0)), allocatable, dimension(:) ::TAU, R, work
        Real (KIND = KIND(0.D0)) :: spin
        Integer :: i, n, m, j, ndim, lwork, info, opn
        Type(Operator) :: Op
        
! setup some test data
        Ndim = 3
        opn = 3
        lwork = 2*opn
        call op_seths()
        Call Op_make(Op, opn)
        Allocate (Uold(opn, opn), Tau(opn), work(lwork), R(opn))
        do i = 1, Op%N
            Op%P(i) = i
            do n = 1,Op%N
            Op%O(i, n) = CMPLX(n+i, n-i, kind(0.D0))
            enddo
        enddo
        CALL Op_Set(Op)
        Op%N_non_zero = 2
        Op%g = 2.D0
        spin =-1.0
        
        do i = 1,Ndim
        do n = 1,Ndim
        matnew(i,n) = CMPLX(i,n, kind(0.D0))
        matold(i,n) = CMPLX(i,n, kind(0.D0))
        enddo
        enddo
        Uold = Op%U
        TAU = Op%U(:, opn)
        Do i = 1, opn - 1
            R(i) = Op%U(i,i)
        ENDDO
        CALL ZUNGQR(Op%N, Op%N, Op%N, Uold, Op%N, TAU, WORK, LWORK, INFO)
        DO i = 1, opn
        Do j = 1, opn
            Uold(i, j) = Uold(i,j) * R(j)
        ENDDO
        ENDDO
        Call Op_mmultL(matnew, Op, spin, Ndim)

! check against old version from Operator_FFA.f90

VH = 0.D0
    do n = 1,Op%N
       Z = exp(Op%g*Op%E(n)*spin)! fixed in contrast to Fakhers old version
       Do m = 1,Op%N
          Z1 = Uold(m,n)* Z
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
          Z1 = conjg(Uold(n,m))
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
    if (Abs(Zre) > MAX(ABS(real(matnew(i,j))), ABS(real(matold(i,j))) )*1D-13) then
    write (*,*) "ERROR in real part", real(matnew(i,j)), real(matold(i,j))
    STOP 2
    endif
    if (Abs(Zim) > MAX(ABS(aimag(matnew(i,j))), ABS(aimag(matold(i,j))) )*1D-13) then
    write (*,*) "ERROR in imag part", aimag(matnew(i,j)), aimag(matold(i,j))
    STOP 3
    endif
    enddo
    enddo
    Deallocate(R, TAU, WORK, Uold)
write (*,*) "success"
end Program OPMULTTEST
