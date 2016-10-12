! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program OPMULTTEST

Use Operator_mod
implicit none

        COMPLEX (KIND=KIND(0.D0)), DIMENSION(3) :: ExpOp, ExpMOp, ExpOpold, ExpMOpold
        Real (KIND = KIND(0.D0)) :: spin
        Integer :: i,n
        Type(Operator) :: Op
        
        Call Op_make(Op, 3)
        
        do i = 1, Op%N
            Op%E(i) = 2*i-3
        enddo
        Op%N_non_zero = 2
        Op%g = 2
        spin =-1.0
        
        Call FillExpOps(ExpOp, ExpMOp, Op, spin)

! check against old version
       do n = 1, Op%N
          ExpOpold(n) = cmplx(1.d0, 0.d0, kind(0.D0))
          If ( n <= OP%N_non_Zero) ExpOpold(n) = exp(Op%g*cmplx(Op%E(n)*spin,0.d0, kind(0.D0)))
          ExpMOpold(n) = cmplx(1.d0, 0.d0, kind(0.D0))
          If ( n <= OP%N_non_Zero) ExpMOpold(n) = exp(-Op%g*cmplx(Op%E(n)*spin,0.d0, kind(0.D0)))
       enddo

do i = 1, 3
if (ExpOpold(i) .ne. ExpOp(i)) then
write (*,*) Expopold(i), expop(i)
STOP 2
endif
if (ExpMOpold(i) .ne. ExpMOp(i)) then
STOP 3
endif
enddo


end Program OPMULTTEST
