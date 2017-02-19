! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a
!
Program FillExpOpsTest
!
      Use Operator_mod
      Implicit None
!
      Complex (Kind=Kind(0.D0)), Dimension (3) :: ExpOp, ExpMOp, &
     & ExpOpold, ExpMOpold
      Real (Kind=Kind(0.D0)) :: spin, diffre, diffim
      Integer :: i, n
      Type (Operator) :: Op
!
      Call Op_make (Op, 3)
!
      Do i = 1, Op%n
         Op%E (i) = 2 * i - 3
      End Do
      Op%N_non_zero = 2
      Op%g = 2
      spin = - 1.0
!
      Call FillExpOps (ExpOp, ExpMOp, Op, spin)
!
! check against old version
      Do n = 1, Op%n
         ExpOpold (n) = cmplx (1.d0, 0.d0, kind(0.D0))
         If (n <= Op%N_non_zero) ExpOpold (n) = Exp &
        & (Op%g*cmplx(Op%E(n)*spin, 0.d0, kind(0.D0)))
         ExpMOpold (n) = cmplx (1.d0, 0.d0, kind(0.D0))
         If (n <= Op%N_non_zero) ExpMOpold (n) = Exp &
        & (-Op%g*cmplx(Op%E(n)*spin, 0.d0, kind(0.D0)))
      End Do
!
      Do i = 1, 3
         If (ExpOpold(i) .Ne. ExpOp(i)) Then
            Write (*,*) "ERROR", ExpOpold (i), ExpOp (i)
            Stop 2
         End If
         diffre = Abs (real(ExpMOpold(i)-ExpMOp(i)))
         diffim = Abs (aimag(ExpMOpold(i)-ExpMOp(i)))
         If (diffre > Max(Abs(real(ExpMOpold(i))), &
        & Abs(real(ExpMOp(i))))*1D-15) Then
            Write (*,*) "ERROR", ExpMOpold (i), ExpMOp (i)
            Stop 3
         End If
         If (diffim > Max(Abs(aimag(ExpMOpold(i))), &
        & Abs(aimag(ExpMOp(i))))*1D-15) Then
            Write (*,*) "ERROR", ExpMOpold (i), ExpMOp (i)
            Stop 4
         End If
      End Do
!
End Program FillExpOpsTest
