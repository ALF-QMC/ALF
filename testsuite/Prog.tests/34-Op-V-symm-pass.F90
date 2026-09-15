! Test that symmetric operator sequences pass the Op_V_is_symmetric check.
! Test 1: palindromic arrangement of non-commuting operators.
! Test 2: commuting on-site operators in non-palindromic order (Hubbard-like).
! Should exit successfully (return code 0).
!
Program Test_Op_V_symm_pass
  Use Operator_mod
  Implicit None

  Type (Operator), allocatable :: Op_V(:,:)
  Integer :: I, Ndim

  Ndim = 4

  ! --- Test 1: palindromic sequence (non-commuting operators, but palindromic order) ---
  allocate(Op_V(5, 1))
  do I = 1, 5
     Call Op_make(Op_V(I, 1), 2)
     Op_V(I, 1)%P(1) = 1
     Op_V(I, 1)%P(2) = 2
     Op_V(I, 1)%type  = 2
     Op_V(I, 1)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  enddo

  ! Op(1) = Op(5): diagonal
  Op_V(1, 1)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.d0))
  Op_V(1, 1)%O(2,2) = cmplx(-1.d0, 0.d0, kind(0.d0))
  Op_V(1, 1)%g      = cmplx(0.5d0, 0.d0, kind(0.d0))
  Op_V(5, 1)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.d0))
  Op_V(5, 1)%O(2,2) = cmplx(-1.d0, 0.d0, kind(0.d0))
  Op_V(5, 1)%g      = cmplx(0.5d0, 0.d0, kind(0.d0))

  ! Op(2) = Op(4): off-diagonal
  Op_V(2, 1)%O(1,2) = cmplx(0.d0, 1.d0, kind(0.d0))
  Op_V(2, 1)%O(2,1) = cmplx(0.d0, -1.d0, kind(0.d0))
  Op_V(2, 1)%g      = cmplx(0.3d0, 0.d0, kind(0.d0))
  Op_V(4, 1)%O(1,2) = cmplx(0.d0, 1.d0, kind(0.d0))
  Op_V(4, 1)%O(2,1) = cmplx(0.d0, -1.d0, kind(0.d0))
  Op_V(4, 1)%g      = cmplx(0.3d0, 0.d0, kind(0.d0))

  ! Op(3): middle
  Op_V(3, 1)%O(1,1) = cmplx(2.d0, 0.d0, kind(0.d0))
  Op_V(3, 1)%O(2,2) = cmplx(2.d0, 0.d0, kind(0.d0))
  Op_V(3, 1)%g      = cmplx(1.d0, 0.d0, kind(0.d0))

  if (.not. Op_V_is_symmetric(Op_V, 1, Ndim)) then
     write(*,*) 'FAIL: palindromic sequence detected as asymmetric'
     STOP 2
  endif
  write(*,*) 'PASS: palindromic sequence'
  deallocate(Op_V)

  ! --- Test 2: commuting, non-palindromic (Hubbard-like on-site operators) ---
  ! Four 1x1 operators on different sites in non-palindromic order [1, 3, 2, 4]
  allocate(Op_V(4, 1))
  do I = 1, 4
     Call Op_make(Op_V(I, 1), 1)
     Op_V(I, 1)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.d0))
     Op_V(I, 1)%alpha   = cmplx(0.d0, 0.d0, kind(0.d0))
     Op_V(I, 1)%type    = 1
  enddo
  Op_V(1, 1)%P(1) = 1
  Op_V(1, 1)%g    = cmplx(0.5d0, 0.d0, kind(0.d0))
  Op_V(2, 1)%P(1) = 3
  Op_V(2, 1)%g    = cmplx(0.7d0, 0.d0, kind(0.d0))
  Op_V(3, 1)%P(1) = 2
  Op_V(3, 1)%g    = cmplx(0.3d0, 0.d0, kind(0.d0))
  Op_V(4, 1)%P(1) = 4
  Op_V(4, 1)%g    = cmplx(0.9d0, 0.d0, kind(0.d0))

  if (.not. Op_V_is_symmetric(Op_V, 1, Ndim)) then
     write(*,*) 'FAIL: commuting non-palindromic detected as asymmetric'
     STOP 2
  endif
  write(*,*) 'PASS: commuting non-palindromic (Hubbard-like)'
  deallocate(Op_V)

  write(*,*) 'success'

end Program Test_Op_V_symm_pass
