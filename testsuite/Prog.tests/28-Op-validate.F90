! Test that Op_set accepts valid Hermitian operators of various kinds.
! This program should exit successfully (return code 0).
!
Program Test_Op_validate
  Use Operator_mod
  Implicit None

  Type (Operator) :: Op
  Integer :: I, J, N

  ! --- Test 1: 1x1 real operator ---
  Call Op_make(Op, 1)
  Op%O(1,1) = cmplx(3.d0, 0.d0, kind(0.d0))
  Op%P(1)   = 1
  Op%g      = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%alpha  = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type   = 1
  Call Op_set(Op)
  Call Op_clear(Op, 1)
  Write(*,*) 'PASS: 1x1 real operator'

  ! --- Test 2: 3x3 real symmetric operator ---
  N = 3
  Call Op_make(Op, N)
  Op%O(1,1) = cmplx( 2.d0, 0.d0, kind(0.d0))
  Op%O(1,2) = cmplx(-1.d0, 0.d0, kind(0.d0))
  Op%O(2,1) = cmplx(-1.d0, 0.d0, kind(0.d0))
  Op%O(2,2) = cmplx( 2.d0, 0.d0, kind(0.d0))
  Op%O(2,3) = cmplx(-1.d0, 0.d0, kind(0.d0))
  Op%O(3,2) = cmplx(-1.d0, 0.d0, kind(0.d0))
  Op%O(3,3) = cmplx( 2.d0, 0.d0, kind(0.d0))
  do I = 1, N
     Op%P(I) = I
  enddo
  Op%g     = cmplx(0.5d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = 1
  Call Op_set(Op)
  Call Op_clear(Op, N)
  Write(*,*) 'PASS: 3x3 real symmetric operator'

  ! --- Test 3: 3x3 complex Hermitian operator ---
  N = 3
  Call Op_make(Op, N)
  Op%O(1,1) = cmplx( 1.d0,  0.d0, kind(0.d0))
  Op%O(1,2) = cmplx( 2.d0,  1.d0, kind(0.d0))
  Op%O(2,1) = cmplx( 2.d0, -1.d0, kind(0.d0))  ! conjg of O(1,2)
  Op%O(1,3) = cmplx( 0.d0,  3.d0, kind(0.d0))
  Op%O(3,1) = cmplx( 0.d0, -3.d0, kind(0.d0))  ! conjg of O(1,3)
  Op%O(2,2) = cmplx( 4.d0,  0.d0, kind(0.d0))
  Op%O(2,3) = cmplx(-1.d0,  0.5d0, kind(0.d0))
  Op%O(3,2) = cmplx(-1.d0, -0.5d0, kind(0.d0)) ! conjg of O(2,3)
  Op%O(3,3) = cmplx( 2.d0,  0.d0, kind(0.d0))
  do I = 1, N
     Op%P(I) = I
  enddo
  Op%g     = cmplx(0.1d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = 2
  Call Op_set(Op)
  Call Op_clear(Op, N)
  Write(*,*) 'PASS: 3x3 complex Hermitian operator'

  ! --- Test 4: 4x4 diagonal operator (all types) ---
  N = 4
  Call Op_make(Op, N)
  do I = 1, N
     Op%O(I,I) = cmplx(real(I, kind(0.d0)), 0.d0, kind(0.d0))
     Op%P(I) = I
  enddo
  Op%g     = cmplx(0.3d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = 3
  Call Op_set(Op)
  Call Op_clear(Op, N)
  Write(*,*) 'PASS: 4x4 diagonal operator (type 3)'

  ! --- Test 5: type 4 operator ---
  N = 2
  Call Op_make(Op, N)
  Op%O(1,1) = cmplx( 1.d0, 0.d0, kind(0.d0))
  Op%O(1,2) = cmplx( 0.5d0, 0.3d0, kind(0.d0))
  Op%O(2,1) = cmplx( 0.5d0,-0.3d0, kind(0.d0))
  Op%O(2,2) = cmplx(-1.d0, 0.d0, kind(0.d0))
  do I = 1, N
     Op%P(I) = I
  enddo
  Op%g     = cmplx(0.d0, 0.2d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Op%type  = 4
  Call Op_set(Op)
  Call Op_clear(Op, N)
  Write(*,*) 'PASS: 2x2 Hermitian operator (type 4)'

  Write(*,*) 'All validation tests passed.'

End Program Test_Op_validate
