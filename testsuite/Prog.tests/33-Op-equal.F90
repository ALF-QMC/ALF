! Test Op_equal: verify it returns .true. for identical operators
! and .false. when any defining field differs.
!
Program Test_Op_equal
  Use Operator_mod
  Implicit None

  Type (Operator) :: Op_a, Op_b
  Integer :: I, N

  N = 3

  ! --- Setup two identical 3x3 Hermitian operators ---
  Call Op_make(Op_a, N)
  Call Op_make(Op_b, N)

  Op_a%O(1,1) = cmplx( 2.d0,  0.d0, kind(0.d0))
  Op_a%O(1,2) = cmplx(-1.d0,  0.5d0, kind(0.d0))
  Op_a%O(2,1) = cmplx(-1.d0, -0.5d0, kind(0.d0))
  Op_a%O(2,2) = cmplx( 3.d0,  0.d0, kind(0.d0))
  Op_a%O(2,3) = cmplx( 0.d0,  1.d0, kind(0.d0))
  Op_a%O(3,2) = cmplx( 0.d0, -1.d0, kind(0.d0))
  Op_a%O(3,3) = cmplx( 1.d0,  0.d0, kind(0.d0))
  do I = 1, N
     Op_a%P(I) = I
  enddo
  Op_a%g     = cmplx(0.5d0, 0.d0, kind(0.d0))
  Op_a%alpha = cmplx(0.1d0, 0.2d0, kind(0.d0))
  Op_a%type  = 2

  ! Copy defining fields to Op_b
  Op_b%O     = Op_a%O
  do I = 1, N
     Op_b%P(I) = Op_a%P(I)
  enddo
  Op_b%g     = Op_a%g
  Op_b%alpha = Op_a%alpha
  Op_b%type  = Op_a%type

  ! --- Test 1: identical operators must be equal ---
  if (.not. Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: identical operators reported as unequal'
     STOP 2
  endif
  write(*,*) 'PASS: identical operators are equal'

  ! --- Test 2: different O matrix ---
  Op_b%O(1,2) = cmplx(99.d0, 0.d0, kind(0.d0))
  if (Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: different O matrix reported as equal'
     STOP 2
  endif
  Op_b%O(1,2) = Op_a%O(1,2)  ! restore
  write(*,*) 'PASS: different O detected'

  ! --- Test 3: different P ---
  Op_b%P(2) = 99
  if (Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: different P reported as equal'
     STOP 2
  endif
  Op_b%P(2) = Op_a%P(2)  ! restore
  write(*,*) 'PASS: different P detected'

  ! --- Test 4: different g ---
  Op_b%g = cmplx(9.d0, 0.d0, kind(0.d0))
  if (Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: different g reported as equal'
     STOP 2
  endif
  Op_b%g = Op_a%g  ! restore
  write(*,*) 'PASS: different g detected'

  ! --- Test 5: different alpha ---
  Op_b%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  if (Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: different alpha reported as equal'
     STOP 2
  endif
  Op_b%alpha = Op_a%alpha  ! restore
  write(*,*) 'PASS: different alpha detected'

  ! --- Test 6: different type ---
  Op_b%type = 3
  if (Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: different type reported as equal'
     STOP 2
  endif
  Op_b%type = Op_a%type  ! restore
  write(*,*) 'PASS: different type detected'

  ! --- Test 7: one has g_t allocated, the other does not ---
  allocate(Op_a%g_t(4))
  Op_a%g_t = cmplx(1.d0, 0.d0, kind(0.d0))
  if (Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: mismatched g_t allocation reported as equal'
     STOP 2
  endif
  write(*,*) 'PASS: mismatched g_t allocation detected'

  ! --- Test 8: both have g_t but different values ---
  allocate(Op_b%g_t(4))
  Op_b%g_t = cmplx(1.d0, 0.d0, kind(0.d0))
  if (.not. Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: identical g_t reported as unequal'
     STOP 2
  endif
  Op_b%g_t(3) = cmplx(99.d0, 0.d0, kind(0.d0))
  if (Op_equal(Op_a, Op_b)) then
     write(*,*) 'FAIL: different g_t values reported as equal'
     STOP 2
  endif
  write(*,*) 'PASS: different g_t values detected'

  write(*,*) 'success'

end Program Test_Op_equal
