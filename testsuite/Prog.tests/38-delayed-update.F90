! The delayed update, composed, must track an immediate one exactly.
!
! Test 36 drives the panel accessors and test 37 drives the panel wrap, each in
! isolation. What neither covers is the two running *together*: a slice in which
! reads, rank-d appends, flushes and vertex conjugations interleave, which is
! what Wrapgr_mod actually does. A Hamiltonian whose Op_V are all diagonal and
! of rank at most 2 -- which is the common case -- never composes the wrap's
! ZSLGEMM arms with an append at all, so this test drives ranks 1 to 4 and both
! diagonalities explicitly.
!
! The arms:
!
!   immediate   an explicit G. Each update is a ZGEMM onto it, each wrap an
!               Op_Wrapup/Op_Wrapdo of it.
!   delayed     G_stale plus panels. Each update is delay_append, each wrap a
!               delay_wrap, and the accessors are read between them.
!
! Both see the same operators, fields, sides and updates in the same order, so
! they must agree: at every step through the accessors, and at delay_close on
! the whole matrix. The wraps are what make this more than test 36 with a bigger
! d -- the identity being checked is
!
!     L*(G_stale + X*Y^T)*R = (L*G_stale*R) + (L*X)*(R^T*Y)^T
!
! applied repeatedly and with appends landing between the conjugations, where a
! panel that was wrapped the wrong number of times, or in the wrong order
! against a flush, gives a plausible wrong answer rather than an obvious one.
!
! Requires ALF_DELAY_K in the environment; CMake sets it. Without it delay_alloc
! is a no-op and every assertion here would pass vacuously, so the test refuses
! to run rather than reporting success.

Program DelayedUpdate

   Use Operator_mod
   Use Fields_mod
   Use delayed_update_mod

   Implicit None

   Integer, Parameter :: Ndim = 10, N_FL = 1, kdepth = 8, dmax = 4

   Complex (Kind=Kind(0.D0)), Allocatable :: GR(:,:,:), g_imm(:,:)
   Complex (Kind=Kind(0.D0)), Allocatable :: xc(:,:), yc(:,:)
   Complex (Kind=Kind(0.D0)), Allocatable :: blk(:,:), rows(:,:), cols(:,:)
   Complex (Kind=Kind(0.D0)) :: one, alpha
   Real    (Kind=Kind(0.D0)) :: err, scale, moved, before
   Integer :: i, j, l, n, step, d, nfail, idiag, iud, nt_f, N_Type
   Integer :: P(dmax)
   Character (Len=1) :: updo
   Type (Operator) :: Op
   Type (Fields)   :: nsigma_single

   one   = cmplx(1.d0, 0.d0, kind(0.D0))
   nfail = 0

   If (delay_depth(Ndim) /= kdepth) Then
      Write (*,*) "ERROR: ALF_DELAY_K must be", kdepth, "for this test; got", &
         &        delay_depth(Ndim)
      Stop 2
   End If

   ! Without this Phi_st is zero, every exponent is exp(0) and the wraps reduce
   ! to the identity -- at which point the comparison would pass on any panel
   ! wrap whatsoever. The power check below is what catches that.
   Call Fields_init()
   Call nsigma_single%make(1,1)

   ! Rank 1 and 2 are the shapes production runs, 3 and 4 the ones only a
   ! foreign Hamiltonian would reach; all four run here because upgrade_mod's
   ! append path branches on d and the wrap's does not.
   Do d = 1, dmax
   Do idiag = 1, 2
   Do iud = 1, 2
   Do N_Type = 1, 2

      updo = 'u'
      if (iud == 2) updo = 'd'
      nt_f = 1 + mod(d, 4)

      Allocate (GR(Ndim,Ndim,N_FL), g_imm(Ndim,Ndim))
      Allocate (xc(Ndim,dmax), yc(Ndim,dmax))
      Allocate (blk(d,d), rows(Ndim,d), cols(Ndim,d))

      ! Support offset from the leading block, so an index computed on the wrong
      ! basis lands somewhere visible rather than aliasing onto the right answer.
      Do i = 1, d
         P(i) = 2*i
      End Do

      Call Op_make (Op, d)
      Do i = 1, Op%n
         Op%P(i) = P(i)
      End Do
      Do i = 1, Op%n
         Do n = 1, Op%n
            if (i == n) then
               Op%O(i,n) = cmplx(0.5d0*dble(i), 0.d0, kind(0.D0))
            elseif (idiag == 2) then
               ! Hermitian off-diagonal, so Op_set's diagonalisation is valid and
               ! Op%diag comes out false -- the wrap branch a production run
               ! never takes, here composed with an append for the first time.
               Op%O(i,n) = cmplx(0.25d0*dble(i+n), 0.1d0*dble(i-n), kind(0.D0))
               Op%O(n,i) = conjg(Op%O(i,n))
            else
               Op%O(i,n) = cmplx(0.d0, 0.d0, kind(0.D0))
            endif
         End Do
      End Do
      Op%type  = nt_f
      Op%g     = cmplx(0.7d0, 0.d0, kind(0.D0))
      Op%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
      Call Op_set (Op)

      Select Case (nt_f)
      Case (1)
         nsigma_single%f(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
      Case (2)
         nsigma_single%f(1,1) = cmplx(2.d0, 0.d0, kind(0.D0))
      Case (3)
         nsigma_single%f(1,1) = cmplx(0.813d0, 0.d0, kind(0.D0))
      Case (4)
         nsigma_single%f(1,1) = cmplx(-1.d0, 0.5d0, kind(0.D0))
      End Select
      nsigma_single%t(1) = Op%type

      Do j = 1, Ndim
         Do i = 1, Ndim
            GR(i,j,1) = cmplx(sin(dble(i+2*j)), cos(dble(3*i-j)), kind(0.D0))
         End Do
      End Do
      g_imm = GR(:,:,1)

      Call delay_alloc(Ndim, N_FL, dmax)
      Call delay_open()

      If (.not. delay_active) Then
         Write (*,*) "ERROR: delay_open did not open a region"
         Stop 2
      End If

      moved = 0.d0

      ! Eleven rounds of (read, append, wrap), and the count is load-bearing.
      ! Every rank must both cross the flush threshold at least once and end
      ! holding a partial panel, or the trailing flush at delay_close goes
      ! untested -- which is exactly the case a whole number of flush periods
      ! hides. With k = 8 the panel empties when the append count is a multiple
      ! of 3 at d = 3 and of 2 at d = 4, so the count has to be odd, not a
      ! multiple of 3, and past 8 for d = 1. Eleven is the smallest such:
      ! d = 1 leaves 3 columns pending, d = 2 leaves 6, d = 3 leaves 6, d = 4
      ! leaves 4, and all four flush along the way.
      Do step = 0, 10

         Call delay_block(1, GR, P, d, blk, d)
         Call delay_row  (1, GR, P, d, rows)
         Call delay_col  (1, GR, P, d, cols)

         scale = maxval(abs(g_imm))

         err = 0.d0
         Do j = 1, d
            Do i = 1, d
               err = max(err, abs(blk(i,j) - g_imm(P(i),P(j))))
            End Do
         End Do
         If (err/scale > 1.d-12) Then
            Write (*,*) "ERROR d", d, "diag", Op%diag, "side ", updo, "N_type", &
               &        N_Type, "step", step, "block rel err", err/scale
            nfail = nfail + 1
         End If

         err = 0.d0
         Do l = 1, d
            Do i = 1, Ndim
               err = max(err, abs(rows(i,l) - g_imm(P(l),i)))
            End Do
         End Do
         If (err/scale > 1.d-12) Then
            Write (*,*) "ERROR d", d, "diag", Op%diag, "side ", updo, "N_type", &
               &        N_Type, "step", step, "row rel err", err/scale
            nfail = nfail + 1
         End If

         err = 0.d0
         Do l = 1, d
            Do i = 1, Ndim
               err = max(err, abs(cols(i,l) - g_imm(i,P(l))))
            End Do
         End Do
         If (err/scale > 1.d-12) Then
            Write (*,*) "ERROR d", d, "diag", Op%diag, "side ", updo, "N_type", &
               &        N_Type, "step", step, "col rel err", err/scale
            nfail = nfail + 1
         End If

         ! One rank-d update, to both arms.
         Do l = 1, d
            Do i = 1, Ndim
               xc(i,l) = cmplx(0.2d0*sin(dble(i+l+step)), &
                  &            0.2d0*cos(dble(i*l)), kind(0.D0))
               yc(i,l) = cmplx(0.2d0*cos(dble(2*i-l)), &
                  &            0.2d0*sin(dble(i+step)), kind(0.D0))
            End Do
         End Do
         alpha = cmplx(0.3d0, 0.1d0, kind(0.D0))

         Call delay_append(1, alpha, xc, yc, d, GR)
         Call ZGEMM('N','T',Ndim,Ndim,d,alpha,xc,Ndim,yc,Ndim,one,g_imm,Ndim)

         ! Then one conjugation, to both arms. This is the composition the other
         ! tests do not make: the panels are wrapped while holding appends that
         ! have not been flushed.
         before = maxval(abs(g_imm))
         ! Paired, as Wrapgr_mod pairs them: delay_wrap conjugates the panels
         ! only, so the stale matrix takes the ordinary wrap beside it. Both
         ! halves of G_stale + X*Y^T have to move or the sum is not conjugated.
         if (updo == 'u') then
            Call Op_Wrapup(GR(:,:,1), Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
            Call delay_wrap(1, Op, nsigma_single%f(1,1), N_Type, 1, updo)
            Call Op_Wrapup(g_imm, Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
         else
            Call Op_Wrapdo(GR(:,:,1), Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
            Call delay_wrap(1, Op, nsigma_single%f(1,1), N_Type, 1, updo)
            Call Op_Wrapdo(g_imm, Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
         endif
         moved = max(moved, abs(maxval(abs(g_imm)) - before)/max(before, 1.d-30))

      End Do

      ! The trailing flush: whatever is still pending has to be paid for, and a
      ! whole number of flush periods would hide it.
      If (delay_pending(1) == 0) Then
         Write (*,*) "NO POWER d", d, ": nothing pending at delay_close"
         nfail = nfail + 1
      End If

      Call delay_close(GR)

      If (delay_active) Then
         Write (*,*) "ERROR: delay_close left the region open"
         nfail = nfail + 1
      End If

      ! The whole matrix, once the panels are flushed into it.
      scale = maxval(abs(g_imm))
      err   = maxval(abs(GR(:,:,1) - g_imm))/max(scale, 1.d-30)
      If (err > 1.d-12) Then
         Write (*,*) "ERROR d", d, "diag", Op%diag, "side ", updo, "N_type", &
            &        N_Type, "flushed rel err", err
         nfail = nfail + 1
      End If

      ! Without this the comparison has no power: a wrap near the identity would
      ! make both arms agree whatever Op_Wrap_panels computed. N_type = 2 on a
      ! diagonal operator is a genuine no-op in both wraps, so it is exempt.
      If (moved < 1.d-9 .and. .not. (N_Type == 2 .and. Op%diag)) Then
         Write (*,*) "NO POWER d", d, "diag", Op%diag, "side ", updo, "N_type", &
            &        N_Type, ": the wraps did not move the matrix"
         nfail = nfail + 1
      End If

      Call delay_dealloc()
      Call Op_clear(Op, d)
      Deallocate (GR, g_imm, xc, yc, blk, rows, cols)

   End Do
   End Do
   End Do
   End Do

   Call nsigma_single%clear()

   If (nfail > 0) Then
      Write (*,*) "FAILURES:", nfail
      Stop 2
   End If

   Write (*,*) "SUCCESS"

End Program DelayedUpdate
