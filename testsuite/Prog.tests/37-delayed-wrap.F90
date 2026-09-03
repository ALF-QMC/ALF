! Op_Wrap_panels must move a factored Green's function exactly as Op_Wrapup and
! Op_Wrapdo move an explicit one.
!
! The delayed update holds G = G_stale + X*Y^T across a time slice. Between two
! vertex updates ALF conjugates the Green's function with the vertex operator,
! and both wraps act as G -> L*G*R with L and R the identity outside Op%P. So
!
!     L*(G_stale + X*Y^T)*R = (L*G_stale*R) + (L*X)*(R^T*Y)^T
!
! must hold exactly, and Op_Wrap_panels is what applies L to X and R^T to Y.
! This test forms both sides and compares them.
!
! Why it exists rather than relying on a simulation: for a Hamiltonian whose
! vertices all have Op%diag = .true. -- which is the common case -- a production
! run never reaches the ZSLGEMM branches or the N_type = 2 branches at all. Half
! of Op_Wrap_panels would otherwise be shipped untested. The sweep below covers
! every combination of type (1..4), rank (1..4), N_type (1,2), side (u,d),
! diagonality and time-dependent coupling, which is all sixteen branches of each
! wrap.
!
! g_t is an axis rather than a detail. When Op%g_t is allocated the wraps take
! the coupling from it per time slice and build the exponential on the spot,
! bypassing the E_Exp table the other arms read; Op_Wrap_panels mirrors that in
! four more branches, which a Hamiltonian that never allocates g_t leaves
! unexecuted. The value below is deliberately not
! Op%g: a panel arm that ignored g_t and read the table instead would otherwise
! agree with the reference for the wrong reason.
!
! It is also the drift guard: Op_Wrap_panels reproduces the branch structure of
! Op_Wrapup/Op_Wrapdo, and if either is edited without the other they diverge
! silently. Nothing else in the suite would notice.

Program DelayedWrap

   Use Operator_mod
   Use Fields_mod

   Implicit None

   Complex (Kind=Kind(0.D0)), Dimension(:,:), Allocatable :: g_ref, g_pan, g_before
   Complex (Kind=Kind(0.D0)), Dimension(:,:), Allocatable :: xp, yp, xr, yr
   Complex (Kind=Kind(0.D0)) :: one, zero
   Real    (Kind=Kind(0.D0)) :: err, scale, moved
   Integer :: Ndim, ncols, i, j, n, opn, nt, N_Type, idiag, iud, igt
   Integer :: nfail
   Character (Len=1) :: updo
   Type (Operator) :: Op
   Type (Fields)   :: nsigma_single

   Ndim  = 6
   ncols = 3
   one   = cmplx(1.d0, 0.d0, kind(0.D0))
   zero  = cmplx(0.d0, 0.d0, kind(0.D0))
   nfail = 0

   ! Without this Phi_st is zero, every exponent is exp(0) = 1, and the wraps
   ! reduce to the identity for the discrete field types -- at which point the
   ! comparison below would pass on any implementation. The power check catches
   ! that, and this is the fix it points to.
   Call Fields_init()
   Call nsigma_single%make(1,1)

   Do nt = 1, 4
   Do opn = 1, 4
   Do idiag = 1, 2
   Do N_Type = 1, 2
   Do iud = 1, 2
   Do igt = 1, 2

      updo = 'u'
      if (iud == 2) updo = 'd'

      Allocate (g_ref(Ndim,Ndim), g_pan(Ndim,Ndim), g_before(Ndim,Ndim))
      Allocate (xp(Ndim,ncols), yp(Ndim,ncols), xr(Ndim,ncols), yr(Ndim,ncols))

      Call Op_make (Op, opn)
      ! Offset from the leading block, so an index computed on the wrong basis
      ! lands somewhere visible instead of aliasing onto the right answer.
      Do i = 1, Op%n
         Op%P(i) = i + 1
      End Do
      Do i = 1, Op%n
         Do n = 1, Op%n
            if (i == n) then
               Op%O(i,n) = cmplx(0.5d0*dble(i), 0.d0, kind(0.D0))
            elseif (idiag == 2) then
               ! Hermitian off-diagonal, so Op_set's diagonalisation is valid and
               ! Op%diag comes out false -- the branch a production run never takes.
               Op%O(i,n) = cmplx(0.25d0*dble(i+n), 0.1d0*dble(i-n), kind(0.D0))
               Op%O(n,i) = conjg(Op%O(i,n))
            else
               Op%O(i,n) = zero
            endif
         End Do
      End Do

      Op%type  = nt
      Op%g     = cmplx(0.7d0, 0.d0, kind(0.D0))
      Op%alpha = zero
      ! Before Op_set, which is what latches g_t_alloc (Operator_mod.F90:267).
      ! Distinct from Op%g so the table arms and the g_t arms cannot agree by
      ! coincidence.
      If (igt == 2) Then
         Allocate (Op%g_t(1))
         Op%g_t(1) = cmplx(0.41d0, 0.23d0, kind(0.D0))
      End If
      Call Op_set (Op)

      Select Case (nt)
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

      ! Entries of order one with no dominant diagonal: nothing here factorises,
      ! and a wrong index lands somewhere visible.
      Do j = 1, Ndim
         Do i = 1, Ndim
            g_ref(i,j) = cmplx(sin(dble(i+2*j)), cos(dble(3*i-j)), kind(0.D0))
         End Do
      End Do
      Do j = 1, ncols
         Do i = 1, Ndim
            xp(i,j) = cmplx(0.3d0*sin(dble(i*j+1)), 0.3d0*cos(dble(i-j)), kind(0.D0))
            yp(i,j) = cmplx(0.3d0*cos(dble(2*i+j)), 0.3d0*sin(dble(i+3*j)), kind(0.D0))
         End Do
      End Do
      g_pan = g_ref
      xr = xp
      yr = yp

      ! Reference: form the explicit Green's function first, then wrap it.
      Call ZGEMM('N','T',Ndim,Ndim,ncols,one,xp,Ndim,yp,Ndim,one,g_ref,Ndim)
      g_before = g_ref
      if (updo == 'u') then
         Call Op_Wrapup(g_ref, Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
      else
         Call Op_Wrapdo(g_ref, Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
      endif

      ! Candidate: wrap the stale matrix and the panels separately, then form it.
      if (updo == 'u') then
         Call Op_Wrapup(g_pan, Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
      else
         Call Op_Wrapdo(g_pan, Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
      endif
      Call Op_Wrap_panels(xr, yr, Op, nsigma_single%f(1,1), Ndim, ncols, N_Type, 1, updo)
      Call ZGEMM('N','T',Ndim,Ndim,ncols,one,xr,Ndim,yr,Ndim,one,g_pan,Ndim)

      scale = maxval(abs(g_ref))
      err   = maxval(abs(g_ref - g_pan))/max(scale, 1.d-30)

      ! How far the wrap moved the matrix at all. Without this the comparison has
      ! no power: an operator whose wrap were near the identity would make both
      ! sides agree whatever Op_Wrap_panels computed. N_type = 2 on a diagonal
      ! operator is a genuine no-op in both wraps, so that case is exempt rather
      ! than expected to move.
      moved = maxval(abs(g_ref - g_before))/max(scale, 1.d-30)
      If (moved < 1.d-6 .and. .not. (N_Type == 2 .and. Op%diag)) Then
         Write (*,*) "NO POWER type", nt, "rank", opn, "N_type", N_Type, &
            &        "diag", Op%diag, "side ", updo, " g_t", igt == 2, &
            &        " moved", moved
         nfail = nfail + 1
      End If

      If (err > 1.d-12) Then
         Write (*,*) "ERROR type", nt, "rank", opn, "N_type", N_Type, &
            &        "diag", Op%diag, "side ", updo, " g_t", igt == 2, &
            &        " rel err", err
         nfail = nfail + 1
      End If

      ! Op_clear deallocates g_t when g_t_alloc is set, so the next iteration
      ! starts from an operator that has none.
      Call Op_clear(Op, opn)
      Deallocate (g_ref, g_pan, g_before, xp, yp, xr, yr)

   End Do
   End Do
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

End Program DelayedWrap