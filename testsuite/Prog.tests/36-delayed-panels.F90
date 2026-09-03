! The panel accessors must reconstruct the Green's function exactly.
!
! delayed_update_mod holds G = G_stale + X*Y^T and never forms it. Three
! accessors serve everything the vertex update needs:
!
!   delay_block(P,d)  the d x d block G(P,P)
!   delay_row (P,d)   d rows,    row(i,l) = G(P(l), i)
!   delay_col (P,d)   d columns, col(i,l) = G(i, P(l))
!
! Each is compared here against an explicitly formed G_stale + X*Y^T. A whole-run
! comparison would catch an error in any of them too -- the chain diverges
! immediately -- but it cannot say which, and it cannot reach the edges: an empty
! panel, a panel one append short of the flush threshold, and the partial panel
! left at delay_close. Those are what this covers.
!
! Requires ALF_DELAY_K in the environment; CMake sets it. Without it delay_alloc
! is a no-op by design and every assertion below would pass vacuously, so the
! test refuses to run rather than reporting success.

Program DelayedPanels

   Use Operator_mod
   Use Fields_mod
   Use delayed_update_mod

   Implicit None

   Integer, Parameter :: Ndim = 8, N_FL = 1, kdepth = 8, dmax = 2

   Complex (Kind=Kind(0.D0)), Allocatable :: GR(:,:,:), g_full(:,:), g_start(:,:)
   Complex (Kind=Kind(0.D0)), Allocatable :: xc(:,:), yc(:,:)
   Complex (Kind=Kind(0.D0)), Allocatable :: blk(:,:), rows(:,:), cols(:,:)
   Complex (Kind=Kind(0.D0)) :: one, alpha
   Real (Kind=Kind(0.D0)) :: err, scale
   Integer :: i, j, l, step, d, nfail, P(dmax)

   one   = cmplx(1.d0, 0.d0, kind(0.D0))
   nfail = 0
   d     = 2
   P     = [3, 6]

   If (delay_depth(Ndim) /= kdepth) Then
      Write (*,*) "ERROR: ALF_DELAY_K must be", kdepth, "for this test; got", &
         &        delay_depth(Ndim)
      Stop 2
   End If

   Allocate (GR(Ndim,Ndim,N_FL), g_full(Ndim,Ndim), g_start(Ndim,Ndim))
   Allocate (xc(Ndim,dmax), yc(Ndim,dmax))
   Allocate (blk(d,d), rows(Ndim,d), cols(Ndim,d))

   Do j = 1, Ndim
      Do i = 1, Ndim
         GR(i,j,1) = cmplx(sin(dble(i+2*j)), cos(dble(3*i-j)), kind(0.D0))
      End Do
   End Do
   g_full  = GR(:,:,1)
   g_start = GR(:,:,1)

   Call delay_alloc(Ndim, N_FL, dmax)
   Call delay_open()

   If (.not. delay_active) Then
      Write (*,*) "ERROR: delay_open did not open a region"
      Stop 2
   End If

   ! Walk past the flush threshold: with d = 2 and k = 8 the panel fills on the
   ! fourth append, so this covers an empty panel (step 1), partial ones, the
   ! flush itself, and a partial panel again afterwards.
   Do step = 0, 5

      Call delay_block(1, GR, P, d, blk, d)
      Call delay_row  (1, GR, P, d, rows)
      Call delay_col  (1, GR, P, d, cols)

      scale = maxval(abs(g_full))

      err = 0.d0
      Do j = 1, d
         Do i = 1, d
            err = max(err, abs(blk(i,j) - g_full(P(i),P(j))))
         End Do
      End Do
      If (err/scale > 1.d-13) Then
         Write (*,*) "ERROR step", step, "delay_block rel err", err/scale
         nfail = nfail + 1
      End If

      err = 0.d0
      Do l = 1, d
         Do i = 1, Ndim
            err = max(err, abs(rows(i,l) - g_full(P(l),i)))
         End Do
      End Do
      If (err/scale > 1.d-13) Then
         Write (*,*) "ERROR step", step, "delay_row rel err", err/scale
         nfail = nfail + 1
      End If

      err = 0.d0
      Do l = 1, d
         Do i = 1, Ndim
            err = max(err, abs(cols(i,l) - g_full(i,P(l))))
         End Do
      End Do
      If (err/scale > 1.d-13) Then
         Write (*,*) "ERROR step", step, "delay_col rel err", err/scale
         nfail = nfail + 1
      End If

      ! Apply one more rank-d update, to the panels and to the reference alike.
      Do l = 1, d
         Do i = 1, Ndim
            xc(i,l) = cmplx(0.2d0*sin(dble(i+l+step)), 0.2d0*cos(dble(i*l)), kind(0.D0))
            yc(i,l) = cmplx(0.2d0*cos(dble(2*i-l)), 0.2d0*sin(dble(i+step)), kind(0.D0))
         End Do
      End Do
      alpha = cmplx(0.3d0, 0.1d0, kind(0.D0))

      Call delay_append(1, alpha, xc, yc, d, GR)
      Call ZGEMM('N','T',Ndim,Ndim,d,alpha,xc,Ndim,yc,Ndim,one,g_full,Ndim)

   End Do

   ! Whatever is still pending must be paid for. Leaving it out would credit the
   ! chain with updates it never applied, and the partial panel here is exactly
   ! the case a whole number of flush periods would hide.
   If (delay_pending(1) == 0) Then
      Write (*,*) "NO POWER: nothing pending at delay_close; the trailing flush is untested"
      nfail = nfail + 1
   End If

   Call delay_close(GR)

   If (delay_active) Then
      Write (*,*) "ERROR: delay_close left the region open"
      nfail = nfail + 1
   End If

   err   = maxval(abs(GR(:,:,1) - g_full))
   scale = maxval(abs(g_full))
   If (err/scale > 1.d-13) Then
      Write (*,*) "ERROR: after delay_close, GR differs by", err/scale
      nfail = nfail + 1
   End If

   ! The updates must have moved the matrix, or none of the above has any power:
   ! if the panels stayed near zero, a reconstruction that ignored them entirely
   ! would still agree with the reference everywhere.
   err = maxval(abs(g_full - g_start))/scale
   If (err < 1.d-3) Then
      Write (*,*) "NO POWER: the updates moved G by only", err, "relative"
      nfail = nfail + 1
   End If

   Call delay_dealloc()
   Deallocate (GR, g_full, g_start, xc, yc, blk, rows, cols)

   If (nfail > 0) Then
      Write (*,*) "FAILURES:", nfail
      Stop 2
   End If

   Write (*,*) "SUCCESS"

End Program DelayedPanels
