! A fixed delay depth must be capped at Ndim.
!
! ALF_DELAY_K takes any positive integer. Past Ndim the panels outgrow G and buy
! nothing, and near huge(0) the panel width kmax + dmax will overflow. Two
! places now cap the depth, and each is checked here.
!
!   delay_resolve  returns min(ALF_DELAY_K, Ndim) and notes a cap in delay_log.
!                  CMake runs this program once per request, passing the depth
!                  expected back and whether a cap should be reported.
!   delay_alloc    caps a depth passed to it directly. It is given huge(0), and
!                  the panels must then flush at exactly Ndim and flush into GR
!                  the same update an explicit rank-1 sum gives.
!
! Usage: 40-delayed-depth-cap <expected k> <1 if capped, else 0>

Program DelayedDepthCap

   Use delayed_update_mod

   Implicit None

   Integer, Parameter :: Ndim = 16, N_FL = 1, dmax = 2

   Complex (Kind=Kind(0.D0)), Allocatable :: GR(:,:,:), g_imm(:,:), xc(:,:), yc(:,:)
   Complex (Kind=Kind(0.D0)) :: one, alpha
   Real (Kind=Kind(0.D0)) :: err, scale
   Integer :: i, j, step, k, k_expected, capped_expected, nfail, u, ios
   Logical :: capped_seen
   Character (Len=256) :: line, arg

   one   = cmplx(1.d0, 0.d0, kind(0.D0)) ! complex identity
   nfail = 0

   If (command_argument_count() /= 2) Then
      Write (*,*) "ERROR: usage: 40-delayed-depth-cap <expected k> <capped 0|1>"
      Stop 2
   End If
   Call get_command_argument(1, arg)
   Read (arg, *) k_expected
   Call get_command_argument(2, arg)
   Read (arg, *) capped_expected

   ! delay_resolve: the depth, and whether delay_log says it was capped.
   k = delay_resolve(Ndim, dmax)
   If (k /= k_expected) Then
      Write (*,*) "ERROR: delay_resolve returned", k, "expected", k_expected
      nfail = nfail + 1
   End If

   Open (newunit=u, status='scratch', action='readwrite', form='formatted')
   Call delay_log(u)
   Rewind (u)
   capped_seen = .false.
   Do
      Read (u, '(a)', iostat=ios) line
      If (ios /= 0) Exit
      If (index(line, 'capped at Ndim') > 0) capped_seen = .true.
   End Do
   Close (u)
   If (capped_seen .neqv. (capped_expected == 1)) Then
      Write (*,*) "ERROR: delay_log reports a cap:", capped_seen, &
      &        "expected:", capped_expected == 1
      nfail = nfail + 1
   End If

   ! delay_alloc: huge(0) passed directly, which overflowed kmax + dmax before
   ! the cap.
   Allocate (GR(Ndim,Ndim,N_FL), g_imm(Ndim,Ndim), xc(Ndim,1), yc(Ndim,1))
   Do j = 1, Ndim
      Do i = 1, Ndim
         GR(i,j,1) = cmplx(sin(dble(i+2*j)), cos(dble(3*i-j)), kind(0.D0))
      End Do
   End Do
   g_imm = GR(:,:,1)

   Call delay_alloc(Ndim, N_FL, dmax, huge(0))
   Call delay_open()
   If (.not. delay_active) Then
      Write (*,*) "ERROR: delay_open did not open a region"
      Stop 2
   End If

   ! Rank-1 appends: Ndim - 1 stay pending, the Ndim-th flushes, and three more
   ! leave a partial panel for delay_close.
   Do step = 1, Ndim + 3
      Do i = 1, Ndim
         xc(i,1) = cmplx(0.2d0*sin(dble(i+step)), 0.2d0*cos(dble(i*step)), kind(0.D0))
         yc(i,1) = cmplx(0.2d0*cos(dble(2*i-step)), 0.2d0*sin(dble(i+3*step)), kind(0.D0))
      End Do
      alpha = cmplx(0.3d0, 0.1d0*step, kind(0.D0))

      Call delay_append(1, alpha, xc, yc, 1, GR)
      Call ZGEMM('N','T',Ndim,Ndim,1,alpha,xc,Ndim,yc,Ndim,one,g_imm,Ndim)

      If (step == Ndim - 1 .and. delay_pending(1) /= Ndim - 1) Then
         Write (*,*) "ERROR: flushed before Ndim; pending", delay_pending(1)
         nfail = nfail + 1
      End If
      If (step == Ndim .and. delay_pending(1) /= 0) Then
         Write (*,*) "ERROR: no flush at Ndim; pending", delay_pending(1)
         nfail = nfail + 1
      End If
   End Do

   Call delay_close(GR)

   scale = maxval(abs(g_imm))
   err   = maxval(abs(GR(:,:,1) - g_imm))/max(scale, 1.d-30)
   If (err > 1.d-12) Then
      Write (*,*) "ERROR: flushed rel err", err
      nfail = nfail + 1
   End If

   Call delay_dealloc()
   Deallocate (GR, g_imm, xc, yc)

   If (nfail > 0) Then
      Write (*,*) "FAILURES:", nfail
      Stop 2
   End If

   Write (*,*) "SUCCESS"

End Program DelayedDepthCap
