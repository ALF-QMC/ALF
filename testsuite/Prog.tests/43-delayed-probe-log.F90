! "auto" must run the probe and log the choice it made, whichever way it went.
!
! With ALF_DELAY_K=auto (CMake sets it) delay_resolve runs the real probe. Check
! that:
!
!   the probe ran           the source is 'probe', or 'formula' on a flat curve
!   a zero depth            only the probe can have chosen it, at this Ndim
!   the immediate update    has exactly one row in the table, as k = 0
!   the depth in force      is marked, once, on its own row: the immediate row
!                           when the depth is zero
!
! It is also the one test that executes delay_probe at all; the others either
! fix the depth or test the decision on made-up curves.

Program DelayedProbeLog

   Use delayed_update_mod

   Implicit None

   ! Ndim wide enough for several candidates, small enough for a quick probe.
   Integer, Parameter :: Ndim = 64, N_FL = 1, dmax = 1

   Character (Len=256) :: line
   Integer :: k, u, ios, n_imm, n_mark, k_row, nfail
   Logical :: mark_on_k

   nfail = 0

   k = delay_resolve(Ndim, dmax)

   If (trim(delay_source) /= 'probe' .and. trim(delay_source) /= 'formula') Then
      Write (*,*) "ERROR: ALF_DELAY_K=auto resolved from '", trim(delay_source), &
         &        "'; is ALF_DELAY_K set?"
      Stop 2
   End If
   If (k == 0 .and. trim(delay_source) /= 'probe') Then
      Write (*,*) "ERROR: depth 0 from '", trim(delay_source), "' at Ndim", Ndim
      nfail = nfail + 1
   End If

   Call delay_alloc(Ndim, N_FL, dmax, k)

   Open (newunit=u, status='scratch', action='readwrite', form='formatted')
   Call delay_log(u)
   Rewind (u)
   n_imm     = 0
   n_mark    = 0
   mark_on_k = .false.
   Do
      Read (u, '(a)', iostat=ios) line
      If (ios /= 0) Exit
      Write (*,'(a)') trim(line)
      ! Table rows open with their depth, the immediate one with 0; header
      ! lines open with a word, and one of them can say "(immediate)" too.
      Read (line, *, iostat=ios) k_row
      If (ios /= 0) Cycle
      If (k_row == 0 .and. index(line, '(immediate') > 0) n_imm = n_imm + 1
      If (index(line, '<- ') > 0) Then
         n_mark = n_mark + 1
         If (k_row == k) mark_on_k = .true.
      End If
   End Do
   Close (u)

   If (n_imm /= 1) Then
      Write (*,*) "ERROR: expected one immediate row in the table, found", n_imm
      nfail = nfail + 1
   End If
   If (n_mark /= 1 .or. .not. mark_on_k) Then
      Write (*,*) "ERROR: expected one marker, on the row for k =", k, &
         &        "; found", n_mark, "marker(s), on k:", mark_on_k
      nfail = nfail + 1
   End If

   Call delay_dealloc()

   If (nfail > 0) Then
      Write (*,*) "FAILURES:", nfail
      Stop 2
   End If

   Write (*,*) "SUCCESS: k =", k, "from ", trim(delay_source)

End Program DelayedProbeLog