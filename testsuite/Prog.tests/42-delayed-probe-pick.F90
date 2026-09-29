! The probe's decision must follow its rules on any cost curve.
!
! delay_probe times kernels and hands the resulting curve, and the time of one
! immediate update, to delay_probe_pick. The timings cannot be asserted on --
! they belong to the machine -- so the decision is tested on its own here,
! against made-up curves over the candidate depths [8, 16, 32, 64, 128, 256]:
!
!   a clear minimum             the largest candidate within the 5% margin of it
!   a flat curve                the closed-form fallback, sqrt(2*Ndim) clamped
!   untimed candidates          ignored, the decision taken on the rest
!   nothing timed               the fallback
!   a zero reading              a failed clock, so the fallback
!
! and against the immediate update, which a rank-d update costs once where the
! best candidate costs d times its per-column cost:
!
!   immediate cheaper           0, no delay
!   within the margin           0: ties go to the immediate update
!   just past the margin        the delay, as without the comparison
!   the rank                    one curve and t_imm, delay at d = 1, not at d = 2
!   flat curve, immediate wins  0: the comparison comes before the fallback
!   a single candidate          compared all the same (8 <= Ndim < 16)
!   t_imm unreadable            the fallback

Program DelayedProbePick

   Use delayed_update_mod

   Implicit None

   Real (Kind=Kind(0.D0)), Parameter :: UNTIMED = huge(1.d0), SLOW = 100.d0
   Real (Kind=Kind(0.D0)), Parameter :: CLEAR(6) = [4.d0, 2.d0, 1.d0, 1.04d0, 1.5d0, 3.d0]
   Real (Kind=Kind(0.D0)), Parameter :: FLAT(6)  = [1.d0, 1.01d0, 1.02d0, 1.d0, 1.03d0, 1.01d0]
   Integer :: nfail

   nfail = 0

   ! The curve alone: the immediate update far slower, so never chosen.
   Call check("clear minimum", CLEAR, SLOW, 1, 256, 64, 'probe')
   Call check("flat curve", FLAT, SLOW, 1, 256, formula(256), 'formula')
   Call check("untimed above Ndim", [2.d0, 1.d0, UNTIMED, UNTIMED, UNTIMED, UNTIMED], &
      &       SLOW, 1, 16, 16, 'probe')
   Call check("nothing timed", [UNTIMED, UNTIMED, UNTIMED, UNTIMED, UNTIMED, UNTIMED], &
      &       SLOW, 1, 256, formula(256), 'formula')
   Call check("zero reading", [0.d0, 2.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
      &       SLOW, 1, 256, formula(256), 'formula')

   ! Against the immediate update; the best candidate of CLEAR costs 1 per column.
   Call check("immediate cheaper", CLEAR, 0.5d0, 1, 256, 0, 'probe')
   Call check("within the margin", CLEAR, 1.04d0, 1, 256, 0, 'probe')
   Call check("just past the margin", CLEAR, 1.06d0, 1, 256, 64, 'probe')
   Call check("rank 1", CLEAR, 1.5d0, 1, 256, 64, 'probe')
   Call check("rank 2", CLEAR, 1.5d0, 2, 256, 0, 'probe')
   Call check("flat curve, immediate wins", FLAT, 0.9d0, 1, 256, 0, 'probe')
   Call check("single candidate, delay wins", [1.d0, UNTIMED, UNTIMED, UNTIMED, UNTIMED, UNTIMED], &
      &       SLOW, 1, 12, formula(12), 'formula')
   Call check("single candidate, immediate wins", [1.d0, UNTIMED, UNTIMED, UNTIMED, UNTIMED, UNTIMED], &
      &       0.5d0, 1, 12, 0, 'probe')
   Call check("t_imm zero", CLEAR, 0.d0, 1, 256, formula(256), 'formula')
   Call check("t_imm untimed", CLEAR, UNTIMED, 1, 256, formula(256), 'formula')

   If (nfail > 0) Then
      Write (*,*) "FAILURES:", nfail
      Stop 2
   End If

   Write (*,*) "SUCCESS"

Contains

   Subroutine check(name, cost, t_imm, d, Ndim, k_expected, source_expected)
      Character (Len=*),      Intent(In) :: name, source_expected
      Real (Kind=Kind(0.D0)), Intent(In) :: cost(:), t_imm
      Integer,                Intent(In) :: d, Ndim, k_expected
      Character (Len=16) :: source
      Integer :: k

      k = delay_probe_pick(cost, t_imm, d, Ndim, source)
      If (k /= k_expected .or. trim(source) /= source_expected) Then
         Write (*,*) "ERROR ", name, ": got k =", k, " from ", trim(source), &
            &        "; expected k =", k_expected, " from ", source_expected
         nfail = nfail + 1
      End If
   End Subroutine check

   ! The closed form, restated rather than called: nint(sqrt(2*Ndim)) in [8, 256].
   Integer Function formula(Ndim)
      Integer, Intent(In) :: Ndim
      formula = min(256, max(8, nint(sqrt(2.d0*real(Ndim, Kind(0.D0))))))
   End Function formula

End Program DelayedProbePick