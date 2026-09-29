!  Copyright (C) 2016 - 2026 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following
!     additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific
!       community. Being part of that community we feel that it is reasonable to
!       require you to give an attribution back to the original authors if you
!       have benefitted from this program. Guidelines for a proper citation can
!       be found on the project's homepage http://alf.physik.uni-wuerzburg.de.
!
!     - We require the preservation of the above copyright notice and this
!       license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source
!       files. To obtain the original source files please visit the homepage
!       http://alf.physik.uni-wuerzburg.de.
!
!     - If you make substantial changes to the program we require you to either
!       consider contributing to the ALF project or to mark your material in a
!       reasonable way as different from the original version.

!-------------------------------------------------------------------------------
!> @brief
!> This module enables us to hold the Green's function in factored form across
!> one imaginary time slice; an accepted update will now append a column pair
!> instead of updating the whole matrix. Thus this is a method of "delaying" the
!> Green's function updates, with the aim of reducing the burden on the memory
!> system of the execution environment. This module is designed only for the
!> sequential auxiliary field sampler.
!>
!> @details
!> When a proposed field change is accepted a rank-d update has to be applied to
!> the entire Green's function (see upgrade_mod). An "instantaneous" scheme
!> requires this rank-d BLAS operation at every such update, which is memory
!> bandwidth bound and hence can cause performance degradation. This is
!> especially true for large matrix dimensions and when the compute environment
!> is fully saturated by many other memory intensive jobs.
!>
!> We can take advantage of the fact that no step of the sequential sweep needs
!> the full Green's function, G. With P = Op%P and d = Op%N_non_zero:
!>
!>    - every proposal, accepted or not, needs only the d x d block G(P,P)
!>      to form the Metropolis ratio;
!>
!>    - an accepted proposal additionally needs the d rows G(P,:) and d
!>      columns G(:,P) of the same vertex, which form the factors of its
!>      rank-d update.
!>
!> Thus we are at liberty to split the Green's function into
!>
!>             G = G_stale + X * Y^T,   X, Y of shape (Ndim, ~k)
!>
!> and pay Ndim**2 only once every k panel columns, i.e. every ~k/d accepted
!> updates. Traffic per accepted field update is then expected to fall from
!> ~2*Ndim**2 to ~d*Ndim*k + 2*d*Ndim**2/k.
!> The matrices X and Y are referred to as "panels". The scheme, including its
!> generalisation to vertices of rank d > 1, follows F. Sun and X. Y. Xu,
!> Phys. Rev. B 109, 235140 (2024); see the "Delayed (rank-k) updates" section
!> of the ALF documentation.
!>
!> The implementation is such that the Green's function is in its "factorised"
!> form only within the sequential vertex loop of a single time slice;
!> thus stabilisation, measurement and global moves always receive the fully
!> flushed G.
!>
!> By default delayed updates are off. The environment variable ALF_DELAY_K set
!> to the appropriate value enables it; see delay_resolve for more details.
!-------------------------------------------------------------------------------

module delayed_update_mod

   use operator_mod
   use runtime_error_mod
   use iso_fortran_env, only: error_unit

   implicit none

   private
   public :: delay_alloc, delay_dealloc, delay_resolve, delay_probe_pick
   public :: delay_assert_inactive, delay_open, delay_close
   public :: delay_block, delay_row, delay_col, delay_append, delay_flush
   public :: delay_wrap, delay_pending, delay_log

   ! Define the panels X and Y with their live column count, ncol, one set per
   ! flavour. The panels are only Ndim*(k+dmax) per flavour hence we can
   ! afford to allocate them once.
   complex (Kind=Kind(0.d0)), private, save, allocatable :: xp(:,:,:), yp(:,:,:)
   integer,                   private, save, allocatable :: ncol(:)

   integer, private, save :: kmax    = 0 ! Flush threshold, the k of the scheme
   integer, private, save :: panel_w = 0 ! Allocated panel width, kmax + dmax
   integer, private, save :: ndim_s  = 0 ! Ndim, recorded for the GR dummies
   integer, private, save :: nfl_s   = 0 ! N_FL, recorded for the GR dummies

   ! Record if a factored region is open
   logical, public, protected, save :: delay_active = .false.

   ! ALF_DELAY_K as parsed by delay_resolve, kept for delay_log: 0 disables the
   ! delay and a positive value fixes the depth, so named requests are negative.
   integer, private, parameter :: K_AUTO    = -1 ! "auto": measure the depth
   integer, private, parameter :: K_FORMULA = -2 ! "formula": k ~ sqrt(2*Ndim)
   integer, private, save      :: k_request = 0

   ! Constants on the depth "auto" may resolve to: a floor and ceiling.
   integer, private, parameter :: K_FLOOR   = 8
   integer, private, parameter :: K_CEILING = 256

   ! How the depth was chosen, for the info file: off, fixed, probe or formula.
   character (Len=16), public, protected, save :: delay_source = 'off'

   ! Depths the delay depth probe times.
   integer, private, parameter :: K_CAND(*) = [8, 16, 32, 64, 128, 256]
   integer, private, parameter :: N_CAND = size(K_CAND)

   ! Delay depth probe timing details.
   real (Kind=Kind(0.d0)), private, parameter :: PROBE_MIN_SECONDS = 5.d-3
   integer,                private, parameter :: PROBE_MAX_REPS    = 4096

   ! Sweep the probe three times to average over the background of the
   ! execution environment.
   integer, private, parameter :: PROBE_SWEEPS = 3

   ! Define a baseline margin in the comparison of the results from the
   ! delay depth probe.
   real (Kind=Kind(0.d0)), private, parameter :: PROBE_MARGIN = 1.05d0

   ! Complex identity "z-one"
   complex (Kind=Kind(0.d0)), private, parameter :: ZONE = (1.d0, 0.d0)

   ! Kernel selector for the flush / panel
   integer, private, parameter :: PROBE_FLUSH = 1 ! stand in for delay_flush
   integer, private, parameter :: PROBE_PANEL = 2 ! stand in for delay_row/col
   integer, private, parameter :: PROBE_IMMEDIATE = 3 ! stand in for Upgrade2's
                                                      ! immediate rank-d update

   ! Details for the delay_log
   real (Kind=Kind(0.d0)), private, save :: probe_cost(N_CAND) = -1.d0
   real (Kind=Kind(0.d0)), private, save :: probe_imm_cost = -1.d0
   real (Kind=Kind(0.d0)), private, save :: probe_seconds = 0.d0
   real (Kind=Kind(0.d0)), private, save :: probe_scratch_mb = 0.d0
   character (Len=64),     private, save :: k_request_text = '<unset>'

   ! GR dummies are explicit shape so GR(1,1,nf) can go to BLAS as a
   ! contiguous matrix; an assumed shape would not guarantee that.

   ! Argument names shared by the routines below:
   !    nf   - flavour
   !    d    - column pairs per update, Op%N_non_zero (or 1)
   !    P(d) - leading entries of Op%P, the rows and columns updated
   !    c    - local; live panel columns, ncol(nf)

contains

!-------------------------------------------------------------------------------
!> @brief
!> Resolve the delay depth for this run; 0 when the delayed update is disabled.
!>
!> @details
!> Reads ALF_DELAY_K, case-insensitively:
!>
!>    - unset, empty, "0", negative, unreadable or over 32 characters: delay
!>      off, with the reason in delay_log
!>    - a positive integer: used verbatim up to Ndim, capped there
!>    - "formula": delay_formula, k ~ sqrt(2*Ndim)
!>    - "auto": delay_probe times the flush and panel costs at this Ndim and
!>      rank dmax against the immediate update. It keeps the immediate update
!>      (k = 0) unless some depth beats it by PROBE_MARGIN, else takes the
!>      largest k within PROBE_MARGIN of the cheapest, falling back to
!>      delay_formula if the measurement is impossible or inconclusive
!>
!> "formula" and "auto" are clamped to [K_FLOOR, K_CEILING], and turn the delay
!> off when Ndim < K_FLOOR. Also sets delay_source and what delay_log reports.
!>
!> Call once and pass the result to delay_alloc: "auto" is a timing, so a
!> second call may answer differently. Under MPI only one rank should call it
!> and broadcast the result, as Wrapgr_delay_alloc does: delay_probe holds an
!> Ndim**2 scratch, and ranks probing at once would contend for the very memory
!> system they are measuring.
!>
!> A measured depth may differ between runs of one chain. This is harmless:
!> different k give the same chain up to rounding, which can only change
!> Metropolis decisions that are near-ties.
!>
!> The depth is deliberately not a simulation parameter, so that a
!> run-to-run choice does not show up as a difference between runs' parameters.
!-------------------------------------------------------------------------------

   integer function delay_resolve(Ndim, dmax)

      implicit none

      integer, intent(in) :: Ndim
      integer, intent(in) :: dmax  ! largest rank one update appends, for "auto"
      character (Len=32) :: text   ! ALF_DELAY_K as the environment gave it
      character (Len=32) :: word   ! the same, trimmed and lower-cased
      integer :: length, status    ! from get_environment_variable
      integer :: value             ! the depth, where the request was a number
      integer :: i

      k_request      = 0
      k_request_text = '<unset>'
      delay_source   = 'off'

      call get_environment_variable("ALF_DELAY_K", text, length, status)
      if (status == -1) then
         ! Longer than text holds: turn delays off, but say so.
         k_request_text = '<too long>'
      else if (status == 0 .and. length > 0) then
         word           = trim(adjustl(text(1:length)))
         k_request_text = word
         ! Capitalisation agnostic
         do i = 1, len_trim(word)
            if (word(i:i) >= 'A' .and. word(i:i) <= 'Z') &
            &  word(i:i) = achar(iachar(word(i:i)) + 32)
         enddo
         select case (word)
          case ('auto')
            k_request = K_AUTO
          case ('formula')
            k_request = K_FORMULA
          case default
            ! Unreadable or negative turns delays off, but say so.
            read (word, *, iostat=status) value
            if (status /= 0) then
               k_request_text = trim(k_request_text)//' (unreadable)'
            else if (value < 0) then
               k_request_text = trim(k_request_text)//' (negative)'
            else
               k_request = value
            endif
         end select
      endif

      select case (k_request)
       case (K_FORMULA)
         delay_resolve = delay_formula(Ndim)
         delay_source  = 'formula'
       case (K_AUTO)
         ! delay_probe sets delay_source itself, since it may have to fall back
         ! to the formula.
         delay_resolve = delay_probe(Ndim, dmax)
       case default
         delay_resolve = k_request
         delay_source  = 'fixed'
      end select
      ! Past Ndim the panels outgrow G and buy nothing, so a depth beyond it is
      ! capped there, and said so. This also keeps kmax + dmax well clear of
      ! overflow in delay_alloc.
      if (delay_resolve > Ndim) then
         write (word, '(i0)') Ndim
         k_request_text = trim(k_request_text)//' (capped at Ndim = '//trim(word)//')'
         delay_resolve  = Ndim
      endif

      ! Whichever path gave it, a zero depth is off -- except when the probe
      ! chose it, which is kept so that the run record says why.
      if (delay_resolve <= 0 .and. trim(delay_source) /= 'probe') delay_source = 'off'

      ! A named request resolves to zero on a matrix below the floor, or when
      ! the probe found the immediate update faster.
      if (delay_resolve <= 0 .and. k_request < 0) then
         if (Ndim < K_FLOOR) then
            write (word, '(i0)') K_FLOOR
            k_request_text = trim(k_request_text)//' (Ndim < '//trim(word)//')'
         else
            k_request_text = trim(k_request_text)//' (probe: immediate faster)'
         endif
      endif
   end function delay_resolve

!-------------------------------------------------------------------------------
!> @brief
!> Log what the delay is set to and the decision pathway.
!>
!> @details
!> Called once at setup, from main under the rank guard. Whenever the probe ran
!> its curve is printed, with the immediate update as k = 0 on the same scale,
!> including when the immediate update won and the delay is off.
!-------------------------------------------------------------------------------

   subroutine delay_log(unit)
      implicit none
      integer, intent(in) :: unit      ! where to write; main passes 6 (stdout)
      integer :: i
      real (Kind=Kind(0.d0)) :: lo          ! best cost, the curve's normaliser
      character (Len=13) :: mark, tag       ! marker for the depth in force
      logical :: timed    ! the probe produced a usable curve
      logical :: placed   ! the depth in force has a row in the table

      write (unit,'(a)')  ' Delayed update:'
      write (unit,'(2a)') '   ALF_DELAY_K            : ', trim(k_request_text)
      if (kmax == 0) then
         write (unit,'(a)') '   status                 : off (immediate)'
         ! The probe chose the immediate update: show the curve it lost to.
         if (trim(delay_source) /= 'probe' .or. probe_cost(1) < 0.d0) return
         write (unit,'(a,i0)') '   Ndim                   : ', ndim_s
         write (unit,'(2a)')   '   chosen by              : ', trim(delay_source)
      else
         write (unit,'(a,i0)') '   Ndim                   : ', ndim_s
         write (unit,'(a,i0)') '   depth k                : ', kmax
         write (unit,'(2a)')   '   chosen by              : ', trim(delay_source)
         write (unit,'(a,i0)') '   panel width (k + dmax) : ', panel_w
         write (unit,'(a,i0,a,i0,a)') '   validated range        : [', &
         & K_FLOOR, ', ', K_CEILING, ']'
      endif

      ! Say why auto fell back to the formula: a curve that is flat to within
      ! PROBE_MARGIN is a result, a probe that gave no reading is not.
      timed = any(probe_cost > 0.d0 .and. probe_cost < huge(1.d0)) .and. &
      &       probe_imm_cost > 0.d0 .and. probe_imm_cost < huge(1.d0)
      if (trim(delay_source) == 'formula' .and. k_request == K_AUTO) then
         if (timed) then
            write (unit,'(a)') '   note: the curve is flat; this is the formula'
         else
            write (unit,'(a)') '   WARN: the probe was refused; this is the formula'
         endif
      endif

      ! Still at its -1 default: the probe never ran, so there is no curve.
      if (probe_cost(1) < 0.d0) return

      write (unit,'(a,f6.3,a,f10.3,a)') '   probe cost             : ', &
      & probe_seconds, ' s, scratch ', probe_scratch_mb, ' MB'
      write (unit,'(a)') '        k   rel. cost   (1.00 = best depth)'

      lo = minval(probe_cost, mask=(probe_cost > 0.d0))

      ! Mark the depth in force; on a fallback it came from the formula.
      mark = '   <- chosen'
      if (trim(delay_source) == 'formula') mark = '   <- formula'

      ! The immediate update, per column like the candidates, as k = 0.
      tag = ''
      if (kmax == 0) tag = mark
      if (probe_imm_cost > 0.d0 .and. probe_imm_cost < huge(1.d0)) then
         write (unit,'(a,i5,f12.3,2a)') '   ', 0, probe_imm_cost/lo, &
         & '  (immediate)', trim(tag)
      else
         write (unit,'(a,i5,a)') '   ', 0, '       --  (immediate, not timed)'
      endif

      ! A formula depth need not be a candidate; it gets its own row, in order.
      placed = kmax == 0 .or. any(K_CAND == kmax)
      do i = 1, N_CAND
         if (.not. placed .and. kmax < K_CAND(i)) then
            write (unit,'(a,i5,2a)') '   ', kmax, '       --  (not timed)', trim(mark)
            placed = .true.
         endif
         tag = ''
         if (K_CAND(i) == kmax) tag = mark
         if (K_CAND(i) > ndim_s) then
            write (unit,'(a,i5,a)') '   ', K_CAND(i), '       --  (above Ndim)'
         else if (probe_cost(i) <= 0.d0) then
            write (unit,'(a,i5,a)') '   ', K_CAND(i), '       --  (not timed)'
         else
            write (unit,'(a,i5,f12.3,a)') '   ', K_CAND(i), probe_cost(i)/lo, &
            & trim(tag)
         endif
      enddo
      if (.not. placed) &
      & write (unit,'(a,i5,2a)') '   ', kmax, '       --  (not timed)', trim(mark)
   end subroutine delay_log

!-------------------------------------------------------------------------------
!> @brief
!> The closed-form depth: nint(sqrt(2*Ndim)), clamped; 0 when Ndim < K_FLOOR.
!>
!> @details
!> Minimises the traffic per accepted update, d*Ndim*k for the panel
!> reconstructions (2*d ZGEMVs against a half-full panel) plus 2*d*Ndim**2/k
!> for the flush; d cancels. It assumes both terms move bytes at the same
!> cost. In practice the panel is often cache-resident and the flush partly
!> compute-bound, so the true optimum tends to lie above this estimate.
!-------------------------------------------------------------------------------

   integer function delay_formula(Ndim)
      implicit none
      integer, intent(in) :: Ndim
      ! Below the floor the matrix is too small for the delay to pay. Above it
      ! the clamp keeps k within the range the delayed path has been exercised
      ! over, which for Ndim >= K_FLOOR is also never wider than the matrix.
      if (Ndim < K_FLOOR) then
         delay_formula = 0
      else
         delay_formula = min(K_CEILING, &
         &               max(K_FLOOR, nint(sqrt(2.d0*real(Ndim, Kind(0.d0))))))
      endif
   end function delay_formula

!-------------------------------------------------------------------------------
!> @brief
!> Pick the delay depth, or no delay, by timing both schemes at this Ndim.
!>
!> @details
!> Per accepted update of rank d the delayed scheme pays a flush,
!> ZGEMM('N','T',Ndim,Ndim,k), once every k/d updates, and 2*d panel ZGEMVs
!> against a panel that is half full on average:
!>
!>     cost(k) = t_gemm(k)*d/k + 2*d*t_gemv(k/2)
!>             = d * [ t_gemm(k)/k + 2*t_gemv(k/2) ]
!>
!> The immediate scheme instead pays one rank-d update of G, t_imm(d): ZGERU
!> for d = 1 and ZGEMM('N','T',Ndim,Ndim,d) otherwise, as Upgrade2 does. d
!> factors out between candidate depths but not against t_imm -- with d small an
!> immediate rank-d update is still about one pass over G -- so the probe is
!> timed at d = dmax, the largest rank of the model. t_gemm, t_gemv and t_imm
!> are probed at runtime.
!>
!> The model leaves out what only the delayed scheme pays -- the d x d block
!> per proposal, the panel wrap per vertex, the partial flush that closes a
!> slice -- so it flatters the delay, and delay_probe_pick breaks ties towards
!> the immediate update.
!>
!> Timed on an otherwise idle node (under MPI one rank probes, the rest wait).
!> Under full load the flush slows more than the cache-resident panel, so the
!> true optimum lies somewhat higher; one reason ties between depths go to the
!> larger k.
!>
!> Falls back to delay_formula when no candidate fits Ndim, the scratch cannot
!> be allocated, or the clock gives no reading; see delay_probe_pick for the
!> rest.
!-------------------------------------------------------------------------------

   integer function delay_probe(Ndim, dmax)

      implicit none

      integer, intent(in) :: Ndim
      integer, intent(in) :: dmax   ! largest rank one update appends

      ! Scratch the kernels run on: g stands in for the Green's function, xs
      ! and ys for the panels (and for the factors of an immediate update), v
      ! for a row of one panel and w for the rebuilt row / column.
      complex (Kind=Kind(0.d0)), allocatable :: g(:,:), xs(:,:), ys(:,:)
      complex (Kind=Kind(0.d0)), allocatable :: v(:), w(:)

      ! cost(i) is the modelled per-column cost at K_CAND(i), built from the
      ! flush time tg and the panel time tv of one reading, this. t_imm is the
      ! time of one immediate rank-dmax update.
      real (Kind=Kind(0.d0)) :: cost(N_CAND), tg, tv, this, t_imm

      ! k is the candidate depth, c the half occupancy the panel is timed at,
      ! kwide the widest candidate that fits Ndim, ncols the width allocated,
      ! stat the allocation status.
      integer :: i, k, c, kwide, ncols, stat, sweep

      ! Wall clock over the whole ladder, which delay_log reports.
      integer (Kind=8) :: wall0, wall1, wall_rate

      ! The fallback
      delay_source = 'formula'
      delay_probe  = delay_formula(Ndim)

      ! No candidate may exceed Ndim. One is enough: a single depth can still be
      ! set against the immediate update, even without a curve to choose from.
      kwide = 0
      do i = 1, N_CAND
         if (K_CAND(i) <= Ndim) kwide = K_CAND(i)
      enddo
      if (kwide < K_CAND(1)) return
      ncols = max(kwide, dmax)

      allocate (g(Ndim,Ndim), xs(Ndim,ncols), ys(Ndim,ncols), &
      & v(kwide), w(Ndim), stat=stat)
      if (stat /= 0) return

      ! Everything just allocated, at 16 bytes per complex entry.
      probe_scratch_mb = 16.d0*(real(Ndim, Kind(0.d0))**2 &
      &                  + 2.d0*real(Ndim, Kind(0.d0))*real(ncols, Kind(0.d0)) &
      &                  + real(kwide + Ndim, Kind(0.d0)))/1048576.d0

      call probe_fill(g,  Ndim*Ndim)
      call probe_fill(xs, Ndim*ncols)
      call probe_fill(ys, Ndim*ncols)
      call probe_fill(v,  kwide)
      call probe_fill(w,  Ndim)

      cost  = huge(1.d0)
      t_imm = huge(1.d0)

      call system_clock(wall0)

      ! The immediate update is timed in every sweep beside the candidates, so
      ! that both schemes see the same cache and background.
      do sweep = 1, PROBE_SWEEPS
         t_imm = min(t_imm, probe_time(PROBE_IMMEDIATE, g, xs, ys, v, w, Ndim, dmax))
         do i = 1, N_CAND
            k = K_CAND(i)
            if (k > kwide) cycle
            c  = k/2
            tg = probe_time(PROBE_FLUSH, g, xs, ys, v, w, Ndim, k)
            tv = probe_time(PROBE_PANEL, g, xs, ys, v, w, Ndim, c)
            this = tg/real(k, Kind(0.d0)) + 2.d0*tv
            ! Use minimum as "best case" scenario
            cost(i) = min(cost(i), this)
         enddo
      enddo

      call system_clock(wall1, wall_rate)

      if (wall_rate > 0) probe_seconds = &
      & real(wall1 - wall0, Kind(0.d0))/real(wall_rate, Kind(0.d0))
      probe_cost = cost
      ! Per column, on the same scale as probe_cost, for delay_log.
      probe_imm_cost = t_imm/real(dmax, Kind(0.d0))

      deallocate (g, xs, ys, v, w)

      delay_probe = delay_probe_pick(cost, t_imm, dmax, Ndim, delay_source)
   end function delay_probe

!-------------------------------------------------------------------------------
!> @brief
!> Choose the depth from a probed cost curve: the decision half of delay_probe.
!>
!> @details
!> Kept apart from the timing so that it can be tested on made-up curves.
!> cost(i) is the modelled per-column cost at K_CAND(i), with huge(1.d0) for a
!> candidate that was not timed; t_imm is the time of one immediate update of
!> rank d, which the best candidate costs d times over. In order:
!>
!>    - no usable reading, of the curve or of t_imm: delay_formula, with source
!>      'formula'
!>    - the immediate update within PROBE_MARGIN of the best candidate, or
!>      cheaper: 0, with source 'probe'. Ties go to the immediate update, the
!>      simpler path, and the one the cost model does not flatter.
!>    - a curve flat to within PROBE_MARGIN: delay_formula, with source
!>      'formula'; any depth then costs about the best, which beats t_imm.
!>    - otherwise the largest candidate within PROBE_MARGIN of the cheapest,
!>      with source 'probe': overshooting tends to give better performance on
!>      average.
!-------------------------------------------------------------------------------

   integer function delay_probe_pick(cost, t_imm, d, Ndim, source)
      implicit none
      real (Kind=Kind(0.d0)), intent(in)  :: cost(N_CAND)
      real (Kind=Kind(0.d0)), intent(in)  :: t_imm   ! one immediate update
      integer,                intent(in)  :: d       ! its rank
      integer,                intent(in)  :: Ndim
      character (Len=*),      intent(out) :: source
      real (Kind=Kind(0.d0)) :: lo, hi   ! bounds of the timed part of the curve

      ! The fallback
      source           = 'formula'
      delay_probe_pick = delay_formula(Ndim)

      ! Nothing timed, or a clock that gave no reading.
      if (count(cost < huge(1.d0)) == 0) return
      if (t_imm <= 0.d0 .or. t_imm >= huge(1.d0)) return

      ! Bounds of the curve; untimed candidates sit at "huge".
      lo = minval(cost, mask=(cost < huge(1.d0)))
      hi = maxval(cost, mask=(cost < huge(1.d0)))
      if (lo <= 0.d0) return

      ! Checked before the flat curve, which would otherwise turn the delay on
      ! without ever asking whether it pays.
      if (t_imm <= PROBE_MARGIN*real(d, Kind(0.d0))*lo) then
         delay_probe_pick = 0
         source           = 'probe'
         return
      endif

      ! When the curve is flat prefer the formula.
      if (hi < PROBE_MARGIN*lo) return

      delay_probe_pick = maxval(K_CAND, mask=(cost <= PROBE_MARGIN*lo))
      source           = 'probe'
   end function delay_probe_pick

!-------------------------------------------------------------------------------
!> @brief
!> Fill probe scratch with finite entries of order one.
!>
!> @details
!> Touches every page before timing starts and keeps NaNs and denormals, which
!> can run slower, out of the kernels. A closed form rather than the RNG, so
!> the probe leaves the Markov chain's random stream untouched.
!-------------------------------------------------------------------------------

   subroutine probe_fill(a, n)
      implicit none
      integer, intent(in) :: n
      complex (Kind=Kind(0.d0)), intent(out) :: a(*)
      integer :: i
      do i = 1, n
         a(i) = cmplx(sin(real(i, Kind(0.d0))), cos(real(3*i, Kind(0.d0))), &
         &            Kind(0.d0))
      enddo
   end subroutine probe_fill

!-------------------------------------------------------------------------------
!> @brief
!> Seconds for one call of a probe kernel, repeated until the clock resolves.
!>
!> @details
!> PROBE_FLUSH times one flush ZGEMM('N','T',Ndim,Ndim,n), PROBE_PANEL one
!> panel ZGEMV against n live columns, and PROBE_IMMEDIATE one immediate update
!> of rank n, as Upgrade2 makes it. A clock reporting no rate returns zero,
!> which delay_probe reads as a failed measurement.
!-------------------------------------------------------------------------------

   real (Kind=Kind(0.d0)) function probe_time(kernel, g, xs, ys, v, w, Ndim, n)
      implicit none
      ! kernel is PROBE_FLUSH, PROBE_PANEL or PROBE_IMMEDIATE; n is the depth k
      ! for the flush, the live column count for the panel and the rank for the
      ! immediate update. Every arm takes every array so that the calls in
      ! delay_probe read alike.
      integer, intent(in) :: kernel, Ndim, n
      complex (Kind=Kind(0.d0)), intent(inout) :: g(Ndim,Ndim), w(Ndim)
      complex (Kind=Kind(0.d0)), intent(in)    :: xs(Ndim,*), ys(Ndim,*), v(*)
      integer :: rep, reps   ! repeats, raised until the clock resolves
      integer (Kind=8) :: c0, c1, rate
      reps = 1
      do
         call system_clock(c0)
         do rep = 1, reps
            if (kernel == PROBE_FLUSH) then
               call ZGEMM('N', 'T', Ndim, Ndim, n, ZONE, xs, Ndim, &
               &          ys, Ndim, ZONE, g, Ndim)
            else if (kernel == PROBE_IMMEDIATE .and. n == 1) then
               call ZGERU(Ndim, Ndim, ZONE, xs, 1, ys, 1, g, Ndim)
            else if (kernel == PROBE_IMMEDIATE) then
               call ZGEMM('N', 'T', Ndim, Ndim, n, ZONE, xs, Ndim, &
               &          ys, Ndim, ZONE, g, Ndim)
            else
               call ZGEMV('N', Ndim, n, ZONE, xs, Ndim, v, 1, &
               &          ZONE, w, 1)
            endif
         enddo
         call system_clock(c1, rate)
         ! No rate means no measurement, at any repeat count: return the zero
         ! delay_probe reads as a failure rather than escalating to
         ! PROBE_MAX_REPS full-size kernels to be told the same thing.
         if (rate <= 0) then
            probe_time = 0.d0
            return
         endif
         probe_time = real(c1 - c0, Kind(0.d0))/real(rate, Kind(0.d0))
         if (probe_time >= PROBE_MIN_SECONDS .or. reps >= PROBE_MAX_REPS) exit
         reps = reps*4
      enddo
      probe_time = probe_time/real(reps, Kind(0.d0))
   end function probe_time

!-------------------------------------------------------------------------------
!> @brief
!> Live column count for one flavour, for tests and diagnostics.
!-------------------------------------------------------------------------------

   integer function delay_pending(nf)
      implicit none
      integer, intent(in) :: nf
      delay_pending = 0
      if (allocated(ncol)) delay_pending = ncol(nf)
   end function delay_pending

!-------------------------------------------------------------------------------
!> @brief
!> Stop if the factored region is open.
!>
!> @details
!> For callers that read, copy or wrap GR as a whole (Wrapgr_PlaceGR,
!> Wrapgr_Random_update): inside the region GR holds only G_stale, so they
!> would silently corrupt the chain.
!-------------------------------------------------------------------------------

   subroutine delay_assert_inactive(where)
      implicit none
      character(len=*), intent(in) :: where
      if (delay_active) then
         write(error_unit,*) 'delayed_update: ', trim(where), &
         & ' reached with a factored Green function open'
         Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
   end subroutine delay_assert_inactive

!-------------------------------------------------------------------------------
!> @brief
!> Allocate the panels, kmax + dmax columns wide; no-op when k <= 0.
!>
!> @param[in] dmax Most columns one update appends, maxval(Op_V%N_non_zero).
!> @param[in] k    The depth, from delay_resolve.
!>
!> @details
!> The extra dmax columns hold the append that crosses kmax before its flush.
!-------------------------------------------------------------------------------

   subroutine delay_alloc(Ndim, N_FL, dmax, k)
      implicit none
      integer, intent(in) :: Ndim, N_FL, dmax, k
      integer :: stat              ! from the panel allocate
      character (Len=256) :: msg   ! the runtime's reason, when it refuses

      ! Recorded even when the delay is off; see the note on the GR dummies.
      ndim_s = Ndim
      nfl_s  = N_FL

      ! Capped at Ndim as delay_resolve caps it, for callers that pass a depth
      ! directly: kmax + dmax must not overflow.
      kmax = min(max(k, 0), Ndim)
      if (kmax == 0) return ! No-op when delays are off

      panel_w = kmax + dmax
      ! The depth was asked for explicitly, so a refusal stops the run rather
      ! than quietly falling back to the immediate update; it happens at setup,
      ! before any work is lost.
      allocate (xp(Ndim, panel_w, N_FL), yp(Ndim, panel_w, N_FL), &
      &         stat=stat, errmsg=msg)
      if (stat /= 0) then
         write(error_unit,'(a,i0,a,i0,a,i0,a,f0.1,a)') &
         & 'delay_alloc: cannot allocate the panels for depth k = ', kmax, &
         & ' (Ndim = ', Ndim, ', N_FL = ', N_FL, ', ', &
         & 32.d0*real(Ndim, Kind(0.d0))*real(panel_w, Kind(0.d0)) &
         & *real(N_FL, Kind(0.d0))/1048576.d0, ' MB)'
         write(error_unit,'(2a)') 'delay_alloc: ALF_DELAY_K = ', trim(k_request_text)
         write(error_unit,'(2a)') 'delay_alloc: ', trim(msg)
         write(error_unit,'(a)')  'delay_alloc: lower ALF_DELAY_K, or unset it'
         Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      allocate (ncol(N_FL))
      ncol   = 0
      delay_active = .false.
   end subroutine delay_alloc

   subroutine delay_dealloc()
      implicit none
      if (allocated(xp)) deallocate (xp)
      if (allocated(yp)) deallocate (yp)
      if (allocated(ncol)) deallocate (ncol)
      ! Back to the state delay_alloc found, so that a second allocation cannot
      ! inherit the shape of the first.
      kmax         = 0
      panel_w      = 0
      ndim_s       = 0
      nfl_s        = 0
      delay_active = .false.
   end subroutine delay_dealloc

!-------------------------------------------------------------------------------
!> @brief
!> Open a factored region.
!-------------------------------------------------------------------------------

   subroutine delay_open()
      implicit none
      if (kmax == 0) return ! No-op when the delay is disabled.
      ncol   = 0
      delay_active = .true.
   end subroutine delay_open

!-------------------------------------------------------------------------------
!> @brief
!> Flush every flavour into GR and close the region.
!-------------------------------------------------------------------------------

   subroutine delay_close(GR)
      implicit none
      complex (kind=kind(0.d0)), intent(inout) :: GR(ndim_s, ndim_s, nfl_s)
      integer :: nf
      if (.not. delay_active) return
      do nf = 1, nfl_s
         call delay_flush(nf, GR)
      enddo
      delay_active = .false.
   end subroutine delay_close

!-------------------------------------------------------------------------------
!> @brief
!> Apply the pending columns of one flavour to GR: G_stale += X*Y^T.
!-------------------------------------------------------------------------------

   subroutine delay_flush(nf, GR)
      implicit none
      integer, intent(in) :: nf
      complex (kind=kind(0.d0)), intent(inout) :: GR(ndim_s, ndim_s, nfl_s)
      if (kmax == 0) return
      if (ncol(nf) == 0) return
      call ZGEMM('N', 'T', ndim_s, ndim_s, ncol(nf), ZONE, xp(1,1,nf), &
      &          ndim_s, yp(1,1,nf), ndim_s, ZONE, GR(1,1,nf), ndim_s)
      ncol(nf) = 0
   end subroutine delay_flush

!-------------------------------------------------------------------------------
!> @brief
!> The d x d block of the current Green's function on the operator's support.
!>
!> @details
!> blk(n,m) = G(P(n), P(m)) = G_stale(P(n),P(m)) + sum_c X(P(n),c)*Y(P(m),c).
!>
!> O(d**2 * c), and paid on every proposal including the rejected ones -- this
!> is the cost the delay adds, and the reason a very large k stops paying.
!-------------------------------------------------------------------------------

   subroutine delay_block(nf, GR, P, d, blk, ldb)
      implicit none
      integer, intent(in) :: nf, d, P(d)
      integer, intent(in) :: ldb ! Leading dimension of blk
      complex (kind=kind(0.d0)), intent(in)  :: GR(ndim_s, ndim_s, nfl_s)
      complex (kind=kind(0.d0)), intent(out) :: blk(ldb,*)
      integer :: n, m   ! row and column within the support
      integer :: c
      c = ncol(nf)
      ! A sum over an empty panel is zero, so this covers c = 0 as it stands.
      do m = 1, d
         do n = 1, d
            blk(n,m) = GR(P(n), P(m), nf) &
            &        + sum(xp(P(n),1:c,nf)*yp(P(m),1:c,nf))
         enddo
      enddo
   end subroutine delay_block

!-------------------------------------------------------------------------------
!> @brief
!> d rows of the current Green's function: rows(i,l) = G(P(l), i).
!>
!> @details
!> Stale row plus Y * X(P(l),:)^T, one ZGEMV each. The stale part is a strided
!> read of a column-major matrix -- Ndim cache lines for Ndim*16 bytes -- which
!> is exactly the traffic the delay exists to amortise.
!-------------------------------------------------------------------------------

   subroutine delay_row(nf, GR, P, d, rows)
      implicit none
      integer, intent(in) :: nf, d, P(d)
      complex (kind=kind(0.d0)), intent(in)  :: GR(ndim_s, ndim_s, nfl_s)
      complex (kind=kind(0.d0)), intent(out) :: rows(ndim_s, d)
      ! One row of X gathered off the panel, the ZGEMV's coefficient vector.
      complex (kind=kind(0.d0)) :: tmp(max(ncol(nf),1))
      integer :: i      ! column of G, i.e. position along the row
      integer :: l, c   ! l indexes the support, as in rows(:,l) = G(P(l),:)
      c = ncol(nf)
      do l = 1, d
         do i = 1, ndim_s
            rows(i,l) = GR(P(l), i, nf)
         enddo
         if (c > 0) then
            tmp(1:c) = xp(P(l), 1:c, nf)
            call ZGEMV('N', ndim_s, c, ZONE, yp(1,1,nf), ndim_s, tmp, 1, &
            &          ZONE, rows(1,l), 1)
         endif
      enddo
   end subroutine delay_row

!-------------------------------------------------------------------------------
!> @brief
!> d columns of the current Green's function: cols(i,l) = G(i, P(l)).
!-------------------------------------------------------------------------------

   subroutine delay_col(nf, GR, P, d, cols)
      implicit none
      integer, intent(in) :: nf, d, P(d)
      complex (kind=kind(0.d0)), intent(in)  :: GR(ndim_s, ndim_s, nfl_s)
      complex (kind=kind(0.d0)), intent(out) :: cols(ndim_s, d)
      ! One row of Y gathered off the panel, the ZGEMV's coefficient vector.
      complex (kind=kind(0.d0)) :: tmp(max(ncol(nf),1))
      integer :: l, c   ! l indexes the support, as in cols(:,l) = G(:,P(l))
      c = ncol(nf)
      do l = 1, d
         call ZCOPY(ndim_s, GR(1, P(l), nf), 1, cols(1,l), 1)
         if (c > 0) then
            tmp(1:c) = yp(P(l), 1:c, nf)
            call ZGEMV('N', ndim_s, c, ZONE, xp(1,1,nf), ndim_s, tmp, 1, &
            &          ZONE, cols(1,l), 1)
         endif
      enddo
   end subroutine delay_col

!-------------------------------------------------------------------------------
!> @brief
!> Apply a rank-d update by appending d column pairs; flush when full.
!>
!> @details
!> The caller's update is G += alpha * xcols * ycols^T. The coefficient is
!> folded into the X column here so the flush stays a plain X*Y^T.
!-------------------------------------------------------------------------------

   subroutine delay_append(nf, alpha, xcols, ycols, d, GR)
      implicit none
      integer, intent(in) :: nf, d
      complex (kind=kind(0.d0)), intent(in)    :: alpha   ! update coefficient
      ! The two factors of the caller's rank-d update, one column pair per rank.
      complex (kind=kind(0.d0)), intent(in)    :: xcols(ndim_s, d)
      complex (kind=kind(0.d0)), intent(in)    :: ycols(ndim_s, d)
      complex (kind=kind(0.d0)), intent(inout) :: GR(ndim_s, ndim_s, nfl_s)
      integer :: l      ! which of the d column pairs is being appended
      integer :: c      ! columns already live, so c+l is where l lands
      c = ncol(nf)
      do l = 1, d
         call ZCOPY(ndim_s, xcols(1,l), 1, xp(1, c+l, nf), 1)
         call ZSCAL(ndim_s, alpha, xp(1, c+l, nf), 1)
         call ZCOPY(ndim_s, ycols(1,l), 1, yp(1, c+l, nf), 1)
      enddo
      ncol(nf) = c + d
      if (ncol(nf) >= kmax) call delay_flush(nf, GR)
   end subroutine delay_append

!-------------------------------------------------------------------------------
!> @brief
!> Conjugate the panels with a vertex operator, mirroring Op_Wrapup/Op_Wrapdo.
!>
!> @details
!> A thin pass-through to Op_Wrap_panels, which lives in Operator_mod beside the
!> routines it mirrors so that an edit to one is visible from the other.
!-------------------------------------------------------------------------------

   subroutine delay_wrap(nf, Op, HS_Field, N_Type, nt, updo)
      implicit none
      ! Op, HS_Field, N_Type and nt are the vertex, its field, the wrap variant
      ! and the time slice, exactly as Op_Wrapup and Op_Wrapdo take them; updo
      ! picks which of the two this call mirrors.
      integer, intent(in) :: nf, N_Type, nt
      Type (Operator), intent(in) :: Op
      complex (kind=kind(0.d0)), intent(in) :: HS_Field
      character(len=1), intent(in) :: updo
      if (.not. delay_active) return
      if (ncol(nf) == 0) return
      call Op_Wrap_panels(xp(1,1,nf), yp(1,1,nf), Op, HS_Field, ndim_s, &
      &                ncol(nf), N_Type, nt, updo)
   end subroutine delay_wrap

end module delayed_update_mod
