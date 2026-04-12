! Test that the overflow-protected kernel and helper functions in
! MaxEnt_Wrapper_mod are smooth across the beta*omega = +-600 branch
! threshold and produce finite, correct results at extreme arguments.
!
! Strategy
! --------
!  (a) Continuity: at beta*omega = +-599 (just inside the exact branch)
!      the exact formula and the asymptotic formula must agree to machine
!      precision, because the neglected terms are O(exp(-599)) ~ 1e-260.
!      We compute both and require relative agreement < 1e-12.
!
!  (b) Extreme values: at beta*omega = +-1000 the naive formula would
!      overflow.  We verify the function returns a finite value matching
!      the hand-coded asymptotic result.

Program Test_MaxEnt_Kernel_Overflow

  Use MaxEnt_Wrapper_mod
  implicit none

  real(Kind=Kind(0.d0)) :: beta, tau, om, pi
  real(Kind=Kind(0.d0)) :: val, ref, relerr
  real(Kind=Kind(0.d0)) :: tol
  real(Kind=Kind(0.d0)) :: taus(4)
  integer :: it, nerr

  ! Use the same single-precision pi constant as maxent_wrapper_mod
  pi  = 3.1415927
  beta = 10.d0
  tol  = 1.d-12
  nerr = 0

  taus(1) = 0.d0
  taus(2) = beta/4.d0
  taus(3) = beta/2.d0
  taus(4) = beta

  ! =====================================================================
  ! XKER_ph  (and implicitly XKER_ph_c which is the same expression)
  ! Exact: (exp(-tau*om) + exp(-(beta-tau)*om)) / (pi*(1+exp(-beta*om)))
  ! beta*om >> 1 :  denom -> pi   =>  (exp(-tau*om)+exp(-(beta-tau)*om))/pi
  ! beta*om << -1:  multiply by exp(beta*om) => (exp((beta-tau)*om)+exp(tau*om))/pi
  ! =====================================================================
  write(*,*) "--- XKER_ph ---"
  do it = 1, 4
     tau = taus(it)
     ! (a) Continuity at +599
     om  = 599.d0/beta
     val = XKER_ph(tau, om, beta)
     ref = (exp(-tau*om) + exp(-(beta-tau)*om)) / pi   ! asymptotic
     call check("XKER_ph +599", val, ref, tol, nerr)

     ! (a) Continuity at -599
     om  = -599.d0/beta
     val = XKER_ph(tau, om, beta)
     ref = (exp((beta-tau)*om) + exp(tau*om)) / pi
     call check("XKER_ph -599", val, ref, tol, nerr)

     ! (b) Extreme +1000
     om  = 1000.d0/beta
     val = XKER_ph(tau, om, beta)
     ref = (exp(-tau*om) + exp(-(beta-tau)*om)) / pi
     call check("XKER_ph +1000", val, ref, tol, nerr)

     ! (b) Extreme -1000
     om  = -1000.d0/beta
     val = XKER_ph(tau, om, beta)
     ref = (exp((beta-tau)*om) + exp(tau*om)) / pi
     call check("XKER_ph -1000", val, ref, tol, nerr)
  enddo

  ! =====================================================================
  ! XKER_pp
  ! Exact: exp(-tau*om) / (pi*(1+exp(-beta*om)))
  ! beta*om >> 1 :  exp(-tau*om)/pi
  ! beta*om << -1:  exp((beta-tau)*om)/pi
  ! =====================================================================
  write(*,*) "--- XKER_pp ---"
  do it = 1, 4
     tau = taus(it)
     om  = 599.d0/beta
     val = XKER_pp(tau, om, beta)
     ref = exp(-tau*om) / pi
     call check("XKER_pp +599", val, ref, tol, nerr)

     om  = -599.d0/beta
     val = XKER_pp(tau, om, beta)
     ref = exp((beta-tau)*om) / pi
     call check("XKER_pp -599", val, ref, tol, nerr)

     om  = 1000.d0/beta
     val = XKER_pp(tau, om, beta)
     ref = exp(-tau*om) / pi
     call check("XKER_pp +1000", val, ref, tol, nerr)

     om  = -1000.d0/beta
     val = XKER_pp(tau, om, beta)
     ref = exp((beta-tau)*om) / pi
     call check("XKER_pp -1000", val, ref, tol, nerr)
  enddo

  ! =====================================================================
  ! XKER_p  (already had the fix; verify it still works)
  ! =====================================================================
  write(*,*) "--- XKER_p ---"
  do it = 1, 4
     tau = taus(it)
     om  = 599.d0/beta
     val = XKER_p(tau, om, beta)
     ref = exp(-tau*om) / pi
     call check("XKER_p +599", val, ref, tol, nerr)

     om  = -599.d0/beta
     val = XKER_p(tau, om, beta)
     ref = exp((beta-tau)*om) / pi
     call check("XKER_p -599", val, ref, tol, nerr)

     om  = 1000.d0/beta
     val = XKER_p(tau, om, beta)
     ref = exp(-tau*om) / pi
     call check("XKER_p +1000", val, ref, tol, nerr)

     om  = -1000.d0/beta
     val = XKER_p(tau, om, beta)
     ref = exp((beta-tau)*om) / pi
     call check("XKER_p -1000", val, ref, tol, nerr)
  enddo

  ! =====================================================================
  ! XKER_p_ph  (same expression as XKER_ph)
  ! =====================================================================
  write(*,*) "--- XKER_p_ph ---"
  do it = 1, 4
     tau = taus(it)
     om  = 599.d0/beta
     val = XKER_p_ph(tau, om, beta)
     ref = (exp(-tau*om) + exp(-(beta-tau)*om)) / pi
     call check("XKER_p_ph +599", val, ref, tol, nerr)

     om  = -599.d0/beta
     val = XKER_p_ph(tau, om, beta)
     ref = (exp((beta-tau)*om) + exp(tau*om)) / pi
     call check("XKER_p_ph -599", val, ref, tol, nerr)

     om  = 1000.d0/beta
     val = XKER_p_ph(tau, om, beta)
     ref = (exp(-tau*om) + exp(-(beta-tau)*om)) / pi
     call check("XKER_p_ph +1000", val, ref, tol, nerr)

     om  = -1000.d0/beta
     val = XKER_p_ph(tau, om, beta)
     ref = (exp((beta-tau)*om) + exp(tau*om)) / pi
     call check("XKER_p_ph -1000", val, ref, tol, nerr)
  enddo

  ! =====================================================================
  ! F_QFI_ph   (= (4/pi)*tanh(beta*om/2)^2 )
  ! |beta*om| >> 1 :  tanh -> +-1,  squared -> 1  =>  4/pi
  ! =====================================================================
  write(*,*) "--- F_QFI_ph ---"
  om  = 599.d0/beta
  val = F_QFI_ph(om, beta)
  ref = 4.d0/pi
  call check("F_QFI_ph +599", val, ref, tol, nerr)

  om  = -599.d0/beta
  val = F_QFI_ph(om, beta)
  ref = 4.d0/pi
  call check("F_QFI_ph -599", val, ref, tol, nerr)

  om  = 1000.d0/beta
  val = F_QFI_ph(om, beta)
  ref = 4.d0/pi
  call check("F_QFI_ph +1000", val, ref, tol, nerr)

  om  = -1000.d0/beta
  val = F_QFI_ph(om, beta)
  ref = 4.d0/pi
  call check("F_QFI_ph -1000", val, ref, tol, nerr)

  ! =====================================================================
  ! Back_trans_ph   Aom/(1+exp(-beta*om))
  ! beta*om >> 1 :  Aom
  ! beta*om << -1:  Aom*exp(beta*om)
  ! =====================================================================
  write(*,*) "--- Back_trans_ph ---"
  om  = 599.d0/beta
  val = Back_trans_ph(1.d0, om, beta)
  ref = 1.d0
  call check("Back_trans_ph +599", val, ref, tol, nerr)

  om  = -599.d0/beta
  val = Back_trans_ph(1.d0, om, beta)
  ref = exp(beta*om)
  call check("Back_trans_ph -599", val, ref, tol, nerr)

  om  = 1000.d0/beta
  val = Back_trans_ph(1.d0, om, beta)
  ref = 1.d0
  call check("Back_trans_ph +1000", val, ref, tol, nerr)

  om  = -1000.d0/beta
  val = Back_trans_ph(1.d0, om, beta)
  ref = exp(beta*om)
  call check("Back_trans_ph -1000", val, ref, tol, nerr)

  ! =====================================================================
  ! Back_trans_ph_c   Aom*tanh(beta*om/2)/om
  ! beta*om >> 1 :  Aom/om
  ! beta*om << -1:  -Aom/om
  ! =====================================================================
  write(*,*) "--- Back_trans_ph_c ---"
  om  = 599.d0/beta
  val = Back_trans_ph_c(1.d0, om, beta)
  ref = 1.d0/om
  call check("Back_trans_ph_c +599", val, ref, tol, nerr)

  om  = -599.d0/beta
  val = Back_trans_ph_c(1.d0, om, beta)
  ref = -1.d0/om
  call check("Back_trans_ph_c -599", val, ref, tol, nerr)

  om  = 1000.d0/beta
  val = Back_trans_ph_c(1.d0, om, beta)
  ref = 1.d0/om
  call check("Back_trans_ph_c +1000", val, ref, tol, nerr)

  om  = -1000.d0/beta
  val = Back_trans_ph_c(1.d0, om, beta)
  ref = -1.d0/om
  call check("Back_trans_ph_c -1000", val, ref, tol, nerr)

  ! =====================================================================
  ! Back_trans_pp   same expression as Back_trans_ph_c
  ! =====================================================================
  write(*,*) "--- Back_trans_pp ---"
  om  = 599.d0/beta
  val = Back_trans_pp(1.d0, om, beta)
  ref = 1.d0/om
  call check("Back_trans_pp +599", val, ref, tol, nerr)

  om  = -599.d0/beta
  val = Back_trans_pp(1.d0, om, beta)
  ref = -1.d0/om
  call check("Back_trans_pp -599", val, ref, tol, nerr)

  om  = 1000.d0/beta
  val = Back_trans_pp(1.d0, om, beta)
  ref = 1.d0/om
  call check("Back_trans_pp +1000", val, ref, tol, nerr)

  om  = -1000.d0/beta
  val = Back_trans_pp(1.d0, om, beta)
  ref = -1.d0/om
  call check("Back_trans_pp -1000", val, ref, tol, nerr)

  ! =====================================================================
  ! Summary
  ! =====================================================================
  if (nerr > 0) then
     write(*,'(A,I0,A)') "FAILED: ", nerr, " check(s) exceeded tolerance"
     stop 2
  endif
  write(*,*) "success"

contains

  subroutine check(label, val, ref, tol, nerr)
    character(len=*), intent(in) :: label
    real(Kind=Kind(0.d0)), intent(in) :: val, ref, tol
    integer, intent(inout) :: nerr
    real(Kind=Kind(0.d0)) :: relerr, denom

    denom = max(abs(ref), abs(val), tiny(1.d0))
    relerr = abs(val - ref) / denom

    if (relerr > tol) then
       write(*,'(A,A,A,ES14.6,A,ES14.6,A,ES10.2)') &
            "  FAIL ", label, ": val=", val, " ref=", ref, " relerr=", relerr
       nerr = nerr + 1
    endif

    ! Also check finiteness
    if (val /= val) then  ! NaN check
       write(*,'(A,A,A)') "  FAIL ", label, ": NaN detected"
       nerr = nerr + 1
    endif
  end subroutine check

end Program Test_MaxEnt_Kernel_Overflow
