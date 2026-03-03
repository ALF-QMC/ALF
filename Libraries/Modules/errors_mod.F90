!  Copyright (C) 2018-2026 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


     MODULE ERRORS
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Comprehensive statistical error analysis toolkit for Monte Carlo simulations.
!
!> @details
!> This module provides a complete suite of error analysis methods for quantum
!> Monte Carlo data, including:
!>
!> **Standard Error Analysis:**
!> - ERRCALC: Standard deviation and error of the mean
!> - Support for real and complex observables
!>
!> **Jackknife Methods:**
!> - ERRCALCJ: Jackknife error estimation (leave-one-out resampling)
!> - Variants with rebinning for correlated data
!> - JS variants: Sign-corrected observables for QMC with sign problem (O = <s·O>/<s>)
!> - Function evaluation on jackknife samples
!>
!> **Bootstrap Resampling:**
!> - Bootstrap: Standard bootstrap error estimation
!> - Bootstrap_fluc: Bootstrap for fluctuations and covariances
!>
!> **Autocorrelation Analysis:**
!> - AUTO_COR: Compute autocorrelation functions and integrated autocorrelation times
!> - Essential for understanding Monte Carlo sampling efficiency
!>
!> **Covariance Calculations:**
!> - COV: Covariance matrices with various binning strategies
!> - Support for complex observables and weighted samples
!>
!> **Numerical Integration:**
!> - INTERGRATE: Trapezoid rule integration for QMC data
!> - INTERGRATE_F: Integration with functional form
!>
!> **Curve Fitting:**
!> - FIT: Linear and polynomial fitting routines
!>
!> @note
!> Jackknife is generally preferred over bootstrap for QMC data due to better
!> bias properties. Rebinning is crucial for properly accounting for autocorrelations.
!>
!> @warning
!> For correlated data (typical in QMC), use rebinned variants or account for
!> autocorrelation time to obtain correct error estimates.
!>
!> @see
!> For theoretical background:
!> - Jackknife: Efron (1982), "The Jackknife, the Bootstrap and Other Resampling Plans"
!> - Autocorrelation: Sokal (1997), "Monte Carlo Methods in Statistical Mechanics"
!--------------------------------------------------------------------
       use runtime_error_mod
       Use MyMats
       Use Random_Wrap
       use iso_fortran_env, only: output_unit, error_unit

       !> Generic interface for standard error calculation (mean and standard deviation)
       !> Supports real (ERRCALC) and complex (ERRCALC_C) data
       INTERFACE ERRCALC
          MODULE PROCEDURE ERRCALC, ERRCALC_C
       END INTERFACE
       
       !> Generic interface for jackknife error estimation
       !> Variants: simple (J), with rebinning (J_REBIN), for ratios (JS), and with functions (JS_F)
       !> Each variant has real and complex versions
       INTERFACE ERRCALCJ
          MODULE PROCEDURE ERRCALC_J,    ERRCALC_J_REBIN,    ERRCALC_JS, ERRCALC_JS_REBIN, &
               &           ERRCALC_J_C,  ERRCALC_J_C_REBIN,  ERRCALC_JS_C, ERRCALC_JS_C_REBIN, &
               &           ERRCALC_JS_F, ERRCALC_JS_REBIN_F, ERRCALC_JS_C_F, ERRCALC_JS_C_REBIN_F
       END INTERFACE
       
       !> Generic interface for covariance calculations
       !> Supports various binning strategies and background subtraction
       INTERFACE COV
          MODULE PROCEDURE COVJ, COVJS, COVJS_C, COVJS_C_REBIN, COVJS_C_BG, COVJS_C_REBIN_BG, &
               &           COVJS_C_BG_Weights,  COVJS_C_BG_Rebin_Weights
       END INTERFACE COV
       
       !> Interface for covariance error estimation
       INTERFACE COV_ERR
          MODULE PROCEDURE COV_ERR
       END INTERFACE
       
       !> Interface for numerical integration with functional form
       INTERFACE INTERGRATE_F
          MODULE PROCEDURE INTER_F
       END INTERFACE
       
       !> Interface for trapezoid rule integration of QMC data
       INTERFACE INTERGRATE
          MODULE PROCEDURE INTER_QMC
       END INTERFACE
       
       !> Interface for curve fitting (linear and polynomial)
       INTERFACE FIT
          MODULE PROCEDURE FIT
       END INTERFACE
       
       !> Interface for autocorrelation function calculation
       INTERFACE AUTO_COR
          MODULE PROCEDURE  AUTO_COR
       END INTERFACE
       
       !> Interface for bootstrap resampling error estimation
       INTERFACE Bootstrap
          MODULE PROCEDURE  Bootstrap
       END INTERFACE
       
       !> Interface for bootstrap estimation of fluctuations and covariances
       INTERFACE Bootstrap_fluc
          MODULE PROCEDURE  BootstrapC_fluc
       END INTERFACE
       
       !> Abstract interface for real-valued functions
       !> Used in error propagation calculations where a function of observables is analyzed
       abstract interface
           function func_r (Z)
               real (Kind=Kind(0.d0)) :: func_r
               real (Kind=Kind(0.d0)), allocatable, intent (in) :: Z(:)
           end function func_r
       end interface
       
       !> Abstract interface for complex-valued functions
       !> Used in error propagation calculations for complex observables
       abstract interface
           function func_c (X)
               COMPLEX (Kind=Kind(0.d0)) :: func_c
               COMPLEX (Kind=Kind(0.d0)), allocatable, intent (in) :: X(:)
           end function func_c
       end interface

       CONTAINS
!***********
         SUBROUTINE ERRCALC(EN,XM,XERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Computes mean and standard error of the mean for real data.
!
!> @details
!> Calculates the sample mean and standard error using:
!> - Mean: μ = Σx_i / N
!> - Variance: σ² = Σ(x_i - μ)² / N
!> - Standard error: SE = σ / √N
!>
!> This assumes independent samples. For correlated Monte Carlo data,
!> use jackknife methods or account for autocorrelation.
!
!> @param[in] EN Input data vector (bins or samples)
!> @param[out] XM Sample mean
!> @param[out] XERR Standard error of the mean
!
!> @see ERRCALC_C (complex version), ERRCALC_J (jackknife version)
!--------------------------------------------------------------------
           IMPLICIT NONE
           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN      ! Input data bins
           REAL (Kind=Kind(0.d0))               ::  XM      ! Mean value
           REAL (Kind=Kind(0.d0))               ::  XERR    ! Standard error of mean
           REAL (Kind=Kind(0.d0))               ::  XSQ     ! Variance
           INTEGER                     ::  NP, NT

           NP = SIZE(EN)  ! Number of bins

           ! Compute mean: μ = (1/N) Σ x_i
           XM  = 0.D0
           DO NT = 1,NP
              XM  = XM  + EN(NT)
           ENDDO
           XM    = XM /DBLE(NP)
           
           ! Compute variance: σ² = (1/N) Σ(x_i - μ)²
           XSQ = 0.D0
           DO NT = 1,NP
              XSQ = XSQ + (EN(NT)-XM)**2
           ENDDO
           XSQ   = XSQ/DBLE(NP)
           
           ! Standard error of mean: SE = σ/√N = √(σ²/N)
           XERR  = XSQ/DBLE(NP)
           IF (XERR.GT.0.D0) THEN
              XERR = SQRT(XERR)
           ELSE
              XERR = 0.D0
           ENDIF

           RETURN
         END SUBROUTINE ERRCALC

         SUBROUTINE ERRCALC_C(EN,ZM,ZERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Computes mean and standard error of the mean for complex data.
!
!> @details
!> Treats real and imaginary parts independently, calculating the error
!> for each component separately using the standard error formula.
!> The complex error is then constructed as: ZERR = XERR_real + i·XERR_imag
!
!> @param[in] EN Input complex data vector
!> @param[out] ZM Complex sample mean
!> @param[out] ZERR Complex standard error (real and imaginary parts independent)
!
!> @see ERRCALC (real version)
!--------------------------------------------------------------------
           IMPLICIT NONE
           Complex (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           Complex (Kind=Kind(0.d0))               ::  ZM, ZERR
           INTEGER                     ::  NP, NT

           ! Local variables
           Real (Kind=Kind(0.d0)), dimension(:), allocatable :: Rhelp
           real (Kind=Kind(0.d0)) :: XM, XERR

           NP = SIZE(EN)
           Allocate (Rhelp(NP))

           ! Process real part
           do nt = 1,np
              Rhelp(nt) = dble(en(nt))
           enddo
           call errcalc(Rhelp, xm, xerr)
           zm   =  cmplx(xm  , 0.d0, kind(0.D0))
           Zerr =  cmplx(xerr, 0.d0, kind(0.D0))

           ! Process imaginary part
           do nt = 1,np
              Rhelp(nt) = aimag(en(nt))
           enddo
           call errcalc(Rhelp, xm, xerr)
           zm   =  zm   + cmplx( 0.d0, xm, kind(0.D0)   )
           Zerr =  Zerr + cmplx( 0.d0, xerr, kind(0.D0) )

           RETURN
         END SUBROUTINE ERRCALC_C

         SUBROUTINE ERRCALC_J(EN,XM,XERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error estimation for real data.
!
!> @details
!> Implements the jackknife (leave-one-out) resampling method:
!> 1. For each bin i, compute mean of all bins except i:
!>    x̄_i = (Σx_j - x_i) / (N-1) for j≠i
!> 2. Calculate variance of these N jackknife estimates
!> 3. Scale by N to obtain error estimate: σ_jack = √(N·Var(x̄_i))
!>
!> The jackknife provides better bias correction than standard error,
!> especially for small samples and derived quantities.
!
!> @param[in] EN Input bins (typically Monte Carlo averages)
!> @param[out] XM Jackknife mean estimate
!> @param[out] XERR Jackknife error estimate
!
!> @note
!> Assumes bins are approximately independent. For correlated data,
!> use ERRCALC_J_REBIN with appropriate bin size.
!
!> @see ERRCALC_J_REBIN, ERRCALC_JS (for ratios)
!--------------------------------------------------------------------
           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, X, Xhelp
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER     :: NP, N, N1

           NP = SIZE(EN)
           ALLOCATE (EN1(NP))

           ! Compute sum of all bins (for efficiency)
           Xhelp = 0.D0
           DO N1 = 1,NP
              Xhelp = Xhelp + EN(N1)
           ENDDO

           ! Build N jackknife estimates (leave-one-out means)
           DO N = 1,NP
              X = Xhelp - EN(N)              ! Sum excluding bin N
              EN1(N) = X / DBLE(NP -1)       ! Mean excluding bin N
           ENDDO
           
           ! Calculate variance of jackknife estimates
           CALL ERRCALC(EN1,XM,XERR)
           
           ! Scale by N for jackknife error formula
           XERR = XERR*DBLE(NP)
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_J

         SUBROUTINE ERRCALC_J_C(EN,ZM,ZERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error estimation for complex data.
!
!> @details
!> Applies the jackknife method to complex observables using leave-one-out
!> resampling. Calculates error using the standard jackknife formula:
!> σ_jack = √(N·Var(z̄_i)) where z̄_i are jackknife estimates.
!
!> @param[in] EN Input complex bins
!> @param[out] ZM Jackknife mean estimate  
!> @param[out] ZERR Jackknife error estimate (complex)
!
!> @see ERRCALC_J (real version), ERRCALC_J_C_REBIN (with rebinning)
!--------------------------------------------------------------------
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           COMPLEX (Kind=Kind(0.d0))               ::  ZM, ZERR, Z, Zhelp
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER     :: NP, N, N1

           NP = SIZE(EN)
           ALLOCATE (EN1(NP))

           ! Compute complex sum of all bins
           Zhelp = CMPLX(0.D0, 0.D0, kind(0.D0))
           DO N1 = 1,NP
              Zhelp = Zhelp + EN(N1)
           ENDDO

           ! Build N jackknife estimates (leave-one-out)
           DO N = 1,NP
              Z =  Zhelp - EN(N)             ! Sum excluding bin N
              EN1(N) = Z / DBLE(NP -1)       ! Mean excluding bin N
           ENDDO
           CALL ERRCALC(EN1,ZM,ZERR)
           ZERR = ZERR*DBLE(NP)              ! Scale by N for jackknife formula
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_J_C

!************
         SUBROUTINE ERRCALC_J_C_REBIN(EN,ZM,ZERR,NREBIN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error estimation with rebinning for complex data.
!
!> @details
!> Combines rebinning with jackknife resampling to properly account for
!> autocorrelations in Monte Carlo data. Workflow:
!> 1. Rebin data: group consecutive bins into larger bins of size NREBIN
!> 2. Apply jackknife method to rebinned data
!>
!> Rebinning reduces autocorrelation by averaging over correlated samples.
!> Choose NREBIN ≈ 2×τ_int where τ_int is the integrated autocorrelation time.
!
!> @param[in] EN Input complex bins (must have size divisible by NREBIN)
!> @param[out] ZM Jackknife mean estimate
!> @param[out] ZERR Jackknife error estimate with rebinning correction
!> @param[in] NREBIN Number of consecutive bins to combine
!
!> @warning Size of EN must be divisible by NREBIN
!
!> @see ERRCALC_J_C, AUTO_COR (for determining τ_int)
!--------------------------------------------------------------------
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           COMPLEX (Kind=Kind(0.d0))               ::  ZM, ZERR, Z
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER     :: NP, N, NP1, NREBIN, NC, NB

           NP = SIZE(EN)
           NP1 = NP/NREBIN                        ! Number of rebinned bins
           ALLOCATE (EN1(NP1))
           
           ! Rebin: average NREBIN consecutive bins
           NC = 0
           DO N = 1,NP1
              Z = CMPLX(0.D0, 0.D0, kind(0.D0))
              DO NB = 1,NREBIN
                 NC = NC + 1
                 Z = Z + EN(NC)
              ENDDO
              Z = Z/DBLE(NREBIN)                  ! Average over rebin window
              EN1(N) = Z
           ENDDO
           
           ! Apply jackknife to rebinned data
           CALL ERRCALC_J_C(EN1,ZM,ZERR)

           DEALLOCATE(EN1)

         END SUBROUTINE ERRCALC_J_C_REBIN

!******************
         SUBROUTINE ERRCALC_J_REBIN(EN,XM,XERR,NREBIN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error estimation with rebinning for real data.
!
!> @details
!> Real-valued version of rebinned jackknife. Workflow:
!> 1. Rebin data by averaging NREBIN consecutive bins
!> 2. Apply jackknife to rebinned data
!>
!> Essential for QMC data with autocorrelations. The rebinning reduces
!> correlations between bins, allowing jackknife to give accurate errors.
!
!> @param[in] EN Input bins (size must be divisible by NREBIN)
!> @param[out] XM Jackknife mean estimate
!> @param[out] XERR Jackknife error estimate with rebinning correction
!> @param[in] NREBIN Number of bins to combine (choose ~ 2×τ_int)
!
!> @see ERRCALC_J, ERRCALC_J_C_REBIN, AUTO_COR
!--------------------------------------------------------------------
           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, X
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER :: NREBIN, NC, N, NB, NP1, NP

           NP = SIZE(EN)
           NP1 = NP/NREBIN                        ! Number of rebinned bins
           ALLOCATE (EN1(NP1))

           ! Rebin: average consecutive bins
           NC = 0
           DO N = 1,NP1
              X = 0.D0
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(NC)
              ENDDO
              X = X/DBLE(NREBIN)                  ! Average over rebin window
              EN1(N) = X
           ENDDO
           
           ! Apply jackknife to rebinned data
           CALL ERRCALC_J(EN1,XM,XERR)

           DEALLOCATE(EN1)
           RETURN
         END SUBROUTINE ERRCALC_J_REBIN

!**********
         SUBROUTINE ERRCALC_JS(EN,SI,XM,XERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error for sign-problematic QMC observables.
!
!> @details
!> Computes jackknife error on O = ⟨s·O⟩/⟨s⟩ where s is the configuration sign.
!> In sign-problematic QMC simulations, the physical observable O cannot be
!> measured directly. Instead, one measures:
!> - EN = s·O (signed observable)
!> - SI = s (configuration sign)
!> 
!> The physical observable is then: O = ⟨s·O⟩ / ⟨s⟩
!>
!> Algorithm:
!> 1. Compute jackknife estimates: O_i = (Σ(s·O)_j - (s·O)_i) / (Σs_j - s_i)
!> 2. Calculate variance of O_i
!> 3. Scale by N for jackknife formula
!>
!> This properly accounts for correlations between numerator and denominator.
!
!> @param[in] EN Signed observable bins (s·O for each configuration)
!> @param[in] SI Configuration signs (s for each configuration)
!> @param[out] XM Mean of physical observable O
!> @param[out] XERR Jackknife error of O
!
!> @note
!> While designed for the sign problem, this method works for any ratio of
!> correlated observables where proper error propagation is needed.
!
!> @see ERRCALC_JS_F (with function evaluation), ERRCALC_JS_REBIN, ERRCALC_JS_C
!--------------------------------------------------------------------
           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN      ! Signed observable: s·O
           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  SI      ! Configuration signs: s
           REAL (Kind=Kind(0.d0))               ::  XM      ! Mean of O = ⟨s·O⟩/⟨s⟩
           REAL (Kind=Kind(0.d0))               ::  XERR    ! Jackknife error
           REAL (Kind=Kind(0.d0))               ::  X,XS    ! Jackknife sums (observable and sign)
           REAL (Kind=Kind(0.d0))               ::  Xhelp   ! Total sum of s·O
           REAL (Kind=Kind(0.d0))               ::  XShelp  ! Total sum of s
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE :: EN1  ! Jackknife estimates
           INTEGER                     ::  N, N1, NP, NP1

           NP = SIZE(EN)
           NP1= SIZE(SI)
           IF (NP1.NE.NP) THEN
              WRITE(error_unit,*) 'Error in Errcalc_JS: EN and SI sizes do not match'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF
           ALLOCATE (EN1(NP))

           ! Compute total sums: Σ(s·O) and Σs
           Xhelp  = 0.D0
           XShelp = 0.D0
           DO N1 = 1,NP
              Xhelp  = Xhelp  + EN(N1)   ! Sum of s·O
              XShelp = XShelp + SI(N1)   ! Sum of s
           ENDDO

           ! Build N jackknife estimates: O_i = [Σ(s·O) - (s·O)_i] / [Σs - s_i]
           DO N = 1,NP
              X  = Xhelp  - EN(N)                 ! Sum excluding bin N: Σ(s·O) - (s·O)_N
              XS = XShelp - SI(N)                 ! Sign sum excluding bin N: Σs - s_N
              EN1(N) = X / XS                     ! Jackknife ratio estimate for bin N
           ENDDO
           
           ! Calculate variance of jackknife estimates and scale by N
           CALL ERRCALC(EN1,XM,XERR)
           XERR = XERR*DBLE(NP)                   ! Jackknife error formula: err_J = N × SE
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_JS

!--------------------------------------------------------------------
         !> @author
         !> ALF-project
         !>
         !> @brief
         !> Jackknife error for functions of sign-corrected QMC observables.
         !>
         !> @details
         !> Computes error on f(O₁, O₂, ...) where each Oᵢ = ⟨s·Oᵢ⟩/⟨s⟩.
         !> Used when you need to apply a function to multiple sign-corrected observables
         !> in QMC with sign problem.
         !>
         !> Algorithm:
         !> 1. Compute jackknife sign-corrected vectors: x_i = (Σ(s·O)_j - (s·O)_i) / (Σs_j - s_i)
         !> 2. Evaluate function: f_i = f(x_i) on each jackknife sample
         !> 3. Calculate variance of f_i and scale by N
         !>
         !> Typical use: structure factors, susceptibilities, or any nonlinear function
         !> of sign-corrected QMC observables (e.g., S(q) computed from multiple correlators).
         !>
         !> @param[in] EN Multi-observable signed data (Nobs × Nbins): s·O₁, s·O₂, ...
         !> @param[in] SI Configuration sign bins (Nbins)
         !> @param[out] XM Mean of f(O₁, O₂, ...)
         !> @param[out] XERR Jackknife error of f(O₁, O₂, ...)
         !> @param[in] f Function pointer: real array -> real (must match func_r interface)
         !>
         !> @note EN has shape (Nobs, Nbins). Each jackknife sample contains Nobs sign-corrected values.
         !>
         !> @see func_r (abstract interface), ERRCALC_JS, ERRCALC_JS_C_F (complex version)
         SUBROUTINE ERRCALC_JS_F(EN, SI, XM, XERR, f)
            IMPLICIT NONE

            ! Declare input variables
            REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: EN
            REAL (Kind=Kind(0.d0)), DIMENSION(:) :: SI
            REAL (Kind=Kind(0.d0)) :: XM, XERR, XS, XShelp
            procedure (func_r), pointer :: f

            ! Declare local variables
            REAL (Kind=Kind(0.d0)), ALLOCATABLE :: EN1(:), Xhelp(:), X(:)
            INTEGER :: N, N1, NP, NP1, Nobs

            ! Get the number of bins and observations
            NP = SIZE(EN,2)
            Nobs = SIZE(EN,1)

            ! Check that the number of SI values matches the number of bins
            NP1 = SIZE(SI)
            IF (NP1.NE.NP) THEN
               WRITE(error_unit,*) 'Error in Errcalc_JS_F: EN and SI sizes do not match'
               Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            ENDIF

            ! Allocate memory for jackknife samples
            ALLOCATE (EN1(NP), Xhelp(nobs), X(nobs))

            ! Compute totals for jackknife (numerator and denominator)
            Xhelp = 0.D0
            XShelp = 0.D0
            DO N1 = 1,NP
               Xhelp = Xhelp + EN(:,N1)
               XShelp = XShelp + SI(N1)
            ENDDO

            ! Build N jackknife function evaluations
            DO N = 1,NP
               XS = XShelp - SI(N)                    ! Denominator excluding bin N
               X = (Xhelp - EN(:,N)) / XS             ! Ratio vector excluding bin N
               EN1(N) = f(X)                          ! Evaluate user function on jackknife sample
            ENDDO

            ! Compute variance and scale for jackknife
            CALL ERRCALC(EN1, XM, XERR)
            XERR = XERR * DBLE(NP)                    ! Jackknife error scaling

            ! Deallocate memory for temporary arrays
            DEALLOCATE (EN1, X, Xhelp)

            RETURN
         END SUBROUTINE ERRCALC_JS_F

!**********
         SUBROUTINE ERRCALC_JS_C(EN,SI,XM,XERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error for complex sign-problematic QMC observables.
!
!> @details
!> Complex version of ERRCALC_JS for observables O = ⟨s·O⟩/⟨s⟩.
!> Used when the physical observable is complex (e.g., Green's functions,
!> dynamical structure factors) but must be measured with sign reweighting
!> in QMC simulations with a sign problem.
!>
!> Computes: O = ⟨s·O_complex⟩ / ⟨s⟩ with proper error propagation.
!
!> @param[in] EN Complex signed observable bins (s·O)
!> @param[in] SI Complex configuration weights (typically real s, but complex for generality)
!> @param[out] XM Mean of physical complex observable
!> @param[out] XERR Jackknife error (complex)
!
!> @see ERRCALC_JS (real version), ERRCALC_JS_C_REBIN
!--------------------------------------------------------------------
           
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN, SI
           COMPLEX (Kind=Kind(0.d0))               ::  XM, XERR, X,XS, Xhelp, XShelp
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE :: EN1
           INTEGER                     ::  N, N1, NP, NP1

           NP = SIZE(EN)
           NP1= SIZE(SI)
           IF (NP1.NE.NP) THEN
              WRITE(error_unit,*) 'Error in Errcalc_JS_C: EN and SI sizes do not match'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF
           ALLOCATE (EN1(NP))

           ! Compute complex sums for jackknife
           Xhelp  = CMPLX(0.D0, 0.D0, kind(0.D0))
           XShelp = CMPLX(0.D0, 0.D0, kind(0.D0))
           DO N1 = 1,NP
              Xhelp  = Xhelp  + EN(N1)
              XShelp = XShelp + SI(N1)
           ENDDO

           ! Build N jackknife complex ratio estimates
           DO N = 1,NP
              X  = Xhelp  - EN(N)                    ! Numerator sum excluding bin N
              XS = XShelp - SI(N)                    ! Denominator sum excluding bin N
              EN1(N) = X / XS                        ! Complex ratio
           ENDDO
           
           ! Calculate variance and scale
           CALL ERRCALC(EN1,XM,XERR)
           XERR = XERR*DBLE(NP)                      ! Jackknife error scaling
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_JS_C

!**********
         SUBROUTINE ERRCALC_JS_C_F(EN,SI,XM,XERR,f)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error for complex functions of sign-corrected QMC observables.
!
!> @details
!> Complex version of ERRCALC_JS_F for f(O₁, O₂, ...) where Oᵢ = ⟨s·Oᵢ⟩/⟨s⟩
!> and observables/function are complex.
!>
!> Use for complex transformations of sign-corrected QMC data such as:
!> - Fourier transforms to momentum space
!> - Spectral functions via analytic continuation  
!> - Complex structure factors
!> - Green's functions in frequency space
!
!> @param[in] EN Multi-observable complex signed data (Nobs × Nbins): s·O₁, s·O₂, ...
!> @param[in] SI Complex weights (typically real signs, but complex for generality)
!> @param[out] XM Mean of f(O₁, O₂, ...)
!> @param[out] XERR Jackknife error of f(O₁, O₂, ...)
!> @param[in] f Function pointer: complex array -> complex (must match func_c)
!
!> @see func_c (abstract interface), ERRCALC_JS_F (real version)
!--------------------------------------------------------------------

           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:) ::  EN
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  SI
           COMPLEX (Kind=Kind(0.d0))               ::  XM, XERR, XS, XShelp
           procedure (func_c), pointer:: f
           COMPLEX (Kind=Kind(0.d0)), ALLOCATABLE  :: EN1(:), Xhelp(:),X(:)
           INTEGER                     ::  N, N1, NP, NP1, Nobs

           NP = SIZE(EN,2)
           Nobs = SIZE(EN,1)
           NP1= SIZE(SI)
           IF (NP1.NE.NP) THEN
              WRITE(error_unit,*) 'Error in Errcalc_JS_C_F: EN and SI sizes do not match'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF
           ALLOCATE (EN1(NP),Xhelp(nobs), X(nobs))
           
           ! Compute complex totals for jackknife
           Xhelp  = CMPLX(0.D0, 0.D0, kind(0.D0))
           XShelp = CMPLX(0.D0, 0.D0, kind(0.D0))
           DO N1 = 1,NP
              Xhelp  = Xhelp  + EN(:,N1)
              XShelp = XShelp + SI(N1)
           ENDDO

           ! Build N jackknife function evaluations
           DO N = 1,NP
              XS = XShelp - SI(N)                       ! Denominator excluding bin N
              X  = (Xhelp  - EN(:,N))/Xs                ! Complex ratio vector
              EN1(N) = f(X)                             ! Evaluate complex function
           ENDDO
           
           ! Calculate variance and scale
           CALL ERRCALC(EN1,XM,XERR)
           XERR = XERR*DBLE(NP)                         ! Jackknife error scaling
           DEALLOCATE  ( EN1, X, Xhelp )

           RETURN
         END SUBROUTINE ERRCALC_JS_C_F



!********
         SUBROUTINE ERRCALC_JS_REBIN(EN,SI,XM,XERR,NREBIN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error for sign-problematic QMC with rebinning.
!
!> @details
!> Combines rebinning with jackknife sign method for autocorrelated data.
!> Computes O = ⟨s·O⟩/⟨s⟩ with proper treatment of autocorrelations.
!>
!> Workflow:
!> 1. Rebin both s·O (EN) and s (SI) independently to reduce autocorrelation
!> 2. Apply jackknife sign method to rebinned data
!>
!> Essential for QMC simulations where both sign problem and autocorrelations
!> are present (typical for fermionic systems).
!
!> @param[in] EN Signed observable bins (s·O, size divisible by NREBIN)
!> @param[in] SI Configuration sign bins (s, same size as EN)
!> @param[out] XM Mean of physical observable O
!> @param[out] XERR Jackknife error with rebinning correction
!> @param[in] NREBIN Number of bins to combine (choose ~ 2×τ_int)
!
!> @see ERRCALC_JS, ERRCALC_JS_C_REBIN, AUTO_COR
!--------------------------------------------------------------------
           
           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN, SI
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, X, Y
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1, SI1
           INTEGER :: NREBIN, NC, N, NB, NP, NP1

           NP = SIZE(EN)
           NP1 = NP/NREBIN                           ! Number of rebinned bins
           ALLOCATE (EN1(NP1))
           ALLOCATE (SI1(NP1))

           ! Rebin both numerator and denominator
           NC = 0
           DO N = 1,NP1
              X = 0.D0; Y = 0.D0
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(NC)                      ! Accumulate numerator
                 Y = Y + SI(NC)                      ! Accumulate denominator
              ENDDO
              X = X/DBLE(NREBIN)                     ! Average numerator
              Y = Y/DBLE(NREBIN)                     ! Average denominator
              EN1(N) = X
              SI1(N) = Y
           ENDDO
           
           ! Apply jackknife ratio to rebinned data
           CALL ERRCALC_JS(EN1,SI1,XM,XERR)

           DEALLOCATE (EN1,SI1)

           RETURN
         END SUBROUTINE ERRCALC_JS_REBIN



!********
         SUBROUTINE ERRCALC_JS_REBIN_F(EN,SI,XM,XERR,NREBIN,f)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error for functions of ratios with rebinning.
!
!> @details
!> Combines rebinning and jackknife for function evaluation on ratios.
!> Workflow:
!> 1. Rebin both EN and SI data
!> 2. Apply jackknife with function evaluation to rebinned data
!>
!> Use for derived quantities computed from correlated QMC observables.
!
!> @param[in] EN Multi-observable numerator (Nobs × Nbins)
!> @param[in] SI Denominator bins (Nbins)
!> @param[out] XM Mean of f(EN/SI)
!> @param[out] XERR Jackknife error with rebinning
!> @param[in] NREBIN Number of bins to combine
!> @param[in] f Function pointer (func_r interface)
!
!> @see ERRCALC_JS_F, ERRCALC_JS_REBIN
!--------------------------------------------------------------------
           
           IMPLICIT NONE

           REAL (Kind=Kind(0.d0))               ::  EN(:,:), SI(:)
           procedure (func_r), pointer:: f
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, Y
           REAL (Kind=Kind(0.d0)), ALLOCATABLE  ::  EN1(:,:), SI1(:), X(:)
           INTEGER :: NREBIN, NC, N, NB, NP, NP1, nobs
 
           NP = SIZE(EN,2)
           Nobs = SIZE(EN,1)
           NP1 = NP/NREBIN
           ALLOCATE (EN1(nobs,NP1), X(nobs))
           ALLOCATE (SI1(NP1))
           
           ! Rebin both observables and weights
           NC = 0
           DO N = 1,NP1
              X = 0.D0; Y = 0.D0
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(:,NC)                    ! Accumulate observable vector
                 Y = Y + SI(NC)                      ! Accumulate weights
              ENDDO
              X = X/DBLE(NREBIN)
              Y = Y/DBLE(NREBIN)
              EN1(:,N) = X
              SI1(N) = Y
           ENDDO
           
           ! Apply jackknife with function
           CALL ERRCALC_JS_F(EN1,SI1,XM,XERR,f)

           DEALLOCATE (EN1,SI1)

           RETURN
         END SUBROUTINE ERRCALC_JS_REBIN_F

!******************
         SUBROUTINE ERRCALC_JS_C_REBIN(EN,SI,XM,XERR,NREBIN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error for complex sign-problematic QMC with rebinning.
!
!> @details
!> Complex version of ERRCALC_JS_REBIN. Computes O = ⟨s·O⟩/⟨s⟩ for complex
!> observables with autocorrelations. Essential for complex QMC observables
!> like Green's functions, dynamical correlators, etc. in sign-problematic
!> simulations.
!
!> @param[in] EN Complex signed observable bins (s·O)
!> @param[in] SI Complex weight bins (typically real signs)
!> @param[out] XM Mean of physical complex observable O
!> @param[out] XERR Jackknife error with rebinning correction
!> @param[in] NREBIN Number of bins to combine
!
!> @see ERRCALC_JS_C, ERRCALC_JS_REBIN (real version)
!--------------------------------------------------------------------
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN, SI
           COMPLEX (Kind=Kind(0.d0))               ::  XM, XERR, X, Y
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1, SI1
           INTEGER :: NREBIN, NC, N, NB, NP, NP1

           NP = SIZE(EN)
           NP1 = NP/NREBIN
           ALLOCATE (EN1(NP1))
           ALLOCATE (SI1(NP1))

           ! Rebin complex data
           NC = 0
           DO N = 1,NP1
              X = cmplx(0.D0,0.d0,kind(0.d0))      ! Initialize numerator
              Y = cmplx(0.D0,0.D0,kind(0.d0))      ! Initialize denominator
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(NC)
                 Y = Y + SI(NC)
              ENDDO
              X = X/DBLE(NREBIN)                   ! Average
              Y = Y/DBLE(NREBIN)
              EN1(N) = X
              SI1(N) = Y
           ENDDO
           
           ! Apply complex jackknife ratio
           CALL ERRCALC_JS_C(EN1,SI1,XM,XERR)

           DEALLOCATE (EN1,SI1)

           RETURN
         END SUBROUTINE ERRCALC_JS_C_REBIN

!******************
         SUBROUTINE ERRCALC_JS_C_REBIN_F(EN,SI,XM,XERR,NREBIN,f)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife error for complex functions of ratios with rebinning.
!
!> @details
!> Most general jackknife variant: combines rebinning, ratio estimation,
!> and function evaluation for complex observables. Use for complex
!> transformations of correlated QMC data.
!
!> @param[in] EN Multi-observable complex data (Nobs × Nbins)
!> @param[in] SI Complex weight bins
!> @param[out] XM Mean of f(EN/SI)
!> @param[out] XERR Jackknife error with rebinning
!> @param[in] NREBIN Number of bins to combine
!> @param[in] f Complex function pointer (func_c interface)
!
!> @see ERRCALC_JS_C_F, ERRCALC_JS_REBIN_F (real version)
!--------------------------------------------------------------------
           
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0))               ::  EN(:,:), SI(:)
           COMPLEX (Kind=Kind(0.d0))               ::  XM, XERR, Y
           procedure (func_c), pointer:: f
           COMPLEX (Kind=Kind(0.d0)), ALLOCATABLE  ::  EN1(:,:), SI1(:), X(:)
           INTEGER :: NREBIN, NC, N, NB, NP, NP1, nobs
 
           Nobs = SIZE(EN,1)
           NP = SIZE(EN,2)
           NP1 = NP/NREBIN
           ALLOCATE (EN1(nobs,NP1))
           ALLOCATE (SI1(NP1), X(nobs))
           
           ! Rebin complex multi-observable data
           NC = 0
           DO N = 1,NP1
              X = cmplx(0.D0,0.d0,kind(0.d0))      ! Initialize observable vector
              Y = cmplx(0.D0,0.D0,kind(0.d0))      ! Initialize weight
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(:,NC)
                 Y = Y + SI(NC)
              ENDDO
              X = X/DBLE(NREBIN)
              Y = Y/DBLE(NREBIN)
              EN1(:,N) = X
              SI1(N) = Y
           ENDDO
           
           ! Apply complex jackknife with function
           CALL ERRCALC_JS_C_F(EN1,SI1,XM,XERR,f)

           DEALLOCATE (EN1,SI1,X)

           RETURN
         END SUBROUTINE ERRCALC_JS_C_REBIN_F

!******************

         SUBROUTINE INTER_QMC(GR, SIGN1, DTAU, RES, ERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Trapezoid rule integration with jackknife error for sign-corrected QMC data.
!
!> @details
!> Integrates imaginary-time dependent observables with sign correction.
!> For sign-problematic QMC, computes: ∫₀^β O(τ) dτ where O(τ) = <s·G(τ)>/<s>
!>
!> Algorithm:
!> 1. Sign-correct each time slice: O(τ) = Σ(s·G)/Σs (jackknife)
!> 2. Integrate using trapezoid rule: ∫ = Δτ Σ (O_i + O_{i+1})/2
!> 3. Calculate jackknife error on integral
!>
!> Typical use: ∫₀^β G(τ) dτ for thermodynamic quantities in sign-problematic QMC.
!
!> @param[in] GR Signed time-dependent data (Ntimes × Nbins): s·O(τ) for each bin
!> @param[in] SIGN1 Configuration signs (Nbins): s for each bin
!> @param[in] DTAU Time step (Δτ)
!> @param[out] RES Integral value
!> @param[out] ERR Jackknife error on integral
!
!> @note Sign is assumed constant across all time slices (typical in QMC).
!
!> @see ERRCALC_JS, INTER_F
!--------------------------------------------------------------------

           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the integral and error
           ! The sign is the same for all Times.
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   ::   SIGN1

           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, Y, Err, Res, DTAU, Xhelp, Yhelp
           INTEGER :: NT, NB, NB1, NTDM, NDATA

           NTDM  = SIZE(GR,1)                      ! Number of time slices
           NDATA = SIZE(GR,2)                      ! Number of bins


           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA) )
           ! HLP1(τ,j) = sign-corrected O(τ) excluding bin j
           ! HLP(j) = integrated value for jackknife sample j
           
           ! Step 1: Build sign-corrected jackknife time series for each bin
           DO NT = 1,NTDM
              Xhelp = 0.d0                         ! Total: Σ_j s·G(τ,j)
              Yhelp = 0.d0                         ! Total: Σ_j s_j
              DO NB1 = 1,NDATA
                  Xhelp = Xhelp + GR(NT,NB1)
                  Yhelp = Yhelp + SIGN1(NB1)
              ENDDO
              ! Create jackknife estimate for each time slice
              DO NB= 1, NDATA
                 X = Xhelp - GR(NT,NB)             ! Σ s·G excluding bin NB
                 Y = Yhelp - SIGN1(NB)             ! Σ s excluding bin NB
                 HLP1(NT,NB) = X/Y                 ! O(τ) = ⟨s·G(τ)⟩/⟨s⟩ for sample NB
              ENDDO
           ENDDO

           ! Step 2: Integrate each jackknife time series using trapezoid rule
           DO NB = 1,NDATA
              X = 0.D0
              DO NT = 1,NTDM-1
                 ! Trapezoid: ∫ ≈ Σ (f_i + f_{i+1})/2 × Δτ
                 X = X + (HLP1(NT,NB) + HLP1(NT+1,NB))*0.5D0
              ENDDO
              HLP(NB) = X * DTAU                   ! Scale by time step Δτ
           ENDDO

           ! Step 3: Calculate jackknife error on the N integrated values
           CALL ERRCALC(HLP, RES, ERR)
           ERR = ERR*DBLE(NDATA)                   ! Jackknife error scaling: N × SE

           DEALLOCATE( HLP, HLP1 )

           RETURN
         END SUBROUTINE INTER_QMC

!******************
         REAL (Kind=Kind(0.d0)) FUNCTION INTER_F(A,B,N,F)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Numerical integration using trapezoid rule for function evaluation.
!
!> @details
!> Integrates a user-supplied function F from A to B using N equally-spaced
!> points with the trapezoid rule:
!> ∫ₐᵇ F(x) dx ≈ h[F(x₀)/2 + F(x₁) + ... + F(x_{N-1}) + F(x_N)/2]
!> where h = (B-A)/N
!
!> @param[in] A Lower integration bound
!> @param[in] B Upper integration bound
!> @param[in] N Number of integration points (more points = higher accuracy)
!> @param[in] F External function to integrate
!> @return Integral value
!
!> @note Accuracy is O(h²) for smooth functions.
!
!> @see INTER_QMC (for QMC time-dependent data)
!--------------------------------------------------------------------

           IMPLICIT NONE

           INTEGER  :: N, I
           REAL (Kind=Kind(0.d0)) ::  A, B, X, X1
           REAL (Kind=Kind(0.d0)), EXTERNAL :: F

           REAL (Kind=Kind(0.d0)) ::  DEL

           DEL = (B-A)/DBLE(N)                     ! Integration step size
           INTER_F = 0.D0
           
           ! Trapezoid rule summation
           DO I = 0, N-1
              X  = A + DBLE(I  )*DEL
              X1 = A + DBLE(I+1)*DEL
              INTER_F = INTER_F + ( F(X) + F(X1) )*0.5D0
           ENDDO
           INTER_F = INTER_F*DEL
         END FUNCTION INTER_F

!****************** Least square fits:
         SUBROUTINE FIT(XDATA,FDATA,ERROR,ARES,CHSQ,F)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Least squares curve fitting using singular value decomposition (SVD).
!
!> @details
!> Fits data to a linear combination of basis functions using SVD.
!> Minimizes χ² = Σ[(FDATA_i - Σ_j A_j·F(j,x_i))² / ERROR_i²]
!>
!> Algorithm:
!> 1. Construct design matrix: A_{ij} = F(j, x_i) / σ_i
!> 2. SVD decomposition: A = U·D·Vᵀ
!> 3. Solve: a = V·D⁻¹·Uᵀ·b (where b = FDATA/ERROR)
!> 4. Compute reduced χ²
!>
!> Typical use: polynomial fits, exponential decay fits in QMC.
!
!> @param[in] XDATA X-coordinates of data points
!> @param[in] FDATA Y-values (data to fit)
!> @param[in] ERROR Uncertainties on FDATA
!> @param[out] ARES Fitted coefficients (A_1, A_2, ..., A_NBASIS)
!> @param[out] CHSQ Reduced chi-squared (goodness of fit)
!> @param[in] F Basis function: F(basis_index, x)
!
!> @note F must return the value of basis function #M at position X.
!> @note ARES size determines number of basis functions.
!
!> @see UDV (SVD decomposition), INV (matrix inversion)
!--------------------------------------------------------------------

           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  XDATA, FDATA, ERROR, ARES
           REAL (Kind=Kind(0.d0))               ::  CHSQ,  X
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:),  ALLOCATABLE :: A, U,V,VINV,V1
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ),  ALLOCATABLE :: B,D
           REAL (Kind=Kind(0.d0)), EXTERNAL  :: F
           INTEGER                     ::  NDATA, NBASIS, I, M, M1, NCON, N

           NDATA = SIZE(XDATA)
           NBASIS= SIZE(ARES)                      ! Number of basis functions

           ALLOCATE (A(NDATA,NBASIS))
           ALLOCATE (U(NDATA,NBASIS))
           ALLOCATE (D(NBASIS))
           ALLOCATE (V   (NBASIS,NBASIS))
           ALLOCATE (V1  (NBASIS,NBASIS))
           ALLOCATE (VINV(NBASIS,NBASIS))
           ALLOCATE (B(NDATA))

           ! Initialize matrices
           A = 0.D0
           U = 0.D0
           D = 0.D0
           V = 0.D0
           VINV = 0.D0
           V1 = 0.D0
           B = 0.D0
           NCON = 0
           
           ! Build design matrix (weighted by errors)
           DO M = 1,NBASIS
              DO I = 1,NDATA
                 A(I,M) = F(M,XDATA(I))/ERROR(I)
              ENDDO
           ENDDO
           
           ! Build right-hand side (weighted data)
           DO I = 1,NDATA
              B(I) = FDATA(I)/ERROR(I)
           ENDDO
           
           ! SVD decomposition: A = U·D·V^T
           CALL UDV(A,U,D,V,NCON)
           
           ! Transpose V for inversion
           DO M = 1,NBASIS
              DO I = 1,NBASIS
                 V1(I,M) = V(M,I)
              ENDDO
           ENDDO
           CALL INV(V1,VINV,X)

           ! Solve for coefficients: ARES = V·D^(-1)·U^T·B
           DO M1 = 1,NBASIS
              X = 0.D0
              DO M = 1,NBASIS
                 DO I = 1,NDATA
                    X = X + B(I)*U(I,M)*VINV(M,M1)/D(M)
                 ENDDO
              ENDDO
              ARES(M1) = X
           ENDDO

           ! Compute reduced chi-squared
           CHSQ = 0.D0
           DO N = 1,NDATA
              X = 0.D0
              DO M = 1,NBASIS
                 X = X + ARES(M)*F(M,XDATA(N))       ! Evaluate fitted function
              ENDDO
              CHSQ = CHSQ + (FDATA(N) - X)**2/ERROR(N)**2
           ENDDO
           CHSQ = CHSQ/DBLE(NDATA)                   ! Normalize by N

           DEALLOCATE (A)
           DEALLOCATE (U)
           DEALLOCATE (D)
           DEALLOCATE (V)
           DEALLOCATE (V1)
           DEALLOCATE (VINV)
           DEALLOCATE (B)

         END SUBROUTINE FIT

         SUBROUTINE COVJ(GR, XCOV, XMEAN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife covariance matrix calculation.
!
!> @details
!> Computes the covariance matrix Cov(O_i, O_j) using jackknife resampling.
!> For time-dependent observables: O_i and O_j at different times.
!> 
!> Covariance: Cov(i,j) = N × ⟨(̅O_i - ⟨O_i⟩)(̅O_j - ⟨O_j⟩)⟩_jackknife
!> where ̅O denotes jackknife estimates.
!
!> @param[in] GR Observable data (Ntimes × Nbins)
!> @param[out] XCOV Covariance matrix (Ntimes × Ntimes)
!> @param[out] XMEAN Mean values at each time (Ntimes)
!
!> @note Diagonal elements are variances; off-diagonal are covariances.
!> @note Essential for correlated error analysis and multi-parameter fits.
!
!> @see COVJS (with sign weighting), COVJS_C (complex version)
!--------------------------------------------------------------------

           IMPLICIT NONE
           !Given GR(Times, Bins)  calculates the mean and the covariance.
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN

           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR, Xhelp
           INTEGER :: NT, NT1, NB, NB1, NTDM, NDATA

           NTDM  = SIZE(GR,1)                      ! Number of observables/times
           NDATA = SIZE(GR,2)                      ! Number of bins

           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(error_unit,*) 'Error in COVJ: XCOV dimension mismatch'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF

           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA) )
           ! HLP1(i,j) = jackknife estimate for observable i, excluding bin j
           ! HLP(j) = temporary storage for jackknife estimates at fixed time
           
           ! Compute jackknife estimates: ̅O_i,j = (Σ_k O_i,k - O_i,j)/(N-1)
           DO NT = 1,NTDM
              Xhelp = 0.0
              ! Sum all bins for observable NT
              DO NB1 = 1,NDATA
                 Xhelp = Xhelp + GR(NT,NB1)
              ENDDO
              ! Build jackknife estimate excluding each bin
              DO NB= 1, NDATA
                 X = Xhelp - GR(NT,NB)                   ! Sum excluding bin NB
                 HLP1(NT,NB) = X/DBLE(NDATA-1)           ! Jackknife average
                 HLP (NB   ) = X/DBLE(NDATA-1)
              ENDDO
              ! Compute mean over jackknife samples
              CALL ERRCALC(HLP,XM ,XERR)
              XMEAN(NT) = XM
           ENDDO


           ! Compute covariance matrix: Cov(i,j) = N × ⟨(̅O_i - ⟨̅O_i⟩)(̅O_j - ⟨̅O_j⟩)⟩
           DO NT = 1,NTDM
              DO NT1= 1,NTDM
                 X = 0.0
                 ! Sum products of deviations across jackknife samples
                 DO NB = 1,NDATA
                    X = X + (HLP1(NT,NB)-XMEAN(NT))*(HLP1(NT1,NB)-XMEAN(NT1))
                 ENDDO
                 X = X/DBLE(NDATA)                       ! Average over samples
                 XCOV(NT,NT1)  = X*DBLE(NDATA)           ! Jackknife scaling factor N
              ENDDO
           ENDDO


           DEALLOCATE( HLP, HLP1 )

           RETURN
         END SUBROUTINE COVJ


         SUBROUTINE COVJS(GR, SIGN1, XCOV, XMEAN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife covariance matrix for sign-problematic QMC.
!
!> @details
!> Computes covariance matrix for observables requiring sign correction:
!> O(τ) = ⟨s·O(τ)⟩ / ⟨s⟩
!>
!> In sign-problematic fermionic QMC, time-dependent observables must be
!> reweighted by the configuration sign. This routine computes the full
!> covariance matrix Cov[O(τ), O(τ')] accounting for:
!> - Sign reweighting correlations
!> - Temporal correlations
!> - Proper error propagation through the ratio
!>
!> Essential for correlated fits of imaginary-time data in systems with sign problem.
!
!> @param[in] GR Signed observable data (Ntimes × Nbins): s·O(τ) for each bin
!> @param[in] SIGN1 Configuration signs (Nbins): s for each bin
!> @param[out] XCOV Covariance matrix of physical observable (Ntimes × Ntimes)
!> @param[out] XMEAN Sign-corrected mean values: O(τ) = ⟨s·O(τ)⟩/⟨s⟩
!
!> @note Sign is assumed constant across all time slices (typical in QMC).
!
!> @see COVJ (without sign problem), COVJS_C (complex version)
!--------------------------------------------------------------------

           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times.
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN, SIGN1

           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR, Y, Xhelp, Yhelp
           INTEGER :: NT, NT1, NB, NTDM, NDATA

           NTDM  = SIZE(GR,1)                      ! Number of observables/times
           NDATA = SIZE(GR,2)                      ! Number of bins

           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(error_unit,*) 'Error in COVJS: XCOV dimension mismatch'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF

           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA) )
           
           ! Compute jackknife ratio estimates
           DO NT = 1,NTDM
              Xhelp=0.d0                           ! Sum of observables
              Yhelp=0.d0                           ! Sum of signs
              DO NB= 1, NDATA
                 Xhelp=Xhelp+GR(NT,NB)
                 Yhelp=Yhelp+SIGN1(NB)
              ENDDO
              DO NB= 1, NDATA
                 X = Xhelp - GR(NT,NB)             ! Observable sum excluding bin NB
                 Y = Yhelp - SIGN1(NB)             ! Sign sum excluding bin NB
                 HLP1(NT,NB) = X/Y                 ! Jackknife ratio estimate
                 HLP (NB   ) = X/Y
              ENDDO
              CALL ERRCALC(HLP,XM ,XERR)
              XMEAN(NT) = XM
           ENDDO


           ! Compute covariance matrix from jackknife estimates
           DO NT = 1,NTDM
              DO NT1= 1,NTDM
                 X = 0.0
                 DO NB = 1,NDATA
                    X = X +  (HLP1(NT,NB)-XMEAN(NT))*(HLP1(NT1,NB)-XMEAN(NT1))
                 ENDDO
                 X = X/DBLE(NDATA)
                 XCOV(NT,NT1)  = X*DBLE(NDATA)     ! Jackknife scaling
              ENDDO
           ENDDO


           DEALLOCATE( HLP, HLP1 )

           RETURN
         END SUBROUTINE COVJS




         SUBROUTINE COVJS_C(GR, SIGN1, XCOV, XMEAN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Jackknife covariance matrix for complex sign-problematic QMC.
!
!> @details
!> Complex version of COVJS. Computes covariance matrix for complex
!> time-dependent observables with sign correction: O(τ) = <s·O(τ)>/<s>
!>
!> Processes real and imaginary parts independently:
!> - Real part: O_re(τ) = <s·Re[O(τ)]>/<s>
!> - Imaginary part: O_im(τ) = <s·Im[O(τ)]>/<s>
!>
!> Essential for:
!> - Complex Green's functions in sign-problematic systems
!> - Dynamical correlation functions
!> - Complex structure factors
!>
!> The covariance matrix includes proper error propagation through sign correction.
!
!> @param[in] GR Complex signed observable (Ntimes × Nbins): s·O(τ) for each bin
!> @param[in] SIGN1 Configuration signs (Nbins): s for each bin
!> @param[out] XCOV Complex covariance matrix (Ntimes × Ntimes)
!> @param[out] XMEAN Sign-corrected complex mean: O(τ) = <s·O(τ)>/<s>
!
!> @note Sign is assumed constant across all time slices.
!
!> @see COVJS (real version), COVJS_C_REBIN
!--------------------------------------------------------------------

           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times.
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           Complex (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)   ::  SIGN1


           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP, XMEAN_R
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR, Y, Xhelp, Yhelp
           INTEGER :: NT, NT1, NB, NTDM, NDATA, Nth
           COMPLEX (Kind=Kind(0.d0)) :: Z

           NTDM  = SIZE(GR,1)                      ! Number of times
           NDATA = SIZE(GR,2)                      ! Number of bins


           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(error_unit,*) 'Error in COVJS_C: XCOV dimension mismatch'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF



           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA), XMEAN_R(NTDM) )
           ! HLP1 = jackknife estimates for current real/imag component
           ! XMEAN_R = real-valued mean for current component
           XMEAN = CMPLX(0.d0, 0.d0, kind(0.D0))
           XCOV  = CMPLX(0.d0, 0.d0, kind(0.D0))

           ! Process real and imaginary parts independently (NTH=1: Re, NTH=2: Im)
           DO NTH = 1,2
              Z = CMPLX(1.D0, 0.D0, kind(0.D0))    ! Z=1: extracts real part
              IF (NTH .EQ. 2 ) Z = CMPLX( 0.D0, -1.D0, kind(0.D0))  ! Z=-i: extracts imag as real
              
              ! Compute sign-corrected jackknife estimates for each time slice
              DO NT = 1,NTDM
                 Xhelp=0.d0                        ! Total: Σ Re[Z·s·O(τ)]
                 Yhelp=0.d0                        ! Total: Σ s
                 DO NB= 1, NDATA
                    Xhelp = Xhelp + DBLE ( Z*GR(NT,NB) )  ! Project to real axis via Z
                    Yhelp = Yhelp + SIGN1(NB)
                 ENDDO
                 ! Build jackknife ratio estimates
                 DO NB= 1, NDATA
                    X = Xhelp - DBLE ( Z*GR(NT,NB) )      ! Sum excluding bin NB
                    Y = Yhelp - SIGN1(NB)
                    HLP1(NT,NB) = X/Y              ! Jackknife: ⟨s·Re[O]⟩_-NB / ⟨s⟩_-NB
                    HLP (NB   ) = X/Y
                 ENDDO
                 CALL ERRCALC(HLP,XM ,XERR)
                 XMEAN(NT) = XMEAN(NT) + CONJG(Z)*XM  ! Accumulate: Re + i·Im
                 XMEAN_R(NT) = XM
              ENDDO

              ! Compute covariance contribution from this component
              DO NT = 1,NTDM
                 DO NT1= 1,NTDM
                    X = 0.d0
                    DO NB = 1,NDATA
                       X = X +  (HLP1(NT,NB)-XMEAN_R(NT))*(HLP1(NT1,NB)-XMEAN_R(NT1))
                    ENDDO
                    X = X/DBLE(NDATA)              ! Average over jackknife samples
                    XCOV(NT,NT1)  = XCOV(NT,NT1) + CONJG(Z)* X *DBLE(NDATA)  ! Scale and accumulate
                 ENDDO
              ENDDO
           ENDDO

           DEALLOCATE( HLP, HLP1, XMEAN_R )

           RETURN
         END SUBROUTINE COVJS_C


!========================

         SUBROUTINE COVJS_C_REBIN(GR, SIGN2, XCOV, XMEAN,NREBIN)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Complex covariance with sign correction and rebinning.
!
!> @details
!> Combines rebinning with complex covariance calculation for
!> sign-problematic QMC data with autocorrelations.
!>
!> Workflow:
!> 1. Rebin both s·O(τ) and s to reduce autocorrelations
!> 2. Compute covariance matrix via COVJS_C on rebinned data
!>
!> Essential for complex time-dependent QMC observables in systems with
!> both sign problem and significant autocorrelation times.
!
!> @param[in] GR Complex signed observable (Ntimes × Nbins): s·O(τ)
!> @param[in] SIGN2 Configuration signs (Nbins)
!> @param[out] XCOV Complex covariance matrix (Ntimes × Ntimes)
!> @param[out] XMEAN Sign-corrected complex mean
!> @param[in] NREBIN Number of bins to combine
!
!> @see COVJS_C, COVJS_C_REBIN_BG
!--------------------------------------------------------------------

           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times.
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           Complex (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)   ::  SIGN2
           INTEGER :: NREBIN


           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE ::  GR1
           REAL    (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  SIGN1

           INTEGER :: NTDM, NDATA, NDATA1, N, NB, NC, NT

           NTDM  = SIZE(GR,1)                      ! Number of times
           NDATA = SIZE(GR,2)                      ! Number of bins

           NDATA1 = NDATA/NREBIN                   ! Number of rebinned bins
           ALLOCATE ( GR1(NTDM,NDATA1), SIGN1(NDATA1) )

           SIGN1 = 0.d0
           GR1   = CMPLX(0.d0,0.d0,kind(0.d0))
           
           ! Rebin: group consecutive bins (average NREBIN bins -> 1 effective bin)
           NC = 0
           DO N = 1,NDATA1
              DO NB = 1,NREBIN
                 NC = NC + 1
                 SIGN1(N) = SIGN1(N) + SIGN2(NC)         ! Accumulate signs
                 DO NT = 1,NTDM
                    GR1(NT,N)  = GR1(NT,N) + GR(NT,NC)   ! Accumulate observables
                 ENDDO
              ENDDO
           ENDDO
           SIGN1 = SIGN1/DBLE (NREBIN)             ! Average signs in each bin
           GR1   = GR1  /CMPLX(DBLE(NREBIN),0.d0,KIND(0.d0))  ! Average observables

           ! Compute covariance on rebinned data (reduces autocorrelation)
           CALL COVJS_C(GR1, SIGN1, XCOV, XMEAN)

           DEALLOCATE ( GR1, SIGN1 )

         END SUBROUTINE COVJS_C_REBIN


!========================
         SUBROUTINE COVJS_C_BG(Gr, SIGN1, XCOV, XMEAN, background)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Complex covariance with sign problem and background subtraction.
!
!> @details
!> Computes sign-averaged covariance matrix with background subtraction
!> for multi-orbital observables. Handles QMC sign problem by computing
!> ⟨s·O⟩/⟨s⟩ and subtracts squared background: ⟨O_bg⟩².
!>
!> Implements jackknife error estimation for:
!>   Trace[Ĝ(τ)] - (⟨O_bg⟩/⟨s⟩)²
!>
!> The background subtraction removes time-independent contributions
!> from multi-orbital correlation functions, important for extracting
!> proper fluctuation spectra in fermionic systems.
!>
!> Real and imaginary parts treated separately for complex arithmetic.
!
!> @param[in] Gr Observable data (Ntimes × Nbins), complex
!> @param[in] SIGN1 Configuration signs (Nbins)
!> @param[out] XCOV Covariance matrix (Ntimes × Ntimes), complex
!> @param[out] XMEAN Background-subtracted mean values (Ntimes), complex
!> @param[in] background Orbital backgrounds (Norb × Nbins), complex
!
!> @see COVJS_C, COVJS_C_BG_Weights
!--------------------------------------------------------------------
           
           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times.

           Complex (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV, background
           Complex (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)   ::  SIGN1


           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP, XMEAN_R
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR, Y, Xhelp, Yhelp
           INTEGER :: NT, NT1, NB, NTDM, NDATA, Nth, Norb, no
           COMPLEX (Kind=Kind(0.d0)) :: Z, tmp
           complex (kind=kind(0.d0)), dimension(:), allocatable :: BGME, BGJ

           NTDM  = SIZE(GR,1)
           NDATA = SIZE(GR,2)
           Norb  = SIZE(background,1)

           allocate ( BGME(Norb), BGJ(Norb) )

           !Write(6,*) 'Errors.F90 ', NTDM, NDATA
           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(error_unit,*) 'Error in COVJS_C'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF



           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA), XMEAN_R(NTDM) )
           XMEAN = CMPLX(0.d0, 0.d0, kind(0.D0))
           XCOV  = CMPLX(0.d0, 0.d0, kind(0.D0))
           BGME  = CMPLX(0.d0, 0.d0, kind(0.D0))

           ! Compute total background sum across all bins: BGME = Σ_b O_bg(b)
           Do NB=1,NDATA
              BGME(:) = BGME(:) + background(:,NB)
           End Do

           ! Process real and imaginary parts independently
           DO NTH = 1,2
              Z = CMPLX(1.D0, 0.D0, kind(0.D0))    ! Z=1: real part
              IF (NTH .EQ. 2 ) Z = CMPLX( 0.D0, -1.D0, kind(0.D0))  ! Z=-i: imag→real
              DO NT = 1,NTDM
                 Xhelp=0.d0                        ! Total: Σ Re[Z·O(τ)]
                 Yhelp=0.d0                        ! Total: Σ s
                 DO NB= 1, NDATA
                    Xhelp = Xhelp + DBLE ( Z*GR(NT,NB) )
                    Yhelp = Yhelp + SIGN1(NB)
                 ENDDO
                 ! Build jackknife with background subtraction
                 DO NB= 1, NDATA
                    Y = Yhelp - SIGN1(NB)          ! ⟨s⟩ excluding bin NB
                    BGJ = BGME - background(:,NB)  ! Background sum excluding bin NB
                    ! Compute background contribution: Σ_n |⟨O_bg,n⟩|²/⟨s⟩
                    tmp = 0.d0
                    DO No = 1,Norb
                       tmp = tmp + BGJ(No)*BGJ(No)/Y
                    enddo
                    ! Observable minus background: ⟨s·O⟩/⟨s⟩ - (⟨O_bg⟩/⟨s⟩)²
                    X = Xhelp - DBLE ( Z*GR(NT,NB) ) - DBLE ( Z*tmp )
                    HLP1(NT,NB) = X/Y
                    HLP (NB   ) = X/Y
                 ENDDO
                 CALL ERRCALC(HLP,XM ,XERR)
                 XMEAN(NT) = XMEAN(NT) + CONJG(Z)*XM  ! Reconstruct complex: Re + i·Im
                 XMEAN_R(NT) = XM
              ENDDO

              ! Compute covariance contribution from this component
              DO NT = 1,NTDM
                 DO NT1= 1,NTDM
                    X = 0.d0
                    DO NB = 1,NDATA
                       X = X +  (HLP1(NT,NB)-XMEAN_R(NT))*(HLP1(NT1,NB)-XMEAN_R(NT1))
                    ENDDO
                    X = X/DBLE(NDATA)
                    XCOV(NT,NT1)  = XCOV(NT,NT1) + CONJG(Z)* X *DBLE(NDATA)
                 ENDDO
              ENDDO
           ENDDO

           DEALLOCATE( HLP, HLP1, XMEAN_R, BGME, BGJ )

           RETURN
         END SUBROUTINE COVJS_C_BG

!========================
         SUBROUTINE COVJS_C_BG_Weights(Gr_in, SIGN1, XCOV, XMEAN, background,L_Back, Weights)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Multi-orbital covariance with sign problem, background, and weights.
!
!> @details
!> Computes sign-averaged covariance matrix for multi-orbital Green's
!> functions with optional orbital weighting and background subtraction.
!> Handles the QMC sign problem via jackknife method on ⟨s·O⟩/⟨s⟩.
!>
!> **Without weights (Weights not present):**
!> Computes trace over orbital indices:
!>   Xmean(τ) = Σₙ ⟨Ô†ₙ(τ)Ôₙ(0)⟩/⟨s⟩ - (⟨O_bg,n⟩/⟨s⟩)²
!>
!> **With weights (Weights present):**
!> Computes weighted combination:
!>   M(τ,b) = Σₙ Weights(n) × Oₙ(τ,b)
!>   Xmean(τ) = ⟨M†(τ)M(0)⟩/⟨s⟩ - (⟨M_bg⟩/⟨s⟩)²
!>
!> Input data structure: Gr_in(τ, n, m, b) represents orbital Green's
!> function ⟨Ô†ₙ(τ,b) Ôₘ(0,b)⟩ for bin b, orbitals n,m, time τ.
!>
!> Real and imaginary parts handled separately for complex arithmetic.
!
!> @param[in] Gr_in Multi-orbital data (Ntimes × Norb × Norb × Nbins)
!> @param[in] SIGN1 Configuration signs (Nbins)
!> @param[out] XCOV Covariance matrix (Ntimes × Ntimes), complex
!> @param[out] XMEAN Background-subtracted mean values (Ntimes), complex
!> @param[in] background Orbital backgrounds (Norb × Nbins), complex
!> @param[in] L_Back Enable/disable background subtraction
!> @param[in] Weights [optional] Orbital weights (Norb), complex
!
!> @see COVJS_C_BG_Rebin_Weights, COVJS_C_BG
!--------------------------------------------------------------------
!!!!===================================================
!!!        On Input   Gr(tau,n,m,b)       =  O^dag_n(tau,b)  O_m(0,b)   with  b  the bin, tau  the imaginary time
!!!                                          and  n  the  orbital index
!!!                   Sign1(b)            =  the  sign  for   each bin          
!!!                   Background(n,b)     =  O_n(b)
!!!                   L_Back              if  true (false)  background is  (not) taken  into account
!!!                   Optional Weights(n) 
!!!       On Output
!!!                  If   Weight  is   not  provided   computes  the  trace
!!!                     Xmean(tau)    =      \sum_n <O^dag_n(tau,:) O_n(0,:)>/<sign(:)>    -     (Background(n,:)>/<sign(:)>)^2
!!!                     Cov  (tau, tau')  =  Covariance matrix              
!!!                  If   Weight  is  provided  
!!!                     Xmean(tau)    =       <M^dag(tau,:) M(0,:)>/<sign(:)>    -     (M_background(:)>/<sign(:)>)^2
!!!                     with   M(tau,b)        =  \sum_n Weights(n) O_n(tau,b)
!!!                     and    M_background(b) =  \sum_n Weights(n) Background(n,b)          
!!!                     Cov  (tau, tau')  =  Covariance matrix
!!!!===================================================
           
           Implicit None
 
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:,:,:), INTENT(IN) ::  GR_IN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)      , INTENT(IN) ::  SIGN1
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:)    , INTENT(IN) ::  Background
           Logical                                      , INTENT(IN) ::  L_Back
           Complex (Kind=Kind(0.d0)), DIMENSION(:)      , INTENT(IN), Optional ::  Weights

           Complex (Kind=Kind(0.d0)), DIMENSION(:)        , INTENT(OUT) ::  XMEAN
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:)      , INTENT(OUT) ::  XCOV


           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP, XMEAN_R
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR, Y, Xhelp, Yhelp
           INTEGER :: NT, NT1, NB, NTDM, NDATA, Nth, Norb, no,  Norb_eff,  no1, no2
           COMPLEX (Kind=Kind(0.d0)) :: Z, tmp
           complex (kind=kind(0.d0)), dimension(:), allocatable :: BGME, BGJ
           complex (kind=kind(0.d0)), dimension(:,:), allocatable :: GR

           NTDM  = SIZE(GR_in,1)
           NDATA = SIZE(GR_in,4)
           Norb  = SIZE(background,1)
           If  (Present(Weights))  then
              Norb_eff = 1
           else
              Norb_eff = Norb
           endif
           allocate ( BGME(Norb_eff), BGJ(Norb_eff) )

           !Write(6,*) 'Errors.F90 ', NTDM, NDATA
           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(error_unit,*) 'Error in COVJS_C'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           ENDIF

           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA), XMEAN_R(NTDM) )
           ALLOCATE( GR(NTDM,NDATA) )
           XMEAN = CMPLX(0.d0, 0.d0, kind(0.D0))
           XCOV  = CMPLX(0.d0, 0.d0, kind(0.D0))
           BGME  = CMPLX(0.d0, 0.d0, kind(0.D0))
           GR    = CMPLX(0.d0,0.d0,kind(0.d0))

           If  (Present(Weights))  then
              do  no1 = 1,Norb
                 do no2  = 1,Norb
                    Z  =   Weights(no2)*  conjg (Weights(no1) ) 
                    do  nt  = 1,NTDM
                       do  nb  =  1,  NDATA
                          GR(nt,nb)   =    GR(nt,nb) +  Z * GR_IN(nt,no1,no2,nb)
                       enddo
                    enddo
                 enddo
              enddo
              do no = 1,Norb
                 Do NB=1,NDATA
                    BGME(1) = BGME(1) + Weights(no)*background(no,NB)
                 End Do
              enddo
           else
              do no = 1,Norb
                 do  nt  = 1,NTDM
                    do  nb  =  1,  NDATA
                       GR(nt,nb)   =    GR(nt,nb) +  GR_IN(nt,no,no,nb)
                    enddo
                 enddo
              enddo
              Do NB=1,NDATA
                 BGME(:) = BGME(:) + background(:,NB)
              End Do
           endif
           
           
           DO NTH = 1,2
              Z = CMPLX(1.D0, 0.D0, kind(0.D0))
              IF (NTH .EQ. 2 ) Z = CMPLX( 0.D0, -1.D0, kind(0.D0))
              DO NT = 1,NTDM
                 Xhelp=0.d0
                 Yhelp=0.d0
                 DO NB= 1, NDATA
                    Xhelp = Xhelp + DBLE ( Z*GR(NT,NB) )
                    Yhelp = Yhelp + SIGN1(NB)
                 ENDDO
                 DO NB= 1, NDATA
                    Y = Yhelp - SIGN1(NB)
                    If  (Present(Weights))  then
                       tmp = cmplx(0.d0,0.d0,kind(0.d0))
                       do  no = 1,Norb
                          tmp =  tmp  +  Weights(no)*background(no,NB)
                       enddo
                       BGJ(1) = BGME(1) - tmp
                    else
                       do no = 1,Norb
                          BGJ(no) = BGME(no) - background(no,NB)
                       enddo
                    endif
                    tmp = cmplx(0.d0,0.d0,kind(0.d0))
                    if (L_Back)  then 
                       DO No = 1,Norb_eff
                          tmp = tmp + BGJ(No)*BGJ(No)/Y
                       enddo
                    endif
                    X = Xhelp - DBLE ( Z*GR(NT,NB) ) - DBLE ( Z*tmp )
                    HLP1(NT,NB) = X/Y
                    HLP (NB   ) = X/Y
                 ENDDO
                 CALL ERRCALC(HLP,XM ,XERR)
                 XMEAN(NT) = XMEAN(NT) + CONJG(Z)*XM
                 XMEAN_R(NT) = XM
                 !if (Nth.eq.2) write(6,*) XM
              ENDDO


              DO NT = 1,NTDM
                 DO NT1= 1,NTDM
                    X = 0.d0
                    DO NB = 1,NDATA
                       X = X +  (HLP1(NT,NB)-XMEAN_R(NT))*(HLP1(NT1,NB)-XMEAN_R(NT1))
                    ENDDO
                    X = X/DBLE(NDATA)
                    XCOV(NT,NT1)  = XCOV(NT,NT1) + CONJG(Z)* X *DBLE(NDATA)
                 ENDDO
              ENDDO
           ENDDO

           DEALLOCATE( HLP, HLP1, XMEAN_R, BGME, BGJ )

           RETURN
         END SUBROUTINE COVJS_C_BG_WEIGHTS

!!!!===================================================
         SUBROUTINE COVJS_C_BG_Rebin_Weights(Gr_in, SIGN, XCOV, XMEAN, background,L_Back, Nrebin,Weights)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Multi-orbital covariance with rebinning, weights, and background.
!
!> @details
!> Combines data rebinning with sign-averaged covariance computation,
!> orbital weighting, and background subtraction for multi-orbital
!> observables. Handles QMC sign problem through jackknife on ⟨s·O⟩/⟨s⟩.
!>
!> **Rebinning step:**
!> Groups consecutive bins (Nrebin bins → 1 effective bin) to reduce
!> autocorrelations before error analysis. All data and signs averaged
!> over Nrebin bins.
!>
!> **Without weights:** Computes orbital trace:
!>   Xmean(τ) = Σₙ ⟨Ô†ₙ(τ)Ôₙ(0)⟩/⟨s⟩ - (⟨O_bg,n⟩/⟨s⟩)²
!>
!> **With weights:** Computes weighted sum:
!>   M(τ,b) = Σₙ Weights(n) × Oₙ(τ,b)
!>   Xmean(τ) = ⟨M†(τ)M(0)⟩/⟨s⟩ - (⟨M_bg⟩/⟨s⟩)²
!>
!> For Nrebin=1, directly calls COVJS_C_BG_Weights without rebinning.
!
!> @param[in] Gr_in Multi-orbital data (Ntimes × Norb × Norb × Nbins)
!> @param[in] SIGN Configuration signs (Nbins)
!> @param[out] XCOV Covariance matrix (Ntimes × Ntimes), complex
!> @param[out] XMEAN Background-subtracted mean values (Ntimes), complex
!> @param[in] background Orbital backgrounds (Norb × Nbins), complex
!> @param[in] L_Back Enable/disable background subtraction
!> @param[in] Nrebin Number of bins to combine (typically 1-10)
!> @param[in] Weights [optional] Orbital weights (Norb), complex
!
!> @warning Nbins must be divisible by Nrebin
!
!> @see COVJS_C_BG_Weights, COVJS_C_REBIN_BG
!--------------------------------------------------------------------
!!!!===================================================
!!!        On Input   Gr(tau,n,m,b)       =  O^dag_n(tau,b)  O_m(0,b)   with  b  the bin, tau  the imaginary time
!!!                                          and  n  the  orbital index
!!!                   Sign1(b)            =  the  sign  for   each bin          
!!!                   Background(n,b)     =  O_n(b)
!!!                   L_Back              if  true (false)  background is  (not) taken  into account
!!!                   N_rebin             Rebinning 
!!!                   Optional Weights(n) 
!!!       On Output
!!!                  If   Weight  is   not  provided   computes  the  trace
!!!                     Xmean(tau)    =      \sum_n <O^dag_n(tau,:) O_n(0,:)>/<sign(:)>    -     (Background(n,:)>/<sign(:)>)^2
!!!                     Cov  (tau, tau')  =  Covariance matrix              
!!!                  If   Weight  is  provided  
!!!                     Xmean(tau)    =       <M^dag(tau,:) M(0,:)>/<sign(:)>    -     (M_background(:)>/<sign(:)>)^2
!!!                     with   M(tau,b)        =  \sum_n Weights(n) O_n(tau,b)
!!!                     and    M_background(b) =  \sum_n Weights(n) Background(n,b)          
!!!                     Cov  (tau, tau')  =  Covariance matrix
!!!!===================================================
           
           Implicit None
 
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:,:,:), INTENT(IN) ::  GR_IN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)      , INTENT(IN) ::  SIGN
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:)    , INTENT(IN) ::  Background
           Integer                                      , INTENT(IN) ::  Nrebin
           Logical                                      , INTENT(IN) ::  L_Back
           Complex (Kind=Kind(0.d0)), DIMENSION(:)      , INTENT(IN), Optional ::  Weights

           Complex (Kind=Kind(0.d0)), DIMENSION(:)        , INTENT(OUT) ::  XMEAN
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:)      , INTENT(OUT) ::  XCOV

           
           
           !Local
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:,:,:), allocatable ::  GR1
           Real    (Kind=Kind(0.d0)), DIMENSION(:)      , allocatable ::  SIGN1
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:)    , allocatable ::  Background1

           Integer ::  Norb, Ntdm, Ndata,  Ndata1, no1,no2,no, N, nc,nb, Nt

           If (NREBIN == 1 )  then
              Call COVJS_C_BG_Weights(Gr_in, SIGN, XCOV, XMEAN, background,L_Back, Weights)
           else

              NTDM  = SIZE(GR_IN,1)
              NDATA = SIZE(GR_IN,4)
              Norb  = SIZE(background,1)
              NDATA1 = NDATA/NREBIN
              ALLOCATE ( GR1(NTDM,Norb,Norb,NDATA1), SIGN1(NDATA1), background1(Norb,NDATA1) )
              
              SIGN1 = 0.d0
              GR1   = CMPLX(0.d0,0.d0,kind(0.d0))
              background1 = CMPLX(0.d0,0.d0,kind(0.d0))
              ! Rebin
              NC = 0
              DO N = 1,NDATA1
                 DO NB = 1,NREBIN
                    NC = NC + 1
                    SIGN1(N) = SIGN1(N) + SIGN(NC)
                    do no1  = 1,Norb
                       Do  no2  = 1,Norb
                          DO NT = 1,NTDM
                             GR1(NT,no1,no2,N)  = GR1(NT,no1,no2,N) + GR_in(NT,no1,no2,NC)
                          ENDDO
                       enddo
                    enddo
                    DO No = 1,Norb
                       background1(No,N) = background1(No,N) + background(No,NC)
                    enddo
                 ENDDO
              ENDDO
              SIGN1 = SIGN1/DBLE (NREBIN)
              GR1   = GR1  /CMPLX(DBLE(NREBIN),0.d0,KIND(0.d0))
              background1 =  background1 /CMPLX(DBLE(NREBIN),0.d0,KIND(0.d0))

              Call COVJS_C_BG_Weights(Gr1, SIGN1, XCOV, XMEAN, background1,L_Back, Weights)
              DEALLOCATE ( GR1, SIGN1, background1 )
           endif

         END SUBROUTINE COVJS_C_BG_REBIN_WEIGHTS

!========================

         SUBROUTINE COVJS_C_REBIN_BG(GR, SIGN2, XCOV, XMEAN,NREBIN, background)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Rebinned covariance with sign problem and background subtraction.
!
!> @details
!> Combines data rebinning with sign-averaged covariance computation
!> and background subtraction. Wrapper around COVJS_C_BG that first
!> rebins data to reduce autocorrelation effects before jackknife
!> error analysis with sign problem treatment.
!>
!> **Rebinning procedure:**
!> Groups consecutive bins (Nrebin bins → 1 effective bin) by
!> averaging all observables and signs over Nrebin consecutive bins.
!> Reduces effective number of bins: Nbins_eff = Nbins / Nrebin.
!>
!> **Error analysis:**
!> After rebinning, computes via COVJS_C_BG:
!>   Xmean(τ) = Trace[Ĝ(τ)] - (⟨O_bg⟩/⟨s⟩)²
!> using jackknife method on rebinned data.
!>
!> Input format simpler than multi-orbital variant: 2D array GR(τ,b)
!> with pre-computed traces or single-orbital data.
!
!> @param[in] GR Observable data (Ntimes × Nbins), complex
!> @param[in] SIGN2 Configuration signs (Nbins)
!> @param[out] XCOV Covariance matrix (Ntimes × Ntimes), complex
!> @param[out] XMEAN Background-subtracted mean values (Ntimes), complex
!> @param[in] NREBIN Number of bins to combine (typically 1-10)
!> @param[in] background Orbital backgrounds (Norb × Nbins), complex
!
!> @warning Nbins must be divisible by Nrebin
!
!> @see COVJS_C_BG, COVJS_C_BG_Rebin_Weights
!--------------------------------------------------------------------

           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times.
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV, background
           Complex (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)   ::  SIGN2
           INTEGER :: NREBIN


           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE ::  GR1, background1
           REAL    (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  SIGN1

           INTEGER :: NTDM, NDATA, NDATA1, N, NB, NC, NT, No, Norb

           NTDM  = SIZE(GR,1)
           NDATA = SIZE(GR,2)
           Norb  = SIZE(background,1)

           NDATA1 = NDATA/NREBIN
           ALLOCATE ( GR1(NTDM,NDATA1), SIGN1(NDATA1), background1(Norb,NDATA1) )

           SIGN1 = 0.d0
           GR1   = CMPLX(0.d0,0.d0,kind(0.d0))
           background1 = CMPLX(0.d0,0.d0,kind(0.d0))
           ! Rebin
           NC = 0
           DO N = 1,NDATA1
              DO NB = 1,NREBIN
                 NC = NC + 1
                 SIGN1(N) = SIGN1(N) + SIGN2(NC)
                 DO NT = 1,NTDM
                    GR1(NT,N)  = GR1(NT,N) + GR(NT,NC)
                 ENDDO
                 DO No = 1,Norb
                    background1(No,N) = background1(No,N) + background(No,NC)
                 enddo
              ENDDO
           ENDDO
           SIGN1 = SIGN1/DBLE (NREBIN)
           GR1   = GR1  /CMPLX(DBLE(NREBIN),0.d0,KIND(0.d0))

           CALL COVJS_C_BG(GR1, SIGN1, XCOV, XMEAN, background1)

           DEALLOCATE ( GR1, SIGN1, background1 )

         END SUBROUTINE COVJS_C_REBIN_BG


         Subroutine COV_ERR(XMEAN, XCOV, ISEED)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Generate random sample consistent with mean and covariance matrix.
!
!> @details
!> Creates synthetic data point by random sampling from a multivariate
!> normal distribution with given mean vector and covariance matrix.
!>
!> Algorithm:
!> 1. Diagonalize covariance matrix: XCOV = U·Λ·Uᵀ
!> 2. Transform mean to eigenspace: μ' = Uᵀ·μ
!> 3. Sample independent Gaussians: x'_i = μ'_i + √λ_i × N(0,1)
!> 4. Transform back: x = U·x'
!>
!> Used for error propagation in fits and Monte Carlo uncertainty
!> quantification.
!
!> @param[in,out] XMEAN On input: mean values. On output: random sample
!> @param[in] XCOV Covariance matrix (Ntimes × Ntimes)
!> @param[in,out] ISEED Random number seed (updated)
!
!> @warning
!> Negative eigenvalues (numerical errors) are replaced with abs(eigenvalue).
!
!> @see COVJ, COVJS (for computing covariance matrices)
!--------------------------------------------------------------------
           !  Given Mean and Cov, diagonalizes the COV and produces a new data set within
           !  the errorbars

           Implicit None
           ! Parameters
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: XCOV
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   :: XMEAN

           Integer :: ntau, I, M, ISeed, ISEED_VEC(1)
           Real (Kind=Kind(0.d0)) :: X

           Real (Kind=Kind(0.d0)), Dimension(:,:),  allocatable ::  UC
           Real (Kind=Kind(0.d0)), Dimension(:),    allocatable ::  XMEAN_1, SIG_1

           ntau = size(Xmean,1)                    ! Dimension of data
           Allocate (UC(ntau,ntau), XMEAN_1(ntau), SIG_1(ntau) )

           ISEED_VEC(1) = ISEED
           CALL RANSET(ISEED_VEC)

           ! Diagonalize covariance: XCOV = UC·Λ·UCᵀ (UC = eigenvectors, SIG_1 = eigenvalues)
           CALL DIAG(XCOV,UC,SIG_1)

           ! Transform mean vector to eigenspace: μ' = UCᵀ·μ
           DO I = 1,NTAU
              X = 0.D0
              DO M = 1,NTAU
                 X  = X + UC(M,I)* XMEAN(M)       ! Matrix-vector product
              ENDDO
              XMEAN_1(I) = X                      ! Mean in eigenspace
           ENDDO
           
           ! Sample independent Gaussians: x'_i ~ N(μ'_i, λ_i)
           DO I = 1,NTAU
              IF (SIG_1(I).LT.0.d0) Then
                  write(6,*) 'Warning in Cov_err: negative eigenvalue', SIG_1(I)
              Endif
              ! Add Gaussian noise scaled by sqrt(eigenvalue)
              XMEAN_1(I) = XMEAN_1(I) + SQRT(ABS(SIG_1(I)))*RANG_WRAP()
           ENDDO
           
           ! Transform random sample back to original space: x = UC·x'
           DO I = 1,NTAU
              X = 0.D0
              DO M = 1,NTAU
                 X  = X + UC(I,M)*XMEAN_1(M)      ! Inverse transform
              ENDDO
              XMEAN(I) = X                        ! Random sample in original space
           ENDDO

           CALL RANGET(ISEED_VEC)
           ISEED  = ISEED_VEC(1)

           Deallocate (UC, XMEAN_1, SIG_1)


         END Subroutine COV_ERR
         

         SUBROUTINE  AUTO_COR(DATA,RES)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Compute normalized autocorrelation function.
!
!> @details
!> Calculates the autocorrelation function:
!> C(τ) = ⟨(x_t - μ)(x_{t+τ} - μ)⟩ / ⟨(x_t - μ)²⟩
!>
!> Essential for determining:
!> - Integrated autocorrelation time: τ_int = 0.5 + Σ_t C(t)
!> - Statistical inefficiency: effective sample size = N / (2τ_int)
!> - Appropriate binning/rebinning sizes
!>
!> For QMC, τ_int determines how many sweeps are needed between
!> independent measurements.
!
!> @param[in] DATA Time series data (typically QMC measurements)
!> @param[out] RES Normalized autocorrelation function C(τ)
!
!> @note
!> - RES(1) = 1.0 by construction
!> - C(τ) decays to 0 for large τ if data is stationary
!> - SIZE(DATA) must be >= SIZE(RES)
!
!> @warning
!> For accurate results, use at least 10×τ_int data points.
!
!> @see ERRCALC_J_REBIN (uses τ_int for rebinning)
!--------------------------------------------------------------------

           Implicit none

           REAL (Kind=Kind(0.d0)),  DIMENSION(:)  :: DATA,RES

           !Local
           Integer  :: nb, nt, ntau, nt1
           Real (Kind=Kind(0.d0)) :: X1, X2, X3

           nb = SIZE(DATA)                         ! Number of data points
           nt = SIZE(RES)                          ! Number of lag times
           if (nb.lt.nt) then
              write(error_unit,*) 'Error in auto_cor: DATA too short for requested lags'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           end if

           ! Compute autocorrelation for each lag time
           DO ntau = 1,  nt
              X1 = 0.0                             ! Numerator: ⟨x(t)x(t+τ)⟩
              X2 = 0.0                             ! Denominator: ⟨x(t)²⟩
              X3 = 0.0                             ! Mean
              
              ! Compute mean for this lag
              DO nt1 = 1, nb - ntau
                 X3 = X3 + DATA(nt1)
              ENDDO
              X3 = X3 / dble(nb - ntau)

              ! Compute covariance and variance
              DO nt1 = 1, nb - ntau
                 X1 = X1 + (DATA(nt1)-x3)*(DATA(nt1 + ntau)-x3)  ! Covariance
                 X2 = X2 + (DATA(nt1)-x3)*(DATA(nt1)-x3)         ! Variance
              ENDDO
              X1 = X1 / dble(nb - ntau)
              X2 = X2 / dble(nb - ntau)

              Res(ntau)  = X1/X2                   ! Normalized autocorrelation

           ENDDO

         END SUBROUTINE AUTO_COR

         SUBROUTINE BOOTSTRAPC_FLUC(A,B,AB,NBOOT,ISEED,ZM,ZERR)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Bootstrap error estimation for connected correlation functions.
!
!> @details
!> Computes bootstrap error on the connected correlator:
!> G_connected = ⟨AB⟩ - ⟨A⟩⟨B⟩
!>
!> This is essential for computing fluctuations, susceptibilities, and
!> connected Green's functions in QMC. Bootstrap properly accounts for
!> correlations between A, B, and AB.
!>
!> Algorithm:
!> 1. Generate NBOOT bootstrap samples by random resampling with replacement
!> 2. For each sample, compute: G = ⟨AB⟩_boot - ⟨A⟩_boot⟨B⟩_boot
!> 3. Calculate mean and standard deviation of G across bootstrap samples
!
!> @param[in] A Observable A bins
!> @param[in] B Observable B bins
!> @param[in] AB Product observable AB bins (same bins as A and B)
!> @param[in] NBOOT Number of bootstrap samples (typically 1000-10000)
!> @param[in,out] ISEED Random seed (updated on exit)
!> @param[out] ZM Mean value of connected correlator
!> @param[out] ZERR Bootstrap error estimate
!
!> @note
!> All input arrays must have the same size (number of bins).
!>
!> @see Bootstrap (simpler version without fluctuations)
!--------------------------------------------------------------------
           !!!  COMPUTES <AB> - <A><B>
           IMPLICIT NONE
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), INTENT(IN) :: A,B,AB
           INTEGER, INTENT(IN)    :: NBOOT
           INTEGER, INTENT(INOUT) :: ISEED
           COMPLEX (Kind=Kind(0.d0)), INTENT(OUT) :: ZM,ZERR


           !Local
           INTEGER :: NP, NB, I, J, ISEED_VEC(1)
           COMPLEX (Kind=Kind(0.d0)) :: Z,  Z1,Z2,Z12


           ISEED_VEC(1) = ISEED
           CALL RANSET(Iseed_vec)
           
           NP = SIZE(A,1)                          ! Number of bins
           ZM   = CMPLX(0.d0,0.d0,Kind=Kind(0.d0))
           ZERR = CMPLX(0.d0,0.d0,Kind=Kind(0.d0))
           
           ! Generate NBOOT bootstrap samples
           DO NB = 1, NBOOT
              Z1  = cmplx(0.d0,0.d0,Kind=Kind(0.d0))  ! Bootstrap mean of A
              Z2  = cmplx(0.d0,0.d0,Kind=Kind(0.d0))  ! Bootstrap mean of B
              Z12 = cmplx(0.d0,0.d0,Kind=Kind(0.d0))  ! Bootstrap mean of AB
              
              ! Resample with replacement
              DO I = 1,NP
                 J = NINT( DBLE(NP)* RANF_WRAP() + 0.5D0 )
                 IF (J == 0) J = 1
                 IF (J > NP) J = NP
                 Z1 = Z1  + A(J)
                 Z2 = Z2  + B(J)
                 Z12 =Z12 + AB(J)
              ENDDO
              Z1 = Z1 /CMPLX(DBLE(NP),0.d0,Kind=Kind(0.d0))
              Z2 = Z2 /CMPLX(DBLE(NP),0.d0,Kind=Kind(0.d0))
              Z12 =Z12/CMPLX(DBLE(NP),0.d0,Kind=Kind(0.d0))

              Z    = Z12 - Z1*Z2                   ! Connected correlator for this sample
              ZM   = ZM   + Z
              ZERR = ZERR + Z*Z
           ENDDO
           
           ! Compute mean and error across bootstrap samples
           ZM   = ZM  /CMPLX(DBLE(NBOOT),0.d0,Kind=Kind(0.d0))
           ZERR = ZERR/CMPLX(DBLE(NBOOT),0.d0,Kind=Kind(0.d0))

           Z = ZERR -  ZM*ZM                       ! Variance
           ZERR = SQRT(Z)                          ! Standard deviation

           CALL RANGET(Iseed_vec)
           ISEED = ISEED_VEC(1)


         END SUBROUTINE BOOTSTRAPC_FLUC

         SUBROUTINE BOOTSTRAP(EN,XM,XERR,NBOOT,ISEED)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Bootstrap resampling for error estimation.
!
!> @details
!> Estimates error using bootstrap resampling (Efron, 1979):
!> 1. Generate NBOOT samples by randomly drawing N points from data
!>    with replacement
!> 2. Calculate mean for each bootstrap sample
!> 3. Compute standard deviation across bootstrap means
!>
!> Bootstrap advantages:
!> - No assumptions about data distribution
!> - Works for any statistic (not just means)
!> - Provides confidence intervals
!>
!> Bootstrap disadvantages for QMC:
!> - Slightly more biased than jackknife
!> - Computationally more expensive
!>
!> For QMC data, jackknife (ERRCALCJ) is generally preferred.
!
!> @param[in] EN Input bins/samples
!> @param[out] XM Bootstrap mean estimate
!> @param[out] XERR Bootstrap error estimate
!> @param[in] NBOOT Number of bootstrap samples (1000-10000 typical)
!> @param[in,out] ISEED Random number seed (updated on exit)
!
!> @note
!> For correlated data, first rebin to reduce autocorrelation.
!
!> @see ERRCALCJ (jackknife - preferred for QMC), BOOTSTRAPC_FLUC
!--------------------------------------------------------------------

           IMPLICIT NONE
           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           REAL (Kind=Kind(0.d0))               ::  XM, XERR,  X
           INTEGER                     ::  NP, NT, NBOOT, NB, I, ISEED

           ! Local
           INTEGER :: ISEED_VEC(1)

           NP = SIZE(EN)                           ! Number of bins
           ISEED_VEC(1) = ISEED
           CALL RANSET(Iseed_vec)


           ! Generate NBOOT bootstrap samples
           XM   = 0.D0
           XERR = 0.D0
           DO NB = 1,NBOOT
              X = 0.D0
              ! Resample with replacement
              DO NT = 1, NP
                 I = NINT( DBLE(NP)* RANF_WRAP() + 0.5D0 )
                 IF (I.EQ.0 .OR. I.GT.NP ) THEN
                    WRITE(error_unit,*) 'ERROR IN BOOTSTRAP: invalid random index'
                    Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
                 ENDIF
                 X = X + EN(I)                     ! Add randomly selected bin
              ENDDO
              X = X/DBLE(NP)                       ! Mean of this bootstrap sample
              XM   = XM + X
              XERR = XERR + X*X
           ENDDO

           ! Compute mean and error across bootstrap samples
           XM   = XM  /DBLE(NBOOT)
           XERR = XERR/DBLE(NBOOT)

           X = XERR - XM*XM                        ! Variance
           XERR = 0.d0
           IF (X.GT.0.d0) XERR = SQRT(X)           ! Standard deviation

           CALL RANGET(Iseed_vec)
           ISEED = ISEED_VEC(1)

         END SUBROUTINE BOOTSTRAP

       END MODULE ERRORS
