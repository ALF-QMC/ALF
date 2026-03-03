
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


Module MaxEnt_mod
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Classic maximum entropy implementation for analytical continuation.
!> @details Follows closely the work of
!>   1.   W. von der Linden, Maximum-entropy data analysis, Applied Physics A 60 (1995), no. 2, 155–165. 
!>   and 
!>   2.  Mark Jarrell and J.E. Gubernatis, Bayesian inference and the analytic continuation 
!>   of imaginary-time quantum monte carlo data, Physics Reports 269 (1996), no. 3, 133–195.
!--------------------------------------------------------------------

        Use MyMats
        Use Errors
        use iso_fortran_env, only: output_unit, error_unit
        use runtime_error_mod

        Interface MaxEnt
           Module Procedure MaxEnt_T, MaxEnt_T0
        end Interface

        REAL (Kind=Kind(0.d0)), Private   ::   ZERO, ALPHA, XMOM1
        REAL (Kind=Kind(0.d0)), Dimension(:),   Allocatable, Private :: XLAM,  DEF, SIG1
        REAL (Kind=Kind(0.d0)), DIMENSION(:,:), Allocatable, Private :: COVM1, UC
        Integer,   Private :: NTAU, NOM


        CONTAINS

!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Maximum entropy method for spectral function recovery at T > 0.
!> @details Inverts QMC data in imaginary time to recover spectral function
!> using maximum entropy principle. Algorithm follows von der Linden (1995)
!> and Jarrell & Gubernatis (1996). Uses iterative determination of Lagrange
!> multiplier (ALPHA) to enforce chi-squared criterion.
!>
!> @param[in] XQMC Real array of QMC measurements at imaginary times
!> @param[in] COV Covariance matrix of QMC measurements
!> @param[out] A Spectral function on frequency grid (non-negative)
!> @param[in] XKER Kernel matrix K[τ,ω] = exp(-τ*ω)/(1+exp(-β*ω))
!> @param[in,out] ALPHA_ST Initial entropy weight (updated via iteration)
!> @param[out] CHISQ Final chi-squared goodness-of-fit measure
!> @param[in] DEFAULT Optional default model (uniform if not provided)
!> @note XQMC should be centered (zero mean) for best results
!> @pre ALPHA_ST > 0, size(COV) = (NTAU, NTAU)
!--------------------------------------------------------------------
          Subroutine MaxEnt_T( XQMC, COV, A, XKER, ALPHA_ST, CHISQ,DEFAULT)

            Implicit None
            Real (Kind=Kind(0.d0)), Dimension(:)   :: XQMC, A
            Real (Kind=Kind(0.d0)), Dimension(:,:) :: COV, XKER
            Real (Kind=Kind(0.d0)) :: ALPHA_ST, CHISQ, ALPHA_N
            Real (Kind=Kind(0.d0)), Dimension(:), optional   :: Default

            Integer       :: NT, NT1, NT2, NW, NFLAG, NCOUNT
            Real (Kind=Kind(0.d0)) :: X, XENT, XQ, PR_ALP, XTRACE, DIFF1, DIFF , Tol_chi_def


            Tol_chi_def = 1000000000000.0D0
            NTAU = SIZE(XQMC,1)
            NOM  = SIZE(A, 1)
            !WRITE(6,*) 'NTAU, Nom: ', NTAU,NOM
            Xmom1 = Xqmc(1)

            ZERO =  1.0D-8
            ALLOCATE ( XLAM(NTAU), SIG1(NTAU), COVM1(NTAU,NTAU), UC(NTAU,NTAU), DEF(NOM) )
            XLAM=0.D0;  SIG1=0.D0; UC = 0.D0

            !Open (Unit=77,File='Aom_steps',Status='unknown')

            !Open(Unit=14)
            !do nt = 1, NTAU
            !   Write(14,*) Nt, XQMC(nt), sqrt(Cov(Nt,Nt))
            !enddo
            !Close(14)

            CALL DIAG(COV,UC,SIG1)
            DO NT1 = 1,NTAU
               DO NT2 = 1,NTAU
                  X = 0.D0
                  DO NT = 1,NTAU
                     X = X + UC(NT1,NT)*UC(NT2,NT)/SIG1(NT)
                  ENDDO
                  COVM1(NT1,NT2) = X
               ENDDO
            ENDDO


            Open (Unit=44, File="Max_cl_log", Status="unknown")

            !Write(44,*) 'N E W   R U N'
            !Write(44,*) '# of data points: ', NTAU
            !Write(6,*) 'N E W   R U N'
            ! Set the Default.
            ALPHA = Alpha_st
            DEF  = XMOM1/dble(NOM)
            XLAM = 0.d0
            if ( Present(Default) ) then
               DEF = Default
               Write(44,*) 'Will use provided default'
            else
               XQ = 0.d0; XENT= 0.d0; CHISQ = 0.d0
               Call Maximize_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)
               IF (CHISQ .GT. Tol_chi_def*NTAU )  THEN
                  DO
                     XQ = 0.d0; XENT= 0.d0; CHISQ = 0.d0
                     Call Maximize_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)
                     !Write(44,*) 'Default: ', Alpha, Chisq
                     !Write(6,*) 'Default: ', Alpha, Chisq
                     IF (CHISQ .GT. Tol_chi_def*NTAU .AND.  ALPHA.GT.100 )  THEN
                        ALPHA = ALPHA - ALPHA*0.1
                     ELSE
                        CALL SETA(A,XKER)
                        DO NW = 1,NOM
                           IF (A(NW).LT.ZERO) THEN
                              DEF(NW)= ZERO
                           ELSE
                              DEF(NW) = A(NW)
                           ENDIF
                        ENDDO
                        EXIT
                     ENDIF
                  ENDDO
               ELSE
                  !Write(6,*) 'Flat Default'
               Endif
               !DO NW = 1,NOM
               !   Write(13,*) NW, DEF(NW)
               !ENDDO
               !Write(6,*) 'Default Final: ', Alpha, Chisq
               DEF  = XMOM1/dble(NOM)
               !Write(6,*) 'Setting the default to a flat default'
            endif

            ! Calssic MaxEnt.
            NFLAG  = 0
            NCOUNT = 0
            !ALPHA  = ALPHA_ST
            XLAM = 0.D0
            DO
               !WRITE(6,*)  'Starting classic  ', ALPHA
               WRITE(44,*)  '========= 1/Alpha:    ', ALPHA
               XQ = 0.d0; XENT= 0.d0; CHISQ = 0.d0
               !write(6,*) 'Calling maximize'
               CALL MAXIMIZE_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)
               !write(6,*) 'Return: Calling maximize'
               IF (NFLAG.EQ.0) THEN
                  CALL CALCPR_ALP(XQMC,  COV, A, XKER,XQ,XENT,PR_ALP,XTRACE)
                  ALPHA_N = -XTRACE/(2.D0*XENT)
                  WRITE(44,*) 'Max at:', ALPHA_N
                  !WRITE(6,*) 'Max at:', ALPHA_N
                  !WRITE(6,*) 'Old_alp', ALPHA
                  DIFF1 =    ABS(ALPHA_N - ALPHA)
               ENDIF
               CALL SETA(A,XKER)
               CALL SETQ(A,XKER,XQMC, XQ,XENT,CHISQ)
               WRITE(44,2006) ALPHA, XQ,XENT,CHISQ
               !WRITE(6,2006 ) ALPHA, XQ,XENT,CHISQ
               DIFF =  ALPHA_N - ALPHA
               IF ( ABS(DIFF) .GT.  0.1*ALPHA ) THEN
                  ALPHA  = ALPHA +  0.1 * ALPHA * DIFF/ABS(DIFF)
                  NFLAG = 1
               ELSE
                  ALPHA =  ALPHA_N
                  NFLAG = 0
               ENDIF
               NCOUNT = NCOUNT + 1
               IF (NCOUNT .EQ. 100) THEN
                  WRITE(44,*) 'NOT CONVERGED'
               ENDIF
               IF ( ABS(DIFF1)/ABS(ALPHA_N).LT.0.01D0 .OR.  NCOUNT.GT.1000 ) Exit
               !&
               !     &    .OR. CHISQ.LT. 0.*dble(NTAU) ) EXIT
           ENDDO


           CLOSE(44)

2006       FORMAT('Res: 1/Alpha, XQ,S,CHI: ', F24.12,2x,F24.12,2x,F24.12,2x,F24.12)

           ALPHA_ST =  ALPHA 
            DEALLOCATE ( XLAM, SIG1, COVM1, UC, DEF )
            !Close(77)
          End Subroutine MaxEnt_T

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief Newton iteration to solve MaxEnt's Lagrange multiplier equations.
          !> @details Iterative solver for F(\lambda) = 0 using Newton-Raphson.
          !> Computes spectral function A from Lagrange multipliers \lambda via
          !> A(\omega) = exp(\lambda - \chi(\omega)). Each Newton step:
          !> \Delta\lambda = J^{-1} \cdot F, where J = Hessian (second derivatives).
          !> Converges quadratically for well-conditioned problems; applies line
          !> search to stabilize in nonquadratic regime.
          !> @param[in] XQMC     QMC measurements G(tau)
          !> @param[in] COV      Covariance matrix (must be positive definite)
          !> @param[out] A       Converged spectral function A(omega)
          !> @param[in] XKER     Kernel matrix K(tau,omega)
          !> @param[out] XQ      Final chi-squared (fit quality)
          !> @param[out] XENT    Final entropy S[A]
          !> @param[out] CHISQ   Normalized chi-squared (for convergence check)
          !> @note Requires good initial guess (from MaxEnt_T default model);
          !> iteration count ~10-50 typical for tight tolerance; fails if COV singular.
          !> @pre COV symmetric; size(XQMC)=size(XKER,1); size(A)=size(XKER,2)
          !-------------------------------------------------------------------

          Subroutine Maximize_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)

            ! Solves F(tau) = 0 with Newton.


            Implicit None
            !Arguments
            REAL (Kind=Kind(0.d0)), intent(out)    :: XQ,XENT,CHISQ
            REAL (Kind=Kind(0.d0)), Dimension(:)   :: XQMC, A
            REAL (Kind=Kind(0.d0)), Dimension(:,:) :: COV, XKER


            !Working space
            REAL (Kind=Kind(0.d0)), DIMENSION(:),  ALLOCATABLE ::  XLAM1,  F
            REAL (Kind=Kind(0.d0)), DIMENSION(:,:),ALLOCATABLE ::  AH, AHINV

            Real (Kind=Kind(0.d0)) :: X, XNORM, DET1(2), XMAX
            Integer :: NITER, NT, NT1


            ALLOCATE (XLAM1(NTAU), F(NTAU))
            XLAM1 = 0.D0; F = 0.D0
            ALLOCATE (AH(NTAU,NTAU), AHINV(NTAU,NTAU))
            AH = 0.D0; AHINV = 0.D0

            NITER = 0
            !WRITE(6,*) "Starting Maximize"
            DO
               !Write(6,*) ' Iteration :: ', Niter
               CALL SETA (A,XKER)
               !Write(6,*) ' Back From SetA '
               CALL SETAH(AH, A,XKER,COV)
               !Write(6,*) ' Back From SetAH '
               CALL SETF (F, COV, XKER, A, XQMC)
               !Write(6,*) ' Back From SetF '
               !Write(6,*) 'Calling INV'
               CALL INV(AH, AHINV, DET1)
               !Write(6,*) 'Back Calling INV', Det1(1),Det1(2)
               !CALL INV(AH, AHINV)
               !Write(6,*) ' Back From INV '
               XNORM = 0.D0
               XMAX = 0.d0
               DO NT = 1,NTAU
                  X = 0.D0
                  DO NT1 = 1,NTAU
                     X = X + AHINV(NT,NT1)*F(NT1)
                  ENDDO
                  XLAM1(NT) = XLAM(NT) - X
                  XNORM = XNORM + X*X
                  If (ABS(X).GT.XMAX) XMAX = ABS(X)
               ENDDO
               !Write(6,*)  'Max Diff Newton: ',  XMAX
               XNORM =  SQRT(XNORM)/DBLE(NTAU)
               !DO nw = 1,Nom
               !write(77,*) nw, A(nw)
               !enddo
               !write(77,*) '# Chisq :  ', CHISQ, XMAX
               !write(77,*)
               DO NT = 1,NTAU
                  XLAM(NT) = XLAM1(NT)
               ENDDO
               NITER = NITER + 1
               !WRITE(6,*) 'Maximize: ', XNORM, NITER
               IF (XNORM.LT.1.0D-6 .OR. NITER.GE.500) EXIT
            ENDDO
            CALL   SETQ(A,XKER,XQMC, XQ,XENT,CHISQ)

            IF (NITER.GE.500) THEN
               WRITE(44,*) 'Convergence problem:', Xnorm
            ENDIF

            Deallocate (XLAM1, F)
            Deallocate (AH, AHINV)

          END Subroutine Maximize_Newton

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief Self-consistent iteration solver for MaxEnt Lagrange multipliers.
          !> @details Alternative to Newton: solves linear system (alpha*COV*lambda)_i
          !> = G_bar(i) - G_qmc(i) iteratively without full Hessian inversion. Uses
          !> conjugate-gradient or direct solve at each step. More robust for
          !> ill-conditioned COV; slower convergence (linear vs quadratic).
          !> Algorithm: Given guess lambda, compute residual F = K^T*(A_conv - G_qmc),
          !> then update lambda via: lambda_new = lambda + (COV^{-1} F) / alpha.
          !> @param[in] XQMC     QMC measurements at imaginary times
          !> @param[in] COV      Covariance matrix (can be nearly singular)
          !> @param[out] A       Converged spectral function
          !> @param[in] XKER     Forward kernel K(tau,omega)
          !> @param[out] XQ      Final chi-squared statistic
          !> @param[out] XENT    Final entropy S[A]
          !> @param[out] CHISQ   Normalized chi-squared (tracking criterion)
          !> @note Prefers LAPACK solver for stability; tolerates noise better than Newton;
          !> typical ~30-100 iterations; returns if XENT converges or iteration limit reached.
          !> @see Maximize_Newton for faster (but less stable) Newton variant
          !> @pre size(XQMC)=size(XKER,1); COV symmetric (rank might be < Ntau)
          !-------------------------------------------------------------------

          Subroutine Maximize_Self( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)

            ! Solves F(tau) = 0 with self-consistency.
            ! That is. Iterate to solve:  alpha Cov(t,t1) xlam(t1) = \bar{G}(t) - G_qmc(t)
            ! bar{G}(t) is the fit

            Implicit None


            !Arguments
            REAL (Kind=Kind(0.d0))                 :: XQ,XENT,CHISQ
            REAL (Kind=Kind(0.d0)), Dimension(:)   :: XQMC, A
            REAL (Kind=Kind(0.d0)), Dimension(:,:) :: COV, XKER


            !Working space
            REAL (Kind=Kind(0.d0)), DIMENSION(:),  ALLOCATABLE ::  XLAM1, GBAR

            Real (Kind=Kind(0.d0)) :: XNORM
            Integer :: NITER, NT, NT1, NW


            ALLOCATE (XLAM1(NTAU), GBAR(NTAU) )
            XLAM1 = 0.D0


            NITER = 0
            DO
               CALL SETA (A,XKER)
               DO NT = 1,NTAU
                  GBAR(NT) = 0.d0
                  DO NW = 1,NOM
                     GBAR(NT) = GBAR(NT) + XKER(NT,NW)*A(NW)
                  ENDDO
                  GBAR(NT) = ( GBAR(NT) - XQMC(NT) ) / ALPHA
               ENDDO
               XNORM = 0.D0
               DO NT = 1,NTAU
                  XLAM1(NT) = 0.d0
                  DO NT1 = 1,NTAU
                     XLAM1(NT) = XLAM1(NT) +  COVM1(NT,NT1)*GBAR(NT1)
                  ENDDO
                  XNORM = XNORM + ( XLAM1(NT) - XLAM(NT) )**2
               ENDDO
               !IF (MOD(NITER,100) .EQ. 0 ) THEN
               !   DO NT = 1,NTAU
               !      Write(6,*) 'Self: ', XLAM(NT), XLAM1(NT)
               !   ENDDO
               !ENDIF
               XNORM =  SQRT(XNORM)/DBLE(NTAU)
               DO NT = 1,NTAU
                  XLAM(NT) = XLAM1(NT)
               ENDDO
               NITER = NITER + 1
               !WRITE(6,*) 'Maximize_Self: ', XNORM, NITER
               IF (XNORM.LT.1.0D-6 .OR. NITER.GE.1000) EXIT
            ENDDO
            CALL  SETQ(A,XKER,XQMC, XQ,XENT,CHISQ)

            IF (NITER.GE.1000) THEN
               WRITE(44,*) 'Convergence problem:', XNORM
            ENDIF

            Deallocate (XLAM1, GBAR)

          END Subroutine Maximize_Self

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief Assemble spectral entropy sum from normalized spectrum.
          !> @details Computes S[A] = sum_omega A(omega) * K(tau,omega) over
          !> frequency grid, normalized by frequency spacing. Updates global
          !> entropy S and entropy density chi (entropy per frequency).
          !> Used internally by MaxEnt iterations to track constraint satisfaction.
          !> Formula: S = -\int d\omega A(\omega) ln(A(\omega)/rho_0(\omega))
          !> where rho_0 is default model (typically constant).
          !> @param[in,out] A      Spectral function (modified in-place if needed)
          !> @param[in] XKER       Kernel matrix (used only for dimension checks)
          !> @note Modifies global XLAM (Lagrange mults), SIG1 (spectrum norm);
          !> assumes A normalized; sum should be 1 for proper spectral weight.
          !> @pre size(A)=size(XKER,2); XLAM, SIG1 allocated elsewhere
          !> @see SETQ for chi-squared computation
          !-------------------------------------------------------------------

          Subroutine SETA(A,XKER)
            Implicit None

            ! Arguments:
            Real (Kind=Kind(0.d0)), Dimension(:) :: A
            Real (Kind=Kind(0.d0)), Dimension(:,:) :: XKER

            Real (Kind=Kind(0.d0)) :: X
            Integer :: Nw, Nt

            DO NW = 1,NOM
               X = 0.D0
               DO NT = 1,NTAU
                  X  = X + XLAM(NT)*XKER(NT,NW)
               ENDDO
               A(NW) = DEF(NW)*EXP(-X)
               !Write(6,*) 'SetA : ',NW, ' ' ,  X, ' ', A(NW)
            ENDDO
          End Subroutine SETA

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief Compute Hessian (2nd derivative matrix) of MaxEnt objective.
          !> @details Calculates AH(\tau,\tau'), the second derivatives of the
          !> free energy functional: AH = \partial^2 F / \partial \lambda_tau \partial \lambda_tau'
          !> Used by Newton solver (Maximize_Newton) for fast convergence.
          !> Formula: AH(\tau,\tau') = K(\tau,\omega) * diag(A(\omega))^2 * K(\omega,\tau')
          !> Hessian must be inverted for Newton step: \Delta\lambda = -H^{-1} * F.
          !> @param[out] AH      Hessian matrix (Ntau x Ntau, symmetric positive-definite)
          !> @param[in] A        Current spectral function A(omega)
          !> @param[in] XKER     Kernel matrix K(tau,omega)
          !> @param[in] COV      Covariance matrix (for dimension info)
          !> @note Hessian used only by Newton method, not self-consistent iteration;
          !> matrix is expensive O(Ntau^2 * Nom) but enables fast convergence for clean data.
          !> @pre A from SETA; size(XKER)=(Ntau,Nom); size(COV)=(Ntau,Ntau)
          !> @see Maximize_Newton (uses this Hessian for step computation)
          !-------------------------------------------------------------------

          Subroutine SETAH(AH, A,XKER,COV)
            Implicit None
            !Given XLAM, A, and alpha, calculates
            !AH(tau,tau1) = \frac{\partial F_tau} {\partial tau1  }

            ! Arguments
            REAL (Kind=Kind(0.d0)), DIMENSION(:,:) ::  AH, COV, XKER
            REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  A

            Integer NT, NT1, NW
            Real (Kind=Kind(0.d0)) :: X

            IF ( SIZE(AH,1).NE.NTAU .OR. SIZE(AH,2).NE.NTAU) THEN
               WRITE(error_unit,*) 'Error in Setah'
               Call Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
            ENDIF

            DO NT  = 1,NTAU
               DO NT1 = 1,NTAU
                  X = 0.D0
                  DO NW = 1,NOM
                     X = X + XKER(NT,NW)*XKER(NT1,NW)*A(NW)
                  ENDDO
                  AH(NT,NT1) =  COV(NT,NT1)*ALPHA + X
               ENDDO
            ENDDO

          End Subroutine SETAH

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief Compute gradient (first derivatives) of MaxEnt objective.
          !> @details Calculates F(\tau), the first derivatives of the free
          !> energy functional: F(\tau) = \partial E_free / \partial \lambda(\tau).
          !> Used by both Newton and self-consistent solvers. Vector F is the
          !> residual in MaxEnt equations: F = K^T * (A_conv - G_qmc).
          !> Formula: F(\tau) = sum_\omega K(\tau,\omega) * (A_conv(\omega) - G_qmc(\tau))
          !> Newton method: solve F = 0 via \Delta\lambda = -H^{-1} * F.
          !> @param[out] F       Gradient/residual vector (length Ntau)
          !> @param[in] COV      Covariance matrix (for dimension checks)
          !> @param[in] XKER     Kernel matrix K(tau,omega)
          !> @param[in] A        Current spectral function A(omega)
          !> @param[in] XQMC     Target measurements G_qmc(tau)
          !> @note Updated every iteration; convergence achieved when ||F|| < tol.
          !> Global ALPHA from line search affects effective weighting.
          !> @pre A normalized; XKER(tau,omega), size(A)=size(XKER,2)
          !> @see Maximize_Newton (uses F for Newton step)
          !-------------------------------------------------------------------

          Subroutine SETF (F,COV,XKER,A,XQMC)
            Implicit None

            !Given XLAM, A, and alpha, calculates F


            !Arguments
            REAL (Kind=Kind(0.d0)), DIMENSION(:) :: F, A, XQMC
            REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: COV, XKER

            REAL (Kind=Kind(0.d0)) :: X, X1
            Integer :: Nt, Nt1, Nw

            IF (SIZE(F,1).NE.NTAU) THEN
               WRITE(error_unit,*) 'Error in Setf'
               Call Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
            ENDIF
            DO NT = 1,NTAU
               X  = 0.D0
               DO NT1 = 1,NTAU
                  X = X + COV(NT,NT1)*XLAM(NT1)
               ENDDO
               X = ALPHA*X
               X1 = 0.D0
               DO NW = 1,NOM
                  X1 = X1 + XKER(NT,NW)*A(NW)
               ENDDO
               F(NT) = X + XQMC(NT) - X1
            ENDDO
          End Subroutine SETF

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief Compute chi-squared and entropy for MaxEnt convergence.
          !> @details Calculates key objective function components:
          !> XQ = chi^2 = sum_{\tau,\tau'} (G_conv(\tau) - G_qmc(\tau)) * COV^{-1}_{\tau\tau'} * (...)
          !> XENT = entropy S[A] = -\int d\omega A(\omega) ln(A(\omega)/rho_0(\omega))
          !> CHISQ = chi^2 / (Ntau - 1) normalized statistic for convergence.
          !> These quantities track: (1) fit quality via XQ, (2) entropy via XENT,
          !> (3) effective degrees of freedom via CHISQ normalization.
          !> @param[in] A        Spectral function A(omega)
          !> @param[in] XKER     Kernel matrix K(tau,omega)
          !> @param[in] XQMC     QMC measurements G_qmc(tau)
          !> @param[out] XQ      Chi-squared (unnormalized)
          !> @param[out] XENT    Entropy value S[A]
          !> @param[out] CHISQ   Normalized chi-squared (expected ~1 at optimum)
          !> @note Convergence criterion: || grad E || < tol; typical loop exits when CHISQ approaches 1.
          !> XENT computed via numerical integration over frequency grid.
          !> @pre A from SETA; global COV diagonal accessible; size(A)=size(XKER,2)
          !> @see CALCPR_ALP (computes alpha derivative); Maximize_Newton/Self (loop control)
          !-------------------------------------------------------------------

          Subroutine SETQ(A,XKER,XQMC, XQ,XENT,CHISQ)
            Implicit None

            !Arguments
            REAL (Kind=Kind(0.d0)), intent(out)     :: XQ, XENT, CHISQ
            Real (Kind=Kind(0.d0)), Dimension(:)   :: A, XQMC
            Real (Kind=Kind(0.d0)), Dimension(:,:) :: XKER

            !Local
            REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE  ::  VHLP
            Integer :: Nw, Nt, Nt1
            Real (Kind=Kind(0.d0)) :: X

            XENT  = 0.D0
            CHISQ = 0.D0
            ALLOCATE (VHLP(NTAU))

            DO NW = 1,NOM
               X  = A(NW)
               IF (A(NW).LT.ZERO) X  = ZERO
               XENT = XENT + X-DEF(NW) - X*log(X/DEF(NW))
            ENDDO

            DO NT = 1,NTAU
               X  = 0.D0
               DO NW = 1,NOM
                  X = X + XKER(NT,NW)*A(NW)
               ENDDO
               VHLP(NT) = XQMC(NT) - X
            ENDDO

            DO NT1= 1,NTAU
               DO NT = 1,NTAU
                  CHISQ = CHISQ + VHLP(NT)*COVM1(NT,NT1)*VHLP(NT1)
               ENDDO
            ENDDO

            XQ = ALPHA*XENT - CHISQ/2.D0

            DEALLOCATE (VHLP)
          End Subroutine SETQ

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief Compute derivative of free energy with respect to entropy weight.
          !> @details Calculates \partial E_free / \partial \alpha, where alpha is
          !> the entropy weight parameter. Used in line search/bracketing to find
          !> optimal alpha value that minimizes free energy E_free = chi^2 - alpha*S.
          !> Formula: PR_ALP = XQ + 0.5*ln(ALPHA)*Nom - 0.5*ln(det(M))
          !> where M = K^T * diag(A) * K + alpha*I (regularized Gram matrix).
          !> Also computes XTRACE = trace(M * M^{-1}) for Hessian information.
          !> @param[in] XQMC     QMC data G_qmc(tau)
          !> @param[in] COV      Data covariance (used for COV^{-1})
          !> @param[in] A        Current spectrum A(omega)
          !> @param[in] XKER     Kernel matrix K(tau,omega)
          !> @param[in] XQ       Current chi-squared value
          !> @param[in] XENT     Current entropy S[A]
          !> @param[out] PR_ALP  Derivative d(E_free)/d(ALPHA)
          !> @param[out] XTRACE  Trace of product (for second derivative approx)
          !> @note Used in all three MaxEnt variants (T, T0, T_Bryan) via line search;
          !> Sign change in PR_ALP indicates optimal alpha (bracketing root).
          !> Inverse COV stored globally (COVM1); matrix operations O(Nom^3).
          !> @pre COVM1 allocated and inverted; size(A)=size(XKER,2)
          !> @see MaxEnt_T, MaxEnt_T0, MaxEnt_T_Bryan (call for alpha optimization)
          !-------------------------------------------------------------------

          SUBROUTINE CALCPR_ALP(XQMC,  COV, A, XKER,XQ,XENT,PR_ALP,XTRACE)
            Implicit None

            Real (Kind=Kind(0.d0)), Dimension(:)   :: XQMC,  A
            Real (Kind=Kind(0.d0)), Dimension(:,:) :: COV, XKER


            ! Arguments
            REAL (Kind=Kind(0.d0))   :: XQ,XENT, PR_ALP,XTRACE


            ! Local
            REAL (Kind=Kind(0.d0)), DIMENSION(:)                :: DET1(2)
            REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE :: XMAT, XMATM1, XKER1

            Integer :: NFLAG,  NW, NT, NT1, NW1
            REAL (Kind=Kind(0.d0)) :: XLDET

            ALLOCATE (XKER1(NTAU,NOM), XMAT(NOM,NOM), XMATM1(NOM,NOM) )
            XKER1 = 0.D0;    XMAT = 0.D0;   XMATM1 = 0.D0
            NFLAG = 0

            IF (NFLAG.EQ.0) THEN

               !WRITE(6,*) 'Hi1'
               XKER1 = 0.D0
               DO NW  = 1,NOM
                  DO NT  = 1,NTAU
                     DO NT1 = 1,NTAU
                        XKER1(NT,NW) = XKER1(NT,NW)+COVM1(NT,NT1)*XKER(NT1,NW)
                     ENDDO
                     XKER1(NT,NW) = XKER1(NT,NW)*SQRT(A(NW))
                  ENDDO
               ENDDO

               DO NW = 1,NOM
                  DO NW1= 1,NOM
                     XMAT(NW,NW1) = 0.D0
                     DO NT = 1,NTAU
                        XMAT(NW,NW1)=XMAT(NW,NW1)+XKER(NT,NW)*XKER1(NT,NW1)
                     ENDDO
                     XMAT(NW,NW1) =  SQRT(A(NW))*XMAT(NW,NW1)
                  ENDDO
               ENDDO

               DO NW = 1,NOM
                  XMAT(NW,NW) = XMAT(NW,NW) + ALPHA
               ENDDO


               CALL INV(XMAT, XMATM1, DET1)

               DO NW = 1,NOM
                  XMAT(NW,NW) = XMAT(NW,NW) - ALPHA
               ENDDO

               !write(6,*) XQ, ALPHA, NOM, DET1(1), DET1(2)
               XLDET = log(DET1(1)) + DET1(2)*log(10.D0)

               PR_ALP = XQ  + 0.5*log(ALPHA)*DBLE(NOM) - 0.5*XLDET

               XTRACE = 0.D0
               DO NW = 1,NOM
                  DO NW1 = 1,NOM
                     XTRACE = XTRACE + XMAT(NW,NW1)*XMATM1(NW1,NW)
                  ENDDO
               ENDDO


            ENDIF

            DEALLOCATE ( XKER1, XMAT, XMATM1 )

            RETURN
          END SUBROUTINE CALCPR_ALP




          !real (Kind=Kind(0.d0))  function f_fit(k,x)
          !  integer k
          !  real (Kind=Kind(0.d0)) x
          !
          !  if ( k.eq.1) f_fit = 1.d0
          !  if ( k.eq.2) f_fit = x
          !
          !  return
          !end function f_fit

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief MaxEnt method at T=0 with automatic noise filtering.
          !> @details Variant of MaxEnt_T that preprocesses noisy early-time
          !> data by truncating measurements where statistical noise dominates.
          !> Optionally applies exponential shift correction via late-time fit.
          !> Algorithm: (1) Filter data by relative error threshold (Rel_err);
          !> (2) Optionally fit late-time decay and shift via f_fit callback;
          !> (3) Call MaxEnt_T on processed subset; (4) Return spectrum A(omega).
          !> References: von der Linden (1995), Jarrell & Gubernatis (1996).
          !> @param[in] XQMC  QMC measurements at imaginary times (tau)
          !> @param[in] COV   Covariance matrix for XQMC (Ntau x Ntau)
          !> @param[out] A     Spectral function A(omega) (Nom real values)
          !> @param[in] XKER   Kernel matrix XKER(tau,omega) at T=0
          !> @param[in] ALPHA_ST  Starting entropy weight (must be > 0)
          !> @param[out] CHISQ    Final goodness-of-fit: chi^2 = ||G-A_conv||^2
          !> @param[in] Rel_err Relative error threshold for keeping data points
          !> @param[in,out] Shft Optional exponential shift (output if auto-fit applied)
          !> @param[in] xtau Optional imaginary times (required if Shft + f_fit given)
          !> @param[in] f_fit Optional external fit function f(x) for late-time tail
          !> @note XQMC should be centered (zero mean); Rel_err ~ 0.1-0.5 typical.
          !> @pre size(XQMC)=size(A)=Nom; size(COV)=(Ntau,Ntau); ALPHA_ST > 0
          !-------------------------------------------------------------------

          Subroutine MaxEnt_T0 ( XQMC,  COV, A, XKER, ALPHA_ST, CHISQ, Rel_err, Shft, xtau, f_fit)

            Implicit None
            Real (Kind=Kind(0.d0)), Dimension(:)   :: XQMC, A
            Real (Kind=Kind(0.d0)), Dimension(:,:) :: COV, XKER
            Real (Kind=Kind(0.d0)) :: ALPHA_ST, CHISQ,  Rel_err
            Real (Kind=Kind(0.d0)), Optional :: Shft
            Real (Kind=Kind(0.d0)), Dimension(:), Optional :: xtau
            Real (Kind=Kind(0.d0)), external,  Optional :: f_fit

            Real (Kind=Kind(0.d0)), Dimension(:)  , Allocatable   :: XQMC_1
            Real (Kind=Kind(0.d0)), Dimension(:,:), Allocatable   :: COV_1, XKER_1

            ! For the fit if requested.
            Real (Kind=Kind(0.d0)) :: chisq_fit,  Ares(2)
            Real (Kind=Kind(0.d0)), Dimension(:), allocatable :: xdata_fit, fdata_fit,  error_fit
            Integer :: Nd_fit
            !real (Kind=Kind(0.d0)), external  :: f_fit

            Integer nt, nt1, ntau_eff, nw
            Real (Kind=Kind(0.d0)) :: X

            ntau = size(xqmc,1)
            Nom  = Size(A,1)
            ntau_eff = 0
            nt = 0
            do
               nt = nt + 1
               X = sqrt( cov(nt,nt) )/ xqmc(nt)
               if ( X.lt.Rel_err)   then
                  ntau_eff = ntau_eff + 1
               else
                  exit
               endif
               if (nt.eq.ntau)  exit
            enddo
            !write(6,*) 'Ntau_eff: ', Ntau_eff

            !Write(6,*) 'Resizing'
            Allocate ( XQMC_1(Ntau_eff), Cov_1(Ntau_eff,Ntau_eff), Xker_1(Ntau_eff,Nom) )
            do nt = 1,Ntau_eff
               xqmc_1(nt) = xqmc(nt)
            enddo
            do nt = 1,Ntau_eff
               do nt1 = 1,Ntau_eff
                  cov_1(nt,nt1) = cov(nt,nt1)
               enddo
            enddo
            do nt = 1,Ntau_eff
               do nw = 1,Nom
                  XKer_1(nt, nw) = XKer(nt, nw)
               enddo
            enddo
            IF ( PRESENT(Shft) .and. PRESENT(xtau) .and. PRESENT(F_FIT) ) Then
               write(6,*) 'The data will be shifted'
               shft = 0.d0
               Nd_fit = Ntau_eff/2
               Allocate   (xdata_fit(Nd_fit), fdata_fit(Nd_fit),  error_fit(Nd_fit) )
               do  nt = 1,Nd_fit
                  xdata_fit(nt) = xtau(nt +  Nd_fit)
                  fdata_fit(nt) = log(xqmc_1(nt +  Nd_fit))
                  error_fit (nt)  = sqrt( cov_1(nt + Nd_fit,nt + Nd_fit) )/xqmc_1(nt + Nd_fit)
               enddo
               call fit(xdata_fit,fdata_fit,error_fit,ares,chisq_fit,f_fit)
               write(6,*) 'The slope is : ', Ares(2)
               shft = -Ares(2)  - 0.2D0
               Deallocate (xdata_fit, fdata_fit,  error_fit )
               do nt = 1,Ntau_eff
                  xqmc_1(nt) = xqmc_1(nt)*exp(xtau(nt)*shft)
               enddo
               do nt = 1,Ntau_eff
                  do nt1 = 1,Ntau_eff
                     cov_1(nt,nt1) = cov_1(nt,nt1)*exp( (xtau(nt) + xtau(nt1))*shft )
                  enddo
               enddo
            else
               write(6,*) 'The data will not  be shifted'
            endif
            Call MaxEnt_T(XQMC_1,  COV_1, A, XKER_1, ALPHA_ST, CHISQ)
            Deallocate ( Xqmc_1, Cov_1, Xker_1 )


          end Subroutine MaxEnt_T0

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief MaxEnt wrapper with automatic Fermionic kernel construction.
          !> @details Convenience routine that constructs the Fermionic kernel
          !> for arbitrary temperatures T > 0, then invokes MaxEnt_T. Kernel:
          !> K(\tau,\omega) = exp(-\tau*\omega) / (1 + exp(-\beta*\omega))
          !> models QMC propagation with proper thermal/Boltzmann weighting.
          !> Handles all kernel construction internally; requires only imaginary
          !> times, frequencies, and temperature (inverse \beta = 1/k_B T).
          !> @param[in] XTAU    Imaginary time points tau (Ntau values)
          !> @param[in] XQMC    QMC measurements at those times
          !> @param[in] COV     Covariance matrix for XQMC
          !> @param[out] A      Recovered spectral function A(omega)
          !> @param[in] XOM     Frequency grid omega (Nom values)
          !> @param[in] Beta    Inverse temperature (1/T in energy units)
          !> @param[in] ALPHA_ST  Starting entropy weight (must be > 0)
          !> @param[out] CHISQ    Goodness-of-fit chi-squared value
          !> @note Fermionic kernel appropriate for electron/fermion problems;
          !> Bosonic variant in MaxEnt_gr_bose. Temperature must be finite (Beta > 0).
          !> @pre size(XTAU)=size(XQMC)=Ntau; size(XOM)=Nom; Beta > 0; size(COV)=(Ntau,Ntau)
          !-------------------------------------------------------------------

          Subroutine MaxEnt_gr(XTAU, XQMC,  COV,  A, XOM,  Beta, ALPHA_ST, CHISQ )
            ! Sets the Kernel for Green functions.
            Implicit none

            Real (Kind=Kind(0.d0)), Dimension(:)   :: XTAU, XQMC, A, XOM
            Real (Kind=Kind(0.d0)), Dimension(:,:) :: COV

            Real (Kind=Kind(0.d0)) :: ALPHA_ST, CHISQ, BETA


            Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: xker

            Integer :: NT,  NW,  NTAU, NOM


            Nom   = Size(Xom ,1)
            Ntau  = Size(Xtau,1)

            Allocate ( Xker(Ntau,Nom) )
            do nt = 1,ntau
               do nw = 1,Nom
                  XKer(nt,nw) = EXP(-xtau(nt)*xom(nw) ) / ( 1.d0 + EXP( -BETA*xom(nw) ) )
               Enddo
            Enddo

            Call MaxEnt_T(XQMC, COV,  A, XKER, ALPHA_ST, CHISQ )

            Deallocate ( Xker )
          End Subroutine MaxEnt_gr

          !-------------------------------------------------------------------
          !> @author ALF-project
          !> @brief MaxEnt with Bryan's method for entropy weight optimization.
          !> @details Alternative to MaxEnt_T that uses Bryan's search strategy
          !> to find optimal entropy weight ALPHA by bracketing the correct range
          !> [ALPHA_ST, ALPHA_EN]. Uses DIAG eigendecomposition + trust-region
          !> iteration with adaptive step control. Applies entropy constraint and
          !> chi-squared normalization (Ntau-1 degrees of freedom).
          !> References: Bryan (1990), Jarrell & Gubernatis (1996).
          !> @param[in] XQMC  QMC measurements (Ntau real values)
          !> @param[in] COV   Covariance matrix (Ntau x Ntau, symmetric)
          !> @param[out] A     Recovered spectral function (Nom real values)
          !> @param[in] XKER   Kernel matrix XKER(tau,omega) for forward problem
          !> @param[in] ALPHA_ST  Lower bound on entropy weight search
          !> @param[in] ALPHA_EN  Upper bound on entropy weight search
          !> @param[out] CHISQ    Final chi-squared statistic (expected ~ Ntau-1)
          !> @note Bryan method automatically determines best ALPHA; assumes
          !> XQMC centered; tolerates numerical noise in eigendecomposition.
          !> @pre size(XQMC,1)=size(COV,1); size(COV,1)=size(COV,2); ALPHA_EN>ALPHA_ST>0
          !-------------------------------------------------------------------

          Subroutine MaxEnt_T_Bryan( XQMC,  COV, A, XKER, ALPHA_ST, ALPHA_EN, CHISQ )

            Implicit None
            Real (Kind=Kind(0.d0)), Dimension(:)   :: XQMC, A
            Real (Kind=Kind(0.d0)), Dimension(:,:) :: COV, XKER
            Real (Kind=Kind(0.d0)) :: ALPHA_ST, ALPHA_N, ALPHA_EN ! , PI
            Real (Kind=Kind(0.D0)), Intent(Out) :: CHISQ

            Integer       :: NT, NT1, NT2, NW, NCOUNT, NTH
!           Integer       :: NFLAG
            Real (Kind=Kind(0.d0)) :: X, XENT, XQ, PR_ALP, XTRACE, DIFF1, DIFF , Tol_chi_def, XNORM, &
                 &           D_ALPHA, ALPHA_OLD, XNORM_TOT

            Real (Kind=Kind(0.d0)), Dimension(:), allocatable   :: A_ME

            Tol_chi_def = 100000000000000.0D0
            NTAU = SIZE(XQMC,1)
            NOM  = SIZE(A, 1)
            ALLOCATE(A_ME(NOM))
            !WRITE(6,*) 'NTAU, Nom: ', NTAU,NOM
!            PI   = ACOS(-1.d0)
            XMOM1= 1.0D0 !PI
            ZERO =  1.0D-8
            ALLOCATE ( XLAM(NTAU), SIG1(NTAU), COVM1(NTAU,NTAU), UC(NTAU,NTAU), DEF(NOM) )
            XLAM=0.D0;  SIG1=0.D0; UC = 0.D0

            !Open (Unit=77,File='Aom_steps',Status='unknown')
            !Open(Unit=14)
            !do nt = 1, NTAU
            !   Write(14,*) Nt, XQMC(nt), sqrt(Cov(Nt,Nt))
            !enddo
            !Close(14)

            CALL DIAG(COV,UC,SIG1)
            DO NT1 = 1,NTAU
               DO NT2 = 1,NTAU
                  X = 0.D0
                  DO NT = 1,NTAU
                     X = X + UC(NT1,NT)*UC(NT2,NT)/SIG1(NT)
                  ENDDO
                  COVM1(NT1,NT2) = X
               ENDDO
            ENDDO


            Open (Unit=44, File="info_Maxent", Status="unknown", position="append")

            Write(44,*) 'N E W   R U N'
            Write(44,*) '# of data points: ', NTAU
            Write(6,*) 'N E W   R U N'
            ! Set the Default.
            ALPHA     = Alpha_st
            DEF       = XMOM1/dble(NOM)
            XLAM      = 0.d0
            Call Maximize_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)
            IF (CHISQ .GT. Tol_chi_def*NTAU )  THEN
               DO
                  Call Maximize_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)
                  Write(44,*) 'Default: ', Alpha, Chisq
                  Write(6,*) 'Default: ', Alpha, Chisq
                  IF (CHISQ .GT. Tol_chi_def*NTAU .AND.  ALPHA.GT.100 )  THEN
                     ALPHA = ALPHA - ALPHA*0.1
                  ELSE
                     CALL SETA(A,XKER)
                     DO NW = 1,NOM
                        IF (A(NW).LT.ZERO) THEN
                           DEF(NW)= ZERO
                        ELSE
                           DEF(NW) = A(NW)
                        ENDIF
                     ENDDO
                     EXIT
                  ENDIF
               ENDDO
            ELSE
               Write(6,*) 'Flat Default'
            Endif
            !DO NW = 1,NOM
            !   Write(13,*) NW, DEF(NW)
            !ENDDO
            Write(6,*) 'Default Final: ', Alpha, Chisq

            DEF  = XMOM1/dble(NOM)
            Write(6,*) 'Setting the default to a flat default'


            ! Classic MaxEnt.
!            NFLAG  = 0
            NCOUNT = 0
            !ALPHA  = ALPHA_ST
            ALPHA_N = ALPHA_EN
            XLAM   = 0.D0
            NTH    = 0
            A_ME   = 0.d0
            XNORM_TOT = 0.d0
            OPEN (Unit=55,File="Tmp",status="unknown")
            DO
               !WRITE(6,*)  'Starting classic  ', ALPHA
               WRITE(44,*)  '========= 1/Alpha:    ', ALPHA
               !write(6,*) 'Calling maximize'
               CALL MAXIMIZE_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)
               !write(6,*) 'Return: Calling maximize'
               !IF (NFLAG.EQ.0) THEN
                  CALL CALCPR_ALP(XQMC,  COV, A, XKER,XQ,XENT,PR_ALP,XTRACE)
                  IF (NTH.EQ.0)   XNORM = EXP(PR_ALP)
                  NTH = NTH + 1
                  !ALPHA_N = -XTRACE/(2.D0*XENT)
                  WRITE(44,*) 'Max at: 1/Alpha', ALPHA_N
                  WRITE(6,*) 'Max at:', ALPHA_N
                  DIFF1 = ABS(ALPHA_N - ALPHA)
               !ENDIF
               CALL SETA(A,XKER)
               CALL SETQ(A,XKER,XQMC, XQ,XENT,CHISQ)
               WRITE(44,2006) ALPHA, XQ,XENT,CHISQ
               WRITE(6,2006 ) ALPHA, XQ,XENT,CHISQ
               DIFF =  ALPHA_N - ALPHA
               ALPHA_OLD = ALPHA
               IF ( ABS(DIFF) .GT.  0.05*ALPHA ) THEN
                  D_alpha =  0.05 * ALPHA
                  ALPHA  = ALPHA +  0.05 * ALPHA * DIFF/ABS(DIFF)
!                  NFLAG = 1
               ELSE
                  D_alpha =  ABS(ALPHA_N - ALPHA)
                  ALPHA =  ALPHA_N
!                  NFLAG = 0
               ENDIF
               NCOUNT = NCOUNT + 1
               IF (NCOUNT .EQ. 100) THEN
                  WRITE(44,*) 'NOT CONVERGED'
               ENDIF
               WRITE(55,*) ALPHA_OLD, EXP(PR_ALP)/XNORM, D_ALPHA
               XNORM_TOT = XNORM_TOT + D_ALPHA*(EXP(PR_ALP)/XNORM)
               do nw = 1, NOM
                  A_ME(nw) = A_ME(nw) + D_ALPHA*A(nw)*(EXP(PR_ALP)/XNORM)
               enddo
               IF ( ABS(DIFF1)/ABS(ALPHA_N).LT.0.01D0 .OR.  NCOUNT.GT.1000  ) EXIT
           ENDDO
           CLOSE(55)

           A_ME  = A_ME/XNORM_TOT
           A = A_ME
           WRITE(44,*) 'Tot Norm:',  XNORM_TOT
           OPEN(Unit=55,File="Tmp", Status="unknown")
           OPEN(Unit=57,File="Pr_alpha", Status="unknown")
           do
              read(55,*,End=10)  ALPHA_OLD, XNORM, D_ALPHA
              XNORM = XNORM/XNORM_TOT
              write(57,*)  ALPHA_OLD, XNORM, D_ALPHA
           enddo
10         continue
           Close(55)
           Close(57)
           CLOSE(44)


2006       FORMAT('Res: Alpha, XQ,S,CHI: ', F14.7,2x,F14.7,2x,F14.7,2x,F14.7)


            DEALLOCATE ( XLAM, SIG1, COVM1, UC, DEF )
            DEALLOCATE ( A_ME )
            !Close(77)
          End Subroutine MaxEnt_T_Bryan

        end Module MaxEnt_mod
