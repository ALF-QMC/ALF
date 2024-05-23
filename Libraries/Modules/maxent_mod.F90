
!  Copyright (C) 2018-2024 The ALF project
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
!> @author
!> ALF-project
!
!> @brief
!> This module provides  an implementation of  the  classic maximum entropy  algorithm 
!> for  analytical  continuation. It  follows  closely the  work of 
!>   1.   W. von der Linden, Maximum-entropy data analysis, Applied Physics A 60 (1995), no. 2, 155–165. 
!>   and 
!>   2.  Mark Jarrell and J.E. Gubernatis, Bayesian inference and the analytic continuation 
!>   of imaginary-time quantum monte carlo data, Physics Reports 269 (1996), no. 3, 133–195.

!
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



          Subroutine Maximize_Newton( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)

            ! Sloves F(tau) = 0 with Newton.


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

          Subroutine Maximize_Self( XQMC,  COV, A, XKER, XQ,XENT,CHISQ)

            ! Sloves F(tau) = 0 with self-consistency.
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

          Subroutine SETAH(AH, A,XKER,COV)
            Implicit None
            !Given XLAM,  A, and alpha,  calcluates
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

          Subroutine SETF (F,COV,XKER,A,XQMC)
            Implicit None

            !Given XLAM,  A, and alpha,  calcluates F


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

          !>  Shft A Shift

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
