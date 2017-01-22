!  Copyright (C) 2016 The ALF project
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

      SUBROUTINE CGR(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>    Computes  GRUP = (1 + UR*DR*VR*VL*DL*UL)^-1
!>    NVAR = 1 Big scales are in DL
!>    NVAR = 2 Big scales are in DR
!
!--------------------------------------------------------------------

        Use UDV_Wrap_mod

        Implicit None
        Integer :: IZAMAX, izmax1
        REAL(Kind=Kind(0.D0)) :: DZSUM1, DZNRM2

	!Arguments.
        COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(IN)   ::  URUP, VRUP, ULUP, VLUP
        COMPLEX(Kind=Kind(0.d0)), Dimension(:),   Intent(IN)   ::  DLUP, DRUP
        COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX(Kind=Kind(0.d0)), Intent(INOUT) :: PHASE
        INTEGER         :: NVAR
 
        !Local
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Allocatable ::  TPUP, RHS
        COMPLEX (Kind=Kind(0.d0)), Dimension(:) , Allocatable ::  DUP
        INTEGER, Dimension(:), Allocatable :: IPVT, VISITED
        COMPLEX (Kind=Kind(0.d0)) ::  alpha, beta, Z
        Integer :: I, J, N_size, NCON, info, LWORK, next, L
        Real (Kind=Kind(0.D0)) :: X, Xmax, sv
        
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK, RWORK
        LOGICAL :: FORWRD
        
        N_size = SIZE(DLUP,1)
        NCON = 0
        alpha = 1.D0
        beta = 0.D0
        Allocate(TPUP(N_size,N_size), VISITED(N_size), RHS(N_size, N_size),&
             &DUP(N_size), IPVT(N_size), TAU(N_size), RWORK(2*N_size))
        !Write(6,*) 'In CGR', N_size
        CALL MMULT(TPUP,VRUP,VLUP)
        DO J = 1,N_size
            TPUP(:,J) = DRUP(:)*TPUP(:,J)*DLUP(J)
        ENDDO
        ! can be inserted again once we are sure that we may assume that UR and UL stem from householder reflectors
!        CALL ZGEMM('C', 'C', N_size, N_size, N_size, alpha, URUP, N_size, ULUP, N_size, alpha, TPUP, N_size)
        CALL ZGEMM('C', 'C', N_size, N_size, N_size, alpha, URUP, N_size, ULUP, N_size, beta, RHS, N_size)
        TPUP = TPUP + RHS
        ! calculate determinant of UR*UL
        PHASE = CONJG(DET_C(RHS, N_size))
        PHASE = PHASE/ABS(PHASE)
        IPVT = 0
        IF (NVAR .NE. 1) THEN
            TPUP = CT(TPUP)
        ENDIF
        ! Query and allocate optimal amount of work space
        call ZGEQP3(N_size, N_size, TPUP, N_size, IPVT, TAU, DUP(1), -1, RWORK, INFO)
        LWORK = INT(DBLE(DUP(1)))
        ALLOCATE(WORK(LWORK))
        ! QR decomposition of TPUP with full column pivoting, AP = QR
        call ZGEQP3(N_size, N_size, TPUP, N_size, IPVT, TAU, WORK, LWORK, RWORK, INFO)
!         CALL PM(TPUP, N_size)
        ! Another attempt at calculating the sign of P
        VISITED = 0
        do i = 1, N_size
            if (VISITED(i) .eq. 0) then
                next = i
                L = 0
                do while (VISITED(next) .eq. 0)
                 L = L + 1
                 VISITED(next) = 1
                 next = IPVT(next)
                enddo
                if(MOD(L, 2) .eq. 0) then
                    PHASE = -PHASE
                endif
            endif
        enddo
        !calculate the determinant of the unitary matrix Q
    DO i = 1, N_size
        Z = TAU(i)
        X = ABS(Z)
        if (Z .ne. CMPLX(0.D0, 0.D0, Kind=Kind(0.D0))) then
         Z = 1.D0 - 2.D0 * (Z/X) * (DBLE(Z)/X)
        IF (NVAR .EQ. 1) THEN
            PHASE = PHASE * Z/ABS(Z)
        ELSE
            PHASE = PHASE * CONJG(Z)/ABS(Z)
        ENDIF
        endif
    ENDDO

        ! consider the triangular right matrix R
        DO i = 1, N_size
        IF(NVAR .EQ. 1) THEN
            PHASE = PHASE * TPUP(i,i)/Abs(TPUP(i,i))
        ELSE
            ! the diagonal should be real, but let's be safe....
            PHASE = PHASE * CONJG(TPUP(i,i))/Abs(TPUP(i,i))
        ENDIF
        enddo
        ! separate off DUP
        do i = 1, N_size
        ! plain diagonal entry
            X = ABS(TPUP(i, i))
!             ! a inf-norm
!             X = TPUP(i, i+izamax(Ndim+1-i, TPUP(i, i), Ndim)-1)
!             ! another inf-norm
!             X = TPUP(i, i-1+izmax1(Ndim+1-i, TPUP(i, i), Ndim))
!             ! 1-norm
!            X = DZSUM1(N_size+1-i, TPUP(i, i), N_size)
            ! 2-norm
!            X = DZNRM2(N_size+1-i, TPUP(i, i), N_size)
            DUP(i) = X
            do j = i, N_size
                TPUP(i, j) = TPUP(i, j) / X
            enddo
        enddo
            
        IF(NVAR .EQ. 1) then
            ! This is supposed to solve the system 
            ! URUP U D V P^dagger ULUP G = 1
            ! initialize the rhs with CT(URUP)
            RHS = CT(URUP)
            ! RHS = U^dagger * RHS
            CALL ZUNMQR('L', 'C', N_size, N_size, N_size, TPUP, N_size, TAU, RHS, N_size, WORK, LWORK, INFO)
            DEALLOCATE(TAU, WORK, RWORK)   
            !apply inverse of D to RHS
        DO J = 1, N_size
           sv = DBLE(DUP(J))
           X = ABS(sv)
           if (J == 1)  Xmax = X
           if ( X  < Xmax ) Xmax = X
           sv = 1.D0/sv
           DO I = 1, N_size
              RHS(I,J) = RHS(I, J) / DUP(I)
           ENDDO
        ENDDO
        ! We solve the equation
        !  A * G = RHS for G with A = R * P^dagger * ULUP
        ! first we solve R *y = RHS. The solution is afterwards in RHS
        CALL ZTRSM('L', 'U', 'N', 'N', N_size, N_size, alpha, TPUP, N_size, RHS, N_size)
        ! apply permutation matrix
            FORWRD = .false.
            CALL ZLAPMR(FORWRD, N_size, N_size, RHS, N_size, IPVT)
        ! perform multiplication with ULUP and store in GRUP
        CALL ZGEMM('C', 'N', N_size, N_size, N_size, alpha, ULUP, N_size, RHS, N_size, beta, GRUP, N_size)
        ELSE
            ! This solves the system G * URUP * P * R^dagger * D * U^dagger * ULUP = 1
            
           ! RHS = ULUP * UUP
           RHS = CT(ULUP)
           CALL ZUNMQR('R', 'N', N_size, N_size, N_size, TPUP, N_size, TAU, RHS, N_size, WORK, LWORK, INFO)
        DEALLOCATE(TAU, WORK, RWORK)

        ! apply D^-1 to RHS
        DO J = 1, N_size
           sv = DBLE(DUP(J))
           X = ABS(sv)
           if (J == 1)  Xmax = X
           if ( X  < Xmax ) Xmax = X
           sv = 1.D0/sv
           DO I = 1, N_size
             RHS(I, J) = RHS(I, J) * sv
           ENDDO
        ENDDO
        
        ! We solve the equation
        ! G * A = RHS for G with A = URUP * P * R^dagger
        ! first we solve y * R^dagger = RHS
        CALL ZTRSM('R', 'U', 'C', 'N', N_size, N_size, alpha, TPUP, N_size, RHS, N_size)
        ! apply inverse permutation matrix
            FORWRD = .false.
            CALL ZLAPMT(FORWRD, N_size, N_size, RHS, N_size, IPVT)
        ! perform multiplication with URUP
        CALL ZGEMM('N', 'C', N_size, N_size, N_size, alpha, RHS, N_size, URUP, N_size, beta, GRUP, N_size)
        ENDIF
        Deallocate(TPUP, DUP, IPVT, VISITED, RHS)
      END SUBROUTINE CGR
