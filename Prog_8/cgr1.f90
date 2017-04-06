!  Copyright (C) 2016, 2017 The ALF project
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

      SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>    Computes  GRUP = (1 + UR*DR*VR*VL*DL*UL)^-1
!>    and      PHASE = det(1 + UR*DR*VR*VL*DL*UL) / abs(det(1 + UR*DR*VR*VL*DL*UL)) 
!>    NVAR = 1 Big scales are in DL
!>    NVAR = 2 Big scales are in DR
!> Implementation note: we calculate the Phase as:
!> NVAR = 1 : Phase = det(URUP * ULUP)/ |det(URUP * ULUP)| * det(P) * det(R) *det(Q)/ |det(R) det(Q)| 
!> NVAR = 2 : Phase = det(URUP * ULUP)/ |det(URUP * ULUP)| * det(P) * det^*(R) *det^*(Q)/ |det(R) det(Q)| 
!
!--------------------------------------------------------------------

        Use UDV_State_mod

#if defined(STAB2) || defined(STAB1)
        Use UDV_Wrap_mod

        Implicit None


	!Arguments.
	CLASS(UDV_State), INTENT(IN) :: udvl, udvr
        COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX(Kind=Kind(0.d0)) :: PHASE
        INTEGER         :: NVAR
 
        !Local
        TYPE(UDV_State) :: udvlocal
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Allocatable :: TPUP, TPUP1, TPUPM1
        INTEGER, Dimension(:), Allocatable :: IPVT
        COMPLEX (Kind=Kind(0.d0)) ::  ZDUP1, ZDDO1, ZDUP2, ZDDO2, Z1, ZUP, ZDO, alpha, beta
        Integer :: I,J, N_size, NCON, info
        Real (Kind=Kind(0.D0)) :: X, Xmax, sv
        
        N_size = udvl%Ndim
        NCON = 0
        alpha = 1.D0
        beta = 0.D0
        Allocate(TPUP(N_size,N_size), TPUP1(N_size,N_size), TPUPM1(N_size,N_size), IPVT(N_size) )
        CALL udvlocal%alloc(N_size)
        !Write(6,*) 'In CGR', N_size
        CALL MMULT(udvlocal%V, udvr%V, udvl%V)
        DO J = 1,N_size
            TPUP(:,J) = udvr%D(:)*udvlocal%V(:,J)*udvl%D(J)
        ENDDO
        CALL ZGEMM('C', 'C', N_size, N_size, N_size, alpha, udvr%U(1,1), N_size, udvl%U(1,1), N_size, alpha, TPUP, N_size)
        !>  Syntax 
        !>  ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        !>  C := alpha*op( A )*op( B ) + beta*C
        !>  TPUP =  (URUP)^(dagger) ULUP^(dagger) + TPUP
        IF (NVAR.EQ.1) THEN
           !WRITE(6,*) 'UDV of U + DR * V * DL'
           CALL UDV_WRAP_Pivot(TPUP,udvlocal%U,udvlocal%D,udvlocal%V,NCON,N_size,N_Size)
           !CALL UDV(TPUP,UUP,DUP,VUP,NCON)
           CALL MMULT(TPUP,udvlocal%V, udvl%U)
           !Do I = 1,N_size
           !   Write(6,*) DLUP(I)
           !enddo
           CALL INV  (TPUP,TPUPM1,ZDUP1)
           CALL MMULT(TPUP1, udvr%U, udvlocal%U)
           CALL ZGETRF(N_size, N_size, TPUP1, N_size, IPVT, info)
           !>  TPUP1 = P * L * U   LU-decomposition
           Z1 = ZDUP1
           Do i = 1, N_size
              IF (IPVT(i) .ne. i) THEN
                 Z1 = -Z1
              endif
              Z1 = Z1 * TPUP1(I, I) 
              
           enddo
        ELSE
           !WRITE(6,*) 'UDV of (U + DR * V * DL)^{*}'
           TPUP1 = CT(TPUP)
           CALL UDV_WRAP_Pivot(TPUP1, udvlocal%U, udvlocal%D, udvlocal%V, NCON,N_size,N_size)
           !CALL UDV(TPUP1,UUP,DUP,VUP,NCON)
           CALL ZGEMM('C', 'N', N_size, N_size, N_size, alpha, udvl%U(1,1), N_size, udvlocal%U, N_size, beta, TPUPM1, N_size)
           CALL ZGEMM('N', 'C', N_size, N_size, N_size, alpha, udvr%U(1,1), N_size, udvlocal%V, N_size, beta, TPUP1, N_size)
           CALL ZGETRF(N_size, N_size, TPUP1, N_size, IPVT, info)
           !>  TPUP1 = P * L * U   LU-decomposition
           ZDUP2 = 1.D0
           do i = 1, N_size
              ZDUP2 = ZDUP2 * TPUP1(I,I)
              IF (IPVT(i) .ne. i) THEN
                 ZDUP2 = -ZDUP2
              endif
           enddo
           TPUP = TPUPM1
           ZDUP1 = DET_C(TPUP, N_size)! Det destroys its argument
           Z1 = ZDUP2/ZDUP1
        ENDIF
        DO J = 1, N_size
           sv = DBLE(udvlocal%D(J))
           X = ABS(sv)
           if (J == 1)  Xmax = X
           if ( X  < Xmax ) Xmax = X
           sv = 1.D0/sv
           DO I = 1, N_size
              udvlocal%U(J, I) = TPUPM1(I, J) * sv
           ENDDO
        ENDDO
        call ZGETRS('T', N_size, N_size, TPUP1, N_size, IPVT, udvlocal%U, N_size, info)
        !> Syntax
        !> ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        !> Op(A) * X = B 
        !> On output  B = X
        GRUP = TRANSPOSE(udvlocal%U)
        PHASE = Z1/ABS(Z1)
        CALL udvlocal%dealloc
        Deallocate(TPUP,TPUP1,TPUPM1, IPVT )

#else

        USE MyMats
        USE QDRP_mod
        Implicit None
	!Arguments.
!         COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(IN)   ::  URUP, VRUP, ULUP, VLUP
!         COMPLEX(Kind=Kind(0.d0)), Dimension(:),   Intent(IN)   ::  DLUP, DRUP
        CLASS(UDV_State), INTENT(IN) :: udvl, udvr
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
        
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK
        LOGICAL :: FORWRD
        
        N_size = udvl%ndim
        NCON = 0
        alpha = 1.D0
        beta = 0.D0
        Allocate(TPUP(N_size,N_size), RHS(N_size, N_size), IPVT(N_size), TAU(N_size), DUP(N_size))
        !Write(6,*) 'In CGR', N_size
        CALL MMULT(TPUP, udvr%V, udvl%V)
        DO J = 1,N_size
            TPUP(:,J) = udvr%D(:) *TPUP(:,J)*udvl%D(J)
        ENDDO
        ! can be inserted again once we are sure that we may assume that UR and UL stem from householder reflectors
!        CALL ZGEMM('C', 'C', N_size, N_size, N_size, alpha, URUP, N_size, ULUP, N_size, alpha, TPUP, N_size)
        CALL ZGEMM('C', 'C', N_size, N_size, N_size, alpha, udvr%U, N_size, udvl%U, N_size, beta, RHS(1, 1), N_size)
        TPUP = TPUP + RHS
        ! calculate determinant of UR*UL
        PHASE = CONJG(DET_C(RHS, N_size))
        PHASE = PHASE/ABS(PHASE)
        IPVT = 0
        IF (NVAR .NE. 1) THEN
            TPUP = CONJG(TRANSPOSE(TPUP))
        ENDIF
        call QDRP_decompose(N_size, TPUP, DUP, IPVT, TAU, WORK, LWORK)
        ALLOCATE(VISITED(N_size))
        ! Calculate the sign of the permutation from the pivoting. Somehow the format used by the QR decomposition of lapack
        ! is different from that of the LU decomposition of lapack
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
        !calculate the determinant of the unitary matrix Q and the upper triangular matrix R
        DO i = 1, N_size
            Z = TAU(i)
            IF(NVAR .EQ. 1) THEN
                PHASE = PHASE * TPUP(i,i)/Abs(TPUP(i,i))
            ELSE
                ! the diagonal should be real, but let's be safe....
                PHASE = PHASE * CONJG(TPUP(i,i))/Abs(TPUP(i,i))
                Z = CONJG(Z) ! conjugate the elementary reflector
            ENDIF
            if (Z .ne. CMPLX(0.D0, 0.D0, Kind=Kind(0.D0))) then
            ! here we calculate the determinant of a single householder reflector: det(1 - tau * v v* ) = 1 - tau * v^* v
            ! In lapack the scalar tau and the vector v are scaled such that |tau|^2 |v|^2 = 2 Re(tau)
            ! The complete determinant det(Q) is the product of all reflectors. See http://www.netlib.org/lapack/lug/node128.html
                X = ABS(Z)
                Z = 1.D0 - 2.D0 * (Z/X) * (DBLE(Z)/X)
                PHASE = PHASE * Z/ABS(Z)
            endif
        enddo
        IF(NVAR .EQ. 1) then
            ! This is supposed to solve the system 
            ! URUP U D V P^dagger ULUP G = 1
            ! initialize the rhs with CT(URUP)
            RHS = CT(udvr%U)
            ! RHS = U^dagger * RHS
            CALL ZUNMQR('L', 'C', N_size, N_size, N_size, TPUP(1, 1), N_size, TAU(1), RHS(1,1), N_size, WORK(1), LWORK, INFO)
            DEALLOCATE(TAU, WORK)
            !apply inverse of D to RHS from the left
            DO J = 1, N_size
                sv = DBLE(DUP(J))
                X = ABS(sv)
                if (J == 1)  Xmax = X
                if ( X  < Xmax ) Xmax = X
                DO I = 1, N_size
                    RHS(I,J) = RHS(I, J) / DUP(I)
                ENDDO
            ENDDO
            ! We solve the equation
            !  A * G = RHS for G with A = R * P^dagger * ULUP
            ! first we solve R *y = RHS. The solution is afterwards in RHS
            CALL ZTRSM('L', 'U', 'N', 'N', N_size, N_size, alpha, TPUP(1,1), N_size, RHS(1,1), N_size)
            ! apply permutation matrix
            FORWRD = .false.
            CALL ZLAPMR(FORWRD, N_size, N_size, RHS(1,1), N_size, IPVT(1))
            ! perform multiplication with ULUP and store in GRUP
            CALL ZGEMM('C', 'N', N_size, N_size, N_size, alpha, udvl%U(1, 1), N_size, RHS(1,1), N_size, beta, GRUP(1,1), N_size)
        ELSE
            ! This solves the system G * URUP * P * R^dagger * D * U^dagger * ULUP = 1
            
            ! RHS = ULUP * UUP
            RHS = CT(udvl%U)
            CALL ZUNMQR('R', 'N', N_size, N_size, N_size, TPUP(1, 1), N_size, TAU(1), RHS(1, 1), N_size, WORK(1), LWORK, INFO)
            DEALLOCATE(TAU, WORK)
            ! apply D^-1 to RHS from the right
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
            CALL ZTRSM('R', 'U', 'C', 'N', N_size, N_size, alpha, TPUP(1, 1), N_size, RHS(1, 1), N_size)
            ! apply inverse permutation matrix
            FORWRD = .false.
            CALL ZLAPMT(FORWRD, N_size, N_size, RHS(1, 1), N_size, IPVT(1))
            ! perform multiplication with URUP
            CALL ZGEMM('N', 'C', N_size, N_size, N_size, alpha, RHS(1, 1), N_size, udvr%U(1,1), N_size, beta, GRUP(1, 1), N_size)
        ENDIF
        Deallocate(TPUP, DUP, IPVT, VISITED, RHS)
#endif
        
      END SUBROUTINE CGR
