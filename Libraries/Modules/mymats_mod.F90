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



!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Wrappers for linear algebra.
!>
!
!--------------------------------------------------------------------


    MODULE MyMats
      use iso_fortran_env, only: output_unit, error_unit
      use runtime_error_mod

       INTERFACE MMULT
          !C = A*B MMULT(C, A, B)
          MODULE PROCEDURE MMULT_R, MMULT_C
       END INTERFACE
       INTERFACE INITD
          MODULE PROCEDURE INITD_R, INITD_C
       END INTERFACE
       INTERFACE COMPARE
          MODULE PROCEDURE COMPARE_R, COMPARE_C
       END INTERFACE
       INTERFACE DET
          MODULE PROCEDURE DET_C
       END INTERFACE DET
       INTERFACE CT
          MODULE PROCEDURE CT
       END INTERFACE CT
       INTERFACE INV
          MODULE PROCEDURE INV_R0, INV_R_Variable, INV_R_VARIABLE_1, INV_R1, INV_R2, INV_C, INV_C1, &
               &        INV_C_Variable
       END INTERFACE
       INTERFACE UDV
          MODULE PROCEDURE UDV1_R, UDV_C
       END INTERFACE
       INTERFACE QR
          MODULE PROCEDURE QR_C
       END INTERFACE QR
       INTERFACE SVD
          MODULE PROCEDURE SVD_C
       END INTERFACE SVD
       INTERFACE DIAG
          MODULE PROCEDURE DIAG_R, DIAG_I
       END INTERFACE
       INTERFACE DIAG_GEN
          MODULE PROCEDURE DIAG_GEN
       END INTERFACE DIAG_GEN
       INTERFACE SECONDS
          MODULE PROCEDURE SECONDS
       END INTERFACE
     CONTAINS

!--------------------------------------------------------------------
       SUBROUTINE DIAG_GEN(Z_MAT,U,W,LR,ICON)
         IMPLICIT NONE
         COMPLEX   (Kind=Kind(0.d0)), INTENT(IN), DIMENSION(:,:) :: Z_MAT
         CHARACTER (LEN=1),  INTENT(IN)  :: LR
         COMPLEX   (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:,:) :: U
         COMPLEX   (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:) :: W
         INTEGER :: ICON

         !!!! Uses Lapack !!!
         ! LR = L  then         U*A   = W*U   Left  eigenvectors
         ! LR = R  then         A*U   = W*U   Right eigenvectors


         !  Local space
         INTEGER :: N, LDA, LDVL, LDVR, INFO, LWORK, I, J, M
         CHARACTER (LEN=1) ::   JOBVL, JOBVR
         COMPLEX (Kind=Kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A, VL, VR
         REAL (Kind=Kind(0.d0))   , ALLOCATABLE, DIMENSION(:) :: RWORK
         COMPLEX (Kind=Kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: WORK

         REAL    (Kind=Kind(0.d0)) :: XMAX, X
         COMPLEX (Kind=Kind(0.d0)) :: Z

         N = SIZE(Z_MAT,1)
         ALLOCATE(A(N,N))
         A = Z_MAT
         LDA = N

         JOBVR  = "N"
         JOBVL  = "N"
         LDVL = 1
         LDVR = 1
         IF (LR == "L" .or. LR == 'l') THEN
            JOBVL ="V"
            LDVL  = N
         ELSEIF (LR == "R" .or. LR == 'r') THEN
            JOBVR ="V"
            LDVR = N
         ELSE
            WRITE(error_unit,*) 'Error in DIAG_GEN'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         ENDIF
         ALLOCATE(VL(LDVL,N),  VR(LDVR,N) )
         LWORK = 2*N
         ALLOCATE (WORK(LWORK), RWORK(LWORK) )

         CALL ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
              &      WORK, LWORK, RWORK, INFO )

         IF (LR == "R" .or. LR == 'r')  THEN
            DO I = 1,N
               DO J = 1,N
                  U(I,J) = VR(I,J)
               ENDDO
            ENDDO
         ELSE
            DO I = 1,N
               DO J = 1,N
                  U(I,J) = CONJG(VL(J,I))
               ENDDO
            ENDDO
         ENDIF

         IF (ICON == 1 ) THEN
            !Test
            XMAX = 0.d0
            DO I = 1,N
               DO J = 1,N
                  IF (LR == "R" .or. LR == 'r')  THEN
                     Z = - W(I)*U(I,J)
                     DO M = 1,N
                        Z = Z + Z_MAT(I,M)*U(M,J)
                     ENDDO
                  ELSE
                     Z = -W(I)*U(I,J)
                     DO M = 1,N
                        Z = Z + U(I,M)*Z_MAT(M,J)
                     ENDDO
                  ENDIF
                  X = ABS(Z)
                  IF ( X > XMAX ) XMAX = X
               ENDDO
            ENDDO
            WRITE(6,*) 'Testing Diag_GEN :', XMAX
            !End Test
         ENDIF

         DEALLOCATE(VL, VR)
         DEALLOCATE(WORK, RWORK)
         DEALLOCATE(A)


       END SUBROUTINE DIAG_GEN
!--------------------------------------------------------------------
       SUBROUTINE MMULT_R(C, A, B)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: A,B,C
         REAL (Kind=Kind(0.d0)) :: ALP, BET
         INTEGER N, M, P, LDA, LDB, LDC
         N = SIZE(A,1) ! Rows in A
         M = SIZE(A,2) ! Columns in A
         P = SIZE(B,2) ! Columns in B
         LDA = N; LDB = SIZE(B,1); LDC = SIZE(C,1)

         ALP = 1.D0
         BET = 0.D0

         CALL DGEMM('n','n',N,P,M,ALP,A,LDA,B,LDB,BET,C,LDC)


! WRITE(6,*) 'In real', N,M,P
! DO I = 1,N
! DO J = 1,P
! X = 0.D0
! DO K = 1,M
! X = X + A(I,K)*B(K,J)
! ENDDO
! C(I,J) = X
! ENDDO
! ENDDO
       END SUBROUTINE MMULT_R

       SUBROUTINE MMULT_C(C, A, B)
         IMPLICIT NONE
         COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:) :: A,B,C
         COMPLEX (Kind=Kind(0.d0)) :: ALP, BET
         INTEGER N, M, P, LDA, LDB, LDC

         N = SIZE(A,1)
         M = SIZE(A,2)
         P = SIZE(B,2)
         LDA = N; LDB = SIZE(B,1); LDC = SIZE(C,1)

         ALP = CMPLX(1.D0,0.D0,Kind=Kind(0d0))
         BET = CMPLX(0.D0,0.D0,Kind=Kind(0d0))

         CALL ZGEMM('n','n',N,P,M,ALP,A,LDA,B,LDB,BET,C,LDC)


         ! WRITE(6,*) 'In complex', N,M,P
         ! DO I = 1,N
         ! DO J = 1,P
         ! X = CMPLX(0.D0,0.D0)
         ! DO K = 1,M
         ! X = X + A(I,K)*B(K,J)
         ! ENDDO
         ! C(I,J) = X
         ! ENDDO
         ! ENDDO

       END SUBROUTINE MMULT_C

!*********
       SUBROUTINE INITD_R(A,X)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: A
         REAL (Kind=Kind(0.d0)) X
         INTEGER I,J, N, M

         N = SIZE(A,1)
         M = SIZE(A,2)

         ! WRITE(6,*) 'In Init1 real', N,M
         DO I = 1,N
         DO J = 1,M
            A(I,J) = 0.D0
         ENDDO
         ENDDO
         DO I = 1,N
            A(I,I) = X
         ENDDO
       END SUBROUTINE INITD_R

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This functions sets the matrix to a diagonal matrix with identical
!> entries on the diagonal.
!
!> @param[inout] A a 2D array constituting the input matrix.
!> @param[in] Xthe scalar that we set the diagonal to.
!--------------------------------------------------------------------
       SUBROUTINE INITD_C(A, X)
         IMPLICIT NONE
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: A
         COMPLEX (KIND=KIND(0.D0)), INTENT(IN) :: X
         INTEGER I, N

         A = (0.D0, 0.D0)
         N = SIZE(A,1)
         DO I = 1,N
            A(I, I) = X
         ENDDO
       END SUBROUTINE INITD_C

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function calculates the LU decomposition and the determinant
!> of the input matrix.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] AINV a 2D array containing the inverse of the matrix A.
!> @param[out] DET the determinant of the input matrix.
!--------------------------------------------------------------------
       SUBROUTINE INV_R0(A, AINV, DET)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), INTENT(OUT) :: DET
         INTEGER I

! Working space.
         REAL (KIND=KIND(0.D0)) :: SGN
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, LDA

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         AINV = A
         CALL DGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         DET = 1.D0
         SGN = 1.D0
         DO i = 1, LDA
         DET = DET * AINV(i,i)
         IF (IPVT(i) .ne. i) THEN
            SGN = -SGN
         ENDIF
         enddo
         DET = SGN * DET
         CALL DGETRI(LDA, AINV, LDA, IPVT, WORK, LDA, INFO)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R0


!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function calculates the LU decomposition and the determinant
!> in a subpart of the input matrix.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] AINV a 2D array containing the inverse of the subpart.
!> @param[out] DET the determinant of the input matrix.
!> @param[in] Ndim The size of the subpart.
!--------------------------------------------------------------------
       SUBROUTINE INV_R_Variable(A, AINV, DET, Ndim)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), INTENT(OUT) :: DET
         INTEGER, INTENT(IN) :: Ndim

! Working space.
         REAL (KIND=KIND(0.D0)) :: SGN
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(LDA) )

         AINV = A
         CALL DGETRF(Ndim, Ndim, AINV, LDA, IPVT, INFO)
         DET = 1.D0
         SGN = 1.D0
         DO i = 1, Ndim
         DET = DET * AINV(i,i)
         IF (IPVT(i) .ne. i) THEN
            SGN = -SGN
         ENDIF
         ENDDO
         DET = SGN * DET
         CALL DGETRI(Ndim, AINV, LDA, IPVT, WORK, LDA, INFO)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R_VARIABLE

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function calculates the LU decomposition and the determinant
!> in a subpart of the input matrix.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] AINV a 2D array containing the inverse of the subpart
!> @param[out] DET The determinant in a kind of Mantissa-Exponent like
!>                 representation. The full determinant is d1*10^d2
!>                 (1.0 <= |d1| < 10.0) OR (d1 == 0.0)
!> @param[in] Ndim The size of the subpart.
!
!> TODO: In a test the best accuracy could be obtained using the log10
!! below. Another possibility would be using the EXPONENT and FRACTION
!! intrinsics. This turned out to be not so good.
!! Currently the case of a singular matrix is not handled, since it was
!! catched in the old linpack version.
!--------------------------------------------------------------------
       SUBROUTINE INV_R_Variable_1(A,AINV,DET,Ndim)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), DIMENSION(2), INTENT(OUT) :: DET
         INTEGER, INTENT(IN) :: Ndim

! Working space.
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(LDA) )

         AINV = A
         CALL DGETRF(Ndim, Ndim, AINV, LDA, IPVT, INFO)
         DET(1) = 1.D0
         DET(2) = 0.D0
!         SGN = 1.0
         DO i = 1, Ndim
         IF (AINV(i, i) < 0.D0) THEN
         DET(1) = -DET(1)
         ENDIF
         DET(2) = DET(2) + LOG10(ABS(AINV(i,i)))
         IF (IPVT(i) .ne. i) THEN
            DET(1) = -DET(1)
         ENDIF
         ENDDO
         CALL DGETRI(Ndim, AINV, LDA, IPVT, WORK, LDA, INFO)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R_VARIABLE_1

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function calculates the LU decomposition and the determinant
!> of the input matrix.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] AINV a 2D array containing the inverse of the matrix A.
!> @param[out] DET The determinant in a kind of Mantissa-Exponent like
!>                 representation. The full determinant is d1*10^d2
!>                 (1.0 <= |d1| < 10.0) OR (d1 == 0.0)
!
!> @todo the same restrictions as in INV_R_VARIABLE_1 apply.
!--------------------------------------------------------------------
       SUBROUTINE INV_R1(A,AINV,DET)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), DIMENSION(2), INTENT(OUT) :: DET

! Working space.
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         AINV = A
         CALL DGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         DET(1) = 1.D0
         DET(2) = 0.D0
!         SGN = 1.0
         DO i = 1, LDA
         IF (AINV(i, i) < 0.D0) THEN
         DET(1) = -DET(1)
         ENDIF
         DET(2) = DET(2) + LOG10(ABS(AINV(i,i)))
         IF (IPVT(i) .ne. i) THEN
            DET(1) = -DET(1)
         ENDIF
         ENDDO
         CALL DGETRI(LDA, AINV, LDA, IPVT, WORK, LDA, INFO)
         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R1

!*************
       SUBROUTINE INV_R2(A,AINV)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: A,AINV

         INTEGER I, J

! Uses Lapack routines.

! Working space.
         REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
         INTEGER INFO, LDA, LWORK

         LDA = SIZE(A,1)

         !Write(6,*) 'Inv_r2:', LDA
         ALLOCATE ( IPIV(LDA) )
         LWORK = LDA
         ALLOCATE ( WORK(LWORK) )
         WORK = 0.0
         IPIV = 0
         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO
         INFO = 0


         CALL DGETRF( LDA, LDA, AINV, LDA, IPIV, INFO )
         CALL DGETRI(LDA, AINV, LDA, IPIV, WORK, LWORK, INFO)





         DEALLOCATE (IPIV)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R2

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function calculates the LU decomposition and the determinant
!> of a complex input matrix.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] AINV a 2D array containing the inverse of the matrix A.
!> @param[out] DET the determinant of the input matrix.
!--------------------------------------------------------------------

       SUBROUTINE INV_C(A,AINV,DET)
         IMPLICIT NONE
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         COMPLEX (KIND=KIND(0.D0)), INTENT(OUT) :: DET
         INTEGER I

! Working space.
         REAL (KIND=KIND(0.D0)) :: SGN
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, LDA

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         AINV = A
         CALL ZGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         DET = (1.D0, 0.D0)
         SGN = 1.D0
         DO i = 1, LDA
         DET = DET * AINV(i,i)
         IF (IPVT(i) .ne. i) THEN
            SGN = -SGN
         ENDIF
         enddo
         DET = SGN * DET
         CALL ZGETRI(LDA, AINV, LDA, IPVT, WORK, LDA, INFO)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_C

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function calculates the LU decomposition and the determinant
!> in a subpart of the complex input matrix.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] AINV a 2D array containing the inverse of the subpart
!> @param[out] DET the determinant of the input matrix.
!> @param[in] Ndim The size of the subpart.
!--------------------------------------------------------------------
       SUBROUTINE INV_C_Variable(A, AINV, DET, Ndim)
         IMPLICIT NONE
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         COMPLEX (KIND=KIND(0.D0)), INTENT(OUT) :: DET
         INTEGER, INTENT(IN) :: Ndim

! Working space.
         REAL (KIND=KIND(0.D0)) :: SGN
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(LDA) )

         AINV = A
         CALL ZGETRF(Ndim, Ndim, AINV, LDA, IPVT, INFO)
         DET = 1.D0
         SGN = 1.D0
         DO i = 1, Ndim
         DET = DET * AINV(i,i)
         IF (IPVT(i) .ne. i) THEN
            SGN = -SGN
         ENDIF
         ENDDO
         DET = SGN * DET
         CALL ZGETRI(Ndim, AINV, LDA, IPVT, WORK, LDA, INFO)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_C_VARIABLE

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function calculates the LU decomposition and the determinant
!> of the complex input matrix.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] AINV a 2D array containing the inverse of the matrix A.
!> @param[out] DET The determinant in a kind of Mantissa-Exponent like
!>                 representation. The full determinant is d1*10^d2
!>                 (1.0 <= |d1| < 10.0) OR (d1 == 0.0)
!
!> @todo the same restrictions as in INV_R_VARIABLE_1 apply.
!--------------------------------------------------------------------
       SUBROUTINE INV_C1(A, AINV, DET)
         IMPLICIT NONE
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(2), INTENT(OUT) :: DET

! Working space.
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, LDA, I
         REAL (KIND=KIND(0.D0)) :: mag

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         AINV = A
         CALL ZGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         DET(1) = (1.D0, 0.D0)
         DET(2) = 0.D0
         DO i = 1, LDA
         mag = ABS(AINV(i, i))
         ! update the phase
         DET(1) = DET(1) * AINV(i, i)/mag
         ! update the magnitude
         DET(2) = DET(2) + LOG10(mag)
         ! consider signs due to permutations
         IF (IPVT(i) .ne. i) THEN
            DET(1) = -DET(1)
         ENDIF
         ENDDO
         CALL ZGETRI(LDA, AINV, LDA, IPVT, WORK, LDA, INFO)
         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_C1
!*****

       SUBROUTINE COMPARE_C(A,B,XMAX,XMEAN)
         IMPLICIT NONE
         COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:), INTENT(in) :: A,B
         REAL (Kind=Kind(0.d0)), INTENT(out) :: XMAX, XMEAN
         INTEGER I,J, N, M

         REAL (Kind=Kind(0.d0)) :: DIFF

         N = SIZE(A,1)
         M = SIZE(A,2)

         XMAX = 0.D0
         XMEAN = 0.D0
         DO I = 1,N
            DO J = 1,M
               DIFF = ABS(A(I,J) - B(I,J))
               IF (DIFF.GT.XMAX) XMAX = DIFF
               XMEAN = XMEAN + DIFF
            ENDDO
         ENDDO
         XMEAN = XMEAN/DBLE(N*M)
       END SUBROUTINE COMPARE_C

       SUBROUTINE COMPARE_R(A,B,XMAX,XMEAN)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)) , INTENT(IN), DIMENSION(:,:) :: A,B
         REAL (Kind=Kind(0.d0)) , INTENT(OUT) :: XMAX, XMEAN
         INTEGER I,J, N, M

         REAL (Kind=Kind(0.d0)) :: DIFF

         N = SIZE(A,1)
         M = SIZE(A,2)

         XMAX = 0.D0
         XMEAN = 0.D0
         DO I = 1,N
            DO J = 1,M
               DIFF = ABS( ( B(I,J) - A(I,J) ) )
               IF (DIFF.GT.XMAX) XMAX = DIFF
               XMEAN = XMEAN + DIFF
            ENDDO
         ENDDO
         XMEAN = XMEAN/DBLE(N*M)
       END SUBROUTINE COMPARE_R

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and Florian Goth
!
!> @brief
!> This function calculates the UDV decomposition using the standard
!> QR algorithm of LaPack.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] U a 2D array containing the left singular vectors.
!> @param[out] D a 1D array containing the sorted singular values.
!> @param[out] V an triangular shaped matrix
!> @param[in] NCON
!--------------------------------------------------------------------
       SUBROUTINE UDV1_R(A,U,D,V,NCON)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:) :: A
         REAL (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
         REAL (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:) :: D
         INTEGER, INTENT(IN) :: NCON

!        The Det of V is not equal to unity.
! Locals:
         INTEGER, DIMENSION(:), ALLOCATABLE :: IVPT, IVPTM1
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: XNORM, VHELP,&
              & TAU, WORK
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: TMP, V1,&
              & TEST, TEST1, TEST2
         REAL (KIND=KIND(0.D0)) :: XMAX, XMEAN, Z
         INTEGER I,J,K, ND1, ND2, NR, IMAX, INFO, LWORK

         ND1 = SIZE(A,1)
         ND2 = SIZE(A,2)



! WRITE(6,*) 'Udv A: ',ND1,ND2
! WRITE(6,*) 'Udv V: ',size(V,1), size(V,2)
! You should now check corresponding sizes for U,V,D.
         IF (SIZE(U,1).NE.ND1 .OR. SIZE(U,2).NE.ND2) THEN
            WRITE(error_unit,*) 'UDV1_R: UDV dim mistake: U'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         ENDIF
         IF (SIZE(D,1).NE.ND2 ) THEN
            WRITE(error_unit,*) 'UDV1_R: UDV dim mistake: D'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         ENDIF
         IF (SIZE(V,1).NE.ND2 .OR. SIZE(V,2).NE.ND2) THEN
            WRITE(error_unit,*) 'UDV1_R: UDV dim mistake: V'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         ENDIF

         ALLOCATE(XNORM (ND2))
         ALLOCATE(VHELP (ND2))
         ALLOCATE(IVPT (ND2))
         ALLOCATE(IVPTM1(ND2))
         ALLOCATE(WORK (ND2))
         ALLOCATE(TAU (ND2))

         ALLOCATE(TMP(ND1,ND2))
         ALLOCATE(V1 (ND2,ND2))

         V1 = 0.D0

         DO I = 1,ND2
            XNORM(I) = 0.D0
            DO NR = 1,ND1
               XNORM(I) = XNORM(I) + ABS(A(NR,I))
            ENDDO
         ENDDO
         VHELP = XNORM

         DO I = 1,ND2
            XMAX = VHELP(1)
            IMAX = 1
            DO J = 2, ND2
               IF (VHELP(J).GT.XMAX) IMAX = J
               IF (VHELP(J).GT.XMAX) XMAX = VHELP(J)
            ENDDO
            VHELP(IMAX) = -1.D0
            IVPTM1(IMAX)=I
            IVPT(I) = IMAX
         ENDDO

         DO I = 1,ND2
            Z = 1.D0/XNORM(IVPT(I))
            K = IVPT(I)
            DO NR = 1,ND1
               TMP(NR,I) = A(NR,K)*Z
            ENDDO
         ENDDO


         !You now want to UDV TMP.
         INFO = 0

        ! Query optimal work space
        CALL DGEQRF(ND1, ND2, TMP, ND1, TAU, WORK, -1, INFO)
        LWORK = INT(WORK(1))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
        CALL DGEQRF(ND1, ND2, TMP, ND1, TAU, WORK, LWORK, INFO)
!         CALL F01QCF(ND1,ND2,TMP,ND1,THETA,INFO)


         !Scale V1 to a unit triangluar matrix.
         DO I = 1,ND2
            D(I) = ABS(TMP(I,I))
         ENDDO
         DO I = 1,ND2
            Z = 1.D0/D(I)
            DO J = I,ND2
               V1(I,J) = TMP(I,J)*Z
            ENDDO
         ENDDO

! Compute U
         INFO = 0
        CALL DORGQR(ND1, ND2, ND2, TMP, ND1, TAU, WORK, LWORK, INFO)
!         CALL F01QEF('Separate', ND1,ND2, ND2, TMP,&
!              & ND1, THETA, WORK, INFO)
        CALL DLACPY('A', ND1, ND2, TMP, ND1, U, Size(U,1))
!
!          DO I = 1,ND1
!             DO J = 1,ND2
!                U(I,J) = TMP(I,J)
!             ENDDO
!          ENDDO
         DEALLOCATE(TMP, WORK, TAU)

! Finish the pivotting.
         DO I = 1,ND2
            D(I) = D(I)*XNORM(IVPT(I))
         ENDDO
         DO I = 1,ND2-1
            Z = 1.D0/XNORM(IVPT(I))
            DO J = I+1,ND2
               V1(I,J) = V1(I,J)*XNORM(IVPT(J))*Z
            ENDDO
         ENDDO

         DO J = 1,ND2
            DO I = 1,ND2
               V(I,J) = V1(I,IVPTM1(J))
            ENDDO
         ENDDO

! Test accuracy.
         IF (NCON.EQ.1) THEN
            ALLOCATE (TEST(ND1,ND2))
            DO J = 1,ND2
               DO I = 1,ND1
                  Z = 0.D0
                  DO NR = 1,ND2
                     Z = Z + U(I,NR)*D(NR)*V(NR,J)
                  ENDDO
                  TEST(I,J) = Z
               ENDDO
            ENDDO
            CALL COMPARE(TEST,A,XMAX,XMEAN)
            WRITE(6,*) 'Accuracy: ',XMAX
            DEALLOCATE (TEST)

            ALLOCATE (TEST (ND2,ND1))
            ALLOCATE (TEST1 (ND2,ND2))
            ALLOCATE (TEST2 (ND2,ND2))
            ! Check orthogonality of U
            DO I = 1,ND1
               DO J = 1,ND2
                  TEST(J,I) = U(I,J)
               ENDDO
            ENDDO
            CALL MMULT(TEST1,TEST,U)
            CALL INITD(TEST2,1.D0)
            CALL COMPARE(TEST1,TEST2,XMAX,XMEAN)
            WRITE(6,*) 'UDV1 orth U: ',XMAX
            DEALLOCATE (TEST )
            DEALLOCATE (TEST1 )
            DEALLOCATE (TEST2 )
         ENDIF


         DEALLOCATE(XNORM )
         DEALLOCATE(VHELP )
         DEALLOCATE(IVPT )
         DEALLOCATE(IVPTM1)
         DEALLOCATE(V1 )

      END SUBROUTINE UDV1_R

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and Florian Goth
!
!> @brief
!> This function calculates a UDV decomposition using the standard
!> QR algorithm of LaPack.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] U a 2D array containing the left singular vectors.
!> @param[out] D a 1D array containing the sorted singular values.
!> @param[out] V an triangular shaped matrix
!> @param[in] NCON
!--------------------------------------------------------------------
      SUBROUTINE UDV_C(A,U,D,V,NCON)
        IMPLICIT NONE
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:) :: D
        INTEGER, INTENT(IN) :: NCON

        !Local
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: TMP, TEST
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: TAU, WORK
        COMPLEX (KIND=KIND(0.D0)) :: Z
        REAL (KIND=KIND(0.D0)) :: DETV, XMDIFF, X
        INTEGER :: NE, LQ, INFO, I, J, NR, LWORK

        LQ = SIZE(A,1)
        NE = SIZE(A,2)

        U = 0.D0 ; V = 0.D0; D = 0.D0
        ALLOCATE (TMP(LQ,NE), TAU(NE), WORK(NE))

        TMP = A

        !You now want to UDV TMP.
        INFO = 0

        ! Query optimal work space. Old style NAG routines use the previously allocated work array
#if !defined(OLDNAG)
#if defined(QRREF)
        CALL ZGEQRF_REF(LQ, NE, TMP, LQ, TAU, WORK, -1, INFO)
#else
        CALL ZGEQRF(LQ, NE, TMP, LQ, TAU, WORK, -1, INFO)
#endif
        LWORK = INT(DBLE(WORK(1)))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
#endif

#if defined(QRREF)
        CALL ZGEQRF_REF(LQ, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#elif !defined(OLDNAG)
        CALL ZGEQRF(LQ, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#else
        CALL F01RCF(LQ,NE,TMP,LQ,TAU,INFO)
#endif
        CALL ZLACPY('U', NE, NE, TMP, LQ, V, Size(V,1))

        DETV = 1.D0
        !V is an NE by NE upper triangular matrix with real diagonal elements.
        DO I = 1,NE
           DETV = DETV * DBLE( TMP(I,I) )
        ENDDO

        !Compute U
! We assume that ZUNGQR and ZGEQRF can work on the same work array.
#if defined(QRREF)
        CALL ZUNGQR_REF(LQ, NE, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#elif !defined(OLDNAG)
        CALL ZUNGQR(LQ, NE, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#else
        CALL F01REF('Separate', LQ,NE, NE, TMP, &
            & LQ, TAU, WORK, INFO)
#endif
        CALL ZLACPY('A', LQ, NE, TMP, LQ, U, Size(U,1))
        DEALLOCATE(TAU, TMP, WORK)
        IF (DBLE(DETV).LT.0.D0) THEN
           DO I = 1,LQ
              U(I,1) = -U(I,1)
           ENDDO
           DO I = 1,NE
              V(1,I) = -V(1,I)
           ENDDO
        ENDIF

        !Scale V1 to a unit triangluar matrix.
        DO I = 1,NE
           X = ABS(DBLE(V(I,I)))
           D(I) = CMPLX(X, 0.D0, kind(0.D0))
           X = 1.D0/X
           DO J = I,NE
              V(I,J) = V(I,J)*X
           ENDDO
        ENDDO

        !Test accuracy.
        IF (NCON.EQ.1) THEN
           ALLOCATE( TEST(LQ,NE) )
           DO J = 1,NE
              DO I = 1,LQ
                 Z = 0.D0
                 DO NR = 1,NE
                    Z = Z + U(I,NR)*D(NR)*V(NR,J)
                 ENDDO
                 TEST(I,J) = Z
              ENDDO
           ENDDO
           XMDIFF = 0.D0
           DO J = 1,LQ
              DO I = 1,NE
                 X = ABS(TEST(J,I)-A(J,I))
                 IF (X.GT.XMDIFF) XMDIFF = X
              ENDDO
           ENDDO
           WRITE(6,*) 'Accuracy, ortho: ',XMDIFF
           DEALLOCATE( TEST )
        ENDIF
        RETURN
      END SUBROUTINE UDV_C

!***************
      SUBROUTINE QR_C(A,U,V,NCON)

        IMPLICIT NONE
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
        INTEGER, INTENT(IN) :: NCON

        !Local
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: TMP, TEST
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: TAU, WORK
        COMPLEX (KIND=KIND(0.D0)) :: Z
        REAL (KIND=KIND(0.D0)) :: DETV, XMDIFF, X
        INTEGER :: NE, LQ, INFO, I, J, LDV, LDU, DU2, DV2, LWORK

        LQ = SIZE(A,1)
        NE = SIZE(A,2)
        LDV = SIZE(V,1)
        LDU = SIZE(U,1)
        DV2 = SIZE(V,2)
        DU2 = SIZE(U,2)
        Z = 0.D0
        call ZLASET('A', LDU, DU2, Z, Z, U, LDU)
        call ZLASET('A', LDV, DV2, Z, Z, V, LDV)
        ALLOCATE (TMP(LQ,NE), TAU(NE), WORK(NE))
        call ZLACPY('A', LQ, NE, A, LQ, TMP, LQ)

        !You now want to UDV TMP. Nag routines.
        INFO = 0

        ! Query optimal work space. Old style NAG routines use the previously allocated work array
#if !defined(OLDNAG)
#if defined(QRREF)
        CALL ZGEQRF_REF(LQ, NE, TMP, LQ, TAU, WORK, -1, INFO)
#else
        CALL ZGEQRF(LQ, NE, TMP, LQ, TAU, WORK, -1, INFO)
#endif
        LWORK = INT(DBLE(WORK(1)))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
#endif
#if defined(QRREF)
        CALL ZGEQRF_REF(LQ, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#elif !defined(OLDNAG)
        CALL ZGEQRF(LQ, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#else
        CALL F01RCF(LQ,NE,TMP,LQ,TAU,INFO)
#endif
        call ZLACPY('U', NE, NE, TMP, LQ, V, LDV)
        DETV = 1.D0
        !V is an NE by NE upper triangular matrix with real diagonal elements.
        DO I = 1,NE
           DETV = DETV * DBLE( TMP(I,I) )
        ENDDO
        !Compute U
! We assume that ZUNGQR and ZGEQRF can work on the same work array.
#if defined(QRREF)
        CALL ZUNGQR_REF(LQ, NE, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#elif !defined(OLDNAG)
        CALL ZUNGQR(LQ, NE, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#else
        CALL F01REF('Separate', LQ,NE, NE, TMP, &
            & LQ, TAU, WORK, INFO)
#endif
        call ZLACPY('A', LQ, NE, TMP, LQ, U, LDU)
        DEALLOCATE(WORK, TAU, TMP)
        IF (DBLE(DETV).LT.0.D0) THEN
           DO I = 1,LQ
              U(I,1) = -U(I,1)
           ENDDO
           DO I = 1,NE
              V(1,I) = -V(1,I)
           ENDDO
        ENDIF

        !Test accuracy.
        IF (NCON.EQ.1) THEN
           ALLOCATE( TEST(LQ,NE) )
           call MMULT(TEST, U, V)
           XMDIFF = 0.D0
           DO J = 1,LQ
              DO I = 1,NE
                 X = ABS(TEST(J,I)-A(J,I))
                 IF (X.GT.XMDIFF) XMDIFF = X
              ENDDO
           ENDDO
           WRITE(6,*) 'Accuracy, QR: ',XMDIFF
           DEALLOCATE( TEST )
        ENDIF
        RETURN
      END SUBROUTINE QR_C

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and Florian Goth
!
!> @brief
!> This function calculates the SVD using the standard QR algorithm
!> of LaPack.
!
!> @note Using the Divide & Conquer algorithm would not yield
!> enough accuracy for using within an auxiliary field type algorithm.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] U a 2D array containing the left singular vectors.
!> @param[out] D a 1D array containing the sorted singular values.
!> @param[out] V a 2D array containing the right singular vectors.
!> @param[in] NCON
!--------------------------------------------------------------------

      SUBROUTINE SVD_C(A,U,D,V,NCON)
        !Uses LaPack Routine
        !#include "machine"

        IMPLICIT NONE
        COMPLEX (KIND=Kind(0.D0)), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=Kind(0.D0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
        COMPLEX (KIND=Kind(0.D0)), INTENT(INOUT), DIMENSION(:) :: D
        INTEGER, INTENT(IN) :: NCON

        !! Local
        REAL    (Kind=Kind(0.D0)), Allocatable :: RWORK(:), S(:)
        COMPLEX (Kind=Kind(0.D0)), Allocatable :: WORK(:), A1(:,:)
        CHARACTER (Len=1):: JOBU,JOBVT
        INTEGER          :: M,N, LDA, LDVT, LDU, LWORK, I, J, I1, INFO
        REAL    (Kind=Kind(0.D0)) :: X, Xmax
        COMPLEX (Kind=Kind(0.D0)) :: Z

        JOBU = "A"
        JOBVT= "A"
        M = SIZE(A,1)
        N = SIZE(A,2)
        Allocate (A1(M,N), S(N))
        A1 = A
        LDA = M
        LDU = M
        LDVT = N
        ALLOCATE( RWORK(5*MIN(M,N)), WORK(10))
! Query optimal amount of memory
        CALL ZGESVD( JOBU, JOBVT, M, N, A1, LDA, S, U, LDU, V, LDVT,&
             &        WORK, -1, RWORK, INFO )
        LWORK = INT(DBLE(WORK(1)))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
        CALL ZGESVD( JOBU, JOBVT, M, N, A1, LDA, S, U, LDU, V, LDVT,&
             &        WORK, LWORK, RWORK, INFO )
        DO I = 1,N
           D(I) = cmplx(S(I), 0.d0, kind(0.D0))
        ENDDO

        IF (NCON ==  1) THEN
           Write(6,*) JobU, JobVT
           Xmax = 0.d0
           DO I = 1,M
              DO I1 = 1,N
                 Z = cmplx(0.d0,0.d0,Kind(0.D0))
                 DO J = 1,N
                    Z  =  Z + U(I,J) *D(J) *V(J,I1)
                 ENDDO
                 X = ABS(Z - A(I,I1))
                 IF (X > Xmax ) Xmax = X
              ENDDO
           ENDDO
           WRITE(6,*) "Success (0), PRE ", INFO, Xmax
        ENDIF


        Deallocate (WORK,RWORK,A1,S)


      END SUBROUTINE SVD_C

!--------------------------------------------------------------------
!> @author
!> Fakher Assaad and  Florian Goth
!
!> @brief
!> This function diagonalizes the input matrix A and returns
!> eigenvalues and vectors using the lapack routine DSYEV.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @param[out] U a 2D array containing the eigen vectors.
!> @param[out] W a 1D array containing the sorted eigenvalues.
!--------------------------------------------------------------------
      SUBROUTINE DIAG_R(A,U,W)
        IMPLICIT NONE
        REAL (KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:) :: A
        REAL (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:,:) :: U
        REAL (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:) :: W


        INTEGER ND1, ND2, IERR, DN
        REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK

         ND1 = SIZE(A,1)
         ND2 = SIZE(A,2)

         IF (ND1.NE.ND2) THEN
            WRITE(error_unit,*) 'DIAG_R: Error in matrix dimension'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         ENDIF

         IERR = 0
         U=A
         W=0
         ! let's just give lapack enough memory
         DN = 3*ND1
         ALLOCATE(WORK(DN))
         CALL DSYEV('V', 'U', ND1, U, ND1, W, WORK, DN, IERR)
         DEALLOCATE(WORK)

      END SUBROUTINE DIAG_R
!*********

      SUBROUTINE DIAG_I(A,U,W)
        ! Uses Lapack
        IMPLICIT NONE
        COMPLEX (Kind=Kind(0.d0)), INTENT(IN)   , DIMENSION(:,:) :: A
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:,:) :: U
        REAL    (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:)   :: W

        CHARACTER (len=1) :: UPLO, JOBZ
        INTEGER :: N, LWORK, INFO
        COMPLEX (Kind=Kind(0.d0)), allocatable :: WORK (:)
        REAL    (Kind=Kind(0.d0)), allocatable :: RWORK(:)
        Logical :: Test
        Integer :: I,J,m
        Complex (Kind=Kind(0.d0)) :: Z
        Real (Kind=Kind(0.d0)) :: X, XMAX

        JOBZ = "V"
        UPLO = "U"
        N = size(A,1)
        U = A
        LWORK = 2*N -1
        Allocate ( WORK(LWORK) )
        Allocate (  RWORK(3*N-2))

        !Write(6,*) 'In Diag'

        Call ZHEEV (JOBZ, UPLO, N, U, N, W, WORK, LWORK, RWORK, INFO)

        Deallocate (WORK, RWORK)

        Test = .false.
        If (Test) then
           XMAX = 0.d0
           DO I = 1,N
              DO J = 1,N
                 Z = cmplx(0.d0,0.d0,Kind=Kind(0.d0))
                 DO m = 1,N
                    Z =  Z + U(I,m)*cmplx(W(m),0.d0, Kind=Kind(0.d0))*Conjg(U(J,m))
                 ENDDO
                 Z = Z - A(I,J)
                 X = ABS(Z)
                 If (X > XMAX ) XMAX = X
              ENDDO
           ENDDO
           write(6,*) ' Test Diag_I: ', XMAX
        endif

      End SUBROUTINE DIAG_I

      SUBROUTINE SECONDS(X)
        IMPLICIT NONE
        REAL (Kind=Kind(0.d0)), INTENT(OUT) :: X

        !DATE_AND_TIME(date, time, zone, values)
        !date_and_time([date][,time][,zone][,values])
        !Subroutine. Die Parameter haben das Attribut intent(out), geben also Werte zurueck.

        ! date: skalare, normale Zeichenvariable von wenigstens 8 Zeichen.
        ! Die linken 8 Zeichen bekommen einen Wert der Form JJJJMMTT . JJJJ Jahr, MM Monat, TT Tag im Monat.
        ! time: skalare, normale Zeichenvariable von wenigstens 10 Zeichen.
        ! Die linken 10 Zeichen bekommen einen Wert der Form hhmmss.sss , wobei hh die Stunde des Tages ist,
        ! mm die Minute innerhalb der Stunde, und ss.sss die Sekunde mit Bruchteilen.
        ! zone: skalare, normale Zeichenvariable von wenigstens 5 Zeichen. Die linken 5 Zeichen bekommen
        ! einen Wert der Form hhmm . hh Stunden, mm Minuten Zeitdifferenz gegenueber der UTC-Weltzeit.
        ! values: Eindimensionales Integer-Feld. Laenge wenigstens 8.
        ! 1 : Jahr, z.B. 1993. 2: Monat. 3: Monatstag. 4: Zeitdifferenz zur Weltzeit in Minuten.
        ! 5: Stunde des Tages. 6: Minute innerhalb der Stunde. 7: Sekunden 8. Millisekunden.

        !character(len=10) :: d,t
        integer,dimension(8) :: V
        !d = ""
        !call date_and_time(date=d,time=t)
        call date_and_time(values=V)

        X = DBLE(V(5)*3600 + V(6)*60 + V(7))

      END SUBROUTINE SECONDS

!====================================================
      Complex (Kind=Kind(0.d0)) Function DET_C(Mat,N)

        Implicit none

        ! Arguments
        Integer, intent(in) :: N
        Complex(Kind=Kind(0.d0)), intent(inout) :: mat(:,:)

        integer :: i, info, LDmat
        integer, allocatable :: ipiv(:)

        integer :: sgn

        allocate(ipiv(N))
        ipiv = 0

        !Lapack LU decomposition
        LDmat=size(mat,1)
        call zgetrf(N, N, mat, LDmat, ipiv, info)

        det_C = cmplx(1.d0, 0.d0, kind(0.d0) )
        do i = 1, N
           det_C = det_C*mat(i, i)
        enddo

        sgn =  1
        do i = 1, N
           if(ipiv(i) /= i)  sgn = -sgn
        enddo
        if (sgn == -1 ) det_C = - det_C

      end function DET_C

!--------------------------------------------------------------------
!> @author
!> F.Assaad
!
!> @brief
!> Returns the determinant as det = \prod_i=1^N d(i).
!> Uses Lapack LU decomposition
!>
!====================================================
      Subroutine DET_C_LU(Mat1,D,N)

        Implicit none

        ! Arguments
        Integer, intent(in) :: N
        Complex(Kind=Kind(0.d0)), intent(in)  :: mat1(N,N)
        Complex(Kind=Kind(0.d0)), intent(out) :: D(N)


        Complex(Kind=Kind(0.d0)) :: mat(N,N)

        integer :: i, info
        integer :: ipiv(N)

        integer :: sgn

        mat = mat1
        ipiv = 0

        !Lapack LU decomposition
        call zgetrf(N, N, mat, N, ipiv, info)

        do i = 1,N
           D(i) = mat(i,i)
        enddo
        sgn =  1
        do i = 1, N
           if(ipiv(i) /= i)  sgn = -sgn
        enddo
        if (sgn == -1 ) D(1) = - D(1)

      end Subroutine DET_C_LU

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief
!> This function returns the Conjugate Transpose, hence the name CT,
!> of a complex input matrix.
!
!> @note The employed variant of chaining a transposition and a conjugation
!> is well optimized by gcc in the sense that it generates a tight inner loop.
!> The original double loop generates quite some additional integer operations.
!
!> @param[in] A a 2D array constituting the input matrix.
!> @return the Conjugate Transpose of A
!--------------------------------------------------------------------
    function ct(a) result(b) ! return the conjugate transpose of a matrix
        complex(kind=kind(0.D0)), dimension(:,:), intent(in) :: a
        complex(kind=kind(0.D0)), dimension(size(a,1),size(a,1)) :: b
        b = conjg(transpose(a))
    end function ct

    END MODULE MyMats
