!  Copyright (C) 2016-2026 The ALF project
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
!> @file mymats_mod.F90
!> @author ALF-project
!> @brief Wrappers for dense linear algebra operations via LAPACK/BLAS.
!> @details This module provides unified interfaces for common linear algebra
!> tasks used in QMC simulations:
!>  - Matrix multiplication (MMULT): generic dispatching to DGEMM/ZGEMM
!>  - Matrix inversion (INV): LU decomposition via DGETRF/ZGETRF + determinant
!>  - Eigenvalue decomposition (DIAG): symmetric/Hermitian via DSYEV/ZHEEV
!>  - Generalized eigenvalue (DIAG_GEN): left/right eigenvectors via ZGEEV
!>  - QR decomposition (QR): Householder via ZGEQRF/ZUNGQR
!>  - SVD (SVD): singular value decomposition via ZGESVD
!>  - UDV decomposition (UDV): pivoted QR for numerical stability
!>  - Determinant (DET): LU-based via ZGETRF and DET_C_LU
!>  - Utilities: matrix initialization (INITD), comparison (COMPARE), timing (SECONDS)
!>
!> Most routines offer real (DGEMM-based) and complex (ZGEMM-based) variants,
!> dispatched via Fortran generic interfaces. Sign-flipping logic handles
!> determinant sign corrections from row pivoting in LU decomposition.
!> @see LAPACK documentation for DGETRF, ZGETRF, DGEMM, ZGEMM, ZHEEV, ZGEEV, etc.
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
!> @brief Generalized eigenvalue decomposition for non-Hermitian matrices.
!> @details Solves A·U = U·W (right) or U·A = W·U (left) for a complex matrix
!> via LAPACK ZGEEV. Computes left eigenvectors (LR='L'), right eigenvectors
!> (LR='R'), or optionally neither. Conjugate transpose applied to left
!> eigenvectors for proper adjoint interpretation.
!>  - Left:  U·A = W·U, computed as CONJG(V^T) where V from ZGEEV
!>  - Right: A·U = U·W, computed directly as V from ZGEEV
!>
!> @param[in] Z_MAT (N×N) complex input matrix to diagonalize
!> @param[out] U (N×N) eigenvector matrix (left if LR='L', right if LR='R')
!> @param[out] W (N) array of eigenvalues (same for left/right)
!> @param[in] LR 'L'/'l'=left eigenvectors, 'R'/'r'=right eigenvectors
!> @param[in] ICON 0=skip validation, 1=test accuracy via residual
!>
!> @note Uses LAPACK ZGEEV generalized eigenvalue solver (QR algorithm)
!> @see ZGEEV documentation and QMC context for non-Hermitian operators
!--------------------------------------------------------------------
       SUBROUTINE DIAG_GEN(Z_MAT,U,W,LR,ICON)
         IMPLICIT NONE
         COMPLEX   (Kind=Kind(0.d0)), INTENT(IN), DIMENSION(:,:) :: Z_MAT
         CHARACTER (LEN=1),  INTENT(IN)  :: LR
         COMPLEX   (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:,:) :: U
         COMPLEX   (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:) :: W
         INTEGER :: ICON

         ! Uses LAPACK ZGEEV for non-Hermitian eigenvalue decomposition
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

         ! Get matrix dimension
         N = SIZE(Z_MAT,1)
         ALLOCATE(A(N,N))
         A = Z_MAT  ! Working copy (ZGEEV modifies input)
         LDA = N

         ! Initialize: no eigenvectors by default
         JOBVR  = "N"
         JOBVL  = "N"
         LDVL = 1
         LDVR = 1
         ! Determine which eigenvectors to compute based on LR flag
         IF (LR == "L" .or. LR == 'l') THEN
            JOBVL ="V"  ! Compute left eigenvectors
            LDVL  = N
         ELSEIF (LR == "R" .or. LR == 'r') THEN
            JOBVR ="V"  ! Compute right eigenvectors
            LDVR = N
         ELSE
            WRITE(error_unit,*) 'Error in DIAG_GEN'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         ENDIF
         ALLOCATE(VL(LDVL,N),  VR(LDVR,N) )
         LWORK = 2*N
         ALLOCATE (WORK(LWORK), RWORK(LWORK) )

         ! Call LAPACK ZGEEV: compute eigenvalues and optionally eigenvectors
         ! W = eigenvalues (complex), VL/VR = left/right eigenvectors (if requested)
         CALL ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
              &      WORK, LWORK, RWORK, INFO )

         ! Copy eigenvectors: right eigenvectors are used directly,
         ! left eigenvectors are conjugate-transposed for proper adjoint form
         IF (LR == "R" .or. LR == 'r')  THEN
            DO I = 1,N
               DO J = 1,N
                  U(I,J) = VR(I,J)  ! Right eigenvectors (column vectors)
               ENDDO
            ENDDO
         ELSE
            DO I = 1,N
               DO J = 1,N
                  U(I,J) = CONJG(VL(J,I))  ! Left eigenvectors as rows (conjugate transpose)
               ENDDO
            ENDDO
         ENDIF

         ! Optional accuracy test: compute residual max|A·U - U·W| (right) or |U·A - W·U| (left)
         IF (ICON == 1 ) THEN
            XMAX = 0.d0
            DO I = 1,N
               DO J = 1,N
                  IF (LR == "R" .or. LR == 'r')  THEN
                     ! Right: residual = (A·U - U·W)(I,J) = sum_M(A(I,M)·U(M,J)) - W(I)·U(I,J)
                     Z = - W(I)*U(I,J)
                     DO M = 1,N
                        Z = Z + Z_MAT(I,M)*U(M,J)
                     ENDDO
                  ELSE
                     ! Left: residual = (U·A - W·U)(I,J) = sum_M(U(I,M)·A(M,J)) - W(I)·U(I,J)
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
         ENDIF

         DEALLOCATE(VL, VR)
         DEALLOCATE(WORK, RWORK)
         DEALLOCATE(A)


       END SUBROUTINE DIAG_GEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!> @brief Real matrix multiplication: C = A·B via BLAS DGEMM.
!> @details Dispatches to DGEMM with default parameters (no transpose,
!> C := 1.0·A·B + 0.0·C). Handles column-major leading dimensions.
!>
!> @param[out] C (N×P) output product matrix
!> @param[in] A (N×M) left input matrix
!> @param[in] B (M×P) right input matrix
!> @note Uses BLAS DGEMM for high performance (DGEMM does not require square matrices)
!> @see DGEMM documentation
!--------------------------------------------------------------------
       SUBROUTINE MMULT_R(C, A, B)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: A,B,C
         REAL (Kind=Kind(0.d0)) :: ALP, BET
         INTEGER N, M, P, LDA, LDB, LDC
         
         ! Extract dimensions: A is N×M, B is M×P, C is N×P
         N = SIZE(A,1)      ! Rows in A
         M = SIZE(A,2)      ! Columns in A = rows in B
         P = SIZE(B,2)      ! Columns in B
         LDA = N; LDB = SIZE(B,1); LDC = SIZE(C,1)  ! Leading dimensions

         ! DGEMM coefficients: C := alpha·A·B + beta·C
         ALP = 1.D0         ! Coefficient for A·B
         BET = 0.D0         ! Coefficient for C (overwrites C)

         ! Call BLAS DGEMM: real matrix multiplication
         CALL DGEMM('n','n',N,P,M,ALP,A,LDA,B,LDB,BET,C,LDC)
       END SUBROUTINE MMULT_R

!--------------------------------------------------------------------
!> @brief Complex matrix multiplication: C = A·B via BLAS ZGEMM.
!> @details Dispatches to ZGEMM with default parameters (no transpose,
!> C := 1.0·A·B + 0.0·C). Handles column-major leading dimensions.
!>
!> @param[out] C (N×P) output product matrix
!> @param[in] A (N×M) left input complex matrix
!> @param[in] B (M×P) right input complex matrix
!> @note Uses BLAS ZGEMM for high performance
!> @see ZGEMM documentation
!--------------------------------------------------------------------
       SUBROUTINE MMULT_C(C, A, B)
         IMPLICIT NONE
         COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:) :: A,B,C
         COMPLEX (Kind=Kind(0.d0)) :: ALP, BET
         INTEGER N, M, P, LDA, LDB, LDC

         ! Extract dimensions: A is N×M, B is M×P, C is N×P
         N = SIZE(A,1)      ! Rows in A
         M = SIZE(A,2)      ! Columns in A = rows in B
         P = SIZE(B,2)      ! Columns in B
         LDA = N; LDB = SIZE(B,1); LDC = SIZE(C,1)  ! Leading dimensions

         ! ZGEMM coefficients: C := alpha·A·B + beta·C
         ALP = CMPLX(1.D0,0.D0,Kind=Kind(0d0))   ! Coefficient for A·B
         BET = CMPLX(0.D0,0.D0,Kind=Kind(0d0))   ! Coefficient for C (overwrites C)

         ! Call BLAS ZGEMM: complex matrix multiplication
         CALL ZGEMM('n','n',N,P,M,ALP,A,LDA,B,LDB,BET,C,LDC)

       END SUBROUTINE MMULT_C

!--------------------------------------------------------------------
!> @brief Initialize real matrix to diagonal with uniform diagonal entries.
!> @details Sets all off-diagonal elements to zero, diagonal to scalar value X.
!> Overwrites input array. Used for identity matrices (X=1.0) and testing.
!>
!> @param[inout] A (N×M) real input matrix to initialize
!> @param[in] X scalar diagonal value assigned to A(i,i) for i=1..min(N,M)
!> @note Matrix is assumed square; behavior for non-square is undefined
!--------------------------------------------------------------------
       SUBROUTINE INITD_R(A,X)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: A
         REAL (Kind=Kind(0.d0)) X
         INTEGER I,J, N, M

         ! Get matrix dimensions
         N = SIZE(A,1)
         M = SIZE(A,2)

         ! Zero all elements
         DO I = 1,N
         DO J = 1,M
            A(I,J) = 0.D0
         ENDDO
         ENDDO
         
         ! Set diagonal to X (valid for min(N,M) diagonal elements)
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

         ! Zero entire matrix with complex zero
         A = (0.D0, 0.D0)
         
         ! Get matrix size (assumes square)
         N = SIZE(A,1)
         
         ! Set diagonal to X
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
       SUBROUTINE INV_R0(A, AINV, DET_VAL)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), INTENT(OUT) :: DET_VAL
         INTEGER I

         ! Working space for LU decomposition
         REAL (KIND=KIND(0.D0)) :: SGN
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT  ! Pivot indices
         INTEGER INFO, LDA

         LDA = SIZE(A,1)
         ! Allocate workspace: IPVT for pivot info, WORK for inverse computation
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         ! Copy input to output (DGETRF/DGETRI modify in-place)
         AINV = A
         
         ! Step 1: LU decomposition via LAPACK DGETRF (no pivoting = in-place LU)
         CALL DGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         
         ! Step 2: Extract determinant from LU diagonal and count row swaps
         DET_VAL = 1.D0    ! Start with det = 1
         SGN = 1.D0    ! Sign correction for row pivots
         DO i = 1, LDA
            ! Determinant is product of diagonal elements from LU
            DET_VAL = DET_VAL * AINV(i,i)
            ! If IPVT(i) ≠ i, a row swap occurred (changes sign)
            IF (IPVT(i) .ne. i) THEN
               SGN = -SGN
            ENDIF
         enddo
         ! Apply sign correction: det(A) = (-1)^(# swaps) × product(diagonal)
         DET_VAL = SGN * DET_VAL
         
         ! Step 3: Compute inverse from LU decomposition via DGETRI
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
       SUBROUTINE INV_R_Variable(A, AINV, DET_VAL, Ndim)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), INTENT(OUT) :: DET_VAL
         INTEGER, INTENT(IN) :: Ndim

         ! Working space for LU decomposition (variable dimension subblock)
         REAL (KIND=KIND(0.D0)) :: SGN
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT  ! Pivot indices
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)  ! Full allocated size
         ! Allocate workspace for Ndim×Ndim subblock
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(LDA) )

         ! Copy full array (invert only upper-left Ndim×Ndim block)
         AINV = A
         
         ! Step 1: LU decomposition of Ndim×Ndim subblock in AINV
         CALL DGETRF(Ndim, Ndim, AINV, LDA, IPVT, INFO)
         
         ! Step 2: Extract determinant from diagonal of Ndim×Ndim LU block
         DET_VAL = 1.D0
         SGN = 1.D0
         DO i = 1, Ndim
            ! Determinant contribution from diagonal
            DET_VAL = DET_VAL * AINV(i,i)
            ! Sign correction from row pivots
            IF (IPVT(i) .ne. i) THEN
               SGN = -SGN
            ENDIF
         ENDDO
         DET_VAL = SGN * DET_VAL
         
         ! Step 3: Compute inverse of Ndim×Ndim subblock via DGETRI
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
       SUBROUTINE INV_R_Variable_1(A,AINV,DET_VAL,Ndim)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), DIMENSION(2), INTENT(OUT) :: DET_VAL
         INTEGER, INTENT(IN) :: Ndim

         ! Working space for LU decomposition (mantissa-exponent format)
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT  ! Pivot indices
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)
         ! Allocate workspace for Ndim×Ndim subblock
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(LDA) )

         ! Copy full array
         AINV = A
         
         ! Step 1: LU decomposition of Ndim×Ndim subblock
         CALL DGETRF(Ndim, Ndim, AINV, LDA, IPVT, INFO)
         
         ! Step 2: Extract determinant in mantissa-exponent form: DET_VAL = DET_VAL(1)·10^DET_VAL(2)
         ! DET_VAL(1) = mantissa (±[1,10)), DET_VAL(2) = log10 of magnitude (allows extreme values)
         DET_VAL(1) = 1.D0       ! Initialize mantissa
         DET_VAL(2) = 0.D0       ! Initialize exponent (log10 scale)
         DO i = 1, Ndim
            ! Extract sign from diagonal element
            IF (AINV(i, i) < 0.D0) THEN
               DET_VAL(1) = -DET_VAL(1)  ! Update sign if diagonal negative
            ENDIF
            ! Accumulate log10 of magnitude to avoid overflow/underflow
            DET_VAL(2) = DET_VAL(2) + LOG10(ABS(AINV(i,i)))
            ! Sign correction from row pivots
            IF (IPVT(i) .ne. i) THEN
               DET_VAL(1) = -DET_VAL(1)
            ENDIF
         ENDDO
         
         ! Step 3: Compute inverse of Ndim×Ndim subblock via DGETRI
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
       SUBROUTINE INV_R1(A,AINV,DET_VAL)
         IMPLICIT NONE
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         REAL (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         REAL (KIND=KIND(0.D0)), DIMENSION(2), INTENT(OUT) :: DET_VAL

         ! Working space for LU decomposition (mantissa-exponent format, full matrix)
         REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT  ! Pivot indices
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)
         ! Allocate workspace for full LDA×LDA matrix
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         ! Copy input matrix (DGETRF/DGETRI modify in-place)
         AINV = A
         
         ! Step 1: LU decomposition via LAPACK DGETRF
         CALL DGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         
         ! Step 2: Extract determinant in mantissa-exponent form: DET_VAL = DET_VAL(1)·10^DET_VAL(2)
         ! Avoids overflow/underflow by storing log10 of magnitude separately
         DET_VAL(1) = 1.D0       ! Initialize mantissa
         DET_VAL(2) = 0.D0       ! Initialize exponent (log10 scale)
         DO i = 1, LDA
            ! Extract sign from diagonal element of LU matrix
            IF (AINV(i, i) < 0.D0) THEN
               DET_VAL(1) = -DET_VAL(1)
            ENDIF
            ! Accumulate log10 of magnitude
            DET_VAL(2) = DET_VAL(2) + LOG10(ABS(AINV(i,i)))
            ! Sign correction from row pivots
            IF (IPVT(i) .ne. i) THEN
               DET_VAL(1) = -DET_VAL(1)
            ENDIF
         ENDDO
         
         ! Step 3: Compute inverse from LU decomposition via DGETRI
         CALL DGETRI(LDA, AINV, LDA, IPVT, WORK, LDA, INFO)
         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R1

!--------------------------------------------------------------------
!> @brief Real matrix inversion (no determinant) via LAPACK DGETRF/DGETRI.
!> @details LU decomposition without determinant extraction, for efficiency
!> when determinant is not needed. Uses in-place computation via DGETRF/DGETRI.
!>
!> @param[in] A (N×N) real input matrix to invert
!> @param[out] AINV (N×N) real output inverse matrix
!> @note Uses LAPACK DGETRF and DGETRI
!--------------------------------------------------------------------
       SUBROUTINE INV_R2(A,AINV)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: A,AINV

         INTEGER I, J

         ! Uses LAPACK LU-based inversion routines

! Working space.
         REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
         INTEGER INFO, LDA, LWORK

         LDA = SIZE(A,1)

         ! Allocate pivot array and workspace
         ALLOCATE ( IPIV(LDA) )
         LWORK = LDA
         ALLOCATE ( WORK(LWORK) )
         WORK = 0.0
         IPIV = 0
         
         ! Transpose copy into AINV (note: column-major ordering)
         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO
         INFO = 0

         ! Step 1: LU decomposition via DGETRF (no determinant extraction)
         CALL DGETRF( LDA, LDA, AINV, LDA, IPIV, INFO )
         
         ! Step 2: Compute inverse from LU via DGETRI
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

       SUBROUTINE INV_C(A,AINV,DET_VAL)
         IMPLICIT NONE
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         COMPLEX (KIND=KIND(0.D0)), INTENT(OUT) :: DET_VAL
         INTEGER I

         ! Working space for complex LU decomposition
         REAL (KIND=KIND(0.D0)) :: SGN
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT  ! Pivot indices
         INTEGER INFO, LDA

         LDA = SIZE(A,1)
         ! Allocate workspace
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         ! Copy input to output (ZGETRF/ZGETRI modify in-place)
         AINV = A
         
         ! Step 1: LU decomposition via LAPACK ZGETRF
         CALL ZGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         
         ! Step 2: Extract determinant from LU diagonal and count row swaps
         DET_VAL = (1.D0, 0.D0)   ! Start with det = 1
         SGN = 1.D0           ! Sign correction for row pivots
         DO i = 1, LDA
            ! Determinant is product of diagonal elements from LU decomposition
            DET_VAL = DET_VAL * AINV(i,i)
            ! If IPVT(i) ≠ i, a row swap occurred (changes sign)
            IF (IPVT(i) .ne. i) THEN
               SGN = -SGN
            ENDIF
         enddo
         ! Apply sign correction: det(A) = (-1)^(# swaps) × product(diagonal)
         DET_VAL = SGN * DET_VAL
         
         ! Step 3: Compute inverse from LU decomposition via ZGETRI
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
       SUBROUTINE INV_C_Variable(A, AINV, DET_VAL, Ndim)
         IMPLICIT NONE
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         COMPLEX (KIND=KIND(0.D0)), INTENT(OUT) :: DET_VAL
         INTEGER, INTENT(IN) :: Ndim

         ! Working space for complex LU decomposition (variable dimension subblock)
         REAL (KIND=KIND(0.D0)) :: SGN
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT  ! Pivot indices
         INTEGER INFO, LDA, I

         LDA = SIZE(A,1)  ! Full allocated size
         ! Allocate workspace for Ndim×Ndim subblock
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(LDA) )

         ! Copy full array (invert only upper-left Ndim×Ndim block)
         AINV = A
         
         ! Step 1: LU decomposition of Ndim×Ndim subblock
         CALL ZGETRF(Ndim, Ndim, AINV, LDA, IPVT, INFO)
         
         ! Step 2: Extract determinant from diagonal of Ndim×Ndim LU block
         DET_VAL = 1.D0
         SGN = 1.D0
         DO i = 1, Ndim
            ! Determinant contribution from diagonal
            DET_VAL = DET_VAL * AINV(i,i)
            ! Sign correction from row pivots
            IF (IPVT(i) .ne. i) THEN
               SGN = -SGN
            ENDIF
         ENDDO
         DET_VAL = SGN * DET_VAL
         
         ! Step 3: Compute inverse of Ndim×Ndim subblock via ZGETRI
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
       SUBROUTINE INV_C1(A, AINV, DET_VAL)
         IMPLICIT NONE
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(IN) :: A
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), INTENT(INOUT) :: AINV
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(2), INTENT(OUT) :: DET_VAL

         ! Working space for complex LU decomposition (mantissa-exponent format)
         ! DET returned as: DET(1) = unit-magnitude phase, DET(2) = log10(magnitude)
         COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT  ! Pivot indices
         INTEGER INFO, LDA, I
         REAL (KIND=KIND(0.D0)) :: mag

         LDA = SIZE(A,1)
         ! Allocate workspace
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )

         ! Copy input to output (ZGETRF/ZGETRI modify in-place)
         AINV = A
         
         ! Step 1: LU decomposition via LAPACK ZGETRF
         CALL ZGETRF(LDA, LDA, AINV, LDA, IPVT, INFO)
         
         ! Step 2: Extract determinant in mantissa-exponent form with phase separation
         ! DET_VAL(1) = phase (unit magnitude, |DET_VAL(1)|=1), DET_VAL(2) = log10(|det|)
         DET_VAL(1) = (1.D0, 0.D0)  ! Initialize phase to 1
         DET_VAL(2) = 0.D0          ! Initialize log10 of magnitude
         DO i = 1, LDA
            mag = ABS(AINV(i, i))  ! Magnitude of diagonal element
            ! Update phase: normalize by magnitude to keep |DET_VAL(1)| = 1
            DET_VAL(1) = DET_VAL(1) * AINV(i, i)/mag
            ! Update log10 of magnitude: accumulate contributions
            DET_VAL(2) = DET_VAL(2) + LOG10(mag)
            ! Sign correction from row pivots (affects phase)
            IF (IPVT(i) .ne. i) THEN
               DET_VAL(1) = -DET_VAL(1)
            ENDIF
         ENDDO
         
         ! Step 3: Compute inverse from LU decomposition via ZGETRI
         CALL ZGETRI(LDA, AINV, LDA, IPVT, WORK, LDA, INFO)
         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_C1
!--------------------------------------------------------------------
!> @brief Compare two complex matrices element-wise for differences.
!> @details Computes maximum and mean absolute differences between A and B.
!> Used for validation of numerical algorithms by comparing against reference.
!>
!> @param[in] A (N×M) complex input matrix
!> @param[in] B (N×M) complex reference matrix
!> @param[out] XMAX maximum absolute difference |A(i,j) - B(i,j)|
!> @param[out] XMEAN mean absolute difference averaged over all N×M elements
!--------------------------------------------------------------------
       SUBROUTINE COMPARE_C(A,B,XMAX,XMEAN)
         IMPLICIT NONE
         COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:), INTENT(in) :: A,B
         REAL (Kind=Kind(0.d0)), INTENT(out) :: XMAX, XMEAN
         INTEGER I,J, N, M

         REAL (Kind=Kind(0.d0)) :: DIFF

         ! Get matrix dimensions
         N = SIZE(A,1)
         M = SIZE(A,2)

         ! Initialize accumulators
         XMAX = 0.D0    ! Maximum absolute difference
         XMEAN = 0.D0   ! Sum of absolute differences (will normalize below)
         
         ! Compare all elements
         DO I = 1,N
            DO J = 1,M
               DIFF = ABS(A(I,J) - B(I,J))  ! Absolute difference at (i,j)
               IF (DIFF.GT.XMAX) XMAX = DIFF  ! Track maximum
               XMEAN = XMEAN + DIFF        ! Accumulate for mean
            ENDDO
         ENDDO
         
         ! Normalize mean by total number of elements
         XMEAN = XMEAN/DBLE(N*M)
       END SUBROUTINE COMPARE_C

!--------------------------------------------------------------------
!> @brief Compare two real matrices element-wise for differences.
!> @details Computes maximum and mean absolute differences between A and B.
!> Used for validation of numerical algorithms by comparing against reference.
!>
!> @param[in] A (N×M) real input matrix
!> @param[in] B (N×M) real reference matrix
!> @param[out] XMAX maximum absolute difference |B(i,j) - A(i,j)|
!> @param[out] XMEAN mean absolute difference averaged over all N×M elements
!--------------------------------------------------------------------
       SUBROUTINE COMPARE_R(A,B,XMAX,XMEAN)
         IMPLICIT NONE
         REAL (Kind=Kind(0.d0)) , INTENT(IN), DIMENSION(:,:) :: A,B
         REAL (Kind=Kind(0.d0)) , INTENT(OUT) :: XMAX, XMEAN
         INTEGER I,J, N, M

         REAL (Kind=Kind(0.d0)) :: DIFF

         ! Get matrix dimensions
         N = SIZE(A,1)
         M = SIZE(A,2)

         ! Initialize accumulators
         XMAX = 0.D0    ! Maximum absolute difference
         XMEAN = 0.D0   ! Sum of absolute differences (will normalize below)
         
         ! Compare all elements: B - A
         DO I = 1,N
            DO J = 1,M
               DIFF = ABS( ( B(I,J) - A(I,J) ) )  ! Absolute difference
               IF (DIFF.GT.XMAX) XMAX = DIFF      ! Track maximum
               XMEAN = XMEAN + DIFF               ! Accumulate for mean
            ENDDO
         ENDDO
         
         ! Normalize mean by total number of elements
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
!--------------------------------------------------------------------
!> @brief QR decomposition of complex matrix: A = Q·R.
!> @details Performs Householder QR factorization via LAPACK ZGEQRF.
!> Returns Q (orthogonal) and R (upper triangular) such that A = Q·R.
!> Uses optional sign flipping to normalize determinant of R.
!>
!> @param[in] A (LQ×NE) complex input matrix to decompose
!> @param[out] U (LQ×NE) complex unitary matrix Q
!> @param[out] V (NE×NE) complex upper triangular matrix R
!> @param[in] NCON 0=skip validation, 1=test accuracy via residual
!> @note Uses LAPACK ZGEQRF and ZUNGQR for Householder QR
!--------------------------------------------------------------------
      SUBROUTINE QR_C(A,U,V,NCON)

        IMPLICIT NONE
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
        INTEGER, INTENT(IN) :: NCON

        ! Local workspace variables for QR decomposition
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: TMP, TEST
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: TAU, WORK
        COMPLEX (KIND=KIND(0.D0)) :: Z
        REAL (KIND=KIND(0.D0)) :: DETV, XMDIFF, X
        INTEGER :: NE, LQ, INFO, I, J, LDV, LDU, DU2, DV2, LWORK

        ! Get matrix dimensions and allocate workspace
        LQ = SIZE(A,1)   ! Number of rows (matrix height)
        NE = SIZE(A,2)   ! Number of columns (matrix width)
        LDV = SIZE(V,1)
        LDU = SIZE(U,1)
        DV2 = SIZE(V,2)
        DU2 = SIZE(U,2)
        
        ! Initialize output matrices to zero
        Z = 0.D0
        call ZLASET('A', LDU, DU2, Z, Z, U, LDU)  ! U = 0
        call ZLASET('A', LDV, DV2, Z, Z, V, LDV)  ! V = 0
        
        ! Allocate workspace for QR factorization
        ALLOCATE (TMP(LQ,NE), TAU(NE), WORK(NE))
        ! Copy input matrix to working array
        call ZLACPY('A', LQ, NE, A, LQ, TMP, LQ)

        ! QR decomposition via LAPACK ZGEQRF (Householder QR)
        INFO = 0

        ! Query optimal work space (for efficiency)
#if !defined(OLDNAG)
#if defined(QRREF)
        CALL ZGEQRF_REF(LQ, NE, TMP, LQ, TAU, WORK, -1, INFO)  ! Query
#else
        CALL ZGEQRF(LQ, NE, TMP, LQ, TAU, WORK, -1, INFO)  ! Query
#endif
        LWORK = INT(DBLE(WORK(1)))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
#endif
        
        ! Perform actual QR decomposition: TMP = Q·R (Q stored as Householder reflectors)
#if defined(QRREF)
        CALL ZGEQRF_REF(LQ, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#elif !defined(OLDNAG)
        CALL ZGEQRF(LQ, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#else
        CALL F01RCF(LQ,NE,TMP,LQ,TAU,INFO)
#endif
        
        ! Extract upper triangular R matrix (stored in TMP(1:NE, 1:NE))
        call ZLACPY('U', NE, NE, TMP, LQ, V, LDV)
        ! Compute determinant of R from diagonal
        DETV = 1.D0
        ! V is NE×NE upper triangular with complex diagonal (product of diagonal = det(R))
        DO I = 1,NE
           DETV = DETV * DBLE( TMP(I,I) )  ! Product of diagonal: det(R)
        ENDDO
        
        ! Compute Q from Householder reflectors stored in TMP via ZUNGQR
        ! We reuse the same work array as ZGEQRF above
#if defined(QRREF)
        CALL ZUNGQR_REF(LQ, NE, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#elif !defined(OLDNAG)
        CALL ZUNGQR(LQ, NE, NE, TMP, LQ, TAU, WORK, LWORK, INFO)
#else
        CALL F01REF('Separate', LQ,NE, NE, TMP, &
            & LQ, TAU, WORK, INFO)
#endif
        ! Copy orthogonal Q matrix to output
        call ZLACPY('A', LQ, NE, TMP, LQ, U, LDU)
        DEALLOCATE(WORK, TAU, TMP)
        
        ! Normalize sign of R: if det(R) < 0, flip first column of Q and first row of R
        IF (DBLE(DETV).LT.0.D0) THEN
           DO I = 1,LQ
              U(I,1) = -U(I,1)      ! Flip first column of U (Q)
           ENDDO
           DO I = 1,NE
              V(1,I) = -V(1,I)      ! Flip first row of V (R)
           ENDDO
        ENDIF

        ! Optional accuracy test: compute max|U·V - A| to verify decomposition
        IF (NCON.EQ.1) THEN
           ALLOCATE( TEST(LQ,NE) )
           call MMULT(TEST, U, V)  ! TEST = U × V = Q × R
           XMDIFF = 0.D0
           ! Compare with original A: max|Q·R - A|
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

!--------------------------------------------------------------------
!> @brief Singular Value Decomposition (SVD) of complex matrix.
!> @details Computes A = U·Σ·V^H where U (left), Σ (singular values),
!> V (right) via LAPACK ZGESVD QR algorithm. Singular values sorted in
!> descending order. Used for numerical stability and rank detection.
!>
!> @param[in] A (M×N) complex input matrix
!> @param[out] U (M×M) complex left singular vectors (orthogonal)
!> @param[out] D (N) complex singular values (real, non-negative, on diagonal)
!> @param[out] V (N×N) complex right singular vectors (orthogonal)
!> @param[in] NCON 0=skip validation, 1=test accuracy via residual
!> @note Uses LAPACK ZGESVD with QR algorithm (not Divide-Conquer for accuracy)
!> @see ZGESVD documentation
!--------------------------------------------------------------------
      SUBROUTINE SVD_C(A,U,D,V,NCON)

        IMPLICIT NONE
        COMPLEX (KIND=Kind(0.D0)), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=Kind(0.D0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
        COMPLEX (KIND=Kind(0.D0)), INTENT(INOUT), DIMENSION(:) :: D
        INTEGER, INTENT(IN) :: NCON

        ! Local workspace variables for SVD decomposition
        REAL    (Kind=Kind(0.D0)), Allocatable :: RWORK(:), S(:)
        COMPLEX (Kind=Kind(0.D0)), Allocatable :: WORK(:), A1(:,:)
        CHARACTER (Len=1):: JOBU,JOBVT
        INTEGER          :: M,N, LDA, LDVT, LDU, LWORK, I, J, I1, INFO
        REAL    (Kind=Kind(0.D0)) :: X, Xmax
        COMPLEX (Kind=Kind(0.D0)) :: Z

        ! Request full U and V matrices (JOBU/JOBVT = 'A')
        JOBU = "A"
        JOBVT= "A"
        M = SIZE(A,1)   ! Matrix height
        N = SIZE(A,2)   ! Matrix width
        Allocate (A1(M,N), S(N))  ! Working copies
        A1 = A          ! ZGESVD modifies input, so use working copy
        LDA = M
        LDU = M
        LDVT = N
        
        ! Allocate workspace for SVD computation
        ALLOCATE( RWORK(5*MIN(M,N)), WORK(10))
        
        ! Query optimal workspace size for efficiency
        CALL ZGESVD( JOBU, JOBVT, M, N, A1, LDA, S, U, LDU, V, LDVT,&
             &        WORK, -1, RWORK, INFO )
        LWORK = INT(DBLE(WORK(1)))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
        
        ! Perform SVD: A = U·Σ·V^H with singular values in S
        CALL ZGESVD( JOBU, JOBVT, M, N, A1, LDA, S, U, LDU, V, LDVT,&
             &        WORK, LWORK, RWORK, INFO )
        
        ! Convert real singular values S to complex D for output
        DO I = 1,N
           D(I) = cmplx(S(I), 0.d0, kind(0.D0))  ! S(i) are real and non-negative
        ENDDO

        ! Optional accuracy test: compute max|U·Σ·V^H - A|
        IF (NCON ==  1) THEN
           Write(6,*) JobU, JobVT
           Xmax = 0.d0
           ! For each element of A, reconstruct via SVD and compute residual
           DO I = 1,M
              DO I1 = 1,N
                 Z = cmplx(0.d0,0.d0,Kind(0.D0))  ! Reconstruct A(I,I1) = sum_J(U(I,J)·D(J)·V(J,I1))
                 DO J = 1,N
                    Z  =  Z + U(I,J) *D(J) *V(J,I1)  ! Sum over singular values
                 ENDDO
                 X = ABS(Z - A(I,I1))  ! Residual at (I,I1)
                 IF (X > Xmax ) Xmax = X  ! Track maximum error
              ENDDO
           ENDDO
           WRITE(6,*) "Success (0), PRE ", INFO, Xmax
        ENDIF

        ! Cleanup workspace
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

        ! Local workspace for symmetric eigenvalue decomposition
        INTEGER ND1, ND2, IERR, DN
        REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: WORK

        ! Check square matrix
        ND1 = SIZE(A,1)
        ND2 = SIZE(A,2)

        IF (ND1.NE.ND2) THEN
           WRITE(error_unit,*) 'DIAG_R: Error in matrix dimension'
           Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        ENDIF

        ! Initialize output: copy input to U (DSYEV modifies in-place)
        IERR = 0
        U = A       ! DSYEV overwrites U with eigenvectors
        W = 0       ! Initialize eigenvalues to zero
        
        ! Allocate sufficient workspace (conservative: 3*N is safe)
        DN = 3*ND1
        ALLOCATE(WORK(DN))
        
        ! Solve via LAPACK DSYEV: diagonalize symmetric real matrix
        ! 'V'=compute eigenvectors, 'U'=use upper triangle
        CALL DSYEV('V', 'U', ND1, U, ND1, W, WORK, DN, IERR)
        DEALLOCATE(WORK)

      END SUBROUTINE DIAG_R
!*********

!--------------------------------------------------------------------
!> @brief Hermitian eigenvalue decomposition for complex matrices.
!> @details Solves A·U = U·W for a complex Hermitian matrix via
!> LAPACK ZHEEV. Returns eigenvalues W (real) and eigenvectors U.
!> Optional validation test checks orthonormality and eigenvalue eq.
!>
!> @param[in] A (N×N) complex Hermitian input matrix
!> @param[out] U (N×N) complex eigenvector matrix (columns are eigenvectors)
!> @param[out] W (N) real eigenvalues (sorted ascending)
!> @note LAPACK ZHEEV solves Hermitian eigenvalue problem via QR
!> @see ZHEEV documentation
!--------------------------------------------------------------------
      SUBROUTINE DIAG_I(A,U,W)
        IMPLICIT NONE
        COMPLEX (Kind=Kind(0.d0)), INTENT(IN)   , DIMENSION(:,:) :: A
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:,:) :: U
        REAL    (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:)   :: W

        ! Variables for Hermitian eigenvalue decomposition via ZHEEV
        CHARACTER (len=1) :: UPLO, JOBZ
        INTEGER :: N, LWORK, INFO
        COMPLEX (Kind=Kind(0.d0)), allocatable :: WORK (:)
        REAL    (Kind=Kind(0.d0)), allocatable :: RWORK(:)
        Logical :: Test
        Integer :: I,J,m
        Complex (Kind=Kind(0.d0)) :: Z
        Real (Kind=Kind(0.d0)) :: X, XMAX

        ! Set options: JOBZ='V' (compute eigenvectors), UPLO='U' (use upper triangle)
        JOBZ = "V"
        UPLO = "U"
        N = size(A,1)  ! Matrix dimension
        
        ! Initialize: copy input (ZHEEV modifies in-place)
        U = A           ! U will be overwritten with eigenvectors
        
        ! Allocate workspace (conservative estimates)
        LWORK = 2*N - 1
        Allocate ( WORK(LWORK) )      ! Complex workspace
        Allocate ( RWORK(3*N-2) )     ! Real workspace

        ! Solve via LAPACK ZHEEV: diagonalize Hermitian complex matrix
        Call ZHEEV (JOBZ, UPLO, N, U, N, W, WORK, LWORK, RWORK, INFO)

        Deallocate (WORK, RWORK)

        ! Optional validation test (disabled by default, set Test=.true. to enable)
        Test = .false.
        If (Test) then
           ! Test reconstruction: A = U·diag(W)·U^H where U^H = conjugate transpose
           XMAX = 0.d0
           DO I = 1,N
              DO J = 1,N
                 ! Reconstruct A(I,J) = sum_m(U(I,m)·W(m)·U(J,m)^*)
                 Z = cmplx(0.d0,0.d0,Kind=Kind(0.d0))
                 DO m = 1,N
                    Z =  Z + U(I,m)*cmplx(W(m),0.d0, Kind=Kind(0.d0))*Conjg(U(J,m))
                 ENDDO
                 ! Compute residual: reconstruction error
                 Z = Z - A(I,J)
                 X = ABS(Z)
                 If (X > XMAX ) XMAX = X
              ENDDO
           ENDDO
           write(6,*) ' Test Diag_I: ', XMAX
        endif

      End SUBROUTINE DIAG_I

!--------------------------------------------------------------------
!> @brief Get current time in seconds since midnight.
!> @details Extracts hour, minute, second from system clock via date_and_time()
!> intrinsic and returns total seconds elapsed since 00:00:00 today.
!> Used for performance timing and profiling.
!>
!> @param[out] X seconds since midnight (real, includes day-to-day discontinuity)
!> @note Discontinuous at midnight; use time differences for accurate elapsed time
!--------------------------------------------------------------------
      SUBROUTINE SECONDS(X)
        IMPLICIT NONE
        REAL (Kind=Kind(0.d0)), INTENT(OUT) :: X

        ! Local workspace: date_and_time() returns 8-element integer array
        ! values(1)=year, values(2)=month, values(3)=day, values(4)=UTC offset (min)
        ! values(5)=hour, values(6)=minute, values(7)=second, values(8)=millisecond
        integer,dimension(8) :: V
        
        ! Query system clock via Fortran intrinsic
        call date_and_time(values=V)

        ! Convert HMS (hour, minute, second) to total seconds since midnight
        ! X = hours×3600 + minutes×60 + seconds (ignore milliseconds for Fortran 2008)
        X = DBLE(V(5)*3600 + V(6)*60 + V(7))

      END SUBROUTINE SECONDS

!====================================================
!> @brief Compute determinant of complex matrix via LU decomposition (function).
!> @details Uses LAPACK ZGETRF to compute LU factorization, then product of
!> diagonal elements (with sign correction from pivots). Modifies input matrix.
!>
!> @param[inout] Mat (N×N) complex input matrix (LU-decomposed in-place)
!> @param[in] N matrix dimension
!> @return complex determinant: product(diagonal) × (-1)^(# row swaps)
!> @note Input matrix is destroyed (replaced by LU decomposition)
!> @see ZGETRF documentation
!====================================================
      Complex (Kind=Kind(0.d0)) Function DET_C(Mat,N)

        Implicit none

        ! Arguments
        Integer, intent(in) :: N
        Complex(Kind=Kind(0.d0)), intent(inout) :: mat(:,:)

        integer :: i, info, LDmat
        integer, allocatable :: ipiv(:)
        integer :: sgn

        ! Allocate pivot array for row swaps
        allocate(ipiv(N))
        ipiv = 0

        ! Step 1: LU decomposition via LAPACK ZGETRF (in-place)
        LDmat=size(mat,1)
        call zgetrf(N, N, mat, LDmat, ipiv, info)

        ! Step 2: Determinant = product of diagonal elements from LU
        det_C = cmplx(1.d0, 0.d0, kind(0.d0) )
        do i = 1, N
           det_C = det_C*mat(i, i)  ! Multiply diagonal elements
        enddo

        ! Step 3: Apply sign correction from row swaps
        ! Each swap (IPVT(i) ≠ i) flips the sign
        sgn =  1
        do i = 1, N
           if(ipiv(i) /= i)  sgn = -sgn
        enddo
        if (sgn == -1 ) det_C = - det_C

        deallocate(ipiv)
      end function DET_C

!--------------------------------------------------------------------
!> @author F.Assaad
!> @brief Return LU diagonal factors from complex matrix LU decomposition.
!> @details Computes LU factorization via ZGETRF and returns diagonal
!> elements D(i) = L(i,i)·U(i,i) from factorization. Product ∏D(i) = det(A)×(-1)^(swaps).
!> Applies sign correction from row pivots to D(1) only.
!>
!> @param[in] Mat1 (N×N) complex input matrix (unchanged)
!> @param[out] D (N) diagonal elements from LU decomposition
!> @param[in] N matrix dimension
!> @note Input matrix Mat1 is not modified; working copy Mat is used
!> @see ZGETRF documentation; determinant = D(1)·D(2)·...·D(N)
!====================================================
      Subroutine DET_C_LU(Mat1,D,N)

        Implicit none

        ! Arguments
        Integer, intent(in) :: N
        Complex(Kind=Kind(0.d0)), intent(in)  :: mat1(N,N)
        Complex(Kind=Kind(0.d0)), intent(out) :: D(N)

        ! Local working copy and pivot array
        Complex(Kind=Kind(0.d0)) :: mat(N,N)
        integer :: i, info
        integer :: ipiv(N)
        integer :: sgn

        ! Copy input (ZGETRF modifies in-place)
        mat = mat1
        ipiv = 0

        ! Step 1: LU decomposition via LAPACK ZGETRF
        call zgetrf(N, N, mat, N, ipiv, info)

        ! Step 2: Extract diagonal elements from LU matrix
        do i = 1,N
           D(i) = mat(i,i)  ! Diagonal of LU = L(i,i)·U(i,i)
        enddo
        
        ! Step 3: Apply sign correction from row swaps to first element only
        sgn =  1
        do i = 1, N
           if(ipiv(i) /= i)  sgn = -sgn  ! Each swap flips sign
        enddo
        if (sgn == -1 ) D(1) = - D(1)  ! Apply accumulated sign to D(1)

      end Subroutine DET_C_LU

!--------------------------------------------------------------------
!> @author Florian Goth
!> @brief Conjugate transpose A^H of a complex matrix.
!> @details Computes A^H = (A*)^T via Fortran intrinsics conjg() and transpose().
!> Equivalent to adjoint operator in quantum mechanics. For Hermitian matrix,
!> A^H = A.
!>
!> @param[in] A (M×N) complex input matrix
!> @return (N×M) complex conjugate transpose (adjoint) of A
!>
!> @note Implementation using chaining conjg(transpose(a)) is compiler-optimized
!> (esp. by gcc) to generate tight inner loops, faster than explicit loops.
!> Fortran's transpose() and conjg() operations compose efficiently.
!--------------------------------------------------------------------
    function ct(a) result(b)
        complex(kind=kind(0.D0)), dimension(:,:), intent(in) :: a
        complex(kind=kind(0.D0)), dimension(size(a,1),size(a,1)) :: b
        ! Conjugate transpose: (A*)^T = conjg(A)^T
        b = conjg(transpose(a))
    end function ct

    END MODULE MyMats
