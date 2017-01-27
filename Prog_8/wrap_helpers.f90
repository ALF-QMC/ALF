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
!> This function updates the UDV matrices with the new matrix stored in TMP:
!
!> @param [in] U
!> @param [inout] D
!> @param [inout] V
!> @param [in] TMP
!> @param [in] TMP1
!> @param [in] Ndim The size of the matrices
!> @param [in] NCON wether we check.
!-------------------------------------------------------------------

SUBROUTINE pm(M, gr)
Implicit NONE
Integer :: gr
COMPLEX (Kind=Kind(0.d0)), intent(in) :: M(gr,gr)
Integer :: i,j
do i = 1, gr
write (*,*) M(i, :)
enddo
write (*,*) "----------------------------"
end subroutine

 SUBROUTINE ul_update_matrices(U, D, V, TMP, TMP1, Ndim, NCON)
        Use UDV_Wrap_mod
        Implicit None
        Integer :: IZAMAX, izmax1
        REAL(Kind=Kind(0.D0)) :: DZSUM1, DZNRM2
        INTEGER, intent(in) :: Ndim, NCON
        COMPLEX (Kind=Kind(0.d0)), intent(in) :: TMP(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)), intent(inout) :: U(Ndim,Ndim), V(Ndim,Ndim), TMP1(Ndim,Ndim), D(Ndim)
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK, RWORK
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        INTEGER, allocatable, Dimension(:) :: IPVT
        INTEGER :: INFO, i, j, t1, LWORK
        LOGICAL :: FORWRD
        REAL(Kind=Kind(0.D0)) :: X

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        ! TMP1 = TMP^dagger * U^dagger
        CALL ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, U(1, 1), Ndim, beta, TMP1, Ndim)
        ! TMP1 = TMP1 * D
        DO i = 1,NDim
            U(:, i) = TMP1(:, i) * D(i)
        ENDDO
        ALLOCATE(TAU(Ndim), RWORK(2*Ndim), IPVT(Ndim))
        IPVT = 0
        ! Query and allocate optimal amount of work space
        call ZGEQP3(Ndim, Ndim, U, Ndim, IPVT, TAU, BETA, -1, RWORK, INFO)
        LWORK = INT(DBLE(BETA))
        ALLOCATE(WORK(LWORK))
        ! QR decomposition of U with full column pivoting, A P = Q R)
        call ZGEQP3(Ndim, Ndim, U, Ndim, IPVT, TAU, WORK, LWORK, RWORK, INFO)
        ! separate off D and calculate the respective row-norms
        do i = 1, Ndim
        ! plain diagonal entry
             X = ABS(U(i, i))
!             ! a inf-norm
!             X = U(i, i+izamax(Ndim+1-i, U(i, i), Ndim)-1)
!             ! another inf-norm
!             X = U(i, i-1+izmax1(Ndim+1-i, U(i, i), Ndim))
!             ! 1-norm
!            X = DZSUM1(Ndim+1-i, U(i, i), Ndim)
            ! 2-norm
!            X = DZNRM2(Ndim+1-i, U(i, i), Ndim)
!             write (*, *) i, ABS(TMP1(i, i)), i+izamax(Ndim+1-i, TMP1(i, i), Ndim)-1, TMP1(i, i+izamax(Ndim+1-i, & 
!             & TMP1(i, i), Ndim)-1), i-1 + izmax1(Ndim+1-i, TMP1(i, i), Ndim), TMP1(i, i-1+izmax1(Ndim+1-i, &
!             & TMP1(i, i), Ndim)), DZSUM1(Ndim+1-i, TMP1(i, i), Ndim), DZNRM2(Ndim+1-i, TMP1(i, i), Ndim)
            D(i) = X
            do j = i, Ndim
                U(i, j) = U(i, j) / X
            enddo
        enddo
        ! Permute V, since we multiply with V from the left we have to permute its columns
        FORWRD = .true.
        CALL ZLAPMT(FORWRD, Ndim, Ndim, V, Ndim, IPVT)
        ! V = V * R^dagger
        CALL ZTRMM('R', 'U', 'C', 'N', Ndim, Ndim, Z_ONE, U, Ndim, V, Ndim)
        ! create explicitly U in the storage already present for it
        CALL ZUNGQR(Ndim, Ndim, Ndim, U, Ndim, TAU, WORK, LWORK, INFO)
        DEALLOCATE(TAU, WORK, RWORK, IPVT)
END SUBROUTINE ul_update_matrices

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function updates the UDV matrices with the new matrix stored in TMP.
!> Essentially we calculate the product TMP * U * D * V
!> For the result we generate a new decomposition in the form U, D, V
!
!> @param [inout] U A unitary matrix in full storage.
!> @param [inout] D The entries of a diagonal matrix.
!> @param [inout] V A full matrix
!> @param [in] TMP A full matrix
!> @param [in] TMP1 temporary storage
!> @param [in] Ndim The size of the matrices
!> @param [in] NCON wether we check.(TODO: currently not used)
!-------------------------------------------------------------------
 SUBROUTINE ur_update_matrices(U, D, V, TMP, TMP1, Ndim, NCON)
        Use UDV_Wrap_mod

        Implicit None
        INTEGER, intent(in) :: Ndim, NCON
        Integer :: IZAMAX, izmax1
        REAL(Kind=Kind(0.D0)) :: DZSUM1, DZNRM2
        COMPLEX (Kind=Kind(0.d0)), intent(in) :: TMP(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)), intent(inout) :: U(Ndim,Ndim), V(Ndim,Ndim), TMP1(Ndim,Ndim), D(Ndim)
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK, RWORK
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        INTEGER :: INFO, i, j, LWORK
        INTEGER, allocatable, Dimension(:) :: IPVT
        LOGICAL :: FORWRD
        REAL(Kind=Kind(0.d0)) :: X!, XMAX, XMEAN
        
        ! QR(TMP * U * D) * V
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        CALL ZGEMM('N', 'N', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, U(1, 1), Ndim, beta, TMP1, Ndim)
        ! TMP1 = TMP1 * D
        DO i = 1,NDim
            U(:, i) = TMP1(:, i)*D(i)
        ENDDO
        ALLOCATE(TAU(Ndim), RWORK(2*Ndim), IPVT(Ndim))
        IPVT = 0
        ! Query and allocate optimal amount of work space
        call ZGEQP3(Ndim, Ndim, U, Ndim, IPVT, TAU, beta, -1, RWORK, INFO)
        LWORK = INT(DBLE(BETA))
        ALLOCATE(WORK(LWORK))
        ! QR decomposition of U with full column pivoting, AP = QR
        call ZGEQP3(Ndim, Ndim, U, Ndim, IPVT, TAU, WORK, LWORK, RWORK, INFO)
        ! separate off D
        do i = 1, Ndim
        ! plain diagonal entry
             X = ABS(U(i, i))
!             ! a inf-norm
!             X = U(i, i+izamax(Ndim+1-i, U(i, i), Ndim)-1)
!             ! another inf-norm
!             X = U(i, i-1+izmax1(Ndim+1-i, U(i, i), Ndim))
            ! 1-norm
!            X = DZSUM1(Ndim+1-i, U(i, i), Ndim)
            ! 2-norm
!            X = DZNRM2(Ndim+1-i, U(i, i), Ndim)
            D(i) = X
            DO j = i, Ndim
                U(i, j) = U(i, j) / X
            ENDDO
        enddo
        ! Permute V. Since we multiply with V from the right we have to permute the rows.
        ! A V = A P P^-1 V = Q R P^-1 V
        FORWRD = .true.
        CALL ZLAPMR(FORWRD, Ndim, Ndim, V, Ndim, IPVT) ! lapack 3.3
        ! V = R * V
        CALL ZTRMM('L', 'U', 'N', 'N', Ndim, Ndim, Z_ONE, U, Ndim, V, Ndim)
        ! Generate explicit U in the previously abused storage of U
        CALL ZUNGQR(Ndim, Ndim, Ndim, U, Ndim, TAU, WORK, LWORK, INFO)
        DEALLOCATE(TAU, WORK, RWORK, IPVT)
END SUBROUTINE ur_update_matrices
