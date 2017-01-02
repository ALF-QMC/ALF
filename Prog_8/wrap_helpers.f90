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
!> This function updates the UDV matrices with the new matrix stored in TMP.
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
write (*,*) real(M(i, :))
enddo
write (*,*) "----------------------------"
end subroutine

 SUBROUTINE ul_update_matrices(U, D, V, V1, TMP, TMP1, Ndim, NCON)
        Use UDV_Wrap_mod
        Implicit None
        INTEGER, intent(in) :: Ndim, NCON
        COMPLEX (Kind=Kind(0.d0)) :: U(Ndim,Ndim), V(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim),TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: TMP2(Ndim, Ndim), TMP3(Ndim, Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim), TAU(Ndim), WORK(2*Ndim), RWORK(2*Ndim)
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        INTEGER :: n, IPVT(Ndim), INFO
        LOGICAL :: FORWRD

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        ! TMP1 = TMP^dagger * U^dagger
        CALL ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, U(1, 1), Ndim, beta, TMP1, Ndim)
        ! TMP1 = TMP1 * D
        DO n = 1,NDim
            TMP1(:, n) = TMP1(:, n) * D(n)
            D(n) = 1.D0
            IPVT(n) = n
        ENDDO
        write (*,*) "INPUT"
        TMP3 =  MATMUL(TMP1, V)
        CALL pm(TMP3, Ndim)
        ! QR Zerlegung von TMP1
!        CALL ZGEQRF(Ndim, Ndim, TMP1, Ndim, TAU, WORK, 2*Ndim, INFO)
        call ZGEQP3(Ndim, Ndim, TMP1, Ndim, IPVT, TAU, WORK, 2*Ndim, RWORK, INFO)
        ! Permute V
!        FORWRD = .true.
!        CALL ZLAPMT(FORWRD, Ndim, Ndim, V, Ndim, IPVT)

        !CALL UDV_WRAP_Pivot(TMP1, TMP, D, V1,NCON,Ndim,Ndim)
        ! V = V * V1^dagger
!        CALL ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_ONE, V(1, 1), Ndim, V1, Ndim, beta, TMP1, Ndim)
!        V = TMP1
write (*,*) "TMP1"
call pm(tmp1, Ndim)
        CALL ZTRMM('R', 'U', 'N', 'N', Ndim, Ndim, Z_ONE, TMP1, Ndim, V, Ndim)
        write (*,*) "New V"
        call pm(V, Ndim)
        CALL ZUNGQR(Ndim, Ndim, Ndim, TMP1, Ndim, TAU, WORK, 2*Ndim, INFO)
        U = TMP1
        WRITE (*,*) "OUTPUT"
        TMP3 = MATMUL(V, U)
        CALL pm(TMP3, Ndim)
        STOP 2
END SUBROUTINE ul_update_matrices

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function updates the UDV matrices with the new matrix stored in TMP.
!
!> @param [inout] U
!> @param [inout] D
!> @param [inout] V
!> @param [in] TMP
!> @param [in] TMP1
!> @param [in] Ndim The size of the matrices
!> @param [in] NCON wether we check.
!-------------------------------------------------------------------
 SUBROUTINE ur_update_matrices(U, D, V, V1, TMP, TMP1, Ndim, NCON)
        Use UDV_Wrap_mod
        Use MyMats
        Implicit None
        INTEGER, intent(in) :: Ndim, NCON
        COMPLEX (Kind=Kind(0.d0)) :: U(Ndim,Ndim), V(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim),TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: TMP2(Ndim, Ndim), TMP3(Ndim, Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim), TAU(Ndim), WORK(2*Ndim), RWORK(2*Ndim)
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta, tempC
        INTEGER :: n, IPVT(Ndim), IPVTINV(Ndim), INFO, IPVT2(Ndim), p(Ndim), tempidx, k
        INTEGER :: i,j
        INTEGER :: visited(Ndim), cycle(Ndim), cycstart, cycend, idx, visitedindices, cyclelen
        LOGICAL :: FORWRD
        REAL(Kind=Kind(0.d0)) :: XMAX, XMEAN
        
        ! QR(TMP * U * D) * V
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        CALL ZGEMM('N', 'N', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, U(1, 1), Ndim, beta, TMP1, Ndim)
        ! TMP1 = TMP1 * D
        DO n = 1,NDim
            TMP1(:, n) = TMP1(:, n)*D(n)
            D(n) = 1.D0
            IPVT(n) = 0
        ENDDO
!                 write (*,*) "COMPARISON"
!                 TMP = MATMUL(TMP1, V)
!                 call pm(TMP, Ndim)
!                write (*,*) "input to QR"
!                call pm(TMP1, Ndim)
        ! QR Zerlegung von TMP1
!        CALL ZGEQRF(Ndim, Ndim, TMP1, Ndim, TAU, WORK, 2*Ndim, INFO)
        call ZGEQP3(Ndim, Ndim, TMP1, Ndim, IPVT, TAU, WORK, 2*Ndim, RWORK, INFO)
        FORWRD = .true.
        CALL ZLAPMR(FORWRD, Ndim, Ndim, V, Ndim, IPVT) ! lapack 3.4.2
        CALL ZTRMM('L', 'U', 'N', 'N', Ndim, Ndim, Z_ONE, TMP1, Ndim, V, Ndim)
        CALL ZUNGQR(Ndim, Ndim, Ndim, TMP1, Ndim, TAU, WORK, 2*Ndim, INFO)
        U = TMP1
!        IF (NCON == 1) THEN
!           !Check the result  A = U D V
!           DO J = 1,N2
!              TMP3(:, J) = D * V(:, J)
!           ENDDO
!           !Write(6,*) 'Here'
!           Call MMULT (TMP3, U, TMP2)
!           !Call MMULT (A2,U,V)
!           CALL COMPARE(A, TMP3, XMAX,XMEAN)
!           Write (6,*) 'Check afer NEW Pivoting', XMAX
!        ENDIF
END SUBROUTINE ur_update_matrices
