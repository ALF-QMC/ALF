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
!> @param [inout] U
!> @param [inout] D
!> @param [inout] V
!> @param [in] TMP
!> @param [in] TMP1
!> @param [in] Ndim The size of the matrices
!> @param [in] NCON wether we check.
!-------------------------------------------------------------------
 SUBROUTINE ul_update_matrices(U, D, V, V1, TMP, TMP1, Ndim, NCON)
        Use UDV_Wrap_mod
        Implicit None
        INTEGER, intent(in) :: Ndim, NCON
        COMPLEX (Kind=Kind(0.d0)) :: U(Ndim,Ndim), V(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim),TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim)
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        INTEGER :: n

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        CALL ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, U(1, 1), Ndim, beta, TMP1, Ndim)
        DO n = 1,NDim
            TMP1(:, n) = TMP1(:, n) * D(n)
        ENDDO
        CALL UDV_WRAP_Pivot(TMP1, TMP, D, V1,NCON,Ndim,Ndim)
        CALL ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_ONE, V(1, 1), Ndim, V1, Ndim, beta, TMP1, Ndim)
        V = TMP1
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
        Implicit None
        INTEGER, intent(in) :: Ndim, NCON
        COMPLEX (Kind=Kind(0.d0)) :: U(Ndim,Ndim), V(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim),TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim)
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        INTEGER :: n

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        CALL ZGEMM('N', 'N', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, U(1, 1), Ndim, beta, TMP1, Ndim)
        DO n = 1,NDim
            TMP1(:, n) = TMP1(:, n)*D(n)
        ENDDO
        CALL UDV_WRAP_Pivot(TMP1, U, D , V1,NCON,Ndim,Ndim)
        CALL ZGEMM('N', 'N', Ndim, Ndim, Ndim, Z_ONE, V1, Ndim, V(1, 1), Ndim, beta, TMP1, Ndim)
        V = TMP1
END SUBROUTINE ur_update_matrices
