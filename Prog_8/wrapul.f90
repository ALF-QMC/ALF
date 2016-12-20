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

      SUBROUTINE WRAPUL(NTAU1, NTAU, UL ,DL, VL)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given    B(LTROT,NTAU1,Nf  ) =  VL, DL, UL
!> Returns  B(LTROT,NTAU, Nf  ) =  VL, DL, UL
!
!--------------------------------------------------------------------


        !NOTE:    NTAU1 > NTAU.
        Use Operator_mod, only : Phi
        Use Hop_mod
        Implicit none

        ! Arguments
        COMPLEX (Kind=Kind(0.d0)) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
        COMPLEX (Kind=Kind(0.d0)) :: DL(Ndim,N_FL)
        Integer :: NTAU1, NTAU


        ! Working space.
        COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:, :) ::  V1, TMP, TMP1
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        Integer :: NT, NCON, n, nf
        Real    (Kind=Kind(0.d0)) ::  X
 
        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.
        Allocate (V1(Ndim,Ndim), TMP(Ndim,Ndim), TMP1(Ndim,Ndim))
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        Do nf = 1, N_FL
           CALL INITD(TMP,Z_ONE)
           DO NT = NTAU1, NTAU+1 , -1
              Do n = Size(Op_V,1),1,-1
                 X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultL(TMP,Op_V(n,nf),X,Ndim)
              enddo
              !CALL MMULT( TMP1,Tmp,Exp_T(:,:,nf) )
              Call  Hop_mod_mmthl (TMP, TMP1,nf)
              TMP = TMP1
           ENDDO
           
           !Carry out U,D,V decomposition.
           CALL ul_update_matrices(UL(:,:,nf), DL(:, nf), VL(:,:,nf), V1, TMP, TMP1, Ndim, NCON)
           UL(:, :, nf) = CONJG(TRANSPOSE(TMP))
        ENDDO
        deallocate(V1, TMP, TMP1)
      END SUBROUTINE WRAPUL
      
