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

#if defined(STAB2) ||  defined(STAB1) 
        !NOTE:    NTAU1 > NTAU.
        Use Operator_mod, only : Phi
        Use Hamiltonian
        Use Hop_mod
        Use UDV_Wrap_mod

        Implicit none

        ! Arguments
        COMPLEX (Kind=Kind(0.d0)) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
        COMPLEX (Kind=Kind(0.d0)) :: DL(Ndim,N_FL)
        Integer :: NTAU1, NTAU


        ! Working space.
        COMPLEX (Kind=Kind(0.d0)) ::  U1(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim), TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) ::  D1(Ndim), Z_ONE, beta
        Integer :: NT, NCON, n, nf
        Real    (Kind=Kind(0.d0)) ::  X
 


        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        Do nf = 1, N_FL
           CALL INITD(TMP,Z_ONE)
           DO NT = NTAU1, NTAU+1 , -1
              Do n = Size(Op_V,1),1,-1
                 X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultL(Tmp,Op_V(n,nf),X,Ndim)
              enddo
              !CALL MMULT( TMP1,Tmp,Exp_T(:,:,nf) )
              Call  Hop_mod_mmthl (Tmp, Tmp1,nf)
              Tmp = Tmp1
           ENDDO
           
           !Carry out U,D,V decomposition.
           CALL ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, UL(1, 1, nf), Ndim, beta, TMP1, Ndim)
           DO n = 1,NDim
              TMP1(:, n) = TMP1(:, n) * DL(n, nf)
           ENDDO
           CALL UDV_WRAP_Pivot(TMP1,U1,D1,V1,NCON,Ndim,Ndim)
           !CALL UDV(TMP,U1,D1,V1,NCON)
           UL(:, :, nf) = CONJG(TRANSPOSE(U1))
           CALL ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_ONE, VL(1, 1, nf), Ndim, V1, Ndim, beta, TMP1, Ndim)
           VL(:, :, nf) = TMP1
           DL(:, nf) = D1
        ENDDO

#else
        !NOTE:    NTAU1 > NTAU.
        Use Operator_mod, only : Phi
        Use Hop_mod
        Use Wrap_helpers
        Implicit none

        ! Arguments
        COMPLEX (Kind=Kind(0.d0)), intent(inout) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
        COMPLEX (Kind=Kind(0.d0)), intent(inout) :: DL(Ndim,N_FL)
        Integer, intent(in) :: NTAU1, NTAU


        ! Working space.
        COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:, :) :: TMP, TMP1
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE
        Integer :: NT, NCON, n, nf
        Real    (Kind=Kind(0.d0)) ::  X
 
        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.
        Allocate (TMP(Ndim,Ndim), TMP1(Ndim,Ndim))
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
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
           CALL ul_update_matrices(UL(:,:,nf), DL(:, nf), VL(:,:,nf), TMP, TMP1, Ndim, NCON)
           UL(:, :, nf) = CONJG(TRANSPOSE(UL(:,:,nf)))
        ENDDO
        deallocate(TMP, TMP1)
#endif
      END SUBROUTINE WRAPUL
      
