
!  Copyright (C) 2016 - 2018 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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
 
     SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief
!> Imaginary time propagation and udv decompostion from Ntau to Ntau1, Ntau1 > Ntau
!>
!> @details
!> On input   B(NTAU ,1)  =  UR*DR*VR \n
!> On output  B(NTAU1,1)  =  UR*DR*VR
!> @param[in] Ntau, Ntau1  Integer
!> @param[inout] UDVR Class(UDV_state)

!-------------------------------------------------------------------

        Use Hop_mod
        Use UDV_State_mod
#if defined(STAB2) || defined(STAB1)         
        Use Hamiltonian_main
        Use UDV_Wrap_mod
        Implicit None

        ! Arguments
        CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: udvr
        Integer,          Intent(IN) :: NTAU1, NTAU


        ! Working space.
        Complex (Kind=Kind(0.d0)) :: Z_ONE
        COMPLEX (Kind=Kind(0.d0)) :: V1(Ndim,Ndim), TMP(Ndim,Ndim), TMP1(Ndim,Ndim)
        Integer ::NCON, NT, I, J, n, nf
        Real (Kind=Kind(0.d0)) :: X

        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        
        Do nf = 1,N_FL
           CALL INITD(TMP,Z_ONE)
           DO NT = NTAU + 1, NTAU1
              !CALL MMULT(TMP1,Exp_T(:,:,nf) ,TMP)
              Call Hop_mod_mmthr(TMP,nf)
!             TMP = TMP1
              Do n = 1,Size(Op_V,1)
!                  X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(Tmp,Op_V(n,nf),nsigma%f(n,nt),'n')
              ENDDO
           ENDDO
           CALL MMULT(TMP1,TMP, udvr(nf)%U)
           if(allocated(udvr(nf)%V)) then
              DO J = 1,NDim
                  DO I = 1,NDim
                    TMP1(I,J) = TMP1(I,J)*udvr(nf)%D(J)
                  ENDDO
              ENDDO
              TMP = udvr(nf)%V
           endif
           CALL UDV_WRAP_Pivot(TMP1(:,1:UDVR(nf)%N_part), udvr(nf)%U, udvr(nf)%D, V1,NCON,Ndim,UDVR(nf)%N_part)
           if(allocated(udvr(nf)%V)) CALL MMULT(udvr(nf)%V, V1, TMP)
        ENDDO
#else
        Implicit None

        ! Arguments
        CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: udvr
        Integer,          Intent(IN) :: NTAU1, NTAU


        ! Working space.
        Integer :: NT, n, nf

        Do nf = 1,N_FL
           DO NT = NTAU + 1, NTAU1
              Call Hop_mod_mmthR(UDVR(nf)%U,nf)
              Do n = 1,Size(Op_V,1)
                 Call Op_mmultR(UDVR(nf)%U,Op_V(n,nf),nsigma%f(n,nt),'n')
              ENDDO
           ENDDO

           CALL UDVR(nf)%decompose
        ENDDO

#endif
        
      END SUBROUTINE WRAPUR
