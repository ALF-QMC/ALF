
!  Copyright (C) 2016 - 2022 The ALF project
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

module wrapur_mod
   implicit none
   contains
 
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
        Integer ::NCON, NT, I, J, n, nf, nf_eff
        Real (Kind=Kind(0.d0)) :: X

        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        
        Do nf_eff = 1,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           CALL INITD(TMP,Z_ONE)
           DO NT = NTAU + 1, NTAU1
              !CALL MMULT(TMP1,Exp_T(:,:,nf) ,TMP)
              Call Hop_mod_mmthr(TMP,nf,nt)
!             TMP = TMP1
              Do n = 1,Size(Op_V,1)
!                  X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(Tmp,Op_V(n,nf),nsigma%f(n,nt),'n',nt)
              ENDDO
           ENDDO
           CALL MMULT(TMP1,TMP, udvr(nf_eff)%U)
           if(allocated(udvr(nf_eff)%V)) then
              DO J = 1,NDim
                  DO I = 1,NDim
                    TMP1(I,J) = TMP1(I,J)*udvr(nf_eff)%D(J)
                  ENDDO
              ENDDO
              TMP = udvr(nf_eff)%V
           endif
           CALL UDV_WRAP_Pivot(TMP1(:,1:UDVR(nf_eff)%N_part), udvr(nf_eff)%U, udvr(nf_eff)%D, V1,NCON,Ndim,UDVR(nf_eff)%N_part)
           if(allocated(udvr(nf_eff)%V)) CALL MMULT(udvr(nf_eff)%V, V1, TMP)
           ! Test if scales in D are appraoching the limit of double precision.
           CALL UDVR(nf_eff)%testscale
        ENDDO
#else
        Implicit None

        ! Arguments
        CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: udvr
        Integer,          Intent(IN) :: NTAU1, NTAU


        ! Working space.
        Integer :: NT, n, nf, nf_eff

        Do nf_eff = 1,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           DO NT = NTAU + 1, NTAU1
              Call Hop_mod_mmthR(UDVR(nf_eff)%U,nf,nt)
              Do n = 1,Size(Op_V,1)
                 Call Op_mmultR(UDVR(nf_eff)%U,Op_V(n,nf),nsigma%f(n,nt),'n',nt)
              ENDDO
           ENDDO

           CALL UDVR(nf_eff)%decompose
           ! Test if scales in D are appraoching the limit of double precision.
           CALL UDVR(nf_eff)%testscale
        ENDDO

#endif
        
      END SUBROUTINE WRAPUR

end module wrapur_mod
