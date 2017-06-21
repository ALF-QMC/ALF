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

  SUBROUTINE WRAPGRUP(GR,NTAU,PHASE)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given the green function matrix GR at time NTAU  the routine   propagates 
!> it to time  NTAU + 1 and carries  out an update of the fields at time NTAU+1
!> NTAU: [0:LTROT-1]
!
!--------------------------------------------------------------------

        Use Hamiltonian
        Use Hop_mod
        Use Wrapgr_mod
        Implicit none

        Interface
           Subroutine Upgrade(GR,N_op,NT,PHASE,Op_dim) 
             Use Hamiltonian 
             Implicit none 
             Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL) 
             Integer, INTENT(IN) :: N_op, Nt, Op_dim
             Complex (Kind=Kind(0.d0)) :: Phase
           End Subroutine Upgrade
        End Interface

        

        ! Arguments
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  GR(:,:,:)
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  PHASE
        INTEGER, INTENT(IN) :: NTAU
        
        !Local 
        Integer :: nf, N_Type, NTAU1,n, m
        Complex (Kind=Kind(0.d0)) :: Mat_TMP(Ndim,Ndim)
        Real    (Kind=Kind(0.d0)) :: X
        Character (Len=1)  :: Direction

        ! Wrap up, upgrade ntau1.  with B^{1}(tau1) 
        NTAU1 = NTAU + 1
        Do nf = 1,N_FL
           CALL HOP_MOD_mmthr( GR(:,:,nf),  MAT_TMP,nf)
           CALL HOP_MOD_mmthl_m1(MAT_TMP,GR(:,:,nf), nf )
           !CALL MMULT ( MAT_TMP,    Exp_T(:,:,nf), GR(:,:,nf)        )
           !CALL MMULT ( GR(:,:,nf), MAT_TMP      , Exp_T_M1(:,:,nf)  )
        Enddo
        Do n = Nt_sequential_start,Nt_sequential_end
           ! Write(6,*) 'Hi'
           Do nf = 1, N_FL
              X = Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
              N_type = 1
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
           enddo
           nf = 1
           Call Upgrade(GR,N,ntau1,PHASE,Op_V(n,nf)%N_non_Zero) 
           do nf = 1,N_FL
              N_type =  2
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
           enddo
        Enddo
        If ( N_Global_tau > 0 ) then 
           direction = "u"
           m         = Nt_sequential_end 
           !Write(6,*) "Calling  Global_tau_mod_Random_update:", m,direction, Size(OP_V,1)
           Call Wrapgr_Random_update(GR,m,ntau1,direction, PHASE )
        Endif
        
        


      END SUBROUTINE WRAPGRUP
