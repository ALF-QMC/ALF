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


      SUBROUTINE WRAPGRDO(GR,NTAU,PHASE)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given the green function matrix GR at time NTAU  the routine   carries out an 
!> update of the fields at time NTAU and  propagates  it to time NTAU-1
!> NTAU: [LTROT:1]
!
!--------------------------------------------------------------------    
    Use Hamiltonian
        Use MyMats
        Use Hop_mod
        Implicit None

        Interface
           Subroutine Upgrade(GR,N_op,NT,PHASE,Op_dim) 
             Use Hamiltonian 
             Implicit none 
             Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL) 
             Integer, INTENT(IN) :: N_op, Nt, Op_dim
             Complex (Kind=Kind(0.d0)) :: Phase
           End Subroutine Upgrade
        End Interface
        
        ! Given GREEN at time NTAU => GREEN at time NTAU - 1,
        ! Upgrade NTAU  [LTROT:1]
        
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: GR(:,:,:)
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
        Integer :: NTAU

        ! Local
        Complex (Kind=Kind(0.d0)) :: Mat_TMP(Ndim,Ndim)
        Integer :: nf, N_Type, n
        real (Kind=Kind(0.d0)) :: spin

        Do n =  size(OP_V,1), 1, -1 
           N_type = 2
           nf = 1
           spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
           do nf = 1,N_FL
              Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), spin, Ndim, N_Type)
           enddo
           !Write(6,*) 'Upgrade : ', ntau,n 
           Call Upgrade(GR,n,ntau,PHASE,Op_V(n,1)%N_non_zero) 
           ! The spin has changed after the upgrade!
           nf = 1
           spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
           N_type = 1
           do nf = 1,N_FL
              Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), spin, Ndim, N_Type )
           enddo
        enddo
        DO nf = 1,N_FL
           Call Hop_mod_mmthl   (GR(:,:,nf), MAT_TMP, nf)
           Call Hop_mod_mmthr_m1(MAT_TMP, GR(:,:,nf), nf)
           !CALL MMULT(MAT_TMP   , GR(:,:,nf)      , Exp_T(:,:,nf) )
           !CALL MMULT(GR(:,:,nf), Exp_T_M1(:,:,nf), MAT_TMP       )
        enddo

!!$        ! Test
!!$        Mat_TMP = cmplx(0.d0,0.d0)
!!$        DO I = 1,Ndim
!!$           Mat_TMP(I,I) = cmplx(1.d0,0.d0) 
!!$        Enddo
!!$        Do n = size(Op_V,1), 1, -1
!!$           N_type = 2
!!$           nf = 1
!!$           spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
!!$           Write(6,*) n, spin
!!$           do nf = 1,N_FL
!!$              Call Op_Wrapdo( Mat_tmp, Op_V(n,nf), spin, Ndim, N_Type)
!!$           enddo
!!$           !Upgrade
!!$           N_type = 1
!!$           do nf = 1,N_FL
!!$              Call Op_Wrapdo( Mat_tmp, Op_V(n,nf), spin, Ndim, N_Type )
!!$           enddo
!!$        enddo
!!$        
!!$        DO I = 1,Ndim
!!$           Do J = 1,NDIM
!!$              WRITE(6,*) I,J, Mat_tmp(I,J)
!!$           ENDDO
!!$        ENDDO
!!$
!!$        STOP

      END SUBROUTINE WRAPGRDO
