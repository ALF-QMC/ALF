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

module gfun_mod
  implicit none
  contains
      
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Computes the Green's function in the projective implementation.
!
!> @param[out] PHASE
!> @param[out] GRUP
!> @param[in] udvr
!> @param[in] udvl
!
!--------------------------------------------------------------------
      SUBROUTINE CGRP(zdet, GRUP, udvr, udvl)
        Use UDV_State_mod
        CLASS(UDV_State), INTENT(IN) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Intent(OUT) :: GRUP
        Complex (Kind=Kind(0.d0)), Intent(OUT)  :: zdet

        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:,:) :: sMat, rMat
        INTEGER, allocatable :: ipiv(:)
        COMPLEX (Kind=Kind(0.d0)) :: alpha, beta, phase
        INTEGER :: Ndim, N_part, info, n

        if((udvl%side .ne. "L") .and. (udvl%side .ne. "l") ) then
          write(*,*) "cgrp: udvl is not of type left"
          write(*,*) "cgrp: actual side is ",udvl%side
        endif
        if((udvr%side .ne. "R") .and. (udvr%side .ne. "r") ) then
          write(*,*) "cgrp: udvr is not of type right"
          write(*,*) "cgrp: actual side is ",udvr%side
        endif

        Ndim = udvl%ndim
        N_part = udvl%n_part
        Allocate(sMat(N_part,N_part), ipiv(N_part), rMat(N_part, Ndim))
        
        ! Gr = Ur (Ul Ur)^-1 Ul
        ! Phase = 1 + Ur (Ul Ur)^-1 Ul
        ! Ul = udvl%U ^dag
        alpha=1.d0
        beta=0.d0
        call ZGEMM('C','N',N_part,N_part,Ndim,alpha,udvl%U(1,1),Ndim,udvr%U(1,1),Ndim,beta,sMat(1,1),N_part)

        ! ZGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(N_part, N_part, sMat, N_part, ipiv, info)
        ! obtain log of det
        zdet  = 0.d0
        phase = 1.d0
        Do n=1,N_part
           if (ipiv(n).ne.n) then
              phase = -phase
           endif
           zdet = zdet + log(sMat(n,n))
        enddo
        zdet = zdet + log(phase)

        rMat = conjg(transpose(udvl%U))
        call zgetrs('N', N_part, Ndim, sMat(1,1), N_part, ipiv, rMat(1,1), N_part, info)
        alpha=-1.d0
        call ZGEMM('N','N',Ndim,Ndim,N_part,alpha,udvr%U(1,1),Ndim,rMat(1,1),N_part,beta,GRUP(1,1),Ndim)
        do n=1,Ndim
          GRUP(n,n)=GRUP(n,n)+cmplx(1.d0, 0.d0, kind(0.d0))
        enddo
        Deallocate(sMat, rMat, ipiv)
      
      END SUBROUTINE CGRP

end module gfun_mod
