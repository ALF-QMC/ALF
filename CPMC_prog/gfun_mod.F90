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
      SUBROUTINE CGRP(phase, GRUP, udvr, udvl)
        Use UDV_State_mod
        CLASS(UDV_State), INTENT(IN) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Intent(OUT) :: GRUP
        COMPLEX (Kind=Kind(0.d0)), Intent(OUT) :: phase

        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:,:) :: sMat, rMat
        INTEGER, allocatable :: ipiv(:)
        COMPLEX (Kind=Kind(0.d0)) :: alpha, beta
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
        phase=1.d0
        Do n=1,N_part
          if (ipiv(n).ne.n) then
            phase = -phase * sMat(n,n)/abs(sMat(n,n))
          else
            phase =  phase * sMat(n,n)/abs(sMat(n,n))
          endif
        enddo
        rMat = conjg(transpose(udvl%U))
        call zgetrs('N', N_part, Ndim, sMat(1,1), N_part, ipiv, rMat(1,1), N_part, info)
        alpha=-1.d0
        call ZGEMM('N','N',Ndim,Ndim,N_part,alpha,udvr%U(1,1),Ndim,rMat(1,1),N_part,beta,GRUP(1,1),Ndim)
        do n=1,Ndim
          GRUP(n,n)=GRUP(n,n)+cmplx(1.d0, 0.d0, kind(0.d0))
        enddo
        Deallocate(sMat, rMat, ipiv)
      
      END SUBROUTINE CGRP

      Subroutine Compute_overlap(Det_Vec, udvr, udvl)

        Use UDV_State_mod
        Use Hamiltonian_main

        Implicit none

        Complex (Kind=Kind(0.d0)), Dimension(:), Intent(OUT)  ::  Det_Vec
        CLASS(UDV_State),  INTENT(INOUT) :: udvr(N_FL_eff)
        CLASS(UDV_State),  INTENT(IN   ) :: udvl(N_FL_eff)

        !> Local variables
        Integer ::  N_size, NCON, I,  J, N_part, info, NSTM, N, nf, nst, nt, nt1, nf_eff
        Integer, allocatable :: ipiv(:)
        COMPLEX (Kind=Kind(0.d0)) :: alpha,beta, Z, Z1
        TYPE(UDV_State) :: udvlocal
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Allocatable ::  TP!, U, V
        COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: D

        N_part=udvr(1)%N_part
        N_size=udvr(1)%ndim
        Det_vec = 0.d0
        Allocate (TP(N_part,N_part), ipiv(N_part))
        do nf_eff=1,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           alpha=1.d0
           beta=0.d0
           CALL ZGEMM('C','N',N_part,N_part,N_size,alpha,udvl(nf_eff)%U(1,1),N_size,udvr(nf_eff)%U(1,1),N_size,beta,TP(1,1),N_part)
           ! ZGETRF computes an LU factorization of a general M-by-N matrix A
           ! using partial pivoting with row interchanges.
           call ZGETRF(N_part, N_part, TP(1,1), N_part, ipiv, info)
           Z=1.d0
           Do J=1,N_part
              if (ipiv(J).ne.J) then
                 Z = -Z
              endif
              Z =  Z * TP(J,J)
           enddo
           Det_vec  (nf) = Det_vec(nf)+log(Z)
        enddo
        Deallocate(TP,ipiv)
        return

      end subroutine Compute_overlap

end module gfun_mod
