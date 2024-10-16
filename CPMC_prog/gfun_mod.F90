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
