! compile with
! gfortran -Wall -std=f2003 -I ../../../Prog/  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.F90 ../../../Prog_8/cgr1.o ../../../Prog/UDV_WRAP.o ../../../Libraries/Modules/modules_90.a -llapack -lblas

Program TESTCGRP
USE UDV_State_mod
use cgr1_mod
implicit none
interface
      SUBROUTINE CGRPold(PHASE, GRUP, udvr, udvl)
        USE UDV_State_mod
        CLASS(UDV_State), intent(in) :: udvr, udvl
        COMPLEX(Kind=Kind(0.D0)), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX(Kind=Kind(0.D0)), Intent(INOUT) :: PHASE
        end subroutine CGRPold
end interface

        TYPE(UDV_State) :: udvr, udvl
        COMPLEX(Kind=Kind(0.D0)), Dimension(:,:), allocatable :: GRUPnew, GRUPold
        COMPLEX(Kind=Kind(0.D0)), Dimension(:), allocatable :: DL, DR
        COMPLEX(Kind=Kind(0.D0)) :: PHASEnew, Phaseold, Zre, Zim
        INTEGER         :: i, j, N_size, N_part
        
        N_size = 5
        N_part = N_size -2
        CALL udvl%alloc(N_size, N_part)
        CALL udvr%alloc(N_size, N_part)
        allocate(GRUPnew(N_size, N_size), GRUPold(N_size, N_size), DL(n_size), DR(n_size))
        ! set up test data
        Phasenew = 1.D0
        Phaseold = 1.D0
        GRUPnew = 0.D0
        GRUPold = 0.D0
        udvr%side = 'R'
        udvl%side = 'L'
        do i = 1, N_size
        do j = 1, N_part
        udvr%U(i, j) = 0!i + j
        udvl%U(i, j) = 0!i + j
        enddo
        if (i <= N_part) then
            DL(i) = 1
            call udvL%setscale(DL(i),i)
            DR(i) = i
            call udvr%setscale(Dr(i),i)
            udvl%U(i,i) = exp(cmplx(0,i, kind=kind(0.D0)))
            udvr%U(i,i) = exp(cmplx(0,0.3*i, kind=kind(0.D0)))
        endif
        enddo
call CGRP(PHASEnew, GRUPnew, udvr, udvl)
! run old code
! write (*,*) GRUPnew

call CGRPold(PHASEold, GRUPold, udvr, udvl)
! write (*,*) GRUPold

! compare GRUP results
! write (*,*) GRUPOLD
! write (*,*) "---------"
! write (*,*) GRUPNEW
    do i = 1, N_size
    do j = 1, N_size
    Zre = DBLE(GRUPnew(i,j)-GRUPold(i,j))
    Zim = aimag(GRUPnew(i,j)-GRUPold(i,j))
    if(Abs(Zre) > 1E-14) then
    if (Abs(Zre) > MAX(ABS(real(GRUPnew(i,j))), ABS(real(GRUPold(i,j))) )*1D-10) then
    write (*,*) "ERROR in real part", i, j, real(GRUPnew(i,j)), real(GRUPold(i,j)), "diff: ", &
    & Abs(Zre), "prec: ", MAX(ABS(real(GRUPnew(i,j))), ABS(real(GRUPold(i,j))) )*1D-7
    STOP 2
    endif
    endif
    if(Abs(Zim) > 1E-14) then
    if (Abs(Zim) > MAX(ABS(aimag(GRUPnew(i,j))), ABS(aimag(GRUPold(i,j))) )*1D-10) then
    write (*,*) "ERROR in imag part", aimag(GRUPnew(i,j)), aimag(GRUPold(i,j))
    STOP 3
    endif
    endif
    enddo
    enddo

    Zre = real(Phasenew-Phaseold)
    Zim = aimag(Phasenew-Phaseold)
    if (Abs(Zre) > MAX(ABS(real(Phasenew)), ABS(real(Phaseold)) )*1D-13) then
    write (*,*) "ERROR in real part", Phasenew, Phaseold, Abs(Phaseold)
    STOP 6
    endif
    if (Abs(Zim) > MAX(ABS(aimag(Phasenew)), ABS(aimag(Phaseold)) )*1D-13) then
    write (*,*) "ERROR in imag part", Phasenew, aimag(Phasenew), aimag(Phaseold)
    STOP 7
    endif
    CALL udvr%dealloc
    CALL udvl%dealloc
deallocate(GRUPnew, GRUPold, Dr, DL)
write (*,*) "success"
end Program TESTCGRP

      SUBROUTINE CGRPOld(PHASE, GRUP, udvr, udvl)
        Use UDV_State_mod
        use MyMats
        CLASS(UDV_State), INTENT(IN) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX (Kind=Kind(0.d0)), Intent(INOUT) :: PHASE
        
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:,:) :: sMat, rMat
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: work 
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
        Allocate(sMat(N_part,N_part),rMat(Ndim,N_part), ipiv(N_part), work(N_part))
        
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
        ! ZGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.do 10,i=1,n
        call ZGETRI(N_part, sMat, N_part, ipiv, work, N_part, info)
        call ZGEMM('N','N',Ndim,N_part,N_part,alpha,udvr%U(1,1),Ndim,sMat(1,1),N_part,beta,rMat(1,1),Ndim)
!         call initd(Grup,alpha)
        alpha=-1.d0
        call ZGEMM('N','C',Ndim,Ndim,N_part,alpha,rMat(1,1),Ndim,udvl%U(1,1),Ndim,beta,GRUP(1,1),Ndim)
        do n=1,Ndim
          Grup(n,n)=Grup(n,n)+cmplx(1.d0, 0.d0, kind(0.d0))
        enddo
        Deallocate(sMat,rMat, ipiv, work)
      
      END SUBROUTINE CGRPOld
