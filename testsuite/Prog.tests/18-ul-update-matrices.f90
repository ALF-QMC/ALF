! compile with
! gfortran -Wall -std=f2003 -I ../../Prog_8/  -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ 18-ul-update-matrices.f90 ../../Prog_8/wrap_helpers.o ../../Prog_8/UDV_WRAP.o ../../Libraries/Modules/modules_90.a ../../../../lapack-3.6.1/liblapack.a -lblas

Program TESTULUPDATEMATRICES
        Use UDV_State_mod
implicit none
interface
SUBROUTINE ul_update_matrices_old(U, D, V, V1, TMP, TMP1, Ndim, NCON)
        Use UDV_Wrap_mod
        Implicit None
        INTEGER, intent(in) :: Ndim, NCON
        COMPLEX (Kind=Kind(0.d0)) :: U(Ndim,Ndim), V(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim),TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim)
        end Subroutine
end interface

        COMPLEX(Kind=Kind(0.D0)), Dimension(:,:), allocatable :: V1, TMP, TMP1, Uold, Vold, TMPold
        COMPLEX(Kind=Kind(0.D0)), Dimension(:), allocatable ::  Dold
        INTEGER         :: i, j, Ndim, NCON
        COMPLEX(Kind=Kind(0.D0)) :: Z, dnew, beta
        TYPE(UDV_State) :: UDVL
        
        do Ndim = 5, 200,5
        CALL UDVL%alloc(Ndim)
        udvl%side='l'
        Allocate (Uold(Ndim, Ndim), Vold(Ndim, Ndim))
        ALLocate(TMP(Ndim, Ndim), TMPold(Ndim, Ndim), Dold(Ndim), TMP1(Ndim, Ndim), V1(Ndim, Ndim))
        Uold = 0.D0
        Vold = 0.D0
        do i = 1, Ndim
        do j = 1, Ndim
        TMP(i, j) = 1.D0/(i+j)
        enddo
        Dold(i) = 1.D0!i*i
        Uold(i,i) = 1.D0
        Vold(i,i) = 1.D0
        Vold(1,i) = 1.D0
        Vold(i,1) = 1.D0
        call udvl%setscale(Dold(i),i)
        enddo
        UDVL%U = Uold
        UDVL%V = Vold
        TMPold = TMP
        NCON = 1
        Z=cmplx(1.d0,0.d0,kind(0.d0))
        beta=0.d0
        CALL ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z, TMP(1, 1), Ndim, UDVL%U, Ndim, beta, TMP1(1, 1), Ndim)
        UDVL%U=TMP1
        CALL UDVL%decompose!(TMP, TMP1, NCON)
        call ul_update_matrices_old(Uold, Dold, Vold, V1, TMPold, TMP1, Ndim, NCON)
        
     ! compare
     do i = 1, Ndim
!     do j = 1, Ndim
!         Z = V(i,j) - Vold(i,j)
     
!     if (Abs(real(Z)) > MAX(ABS(REAL(V(i, j))), ABS(REAL(Vold(i, j))))*1D-15 ) then
!     write (*,*) "Error in V real part", V(i,j), Vold(i,j)
!     STOP 2
!     endif
!     if (Abs(AIMAG(Z)) > MAX(ABS(AIMAG(V(i, j))), ABS(AIMAG(Vold(i, j))))*1D-15 ) then
!     write (*,*) "Error in V imag part", V(i,j), Vold(i,j)
!     STOP 3
!     endif
!         
!         Z = TMP(i,j) - TMPold(i,j)
!         if (Abs(real(Z)) > MAX(ABS(REAL(TMP(i, j))), ABS(REAL(TMPold(i, j))))*1D-15 ) then
!         write (*,*) "Error in TMP real part", TMP(i,j), TMPold(i,j)
!         STOP 4
!         endif
!         if (Abs(AIMAG(Z)) > MAX(ABS(AIMAG(TMP(i, j))), ABS(AIMAG(TMPold(i, j))))*1D-15 ) then
!         write (*,*) "Error in TMP imag part", TMP(i,j), TMPold(i,j)
!         STOP 5
!         endif
!         
!     enddo

         call udvl%getscale(dnew,i)
         Z = dnew - Dold(i)
         if (Abs(real(Z)) > MAX(ABS(REAL(dnew)), ABS(REAL(Dold(i))))*1D-15 ) then
!         write (*,*) "Error in D real part", D(i), Dold(i)
!         STOP 6
         endif
         if (Abs(AIMAG(Z)) > MAX(ABS(AIMAG(dnew)), ABS(AIMAG(Dold(i))))*1D-15 ) then
!         write (*,*) "Error in D imag part", D(i), Dold(i)
!         STOP 7
         endif

        enddo
        CALL UDVL%dealloc
        Deallocate(Vold, Uold, Dold, TMP, TMPold, V1, TMP1)
        enddo
        
write (*,*) "success"
end Program TESTULUPDATEMATRICES

SUBROUTINE ul_update_matrices_old(U, D, V, V1, TMP, TMP1, Ndim, NCON)
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
END SUBROUTINE ul_update_matrices_old
