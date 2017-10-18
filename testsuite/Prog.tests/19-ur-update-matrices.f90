Program TESTURUPDATEMATRICES
Use UDV_State_mod
implicit none
interface
SUBROUTINE ur_update_matrices_old(U, D, V, V1, TMP, TMP1, Ndim, NCON)
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
        TYPE(UDV_State) :: udvl 
        COMPLEX(Kind=Kind(0.D0)) :: Z, dnew
        
        do Ndim = 10, 200,10
        CALL udvl%alloc(Ndim)
        Allocate (Uold(Ndim, Ndim), Vold(Ndim, Ndim))
        ALLocate(TMP(Ndim, Ndim), TMPold(Ndim, Ndim), Dold(Ndim), TMP1(Ndim, Ndim), V1(Ndim, Ndim))
        Uold = 0.D0
        Vold = 0.D0
        do i = 1, Ndim
        do j = 1, Ndim
        TMP(i, j) = 1.D0/(i+j)
        enddo
        Dold(i) = i*i
        Uold(i,i) = 1.D0
        Vold(i,i) = 1.D0
        call udvl%setscale(Dold(i),i)
        enddo
        udvl%U = Uold
        udvl%V = Vold
        TMPold = TMP
        NCON = 1
        CALL UDVL%matmultleft(TMP, TMP1, NCON)
        call ur_update_matrices_old(Uold, Dold, Vold, V1, TMPold, TMP1, Ndim, NCON)
        
        ! compare
        do i = 1, Ndim
        do j = 1, Ndim
        Z = udvl%V(i,j) - Vold(i,j)
        
        if (Abs(real(Z)) > MAX(ABS(REAL(udvl%V(i, j))), ABS(REAL(Vold(i, j))))*1D-15 ) then
!        write (*,*) "Error in V real part", V(i,j), Vold(i,j)
!        STOP 2
        endif
        if (Abs(AIMAG(Z)) > MAX(ABS(AIMAG(udvl%V(i, j))), ABS(AIMAG(Vold(i, j))))*1D-15 ) then
!        write (*,*) "Error in V imag part", V(i,j), Vold(i,j)
!        STOP 3
        endif
        
        Z = udvl%U(i,j) - Uold(i,j)
        if (Abs(real(Z)) > MAX(ABS(REAL(udvl%U(i, j))), ABS(REAL(Uold(i, j))))*1D-15 ) then
!        write (*,*) "Error in U real part", U(i,j), Uold(i,j)
!        STOP 4
        endif
        if (Abs(AIMAG(Z)) > MAX(ABS(AIMAG(TMP(i, j))), ABS(AIMAG(TMPold(i, j))))*1D-15 ) then
!        write (*,*) "Error in TMP imag part", TMP(i,j), TMPold(i,j)
!        STOP 5
        endif
        
        enddo
        
        call udvl%getscale(dnew,i)
        Z = dnew - Dold(i)
        if (Abs(real(Z)) > MAX(ABS(REAL(dnew)), ABS(REAL(Dold(i))))*1D-15) then
!        write (*,*) "Error in D real part", D(i), Dold(i)
!        STOP 6
        endif
        if (Abs(AIMAG(Z)) > MAX(ABS(AIMAG(dnew)), ABS(AIMAG(Dold(i))))*1D-15 ) then
!        write (*,*) "Error in D imag part", D(i), Dold(i)
!        STOP 7
        endif
        
        enddo
        CALL udvl%dealloc
        Deallocate(Vold, Uold, Dold, TMP, TMPold, V1, TMP1)
        enddo
        
write (*,*) "success"
end Program TESTURUPDATEMATRICES

SUBROUTINE ur_update_matrices_old(U, D, V, V1, TMP, TMP1, Ndim, NCON)
        Use UDV_Wrap_mod
        Implicit None
        INTEGER, intent(in) :: Ndim, NCON
        COMPLEX (Kind=Kind(0.d0)) :: U(Ndim,Ndim), V(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim),TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) :: D(Ndim)
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        INTEGER :: n

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        CALL ZGEMM('N', 'N', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, U(1, 1), Ndim, beta, TMP1, Ndim)
        DO n = 1,NDim
            TMP1(:, n) = TMP1(:, n)*D(n)
        ENDDO
        CALL UDV_WRAP_Pivot(TMP1, U, D , V1,NCON,Ndim,Ndim)
        CALL ZGEMM('N', 'N', Ndim, Ndim, Ndim, Z_ONE, V1, Ndim, V(1, 1), Ndim, beta, TMP1, Ndim)
        V = TMP1
END SUBROUTINE ur_update_matrices_old
