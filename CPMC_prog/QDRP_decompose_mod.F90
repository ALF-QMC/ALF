!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This constructs a decompostion Mat = Q D R P^* using a pivoted QR decomposition
!--------------------------------------------------------------------
Module QDRP_mod

Contains

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> ! This constructs a decompostion Mat = Q D R P^* using a pivoted QR decomposition
!
!> @param Ndim[in] The size of the involved matrices
!> @param Mat[inout] The matrix that we want to decompose. The Householder reflectors and the upper
!>            triangular matrix are returned in Mat on exit
!> @param D[inout] The diagonal elements
!> @param IPVT[inout] A Pivoting vector that collects the permutations. Can be used by ?LAPMR and ?LAPMT
!> @param TAU[inout] The scalar factors for the Householder decomposition
!> @param WORK[inout] work memory. We query and allocate it in this routine. Needs to be deallocated outside.
!> @param LWORK[inout] optimal size of the work memory.
!--------------------------------------------------------------------
    SUBROUTINE QDRP_decompose(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK)
      Implicit None
      Integer, intent(in) :: Ndim
      Integer, intent(in) :: N_part
      Integer, intent(inout) :: LWORK
      Integer, Dimension(:), intent(inout), Allocatable :: IPVT
      COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(inout) :: Mat
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(inout) :: D
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(inout), Allocatable :: TAU
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(INOUT), Allocatable :: WORK
      
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Allocatable :: RWORK
      COMPLEX(Kind=Kind(0.d0)) :: Z
      Integer :: info, i, j
      Real(Kind=Kind(0.d0)) :: X
      
      ALLOCATE(RWORK(2*Ndim))
      ! Query optimal amount of memory
      call ZGEQP3(Ndim, N_part, Mat(1,1), Ndim, IPVT, TAU(1), Z, -1, RWORK(1), INFO)
      LWORK = INT(DBLE(Z))
      ALLOCATE(WORK(LWORK))
      ! QR decomposition of Mat with full column pivoting, Mat * P = Q * R
      call ZGEQP3(Ndim, N_part, Mat(1,1), Ndim, IPVT, TAU(1), WORK(1), LWORK, RWORK(1), INFO)
      DEALLOCATE(RWORK)
      ! separate off D
      do i = 1, N_part
         ! plain diagonal entry
         X = ABS(Mat(i, i))
         !             ! a inf-norm
         !             X = TPUP(i, i+izamax(Ndim+1-i, TPUP(i, i), Ndim)-1)
         !             ! another inf-norm
         !             X = TPUP(i, i-1+izmax1(Ndim+1-i, TPUP(i, i), Ndim))
         !             ! 1-norm
         !            X = DZSUM1(N_size+1-i, TPUP(i, i), N_size)
         ! 2-norm
         !            X = DZNRM2(N_size+1-i, TPUP(i, i), N_size)
         D(i) = X
         do j = i, N_part
            Mat(i, j) = Mat(i, j) / X
         enddo
      enddo
    END SUBROUTINE QDRP_decompose
    
    SUBROUTINE Pivot_phase(Phase, IPVT, N_size)
      Implicit none
      COMPLEX(kind=kind(0.d0)), Intent(INOUT) :: Phase
      Integer, Dimension(:), Intent(IN)       :: IPVT
      Integer,               Intent(IN)       :: N_size
      
      Integer:: i, next, L, VISITED(N_size)
      
      VISITED=0
      do i = 1, N_size
         if (VISITED(i) .eq. 0) then
            next = i
            L = 0
            do while (VISITED(next) .eq. 0)
               L = L + 1
               VISITED(next) = 1
               next = IPVT(next)
            enddo
            if(MOD(L, 2) .eq. 0) then
               PHASE = -PHASE
            endif
         endif
      enddo
    END SUBROUTINE Pivot_phase
    
  End Module QDRP_mod
