!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This constructs a decompostion Mat = Q D R P^* using a pivoted QR decomposition
!--------------------------------------------------------------------
module QDRP_mod

contains

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
   subroutine QDRP_decompose(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK)
      implicit none
      integer, intent(in) :: Ndim
      integer, intent(in) :: N_part
      integer, intent(inout) :: LWORK
      integer, dimension(:), intent(inout), allocatable :: IPVT
      complex(Kind=kind(0.d0)), dimension(:, :), intent(inout) :: Mat
      complex(Kind=kind(0.d0)), dimension(:), intent(inout) :: D
      complex(Kind=kind(0.d0)), dimension(:), intent(inout), allocatable :: TAU
      complex(Kind=kind(0.d0)), dimension(:), intent(INOUT), allocatable :: WORK

      complex(Kind=kind(0.d0)), dimension(:), allocatable :: RWORK
      complex(Kind=kind(0.d0)) :: Z
      integer :: info, i, j
      real(Kind=kind(0.d0)) :: X

      allocate (RWORK(2*Ndim))
      ! Query optimal amount of memory
      call ZGEQP3(Ndim, N_part, Mat(1, 1), Ndim, IPVT, TAU(1), Z, -1, RWORK(1), INFO)
      LWORK = int(dble(Z))
      allocate (WORK(LWORK))
      ! QR decomposition of Mat with full column pivoting, Mat * P = Q * R
      call ZGEQP3(Ndim, N_part, Mat(1, 1), Ndim, IPVT, TAU(1), WORK(1), LWORK, RWORK(1), INFO)
      deallocate (RWORK)
      ! separate off D
      do i = 1, N_part
         ! plain diagonal entry
         X = abs(Mat(i, i))
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
            Mat(i, j) = Mat(i, j)/X
         end do
      end do
   end subroutine QDRP_decompose

   subroutine Pivot_phase(Phase, IPVT, N_size)
      implicit none
      complex(kind=kind(0.d0)), intent(INOUT) :: Phase
      integer, dimension(:), intent(IN)       :: IPVT
      integer, intent(IN)       :: N_size

      integer:: i, next, L, VISITED(N_size)

      VISITED = 0
      do i = 1, N_size
         if (VISITED(i) .eq. 0) then
            next = i
            L = 0
            do while (VISITED(next) .eq. 0)
               L = L + 1
               VISITED(next) = 1
               next = IPVT(next)
            end do
            if (mod(L, 2) .eq. 0) then
               PHASE = -PHASE
            end if
         end if
      end do
   end subroutine Pivot_phase

end module QDRP_mod
