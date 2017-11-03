!  Copyright (C) 2017 The ALF project
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
        call ZGEQP3(Ndim, N_part, Mat, Ndim, IPVT, TAU(1), Z, -1, RWORK(1), INFO)
        LWORK = INT(DBLE(Z))
        ALLOCATE(WORK(LWORK))
        ! QR decomposition of Mat with full column pivoting, Mat * P = Q * R
        call ZGEQP3(Ndim, N_part, Mat, Ndim, IPVT, TAU(1), WORK(1), LWORK, RWORK(1), INFO)
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
END SUBROUTINE

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
END SUBROUTINE

End Module QDRP_mod
