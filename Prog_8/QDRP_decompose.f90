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

! This constructs a decompostion Mat = Q D R P^*

Module QDRP_mod

Contains

SUBROUTINE QDRP_decompose(Ndim, Mat, D, IPVT, TAU, WORK, LWORK)
Implicit None
Integer, intent(in) :: Ndim
Integer, intent(inout) :: LWORK
Integer, Dimension(:), intent(inout) :: IPVT
COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(inout) :: Mat
COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(inout) :: D, TAU
COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(INOUT), Allocatable :: WORK

COMPLEX(Kind=Kind(0.d0)), Dimension(:), Allocatable :: RWORK
COMPLEX(Kind=Kind(0.d0)) :: Z
Integer :: info, i, j
Real(Kind=Kind(0.d0)) :: X

        ALLOCATE(RWORK(2*Ndim))
        call ZGEQP3(Ndim, Ndim, Mat(1, 1), Ndim, IPVT, TAU(1), Z, -1, RWORK(1), INFO)
        LWORK = INT(DBLE(Z))
        ALLOCATE(WORK(LWORK))
        ! QR decomposition of TPUP with full column pivoting, AP = QR
        call ZGEQP3(Ndim, Ndim, Mat(1, 1), Ndim, IPVT, TAU(1), WORK(1), LWORK, RWORK(1), INFO)
        DEALLOCATE(RWORK)
        ! separate off D3
        do i = 1, Ndim
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
            do j = i, Ndim
                Mat(i, j) = Mat(i, j) / X
            enddo
        enddo
END SUBROUTINE

End Module QDRP_mod
