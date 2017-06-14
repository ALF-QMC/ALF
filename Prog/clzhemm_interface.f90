!  Copyright (C) 2017 The ALF project
! 
!  This file is part of the ALF project.
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

module clALF
    interface
    subroutine initopenclandclblas(info) bind(c)
    use iso_c_binding
    IMPLICIT NONE
    INTEGER(c_int32_t), intent(out) :: info 
    end subroutine

    subroutine clalfzhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc, info) bind(c)
    use iso_c_binding
    IMPLICIT NONE
    INTEGER(c_int32_t), intent(out) :: info 
    CHARACTER, intent(in) :: side, uplo
    integer(c_int32_t), intent(in) :: m,n, lda, ldb, ldc
    complex (Kind=Kind(0.d0)), intent(in) :: alpha, beta
    COMPLEX(kind = kind(0.D0)), intent(in) :: A(LDA,*),B(LDB,*)
    COMPLEX(kind = kind(0.D0)), intent(inout) ::C(LDC,*)
    end subroutine

    subroutine teardown(info) bind(c)
    use iso_c_binding
    IMPLICIT NONE
    INTEGER(c_int32_t), intent(out) :: info 
    end subroutine
    end interface

END MODULE clALF
