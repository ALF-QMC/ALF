!  Copyright (C) 2020 The ALF project
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
!       https://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage https://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


module matTypes_mod
    use ContainerElementBase_mod
    implicit none

    private
    public :: RealMat, CmplxMat

    type, extends(ContainerElementBase) :: RealMat
        Real(kind=kind(0.d0)), allocatable, dimension(:,:) :: mat
        Integer :: m, n
    contains
        procedure :: init => RealMat_init
        procedure :: simt => RealMat_simt
        procedure :: rmult => RealMat_rmult
        procedure :: lmult => RealMat_lmult
        procedure :: rmultinv => RealMat_rmultinv
        procedure :: lmultinv => RealMat_lmultinv
    end type RealMat

    type, extends(ContainerElementBase) :: CmplxMat
        Complex(kind=kind(0.d0)),allocatable, dimension(:,:):: mat
        Integer :: m, n
    contains
        procedure :: init => CmplxMat_init
        procedure :: simt => CmplxMat_simt
        procedure :: rmult => CmplxMat_rmult
        procedure :: lmult => CmplxMat_lmult
        procedure :: rmultinv => CmplxMat_rmultinv
        procedure :: lmultinv => CmplxMat_lmultinv
    end type CmplxMat

contains
    subroutine RealMat_init(this, arg)
        class(RealMat) :: this
        Real(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
        Integer :: i,j
        this%mat = arg !copy argument to local storage
        this%m = size(arg, 1)
        this%n = size(arg, 2)

!                 do i = 1, size(arg,1)
! write (*,*) (this%mat(i,j), j = 1, size(arg,2) )
! enddo
    end subroutine
    
    subroutine RealMat_simt(this, arg)
        class(RealMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine
    
    subroutine RealMat_rmult(this, arg)
        class(RealMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: out
        Real(kind=kind(0.D0)), allocatable, dimension(:) :: rwork
        Integer :: i, j, sz1, sz2
        
        sz1 = size(arg, 1)
        sz2 = size(arg, 2)
        allocate(out(sz1, sz2), rwork(2*sz1*sz2))
        call zlacrm(sz1, sz2, arg, sz1, this%mat, this%m, out, this%m, rwork) ! zlarcm assumes mat to be square
        arg = out
        deallocate(out, rwork)
    end subroutine
    
    subroutine RealMat_rmultinv(this, arg)
        class(RealMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: out
        Real(kind=kind(0.D0)), allocatable, dimension(:) :: rwork
        Integer :: i, j, sz1, sz2
        
    end subroutine
    
    subroutine RealMat_lmult(this, arg)
        class(RealMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine

    subroutine RealMat_lmultinv(this, arg)
        class(RealMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine
    
    subroutine CmplxMat_init(this, arg)
        class(CmplxMat) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
        Integer :: i,j
!         do i = 1, size(arg,1)
! write (*,*) (this%mat(i,j), j = 1, size(arg,2) )
! enddo
        this%m = size(arg, 1)
        this%n = size(arg, 2)
        if (allocated(this%mat)) deallocate(this%mat)
        allocate(this%mat(this%m, this%n))
        do i = 1, this%m
        do j = 1, this%n
        this%mat(i, j) = arg(i, j)
        enddo
        enddo
!        this%mat = arg !copy argument to local storage
        
    end subroutine

    subroutine CmplxMat_simt(this, arg)
        class(CmplxMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
    end subroutine
    
    subroutine CmplxMat_rmult(this, arg)
        class(CmplxMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: out
        Complex(kind=kind(0.d0)) :: alpha, zero
        Integer :: i, j, sz1, sz2
        
        alpha = 1.0
        zero = 0
        sz1 = size(arg, 1)
        sz2 = size(arg, 2)
        allocate(out(sz1, sz2))
        call zhemm('R', 'U', sz1, sz2, alpha, arg, sz1, this%mat, this%m, zero, out, sz1)
        arg = out
        deallocate(out)
    end subroutine
    
    subroutine CmplxMat_rmultinv(this, arg)
        class(CmplxMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: out
        Complex(kind=kind(0.d0)) :: alpha, zero
        Integer :: i, j, sz1, sz2
        
    end subroutine
    
    subroutine CmplxMat_lmult(this, arg)
        class(CmplxMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
    end subroutine

    subroutine CmplxMat_lmultinv(this, arg)
        class(CmplxMat), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
    end subroutine
end module matTypes_mod
