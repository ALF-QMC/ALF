
!  Copyright (C) 2018-2026 The ALF project
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
     MODULE Matrix

!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Lightweight matrix container types and lifecycle helpers.
!
!> @details
!> Defines two simple square matrix containers with explicit memory
!> management:
!> - `Mat_C`: complex matrix
!> - `Mat_R`: real matrix
!>
!> Memory is allocated via `Make_Mat` and released via `Clear_Mat`.
!> Generic interfaces dispatch to complex/real implementations.
!> 
!--------------------------------------------------------------------

       
  !--------------------------------------------------------------------
  !> @brief Complex square matrix container.
  !>
  !> @details
  !> Stores a dynamically allocated complex matrix `el(dim,dim)` and the
  !> corresponding dimension.
  !--------------------------------------------------------------------
       Type Mat_C
         complex (Kind=Kind(0.d0)), pointer :: el(:,:)  ! Matrix elements
         Integer :: dim                                  ! Matrix dimension
       end Type Mat_C

  !--------------------------------------------------------------------
  !> @brief Real square matrix container.
  !>
  !> @details
  !> Stores a dynamically allocated real matrix `el(dim,dim)` and the
  !> corresponding dimension.
  !--------------------------------------------------------------------
       Type Mat_R
         Real (Kind=Kind(0.d0)), pointer :: el(:,:)  ! Matrix elements
         Integer :: dim                               ! Matrix dimension
       end Type Mat_R

       !> Generic constructor for matrix containers (real/complex)
       Interface Make_Mat
          module procedure constructor_C, constructor_R
       end Interface
       !> Generic destructor for matrix containers (real/complex)
       Interface Clear_Mat
          module procedure Destroy_C, Destroy_R
       end Interface

       Contains

  !--------------------------------------------------------------------
  !> @brief Allocate and initialize complex matrix container.
  !>
  !> @param[out] Mat Complex matrix container to initialize
  !> @param[in] N Matrix dimension (allocates N x N)
  !>
  !> @details
  !> Allocates `Mat%el(N,N)`, initializes all entries to complex zero, and
  !> stores `Mat%dim = N`.
  !--------------------------------------------------------------------
         subroutine constructor_C(Mat,N)
          type (Mat_C) :: Mat
          Integer :: N
          ! Allocate complex matrix storage
           allocate (Mat%el(N,N))
          ! Initialize all entries to zero
           Mat%el = cmplx(0.D0,0.D0, kind(0.D0))
          ! Store matrix dimension
           Mat%dim = N
         end subroutine constructor_C

  !--------------------------------------------------------------------
  !> @brief Allocate and initialize real matrix container.
  !>
  !> @param[out] Mat Real matrix container to initialize
  !> @param[in] N Matrix dimension (allocates N x N)
  !>
  !> @details
  !> Allocates `Mat%el(N,N)`, initializes all entries to zero, and stores
  !> `Mat%dim = N`.
  !--------------------------------------------------------------------
         subroutine constructor_R(Mat,N)
          type (Mat_R) :: Mat
          Integer :: N
          ! Allocate real matrix storage
           allocate (Mat%el(N,N))
          ! Initialize all entries to zero
           Mat%el = 0.d0
          ! Store matrix dimension
           Mat%dim = N
         end subroutine constructor_R

  !--------------------------------------------------------------------
  !> @brief Deallocate complex matrix container storage.
  !>
  !> @param[in,out] Mat Complex matrix container
  !
  !> @details
  !> Releases allocated memory for `Mat%el`.
  !--------------------------------------------------------------------
         subroutine Destroy_C(Mat)
           type (Mat_C) :: Mat
          ! Release matrix storage
           deallocate (Mat%el)
         end subroutine Destroy_C

  !--------------------------------------------------------------------
  !> @brief Deallocate real matrix container storage.
  !>
  !> @param[in,out] Mat Real matrix container
  !
  !> @details
  !> Releases allocated memory for `Mat%el`.
  !--------------------------------------------------------------------
         subroutine Destroy_R(Mat)
           type (Mat_R) :: Mat
          ! Release matrix storage
           deallocate (Mat%el)
         end subroutine Destroy_R
     end MODULE Matrix



