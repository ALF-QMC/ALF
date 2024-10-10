!  Copyright (C) 2016 - 2020 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module provides a set of predefined lattices.
!>
!
!--------------------------------------------------------------------

module Predefined_Lattices

   use runtime_error_mod
   use Lattices_v3
   use iso_fortran_env, only: output_unit, error_unit
   implicit none

contains
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Definition of  a set of lattices: Square,  Honeycomb, Pi_Flux
!>
!> @param [in]  Latttice_type
!>\verbatim
!> Character(64)
!> Can take the values
!> Square,  Honeycomb, Pi_Flux
!> \endverbatim
!> @param [in]  L1, L2
!>\verbatim
!>    Integer
!>    Size of the lattice in units of the lattice constants
!>\endverbatim
!> @param [out]  Latt_unit
!>\verbatim
!>    Type (Unit_cell)
!>    The unit cell. Contains Norb, N_coord and positions of orbitals.
!>\endverbatim
!> @param [out]  Ndim
!>\verbatim
!>    Integer
!>    Number of orbitals
!>\endverbatim
!> @param [out]  List, Invlist
!>\verbatim
!>    Integer(:,:)
!>    List(I=1.. Ndim,1)    =   Unit cell of site I
!>    List(I=1.. Ndim,2)    =   Orbital index  of site I
!>    Invlist(1..Unit_cell,1..Orbital) = site I
!>\endverbatim
!> @param [out]  Latt
!>\verbatim
!>    Type(Lattice)
!>    Sets the lattice
!>\endverbatim
!> @param [out]  Latt_unit
!>\verbatim
!>    Type(Unit_cell)
!>    Sets the lattice
!>\endverbatim
!>
!-------------------------------------------------------------------
   subroutine Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)

      implicit none

      !Set the lattice
      character(len=64), intent(IN)                     :: Lattice_type
      integer, intent(IN)                                :: L1, L2
      integer, intent(OUT)                               :: Ndim
      integer, intent(OUT), dimension(:, :), allocatable  :: List, Invlist
      type(Unit_cell), intent(Out)                       :: Latt_Unit
      type(Lattice), intent(Out)                         :: Latt
      real(Kind=kind(0.d0))  :: A1_p(2), a2_p(2), L1_p(2), L2_p(2)
      integer :: I, nc, no, n

      select case (Lattice_type)
      case ("Square")
         if (L2 == 1 .and. L1 > 1) then
            Latt_Unit%N_coord = 1
         elseif (L2 > 1 .and. L1 > 1) then
            Latt_Unit%N_coord = 2
         else
            write (error_unit, *) 'For one-dimensional lattices set L2=1.'
            write (error_unit, *) 'You can also use use n_leg_ladder with n=1'
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if
         Latt_Unit%Norb = 1
         allocate (Latt_unit%Orb_pos_p(1, 2))
         Latt_Unit%Orb_pos_p(1, :) = 0.d0
         a1_p(1) = 1.0; a1_p(2) = 0.d0
         a2_p(1) = 0.0; a2_p(2) = 1.d0
         L1_p = dble(L1)*a1_p
         L2_p = dble(L2)*a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)
      case ('square_anisotropic')
         Latt_Unit%Norb = 1
         allocate (Latt_unit%Orb_pos_p(1, 2))
         Latt_Unit%Orb_pos_p(1, :) = 0.d0
         a1_p(1) = 1.0; a1_p(2) = 0.d0
         a2_p(1) = 0.0; a2_p(2) = 1.d0
         L1_p = dble(L1)*a1_p
         L2_p = dble(L2)*a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)
      case ("N_leg_ladder")
         a1_p(1) = 1.0; a1_p(2) = 0.d0
         a2_p(1) = 0.0; a2_p(2) = 1.d0
         L1_p = dble(L1)*a1_p
         L2_p = a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)

         Latt_Unit%Norb = L2
         Latt_Unit%N_coord = 1
         allocate (Latt_unit%Orb_pos_p(L2, 2))
         do no = 1, L2
            Latt_Unit%Orb_pos_p(no, 1) = 0.d0
            Latt_Unit%Orb_pos_p(no, 2) = real(no - 1, kind(0.d0))
         end do
      case ("Bilayer_square")
         a1_p(1) = 1.0; a1_p(2) = 0.d0
         a2_p(1) = 0.0; a2_p(2) = 1.d0
         L1_p = dble(L1)*a1_p
         L2_p = dble(L2)*a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)

         Latt_Unit%Norb = 2
         Latt_Unit%N_coord = 2
         allocate (Latt_unit%Orb_pos_p(2, 3))
         do no = 1, 2
            Latt_Unit%Orb_pos_p(no, 1) = 0.d0
            Latt_Unit%Orb_pos_p(no, 2) = 0.d0
            Latt_Unit%Orb_pos_p(no, 3) = real(1 - no, kind(0.d0))
         end do

      case ("Honeycomb")
         if (L1 == 1 .or. L2 == 1) then
            write (error_unit, *) 'The Honeycomb lattice cannot be one-dimensional.'
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if
         Latt_Unit%Norb = 2
         Latt_Unit%N_coord = 3
         a1_p(1) = 1.d0; a1_p(2) = 0.d0
         a2_p(1) = 0.5d0; a2_p(2) = sqrt(3.d0)/2.d0
         allocate (Latt_Unit%Orb_pos_p(2, 2))
         Latt_Unit%Orb_pos_p(1, :) = 0.d0
         Latt_Unit%Orb_pos_p(2, :) = (a2_p(:) - 0.5d0*a1_p(:))*2.d0/3.d0
         L1_p = dble(L1)*a1_p
         L2_p = dble(L2)*a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)
      case ("Bilayer_honeycomb")
         a1_p(1) = 1.d0; a1_p(2) = 0.d0
         a2_p(1) = 0.5d0; a2_p(2) = sqrt(3.d0)/2.d0
         L1_p = dble(L1)*a1_p
         L2_p = dble(L2)*a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)

         Latt_Unit%Norb = 4
         Latt_Unit%N_coord = 3
         allocate (Latt_unit%Orb_pos_p(4, 3))
         Latt_unit%Orb_pos_p = 0.d0
         do n = 1, 2
            Latt_Unit%Orb_pos_p(1, n) = 0.d0
            Latt_Unit%Orb_pos_p(2, n) = (a2_p(n) - 0.5d0*a1_p(n))*2.d0/3.d0
            Latt_Unit%Orb_pos_p(3, n) = 0.d0
            Latt_Unit%Orb_pos_p(4, n) = (a2_p(n) - 0.5d0*a1_p(n))*2.d0/3.d0
         end do
         Latt_Unit%Orb_pos_p(3, 3) = -1.d0
         Latt_Unit%Orb_pos_p(4, 3) = -1.d0
      case ("Pi_Flux")
         if (L1 == 1 .or. L2 == 1) then
            write (error_unit, *) 'The Pi Flux lattice cannot be one-dimensional.'
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if
         Latt_Unit%Norb = 2
         Latt_Unit%N_coord = 4
         a1_p(1) = 1.d0; a1_p(2) = 1.d0
         a2_p(1) = 1.d0; a2_p(2) = -1.d0
         allocate (Latt_Unit%Orb_pos_p(2, 2))
         Latt_Unit%Orb_pos_p(1, :) = 0.d0
         Latt_Unit%Orb_pos_p(2, :) = (a1_p(:) - a2_p(:))/2.d0
         L1_p = dble(L1)*(a1_p - a2_p)/2.d0
         L2_p = dble(L2)*(a1_p + a2_p)/2.d0
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)
      case default
         write (error_unit, *) "Predefined_Latt: Lattice not yet implemented!"
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end select
      ! Call Print_latt(Latt)
      ! This is for the orbital structure.

      Ndim = Latt%N*Latt_Unit%Norb
      allocate (List(Ndim, 2), Invlist(Latt%N, Latt_Unit%Norb))
      nc = 0
      do I = 1, Latt%N
         do no = 1, Latt_Unit%Norb
            ! For the Honeycomb and pi-flux lattices no = 1,2 corresponds to the A,and B sublattice.
            nc = nc + 1
            List(nc, 1) = I
            List(nc, 2) = no
            Invlist(I, no) = nc
         end do
      end do

   end subroutine Predefined_Latt

end module Predefined_Lattices
