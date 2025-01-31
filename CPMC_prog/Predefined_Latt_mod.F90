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
      real(Kind=kind(0.d0))  :: A1_p(2), a2_p(2), a3_p(2), L1_p(2), L2_p(2)
      integer :: I, nc, no, n, lx, ly, no_y, no_x

      select case (Lattice_type)
      case ("Pi_Flux")
         Latt_Unit%Norb = 2
         Latt_Unit%N_coord = 4
         a1_p(1) = 1.d0; a1_p(2) = 0.d0
         a2_p(1) = 0.d0; a2_p(2) = 1.d0
         allocate (Latt_Unit%Orb_pos_p(2, 2))
         Latt_Unit%Orb_pos_p(1, :) = 0.d0
         Latt_Unit%Orb_pos_p(2, :) = 0.5d0
         L1_p = dble(L1)*a1_p
         L2_p = dble(L2)*a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)
      case ("Pi_Flux_ob")
         a1_p(1) = 1.d0; a1_p(2) = 0.d0
         a2_p(1) = 0.d0; a2_p(2) = 1.d0
         !! for sublattice vec
         a3_p(1) = 0.5d0; a3_p(2) = 0.5d0

         L1_p = dble(L1)*a1_p
         L2_p = dble(1)*a2_p
         call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)

         Latt_Unit%Norb = 2*L2
         Latt_Unit%N_coord = 2
         allocate (Latt_Unit%Orb_pos_p(2*L2, 2))
         Latt_unit%Orb_pos_p = 0.d0
         do nc = 1, 2*L2
            ly = (nc - 1)/2+1
            Latt_Unit%Orb_pos_p(nc, :) = 0.d0 + dble(ly-1)*a2_p(:) 
            if (mod(nc, 2) .eq. 0) then
               Latt_Unit%Orb_pos_p(nc, :) = Latt_Unit%Orb_pos_p(nc, :) + a3_p(:)
            end if
         end do
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
