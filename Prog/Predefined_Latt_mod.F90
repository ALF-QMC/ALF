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


    Module Predefined_Lattices

      use runtime_error_mod
      Use Lattices_v3
      Use files_mod
      use iso_fortran_env, only: output_unit, error_unit
      Implicit none


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
      Subroutine Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit )

        Implicit none

        !Set the lattice
        Character (len=64), Intent(IN)                     :: Lattice_type
        Integer, Intent(IN)                                :: L1,L2
        Integer, Intent(OUT)                               :: Ndim
        Integer, Intent(OUT), Dimension(:,:), allocatable  :: List, Invlist
        Type(Unit_cell), Intent(Out)                       :: Latt_Unit
        Type(Lattice), Intent(Out)                         :: Latt
        ! Local variables: primitive lattice vectors (a1_p, a2_p) and simulation-cell
        ! boundary vectors (L1_p, L2_p) are all in Cartesian coordinates (Angstrom-like units).
        Real (Kind=Kind(0.d0))  :: A1_p(2), a2_p(2), L1_p(2), L2_p(2)
        Integer :: I, nc, no, n  ! nc = running site/orbital index; no = orbital index within unit cell

        ! Dispatch to the appropriate lattice geometry.  Each branch sets:
        !   Latt_Unit%Norb    = number of orbitals per unit cell
        !   Latt_Unit%N_coord = coordination number (relevant for checkerboard)
        !   Latt_Unit%Orb_pos_p(no, 1:2 or 1:3) = orbital positions in Cartesian coords
        !     (third coordinate, if present, encodes a layer index for bilayer geometries)
        !   a1_p, a2_p = primitive lattice vectors
        !   L1_p, L2_p = simulation-cell boundary vectors
        select case (str_to_upper(Lattice_type))
        case("SQUARE")
           ! Single-orbital 2D square lattice with orthogonal primitive vectors.
           ! N_coord = 2 (two distinct bond directions: along a1 and a2) for 2D,
           ! reduced to N_coord = 1 for a chain (L2 = 1).
           If ( L2==1 .and. L1 > 1 ) then
              Latt_Unit%N_coord   = 1
           elseif (L2 >1 .and. L1 > 1) then
              Latt_Unit%N_coord   = 2
           else
              Write(error_unit,*) 'For one-dimensional lattices set L2=1.'
              Write(error_unit,*) 'You can also use use n_leg_ladder with n=1'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           Latt_Unit%Norb      = 1
           Allocate (Latt_unit%Orb_pos_p(1,2))
           Latt_Unit%Orb_pos_p(1,:) = 0.d0
           a1_p(1) =  1.0  ; a1_p(2) =  0.d0
           a2_p(1) =  0.0  ; a2_p(2) =  1.d0
           L1_p    =  dble(L1)*a1_p
           L2_p    =  dble(L2)*a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case("N_LEG_LADDER")
           ! N-leg ladder: L2 legs of length L1 rungs.  The L2 orbitals per unit cell
           ! are stacked along the a2 direction (rung direction); lattice spans only
           ! one unit cell in a2 (L2_p = a2_p), so the full ladder is encoded via
           ! the multi-orbital structure rather than replicated unit cells.
           a1_p(1) =  1.0  ; a1_p(2) =  0.d0
           a2_p(1) =  0.0  ; a2_p(2) =  1.d0
           L1_p    =  dble(L1)*a1_p
           L2_p    =           a2_p   ! single unit cell in rung direction
           Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

           Latt_Unit%Norb     = L2   ! one orbital per leg
           Latt_Unit%N_coord  = 1
           Allocate (Latt_unit%Orb_pos_p(L2,2))
           do no = 1,L2
              Latt_Unit%Orb_pos_p(no,1) = 0.d0
              Latt_Unit%Orb_pos_p(no,2) = real(no-1,kind(0.d0))  ! offset each leg by one unit
           enddo
        case("BILAYER_SQUARE")
           ! Two-layer square lattice.  Orbitals 1 and 2 sit at the same (x,y) position
           ! but at z = 0 and z = -1 respectively (encoded in Orb_pos_p(:,3)).
           a1_p(1) =  1.0  ; a1_p(2) =  0.d0
           a2_p(1) =  0.0  ; a2_p(2) =  1.d0
           L1_p    =  dble(L1)*a1_p
           L2_p    =  dble(L2)*a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

           Latt_Unit%Norb     = 2
           Latt_Unit%N_coord  = 2
           Allocate (Latt_unit%Orb_pos_p(2,3))  ! third coord encodes layer
           do no = 1,2
              Latt_Unit%Orb_pos_p(no,1) = 0.d0
              Latt_Unit%Orb_pos_p(no,2) = 0.d0
              Latt_Unit%Orb_pos_p(no,3) = real(1-no,kind(0.d0))  ! layer 1: z=0, layer 2: z=-1
           enddo
        case("TRIANGULAR")
           ! Single-orbital 2D triangular lattice with primitive vectors
           ! a1 = (1, 0) and a2 = (1/2, sqrt(3)/2).  The non-orthogonal a2
           ! gives the 60-degree bond and 3 distinct NN directions (N_coord=3).
           If (L1==1 .or. L2==1 ) then
              Write(error_unit,*) 'The triangular lattice cannot be one-dimensional.'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           Latt_Unit%Norb    = 1
           Latt_Unit%N_coord = 3
           a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
           a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0  ! 60-degree angle
           Allocate (Latt_Unit%Orb_pos_p(1,2))
           Latt_Unit%Orb_pos_p = 0.d0
           L1_p    =  dble(L1) * a1_p
           L2_p    =  dble(L2) * a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         case("KAGOME")
            ! Kagome lattice: 3 orbitals per unit cell placed at the midpoints of the
            ! triangular-lattice bonds (a1/2 and a2/2 from the origin).  N_coord=4
            ! because the unit cell has 4 distinct NN bond types.
            If (L1==1 .or. L2==1 ) then
               Write(error_unit,*) 'The kagome lattice cannot be one-dimensional.'
               CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            endif
            Latt_Unit%Norb    = 3
            Latt_Unit%N_coord = 4
            a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
            a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0  ! shared with triangular
            Allocate (Latt_Unit%Orb_pos_p(3,2))
            Latt_Unit%Orb_pos_p = 0.d0
            Latt_Unit%Orb_pos_p(2,:) = a1_p(:)/2.d0  ! midpoint of a1 bond
            Latt_Unit%Orb_pos_p(3,:) = a2_p(:)/2.d0  ! midpoint of a2 bond
            L1_p    =  dble(L1) * a1_p
            L2_p    =  dble(L2) * a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case("HONEYCOMB")
           ! Honeycomb lattice: 2-orbital unit cell with A (orbital 1) at origin and
           ! B (orbital 2) displaced by (2/3)(a2 - a1/2).  This places B at the
           ! nearest-neighbour position within the cell.  Same primitive vectors as
           ! triangular; sqrt(3)/2 gives the correct in-plane geometry.
           If (L1==1 .or. L2==1 ) then
              Write(error_unit,*) 'The Honeycomb lattice cannot be one-dimensional.'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           Latt_Unit%Norb    = 2
           Latt_Unit%N_coord = 3
           a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
           a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0
           Allocate (Latt_Unit%Orb_pos_p(2,2))
           Latt_Unit%Orb_pos_p(1,:) = 0.d0
           Latt_Unit%Orb_pos_p(2,:) = (a2_p(:) - 0.5D0*a1_p(:) ) * 2.D0/3.D0  ! B-sublattice position
           L1_p    =  dble(L1) * a1_p
           L2_p    =  dble(L2) * a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case("BILAYER_HONEYCOMB")
           ! Bilayer honeycomb: 4 orbitals per unit cell, two per layer.
           !   Orbitals 1,2 = layer 1 (A/B sublattice, z = 0)
           !   Orbitals 3,4 = layer 2 (A/B sublattice, z = -1)
           ! The in-plane positions of orbitals 3,4 mirror those of 1,2;
           ! the third column of Orb_pos_p encodes the layer separation.
           a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
           a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0
           L1_p    =  dble(L1)*a1_p
           L2_p    =  dble(L2)*a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

           Latt_Unit%Norb     = 4
           Latt_Unit%N_coord  = 3
           Allocate (Latt_unit%Orb_pos_p(4,3))  ! third coord = layer offset
           Latt_unit%Orb_pos_p = 0.d0
           do n = 1,2
              Latt_Unit%Orb_pos_p(1,n) = 0.d0
              Latt_Unit%Orb_pos_p(2,n) = (a2_p(n) - 0.5D0*a1_p(n) ) * 2.D0/3.D0  ! B position, layer 1
              Latt_Unit%Orb_pos_p(3,n) = 0.d0
              Latt_Unit%Orb_pos_p(4,n) = (a2_p(n) - 0.5D0*a1_p(n) ) * 2.D0/3.D0  ! B position, layer 2
           Enddo
           Latt_Unit%Orb_pos_p(3,3) = -1.d0  ! layer 2 at z = -1
           Latt_Unit%Orb_pos_p(4,3) = -1.d0
        case("PI_FLUX")
           ! Pi-flux square lattice: a pi flux is threaded per plaquette by choosing
           ! non-orthogonal primitive vectors a1=(1,1) and a2=(1,-1) and then defining
           ! the simulation-cell vectors as L1_p = L1*(a1-a2)/2 = L1*(0,1) and
           ! L2_p = L2*(a1+a2)/2 = L2*(1,0).  The 2-orbital unit cell has A at the
           ! origin and B displaced by (a1-a2)/2.
           If (L1==1 .or. L2==1 ) then
              Write(error_unit, *) 'The Pi Flux lattice cannot be one-dimensional.'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           Latt_Unit%Norb    = 2
           Latt_Unit%N_coord = 4  ! 4 NN bonds (all equivalent by symmetry), no checkerboard phase
           a1_p(1) =  1.D0   ; a1_p(2) =   1.d0
           a2_p(1) =  1.D0   ; a2_p(2) =  -1.d0
           Allocate (Latt_Unit%Orb_pos_p(2,2))
           Latt_Unit%Orb_pos_p(1,:) = 0.d0
           Latt_Unit%Orb_pos_p(2,:) = (a1_p(:) - a2_p(:))/2.d0  ! B sublattice: (0,1)
           L1_p    =  dble(L1) * (a1_p - a2_p)/2.d0  ! = L1*(0,1)
           L2_p    =  dble(L2) * (a1_p + a2_p)/2.d0  ! = L2*(1,0)
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case default
           Write(error_unit,*) "Predefined_Latt: Lattice not yet implemented!"
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        end select
        ! Call Print_latt(Latt)

        ! Build the site-to-(unit-cell, orbital) mapping tables.
        ! After this loop:
        !   List(nc, 1) = unit-cell index I    for composite site nc
        !   List(nc, 2) = orbital index    no   for composite site nc
        !   Invlist(I, no) = nc   (the inverse mapping)
        ! nc runs from 1 to Ndim = Latt%N * Latt_Unit%Norb.
        ! For Honeycomb and Pi_Flux: orbital 1 = A sublattice, orbital 2 = B sublattice.
        Ndim = Latt%N*Latt_Unit%Norb
        Allocate (List(Ndim,2), Invlist(Latt%N,Latt_Unit%Norb))
        nc = 0
        Do I = 1,Latt%N
           Do no = 1,Latt_Unit%Norb
              ! For the Honeycomb and pi-flux lattices no = 1,2 corresponds to the A,and B sublattice.
              nc = nc + 1
              List(nc,1) = I
              List(nc,2) = no
              Invlist(I,no) = nc
           Enddo
        Enddo

      end Subroutine Predefined_Latt


    end Module Predefined_Lattices
