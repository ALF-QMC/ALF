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
        Real (Kind=Kind(0.d0))  :: A1_p(2), a2_p(2), L1_p(2), L2_p(2)
        Integer :: I, nc, no,n, lx, ly, no_y, no_x

        select case (Lattice_type)
        case("Square")
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
        case("N_leg_ladder")
           a1_p(1) =  1.0  ; a1_p(2) =  0.d0
           a2_p(1) =  0.0  ; a2_p(2) =  1.d0
           L1_p    =  dble(L1)*a1_p
           L2_p    =           a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

           Latt_Unit%Norb     = L2
           Latt_Unit%N_coord  = 1
           Allocate (Latt_unit%Orb_pos_p(L2,2))
           do no = 1,L2
              Latt_Unit%Orb_pos_p(no,1) = 0.d0
              Latt_Unit%Orb_pos_p(no,2) = real(no-1,kind(0.d0))
           enddo
        case("Bilayer_square")
           a1_p(1) =  1.0  ; a1_p(2) =  0.d0
           a2_p(1) =  0.0  ; a2_p(2) =  1.d0
           L1_p    =  dble(L1)*a1_p
           L2_p    =  dble(L2)*a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

           Latt_Unit%Norb     = 2
           Latt_Unit%N_coord  = 2
           Allocate (Latt_unit%Orb_pos_p(2,3))
           do no = 1,2
              Latt_Unit%Orb_pos_p(no,1) = 0.d0
              Latt_Unit%Orb_pos_p(no,2) = 0.d0
              Latt_Unit%Orb_pos_p(no,3) = real(1-no,kind(0.d0))
           enddo

        case("Honeycomb")
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
           Latt_Unit%Orb_pos_p(2,:) = (a2_p(:) - 0.5D0*a1_p(:) ) * 2.D0/3.D0
           L1_p    =  dble(L1) * a1_p
           L2_p    =  dble(L2) * a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case("Bilayer_honeycomb")
           a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
           a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0
           L1_p    =  dble(L1)*a1_p
           L2_p    =  dble(L2)*a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

           Latt_Unit%Norb     = 4
           Latt_Unit%N_coord  = 3
           Allocate (Latt_unit%Orb_pos_p(4,3))
           Latt_unit%Orb_pos_p = 0.d0
           do n = 1,2
              Latt_Unit%Orb_pos_p(1,n) = 0.d0
              Latt_Unit%Orb_pos_p(2,n) = (a2_p(n) - 0.5D0*a1_p(n) ) * 2.D0/3.D0
              Latt_Unit%Orb_pos_p(3,n) = 0.d0
              Latt_Unit%Orb_pos_p(4,n) = (a2_p(n) - 0.5D0*a1_p(n) ) * 2.D0/3.D0
           Enddo
           Latt_Unit%Orb_pos_p(3,3) = -1.d0
           Latt_Unit%Orb_pos_p(4,3) = -1.d0
        case("Pi_Flux")
           Latt_Unit%Norb    = 2
           Latt_Unit%N_coord = 4
           a1_p(1) = 1.d0; a1_p(2) = 0.d0
           a2_p(1) = 0.d0; a2_p(2) = 1.d0
           Allocate (Latt_Unit%Orb_pos_p(2,2))
           Latt_Unit%Orb_pos_p(1,:) = 0.d0
           Latt_Unit%Orb_pos_p(2,:) = 0.5d0
           L1_p    =  dble(L1) * a1_p
           L2_p    =  dble(L2) * a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case("Pi_Flux_ob")
           a1_p(1) = 1.d0; a1_p(2) = 0.d0
           a2_p(1) = 0.d0; a2_p(2) = 1.d0
           
           L1_p    =  dble(L1) * a1_p
           L2_p    =  dble( 1) * a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )

           Latt_Unit%Norb    = 2*L2
           Latt_Unit%N_coord = 2
           Allocate (Latt_Unit%Orb_pos_p(2*L2,2))
           Latt_unit%Orb_pos_p = 0.d0
           do nc = 1, 2*L2
              ly = (nc-1)/2+1
              Latt_Unit%Orb_pos_p(nc,:) = 0.d0 + (ly-1)*a2_p(:)
              if ( mod(nc,2) .eq. 0 ) then
                  Latt_Unit%Orb_pos_p(nc,:) = Latt_Unit%Orb_pos_p(nc,:)+0.5d0
              endif
           enddo
        case default
           Write(error_unit,*) "Predefined_Latt: Lattice not yet implemented!"
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        end select
        ! Call Print_latt(Latt)
        ! This is for the orbital structure.


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
