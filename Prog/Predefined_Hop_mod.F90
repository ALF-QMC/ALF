!  Copyright (C) 2016 - 2023 The ALF project
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
!> This module provides a set of predefined hopping matrices for the lattice types
!> supported by ALF (square, triangular, kagome, honeycomb, bilayer variants, N-leg ladder)
!> as well as a general framework for specifying the hopping Hamiltonian of a
!> translation-invariant multi-orbital system.
!>
!> The central data structure is Hopping_Matrix_type, which stores the bond list, hopping
!> amplitudes, flux parameters, and (optionally) the checkerboard decomposition families
!> needed for an efficient Suzuki-Trotter factorisation.
!>
!> Typical usage:
!> 1. Call one of the \c Set_Default_hopping_parameters_* routines to populate
!>    an allocatable array of type(Hopping_Matrix_type).
!> 2. Call Predefined_Hoppings_set_OPT to convert the hopping matrix into the
!>    array of Operator objects consumed by the QMC engine.
!>
!> @see Operator_mod
!> @see Lattices_v3
!--------------------------------------------------------------------


    Module Predefined_Hoppings

      Use runtime_error_mod
      Use Lattices_v3
      Use Operator_mod
      Use WaveFunction_mod
      Use MyMats
      use iso_fortran_env, only: output_unit, error_unit
      use Hamiltonian_main
      Implicit none

      Logical, private, save :: pinning_notice_issued = .false.
      Logical, private, save :: first_pinning_notice_issued = .true.

      Type Hopping_Matrix_type
         Integer                   :: N_bonds  !< Number of inequivalent bond types in the unit cell.
         Complex (Kind=Kind(0.d0)), pointer :: T    (:)  !< Hopping amplitude T(N_b) for each bond type; does not include on-site terms.
         Complex (Kind=Kind(0.d0)), pointer :: T_loc(:)  !< On-site (diagonal) term for each orbital, e.g. chemical potential.
         Integer                  , pointer :: list(:,:) !< Bond definition table: List(N_b,1)=orbital_from, List(N_b,2)=orbital_to, List(N_b,3)=n_1, List(N_b,4)=n_2.
         ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)
         Integer                   :: N_Phi    !< Number of magnetic flux quanta \f$N_\Phi\f$ threading the lattice (Landau gauge).
         Real    (Kind=Kind(0.d0)) :: Phi_X    !< Twist / Aharonov-Bohm flux along the a1 direction.
         Real    (Kind=Kind(0.d0)) :: Phi_Y    !< Twist / Aharonov-Bohm flux along the a2 direction.
         Logical                   :: Bulk     !< If .true., twist is inserted as a vector potential in the bulk; if .false., only at the boundary.

         ! For Checkerboard decomposition
         Integer                            :: N_Fam      !< Number of commuting families in the checkerboard decomposition.
         Integer                  , pointer :: L_Fam(:)   !< L_Fam(n): number of bonds in family n.
         Integer                  , pointer :: List_Fam(:,:,:) !< List_Fam(n,l,1:2): unit cell and bond index of the l-th bond in family n.
         Real    (Kind=Kind(0.d0)), pointer :: Prop_Fam(:) !< Trotter weight for family n; 1 for standard, 0.5 for symmetric satellite families.

         Integer, private         , allocatable :: Multiplicity(:) !< Number of times each orbital appears in the bond list; used when distributing on-site terms across checkerboard operators.
      End type Hopping_Matrix_Type



    contains
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Deallocates this
!>
!
!--------------------------------------------------------------------
      Subroutine Predefined_hoppings_clear(this)

        Implicit none

        Type  (Hopping_Matrix_type)   , allocatable    :: this(:)
        Integer :: n
        If ( allocated(this) ) then
           do n = 1, size(This,1)
              deallocate (this(n)%T,this(n)%T_loc,this(n)%list)
           enddo
           deallocate (this(1)%L_Fam, this(1)%List_Fam, this(1)%Prop_Fam )
           if( allocated(this(1)%Multiplicity) ) deallocate(this(1)%Multiplicity)
        endif

      end Subroutine Predefined_hoppings_clear
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!>  Checks if the Hopping is
!>  zero     -->  inquire_hop = 0
!>  diagonal -->  inquire_hop = 1
!>  full     -->  inquire_hop = 2
!
!--------------------------------------------------------------------
      Integer function  inquire_hop(this)

        Implicit none
        Type  (Hopping_Matrix_type), Intent(In)  :: this(:)

        Real (Kind=Kind(0.d0)) :: Xmax_loc, Xmax_hop, Zero, X
        Integer :: nc, nf

        Zero     =  1.D-10
        Xmax_loc =  0.d0
        Xmax_hop =  0.d0

        do nf = 1,size(this,1)
           do nc = 1, size(this(1)%T_Loc,1)
              X = sqrt(Real(this(nf)%T_Loc(nc)*conjg(this(nf)%T_Loc(nc)),kind(0.d0)))
              If ( X  > Xmax_loc)   Xmax_loc =  X
           enddo
           do nc = 1,this(1)%N_bonds
              X = sqrt( Real(this(nf)%T(nc)*conjg( this(nf)%T(nc)),kind(0.d0) )  )
              If ( X  > Xmax_hop )   Xmax_hop =  X
           enddo
        enddo

        If (     Xmax_loc < Zero  .and.  Xmax_hop < Zero )  then
           inquire_hop = 0     !  Zero
        elseif ( Xmax_loc > Zero  .and.  Xmax_hop < Zero )  then
           inquire_hop = 1     !  Diagonal
        else
           inquire_hop = 2     !  Full
        endif

      end function inquire_hop

!--------------------------------------------------------------------
!> @author
!> Francesco Parisen Toldin
!>
!> @brief
!> Check the consistency of the checkerboard decomposition.
!> The following tests are performed:
!>   - The allocated second dimension of \c List_Fam must be at least as large as the
!>     maximum family size (a smaller allocation is a fatal error; a larger one issues a warning).
!>   - Every bond \c (unit_cell, bond_index) must appear exactly once across all families.
!>   - Within each family, no two bonds may share a site (commutativity check).
!>
!> @param [in]  this
!> \verbatim
!>   Type(Hopping_Matrix_type). The hopping object whose checkerboard decomposition
!>   (N_FAM, L_Fam, List_Fam, List) is to be verified.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry; provides N (number of unit cells) and nnlist.
!> \endverbatim
!> @param [in]  inv_list
!> \verbatim
!>   Integer(:,:). Inverse site map: inv_list(unit_cell, orbital) = global site index.
!> \endverbatim
!> @return  .true. if all tests pass; .false. (with messages to error_unit) otherwise.
!
!--------------------------------------------------------------------
      Logical Function test_checkerboard_decomposition(this, Latt, inv_list)
        Implicit none

        Type(Hopping_Matrix_type), intent(IN) :: this
        Type(Lattice), intent(IN)             :: Latt
        Integer, intent(IN)                   :: inv_list(:, :)

        ! Local variables
        Logical, allocatable :: all_bonds(:, :)
        Integer              :: maxl, i, j, n1, n2, unit1, bond1, site1a, site1b, unit2, bond2, site2a, site2b

        test_checkerboard_decomposition = .true.
        allocate(all_bonds(Latt%N, this%N_bonds))
        all_bonds = .false.

        ! Check size of families
        maxl = this%L_Fam(1)
        do i = 2, this%N_Fam
           if (maxl < this%L_Fam(i)) maxl = this%L_Fam(i)
        end do
        if (maxl > size(this%List_Fam, 2)) then
           write(error_unit, *) 'Error in the length of families. Maximum length found is ', maxl, ' allocated size is ', size(this%List_Fam, 2)
           test_checkerboard_decomposition = .false.
        else if (maxl < size(this%List_Fam, 2)) then
           write(error_unit, *) 'Warning: the maximum family length is ', maxl, ' allocated size is ', size(this%List_Fam, 2)
        end if
        
        ! Check duplicates
        do i = 1, this%N_Fam
           do j = 1, this%L_Fam(i)
              if (all_bonds(this%List_Fam(i, j, 1), this%List_Fam(i, j, 2))) then
                 write(error_unit, *) 'Error in decomposition: bond at List_Fam(', i, ' ', j, ') is present twice'
                 test_checkerboard_decomposition = .false.
              else
                 all_bonds(this%List_Fam(i, j, 1), this%List_Fam(i, j, 2)) = .true.
              end if
           end do
        end do

        ! Check that all bonds are present in the decomposition
        do i = 1, Latt%N
           do j = 1, this%N_bonds
              if (.not.(all_bonds(i, j))) then
                 write(error_unit, *) 'Error: bonds at Nunit_cell = ', i, ' bond no. ', j, ' is missing'
                 test_checkerboard_decomposition = .false.
              end if
           end do
        end do

        ! Check commutativity
        do i = 1, this%N_Fam
           do n1 = 1, this%L_Fam(i) - 1
              ! Sites of the first bond
              unit1 = this%List_Fam(i, n1, 1)
              bond1 = this%List_Fam(i, n1, 2)
              site1a = inv_list(unit1, this%List(bond1, 1))
              site1b = inv_list(Latt%nnlist(unit1, this%List(bond1, 3), this%List(bond1, 4)), this%List(bond1, 2))
              do n2 = n1 + 1, this%L_Fam(i)
                 unit2 = this%List_Fam(i, n2, 1)
                 bond2 = this%List_Fam(i, n2, 2)
                 site2a = inv_list(unit2, this%List(bond2, 1))
                 site2b = inv_list(Latt%nnlist(unit2, this%List(bond2, 3), this%List(bond2, 4)), this%List(bond2, 2))

                 if ((site1a == site2a) .or. (site1a == site2b) .or. (site1b == site2a) .or. (site1b == site2b)) then
                    write(error_unit, *) 'Error: non-communting hoppings at family ', i, ' n1 = ', n1, ' List_Fam(i, n1) = ', &
                         & this%List_Fam(i, n1, 1), ' ', this%List_Fam(i, n1, 2), ' site1a = ', site1a, ' site1b = ', site1b, &
                         & ' n2 = ', n2, ' List_Fam(i, n2) = ', &
                         & this%List_Fam(i, n2, 1), ' ', this%List_Fam(i, n2, 2), ' site2a = ', site2a, ' site2b = ', site2b
                    test_checkerboard_decomposition = .false.
                 end if
              end do
           end do
        end do
        
        deallocate(all_bonds)
      end Function test_checkerboard_decomposition
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Set the default hopping matrix for the square lattice.
!> \c Ham_T is the nearest-neighbour hopping amplitude and \c Ham_Chem the on-site
!> chemical potential.  The checkerboard decomposition uses 4 families:
!> sites with \c mod(List(I,1)+List(I,2),2)==0 (sublattice A) contribute bonds 1 and 2
!> to families 1 and 2; the remaining sites (sublattice B) contribute to families 3 and 4.
!> If \c L_2==1 the routine delegates to Set_Default_hopping_parameters_N_Leg_Ladder.
!>
!> @param [out]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!>   Allocated and filled by this routine.
!> \endverbatim
!> @param [in]  Ham_T_vec
!> \verbatim
!>   Real(:). Nearest-neighbour hopping amplitude for each flavor.
!> \endverbatim
!> @param [in]  Ham_Chem_vec
!> \verbatim
!>   Real(:). Chemical potential for each flavor.
!> \endverbatim
!> @param [in]  Phi_X_vec, Phi_Y_vec
!> \verbatim
!>   Real(:). Twist boundary conditions / Aharonov-Bohm flux along a1 and a2 for each flavor.
!> \endverbatim
!> @param [in]  Bulk
!> \verbatim
!>   Logical. If .true. the twist is applied as a gauge field in the bulk;
!>   if .false. it is applied only at the boundary.
!> \endverbatim
!> @param [in]  N_Phi_vec
!> \verbatim
!>   Integer(:). Number of magnetic flux quanta threading the lattice for each flavor.
!> \endverbatim
!> @param [in]  N_FL
!> \verbatim
!>   Integer. Number of fermion flavors.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!>   List(I,1) = unit cell, List(I,2) = orbital; Invlist(uc,orb) = site index.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information (Norb, orbital positions).
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_square(this, Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
           &                                           List, Invlist, Latt, Latt_unit )

        Implicit none

        Type  (Hopping_Matrix_type), allocatable     :: this(:)
        Real (Kind=Kind(0.d0)), Intent(IN),Dimension(:)   :: Ham_T_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
        Integer, Intent(IN),Dimension(:)                  :: N_Phi_vec
        Integer, Intent(IN)                               :: N_FL
        Logical, Intent(IN)                               :: Bulk
        Integer, Intent(IN), Dimension(:,:)               :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit


        ! Local
        Integer :: nf,N_Bonds, nc, I, I1
        Real (Kind = Kind(0.d0) ) :: Zero = 1.0E-8,  Ham_T_max
        Real (Kind = Kind(0.d0) ), allocatable :: Ham_T_perp_vec(:)

        
        
        If ( Xnorm(Latt%L2_p - Latt%a2_p)  < Zero )  then
           Allocate( Ham_T_perp_vec(N_FL) )
           Ham_T_perp_vec = 0.d0
           Call Set_Default_hopping_parameters_N_Leg_Ladder(this,Ham_T_vec, Ham_T_perp_vec, Ham_Chem_vec, Phi_X_vec, &
                &                                           Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
                &                                           List, Invlist, Latt, Latt_unit )
           Deallocate ( Ham_T_perp_vec )
        else
           If (  mod(nint(latt%L1_p(1)),2)  /=  0   .or.   mod(nint(latt%L2_p(2)),2)  /= 0 )  then
              Write(error_unit,*) '*** For  the  square  lattice,  our  implementation of the checkerborad '
              Write(error_unit,*) 'decomposition  requires even  values of L_1  and L_2  ***'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           Allocate( this(N_FL) )

           Ham_T_max = 0.d0
           Do nf = 1,N_FL
              If ( Abs(Ham_T_vec(nf))   >  Ham_T_max )  Ham_T_max = Abs(Ham_T_vec(nf))
           Enddo

           ! Bond-list convention (applies to all lattices):
           !   List(nc, 1) = no_from  : orbital index of the source site in the unit cell
           !   List(nc, 2) = no_to    : orbital index of the target site in the unit cell
           !   List(nc, 3) = n1       : unit-cell offset along lattice vector a1
           !   List(nc, 4) = n2       : unit-cell offset along lattice vector a2
           ! The encoded hopping matrix element is:
           !   H[ (i, no_from), (i + n1*a1 + n2*a2, no_to) ] = T(nc)  where T(nc) = -Ham_T
           do nf = 1,N_FL
              this(nf)%N_bonds = 0
              if ( abs(Ham_T_max) > Zero)  then
                 this(nf)%N_bonds = 2
                 Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                      &    this(nf)%T(this(nf)%N_bonds) )
                 nc = 0
                 ! Bond 1: orb1 -> orb1, offset (n1,n2) = (0,+1)  =>  NN hop along a2 (y-direction)
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 1
                 this(nf)%List(nc,2) = 1
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 1

                 ! Bond 2: orb1 -> orb1, offset (n1,n2) = (+1,0)  =>  NN hop along a1 (x-direction)
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 1
                 this(nf)%List(nc,2) = 1
                 this(nf)%List(nc,3) = 1
                 this(nf)%List(nc,4) = 0
              Endif
              Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
              do nc = 1,Latt_Unit%Norb
                 this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
              enddo
              this(nf)%N_Phi =  N_Phi_vec(nf)
              this(nf)%Phi_X =  Phi_X_vec(nf)
              this(nf)%Phi_Y =  Phi_Y_vec(nf)
              this(nf)%Bulk  =  Bulk
           enddo

           !Set Checkerboard
           ! The two bond types are split across 4 families using the checkerboard sublattice
           ! parity mod(n1+n2, 2) of each unit cell, so that all bonds within a family connect
           ! disjoint site pairs and their single-bond propagators exp(-dtau*h_bond) commute.
           !   Sublattice A  (mod(n1+n2,2)==0):  Family 1 -> bond 1 (a2-hop)
           !                                     Family 2 -> bond 2 (a1-hop)
           !   Sublattice B  (mod(n1+n2,2)==1):  Family 3 -> bond 1 (a2-hop)
           !                                     Family 4 -> bond 2 (a1-hop)
           if ( Ham_T_max   > Zero ) then
              this(1)%N_FAM  = 4
              Allocate (this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM))
              this(1)%L_FAM  = Latt%N/2
              this(1)%Prop_Fam= 1.d0
              Allocate (this(1)%List_Fam(this(1)%N_FAM,this(1)%L_Fam(1),2))
              this(1)%L_FAM  = 0
              do I = 1,Latt%N
                 if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                    Nf = 1
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1 ! The bond (See above)
                    Nf = 2
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2
                 else
                    Nf = 3
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1
                    Nf = 4
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2
                 endif
              enddo
           endif
        endif
      end Subroutine Set_Default_hopping_parameters_square

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Set the default hopping matrix for the triangular lattice.
!> \c Ham_T is the nearest-neighbour hopping amplitude (three bond directions) and
!> \c Ham_Chem the on-site chemical potential.
!> The checkerboard decomposition uses 6 families:
!> bonds 1 and 2 (along a1 and a2) are split by the parity of
!> \c mod(List(I,1)+List(I,2),2): even sites -> families 1 and 2; odd sites -> families 4 and 5.
!> Bond 3 (diagonal, along a2-a1) is split by \c mod(List(I,1),2):
!> even -> family 3; odd -> family 6.
!>
!> @param [out]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!> \endverbatim
!> @param [in]  Ham_T_vec
!> \verbatim
!>   Real(:). Nearest-neighbour hopping amplitude for each flavor.
!> \endverbatim
!> @param [in]  Ham_Chem_vec
!> \verbatim
!>   Real(:). Chemical potential for each flavor.
!> \endverbatim
!> @param [in]  Phi_X_vec, Phi_Y_vec
!> \verbatim
!>   Real(:). Twist boundary conditions / flux along a1 and a2 for each flavor.
!> \endverbatim
!> @param [in]  Bulk
!> \verbatim
!>   Logical. Twist applied in bulk (.true.) or at boundary (.false.).
!> \endverbatim
!> @param [in]  N_Phi_vec
!> \verbatim
!>   Integer(:). Number of magnetic flux quanta for each flavor.
!> \endverbatim
!> @param [in]  N_FL
!> \verbatim
!>   Integer. Number of fermion flavors.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information.
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_triangular(this, Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
         &                                           List, Invlist, Latt, Latt_unit )

      Implicit none

      Type  (Hopping_Matrix_type), allocatable     :: this(:)
      Real (Kind=Kind(0.d0)), Intent(IN),Dimension(:)   :: Ham_T_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      Integer, Intent(IN),Dimension(:)                  :: N_Phi_vec
      Integer, Intent(IN)                               :: N_FL
      Logical, Intent(IN)                               :: Bulk
      Integer, Intent(IN), Dimension(:,:)               :: List, Invlist
      Type(Lattice),  Intent(in)            :: Latt
      Type(Unit_cell),Intent(in)            :: Latt_unit


      ! Local
      Integer :: nf,N_Bonds, nc, I, I1
      Real (Kind = Kind(0.d0) ) :: Zero = 1.0E-8,  Ham_T_max, x_p(2)

         !Write(6,*)   Iscalar(Latt%L1_p,Latt%BZ1_p)/(2.d0*acos(-1.d0)),Iscalar(Latt%L2_p,Latt%BZ2_p)/(2.d0*acos(-1.d0)) 
      
         If (  mod(nint(Iscalar(Latt%L1_p,Latt%BZ1_p)/(2.d0*acos(-1.d0))),2)  /=  0   .or. & 
             & mod(nint(Iscalar(Latt%L2_p,Latt%BZ2_p)/(2.d0*acos(-1.d0))),2)  /= 0 )  then
            Write(error_unit,*) '*** For  the  triangular  lattice,  our  implementation of the checkerborad '
            Write(error_unit,*) 'decomposition  requires even  values of L_1  and L_2  ***'
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif
         Allocate( this(N_FL) )

         Ham_T_max = 0.d0
         Do nf = 1,N_FL
            If ( Abs(Ham_T_vec(nf))   >  Ham_T_max )  Ham_T_max = Abs(Ham_T_vec(nf))
         Enddo

         do nf = 1,N_FL
            this(nf)%N_bonds = 0
            if ( abs(Ham_T_max) > Zero)  then
               this(nf)%N_bonds = 3
               Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                    &    this(nf)%T(this(nf)%N_bonds) )
               nc = 0
               ! Bond 1: orb1 -> orb1, offset (n1,n2) = (+1, 0)  =>  NN along a1
               nc = nc + 1
               this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
               this(nf)%List(nc,1) = 1
               this(nf)%List(nc,2) = 1
               this(nf)%List(nc,3) = 1
               this(nf)%List(nc,4) = 0

               ! Bond 2: orb1 -> orb1, offset (n1,n2) = ( 0,+1)  =>  NN along a2
               nc = nc + 1
               this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
               this(nf)%List(nc,1) = 1
               this(nf)%List(nc,2) = 1
               this(nf)%List(nc,3) = 0
               this(nf)%List(nc,4) = 1

               ! Bond 3: orb1 -> orb1, offset (n1,n2) = (-1,+1)  =>  NN along a2-a1 (third triangular direction)
               nc = nc + 1
               this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
               this(nf)%List(nc,1) = 1
               this(nf)%List(nc,2) = 1
               this(nf)%List(nc,3) =-1
               this(nf)%List(nc,4) = 1
            Endif
            Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
            do nc = 1,Latt_Unit%Norb
               this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
            enddo
            this(nf)%N_Phi =  N_Phi_vec(nf)
            this(nf)%Phi_X =  Phi_X_vec(nf)
            this(nf)%Phi_Y =  Phi_Y_vec(nf)
            this(nf)%Bulk  =  Bulk
         enddo

         ! Do I = 1,Latt%N 
         !    write(10,*)  Latt%List(I,1)*Latt%a1_p  +Latt%List(I,2)*Latt%a2_p  
         !    if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
         !       write(11,*)  Latt%List(I,1)*Latt%a1_p  +Latt%List(I,2)*Latt%a2_p  
         !    endif
         !    if ( mod(Latt%List(I,1),2) == 0 ) then
         !       write(12,*)  Latt%List(I,1)*Latt%a1_p  +Latt%List(I,2)*Latt%a2_p  
         !    endif
         ! enddo

         !Set Checkerboard
         ! The triangular lattice has 3 bond directions requiring 6 families total.
         !
         ! Bonds 1 (a1-hop, List(nc,3)=+1,4=0) and 2 (a2-hop, List(nc,3)=0,4=+1) connect
         ! sites of OPPOSITE checkerboard parity mod(n1+n2,2): the source is at (n1,n2)
         ! and the target at (n1+1,n2) or (n1,n2+1), both with opposite parity.  Splitting
         ! by mod(n1+n2,2) therefore guarantees disjoint site pairs within each family:
         !   Sublattice A (mod(n1+n2,2)==0):  Family 1 -> bond 1 (a1-hop)
         !                                    Family 2 -> bond 2 (a2-hop)
         !   Sublattice B (mod(n1+n2,2)==1):  Family 4 -> bond 1 (a1-hop)
         !                                    Family 5 -> bond 2 (a2-hop)
         !
         ! Bond 3 (diagonal, List(nc,3)=-1,4=+1) connects (n1,n2) -> (n1-1,n2+1).
         ! Its target has parity mod((n1-1)+(n2+1),2) = mod(n1+n2,2): the SAME as the
         ! source.  Splitting by mod(n1+n2,2) would therefore put both endpoints of
         ! different bond-3 instances in the same sublattice, causing site sharing
         ! (e.g., bond-3 at (1,1) targets (0,2), which is the source of bond-3 at (0,2)).
         ! Instead, bond 3 is split by mod(n1,2): the source has n1 and the target n1-1,
         ! so sources and targets always have opposite n1-parity, ensuring disjoint pairs:
         !   Even n1 (mod(n1,2)==0):          Family 3 -> bond 3
         !   Odd  n1 (mod(n1,2)==1):          Family 6 -> bond 3
         if ( Ham_T_max   > Zero ) then
            this(1)%N_FAM  = 6
            Allocate (this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM))
            this(1)%L_FAM  = Latt%N/2
            this(1)%Prop_Fam= 1.d0
            Allocate (this(1)%List_Fam(this(1)%N_FAM,this(1)%L_Fam(1),2))
            this(1)%L_FAM  = 0
            do I = 1,Latt%N
               if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                  Nf = 1
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1 ! The bond (See above)
                  Nf = 2
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2
               else
                  Nf = 4
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1
                  Nf = 5
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2
               endif
               if ( mod(Latt%List(I,1), 2) == 0 ) then
                  Nf = 3
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 3
               else
                  Nf = 6
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                  this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 3
               endif 
            enddo
         endif

    end Subroutine Set_Default_hopping_parameters_triangular
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Set the default hopping matrix for the kagome lattice.
!> \c Ham_T is the nearest-neighbour hopping amplitude (six bond types,
!> three orbitals per unit cell) and \c Ham_Chem the on-site chemical potential.
!> The checkerboard decomposition uses 6 families that map one-to-one onto the
!> 6 bond types: every unit cell contributes its bond \c Nf to family \c Nf,
!> so all bonds within a family connect disjoint site pairs and hence commute.
!>
!> @param [out]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!> \endverbatim
!> @param [in]  Ham_T_vec
!> \verbatim
!>   Real(:). Nearest-neighbour hopping amplitude for each flavor.
!> \endverbatim
!> @param [in]  Ham_Chem_vec
!> \verbatim
!>   Real(:). Chemical potential for each flavor.
!> \endverbatim
!> @param [in]  Phi_X_vec, Phi_Y_vec
!> \verbatim
!>   Real(:). Twist boundary conditions / flux along a1 and a2 for each flavor.
!> \endverbatim
!> @param [in]  Bulk
!> \verbatim
!>   Logical. Twist applied in bulk (.true.) or at boundary (.false.).
!> \endverbatim
!> @param [in]  N_Phi_vec
!> \verbatim
!>   Integer(:). Number of magnetic flux quanta for each flavor.
!> \endverbatim
!> @param [in]  N_FL
!> \verbatim
!>   Integer. Number of fermion flavors.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information.
!> \endverbatim
!
!--------------------------------------------------------------------
    Subroutine Set_Default_hopping_parameters_kagome(this, Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
      &                                           List, Invlist, Latt, Latt_unit )

   Implicit none

   Type  (Hopping_Matrix_type), allocatable     :: this(:)
   Real (Kind=Kind(0.d0)), Intent(IN),Dimension(:)   :: Ham_T_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
   Integer, Intent(IN),Dimension(:)                  :: N_Phi_vec
   Integer, Intent(IN)                               :: N_FL
   Logical, Intent(IN)                               :: Bulk
   Integer, Intent(IN), Dimension(:,:)               :: List, Invlist
   Type(Lattice),  Intent(in)            :: Latt
   Type(Unit_cell),Intent(in)            :: Latt_unit


   ! Local
   Integer :: nf,N_Bonds, nc, I, I1
   Real (Kind = Kind(0.d0) ) :: Zero = 1.0E-8,  Ham_T_max, x_p(2)

      !Write(6,*)   Iscalar(Latt%L1_p,Latt%BZ1_p)/(2.d0*acos(-1.d0)),Iscalar(Latt%L2_p,Latt%BZ2_p)/(2.d0*acos(-1.d0)) 

      Allocate( this(N_FL) )

      Ham_T_max = 0.d0
      Do nf = 1,N_FL
         If ( Abs(Ham_T_vec(nf))   >  Ham_T_max )  Ham_T_max = Abs(Ham_T_vec(nf))
      Enddo

      ! List(N_b,1) = no_1
      ! List(N_b,2) = no_2
      ! List(N_b,3) = n_1
      ! List(N_b,4) = n_2
      ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)   
      !
      ! Kagome unit cell: 3 orbitals per cell (orb1, orb2, orb3 forming a corner-sharing triangle).
      ! The 6 bond types cover all nearest-neighbour pairs; 3 are intra-cell (n1=n2=0)
      ! and 3 connect to neighbouring cells.
      do nf = 1,N_FL
         this(nf)%N_bonds = 0
         if ( abs(Ham_T_max) > Zero)  then
            this(nf)%N_bonds = 6
            Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                 &    this(nf)%T(this(nf)%N_bonds) )
            nc = 0
            ! Bond 1: orb1 -> orb2, (0, 0)  intra-cell (short bond along the a1-base)
            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 1
            this(nf)%List(nc,2) = 2
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 0

            ! Bond 2: orb1 -> orb3, (0, 0)  intra-cell
            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 1
            this(nf)%List(nc,2) = 3
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 0

            ! Bond 3: orb2 -> orb3, (0, 0)  intra-cell
            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 2
            this(nf)%List(nc,2) = 3
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 0

            ! Bond 4: orb3 -> orb1, ( 0,+1)  inter-cell, +a2
            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 3
            this(nf)%List(nc,2) = 1
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 1

            ! Bond 5: orb3 -> orb2, (-1,+1)  inter-cell, along a2-a1
            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 3
            this(nf)%List(nc,2) = 2
            this(nf)%List(nc,3) = -1
            this(nf)%List(nc,4) = 1

            ! Bond 6: orb1 -> orb2, (-1, 0)  inter-cell, along -a1
            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 1
            this(nf)%List(nc,2) = 2
            this(nf)%List(nc,3) = -1
            this(nf)%List(nc,4) = 0
         Endif
         Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
         do nc = 1,Latt_Unit%Norb
            this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
         enddo
         this(nf)%N_Phi =  N_Phi_vec(nf)
         this(nf)%Phi_X =  Phi_X_vec(nf)
         this(nf)%Phi_Y =  Phi_Y_vec(nf)
         this(nf)%Bulk  =  Bulk
      enddo

      !Set Checkerboard
      ! The kagome lattice uses one family per bond type (6 families total).
      ! Every unit cell contributes exactly one bond of each type; because the two orbital
      ! endpoints of any given bond type are distinct for every unit cell and are not shared
      ! with the same bond type from any other unit cell, all bonds in a family connect
      ! disjoint site pairs and their propagators commute trivially.
      ! Note: bonds 1 and 6 both connect orb1->orb2 but with different cell offsets
      ! ((0,0) vs (-1,0)); they must be in separate families because they share the orb1
      ! site at the source unit cell.
      !   Family 1 (bond 1): orb1 -> orb2, ( 0, 0)  intra-cell, midpoint of a1-edge
      !   Family 2 (bond 2): orb1 -> orb3, ( 0, 0)  intra-cell, midpoint of a2-edge
      !   Family 3 (bond 3): orb2 -> orb3, ( 0, 0)  intra-cell, closes the triangle
      !   Family 4 (bond 4): orb3 -> orb1, ( 0,+1)  inter-cell, bridge along +a2
      !   Family 5 (bond 5): orb3 -> orb2, (-1,+1)  inter-cell, bridge along a2-a1
      !   Family 6 (bond 6): orb1 -> orb2, (-1, 0)  inter-cell, bridge along -a1
      if ( Ham_T_max   > Zero ) then
         this(1)%N_FAM  = 6
         Allocate (this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM))
         this(1)%L_FAM  = Latt%N
         this(1)%Prop_Fam= 1.d0
         Allocate (this(1)%List_Fam(this(1)%N_FAM,this(1)%L_Fam(1),2))
         this(1)%L_FAM  = 0
         do I = 1,Latt%N
            Do Nf = 1,this(1)%N_FAM 
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I  ! Unit cell
               this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = Nf ! The bond (See above)
            enddo
         enddo
      endif

 end Subroutine Set_Default_hopping_parameters_kagome
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Set the default hopping matrix for an N-leg ladder.
!> \c Ham_T is the hopping along each chain (nearest neighbours along a1),
!> \c Ham_T_perp is the inter-rung hopping (nearest neighbours connecting orbital
!> \c n to \c n+1 within the same unit cell), and \c Ham_Chem is the on-site
!> chemical potential.  When \c Latt%N==1 the routine builds an open-boundary
!> one-dimensional chain; otherwise periodic boundary conditions along a1 are used.
!> The checkerboard decomposition is based on the parity of the unit-cell index
!> along a1 (families 1/2) and the parity of the rung bond index (families 3/4).
!>
!> @param [out]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!> \endverbatim
!> @param [in]  Ham_T_vec
!> \verbatim
!>   Real(:). Hopping amplitude along each chain for each flavor.
!> \endverbatim
!> @param [in]  Ham_T_perp_vec
!> \verbatim
!>   Real(:). Inter-rung hopping amplitude for each flavor.
!> \endverbatim
!> @param [in]  Ham_Chem_vec
!> \verbatim
!>   Real(:). Chemical potential for each flavor.
!> \endverbatim
!> @param [in]  Phi_X_vec, Phi_Y_vec
!> \verbatim
!>   Real(:). Twist boundary conditions / flux along a1 and a2 for each flavor.
!> \endverbatim
!> @param [in]  Bulk
!> \verbatim
!>   Logical. Twist applied in bulk (.true.) or at boundary (.false.).
!> \endverbatim
!> @param [in]  N_Phi_vec
!> \verbatim
!>   Integer(:). Number of magnetic flux quanta for each flavor.
!> \endverbatim
!> @param [in]  N_FL
!> \verbatim
!>   Integer. Number of fermion flavors.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information (Norb = number of legs).
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_N_Leg_Ladder  &
           &               (this, Ham_T_vec, Ham_T_perp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk, &
           &                N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)


        Implicit none

        type   (Hopping_Matrix_type), allocatable :: this(:)
        Integer, Intent(IN)                   :: N_FL
        Real (Kind=Kind(0.d0)), Intent(IN), Dimension(:)    :: Ham_T_vec, Ham_T_perp_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
        Integer, Intent(IN), Dimension(:)                   :: N_Phi_vec
        Logical, Intent(IN)                                 :: Bulk
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit


        ! Local
        Integer :: nf,N_Bonds, nc, I, I1, n, no
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8

        
        select case (Latt%N)
        case(1)   !  Here the length of the  N_leg_ladder is unity such  that it
                  !  effectivley maps onto a one-dimensional chain with open boundary conditions.
           ! With only one unit cell along a1 there are no leg hops; only the
           ! N_legs-1 intra-cell rung bonds between adjacent legs are included.
           
           
           Allocate( this(N_FL) )
           do nf = 1,N_FL
              this(nf)%N_bonds =  Latt_unit%Norb - 1
              Allocate (this(nf)%List( this(nf)%N_bonds,4 ), this(nf)%T( this(nf)%N_bonds ) )
              nc = 0
              do n = 1,Latt_unit%Norb  -1
                 ! Bond n: orb_n -> orb_{n+1}, (0,0)  rung hop between legs n and n+1
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T_perp_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = n
                 this(nf)%List(nc,2) = n + 1
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 0
              enddo
              
              Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
              do nc = 1,Latt_Unit%Norb
                 this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
              enddo
              this(nf)%N_Phi =  N_Phi_vec(nf)
              this(nf)%Phi_X =  Phi_X_vec(nf)
              this(nf)%Phi_Y =  Phi_Y_vec(nf)
              this(nf)%Bulk =   Bulk
           enddo
           
           ! Set Checkerboard (single unit cell / open chain)
           ! Rung bonds are split into two families by the parity of the bond index n,
           ! so bonds within a family share no sites and their propagators commute:
           !   Family 1: odd-indexed rung bonds  (n = 1, 3, 5, ...)
           !   Family 2: even-indexed rung bonds (n = 2, 4, 6, ...)
           If  ( Latt_Unit%Norb  <=  2 ) then
              this(1)%N_FAM        = 1
           else
              this(1)%N_FAM        = 2
           endif
           Allocate ( this(1)%L_Fam( this(1)%N_FAM ),  this(1)%Prop_Fam( this(1)%N_FAM ) )
           this(1)%L_Fam    = Latt_unit%Norb/2
           this(1)%Prop_Fam = 1.d0
           Allocate ( this(1)%List_Fam(this(1)%N_FAM,this(1)%L_Fam(1),2) )
           
           
           this(1)%L_FAM  = 0
           do no = 1,Latt_unit%Norb - 1
              if (mod(no,2) == 1 ) then
                 Nf = 1
                 !Write(6,*)  NF, no 
                 do I = 1,Latt%N
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no 
                 enddo
              else
                 Nf = 2
                 !Write(6,*)  NF, no 
                 do I = 1,Latt%N
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no 
                 enddo
              endif
           enddo

           
        case default
           If (  mod(nint(latt%L1_p(1)),2)  /=  0  )  then
              Write(error_unit,*) '*** For  the N_leg_ladder  lattice,  our  implementation of the checkerborad '
              Write(error_unit,*) 'decomposition  requires L_1 = 1 or  L_1   even ***'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif

           !Write(6,*) Ham_T_vec,  Ham_T_perp_vec, Ham_chem_vec
           Allocate( this(N_FL) )
           do nf = 1,N_FL
              ! Bonds 1..N_legs: leg hops along a1 for each orbital (leg)
              !   Bond n: orb_n -> orb_n, (1,0)
              ! Bonds N_legs+1..2*N_legs-1: rung hops within the unit cell
              !   Bond N_legs+n: orb_n -> orb_{n+1}, (0,0)
              this(nf)%N_bonds = Latt_unit%Norb +  (Latt_unit%Norb - 1 )
              Allocate (this(nf)%List( this(nf)%N_bonds,4 ), &
                   &    this(nf)%T( this(nf)%N_bonds ) )
              nc = 0
              do n = 1,Latt_unit%Norb
                 ! Leg hop for orbital (leg) n along a1
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = n
                 this(nf)%List(nc,2) = n
                 this(nf)%List(nc,3) = 1
                 this(nf)%List(nc,4) = 0
              enddo
              
              do n = 1,Latt_unit%Norb -1
                 ! Rung hop between legs n and n+1 within the same unit cell
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T_perp_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = n
                 this(nf)%List(nc,2) = n + 1
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 0
              enddo
              
              Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
              do nc = 1,Latt_Unit%Norb
                 this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
              enddo
              this(nf)%N_Phi =  N_Phi_vec(nf)
              this(nf)%Phi_X =  Phi_X_vec(nf)
              this(nf)%Phi_Y =  Phi_Y_vec(nf)
              this(nf)%Bulk =   Bulk
           enddo
           
           ! Write(6,*) Latt_unit%Norb
           ! Set Checkerboard (full N-leg ladder)
           ! Leg bonds (1..N_legs) are split into families 1 and 2 by the parity of
           ! the unit-cell index along a1: mod(n1,2)==0 -> family 1, else -> family 2.
           ! Rung bonds (N_legs+1..2*N_legs-1) are split into families 3 and 4 by the
           ! parity of the rung bond index n: odd n -> family 3, even n -> family 4.
           If     ( Latt_Unit%Norb  == 1 ) then
              this(1)%N_FAM        = 2    ! Only leg bonds; no rung bonds
           elseif ( Latt_Unit%Norb  == 2 ) then
              this(1)%N_FAM        = 3    ! One rung bond -> needs only one extra family
           else
              this(1)%N_FAM        = 4
           endif
           Allocate ( this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM) )
           this(1)%L_Fam    = Latt%N*Latt_unit%Norb/2
           this(1)%Prop_Fam= 1.d0
           Allocate ( this(1)%List_Fam(this(1)%N_FAM,this(1)%L_Fam(1),2) )
           
           
           this(1)%L_FAM  = 0
           do I = 1,Latt%N
              if ( mod(Latt%List(I,1),2) == 0 ) then
                 Nf = 1   ! Even unit-cell index along a1: leg bonds go to family 1
                 do no = 1,Latt_unit%Norb
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no ! The bond (See above)
                 enddo
              else
                 Nf = 2   ! Odd unit-cell index along a1: leg bonds go to family 2
                 do no = 1,Latt_unit%Norb
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no
                 enddo
              endif
           enddo
           do no = 1,Latt_unit%Norb - 1
              if (mod(no,2) == 1 ) then
                 Nf = 3   ! Odd rung index: rung bond no+N_legs goes to family 3
                 !Write(6,*)  NF, no + Latt_unit%Norb
                 do I = 1,Latt%N
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no + Latt_unit%Norb
                 enddo
              else
                 Nf = 4   ! Even rung index: rung bond no+N_legs goes to family 4
                 !Write(6,*)  NF, no + Latt_unit%Norb
                 do I = 1,Latt%N
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no + Latt_unit%Norb
                 enddo
              endif
           enddo
        end select
        
      end Subroutine Set_Default_hopping_parameters_N_Leg_Ladder



      
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Set the default hopping matrix for the honeycomb lattice.
!> \c Ham_T is the nearest-neighbour hopping amplitude (three bond types connecting
!> the two sublattice orbitals), \c Lambda is the Kane-Mele spin-orbit coupling
!> (**not yet implemented** — a non-zero value triggers a runtime error), and
!> \c Ham_Chem is the on-site chemical potential.
!> The checkerboard decomposition uses 3 families that map one-to-one onto the 3
!> nearest-neighbour bond types; each family contains all \c Latt%N bonds of that type.
!>
!> @param [out]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!> \endverbatim
!> @param [in]  Ham_T_vec
!> \verbatim
!>   Real(:). Nearest-neighbour hopping amplitude for each flavor.
!> \endverbatim
!> @param [in]  Ham_Lambda_vec
!> \verbatim
!>   Real(:). Kane-Mele spin-orbit coupling for each flavor.
!>   Currently not implemented; must be zero.
!> \endverbatim
!> @param [in]  Ham_Chem_vec
!> \verbatim
!>   Real(:). Chemical potential for each flavor.
!> \endverbatim
!> @param [in]  Phi_X_vec, Phi_Y_vec
!> \verbatim
!>   Real(:). Twist boundary conditions / flux along a1 and a2 for each flavor.
!> \endverbatim
!> @param [in]  Bulk
!> \verbatim
!>   Logical. Twist applied in bulk (.true.) or at boundary (.false.).
!> \endverbatim
!> @param [in]  N_Phi_vec
!> \verbatim
!>   Integer(:). Number of magnetic flux quanta for each flavor.
!> \endverbatim
!> @param [in]  N_FL
!> \verbatim
!>   Integer. Number of fermion flavors.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information.
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_honeycomb(this,Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, &
           &                                              Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
           &                                              List, Invlist, Latt, Latt_unit)

        Implicit none

        type (Hopping_Matrix_type), allocatable            :: this(:)
        Real (Kind=Kind(0.d0)), Intent(IN),Dimension(:), allocatable  :: Ham_T_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec, Ham_Lambda_vec
        Integer, Intent(IN),Dimension(:), allocatable                 :: N_Phi_vec
        Integer, Intent(IN)                                           :: N_FL
        Logical, Intent(IN)                                           :: Bulk
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit

        ! Local
        Integer :: nf,N_Bonds, nc, I, I1, n, no
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8, Ham_Lambda_Max

        !Write(6,*) Ham_T_vec, Ham_Chem_vec
        Ham_Lambda_Max = 0.d0
        do nf = 1,N_FL
           if ( Abs(Ham_Lambda_vec(nf)) > Ham_Lambda_Max ) Ham_Lambda_Max =  Abs(Ham_Lambda_vec(nf))
        enddo
        If (abs(Ham_Lambda_max) > 0 ) then
           Write(error_unit,*) 'Kane Mele term is not yet implemented'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        Allocate( this(N_FL) )
        do nf = 1,N_FL
           ! Honeycomb unit cell: 2 orbitals per cell (A-sublattice = orb1, B-sublattice = orb2).
           ! The 3 nearest-neighbour bond types cover all A-B connections.
           this(nf)%N_bonds =  3
           Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                &    this(nf)%T(this(nf)%N_bonds) )
           nc = 0
           ! Bond 1: A(orb1) -> B(orb2), ( 0,  0)  delta1 bond (same unit cell)
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  1
           this(nf)%List(nc,2) =  2
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  0

           ! Bond 2: B(orb2) -> A(orb1), ( 0, +1)  delta2 bond (+a2)
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  2
           this(nf)%List(nc,2) =  1
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  1

           ! Bond 3: A(orb1) -> B(orb2), (+1, -1)  delta3 bond (a1-a2)
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  1
           this(nf)%List(nc,2) =  2
           this(nf)%List(nc,3) =  1
           this(nf)%List(nc,4) = -1

           Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
           do nc = 1,Latt_Unit%Norb
              this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
           enddo
           this(nf)%N_Phi =  N_Phi_vec(nf)
           this(nf)%Phi_X =  Phi_X_vec(nf)
           this(nf)%Phi_Y =  Phi_Y_vec(nf)
           this(nf)%Bulk =   Bulk
        enddo

        ! Set Checkerboard
        ! The honeycomb lattice uses one family per bond type (3 families total).
        ! Each family consists of all Latt%N bonds of that type, one per unit cell.
        ! Since each bond type connects A and B sites that are not shared between
        ! different unit cells for the same bond type, the bonds within a family
        ! connect disjoint site pairs and their propagators commute trivially.
        !   Family 1: all bond-1 (delta1, same-cell A->B) bonds
        !   Family 2: all bond-2 (delta2, +a2 B->A) bonds
        !   Family 3: all bond-3 (delta3, a1-a2 A->B) bonds
        this(1)%N_FAM  = 3
        Allocate (this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM))
        this(1)%L_FAM  = Latt%N
        this(1)%Prop_Fam= 1.d0
        Allocate (this(1)%List_Fam(this(1)%N_FAM,this(1)%L_Fam(1),2))
        do I = 1,Latt%N
           Do  nf = 1,this(1)%N_FAM
              this(1)%List_Fam(nf,I,1) = I  ! Unit cell
              this(1)%List_Fam(nf,I,2) = nf ! The bond (See above)
           Enddo
        enddo

      end Subroutine Set_Default_hopping_parameters_honeycomb

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Set the default hopping matrix for the bilayer square lattice.
!> \c Ham_T1 and \c Ham_T2 are the nearest-neighbour hopping amplitudes on layer 1 and
!> layer 2, respectively (orbitals 1 and 2 per unit cell), and \c Ham_Tperp is the
!> interlayer hopping amplitude (vertical rung between orbital 1 and orbital 2).
!> The checkerboard decomposition is based on the parity of \c mod(List(I,1)+List(I,2),2)
!> for the in-plane bonds and adds a separate family for the interlayer bonds.
!>
!> @param [out]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!> \endverbatim
!> @param [in]  Ham_T1_vec
!> \verbatim
!>   Real(:). In-plane hopping on layer 1 (orbital 1) for each flavor.
!> \endverbatim
!> @param [in]  Ham_T2_vec
!> \verbatim
!>   Real(:). In-plane hopping on layer 2 (orbital 2) for each flavor.
!> \endverbatim
!> @param [in]  Ham_Tperp_vec
!> \verbatim
!>   Real(:). Interlayer hopping for each flavor.
!> \endverbatim
!> @param [in]  Ham_Chem_vec
!> \verbatim
!>   Real(:). Chemical potential for each flavor.
!> \endverbatim
!> @param [in]  Phi_X_vec, Phi_Y_vec
!> \verbatim
!>   Real(:). Twist boundary conditions / flux along a1 and a2 for each flavor.
!> \endverbatim
!> @param [in]  Bulk
!> \verbatim
!>   Logical. Twist applied in bulk (.true.) or at boundary (.false.).
!> \endverbatim
!> @param [in]  N_Phi_vec
!> \verbatim
!>   Integer(:). Number of magnetic flux quanta for each flavor.
!> \endverbatim
!> @param [in]  N_FL
!> \verbatim
!>   Integer. Number of fermion flavors.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information (Norb=2 for the two layers).
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_Bilayer_square(this,Ham_T1_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
           &                                                   Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
           &                                                   List, Invlist, Latt, Latt_unit )

        Implicit none

        type  (Hopping_Matrix_type), allocatable            :: this(:)
        Real (Kind=Kind(0.d0)), Intent(IN),dimension(:)  :: Ham_T1_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
        Integer, Intent(IN),dimension(:)                 :: N_Phi_vec
        Logical, Intent(IN)                   :: Bulk
        Integer, Intent(IN)                   :: N_FL
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit

        

        ! Local
        Integer :: nf,N_Bonds, nc, I, I1, No_Shift, n, nb
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8
        Logical :: Test=.false.
        Real (Kind=Kind(0.d0))                :: Ham_T1_max, Ham_T2_max, Ham_Tperp_max

        Ham_T1_max    = 0.d0
        Ham_T2_max    = 0.d0
        Ham_Tperp_max = 0.d0
        do nf = 1,N_FL
           if (abs(Ham_T1_vec   (nf)) > Ham_T1_max    ) Ham_T1_max    = abs(Ham_T1_vec(nf)   )
           if (abs(Ham_T2_vec   (nf)) > Ham_T2_max    ) Ham_T2_max    = abs(Ham_T2_vec(nf)   )
           if (abs(Ham_Tperp_vec(nf)) > Ham_Tperp_max ) Ham_Tperp_max = abs(Ham_Tperp_vec(nf))
        enddo



        If ( nint( Latt%L2_p(2) )   == 1  )  then
           If (  mod(nint(latt%L1_p(1)),2)  /=  0 )  then
              Write(error_unit,*) '*** For  the Bilayer square lattice,  our  implementation of the checkerborad '
              Write(error_unit,*) 'decomposition  requires L_1  to be  even ***'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           Allocate( this(N_FL) )
           do nf = 1,N_FL
              ! Bilayer square, L2==1 (quasi-1D ladder geometry).
              ! orb1 = layer 1, orb2 = layer 2.  Bond counting:
              !   Bond 1 (always):        layer-1 leg hop along a1, orb1->orb1, (1,0)
              !   Bond 2 (if Tperp /= 0): interlayer rung,           orb1->orb2, (0,0)
              !   Bond 3 (if T2    /= 0): layer-2 leg hop along a1,  orb2->orb2, (1,0)
              ! No_Shift=1 when Tperp /= 0: the interlayer bond occupies index 2,
              !   so the layer-2 leg bond is shifted to index 3 (2+No_Shift).
              N_bonds = 0
              N_bonds = N_bonds + 1
              if (abs(Ham_Tperp_max) > Zero )  N_bonds = N_bonds + 1
              if (abs(Ham_T2_max)    > Zero )  N_bonds = N_bonds + 1
              this(nf)%N_bonds = N_bonds
              Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                   &    this(nf)%T(this(nf)%N_bonds) )
              nc = 0
              ! Bond 1: orb1 -> orb1, (+1, 0)  layer-1 leg hop along a1
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) = 1
              this(nf)%List(nc,2) = 1
              this(nf)%List(nc,3) = 1
              this(nf)%List(nc,4) = 0
              
              If (abs(Ham_Tperp_max) > Zero ) Then
                 ! Bond 2: orb1 -> orb2, (0, 0)  interlayer rung hop
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 1
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 0
              endif
              
              If (abs(Ham_T2_max) > Zero ) Then
                 ! Bond 2+No_Shift: orb2 -> orb2, (+1, 0)  layer-2 leg hop along a1
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 2
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 1
                 this(nf)%List(nc,4) = 0
              endif
              
              Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
              do nc = 1,Latt_Unit%Norb
                 this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
              enddo
              If (Abs(Ham_T2_max) < Zero .and. Abs(Ham_Tperp_max) < Zero ) this(nf)%T_Loc(2)  = cmplx(0.0,0.d0,kind(0.d0))
              this(nf)%N_Phi =  N_Phi_vec(nf)
              this(nf)%Phi_X =  Phi_X_vec(nf)
              this(nf)%Phi_Y =  Phi_Y_vec(nf)
              this(nf)%Bulk =   Bulk
           enddo

           ! Set Checkerboard (L2==1 quasi-1D branch)
           ! Leg bonds are split by the checkerboard sublattice parity mod(n1+n2,2)
           ! of the unit cell, keeping layer-1 and layer-2 leg bonds in the same family
           ! so each family has only disjoint site pairs.
           !   Sublattice A (mod(n1+n2,2)==0):  Family 1 -> bond 1 [+ bond 2+No_Shift if T2 active]
           !   Sublattice B (mod(n1+n2,2)==1):  Family 2 -> bond 1 [+ bond 2+No_Shift if T2 active]
           !   Family 3 (if Tperp /= 0):        all interlayer rung bonds (bond 2), one per unit cell
           ! No_Shift=1 when Tperp /= 0: the layer-2 bond index is 2+No_Shift to skip the rung bond.
           this(1)%N_FAM  = 2
           if (abs(Ham_Tperp_max) > Zero )  this(1)%N_FAM=3

           Allocate (this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM))
           this(1)%Prop_Fam= 1.d0
           
           No_Shift = 0
           If (abs(Ham_Tperp_max) > Zero ) No_Shift=1
           If     ( abs(Ham_T2_max)   <  Zero  .and. abs(Ham_Tperp_max) < Zero)    then
              this(1)%L_FAM  = Latt%N/2
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N/2,2))
           elseif ( abs(Ham_T2_max)   <  Zero  .and. abs(Ham_Tperp_max) > Zero)    then
              this(1)%L_FAM    = Latt%N/2
              this(1)%L_Fam(3) = Latt%N
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N,2))
           elseif ( abs(Ham_T2_max)   >  Zero  .and. abs(Ham_Tperp_max) < Zero)    then
              this(1)%L_FAM    = Latt%N
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N,2))
           elseif ( abs(Ham_T2_max)   >  Zero  .and. abs(Ham_Tperp_max) > Zero)    then
              this(1)%L_FAM    = Latt%N
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N,2))
           endif
           this(1)%L_FAM  = 0
           do I = 1,Latt%N
              if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                 Nf = 1
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1 ! The bond (See above)
                 If (Abs(Ham_T2_max) > Zero) then
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2 + No_Shift ! The bond (See above)
                 endif
              else
                 Nf = 2
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1
                 If (Abs(Ham_T2_max) > Zero) then
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2 + No_Shift ! The bond (See above)
                 endif
              endif
              If (Abs(Ham_Tperp_max) > Zero) then
                 Nf = 3
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2
              Endif
           enddo
        Else
           If (  mod(nint(latt%L1_p(1)),2)  /=  0 .or.  mod(nint(latt%L2_p(2)),2)  /=  0  )  then
              Write(error_unit,*) '*** For  the Bilayer square lattice,  our  implementation of the checkerborad '
              Write(error_unit,*) 'decomposition  requires L_1 and  L_2 to be  even ***'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           
           Allocate( this(N_FL) )
           do nf = 1,N_FL
              ! Bilayer square, 2D geometry (L2 > 1).  orb1 = layer 1, orb2 = layer 2.
              ! Bond counting:
              !   Bond 1 (always):        layer-1 hop along a2, orb1->orb1, (0,+1)
              !   Bond 2 (always):        layer-1 hop along a1, orb1->orb1, (+1,0)
              !   Bond 3 (if Tperp /= 0): interlayer rung,       orb1->orb2, (0,0)
              !   Bond 4 (if T2    /= 0): layer-2 hop along a2,  orb2->orb2, (0,+1)
              !   Bond 5 (if T2    /= 0): layer-2 hop along a1,  orb2->orb2, (+1,0)
              ! No_Shift=1 when Tperp /= 0: layer-2 bond indices are shifted by 1
              !   (4+No_Shift and 5+No_Shift would be wrong; here the pattern is 3+No_Shift, 4+No_Shift).
              N_bonds = 0
              N_bonds = N_bonds + 2
              if (abs(Ham_Tperp_max) > Zero )  N_bonds = N_bonds + 1
              if (abs(Ham_T2_max)    > Zero )  N_bonds = N_bonds + 2
              this(nf)%N_bonds = N_bonds
              Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                   &    this(nf)%T(this(nf)%N_bonds) )
              nc = 0
              ! Bond 1: orb1 -> orb1, ( 0,+1)  layer-1 hop along a2
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) = 1
              this(nf)%List(nc,2) = 1
              this(nf)%List(nc,3) = 0
              this(nf)%List(nc,4) = 1
              
              ! Bond 2: orb1 -> orb1, (+1, 0)  layer-1 hop along a1
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) = 1
              this(nf)%List(nc,2) = 1
              this(nf)%List(nc,3) = 1
              this(nf)%List(nc,4) = 0
              
              If (abs(Ham_Tperp_max) > Zero ) Then
                 ! Bond 3: orb1 -> orb2, (0, 0)  interlayer rung hop
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 1
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 0
              endif
              
              If (abs(Ham_T2_max) > Zero ) Then
                 ! Bond 3+No_Shift: orb2 -> orb2, ( 0,+1)  layer-2 hop along a2
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 2
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 1
                 
                 ! Bond 4+No_Shift: orb2 -> orb2, (+1, 0)  layer-2 hop along a1
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 2
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 1
                 this(nf)%List(nc,4) = 0
              endif
              
              
              Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
              do nc = 1,Latt_Unit%Norb
                 this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
              enddo
              If (Abs(Ham_T2_max) < Zero .and. Abs(Ham_Tperp_max) < Zero ) this(nf)%T_Loc(2)  = cmplx(0.0,0.d0,kind(0.d0))
              this(nf)%N_Phi =  N_Phi_vec(nf)
              this(nf)%Phi_X =  Phi_X_vec(nf)
              this(nf)%Phi_Y =  Phi_Y_vec(nf)
              this(nf)%Bulk =   Bulk
           enddo
           
           ! Set Checkerboard (2D bilayer square branch)
           ! Both layer-1 and layer-2 in-plane bonds are split by the checkerboard
           ! sublattice parity mod(n1+n2,2), and layer-1 and layer-2 bonds of the
           ! same bond direction are placed in the same family.  This works because
           ! the two layers are independent in real space and share no lattice sites.
           !   Sublattice A (mod(n1+n2,2)==0):  Family 1 -> bond 1 (a2) [+ T2 a2]
           !                                    Family 2 -> bond 2 (a1) [+ T2 a1]
           !   Sublattice B (mod(n1+n2,2)==1):  Family 3 -> bond 1 (a2) [+ T2 a2]
           !                                    Family 4 -> bond 2 (a1) [+ T2 a1]
           !   Family 5 (if Tperp /= 0):        all interlayer rung bonds, one per unit cell
           ! No_Shift=1 when Tperp /= 0: layer-2 bond indices are 3+No_Shift and 4+No_Shift.
           this(1)%N_FAM  = 4
           if (abs(Ham_Tperp_max) > Zero )  this(1)%N_FAM=5
           
           Allocate (this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM))
           this(1)%Prop_Fam= 1.d0
           
           No_Shift = 0
           If (abs(Ham_Tperp_max) > Zero ) No_Shift=1
           
           If     ( abs(Ham_T2_max)   <  Zero  .and. abs(Ham_Tperp_max) < Zero)    then
              this(1)%L_FAM  = Latt%N/2
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N/2,2))
           elseif ( abs(Ham_T2_max)   <  Zero  .and. abs(Ham_Tperp_max) > Zero)    then
              this(1)%L_FAM    = Latt%N/2
              this(1)%L_Fam(5) = Latt%N
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N,2))
           elseif ( abs(Ham_T2_max)   >  Zero  .and. abs(Ham_Tperp_max) < Zero)    then
              this(1)%L_FAM    = Latt%N
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N,2))
           elseif ( abs(Ham_T2_max)   >  Zero  .and. abs(Ham_Tperp_max) > Zero)    then
              this(1)%L_FAM    = Latt%N
              Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N,2))
              No_Shift     = 1
           endif
           this(1)%L_FAM  = 0
           do I = 1,Latt%N
              if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                 Nf = 1
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1 ! The bond (See above)
                 If (Abs(Ham_T2_max) > Zero) then
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 3 + No_Shift ! The bond (See above)
                 endif
                 Nf = 2
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2
                 If (Abs(Ham_T2_max) > Zero) then
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 4 + No_Shift ! The bond (See above)
                 endif
              else
                 Nf = 3
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 1
                 If (Abs(Ham_T2_max) > Zero) then
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 3 + No_Shift ! The bond (See above)
                 endif
                 Nf = 4
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 2
                 If (Abs(Ham_T2_max) > Zero) then
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 4 + No_Shift ! The bond (See above)
                 endif
              endif
              If (Abs(Ham_Tperp_max) > Zero) then
                 Nf = 5
                 this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                 this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = 3
              Endif
           enddo
        endif
        ! Test
        If (Test) then
           Write(6,*)  this(1)%N_FAM,  this(1)%L_FAM
           Write(6,*)  Ham_T1_max,Ham_T2_max, Ham_Tperp_max
           Do nf = 1,this(1)%N_FAM
              Do n = 1,this(1)%L_Fam(nf)
                 I =  this(1)%List_Fam(Nf,n,1)
                 nb = this(1)%List_Fam(Nf,n,2)
                 Write(6,"(I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,F6.3)")   Latt%list(I,1), Latt%list(I,2), this(1)%List(nb,1),this(1)%List(nb,2), &
                      &this(1)%List(nb,3), this(1)%List(nb,4), real(this(1)%T(nb))
              enddo
              Write(6,*)
           enddo
        endif



      end Subroutine Set_Default_hopping_parameters_Bilayer_square

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for the bilayer honeycomb lattice. \c Ham_T1 and \c Ham_T2 are the nearest-neighbour
!> hopping amplitudes on the first and second honeycomb layer, respectively, and
!> \c Ham_T_perp is the interlayer hopping amplitude.
!> Each honeycomb layer uses the same three-bond checkerboard decomposition as
!> Set_Default_hopping_parameters_honeycomb.
!>
!> @see Set_Default_hopping_parameters_honeycomb
!>
!> @param [out]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!> \endverbatim
!> @param [in]  Ham_T1_vec
!> \verbatim
!>   Real(:). In-plane nearest-neighbour hopping on layer 1 (orbitals 1-2) for each flavor.
!> \endverbatim
!> @param [in]  Ham_T2_vec
!> \verbatim
!>   Real(:). In-plane nearest-neighbour hopping on layer 2 (orbitals 3-4) for each flavor.
!> \endverbatim
!> @param [in]  Ham_Tperp_vec
!> \verbatim
!>   Real(:). Interlayer hopping amplitude for each flavor.
!> \endverbatim
!> @param [in]  Ham_Chem_vec
!> \verbatim
!>   Real(:). Chemical potential for each flavor.
!> \endverbatim
!> @param [in]  Phi_X_vec, Phi_Y_vec
!> \verbatim
!>   Real(:). Twist boundary conditions / flux along a1 and a2 for each flavor.
!> \endverbatim
!> @param [in]  Bulk
!> \verbatim
!>   Logical. Twist applied in bulk (.true.) or at boundary (.false.).
!> \endverbatim
!> @param [in]  N_Phi_vec
!> \verbatim
!>   Integer(:). Number of magnetic flux quanta for each flavor.
!> \endverbatim
!> @param [in]  N_FL
!> \verbatim
!>   Integer. Number of fermion flavors.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information (Norb=4 for two honeycomb layers).
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_Bilayer_honeycomb(this,Ham_T1_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
           &                                                      Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
           &                                                      List, Invlist, Latt, Latt_unit)

        Implicit none

        type (Hopping_Matrix_type), allocatable           :: this(:)
        Real (Kind=Kind(0.d0)), Intent(IN),dimension(:)  :: Ham_T1_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
        Integer, Intent(IN),dimension(:)                 :: N_Phi_vec
        Integer, Intent(IN)                   :: N_FL
        Logical, Intent(IN)                   :: Bulk
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit

        Real (Kind=Kind(0.d0))                :: Ham_T1_max, Ham_T2_max, Ham_Tperp_max

        ! Local
        Integer :: nf,N_Bonds, nc, I, I1, No_Shift, n, nb, no
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8
        Logical :: Test=.false.

        Ham_T1_max    = 0.d0
        Ham_T2_max    = 0.d0
        Ham_Tperp_max = 0.d0
        do nf = 1,N_FL
           if (abs(Ham_T1_vec   (nf)) > Ham_T1_max    ) Ham_T1_max    = abs(Ham_T1_vec(nf)   )
           if (abs(Ham_T2_vec   (nf)) > Ham_T2_max    ) Ham_T2_max    = abs(Ham_T2_vec(nf)   )
           if (abs(Ham_Tperp_vec(nf)) > Ham_Tperp_max ) Ham_Tperp_max = abs(Ham_Tperp_vec(nf))
        enddo

!!$        If (abs(Ham_T1_max) < Zero ) Then
!!$           Write(error_unit,*) 'At least Ham_T1 has to be bigger than zero'
!!$           error stop 1
!!$        endif


        Allocate( this(N_FL) )
        do nf = 1,N_FL
           N_bonds = 0
           N_bonds = N_bonds + 3
           if (abs(Ham_Tperp_max) > Zero )  N_bonds = N_bonds + 2
           if (abs(Ham_T2_max)    > Zero )  N_bonds = N_bonds + 3
           this(nf)%N_bonds =  N_Bonds
           Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                &    this(nf)%T(this(nf)%N_bonds) )
           nc = 0
           ! Bilayer honeycomb unit cell: 4 orbitals (orb1=A-layer1, orb2=B-layer1,
           !   orb3=A-layer2, orb4=B-layer2).  The 1+2 = 3 and 2+2 = 4 notation
           !   below simply means orb3 and orb4 (layer-2 counterparts of orb1 and orb2).
           !
           ! Bonds 1-3: layer-1 honeycomb bonds (same geometry as standalone honeycomb)
           ! Bond 1: A1(orb1) -> B1(orb2), ( 0,  0)  delta1 bond (same unit cell)
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  1
           this(nf)%List(nc,2) =  2
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  0

           ! Bond 2: B1(orb2) -> A1(orb1), ( 0, +1)  delta2 bond (+a2)
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  2
           this(nf)%List(nc,2) =  1
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  1

           ! Bond 3: A1(orb1) -> B1(orb2), (+1, -1)  delta3 bond (a1-a2)
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  1
           this(nf)%List(nc,2) =  2
           this(nf)%List(nc,3) =  1
           this(nf)%List(nc,4) = -1

           If (abs(Ham_Tperp_Max) > Zero )  then
              ! Bonds 4-5: interlayer (vertical) rung bonds connecting the two layers
              ! Bond 4: A1(orb1) -> A2(orb3), (0, 0)  A-sublattice interlayer rung
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  1
              this(nf)%List(nc,2) =  3
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  0

              ! Bond 5: B1(orb2) -> B2(orb4), (0, 0)  B-sublattice interlayer rung
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  2
              this(nf)%List(nc,2) =  4
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  0
           endif
           If (abs(Ham_T2_Max) > Zero )  then
              ! Bonds 3+No_Shift+1 to 3+No_Shift+3: layer-2 honeycomb bonds
              ! (orb3 = 1+2, orb4 = 2+2; No_Shift=2 when Tperp /= 0 to skip bonds 4-5)
              ! Bond 3+No_Shift+1: A2(orb3) -> B2(orb4), ( 0,  0)  delta1 bond
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  1 + 2
              this(nf)%List(nc,2) =  2 + 2
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  0

              ! Bond 3+No_Shift+2: B2(orb4) -> A2(orb3), ( 0, +1)  delta2 bond (+a2)
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  2 + 2
              this(nf)%List(nc,2) =  1 + 2
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  1

              ! Bond 3+No_Shift+3: A2(orb3) -> B2(orb4), (+1, -1)  delta3 bond (a1-a2)
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  1 + 2
              this(nf)%List(nc,2) =  2 + 2
              this(nf)%List(nc,3) =  1
              this(nf)%List(nc,4) = -1
           endif
           Allocate ( this(nf)%T_Loc(Latt_Unit%Norb) )
           do nc = 1,Latt_Unit%Norb
              this(nf)%T_Loc(nc)  = cmplx(-Ham_Chem_vec(nf),0.d0,kind(0.d0))
           enddo
           If (abs(Ham_Tperp_Max) < Zero .and. abs(Ham_T2_Max) < Zero ) then
              this(nf)%T_Loc(3) = cmplx(0.d0,0.d0,kind(0.d0))
              this(nf)%T_Loc(4) = cmplx(0.d0,0.d0,kind(0.d0))
           Endif
           this(nf)%N_Phi =  N_Phi_vec(nf)
           this(nf)%Phi_X =  Phi_X_vec(nf)
           this(nf)%Phi_Y =  Phi_Y_vec(nf)
           this(nf)%Bulk =   Bulk

        enddo

        ! Set Checkerboard (bilayer honeycomb)
        ! Each layer uses its own set of 3 honeycomb families (one per bond type).
        ! Since each bond type in one layer involves sites from only that layer,
        ! layer-1 and layer-2 bonds of the same type are placed in the same family.
        ! The interlayer rung bonds (if present) form a separate 4th family.
        !   Family 1: bond type 1 (delta1) from both layers -> all Latt%N [or 2*Latt%N] bonds
        !   Family 2: bond type 2 (delta2) from both layers
        !   Family 3: bond type 3 (delta3) from both layers
        !   Family 4 (if Tperp /= 0): both interlayer rung bond types (bonds 4 and 5)
        ! No_Shift=2 when Tperp /= 0: layer-2 bonds start at index 3+No_Shift=6.
        this(1)%N_FAM  = 3
        If ( abs(Ham_Tperp_Max) > Zero ) this(1)%N_FAM = 4
        Allocate (this(1)%L_Fam(this(1)%N_FAM),  this(1)%Prop_Fam(this(1)%N_FAM))
        this(1)%Prop_Fam= 1.d0

        No_Shift = 0
        If (abs(Ham_Tperp_Max) > Zero ) No_Shift=2

        If     ( abs(Ham_T2_Max)   <  Zero  .and. abs(Ham_Tperp_Max) < Zero)    then
           this(1)%L_FAM  = Latt%N
           Allocate (this(1)%List_Fam(this(1)%N_FAM,Latt%N,2))
        elseif ( abs(Ham_T2_Max)   <  Zero  .and. abs(Ham_Tperp_Max) > Zero)    then
           this(1)%L_FAM    =   Latt%N
           this(1)%L_Fam(4) = 2*Latt%N
           Allocate (this(1)%List_Fam(this(1)%N_FAM,2*Latt%N,2))
        elseif ( abs(Ham_T2_Max)   >  Zero  .and. abs(Ham_Tperp_Max) < Zero)    then
           this(1)%L_FAM    = 2*Latt%N
           Allocate (this(1)%List_Fam(this(1)%N_FAM,2*Latt%N,2))
        elseif ( abs(Ham_T2_Max)   >  Zero  .and. abs(Ham_Tperp_Max) > Zero)    then
           this(1)%L_FAM    = 2*Latt%N
           Allocate (this(1)%List_Fam(this(1)%N_FAM,2*Latt%N,2))
           No_Shift     = 2
        endif

        do I = 1,Latt%N
           Do  nf = 1,this(1)%N_FAM
              this(1)%List_Fam(nf,I,1) = I  ! Unit cell
              this(1)%List_Fam(nf,I,2) = nf ! The bond (See above)
           Enddo
        enddo
        if (abs(Ham_T2_Max)   >  Zero ) Then
           do I = 1,Latt%N
              Do  nf = 1,this(1)%N_FAM
                 this(1)%List_Fam(nf,I + Latt%N,1) = I                   ! Unit cell
                 this(1)%List_Fam(nf,I + Latt%N,2) = nf + 3 +  No_Shift  ! The bond (See above)
              Enddo
           enddo
        endif
        if (abs(Ham_Tperp_Max)   >  Zero ) Then
           do no = 0,1
              do I = 1,Latt%N
                 this(1)%List_Fam(4,I + no*Latt%N,1) = I       ! Unit cell
                 this(1)%List_Fam(4,I + no*Latt%N,2) = 4 + no  ! The bond (See above)
              Enddo
           enddo
        endif
        ! Test
        If (Test) then
           Write(6,*)  this(1)%N_FAM,  this(1)%L_FAM
           Write(6,*)  Ham_T1_Max,Ham_T2_Max, Ham_Tperp_Max
           Do nf = 1,this(1)%N_FAM
              Do n = 1,this(1)%L_Fam(nf)
                 I =  this(1)%List_Fam(Nf,n,1)
                 nb = this(1)%List_Fam(Nf,n,2)
                 Write(6,"(I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,F6.3)")   Latt%list(I,1), Latt%list(I,2), this(1)%List(nb,1),this(1)%List(nb,2), &
                      &this(1)%List(nb,3), this(1)%List(nb,4), Real(this(1)%T(nb))
              enddo
              Write(6,*)
           enddo
        endif

      end Subroutine Set_Default_hopping_parameters_Bilayer_honeycomb

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Given the checkerboard decomposition set by a \c Set_Default_hopping_parameters_* routine,
!> this routine generates the data for a symmetric Trotter decomposition.
!> The number of families is expanded from \c N_FAM to \c 2*N_FAM-1 by appending the original
!> sequence in reverse order.  The largest family is always placed in the centre slot so that
!> it receives the full time-step weight (\c Prop_Fam=1), while all satellite families receive
!> half the time-step weight (\c Prop_Fam=0.5).  This preserves the hermitian properties of the
!> hopping propagator when \c Symm=.true. is passed to Predefined_Hoppings_set_OPT.
!>
!> @param [inout] this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type (N_FL flavors).  On entry the standard
!>   N_FAM-family decomposition is set; on exit the symmetric 2*N_FAM-1 form is stored.
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Symmetrize_Families(this)
        implicit none


        type  (Hopping_Matrix_type), allocatable         :: this(:)
        ! In Families.  Out Symmetrized Families.

        !  Make a copy  of the unsymmetrized forms
        Integer                              ::  N_FAM_C
        Integer, allocatable                 ::  L_Fam_C(:),  List_Fam_C(:,:,:)
        Real (Kind=Kind(0.d0)), allocatable  ::  Prop_Fam_C(:)

        Integer :: n,n1,n2, n_f_max, n_l_max, nc
        Integer, allocatable ::  list_Fam_tmp(:)

        ! Copy
        N_FAM_C = this(1)%N_FAM
        Allocate(L_FAM_C(N_FAM_C))
        n2 = Size(this(1)%List_Fam,2)
        Allocate ( List_Fam_C(N_FAM_C,n2,2), Prop_Fam_C(N_FAM_C) )
        L_FAM_C    = this(1)%L_FAM
        List_Fam_C = this(1)%List_Fam
        Prop_Fam_C = this(1)%Prop_Fam

        ! Re-allocate
        this(1)%N_FAM  =  2*N_FAM_C - 1
        Deallocate (this(1)%L_Fam, this(1)%List_Fam, this(1)%Prop_Fam)
        Allocate   (this(1)%L_Fam(this(1)%N_FAM), this(1)%List_Fam(this(1)%N_FAM,n2,2), this(1)%Prop_Fam(this(1)%N_FAM) )

        ! Symmetrize
        ! Find the longest family.
        n_l_max = 0
        n_f_max = 0
        do n = 1, N_FAM_C
           if (L_FAM_C(n) > n_l_max ) then
              n_l_max = L_FAM_C(n)
              n_f_max = n
           endif
        enddo
        !Write(6,*) 'N_f_max' , n_f_max
        Allocate( list_Fam_tmp(this(1)%N_FAM) )
        nc = 0
        Do n = 1, N_FAM_C
           nc = nc + 1
           list_Fam_tmp(nc) = n
        Enddo
        Do n = N_FAM_C-1,1,-1
           nc = nc + 1
           list_Fam_tmp(nc) = n
        Enddo

        ! Place the largest familly in the middle and set the time step.
        this(1)%Prop_Fam         = 0.5D0
        this(1)%Prop_Fam(N_FAM_C) = 1.D0
        If (N_F_Max .ne. N_FAM_C )  then
           list_Fam_tmp(N_FAM_C)        = n_f_max
           list_Fam_tmp(1)              = N_FAM_C
           list_Fam_tmp(this(1)%N_FAM ) = N_Fam_C
        endif

        do n = 1,this(1)%N_FAM
           n1 = list_Fam_tmp(n)
           this(1)%L_Fam(n)        = L_FAM_C(n1)
           this(1)%List_Fam(n,:,:) = List_Fam_C(n1,:,:)
        enddo

        ! Clean
        Deallocate( L_FAM_C, List_Fam_C, Prop_Fam_C, List_Fam_tmp )

        !Write(6,*)  this(1)%N_FAM
        !Write(6,*)  this(1)%L_FAM
        !Write(6,*)  this(1)%Prop_Fam

      end Subroutine Symmetrize_Families

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Given a bond (I, J) or a site I=J, search whether it belongs to the set of
!> pinned vertices.  Returns the 1-based index into \c pinned_vertices if found,
!> or 0 if not found.  The search is symmetric: the pair (I,J) and (J,I) both match.
!>
!> @param [in]  I, J
!> \verbatim
!>   Integer. Global site indices of the bond endpoint (I=J for on-site terms).
!> \endverbatim
!> @param [in]  N_pinned_vertices
!> \verbatim
!>   Integer. Number of pinned vertices (first dimension of pinned_vertices).
!> \endverbatim
!> @param [in]  pinned_vertices
!> \verbatim
!>   Integer(N_pinned_vertices, 2). Each row stores the two site indices of a pinned bond/site.
!> \endverbatim
!> @return  1-based index into pinned_vertices if the bond is pinned, 0 otherwise.
!--------------------------------------------------------------------
      integer pure function get_i_pinned_vertex(I, J, N_pinned_vertices, pinned_vertices)
         integer, intent(in) :: I, J, N_pinned_vertices, pinned_vertices(N_pinned_vertices, 2)
         integer :: n

         get_i_pinned_vertex = 0
         do n = 1, N_pinned_vertices
            if ((I == pinned_vertices(n,1) .and. J == pinned_vertices(n,2)) .or. &
                (J == pinned_vertices(n,1) .and. I == pinned_vertices(n,2))) then
               get_i_pinned_vertex = n
               exit
            endif
         enddo
      end function get_i_pinned_vertex


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Given a bond (I, J) or a site I=J, return the multiplicative pinning factor for that bond
!> in flavor \c nf.  If the bond is not in the pinning list the routine returns 1.0 (no modification).
!> Sets the module-level flag \c pinning_notice_issued if the returned factor differs from 1.
!>
!> @param [in]  I, J
!> \verbatim
!>   Integer. Global site indices of the bond endpoints.
!> \endverbatim
!> @param [in]  N_pinned_vertices
!> \verbatim
!>   Integer. Number of pinned bonds/sites.
!> \endverbatim
!> @param [in]  pinned_vertices
!> \verbatim
!>   Integer(N_pinned_vertices, 2). Site-index pairs of pinned bonds/sites.
!> \endverbatim
!> @param [in]  pinning_factor
!> \verbatim
!>   Complex(:,:). Shape (N_pinned_vertices, N_FL). Multiplicative factor for each pinned bond and flavor.
!> \endverbatim
!> @param [in]  nf
!> \verbatim
!>   Integer. Flavor index (1..N_FL).
!> \endverbatim
!> @return  The pinning factor for bond (I,J) in flavor nf, or cmplx(1,0) if not pinned.
!--------------------------------------------------------------------
      complex(Kind=Kind(0.d0)) function get_pinning_factor(I, J, N_pinned_vertices, pinned_vertices, pinning_factor, nf)
         integer, intent(in) :: I, J, N_pinned_vertices, pinned_vertices(N_pinned_vertices, 2), nf
         complex(Kind=Kind(0.d0)), Intent(IN) :: pinning_factor(:,:)

         integer :: i_pinned_vertex

         i_pinned_vertex = get_i_pinned_vertex(I, J, N_pinned_vertices, pinned_vertices)
         if(i_pinned_vertex .ne. 0) then
            get_pinning_factor = pinning_factor(i_pinned_vertex, nf)
         else
            get_pinning_factor = cmplx(1.d0,0.d0, kind(0.d0))
         endif
         if( abs(get_pinning_factor  - cmplx(1.d0,0.d0, kind(0.d0))) >= 1.d-10 ) pinning_notice_issued = .true.
      end function get_pinning_factor

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Given a bond (I, J) or a site I=J, return the additive pinning offset for that bond
!> in flavor \c nf.  If the bond is not in the pinning list the routine returns 0.0 (no modification).
!> Sets the module-level flag \c pinning_notice_issued if the returned offset is non-zero.
!>
!> @param [in]  I, J
!> \verbatim
!>   Integer. Global site indices of the bond endpoints.
!> \endverbatim
!> @param [in]  N_pinned_vertices
!> \verbatim
!>   Integer. Number of pinned bonds/sites.
!> \endverbatim
!> @param [in]  pinned_vertices
!> \verbatim
!>   Integer(N_pinned_vertices, 2). Site-index pairs of pinned bonds/sites.
!> \endverbatim
!> @param [in]  pinning_offset
!> \verbatim
!>   Complex(:,:). Shape (N_pinned_vertices, N_FL). Additive offset for each pinned bond and flavor.
!> \endverbatim
!> @param [in]  nf
!> \verbatim
!>   Integer. Flavor index (1..N_FL).
!> \endverbatim
!> @return  The pinning offset for bond (I,J) in flavor nf, or cmplx(0,0) if not pinned.
!--------------------------------------------------------------------
      complex(Kind=Kind(0.d0)) function get_pinning_offset(I, J, N_pinned_vertices, pinned_vertices, pinning_offset, nf)
         integer, intent(in) :: I, J, N_pinned_vertices, pinned_vertices(N_pinned_vertices, 2), nf
         complex(Kind=Kind(0.d0)), Intent(IN) :: pinning_offset(:,:)

         integer :: i_pinned_vertex

         i_pinned_vertex = get_i_pinned_vertex(I, J, N_pinned_vertices, pinned_vertices)
         if(i_pinned_vertex .ne. 0) then
            get_pinning_offset = pinning_offset(i_pinned_vertex, nf)
         else
            get_pinning_offset = cmplx(0.d0,0.d0, kind(0.d0))
         endif
         if( abs(get_pinning_offset) >= 1.d-10 ) pinning_notice_issued = .true.
      end function get_pinning_offset

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This routine issues a warning if  pinning is used in the simulation. In particular it notifies the user that translation
!> symmetry is broken. 
!--------------------------------------------------------------------
      subroutine issue_pinning_notice()
         Character (len=64) :: file_info
         integer :: unit_info
#ifdef MPI
         integer :: Ierr, Isize, Irank, irank_g, isize_g, igroup
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
#endif

         If (first_pinning_notice_issued) then 
#if defined(TEMPERING)
            write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
            If (Irank_g == 0 ) write(error_unit, '(A,I0,A)') &
               'Warning: you are using pinning on parameter set ', igroup, &
               '. Results will not have translation symmetry.'
#else
            file_info = "info"
#endif

#ifdef MPI
            If (Irank_g == 0) &
#endif
               write(error_unit, *) 'Warning: you are using pinning, results will not have translation symmetry.'

#if defined(MPI)
            If (Irank_g == 0 ) then
#endif
               Open (newunit=unit_info, file=file_info, status="unknown", position="append")
               Write(unit_info,*) ' Pinning is used. Results will not have translation symmetry.'
               close(unit_info)
#if defined(MPI)
            endif
#endif
            first_pinning_notice_issued = .false.
         endif
      end subroutine issue_pinning_notice

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Allocate and set the hopping propagator array \c Op_T from the pre-built hopping matrix \c this.
!> Three propagation modes are supported depending on \c Checkerboard and \c Symm:
!>   - **Zero / Diagonal** hopping: one operator per site or a single diagonal block.
!>   - **Full matrix without checkerboard**: one dense operator block per flavor.
!>   - **Checkerboard decomposition**: one operator per checkerboard family per flavor;
!>     if \c Symm=.true., Symmetrize_Families is called first to produce a symmetric
!>     Trotter decomposition \f$ e^{-\Delta\tau H_t} \approx \prod_{n=1}^{2N-1} e^{-\Delta\tau_n H_t(n)} \f$.
!>
!> @param [inout]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!>   If Symm=.true. and Checkerboard=.true. its N_FAM is expanded by Symmetrize_Families.
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information.
!> \endverbatim
!> @param [in]  Dtau
!> \verbatim
!>   Real. Imaginary-time step.
!> \endverbatim
!> @param [in]  Checkerboard
!> \verbatim
!>   Logical. If .true., use the checkerboard decomposition stored in this.
!> \endverbatim
!> @param [in]  Symm
!> \verbatim
!>   Logical. If .true. and Checkerboard=.true., apply a symmetric Trotter decomposition.
!> \endverbatim
!> @param [out]  Op_T
!> \verbatim
!>   Type(Operator)(:,:), allocatable. Shape (N_ops, N_FL).
!>   N_ops = 1 for full-matrix mode, Ndim for diagonal, or sum of L_Fam for checkerboard.
!> \endverbatim
!> @param [in]  pinned_vertices  (optional)
!> \verbatim
!>   Integer(:,:). Shape (N_pinned_vertices, 2). Site-index pairs of bonds/sites whose
!>   hopping matrix elements are to be selectively modified.
!>   Must be supplied together with pinning_factor and pinning_offset.
!> \endverbatim
!> @param [in]  pinning_factor  (optional)
!> \verbatim
!>   Complex(:,:). Shape (N_pinned_vertices, N_FL). Multiplicative factor applied to
!>   each pinned bond matrix element per flavor.
!> \endverbatim
!> @param [in]  pinning_offset  (optional)
!> \verbatim
!>   Complex(:,:). Shape (N_pinned_vertices, N_FL). Additive offset applied to
!>   each pinned bond matrix element per flavor.
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine Predefined_Hoppings_set_OPT(this,List,Invlist,Latt,  Latt_unit,  Dtau,Checkerboard, Symm, OP_T,  & 
                                             & pinned_vertices, pinning_factor, pinning_offset)

        Implicit none

        type (Hopping_Matrix_type), allocatable             :: this(:)
        Integer, Intent(IN), Dimension(:,:)                 :: List, Invlist
        Type(Lattice),  Intent(in)                          :: Latt
        Type(Unit_cell),Intent(in)                          :: Latt_unit
        Real (Kind=Kind(0.d0)), Intent(In)                  :: Dtau
        Logical, Intent(IN)                                 :: Checkerboard, Symm
        
        Type(Operator), Intent(Out),  dimension(:,:), allocatable  :: Op_T
        
        ! Indices of pinned vertices. Shape [N_pinned_vertices, 2]
        Integer, Intent(IN), optional                       :: pinned_vertices(:,:)
        ! Factor, by which the vertex matrix elements will get multiplied. Shape [N_pinned_vertices, N_FL]
        complex(Kind=Kind(0.d0)), Intent(IN), optional      :: pinning_factor(:,:)
        ! Offset of which the vertex matrix element Shape [N_pinned_vertices, N_FL]
        complex(Kind=Kind(0.d0)), Intent(IN), optional      :: pinning_offset(:,:)
        
        
        ! Local
        Integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
        Integer                           :: n_1, n_2, Nb, n_f,l_f, n_l, N, nc, orb
        Real   (Kind=Kind(0.d0))          :: Ham_T, Ham_Chem,  Phi_X, Phi_Y
        Logical                           :: Bulk
        Complex(Kind=Kind(0.d0))          :: Z,  Z1, Z2

        Integer                           :: N_pinned_vertices, i_pinned_vertex

        
        N_FL =  size(this,1)
        !Write(6,*)  'N_FL ', N_FL
        Ndim =  Latt%N * Latt_Unit%Norb

        if( (present(pinned_vertices) .and. .not.(present(pinning_factor ).and. present(pinning_offset))) .or. &
          & (present(pinning_factor)  .and. .not.(present(pinned_vertices).and. present(pinning_offset))) .or. &
          & (present(pinning_offset)  .and. .not.(present(pinned_vertices).and. present(pinning_factor))) ) then
           write(error_unit, *) 'All pinned_vertices, pinning_factor and pinning_offset need to be supplied for pinning.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        if(present(pinned_vertices)) N_pinned_vertices = size(pinned_vertices, 1)
        if(present(pinned_vertices)) then
           if(size(pinned_vertices, 2) .ne. 2) then
             write(error_unit, *) 'Second dimension of pinned_vertices has to be 2.'
             CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif

           if(size(pinning_factor, 1) .ne. N_pinned_vertices) then
            write(error_unit, *) 'First dimension of pinning_factor has to be equal to the number of pinned vertices.'
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          endif

          if(size(pinning_factor, 2) .ne. N_FL) then
            write(error_unit, *) 'Second dimension of pinning_factor has to be equal to the number of flavors N_FL.'
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          endif
        endif

        ! Test of correctness of checkerboard decomposition
        If (checkerboard) then
           if (.not.(test_checkerboard_decomposition(this(1), Latt, invlist)))  CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        end If

        select case (inquire_hop(this))
        case(0)  !  Zero
           allocate(Op_T(1,N_FL))
           do nf = 1,N_FL
              Call Op_make(Op_T(1,nf),1)
              Op_T(1,nf)%P(1)   = 1
              Op_T(1,nf)%O(1,1) = cmplx(0.d0,0.d0, kind(0.d0))
              Op_T(1,nf)%g      = 0.d0
              Op_T(1,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
              Call Op_set(Op_T(1,nf))
           enddo
        case(1)  ! Diagonal
           allocate(Op_T(Ndim,N_FL))
           do nf = 1,N_FL
              do n = 1,ndim
                  Call Op_make(Op_T(n,nf),1)
                  Op_T(n,nf)%P(1)   = n
                  Z = this(nf)%T_Loc(list(n,2))
                  if(present(pinned_vertices)) then
                     Z = Z*get_pinning_factor(n, n, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + &
                        &  get_pinning_offset(n, n, N_pinned_vertices, pinned_vertices, pinning_offset, nf)
                  endif
                  Op_T(n,nf)%O(1,1) =  Z
                  Op_T(n,nf)%g      = -Dtau
                  Op_T(n,nf)%alpha  =  cmplx(0.d0,0.d0, kind(0.D0))
                  Call Op_set(Op_T(n,nf))
              enddo
           enddo
        case default
           If ( .not. Checkerboard) then
              allocate(Op_T(1,N_FL))
              do nf = 1,N_FL
                 !Write(6,*)
                 Call Op_make(Op_T(1,nf),Ndim)   ! This is too restrictive for the Kondo type models. The hopping only occurs on one subsystem.
                 N_Phi     = this(nf)%N_Phi
                 Phi_X     = this(nf)%Phi_X
                 Phi_Y     = this(nf)%Phi_Y
                 Bulk      = this(nf)%Bulk
                 !Write(6,*) N_Phi, Phi_X,Phi_Y, Bulk
                 !Write(6,*) This(nf)%list
                 DO I = 1, Latt%N
                    do Nb = 1, this(nf)%N_bonds
                       no_I = this(nf)%list(Nb,1)
                       no_J = this(nf)%list(Nb,2)
                       n_1  = this(nf)%list(Nb,3)
                       n_2  = this(nf)%list(Nb,4)
                       J    = Latt%nnlist(I,n_1,n_2)
                       Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                       I1   = Invlist(I,no_I)
                       J1   = Invlist(J,no_J)
                       if(present(pinned_vertices)) then
                          Z = Z*get_pinning_factor(I1, J1, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + & 
                              & get_pinning_offset(I1, J1, N_pinned_vertices, pinned_vertices, pinning_offset, nf)
                       endif
                       Op_T(1,nf)%O(I1,J1) = this(nf)%T(Nb)*Z
                       Op_T(1,nf)%O(J1,I1) = Conjg(this(nf)%T(Nb)*Z)
                    enddo
                    ! T(N_b=1..N_bonds)
                    ! List(N_b,1) = no_1
                    ! List(N_b,2) = no_2
                    ! List(N_b,3) = n_1
                    ! List(N_b,4) = n_2
                    ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)
                    Do no_I = 1, Latt_Unit%Norb
                        I1   = Invlist(I,no_I)
                        Z = this(nf)%T_Loc(no_I) 
                        if(present(pinned_vertices)) then
                           Z = Z*get_pinning_factor(I1, I1, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + &
                              &  get_pinning_offset(I1, I1, N_pinned_vertices, pinned_vertices, pinning_offset, nf)
                        endif
                        Op_T(1,nf)%O(I1,I1) = Z
                    Enddo
                 enddo
                 Do I = 1,Ndim
                    Op_T(1,nf)%P(i) = i
                 Enddo
                 Op_T(1,nf)%g = -Dtau
                 Op_T(1,nf)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                 Call Op_set(Op_T(1,nf))
                 !Do I = 1,Size(Op_T(1,nf)%E,1)
                 !   Write(6,*) Op_T(1,nf)%E(I)
                 !Enddo
              Enddo
           Elseif (Checkerboard) then
              If (Symm) Call Symmetrize_families(this)
              N = 0
              do n_f = 1,this(1)%N_FAM
                 N = N +  this(1)%L_Fam(n_f)
              enddo
              allocate(Op_T(N,N_FL))
              do nf = 1,N_FL
                 ! Compute Multiplicity: count how many checkerboard bond types contain each orbital.
                 ! An orbital can appear as source (List(i,1)) or target (List(i,2)) of multiple bond
                 ! types, so Multiplicity(orb) >= 1.  The on-site term T_Loc(orb) is later divided by
                 ! Multiplicity(orb) so that, when all checkerboard operators that touch orb are summed,
                 ! the full T_Loc(orb) is recovered without double-counting.
                 allocate(this(nf)%Multiplicity(Latt_Unit%Norb))
                 this(nf)%Multiplicity = 0
                 do i = 1, size(this(nf)%List, 1)
                    orb = this(nf)%List(i, 1)
                    this(nf)%Multiplicity(orb) = this(nf)%Multiplicity(orb) + 1
                    orb = this(nf)%List(i, 2)
                    this(nf)%Multiplicity(orb) = this(nf)%Multiplicity(orb) + 1
                 end do

                 N_Phi     = this(nf)%N_Phi
                 Phi_X     = this(nf)%Phi_X
                 Phi_Y     = this(nf)%Phi_Y
                 Bulk      = this(nf)%Bulk
                 do nc = 1, Size(Op_T,1)
                    Call Op_make(Op_T(nc,nf),2)
                 enddo
                 nc = 0
                 Do n_f = 1, this(1)%N_FAM
                    Do l_f = 1, this(1)%L_Fam(n_f)
                       I  = this(1)%List_Fam(n_f,l_f,1)
                       nb = this(1)%List_Fam(n_f,l_f,2)
                       no_I = this(nf)%list(Nb,1)
                       no_J = this(nf)%list(Nb,2)
                       n_1  = this(nf)%list(Nb,3)
                       n_2  = this(nf)%list(Nb,4)
                       J    = Latt%nnlist(I,n_1,n_2)
                       I1   = Invlist(I,no_I)
                       J1   = Invlist(J,no_J)
                       ! Z: Peierls gauge factor for the hop (i,no_I) -> (j,no_J); see Generic_hopping.
                       Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                       ! Z1, Z2: the fraction of the on-site (chemical-potential) energy assigned to
                       ! this 2x2 block for sites I1 and J1 respectively.  Dividing by Multiplicity
                       ! ensures T_Loc is not over-counted when multiple checkerboard families share
                       ! the same orbital.
                       Z1   = this(nf)%T_loc(no_I)/this(1)%Multiplicity(no_I)
                       Z2   = this(nf)%T_loc(no_J)/this(1)%Multiplicity(no_J)
                       if(present(pinned_vertices)) then
                          Z = Z*get_pinning_factor(I1, J1, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + & 
                              & get_pinning_offset(I1, J1, N_pinned_vertices, pinned_vertices, pinning_offset, nf)
                          Z1 = Z1*get_pinning_factor(I1, I1, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + &
                              &   get_pinning_offset(I1, I1, N_pinned_vertices, pinned_vertices, pinning_offset, nf)/this(1)%Multiplicity(no_I)
                          Z2 = Z2*get_pinning_factor(J1, J1, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + &
                              &   get_pinning_offset(J1, J1, N_pinned_vertices, pinned_vertices, pinning_offset, nf)/this(1)%Multiplicity(no_J)
                       endif
                       nc = nc + 1
                       Op_T(nc,nf)%P(1) = I1           ! global index of source site
                       Op_T(nc,nf)%P(2) = J1           ! global index of target site
                       ! The 2x2 operator block in the basis (I1, J1):
                       Op_T(nc,nf)%O(1,2) = this(nf)%T(Nb)*Z      ! off-diagonal: hopping I1->J1 (with Peierls phase)
                       Op_T(nc,nf)%O(2,1) = Conjg(this(nf)%T(Nb)*Z) ! conjugate: hopping J1->I1 (hermiticity)
                       Op_T(nc,nf)%O(1,1) = Z1  ! on-site energy at I1, weighted by 1/Multiplicity(no_I)
                       Op_T(nc,nf)%O(2,2) = Z2  ! on-site energy at J1, weighted by 1/Multiplicity(no_J)
                       Op_T(nc,nf)%g = -Dtau*this(1)%Prop_Fam(n_f)
                       Op_T(nc,nf)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                       Call Op_set(Op_T(nc,nf))
                    Enddo
                 enddo
              enddo
           endif
        end select

        If (pinning_notice_issued) call issue_pinning_notice()

      end Subroutine Predefined_Hoppings_set_OPT

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Compute the kinetic energy \f$ \langle H_t \rangle \f$ from the equal-time
!> Green's function \c GRC using the generic hopping stored in \c this.
!> Handles zero, diagonal, and full (non-checkerboard) hopping; the checkerboard
!> storage is not used here since kinetic energy only requires the bond list.
!>
!> @param [in]  this
!> \verbatim
!>   Allocatable array of Hopping_Matrix_type, dimension(N_FL).
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>   Integer(:,:). Site-to-(unit-cell, orbital) map and its inverse.
!> \endverbatim
!> @param [in]  Latt
!> \verbatim
!>   Type(Lattice). Lattice geometry.
!> \endverbatim
!> @param [in]  Latt_unit
!> \verbatim
!>   Type(Unit_cell). Unit-cell information.
!> \endverbatim
!> @param [in]  GRC
!> \verbatim
!>   Complex(:,:,:). Equal-time Green's function GRC(I,J,nf) = <c_I c_J^dag>_nf.
!> \endverbatim
!> @param [out]  Z_Kin
!> \verbatim
!>   Complex. Kinetic energy summed over all sites and flavors.
!> \endverbatim
!> @param [in]  pinned_vertices  (optional)
!> \verbatim
!>   Integer(:,:). Shape (N_pinned_vertices, 2). Pinned bond site-index pairs.
!> \endverbatim
!> @param [in]  pinning_factor  (optional)
!> \verbatim
!>   Complex(:,:). Multiplicative pinning factors, shape (N_pinned_vertices, N_FL).
!> \endverbatim
!> @param [in]  pinning_offset  (optional)
!> \verbatim
!>   Complex(:,:). Additive pinning offsets, shape (N_pinned_vertices, N_FL).
!> \endverbatim
!
!--------------------------------------------------------------------
      Subroutine  Predefined_Hoppings_Compute_Kin(this,List,Invlist, Latt, Latt_unit, GRC, Z_Kin, & 
                                             &    pinned_vertices, pinning_factor, pinning_offset)

        Implicit none

        type (Hopping_Matrix_type), allocatable  :: this(:)
        Integer, Intent(IN), Dimension(:,:)                 :: List, Invlist
        Type(Lattice),  Intent(in)                          :: Latt
        Type(Unit_cell),Intent(in)                          :: Latt_unit   
        Complex (Kind=Kind(0.d0)), intent(in), Dimension(:,:,:) :: GRC(:,:,:)
        Complex (Kind=Kind(0.d0)),  intent(out) :: Z_kin
        ! Indices of pinned vertices. Shape [N_pinned_vertices, 2]
        Integer, Intent(IN), optional                       :: pinned_vertices(:,:)
        ! Factor, by which the vertex matrix elements will get multiplied. Shape [N_pinned_vertices, N_FL]
        complex(Kind=Kind(0.d0)), Intent(IN), optional      :: pinning_factor(:,:)
        ! Offset of which the vertex matrix element Shape [N_pinned_vertices, N_FL]
        complex(Kind=Kind(0.d0)), Intent(IN), optional      :: pinning_offset(:,:)
        

        !Local
        Integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
        Integer                           :: n_1, n_2, Nb, n_f,l_f, n_l, N, nc
        Real   (Kind=Kind(0.d0))          :: Ham_T, Ham_Chem,  Phi_X, Phi_Y
        Logical                           :: Bulk
        Complex(Kind=Kind(0.d0))          :: Z
        Integer                           :: N_pinned_vertices, i_pinned_vertex

        if( (present(pinned_vertices) .and. .not.(present(pinning_factor ).and. present(pinning_offset))) .or. &
          & (present(pinning_factor)  .and. .not.(present(pinned_vertices).and. present(pinning_offset))) .or. &
          & (present(pinning_offset)  .and. .not.(present(pinned_vertices).and. present(pinning_factor))) ) then
           write(error_unit, *) 'All pinned_vertices, pinning_factor and pinning_offset need to be supplied for pinning.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif


        if(present(pinned_vertices)) N_pinned_vertices = size(pinned_vertices, 1)
        if(present(pinned_vertices)) then
            if(size(pinned_vertices, 2) .ne. 2) then
               write(error_unit, *) 'Second dimension of pinned_vertices has to be 2.'
               call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            endif

           if(size(pinning_factor, 1) .ne. N_pinned_vertices) then
               write(error_unit, *) 'First dimension of pinning_factor has to be equal to the number of pinned vertices.'
               call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif

           if(size(pinning_offset, 1) .ne. N_pinned_vertices) then
               write(error_unit, *) 'First dimension of pinning_offset has to be equal to the number of pinned vertices.'
               call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif

           if(size(pinning_factor, 2) .ne. N_FL) then
               write(error_unit, *) 'Second dimension of pinning_factor has to be equal to the number of flavors N_FL.'
               call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
        endif

        select case (inquire_hop(this))
        case(0)  !  Zero
           Z_Kin = cmplx(0.d0,0.d0,Kind(0.d0))
        case(1)
           Z_Kin = cmplx(0.d0,0.d0,Kind(0.d0))
           N_FL  =  Size(GRC,3)
           do nf = 1,N_FL
              do I = 1, Latt%N
                 Do no_I = 1, Latt_Unit%Norb
                    I1   = Invlist(I,no_I)
                    Z_Kin = Z_Kin   +  this(nf)%T_Loc(no_I)*GRC(I1,I1,nf)
                 Enddo
              enddo
           enddo
        case default
           N_FL  =  Size(GRC,3)
           Z_Kin = cmplx(0.d0,0.d0,Kind(0.d0))
           do nf = 1,N_FL
              N_Phi     = this(nf)%N_Phi
              Phi_X     = this(nf)%Phi_X
              Phi_Y     = this(nf)%Phi_Y
              Bulk      = this(nf)%Bulk
              DO I = 1, Latt%N
                 do Nb = 1, this(nf)%N_bonds
                    no_I = this(nf)%list(Nb,1)
                    no_J = this(nf)%list(Nb,2)
                    n_1  = this(nf)%list(Nb,3)
                    n_2  = this(nf)%list(Nb,4)
                    J    = Latt%nnlist(I,n_1,n_2)
                    Z    = this(nf)%T(Nb)*Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                    I1   = Invlist(I,no_I)
                    J1   = Invlist(J,no_J)
                    if(present(pinned_vertices)) then
                          Z = Z*get_pinning_factor(I1, J1, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + & 
                              & get_pinning_offset(I1, J1, N_pinned_vertices, pinned_vertices, pinning_offset, nf)
                    endif
                    Z_Kin = Z_Kin + Z * GRC(I1,J1,nf) + conjg(Z)*GRC(J1,I1,nf)
                 enddo
                 Do no_I = 1, Latt_Unit%Norb
                    I1   = Invlist(I,no_I)
                    Z = this(nf)%T_Loc(no_I)
                    if(present(pinned_vertices)) then
                          Z = Z*get_pinning_factor(I1, I1, N_pinned_vertices, pinned_vertices, pinning_factor, nf) + & 
                              & get_pinning_offset(I1, I1, N_pinned_vertices, pinned_vertices, pinning_offset, nf)
                    endif
                    Z_Kin = Z_Kin   +  Z*GRC(I1,I1,nf)
                 Enddo
              enddo
           enddo
        end select

      end Subroutine Predefined_Hoppings_Compute_Kin
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!>    Hopping, with or without checkerboard
!>    Per flavor, the  hopping is given by
!>    \f[  e^{ - \Delta \tau  H_t  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n) }   \f]
!>    If  _Symm_ is set to true and if  _Checkeboard_ is on, then  one will carry out a
!>    symmetric decomposition so as to preserve  the hermitian properties of the hopping.
!>    Thereby   OP_T has dimension OP_T(N,N_FL)
!> @param [in]  Latttice_type
!>    Character(64)
!>\verbatim
!>     Square,  Honeycomb, Pi_Flux
!>\endverbatim
!> @param [in]  Latt_unit
!>    Type(Unit_cell)
!> \verbatim
!>     Contains number of orbitals per unit cell and positions, as well as coordination number
!> \endverbatim
!> @param [in]  Ndim
!>    Integer
!> \verbatim
!>     Number of orbitals
!> \endverbatim
!> @param [in]  List, Invlist
!>    Integer(:,:)
!> \verbatim
!>      List(I=1.. Ndim,1)    =   Unit cell of site I
!>      List(I=1.. Ndim,2)    =   Orbital index  of site I
!>      Invlist(Unit_cell,Orbital) = site I
!> \endverbatim
!> @param [in]    Latt
!>    Type(Lattice)
!> \verbatim
!>      The Lattice
!> \endverbatim
!> @param [in]  Dtau
!>    Real
!> \verbatim
!>      Imaginary time step
!> \endverbatim
!> @param [in]  Ham_T
!>    Real
!> \verbatim
!>      Hopping matrix element
!> \endverbatim
!> @param [in]  Ham_Chem
!>    Real
!> \verbatim
!>      Chemical potential
!> \endverbatim
!> @param [in]  XB_X, YB_Y
!>    Real
!> \verbatim
!>      X, Y  Boundary conditions
!> \endverbatim
!> @param [in]  Phi_X, Phi_Y
!>    Real
!> \verbatim
!>      X, Y  Fluxes
!> \endverbatim
!> @param [in]  N_FL
!>    Integer
!> \verbatim
!>      Flavors
!> \endverbatim
!> @param [in]  Checkerboard
!>    Logical
!> \verbatim
!>      Allows for checkerboard decomposition
!> \endverbatim
!> @param [in]  Symm
!>    Logical
!> \verbatim
!>      Allows for symmetric checkerboard decomposition
!> \endverbatim
!> @param [out]  OP_T
!>    Type(operator)(N,N_FL)
!> \verbatim
!>      Hopping
!> \endverbatim
!> @param [in]  Dimer
!>    Real, Optional.  Modulation of hopping that breaks lattice symmetries so as to generate a unique
!>    ground state for the half-filled case.  This option is  effective only  if the checkerboard
!>    decomposition is not used. It is presently implemented for the square and one-dimensional lattices.
!> \verbatim
!>      Hopping
!> \endverbatim
!>

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This  function provides generic hopping.
!>
!--------------------------------------------------------------------

      complex  (Kind=kind(0.d0)) function Generic_hopping(i,no_i, del_1, del_2, no_j, N_Phi, Flux_1,Flux_2, Bulk, Latt, Latt_Unit)

        Use Lattices_v3
        Implicit none


        Integer        ,  Intent(In) :: N_Phi, i, no_i, del_1, del_2, no_j
        Type(Unit_cell),  Intent(In) :: Latt_Unit
        Type(Lattice)  ,  Intent(In) :: Latt
        Real (Kind = Kind(0.d0)), intent(In) :: Flux_1,Flux_2
        Logical        ,  Intent(In) :: Bulk


        !Local
        Integer                   :: j, N1, N2,n
        real (Kind=Kind(0.d0))    :: xj_p(2), xi_p(2), xjp_p(2), del_p(2), A_p(2), pi, XB_p(2), V, B, Zero, x_p(2), x1_p(2)

        Complex (Kind=Kind(0.d0)) :: Z_hop


        ! Initialize the gauge factor to 1; three independent contributions are multiplied below.
        ! The final result is Z_hop = exp(i * total_phase), where the phase accumulates from:
        !   1. Aharonov-Bohm / twist boundary conditions  (Flux_1, Flux_2)
        !   2. Orbital magnetic field in Landau gauge      (N_Phi flux quanta through the system)
        !   3. Boundary-crossing correction for the Landau gauge when the hop wraps around the
        !      periodic boundary (N1 /= 0 and/or N2 /= 0 winding numbers)
        Z_hop = cmplx(1.d0,0.d0,kind(0.d0))

        ! Compute the real-space position of the target site j = i + del_1*a1 + del_2*a2 (before PBC).
        xj_p =  real(latt%list(i,1) + del_1 ,kind(0.d0)) * latt%a1_p  +  real(latt%list(i,2) + del_2 ,kind(0.d0)) * latt%a2_p
        ! Fold xj_p back into the simulation cell: xj_p = xjp_p + N1*L1_p + N2*L2_p,
        ! where xjp_p is the image of j inside the simulation cell and (N1,N2) are the winding numbers.
        N1 = 0; N2 = 0
        Call npbc(xjp_p, xj_p, Latt%L1_p, Latt%L2_p,  N1, N2)
        XB_p = real(N1,kind(0.d0))*Latt%L1_p  +  real(N2,kind(0.d0))*Latt%L2_p
        Do n = 1,2
           xj_p (n) = xj_p (n) + Latt_unit%Orb_pos_p(no_j,n)    ! add orbital offset at target
           xjp_p(n) = xjp_p(n) + Latt_unit%Orb_pos_p(no_j,n)
        enddo
        xi_p    = real(latt%list(i,1), kind(0.d0)) * latt%a1_p  +  real(latt%list(i,2),kind(0.d0)) * latt%a2_p
        Do n = 1,2
           xi_p(n) = xi_p(n) +  Latt_unit%Orb_pos_p(no_i,n)     ! add orbital offset at source
        Enddo
        !!Check that  xjp_p(:) + XB_p(:) =  xj_p(:)
        !!x1_p(:) = xjp_p(:) + XB_p
        !!Write(6,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,I2,3x,I2)")  x1_p(1),x1_p(2), xj_p(1), xj_p(2), N1,N2
        !! -->  i + del_1*a_1 + del_2* a_2  =  i' + N1*L1_p + N2*L2_p  with i' in the set of lattice points.


        ! del_p: the real-space displacement vector of the hop (including orbital offsets).
        del_p  =  xj_p - xi_p

        ! ---------- Contribution 1: Aharonov-Bohm / twist boundary condition ----------
        ! A_p is the vector potential that implements twist phases Flux_1 and Flux_2 along
        ! the two lattice directions.  bZ1_p, bZ2_p are the reciprocal-lattice vectors dual
        ! to L1_p, L2_p, so A_p · L_{1,2} = 2*pi*Flux_{1,2}.
        pi = acos(-1.d0)
        A_p(:)  =   Flux_1 * Xnorm(Latt%a1_p) * latt%bZ1_p(:)  /  Xnorm(Latt%L1_p) + &
             &      Flux_2 * Xnorm(Latt%a2_p) * latt%bZ2_p(:)  /  Xnorm(Latt%L2_p)

        if (Bulk) then
           ! Bulk (distributed) twist: phase = A_p · del_p (line integral along the bond).
           Z_hop = Z_hop * exp(cmplx(0.d0,Iscalar(A_p,del_p),Kind(0.d0)))
        else
           ! Boundary twist: phase = A_p · XB_p (only non-zero when the hop crosses the boundary).
           Z_hop = Z_hop * exp(cmplx(0.d0,Iscalar(A_p,XB_p ),Kind(0.d0)))
        endif

        ! ---------- Contribution 2: orbital magnetic field in the Landau gauge ----------
        ! B = N_Phi / V is the uniform magnetic field (N_Phi flux quanta, V = cell area).
        ! The Landau-gauge phase for a hop (xi_p -> xj_p) is -2*pi*B * del_p(1) * (xj_p(2)+xi_p(2))/2,
        ! i.e. a symmetric midpoint convention along the y-coordinate (Peierls substitution).
        Zero =  1.0E-8
        V  =  abs(Latt%L1_p(1) * Latt%L2_p(2)  -  Latt%L1_p(2) * Latt%L2_p(1) )
        If ( V > Zero )  then
           B = real(N_Phi,kind(0.d0))/V
           Z_hop = Z_hop*exp(cmplx(0.d0, -2.d0*pi* B * del_p(1) *  ( xj_p(2) + xi_p(2) )/2.d0,kind(0.d0) ) )
           ! ---------- Contribution 3: boundary-crossing correction for the Landau gauge ----------
           ! When the hop wraps around the periodic boundary (N1 or N2 /= 0), the Landau-gauge
           ! vector potential is not single-valued; Chi() computes the resulting correction phase.
           x_p   =  Real(N2,Kind(0.d0))*Latt%L2_p
           x1_p  =  Xjp_p + Real(N1,Kind(0.d0))*Latt%L1_p
           Z_hop =  Z_hop  *  exp(cmplx( 0.d0, -Chi(x_p, x1_p,B,pi),kind(0.d0)))
           x_p   =  Real(N1,Kind(0.d0))*Latt%L1_p
           x1_p  =  Xjp_p
           Z_hop =  Z_hop  *  exp(cmplx( 0.d0, -Chi(x_p, x1_p,B,pi),kind(0.d0)))
        endif

        Generic_hopping =  Z_hop

      end function GENERIC_HOPPING

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Periodic boundary conditions for Landau gauge: c_{i+L} = e{-i Chi(L,i)} c_{i}
!>
!--------------------------------------------------------------------
      Real (Kind=kind(0.d0)) function Chi(L_p,X_p,B,pi)
        Implicit none

        Real (Kind=Kind(0.d0)), Intent(In) :: L_p(2), X_p(2), B, pi

        Chi =  - 2.d0 * pi *B * L_p(2) * X_p(1)
      end function Chi

    end Module Predefined_Hoppings
