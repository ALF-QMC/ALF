!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module provides a set of predefined hoppings as well as a
!> general framework to specify the hopping matrix for translation invariant
!> multi-orbital systems.
!>
!--------------------------------------------------------------------

module Predefined_Hoppings

   use runtime_error_mod
   use Lattices_v3
   use Operator_mod
   use WaveFunction_mod
   use MyMats
   use iso_fortran_env, only: output_unit, error_unit
   implicit none

   type Hopping_Matrix_type
      integer                   :: N_bonds
      complex(Kind=kind(0.d0)), pointer :: T(:)    !  This does not include  local terms.
      complex(Kind=kind(0.d0)), pointer :: T_loc(:)    !  This is just for the local matrix elements such as chemical potential.
      integer, pointer :: list(:, :)
      ! T(N_b=1..N_bonds)
      ! List(N_b,1) = no_1
      ! List(N_b,2) = no_2
      ! List(N_b,3) = n_1
      ! List(N_b,4) = n_2
      ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)
      integer                   :: N_Phi
      real(Kind=kind(0.d0)) :: Phi_X, Phi_Y
      logical                   :: Bulk
      ! N_Phi         = #  of flux quanta  piercieng the lattice
      ! Phi_X, Phi_Y  =  Twist
      ! Bulk          =  Twist as boundary condtion (Bulk=.F.)
      !               =  Twist as vector potential  (Bulk=.T.)

      ! For Checkerboard decomposition
      integer                            :: N_Fam
      integer, pointer :: L_Fam(:), List_Fam(:, :, :)
      real(Kind=kind(0.d0)), pointer :: Prop_Fam(:)

      integer, private, allocatable :: Multiplicity(:) !> Numer of times a given orbital occurs in the list of bonds, automatically computed
   end type Hopping_Matrix_Type

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
   subroutine Predefined_hoppings_clear(this)

      implicit none

      type(Hopping_Matrix_type), allocatable    :: this(:)
      integer :: n
      if (allocated(this)) then
         do n = 1, size(This, 1)
            deallocate (this(n)%T, this(n)%T_loc, this(n)%list)
         end do
         deallocate (this(1)%L_Fam, this(1)%List_Fam, this(1)%Prop_Fam)
         if (allocated(this(1)%Multiplicity)) deallocate (this(1)%Multiplicity)
      end if

   end subroutine Predefined_hoppings_clear
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
   integer function inquire_hop(this)

      implicit none
      type(Hopping_Matrix_type), intent(In)  :: this(:)

      real(Kind=kind(0.d0)) :: Xmax_loc, Xmax_hop, Zero, X
      integer :: nc, nf

      Zero = 1.d-10
      Xmax_loc = 0.d0
      Xmax_hop = 0.d0

      do nf = 1, size(this, 1)
         do nc = 1, size(this(1)%T_Loc, 1)
            X = sqrt(real(this(nf)%T_Loc(nc)*conjg(this(nf)%T_Loc(nc)), kind(0.d0)))
            if (X > Xmax_loc) Xmax_loc = X
         end do
         do nc = 1, this(1)%N_bonds
            X = sqrt(real(this(nf)%T(nc)*conjg(this(nf)%T(nc)), kind(0.d0)))
            if (X > Xmax_hop) Xmax_hop = X
         end do
      end do

      if (Xmax_loc < Zero .and. Xmax_hop < Zero) then
         inquire_hop = 0     !  Zero
      elseif (Xmax_loc > Zero .and. Xmax_hop < Zero) then
         inquire_hop = 1     !  Diagonal
      else
         inquire_hop = 2     !  Full
      end if

   end function inquire_hop

!--------------------------------------------------------------------
!> @author
!> Francesco Parisen Toldin
!>
!> @brief
!> Check the consistency of the checkerboard decomposition.
!> Following tests are done:
!>   -it checks that the allocated size of second index of List_Fam is
!>    at least maximum size of families. If the allocated size exceeds
!>    the required size, it issues a warning.
!>   -it checks that all bonds enter once and only once in the decomposition
!>   -it checks that every pair of bonds in each family commute
!
!--------------------------------------------------------------------
   logical function test_checkerboard_decomposition(this, Latt, inv_list)
      implicit none

      type(Hopping_Matrix_type), intent(IN) :: this
      type(Lattice), intent(IN)             :: Latt
      integer, intent(IN)                   :: inv_list(:, :)

      ! Local variables
      logical, allocatable :: all_bonds(:, :)
      integer              :: maxl, i, j, n1, n2, unit1, bond1, site1a, site1b, unit2, bond2, site2a, site2b

      test_checkerboard_decomposition = .true.
      allocate (all_bonds(Latt%N, this%N_bonds))
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
         write (error_unit, *) 'Warning: the maximum family length is ', maxl, ' allocated size is ', size(this%List_Fam, 2)
      end if

      ! Check duplicates
      do i = 1, this%N_Fam
         do j = 1, this%L_Fam(i)
            if (all_bonds(this%List_Fam(i, j, 1), this%List_Fam(i, j, 2))) then
               write (error_unit, *) 'Error in decomposition: bond at List_Fam(', i, ' ', j, ') is present twice'
               test_checkerboard_decomposition = .false.
            else
               all_bonds(this%List_Fam(i, j, 1), this%List_Fam(i, j, 2)) = .true.
            end if
         end do
      end do

      ! Check that all bonds are present in the decomposition
      do i = 1, Latt%N
         do j = 1, this%N_bonds
            if (.not. (all_bonds(i, j))) then
               write (error_unit, *) 'Error: bonds at Nunit_cell = ', i, ' bond no. ', j, ' is missing'
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
                  write (error_unit, *) 'Error: non-communting hoppings at family ', i, ' n1 = ', n1, ' List_Fam(i, n1) = ', &
                       & this%List_Fam(i, n1, 1), ' ', this%List_Fam(i, n1, 2), ' site1a = ', site1a, ' site1b = ', site1b, &
                       & ' n2 = ', n2, ' List_Fam(i, n2) = ', &
                       & this%List_Fam(i, n2, 1), ' ', this%List_Fam(i, n2, 2), ' site2a = ', site2a, ' site2b = ', site2b
                  test_checkerboard_decomposition = .false.
               end if
            end do
         end do
      end do

      deallocate (all_bonds)
   end function test_checkerboard_decomposition

   subroutine set_hopping_parameters_n_ladder_anisotropic(this, Ham_tx_vec, ham_ty_vec, Ham_Chem_vec, &
           & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable          :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:) :: Ham_tx_vec, Ham_ty_vec, Ham_Chem_vec
      real(Kind=kind(0.d0)), intent(IN), dimension(:) :: Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)               :: N_Phi_vec
      integer, intent(IN)                             :: N_FL
      logical, intent(IN)                             :: Bulk
      integer, intent(IN), dimension(:, :)            :: List, Invlist
      type(Lattice), intent(in)   :: Latt
      type(Unit_cell), intent(in) :: Latt_unit

      ! Local
      integer :: n, nf, N_Bonds, nc, I, I1, no
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Ham_T_max
      real(Kind=kind(0.d0)), allocatable :: Ham_T_perp_vec(:)

      allocate (this(N_FL))

      Ham_T_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_Tx_vec(nf)) > Ham_T_max) Ham_T_max = abs(Ham_Tx_vec(nf))
      end do

      do nf = 1, N_FL
         this(nf)%N_bonds = latt_unit%norb + (latt_unit%norb - 1)

         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
         nc = 0

         do n = 1, Latt_unit%Norb
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-ham_tx_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = n
            this(nf)%List(nc, 2) = n
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0
         end do

         do n = 1, Latt_unit%Norb - 1
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-ham_ty_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = n
            this(nf)%List(nc, 2) = n + 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0
         end do

         allocate (this(nf)%T_Loc(Latt_Unit%Norb))
         do nc = 1, Latt_Unit%Norb
            this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
         end do
         this(nf)%N_Phi = N_Phi_vec(nf)
         this(nf)%Phi_X = Phi_X_vec(nf)
         this(nf)%Phi_Y = Phi_Y_vec(nf)
         this(nf)%Bulk = Bulk
      end do

      !Set Checkerboard
      if (Latt_Unit%Norb == 1) then
         this(1)%N_FAM = 2
      elseif (Latt_Unit%Norb == 2) then
         this(1)%N_FAM = 3
      else
         this(1)%N_FAM = 4
      end if
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%L_Fam = Latt%N*Latt_unit%Norb/2
      this(1)%Prop_Fam = 1.d0
      allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))

      this(1)%L_FAM = 0
      do I = 1, Latt%N
         if (mod(Latt%List(I, 1), 2) == 0) then
            Nf = 1
            do no = 1, Latt_unit%Norb
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no ! The bond (See above)
            end do
         else
            Nf = 2
            do no = 1, Latt_unit%Norb
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no
            end do
         end if
      end do
      do no = 1, Latt_unit%Norb - 1
         if (mod(no, 2) == 1) then
            Nf = 3
            !Write(6,*)  NF, no + Latt_unit%Norb
            do I = 1, Latt%N
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no + Latt_unit%Norb
            end do
         else
            Nf = 4
            !Write(6,*)  NF, no + Latt_unit%Norb
            do I = 1, Latt%N
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no + Latt_unit%Norb
            end do
         end if
      end do

   end subroutine set_hopping_parameters_n_ladder_anisotropic

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for the square lattice.  Ham_T is the nearest
!> neighbour hopping and Ham_Chem the chemical potential.
!>
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_square(this, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
     &                                           List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable     :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:)   :: Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                  :: N_Phi_vec
      integer, intent(IN)                               :: N_FL
      logical, intent(IN)                               :: Bulk
      integer, intent(IN), dimension(:, :)               :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, I1
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Ham_T_max, Ham_T2_max
      real(Kind=kind(0.d0)), allocatable :: Ham_T_perp_vec(:)

      if (Xnorm(Latt%L2_p - Latt%a2_p) < Zero) then
         allocate (Ham_T_perp_vec(N_FL))
         Ham_T_perp_vec = 0.d0
         call Set_Default_hopping_parameters_N_Leg_Ladder(this, Ham_T_vec, Ham_T_perp_vec, Ham_Chem_vec, Phi_X_vec, &
              &                                           Phi_Y_vec, Bulk, N_Phi_vec, N_FL, &
              &                                           List, Invlist, Latt, Latt_unit)
         deallocate (Ham_T_perp_vec)
      else
         if (mod(nint(latt%L1_p(1)), 2) /= 0 .or. mod(nint(latt%L2_p(2)), 2) /= 0) then
            write (error_unit, *) '*** For  the  square  lattice,  our  implementation of the checkerborad '
            write (error_unit, *) 'decomposition  requires even  values of L_1  and L_2  ***'
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if
         allocate (this(N_FL))

         Ham_T_max = 0.d0
         Ham_T2_max = 0.d0
         do nf = 1, N_FL
            if (abs(Ham_T_vec(nf)) > Ham_T_max) Ham_T_max = abs(Ham_T_vec(nf))
            if (abs(Ham_T2_vec(nf)) > Ham_T2_max) Ham_T2_max = abs(Ham_T2_vec(nf))
         end do

         do nf = 1, N_FL
            N_bonds = 0
            N_bonds = N_bonds + 2
            if (abs(Ham_T2_max) > Zero) N_bonds = N_bonds + 2
            this(nf)%N_bonds = N_bonds
            if (this(nf)%N_bonds .gt. 0) then
               allocate (this(nf)%List(this(nf)%N_bonds, 4), &
                    &    this(nf)%T(this(nf)%N_bonds))
               nc = 0

               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1
               this(nf)%List(nc, 2) = 1
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 1

               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1
               this(nf)%List(nc, 2) = 1
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0

               if (abs(Ham_T2_max) > Zero) then
                  nc = nc + 1
                  this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
                  this(nf)%List(nc, 1) = 1
                  this(nf)%List(nc, 2) = 1
                  this(nf)%List(nc, 3) = 1
                  this(nf)%List(nc, 4) = 1

                  nc = nc + 1
                  this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
                  this(nf)%List(nc, 1) = 1
                  this(nf)%List(nc, 2) = 1
                  this(nf)%List(nc, 3) = -1
                  this(nf)%List(nc, 4) = 1
               end if

            end if
            allocate (this(nf)%T_Loc(Latt_Unit%Norb))
            do nc = 1, Latt_Unit%Norb
               this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
            end do
            this(nf)%N_Phi = N_Phi_vec(nf)
            this(nf)%Phi_X = Phi_X_vec(nf)
            this(nf)%Phi_Y = Phi_Y_vec(nf)
            this(nf)%Bulk = Bulk
         end do

         !Set Checkerboard
         this(1)%N_FAM = 4
         if (abs(Ham_T2_max) > Zero) this(1)%N_FAM = 8

         allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
         this(1)%L_FAM = Latt%N/2
         this(1)%Prop_Fam = 1.d0
         allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))

         this(1)%L_FAM = 0
         do I = 1, Latt%N
            if (mod(Latt%List(I, 1) + Latt%List(I, 2), 2) == 0) then
               Nf = 1
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1
               Nf = 2
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2
            else
               Nf = 3
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1
               Nf = 4
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2
            end if
         end do

         if (abs(Ham_T2_Max) > Zero) then
            do I = 1, Latt%N
               if ((mod(Latt%List(I, 1) + Latt%List(I, 2), 2) == 0) .and. &
                   & (mod(Latt%List(I, 1), 2) == 0)) then
                  Nf = 5
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3

                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%nnlist(I, 1, 0)
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4
               elseif ((mod(Latt%List(I, 1) + Latt%List(I, 2), 2) == 0) .and. &
                   & (mod(Latt%List(I, 1), 2) .ne. 0)) then
                  Nf = 6
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3

                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%nnlist(I, 1, 0)
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4
               elseif ((mod(Latt%List(I, 1) + Latt%List(I, 2), 2) .ne. 0) .and. &
                   & (mod(Latt%List(I, 1), 2) == 0)) then
                  Nf = 7
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3

                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%nnlist(I, 1, 0)
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4
               elseif ((mod(Latt%List(I, 1) + Latt%List(I, 2), 2) .ne. 0) .and. &
                   & (mod(Latt%List(I, 1), 2) .ne. 0)) then
                  Nf = 8
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3

                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%nnlist(I, 1, 0)
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4
               end if
            end do
         end if

      end if
   end subroutine Set_Default_hopping_parameters_square

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for n-leg-ladder.  Ham_T is the nearest neighbour hopping along the chain,  Ham_T_perp  the
!> interrung hopping and H_chem the chemical potential.
!>
!
!--------------------------------------------------------------------
   subroutine Set_Default_hopping_parameters_N_Leg_Ladder  &
        &               (this, Ham_T_vec, Ham_T_perp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk, &
        &                N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable :: this(:)
      integer, intent(IN)                   :: N_FL
      real(Kind=kind(0.d0)), intent(IN), dimension(:)    :: Ham_T_vec, Ham_T_perp_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                   :: N_Phi_vec
      logical, intent(IN)                                 :: Bulk
      integer, intent(IN), dimension(:, :)   :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, I1, n, no
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8

      select case (Latt%N)
      case (1)   !  Here the length of the  N_leg_ladder is unity such  that it
         !  effectivley maps onto a one-dimensional chain with open boundary conditions.

         allocate (this(N_FL))
         do nf = 1, N_FL
            this(nf)%N_bonds = Latt_unit%Norb - 1
            allocate (this(nf)%List(this(nf)%N_bonds, 4), this(nf)%T(this(nf)%N_bonds))
            nc = 0
            do n = 1, Latt_unit%Norb - 1
               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T_perp_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = n
               this(nf)%List(nc, 2) = n + 1
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0
            end do

            allocate (this(nf)%T_Loc(Latt_Unit%Norb))
            do nc = 1, Latt_Unit%Norb
               this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
            end do
            this(nf)%N_Phi = N_Phi_vec(nf)
            this(nf)%Phi_X = Phi_X_vec(nf)
            this(nf)%Phi_Y = Phi_Y_vec(nf)
            this(nf)%Bulk = Bulk
         end do

         ! Set Checkerboard
         if (Latt_Unit%Norb <= 2) then
            this(1)%N_FAM = 1
         else
            this(1)%N_FAM = 2
         end if
         allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
         this(1)%L_Fam = Latt_unit%Norb/2
         this(1)%Prop_Fam = 1.d0
         allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))

         this(1)%L_FAM = 0
         do no = 1, Latt_unit%Norb - 1
            if (mod(no, 2) == 1) then
               Nf = 1
               !Write(6,*)  NF, no
               do I = 1, Latt%N
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no
               end do
            else
               Nf = 2
               !Write(6,*)  NF, no
               do I = 1, Latt%N
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no
               end do
            end if
         end do

      case default
         if (mod(nint(latt%L1_p(1)), 2) /= 0) then
            write (error_unit, *) '*** For  the N_leg_ladder  lattice,  our  implementation of the checkerborad '
            write (error_unit, *) 'decomposition  requires L_1 = 1 or  L_1   even ***'
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if

         !Write(6,*) Ham_T_vec,  Ham_T_perp_vec, Ham_chem_vec
         allocate (this(N_FL))
         do nf = 1, N_FL
            this(nf)%N_bonds = Latt_unit%Norb + (Latt_unit%Norb - 1)
            allocate (this(nf)%List(this(nf)%N_bonds, 4), &
                 &    this(nf)%T(this(nf)%N_bonds))
            nc = 0
            do n = 1, Latt_unit%Norb
               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = n
               this(nf)%List(nc, 2) = n
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end do

            do n = 1, Latt_unit%Norb - 1
               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T_perp_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = n
               this(nf)%List(nc, 2) = n + 1
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0
            end do

            allocate (this(nf)%T_Loc(Latt_Unit%Norb))
            do nc = 1, Latt_Unit%Norb
               this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
            end do
            this(nf)%N_Phi = N_Phi_vec(nf)
            this(nf)%Phi_X = Phi_X_vec(nf)
            this(nf)%Phi_Y = Phi_Y_vec(nf)
            this(nf)%Bulk = Bulk
         end do

         ! Write(6,*) Latt_unit%Norb
         ! Set Checkerboard
         if (Latt_Unit%Norb == 1) then
            this(1)%N_FAM = 2
         elseif (Latt_Unit%Norb == 2) then
            this(1)%N_FAM = 3
         else
            this(1)%N_FAM = 4
         end if
         allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
         this(1)%L_Fam = Latt%N*Latt_unit%Norb/2
         this(1)%Prop_Fam = 1.d0
         allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))

         this(1)%L_FAM = 0
         do I = 1, Latt%N
            if (mod(Latt%List(I, 1), 2) == 0) then
               Nf = 1
               do no = 1, Latt_unit%Norb
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no ! The bond (See above)
               end do
            else
               Nf = 2
               do no = 1, Latt_unit%Norb
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no
               end do
            end if
         end do
         do no = 1, Latt_unit%Norb - 1
            if (mod(no, 2) == 1) then
               Nf = 3
               !Write(6,*)  NF, no + Latt_unit%Norb
               do I = 1, Latt%N
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no + Latt_unit%Norb
               end do
            else
               Nf = 4
               !Write(6,*)  NF, no + Latt_unit%Norb
               do I = 1, Latt%N
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = no + Latt_unit%Norb
               end do
            end if
         end do
      end select

   end subroutine Set_Default_hopping_parameters_N_Leg_Ladder

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for Honeycomb lattice.  Ham_T is the nearest neighbour hopping along the chain, Lambda is  the
!> Kane-Mele term and H_chem the chemical potential.  **Note**  The Kane-Mele term is not yet implemented.
!>
!
!--------------------------------------------------------------------
   subroutine Set_Default_hopping_parameters_honeycomb(this, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, &
        &                                              Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
        &                                              List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable            :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:), allocatable  :: Ham_T_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec, Ham_Lambda_vec
      integer, intent(IN), dimension(:), allocatable                 :: N_Phi_vec
      integer, intent(IN)                                           :: N_FL
      logical, intent(IN)                                           :: Bulk
      integer, intent(IN), dimension(:, :)   :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, I1, n, no
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Ham_Lambda_Max

      !Write(6,*) Ham_T_vec, Ham_Chem_vec
      Ham_Lambda_Max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_Lambda_vec(nf)) > Ham_Lambda_Max) Ham_Lambda_Max = abs(Ham_Lambda_vec(nf))
      end do
      if (abs(Ham_Lambda_max) > 0) then
         write (error_unit, *) 'Kane Mele term is not yet implemented'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      allocate (this(N_FL))
      do nf = 1, N_FL
         this(nf)%N_bonds = 3
         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
         nc = 0
         nc = nc + 1
         this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 1
         this(nf)%List(nc, 2) = 2
         this(nf)%List(nc, 3) = 0
         this(nf)%List(nc, 4) = 0

         nc = nc + 1
         this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 2
         this(nf)%List(nc, 2) = 1
         this(nf)%List(nc, 3) = 0
         this(nf)%List(nc, 4) = 1

         nc = nc + 1
         this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 1
         this(nf)%List(nc, 2) = 2
         this(nf)%List(nc, 3) = 1
         this(nf)%List(nc, 4) = -1

         allocate (this(nf)%T_Loc(Latt_Unit%Norb))
         do nc = 1, Latt_Unit%Norb
            this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
         end do
         this(nf)%N_Phi = N_Phi_vec(nf)
         this(nf)%Phi_X = Phi_X_vec(nf)
         this(nf)%Phi_Y = Phi_Y_vec(nf)
         this(nf)%Bulk = Bulk
      end do

      ! Set Checkerboard
      this(1)%N_FAM = 3
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%L_FAM = Latt%N
      this(1)%Prop_Fam = 1.d0
      allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))
      do I = 1, Latt%N
         do nf = 1, this(1)%N_FAM
            this(1)%List_Fam(nf, I, 1) = I  ! Unit cell
            this(1)%List_Fam(nf, I, 2) = nf ! The bond (See above)
         end do
      end do

   end subroutine Set_Default_hopping_parameters_honeycomb

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for a bilayer square. Ham_T1, Ham_T2, are the nearest neighbour hopping on the first and second layer and
!> Ham_T_perp   is the interlayer  hopping.
!>
!
!--------------------------------------------------------------------
   subroutine Set_Default_hopping_parameters_Bilayer_square(this, Ham_T1_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, &
        &                                                   Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
        &                                                   List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable            :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:)  :: Ham_T1_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                 :: N_Phi_vec
      logical, intent(IN)                   :: Bulk
      integer, intent(IN)                   :: N_FL
      integer, intent(IN), dimension(:, :)   :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, I1, No_Shift, n, nb
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8
      logical :: Test = .false.
      real(Kind=kind(0.d0))                :: Ham_T1_max, Ham_T2_max, Ham_Tperp_max

      Ham_T1_max = 0.d0
      Ham_T2_max = 0.d0
      Ham_Tperp_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_T1_vec(nf)) > Ham_T1_max) Ham_T1_max = abs(Ham_T1_vec(nf))
         if (abs(Ham_T2_vec(nf)) > Ham_T2_max) Ham_T2_max = abs(Ham_T2_vec(nf))
         if (abs(Ham_Tperp_vec(nf)) > Ham_Tperp_max) Ham_Tperp_max = abs(Ham_Tperp_vec(nf))
      end do

      if (nint(Latt%L2_p(2)) == 1) then
         if (mod(nint(latt%L1_p(1)), 2) /= 0) then
            write (error_unit, *) '*** For  the Bilayer square lattice,  our  implementation of the checkerborad '
            write (error_unit, *) 'decomposition  requires L_1  to be  even ***'
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if
         allocate (this(N_FL))
         do nf = 1, N_FL
            N_bonds = 0
            N_bonds = N_bonds + 1
            if (abs(Ham_Tperp_max) > Zero) N_bonds = N_bonds + 1
            if (abs(Ham_T2_max) > Zero) N_bonds = N_bonds + 1
            this(nf)%N_bonds = N_bonds
            allocate (this(nf)%List(this(nf)%N_bonds, 4), &
                 &    this(nf)%T(this(nf)%N_bonds))
            nc = 0
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T1_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            if (abs(Ham_Tperp_max) > Zero) then
               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_Tperp_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1
               this(nf)%List(nc, 2) = 2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0
            end if

            if (abs(Ham_T2_max) > Zero) then
               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2
               this(nf)%List(nc, 2) = 2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if

            allocate (this(nf)%T_Loc(Latt_Unit%Norb))
            do nc = 1, Latt_Unit%Norb
               this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
            end do
            if (abs(Ham_T2_max) < Zero .and. abs(Ham_Tperp_max) < Zero) this(nf)%T_Loc(2) = cmplx(0.0, 0.d0, kind(0.d0))
            this(nf)%N_Phi = N_Phi_vec(nf)
            this(nf)%Phi_X = Phi_X_vec(nf)
            this(nf)%Phi_Y = Phi_Y_vec(nf)
            this(nf)%Bulk = Bulk
         end do

         ! Set Checkerboard
         this(1)%N_FAM = 2
         if (abs(Ham_Tperp_max) > Zero) this(1)%N_FAM = 3

         allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
         this(1)%Prop_Fam = 1.d0

         No_Shift = 0
         if (abs(Ham_Tperp_max) > Zero) No_Shift = 1
         if (abs(Ham_T2_max) < Zero .and. abs(Ham_Tperp_max) < Zero) then
            this(1)%L_FAM = Latt%N/2
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N/2, 2))
         elseif (abs(Ham_T2_max) < Zero .and. abs(Ham_Tperp_max) > Zero) then
            this(1)%L_FAM = Latt%N/2
            this(1)%L_Fam(3) = Latt%N
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N, 2))
         elseif (abs(Ham_T2_max) > Zero .and. abs(Ham_Tperp_max) < Zero) then
            this(1)%L_FAM = Latt%N
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N, 2))
         elseif (abs(Ham_T2_max) > Zero .and. abs(Ham_Tperp_max) > Zero) then
            this(1)%L_FAM = Latt%N
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N, 2))
         end if
         this(1)%L_FAM = 0
         do I = 1, Latt%N
            if (mod(Latt%List(I, 1) + Latt%List(I, 2), 2) == 0) then
               Nf = 1
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1 ! The bond (See above)
               if (abs(Ham_T2_max) > Zero) then
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2 + No_Shift ! The bond (See above)
               end if
            else
               Nf = 2
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1
               if (abs(Ham_T2_max) > Zero) then
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2 + No_Shift ! The bond (See above)
               end if
            end if
            if (abs(Ham_Tperp_max) > Zero) then
               Nf = 3
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2
            end if
         end do
      else
         if (mod(nint(latt%L1_p(1)), 2) /= 0 .or. mod(nint(latt%L2_p(2)), 2) /= 0) then
            write (error_unit, *) '*** For  the Bilayer square lattice,  our  implementation of the checkerborad '
            write (error_unit, *) 'decomposition  requires L_1 and  L_2 to be  even ***'
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if

         allocate (this(N_FL))
         do nf = 1, N_FL
            N_bonds = 0
            N_bonds = N_bonds + 2
            if (abs(Ham_Tperp_max) > Zero) N_bonds = N_bonds + 1
            if (abs(Ham_T2_max) > Zero) N_bonds = N_bonds + 2
            this(nf)%N_bonds = N_bonds
            allocate (this(nf)%List(this(nf)%N_bonds, 4), &
                 &    this(nf)%T(this(nf)%N_bonds))
            nc = 0
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T1_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T1_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            if (abs(Ham_Tperp_max) > Zero) then
               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_Tperp_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1
               this(nf)%List(nc, 2) = 2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0
            end if

            if (abs(Ham_T2_max) > Zero) then
               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2
               this(nf)%List(nc, 2) = 2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 1

               nc = nc + 1
               this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2
               this(nf)%List(nc, 2) = 2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if

            allocate (this(nf)%T_Loc(Latt_Unit%Norb))
            do nc = 1, Latt_Unit%Norb
               this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
            end do
            if (abs(Ham_T2_max) < Zero .and. abs(Ham_Tperp_max) < Zero) this(nf)%T_Loc(2) = cmplx(0.0, 0.d0, kind(0.d0))
            this(nf)%N_Phi = N_Phi_vec(nf)
            this(nf)%Phi_X = Phi_X_vec(nf)
            this(nf)%Phi_Y = Phi_Y_vec(nf)
            this(nf)%Bulk = Bulk
         end do

         ! Set Checkerboard
         this(1)%N_FAM = 4
         if (abs(Ham_Tperp_max) > Zero) this(1)%N_FAM = 5

         allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
         this(1)%Prop_Fam = 1.d0

         No_Shift = 0
         if (abs(Ham_Tperp_max) > Zero) No_Shift = 1

         if (abs(Ham_T2_max) < Zero .and. abs(Ham_Tperp_max) < Zero) then
            this(1)%L_FAM = Latt%N/2
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N/2, 2))
         elseif (abs(Ham_T2_max) < Zero .and. abs(Ham_Tperp_max) > Zero) then
            this(1)%L_FAM = Latt%N/2
            this(1)%L_Fam(5) = Latt%N
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N, 2))
         elseif (abs(Ham_T2_max) > Zero .and. abs(Ham_Tperp_max) < Zero) then
            this(1)%L_FAM = Latt%N
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N, 2))
         elseif (abs(Ham_T2_max) > Zero .and. abs(Ham_Tperp_max) > Zero) then
            this(1)%L_FAM = Latt%N
            allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N, 2))
            No_Shift = 1
         end if
         this(1)%L_FAM = 0
         do I = 1, Latt%N
            if (mod(Latt%List(I, 1) + Latt%List(I, 2), 2) == 0) then
               Nf = 1
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1 ! The bond (See above)
               if (abs(Ham_T2_max) > Zero) then
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3 + No_Shift ! The bond (See above)
               end if
               Nf = 2
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2
               if (abs(Ham_T2_max) > Zero) then
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4 + No_Shift ! The bond (See above)
               end if
            else
               Nf = 3
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1
               if (abs(Ham_T2_max) > Zero) then
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3 + No_Shift ! The bond (See above)
               end if
               Nf = 4
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2
               if (abs(Ham_T2_max) > Zero) then
                  this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
                  this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4 + No_Shift ! The bond (See above)
               end if
            end if
            if (abs(Ham_Tperp_max) > Zero) then
               Nf = 5
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3
            end if
         end do
      end if
      ! Test
      if (Test) then
         write (6, *) this(1)%N_FAM, this(1)%L_FAM
         write (6, *) Ham_T1_max, Ham_T2_max, Ham_Tperp_max
         do nf = 1, this(1)%N_FAM
            do n = 1, this(1)%L_Fam(nf)
               I = this(1)%List_Fam(Nf, n, 1)
               nb = this(1)%List_Fam(Nf, n, 2)
    Write(6,"(I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,F6.3)")   Latt%list(I,1), Latt%list(I,2), this(1)%List(nb,1),this(1)%List(nb,2), &
                               &this(1)%List(nb, 3), this(1)%List(nb, 4), real(this(1)%T(nb))
            end do
            write (6, *)
         end do
      end if

   end subroutine Set_Default_hopping_parameters_Bilayer_square

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for a bilayer square. Ham_T1, Ham_T2, are the nearest neighbour hopping on the first and second layer and
!> Ham_T_perp   is the interlayer  hopping.
!>
!
!--------------------------------------------------------------------
   subroutine Set_Default_hopping_parameters_Bilayer_honeycomb(this, Ham_T1_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, &
        &                                                      Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL,&
        &                                                      List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable           :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:)  :: Ham_T1_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                 :: N_Phi_vec
      integer, intent(IN)                   :: N_FL
      logical, intent(IN)                   :: Bulk
      integer, intent(IN), dimension(:, :)   :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      real(Kind=kind(0.d0))                :: Ham_T1_max, Ham_T2_max, Ham_Tperp_max

      ! Local
      integer :: nf, N_Bonds, nc, I, I1, No_Shift, n, nb, no
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8
      logical :: Test = .false.

      Ham_T1_max = 0.d0
      Ham_T2_max = 0.d0
      Ham_Tperp_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_T1_vec(nf)) > Ham_T1_max) Ham_T1_max = abs(Ham_T1_vec(nf))
         if (abs(Ham_T2_vec(nf)) > Ham_T2_max) Ham_T2_max = abs(Ham_T2_vec(nf))
         if (abs(Ham_Tperp_vec(nf)) > Ham_Tperp_max) Ham_Tperp_max = abs(Ham_Tperp_vec(nf))
      end do

!!$        If (abs(Ham_T1_max) < Zero ) Then
!!$           Write(error_unit,*) 'At least Ham_T1 has to be bigger than zero'
!!$           error stop 1
!!$        endif

      allocate (this(N_FL))
      do nf = 1, N_FL
         N_bonds = 0
         N_bonds = N_bonds + 3
         if (abs(Ham_Tperp_max) > Zero) N_bonds = N_bonds + 2
         if (abs(Ham_T2_max) > Zero) N_bonds = N_bonds + 3
         this(nf)%N_bonds = N_Bonds
         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
         nc = 0
         nc = nc + 1
         this(nf)%T(nc) = cmplx(-Ham_T1_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 1
         this(nf)%List(nc, 2) = 2
         this(nf)%List(nc, 3) = 0
         this(nf)%List(nc, 4) = 0

         nc = nc + 1
         this(nf)%T(nc) = cmplx(-Ham_T1_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 2
         this(nf)%List(nc, 2) = 1
         this(nf)%List(nc, 3) = 0
         this(nf)%List(nc, 4) = 1

         nc = nc + 1
         this(nf)%T(nc) = cmplx(-Ham_T1_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 1
         this(nf)%List(nc, 2) = 2
         this(nf)%List(nc, 3) = 1
         this(nf)%List(nc, 4) = -1

         if (abs(Ham_Tperp_Max) > Zero) then
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_Tperp_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 3
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_Tperp_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 4
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0
         end if
         if (abs(Ham_T2_Max) > Zero) then
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1 + 2
            this(nf)%List(nc, 2) = 2 + 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + 2
            this(nf)%List(nc, 2) = 1 + 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1 + 2
            this(nf)%List(nc, 2) = 2 + 2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = -1
         end if
         allocate (this(nf)%T_Loc(Latt_Unit%Norb))
         do nc = 1, Latt_Unit%Norb
            this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
         end do
         if (abs(Ham_Tperp_Max) < Zero .and. abs(Ham_T2_Max) < Zero) then
            this(nf)%T_Loc(3) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%T_Loc(4) = cmplx(0.d0, 0.d0, kind(0.d0))
         end if
         this(nf)%N_Phi = N_Phi_vec(nf)
         this(nf)%Phi_X = Phi_X_vec(nf)
         this(nf)%Phi_Y = Phi_Y_vec(nf)
         this(nf)%Bulk = Bulk

      end do

      ! Set Checkerboard
      this(1)%N_FAM = 3
      if (abs(Ham_Tperp_Max) > Zero) this(1)%N_FAM = 4
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%Prop_Fam = 1.d0

      No_Shift = 0
      if (abs(Ham_Tperp_Max) > Zero) No_Shift = 2

      if (abs(Ham_T2_Max) < Zero .and. abs(Ham_Tperp_Max) < Zero) then
         this(1)%L_FAM = Latt%N
         allocate (this(1)%List_Fam(this(1)%N_FAM, Latt%N, 2))
      elseif (abs(Ham_T2_Max) < Zero .and. abs(Ham_Tperp_Max) > Zero) then
         this(1)%L_FAM = Latt%N
         this(1)%L_Fam(4) = 2*Latt%N
         allocate (this(1)%List_Fam(this(1)%N_FAM, 2*Latt%N, 2))
      elseif (abs(Ham_T2_Max) > Zero .and. abs(Ham_Tperp_Max) < Zero) then
         this(1)%L_FAM = 2*Latt%N
         allocate (this(1)%List_Fam(this(1)%N_FAM, 2*Latt%N, 2))
      elseif (abs(Ham_T2_Max) > Zero .and. abs(Ham_Tperp_Max) > Zero) then
         this(1)%L_FAM = 2*Latt%N
         allocate (this(1)%List_Fam(this(1)%N_FAM, 2*Latt%N, 2))
         No_Shift = 2
      end if

      do I = 1, Latt%N
         do nf = 1, this(1)%N_FAM
            this(1)%List_Fam(nf, I, 1) = I  ! Unit cell
            this(1)%List_Fam(nf, I, 2) = nf ! The bond (See above)
         end do
      end do
      if (abs(Ham_T2_Max) > Zero) then
         do I = 1, Latt%N
            do nf = 1, this(1)%N_FAM
               this(1)%List_Fam(nf, I + Latt%N, 1) = I                   ! Unit cell
               this(1)%List_Fam(nf, I + Latt%N, 2) = nf + 3 + No_Shift  ! The bond (See above)
            end do
         end do
      end if
      if (abs(Ham_Tperp_Max) > Zero) then
         do no = 0, 1
            do I = 1, Latt%N
               this(1)%List_Fam(4, I + no*Latt%N, 1) = I       ! Unit cell
               this(1)%List_Fam(4, I + no*Latt%N, 2) = 4 + no  ! The bond (See above)
            end do
         end do
      end if
      ! Test
      if (Test) then
         write (6, *) this(1)%N_FAM, this(1)%L_FAM
         write (6, *) Ham_T1_Max, Ham_T2_Max, Ham_Tperp_Max
         do nf = 1, this(1)%N_FAM
            do n = 1, this(1)%L_Fam(nf)
               I = this(1)%List_Fam(Nf, n, 1)
               nb = this(1)%List_Fam(Nf, n, 2)
    Write(6,"(I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,F6.3)")   Latt%list(I,1), Latt%list(I,2), this(1)%List(nb,1),this(1)%List(nb,2), &
                               &this(1)%List(nb, 3), this(1)%List(nb, 4), real(this(1)%T(nb))
            end do
            write (6, *)
         end do
      end if

   end subroutine Set_Default_hopping_parameters_Bilayer_honeycomb

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for a Pi_Flux square. Ham_T1 are the nearest neighbour hopping
!>
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_pi_flux(this, Ham_T_vec, Ham_T2_vec, Ham_T3_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
     &                                           List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable          :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:)   :: Ham_T_vec, Ham_T2_vec, Ham_T3_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                  :: N_Phi_vec
      integer, intent(IN)                               :: N_FL
      logical, intent(IN)                               :: Bulk
      integer, intent(IN), dimension(:, :)               :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, I1, amx, amy, i0, ix, iy, dx, dy, l1, l2
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Pi = acos(-1.d0), Ham_T_max, Ham_T2_max, Ham_T3_max
      real(Kind=kind(0.d0)), allocatable :: Ham_T_perp_vec(:)

      allocate (this(N_FL))

      l1 = int(Latt%L1_p(1)/Latt%a1_p(1))
      l2 = int(Latt%L2_p(2)/Latt%a2_p(2))
      amx = mod(l1, 2)
      amy = mod(l2, 2)

      Ham_T_max = 0.d0
      Ham_T2_max = 0.d0
      Ham_T3_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_T_vec(nf)) > Ham_T_max) Ham_T_max = abs(Ham_T_vec(nf))
         if (abs(Ham_T2_vec(nf)) > Ham_T2_max) Ham_T2_max = abs(Ham_T2_vec(nf))
         if (abs(Ham_T3_vec(nf)) > Ham_T3_max) Ham_T3_max = abs(Ham_T3_vec(nf))
      end do

      do nf = 1, N_FL
         this(nf)%N_bonds = 0
         if (abs(Ham_T_max) > Zero) this(nf)%N_bonds = this(nf)%N_bonds + 4
         if (abs(Ham_T2_max) > Zero) this(nf)%N_bonds = this(nf)%N_bonds + 4
         if (abs(Ham_T3_max) > Zero) this(nf)%N_bonds = this(nf)%N_bonds + 4
         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
      end do

      do nf = 1, N_FL
         nc = 0
         if (abs(Ham_T_max) > Zero) then ! hop t1
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, pi/4, kind(0.d0)))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, pi/4, kind(0.d0)))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, pi/4, kind(0.d0)))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, -pi/4, kind(0.d0)))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 1
         end if

         if (abs(Ham_T2_max) > Zero) then ! hop t2
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1
         end if

         if (abs(Ham_T3_max) > Zero) then ! hop t3
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = -1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = -1
         end if

         ! local
         allocate (this(nf)%T_Loc(Latt_Unit%Norb))
         do nc = 1, Latt_Unit%Norb
            this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
         end do
         this(nf)%N_Phi = N_Phi_vec(nf)
         this(nf)%Phi_X = Phi_X_vec(nf)
         this(nf)%Phi_Y = Phi_Y_vec(nf)
         this(nf)%Bulk = Bulk
      end do

      !Set Checkerboard
      this(1)%N_FAM = 0
      if (Ham_T_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4
      if (Ham_T2_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4 + amx + amy
      if (Ham_T3_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4 + amy + amy
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%Prop_Fam = 1.d0
      this(1)%L_FAM = Latt%N
      allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))
      this(1)%L_FAM = 0

      if (Ham_T_max > Zero) then
         do I = 1, Latt%N
            Nf = 1
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1 ! The bond (See above)

            Nf = 2
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2

            Nf = 3
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3

            Nf = 4
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4
         end do
      end if

      if (Ham_T2_max > Zero) then
         i0 = 1
         do dx = 1, L1 - amx
         do dy = 1, L2
            if (mod(latt%list(i0, 1), 2) == 0) then
               Nf = 5
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 5 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 6 ! The bond (See above)
            end if

            if (mod(latt%list(i0, 1), 2) .ne. 0) then
               Nf = 6
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 5 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 6 ! The bond (See above)
            end if
            i0 = latt%nnlist(I0, 0, 1)
         end do
         i0 = latt%nnlist(I0, 1, 0)
         end do

         do dx = 1, amx
         do dy = 1, L2

            Nf = 13
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 5 ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 6 ! The bond (See above)

            i0 = latt%nnlist(I0, 0, 1)
         end do
         i0 = latt%nnlist(I0, 1, 0)
         end do

         i0 = 1
         do dy = 1, L2 - amy
         do dx = 1, L1
            if (mod(Latt%List(i0, 2), 2) == 0) then
               Nf = 7
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 7 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 8 ! The bond (See above)
            end if

            if (mod(Latt%List(i0, 2), 2) .ne. 0) then
               Nf = 8
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 7 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 8 ! The bond (See above)
            end if
            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

         do dy = 1, amy
         do dx = 1, L1
            Nf = 14
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 7 ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 8 ! The bond (See above)

            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

      end if

      if (Ham_T3_max > Zero) then
         i0 = 1
         do dy = 1, L2 - amy
         do dx = 1, L1
            if (mod(Latt%List(i0, 2), 2) == 0) then
               Nf = 9
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 9  ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 10 ! The bond (See above)
            end if

            if (mod(Latt%List(i0, 2), 2) .ne. 0) then
               Nf = 10
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 9  ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 10 ! The bond (See above)
            end if

            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

         do dy = 1, amy
         do dx = 1, L1
            Nf = 15
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 9  ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 10 ! The bond (See above)

            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

         i0 = 1
         do dy = 1, amy
         do dx = 1, L1
            Nf = 16
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 11 ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 12 ! The bond (See above)

            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

         do dy = 1, L2 - amy
         do dx = 1, L1
            if (mod(Latt%List(i0, 2), 2) == 0) then
               Nf = 11
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 11 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 12 ! The bond (See above)
            end if

            if (mod(Latt%List(i0, 2), 2) .ne. 0) then
               Nf = 12
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 11 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 12 ! The bond (See above)
            end if

            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

      end if

   end subroutine Set_Default_hopping_parameters_pi_flux

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Default hopping for a Pi_Flux square. Ham_T1 are the nearest neighbour hopping
!>
!
!--------------------------------------------------------------------
      Subroutine Set_Default_hopping_parameters_pi_flux_ob(this, Ham_T_vec, Ham_T2_vec, Ham_T3_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
     &                                           List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable          :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:)   :: Ham_T_vec, Ham_T2_vec, Ham_T3_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                  :: N_Phi_vec
      integer, intent(IN)                               :: N_FL
      logical, intent(IN)                               :: Bulk
      integer, intent(IN), dimension(:, :)               :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, J, I1, N_bonds_tmp
      integer :: lx, ly, ix, iy, iy_no1, iy_no2, ix_no1, ix_no2, amx, amy
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Pi = acos(-1.d0), Ham_T_max, Ham_T2_max, Ham_T3_max
      real(Kind=kind(0.d0)), allocatable :: Ham_T_perp_vec(:)
      integer, allocatable :: ly_bond(:), inv_ly_bond(:, :)

      allocate (this(N_FL))

      ly = Latt_Unit%Norb/2
      lx = Latt%N

      amx = mod(lx, 2)
      amy = mod(ly, 2)

      Ham_T_max = 0.d0
      Ham_T2_max = 0.d0
      Ham_T3_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_T_vec(nf)) > Ham_T_max) Ham_T_max = abs(Ham_T_vec(nf))
         if (abs(Ham_T2_vec(nf)) > Ham_T2_max) Ham_T2_max = abs(Ham_T2_vec(nf))
         if (abs(Ham_T3_vec(nf)) > Ham_T3_max) Ham_T3_max = abs(Ham_T3_vec(nf))
      end do

      do nf = 1, N_FL
         this(nf)%N_bonds = 0
         N_bonds_tmp = 0
         if (abs(Ham_T_max) > Zero) N_bonds_tmp = N_bonds_tmp + 4
         if (abs(Ham_T2_max) > Zero) N_bonds_tmp = N_bonds_tmp + 4
         if (abs(Ham_T3_max) > Zero) N_bonds_tmp = N_bonds_tmp + 4
         this(nf)%N_bonds = N_bonds_tmp*ly
         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
      end do
      allocate (ly_bond(this(1)%N_bonds), inv_ly_bond(ly, N_bonds_tmp))

      ! y direction unit cells
      do nf = 1, N_FL
         nc = 0
         do iy = 1, ly
            if (abs(Ham_T_max) > Zero) then ! hop t1
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 1) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, pi/4, kind(0.d0)))
               this(nf)%List(nc, 1) = 1 + (iy - 1)*2
               this(nf)%List(nc, 2) = 2 + (iy - 1)*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 2) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, pi/4, kind(0.d0)))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         do iy = 1, ly - 1
            if (abs(Ham_T_max) > Zero) then ! hop t1
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 3) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, pi/4, kind(0.d0)))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + iy*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 4) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))*exp(cmplx(0.d0, -pi/4, kind(0.d0)))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + iy*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         !--------------------------------!
         ! Periodic
         !--------------------------------!
         if (abs(Ham_T_max) > Zero) then ! hop t2
            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 3) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))*exp(cmplx(0.d0,pi/4,kind(0.d0)))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + (ly - 1)*2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 4) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))*exp(cmplx(0.d0,-pi/4,kind(0.d0)))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + (ly - 1)*2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0
         end if
         !--------------------------------!

         do iy = 1, ly
            if (abs(Ham_T2_max) > Zero) then ! hop t2
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 5) = nc
               this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 6) = nc
               this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 2 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         do iy = 1, ly - 1
            if (abs(Ham_T2_max) > Zero) then ! hop t2
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 7) = nc
               this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 2 + iy*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 8) = nc
               this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + iy*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0
            end if
         end do

         !--------------------------------!
         ! Periodic
         !--------------------------------!
         if (abs(Ham_T2_max) > Zero) then ! hop t2
            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 7) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + (ly - 1)*2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 8) = nc
            !this(nf)%T(nc)    = cmplx( Ham_T2_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1 + (ly - 1)*2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0
         end if
         !--------------------------------!

         do iy = 1, ly - 1
            if (abs(Ham_T3_max) > Zero) then ! hop t3
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 9) = nc
               this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + iy*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 10) = nc
               this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 2 + iy*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         do iy = 1, ly - 1
            if (abs(Ham_T3_max) > Zero) then ! hop t3
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 11) = nc
               this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1 + iy*2
               this(nf)%List(nc, 2) = 1 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 12) = nc
               this(nf)%T(nc) = cmplx(-Ham_T3_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + iy*2
               this(nf)%List(nc, 2) = 2 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         !--------------------------------!
         ! Periodic
         !--------------------------------!
         if (abs(Ham_T3_max) > Zero) then ! hop t3
            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 9) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T3_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1 + (ly - 1)*2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 10) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T3_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + (ly - 1)*2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 11) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T3_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1 + (ly - 1)*2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 12) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T3_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 2 + (ly - 1)*2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0
         end if
         !--------------------------------!

         ! local
         allocate (this(nf)%T_Loc(Latt_Unit%Norb))
         do nc = 1, Latt_Unit%Norb
            this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
         end do
         this(nf)%N_Phi = N_Phi_vec(nf)
         this(nf)%Phi_X = Phi_X_vec(nf)
         this(nf)%Phi_Y = Phi_Y_vec(nf)
         this(nf)%Bulk = Bulk
      end do

      !Set Checkerboard
      this(1)%N_FAM = 0
      if (Ham_T_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4
      if (Ham_T2_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4 + amx + amy
      if (Ham_T3_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4 + 2*amy
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%Prop_Fam = 1.d0
      this(1)%L_FAM = lx*ly
      allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))
      this(1)%L_FAM = 0

      if (Ham_T_max > Zero) then
         do I = 1, Latt%N
         do iy = 1, ly
            Nf = 1
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 1) ! The bond (See above)

            Nf = 2
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 2)

            Nf = 3
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 3)

            Nf = 4
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 4)
         end do
         end do
      end if

      if (Ham_T2_max > Zero) then
         do I = 1, Latt%N - amx
         do iy = 1, ly
            if (mod(Latt%List(I, 1), 2) == 0) then
               Nf = 5
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 5) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 6) ! The bond (See above)
            end if

            if (mod(Latt%List(I, 1), 2) .ne. 0) then
               Nf = 6
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 5) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 6) ! The bond (See above)
            end if
         end do
         end do

         !--------------------------------!
         do I = 1, amx
         do iy = 1, ly
            Nf = 13
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%N ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 5) ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%N
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 6)
         end do
         end do
         !--------------------------------!

         do I = 1, Latt%N
         do iy = 1, ly - amy
            if (mod(iy, 2) == 0) then
               Nf = 7
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 7) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 8) ! The bond (See above)
            end if

            if (mod(iy, 2) .ne. 0) then
               Nf = 8
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 7) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 8) ! The bond (See above)
            end if
         end do
         end do

         !--------------------------------!
         do iy = 1, amy
         do I = 1, Latt%N
            Nf = 14
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 7) ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 8)
         end do
         end do
         !--------------------------------!

      end if

      if (Ham_T3_max > Zero) then
         do I = 1, Latt%N
         do iy = 1, ly - amy
            if (mod(iy, 2) == 0) then
               Nf = 9
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 9) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I  ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 10) ! The bond (See above)
            end if

            if (mod(iy, 2) .ne. 0) then
               Nf = 10
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 9) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I  ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 10) ! The bond (See above)
            end if

            if (mod(iy, 2) == 0) then
               Nf = 11
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I  ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 11) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I  ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 12) ! The bond (See above)
            end if

            if (mod(iy, 2) .ne. 0) then
               Nf = 12
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I  ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 11) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I  ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 12) ! The bond (See above)
            end if
         end do
         end do

         !--------------------------------!
         do I = 1, Latt%N
         do J = 1, amy
            Nf = 15
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 9) ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 10)

            Nf = 16
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 11)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 12)
         end do
         end do
         !--------------------------------!

      end if

      deallocate (ly_bond, inv_ly_bond)

   end subroutine Set_Default_hopping_parameters_pi_flux_ob

      Subroutine Set_Default_hopping_parameters_pi_flux_qbt(this, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
     &                                           List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable          :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:)   :: Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                  :: N_Phi_vec
      integer, intent(IN)                               :: N_FL
      logical, intent(IN)                               :: Bulk
      integer, intent(IN), dimension(:, :)               :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, I1, amx, amy, i0, ix, iy, dx, dy, l1, l2
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Pi = acos(-1.d0), Ham_T_max, Ham_T2_max, Ham_T3_max
      real(Kind=kind(0.d0)), allocatable :: Ham_T_perp_vec(:)

      allocate (this(N_FL))

      l1 = int(Latt%L1_p(1)/Latt%a1_p(1))
      l2 = int(Latt%L2_p(2)/Latt%a2_p(2))
      amx = mod(l1, 2)
      amy = mod(l2, 2)

      Ham_T_max = 0.d0
      Ham_T2_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_T_vec(nf)) > Ham_T_max) Ham_T_max = abs(Ham_T_vec(nf))
         if (abs(Ham_T2_vec(nf)) > Ham_T2_max) Ham_T2_max = abs(Ham_T2_vec(nf))
      end do

      do nf = 1, N_FL
         this(nf)%N_bonds = 0
         if (abs(Ham_T_max) > Zero) this(nf)%N_bonds = this(nf)%N_bonds + 4
         if (abs(Ham_T2_max) > Zero) this(nf)%N_bonds = this(nf)%N_bonds + 4
         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
      end do

      do nf = 1, N_FL
         nc = 0
         if (abs(Ham_T_max) > Zero) then ! hop t1
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 1
         end if

         if (abs(Ham_T2_max) > Zero) then ! hop t2
            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1

            nc = nc + 1
            this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 1
         end if

         ! local
         allocate (this(nf)%T_Loc(Latt_Unit%Norb))
         do nc = 1, Latt_Unit%Norb
            this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
         end do
         this(nf)%N_Phi = N_Phi_vec(nf)
         this(nf)%Phi_X = Phi_X_vec(nf)
         this(nf)%Phi_Y = Phi_Y_vec(nf)
         this(nf)%Bulk = Bulk
      end do

      !Set Checkerboard
      this(1)%N_FAM = 0
      if (Ham_T_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4
      if (Ham_T2_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4 + amx + amy
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%Prop_Fam = 1.d0
      this(1)%L_FAM = Latt%N
      allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))
      this(1)%L_FAM = 0

      if (Ham_T_max > Zero) then
         do I = 1, Latt%N
            Nf = 1
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1 ! The bond (See above)

            Nf = 2
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 2

            Nf = 3
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 3

            Nf = 4
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 4
         end do
      end if

      if (Ham_T2_max > Zero) then
         i0 = 1
         do dx = 1, L1 - amx
         do dy = 1, L2
            if (mod(latt%list(i0, 1), 2) == 0) then
               Nf = 5
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 5 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 6 ! The bond (See above)
            end if

            if (mod(latt%list(i0, 1), 2) .ne. 0) then
               Nf = 6
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 5 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 6 ! The bond (See above)
            end if
            i0 = latt%nnlist(I0, 0, 1)
         end do
         i0 = latt%nnlist(I0, 1, 0)
         end do

         do dx = 1, amx
         do dy = 1, L2

            Nf = 9
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 5 ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 6 ! The bond (See above)

            i0 = latt%nnlist(I0, 0, 1)
         end do
         i0 = latt%nnlist(I0, 1, 0)
         end do

         i0 = 1
         do dy = 1, L2 - amy
         do dx = 1, L1
            if (mod(Latt%List(i0, 2), 2) == 0) then
               Nf = 7
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 7 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 8 ! The bond (See above)
            end if

            if (mod(Latt%List(i0, 2), 2) .ne. 0) then
               Nf = 8
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 7 ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 8 ! The bond (See above)
            end if
            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

         do dy = 1, amy
         do dx = 1, L1
            Nf = 10
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 7 ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I0 ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 8 ! The bond (See above)

            i0 = latt%nnlist(I0, 1, 0)
         end do
         i0 = latt%nnlist(I0, 0, 1)
         end do

      end if

   end subroutine Set_Default_hopping_parameters_pi_flux_qbt

      Subroutine Set_Default_hopping_parameters_pi_flux_qbt_ob(this, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
     &                                           List, Invlist, Latt, Latt_unit)

      implicit none

      type(Hopping_Matrix_type), allocatable          :: this(:)
      real(Kind=kind(0.d0)), intent(IN), dimension(:)   :: Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_x_vec, Phi_y_vec
      integer, intent(IN), dimension(:)                  :: N_Phi_vec
      integer, intent(IN)                               :: N_FL
      logical, intent(IN)                               :: Bulk
      integer, intent(IN), dimension(:, :)               :: List, Invlist
      type(Lattice), intent(in)            :: Latt
      type(Unit_cell), intent(in)            :: Latt_unit

      ! Local
      integer :: nf, N_Bonds, nc, I, J, I1, N_bonds_tmp
      integer :: lx, ly, ix, iy, iy_no1, iy_no2, ix_no1, ix_no2, amx, amy
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Pi = acos(-1.d0), Ham_T_max, Ham_T2_max, Ham_T3_max
      real(Kind=kind(0.d0)), allocatable :: Ham_T_perp_vec(:)
      integer, allocatable :: ly_bond(:), inv_ly_bond(:, :)

      allocate (this(N_FL))

      ly = Latt_Unit%Norb/2
      lx = Latt%N

      amx = mod(lx, 2)
      amy = mod(ly, 2)

      Ham_T_max = 0.d0
      Ham_T2_max = 0.d0
      Ham_T3_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_T_vec(nf)) > Ham_T_max) Ham_T_max = abs(Ham_T_vec(nf))
         if (abs(Ham_T2_vec(nf)) > Ham_T2_max) Ham_T2_max = abs(Ham_T2_vec(nf))
      end do

      do nf = 1, N_FL
         this(nf)%N_bonds = 0
         N_bonds_tmp = 0
         if (abs(Ham_T_max) > Zero) N_bonds_tmp = N_bonds_tmp + 4
         if (abs(Ham_T2_max) > Zero) N_bonds_tmp = N_bonds_tmp + 4
         this(nf)%N_bonds = N_bonds_tmp*ly
         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
      end do
      allocate (ly_bond(this(1)%N_bonds), inv_ly_bond(ly, N_bonds_tmp))

      ! y direction unit cells
      do nf = 1, N_FL
         nc = 0
         do iy = 1, ly
            if (abs(Ham_T_max) > Zero) then ! hop t1
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 1) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1 + (iy - 1)*2
               this(nf)%List(nc, 2) = 2 + (iy - 1)*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 2) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         do iy = 1, ly - 1
            if (abs(Ham_T_max) > Zero) then ! hop t1
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 3) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + iy*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 4) = nc
               this(nf)%T(nc) = cmplx(-Ham_T_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + iy*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         !--------------------------------!
         ! Periodic
         !--------------------------------!
         if (abs(Ham_T_max) > Zero) then ! hop t2
            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 3) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + (ly - 1)*2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 4) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + (ly - 1)*2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 1
            this(nf)%List(nc, 4) = 0
         end if
         !--------------------------------!

         do iy = 1, ly
            if (abs(Ham_T2_max) > Zero) then ! hop t2
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 5) = nc
               this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 6) = nc
               this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 2 + (iy - 1)*2
               this(nf)%List(nc, 3) = 1
               this(nf)%List(nc, 4) = 0
            end if
         end do

         do iy = 1, ly - 1
            if (abs(Ham_T2_max) > Zero) then ! hop t2
               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 7) = nc
               this(nf)%T(nc) = cmplx(-Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 2 + (iy - 1)*2
               this(nf)%List(nc, 2) = 2 + iy*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0

               nc = nc + 1
               ly_bond(nc) = iy
               inv_ly_bond(iy, 8) = nc
               this(nf)%T(nc) = cmplx(Ham_T2_vec(nf), 0.d0, kind(0.d0))
               this(nf)%List(nc, 1) = 1 + (iy - 1)*2
               this(nf)%List(nc, 2) = 1 + iy*2
               this(nf)%List(nc, 3) = 0
               this(nf)%List(nc, 4) = 0
            end if
         end do

         !--------------------------------!
         ! Periodic
         !--------------------------------!
         if (abs(Ham_T2_max) > Zero) then ! hop t2
            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 7) = nc
            !this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 2 + (ly - 1)*2
            this(nf)%List(nc, 2) = 2
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0

            nc = nc + 1
            ly_bond(nc) = ly
            inv_ly_bond(ly, 8) = nc
            !this(nf)%T(nc)    = cmplx( Ham_T2_vec(nf),0.d0,kind(0.d0))
            this(nf)%T(nc) = cmplx(0.d0, 0.d0, kind(0.d0))
            this(nf)%List(nc, 1) = 1 + (ly - 1)*2
            this(nf)%List(nc, 2) = 1
            this(nf)%List(nc, 3) = 0
            this(nf)%List(nc, 4) = 0
         end if
         !--------------------------------!

         ! local
         allocate (this(nf)%T_Loc(Latt_Unit%Norb))
         do nc = 1, Latt_Unit%Norb
            this(nf)%T_Loc(nc) = cmplx(-Ham_Chem_vec(nf), 0.d0, kind(0.d0))
         end do
         this(nf)%N_Phi = N_Phi_vec(nf)
         this(nf)%Phi_X = Phi_X_vec(nf)
         this(nf)%Phi_Y = Phi_Y_vec(nf)
         this(nf)%Bulk = Bulk
      end do

      !Set Checkerboard
      this(1)%N_FAM = 0
      if (Ham_T_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4
      if (Ham_T2_max > Zero) this(1)%N_FAM = this(1)%N_FAM + 4 + amx + amy
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%Prop_Fam = 1.d0
      this(1)%L_FAM = lx*ly
      allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))
      this(1)%L_FAM = 0

      if (Ham_T_max > Zero) then
         do I = 1, Latt%N
         do iy = 1, ly
            Nf = 1
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 1) ! The bond (See above)

            Nf = 2
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 2)

            Nf = 3
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 3)

            Nf = 4
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 4)
         end do
         end do
      end if

      if (Ham_T2_max > Zero) then
         do I = 1, Latt%N - amx
         do iy = 1, ly
            if (mod(Latt%List(I, 1), 2) == 0) then
               Nf = 5
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 5) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 6) ! The bond (See above)
            end if

            if (mod(Latt%List(I, 1), 2) .ne. 0) then
               Nf = 6
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 5) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 6) ! The bond (See above)
            end if
         end do
         end do

         !--------------------------------!
         do I = 1, amx
         do iy = 1, ly
            Nf = 9
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%N ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 5) ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = Latt%N
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 6)
         end do
         end do
         !--------------------------------!

         do I = 1, Latt%N
         do iy = 1, ly - amy
            if (mod(iy, 2) == 0) then
               Nf = 7
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 7) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 8) ! The bond (See above)
            end if

            if (mod(iy, 2) .ne. 0) then
               Nf = 8
               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 7) ! The bond (See above)

               this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
               this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(iy, 8) ! The bond (See above)
            end if
         end do
         end do

         !--------------------------------!
         do iy = 1, amy
         do I = 1, Latt%N
            Nf = 10
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 7) ! The bond (See above)

            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = inv_ly_bond(ly, 8)
         end do
         end do
         !--------------------------------!

      end if

      deallocate (ly_bond, inv_ly_bond)

   end subroutine Set_Default_hopping_parameters_pi_flux_qbt_ob

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Given the checkerbord decompostion
!> the routine allocates and sets OP_T Set_Default_hopping_parameters_"Lattice" routine,
!> this routine generates the data for the symmetric decomposition.
!
!--------------------------------------------------------------------
   subroutine Symmetrize_Families(this)
      implicit none

      type(Hopping_Matrix_type), allocatable         :: this(:)
      ! In Families.  Out Symmetrized Families.

      !  Make a copy  of the unsymmetrized forms
      integer                              ::  N_FAM_C
      integer, allocatable                 ::  L_Fam_C(:), List_Fam_C(:, :, :)
      real(Kind=kind(0.d0)), allocatable  ::  Prop_Fam_C(:)

      integer :: n, n1, n2, n_f_max, n_l_max, nc
      integer, allocatable ::  list_Fam_tmp(:)

      ! Copy
      N_FAM_C = this(1)%N_FAM
      allocate (L_FAM_C(N_FAM_C))
      n2 = size(this(1)%List_Fam, 2)
      allocate (List_Fam_C(N_FAM_C, n2, 2), Prop_Fam_C(N_FAM_C))
      L_FAM_C = this(1)%L_FAM
      List_Fam_C = this(1)%List_Fam
      Prop_Fam_C = this(1)%Prop_Fam

      ! Re-allocate
      this(1)%N_FAM = 2*N_FAM_C - 1
      deallocate (this(1)%L_Fam, this(1)%List_Fam, this(1)%Prop_Fam)
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%List_Fam(this(1)%N_FAM, n2, 2), this(1)%Prop_Fam(this(1)%N_FAM))

      ! Symmetrize
      ! Find the longest family.
      n_l_max = 0
      n_f_max = 0
      do n = 1, N_FAM_C
         if (L_FAM_C(n) > n_l_max) then
            n_l_max = L_FAM_C(n)
            n_f_max = n
         end if
      end do
      !Write(6,*) 'N_f_max' , n_f_max
      allocate (list_Fam_tmp(this(1)%N_FAM))
      nc = 0
      do n = 1, N_FAM_C
         nc = nc + 1
         list_Fam_tmp(nc) = n
      end do
      do n = N_FAM_C - 1, 1, -1
         nc = nc + 1
         list_Fam_tmp(nc) = n
      end do

      ! Place the largest familly in the middle and set the time step.
      this(1)%Prop_Fam = 0.5d0
      this(1)%Prop_Fam(N_FAM_C) = 1.d0
      if (N_F_Max .ne. N_FAM_C) then
         list_Fam_tmp(N_FAM_C) = n_f_max
         list_Fam_tmp(1) = N_FAM_C
         list_Fam_tmp(this(1)%N_FAM) = N_Fam_C
      end if

      do n = 1, this(1)%N_FAM
         n1 = list_Fam_tmp(n)
         this(1)%L_Fam(n) = L_FAM_C(n1)
         this(1)%List_Fam(n, :, :) = List_Fam_C(n1, :, :)
      end do

      ! Clean
      deallocate (L_FAM_C, List_Fam_C, Prop_Fam_C, List_Fam_tmp)

      !Write(6,*)  this(1)%N_FAM
      !Write(6,*)  this(1)%L_FAM
      !Write(6,*)  this(1)%Prop_Fam

   end subroutine Symmetrize_Families

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Given the Hopping-matrix, and if required the checkerboard decomposion (i.e. private data of this module)
!> the routine allocates and sets OP_T.
!
!--------------------------------------------------------------------
   subroutine Predefined_Hoppings_set_OPT(this, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_T)

      implicit none

      type(Hopping_Matrix_type), allocatable             :: this(:)
      integer, intent(IN), dimension(:, :)                 :: List, Invlist
      type(Lattice), intent(in)                          :: Latt
      type(Unit_cell), intent(in)                          :: Latt_unit
      real(Kind=kind(0.d0)), intent(In)                  :: Dtau
      logical, intent(IN)                                 :: Checkerboard, Symm

      type(operator), intent(Out), dimension(:, :), allocatable  :: Op_T

      ! Local
      integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
      integer                           :: n_1, n_2, Nb, n_f, l_f, n_l, N, nc, orb
      real(Kind=kind(0.d0))          :: Ham_T, Ham_Chem, Phi_X, Phi_Y
      logical                           :: Bulk
      complex(Kind=kind(0.d0))          :: Z

      N_FL = size(this, 1)
      !Write(6,*)  'N_FL ', N_FL
      Ndim = Latt%N*Latt_Unit%Norb

      ! Test of correctness of checkerboard decomposition
      if (checkerboard) then
     if (.not. (test_checkerboard_decomposition(this(1), Latt, invlist))) call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      select case (inquire_hop(this))
      case (0)  !  Zero
         allocate (Op_T(1, N_FL))
         do nf = 1, N_FL
            call Op_make(Op_T(1, nf), 1)
            Op_T(1, nf)%P(1) = 1
            Op_T(1, nf)%O(1, 1) = cmplx(0.d0, 0.d0, kind(0.d0))
            Op_T(1, nf)%g = 0.d0
            Op_T(1, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
            call Op_set(Op_T(1, nf))
         end do
      case (1)  ! Diagonal
         allocate (Op_T(Ndim, N_FL))
         do nf = 1, N_FL
            do n = 1, ndim
               call Op_make(Op_T(n, nf), 1)
               Op_T(n, nf)%P(1) = n
               Op_T(n, nf)%O(1, 1) = this(nf)%T_Loc(list(n, 2))
               Op_T(n, nf)%g = -Dtau
               Op_T(n, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
               call Op_set(Op_T(n, nf))
            end do
         end do
      case default
         if (.not. Checkerboard) then
            allocate (Op_T(1, N_FL))
            do nf = 1, N_FL
               !Write(6,*)
               call Op_make(Op_T(1, nf), Ndim)   ! This is too restrictive for the Kondo type models. The hopping only occurs on one subsystem.
               N_Phi = this(nf)%N_Phi
               Phi_X = this(nf)%Phi_X
               Phi_Y = this(nf)%Phi_Y
               Bulk = this(nf)%Bulk
               !Write(6,*) N_Phi, Phi_X,Phi_Y, Bulk
               !Write(6,*) This(nf)%list
               do I = 1, Latt%N
                  do Nb = 1, this(nf)%N_bonds
                     no_I = this(nf)%list(Nb, 1)
                     no_J = this(nf)%list(Nb, 2)
                     n_1 = this(nf)%list(Nb, 3)
                     n_2 = this(nf)%list(Nb, 4)
                     J = Latt%nnlist(I, n_1, n_2)
                     Z = Generic_hopping(I, no_I, n_1, n_2, no_J, N_Phi, Phi_x, Phi_y, Bulk, Latt, Latt_Unit)
                     I1 = Invlist(I, no_I)
                     J1 = Invlist(J, no_J)
                     Op_T(1, nf)%O(I1, J1) = this(nf)%T(Nb)*Z
                     Op_T(1, nf)%O(J1, I1) = conjg(this(nf)%T(Nb)*Z)
                  end do
                  ! T(N_b=1..N_bonds)
                  ! List(N_b,1) = no_1
                  ! List(N_b,2) = no_2
                  ! List(N_b,3) = n_1
                  ! List(N_b,4) = n_2
                  ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)
                  do no_I = 1, Latt_Unit%Norb
                     I1 = Invlist(I, no_I)
                     Op_T(1, nf)%O(I1, I1) = this(nf)%T_Loc(no_I)
                  end do
               end do
               do I = 1, Ndim
                  Op_T(1, nf)%P(i) = i
               end do
               Op_T(1, nf)%g = -Dtau
               Op_T(1, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
               call Op_set(Op_T(1, nf))
               !Do I = 1,Size(Op_T(1,nf)%E,1)
               !   Write(6,*) Op_T(1,nf)%E(I)
               !Enddo
            end do
         elseif (Checkerboard) then
            if (Symm) call Symmetrize_families(this)
            N = 0
            do n_f = 1, this(1)%N_FAM
               N = N + this(1)%L_Fam(n_f)
            end do
            allocate (Op_T(N, N_FL))
            do nf = 1, N_FL
               ! Compute Multiplicity
               allocate (this(nf)%Multiplicity(Latt_Unit%Norb))
               this(nf)%Multiplicity = 0
               do i = 1, size(this(nf)%List, 1)
                  orb = this(nf)%List(i, 1)
                  this(nf)%Multiplicity(orb) = this(nf)%Multiplicity(orb) + 1
                  orb = this(nf)%List(i, 2)
                  this(nf)%Multiplicity(orb) = this(nf)%Multiplicity(orb) + 1
               end do

               N_Phi = this(nf)%N_Phi
               Phi_X = this(nf)%Phi_X
               Phi_Y = this(nf)%Phi_Y
               Bulk = this(nf)%Bulk
               do nc = 1, size(Op_T, 1)
                  call Op_make(Op_T(nc, nf), 2)
               end do
               nc = 0
               do n_f = 1, this(1)%N_FAM
                  do l_f = 1, this(1)%L_Fam(n_f)
                     I = this(1)%List_Fam(n_f, l_f, 1)
                     nb = this(1)%List_Fam(n_f, l_f, 2)
                     no_I = this(nf)%list(Nb, 1)
                     no_J = this(nf)%list(Nb, 2)
                     n_1 = this(nf)%list(Nb, 3)
                     n_2 = this(nf)%list(Nb, 4)
                     J = Latt%nnlist(I, n_1, n_2)
                     Z = Generic_hopping(I, no_I, n_1, n_2, no_J, N_Phi, Phi_x, Phi_y, Bulk, Latt, Latt_Unit)
                     I1 = Invlist(I, no_I)
                     J1 = Invlist(J, no_J)
                     nc = nc + 1
                     Op_T(nc, nf)%P(1) = I1
                     Op_T(nc, nf)%P(2) = J1
                     Op_T(nc, nf)%O(1, 2) = this(nf)%T(Nb)*Z
                     Op_T(nc, nf)%O(2, 1) = conjg(this(nf)%T(Nb)*Z)
                     Op_T(nc, nf)%O(1, 1) = this(nf)%T_loc(no_I)/this(1)%Multiplicity(no_I)
                     Op_T(nc, nf)%O(2, 2) = this(nf)%T_loc(no_J)/this(1)%Multiplicity(no_J)
                     Op_T(nc, nf)%g = -Dtau*this(1)%Prop_Fam(n_f)
                     Op_T(nc, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
                     call Op_set(Op_T(nc, nf))
                  end do
               end do
            end do
         end if
      end select

   end subroutine Predefined_Hoppings_set_OPT

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> The subroutine computes the kinetic energy based on the generic form of the
!> the hopping matrix.
!>
!--------------------------------------------------------------------
   subroutine Predefined_Hoppings_Compute_Kin(this, List, Invlist, Latt, Latt_unit, GRC, Z_Kin)

      implicit none

      type(Hopping_Matrix_type), allocatable  :: this(:)
      integer, intent(IN), dimension(:, :)                 :: List, Invlist
      type(Lattice), intent(in)                          :: Latt
      type(Unit_cell), intent(in)                          :: Latt_unit
      complex(Kind=kind(0.d0)), intent(in), dimension(:, :, :) :: GRC(:, :, :)
      complex(Kind=kind(0.d0)), intent(out) :: Z_kin

      !Local
      integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
      integer                           :: n_1, n_2, Nb, n_f, l_f, n_l, N, nc
      real(Kind=kind(0.d0))          :: Ham_T, Ham_Chem, Phi_X, Phi_Y
      logical                           :: Bulk
      complex(Kind=kind(0.d0))          :: Z

      select case (inquire_hop(this))
      case (0)  !  Zero
         Z_Kin = cmplx(0.d0, 0.d0, kind(0.d0))
      case (1)
         Z_Kin = cmplx(0.d0, 0.d0, kind(0.d0))
         N_FL = size(GRC, 3)
         do nf = 1, N_FL
            do I = 1, Latt%N
               do no_I = 1, Latt_Unit%Norb
                  I1 = Invlist(I, no_I)
                  Z_Kin = Z_Kin + this(nf)%T_Loc(no_I)*GRC(I1, I1, nf)
               end do
            end do
         end do
      case default
         N_FL = size(GRC, 3)
         Z_Kin = cmplx(0.d0, 0.d0, kind(0.d0))
         do nf = 1, N_FL
            N_Phi = this(nf)%N_Phi
            Phi_X = this(nf)%Phi_X
            Phi_Y = this(nf)%Phi_Y
            Bulk = this(nf)%Bulk
            do I = 1, Latt%N
               do Nb = 1, this(nf)%N_bonds
                  no_I = this(nf)%list(Nb, 1)
                  no_J = this(nf)%list(Nb, 2)
                  n_1 = this(nf)%list(Nb, 3)
                  n_2 = this(nf)%list(Nb, 4)
                  J = Latt%nnlist(I, n_1, n_2)
                  Z = Generic_hopping(I, no_I, n_1, n_2, no_J, N_Phi, Phi_x, Phi_y, Bulk, Latt, Latt_Unit)
                  I1 = Invlist(I, no_I)
                  J1 = Invlist(J, no_J)
                  Z_Kin = Z_Kin + this(nf)%T(Nb)*Z*GRC(I1, J1, nf) + conjg(this(nf)%T(Nb)*Z)*GRC(J1, I1, nf)
               end do
               do no_I = 1, Latt_Unit%Norb
                  I1 = Invlist(I, no_I)
                  Z_Kin = Z_Kin + this(nf)%T_Loc(no_I)*GRC(I1, I1, nf)
               end do
            end do
         end do
      end select

   end subroutine Predefined_Hoppings_Compute_Kin
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

   complex(Kind=kind(0.d0)) function Generic_hopping(i, no_i, del_1, del_2, no_j, N_Phi, Flux_1, Flux_2, Bulk, Latt, Latt_Unit)

      use Lattices_v3
      implicit none

      integer, intent(In) :: N_Phi, i, no_i, del_1, del_2, no_j
      type(Unit_cell), intent(In) :: Latt_Unit
      type(Lattice), intent(In) :: Latt
      real(Kind=kind(0.d0)), intent(In) :: Flux_1, Flux_2
      logical, intent(In) :: Bulk

      !Local
      integer                   :: j, N1, N2, n
      real(Kind=kind(0.d0))    :: xj_p(2), xi_p(2), xjp_p(2), del_p(2), A_p(2), pi, XB_p(2), V, B, Zero, x_p(2), x1_p(2)

      complex(Kind=kind(0.d0)) :: Z_hop

      Z_hop = cmplx(1.d0, 0.d0, kind(0.d0))

      xj_p = real(latt%list(i, 1) + del_1, kind(0.d0))*latt%a1_p + real(latt%list(i, 2) + del_2, kind(0.d0))*latt%a2_p
      ! Check if you have crossed the boundary:  xj_p  = xjp_p + N1*L1_p  + N2*L2_p  with  xjp_p  in the set of lattice sites.
      N1 = 0; N2 = 0
      call npbc(xjp_p, xj_p, Latt%L1_p, Latt%L2_p, N1, N2)
      XB_p = real(N1, kind(0.d0))*Latt%L1_p + real(N2, kind(0.d0))*Latt%L2_p
      do n = 1, 2
         xj_p(n) = xj_p(n) + Latt_unit%Orb_pos_p(no_j, n)
         xjp_p(n) = xjp_p(n) + Latt_unit%Orb_pos_p(no_j, n)
      end do
      xi_p = real(latt%list(i, 1), kind(0.d0))*latt%a1_p + real(latt%list(i, 2), kind(0.d0))*latt%a2_p
      do n = 1, 2
         xi_p(n) = xi_p(n) + Latt_unit%Orb_pos_p(no_i, n)
      end do
        !!Check that  xjp_p(:) + XB_p(:) =  xj_p(:)
        !!x1_p(:) = xjp_p(:) + XB_p
        !!Write(6,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,I2,3x,I2)")  x1_p(1),x1_p(2), xj_p(1), xj_p(2), N1,N2
        !! -->  i + del_1*a_1 + del_2* a_2  =  i' + N1*L1_p + N2*L2_p  with i' in the set of lattice points.

      ! The hopping increment.
      del_p = xj_p - xi_p

      !Twist
      pi = acos(-1.d0)
      A_p(:) = Flux_1*Xnorm(Latt%a1_p)*latt%bZ1_p(:)/Xnorm(Latt%L1_p) + &
           &      Flux_2*Xnorm(Latt%a2_p)*latt%bZ2_p(:)/Xnorm(Latt%L2_p)

      if (Bulk) then
         !Twist in bulk
         Z_hop = Z_hop*exp(cmplx(0.d0, Iscalar(A_p, del_p), kind(0.d0)))
      else
         !Twist as boundary
         Z_hop = Z_hop*exp(cmplx(0.d0, Iscalar(A_p, XB_p), kind(0.d0)))
      end if

      !Orbital magnetic field (Landau gauge)
      Zero = 1.0e-8
      V = abs(Latt%L1_p(1)*Latt%L2_p(2) - Latt%L1_p(2)*Latt%L2_p(1))
      if (V > Zero) then
         B = real(N_Phi, kind(0.d0))/V
         Z_hop = Z_hop*exp(cmplx(0.d0, -2.d0*pi*B*del_p(1)*(xj_p(2) + xi_p(2))/2.d0, kind(0.d0)))
         ! Boundary
         x_p = real(N2, kind(0.d0))*Latt%L2_p
         x1_p = Xjp_p + real(N1, kind(0.d0))*Latt%L1_p
         Z_hop = Z_hop*exp(cmplx(0.d0, -Chi(x_p, x1_p, B, pi), kind(0.d0)))
         x_p = real(N1, kind(0.d0))*Latt%L1_p
         x1_p = Xjp_p
         Z_hop = Z_hop*exp(cmplx(0.d0, -Chi(x_p, x1_p, B, pi), kind(0.d0)))
      end if

      Generic_hopping = Z_hop

   end function GENERIC_HOPPING

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Periodic boundary conditions for Landau gauge: c_{i+L} = e{-i Chi(L,i)} c_{i}
!>
!--------------------------------------------------------------------
   real(Kind=kind(0.d0)) function Chi(L_p, X_p, B, pi)
      implicit none

      real(Kind=kind(0.d0)), intent(In) :: L_p(2), X_p(2), B, pi

      Chi = -2.d0*pi*B*L_p(2)*X_p(1)
   end function Chi

end module Predefined_Hoppings
