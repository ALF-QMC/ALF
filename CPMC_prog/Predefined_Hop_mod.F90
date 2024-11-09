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
   
   subroutine set_hopping_parameters_square_anisotropic(this, Ham_tx_vec, ham_ty_vec, Ham_Chem_vec, &
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
      integer :: nf, N_Bonds, nc, I, I1
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, Ham_T_max
      real(Kind=kind(0.d0)), allocatable :: Ham_T_perp_vec(:)

      allocate (this(N_FL))

      Ham_T_max = 0.d0
      do nf = 1, N_FL
         if (abs(Ham_Tx_vec(nf)) > Ham_T_max) Ham_T_max = abs(Ham_Tx_vec(nf))
      end do

      do nf = 1, N_FL
         this(nf)%N_bonds = 2

         allocate (this(nf)%List(this(nf)%N_bonds, 4), &
              &    this(nf)%T(this(nf)%N_bonds))
         nc = 0

         nc = nc + 1
         this(nf)%T(nc) = cmplx(-Ham_tx_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 1
         this(nf)%List(nc, 2) = 1
         this(nf)%List(nc, 3) = 1
         this(nf)%List(nc, 4) = 0
         nc = nc + 1

         this(nf)%T(nc) = cmplx(-Ham_ty_vec(nf), 0.d0, kind(0.d0))
         this(nf)%List(nc, 1) = 1
         this(nf)%List(nc, 2) = 1
         this(nf)%List(nc, 3) = 0
         this(nf)%List(nc, 4) = 1

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
      allocate (this(1)%L_Fam(this(1)%N_FAM), this(1)%Prop_Fam(this(1)%N_FAM))
      this(1)%L_FAM = Latt%N/2
      this(1)%Prop_Fam = 1.d0
      allocate (this(1)%List_Fam(this(1)%N_FAM, this(1)%L_Fam(1), 2))
      this(1)%L_FAM = 0
      do I = 1, Latt%N
         if (mod(Latt%List(I, 1) + Latt%List(I, 2), 2) == 0) then
            Nf = 1
            this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 1) = I ! Unit cell
            this(1)%List_Fam(Nf, this(1)%L_Fam(Nf), 2) = 1 ! The bond (See above)
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

   end subroutine set_hopping_parameters_square_anisotropic

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
