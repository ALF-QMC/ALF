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
!> This module provides a set of predefined hoppings as well as a
!> general framework to specify the hopping matrix for translation invariant
!> multi-orbital systems.
!>
!--------------------------------------------------------------------


    Module Predefined_Hoppings

      Use runtime_error_mod
      Use Lattices_v3
      Use Operator_mod
      Use WaveFunction_mod
      Use MyMats
      use iso_fortran_env, only: output_unit, error_unit
      Implicit none

      Type Hopping_Matrix_type
         Integer                   :: N_bonds
         Complex (Kind=Kind(0.d0)), pointer :: T    (:)    !  This does not include  local terms.
         Complex (Kind=Kind(0.d0)), pointer :: T_loc(:)    !  This is just for the local matrix elements such as chemical potential.
         Integer                  , pointer :: list(:,:)
         ! T(N_b=1..N_bonds)
         ! List(N_b,1) = no_1
         ! List(N_b,2) = no_2
         ! List(N_b,3) = n_1
         ! List(N_b,4) = n_2
         ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)
         Integer                   :: N_Phi
         Real    (Kind=Kind(0.d0)) :: Phi_X, Phi_Y
         Logical                   :: Bulk
         ! N_Phi         = #  of flux quanta  piercieng the lattice
         ! Phi_X, Phi_Y  =  Twist
         ! Bulk          =  Twist as boundary condtion (Bulk=.F.)
         !               =  Twist as vector potential  (Bulk=.T.)

         ! For Checkerboard decomposition
         Integer                            :: N_Fam
         Integer                  , pointer :: L_Fam(:),  List_Fam(:,:,:)
         Real    (Kind=Kind(0.d0)), pointer :: Prop_Fam(:)

         Integer, private         , allocatable :: Multiplicity(:) !> Numer of times a given orbital occurs in the list of bonds, automatically computed
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
!> Following tests are done:
!>   -it checks that the allocated size of second index of List_Fam is
!>    at least maximum size of families. If the allocated size exceeds
!>    the required size, it issues a warning.
!>   -it checks that all bonds enter once and only once in the decomposition
!>   -it checks that every pair of bonds in each family commute
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
!> Default hopping for the square lattice.  Ham_T is the nearest
!> neighbour hopping and Ham_Chem the chemical potential.
!>
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

           do nf = 1,N_FL
              this(nf)%N_bonds = 0
              if ( abs(Ham_T_max) > Zero)  then
                 this(nf)%N_bonds = 2
                 Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                      &    this(nf)%T(this(nf)%N_bonds) )
                 nc = 0
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 1
                 this(nf)%List(nc,2) = 1
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 1

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
!> Default hopping for the triangular lattice.  Ham_T is the nearest
!> neighbour hopping and Ham_Chem the chemical potential.
!>
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
               nc = nc + 1
               this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
               this(nf)%List(nc,1) = 1
               this(nf)%List(nc,2) = 1
               this(nf)%List(nc,3) = 1
               this(nf)%List(nc,4) = 0

               nc = nc + 1
               this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
               this(nf)%List(nc,1) = 1
               this(nf)%List(nc,2) = 1
               this(nf)%List(nc,3) = 0
               this(nf)%List(nc,4) = 1

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
!> Default hopping for the kagome lattice.  Ham_T is the nearest
!> neighbour hopping and Ham_Chem the chemical potential.
!>
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
      do nf = 1,N_FL
         this(nf)%N_bonds = 0
         if ( abs(Ham_T_max) > Zero)  then
            this(nf)%N_bonds = 6
            Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                 &    this(nf)%T(this(nf)%N_bonds) )
            nc = 0
            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 1
            this(nf)%List(nc,2) = 2
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 0

            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 1
            this(nf)%List(nc,2) = 3
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 0

            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 2
            this(nf)%List(nc,2) = 3
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 0

            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 3
            this(nf)%List(nc,2) = 1
            this(nf)%List(nc,3) = 0
            this(nf)%List(nc,4) = 1

            nc = nc + 1
            this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
            this(nf)%List(nc,1) = 3
            this(nf)%List(nc,2) = 2
            this(nf)%List(nc,3) = -1
            this(nf)%List(nc,4) = 1

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
!> Default hopping for n-leg-ladder.  Ham_T is the nearest neighbour hopping along the chain,  Ham_T_perp  the
!> interrung hopping and H_chem the chemical potential.
!>
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
           
           
           Allocate( this(N_FL) )
           do nf = 1,N_FL
              this(nf)%N_bonds =  Latt_unit%Norb - 1
              Allocate (this(nf)%List( this(nf)%N_bonds,4 ), this(nf)%T( this(nf)%N_bonds ) )
              nc = 0
              do n = 1,Latt_unit%Norb  -1
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
           
           ! Set Checkerboard
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
              this(nf)%N_bonds = Latt_unit%Norb +  (Latt_unit%Norb - 1 )
              Allocate (this(nf)%List( this(nf)%N_bonds,4 ), &
                   &    this(nf)%T( this(nf)%N_bonds ) )
              nc = 0
              do n = 1,Latt_unit%Norb
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = n
                 this(nf)%List(nc,2) = n
                 this(nf)%List(nc,3) = 1
                 this(nf)%List(nc,4) = 0
              enddo
              
              do n = 1,Latt_unit%Norb -1
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
           ! Set Checkerboard
           If     ( Latt_Unit%Norb  == 1 ) then
              this(1)%N_FAM        = 2
           elseif ( Latt_Unit%Norb  == 2 ) then
              this(1)%N_FAM        = 3
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
                 Nf = 1
                 do no = 1,Latt_unit%Norb
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I ! Unit cell
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no ! The bond (See above)
                 enddo
              else
                 Nf = 2
                 do no = 1,Latt_unit%Norb
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no
                 enddo
              endif
           enddo
           do no = 1,Latt_unit%Norb - 1
              if (mod(no,2) == 1 ) then
                 Nf = 3
                 !Write(6,*)  NF, no + Latt_unit%Norb
                 do I = 1,Latt%N
                    this(1)%L_Fam(Nf) = this(1)%L_Fam(Nf) + 1
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),1) = I
                    this(1)%List_Fam(Nf,this(1)%L_Fam(Nf),2) = no + Latt_unit%Norb
                 enddo
              else
                 Nf = 4
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
!> Default hopping for Honeycomb lattice.  Ham_T is the nearest neighbour hopping along the chain, Lambda is  the
!> Kane-Mele term and H_chem the chemical potential.  **Note**  The Kane-Mele term is not yet implemented.
!>
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
           this(nf)%N_bonds =  3
           Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                &    this(nf)%T(this(nf)%N_bonds) )
           nc = 0
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  1
           this(nf)%List(nc,2) =  2
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  0

           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  2
           this(nf)%List(nc,2) =  1
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  1

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
!> Default hopping for a bilayer square. Ham_T1, Ham_T2, are the nearest neighbour hopping on the first and second layer and
!> Ham_T_perp   is the interlayer  hopping.
!>
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
              N_bonds = 0
              N_bonds = N_bonds + 1
              if (abs(Ham_Tperp_max) > Zero )  N_bonds = N_bonds + 1
              if (abs(Ham_T2_max)    > Zero )  N_bonds = N_bonds + 1
              this(nf)%N_bonds = N_bonds
              Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                   &    this(nf)%T(this(nf)%N_bonds) )
              nc = 0
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) = 1
              this(nf)%List(nc,2) = 1
              this(nf)%List(nc,3) = 1
              this(nf)%List(nc,4) = 0
              
              If (abs(Ham_Tperp_max) > Zero ) Then
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 1
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 0
              endif
              
              If (abs(Ham_T2_max) > Zero ) Then
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

           ! Set Checkerboard
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
              N_bonds = 0
              N_bonds = N_bonds + 2
              if (abs(Ham_Tperp_max) > Zero )  N_bonds = N_bonds + 1
              if (abs(Ham_T2_max)    > Zero )  N_bonds = N_bonds + 2
              this(nf)%N_bonds = N_bonds
              Allocate (this(nf)%List(this(nf)%N_bonds,4), &
                   &    this(nf)%T(this(nf)%N_bonds) )
              nc = 0
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) = 1
              this(nf)%List(nc,2) = 1
              this(nf)%List(nc,3) = 0
              this(nf)%List(nc,4) = 1
              
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) = 1
              this(nf)%List(nc,2) = 1
              this(nf)%List(nc,3) = 1
              this(nf)%List(nc,4) = 0
              
              If (abs(Ham_Tperp_max) > Zero ) Then
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 1
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 0
              endif
              
              If (abs(Ham_T2_max) > Zero ) Then
                 nc = nc + 1
                 this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
                 this(nf)%List(nc,1) = 2
                 this(nf)%List(nc,2) = 2
                 this(nf)%List(nc,3) = 0
                 this(nf)%List(nc,4) = 1
                 
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
           
           ! Set Checkerboard
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
!> Default hopping for a bilayer square. Ham_T1, Ham_T2, are the nearest neighbour hopping on the first and second layer and
!> Ham_T_perp   is the interlayer  hopping.
!>
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
           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  1
           this(nf)%List(nc,2) =  2
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  0

           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  2
           this(nf)%List(nc,2) =  1
           this(nf)%List(nc,3) =  0
           this(nf)%List(nc,4) =  1

           nc = nc + 1
           this(nf)%T(nc)    = cmplx(-Ham_T1_vec(nf),0.d0,kind(0.d0))
           this(nf)%List(nc,1) =  1
           this(nf)%List(nc,2) =  2
           this(nf)%List(nc,3) =  1
           this(nf)%List(nc,4) = -1

           If (abs(Ham_Tperp_Max) > Zero )  then
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  1
              this(nf)%List(nc,2) =  3
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  0

              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_Tperp_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  2
              this(nf)%List(nc,2) =  4
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  0
           endif
           If (abs(Ham_T2_Max) > Zero )  then
              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  1 + 2
              this(nf)%List(nc,2) =  2 + 2
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  0

              nc = nc + 1
              this(nf)%T(nc)    = cmplx(-Ham_T2_vec(nf),0.d0,kind(0.d0))
              this(nf)%List(nc,1) =  2 + 2
              this(nf)%List(nc,2) =  1 + 2
              this(nf)%List(nc,3) =  0
              this(nf)%List(nc,4) =  1

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

        ! Set Checkerboard
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
!> Given the checkerbord decompostion
!> the routine allocates and sets OP_T Set_Default_hopping_parameters_"Lattice" routine,
!> this routine generates the data for the symmetric decomposition.
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
!> Given the Hopping-matrix, and if required the checkerboard decomposion (i.e. private data of this module)
!> the routine allocates and sets OP_T.
!
!--------------------------------------------------------------------
      Subroutine Predefined_Hoppings_set_OPT(this,List,Invlist,Latt,  Latt_unit,  Dtau,Checkerboard, Symm,  OP_T )

        Implicit none

        type (Hopping_Matrix_type), allocatable             :: this(:)
        Integer, Intent(IN), Dimension(:,:)                 :: List, Invlist
        Type(Lattice),  Intent(in)                          :: Latt
        Type(Unit_cell),Intent(in)                          :: Latt_unit
        Real (Kind=Kind(0.d0)), Intent(In)                  :: Dtau
        Logical, Intent(IN)                                 :: Checkerboard, Symm

        Type(Operator), Intent(Out),  dimension(:,:), allocatable  :: Op_T



        ! Local
        Integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
        Integer                           :: n_1, n_2, Nb, n_f,l_f, n_l, N, nc, orb
        Real   (Kind=Kind(0.d0))          :: Ham_T, Ham_Chem,  Phi_X, Phi_Y
        Logical                           :: Bulk
        Complex(Kind=Kind(0.d0))          :: Z

        N_FL =  size(this,1)
        !Write(6,*)  'N_FL ', N_FL
        Ndim =  Latt%N * Latt_Unit%Norb

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
                 Op_T(n,nf)%O(1,1) =  this(nf)%T_Loc(list(n,2))
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
                       Op_T(1,nf)%O(I1,I1) = this(nf)%T_Loc(no_I)
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
                 ! Compute Multiplicity
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
                       Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                       I1   = Invlist(I,no_I)
                       J1   = Invlist(J,no_J)
                       nc = nc + 1
                       Op_T(nc,nf)%P(1) = I1
                       Op_T(nc,nf)%P(2) = J1
                       Op_T(nc,nf)%O(1,2) = this(nf)%T(Nb)*Z
                       Op_T(nc,nf)%O(2,1) = Conjg(this(nf)%T(Nb)*Z)
                       Op_T(nc,nf)%O(1,1) = this(nf)%T_loc(no_I)/this(1)%Multiplicity(no_I)
                       Op_T(nc,nf)%O(2,2) = this(nf)%T_loc(no_J)/this(1)%Multiplicity(no_J)
                       Op_T(nc,nf)%g = -Dtau*this(1)%Prop_Fam(n_f)
                       Op_T(nc,nf)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                       Call Op_set(Op_T(nc,nf))
                    Enddo
                 enddo
              enddo
           endif
        end select

      end Subroutine Predefined_Hoppings_set_OPT

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> The subroutine computes the kinetic energy based on the generic form of the
!> the hopping matrix.
!>
!--------------------------------------------------------------------
      Subroutine  Predefined_Hoppings_Compute_Kin(this,List,Invlist, Latt, Latt_unit, GRC, Z_Kin)

        Implicit none

        type (Hopping_Matrix_type), allocatable  :: this(:)
        Integer, Intent(IN), Dimension(:,:)                 :: List, Invlist
        Type(Lattice),  Intent(in)                          :: Latt
        Type(Unit_cell),Intent(in)                          :: Latt_unit
        Complex (Kind=Kind(0.d0)), intent(in), Dimension(:,:,:) :: GRC(:,:,:)
        Complex (Kind=Kind(0.d0)),  intent(out) :: Z_kin

        !Local
        Integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
        Integer                           :: n_1, n_2, Nb, n_f,l_f, n_l, N, nc
        Real   (Kind=Kind(0.d0))          :: Ham_T, Ham_Chem,  Phi_X, Phi_Y
        Logical                           :: Bulk
        Complex(Kind=Kind(0.d0))          :: Z


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
                    Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                    I1   = Invlist(I,no_I)
                    J1   = Invlist(J,no_J)
                    Z_Kin = Z_Kin + this(nf)%T(Nb)*Z * GRC(I1,J1,nf) + conjg(this(nf)%T(Nb)*Z)*GRC(J1,I1,nf)
                 enddo
                 Do no_I = 1, Latt_Unit%Norb
                    I1   = Invlist(I,no_I)
                    Z_Kin = Z_Kin   +  this(nf)%T_Loc(no_I)*GRC(I1,I1,nf)
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


        Z_hop = cmplx(1.d0,0.d0,kind(0.d0))

        xj_p =  real(latt%list(i,1) + del_1 ,kind(0.d0)) * latt%a1_p  +  real(latt%list(i,2) + del_2 ,kind(0.d0)) * latt%a2_p
        ! Check if you have crossed the boundary:  xj_p  = xjp_p + N1*L1_p  + N2*L2_p  with  xjp_p  in the set of lattice sites.
        N1 = 0; N2 = 0
        Call npbc(xjp_p, xj_p, Latt%L1_p, Latt%L2_p,  N1, N2)
        XB_p = real(N1,kind(0.d0))*Latt%L1_p  +  real(N2,kind(0.d0))*Latt%L2_p
        Do n = 1,2
           xj_p (n) = xj_p (n) + Latt_unit%Orb_pos_p(no_j,n)
           xjp_p(n) = xjp_p(n) + Latt_unit%Orb_pos_p(no_j,n)
        enddo
        xi_p    = real(latt%list(i,1), kind(0.d0)) * latt%a1_p  +  real(latt%list(i,2),kind(0.d0)) * latt%a2_p
        Do n = 1,2
           xi_p(n) = xi_p(n) +  Latt_unit%Orb_pos_p(no_i,n)
        Enddo
        !!Check that  xjp_p(:) + XB_p(:) =  xj_p(:)
        !!x1_p(:) = xjp_p(:) + XB_p
        !!Write(6,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,I2,3x,I2)")  x1_p(1),x1_p(2), xj_p(1), xj_p(2), N1,N2
        !! -->  i + del_1*a_1 + del_2* a_2  =  i' + N1*L1_p + N2*L2_p  with i' in the set of lattice points.


        ! The hopping increment.
        del_p  =  xj_p - xi_p

        !Twist
        pi = acos(-1.d0)
        A_p(:)  =   Flux_1 * Xnorm(Latt%a1_p) * latt%bZ1_p(:)  /  Xnorm(Latt%L1_p) + &
             &      Flux_2 * Xnorm(Latt%a2_p) * latt%bZ2_p(:)  /  Xnorm(Latt%L2_p)

        if (Bulk) then
           !Twist in bulk
           Z_hop = Z_hop * exp(cmplx(0.d0,Iscalar(A_p,del_p),Kind(0.d0)))
        else
           !Twist as boundary
           Z_hop = Z_hop * exp(cmplx(0.d0,Iscalar(A_p,XB_p ),Kind(0.d0)))
        endif

        !Orbital magnetic field (Landau gauge)
        Zero =  1.0E-8
        V  =  abs(Latt%L1_p(1) * Latt%L2_p(2)  -  Latt%L1_p(2) * Latt%L2_p(1) )
        If ( V > Zero )  then
           B = real(N_Phi,kind(0.d0))/V
           Z_hop = Z_hop*exp(cmplx(0.d0, -2.d0*pi* B * del_p(1) *  ( xj_p(2) + xi_p(2) )/2.d0,kind(0.d0) ) )
           ! Boundary
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
