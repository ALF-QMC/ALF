!  Copyright (C) 2016 - 2018 The ALF project
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
!> This module provides a set of predefined lattices, hoppings, interactions,
!> trial wave functions as well as observables.
!>       
!
!--------------------------------------------------------------------

    Module Predefined_structures
      
      Use Lattices_v3
      Use Operator_mod
      Use WaveFunction_mod
      Use MyMats
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
        Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
        Integer :: I, nc, no
        
        select case (Lattice_type)
        case("Square")
           If ( L2==1 .and. L1 > 1 ) then
              Latt_Unit%N_coord   = 1
           elseif (L2 >1 .and. L1 > 1) then
              Latt_Unit%N_coord   = 2
           else
              Write(6,*) 'For one-dimnesional lattices set L2=1'
              Stop
           endif
           Latt_Unit%Norb      = 1
           Allocate (Latt_unit%Orb_pos_p(1,2))
           Latt_Unit%Orb_pos_p(1,:) = 0.d0 
           a1_p(1) =  1.0  ; a1_p(2) =  0.d0
           a2_p(1) =  0.0  ; a2_p(2) =  1.d0
           L1_p    =  dble(L1)*a1_p
           L2_p    =  dble(L2)*a2_p
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case("Honeycomb")
           If (L1==1 .or. L2==1 ) then
              Write(6,*) 'For one-dimensional lattices set : L2 = 1'
              stop
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
        case("Pi_Flux")
           If (L1==1 .or. L2==1 ) then
              Write(6,*) 'For one-dimensional lattices set : L2 = 1'
              stop
           endif
           Latt_Unit%Norb    = 2
           Latt_Unit%N_coord = 4
           a1_p(1) =  1.D0   ; a1_p(2) =   1.d0
           a2_p(1) =  1.D0   ; a2_p(2) =  -1.d0
           Allocate (Latt_Unit%Orb_pos_p(2,2))
           Latt_Unit%Orb_pos_p(1,:) = 0.d0 
           Latt_Unit%Orb_pos_p(2,:) = (a1_p(:) - a2_p(:))/2.d0 
           L1_p    =  dble(L1) * (a1_p - a2_p)/2.d0
           L2_p    =  dble(L2) * (a1_p + a2_p)/2.d0
           Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case default 
           Write(6,*) "Lattice not yet implemented!"
           Stop
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
      
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!>    Hopping, with ot without checkerboard  
!>    Per flavor, the  hopping is given by
!>    \f[  e^{ - \Delta \tau  H_t  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n) }   \f]
!>    If  _Symm_ is set to true and if  _Checkeborad_ is on, then  one will carry out a
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
!------------------------------------------------------------------
      Subroutine Predefined_Hopping(Lattice_type, Ndim, List,Invlist,Latt,  Latt_unit,  &
           &                        Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, &
           &                        N_FL,  Checkerboard, Symm, OP_T, Dimer )

        Implicit none

        Character (len=64), Intent(IN)               :: Lattice_type
        Integer, Intent(IN)                          :: Ndim, N_FL
        Integer, Intent(IN), Dimension(:,:)          :: List, Invlist
        Type(Lattice),  Intent(in)                   :: Latt
        Type(Unit_cell),Intent(in)                   :: Latt_unit
        Real (Kind=Kind(0.d0)), Intent(In)           :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y
        Logical                                      :: Checkerboard, Symm
        Real(Kind=Kind(0.d0)), Intent(IN), Optional  :: Dimer
        
        Type(Operator), Intent(Out),  dimension(:,:), allocatable  :: Op_T 


        !Local
        Integer :: I, I1, J1, I2, n, Ncheck,nc, nc1, no, N_Fam, L_FAM, Ix, Iy, n1
        Complex (Kind=Kind(0.d0)) :: ZX, ZY
        Real    (Kind=Kind(0.d0)) :: del_p(2), X, g
        
        If ( .not. Checkerboard) then
           allocate(Op_T(1,N_FL))
           do n = 1,N_FL
              Call Op_make(Op_T(1,n),Ndim)
              nc = 1
              Select case (Lattice_type)
              Case ("Square")
                 If ( Latt_unit%N_coord == 2 ) then 
                    ZX  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/Xnorm(Latt%L1_p), kind=kind(0.d0) ) )
                    ZY  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_Y/Xnorm(Latt%L2_p), kind=kind(0.d0) ) )
                    DO I = 1, Latt%N
                       I1 = Latt%nnlist(I,1,0)
                       I2 = Latt%nnlist(I,0,1)
                       If ( Latt%list(I,1) == 0 ) then
                          Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*ZX
                          Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*conjg(ZX)
                       else
                          Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*ZX
                          Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(ZX)
                       endif
                       If ( Latt%list(I,2) == 0 ) then
                          Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T*XB_Y,    0.d0, kind(0.D0))*ZY
                          Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T*XB_Y,    0.d0, kind(0.D0))*conjg(ZY)
                       else
                          Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T     ,    0.d0, kind(0.D0))*ZY
                          Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T     ,    0.d0, kind(0.D0))*conjg(ZY)
                       endif
                       Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    Enddo
                    if (Present(Dimer)) then
                       Do I = 1, Latt%N
                          Ix = Latt%list(I,1); Iy = Latt%list(I,2)
                          If ( mod(Ix + Iy,2)  == 0 ) then
                             I1 = I
                             I2 = Latt%nnlist(I,1,0)
                             Op_T(nc,n)%O(I1,I2)  = Op_T(nc,n)%O(I1,I2) + Dimer
                             Op_T(nc,n)%O(I2,I1)  = Op_T(nc,n)%O(I2,I1) + Dimer
                          endif
                       enddo
                    endif
                 else
                    ZX  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/Xnorm(Latt%L1_p), kind=kind(0.d0) ) )
                    DO I = 1, Latt%N
                       I1 = Latt%nnlist(I,1,0)
                       If ( Latt%list(I,1) == 0 ) then
                          Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*ZX
                          Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*conjg(ZX)
                       else
                          Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*ZX
                          Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(ZX)
                       endif
                       Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    Enddo
                    if (Present(Dimer)) then
                       Do I = 1, Latt%N
                          Ix = Latt%list(I,1); Iy = Latt%list(I,2)
                          If ( mod(Ix + Iy,2)  == 0 ) then
                             I1 = I
                             I2 = Latt%nnlist(I,1,0)
                             Op_T(nc,n)%O(I1,I2)  = Op_T(nc,n)%O(I1,I2) + Dimer
                             Op_T(nc,n)%O(I2,I1)  = Op_T(nc,n)%O(I2,I1) + Dimer
                          endif
                       enddo
                    endif
                 endif
              Case ("Honeycomb")
                 X = 2.d0 * acos(-1.d0)*Phi_X / ( Xnorm(Latt%L1_p) * (Xnorm(Latt%a1_p)**2)  )
                 DO I = 1, Latt%N
                    do no = 1,Latt_unit%Norb
                       I1 = Invlist(I,no)
                       Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    enddo
                    I1 = Invlist(I,1)
                    J1 = I1
                    Do nc1 = 1,Latt_unit%N_coord
                       select case (nc1)
                       case (1)
                          J1 = invlist(I,2)
                          del_p(:)  =  Latt_unit%Orb_pos_p(2,:) 
                       case (2)
                          J1 = invlist(Latt%nnlist(I,1,-1),2)
                          del_p(:)   =  Latt%a1_p(:) - Latt%a2_p(:)  + Latt_unit%Orb_pos_p(2,:)
                       case (3)
                          J1 = invlist(Latt%nnlist(I,0,-1),2) 
                          del_p(:)   =  - Latt%a2_p(:) +  Latt_unit%Orb_pos_p(2,:) 
                       case default
                          Write(6,*) ' Error in  Ham_Hop '  
                          Stop
                       end select
                       ZX = exp( cmplx(0.d0, X*Iscalar(Latt%a1_p,del_p), kind(0.D0) ) )
                       Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) * ZX
                       Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) * CONJG(ZX) 
                    Enddo
                 Enddo
              case("Pi_Flux")
                 DO I = 1, Latt%N
                    do no = 1,Latt_unit%Norb
                       I1 = Invlist(I,no)
                       Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    enddo
                    I1 = Invlist(I,1)
                    J1 = I1
                    Do nc1 = 1,Latt_unit%N_coord
                       select case (nc1)
                       case (1)
                          J1 = invlist(I,2) 
                       case (2)
                          J1 = invlist(Latt%nnlist(I,0, 1),2) 
                       case (3)
                          J1 = invlist(Latt%nnlist(I,-1,1),2) 
                       case (4)
                          J1 = invlist(Latt%nnlist(I,-1,0),2) 
                       case default
                          Write(6,*) ' Error in  Ham_Hop '  
                          Stop
                       end select
                       if (nc1 == 1 ) then
                          Op_T(nc,n)%O(I1,J1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                          Op_T(nc,n)%O(J1,I1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                       Else
                          Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                          Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                       endif
                    Enddo
                 Enddo
              case default 
                 Write(6,*) "Lattice not yet implemented!"
                 Stop
              end select
              Do I = 1,Ndim
                 Op_T(nc,n)%P(i) = i 
              Enddo
              if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                 Op_T(nc,n)%g = 0.d0
              else
                 Op_T(nc,n)%g = -Dtau
              endif
              Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
              Call Op_set(Op_T(nc,n))
              !Do I = 1,Size(Op_T(nc,n)%E,1)
              !   Write(6,*) Op_T(nc,n)%E(I)
              !Enddo
           Enddo
        else
           select case (Lattice_type)
           case ("Square")
              If (Symm) then
                 N_Fam = 2*(2*Latt_Unit%N_coord) -1  !  Number of Families = 4
                 L_Fam = Latt%N/2          !  Length of Families = LQ/4
                 Allocate(Op_T(N_FAM*L_FAM,N_FL))
                 do n = 1,N_FL
                    do I  =  1, N_FAM*L_FAM
                       call Op_make(Op_T(I,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,N_Fam
                       do I = 1,Latt%N
                          if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                             I1 = I
                             nc = nc + 1
                             select case (nc1)
                             case(1)
                                I2 = latt%nnlist(I1, 1, 0) 
                                g = -Dtau/2.d0
                             case(2)
                                I2 = latt%nnlist(I1, 0, 1)
                                g = -Dtau/2.d0
                             case(3)
                                I2 = latt%nnlist(I1,-1, 0)
                                g = -Dtau/2.d0
                             case(4)
                                I2 = latt%nnlist(I1, 0,-1)
                                g = -Dtau
                             case(5)
                                I2 = latt%nnlist(I1,-1, 0)
                                g = -Dtau/2.d0
                             case(6)
                                I2 = latt%nnlist(I1, 0, 1)
                                g = -Dtau/2.d0
                             case(7)
                                I2 = latt%nnlist(I1, 1, 0) 
                                g = -Dtau/2.d0
                             end select
                             !Write(6,*) nc,nc1, Latt%List(I1,1), Latt%List(I1,2),Latt%List(I2,1), Latt%List(I2,2), I1, I2
                             Op_T(nc,n)%P(1) = I1
                             Op_T(nc,n)%P(2) = I2
                             Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) 
                             Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0))
                             Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                             Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                             if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                                Op_T(nc,n)%g = 0.d0
                             else
                                Op_T(nc,n)%g = g
                             endif
                             Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                             Call Op_set( Op_T(nc,n) )
                          endif
                       Enddo
                    Enddo
                 Enddo
              Else
                 N_Fam = 2*Latt_unit%N_coord !  Number of Families = 4
                 L_Fam = Latt%N/2  !  Length of Families = LQ/4
                 Allocate(Op_T(N_FAM*L_FAM,N_FL))
                 do n = 1,N_FL
                    do I  =  1, N_FAM*L_FAM
                       call Op_make(Op_T(I,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,N_Fam
                       do I = 1,Latt%N
                          if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                             I1 = I
                             nc = nc + 1
                             if (nc1 == 1 ) I2 = latt%nnlist(I1, 1, 0) 
                             if (nc1 == 2 ) I2 = latt%nnlist(I1, 0, 1)
                             if (nc1 == 3 ) I2 = latt%nnlist(I1,-1, 0)
                             if (nc1 == 4 ) I2 = latt%nnlist(I1, 0,-1)
                             !Write(6,*) nc,nc1, Latt%List(I1,1), Latt%List(I1,2),Latt%List(I2,1), Latt%List(I2,2), I1, I2
                             Op_T(nc,n)%P(1) = I1
                             Op_T(nc,n)%P(2) = I2
                             Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) 
                             Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0))
                             Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                             Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                             if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                                Op_T(nc,n)%g = 0.d0
                             else
                                Op_T(nc,n)%g = -Dtau
                             endif
                             Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                             Call Op_set( Op_T(nc,n) )
                          endif
                       Enddo
                    Enddo
                 Enddo
              Endif
           case ("Pi_Flux")
              If (Symm) Then
                 Write(6,*) 'Pi-Flux Symm is not yet implemented'
                 Stop
              else
                 Allocate(Op_T(Latt_unit%N_coord*Latt%N,N_FL))
                 do n = 1,N_FL
                    !Write(6,*) 'N_coord, Latt%N ',  N_coord, Latt%N
                    do i  =  1, Latt_unit%N_coord*Latt%N
                       call Op_make(Op_T(i,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,Latt_unit%N_coord
                       do I = 1,Latt%N
                          I1 = Invlist(I,1)
                          nc = nc + 1
                          If (nc1 == 1 )  I2 = invlist(I,2)
                          If (nc1 == 2 )  I2 = invlist(Latt%nnlist(I,0, 1),2) 
                          If (nc1 == 3 )  I2 = invlist(Latt%nnlist(I,-1,1),2)
                          If (nc1 == 4 )  I2 = invlist(Latt%nnlist(I,-1,0),2) 
                          Op_T(nc,n)%P(1) = I1
                          Op_T(nc,n)%P(2) = I2
                          if (nc1 == 1 ) then
                             Op_T(nc,n)%O(1,2) = cmplx( Ham_T,    0.d0, kind(0.D0))
                             Op_T(nc,n)%O(2,1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                          Else
                             Op_T(nc,n)%O(1,2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                             Op_T(nc,n)%O(2,1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                          endif
                          Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                          if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                             Op_T(nc,n)%g = 0.d0
                          else
                             Op_T(nc,n)%g = -Dtau
                          endif
                          Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                          Call Op_set( Op_T(nc,n) )
                       Enddo
                    Enddo
                 Enddo
              Endif
           case ("Honeycomb")
              If (Symm) then
                 Allocate(Op_T((2*Latt_unit%N_coord-1)*Latt%N,N_FL))
                 do n = 1,N_FL
                    do i  =  1, (2*Latt_unit%N_coord-1)*Latt%N
                       call Op_make(Op_T(i,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,(2*Latt_unit%N_coord -1) 
                       do I = 1,Latt%N
                          I1 = invlist(I,1)
                          nc = nc + 1
                          select case (nc1)
                          case(1)
                             I2 = invlist(I,2)
                             g = -Dtau/2.d0
                          case(2)
                             I2 = invlist(latt%nnlist(I,1,-1),2)
                             g = -Dtau/2.d0
                          case(3)
                             I2 = invlist(latt%nnlist(I,0,-1),2)
                             g = -Dtau
                          case(4)
                             I2 = invlist(latt%nnlist(I,1,-1),2)
                             g = -Dtau/2.d0
                          case(5)
                             I2 = invlist(I,2)
                             g = -Dtau/2.d0
                          end select
                          Op_T(nc,n)%P(1) = I1
                          Op_T(nc,n)%P(2) = I2
                          Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0))
                          Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0))
                          if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                             Op_T(nc,n)%g = 0.d0
                          else
                             Op_T(nc,n)%g = g
                          endif
                          Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                          Call Op_set( Op_T(nc,n) )
                       Enddo
                    Enddo
                 Enddo
              Else
                 Allocate(Op_T(Latt_unit%N_coord*Latt%N,N_FL))
                 do n = 1,N_FL
                    do i  =  1, Latt_unit%N_coord*Latt%N
                       call Op_make(Op_T(i,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,Latt_unit%N_coord
                       do I = 1,Latt%N
                          I1 = invlist(I,1)
                          nc = nc + 1
                          if (nc1 == 1 ) I2 = invlist(I,2)
                          if (nc1 == 2 ) I2 = invlist(latt%nnlist(I,1,-1),2)
                          if (nc1 == 3 ) I2 = invlist(latt%nnlist(I,0,-1),2)
                          Op_T(nc,n)%P(1) = I1
                          Op_T(nc,n)%P(2) = I2
                          Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0))
                          Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0))
                          if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                             Op_T(nc,n)%g = 0.d0
                          else
                             Op_T(nc,n)%g = -Dtau
                          endif
                          Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                          Call Op_set( Op_T(nc,n) )
                       Enddo
                    Enddo
                 Enddo
              Endif
              case default
                 Write(6,*) "Checkeboard is not implemented for this lattice"
              Stop
           End select
        endif
        
      End Subroutine Predefined_Hopping

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!>    Sets the trial wave function corresponding to the solution of the non-interacting
!>    tight binding Hamiltonian on the given lattice. Twisted boundary conditions (Phi_X=0.01)
!>    are implemented so as to generate a non-degenerate trial wave functions. 
!> @param [in]  Lattice_type
!>    Character(64)
!> \verbatim
!>    Square,  Honeycomb, Pi_Flux
!> \endverbatim
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
!>    List(I=1.. Ndim,1)    =   Unit cell of site I    
!>    List(I=1.. Ndim,2)    =   Orbital index  of site I    
!>    Invlist(Unit_cell,Orbital) = site I    
!> \endverbatim
!> @param [in]    Latt
!>    Type(Lattice)
!> \verbatim
!>    The Lattice
!> \endverbatim
!> @param [in]  N_part
!>    Integer
!> \verbatim
!>    Particle number for each flavor
!> \endverbatim
!> @param [in]  N_FL
!>    Integer
!> \verbatim
!>    Flavor
!> \endverbatim
!> @param [out]  WF_L, WF_R
!>    Type(Wavefunction)(N_FL)
!> \verbatim
!>    Wavefunction
!>    Also sets the degeneracy:  E(N_part + 1) - E(N_part). Energy eigenvalues are ordered in ascending order.
!> \endverbatim
!>       
!------------------------------------------------------------------       
      Subroutine Predefined_TrialWaveFunction(Lattice_type, Ndim,  List,Invlist,Latt, Latt_unit, &
           &                                  N_part, N_FL,  WF_L, WF_R) 



        Implicit none
        Character (len=64), Intent(IN)                :: Lattice_type
        Integer, Intent(IN)                           :: Ndim, N_FL, N_part
        Integer, Intent(IN), Dimension(:,:)           :: List, Invlist
        Type(Lattice),   Intent(in)                   :: Latt
        Type(Unit_cell), Intent(in)                   :: Latt_Unit
        Type(WaveFunction), Intent(out), Dimension(:), allocatable :: WF_L, WF_R

        
        Type(Operator),  dimension(:,:), allocatable  :: OP_tmp
        Real (Kind=Kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer
        Logical                                       :: Checkerboard, Symm, Kekule_Trial

        Type (Lattice)                                :: Latt_Kekule
        Real (Kind=Kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2),x1_p(2), hop(3), del_p(2)
        Real (Kind=Kind(0.d0))  :: delta = 0.01, Ham_T1
        
        Integer :: N, nf, I, I1, I2, nc, nc1, IK_u, I_u, J1, lp, J
        Logical :: Test=.false.
        Complex (Kind=Kind(0.d0)) :: Z_norm
        
        
        Allocate(WF_L(N_FL),WF_R(N_FL))
        do n=1,N_FL
           Call WF_alloc(WF_L(n),Ndim,N_part)
           Call WF_alloc(WF_R(n),Ndim,N_part)
        enddo



        Kekule_Trial = .false.
        Select case (Lattice_type)
           
        case ("Honeycomb")
           If (Kekule_Trial) then
              !  Kekule Mass term to avoid  degeneracy at half-filling.
              Allocate(Op_Tmp(1,N_FL))
              do n = 1,N_FL
                 Call Op_make(Op_Tmp(1,n),Ndim)
              enddo
              
              If (test) then
                 Open (Unit=31,status="Unknown", file="Tmp1_latt") 
                 Open (Unit=32,status="Unknown", file="Tmp2_latt") 
                 Open (Unit=33,status="Unknown", file="Tmp3_latt")
              endif
              A1_p = 2.d0 * Latt%a1_p  - Latt%a2_p
              A2_p =        Latt%a1_p  + Latt%a2_p
              L1_p = Latt%L1_p
              L2_p = Latt%L2_p
              Call Make_Lattice( L1_p, L2_p, A1_p,  A2_p, Latt_Kekule)
              Call Print_latt(Latt_Kekule)
              
              DO I = 1, Latt_Kekule%N
                 x_p = dble(Latt_Kekule%list(I,1))*Latt_Kekule%a1_p + dble(Latt_Kekule%list(I,2))*Latt_Kekule%a2_p
                 IK_u   = Inv_R(x_p,Latt)
                 do nc  = 1, 3
                    select case (nc)
                    case (1)
                       I_u    =  IK_u 
                       hop(1) =  1.d0 + delta
                       hop(2) =  1.d0 - delta
                       hop(3) =  1.d0
                    case (2)
                       I_u    = Latt%nnlist(IK_u,0,1) 
                       hop(1) =  1.d0 
                       hop(2) =  1.d0 + delta
                       hop(3) =  1.d0 - delta
                    case (3)
                       I_u     = Latt%nnlist(IK_u,1,0) 
                       hop(1) =  1.d0 - delta
                       hop(2) =  1.d0 
                       hop(3) =  1.d0 + delta
                    end select
                    x_p = dble(Latt%list(I_u,1))*Latt%a1_p + dble(Latt%list(I_u,2))*Latt%a2_p
                    I1 = invlist(I_u,1)
                    do nc1 = 1,3
                       select case (nc1)
                       case (1)
                          J1 = invlist(I_u,2)
                          del_p(:)  =  Latt_unit%Orb_pos_p(2,:) 
                       case (2)
                          J1 = invlist(Latt%nnlist(I_u,1,-1),2)
                          del_p(:)   =  Latt%a1_p(:) - Latt%a2_p(:)  + Latt_unit%Orb_pos_p(2,:)
                       case (3)
                          J1 = invlist(Latt%nnlist(I_u,0,-1),2) 
                          del_p(:)   =  - Latt%a2_p(:) +  Latt_unit%Orb_pos_p(2,:) 
                       end select
                       
                       x1_p = X_p + del_p
                       lp = 32
                       if (hop(nc1) > 1.d0 ) lp = 33
                       if (hop(nc1) < 1.d0 ) lp = 31
                       If (test) then
                          Write(lp,"(F14.7,2x,F14.7)")  x_p(1), x_p(2)
                          Write(lp,"(F14.7,2x,F14.7)")  x1_p(1), x1_p(2)
                          Write(lp,*)
                       endif
                       do n = 1,N_FL
                          Op_Tmp(1,n)%O(I1,J1) =   cmplx( - hop(nc1),    0.d0, kind(0.D0))
                          Op_Tmp(1,n)%O(J1,I1) =   cmplx( - hop(nc1),    0.d0, kind(0.D0))
                       enddo
                    enddo
                 enddo
              Enddo
              do n = 1,N_FL
                 Do I = 1,Ndim
                    Op_Tmp(1,n)%P(i) = i 
                 Enddo
                 Op_Tmp(1,n)%g    = cmplx(1.d0, 0.d0,kind(0.d0))
                 Op_Tmp(1,n)%alpha= cmplx(0.d0,0.d0, kind(0.D0))
                 Call Op_set(Op_Tmp(1,n))
              Enddo
              If (test) then
                 Close(31)
                 Close(32)
                 Close(33)
              endif
           else
              Ham_T = 1.d0
              Ham_T1 = -delta*Ham_T
              Allocate(Op_Tmp(1,N_FL))
              do n = 1,N_FL
                 Call Op_make(Op_Tmp(1,n),Ndim)
                 Do I = 1,Latt%N
                    I1 = Invlist(I,1)
                    Do nc1 = 1,Latt_unit%N_coord
                       select case (nc1)
                       case (1)
                          J1 = invlist(I,2)
                       case (2)
                          J1 = invlist(Latt%nnlist(I,1,-1),2)
                       case (3)
                          J1 = invlist(Latt%nnlist(I,0,-1),2) 
                       case default
                          Write(6,*) ' Error in  Ham_Hop '  
                          Stop
                       end select
                       Op_Tmp(1,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) 
                       Op_Tmp(1,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) 
                    Enddo
                    I1 = invlist(Latt%nnlist(I,1,-1),2)
                    J1 = invlist(Latt%nnlist(I,0, 1),1)
                    Op_Tmp(1,n)%O(I1,J1) = cmplx(-Ham_T1,    0.d0, kind(0.D0)) 
                    Op_Tmp(1,n)%O(J1,I1) = cmplx(-Ham_T1,    0.d0, kind(0.D0))
                 enddo
                 do I = 1,Ndim
                    Op_Tmp(1,n)%P(i) = i 
                 Enddo
                 Op_Tmp(1,n)%g    = cmplx(1.d0, 0.d0,kind(0.d0))
                 Op_Tmp(1,n)%alpha= cmplx(0.d0,0.d0, kind(0.D0))
                 Call Op_set(Op_Tmp(1,n))
              Enddo
           endif
        case default 
           Dtau     = 1.d0
           Ham_T    = 1.d0
           Ham_Chem = 0.d0
           XB_X     = 1.d0
           XB_Y     = 1.d0
           Phi_X    = 0.00
           Phi_Y    = 0.d0
           Dimer    = 0.d0
           Checkerboard  = .false.
           Symm          = .false.
           !If (Lattice_type == "Square"  ) then 
           !   Dimer    = 0.001d0
           !else
           Phi_X    = 0.01
           !endif
           
           Call Predefined_Hopping(Lattice_type ,Ndim, List,Invlist,Latt, Latt_Unit, &
                &                    Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, &
                &                    N_FL,  Checkerboard, Symm,  OP_tmp, Dimer )
           
           
           
        end Select
        
        
        Do nf = 1,N_FL
           Call Diag(Op_tmp(1,nf)%O,Op_tmp(1,nf)%U,Op_tmp(1,nf)%E)
           do I2=1,N_part
              do I1=1,Ndim
                 WF_L(nf)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
                 WF_R(nf)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
              enddo
           enddo
           WF_L(nf)%Degen = Op_tmp(1,nf)%E(N_part+1) - Op_tmp(1,nf)%E(N_part)
           WF_R(nf)%Degen = Op_tmp(1,nf)%E(N_part+1) - Op_tmp(1,nf)%E(N_part)
        enddo

        Do nf = 1,N_FL
           Call WF_overlap(WF_L(nf), WF_R(nf), Z_norm)
           Write(6,*) " Z_norm ", Z_norm
        enddo
        
        If (test) then
           DO  I = 1,NDim
              Write(6,*) Op_tmp(1,1)%E(I)
           enddo
           Do I = 1,Ndim
              do J = 1,Ndim
                 Write(6,*) Op_tmp(1,1)%O(I,J)
              enddo
           enddo
        endif
        Do nf = 1,N_FL
           Call Op_clear(OP_tmp(1,nf),Ndim)
        enddo
        Deallocate (OP_tmp)

       end Subroutine Predefined_TrialWaveFunction

      
      
    end Module Predefined_structures
