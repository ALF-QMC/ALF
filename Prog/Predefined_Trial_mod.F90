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
!> This module provides a set of predefined trial wave functions.
!>
!
!--------------------------------------------------------------------


    Module Predefined_Trial

      Use runtime_error_mod
      Use Lattices_v3
      Use Operator_mod
      Use WaveFunction_mod
      Use MyMats
      Use Predefined_Hoppings
      use iso_fortran_env, only: output_unit, error_unit
      use Files_mod

      Implicit none


    contains

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
        Type (Hopping_Matrix_type), allocatable       :: Hopping_Matrix_tmp(:)
        Real (Kind=Kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer
        Logical                                       :: Checkerboard, Symm, Kekule_Trial

        Type (Lattice)                                :: Latt_Kekule
        Real (Kind=Kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2),x1_p(2), hop(3), del_p(2)
        Real (Kind=Kind(0.d0))  :: delta = 0.01, Ham_T1, Ham_T2, Ham_Tperp, dom, om

        Integer :: N, nf, I, I1, I2, nc, nc1, IK_u, I_u, J1, lp, J, N_Phi,  den_file, Nom = 200 , nw
        Logical :: Test=.false. ,  Bulk =.true.
        Complex (Kind=Kind(0.d0)) :: Z_norm, Z

        Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                Ham_T2_vec(:)
        Integer, allocatable ::   N_Phi_vec(:)



        Allocate(WF_L(N_FL),WF_R(N_FL))
        do n=1,N_FL
           Call WF_alloc(WF_L(n),Ndim,N_part)
           Call WF_alloc(WF_R(n),Ndim,N_part)
        enddo



        Allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
             &                                    N_Phi_vec(N_FL) )


        Checkerboard  = .false.
        Kekule_Trial  = .false.
        Symm          = .false.


        N_Phi       = 0
        Phi_X       = 0.d0
        Phi_Y       = 0.d0
        Bulk        = .false.
        Ham_T       = 1.d0
        Ham_T2      = 0.d0
        Ham_Tperp   = 0.d0
        Ham_Chem    = 0.d0
        Dtau        = 1.d0

        N_Phi_vec      = N_Phi
        Phi_X_vec      = Phi_X
        Phi_Y_vec      = Phi_Y
        Ham_T_vec      = Ham_T
        Ham_Tperp_vec  = Ham_Tperp
        Ham_T2_vec     = Ham_T2
        Ham_Chem_vec   = Ham_Chem


        Select case (str_to_upper(Lattice_type))

        case ("HONEYCOMB")
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
                          Write(error_unit,*) 'Error in  Predefined_TrialWaveFunction'
                          CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
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
        Case ("SQUARE")
           Phi_X_vec    = 0.01
           Call  Set_Default_hopping_parameters_square(Hopping_Matrix_tmp,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                &                                      Bulk,  N_Phi_vec, N_FL, &
                &                                      List, Invlist, Latt, Latt_unit )
        Case ("TRIANGULAR")
           Call  Set_Default_hopping_parameters_triangular(Hopping_Matrix_tmp,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                &                                      Bulk,  N_Phi_vec, N_FL, &
                &                                      List, Invlist, Latt, Latt_unit )
        Case ("KAGOME")
           Call  Set_Default_hopping_parameters_kagome(Hopping_Matrix_tmp,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                &                                      Bulk,  N_Phi_vec, N_FL, &
                &                                      List, Invlist, Latt, Latt_unit )
        Case ("N_LEG_LADDER")
           Ham_T_vec     = 1.d0
           Ham_Tperp_vec = 1.d0
           Phi_X_vec     = 0.01
           Call  Set_Default_hopping_parameters_n_leg_ladder(Hopping_Matrix_tmp, Ham_T_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
                &                                            Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
                &                                            List, Invlist, Latt, Latt_unit )
           !Case ("Honeycomb")
           !   Ham_Lambda = 0.d0
           !   Call  Set_Default_hopping_parameters_honeycomb(Hopping_Matrix_tmp, Ham_T, Ham_Lambda, Ham_Chem, Phi_X, Phi_Y, Bulk,  N_Phi, N_FL, &
           !        &                                       List, Invlist, Latt, Latt_unit )
        Case ("BILAYER_SQUARE")
           Ham_T_vec     = 1.d0
           Ham_T2_vec    = 0.d0
           Ham_Tperp_vec = 1.d0
           Phi_X_vec     = 0.00
           Call  Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix_tmp,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
                &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                &                                              List, Invlist, Latt, Latt_unit )
        Case ("BILAYER_HONEYCOMB")
           Ham_T_vec     = 1.d0
           Ham_T2_vec    = 0.d0
           Ham_Tperp_vec = 1.d0
           Phi_X_vec     = 0.00
           Call  Set_Default_hopping_parameters_Bilayer_honeycomb(Hopping_Matrix_tmp,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
                &                                                 Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                &                                                 List, Invlist, Latt, Latt_unit )
        case default
           Write(error_unit,*) 'No predefined trial wave function for this lattice.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        end Select

        If (str_to_upper(Lattice_type) .ne. "HONEYCOMB" )   &
             &     Call  Predefined_Hoppings_set_OPT(Hopping_Matrix_tmp,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_tmp )


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
           !Write(6,*) " Z_norm ", Z_norm
        enddo

        If (test) then
           DO  I = 1,NDim
              Write(6,*) Op_tmp(1,1)%E(I)
           enddo
           Open(newunit=den_file, file="Den_H0", status="unknown")
           delta = 2.d0 * Acos(-1.d0)/(Iscalar(Latt%L1_p,Latt%BZ1_p))
           dom =(Op_tmp(1,1)%E(ndim) - Op_tmp(1,1)%E(1)) / dble(Nom)
           om = Op_tmp(1,1)%E(1)
           do nw = 1,Nom + 1
               Z = cmplx(0.d0,0.d0,kind(0.d0))
               Do I  = 1,NDim 
                  Z  =  Z  + 1.d0/cmplx(om - Op_tmp(1,1)%E(I),delta, kind=kind(0.d0))
               enddo
               Write(den_file,"(F14.7,2x,F14.7)") om, - Aimag(Z)/(dble(Ndim)*acos(-1.d0))
               om = om + dom
           enddo  
           close(den_file)
        endif
        Do nf = 1,N_FL
           Call Op_clear(OP_tmp(1,nf),Ndim)
        enddo
        Deallocate (OP_tmp)
        Call Predefined_hoppings_clear(Hopping_Matrix_tmp)

        Deallocate (Ham_T_vec, Ham_Tperp_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec,  N_Phi_vec )

      end Subroutine Predefined_TrialWaveFunction



     end Module Predefined_Trial
