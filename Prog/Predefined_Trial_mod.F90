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

module Predefined_Trial

   use runtime_error_mod
   use Lattices_v3
   use Operator_mod
   use WaveFunction_mod
   use MyMats
   use Predefined_Hoppings
   use random_wrap
   use iso_fortran_env, only: output_unit, error_unit

   implicit none

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
   subroutine Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
        &                                  N_part, alpha, N_FL, WF_L, WF_R)

      implicit none
      character(len=64), intent(IN)                :: Lattice_type
      integer, intent(IN)                          :: Ndim, N_FL, N_part
      integer, intent(IN), dimension(:, :)         :: List, Invlist
      real(Kind=kind(0.d0)), intent(in)            :: alpha
      type(Lattice), intent(in)                    :: Latt
      type(Unit_cell), intent(in)                  :: Latt_Unit
      type(WaveFunction), intent(out), dimension(:), allocatable :: WF_L, WF_R

      type(operator), dimension(:, :), allocatable :: OP_tmp
      type(Hopping_Matrix_type), allocatable       :: Hopping_Matrix_tmp(:)
      real(Kind=kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer
      logical                                      :: Checkerboard, Symm, Kekule_Trial

      type(Lattice)                                :: Latt_Kekule
      real(Kind=kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2), x1_p(2), hop(3), del_p(2)
      real(Kind=kind(0.d0))  :: delta = 0.01, ham_tx, ham_ty, Ham_T1, Ham_T2, Ham_Tperp, stag_mass, stag_sgn

      integer :: N, nf, I, I1, I2, nc, nc1, IK_u, I_u, J1, lp, J, k, k1, N_Phi, no, j2, k2
      logical :: Test = .false., Bulk = .true.
      complex(Kind=kind(0.d0)) :: Z_norm

      real(Kind=kind(0.d0)), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
             &  ham_tx_vec(:), ham_ty_vec(:), Ham_T2_vec(:)
      integer, allocatable ::   N_Phi_vec(:)

      allocate (WF_L(N_FL), WF_R(N_FL))
      do n = 1, N_FL
         call WF_alloc(WF_L(n), Ndim, N_part)
         call WF_alloc(WF_R(n), Ndim, N_part)
      end do

      allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
           &    ham_tx_vec(N_FL), ham_ty_vec(N_FL), N_Phi_vec(N_FL))

      Checkerboard = .false.
      Kekule_Trial = .false.
      Symm = .false.

      N_Phi = 0
      Phi_X = 0.d0
      Phi_Y = 0.d0
      Bulk = .false.
      Ham_T = 1.d0
      Ham_T2 = 0.d0
      Ham_Tperp = 0.d0
      Ham_Chem = 0.d0
      Dtau = 1.d0

      N_Phi_vec = N_Phi
      Phi_X_vec = Phi_X
      Phi_Y_vec = Phi_Y
      Ham_T_vec = Ham_T
      Ham_Tperp_vec = Ham_Tperp
      Ham_T2_vec = Ham_T2
      Ham_Chem_vec = Ham_Chem

      select case (Lattice_type)

      case ('bilayer_square')
         allocate (op_tmp(1, n_fl))
         do n = 1, n_fl
            call op_make(op_tmp(1, n), ndim)
            do i = 1, ndim
               op_tmp(1, n)%P(i) = i
            end do
            op_tmp(1, n)%g = cmplx(1.d0, 0.d0, kind(0.d0))
            op_tmp(1, n)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
         end do

         do i = 1, latt%n
            delta = 0.005d0*ranf_wrap()
            ham_t = 1.d0 + delta

            i1 = invlist(i, 1)
            j1 = invlist(latt%nnlist(i, 1, 0), 1)
            k1 = invlist(latt%nnlist(i, 0, 1), 1)

            i2 = invlist(i, 2)
            j2 = invlist(latt%nnlist(i, 1, 0), 2)
            k2 = invlist(latt%nnlist(i, 0, 1), 2)

            ham_tx_vec(1) = ham_t; 
            ham_ty_vec(1) = ham_t*alpha; 
            do n = 1, n_fl
               op_tmp(1, n)%o(i1, j1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(j1, i1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(i1, k1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(k1, i1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))

               op_tmp(1, n)%o(i2, j2) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(j2, i2) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(i2, k2) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(k2, i2) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
            end do
         end do

         do n = 1, n_fl
            call op_set(op_tmp(1, n))
         end do

      case ('square_anisotropic')
         allocate (op_tmp(1, n_fl))
         do n = 1, n_fl
            call op_make(op_tmp(1, n), ndim)
            do i = 1, ndim
               op_tmp(1, n)%P(i) = i
            end do
            op_tmp(1, n)%g = cmplx(1.d0, 0.d0, kind(0.d0))
            op_tmp(1, n)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
         end do

         do i = 1, latt%n
            delta = 0.005d0*ranf_wrap()
            ham_t = 1.d0!+delta

            i1 = invlist(i, 1)
            j1 = invlist(latt%nnlist(i, 1, 0), 1)
            k1 = invlist(latt%nnlist(i, 0, 1), 1)

            ham_tx_vec(1) = ham_t; 
            ham_tx_vec(2) = ham_t*alpha; 
            ham_ty_vec(1) = ham_t*alpha; 
            ham_ty_vec(2) = ham_t; 
            do n = 1, n_fl
               op_tmp(1, n)%o(i1, j1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(j1, i1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(i1, k1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(k1, i1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
            end do
         end do

         do n = 1, n_fl
            call op_set(op_tmp(1, n))
         end do

      case ('N_leg_ladder')
         allocate (op_tmp(1, n_fl))
         do n = 1, n_fl
            call op_make(op_tmp(1, n), ndim)
            do i = 1, ndim
               op_tmp(1, n)%P(i) = i
            end do
            op_tmp(1, n)%g = cmplx(1.d0, 0.d0, kind(0.d0))
            op_tmp(1, n)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
         end do

         do i = 1, latt%n
            do no = 1, Latt_unit%Norb - 1
               delta = 0.005d0*ranf_wrap()
               ham_t = 1.d0 + delta

               i1 = invlist(i, no)
               j1 = invlist(latt%nnlist(i, 1, 0), no)
               k1 = invlist(i, no + 1)

               ham_tx_vec(1) = ham_t; 
               ham_tx_vec(2) = ham_t*alpha; 
               ham_ty_vec(1) = ham_t*alpha; 
               ham_ty_vec(2) = ham_t; 
               do n = 1, n_fl
                  op_tmp(1, n)%o(i1, j1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
                  op_tmp(1, n)%o(j1, i1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
                  op_tmp(1, n)%o(i1, k1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
                  op_tmp(1, n)%o(k1, i1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
               end do
            end do
            delta = 0.005d0*ranf_wrap()
            ham_t = 1.d0 + delta
            no = Latt_unit%Norb
            i1 = invlist(i, no)
            j1 = invlist(latt%nnlist(i, 1, 0), no)
            ham_tx_vec(1) = ham_t; 
            ham_tx_vec(2) = ham_t*alpha; 
            do n = 1, n_fl
               op_tmp(1, n)%o(i1, j1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
               op_tmp(1, n)%o(j1, i1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
            end do
         end do

         do n = 1, n_fl
            call op_set(op_tmp(1, n))
         end do

      case ("qbt")
         Ham_T_vec = 1.d0
         Ham_T2_vec = 0.5d0
         !Phi_X_vec     = 0.01
         call set_hopping_parameters_pi_flux_qbt(Hopping_Matrix_tmp, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, &
             & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)

         call Predefined_Hoppings_set_OPT(Hopping_Matrix_tmp, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_tmp)

           !! add stagger mass to avoid the degeneracy of qbt
         stag_mass = 0.01
         do nf = 1, N_FL
            do I = 1, Latt%N
            do no = 1, Latt_unit%norb
               stag_sgn = 1.d0
               if (mod(no, 2) .eq. 0) stag_sgn = -1.d0
               I1 = invlist(I, no)
                 !! onsite sublattice mass
               op_tmp(1, nf)%o(I1, I1) = stag_sgn*stag_mass
            end do
            end do
         end do

      case default
         write (error_unit, *) 'No predefined trial wave function for this lattice.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end select

      do nf = 1, N_FL
         call Diag(Op_tmp(1, nf)%O, Op_tmp(1, nf)%U, Op_tmp(1, nf)%E)
         do I2 = 1, N_part
            do I1 = 1, Ndim
               WF_L(nf)%P(I1, I2) = Op_tmp(1, nf)%U(I1, I2)
               WF_R(nf)%P(I1, I2) = Op_tmp(1, nf)%U(I1, I2)
            end do
         end do
         WF_L(nf)%Degen = Op_tmp(1, nf)%E(N_part + 1) - Op_tmp(1, nf)%E(N_part)
         WF_R(nf)%Degen = Op_tmp(1, nf)%E(N_part + 1) - Op_tmp(1, nf)%E(N_part)
      end do

      do nf = 1, N_FL
         call WF_overlap(WF_L(nf), WF_R(nf), Z_norm)
      end do

      do nf = 1, N_FL
         call Op_clear(OP_tmp(1, nf), Ndim)
      end do
      deallocate (OP_tmp)
      call Predefined_hoppings_clear(Hopping_Matrix_tmp)

      deallocate (Ham_T_vec, Ham_Tperp_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec)
      deallocate (ham_tx_vec, ham_ty_vec)

   end subroutine Predefined_TrialWaveFunction

end module Predefined_Trial
