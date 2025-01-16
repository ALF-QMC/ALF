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
   use Random_wrap
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
        &                                  N_part, alpha, N_FL, N_slat, WF_L, WF_R)

      implicit none
      character(len=64), intent(IN)                :: Lattice_type
      integer, intent(IN)                           :: Ndim, N_FL, N_slat, N_part
      integer, intent(IN), dimension(:, :)           :: List, Invlist
      real(Kind=kind(0.d0)), intent(in)             :: alpha
      type(Lattice), intent(in)                   :: Latt
      type(Unit_cell), intent(in)                   :: Latt_Unit
      type(WaveFunction), intent(out), dimension(:, :), allocatable :: WF_L, WF_R

      type(operator), dimension(:, :), allocatable  :: OP_tmp
      type(Hopping_Matrix_type), allocatable       :: Hopping_Matrix_tmp(:)
      real(Kind=kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer, mass
      real(Kind=kind(0.d0))                        :: ham_U, sgn_i, sgn_updn, rmu
      logical                                       :: Checkerboard, Symm, Kekule_Trial, hatree_fock, l_cmplx_trial

      type(Lattice)                                :: Latt_Kekule
      real(Kind=kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2), x1_p(2), hop(3), del_p(2)
      real(Kind=kind(0.d0))  :: delta = 0.01, Ham_T1, Ham_T2, Ham_Tperp

      integer :: N, nf, I, I1, I2, nc, nc1, IK_u, I_u, J1, lp, J, N_Phi, ns, no, k1, N_part_tot, Np_arr(2), i0, j0
      integer :: j2, k2
      logical :: Test = .false., Bulk = .true.
      complex(Kind=kind(0.d0)) :: Z_norm

      real(Kind=kind(0.d0)), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
             & ham_tx_vec(:), ham_ty_vec(:), Ham_T2_vec(:), Ham_lambda_vec(:)
      integer, allocatable ::   N_Phi_vec(:)
      real(Kind=kind(0.d0)), allocatable :: eig_sort_arr(:,:)

      allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
           &    ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_lambda_vec(N_FL), N_Phi_vec(N_FL))

      Checkerboard = .false.
      Kekule_Trial = .false.
      hatree_fock = .false.
      l_cmplx_trial = .false.
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
      Ham_lambda_vec = 0.d0

      select case (Lattice_type)
      
      case ('bilayer_square')
           allocate(op_tmp(1,n_fl))
           do n = 1,n_fl
              call op_make(op_tmp(1,n),ndim)
              do i = 1,ndim
                 op_tmp(1,n)%P(i) = i
              enddo
              op_tmp(1,n)%g    = cmplx(1.d0,0.d0,kind(0.d0))
              op_tmp(1,n)%alpha= cmplx(0.d0,0.d0,kind(0.D0))
           enddo

           do i = 1,latt%n
              delta = 0.005d0*ranf_wrap()
              ham_t = 1.d0+delta

              i1 = invlist(i,1)
              j1 = invlist(latt%nnlist(i,1,0),1)
              k1 = invlist(latt%nnlist(i,0,1),1)
              
              i2 = invlist(i,2)
              j2 = invlist(latt%nnlist(i,1,0),2)
              k2 = invlist(latt%nnlist(i,0,1),2)
              
              ham_tx_vec(1) = ham_t;
              ham_ty_vec(1) = ham_t*alpha;
              
              do n = 1, n_fl
                 op_tmp(1,n)%o(i1,j1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(j1,i1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(i1,k1) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(k1,i1) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
                 
                 op_tmp(1,n)%o(i2,j2) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(j2,i2) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(i2,k2) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(k2,i2) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
              enddo
           enddo

           do n = 1, n_fl
              call op_set(op_tmp(1,n))
           enddo

      case default
         write (error_unit, *) 'No predefined trial wave function for this lattice.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end select

      if (l_cmplx_trial )   &
      &     call Predefined_Hoppings_set_OPT(Hopping_Matrix_tmp, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_tmp)

      do nf = 1,N_FL
         call diag(op_tmp(1,nf)%o,op_tmp(1,nf)%u,op_tmp(1,nf)%e)
         do I2=1,N_part
            do I1=1,Ndim
               wf_l(nf,1)%p(I1,I2)=op_tmp(1,nf)%u(I1,I2)
               wf_r(nf,1)%p(I1,I2)=op_tmp(1,nf)%u(I1,I2)
            enddo
         enddo
         wf_l(nf,1)%degen = op_tmp(1,nf)%e(N_part+1) - op_tmp(1,nf)%e(N_part)
         wf_r(nf,1)%degen = op_tmp(1,nf)%e(N_part+1) - op_tmp(1,nf)%e(N_part)
      enddo

      do nf = 1, n_fl
         call wf_overlap(wf_l(nf, 1), wf_r(nf, 1), z_norm)
         wf_l(nf, 1)%p(:, :) = wf_l(nf, 1)%p(:, :)/sqrt(dble(N_slat))
         wf_r(nf, 1)%p(:, :) = wf_r(nf, 1)%p(:, :)/sqrt(dble(N_slat))
      end do

      do nf = 1, n_fl
      do ns = 2, n_slat
         wf_l(nf, ns)%degen = wf_l(nf, 1)%degen
         wf_r(nf, ns)%degen = wf_r(nf, 1)%degen
         wf_l(nf, ns)%p(:, :) = wf_l(nf, 1)%p(:, :)
         wf_r(nf, ns)%p(:, :) = wf_r(nf, 1)%p(:, :)
      end do
      end do

      do nf = 1, N_FL
         call op_clear(op_tmp(1, nf), ndim)
      end do
      deallocate (op_tmp)
      call Predefined_hoppings_clear(Hopping_Matrix_tmp)

      deallocate (Ham_T_vec, Ham_Tperp_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec)
      deallocate (ham_tx_vec, ham_ty_vec)

   end subroutine Predefined_TrialWaveFunction

end module Predefined_Trial
