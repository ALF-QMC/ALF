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

   subroutine Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
        &                                  N_part, N_FL, N_slat, WF_L, WF_R)

      implicit none
      character(len=64), intent(IN)                :: Lattice_type
      integer, intent(IN)                           :: Ndim, N_FL, N_slat, N_part
      integer, intent(IN), dimension(:, :)           :: List, Invlist
      type(Lattice), intent(in)                   :: Latt
      type(Unit_cell), intent(in)                   :: Latt_Unit
      type(WaveFunction), intent(out), dimension(:, :), allocatable :: WF_L, WF_R

      type(operator), dimension(:, :), allocatable  :: OP_tmp
      type(Hopping_Matrix_type), allocatable       :: Hopping_Matrix_tmp(:)
      real(Kind=kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer, mass
      real(Kind=kind(0.d0))                        :: ham_U, sgn_i, sgn_updn, rmu, stag_sgn, stag_mass
      logical                                       :: Checkerboard, Symm, Kekule_Trial, hatree_fock, l_cmplx_trial

      type(Lattice)                                :: Latt_Kekule
      real(Kind=kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2), x1_p(2), hop(3), del_p(2)
      real(Kind=kind(0.d0))  :: delta = 0.01, Ham_T1, Ham_T2, Ham_Tperp, v1_eff, v2_eff, dtmp

      integer :: N, nf, I, I1, I2, nc, nc1, IK_u, I_u, J1, lp, J, N_Phi, ns, no, k1, N_part_tot, Np_arr(2), i0, j0
      integer :: j2, k2, ix, l_width, ierr, nx, ny, nn_bond(4,2), nnn_bond(4,2)
      logical :: Test = .false., Bulk = .true., lconf
      complex(Kind=kind(0.d0)) :: Z_norm

      real(Kind=kind(0.d0)), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
             & ham_tx_vec(:), ham_ty_vec(:), Ham_T2_vec(:), Ham_lambda_vec(:), Ham_T3_vec(:)
      integer, allocatable ::   N_Phi_vec(:)
      real(Kind=kind(0.d0)), allocatable :: ni_in(:)
      character(len=64) :: file_tg

      allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
           &    ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_lambda_vec(N_FL), N_Phi_vec(N_FL), Ham_T3_vec(N_FL))

      allocate (wf_l(n_fl, n_slat), wf_r(n_fl, n_slat))
      do ns = 1, n_slat
      do n = 1, n_fl
         call wf_alloc(wf_l(n, ns), ndim, n_part)
         call wf_alloc(wf_r(n, ns), ndim, n_part)
      end do
      end do

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
      Ham_T3_vec = 0.d0
      Ham_Chem_vec = Ham_Chem
      Ham_lambda_vec = 0.d0

      nx = int(latt%l1_p(1)/latt%a1_p(1))
      ny = int(latt%l2_p(2)/latt%a2_p(2))
      file_tg = "u_eff_in.dat"
      inquire (file=file_tg, exist=lconf)
      allocate(ni_in(ndim))
      !!if (lconf) then
      !!   open (unit=5, file=file_tg, status='old', action='read', iostat=ierr)
      !!   read (5, *) v1_eff, v2_eff
      !!   i1 = 1
      !!   do i = 1, nx
      !!   do j = 1, ny
      !!       !! a sublattice
      !!      read (5, *) dtmp
      !!      i2 = invlist(i1,1)
      !!      ni_in(i2) = dtmp
      !!      
      !!       !! b sublattice
      !!      read (5, *) dtmp
      !!      i2 = invlist(i1,2)
      !!      ni_in(i2) = dtmp
      !!      
      !!      i1 = latt%nnlist(i1, 0, 1)
      !!   enddo
      !!   i1 = latt%nnlist(i1, 1, 0)
      !!   enddo
      !!   close(5)
      !!endif

      select case (Lattice_type)

      case ("Pi_Flux_ob")
          !! real trial
         Ham_T_vec = 1.d0
         Ham_T2_vec = 0.5d0
         call set_hopping_parameters_pi_flux_qbt_ob(Hopping_Matrix_tmp, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, &
             & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)

          !! cmplx trial
         !!ham_t3_vec = 0.0d0
         !!call Set_hopping_parameters_pi_flux_ob(Hopping_Matrix_tmp, Ham_T_vec, Ham_T2_vec, Ham_T3_vec, Ham_Chem_vec, &
         !!    & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
      
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

      case ("Pi_Flux")
          !! real trial
         Ham_T_vec = 1.d0
         Ham_T2_vec = 0.5d0
         !Phi_Y_vec = 0.01
         call set_hopping_parameters_pi_flux_qbt(Hopping_Matrix_tmp, Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, &
             & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
      
         call Predefined_Hoppings_set_OPT(Hopping_Matrix_tmp, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_tmp)
      
         !!if (lconf) then
         !!   do nf = 1, N_FL
         !!   do I = 1, latt%N

         !!       i1 = invlist(i, 1) !! a sublattice
         !!       i2 = invlist(i, 2) !! b sublattice

         !!       nn_bond(:,2) = i2; 
         !!       nn_bond(1,1) = i1; 
         !!       nn_bond(2,1) = invlist(latt%nnlist(i,1,0),1);
         !!       nn_bond(3,1) = invlist(latt%nnlist(i,0,1),1); 
         !!       nn_bond(4,1) = invlist(latt%nnlist(i,1,1),1);

         !!       do k1 = 1, 4
         !!          op_tmp(1,nf)%o(nn_bond(k1,1),nn_bond(k1,1)) = & 
         !!              & op_tmp(1,nf)%o(nn_bond(k1,1),nn_bond(k1,1)) & 
         !!              & + v1_eff*(ni_in(nn_bond(k1,2))-0.5d0)
         !!          
         !!          op_tmp(1,nf)%o(nn_bond(k1,2),nn_bond(k1,2)) = & 
         !!              & op_tmp(1,nf)%o(nn_bond(k1,2),nn_bond(k1,2)) & 
         !!              & + v1_eff*(ni_in(nn_bond(k1,1))-0.5d0)
         !!       enddo
         !!       
         !!       nnn_bond(1,1) = i1; nnn_bond(1,2) = invlist(latt%nnlist(i,1,0),1);
         !!       nnn_bond(2,1) = i2; nnn_bond(2,2) = invlist(latt%nnlist(i,1,0),2);
         !!       nnn_bond(3,1) = i1; nnn_bond(3,2) = invlist(latt%nnlist(i,0,1),1);
         !!       nnn_bond(4,1) = i2; nnn_bond(4,2) = invlist(latt%nnlist(i,0,1),2);
         !!       
         !!       do k1 = 1, 4
         !!          op_tmp(1,nf)%o(nnn_bond(k1,1),nnn_bond(k1,1)) = & 
         !!              & op_tmp(1,nf)%o(nnn_bond(k1,1),nnn_bond(k1,1)) & 
         !!              & + v2_eff*(ni_in(nnn_bond(k1,2))-0.5d0)
         !!          
         !!          op_tmp(1,nf)%o(nnn_bond(k1,2),nnn_bond(k1,2)) = & 
         !!              & op_tmp(1,nf)%o(nnn_bond(k1,2),nnn_bond(k1,2)) & 
         !!              & + v2_eff*(ni_in(nnn_bond(k1,1))-0.5d0)
         !!       enddo

         !!   enddo
         !!   enddo

         !!else

            !! add stagger mass to avoid the degeneracy of qbt

            stag_mass = 0.001
            do nf = 1, N_FL
               do I = 1, latt%N
                  do no = 1, Latt_unit%norb
                     stag_sgn = -1.d0
                     if (mod(no, 2) .eq. 0) stag_sgn = 1.d0
                     I1 = invlist(I, no)
                     !! onsite sublattice mass
                     op_tmp(1, nf)%o(I1, I1) = stag_sgn*stag_mass
                  end do
               end do
            end do
            
            !!l_width = int(latt%l2_p(2)/latt%a2_p(2))
            !!stag_mass = 0.005
            !!do nf = 1, N_FL
            !!   I = 1
            !!   do J = 1, 1!l_width
            !!      do no = 1, Latt_unit%norb
            !!         stag_sgn = 1.d0
            !!         if (mod(no, 2) .eq. 0) stag_sgn = -1.d0
            !!         I1 = invlist(I, no)
            !!         !! onsite sublattice mass
            !!         op_tmp(1, nf)%o(I1, I1) = stag_sgn*stag_mass
            !!      end do
            !!      I = latt%nnlist(I,0,1)
            !!   end do
            !!end do

            !!!! pinning field
            !!stag_mass = 0.005
            !!do nf = 1, N_FL
            !!   I = 1
            !!   I = latt%nnlist(I,1,1)
            !!   do J = 1, 1!l_width
            !!      do no = 1, Latt_unit%norb
            !!          I1 = invlist(I, no)
            !!          J1 = invlist(latt%nnlist(I,1,0), no)
            !!          K1 = invlist(latt%nnlist(I,0,1), no)
            !!          !! Hopping amplitude
            !!          stag_sgn = -1.d0
            !!          if (mod(no, 2) .eq. 0) stag_sgn = 1.d0
            !!          op_tmp(1, nf)%o(I1, J1) = op_tmp(1, nf)%o(I1, J1) + &
            !!              & stag_sgn*stag_mass
            !!          op_tmp(1, nf)%o(J1, I1) = op_tmp(1, nf)%o(J1, I1) + &
            !!              & stag_sgn*stag_mass
            !!          op_tmp(1, nf)%o(I1, K1) = op_tmp(1, nf)%o(I1, K1) + &
            !!              & stag_sgn*stag_mass
            !!          op_tmp(1, nf)%o(K1, I1) = op_tmp(1, nf)%o(K1, I1) + &
            !!              & stag_sgn*stag_mass
            !!      enddo
            !!      I = latt%nnlist(I,0,1)
            !!   enddo
            !!enddo

         !!endif

      case default
         write (error_unit, *) 'No predefined trial wave function for this lattice.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end select

      do nf = 1, N_FL
         call diag(op_tmp(1, nf)%o, op_tmp(1, nf)%u, op_tmp(1, nf)%e)
         do I2 = 1, N_part
            do I1 = 1, Ndim
               wf_l(nf, 1)%p(I1, I2) = op_tmp(1, nf)%u(I1, I2)
               wf_r(nf, 1)%p(I1, I2) = op_tmp(1, nf)%u(I1, I2)
            end do
         end do
         wf_l(nf, 1)%degen = op_tmp(1, nf)%e(N_part + 1) - op_tmp(1, nf)%e(N_part)
         wf_r(nf, 1)%degen = op_tmp(1, nf)%e(N_part + 1) - op_tmp(1, nf)%e(N_part)
      end do

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
      deallocate (ham_tx_vec, ham_ty_vec, ham_t3_vec)
      deallocate (ni_in)

   end subroutine Predefined_TrialWaveFunction

end module Predefined_Trial
