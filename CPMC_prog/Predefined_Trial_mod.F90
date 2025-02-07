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
        &                                  N_part, alpha, N_FL, N_hfb, WF_L, WF_R)

      implicit none
      character(len=64), intent(IN)                :: Lattice_type
      integer, intent(IN)                           :: Ndim, N_FL, N_hfb, N_part
      integer, intent(IN), dimension(:, :)           :: List, Invlist
      real(Kind=kind(0.d0)), intent(in)             :: alpha
      type(Lattice), intent(in)                   :: Latt
      type(Unit_cell), intent(in)                   :: Latt_Unit
      type(WaveFunction), intent(out), dimension(:), allocatable :: wf_l, wf_r

      !! local
      type(operator), dimension(:, :), allocatable  :: op_tmp
      type(Hopping_Matrix_type), allocatable       :: Hopping_Matrix_tmp(:)
      real(Kind=kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, Phi_X, Phi_Y, Dimer, mass
      logical                                       :: Checkerboard, Symm, Kekule_Trial, hatree_fock, l_cmplx_trial
      real(Kind=kind(0.d0))  :: delta = 0.01, Ham_T1, Ham_T2, Ham_Tperp, kx, ky
      real(Kind=kind(0.d0))  :: sc_modu, epsilon_k, xi_k, k_p(2), r_p(2), r1_p(2), dr_p(2), ang
      integer :: N, nf, I, I1, I2, nc, nc1, J1, lp, J, K, N_Phi, ns, no, k1, N_part_tot, i0, j0
      logical :: Test = .false., Bulk = .true.
      complex(Kind=kind(0.d0)) :: Z_norm, ztmp
      real(Kind=kind(0.d0)), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
             & ham_tx_vec(:), ham_ty_vec(:), Ham_T2_vec(:), Ham_lambda_vec(:)
      complex(Kind=kind(0.d0)), allocatable :: uk_bg_vec(:), vk_bg_vec(:), F_ssc(:,:)
      integer, allocatable ::   N_Phi_vec(:)

      allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
           &    ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_lambda_vec(N_FL), N_Phi_vec(N_FL))
       
      allocate(wf_l(n_hfb),wf_r(n_fl))

      !!============================================================================================!!
      !! Set up Slater determinant
      !!============================================================================================!!
      !! Slater determinant
      do n=1,n_fl
         call wf_alloc(wf_r(n),ndim,n_part)
      enddo
      
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
      
      case ('square_anisotropic')
           Allocate(op_tmp(1,n_fl))
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
              
              ham_tx_vec(1) = ham_t;
              ham_tx_vec(2) = ham_t*alpha;
              ham_ty_vec(1) = ham_t*alpha;
              ham_ty_vec(2) = ham_t;
              
              do n = 1, n_fl
                 op_tmp(1,n)%o(i1,j1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(j1,i1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(i1,k1) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
                 op_tmp(1,n)%o(k1,i1) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
              enddo
           enddo

           do n = 1, n_fl
              Call op_set(op_tmp(1,n))
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
               wf_r(nf)%p(I1,I2)=op_tmp(1,nf)%u(I1,I2)
            enddo
         enddo
         wf_r(nf)%degen = op_tmp(1,nf)%e(N_part+1) - op_tmp(1,nf)%e(N_part)
      enddo
      
      !! deallocate tmp matrix
      do nf = 1, N_FL
         call op_clear(op_tmp(1, nf), ndim)
      end do
      deallocate (op_tmp)
      call Predefined_hoppings_clear(Hopping_Matrix_tmp)
      !!============================================================================================!!
      
      !!============================================================================================!!
      !! Set up BCS trial
      !!============================================================================================!!
      !! BCS trial
      !! isotropic s wave trial input
      !! only store Z_{\up\dn} sub block
      do ns=1,n_hfb
         call wf_alloc(wf_l(ns),ndim,ndim)
      enddo
      allocate(uk_bg_vec(latt%n), vk_bg_vec(latt%n))
      allocate(F_ssc(latt%n,latt%n))
    
      do i = 1, latt%n
          sc_modu = 0.5d0
          
          k_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p
          kx = k_p(1); ky = k_p(2)

          epsilon_k = -2.d0*ham_t*(cos(kx)+cos(ky))
          xi_k = sqrt(epsilon_k**2 + sc_modu**2)
          
          !!!u_k
          uk_bg_vec(i) = cmplx(sqrt(0.5d0*(1.d0+epsilon_k/xi_k)), 0.d0, kind(0.d0))
          !!!v_k
          vk_bg_vec(i) = cmplx(sqrt(0.5d0*(1.d0-epsilon_k/xi_k)), 0.d0, kind(0.d0))
      enddo
      
      f_ssc(:,:) = 0.d0
      do i = 1, latt%n
      do j = 1, latt%n
         r_p  = dble(Latt%listk(i,1))*Latt%a1_p + dble(Latt%listk(i,2))*Latt%a2_p
         r1_p = dble(Latt%listk(j,1))*Latt%a1_p + dble(Latt%listk(j,2))*Latt%a2_p
         dr_p = r_p-r1_p

         ztmp = cmplx(0.d0,0.d0,kind(0.d0))
         !!! F_{r,r'}=\frac{1}{L}\sum_k e^{-ik*(r-r')}u_k/v_k
         do k = 1, latt%n
            k_p = dble(Latt%listk(k,1))*Latt%b1_p + dble(Latt%listk(k,2))*Latt%b2_p
            ang = -(k_p(1)*dr_p(1)+k_p(2)*dr_p(2))
            ztmp = ztmp + cmplx(cos(ang), sin(ang), kind(0.d0))*vk_bg_vec(k)/uk_bg_vec(k)
         enddo
         ztmp = ztmp/dble(latt%n)
         F_ssc(i,j) = ztmp
      enddo
      enddo

      !! In our model ndim = latt%n
      do ns=1,n_hfb
         wf_l(ns)%p(:,:) = F_ssc(:,:)
      enddo

      deallocate(F_ssc, uk_bg_vec, vk_bg_vec)
      
      !!============================================================================================!!

      call wf_overlap(wf_l, wf_r, z_norm)

      deallocate (Ham_T_vec, Ham_Tperp_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec)
      deallocate (ham_tx_vec, ham_ty_vec)

   end subroutine Predefined_TrialWaveFunction

end module Predefined_Trial
