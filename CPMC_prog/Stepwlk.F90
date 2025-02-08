module stepwlk_mod
   use Hamiltonian_main
   use Operator_mod
   use Control
   use Hop_mod
   use UDV_State_mod
   use gfun_mod
   use upgrade_mod

contains

   subroutine stepwlk_move(phi_trial, phi_0, gr, kappa, kappa_bar, ntau_bp)

      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(in   ) :: phi_trial
      class(udv_state), dimension(:, :), allocatable, intent(inout) :: phi_0
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(inout) :: gr
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(inout) :: kappa, kappa_bar
      integer, intent(in) :: ntau_bp

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op
      real(kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(kind=kind(0.d0)) :: zero = 1.0e-12, sign_w

      do i_wlk = 1, N_wlk
         ! update weight by fac_norm
         sign_w = cos(aimag(weight_k(i_wlk)))
         if (sign_w .gt. zero) weight_k(i_wlk) = weight_k(i_wlk) + fac_norm
      end do

      call half_K_propagation(phi_trial, phi_0)

          !! propagate with interaction
      call int_propagation(phi_0, ntau_bp)

          !! Kinetic part exp(-/Delta/tau T/2)
      call half_K_propagation(phi_trial, phi_0)
      
          !! update weight, Green's function and force bias after propagation 
      call update_weight_and_overlap(phi_trial, phi_0, gr, kappa, kappa_bar, ntau_bp)

   end subroutine stepwlk_move

   subroutine int_propagation(phi_0, ntau_bp)

      use random_wrap
      implicit none

      class(udv_state), intent(inout), allocatable :: phi_0(:, :)
      integer, intent(in) :: ntau_bp

      ! local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed
      integer :: n_op, ns, i_grc
      complex(Kind=kind(0.d0)) :: Z, z_alpha, cspin
      complex(Kind=kind(0.d0)) :: det_Vec(N_FL)
      real(Kind=kind(0.d0)) :: sign_w, spin, zero = 1.0e-12

      n_op = size(op_v, 1)

      do i_wlk = 1, n_wlk

         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then

            ! upgrade Green's function
            do n = 1, n_op

                !! sample gaussian field
               spin = rang_wrap()
               cspin = spin - x_local(n,i_wlk)
               nsigma_bp(i_wlk)%f(n, ntau_bp) = cspin

               ! propagate slater determinant
               do nf = 1, N_FL
                  call op_mmultr(phi_0(nf, i_wlk)%u, op_v(n, nf), nsigma_bp(i_wlk)%f(n, ntau_bp), 'n', 1)
               end do

            end do

         end if

      end do

   end subroutine int_propagation

   subroutine half_K_propagation(phi_trial, phi_0)

      implicit none

      class(udv_state), intent(in), allocatable :: phi_trial(:, :)
      class(udv_state), intent(inout), allocatable :: phi_0(:, :)

      ! local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed, ns, i_grc
      complex(Kind=kind(0.d0)) :: overlap_old, overlap_new, Z, sum_o_new, sum_o_old
      real(Kind=kind(0.d0)) :: overlap_ratio, re_overlap, re_o_max
      real(Kind=kind(0.d0)) :: zero = 1.0e-12, sign_w

      do i_wlk = 1, N_wlk

         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then

             !! multi exp(-\Delta\tau T/2)
            do nf = 1, n_fl
               call Hop_mod_mmthr_1D2(phi_0(nf, i_wlk)%u, nf, 1)
            end do

         end if

      end do

   end subroutine half_K_propagation
      
   subroutine update_weight_and_overlap(phi_trial, phi_0, gr, kappa, kappa_bar, ntau_bp)

      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(in   ) :: phi_trial
      class(udv_state), dimension(:, :), allocatable, intent(inout) :: phi_0
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(inout) :: gr
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(inout) :: kappa, kappa_bar
      integer, intent(in) :: ntau_bp
      
      !Local
      integer :: nf, n_type, ntau1, n, ns, m, nt, NVAR, i_wlk, N_op, i_st, i_ed, i_grc
      real(kind=kind(0.d0)) :: zero = 1.0e-12, sign_w, spin, hs_new, costheta, pi = acos(-1.d0)
      real(Kind=kind(0.d0)) :: re_overlap, re_o_max, logcostheta
      complex(kind=kind(0.d0)) :: gauss_spin, detO, sum_o_new, sum_o_old, s_d_hs, x_bar, overlap_ratio
      complex(Kind=kind(0.d0)) :: det_Vec(n_fl), log_o_new(n_hfb), log_o_old(n_hfb), c_log_I, z_alpha

      n_op = size(op_v, 1)

      do i_wlk = 1, N_wlk

         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then

             re_o_max = 0.d0
             
             !! < BCS | phi^n_k >
             do ns = 1, n_hfb
                i_grc = ns + (i_wlk - 1)*n_hfb
                log_o_old(ns) = overlap(i_grc)
                if ( dble(log_o_old(ns)) .gt. re_o_max ) re_o_max = dble(log_o_old(ns))
             enddo

             !! < BCS | phi^{n+1}_k >
             !! upgrade Green's function
             do ns = 1, N_hfb
                i_grc = ns + (i_wlk - 1)*N_hfb

                call cgrp(detO, gr(:, :, :, i_grc), kappa(:, :, :, i_grc), kappa_bar(:, :, :, i_grc), &
                    &     phi_0(1, i_wlk), phi_0(2, i_wlk), phi_trial(1, ns))
                overlap(i_grc) = detO*n_sun
                log_o_new(ns) = overlap(i_grc)
                if ( dble(log_o_new(ns)) .gt. re_o_max ) re_o_max = dble(log_o_new(ns))
            enddo

            !! < BCS | phi^{n+1}_k >/< BCS | phi^{n}_k >
            sum_o_old = 0.d0
            sum_o_new = 0.d0
            do ns = 1, N_hfb
               sum_o_old = sum_o_old + exp(log_o_old(ns)-cmplx(re_o_max,0.d0,kind(0.d0)))
               sum_o_new = sum_o_new + exp(log_o_new(ns)-cmplx(re_o_max,0.d0,kind(0.d0)))
            enddo
            overlap_ratio = sum_o_new/sum_o_old
            
            !! \prod_i \frac{p(x(i)-\bar{x}(i))}{p(x(i))}
            s_d_hs = 0.d0
            do n = 1, n_op
               x_bar = nsigma_bp(i_wlk)%phi(n,ntau_bp) + x_local(n,i_wlk)
               s_d_hs = s_d_hs + x_bar*x_local(n,i_wlk) - x_local(n,i_wlk)*x_local(n,i_wlk)/2.d0
            enddo
            
            !! alpha factor from operators
            z_alpha = 0.d0
            do n = 1, n_op
               do nf = 1, n_fl
                  z_alpha = z_alpha + &
                      & op_v(n, nf)%g*op_v(n, nf)%alpha*nsigma_bp(i_wlk)%phi(n, ntau_bp)
               enddo
            enddo

            !! logarithmic of I = < BCS | phi^{n+1}_k >/< BCS | phi^{n}_k >*\prod_i \frac{p(x(i)-\bar{x}(i))}{p(x(i))}
            c_log_I = log(overlap_ratio) + s_d_hs + z_alpha
            costheta = cos(aimag(c_log_I))
            if ( costheta .gt. zero ) then
               logcostheta = log(costheta)
               weight_k(i_wlk) = weight_k(i_wlk) + dble(c_log_I) + logcostheta
            else
               weight_k(i_wlk) = cmplx(0.d0,pi,kind(0.d0))
            endif

         endif

      enddo

      call ham%set_xloc(gr, kappa, kappa_bar)
   
   end subroutine update_weight_and_overlap

   subroutine re_orthonormalize_walkers(Phi_0)
      
      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(INOUT) :: phi_0

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, N_size, I, i_wlk
      integer :: ndistance, i_st, i_ed, ns, nrs, i_slat, i_grc
      real(Kind=kind(0.d0)) :: Overlap_ratio, zero = 1.0e-12, pi = acos(-1.d0)
      real(Kind=kind(0.d0)) :: re_o_am, re_o_ph, sign_w
      complex(Kind=kind(0.d0)) :: det_D(N_FL)

      do i_wlk = 1, n_wlk

         sign_w = cos(aimag(weight_k(i_wlk)))
         if (sign_w .gt. zero) then

            !Carry out U,D,V decomposition.
            do nf = 1, N_FL
               call Phi_0(nf, i_wlk)%decompose
            end do

            Det_D = cmplx(0.d0, 0.d0, kind(0.d0))

            do nf = 1, N_FL
               N_size = phi_0(nf, 1)%n_part
               do I = 1, N_size
                  det_D(nf) = det_D(nf) + log(Phi_0(nf, i_wlk)%D(I))
               end do
            end do

            !! update overlap
            do ns = 1, n_hfb
                i_grc = ns + (i_wlk - 1)*N_hfb
                overlap(i_grc) = overlap(i_grc) - sum(det_D(:))
            enddo

            do nf = 1, N_FL
               Phi_0(nf, i_wlk)%D(:) = cmplx(1.d0, 0.d0, kind(0.d0))
            end do

         end if

      end do

   end subroutine re_orthonormalize_walkers
   
   subroutine population_control(phi_0, phi_bp_r)
#ifdef MPI
      use mpi
#endif
      use Random_wrap

      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(inout) :: phi_0, phi_bp_r

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, it_wlk, n_exc, pop_exc(N_wlk_mpi, 4)
      integer :: j, it, i_t, i_st, i_ed, nu_wlk, i_src, i_wlk, j_src, j_wlk, n1, n2, nrg, nfrg, ilabel, ncslat
      integer :: i_llim, i_rlim, j_llim, j_rlim, n_part, ii
      real(Kind=kind(0.d0)) :: zero = 1.0e-12, d_scal, sum_w, w_count, re_weight(N_wlk_mpi), ang_w, re_lw
      real(Kind=kind(0.d0)) :: exp_o_abs(n_hfb), exp_o_phase(n_hfb), dz2, max_re_w
      complex(Kind=kind(0.d0)), allocatable :: w_arr(:), weight_mpi(:)
      complex(Kind=kind(0.d0)) :: overlap_tmp(n_grc)
      complex(Kind=kind(0.d0)) :: Z1, Z2, Z3, Z_s_array(N_hfb), Z_r_array(N_hfb), zp
      type(Fields), dimension(:), allocatable :: nsigma_store
      class(UDV_State), dimension(:, :), allocatable :: phi_0_m, phi_bp_m

#ifdef MPI
      integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

          !! temporary store
          !! store slater determinant
      pop_exc = 0
      allocate (phi_0_m (N_FL, N_wlk))
      allocate (phi_bp_m(N_FL, N_wlk))
      do i_wlk = 1, N_wlk
         do nf = 1, N_FL
            n_part = phi_0(nf, 1)%n_part
            call phi_0_m (nf, i_wlk)%alloc(ndim, n_part)
            call phi_bp_m(nf, i_wlk)%alloc(ndim, n_part)
         end do
      end do

      do nf = 1, N_FL
         do i_wlk = 1, N_wlk
            phi_0_m (nf, i_wlk) = phi_0   (nf, i_wlk)
            phi_bp_m(nf, i_wlk) = phi_bp_r(nf, i_wlk)
         end do
      end do

          !! store fields
      n1 = size(nsigma_bp(1)%f, 1)
      n2 = size(nsigma_bp(1)%f, 2)
      allocate (nsigma_store(N_wlk))
      do i_wlk = 1, N_wlk
         call nsigma_store(i_wlk)%make(n1, n2)
         nsigma_store(i_wlk)%f = nsigma_bp(i_wlk)%f
         nsigma_store(i_wlk)%t = nsigma_bp(i_wlk)%t
      end do

          !! store overlap
      overlap_tmp(:) = overlap(:)

      if ( irank_g .eq. 0 ) then
          allocate(weight_mpi(N_wlk_mpi))
          allocate(w_arr     (N_wlk    ))
          i_st = 1
          i_ed = N_wlk
          weight_mpi(i_st:i_ed) = weight_k(:)
      endif

      if ( irank_g .eq. 0 ) then
          do ii = 1, isize_g - 1
             i_st = ii*N_wlk + 1
             i_ed = (ii + 1)*N_wlk
             call mpi_recv(w_arr, n_wlk, MPI_COMPLEX16, ii, 1, group_comm, status, ierr)
             weight_mpi(i_st:i_ed) = w_arr
          enddo
      else
          call mpi_send(weight_k, n_wlk, MPI_COMPLEX16, 0, 1, group_comm, ierr)
      endif

      if ( irank_g .eq. 0 ) then
          max_re_w = maxval(dble(weight_mpi))
          weight_mpi(:) = weight_mpi(:) - max_re_w
          do i_wlk = 1, n_wlk_mpi
              ang_w = aimag(weight_mpi(i_wlk)) 
              re_lw = dble (weight_mpi(i_wlk))
              re_weight(i_wlk) = 0.d0
              if ( cos(ang_w) .gt. 0.d0 ) re_weight(i_wlk) = exp(re_lw)*cos(ang_w)
          enddo
          deallocate(weight_mpi, w_arr)
      endif

      call MPI_BCAST(re_weight, N_wlk_mpi, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      d_scal = dble(N_wlk_mpi)/sum(re_weight)
      if (irank_g == 0) sum_w = -ranf_wrap()
      call MPI_BCAST(sum_w, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      nu_wlk = 0

      n_exc = 0
      do it_wlk = 1, N_wlk_mpi
         i_src = (it_wlk - 1)/N_wlk
         i_wlk = it_wlk - N_wlk*i_src

         sum_w = sum_w + re_weight(it_wlk)*d_scal
         n = ceiling(sum_w); 
         do j = (nu_wlk + 1), n
            j_src = (j - 1)/N_wlk
            j_wlk = j - N_wlk*j_src

            if (j_src .ne. i_src) then
               n_exc = n_exc + 1
               pop_exc(n_exc, 1) = i_src
               pop_exc(n_exc, 2) = i_wlk
               pop_exc(n_exc, 3) = j_src
               pop_exc(n_exc, 4) = j_wlk
            else
               if (irank_g .eq. i_src) then
                  do nf = 1, N_FL
                     phi_0_m (nf, j_wlk) = phi_0   (nf, i_wlk)
                     phi_bp_m(nf, j_wlk) = phi_bp_r(nf, i_wlk)
                  end do
                  nsigma_store(j_wlk)%f = nsigma_bp(i_wlk)%f

                  i_llim = 1 + (i_wlk - 1)*N_hfb; i_rlim = i_wlk*N_hfb
                  j_llim = 1 + (j_wlk - 1)*N_hfb; j_rlim = j_wlk*N_hfb
                  overlap_tmp(j_llim:j_rlim) = overlap(i_llim:i_rlim)
               end if
            end if
         end do
         nu_wlk = n
      end do

      nrg = 2 + N_fl*6

      do it = 1, n_exc
         i_src = pop_exc(it, 1); i_wlk = pop_exc(it, 2)
         j_src = pop_exc(it, 3); j_wlk = pop_exc(it, 4)
         if (irank_g .eq. i_src) then
            i_llim = 1 + (i_wlk - 1)*N_hfb; i_rlim = i_wlk*N_hfb
            Z_s_array(:) = overlap(i_llim:i_rlim)

            ilabel = (it - 1)*nrg
            call mpi_send(Z_s_array, N_hfb, MPI_COMPLEX16, j_src, ilabel, Group_comm, IERR)

            do nf = 1, N_FL
               ilabel = (it - 1)*nrg + (nf - 1)*6 + 1
               call phi_0(nf, i_wlk)%MPI_send_general(j_src, ilabel, ierr)
               ilabel = (it - 1)*nrg + (nf - 1)*6 + 4
               call phi_bp_r(nf, i_wlk)%MPI_send_general(j_src, ilabel, ierr)
            end do

            ilabel = (it - 1)*nrg + (n_fl)*6 + 1
            call mpi_send(nsigma_bp(i_wlk)%f, n1*n2, MPI_REAL8, j_src, ilabel, Group_comm, IERR)
         end if
         if (irank_g .eq. j_src) then
            ilabel = (it - 1)*nrg
            call mpi_recv(Z_r_array, N_hfb, MPI_COMPLEX16, i_src, ilabel, Group_comm, STATUS, IERR)

            j_llim = 1 + (j_wlk - 1)*N_hfb; j_rlim = j_wlk*N_hfb
            overlap_tmp(j_llim:j_rlim) = Z_r_array(:)

            do nf = 1, N_FL
               ilabel = (it - 1)*nrg + (nf - 1)*6 + 1
               call phi_0_m(nf, j_wlk)%MPI_recv_general(i_src, ilabel, status, ierr)
               ilabel = (it - 1)*nrg + (nf - 1)*6 + 4
               call phi_bp_m(nf, j_wlk)%MPI_recv_general(i_src, ilabel, status, ierr)
            end do

            ilabel = (it - 1)*nrg + (n_fl)*6 + 1
            call mpi_recv(nsigma_store(j_wlk)%f, n1*n2, MPI_REAL8, i_src, ilabel, Group_comm, STATUS, IERR)
         end if
      end do

      overlap(:) = overlap_tmp(:)
      ! reset weight
      weight_k(:) = cmplx(0.d0, 0.d0, kind(0.d0))
      do nf = 1, N_FL
         do i_wlk = 1, N_wlk
            phi_0   (nf, i_wlk) = phi_0_m (nf, i_wlk)
            phi_bp_r(nf, i_wlk) = phi_bp_m(nf, i_wlk)
         end do
      end do

      do i_wlk = 1, N_wlk
         nsigma_bp(i_wlk)%f = nsigma_store(i_wlk)%f
         nsigma_bp(i_wlk)%t = nsigma_store(i_wlk)%t
      end do

      ! deallocate tmp udv class
      do i_wlk = 1, N_wlk
         do nf = 1, N_FL
            call phi_0_m (nf, i_wlk)%dealloc
            call phi_bp_m(nf, i_wlk)%dealloc
         end do
      end do

      do i_wlk = 1, N_wlk
         call nsigma_store(i_wlk)%clear
      end do

      deallocate (phi_0_m, phi_bp_m)
      deallocate (nsigma_store)

   end subroutine population_control

   subroutine store_phi(phi_0, phi_bp_r)

      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0, phi_bp_r

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk

      do nf = 1, N_FL
         do i_wlk = 1, N_wlk
            phi_bp_r(nf, i_wlk) = phi_0(nf, i_wlk)
                !! phi_bp_r is not used in propagation
            call phi_bp_r(nf, i_wlk)%decompose
         end do
      end do

   end subroutine store_phi

   subroutine backpropagation(phi_bp_l, phi_bp_r, udvst, stab_nt, ltau)
#ifdef MPI
      use mpi
#endif
      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(inout) :: phi_bp_l, phi_bp_r
      class(udv_state), dimension(:, :, :), allocatable, intent(inout) :: udvst
      integer, dimension(:), allocatable, intent(in) :: stab_nt
      integer, intent(in) :: ltau

      !Local
      complex(Kind=kind(0.d0)) :: gr_bp(NDIM, NDIM, N_FL, n_grc)
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, ltrot_bp, N_op, nstm, nst, ntau
      integer :: i_grc, i_st, i_ed, ns
      complex(Kind=kind(0.d0)) :: z, z_weight, z_sum_overlap, exp_overlap(N_hfb)
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(Kind=kind(0.d0))    :: zero = 1.0e-12, sign_w, ang_w, re_lw, re_we

#ifdef MPI
      integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

!      N_op = size(OP_V, 1)
!      ltrot_bp = size(nsigma_bp(1)%f, 2)
!      nstm = size(udvst, 1)
!
!          !! initialization
!      do i_grc = 1, N_grc
!         i_wlk = (i_grc - 1)/N_hfb + 1     ! index for walker
!         ns = i_grc - (i_wlk - 1)*N_hfb ! index for slater det
!         do nf = 1, N_FL
!            call phi_bp_l(nf, i_grc)%reset('l', wf_l(nf, ns)%P)
!            call udvst(nstm, nf, i_grc)%reset('l', wf_l(nf, ns)%P)
!         end do
!      end do
!
!          !! backpropagation
!      nst = nstm - 1
!      do nt = ltrot_bp, 1, -1
!         ntau = nt - 1
!
!         do i_wlk = 1, N_wlk
!
!            sign_w = cos(aimag(weight_k(i_wlk)))
!            if ( sign_w .gt. zero ) then
!             
!                do ns = 1, N_hfb
!                   i_grc = ns + (i_wlk - 1)*N_hfb
!                   do nf = 1, N_FL
!
!                      call Hop_mod_mmthlc_1D2(phi_bp_l(nf, i_grc)%U, nf, 1)
!
!                      do n = N_op, 1, -1
!                         call Op_mmultR(phi_bp_l(nf, i_grc)%U, Op_V(n, nf), nsigma_bp(i_wlk)%f(n, nt), 'c', 1)
!                      end do
!
!                      call Hop_mod_mmthlc_1D2(phi_bp_l(nf, i_grc)%U, nf, 1)
!
!                   end do
!                end do
!
!            end if
!
!         end do
!
!         if (ntau .eq. stab_nt(nst) .and. ntau .ne. 0) then
!            call re_orthonormalize_walkers(phi_bp_l, 'N')
!            do i_grc = 1, N_grc
!            do nf = 1, N_FL
!               udvst(nst, nf, i_grc) = phi_bp_l(nf, i_grc)
!            end do
!            end do
!            nst = nst - 1
!         end if
!
!      end do
!
!          !! svd at tau = 0
!      call re_orthonormalize_walkers(phi_bp_l, 'N')
!
!          !! compute the total weight
!      call ham%sum_weight(z_weight)
!
!          !! equal time measurement
!      if ( irank .eq. 0 ) call ham%count_obs
!
!      do i_wlk = 1, N_wlk
!         
!         sign_w = cos(aimag(weight_k(i_wlk)))
!         if ( sign_w .gt. zero ) then 
!         
!             ang_w = aimag(weight_k(i_wlk)) 
!             re_lw = dble (weight_k(i_wlk))
!             re_we = exp(re_lw)*cos(ang_w)
!             
!             i_st = 1 + (i_wlk - 1)*N_hfb; i_ed = i_wlk*N_hfb
!             exp_overlap(:) = exp(overlap(i_st:i_ed))
!             z_sum_overlap = sum(exp_overlap(:))
!             do ns = 1, N_hfb
!                i_grc = ns + (i_wlk - 1)*N_hfb
!                do nf = 1, N_Fl
!                   call CGRP(Z, GR_bp(:, :, nf, i_grc), phi_bp_r(nf, i_wlk), phi_bp_l(nf, i_grc))
!                end do
!                call ham%obser(gr_bp(:, :, :, i_grc), gr_mix(:, :, :, i_grc), i_grc, re_we, z_weight, z_sum_overlap)
!             end do
!     
!         endif
!
!      end do

   end subroutine backpropagation

end module stepwlk_mod
