module stepwlk_mod
   use Hamiltonian_main
   use Operator_mod
   use Control
   use Hop_mod
   use UDV_State_mod
   use gfun_mod
   use upgrade_mod

contains

   subroutine stepwlk_move(phi_trial, phi_0, GR, ntau_bp)

      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(IN)    :: phi_trial
      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(INOUT) :: GR
      integer, intent(IN) :: ntau_bp

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(Kind=kind(0.d0)) :: zero = 1.0e-12, sign_w

      N_op = size(OP_V, 1)

      do i_wlk = 1, N_wlk
         ! update weight by fac_norm
         sign_w = cos(aimag(weight_k(i_wlk)))
         if (sign_w .gt. zero) weight_k(i_wlk) = weight_k(i_wlk) + fac_norm
      end do

      call half_K_propagation(phi_trial, phi_0, GR)

          !! propagate with interaction
      call int_propagation(phi_0, GR, ntau_bp)

          !! Kinetic part exp(-/Delta/tau T/2)
      call half_K_propagation(phi_trial, phi_0, GR)

          !! rescale overlap after each step
      call rescale_overlap(overlap)

   end subroutine stepwlk_move

   subroutine int_propagation(phi_0, GR, ntau_bp)

      implicit none

      class(udv_state), intent(inout), allocatable :: phi_0(:, :)
      complex(Kind=kind(0.d0)), intent(inout), allocatable :: GR(:, :, :, :)
      integer, intent(in) :: ntau_bp

      ! local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed
      integer :: n_op, ns, i_grc
      complex(Kind=kind(0.d0)) :: Z, z_alpha
      complex(Kind=kind(0.d0)) :: det_Vec(N_FL)
      real(Kind=kind(0.d0)) :: sign_w, spin, zero = 1.0e-12

      N_op = size(OP_V, 1)

      do i_wlk = 1, N_wlk

         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then

            ! upgrade Green's function
            z_alpha = 0.d0
            do n = 1, N_op

               N_type = 2
               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf = 1, N_FL
                     call op_wrapdo(gr(:, :, nf, i_grc), op_v(n, nf), 1.d0, ndim, n_Type, 1)
                  end do
               end do

               call upgrade(gr, n, spin, i_wlk)
               nsigma_bp(i_wlk)%f(n, ntau_bp) = spin

               N_type = 2
               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf = 1, N_FL
                     call op_wrapup(gr(:, :, nf, i_grc), op_v(n, nf), 1.d0, ndim, n_Type, 1)
                  end do
               end do

               ! propagate slater determinant
               do nf = 1, N_FL
                  call op_mmultr(phi_0(nf, i_wlk)%u, op_v(n, nf), spin, 'n', 1)
                  ! store alpha factor in z_alpha
                  z_alpha = z_alpha + &
                      & op_v(n_op, nf)%g*op_v(n_op, nf)%alpha*nsigma_bp(i_wlk)%Phi(n, ntau_bp)
               end do

            end do

            !!! comment this since it would cause instability when z_alpha is real 
            !!! rescale U matrix with z_alpha factor
            !!do nf = 1, N_FL
            !!   phi_0(nf, i_wlk)%U(:, :) = exp(z_alpha)*phi_0(nf, i_wlk)%U(:, :)
            !!end do

         end if

      end do

   end subroutine int_propagation

   subroutine half_K_propagation(phi_trial, phi_0, GR)

      implicit none

      class(udv_state), intent(in), allocatable :: phi_trial(:, :)
      class(udv_state), intent(inout), allocatable :: phi_0(:, :)
      complex(Kind=kind(0.d0)), intent(inout), allocatable :: GR(:, :, :, :)

      ! local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed, ns, i_grc
      complex(Kind=kind(0.d0)) :: Overlap_old, Overlap_new, Z, sum_o_new, sum_o_old
      complex(Kind=kind(0.d0)) :: det_Vec(N_FL), log_o_new(n_slat), log_o_old(n_slat)
      real(Kind=kind(0.d0)) :: overlap_ratio, re_overlap, re_o_max
      real(Kind=kind(0.d0)) :: zero = 1.0e-12, sign_w, pi = acos(-1.d0)

      do i_wlk = 1, N_wlk

         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then

            re_o_max = 0.d0

            do ns = 1, N_slat
               i_grc = ns + (i_wlk - 1)*N_slat
               ! Update Green's function
               do nf = 1, N_Fl
                  call CGRP(Z, GR(:, :, nf, i_grc), phi_0(nf, i_wlk), phi_trial(nf, ns))
                  det_vec(nf) = Z
               end do
               det_Vec(:) = det_Vec(:)*N_SUN
               log_o_old(ns) = sum(det_Vec)
               if ( dble(log_o_old(ns)) .gt. re_o_max ) re_o_max = dble(log_o_old(ns))
            end do

             !! multi exp(-\Delta\tau T/2)
            do nf = 1, N_FL
               call Hop_mod_mmthr_1D2(phi_0(nf, i_wlk)%U, nf, 1)
            end do

            do ns = 1, N_slat
               i_grc = ns + (i_wlk - 1)*N_slat
               ! Update Green's function
               do nf = 1, N_Fl
                  call CGRP(Z, GR(:, :, nf, i_grc), phi_0(nf, i_wlk), phi_trial(nf, ns))
                  det_vec(nf) = Z
               end do
               det_Vec(:) = det_Vec(:)*N_SUN
               log_o_new(ns) = sum(det_Vec)
               if ( dble(log_o_new(ns)) .gt. re_o_max ) re_o_max = dble(log_o_new(ns))
            end do

            sum_o_old = 0.d0
            sum_o_new = 0.d0
            do ns = 1, N_slat
               sum_o_old = sum_o_old + exp(log_o_old(ns)-cmplx(re_o_max,0.d0,kind(0.d0)))
               sum_o_new = sum_o_new + exp(log_o_new(ns)-cmplx(re_o_max,0.d0,kind(0.d0)))
            enddo

            overlap_ratio = sum_o_new/sum_o_old
            re_overlap = dble(overlap_ratio)
            if (re_overlap .gt. zero) then
                 !! upgrade overlap
               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  overlap(i_grc) = overlap(i_grc) + &
                      & (log_o_new(ns) - log_o_old(ns))
               end do
               weight_k(i_wlk) = weight_k(i_wlk) + log(re_overlap)
            else
               weight_k(i_wlk) = cmplx(0.d0,pi,kind(0.d0))
            end if

         end if

      end do

   end subroutine half_K_propagation

   subroutine re_orthonormalize_walkers(Phi_0, cop)

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0
      character, intent(IN)    :: cop

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, N_size, I, i_wlk
      integer :: ndistance, i_wlk_eff, i_st, i_ed, n_wlk_eff, ns, nrs, i_slat
      real(Kind=kind(0.d0)) :: Overlap_ratio, zero = 1.0e-12, pi = acos(-1.d0)
      real(Kind=kind(0.d0)) :: re_o_am, re_o_ph
      complex(Kind=kind(0.d0)) :: det_D(N_FL)

      n_wlk_eff = size(phi_0, 2)
      ndistance = n_wlk_eff/n_wlk

      do i_wlk_eff = 1, n_wlk_eff

         i_wlk = (i_wlk_eff - 1)/ndistance + 1     ! index for walker
         ns = i_wlk_eff - (i_wlk - 1)*ndistance ! index for slater det

         sign_w = cos(aimag(weight_k(i_wlk)))
         if (sign_w .gt. zero) then

            !Carry out U,D,V decomposition.
            do nf = 1, N_FL
               call Phi_0(nf, i_wlk_eff)%decompose
            end do

            Det_D = cmplx(0.d0, 0.d0, kind(0.d0))

            do nf = 1, N_FL
               N_size = phi_0(nf, 1)%n_part
               do I = 1, N_size
                  det_D(nf) = det_D(nf) + log(Phi_0(nf, i_wlk_eff)%D(I))
               end do
            end do

            do nf = 1, N_FL
               Phi_0(nf, i_wlk_eff)%D(:) = cmplx(1.d0, 0.d0, kind(0.d0))
            end do

                 !!!update the overlap when normal propagation
            !if (cop == 'U') then
            !    i_slat = (i_wlk_eff-1)*N_slat
            !    do nrs = 1, N_slat
            !       i_slat = i_slat + 1
            !       overlap(i_slat)=overlap(i_slat)-sum(Det_D)
            !       re_o_am = dble( overlap(i_slat) )
            !       re_o_ph = mod(aimag( overlap(i_slat) ), 2*pi)
            !       overlap(i_slat) = dcmplx(re_o_am, re_o_ph)
            !    enddo
            !endif

         end if

      end do

   end subroutine re_orthonormalize_walkers
   
   subroutine population_control(phi_0, phi_bp_r)
#ifdef MPI
      use mpi
#endif
      use Random_wrap

      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0, phi_bp_r

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, it_wlk, n_exc, pop_exc(N_wlk_mpi, 4)
      integer :: j, it, i_t, i_st, i_ed, nu_wlk, i_src, i_wlk, j_src, j_wlk, n1, n2, nrg, nfrg, ilabel, ncslat
      integer :: i_llim, i_rlim, j_llim, j_rlim, n_part, ii
      real(Kind=kind(0.d0)) :: zero = 1.0e-12, d_scal, sum_w, w_count, re_weight(N_wlk_mpi), ang_w, re_lw
      real(Kind=kind(0.d0)) :: exp_o_abs(n_slat), exp_o_phase(n_slat), dz2, max_re_w
      complex(Kind=kind(0.d0)), allocatable :: w_arr(:), weight_mpi(:)
      complex(Kind=kind(0.d0)) :: overlap_tmp(N_grc)
      complex(Kind=kind(0.d0)) :: Z1, Z2, Z3, Z_s_array(N_slat), Z_r_array(N_slat), zp
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
              ang_w = aimag(weight_k(i_wlk)) 
              re_lw = dble (weight_k(i_wlk))
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

                  i_llim = 1 + (i_wlk - 1)*N_slat; i_rlim = i_wlk*N_slat
                  j_llim = 1 + (j_wlk - 1)*N_slat; j_rlim = j_wlk*N_slat
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
            i_llim = 1 + (i_wlk - 1)*N_slat; i_rlim = i_wlk*N_slat
            Z_s_array(:) = overlap(i_llim:i_rlim)

            ilabel = (it - 1)*nrg
            call mpi_send(Z_s_array, N_slat, MPI_COMPLEX16, j_src, ilabel, Group_comm, IERR)

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
            call mpi_recv(Z_r_array, N_slat, MPI_COMPLEX16, i_src, ilabel, Group_comm, STATUS, IERR)

            j_llim = 1 + (j_wlk - 1)*N_slat; j_rlim = j_wlk*N_slat
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

   subroutine backpropagation(GR_mix, phi_bp_l, phi_bp_r, udvst, Stab_nt, ltau)
#ifdef MPI
      use mpi
#endif
      implicit none

      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(INOUT) :: GR_mix
      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_bp_l, phi_bp_r
      class(UDV_State), dimension(:, :, :), allocatable, intent(INOUT) :: udvst
      integer, dimension(:), allocatable, intent(IN) :: Stab_nt
      integer, intent(IN) :: ltau

      !Local
      complex(Kind=kind(0.d0)) :: GR_bp(NDIM, NDIM, N_FL, N_grc)
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, ltrot_bp, N_op, nstm, nst, ntau
      integer :: i_grc, i_st, i_ed, ns
      complex(Kind=kind(0.d0)) :: z, z_weight, z_sum_overlap, exp_overlap(N_slat)
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

      N_op = size(OP_V, 1)
      ltrot_bp = size(nsigma_bp(1)%f, 2)
      nstm = size(udvst, 1)

          !! initialization
      do i_grc = 1, N_grc
         i_wlk = (i_grc - 1)/N_slat + 1     ! index for walker
         ns = i_grc - (i_wlk - 1)*N_slat ! index for slater det
         do nf = 1, N_FL
            call phi_bp_l(nf, i_grc)%reset('l', wf_l(nf, ns)%P)
            call udvst(nstm, nf, i_grc)%reset('l', wf_l(nf, ns)%P)
         end do
      end do

          !! backpropagation
      nst = nstm - 1
      do nt = ltrot_bp, 1, -1
         ntau = nt - 1

         do i_wlk = 1, N_wlk

            sign_w = cos(aimag(weight_k(i_wlk)))
            if ( sign_w .gt. zero ) then
             
                do ns = 1, N_slat
                   i_grc = ns + (i_wlk - 1)*N_slat
                   do nf = 1, N_FL

                      call Hop_mod_mmthlc_1D2(phi_bp_l(nf, i_grc)%U, nf, 1)

                      do n = N_op, 1, -1
                         call Op_mmultR(phi_bp_l(nf, i_grc)%U, Op_V(n, nf), nsigma_bp(i_wlk)%f(n, nt), 'c', 1)
                      end do

                      call Hop_mod_mmthlc_1D2(phi_bp_l(nf, i_grc)%U, nf, 1)

                   end do
                end do

            end if

         end do

         if (ntau .eq. stab_nt(nst) .and. ntau .ne. 0) then
            call re_orthonormalize_walkers(phi_bp_l, 'N')
            do i_grc = 1, N_grc
            do nf = 1, N_FL
               udvst(nst, nf, i_grc) = phi_bp_l(nf, i_grc)
            end do
            end do
            nst = nst - 1
         end if

      end do

          !! svd at tau = 0
      call re_orthonormalize_walkers(phi_bp_l, 'N')

          !! compute the total weight
      call ham%sum_weight(z_weight)

          !! equal time measurement
      if ( irank .eq. 0 ) call ham%count_obs

      do i_wlk = 1, N_wlk
         
         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then 
         
             ang_w = aimag(weight_k(i_wlk)) 
             re_lw = dble (weight_k(i_wlk))
             re_we = exp(re_lw)*cos(ang_w)
             
             i_st = 1 + (i_wlk - 1)*N_slat; i_ed = i_wlk*N_slat
             exp_overlap(:) = exp(overlap(i_st:i_ed))
             z_sum_overlap = sum(exp_overlap(:))
             do ns = 1, N_slat
                i_grc = ns + (i_wlk - 1)*N_slat
                do nf = 1, N_Fl
                   call CGRP(Z, GR_bp(:, :, nf, i_grc), phi_bp_r(nf, i_wlk), phi_bp_l(nf, i_grc))
                end do
                call ham%obser(gr_bp(:, :, :, i_grc), gr_mix(:, :, :, i_grc), i_grc, re_we, z_weight, z_sum_overlap)
             end do
     
         endif

      end do

           !! time dependence measurement if not call metropolis
      if (ltau .eq. 1) call bp_measure_tau(phi_bp_l, phi_bp_r, udvst, stab_nt)

   end subroutine backpropagation

   subroutine bp_measure_tau(phi_bp_l, phi_bp_r, udvst, stab_nt)
#ifdef MPI
      use mpi
#endif
      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_bp_l, phi_bp_r
      class(UDV_State), dimension(:, :, :), allocatable, intent(IN) :: udvst
      integer, dimension(:), allocatable, intent(IN) :: Stab_nt

      !Local
      complex(Kind=kind(0.d0)) :: Temp(NDIM, NDIM), GR(NDIM, NDIM, N_FL), GRC(NDIM, NDIM, N_FL)
      complex(Kind=kind(0.d0)) :: GT0(NDIM, NDIM, N_FL, N_grc), G00(NDIM, NDIM, N_FL, N_grc)
      complex(Kind=kind(0.d0)) :: GTT(NDIM, NDIM, N_FL, N_grc), G0T(NDIM, NDIM, N_FL, N_grc)
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
      integer :: ns, i_grc, i_st, i_ed
      complex(Kind=kind(0.d0)) :: Z, Z_weight, DETZ, z_sum_overlap, exp_overlap(N_slat)
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(Kind=kind(0.d0)) :: zero = 1.0e-12, sign_w, ang_w, re_lw, re_we

#ifdef MPI
      integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

          !! initialization
      N_op = size(OP_V, 1)
      nstm = size(udvst, 1)

          !! compute the total weight
      call ham%sum_weight(z_weight)

      do i_grc = 1, N_grc
         i_wlk = (i_grc - 1)/N_slat + 1     ! index for walker
         ns = i_grc - (i_wlk - 1)*N_slat ! index for slater det

         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then
            do nf = 1, N_FL
               call CGRP(DetZ, gr(:, :, nf), phi_bp_r(nf, i_wlk), phi_bp_l(nf, i_grc))
            end do
            gtt(:, :, :, i_grc) = gr
         endif
      end do

      G00 = GTT
      GT0 = GTT
      G0T = GTT
      do i_grc = 1, N_grc
      do nf = 1, N_FL
         do I = 1, Ndim
            G0T(I, I, nf, i_grc) = G0T(I, I, nf, i_grc) - 1.d0
         end do
      end do
      end do

      ntau = 0
      do i_wlk = 1, N_wlk
         
         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then
         
             ang_w = aimag(weight_k(i_wlk)) 
             re_lw = dble (weight_k(i_wlk))
             re_we = exp(re_lw)*cos(ang_w)
             
             i_st = 1 + (i_wlk - 1)*N_slat; i_ed = i_wlk*N_slat
             exp_overlap(:) = exp(overlap(i_st:i_ed))
             z_sum_overlap = sum(exp_overlap(:))
             do ns = 1, N_slat
                i_grc = ns + (i_wlk - 1)*N_slat
                call ham%obsert(ntau, gt0(:, :, :, i_grc), g0t(:, :, :, i_grc), g00(:, :, :, i_grc), &
                    & gtt(:, :, :, i_grc), i_grc, re_we, z_weight, z_sum_overlap)
             end do

         endif

      end do

      NST = 1
      do ntau = 1, ltrot

         do i_wlk = 1, N_wlk
                !! Propagate wave function
            sign_w = cos(aimag(weight_k(i_wlk)))
            if ( sign_w .gt. zero ) then

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat

                   !! Propagate Green's function
                  call PROPR(GT0(:, :, :, i_grc), ntau, i_wlk)
                  call PROPRM1(G0T(:, :, :, i_grc), ntau, i_wlk)
                  call PROPRM1(GTT(:, :, :, i_grc), ntau, i_wlk)
                  call PROPR(GTT(:, :, :, i_grc), ntau, i_wlk)

               end do

               do nf = 1, N_FL

                  call Hop_mod_mmthr_1D2(phi_bp_r(nf, i_wlk)%U, nf, 1)

                  do n = 1, N_op
                     call Op_mmultR(phi_bp_r(nf, i_wlk)%U, op_v(n, nf), nsigma_bp(i_wlk)%f(n, ntau), 'n', 1)
                  end do

                  call Hop_mod_mmthr_1D2(phi_bp_r(nf, i_wlk)%U, nf, 1)

               end do

            end if

            !! measure
            if ( sign_w .gt. zero ) then

                 ang_w = aimag(weight_k(i_wlk)) 
                 re_lw = dble (weight_k(i_wlk))
                 re_we = exp(re_lw)*cos(ang_w)

                 i_st = 1 + (i_wlk - 1)*N_slat; i_ed = i_wlk*N_slat
                 exp_overlap(:) = exp(overlap(i_st:i_ed))
                 z_sum_overlap = sum(exp_overlap(:))
                 do ns = 1, N_slat
                    i_grc = ns + (i_wlk - 1)*N_slat
                    call ham%obsert(ntau, gt0(:, :, :, i_grc), g0t(:, :, :, i_grc), g00(:, :, :, i_grc), &
                        & gtt(:, :, :, i_grc), i_grc, re_we, z_weight, z_sum_overlap)
                 end do

            endif

         end do

             !! call svd
         if (ntau .eq. stab_nt(nst)) then
            call re_orthonormalize_walkers(phi_bp_r, 'N')

            do i_wlk = 1, N_wlk
            
               sign_w = cos(aimag(weight_k(i_wlk)))
               if ( sign_w .gt. zero ) then

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf = 1, N_FL
                     phi_bp_l(nf, i_grc) = udvst(nst, nf, i_grc)
                     call CGRP(DetZ, gr(:, :, nf), phi_bp_r(nf, i_wlk), phi_bp_l(nf, i_grc))
                     call Control_Precision_tau(gtt(:, :, nf, i_grc), GR(:, :, nf), Ndim)
                  end do
                  GTT(:, :, :, i_grc) = GR

                  GRC = -GR
                  do nf = 1, N_FL
                     do I = 1, Ndim
                        GRC(I, I, nf) = GRC(I, I, nf) + 1.d0
                     end do
                  end do

                  do nf = 1, N_FL
                     call MMULT(TEMP, GR(:, :, nf), GT0(:, :, nf, i_grc))
                     GT0(:, :, nf, i_grc) = TEMP
                     call MMULT(TEMP, G0T(:, :, nf, i_grc), GRC(:, :, nf))
                     G0T(:, :, nf, i_grc) = TEMP
                  end do

               end do

               endif

            end do
            nst = nst + 1
         end if

      end do

   end subroutine bp_measure_tau

   subroutine rescale_overlap(overlap_in)

      implicit none

      complex(Kind=kind(0.d0)), dimension(:), allocatable, intent(inout) :: overlap_in

      integer :: nf, n, m, nt, i_wlk, i_grc, ns
      integer :: i_st, i_ed, ncslat
      real(Kind=kind(0.d0)) :: log_o_abs(n_slat), log_o_phase(n_slat), dz2
      real(kind=kind(0.d0)) :: pi = acos(-1.d0), dre_o, zero = 1.0e-12, sign_w
      complex(Kind=kind(0.d0)) :: z1, zp

      do i_wlk = 1, N_wlk

         sign_w = cos(aimag(weight_k(i_wlk)))
         if ( sign_w .gt. zero ) then

            i_st = 1 + (i_wlk - 1)*N_slat
            i_ed = i_wlk*N_slat

            ncslat = 0
            do i_grc = i_st, i_ed
               ncslat = ncslat + 1
               log_o_abs(ncslat) = dble(overlap_in(i_grc))
               log_o_phase(ncslat) = mod(aimag(overlap_in(i_grc)), 2.d0*pi)
            end do

            dz2 = maxval(log_o_abs(:))

            ncslat = 0
            do i_grc = i_st, i_ed
               ncslat = ncslat + 1
               dre_o = log_o_abs(ncslat) - dz2
               overlap_in(i_grc) = cmplx(dre_o, log_o_phase(ncslat), kind(0.d0))
            end do

         end if
      end do

   end subroutine rescale_overlap

   subroutine PROPR(AIN, NT, i_wlk)

      ! Ain =       B(NT-1, NT1)
      ! Aout= Ain = B(NT  , NT1)

      implicit none
      complex(Kind=kind(0.d0)), intent(INOUT) :: Ain(Ndim, Ndim, N_FL)
      integer, intent(IN) :: NT, i_wlk

      !Locals
      integer :: nf, n

      do nf = 1, N_FL
         call Hop_mod_mmthr_1D2(Ain(:, :, nf), nf, nt)
         do n = 1, size(Op_V, 1)
            call Op_mmultR(Ain(:, :, nf), Op_V(n, nf), nsigma_bp(i_wlk)%f(n, nt), 'n', nt)
         end do
         call Hop_mod_mmthr_1D2(Ain(:, :, nf), nf, nt)
      end do

   end subroutine PROPR

   subroutine PROPRM1(AIN, NT, i_wlk)

      !Ain = B^{-1}(NT-1, NT1)
      !Aout= B^{-1}(NT  , NT1)

      implicit none

      !Arguments
      complex(Kind=kind(0.d0)), intent(Inout) ::  AIN(Ndim, Ndim, N_FL)
      integer :: NT, i_wlk

      ! Locals
      integer :: nf, n

      do nf = 1, N_FL
         call Hop_mod_mmthl_m1_1D2(Ain(:, :, nf), nf, nt)
         do n = 1, size(Op_V, 1)
            call Op_mmultL(Ain(:, :, nf), Op_V(n, nf), -nsigma_bp(i_wlk)%f(n, nt), 'n', nt)
         end do
         call Hop_mod_mmthl_m1_1D2(Ain(:, :, nf), nf, nt)
      end do

   end subroutine PROPRM1

end module stepwlk_mod
