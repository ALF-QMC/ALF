module stepwlk_mod
   use Hamiltonian_main
   use Operator_mod
   use Control
   use Hop_mod
   use UDV_State_mod
   use gfun_mod
   use upgrade_mod

contains

   subroutine initial_wlk(phi_trial, phi_0, phi_bp_l, phi_bp_r, udvst, stab_nt, GR, nwrap)
#ifdef MPI
      use mpi
#endif
      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(INOUT) :: phi_trial, phi_0, phi_bp_l, phi_bp_r
      class(udv_state), dimension(:, :, :), allocatable, intent(INOUT) :: udvst
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(INOUT) :: GR
      integer, dimension(:), allocatable, intent(INOUT) :: stab_nt
      integer, intent(IN) :: nwrap

      !Local
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_grc, NSTM, NST, ltrot_bp, ns
      integer :: i_st, i_ed, ncslat
      complex(Kind=kind(0.d0)) :: overlap_old, overlap_new, Z, Z1, Z2, tot_ene, ZP
      complex(Kind=kind(0.d0)) :: tot_c_weight, el_tmp
      complex(Kind=kind(0.d0)) :: det_Vec(N_FL)
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio, X1, wtmp
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, tot_re_weight, dz2
      character(LEN=64) :: FILE_TG, FILE_seeds, file_inst, file_antiinst
      logical ::   LCONF, LCONF_H5, lconf_inst, lconf_antiinst

#ifdef MPI
      integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

      nstm = size(udvst, 1)
      ltrot_bp = size(nsigma_bp(1)%f, 2)
      tot_ene = cmplx(0.d0, 0.d0, kind(0.d0))
      tot_re_weight = 0.d0

      allocate (Stab_nt(0:nstm))
      Stab_nt(0) = 0
      do n = 1, Nstm - 1
         Stab_nt(n) = nwrap*n
      end do
      Stab_nt(Nstm) = ltrot_bp

      do i_wlk = 1, N_wlk
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call phi_0(nf_eff, i_wlk)%init(ndim, 'r', WF_R(nf, 1)%P)
            call phi_bp_r(nf_eff, i_wlk)%init(ndim, 'r', WF_R(nf, 1)%P)
         end do
      end do

      do ns = 1, N_slat
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call phi_trial(nf_eff, ns)%init(ndim, 'l', WF_L(nf, ns)%P)
         end do
         do i_wlk = 1, N_wlk
            i_grc = ns + (i_wlk - 1)*N_slat
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               do n = 1, nstm
                  call udvst(n, nf_eff, i_grc)%alloc(ndim)
               end do
               call phi_bp_l(nf_eff, i_grc)%init(ndim, 'l', WF_L(nf, ns)%P)
            end do
         end do
      end do

      file_inst = 'trial_inst.h5'
      file_antiinst = 'trial_antiinst.h5'
      inquire (file=file_inst, exist=lconf_inst)
      inquire (file=file_antiinst, exist=lconf_antiinst)
      if (lconf_inst .and. lconf_antiinst) then
         if (irank_g .eq. 0) write (*, *) "read inst-anti inst config as trial wave function"
         call trial_in_hdf5(phi_0, phi_trial, file_inst, file_antiinst)
      end if

      file_tg = "phiin_0.h5"
      inquire (FILE=file_tg, EXIST=LCONF_H5)
      if (LCONF_H5) then
         if (irank_g .eq. 0) write (*, *) "read input walkers"
         call wavefunction_in_hdf5(phi_0, file_tg)
      end if

          !! initial overlap and green's function
      do i_wlk = 1, N_wlk

         if (weight_k(i_wlk) .le. 0.d0) weight_k(i_wlk) = 0.d0

         do ns = 1, N_slat
            i_grc = ns + (i_wlk - 1)*N_slat

            do nf_eff = 1, N_Fl_eff
               nf = Calc_Fl_map(nf_eff)
               call cgrp(Z, GR(:, :, nf, i_grc), phi_0(nf_eff, i_wlk), phi_trial(nf_eff, ns))
               det_vec(nf_eff) = Z
            end do
            det_Vec(:) = det_Vec(:)*N_SUN
            if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
            if (.not. LCONF_H5) overlap(i_grc) = sum(det_vec)
         end do
      end do

          !! rescale overlap
      call rescale_overlap(overlap)

          !! initial energy
      call ham%update_fac_norm(GR, 0)

      file_seeds = "seedvec_in"
      inquire (FILE=file_seeds, EXIST=LCONF)
      if (LCONF) then
         if (irank_g .eq. 0) write (*, *) "read input seeds"
         call seed_vec_in(file_seeds)
      end if

   end subroutine initial_wlk

   subroutine stepwlk_move(phi_trial, phi_0, GR, ntau_bp)

      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(IN)    :: phi_trial
      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(INOUT) :: GR
      integer, intent(IN) :: ntau_bp

      !Local
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8

      N_op = size(OP_V, 1)

      do i_wlk = 1, N_wlk
         ! update weight by fac_norm
         if (weight_k(i_wlk) .gt. Zero) weight_k(i_wlk) = weight_k(i_wlk)*exp(fac_norm)
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
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed
      integer :: n_op, ns, i_grc
      complex(Kind=kind(0.d0)) :: Z, z_alpha
      complex(Kind=kind(0.d0)) :: det_Vec(N_FL)
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, spin

      N_op = size(OP_V, 1)

      do i_wlk = 1, N_wlk

         if (weight_k(i_wlk) .gt. Zero) then

            ! upgrade Green's function
            z_alpha = 0.d0
            do n = 1, N_op

               N_type = 2
               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     call Op_Wrapdo(GR(:, :, nf, i_grc), Op_V(n, nf), 1.d0, Ndim, N_Type, 1)
                  end do
               end do

               call Upgrade(GR, n, spin, i_wlk)
               nsigma_bp(i_wlk)%f(n, ntau_bp) = spin

               N_type = 2
               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     call Op_Wrapup(Gr(:, :, nf, i_grc), Op_V(n, nf), 1.d0, Ndim, N_Type, 1)
                  end do
               end do

               ! propagate slater determinant
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  call Op_mmultR(phi_0(nf_eff, i_wlk)%U, Op_V(n, nf), spin, 'n', 1)
                  ! store alpha factor in z_alpha
                  z_alpha = z_alpha + &
                      & op_v(n_op, nf)%g*op_v(n_op, nf)%alpha*nsigma_bp(i_wlk)%Phi(n, ntau_bp)
               end do

            end do

            !! comment this since it would cause instability when z_alpha is real 
            !! rescale U matrix with z_alpha factor
            !do nf_eff = 1, N_FL_eff
            !   nf = Calc_Fl_map(nf_eff)
            !   phi_0(nf_eff, i_wlk)%U(:, :) = exp(z_alpha)*phi_0(nf_eff, i_wlk)%U(:, :)
            !end do

         end if

      end do

   end subroutine int_propagation

   subroutine half_K_propagation(phi_trial, phi_0, GR)

      implicit none

      class(udv_state), intent(in), allocatable :: phi_trial(:, :)
      class(udv_state), intent(inout), allocatable :: phi_0(:, :)
      complex(Kind=kind(0.d0)), intent(inout), allocatable :: GR(:, :, :, :)

      ! local
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed, ns, i_grc
      complex(Kind=kind(0.d0)) :: Overlap_old, Overlap_new, Z, sum_o_new, sum_o_old
      complex(Kind=kind(0.d0)) :: det_Vec(N_FL), log_o_new(n_slat), log_o_old(n_slat)
      real(Kind=kind(0.d0)) :: overlap_ratio, re_overlap
      real(Kind=kind(0.d0)) :: zero = 1.0e-8

      do i_wlk = 1, N_wlk

         if (weight_k(i_wlk) .gt. zero) then

            sum_o_old = 0.d0
            do ns = 1, N_slat
               i_grc = ns + (i_wlk - 1)*N_slat
               ! Update Green's function
               do nf_eff = 1, N_Fl_eff
                  nf = Calc_Fl_map(nf_eff)
                  call CGRP(Z, GR(:, :, nf, i_grc), phi_0(nf_eff, i_wlk), phi_trial(nf_eff, ns))
                  det_vec(nf) = Z
               end do
               det_Vec(:) = det_Vec(:)*N_SUN
               if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
               log_o_old(ns) = sum(det_Vec)
               sum_o_old = sum_o_old + exp(log_o_old(ns))
            end do

             !! multi exp(-\Delta\tau T/2)
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call Hop_mod_mmthr_1D2(phi_0(nf_eff, i_wlk)%U, nf, 1)
            end do

            sum_o_new = 0.d0
            do ns = 1, N_slat
               i_grc = ns + (i_wlk - 1)*N_slat
               ! Update Green's function
               do nf_eff = 1, N_Fl_eff
                  nf = Calc_Fl_map(nf_eff)
                  call CGRP(Z, GR(:, :, nf, i_grc), phi_0(nf_eff, i_wlk), phi_trial(nf_eff, ns))
                  det_vec(nf) = Z
               end do
               det_Vec(:) = det_Vec(:)*N_SUN
               if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
               log_o_new(ns) = sum(det_Vec)
               sum_o_new = sum_o_new + exp(log_o_new(ns))
            end do

            overlap_ratio = sum_o_new/sum_o_old
            re_overlap = dble(overlap_ratio)
            if (re_overlap .gt. zero) then
                 !! upgrade overlap
               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  overlap(i_grc) = overlap(i_grc) + &
                      & (log_o_new(ns) - log_o_old(ns))
               end do
               weight_k(i_wlk) = weight_k(i_wlk)*re_overlap
            else
               weight_k(i_wlk) = 0.d0
            end if

         end if

      end do

   end subroutine half_K_propagation

   subroutine re_orthonormalize_walkers(Phi_0, cop)

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0
      character, intent(IN)    :: cop

      !Local
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, N_size, I, i_wlk
      integer :: ndistance, i_wlk_eff, i_st, i_ed, n_wlk_eff, ns, nrs, i_slat
      real(Kind=kind(0.d0)) :: Overlap_ratio, Zero = 1.0e-8, pi = acos(-1.d0)
      real(Kind=kind(0.d0)) :: re_o_am, re_o_ph
      complex(Kind=kind(0.d0)) :: det_D(N_FL)

      n_wlk_eff = size(phi_0, 2)
      ndistance = n_wlk_eff/n_wlk

      do i_wlk_eff = 1, n_wlk_eff

         i_wlk = (i_wlk_eff - 1)/ndistance + 1     ! index for walker
         ns = i_wlk_eff - (i_wlk - 1)*ndistance ! index for slater det

         if (weight_k(i_wlk) .gt. Zero) then

            N_size = phi_0(1, 1)%n_part

            !Carry out U,D,V decomposition.
            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               call Phi_0(nf_eff, i_wlk_eff)%decompose
            end do

            Det_D = cmplx(0.d0, 0.d0, kind(0.d0))

            do nf_eff = 1, N_FL_eff
               nf = Calc_Fl_map(nf_eff)
               do I = 1, N_size
                  det_D(nf) = det_D(nf) + log(Phi_0(nf_eff, i_wlk_eff)%D(I))
               end do
            end do

            if (reconstruction_needed) call ham%weight_reconstruction(Det_D)

            do nf_eff = 1, N_FL_eff
               Phi_0(nf_eff, i_wlk_eff)%D(:) = cmplx(1.d0, 0.d0, kind(0.d0))
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
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, it_wlk, n_exc, pop_exc(N_wlk_mpi, 4)
      integer :: j, it, i_t, i_st, i_ed, nu_wlk, i_src, i_wlk, j_src, j_wlk, n1, n2, nrg, nfrg, ilabel, ncslat
      integer :: i_llim, i_rlim, j_llim, j_rlim
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8, d_scal, sum_w, w_count, w_tmp(N_wlk_mpi), weight_mpi(N_wlk_mpi)
      real(Kind=kind(0.d0)) :: exp_o_abs(n_slat), exp_o_phase(n_slat), dz2
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
      allocate (phi_0_m(N_FL_eff, N_wlk))
      allocate (phi_bp_m(N_FL_eff, N_wlk))
      do i_wlk = 1, N_wlk
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call phi_0_m(nf_eff, i_wlk)%alloc(ndim, phi_0(1, 1)%n_part)
            call phi_bp_m(nf_eff, i_wlk)%alloc(ndim, phi_0(1, 1)%n_part)
         end do
      end do

      do nf_eff = 1, N_FL_eff
         do i_wlk = 1, N_wlk
            phi_0_m(nf_eff, i_wlk) = phi_0(nf_eff, i_wlk)
            phi_bp_m(nf_eff, i_wlk) = phi_bp_r(nf_eff, i_wlk)
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

      ! population control
      weight_mpi(:) = 0.d0
      i_st = irank_g*N_wlk + 1
      i_ed = (irank_g + 1)*N_wlk

      weight_mpi(i_st:i_ed) = weight_k(:)

      call MPI_REDUCE(weight_mpi, w_tmp, N_wlk_mpi, MPI_REAL8, MPI_SUM, 0, Group_comm, IERR)
      if (irank_g == 0) weight_mpi = w_tmp
      call MPI_BCAST(weight_mpi, N_wlk_mpi, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      d_scal = dble(N_wlk_mpi)/sum(weight_mpi)
      if (irank_g == 0) sum_w = -ranf_wrap()
      call MPI_BCAST(sum_w, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      nu_wlk = 0

      n_exc = 0
      do it_wlk = 1, N_wlk_mpi
         i_src = (it_wlk - 1)/N_wlk
         i_wlk = it_wlk - N_wlk*i_src

         sum_w = sum_w + weight_mpi(it_wlk)*d_scal
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
                  do nf_eff = 1, N_FL_eff
                     phi_0_m(nf_eff, j_wlk) = phi_0(nf_eff, i_wlk)
                     phi_bp_m(nf_eff, j_wlk) = phi_bp_r(nf_eff, i_wlk)
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

      nrg = 2 + N_fl_eff*6

      do it = 1, n_exc
         i_src = pop_exc(it, 1); i_wlk = pop_exc(it, 2)
         j_src = pop_exc(it, 3); j_wlk = pop_exc(it, 4)
         if (irank_g .eq. i_src) then
            i_llim = 1 + (i_wlk - 1)*N_slat; i_rlim = i_wlk*N_slat
            Z_s_array(:) = overlap(i_llim:i_rlim)

            ilabel = (it - 1)*nrg
            call mpi_send(Z_s_array, N_slat, MPI_COMPLEX16, j_src, ilabel, Group_comm, IERR)

            do nf_eff = 1, N_FL_eff
               ilabel = (it - 1)*nrg + (nf_eff - 1)*6 + 1
               call phi_0(nf_eff, i_wlk)%MPI_send_general(j_src, ilabel, ierr)
               ilabel = (it - 1)*nrg + (nf_eff - 1)*6 + 4
               call phi_bp_r(nf_eff, i_wlk)%MPI_send_general(j_src, ilabel, ierr)
            end do

            ilabel = (it - 1)*nrg + (n_fl_eff)*6 + 1
            call mpi_send(nsigma_bp(i_wlk)%f, n1*n2, MPI_REAL8, j_src, ilabel, Group_comm, IERR)
         end if
         if (irank_g .eq. j_src) then
            ilabel = (it - 1)*nrg
            call mpi_recv(Z_r_array, N_slat, MPI_COMPLEX16, i_src, ilabel, Group_comm, STATUS, IERR)

            j_llim = 1 + (j_wlk - 1)*N_slat; j_rlim = j_wlk*N_slat
            overlap_tmp(j_llim:j_rlim) = Z_r_array(:)

            do nf_eff = 1, N_FL_eff
               ilabel = (it - 1)*nrg + (nf_eff - 1)*6 + 1
               call phi_0_m(nf_eff, j_wlk)%MPI_recv_general(i_src, ilabel, status, ierr)
               ilabel = (it - 1)*nrg + (nf_eff - 1)*6 + 4
               call phi_bp_m(nf_eff, j_wlk)%MPI_recv_general(i_src, ilabel, status, ierr)
            end do

            ilabel = (it - 1)*nrg + (n_fl_eff)*6 + 1
            call mpi_recv(nsigma_store(j_wlk)%f, n1*n2, MPI_REAL8, i_src, ilabel, Group_comm, STATUS, IERR)
         end if
      end do

      overlap(:) = overlap_tmp(:)
      ! reset weight
      weight_k(:) = 1.d0
      do nf_eff = 1, N_FL_eff
         do i_wlk = 1, N_wlk
            phi_0(nf_eff, i_wlk) = phi_0_m(nf_eff, i_wlk)
            phi_bp_r(nf_eff, i_wlk) = phi_bp_m(nf_eff, i_wlk)
         end do
      end do

          !!! rescale overlap
      !call rescale_overlap(overlap)

      do i_wlk = 1, N_wlk
         nsigma_bp(i_wlk)%f = nsigma_store(i_wlk)%f
         nsigma_bp(i_wlk)%t = nsigma_store(i_wlk)%t
      end do

      ! deallocate tmp udv class
      do i_wlk = 1, N_wlk
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call phi_0_m(nf_eff, i_wlk)%dealloc
            call phi_bp_m(nf_eff, i_wlk)%dealloc
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
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk

      do nf_eff = 1, N_FL_eff
         do i_wlk = 1, N_wlk
            phi_bp_r(nf_eff, i_wlk) = phi_0(nf_eff, i_wlk)
                !! phi_bp_r is not used in propagation
            call phi_bp_r(nf_eff, i_wlk)%decompose
         end do
      end do

   end subroutine store_phi

   subroutine backpropagation(GR_mix, phi_bp_l, phi_bp_r, udvst, Stab_nt, ltau, lmetropolis, nsweep, nwarmup)
#ifdef MPI
      use mpi
#endif
      implicit none

      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(INOUT) :: GR_mix
      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_bp_l, phi_bp_r
      class(UDV_State), dimension(:, :, :), allocatable, intent(INOUT) :: udvst
      integer, dimension(:), allocatable, intent(IN) :: Stab_nt
      integer, intent(IN) :: ltau, lmetropolis, nsweep, nwarmup

      !Local
      complex(Kind=kind(0.d0)) :: GR_bp(NDIM, NDIM, N_FL, N_grc)
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, ltrot_bp, N_op, nstm, nst, ntau
      integer :: i_grc, i_st, i_ed, act_mea, ns
      complex(Kind=kind(0.d0)) :: z, z_weight, z_sum_overlap, exp_overlap(N_slat)
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(Kind=kind(0.d0))    :: Zero = 1.0e-8

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
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call phi_bp_l(nf_eff, i_grc)%reset('l', wf_l(nf, ns)%P)
            call udvst(nstm, nf_eff, i_grc)%reset('l', wf_l(nf, ns)%P)
         end do
      end do

          !! backpropagation
      nst = nstm - 1
      do nt = ltrot_bp, 1, -1
         ntau = nt - 1

         do i_wlk = 1, N_wlk

            if (weight_k(i_wlk) .gt. Zero) then

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf_eff = 1, N_FL_eff

                     nf = Calc_Fl_map(nf_eff)
                     call Hop_mod_mmthlc_1D2(phi_bp_l(nf_eff, i_grc)%U, nf, 1)

                     do n = N_op, 1, -1
                        call Op_mmultR(phi_bp_l(nf_eff, i_grc)%U, Op_V(n, nf), nsigma_bp(i_wlk)%f(n, nt), 'c', 1)
                     end do

                     call Hop_mod_mmthlc_1D2(phi_bp_l(nf_eff, i_grc)%U, nf, 1)

                  end do
               end do

            end if

         end do

         if (ntau .eq. stab_nt(nst) .and. ntau .ne. 0) then
            call re_orthonormalize_walkers(phi_bp_l, 'N')
            do i_grc = 1, N_grc
            do nf_eff = 1, N_FL_eff
               udvst(nst, nf_eff, i_grc) = phi_bp_l(nf_eff, i_grc)
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
      act_mea = 0 + irank
      do i_wlk = 1, N_wlk
         i_st = 1 + (i_wlk - 1)*N_slat; i_ed = i_wlk*N_slat
         exp_overlap(:) = exp(overlap(i_st:i_ed))
         z_sum_overlap = sum(exp_overlap(:))
         do ns = 1, N_slat
            i_grc = ns + (i_wlk - 1)*N_slat
            do nf_eff = 1, N_Fl_eff
               nf = Calc_Fl_map(nf_eff)
               call CGRP(Z, GR_bp(:, :, nf, i_grc), phi_bp_r(nf_eff, i_wlk), phi_bp_l(nf_eff, i_grc))
            end do
            if (reconstruction_needed) call ham%GR_reconstruction(GR_bp(:, :, :, i_grc))
            if (reconstruction_needed) call ham%GR_reconstruction(GR_mix(:, :, :, i_grc))
            call ham%Obser(GR_bp(:, :, :, i_grc), GR_mix(:, :, :, i_grc), i_wlk, i_grc, z_weight, z_sum_overlap, act_mea)
            act_mea = act_mea + 1
         end do
      end do

      if (lmetropolis .eq. 1) then
         ! metripolis sampling
         call metropolis(phi_bp_r, phi_bp_l, udvst, stab_nt, nsweep, nwarmup, ltau)

         if (ltau .eq. 1) then
            act_mea = 0 + irank
            do i_wlk = 1, N_wlk
               call ham%bp_obsert(i_wlk, i_grc, z_weight, act_mea)
               act_mea = act_mea + 1
            end do
         end if
      else
              !! time dependence measurement if not call metropolis
         if (ltau .eq. 1) call bp_measure_tau(phi_bp_l, phi_bp_r, udvst, stab_nt)
      end if

   end subroutine backpropagation

   subroutine metropolis(phi_bp_r, phi_l_m, udvst, stab_nt, nsweep, nwarmup, ltau)

      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(in)    :: phi_bp_r
      class(udv_state), dimension(:, :), allocatable, intent(inout) :: phi_l_m
      class(udv_state), dimension(:, :, :), allocatable, intent(inout) :: udvst
      integer, dimension(:), allocatable, intent(in) :: stab_nt
      integer, intent(in) :: nsweep, nwarmup, ltau

      !Local
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
      integer :: nsw, ltrot_bp, nmea, ltrot_eff, i_slat, ns, i_grc
      complex(Kind=kind(0.d0)) :: z, z_weight, detz, z1, z2, zp, ztmp, z_avg, z_sgn_avg, ener_tmp
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, overlap_ratio, hs_field
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8
      complex(Kind=kind(0.d0)) :: gr(ndim, ndim, n_fl)
      complex(Kind=kind(0.d0)) :: det_vec(n_fl), zph1, zph2
      class(udv_state), dimension(:, :), allocatable :: phi_r_m
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable :: gtt
      complex(Kind=kind(0.d0)), dimension(:), allocatable :: overlap_mc

#ifdef MPI
      integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
      integer        :: STATUS(MPI_STATUS_SIZE)
      call mpi_comm_size(mpi_comm_world, isize, ierr)
      call mpi_comm_rank(mpi_comm_world, irank, ierr)
      call mpi_comm_rank(group_comm, irank_g, ierr)
      call mpi_comm_size(group_comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

          !! initialization
      n_op = size(op_v, 1)
      nstm = size(udvst, 1)
      ltrot_bp = size(nsigma_bp(1)%f, 2)
      ltrot_eff = ltrot

      allocate (overlap_mc(n_grc))
      overlap_mc(:) = overlap(:)

          !! init observables
      if (ltau .eq. 1) call ham%init_obs_mc

          !! allocate tmp wavefunction
      allocate (phi_r_m(N_FL_eff, N_wlk))
      allocate (gtt(ndim, ndim, n_fl, n_grc))

          !! compute the total weight
      call ham%sum_weight(z_weight)

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do i_wlk = 1, N_wlk
            call phi_r_m(nf_eff, i_wlk)%init(ndim, 'r', phi_bp_r(nf_eff, i_wlk)%U)
         end do
         do i_grc = 1, N_grc
            i_wlk = (i_grc - 1)/N_slat + 1      ! index for walker
            i_slat = i_grc - (i_wlk - 1)*N_slat  ! index for slater det
                !! init Green's function
            call cgrp(z, gtt(:, :, nf, i_grc), phi_r_m(nf_eff, i_wlk), phi_l_m(nf_eff, i_grc))
         end do
      end do

      do nsw = 1, nsweep

         nst = 1
         do ntau = 1, ltrot_bp

            do i_wlk = 1, N_wlk

               if (weight_k(i_wlk) .gt. zero) then

                 !! Propagate Green's function
                  do ns = 1, N_slat
                     i_grc = ns + (i_wlk - 1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        call hop_mod_mmthr_1d2(gtt(:, :, nf, i_grc), nf, 1)
                        call hop_mod_mmthl_m1_1d2(gtt(:, :, nf, i_grc), nf, 1)
                     end do
                  end do

                  do n = 1, n_op

                     n_type = 1
                     do ns = 1, N_slat
                        i_grc = ns + (i_wlk - 1)*N_slat
                        do nf_eff = 1, N_FL_eff
                           nf = Calc_Fl_map(nf_eff)
                           hs_field = nsigma_bp(i_wlk)%f(n, ntau)
                           call op_wrapup(gtt(:, :, nf, i_grc), op_v(n, nf), hs_field, ndim, n_type, ntau)
                        end do
                     end do

                     ! metropolis update
                     hs_new = nsigma_bp(i_wlk)%flip(n, ntau)
                     if (ntau .le. ltrot_eff) then
                        call upgrade_mc(gtt, n, ntau, hs_new, i_wlk, overlap_mc)
                     end if

                     n_type = 2
                     do ns = 1, N_slat
                        i_grc = ns + (i_wlk - 1)*N_slat
                        do nf_eff = 1, N_FL_eff
                           nf = Calc_Fl_map(nf_eff)
                           hs_field = nsigma_bp(i_wlk)%f(n, ntau)
                           call op_wrapup(gtt(:, :, nf, i_grc), op_v(n, nf), hs_field, ndim, n_type, ntau)
                        end do
                     end do

                  end do

                  do ns = 1, N_slat
                     i_grc = ns + (i_wlk - 1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        call hop_mod_mmthr_1d2(gtt(:, :, nf, i_grc), nf, 1)
                        call hop_mod_mmthl_m1_1d2(gtt(:, :, nf, i_grc), nf, 1)
                     end do
                  end do

                 !! Propagate wave function
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     call hop_mod_mmthr_1d2(phi_r_m(nf_eff, i_wlk)%U, nf, 1)
                     do n = 1, N_op
                        call op_mmultr(phi_r_m(nf_eff, i_wlk)%U, op_v(n, nf), nsigma_bp(i_wlk)%f(n, ntau), 'n', 1)
                     end do
                     call hop_mod_mmthr_1d2(phi_r_m(nf_eff, i_wlk)%U, nf, 1)
                  end do

               end if

            end do

             !! call svd
            if (ntau .eq. stab_nt(nst)) then
               call re_orthonormalize_walkers(phi_r_m, 'N')

               do i_wlk = 1, N_wlk

                  if (weight_k(i_wlk) .gt. zero) then

                     do ns = 1, N_slat
                        i_grc = ns + (i_wlk - 1)*N_slat
                        do nf_eff = 1, N_FL_eff
                           nf = Calc_Fl_map(nf_eff)
                           phi_l_m(nf_eff, i_grc) = udvst(nst, nf_eff, i_grc)
                           !! compute Green's function
                           call cgrp(z, gr(:, :, nf), phi_r_m(nf_eff, i_wlk), phi_l_m(nf_eff, i_grc))
                           det_vec(nf_eff) = z
                           call control_precisiong(gr(:, :, nf), gtt(:, :, nf, i_grc), Ndim)
                           gtt(:, :, nf, i_grc) = gr(:, :, nf)
                        end do
                        det_Vec(:) = det_Vec(:)*N_SUN
                        zph1 = exp(dcmplx(0.d0, aimag(sum(det_vec))))
                        zph2 = exp(dcmplx(0.d0, aimag(overlap_mc(i_grc))))
                        if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                        call control_Precisionp(zph1, zph2)
                        overlap_mc(i_grc) = sum(det_vec)
                     end do

                     !! store on the first slat det storage
                     i_grc = 1 + (i_wlk - 1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        udvst(nst, nf_eff, i_grc) = phi_r_m(nf_eff, i_wlk)
                     end do

                  end if

               end do
               call rescale_overlap(overlap_mc)
               nst = nst + 1
            end if

         end do

         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            do i_grc = 1, N_grc
               i_wlk = (i_grc - 1)/N_slat + 1      ! index for walker
               i_slat = i_grc - (i_wlk - 1)*N_slat  ! index for slater det
               call phi_l_m(nf_eff, i_grc)%reset('l', wf_l(nf, i_slat)%P)
            end do
         end do

         nst = nstm - 1
         do ntau = ltrot_bp, 1, -1
            ntau1 = ntau - 1

            do i_wlk = 1, N_wlk

               if (weight_k(i_wlk) .gt. Zero) then

                 !! Propagate Green's function
                  do ns = 1, N_slat
                     i_grc = ns + (i_wlk - 1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        call hop_mod_mmthl_1d2(gtt(:, :, nf, i_grc), nf, 1)
                        call hop_mod_mmthr_m1_1d2(gtt(:, :, nf, i_grc), nf, 1)
                     end do
                  end do

                  do n = n_op, 1, -1

                     n_type = 2
                     do ns = 1, N_slat
                        i_grc = ns + (i_wlk - 1)*N_slat
                        do nf_eff = 1, N_FL_eff
                           nf = Calc_Fl_map(nf_eff)
                           hs_field = nsigma_bp(i_wlk)%f(n, ntau)
                           call op_wrapdo(gtt(:, :, nf, i_grc), op_v(n, nf), hs_field, ndim, N_Type, ntau)
                        end do
                     end do

                     ! metropolis update
                     hs_new = nsigma_bp(i_wlk)%flip(n, ntau)
                     if (ntau .le. ltrot_eff) then
                        call upgrade_mc(gtt, n, ntau, hs_new, i_wlk, overlap_mc)
                     end if

                     n_type = 1
                     hs_field = nsigma_bp(i_wlk)%f(n, ntau)
                     do ns = 1, N_slat
                        i_grc = ns + (i_wlk - 1)*N_slat
                        do nf_eff = 1, N_FL_eff
                           nf = Calc_Fl_map(nf_eff)
                           call op_wrapdo(gtt(:, :, nf, i_grc), op_v(n, nf), hs_field, ndim, n_Type, ntau)
                        end do
                     end do

                  end do

                  do ns = 1, N_slat
                     i_grc = ns + (i_wlk - 1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        call hop_mod_mmthl_1d2(gtt(:, :, nf, i_grc), nf, 1)
                        call hop_mod_mmthr_m1_1d2(gtt(:, :, nf, i_grc), nf, 1)
                     end do
                  end do

                 !! Propagate wave function
                  do ns = 1, N_slat
                     i_grc = ns + (i_wlk - 1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        call hop_mod_mmthlc_1d2(phi_l_m(nf_eff, i_grc)%U, nf, 1)
                        do n = N_op, 1, -1
                           call op_mmultr(phi_l_m(nf_eff, i_grc)%U, op_v(n, nf), nsigma_bp(i_wlk)%f(n, ntau), 'c', 1)
                        end do
                        call hop_mod_mmthlc_1d2(phi_l_m(nf_eff, i_grc)%U, nf, 1)
                     end do
                  end do

               end if

            end do

             !! call svd
            if (ntau1 .eq. stab_nt(nst) .and. ntau1 .ne. 0) then
               call re_orthonormalize_walkers(phi_l_m, 'N')

               do i_wlk = 1, N_wlk
                  if (weight_k(i_wlk) .gt. Zero) then

                    !! read phi_r from the 1st slat det storage
                     i_grc = 1 + (i_wlk - 1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        phi_r_m(nf_eff, i_wlk) = udvst(nst, nf_eff, i_grc)
                     end do

                     do ns = 1, N_slat
                        i_grc = ns + (i_wlk - 1)*N_slat
                        do nf_eff = 1, N_FL_eff
                           nf = Calc_Fl_map(nf_eff)
                          !! compute Green's function
                           call cgrp(z, gr(:, :, nf), phi_r_m(nf_eff, i_wlk), phi_l_m(nf_eff, i_grc))
                           det_vec(nf_eff) = z
                           call control_precisiong(gr(:, :, nf), gtt(:, :, nf, i_grc), Ndim)
                           gtt(:, :, nf, i_grc) = gr(:, :, nf)
                           udvst(nst, nf_eff, i_grc) = phi_l_m(nf_eff, i_grc)
                        end do
                        det_Vec(:) = det_Vec(:)*N_SUN
                        zph1 = exp(dcmplx(0.d0, aimag(sum(det_vec))))
                        zph2 = exp(dcmplx(0.d0, aimag(overlap_mc(i_grc))))
                        if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                        call control_Precisionp(zph1, zph2)
                        overlap_mc(i_grc) = sum(det_vec)
                     end do

                  end if
               end do
               call rescale_overlap(overlap_mc)
               nst = nst - 1
            end if

         end do

         call re_orthonormalize_walkers(phi_l_m, 'N')

         do i_wlk = 1, N_wlk

            if (weight_k(i_wlk) .gt. zero) then

               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  call phi_r_m(nf_eff, i_wlk)%reset('r', phi_bp_r(nf_eff, i_wlk)%U)
               end do

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                   !! compute Green's function
                     call cgrp(z, gr(:, :, nf), phi_r_m(nf_eff, i_wlk), phi_l_m(nf_eff, i_grc))
                     det_vec(nf_eff) = z
                     call control_precisiong(gr(:, :, nf), gtt(:, :, nf, i_grc), Ndim)
                     gtt(:, :, nf, i_grc) = gr(:, :, nf)
                  end do
                  det_Vec(:) = det_Vec(:)*N_SUN
                  zph1 = exp(dcmplx(0.d0, aimag(sum(det_vec))))
                  zph2 = exp(dcmplx(0.d0, aimag(overlap_mc(i_grc))))
                  if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                  call control_Precisionp(zph1, zph2)
                  overlap_mc(i_grc) = sum(det_vec)
               end do

            end if
         end do
         call rescale_overlap(overlap_mc)

         do ns = 1, N_slat
            do i_wlk = 1, N_wlk
               i_grc = ns + (i_wlk - 1)*N_slat
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  call udvst(nstm, nf_eff, i_grc)%reset('l', wf_l(nf, ns)%P)
               end do
            end do
         end do

         if ((ltau .eq. 1) .and. (nsw .gt. nwarmup)) call mc_measure_dyn(udvst, gtt, phi_bp_r, stab_nt, overlap_mc)

      end do

      ! deallocate tmp udv class
      do i_wlk = 1, N_wlk
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call phi_r_m(nf_eff, i_wlk)%dealloc
         end do
      end do

      deallocate (phi_r_m)
      deallocate (gtt)
      deallocate (overlap_mc)

   end subroutine metropolis

   subroutine mc_measure_dyn(udvst, gr_in, phi_r_in, stab_nt, overlap_in)
#ifdef MPI
      use mpi
#endif
      implicit none

      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(in) :: gr_in
      complex(Kind=kind(0.d0)), dimension(:), allocatable, intent(in) :: overlap_in
      class(udv_state), dimension(:, :), allocatable, intent(in) :: phi_r_in
      class(udv_state), dimension(:, :, :), allocatable, intent(in) :: udvst
      integer, dimension(:), allocatable, intent(in) :: stab_nt

      !Local
      complex(Kind=kind(0.d0)) :: temp(ndim, ndim), gr(ndim, ndim, n_fl), grc(ndim, ndim, n_fl)
      complex(Kind=kind(0.d0)) :: gt0(ndim, ndim, n_fl, n_grc), g00(ndim, ndim, n_fl, n_grc)
      complex(Kind=kind(0.d0)) :: gtt(ndim, ndim, n_fl, n_grc), g0t(ndim, ndim, n_fl, n_grc)
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
      integer :: ns, i_grc, act_mea, i_st, i_ed, n_part
      complex(Kind=kind(0.d0)) :: Z, Z_weight, DETZ, z_sum_overlap, exp_overlap(N_slat)
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8

      class(udv_state), dimension(:, :), allocatable :: phi_r_mea, phi_l_mea

#ifdef MPI
      integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
      integer        :: STATUS(MPI_STATUS_SIZE)
      call mpi_comm_size(MPI_COMM_WORLD, ISIZE, IERR)
      call mpi_comm_rank(MPI_COMM_WORLD, IRANK, IERR)
      call mpi_comm_rank(Group_Comm, irank_g, ierr)
      call mpi_comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

          !! initialization
      n_op = size(op_v, 1)
      nstm = size(udvst, 1)

          !! allocate tmp wavefunction
      allocate (phi_r_mea(N_FL_eff, N_wlk))
      allocate (phi_l_mea(N_FL_eff, N_grc))

      n_part = phi_r_in(1, 1)%n_part

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do i_wlk = 1, N_wlk
            call phi_r_mea(nf_eff, i_wlk)%alloc(ndim, n_part)
            phi_r_mea(nf_eff, i_wlk) = phi_r_in(nf_eff, i_wlk)
            do ns = 1, n_slat
               i_grc = ns + (i_wlk - 1)*N_slat
               call phi_l_mea(nf_eff, i_grc)%alloc(ndim, n_part)
            end do
         end do
      end do

      gtt = gr_in
      g00 = gtt
      gt0 = gtt
      g0t = gtt
      do i_grc = 1, N_grc
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do I = 1, Ndim
            g0t(I, I, nf, i_grc) = g0t(I, I, nf, i_grc) - 1.d0
         end do
      end do
      end do

      ntau = 0
          !! call reconstruction of non-calculated flavor blocks
      do i_grc = 1, n_grc
         if (reconstruction_needed) then
            call ham%gr_reconstruction(g00(:, :, :, i_grc))
            call ham%gr_reconstruction(gtt(:, :, :, i_grc))
            call ham%grt_reconstruction(gt0(:, :, :, i_grc), g0t(:, :, :, i_grc))
         end if
      end do
      call ham%obsert_mc(ntau, gt0, g0t, g00, gtt, overlap_in)

      NST = 1
      do ntau = 1, ltrot

         do i_wlk = 1, N_wlk
                !! Propagate wave function
            if (weight_k(i_wlk) .gt. Zero) then

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat

                   !! Propagate Green's function
                  call propr(gt0(:, :, :, i_grc), ntau, i_wlk)
                  call proprm1(g0t(:, :, :, i_grc), ntau, i_wlk)
                  call proprm1(gtt(:, :, :, i_grc), ntau, i_wlk)
                  call propr(gtt(:, :, :, i_grc), ntau, i_wlk)

               end do

               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  call Hop_mod_mmthr_1D2(phi_r_mea(nf_eff, i_wlk)%U, nf, 1)
                  do n = 1, N_op
                     call Op_mmultR(phi_r_mea(nf_eff, i_wlk)%U, Op_V(n, nf), nsigma_bp(i_wlk)%f(n, ntau), 'n', 1)
                  end do
                  call Hop_mod_mmthr_1D2(phi_r_mea(nf_eff, i_wlk)%U, nf, 1)
               end do

            end if

                !! call reconstruction of non-calculated flavor blocks
            do i_grc = 1, n_grc
               if (reconstruction_needed) then
                  call ham%gr_reconstruction(g00(:, :, :, i_grc))
                  call ham%gr_reconstruction(gtt(:, :, :, i_grc))
                  call ham%grt_reconstruction(gt0(:, :, :, i_grc), g0t(:, :, :, i_grc))
               end if
            end do
            call ham%obsert_mc(ntau, gt0, g0t, g00, gtt, overlap_in)

         end do

             !! call svd
         if (ntau .eq. stab_nt(nst)) then
            call re_orthonormalize_walkers(phi_r_mea, 'N')

            do i_wlk = 1, N_wlk

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     phi_l_mea(nf_eff, i_grc) = udvst(nst, nf_eff, i_grc)
                     call cgrp(detz, gr(:, :, nf), phi_r_mea(nf_eff, i_wlk), phi_l_mea(nf_eff, i_grc))
                     call control_precision_tau(gtt(:, :, nf, i_grc), gr(:, :, nf), Ndim)
                  end do
                  gtt(:, :, :, i_grc) = gr

                  grc = -gr
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     do I = 1, Ndim
                        grc(I, I, nf) = grc(I, I, nf) + 1.d0
                     end do
                  end do

                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     call mmult(temp, gr(:, :, nf), gt0(:, :, nf, i_grc))
                     gt0(:, :, nf, i_grc) = temp
                     call mmult(temp, g0t(:, :, nf, i_grc), grc(:, :, nf))
                     g0t(:, :, nf, i_grc) = temp
                  end do

               end do

            end do
            nst = nst + 1
         end if

      end do

      ! deallocate tmp udv class
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do i_wlk = 1, N_wlk
            call phi_r_mea(nf_eff, i_wlk)%dealloc
            do ns = 1, n_slat
               i_grc = ns + (i_wlk - 1)*N_slat
               call phi_l_mea(nf_eff, i_grc)%dealloc
            end do
         end do
      end do

      deallocate (phi_r_mea, phi_l_mea)

   end subroutine mc_measure_dyn

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
      integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
      integer :: ns, i_grc, act_mea, i_st, i_ed
      complex(Kind=kind(0.d0)) :: Z, Z_weight, DETZ, z_sum_overlap, exp_overlap(N_slat)
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
      real(Kind=kind(0.d0)) :: Zero = 1.0e-8

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
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            call CGRP(DetZ, GR(:, :, nf), phi_bp_r(nf_eff, i_wlk), phi_bp_l(nf_eff, i_grc))
         end do
         GTT(:, :, :, i_grc) = GR
      end do

      G00 = GTT
      GT0 = GTT
      G0T = GTT
      do i_grc = 1, N_grc
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do I = 1, Ndim
            G0T(I, I, nf, i_grc) = G0T(I, I, nf, i_grc) - 1.d0
         end do
      end do
      end do

      ntau = 0
      act_mea = 0 + irank
      do i_wlk = 1, N_wlk
         i_st = 1 + (i_wlk - 1)*N_slat; i_ed = i_wlk*N_slat
         exp_overlap(:) = exp(overlap(i_st:i_ed))
         z_sum_overlap = sum(exp_overlap(:))
         do ns = 1, N_slat
            i_grc = ns + (i_wlk - 1)*N_slat
                !! call reconstruction of non-calculated flavor blocks
            if (reconstruction_needed) then
               call ham%GR_reconstruction(G00(:, :, :, i_grc))
               call ham%GR_reconstruction(GTT(:, :, :, i_grc))
               call ham%GRT_reconstruction(GT0(:, :, :, i_grc), G0T(:, :, :, i_grc))
            end if
            call ham%obserT(ntau, GT0(:, :, :, i_grc), G0T(:, :, :, i_grc), G00(:, :, :, i_grc), &
                & GTT(:, :, :, i_grc), i_wlk, i_grc, z_weight, z_sum_overlap, act_mea)
            act_mea = act_mea + 1
         end do
      end do

      NST = 1
      do ntau = 1, ltrot

         do i_wlk = 1, N_wlk
                !! Propagate wave function
            if (weight_k(i_wlk) .gt. Zero) then

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat

                   !! Propagate Green's function
                  call PROPR(GT0(:, :, :, i_grc), ntau, i_wlk)
                  call PROPRM1(G0T(:, :, :, i_grc), ntau, i_wlk)
                  call PROPRM1(GTT(:, :, :, i_grc), ntau, i_wlk)
                  call PROPR(GTT(:, :, :, i_grc), ntau, i_wlk)

               end do

               do nf_eff = 1, N_FL_eff

                  nf = Calc_Fl_map(nf_eff)
                  call Hop_mod_mmthr_1D2(phi_bp_r(nf_eff, i_wlk)%U, nf, 1)

                  do n = 1, N_op
                     call Op_mmultR(phi_bp_r(nf_eff, i_wlk)%U, Op_V(n, nf), nsigma_bp(i_wlk)%f(n, ntau), 'n', 1)
                  end do

                  call Hop_mod_mmthr_1D2(phi_bp_r(nf_eff, i_wlk)%U, nf, 1)

               end do

            end if

                !! call reconstruction of non-calculated flavor blocks
            i_st = 1 + (i_wlk - 1)*N_slat; i_ed = i_wlk*N_slat
            exp_overlap(:) = exp(overlap(i_st:i_ed))
            z_sum_overlap = sum(exp_overlap(:))
            do ns = 1, N_slat
               i_grc = ns + (i_wlk - 1)*N_slat
               if (reconstruction_needed) then
                  call ham%GR_reconstruction(G00(:, :, :, i_grc))
                  call ham%GR_reconstruction(GTT(:, :, :, i_grc))
                  call ham%GRT_reconstruction(GT0(:, :, :, i_grc), G0T(:, :, :, i_grc))
               end if
               call ham%obserT(ntau, GT0(:, :, :, i_grc), G0T(:, :, :, i_grc), G00(:, :, :, i_grc), &
                   & GTT(:, :, :, i_grc), i_wlk, i_grc, z_weight, z_sum_overlap, act_mea)
               act_mea = act_mea + 1
            end do
         end do

             !! call svd
         if (ntau .eq. stab_nt(nst)) then
            call re_orthonormalize_walkers(phi_bp_r, 'N')

            do i_wlk = 1, N_wlk

               do ns = 1, N_slat
                  i_grc = ns + (i_wlk - 1)*N_slat
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     phi_bp_l(nf_eff, i_grc) = udvst(nst, nf_eff, i_grc)
                     call CGRP(DetZ, GR(:, :, nf), phi_bp_r(nf_eff, i_wlk), phi_bp_l(nf_eff, i_grc))
                     call Control_Precision_tau(GTT(:, :, nf, i_grc), GR(:, :, nf), Ndim)
                  end do
                  GTT(:, :, :, i_grc) = GR

                  GRC = -GR
                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     do I = 1, Ndim
                        GRC(I, I, nf) = GRC(I, I, nf) + 1.d0
                     end do
                  end do

                  do nf_eff = 1, N_FL_eff
                     nf = Calc_Fl_map(nf_eff)
                     call MMULT(TEMP, GR(:, :, nf), GT0(:, :, nf, i_grc))
                     GT0(:, :, nf, i_grc) = TEMP
                     call MMULT(TEMP, G0T(:, :, nf, i_grc), GRC(:, :, nf))
                     G0T(:, :, nf, i_grc) = TEMP
                  end do

               end do

            end do
            nst = nst + 1
         end if

      end do

   end subroutine bp_measure_tau

   subroutine rescale_overlap(overlap_in)

      implicit none

      complex(Kind=kind(0.d0)), dimension(:), allocatable, intent(inout) :: overlap_in

      integer :: nf, nf_eff, n, m, nt, i_wlk, i_grc, ns
      integer :: i_st, i_ed, ncslat
      real(Kind=kind(0.d0)) :: log_o_abs(n_slat), log_o_phase(n_slat), dz2
      real(kind=kind(0.d0)) :: pi = acos(-1.d0), dre_o, zero = 1.0e-8
      complex(Kind=kind(0.d0)) :: z1, zp

      do i_wlk = 1, N_wlk

         if (weight_k(i_wlk) .gt. zero) then

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
      integer :: nf, nf_eff, n

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
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
      integer :: nf, nf_eff, n

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call Hop_mod_mmthl_m1_1D2(Ain(:, :, nf), nf, nt)
         do n = 1, size(Op_V, 1)
            call Op_mmultL(Ain(:, :, nf), Op_V(n, nf), -nsigma_bp(i_wlk)%f(n, nt), 'n', nt)
         end do
         call Hop_mod_mmthl_m1_1D2(Ain(:, :, nf), nf, nt)
      end do

   end subroutine PROPRM1

   subroutine seed_vec_in(file_tg)
#ifdef MPI
      use mpi
#endif
      implicit none

      character(LEN=64), intent(IN)  :: FILE_TG

      ! LOCAL
      integer             :: K, ii, intseed(2), seed_tmp(2)
      character(LEN=64)  :: filename
      integer, allocatable :: SEED_VEC(:), seed_output(:, :), s_tmp(:)

#if defined(MPI)
      integer        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

      filename = file_tg

      call GET_SEED_LEN(K)
      allocate (SEED_VEC(K))
      allocate (s_tmp(K))

      if (irank_g .eq. 0) then
         allocate (seed_output(K, isize_g))
         open (UNIT=10, FILE=filename, STATUS='OLD', ACTION='Read')
         do ii = 1, isize_g
            read (10, *) seed_output(:, ii)
         end do
         close (10)
      end if

      if (irank_g .ne. 0) then
         call mpi_recv(seed_vec, K, mpi_integer, 0, 1, MPI_COMM_WORLD, STATUS, IERR)
         call RANSET(SEED_VEC)
      else
         do ii = 1, isize_g - 1
            s_tmp(:) = seed_output(:, ii + 1)
            call mpi_send(s_tmp, K, mpi_integer, ii, 1, MPI_COMM_WORLD, IERR)
         end do
      end if

      if (irank_g .eq. 0) then
         SEED_VEC = seed_output(:, 1)
         call RANSET(SEED_VEC)
         deallocate (seed_output)
      end if

      deallocate (seed_vec)
      deallocate (s_tmp)

   end subroutine seed_vec_in

   subroutine seed_vec_out
#ifdef MPI
      use mpi
#endif
      implicit none

      ! LOCAL
      integer             :: K, ii, intseed(2), seed_tmp(2)
      character(LEN=64)  :: FILE_TG, filename
      integer, allocatable :: SEED_VEC(:), seed_output(:, :), s_tmp(:)

#if defined(MPI)
      integer        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif
      filename = "seedvec_out"

      call GET_SEED_LEN(K)
      allocate (SEED_VEC(K))
      allocate (s_tmp(K))
      call RANGET(SEED_VEC)

      if (irank_g .eq. 0) then
         allocate (seed_output(K, isize_g))
         seed_output(:, 1) = SEED_VEC(:)
      end if

      if (irank_g .ne. 0) then
         call mpi_send(seed_vec, K, mpi_integer, 0, 1, MPI_COMM_WORLD, IERR)
      else
         do ii = 1, isize_g - 1
            call mpi_recv(s_tmp, K, mpi_integer, ii, 1, MPI_COMM_WORLD, STATUS, IERR)
            seed_output(:, ii + 1) = s_tmp(:)
         end do
      end if

      if (irank_g .eq. 0) then
         open (UNIT=10, FILE=filename, STATUS='UNKNOWN', ACTION='WRITE')
         do ii = 1, isize_g
            write (10, *) seed_output(:, ii)
         end do
         close (10)
         deallocate (seed_output)
      end if

      deallocate (seed_vec)
      deallocate (s_tmp)

   end subroutine seed_vec_out

#if defined HDF5
   subroutine wavefunction_out_hdf5(phi_0)
#ifdef MPI
      use mpi
#endif
#if defined HDF5
      use hdf5
      use h5lt
#endif
      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(IN) :: phi_0

      ! LOCAL
      character(LEN=64) :: FILE_TG, filename
      complex(Kind=kind(0.d0)), pointer :: phi0_out(:, :, :, :)
      complex(Kind=kind(0.d0)), pointer :: overlap_out(:)
      real(Kind=kind(0.d0)), pointer :: weight_out(:)
      complex(Kind=kind(0.d0)), allocatable :: otphi_tmp(:), p0_tmp(:, :, :, :), p1_tmp(:, :, :, :)
      real(Kind=kind(0.d0)), allocatable :: wt_tmp(:)

      integer             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii
      integer             :: i_st2, i_ed2
      integer(HSIZE_T), allocatable :: dims(:), dimsc(:)
      logical             :: file_exists
      integer(HID_T)      :: file_id, crp_list, space_id, dset_id, dataspace
      character(len=64)  :: dset_name
      type(c_ptr)         :: dat_ptr

#if defined(MPI)
      integer        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif
      filename = "phiout_0"

      write (filename, '(A,A)') trim(filename), ".h5"

      n_part = phi_0(1, 1)%n_part

      if (irank_g .eq. 0) then
         allocate (phi0_out(ndim, n_part, n_fl_eff, n_wlk_mpi))
         allocate (weight_out(N_wlk_mpi), overlap_out(N_grc_mpi))
      end if

      allocate (p0_tmp(ndim, n_part, n_fl_eff, n_wlk))
      allocate (p1_tmp(ndim, n_part, n_fl_eff, n_wlk))
      allocate (wt_tmp(N_wlk))
      allocate (otphi_tmp(N_grc))

      do nf = 1, N_FL_eff
      do nw = 1, N_wlk
         p0_tmp(:, :, nf, nw) = phi_0(nf, nw)%U(:, :)
      end do
      end do

      if (irank_g .ne. 0) then
         call mpi_send(overlap, N_grc, mpi_complex16, 0, 0, MPI_COMM_WORLD, IERR)
         call mpi_send(weight_k, N_wlk, mpi_real8, 0, 1, MPI_COMM_WORLD, IERR)
         Ndt = N_FL_eff*N_wlk*ndim*n_part
         call mpi_send(p0_tmp, Ndt, mpi_complex16, 0, 2, MPI_COMM_WORLD, IERR)
      else
         do ii = 1, isize_g - 1
            i_st = ii*N_wlk + 1; i_st2 = ii*N_grc + 1
            i_ed = (ii + 1)*N_wlk; i_ed2 = (ii + 1)*N_grc

            call mpi_recv(otphi_tmp, N_grc, mpi_complex16, ii, 0, MPI_COMM_WORLD, STATUS, IERR)
            overlap_out(i_st2:i_ed2) = otphi_tmp(:)
            call mpi_recv(wt_tmp, N_wlk, mpi_real8, ii, 1, MPI_COMM_WORLD, STATUS, IERR)
            weight_out(i_st:i_ed) = wt_tmp(:)
            Ndt = N_FL_eff*N_wlk*ndim*n_part
            call mpi_recv(p1_tmp, Ndt, mpi_complex16, ii, 2, MPI_COMM_WORLD, STATUS, IERR)
            phi0_out(:, :, :, i_st:i_ed) = p1_tmp
         end do
      end if

      if (irank_g .eq. 0) then
         i_st = 1; i_st2 = 1
         i_ed = N_wlk; i_ed2 = N_grc
         overlap_out(i_st2:i_ed2) = overlap(:)
         weight_out(i_st:i_ed) = weight_k(:)
         phi0_out(:, :, :, i_st:i_ed) = p0_tmp
      end if

      if (irank .eq. 0) then

         inquire (file=filename, exist=file_exists)
         if (.not. file_exists) then
            call h5open_f(ierr)
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferr)

            !Create and write dataset for field overlap
            dset_name = "overlap"
            rank = 2
            allocate (dims(2), dimsc(2))
            dims = [2, N_grc_mpi]
            dimsc = dims
            call h5screate_simple_f(rank, dims, space_id, hdferr)
            call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
            call h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
#if defined HDF5_ZLIB
            ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
            call h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
            !Create a dataset using cparms creation properties.
            call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                             dset_id, hdferr, crp_list)
            dat_ptr = c_loc(overlap_out(1))
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !Close objects
            deallocate (dims, dimsc)
            call h5sclose_f(space_id, hdferr)
            call h5pclose_f(crp_list, hdferr)
            call h5dclose_f(dset_id, hdferr)

            !Create and write dataset for real part of weight
            dset_name = "weight_re"
            rank = 2
            allocate (dims(2), dimsc(2))
            dims = [1, N_wlk_mpi]
            dimsc = dims
            call h5screate_simple_f(rank, dims, space_id, hdferr)
            call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
            call h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
#if defined HDF5_ZLIB
            ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
            call h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
            !Create a dataset using cparms creation properties.
            call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                             dset_id, hdferr, crp_list)
            dat_ptr = c_loc(weight_out(1))
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !Close objects
            deallocate (dims, dimsc)
            call h5sclose_f(space_id, hdferr)
            call h5pclose_f(crp_list, hdferr)
            call h5dclose_f(dset_id, hdferr)

            !Create and write dataset for wave function
            dset_name = "wavefunction"
            rank = 5
            allocate (dims(5), dimsc(5))
            dims = [2, ndim, n_part, N_FL_eff, N_wlk_mpi]
            dimsc = dims
            call h5screate_simple_f(rank, dims, space_id, hdferr)
            call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
            call h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
#if defined HDF5_ZLIB
            ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
            call h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
            !Create a dataset using cparms creation properties.
            call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                             dset_id, hdferr, crp_list)
            dat_ptr = c_loc(phi0_out(1, 1, 1, 1))
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !Close objects
            deallocate (dims, dimsc)
            call h5sclose_f(space_id, hdferr)
            call h5pclose_f(crp_list, hdferr)
            call h5dclose_f(dset_id, hdferr)

            !close file
            call h5fclose_f(file_id, hdferr)

         else
            !open file
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdferr)

            !open and write field overlap
            dset_name = "overlap"
            !Open the  dataset.
            call h5dopen_f(file_id, dset_name, dset_id, hdferr)
            dat_ptr = c_loc(overlap_out(1))
            !Write data
            call H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !close objects
            call h5dclose_f(dset_id, hdferr)

            !open and write real weight
            dset_name = "weight_re"
            !Open the  dataset.
            call h5dopen_f(file_id, dset_name, dset_id, hdferr)
            dat_ptr = c_loc(weight_out(1))
            !Write data
            call H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !close objects
            call h5dclose_f(dset_id, hdferr)

            !open and write real weight
            dset_name = "wavefunction"
            !Open the  dataset.
            call h5dopen_f(file_id, dset_name, dset_id, hdferr)
            dat_ptr = c_loc(phi0_out(1, 1, 1, 1))
            !Write data
            call H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !close objects
            call h5dclose_f(dset_id, hdferr)

            !close file
            call h5fclose_f(file_id, hdferr)
         end if

      end if !irank 0

      if (irank_g .eq. 0) then
         deallocate (phi0_out, weight_out, overlap_out)
      end if
      deallocate (p0_tmp, p1_tmp, wt_tmp, otphi_tmp)

   end subroutine wavefunction_out_hdf5

   subroutine wavefunction_in_hdf5(phi_0, file_tg)
#ifdef MPI
      use mpi
#endif
#if defined HDF5
      use hdf5
      use h5lt
#endif

      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0
      character(LEN=64), intent(IN)  :: FILE_TG

      ! LOCAL
      character(LEN=64) :: filename
      complex(Kind=kind(0.d0)), pointer :: phi0_out(:, :, :, :)
      complex(Kind=kind(0.d0)), pointer :: overlap_out(:)
      real(Kind=kind(0.d0)), pointer :: weight_out(:)
      complex(Kind=kind(0.d0)), allocatable :: otphi_tmp(:), p0_tmp(:, :, :, :), p1_tmp(:, :, :, :)
      real(Kind=kind(0.d0)), allocatable :: wt_tmp(:)

      integer             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii
      integer             :: nwalk_in, ngrc_in, i_st2, i_ed2
      integer             :: nf_eff
      integer(HSIZE_T), allocatable :: dims(:), dimsc(:), maxdims(:)
      logical             :: file_exists
      integer(HID_T)      :: file_id, crp_list, space_id, dset_id, dataspace
      character(len=64)  :: dset_name
      type(c_ptr)         :: dat_ptr

#if defined(MPI)
      integer        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif
      filename = file_tg

      n_part = phi_0(1, 1)%n_part

      allocate (p0_tmp(ndim, n_part, n_fl_eff, n_wlk))
      allocate (p1_tmp(ndim, n_part, n_fl_eff, n_wlk))
      allocate (wt_tmp(N_wlk))
      allocate (otphi_tmp(N_grc))

      if (irank .eq. 0) then

         !open file
         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdferr)

         !open and read overlap
         dset_name = "overlap"
         !Open the  dataset.
         call h5dopen_f(file_id, dset_name, dset_id, hdferr)
         !Get dataset's dataspace handle.
         call h5dget_space_f(dset_id, dataspace, hdferr)
         !Get dataspace's rank.
         call h5sget_simple_extent_ndims_f(dataspace, rank, ierr)
         allocate (dims(rank), maxdims(rank))
         !Get dataspace's dimensions.
         call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, ierr)
         ngrc_in = dims(rank)
              !! allocate !!
         allocate (overlap_out(ngrc_in))
              !!-----------!!
         dat_ptr = c_loc(overlap_out(1))
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close dataspace
         call h5sclose_f(dataspace, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)
         deallocate (dims, maxdims)

         !open and read weight
         dset_name = "weight_re"
         !Open the  dataset.
         call h5dopen_f(file_id, dset_name, dset_id, hdferr)
         !Get dataset's dataspace handle.
         call h5dget_space_f(dset_id, dataspace, hdferr)
         !Get dataspace's rank.
         call h5sget_simple_extent_ndims_f(dataspace, rank, ierr)
         allocate (dims(rank), maxdims(rank))
         !Get dataspace's dimensions.
         call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, ierr)
         nwalk_in = dims(rank)
              !! allocate !!
         allocate (phi0_out(ndim, n_part, n_fl_eff, nwalk_in))
         allocate (weight_out(nwalk_in))
              !!-----------!!
         dat_ptr = c_loc(weight_out(1))
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close dataspace
         call h5sclose_f(dataspace, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)
         deallocate (dims, maxdims)

         !open and read real weight
         dset_name = "wavefunction"
         !Open the  dataset.
         call h5dopen_f(file_id, dset_name, dset_id, hdferr)
         dat_ptr = c_loc(phi0_out(1, 1, 1, 1))
         !Write data
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)

         !close file
         call h5fclose_f(file_id, hdferr)

      end if !irank 0

      if (irank_g .eq. 0) then
         do ii = 1, isize_g - 1
            i_st = ii*N_wlk + 1; i_st2 = ii*N_grc + 1
            i_ed = (ii + 1)*N_wlk; i_ed2 = ii*N_grc + 1

            otphi_tmp(:) = overlap_out(i_st2:i_ed2)
            call mpi_send(otphi_tmp, N_grc, mpi_complex16, ii, 0, MPI_COMM_WORLD, IERR)
            wt_tmp(:) = weight_out(i_st:i_ed)
            call mpi_send(wt_tmp, N_wlk, mpi_real8, ii, 1, MPI_COMM_WORLD, IERR)
            Ndt = N_FL_eff*N_wlk*ndim*n_part
            p1_tmp = phi0_out(:, :, :, i_st:i_ed)
            call mpi_send(p1_tmp, Ndt, mpi_complex16, ii, 2, MPI_COMM_WORLD, IERR)
         end do
      else
         call mpi_recv(otphi_tmp, N_grc, mpi_complex16, 0, 0, MPI_COMM_WORLD, STATUS, IERR)
         overlap(:) = otphi_tmp(:)
         call mpi_recv(wt_tmp, N_wlk, mpi_real8, 0, 1, MPI_COMM_WORLD, STATUS, IERR)
         weight_k(:) = wt_tmp(:)
         Ndt = N_FL_eff*N_wlk*ndim*n_part
         call mpi_recv(p0_tmp, Ndt, mpi_complex16, 0, 2, MPI_COMM_WORLD, STATUS, IERR)
         do nf = 1, N_FL_eff
         do nw = 1, N_wlk
            phi_0(nf, nw)%U(:, :) = p0_tmp(:, :, nf, nw)
         end do
         end do
      end if

      if (irank_g .eq. 0) then
         i_st = 1; i_st2 = 1
         i_ed = N_wlk; i_ed2 = N_grc
         p0_tmp = phi0_out(:, :, :, i_st:i_ed)

         weight_k(:) = weight_out(i_st:i_ed)
         overlap(:) = overlap_out(i_st2:i_ed2)
         do nf = 1, N_FL_eff
         do nw = 1, N_wlk
            phi_0(nf, nw)%U(:, :) = p0_tmp(:, :, nf, nw)
         end do
         end do
      end if

      if (irank_g .eq. 0) then
         deallocate (phi0_out, weight_out, overlap_out)
      end if
      deallocate (p0_tmp, p1_tmp, wt_tmp, otphi_tmp)

   end subroutine wavefunction_in_hdf5

   subroutine trial_in_hdf5(phi_0_r, phi_0_l, file_inst, file_antiinst)
#ifdef MPI
      use mpi
#endif
#if defined HDF5
      use hdf5
      use h5lt
#endif

      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0_r, phi_0_l
      character(LEN=64), intent(IN)  :: file_inst, file_antiinst

      ! LOCAL
      character(LEN=64) :: filename
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), pointer :: phi0_in(:, :, :)
      complex(kind=kind(0.d0)), allocatable :: p1_tmp(:, :, :), p2_tmp(:, :, :), log_zdet(:)
      complex(kind=kind(0.d0)), allocatable, dimension(:, :) :: sMat
      complex(kind=kind(0.d0)) :: alpha, beta, zdet, phase, t_overlap, z_norm, c1, ctmp
      real(kind=kind(0.d0)) :: d_norm

      integer             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii, nwalk_in
      integer             :: nf_eff, n_fl_in, ns, i_wlk, info, n
      integer(HSIZE_T), allocatable :: dims(:), dimsc(:), maxdims(:)
      logical             :: file_exists
      integer(HID_T)      :: file_id, crp_list, space_id, dset_id, dataspace
      character(len=64)  :: dset_name
      type(c_ptr)         :: dat_ptr

#if defined(MPI)
      integer        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
      integer        :: STATUS(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif
      n_part = phi_0_l(1, 1)%n_part

      n_fl_in = 1

      allocate (p1_tmp(ndim, n_part, n_fl))
      allocate (p2_tmp(ndim, n_part, n_fl))

      if (irank .eq. 0) then

              !! inst
         !open file
         call h5fopen_f(file_inst, H5F_ACC_RDWR_F, file_id, hdferr)

              !! allocate !!
         allocate (phi0_in(ndim, n_part, n_fl_in))
              !!-----------!!
         !open and read real weight
         dset_name = "wavefunction"
         !Open the  dataset.
         call h5dopen_f(file_id, dset_name, dset_id, hdferr)
         dat_ptr = c_loc(phi0_in(1, 1, 1))
         !Write data
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)

         !close file
         call h5fclose_f(file_id, hdferr)

              !! store input
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            p1_tmp(:, :, nf) = phi0_in(:, :, 1)
         end do

              !! anti-inst
         !open file
         call h5fopen_f(file_antiinst, H5F_ACC_RDWR_F, file_id, hdferr)

              !!-----------!!
         !open and read real weight
         dset_name = "wavefunction"
         !Open the  dataset.
         call h5dopen_f(file_id, dset_name, dset_id, hdferr)
         dat_ptr = c_loc(phi0_in(1, 1, 1))
         !Write data
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)

         !close file
         call h5fclose_f(file_id, hdferr)

              !! store input
         do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            p2_tmp(:, :, nf) = phi0_in(:, :, 1)
         end do

      end if !irank 0

      if (irank_g .eq. 0) then
         do ii = 1, isize_g - 1
            Ndt = N_FL*ndim*n_part
            call mpi_send(p1_tmp, Ndt, mpi_complex16, ii, 2, MPI_COMM_WORLD, IERR)
            call mpi_send(p2_tmp, Ndt, mpi_complex16, ii, 3, MPI_COMM_WORLD, IERR)
         end do
      else
         Ndt = N_FL*ndim*n_part
         call mpi_recv(p1_tmp, Ndt, mpi_complex16, 0, 2, MPI_COMM_WORLD, STATUS, IERR)
         call mpi_recv(p2_tmp, Ndt, mpi_complex16, 0, 3, MPI_COMM_WORLD, STATUS, IERR)
      end if

         !! allocate tmp matrix
      allocate (smat(n_part, n_part), ipiv(N_part), log_zdet(N_slat))

         !! test whether overlap has minus sign
      alpha = 1.d0
      beta = 0.d0
      ctmp = 0.d0
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call zgemm('C','N',N_part,N_part,Ndim,alpha,p1_tmp(1,1,nf),Ndim,phi_0_r(nf_eff,1)%U(1,1),ndim,beta,smat(1,1),N_part)
         ! ZGETRF computes an LU factorization of a general M-by-N matrix A
         ! using partial pivoting with row interchanges.
         call ZGETRF(N_part, N_part, smat, N_part, ipiv, info)
         ! obtain log of det
         zdet = 0.d0
         phase = 1.d0
         do n = 1, N_part
            if (ipiv(n) .ne. n) then
               phase = -phase
            end if
            zdet = zdet + log(smat(n, n))
         end do
         zdet = zdet + log(phase)

         ctmp = ctmp + zdet
      end do
      z_norm = exp(ctmp)
      d_norm = dble(z_norm)
      if (d_norm .lt. 0.d0) then
         c1 = cmplx(0.d0, 1.d0, kind(1.d0))
      else
         c1 = cmplx(1.d0, 0.d0, kind(1.d0))
      end if

         !! combine the trial wave function by
         !! | \Psi_T > = c1*| \Psi_T^1 > + c1^{\dagger}*| \Psi_T^2 >
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         phi_0_l(nf_eff, 1)%U(:, :) = p1_tmp(:, :, nf)*c1
         phi_0_l(nf_eff, 2)%U(:, :) = p2_tmp(:, :, nf)*conjg(c1)
         WF_L(nf, 1)%P(:, :) = p1_tmp(:, :, nf)*c1
         WF_L(nf, 2)%P(:, :) = p2_tmp(:, :, nf)*conjg(c1)
      end do

         !! normalization of overlap <\Psi_T | \phi_k^0>

      alpha = 1.d0
      beta = 0.d0
      log_zdet(:) = 0.d0
      do ns = 1, n_slat
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
     call zgemm('C','N',N_part,N_part,Ndim,alpha,phi_0_l(nf_eff,ns)%U(1,1),Ndim,phi_0_r(nf_eff,1)%U(1,1),ndim,beta,smat(1,1),N_part)
         ! ZGETRF computes an LU factorization of a general M-by-N matrix A
         ! using partial pivoting with row interchanges.
         call ZGETRF(N_part, N_part, smat, N_part, ipiv, info)
         ! obtain log of det
         zdet = 0.d0
         phase = 1.d0
         do n = 1, N_part
            if (ipiv(n) .ne. n) then
               phase = -phase
            end if
            zdet = zdet + log(smat(n, n))
         end do
         zdet = zdet + log(phase)

         log_zdet(ns) = log_zdet(ns) + zdet
      end do
      end do

      z_norm = exp(log_zdet(1)) + exp(log_zdet(2))
      if (irank_g .eq. 0) write (*, *) 'overlap for input slater', z_norm
      !d_norm = dble(z_norm)
      !z_norm = (1.d0/d_norm)**(1.d0/dble(2*N_part))
      z_norm = (1.d0/z_norm)**(1.d0/dble(2*N_part))
      if (irank_g .eq. 0) write (*, *) 'renormalized factor for input slater z_norm', z_norm

      do i_wlk = 1, n_wlk
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         phi_0_r(nf_eff, i_wlk)%U(:, :) = phi_0_r(nf_eff, i_wlk)%U(:, :)*z_norm
      end do
      end do

      if (irank_g .eq. 0) then
         deallocate (phi0_in)
      end if
      deallocate (p1_tmp)
      deallocate (p2_tmp)
      deallocate (smat, log_zdet, ipiv)

   end subroutine trial_in_hdf5
#endif

end module stepwlk_mod
