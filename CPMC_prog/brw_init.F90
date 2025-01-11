module BRW_init_mod

   use hamiltonian_main
   use udv_state_mod
   use gfun_mod
   use set_random
   use fields_mod
   use operator_mod
#ifdef MPI
   use mpi
#endif
   
contains
   
   subroutine initial_setup(ltrot_bp, nwrap, ltau)
      
       implicit none
       integer, intent(in) :: ltrot_bp, nwrap, ltau

       !local
       integer :: n_op, i, seed_in, i_wlk, n
       character(len=64) :: file_seeds
   
       call Fields_init()

       n_op = size(op_v, 1)

       File_seeds = "seeds"
       call Set_Random_number_Generator(File_seeds, Seed_in)

       !! To do: modify nsigma structure
       allocate (nsigma_bp(n_wlk))
       if (ltau .eq. 1) ltrot_bp = ltrot_bp + ltrot
       do i_wlk = 1, N_wlk
          call nsigma_bp(i_wlk)%make(n_op, ltrot_bp)
          do n = 1, N_op
             nsigma_bp(i_wlk)%t(n) = op_v(n, 1)%type
          end do
       end do

       call Hop_mod_init

       !! init log of weight
       weight_k(:) = cmplx(0.d0, 0.d0, kind(0.d0))
   
   end subroutine initial_setup

   subroutine initial_wlk(phi_trial, phi_0, phi_bp_l, phi_bp_r, udvst, stab_nt, GR, nwrap)
      
      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(INOUT) :: phi_trial, phi_0, phi_bp_l, phi_bp_r
      class(udv_state), dimension(:, :, :), allocatable, intent(INOUT) :: udvst
      complex(Kind=kind(0.d0)), dimension(:, :, :, :), allocatable, intent(INOUT) :: GR
      integer, dimension(:), allocatable, intent(INOUT) :: stab_nt
      integer, intent(IN) :: nwrap

      !Local
      integer :: nf, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_grc, NSTM, NST, ltrot_bp, ns
      integer :: i_st, i_ed, ncslat
      complex(Kind=kind(0.d0)) :: overlap_old, overlap_new, Z, Z1, Z2, tot_ene, ZP
      complex(Kind=kind(0.d0)) :: tot_c_weight, el_tmp
      complex(Kind=kind(0.d0)) :: det_Vec(N_FL)
      real(Kind=kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio, X1, wtmp, sign_w
      real(Kind=kind(0.d0)) :: zero = 1.0e-12, tot_re_weight, dz2
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
         do nf = 1, N_FL
            call phi_0   (nf, i_wlk)%init(ndim, 'r', WF_R(nf, 1)%P)
            call phi_bp_r(nf, i_wlk)%init(ndim, 'r', WF_R(nf, 1)%P)
         end do
      end do

      do ns = 1, N_slat
         do nf = 1, N_FL
            call phi_trial(nf, ns)%init(ndim, 'l', WF_L(nf, ns)%P)
         end do
         do i_wlk = 1, N_wlk
            i_grc = ns + (i_wlk - 1)*N_slat
            do nf = 1, N_FL
               do n = 1, nstm
                  call udvst(n, nf, i_grc)%alloc(ndim)
               end do
               call phi_bp_l(nf, i_grc)%init(ndim, 'l', wf_l(nf, ns)%P)
            end do
         end do
      end do

      file_tg = 'trial_0.h5'
      inquire (file=file_tg, exist=lconf_h5)
      if (lconf_h5) then
         if (irank_g .eq. 0) write (*, *) "read input trial wave function"
         call trial_in_hdf5(phi_0, phi_trial, file_tg)
      end if

      file_tg = "phiin_0.h5"
      inquire (FILE=file_tg, EXIST=LCONF_H5)
      if (LCONF_H5) then
         if (irank_g .eq. 0) write (*, *) "read input walkers"
         call wavefunction_in_hdf5(phi_0, file_tg)
      end if

          !! initial overlap and green's function
      do i_wlk = 1, N_wlk

         do ns = 1, N_slat
            i_grc = ns + (i_wlk - 1)*N_slat

            do nf = 1, N_Fl
               call cgrp(Z, gr(:, :, nf, i_grc), phi_0(nf, i_wlk), phi_trial(nf, ns))
               det_vec(nf) = Z
            end do
            det_vec(:) = det_Vec(:)*N_SUN
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

#if defined HDF5
   subroutine wavefunction_out_hdf5(phi_0)
#if defined HDF5
      use hdf5
      use h5lt
#endif
      implicit none

      class(udv_state), dimension(:, :), allocatable, intent(IN) :: phi_0

      ! LOCAL
      character(LEN=64) :: FILE_TG, filename
      complex(Kind=kind(0.d0)), pointer :: phi0_up_out(:, :, :), phi0_dn_out(:, :, :)
      complex(Kind=kind(0.d0)), pointer :: overlap_out(:)
      complex(Kind=kind(0.d0)), pointer :: weight_out(:)
      complex(Kind=kind(0.d0)), allocatable :: otphi_tmp(:), p0_tmp(:, :, :), p1_tmp(:, :, :)
      complex(Kind=kind(0.d0)), allocatable :: wt_tmp(:)

      integer             :: K, hdferr, rank, nf, nw, i0, i1, i2, i_st, i_ed, Ndt, ii
      integer             :: i_st2, i_ed2, n_part_1, n_part_2
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
      
      n_part_1 = phi_0(1, 1)%n_part
      n_part_2 = phi_0(2, 1)%n_part

      if (irank_g .eq. 0) then
         allocate (phi0_up_out(ndim, n_part_1, n_wlk_mpi))
         allocate (phi0_dn_out(ndim, n_part_2, n_wlk_mpi))
         allocate (weight_out(N_wlk_mpi), overlap_out(N_grc_mpi))
      end if

      allocate (p0_tmp(ndim, n_part_1, n_wlk))
      allocate (p1_tmp(ndim, n_part_2, n_wlk))
      allocate (wt_tmp(N_wlk))
      allocate (otphi_tmp(N_grc))

      if (irank_g .ne. 0) then
          do nw = 1, N_wlk
             p0_tmp(:, :, nw) = phi_0(1, nw)%U(:, :)
             p1_tmp(:, :, nw) = phi_0(2, nw)%U(:, :)
          end do
      endif

      if (irank_g .ne. 0) then
         call mpi_send(overlap, N_grc, mpi_complex16, 0, 0, MPI_COMM_WORLD, IERR)
         call mpi_send(weight_k, N_wlk, mpi_real8, 0, 1, MPI_COMM_WORLD, IERR)
         Ndt = N_wlk*ndim*n_part_1
         call mpi_send(p0_tmp, Ndt, mpi_complex16, 0, 2, MPI_COMM_WORLD, IERR)
         Ndt = N_wlk*ndim*n_part_2
         call mpi_send(p1_tmp, Ndt, mpi_complex16, 0, 3, MPI_COMM_WORLD, IERR)
      else
         do ii = 1, isize_g - 1
            i_st = ii*N_wlk + 1; i_st2 = ii*N_grc + 1
            i_ed = (ii + 1)*N_wlk; i_ed2 = (ii + 1)*N_grc

            call mpi_recv(otphi_tmp, N_grc, mpi_complex16, ii, 0, MPI_COMM_WORLD, STATUS, IERR)
            overlap_out(i_st2:i_ed2) = otphi_tmp(:)
            call mpi_recv(wt_tmp, N_wlk, mpi_complex16, ii, 1, MPI_COMM_WORLD, STATUS, IERR)
            weight_out(i_st:i_ed) = wt_tmp(:)
            Ndt = N_wlk*ndim*n_part_1
            call mpi_recv(p0_tmp, Ndt, mpi_complex16, ii, 2, MPI_COMM_WORLD, STATUS, IERR)
            phi0_up_out(:, :, i_st:i_ed) = p0_tmp
            Ndt = N_wlk*ndim*n_part_2
            call mpi_recv(p1_tmp, Ndt, mpi_complex16, ii, 3, MPI_COMM_WORLD, STATUS, IERR)
            phi0_dn_out(:, :, i_st:i_ed) = p1_tmp
         end do
      end if

      if (irank_g .eq. 0) then
          i_st = 1; i_st2 = 1
          i_ed = N_wlk; i_ed2 = N_grc
          overlap_out(i_st2:i_ed2) = overlap(:)
          weight_out(i_st:i_ed) = weight_k(:)

          do nw = 1, N_wlk
             p0_tmp(:, :, nw) = phi_0(1, nw)%U(:, :)
             p1_tmp(:, :, nw) = phi_0(2, nw)%U(:, :)
          end do
          phi0_up_out(:, :, i_st:i_ed) = p0_tmp
          phi0_dn_out(:, :, i_st:i_ed) = p1_tmp
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
!!#if defined HDF5_ZLIB
!!            ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
!!            call h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
!!#endif
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
            dims = [2, N_wlk_mpi]
            dimsc = dims
            call h5screate_simple_f(rank, dims, space_id, hdferr)
            call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
            call h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
!!#if defined HDF5_ZLIB
!!            ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
!!            call h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
!!#endif
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
            dset_name = "phi_trial_up"
            rank = 4
            allocate (dims(4), dimsc(4))
            dims = [2, ndim, n_part_1, N_wlk_mpi]
            dimsc = dims
            call h5screate_simple_f(rank, dims, space_id, hdferr)
            call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
            call h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
!!#if defined HDF5_ZLIB
!!            ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
!!            call h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
!!#endif
            !Create a dataset using cparms creation properties.
            call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                             dset_id, hdferr, crp_list)
            dat_ptr = c_loc(phi0_up_out(1, 1, 1))
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !Close objects
            deallocate (dims, dimsc)
            call h5sclose_f(space_id, hdferr)
            call h5pclose_f(crp_list, hdferr)
            call h5dclose_f(dset_id, hdferr)
            
            dset_name = "phi_trial_dn"
            rank = 4
            allocate (dims(4), dimsc(4))
            dims = [2, ndim, n_part_2, N_wlk_mpi]
            dimsc = dims
            call h5screate_simple_f(rank, dims, space_id, hdferr)
            call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
            call h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
!!#if defined HDF5_ZLIB
!!            ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
!!            call h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
!!#endif
            !Create a dataset using cparms creation properties.
            call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                             dset_id, hdferr, crp_list)
            dat_ptr = c_loc(phi0_dn_out(1, 1, 1))
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

            !open and write wave function
            dset_name= "phi_trial_up"
            !Open the  dataset.
            call h5dopen_f(file_id, dset_name, dset_id, hdferr)
            dat_ptr = c_loc(phi0_up_out(1, 1, 1))
            !Write data
            call H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !close objects
            call h5dclose_f(dset_id, hdferr)
            
            dset_name= "phi_trial_dn"
            !Open the  dataset.
            call h5dopen_f(file_id, dset_name, dset_id, hdferr)
            dat_ptr = c_loc(phi0_dn_out(1, 1, 1))
            !Write data
            call H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
            !close objects
            call h5dclose_f(dset_id, hdferr)

            !close file
            call h5fclose_f(file_id, hdferr)
         end if

      end if !irank 0

      if (irank_g .eq. 0) then
         deallocate (phi0_up_out, phi0_dn_out, weight_out, overlap_out)
      end if
      deallocate (p0_tmp, p1_tmp, wt_tmp, otphi_tmp)

   end subroutine wavefunction_out_hdf5

   subroutine wavefunction_in_hdf5(phi_0, file_tg)
#if defined HDF5
      use hdf5
      use h5lt
#endif

      implicit none

      class(UDV_State), dimension(:, :), allocatable, intent(INOUT) :: phi_0
      character(LEN=64), intent(IN)  :: FILE_TG

      ! LOCAL
      character(LEN=64) :: filename
      complex(Kind=kind(0.d0)), pointer :: phi0_up_out(:, :, :), phi0_dn_out(:, :, :)
      complex(Kind=kind(0.d0)), pointer :: overlap_out(:)
      complex(Kind=kind(0.d0)), pointer :: weight_out(:)
      complex(Kind=kind(0.d0)), allocatable :: otphi_tmp(:), p0_tmp(:, :, :), p1_tmp(:, :, :)
      complex(Kind=kind(0.d0)), allocatable :: wt_tmp(:)

      integer             :: K, hdferr, rank, nf, nw, i0, i1, i2, i_st, i_ed, Ndt, ii
      integer             :: nwalk_in, ngrc_in, i_st2, i_ed2, n_part_1, n_part_2
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

      n_part_1 = phi_0(1, 1)%n_part
      n_part_2 = phi_0(2, 1)%n_part

      allocate (p0_tmp(ndim, n_part_1, n_wlk))
      allocate (p1_tmp(ndim, n_part_2, n_wlk))
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
         allocate (phi0_up_out(ndim, n_part_1, nwalk_in))
         allocate (phi0_dn_out(ndim, n_part_2, nwalk_in))
         allocate (weight_out (nwalk_in))
              !!-----------!!
         dat_ptr = c_loc(weight_out(1))
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close dataspace
         call h5sclose_f(dataspace, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)
         deallocate (dims, maxdims)

         !open and read trial wave function
         dset_name= "phi_trial_up"
         !Open the  dataset.
         call h5dopen_f(file_id, dset_name, dset_id, hdferr)
         dat_ptr = c_loc(phi0_up_out(1, 1, 1))
         !Write data
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)
         
         dset_name= "phi_trial_dn"
         !Open the  dataset.
         call h5dopen_f(file_id, dset_name, dset_id, hdferr)
         dat_ptr = c_loc(phi0_dn_out(1, 1, 1))
         !Write data
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
         !close objects
         call h5dclose_f(dset_id, hdferr)

         !close file
         call h5fclose_f(file_id, hdferr)

      end if !irank 0

      if (irank_g .eq. 0) then
         do ii = 1, isize_g - 1
            i_st = ii*N_wlk+1  ; i_st2=ii*N_grc+1
            i_ed = (ii+1)*N_wlk; i_ed2=(ii+1)*N_grc

            otphi_tmp(:) = overlap_out(i_st2:i_ed2)
            call mpi_send(otphi_tmp, N_grc, mpi_complex16, ii, 0, MPI_COMM_WORLD, IERR)
            wt_tmp(:) = weight_out(i_st:i_ed)
            call mpi_send(wt_tmp, N_wlk, mpi_complex16, ii, 1, MPI_COMM_WORLD, IERR)
            Ndt = N_wlk*ndim*n_part_1
            p0_tmp = phi0_up_out(:, :, i_st:i_ed)
            call mpi_send(p0_tmp, Ndt, mpi_complex16, ii, 2, MPI_COMM_WORLD, IERR)
            Ndt = N_wlk*ndim*n_part_2
            p1_tmp = phi0_dn_out(:, :, i_st:i_ed)
            call mpi_send(p1_tmp, Ndt, mpi_complex16, ii, 3, MPI_COMM_WORLD, IERR)
         end do
      else
         call mpi_recv(otphi_tmp, N_grc, mpi_complex16, 0, 0, MPI_COMM_WORLD, STATUS, IERR)
         overlap(:) = otphi_tmp(:)
         call mpi_recv(wt_tmp, N_wlk, mpi_complex16, 0, 1, MPI_COMM_WORLD, STATUS, IERR)
         weight_k(:) = wt_tmp(:)
         Ndt = N_wlk*ndim*n_part_1
         call mpi_recv(p0_tmp, Ndt, mpi_complex16, 0, 2, MPI_COMM_WORLD, STATUS, IERR)
         Ndt = N_wlk*ndim*n_part_2
         call mpi_recv(p1_tmp, Ndt, mpi_complex16, 0, 3, MPI_COMM_WORLD, STATUS, IERR)
         do nw = 1, N_wlk
            phi_0(1, nw)%U(:, :) = p0_tmp(:, :, nw)
            phi_0(2, nw)%U(:, :) = p1_tmp(:, :, nw)
         end do
      end if

      if (irank_g .eq. 0) then
         i_st = 1    ; i_st2 = 1
         i_ed = N_wlk; i_ed2 = N_grc
         p0_tmp = phi0_up_out(:, :, i_st:i_ed)
         p1_tmp = phi0_dn_out(:, :, i_st:i_ed)

         weight_k(:) = weight_out(i_st:i_ed)
         overlap(:) = overlap_out(i_st2:i_ed2)
         do nw = 1, N_wlk
            phi_0(1, nw)%U(:, :) = p0_tmp(:, :, nw)
            phi_0(2, nw)%U(:, :) = p1_tmp(:, :, nw)
         end do
      end if

      if (irank_g .eq. 0) then
         deallocate (phi0_up_out, phi0_dn_out, weight_out, overlap_out)
      end if
      deallocate (p0_tmp, p1_tmp, wt_tmp, otphi_tmp)

   end subroutine wavefunction_in_hdf5

   subroutine trial_in_hdf5( phi_0_r, phi_0_l, file_tg )
#if defined HDF5
      use hdf5
      use h5lt
#endif
          
      implicit none
     
      class(udv_state), dimension(:,:), allocatable, intent(inout) :: phi_0_r, phi_0_l
      character (LEN=64), intent(in)  :: file_tg

      ! LOCAL
      character (LEN=64) :: filename
      complex (Kind=Kind(0.d0)), pointer :: phi0_up_out(:,:), phi0_dn_out(:,:)
      complex (Kind=Kind(0.d0)), allocatable :: p0_tmp(:,:), p1_tmp(:,:), p2_tmp(:,:), p3_tmp(:,:)
      complex (kind=kind(0.d0)), allocatable, dimension(:,:) :: smat_up, smat_dn
      complex (kind=kind(0.d0)) :: alpha, beta, zdet, z_norm_up, z_norm_dn, z_norm
      complex (kind=kind(0.d0)) :: phase
      integer, allocatable :: ipiv_up(:), ipiv_dn(:)

      INTEGER             :: K, hdferr, rank, nf, nw, i0, i1, i2, i_st, i_ed, Ndt, ii, nwalk_in
      Integer             :: n_part_1, n_part_2, n, info
      INTEGER(HSIZE_T), allocatable :: dims(:), dimsc(:), maxdims(:)
      Logical             :: file_exists
      INTEGER(HID_T)      :: file_id, crp_list, space_id, dset_id, dataspace
      Character (len=64)  :: dset_name
      TYPE(C_PTR)         :: dat_ptr

#if defined(MPI)
      INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
      Integer        :: STATUS(MPI_STATUS_SIZE)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup           = irank/isize_g
#endif
      filename = file_tg

      n_part_1 = phi_0_l(1, 1)%n_part
      n_part_2 = phi_0_l(2, 1)%n_part
      
      allocate(p0_tmp(ndim,n_part_1))
      allocate(p1_tmp(ndim,n_part_2))
      allocate(p2_tmp(ndim,n_part_1))
      allocate(p3_tmp(ndim,n_part_2))
      
      if ( irank .eq. 0 ) then

          !open file
          CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

          !! allocate !!
          allocate(phi0_up_out(ndim,n_part_1))
          allocate(phi0_dn_out(ndim,n_part_2))
          !!-----------!!
          !open and read wave function
          dset_name= "phi_trial_up"
          !Open the  dataset.
          call h5dopen_f(file_id, dset_name, dset_id, hdferr)
          dat_ptr = c_loc(phi0_up_out(1,1))
          !Write data
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
          !close objects
          call h5dclose_f(dset_id,   hdferr)
          
          dset_name= "phi_trial_dn"
          !Open the  dataset.
          call h5dopen_f(file_id, dset_name, dset_id, hdferr)
          dat_ptr = c_loc(phi0_dn_out(1,1))
          !Write data
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
          !close objects
          call h5dclose_f(dset_id,   hdferr)
            
          !close file
          call h5fclose_f(file_id, hdferr)

          p0_tmp(:,:)=phi0_up_out(:,:)
          p1_tmp(:,:)=phi0_dn_out(:,:)

      endif !irank 0

      if ( irank_g .eq. 0 ) then
         do ii = 1, isize_g-1
            Ndt=ndim*n_part_1
            call mpi_send(p0_tmp,  Ndt,mpi_complex16, ii, 2, MPI_COMM_WORLD,IERR)
            Ndt=ndim*n_part_2
            call mpi_send(p1_tmp,  Ndt,mpi_complex16, ii, 3, MPI_COMM_WORLD,IERR)
         ENDDO
      else
         Ndt=ndim*n_part_1
         call mpi_recv(p0_tmp,  Ndt,mpi_complex16, 0, 2, MPI_COMM_WORLD,STATUS,IERR)
         Ndt=ndim*n_part_2
         call mpi_recv(p1_tmp,  Ndt,mpi_complex16, 0, 3, MPI_COMM_WORLD,STATUS,IERR)
      ENDIF

      p2_tmp = p0_tmp
      p3_tmp = p1_tmp

      !! allocate tmp matrix
      allocate(smat_up(n_part_1,n_part_1), smat_dn(n_part_2,n_part_2))
      allocate(ipiv_up(n_part_1), ipiv_dn(n_part_2) )

      !! compute overlap
      alpha = 1.d0
      beta  = 0.d0
      call zgemm('C','N',n_part_1,n_part_1,ndim,alpha,p0_tmp(1,1), ndim, p2_tmp(1,1), &
          & ndim,beta,smat_up(1,1), n_part_1)
      ! ZGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call zgetrf(n_part_1, n_part_1, smat_up, n_part_1, ipiv_up, info)
      ! obtain log of det
      zdet  = 0.d0
      phase = 1.d0
      do n = 1, n_part_1
         if (ipiv_up(n).ne.n) then
            phase = -phase
         endif
         zdet = zdet + log(smat_up(n,n))
      enddo
      zdet = zdet + log(phase)
      z_norm_up = exp(zdet)
      if (irank_g .eq. 0 ) write(*,*) '======================================================================'
      if (irank_g .eq. 0 ) write(*,*) 'overlap for input slater spin up', z_norm_up
      z_norm_up = (1.d0/z_norm_up)**(1.d0/dble(N_part_1))
      if (irank_g .eq. 0 ) write(*,*) 'renormalized factor for input slater z_norm spin up', z_norm_up
      if (irank_g .eq. 0 ) write(*,*) '======================================================================'

      call zgemm('C','N',n_part_2,n_part_2,ndim,alpha,p1_tmp(1,1), ndim, p3_tmp(1,1), &
          & ndim, beta, smat_dn(1,1), n_part_2)
      ! ZGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call zgetrf(n_part_2, n_part_2, smat_dn, n_part_2, ipiv_dn, info)
      ! obtain log of det
      zdet  = 0.d0
      phase = 1.d0
      do n=1,n_part_2
         if (ipiv_dn(n).ne.n) then
            phase = -phase
         endif
         zdet = zdet + log(smat_dn(n,n))
      enddo
      zdet = zdet + log(phase)
      z_norm_dn = exp(zdet)
      if (irank_g .eq. 0 ) write(*,*) '======================================================================'
      if (irank_g .eq. 0 ) write(*,*) 'overlap for input slater spin dn', z_norm_dn
      z_norm_dn = (1.d0/z_norm_dn)**(1.d0/dble(n_part_2))
      if (irank_g .eq. 0 ) write(*,*) 'renormalized factor for input slater z_norm spin dn', z_norm_dn
      if (irank_g .eq. 0 ) write(*,*) '======================================================================'

      do nw = 1, N_wlk
          phi_0_r(1,nw)%U(:,:)=p2_tmp(:,:)*z_norm_up
          phi_0_r(2,nw)%U(:,:)=p3_tmp(:,:)*z_norm_dn
      enddo
      !! here we assume single SD as trial wave function
      WF_L(1,1)%P(:,:)=p0_tmp(:,:)
      WF_R(1,1)%P(:,:)=p0_tmp(:,:)
      phi_0_l(1,1)%U(:,:)=p0_tmp(:,:)

      WF_L(2,1)%P(:,:)=p1_tmp(:,:)
      WF_R(2,1)%P(:,:)=p1_tmp(:,:)
      phi_0_l(2,1)%U(:,:)=p1_tmp(:,:)

      if (irank_g .eq. 0 ) then
          deallocate(phi0_up_out, phi0_dn_out)
      endif
      deallocate(p0_tmp, p1_tmp, p2_tmp, p3_tmp)
      deallocate(ipiv_up, ipiv_dn)
      deallocate(smat_up, smat_dn)

   end subroutine trial_in_hdf5
#endif

   subroutine seed_vec_in(file_tg)
      
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

end module BRW_init_mod
