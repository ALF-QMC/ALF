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
!
!> @brief
!> This module defines the Obser_Vec and Obser_Latt types and provides
!> routine to initialize them and to print out the bins
!
!--------------------------------------------------------------------

module Observables

#if !defined HDF5 && !defined OBS_LEGACY
#define OBS_LEGACY 1
#endif
   use runtime_error_mod
   use Lattices_v3, only: Unit_cell, Lattice
   use iso_fortran_env, only: output_unit, error_unit

   type :: Obser_Vec
!>  Data structure for
!>  < O_n >  n : =1, size(Obs,1)
      !private
      integer                     :: N                    ! Number of measurements
      real(Kind=kind(0.d0)) :: Ave_Sign             ! Averarge sign
      complex(Kind=kind(0.d0)), pointer :: Obs_vec(:)  ! Vector of observables
      character(len=64) :: File_Vec                      ! Name of file in which the bins will be written out
      character(len=64) :: analysis_mode                 ! How to analyze the observable
      character(len=64), allocatable :: description(:)   ! Optional short description
   contains
      !procedure :: make        => Obser_vec_make
      procedure :: init => Obser_vec_init
      procedure :: print_bin => print_bin_vec
      procedure :: measure => Obser_vec_measure
   end type Obser_Vec

   type :: Obser_Latt
!>  Data structure for
!>  < O^{dagger}(i,tau)_n O(j,0)_m>  - < O_n> <O_m>
!>  where it is assumed that translation symmetry as specified by the lattice Latt is present.
!>  Obs_Latt(i-j,tau,n,m) = < O^{dagger}(i,tau)_n O(j,0)_m>
!>  Obs_Latt0(n) = < O_n>
!>  For equal   time correlation functions, tau runs from 1,1
!>  For unequal time correlation functions, tau runs from 1,Ltrot+1
      integer            :: N                                    ! Number of measurements
      real(Kind=kind(0.d0)) :: Ave_Sign                    ! Averarge sign
      complex(Kind=kind(0.d0)), pointer :: Obs_Latt(:, :, :, :) ! i-j, tau, norb, norb
      complex(Kind=kind(0.d0)), pointer :: Obs_Latt0(:)       ! norb
      character(len=64) :: File_Latt                            ! Name of file in which the bins will be written out
      type(Lattice), pointer :: Latt                      ! Pointer to Bravais lattice
      type(Unit_cell), pointer :: Latt_unit                 ! Pointer to unit cell
      real(Kind=kind(0.d0))   :: dtau                      ! Imaginary time step
      character(len=:), allocatable  :: Channel    ! Type of observable. Possible values:
      ! - T0  : zero temperature
      ! - P   : finite temperature particle
      ! - P_PH: finite temperature particle   with particle-hole  symmetry
      ! - PH  : finite temperature particle-hole
      ! - PP  : finite temperature particle-particle
   contains
      !procedure :: make        => Obser_latt_make
      procedure :: init => Obser_latt_init
      procedure :: print_bin => print_bin_latt
   end type Obser_Latt

contains

   subroutine Obser_Latt_make(Obs, Nt, Filename, Latt, Latt_unit, Channel, dtau)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Create lattice type observable
!>
!> @details
!> Create lattice type observable. Be aware that Latt and Latt_unit don't get copied
!> but linked, meaning changing them after making the observable still affects the
!> observable.
!>
!> @param [INOUT] Obs, Type(Obser_Latt)
!> \verbatim
!>  Observable to define
!> \endverbatim
!> @param [IN] Nt, Integer
!> \verbatim
!>  Number of imaginary time points, set to 1 for equal time correlators.
!> \endverbatim
!> @param [IN] Filename, Character(len=64)
!> \verbatim
!>  Name of file in which the bins will be written out.
!> \endverbatim
!> @param [IN] Latt, Type(Lattice)
!> \verbatim
!>  Bravais lattice. Only gets linked, needs attribute target or pointer.
!> \endverbatim
!> @param [IN] Latt_unit, Type(Unit_cell)
!> \verbatim
!>  Unit cell. Only gets linked, needs attribute target or pointer.
!> \endverbatim
!> @param [IN] Channel, Character(len=*)
!> \verbatim
!>  MaxEnt channel. Only relevant for time displaced observables.
!> \endverbatim
!> @param [IN] dtau, Real(Kind=Kind(0.d0))
!> \verbatim
!>  Imaginary time step. Only relevant for time displaced observables.
!> \endverbatim
!-------------------------------------------------------------------
      implicit none
      type(Obser_Latt), intent(INOUT)      :: Obs
      integer, intent(IN)         :: Nt
      character(len=64), intent(IN)         :: Filename
      type(Lattice), intent(IN), target :: Latt
      type(Unit_cell), intent(IN), target :: Latt_unit
      character(len=*), intent(IN)         :: Channel
      real(Kind=kind(0.d0)), intent(IN)    :: dtau
      allocate (Obs%Obs_Latt(Latt%N, Nt, Latt_unit%Norb, Latt_unit%Norb))
      allocate (Obs%Obs_Latt0(Latt_unit%Norb))
      Obs%File_Latt = Filename
      Obs%Latt => Latt
      Obs%Latt_unit => Latt_unit
      Obs%Channel = Channel
      Obs%dtau = dtau
   end subroutine Obser_Latt_make
!--------------------------------------------------------------------

   subroutine Obser_Latt_Init(Obs)
      implicit none
      class(Obser_Latt), intent(INOUT) :: Obs
      Obs%Obs_Latt = cmplx(0.d0, 0.d0, kind(0.d0))
      Obs%Obs_Latt0 = cmplx(0.d0, 0.d0, kind(0.d0))
      Obs%N = 0
      Obs%Ave_Sign = 0.d0
   end subroutine Obser_Latt_Init

!--------------------------------------------------------------------

   subroutine Obser_Vec_make(Obs, N, Filename, analysis_mode, description)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Create scalar type observable
!>
!> @param [INOUT] Obs, Type(Obser_vec)
!> \verbatim
!>  Observable to define
!> \endverbatim
!> @param [IN] N, Integer
!> \verbatim
!>  Number of scalars in this observable.
!> \endverbatim
!> @param [IN] Filename, Character(len=64)
!> \verbatim
!>  Name of file in which the bins will be written out.
!> \endverbatim
!> @param [IN] analysis_mode, Character(len=64), optional
!> \verbatim
!>  How to analyze the observable.
!> \endverbatim
!> @param [IN] description(:), Character(len=64), optional
!> \verbatim
!>  Optional array to describe observable.
!> \endverbatim
!-------------------------------------------------------------------
      implicit none
      type(Obser_vec), intent(INOUT) :: Obs
      integer, intent(IN)             :: N
      character(len=64), intent(IN)  :: Filename
      character(len=64), intent(IN), optional :: analysis_mode
      character(len=64), intent(IN), optional :: description(:)

      allocate (Obs%Obs_vec(N))
      Obs%File_Vec = Filename
      if (present(analysis_mode)) then
         Obs%analysis_mode = analysis_mode
      else
         Obs%analysis_mode = 'identity'
      end if
      if (present(description)) then
         allocate (Obs%description(size(description, 1)))
         Obs%description = description
      end if
   end subroutine Obser_Vec_make
!--------------------------------------------------------------------

   subroutine Obser_Vec_Init(Obs)
      implicit none
      class(Obser_vec), intent(INOUT) :: Obs
      Obs%Obs_vec = cmplx(0.d0, 0.d0, kind(0.d0))
      Obs%N = 0
      Obs%Ave_Sign = 0.d0
   end subroutine Obser_Vec_Init

!--------------------------------------------------------------------

   subroutine Obser_vec_measure(obs, value, Phase)
      implicit none

      class(Obser_vec), intent(Inout) :: Obs
      complex(Kind=kind(0.d0)), intent(In)    :: value(:)  ! Vector of observables
      complex(Kind=kind(0.d0)), intent(IN), optional    :: Phase
      !Local
      complex(Kind=kind(0.d0)) :: ZP, ZS

      obs%N = obs%N + 1

      if (present(Phase)) then
         ZP = PHASE/real(Phase, kind(0.d0))
         ZS = real(Phase, kind(0.d0))/abs(real(Phase, kind(0.d0)))

         obs%Ave_sign = obs%Ave_sign + real(ZS, kind(0.d0))
         obs%obs_vec = obs%obs_vec + value*ZS*ZP
      else
         obs%Ave_sign = obs%Ave_sign + 1.d0
         obs%obs_vec = obs%obs_vec + value
      end if

   end subroutine Obser_vec_measure

!--------------------------------------------------------------------

   subroutine Print_bin_Latt(Obs, Group_Comm)
      use Lattices_v3
#if defined HDF5
      use hdf5
      use alf_hdf5
#endif
#if defined MPI
      use mpi
#endif
      implicit none

      class(Obser_Latt), intent(Inout)   :: Obs
      integer, intent(In)      :: Group_Comm

      ! Local
      integer :: Ns, Nt, no, no1, I, Ntau
      complex(Kind=kind(0.d0)), allocatable, target :: Tmp(:, :, :, :)
      real(Kind=kind(0.d0))              :: x_p(2)
      complex(Kind=kind(0.d0))              :: Sign_bin
      character(len=64) :: File_pr, File_suff, File_aux, tmp_str
      logical            :: File_exists
#ifdef HDF5
      character(len=7), parameter  :: File_h5 = "data.h5"
      character(len=64)            :: filename, groupname, obs_dsetname, bak_dsetname, sgn_dsetname
      integer(HID_T)                :: file_id, group_id
      logical                       :: link_exists
      integer                       :: hdferr
      integer(HSIZE_T), allocatable :: dims(:)
      type(c_ptr)                   :: dat_ptr
      real(Kind=kind(0.d0)), target :: sgn
#endif
#ifdef MPI
      complex(Kind=kind(0.d0)), allocatable :: Tmp1(:)
      complex(Kind=kind(0.d0)) :: Z
      real(Kind=kind(0.d0)) :: X
      integer         :: Ierr, Isize, Irank
      integer         :: irank_g, isize_g, igroup

      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

      Ns = size(Obs%Obs_Latt, 1)
      Ntau = size(Obs%Obs_Latt, 2)
      if (.not. (Obs%Latt%N == Ns)) then
         write (error_unit, *) 'Error in Print_bin_Latt'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      if (Ntau == 1) then
         File_suff = "_eq"
      else
         File_suff = "_tau"
      end if
      write (File_pr, '(A,A)') trim(Obs%File_Latt), trim(File_suff)
#if defined HDF5
      groupname = File_pr
      filename = File_h5
#endif
      allocate (Tmp(Ns, Ntau, Obs%Latt_unit%Norb, Obs%Latt_unit%Norb))
      Obs%Obs_Latt = Obs%Obs_Latt/dble(Obs%N)
      Obs%Obs_Latt0 = Obs%Obs_Latt0/dble(Obs%N*Ns*Ntau)
      Obs%Ave_sign = Obs%Ave_Sign/dble(Obs%N)

#if defined(MPI)
      I = Obs%Latt%N*Ntau*Obs%Latt_unit%Norb*Obs%Latt_unit%Norb
      Tmp = cmplx(0.d0, 0.d0, kind(0.d0))
      call MPI_REDUCE(Obs%Obs_Latt, Tmp, I, MPI_COMPLEX16, MPI_SUM, 0, Group_Comm, IERR)
      Obs%Obs_Latt = Tmp/dble(ISIZE_g)

      I = 1
      X = 0.d0
      call MPI_REDUCE(Obs%Ave_sign, X, I, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      Obs%Ave_sign = X/dble(ISIZE_g)

      I = Obs%Latt_unit%Norb
      allocate (Tmp1(Obs%Latt_unit%Norb))
      Tmp1 = cmplx(0.d0, 0.d0, kind(0.d0))
      call MPI_REDUCE(Obs%Obs_Latt0, Tmp1, I, MPI_COMPLEX16, MPI_SUM, 0, Group_Comm, IERR)
      Obs%Obs_Latt0 = Tmp1/dble(ISIZE_g)
      deallocate (Tmp1)

      if (Irank_g == 0) then
#endif
#if defined TEMPERING
         write (File_pr, '(A,I0,A,A,A)') "Temp_", igroup, "/", trim(Obs%File_Latt), trim(File_suff)
#if defined HDF5
         write (filename, '(A,I0,A,A)') "Temp_", igroup, "/", trim(File_h5)
#endif
#endif

         do nt = 1, Ntau
            do no = 1, Obs%Latt_unit%Norb
               do no1 = 1, Obs%Latt_unit%Norb
                  call Fourier_R_to_K(Obs%Obs_Latt(:, nt, no, no1), Tmp(:, nt, no, no1), Obs%Latt)
               end do
            end do
         end do

#if defined OBS_LEGACY
         write (File_aux, '(A,A)') trim(File_pr), "_info"
         inquire (file=File_aux, exist=File_exists)
         if (.not. File_exists) then
11          format(A20, ': ', A)
12          format(A20, ': ', I10)
13          format(A20, ': ', *(E26.17E3))
14          format(A20, ': ')
15          format((E26.17E3))
            open (10, file=File_aux, status='new')
            write (tmp_str, '(A, A)') trim(Obs%File_Latt), trim(File_suff)
            write (10, 11) 'Observable', trim(tmp_str)
            write (10, 11) 'Channel', trim(Obs%Channel)
            write (10, 12) 'Ntau', Ntau
            write (10, 13) 'dtau', Obs%dtau
            write (10, '(A)') '       ====== Bravais Lattice ======'
            write (10, 12) 'Unit cells', Obs%Latt%N
            write (10, 13) 'L1', Obs%Latt%L1_p
            write (10, 13) 'L2', Obs%Latt%L2_p
            write (10, 13) 'a1', Obs%Latt%a1_p
            write (10, 13) 'a2', Obs%Latt%a2_p
            write (10, '(A)') '       ========= Unit cell ========='
            write (10, 12) 'Coordination number', Obs%Latt_unit%N_coord
            write (10, 12) 'Number of orbitals', Obs%Latt_unit%Norb
            write (10, 12) 'Ndim', size(Obs%Latt_unit%Orb_pos_p, 2)
            do no = 1, Obs%Latt_unit%Norb
               write (tmp_str, '("Orbital ",I0)') no
               write (10, 14, advance='no') trim(tmp_str)
               do i = 1, size(Obs%Latt_unit%Orb_pos_p, 2)
                  write (10, 15, advance='no') Obs%Latt_unit%Orb_pos_p(no, i)
               end do
               write (10, *)
            end do
            close (10)
         end if

         open (Unit=10, File=File_pr, status="unknown", position="append")
         if (Ntau == 1) then
            write (10, '(E25.17E3, 2(I11))') Obs%Ave_sign, Obs%Latt_unit%Norb, Obs%Latt%N
         else
            write (10, '(E25.17E3, 3(I11), E26.17E3)') Obs%Ave_sign, Obs%Latt_unit%Norb, Obs%Latt%N, Ntau, Obs%dtau
         end if
         do no = 1, Obs%Latt_unit%Norb
            write (10, '("(", E25.17E3, ",", E25.17E3, ")")') Obs%Obs_Latt0(no)
         end do
         do I = 1, Obs%Latt%N
            x_p = dble(Obs%Latt%listk(i, 1))*Obs%Latt%b1_p + dble(Obs%Latt%listk(i, 2))*Obs%Latt%b2_p
            write (10, '(E25.17E3, 1x, E25.17E3)') X_p(1), X_p(2)
            do nt = 1, Ntau
               do no = 1, Obs%Latt_unit%Norb
                  do no1 = 1, Obs%Latt_unit%Norb
                     write (10, '("(", E25.17E3, ",", E25.17E3, ")")') tmp(I, nt, no, no1)
                  end do
               end do
            end do
         end do
         close (10)
#endif

#if defined HDF5
         write (obs_dsetname, '(A,A,A)') trim(groupname), "/obser"
         write (bak_dsetname, '(A,A,A)') trim(groupname), "/back"
         write (sgn_dsetname, '(A,A,A)') trim(groupname), "/sign"

         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdferr)

         !Writes the lattice to HDF5 file if it doesn't already exist
         call write_latt(file_id, Obs%Latt, Obs%Latt_Unit)

         call h5lexists_f(file_id, groupname, link_exists, hdferr)
         if (.not. link_exists) then
            !Create Group for observable
            call h5gcreate_f(file_id, groupname, group_id, hdferr)
            call write_attribute(group_id, '.', "dtau", Obs%dtau, hdferr)
            call write_attribute(group_id, '.', "Channel", Obs%Channel, hdferr)
            call write_latt(group_id, Obs%Latt, Obs%Latt_Unit)
            call h5gclose_f(group_id, hdferr)

            !Create Dataset for data
            allocate (dims(6))
            dims = [2, Obs%Latt%N, Ntau, Obs%Latt_unit%Norb, Obs%Latt_unit%Norb, 0]
            call init_dset(file_id, obs_dsetname, dims, .true.)
            deallocate (dims)

            !Create Dataset for background
            allocate (dims(3))
            dims = [2, Obs%Latt_unit%Norb, 0]
            call init_dset(file_id, bak_dsetname, dims, .true.)
            deallocate (dims)

            !Create Dataset for sign
            allocate (dims(1))
            dims = [0]
            call init_dset(file_id, sgn_dsetname, dims, .false.)
            deallocate (dims)
         end if

         !Write data
         dat_ptr = c_loc(tmp(1, 1, 1, 1))
         call append_dat(file_id, obs_dsetname, dat_ptr)

         !Write background
         dat_ptr = c_loc(Obs%Obs_Latt0(1))
         call append_dat(file_id, bak_dsetname, dat_ptr)

         !Write sign
         sgn = Obs%Ave_sign
         dat_ptr = c_loc(sgn)
         call append_dat(file_id, sgn_dsetname, dat_ptr)

         call h5fclose_f(file_id, hdferr)
#endif
#if defined MPI
      end if
#endif

      deallocate (Tmp)

   end subroutine Print_bin_Latt

!--------------------------------------------------------------------

   subroutine Print_bin_Vec(Obs, Group_Comm)
#if defined MPI
      use mpi
#endif
#if defined HDF5
      use hdf5
      use alf_hdf5
#endif
      implicit none

      class(Obser_vec), intent(Inout) :: Obs
      integer, intent(IN)  :: Group_Comm

      ! Local
      integer :: I
      character(len=64) :: File_pr, File_suff, File_aux
      logical            :: File_exists

#if defined HDF5
      character(len=7), parameter  :: File_h5 = "data.h5"
      character(len=64)            :: filename, groupname, obs_dsetname, sgn_dsetname
      integer(HID_T)                :: file_id, group_id
      logical                       :: link_exists
      integer                       :: hdferr
      integer(HSIZE_T), allocatable :: dims(:)
      type(c_ptr)                   :: dat_ptr
      real(Kind=kind(0.d0)), target :: sgn
#endif
#if defined MPI
      integer        :: Ierr, Isize, Irank, No
      integer        :: irank_g, isize_g, igroup
      complex(Kind=kind(0.d0)), allocatable :: Tmp(:)
      real(Kind=kind(0.d0)) :: X

      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif
      Obs%Obs_vec = Obs%Obs_vec/dble(Obs%N)
      Obs%Ave_sign = Obs%Ave_sign/dble(Obs%N)
      write (File_pr, '(A,A)') trim(Obs%File_Vec), "_scal"
#if defined HDF5
      groupname = File_pr
      filename = File_h5
#endif

#if defined MPI
      No = size(Obs%Obs_vec, 1)
      allocate (Tmp(No))
      Tmp = cmplx(0.d0, 0.d0, kind(0.d0))
      call MPI_REDUCE(Obs%Obs_vec, Tmp, No, MPI_COMPLEX16, MPI_SUM, 0, Group_Comm, IERR)
      Obs%Obs_vec = Tmp/dble(ISIZE_g)
      deallocate (Tmp)

      I = 1
      X = 0.d0
      call MPI_REDUCE(Obs%Ave_sign, X, I, MPI_REAL8, MPI_SUM, 0, Group_comm, IERR)
      Obs%Ave_sign = X/dble(ISIZE_g)

      if (Irank_g == 0) then
#endif

#if defined TEMPERING
         write (File_pr, '(A,I0,A,A,A)') "Temp_", igroup, "/", trim(Obs%File_Vec), "_scal"
#if defined HDF5
         write (filename, '(A,I0,A,A)') "Temp_", igroup, "/", trim(File_h5)
#endif
#endif

#if defined OBS_LEGACY
         write (File_aux, '(A,A)') trim(File_pr), "_info"
         inquire (file=File_aux, exist=File_exists)
         if (.not. File_exists) then
            open (10, file=File_aux, status='new')
            write (10, '(A)') '====== Analysis Mode ======'
            write (10, '(A)') trim(Obs%analysis_mode)
            if (allocated(Obs%description)) then
               write (10, '(A)') '====== Description ======'
               do i = 1, size(Obs%description, 1)
                  write (10, '(A)') trim(Obs%description(i))
               end do
            end if
            close (10)
         end if
         open (Unit=10, File=File_pr, status="unknown", position="append")
         !WRITE(10,*) size(Obs%Obs_vec,1)+1, (Obs%Obs_vec(I), I=1,size(Obs%Obs_vec,1)), Obs%Ave_sign
         write (10, '(I10)', advance='no') size(Obs%Obs_vec, 1) + 1
         do I = 1, size(Obs%Obs_vec, 1)
            write (10, '(" (",E25.17E3,",",E25.17E3,")")', advance='no') Obs%Obs_vec(I)
         end do
         write (10, '(E26.17E3)') Obs%Ave_sign
         close (10)
#endif

#if defined HDF5
         write (obs_dsetname, '(A,A,A)') trim(groupname), "/obser"
         write (sgn_dsetname, '(A,A,A)') trim(groupname), "/sign"

         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdferr)

         !Check if observable already exists in hdf5 file
         call h5lexists_f(file_id, groupname, link_exists, hdferr)

         if (.not. link_exists) then
            !Create Group for observable and write auxiliary info
            call h5gcreate_f(file_id, groupname, group_id, hdferr)
            call write_attribute(group_id, '.', "analysis_mode", Obs%analysis_mode, hdferr)
            if (allocated(Obs%description)) then
               call write_comment(group_id, '.', "description", Obs%description, hdferr)
            end if
            call h5gclose_f(group_id, hdferr)

            !Create Dataset for data
            allocate (dims(3))
            dims = [2, size(Obs%Obs_vec, 1), 0]
            call init_dset(file_id, obs_dsetname, dims, .true.)
            deallocate (dims)

            !Create Dataset for sign
            allocate (dims(1))
            dims = [0]
            call init_dset(file_id, sgn_dsetname, dims, .false.)
            deallocate (dims)
         end if

         !Write data
         dat_ptr = c_loc(Obs%Obs_vec(1))
         call append_dat(file_id, obs_dsetname, dat_ptr)

         !Write sign
         sgn = Obs%Ave_sign
         dat_ptr = c_loc(sgn)
         call append_dat(file_id, sgn_dsetname, dat_ptr)

         call h5fclose_f(file_id, hdferr)
#endif
#if defined MPI
      end if
#endif

   end subroutine Print_bin_Vec

end module Observables
