!  Copyright (C) 2019 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

   module ana_mod
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Collection of routines for postprocessing the Monte Carlo bins
!
!--------------------------------------------------------------------
      use iso_fortran_env, only: output_unit, error_unit
      use Files_mod
      Use Errors
      Use MyMats
      Use Matrix
      Use Lattices_v3, only: Unit_cell, Lattice, Make_lattice, Fourier_K_to_R, Inv_K
      Use Predefined_Lattices, only: Predefined_Latt
#ifdef HDF5
      use hdf5
      use h5lt
      Use alf_hdf5
#endif

   contains

   Subroutine read_vec(file, sgn, bins, analysis_mode)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in bins of scalar observables from file
!>
!> @param [IN] file Character(len=*)
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!> @param [OUT] sgn Real(:)
!> \verbatimam
!>  Sign of bins
!> \endverbatim
!> @param [OUT] bins Complex(:,:)
!> \verbatim
!>  Monte Carlo bins
!> \endverbatim
!> @param [OUT] analysis_mode Character(len=64)
!> \verbatim
!>  How to analyze the observable
!> \endverbatim
!-------------------------------------------------------------------
      Implicit none
      Character (len=*), intent(in) :: file
      Real    (Kind=Kind(0.d0)), allocatable, intent(out) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer, intent(out) :: bins(:,:)
      Character (len=64), intent(out) :: analysis_mode

      Integer :: N, N1, I, Nobs, Nbins, stat
      Real    (Kind=Kind(0.d0)) :: X
      Complex (Kind=Kind(0.d0)), Allocatable  :: tmp(:)
      Character (len=64) :: file_aux
      logical :: file_exists

      write(file_aux, '(A,A)') trim(file), "_info"
      inquire(file=file_aux, exist=file_exists)
      if(file_exists) then
        open(Unit=10, File=file_aux, status="old", action='read')
        read(10, *)
        read(10, '(A)') analysis_mode
        close(10)
      else
        analysis_mode = 'identity'
      endif

      open(Unit=10, File=file, status="old", action='read')
      read(10,*) NOBS
      NOBS = NOBS-1
      allocate(tmp(NOBS))
      rewind(10)

      Nbins = 0
      do
         read(10, *, iostat=stat) N1, (tmp(I), I=1,size(Tmp,1)),X
         if (stat /= 0) exit
         Nbins = Nbins + 1
      enddo
      rewind(10)

      allocate(bins(NOBS,Nbins))
      allocate(sgn(Nbins))

      do N = 1,Nbins
         read(10,*) N1, (tmp(I), I=1,size(Tmp,1)),X

         bins(:,N) = Tmp(:)
         sgn(N)  = X
      enddo
      close(10)
      deallocate(tmp)
   End Subroutine read_vec

!==============================================================================

#ifdef HDF5
   Subroutine read_vec_hdf5(filename, groupname, sgn, bins, analysis_mode)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in bins of scalar observables from HDF5 file
!>
!> @param [IN] filename Character(len=*)
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!> @param [IN] groupname Character(len=*)
!> \verbatim
!>  Name of observable that gets read in
!> \endverbatim
!> @param [OUT] sgn Real(:)
!> \verbatimam
!>  Sign of bins
!> \endverbatim
!> @param [OUT] bins Complex(:,:)
!> \verbatim
!>  Monte Carlo bins
!> \endverbatim
!> @param [OUT] analysis_mode Character(len=64)
!> \verbatim
!>  How to analyze the observable
!> \endverbatim
!-------------------------------------------------------------------
      Implicit none
      Character (len=*), intent(in) :: filename
      Character (len=*), intent(in) :: groupname
      Real    (Kind=Kind(0.d0)), allocatable, intent(out) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer, intent(out) :: bins(:,:)
      Character (len=64), intent(out) :: analysis_mode

      Integer :: Nobs, Nbins

      Character (len=64) :: file, obs_dsetname, sgn_dsetname
      INTEGER                       :: rank, hdferr
      INTEGER(HSIZE_T)              :: mem_dims(1)
      INTEGER(HSIZE_T), allocatable :: dims(:), maxdims(:)
      INTEGER(HID_T)                :: file_id, dset_id, dataspace, memspace
      TYPE(C_PTR)                   :: dat_ptr

      file = 'data.h5'
      write(obs_dsetname,'(2A)') trim(groupname), "/obser"
      write(sgn_dsetname,'(2A)') trim(groupname), "/sign"

      CALL h5open_f(hdferr)

      CALL h5fopen_f (File, H5F_ACC_RDONLY_F, file_id, hdferr)
      call read_attribute(file_id, groupname, "analysis_mode", analysis_mode, hdferr)
      !Open the  dataset.
      CALL h5dopen_f(file_id, obs_dsetname, dset_id, hdferr)

      !Get dataset's dataspace handle.
      CALL h5dget_space_f(dset_id, dataspace, hdferr)

      !Get dataspace's rank.
      CALL h5sget_simple_extent_ndims_f(dataspace, rank, hdferr)
      allocate( dims(rank), maxdims(rank) )

      !Get dataspace's dimensions.
      CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, hdferr)

      Nobs  = int( dims(2) )
      Nbins = int( dims(3) )
      deallocate( dims )

      Allocate ( Bins(Nobs,Nbins), sgn(Nbins))

      dat_ptr = C_LOC(bins(1,1))
      CALL h5dopen_f(file_id, obs_dsetname, dset_id, hdferr)
      CALL h5dget_space_f(dset_id, dataspace, hdferr)
      CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr, &
                  mem_space_id = dataspace, file_space_id = dataspace)

      allocate( dims(1) )
      dims(1) = Nbins
      CALL h5dopen_f(file_id, sgn_dsetname, dset_id, hdferr)
      CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, sgn, dims, hdferr)
      deallocate( dims )

   End Subroutine read_vec_hdf5
#endif

!==============================================================================


   Subroutine read_latt(file, sgn, bins, bins0, Latt, Latt_unit, dtau, Channel)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in bins of lattice-type observables (both equal time and timedisplaced) from file
!>
!> @param [IN] file Character(len=64)
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!> @param [OUT] sgn Real(:)
!> \verbatimam
!>  Sign of bins
!> \endverbatim
!> @param [OUT] bins Complex(:,:,:,:,:)
!> \verbatim
!>  Monte Carlo bins of correlation
!> @param [OUT] bins0 Complex(:,:)
!> \verbatim
!>  Monte Carlo bins of background
!> \endverbatim
!> @param [OUT] Latt Type(Lattice)
!> \verbatim
!>  Bravais lattice
!> \endverbatim
!> @param [OUT] Latt_unit Type(Unit_cell)
!> \verbatim
!>  Unit cell
!> \endverbatim
!> @param [OUT] dtau Real
!> \verbatim
!>  Size of imaginary time step
!> \endverbatim
!> @param [OUT] file Character(len=2)
!> \verbatim
!>  MaxEnt Channel. Relevant for timedisplaced correlation.
!> \endverbatim
!-------------------------------------------------------------------
      Implicit none
      Character (len=*), intent(in) :: file
      Real    (Kind=Kind(0.d0)), allocatable, intent(out) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer    , intent(out) :: bins(:,:,:,:,:)
      Complex (Kind=Kind(0.d0)), pointer    , intent(out) :: bins0(:,:)
      Type (Lattice)                        , intent(out) :: Latt
      Type (Unit_cell)                      , intent(out) :: Latt_unit
      Real    (Kind=Kind(0.d0))             , intent(out) :: dtau
      Character (len=:), allocatable        , intent(out) :: Channel
      

      Character (len=64) :: file_aux, str_temp1,  str_temp2
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      Integer :: no, no1, n, nt, nb, Ntau, Ndim, Nbins, stat, Ndim_unit
      Real(Kind=Kind(0.d0)) :: X
      Real(Kind=Kind(0.d0)), allocatable :: Xk_p(:,:)
      Real(Kind=Kind(0.d0)) :: x_p(2), a1_p(2), a2_p(2), L1_p(2), L2_p(2)
      logical            :: file_exists

      Integer             :: L1, L2
      Character (len=64)  :: Model, Lattice_type
      NAMELIST /VAR_Lattice/ L1, L2, Lattice_type, Model


      write(file_aux, '(A,A)') trim(file), "_info"
      inquire(file=file_aux, exist=file_exists)
      if(file_exists) then
        open(Unit=10, File=file_aux, status="old", action='read')
        11 format(A22, A)
        12 format(A22, I10)
        13 format(A22, *(E26.17E3))
        read(10, *)
        read(10, 11) str_temp1, str_temp2
        Channel  = trim(str_temp2)
        read(10, 12) str_temp1, Ntau
        read(10, 13) str_temp1, dtau
        read(10, *)
        read(10, 12) str_temp1, Latt%N
        read(10, 13) str_temp1, L1_p
        read(10, 13) str_temp1, L2_p
        read(10, 13) str_temp1, a1_p
        read(10, 13) str_temp1, a2_p
        read(10, *)
        read(10, 12) str_temp1, Latt_unit%N_coord
        read(10, 12) str_temp1, Latt_unit%Norb
        read(10, 12) str_temp1, Ndim_unit
        allocate(Latt_unit%Orb_pos_p(Latt_unit%Norb, Ndim_unit))
        do no = 1, Latt_unit%Norb
          read(10, 13) str_temp1, Latt_unit%Orb_pos_p(no,:)
        enddo
        close(10)
        Call Make_Lattice(L1_p, L2_p, a1_p, a2_p, Latt)
        Ndim = Latt%N*Latt_Unit%Norb
      else
        Channel = '--'
        open(Unit=10, File='parameters', status="old", action='read')
        read(10, NML=VAR_lattice)
        close(10)
        Call Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)
        open(Unit=10, File=file, status="old", action='read')
        Read(10, *, iostat=stat) X, Latt_unit%Norb, Latt%N, Ntau, dtau
        if (stat /= 0) then
           rewind(10)
           Ntau = 1
           dtau = -1.d0
           Read(10, *) X, Latt_unit%Norb, Latt%N
        endif
        close(10)
      endif

      ! Determine the number of bins.
      open(Unit=10, File=file, status="old", action='read')
      nbins = 0
      do
         if(Ntau == 1) then
            read(10, *, iostat=stat) X, Latt_unit%Norb, Latt%N
         else
            read(10, *, iostat=stat) X, Latt_unit%Norb, Latt%N, Ntau, dtau
         endif
         if (stat /= 0) exit
         Do no = 1, Latt_unit%Norb
            Read(10,*)
         enddo
         do n = 1, Latt%N
            Read(10,*)
            do nt = 1, Ntau
               do no = 1, Latt_unit%Norb
                  do no1 = 1, Latt_unit%Norb
                     read(10,*)
                  enddo
               enddo
            enddo
         enddo
         nbins = nbins + 1
      enddo
      rewind(10)

      ! Allocate  space
      Allocate(bins(Latt%N, Ntau, Latt_unit%Norb, Latt_unit%Norb, Nbins))
      Allocate(sgn(Nbins), Xk_p(2, Latt%N), bins0(Latt_unit%Norb, Nbins))

      do nb = 1, nbins
         if(Ntau == 1) then
            read(10, *) sgn(nb), Latt_unit%Norb, Latt%N
         else
            read(10, *) sgn(nb), Latt_unit%Norb, Latt%N, Ntau, dtau
         endif
         Do no = 1, Latt_unit%Norb
            Read(10,*) bins0(no,nb)
         Enddo
         do n = 1, Latt%N
            Read(10,*) Xk_p(1,n), Xk_p(2,n)
            do nt = 1, Ntau
               do no = 1, Latt_unit%Norb
                  do no1 = 1, Latt_unit%Norb
                     read(10,*) bins(n,nt,no,no1,nb)
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(10)

      do n = 1, Latt%N
         x_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
         x   = (x_p(1)-Xk_p(1,n))**2 + (x_p(2)-Xk_p(2,n))**2
         if ( x > 0.00001 ) then
            Write(error_unit,*) "Error in read_latt: momenta do not fit", x, x_p, Xk_p(1,n)
            error stop
         endif
      enddo

   End Subroutine read_latt

!==============================================================================

#ifdef HDF5
Subroutine read_latt_hdf5(filename, name, sgn, bins, bins0, Latt, Latt_unit, dtau, Channel)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in bins of lattice-type observables (both equal time and timedisplaced)
!> from HDF5 file data.h5
!>
!> @param [IN] name Character(len=64)
!> \verbatim
!>  Name of observable that gets read in
!> \endverbatim
!> @param [OUT] sgn Real(:)
!> \verbatimam
!>  Sign of bins
!> \endverbatim
!> @param [OUT] bins Complex(:,:,:,:,:)
!> \verbatim
!>  Monte Carlo bins of correlation
!> @param [OUT] bins0 Complex(:,:)
!> \verbatim
!>  Monte Carlo bins of background
!> \endverbatim
!> @param [OUT] Latt Type(Lattice)
!> \verbatim
!>  Bravais lattice
!> \endverbatim
!> @param [OUT] Latt_unit Type(Unit_cell)
!> \verbatim
!>  Unit cell
!> \endverbatim
!> @param [OUT] dtau Real
!> \verbatim
!>  Size of imaginary time step
!> \endverbatim
!> @param [OUT] file Character(len=2)
!> \verbatim
!>  MaxEnt Channel. Relevant for timedisplaced correlation.
!> \endverbatim
!-------------------------------------------------------------------
      Implicit none
      Character (len=*), intent(in) :: filename
      Character (len=*), intent(in) :: name
      Real    (Kind=Kind(0.d0)), allocatable, intent(out) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer    , intent(out) :: Bins(:,:,:,:,:)
      Complex (Kind=Kind(0.d0)), pointer    , intent(out) :: Bins0(:,:)
      Type (Lattice)                        , intent(out) :: Latt
      Type (Unit_cell)                      , intent(out) :: Latt_unit
      Real    (Kind=Kind(0.d0))             , intent(out) :: dtau
      Character (len=*)                     , intent(out) :: Channel

      Integer    :: Nbins, Norb

      Character (len=64) :: obs_dsetname, bak_dsetname, sgn_dsetname, par_dsetname, attr_name
      INTEGER                       :: ierr, rank, Nunit, Ntau, Ndim, no
      INTEGER(HSIZE_T), allocatable :: dims(:), maxdims(:)
      INTEGER(HID_T)                :: file_id, dset_id, grp_id, dataspace
      TYPE(C_PTR)                   :: dat_ptr
      Real (Kind=Kind(0.d0))        :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)

      write(obs_dsetname,'(2A)') trim(name), "/obser"
      write(bak_dsetname,'(2A)') trim(name), "/back"
      write(sgn_dsetname,'(2A)') trim(name), "/sign"

      !Initialize HDF5
      CALL h5open_f(ierr)

      !Open HDF5 file
      CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, ierr)

      !Open the dataset.
      CALL h5dopen_f(file_id, obs_dsetname, dset_id, ierr)

      !Get dataset's dataspace handle.
      CALL h5dget_space_f(dset_id, dataspace, ierr)

      !Get dataspace's rank.
      CALL h5sget_simple_extent_ndims_f(dataspace, rank, ierr)
      allocate( dims(rank), maxdims(rank) )

      !Get dataspace's dimensions.
      CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, ierr)
      if ( rank == 5 ) then
         Nunit = int( dims(2) )
         Ntau  = 1
         Norb  = int( dims(3) )
         Nbins = int( dims(5) )
      else
         Nunit = int( dims(2) )
         Ntau  = int( dims(3) )
         Norb  = int( dims(4) )
         Nbins = int( dims(6) )
      endif
      deallocate( dims )

      CALL h5gopen_f(file_id, name, grp_id, ierr)
      call read_attribute(grp_id, '.', "dtau", dtau, ierr)
      call read_attribute(grp_id, '.', "Channel", Channel, ierr)
      call h5gclose_f(grp_id, ierr)

      par_dsetname = "lattice"
      write(par_dsetname, '(A, "/lattice")') trim(name)
      CALL h5gopen_f(file_id, par_dsetname, grp_id, ierr)
      call h5ltget_attribute_double_f(grp_id, '.', "a1", a1_p , ierr)
      call h5ltget_attribute_double_f(grp_id, '.', "a2", a2_p , ierr)
      call h5ltget_attribute_double_f(grp_id, '.', "L1", L1_p , ierr)
      call h5ltget_attribute_double_f(grp_id, '.', "L2", L2_p , ierr)
      Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )
      
      attr_name = "N_coord"
      call read_attribute(grp_id, '.', attr_name, Latt_unit%Norb, ierr)
      attr_name = "Norb"
      call read_attribute(grp_id, '.', attr_name, Latt_unit%N_coord, ierr)
      attr_name = "Ndim"
      call read_attribute(grp_id, '.', attr_name, Ndim, ierr)
      allocate(Latt_unit%Orb_pos_p(Latt_unit%Norb, Ndim))
      
      do no = 1, Latt_unit%Norb
         write(attr_name, '("Orbital", I0)') no
         call  h5ltget_attribute_double_f(grp_id, '.', attr_name, Latt_unit%Orb_pos_p(no,:), ierr )
      enddo

      Allocate ( Bins(Nunit,Ntau,Norb,Norb,Nbins), Bins0(Norb,Nbins), sgn(Nbins) )

      dat_ptr = C_LOC(Bins(1,1,1,1,1))
      CALL h5dopen_f(file_id, obs_dsetname, dset_id, ierr)
      CALL h5dget_space_f(dset_id, dataspace, ierr)
      CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, ierr, &
                     mem_space_id = dataspace, file_space_id = dataspace)
      call h5dclose_f(dset_id, ierr)

      dat_ptr = C_LOC(Bins0(1,1))
      CALL h5dopen_f(file_id, bak_dsetname, dset_id, ierr)
      CALL h5dget_space_f(dset_id, dataspace, ierr)
      CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, ierr, &
                     mem_space_id = dataspace, file_space_id = dataspace)
      call h5dclose_f(dset_id, ierr)

      allocate( dims(1) )
      dims(1) = Nbins
      CALL h5dopen_f(file_id, sgn_dsetname, dset_id, ierr)
      CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, sgn, dims, ierr)
      call h5dclose_f(dset_id, ierr)
      deallocate( dims )

      call h5fclose_f(file_id, ierr)

   End Subroutine read_latt_hdf5
#endif

!==============================================================================

   Subroutine jack(func, data, N_skip, N_rebin, best, err)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Performs jackknife error Analysis
!>
!-------------------------------------------------------------------
      Implicit none

      Complex (Kind=Kind(0.d0)), External                :: func
      Complex (Kind=Kind(0.d0)), intent(in)              :: data(:,:)
      Integer, intent(in)                                :: N_skip, N_rebin

      Complex (Kind=Kind(0.d0)), intent(out) :: best, err

      Complex (Kind=Kind(0.d0)), allocatable :: data_r(:,:), j(:), X(:), X2(:)
      Real (Kind=Kind(0.d0)), allocatable :: Rhelp(:)
      real (Kind=Kind(0.d0)) :: XM, XERR
      Integer :: Nobs, Nbins, Nbins_r
      Integer :: N_r, N1, N


      Nobs  = size(data,1)
      Nbins = size(data,2)

      ! Skipping and Rebinning
      Nbins_r = (Nbins-N_skip)/N_rebin
      Allocate( data_r(Nobs,Nbins_r), X(Nobs) )
      N = N_skip
      Do N_r = 1, Nbins_r
         X(:) = cmplx(0.D0,0.d0,kind(0.d0))
         Do N1 = 1, N_rebin
            N = N + 1
            X(:) = X(:) + data(:,N)
         enddo
         data_r(:,N_r) = X(:)/dble(N_rebin)
      enddo

      ! Creating jackknife bins
      Allocate( j(Nbins_r), X2(Nobs) )
      X(:) = cmplx(0.D0,0.d0,kind(0.d0))
      do N = 1, Nbins_r
         X(:) = X(:) + data_r(:,N)
      enddo
      do N = 1, Nbins_r
         X2(:) = (X(:) - data_r(:,N))/dble(Nbins_r-1)
         j(N) = func(X2)
      enddo

      ! Calculate standard deviation of jackknife bins
      Allocate (Rhelp(Nbins_r))
      do N = 1, Nbins_r
         Rhelp(N) = dble(j(N))
      enddo
      call errcalc(Rhelp, xm, xerr)
      best =  cmplx(xm  , 0.d0, kind(0.D0))
      err  =  cmplx(xerr, 0.d0, kind(0.D0))

      do N = 1, Nbins_r
         Rhelp(N) = aimag(j(N))
      enddo
      call errcalc(Rhelp, xm, xerr)
      best =  best + cmplx( 0.d0, xm, kind(0.D0)   )
      err  =  err  + cmplx( 0.d0, xerr, kind(0.D0) )

      err = err * dble(Nbins_r)

      deallocate( data_r, X, j, X2, Rhelp )

   End Subroutine jack

!==============================================================================

   subroutine auto(func, data, N_skip, res)
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Calculates autocorrelation
!-------------------------------------------------------------------
      implicit none

      complex (Kind=Kind(0.d0)), External   :: func
      complex (Kind=Kind(0.d0)), intent(in) :: data(:,:)
      Integer, intent(in)                   :: N_skip

      Real (Kind=Kind(0.d0)), intent(inout) :: res(:)

      Integer                :: N_obs, N_bins, N_auto, ntau, nt
      complex (Kind=Kind(0.d0)), allocatable :: data1(:), data2(:)
      complex (Kind=Kind(0.d0)) :: mean, X1, X2

      N_obs  = size(data,1)
      N_bins = size(data,2) - N_skip
      N_auto = size(res)
      allocate( data1(N_obs), data2(N_obs) )

      do ntau = 1, N_auto
         X1 = 0.0
         X2 = 0.0
         mean = 0.0
         do nt = 1, N_bins - ntau
            data1(:) = data(:,nt+N_skip)
            mean = mean + func(data1)
         enddo
         mean = mean / dble(N_bins - ntau)
         mean = func(mean)

         do nt = 1, N_bins - ntau
            data1(:) = data(:,nt+N_skip)
            data2(:) = data(:,nt+N_skip+ntau)
            X1 = X1 + (func(data1)-mean)*(func(data2)-mean)
            X2 = X2 + (func(data1)-mean)*(func(data1)-mean)
         enddo
      !     X1 = X1 / dble(N_bins - ntau)
      !     X2 = X2 / dble(N_bins - ntau)

         Res(ntau) = dble( X1/X2 )
      enddo


   end subroutine auto

!==============================================================================


   subroutine Cov_tau(name_obs, filename_h5)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Analysis of imaginary time displaced  correlation functions.
!>
!> @param [IN] name_obs Character(len=64)
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!> @param [IN] filename_h5 Character(len=64), optional
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!
!--------------------------------------------------------------------

      Implicit none
      Character (len=*), intent(in) :: name_obs
      Character (len=*), intent(in), optional :: filename_h5

      Character (len=64) :: name_obs2
      Real    (Kind=Kind(0.d0)), allocatable :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer :: Bins_raw(:,:,:,:,:), Bins0_raw(:,:)
      Type (Lattice)   :: Latt
      Type (Unit_cell) :: Latt_unit
      real    (Kind=Kind(0.d0)):: dtau
      Character (len=:), allocatable :: Channel 
      Integer :: i

      i = len(trim(name_obs)) - 4
      name_obs2 = name_obs(:i)

      if( present(filename_h5) ) then
#ifdef HDF5
        call read_latt_hdf5(filename_h5, name_obs, sgn, bins_raw, bins0_raw, Latt, Latt_unit, dtau, Channel)
#endif
      else
        call read_latt(name_obs, sgn, bins_raw, bins0_raw, Latt, Latt_unit, dtau, Channel)
      endif

      Write(6,*)  'Channel is ', Channel

      call ana_tau(name_obs2, sgn, bins_raw, bins0_raw, Latt, Latt_unit, dtau, Channel)
   end subroutine Cov_tau

!==============================================================================

   Subroutine ana_tau(name_obs, sgn, bins_raw, bins0_raw, Latt, Latt_unit, dtau, Channel)
      Implicit none
      Character (len=64), intent(in) :: name_obs
      Real    (Kind=Kind(0.d0)), allocatable, intent(in) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer    , intent(in) :: Bins_raw(:,:,:,:,:)
      Complex (Kind=Kind(0.d0)), pointer    , intent(in) :: Bins0_raw(:,:)
      Type (Lattice)                        , intent(in) :: Latt
      Type (Unit_cell)                      , intent(in) :: Latt_unit
      Real    (Kind=Kind(0.d0))             , intent(in) :: dtau
      Character (len=*)                     , intent(in) :: Channel

      Logical :: PartHole,  L_Back
      Character (len=64) :: File_out, command
      Real    (Kind=Kind(0.d0)), parameter :: Zero=1.D-8
      Integer :: N_skip, N_rebin, N_Cov, N_Back, N_auto
      Integer :: Nbins, LT, Lt_eff,  n_mk
      Integer :: nb, no, no1, no2, n, nt, nt1, ierr, Norb
      Complex (Kind=Kind(0.d0)) :: Z, Zmean, Zerr
      Real    (Kind=Kind(0.d0)), allocatable :: Phase(:)
      Complex (Kind=Kind(0.d0)), allocatable :: V_help_loc(:,:,:,:), Bins_help(:,:,:,:) 
      Real    (Kind=Kind(0.d0)), allocatable :: Xk_p(:,:),  Xk_p1(:) 
      Complex (Kind=Kind(0.d0)), allocatable :: V_help_suscep(:,:,:,:), Weights(:) 
      Complex (Kind=Kind(0.d0)), allocatable :: Background(:,:)
      Complex (Kind=Kind(0.d0)), allocatable :: Xmean(:), Xcov(:,:),  Xmean_st(:),  Xerr_st(:)


      NAMELIST /VAR_errors/ n_skip, N_rebin, N_Cov, N_Back, N_auto

      PartHole = .false.
      if(str_to_upper(Channel) == 'PH' .or. str_to_upper(Channel) ==  'P_PH') PartHole = .true.
      
      N_skip = 1
      N_rebin = 1
      N_Back = 1
      N_auto = 0
      N_Cov  = 0
      OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
      IF (ierr /= 0) THEN
         Write(error_unit,*) 'unable to open <parameters>',ierr
         error stop 1
      END IF
      READ(5,NML=VAR_errors)
      CLOSE(5)

      Nbins = size(bins_raw,5)
      LT    = size(bins_raw,2)

      Write(6, '(A22, I0)') "# of bins: ", Nbins
      nbins = Nbins - n_skip
      Write(6, '(A22, I0)') "Effective # of bins: ", Nbins/N_rebin
      if(Nbins/N_rebin < 2) then
         Write(error_unit,*) "Effective # of bins smaller than 2. Analysis impossible!"
         error stop 1
      endif

      if (PartHole .and. mod(Lt-1,2) == 0 ) then
         Lt_eff = (Lt -1 ) /2   + 1
      elseif (PartHole) then
         Lt_eff = Lt/2
      else
         Lt_eff = Lt
      endif

      ! Allocate  space
      Norb = Latt_unit%Norb
      Allocate (  Phase(Nbins), Xk_p(2,Latt%N), Xk_p1(2),  &
            &     V_help_loc(Lt_eff,Norb,Norb,Nbins), &
            &     background(Norb,Nbins))
      Allocate ( Bins_help(Lt_eff,Norb,Norb,Nbins)  )
      Allocate ( Xmean(Lt_eff), Xcov(Lt_eff,Lt_eff) )

      do n = 1,Latt%N
         xk_p(:,n) = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
      enddo
      do nb = 1, nbins
         Phase (nb) = sgn(nb+n_skip)
      enddo

      ! Do timedisplaced k-resolved ===============================
      do nb = 1, nbins
         do no = 1,Norb
            background(no,nb) = bins0_raw(no,nb+n_skip)*sqrt(dble(Latt%N))
         enddo
      enddo
      !!  Normalization:  \sum_q e^{iqr} ( <O_n(r,t) O_m(0)> - <O_n><O_m> ) =
      !!                  \sum_q e^{iqr}  <O_n(r,t) O_m(0)>  - N \delta_{q,0} <O_n><O_m>
      V_help_loc = cmplx(0.d0,0.d0,kind(0.d0))
      do n = 1,Latt%N
         XK_p1  =  - xk_p(:,n)
         n_mk   =   Inv_K(XK_P1,Latt)
         bins_help  = cmplx(0.d0,0.d0,Kind(0.d0))
         do nb  =  1, nbins 
            do no2 = 1,Latt_unit%Norb
               do no1 = 1,Latt_unit%Norb
                  do nt = 1,Lt_eff
                     if (PartHole) then
                        bins_help(nt,no1,no2,nb) = &
                             &  ( bins_raw(n,nt,no1,no2,nb+n_skip) + conjg(bins_raw(n_mk,Lt-nt+1,no1,no2,nb+n_skip)) ) &
                             &   / cmplx(2.d0,0.d0,Kind(0.d0))
                     else
                        bins_help(nt,no1,no2,nb) =  bins_raw(n,nt,no1,no2,nb+n_skip)
                     endif
                     V_help_loc(nt,no1,no2,nb) =  V_help_loc(nt,no1,no2,nb) +  bins_help(nt,no1,no2,nb)  !  For local 
                  enddo
               enddo
            enddo
         enddo
         
         if (  Xk_p(1,n) >= -zero .and. XK_p(2,n) >= -zero ) then
            if ( sqrt(Xk_p(1,n)**2 + Xk_p(2,n)**2) < 1.D-6 .and. N_Back == 1 ) then
               L_back  = .true. 
               call COV(bins_help, phase, Xcov, Xmean, background, L_back, N_rebin  )
            else
               L_back  = .false.
               call COV(bins_help, phase, Xcov, Xmean, background, L_back, N_rebin )
            endif
            write(File_out,'(A,"_",F4.2,"_",F4.2,"/g_dat")') trim(name_obs), Xk_p(1,n), Xk_p(2,n)
            write(command, '("mkdir -p ",A,"_",F4.2,"_",F4.2)') trim(name_obs), Xk_p(1,n), Xk_p(2,n)
            CALL EXECUTE_COMMAND_LINE(command)
            Open (Unit=10, File=File_out, status="unknown")
            Write(10, '(2(I11), E26.17E3, I11, " ", A)') &
                  & Lt_eff, nbins/N_rebin, real(lt-1,kind(0.d0))*dtau, Latt_unit%Norb, trim(Channel)
            do nt = 1, LT_eff
               Write(10, '(3(E26.17E3))') &
                     & dble(nt-1)*dtau,  dble(Xmean(nt)), sqrt(abs(dble(Xcov(nt,nt))))
            enddo
            If (N_cov == 1) Then ! print covarariance
               Do nt = 1,LT_eff
                  Do nt1 = 1,LT_eff
                     Write(10, '(E25.17E3)') dble(Xcov(nt,nt1))
                  Enddo
               Enddo
            Endif
            close(10)
         endif
      enddo


      ! Do timedisplaced r=0 ===============================
      do nb = 1, nbins
         do no = 1,Norb
            background(no,nb) = bins0_raw(no,nb+n_skip)
         enddo
      enddo
      !!  Normalization:   <O_n(0,t) O_m(0)> - <O_n><O_m> 
      V_help_loc =  V_help_loc/dble(Latt%N)
      if (N_Back == 1) then
         L_back  = .true. 
         call COV(V_help_loc, phase, Xcov, Xmean, background, L_back, N_rebin  )
      else
         L_back  = .false. 
         call COV(V_help_loc, phase, Xcov, Xmean, background, L_back, N_rebin  )
      endif
      write(File_out,'(A,"_R0/g_dat")') trim(name_obs)
      write(command, '("mkdir -p ",A,"_R0")') trim(name_obs)
      CALL EXECUTE_COMMAND_LINE(command)
      Open (Unit=10,File=File_out,status="unknown")
      Write(10, '(2(I11), E26.17E3, I11, " ", A)') &
            & LT_eff, nbins/N_rebin, real(lt-1,kind(0.d0))*dtau, Latt_unit%Norb, trim(Channel)
      do nt = 1, LT_eff
         Write(10, '(3(E26.17E3))') &
               & dble(nt-1)*dtau,  dble(Xmean(nt)), sqrt(abs(dble(Xcov(nt,nt))))
      enddo
      If (N_cov == 1) Then ! Print  covariance
         Do nt = 1,LT_eff
            Do nt1 = 1,LT_eff
               Write(10, '(E25.17E3)') dble(Xcov(nt,nt1))
            Enddo
         Enddo
      Endif
      close(10)
      Deallocate( Xmean, Xcov, V_help_loc) 

      ! Do susceptibilities ===============================
      Allocate(V_help_suscep(1,Norb,Norb,Nbins), Weights(Norb) )
      Allocate(Xmean(1), Xcov(1,1), Xmean_st(Latt%N), Xerr_st(Latt%N) ) 
      Weights=cmplx(1.d0,0.d0,kind(0.d0))
      Z = Latt%N*(Lt_eff -1)
      if (PartHole) Z = Z*cmplx(2.d0,0.d0,Kind(0.d0))
      Z=sqrt(Z)
      do nb = 1, nbins
         do no = 1,Norb
            background(no,nb) = bins0_raw(no,nb+n_skip)*Z
         enddo
      enddo
      !!  Normalization:  \int_0^beta \sum_q e^{iqr} ( <O_n(r,t) O_m(0)> - <O_n><O_m> ) =
      !!                  \int_0^beta \sum_q e^{iqr}  <O_n(r,t) O_m(0)>  - N \beta \delta_{q,0} <O_n><O_m>
      do n = 1,Latt%N
         V_help_suscep =   cmplx(0.d0,0.d0,Kind(0.d0))
         do nb = 1, nbins
            Z = cmplx(0.d0,0.d0,kind(0.d0))
            Do nt = 1,Lt_eff -1
               do no = 1,Latt_unit%Norb
                  do no1 = 1,Latt_unit%Norb
                     Z = Z + cmplx(0.5d0,0.d0,Kind(0.d0)) * ( bins_raw(n,nt,no,no1,nb+n_skip) + bins_raw(n,nt+1,no,no1,nb+n_skip) )
                     V_help_suscep(1,no,no1,nb) =  V_help_suscep(1,no,no1,nb) +  &
                          &   cmplx(0.5d0,0.d0,Kind(0.d0)) * ( bins_raw(n,nt,no,no1,nb+n_skip) + bins_raw(n,nt+1,no,no1,nb+n_skip) ) 
                  enddo
               enddo
            enddo
         enddo
         if (PartHole)  V_help_suscep   =  V_help_suscep * cmplx(2.d0,0.d0,Kind(0.d0))
         if ( sqrt(Xk_p(1,n)**2 + Xk_p(2,n)**2) < 1.D-6 .and. N_Back == 1 ) then
            L_back  = .true. 
            call COV(V_help_suscep, phase, Xcov, Xmean, background, L_back, N_rebin, Weights  )
         else
            L_back  = .false. 
            call COV(V_help_suscep, phase, Xcov, Xmean, background, L_back, N_rebin, Weights  )
         endif
         Xmean_st(n) = Xmean(1)*dtau
         Xerr_st(n)  =  Sqrt(Xcov(1,1))*dtau
      enddo
      write(File_out,'(A,"_tauJK")') trim(name_obs)
      Open (Unit=33,File=File_out ,status="unknown")
      Do n = 1,Latt%N
         Write(33, '(6(E26.17E3))') &
               &   Xk_p(1,n), Xk_p(2,n), dble(XMean_st(n)), dble (Xerr_st(n)), &
               &                        aimag(XMean_st(n)), aimag(Xerr_st(n))
         
      enddo
      Close(33)
      deallocate(Xmean, Xcov,Xmean_st, Xerr_st,V_Help_suscep, Weights) 
      ! EndDo susceptibilities ===============================


      ! Deallocate space ===============================
      Deallocate ( Phase, Xk_p,  Xk_p1, Bins_help)
      Deallocate ( background)

   end Subroutine ana_tau

!==============================================================================

   subroutine Cov_eq(name_obs, filename_h5)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Analysis program for equal time observables.
!>
!
!--------------------------------------------------------------------

      Implicit none
      Character (len=*), intent(in) :: name_obs
      Character (len=*), intent(in), optional :: filename_h5

      Real    (Kind=Kind(0.d0)), allocatable :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer :: Bins_raw(:,:,:,:,:), Bins0_raw(:,:)
      Type (Lattice)   :: Latt
      Type (Unit_cell) :: Latt_unit
      Real (Kind=Kind(0.d0)) :: dtau
      Character (len=:), allocatable      :: Channel
      
      if( present(filename_h5) ) then
#ifdef HDF5
         call read_latt_hdf5(filename_h5, name_obs, sgn, bins_raw, bins0_raw, Latt, Latt_unit, dtau, Channel)
#endif
      else
         call read_latt(name_obs, sgn, bins_raw, bins0_raw, Latt, Latt_unit, dtau, Channel)
      endif
      call ana_eq(name_obs, sgn, bins_raw, bins0_raw, Latt, Latt_unit)

   end subroutine Cov_eq

!==============================================================================

   Subroutine ana_eq(name, sgn, bins_raw, bins0_raw, Latt, Latt_unit)
     Implicit none
     Character (len=*)                     , intent(in) :: name
     Real    (Kind=Kind(0.d0)), allocatable, intent(in) :: sgn(:)
     Complex (Kind=Kind(0.d0)), pointer    , intent(in) :: Bins_raw(:,:,:,:,:)
     Complex (Kind=Kind(0.d0)), pointer    , intent(in) :: Bins0_raw(:,:)
     Type (Lattice)                        , intent(in) :: Latt
     Type (Unit_cell)                      , intent(in) :: Latt_unit
     
     Character (len=64) :: File_out
     Integer :: N_skip, N_rebin, N_Cov, N_Back, N_auto
     Integer :: Nbins
     Integer :: i, n, nb, no, no1, ierr
     Type  (Mat_C), allocatable :: Bins (:,:), Bins_R(:,:)
     Complex (Kind=Kind(0.d0)), allocatable :: Phase(:)
     Complex (Kind=Kind(0.d0)), allocatable :: V_help(:,:), V_help_R(:,:)
     Complex (Kind=Kind(0.d0)), allocatable :: Bins0(:,:)
     Real (Kind=Kind(0.d0)), allocatable :: Xk_p_s(:,:)
     Real (Kind=Kind(0.d0)), allocatable :: AutoCorr(:),En(:)
     Real    (Kind=Kind(0.d0)) :: Xk_p(2), Xr_p(2)
     Complex (Kind=Kind(0.d0)) :: Z, Xmean, Xerr, Xmean_r, Xerr_r
     Real (Kind=Kind(0.d0)) :: Xm,Xe
     procedure (func_c), pointer :: f_ptr => Background_eq
     
      NAMELIST /VAR_errors/ N_skip, N_rebin, N_Cov, N_Back, N_auto

      N_skip = 1
      N_rebin = 1
      N_Back = 1
      N_auto = 0
      N_Cov  = 0
      OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
      IF (ierr /= 0) THEN
         Write(error_unit,*) 'unable to open <parameters>',ierr
         error stop 1
      END IF
      READ(5,NML=VAR_errors)
      CLOSE(5)

      Nbins = size(bins_raw,5)

      Write(6, '(A22, I0)') "# of bins: ", Nbins
      Nbins  = Nbins - n_skip
      Write(6, '(A22, I0)') "Effective # of bins: ", Nbins/N_rebin
      N_auto=min(N_auto,Nbins/3)
      if(Nbins/N_rebin < 2) then
         Write(error_unit,*) "Effective # of bins smaller than 2. Analysis impossible!"
         error stop 1
      endif

      ! Allocate  space
      Allocate( bins(Latt%N,Nbins), bins_r(Latt%N,Nbins), Phase(Nbins) )
      Allocate( V_help(3,Nbins), V_help_R(3,Nbins), Bins0(Nbins, Latt_unit%Norb) )
      Do n = 1,Latt%N
         do nb = 1,nbins
            Call Make_Mat(bins  (n,nb), Latt_unit%Norb)
            Call Make_Mat(bins_r(n,nb), Latt_unit%Norb)
            bins_r(n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
            bins  (n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
         Enddo
      Enddo

      allocate( Xk_p_s(2, Latt%N) )
      do n = 1,Latt%N
         Xk_p_s(:,n) = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
      enddo

      Bins0 = cmplx(0.d0,0.d0,kind(0.d0))
      do nb = 1, nbins + n_skip
         if (nb > n_skip ) then
            Phase(nb-n_skip) = cmplx(sgn(nb),0.d0,kind(0.d0))
            Do no = 1,Latt_unit%Norb
               if (N_Back == 1 ) Bins0(nb-n_skip,no) = Bins0_raw(no,nb)
            enddo
            do n = 1,Latt%N
               do no = 1,Latt_unit%Norb
                  do no1 = 1,Latt_unit%Norb
                     bins(n,nb-n_skip)%el(no,no1) = Bins_raw(n,1,no,no1,nb)
                  enddo
               enddo
!!$               FFA:  Legacy                
!!$               Xk_p(:) = Xk_p_s(:,n)
!!$               if ( sqrt(Xk_p(1)**2 + Xk_p(2)**2) < 1.D-6 ) then
!!$                  do no = 1,Latt_unit%Norb
!!$                     do no1 = 1,Latt_unit%Norb
!!$                        bins(n,nb-n_skip)%el(no,no1)  =  bins(n,nb-n_skip)%el(no,no1) !-  &
!!$                        ! &        cmplx(dble(Latt%N),0.d0,kind(0.d0))*Bins0(nb-n_skip,no)*Bins0(nb-n_skip,no1) &
!!$                        ! &        /Phase(nb-n_skip)
!!$                     enddo
!!$                  enddo
!!$               endif
            enddo
         endif
      enddo
      N_auto=min(N_auto,Nbins/3)


      Call Fourier_K_to_R(bins,bins_r,Latt)
      
!!$#ifdef test
!!$      ! Setup symmetries for square lattice.
!!$      do n = 1,Latt%N
!!$         n1 = n
!!$         Write(6, "(2(E26.17E3))") Xk_p(1,n1), Xk_p(2,n1)
!!$         do m = 1,4
!!$            n1 = Rot90(n1, Xk_p, Latt%N)
!!$            Write(6, "(I11, 2(E26.17E3))") n1, Xk_p(1,n1), Xk_p(2,n1)
!!$         enddo
!!$         Write(6,*)
!!$      enddo
!!$#endif
      write(File_out,'(A,A)') trim(name), "JK"
      Open (Unit=33,File=File_out ,status="unknown")
      write(File_out,'(A,A)') trim(name), "JR"
      Open (Unit=34,File=File_out ,status="unknown")
      Do n = 1,Latt%N
         Xk_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
         Xr_p = dble(Latt%list (n,1))*Latt%a1_p + dble(Latt%list (n,2))*Latt%a2_p
         Write(33, '(2(E26.17E3))')  Xk_p(1), Xk_p(2)
         Write(34, '(2(E26.17E3))')  Xr_p(1), Xr_p(2)
         Do no = 1,Latt_unit%Norb
            do no1 = 1,Latt_unit%Norb
               do nb = 1,Nbins
                  V_help(1,nb) = bins  (n,nb)%el(no,no1)
               enddo
               if ( sqrt(Xk_p(1)**2 + Xk_p(2)**2) < 1.D-6 ) then
                  do nb = 1,Nbins
                     V_help(2,nb) = Bins0(nb,no)*Latt%N
                     V_help(3,nb) = Bins0(nb,no1)
                  enddo
               else
                  do nb = 1,Nbins
                     V_help(2,nb) = 0.0d0
                     V_help(3,nb) = 0.0d0
                  enddo
               endif
               call ERRCALCJ( V_help, Phase,XMean, XERR, N_rebin, f_ptr )
               Write(33, "(2(I11), 4(E26.17E3))") &
                     &  no, no1, dble(XMean), dble(XERR), aimag(XMean), aimag(XERR)
               do nb = 1,Nbins
                  V_help_r(1,nb) = bins_r(n,nb)%el(no,no1)
                  V_help_r(2,nb) = Bins0(nb,no)
                  V_help_r(3,nb) = Bins0(nb,no1)
               enddo
               call ERRCALCJ( V_help_R,Phase, XMean_r, XERR_r, N_rebin, f_ptr )
               Write(34, "(2(I11), 4(E26.17E3))") &
                     &  no, no1, dble(XMean_r), dble(XERR_r), aimag(XMean_r), aimag(XERR_r)
            enddo
         enddo
      enddo

      Close(33)
      Close(34)

      if ( N_auto > 0 ) then
         ALLOCATE(AutoCorr(N_auto))
         ALLOCATE(EN(Nbins))
         Do n = 1,Latt%N
            Xk_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
            if (Xk_p(1) >= -1.d-8 .and. XK_p(2) >= -1.d-8) then
               write(File_out,'(A,"_Auto_Tr_",F4.2,"_",F4.2)') trim(name), Xk_p(1), Xk_p(2)
               OPEN (UNIT=21, FILE=File_out, STATUS='unknown')
               WRITE(21,*)
               do nb = 1,Nbins
                  Z=0
                  do no = 1,Latt_unit%Norb
                     Z = Z+bins  (n,nb)%el(no,no)
                  enddo
                  En(nb)=dble(Z)
               enddo
               Call AUTO_COR(En,AutoCorr)
               do i = 1,N_auto
                  CALL ERRCALCJ(En,XM, XE,i)
                  write(21, "(I11, 2(E26.17E3))") i, AutoCorr(i), Xe
               enddo
               CLOSE(21)
            endif
         enddo
         DEALLOCATE(AutoCorr)
      endif
      
    end Subroutine ana_eq


    subroutine Cov_vec(name_obs, filename_h5)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Analysis program for scalar observables
!
!--------------------------------------------------------------------

      Implicit none
      Character (len=*), intent(in) :: name_obs
      Character (len=*), intent(in), optional :: filename_h5
      
      Real    (Kind=Kind(0.d0)), allocatable :: sgn_raw(:)
      Complex (Kind=Kind(0.d0)), pointer     :: Bins_raw(:,:)
      Character (len=64)                     :: analysis_mode
      
      if( present(filename_h5) ) then
#ifdef HDF5
         call read_vec_hdf5(filename_h5, name_obs, sgn_raw, bins_raw, analysis_mode)
#endif
      else
         call read_vec(name_obs, sgn_raw, bins_raw, analysis_mode)
      endif
      call ana_vec(name_obs, sgn_raw, bins_raw, analysis_mode)
      
    END subroutine Cov_vec
    
!==============================================================================

    subroutine ana_vec(name, sgn_raw, bins_raw, analysis_mode)
      Implicit none
      Character (len=*), intent(in) :: name
      Real    (Kind=Kind(0.d0)), allocatable, intent(inout) :: sgn_raw(:)
      Complex (Kind=Kind(0.d0)), pointer,     intent(inout) :: bins_raw(:,:)
      Character (len=64),                     intent(in)    :: analysis_mode
      
      REAL    (Kind=Kind(0.d0)), DIMENSION(:),   ALLOCATABLE :: EN, sgn
      REAL    (Kind=Kind(0.d0)), DIMENSION(:,:),   ALLOCATABLE :: EN_f_arg
      REAL    (Kind=Kind(0.d0)) :: XM, XERR
      
      Complex (Kind=Kind(0.d0)), Allocatable  :: Bins(:,:)
      REAL    (Kind=Kind(0.d0)), Allocatable  :: AutoCorr(:)
      Integer :: Nobs, Nobs_output, data_range
      Integer :: Nbins, Nbins_eff, I, IOBS, N_Back

      Integer :: N_skip, N_rebin, N_Cov, ierr, N_auto
      Character (len=64) :: File_out
      NAMELIST /VAR_errors/   N_skip, N_rebin, N_Cov, N_Back, N_auto

      !New Stuff for Autocorrelation
      REAL(Kind=Kind(0.d0)), DIMENSION(:)  , ALLOCATABLE :: vec, vec_err
      
!       abstract interface
!          function func (X)
!             real (Kind=Kind(0.d0)) :: func
!             real (Kind=Kind(0.d0)), allocatable, intent (in) :: X(:)
!          end function func
!       end interface
      

      procedure (func_r), pointer :: f_ptr => null ()

      N_skip = 1
      N_rebin = 1
      N_Back = 1
      N_auto = 0
      N_Cov  = 0
      OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
      IF (ierr /= 0) THEN
         Write(error_unit,*) 'unable to open <parameters>',ierr
         error stop 1
      END IF
      READ(5,NML=VAR_errors)
      CLOSE(5)

      Nobs  = size(bins_raw, 1)
      Nbins = size(bins_raw, 2)
      
      if (str_to_upper(analysis_mode) == 'IDENTITY') then
         f_ptr => identity
         Nobs_output = Nobs
         data_range  = 0
      elseif(str_to_upper(analysis_mode) == 'RENYI_ENTROPIE') then
         f_ptr => entanglement
         Nobs_output = Nobs
         data_range  = 0
      elseif(str_to_upper(analysis_mode) == 'MUTUAL_INFORMATION') then
         if (Nobs .ne. 3) then
            Write(error_unit,*) 'Evaluating the mutual information between A and B requires the &
                 &   entanglement entropies of A, B and the union of A and B, i.e. Nobs=4 (3 + 1 for the phase)'
            error stop 1
         endif
         f_ptr => mutinf
         Nobs_output = 1
         data_range  = 2
      else
         Write(error_unit,*) 'Unknown observable function! Continue with identity operation.'
         f_ptr => identity
         Nobs_output = Nobs
         data_range  = 0
      endif
      
      Write(6, '(A22, I0)') "# of bins: ", Nbins
      Nbins_eff  = Nbins - n_skip
      Write(6, '(A22, I0)') "Effective # of bins: ", Nbins_eff/N_rebin
      N_auto=min(N_auto,Nbins_eff/3)
      if(Nbins_eff/N_rebin < 2) then
         Write(error_unit,*) "Effective # of bins smaller than 2. Analysis impossible!"
         error stop 1
      endif
      
      ! Allocate  space
      Allocate ( Bins(Nobs,Nbins_eff), sgn(Nbins_eff) )
      
      do i =1,Nbins_eff
         Bins(:,i) = Bins_raw(:,i+n_skip)
         sgn(i) = sgn_raw(i+n_skip)
      enddo
      
      write(File_out,'(A,A)') trim(name), "J"
      OPEN (UNIT=21, FILE=File_out, STATUS='unknown')
      WRITE(21,*) 'Effective number of bins, and bins: ', Nbins_eff/N_rebin, Nbins
      ALLOCATE (EN(Nbins_eff), EN_f_arg(data_range+1,Nbins_eff), vec(NOBS), vec_err(NOBS))
      DO IOBS = 1,Nobs_output
         EN(:) = Real(Bins(IOBS,:), kind(0.d0)) ! not used any more, too be deleted
         EN_f_arg(:,:) = Real(Bins(IOBS:IOBS+data_range,:), kind(0.d0)) !+data_range
         CALL ERRCALCJ(EN_f_arg,sgn,XM,XERR,N_Rebin,f_ptr)
!          CALL ERRCALCJ(EN,sgn,XM,XERR,N_Rebin)
         vec    (IOBS) = XM
         vec_err(IOBS) = XERR
         WRITE(21,*)
         WRITE(21, "(I11, 2(E26.17E3))") IOBS, XM,  XERR
      ENDDO
      CALL ERRCALCJ(sgn, XM,XERR,N_Rebin)
      WRITE(21,*)
      WRITE(21, "(I11, 2(E26.17E3))") Nobs_output+1, XM,  XERR
      CLOSE(21)
      
      if(N_auto>0) then
         ALLOCATE(AutoCorr(N_auto))
         DO IOBS = 1,NOBS
            write(File_out,'(A,A,I1.1)') trim(name), '_Auto_', iobs
            write(*,*) File_out
            OPEN (UNIT=21, FILE=File_out, STATUS='unknown')
            WRITE(21,*)
            EN(:) = Real(Bins(IOBS,:), kind(0.d0))
            Call AUTO_COR(EN,AutoCorr)
            do i = 1,N_auto
               CALL ERRCALCJ(EN,XM,XERR,i)
               write(21, "(I11, 2(E26.17E3))") i, AutoCorr(i), Xerr
            enddo
            CLOSE(21)
         ENDDO
         DEALLOCATE(AutoCorr)
      endif
      
      DEALLOCATE (EN, EN_f_arg,vec,vec_err,sgn_raw,sgn,Bins_raw,Bins)
      
    END subroutine ana_vec
    
    complex (Kind=Kind(0.d0)) function Background_eq(X)
      
      Implicit None
      complex (Kind=Kind(0.d0)), allocatable, intent (in) :: X(:)
      
      Background_eq = X(1) - X(2)*X(3)

    end function Background_eq
   
    complex (Kind=Kind(0.d0)) function Background_sus(X)
     
      Implicit None
      complex (Kind=Kind(0.d0)), allocatable, intent (in) :: X(:)
      
      integer :: Norb, no, no1
      
      Norb = size(X,1)-1
      
      Background_sus = X(1)
      do no = 1,Norb
         do no1 = 1,Norb
            Background_sus = Background_sus - X(1+no)*X(1+no1)
         enddo
      enddo
      
    end function Background_sus
   
    Real (Kind=Kind(0.d0)) function mutinf(X)
      
      Implicit None
      Real (Kind=Kind(0.d0)), allocatable, intent (in) :: X(:)
      
      mutinf = log(X(3)/(X(1)*X(2)))
      
    end function mutinf
    
    Real (Kind=Kind(0.d0)) function identity(X)
      
      Implicit None
      Real (Kind=Kind(0.d0)), allocatable, intent (in) :: X(:)
      
      identity = X(1)
      
    end function identity
    
    Real (Kind=Kind(0.d0)) function entanglement(X)
      
      Implicit None
      Real (Kind=Kind(0.d0)), allocatable, intent (in) :: X(:)
      
      entanglement = -log(X(1))
      
    end function entanglement
  end module ana_mod
