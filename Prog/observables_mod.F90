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


     Module Observables

#if !defined HDF5 && !defined OBS_LEGACY
#define OBS_LEGACY 1
#endif
       Use runtime_error_mod
       Use Lattices_v3, only: Unit_cell, Lattice
       use iso_fortran_env, only: output_unit, error_unit

       Type :: Obser_Vec
!>  Data structure for
!>  < O_n >  n : =1, size(Obs,1)
          !private
          Integer                     :: N                    ! Number of measurements
          real      (Kind=Kind(0.d0)) :: Ave_Sign             ! Averarge sign
          complex   (Kind=Kind(0.d0)), pointer :: Obs_vec(:)  ! Vector of observables
          Character (len=64) :: File_Vec                      ! Name of file in which the bins will be written out
          Character (len=64) :: analysis_mode                 ! How to analyze the observable
          Character (len=64), allocatable :: description(:)   ! Optional short description
       contains
          !procedure :: make        => Obser_vec_make
          procedure :: init        => Obser_vec_init
          procedure :: print_bin   => print_bin_vec
          procedure :: measure     => Obser_vec_measure
       end type Obser_Vec


       Type :: Obser_Latt
!>  Data structure for
!>  < O^{dagger}(i,tau)_n O(j,0)_m>  - < O_n> <O_m>
!>  where it is assumed that translation symmetry as specified by the lattice Latt is present.
!>  Obs_Latt(i-j,tau,n,m) = < O^{dagger}(i,tau)_n O(j,0)_m>
!>  Obs_Latt0(n) = < O_n>
!>  For equal   time correlation functions, tau runs from 1,1
!>  For unequal time correlation functions, tau runs from 1,Ltrot+1
          Integer            :: N                                    ! Number of measurements
          Real      (Kind=Kind(0.d0)) :: Ave_Sign                    ! Averarge sign
          complex   (Kind=Kind(0.d0)), pointer :: Obs_Latt (:,:,:,:) ! i-j, tau, norb, norb
          complex   (Kind=Kind(0.d0)), pointer :: Obs_Latt0(:)       ! norb
          Character (len=64) :: File_Latt                            ! Name of file in which the bins will be written out
          Type (Lattice),       pointer :: Latt                      ! Pointer to Bravais lattice
          Type (Unit_cell),     pointer :: Latt_unit                 ! Pointer to unit cell
          Real      (Kind=Kind(0.d0))   :: dtau                      ! Imaginary time step
          Character (len=2)  :: Channel    ! Type of observable. Possible values:
                                           ! - T0: zero temperature
                                           ! - P:  finite temperature particle
                                           ! - PH: finite temperature particle-hole
                                           ! - PP: finite temperature particle-particle
       contains
          !procedure :: make        => Obser_latt_make
          procedure :: init        => Obser_latt_init
          procedure :: print_bin   => print_bin_latt
       end type Obser_Latt

       Contains

         Subroutine Obser_Latt_make(Obs, Nt, Filename, Latt, Latt_unit, Channel, dtau)
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
!> @param [IN] Channel, Character(len=2)
!> \verbatim
!>  MaxEnt channel. Only relevant for time displaced observables.
!> \endverbatim
!> @param [IN] dtau, Real(Kind=Kind(0.d0))
!> \verbatim
!>  Imaginary time step. Only relevant for time displaced observables.
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           type(Obser_Latt), Intent(INOUT)      :: Obs
           Integer,           Intent(IN)         :: Nt
           Character(len=64), Intent(IN)         :: Filename
           Type(Lattice),     Intent(IN), target :: Latt
           Type(Unit_cell),   Intent(IN), target :: Latt_unit
           Character(len=2),  Intent(IN)         :: Channel
           Real(Kind=Kind(0.d0)),  Intent(IN)    :: dtau
           Allocate (Obs%Obs_Latt(Latt%N, Nt, Latt_unit%Norb, Latt_unit%Norb))
           Allocate (Obs%Obs_Latt0(Latt_unit%Norb))
           Obs%File_Latt = Filename
           Obs%Latt => Latt
           Obs%Latt_unit => Latt_unit
           Obs%Channel = Channel
           Obs%dtau = dtau
         end subroutine Obser_Latt_make
!--------------------------------------------------------------------

         Subroutine Obser_Latt_Init(Obs)
           Implicit none
           class(Obser_Latt), intent(INOUT) :: Obs
           Obs%Obs_Latt  = cmplx(0.d0,0.d0,kind(0.d0))
           Obs%Obs_Latt0 = cmplx(0.d0,0.d0,kind(0.d0))
           Obs%N         = 0
           Obs%Ave_Sign  = 0.d0
         end subroutine Obser_Latt_Init

!--------------------------------------------------------------------

         Subroutine Obser_Vec_make(Obs, N, Filename, analysis_mode, description)
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
           Implicit none
           type(Obser_vec), intent(INOUT) :: Obs
           Integer, Intent(IN)             :: N
           Character (len=64), Intent(IN)  :: Filename
           Character (len=64), Intent(IN), optional :: analysis_mode
           Character (len=64), Intent(IN), optional :: description(:)
           
           Allocate (Obs%Obs_vec(N))
           Obs%File_Vec = Filename
           if(present(analysis_mode)) then
             Obs%analysis_mode = analysis_mode
           else
             Obs%analysis_mode = 'identity'
           endif
           if(present(description)) then
             allocate(Obs%description(size(description, 1)))
             Obs%description = description
           endif
         end subroutine Obser_Vec_make
!--------------------------------------------------------------------

         Subroutine Obser_Vec_Init(Obs)
           Implicit none
           class(Obser_vec), intent(INOUT) :: Obs
           Obs%Obs_vec = cmplx(0.d0,0.d0,kind(0.d0))
           Obs%N       = 0
           Obs%Ave_Sign= 0.d0
         end subroutine Obser_Vec_Init

!--------------------------------------------------------------------

         Subroutine Obser_vec_measure(obs, value, Phase)
           Implicit none

           class (Obser_vec),        Intent(Inout) :: Obs
           complex(Kind=Kind(0.d0)), Intent(In)    :: value(:)  ! Vector of observables
           complex(Kind=Kind(0.d0)), Intent(IN), optional    :: Phase
            !Local
           Complex (Kind=Kind(0.d0)) :: ZP, ZS

           obs%N = obs%N + 1

           if ( present(Phase) ) then
              ZP = PHASE/Real(Phase, kind(0.D0))
              ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

              obs%Ave_sign  = obs%Ave_sign + real(ZS, kind(0.D0))
              obs%obs_vec   = obs%obs_vec + value *ZS*ZP
           else
              obs%Ave_sign  = obs%Ave_sign + 1.d0
              obs%obs_vec   = obs%obs_vec  + value
           endif

         end Subroutine  Obser_vec_measure

!--------------------------------------------------------------------

         Subroutine  Print_bin_Latt(Obs, Group_Comm)
           Use Lattices_v3
#if defined HDF5
           Use hdf5
           Use alf_hdf5
#endif
#if defined MPI
           Use mpi
#endif
           Implicit none

           class(Obser_Latt),        Intent(Inout)   :: Obs
           Integer,                  Intent(In)      :: Group_Comm

           ! Local
           Integer :: Ns, Nt, no, no1, I, Ntau
           Complex (Kind=Kind(0.d0)), pointer     :: Tmp(:,:,:,:)
           Real    (Kind=Kind(0.d0))              :: x_p(2)
           Complex (Kind=Kind(0.d0))              :: Sign_bin
           Character (len=64) :: File_pr,  File_suff, File_aux, tmp_str
           logical            :: File_exists
#ifdef HDF5
           Character (len=7), parameter  :: File_h5 = "data.h5"
           Character (len=64)            :: filename, groupname, obs_dsetname, bak_dsetname, sgn_dsetname
           INTEGER(HID_T)                :: file_id, group_id
           logical                       :: link_exists
           INTEGER                       :: hdferr
           INTEGER(HSIZE_T), allocatable :: dims(:)
           TYPE(C_PTR)                   :: dat_ptr
           real(Kind=Kind(0.d0)), target :: sgn
#endif
#ifdef MPI
           Complex (Kind=Kind(0.D0)), allocatable :: Tmp1(:)
           Complex (Kind=Kind(0.d0)) :: Z
           Real    (Kind=Kind(0.d0)) :: X
           Integer         :: Ierr, Isize, Irank
           INTEGER         :: irank_g, isize_g, igroup

           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
           call MPI_Comm_rank(Group_Comm, irank_g, ierr)
           call MPI_Comm_size(Group_Comm, isize_g, ierr)
           igroup           = irank/isize_g
#endif

           Ns    = Size(Obs%Obs_Latt,1)
           Ntau  = Size(Obs%Obs_Latt,2)
           if ( .not. (Obs%Latt%N == Ns ) ) then
              Write(error_unit,*) 'Error in Print_bin_Latt'
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           If (Ntau == 1) then
              File_suff = "_eq"
           else
              File_suff = "_tau"
           endif
           write(File_pr, '(A,A)') trim(Obs%File_Latt), Trim(File_suff)
#if defined HDF5
           groupname = File_pr
           filename = File_h5
#endif
           Allocate (Tmp(Ns, Ntau, Obs%Latt_unit%Norb, Obs%Latt_unit%Norb))
           Obs%Obs_Latt  = Obs%Obs_Latt /dble(Obs%N   )
           Obs%Obs_Latt0 = Obs%Obs_Latt0/dble(Obs%N*Ns*Ntau)
           Obs%Ave_sign  = Obs%Ave_Sign /dble(Obs%N   )

#if defined(MPI)
           I = Obs%Latt%N * Ntau * Obs%Latt_unit%Norb * Obs%Latt_unit%Norb
           Tmp = cmplx(0.d0, 0.d0, kind(0.D0))
           CALL MPI_REDUCE(Obs%Obs_Latt,Tmp,I,MPI_COMPLEX16,MPI_SUM, 0,Group_Comm,IERR)
           Obs%Obs_Latt = Tmp/DBLE(ISIZE_g)

           I = 1
           X = 0.d0
           CALL MPI_REDUCE(Obs%Ave_sign,X,I,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
           Obs%Ave_sign = X/DBLE(ISIZE_g)

           I = Obs%Latt_unit%Norb
           Allocate(Tmp1(Obs%Latt_unit%Norb))
           Tmp1 = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Obs%Obs_Latt0,Tmp1,I,MPI_COMPLEX16,MPI_SUM, 0,Group_Comm,IERR)
           Obs%Obs_Latt0 = Tmp1/DBLE(ISIZE_g)
           Deallocate(Tmp1)

           If (Irank_g == 0 ) then
#endif
#if defined TEMPERING
              write(File_pr ,'(A,I0,A,A,A)') "Temp_",igroup,"/",trim(Obs%File_Latt),trim(File_suff )
#if defined HDF5
              write(filename ,'(A,I0,A,A)') "Temp_",igroup,"/",trim(File_h5)
#endif
#endif

              do nt = 1, Ntau
                 do no = 1, Obs%Latt_unit%Norb
                    do no1 = 1, Obs%Latt_unit%Norb
                       Call Fourier_R_to_K(Obs%Obs_Latt(:,nt,no,no1), Tmp(:,nt,no,no1), Obs%Latt)
                    enddo
                 enddo
              enddo

#if defined OBS_LEGACY
              write(File_aux, '(A,A)') trim(File_pr), "_info"
              inquire(file=File_aux, exist=File_exists)
              if (.not.File_exists) then
                 11 format(A20, ': ', A)
                 12 format(A20, ': ', I10)
                 13 format(A20, ': ', *(E26.17E3))
                 14 format(A20, ': ')
                 15 format((E26.17E3))
                 open(10, file=File_aux, status='new')
                 write(tmp_str, '(A, A)') trim(Obs%File_Latt), trim(File_suff)
                 write(10, 11) 'Observable', trim(tmp_str)
                 write(10, 11) 'Channel', trim(Obs%Channel)
                 write(10, 12) 'Ntau', Ntau
                 write(10, 13) 'dtau', Obs%dtau
                 write(10, '(A)') '       ====== Bravais Lattice ======'
                 write(10, 12) 'Unit cells', Obs%Latt%N
                 write(10, 13) 'L1', Obs%Latt%L1_p
                 write(10, 13) 'L2', Obs%Latt%L2_p
                 write(10, 13) 'a1', Obs%Latt%a1_p
                 write(10, 13) 'a2', Obs%Latt%a2_p
                 write(10, '(A)') '       ========= Unit cell ========='
                 write(10, 12) 'Coordination number', Obs%Latt_unit%N_coord
                 write(10, 12) 'Number of orbitals', Obs%Latt_unit%Norb
                 write(10, 12) 'Ndim', size(Obs%Latt_unit%Orb_pos_p, 2)
                 do no = 1, Obs%Latt_unit%Norb
                    write(tmp_str, '("Orbital ",I0)') no
                    write(10, 14, advance='no') trim(tmp_str)
                    do i = 1, size(Obs%Latt_unit%Orb_pos_p, 2)
                       write(10, 15, advance='no') Obs%Latt_unit%Orb_pos_p(no,i)
                    enddo
                    write(10, *)
                 enddo
                 close(10)
              endif

              Open(Unit=10, File=File_pr, status="unknown",  position="append")
              If ( Ntau == 1 ) then
                 Write(10, '(E25.17E3, 2(I11))') Obs%Ave_sign, Obs%Latt_unit%Norb, Obs%Latt%N
              else
                 Write(10, '(E25.17E3, 3(I11), E26.17E3)') Obs%Ave_sign, Obs%Latt_unit%Norb, Obs%Latt%N, Ntau, Obs%dtau
              endif
              Do no = 1, Obs%Latt_unit%Norb
                 Write(10, '("(", E25.17E3, ",", E25.17E3, ")")')  Obs%Obs_Latt0(no)
              enddo
              do I = 1, Obs%Latt%N
                 x_p = dble(Obs%Latt%listk(i,1))*Obs%Latt%b1_p + dble(Obs%Latt%listk(i,2))*Obs%Latt%b2_p
                 Write(10, '(E25.17E3, 1x, E25.17E3)') X_p(1), X_p(2)
                 Do nt = 1, Ntau
                    do no = 1, Obs%Latt_unit%Norb
                       do no1 = 1, Obs%Latt_unit%Norb
                          Write(10, '("(", E25.17E3, ",", E25.17E3, ")")') tmp(I,nt,no,no1)
                       enddo
                    enddo
                 enddo
              enddo
              close(10)
#endif

#if defined HDF5
              write(obs_dsetname,'(A,A,A)') trim(groupname), "/obser"
              write(bak_dsetname,'(A,A,A)') trim(groupname), "/back"
              write(sgn_dsetname,'(A,A,A)') trim(groupname), "/sign"

              CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

              !Writes the lattice to HDF5 file if it doesn't already exist
              call write_latt(file_id, Obs%Latt, Obs%Latt_Unit)

              CALL h5lexists_f(file_id, groupname, link_exists, hdferr)
              if ( .not. link_exists ) then
                !Create Group for observable
                CALL h5gcreate_f (file_id, groupname, group_id, hdferr)
                call write_attribute(group_id, '.', "dtau", Obs%dtau, hdferr)
                call write_attribute(group_id, '.', "Channel", Obs%Channel, hdferr)
                call write_latt(group_id, Obs%Latt, Obs%Latt_Unit)
                CALL h5gclose_f (group_id, hdferr)

                !Create Dataset for data
                allocate( dims(6) )
                dims = [2, Obs%Latt%N, Ntau, Obs%Latt_unit%Norb, Obs%Latt_unit%Norb, 0]
                CALL init_dset(file_id, obs_dsetname, dims, .true.)
                deallocate( dims )

                !Create Dataset for background
                allocate( dims(3) )
                dims = [2, Obs%Latt_unit%Norb, 0]
                CALL init_dset(file_id, bak_dsetname, dims, .true.)
                deallocate( dims )

                !Create Dataset for sign
                allocate( dims(1) )
                dims = [0]
                CALL init_dset(file_id, sgn_dsetname, dims, .False.)
                deallocate( dims )
              endif

              !Write data
              dat_ptr = C_LOC(tmp(1,1,1,1))
              CALL append_dat(file_id, obs_dsetname, dat_ptr)

              !Write background
              dat_ptr = C_LOC(Obs%Obs_Latt0(1))
              CALL append_dat(file_id, bak_dsetname, dat_ptr)

              !Write sign
              sgn = Obs%Ave_sign
              dat_ptr = C_LOC(sgn)
              CALL append_dat(file_id, sgn_dsetname, dat_ptr)

              CALL h5fclose_f(file_id, hdferr)
#endif
#if defined MPI
           Endif
#endif

           deallocate (Tmp)

         End Subroutine Print_bin_Latt

!--------------------------------------------------------------------

         Subroutine  Print_bin_Vec(Obs,Group_Comm)
#if defined MPI
           Use mpi
#endif
#if defined HDF5
           Use hdf5
           Use alf_hdf5
#endif
           Implicit none

           class(Obser_vec), intent(Inout) :: Obs
           Integer, INTENT(IN)  :: Group_Comm

           ! Local
           Integer :: I
           Character (len=64) :: File_pr, File_suff, File_aux
           logical            :: File_exists

#if defined HDF5
           Character (len=7), parameter  :: File_h5 = "data.h5"
           Character (len=64)            :: filename, groupname, obs_dsetname, sgn_dsetname
           INTEGER(HID_T)                :: file_id, group_id
           logical                       :: link_exists
           INTEGER                       :: hdferr
           INTEGER(HSIZE_T), allocatable :: dims(:)
           TYPE(C_PTR)                   :: dat_ptr
           real(Kind=Kind(0.d0)), target :: sgn
#endif
#if defined MPI
           Integer        :: Ierr, Isize, Irank, No
           INTEGER        :: irank_g, isize_g, igroup
           Complex  (Kind=Kind(0.d0)), allocatable :: Tmp(:)
           Real     (Kind=Kind(0.d0)) :: X


           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
           call MPI_Comm_rank(Group_Comm, irank_g, ierr)
           call MPI_Comm_size(Group_Comm, isize_g, ierr)
           igroup           = irank/isize_g
#endif
           Obs%Obs_vec  = Obs%Obs_vec /dble(Obs%N)
           Obs%Ave_sign = Obs%Ave_sign/dble(Obs%N)
           write(File_pr, '(A,A)') trim(Obs%File_Vec), "_scal"
#if defined HDF5
           groupname = File_pr
           filename = File_h5
#endif

#if defined MPI
           No = size(Obs%Obs_vec, 1)
           Allocate (Tmp(No) )
           Tmp = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Obs%Obs_vec,Tmp,No,MPI_COMPLEX16,MPI_SUM, 0,Group_Comm,IERR)
           Obs%Obs_vec = Tmp/DBLE(ISIZE_g)
           deallocate (Tmp )

           I = 1
           X = 0.d0
           CALL MPI_REDUCE(Obs%Ave_sign,X,I,MPI_REAL8,MPI_SUM, 0,Group_comm,IERR)
           Obs%Ave_sign = X/DBLE(ISIZE_g)

           if (Irank_g == 0 ) then
#endif

#if defined TEMPERING
              write(File_pr,'(A,I0,A,A,A)') "Temp_",igroup,"/",trim(Obs%File_Vec), "_scal"
#if defined HDF5
              write(filename ,'(A,I0,A,A)') "Temp_",igroup,"/",trim(File_h5)
#endif
#endif

#if defined OBS_LEGACY
              write(File_aux, '(A,A)') trim(File_pr), "_info"
              inquire(file=File_aux, exist=File_exists)
              if (.not.File_exists) then
                 open(10, file=File_aux, status='new')
                 write(10, '(A)') '====== Analysis Mode ======'
                 write(10, '(A)') trim(Obs%analysis_mode)
                 if(allocated(Obs%description)) then
                   write(10, '(A)') '====== Description ======'
                   do i=1, size(Obs%description, 1)
                     write(10, '(A)') trim(Obs%description(i))
                   enddo
                 endif
                 close(10)
              endif
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              !WRITE(10,*) size(Obs%Obs_vec,1)+1, (Obs%Obs_vec(I), I=1,size(Obs%Obs_vec,1)), Obs%Ave_sign
              write(10, '(I10)', advance='no') size(Obs%Obs_vec,1)+1
              do I=1,size(Obs%Obs_vec,1)
                 write(10, '(" (",E25.17E3,",",E25.17E3,")")', advance='no') Obs%Obs_vec(I)
              enddo
              write(10, '(E26.17E3)') Obs%Ave_sign
              close(10)
#endif

#if defined HDF5
              write(obs_dsetname,'(A,A,A)') trim(groupname), "/obser"
              write(sgn_dsetname,'(A,A,A)') trim(groupname), "/sign"

              CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

              !Check if observable already exists in hdf5 file
              CALL h5lexists_f(file_id, groupname, link_exists, hdferr)

              if ( .not. link_exists ) then
                !Create Group for observable and write auxiliary info
                CALL h5gcreate_f (file_id, groupname, group_id, hdferr)
                call write_attribute(group_id, '.', "analysis_mode", Obs%analysis_mode, hdferr)
                if(allocated(Obs%description)) then
                  call write_comment(group_id, '.', "description", Obs%description, hdferr)
                endif
                CALL h5gclose_f (group_id, hdferr)

                !Create Dataset for data
                allocate( dims(3) )
                dims = [2, size(Obs%Obs_vec,1), 0]
                CALL init_dset(file_id, obs_dsetname, dims, .true.)
                deallocate( dims )

                !Create Dataset for sign
                allocate( dims(1) )
                dims = [0]
                CALL init_dset(file_id, sgn_dsetname, dims, .false.)
                deallocate( dims )
              endif

              !Write data
              dat_ptr = C_LOC(Obs%Obs_vec(1))
              CALL append_dat(file_id, obs_dsetname, dat_ptr)

              !Write sign
              sgn = Obs%Ave_sign
              dat_ptr = C_LOC(sgn)
              CALL append_dat(file_id, sgn_dsetname, dat_ptr)

              CALL h5fclose_f(file_id, hdferr)
#endif
#if defined MPI
           endif
#endif

         End Subroutine Print_bin_Vec

       end Module Observables
