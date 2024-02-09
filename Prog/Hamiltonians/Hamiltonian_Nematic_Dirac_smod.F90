!  Copyright (C) 2021-2022 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version

!--------------------------------------------------------------------
!> @author
!> J. Schwab
!>
!> @brief
!> This module defines the  Hamiltonian and observables of several models
!> that correcond to Dirac systems undergoing a nematic quantum phase transition.
!>
!--------------------------------------------------------------------
   submodule (Hamiltonian_main) ham_Nematic_Dirac_mod
     
#ifdef MPI
      use mpi
#endif
      Use Operator_mod
      Use Observables
      Use Lattices_v3
      Use Random_Wrap

      Implicit none

      type, extends(ham_base) :: ham_Nematic_Dirac
      contains
         procedure, nopass :: Ham_Set
         procedure, nopass :: S0
         procedure, nopass :: Delta_S0_global
         procedure, nopass :: Alloc_obs
         procedure, nopass :: Obser
         procedure, nopass :: ObserT
         procedure, nopass :: Global_move
         procedure, nopass :: Hamiltonian_set_nsigma
         procedure, nopass :: Pr_obs
         procedure, nopass :: Init_obs
#ifdef HDF5
         procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Nematic_Dirac


      type :: Obser_hist
!>  Defines histogram observable
         Integer                     :: N        ! Number of measurements
         real      (Kind=Kind(0.d0)) :: Ave_Sign ! Averarge sign
         Character (len=64)          :: name     ! Name of file in which the bins will be written out

         !real      (Kind=Kind(0.d0)) :: val ! temporary storage for value

         integer                     :: N_classes        ! Number of classes, in which the observable gets counted
         real      (Kind=Kind(0.d0)) :: upper, lower     ! Range, in which the observable should lie.
         real      (Kind=Kind(0.d0)) :: d_class          ! width of one class

         integer, allocatable        :: counts(:)        !
         integer                     :: N_below, N_above ! Counts occurences outside of range

      contains
         procedure :: make        => Obser_hist_make
         procedure :: init        => Obser_hist_init
         procedure :: print_bin   => print_bin_hist
         procedure :: measure     => Obser_hist_measure
      end type Obser_hist
      
      integer, parameter :: dp=kind(0.d0)  ! double precision

      !#PARAMETERS START# VAR_Nematic_Dirac
      !Integer :: N_SUN = 2       ! SU(N) symmetry
      real(dp) :: dtau = 0.1d0    ! Imaginary time step size
      Character (len=64) :: Global_type = '' ! Type of global update. Possible values: 'Wolff', 'Geo', 'switch', 'flip'
      Integer  :: L1 = 4          ! Size of lattice in a1 direction
      Integer  :: L2 = 4          ! Size of lattice in a2 direction
      Integer  :: Model_vers = 1  ! Version of model. 1: C_2v model, 2: C_4v model
      real(dp) :: ham_t = 1.d0    ! Hopping amplitude of fermions
      real(dp) :: beta = 10.d0    ! Reciprocal temperature
      real(dp) :: Ham_h = 3.d0    ! Ising transverse field
      real(dp) :: Ham_J = 1.d0    ! Ferromagnetic Ising interaction
      real(dp) :: Ham_xi = 1.d0   ! Coupling strength Ising spins <-> fermions
      real(dp) :: Ham_xi2 = 0.d0  ! Static fermion hopping "distortion"
      real(dp) :: Ham_chem = 0.d0 ! Chemical potential
      real(dp) :: Global_J = 1.d0 ! J for proposing global updates
      real(dp) :: Global_h = 3.d0 ! h for proposing global updates
      real(dp) :: Phi_1 = 0.d0    ! Twisted boundary in a1 direction
      real(dp) :: Phi_2 = 0.d0    ! Twisted boundary in a2 direction
      Character (len=64) :: init_type = 'random' ! How to initialize Ising field. Possile values: 'random', 'up', 'down', 'updown'
      !#PARAMETERS END#

#ifdef MPI
      Integer        :: Isize, Irank, irank_g, isize_g, igroup
      Integer        :: STATUS(MPI_STATUS_SIZE)
#endif
      
     Integer                  :: N_coord, Norb
     Type (Lattice),   target :: Latt
     Type (Unit_cell), target :: Latt_unit
     Integer, allocatable     :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell


     real (Kind=Kind(0.d0)) :: Model_sign

     !>    Private variables for observing Z_x_ising
     Real (Kind=Kind(0.d0)) :: eq_x_ising, neq_x_ising
     Integer                         :: nBlub, nBlub2

     !>    Storage for the Ising action
     Real (Kind=Kind(0.d0)) :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
     Integer, allocatable   :: Ising_nnlist(:,:)

     !>    Variables for the Wolff cluster update
     Real (Kind=Kind(0.d0)) :: Wolff_addProb_space, Wolff_addProb_tau
     Integer                :: N_ising
     !>    Variables for the Geometric cluster update
     Real (Kind=Kind(0.d0)) :: Geo_AddProb_space, Geo_AddProb_tau
     Integer :: R_init(2)

     !>    Experimenting
     Integer :: n_global = 0

     !> container for histogram observables
     Type (Obser_hist), dimension(:), allocatable :: Obs_hist

     !>    For measuring histograms once per sweep
     logical :: measure_hist_on_next_ntau1 = .true.

   contains

   module subroutine Ham_Alloc_Nematic_Dirac()
     allocate(ham_Nematic_Dirac::ham)
   end subroutine Ham_Alloc_Nematic_Dirac

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Nematic_Dirac_read_write_parameters.F90"


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Prints out the bins.  Modified from default to include histograms
!-------------------------------------------------------------------
      Subroutine  Pr_obs(LTAU)

         Implicit none

         Integer,  Intent(In) ::  Ltau

         !Local
         Integer :: I


         if ( allocated(Obs_scal) ) then
           Do I = 1,Size(Obs_scal,1)
              Call Print_bin_Vec(Obs_scal(I), Group_Comm)
           enddo
         endif
         if ( allocated(Obs_hist) ) then
            Do I = 1,Size(Obs_hist,1)
               Call  Obs_hist(I)%print_bin(Group_Comm)
            enddo
         endif
         if ( allocated(Obs_eq) ) then
           Do I = 1,Size(Obs_eq,1)
              Call Print_bin_Latt(Obs_eq(I), Group_Comm)
           enddo
         endif
         if ( allocated(Obs_tau) ) then
           Do I = 1,Size(Obs_tau,1)
              Call Print_bin_Latt(Obs_tau(I), Group_Comm)
           enddo
         endif

      end Subroutine Pr_obs


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Initializes observables to zero before each bins.
!> Modified from default to include histograms
!-------------------------------------------------------------------
      Subroutine  Init_obs(Ltau)

         Implicit none
         Integer, Intent(In) :: Ltau

         ! Local
         Integer :: I

         if ( allocated(Obs_scal) ) then
           Do I = 1,Size(Obs_scal,1)
              Call Obser_vec_Init(Obs_scal(I))
           Enddo
         endif

         if ( allocated(Obs_hist) ) then
            Do I = 1,Size(Obs_hist,1)
               Call Obs_hist(I)%init()
            Enddo
         endif

         if ( allocated(Obs_eq) ) then
           Do I = 1,Size(Obs_eq,1)
              Call Obser_Latt_Init(Obs_eq(I))
           Enddo
         endif

         if ( allocated(Obs_tau) ) then
           Do I = 1,Size(Obs_tau,1)
              Call Obser_Latt_Init(Obs_tau(I))
           Enddo
         Endif

      end Subroutine Init_obs

!--------------------------------------------------------------------

          Subroutine Obser_hist_make(obs, N_classes, upper, lower ,name)
            Implicit none
            class (Obser_hist)         , intent(INOUT) :: obs
            Integer                    , Intent(IN)    :: N_classes
            real      (Kind=Kind(0.d0)), Intent(IN)    :: upper, lower
            Character (len=64)         , Intent(IN)    :: name

            obs%N_classes = N_classes
            obs%upper = upper
            obs%lower = lower

            obs%d_class = (upper - lower) / real(N_classes, Kind=Kind(0.d0))

            Allocate ( obs%counts(N_classes) )
            Obs%name = name
          end subroutine Obser_hist_make
!--------------------------------------------------------------------

          Subroutine Obser_hist_Init(obs)
            Implicit none
            class (Obser_hist), intent(INOUT) :: obs
            obs%N       = 0
            obs%Ave_Sign= 0.d0

            obs%counts(:) = 0
            obs%N_below = 0
            obs%N_above = 0
          end subroutine Obser_hist_Init

!--------------------------------------------------------------------

          Subroutine  Obser_hist_measure(obs, value, sign_in)
            Implicit none

            class (Obser_hist),       Intent(Inout)   :: Obs
            Real (Kind=Kind(0.d0)),   Intent(In)      :: value, sign_in

            Integer :: n_x

            obs%N = obs%N + 1
            obs%Ave_sign  =  obs%Ave_sign + sign_in

            n_x = ceiling( (value - obs%lower) / obs%d_class )
            if ( n_x < 1 ) then
              obs%N_below = obs%N_below + 1
            elseif ( n_x > obs%N_classes ) then
              obs%N_above = obs%N_above + 1
            else
              obs%counts(n_x) = obs%counts(n_x) + 1
            endif

          end Subroutine  Obser_hist_measure

!--------------------------------------------------------------------

          Subroutine Print_bin_hist(obs,Group_Comm)
#if defined MPI
            Use mpi
#endif
#if defined HDF5
            Use hdf5
            Use alf_hdf5
#endif
            Implicit none

            class (Obser_hist), Intent(Inout) :: obs
            Integer           , Intent(In)    :: Group_Comm

            ! Local
            Integer :: I
            Character (len=64)           :: File_pr
            real(Kind=Kind(0.d0)), pointer :: counts_out(:)
            real(Kind=Kind(0.d0)), target :: N_above_out, N_below_out

#if defined HDF5
            Character (len=7), parameter  :: File_h5 = "data.h5"
            Character (len=64)            :: filename, groupname, obs_dsetname, sgn_dsetname, above_dsetname, below_dsetname
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
            real     (Kind=Kind(0.d0)), allocatable :: Tmp(:)
            Real     (Kind=Kind(0.d0)) :: X

            CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
            CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
            call MPI_Comm_rank(Group_Comm, irank_g, ierr)
            call MPI_Comm_size(Group_Comm, isize_g, ierr)
            igroup           = irank/isize_g
#endif
            allocate( counts_out(obs%N_classes) )
            counts_out  = dble(Obs%counts ) / dble(Obs%N)
            N_above_out = dble(Obs%N_above) / dble(Obs%N)
            N_below_out = dble(Obs%N_below) / dble(Obs%N)
            Obs%Ave_sign = Obs%Ave_sign/dble(Obs%N)
            write(File_pr, '(A,A)') trim(Obs%name), "_hist"
#if defined HDF5
            groupname = File_pr
            filename = File_h5
#endif

#if defined MPI
            No = size(counts_out, 1)
            Allocate (Tmp(No) )
            Tmp = cmplx(0.d0,0.d0,kind(0.d0))
            CALL MPI_REDUCE(counts_out,Tmp,No,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
            counts_out = Tmp/DBLE(ISIZE_g)
            deallocate (Tmp )

            I = 1
            X = 0.d0
            CALL MPI_REDUCE(Obs%Ave_sign,X,I,MPI_REAL8,MPI_SUM, 0,Group_comm,IERR)
            Obs%Ave_sign = X/DBLE(ISIZE_g)

            X = 0.d0
            CALL MPI_REDUCE(N_above_out,X,I,MPI_REAL8,MPI_SUM, 0,Group_comm,IERR)
            N_above_out = X/DBLE(ISIZE_g)

            X = 0.d0
            CALL MPI_REDUCE(N_below_out,X,I,MPI_REAL8,MPI_SUM, 0,Group_comm,IERR)
            N_below_out = X/DBLE(ISIZE_g)

            if (Irank_g == 0 ) then
#endif

#if defined TEMPERING
              write(File_pr,'(A,I0,A,A,A)') "Temp_",igroup,"/",trim(Obs%name), "_hist"
#if defined HDF5
              write(filename ,'(A,I0,A,A)') "Temp_",igroup,"/",trim(File_h5)
#endif
#endif

#if defined HDF5
              write(obs_dsetname,'(A,A,A)') trim(groupname), "/obser"
              write(sgn_dsetname,'(A,A,A)') trim(groupname), "/sign"
              write(above_dsetname,'(A,A,A)') trim(groupname), "/above"
              write(below_dsetname,'(A,A,A)') trim(groupname), "/below"

              CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

              !Check if observable already exists in hdf5 file
              CALL h5lexists_f(file_id, groupname, link_exists, hdferr)

              if ( .not. link_exists ) then
                !Create Group for observable
                CALL h5gcreate_f (file_id, groupname, group_id, hdferr)
                call write_attribute(group_id, '.', "N_classes", obs%N_classes, hdferr)
                call write_attribute(group_id, '.', "upper"    , obs%upper    , hdferr)
                call write_attribute(group_id, '.', "lower"    , obs%lower    , hdferr)
                CALL h5gclose_f  (group_id, hdferr)

                !Create Dataset for data
                allocate( dims(2) )
                dims = (/size(counts_out,1), 0/)
                CALL init_dset(file_id, obs_dsetname, dims, .false.)
                deallocate( dims )

                !Create Dataset for sign
                allocate( dims(1) )
                dims = (/0/)
                CALL init_dset(file_id, sgn_dsetname, dims, .false.)
                deallocate( dims )

                !Create Dataset
                allocate( dims(1) )
                dims = (/0/)
                CALL init_dset(file_id, above_dsetname, dims, .false.)
                deallocate( dims )

                !Create Dataset
                allocate( dims(1) )
                dims = (/0/)
                CALL init_dset(file_id, below_dsetname, dims, .false.)
                deallocate( dims )
              endif

              !Write data
              dat_ptr = C_LOC(counts_out(1))
              CALL append_dat(file_id, obs_dsetname, dat_ptr)

              !Write sign
              sgn = Obs%Ave_sign
              dat_ptr = C_LOC(sgn)
              CALL append_dat(file_id, sgn_dsetname, dat_ptr)

              !Write
              dat_ptr = C_LOC(N_above_out)
              CALL append_dat(file_id, above_dsetname, dat_ptr)

              !Write
              dat_ptr = C_LOC(N_below_out)
              CALL append_dat(file_id, below_dsetname, dat_ptr)

              CALL h5fclose_f(file_id, hdferr)
#else
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              WRITE(10,*) obs%N_classes, obs%upper, obs%lower, N_above_out, N_below_out, &
                          & (counts_out(I), I=1,size(counts_out,1)), Obs%Ave_sign
              close(10)
#endif
#if defined MPI
            endif
#endif
            deallocate( counts_out )


          end Subroutine Print_bin_hist

          Subroutine Obser_Latt_make_norb(Obs, Nt, Norb, Filename, Latt, Channel, dtau)
             Implicit none
             class(Obser_Latt), Intent(INOUT)      :: Obs
             Integer,           Intent(IN)         :: Nt
             Integer,           Intent(IN)         :: Norb
             Character(len=64), Intent(IN)         :: Filename
             Type(Lattice),     Intent(IN), target :: Latt
             Character(len=2),  Intent(IN)         :: Channel
             Real(Kind=Kind(0.d0)),  Intent(IN)    :: dtau
            
             allocate(Obs%Latt_unit)
             Obs%Latt_unit%N_coord = 0
             Obs%Latt_unit%Norb = Norb
             allocate(Obs%Latt_unit%Orb_pos_p(Norb, 2))
             Obs%Latt_unit%Orb_pos_p = 0.d0
             
             Allocate (Obs%Obs_Latt(Latt%N, Nt, Obs%Latt_unit%Norb, Obs%Latt_unit%Norb))
             Allocate (Obs%Obs_Latt0(Obs%Latt_unit%Norb))
             Obs%File_Latt = Filename
             Obs%Latt => Latt
             Obs%Channel = Channel
             Obs%dtau = dtau
          end subroutine Obser_Latt_make_norb

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set()
         Implicit none
         Character (len=64) :: file_info
         integer            :: unit_info, ierr
         
#ifdef MPI
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
#endif

           ! From dynamically generated file "Hamiltonian_Nematic_Dirac_read_write_parameters.F90"
           call read_parameters()

#ifdef MPI
           If (Irank_g == 0 ) then
#endif
#if defined(TEMPERING)
           write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
           file_info = "info"
#endif
           OPEN(newunit=unit_info, file=file_info, status="unknown", position="append")
           Write(unit_info,*) '====================================='
           Write(unit_info,*) 'Model is            : ', 'Nematic_Dirac', Model_vers
           Write(unit_info,*) 'Global_type         : ', Global_type
           Write(unit_info,*) 'L1                  : ', L1
           Write(unit_info,*) 'L2                  : ', L2
           Write(unit_info,*) 'N_SUN               : ', N_SUN
           Write(unit_info,*) 'ham_t               : ', ham_t
           Write(unit_info,*) 'dtau                : ', dtau
           Write(unit_info,*) 'beta                : ', beta
           Write(unit_info,*) 'Ham_h               : ', Ham_h
           Write(unit_info,*) 'Ham_J               : ', Ham_J
           Write(unit_info,*) 'Ham_xi              : ', Ham_xi
           Write(unit_info,*) 'Ham_xi2             : ', Ham_xi2
           Write(unit_info,*) 'Ham_chem            : ', Ham_chem
           Write(unit_info,*) 'Global_J            : ', Global_J
           Write(unit_info,*) 'Global_h            : ', Global_h
           Write(unit_info,*) 'Phi_1               : ', Phi_1
           Write(unit_info,*) 'Phi_2               : ', Phi_2
           Write(unit_info,*) 'init_type           : ', init_type
           close(unit_info)
#ifdef MPI
          endif
#endif

          Call Ham_latt()

          N_FL = 1

          Call Ham_hop()
          Ltrot = nint(beta/dtau)
          Projector = .false.
          Thtrot = 0
          Symm = .false.

          Call Setup_Ising_action()
          call Ham_V()
      End Subroutine Ham_Set

!=============================================================================
        Subroutine Ham_Latt()
          Implicit none
          !Set the lattice

          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no

          select case( Model_vers )
          case( 0, 1, 3 )
            Norb = 2
            N_coord   = 4
            a1_p(1) =  1.D0/sqrt(2.D0)  ; a1_p(2) =  1.D0/sqrt(2.D0)
            a2_p(1) =  1.D0/sqrt(2.D0)  ; a2_p(2) = -1.D0/sqrt(2.D0)
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
            If ( L1 == 1 .or. L2 == 1 ) then
              Write(6,*) ' One dimensional systems not implemented '
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            endif
          case( 2 )
            Norb = 2
            N_coord   = 4
            a1_p(1) = 1.D0  ; a1_p(2) = 0.D0
            a2_p(1) = 0.D0  ; a2_p(2) = 1.D0
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
            If ( L1 == 1 .or. L2 == 1 ) then
              Write(6,*) ' One dimensional systems not implemented '
              CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            endif
          case default
            Write(6,*) "Lattice not yet implemented!"
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end select

          ! This is for the orbital structure.
          Ndim = Latt%N*Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc
             Enddo
          Enddo

        end Subroutine Ham_Latt
        
        
        pure function twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, latt, L1, L2) result(f)
          Implicit none

          Integer      , intent(in) :: i, del_1, del_2, L1, L2
          real(dp)     , intent(in) :: phi_1, phi_2
          type(Lattice), intent(in) :: Latt
          complex(dp) :: f
          
          real(Kind=Kind(0.d0)), parameter :: pi = dacos(-1.d0)
          
          f = cmplx(1.d0,0.d0, kind(0.D0))
          if((latt%list(i,1)+del_1) > L1) then
            f = f * exp(cmplx(0.d0, 2.d0*pi*phi_1, kind(0.d0)))
          endif
          if((latt%list(i,1)+del_1) < 0) then
            f = f * exp(cmplx(0.d0, -2.d0*pi*phi_1, kind(0.d0)))
          endif
          if((latt%list(i,2)+del_2) > L2) then
            f = f * exp(cmplx(0.d0, 2.d0*pi*phi_2, kind(0.d0)))
          endif
          if((latt%list(i,2)+del_2) < 0) then
            f = f * exp(cmplx(0.d0, -2.d0*pi*phi_2, kind(0.d0)))
          endif
        end function twisted_boundary_f

!===================================================================================
        Subroutine Ham_hop()
          Implicit none

          Integer :: I, I1, J1, n, Ncheck, nc, nc1, del_1, del_2
          complex (Kind=Kind(0.d0)) :: t_temp

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
            Do nc = 1,Ncheck
              Call Op_make(Op_T(nc,n),Ndim)
              select case( Model_vers )
              case( 0 )
                DO I = 1, Latt%N
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      J1 = invlist(I,2)
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case (2)
                      J1 = invlist(Latt%nnlist(I, 0,1),2)
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case (3)
                      J1 = invlist(Latt%nnlist(I,-1,1),2)
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case (4)
                      J1 = invlist(Latt%nnlist(I,-1,0),2)
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '
                    end select
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              case( 1, 3 )
                DO I = 1, Latt%N
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      del_1 = 0
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case (2)
                      del_1 = 0
                      del_2 = 1
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case (3)
                      del_1 = -1
                      del_2 = 1
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case (4)
                      del_1 = -1
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '
                    end select
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    J1 = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              case( 2 )
                DO I = 1, Latt%N
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      del_1 = 1
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (2)
                      del_1 = -1
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (3)
                      del_1 = 0
                      del_2 = 1
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case (4)
                      del_1 = 0
                      del_2 = -1
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '
                    end select
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    J1 = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              case default
                Write(6,*) ' This hopping is not yet implemented '
                CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
              end select

                Do I = 1,Ndim
                   Op_T(nc,n)%P(i) = i
                   Op_T(nc,n)%O(i,i) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                Enddo
                Op_T(nc,n)%g = -Dtau
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n))
             enddo
          enddo
        end Subroutine Ham_hop

!===================================================================================
        Subroutine Ham_V()

          Implicit none

          Integer :: nf, I, nc1, del_1, del_2
          complex (Kind=Kind(0.d0)) :: t_temp

          select case( Model_vers )
          case( 0 )
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I,nf),5)
                Op_V(I,nf)%P(1) = Invlist(I,1)
                Do nc1 = 1,N_coord
                  select case (nc1)
                  case (1)
                    Op_V(I,nf)%P(nc1+1) = invlist(I,2)
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (2)
                    Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I, 0,1),2)
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case (3)
                    Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I,-1,1),2)
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (4)
                    Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I,-1,0),2)
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case default
                    Write(6,*) ' Error in  Ham_V '
                  end select
                  Op_V(I,nf)%O(1   ,nc1+1) = t_temp
                  Op_V(I,nf)%O(nc1+1,1   ) = conjg(t_temp)
                Enddo
                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
              Enddo
            Enddo
          case( 1 )
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I,nf),5)
                Op_V(I,nf)%P(1) = Invlist(I,1)
                Do nc1 = 1,N_coord
                  select case (nc1)
                  case (1)
                    del_1 = 0
                    del_2 = 0
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (2)
                    del_1 = 0
                    del_2 = 1
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case (3)
                    del_1 = -1
                    del_2 = 1
                    t_temp =  cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (4)
                    del_1 = -1
                    del_2 = 0
                    t_temp = -cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case default
                    Write(6,*) ' Error in  Ham_V '
                  end select
                  t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                       & latt, L1, L2)
                  Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                  Op_V(I,nf)%O(1   ,nc1+1) = t_temp
                  Op_V(I,nf)%O(nc1+1,1   ) = conjg(t_temp)
                Enddo
                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
              Enddo
            Enddo
          case( 3 )
            Allocate(Op_V(2*Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I         ,nf),3)
                call Op_make(Op_V(I+Latt%N,nf),3)
                Op_V(I       ,nf)%P(1) = Invlist(I,1)
                Op_V(I+Latt%N,nf)%P(1) = Invlist(I,1)
                Do nc1 = 1,4
                  select case (nc1)
                  case (1)
                    del_1 = 0
                    del_2 = 0
                    Op_V(I,nf)%P(2) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I,nf)%O(1,2) = t_temp
                    Op_V(I,nf)%O(2,1) = conjg(t_temp)
                  case (2)
                    del_1 = 0
                    del_2 = 1
                    Op_V(I+Latt%N,nf)%P(2) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I+Latt%N,nf)%O(1,2) = t_temp
                    Op_V(I+Latt%N,nf)%O(2,1) = conjg(t_temp)
                  case (3)
                    del_1 = -1
                    del_2 = 1
                    Op_V(I,nf)%P(3) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp =  cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I,nf)%O(1,3) = t_temp
                    Op_V(I,nf)%O(3,1) = conjg(t_temp)
                  case (4)
                    del_1 = -1
                    del_2 = 0
                    Op_V(I+Latt%N,nf)%P(3) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp = -cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I+Latt%N,nf)%O(1,3) = t_temp
                    Op_V(I+Latt%N,nf)%O(3,1) = conjg(t_temp)
                  case default
                    Write(6,*) ' Error in  Ham_V '
                  end select
                Enddo
                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
                Op_V(I+Latt%N,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I+Latt%N,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I+Latt%N,nf)%type   = 1
                Call Op_set( Op_V(I+Latt%N,nf) )
              Enddo
            Enddo
          case( 2 )
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I,nf),2)
                Op_V(I,nf)%P(1) = invlist(I,1)
                Op_V(I,nf)%P(2) = invlist(I,2)

                Op_V(I,nf)%O(1,2) = cmplx(0.d0, 1.d0, kind(0.D0))
                Op_V(I,nf)%O(2,1) = cmplx(0.d0,-1.d0, kind(0.D0))

                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
              Enddo
            Enddo
          case default
            Write(6,*) ' This interaction is not yet implemented '
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end select
        end Subroutine Ham_V

!===================================================================================
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Single spin flip S0 ratio
!> @details
!> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
!> a spin flip of Operator n on time slice nt
!> @details
!--------------------------------------------------------------------
        Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)
          Implicit none
          !> Operator index
          Integer, Intent(IN) :: n
          !> Time slice
          Integer, Intent(IN) :: nt
          !> New local field on time slice nt and operator index n
          complex (Kind=Kind(0.d0)), Intent(In) :: Hs_new

          Integer :: nt1,I
          S0 = 1.d0
          If ( Op_V(n,1)%type == 1 ) then
             do i = 1,4
                S0 = S0*DW_Ising_space(nsigma%i(n,nt)*nsigma%i(Ising_nnlist(n,i),nt))
             enddo
             nt1 = nt +1
             if (nt1 > Ltrot) nt1 = 1
             S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
             nt1 = nt - 1
             if (nt1 < 1  ) nt1 = Ltrot
             S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
             If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
          endif

        end function S0

!===================================================================================
        Subroutine  Hamiltonian_set_nsigma(Initial_field)

          ! The user can set the initial configuration

          Implicit none
          complex (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(INOUT) :: Initial_field

          Integer :: I, nt, ierr

          !write(*,*) init_type
          Allocate( Initial_field(Size(OP_V,1), Ltrot) )

          select case (init_type)
          case ('random')
            deallocate( Initial_field )
          case ('up')
            Do nt = 1,Ltrot
              Do I = 1,Size(OP_V,1)
                Initial_field(I,nt)  = 1
              enddo
            enddo
          case ('down')
            Do nt = 1,Ltrot
              Do I = 1,Size(OP_V,1)
                Initial_field(I,nt)  = -1
              enddo
            enddo
          case ('updown')
            Do nt = 1,Ltrot
              Do I = 1,Size(OP_V,1)
                Initial_field(I,nt)  = 1
                If ( I > Latt%N ) Initial_field(I,nt)  = -1
              enddo
            enddo
#ifdef MPI
          case ('part_rand_up')
            If ( mod(Irank_g,2) == 0 ) then
              deallocate( Initial_field )
            else
              Initial_field(:,:)  = 1.d0
            endif
#endif
          case default
            Write(6,*) ' Error in  Hamiltonian_set_random_nsigma'
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end select

        end Subroutine Hamiltonian_set_nsigma

!===================================================================================
        Subroutine Setup_Ising_action()

          ! This subroutine sets up lists and arrays so as to enable an
          ! an efficient calculation of  S0(n,nt)

          Implicit none
          Integer :: I

          select case( Model_vers )
          case( 0,1 )
            allocate(Ising_nnlist(Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
            enddo
          case( 3 )
            allocate(Ising_nnlist(2*Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
               Ising_nnlist(I+Latt%N,1) = Latt%nnlist(I, 1, 0) + Latt%N
               Ising_nnlist(I+Latt%N,2) = Latt%nnlist(I, 0, 1) + Latt%N
               Ising_nnlist(I+Latt%N,3) = Latt%nnlist(I,-1, 0) + Latt%N
               Ising_nnlist(I+Latt%N,4) = Latt%nnlist(I, 0,-1) + Latt%N
            enddo
          case( 2 )
            allocate(Ising_nnlist(Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
            enddo
          case default
            Write(6,*) ' Error in Setup_Ising_action '
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end select

          !print*, Global_J, Global_h

          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J)
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J)

          Wolff_addProb_space = 1 - exp(-2.d0*Dtau*Global_J)
          Wolff_addProb_tau   = 1 - tanh(Dtau*Global_h)

          Geo_addProb_space = 1 - exp(-4.d0*Dtau*Ham_J)
          Geo_addProb_tau   = 1 - tanh(Dtau*Ham_h)**2
        End Subroutine Setup_Ising_action
!===================================================================================
        Subroutine  Alloc_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N,Nt, No
          Character (len=64) ::  name
          Character (len=2)  ::  Channel

          nBlub2 = 100

          ! Scalar observables
          If ( Model_vers == 3 ) then
            Allocate ( Obs_scal(11) )
          else
            Allocate ( Obs_scal(9) )
          endif
          Do I = 1,Size(Obs_scal,1)
            select case (I)
            case (1)
              N = 1;   name ="Part"
            case (2)
              N = 1;   name ="ising_z"
            case (3)
              N = 3;   name ="Kin_Pot_E"
            case (4)
              N = 1;   name ="ising_x"
              eq_x_ising  = tanh(Dtau*Ham_h)
              neq_x_ising = 1/tanh(Dtau*Ham_h)
            case (5)
              N = 3;   name ="m"
            case (6)
              N = 2;   name ="chi2"

            case (7)
              N = 1;   name ="ising_z_alt"
            case (8)
              N = 1;   name ="ising_x_alt"
            case (9)
              N = 3;   name ="m_alt"

            case (10)
              N = 1; name ="d_Z"
            case (11)
              N = 1; name ="d_X"
            case default
              Write(6,*) ' Error in Alloc_obs '
            end select
            Call Obser_vec_make(Obs_scal(I), N,name)
          enddo

          ! Histogram observables
          If ( Model_vers == 3 ) then
            Allocate ( Obs_hist(5) )
          else
            Allocate ( Obs_hist(3) )
          endif
          Do I = 1,Size(Obs_hist,1)
            select case (I)
            case (1)
              N = 400; name = "X"
            case (2)
              N = 400; name = "m"
            case (3)
              N = 400; name = "B"
            case (4)
              N = 400; name = "d_Z"
            case (5)
              N = 400; name = "d_X"
            ! case (5)
            !   N = 1; name ="m12"
            ! case (6)
            !   N = 1; name ="X12"
            case default
              Write(6,*) ' Error in Alloc_obs '
            end select
            Call Obs_hist(I)%make(N, 1.5d0, -0.5d0, name)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(5) )
          Do I = 1,Size(Obs_eq,1)
            select case (I)
            case (1)
              if ( Model_vers == 3 ) then
                No = 2
              else
                No = 1
              endif
              name ="IsingX"
              Nt = 1
            case (2)
              if ( Model_vers == 3 ) then
                No = 2
              else
                No = 1
              endif
              name ="IsingZ"
              Nt = 1
            case (3)
              No = Norb; name ="Green"
              Nt = 1
            case (4)
              No = 1;     name ="IsingZT"
              Nt = Ltrot+1
            case (5)
              No = 1; name ="IsingXT"
              Nt = Ltrot+1
            case default
              Write(6,*) ' Error in Alloc_obs '
            end select
            !Nt = 1
            Channel = '--'
            Call Obser_Latt_make_norb(Obs_eq(I), Nt, No, name, Latt, Channel, dtau)
          enddo

          If (Ltau == 1) then
             ! Equal time correlators
             Allocate ( Obs_tau(1) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   No = Norb;  name ="Green"
                case (2)
                   No = Norb;  name ="GreenUp"
                case (3)
                   No = Norb;  name ="GreenDow"
                case (4)
                   No = Norb;  name ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1-2*Thtrot
                Channel = '--'
                Call Obser_Latt_make_norb(Obs_tau(I), Nt, No, name, Latt, Channel, dtau)
             enddo
          endif

        end Subroutine Alloc_obs

!========================================================================
        ! Functions for Global moves.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)

          !>  The input is the field nsigma declared in this module. This routine generates a
          !>  global update with  and returns the propability
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
          !>
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
      !     Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
          type(fields), intent(in) :: nsigma_old
          Integer :: size_cluster
      
          Integer :: i, nt
      
          n_global = n_global + 1
      
          select case( Global_type )
          case ( 'Wolff' )
            call Wolff_cluster_start(size_cluster, T0_Proposal_ratio, nsigma_old)
          case ( 'Geo' )
            call Geo_cluster_start(size_cluster)
          case ( 'switch' )
            size_clust = latt%N*2
            T0_Proposal_ratio = 1
            Do I = 1, latt%N
              Do nt = 1,Ltrot
                nsigma%f(I,nt) = nsigma_old%f(I+latt%N,nt)
                nsigma%f(I+latt%N,nt) = nsigma_old%f(I,nt)
              enddo
            enddo
            if ( mod(n_global,2) == 0 ) nsigma%f(:,:) = -nsigma_old%f(:,:)
          case ( 'flip' )
            T0_Proposal_ratio = 1
            nsigma%f(:,:) = -nsigma_old%f(:,:)
          case default
            write(*,*) "Error in Global_move: Unknown Global_type:", Global_type
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end select
      
          write(*,*) Delta_S0_global(nsigma_old) * T0_Proposal_ratio
      
          !Write(6,*) "cluster finished. Size,Size/N_Spins:", size_cluster, dble(size_cluster)/dble(N_ising*Ltrot)
          size_clust = dble(size_cluster)/dble(N_ising*Ltrot)
          !T0_Proposal_ratio = 1
        End Subroutine Global_move
      
        Subroutine Wolff_cluster_start(size_cluster, T0_Proposal_ratio, nsigma_old)
          implicit none
          Integer, intent(out) :: size_cluster
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio
      !     Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
          type(fields), intent(in) :: nsigma_old
      
          Integer :: n0, nt0, sigma_clust
          Integer, dimension(:), allocatable :: list_R, list_t
          !Real (Kind=Kind(0.d0)) :: Delta_S0
      
          Allocate( list_R(N_ising*Ltrot), list_t(N_ising*Ltrot) )
      
          size_cluster = 0
          n0  = ceiling( dble(N_ising) * RANF_WRAP() )
          nt0 = ceiling( dble(Ltrot)   * RANF_WRAP() )
          sigma_clust = nsigma%i(n0,nt0)
      
          size_cluster = size_cluster + 1
          list_R(size_cluster) = n0
          list_t(size_cluster) = nt0
          call Wolff_cluster_add(n0, nt0, sigma_clust, size_cluster, list_R, list_t)
      
          T0_Proposal_ratio = Wolff_T0_Proposal_ratio(size_cluster, sigma_clust, nsigma_old, list_R, list_t)
          !Delta_S0 = Delta_S0_global(nsigma_old)
          !Write(6,*) T0_Proposal_ratio*Delta_S0
      
        End Subroutine Wolff_cluster_start
      
        Function Wolff_T0_Proposal_ratio(size_cluster, sigma_clust, nsigma_old, list_R, list_t)
          implicit none
          Real (Kind=Kind(0.d0)) :: Wolff_T0_Proposal_ratio
          Integer, intent(in) :: size_cluster, sigma_clust
      !     Integer, dimension(:,:), allocatable, intent(in) :: nsigma_old
          type(fields), intent(in) :: nsigma_old
          Integer, dimension(:), allocatable, intent(in) :: list_R, list_t
          Integer :: x, n, nt, i, n1, nt1, n_bond_space, n_bond_tau
          n_bond_space = 0
          n_bond_tau = 0
      
          do x = 1, size_cluster
            n  = list_R(x)
            nt = list_t(x)
      
            do i = 1,4
              n1 = Ising_nnlist(n,i)
              If ( nsigma%i(n1,nt) == nsigma_old%i(n1,nt) ) then
                If ( nsigma%i(n1,nt) == sigma_clust ) then
                  n_bond_space = n_bond_space+1
                else
                  n_bond_space = n_bond_space-1
                endif
              Endif
            enddo
      
            nt1 = nt +1
            if (nt1 > Ltrot) nt1 = 1
            If ( nsigma%i(n,nt1) == nsigma_old%i(n,nt1) ) then
              if ( nsigma%i(n,nt1) == sigma_clust ) then
                n_bond_tau = n_bond_tau+1
              else
                n_bond_tau = n_bond_tau-1
              endif
            Endif
      
            nt1 = nt - 1
            if (nt1 < 1  ) nt1 = Ltrot
            If ( nsigma%i(n,nt1) == nsigma_old%i(n,nt1) ) then
              if ( nsigma%i(n,nt1) == sigma_clust ) then
                n_bond_tau = n_bond_tau+1
              else
                n_bond_tau = n_bond_tau-1
              endif
            Endif
          enddo
      
          Wolff_T0_Proposal_ratio = (1-Wolff_addProb_space)**(-n_bond_space) * (1-Wolff_addProb_tau)**(-n_bond_tau)
      
        End Function Wolff_T0_Proposal_ratio
      
        Recursive Subroutine Wolff_cluster_add(n, nt, sigma_clust, size_cluster, list_R, list_t)
          Implicit none
          Integer, intent(in) :: n, nt, sigma_clust
          Integer :: i, n1, nt1
          Integer, intent(inout) :: size_cluster
          Integer, dimension(:), allocatable, intent(inout) :: list_R, list_t
      
          nsigma%f(n,nt) = -nsigma%f(n,nt)
      
          do i = 1,4
            n1 = Ising_nnlist(n,i)
            If ( ( nsigma%i(n1,nt) == sigma_clust ) .and. ( Wolff_addProb_space > RANF_WRAP() ) ) then
              size_cluster = size_cluster + 1
              list_R(size_cluster) = n1
              list_t(size_cluster) = nt
              call Wolff_cluster_add(n1, nt, sigma_clust, size_cluster, list_R, list_t)
            Endif
          enddo
      
          nt1 = nt +1
          if (nt1 > Ltrot) nt1 = 1
          If ( ( nsigma%i(n,nt1) == sigma_clust ) .and. ( Wolff_addProb_tau > RANF_WRAP() ) ) then
            size_cluster = size_cluster + 1
            list_R(size_cluster) = n
            list_t(size_cluster) = nt1
            call Wolff_cluster_add(n, nt1, sigma_clust, size_cluster, list_R, list_t)
          Endif
      
          nt1 = nt - 1
          if (nt1 < 1  ) nt1 = Ltrot
          If ( ( nsigma%i(n,nt1) == sigma_clust ) .and. ( Wolff_addProb_tau > RANF_WRAP() ) ) then
            size_cluster = size_cluster + 1
            list_R(size_cluster) = n
            list_t(size_cluster) = nt1
            call Wolff_cluster_add(n, nt1, sigma_clust, size_cluster, list_R, list_t)
          Endif
        End Subroutine Wolff_cluster_add
      
      
        Subroutine Geo_cluster_start(size_cluster)
          Implicit none
          Integer, intent(out) :: size_cluster
          Integer :: i1, i1t, i2, i2t, sigma_i1, sigma_i2
          Integer :: j1, j1t, j2, j2t, i
          Logical, allocatable, dimension(:,:) :: Geo_cluster
      
          size_cluster = 0
      
          ! Start by randomly selecting two space-time-points (i1, i1t) and (i2, i2t)
          ! The center of these two points becomes the symmetry center of the move
          i1  = ceiling( dble(latt%N) * RANF_WRAP() )
          i1t = ceiling( dble(Ltrot)  * RANF_WRAP() )
          i2  = ceiling( dble(latt%N) * RANF_WRAP() )
          i2t = ceiling( dble(Ltrot)  * RANF_WRAP() )
      
          sigma_i1 = nsigma%i(i1, i1t)
          sigma_i2 = nsigma%i(i2, i2t)
          If(sigma_i1 == sigma_i2) return
      
          Allocate (Geo_cluster(latt%N,Ltrot))
          !Allocate (Geo_cluster_test(latt%N,Ltrot))
          Geo_cluster(:,:) = .false.
          !Geo_cluster_test(:,:) = .false.
      
          !Defining symmetry
          R_init(1) = latt%List(i1,1) + latt%List(i2,1)
          R_init(2) = latt%List(i1,2) + latt%List(i2,2)
          if ( R_init(1) < -(L1-1)/2  ) R_init(1) = R_init(1) + L1
          if ( R_init(1) >   L1/2     ) R_init(1) = R_init(1) - L1
          if ( R_init(2) < -(L2-1)/2  ) R_init(2) = R_init(2) + L2
          if ( R_init(2) >   L2/2     ) R_init(2) = R_init(2) - L2
          Write(6,*) "Symmetry-defining vector", R_init(1), R_init(2)
      
          nsigma%f(i1, i1t) = sigma_i2
          nsigma%f(i2, i2t) = sigma_i1
          Geo_cluster(i1, i1t) = .true.
          Geo_cluster(i2, i2t) = .true.
          size_cluster = 2
      
          do i = 1,4
            j1 = Ising_nnlist(i1,i)
            j2 = Find_geo_partner_space(j1)
            call Geo_cluster_tryadd(sigma_i1, j1, i1t, j2, i2t, .false., size_cluster, Geo_cluster)
          enddo
      
          j1t = i1t + 1
          if (j1t > Ltrot) j1t = 1
          j2t = i2t - 1
          if (j2t < 1) j2t = Ltrot
          call Geo_cluster_tryadd(sigma_i1, i1, j1t, i2, j2t, .true., size_cluster, Geo_cluster)
      
          j1t = i1t - 1
          if (j1t < 1) j1t = Ltrot
          j2t = i2t + 1
          if (j2t > Ltrot) j2t =1
          call Geo_cluster_tryadd(sigma_i1, i1, j1t, i2, j2t, .true., size_cluster, Geo_cluster)
      
          deallocate(Geo_cluster)
        End Subroutine Geo_cluster_start
      
        Integer Function Find_geo_partner_space(j1)
          Implicit none
          Integer, intent(in)  :: j1
          Integer :: R_j2(2)
          Integer :: test(2)
      
          R_j2(1) = R_init(1) - latt%List(j1,1)
          R_j2(2) = R_init(2) - latt%List(j1,2)
          if ( R_j2(1) < -(L1-1)/2  ) R_j2(1) = R_j2(1) + L1
          if ( R_j2(1) >   L1/2     ) R_j2(1) = R_j2(1) - L1
          if ( R_j2(2) < -(L2-1)/2  ) R_j2(2) = R_j2(2) + L2
          if ( R_j2(2) >   L2/2     ) R_j2(2) = R_j2(2) - L2
          Find_geo_partner_space = latt%invlist(R_j2(1),R_j2(2))
      
          test(1) = R_init(1) - (latt%list(j1,1) + latt%list(Find_geo_partner_space,1))
          If ( mod(test(1),L1) .ne. 0 ) then
            Write(6,*) "ERROR: this is not zero:", test(1)
          endif
          test(2) = R_init(2) - (latt%list(j1,2) + latt%list(Find_geo_partner_space,2))
          If (  mod(test(2),L2) .ne. 0 ) then
            Write(6,*) "ERROR: this is not zero:", test(2)
          endif
        End Function Find_geo_partner_space
      
      
        Recursive Subroutine Geo_cluster_tryadd(sigma_i1, j1, j1t, j2, j2t, in_tau, size_cluster, Geo_cluster)
          Implicit none
          Integer, intent(in) :: sigma_i1, j1, j1t, j2, j2t
          Logical, intent(in) :: in_tau
          Integer, intent(inout) :: size_cluster
          Logical, allocatable, dimension(:,:), intent(inout) :: Geo_cluster
          Integer :: sigma_j1, sigma_j2
          Integer :: i, k1, k1t, k2, k2t
      
          If( Geo_cluster(j1, j1t) ) return
      
          sigma_j1 = nsigma%i(j1, j1t)
          sigma_j2 = nsigma%i(j2, j2t)
          If(sigma_j1 == sigma_j2) return
      
          If( in_tau ) then
            If ( RANF_WRAP() > Geo_addProb_tau   ) return
          else
            If ( RANF_WRAP() > Geo_addProb_space ) return
          Endif
          size_cluster = size_cluster + 2
      
          nsigma%f(j1, j1t) = sigma_j2
          nsigma%f(j2, j2t) = sigma_j1
          Geo_cluster(j1, j1t) = .true.
          Geo_cluster(j2, j2t) = .true.
      
          Do i = 1,4
            k1 = Ising_nnlist(j1,i)
            k2 = Find_geo_partner_space(k1)
            call Geo_cluster_tryadd(sigma_j1, k1, j1t, k2, j2t, .false., size_cluster, Geo_cluster)
          enddo
      
          k1t = j1t + 1
          if (k1t > Ltrot) k1t = 1
          k2t = j2t - 1
          if (k2t < 1) k2t = Ltrot
          call Geo_cluster_tryadd(sigma_j1, j1, k1t, j2, k2t, .true., size_cluster, Geo_cluster)
      
          k1t = j1t - 1
          if (k1t < 1) k1t = Ltrot
          k2t = j2t + 1
          if (k2t > Ltrot) k2t =1
          call Geo_cluster_tryadd(sigma_j1, j1, k1t, j2, k2t, .true., size_cluster, Geo_cluster)
      
        End Subroutine Geo_cluster_tryadd
      !========================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(nsigma_old)
      
          !>  This function computes the ratio:  e^{-S0(nsigma%f)}/e^{-S0(nsigma_old)}
          Implicit none
      
          !> Arguments
          type(fields), intent(in) :: nsigma_old
      !     Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
          !> Local
          Integer :: I, nt, nt1, I1, I2, nc_J, nc_h_p, nc_h_m, N
      
          Delta_S0_global = 1.D0
          nc_J = 0
          nc_h_p = 0
          nc_h_m = 0
      
          select case( Model_vers )
          case( 0,1 )
            N = latt%N
          case( 3 )
            N = latt%N*2
          case( 2 )
            N = latt%N
          case default
            Write(6,*) "Error in Delta_S0_global: Model not yet implemented!"
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end select
      
          Do I = 1,N
            If ( I > latt%N ) then
              I1 = latt%nnlist(I - latt%N,1,0) + latt%N
              I2 = latt%nnlist(I - latt%N,0,1) + latt%N
            else
              I1 = latt%nnlist(I,1,0)
              I2 = latt%nnlist(I,0,1)
            endif
            Do nt = 1,Ltrot
                nt1 = nt + 1
                if (nt == Ltrot) nt1 = 1
                if (nsigma%i(I,nt) == nsigma%i(I,nt1) ) then
                  nc_h_p = nc_h_p + 1
                else
                  nc_h_m = nc_h_m + 1
                endif
                if (nsigma_old%i(I,nt) == nsigma_old%i(I,nt1) ) then
                  nc_h_p = nc_h_p - 1
                else
                  nc_h_m = nc_h_m - 1
                endif
      
                nc_J = nc_J + nsigma%i(I,nt)*nsigma%i(I1,nt) &
                    &      + nsigma%i(I,nt)*nsigma%i(I2,nt) &
                    &      - nsigma_old%i(I,nt)*nsigma_old%i(I1,nt) &
                    &      - nsigma_old%i(I,nt)*nsigma_old%i(I2,nt)
            enddo
          enddo
      
          Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) &
                  &         * exp( Dtau * Ham_J*real(nc_J,kind(0.d0)))
      
        end Function Delta_S0_global
      
!========================================================================


        function E_kin(GRC)
          Implicit none
          Complex (Kind=Kind(0.d0)), intent(in) :: GRC(Ndim,Ndim,N_FL)
       
          Complex (Kind=Kind(0.d0)) :: E_kin
       
          Integer :: n, nf, I, J
       
          E_kin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do n  = 1,Size(Op_T,1)
             Do nf = 1,N_FL
                Do I = 1,Size(Op_T(n,nf)%O,1)
                   Do J = 1,Size(Op_T(n,nf)%O,2)
                      E_kin = E_kin + Op_T(n,nf)%O(i, j)*Grc( Op_T(n,nf)%P(I), Op_T(n,nf)%P(J), nf )
                   ENddo
                Enddo
             Enddo
          Enddo
          E_kin = E_kin * dble(N_SUN)
       end function E_kin
       
       
       Subroutine Obser(GR,Phase,Ntau, Mc_step_weight)
         Implicit none
       
         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
         Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
         Integer, INTENT(IN)          :: Ntau
         Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
         !Local
         Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZP, ZS, Z
       
         Complex (Kind=Kind(0.d0)) :: Z_z_ising1, Z_x_ising1, Z_m1, Z_chi
         Complex (Kind=Kind(0.d0)) :: Z_z_ising2, Z_x_ising2, Z_m2
         Complex (Kind=Kind(0.d0)) :: Zkin, Zpot, Zrho
         Complex (Kind=Kind(0.d0)) :: d_Z, d_X
         Integer :: nf, I, J, I1, nc1, imj, J1, no_I, no_J, nt1, nt, dnt, Ntau1, N_ising, n
       
         ZP = PHASE/Real(Phase, kind(0.D0))
         ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
       
       !   Write(6,*) "Ntau =", Ntau
       !   Write(6,*) "i_sweep =", i_sweep, ZP, ZS

         if (Ntau == 1) then
            if (measure_hist_on_next_ntau1) then
              CALL measure_hist(PHASE)
              measure_hist_on_next_ntau1 = .false.
            else
              measure_hist_on_next_ntau1 = .true.
            endif
         endif
       
         Do nf = 1,N_FL
           Do I = 1,Ndim
             Do J = 1,Ndim
               GRC(I, J, nf) = -GR(J, I, nf)
             Enddo
             GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
           Enddo
         Enddo
         ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >
       
         Z_z_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_x_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_m1       = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_z_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_x_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_m2       = cmplx(0.d0, 0.d0, kind(0.D0))
       
         nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
       
         Zkin = E_kin(GRC)
       
         Zpot = cmplx(0.d0, 0.d0, kind(0.D0))
         do I = 1, size(Op_V, 1)
           I1 = Op_V(I,1)%P(1)
           do nc1 = 2, size( Op_V(I,1)%O, 1 )
             Zpot = Zpot + nsigma%f(I,Ntau) * GRC(I1, Op_V(I,1)%P(nc1) ,1) * Op_V(I,1)%O(1,nc1)
             Zpot = Zpot + nsigma%f(I,Ntau) * GRC(Op_V(I,1)%P(nc1), I1 ,1) * Op_V(I,1)%O(nc1,1)
           enddo
         enddo
         Zpot = Zpot * dble(N_SUN) * ham_t * ham_xi
       
         call Obs_scal(3)%measure( [Zkin, Zpot, Zkin+Zpot], Phase )
       
         Zrho = cmplx(0.d0, 0.d0, kind(0.D0))
         do I = 1, Ndim
           Zrho = Zrho + Grc(i, i, 1)
         enddo
         Zrho = Zrho * dble(N_SUN)
       
         call Obs_scal(1)%measure( [Zrho], Phase )
       
         do I = 1,Latt%N
           Z_z_ising1 = Z_z_ising1 + nsigma%f(I,Ntau)
           Z_m1 = Z_m1 + nsigma%f(1,Ntau)*nsigma%f(I,Ntau)
       
           if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
             Z_x_ising1 = Z_x_ising1 + eq_x_ising
           else
             Z_x_ising1 = Z_x_ising1 + neq_x_ising
           endif
         enddo
         Z_z_ising1 = Z_z_ising1/Latt%N
         Z_x_ising1 = Z_x_ising1/Latt%N
         Z_m1 = Z_m1/Latt%N
       
         call Obs_scal(2)%measure( [Z_z_ising1], Phase )
         call Obs_scal(4)%measure( [Z_x_ising1], Phase )
         call Obs_scal(5)%measure( [Z_m1, Z_m1**2, Z_m1**4], Phase )
       
         If ( Model_vers == 3 ) then
           do I = 1+Latt%N, 2*Latt%N
             Z_z_ising2 = Z_z_ising2 + nsigma%f(I,Ntau)
             Z_m2 = Z_m2 + nsigma%f(1+Latt%N,Ntau)*nsigma%f(I,Ntau)
       
             if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
               Z_x_ising2 = Z_x_ising2 + eq_x_ising
             else
               Z_x_ising2 = Z_x_ising2 + neq_x_ising
             endif
           enddo
           Z_z_ising2 = Z_z_ising2/Latt%N
           Z_x_ising2 = Z_x_ising2/Latt%N
           Z_m2 = Z_m2/Latt%N
       
           call Obs_scal(2)%measure( [Z_z_ising2], Phase )
           call Obs_scal(4)%measure( [Z_x_ising2], Phase )
           call Obs_scal(5)%measure( [Z_m2, Z_m2**2, Z_m2**4], Phase )
         endif
       
       
         ! counting up correlation functions
       !   DO I = 1,Size(Obs_eq,1)
         DO I = 1,3
           Obs_eq(I)%N        = Obs_eq(I)%N + 1
           Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
         ENDDO
       
         ! Compute Ising X-X and Z-Z correlation functions
         nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
         Do I = 1,Latt%N
           do no_I = 1, size(Obs_eq(1)%Obs_Latt0,1)
             I1 = I + (no_I-1)*Latt%N
             if ( nsigma%f(I1,nt1) == nsigma%f(I1,Ntau) ) then
               Z_x_ising1 = eq_x_ising
             else
               Z_x_ising1 = neq_x_ising
             endif
             Obs_eq(1)%Obs_Latt0(no_I) = Obs_eq(1)%Obs_Latt0(no_I) + Z_x_ising1      * ZP*ZS
         !     Obs_eq(2)%Obs_Latt0(1) = Obs_eq(2)%Obs_Latt0(1) + nsigma%f(I,Ntau) * ZP*ZS
             Do J = 1,Latt%N
               do no_J = 1, size(Obs_eq(1)%Obs_Latt0,1)
                 J1 = J + (no_J-1)*Latt%N
                 imj = Latt%imj(I,J)
                 Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + nsigma%f(I1,Ntau) * nsigma%f(J1,Ntau) * ZP*ZS
                 if ( nsigma%f(J1,nt1) == nsigma%f(J1,Ntau) ) then
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + Z_x_ising1 * eq_x_ising  * ZP*ZS
                 else
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + Z_x_ising1 * neq_x_ising * ZP*ZS
                 endif
               enddo
             enddo
           enddo
         enddo
       
       !   If ( Model_vers == 3 ) then
       !     Do I = 1+Latt%N,2*Latt%N
       !       if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
       !         Z_x_ising1 = eq_x_ising
       !       else
       !         Z_x_ising1 = neq_x_ising
       !       endif
       !       Obs_eq(1)%Obs_Latt0(1) = Obs_eq(1)%Obs_Latt0(1) + Z_x_ising1      * ZP*ZS
       ! !       Obs_eq(2)%Obs_Latt0(1) = Obs_eq(2)%Obs_Latt0(1) + nsigma%f(I,Ntau) * ZP*ZS
       !       Do J = 1+Latt%N,2*Latt%N
       !         imj = Latt%imj(I-Latt%N,J-Latt%N)
       !         Obs_eq(2)%Obs_Latt(imj,1,1,1) = Obs_eq(2)%Obs_Latt(imj,1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,Ntau) * ZP*ZS
       !         if ( nsigma%f(J,nt1) == nsigma%f(J,Ntau) ) then
       !           Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising1 * eq_x_ising  * ZP*ZS
       !         else
       !           Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising1 * neq_x_ising * ZP*ZS
       !         endif
       !       enddo
       !     enddo
       !   endif
       
         ! Computing time-displaced X-X and Z-Z correlation functions and chi
       
         nBlub2 = nBlub2 + 1
         if ( nBlub2 > 48 ) then
           nBlub2 = 0
       
         DO I = 4,5
           Obs_eq(I)%N        = Obs_eq(I)%N + 1
           Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
         ENDDO
       
         Z_chi     = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_x_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
         Ntau1 = Ntau+1; if ( Ntau == Ltrot ) Ntau1 = 1
       
         do dnt= 0, Ltrot
           nt = Ntau + dnt
           if ( nt > Ltrot ) nt = nt - Ltrot
           nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
           do I = 1,Latt%N
             if ( nsigma%f(I,nt1) == nsigma%f(I,nt) ) then
               Z_x_ising1 = eq_x_ising
             else
               Z_x_ising1 = neq_x_ising
             endif
             Z_x_ising1 = Z_x_ising1 + Z_x_ising1
             Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + Z_x_ising1    * ZP*ZS
       !       Obs_eq(4)%Obs_Latt0(1) = Obs_eq(4)%Obs_Latt0(1) + nsigma%f(I,nt) * ZP*ZS
             do J = 1,Latt%N
               imj = Latt%imj(I,J)
               Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,nt) * ZP*ZS
               if ( nt == Ntau .and. I == J ) then
                 Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + 1  * ZP*ZS
                 Z_chi = Z_chi + 1
               elseif ( nsigma%f(J,Ntau1) == nsigma%f(J,Ntau) ) then
                 Z_chi = Z_chi + Z_x_ising1 * eq_x_ising
                 Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * eq_x_ising  * ZP*ZS
               else
                 Z_chi = Z_chi + Z_x_ising1 * neq_x_ising
                 Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * neq_x_ising * ZP*ZS
               endif
             enddo
           enddo
         enddo
       
         If ( Model_vers == 3 ) then
           do dnt= 0, Ltrot
             nt = Ntau + dnt
             if ( nt > Ltrot ) nt = nt - Ltrot
             nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
             do I = 1+Latt%N,2*Latt%N
               if ( nsigma%f(I,Ntau1) == nsigma%f(I,Ntau) ) then
                 Z_x_ising1 = eq_x_ising
               else
                 Z_x_ising1 = neq_x_ising
               endif
               Z_x_ising1 = Z_x_ising1 + Z_x_ising1
               Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + Z_x_ising1    * ZP*ZS
       !         Obs_eq(4)%Obs_Latt0(1) = Obs_eq(4)%Obs_Latt0(1) + nsigma%f(I,nt) * ZP*ZS
               do J = 1+Latt%N,2*Latt%N
                 imj = Latt%imj(I-Latt%N,J-Latt%N)
                 Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,nt) * ZP*ZS
                 if ( nt == Ntau .and. I == J ) then
                   Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + 1  * ZP*ZS
                   Z_chi = Z_chi + 1
                 elseif ( nsigma%f(J,nt) == nsigma%f(J,Ntau) ) then
                   Z_chi = Z_chi + Z_x_ising1 * eq_x_ising
                   Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * eq_x_ising  * ZP*ZS
                 else
                   Z_chi = Z_chi + Z_x_ising1 * neq_x_ising
                   Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * neq_x_ising * ZP*ZS
                 endif
               enddo
             enddo
           enddo
         endif
       
         If ( Model_vers == 3 ) then
           N_ising = 2*Latt%N
         else
           N_ising = Latt%N
         endif
       
         call Obs_scal(6)%measure( [ Z_chi/(N_ising**2), Z_x_ising1/N_ising/Ltrot ], Phase )
       
         endif
       
         ! Compute Green-function
         Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
         Do I1 = 1,Ndim
           I    = List(I1,1)
           no_I = List(I1,2)
           Do J1 = 1,Ndim
             J    = List(J1,1)
             no_J = List(J1,2)
             imj = Latt%imj(I,J)
             Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + Z * GRC(I1,J1,1) *  ZP*ZS
           enddo
         enddo
       
       end Subroutine Obser
       
       
       Subroutine measure_hist(Phase)
         Implicit none
         Complex (Kind=Kind(0.d0)), Intent(IN)    :: Phase
       
         Complex (Kind=Kind(0.d0)) :: ZP, ZS
       
         Complex (Kind=Kind(0.d0)) :: Z_z_ising1, Z_x_ising1, Z_m1, Z_m1_sq, Z_m1_quad
         Complex (Kind=Kind(0.d0)) :: Z_z_ising2, Z_x_ising2, Z_m2, Z_m2_sq, Z_m2_quad
         Complex (Kind=Kind(0.d0)) :: Z_m_temp
         real    (Kind=Kind(0.d0)) :: B
         Complex (Kind=Kind(0.d0)) :: d_Z, d_X
         Integer :: I, nt1, nt
       
       
         ZP = PHASE/Real(Phase, kind(0.D0))
         ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
       
       
         Z_z_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_x_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_m1       = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_m1_sq    = cmplx(0.d0, 0.d0, kind(0.D0))
         Z_m1_quad  = cmplx(0.d0, 0.d0, kind(0.D0))
       
         do nt = 1, Ltrot
           nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
           Z_m_temp = cmplx(0.d0, 0.d0, kind(0.D0))
           do I = 1,Latt%N
             Z_z_ising1 = Z_z_ising1 + nsigma%f(I,nt)
             Z_m_temp = Z_m_temp + nsigma%f(1,nt)*nsigma%f(I,nt)
       
             if ( nsigma%f(I,nt1) == nsigma%f(I,nt) ) then
               Z_x_ising1 = Z_x_ising1 + eq_x_ising
             else
               Z_x_ising1 = Z_x_ising1 + neq_x_ising
             endif
           enddo
           Z_m_temp = Z_m_temp/Latt%N
           Z_m1 = Z_m1 + Z_m_temp
           Z_m1_sq = Z_m1_sq + Z_m_temp**2
           Z_m1_quad = Z_m1_quad + Z_m_temp**4
         enddo
       
         Z_z_ising1 = Z_z_ising1/(Latt%N*Ltrot)
         Z_x_ising1 = Z_x_ising1/(Latt%N*Ltrot)
         Z_m1 = Z_m1/Ltrot
         Z_m1_sq = Z_m1_sq/Ltrot
         Z_m1_quad = Z_m1_quad/Ltrot
       
         call Obs_scal(7)%measure( [Z_z_ising1], Phase )
         call Obs_scal(8)%measure( [Z_x_ising1], Phase )
         call Obs_scal(9)%measure( [Z_m1, Z_m1**2, Z_m1**4], Phase )
       
         call Obs_hist(1)%measure( dble(Z_x_ising1), dble(ZS) )
         call Obs_hist(2)%measure( dble(Z_m1), dble(ZS) )
         B = (3.d0 - dble(Z_m1_quad/Z_m1_sq**2))/2.d0
         call Obs_hist(3)%measure( B, dble(ZS) )
       
         If ( Model_vers == 3 ) then
           Z_z_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
           Z_x_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
           Z_m2       = cmplx(0.d0, 0.d0, kind(0.D0))
           Z_m2_sq    = cmplx(0.d0, 0.d0, kind(0.D0))
           Z_m2_quad  = cmplx(0.d0, 0.d0, kind(0.D0))
       
           do nt = 1, Ltrot
             nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
             Z_m_temp = cmplx(0.d0, 0.d0, kind(0.D0))
             do I = 1+Latt%N, 2*Latt%N
               Z_z_ising2 = Z_z_ising2 + nsigma%f(I,nt)
               Z_m_temp = Z_m_temp + nsigma%f(1+Latt%N,nt)*nsigma%f(I,nt)
       
               if ( nsigma%f(I,nt1) == nsigma%f(I,nt) ) then
                 Z_x_ising2 = Z_x_ising2 + eq_x_ising
               else
                 Z_x_ising2 = Z_x_ising2 + neq_x_ising
               endif
             enddo
             Z_m_temp = Z_m_temp/Latt%N
             Z_m2 = Z_m2 + Z_m_temp
             Z_m2_sq = Z_m2_sq + Z_m_temp**2
             Z_m2_quad = Z_m2_quad + Z_m_temp**4
           enddo
           
           Z_z_ising2 = Z_z_ising2/(Latt%N*Ltrot)
           Z_x_ising2 = Z_x_ising2/(Latt%N*Ltrot)
           Z_m2 = Z_m2/Ltrot
           Z_m2_sq = Z_m2_sq/Ltrot
           Z_m2_quad = Z_m2_quad/Ltrot
       
           d_Z = (Z_z_ising1 - Z_z_ising2)**2
           d_X = (Z_x_ising1 - Z_x_ising2)**2
       
           call Obs_scal(7)%measure( [Z_z_ising2], Phase )
           call Obs_scal(8)%measure( [Z_x_ising2], Phase )
           call Obs_scal(9)%measure( [Z_m2, Z_m2**2, Z_m2**4], Phase )
       
           call Obs_scal(10)%measure( [d_Z], Phase )
           call Obs_scal(11)%measure( [d_X], Phase )
       
           call Obs_hist(1)%measure( dble(Z_x_ising2), dble(ZS) )
           call Obs_hist(2)%measure( dble(Z_m2), dble(ZS) )
           call Obs_hist(3)%measure( (3.d0 - dble(Z_m2_quad/Z_m2_sq**2))/2.d0, dble(ZS) )
       
           call Obs_hist(4)%measure( dble(d_Z), dble(ZS) )
           call Obs_hist(5)%measure( dble(d_X), dble(ZS) )
         endif
       
       end Subroutine measure_hist
       
!=====================================================
       Subroutine ObserT(NT, GT0,G0T,G00,GTT, PHASE, Mc_step_weight)
        Implicit none

        Integer         , INTENT(IN) :: NT
        Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
        Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
        Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

        !Locals
        Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
        Integer :: IMJ, I, J, I1, J1, no_I, no_J

        real (Kind=Kind(0.d0)) :: Z_z_ising
        Integer                   :: Ntau

!           Write(6,*) "ObserT NT =" , NT

        Z_z_ising = 0.d0
        do I = 1,Latt%N
          do Ntau = 1, Ltrot
            Z_z_ising = Z_z_ising + real(nsigma%f(I,Ntau))
          enddo
        enddo

        ZP = PHASE/Real(Phase, kind(0.D0))
        ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
        If (NT == 0 ) then
           Obs_tau(1)%N = Obs_tau(1)%N + 1
           Obs_tau(1)%Ave_sign = Obs_tau(1)%Ave_sign + Real(ZS,kind(0.d0))
!              if ( Z_z_ising > 0 ) then
!                 Obs_tau(2)%N = Obs_tau(2)%N + 1
!                 Obs_tau(2)%Ave_sign = Obs_tau(2)%Ave_sign + Real(ZS,kind(0.d0))
!              else
!                 Obs_tau(3)%N = Obs_tau(3)%N + 1
!                 Obs_tau(3)%Ave_sign = Obs_tau(3)%Ave_sign + Real(ZS,kind(0.d0))
!              endif
        endif
           Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
           Do I1 = 1,Ndim
              I    = List(I1,1)
              no_I = List(I1,2)
              Do J1 = 1,Ndim
                 J    = List(J1,1)
                 no_J = List(J1,2)
                 imj = latt%imj(I,J)
                 ! Green
                 Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                      & +  Z * GT0(I1,J1,1) * ZP* ZS
!                    if ( Z_z_ising > 0 ) then
!                       Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                            & +  Z * GT0(I1,J1,1) * ZP* ZS
!                    else
!                       Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                            & +  Z * GT0(I1,J1,1) * ZP* ZS
!                    endif
              Enddo
           Enddo

      end Subroutine OBSERT

!==========================================================

   end submodule Ham_Nematic_Dirac_mod
