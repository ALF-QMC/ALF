!  Copyright (C) 2016 - 2023 The ALF project
!
!  This file is part of the ALF project.
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
!>
!> @brief
!> Handles Hubbard Stratonovitch fields.
!>
!> @details
!> A general operator has the form: \f$ \gamma_{n,\tau} e^{ \phi_{n,\tau} g \hat{O}_{n,\tau} }  \f$.
!>
!> For  type=1 the fields, f, take two  integer values, \f$\pm 1 \f$ and  \f$ \gamma_{n,\tau}(f) = 1,  \phi_{n,\tau}(f) = f \f$
!>
!> For  type=2 the fields, f, take four integer values \f$\pm 1, \pm 2 \f$ and
!>     \f[ \gamma_{n,\tau}(\pm 1)  = 1 + \sqrt{6}/3,
!>      \gamma_{n,\tau}(\pm 2)  = 1 - \sqrt{6}/3,  \phi_{n,\tau}(\pm 1) = \pm \sqrt{2  ( 3 - \sqrt{6} ) },
!>       \phi_{n,\tau}(\pm 2) = \pm \sqrt{2  ( 3 + \sqrt{6} ) }  \f]
!> For  type=3 the fields, f, are real and  \f$ \gamma_{n,\tau}(f)  = 1, \phi_{n,\tau}(f) = f \f$
!>
!> For  type=4   is  for  two  HS  fields  per  vertex,  encoded in the  real  and imaginary parts of f.
!>                Re(f) = \pm 1, pm 2,  and   Im(f) is  real    with    
!>                 gamma = gamma(Re(f))   and   phi = \sqrt{1 + aimag(f) } eta(real(f))  
!--------------------------------------------------------------------


     Module Fields_mod
       
#ifdef MPI
       Use mpi
#endif
#if defined HDF5
       Use hdf5
       use h5lt
#endif
       Use runtime_error_mod
       Use Random_Wrap
       use iso_fortran_env, only: output_unit, error_unit

       Public Fields
       Public Fields_init

       Private
       Real (Kind=Kind(0.d0))  :: Phi_st(-2:2,2),  Gama_st(-2:2,2)
       Real (Kind=Kind(0.d0))  :: FLIP_st(-2:2,3)
       Real (Kind=Kind(0.d0))  :: Amplitude

       Type Fields
          Complex (Kind=Kind(0.d0)), allocatable    :: f(:,:)             ! Field
          Integer                  , allocatable    :: t(:)               ! Type
          Integer                  , allocatable    :: Flip_protocol(:)   ! Flip_protocol   Allows  for different  flip  protocols.
                                                                          ! Used only  for  type  4.
        CONTAINS
          procedure  :: make  => Fields_make
          procedure  :: clear => Fields_clear
          procedure  :: set   => Fields_set
          procedure  :: out   => Fields_out
          procedure  :: in    => Fields_in
          procedure  :: i     => Fields_get_i
          procedure  :: Phi   => Fields_Phi
          procedure  :: Gama  => Fields_Gama
          procedure  :: Flip  => Fields_Flip
          procedure, private  :: read_conf    => Fields_read_conf
#if defined HDF5
          procedure, private  :: read_conf_h5 => Fields_read_conf_h5
#endif
          procedure, private  :: write_conf   => Fields_write_conf
       END TYPE Fields

    Contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Returns Phi of the field this(n_op,n_tau)
!>
!-------------------------------------------------------------------

      Complex (Kind=Kind(0.d0)) function  Fields_Phi(this,n_op,n_tau)

        Implicit none
        Class (Fields) :: this
        Integer, INTENT(IN) ::  n_op, n_tau


        select case (this%t(n_op))
        case(1)
           Fields_Phi = cmplx(Phi_st(Nint(real(this%f(n_op,n_tau))),1), 0.d0,kind(0.d0))
        case(2)
           Fields_Phi = cmplx(Phi_st(Nint(real(this%f(n_op,n_tau))),2), 0.d0,kind(0.d0))
        case(3)
           Fields_Phi = cmplx(real(this%f(n_op,n_tau),kind(0.d0))     , 0.d0,kind(0.d0))
        case(4)
           Fields_Phi = cmplx(Phi_st(Nint(real(this%f(n_op,n_tau))),2),0.d0,kind(0.d0)) * &
                &       sqrt(cmplx( 1.d0 +  aimag(this%f(n_op,n_tau)), 0.d0,kind(0.d0)) )
        case default
           Write(error_unit,*) 'Error in Fields_Phi'
           CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
        end select
      end function Fields_Phi

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Returns Gamma of the field this(n_op,n_tau)
!>
!-------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function Fields_Gama(this,n_op,n_tau)

        Implicit none
        Class (Fields) :: this
        Integer, INTENT(IN) ::  n_op, n_tau

        select case (this%t(n_op))
        case(1)
           Fields_GAMA = 1.d0
        case(2)
           Fields_GAMA = GAMA_st(Nint(Real(this%f(n_op,n_tau))),2)
        case(3)
           Fields_GAMA = 1.d0
        case(4)
           Fields_GAMA = GAMA_st(Nint(Real(this%f(n_op,n_tau))),2)
        case default
           Write(error_unit,*) 'Error in Fields_GAMA'
           CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
        end select

      end function Fields_Gama

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Flips the field this(n_op,n_tau)
!>
!-------------------------------------------------------------------
      Complex (Kind=Kind(0.d0)) function  Fields_flip(this,n_op,n_tau)

        Implicit none
        Class (Fields)      :: this
        Integer, INTENT(IN) :: n_op, n_tau

       
        select case (this%t(n_op))
        case(1)
           Fields_flip = - this%f(n_op,n_tau)
        case (2)
           Fields_flip =   cmplx(Flip_st( nint(real(this%f(n_op,n_tau))),nranf(3)),  0.d0,  Kind(0.d0) )
        case (3)
           Fields_flip =   cmplx(real(this%f(n_op,n_tau)) + Amplitude*( ranf_wrap() - 0.5D0), 0.d0,  Kind(0.d0))
        case (4)
           Select case (this%Flip_Protocol(n_op))
           case(1) ! Flip one of  the  two  fields randomly,  this is the  default
              If (ranf_wrap() > 0.5D0 )  then 
                 Fields_flip =   cmplx( Flip_st( nint(real(this%f(n_op,n_tau))),nranf(3))             , &
                      &                 aimag(this%f(n_op,n_tau)), Kind(0.d0))
              else
                 Fields_flip =   cmplx( real(this%f(n_op,n_tau)) , &
                      &                 aimag(this%f(n_op,n_tau)) +  Amplitude*( ranf_wrap() - 0.5D0) , Kind(0.d0))
              endif
           case(2) ! Flip both  fields
              Fields_flip =   cmplx( Flip_st( nint(real(this%f(n_op,n_tau))),nranf(3))                , &
                   &                 aimag(this%f(n_op,n_tau)) +  Amplitude*( ranf_wrap() - 0.5D0) , Kind(0.d0))
              
           case(3) ! Flip only  the  real  part
              Fields_flip =   cmplx( Flip_st( nint(real(this%f(n_op,n_tau))),nranf(3))             , &
                   &                 aimag(this%f(n_op,n_tau)), Kind(0.d0))
           case(4) ! Flip only  the   impaginary  part
              Fields_flip =   cmplx( real(this%f(n_op,n_tau)) , &
                   &                 aimag(this%f(n_op,n_tau)) +  Amplitude*( ranf_wrap() - 0.5D0) , Kind(0.d0))
           case default
              Write(error_unit,*) 'No flip protocol provided for  field. '
              CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
           end Select
           
        case default
           Write(error_unit,*) 'Error in Fields. '
           CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
        end select
        
      end function Fields_flip
      
!-------------------------------------------------------------------

      Integer function Fields_get_i(this,n_op,n_tau)

        Implicit none
        Class (Fields) :: this
        Integer, INTENT(IN) ::  n_op, n_tau

        if ( this%t(n_op) == 1 .or.   this%t(n_op) == 2 ) then
           Fields_get_i = NINT(Real(this%f(n_op,n_tau)))
        else
           Write(error_unit,*) "Error in fields"
           CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
        endif

      end function Fields_get_i
!-------------------------------------------------------------------
      Subroutine Fields_make(this,N_OP,N_tau)
        Implicit none
        Class (Fields), INTENT(INOUT)  :: this
        Integer, INTENT(IN)            :: N_OP, N_tau

        !Write(6,*) "Allocating  fields: ", N_op, N_tau
        allocate (this%f(N_OP,N_tau), this%t(N_OP), this%Flip_protocol(N_OP) )

        
        this%f = cmplx(0.d0,0.d0,kind(0.d0)) ;  this%t = 0; this%flip_protocol = 1

      end Subroutine Fields_make
!-------------------------------------------------------------------
      Subroutine Fields_clear(this)
        Implicit none
        Class (Fields) :: this

        deallocate (this%f, this%t,  this%Flip_protocol )
      end Subroutine Fields_clear

!-------------------------------------------------------------------

      Subroutine Fields_init(Amplitude_in)

        Implicit none

        Real  (Kind=Kind(0.d0)), Optional, Intent(IN) :: Amplitude_in

        !Local
        Integer :: n

        Amplitude = 1.d0
        If (Present(Amplitude_in)) Amplitude = Amplitude_in

        Phi_st = 0.d0
        do n = -2,2
           Phi_st(n,1) = real(n,Kind=Kind(0.d0))
        enddo
        Phi_st(-2,2) = - SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
        Phi_st(-1,2) = - SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
        Phi_st( 1,2) =   SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
        Phi_st( 2,2) =   SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )

        Do n = -2,2
           gama_st(n,1) = 1.d0
        Enddo
        GAMA_st(-2,2) = 1.D0 - SQRT(6.D0)/3.D0
        GAMA_st( 2,2) = 1.D0 - SQRT(6.D0)/3.D0
        GAMA_st(-1,2) = 1.D0 + SQRT(6.D0)/3.D0
        GAMA_st( 1,2) = 1.D0 + SQRT(6.D0)/3.D0

        FLIP_st(-2,1) = -1.d0
        FLIP_st(-2,2) =  1.d0
        FLIP_st(-2,3) =  2.d0

        FLIP_st(-1,1) =  1.d0
        FLIP_st(-1,2) =  2.d0
        FLIP_st(-1,3) = -2.d0

        FLIP_st( 1,1) =  2.d0
        FLIP_st( 1,2) = -2.d0
        FLIP_st( 1,3) = -1.d0

        FLIP_st( 2,1) = -2.d0
        FLIP_st( 2,2) = -1.d0
        FLIP_st( 2,3) =  1.d0

      end Subroutine Fields_init



!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Reads in field configuration
!>
!> @details
!> Reads in the field configuration and seeds if present so as to
!> pursue a run. If  the configuration is not present  the
!> routine will generate one randomly. Note that the random number generator is
!> initialized from the  seeds file in the routine Set_Random_number_Generator  of the
!> module random_wrap_mod.F90
!>
!> @param [INOUT] this
!> \verbatim
!> Type Fields
!> On input  test%t(:) is set. The operator types are time independent.
!> On output test%f(:,:) is initialized \endverbatim
!>
!> @param [IN] Group_Comm
!> \verbatim
!> Type Integer
!> Communicator for MPI \endverbatim
!>
!> @param [Optional]  Initial_field
!> \verbatim
!> Type Complex
!> Initial field \endverbatim
!--------------------------------------------------------------------
      Subroutine Fields_in(this,Group_Comm,Initial_field)

        Implicit none

        Class (Fields)        , INTENT(INOUT) :: this
        Integer               , INTENT(IN   ) :: Group_Comm
        Complex (Kind=Kind(0.d0)), Dimension(:,:), Optional   :: Initial_field

        ! LOCAL
        Integer                 :: I, I1, IERR, SEED_IN, K, NT
        Real (Kind=Kind(0.d0) ) :: X
        Integer, DIMENSION(:), ALLOCATABLE :: SEED_VEC
        Logical ::   LCONF, LCONF_H5
        Character (LEN=64) :: FILE_SR, FILE_TG, FILE_seeds, FILE_info, File1, FILE_TG_H5, File1_h5

#ifdef MPI
        INTEGER        :: STATUS(MPI_STATUS_SIZE), irank_g, isize_g, igroup, ISIZE, IRANK
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif


#if defined(MPI)
#if defined(TEMPERING)
            write(FILE1,'(A,I0,A)')      "Temp_",igroup,"/confin_0"
            write(FILE_TG,'(A,I0,A,I0)') "Temp_",igroup,"/confin_",irank_g
            write(FILE_info,'(A,I0,A)')  "Temp_",igroup,"/info"
#else
            File1 = "confin_0"
            write(FILE_TG,'(A,I0)') "confin_",irank_g
            FILE_info="info"
#endif
#else
            File1   = "confin_0"
            FILE_TG = "confin_0"
            FILE_info="info"
#endif
            FILE_seeds="seeds"
            write(FILE1_H5,  '(A,A)') trim(FILE1)  , ".h5"
            write(FILE_TG_H5,'(A,A)') trim(FILE_TG), ".h5"

            INQUIRE (FILE=File1, EXIST=LCONF)
            INQUIRE (FILE=File1_h5, EXIST=LCONF_H5)
            IF (LCONF .and. LCONF_H5) THEN
               write(error_unit,*) "ERROR: Both plain text configuration file confin_0 and HDF5 file confin_0.h5 exists."
               write(error_unit,*) "   Delete one of them to remove ambiguity, program aborted!"
               CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            ENDIF
#if defined HDF5
            CALL h5open_f(ierr)
            IF (LCONF_H5) THEN
               CALL this%read_conf_h5(FILE_TG_H5)
            elseif (LCONF) THEN
               write(error_unit,*) "WARNING: READING plain text configuration file confin_0 even though HDF5 support is compiled in."
               write(error_unit,*) "   This feature is only for switching from plain text to HDF5, please do not mix up both formats!"
               CALL this%read_conf(FILE_TG)
            ELSE
#else
            IF (LCONF_H5) THEN
               write(error_unit,*) "ERROR: HDF5 configuration file confin_0.h5 exists, even though program is compiled"
               write(error_unit,*) "   without HDF5! You cannot mix up HDF5 runs with non-HDF5 runs, program aborted!"
               CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            ENDIF
            IF (LCONF) THEN
               CALL this%read_conf(FILE_TG)
            ELSE
#endif
#if defined MPI
               IF (IRANK == 0) THEN
#endif
                  WRITE(6,*) 'No initial configuration'
                  OPEN(UNIT=5,FILE=FILE_seeds,STATUS='OLD',ACTION='READ',IOSTAT=IERR)
                  IF (IERR /= 0) THEN
                     WRITE(error_unit,*) 'Fields_in: unable to open <seeds>',IERR
                     CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
                  END IF
#if defined MPI
                  DO I = ISIZE-1,1,-1
                     READ (5,*) SEED_IN
                     CALL MPI_SEND(SEED_IN,1,MPI_INTEGER, I, I+1024, MPI_COMM_WORLD,IERR)
                  ENDDO
                  READ(5,*) SEED_IN
                  CLOSE(5)
               ELSE
                  CALL MPI_RECV(SEED_IN, 1, MPI_INTEGER,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
               ENDIF
#else
               READ (5,*) SEED_IN
               CLOSE(5)
#endif
               ALLOCATE (SEED_VEC(1))
               SEED_VEC(1) = SEED_IN
               CALL RANSET(SEED_VEC)
               DEALLOCATE (SEED_VEC)
               If (Present(Initial_field)) then
                  this%f = Initial_field
               else
                  Call  this%set()
               endif
#if defined MPI
               if (irank_g == 0) then
#endif
                  Open (Unit = 50,file=FILE_info,status="unknown",position="append")
                  WRITE(50,*) 'No initial configuration, Seed_in', SEED_IN
                  Close(50)
#if defined MPI
               endif
#endif
            ENDIF
            call Fields_test(this)

       end Subroutine Fields_in
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Writes out the field configuration
!>
!> @details
!
!> @param [IN] this
!> \verbatim
!> Type Fields
!>
!> @param [IN] Group_Comm
!> \verbatim
!> Type Integer
!> Communicator for MPI \endverbatim
!>
!--------------------------------------------------------------------

       SUBROUTINE Fields_out(this,Group_Comm)

         IMPLICIT NONE

         Class (Fields), INTENT(INOUT) :: this
         Integer,        INTENT(IN   ) :: Group_Comm

         ! LOCAL
         CHARACTER (LEN=64) :: FILE_TG

#if defined(MPI)
         INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
         !Write(6,*) "Group, rank :", igroup, irank_g
#if defined(TEMPERING)
         write(FILE_TG,'(A,I0,A,I0)') "Temp_",igroup,"/confout_",irank_g
#else
         write(FILE_TG,'(A,I0)') "confout_",irank_g
#endif
#else
         FILE_TG = "confout_0"
#endif

#if defined HDF5
         write(FILE_TG,'(A,A)') trim(FILE_TG), ".h5"
#endif
         call this%write_conf(FILE_TG)

       END SUBROUTINE Fields_out


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Test if the field values are consistent with the field types.
!>
!> @details
!>
!> @param [INOUT] this
!> \verbatim
!> Type Fields
!> \endverbatim
!--------------------------------------------------------------------
       Subroutine  Fields_test(this)

         Implicit none

         Class (Fields), INTENT(INOUT) :: this

         Integer :: nt, I, I1

         !Write(6,*) "Fields_set", size(this%f,1), size(this%f,2)
         Do nt = 1,size(this%f,2)
            Do I = 1,size(this%f,1)
               Select Case(this%t(i))
               Case(1)
                  ! Field should be 1 or -1
                  if (abs(this%f(I,nt))-1.d0 > 1e-8 .or. abs(aimag(this%f(I,nt))) > 1e-8) then
                     write(error_unit, '("A","I0","A","I0","A","I0","A","I0")') &
                       & "Field I=", I, " nt=", nt, "value=", this%f(I,nt), " does not fit to type ", this%t(i)
                     CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
                  endif
               Case(2)
                  ! Field should be 1, -1, 2, or -2
                  if ((abs(this%f(I,nt))-1.d0 > 1e-8 .and. abs(this%f(I,nt))-2.d0 > 1e-8) .or. abs(aimag(this%f(I,nt))) > 1e-8) then
                     write(error_unit, '("A","I0","A","I0","A","I0","A","I0")') &
                       & "Field I=", I, " nt=", nt, "value=", this%f(I,nt), " does not fit to type ", this%t(i)
                     CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
                  endif
               Case(3)
                  ! Field should be real
                  if (abs(aimag(this%f(I,nt))) > 0) then
                     write(error_unit, '("A","I0","A","I0","A","I0","A","I0")') &
                       & "Field I=", I, " nt=", nt, "value=", this%f(I,nt), " does not fit to type ", this%t(i)
                     CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
                  endif
               Case(4)
                  ! Real part of Field should be 1, -1, 2, or -2
                  if (abs(real(this%f(I,nt)))-1.d0 > 1e-8 .and. abs(real(this%f(I,nt)))-2.d0 > 1e-8) then
                     write(error_unit, '("A","I0","A","I0","A","I0","A","I0")') &
                       & "Real part of field I=", I, " nt=", nt, "value=", this%f(I,nt), " does not fit to type ", this%t(i)
                     CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
                  endif
               Case default
                  write(error_unit, '("A","I0","A","I0")') "Field ", I, " has unrecongnized type ", this%t(i)
                  CALL Terminate_on_error(ERROR_FIELDS,__FILE__,__LINE__)
               end Select
            enddo
         enddo

       end Subroutine Fields_test


!--------------------------------------------------------------------
       !> @author
!> ALF-project
!>
!> @brief
!> Sets the field.
!>
!> @details
!>
!> @param [INOUT] this
!> \verbatim
!> Type Fields
!> On input the size if this%f is used test%t is set.
!> On output this%f is  initialized to a random configuration \endverbatim
!--------------------------------------------------------------------
       Subroutine  Fields_set(this)

         Implicit none

         Class (Fields), INTENT(INOUT) :: this

         Integer :: nt, I, I1

         !Write(6,*) "Fields_set", size(this%f,1), size(this%f,2)
         Do nt = 1,size(this%f,2)
            Do I = 1,size(this%f,1)
               if (this%t(i)  < 4 ) then
                  this%f(I,nt)  = cmplx(1.d0,0.d0,kind(0.d0))
                  if ( ranf_wrap() > 0.5D0 ) this%f(I,nt) = cmplx(-1.d0,0.d0,kind(0.d0))
               else 
                  I1 = 1
                  if ( ranf_wrap() > 0.5D0 ) I1 = -1
                  this%f(I,nt)  = cmplx(dble(I1),  Amplitude*( ranf_wrap() - 0.5D0) ,Kind(0.d0))
               endif
            enddo
         enddo

       end Subroutine Fields_set

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Reads in field configuration for single process, private Subroutine.
!>
!
!> @param [INOUT] this
!> \verbatim
!> Type Fields
!> On input  test%t(:) is set. The operator types are time independent.
!> On output test%f(:,:) is initialized \endverbatim
!>
!> @param [IN] filename
!> \verbatim
!> Type CHARACTER (LEN=64)
!> Name of file from which to read configuration and random seed \endverbatim
!--------------------------------------------------------------------
        SUBROUTINE Fields_read_conf(this, filename)
            IMPLICIT NONE

            Class (Fields)    , INTENT(INOUT) :: this
            CHARACTER (LEN=64), intent(in)    :: filename

            INTEGER                :: K, I, NT, I1
            INTEGER,ALLOCATABLE    :: SEED_VEC(:)
            Real (Kind=Kind(0.d0)) :: X

            CALL GET_SEED_LEN(K)
            ALLOCATE(SEED_VEC(K))

            OPEN (UNIT = 10, FILE=filename, STATUS='OLD', ACTION='READ')
            READ(10,*) SEED_VEC
            CALL RANSET(SEED_VEC)
            DO NT = 1,SIZE(this%f,2)
               DO I = 1,SIZE(this%f,1)
                  IF (this%t(I) == 1 .or.  this%t(I) == 2) then
                     Read(10,*)  I1
                     this%f(I,NT) = cmplx(real(I1,kind(0.d0)), 0.d0,kind(0.d0))
                  elseif (this%t(I) == 3)   then
                     Read(10,*)  X
                     this%f(I,NT) = cmplx(X,0.d0,kind(0.d0))
                  elseif (this%t(I) == 4)   then
                     Read(10,*)  this%f(I,NT)
                  Endif
               ENDDO
            ENDDO
            CLOSE(10)
            DEALLOCATE(SEED_VEC)
        END SUBROUTINE Fields_read_conf
!--------------------------------------------------------------------
#if defined HDF5
        SUBROUTINE Fields_read_conf_h5(this, filename)
            IMPLICIT NONE

            Class (Fields)    , INTENT(INOUT) :: this
            CHARACTER (LEN=64), intent(in)    :: filename

            INTEGER             :: K, hdferr
            INTEGER,ALLOCATABLE :: SEED_VEC(:)
            INTEGER(HID_T)      :: file_id, dset_id, dataspace
            Character (len=64)  :: dset_name
            INTEGER(HSIZE_T), allocatable :: dims(:)

            INTEGER :: rank
            real    (Kind=Kind(0.d0)), allocatable :: f_tmp_real(:,:)
            Complex (Kind=Kind(0.d0)), allocatable, target :: f_tmp_cplx(:,:)
            TYPE(C_PTR)                   :: dat_ptr

            CALL GET_SEED_LEN(K)
            ALLOCATE(SEED_VEC(K))

            CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr)

            !Open and read random seed dataset
            dset_name = "seed"
            allocate( dims(1) )
            dims(1) = K
            CALL h5ltread_dataset_int_f(file_id, dset_name, SEED_VEC, dims, hdferr)
            CALL RANSET(SEED_VEC)
            deallocate( dims, SEED_VEC )

            if (SIZE(this%f) > 0) then
               !Open and read configuration dataset
               dset_name = "configuration"
               !Open the  dataset.
               CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
         
               !Get dataset's dataspace handle.
               CALL h5dget_space_f(dset_id, dataspace, hdferr)
         
               !Get dataspace's rank.
               CALL h5sget_simple_extent_ndims_f(dataspace, rank, hdferr)
               if (rank == 2) then
                  ! rank=2 -> Read real config values
                  allocate( dims(2) )
                  dims = [SIZE(this%f,1), SIZE(this%f,2)]
                  allocate(f_tmp_real(dims(1), dims(2)))
                  CALL h5ltread_dataset_double_f(file_id, dset_name, f_tmp_real, dims, hdferr)
                  this%f(:,:) = f_tmp_real(:,:)
                  deallocate(dims, f_tmp_real)
               else
                  ! Read complex config values
                  allocate( dims(3) )
                  dims = [2, SIZE(this%f,1), SIZE(this%f,2)]
                  allocate(f_tmp_cplx(SIZE(this%f,1), SIZE(this%f,2)))
                  dat_ptr = C_LOC(f_tmp_cplx(1,1))
                  CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
                  this%f(:,:) = f_tmp_cplx(:,:)
                  deallocate(dims, f_tmp_cplx)
               endif

               CALL h5sclose_f(dataspace, hdferr)
               CALL h5dclose_f(dset_id, hdferr)
            endif
            CALL h5fclose_f(file_id, hdferr)
        END SUBROUTINE Fields_read_conf_h5
#endif


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Writes out field configuration for single process, private Subroutine.
!>
!
!> @param [INOUT] this
!> \verbatim
!> Type Fields
!>
!> @param [IN] filename
!> \verbatim
!> Type CHARACTER (LEN=64)
!> Name of file in which to write configuration and random seed \endverbatim
!--------------------------------------------------------------------
        SUBROUTINE Fields_write_conf(this, filename)
#if !defined HDF5
            IMPLICIT NONE

            Class (Fields)    , INTENT(INOUT) :: this
            CHARACTER (LEN=64), intent(in)    :: filename

            INTEGER             :: K, I, NT
            INTEGER,ALLOCATABLE :: SEED_VEC(:)

            CALL GET_SEED_LEN(K)
            ALLOCATE(SEED_VEC(K))
            CALL RANGET(SEED_VEC)

            OPEN (UNIT = 10, FILE=filename, STATUS='UNKNOWN', ACTION='WRITE')
            WRITE(10,*) SEED_VEC
            DO NT = 1,size(this%f,2)
               DO I = 1,size(this%f,1)
                  if (this%t(i) ==  3 ) then
                     WRITE(10,*) real(this%f(I,NT))
                  elseif ( this%t(i) ==  1  .or. this%t(i) ==  2 ) then
                     WRITE(10,*) nint(real(this%f(I,NT)))
                  elseif ( this%t(i) ==  4  ) then
                     WRITE(10,*) this%f(I,NT) 
                  endif
               ENDDO
            ENDDO
            CLOSE(10)
            DEALLOCATE(SEED_VEC)
#else
            IMPLICIT NONE

            Class (Fields)    , INTENT(INOUT) :: this
            CHARACTER (LEN=64), intent(in)    :: filename

            INTEGER             :: K, hdferr, rank
            INTEGER(HSIZE_T), allocatable :: dims(:)
            Logical             :: file_exists
            INTEGER,ALLOCATABLE :: SEED_VEC(:)
            INTEGER(HID_T)      :: file_id, crp_list, space_id, dset_id
            Character (len=64)  :: dset_name
            Complex (Kind=Kind(0.d0)), allocatable, target :: f_tmp(:,:)
            TYPE(C_PTR) :: dat_ptr


            CALL GET_SEED_LEN(K)
            ALLOCATE(SEED_VEC(K))
            CALL RANGET(SEED_VEC)

            inquire (file=filename, exist=file_exists)
            IF (.not. file_exists) THEN
                CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferr)

                !Create and write dataset for random seed
                dset_name = "seed"
                rank = 1
                allocate( dims(1) )
                dims(1) = K
                CALL  h5ltmake_dataset_int_f(file_id, dset_name, rank, dims, SEED_VEC, hdferr)
                deallocate( dims )

                
                if (SIZE(this%f) > 0) then
                  !Create and write dataset for configuration
                  dset_name = "configuration"
                  rank = 3
                  allocate( dims(rank) )
                  dims = [2, SIZE(this%f,1), SIZE(this%f,2)]
                  !Create Dataspace
                  CALL h5screate_simple_f(rank, dims, space_id, hdferr)
                  !Modify dataset creation properties, i.e. enable chunking
                  CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
                  CALL h5pset_chunk_f(crp_list, rank, dims, hdferr)
#if defined HDF5_ZLIB
                  ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
                  CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
                  !Create a dataset using cparms creation properties.
                  CALL h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                                 dset_id, hdferr, crp_list )
                  !Write configuration
                  allocate(f_tmp(SIZE(this%f,1), SIZE(this%f,2)))
                  f_tmp(:,:) = this%f(:,:)
                  dat_ptr = C_LOC(f_tmp(1,1))
                  CALL H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
                  !Close objects
                  deallocate(dims, f_tmp)
                  CALL h5sclose_f(space_id, hdferr)
                  CALL h5pclose_f(crp_list, hdferr)
                  CALL h5dclose_f(dset_id, hdferr)
                endif
                CALL h5fclose_f(file_id, hdferr)
            else
                !open file
                CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

                !open and write random seed dataset
                dset_name = "seed"
                CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
                allocate( dims(1) )
                dims(1) = K
                CALL H5dwrite_f(dset_id, H5T_NATIVE_INTEGER, SEED_VEC, dims, hdferr)
                deallocate( dims )
                CALL h5dclose_f(dset_id,   hdferr)


                if (SIZE(this%f) > 0) then
                  !open and write configuration dataset
                  dset_name = "configuration"
                  CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
                  allocate(f_tmp(SIZE(this%f,1), SIZE(this%f,2)))
                  f_tmp(:,:) = this%f(:,:)
                  dat_ptr = C_LOC(f_tmp(1,1))
                  CALL H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
                  !Close objects
                  deallocate(f_tmp)
                  CALL h5dclose_f(dset_id,   hdferr)
                endif

                CALL h5fclose_f(file_id, hdferr)
            endif
            DEALLOCATE(SEED_VEC)
#endif
            END SUBROUTINE Fields_write_conf
     end Module Fields_Mod
