!  Copyright (C) 2016 - 2018 The ALF project
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
!--------------------------------------------------------------------

     Module Fields_mod

       Use Random_Wrap
       use iso_fortran_env, only: output_unit, error_unit

       Public Fields
       Public Fields_init

       Private
       Real (Kind=Kind(0.d0))  :: Phi_st(-2:2,2),  Gama_st(-2:2,2)
       Real (Kind=Kind(0.d0))  :: Del, FLIP_st(-2:2,3)
       Real (Kind=Kind(0.d0))  :: Amplitude=5.d0

       Type Fields
          Real    (Kind=Kind(0.d0)), allocatable    :: f(:,:)
          Integer                  , allocatable    :: t(:)
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

      Real (Kind=Kind(0.d0)) function  Fields_Phi(this,n_op,n_tau)

        Implicit none
        Class (Fields) :: this
        Integer, INTENT(IN) ::  n_op, n_tau


        select case (this%t(n_op))
        case(1)
           Fields_Phi = Phi_st(Nint(this%f(n_op,n_tau)),1)
        case(2)
           Fields_Phi = Phi_st(Nint(this%f(n_op,n_tau)),2)
        case(3)
           Fields_Phi = this%f(n_op,n_tau)
        case default
           Write(error_unit,*) 'Error in Fields_Phi'
           error stop 1
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
           Fields_GAMA = GAMA_st(Nint(this%f(n_op,n_tau)),2)
        case(3)
           Fields_GAMA = 1.d0
        case default
           Write(error_unit,*) 'Error in Fields_GAMA'
           error stop 1
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

      Real (Kind=Kind(0.d0)) function Fields_flip(this,n_op,n_tau)

        Implicit none
        Class (Fields) :: this
        Integer, INTENT(IN) ::  n_op, n_tau

        select case (this%t(n_op))
        case(1)
           Fields_flip = - this%f(n_op,n_tau)
        case (2)
           Fields_flip =   Flip_st( nint(this%f(n_op,n_tau)),nranf(3))
        case (3)
           Fields_flip =   this%f(n_op,n_tau) + Amplitude*( ranf_wrap() - 0.5D0)
        case default
           Write(error_unit,*) 'Error in Fields. '
           error stop 1
        end select

      end function Fields_Flip


      Integer function Fields_get_i(this,n_op,n_tau)

        Implicit none
        Class (Fields) :: this
        Integer, INTENT(IN) ::  n_op, n_tau

        if ( this%t(n_op) == 1 .or.   this%t(n_op) == 2 ) then
           Fields_get_i = NINT(this%f(n_op,n_tau))
        else
           Write(error_unit,*) "Error in fields"
           error stop 1
        endif

      end function Fields_get_i


      Subroutine Fields_make(this,N_OP,N_tau)
        Implicit none
        Class (Fields), INTENT(INOUT)  :: this
        Integer, INTENT(IN)            :: N_OP, N_tau

        !Write(6,*) "Allocating  fields: ", N_op, N_tau
        allocate (this%f(N_OP,N_tau), this%t(N_OP) )

        this%f = 0.d0;  this%t = 0

      end Subroutine Fields_make

      Subroutine Fields_clear(this)
        Implicit none
        Class (Fields) :: this

        deallocate (this%f, this%t )
      end Subroutine Fields_clear

      Subroutine Fields_init(Delta_X)

        Implicit none

        Real  (Kind=Kind(0.d0)), Optional, Intent(IN) :: Delta_X

        !Local
        Integer :: n

        Del = 0.d0
        If (Present(Delta_X)) Del = Delta_X

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
!> routine will generate one based on the seeds read in from the file
!> seeds
!
!> @param [INOUT] this
!> \verbatim
!> Type Fields
!> On input  test%t(:) is set.  The operator types are time independent.
!> On output test%f(:,:) is initialized \endverbatim
!>
!> @param [IN] Group_Comm
!> \verbatim
!> Type Integer
!> Communicator for MPI \endverbatim
!>
!> @param [Optional]  Initial_field
!> \verbatim
!> Type Real
!> Initial field \endverbatim
!--------------------------------------------------------------------
      Subroutine Fields_in(this,Group_Comm,Initial_field)

#ifdef MPI
        Use mpi
#endif

        Implicit none

        Class (Fields)        , INTENT(INOUT) :: this
        Integer               , INTENT(IN   ) :: Group_Comm
        Real (Kind=Kind(0.d0)), Dimension(:,:), Optional   :: Initial_field

        ! LOCAL
        Integer                 :: I, I1, IERR, SEED_IN, K, NT
        Real (Kind=Kind(0.d0) ) :: X
        Integer, DIMENSION(:), ALLOCATABLE :: SEED_VEC
        Logical ::   LCONF
        Character (LEN=64) :: FILE_SR, FILE_TG, FILE_seeds, FILE_info, File1

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
            write(FILE1,'(A,I0,A)') "Temp_",igroup,"/confin_0"
#else
            File1 = "confin_0"
#endif
            INQUIRE (FILE=File1, EXIST=LCONF)
            IF (LCONF) THEN
#if defined(TEMPERING)
               write(FILE_TG,'(A,I0,A,I0)') "Temp_",igroup,"/confin_",irank_g
#else
               write(FILE_TG,'(A,I0)') "confin_",irank_g
#endif
               CALL GET_SEED_LEN(K)
               ALLOCATE(SEED_VEC(K))
               OPEN (UNIT = 10, FILE=FILE_TG, STATUS='OLD', ACTION='READ')
               READ(10,*) SEED_VEC
               CALL RANSET(SEED_VEC)
               DO NT = 1,SIZE(this%f,2)
                  DO I = 1,SIZE(this%f,1)
                     IF (this%t(I) == 1 .or.  this%t(I) == 2) then
                        Read(10,*)  I1
                        this%f(I,NT) = real(I1,kind(0.d0))
                     else
                        Read(10,*)  this%f(I,NT)
                     Endif
                  Enddo
               Enddo
               CLOSE(10)
               DEALLOCATE(SEED_VEC)
            ELSE
               IF (IRANK == 0) THEN
                  WRITE(6,*) 'No initial configuration'
                  OPEN(UNIT=5,FILE='seeds',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
                  IF (IERR /= 0) THEN
                     WRITE(error_unit,*) 'Fields_in: unable to open <seeds>',IERR
                     error stop 1
                  END IF
                  DO I = ISIZE-1,1,-1
                     READ (5,*) SEED_IN
                     CALL MPI_SEND(SEED_IN,1,MPI_INTEGER, I, I+1024, MPI_COMM_WORLD,IERR)
                  ENDDO
                  READ(5,*) SEED_IN
                  CLOSE(5)
               ELSE
                  CALL MPI_RECV(SEED_IN, 1, MPI_INTEGER,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
               ENDIF
               ALLOCATE (SEED_VEC(1))
               SEED_VEC(1) = SEED_IN
               CALL RANSET(SEED_VEC)
               DEALLOCATE (SEED_VEC)
               If (Present(Initial_field)) then
                  this%f = Initial_field
               else
                  Call  this%set()
               endif
               if (irank_g == 0) then
#if defined(TEMPERING)
                  write(FILE_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
                  FILE_info="info"
#endif
                  Open (Unit = 50,file=FILE_info,status="unknown",position="append")
                  WRITE(50,*) 'No initial configuration, Seed_in', SEED_IN
                  Close(50)
               endif
         ENDIF

#else
         FILE_TG = "confin_0"
         INQUIRE (FILE=FILE_TG, EXIST=LCONF)
         IF (LCONF) THEN
            CALL GET_SEED_LEN(K)
            ALLOCATE(SEED_VEC(K))
            OPEN (UNIT = 10, FILE=FILE_TG, STATUS='OLD', ACTION='READ')
            READ(10,*) SEED_VEC
            CALL RANSET(SEED_VEC)
            DO NT = 1,SIZE(this%f,2)
               DO I = 1,SIZE(this%f,1)
                  IF (this%t(I) == 1 .or.  this%t(I) == 2) then
                     Read(10,*)  I1
                     this%f(I,NT) = real(I1,kind(0.d0))
                  else
                     Read(10,*)  this%f(I,NT)
                  Endif
               Enddo
            Enddo
            DEALLOCATE(SEED_VEC)
         ELSE
            FILE_seeds="seeds"
            OPEN(UNIT=5,FILE=FILE_seeds,STATUS='OLD',ACTION='READ',IOSTAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'Fields_in: unable to open <seeds>',IERR
               error stop 1
            END IF
            READ (5,*) SEED_IN
            CLOSE(5)
            FILE_info="info"
            Open (Unit = 50,file=FILE_info,status="unknown",position="append")
            WRITE(50,*) 'No initial configuration, Seed_in', SEED_IN
            Close(50)

            ALLOCATE(SEED_VEC(1))
            SEED_VEC(1) = SEED_IN
            CALL RANSET (SEED_VEC)
            DEALLOCATE  (SEED_VEC)
            If (Present(Initial_field)) then
               this%f = Initial_field
            else
               Call  this%set()
            endif
         ENDIF
#endif
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

#ifdef MPI
         Use mpi
#endif
         IMPLICIT NONE

         Class (Fields), INTENT(INOUT) :: this
         Integer,        INTENT(IN   ) :: Group_Comm

         ! LOCAL
         INTEGER        :: I, K, NT
         INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_VEC
         CHARACTER (LEN=64) :: FILE_TG

#if defined(MPI)
         INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
         !Write(6,*) "Group, rank :", igroup, irank_g

         CALL GET_SEED_LEN(K)
         ALLOCATE(SEED_VEC(K))
         CALL RANGET(SEED_VEC)
#if defined(TEMPERING)
         write(FILE_TG,'(A,I0,A,I0)') "Temp_",igroup,"/confout_",irank_g
#else
         write(FILE_TG,'(A,I0)') "confout_",irank_g
#endif
         OPEN (UNIT = 10, FILE=FILE_TG, STATUS='UNKNOWN', ACTION='WRITE')
         WRITE(10,*) SEED_VEC
         DO NT = 1,size(this%f,2)
            DO I = 1,size(this%f,1)
               if (this%t(i) ==  3 ) then
                  WRITE(10,*) this%f(I,NT)
               else
                  WRITE(10,*) nint(this%f(I,NT))
               endif
            ENDDO
         ENDDO
         CLOSE(10)
         DEALLOCATE(SEED_VEC)
#else
         CALL GET_SEED_LEN(K)
         ALLOCATE(SEED_VEC(K))
         CALL RANGET(SEED_VEC)
         FILE_TG = "confout_0"
         OPEN (UNIT = 10, FILE=FILE_TG, STATUS='UNKNOWN', ACTION='WRITE')
         WRITE(10,*) SEED_VEC
         DO NT = 1,size(this%f,2)
            DO I = 1,size(this%f,1)
               if (this%t(i) ==  3 ) then
                  WRITE(10,*) this%f(I,NT)
               else
                  WRITE(10,*) nint(this%f(I,NT))
               endif
            ENDDO
         ENDDO
         CLOSE(10)
         DEALLOCATE(SEED_VEC)
#endif

       END SUBROUTINE Fields_out


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

         Integer :: nt, I

         !Write(6,*) "Fields_set", size(this%f,1), size(this%f,2)
         Do nt = 1,size(this%f,2)
            Do I = 1,size(this%f,1)
               if (this%t(i)  < 4 ) then
                  this%f(I,nt)  = 1.d0
                  if ( ranf_wrap() > 0.5D0 ) this%f(I,nt) = -1.d0
               else
                  this%f(I,nt)  = ranf_wrap() - 0.5d0
               endif
            enddo
         enddo

       end Subroutine Fields_set
!---------------------------------------------------------------------


     end Module Fields_Mod
