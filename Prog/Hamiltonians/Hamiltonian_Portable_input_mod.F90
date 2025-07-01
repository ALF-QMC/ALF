!  Copyright (C) 2016 - 2022 The ALF project
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
!> Reads in the parameters for the construction of the Hamiltonian from the files geometry.txt (Ham_latt), hoppings.txt (Ham_hop) and potentials.txt (Ham_V) and for the observables.
!>
!--------------------------------------------------------------------

module Hamiltonian_Portable_input_mod

   Use runtime_error_mod

   implicit none

   type operator_matrix
      ! file      : only for observables , file name in the form of <file>_scal, <file>_eq or <file>_tau
      ! list      :                        list that contains information for shift of unit cells, orbitals, and spin
      !                                    i) hoppings.txt, obs_scal_ph.txt:
      !                                       list(:,1) = orbital 1
      !                                       list(:,2) = spin 1
      !                                       list(:,3) = shift of unit cell with vector a_1
      !                                       list(:,4) = shift of unit cell with vector a_2
      !                                       list(:,5) = orbital 2
      !                                       list(:,6) = spin 2
      !                                    i) potentials.txt, obs_corr_ph.txt:
      !                                       list(:,1) = shift 1 of unit cell with vector a_1
      !                                       list(:,2) = shift 1 of unit cell with vector a_1
      !                                       list(:,3) = orbital 1
      !                                       list(:,4) = spin 1
      !                                       list(:,5) = shift 2 of unit cell with vector a_1
      !                                       list(:,6) = shift 2 of unit cell with vector a_1
      !                                       list(:,7) = orbital 2
      !                                       list(:,8) = spin 2
      ! g         :                        matrix elements for hopping, interaction, observables
      ! N_orbitals: only for interactions, number of interacting orbitals per interaction term, determines the effective dimension
      !                                    of the operator in Op_make(OP_V, N_orbitals)
      ! u         : only for interactions, prefactor U_k of interaction terms
      ! alpha     : only for interactions, complex shift \alpha_k in interaction term
      Character (len=64)        :: File
      integer, allocatable                   :: list(:,:)
      complex (kind=kind(0.d0)), allocatable :: g(:)
      integer                   :: N_orbitals
      real (kind=kind(0.d0))    :: u
      complex (kind=kind(0.d0)) :: alpha
   end type operator_matrix

   contains

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in  the  Lattice parameters
!--------------------------------------------------------------------

        Subroutine read_latt(L1, l2, Norb, a_p, Orb_pos, Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          Real (Kind=Kind(0.d0)), intent(out)  :: a_p(3,3)
          real (kind=kind(0.d0)), allocatable, intent(out) :: Orb_pos(:,:)
          integer, intent(out) :: L1, L2, Norb
          integer, intent(in)  :: Group_Comm
 
          integer                :: ierr, unit_latt, no, i
          Character (len=64)     :: file_latt
          real (kind=kind(0.d0)) :: x, y, z

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

          If (Irank_g == 0) then
#endif
             File_latt = "geometry.txt"
#if defined(TEMPERING)
             write(File_latt,'(A,I0,A)') "Temp_",igroup,"/geometry.txt"
#endif
             Open(newunit=unit_latt, file=file_latt, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <geometry.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF

             read(unit_latt,*) L1, L2
             do i = 1, size(a_p,1)
                read(unit_latt,*) x, y, z
                a_p(i,1) = x; a_p(i,2) = y; a_p(i,3) = z
             enddo
             read(unit_latt,*) norb
             allocate(Orb_pos(Norb,3))
             do no = 1, Norb
                read(unit_latt,*) i, x, y, z
                Orb_pos(no,1) = x; Orb_pos(no,2) = y; Orb_pos(no,3) = z
             enddo

             Close(unit_latt)
#ifdef MPI
          Endif

          CALL MPI_BCAST(L1       ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(L2       ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(a_p      ,  size(a_p)   ,MPI_REAL8    ,0,Group_Comm,ierr)
          CALL MPI_BCAST(Norb     ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate(Orb_pos(Norb,3))
          CALL MPI_BCAST(Orb_pos,  size(Orb_pos) ,MPI_REAL8    ,0,Group_Comm,ierr)
#endif

          If (L1 < 1 .or. L2 < 1) then
             Write(error_unit,*) 'The lattice sizes L1 and L2 have to be positive integers.'
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
          If (Norb < 1) then
             Write(error_unit,*) 'The number of orbitals per unit cell Norb has to be a positive integer.'
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

        end subroutine read_latt

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads  the  hopping parameters
!--------------------------------------------------------------------

        Subroutine read_hop(this, Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          type (operator_matrix), intent(out) :: this
          integer, intent(in)      :: Group_Comm

          integer                :: ierr, unit_hop, n_hop, nh, i
          integer                :: no1, no2, s1, s2, n1, n2
          Character (len=64)     :: file_hop
          real (kind=kind(0.d0)) :: x, y

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

          If (Irank_g == 0) then
#endif
             File_hop = "hoppings.txt"
#if defined(TEMPERING)
             write(File_hop,'(A,I0,A)') "Temp_",igroup,"/hoppings.txt"
#endif
             Open(newunit=unit_hop, file=file_hop, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <hoppings.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF

             read(unit_hop,*) n_hop
             
             allocate( this%g(n_hop), this%list(n_hop,6) )
      
             do nh = 1, n_hop
                read(unit_hop,*) i, no1, s1, n1, n2, no2, s2, x, y
                this%list(nh,1) = no1
                this%list(nh,2) = s1
                this%list(nh,3) = n1
                this%list(nh,4) = n2
                this%list(nh,5) = no2
                this%list(nh,6) = s2
                this%g(nh)      = cmplx( x, y, kind(0.d0))
             enddo

             Close(unit_hop)
#ifdef MPI
          Endif

          CALL MPI_BCAST(n_hop       ,  1             ,MPI_INTEGER  ,0,Group_Comm,ierr)

          if ( irank_g /= 0 ) allocate( this%g(n_hop), this%list(n_hop,6) )
          CALL MPI_BCAST(this%list    ,  size(this%list),MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(this%g       ,  size(this%g)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
 
#endif

        end subroutine read_hop

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in   the  interaction parameters
!--------------------------------------------------------------------

        Subroutine read_int(this,Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          type (operator_matrix), allocatable, intent(out) :: this(:)
          integer, intent(in)    :: Group_Comm
 
          integer                :: ierr, unit_int, mk, no, n, i1, i2
          integer                :: no1, s1, j1, j2, no2, s2, i, N_ops
          Character (len=64)     :: file_int
          real (kind=kind(0.d0)) :: u, x, y

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
   
          If (Irank_g == 0) then
#endif
             File_int = "potentials.txt"
#if defined(TEMPERING)
             write(File_int,'(A,I0,A)') "Temp_",igroup,"/potentials.txt"
#endif
             Open(newunit=unit_int, file=file_int, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <potentials.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF

             read(unit_int,*) N_ops

             allocate( this(N_ops) )

             do no = 1, N_ops
                read(unit_int,*) mk, u, x, y
                this(no)%u     = u
                this(no)%alpha = cmplx( x, y, kind(0.d0) )

                allocate( this(no)%list(mk,8), this(no)%g(mk) )

                do n = 1, mk
                   read(unit_int,*) i, i1, i2, no1, s1, j1, j2, no2, s2, x, y
                   this(no)%list(n,1) = i1
                   this(no)%list(n,2) = i2
                   this(no)%list(n,3) = no1
                   this(no)%list(n,4) = s1
                   this(no)%list(n,5) = j1
                   this(no)%list(n,6) = j2
                   this(no)%list(n,7) = no2
                   this(no)%list(n,8) = s2
                   this(no)%g(n)      = cmplx( x, y, kind(0.d0) )
                enddo
             enddo

             Close(unit_int)

#ifdef MPI
          Endif

          CALL MPI_BCAST(n_ops       ,  1              ,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate( this(N_ops) )
      
          do no = 1, N_ops
             CALL MPI_BCAST(this(no)%alpha  ,     1    ,MPI_COMPLEX16,0,Group_Comm,ierr)
             CALL MPI_BCAST(this(no)%u      ,     1    ,MPI_REAL8    ,0,Group_Comm,ierr)
             if (irank_g == 0) mk = size(this(no)%list,1)
             CALL MPI_BCAST(mk       ,  1              ,MPI_INTEGER  ,0,Group_Comm,ierr)
             if (irank_g /= 0) allocate( this(no)%list(mk,8), this(no)%g(mk) )
             CALL MPI_BCAST(this(no)%list   ,size(this(no)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(this(no)%g      ,   size(this(no)%g),MPI_COMPLEX16,0,Group_Comm,ierr)
          enddo

#endif

        end subroutine read_int

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in scalar observables in the particle-hole channel
!--------------------------------------------------------------------

        Subroutine read_obs_scal_ph(this,Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          type (operator_matrix), allocatable, intent(out) :: this(:)
          integer, intent(in) :: Group_Comm
 
          integer                :: ierr, unit_obs, mk, no, n, N_obs
          integer                :: i, no1, s1, n1, n2, no2, s2
          Character (len=64)     :: file_obs, file
          real (kind=kind(0.d0)) :: x, y

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

          If (Irank_g == 0) then
#endif
             File_obs = "obs_scal_ph.txt"
#if defined(TEMPERING)
             write(File_obs,'(A,I0,A)') "Temp_",igroup,"/obs_scal_ph.txt"
#endif
             Open(newunit=unit_obs, file=file_obs, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <obs_scal_ph.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF

             read(unit_obs,*) N_obs
             allocate(this(N_obs))
      
             do no = 1, N_obs
                read(unit_obs,*) file, mk
                this(no)%file = file
                allocate( this(no)%list(mk,6), this(no)%g(mk) )
                
                do n = 1, mk
                   read(unit_obs,*) i, no1, s1, n1, n2, no2, s2, x, y
                   this(no)%list(n,1) = no1
                   this(no)%list(n,2) = s1
                   this(no)%list(n,3) = n1
                   this(no)%list(n,4) = n2
                   this(no)%list(n,5) = no2
                   this(no)%list(n,6) = s2
                   this(no)%g(n)     = cmplx(x, y, kind(0.d0))
                enddo
             enddo

             Close(unit_obs)

#ifdef MPI
          Endif

          CALL MPI_BCAST(N_obs                     ,  1,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate(this(N_obs))
          do no = 1, N_obs
             CALL MPI_BCAST(this(no)%file   , 64,MPI_CHARACTER,0,Group_Comm,ierr)
             if (irank_g == 0) mk = size(this(no)%list,1)
             CALL MPI_BCAST(mk              ,  1                  ,MPI_INTEGER  ,0,Group_Comm,ierr)
             if (irank_g /= 0) allocate( this(no)%list(mk,6), this(no)%g(mk) )
             CALL MPI_BCAST(this(no)%list   ,  size(this(no)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(this(no)%g      ,  size(this(no)%g)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
          enddo

#endif

        end subroutine read_obs_scal_ph

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in correlation functions in the particle-hole channel
!--------------------------------------------------------------------

        Subroutine read_obs_corr(this,Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          type (operator_matrix), allocatable, intent(out) :: this(:)
          integer, intent(in) :: Group_Comm
 
          integer                :: ierr, unit_obs, mk, no, n, N_obs
          integer                :: i, no1, s1, i1, i2, j1, j2, no2, s2
          Character (len=64)     :: file_obs, file
          real (kind=kind(0.d0)) :: x, y

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

   
          If (Irank_g == 0) then
#endif
             File_obs = "obs_corr_ph.txt"
#if defined(TEMPERING)
             write(File_obs,'(A,I0,A)') "Temp_",igroup,"/obs_corr_ph.txt"
#endif
             Open(newunit=unit_obs, file=file_obs, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <obs_corr_ph.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF
         
             !local correlation functions
             read(unit_obs,*) N_obs
             allocate(this(N_obs))

             do no = 1, N_obs
                read(unit_obs,*) file, mk
                this(no)%file  = file
                allocate( this(no)%list(mk,8), this(no)%g(mk) )

                do n = 1, mk
                   read(unit_obs,*) i, i1, i2, no1, s1, j1, j2, no2, s2, x, y
                   this(no)%list(n,1) = i1
                   this(no)%list(n,2) = i2
                   this(no)%list(n,3) = no1
                   this(no)%list(n,4) = s1
                   this(no)%list(n,5) = j1
                   this(no)%list(n,6) = j2
                   this(no)%list(n,7) = no2
                   this(no)%list(n,8) = s2
                   this(no)%g(n)      = cmplx( x, y, kind(0.d0) )
                enddo
             enddo

             Close(unit_obs)

#ifdef MPI
          Endif

          CALL MPI_BCAST(N_obs             ,  1,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate(this(N_obs))
          do no = 1, N_obs
             CALL MPI_BCAST(this(no)%file  , 64,MPI_CHARACTER,0,Group_Comm,ierr)
             if (irank_g == 0) mk = size(this(no)%list,1)
             CALL MPI_BCAST(mk             ,  1                  ,MPI_INTEGER  ,0,Group_Comm,ierr)
             if (irank_g /= 0) allocate( this(no)%list(mk,8), this(no)%g(mk) )
             CALL MPI_BCAST(this(no)%list  ,  size(this(no)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(this(no)%g     ,  size(this(no)%g)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
          enddo

#endif

        end subroutine read_obs_corr

end module Hamiltonian_Portable_input_mod
