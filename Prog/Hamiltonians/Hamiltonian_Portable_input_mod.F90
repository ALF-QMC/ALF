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
   use iso_fortran_env, only: output_unit, error_unit

   implicit none

   type operator_matrix
      ! a) hopping <hoppings.txt>:
      ! \hat{H}_t = - \sum_R \sum_s \sum_{\sigma} \sum_n ( \hat{c}^{\dagger}_{x,s,\sigma} T_{n,s} \hat{c}_{y,s,\sigma} + h.c. )
      ! with x = (R, \mu_1^{n,s}) and y = (R + n_1^{n,s} a_1 + n_2^{n,s} a_2, \mu_2^{n,s})
      !
      ! b) interaction <potentials.txt>:
      ! \hat{H}_I = \sum_R \sum_k U_k { \sum_s \sum_{\sigma} [ \sum_n ( \hat{c}^{\dagger}_{x,s,\sigma} g_{k,n,s} \hat{c}_{y,s,\sigma} + h.c. ) + \alpha_{k,s} ] }^2
      ! with x = (R + n_{1,1}^{k,n,s} a_1 + n_{1,2}^{k,n,s} a_2, \mu_1^{k,n,s}) and y = (R + n_{2,1}^{k,n,s} a_1 + n_{2,2}^{k,n,s} a_2, \mu_2^{k,n,s})
      !
      ! c) scalar observable in particle-hole channel <obs_scal_ph.txt>:
      ! \sum_R \langle \hat{O}_R \rangle
      ! with \hat{O}_R =  \sum_n \sum_sigma ( \hat{c}^{\dagger}_{x,s_n,\sigma} O_{n} \hat{c}_{y,s_n,\sigma} )
      ! with x = (R + n_{1,1}^{n} a_1 + n_{1,2}^{n} a_2, \mu_1^{n}) and y = (R + n_{2,1}^{n} a_1 + n_{2,2}^{n} a_2, \mu_2^{n})
      !
      ! d) correlation function in particle-hole channel <obs_corr_ph.txt>:
      ! \sum_{R_1,R_2} ( \langle \hat{O}_{R_1} \hat{O}_{R_2} \rangle - \langle \hat{O}_{R_1} \rangle \langle \hat{O}_{R_2} \rangle )
      ! with \hat{O}_R =  \sum_n ( \hat{c}^{\dagger}_{x,s_{n,1},\sigma_{n,1}} O_{n} \hat{c}_{y,s_{n,2},\sigma_{n,2}} )
      ! with x and y as defined in c)
      !
      ! list      : list that contains information for shift of unit cells, orbitals, and spin
      !             i)   hoppings.txt:
      !                  list(n,1)  = orbital 1 \mu_1^{n,s}
      !                  list(n,2)  = shift of unit cell with vector a_1 n_1^{n,s}
      !                  list(n,3)  = shift of unit cell with vector a_2 n_2^{n,s}
      !                  list(n,4)  = orbital 2 \mu_2^{n,s}
      !             ii)  potentials.txt:
      !                  list(n,1)  = shift 1 of unit cell with vector a_1 n_{1,1}^{k,n,s}
      !                  list(n,2)  = shift 1 of unit cell with vector a_2 n_{1,2}^{k,n,s}
      !                  list(n,3)  = orbital 1 \mu_1^{k,n,s}
      !                  list(n,4)  = shift 2 of unit cell with vector a_1 n_{2,1}^{k,n,s}
      !                  list(n,5)  = shift 2 of unit cell with vector a_2 n_{2,2}^{k,n,s}
      !                  list(n,6)  = orbital 2 \mu_2^{k,n,s}
      !             iii) obs_scal_ph.txt:
      !                  list(n,1)  = shift 1 of unit cell with vector a_1 n_{1,1}^{n}
      !                  list(n,2)  = shift 1 of unit cell with vector a_2 n_{1,2}^{n}
      !                  list(n,3)  = orbital 1 \mu_1^{n}
      !                  list(n,4)  = flavor 1 s_n
      !                  list(n,5)  = shift 2 of unit cell with vector a_1 n_{2,1}^{n}
      !                  list(n,6)  = shift 2 of unit cell with vector a_2 n_{2,2}^{n}
      !                  list(n,7)  = orbital 2 \mu_2^{n}
      !             iv)  obs_corr_ph.txt:
      !                  list(n,1)  = shift 1 of unit cell with vector a_1 n_{1,1}^{n}
      !                  list(n,2)  = shift 1 of unit cell with vector a_2 n_{1,2}^{n}
      !                  list(n,3)  = orbital 1 \mu_1^{n}
      !                  list(n,4)  = flavor 1 s_{n,1}
      !                  list(n,5)  = color 1 \sigma_{n,1}
      !                  list(n,6)  = shift 2 of unit cell with vector a_1 n_{2,1}^{n}
      !                  list(n,7)  = shift 2 of unit cell with vector a_2 n_{2,2}^{n}
      !                  list(n,8)  = orbital 2 \mu_2^{n}
      !                  list(n,9)  = flavor 2 s_{n,2}
      !                  list(n,10) = color 2 \sigma_{n,2}
      ! g         : matrix elements for hopping (T_{n,s}), interaction (g_{k,n,s}), observables (O_{n})
      ! u         : only for interactions b)
      !             prefactor U_k of interaction terms
      ! alpha     : only for interactions b)
      !             complex shift \alpha_{k,s} in interaction term
      integer, allocatable                   :: list(:,:)
      complex (kind=kind(0.d0)), allocatable :: g(:)
      real (kind=kind(0.d0))    :: u
      complex (kind=kind(0.d0)) :: alpha
   end type operator_matrix

   type obs_orbitals
      integer, allocatable :: orb_list(:,:,:)
      integer, allocatable :: nonzero_corr(:,:,:), nonzero_back(:)
   end type obs_orbitals

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
          CALL MPI_BCAST(Orb_pos  ,size(Orb_pos) ,MPI_REAL8    ,0,Group_Comm,ierr)
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

          type (operator_matrix), intent(out) :: this(:)
          integer, intent(in)      :: Group_Comm

          integer                :: ierr, unit_hop, n_hop, nh, i
          integer                :: no1, no2, n1, n2, nf
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

             do nf = 1, size(this,1)

                read(unit_hop,*) n_hop
             
                allocate( this(nf)%g(n_hop), this(nf)%list(n_hop,4) )
      
                do nh = 1, n_hop
                   read(unit_hop,*) i, no1, n1, n2, no2, x, y
                   this(nf)%list(nh,1) = no1
                   this(nf)%list(nh,2) = n1
                   this(nf)%list(nh,3) = n2
                   this(nf)%list(nh,4) = no2
                   this(nf)%g(nh)      = cmplx( x, y, kind(0.d0))
                enddo

             enddo

             Close(unit_hop)
#ifdef MPI
          Endif

          do nf = 1, size(this,1)
             if ( irank_g == 0 ) n_hop = size(this(nf)%list,1)
             CALL MPI_BCAST(n_hop        ,  1                ,MPI_INTEGER  ,0,Group_Comm,ierr)

             if ( irank_g /= 0 ) allocate( this(nf)%g(n_hop), this(nf)%list(n_hop,4) )
             CALL MPI_BCAST(this(nf)%list,size(this(nf)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(this(nf)%g   ,size(this(nf)%g)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
          enddo
 
#endif

        end subroutine read_hop

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in   the  interaction parameters
!--------------------------------------------------------------------

        Subroutine read_int(this,N_fl,Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          type (operator_matrix), allocatable, intent(out) :: this(:,:)
          integer, intent(in)    :: Group_Comm, N_FL
 
          integer                :: ierr, unit_int, mk, no, n, i1, i2
          integer                :: no1, j1, j2, no2, i, N_ops, nf
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

             allocate( this(N_ops,N_Fl) )

             do no = 1, N_ops
                read(unit_int,*) u
                this(no,1)%u   = u
                do nf = 1, N_Fl
                   read(unit_int,*) mk, x, y
                   this(no,nf)%alpha = cmplx( x, y, kind(0.d0) )

                   allocate( this(no,nf)%list(mk,6), this(no,nf)%g(mk) )

                   do n = 1, mk
                      read(unit_int,*) i, i1, i2, no1, j1, j2, no2, x, y
                      this(no,nf)%list(n,1) = i1
                      this(no,nf)%list(n,2) = i2
                      this(no,nf)%list(n,3) = no1
                      this(no,nf)%list(n,4) = j1
                      this(no,nf)%list(n,5) = j2
                      this(no,nf)%list(n,6) = no2
                      this(no,nf)%g(n)      = cmplx( x, y, kind(0.d0) )
                   enddo
                enddo
             enddo

             Close(unit_int)

#ifdef MPI
          Endif

          CALL MPI_BCAST(n_ops          ,  1  ,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate( this(N_ops,N_Fl) )
      
          do no = 1, N_ops
             CALL MPI_BCAST(this(no,1)%u,  1  ,MPI_REAL8    ,0,Group_Comm,ierr)
             do nf = 1, N_Fl
                CALL MPI_BCAST(this(no,nf)%alpha,  1                   ,MPI_COMPLEX16,0,Group_Comm,ierr)
             
                if (irank_g == 0) mk = size(this(no,nf)%list,1)
                CALL MPI_BCAST(mk               ,  1                   ,MPI_INTEGER  ,0,Group_Comm,ierr)
                if (irank_g /= 0) allocate( this(no,nf)%list(mk,6), this(no,nf)%g(mk) )
                CALL MPI_BCAST(this(no,nf)%list ,size(this(no,nf)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
                CALL MPI_BCAST(this(no,nf)%g    ,size(this(no,nf)%g)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
             enddo
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

        Subroutine read_obs_scal_ph(this,filenames,Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          type (operator_matrix), allocatable, intent(out) :: this(:)
          integer, intent(in) :: Group_Comm
          character (len=64), allocatable, intent(out) :: filenames(:)
 
          integer                :: ierr, unit_obs, mk, no, n, N_obs
          integer                :: i, no1, no2, nf, i1, i2, j1, j2
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
             allocate(this(N_obs),filenames(N_obs))
      
             do no = 1, N_obs
                read(unit_obs,*) file, mk
                filenames(no) = file
                allocate( this(no)%list(mk,7), this(no)%g(mk) )
                
                do n = 1, mk
                   read(unit_obs,*) i, i1, i2, no1, nf, j1, j2, no2, x, y
                   this(no)%list(n,1) = i1
                   this(no)%list(n,2) = i2
                   this(no)%list(n,3) = no1

                   this(no)%list(n,4) = nf

                   this(no)%list(n,5) = j1
                   this(no)%list(n,6) = j2
                   this(no)%list(n,7) = no2

                   this(no)%g(n) = cmplx(x, y, kind(0.d0))
                enddo
             enddo

             Close(unit_obs)

#ifdef MPI
          Endif

          CALL MPI_BCAST(N_obs                 ,  1,MPI_INTEGER   ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate(this(N_obs),filenames(N_obs))
          call MPI_BCAST(filenames             ,64*size(filenames), MPI_CHARACTER,0,Group_Comm,ierr)
          do no = 1, N_obs
             if (irank_g == 0) mk = size(this(no)%list,1)
             CALL MPI_BCAST(mk              ,  1                   ,MPI_INTEGER  ,0,Group_Comm,ierr)
             if (irank_g /= 0) allocate( this(no)%list(mk,7), this(no)%g(mk) )
             CALL MPI_BCAST(this(no)%list,size(this(no)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(this(no)%g   ,size(this(no)%g)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
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

        Subroutine read_obs_corr(this,filenames,Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          type (operator_matrix), allocatable, intent(out) :: this(:)
          integer, intent(in) :: Group_Comm
          character (len=64), allocatable, intent(out) :: filenames(:)
 
          integer                :: ierr, unit_obs, mk, no, n, N_obs
          integer                :: i, no1, i1, i2, j1, j2, no2, ns1, ns2
          integer                :: nf1, nf2
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

             read(unit_obs,*) N_obs
             allocate(this(N_obs),filenames(N_obs))

             do no = 1, N_obs
                read(unit_obs,*) file, mk
                filenames(no)  = file
                allocate( this(no)%list(mk,10), this(no)%g(mk) )

                do n = 1, mk
                   read(unit_obs,*) i, i1, i2, no1, nf1, ns1, j1, j2, no2, nf2, ns2, x, y
                   this(no)%list(n,1)  = i1
                   this(no)%list(n,2)  = i2
                   this(no)%list(n,3)  = no1
                   this(no)%list(n,4)  = nf1
                   this(no)%list(n,5)  = ns1

                   this(no)%list(n,6)  = j1
                   this(no)%list(n,7)  = j2
                   this(no)%list(n,8)  = no2
                   this(no)%list(n,9)  = nf2
                   this(no)%list(n,10) = ns2

                   this(no)%g(n) = cmplx(x, y, kind(0.d0))
                enddo
             enddo

             Close(unit_obs)

#ifdef MPI
          Endif

          CALL MPI_BCAST(N_obs                 ,   1,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate(this(N_obs), filenames(N_obs))
          CALL MPI_BCAST(filenames             ,64*size(filenames),MPI_CHARACTER,0,Group_Comm,ierr)
          do no = 1, N_obs
             if (irank_g == 0) mk = size(this(no)%list,1)
             CALL MPI_BCAST(mk              ,  1                     ,MPI_INTEGER  ,0,Group_Comm,ierr)
             if (irank_g /= 0) allocate( this(no)%list(mk,10), this(no)%g(mk) )
             CALL MPI_BCAST(this(no)%list,  size(this(no)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(this(no)%g   ,  size(this(no)%g)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
          enddo

#endif

        end subroutine read_obs_corr

end module Hamiltonian_Portable_input_mod
