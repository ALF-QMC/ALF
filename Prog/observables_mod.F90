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

       Use Lattices_v3, only: Unit_cell, Lattice
       use iso_fortran_env, only: output_unit, error_unit

       Type Obser_Vec
!>  Data structure for
!>  < O_n >  n : =1, size(Obs,1)
          Integer                     :: N                    ! Number of measurements
          real      (Kind=Kind(0.d0)) :: Ave_Sign             ! Averarge sign
          complex   (Kind=Kind(0.d0)), pointer :: Obs_vec(:)  ! Vector of observables
          Character (len=64) :: File_Vec                      ! Name of file in which the bins will be written out
          Character (len=64) :: analysis_mode                 ! How to analyze the observable
          Character (len=64), allocatable :: description(:)   ! Optional short description
       end type Obser_Vec


       Type Obser_Latt
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
       end type Obser_Latt



       Contains

         Subroutine Obser_Latt_make_eq(Obs, Filename, Latt, Latt_unit)
           Implicit none
           Type(Obser_Latt),  Intent(INOUT)      :: Obs
           Character(len=64), Intent(IN)         :: Filename
           Type(Lattice),     Intent(IN), target :: Latt
           Type(Unit_cell),   Intent(IN), target :: Latt_unit
           
           Integer               :: Nt
           Character(len=2)      :: Channel
           Real(Kind=Kind(0.d0)) :: dtau
           
           Nt = 1
           Channel = '--'
           dtau = -1.d0
           
           call Obser_Latt_make(Obs, Nt, Filename, Latt, Latt_unit, Channel, dtau)
         end subroutine Obser_Latt_make_eq

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
           Type(Obser_Latt),  Intent(INOUT)      :: Obs
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
           Type (Obser_Latt), intent(INOUT) :: Obs
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
           Type (Obser_vec), intent(INOUT) :: Obs
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
           Type (Obser_vec), intent(INOUT) :: Obs
           Obs%Obs_vec = cmplx(0.d0,0.d0,kind(0.d0))
           Obs%N       = 0
           Obs%Ave_Sign= 0.d0
         end subroutine Obser_Vec_Init

!--------------------------------------------------------------------

         Subroutine  Print_bin_Latt(Obs, Group_Comm)
           Use Lattices_v3
#ifdef MPI
           Use mpi
#endif
           Implicit none

           Type (Obser_Latt),        Intent(Inout)   :: Obs
           Integer,                  Intent(In)      :: Group_Comm

           ! Local
           Integer :: Ns, Nt, no, no1, I, Ntau
           Complex (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:,:)
           Real    (Kind=Kind(0.d0))              :: x_p(2)
           Complex (Kind=Kind(0.d0))              :: Sign_bin
           Character (len=64) :: File_pr,  File_suff, File_aux, tmp_str
           logical            :: File_exists
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
              error stop 1
           endif
           If (Ntau == 1) then
              File_suff = "_eq"
           else
              File_suff = "_tau"
           endif
           write(File_pr, '(A,A)') trim(Obs%File_Latt), Trim(File_suff)
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
#if defined(TEMPERING)
              write(File_pr ,'(A,I0,A,A,A)') "Temp_",igroup,"/",trim(Obs%File_Latt),trim(File_suff)
#endif
              write(File_aux, '(A,A)') trim(File_pr), "_info"
              inquire(file=File_aux, exist=File_exists)
              if (.not.File_exists) then
                 11 format(A20, ': ', A)
                 12 format(A20, ': ', I10)
                 13 format(A20, ': ', *(E26.17E3))
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
                    write(10, 13) trim(tmp_str), Obs%Latt_unit%Orb_pos_p(no,:)
                 enddo
                 close(10)
              endif

              do nt = 1, Ntau
                 do no = 1, Obs%Latt_unit%Norb
                    do no1 = 1, Obs%Latt_unit%Norb
                       Call Fourier_R_to_K(Obs%Obs_Latt(:,nt,no,no1), Tmp(:,nt,no,no1), Obs%Latt)
                    enddo
                 enddo
              enddo
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
#if defined(MPI)
           Endif
#endif

           deallocate (Tmp)

         End Subroutine Print_bin_Latt

!--------------------------------------------------------------------

         Subroutine  Print_bin_Vec(Obs,Group_Comm)
#ifdef MPI
           Use mpi
#endif
           Implicit none

           Type (Obser_vec), intent(Inout) :: Obs
           Integer, INTENT(IN)  :: Group_Comm

           ! Local
           Integer :: I
           Character (len=64) :: File_pr, File_suff, File_aux
           logical            :: File_exists
#ifdef MPI
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
           File_suff = "_scal"
           write(File_pr, '(A,A)') trim(Obs%File_Vec), Trim(File_suff)

#if defined(MPI)
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
#if defined(TEMPERING)
              write(File_pr,'(A,I0,A,A,A)') "Temp_",igroup,"/",trim(Obs%File_Vec),trim(File_suff)
#endif
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
#if defined(MPI)
           endif
#endif

         End Subroutine Print_bin_Vec

       end Module Observables
