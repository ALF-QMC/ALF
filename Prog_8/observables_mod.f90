!  Copyright (C) 2016 The ALF project
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

     Module Observables

       Use Files_mod
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief 
!> This module defines the Obser_Vec and Obser_Latt types and provides  
!> routine to initialize them and to print out the bins
!
!--------------------------------------------------------------------
    
       Type Obser_Vec
!>  Data structure for 
!>  < O_n >  n : =1, size(Obs,1)  
          Integer            :: N                    ! Number of measurements
          real      (Kind=Kind(0.d0)) :: Ave_Sign             ! Averarge sign
          complex   (Kind=Kind(0.d0)), pointer :: Obs_vec(:)  ! Vector of observables
          Character (len=64) :: File_Vec             ! Name of file in which the bins will be written out
       end type Obser_Vec
       

       Type Obser_Latt
!>  Data structure for 
!>  < O^{dagger}(i,tau)_n O(j,0)_m>  - < O_n> <O_m> 
!>  where it is assumed that translation symmetry as specified by the lattice Latt is present.  
!>  Obs_Latt(i-j,tau,n,m) = < O^{dagger}(i,tau)_n O(j,0)_m>  
!>  Obs_Latt0(n) = < O_n>
!>  For equal   time correlation functions, tau runs from 1,1 
!>  For unequal time correlation functions, tau runs from 1,Ltrot+1  
          Integer            :: N                           ! Number of measurements
          Real      (Kind=Kind(0.d0)) :: Ave_Sign                    ! Averarge sign
          complex   (Kind=Kind(0.d0)), pointer :: Obs_Latt (:,:,:,:) ! i-j, tau, norb, norb  
          complex   (Kind=Kind(0.d0)), pointer :: Obs_Latt0(:)       ! norb 
          Character (len=64) :: File_Latt                   ! Name of file in which the bins will be written out
       end type Obser_Latt
       


       Contains

         Subroutine Obser_Latt_make(Obs,Ns,Nt,No,Filename)
           Implicit none
           Type (Obser_Latt), intent(INOUT) :: Obs
           Integer, Intent(IN)             :: Ns,Nt,No
           Character (len=64), Intent(IN)  :: Filename 
           Allocate (Obs%Obs_Latt (Ns,Nt,No,No))
           Allocate (Obs%Obs_Latt0(No)         )
           Obs%File_Latt = Filename
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

         Subroutine Obser_Vec_make(Obs,N,Filename)
           Implicit none
           Type (Obser_vec), intent(INOUT) :: Obs
           Integer, Intent(IN)             :: N
           Character (len=64), Intent(IN)  :: Filename 
           Allocate (Obs%Obs_vec(N))
           Obs%File_Vec = Filename
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
         
         Subroutine  Print_bin_Latt(Obs,Latt,dtau)
           Use Lattices_v3
           Implicit none

#ifdef MPI
           include 'mpif.h'
#endif   

           Type (Obser_Latt),        Intent(Inout)   :: Obs
           Type (Lattice),           Intent(In)      :: Latt
           Real (Kind=Kind(0.d0)),            Intent(In)      :: dtau

           ! Local
           Integer :: Ns,Nt, Norb, no, no1, I , Ntau
           Complex (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:,:), Tmp1(:)
           Real    (Kind=Kind(0.d0))              :: x_p(2) 
           Complex (Kind=Kind(0.d0))              :: Sign_bin
           Character (len=64)            :: File_pr, File_suff
#ifdef MPI
           Complex (Kind=Kind(0.d0)):: Z
           Real    (Kind=Kind(0.d0)):: X
           Integer         :: Ierr, Isize, Irank
           INTEGER         :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

           Ns    = Size(Obs%Obs_Latt,1)
           Ntau  = Size(Obs%Obs_Latt,2)
           Norb  = Size(Obs%Obs_Latt,3)
           if ( .not. (Latt%N  == Ns ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           If (Ntau == 1) then
              File_suff ="_eq"
           else
              File_suff ="_tau"
           endif
           File_pr = file_add(Obs%File_Latt,File_suff)
           Allocate (Tmp(Ns,Ntau,Norb,Norb), Tmp1(Norb) )
           Obs%Obs_Latt  =   Obs%Obs_Latt /dble(Obs%N   )
           Obs%Obs_Latt0 =   Obs%Obs_Latt0/dble(Obs%N*Ns*Ntau)
           Obs%Ave_sign  =   Obs%Ave_Sign /dble(Obs%N   )

#ifdef MPI
           I = Ns*Ntau*Norb*Norb
           Tmp = cmplx(0.d0, 0.d0, kind(0.D0))
           CALL MPI_REDUCE(Obs%Obs_Latt,Tmp,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Obs_Latt = Tmp/DBLE(ISIZE)

           I = 1
           X = 0.d0
           CALL MPI_REDUCE(Obs%Ave_sign,X,I,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Ave_sign = X/DBLE(ISIZE)

           I = Norb
           Tmp1 = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Obs%Obs_Latt0,Tmp1,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Obs_Latt0 = Tmp1/DBLE(ISIZE)

           If (Irank == 0 ) then
#endif
              do nt = 1,Ntau
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Call  Fourier_R_to_K(Obs%Obs_Latt(:,nt,no,no1), Tmp(:,nt,no,no1), Latt)
                    enddo
                 enddo
              enddo
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              If ( Ntau == 1 ) then
                 Write(10,*) Obs%Ave_sign,Norb,Latt%N
              else
                 Write(10,*) Obs%Ave_sign,Norb,Latt%N, Ntau, dtau
              endif
              Do no = 1,Norb
                 Write(10,*)  Obs%Obs_Latt0(no)
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
                 Write(10,*) X_p(1), X_p(2)
                 Do nt = 1,Ntau
                    do no = 1,Norb
                       do no1 = 1,Norb
                          Write(10,*) tmp(I,nt,no,no1)
                       enddo
                    enddo
                 enddo
              enddo
              close(10)
#ifdef MPI
           Endif
#endif
              
           deallocate (Tmp, tmp1 )
          

         End Subroutine Print_bin_Latt
         
!--------------------------------------------------------------------

         Subroutine  Print_bin_Vec(Obs)
           
           Implicit none

#ifdef MPI
           include 'mpif.h'
#endif   
           
           Type (Obser_vec), intent(Inout) :: Obs
           
           ! Local
           Integer :: No,I
           Character (len=64)             :: File_pr, File_suff
           
           
#ifdef MPI
           Integer        :: Ierr, Isize, Irank
           INTEGER        :: STATUS(MPI_STATUS_SIZE)
           Complex  (Kind=Kind(0.d0)), allocatable :: Tmp(:)
           Real     (Kind=Kind(0.d0)) :: X
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
           
           No = size(Obs%Obs_vec,1)
           Obs%Obs_vec  = Obs%Obs_vec /dble(Obs%N)
           Obs%Ave_sign = Obs%Ave_sign/dble(Obs%N)
           File_suff ="_scal"
           File_pr = file_add(Obs%File_Vec,File_suff)


#ifdef MPI
           Allocate (Tmp(No) )
           Tmp = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Obs%Obs_vec,Tmp,No,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Obs_vec = Tmp/DBLE(ISIZE)
           deallocate (Tmp )

           I = 1
           X = 0.d0
           CALL MPI_REDUCE(Obs%Ave_sign,X,I,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Ave_sign = X/DBLE(ISIZE)

           if (Irank == 0 ) then
#endif
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              WRITE(10,*) size(Obs%Obs_vec,1)+1, (Obs%Obs_vec(I), I=1,size(Obs%Obs_vec,1)), Obs%Ave_sign
              close(10)
#ifdef MPI
           endif
#endif
           
         End Subroutine Print_bin_Vec

         
       end Module Observables
