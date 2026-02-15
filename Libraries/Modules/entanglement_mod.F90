!  Copyright (C) 2020 The ALF project
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
 
Module entanglement_mod

#ifndef MPI
#warning "You are compiling entanglement without MPI. No Renyi entropy results possible, all other observables still work!"
#endif

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This module generates one and two dimensional Bravais lattices.
!
!--------------------------------------------------------------------
      ! Used for MPI
      INTEGER, save, private :: ENTCOMM, ENT_RANK, ENT_SIZE=0, Norm, group
      Real (kind=kind(0.d0)), save, private :: weight

      INTERFACE Calc_Renyi_Ent
        !> Interface to Calc_Renyi_Ent function.
        MODULE PROCEDURE Calc_Renyi_Ent_gen_all, Calc_Renyi_Ent_indep, Calc_Renyi_Ent_gen_fl
      END INTERFACE
      INTERFACE Calc_Mutual_Inf
        !> Interface to Calc_Mutual_Inf Subroutine.
        MODULE PROCEDURE Calc_Mutual_Inf_indep, Calc_Mutual_Inf_gen_fl, Calc_Mutual_Inf_gen_all
      END INTERFACE
      Contains
!========================================================================

        Subroutine Init_Entanglement_replicas(Group_Comm)
#ifdef MPI
          Use mpi
#endif
          Implicit none
          Integer, INTENT(IN)               :: Group_Comm
          
          Integer                           :: ISIZE, IRANK, IERR
#ifdef MPI
          ! Create subgroups of two replicas each
          CALL MPI_COMM_SIZE(Group_Comm,ISIZE,IERR)
          CALL MPI_COMM_RANK(Group_Comm,IRANK,IERR)
          CALL MPI_COMM_SPLIT(Group_Comm, IRANK / 2, IRANK, ENTCOMM, IERR)
          CALL MPI_COMM_RANK(ENTCOMM,ENT_RANK,IERR)
          CALL MPI_COMM_SIZE(ENTCOMM,ENT_SIZE,IERR)
          Norm=ISIZE/2 ! number of pairs
          norm=2*norm ! but each task of pair contributes
          group=Group_Comm
          weight=dble(ISIZE)/dble(Norm)
#endif
          
        end Subroutine Init_Entanglement_replicas
          
!========================================================================
        ! Calculation of the Renyi entanglement entropy
        ! The algorithm works only for an MPI program
        ! We partition the nodes into groups of 2 replicas:
        ! ! (n, n+1), with n=0,2,...
        Subroutine Calc_Mutual_Inf_indep(GRC,List_A,Nsites_A,List_B,Nsites_B,N_SUN,Renyi_A,Renyi_B,Renyi_AB)

          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:), INTENT(IN) :: List_A, List_B
          Integer, INTENT(IN)               :: Nsites_A ,Nsites_B, N_SUN
          Complex (kind=kind(0.d0)), INTENT(OUT)   :: Renyi_A, Renyi_B, Renyi_AB

          Integer, Dimension(:), Allocatable :: List_AB
          Integer          :: I, J, IERR, INFO, Nsites_AB

          Nsites_AB=Nsites_A+Nsites_B
          
          allocate(List_AB(Nsites_AB))
          
          DO I = 1, Nsites_A
             List_AB(I) = List_A(I)
          END DO
          DO I = 1, Nsites_B
             List_AB(I+Nsites_A) = List_B(I)
          END DO
          
          Renyi_A  = Calc_Renyi_Ent_indep(GRC,List_A,Nsites_A,N_SUN)
          Renyi_B  = Calc_Renyi_Ent_indep(GRC,List_B,Nsites_B,N_SUN)
          Renyi_AB = Calc_Renyi_Ent_indep(GRC,List_AB,Nsites_AB,N_SUN)
          
          deallocate(List_AB)
          
        End Subroutine Calc_Mutual_Inf_indep
          
!========================================================================
        ! Calculation of the Renyi entanglement entropy
        ! The algorithm works only for an MPI program
        ! We partition the nodes into groups of 2 replicas:
        ! ! (n, n+1), with n=0,2,...
        Subroutine Calc_Mutual_Inf_gen_fl(GRC,List_A,Nsites_A,List_B,Nsites_B,N_SUN,Renyi_A,Renyi_B,Renyi_AB)

          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)    :: GRC(:,:,:)
          Integer, Dimension(:,:), INTENT(IN)      :: List_A, List_B
          Integer, Dimension(:), INTENT(IN)        :: Nsites_A ,Nsites_B, N_SUN
          Complex (kind=kind(0.d0)), INTENT(OUT)   :: Renyi_A, Renyi_B, Renyi_AB

          Integer, Allocatable :: List_AB(:,:), Nsites_AB(:)
          Integer          :: I, J, IERR, INFO, N_FL, Nsites_AB_max

          N_FL=size(GRC,3)
          Allocate(Nsites_AB(N_FL))
          Nsites_AB_max=0
          do I=1,N_FL
            Nsites_AB(I)=Nsites_A(I)+Nsites_B(I)
            if (Nsites_AB(I)>Nsites_AB_max) Nsites_AB_max=Nsites_AB(I)
          enddo
          
          allocate(List_AB(Nsites_AB_max,N_FL))
          
          do J=1,N_FL
            DO I = 1, Nsites_A(J)
              List_AB(I,J) = List_A(I,J)
            END DO
            DO I = 1, Nsites_B(J)
              List_AB(I+Nsites_A(J),J) = List_B(I,J)
            END DO
          enddo
          
          Renyi_A  = Calc_Renyi_Ent_gen_fl(GRC,List_A,Nsites_A,N_SUN)
          Renyi_B  = Calc_Renyi_Ent_gen_fl(GRC,List_B,Nsites_B,N_SUN)
          Renyi_AB = Calc_Renyi_Ent_gen_fl(GRC,List_AB,Nsites_AB,N_SUN)
          
          deallocate(List_AB,Nsites_AB)
          
        End Subroutine Calc_Mutual_Inf_gen_fl
          
!========================================================================
        ! Calculation of the Renyi entanglement entropy
        ! The algorithm works only for an MPI program
        ! We partition the nodes into groups of 2 replicas:
        ! ! (n, n+1), with n=0,2,...
        Subroutine Calc_Mutual_Inf_gen_all(GRC,List_A,Nsites_A,List_B,Nsites_B,Renyi_A,Renyi_B,Renyi_AB)

          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)    :: GRC(:,:,:)
          Integer, Dimension(:,:,:), INTENT(IN)      :: List_A, List_B
          Integer, Dimension(:,:), INTENT(IN)        :: Nsites_A ,Nsites_B
          Complex (kind=kind(0.d0)), INTENT(OUT)   :: Renyi_A, Renyi_B, Renyi_AB

          Integer, Allocatable :: List_AB(:,:,:), Nsites_AB(:,:)
          Integer          :: I, J, IERR, INFO, N_FL, Nsites_AB_max, nc, num_nc

          N_FL=size(GRC,3)
          num_nc=size(List_B,3)
          Allocate(Nsites_AB(N_FL,num_nc))
          Nsites_AB_max=0
          do nc=1,num_nc
            do I=1,N_FL
              Nsites_AB(I,nc)=Nsites_A(I,nc)+Nsites_B(I,nc)
              if (Nsites_AB(I,nc)>Nsites_AB_max) Nsites_AB_max=Nsites_AB(I,nc)
            enddo
          enddo
          
          allocate(List_AB(Nsites_AB_max,N_FL,num_nc))
          
          do nc=1, num_nc
            do J=1,N_FL
              DO I = 1, Nsites_A(J,nc)
                List_AB(I,J,nc) = List_A(I,J,nc)
              END DO
              DO I = 1, Nsites_B(J,nc)
                List_AB(I+Nsites_A(J,nc),J,nc) = List_B(I,J,nc)
              END DO
            enddo
          enddo
          
          Renyi_A  = Calc_Renyi_Ent_gen_all(GRC,List_A,Nsites_A)
          Renyi_B  = Calc_Renyi_Ent_gen_all(GRC,List_B,Nsites_B)
          Renyi_AB = Calc_Renyi_Ent_gen_all(GRC,List_AB,Nsites_AB)
          
          deallocate(List_AB,Nsites_AB)
          
        End Subroutine Calc_Mutual_Inf_gen_all

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Compute the Renyi entropy for patches of the same N sites for each 
!> flavor and color degree of freedom. Returns expontiated Renyi entropy.
!> @details
!> @param [IN] GRC Complex(:,:,:)
!> \verbatim
!>  Greens function
!> \endverbatim
!> @param [IN] List  Integer(:)
!> \verbatim
!>  List of sites that each patch will use.
!> \endverbatim
!> @param [IN] Nsites   Integer
!> \verbatim
!>  Number of sites in each patch.
!> \endverbatim
!> @param [IN] N_SUN   Integer
!> \verbatim
!>  Number of color degrees of freedom.
!> \endverbatim
!-------------------------------------------------------------------

          Complex (kind=kind(0.d0)) function Calc_Renyi_Ent_indep(GRC,List,Nsites,N_SUN)
#ifdef MPI
          Use mpi
#endif
          use iso_fortran_env, only: error_unit
  
          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)      :: GRC(:,:,:)
          Integer, INTENT(IN)               :: List(:)
          Integer, INTENT(IN)               :: Nsites, N_SUN

          Complex (kind=kind(0.d0)), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          ! Integer, Dimension(:), Allocatable :: PIVOT
          Complex (kind=kind(0.d0)) :: DET, PRODDET, alpha, beta
          Integer          :: J, IERR, INFO, N_FL, nf, N_FL_half
          Integer         , Dimension(:,:), Allocatable :: List_tmp
          Integer         , Dimension(2)              :: Nsites_tmp,nf_list,N_SUN_tmp
          Logical, save :: First_call = .True.
          
          EXTERNAL ZGEMM
          EXTERNAL ZGETRF
          
          N_FL = size(GRC,3)
          N_FL_half = N_FL/2
            
          Calc_Renyi_Ent_indep=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))

          if (Nsites==0) then
            Calc_Renyi_Ent_indep=0.d0
            return
          endif
            
#ifdef MPI
          ! Check if entanglement replica group is of size 2 such that the second renyi entropy can be calculated
          if(ENT_SIZE==2) then
          
            allocate(List_tmp(NSITES,2))
            Allocate(GreenA(Nsites,2*Nsites),GreenA_tmp(Nsites,2*Nsites),IDA(Nsites,Nsites)) ! new

            DO J = 1, 2
              List_tmp(:,J)=List(:)
              Nsites_tmp(J) = Nsites
              N_SUN_tmp(J) = N_SUN
            enddo

            DO nf=1,N_FL_half
              
              DO J = 1, 2
                nf_list(J) = 2*nf-2+J
              enddo
              PRODDET = Calc_Renyi_Ent_pair(GRC,List_tmp,Nsites_tmp,nf_list,N_SUN_tmp,GreenA, GreenA_tmp, IDA)
              Calc_Renyi_Ent_indep = Calc_Renyi_Ent_indep * PRODDET
              
            Enddo
              
            if (N_FL/=2*N_FL_half) then
              PRODDET = Calc_Renyi_Ent_single(GRC,List_tmp(:,1),Nsites,N_fl,N_SUN,GreenA, GreenA_tmp, IDA)
              Calc_Renyi_Ent_indep = Calc_Renyi_Ent_indep * PRODDET
            
            endif
            
            Deallocate(GreenA,GreenA_tmp,IDA,List_tmp)
            
          else
            ! if there had been an odd number of task in tempering group / world, set renyi to 0
            Calc_Renyi_Ent_indep=CMPLX(0.d0,0.d0,kind(0.d0))
          endif
          
          ! average over all pairs of replicas, the single task contributes nothing even so it takes part in the call
          Calc_Renyi_Ent_indep=Calc_Renyi_Ent_indep*weight
          ! At this point, each task of the temepering group / world returns the same averaged value of the pairs, including the possible "free"/ unpaired one.
          ! This mechanisms leads to some syncronization, but I (Johannes) am lacking a better way to treat odd number of tasks.
#else
          Calc_Renyi_Ent_indep=0.0d0
          if (First_call) then
            write(error_unit,*) "Entanglement module compiled without MPI, no Renyi entropy results possible!"
            First_call = .false.
          endif
#endif
              
          End function Calc_Renyi_Ent_indep

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Compute the Renyi entropy for a list of site patches that may differ for
!> each flavor, but each given flavor patch f will be the same across
!> NSUN(f) color degrees of freedom. Returns expontiated Renyi entropy.
!> @details
!> @param [IN] GRC Complex(:,:,:)
!> \verbatim
!>  Greens function
!> \endverbatim
!> @param [IN] List  Integer(:,:)
!> \verbatim
!>  List(:,f) gives the list of sites that the flavor patch i contains.
!> \endverbatim
!> @param [IN] Nsites(:)   Integer
!> \verbatim
!>  Number of sites in each flavor patch.
!> \endverbatim
!> @param [IN] N_SUN(:)   Integer
!> \verbatim
!>  Number of color degrees of freedom for each flavor patch.
!> \endverbatim
!-------------------------------------------------------------------        
        
        Complex (kind=kind(0.d0)) function Calc_Renyi_Ent_gen_fl(GRC,List,Nsites,N_SUN)
#ifdef MPI
          Use mpi
#endif
          use iso_fortran_env, only: error_unit

          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)      :: GRC(:,:,:)
          !Integer, Dimension(:,:), INTENT(IN) :: List ! new
          Integer, INTENT(IN) :: List(:,:)
          Integer, INTENT(IN)               :: Nsites(:), N_SUN(:) ! new

          Complex (kind=kind(0.d0)), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          ! Integer, Dimension(:), Allocatable :: PIVOT
          Complex (kind=kind(0.d0)) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, nf_eff, start_flav
          Integer         , Dimension(:), Allocatable :: SortedFlavors ! new
          Integer         , Dimension(:,:), Allocatable :: List_tmp
          Integer         , Dimension(2)              :: Nsites_tmp,nf_list,N_SUN_tmp
          Logical, save :: First_call = .True.
          
          EXTERNAL ZGEMM
          EXTERNAL ZGETRF
          
          N_FL = size(GRC,3)

          Allocate(SortedFlavors(N_FL)) ! new

          ! insertion sort for small number of elements, SortedFlavors lists flavors in order of size
          start_flav=0
          if (Nsites(1)==0) start_flav = 1
          SortedFlavors(1) = 1
          DO I=2,N_FL
            x = Nsites(I)
            if (x==0) start_flav = start_flav + 1
            J = I-1
            DO while(J >= 1)
              if(Nsites(J) <= x) exit
              SortedFlavors(J+1) = J 
              J = J - 1
            end do
            SortedFlavors(J+1) = I
          END DO
          if(start_flav==N_FL) then
            Calc_Renyi_Ent_gen_fl=0.0d0
            return
          endif
          N_FL_half = (N_FL-start_flav)/2
          
          Calc_Renyi_Ent_gen_fl=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
          
#ifdef MPI
          ! Check if entanglement replica group is of size 2 such that the second reny entropy can be calculated
          if(ENT_SIZE==2) then
          
            !Allocate(GreenA(Nsites,2*Nsites),GreenA_tmp(Nsites,2*Nsites),IDA(Nsites,Nsites))
            dim = Nsites(SortedFlavors(N_FL)) ! new
            allocate(List_tmp(dim,2))
            Allocate(GreenA(dim,2*dim),GreenA_tmp(dim,2*dim),IDA(dim,dim)) ! new

            DO nf=1,N_FL_half

              DO J = 1, 2
                nf_eff = SortedFlavors(start_flav+2*nf-2+J)
                Nsites_tmp(J)=Nsites(nf_eff)
                List_tmp(:,J)=List(1:dim,nf_eff)
                N_sun_tmp(J)=N_SUN(nf_eff)
                nf_list(J)=nf_eff
              enddo
              PRODDET = Calc_Renyi_Ent_pair(GRC,List_tmp,Nsites_tmp,nf_list,N_SUN_tmp,GreenA, GreenA_tmp, IDA)
              Calc_Renyi_Ent_gen_fl = Calc_Renyi_Ent_gen_fl * PRODDET
              
            Enddo
              
            if (N_FL/=2*N_FL_half+start_flav) then

              nf_eff = SortedFlavors(N_fl)
              List_tmp(:,1)=List(1:dim,nf_eff)
            
              PRODDET = Calc_Renyi_Ent_single(GRC,List_tmp(:,1),Nsites(nf_eff),nf_eff,N_SUN(nf_eff),GreenA, GreenA_tmp, IDA)
              Calc_Renyi_Ent_gen_fl = Calc_Renyi_Ent_gen_fl * PRODDET
            
            endif
            
            Deallocate(GreenA,GreenA_tmp,IDA,List_tmp)
            
          else
            ! if there had been an odd number of task in tempering group / world, set renyi to 0
            Calc_Renyi_Ent_gen_fl=CMPLX(0.d0,0.d0,kind(0.d0))
          endif
          
          ! average over all pairs of replicas, the single task contributes nothing even so it takes part in the call
          Calc_Renyi_Ent_gen_fl=Calc_Renyi_Ent_gen_fl*weight
          ! At this point, each task of the temepering group / world returns the same averaged value of the pairs, including the possible "free"/ unpaired one.
          ! This mechanisms leads to some syncronization, but I (Johannes) am lacking a better way to treat odd number of tasks.
#else
          Calc_Renyi_Ent_gen_fl=0.0d0
          if (First_call) then
            write(error_unit,*) "Entanglement module compiled without MPI, no Renyi entropy results possible!"
            First_call = .false.
          endif
#endif
            
        End function Calc_Renyi_Ent_gen_fl


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Compute the Renyi entropy for a generic set of site patches that may differ
!> for each flavor and color degree of freedom combination.
!> Returns expontiated Renyi entropy.
!> @details
!> @param [IN] GRC Complex(:,:,:)
!> \verbatim
!>  Greens function
!> \endverbatim
!> @param [IN] List  Integer(:,:,:)
!> \verbatim
!>  List(:,f,c) gives the list of sites that the patch with flavor index f and
!>  color index c contains.
!> \endverbatim
!> @param [IN] Nsites(:,:)   Integer
!> \verbatim
!>  Nsites(f,c) gives the number of sites in the patch for flavor f and color c.
!> \endverbatim
!-------------------------------------------------------------------   

        Complex (kind=kind(0.d0)) function Calc_Renyi_Ent_gen_all(GRC,List,Nsites)
#ifdef MPI
          Use mpi
#endif
          use iso_fortran_env, only: error_unit

          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:,:,:), INTENT(IN) :: List
          Integer, INTENT(IN)               :: Nsites(:,:)

          Complex (kind=kind(0.d0)), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          ! Integer, Dimension(:), Allocatable :: PIVOT
          Complex (kind=kind(0.d0)) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, nf_eff, start_flav
          Integer          :: nc, num_nc
          Integer         , Dimension(:), Allocatable :: SortedFlavors,N_SUN_fl,df_list
          Integer         , Dimension(:,:), Allocatable :: List_tmp, eff_ind, eff_ind_inv
          Integer         , Dimension(2)              :: Nsites_tmp,nf_list,N_SUN_tmp
          Logical, save :: First_call = .True.
          
          EXTERNAL ZGEMM
          EXTERNAL ZGETRF
          
          N_FL = size(GRC,3)
          num_nc = size(List,3)
          Allocate(SortedFlavors(num_nc*N_FL),N_SUN_fl(num_nc*N_FL),eff_ind(N_FL,num_nc),eff_ind_inv(2,N_FL*num_nc))

          I=0
          do nc=1,num_nc
            do nf=1,N_FL
              I=I+1
              eff_ind(nf,nc)=i
              eff_ind_inv(1,I)=nf
              eff_ind_inv(2,I)=nc
            enddo
          enddo
          N_SUN_fl=1

          ! insertion sort for small number of elements, SortedFlavors lists flavors in order of size
          start_flav=0
          if (Nsites(1,1)==0) start_flav = 1
          SortedFlavors(1) = 1
          ! might have an update in the future to exchange color and flavor loops--optimization
          DO I=2,N_FL*num_nc
            x = Nsites(eff_ind_inv(1,I),eff_ind_inv(2,I))
            if (x==0) start_flav = start_flav + 1
            J = I-1
            DO while(J >= 1)
              if(Nsites(eff_ind_inv(1,J),eff_ind_inv(2,J)) <= x) exit
              SortedFlavors(J+1) = J 
              J = J - 1
            end do
            SortedFlavors(J+1) = I
          END DO

          if(start_flav==N_FL*num_nc) then
            Calc_Renyi_Ent_gen_all=0.0d0
            return
          endif

          N_FL_half = (N_FL*num_nc-start_flav)/2
          
          Calc_Renyi_Ent_gen_all=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
            
#ifdef MPI
          ! Check if entanglement replica group is of size 2 such that the second reny entropy can be calculated
          if(ENT_SIZE==2) then
          
            !Allocate(GreenA(Nsites,2*Nsites),GreenA_tmp(Nsites,2*Nsites),IDA(Nsites,Nsites))
            nf=eff_ind_inv(1,SortedFlavors(N_FL*num_nc))
            nc=eff_ind_inv(2,SortedFlavors(N_FL*num_nc))
            dim = Nsites(nf,nc) ! new
            allocate(List_tmp(dim,2))
            Allocate(GreenA(dim,2*dim),GreenA_tmp(dim,2*dim),IDA(dim,dim)) ! new

            DO I=1,N_FL_half

              DO J = 1, 2
                nf=eff_ind_inv(1,SortedFlavors(start_flav+2*I-2+J))
                nc=eff_ind_inv(2,SortedFlavors(start_flav+2*I-2+J))
                Nsites_tmp(J)=Nsites(nf,nc)
                List_tmp(:,J)=List(1:dim,nf,nc)
                N_sun_tmp(J)=1
                nf_list(J)=nf
              enddo
              PRODDET = Calc_Renyi_Ent_pair(GRC,List_tmp,Nsites_tmp,nf_list,N_SUN_tmp,GreenA, GreenA_tmp, IDA)
              Calc_Renyi_Ent_gen_all = Calc_Renyi_Ent_gen_all * PRODDET
              
            Enddo
              
            if (N_FL*num_nc/=2*N_FL_half+start_flav) then

              nf=eff_ind_inv(1,SortedFlavors(N_FL*num_nc))
              nc=eff_ind_inv(2,SortedFlavors(N_FL*num_nc))
              List_tmp(:,1)=List(1:dim,nf,nc)
            
              proddet = Calc_Renyi_Ent_single(GRC,List_tmp(:,1),Nsites(nf,nc),nf,1,GreenA, GreenA_tmp, IDA)
              Calc_Renyi_Ent_gen_all = Calc_Renyi_Ent_gen_all * PRODDET
            
            endif
            
            Deallocate(GreenA,GreenA_tmp,IDA,List_tmp)
              
          else
            ! if there had been an odd number of task in tempering group / world, set renyi to 0
            Calc_Renyi_Ent_gen_all=CMPLX(0.d0,0.d0,kind(0.d0))
          endif
            
          ! average over all pairs of replicas, the single task contributes nothing even so it takes part in the call
          Calc_Renyi_Ent_gen_all=Calc_Renyi_Ent_gen_all*weight

          ! At this point, each task of the temepering group / world returns the same averaged value of the pairs, including the possible "free"/ unpaired one.
          ! This mechanisms leads to some syncronization, but I (Johannes) am lacking a better way to treat odd number of tasks.
#else
          Calc_Renyi_Ent_gen_all=0.0d0
          if (First_call) then
            write(error_unit,*) "Entanglement module compiled without MPI, no Renyi entropy results possible!"
            First_call = .false.
          endif
#endif

            
        End function Calc_Renyi_Ent_gen_all
       
       
#ifdef MPI
! The following two helper function cannot do meaningfull calculatios without MPI.
! Hence they are only declared and defined when MPI is present.
! They are designed as helper functions within this module and won't be called here if MPI isn't defined
! However, they are currently still public such that a user may call them directly,
! with the neccessary MPI protection clauses.


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Compute the Renyi entropy for a generic pair of patches.
!> Returns expontiated Renyi entropy.
!> @details
!> @param [IN] GRC Complex(:,:,:)
!> \verbatim
!>  Greens function
!> \endverbatim
!> @param [IN] List  Integer(:,:)
!> \verbatim
!>  List(:,1:2) gives the list of sites of the first/second patch of the pair.
!> \endverbatim
!> @param [IN] Nsites(2)   Integer
!> \verbatim
!>  Nsites(1:2) gives the number of sites in the patches of the pair.
!> \endverbatim
!> @param [IN] nf_list(2)   Integer
!> \verbatim
!>  nf_list(1:2) gives the flavor of the patches in the pair.
!> \endverbatim
!> @param [IN] N_SUN(2)   Integer
!> \verbatim
!>  N_SUN(1:2) gives the number of colors of the patches in the pair.
!> \endverbatim
!> @param [IN] GreenA Complex(:,:)
!> \verbatim
!>  temporary Greens function storage of dimension (dim,2*dim).
!>  Dim has to be larger or equal to the larger patch size.
!> \endverbatim
!> @param [IN] GreenA_tmp Complex(:,:)
!> \verbatim
!>  same as GreenA.
!> \endverbatim
!> @param [IN] IDA Complex(:,:)
!> \verbatim
!>  temporary Greens function storage of dimension (dim,dim).
!>  Dim has to be the same as in GreenA and GreenA_tmp.
!> \endverbatim
!-------------------------------------------------------------------   
        Complex (kind=kind(0.d0)) function Calc_Renyi_Ent_pair(GRC,List,Nsites,nf_list,N_SUN,GreenA, GreenA_tmp, IDA)
          Use mpi

          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:,:), INTENT(IN) :: List ! new
          Integer, INTENT(IN)               :: Nsites(2), N_SUN(2),nf_list(2) ! new
          Complex (kind=kind(0.d0)), INTENT(OUT), Dimension(:,:) :: GreenA, GreenA_tmp, IDA

          Integer, Dimension(:), Allocatable :: PIVOT
          Complex (kind=kind(0.d0)) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, nf_eff, start_flav, dim_sq
          Integer         , Dimension(:), Allocatable :: SortedFlavors ! new


          Calc_Renyi_Ent_pair=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
          
          dim=size(IDA,1)
          dim_sq=dim*dim
          ! We store the reduced Green's function in GreenA_c_tmp (c electrons)
          ! and GreenA_f_tmp (f electrons)
          ! The first Nsites columns contains the spin up sector,
          ! the last Nsites column the spin down sector
          nf_eff = nf_list(1)

          DO J = 1, Nsites(1)
            Do I = 1, Nsites(1)
              GreenA_tmp(I,J) = GRC(List(I,1),List(J,1),nf_eff)
            END DO
          END Do

          nf_eff = nf_list(2)
          DO J = 1, Nsites(2)
            DO I = 1, Nsites(2)
                GreenA_tmp(I,J+dim) = GRC(List(I,2), List(J,2), nf_eff)
            END DO
          END DO
          ! This exchange the last Nsites columns of GreenA_c_tmp between the two replicas
          ! such that GreenA contains the reduced Green's function for two replicas
          ! and a fixed spin sector.
          CALL MPI_ALLTOALL(GreenA_tmp, dim_sq, MPI_COMPLEX16, GreenA, dim_sq, MPI_COMPLEX16, ENTCOMM, IERR)

          ! Compute Identity - GreenA(replica=1) - GreenA(replica=2) + 2 GreenA(replica=1) * GreenA(replica=2)
          dim_eff = Nsites(1+ENT_RANK)
          IDA(1:dim_eff,1:dim_eff) = - GreenA(1:dim_eff, 1:dim_eff) - GreenA(1:dim_eff, dim+1:dim+dim_eff)
          DO I = 1, dim_eff
              IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
          END DO
          CALL ZGEMM('n', 'n', dim_eff, dim_eff, dim_eff, alpha, GreenA(1, 1), &
              & dim, GreenA(1, dim+1), dim, beta, IDA, dim)
          ! Compute determinant
          SELECT CASE(dim_eff)
          CASE (1)
            DET = IDA(1,1)
          CASE (2)
            DET = IDA(1,1) * IDA(2,2) - IDA(1,2) * IDA(2,1)
          CASE DEFAULT
            Allocate(PIVOT(dim_eff))
            CALL ZGETRF(dim_eff, dim_eff, IDA, dim, PIVOT, INFO)
            DET = cmplx(1.D0,0.D0,KIND(0.D0))
            DO I = 1, dim_eff
                IF (PIVOT(I).NE.I) THEN
                  DET = -DET * IDA(I,I)
                ELSE
                  DET = DET * IDA(I,I)
                END IF
            ENDDO
            Deallocate(PIVOT)
          END SELECT
          DET=DET**N_SUN(1+ENT_RANK)
          ! Compute the product of determinants for up and down spin sectors.
          CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
          ! Now each thread contains in PRODDET the full determinant, as obtained by
          ! a pair of replicas.
          Calc_Renyi_Ent_pair = Calc_Renyi_Ent_pair * PRODDET
        end function Calc_Renyi_Ent_pair

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Compute the Renyi entropy for a single unpaired patches.
!> Returns expontiated Renyi entropy.
!> @details
!> @param [IN] GRC Complex(:,:,:)
!> \verbatim
!>  Greens function
!> \endverbatim
!> @param [IN] List  Integer(:)
!> \verbatim
!>  List(:,1:2) gives the list of sites of the patch.
!> \endverbatim
!> @param [IN] Nsites   Integer
!> \verbatim
!>  Nsites gives the number of sites in the patch.
!> \endverbatim
!> @param [IN] nf_list   Integer
!> \verbatim
!>  nf_list gives the flavor of the patch.
!> \endverbatim
!> @param [IN] N_SUN   Integer
!> \verbatim
!>  N_SUN gives the number of colors of the patch.
!> \endverbatim
!> @param [IN] GreenA Complex(:,:)
!> \verbatim
!>  temporary Greens function storage of dimension (dim,2*dim).
!>  Dim has to be larger or equal to the larger patch size.
!> \endverbatim
!> @param [IN] GreenA_tmp Complex(:,:)
!> \verbatim
!>  same as GreenA.
!> \endverbatim
!> @param [IN] IDA Complex(:,:)
!> \verbatim
!>  temporary Greens function storage of dimension (dim,dim).
!>  Dim has to be the same as in GreenA and GreenA_tmp.
!> \endverbatim
!-------------------------------------------------------------------   
        Complex (Kind=8) function Calc_Renyi_Ent_single(GRC,List,Nsites,nf_eff,N_SUN,GreenA, GreenA_tmp, IDA)
          Use mpi
          
          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:), INTENT(IN) :: List ! new
          Integer, INTENT(IN)               :: Nsites, N_SUN,nf_eff ! new
          Complex (kind=kind(0.d0)), INTENT(OUT), Dimension(:,:) :: GreenA, GreenA_tmp, IDA

          Integer, Dimension(:), Allocatable :: PIVOT
          Complex (kind=kind(0.d0)) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, start_flav, dim_sq
          Integer         , Dimension(:), Allocatable :: SortedFlavors ! new

          Calc_Renyi_Ent_single=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
          
          dim=size(IDA,1)
          dim_sq=dim*dim
          ! We store the reduced Green's function in GreenA_c_tmp (c electrons)
          ! and GreenA_f_tmp (f electrons)
          ! The first Nsites columns contains the last flavour sector,
          ! the last Nsites column are old (Nf>2) or random, but they won't be used
          DO J = 1, Nsites
            DO I = 1, Nsites
                GreenA_tmp(I,J) = GRC(List(I), List(J), nf_eff)
            END DO
          END DO
          
          ! This exchange the last Nsites columns of GreenA_c_tmp between the two replicas
          ! such that GreenA contains the reduced Green's function for two replicas
          ! and a fixed spin sector.
          CALL MPI_ALLTOALL(GreenA_tmp, dim_sq, MPI_COMPLEX16, GreenA, dim_sq, MPI_COMPLEX16, ENTCOMM, IERR)
          
          DET = cmplx(1.D0,0.D0,KIND(0.D0))
          if(ENT_RANK==0) then
            dim_eff = NSites
            ! Compute Identity - GreenA(replica=1) - GreenA(replica=2) + 2 GreenA(replica=1) * GreenA(replica=2)
            IDA = - GreenA(1:dim_eff, 1:dim_eff) - GreenA(1:dim_eff, dim+1:dim+dim_eff)
            DO I = 1, dim_eff
                IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
            END DO
            CALL ZGEMM('n', 'n', dim_eff, dim_eff, dim_eff, alpha, GreenA(1, 1), &
                & dim, GreenA(1, dim+1), dim, beta, IDA, dim)
            ! Compute determinant
            SELECT CASE(dim_eff)
            CASE (1)
              DET = IDA(1,1)
            CASE (2)
              DET = IDA(1,1) * IDA(2,2) - IDA(1,2) * IDA(2,1)
            CASE DEFAULT
              Allocate(PIVOT(dim_eff))
              CALL ZGETRF(dim_eff, dim_eff, IDA, dim, PIVOT, INFO)
              DO I = 1, dim_eff
                  IF (PIVOT(I).NE.I) THEN
                    DET = -DET * IDA(I,I)
                  ELSE
                    DET = DET * IDA(I,I)
                  END IF
              ENDDO
              Deallocate(PIVOT)
            END SELECT
            Det=Det**N_SUN
          endif
          
          ! Compute the product of determinants for up and down spin sectors.
          CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
          ! Now each thread contains in PRODDET the full determinant, as obtained by
          ! a pair of replicas.
          Calc_Renyi_Ent_single = Calc_Renyi_Ent_single * PRODDET
        end function Calc_Renyi_Ent_single
#endif
        
      end Module entanglement_mod
      
