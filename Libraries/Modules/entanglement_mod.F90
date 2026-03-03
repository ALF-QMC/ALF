!  Copyright (C) 2020-2026 The ALF project
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
!> @author ALF-project
!> @brief Calculation of Renyi entanglement entropy and mutual information for quantum many-body systems.
!
!> @details
!> This module provides functionality for computing quantum entanglement measures
!> using the replica trick in Monte Carlo simulations. The key observables are:
!> - Renyi entanglement entropy: Measures quantum entanglement of subsystems
!> - Mutual information: Quantifies correlations between separated regions
!>
!> The module implements the second Renyi entropy calculation via the replica trick,
!> which requires pairing MPI replicas. The algorithm works by:
!> 1. Partitioning MPI tasks into pairs of replicas
!> 2. Computing reduced Green's functions on specified site patches
!> 3. Exchanging and combining Green's functions between replica pairs
!> 4. Calculating determinants to obtain the entanglement entropy
!>
!> The module supports three increasingly general patch specifications:
!> - Independent: Same patch for all flavors and colors
!> - Flavor-dependent: Different patches per flavor, same across colors
!> - Fully general: Independent patches for each flavor-color combination
!
!> @note
!> MPI is required for Renyi entropy calculations. Without MPI, only other
!> observables from the Green's function are computed. An even number of MPI
!> tasks is optimal; odd numbers result in one unpaired task.
!
!> @warning
!> The replica trick algorithm modifies the global MPI communicator structure
!> via Init_Entanglement_replicas(), which must be called before computing
!> any entanglement observables.
!
!> @see
!> For theoretical background: Phys. Rev. Lett. 104, 157201 (2010)
!--------------------------------------------------------------------
      ! Module-level MPI state variables for replica pairing
      INTEGER, save, private :: ENTCOMM      !< MPI communicator for replica pairs
      INTEGER, save, private :: ENT_RANK     !< Rank within ENTCOMM (0 or 1)
      INTEGER, save, private :: ENT_SIZE=0   !< Size of ENTCOMM (should be 2)
      INTEGER, save, private :: Norm         !< Normalization factor for averaging
      INTEGER, save, private :: group        !< Global MPI group communicator
      Real (kind=kind(0.d0)), save, private :: weight  !< Weight factor for replica averaging

      !> Generic interface for computing Renyi entanglement entropy
      !> Supports independent, flavor-dependent, and fully general patch specifications
      INTERFACE Calc_Renyi_Ent
        MODULE PROCEDURE Calc_Renyi_Ent_gen_all, Calc_Renyi_Ent_indep, Calc_Renyi_Ent_gen_fl
      END INTERFACE
      
      !> Generic interface for computing mutual information between regions
      !> Returns I(A:B) = S_A + S_B - S_AB where S denotes Renyi entropy
      INTERFACE Calc_Mutual_Inf
        MODULE PROCEDURE Calc_Mutual_Inf_indep, Calc_Mutual_Inf_gen_fl, Calc_Mutual_Inf_gen_all
      END INTERFACE
      Contains
!========================================================================

        Subroutine Init_Entanglement_replicas(Group_Comm)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Initializes MPI communicator structure for replica-based entanglement calculations.
!
!> @details
!> This subroutine partitions the MPI tasks into pairs of replicas for the
!> replica trick algorithm. Each pair forms a sub-communicator ENTCOMM with
!> 2 members (ranks 0 and 1). The pairing scheme is:
!> - Tasks (0,1), (2,3), (4,5), ..., (2n, 2n+1)
!>
!> The subroutine also computes normalization factors needed for averaging
!> entanglement observables across all replica pairs.
!
!> @param[in] Group_Comm MPI communicator containing all replicas (typically MPI_COMM_WORLD)
!
!> @note
!> Must be called before any Calc_Renyi_Ent or Calc_Mutual_Inf functions.
!> If the number of tasks is odd, one task remains unpaired and contributes
!> zero to the entanglement observable.
!
!> @warning
!> Only available when compiled with MPI support. Without MPI, this is a no-op.
!--------------------------------------------------------------------
#ifdef MPI
          Use mpi
#endif
          Implicit none
          Integer, INTENT(IN)               :: Group_Comm
          
          Integer                           :: ISIZE, IRANK, IERR
#ifdef MPI
          ! Query global communicator properties
          CALL MPI_COMM_SIZE(Group_Comm,ISIZE,IERR)
          CALL MPI_COMM_RANK(Group_Comm,IRANK,IERR)
          
          ! Create subgroups of two replicas each using integer division
          ! Tasks 2n and 2n+1 get the same color (n) and thus share ENTCOMM
          CALL MPI_COMM_SPLIT(Group_Comm, IRANK / 2, IRANK, ENTCOMM, IERR)
          CALL MPI_COMM_RANK(ENTCOMM,ENT_RANK,IERR)
          CALL MPI_COMM_SIZE(ENTCOMM,ENT_SIZE,IERR)
          
          ! Compute normalization factors for averaging over replica pairs
          Norm=ISIZE/2        ! Number of replica pairs
          norm=2*norm         ! Each task in a pair contributes independently
          group=Group_Comm    ! Store global communicator
          weight=dble(ISIZE)/dble(Norm)  ! Weight for proper averaging
#endif
          
        end Subroutine Init_Entanglement_replicas
          
!========================================================================
        Subroutine Calc_Mutual_Inf_indep(GRC,List_A,Nsites_A,List_B,Nsites_B,N_SUN,Renyi_A,Renyi_B,Renyi_AB)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Computes mutual information I(A:B) for two regions using independent patches.
!
!> @details
!> The mutual information is computed as I(A:B) = S_A + S_B - S_AB where
!> S denotes the exponentiated Renyi entropy. This measures the quantum
!> correlations between regions A and B.
!>
!> This version uses the same patch specification (List, Nsites, N_SUN) for
!> all flavor and color degrees of freedom.
!
!> @param[in] GRC Green's function array
!> @param[in] List_A Site indices for region A
!> @param[in] Nsites_A Number of sites in region A
!> @param[in] List_B Site indices for region B
!> @param[in] Nsites_B Number of sites in region B
!> @param[in] N_SUN Number of color degrees of freedom
!> @param[out] Renyi_A Exponentiated Renyi entropy of region A
!> @param[out] Renyi_B Exponentiated Renyi entropy of region B
!> @param[out] Renyi_AB Exponentiated Renyi entropy of combined region A∪B
!
!> @see Calc_Renyi_Ent_indep
!--------------------------------------------------------------------
          Implicit none
          
          Complex (kind=kind(0.d0)), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:), INTENT(IN) :: List_A, List_B
          Integer, INTENT(IN)               :: Nsites_A ,Nsites_B, N_SUN
          Complex (kind=kind(0.d0)), INTENT(OUT)   :: Renyi_A, Renyi_B, Renyi_AB

          Integer, Dimension(:), Allocatable :: List_AB
          Integer          :: I, J, IERR, INFO, Nsites_AB

          ! Combine regions A and B into a single list
          Nsites_AB=Nsites_A+Nsites_B
          allocate(List_AB(Nsites_AB))
          
          ! Copy region A sites
          DO I = 1, Nsites_A
             List_AB(I) = List_A(I)
          END DO
          
          ! Append region B sites
          DO I = 1, Nsites_B
             List_AB(I+Nsites_A) = List_B(I)
          END DO
          
          Renyi_A  = Calc_Renyi_Ent_indep(GRC,List_A,Nsites_A,N_SUN)
          Renyi_B  = Calc_Renyi_Ent_indep(GRC,List_B,Nsites_B,N_SUN)
          Renyi_AB = Calc_Renyi_Ent_indep(GRC,List_AB,Nsites_AB,N_SUN)
          
          deallocate(List_AB)
          
        End Subroutine Calc_Mutual_Inf_indep
          
!========================================================================
        Subroutine Calc_Mutual_Inf_gen_fl(GRC,List_A,Nsites_A,List_B,Nsites_B,N_SUN,Renyi_A,Renyi_B,Renyi_AB)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Computes mutual information with flavor-dependent patches.
!
!> @details
!> Similar to Calc_Mutual_Inf_indep, but allows different patch specifications
!> for each flavor. Within each flavor, all color degrees of freedom use the
!> same patch.
!
!> @param[in] GRC Green's function array
!> @param[in] List_A Site lists for region A, List_A(:,f) for flavor f
!> @param[in] Nsites_A Number of sites per flavor for region A
!> @param[in] List_B Site lists for region B, List_B(:,f) for flavor f
!> @param[in] Nsites_B Number of sites per flavor for region B
!> @param[in] N_SUN Number of colors per flavor
!> @param[out] Renyi_A Exponentiated Renyi entropy of region A
!> @param[out] Renyi_B Exponentiated Renyi entropy of region B
!> @param[out] Renyi_AB Exponentiated Renyi entropy of combined region A∪B
!
!> @see Calc_Renyi_Ent_gen_fl
!--------------------------------------------------------------------
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
        Subroutine Calc_Mutual_Inf_gen_all(GRC,List_A,Nsites_A,List_B,Nsites_B,Renyi_A,Renyi_B,Renyi_AB)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Computes mutual information with fully general patch specifications.
!
!> @details
!> Most general version allowing independent patch specifications for every
!> combination of flavor and color degrees of freedom. List_A(:,f,c) specifies
!> the patch for flavor f and color c.
!
!> @param[in] GRC Green's function array
!> @param[in] List_A Site lists for region A, List_A(:,f,c) for flavor f, color c
!> @param[in] Nsites_A Number of sites, Nsites_A(f,c) for flavor f, color c
!> @param[in] List_B Site lists for region B, List_B(:,f,c) for flavor f, color c
!> @param[in] Nsites_B Number of sites, Nsites_B(f,c) for flavor f, color c
!> @param[out] Renyi_A Exponentiated Renyi entropy of region A
!> @param[out] Renyi_B Exponentiated Renyi entropy of region B
!> @param[out] Renyi_AB Exponentiated Renyi entropy of combined region A∪B
!
!> @see Calc_Renyi_Ent_gen_all
!--------------------------------------------------------------------
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
!> Computes second Renyi entropy for a single patch specification across all flavors.
!>
!> @details
!> Returns exp(-S_2) where S_2 is the second Renyi entropy. The same patch
!> (List, Nsites) is used for all flavor and color degrees of freedom.
!>
!> The algorithm:
!> 1. Pairs replicas via MPI communicators (set up by Init_Entanglement_replicas)
!> 2. Extracts reduced Green's functions for the specified sites
!> 3. Exchanges Green's functions between replica pairs
!> 4. Computes determinants of (1-G₁-G₂+2G₁G₂) matrices
!> 5. Multiplies determinants across flavors and averages over replica pairs
!>
!> @param[in] GRC Green's function, GRC(i,j,f) for sites i,j and flavor f
!> @param[in] List Site indices defining the entanglement patch
!> @param[in] Nsites Number of sites in the patch
!> @param[in] N_SUN Number of color degrees of freedom
!
!> @return Complex value of exp(-S_2), or 0 if MPI not available or odd number of replicas
!
!> @note
!> Requires ENT_SIZE=2 (replica pairs). Returns 0 for unpaired replicas.
!> Without MPI compilation, always returns 0 and prints a warning on first call.
!
!> @warning
!> Init_Entanglement_replicas() must be called before using this function.
!
!> @see Init_Entanglement_replicas, Calc_Renyi_Ent_pair, Calc_Renyi_Ent_single
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
!> Computes Renyi entropy with flavor-dependent patches.
!>
!> @details
!> Extension of Calc_Renyi_Ent_indep allowing different patches for each flavor.
!> Within a given flavor f, all N_SUN(f) color degrees of freedom use the same
!> patch List(:,f) with Nsites(f) sites.
!>
!> The algorithm sorts flavors by patch size for computational efficiency,
!> processing larger patches first to optimize memory allocation.
!>
!> @param[in] GRC Green's function array
!> @param[in] List Site lists, List(:,f) for flavor f
!> @param[in] Nsites Number of sites per flavor, Nsites(f) for flavor f
!> @param[in] N_SUN Number of color degrees of freedom per flavor
!
!> @return Complex value of exp(-S_2), averaged over replica pairs
!
!> @note
!> Flavors with Nsites(f)=0 are skipped. The function uses insertion sort
!> to order patches by size (efficient for small numbers of flavors).
!
!> @see Calc_Renyi_Ent_indep, Calc_Renyi_Ent_gen_all
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
!> Computes Renyi entropy with fully general patch specifications.
!>
!> @details
!> Most general version allowing independent patch List(:,f,c) and size Nsites(f,c)
!> for every combination of flavor f and color c.
!>
!> The algorithm treats each (flavor, color) combination as an independent degree
!> of freedom and sorts all N_FL×num_nc patches by size. This ensures optimal
!> memory usage during the calculation.
!>
!> @param[in] GRC Green's function array
!> @param[in] List Site lists, List(:,f,c) for flavor f and color c
!> @param[in] Nsites Number of sites, Nsites(f,c) for flavor f and color c
!
!> @return Complex value of exp(-S_2), averaged over replica pairs
!
!> @note
!> This is the most flexible but also most memory-intensive version.
!> Use flavor-dependent or independent versions when possible for efficiency.
!
!> @see Calc_Renyi_Ent_indep, Calc_Renyi_Ent_gen_fl
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
!> Helper function computing Renyi entropy contribution from a pair of patches.
!>
!> @details
!> This internal helper function processes a pair of patches (one per replica)
!> that may differ in size, flavor, and color count. It:
!> 1. Extracts reduced Green's functions for both patches
!> 2. Exchanges Green's functions between replica pairs via MPI_ALLTOALL
!> 3. Computes matrix M = 1 - G₁ - G₂ + 2G₁G₂
!> 4. Calculates det(M)^N_SUN for appropriate N_SUN
!> 5. Multiplies determinants from both replicas
!>
!> @param[in] GRC Green's function array
!> @param[in] List Site lists, List(:,1:2) for replica 1 and 2
!> @param[in] Nsites Number of sites, Nsites(1:2) for each patch
!> @param[in] nf_list Flavor indices, nf_list(1:2) for each patch
!> @param[in] N_SUN Number of colors, N_SUN(1:2) for each patch
!> @param[out] GreenA Work array (dim,2*dim) for Green's function storage
!> @param[out] GreenA_tmp Work array (dim,2*dim) for MPI exchange buffer
!> @param[out] IDA Work array (dim,dim) for matrix operations
!
!> @return Product of determinants from both replicas: det₁^N_SUN(1) × det₂^N_SUN(2)
!
!> @note
!> This is an internal helper function. Use the public Calc_Renyi_Ent_* interfaces.
!> Requires dim ≥ max(Nsites(1), Nsites(2)).
!
!> @warning
!> Only available with MPI. Work arrays must be pre-allocated by caller.
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
          
          ! Extract reduced Green's function for patch 1 (first replica)
          ! Store in first Nsites(1) columns of GreenA_tmp
          nf_eff = nf_list(1)

          DO J = 1, Nsites(1)
            Do I = 1, Nsites(1)
              GreenA_tmp(I,J) = GRC(List(I,1),List(J,1),nf_eff)
            END DO
          END Do

          ! Extract reduced Green's function for patch 2 (second replica)
          ! Store in columns dim+1 to dim+Nsites(2) of GreenA_tmp
          nf_eff = nf_list(2)
          DO J = 1, Nsites(2)
            DO I = 1, Nsites(2)
                GreenA_tmp(I,J+dim) = GRC(List(I,2), List(J,2), nf_eff)
            END DO
          END DO
          
          ! Exchange Green's functions between replica pairs via MPI
          ! After this, GreenA contains cross-replica Green's functions:
          ! Columns 1:dim from partner, columns dim+1:2*dim from own replica
          CALL MPI_ALLTOALL(GreenA_tmp, dim_sq, MPI_COMPLEX16, GreenA, dim_sq, MPI_COMPLEX16, ENTCOMM, IERR)

          ! Compute the matrix: M = I - G₁ - G₂ + 2G₁G₂
          ! This is the core of the replica trick calculation
          dim_eff = Nsites(1+ENT_RANK)
          ! Start with -G₁ - G₂
          IDA(1:dim_eff,1:dim_eff) = - GreenA(1:dim_eff, 1:dim_eff) - GreenA(1:dim_eff, dim+1:dim+dim_eff)
          ! Add identity matrix
          DO I = 1, dim_eff
              IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
          END DO
          ! Add 2G₁G₂ term using matrix multiplication
          CALL ZGEMM('n', 'n', dim_eff, dim_eff, dim_eff, alpha, GreenA(1, 1), &
              & dim, GreenA(1, dim+1), dim, beta, IDA, dim)
              
          ! Compute determinant of M using different algorithms based on size
          SELECT CASE(dim_eff)
          CASE (1)
            ! Trivial case: det = M₁₁
            DET = IDA(1,1)
          CASE (2)
            ! Analytical formula for 2×2 matrix
            DET = IDA(1,1) * IDA(2,2) - IDA(1,2) * IDA(2,1)
          CASE DEFAULT
            ! Use LU decomposition for larger matrices
            Allocate(PIVOT(dim_eff))
            CALL ZGETRF(dim_eff, dim_eff, IDA, dim, PIVOT, INFO)
            DET = cmplx(1.D0,0.D0,KIND(0.D0))
            ! Compute determinant from LU factors, accounting for pivoting
            DO I = 1, dim_eff
                IF (PIVOT(I).NE.I) THEN
                  DET = -DET * IDA(I,I)  ! Pivoting changes sign
                ELSE
                  DET = DET * IDA(I,I)
                END IF
            ENDDO
            Deallocate(PIVOT)
          END SELECT
          
          ! Raise to power N_SUN to account for color degrees of freedom
          DET=DET**N_SUN(1+ENT_RANK)
          
          ! Multiply determinants from both replicas in the pair
          CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
          
          ! Each replica now has the complete product det₁ × det₂
          Calc_Renyi_Ent_pair = Calc_Renyi_Ent_pair * PRODDET
        end function Calc_Renyi_Ent_pair

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Helper function for handling unpaired replica when total number of replicas is odd.
!>
!> @details
!> When the number of MPI tasks is odd, one replica remains unpaired. This function
!> handles the Renyi entropy calculation for the unpaired patch. Only the replica
!> with ENT_RANK=0 performs the actual calculation; the other replica contributes
!> det=1.
!>
!> The algorithm follows the same steps as Calc_Renyi_Ent_pair but with only one
!> active patch.
!>
!> @param[in] GRC Green's function array
!> @param[in] List Site list for the patch
!> @param[in] Nsites Number of sites in the patch
!> @param[in] nf_eff Flavor index of the patch
!> @param[in] N_SUN Number of color degrees of freedom
!> @param[out] GreenA Work array (dim,2*dim) for Green's function storage
!> @param[out] GreenA_tmp Work array (dim,2*dim) for MPI exchange buffer
!> @param[out] IDA Work array (dim,dim) for matrix operations
!
!> @return det^N_SUN from the unpaired patch (or 1 from the paired partner)
!
!> @note
!> This is an internal helper function called automatically when N_FL is odd.
!> The contribution is averaged with weight accounting for the unpaired status.
!
!> @warning
!> Only available with MPI. Only rank 0 in ENTCOMM performs the calculation.
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
          
          ! Extract reduced Green's function for the unpaired patch
          ! Store in first Nsites columns
          DO J = 1, Nsites
            DO I = 1, Nsites
                GreenA_tmp(I,J) = GRC(List(I), List(J), nf_eff)
            END DO
          END DO
          
          ! Exchange with partner (even though partner won't use it)
          ! This maintains MPI communication symmetry
          CALL MPI_ALLTOALL(GreenA_tmp, dim_sq, MPI_COMPLEX16, GreenA, dim_sq, MPI_COMPLEX16, ENTCOMM, IERR)
          
          DET = cmplx(1.D0,0.D0,KIND(0.D0))
          
          ! Only rank 0 in ENTCOMM performs the actual calculation
          ! The partner (rank 1) contributes DET=1
          if(ENT_RANK==0) then
            dim_eff = NSites
            ! Compute M = I - G₁ - G₂ + 2G₁G₂ (same as paired case)
            IDA = - GreenA(1:dim_eff, 1:dim_eff) - GreenA(1:dim_eff, dim+1:dim+dim_eff)
            DO I = 1, dim_eff
                IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
            END DO
            CALL ZGEMM('n', 'n', dim_eff, dim_eff, dim_eff, alpha, GreenA(1, 1), &
                & dim, GreenA(1, dim+1), dim, beta, IDA, dim)
                
            ! Compute determinant using appropriate method
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
            
            ! Raise to power N_SUN for color degrees of freedom
            Det=Det**N_SUN
          endif
          
          ! Broadcast result from rank 0 to rank 1 (rank 1 contributes DET=1)
          CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
          
          ! Both replicas now have the result
          Calc_Renyi_Ent_single = Calc_Renyi_Ent_single * PRODDET
        end function Calc_Renyi_Ent_single
#endif
        
      end Module entanglement_mod
      
