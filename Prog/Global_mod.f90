!  Copyright (C) 2016, 2017 The ALF project
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



Module Global_mod

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles global updates.
!> Handles parrallel tempering
!
!--------------------------------------------------------------------

      Use Hamiltonian
      Use MyMats 
      Use Operator_mod
      Use Control
      Use Observables
      
      Implicit none

      Type (Obser_Vec ),  private  ::   Tempering_acceptance

      
    Contains
#if defined(TEMPERING)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles parrallel tempering
!> This subroutine is called only if the tempering flag is switched on. In this 
!> case the MPI flag is also switched on. 
!> 
!--------------------------------------------------------------------
      Subroutine Exchange_Step(Phase,GR, udvr, udvl, Stab_nt, udvst, N_exchange_steps, Tempering_calc_det)
        Use UDV_State_mod
        Implicit none
        
        include 'mpif.h'
        
        Interface
           SUBROUTINE WRAPUL(NTAU1, NTAU, udvl)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout) :: udvl(N_FL)
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: udvl, udvr
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface
        
        !>  Arguments
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT)                   :: Phase
        CLASS(UDV_State), intent(inout), allocatable, Dimension(:) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
        CLASS(UDV_State), intent(inout), allocatable, Dimension(:, :) :: udvst
        INTEGER, dimension(:),     INTENT   (IN), allocatable      :: Stab_nt
        !>  On entry and on exit the left storage is full, and the Green function is on time slice 0 and the phase is set.
        
        
        !>  Local variables.
        Integer :: NST, NSTM, NF, NT, NT1, NVAR,N, N1,N2, I, NC, I_Partner, n_step, N_exchange_steps, N_count
        Integer, Dimension(:,:),  allocatable :: nsigma_old
        Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, Weight, Weight1
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Ratiotot, Ratiotot_p, Phase_old, Phase_new
        Complex (Kind=Kind(0.d0)), allocatable :: Det_vec_old(:,:), Det_vec_new(:,:), Phase_Det_new(:), Phase_Det_old(:)
        Complex (Kind=Kind(0.d0)) :: Ratio(2), Ratio_p(2)
        Logical :: TOGGLE, L_Test
        
        Integer, allocatable :: List_partner(:), List_masters(:)
        
        !> Additional variables for running without Fermion weight
        Logical :: Tempering_calc_det
        Integer        :: nsigma_irank, nsigma_old_irank, nsigma_irank_temp ! Keeps track of where the configuration originally comes from
        Integer        :: n_GR
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), allocatable :: GR_new
        !Integer, Dimension(:,:),  allocatable :: nsigma_orig, nsigma_test
        !Integer :: I1, I2
        
        Integer        :: Isize, Irank, Ierr, irank_g, isize_g, igroup
        Integer        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
        nsigma_irank = irank
        
        
        
        n1 = size(nsigma,1)
        n2 = size(nsigma,2)
        NSTM = Size(udvst, 1)
        Allocate ( nsigma_old(n1,n2) )
        if (Tempering_calc_det) then
           Allocate ( Det_vec_old(NDIM,N_FL), Det_vec_new(NDIM,N_FL) ) 
           Allocate ( Phase_Det_new(N_FL), Phase_Det_old(N_FL) )
        endif
        Allocate ( List_partner(0:Isize-1), List_masters(Isize/2) )
        
        !         Allocate ( nsigma_orig(n1,n2) )
        !         nsigma_orig = nsigma
        
        !>  Compute for each core the old weights.     
        L_test = .false.
        if (Tempering_calc_det) then
           ! Set old weight. 
           Phase_old =cmplx(1.d0,0.d0,kind(0.d0))
           do nf = 1,N_Fl
              Call Compute_Fermion_Det(Z,Det_Vec_old(:,nf), udvl(nf))
              Phase_det_old(nf) = Z
              Phase_old = Phase_old*Z
           Enddo
           call Op_phase(Phase_old,OP_V,Nsigma,N_SUN) 
        endif
        !> Store old configuration
        nsigma_old = nsigma
        nsigma_old_irank = nsigma_irank
        ! Setup the list of masters
        nc = 0
        Do I  = 0,Isize-1,2*Isize_g
           do n = 0,isize_g-1
              nc = nc + 1
              List_masters(nc)  = I + n
           enddo
        Enddo
        if (L_Test) then
           if (Irank == 0)  then
              Write(6,*) 'List  of masters: '
              Do nc = 1,Isize/2
                 Write(6,*) nc, list_masters(nc)
              Enddo
           endif
        endif
        DO N_count = 1,N_exchange_steps
           
           !>  Set the partner rank on each core
           If (Irank == 0 ) then
              n_step = isize_g
              if (  ranf_wrap() > 0.5d0 ) n_step = -isize_g
              Do I = 0,Isize-1,2*isize_g
                 do n = 0,isize_g-1
                    List_partner(npbc_tempering(I  + n             ,Isize)) =  npbc_tempering(I +  n   + n_step,Isize)
                    List_partner(npbc_tempering(I  + n  + n_step  , Isize)) =  npbc_tempering(I + n            ,Isize)
                 enddo
              enddo
           endif

           CALL MPI_BCAST(List_partner, Isize  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
           
           If (L_test) then
              if (Irank == 0 ) then
                 Write(6,*) 'Testing global : '
                 do  I = 0,Isize -1 
                    Write(6,*)  I, List_partner(I) !, Phase_old, Phase
                    Write(11,*)  I, List_partner(I) !, Phase_old, Phase
                 enddo
                 Write(11,*)
              endif
           Endif
           !call MPI_Barrier(MPI_COMM_WORLD)
           !Write(6,*) '---------'
!!$    select case (IRANK)
!!$    case(0)
!!$       nsigma_old(1,1) =  1;  nsigma_old(2,1) = 1
!!$    case(1)
!!$       nsigma_old(1,1) = -1;  nsigma_old(2,1) = 1
!!$    case(2)
!!$       nsigma_old(1,1) =  1;  nsigma_old(2,1) = -1
!!$    case(3)
!!$       nsigma_old(1,1) = -1;  nsigma_old(2,1) = -1
!!$    case default
!!$    end select
    

           !>  Exchange configurations
           n = size(nsigma_old,1)*size(nsigma_old,2)
           CALL MPI_Sendrecv(nsigma_old      , n, MPI_INTEGER, List_partner(IRANK), 0, &
                    &        nsigma          , n, MPI_INTEGER, List_partner(IRANK), 0, MPI_COMM_WORLD,STATUS,IERR)
           CALL MPI_Sendrecv(nsigma_old_irank, 1, MPI_INTEGER, List_partner(IRANK), 0, &
                    &        nsigma_irank    , 1, MPI_INTEGER, List_partner(IRANK), 0, MPI_COMM_WORLD,STATUS,IERR)
           
           !>  Each node now has a new configuration nsigma
           
!!$    If (L_test) then
!!$       Write(6,*) 'Testing global : ', Irank,List_partner(IRANK), nsigma_old(1,1),  nsigma_old(2,1), nsigma(1,1),  nsigma(2,1) 
!!$    Endif
           
           
           
    if (Tempering_calc_det) then
           !>  Compute ratio on weights one each rank
           DO nf = 1,N_FL
              CALL udvl(nf)%reset
           ENDDO
           DO NST = NSTM-1,1,-1
              NT1 = Stab_nt(NST+1)
              NT  = Stab_nt(NST  )
              !Write(6,*) NT1,NT, NST
              CALL WRAPUL(NT1,NT, udvl)
              Do nf = 1,N_FL
                 udvst(NST, nf) = udvl(nf)
              ENDDO
           ENDDO
           NT1 = stab_nt(1)
           CALL WRAPUL(NT1,0, udvl)
           Phase_new = cmplx(1.d0,0.d0,kind(0.d0))
           do nf = 1,N_Fl
              Call Compute_Fermion_Det(Z,Det_Vec_new(:,nf), udvl(nf))
              Phase_det_new(nf) = Z
              Phase_new = Phase_new*Z
           Enddo
           call Op_phase(Phase_new,OP_V,Nsigma,N_SUN) 
           
           T0_Proposal_ratio = 1.d0
           Ratiotot = Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
                &                               Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio,Ratio) 
           
           If (L_Test) Write(6,*) 'Ratio_global: Irank, Partner',Irank,List_partner(Irank), Ratiotot, Ratio(1)*exp(Ratio(2))
        else
           Ratiotot = Delta_S0_global(Nsigma_old)
           Ratio(1) = Ratiotot
           Ratio(2) = 0
        endif
    
           !>  Acceptace/rejection decision taken on master node after receiving information from slave
           Do nc = 1,Isize/2 ! Loop over masters
              I = List_masters(nc)
              If (Irank == I ) Then
                 !CALL MPI_SEND(Ratiotot,1, MPI_COMPLEX16  , List_partner(I), I+512 , MPI_COMM_WORLD,IERR)
                 CALL MPI_SEND(Ratio   ,2, MPI_COMPLEX16  , List_partner(I), I+1024, MPI_COMM_WORLD,IERR)
              else if (IRANK == List_Partner(I) ) Then
                 !CALL MPI_RECV(Ratiotot_p , 1, MPI_COMPLEX16,  I, I+512 , MPI_COMM_WORLD,STATUS,IERR)
                 CALL MPI_RECV(Ratio_p    , 2, MPI_COMPLEX16,  I, I+1024, MPI_COMM_WORLD,STATUS,IERR)
                 !Weight = abs(Ratiotot_p*Ratiotot)
                 Weight= abs(Ratio(1) * Ratio_p(1) * exp( Ratio_p(2) + Ratio(2)  ) )
                 TOGGLE = .false. 
                 if ( Weight > ranf_wrap() )  Toggle =.true.
                 If (L_Test) Write(6,*) 'Master : ', List_Partner(I), I, Weight, Toggle
              endif
           enddo
           
           !>  Send result of acceptance/rejection decision form master to slave
           Do nc =  1,Isize/2 ! Loop over masters
              I = List_masters(nc)
              If (Irank == List_Partner(I) ) Then
                 ! Write(6,*) 'Send from ', List_Partner(I), 'to, ', I, I + 512
                 CALL MPI_SEND(Toggle, 1, MPI_LOGICAL, I, I+512, MPI_COMM_WORLD,IERR)
              else if (IRANK == I ) Then
                 CALL MPI_RECV(Toggle , 1, MPI_LOGICAL,   List_partner(I), I+512 ,MPI_COMM_WORLD,STATUS,IERR)
                 If (L_Test) Write(6,*) 'Slave : ', Irank,  Toggle
              endif
           enddo
           
           If (L_Test) Write(6,*) 'Final: ',  Irank, List_partner(Irank), toggle
           Call Global_Tempering_obser(Toggle)
           Call Control_upgrade_Temp  (Toggle) 
           If (toggle)  then
              !>     Move has been accepted
    if (Tempering_calc_det) then
              Phase_old     = Phase_new
              Phase_det_old = Phase_det_new
              Det_vec_old   = Det_vec_new
    endif
              nsigma_old       = nsigma
              nsigma_old_irank = nsigma_irank
           else
              nsigma       = nsigma_old
              nsigma_irank = nsigma_old_irank
           endif
        enddo
        
        !> Finalize
    if (Tempering_calc_det) then
        !> If move has been accepted, no use to recomute storage
        If (.not.TOGGLE) then
           DO nf = 1,N_FL
              CALL udvl(nf)%reset
           ENDDO
           DO NST = NSTM-1,1,-1
              NT1 = Stab_nt(NST+1)
              NT  = Stab_nt(NST  )
              !Write(6,*) NT1,NT, NST
              CALL WRAPUL(NT1,NT, udvl)
              Do nf = 1,N_FL
                 udvst(NST, nf) = udvl(nf)
              ENDDO
           ENDDO
           NT1 = stab_nt(1)
           CALL WRAPUL(NT1,0, udvl)
        Endif
        !> Compute the Green functions so as to provide correct starting point for the sequential updates.
        NVAR  = 1
        Phase = cmplx(1.d0,0.d0,kind(0.d0))
        do nf = 1,N_Fl
           CALL CGR(Z, NVAR, GR(:,:,nf), udvr(nf),  udvl(nf))
           Phase = Phase*Z
        Enddo
        call Op_phase(Phase,OP_V,Nsigma,N_SUN)     
    else
        !> Send >>Phase, GR, udvr, udvl, udvst<< to new node 
        !  First step: Each node sends to IRANK=0 its value nsigma_irank, which is the node where its new Phase, GR, udvr, udvl, udvst is stored
        !              This node then tells each node where to send its now old Phase, GR, udvr, udvl, udvst
        !              Finally, the variables get submitted
        If (Irank == 0) then
           Do I = 1,Isize-1
              CALL MPI_RECV(nsigma_irank_temp , 1, MPI_INTEGER, I, 0, MPI_COMM_WORLD,STATUS,IERR)
              If ( nsigma_irank_temp == 0) then
                 nsigma_old_irank = I
              else
                 CALL MPI_SEND(I , 1, MPI_INTEGER, nsigma_irank_temp, 0, MPI_COMM_WORLD,IERR)
              endif
           enddo
           If ( nsigma_irank /= 0 ) then
              CALL MPI_SEND(0 , 1, MPI_INTEGER, nsigma_irank, 0, MPI_COMM_WORLD,IERR)
           endif
        else
           CALL MPI_SEND(nsigma_irank     , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD,IERR)
           CALL MPI_RECV(nsigma_old_irank , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD,STATUS,IERR)
        endif
        
        if ( nsigma_irank /= irank ) then
           CALL MPI_Sendrecv(Phase,     1, MPI_COMPLEX16, nsigma_old_irank, 0, &
                    &        Phase_new, 1, MPI_COMPLEX16, nsigma_irank    , 0, MPI_COMM_WORLD,STATUS,IERR)
           Phase = Phase_new
           
           n_GR = size(GR,1)*size(GR,2)*size(GR,3)
           Allocate ( GR_new(size(GR,1),size(GR,2),size(GR,3)) )
           CALL MPI_Sendrecv(GR,     n_GR, MPI_COMPLEX16, nsigma_old_irank, 0, &
                    &        GR_new, n_GR, MPI_COMPLEX16, nsigma_irank    , 0, MPI_COMM_WORLD,STATUS,IERR)
           GR = GR_new
           Deallocate ( GR_new )
           
           do nf = 1,N_Fl
              CALL udvr(nf)%MPI_Sendrecv(nsigma_old_irank, 0, nsigma_irank, 0, STATUS, IERR)
           enddo
           do nf = 1,N_Fl
              CALL udvl(nf)%MPI_Sendrecv(nsigma_old_irank, 0, nsigma_irank, 0, STATUS, IERR)
           enddo
           do NST = 1, NSTM
              do nf = 1,N_Fl
                 CALL udvst(NST, nf)%MPI_Sendrecv(nsigma_old_irank, 0, nsigma_irank, 0, STATUS, IERR)
              enddo
           enddo
        endif   
    endif
        
        Deallocate ( nsigma_old )
    if (Tempering_calc_det) then
        Deallocate ( Det_vec_old, Det_vec_new ) 
        Deallocate ( Phase_Det_new, Phase_Det_old )
    endif
        Deallocate ( List_partner, List_masters )
        
      end Subroutine Exchange_Step
#endif
!---------------------------------------------------------------------
      Subroutine Global_Updates(Phase,GR, udvr, udvl, Stab_nt, udvst,N_Global)
        Use UDV_State_mod
        Implicit none
        
        Interface
           SUBROUTINE WRAPUL(NTAU1, NTAU, udvl)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout) :: udvl(N_FL)
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: udvl, udvr
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface
        
        !>  Arguments
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT)                   :: Phase
        CLASS   (UDV_State), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
        INTEGER, dimension(:),     INTENT   (IN), allocatable      :: Stab_nt
        Integer, INTENT(IN) :: N_Global
        !>  On entry and on exit the left storage is full, and the Green function is on time slice 0 and the phase is set.
        
        
        !>  Local variables.
        Integer :: NST, NSTM, NF, NT, NT1, NVAR,N, N1,N2, I, NC
        Integer, Dimension(:,:),  allocatable :: nsigma_old
        Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, Weight
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Ratiotot, Phase_old, Phase_new
        Complex (Kind=Kind(0.d0)), allocatable :: Det_vec_old(:,:), Det_vec_new(:,:), Phase_Det_new(:), Phase_Det_old(:)
        Complex (Kind=Kind(0.d0)) :: Ratio(2)
        Logical :: TOGGLE, L_Test
        Real    (Kind=Kind(0.d0)) :: size_clust
        
        
        
        !>  On entry and on exit the left storage is full, and the Green function is on time slice 0 and the phase is set.
        
        n1 = size(nsigma,1)
        n2 = size(nsigma,2)
        NSTM = Size(udvst, 1)
        Allocate ( nsigma_old(n1,n2) )
        Allocate ( Det_vec_old(NDIM,N_FL), Det_vec_new(NDIM,N_FL) ) 
        Allocate ( Phase_Det_new(N_FL), Phase_Det_old(N_FL) )
        
        
        L_test = .false.
        ! Write(6,*)
        ! Set old weight. 
        Phase_old =cmplx(1.d0,0.d0,kind(0.d0))
        do nf = 1,N_Fl
           Call Compute_Fermion_Det(Z,Det_Vec_old(:,nf), udvl(nf))
           Phase_det_old(nf) = Z
           Phase_old = Phase_old*Z
        Enddo
        call Op_phase(Phase_old,OP_V,Nsigma,N_SUN) 
        If (L_test) then
           Write(6,*) 'Testing global : ',  Phase_old, Phase
        Endif
        
        
        If (L_test) then 
           ! Testing    
           Do nf = 1,N_FL
              CALL udvr(nf)%reset
           Enddo
           NVAR = 1
           Phase = cmplx(1.d0,0.d0,kind(0.d0))
           do nf = 1,N_Fl
              CALL CGR(Z, NVAR, GR(:,:,nf), udvr(nf), udvl(nf))
              Phase = Phase*Z
           Enddo
           call Op_phase(Phase,OP_V,Nsigma,N_SUN) 
           Do Nf = 1,N_FL
              Call DET_C_LU(GR(:,:,nf),Det_vec_new(:,nf),Ndim)
              Z = Phase_det_old(nf)
              DO I = 1,Ndim
                 Z = Z*Det_vec_new(I,nf)*Det_vec_old(I,nf)
              Enddo
              Write(6,*) 'Testing weight: ', Z
           Enddo
        Endif
        
        !> Store old configuration
        nsigma_old = nsigma 
        !> Phase_old, Phase_det_old and Det_vec_old  are all set. 
        NC = 0
        Do n = 1,N_Global
           !> Draw a new spin configuration. This is provided by the user in the Hamiltonian module
           !> Note that nsigma is a variable in the module Hamiltonian
           Call Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
           If (T0_Proposal_ratio > 1.D-24) then
              NC = NC + 1
              !> Compute the new Green function
              DO nf = 1,N_FL
                 CALL udvl(nf)%reset
              ENDDO
              DO NST = NSTM-1,1,-1
                 NT1 = Stab_nt(NST+1)
                 NT  = Stab_nt(NST  )
                 !Write(6,*) NT1,NT, NST
                 CALL WRAPUL(NT1,NT,udvl)
                 Do nf = 1,N_FL
                    udvst(NST, nf) = udvl(nf)
                 ENDDO
              ENDDO
              NT1 = stab_nt(1)
              CALL WRAPUL(NT1,0, udvl)
              !You could now compute the det directly here.
              Phase_new = cmplx(1.d0,0.d0,kind(0.d0))
              do nf = 1,N_Fl
                 Call Compute_Fermion_Det(Z,Det_Vec_new(:,nf),udvl(nf))
                 Phase_det_new(nf) = Z
                 Phase_new = Phase_new*Z
              Enddo
              call Op_phase(Phase_new,OP_V,Nsigma,N_SUN) 
              
              Ratiotot = Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
                   &                               Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio, Ratio) 
              
              !Write(6,*) 'Ratio_global: ', Ratiotot
              
              Weight = abs(  real(Phase_old * Ratiotot, kind=Kind(0.d0))/real(Phase_old,kind=Kind(0.d0)) )
              
              Z = Phase_old * Ratiotot/ABS(Ratiotot)
              Call Control_PrecisionP_Glob(Z,Phase_new)
              !Write(6,*) Z, Phase_new
              
              
              TOGGLE = .false. 
              if ( Weight > ranf_wrap() )  Then
                 TOGGLE = .true.
                 Phase_old     = Phase_new
                 Phase_det_old = Phase_det_new
                 nsigma_old    = nsigma
                 Det_vec_old   = Det_vec_new
              else
                 nsigma = nsigma_old
              endif
              Call Control_upgrade_Glob(TOGGLE,size_clust)
           endif
        Enddo
        
        If (NC > 0 ) then
           If (.not.TOGGLE) then
              DO nf = 1,N_FL
                 CALL udvl(nf)%reset
              ENDDO
              DO NST = NSTM-1,1,-1
                 NT1 = Stab_nt(NST+1)
                 NT  = Stab_nt(NST  )
                 !Write(6,*) NT1,NT, NST
                 CALL WRAPUL(NT1,NT, udvl)
                 Do nf = 1,N_FL
                    udvst(NST, nf) = udvl(nf)
                 ENDDO
              ENDDO
              NT1 = stab_nt(1)
              CALL WRAPUL(NT1,0, udvl)
           Endif
           !Compute the Green functions so as to provide correct starting point for the sequential updates.
           NVAR  = 1
           Phase = cmplx(1.d0,0.d0,kind(0.d0))
           do nf = 1,N_Fl
              CALL CGR(Z, NVAR, GR(:,:,nf), udvr(nf), udvl(nf))
              Phase = Phase*Z
           Enddo
           call Op_phase(Phase,OP_V,Nsigma,N_SUN) 
        endif
        
        
        Deallocate ( nsigma_old)
        Deallocate ( Det_vec_old  , Det_vec_new  ) 
        Deallocate ( Phase_Det_new, Phase_Det_old )
        
        
      End Subroutine Global_Updates
      
      
      
!--------------------------------------------------------------------
      Complex (Kind=Kind(0.d0)) Function Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
           &                    Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio,Ratio)
!--------------------------------------------------------------------
!> @author
!> Fakher Assaad 
!>
!> @brief 
!> This fucntion computes ratio of weights  T0(nsigma--> nsigma_old) W(nsigma)/ 
!>                                          T0(nsigma_old--> nsigma) W(nsigma_old) =
!>                                          T0_Proposal_ratio   W(nsigma)/  W(nsigma_old)
!> 
!> Note that the new configuration, nsigma, is contained in the Hamiltonian moddule
!> The fermionic determinant stems from the routine Compute_Fermion_Det. 
!> Since the ratio can be a very large number, it is encoded as Ratio(1)*exp(Ratio(2))
!--------------------------------------------------------------------

    
        Implicit none
        
        !> Arguments
        Complex (Kind=Kind(0.d0)), allocatable, INTENT(IN) :: Phase_Det_old(:), Phase_Det_new(:), &
             &                                                Det_vec_old(:,:), Det_vec_new(:,:)
        Real    (Kind=Kind(0.d0)) :: T0_proposal_ratio 
        Integer, allocatable      :: nsigma_old(:,:)
        Complex (Kind=Kind(0.d0)), INTENT(out) :: Ratio(2)
        
        !> Local 
        Integer                                :: Nf, i, nt
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Real    (Kind=Kind(0.d0)) :: X, Ratio_2
        
        Ratio = cmplx(0.d0,0.d0,kind(0.d0))
        Ratio_2 = 0.d0
        !X = 1.d0
        Do nf = 1,N_Fl
           DO I = 1,Ndim
              !X= X*real(Det_vec_new(I,nf),kind(0.d0)) / Real(Det_vec_old(I,nf),kind(0.d0) )
              Ratio_2 = Ratio_2 +  log(real(Det_vec_new(I,nf),kind(0.d0))) - log( Real(Det_vec_old(I,nf),kind(0.d0) ) )
           enddo
        enddo
        !Z = cmplx(X,0.d0,kind(0.d0))
        Ratio(1) = cmplx(1.d0,0.d0,kind(0.d0))
        Do nf = 1,N_FL
           !Z = Z*Phase_Det_new(nf)/Phase_Det_old(nf)
           Ratio(1) = Ratio(1) *  Phase_Det_new(nf)/Phase_Det_old(nf)
        enddo
        !Z = Z**N_SUN 
        Ratio(1) = Ratio(1)**N_SUN
        Ratio_2 = real(N_SUN,kind(0.d0))*Ratio_2
        
        Do I = 1,Size(Op_V,1)
           If (Op_V(i,1)%type == 2) then 
              X = 0.d0
              Do nt = 1,Ltrot
                 if ( nsigma(i,nt) /= nsigma_old(i,nt) )  then 
                    !Z = Z * cmplx( Gaml(nsigma(i,nt),2)/Gaml(nsigma_old(i,nt),2),0.d0,kind(0.d0) ) 
                    Ratio(1) = Ratio(1) * cmplx( Gaml(nsigma(i,nt),2)/Gaml(nsigma_old(i,nt),2),0.d0,kind(0.d0) )
                    X = X + Phi(nsigma(i,nt),2) - Phi(nsigma_old(i,nt),2)
                 endif
              Enddo
              Do nf = 1,N_FL
                 !Z = Z * exp(cmplx( X*Real(N_SUN,Kind(0.d0)), 0.d0,kind(0.d0)) * Op_V(i,nf)%g * Op_V(i,nf)%alpha )
                 Ratio(1) = Ratio(1) * exp(cmplx( X*Real(N_SUN,Kind(0.d0)), 0.d0,kind(0.d0)) * Op_V(i,nf)%g * Op_V(i,nf)%alpha )
              Enddo
           endif
        Enddo
        !Z =  Z * cmplx( Delta_S0_global(Nsigma_old),0.d0,kind(0.d0) )
        !Z =  Z * cmplx( T0_Proposal_ratio, 0.d0,kind(0.d0))
        Ratio_2 = Ratio_2 + log(Delta_S0_global(Nsigma_old)) + log(T0_Proposal_ratio)
        
        Ratio(2) = Ratio_2
        Compute_Ratio_Global = Ratio(1)*exp(Ratio(2))
        
        
      end Function Compute_Ratio_Global
      

!--------------------------------------------------------------------
      subroutine Compute_Fermion_Det(Phase,Det_Vec, udvl)
!--------------------------------------------------------------------
!> @author 
!> Fakher Assaad 
!>
!> @brief 
!> Computes det( 1 +  VL*DL*UL)   = Phase * Det_vec(1)*..*Det_vec(Ndim)
!> Note that Phase is a unit complex number and the Det_vec contains only 
!> positive, implying real, numbers.
!>  
!--------------------------------------------------------------------

        Use  UDV_Wrap_mod
        Use UDV_State_mod
        Implicit none
        
        CLASS(UDV_State), INTENT(INOUT) :: udvl
        Complex (Kind=Kind(0.d0)), Dimension(:), Intent(OUT)  ::  Det_Vec
        Complex (Kind=Kind(0.d0)) :: Phase
        
        !> Local variables
        Integer ::  N_size, NCON, J
        COMPLEX (Kind=Kind(0.d0)) :: alpha,beta, Z, Z1
        TYPE(UDV_State) :: udvlocal
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Allocatable ::  TP!, U, V
        !    COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: D
        
        !    N_size = SIZE(DL,1)
        N_size = udvl%ndim
        NCON  = 0
        alpha = cmplx(1.d0,0.d0,kind(0.d0))
        beta  = cmplx(0.d0,0.d0,kind(0.d0))
        Allocate (TP(N_Size,N_Size))
        TP = CT(udvl%U)
        DO J = 1,N_size
           TP(:,J) = TP(:,J) +  udvl%V(:,J)*exp(udvl%L(J))
        ENDDO
        CALL udvlocal%alloc(N_size)
        Call  UDV_WRAP_Pivot(TP,udvlocal%U, udvlocal%D, udvlocal%V, NCON,N_size,N_Size)
        Z  = DET_C(udvlocal%V, N_size) ! Det destroys its argument
        Call MMULT(TP, udvl%U, udvlocal%U)
        Z1 = Det_C(TP, N_size) 
        Deallocate (TP)
        Phase   = Z*Z1/ABS(Z*Z1)
        Det_vec = udvlocal%D
        Det_vec(1) = Det_vec(1)*ABS(Z*Z1)
        
        CALL udvlocal%dealloc
        
      end subroutine Compute_Fermion_Det

      
!--------------------------------------------------------------------
      Integer function  npbc_tempering(n,Isize)
!--------------------------------------------------------------------
!> @author 
!> Fakher Assaad 
!>
!> @brief 
!> Periodic boundary conditions required to defined master and slave for the  
!> tempering
!--------------------------------------------------------------------
        implicit none
        Integer,  INTENT(IN)   :: Isize,n
        
        npbc_tempering = n
        if (  npbc_tempering < 0       ) npbc_tempering = npbc_tempering +  Isize
        if (  npbc_tempering > Isize -1) npbc_tempering = npbc_tempering -  Isize
        
      end function npbc_tempering


!--------------------------------------------------------------------
!> @author 
!> Fakher Assaad 
!>
!> @brief 
!> The following routine monitors the acceptance of tempering moves
!> 
!--------------------------------------------------------------------
      Subroutine Global_Tempering_setup
        Integer    ::   N
        Character (len=64) ::  Filename
        N = 1
        Filename = "Acc_Temp"
        Call Obser_Vec_make(Tempering_acceptance,N,Filename)
      End Subroutine Global_Tempering_setup

!--------------------------------------------------------------------

      Subroutine Global_Tempering_init_obs
        Call Obser_vec_Init( Tempering_acceptance )
      end Subroutine Global_Tempering_init_obs

!--------------------------------------------------------------------

      Subroutine Global_Tempering_obser(toggle)
        Implicit none

        Logical, intent(in) :: toggle
        Tempering_acceptance%N         =  Tempering_acceptance%N + 1
        Tempering_acceptance%Ave_sign  =  Tempering_acceptance%Ave_sign + 1.d0
        if (toggle) Tempering_acceptance%Obs_vec(1) = Tempering_acceptance%Obs_vec(1) +  cmplx(1.d0,0.d0,kind(0.d0))
        
      end Subroutine Global_Tempering_obser
!--------------------------------------------------------------------
      Subroutine Global_Tempering_Pr
        Implicit none
        Call  Print_bin_Vec(Tempering_acceptance,Group_Comm)
      end Subroutine Global_Tempering_Pr
      
      
end Module Global_mod
