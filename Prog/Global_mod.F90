!  Copyright (C) 2016 - 2022 The ALF project
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
!
!> @brief
!> Handles global updates and parallel tempering
!
!--------------------------------------------------------------------

Module Global_mod

      Use runtime_error_mod
      Use Hamiltonian_main
      Use MyMats
      Use Operator_mod
      Use Control
      Use Observables
      Use Fields_mod
      Use Random_Wrap
      use iso_fortran_env, only: output_unit, error_unit

      Implicit none

      Type (Obser_Vec ),  private  ::   Tempering_acceptance


    Contains
#if defined(TEMPERING)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Handles parrallel tempering \n
!> This subroutine is called only if the tempering flag is switched on \n
!> Requires MPI
!>
!> @details
!> On entry and on exit the left storage is full, the Green function is on time slice 0 and the phase is set.
!> The same holds on exit but with updated quantities.
!>
!> @param[inout] Phase  Complex
!> \verbatim
!>  Is updated upon acceptance
!> \endverbatim
!> @param[inout] udvl, udvr Class(UDV_State)
!> \verbatim
!>  On entry udvl contains the last udv decomposition such that G=(1 + VL*DL*UL^{dag})^{-1}
!>  udvr  is used as working space
!> \endverbatim
!> @param[inout] GR Complex
!> \verbatim
!>  Green function. Is updated upon acceptance.
!> \endverbatim
!> @param[inout] udvst Class(UDV_State)
!> \verbatim
!>  Storage. Is updated with left propagation upon acceptance
!> \endverbatim
!> @param[in] Stab_nt Integer
!> \verbatim
!>  List of time slices for stabilization.
!> \endverbatim
!> @param[in] N_exchange_steps Integer
!> \verbatim
!> Number of exchange steps
!> \endverbatim
!> @param[in] Tempering_calc_det Logical
!> \verbatim
!>  If true then the fermion determinant is computed in the calulcation of the exchange weight.
!> \endverbatim
!>
!--------------------------------------------------------------------
      Subroutine Exchange_Step(Phase,GR, udvr, udvl, Stab_nt, udvst, N_exchange_steps, Tempering_calc_det)
        Use UDV_State_mod
        Use mpi
        use wrapul_mod
        use cgr1_mod
        Implicit none

        !  Arguments
        COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT)                                :: Phase
        CLASS(UDV_State), intent(inout), allocatable, Dimension(:)              :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
        CLASS(UDV_State), intent(inout), allocatable, Dimension(:, :) :: udvst
        INTEGER, dimension(:),     INTENT   (IN), allocatable         :: Stab_nt
        INTEGER, INTENT(IN) :: N_exchange_steps
        Logical, INTENT(IN) :: Tempering_calc_det


        !>  Local variables.
        Integer :: NST, NSTM, NF, nf_eff, NT, NT1, NVAR,N, N1,N2, I, NC, I_Partner, n_step,  N_count, N_part
        Type    (Fields)           :: nsigma_old
        Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, Weight, Weight1, delta_S0_log, exp_delta_S0
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Ratiotot, Ratiotot_p, Phase_old, Phase_new
        Real    (Kind=Kind(0.d0)), allocatable :: Det_vec_old(:,:), Det_vec_new(:,:)
        Complex (Kind=Kind(0.d0)), allocatable :: Phase_Det_new(:), Phase_Det_old(:)
        Complex (Kind=Kind(0.d0)) :: Ratio(2), Ratio_p(2),Phase_array(N_FL)
        Logical :: TOGGLE, L_Test
        Integer, allocatable :: List_partner(:), List_masters(:)

        ! Keep track of where the configuration originally came from
        Integer        :: nsigma_irank, nsigma_old_irank, nsigma_irank_temp
        Integer        :: n_GR

        !Integer, Dimension(:,:),  allocatable :: nsigma_orig, nsigma_test
        !Integer :: I1, I2

        Integer        :: Isize, Irank, Ierr, irank_g, isize_g, igroup
        Integer        :: STATUS(MPI_STATUS_SIZE)

        Character (Len=64)  :: storage


        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
        nsigma_irank     = irank

        n1 = size(nsigma%f,1)
        n2 = size(nsigma%f,2)
        NSTM = Size(udvst, 1)
        call nsigma_old%make(n1, n2)
        if (Tempering_calc_det) then
           Allocate ( Det_vec_old(NDIM,N_FL), Det_vec_new(NDIM,N_FL) )
           Allocate ( Phase_Det_new(N_FL), Phase_Det_old(N_FL) )
        endif
        Allocate ( List_partner(0:Isize-1), List_masters(Isize/2) )

        !         Allocate ( nsigma_orig(n1,n2) )
        !         nsigma_orig = nsigma

        !  Compute for each core the old weights.
        L_test = .false.
        if (Tempering_calc_det) then
           ! Set old weight.
           storage = "Full"
           Call Compute_Fermion_Det(Phase_det_old,Det_Vec_old, udvl, udvst, Stab_nt, storage)
           Phase_array=1.0d0
           Do nf_eff = 1,N_Fl_eff
              nf=Calc_Fl_map(nf_eff)
              Phase_array(nf) = Phase_det_old(nf)
              Call Op_phase(Phase_array(nf),OP_V,Nsigma,nf)
           Enddo
           if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
           Phase_old=product(Phase_array)
           Phase_old=Phase_old**N_SUN
        endif
        !> Store old configuration
        nsigma_old%f = nsigma%f
        nsigma_old%t = nsigma%t
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
           !  Set the partner rank on each core
           If (Irank == 0 ) then
              n_step = isize_g
              if (  ranf_wrap() > 0.5d0 ) n_step = -isize_g
              Do I = 0,Isize-1,2*isize_g
                 do n = 0,isize_g-1
                    List_partner(npbc_tempering(I  + n             ,Isize)) =  npbc_tempering(I + n   + n_step ,Isize)
                    List_partner(npbc_tempering(I  + n  + n_step  , Isize)) =  npbc_tempering(I + n            ,Isize)
                 enddo
              enddo
           endif

           CALL MPI_BCAST(List_partner, Isize  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)

           If (L_test) then
              if (Irank == 0 ) then
                 Write(6,*) 'Testing global : '
                 do  I = 0,Isize -1
                    Write(6,*)   I, List_partner(I) !, Phase_old, Phase
                    Write(11,*)  I, List_partner(I) !, Phase_old, Phase
                 enddo
                 Write(11,*)
              endif
           Endif


           !  Exchange configurations
           !  The types do not change --> no need to exchange them
           n = size(nsigma_old%f,1)*size(nsigma_old%f,2)
           CALL MPI_Sendrecv(nsigma_old%f    , n, MPI_COMPLEX16, List_partner(IRANK), 0, &
                    &        nsigma%f        , n, MPI_COMPLEX16, List_partner(IRANK), 0, MPI_COMM_WORLD,STATUS,IERR)
           CALL MPI_Sendrecv(nsigma_old_irank, 1, MPI_INTEGER, List_partner(IRANK), 0, &
                    &        nsigma_irank    , 1, MPI_INTEGER, List_partner(IRANK), 0, MPI_COMM_WORLD,STATUS,IERR)

           !  Each node now has a new configuration nsigma


           if (Tempering_calc_det) then
              !  HERE
              !  Compute ratio on weights one each rank
              storage = "Empty"
              Call Compute_Fermion_Det(Phase_det_new,Det_Vec_new, udvl, udvst, Stab_nt, storage)

              Phase_array = cmplx(1.d0,0.d0,kind=kind(0.d0))
              Do nf_eff = 1,N_Fl_eff
                 nf=Calc_FL_map(nf_eff)
                 Phase_array(nf) = Phase_det_new(nf)
                 Call Op_phase(Phase_array(nf),OP_V,Nsigma,nf)
              Enddo
              if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
              Phase_new=product(Phase_array)
              Phase_new=Phase_new**N_SUN

              T0_Proposal_ratio = 1.d0
              Ratiotot = Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
                   &            Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio,Ratio)

              If (L_Test) Write(6,*) 'Ratio_global: Irank, Partner',Irank,List_partner(Irank), &
                   &                  Ratiotot, Ratio(1)*exp(Ratio(2))
           else
              delta_S0_log = ham%Get_Delta_S0_global(Nsigma_old)
              Ratio(1) = 1.0d0 !!! @FAKHER @JONAS double check please
              Ratio(2) = delta_S0_log
           endif

           !  Acceptace/rejection decision taken on master node after receiving information from slave
           Do nc = 1,Isize/2 ! Loop over masters
              I = List_masters(nc)
              If (Irank == I ) Then
                 CALL MPI_SEND(Ratio   ,2, MPI_COMPLEX16  , List_partner(I), I+1024, MPI_COMM_WORLD,IERR)
              else if (IRANK == List_Partner(I) ) Then
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
                 CALL MPI_SEND(Toggle, 1, MPI_LOGICAL, I, I+1024, MPI_COMM_WORLD,IERR)
              else if (IRANK == I ) Then
                 CALL MPI_RECV(Toggle , 1, MPI_LOGICAL,   List_partner(I), I+1024 ,MPI_COMM_WORLD,STATUS,IERR)
                 If (L_Test) Write(6,*) 'Slave : ', Irank,  Toggle
              endif
           enddo

           If (L_Test) Write(6,*) 'Final: ',  Irank, List_partner(Irank), toggle
           Call Global_Tempering_obser(Toggle)
           Call Control_upgrade_Temp  (Toggle)
           If (toggle)  then
              !     Move has been accepted
              if (Tempering_calc_det) then
                 Phase_old     = Phase_new
                 Phase_det_old = Phase_det_new
                 Det_vec_old   = Det_vec_new
              endif
              nsigma_old%f       = nsigma%f
              nsigma_old%t       = nsigma%t
              nsigma_old_irank = nsigma_irank
           else
              nsigma%f       = nsigma_old%f
              nsigma%t       = nsigma_old%t
              nsigma_irank = nsigma_old_irank
           endif
        enddo

        ! Finalize
        if (Tempering_calc_det) then
           ! If move has been accepted, no use to recomute storage
           If (.not.TOGGLE) then
              DO nf_eff = 1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 if (Projector) then
                    CALL udvl(nf_eff)%reset('l',WF_L(nf)%P)
                 else
                    CALL udvl(nf_eff)%reset('l')
                 endif
              ENDDO
              DO NST = NSTM-1,1,-1
                 NT1 = Stab_nt(NST+1)
                 NT  = Stab_nt(NST  )
                 !Write(6,*) NT1,NT, NST
                 CALL WRAPUL(NT1,NT, udvl)
                 Do nf_eff = 1,N_FL_eff
                    udvst(NST, nf_eff) = udvl(nf_eff)
                 ENDDO
              ENDDO
              NT1 = stab_nt(1)
              CALL WRAPUL(NT1,0, udvl)
           Endif
           ! Compute the Green functions so as to provide correct starting point for the sequential updates.
           NVAR  = 1
           Phase_array = cmplx(1.d0,0.d0,kind(0.d0))
           do nf_eff = 1,N_Fl_eff
              nf=Calc_Fl_map(nf_eff)
              CALL CGR(Z, NVAR, GR(:,:,nf), udvr(nf_eff),  udvl(nf_eff))
              call Op_phase(Z,OP_V,Nsigma,nf)
              Phase_array(nf)=Z
           Enddo
           if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
           Phase=product(Phase_array)
           Phase=Phase**N_SUN
        else
           !  Send >> Phase, GR, udvr, udvl, udvst << to new node
           !  First step: Each node sends to IRANK=0 its value nsigma_irank,
           !  which is the node where its new Phase, GR, udvr, udvl, udvst is stored
           !  This node then tells each node where to send its now old Phase, GR, udvr, udvl, udvst
           !  Finally, the variables get submitted
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
              CALL MPI_Sendrecv_Replace(Phase, 1, MPI_COMPLEX16, nsigma_old_irank, 0, &
                   &        nsigma_irank, 0, MPI_COMM_WORLD, STATUS, IERR)

              n_GR = size(GR,1)*size(GR,2)*size(GR,3)
              CALL MPI_Sendrecv_Replace(GR, n_GR, MPI_COMPLEX16, nsigma_old_irank, 0, &
                   &        nsigma_irank, 0, MPI_COMM_WORLD, STATUS, IERR)

              do nf_eff = 1,N_Fl_eff
                 CALL udvr(nf_eff)%MPI_Sendrecv(nsigma_old_irank, 0, nsigma_irank, 0, STATUS, IERR)
              enddo
              do nf_eff = 1,N_Fl_eff
                 CALL udvl(nf_eff)%MPI_Sendrecv(nsigma_old_irank, 0, nsigma_irank, 0, STATUS, IERR)
              enddo
              do NST = 1, NSTM
                 do nf_eff = 1,N_Fl_eff
                    CALL udvst(NST, nf_eff)%MPI_Sendrecv(nsigma_old_irank, 0, nsigma_irank, 0, STATUS, IERR)
                 enddo
              enddo
           endif
        endif


        call nsigma_old%clear

        if (Tempering_calc_det) then
           Deallocate ( Det_vec_old, Det_vec_new )
           Deallocate ( Phase_Det_new, Phase_Det_old )
        endif
        Deallocate ( List_partner, List_masters )

      end Subroutine Exchange_Step
#endif


!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Carries out N_global global updates as defeined in the Global_move routine of the Hamiltonian module
!>
!> @details
!> On entry and on exit the left storage is full, the Green function is on time slice 0 and the phase is set.
!> The same holds on exit but with updated quantities.
!>
!> @param[inout] Phase  Complex
!> \verbatim
!>  Is updated upon acceptance
!> \endverbatim
!> @param[inout] udvl, udvr Class(UDV_State)
!> \verbatim
!>  On entry udvl contains the last udv decomposition such that G=(1 + VL*DL*UL^{dag})^{-1}
!>  udvr  is used as working space
!> \endverbatim
!> @param[inout] GR Complex
!> \verbatim
!>  Green function. Is updated upon acceptance.
!> \endverbatim
!> @param[inout] udvst Class(UDV_State)
!> \verbatim
!>  Storage. Is updated with left propagation upon acceptance
!> \endverbatim
!> @param[in] Stab_nt Integer
!> \verbatim
!>  List of time slices for stabilization.
!> \endverbatim
!> @param[in] N_Global Integer
!> \verbatim
!>  Number of global moves that will be caried out
!> \endverbatim
!>
!--------------------------------------------------------------------
      Subroutine Global_Updates(Phase,GR, udvr, udvl, Stab_nt, udvst,N_Global)

        Use UDV_State_mod
        use wrapul_mod
        use cgr1_mod
        Implicit none

        !  Arguments
        COMPLEX (Kind=Kind(0.d0)),                                INTENT(INOUT) :: Phase
        CLASS   (UDV_State), DIMENSION(:), ALLOCATABLE,           INTENT(INOUT) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:),  allocatable,INTENT(INOUT) :: GR
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE,            INTENT(INOUT) :: udvst
        INTEGER, dimension(:),   allocatable,                     INTENT(IN)    :: Stab_nt
        Integer,                                                  INTENT(IN)    :: N_Global


        !  Local variables.
        Integer :: NST, NSTM, NF, NT, NT1, NVAR,N, N1,N2, I, NC, N_part,j, nf_eff
        Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, Weight
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Ratiotot, Phase_old, Phase_new
        Complex (Kind=Kind(0.d0)), allocatable :: Det_vec_test(:,:), Phase_Det_new(:), Phase_Det_old(:)
        Real    (Kind=Kind(0.d0)), allocatable :: Det_vec_old(:,:), Det_vec_new(:,:)
        Type   (Fields)   :: nsigma_old

        Complex (Kind=Kind(0.d0)) :: Ratio(2), Phase_array(N_Fl)
        Logical :: TOGGLE, L_Test
        Real    (Kind=Kind(0.d0)) :: size_clust
        Real    (Kind=Kind(0.d0)) :: ratio_2_test
        Character (Len=64)  :: storage



        n1 = size(nsigma%f,1)
        n2 = size(nsigma%f,2)
        NSTM = Size(udvst, 1)
        call nsigma_old%make(n1, n2)

        Allocate ( Det_vec_old(NDIM,N_FL_eff), Det_vec_new(NDIM,N_FL_eff), Det_vec_test(NDIM,N_FL_eff) )
        Allocate ( Phase_Det_new(N_FL_eff), Phase_Det_old(N_FL_eff) )

        L_test = .false.
        ! Set old weight.
        storage = "Full"
        Call Compute_Fermion_Det(Phase_det_old,Det_Vec_old, udvl, udvst, Stab_nt, storage)
        Phase_array = cmplx(1.d0,0.d0,kind=kind(0.d0))
        Do nf_eff = 1,N_Fl_eff
           nf=Calc_FL_map(nf_eff)
           Phase_array(nf) = Phase_det_old(nf_eff)
           Call Op_phase(Phase_array(nf),OP_V,Nsigma,nf)
        Enddo
        if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
        Phase_old=product(Phase_array)
        Phase_old=Phase_old**N_SUN

        If (L_test) then
           ! Testing
           Do nf_eff = 1,N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              if (Projector) then
                 CALL udvr(nf_eff)%reset('r',WF_R(nf)%P)
              else
                 CALL udvr(nf_eff)%reset('r')
              endif
           Enddo
           NVAR = 1
           Phase_array = cmplx(1.d0,0.d0,kind(0.d0))
           do nf_eff = 1,N_Fl_eff
              nf=Calc_Fl_map(nf_eff)
              CALL CGR(Z, NVAR, GR(:,:,nf), udvr(nf_eff), udvl(nf_eff))
              call Op_phase(Z,OP_V,Nsigma,nf)
              Phase_array(nf)=Z
           Enddo
           if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
           Phase=product(Phase_array)
           Phase=Phase**N_SUN
           Do nf_eff = 1,N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              Call DET_C_LU(GR(:,:,nf),Det_vec_test(:,nf_eff),Ndim)
              !!!ATTENTION HOW DOES BELOW DEPEND ON NF BLOCK symm
              Z = Phase_det_old(nf_eff)
              ratio_2_test=0.d0
              DO I = 1,Ndim
                 Z = Z*Det_vec_test(I,nf_eff)/ABS(Det_vec_test(I,nf_eff))
                 ratio_2_test=ratio_2_test+log(ABS(Det_vec_test(I,nf_eff)))+Det_vec_old(I,nf_eff)
              Enddo
              Z=Z*cmplx(exp(ratio_2_test),0.d0,kind(0.d0))
              Write(6,*) 'Testing weight: ', Z
              if(abs(ratio_2_test)>650) write(6,*) "Weight is about to reach double underflow!"
           Enddo
        Endif

        ! Store old configuration
        nsigma_old%f = nsigma%f
        nsigma_old%t = nsigma%t
        ! Phase_old, Phase_det_old and Det_vec_old  are all set.
        NC = 0
        Do n = 1,N_Global
           ! Draw a new spin configuration. This is provided by the user in the Hamiltonian module
           ! Note that nsigma is a variable in the module Hamiltonian
           Call ham%Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
           If (T0_Proposal_ratio > 0.d0) then
              NC = NC + 1
              ! Compute the new Green function
              storage = "Empty"
              Call Compute_Fermion_Det(Phase_det_new,Det_Vec_new, udvl, udvst, Stab_nt, storage)

              Phase_array = cmplx(1.d0,0.d0,kind=kind(0.d0))
              Do nf_eff = 1,N_Fl_eff
                 nf=Calc_FL_map(nf_eff)
                 Phase_array(nf) = Phase_det_new(nf)
                 Call Op_phase(Phase_array(nf),OP_V,Nsigma,nf)
              Enddo
              if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
              Phase_new=product(Phase_array)
              Phase_new=Phase_new**N_SUN

              Ratiotot = Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
                   &                          Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio, Ratio)

              !Write(6,*) 'Ratio_global: ', Ratiotot

              Weight = abs(  real( Phase_old * Ratiotot, kind=Kind(0.d0))/real(Phase_old,kind=Kind(0.d0)) )

              Z = Phase_old * Ratiotot/ABS(Ratiotot)
              Call Control_PrecisionP_Glob(Z,Phase_new)
              !Write(6,*) Z, Phase_new


              TOGGLE = .false.
              if ( Weight > ranf_wrap() )  Then
                 TOGGLE = .true.
                 Phase_old     = Phase_new
                 Phase_det_old = Phase_det_new
                 nsigma_old%t    = nsigma%t
                 nsigma_old%f    = nsigma%f
                 Det_vec_old     = Det_vec_new
              else
                 nsigma%t = nsigma_old%t
                 nsigma%f = nsigma_old%f
              endif
              Call Control_upgrade_Glob(TOGGLE,size_clust)
           else
              nsigma%t = nsigma_old%t
              nsigma%f = nsigma_old%f
           endif
        Enddo

        If (NC > 0 ) then
           If (.not.TOGGLE) then
              DO nf_eff = 1,N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                if (Projector) then
                  CALL udvl(nf_eff)%reset('l',WF_L(nf)%P)
                else
                  CALL udvl(nf_eff)%reset('l')
                endif
              ENDDO
              DO NST = NSTM-1,1,-1
                 NT1 = Stab_nt(NST+1)
                 NT  = Stab_nt(NST  )
                 !Write(6,*) NT1,NT, NST
                 CALL WRAPUL(NT1,NT, udvl)
                 Do nf_eff = 1,N_FL_eff
                    udvst(NST, nf_eff) = udvl(nf_eff)
                 ENDDO
              ENDDO
              NT1 = stab_nt(1)
              CALL WRAPUL(NT1,0, udvl)
           Endif
           !Compute the Green functions so as to provide correct starting point for the sequential updates.
           NVAR  = 1
           Phase_array = cmplx(1.d0,0.d0,kind(0.d0))
           do nf_eff = 1,N_Fl_eff
              nf=Calc_Fl_map(nf_eff)
              CALL CGR(Z, NVAR, GR(:,:,nf), udvr(nf_eff), udvl(nf_eff))
              call Op_phase(Z,OP_V,Nsigma,nf)
              Phase_array(nf)=Z
           Enddo
           if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
           Phase=product(Phase_array)
           Phase=Phase**N_SUN
        endif


        call nsigma_old%clear
        Deallocate ( Det_vec_old  , Det_vec_new, Det_vec_test  )
        Deallocate ( Phase_Det_new, Phase_Det_old )


      End Subroutine Global_Updates


!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> This fucntion computes ratio of weights for global moves
!>
!> @details
!> \verbatim
!>  Ratio_Global = T0(nsigma--> nsigma_old) W(nsigma)/T0(nsigma_old--> nsigma) W(nsigma_old) =
!>                 T0_Proposal_ratio   W(nsigma)/  W(nsigma_old)
!>  On entry the old and new fermion determinants read: Phase_det*e^{ \sum_{n=1}^{ndim} Det_vec(n) }
!>  as obtained from the Compute_Fermion_det routine.
!>  On exit Ratio_Global = Ratio(1)*exp(Ratio(2))
!> \endverbatim
!>
!> @param[in]  Phase_det_new  Complex, Dimension(N_FL)
!> @param[in]  Phase_det_old  Complex, Dimension(N_FL)
!> @param[in]  Det_vec_new  Real, Dimension(:,N_FL)
!> @param[in]  Det_vec_old  Real, Dimension(:,N_FL)
!> @param[in]  T0_proposal_ratio   Real
!> @param[in]  nsigma_old Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma. nsigma is a globale variable
!>  contained in the Hamiltonian module.
!> \endverbatim
!--------------------------------------------------------------------
      Complex (Kind=Kind(0.d0)) Function Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
           &                    Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio,Ratio)


        Implicit none

        ! Arguments
        Complex (Kind=Kind(0.d0)), allocatable, INTENT(IN) :: Phase_Det_old(:), Phase_Det_new(:)
        REAL    (Kind=Kind(0.d0)), allocatable, INTENT(IN) :: Det_vec_old(:,:), Det_vec_new(:,:)
        Real    (Kind=Kind(0.d0)),    INTENT(IN)  :: T0_proposal_ratio
        Type    (Fields),             INTENT(IN)  :: nsigma_old
        Complex (Kind=Kind(0.d0)),    INTENT(out) :: Ratio(2)

        ! Local
        Integer  :: Nf, i, nt, nf_eff
        Complex (Kind=Kind(0.d0)) :: Z, Z1, Ratio_1_array(N_FL), Ratio_2_array(N_FL), g_loc
        Real    (Kind=Kind(0.d0)) :: X, Ratio_2, delta, log_delta

        Ratio = cmplx(0.d0,0.d0,kind(0.d0))
        Ratio_2_array = 0.d0
        !X = 1.d0
        Do nf_eff = 1,N_Fl_eff
           nf=Calc_Fl_map(nf_eff)
           DO I = 1,Ndim
              !X= X*real(Det_vec_new(I,nf),kind(0.d0)) / Real(Det_vec_old(I,nf),kind(0.d0) )
              Ratio_2_array(nf) = Ratio_2_array(nf) +  Det_vec_new(I,nf_eff) - Det_vec_old(I,nf_eff)
           enddo
           !Z = Z**N_SUN
           Ratio_2_array(nf) = real(N_SUN,kind(0.d0))*Ratio_2_array(nf)
        enddo
        !Z = cmplx(X,0.d0,kind(0.d0))
        Ratio_1_array = cmplx(1.d0,0.d0,kind(0.d0))
        Do nf_eff = 1,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           !Z = Z*Phase_Det_new(nf)/Phase_Det_old(nf)
           Ratio_1_array(nf_eff) = Ratio_1_array(nf) *  Phase_Det_new(nf_eff)/Phase_Det_old(nf_eff)
           !Z = Z**N_SUN
           Ratio_1_array(nf) = Ratio_1_array(nf)**N_SUN
        enddo

        Ratio(1)=1.d0
        Do I = 1,Size(Op_V,1)
           X = 0.d0
           Do nt = 1,Ltrot
              !Z = Z * cmplx( Gaml(nsigma(i,nt),2)/Gaml(nsigma_old(i,nt),2),0.d0,kind(0.d0) )
              Ratio(1) = Ratio(1) * cmplx( nsigma%Gama(i,nt)/nsigma_old%Gama(i,nt),0.d0,kind(0.d0) ) !You could put this in Ratio_2
              !X = X + nsigma%Phi(i,nt) - nsigma_old%Phi(i,nt)
              Do nf_eff = 1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 Z =  nsigma%Phi(i,nt) - nsigma_old%Phi(i,nt)
                 g_loc = Op_V(i,nf)%g
                 if ( Op_V(i,nf)%get_g_t_alloc()) g_loc = Op_V(i,nf)%g_t(nt)
                 !Z = Z * exp(cmplx( X*Real(N_SUN,Kind(0.d0)), 0.d0,kind(0.d0)) * Op_V(i,nf)%g * Op_V(i,nf)%alpha )
                 Ratio_1_array(nf) = Ratio_1_array(nf) * exp(Z*cmplx( Real(N_SUN,Kind(0.d0)), 0.d0,kind(0.d0)) * g_loc * Op_V(i,nf)%alpha )
              Enddo
           Enddo
        Enddo
        if (reconstruction_needed) then
            call ham%weight_reconstruction(Ratio_1_array)
            call ham%weight_reconstruction(Ratio_2_array)
        endif
        Ratio(1)=Ratio(1)*product(Ratio_1_array)
        !Z =  Z * cmplx( ham%Delta_S0_global(Nsigma_old),0.d0,kind(0.d0) )
        !Z =  Z * cmplx( T0_Proposal_ratio, 0.d0,kind(0.d0))
        Ratio(2) = sum(Ratio_2_array)
        log_delta = ham%Get_Delta_S0_global(Nsigma_old)
        Ratio(2) = Ratio(2) + log_delta + log(T0_Proposal_ratio)

        Compute_Ratio_Global = Ratio(1)*exp(Ratio(2))


      end Function Compute_Ratio_Global

!--------------------------------------------------------------------
!> @author
!>
!> @brief
!> Computes the fermion determinant from scratch or from the storage.
!>
!> @details
!> \verbatim
!> Finite temperature: If the storage is full then computes 
!>                     det( 1 +  VL*DL*UL^{dag} )= Phase_det*e^{ \sum_{n=1}^{ndim} Det_vec(n)}
!>                     Note that  Udvl%U = UL, Udvl%D = DL and  Udvl%V = VL^{dag}
!>                     If storage is empty one first computes  Udvl and in doing so fills up the storage.
!> Zero temperature:   If storage is full then computes  det \Psi_L U(\theta_tot,0) \Psi_R  =  Phase_det*e^{ \sum_{n=1}^{N_part} Det_vec(n) }
!>                     For the projective code the required scales which one can throw away for the Green function are in udvst.
!>                     If the storage is not full, then things will be computed from scratch and intermediate results will be stored in udvst.
!> \endverbatim
!>
!> @param[inout] udvl, udvr Class(UDV_State)
!> \verbatim
!>  If storage=Full  : On entry and exit udvl contains the last udv decomposition such that G=(1 + VL*DL*UL^{dag})^{-1}
!>  If storage=Empty : On exit udv1 will contain the last udv decomposition such that G=(1 + VL*DL*UL^{dag})^{-1}
!> \endverbatim
!> @param[out] Phase_det  Complex, Dimension(N_FL)
!> \verbatim
!>  The phase, per flavor
!> \endverbatim
!> @param[out] Det_vec  Real, Dimension(:,N_FL)
!> @param[in] Stab_nt  Integer, Dimension(:)
!> \verbatim
!>  List of time slices for stabilization.
!> \endverbatim
!> @param[in] Storage Character
!> \verbatim
!>  storage = "Full".  Compute the fermion det  with the use of the storage.
!>  storage = "Empty". Compute the fermion det from scratch and  in doing so, fill the storage.
!> \endverbatim
!>
!--------------------------------------------------------------------
      Subroutine Compute_Fermion_Det(Phase_det,Det_Vec, udvl, udvst, Stab_nt, storage)
!--------------------------------------------------------------------

        Use UDV_Wrap_mod
        Use UDV_State_mod
        use wrapul_mod

        Implicit none

        REAL    (Kind=Kind(0.d0)), Dimension(:,:), Intent(OUT)  ::  Det_Vec
        Complex (Kind=Kind(0.d0)), Dimension(:)  , Intent(OUT)  ::  Phase_det
        CLASS(UDV_State), DIMENSION(:)  , ALLOCATABLE,  INTENT(INOUT) :: udvl
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE,  INTENT(INOUT) :: udvst
        INTEGER,          dimension(:),   allocatable,    INTENT(IN)  :: Stab_nt
        Character (Len=64), Intent(IN) :: storage


        !> Local variables
        Integer ::  N_size, NCON, I,  J, N_part, info, NSTM, N, nf, nst, nt, nt1, nf_eff
        Integer, allocatable :: ipiv(:)
        COMPLEX (Kind=Kind(0.d0)) :: alpha,beta, Z, Z1
        TYPE(UDV_State) :: udvlocal
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Allocatable ::  TP!, U, V
        COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: D

        !! TODO adapt to flavor symmetry

        if(udvl(1)%side .ne. "L" .and. udvl(1)%side .ne. "l" ) then
           write(error_unit,*) "Compute_Fermion_Det: calling wrong decompose"
           CALL Terminate_on_error(ERROR_GLOBAL_UPDATES,__FILE__,__LINE__)
        endif

        NSTM = Size(udvst, 1)
        If (str_to_upper(storage) == "EMPTY" ) then
           DO nf_eff = 1,N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              if (Projector) then
                 CALL udvl(nf_eff)%reset('l',WF_L(nf)%P)
              else
                 CALL udvl(nf_eff)%reset('l')
              endif
           ENDDO
           DO NST = NSTM-1,1,-1
              NT1 = Stab_nt(NST+1)
              NT  = Stab_nt(NST  )
              !Write(6,*) 'Call wrapul', NT1,NT, udvl(1)%d(1), udvl(1)%N_part, udvl(1)%Ndim
              CALL WRAPUL(NT1,NT, udvl)
              Do nf_eff = 1,N_FL_eff
                 udvst(NST, nf_eff) = udvl(nf_eff)
              ENDDO
           ENDDO
           NT1 = stab_nt(1)
           CALL WRAPUL(NT1,0, udvl)
        Endif

        if (Projector) then
           N_part=udvst(1,1)%N_part
           N_size=udvl(1)%ndim
           Det_vec = 0.d0
           Allocate (TP(N_part,N_part), ipiv(N_part))
           do nf_eff=1,N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              N_part=udvst(1,nf)%N_part
              do i=1,NSTM-1
                 do n=1,N_part
#if !defined(STABLOG)
                    Det_Vec(n,nf_eff)=Det_Vec(n,nf_eff)+log(dble(udvst(i,nf)%D(n)))
#else
                    Det_Vec(n,nf_eff)=Det_Vec(n,nf_eff)+udvst(i,nf)%L(n)
#endif
                 enddo
              enddo
              do n=1,N_part
#if !defined(STABLOG)
                 Det_Vec(n,nf_eff)=Det_Vec(n,nf_eff)+log(dble(udvl(nf_eff)%D(n)))
#else
                 Det_Vec(n,nf_eff)=Det_Vec(n,nf_eff)+udvl(nf_eff)%L(n)
#endif
              enddo
              alpha=1.d0
              beta=0.d0
              CALL ZGEMM('C','N',N_part,N_part,N_size,alpha,udvl(nf_eff)%U(1,1),N_size,WF_R(nf)%P(1,1),N_size,beta,TP(1,1),N_part)
              ! ZGETRF computes an LU factorization of a general M-by-N matrix A
              ! using partial pivoting with row interchanges.
              call ZGETRF(N_part, N_part, TP(1,1), N_part, ipiv, info)
              Z=1.d0
              Do J=1,N_part
                 if (ipiv(J).ne.J) then
                    Z = -Z
                 endif
                 Z =  Z * TP(J,J)
              enddo
              Phase_det(nf_eff) = Z/abs(Z)
              Det_vec(1,nf) = Det_vec(1,nf)+log(abs(Z))
           enddo
           Deallocate(TP,ipiv)
           return
        endif

        !    N_size = SIZE(DL,1)
        N_size = udvl(1)%ndim
        NCON  = 0
        alpha = cmplx(1.d0,0.d0,kind(0.d0))
        beta  = cmplx(0.d0,0.d0,kind(0.d0))
        Allocate (TP(N_Size,N_Size),D(N_size))
        CALL udvlocal%alloc(N_size)
        Do nf_eff = 1,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           TP = udvl(nf_eff)%U !udvl stores U^dag instead of U !CT(udvl%U)
#if !defined(STABLOG)
#if !defined(STAB3)
           DO J = 1,N_size
              TP(:,J) = TP(:,J) +  udvl(nf_eff)%V(:,J)*udvl(nf_eff)%D(J)
           ENDDO
#else
           DO J = 1,N_size
              if ( dble(udvl(nf_eff)%D(J)) <= 1.d0 ) then
                 TP(:,J) = TP(:,J) +  udvl(nf_eff)%V(:,J)*udvl(nf_eff)%D(J)
              else
                 TP(:,J) = TP(:,J)/udvl(nf_eff)%D(J) +  udvl(nf_eff)%V(:,J)
              endif
           ENDDO
#endif
#else
           DO J = 1,N_size
              if ( udvl(nf_eff)%L(J) <= 0.d0 ) then
                 TP(:,J) = TP(:,J) +  udvl(nf_eff)%V(:,J)*cmplx(exp(udvl(nf_eff)%L(J)),0.d0,kind(0.d0))
              else
                 TP(:,J) = TP(:,J)*cmplx(exp(-udvl(nf_eff)%L(J)),0.d0,kind(0.d0)) +  udvl(nf_eff)%V(:,J)
              endif
           ENDDO
#endif

           !Call  UDV_WRAP_Pivot(TP,udvlocal%U, D, udvlocal%V, NCON,N_size,N_Size)
           Call  UDV_WRAP(TP,udvlocal%U, D, udvlocal%V, NCON )
           Z  = DET_C(udvlocal%V, N_size) ! Det destroys its argument
           !         Call MMULT(TP, udvl%U, udvlocal%U)
           CALL ZGEMM('C', 'N', N_size, N_size, N_size, alpha, udvl(nf_eff)%U(1,1), N_size, udvlocal%U(1,1), N_size, beta, TP, N_size)
           Z1 = Det_C(TP, N_size)
           Phase_det(nf_eff)   = Z*Z1/ABS(Z*Z1)
#if !defined(STABLOG)
#if !defined(STAB3)
           Det_vec(:,nf_eff) = log(real(D(:)))
           Det_vec(1,nf_eff) = log(real(D(1))) + log(ABS(Z*Z1))
#else
           Det_vec(1,nf_eff) = log(real(D(1))*ABS(Z*Z1))
           if (dble(udvl(nf_eff)%D(1)) > 1.d0) Det_vec(1,nf_eff)=Det_Vec(1,nf_eff)+log(dble(udvl(nf_eff)%D(1)))
           Do J=2,Ndim
              if (dble(udvl(nf_eff)%D(J))<=1.d0) then
                 Det_vec(J,nf_eff) = log(real(D(J)))
              else
                 Det_vec(J,nf_eff) = log(real(D(J)))+log(dble(udvl(nf_eff)%D(J)))
              endif
           enddo
#endif
#else
           Det_vec(:,nf_eff) = log(real(D(:)))
           Det_vec(1,nf_eff) = log(real(D(1))*ABS(Z*Z1))
           if (udvl(nf_eff)%L(1) > 0.d0) Det_vec(1,nf_eff)=Det_Vec(1,nf_eff)+udvl(nf_eff)%L(1)
           Do J=2,Ndim
              if (udvl(nf_eff)%L(J)<=0.d0) then
                 Det_vec(J,nf_eff) = log(real(D(J)))
              else
                 Det_vec(J,nf_eff) = log(real(D(J)))+udvl(nf_eff)%L(J)
              endif
           enddo
#endif
        Enddo
        Deallocate (TP)
        CALL udvlocal%dealloc

      end subroutine Compute_Fermion_Det


!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Periodic boundary conditions required to defined master and slave for the tempering
!--------------------------------------------------------------------
      Integer function  npbc_tempering(n,Isize)
        implicit none
        Integer,  INTENT(IN)   :: Isize,n

        npbc_tempering = n
        if (  npbc_tempering < 0       ) npbc_tempering = npbc_tempering +  Isize
        if (  npbc_tempering > Isize -1) npbc_tempering = npbc_tempering -  Isize

      end function npbc_tempering


!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Allocates the Tempering_acceptance (Type Obser_vec) variable  that monitors the exchange acceptance
!--------------------------------------------------------------------
      Subroutine Global_Tempering_setup
        Integer    ::   N
        Character (len=64) ::  Filename
        N = 1
        Filename = "Acc_Temp"
        Call Obser_Vec_make(Tempering_acceptance,N,Filename)
      End Subroutine Global_Tempering_setup

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Initializes  Tempering_acceptance
!--------------------------------------------------------------------

      Subroutine Global_Tempering_init_obs
        Call Obser_vec_Init( Tempering_acceptance )
      end Subroutine Global_Tempering_init_obs

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Measures Tempering_acceptance
!--------------------------------------------------------------------

      Subroutine Global_Tempering_obser(toggle)
        Implicit none

        Logical, intent(in) :: toggle
        Tempering_acceptance%N         =  Tempering_acceptance%N + 1
        Tempering_acceptance%Ave_sign  =  Tempering_acceptance%Ave_sign + 1.d0
        if (toggle) Tempering_acceptance%Obs_vec(1) = Tempering_acceptance%Obs_vec(1) +  cmplx(1.d0,0.d0,kind(0.d0))

      end Subroutine Global_Tempering_obser
!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Prints Tempering_acceptance
!--------------------------------------------------------------------
      Subroutine Global_Tempering_Pr
        Implicit none
        Call  Print_bin_Vec(Tempering_acceptance,Group_Comm)
      end Subroutine Global_Tempering_Pr


end Module Global_mod
