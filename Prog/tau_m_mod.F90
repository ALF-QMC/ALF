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
!> This module  handles calculation of imagimary time displaced Green functions and  
!> calls the routine ObserT.F90 in the Hamiltonian module,  so as to compute the 
!> defined  time dispalced correlations functions. This modules is for the finite temperature code.
!> 
!
!--------------------------------------------------------------------
     Module Tau_m_mod
       Use Hamiltonian_main 
       Use Operator_mod
       Use Control
       Use Hop_mod
       Use UDV_State_mod
       Use Langevin_HMC_mod
       use wrapur_mod
       use cgr2_2_mod
       

       Contains

         SUBROUTINE TAU_M(udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN) 
           
           Implicit none
     
           Integer, Intent(In) :: NSTM, NWRAP
           CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN) :: udvst
           Complex (Kind=Kind(0.d0)), Intent(in) :: GR(NDIM,NDIM,N_FL),  Phase
           Integer, Intent(In) :: STAB_NT(0:NSTM)         
           Integer, Intent(In) :: LOBS_ST, LOBS_EN

           ! Local 
           ! This could be placed as  private for the module 
           Complex (Kind=Kind(0.d0))  :: GT0(NDIM,NDIM,N_FL),  G00(NDIM,NDIM,N_FL), GTT(NDIM,NDIM,N_FL), G0T(NDIM,NDIM,N_FL)
           Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable  :: GT0_T,  G00_T, GTT_T, G0T_T
           CLASS(UDV_State), DIMENSION(:), ALLOCATABLE :: udvr
           Complex (Kind=Kind(0.d0))  :: HLP4(Ndim,Ndim), HLP5(Ndim,Ndim), HLP6(Ndim,Ndim)
           
           Complex (Kind=Kind(0.d0))  ::  Z
           Integer  ::  I, J, nf, nf_eff, NT, NT1, NTST, NST, N,  N_type
           Real (Kind=Kind(0.d0))  ::  spin,  Mc_step_Weight
           
           If (Symm) Then
              Allocate ( G00_T(Ndim,Ndim,N_FL), G0T_T(Ndim,Ndim,N_FL), GT0_T(Ndim,Ndim,N_FL),  GTT_T(Ndim,Ndim,N_FL) )
           endif

           
           Mc_step_Weight  = 1.d0
           !needed for integration weight in Langevin, dt is taken care of by acceptance ratio for HMC
           if ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" ) Mc_step_weight = Langevin_HMC%get_Delta_t_running()
           
           !Tau = 0
           Do nf_eff = 1, N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              DO J = 1,Ndim 
                 DO I = 1,Ndim
                    Z = cmplx(0.d0, 0.d0, kind(0.D0))
                    if (I == J ) Z = cmplx(1.d0,0.d0,kind(0.d0))
                    G00(I,J,nf) = GR(I,J,nf)
                    GT0(I,J,nf) = GR(I,J,nf)
                    GTT(I,J,nf) = GR(I,J,nf)
                    G0T(I,J,nf) = -(Z - GR(I,J,nf))
                 ENDDO
              ENDDO
           Enddo
           NT = 0
           ! In Module Hamiltonian
           If (Symm) then
              Call Hop_mod_Symm(G00_T,G00,nt)
              Call Hop_mod_Symm(GTT_T,GTT,nt)
              Call Hop_mod_Symm(G0T_T,G0T,nt)
              Call Hop_mod_Symm(GT0_T,GT0,nt)
              !call reconstruction of non-calculated flavor blocks
              If (reconstruction_needed) then
                  Call ham%GR_reconstruction( G00_T )
                  Call ham%GR_reconstruction( GTT_T )
                  Call ham%GRT_reconstruction( GT0_T, G0T_T )
              endif
              CALL ham%OBSERT(NT,  GT0_T,G0T_T,G00_T,GTT_T, PHASE, Mc_step_Weight)
           Else
              !call reconstruction of non-calculated flavor blocks
              If (reconstruction_needed) then
                  Call ham%GR_reconstruction( G00 )
                  Call ham%GR_reconstruction( GTT )
                  Call ham%GRT_reconstruction( GT0, G0T )
              endif
              CALL ham%OBSERT(NT,  GT0,G0T,G00,GTT, PHASE, Mc_step_Weight)
           Endif
           
           ALLOCATE(udvr(N_FL_eff))
           Z = cmplx(1.d0,0.d0,kind(0.d0))
           Do nf_eff = 1, N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              if (Projector) then
                CALL udvr(nf_eff)%init(ndim,'r',WF_R(nf)%P)
              else
                CALL udvr(nf_eff)%init(ndim,'r')
              endif
           enddo

           NST = 1
           DO NT = 0,LTROT - 1
              ! Now wrapup:
              NT1 = NT + 1
              CALL PROPR   (GT0,NT1)
              CALL PROPRM1 (G0T,NT1)
              If  ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
              & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" ) then
                 Call Langevin_HMC%Wrap_Forces(GTT,NT1)
              else
                 CALL PROPRM1 (GTT,NT1)
                 CALL PROPR   (GTT,NT1)
              endif
              ! In Module Hamiltonian
              If (Symm) then
                 Call Hop_mod_Symm(G00_T,G00,0)
                 Call Hop_mod_Symm(GTT_T,GTT,nt1)
                 Call Hop_mod_Symm(G0T_T,G0T,Ltrot,nt1)
                 Call Hop_mod_Symm(GT0_T,GT0,nt1,Ltrot)
                 !call reconstruction of non-calculated flavor blocks
                 If (reconstruction_needed) then
                     Call ham%GR_reconstruction( G00_T )
                     Call ham%GR_reconstruction( GTT_T )
                     Call ham%GRT_reconstruction( GT0_T, G0T_T )
                 endif
                 CALL ham%OBSERT(NT1, GT0_T,G0T_T,G00_T,GTT_T,PHASE, Mc_step_weight)
                 If (( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN"  &
                      & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" ) &
                      &  .and. NT1.ge.LOBS_ST .and. NT1.le.LOBS_EN ) CALL ham%Obser( GTT_T, PHASE, NT1, Mc_step_weight )
              Else
                 !call reconstruction of non-calculated flavor blocks
                 If (reconstruction_needed) then 
                     Call ham%GR_reconstruction( G00 )
                     Call ham%GR_reconstruction( GTT )
                     Call ham%GRT_reconstruction( GT0, G0T )
                 endif
                 CALL ham%OBSERT(NT1, GT0,G0T,G00,GTT,PHASE, Mc_step_weight)
                 If ( ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
                      & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" ) &
                      & .and. NT1.ge.LOBS_ST .and. NT1.le.LOBS_EN ) CALL ham%Obser( GTT, PHASE, NT1, Mc_step_weight )
              Endif
              
              IF ( Stab_nt(NST) == NT1 .AND.  NT1 .NE. LTROT ) THEN
                 !NTST = NT1 - NWRAP
                 !NST  = NT1/(NWRAP)
                 NTST = Stab_nt(NST-1)
                 ! WRITE(6,*) 'NT1, NST: ', NT1,NST
                 CALL WRAPUR(NTST, NT1, udvr)
                 DO nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    HLP4(:,:) = GTT(:,:,nf)
                    HLP5(:,:) = GT0(:,:,nf)
                    HLP6(:,:) = G0T(:,:,nf)
                    Call CGR2_2(GT0(:,:,nf), G00(:,:,nf), GTT(:,:,nf), G0T(:,:,nf), &
                         & udvr(nf_eff), udvst(NST, nf_eff), NDIM)
                    Call Control_Precision_tau(GR(:,:,nf), G00(:,:,nf), Ndim)
                    Call Control_Precision_tau(HLP4      , GTT(:,:,nf), Ndim)
                    Call Control_Precision_tau(HLP5      , GT0(:,:,nf), Ndim)
                    Call Control_Precision_tau(HLP6      , G0T(:,:,nf), Ndim)
                 Enddo
                 NST = NST + 1
              Endif
           ENDDO
           
           DO nf_eff = 1, N_Fl_eff
              CALL udvr(nf_eff)%dealloc
           ENDDO
           DEALLOCATE(udvr)
           If (Symm) Then
              Deallocate ( G00_T, G0T_T, GT0_T,  GTT_T )
           endif

         END SUBROUTINE TAU_M

!--------------------------------------------------------------------

         SUBROUTINE PROPR(AIN,NT) 

           ! Ain =       B(NT-1, NT1) 
           ! Aout= Ain = B(NT  , NT1)

           Implicit none
           Complex (Kind=Kind(0.D0)), intent(INOUT) :: Ain(Ndim,Ndim,N_FL)
           Integer, INTENT(IN) :: NT

           !Locals
           Integer :: nf,nf_eff,n 

           Do nf_eff = 1,N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              Call Hop_mod_mmthr(Ain(:,:,nf),nf,nt)
              Do n = 1,Size(Op_V,1)
!                  X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(Ain(:,:,nf),Op_V(n,nf),nsigma%f(n,nt),'n',nt)
              ENDDO
           Enddo

         end SUBROUTINE PROPR

!--------------------------------------------------------------------

         SUBROUTINE PROPRM1(AIN,NT)

           !Ain = B^{-1}(NT-1, NT1) 
           !Aout= B^{-1}(NT  , NT1)

           Implicit none

           !Arguments 
           Complex (Kind=Kind(0.D0)), intent(Inout) ::  AIN(Ndim, Ndim, N_FL) 
           Integer :: NT, sign = - 1

           ! Locals 
           Integer :: nf,nf_eff, n 

           do nf_eff = 1,N_FL_eff
              nf=Calc_Fl_map(nf_eff)
              !Call MMULT(HLP4,Ain(:,:,nf),Exp_T_M1(:,:,nf) )
              Call Hop_mod_mmthl_m1(Ain(:,:,nf),nf,nt)
              Do n =1,Size(Op_V,1)
                 Call Op_mmultL(Ain(:,:,nf),Op_V(n,nf),nsigma%f(n,nt),'n',nt,sign)
              Enddo
           enddo

         END SUBROUTINE PROPRM1

       end Module Tau_m_mod
