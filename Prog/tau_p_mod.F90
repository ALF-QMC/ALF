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

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief 
!> This module handles calculation of imaginary-time-displaced Green functions and  
!> calls the routine ObserT.F90 in the Hamiltonian module, so as to compute  user
!> defined time-displaced correlations functions. This module is for the projector code.
!--------------------------------------------------------------------

     Module Tau_p_mod
       Use Hamiltonian_main
       Use Operator_mod
       Use Control
       Use Hop_mod
       Use UDV_State_mod
       Use Langevin_HMC_mod
       Use tau_m_mod  !, only propr, proprm1
       use cgr1_mod, only: cgrp
       use wrapur_mod

     Contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief      This routine computes the time displaced  zero termperature Green functions and call  obserT. 
!> On input:   a) GR,  the equal time Green function,  as  well as udvl, udvr are on time slice 
!>                nt_in= stab_nt(nst)   with   stab_nt(NST) <= THTROT+1  and  stab_nt( NST +1 )  > THTROT+1.
!>             b) The storage, udvst, is full with left propagations from  Ltrot to    stab_nt( NST +1 ).                    
!> On_input    a) GR,  the equal time Green function,  as  well as udvl, udvr are on time slice 
!>                nt_in = = Stab_nt(NST_IN)  with nt_in <= THTROT+1.
!>             b) The storage is full with left propagations for all n's with stab_nt(n) > nt_in.
!>             c) udvl and udvr are on time slice nt_in  such that  a call to CGR with udvl and udvr will
!>                produce Gr.
!>
!> On_output   a) The time displaced Green functions have been computed and measurements carried out.
!>             b) If Langevin then 1)  nt_in = 0, 2) forces are computed  3)  time displaced and equal time
!>                observables are  measured.      
!--------------------------------------------------------------------

       SUBROUTINE Tau_p(udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST_IN, LOBS_ST, LOBS_EN )

         Implicit none

        ! Storage is full with U^{<} (left)  propagations.


        Integer, Intent(In) :: NSTM, NST_IN
        CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(IN) :: udvl, udvr
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN) :: udvst
        Complex (Kind=Kind(0.d0)), Intent(in) :: GR(NDIM,NDIM,N_FL),  Phase
        Integer, Intent(In) :: STAB_NT(0:NSTM)
        Integer, Intent(In) :: LOBS_ST, LOBS_EN
        
!       Local.
        CLASS(UDV_State), Dimension(:), ALLOCATABLE :: udvr_local
        Complex (Kind=Kind(0.d0)) :: DETZ, ZK, DET1(2)
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable  ::  GRUPB, GRUP
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable  ::  G00UP, G0TUP, GT0UP,  GTTUP
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable  ::  G00UP_T, G0TUP_T, GT0UP_T,  GTTUP_T
        Complex (Kind=Kind(0.d0)), allocatable  :: TEMP(:,:), TMPUP(:,:)

        Real    (Kind=kind(0.d0))  :: XMEAN_DYN, XMAX_DYN

        Integer :: NTAUIN,  NTDM,  LFAM, NFAM, N_Part,  LQ , I, NCON, NF, nf_eff, NFLAG, NL, NT1, NT_ST, NT, NTAU, NTAU1,n

        Real (Kind=Kind(0.d0)) :: XMEAN, XMAX
        Real (Kind=Kind(0.d0)) :: Mc_step_weight
        
        LQ = ndim

        Mc_step_weight = 1.d0
        if ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" ) Mc_step_weight =  Langevin_HMC%get_Delta_t_running()

        
        ALLOCATE (  GRUPB(LQ,LQ,N_FL), GRUP(LQ,LQ,N_FL), G00UP(LQ,LQ,N_FL), G0TUP(LQ,LQ,N_FL), &
             &      GT0UP(LQ,LQ,N_FL),  GTTUP(LQ,LQ,N_FL), TEMP(LQ,LQ), udvr_local(N_FL_eff) )

        If (Symm) Then
           Allocate ( G00UP_T(LQ,LQ,N_FL), G0TUP_T(LQ,LQ,N_FL), GT0UP_T(LQ,LQ,N_FL),  GTTUP_T(LQ,LQ,N_FL) )
        endif

        do nf_eff=1,N_FL_eff
           call udvr_local(nf_eff)%alloc(ndim,udvr(nf_eff)%N_part)
           udvr_local(nf_eff)=udvr(nf_eff)
        enddo

        GTTUP = GR ! On time slice Stab_nt(NST_IN)
        NT_ST = NST_IN
        do NT = Stab_nt(NT_ST)+1, Thtrot + 1
           If  ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
              & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" ) then            
              Call Langevin_HMC%Wrap_Forces(GTTUP,NT)
           else
              CALL PROPRM1 (GTTUP,NT)
              CALL PROPR   (GTTUP,NT)
           endif
           If (( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN"  &
              & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" ) &
              .and. NT .ge. LOBS_ST .and. NT .le. LOBS_EN ) then
              If (Symm) then
                 Call Hop_mod_Symm(GTTUP_T,GTTUP,nt)
                 !call reconstruction of non-calculated flavor blocks
                 If (reconstruction_needed) Call ham%GR_reconstruction( GTTUP_T )
                 CALL ham%Obser( GTTUP_T, PHASE, NT, Mc_step_weight )
              else
               !call reconstruction of non-calculated flavor blocks
                 If (reconstruction_needed) Call ham%GR_reconstruction( GTTUP )
                 CALL ham%Obser( GTTUP, PHASE, NT, Mc_step_weight )
              endif
           endif
           IF ( NT .EQ. STAB_NT(NT_ST+1) ) THEN
              Call Wrapur(STAB_NT(NT_ST), STAB_NT(NT_ST+1), UDVR_local)
              do nf_eff=1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 CALL CGRP(DetZ, GRUP(:,:,nf), udvr_local(nf_eff), udvst(nt_st + 1,nf_eff))
              enddo
              Do nf_eff = 1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 Call Control_Precision_tau(GTTUP(:,:,nf), GRUP(:,:,nf), Ndim)
              enddo
              GTTUP = GRUP
              NT_ST = NT_ST+1
           endif
        enddo

        GRUPB = GTTUP
        do nf_eff=1,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           do I=1,Ndim
              GRUPB(I,I,nf)=GRUPB(I,I,nf)-1.d0
           enddo
        enddo

        G00UP = GTTUP
        !GTTUP = GTTUP
        GT0UP = GTTUP
        G0TUP = GRUPB
        NTAU   = 0
        If (Symm) then
           Call Hop_mod_Symm(G00UP_T,G00UP,Thtrot + 1)
           Call Hop_mod_Symm(GTTUP_T,GTTUP,Thtrot + 1)
           Call Hop_mod_Symm(G0TUP_T,G0TUP,Thtrot + 1)
           Call Hop_mod_Symm(GT0UP_T,GT0UP,Thtrot + 1)
           !call reconstruction of non-calculated flavor blocks
           If (reconstruction_needed) then
               Call ham%GR_reconstruction( G00UP_T )
               Call ham%GR_reconstruction( GTTUP_T )
               Call ham%GRT_reconstruction( GT0UP_T, G0TUP_T )
           endif
           CALL ham%OBSERT (NTAU,GT0UP_T,G0TUP_T,G00UP_T,GTTUP_T,PHASE, Mc_step_Weight)
        else
           !call reconstruction of non-calculated flavor blocks
           If (reconstruction_needed) then
              Call ham%GR_reconstruction( G00UP )
              Call ham%GR_reconstruction( GTTUP )
              Call ham%GRT_reconstruction( GT0UP, G0TUP )
           endif
           CALL ham%OBSERT (NTAU,GT0UP,G0TUP,G00UP,GTTUP,PHASE, Mc_step_Weight)
        endif
        DO NT = THTROT+1, Ltrot-THTROT
           ! UR is on time slice NT
           NTAU = NT - THTROT -1
           IF ( NT .EQ. STAB_NT(NT_ST+1) .and. NTAU /= 0) THEN
              Call Wrapur(STAB_NT(NT_ST), STAB_NT(NT_ST+1), UDVR_local)
              do nf_eff=1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 CALL CGRP(DetZ, GRUP(:,:,nf), udvr_local(nf_eff), udvst(nt_st + 1,nf_eff))
              enddo
              NT_ST = NT_ST+1
              Do nf_eff = 1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 Call Control_Precision_tau(GTTUP(:,:,nf), GRUP(:,:,nf), Ndim)
              enddo
              GTTUP = GRUP

              GRUPB = -GRUP
              do nf_eff=1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 do I=1,Ndim
                    GRUPB(I,I,nf)=GRUPB(I,I,nf)+1.d0
                 enddo
              enddo
              do nf_eff=1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 CALL MMULT(TEMP,GRUP(:,:,nf),GT0UP(:,:,nf))
                 GT0UP(:,:,nf) = TEMP
                 CALL MMULT(TEMP,G0TUP(:,:,nf),GRUPB(:,:,nf))
                 G0TUP(:,:,nf) = TEMP
              enddo
           ENDIF
           NT1 = NT + 1
           CALL PROPR  (GT0UP,NT1)
           CALL PROPRM1(G0TUP,NT1)
           If  ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
            & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") then
              Call Langevin_HMC%Wrap_Forces(GTTUP,NT1)
           else
              CALL PROPRM1 (GTTUP,NT1)
              CALL PROPR   (GTTUP,NT1)
           endif

           NTAU1 = NTAU + 1
           If (Symm) then
              Call Hop_mod_Symm(G00UP_T,G00UP,Thtrot + 1)
              Call Hop_mod_Symm(GTTUP_T,GTTUP,nt1)
              Call Hop_mod_Symm(G0TUP_T,G0TUP,Thtrot + 1,nt1)
              Call Hop_mod_Symm(GT0UP_T,GT0UP,nt1,Thtrot + 1)
              !call reconstruction of non-calculated flavor blocks
              If (reconstruction_needed) then
                  Call ham%GR_reconstruction( G00UP_T )
                  Call ham%GR_reconstruction( GTTUP_T )
                  Call ham%GRT_reconstruction( GT0UP_T, G0TUP_T )
              endif
              Call ham%OBSERT (NTAU1,GT0UP_T,G0TUP_T,G00UP_T,GTTUP_T,PHASE,Mc_step_weight)
              If (( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
                   & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" )&
                   &.and. NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN ) CALL ham%Obser( GTTUP_T, PHASE, NT1, Mc_step_weight )
           else
              !call reconstruction of non-calculated flavor blocks
              If (reconstruction_needed) then
                  Call ham%GR_reconstruction( G00UP )
                  Call ham%GR_reconstruction( GTTUP )
                  Call ham%GRT_reconstruction( GT0UP, G0TUP )
              endif
              Call ham%OBSERT (NTAU1,GT0UP,G0TUP,G00UP,GTTUP,PHASE,Mc_step_weight)
              If ((str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
                   & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") &
                   & .and. NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN ) CALL ham%Obser( GTTUP, PHASE, NT1, Mc_step_weight )
           endif

        ENDDO

        If (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" &
           & .or. str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") then   ! Finish calculating the forces
           DO NT = Ltrot-THTROT + 1, Ltrot - 1
              ! UR is on time slice NT
              IF ( NT .EQ. STAB_NT(NT_ST+1) ) THEN
                 Call Wrapur(STAB_NT(NT_ST), STAB_NT(NT_ST+1), UDVR_local)
                 do nf_eff=1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    CALL CGRP(DetZ, GRUP(:,:,nf), udvr_local(nf_eff), udvst(nt_st + 1,nf_eff))
                 enddo
                 NT_ST = NT_ST+1
                 Do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    Call Control_Precision_tau(GTTUP(:,:,nf), GRUP(:,:,nf), Ndim)
                 enddo
                 GTTUP = GRUP
              ENDIF
              NT1 = NT + 1
              !TODO reconstruction required?
              Call Langevin_HMC%Wrap_Forces(GTTUP,NT1)
              If (NT1 .ge. LOBS_ST .and. NT1 .le. LOBS_EN ) Then
                 If (Symm) then
                    Call Hop_mod_Symm(GTTUP_T,GTTUP,nt1)
                    !call reconstruction of non-calculated flavor blocks
                    If (reconstruction_needed) Call ham%GR_reconstruction( GTTUP_T )
                    CALL ham%Obser( GTTUP_T, PHASE, NT1, Mc_step_weight )
                 else
                    !call reconstruction of non-calculated flavor blocks
                    If (reconstruction_needed) Call ham%GR_reconstruction( GTTUP )
                    CALL ham%Obser( GTTUP, PHASE, NT1, Mc_step_weight )
                 endif
              endif
           Enddo
        endif
        
        Do nf_eff=1,N_FL_eff
           call udvr_local(nf_eff)%dealloc
        enddo
        DEALLOCATE (GRUPB, GRUP, G00UP, G0TUP, GT0UP,  GTTUP, TEMP, udvr_local)
        If (Symm) Then
           Deallocate ( G00UP_T, G0TUP_T, GT0UP_T,  GTTUP_T )
        endif
        
      END SUBROUTINE Tau_p
      
    End Module Tau_p_mod
