!  Copyright (C) 2016 - 2018 The ALF project
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
!> defined time-displaced correlations functions. This modules is for the projector code.
!> 
!
!--------------------------------------------------------------------

     Module Tau_p_mod
       Use Hamiltonian
       Use Operator_mod
       Use Control
       Use Hop_mod
       Use UDV_State_mod
       Use tau_m_mod !, only propr, proprm1

     Contains

       SUBROUTINE Tau_p(udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST )

         Implicit none

         Interface
            SUBROUTINE CGRP(PHASE, GRUP, udvr, udvl)
              Use UDV_State_mod
              use MyMats
              CLASS(UDV_State), INTENT(IN) :: udvl, udvr
              COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Intent(INOUT) :: GRUP
              COMPLEX (Kind=Kind(0.d0)), Intent(INOUT) :: PHASE
            End Subroutine CGRP
            SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
              Use Hamiltonian
              Use UDV_Wrap_mod
              Use UDV_State_mod
              Implicit None
              CLASS(UDV_State), intent(inout) :: UDVR(N_FL)
              Integer :: NTAU1, NTAU
            END SUBROUTINE WRAPUR
         End Interface

        ! Storage is full with U^{<} (left)  propagations.


        Integer, Intent(In) :: NSTM, NST
        CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(IN) :: udvl, udvr
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN) :: udvst
        Complex (Kind=Kind(0.d0)), Intent(in) :: GR(NDIM,NDIM,N_FL),  Phase
        Integer, Intent(In) :: STAB_NT(0:NSTM)

!       Local.
        CLASS(UDV_State), Dimension(:), ALLOCATABLE :: udvr_local
        Complex (Kind=Kind(0.d0)) :: DETZ, ZK, DET1(2)
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable  ::  GRUPB, GRUP
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable  ::  G00UP, G0TUP, GT0UP,  GTTUP
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable  ::  G00UP_T, G0TUP_T, GT0UP_T,  GTTUP_T

        Complex (Kind=Kind(0.d0)), allocatable  :: TEMP(:,:), TMPUP(:,:)

        Real    (Kind=kind(0.d0))  :: XMEAN_DYN, XMAX_DYN

        Integer :: NTAUIN,  NTDM,  LFAM, NFAM, N_Part,  LQ , I, NCON, NF, NFLAG, NL, NT1, NT_ST, NT, NTAU, NTAU1,n

        Real (Kind=Kind(0.d0)):: XMEAN, XMAX

        LQ = ndim
        nt_st=0

        do while(STAB_NT(NT_ST)<=THTROT+1)
           nt_st=nt_st+1
        enddo

        ALLOCATE (  GRUPB(LQ,LQ,N_FL), GRUP(LQ,LQ,N_FL), G00UP(LQ,LQ,N_FL), G0TUP(LQ,LQ,N_FL), &
             &      GT0UP(LQ,LQ,N_FL),  GTTUP(LQ,LQ,N_FL), TEMP(LQ,LQ), udvr_local(N_FL) )

        If (Symm) Then
           Allocate ( G00UP_T(LQ,LQ,N_FL), G0TUP_T(LQ,LQ,N_FL), GT0UP_T(LQ,LQ,N_FL),  GTTUP_T(LQ,LQ,N_FL) )
        endif

        do nf=1,N_FL
           call udvr_local(nf)%alloc(ndim,udvr(nf)%N_part)
           udvr_local(nf)=udvr(nf)
        enddo


        Call Wrapur(STAB_NT(NST), STAB_NT(NT_ST), UDVR_local)

        GRUP=GR
        do nt=stab_nt(nst)+1,THTROT+1
           CALL PROPRM1(GRUP,nt)
           CALL PROPR  (GRUP,nt)
        enddo

        GRUPB = GRUP
        do nf=1,N_FL
           do I=1,Ndim
              GRUPB(I,I,nf)=GRUPB(I,I,nf)-1.d0
           enddo
        enddo

        G00UP = GRUP
        GTTUP = GRUP
        GT0UP = GRUP
        G0TUP = GRUPB
        NT1   = 0
        If (Symm) then
           Call Hop_mod_Symm(G00UP_T,G00UP)
           Call Hop_mod_Symm(GTTUP_T,GTTUP)
           Call Hop_mod_Symm(G0TUP_T,G0TUP)
           Call Hop_mod_Symm(GT0UP_T,GT0UP)
           CALL OBSERT (NT1,GT0UP_T,G0TUP_T,G00UP_T,GTTUP_T,PHASE)
        else
           CALL OBSERT (NT1,GT0UP,G0TUP,G00UP,GTTUP,PHASE)
        endif
        ! WRITE(6,*) 'Starting Dyn'
        DO NT = THTROT+1, Ltrot-THTROT
           ! UR is on time slice NT
           NTAU = NT - THTROT -1
           !#ifdef test
           !            WRITE(6,*) 'Ntau: ', NTAU, "on slice",NT
           !#endif
           IF ( NT.EQ.STAB_NT(NT_ST) .and. NTAU/=0) THEN
              !               write(*,*) "stabilization at ",stab_nt(NT_st)
              do nf=1,N_FL
                 !                   write(*,*) "udvl side:",udvst(nt_st,nf)%side
                 CALL CGRP(DetZ, GRUP(:,:,nf), udvr_local(nf), udvst(nt_st,nf))
              enddo
              Call Wrapur(STAB_NT(NT_ST), STAB_NT(NT_ST+1), UDVR_local)
              NT_ST=NT_ST+1

              GRUPB = -GRUP
              do nf=1,N_FL
                 do I=1,Ndim
                    GRUPB(I,I,nf)=GRUPB(I,I,nf)+1.d0
                 enddo
              enddo

              Do nf = 1,N_FL
                 Call Control_Precision_tau(GTTUP(:,:,nf), GRUP(:,:,nf), Ndim)
              enddo

!!$              XMAX      = 0.D0
!!$              XMEAN     = 0.D0
!!$              XMAX_DYN  = 0.D0
!!$              XMEAN_DYN = 0.D0
!!$              do nf=1,N_FL
!!$                 CALL COMPARE (GTTUP(:,:,nf), GRUP(:,:,nf), XMAX,XMEAN)
!!$                 IF (XMAX.GT.XMAX_DYN) XMAX_DYN = XMAX
!!$                 XMEAN_DYN = XMEAN_DYN + XMEAN
!!$              enddo
!!$              !#ifdef test
!!$                                WRITE(6,*) 'Compare up: ',XMEAN/Real(N_FL*LQ*LQ,Kind(0.d0)), XMAX
!!$              !#endif

              GTTUP = GRUP

              do nf=1,N_FL
                 CALL MMULT(TEMP,GRUP(:,:,nf),GT0UP(:,:,nf))
                 GT0UP(:,:,nf) = TEMP
                 CALL MMULT(TEMP,G0TUP(:,:,nf),GRUPB(:,:,nf))
                 G0TUP(:,:,nf) = TEMP
              enddo
           ENDIF                ! Ortho.
           ! Now propagate to Ntau + 1 and call OBSERT.
           NT1 = NT + 1
           !Write(6,*) "CALL PROPR  (GT0UP,NT1)"
           CALL PROPR  (GT0UP,NT1)
           !Write(6,*) "Ret"
           !Write(6,*) " CALL PROPRM1(G0TUP,NT1)"
           CALL PROPRM1(G0TUP,NT1)
           !Write(6,*) "Ret"

           !Write(6,*) "CALL PROPRM1(GTTUP,NT1)"
           CALL PROPRM1(GTTUP,NT1)
           !Write(6,*) "Ret"
           !Write(6,*) "CALL PROPR  (GTTUP,NT1)"
           CALL PROPR  (GTTUP,NT1)
           !Write(6,*) "Ret"


           NTAU1 = NTAU + 1
           If (Symm) then
              Call Hop_mod_Symm(G00UP_T,G00UP)
              Call Hop_mod_Symm(GTTUP_T,GTTUP)
              Call Hop_mod_Symm(G0TUP_T,G0TUP)
              Call Hop_mod_Symm(GT0UP_T,GT0UP)
              Call OBSERT (NTAU1,GT0UP_T,G0TUP_T,G00UP_T,GTTUP_T,PHASE)
           else
              Call OBSERT (NTAU1,GT0UP,G0TUP,G00UP,GTTUP,PHASE)
           endif

        ENDDO

        Do nf=1,N_FL
           call udvr_local(nf)%dealloc
        enddo
        DEALLOCATE (GRUPB, GRUP, G00UP, G0TUP, GT0UP,  GTTUP, TEMP, udvr_local)
        If (Symm) Then
           Deallocate ( G00UP_T, G0TUP_T, GT0UP_T,  GTTUP_T )
        endif

        !         DEALLOCATE (TMPUP, V, D  )
        !         DEALLOCATE (ULR , ULRINV )
        
        RETURN
      END SUBROUTINE Tau_p

    End Module Tau_p_mod
