!  Copyright (C) 2016 - 2019 The ALF project
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


      Module Langevin_mod
        
        Use Hamiltonian
        Use UDV_State_mod
        Use Control
        Use Hop_mod

        
        Implicit none
        
        
        Complex (Kind=Kind(0.d0)),  allocatable, private ::  Forces(:,:)



      Contains
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles a  Langevin sweep.
!>   On input GR is on the first time slice and  the storage is full with
!> left propagations.   Udvr  and Udvl are on time slice 1.
!>   On output. The  field configuration is  updated.  GR, Udvr,  Udvl and Udvst are as on input but with the
!> updated configuration.  
!> Equal time measurments as well as time displaced ones  is projector is true are also carried out. 
!> 
!--------------------------------------------------------------------

      SUBROUTINE  Langevin_update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN)
        Implicit none
        
        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: UDVL, UDVR
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:)   :: Test
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR, GR_Tilde
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
        Integer, intent(in) :: LOBS_ST, LOBS_EN


        !Local
        Integer :: NSTM, n, nf, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Real    (Kind=Kind(0.d0)) :: spin
        
        NSTM = Size(Stab_nt,1) - 1 
        !Do  n = 0,NSTM
        !   Write(6,*)  n, Stab_nt(n)
        !Enddo
        
        Forces = cmplx(0.d0,0.d0,Kind(0.d0))
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        Enddo
        NST = 1
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           
           Do nf = 1,N_FL
              CALL HOP_MOD_mmthr   (GR(:,:,nf), nf )
              CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf )
           Enddo
           Do n = 1, size(OP_V,1) 
              Do nf = 1, N_FL
                 spin = nsigma%f(n,ntau1) ! Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
                 N_type = 1
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
                 N_type =  2
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
              enddo
              ! Compute forces here.
              if (OP_V(n,1)%type == 3 ) then
                 Do nf = 1, N_Fl
                    Z = cmplx(0.d0,0.d0,Kind(0.d0))
                    do I = 1,size(OP_V(n,nf)%P,1)
                       do J = 1,size(OP_V(n,nf)%P,1)
                          Z1 =  cmplx(0.d0,0.d0,Kind(0.d0))
                          if ( I == J ) Z1 = cmplx(1.d0,0.d0,Kind(0.d0))
                          Z  = Z +    Op_V(n,nf)%O(I,J) * ( Z1 - Gr(Op_V(n,nf)%P(J),Op_V(n,nf)%P(I), nf) )
                       Enddo
                    Enddo
                    Forces(n,ntau1) = Forces(n,ntau1)  - Op_V(n,nf)%g * Z *  cmplx(real(N_SUN,Kind(0.d0)), 0.d0, Kind(0.d0)) 
                 Enddo
              endif
           enddo

           
           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              CALL WRAPUR(NT1, NTAU1, udvr)
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf = 1, N_FL
                 ! Read from storage left propagation from LTROT to  NTAU1
                 udvl(nf) = udvst(NST, nf)
                 ! Write in storage right prop from 1 to NTAU1
                 udvst(NST, nf) = udvr(nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
                 Z = Z*Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST + 1
           ENDIF
           
           IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
              !Call  Global_tau_mod_Test(Gr,ntau1)
              !Stop
              !write(*,*) "GR before obser sum: ",sum(GR(:,:,1))
              !write(*,*) "Phase before obser : ",phase
              If (Symm) then
                 Call Hop_mod_Symm(GR_Tilde,GR)
                 CALL Obser( GR_Tilde, PHASE, Ntau1 )
              else
                 CALL Obser( GR, PHASE, Ntau1 )
              endif
           endif
        enddo

        Call Control_Langevin(Forces)
        Call Ham_Langevin_update( Forces )
        !DO NTAU1 = 1, LTROT
        !   Do n = 1, size(OP_V,1)
        !      nsigma%f(n,ntau1) = nsigma%flip(n,ntau1)
        !   enddo
        !enddo
        ! Reset Left storage with new fields
        
        Do nf = 1,N_FL
           if (Projector) then
              CALL udvl(nf)%reset('l',WF_L(nf)%P)
              CALL udvst(NSTM, nf)%reset('l',WF_L(nf)%P)
           else
              CALL udvl(nf)%reset('l')
              CALL udvst(NSTM, nf)%reset('l')
           endif
        ENDDO
        
        DO NST = NSTM-1,1,-1
           NT1 = Stab_nt(NST+1)
           NT  = Stab_nt(NST  )
           !Write(6,*)'Hi', NT1,NT, NST
           CALL WRAPUL(NT1, NT, UDVL)
           Do nf = 1,N_FL
              UDVST(NST, nf) = UDVL(nf)
           ENDDO
        ENDDO
        NT1 = stab_nt(1)
        CALL WRAPUL(NT1, 0, UDVL)
        
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        ENDDO
        
        NVAR = 1
        Phase = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf = 1,N_Fl
           CALL CGR(Z, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
           Phase = Phase*Z
        Enddo
        call Op_phase(Phase,OP_V,Nsigma,N_SUN)

        
      END SUBROUTINE LANGEVIN_UPDATE
     

      SUBROUTINE  Langevin_setup
        Implicit none

        Integer :: Nr,Nt
        
        Nr = size(nsigma%f,1)
        Nt = size(nsigma%f,2)
        Allocate ( Forces(Nr,Nt) )
        !Write(6,*) "Langevin: ", Nr,Nt
      end SUBROUTINE Langevin_setup


      SUBROUTINE  Langevin_clear
        Implicit none
        Deallocate ( Forces )
      end SUBROUTINE Langevin_clear

    end Module Langevin_mod
