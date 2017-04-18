!  Copyright (C) 2016, 2017 The ALF project
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

     Module Tau_m_mod
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This routine  handles calculation of imagimary time displaced Green functions and  
!> calls the routine ObserT.f90 in the Hamiltonian module,  so as to compute the user
!> defined  correlations. 
!
!--------------------------------------------------------------------
       Use Hamiltonian 
       Use Operator_mod
       Use Control
       Use Hop_mod
       Use UDV_State_mod

       Contains

         SUBROUTINE TAU_M(udvst, GR, PHASE, NSTM, NWRAP, STAB_NT  ) 
           
           Implicit none

           Interface
              SUBROUTINE WRAPUR(NTAU1, NTAU, udvr)
                Use Hamiltonian
                Use UDV_State_mod
                Implicit none
                CLASS(UDV_State), intent(inout) :: UDVr(N_FL)
                Integer :: NTAU1, NTAU
              END SUBROUTINE WRAPUR
              SUBROUTINE CGR2_2(GRT0, GR00, GRTT, GR0T, udv2, udv1, LQ)
                Use MyMats
                Use UDV_WRAP_mod
                Use UDV_State_mod
                Implicit none
                !  Arguments
                Integer,  intent(in) :: LQ
                CLASS(UDV_State), intent(in) :: udv1, udv2
                Complex (Kind=Kind(0.d0)), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)
              end SUBROUTINE CGR2_2
           end Interface
     
           Integer, Intent(In) :: NSTM, NWRAP
           CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN) :: udvst
           Complex (Kind=Kind(0.d0)), Intent(in) :: GR(NDIM,NDIM,N_FL),  Phase
           Integer, Intent(In) :: STAB_NT(0:NSTM)         

           

           ! Local 
           ! This could be placed as  private for the module 
           Complex (Kind=Kind(0.d0))  :: GT0(NDIM,NDIM,N_FL),  G00(NDIM,NDIM,N_FL), GTT(NDIM,NDIM,N_FL), G0T(NDIM,NDIM,N_FL)
           CLASS(UDV_State), DIMENSION(:), ALLOCATABLE :: udvr
           Complex (Kind=Kind(0.d0))  :: HLP4(Ndim,Ndim), HLP5(Ndim,Ndim), HLP6(Ndim,Ndim)
           
           Complex (Kind=Kind(0.d0))  ::  Z
           Integer  ::  I, J, nf, NT, NT1, NTST, NST, NVAR
           
           !Tau = 0
           Do nf = 1, N_FL
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
           CALL OBSERT(NT,  GT0,G0T,G00,GTT, PHASE)
           
           ALLOCATE(udvr(N_FL))
           Z = cmplx(1.d0,0.d0,kind(0.d0))
           Do nf = 1, N_FL
              CALL udvr(nf)%init(ndim)
           enddo
              
           NST = 1
           DO NT = 0,LTROT - 1
              ! Now wrapup:
              NT1 = NT + 1
              CALL PROPR   (GT0,NT1)
              CALL PROPRM1 (G0T,NT1)
              CALL PROPRM1 (GTT,NT1)
              CALL PROPR   (GTT,NT1)
              ! In Module Hamiltonian
              CALL OBSERT(NT1, GT0,G0T,G00,GTT,PHASE)
              
              IF ( Stab_nt(NST) == NT1 .AND.  NT1 .NE. LTROT ) THEN
                 !NTST = NT1 - NWRAP
                 !NST  = NT1/(NWRAP)
                 NTST = Stab_nt(NST-1)
                 ! WRITE(6,*) 'NT1, NST: ', NT1,NST
                 CALL WRAPUR(NTST, NT1, udvr)
                 DO nf = 1,N_FL
                    HLP4(:,:) = GTT(:,:,nf)
                    HLP5(:,:) = GT0(:,:,nf)
                    HLP6(:,:) = G0T(:,:,nf)
                    NVAR = 1
                    IF (NT1  >  LTROT/2) NVAR = 2
                    Call CGR2_2(GT0(:,:,nf), G00(:,:,nf), GTT(:,:,nf), G0T(:,:,nf), &
                         & udvr(nf), udvst(NST, nf), NDIM)
                    Call Control_Precision_tau(GR(:,:,nf), G00(:,:,nf), Ndim)
                    Call Control_Precision_tau(HLP4      , GTT(:,:,nf), Ndim)
                    Call Control_Precision_tau(HLP5      , GT0(:,:,nf), Ndim)
                    Call Control_Precision_tau(HLP6      , G0T(:,:,nf), Ndim)
                 Enddo
                 NST = NST + 1
              Endif
           ENDDO
           
           DO nf = 1, N_Fl
              CALL udvr(nf)%dealloc
           ENDDO
           DEALLOCATE(udvr)
         END SUBROUTINE TAU_M

!--------------------------------------------------------------------
         
         SUBROUTINE PROPR(AIN,NT) 

           ! Ain =       B(NT-1, NT1) 
           ! Aout= Ain = B(NT  , NT1)
           
           Implicit none
           Complex (Kind=Kind(0.D0)), intent(INOUT) :: Ain(Ndim,Ndim,N_FL)
           Integer, INTENT(IN) :: NT

           !Locals
           Integer :: J,I,nf,n 
           Complex (Kind=Kind(0.D0)), allocatable, Dimension(:, :) :: HLP4
           Real    (Kind=Kind(0.D0)) :: X
           
           Allocate(HLP4(Ndim, Ndim))

           Do nf = 1,N_FL
              Call Hop_mod_mmthr(Ain(:,:,nf),HLP4,nf)
              Do n = 1,Size(Op_V,1)
                 X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(HLP4,Op_V(n,nf),X,Ndim)
              ENDDO
              Call ZLACPY('A', Ndim, Ndim, HLP4, Ndim, Ain(1,1, nf), Ndim)
           Enddo
           Deallocate(HLP4)
           
         end SUBROUTINE PROPR

!--------------------------------------------------------------------

         SUBROUTINE PROPRM1(AIN,NT)

           !Ain = B^{-1}(NT-1, NT1) 
           !Aout= B^{-1}(NT  , NT1)
           
           
           Implicit none
           
           !Arguments 
           Complex (Kind=Kind(0.D0)), intent(Inout) ::  AIN(Ndim, Ndim, N_FL) 
           Integer :: NT
           
           ! Locals 
           Integer :: J,I,nf,n 
           Complex (Kind=Kind(0.D0)), allocatable, Dimension(:, :) :: HLP4
           Real    (Kind=Kind(0.D0)) :: X
           
           Allocate(HLP4(Ndim, Ndim))
           
           do nf = 1,N_FL
              !Call MMULT(HLP4,Ain(:,:,nf),Exp_T_M1(:,:,nf) )
              Call Hop_mod_mmthl_m1(Ain(:,:,nf),HLP4,nf)
              Do n =1,Size(Op_V,1)
                 X = -Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultL(HLP4,Op_V(n,nf),X,Ndim)
              Enddo
              Call ZLACPY('A', Ndim, Ndim, HLP4, Ndim, Ain(1,1, nf), Ndim)
           enddo
           Deallocate(HLP4)
           
         END SUBROUTINE PROPRM1

       end Module Tau_m_mod
