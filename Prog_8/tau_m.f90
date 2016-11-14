     Module Tau_m_mod

       Use Hamiltonian 
       Use Operator_mod
       Use Precdef
       Use Control
       Use Hop_mod

       Contains

         SUBROUTINE TAU_M( UST,DST,VST, GR, PHASE, NSTM, NWRAP, STAB_NT  ) 
           
           Implicit none

           Interface
              SUBROUTINE WRAPUL(NTAU1, NTAU, UL ,DL, VL)
                Use Hamiltonian
                Implicit none
                COMPLEX (KIND=8) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
                COMPLEX (KIND=8) :: DL(Ndim,N_FL)
                Integer :: NTAU1, NTAU
              END SUBROUTINE WRAPUL
              SUBROUTINE CGR2_2(GRT0, GR00, GRTT, GR0T, U2, D2, V2, U1, D1, V1, LQ)
                Use Precdef
                Use MyMats
                Use UDV_WRAP_mod
                Implicit none
                
                !  Arguments
                Integer,  intent(in) :: LQ
                Complex (Kind=double), intent(in)    :: U1(LQ,LQ), V1(LQ,LQ), U2(LQ,LQ), V2(LQ,LQ)
                Complex (Kind=double), intent(in)    :: D2(LQ), D1(LQ)
                Complex (Kind=double), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)
              end SUBROUTINE CGR2_2
              SUBROUTINE CGR2_1(GRT0, GR00, GRTT, GR0T, U2, D2, V2, U1, D1, V1, LQ, NVAR)
                Use Precdef
                Use MyMats
                USe UDV_Wrap_mod
                Implicit none
                !  Arguments
                Integer,  intent(in) :: LQ, NVAR
                Complex (Kind=double), intent(in)    :: U1(LQ,LQ), V1(LQ,LQ), U2(LQ,LQ), V2(LQ,LQ)
                Complex (Kind=double), intent(in)    :: D2(LQ), D1(LQ)
                Complex (Kind=double), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)
              end SUBROUTINE CGR2_1
              SUBROUTINE CGR2(GRT0, GR00, GRTT, GR0T, U2, D2, V2, U1, D1, V1, LQ)
          
                !       B2 = U2*D2*V2
                !       B1 = V1*D1*U1
                !Calc:      (  1   B1 )^-1   i.e. 2*LQ \times 2*LQ matrix
                !           (-B2   1  )
                
                
                Use Precdef
                Use UDV_WRAP_mod
                Use MyMats
                
                Implicit none
                
                !  Arguments
                Integer :: LQ
                Complex (Kind=double), intent(in)    :: U1(LQ,LQ), V1(LQ,LQ), U2(LQ,LQ), V2(LQ,LQ)
                Complex (Kind=double), intent(in)    :: D2(LQ), D1(LQ)
                Complex (Kind=double), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)
              end SUBROUTINE CGR2
           end Interface
     
           Integer, Intent(In) :: NSTM, NWRAP
           Complex (Kind=double), Intent(in) :: UST(NDIM,NDIM,NSTM,N_FL), VST(NDIM,NDIM,NSTM,N_FL), DST(NDIM,NSTM,N_FL) 
           Complex (Kind=double), Intent(in) :: GR(NDIM,NDIM,N_FL),  Phase
           Integer, Intent(In) :: STAB_NT(0:NSTM)         

           

           ! Local 
           ! This could be placed as  private for the module 
           Complex (Kind=double)  :: GT0(NDIM,NDIM,N_FL),  G00(NDIM,NDIM,N_FL), GTT(NDIM,NDIM,N_FL), G0T(NDIM,NDIM,N_FL)
           Complex (Kind=double)  :: UL(Ndim,Ndim,N_FL), DL(Ndim,N_FL), VL(Ndim,Ndim,N_FL) 
           Complex (Kind=double)  :: UR(Ndim,Ndim,N_FL), DR(Ndim,N_FL), VR(Ndim,Ndim,N_FL) 
           Complex (Kind=double)  :: HLP4(Ndim,Ndim), HLP5(Ndim,Ndim), HLP6(Ndim,Ndim)
           
           Complex (Kind=double)  ::  Z
           Integer  ::  I, J, nf, NT, NT1, NTST, NST, NVAR
           
           !Tau = 0
           Do nf = 1, N_FL
              DO J = 1,Ndim 
                 DO I = 1,Ndim
                    Z = cmplx(0.d0, 0.d0, kind(0.D0))
                    if (I == J ) Z = cone
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
           
           Do nf = 1, N_FL
              CALL INITD(UR(:,:,nf),cone)
              CALL INITD(VR(:,:,nf),cone)
           enddo
           DR = cone
              
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
                 CALL WRAPUR(NTST, NT1,UR, DR, VR)
                 DO nf = 1,N_FL
                    UL(:,:,nf) = UST(:,:,NST,nf)
                    VL(:,:,nf) = VST(:,:,NST,nf)
                    DL(:  ,nf) = DST(:  ,NST,nf)
                 Enddo
                 Do nf = 1,N_FL
                    Do J = 1,Ndim
                       DO I = 1,Ndim
                          HLP4(I,J) = GTT(I,J,nf)
                          HLP5(I,J) = GT0(I,J,nf)
                          HLP6(I,J) = G0T(I,J,nf)
                       Enddo
                    Enddo
                    NVAR = 1
                    IF (NT1  >  LTROT/2) NVAR = 2
                    !DO I = 1,Ndim
                    !   Write(6,*) DL(I,nf)*DR(I,nf)
                    !enddo
                    !Write(6,*) 'Call CGR2'
                    Call CGR2_2(GT0(:,:,nf), G00(:,:,nf), GTT(:,:,nf), G0T(:,:,nf), &
                         &      UR(:,:,nf),DR(:,nf),VR(:,:,nf), UL(:,:,nf),DL(:,nf),VL(:,:,nf),NDIM)
                    !Call CGR2   (GT0(:,:,nf), G00(:,:,nf), GTT(:,:,nf), G0T(:,:,nf), &
                    !     &      UR(:,:,nf),DR(:,nf),VR(:,:,nf), UL(:,:,nf),DL(:,nf),VL(:,:,nf),NDIM)
                    
                    !Call CGR2_1(GT0(:,:,nf), G00(:,:,nf), GTT(:,:,nf), G0T(:,:,nf), &
                    !      &    UR(:,:,nf),DR(:,nf),VR(:,:,nf), UL(:,:,nf),DL(:,nf),VL(:,:,nf),NDIM,NVAR)
                    
                    !Write(6,*) 'End Call CGR2'
                    !Write(6,*) ' Tau ', NT1
                    !Write(6,*) ' G00 '
                    Call Control_Precision_tau(GR(:,:,nf), G00(:,:,nf), Ndim)
                    !Write(6,*) ' GTT '
                    Call Control_Precision_tau(HLP4      , GTT(:,:,nf), Ndim)
                    !Write(6,*) ' GT0 '
                    Call Control_Precision_tau(HLP5      , GT0(:,:,nf), Ndim)
                    !Write(6,*) ' G0T '
                    Call Control_Precision_tau(HLP6      , G0T(:,:,nf), Ndim)
                 Enddo
                 NST = NST + 1
              Endif
           ENDDO
           
         END SUBROUTINE TAU_M
            
!==============================================================
         
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
              !CALL MMULT(HLP4,Exp_T(:,:,nf) ,Ain(:,:,nf))
              Call Hop_mod_mmthr(Ain(:,:,nf),HLP4,nf)
              Do n = 1,Size(Op_V,1)
                 X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(HLP4,Op_V(n,nf),X,Ndim)
              ENDDO
              Call ZLACPY('A', Ndim, Ndim, HLP4, Ndim, Ain(1,1, nf), Ndim)
           Enddo
           Deallocate(HLP4)
           
         end SUBROUTINE PROPR
!==============================================================
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
!==============================================================
       end Module Tau_m_mod
