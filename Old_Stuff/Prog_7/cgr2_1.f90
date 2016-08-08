      SUBROUTINE CGR2_1(GRT0, GR00, GRTT, GR0T, U2, D2, V2, U1, D1, V1, LQ, NVAR)

        !       B2 = U2*D2*V2 is right (i.e. from time slice 0 to tau) propagation to time tau 
        !       B1 = V1*D1*U1 is left  (i.e. from time slice Ltrot to tau) propagation to time tau
        !Calc:      (  1   B1 )^-1       ( G00    G0T )   
        !           (-B2   1  )     ==   ( GT0    GTT )
        !            
        !       G00 = (1 + B1*B2)^-1   G0T = -(1 -  G00 )*B2^-1
        !       GT0  =   B2 * G00      GTT =  (1 + B2*B1)^-1
 
        !  Here you want to compute  G00, G0T, GT0 and GTT just by involving LQ x LQ matrix operations. 
        !  If NVAR == 1  then the  large scales are in D1 
        !  If NVAR == 2  then the  large scales are in D2
        Use Precdef
        Use MyMats
        USe UDV_Wrap_mod

        Implicit none
        
        Interface
           SUBROUTINE CGR(Z,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)
             COMPLEX(Kind=8), Dimension(:,:), Intent(In) ::  URUP, VRUP, ULUP, VLUP
             COMPLEX(Kind=8), Dimension(:), Intent(In)   ::  DLUP, DRUP
             COMPLEX(Kind=8), Dimension(:,:), Intent(INOUT) :: GRUP
        
             COMPLEX(Kind=8) :: Z
           END SUBROUTINE CGR
        end Interface


        !  Arguments
        Integer,  intent(in) :: LQ, NVAR
        Complex (Kind=double), intent(in)    :: U1(LQ,LQ), V1(LQ,LQ), U2(LQ,LQ), V2(LQ,LQ)
        Complex (Kind=double), intent(in)    :: D2(LQ), D1(LQ)
        Complex (Kind=double), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)


        ! Local::
        Complex  (Kind=double) :: HLP1(LQ,LQ), HLP2(LQ,LQ), U(LQ,LQ), D(LQ), V(LQ,LQ)
        Complex  (Kind=double) :: Z, Z1, Z2
        Real     (Kind=double) :: X, Xmax, Xmin, X1, X2, Xmax1, Xmax2, Xmean
        Integer                :: I, J, M, NCON, NVAR1

        Complex  (Kind=double) :: V2inv(LQ,LQ), V1inv(LQ,LQ)
        

        NCON = 0

        Call INV( V2, V2inv, Z2)
        CALL INV( V1, V1inv, Z1)

        
        DO I = 1,LQ
           DO J = 1,LQ
              HLP1(I,J) = CONJG( U1(J,I) ) 
           ENDDO
        ENDDO
        CALL MMULT(HLP2,HLP1,U1)
        HLP1 = cmplx(0.d0,0.d0,kind=8)
        DO I = 1,LQ
           HLP1(I,I) =  cmplx(1.d0,0.d0,kind=8)
        ENDDO
        Xmax = 0.d0
        CALL COMPARE(HLP1, HLP2, XMAX, XMEAN)
        
        DO I = 1,LQ
           DO J = 1,LQ
              HLP1(I,J) = CONJG( U2(J,I) ) 
           ENDDO
        ENDDO
        CALL MMULT(HLP2,HLP1,U2)
        HLP1 = cmplx(0.d0,0.d0,kind=8)
        DO I = 1,LQ
           HLP1(I,I) =  cmplx(1.d0,0.d0,kind=8)
        ENDDO
        Xmax1 = 0.d0
        CALL COMPARE(HLP1, HLP2, XMAX1, XMEAN)
        Write(77,*) "Cgr2_1  V1inv V2inv : ", Xmax, Xmax1 

!!$        Xmax = 0.d0
!!$        do  I = 1,LQ
!!$           do j = 1,LQ
!!$              X = sqrt(dble(V1(i,j)*conjg(V1(i,j))))
!!$              if (X > Xmax) Xmax = X
!!$           enddo
!!$        enddo
!!$        Write(77,*) 'In cgr2_1 Xmax V1: ', Xmax, Z2
!!$        do  I = 1,LQ
!!$           do j = 1,LQ
!!$              X = sqrt(dble(V2(i,j)*conjg(V2(i,j))))
!!$              if (X > Xmax) Xmax = X
!!$           enddo
!!$        enddo
!!$        Write(77,*) 'In cgr2_1 Xmax V2: ', Xmax, Z1

        ! Compute G00
        ! G00 = (1 + B1*B2)^-1 = (1 + V1 D1 U1 U2 D2 V2 )^-1 = 
        !                      = ( V1  (   V1^-1 V2^-1 +  D1 U1 U2 D2 ) V2 )^-1 =
        !                      = V2^-1 ( (V2 V1)^-1 + D1 U1 U2 D2     )^-1 V1^-1
        Call MMULT(HLP1,V1inv,V2inv)
        Call MMULT(HLP2,U1,U2)
        DO J = 1,LQ
           DO I = 1,LQ
              HLP2(I,J) = D1(I)*HLP2(I,J)*D2(J)  + HLP1(I,J)
           ENDDO
        ENDDO
        Xmax1 = dble( D1(1) )
        Xmax2 = dble( D2(1) )
        DO I = 2,LQ
           If (  dble( D1(I) ) > Xmax1 ) Xmax1 =  dble( D1(I) )
           If (  dble( D2(I) ) > Xmax2 ) Xmax2 =  dble( D2(I) )
        Enddo
        Nvar1 = 1
        If ( Xmax2 > Xmax1) Nvar1 = 2
        If (Nvar1 == 1) then
           !  V2^-1 (UDV  )^-1 V1^-1 =  V2^-1 V^-1 D^-1 U^-1 V1^-1
           Call UDV_WRAP(HLP2, U, D, V, Ncon) 
           CALL INV  (V,HLP2 ,Z   )
           CALL MMULT(V,V2inv,HLP2)
           DO J = 1,LQ
              DO I = 1,LQ
                 V(I,J) = V(I,J)/D(J)
              ENDDO
           ENDDO
           DO I = 1,LQ
              DO J = 1,LQ
                 HLP1(I,J) = Conjg(U(J,I))
              ENDDO
           ENDDO
           CALL MMULT( HLP2, HLP1,V1inv)
           CALL MMULT (GR00, V, HLP2)
        else
           !  V2^-1 (UDV  )^(-1,*) V1^-1 =  V2^-1 U  D^-1 V^(-1,*) V1^-1
           DO J = 1,LQ
              DO I = 1,LQ
                 HLP1(I,J) = conjg(HLP2(J,I))
              ENDDO
           ENDDO
           Call UDV_WRAP(HLP1, U, D, V, Ncon) 
           Call MMULT(HLP1, V2inv, U)
           DO J = 1,LQ
              DO I = 1,LQ
                 HLP1(I,J) = HLP1(I,J)/D(J)
              ENDDO
           ENDDO
           CALL INV (V, HLP2, Z)
           DO J = 1,LQ
              DO I = 1,LQ
                 V(I,J) = CONJG(HLP2(J,I))
              ENDDO
           ENDDO
           CALL MMULT(HLP2,V,V1inv)
           CALL MMULT(GR00,HLP1,HLP2)
        endif
        
        ! Compute G0T
        ! G00 = (1 + B1*B2)^-1   G0T = -(1 -  (1 + B1*B2)^-1 )*B2^-1              =
        !                            = -( 1 +  B1*B2 - 1)  (1 + B1*B2)^-1 * B2^-1 =
        !                            = - B1 * B2 (1 + B1*B2)^-1 * B2^-1           = 
        !                            = -( B2 ( 1+ B1*B2) * B2^-1 B1^-1)^-1        =
        !                            = -( B2 ( 1+ B1*B2) * (B1 B2)^-1)^-1         = 
        !                            = -( B2 ( (B1 B2)^-1  + 1 ) )^-1             =
        !                            = -( B1^-1   +   B2)^-1                      =
        !                       -G0T*=  ( B1*^-1  +  B2*)^-1                      =
        !                            =  ( V1*^-1 D1*^-1 U1  + V2* D2* U2*)^-1     =
        !                            =  ( V1*^-1 (  D1*^-1 U1 U2 + V1* V2* D2* ) U2* )^-1 =
        !                            =  U2 ( D1*^-1 (U1 U2) + ( V2 V1)* D2* )^-1 V1*
        !                            =  U2 ( D1*^-1 (U1 U2) + ( V2 V1)* D2* )^-1 V1*
        !       B2 = U2*D2*V2
        !       B1 = V1*D1*U1
        Xmax2 = dble(cmplx(1.d0,0.d0)/D1(1))
        Xmax1 = dble(D2(1))
        Do I = 2,LQ
           X2 = dble(cmplx(1.d0,0.d0)/D1(I))
           X1 = dble(D2(I))
           If ( X2 > Xmax2 ) Xmax2 = X2
           If ( X1 > Xmax1 ) Xmax1 = X1
        ENDDO
        NVAR1 = 1
        If (Xmax2  > Xmax1)  Nvar1 = 2
        Call  MMULT(HLP1,U1,U2)
        DO J = 1,LQ
           DO I =1,LQ
              HLP1(I,J) =  HLP1(I,J)/conjg(D1(I))
           ENDDO
        ENDDO
        Call MMULT(V,V2,V1)
        DO J = 1,LQ
           DO I = 1,LQ
              HLP2(I,J) = Conjg(V(J,I))
           ENDDO
        ENDDO
        DO J = 1,LQ
           DO I =1,LQ
              HLP2(I,J) = HLP1(I,J) +  HLP2(I,J)*conjg(D2(J))
           ENDDO
        ENDDO
        NCON = 0
        IF ( NVAR1 == 1 ) Then
           !  UDV of HLP2
           !  -G0T*= U2 V^-1 D^-1 U* V1*
           CALL UDV_WRAP(HLP2,U,D,V,NCON)
           CALL MMULT (HLP1, V1, U)
           DO I = 1,LQ
              DO J = 1,LQ
                 U(I,J) = conjg(HLP1(J,I))
              ENDDO
           ENDDO
           CALL INV(V,HLP2,Z)
           Call MMULT(HLP1,U2,HLP2)
           DO J = 1,LQ
              Z = cmplx(1.d0,0.d0,kind=8)/D(J)
              DO I = 1,LQ
                 HLP1(I,J) = HLP1(I,J)*Z
              ENDDO
           ENDDO
           Call MMULT (HLP2,HLP1,U) 
           DO I = 1,LQ
              DO J = 1,LQ
                 GR0T(I,J) = -conjg(HLP2(J,I))
              ENDDO
           ENDDO
        ELSE
           !  UDV of HLP2*
           !  -G0T*= U2 (U D V)*^-1 V1* =  U2 U D*^-1 V*^-1 V1*
           DO I = 1,LQ
              DO J = 1,LQ
                 HLP1(I,J) = conjg(HLP2(J,I))
              ENDDO
           ENDDO
           CALL UDV_WRAP(HLP1,U,D,V,NCON)
           CALL MMULT (HLP1, U2, U)
           DO J = 1,LQ
              Z = cmplx(1.d0,0.d0,kind=8)/D(J)
              DO I = 1,LQ
                 HLP1(I,J) = HLP1(I,J)*Z
              ENDDO
           ENDDO
           CALL INV(V,HLP2,Z)
           Call MMULT(V,V1,HLP2)
           DO I = 1,LQ
              DO J = 1,LQ
                 HLP2(I,J) = conjg(V(J,I))
              ENDDO
           ENDDO
           Call MMULT (V,HLP1,HLP2)
           DO I = 1,LQ
              DO J = 1,LQ
                 GR0T(I,J) = -conjg(V(J,I))
              ENDDO
           ENDDO
        ENDIF




        ! Compute GT0
        ! GT0  =   B2 * G00   =   (  (   1   + B1* B2) * B2^-1  )^-1  =   ( B2^-1   + B1)^-1 =
        !                     =   (V2^-1 D2^-1 U2^-1  + V1 D1 U1)^-1 =  
        !                     =   ( (V2^-1 D2^-1 U2^-1 U1^-1  + V1 D1 ) U1  )^-1 = 
        !                     =   U1^-1 (  ( D2^-1 (U1 U2)^-1  + V2*V1 D1 )   )^-1 V2 
        Xmax2 = dble(cmplx(1.d0,0.d0)/D2(1))
        Xmax1 = dble(D1(1))
        Do I = 2,LQ
           X2 = dble(cmplx(1.d0,0.d0)/D2(I))
           X1 = dble(D1(I))
           If ( X2 > Xmax2 ) Xmax2 = X2
           If ( X1 > Xmax1 ) Xmax1 = X1
        ENDDO
        NVAR1 = 1
        If (Xmax2 > Xmax1 ) NVAR1 = 2
        !Write(6,*) "CGR2_1: NVAR,NVAR1 ", NVAR, NVAR1
        Call  MMULT(HLP2,U1,U2)
        DO J = 1,LQ
           DO I = 1,LQ
              HLP1(I,J) = Conjg(HLP2(J,I))
           ENDDO
        ENDDO
        DO J = 1,LQ
           DO I =1,LQ
              HLP1(I,J) =  HLP1(I,J)/D2(I) 
           ENDDO
        ENDDO
        Call MMULT(HLP2,V2,V1)
        DO J = 1,LQ
           DO I =1,LQ
              HLP2(I,J) = HLP1(I,J) +  HLP2(I,J)*D1(J) 
           ENDDO
        ENDDO
        NCON = 0
        IF ( NVAR1 == 1 ) Then
           !  UDV of HLP2
           CALL UDV_WRAP(HLP2,U,D,V,NCON)
           CALL MMULT (HLP1, V, U1) 
           CALL INV(HLP1,HLP2,Z)
           DO J = 1,LQ
              Z = cmplx(1.d0,0.d0)/D(J)
              DO I = 1,LQ
                 HLP2(I,J) = HLP2(I,J)*Z
              ENDDO
           ENDDO
           DO I = 1,LQ
              DO J = 1,LQ
                 HLP1(I,J) = Conjg(U(J,I))
              ENDDO
           ENDDO
           CALL MMULT(U,HLP1,V2)
           Call MMULT (GRT0, HLP2,U)
        ELSE
           !UDV of HLP2^*
           DO J = 1,LQ
              DO I =1,LQ
                 HLP1(I,J) = Conjg(HLP2(J,I))
              ENDDO
           ENDDO
           CALL UDV_WRAP(HLP1,U,D,V,NCON) 
           DO I = 1,LQ
              DO J = 1,LQ
                 HLP1(I,J) = conjg(U1(J,I))
              ENDDO
           ENDDO
           CALL MMULT( HLP2, HLP1,U)
           DO J = 1,LQ
              DO I = 1,LQ
                 HLP2(I,J) = HLP2(I,J)/Conjg(D(J))
              ENDDO
           ENDDO
           DO I = 1,LQ
              DO J = 1,LQ
                 HLP1(I,J) = conjg(V(J,I))
              ENDDO
           ENDDO
           CALL INV(HLP1,V,Z)
           CALL MMULT(U,V,V2)
           Call MMULT (GRT0, HLP2,U)
        ENDIF
        Xmin = abs(dble(D(1)))
        DO I = 1,LQ
           if (abs(dble(D(I))) <   Xmin ) Xmin = abs(dble(D(I)))
        ENDDO
        Write(6,*) 'Cgr2_1 T0, Xmin: ', Xmin


        !Compute GRTT
        Z  = cmplx(1.d0,0.d0,kind=8)
        Z1 = cmplx(1.d0,0.d0,kind=8)
        CALL CGR(Z,NVAR,GRTT, U2,D2,V2, U1,D1,V1)

 
      END SUBROUTINE CGR2_1


!!$        ! Compute G0T
!!$        ! G00 = (1 + B1*B2)^-1   G0T = -(1 -  (1 + B1*B2)^-1 )*B2^-1              =
!!$        !                            = -( 1 +  B1*B2 - 1)  (1 + B1*B2)^-1 * B2^-1 =
!!$        !                            = - B1 * B2 (1 + B1*B2)^-1 * B2^-1           = 
!!$        !                            = -( B2 ( 1+ B1*B2) * B2^-1 B1^-1)^-1        =
!!$        !                            = -( B2 ( 1+ B1*B2) * (B1 B2)^-1)^-1         = 
!!$        !                            = -( B2 ( (B1 B2)^-1  + 1 ) )^-1             =
!!$        !                            = -( B1^-1   +   B2)^-1                      =
!!$        !                            = -( U1^-1 D1^-1 V1^-1  + U2 D2 V2)^-1       =
!!$        !                            = -( ( U1^-1 D1^-1 V1^-1 V2^-1 + U2 D2 ) V2 )^-1 =
!!$        !                            = -( U1^-1( D1^-1 (V2 V1)^-1 + U1 U2 D2) V2 )^-1 = 
!!$        !                            = -  V2^-1( D1^-1 (V2 V1)^-1 + U1 U2 D2)^-1  U1 
!!$        !       B2 = U2*D2*V2
!!$        !       B1 = V1*D1*U1
!!$        Call MMULT (HLP2, V1inv,V2inv)
!!$        DO J = 1,LQ
!!$           DO I = 1,LQ
!!$              HLP2(I,J) =   HLP2(I,J)/D1(I) 
!!$           ENDDO
!!$        ENDDO
!!$        Call MMULT (HLP1, U1,U2)
!!$        DO J = 1,LQ
!!$           DO I = 1,LQ
!!$              HLP2(I,J) =   HLP2(I,J)  + HLP1(I,J)*D2(J)
!!$           ENDDO
!!$        ENDDO
!!$        Xmax2 = dble(cmplx(1.d0,0.d0,Kind=8)/D1(1))
!!$        Xmax1 = dble(D2(1))
!!$        Do I = 2,LQ
!!$           X2 = dble(cmplx(1.d0,0.d0,Kind=8)/D1(I))
!!$           X1 = dble(D2(I))
!!$           If ( X2 > Xmax2 ) Xmax2 = X2
!!$           If ( X1 > Xmax1 ) Xmax1 = X1
!!$        ENDDO
!!$        NVAR1 = 1
!!$        If (Xmax2 > Xmax1 ) NVAR1 = 2
!!$        IF (NVAR1 == 1) Then 
!!$           ! UDV of HLP2
!!$           != -  V2^-1( U D V )^-1 ) U1  = 
!!$           != - V2^-1 V^-1 D^-1 U^-1 U1 = - (V V2)^-1 D^-1  U^-1 U1
!!$           CALL UDV(HLP2,U,D,V,NCON) 
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 HLP1(I,J) = conjg(U(J,I))
!!$              ENDDO
!!$           ENDDO
!!$           CALL MMULT( U, HLP1, U1 )
!!$           CALL MMULT(HLP1,V, V2)
!!$           CALL INV  (HLP1,V  ,Z)
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 HLP1(I,J) = - V(I,J)/D(J)
!!$              ENDDO
!!$           ENDDO
!!$           CALL MMULT(GR0T, HLP1, U)           
!!$        Else
!!$           ! UDV of HLP2^*
!!$           != -  V2^-1( U D V)^*,-1 ) U1  = 
!!$           != -  V2^-1 U  D^-1 V^*,-1 U1 
!!$           DO I = 1,LQ
!!$              DO J = 1,LQ
!!$                 HLP1(J,I) = Conjg(HLP2(I,J))
!!$              ENDDO
!!$           ENDDO
!!$           CALL UDV(HLP1,U,D,V,NCON)
!!$           CALL INV(V2,HLP1,Z)
!!$           CALL MMULT(HLP2,HLP1,U)
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 HLP1(I,J) = -HLP2(I,J)/conjg(D(J))
!!$              ENDDO
!!$           ENDDO
!!$           CALL INV(V,HLP2,Z)
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 V(I,J) = Conjg(HLP2(J,I))
!!$              ENDDO
!!$           ENDDO
!!$           CALL MMULT(HLP2,V,U1)
!!$           CALL MMULT(GR0T, HLP1,HLP2)           
!!$        endif
!!$        Xmin = abs(dble(D(1)))
!!$        DO I = 1,LQ
!!$           if (abs(dble(D(I))) <   Xmin ) Xmin = abs(dble(D(I)))
!!$        ENDDO
!!$        Write(6,*) 'Cgr2_1 0T, Xmin: ', Xmin




!!$        ! Compute G0T
!!$        ! G00 = (1 + B1*B2)^-1   G0T = -(1 -  (1 + B1*B2)^-1 )*B2^-1              =
!!$        !                            = -( 1 +  B1*B2 - 1)  (1 + B1*B2)^-1 * B2^-1 =
!!$        !                            = - B1 * B2 (1 + B1*B2)^-1 * B2^-1           = 
!!$        !                            = -( B2 ( 1+ B1*B2) * B2^-1 B1^-1)^-1        =
!!$        !                            = -( B2 ( 1+ B1*B2) * (B1 B2)^-1)^-1         = 
!!$        !                            = -( B2 ( (B1 B2)^-1  + 1 ) )^-1             =
!!$        !                            = -( B1^-1   +   B2)^-1                      =
!!$        !                            = -( U1^-1 D1^-1 V1^-1  + U2 D2 V2)^-1       =
!!$        !                            = -(U2 (U2^-1 U1^-1 D1^-1  +  D2 V2 V1 ) V1^-1  )^-1   =
!!$        !                            = - V1 ( (U1 U2)^-1 D1^-1  +  D2 V2 V1 )^-1 U2^-1      
!!$        !       B2 = U2*D2*V2
!!$        !       B1 = V1*D1*U1
!!$        Call MMULT (HLP1, U1,U2)
!!$        DO J = 1,LQ
!!$           DO I = 1,LQ 
!!$              HLP2(I,J) =   Conjg(HLP1(J,I))
!!$           ENDDO
!!$        ENDDO
!!$        DO J = 1,LQ
!!$           DO I = 1,LQ 
!!$              HLP2(I,J) =    HLP2(I,J) / D1(J)
!!$           ENDDO
!!$        ENDDO
!!$
!!$        Call MMULT (HLP1, V2,V1)
!!$        DO J = 1,LQ
!!$           DO I = 1,LQ
!!$              HLP2(I,J) =   HLP2(I,J)  + D2(I)*HLP1(I,J)
!!$           ENDDO
!!$        ENDDO
!!$        Xmax2 = dble(cmplx(1.d0,0.d0)/D1(1))
!!$        Xmax1 = dble(D2(1))
!!$        Do I = 2,LQ
!!$           X2 = dble(cmplx(1.d0,0.d0)/D1(I))
!!$           X1 = dble(D2(I))
!!$           If ( X2 > Xmax2 ) Xmax2 = X2
!!$           If ( X1 > Xmax1 ) Xmax1 = X1
!!$        ENDDO
!!$        NVAR1 = 1
!!$        If (Xmax1 > Xmax2 ) NVAR1 = 2
!!$        IF (NVAR1 == 1) Then 
!!$           ! UDV of HLP2
!!$           != - V1 ( U D V)^-1 U2^-1      
!!$           != - V1 V^-1 D^-1 U^-1 U2^-1 = - V1 V^-1 D^-1 (U2 U)^-1
!!$           CALL UDV(HLP2,U,D,V,NCON) 
!!$           CALL MMULT( HLP2, U2, U )
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 HLP1(I,J) = conjg(HLP2(J,I))
!!$              ENDDO
!!$           ENDDO
!!$           CALL INV  (V, HLP2  ,Z)
!!$           CALL MMULT(V, V1, HLP2)
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 HLP2(I,J) = - V(I,J)/D(J)
!!$              ENDDO
!!$           ENDDO
!!$           CALL MMULT(GR0T, HLP2, HLP1)           
!!$        Else
!!$           ! UDV of HLP2^*
!!$           != - V1 ( U D V)^(*,-1) U2^-1      
!!$           != - V1 U D^(*,-1) V^(*,-1)  U2^-1 
!!$           DO I = 1,LQ
!!$              DO J = 1,LQ
!!$                 HLP1(J,I) = Conjg(HLP2(I,J))
!!$              ENDDO
!!$           ENDDO
!!$           CALL UDV(HLP1,U,D,V,NCON)
!!$           CALL MMULT(HLP2,V1,U)
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 HLP2(I,J) = -HLP2(I,J)/conjg(D(J))
!!$              ENDDO
!!$           ENDDO
!!$
!!$           CALL INV(V,HLP1,Z)
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 V(I,J) = Conjg(HLP1(J,I))
!!$              ENDDO
!!$           ENDDO
!!$           DO J = 1,LQ
!!$              DO I = 1,LQ
!!$                 U(I,J) = Conjg(U2(J,I))
!!$              ENDDO
!!$           ENDDO
!!$           CALL MMULT(HLP1,V,U)
!!$           
!!$           CALL MMULT(GR0T, HLP2,HLP1)           
!!$        endif
!!$        Xmin = abs(dble(D(1)))
!!$        DO I = 1,LQ
!!$           if (abs(dble(D(I))) <   Xmin ) Xmin = abs(dble(D(I)))
!!$        ENDDO
!!$        Write(6,*) 'Cgr2_1 0T, Xmin: ', Xmin, NVAR1
