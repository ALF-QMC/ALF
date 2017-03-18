!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief
!> This function scales every Row of M with an entry from the vector D
!>
!> @note
!> GCC does not by itself pull out the division by D from the inner loop
!> in every other expression that I tried, therefore this explicit loop construct
!> is necessary.
!>
!> @param[inout] M The matrix whose rows we scale
!> @param[in] D A vector with the scaling factors
!> @param[in] A boolean variable wether we do a ccomplex conjugate on D
!> @param[in] LQ The length of the vectors and the dimension of M
!--------------------------------------------------------------------
        Subroutine scalematrix(M, D, c, LQ)
        Integer, Intent(In) :: LQ
        Integer :: J
        Complex (Kind = Kind(0.D0)), Intent(Inout) :: M(LQ,LQ)
        Complex (Kind = Kind(0.D0)), Intent(In) :: D(LQ)
        Complex (Kind = Kind(0.D0)) :: Z
        Logical, Intent(In) :: c

        DO J = 1,LQ
            Z = 1.D0/D(J)
            If (c) Z = Conjg(z)
            M(:,J) = M(:,J)*Z
        ENDDO

        end Subroutine

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
        Use MyMats
        USe UDV_Wrap_mod

        Implicit none
        
        Interface
           SUBROUTINE CGR(Z,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(In) ::  URUP, VRUP, ULUP, VLUP
             COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(In)   ::  DLUP, DRUP
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(INOUT) :: GRUP
        
             COMPLEX(Kind=Kind(0.d0)) :: Z
           END SUBROUTINE CGR
        end Interface


        !  Arguments
        Integer,  intent(in) :: LQ, NVAR
        Complex (Kind=Kind(0.d0)), intent(in)    :: U1(LQ,LQ), V1(LQ,LQ), U2(LQ,LQ), V2(LQ,LQ)
        Complex (Kind=Kind(0.d0)), intent(in)    :: D2(LQ), D1(LQ)
        Complex (Kind=Kind(0.d0)), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)


        ! Local::
        Complex  (Kind=Kind(0.d0)) :: HLP1(LQ,LQ), HLP2(LQ,LQ), U(LQ,LQ), D(LQ), V(LQ,LQ)
        Complex  (Kind=Kind(0.d0)) :: Z, Z1, Z2
        Real     (Kind=Kind(0.d0)) :: Xmax, Xmin, Xmax1, Xmax2, Xmean
        Integer                :: I, J, NCON, NVAR1

        Complex  (Kind=Kind(0.d0)) :: V2inv(LQ,LQ), V1inv(LQ,LQ), alpha, beta
        

        NCON = 0

        Call INV( V2, V2inv, Z2)
        CALL INV( V1, V1inv, Z1)
        alpha = 1.D0
        beta = 0.D0

        CALL ZGEMM('C', 'N', LQ, LQ, LQ, alpha, U1, LQ, U1, LQ, beta, HLP2, LQ)
        CALL ZLASET('A', LQ, LQ, alpha, beta, HLP1, LQ)
        Xmax = 0.d0
        CALL COMPARE(HLP1, HLP2, XMAX, XMEAN)

        CALL ZGEMM('C', 'N', LQ, LQ, LQ, alpha, U2, LQ, U2, LQ, beta, HLP2, LQ)
        CALL ZLASET('A', LQ, LQ, alpha, beta, HLP1, LQ)
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
            HLP2(:,J) = D1(:) * HLP2(:, J) * D2(J) + HLP1(:,J)
        ENDDO
        
        Xmax1 = maxval(dble(D1))
        Xmax2 = maxval(dble(D2))
        Nvar1 = 1
        If ( Xmax2 > Xmax1) Nvar1 = 2
        If (Nvar1 == 1) then
           !  V2^-1 (UDV  )^-1 V1^-1 =  V2^-1 V^-1 D^-1 U^-1 V1^-1
           Call UDV_WRAP(HLP2, U, D, V, Ncon) 
           CALL INV  (V,HLP2 ,Z   )
           CALL MMULT(HLP1,V2inv,HLP2)
           CALL ZGEMM('C', 'N', LQ, LQ, LQ, alpha, U, LQ, V1Inv, LQ, beta, HLP2, LQ)
        else
           !  V2^-1 (UDV  )^(-1,*) V1^-1 =  V2^-1 U  D^-1 V^(-1,*) V1^-1
           HLP1 = CT(HLP2)
           Call UDV_WRAP(HLP1, U, D, V, Ncon) 
           CALL INV (V, HLP1, Z)
           CALL ZGEMM('C', 'N', LQ, LQ, LQ, alpha, HLP1, LQ, V1Inv, LQ, beta, HLP2, LQ)
           Call MMULT(HLP1, V2inv, U)
        endif
        Call SCALEMATRIX(HLP1, D, .FALSE., LQ)
        CALL MMULT(GR00,HLP1,HLP2)
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
        Xmax2 = maxval(dble(1.D0/D1))
        Xmax1 = maxval(dble(D2))
        NVAR1 = 1
        If (Xmax2  > Xmax1)  Nvar1 = 2
        Call  MMULT(HLP1,U1,U2)
!TODO: Consider benchmarking wether it is beneficial to interchange loops. Here it is saving divisions vs. proper mem access
        DO J = 1,LQ
           HLP1(:, J) =  HLP1(:,J)/conjg(D1(:))
        ENDDO
        CALL ZGEMM('C', 'C', LQ, LQ, LQ, alpha, V1, LQ, V2, LQ, beta, HLP2, LQ)
        DO J = 1,LQ
            HLP2(:, J) = HLP1(:,J) + HLP2(:,J) * conjg(D2(J))
        ENDDO
        NCON = 0
        IF ( NVAR1 == 1 ) Then
           !  UDV of HLP2
           !  -G0T*= U2 V^-1 D^-1 U* V1*
           CALL UDV_WRAP(HLP2,U,D,V,NCON)
           CALL INV(V,HLP2,Z)
           CALL MMULT (V, V1, U) ! V = V1 * U, reuse of the variable V as temporary storage
           Call MMULT(HLP1,U2,HLP2)
        ELSE
           !  UDV of HLP2*
           !  -G0T*= U2 (U D V)*^-1 V1* =  U2 U D*^-1 V*^-1 V1*
           HLP1 = CT(HLP2)
           CALL UDV_WRAP(HLP1,U,D,V,NCON)
           CALL INV(V,HLP2,Z)
           Call MMULT(V,V1,HLP2)
           CALL MMULT (HLP1, U2, U)
        ENDIF
        CALL SCALEMATRIX(HLP1, (-D), .FALSE., LQ)
        CALL ZGEMM('N', 'C', LQ, LQ, LQ, alpha, V, LQ, HLP1, LQ, beta, GR0T, LQ)


        ! Compute GT0
        ! GT0  =   B2 * G00   =   (  (   1   + B1* B2) * B2^-1  )^-1  =   ( B2^-1   + B1)^-1 =
        !                     =   (V2^-1 D2^-1 U2^-1  + V1 D1 U1)^-1 =  
        !                     =   ( (V2^-1 D2^-1 U2^-1 U1^-1  + V1 D1 ) U1  )^-1 = 
        !                     =   U1^-1 (  ( D2^-1 (U1 U2)^-1  + V2*V1 D1 )   )^-1 V2 
        Xmax2 = maxval(dble(1.D0/D2))
        Xmax1 = maxval(dble(D1))
        NVAR1 = 1
        If (Xmax2 > Xmax1 ) NVAR1 = 2
        !Write(6,*) "CGR2_1: NVAR,NVAR1 ", NVAR, NVAR1
        CALL ZGEMM('C', 'C', LQ, LQ, LQ, alpha, U2, LQ, U1, LQ, beta, HLP1, LQ)
!TODO: Consider benchmarking wether it is beneficial to interchange loops. Here it is saving divisions vs. proper mem access
        DO J = 1,LQ
           HLP1(:, J) =  HLP1(:,J)/D2(:)
        ENDDO
        Call MMULT(HLP2,V2,V1)
        DO J = 1,LQ
            HLP2(:, J) = HLP1(:,J) + HLP2(:,J) * D1(J)
        ENDDO
        NCON = 0
        IF ( NVAR1 == 1 ) Then
           !  UDV of HLP2
           CALL UDV_WRAP(HLP2,U,D,V,NCON)
           CALL MMULT (HLP1, V, U1) 
           CALL INV(HLP1,HLP2,Z)
           CALL ZGEMM('C', 'N', LQ, LQ, LQ, alpha, U, LQ, V2, LQ, beta, V, LQ)
        ELSE
           !UDV of HLP2^*
           HLP1 = CT(HLP2)
           CALL UDV_WRAP(HLP1,U,D,V,NCON)
           CALL ZGEMM('C', 'N', LQ, LQ, LQ, alpha, U1, LQ, U, LQ, beta, HLP2, LQ)
           CALL INV(V,HLP1,Z)
           CALL ZGEMM('C', 'N', LQ, LQ, LQ, alpha, HLP1, LQ, V2, LQ, beta, V, LQ)
        ENDIF
        CALL SCALEMATRIX(HLP2, D, NVAR .ne. 1, LQ)
        Xmin = minval(abs(dble(D)))
        Call MMULT (GRT0, HLP2, V)
        Write(6,*) 'Cgr2_1 T0, Xmin: ', Xmin

        !Compute GRTT
        Z  = cmplx(1.d0,0.d0,Kind=Kind(0.d0))
        Z1 = cmplx(1.d0,0.d0,Kind=Kind(0.d0))
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
!!$        Xmax2 = dble(cmplx(1.d0,0.d0,Kind=Kind(0.d0))/D1(1))
!!$        Xmax1 = dble(D2(1))
!!$        Do I = 2,LQ
!!$           X2 = dble(cmplx(1.d0,0.d0,Kind=Kind(0.d0))/D1(I))
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
