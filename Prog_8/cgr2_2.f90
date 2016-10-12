      SUBROUTINE CGR2_2(GRT0, GR00, GRTT, GR0T, U2, D2, V2, U1, D1, V1, LQ)

      
        !       B2 = U2*D2*V2 is right (i.e. from time slice 0 to tau) propagation to time tau 
        !       B1 = V1*D1*U1 is left  (i.e. from time slice Ltrot to tau) propagation to time tau
        !Calc:      (  1   B1 )^-1       ( G00    G0T )   
        !           (-B2   1  )     ==   ( GT0    GTT )
        !            
        !       G00 = (1 + B1*B2)^-1   G0T = -(1 -  G00)*B2^-1
        !       GT0  =   B2 * G00      GTT = (1 + B2*B1)^-1
 
        !(  1   V1*D1*U1 )^-1      ( ( V1   0 )   ( V1^-1     D1*U1 ) )^-1
        !(-U2*D2*V2   1  )      == ( ( 0   U2 ) * (-D2*V2     U2^-1 ) )       == I
        !  You should transpose before carrying out the singular value decomposition 
        !
        !     
        !       ( ( V1   0 )   ( V1^-1     D1*U1 )^*^* )^-1                      ( V1^-1    0 ) 
        ! I ==  ( ( 0   U2 ) * (-D2*V2     U2^-1 )     )      =   (UDV^*)^(-1) * ( 0     U2^-1) =
        !
        !                                      ( V1^-1    0 )
        !  ==   U * D^(*,-1) * V^(*,-1)     *  ( 0     U2^-1)

        ! Let's see if this could work.
        Use Precdef
        Use MyMats
        Use UDV_WRAP_mod
        Implicit none

        !  Arguments
        Integer,  intent(in) :: LQ
        Complex (Kind=double), intent(in)    :: U1(LQ,LQ), V1(LQ,LQ), U2(LQ,LQ), V2(LQ,LQ)
        Complex (Kind=double), intent(in)    :: D2(LQ), D1(LQ)
        Complex (Kind=double), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)


        ! Local::
        Complex  (Kind=double) :: U3B(2*LQ,2*LQ), V3B(2*LQ,2*LQ), HLPB1(2*LQ,2*LQ), HLPB2(2*LQ,2*LQ), &
             &                   V1INV(LQ,LQ)
        Complex  (Kind=double) :: D3B(2*LQ)
        Complex  (Kind=double) :: Z, alpha, beta
        Complex(Kind = Kind(0.D0)), allocatable, Dimension(:) :: TMPVEC
        Complex(Kind = Kind(0.D0)), allocatable, Dimension(:, :) :: MYU2
        INTEGER, Dimension(:), Allocatable :: IPVT

        Integer :: LQ2, I,J, NCON, info
        
        LQ2 = LQ*2
        NCON = 0
        alpha = 1.D0
        beta = 0.D0
        ALLOCATE(TMPVEC(LQ2), MYU2(LQ, LQ), IPVT(LQ2))
        MYU2 = CONJG(TRANSPOSE(U2))
        CALL INV(V1,V1INV,Z)
        If (dble(D1(1)) >  dble(D2(1)) ) Then 
           !Write(6,*) "D1(1) >  D2(1)", dble(D1(1)), dble(D2(1))
           DO J = 1,LQ
              DO I = 1,LQ
                 HLPB2(I   , J    ) =  V1INV(I,J)
                 HLPB2(I   , J+LQ ) =  D1(I)*U1(I,J)
                 HLPB2(I+LQ, J+LQ ) =  MYU2(I, J)
                 HLPB2(I+LQ, J    ) = -D2(I)*V2(I,J)
              ENDDO
           ENDDO
           HLPB1 = CT(HLPB2)

           !CALL UDV_wrap(HLPB1,U3B,D3B,V3B,NCON)
           CALL UDV_wrap_Pivot(HLPB1,U3B,D3B,V3B,NCON,LQ2,LQ2)
           TMPVEC = conjg(1.D0/D3B)
!!$!!!!!!!!!!!!!  Tests
!!$        Xmax = 0.d0
!!$        DO I = 1,LQ2
!!$           DO J = 1,LQ2
!!$              Z = cmplx(0.d0,0.d0)
!!$              DO N = 1,LQ2
!!$                 Z = Z + U3B(I,N) *conjg(U3B(J,N))
!!$              ENDDO
!!$              if (I == J)  Z = Z - cmplx(1.d0,0.d0)
!!$              X = real(SQRT( Z* conjg(Z)),kind=8)
!!$              if (X > Xmax) Xmax = X
!!$           ENDDO
!!$        ENDDO
!!$        !Write(6,*) 'Cgr2_2, ortho: ', Xmax
!!$        DO I = 1,LQ2
!!$           Z =  D3B(I)
!!$           if (I == 1)  Xmax = real(SQRT( Z* conjg(Z)),kind=8) 
!!$           if ( real(SQRT( Z* conjg(Z)),kind=8)  < Xmax ) Xmax = &
!!$                & real(SQRT( Z* conjg(Z)),kind=8)
!!$        ENDDO
!!$        !Write(6,*) 'Cgr2_2, Cutoff: ', Xmax
!!$!!!!!!!!!!!!! End Tests
           HLPB1 = cmplx(0.d0,0.d0,double)
           DO I = 1,LQ
              DO J = 1,LQ
                 HLPB1(I   , J    ) =  V1INV(I,J)
                 HLPB1(I+LQ, J+LQ ) =  MYU2(I, J)
              ENDDO
           ENDDO
           CALL ZGETRF(LQ2, LQ2, V3B, LQ2, IPVT, info)
           CALL ZGETRS('C', LQ2, LQ2, V3B, LQ2, IPVT, HLPB1, LQ2, info)! Block structure of HLPB1 is not exploited
           DO J = 1,LQ2
              DO I = 1,LQ2
                 HLPB1(I,J)  = TMPVEC(I)*HLPB1(I,J)
              ENDDO
           ENDDO
           CALL get_blocks_of_prod(GR00, GR0T, GRT0, GRTT, U3B, HLPB1, LQ)
        Else
           !Write(6,*) "D1(1) <  D2(1)", dble(D1(1)), dble(D2(1))
           DO J = 1,LQ
              DO I = 1,LQ
                 HLPB2(I   , J    ) =  MYU2(I, J)
                 HLPB2(I   , J+LQ ) = -D2(I)*V2(I,J)
                 HLPB2(I+LQ, J+LQ ) =  V1INV(I,J)
                 HLPB2(I+LQ, J    ) =  D1(I)*U1(I,J)
              ENDDO
           ENDDO
           HLPB1 = CT(HLPB2)
           
           !CALL UDV_wrap(HLPB1,U3B,D3B,V3B,NCON)
           CALL UDV_wrap_Pivot(HLPB1,U3B,D3B,V3B,NCON,LQ2,LQ2)
           TMPVEC = conjg(1.D0/D3B)
           HLPB1 = cmplx(0.d0,0.d0,double)
           DO I = 1,LQ
              DO J = 1,LQ
                 HLPB1(I   , J    ) =  MYU2(I, J)
                 HLPB1(I+LQ, J+LQ ) =  V1INV(I,J)
              ENDDO
           ENDDO
           CALL ZGETRF(LQ2, LQ2, V3B, LQ2, IPVT, info)
           CALL ZGETRS('C', LQ2, LQ2, V3B, LQ2, IPVT, HLPB1, LQ2, info)! Block structure of HLPB1 is not exploited
           DO J = 1,LQ2
              DO I = 1,LQ2
                 HLPB1(I,J)  = TMPVEC(I)*HLPB1(I,J)
              ENDDO
           ENDDO
           call get_blocks_of_prod(GRTT, GRT0, GR0T, GR00, U3B, HLPB1, LQ)
        Endif
        DEALLOCATE(TMPVEC, MYU2, IPVT)
      END SUBROUTINE CGR2_2
