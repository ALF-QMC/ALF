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


        ! Local::
        Complex (Kind=double) :: U3B(2*LQ,2*LQ), V3B(2*LQ,2*LQ), HLPB1(2*LQ,2*LQ), HLPB2(2*LQ,2*LQ), &
             &                   V2INV(LQ,LQ), V1INV(LQ,LQ), HLP2(LQ,LQ)
        Complex  (Kind=double) :: D3B(2*LQ)
        Complex  (Kind=double) :: Z, alpha, beta
        Integer, dimension(:), allocatable :: IPVT
        
        Integer :: LQ2, I,J, M, ILQ, JLQ, NCON
        
        LQ2 = LQ*2
        alpha = 1.D0
        beta = 0.D0
        HLPB1 = cmplx(0.D0,0.d0,double)
        DO I = 1,LQ
           HLPB1(I   , I + LQ ) =  D1(I)
           HLPB1(I+LQ, I      ) = -D2(I)
        ENDDO
        CALL INV(V2,V2INV,Z)
        CALL INV(V1,V1INV,Z)
        CALL MMULT(HLP2,V1INV,V2INV)
        DO J = 1,LQ
           DO I = 1,LQ
              HLPB1(I,J) = HLP2(I,J)
           ENDDO
        ENDDO
        CALL MMULT(HLP2,U1,U2)
        DO I = 1,LQ
           ILQ = I+LQ
           DO J = 1,LQ
              JLQ = J + LQ 
              HLPB1(ILQ,JLQ) = conjg( HLP2(J,I) ) ! = (U1*U2)^T
           ENDDO
        ENDDO
        NCON = 0
        CALL UDV_wrap(HLPB1,U3B,D3B,V3B,NCON)
        
        !       Multiplication:
        !	( V2INV   0    )*(V3B)^{-1}   = V3B
        !	( 0       U1^T )
        
        HLPB2 = cmplx(0.d0,0.d0,double)
        call ZLACPY('A', LQ, LQ, V2INV, LQ, HLPB2, LQ2)
        DO I = 1,LQ
           ILQ = I + LQ
           DO J = 1, LQ
              JLQ = J + LQ
              HLPB2(ILQ,JLQ) = conjg(U1(J,I))
           ENDDO
        ENDDO
        allocate(IPVT(LQ2))
        call ZGESV(LQ2, LQ2, V3B, LQ2, IPVT, HLPB2, LQ2, I)
        call ZLACPY('A', LQ2, LQ2, HLPB2, LQ2, V3B, LQ2)
        deallocate(IPVT)
        !       Multiplication:
        !	U3B^T * ( V1INV  0   )   = U3B
        !	        ( 0      U2^T )
        
        HLPB2 = cmplx(0.D0,0.d0,double)
        call ZLACPY('A', LQ, LQ, V1INV, LQ, HLPB2, LQ2)
        DO I = 1,LQ
           ILQ = I + LQ
           DO J = 1,LQ
              JLQ = J + LQ
              HLPB2(ILQ,JLQ) = conjg(U2(J,I))
           ENDDO
        ENDDO
        CALL ZGEMM('C', 'N', LQ2, LQ2, LQ2, alpha, U3B, LQ2, HLPB2, LQ2, beta, HLPB1, LQ2)
        ! G = V3B * D3B^{-1}* U3B
        DO M = 1,LQ2
           Z = cone/D3B(M)
           DO J = 1,LQ2
              HLPB1(M,J) =   Z * HLPB1(M,J) 
           ENDDO
        ENDDO
        call MMULT(HLPB2, V3B, HLPB1)
        call get_blocks(GR00, GR0T, GRT0, GRTT, HLPB1, LQ)
 
    END SUBROUTINE CGR2
