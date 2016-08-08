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
        Complex  (Kind=double) :: Z
        
        Integer :: LQ2, I,J, M, ILQ, JLQ, NCON, I1, J1
        
        LQ2 = LQ*2
        
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
        !	U3B^T * ( V1INV  0   )   = U3B
        !	        ( 0      U2^T )
        
        DO I = 1,LQ2
           DO J = 1,LQ2
              HLPB1(I,J) = conjg(U3B(J,I))
           ENDDO
        ENDDO
        HLPB2 = cmplx(0.D0,0.d0,double)
        DO I = 1,LQ
           DO J = 1,LQ
              HLPB2(I,J) = V1INV(I,J)
           ENDDO
        ENDDO
        DO I = 1,LQ
           ILQ = I + LQ
           DO J = 1,LQ
              JLQ = J + LQ
              HLPB2(ILQ,JLQ) = conjg(U2(J,I))
           ENDDO
        ENDDO
        CALL MMULT(U3B,HLPB1,HLPB2)
        
        
        !       Multiplication:
        !	( V2INV   0    )*(V3B)^{-1}   = V3B
        !	( 0       U1^T )
        
        CALL INV(V3B,HLPB1,Z)
        HLPB2 = cmplx(0.d0,0.d0,double)
        DO I = 1,LQ
           DO J = 1,LQ
              HLPB2(I,J) = V2INV(I,J)
           ENDDO
        ENDDO
        DO I = 1,LQ
           ILQ = I + LQ
           DO J = 1, LQ
              JLQ = J + LQ
              HLPB2(ILQ,JLQ) = conjg(U1(J,I))
           ENDDO
        ENDDO
        CALL MMULT(V3B,HLPB2,HLPB1)
        
        
        ! G = V3B * D3B^{-1}* U3B
        DO M = 1,LQ2
           Z = cone/D3B(M)
           DO J = 1,LQ2
              U3B(M,J) =   Z * U3B(M,J) 
           ENDDO
        ENDDO
        CALL MMULT(HLPB2, V3B, U3B)
        DO I = 1,LQ
           I1 = I+LQ
           DO J = 1,LQ
              J1 = J + LQ
              GR00(I,J) = HLPB2(I ,J )
              GRTT(I,J) = HLPB2(I1,J1)
              GRT0(I,J) = HLPB2(I1,J )
              GR0T(I,J) = HLPB2(I,J1 )
           ENDDO
	ENDDO
 
    END SUBROUTINE CGR2
