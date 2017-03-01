   Program UDV_Test

     Use MyMats
     Use UDV_Wrap_mod
     Use Random_Wrap

     COMPLEX (KIND=KIND(0.D0)), allocatable :: U(:,:), D(:), V(:,:), A(:,:)
     Integer :: i, j, LQ, LQ1,I1, J1
     Integer :: Iseed(1)
     
     
     ISEED(1) = 74392
     CALL RANSET(ISEED)
     LQ = 100; LQ1=LQ
     ALLOCATE(A(LQ,LQ1),U(LQ,LQ1), D(LQ1), V(LQ1,LQ1))
     NCON = 1
     DO I = 1,LQ
        DO J = 1,LQ1
           A(I,J) = CMPLX(RANF_WRAP(),RANF_WRAP(),KIND(0.d0))
        ENDDO
     ENDDO
     !WRITE(6,*) A
     CALL UDV_WRAP_Pivot(A,U,D,V,NCON,LQ,LQ1)

   end Program UDV_TEST
