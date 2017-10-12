! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a
!
PROGRAM SLMATMULTS

        IMPLICIT NONE
        INTEGER                   :: N, M, M2, I, J, K
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: A
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: Mat, MatL, MatR
        INTEGER                  , DIMENSION(:)  , ALLOCATABLE :: P
        COMPLEX (KIND=KIND(0.D0)) :: Z
        
        M=400
        M2=600
        
        ALLOCATE(Mat(M,M2),MatL(M,M2),MatR(M,M2))
        
        Do N=1,6
          
          ALLOCATE(A(N,N),P(N))
          
          P(1)=1
          if (N>1) P(2)=2
          if (N>2) P(3)=6
          if (N>3) P(4)=5
          if (N>4) P(5)=M-1
          if (N>5) P(6)=M
          
          !Set up small matrix a as a hermitian one
          DO I=1,N
          DO J=1,N
            A(I,J) = CMPLX(DBLE(i+j),DBLE(i-j),kind(0.d0))
          ENDDO
          ENDDO
          
          !test left-mult using upper
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatL(:,:)=Mat(:,:)
          DO I=1,N
          DO J=1,M2
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + A(I,K)*Mat(P(K),J) 
          ENDDO
          MatL(P(I),J) = Z
          ENDDO
          ENDDO
          
          CALL ZSLHEMM('L','u', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatL))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLHEMM left upper"
            write(*,*) sum(abs(Mat-MatL))/dble(M*M2), " is larger then 1E-14!"
            stop 2
          endif
          
          !test left-mult using lower
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatL(:,:)=Mat(:,:)
          DO I=1,N
          DO J=1,M2
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + A(I,K)*Mat(P(K),J) 
          ENDDO
          MatL(P(I),J) = Z
          ENDDO
          ENDDO
          
          CALL ZSLHEMM('L','l', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatL))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLHEMM left lower"
            write(*,*) sum(abs(Mat-MatL))/dble(M*M2), " is larger then 1E-14!"
            stop 3
          endif
          
          !test right-mult using upper
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatR(:,:)=Mat(:,:)
          DO I=1,M
          DO J=1,N
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + Mat(I,P(K)) * A(k,j)
          ENDDO
          MatR(I,P(J)) = Z
          ENDDO
          ENDDO
          
          CALL ZSLHEMM('R','u', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatR))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLHEMM right upper"
            write(*,*) sum(abs(Mat-MatR))/dble(M*M2), " is larger then 1E-14!"
            stop 4
          endif
          
          !test right-mult using lower
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatR(:,:)=Mat(:,:)
          DO I=1,M
          DO J=1,N
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + Mat(I,P(K)) * A(k,j)
          ENDDO
          MatR(I,P(J)) = Z
          ENDDO
          ENDDO
          
          CALL ZSLHEMM('R','l', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatR))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLHEMM right lower"
            write(*,*) sum(abs(Mat-MatR))/dble(M*M2), " is larger then 1E-14!"
            stop 5
          endif
          
          !test left-mult op='n'
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatL(:,:)=Mat(:,:)
          DO I=1,N
          DO J=1,M2
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + A(i,k)*Mat(P(K),J) 
          ENDDO
          MatL(P(I),J) = Z
          ENDDO
          ENDDO
          
          CALL ZSLGEMM('L','n', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatL))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLGEMM left no-transpose"
            write(*,*) sum(abs(Mat-MatL))/dble(M*M2), " is larger then 1E-14!"
            stop 6
          endif
          
          !test left-mult op='t'
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatL(:,:)=Mat(:,:)
          DO I=1,N
          DO J=1,M2
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + A(k,i)*Mat(P(K),J) 
          ENDDO
          MatL(P(I),J) = Z
          ENDDO
          ENDDO
          
          CALL ZSLGEMM('L','t', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatL))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLGEMM left transpose"
            write(*,*) sum(abs(Mat-MatL))/dble(M*M2), " is larger then 1E-14!"
            stop 7
          endif
          
          !test left-mult op='c'
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatL(:,:)=Mat(:,:)
          DO I=1,N
          DO J=1,M2
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + conjg(A(k,i))*Mat(P(K),J) 
          ENDDO
          MatL(P(I),J) = Z
          ENDDO
          ENDDO
          
          CALL ZSLGEMM('L','c', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatL))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLGEMM left conjg-transpose"
            write(*,*) sum(abs(Mat-MatL))/dble(M*M2), " is larger then 1E-14!"
            stop 8
          endif
          
          !test right-mult op='n'
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatR(:,:)=Mat(:,:)
          DO I=1,M
          DO J=1,N
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + Mat(I,P(K)) * A(k,j)
          ENDDO
          MatR(I,P(J)) = Z
          ENDDO
          ENDDO
          
          CALL ZSLGEMM('R','n', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatR))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLGEMM right no-transpose"
            write(*,*) sum(abs(Mat-MatR))/dble(M*M2), " is larger then 1E-14!"
            stop 9
          endif
          
          !test right-mult op='t'
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatR(:,:)=Mat(:,:)
          DO I=1,M
          DO J=1,N
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + Mat(I,P(K)) * A(j,k)
          ENDDO
          MatR(I,P(J)) = Z
          ENDDO
          ENDDO
          
          CALL ZSLGEMM('R','T', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatR))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLGEMM right transpose"
            write(*,*) sum(abs(Mat-MatR))/dble(M*M2), " is larger then 1E-14!"
            stop 10
          endif
          
          !test right-mult op='c'
          DO I=1,M
          DO J=1,M2
            Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
          ENDDO
          ENDDO
          MatR(:,:)=Mat(:,:)
          DO I=1,M
          DO J=1,N
          Z = CMPLX(0.d0,0.d0,kind(0.d0))
          DO K=1,N
            Z = Z + Mat(I,P(K)) * conjg(A(j,k))
          ENDDO
          MatR(I,P(J)) = Z
          ENDDO
          ENDDO
          
          CALL ZSLGEMM('R','c', N, M, M2, A, P, Mat)
          IF(sum(abs(Mat-MatR))/dble(M*M2)>1E-14) then
            write(*,*) "ERROR in ZSLGEMM right conjg-transpose"
            write(*,*) sum(abs(Mat-MatR))/dble(M*M2), " is larger then 1E-14!"
            stop 11
          endif
          
          Deallocate(A,P)
        enddo ! N loop
        Deallocate(Mat,MatL,MatR)
        write(*,*) "SUCCESS"
        
END PROGRAM SLMATMULTS
