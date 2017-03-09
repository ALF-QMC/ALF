   Module UDV_Wrap_mod
     Use MyMats
     Use Files_mod
     
   Contains

!***************************************************************
     Subroutine UDV_Wrap_Pivot(A,U,D,V,NCON,N1,N2)
       
       Implicit NONE
       COMPLEX (KIND=8), INTENT(IN),    DIMENSION(:,:) :: A
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:) :: D
       INTEGER, INTENT(IN) :: NCON
       INTEGER, INTENT(IN) :: N1,N2
       
       ! Locals
       REAL (Kind=8) :: VHELP(N2), XNORM(N2), XMAX, XMEAN
       INTEGER :: IVPT(N2), IVPTM1(N2), I, J, K, IMAX
       COMPLEX (KIND=8)  :: A1(N1,N2), A2(N1,N2), V1(N2,N2), U1(N2,N2), Z
       
       DO I = 1,N2
          XNORM(I) = 0.D0
          DO J = 1,N1
             XNORM(I) = XNORM(I) + DBLE( A(J,I) * CONJG( A(J,I) ) )
          ENDDO
       ENDDO
       VHELP = XNORM
      
       DO I = 1,N2
          XMAX = VHELP(1)
          IMAX = 1
          DO J = 2,N2
             IF (VHELP(J).GT.XMAX)  THEN
                IMAX = J
                XMAX = VHELP(J)
             ENDIF
          ENDDO
          VHELP(IMAX) = -1.D0
          IVPTM1(IMAX)=  I
          IVPT(I) = IMAX
       ENDDO
       DO I = 1,N2
          K = IVPT(I)
          A1(:, I) = A(:, K)/CMPLX(XNORM(K),0.d0,Kind(0.d0))
       ENDDO
       
       CALL UDV(A1,U,D,V1,NCON)
       

       ! Finish the pivotting.
       DO I = 1,N2
          D(I) = D(I)*XNORM(IVPT(I))
       ENDDO
       DO I = 1,N2-1
          Z = 1.D0/XNORM(IVPT(I))
          DO J = I+1,N2
             V1(I,J) = V1(I,J)*XNORM(IVPT(J))*Z
          ENDDO
       ENDDO

       DO J = 1,N2
          DO I = 1,N2
             V(I,J) = V1(I,IVPTM1(J))
          ENDDO
       ENDDO
       

       IF (NCON == 1) THEN
          !Check the result  A = U D V
          DO J = 1,N2
             A1(:, J) = D * V(:, J)
          ENDDO
          !Write(6,*) 'Here'
          Call MMULT (A2,U,A1)
          !Call MMULT (A2,U,V)
          CALL COMPARE(A,A2,XMAX,XMEAN)
          Write (6,*) 'Check afer  Pivoting', XMAX
       ENDIF

       
       
     End Subroutine UDV_Wrap_Pivot

!***************************************************************

     Subroutine UDV_Wrap_Pivot2(A,U,D,V,NCON,N1,N2)
       
       Implicit NONE
       COMPLEX (KIND=8), INTENT(IN),    DIMENSION(:,:) :: A
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:) :: D
       INTEGER, INTENT(IN) :: NCON
       INTEGER, INTENT(IN) :: N1,N2
       
       ! Locals
       REAL (Kind=8) :: VHELP(N2), XNORM(N2), XMAX, XMEAN
       INTEGER :: IVPT(N2), IVPTM1(N2), I, J, K, IMAX
       COMPLEX (KIND=8)  :: A1(N1,N2), A2(N1,N2), V1(N2,N2), U1(N2,N2)
       
       DO I = 1,N2
          XNORM(I) = 0.D0
          DO J = 1,N1
             XNORM(I) = XNORM(I) + DBLE( A(J,I) * CONJG( A(J,I) ) )
          ENDDO
       ENDDO
       VHELP = XNORM
      
       DO I = 1,N2
          XMAX = VHELP(1)
          IMAX = 1
          DO J = 2,N2
             IF (VHELP(J).GT.XMAX)  THEN
                IMAX = J
                XMAX = VHELP(J)
             ENDIF
          ENDDO
          VHELP(IMAX) = -1.D0
          IVPTM1(IMAX)=  I
          IVPT(I) = IMAX
       ENDDO
       DO I = 1,N2
          K = IVPT(I)
          A1(:, I) = A(:, K)/CMPLX(XNORM(K),0.d0,Kind(0.d0))
       ENDDO
       

       CALL QR(A1,U,V1,NCON)

       V = cmplx(0.d0,0.d0,Kind(0.d0))
       DO I = 1,N2
          K = IVPTM1(I)
          Do J = 1,K
             V(J, I) = V1(J,K)*CMPLX(XNORM(I),0.d0,Kind(0.d0))
          ENDDO
       ENDDO
       
       A1 = V
       CALL SVD(A1,U1,D,V,NCON)
       Call MMULT(A1,U,U1)
       U = A1

       IF (NCON == 1) THEN
          !Check the result  A = U D V
          DO J = 1,N2
             A1(:, J) = D * V(:, J)
          ENDDO
          !Write(6,*) 'Here'
          Call MMULT (A2,U,A1)
          !Call MMULT (A2,U,V)
          CALL COMPARE(A,A2,XMAX,XMEAN)
          Write (6,*) 'Check afer  Pivoting', XMAX
       ENDIF

       
       
     End Subroutine UDV_Wrap_Pivot2

!***************************************************************

     Subroutine UDV_Wrap_Pivot1(A,U,D,V,NCON,N1,N2)
       
       Implicit NONE
       COMPLEX (KIND=8), INTENT(IN),    DIMENSION(:,:) :: A
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:) :: D
       INTEGER, INTENT(IN) :: NCON
       INTEGER, INTENT(IN) :: N1,N2
       
       ! Locals
       REAL (Kind=8) :: VHELP(N2), XNORM(N2), XMAX, XMEAN
       INTEGER :: IVPT(N2), IVPTM1(N2), I, J, K, IMAX
       COMPLEX (KIND=8)  :: A1(N1,N2), A2(N1,N2), V1(N2,N2)
       
       DO I = 1,N2
          XNORM(I) = 0.D0
          DO J = 1,N1
             XNORM(I) = XNORM(I) + DBLE( A(J,I) * CONJG( A(J,I) ) )
          ENDDO
       ENDDO
       VHELP = XNORM
      
       DO I = 1,N2
          XMAX = VHELP(1)
          IMAX = 1
          DO J = 2,N2
             IF (VHELP(J).GT.XMAX)  THEN
                IMAX = J
                XMAX = VHELP(J)
             ENDIF
          ENDDO
          VHELP(IMAX) = -1.D0
          IVPTM1(IMAX)=  I
          IVPT(I) = IMAX
       ENDDO
       DO I = 1,N2
          K = IVPT(I)
          A1(:, I) = A(:, K)
       ENDDO
       
       CALL UDV_Wrap(A1,U,D,V,NCON)
       
       V1 = V
       DO I = 1,N2
          K = IVPTM1(I)
          V(:, I) = V1(:,K)
       ENDDO
       

       IF (NCON == 1) THEN
          !Check the result  A = U D V
          DO J = 1,N2
             A1(:, J) = D * V(:, J)
          ENDDO
          Write(6,*) 'Here'
          Call MMULT (A2,U,A1)
          CALL COMPARE(A,A2,XMAX,XMEAN)
          Write (6,*) 'Check afer  Pivoting', XMAX
       ENDIF

       
       
     End Subroutine UDV_Wrap_Pivot1
!***************************************************************
     Subroutine UDV_Wrap(A,U,D,V,NCON)
       Implicit None
#ifdef MPI            
       INCLUDE 'mpif.h'
#endif
       COMPLEX (KIND=8), INTENT(IN),    DIMENSION(:,:) :: A
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
       COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:) :: D
       INTEGER, INTENT(IN) :: NCON
       
       !Local 
       Complex (Kind=8), Allocatable ::  A1(:,:),U1(:,:)
       Integer :: N,I,J
       character (len=64) :: file_sr, File
#ifdef MPI  
       INTEGER :: STATUS(MPI_STATUS_SIZE)
       INTEGER :: Isize, Irank,Ierr
            
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
   
       File_sr = "SDV"
#ifdef MPI 
       File = File_i(File_sr, Irank)
#else
       File = File_sr
#endif 
       !Open (Unit = 78,File=File, Status='UNKNOWN', action="write", position="append")
       !Write(78,*) 'Call QR'
       !Close(78)
       CALL QR(A,U,V,NCON)
       !Open (Unit = 78,File=File, Status='UNKNOWN', action="write", position="append")
       !Write(78,*) 'End call QR'
       !Close(78)
       N = Size(A,1)
       Allocate (A1(N,N),U1(N,N))
       DO I = 1,N
          DO J = 1,N
             WRITE(6,*) I,J,V(I,J)
          ENDDO
       ENDDO
       A1 = V
       !Open (Unit = 78,File=File, Status='UNKNOWN')
       !Write(78,*) 'Call SVD'
       !DO I = 1,N
       !   Write(78,*) Real(V(I,I))
       !ENDDO
       !Close(78)
       CALL SVD(A1,U1,D,V,NCON)
       !Open (Unit = 78,File=File, Status='UNKNOWN', action="write", position="append")
       !Write(78,*) 'End call SVD'
       !Close(78)
       Call MMULT(A1,U,U1)
       U = A1
       
     End Subroutine UDV_Wrap
     
   End Module UDV_Wrap_mod
   
