      SUBROUTINE F01QCF(M,N,A,LDA,ZETA,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QCF  finds  the  QR factorization  of the real  m by n,  m .ge. n,
C  matrix A,  so that  A is reduced to upper triangular form by means of
C  orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )   when   m.gt.n,
C           ( 0 )
C
C     A = Q*R       when   m = n,
C
C  where  Q  is an  m by m orthogonal matrix and  R  is an  n by n upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used  to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k ) and z( k )  are chosen to annhilate the elements  below the
C  triangular part of  A.
C
C  The vector  u( k ) is returned in the kth element of  ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of  z( k ) are in  a( k + 1, k ), ..., a( m, k ).  The elements of  R
C  are returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( n )*Q( n - 1 )*...*Q( 1 ) )'.
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at least  n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be  at  least zero. When  N = 0  then an immediate return is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the upper triangular matrix R and the  M by N strictly lower
C           triangular  part   of   A   will  contain  details   of  the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of DIMENSION at least ( n ).
C
C           On exit,  ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then  ZETA( k ) = 0.0,
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and  zeta( k ) is always in the range  ( 1.0, sqrt( 2.0 ) ).
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On successful  exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to  -1  indicating that an  input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        M   .lt. N
C        N   .lt. 0
C        LDA .lt. M
C
C  If  on  entry,  IFAIL  was  either  -1 or 0  then  further diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C        B := Q*B   and   B := Q'*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  NAG Library routine  F01QDF. The  operation  B := Q*B can be obtained
C  by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  and  B := Q'*B  can be obtained by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine F01QDF for further
C  details.
C
C  The first k columns of the orthogonal matrix Q can either be obtained
C  by setting  B to the first k columns of the unit matrix and using the
C  first of the above two calls,  or by calling the  NAG Library routine
C  F01QEF, which overwrites the k columns of Q on the first k columns of
C  the array A.  Q is obtained by the call:
C
C     CALL F01QEF( 'Separate', M, N, K, A, LDA, ZETA, WORK, IFAIL )
C
C  As above WORK must be a k element array.  If K is larger than N, then
C  A must have been declared to have at least K columns.
C
C  Operations involving the matrix  R  can readily  be performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV  (see Chapter F06), but note
C  that no test for  near singularity  of  R  is incorporated in DTRSV .
C  If  R  is singular,  or nearly singular then the  NAG Library routine
C  F02WUF  can be  used to  determine  the  singular value decomposition
C  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 21-December-1985.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ZETA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           IERR, K, LA
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, F06FRF, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF (M.LT.N) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the factorization.
C
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      LA = LDA
      DO 20 K = 1, MIN(M-1,N)
C
C        Use a  Householder reflection  to  zero the  kth column  of  A.
C        First set up the reflection.
C
         CALL F06FRF(M-K,A(K,K),A(K+1,K),1,ZERO,ZETA(K))
         IF ((ZETA(K).GT.ZERO) .AND. (K.LT.N)) THEN
            IF ((K+1).EQ.N) LA = M - K + 1
C
C           Temporarily  store  beta and  put  zeta( k )  in  a( k, k ).
C
            TEMP = A(K,K)
            A(K,K) = ZETA(K)
C
C           We now perform the operation  A := Q( k )*A.
C
C           Let  B  denote  the bottom  ( m - k + 1 ) by ( n - k )  part
C           of  A.
C
C           First form   work = B'*u.  ( work  is stored in the elements
C           ZETA( k + 1 ), ..., ZETA( n ). )
C
            CALL DGEMV('Transpose',M-K+1,N-K,ONE,A(K,K+1),LA,A(K,K),1,
     *                 ZERO,ZETA(K+1),1)
C
C           Now form  B := B - u*work'.
C
            CALL DGER(M-K+1,N-K,-ONE,A(K,K),1,ZETA(K+1),1,A(K,K+1),LA)
C
C           Restore beta.
C
            A(K,K) = TEMP
         END IF
   20 CONTINUE
C
C     Set the final  ZETA  when  m.eq.n.
C
      IF (M.EQ.N) ZETA(N) = ZERO
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QCF. ( SGEQR )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
