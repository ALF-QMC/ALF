      SUBROUTINE F01RCF(M,N,A,LDA,THETA,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-732 (DEC 1989).
C
C  1. Purpose
C     =======
C
C  F01RCF  finds the  QR factorization of the complex m by n,  m .ge. n,
C  matrix A,  so that  A is reduced to upper triangular form by means of
C  unitary transformations.
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
C  where  Q  is an  m by m  unitary matrix  and  R  is an  n by n  upper
C  triangular matrix with real diagonal elements.
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
C     T( k ) = I - gamma( k )*u( k )*conjg( u( k )' ),
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  gamma( k ) is a scalar for which  real( gamma( k ) ) = 1.0, zeta( k )
C  is  a  real  scalar  and  z( k )  is  an  ( m - k )  element  vector.
C  gamma( k ), zeta( k ) and z( k ) are chosen to annhilate the elements
C  below  the  triangular part of  A  and to make  the diagonal elements
C  real.
C
C  The scalar  gamma( k ) and the vector  u( k ) are returned in the kth
C  element of  THETA and in the kth column of  A, such that  theta( k ),
C  given by
C
C     theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ),
C
C  is in  THETA( k )  and the elements of  z( k ) are in  a( k + 1, k ),
C  ..., a( m, k ).   The  elements  of  R  are  returned  in  the  upper
C  triangular part of  A.
C
C  Q is given by
C
C     Q = conjg( ( Q( n )*Q( n - 1 )*...*Q( 1 ) )' ).
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
C  A      - COMPLEX array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the upper triangular matrix  R, with the  imaginary parts of
C           the diagonal elements set to zero, and the  M by N  strictly
C           lower triangular part of  A  will  contain  details  of  the
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
C  THETA  - COMPLEX array of DIMENSION at least ( n ).
C
C           On exit, THETA( k )  contains the scalar  theta( k ) for the
C           kth  transformation. If  T( k ) = I  then  THETA( k ) = 0.0,
C           if
C
C              T( k ) = ( alpha  0 ),   real( alpha ) .lt. 0.0,
C                       (   0    I )
C
C           then   THETA( k ) = alpha,   otherwise  THETA( k )  contains
C           theta( k )  as  described above  and  real( theta( k ) )  is
C           always in the range ( 1.0, sqrt( 2.0 ) ).
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
C        B := Q*B   and   B := conjg( Q' )*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  NAG Library routine  F01RDF. The  operation  B := Q*B can be obtained
C  by the call:
C
C     IFAIL = 0
C     CALL F01RDF( 'No conjugate', 'Separate', M, N, A, LDA, THETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  and  B := conjg( Q' )*B  can be obtained by the call:
C
C     IFAIL = 0
C     CALL F01RDF( 'Conjugate', 'Separate', M, N, A, LDA, THETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine F01RDF for further
C  details.
C
C  The first  k columns of the unitary matrix  Q  can either be obtained
C  by setting  B to the first k columns of the unit matrix and using the
C  first of the above two calls,  or by calling the  NAG Library routine
C  F01REF, which overwrites the k columns of Q on the first k columns of
C  the array A.  Q is obtained by the call:
C
C     CALL F01REF( 'Separate', M, N, K, A, LDA, THETA, WORK, IFAIL )
C
C  As above WORK must be a k element array.  If K is larger than N, then
C  A must have been declared to have at least K columns.
C
C  Operations involving the matrix  R  can readily  be performed by  the
C  Level 2 BLAS  routines  CTRSV  and CTRMV  (see Chapter F06), but note
C  that no test for  near singularity  of  R  is incorporated in CTRSV .
C  If  R  is singular,  or nearly singular then the  NAG Library routine
C  F02XUF  can be  used to  determine  the  singular value decomposition
C  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 21-December-1985.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01RCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), THETA(*)
C     .. Local Scalars ..
      COMPLEX*16        GAMMA
      DOUBLE PRECISION  TEMP
      INTEGER           IERR, K, LA
C     .. Local Arrays ..
      COMPLEX*16        DUMMY(1)
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERC, ZSCAL, F06HRF, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, DCMPLX, MIN, DREAL
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
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
      LA = LDA
      DO 20 K = 1, MIN(M-1,N)
C
C        Use a  Householder reflection  to  zero the  kth column  of  A.
C        First set up the reflection.
C
         CALL F06HRF(M-K,A(K,K),A(K+1,K),1,DREAL(ZERO),THETA(K))
         IF ((DREAL(THETA(K)).GT.DREAL(ZERO)) .AND. (K.LT.N)) THEN
            IF ((K+1).EQ.N) LA = M - K + 1
C
C           Temporarily store  beta,  put  zeta( k )  in  a( k, k )  and
C           form  gamma( k ).
C
            TEMP = A(K,K)
            A(K,K) = DREAL(THETA(K))
            GAMMA = DCMPLX(DREAL(ONE),DIMAG(THETA(K)))
C
C           We now perform the operation  A := Q( k )*A.
C
C           Let  B  denote  the bottom  ( m - k + 1 ) by ( n - k )  part
C           of  A.
C
C           First form   work = conjg( B' )*u.  ( work  is stored in the
C           elements  THETA( k + 1 ), ..., THETA( n ). )
C
            CALL ZGEMV('Conjugate',M-K+1,N-K,ONE,A(K,K+1),LA,A(K,K),1,
     *                 ZERO,THETA(K+1),1)
C
C           Now form  B := B - gamma( k )*u*conjg( work' ).
C
            CALL ZGERC(M-K+1,N-K,-GAMMA,A(K,K),1,THETA(K+1),1,A(K,K+1),
     *                 LA)
C
C           Restore beta.
C
            A(K,K) = TEMP
         ELSE IF (DIMAG(THETA(K)).NE.DREAL(ZERO)) THEN
            CALL ZSCAL(N-K,THETA(K),A(K,K+1),LDA)
         END IF
   20 CONTINUE
C
C     Find the final  THETA  when  m.eq.n.  This ensures  that the  last
C     diagonal element of  R  is real.
C
      IF (M.EQ.N) CALL F06HRF(0,A(N,N),DUMMY,1,DREAL(ZERO),THETA(N))
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C     End of F01RCF. ( CGEQR )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
