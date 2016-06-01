      SUBROUTINE F01QDF(TRANS,WHERET,M,N,A,LDA,ZETA,NCOLB,B,LDB,WORK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QDF performs one of the transformations
C
C     B := Q*B   or   B := Q'*B,
C
C  where  B is an  m by ncolb real matrix and  Q is an m by m orthogonal
C  matrix, given as the product of  Householder transformation matrices.
C
C  This  routine  is  intended  for  use  following  NAG Fortran Library
C  routine  F01QCF.
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = ( Q( n )*Q( n - 1 )*...*Q( 1 ) )',
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )'
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  a( k + 1, k ), ..., a( m, k )  and  zeta( k ) must be supplied either
C  in  a( k, k ) or in zeta( k ),  depending upon the parameter  WHERET.
C
C  To obtain Q explicitly B may be set to I and premultiplied by Q. This
C  is more efficient than obtaining Q'.
C
C  3. Parameters
C     ==========
C
C  TRANS  - CHARACTER*1.
C
C           On entry, TRANS  specifies the operation to be performed  as
C           follows.
C
C           TRANS = 'N' or 'n'  ( No transpose )
C
C              Perform the operation  B := Q*B.
C
C           TRANS = 'T' or 't' or 'C' or 'c'  ( Transpose )
C
C              Perform the operation  B := Q'*B.
C
C           Unchanged on exit.
C
C  WHERET - CHARACTER*1.
C
C           On entry,  WHERET  specifies where the elements of  zeta are
C           to be found as follows.
C
C           WHERET = 'I' or 'i'   ( In A )
C
C              The elements of zeta are in A.
C
C           WHERET = 'S' or 's'   ( Separate )
C
C              The elements of zeta are separate from A, in ZETA.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHERET = 'I' or 'i'  then  the  diagonal
C           elements of A must contain the elements of zeta as described
C           under the argument  ZETA  below.
C
C           When  WHERET = 'S' or 's'  then the diagonal elements of the
C           array  A  are referenced, since they are used temporarily to
C           store the  zeta( k ), but they contain their original values
C           on return.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least m.
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of  DIMENSION  at least  ( n ),  when
C           WHERET = 'S' or 's'.
C
C           Before entry with  WHERET = 'S' or 's', the array  ZETA must
C           contain  the  elements  of  zeta.  If  ZETA( k ) = 0.0  then
C           T( k )  is assumed  to be  I otherwise  ZETA( k ) is assumed
C           to contain zeta( k ).
C
C           When WHERET = 'I' or 'i', the array  ZETA is not referenced.
C
C           Unchanged on exit.
C
C  NCOLB  - INTEGER.
C
C           On  entry, NCOLB  must specify  the number of columns of  B.
C           NCOLB  must  be  at  least  zero.  When  NCOLB = 0  then  an
C           immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL             array of DIMENSION ( LDB, ncolb ).
C
C           Before entry, the leading  M by NCOLB  part of  the array  B
C           must  contain  the matrix to be  transformed.
C
C           On  exit,  B  is  overwritten  by  the  transformed  matrix.
C
C  LDB    - INTEGER.
C
C           On  entry, LDB  must specify  the  leading dimension of  the
C           array  B as declared in the calling (sub) program. LDB  must
C           be at least m.
C
C           Unchanged on exit.
C
C  WORK   - REAL             array of DIMENSION at least ( ncolb ).
C
C           Used as internal workspace.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to   -1  indicating that an input parameter has
C           been  incorrectly  set. See  the  next  section  for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        TRANS  .ne. 'N' or 'n' or 'T' or 't' or 'C' or 'c'
C        WHERET .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        LDA    .lt. M
C        NCOLB  .lt. 0
C        LDB    .lt. M
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QDF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDB, M, N, NCOLB
      CHARACTER*1       TRANS, WHERET
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), WORK(*), ZETA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, ZETAK
      INTEGER           IERR, K, KK, LB
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF ((TRANS.NE.'N') .AND. (TRANS.NE.'n') .AND. (TRANS.NE.'T')
     *    .AND. (TRANS.NE.'t') .AND. (TRANS.NE.'C') .AND. (TRANS.NE.'c')
     *    ) CALL P01ABW(TRANS,'TRANS',IFAIL,IERR,SRNAME)
      IF ((WHERET.NE.'I') .AND. (WHERET.NE.'i') .AND. (WHERET.NE.'S')
     *     .AND. (WHERET.NE.'s')) CALL P01ABW(WHERET,'WHERET',IFAIL,
     *    IERR,SRNAME)
      IF (M.LT.N) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (NCOLB.LT.0) CALL P01ABY(NCOLB,'NCOLB',IFAIL,IERR,SRNAME)
      IF (LDB.LT.M) CALL P01ABY(LDB,'LDB',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the transformation.
C
      IF (MIN(N,NCOLB).EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      LB = LDB
      DO 20 KK = 1, N
         IF ((TRANS.EQ.'T') .OR. (TRANS.EQ.'t') .OR. (TRANS.EQ.'C')
     *       .OR. (TRANS.EQ.'c')) THEN
C
C           Q'*B = Q( n )*...*Q( 2 )*Q( 1 )*B,
C
            K = KK
         ELSE
C
C           Q*B  = Q( 1 )'*Q( 2 )'*...*Q( n )'*B,
C
            K = N + 1 - KK
         END IF
         IF ((WHERET.EQ.'S') .OR. (WHERET.EQ.'s')) THEN
            ZETAK = ZETA(K)
         ELSE
            ZETAK = A(K,K)
         END IF
         IF (ZETAK.GT.ZERO) THEN
            TEMP = A(K,K)
            A(K,K) = ZETAK
            IF (NCOLB.EQ.1) LB = M - K + 1
C
C           Let C denote the bottom ( m - k + 1 ) by ncolb part of B.
C
C           First form  work = C'*u.
C
            CALL DGEMV('Transpose',M-K+1,NCOLB,ONE,B(K,1),LB,A(K,K),1,
     *                 ZERO,WORK,1)
C
C           Now form  C := C - u*work'.
C
            CALL DGER(M-K+1,NCOLB,-ONE,A(K,K),1,WORK,1,B(K,1),LB)
C
C           Restore the diagonal element of A.
C
            A(K,K) = TEMP
         END IF
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QDF. ( SGEAPQ )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
