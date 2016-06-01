      SUBROUTINE F01REF(WHERET,M,N,NCOLQ,A,LDA,THETA,WORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-733 (DEC 1989).
C     MARK 14C REVISED. IER-886 (NOV 1990).
C
C  1. Purpose
C     =======
C
C  F01REF  returns the first  ncolq columns of the m by m unitary matrix
C  Q,  where  Q  is given  as the product of  Householder transformation
C  matrices.
C
C  This  routine  is  intended  for  use  following  NAG Fortran Library
C  routine F01RCF.
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = conjg( ( Q( n )*Q( n - 1 )*...*Q( 1 ) )' ),
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - gamma( k )*u( k )*conjg( u( k )' )
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  gamma( k ) is a scalar for which  real( gamma( k ) ) = 1.0, zeta( k )
C  is a real scalar and z( k ) is an ( m - k ) element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  a( k + 1, k ), ..., a( m, k ) and theta( k ), given by
C
C     theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ),
C
C  must be supplied either in a( k, k ) or in theta( k ), depending upon
C  the parameter WHERET.
C
C  3. Parameters
C     ==========
C
C  WHERET - CHARACTER*1.
C
C           On entry, WHERET  specifies where the elements of  theta are
C           to be found as follows.
C
C           WHERET = 'I' or 'i'   ( In A )
C
C              The elements of theta are in A.
C
C           WHERET = 'S' or 's'   ( Separate )
C
C              The elements of theta are separate from A, in THETA.
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
C           be at least zero.
C
C           Unchanged on exit.
C
C  NCOLQ  - INTEGER.
C
C           On entry, NCOLQ  must specify the required number of columns
C           of Q.  NCOLQ must be at least zero and not be larger than m.
C           When   NCOLQ = 0  then  an  immediate  return  is  effected.
C
C           Unchanged on exit.
C
C  A      - COMPLEX array of DIMENSION ( LDA, nca ),  where  nca must be
C           at least  max( n, ncolq ).
C
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHERET = 'I' or 'i'  then  the  diagonal
C           elements  of  A  must  contain  the  elements  of  theta  as
C           described under the argument  THETA  below.
C
C           On  exit, the  first  NCOLQ  columns  of  the  array  A  are
C           overwritten by the first ncolq columns of the m by m unitary
C           matrix Q.
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
C  THETA  - COMPLEX array of DIMENSION at least ( n ), when WHERET = 'S'
C           or 's'.
C
C           Before entry with  WHERET = 'S' or 's', the array THETA must
C           contain  the elements of  theta.  If  THETA( k ) = 0.0  then
C           T( k )  is assumed  to be  I,  if  THETA( k ) = alpha,  with
C           real( alpha ) .lt. 0.0  then  T( k )  is assumed  to  be  of
C           the form
C
C              T( k ) = ( alpha  0 ),
C                       (   0    I )
C
C           otherwise  THETA( k ) is assumed to contain theta( k ) given
C           by  theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ).
C
C           When WHERET = 'I' or 'i', the array THETA is not referenced.
C
C           Unchanged on exit.
C
C  WORK   - COMPLEX array of DIMENSION at least ( ncolq ).
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
C        WHERET .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        NCOLQ  .lt. 0  .or.  NCOLQ .gt. M
C        LDA    .lt. M
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
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01REF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N, NCOLQ
      CHARACTER*1       WHERET
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), THETA(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        GAMMA, THETAK
      INTEGER           IERR, K, NCQ, P
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERC, ZSCAL, F06HBF, F06THF, P01ABW,
     *                  P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, DCMPLX, DCONJG, MIN, DREAL
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IF (NCOLQ.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      IERR = 0
      IF ((WHERET.NE.'I') .AND. (WHERET.NE.'i') .AND. (WHERET.NE.'S')
     *     .AND. (WHERET.NE.'s')) CALL P01ABW(WHERET,'WHERET',IFAIL,
     *    IERR,SRNAME)
      IF (M.LT.N) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF ((NCOLQ.LT.0) .OR. (NCOLQ.GT.M)) CALL P01ABY(NCOLQ,'NCOLQ',
     *    IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Start to form Q. First set the elements above the leading diagonal
C     to zero.
C
      P = MIN(N,NCOLQ)
      IF (P.GT.1) CALL F06THF('Upper',P-1,P-1,ZERO,ZERO,A(1,2),LDA)
      IF (NCOLQ.GT.N) THEN
         NCQ = NCOLQ - N
C
C        Set the last  ( ncolq - n ) columns of  Q  to those of the unit
C        matrix.
C
         CALL F06THF('General',N,NCOLQ-N,ZERO,ZERO,A(1,N+1),LDA)
         CALL F06THF('General',M-N,NCOLQ-N,ZERO,ONE,A(N+1,N+1),LDA)
      ELSE
         NCQ = 0
      END IF
      DO 20 K = P, 1, -1
C
C        Q*E( ncolq ) =
C           conjg( Q( 1 )' )*...*conjg( Q( p )' )*E( ncolq ),
C        where  E( ncolq )  is the  matrix  containing  the first  ncolq
C        columns of  I.
C
         IF ((WHERET.EQ.'S') .OR. (WHERET.EQ.'s')) THEN
            THETAK = THETA(K)
         ELSE
            THETAK = A(K,K)
         END IF
C
C        If  real( THETA( k ) ) .le. zero  then Q( k ) is special.
C
         IF (DREAL(THETAK).GT.DREAL(ZERO)) THEN
            A(K,K) = DREAL(THETAK)
            GAMMA = DCMPLX(DREAL(ONE),-DIMAG(THETAK))
C
C           Let C denote the bottom ( m - k + 1 ) by ncq part of Q.
C
C           First form  work = conjg( C' )*u.
C
            IF ((K.LT.M) .AND. (NCQ.GT.0)) THEN
               CALL ZGEMV('Conjugate',M-K+1,NCQ,ONE,A(K,K+1),LDA,A(K,K),
     *                    1,ZERO,WORK,1)
C
C              Now form  C := C - gamma( k )*u*conjg( work' ).
C
               CALL ZGERC(M-K+1,NCQ,-GAMMA,A(K,K),1,WORK,1,A(K,K+1),LDA)
            END IF
C
C           Now form the kth column of Q.
C
            CALL ZSCAL(M-K+1,-GAMMA*DREAL(THETAK),A(K,K),1)
            A(K,K) = ONE + A(K,K)
         ELSE
            IF (DIMAG(THETAK).EQ.DREAL(ZERO)) THEN
               A(K,K) = ONE
            ELSE
               A(K,K) = DCONJG(THETAK)
            END IF
            IF (K.LT.M) CALL F06HBF(M-K,ZERO,A(K+1,K),1)
         END IF
         NCQ = NCQ + 1
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C     End of F01REF. ( CGEFQ  )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
