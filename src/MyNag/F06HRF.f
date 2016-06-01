      SUBROUTINE F06HRF( N, ALPHA, X, INCX, TOL, THETA )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, THETA
      DOUBLE PRECISION   TOL
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06HRF generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   conjg( P' )*P = I,  aimag( beta ) = 0.0.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - gamma*( zeta )*( zeta  conjg( z' ) ),
C                   (   z  )
C
C  where z is an n element vector, gamma is a scalar such that
C
C     real ( gamma ) = 1.0,
C     aimag( gamma ) = aimag( alpha )/( beta - real( alpha ) )
C
C  and zeta is a real scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  Note that when alpha is real then gamma = 1.0.
C
C  gamma and zeta are returned in THETA as
C
C     THETA = ( zeta, aimag( gamma ) )
C
C  unless x is such that
C
C     max( abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ) .le.
C     max( tol, eps*max( abs( real( alpha  ) ),
C                        abs( aimag( alpha  ) ) ) ),
C
C  where eps is the relative machine precision and tol is the user
C  supplied tolerance  TOL,  in which case  THETA is returned as 0.0, or
C  THETA  is such that  real( THETA ) .le. 0.0, in which case  P  can be
C  taken to be
C
C     P = I              when   THETA = 0.0,
C
C     P = ( THETA  0 )   when   real( THETA ) .le. 0.0, THETA .ne. 0.0.
C         (   0    I )
C
C  beta is overwritten on alpha with the imaginary part of alpha set to
C  zero and z is overwritten on x.
C
C  The routine may be called with  n = 0.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         GAMMA
      DOUBLE PRECISION   BETA, EPS, SCALE, SSQ, ZETA
      LOGICAL            FIRST
C     .. Local Arrays ..
      COMPLEX*16         WORK( 1 )
C     .. External Functions ..
      DOUBLE PRECISION   X02AJF
      EXTERNAL           X02AJF
C     .. External Subroutines ..
      EXTERNAL           ZSCAL, ZDSCAL, F06KJF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DIMAG, DCMPLX, DCONJG, MAX, DREAL, SIGN,
     $                   SQRT
C     .. Save statement ..
      SAVE               EPS, FIRST
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         IF( DIMAG( ALPHA ).EQ.DREAL( ZERO ) )THEN
            THETA = ZERO
         ELSE
            BETA  = -SIGN  ( ABS( ALPHA ), DREAL( ALPHA ) )
            THETA =  DCONJG( ALPHA )/BETA
            ALPHA =  BETA
         END IF
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         IF( DIMAG( ALPHA ).EQ.DREAL( ZERO ) )THEN
            THETA = ZERO
         ELSE
            BETA  = -SIGN  ( ABS( ALPHA ), DREAL( ALPHA ) )
            THETA =  DCONJG( ALPHA )/BETA
            ALPHA =  BETA
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            EPS   =  X02AJF( )
         END IF
C
         SSQ   = ONE
         SCALE = DREAL( ZERO )
         CALL F06KJF( N, X, INCX, SCALE, SSQ )
C
C        Treat cases where  SCALE = zero,  SCALE is negligible
C        and  ALPHA = zero  specially.
C        Note that
C        SCALE = max( abs( real( X( i ) ) ), abs( aimag( X( i ) ) ) ).
C
         IF( ( SCALE.EQ.DREAL( ZERO ) ).OR.
     $       ( SCALE.LE.MAX( TOL,
     $                       EPS*MAX( ABS( DREAL ( ALPHA ) ),
     $                                ABS( DIMAG( ALPHA ) )  ) ) ) )THEN
            IF( DIMAG( ALPHA ).EQ.DREAL( ZERO ) )THEN
               THETA =  ZERO
            ELSE
               BETA  = -SIGN  ( ABS( ALPHA ), DREAL( ALPHA ) )
               THETA =  DCONJG( ALPHA )/BETA
               ALPHA =  BETA
            END IF
         ELSE IF( ALPHA.EQ.ZERO )THEN
            THETA = ONE
            BETA  = SCALE*SQRT( SSQ )
            CALL ZDSCAL( N, -1/BETA, X, INCX )
            ALPHA = BETA
         ELSE
            WORK( 1 ) = ALPHA
            CALL F06KJF( 1, WORK, 1, SCALE, SSQ )
            BETA = SCALE*SQRT( SSQ )
            ZETA = SQRT( ( BETA + ABS( DREAL( ALPHA ) ) )/BETA )
            IF( DREAL( ALPHA ).GT.DREAL( ZERO ) )
     $         BETA = -BETA
            IF( DIMAG( ALPHA ).EQ.DREAL( ZERO ) )THEN
               CALL ZDSCAL( N, -1/( ZETA*BETA ), X, INCX )
               THETA = ZETA
               ALPHA = BETA
            ELSE
               GAMMA = DCMPLX( ONE,
     $                         DIMAG( ALPHA )/( BETA - DREAL( ALPHA )))
               CALL ZSCAL( N, -1/( DCONJG( GAMMA )*ZETA*BETA ),
     $                     X, INCX )
               THETA = DCMPLX( ZETA, DIMAG( GAMMA ) )
               ALPHA = BETA
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06HRF. ( CGRFG )
C
      END
