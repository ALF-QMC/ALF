      SUBROUTINE F06FRF( N, ALPHA, X, INCX, TOL, ZETA )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, TOL, ZETA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FRF generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA, EPS, SCALE, SSQ
      LOGICAL            FIRST
C     .. External Functions ..
      DOUBLE PRECISION   X02AJF
      EXTERNAL           X02AJF
C     .. External Subroutines ..
      EXTERNAL           F06FJF, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
C     .. Save statement ..
      SAVE               EPS, FIRST
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         ZETA = ZERO
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         ZETA = ZERO
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            EPS   =  X02AJF( )
         END IF
C
C        Treat case where P is a 2 by 2 matrix specially.
C
         IF( N.EQ.1 )THEN
C
C           Deal with cases where  ALPHA = zero  and
C           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
C
            IF( ALPHA.EQ.ZERO )THEN
               ZETA   =  ONE
               ALPHA  =  ABS ( X( 1 ) )
               X( 1 ) = -SIGN( ONE, X( 1 ) )
            ELSE IF( ABS( X( 1 ) ).LE.MAX( EPS*ABS( ALPHA ), TOL ) )THEN
               ZETA   =  ZERO
            ELSE
               IF( ABS( ALPHA ).GE.ABS( X( 1 ) ) )THEN
                  BETA = ABS( ALPHA ) *SQRT( 1 + ( X( 1 )/ALPHA )**2 )
               ELSE
                  BETA = ABS( X( 1 ) )*SQRT( 1 + ( ALPHA/X( 1 ) )**2 )
               END IF
               ZETA = SQRT( ( ABS( ALPHA ) + BETA )/BETA )
               IF( ALPHA.GE.ZERO )
     $            BETA = -BETA
               X( 1 ) = -X( 1 )/( ZETA*BETA )
               ALPHA  = BETA
            END IF
         ELSE
C
C           Now P is larger than 2 by 2.
C
            SSQ   = ONE
            SCALE = ZERO
            CALL F06FJF( N, X, INCX, SCALE, SSQ )
C
C           Treat cases where  SCALE = zero,
C           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
C           ALPHA = zero  specially.
C           Note that  SCALE = max( abs( X( i ) ) ).
C
            IF( ( SCALE.EQ.ZERO ).OR.
     $          ( SCALE.LE.MAX( EPS*ABS( ALPHA ), TOL ) ) )THEN
               ZETA  = ZERO
            ELSE IF( ALPHA.EQ.ZERO )THEN
               ZETA  = ONE
               ALPHA = SCALE*SQRT( SSQ )
               CALL DSCAL( N, -1/ALPHA, X, INCX )
            ELSE
               IF( SCALE.LT.ABS( ALPHA ) )THEN
                  BETA = ABS( ALPHA )*SQRT( 1 + SSQ*( SCALE/ALPHA )**2 )
               ELSE
                  BETA = SCALE       *SQRT( SSQ +   ( ALPHA/SCALE )**2 )
               END IF
               ZETA = SQRT( ( BETA + ABS( ALPHA ) )/BETA )
               IF( ALPHA.GT.ZERO )
     $            BETA = -BETA
               CALL DSCAL( N, -1/( ZETA*BETA ), X, INCX )
               ALPHA = BETA
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FRF. ( SGRFG )
C
      END
