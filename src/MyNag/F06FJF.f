      SUBROUTINE F06FJF( N, X, INCX, SCALE, SUMSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FJF returns the values scl and smsq such that
C
C     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
C  to be at least unity and the value of smsq will then satisfy
C
C     1.0 .le. smsq .le. ( sumsq + n ) .
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( x( i ) ) ) .
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABSXI
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ +       ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
C
C     End of F06FJF. ( SSSQ )
C
      END
