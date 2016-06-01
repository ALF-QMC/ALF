      SUBROUTINE F06KJF( N, X, INCX, SCALE, SUMSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06KJF returns the values scl and ssq such that
C
C     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
C  assumed to be at least unity and the value of ssq will then satisfy
C
C     1.0 .le. ssq .le. ( sumsq + 2*n ).
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
C            i
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  SCALE and SUMSQ are overwritten by scl and ssq respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 basic linear algebra routine.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-April-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DIMAG, DREAL
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( DREAL( X( IX ) ).NE.ZERO )THEN
               TEMP1 = ABS( DREAL( X( IX ) ) )
               IF( SCALE.LT.TEMP1 )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ +       ( TEMP1/SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO )THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ +       ( TEMP1/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
C
      RETURN
C
C     End of F06KJF. ( SCSSQ )
C
      END
