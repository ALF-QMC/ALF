      SUBROUTINE F06HBF( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         CONST
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06HBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06HBF. ( CLOAD )
C
      END
