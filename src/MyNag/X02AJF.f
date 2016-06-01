      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
C     .. Local Scalars ..
      DOUBLE PRECISION Z
      DATA            Z/0.222044604925031336E-15/
C     .. Executable Statements ..
      X02AJF = Z
      RETURN
      END
