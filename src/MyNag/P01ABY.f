      SUBROUTINE P01ABY(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABY increases the value of IERR by 1 and, if
C
C        ( mod( INFORM, 10 ).ne.1 ).or.( mod( INFORM/10, 10 ).ne.0 )
C
C     writes a message on the current error message channel giving the
C     value of N, a message to say that N is invalid and the strings
C     NAME and SRNAME.
C
C     NAME must be the name of the actual argument for N and SRNAME must
C     be the name of the calling routine.
C
C     This routine is intended for use when N is an invalid input
C     parameter to routine SRNAME. For example
C
C        IERR = 0
C        IF( N.LT.1 )CALL P01ABY( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 23-February-1984.  Sven.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM, N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(2)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = IERR + 1
      IF ((MOD(INFORM,10).NE.1) .OR. (MOD(INFORM/10,10).NE.0)) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) NAME, SRNAME, N
         CALL X04BAF(NERR,' ')
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
      END IF
      RETURN
C
C
C     End of P01ABY.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is ',I6)
      END
