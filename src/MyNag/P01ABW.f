      SUBROUTINE P01ABW(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABW increases the value of IERR by 1 and, if
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
C        IF( N.NE.'Valid value' )
C     $     CALL P01ABW( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM
      CHARACTER*(*)     N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(3)
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
         CALL X04BAF(NERR,REC(3))
      END IF
      RETURN
C
C
C     End of P01ABW.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is',/8X,A)
      END
