      SUBROUTINE F06AAZ ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 15 REVISED. IER-915 (APR 1991).
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IERR, IFAIL
      CHARACTER*4        VARBNM
C     .. Local Arrays ..
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ACF
      EXTERNAL           P01ACF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IF (SRNAME(1:3).EQ.'F06') THEN
         IERR = -1
         VARBNM = '    '
      ELSE
         IERR = -INFO
         VARBNM = 'INFO'
      END IF
      IFAIL = 0
      IFAIL = P01ACF (IFAIL, IERR, SRNAME(1:6), VARBNM, 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
