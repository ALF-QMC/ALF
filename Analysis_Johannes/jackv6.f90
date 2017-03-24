      Program enerJ

        
        Use ERRORS
        Implicit none

        Complex (KIND=8), DIMENSION(:,:), ALLOCATABLE :: OBS
        Complex (KIND=8), DIMENSION(:),   ALLOCATABLE :: EN, SIGN
        Complex (KIND=8) :: XM, XERR,x,x1,y,y1

        Complex (Kind=8) Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10, Z11
        Integer :: NST, NS, NS1, NS2, NSTEP, NC, NP, NOBS, Nbins, NP_EFF, ISEED, I, IOBS
        Integer :: N, NBIN

        ! Count the number of bins
        Open (Unit=10, File="ener", status="unknown") 
        !Open (Unit=12, File="ener_hist", status="unknown") 
        nbins = 0
        do
            read(10,*,End=10)  Z1, Z2, Z3, Z4, Z5, Z5, Z5, Z5, Z5, Z5!, Z11
            nbins = nbins + 1
         enddo
10       continue
         Write(6,*) "# of bins: ", Nbins
         Close(10) 
         !Close(12)
         
         NP = NBINS
         NOBS = 10
         
         ALLOCATE(OBS(NP,NOBS))
         !	Error on energy

         !Open (Unit=25, File="statdat1", status="unknown") 
         !read(25,*) NST, NS1, NS2, NSTEP
         !Close(25)
         NST = 0; NS1 = 5; NS2 = 30; NSTEP = 5
         !If ( L == 15 ) NST = 10
         !If ( L == 12 ) NST = 8 
         !If ( L == 9  ) NST = 3 
         !If ( L == 6  ) NST = 2 
         !If ( L == 3  ) NST = 2 
         OPEN (UNIT=20, FILE='ener', STATUS='old')
         Open ( Unit=55, File="U1_timeline", status="unknown" )
         x1=0.d0
         y1=0.d0
         NC = 0
         DO N = 1,NP
            IF (N.GE.NST) THEN
               NC = NC + 1
               READ(20,*) Z1,Z2,Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10!, Z11
               OBS(NC,1) = Z1 
               OBS(NC,2) = Z2 
               OBS(NC,3) = Z3
               OBS(NC,4) = Z4 
               OBS(NC,5) = Z5
               OBS(NC,6) = Z6
               OBS(NC,7) = Z7
               OBS(NC,8) = Z8
               OBS(NC,9) = Z9
               OBS(NC,10) = Z10
!                OBS(NC,11) = Z11
               X=Obs(nc,8)
               X1=x1+X
               y=Obs(nc,9)
               y1=y1+y
               write(55,*) dble(X), aimag(X), dble(x1/(nc)), aimag(x1/(nc)), dble(y), aimag(y), dble(y1/(nc)), aimag(y1/(nc))
            ELSE
               READ(20,*) Z1,Z2,Z3, Z4, Z5, Z6, Z5, Z6, Z5, Z5!, Z11
            ENDIF
         ENDDO
         CLOSE(20)
2100     FORMAT(I6,2X,F16.8)
         
         OPEN (UNIT=21, FILE='enerJ', STATUS='unknown')
         WRITE(21,*) 'Effective number of bins, and bins: ', NC, NP
         NP_EFF = NC
         ALLOCATE (EN(NP_EFF), SIGN(NP_EFF))
         DO IOBS = 1,NOBS
            WRITE(21,*)
            DO I = 1,NP_EFF
               EN  (I) = OBS(I,IOBS)
               SIGN(I) = OBS(I,NOBS)
            ENDDO
            IF (IOBS.EQ.1) WRITE(21,*) ' rho           '
            IF (IOBS.EQ.2) WRITE(21,*) ' kin           '
            IF (IOBS.EQ.3) WRITE(21,*) ' pot           '
            IF (IOBS.EQ.4) WRITE(21,*) ' Energy        '
            IF (IOBS.EQ.5) WRITE(21,*) ' dF/dt         '
            IF (IOBS.EQ.6) WRITE(21,*) ' dF/dU         '
            IF (IOBS.EQ.7) WRITE(21,*) ' dF/dl         '
            IF (IOBS.EQ.8) WRITE(21,*) ' U1            '
            IF (IOBS.EQ.9) WRITE(21,*) ' U1xG          '
!             IF (IOBS.EQ.10) WRITE(21,*) ' U1yG          '
            IF (IOBS.EQ.NOBS) WRITE(21,*) ' phase         '
            DO NBIN = NS1, NS2, NSTEP
               if (NBIN.gt.0) then
                  IF (IOBS.EQ.NOBS) then 
                     CALL ERRCALCJ(EN,XM,XERR,NBIN)
                  else
                     CALL ERRCALCJ(EN,SIGN,XM,XERR,NBIN)
                  endif
                  WRITE(21,2001) IOBS, dble(XM),  dble(XERR), aimag(XM),  aimag(XERR)
                  ! Test
                  ! NBOOT = 40
                  ! CALL BOOTSTRAP( EN,XM_BS,XERR_BS,NBOOT,ISEED)
                  ! WRITE(21,2001) IOBS, XM_BS,  XERR_BS
                  ! IF (IOBS == 4) Write(22,"(F14.7,2x,F14.7)")  XM/dble(L*L), XERR/dble(L*L)
               endif
            ENDDO
         ENDDO
         CLOSE(21)
2001     FORMAT('OBS : ', I4,4x,F12.6,2X, F12.6,2x,F12.6,2X, F12.6)
         
         DEALLOCATE (EN,SIGN,OBS)
         
       END Program enerJ
       
