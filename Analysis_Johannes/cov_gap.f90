      PROGRAM COV_GAP

        
        Use Errors
        Use MyMats
        Use Random_wrap
        IMPLICIT REAL (KIND=8) (A-G,O-Z)

        REAL (KIND=8), DIMENSION(:,:), ALLOCATABLE ::  HELP, COV1, UC1
        REAL (KIND=8), DIMENSION(:  ), ALLOCATABLE ::  SIG1, GRT0, XTAU, XDATA, FDATA, FDATA1, VHELP, &
             &                                         ERROR, ARES, AME1, AME, ASQ, ERRA
        REAL (KIND=8)                              ::  CHISQ, CHISQ1
        REAL (KIND=8), EXTERNAL   :: F
        Integer, allocatable :: Seed_start(:)

!        Nred=3
!	NSETS = 100
        OPEN (UNIT=55, FILE='stat_gap', STATUS='unknown')
        READ(55,*) BETAST, BETAEN
        READ(55,*) NCOV
        READ(55,*) Nred
        READ(55,*) NSETS
        CLOSE(55)
!        print*,'Nred=',Nred

        OPEN (UNIT=20,FILE='g_dat',STATUS='unknown')
        READ(20,*) NTDM
        ALLOCATE (GRT0(NTDM/Nred), XTAU(NTDM/Nred),  HELP(NTDM/Nred,NTDM/Nred) )
        Help=0.d0
        DO NT = 1,NTDM
           READ(20,*) helpx, helpy, helperr, X
           if (NT.ge.Nred) then
           XTAU(NT/Nred)=helpx
           GRT0(NT/Nred)=helpy
           help(NT/Nred,NT/Nred)=helperr*helperr
           endif
        ENDDO
        if (Ncov==1) then
        DO NT1 = 1,NTDM
           DO NT2 = 1,NTDM
              READ(20,*) HELPx
              if ((NT1.ge.Nred).and.(NT2.ge.Nred)) then
               HELP(NT1/Nred,NT2/Nred)=HELPx
              endif
           ENDDO
        ENDDO
        endif
        
        DTAU = XTAU(2) - XTAU(1)
        Open(Unit=25,File="Test_out",status="unknown")
           do nt = 1,Ntdm/Nred
              Write(25,"(F14.7,2x,F24.12,2x,F24.12)") xtau(nt), Grt0(nt), sqrt(help(nt,nt))
           enddo
        Close(25)


	      
        NTST = NINT(BETAST/DTAU) 
        NTEN = NINT(BETAEN/DTAU)
        NDATA= NTEN - NTST 
        OPEN (UNIT=25, FILE='Gap_res', STATUS='unknown')
        WRITE(25,*) '# Fit from: ', FLOAT(NTST)*DTAU, FLOAT(NTEN)*DTAU


        ALLOCATE( COV1(NDATA,NDATA), UC1(NDATA,NDATA), SIG1(NDATA) )
        ALLOCATE( XDATA(NDATA), FDATA(NDATA), FDATA1(NDATA), ERROR(NDATA), VHELP(NDATA) )
        DO NT1 = 1,NDATA
           DO NT2 = 1,NDATA
              COV1(NT1,NT2) = HELP(NT1+NTST,NT2+NTST)
           ENDDO
        ENDDO
        IF (NCOV.EQ.1) THEN
           print*, 'Cov taken into account'
           CALL DIAG(COV1,UC1,SIG1)
        ELSE
           print*, 'Cov not taken into account'
           UC1 = 0.D0
           DO NT = 1,NDATA
              SIG1(NT) = COV1(NT,NT)
              UC1(NT,NT) = 1.D0
           ENDDO
        ENDIF
!        Do NT=1,NDATA
!         print*,'NT=',NT,'SIG1=',SIG1(NT)
!        Enddo
        
        !Calculate the mean, no errors.
        DO I = 1,NDATA
           NT = I + NTST
           XDATA(I) = XTAU(NT)
           FDATA(I) = LOG(GRT0(NT))
           ERROR(I) = SQRT(COV1(I,I))/GRT0(NT)
        ENDDO
        NBASIS = 2
        ALLOCATE (ARES(NBASIS), AME1(NBASIS), AME(NBASIS), ASQ(NBASIS), ERRA(NBASIS) ) 
        CALL FIT(XDATA,FDATA,ERROR,ARES,CHSQ1,F)
        ARES(1) = EXP(ARES(1))
        AME1(1) = ARES(1)
        AME1(2) = ARES(2)
        WRITE(25,2012)  ARES(1),ARES(2)
2012    FORMAT(' # Mean of fit:', F16.8,2x,F16.8)
        !Calculate errors.
        ISEED = 73288979
        ASQ = 0.D0; AME = 0.D0
        DO I = 1,NDATA
           X = 0.D0
           DO M = 1,NDATA
              X  = X + UC1(M,I)* GRT0(M + NTST)
           ENDDO
           VHELP(I) = X
        ENDDO
  
        call Get_seed_Len(nK)   
        Allocate         (SEED_start(nK) )
        Seed_start=0
        call Ranset(seed_start)
        DO NTH = 1,NSETS
           DO I = 1,NDATA
              FDATA(I) = VHELP(I) + SQRT(SIG1(I))*rang_wrap(ISEED)
           ENDDO
           DO I = 1,NDATA
              X = 0.D0
              DO M = 1,NDATA
                 X  = X + UC1(I,M)*FDATA(M)
              ENDDO
              FDATA1(I) = X
           ENDDO
           DO I = 1,NDATA
              NT = I + NTST
              XDATA(I) = XTAU(NT)
              ERROR(I) = SQRT(COV1(I,I))/FDATA1(I)
              FDATA1(I)= LOG(FDATA1(I))
           ENDDO
           CALL FIT(XDATA,FDATA1,ERROR,ARES,CHSQ,F)
              
           ARES(1) = EXP(ARES(1))
           DO I = 1,NBASIS
!           print*,I,'ARES',ARES(I)
              ASQ(I) = ASQ(I) + ARES(I)*ARES(I)
              AME(I) = AME(I) + ARES(I)
           ENDDO
        ENDDO
   
        DO I = 1,NBASIS
           ASQ(I) = ASQ(I)/FLOAT(NSETS)
           AME(I) = AME(I)/FLOAT(NSETS)
           ERRA(I)= 0.D0
           TMP = ASQ(I) - AME(I)*AME(I)
!           print*,'TMP=',TMP
           IF (TMP.GT.0) ERRA(I) = SQRT( TMP )
        ENDDO

        WRITE(25,*) '# Chisquare: ', CHSQ1
        WRITE(25,2006)  AME1(1) ,ERRA(1)
        WRITE(25,2006) -AME1(2) ,ERRA(2)
        
!2001    FORMAT(F16.8,2x,F16.8,2x,F16.8)
!2005    FORMAT(F16.8,2x,F16.8,2x,F16.8)
2006    FORMAT(F16.8,2x,F16.8)
        
        CLOSE(25)
        
        STOP
      END PROGRAM COV_GAP
      
!**********
      REAL (KIND=8)  FUNCTION F(K,X)

        INTEGER K
        REAL (KIND=8) X

        IF ( K.EQ.1) F = 1.D0
        IF ( K.EQ.2) F = X

        RETURN
      END FUNCTION F
!***********
