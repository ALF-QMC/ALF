      SUBROUTINE CGR(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)

        Use UDV_Wrap_mod

        Implicit None

	!!!  GRUP = (1 + UR*DR*VR*VL*DL*UL)^-1
	!!!  NVAR = 1 Big scales are in DL
	!!!  NVAR = 2 Big scales are in DR

	!Arguments.
        COMPLEX(Kind=8), Dimension(:,:), Intent(IN)   ::  URUP, VRUP, ULUP, VLUP
        COMPLEX(Kind=8), Dimension(:),   Intent(In)   ::  DLUP, DRUP
        COMPLEX(Kind=8), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX(Kind=8) :: PHASE
        INTEGER         :: NVAR
 
	!Local
	COMPLEX (Kind=8), Dimension(:,:), Allocatable ::  UUP,  VUP, TPUP, TPUP1, &
             &	                                          TPUPM1,TPUP1M1,  UUPM1, VUP1
        COMPLEX (Kind=8), Dimension(:) , Allocatable ::  DUP
	COMPLEX (Kind=8) ::  ZDUP1, ZDDO1, ZDUP2, ZDDO2, Z1, ZUP, ZDO, Z
        Integer :: I,J, N_size, NCON, NR, NT, N
        Real (Kind=8) :: X, Xmax
        
        N_size = SIZE(DLUP,1)
	NCON = 0

        Allocate( UUP(N_size,N_size),  VUP(N_size,N_size), TPUP(N_size,N_size), TPUP1(N_size,N_size), &
             & 	  TPUPM1(N_size,N_size),TPUP1M1(N_size,N_size),  UUPM1(N_size,N_size), VUP1(N_size,N_size), DUP(N_size) )
	
        !Write(6,*) 'In CGR', N_size
        CALL MMULT(VUP,VRUP,VLUP)
        DO J = 1,N_size
        DO I = 1,N_size
           TPUP(I,J) = DRUP(I)*VUP(I,J)*DLUP(J)
	ENDDO
	ENDDO
        CALL MMULT(UUP,ULUP,URUP)
        DO J = 1,N_size
        DO I = 1,N_size
           UUPM1(I,J) =  CONJG(UUP(J,I))
	ENDDO
	ENDDO
        DO J = 1,N_size
        DO I = 1,N_size
           TPUP(I,J) = TPUP(I,J) + UUPM1(I,J)
	ENDDO
	ENDDO
        IF (NVAR.EQ.1) THEN
           !WRITE(6,*) 'UDV of U + DR * V * DL'
           CALL UDV_WRAP(TPUP,UUP,DUP,VUP,NCON)
           !CALL UDV(TPUP,UUP,DUP,VUP,NCON)
           CALL MMULT(TPUP,VUP,ULUP)
           !Do I = 1,N_size
           !   Write(6,*) DLUP(I)
           !enddo
           CALL INV(TPUP,TPUPM1,ZDUP1)
           !WRITE(6,*) 'End called Inv'
           CALL MMULT(TPUP1,URUP,UUP)
           CALL INV(TPUP1,TPUP1M1,ZDUP2)
           Z1 = ZDUP1*ZDUP2
        ELSEIF (NVAR.EQ.2) THEN
           !WRITE(6,*) 'UDV of (U + DR * V * DL)^{*}'
           DO J = 1,N_size
           DO I = 1,N_size
              TPUP1(I,J) = CONJG( TPUP(J,I) )
	   ENDDO
	   ENDDO
           CALL UDV_WRAP(TPUP1,UUP,DUP,VUP,NCON)
           !CALL UDV(TPUP1,UUP,DUP,VUP,NCON)
           DO J = 1,N_size
           DO I = 1,N_size
              TPUP(I,J) = CONJG( ULUP(J,I) )
	   ENDDO
	   ENDDO
           CALL MMULT(TPUPM1,TPUP,UUP)
           DO J = 1,N_size
           DO I = 1,N_size
              VUP1(I,J) = CONJG( VUP(J,I) )
	   ENDDO
	   ENDDO
           CALL MMULT(TPUP1,URUP,VUP1)
           CALL INV(TPUP1,TPUP1M1,ZDUP2)
           CALL INV(TPUPM1, TPUP, ZDUP1)
           Z1 = ZDUP2/ZDUP1
        ENDIF
        DO I = 1,N_size
           Z =  DUP(I)
           if (I == 1)  Xmax = real(SQRT( Z* conjg(Z)),kind=8) 
           if ( real(SQRT( Z* conjg(Z)),kind=8)  < Xmax ) Xmax = &
                & real(SQRT( Z* conjg(Z)),kind=8)
        ENDDO
        !Write(6,*) 'Cgr1, Cutoff: ', Xmax


        DO J = 1,N_size
        DO I = 1,N_size
           ZUP = CMPLX(0.D0,0.D0)
           DO NR = 1,N_size
             ZUP = ZUP + TPUPM1(I,NR)*TPUP1M1(NR,J)/DUP(NR)
	   ENDDO
           GRUP(I,J) = ZUP
	ENDDO
	ENDDO
        PHASE = Z1/SQRT( Z1* CONJG(Z1) )

        Deallocate(UUP, VUP, TPUP,TPUP1,TPUPM1, TPUP1M1, UUPM1, VUP1, DUP )

      END SUBROUTINE CGR
