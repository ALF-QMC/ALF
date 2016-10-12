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
        COMPLEX (Kind=8), Dimension(:,:), Allocatable ::  UUP, VUP, TPUP, TPUP1, TPUPM1, TPUP1M1, UUPM1, VUP1
        COMPLEX (Kind=8), Dimension(:) , Allocatable ::  DUP
        COMPLEX (Kind=8) ::  ZDUP1, ZDDO1, ZDUP2, ZDDO2, Z1, ZUP, ZDO
        Integer :: I,J, N_size, NCON
        Real (Kind=8) :: X, Xmax
        
        N_size = SIZE(DLUP,1)
        NCON = 0

        Allocate( UUP(N_size,N_size),  VUP(N_size,N_size), TPUP(N_size,N_size), TPUP1(N_size,N_size), &
             & TPUPM1(N_size,N_size),TPUP1M1(N_size,N_size),  UUPM1(N_size,N_size), VUP1(N_size,N_size), DUP(N_size) )
	
        !Write(6,*) 'In CGR', N_size
        CALL MMULT(VUP,VRUP,VLUP)
        DO J = 1,N_size
            TPUP(:,J) = DRUP(:)*VUP(:,J)*DLUP(J)
        ENDDO
        CALL MMULT(UUP,ULUP,URUP)
! Fusing the CT and the Matrix Addition breaks the vectorization on GCC. Hence only benchmarks can decide.
        UUPM1 = CT(UUP)
        TPUP=TPUP + UUPM1
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
        ELSE
           !WRITE(6,*) 'UDV of (U + DR * V * DL)^{*}'
           TPUP1 = CT(TPUP)
           CALL UDV_WRAP(TPUP1,UUP,DUP,VUP,NCON)
           !CALL UDV(TPUP1,UUP,DUP,VUP,NCON)
           TPUP = CT(ULUP)
           CALL MMULT(TPUPM1,TPUP,UUP)
           VUP1 = CT(VUP)
           CALL MMULT(TPUP1,URUP,VUP1)
           CALL INV(TPUP1,TPUP1M1,ZDUP2)
           CALL INV(TPUPM1, TPUP, ZDUP1)
           Z1 = ZDUP2/ZDUP1
        ENDIF
        DO I = 1,N_size
           X = ABS(DUP(I))
           if (I == 1)  Xmax = X
           if ( X  < Xmax ) Xmax = X
        ENDDO
        !Write(6,*) 'Cgr1, Cutoff: ', Xmax

        DO J = 1,N_size
        DO I = 1,N_size
           GRUP(I, J) = Sum(TPUPM1(I,:) * TPUP1M1(:,J)/DUP)
        ENDDO
        ENDDO
        PHASE = Z1/ABS(Z1)

        Deallocate(UUP, VUP, TPUP,TPUP1,TPUPM1, TPUP1M1, UUPM1, VUP1, DUP )

      END SUBROUTINE CGR
