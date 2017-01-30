! compile with
! gfortran -Wall -std=f2003 -I ../../../Prog_8/  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/cgr1.o ../../../Prog_8/UDV_WRAP.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program TESTCGR
implicit none
interface
      SUBROUTINE CGRold(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)
        COMPLEX(Kind=8), Dimension(:,:), Intent(IN)   ::  URUP, VRUP, ULUP, VLUP
        COMPLEX(Kind=8), Dimension(:),   Intent(In)   ::  DLUP, DRUP
        COMPLEX(Kind=8), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX(Kind=8) :: PHASE
        INTEGER, Intent(In)         :: NVAR
        end subroutine CGRold
        
      SUBROUTINE CGR(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)
        COMPLEX(Kind=8), Dimension(:,:), Intent(IN)   ::  URUP, VRUP, ULUP, VLUP
        COMPLEX(Kind=8), Dimension(:),   Intent(In)   ::  DLUP, DRUP
        COMPLEX(Kind=8), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX(Kind=8) :: PHASE
        INTEGER, Intent(In)         :: NVAR
        end subroutine CGR
end interface

        COMPLEX(Kind=Kind(0.D0)), Dimension(:,:), allocatable ::  URUP, VRUP, ULUP, VLUP
        COMPLEX(Kind=Kind(0.D0)), Dimension(:), allocatable ::  DLUP, DRUP
        COMPLEX(Kind=Kind(0.D0)), Dimension(:,:), allocatable :: GRUPnew, GRUPold
        COMPLEX(Kind=Kind(0.D0)) :: PHASEnew, Phaseold, Zre, Zim
        INTEGER         :: NVAR, i, j, N_size
        
        N_size = 5
        
        do NVAR = 1,2
        allocate(URUP(N_size, N_size), VRUP(N_size, N_size), ULUP(N_size, N_size), VLUP(N_size, N_size))
        allocate(GRUPnew(N_size, N_size), GRUPold(N_size, N_size), DLUP(N_size), DRUP(N_size))
        ! set up test data
        Phasenew = 1.D0
        Phaseold = 1.D0
        GRUPnew = 0.D0
        GRUPold = 0.D0
        do i = 1, N_size
        do j = 1, N_size
        URUP(i, j) = 0!i + j
        ULUP(i, j) = 0!i + j
        VRUP(i, j) = CMPLX(i,j, kind=kind(0.D0))
        VLUP(i, j) = i + j
        enddo
        DLUP(i) = 1
        ULUP(i,i) = 1
        URUP(i,i) = 1
        VRUP(i, i) = CMPLX(1,1)
        VLUP(i, i) = CMPLX(1,1)
        DRUP(i) = i
        enddo
call CGR(PHASEnew, NVAR, GRUPnew, URUP, DRUP, VRUP, ULUP, DLUP, VLUP)
! run old code

call CGRold(PHASEold, NVAR, GRUPold, URUP, DRUP, VRUP, ULUP, DLUP, VLUP)

! compare GRUP results
! write (*,*) GRUPOLD
! write (*,*) "---------"
! write (*,*) GRUPNEW
    do i = 1, N_size
    do j = 1, N_size
    Zre = DBLE(GRUPnew(i,j)-GRUPold(i,j))
    Zim = aimag(GRUPnew(i,j)-GRUPold(i,j))
    if (Abs(Zre) > MAX(ABS(real(GRUPnew(i,j))), ABS(real(GRUPold(i,j))) )*1D-10) then
    write (*,*) "ERROR in real part", i, j, real(GRUPnew(i,j)), real(GRUPold(i,j)), "diff: ", &
    & Abs(Zre), "prec: ", MAX(ABS(real(GRUPnew(i,j))), ABS(real(GRUPold(i,j))) )*1D-7
    STOP 2
    endif
    if (Abs(Zim) > MAX(ABS(aimag(GRUPnew(i,j))), ABS(aimag(GRUPold(i,j))) )*1D-10) then
    write (*,*) "ERROR in imag part", aimag(GRUPnew(i,j)), aimag(GRUPold(i,j))
    STOP 3
    endif
    enddo
    enddo

    Zre = real(Phasenew-Phaseold)
    Zim = aimag(Phasenew-Phaseold)
    if (Abs(Zre) > MAX(ABS(real(Phasenew)), ABS(real(Phaseold)) )*1D-13) then
    write (*,*) "ERROR in real part", Phasenew, Phaseold, Abs(Phaseold)
    STOP 6
    endif
    if (Abs(Zim) > MAX(ABS(aimag(Phasenew)), ABS(aimag(Phaseold)) )*1D-13) then
    write (*,*) "ERROR in imag part", Phasenew, aimag(Phasenew), aimag(Phaseold)
    STOP 7
    endif
deallocate(URUP, VRUP, ULUP, VLUP, GRUPnew, GRUPold, DLUP, DRUP)
enddo
write (*,*) "success"
end Program TESTCGR

      SUBROUTINE CGRold(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)

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
        INTEGER, Intent(In)         :: NVAR
 
        !Local
        COMPLEX (Kind=8), Dimension(:,:), Allocatable ::  UUP, VUP, TPUP, TPUP1, TPUPM1, TPUP1M1, UUPM1, VUP1
        COMPLEX (Kind=8), Dimension(:) , Allocatable ::  DUP
        COMPLEX (Kind=8) ::  ZDUP1, ZDUP2, Z1
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
        UUPM1 = Conjg(Transpose((UUP)))
        TPUP=TPUP + UUPM1
        IF (NVAR.EQ.1) THEN
           !WRITE(6,*) 'UDV of U + DR * V * DL'
           CALL UDV_WRAP_Pivot(TPUP,UUP,DUP,VUP,NCON, N_size, N_size)
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
           TPUP1 = Conjg(transpose((TPUP)))
           CALL UDV_WRAP_Pivot(TPUP1,UUP,DUP,VUP,NCON, N_size, N_size)
           !CALL UDV(TPUP1,UUP,DUP,VUP,NCON)
           TPUP = Conjg(Transpose((ULUP)))
           CALL MMULT(TPUPM1,TPUP,UUP)
           VUP1 = Conjg(Transpose((VUP)))
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

      END SUBROUTINE CGRold
