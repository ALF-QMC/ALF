! compile with
! gfortran -Wall -std=f2003  -I ../../Libraries/Modules/ -I ../../Prog_8/ 
! -L ../../Libraries/Modules/ -L ../../Prog_8/  17-solve-extended-system.f90 
! ../../Prog_8/cgr2_2.o ../../Prog_8/UDV_WRAP.o ../../Libraries/Modules/modules_90.a -llapack -lblas


Program SOLVEXTENDEDSYSTEM
Use MyMats
Use QDRP_mod
implicit none

        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), allocatable :: HLPo, HLPn, HLPB1
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), allocatable :: UCT, VINV, U3, V3, V3t, testn, testo, INPUT
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), allocatable :: D3, TMPVEC, TAU, WORK
        Integer :: I, J, LQ, LQ2, info, LWORK
        COMPLEX (KIND=KIND(0.D0)) :: alpha, beta
        INTEGER, Dimension(:), Allocatable :: IPVT
        LOGICAL :: FORWRD
        
        alpha = 1.D0
        beta = 0.D0
        
do LQ = 1, 3
LQ2 = 2*LQ
allocate(HLPo(LQ2, LQ2), HLPn(LQ2, LQ2), UCT(LQ, LQ), VINV(LQ, LQ), U3(LQ2, LQ2), V3(LQ2, LQ2), D3(LQ2), TAU(LQ2))
allocate(HLPB1(LQ2, LQ2), TMPVEC(LQ2), V3t(LQ2, LQ2), IPVT(LQ2), testn(LQ2, LQ2), testo(LQ2, LQ2), input(LQ2,LQ2))


! generate input data
        do i=1,LQ2
            do j=1,LQ2
                input(i,j) = 1.D0/(i+j) !Hilbert Matrix
            enddo
        enddo
        do i=1,LQ
            do j=1,LQ
                VINV(i, j) = (i+j)
                UCT(i,j) = CMPLX(1,0, kind(0.D0))
            enddo
        enddo        
IPVT = 0
call QDRP_decompose(LQ2, input, D3, IPVT, TAU, WORK, LWORK)
U3 = input
CALL ZUNGQR(LQ2, LQ2, LQ2, U3, LQ2, TAU, WORK, LWORK, INFO)
call ZLACPY('U', LQ2, LQ2, input, LQ2, V3, LQ2)
FORWRD = .false.
CALL ZLAPMT(FORWRD, LQ2, LQ2, V3, LQ2, IPVT)
call solve_extended_System(HLPn, UCT, VINV, input, D3, TAU, IPVT, LQ, WORK, LWORK)

! check old version
! by using the testo matrices that can be generated after the matrix inversion it was found 
! that the variant using ZGETRS is about 2-3 orders of magnitude more precise than Fakhers
! original implementaion(and faster, too!)
! Hence the baseline is given by this more accurate implementation.

TMPVEC = conjg(1.D0/D3)
HLPB1 = cmplx(0.d0,0.d0,kind(0.D0))
           DO I = 1,LQ
              DO J = 1,LQ
                 HLPB1(I   , J    ) =  UCT(I, J)
                 HLPB1(I+LQ, J+LQ ) =  VINV(I,J)
              ENDDO
           ENDDO
           V3t = V3
           HLPo = HLPB1
           CALL ZGETRF(LQ2, LQ2, V3t, LQ2, IPVT, info)
           CALL ZGETRS('C', LQ2, LQ2, V3t, LQ2, IPVT, HLPB1, LQ2, info)
           
           DO J = 1,LQ2
              DO I = 1,LQ2
                 HLPB1(I,J)  = TMPVEC(I)*HLPB1(I,J)
              ENDDO
           ENDDO
           CALL MMULT(HLPo,U3,HLPB1)

    do i = 1,LQ2
        do j = 1,LQ2
        if (Abs(HLPo(i, j) - HLPn(i, j)) > MAX(ABS(HLPn(i, j)), ABS(HLPo(i, j))) * 1E-01) then
        write (*,*) "LQ = ", LQ, " (i,j) = ", i,j
        write (*,*) HLPo(i, j), HLPn(i, j), HLPo(i, j) - HLPn(i, j)
!         write (*,*) "===================== old(Fakher)"
!         write (*,*) testo
!         write (*,*) "===================== new(LU)"
!         write (*,*) testn
!        STOP 1
        endif
        enddo
    enddo
deallocate(HLPo, HLPn, UCT, VINV, U3, V3, D3, HLPB1, TMPVEC, V3t, IPVT, testo, testn, WORK, TAU, input)
enddo
end Program SOLVEXTENDEDSYSTEM
