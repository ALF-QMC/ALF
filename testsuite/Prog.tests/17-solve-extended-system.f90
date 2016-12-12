! compile with
! gfortran -Wall -std=f2003  -I ../../Libraries/Modules/ -I ../../Prog_8/ -L ../../Libraries/Modules/ -L ../../Prog_8/  17-solve-extended-system.f90 ../../Prog_8/cgr2_2.o ../../Prog_8/UDV_WRAP.o ../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyNag/libnag.a


Program SOLVEXTENDEDSYSTEM
Use MyMats
implicit none

        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), allocatable :: HLPo, HLPn, HLPB1
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), allocatable :: UCT, VINV, U3, V3, V3t, testn, testo
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:), allocatable :: D3, TMPVEC
        Integer :: I, J, LQ, LQ2, info
        COMPLEX (KIND=KIND(0.D0)) :: alpha, beta
        INTEGER, Dimension(:), Allocatable :: IPVT
        
        alpha = 1.D0
        beta = 0.D0
        
do LQ = 1, 6
LQ2 = 2*LQ
allocate(HLPo(LQ2, LQ2), HLPn(LQ2, LQ2), UCT(LQ, LQ), VINV(LQ, LQ), U3(LQ2, LQ2), V3(LQ2, LQ2), D3(LQ2))
allocate(HLPB1(LQ2, LQ2), TMPVEC(LQ2), V3t(LQ2, LQ2), IPVT(LQ2), testn(LQ2, LQ2), testo(LQ2, LQ2))


! generate input data
        do i=1,LQ2
            do j=1,LQ2
                V3(i, j) = 1.D0/(i+j) !Hilbert Matrix
                U3(i, j) = CMPLX(i, j, kind(0.D0))
            enddo
            D3(i) = i -LQ + 0.5
        enddo
        do i=1,LQ
            do j=1,LQ
                VINV(i, j) = (i+j)
                UCT(i,j) = CMPLX(i, j, kind(0.D0))
            enddo
        enddo        

! insert example code for comparison here



! TMPVEC = conjg(1.D0/D3)
! HLPB1 = cmplx(0.d0,0.d0,kind(0.D0))
!            DO I = 1,LQ
!               DO J = 1,LQ
!                  HLPB1(I   , J    ) =  UCT(I, J)
!                  HLPB1(I+LQ, J+LQ ) =  VINV(I,J)
!               ENDDO
!            ENDDO
!            V3t = V3
!            HLPn = HLPB1
! CALL ZGETRF(LQ2, LQ2, V3t, LQ2, IPVT, info)
!            CALL ZGETRS('C', LQ2, LQ2, V3t, LQ2, IPVT, HLPB1, LQ2, info)! Block structure of HLPB1 is not exploited
!            
!            
!             CALL ZGEMM('C', 'N', LQ2, LQ2, LQ2, alpha, V3, LQ2, HLPB1, LQ2, beta, testn, LQ2) ! Block structure of HLPB1 is not exploited
!             testn = testn - HLPn
!            
!            
!            DO J = 1,LQ2
!               DO I = 1,LQ2
!                  HLPB1(I,J)  = TMPVEC(I)*HLPB1(I,J)
!               ENDDO
!            ENDDO
!            CALL MMULT(HLPn,U3,HLPB1)
        V3t = V3
        call solve_extended_system(HLPn, UCT, VINV, U3, D3, V3t, LQ)
        
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
           
           !create test matrix for judging the quality of the inversion
            CALL ZGEMM('C', 'N', LQ2, LQ2, LQ2, alpha, V3, LQ2, HLPB1, LQ2, beta, testo, LQ2)
            testo = testo - HLPo
           
           
           DO J = 1,LQ2
              DO I = 1,LQ2
                 HLPB1(I,J)  = TMPVEC(I)*HLPB1(I,J)
              ENDDO
           ENDDO
           CALL MMULT(HLPo,U3,HLPB1)

    do i = 1,LQ2
        do j = 1,LQ2
        if (Abs(HLPo(i, j) - HLPn(i, j)) > MAX(ABS(HLPn(i, j)), ABS(HLPo(i, j))) * 1E-15) then
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
deallocate(HLPo, HLPn, UCT, VINV, U3, V3, D3, HLPB1, TMPVEC, V3t, IPVT, testo, testn)
enddo
end Program SOLVEXTENDEDSYSTEM
