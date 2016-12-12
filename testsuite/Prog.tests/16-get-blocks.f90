! compile with
!gfortran -Wall -std=f2003  -I ../../../Libraries/Modules/ -I ../../../Prog_8/ -L ../../../Libraries/Modules/ -L ../../../Prog_8/  main.f90 ../../../Prog_8/cgr2.o ../../../Prog_8/UDV_WRAP.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program TESTGETMATRIXSUBBLOCKS

Use MyMats


        COMPLEX (KIND=KIND(0.D0)), DIMENSION(10,10) :: V
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(5,5) :: GR00n, GRTTn, GRT0n, GR0Tn
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(5,5) :: GR00o, GRTTo, GRT0o, GR0To
        Integer :: i, j, LQ, I1, J1
        
        LQ=5
        do i=1,2*LQ
            do j=1,2*LQ
                V(i, j) = cmplx(j, i, kind(0.D0))
            enddo
        enddo
        
        call get_blocks(GR00n, GR0Tn, GRT0n, GRTTn, V, LQ)
        
! check old version
        DO I = 1,LQ
          I1 = I+LQ
          DO J = 1,LQ
              J1 = J + LQ
              GR00o(I,J) = V(I ,J )
              GRTTo(I,J) = V(I1,J1)
              GRT0o(I,J) = V(I1,J )
              GR0To(I,J) = V(I,J1 )
          ENDDO
        ENDDO

    do i = 1,LQ
        do j = 1,LQ
        if (GR00o(i,j) .ne. GR00n(i,j)) then
        STOP 1
        endif
        
        if (GR0To(i,j) .ne. GR0Tn(i,j)) then
        STOP 2
        endif
        
        if (GRT0o(i,j) .ne. GRT0n(i,j)) then
        STOP 3
        endif

        if (GRTTo(i,j) .ne. GRTTn(i,j)) then
        STOP 4
        endif
        enddo
    enddo

end Program TESTGETMATRIXSUBBLOCKS
