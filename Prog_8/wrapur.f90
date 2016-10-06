      SUBROUTINE WRAPUR(NTAU, NTAU1, UR, DR, VR)

        ! Given    B(NTAU,  1 ) =  UR, DR, VR
        ! Returns  B(NTAU1, 1 ) =  UR, DR, VR
        ! NOTE:    NTAU1 > NTAU.
        
        Use Hamiltonian
        Use UDV_Wrap_mod
        Use Hop_mod
        Implicit None

        ! Arguments
        COMPLEX (KIND=8) :: UR(Ndim,Ndim,N_FL), VR(Ndim,Ndim,N_FL)
        COMPLEX (KIND=8) :: DR(Ndim,N_FL)
        Integer :: NTAU1, NTAU


        ! Working space.
        Complex (Kind=8) :: Z_ONE
        COMPLEX (Kind=8) :: V1(Ndim,Ndim), TMP(Ndim,Ndim), TMP1(Ndim,Ndim)
        Integer ::NCON, NT, I, J, n, nf
        Real (Kind=8) :: X

        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        
        Do nf = 1,N_FL
           CALL INITD(TMP,Z_ONE)
           DO NT = NTAU + 1, NTAU1
              !CALL MMULT(TMP1,Exp_T(:,:,nf) ,TMP)
              Call Hop_mod_mmthr(TMP,TMP1,nf)
              TMP = TMP1
              Do n = 1,Size(Op_V,1)
                 X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(Tmp,Op_V(n,nf),X,Ndim)
              ENDDO
           ENDDO
           CALL MMULT(TMP1,TMP,UR(:,:,nf))
           DO J = 1,NDim
              DO I = 1,NDim
                 TMP1(I,J) = TMP1(I,J)*DR(J,nf)
                 TMP(I,J)  = VR(I,J,nf)
              ENDDO
           ENDDO
           CALL UDV_WRAP(TMP1,UR(:,:,nf),DR(:,nf),V1,NCON)
           CALL MMULT(VR(:,:,nf),V1,TMP)
        ENDDO
      END SUBROUTINE WRAPUR
