      SUBROUTINE WRAPUL(NTAU1, NTAU, UL ,DL, VL)

        !Given    B(LTROT,NTAU1,Nf  ) =  VLUP, DLUP, ULUP
        !Returns  B(LTROT,NTAU, Nf  ) =  VLUP, DLUP, ULUP


        !NOTE:    NTAU1 > NTAU.
        ! Does this for all replicas
        Use Hamiltonian
        Use Hop_mod
        Use UDV_Wrap_mod

        Implicit none

        ! Arguments
        COMPLEX (KIND=8) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
        COMPLEX (KIND=8) :: DL(Ndim,N_FL)
        Integer :: NTAU1, NTAU


        ! Working space.
        COMPLEX (Kind=8) ::  U1(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim), TMP1(Ndim,Ndim)
        COMPLEX (Kind=8) ::  D1(Ndim), Z_ONE, beta
        Integer :: NT, NCON, n, nf
        Real    (Kind=8) ::  X
 


        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        Do nf = 1, N_FL
           CALL INITD(TMP,Z_ONE)
           DO NT = NTAU1, NTAU+1 , -1
              Do n = Size(Op_V,1),1,-1
                 X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultL(Tmp,Op_V(n,nf),X,Ndim)
              enddo
              !CALL MMULT( TMP1,Tmp,Exp_T(:,:,nf) )
              Call  Hop_mod_mmthl (Tmp, Tmp1,nf)
              Tmp = Tmp1
           ENDDO
           
           !Carry out U,D,V decomposition.
!           TMP1 = CONJG(TRANSPOSE(TMP))
           CALL ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_ONE, TMP, Ndim, UL(:, :, nf), Ndim, beta, TMP1, Ndim)
           DO n = 1,NDim
              TMP1(:, n) = TMP1(:, n) * DL(n, nf)
           ENDDO
           CALL UDV_WRAP(TMP1,U1,D1,V1,NCON)
           !CALL UDV(TMP,U1,D1,V1,NCON)
           UL(:, :, nf) = CONJG(TRANSPOSE(U1))
           CALL ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_ONE, VL(:, :, nf), Ndim, V1, Ndim, beta, TMP1, Ndim)
           VL(:, :, nf) = TMP1
           DL(:, nf) = D1
        ENDDO
        
      END SUBROUTINE WRAPUL
      
