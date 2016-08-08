      SUBROUTINE WRAPGRUP(GR,NTAU,PHASE)

        Use Hamiltonian
        Use Hop_mod
        Implicit none

        Interface
           Subroutine Upgrade(GR,N_op,NT,PHASE,Op_dim) 
             Use Hamiltonian 
             Implicit none 
             Complex (Kind=8) :: GR(Ndim,Ndim, N_FL) 
             Integer, INTENT(IN) :: N_op, Nt, Op_dim
             Complex (Kind=8) :: Phase
           End Subroutine Upgrade
        End Interface

        ! Given GRUP at time NTAU => GRUP at time NTAU + 1.
	! Upgrade NTAU + 1     NTAU: [0:LTROT-1]

        ! Arguments
	COMPLEX (Kind=8), INTENT(INOUT) ::  GR(Ndim,Ndim,N_FL)
        COMPLEX (Kind=8), INTENT(INOUT) ::  PHASE
        INTEGER, INTENT(IN) :: NTAU
        
        !Local 
        Integer :: nf, N_Type, NTAU1,n 
        Complex (Kind=8) :: Mat_TMP(Ndim,Ndim)
        Real    (Kind=8) :: X

        ! Wrap up, upgrade ntau1.  with B^{1}(tau1) 
        NTAU1 = NTAU + 1
        Do nf = 1,N_FL
           CALL HOP_MOD_mmthr( GR(:,:,nf),  MAT_TMP,nf)
           CALL HOP_MOD_mmthl_m1(MAT_TMP,GR(:,:,nf), nf )
           !CALL MMULT ( MAT_TMP,    Exp_T(:,:,nf), GR(:,:,nf)        )
           !CALL MMULT ( GR(:,:,nf), MAT_TMP      , Exp_T_M1(:,:,nf)  )
        Enddo
        Do n = 1,Size(Op_V,1)
           Do nf = 1, N_FL
              X = Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
              N_type = 1
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
           enddo
           nf = 1
           !Write(6,*) 'Upgrade: ', ntau1,n
           Call Upgrade(GR,N,ntau1,PHASE,Op_V(n,nf)%N) 
           do nf = 1,N_FL
              N_type =  2
              Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
           enddo
        Enddo
        
      END SUBROUTINE WRAPGRUP
