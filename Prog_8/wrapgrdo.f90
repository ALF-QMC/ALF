      SUBROUTINE WRAPGRDO(GR,NTAU,PHASE)

        Use Hamiltonian
        Use MyMats
        Use Hop_mod
        Implicit None

        Interface
           Subroutine Upgrade(GR,N_op,NT,PHASE,Op_dim) 
             Use Hamiltonian 
             Implicit none 
             Complex (Kind=8) :: GR(Ndim,Ndim, N_FL) 
             Integer, INTENT(IN) :: N_op, Nt, Op_dim
             Complex (Kind=8) :: Phase
           End Subroutine Upgrade
        End Interface
        
        ! Given GREEN at time NTAU => GREEN at time NTAU - 1,
        ! Upgrade NTAU  [LTROT:1]
        
        COMPLEX (Kind=8), INTENT(INOUT) :: GR(Ndim,Ndim,N_FL)
        COMPLEX (Kind=8), INTENT(INOUT) :: PHASE
        Integer :: NTAU

        ! Local
        Complex (Kind=8) :: Mat_TMP(Ndim,Ndim)
        Integer :: nf, N_Type, n
        real (Kind=8) :: spin

        Do n = size(Op_V,1), 1, -1
           N_type = 2
           nf = 1
           spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
           do nf = 1,N_FL
              Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), spin, Ndim, N_Type)
           enddo
           !Write(6,*) 'Upgrade : ', ntau,n 
           Call Upgrade(GR,n,ntau,PHASE,Op_V(n,1)%N_non_zero) 
           ! The spin has changed after the upgrade!
           nf = 1
           spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
           N_type = 1
           do nf = 1,N_FL
              Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), spin, Ndim, N_Type )
           enddo
        enddo
        DO nf = 1,N_FL
           Call Hop_mod_mmthl   (GR(:,:,nf), MAT_TMP, nf)
           Call Hop_mod_mmthr_m1(MAT_TMP, GR(:,:,nf), nf)
           !CALL MMULT(MAT_TMP   , GR(:,:,nf)      , Exp_T(:,:,nf) )
           !CALL MMULT(GR(:,:,nf), Exp_T_M1(:,:,nf), MAT_TMP       )
        enddo

!!$        ! Test
!!$        Mat_TMP = cmplx(0.d0,0.d0)
!!$        DO I = 1,Ndim
!!$           Mat_TMP(I,I) = cmplx(1.d0,0.d0) 
!!$        Enddo
!!$        Do n = size(Op_V,1), 1, -1
!!$           N_type = 2
!!$           nf = 1
!!$           spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
!!$           Write(6,*) n, spin
!!$           do nf = 1,N_FL
!!$              Call Op_Wrapdo( Mat_tmp, Op_V(n,nf), spin, Ndim, N_Type)
!!$           enddo
!!$           !Upgrade
!!$           N_type = 1
!!$           do nf = 1,N_FL
!!$              Call Op_Wrapdo( Mat_tmp, Op_V(n,nf), spin, Ndim, N_Type )
!!$           enddo
!!$        enddo
!!$        
!!$        DO I = 1,Ndim
!!$           Do J = 1,NDIM
!!$              WRITE(6,*) I,J, Mat_tmp(I,J)
!!$           ENDDO
!!$        ENDDO
!!$
!!$        STOP

      END SUBROUTINE WRAPGRDO
