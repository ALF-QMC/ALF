      Subroutine Gperp_sub( G, Gperp, Ndim,Irank) 
        
        Use Precdef
        Use MyMats
        Implicit none
        
        ! Arguments
        Integer, Intent(In) :: Ndim, Irank
        Complex (kind=double), Intent(In)    ::  G(ndim,ndim)
        Complex (kind=double), Intent(InOut) ::  Gperp(ndim,ndim)

        ! Local space
        Complex   (Kind=double) :: A(ndim,ndim), W(ndim), VL(Ndim,ndim), VR(Ndim,ndim)
        Character (len=1)       :: JOBVL, JOBVR
        Integer                 :: INFO, LDA, LDVL, LDVR, N, lp, LWORK, N_c,m, i, j, NCon
        Complex (Kind=double)   :: WORK(2*Ndim), U(Ndim,Ndim/2), Vec(Ndim),Z 
        Real    (Kind=double)   :: RWORK(2*ndim), X, Xmax, Xmean
        Complex (Kind=double)   :: U1(Ndim,Ndim/2),  V(Ndim/2,Ndim/2),  D(Ndim/2)
        

        A = G
        JOBVL = "N"
        JOBVR = "V"
        LDA   =  Ndim
        LWORK =  2*Ndim 
        LDVL  = Ndim
        LDVR  = Ndim
        N     = Ndim

        Call  ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
        
        !lp = 70 + Irank
        !Write(lp,*) "Info: ", INFO
        N_c = 0
        do n = 1,Ndim
           !Write(lp,*) n,  W(n)
           if ( abs( dble(W(n)) ) <  0.00001 ) then
              N_c = N_c + 1
              U1(:, N_c) = VR(:, n)
           endif
        enddo
        !Write(6,*) "N_c ", N_c 
        NCON = 0
        Call UDV  (U1,U,D,V,NCON)
        
        ! Setup G_perp
        Do i = 1,Ndim
           do j = 1,Ndim
            Gperp(i, j) = dot_product(U(j, :), U(i, :))
           enddo
        enddo

#ifdef Test_gperp
        X = 0.05
        A  =  cmplx(1.d0-X,0.d0, kind(0.D0)) * G    + cmplx(x,0.d0, kind(0.D0))*Gperp
        Call Inv(A,VR,Z) 
        Write(lp,*) "Det is ", Z 
        Call MMult(VL,A,VR)
        VR = cmplx(0.d0,0.d0, kind (0.D0))
        do i = 1,Ndim
           VR(I,I)  = cmplx(1.d0,0.d0, kind(0.D0))
        enddo
        Call Compare(VL,VR,Xmax,Xmean) 
        Write(lp,*) 'Compare: ', Xmax, Xmean

        ! This is for testing
        do n = 1,N_c
           Vec  = cmplx(0.d0,0.d0, kind(0.D0))
           do i = 1,Ndim
              do j = 1,Ndim
                 Vec(i)  = Vec(i) + G(i,j) * U(j,n)
              enddo
           enddo
           X = 0.d0
           do i  = 1,Ndim
              X = X  + dble( Vec(i) * conjg(Vec(i)))
           enddo
           X = sqrt(x)
           Write(lp,*) 'n, G*v = ', n, X
        enddo

        do n = 1,N_c
           do m = n,N_c
              Z = cmplx(0.d0,0.d0, kind(0.D0))
              do j = 1,Ndim
                 Z = Z + Conjg(U(j,m)) * U(j,n)
              enddo
              Write(lp,*) "n,m,z ", n,m,z
           enddo
        enddo
#endif

      end Subroutine Gperp_sub
