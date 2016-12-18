      Subroutine Upgrade(GR,N_op,NT,PHASE,Op_dim) 

        Use Hamiltonian
        Use Random_wrap
        Use Control
        Implicit none 
        
        Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL) 
        Integer, INTENT(IN) :: N_op, Nt, Op_dim
        Complex (Kind=Kind(0.d0)) :: Phase

        ! Local ::
        Complex (Kind=Kind(0.d0)) :: Mat(Op_dim,Op_Dim), Delta(Op_dim,N_FL)
        Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Ratiotot, Z1 
        Integer :: ns_new, ns_old, n,m,nf, i,j
        Complex (Kind=Kind(0.d0)) :: ZK, Z, D_Mat
        Integer, external :: nranf
        
        Real    (Kind =Kind(0.d0)) :: Weight
        Complex (Kind =Kind(0.d0)) :: u(Ndim,Op_dim), v(Ndim,Op_dim) 
        Complex (Kind =Kind(0.d0)) :: x_v(Ndim,Op_dim), y_v(Ndim,Op_dim),   xp_v(Ndim,Op_dim)
        Complex (Kind =Kind(0.d0)) :: s_xv(Op_dim), s_yu(Op_dim)
        
        Logical :: Log

        
        if ( sqrt(dble(OP_V(n_op,1)%g*conjg(OP_V(n_op,1)%g))) < 1.D-6 ) return
        
        ! Compute the ratio
        nf = 1
        ns_old = nsigma(n_op,nt)
        If ( Op_V(n_op,nf)%type == 1) then
           ns_new = -ns_old
        else
           ns_new = NFLIPL(Ns_old,nranf(3))
        endif
        Do nf = 1,N_FL
           Z1 = Op_V(n_op,nf)%g * cmplx( Phi(ns_new,Op_V(n_op,nf)%type) -  Phi(ns_old,Op_V(n_op,nf)%type), 0.d0)  
           Do m = 1,Op_V(n_op,nf)%N_non_zero
              Z =  exp( Z1* Op_V(n_op,nf)%E(m) ) - cmplx(1.d0,0.d0)
              Delta(m,nf) = Z
              do n = 1,Op_V(n_op,nf)%N_non_zero
                 ZK = cmplx(0.d0,0.d0)
                 If (n == m ) ZK = cmplx(1.d0,0.d0)
                 Mat(n , m )  = ZK  +  ( ZK - GR( Op_V(n_op,nf)%P(n), Op_V(n_op,nf)%P(m),nf )) * Z
              Enddo
           Enddo
           If (Size(Mat,1) == 1 ) then
              D_mat = Mat(1,1)
           elseif (Size(Mat,1) == 2 ) then
              D_mat = Mat(1,1)*Mat(2,2) - Mat(2,1)*Mat(1,2)
           else
              D_mat = Det(Mat,Size(Mat,1))
           endif
           Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
        Enddo
        
        Ratiotot = cmplx(1.d0,0.d0)
        Do nf = 1,N_FL
           Ratiotot = Ratiotot * Ratio(nf) 
        enddo
        nf = 1
        Ratiotot = (Ratiotot**dble(N_SUN)) * cmplx(Gaml(ns_new, Op_V(n_op,nf)%type)/Gaml(ns_old, Op_V(n_op,nf)%type),0.d0)
        Ratiotot = Ratiotot*cmplx(S0(n_op,nt),0.d0)
        

        !Write(6,*) Ratiotot
        
        Weight = abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )
      
        Log = .false. 
        if ( Weight > ranf() )  Then
           Log = .true.
           Phase = Phase * Ratiotot/cmplx(weight,0.d0)
           !Write(6,*) 'Accepted : ', Ratiotot

           Do nf = 1,N_FL
              ! Setup u(i,n), v(n,i) 
              u = cmplx(0.d0,0.d0)
              v = cmplx(0.d0,0.d0)
              do n = 1,Op_V(n_op,nf)%N_non_zero
                 u( Op_V(n_op,nf)%P(n), n) = Delta(n,nf)
                 do i = 1,Ndim
                    v(i,n) = - GR( Op_V(n_op,nf)%P(n), i, nf )
                 enddo
                 v(Op_V(n_op,nf)%P(n), n)  = cmplx(1.d0,0.d0) - GR( Op_V(n_op,nf)%P(n),  Op_V(n_op,nf)%P(n), nf)
              enddo

              
              x_v = cmplx(0.d0,0.d0)
              y_v = cmplx(0.d0,0.d0)
              i = Op_V(n_op,nf)%P(1)
              x_v(i,1) = u(i,1)/(cmplx(1.d0,0.d0) + v(i,1)*u(i,1) )
              Do i = 1,Ndim
                 y_v(i,1) = v(i,1)
              enddo
              do n = 2,Op_V(n_op,nf)%N_non_zero
                 s_yu = cmplx(0.d0,0.d0)
                 s_xv = cmplx(0.d0,0.d0)
                 do m = 1,n-1
                    do i = 1,Ndim
                       s_yu(m) = s_yu(m) + y_v(i,m)*u(i,n)
                       s_xv(m) = s_xv(m) + x_v(i,m)*v(i,n)
                    enddo
                 enddo
                 Do i = 1,Ndim
                    x_v(i,n) = u(i,n)
                    y_v(i,n) = v(i,n)
                 enddo
                 Z = cmplx(1.d0,0.d0) +  u( Op_V(n_op,nf)%P(n), n)*v(Op_V(n_op,nf)%P(n),n)
                 do m = 1,n-1
                    Z = Z - s_xv(m)*s_yu(m)
                    Do i = 1,Ndim
                       x_v(i,n) = x_v(i,n) - x_v(i,m)*s_yu(m)
                       y_v(i,n) = y_v(i,n) - y_v(i,m)*s_xv(m)
                    enddo
                 enddo
                 Do i = 1,Ndim
                    x_v(i,n) = x_v(i,n)/Z
                 Enddo
              enddo
              xp_v = cmplx(0.d0,0.d0)
              do n = 1,Op_dim
                 do m = 1,Op_dim
                    j = Op_V(n_op,nf)%P(m)
                    do i = 1,Ndim
                       xp_v(i,n) = xp_v(i,n) + gr(i,j,nf)*x_v(j,n)
                    enddo
                 enddo
              enddo
              
              do n = 1,Op_dim
                 do j = 1,Ndim
                    do i = 1,Ndim
                       gr(i,j,nf) = gr(i,j,nf) - xp_v(i,n)*y_v(j,n)
                    enddo
                 enddo
              enddo
           enddo
           
           ! Flip the spin
           nsigma(n_op,nt) = ns_new
        endif

        Call Control_upgrade(Log)
        

      End Subroutine Upgrade
