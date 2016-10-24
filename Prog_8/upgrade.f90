      Subroutine Upgrade(GR,N_op,NT,PHASE,Op_dim) 

!!!!!! This version of  Upgrade.f90  contains optimization carried out by Johannes Hofmann
!!!!!! The original version of this routine can be found in upgrade_FFA.f90 
!!!!!! Both versions  must give the same results

       
        Use Hamiltonian
        Use Random_wrap
        Use Control
        Use Precdef
        Implicit none 
        
        Complex (Kind=double) :: GR(Ndim,Ndim, N_FL)
        Integer, INTENT(IN) :: N_op, Nt, Op_dim
        Complex (Kind=double) :: Phase

        ! Local ::
        Complex (Kind=double) :: Mat(Op_dim,Op_Dim), Delta(Op_dim,N_FL)
        Complex (Kind=double) :: Ratio(N_FL), Ratiotot, Z1 
        Integer :: ns_new, ns_old, n,m,nf, i,j
        Complex (Kind= double) :: ZK, Z, D_Mat, Z2, myexp, tmpsxv, tmpsyu
        Integer, external :: nranf
        
        Real (Kind = double) :: Weight, reZ, imZ
        Complex (Kind = double) :: u(Ndim,Op_dim), v(Ndim,Op_dim) ,alpha, beta
        Complex (Kind = double) :: y_v(Ndim,Op_dim), xp_v(Ndim,Op_dim)
        Complex (Kind = double) :: x_v(Ndim,Op_dim)
        Logical :: Log
        Complex (Kind = Kind(0.D0)), Dimension(:, :), Allocatable :: Zarr, grarr
        Complex (Kind = Kind(0.D0)), Dimension(:), Allocatable :: sxv, syu

        if ( abs(OP_V(n_op,1)%g) < 1.D-6 ) return

        ! Compute the ratio
        nf = 1
        ns_old = nsigma(n_op,nt)
        If ( Op_V(n_op,nf)%type == 1) then
           ns_new = -ns_old
        else
           ns_new = NFLIPL(Ns_old,nranf(3))
        endif
        Do nf = 1,N_FL
           Z1 = Op_V(n_op,nf)%g * ( Phi(ns_new,Op_V(n_op,nf)%type) -  Phi(ns_old,Op_V(n_op,nf)%type))
           Do m = 1,Op_V(n_op,nf)%N_non_zero
              myexp = exp( Z1* Op_V(n_op,nf)%E(m) )
              Z = myexp - 1.d0
              Delta(m,nf) = Z
              do n = 1,Op_V(n_op,nf)%N_non_zero
              Mat(n,m) = - Z * GR( Op_V(n_op,nf)%P(n), Op_V(n_op,nf)%P(m),nf )
              Enddo
              Mat(m,m) = myexp + Mat(m,m)
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
        
        Ratiotot = Product(Ratio)
        nf = 1
        Ratiotot = (Ratiotot**dble(N_SUN)) * Gaml(ns_new, Op_V(n_op,nf)%type)/Gaml(ns_old, Op_V(n_op,nf)%type)
        Ratiotot = Ratiotot * real(S0(n_op,nt), kind(0.D0))!Just to be save since S0 seems to be user supplied
        

        !Write(6,*) Ratiotot
        
        Weight = abs(  real(Phase * Ratiotot, kind=double)/real(Phase,kind=double) )
      
        Log = .false. 
        if ( Weight > ranf() )  Then
           Log = .true.
           Phase = Phase * Ratiotot/weight
           !Write(6,*) 'Accepted : ', Ratiotot

           Do nf = 1,N_FL
              ! Setup u(i,n), v(n,i) 
              beta = 0.D0
              call zlaset('N', Ndim, Op_dim, beta, beta, u, size(u, 1))
              call zlaset('N', Ndim, Op_dim, beta, beta, v, size(v, 1))
              do n = 1,Op_V(n_op,nf)%N_non_zero
                 u( Op_V(n_op,nf)%P(n), n) = Delta(n,nf)
                 do i = 1,Ndim
                    v(i,n) = - GR( Op_V(n_op,nf)%P(n), i, nf )
                 enddo
                 v(Op_V(n_op,nf)%P(n), n)  = 1.d0 - GR( Op_V(n_op,nf)%P(n),  Op_V(n_op,nf)%P(n), nf)
              enddo

              call zlaset('N', Ndim, Op_dim, beta, beta, x_v, size(x_v, 1))
              call zlaset('N', Ndim, Op_dim, beta, beta, y_v, size(y_v, 1))
              i = Op_V(n_op,nf)%P(1)
              x_v(i, 1) = u(i, 1)/(1.d0 + v(i,1)*u(i,1) )
              call zcopy(Ndim, v(:, 1), 1, y_v(:, 1), 1)
              do n = 2,Op_V(n_op,nf)%N_non_zero
                 call zcopy(Ndim, u(:, n), 1, x_v(:, n), 1)
                 call zcopy(Ndim, v(:, n), 1, y_v(:, n), 1)
                 Z = 1.d0 + u( Op_V(n_op,nf)%P(n), n)*v(Op_V(n_op,nf)%P(n),n)
                 alpha = -1.D0
                 Allocate(syu(n), sxv(n))
                 call zgemv('T', NDim, n-1, alpha, y_v, Ndim, u(1,n), 1, beta , syu, 1)
                 call zgemv('T', NDim, n-1, alpha, x_v, Ndim, v(1,n), 1, beta , sxv, 1)
                 do m = 1,n-1
                    Z = Z - syu(m)*sxv(m)
                    call zaxpy(Ndim, syu(m), x_v(1, m), 1, x_v(1, n), 1)
                    call zaxpy(Ndim, sxv(m), y_v(1, m), 1, y_v(1, n), 1)
                 enddo
                 Z = 1.D0/Z
                 call zscal(Ndim, Z, x_v(1, n), 1)
                 Deallocate(syu, sxv)
              enddo
              Allocate (Zarr(Op_dim,Op_dim), grarr(NDim, Op_dim))
              Zarr = x_v(Op_V(n_op,nf)%P, :)
              grarr = gr(:, Op_V(n_op,nf)%P, nf)
              alpha = 1.D0
              CALL ZGEMM('N', 'N', NDim, Op_Dim, Op_Dim, alpha, grarr, size(grarr,1), Zarr, size(Zarr,1), beta, xp_v, size(xp_v,1))
              Deallocate(Zarr, grarr)
              !do n = 1,Op_dim
              !   do j = 1,Ndim
              !      do i = 1,Ndimop
              !         gr(i,j,nf) = gr(i,j,nf) - xp_v(i,n)*y_v(j,n)
              !      enddo
              !   enddo
              !enddo
              ! gr(:,:,nf) -= xp_v(:,:) * y_v(:,:)^T
              ! Replace by Zgemm 
              alpha = cmplx (-1.0d0, 0.0d0, kind(0.D0))
              beta  = cmplx ( 1.0d0, 0.0d0, kind(0.D0))
              CALL ZGEMM('N','T',Ndim,Ndim,Op_dim,alpha,xp_v,size(xp_v,1),y_v,size(y_v,1),beta,gr(1,1,nf),size(gr,1))


!!!!!         Requires additional space
!             Complex (Kind = double) ::  tmpMat(Ndim,Ndim), tmp
!           
! 	      !$OMP PARALLEL DO PRIVATE(tmp)
!               do j = 1,Ndim
!                  do i = 1,Ndim
! 		    tmp=cmplx(0.d0,0.d0)
! 		    do n = 1,Op_dim
! 		       tmp = tmp - xp_v(i,n)*y_v(j,n)
!                     enddo
!                     if (abs(tmpMat(i,j)-tmp) >= 0.00001) then
!                       write(*,*) tmpMat(i,j), tmp, abs(tmpMat(i,j)-tmp)
!                     else
!                       write(*,*) "OK"
!                     endif
! 		    gr(i,j,nf) = gr(i,j,nf) + tmp
!                  enddo
!               enddo
! 	      !$OMP END PARALLEL DO

           enddo
           
           ! Flip the spin
           nsigma(n_op,nt) = ns_new
        endif
        Call Control_upgrade(Log)

      End Subroutine Upgrade
