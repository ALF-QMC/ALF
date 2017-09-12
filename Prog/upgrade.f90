!  Copyright (C) 2016, 2017 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

      Subroutine Upgrade(GR,N_op,NT,PHASE,Op_dim,Propose_S0) 

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This routine updates the field associated to the operator N_op on time 
!> slice NT. If  the local flip is accepted, the Green function is updated.
!
!--------------------------------------------------------------------
       
        Use Hamiltonian
        Use Random_wrap
        Use Control
        Implicit none 
        
        Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL)
        Integer, INTENT(IN)       :: N_op, Nt, Op_dim
        Complex (Kind=Kind(0.d0)) :: Phase
        LOGICAL, INTENT(IN)       :: Propose_S0


        ! Local ::
        Complex (Kind=Kind(0.d0)) :: Mat(Op_dim,Op_Dim), Delta(Op_dim,N_FL)
        Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Ratiotot, Z1 
        Integer :: ns_new, ns_old, n,m,nf, i,j
        Complex (Kind=Kind(0.d0)) :: ZK, Z, D_Mat, Z2, myexp, s1, s2
        
        Real    (Kind=Kind(0.d0)) :: Weight, reZ, imZ
        Complex (Kind=Kind(0.d0)) :: u(Ndim,Op_dim), v(Ndim,Op_dim) ,alpha, beta
        Complex (Kind=Kind(0.d0)) :: y_v(Ndim,Op_dim), xp_v(Ndim,Op_dim)
        Complex (Kind=Kind(0.d0)) :: x_v(Ndim,Op_dim)
        Logical :: toggle
        Complex (Kind=Kind(0.D0)), Dimension(:, :), Allocatable :: Zarr, grarr
        Complex (Kind=Kind(0.D0)), Dimension(:), Allocatable :: sxv, syu

        if ( abs(OP_V(n_op,1)%g) < 1.D-6 ) return

        ! Compute the ratio
        nf = 1
        ns_old = nsigma(n_op,nt)
        If ( Op_V(n_op,nf)%type == 1) then
           if ( Propose_S0 ) then
              Weight = 1.d0 - 1.d0/(1.d0+S0(n_op,nt))
              If ( Weight < ranf_wrap() ) then
                 Call Control_upgrade_eff(.false.)
                 Return
              endif
           endif
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
              s1 = Mat(1,1)*Mat(2,2)
              s2 = Mat(2,1)*Mat(1,2)
              If (Abs(s1) > Abs(s2)) then
                D_mat = s1*(1.D0 - s2/s1)
              else
                D_mat = s2*(s1/s2 - 1.D0)
              Endif
!              D_mat =  - Mat(2,1)*Mat(1,2)
           else
              D_mat = Det(Mat,Size(Mat,1))
           endif
           Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
        Enddo
        
        Ratiotot = Product(Ratio)
        nf = 1
        Ratiotot = (Ratiotot**dble(N_SUN)) * Gaml(ns_new, Op_V(n_op,nf)%type)/Gaml(ns_old, Op_V(n_op,nf)%type)
        if ( .not. Propose_S0 ) &
             &  Ratiotot = Ratiotot * real(S0(n_op,nt), kind(0.D0))  ! Just to be safe since S0 seems to be user supplied
        

        !Write(6,*) Ratiotot
        
        Weight = abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )
      
        toggle = .false. 
        if ( Weight > ranf_wrap() )  Then
           toggle = .true.
           Phase = Phase * Ratiotot/abs(Ratiotot)
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
                 alpha = 1.D0
                 call zgemv('N', NDim, n-1, alpha, x_v, Ndim, syu, 1, alpha, x_v(1, n), 1)
                 call zgemv('N', NDim, n-1, alpha, y_v, Ndim, sxv, 1, alpha, y_v(1, n), 1)
                 do m = 1,n-1
                    Z = Z - syu(m)*sxv(m)
                 enddo
                 Z = 1.D0/Z
                 call zscal(Ndim, Z, x_v(1, n), 1)
                 Deallocate(syu, sxv)
              enddo
              IF (Op_dim == 1) THEN
                CALL ZCOPY(Ndim, gr(1, Op_V(n_op,nf)%P(1), nf), 1, xp_v(1, 1), 1)
                CALL ZGERU(Ndim, Ndim, -x_v(Op_V(n_op,nf)%P(1), 1), xp_v(1,1), 1, y_v(1, 1), 1, gr(1,1,nf), Ndim)
              ELSE
                Allocate (Zarr(Op_dim,Op_dim), grarr(NDim, Op_dim))
                Zarr = x_v(Op_V(n_op,nf)%P, :)
                grarr = gr(:, Op_V(n_op,nf)%P, nf)
                alpha = 1.D0
                CALL ZGEMM('N', 'N', NDim, Op_Dim, Op_Dim, alpha, grarr, size(grarr,1), Zarr, Op_Dim, beta, xp_v, size(xp_v,1))
                Deallocate(Zarr, grarr)
                beta  = cmplx ( 1.0d0, 0.0d0, kind(0.D0))
                alpha = -1.D0
                CALL ZGEMM('N','T',Ndim,Ndim,Op_dim,alpha,xp_v, Ndim,y_v, Ndim,beta,gr(1,1,nf), Ndim)
              ENDIF
              !do n = 1,Op_dim
              !   do j = 1,Ndim
              !      do i = 1,Ndimop
              !         gr(i,j,nf) = gr(i,j,nf) - xp_v(i,n)*y_v(j,n)
              !      enddo
              !   enddo
              !enddo
              ! gr(:,:,nf) -= xp_v(:,:) * y_v(:,:)^T
              ! Replace by Zgemm 

!!!!!         Requires additional space
!             Complex (Kind =Kind(0.d0)) ::  tmpMat(Ndim,Ndim), tmp
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

        Call Control_upgrade(toggle)
        Call Control_upgrade_eff(toggle)

      End Subroutine Upgrade

!--------------------------------------------------------------------
      Subroutine Upgrade2(GR,N_op,NT,PHASE,Op_dim,ns_new, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, toggle,  mode) 
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This routine updates the field associated to the operator N_op on time 
!> slice NT to the value ns_new
!> if mode = final  the move is  accepted according to T0_proposal_ratio*S0_ratio*Prev_Ratiotot*ratio  and the Green function is updated
!> if mode = intermediate the move is carried our deterministically  and the Green function updated.  Also Prev_Ratio = Prev_Ration*ratio
!> The ratio is computed in the routine.
!--------------------------------------------------------------------
       
        Use Hamiltonian
        Use Random_wrap
        Use Control
        Implicit none 
        
        Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL)
        Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot
        Integer, INTENT(IN)       :: N_op, Nt, Op_dim
        Complex (Kind=Kind(0.d0)) :: Phase
        Integer                   :: ns_new
        Real    (Kind=Kind(0.d0)) :: S0_ratio, T0_proposal_ratio
        Character (Len=64)        :: Mode
        Logical                   :: toggle

        
        ! Local ::
        Complex (Kind=Kind(0.d0)) :: Mat(Op_dim,Op_Dim), Delta(Op_dim,N_FL)
        Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Ratiotot, Z1 
        Integer :: ns_old, n,m,nf, i,j
        Complex (Kind=Kind(0.d0)) :: ZK, Z, D_Mat, Z2, myexp, s1, s2
        
        Real    (Kind=Kind(0.d0)) :: Weight, reZ, imZ
        Complex (Kind=Kind(0.d0)) :: u(Ndim,Op_dim), v(Ndim,Op_dim) ,alpha, beta
        Complex (Kind=Kind(0.d0)) :: y_v(Ndim,Op_dim), xp_v(Ndim,Op_dim)
        Complex (Kind=Kind(0.d0)) :: x_v(Ndim,Op_dim)
        Complex (Kind=Kind(0.D0)), Dimension(:, :), Allocatable :: Zarr, grarr
        Complex (Kind=Kind(0.D0)), Dimension(:), Allocatable :: sxv, syu

        toggle = .false.
        ! if ( abs(OP_V(n_op,1)%g) < 1.D-12 )   return

        ! Compute the ratio
        nf = 1
        ns_old = nsigma(n_op,nt)
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
              s1 = Mat(1,1)*Mat(2,2)
              s2 = Mat(2,1)*Mat(1,2)
              If (Abs(s1) > Abs(s2)) then
                D_mat = s1*(1.D0 - s2/s1)
              else
                D_mat = s2*(s1/s2 - 1.D0)
              Endif
              !  D_mat =  - Mat(2,1)*Mat(1,2)
           else
              D_mat = Det(Mat,Size(Mat,1))
           endif
           Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
        Enddo
        
        Ratiotot = Product(Ratio)
        nf = 1
        Ratiotot = (Ratiotot**dble(N_SUN)) * Gaml(ns_new, Op_V(n_op,nf)%type)/Gaml(ns_old, Op_V(n_op,nf)%type)

        !Write(6,*) Ratiotot
        
        If      (mode  == "Final"       ) Then
           !Write(6,*) "Up_I: ", Ratiotot
           Ratiotot = Ratiotot * Prev_Ratiotot
           !Write(6,*) "Up_f: ", Ratiotot
           Weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )
           !Write(6,*) Phase, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, ns_old,ns_new
        elseif  (mode == "Intermediate" ) Then
           Weight = 1.5
           !Write(6,*) "Up_I: ", Ratiotot
           Prev_Ratiotot = Prev_Ratiotot*Ratiotot
        else
           Write(6,*) 'Error'
           stop
        endif

        toggle = .false. 
        if ( Weight > ranf_wrap() )  Then
           toggle = .true.
           If (mode == "Final"  )  Phase = Phase * Ratiotot/sqrt(Ratiotot*conjg(Ratiotot))
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
                 alpha = 1.D0
                 call zgemv('N', NDim, n-1, alpha, x_v, Ndim, syu, 1, alpha, x_v(1, n), 1)
                 call zgemv('N', NDim, n-1, alpha, y_v, Ndim, sxv, 1, alpha, y_v(1, n), 1)
                 do m = 1,n-1
                    Z = Z - syu(m)*sxv(m)
                 enddo
                 Z = 1.D0/Z
                 call zscal(Ndim, Z, x_v(1, n), 1)
                 Deallocate(syu, sxv)
              enddo
              IF (size(Op_V(n_op,nf)%P, 1) == 1) THEN
                CALL ZCOPY(Ndim, gr(1, Op_V(n_op,nf)%P(1), nf), 1, xp_v(1, 1), 1)
                Z = -x_v(Op_V(n_op,nf)%P(1), 1)
                CALL ZGERU(Ndim, Ndim, Z, xp_v(1,1), 1, y_v(1, 1), 1, gr(1,1,nf), Ndim)
              ELSE
                Allocate (Zarr(size(Op_V(n_op,nf)%P, 1), Op_dim), grarr(NDim, Op_dim))
                Zarr = x_v(Op_V(n_op,nf)%P, :)
                grarr = gr(:, Op_V(n_op,nf)%P, nf)
                alpha = 1.D0
                CALL ZGEMM('N', 'N', NDim, Op_Dim, Op_Dim, alpha, grarr, Ndim, Zarr, size(Op_V(n_op,nf)%P, 1), beta, xp_v, Ndim)
                Deallocate(Zarr, grarr)
                beta  = cmplx ( 1.0d0, 0.0d0, kind(0.D0))
                alpha = -1.D0
                CALL ZGEMM('N','T',Ndim,Ndim,Op_dim,alpha,xp_v, Ndim, y_v, Ndim,beta,gr(1,1,nf), Ndim)
              ENDIF

           enddo
           
           ! Flip the spin
           nsigma(n_op,nt) = ns_new
        endif

        If ( mode == "Final" )  then
           Call Control_upgrade(toggle)
           Call Control_upgrade_eff(toggle)
        endif

      End Subroutine Upgrade2
