module upgrade_mod
   use runtime_error_mod
   implicit none
   contains

      Subroutine Upgrade(GR,N_op,Hs_new,i_wlk)

        Use Hamiltonian_main
        Use Random_wrap
        Use Control
        Use Fields_mod
        Use Operator_mod
        use iso_fortran_env, only: output_unit, error_unit
        Implicit none

        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim,N_FL,N_grc)
        Integer                  , INTENT(IN)    :: N_op, i_wlk
        Real    (Kind=Kind(0.d0)), INTENT(OUT)   :: Hs_new

        ! Local ::
        Type   (Fields)   ::  nsigma_new
        Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Ratio_f(N_FL), ratiotot, Z1, exp_o_new, exp_o_old
        Integer ::  n,ns,m,nf, nf_eff, i, Op_dim, op_dim_nf, nu_spin, nu_c, n_prop, i_grc
        Complex (Kind=Kind(0.d0)) :: Z, D_Mat, myexp, s1, s2, ratioD
        Real    (Kind=Kind(0.d0)) :: S0_ratio

        Real    (Kind=Kind(0.d0)) :: Weight, st_r, ed_r, sum_ratio, rand_nu
        Complex (Kind=Kind(0.d0)) :: alpha, beta, g_loc
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: Mat, Delta
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: u, v
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: y_v, xp_v
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: x_v
        Complex (Kind=Kind(0.D0)), Dimension(:, :), Allocatable :: Zarr, grarr
        Complex (Kind=Kind(0.D0)), Dimension(:), Allocatable :: sxv, syu
        
        Real    (Kind=Kind(0.D0)), Dimension(:), Allocatable :: ratio_field, field_list
        Complex (Kind=Kind(0.D0)), Dimension(:,:), Allocatable :: ratio_o

        Call nsigma_new%make(1,1)
        
        op_dim=Op_V(n_op,Calc_FL_map(1))%N_non_zero
        Do nf_eff = 2,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           if (op_dim<Op_V(n_op,nf)%N_non_zero) op_dim=Op_V(n_op,nf)%N_non_zero
        Enddo
        
        if (op_dim > 0) then
           Allocate ( Mat(Op_dim,Op_Dim), Delta(Op_dim,N_FL_eff), u(Ndim,Op_dim), v(Ndim,Op_dim) )
           Allocate ( y_v(Ndim,Op_dim), xp_v(Ndim,Op_dim), x_v(Ndim,Op_dim) )
        endif

        ! Compute the ratio for each posibility
        nf = 1
        nsigma_new%t(1)    = Op_V(n_op,nf)%Type
        if ( Op_V(n_op,nf)%Type .eq. 1 ) then
            allocate(field_list(2))
            allocate(ratio_field(2))
            allocate(ratio_o(2,N_slat))
            field_list(1)=1.d0; 
            field_list(2)=2.d0;
            nu_spin=2
        elseif ( Op_V(n_op,nf)%Type .eq. 2 ) then
            allocate(field_list(4))
            allocate(ratio_field(4))
            allocate(ratio_o(4,N_slat))
            field_list(1)= 1.d0; 
            field_list(2)=-1.d0; 
            field_list(3)= 2.d0; 
            field_list(4)=-2.d0; 
            nu_spin=4
        endif

        do nu_c = 1, nu_spin

            nsigma_new%f(1,1)  = field_list(nu_c) !real(ns_new,kind=kind(0.d0))
            S0_ratio           = ham%S0(n_op,1,dble(field_list(nu_c)))
            
            exp_o_new = 0.d0
            exp_o_old = 0.d0
            do ns = 1, N_slat
               i_grc = ns+(i_wlk-1)*N_slat
               
               do nf_eff = 1,N_FL_eff
                  nf=Calc_Fl_map(nf_eff)
                  g_loc = Op_V(n_op,nf)%g
                  Z1 = g_loc * ( nsigma_new%Phi(1,1) )
                  op_dim_nf = Op_V(n_op,nf)%N_non_zero
                  Do m = 1,op_dim_nf
                     myexp = exp( Z1* Op_V(n_op,nf)%E(m) )
                     Z = myexp - 1.d0
                     Delta(m,nf_eff) = Z
                     do n = 1,op_dim_nf
                        Mat(n,m) = - Z * GR( Op_V(n_op,nf)%P(n), Op_V(n_op,nf)%P(m),nf, i_grc )
                     Enddo
                     Mat(m,m) = myexp + Mat(m,m)
                  Enddo
                  If (op_dim_nf == 0 ) then
                     D_mat = 1.0d0
                  elseIf (op_dim_nf == 1 ) then
                     D_mat = Mat(1,1)
                  elseif (op_dim_nf == 2 ) then
                     s1 = Mat(1,1)*Mat(2,2)
                     s2 = Mat(2,1)*Mat(1,2)
                     If (Abs(s1) > Abs(s2)) then
                       D_mat = s1*(1.D0 - s2/s1)
                     else
                       D_mat = s2*(s1/s2 - 1.D0)
                     Endif
                  else
                     D_mat = Det(Mat,op_dim_nf)
                  endif
                  Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
               Enddo

               !call reconstruct weight subroutine to fill the non-calculated blocks
               if (reconstruction_needed) call ham%weight_reconstruction(Ratio)
            
               ratiotot = Product(Ratio)**dble(N_SUN)
               ratio_o(nu_c,ns) = ratiotot

               exp_o_new = exp_o_new + exp(overlap(i_grc))*ratio_o(nu_c,ns)
               exp_o_old = exp_o_old + exp(overlap(i_grc))

            enddo

            ratiotot = exp_o_new/exp_o_old
            weight = S0_ratio * dble(ratiotot*nsigma_new%px0(1,1))
            if ( weight .le. 0.d0 ) weight = 0.d0
            ratio_field(nu_c) = weight

        enddo

        sum_ratio = sum(ratio_field)

        if ( sum_ratio .le. 0.d0  ) then
            weight_k(i_wlk) = 0.d0
            Hs_new = field_list(1) ! randomly set up a Hs_new for output
        else
            !! Decide the field in the next propagation
            n_prop = -1
            rand_nu = ranf_wrap()
            st_r = 0.d0
            do nu_c = 1, nu_spin
                ed_r = st_r + ratio_field(nu_c)/sum_ratio
                if ( (rand_nu .ge. st_r ) .and. (rand_nu .lt. ed_r )  ) then
                    n_prop = nu_c
                    exit
                endif
                st_r = st_r + ratio_field(nu_c)/sum_ratio
            enddo
            nsigma_new%f(1,1)  = field_list(n_prop)
            
            !! update weight
            weight_k(i_wlk) = weight_k(i_wlk)*sum_ratio

            !! Update Green's function
            ! update delta
            Do nf_eff = 1,N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               g_loc = Op_V(n_op,nf)%g
               Z1 = g_loc * ( nsigma_new%Phi(1,1) )
               op_dim_nf = Op_V(n_op,nf)%N_non_zero
               Do m = 1,op_dim_nf
                  myexp = exp( Z1* Op_V(n_op,nf)%E(m) )
                  Z = myexp - 1.d0
                  Delta(m,nf_eff) = Z
               Enddo
            Enddo

            do ns = 1, N_slat
               
               i_grc = ns+(i_wlk-1)*N_slat
               !! update overlap
               overlap (i_grc) = overlap (i_grc) + log(ratio_o(n_prop,ns))

               Do nf_eff = 1,N_FL_eff
                  nf=Calc_Fl_map(nf_eff)
                  ! Setup u(i,n), v(n,i)
                  op_dim_nf = Op_V(n_op,nf)%N_non_zero
                  if (op_dim_nf > 0) then
                    beta = 0.D0
                    call zlaset('N', Ndim, op_dim_nf, beta, beta, u, size(u, 1))
                    call zlaset('N', Ndim, op_dim_nf, beta, beta, v, size(v, 1))
                    do n = 1,op_dim_nf
                        u( Op_V(n_op,nf)%P(n), n) = Delta(n,nf_eff)
                        do i = 1,Ndim
                            v(i,n) = - GR( Op_V(n_op,nf)%P(n), i, nf, i_grc )
                        enddo
                        v(Op_V(n_op,nf)%P(n), n)  = 1.d0 - GR( Op_V(n_op,nf)%P(n),  Op_V(n_op,nf)%P(n), nf, i_grc)
                    enddo

                    call zlaset('N', Ndim, op_dim_nf, beta, beta, x_v, size(x_v, 1))
                    call zlaset('N', Ndim, op_dim_nf, beta, beta, y_v, size(y_v, 1))
                    i = Op_V(n_op,nf)%P(1)
                    x_v(i, 1) = u(i, 1)/(1.d0 + v(i,1)*u(i,1) )
                    call zcopy(Ndim, v(:, 1), 1, y_v(:, 1), 1)
                    do n = 2,op_dim_nf
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
                        CALL ZCOPY(Ndim, gr(1, Op_V(n_op,nf)%P(1), nf, i_grc), 1, xp_v(1, 1), 1)
                        Z = -x_v(Op_V(n_op,nf)%P(1), 1)
                        CALL ZGERU(Ndim, Ndim, Z, xp_v(1,1), 1, y_v(1, 1), 1, gr(1,1,nf,i_grc), Ndim)
                    ELSE
                        Allocate (zarr(op_dim_nf, op_dim_nf), grarr(NDim, op_dim_nf))
                        Zarr = x_v(Op_V(n_op,nf)%P(1:op_dim_nf), :)
                        grarr = gr(:, Op_V(n_op,nf)%P(1:op_dim_nf), nf, i_grc)
                        beta  = 0.d0
                        alpha = 1.D0
                        CALL ZGEMM('N', 'N', NDim, op_dim_nf, op_dim_nf, alpha, grarr, Ndim, Zarr, op_dim_nf, beta, xp_v, Ndim)
                        Deallocate(Zarr, grarr)
                        beta  = cmplx ( 1.0d0, 0.0d0, kind(0.D0))
                        alpha = -1.D0
                        CALL ZGEMM('N','T',Ndim,Ndim,op_dim_nf,alpha,xp_v, Ndim, y_v, Ndim,beta,gr(1,1,nf,i_grc), Ndim)
                    ENDIF
                  endif
               enddo

            Enddo

            ! Output the field
            Hs_new = nsigma_new%f(1,1) ! real(ns_new,Kind=kind(0.d0))

        endif
        
        if (op_dim > 0) then
           deallocate ( Mat, Delta, u, v )
           deallocate ( y_v, xp_v, x_v )
           deallocate ( ratio_field, field_list, ratio_o )
        endif

        Call nsigma_new%clear()

      End Subroutine Upgrade

end module upgrade_mod
