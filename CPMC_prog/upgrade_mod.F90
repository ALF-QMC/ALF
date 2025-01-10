module upgrade_mod
   use runtime_error_mod
   implicit none
contains

   subroutine upgrade_ph(GR, N_op, PHASE, Hs_new, i_wlk, log_I_ratio)

      use Hamiltonian_main
      use Random_wrap
      use Control
      use Fields_mod
      use Operator_mod
      use iso_fortran_env, only: output_unit, error_unit
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT) :: GR(Ndim, Ndim, N_FL)
      integer, intent(IN)    :: N_op, i_wlk
      complex(Kind=kind(0.d0)), intent(INOUT) :: Phase
      complex(Kind=kind(0.d0)), intent(OUT)   :: Hs_new
      complex(Kind=kind(0.d0)), intent(INOUT) :: log_I_ratio

      ! Local ::
      type(Fields)   ::  nsigma_new
      complex(Kind=kind(0.d0)) :: Ratio(N_FL), Ratio_f(N_FL), Ratiotot, Z1, Phase_a_array(N_FL)
      integer ::  n, m, nf, nf_eff, i, Op_dim, op_dim_nf, nu_spin, nu_c, n_prop, ratioD
      complex(Kind=kind(0.d0)) :: Z, D_Mat, myexp, s1, s2, gamma_ratio
      real(Kind=kind(0.d0)) :: S0_ratio

      real(Kind=kind(0.d0)) :: Weight, st_r, ed_r, sum_ratio, rand_nu, x_field
      complex(Kind=kind(0.d0)) :: alpha, beta, g_loc
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: Mat, Delta
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: u, v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: y_v, xp_v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: x_v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: Zarr, grarr
      complex(Kind=kind(0.d0)), dimension(:), allocatable :: sxv, syu

      real(Kind=kind(0.d0)), dimension(:), allocatable :: field_list
      complex(Kind=kind(0.d0)) :: ratio_field, ratio_O, ratio_G

      call nsigma_new%make(1, 1)

      op_dim = op_v(n_op, calc_fl_map(1))%n_non_zero
      do nf_eff = 2, n_fl_eff
         nf = calc_fl_map(nf_eff)
         if (op_dim < op_v(n_op, nf)%n_non_zero) op_dim = op_v(n_op, nf)%n_non_zero
      end do

      if (op_dim > 0) then
         allocate (mat(op_dim, op_dim), delta(op_dim, n_fl_eff), u(ndim, op_dim), v(ndim, op_dim))
         allocate (y_v(ndim, op_dim), xp_v(ndim, op_dim), x_v(ndim, op_dim))
      end if

      ! Compute the ratio for each posibility
      nf = 1
      nsigma_new%t(1) = op_v(n_op, nf)%type

      !--------------------------------------------------!
      !! Decide the field in the next propagation
      x_field = rang_wrap()
      gamma_ratio = cmplx(x_field, 0.d0, kind(0.d0))*x_local(N_op) - 0.5d0*x_local(N_op)**2
      nsigma_new%f(1, 1) = cmplx(x_field, 0.d0, kind(0.d0)) - x_local(N_op)

      ! Output the field
      Hs_new = nsigma_new%f(1, 1)

      S0_ratio = ham%S0(n_op, 1, Hs_new)

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         g_loc = Op_V(n_op, nf)%g
         Z1 = g_loc*nsigma_new%Phi(1, 1)
         op_dim_nf = Op_V(n_op, nf)%N_non_zero
         do m = 1, op_dim_nf
            myexp = exp(Z1*Op_V(n_op, nf)%E(m))
            Z = myexp - 1.d0
            Delta(m, nf_eff) = Z
            do n = 1, op_dim_nf
               Mat(n, m) = -Z*GR(Op_V(n_op, nf)%P(n), Op_V(n_op, nf)%P(m), nf)
            end do
            Mat(m, m) = myexp + Mat(m, m)
         end do
         if (op_dim_nf == 0) then
            D_mat = 1.0d0
         elseif (op_dim_nf == 1) then
            D_mat = Mat(1, 1)
         elseif (op_dim_nf == 2) then
            s1 = Mat(1, 1)*Mat(2, 2)
            s2 = Mat(2, 1)*Mat(1, 2)
            if (abs(s1) > abs(s2)) then
               D_mat = s1*(1.d0 - s2/s1)
            else
               D_mat = s2*(s1/s2 - 1.d0)
            end if
         else
            D_mat = Det(Mat, op_dim_nf)
         end if
         Ratio(nf) = D_Mat*exp(Z1*Op_V(n_op, nf)%alpha)
         Ratio_f(nf) = D_Mat
      end do

      !call reconstruct weight subroutine to fill the non-calculated blocks
      if (reconstruction_needed) call ham%weight_reconstruction(Ratio)

      Ratiotot = product(Ratio)
      Ratiotot = (Ratiotot**dble(N_SUN))

      ratio_O = Ratiotot
      ratio_G = product(Ratio_f)**dble(N_SUN)
      !--------------------------------------------------!

        !! Compute the phase of the weight
      !Phase = Phase * Ratio_g/sqrt(Ratio_g*conjg(Ratio_g))

        !! update ratio and overlap
      Overlap(i_wlk) = Overlap(i_wlk)*ratio_O
      log_I_ratio = log_I_ratio + log(ratio_O) + gamma_ratio

        !! Update Green's function
      ! update delta
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         g_loc = Op_V(n_op, nf)%g
         Z1 = g_loc*(nsigma_new%Phi(1, 1))
         op_dim_nf = Op_V(n_op, nf)%N_non_zero
         do m = 1, op_dim_nf
            myexp = exp(Z1*Op_V(n_op, nf)%E(m))
            Z = myexp - 1.d0
            Delta(m, nf_eff) = Z
         end do
      end do

      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         ! Setup u(i,n), v(n,i)
         op_dim_nf = Op_V(n_op, nf)%N_non_zero
         if (op_dim_nf > 0) then
            beta = 0.d0
            call zlaset('N', Ndim, op_dim_nf, beta, beta, u, size(u, 1))
            call zlaset('N', Ndim, op_dim_nf, beta, beta, v, size(v, 1))
            do n = 1, op_dim_nf
               u(Op_V(n_op, nf)%P(n), n) = Delta(n, nf_eff)
               do i = 1, Ndim
                  v(i, n) = -GR(Op_V(n_op, nf)%P(n), i, nf)
               end do
               v(Op_V(n_op, nf)%P(n), n) = 1.d0 - GR(Op_V(n_op, nf)%P(n), Op_V(n_op, nf)%P(n), nf)
            end do

            call zlaset('N', Ndim, op_dim_nf, beta, beta, x_v, size(x_v, 1))
            call zlaset('N', Ndim, op_dim_nf, beta, beta, y_v, size(y_v, 1))
            i = Op_V(n_op, nf)%P(1)
            x_v(i, 1) = u(i, 1)/(1.d0 + v(i, 1)*u(i, 1))
            call zcopy(Ndim, v(:, 1), 1, y_v(:, 1), 1)
            do n = 2, op_dim_nf
               call zcopy(Ndim, u(:, n), 1, x_v(:, n), 1)
               call zcopy(Ndim, v(:, n), 1, y_v(:, n), 1)
               Z = 1.d0 + u(Op_V(n_op, nf)%P(n), n)*v(Op_V(n_op, nf)%P(n), n)
               alpha = -1.d0
               allocate (syu(n), sxv(n))
               call zgemv('T', NDim, n - 1, alpha, y_v, Ndim, u(1, n), 1, beta, syu, 1)
               call zgemv('T', NDim, n - 1, alpha, x_v, Ndim, v(1, n), 1, beta, sxv, 1)
               alpha = 1.d0
               call zgemv('N', NDim, n - 1, alpha, x_v, Ndim, syu, 1, alpha, x_v(1, n), 1)
               call zgemv('N', NDim, n - 1, alpha, y_v, Ndim, sxv, 1, alpha, y_v(1, n), 1)
               do m = 1, n - 1
                  Z = Z - syu(m)*sxv(m)
               end do
               Z = 1.d0/Z
               call zscal(Ndim, Z, x_v(1, n), 1)
               deallocate (syu, sxv)
            end do
            if (size(Op_V(n_op, nf)%P, 1) == 1) then
               call ZCOPY(Ndim, gr(1, Op_V(n_op, nf)%P(1), nf), 1, xp_v(1, 1), 1)
               Z = -x_v(Op_V(n_op, nf)%P(1), 1)
               call ZGERU(Ndim, Ndim, Z, xp_v(1, 1), 1, y_v(1, 1), 1, gr(1, 1, nf), Ndim)
            else
               allocate (zarr(op_dim_nf, op_dim_nf), grarr(NDim, op_dim_nf))
               Zarr = x_v(Op_V(n_op, nf)%P(1:op_dim_nf), :)
               grarr = gr(:, Op_V(n_op, nf)%P(1:op_dim_nf), nf)
               beta = 0.d0
               alpha = 1.d0
               call ZGEMM('N', 'N', NDim, op_dim_nf, op_dim_nf, alpha, grarr, Ndim, Zarr, op_dim_nf, beta, xp_v, Ndim)
               deallocate (Zarr, grarr)
               beta = cmplx(1.0d0, 0.0d0, kind(0.d0))
               alpha = -1.d0
               call ZGEMM('N', 'T', Ndim, Ndim, op_dim_nf, alpha, xp_v, Ndim, y_v, Ndim, beta, gr(1, 1, nf), Ndim)
            end if
         end if

      end do

      if (op_dim > 0) then
         deallocate (Mat, Delta, u, v)
         deallocate (y_v, xp_v, x_v)
      end if

      call nsigma_new%clear()

   end subroutine Upgrade_ph

end module upgrade_mod
