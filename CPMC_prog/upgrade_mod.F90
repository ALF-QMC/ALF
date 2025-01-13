module upgrade_mod
   use runtime_error_mod
   implicit none
contains

   subroutine Upgrade(GR, N_op, Hs_new, i_wlk)

      use Hamiltonian_main
      use Random_wrap
      use Control
      use Fields_mod
      use Operator_mod
      use iso_fortran_env, only: output_unit, error_unit
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT) :: GR(Ndim, Ndim, N_FL, N_grc)
      integer, intent(IN)    :: N_op, i_wlk
      real(Kind=kind(0.d0)), intent(OUT)   :: Hs_new

      ! Local ::
      type(Fields)   ::  nsigma_new
      complex(Kind=kind(0.d0)) :: Ratio(N_FL), Ratio_f(N_FL), ratiotot, Z1, exp_o_new, exp_o_old
      integer ::  n, ns, m, nf, i, Op_dim, op_dim_nf, nu_spin, nu_c, n_prop, i_grc
      complex(Kind=kind(0.d0)) :: Z, D_Mat, myexp, s1, s2, ratioD
      real(Kind=kind(0.d0)) :: S0_ratio, pi = acos(-1.d0)

      real(Kind=kind(0.d0)) :: Weight, st_r, ed_r, sum_ratio, rand_nu
      complex(Kind=kind(0.d0)) :: alpha, beta, g_loc
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: Mat, Delta
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: u, v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: y_v, xp_v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: x_v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: Zarr, grarr
      complex(Kind=kind(0.d0)), dimension(:), allocatable :: sxv, syu

      real(Kind=kind(0.d0)), dimension(:), allocatable :: ratio_field, field_list
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: ratio_o

      call nsigma_new%make(1, 1)

      op_dim = op_v(n_op, 1)%N_non_zero
      do nf = 2, N_FL
         if (op_dim < op_V(n_op, nf)%N_non_zero) op_dim = op_v(n_op, nf)%N_non_zero
      end do

      if (op_dim > 0) then
         allocate (Mat(Op_dim, Op_Dim), Delta(Op_dim, N_FL), u(Ndim, Op_dim), v(Ndim, Op_dim))
         allocate (y_v(Ndim, Op_dim), xp_v(Ndim, Op_dim), x_v(Ndim, Op_dim))
      end if

      ! Compute the ratio for each posibility
      nf = 1
      nsigma_new%t(1) = Op_V(n_op, nf)%type
      if (Op_V(n_op, nf)%type .eq. 1) then
         allocate (field_list(2))
         allocate (ratio_field(2))
         allocate (ratio_o(2, N_slat))
         field_list(1) = 1.d0; 
         field_list(2) = 2.d0; 
         nu_spin = 2
      elseif (Op_V(n_op, nf)%type .eq. 2) then
         allocate (field_list(4))
         allocate (ratio_field(4))
         allocate (ratio_o(4, N_slat))
         field_list(1) = 1.d0; 
         field_list(2) = -1.d0; 
         field_list(3) = 2.d0; 
         field_list(4) = -2.d0; 
         nu_spin = 4
      end if

      do nu_c = 1, nu_spin

         nsigma_new%f(1, 1) = field_list(nu_c) !real(ns_new,kind=kind(0.d0))
         S0_ratio = ham%S0(n_op, 1, dble(field_list(nu_c)))

         exp_o_new = 0.d0
         exp_o_old = 0.d0
         do ns = 1, N_slat
            i_grc = ns + (i_wlk - 1)*N_slat

            do nf = 1, N_FL
               g_loc = Op_V(n_op, nf)%g
               Z1 = g_loc*(nsigma_new%Phi(1, 1))
               op_dim_nf = Op_V(n_op, nf)%N_non_zero
               do m = 1, op_dim_nf
                  myexp = exp(Z1*Op_V(n_op, nf)%E(m))
                  Z = myexp - 1.d0
                  Delta(m, nf) = Z
                  do n = 1, op_dim_nf
                     Mat(n, m) = -Z*GR(Op_V(n_op, nf)%P(n), Op_V(n_op, nf)%P(m), nf, i_grc)
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
            end do

            ratiotot = product(Ratio)**dble(N_SUN)
            ratio_o(nu_c, ns) = ratiotot

            exp_o_new = exp_o_new + exp(overlap(i_grc))*ratio_o(nu_c, ns)
            exp_o_old = exp_o_old + exp(overlap(i_grc))

         end do

         ratiotot = exp_o_new/exp_o_old
         weight = S0_ratio*dble(ratiotot*nsigma_new%px0(1, 1))
         if (weight .le. 0.d0) weight = 0.d0
         ratio_field(nu_c) = weight

      end do

      sum_ratio = sum(ratio_field)

      if (sum_ratio .le. 0.d0) then
         weight_k(i_wlk) = cmplx(0.d0,pi,kind(0.d0))
         Hs_new = field_list(1) ! randomly set up a Hs_new for output
      else
            !! Decide the field in the next propagation
         n_prop = -1
         rand_nu = ranf_wrap()
         st_r = 0.d0
         do nu_c = 1, nu_spin
            ed_r = st_r + ratio_field(nu_c)/sum_ratio
            if ((rand_nu .ge. st_r) .and. (rand_nu .lt. ed_r)) then
               n_prop = nu_c
               exit
            end if
            st_r = st_r + ratio_field(nu_c)/sum_ratio
         end do
         nsigma_new%f(1, 1) = field_list(n_prop)

            !! update weight
         weight_k(i_wlk) = weight_k(i_wlk) + sum_ratio

            !! Update Green's function
         ! update delta
         do nf = 1, N_FL
            g_loc = Op_V(n_op, nf)%g
            Z1 = g_loc*(nsigma_new%Phi(1, 1))
            op_dim_nf = Op_V(n_op, nf)%N_non_zero
            do m = 1, op_dim_nf
               myexp = exp(Z1*Op_V(n_op, nf)%E(m))
               Z = myexp - 1.d0
               Delta(m, nf) = Z
            end do
         end do

         do ns = 1, N_slat

            i_grc = ns + (i_wlk - 1)*N_slat
               !! update overlap
            overlap(i_grc) = overlap(i_grc) + log(ratio_o(n_prop, ns))

            do nf = 1, N_FL
               ! Setup u(i,n), v(n,i)
               op_dim_nf = Op_V(n_op, nf)%N_non_zero
               if (op_dim_nf > 0) then
                  beta = 0.d0
                  call zlaset('N', Ndim, op_dim_nf, beta, beta, u, size(u, 1))
                  call zlaset('N', Ndim, op_dim_nf, beta, beta, v, size(v, 1))
                  do n = 1, op_dim_nf
                     u(Op_V(n_op, nf)%P(n), n) = Delta(n, nf)
                     do i = 1, Ndim
                        v(i, n) = -GR(Op_V(n_op, nf)%P(n), i, nf, i_grc)
                     end do
                     v(Op_V(n_op, nf)%P(n), n) = 1.d0 - GR(Op_V(n_op, nf)%P(n), Op_V(n_op, nf)%P(n), nf, i_grc)
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
                     call ZCOPY(Ndim, gr(1, Op_V(n_op, nf)%P(1), nf, i_grc), 1, xp_v(1, 1), 1)
                     Z = -x_v(Op_V(n_op, nf)%P(1), 1)
                     call ZGERU(Ndim, Ndim, Z, xp_v(1, 1), 1, y_v(1, 1), 1, gr(1, 1, nf, i_grc), Ndim)
                  else
                     allocate (zarr(op_dim_nf, op_dim_nf), grarr(NDim, op_dim_nf))
                     Zarr = x_v(Op_V(n_op, nf)%P(1:op_dim_nf), :)
                     grarr = gr(:, Op_V(n_op, nf)%P(1:op_dim_nf), nf, i_grc)
                     beta = 0.d0
                     alpha = 1.d0
                     call ZGEMM('N', 'N', NDim, op_dim_nf, op_dim_nf, alpha, grarr, Ndim, Zarr, op_dim_nf, beta, xp_v, Ndim)
                     deallocate (Zarr, grarr)
                     beta = cmplx(1.0d0, 0.0d0, kind(0.d0))
                     alpha = -1.d0
                     call ZGEMM('N', 'T', Ndim, Ndim, op_dim_nf, alpha, xp_v, Ndim, y_v, Ndim, beta, gr(1, 1, nf, i_grc), Ndim)
                  end if
               end if
            end do

         end do

         ! Output the field
         Hs_new = nsigma_new%f(1, 1) ! real(ns_new,Kind=kind(0.d0))

      end if

      if (op_dim > 0) then
         deallocate (Mat, Delta, u, v)
         deallocate (y_v, xp_v, x_v)
         deallocate (ratio_field, field_list, ratio_o)
      end if

      call nsigma_new%clear()

   end subroutine upgrade

   subroutine upgrade_mc(gr, n_op, nt, hs_new, i_wlk, overlap_mc)

      use Hamiltonian_main
      use Random_wrap
      use Control
      use Fields_mod
      use Operator_mod
      use iso_fortran_env, only: output_unit, error_unit
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT) :: gr(ndim, ndim, n_fl, n_grc)
      integer, intent(IN)    :: n_op, i_wlk, nt
      complex(Kind=kind(0.d0)), intent(INOUT) :: overlap_mc(n_grc)
      real(Kind=kind(0.d0)), intent(IN)    :: hs_new

      ! Local ::
      type(Fields)   ::  nsigma_new
      complex(Kind=kind(0.d0)) :: Ratio(N_FL), Ratio_f(N_FL), Ratiotot, Z1, Phase_a_array(N_FL)
      integer ::  n, m, nf, i, Op_dim, op_dim_nf, nu_spin, nu_c, n_prop, ratioD, ns, i_grc
      complex(Kind=kind(0.d0)) :: Z, D_Mat, myexp, s1, s2
      logical                   :: toggle

      real(Kind=kind(0.d0)) :: Weight, st_r, ed_r, sum_ratio, rand_nu
      real(Kind=kind(0.d0)) :: S0_ratio
      complex(Kind=kind(0.d0)) :: alpha, beta, g_loc
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: Mat, Delta
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: u, v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: y_v, xp_v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: x_v
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: Zarr, grarr
      complex(Kind=kind(0.d0)), dimension(:), allocatable :: sxv, syu

      complex(kind=kind(0.d0)), dimension(:), allocatable :: ratio_o
      complex(kind=kind(0.d0)) :: exp_o_new, exp_o_old
      real(Kind=kind(0.d0)) :: dexp_o_new, dexp_o_old

      allocate (ratio_o(N_slat))
      call nsigma_new%make(1, 1)

      op_dim = op_v(n_op, 1)%N_non_zero
      do nf = 2, N_FL
         if (op_dim < Op_V(n_op, nf)%N_non_zero) op_dim = op_v(n_op, nf)%N_non_zero
      end do

      if (op_dim > 0) then
         allocate (Mat(Op_dim, Op_Dim), Delta(Op_dim, N_FL), u(Ndim, Op_dim), v(Ndim, Op_dim))
         allocate (y_v(Ndim, Op_dim), xp_v(Ndim, Op_dim), x_v(Ndim, Op_dim))
      end if

      S0_ratio = ham%S0(n_op, nt, Hs_new)

      ! Compute the ratio for each posibility
      nf = 1
      nsigma_new%f(1, 1) = Hs_new
      nsigma_new%t(1) = Op_V(n_op, nf)%type
      exp_o_new = 0.d0
      exp_o_old = 0.d0

      do ns = 1, N_slat
         i_grc = ns + (i_wlk - 1)*N_slat

         do nf = 1, N_FL
            g_loc = op_v(n_op, nf)%g
            Z1 = g_loc*(nsigma_new%Phi(1, 1) - nsigma_bp(i_wlk)%Phi(n_op, nt))
            op_dim_nf = Op_V(n_op, nf)%N_non_zero
            do m = 1, op_dim_nf
               myexp = exp(Z1*Op_V(n_op, nf)%E(m))
               Z = myexp - 1.d0
               Delta(m, nf) = Z
               do n = 1, op_dim_nf
                  Mat(n, m) = -Z*GR(Op_V(n_op, nf)%P(n), Op_V(n_op, nf)%P(m), nf, i_grc)
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
         end do

         ratiotot = product(Ratio)**dble(N_SUN)*nsigma_new%Gama(1, 1)/nsigma_bp(i_wlk)%Gama(n_op, nt)
         ratio_o(ns) = ratiotot

         exp_o_new = exp_o_new + exp(overlap_mc(i_grc))*ratio_o(ns)
         exp_o_old = exp_o_old + exp(overlap_mc(i_grc))

      end do

      dexp_o_new = dble(exp_o_new)
      dexp_o_old = dble(exp_o_old)

      weight = S0_ratio*abs(dexp_o_new/dexp_o_old)

      toggle = .false.
      if (weight > ranf_wrap()) then
         toggle = .true.

         do ns = 1, N_slat

            i_grc = ns + (i_wlk - 1)*N_slat
              !! update overlap
            overlap_mc(i_grc) = overlap_mc(i_grc) + log(ratio_o(ns))

            do nf = 1, N_FL
               op_dim_nf = Op_V(n_op, nf)%N_non_zero
               if (op_dim_nf > 0) then
                  beta = 0.d0
                  call zlaset('N', Ndim, op_dim_nf, beta, beta, u, size(u, 1))
                  call zlaset('N', Ndim, op_dim_nf, beta, beta, v, size(v, 1))
                  do n = 1, op_dim_nf
                     u(Op_V(n_op, nf)%P(n), n) = Delta(n, nf)
                     do i = 1, Ndim
                        v(i, n) = -gr(Op_V(n_op, nf)%P(n), i, nf, i_grc)
                     end do
                     v(Op_V(n_op, nf)%P(n), n) = 1.d0 - gr(Op_V(n_op, nf)%P(n), Op_V(n_op, nf)%P(n), nf, i_grc)
                  end do

                  call zlaset('N', Ndim, op_dim_nf, beta, beta, x_v, size(x_v, 1))
                  call zlaset('N', Ndim, op_dim_nf, beta, beta, y_v, size(y_v, 1))
                  i = op_V(n_op, nf)%P(1)
                  x_v(i, 1) = u(i, 1)/(1.d0 + v(i, 1)*u(i, 1))
                  call zcopy(ndim, v(:, 1), 1, y_v(:, 1), 1)
                  do n = 2, op_dim_nf
                     call zcopy(ndim, u(:, n), 1, x_v(:, n), 1)
                     call zcopy(ndim, v(:, n), 1, y_v(:, n), 1)
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
                     call ZCOPY(Ndim, gr(1, Op_V(n_op, nf)%P(1), nf, i_grc), 1, xp_v(1, 1), 1)
                     z = -x_v(Op_V(n_op, nf)%P(1), 1)
                     call ZGERU(Ndim, Ndim, Z, xp_v(1, 1), 1, y_v(1, 1), 1, gr(1, 1, nf, i_grc), Ndim)
                  else
                     allocate (zarr(op_dim_nf, op_dim_nf), grarr(NDim, op_dim_nf))
                     Zarr = x_v(Op_V(n_op, nf)%P(1:op_dim_nf), :)
                     grarr = gr(:, Op_V(n_op, nf)%P(1:op_dim_nf), nf, i_grc)
                     beta = 0.d0
                     alpha = 1.d0
                     call zgemm('N', 'N', NDim, op_dim_nf, op_dim_nf, alpha, grarr, Ndim, Zarr, op_dim_nf, beta, xp_v, Ndim)
                     deallocate (Zarr, grarr)
                     beta = cmplx(1.0d0, 0.0d0, kind(0.d0))
                     alpha = -1.d0
                     call ZGEMM('N', 'T', Ndim, Ndim, op_dim_nf, alpha, xp_v, Ndim, y_v, Ndim, beta, gr(1, 1, nf, i_grc), Ndim)
                  end if
               end if
            end do

         end do

         ! Flip the spin
         nsigma_bp(i_wlk)%f(n_op, nt) = nsigma_new%f(1, 1)
      end if

      if (op_dim > 0) then
         deallocate (Mat, Delta, u, v)
         deallocate (y_v, xp_v, x_v)
      end if

      deallocate (ratio_o)

      call control_upgrade(toggle)

      call nsigma_new%clear()

   end subroutine Upgrade_mc

end module upgrade_mod
