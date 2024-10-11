!  Copyright (C) 2016 - 2022 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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


module upgrade_mod
   use runtime_error_mod
   implicit none
   contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This routine updates the field associated to the operator N_op on time
!> slice NT to the value ns_new.
!> @details
!> @param [inout] Gr(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  Green function.  In the move is accepted, then the Green function will be updated
!> \endverbatim
!> @param [in] N_op
!> Integer
!> \verbatim
!>  Operator number
!> \endverbatim
!> @param [in] Nt
!> Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!> @param [in] Op_dim
!> Integer
!> \verbatim
!>  Number of non-zero eigenvalues of the Operator OP_V(N_op,Nt)
!> \endverbatim
!> @param [in] HS_new
!> Real
!> \verbatim
!>  HS_new contains the new field for nsigma%f(n_op,nt)
!> \endverbatim
!> @param [in] Mode
!>  Character (Len=64)
!> \verbatim
!>  If  Mode=Intermediate  then the routine cumulates the fermion weight ratio in the variable
!>                         Prev_Ratiotot and updates the green function.
!>  If  Mode=Final         then the routine carries out the move and updates the Green function if accepted
!> \endverbatim
!> @param [inout] Prev_Ratiotot
!>  Complex
!> \verbatim
!>  Prev_Ratiotot = Prev_Ratiotot*Ratiotot.  Ratiotot corresponds to the ratio of the fermion weights.
!> \endverbatim
!> @param [in] S0_ratio,  T0_proposal_ratio
!>  Real
!> \verbatim
!>  If the mode is set to final then Weight = S0_Ratio*T0_proposal_ratio*Abs(Real(Ratiotot))
!> \endverbatim
!>  *  S0_ratio corresponds to ratio of the bosonic part of the action \f$ \frac{e^{-S_0(C')}}{e^{-S_0(C)}} \f$
!>  *  T0_proposal_ratio corresponds to the proposal probability ratio
!>   \f$ \frac{T_0 ( C'\rightarrow C)}{T_0 (C\rightarrow C')} \f$
!> @param [inout] toggle
!>  Logical
!> \verbatim
!>  Returns true if the move is accepted
!> \endverbatim
!> @param [inout] Phase
!>  Complex
!> \verbatim
!>  If mode=final and the move is accepted then phase is updated.
!> \endverbatim
!>  *  The phase is defeined as  \f$ e^{i \phi'}  = \frac{W(C')}{| W(C')|} \f$  where \f$ W(C') \f$ is the full fermion weight
!>  *  The weight reads \f$ W(C)  =   \left[ \left( \prod_{n,\tau}  \exp \left[ g(n,\tau) \alpha(n,\tau) \phi(\sigma(n,\tau)) \right] \right) \det(M(C))\right]^{N_{SUN}} \prod_{n,\tau }\gamma(\sigma(n,\tau)) \f$
!--------------------------------------------------------------------

      Subroutine Upgrade2(GR,N_op,NT,PHASE,Hs_new, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, toggle,  mode)
!--------------------------------------------------------------------

        Use Hamiltonian_main
        Use Random_wrap
        Use Control
        Use Fields_mod
        use iso_fortran_env, only: output_unit, error_unit
        Implicit none

        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim, N_FL)
        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Prev_Ratiotot
        Integer                  , INTENT(IN)    :: N_op, Nt
        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Phase
        Complex (Kind=Kind(0.d0)), INTENT(IN)    :: Hs_new
        Real    (Kind=Kind(0.d0)), INTENT(IN)    :: S0_ratio, T0_proposal_ratio
        Character (Len=64)       , INTENT(IN)    :: Mode
        Logical                  , INTENT(INOUT) :: toggle



        ! Local ::
        Type   (Fields)   ::  nsigma_new
        Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Ratiotot, Z1
        Integer ::  n,m,nf, nf_eff, i, Op_dim, op_dim_nf
        Complex (Kind=Kind(0.d0)) :: Z, D_Mat, myexp, s1, s2

        Real    (Kind=Kind(0.d0)) :: Weight, tmp_r
        Complex (Kind=Kind(0.d0)) :: alpha, beta, g_loc
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: Mat, Delta
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: u, v
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: y_v, xp_v
        Complex (Kind=Kind(0.d0)), Dimension(:, :), Allocatable :: x_v
        Complex (Kind=Kind(0.D0)), Dimension(:, :), Allocatable :: Zarr, grarr
        Complex (Kind=Kind(0.D0)), Dimension(:), Allocatable :: sxv, syu

        toggle = .false.
        ! if ( abs(OP_V(n_op,1)%g) < 1.D-12 )   return

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

        ! Compute the ratio
        nf = 1
        nsigma_new%f(1,1)  = Hs_new
        nsigma_new%t(1)    = Op_V(n_op,nf)%Type
        Do nf_eff = 1,N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           !Z1 = Op_V(n_op,nf)%g * ( Phi_st(ns_new,Op_V(n_op,nf)%type) -  Phi_st(ns_old,Op_V(n_op,nf)%type))
           g_loc = Op_V(n_op,nf)%g
           if (op_v(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)
           Z1 = g_loc * ( nsigma_new%Phi(1,1) -  nsigma%Phi(n_op,nt) )
           op_dim_nf = Op_V(n_op,nf)%N_non_zero
           Do m = 1,op_dim_nf
              myexp = exp( Z1* Op_V(n_op,nf)%E(m) )
              Z = myexp - 1.d0
              Delta(m,nf_eff) = Z
              do n = 1,op_dim_nf
                 Mat(n,m) = - Z * GR( Op_V(n_op,nf)%P(n), Op_V(n_op,nf)%P(m),nf )
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
              !  D_mat =  - Mat(2,1)*Mat(1,2)
           else
              D_mat = Det(Mat,op_dim_nf)
           endif
           Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
        Enddo

        !call reconstruct weight subroutine to fill the non-calculated blocks
        if (reconstruction_needed) call ham%weight_reconstruction(Ratio)

        Ratiotot = Product(Ratio)
        nf = 1
        !Ratiotot = (Ratiotot**dble(N_SUN)) * Gama_st(ns_new, Op_V(n_op,nf)%type)/Gama_st(ns_old, Op_V(n_op,nf)%type)
        Ratiotot = (Ratiotot**dble(N_SUN)) * nsigma_new%Gama(1,1)/nsigma%Gama(n_op,nt)

        !Write(6,*) Ratiotot

        If      ( str_to_upper(mode) == "FINAL"       ) Then
           !Write(6,*) "Up_I: ", Ratiotot
           Ratiotot = Ratiotot * Prev_Ratiotot
           !Write(6,*) "Up_f: ", Ratiotot
           Weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )
           !Write(6,*) Phase, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, ns_old,ns_new
        elseif  ( str_to_upper(mode) == "INTERMEDIATE" ) Then
           Weight = 1.5D0
           !Write(6,*) "Up_I: ", Ratiotot
           Prev_Ratiotot = Prev_Ratiotot*Ratiotot
        else
           Write(error_unit,*) 'Error in upgrade2'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        toggle = .false.
        if ( Weight > ranf_wrap() )  Then
           toggle = .true.
           If ( str_to_upper(mode) == "FINAL"  )  Phase = Phase * Ratiotot/sqrt(Ratiotot*conjg(Ratiotot))
           !Write(6,*) 'Accepted : ', Ratiotot

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
                        v(i,n) = - GR( Op_V(n_op,nf)%P(n), i, nf )
                    enddo
                    v(Op_V(n_op,nf)%P(n), n)  = 1.d0 - GR( Op_V(n_op,nf)%P(n),  Op_V(n_op,nf)%P(n), nf)
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
                    CALL ZCOPY(Ndim, gr(1, Op_V(n_op,nf)%P(1), nf), 1, xp_v(1, 1), 1)
                    Z = -x_v(Op_V(n_op,nf)%P(1), 1)
                    CALL ZGERU(Ndim, Ndim, Z, xp_v(1,1), 1, y_v(1, 1), 1, gr(1,1,nf), Ndim)
                ELSE
                    Allocate (zarr(op_dim_nf, op_dim_nf), grarr(NDim, op_dim_nf))
                    Zarr = x_v(Op_V(n_op,nf)%P(1:op_dim_nf), :)
                    grarr = gr(:, Op_V(n_op,nf)%P(1:op_dim_nf), nf)
                    beta  = 0.d0
                    alpha = 1.D0
                    CALL ZGEMM('N', 'N', NDim, op_dim_nf, op_dim_nf, alpha, grarr, Ndim, Zarr, op_dim_nf, beta, xp_v, Ndim)
                    Deallocate(Zarr, grarr)
                    beta  = cmplx ( 1.0d0, 0.0d0, kind(0.D0))
                    alpha = -1.D0
                    CALL ZGEMM('N','T',Ndim,Ndim,op_dim_nf,alpha,xp_v, Ndim, y_v, Ndim,beta,gr(1,1,nf), Ndim)
                ENDIF
              endif

           enddo

           ! Flip the spin
           nsigma%f(n_op,nt) = nsigma_new%f(1,1) ! real(ns_new,Kind=kind(0.d0))
        endif
        
        if (op_dim > 0) then
           deallocate ( Mat, Delta, u, v )
           deallocate ( y_v, xp_v, x_v )
        endif

        If ( str_to_upper(mode) == "FINAL" )  then
           Call Control_upgrade(toggle)
           Call Control_upgrade_eff(toggle)
        endif

        Call nsigma_new%clear()

      End Subroutine Upgrade2

end module upgrade_mod
