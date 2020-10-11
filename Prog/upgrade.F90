!  Copyright (C) 2016 - 2018 The ALF project
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

      Subroutine Upgrade2(GR,N_op,NT,PHASE,Op_dim,Hs_new, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, toggle,  mode)
!--------------------------------------------------------------------

        Use Hamiltonian
        Use Random_wrap
        Use Control
        Use Fields_mod
        use iso_fortran_env, only: output_unit, error_unit
        Implicit none

        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim, N_FL)
        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Prev_Ratiotot
        Integer                  , INTENT(IN)    :: N_op, Nt, Op_dim
        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Phase
        Real    (Kind=Kind(0.d0)), INTENT(IN)    :: Hs_new
        Real    (Kind=Kind(0.d0)), INTENT(IN)    :: S0_ratio, T0_proposal_ratio
        Character (Len=64)       , INTENT(IN)    :: Mode
        Logical                  , INTENT(INOUT) :: toggle



        ! Local ::
        Type   (Fields)   ::  nsigma_new
        Complex (Kind=Kind(0.d0)) :: Mat(Op_dim,Op_Dim), Delta(Op_dim,N_FL)
        Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Ratiotot, Z1
        Integer ::  n,m,nf, i
        Complex (Kind=Kind(0.d0)) :: Z, D_Mat, myexp, s1, s2

        Real    (Kind=Kind(0.d0)) :: Weight
        Complex (Kind=Kind(0.d0)) :: u(Ndim,Op_dim), v(Ndim,Op_dim) ,alpha, beta
        Complex (Kind=Kind(0.d0)) :: y_v(Ndim,Op_dim), xp_v(Ndim,Op_dim)
        Complex (Kind=Kind(0.d0)) :: x_v(Ndim,Op_dim)
        Complex (Kind=Kind(0.D0)), Dimension(:, :), Allocatable :: Zarr, grarr
        Complex (Kind=Kind(0.D0)), Dimension(:), Allocatable :: sxv, syu

        toggle = .false.
        ! if ( abs(OP_V(n_op,1)%g) < 1.D-12 )   return

        Call nsigma_new%make(1,1)

        ! Compute the ratio
        nf = 1
        nsigma_new%f(1,1)  = Hs_new !real(ns_new,kind=kind(0.d0))
        nsigma_new%t(1)    = Op_V(n_op,nf)%Type
        Do nf = 1,N_FL
           !Z1 = Op_V(n_op,nf)%g * ( Phi_st(ns_new,Op_V(n_op,nf)%type) -  Phi_st(ns_old,Op_V(n_op,nf)%type))
           Z1 = Op_V(n_op,nf)%g * ( nsigma_new%Phi(1,1) -  nsigma%Phi(n_op,nt) )
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
        !Ratiotot = (Ratiotot**dble(N_SUN)) * Gama_st(ns_new, Op_V(n_op,nf)%type)/Gama_st(ns_old, Op_V(n_op,nf)%type)
        Ratiotot = (Ratiotot**dble(N_SUN)) * nsigma_new%Gama(1,1)/nsigma%Gama(n_op,nt)

        !Write(6,*) Ratiotot

        If      (mode  == "Final"       ) Then
           !Write(6,*) "Up_I: ", Ratiotot
           Ratiotot = Ratiotot * Prev_Ratiotot
           !Write(6,*) "Up_f: ", Ratiotot
           Weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )
           !Write(6,*) Phase, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, ns_old,ns_new
        elseif  (mode == "Intermediate" ) Then
           Weight = 1.5D0
           !Write(6,*) "Up_I: ", Ratiotot
           Prev_Ratiotot = Prev_Ratiotot*Ratiotot
        else
           Write(error_unit,*) 'Error in upgrade2'
           error stop 1
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
                Allocate (zarr(Op_dim, Op_dim), grarr(NDim, Op_dim))
                Zarr = x_v(Op_V(n_op,nf)%P(1:Op_dim), :)
                grarr = gr(:, Op_V(n_op,nf)%P(1:Op_dim), nf)
                beta  = 0.d0
                alpha = 1.D0
                CALL ZGEMM('N', 'N', NDim, Op_Dim, Op_Dim, alpha, grarr, Ndim, Zarr, Op_dim, beta, xp_v, Ndim)
                Deallocate(Zarr, grarr)
                beta  = cmplx ( 1.0d0, 0.0d0, kind(0.D0))
                alpha = -1.D0
                CALL ZGEMM('N','T',Ndim,Ndim,Op_dim,alpha,xp_v, Ndim, y_v, Ndim,beta,gr(1,1,nf), Ndim)
              ENDIF

           enddo

           ! Flip the spin
           nsigma%f(n_op,nt) = nsigma_new%f(1,1) ! real(ns_new,Kind=kind(0.d0))
        endif

        If ( mode == "Final" )  then
           Call Control_upgrade(toggle)
           Call Control_upgrade_eff(toggle)
        endif

        Call nsigma_new%clear()

      End Subroutine Upgrade2
