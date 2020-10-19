
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Prints out the bins.  No need to change this routine.
!-------------------------------------------------------------------
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau

          !Local
          Integer :: I


          Do I = 1,Size(Obs_scal,1)
             Call Print_bin_Vec(Obs_scal(I), Group_Comm)
          enddo
          Do I = 1,Size(Obs_eq,1)
             Call Print_bin_Latt(Obs_eq(I), Group_Comm)
          enddo
          If (Ltau == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call Print_bin_Latt(Obs_tau(I), Group_Comm)
             enddo
          endif

        end Subroutine Pr_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Initializes observables to zero before each bins.  No need to change
!> this routine.
!-------------------------------------------------------------------
        Subroutine  Init_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau

          ! Local
          Integer :: I

          Do I = 1,Size(Obs_scal,1)
             Call Obser_vec_Init(Obs_scal(I))
          Enddo

          Do I = 1,Size(Obs_eq,1)
             Call Obser_Latt_Init(Obs_eq(I))
          Enddo

          If (Ltau == 1) then
             Do I = 1,Size(Obs_tau,1)
                Call Obser_Latt_Init(Obs_tau(I))
             Enddo
          Endif

        end Subroutine Init_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specify a global move on a given time slice tau.
!>
!> @details
!> @param[in] ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!> @param[out] T0_Proposal_ratio, Real
!> \verbatim
!>  T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)
!> \endverbatim
!> @param[out] S0_ratio, Real
!> \verbatim
!>  S0_ratio = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
!> \endverbatim
!> @param[out] Flip_length  Integer
!> \verbatim
!>  Number of flips stored in the first  Flip_length entries of the array Flip_values.
!>  Has to be smaller than NDIM
!> \endverbatim
!> @param[out] Flip_list  Integer(Ndim)
!> \verbatim
!>  List of spins to be flipped: nsigma%f(Flip_list(1),ntau) ... nsigma%f(Flip_list(Flip_Length),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!> @param[out] Flip_value  Real(Ndim)
!> \verbatim
!>  Flip_value(:)= nsigma%flip(Flip_list(:),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!--------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)


          Implicit none
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: T0_Proposal_ratio,  S0_ratio
          Integer                , INTENT(OUT) :: Flip_list(:)
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: Flip_value(:)
          Integer, INTENT(OUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          ! Local
          Integer :: n_op, n, ns
          Real (Kind=Kind(0.d0)) :: T0_proposal

          Flip_length = nranf(4)
          do n = 1,flip_length
             n_op = nranf(size(OP_V,1))
             Flip_list(n)  = n_op
             Flip_value(n) = nsigma%flip(n_op,ntau)
             If ( OP_V(n_op,1)%type == 1 ) then
                S0_ratio          =   S0(n_op,ntau,Flip_value(n))
                T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0_ratio) ! No move prob
                If ( T0_Proposal > Ranf_wrap() ) then
                   T0_Proposal_ratio =  1.d0 / S0_ratio
                else
                   T0_Proposal_ratio = 0.d0
                endif
             else
                T0_Proposal_ratio = 1.d0
                S0_ratio          = 1.d0
             endif
          Enddo

        end Subroutine Global_move_tau

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> The user can set the initial field.
!>
!> @details
!> @param[OUT] Initial_field Real(:,:)
!> \verbatim
!>  Upon entry Initial_field is not allocated. If alloacted then it will contain the
!>  the initial field
!> \endverbatim
!--------------------------------------------------------------------
     Subroutine  Hamiltonian_set_nsigma(Initial_field)
        Implicit none

        Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(OUT) :: Initial_field


      end Subroutine Hamiltonian_set_nsigma

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This routine allows to user to  determine the global_tau sampling parameters at run time
!> It is especially usefull if these parameters are dependent on other parameters.
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Subroutine Overide_global_tau_sampling_parameters(Nt_sequential_start,Nt_sequential_end,N_Global_tau)

        Implicit none
        Integer, Intent(INOUT) :: Nt_sequential_start,Nt_sequential_end, N_Global_tau
      end Subroutine Overide_global_tau_sampling_parameters

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Single spin flip S0 ratio
!> @details
!> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
!> a spin flip of Operator n on time slice nt
!> @details
!--------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)
        Implicit none
        !> Operator index
        Integer, Intent(IN) :: n
        !> Time slice
        Integer, Intent(IN) :: nt
        !> New local field on time slice nt and operator index n
        Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new

        Integer :: nt1,I
        !Write(6,*) "Hi1"

        S0 = 1.d0

      end function S0
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Global moves
!>
!> @details
!>  This routine generates a
!>  global update  and returns the propability T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!> @param [OUT]  T0_Proposal_ratio Real
!> \verbatimam
!>  T0_Proposal_ratio  =  T0( sigma_new -> sigma_old ) /  T0( sigma_old -> sigma_new)
!> \endverbatim
!> @param [OUT]  Size_clust Real
!> \verbatim
!>  Size of cluster that will be flipped.
!> \endverbatim
!-------------------------------------------------------------------
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)

          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Type (Fields),  Intent(IN)  :: nsigma_old

          ! Local
          Integer :: N_op, N_tau, n1,n2, n

          T0_Proposal_ratio  = 1.d0
          size_clust         = 3.d0

          nsigma%f = nsigma_old%f
          nsigma%t = nsigma_old%t
          N_op  = size(nsigma_old%f,1)
          N_tau = size(nsigma_old%f,2)
          Do n  = 1, Nint(size_clust)
             n1 = nranf(N_op)
             n2 = nranf(N_tau)
             nsigma%f(n1,n2) = nsigma_old%flip(n1,n2)
          enddo

        End Subroutine Global_move
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes the ratio exp(S0(new))/exp(S0(old))
!>
!> @details
!> This function computes the ratio \verbatim  e^{-S0(nsigma)}/e^{-S0(nsigma_old)} \endverbatim
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!-------------------------------------------------------------------
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none

          ! Arguments
          Type (Fields),  INTENT(IN) :: nsigma_old

          ! Local
          Integer :: I,n,n1,n2,n3,n4,nt,nt1, nc_F, nc_J, nc_h_p, nc_h_m


          Delta_S0_global = 1.d0


        end Function Delta_S0_global
