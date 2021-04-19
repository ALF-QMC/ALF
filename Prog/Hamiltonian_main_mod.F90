!  Copyright (C) 2016 - 2020 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module defines the  Hamiltonian and observables.  Here, we have included a
!> set of predefined Hamiltonians. They include the Hubbard and SU(N) tV models
!> on honeycomb, pi-flux and square lattices.

!> @details
!> The public variables of this module are the following
!>
!>
!> @param [public] OP_V
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> List of operators of type=1,2 and 3 describing the sequence of interactions on a time slice.
!> The first index runs over this sequence. The second corresponds to the flavor index.  \endverbatim
!>
!> @param [public] OP_T
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> Sequence of  operators  accounting for the  hopping on a  time slice. This can include  various
!> checkerboard decompositions. The first index runs over this sequence. The second corresponds to
!> the flavor index. \endverbatim
!> *  The progagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n}  \f$.  That is
!> first the hopping and then the potential energy.
!>
!>@param [public] WF_L
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma
!> \verbatim Type(Fields)
!> Contains all auxiliary fields in the variable f(:,:). The first index runs through the operator
!> sequence. The second through the time slices.   \endverbatim
!
!> @param [public]  Ndim
!> \verbatim Integer
!> Total number of orbitals. e.g. # unit cells * # orbitals per unit cell.  \endverbatim
!
!> @param [public]  N_FL
!> \verbatim Integer
!> # of flavors.  Propagation is block diagonal in flavors.  \endverbatim
!
!> @param [public]  N_SUN
!> \verbatim Integer
!> # of colors.  Propagation is color independent.  \endverbatim
!>
!> @param [public] Ltrot
!> \verbatim Integer
!> Available measurment interval in units of Delta Tau. \endverbatim
!>
!> @param [public] Thtrot
!>  \verbatim Integer
!> Effective projection parameter in units of Delta Tau.  (Only relevant if projective option is turned on) \endverbatim
!>
!> @param [public] Projector
!> \verbatim Logical
!> Flag for projector. If true then the total number of time slices will correspond to Ltrot + 2*Thtrot \endverbatim
!>
!> @param [public] Group_Comm
!> \verbatim Integer
!> Defines MPI communicator  \endverbatim
!
!> @param [public] Symm
!> \verbatim Logical  \endverbatim
!> If set to true then the green functions will be symmetrized
!> before being  sent to the Obser, ObserT subroutines.
!> In particular, the transformation,  \f$ \tilde{G} =  e^{-\Delta \tau T /2 } G e^{\Delta \tau T /2 } \f$
!> will be carried out  and \f$ \tilde{G} \f$  will be sent to the Obser and ObserT subroutines.  Note that
!> if you want to use this  feature, then you have to be sure the hopping and interaction terms are decomposed
!> symmetrically. If Symm is true, the propagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=N_T}^{1}e^{T_n/2} \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n/2}  \f$
!>
!>
!> You still have to add some docu for the other private variables in this module.
!>
!--------------------------------------------------------------------

    Module Hamiltonian_main
      Use Operator_mod, only: Operator
      Use WaveFunction_mod, only: WaveFunction
      Use Observables
      Use Fields_mod, only: Fields
      use iso_fortran_env, only: output_unit, error_unit
      
    Implicit none
    
    private
    public :: Alloc_Ham, ham_base, ham
#ifdef __PGI
    public :: Obs_scal, Obs_eq, Obs_tau
#endif
      type ham_base
      contains
        procedure, nopass :: ham_set => ham_set_base
        procedure, nopass :: Alloc_obs => Alloc_obs_base
        procedure, nopass :: Obser => Obser_base
        procedure, nopass :: ObserT => ObserT_base
        procedure, nopass :: Pr_obs => Pr_obs_base
        procedure, nopass :: Init_obs => Init_obs_base
        procedure, nopass :: Global_move_tau => Global_move_tau_base
        procedure, nopass :: Hamiltonian_set_nsigma => Hamiltonian_set_nsigma_base
        procedure, nopass :: Overide_global_tau_sampling_parameters => Overide_global_tau_sampling_parameters_base
        procedure, nopass :: Global_move => Global_move_base
        procedure, nopass :: Delta_S0_global => Delta_S0_global_base
        procedure, nopass :: S0 => S0_base
        procedure, nopass :: Ham_Langevin_HMC_S0 => Ham_Langevin_HMC_S0_base
      end type ham_base
      
      class(ham_base), allocatable :: ham

      Type (Operator),     dimension(:,:), allocatable, public :: Op_V
      Type (Operator),     dimension(:,:), allocatable, public :: Op_T
      Type (WaveFunction), dimension(:),   allocatable, public :: WF_L
      Type (WaveFunction), dimension(:),   allocatable, public :: WF_R
      Type (Fields), public        :: nsigma
      Integer      , public        :: Ndim
      Integer      , public        :: N_FL
      Integer      , public        :: N_SUN
      Integer      , public        :: Ltrot
      Integer      , public        :: Thtrot
      Logical      , public        :: Projector
      Integer      , public        :: Group_Comm
      Logical      , public        :: Symm


      !>    Privat Observables
      Type (Obser_Vec ), dimension(:), allocatable :: Obs_scal
      Type (Obser_Latt), dimension(:), allocatable :: Obs_eq
      Type (Obser_Latt), dimension(:), allocatable :: Obs_tau


#include "Hamiltonians_interface.h"
!!$      This file will  be dynamically generated and appended
!!$      interface
!!$         module subroutine Ham_Alloc_Kondo()
!!$         end subroutine Ham_Alloc_Kondo
!!$         module subroutine Ham_Alloc_Hubbard()
!!$         end subroutine Ham_Alloc_Hubbard
!!$         module subroutine Ham_Alloc_Hubbard_Plain_Vanilla()
!!$         end subroutine Ham_Alloc_Hubbard_Plain_Vanilla
!!$         module subroutine Ham_Alloc_tV()
!!$         end subroutine Ham_Alloc_tV
!!$         module subroutine Ham_Alloc_LRC()
!!$         end subroutine Ham_Alloc_LRC
!!$         module subroutine Ham_Alloc_Z2_Matter()
!!$         end subroutine Ham_Alloc_Z2_Matter
!!$      end interface
      
    contains

    subroutine Alloc_Ham()
       Implicit none
       Integer :: ierr
       Character (len=64) :: ham_name
       NAMELIST /VAR_HAM_NAME/ ham_name
       
       OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(error_unit,*) 'Alloc_Ham: unable to open <parameters>',ierr
          error stop 1
       END IF
       READ(5,NML=VAR_HAM_NAME)
       CLOSE(5)

       Select Case (ham_name)
#include "Hamiltonians_case.h"
!!$       This file will be dynamically generated and appended
!!$       Case ('Kondo')
!!$          call Ham_Alloc_Kondo()
!!$       Case ('Hubbard')
!!$          call Ham_Alloc_Hubbard()
!!$       Case ('Hubbard_Plain_Vanilla')
!!$          call Ham_Alloc_Hubbard_Plain_Vanilla()
!!$       Case ('tV')
!!$          call Ham_Alloc_tV()
!!$       Case ('LRC')
!!$          call Ham_Alloc_LRC()
!!$       Case ('Z2_Matter')
!!$          call Ham_Alloc_Z2_Matter()
       Case default
          write(error_unit, '("A","A","A")') 'Hamiltonian ', ham_name, ' not yet implemented!'
          error stop 1
       end Select
    end subroutine Alloc_Ham
    
    !--------------------------------------------------------------------
    !> @brief
    !> Sets the Hamiltonian.
    !> @details
    !> This has to be overloaded in the Hamiltonian submodule.
    !--------------------------------------------------------------------
    subroutine Ham_Set_base()
      implicit none
      write(error_unit, *) 'Ham_set not defined!'
      error stop 1
    end subroutine Ham_Set_base
    
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
          Real (Kind=Kind(0.d0)) function S0_base(n,nt,Hs_new)
             Implicit none
             !> Operator index
             Integer, Intent(IN) :: n
             !> Time slice
             Integer, Intent(IN) :: nt
             !> New local field on time slice nt and operator index n
             Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new

             S0_base = 1.d0
             If ( Op_V(n,1)%type == 1 ) then
               write(error_unit, *) 'function S0 not implemented'
               error stop 1
             endif

          end function S0_base


    !--------------------------------------------------------------------
    !> @author
    !> ALF Collaboration
    !>
    !> @brief
    !> Specifiy the equal time and time displaced observables
    !> @details
    !--------------------------------------------------------------------
          Subroutine  Alloc_obs_base(Ltau)

             Implicit none
             !>  Ltau=1 if time displaced correlations are considered.
             Integer, Intent(In) :: Ltau
             write(error_unit, *) "Warning: Alloc_obs not implemented."
          End Subroutine Alloc_obs_base


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
          Subroutine Global_move_base(T0_Proposal_ratio, nsigma_old, size_clust)

             Implicit none
             Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
             Type (Fields),  Intent(IN)  :: nsigma_old

             write(error_unit, *) 'Global_move not implemented'
             error stop 1

          End Subroutine Global_move_base


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
          Real (Kind=kind(0.d0)) Function Delta_S0_global_base(Nsigma_old)

             !  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
             Implicit none

             ! Arguments
             Type (Fields),  INTENT(IN) :: nsigma_old

             Delta_S0_global_base = 1.d0

          end Function Delta_S0_global_base


    !--------------------------------------------------------------------
    !> @author
    !> ALF Collaboration
    !>
    !> @brief
    !> Computes equal time observables
    !> @details
    !> @param [IN] Gr   Complex(:,:,:)
    !> \verbatim
    !>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
    !> \endverbatim
    !> @param [IN] Phase   Complex
    !> \verbatim
    !>  Phase
    !> \endverbatim
    !> @param [IN] Ntau Integer
    !> \verbatim
    !>  Time slice
    !> \endverbatim
    !-------------------------------------------------------------------
          subroutine Obser_base(GR,Phase,Ntau, Mc_step_weight)

             Implicit none

             Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
             Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
             Integer, INTENT(IN)          :: Ntau
             Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
             
             write(error_unit, *) "Warning: Obser not implemented."

          end Subroutine Obser_base


    !--------------------------------------------------------------------
    !> @author
    !> ALF Collaboration
    !>
    !> @brief
    !> Computes time displaced  observables
    !> @details
    !> @param [IN] NT, Integer
    !> \verbatim
    !>  Imaginary time
    !> \endverbatim
    !> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
    !> \verbatim
    !>  Green functions:
    !>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
    !>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
    !>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
    !>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
    !> \endverbatim
    !> @param [IN] Phase   Complex
    !> \verbatim
    !>  Phase
    !> \endverbatim
    !-------------------------------------------------------------------
          Subroutine ObserT_base(NT, GT0, G0T, G00, GTT, PHASE, Mc_step_weight)
             Implicit none
    
             Integer         , INTENT(IN) :: NT
             Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL), G0T(Ndim,Ndim,N_FL)
             Complex (Kind=Kind(0.d0)), INTENT(IN) :: G00(Ndim,Ndim,N_FL), GTT(Ndim,Ndim,N_FL)
             Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
             Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
             
             write(error_unit, *) "Warning: ObserT not implemented."
    
          end Subroutine ObserT_base


    !--------------------------------------------------------------------
    !> @author
    !> ALF Collaboration
    !>
    !> @brief
    !> Prints out the bins.  No need to change this routine.
    !-------------------------------------------------------------------
          Subroutine  Pr_obs_base(LTAU)
    
             Implicit none
    
             Integer,  Intent(In) ::  Ltau
    
             !Local
             Integer :: I
    
    
             if ( allocated(Obs_scal) ) then
               Do I = 1,Size(Obs_scal,1)
                  Call Print_bin_Vec(Obs_scal(I), Group_Comm)
               enddo
             endif
             if ( allocated(Obs_eq) ) then
               Do I = 1,Size(Obs_eq,1)
                  Call Print_bin_Latt(Obs_eq(I), Group_Comm)
               enddo
             endif
             if ( allocated(Obs_tau) ) then
               Do I = 1,Size(Obs_tau,1)
                  Call Print_bin_Latt(Obs_tau(I), Group_Comm)
               enddo
             endif
    
          end Subroutine Pr_obs_base


    !--------------------------------------------------------------------
    !> @author
    !> ALF Collaboration
    !>
    !> @brief
    !> Initializes observables to zero before each bins.  No need to change
    !> this routine.
    !-------------------------------------------------------------------
          Subroutine  Init_obs_base(Ltau)
    
             Implicit none
             Integer, Intent(In) :: Ltau
    
             ! Local
             Integer :: I
    
             if ( allocated(Obs_scal) ) then
               Do I = 1,Size(Obs_scal,1)
                  Call Obser_vec_Init(Obs_scal(I))
               Enddo
             endif
    
             if ( allocated(Obs_eq) ) then
               Do I = 1,Size(Obs_eq,1)
                  Call Obser_Latt_Init(Obs_eq(I))
               Enddo
             endif
    
             if ( allocated(Obs_tau) ) then
               Do I = 1,Size(Obs_tau,1)
                  Call Obser_Latt_Init(Obs_tau(I))
               Enddo
             Endif
    
          end Subroutine Init_obs_base


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
          Subroutine Global_move_tau_base(T0_Proposal_ratio, S0_ratio, &
                &                     Flip_list, Flip_length,Flip_value,ntau)

             Implicit none
             Real (Kind = Kind(0.d0)),INTENT(OUT) :: T0_Proposal_ratio,  S0_ratio
             Integer                , INTENT(OUT) :: Flip_list(:)
             Real (Kind = Kind(0.d0)),INTENT(OUT) :: Flip_value(:)
             Integer, INTENT(OUT) :: Flip_length
             Integer, INTENT(IN)  :: ntau
             
             write(error_unit, *) 'Global_move_tau not implemented'
             error stop 1
          end Subroutine Global_move_tau_base


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
          Subroutine  Hamiltonian_set_nsigma_base(Initial_field)
             Implicit none

             Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(INOUT) :: Initial_field

          end Subroutine Hamiltonian_set_nsigma_base


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
          Subroutine Overide_global_tau_sampling_parameters_base(Nt_sequential_start,Nt_sequential_end,N_Global_tau)

             Implicit none
             Integer, Intent(INOUT) :: Nt_sequential_start,Nt_sequential_end, N_Global_tau
          end Subroutine Overide_global_tau_sampling_parameters_base


  !--------------------------------------------------------------------
  !> @author 
  !> ALF Collaboration
  !>
  !> @brief 
  !>   Forces_0  = \partial S_0 / \partial s  are calculated and returned to  main program.
  !> 
  !-------------------------------------------------------------------
          Subroutine Ham_Langevin_HMC_S0_base(Forces_0)

            Implicit none

            Real (Kind=Kind(0.d0)), Intent(inout), allocatable :: Forces_0(:,:)

            Forces_0  = 0.d0
            
          end Subroutine Ham_Langevin_HMC_S0_base


    end Module Hamiltonian_main
