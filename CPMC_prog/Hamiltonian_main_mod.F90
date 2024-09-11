!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module defines the interface between the Hamiltonians (= model and observables definition) and
!> the Monte Carlo core. Hamiltonians are defined as submodules of this module. The Monte Carlo core
!> has only access to public members of this module. For defining a new Hamiltonian named <new_ham_name>, 
!> the user has to create the file Hamiltonians/Hamiltonian_<new_ham_name>_smod.F90 and add the line
!> <new_ham_name> to Hamiltonians.list.

!> @details
!> The public variables of this module are the following
!>
!>
!> @param [public] ham
!> \verbatim
!> class(ham_base), allocatable
!> Object that contains all the Hamiltonian-specific procedures needed by the Monte Carlo core.
!> The procedures can be overloaded in the Hamiltonians. \endverbatim
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
      Use runtime_error_mod
      Use Operator_mod, only: Operator
      Use WaveFunction_mod, only: WaveFunction
      Use Observables
      Use Fields_mod, only: Fields
      use iso_fortran_env, only: output_unit, error_unit
      
    Implicit none
    
    private
    public :: Alloc_Ham, ham_base, ham
#ifdef __PGI
    public :: Obs_scal
#endif
      type ham_base
      contains
        procedure, nopass :: ham_set => ham_set_base
        procedure, nopass :: Alloc_obs => Alloc_obs_base
        procedure, nopass :: Obser => Obser_base
        procedure, nopass :: ObserT => ObserT_base
        procedure, nopass :: E0_local => E0_local_base
        procedure, nopass :: localmeas => localmeas_base
        procedure, nopass :: sum_weight => sum_weight_base
        procedure, nopass :: update_fac_norm => update_fac_norm_base
        procedure, nopass :: Pr_obs => Pr_obs_base
        procedure, nopass :: Init_obs => Init_obs_base
        procedure, nopass :: S0 => S0_base
        procedure, nopass :: weight_reconstruction => weight_reconstruction_base
        procedure, nopass :: GR_reconstruction => GR_reconstruction_base
        procedure, nopass :: GRT_reconstruction => GRT_reconstruction_base
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5 => write_parameters_hdf5_base
#endif
      end type ham_base
      
      class(ham_base), allocatable :: ham

      Type (Operator),     dimension(:,:), allocatable, public :: Op_V
      Type (Operator),     dimension(:,:), allocatable, public :: Op_T
      Type (WaveFunction), dimension(:,:), allocatable, public :: WF_L
      Type (WaveFunction), dimension(:,:), allocatable, public :: WF_R
      Logical            , dimension(:),   allocatable, public :: Calc_Fl
      Integer            , dimension(:),   allocatable, public :: Calc_Fl_map
      Type (Fields)      , dimension(:),   allocatable, public :: nsigma_bp, nsigma_qr
      Integer      , public        :: Ndim
      Integer      , public        :: N_FL, N_FL_eff
      Integer      , public        :: N_SUN
      Integer      , public        :: N_wlk
      Integer      , public        :: N_slat
      Integer      , public        :: N_grc
      Integer      , public        :: N_wlk_mpi
      Integer      , public        :: N_grc_mpi
      Integer      , public        :: ltrot
      Integer      , public        :: Group_Comm
      Logical      , public        :: Symm
      Logical      , public        :: reconstruction_needed
      
      Real    (Kind=Kind(0.d0)), public :: fac_norm
      Real    (Kind=Kind(0.d0)), dimension(:), allocatable, public :: weight_k
      Complex (Kind=Kind(0.d0)), dimension(:), allocatable, public :: overlap

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

    subroutine Alloc_Ham(ham_name)
       Implicit none
       Character (len=64), intent(in) :: ham_name

       Select Case (ham_name)
#include "Hamiltonians_case.h"
!!$       This file will be dynamically generated and appended
       Case default
          write(error_unit, '("A","A","A")') 'Hamiltonian ', ham_name, ' not yet implemented!'
          CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
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
      CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
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
               CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
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
    !> Computes equal time observables
    !> @details
    !> @param [IN] Gr   Complex(:,:,:)
    !> \verbatim
    !>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
    !> \endverbatim
    !> @param [IN] Ntau Integer
    !> \verbatim
    !>  Time slice
    !> \endverbatim
    !-------------------------------------------------------------------
          subroutine Obser_base(GR,GR_mix,i_wlk,sum_w)

             Implicit none

             Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
             Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR_mix(Ndim,Ndim,N_FL)
             Complex (Kind=Kind(0.d0)), Intent(IN) :: sum_w
             Integer                  , Intent(IN) :: i_wlk
             Logical, save              :: first_call=.True.
             
             If  (first_call)    then 
                write(error_unit, *) "Warning: Obser not implemented."
                first_call=.false.
             endif
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
    !-------------------------------------------------------------------
          Subroutine ObserT_base(NT, GT0, G0T, G00, GTT, i_wlk,sum_w)
             Implicit none
    
             Integer         , INTENT(IN) :: NT, i_wlk
             Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL), G0T(Ndim,Ndim,N_FL)
             Complex (Kind=Kind(0.d0)), INTENT(IN) :: G00(Ndim,Ndim,N_FL), GTT(Ndim,Ndim,N_FL)
             Complex (Kind=Kind(0.d0)), INTENT(IN) :: sum_w
             Logical, save              :: first_call=.True.
             
             If  (first_call)    then 
                write(error_unit, *) "Warning: ObserT not implemented."
                first_call=.false.
             endif
    
          end Subroutine ObserT_base

          complex (Kind=Kind(0.d0)) function E0_local_base(GR)
            Implicit none
             
            Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
            
            E0_local_base = cmplx(0.d0, 0.d0, kind(0.d0))
            
          end function E0_local_base
          
          complex (Kind=Kind(0.d0)) function localmeas_base(GR)
            Implicit none
             
            Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
            
            localmeas_base = cmplx(0.d0, 0.d0, kind(0.d0))
            
          end function localmeas_base

          complex (Kind=Kind(0.d0)) function sum_weight_base
            Implicit none
            
            sum_weight_base = cmplx(0.d0, 0.d0, kind(0.d0))
            
          end function sum_weight_base

          Subroutine update_fac_norm_base(GR, ntw)
            Implicit none
             
            complex (Kind=Kind(0.d0)), intent(in) :: GR(Ndim,Ndim,N_FL,N_wlk)
            integer                  , intent(in) :: ntw
            
          end subroutine update_fac_norm_base

    !--------------------------------------------------------------------
    !> @author
    !> ALF Collaboration
    !>
    !> @brief
    !> Prints out the bins.  No need to change this routine.
    !-------------------------------------------------------------------
          Subroutine  Pr_obs_base(LTAU)
    
             Implicit none
             
             Integer, Intent(In) :: Ltau
    
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
!> @brief
!> Reconstructs dependent flavors of the configuration's weight.
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!--------------------------------------------------------------------
          subroutine weight_reconstruction_base(weight)
            implicit none
            complex (Kind=Kind(0.d0)), Intent(inout) :: weight(:)

            write(error_unit, *) 'weight_reconstruction not defined!'
            CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
            
          end subroutine weight_reconstruction_base


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of equal time Greens function
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!-------------------------------------------------------------------
          subroutine GR_reconstruction_base(GR)
            
            Implicit none
            
            Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim,N_FL)
            
            write(error_unit, *) "Warning: GR_reconstruction not implemented."
            CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          end Subroutine GR_reconstruction_base

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of time displaced Greens function G0T and GT0
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] GT0, G0T,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!-------------------------------------------------------------------
         Subroutine GRT_reconstruction_base(GT0, G0T)
           Implicit none
           
           Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GT0(Ndim,Ndim,N_FL), G0T(Ndim,Ndim,N_FL)
           
           write(error_unit, *) "Warning: GRT_reconstruction not implemented."
           CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
         end Subroutine GRT_reconstruction_base
         
#ifdef HDF5
         subroutine write_parameters_hdf5_base(filename)
           implicit none
           
           Character (len=64), intent(in) :: filename

           write(output_unit,*)
           write(output_unit,*) "ATTENTION:     Base implementation of write_parameters_hdf5 is getting calling!"
           write(output_unit,*) "This routine does not actually store the parameters in the hdf5 file."
           write(output_unit,*) "Usually, a python script is generating a model specific routine, granted that the Hamiltonian"
           write(output_unit,*) "respects the design rule. Without parameters, the HDF5 file's compatibility with pyALF is limited."
           write(output_unit,*)
           
         end subroutine write_parameters_hdf5_base
#endif

    end Module Hamiltonian_main
