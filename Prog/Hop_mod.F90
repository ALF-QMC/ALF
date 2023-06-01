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
!>
!> @brief
!> This module computes and stores the exponential of the hopping matrix.  It also provide routines to carry out multiplications  with
!> \f$  e^{ - \Delta \tau  H_t  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n) }   \f$,
!> \f$  e^{   \Delta \tau  H_t  }   = \prod_{n=N}^{1} e^{   \Delta \tau_n  H_t(n) }   \f$,
!> \f$  e^{ - \Delta \tau  H_t/2  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n)/2 }   \f$, and
!> \f$  e^{   \Delta \tau  H_t/2  }   = \prod_{n=N}^{1} e^{   \Delta \tau_n  H_t(n)/2 }   \f$.
!> The last equality are important for the symmetric Trotter option. (See variable: Symm in the Hamiltonian module)
!
!>
!--------------------------------------------------------------------


    Module Hop_mod

      Use Hamiltonian_main
      Use Operator_mod
      Use Random_wrap
      Use DynamicMatrixArray_mod
      Use ContainerElementBase_mod
      Use OpTTypes_mod
      use iso_fortran_env, only: output_unit, error_unit

      ! Private variables
      Type(DynamicMatrixArray), private, allocatable :: ExpOpT_vec(:) ! for now we have for simplicity for each flavour a vector
      Integer, private, save ::  Ncheck
      Real (Kind=Kind(0.d0)), private, save  :: Zero

      Contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function serves as a central entry point to collect the 
!> processing that occurs in mapping an OpT input matrix to the internal
!> matrix-like data structure.
!
!> @param ExpOpT_vec[inout] a DynamicMatrixArray structure to which we append new elements.
!> @param op[in] an Operator that describes an OpT hopping matrix.
!
!--------------------------------------------------------------------
        subroutine OpT_postprocess(ExpOpT_vec, op)
            use Operator_mod
            implicit none
            
            Type(DynamicMatrixArray), intent(inout) :: ExpOpT_vec
            Type(Operator), intent(in) :: op
            
            Class(CmplxExpOpT), pointer :: cmplxexp => null()
            Class(RealExpOpT), pointer :: realexp => null()
            
            if (Op_is_real(op)) then
                ! branch for real operators
                    allocate(realexp) ! Yep, this is a manifest memory leak. Using the ptr we can allocate onto the same variable
                    call realexp%init(op)
                    call ExpOpT_vec%pushback(realexp)
                else
                ! branch for complex operators
                    allocate(cmplxexp)
                    call cmplxexp%init(op)
                    call ExpOpT_vec%pushback(cmplxexp)
                endif
        end subroutine
        
      
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This functions sets up the exponentiated matrices.
!> We symmetrize the upper part of those matrices.
!
!--------------------------------------------------------------------
        subroutine Hop_mod_init
          use runtime_error_mod
          Implicit none

          Integer :: nc, nf

          Ncheck = size(Op_T,1)
          If ( size(Op_T,2) /= N_FL ) then
             Write(error_unit,*) 'Hop_mod_init: Error in the number of flavors.'
             CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          Endif

          allocate(ExpOpT_vec(N_FL))

          do nf = 1,N_FL
             call ExpOpT_vec(nf)%init()
             do nc = 1,Ncheck
                call OpT_postprocess(ExpOpT_vec(nf), Op_T(nc, nf))
             enddo
          enddo

          Zero = 1.E-12
        end subroutine Hop_mod_init

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthr(In,nf)


          ! InOut:  In = e^{ -dtau T }.IN
          Implicit none

          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer, intent(IN) :: nf

          !Local
          Integer :: nc
          class(ContainerElementBase), pointer :: dummy

          do nc =  Ncheck,1,-1
            dummy => ExpOpT_vec(nf)%at(nc)
            call dummy%lmult(In)
          Enddo
        end Subroutine Hop_mod_mmthr

        Subroutine Hop_mod_mmthr_m1(In,nf)


          ! InOut:  In = e^{  dtau T }.IN
          Implicit none

          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf

          !Local
          Integer :: nc
          class(ContainerElementBase), pointer :: dummy

          do nc =  1,Ncheck
            dummy => ExpOpT_vec(nf)%at(nc)
            call dummy%lmultinv(In)
          Enddo

        end Subroutine Hop_mod_mmthr_m1

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthl (In,nf)


          ! InOut:  In = IN * e^{ -dtau T }
          Implicit none

          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf

          !Local
          Integer :: nc
          class(ContainerElementBase), pointer :: dummy

          do nc =  1, Ncheck
            dummy => ExpOpT_vec(nf)%at(nc)
            call dummy%rmult(In)
          Enddo

        end Subroutine Hop_mod_mmthl

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthlc (In,nf)


          ! InOut:  In = IN * e^{ -dtau T }
          Implicit none

          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf

          !Local
          Integer :: nc
          class(ContainerElementBase), pointer :: dummy

          do nc =  1, Ncheck
            dummy => ExpOpT_vec(nf)%at(nc)
            call dummy%lmult(In)
          Enddo

        end Subroutine Hop_mod_mmthlc

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthl_m1 (In, nf)


          ! InOut:  In = IN * e^{ dtau T }
          Implicit none

          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf

          !Local
          Integer :: nc
          class(ContainerElementBase), pointer :: dummy

          do nc =  Ncheck,1,-1
            dummy => ExpOpT_vec(nf)%at(nc)
            call dummy%rmultinv(In)
          Enddo

        end Subroutine Hop_mod_mmthl_m1

!!$        Subroutine  Hop_mod_test
!!$
!!$          Implicit none
!!$
!!$          Complex (Kind=Kind(0.d0)) ::  IN(Ndim,Ndim),Out(Ndim,Ndim)
!!$          Complex (Kind=Kind(0.d0)) ::  Test(Ndim,Ndim)
!!$
!!$          Integer :: I,J
!!$
!!$          DO I = 1,Ndim
!!$             DO J = 1,Ndim
!!$                IN(J,I) = cmplx(Ranf(),Ranf())
!!$             ENDDO
!!$          ENDDO
!!$
!!$          !Write(6,*) IN
!!$        end Subroutine Hop_mod_test

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Out = exp(-Delta tau /2 T) In exp( Delta tau /2 T)
!>
!
!--------------------------------------------------------------------
        Subroutine Hop_mod_Symm(Out,In)

          Implicit none
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Intent(Out):: Out
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Intent(IN):: In

          Integer :: nf, nc, nf_eff
          class(ContainerElementBase), pointer :: dummy

          Out = In
          Do nf_eff = 1, N_FL_eff !size(In,3)
             nf=Calc_Fl_map(nf_eff)
             do nc =  Ncheck,1,-1
                dummy => ExpOpT_vec(nf)%at(nc)
                call dummy%adjointaction(Out(:, :, nf))
             enddo
          enddo

        End Subroutine Hop_mod_Symm
      end Module Hop_mod
