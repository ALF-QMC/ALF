!  Copyright (C) 2016 - 2019 The ALF project
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

      Use Hamiltonian
      Use Random_wrap
      
      ! Private variables
      Complex (Kind=Kind(0.d0)), allocatable, private :: Exp_T(:,:,:,:), Exp_T_M1(:,:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: Exp_T_1D2(:,:,:,:), Exp_T_M1_1D2(:,:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: U_HLP(:,:), U_HLP1(:,:),  V_HLP(:,:), V_HLP1(:,:)
      Integer, private, save ::  Ncheck, Ndim_hop
      Real (Kind=Kind(0.d0)), private, save  :: Zero

      Contains
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

          Implicit none
          
          Integer :: nc, nf, i,j
          Complex (Kind=Kind(0.d0)) :: g

          Ncheck = size(Op_T,1)
          If ( size(Op_T,2) /= N_FL ) then 
             Write(6,*) 'Error in the number of flavors.'
             Stop
          Endif
          Ndim_hop = Op_T(1,1)%N
          !Write(6,*) 'In Hop_mod: ', Ndim, Ndim_hop, Ncheck
          Do nc = 1, Ncheck
             do nf = 1,N_FL
                if ( Ndim_hop /= Op_T(nc,nf)%N ) Then 
                   Write(6,*)  'Different size of Hoppings not  implemented '
                   Stop
                endif
             enddo
          enddo

          Allocate ( Exp_T       (Ndim_hop,Ndim_hop,Ncheck,N_FL) ) 
          Allocate ( Exp_T_M1    (Ndim_hop,Ndim_hop,Ncheck,N_FL) ) 
          Allocate ( Exp_T_1D2   (Ndim_hop,Ndim_hop,Ncheck,N_FL) ) 
          Allocate ( Exp_T_M1_1D2(Ndim_hop,Ndim_hop,Ncheck,N_FL) ) 

          Allocate ( V_Hlp(Ndim_hop,Ndim) )
          Allocate ( V_Hlp1(Ndim_hop,Ndim) )
          Allocate ( U_Hlp (Ndim, Ndim_hop) )
          Allocate ( U_Hlp1(Ndim, Ndim_hop) )
          
          Exp_T = cmplx(0.d0, 0.d0, kind(0.D0))
          Exp_T_M1 = cmplx(0.d0, 0.d0, kind(0.D0))
          do nf = 1,N_FL
             do nc = 1,Ncheck
                g = Op_T(nc,nf)%g
                Call  Op_exp(g,Op_T(nc,nf),Exp_T(:,:,nc,nf))
                g = -Op_T(nc,nf)%g
                Call  Op_exp(g,Op_T(nc,nf),Exp_T_M1(:,:,nc,nf))
                g = Op_T(nc,nf)%g/2.d0
                Call  Op_exp(g,Op_T(nc,nf),Exp_T_1D2(:,:,nc,nf))
                g = -Op_T(nc,nf)%g/2.d0
                Call  Op_exp(g,Op_T(nc,nf),Exp_T_M1_1D2(:,:,nc,nf))
                ! symmetrize the upper part of Exp_T and Exp_T_M1
                DO i = 1, Ndim_hop
                   DO j = i, Ndim_hop
                      Exp_T(i, j, nc, nf) = (Exp_T(i, j, nc, nf) + Conjg(Exp_T(j, i, nc, nf)))/2.D0
                      Exp_T_M1(i, j, nc, nf) = (Exp_T_M1(i, j, nc, nf) + Conjg(Exp_T_M1(j, i, nc, nf)))/2.D0
                   ENDDO
                ENDDO
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
          Integer :: nc, N1, N2
          
          N1=size(In,1)
          N2=size(In,2)
          
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('L','U',Ndim_hop,N1,N2,Exp_T(:,:,nc,nf),Op_T(nc,nf)%P,In)
             Endif
          Enddo
        end Subroutine Hop_mod_mmthr

!--------------------------------------------------------------------
        Subroutine Hop_mod_mmthr_1D2(In,nf)
          

          ! InOut:  In = e^{ -dtau T /2 }.IN
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer, intent(IN) :: nf
          
          !Local 
          Integer :: nc, N1, N2
          
          N1=size(In,1)
          N2=size(In,2)
          
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('L','U',Ndim_hop,N1,N2,Exp_T_1D2(:,:,nc,nf),Op_T(nc,nf)%P,In)
             Endif
          Enddo
        end Subroutine Hop_mod_mmthr_1D2

!--------------------------------------------------------------------
        
        Subroutine Hop_mod_mmthr_m1(In,nf)
          

          ! InOut:  In = e^{  dtau T }.IN
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf
          
          !Local 
          Integer :: nc , N1, N2
          
          N1=size(In,1)
          N2=size(In,2)
          
          do nc =  1,Ncheck
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('L','U',Ndim_hop,N1,N2,Exp_T_m1(:,:,nc,nf),Op_T(nc,nf)%P,In)
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthr_m1

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthl (In,nf)
          

          ! InOut:  In = IN * e^{ -dtau T }
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf
          
          !Local 
          Integer :: nc, N1, N2
          
          N1=size(In,1)
          N2=size(In,2)
          
          do nc =  1, Ncheck
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('R','U',Ndim_hop,N1,N2,Exp_T(:,:,nc,nf),Op_T(nc,nf)%P,In)
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthl

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthlc (In,nf)
          

          ! InOut:  In = IN * e^{ -dtau T }
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf
          
          !Local 
          Integer :: nc, N1, N2
          
          N1=size(In,1)
          N2=size(In,2)
          
          do nc =  1, Ncheck
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('L','U',Ndim_hop,N1,N2,Exp_T(:,:,nc,nf),Op_T(nc,nf)%P,In)
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthlc

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthl_m1 (In, nf)
          

          ! InOut:  In = IN * e^{ dtau T }
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf
          
          !Local 
          Integer :: nc, N1, N2
          
          N1=size(In,1)
          N2=size(In,2)
          
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('R','U',Ndim_hop,N1,N2,Exp_T_m1(:,:,nc,nf),Op_T(nc,nf)%P,In)
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthl_m1


!--------------------------------------------------------------------
        
        Subroutine Hop_mod_mmthl_m1_1D2(In, nf)
          

          ! InOut:  In = IN * e^{ dtau T/2 }
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(INOUT)  :: IN(:,:)
          Integer :: nf
          
          !Local 
          Integer :: nc, N1, N2
          
          N1=size(In,1)
          N2=size(In,2)
          
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('R','U',Ndim_hop,N1,N2,Exp_T_m1_1D2(:,:,nc,nf),Op_T(nc,nf)%P,In)
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthl_m1_1D2


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

          Integer :: nf

          Out = In
          Do nf = 1, size(In,3)
             Call Hop_mod_mmthr_1D2   (Out(:,:,nf), nf )
             Call Hop_mod_mmthl_m1_1D2(Out(:,:,nf), nf )
          enddo
          
        End Subroutine Hop_mod_Symm
      end Module Hop_mod
