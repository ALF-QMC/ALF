!  Copyright (C) 2016 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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


    Module Hop_mod

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This module computes and store the exponential of the hopping matrix. 
!> It also provide routines to carry out multiplication  with e^{ -dtau T }       
!
!--------------------------------------------------------------------

      Use Hamiltonian
      Use Random_wrap
      Use MyMats 
      
      ! Private variables
      Complex (Kind=Kind(0.d0)), allocatable, private :: Exp_T(:,:,:,:), Exp_T_M1(:,:,:,:)
      Complex (Kind=Kind(0.d0)), allocatable, private :: U_HLP(:,:), U_HLP1(:,:),  V_HLP(:,:), V_HLP1(:,:)
      Integer, private, save ::  Ncheck, Ndim_hop
      Real (Kind=Kind(0.d0)), private, save  :: Zero

      Contains
        
        subroutine Hop_mod_init

          Implicit none
          
          Integer :: nc, nf
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

          Allocate ( Exp_T   (Ndim_hop,Ndim_hop,Ncheck,N_FL) ) 
          Allocate ( Exp_T_M1(Ndim_hop,Ndim_hop,Ncheck,N_FL) ) 
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
             enddo
          enddo
          
          Zero = 1.E-12
          
        end subroutine Hop_mod_init

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthr(In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = e^{ -dtau T }.IN
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=Kind(0.d0)), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer, intent(IN) :: nf
          
          !Local 
          Complex (Kind=Kind(0.D0)) :: alpha, beta
          Integer :: nc, n
          
          Out = In
          alpha = 1.D0
          beta = 0.D0
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                do n = 1,Ndim_hop
                    call ZCOPY(Ndim, Out(Op_T(nc,nf)%P(n),1), NDim, U_Hlp(1, n), 1)
                enddo
                Call ZGEMM('N', 'T', Ndim, Ndim_hop, Ndim_hop, alpha, U_Hlp, NDim, Exp_T(:,:,nc,nf), Ndim_hop, beta, U_HLP1, Ndim)
                do n = 1,Ndim_hop
                    call ZCOPY(Ndim, U_hlp1(1,n), 1, OUT(OP_T(nc,nf)%P(n),1), Ndim)
                Enddo
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthr

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthr_m1(In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = e^{  dtau T }.IN
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=Kind(0.d0)), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer :: nf
          
          !Local 
          Integer :: nc, I, n 
          
          
          Out = In
          do nc =  1,Ncheck
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                do I = 1,Ndim
                   do n = 1,Ndim_hop
                      V_Hlp(n,I) = Out(Op_T(nc,nf)%P(n),I)
                   enddo
                enddo
                Call mmult(V_HLP1,Exp_T_m1(:,:,nc,nf),V_Hlp)
                DO I = 1,Ndim
                   do n = 1,Ndim_hop
                      OUT(OP_T(nc,nf)%P(n),I) = V_hlp1(n,I)
                   Enddo
                Enddo
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthr_m1

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthl (In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = IN * e^{ -dtau T }
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=Kind(0.d0)), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer :: nf
          
          !Local 
          Integer :: nc, I, n  
          
          Out = In
          do nc =  1, Ncheck
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                do n = 1,Ndim_hop
                  call zcopy(Ndim, Out(1, Op_T(nc,nf)%P(n)), 1, U_Hlp(1, n), 1)
                enddo
                Call mmult(U_Hlp1,U_Hlp,Exp_T(:,:,nc,nf))
                do n = 1,Ndim_hop
                  call zcopy(Ndim, U_hlp1(1, n), 1, OUT(1,OP_T(nc,nf)%P(n)), 1)
                Enddo
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthl

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthl_m1 (In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = IN * e^{ dtau T }
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=Kind(0.d0)), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer :: nf
          
          !Local 
          Integer :: nc, I, n  
          
          Out = In
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                do n = 1,Ndim_hop
                   call zcopy(Ndim, Out(1, Op_T(nc,nf)%P(n)), 1, U_Hlp(1, n), 1)
                enddo
                Call mmult(U_Hlp1,U_Hlp,Exp_T_M1(:,:,nc,nf))
                do n = 1,Ndim_hop
                   call zcopy(Ndim, U_Hlp1(1, n), 1, Out(1, Op_T(nc,nf)%P(n)), 1)
                Enddo
             Endif
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
       
      end Module Hop_mod
