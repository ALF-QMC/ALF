!  Copyright (C) 2016, 2017 The ALF project
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
      
      ! Private variables
      Complex (Kind=Kind(0.d0)), allocatable, private :: Exp_T(:,:,:,:), Exp_T_M1(:,:,:,:)
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

        Subroutine Hop_mod_mmthr(In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = e^{ -dtau T }.IN
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=Kind(0.d0)), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer, intent(IN) :: nf
          
          !Local 
          Complex (Kind=Kind(0.D0)) :: alpha, beta
          Complex(Kind = Kind(0.D0)), allocatable, dimension(:,:) :: tmp
          Integer :: nc, n
          
          Out = In
          alpha = 1.D0
          beta = 0.D0
          Allocate(tmp(Ndim_hop, Ndim_hop))
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('L','U',Ndim_hop,Ndim,Ndim,Exp_T(:,:,nc,nf),Op_T(nc,nf)%P,Out)
!                 do n = 1,Ndim_hop
!                     call ZCOPY(Ndim, Out(Op_T(nc,nf)%P(n),1), NDim, V_Hlp(n, 1), Ndim_hop)
!                 enddo
!                 CALL ZHEMM('L', 'U', Ndim_hop, Ndim, alpha, Exp_T(:,:,nc,nf),Ndim_hop,V_hlp(1,1),NDim_hop,beta,V_HLP1(1,1),Ndim_hop)
!                 do n = 1,Ndim_hop
!                     call ZCOPY(Ndim, V_hlp1(n,1), Ndim_hop, OUT(OP_T(nc,nf)%P(n),1), Ndim)
!                 Enddo
             Endif
          Enddo
          deallocate(tmp)
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
          Complex (Kind=Kind(0.D0)) :: a, b
          a = 1.D0
          b = 0.D0
          
          Out = In
          do nc =  1,Ncheck
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('L','U',Ndim_hop,Ndim,Ndim,Exp_T_m1(:,:,nc,nf),Op_T(nc,nf)%P,Out)
!                 do I = 1,Ndim
!                    do n = 1,Ndim_hop
!                       V_Hlp(n,I) = Out(Op_T(nc,nf)%P(n),I)
!                    enddo
!                 enddo
!                 CALL ZHEMM('L','U', Ndim_hop, Ndim, a, Exp_T_m1(:, :,nc,nf),Ndim_hop, V_hlp(1,1), NDim_hop, b, V_HLP1(1,1),Ndim_hop)
!                 DO I = 1,Ndim
!                    do n = 1,Ndim_hop
!                       OUT(OP_T(nc,nf)%P(n),I) = V_hlp1(n,I)
!                    Enddo
!                 Enddo
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
          Integer :: nc, n
          Complex (Kind=Kind(0.D0)) :: alpha, beta
          
          alpha = 1.D0
          beta = 0.D0
          Out = In
          do nc =  1, Ncheck
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('R','U',Ndim_hop,Ndim,Ndim,Exp_T(:,:,nc,nf),Op_T(nc,nf)%P,Out)
!                 do n = 1,Ndim_hop
!                   call zcopy(Ndim, Out(1, Op_T(nc,nf)%P(n)), 1, U_Hlp(1, n), 1)
!                 enddo
!                 CALL ZHEMM('R', 'U', Ndim, Ndim_hop, alpha, Exp_T(:, :, nc, nf), Ndim_hop, U_hlp(1,1), NDim, beta, U_HLP1(1,1),Ndim)
!                 do n = 1,Ndim_hop
!                   call zcopy(Ndim, U_hlp1(1, n), 1, OUT(1,OP_T(nc,nf)%P(n)), 1)
!                 Enddo
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthl

!--------------------------------------------------------------------

        Subroutine Hop_mod_mmthl_m1 (In, Out, nf)
          

          ! In:   IN
          ! Out:  OUT = IN * e^{ dtau T }
          Implicit none
          
          Complex (Kind=Kind(0.d0)), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=Kind(0.d0)), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer :: nf
          
          !Local 
          Integer :: nc, n
          Complex (Kind=Kind(0.D0)) :: a, b
          a = 1.D0
          b = 0.D0
          
          Out = In
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                call ZSLHEMM('R','U',Ndim_hop,Ndim,Ndim,Exp_T_m1(:,:,nc,nf),Op_T(nc,nf)%P,Out)
!                 do n = 1,Ndim_hop
!                    call zcopy(Ndim, Out(1, Op_T(nc,nf)%P(n)), 1, U_Hlp(1, n), 1)
!                 enddo
!                 CALL ZHEMM('R','U', Ndim, Ndim_hop, a, Exp_T_m1(:, :,nc,nf),Ndim_hop, U_hlp(1,1), NDim, b, U_HLP1(1,1),Ndim)
!                 do n = 1,Ndim_hop
!                    call zcopy(Ndim, U_Hlp1(1, n), 1, Out(1, Op_T(nc,nf)%P(n)), 1)
!                 Enddo
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
