!  This is for  the Kondo project with tarun.
    Module Hop_mod


      Use Hamiltonian
      Use Random_wrap
      Use MyMats 
      
      ! Private variables
      Complex (Kind=8), allocatable, private :: Exp_T(:,:,:,:), Exp_T_M1(:,:,:,:)
      Complex (Kind=8), allocatable, private :: U_HLP(:,:), U_HLP1(:,:),  V_HLP(:,:), V_HLP1(:,:)
      Integer, private, save ::  Ncheck, Ndim_hop
      Real (Kind=8), private, save  :: Zero

      Contains
        
        subroutine Hop_mod_init

          Implicit none
          
          Integer :: nc, nf
          Complex (Kind=8) :: g

          Ncheck = size(Op_T,1)
          If ( size(Op_T,2) /= N_FL ) then 
             Write(6,*) 'Error in the number of flavors.'
             Stop
          Endif
          Ndim_hop = Op_T(1,1)%N
          Write(6,*) 'In Hop_mod: ', Ndim, Ndim_hop, Ncheck
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

!============================================================================
        Subroutine Hop_mod_mmthr(In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = e^{ -dtau T }.IN
          Implicit none
          
          Complex (Kind=8), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=8), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer :: nf
          
          !Local 
          Integer :: nc, n
          
          Out = In

          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                   do n = 1,Ndim_hop
                      V_Hlp(n,:) = Out(Op_T(nc,nf)%P(n),:)
                   enddo
                Call mmult(V_HLP1,Exp_T(:,:,nc,nf),V_Hlp)
                do n = 1,Ndim_hop
                    OUT(OP_T(nc,nf)%P(n),:) = V_hlp1(n,:)
                Enddo
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthr

!============================================================================
        Subroutine Hop_mod_mmthr_m1(In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = e^{  dtau T }.IN
          Implicit none
          
          Complex (Kind=8), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=8), intent(INOUT) :: OUT(Ndim,Ndim)
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

!============================================================================
        Subroutine Hop_mod_mmthl (In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = IN * e^{ -dtau T }
          Implicit none
          
          Complex (Kind=8), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=8), intent(INOUT) :: OUT(Ndim,Ndim)
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
!============================================================================
        Subroutine Hop_mod_mmthl_m1 (In, Out,nf)
          

          ! In:   IN
          ! Out:  OUT = IN * e^{ dtau T }
          Implicit none
          
          Complex (Kind=8), intent(IN)  :: IN(Ndim,Ndim)
          Complex (Kind=8), intent(INOUT) :: OUT(Ndim,Ndim)
          Integer :: nf
          
          !Local 
          Integer :: nc, I, n  
          
          Out = In
          do nc =  Ncheck,1,-1
             If ( dble( Op_T(nc,nf)%g*conjg(Op_T(nc,nf)%g) ) > Zero ) then
                do n = 1,Ndim_hop
                   do I = 1,Ndim
                      U_Hlp(I,n) = Out(I,Op_T(nc,nf)%P(n))
                   enddo
                enddo
                Call mmult(U_Hlp1,U_Hlp,Exp_T_M1(:,:,nc,nf))
                do n = 1,Ndim_hop
                   DO I = 1,Ndim
                      OUT(I,OP_T(nc,nf)%P(n)) = U_hlp1(I,n)
                   Enddo
                Enddo
             Endif
          Enddo
          
        end Subroutine Hop_mod_mmthl_m1
        
!============================================================================
!!$        Subroutine  Hop_mod_test
!!$          
!!$          Implicit none
!!$
!!$          Complex (Kind=8) ::  IN(Ndim,Ndim),Out(Ndim,Ndim)
!!$          Complex (Kind=8) ::  Test(Ndim,Ndim)
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
