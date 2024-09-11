!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Defines the wavefunction type
!> 
!
!--------------------------------------------------------------------

Module WaveFunction_mod

  Use MyMats

  Implicit none
  
  Type WaveFunction
     !> P is an Ndim x N_part matrix. N_part is the number of particles
     complex (Kind=Kind(0.d0)), allocatable :: P(:,:)
     !> Degeneracy of trial wave function 
     Real (Kind=Kind(0.d0)) :: Degen
  end type WaveFunction
  
Contains

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief
!> If the left and right trial wave functions are different, then one of the trial wave functions has to
!> be multiple by a phase such that  <Psi_L | Psi_R > = 1  otherwise, the program will produce a "fake"
!> negative sign problem.
!> 
!
!--------------------------------------------------------------------

  Subroutine WF_overlap(WF_L, WF_R, Z_norm)
    Implicit none
    Type (WaveFunction), intent(IN   )     :: WF_L 
    Type (WaveFunction), intent(INOUT)     :: WF_R
    Complex (Kind=Kind(0.d0)), intent(OUT) :: Z_norm
    
    
    ! Local
    Integer :: N_Part, Ndim, n,ne
    Complex (Kind=Kind(0.d0)), allocatable ::  mat(:,:)
    Complex (Kind=Kind(0.d0)) :: alpha, beta

    N_part = Size(WF_R%P,2)
    Ndim  = Size(WF_R%P,1)
    Allocate  (Mat(N_part,N_Part)) 

    alpha=1.d0
    beta=0.d0
    call ZGEMM('C','N',N_part,N_part,Ndim,alpha,WF_L%P(1,1),Ndim,WF_R%P(1,1),Ndim,beta,Mat(1,1),N_part)
    ! Mat = (WL_L%P)^{dagger} WL_L%R 
    
    Z_norm =  Det(Mat,N_part)

    Z_norm  = (cmplx(1.d0,0.d0,Kind(0.d0))/Z_norm)**(1.d0/Real(N_part,Kind(0.d0)))
    
    WF_R%P  = Z_norm * WF_R%P
    
    Deallocate  (Mat) 
    
  end subroutine WF_overlap
  
!--------------------------------------------------------------------

  Pure subroutine WF_alloc(WF, Ndim, N_part)
    Implicit none
    Type (WaveFunction), intent(INOUT) :: WF
    Integer, Intent(IN) :: Ndim, N_part
    Allocate (WF%P(Ndim, N_part))
    WF%P = cmplx(0.d0, 0.d0, kind(0.D0))
  end subroutine WF_alloc

!--------------------------------------------------------------------

  Pure subroutine WF_clear(WF)
    Implicit none
    Type (WaveFunction), intent(INOUT) :: WF
    Deallocate (WF%P)
  end subroutine WF_clear 

end Module WaveFunction_mod
