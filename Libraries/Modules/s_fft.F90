!--------------------------------------------------------------------
!> @author 
!> Zi Hong Liu
!
!> @brief 
!> Do the discrete fast fourier transformation with MKL
!
!--------------------------------------------------------------------

!==============================================================================!
  subroutine onedimension_fft( ndim, X )
!==============================================================================!
! 1D complex to complex
  Use MKL_DFTI
  implicit none
  integer, intent(in) :: ndim
  Complex(Kind=Kind(0.d0)), intent(inout) :: X(ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
  Integer :: Status
  !...put input data into X(1),...,X(32); Y(1),...,Y(32)

  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
     DFTI_COMPLEX, 1, ndim )
  Status = DftiCommitDescriptor( My_Desc1_Handle )
  Status = DftiComputeForward( My_Desc1_Handle, X )
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by {X(1),X(2),...,X(32)}

  end subroutine onedimension_fft
!==============================================================================!

!==============================================================================!
  subroutine onedimension_invfft( ndim, X )
!==============================================================================!
! 1D complex to complex
  Use MKL_DFTI
  implicit none
  integer, intent(in) :: ndim
  Complex(Kind=Kind(0.d0)), intent(inout) :: X(ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
  Integer :: Status
  !...put input data into X(1),...,X(32); Y(1),...,Y(32)

  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
     DFTI_COMPLEX, 1, ndim )
  Status = DftiCommitDescriptor( My_Desc1_Handle )
  Status = DftiComputeBackward( My_Desc1_Handle, X )
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by {X(1),X(2),...,X(32)}

  end subroutine onedimension_invfft
!================================================================================!

!================================================================================!
  subroutine Discrete_Forward_FFT( L1, L2, Ndim, X )
!================================================================================!
! 2D complex to complex discrete fourier transform
  Use MKL_DFTI
  implicit none
  integer, intent(in) :: L1, L2, Ndim
  Complex (Kind=Kind(0.d0)), intent(inout) :: X(Ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
  Integer :: Status, LN(2)
  !...put input data into X_2D(j,k), Y_2D(j,k), 1<=j=32,1<=k<=100
  !...set LN(1) = ndim, LN(2) = ndim
  !...the transform is a ndim-by-ndim
  
  LN(1) = L1; LN(2) = L2
   
  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
            DFTI_COMPLEX, 2, LN)
  Status = DftiCommitDescriptor( My_Desc1_Handle)
  Status = DftiComputeForward( My_Desc1_Handle, X)
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by X_2D(j,k), 1<=j<=32, 1<=k<=100
  
  end subroutine Discrete_Forward_FFT
!==================================================================================!

!================================================================================!
  subroutine Discrete_Backward_FFT( L1, L2, Ndim, X )
!================================================================================!
! 2D complex to complex discrete fourier transform
  Use MKL_DFTI
  implicit none
  integer, intent(in) :: L1, L2, Ndim
  Complex (Kind=Kind(0.d0)), intent(inout) :: X(Ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
  Integer :: Status, LN(2)
  !...put input data into X_2D(j,k), Y_2D(j,k), 1<=j=32,1<=k<=100
  !...set LN(1) = ndim, LN(2) = ndim
  !...the transform is a ndim-by-ndim
  
  LN(1) = L1; LN(2) = L2
   
  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
            DFTI_COMPLEX, 2, LN)
  Status = DftiCommitDescriptor( My_Desc1_Handle)
  Status = DftiComputeBackward( My_Desc1_Handle, X)
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by X_2D(j,k), 1<=j<=32, 1<=k<=100
  
  end subroutine Discrete_Backward_FFT
!==================================================================================!

!================================================================================!
  subroutine Discrete_Forward_FFT_3D( lx, ly, lz, Ndim, X )
!================================================================================!
! 3D complex to complex discrete fourier transform
  Use MKL_DFTI
  implicit none
  integer, intent(in) :: lx, ly, lz, Ndim
  Complex (Kind=Kind(0.d0)), intent(inout) :: X(Ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
  Integer :: Status, LN(3)
  
  LN(1) = lx; LN(2) = ly; LN(3) = lz
   
  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
            DFTI_COMPLEX, 3, LN)
  Status = DftiCommitDescriptor( My_Desc1_Handle)
  Status = DftiComputeForward( My_Desc1_Handle, X)
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by X
  
  end subroutine Discrete_Forward_FFT_3D
!==================================================================================!

!================================================================================!
  subroutine Discrete_Backward_FFT_3D( lx, ly, lz, Ndim, X )
!================================================================================!
! 3D complex to complex discrete backward fourier transform
  Use MKL_DFTI
  implicit none
  integer, intent(in) :: lx, ly, lz, Ndim
  Complex (Kind=Kind(0.d0)), intent(inout) :: X(Ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
  Integer :: Status, LN(2)
  
  LN(1) = lx; LN(2) = ly; LN(3) = lz
   
  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
            DFTI_COMPLEX, 3, LN)
  Status = DftiCommitDescriptor( My_Desc1_Handle)
  Status = DftiComputeBackward( My_Desc1_Handle, X)
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by X
  
  end subroutine Discrete_Backward_FFT_3D
!==================================================================================!
