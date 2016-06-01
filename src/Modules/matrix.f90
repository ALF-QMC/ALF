     MODULE Matrix



       Type Mat_C
          complex (Kind=8), pointer :: el(:,:)
          Integer :: dim
       end Type Mat_C

       Type Mat_R
          Real (Kind=8), pointer :: el(:,:)
          Integer :: dim
       end Type Mat_R

       Interface Make_Mat
          module procedure constructor_C, constructor_R
       end Interface
       Interface Clear_Mat
          module procedure Destroy_C, Destroy_R
       end Interface

       Contains
         subroutine constructor_C(Mat,N)
           type (Mat_C) :: Mat
           Integer :: N
           allocate (Mat%el(N,N))
           Mat%el = cmplx(0.0,0.0)
           Mat%dim = N
         end subroutine constructor_C

         subroutine constructor_R(Mat,N)
           type (Mat_R) :: Mat
           Integer :: N
           allocate (Mat%el(N,N))
           Mat%el = 0.0
           Mat%dim = N
         end subroutine constructor_R

         subroutine Destroy_C(Mat)
           type (Mat_C) :: Mat
           deallocate (Mat%el)
         end subroutine Destroy_C

         subroutine Destroy_R(Mat)
           type (Mat_R) :: Mat
           deallocate (Mat%el)
         end subroutine Destroy_R
     end MODULE Matrix



!!!!!!!!!!!!! Would be nice to implement one day.... !!!!!!!!!!!!!!!!!!!!!
! Use MyMats
!
! interface assignment(=)
! module procedure Equal_C
! end interface
! interface operator(*)
! module procedure Mat_mult_C
! end interface
! subroutine Equal_C(Z_out, Z_in)
! type (Mat_C), intent(in) :: Z_in
! type (MAT_C), intent(out) :: Z_out
! Z_out%A = Z_in%A
! end subroutine Equal_C
!
! function Mat_mult_C(Z1, Z2) result(Z3)
!
! type (Mat_C) , intent(in) :: Z1 , Z2
! type (Mat_C) :: Z3
! integer N
!
! N = size(Z1%A,1)
! Call Construct_Mat(Z3,N)
!
! Call MMULT(Z3%A, Z1%A, Z2%A)
!
! end function Mat_mult_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
