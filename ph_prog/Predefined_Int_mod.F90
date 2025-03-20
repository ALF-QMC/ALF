!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module provides a set of predefined interactions.
!>
!
!--------------------------------------------------------------------

module Predefined_Int

   use Operator_mod

   implicit none

contains
!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Hubbard U, SU(N)
!> U/N_{N_SUN}  [ \sum_{sigma=1}^{N_SUN}( n_{i,sigma} - 1/2 ) ]^2
!> This is for an  SU(N) code. In this case the Hubbard U  is defined only by one operator.
!--------------------------------------------------------------------
   subroutine Predefined_Int_U_SUN(OP, I, N_SUN, DTAU, U)

      implicit none
      integer, intent(In) :: I, N_SUN
      real(Kind=kind(0.d0)), intent(IN) :: Dtau, U
      type(operator), intent(Out) :: OP

      call OP_Make(Op, 1)

      Op%P(1) = I
      Op%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op%alpha = cmplx(-0.5d0, 0.d0, kind(0.d0))
      Op%g = sqrt(cmplx(-DTAU*U/(dble(N_SUN)), 0.d0, kind(0.d0)))
      Op%type = 2

      call Op_set(Op)

   end subroutine Predefined_Int_U_SUN

!-------------------------------------------------------------------
!> @Author
!> ALF-project
!
!> @brief
!> Hubbard U, SU(N)
!>   U/N_{N_SUN}  [ \sum_{sigma=1}^{N_SUN}( n_{i,sigma} - 1/2 ) ]^2
!>   The routine uses continuous   HS  transformation
!>   e^{A^2} = \int dx e^{ - x^2/2 + \sqrt{2} A x }
!>   Here A^2 =  -\Dtau U/N_{N_SUN}  [ \sum_{sigma=1}^{N_SUN}( n_{i,sigma} - 1/2 ) ]^2
!>
!--------------------------------------------------------------------
   subroutine Predefined_Int_U_SUN_continuous_HS(OP, I, N_SUN, DTAU, U)

      implicit none
      integer, intent(In) :: I, N_SUN
      real(Kind=kind(0.d0)), intent(IN) :: Dtau, U
      type(operator), intent(Out) :: OP

      call OP_Make(Op, 1)

      Op%P(1) = I
      Op%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op%alpha = cmplx(-0.5d0, 0.d0, kind(0.d0))
      Op%g = sqrt(cmplx(-DTAU*U*2.d0/(dble(N_SUN)), 0.d0, kind(0.d0)))
      Op%type = 3

      call Op_set(Op)

   end subroutine Predefined_Int_U_SUN_Continuous_HS

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets  Hubbard U Mz
!>   - U/2  ( (n_{i,up} -1/2) - (n_{i,do} - 1/2)  ) ^2
!>   Here N_FL = 2 and  the routine sets both the up and down operators.
!>   The routine uses a  descrete  HS transformation  with four fields
!>   (see documentation)
!--------------------------------------------------------------------
   subroutine Predefined_Int_U_MZ(OP_up, Op_do, I, DTAU, U)

      implicit none
      integer, intent(In) :: I
      real(Kind=kind(0.d0)), intent(IN) ::  Dtau, U
      type(operator), intent(Out) :: OP_up, Op_do

      call OP_Make(Op_up, 1)
      call OP_Make(Op_do, 1)

      Op_up%P(1) = I
      Op_up%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op_up%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op_up%g = sqrt(cmplx(DTAU*U/2.d0, 0.d0, kind(0.d0)))
      Op_up%alpha = cmplx(-0.5d0, 0.d0, kind(0.d0))
      Op_up%type = 2

      Op_do%P(1) = I
      Op_do%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op_do%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op_do%g = -sqrt(cmplx(DTAU*U/2.d0, 0.d0, kind(0.d0)))
      Op_do%alpha = cmplx(-0.5d0, 0.d0, kind(0.d0))
      Op_do%type = 2

      call Op_set(Op_up)
      call Op_set(Op_do)

   end subroutine Predefined_Int_U_MZ
!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets  Hubbard U Mz
!>   - U/2  ( n_{i,up} - n_{i,do} ) ^2
!>   Here N_FL = 2 and  the routine sets both the up and down operators.
!>   The routine uses continuous   HS  transformation
!>   e^{A^2} = \int dx e^{ - x^2/2 + \sqrt{2} A x }
!>   Here A^2 =  \Dtau U /2  ( n_{i,up} - n_{i,do} ) ^2
!--------------------------------------------------------------------
   subroutine Predefined_Int_U_MZ_continuous_HS(OP_up, Op_do, I, DTAU, U)

      implicit none
      integer, intent(In) :: I
      real(Kind=kind(0.d0)), intent(IN) ::  Dtau, U
      type(operator), intent(Out) :: OP_up, Op_do

      call OP_Make(Op_up, 1)
      call OP_Make(Op_do, 1)

      Op_up%P(1) = I
      Op_up%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op_up%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op_up%g = sqrt(cmplx(DTAU*U, 0.d0, kind(0.d0)))
      Op_up%type = 3

      Op_do%P(1) = I
      Op_do%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op_do%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op_do%g = -sqrt(cmplx(DTAU*U, 0.d0, kind(0.d0)))
      Op_do%type = 3

      call Op_set(Op_up)
      call Op_set(Op_do)

   end subroutine Predefined_Int_U_MZ_Continuous_HS

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets
!>   - V/N_SUN  [ \sum_{s=1}^{N_SUN}( c^{dag}_{i,s} c_{j,s} + c^{dag}_{j,s} c_{i,s} ) ]^2
!>
!>   This is for an SU(N) symmetric code such that the above  interaction corresponds to just one operator
!>
!--------------------------------------------------------------------
   subroutine Predefined_Int_V_SUN(OP, I, J, N_SUN, DTAU, V)

      implicit none
      integer, intent(In) :: I, J, N_SUN
      real(Kind=kind(0.d0)), intent(IN) ::  Dtau, V
      type(operator), intent(Out) :: OP

      call OP_Make(Op, 2)

      Op%P(1) = I
      Op%P(2) = J
      Op%O(1, 2) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op%O(2, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op%g = sqrt(cmplx(DTAU*V/real(N_SUN, kind(0.d0)), 0.d0, kind(0.d0)))
      Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op%type = 2

      call Op_set(Op)

   end subroutine Predefined_Int_V_SUN

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets
!>   - V/N_SUN  [ \sum_{s=1}^{N_SUN}( ic^{dag}_{i,s} c_{j,s}  - ic^{dag}_{j,s} c_{i,s} ) ]^2
!>
!>   This is for an SU(N) symmetric code such that the above  interaction corresponds to just one operator
!>
!--------------------------------------------------------------------
   subroutine Predefined_Int_VJ_SUN(OP, I, J, N_SUN, DTAU, V)

      implicit none
      integer, intent(In) :: I, J, N_SUN
      real(Kind=kind(0.d0)), intent(IN) ::  Dtau, V
      type(operator), intent(Out) :: OP

      call OP_Make(Op, 2)

      Op%P(1) = I
      Op%P(2) = J
      Op%O(1, 2) = cmplx(0.d0, 1.d0, kind(0.d0))
      Op%O(2, 1) = cmplx(0.d0, -1.d0, kind(0.d0))
      Op%g = sqrt(cmplx(DTAU*V/real(N_SUN, kind(0.d0)), 0.d0, kind(0.d0)))
      Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op%type = 2

      call Op_set(Op)

   end subroutine Predefined_Int_VJ_SUN

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets
!>     \hat{Z}_{i,j} \xi  \sum_{s=1}^{N_SUN}( c^{dag}_{i,s} c_{j,s} + c^{dag}_{j,s} c_{i,s} )
!>
!>   This is for an SU(N) symmetric code, such that the above corresponds to one operator.
!>
!--------------------------------------------------------------------
   subroutine Predefined_Int_Ising_SUN(OP, I, J, DTAU, Xi)

      implicit none
      integer, intent(In) :: I, J
      real(Kind=kind(0.d0)), intent(IN) ::  Dtau, Xi
      type(operator), intent(Out) :: OP

      call OP_Make(Op, 2)

      Op%P(1) = I
      Op%P(2) = J
      Op%O(1, 2) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op%O(2, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op%g = cmplx(-dtau*xi, 0.d0, kind(0.d0))
      Op%alpha = cmplx(0d0, 0.d0, kind(0.d0))
      Op%type = 1

      call Op_set(Op)

   end subroutine Predefined_Int_Ising_SUN

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets
!>       \sum_{s=1}^{N_SUN} \Phi i ( c^{dag}_{i,s} c_{i,s} - 0.5 )
!>
!>   For the long range Coulomb repulsion.
!>
!--------------------------------------------------------------------
   subroutine Predefined_Int_LRC(OP, I, DTAU)

      implicit none
      integer, intent(In) :: I
      real(Kind=kind(0.d0)), intent(IN) ::  Dtau
      type(operator), intent(Out) :: OP

      call OP_Make(Op, 1)

      Op%P(1) = I
      Op%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op%alpha = cmplx(-0.5d0, 0.d0, kind(0.d0))
      Op%g = cmplx(0.d0, Dtau, kind(0.d0))
      Op%type = 3

      call Op_set(Op)

   end subroutine Predefined_Int_LRC

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets  Jz-Jz interaction
!>     - |J_z|/2  ( S^{z}_i - sign|J_z| S^{z}_j ) ^2 =
!>       J_z  S^{z}_i  S^{z}_j  -   |J_z|/2   (S^{z}_i)^2 - |J_z|/2   (S^{z}_j)^2
!>   Here N_FL = 2 and the routine sets both the up and down operators.
!>   If  particle fluctuations are frozen on the i and j sites then (S^{z}_i)^2 = 1/4
!>   and the interactions corresponds to a Jz-Jz ferro or antiferro coupling.
!--------------------------------------------------------------------
   subroutine Predefined_Int_Jz(OP_up, Op_do, I, J, DTAU, Jz)

      implicit none
      integer, intent(In) :: I, J
      real(Kind=kind(0.d0)), intent(IN) ::  Dtau, Jz
      type(operator), intent(Out) :: OP_up, Op_do

      call OP_Make(Op_up, 2)
      call OP_Make(Op_do, 2)

      Op_up%P(1) = I
      Op_up%P(2) = J
      Op_up%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op_up%O(2, 2) = cmplx(-Jz/abs(Jz), 0.d0, kind(0.d0))
      Op_up%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op_up%g = sqrt(cmplx(DTAU*Jz/8.d0, 0.d0, kind(0.d0)))
      Op_up%type = 2

      Op_do%P(1) = I
      Op_do%P(2) = J
      Op_do%O(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
      Op_do%O(2, 2) = cmplx(-Jz/abs(Jz), 0.d0, kind(0.d0))
      Op_do%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op_do%g = -sqrt(cmplx(DTAU*Jz/8.d0, 0.d0, kind(0.d0)))
      Op_do%type = 2

      call Op_set(Op_up)
      call Op_set(Op_do)

   end subroutine Predefined_Int_Jz

!!$       Still has to be done when implementing Hamiltoninas.
!!$      Subroutine Predefined_Int_J_SUN( )
!!$        Implicit none
!!$      end Subroutine Predefined_Int_J_SUN

end module Predefined_Int
