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
!> This module provides a set of predefined interactions.
!>       
!
!--------------------------------------------------------------------

    Module Predefined_Int
      
      Use Operator_mod
      
      Implicit none
      
      
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
      Subroutine Predefined_Int_U_SUN( OP, I, N_SUN, DTAU, U  )
        
        Implicit none
        Integer,  Intent(In) :: I, N_SUN 
        Real (Kind=Kind(0.d0)), Intent(IN) :: Dtau, U
        Type(Operator), Intent(Out) :: OP

        Call OP_Make( Op,1 )

        Op%P(1)   = I
        Op%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
        Op%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
        Op%g      = SQRT(CMPLX(-DTAU*U/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
        Op%type   = 2
        
        Call Op_set( Op )
        
      end Subroutine Predefined_Int_U_SUN

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
      Subroutine Predefined_Int_U_SUN_continuous_HS( OP, I, N_SUN, DTAU, U  )
        
        Implicit none
        Integer,  Intent(In) :: I, N_SUN 
        Real (Kind=Kind(0.d0)), Intent(IN) :: Dtau, U
        Type(Operator), Intent(Out) :: OP

        Call OP_Make( Op,1 )

        Op%P(1)   = I
        Op%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
        Op%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
        Op%g      = SQRT(CMPLX(-DTAU*U*2.d0/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
        Op%type   = 3
        
        Call Op_set( Op )
        
      end Subroutine Predefined_Int_U_SUN_Continuous_HS

      
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
      Subroutine Predefined_Int_U_MZ (  OP_up,  Op_do, I, DTAU, U  ) 

        Implicit none
        Integer,  Intent(In) :: I
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau, U
        Type(Operator), Intent(Out) :: OP_up, Op_do

        Call OP_Make( Op_up,1 )
        Call OP_Make( Op_do,1 )
        
        Op_up%P(1)   =  I
        Op_up%O(1,1) =  cmplx(1.d0  ,0.d0, kind(0.D0))
        Op_up%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_up%g      =  SQRT(CMPLX(DTAU*U/2.d0, 0.D0, kind(0.D0))) 
        Op_up%alpha  =  cmplx(-0.5d0,0.d0, kind(0.D0))
        Op_up%type   =  2

        Op_do%P(1)   =  I
        Op_do%O(1,1) =  cmplx(1.d0  ,0.d0, kind(0.D0))
        Op_do%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_do%g      = -SQRT(CMPLX(DTAU*U/2.d0, 0.D0, kind(0.D0))) 
        Op_do%alpha  =  cmplx(-0.5d0,0.d0, kind(0.D0))
        Op_do%type   =  2

        Call Op_set( Op_up )
        Call Op_set( Op_do )
        
      end Subroutine Predefined_Int_U_MZ
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
      Subroutine Predefined_Int_U_MZ_continuous_HS(  OP_up,  Op_do, I, DTAU, U  ) 

        Implicit none
        Integer,  Intent(In) :: I
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau, U
        Type(Operator), Intent(Out) :: OP_up, Op_do

        Call OP_Make( Op_up,1 )
        Call OP_Make( Op_do,1 )
        
        Op_up%P(1)   =  I
        Op_up%O(1,1) =  cmplx(1.d0  ,0.d0, kind(0.D0))
        Op_up%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_up%g      =  SQRT(CMPLX(DTAU*U, 0.D0, kind(0.D0))) 
        Op_up%type   =  3

        Op_do%P(1)   =  I
        Op_do%O(1,1) =  cmplx(1.d0  ,0.d0, kind(0.D0))
        Op_do%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_do%g      = -SQRT(CMPLX(DTAU*U, 0.D0, kind(0.D0))) 
        Op_do%type   =  3

        Call Op_set( Op_up )
        Call Op_set( Op_do )
        
      end Subroutine Predefined_Int_U_MZ_Continuous_HS

      
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
      Subroutine Predefined_Int_V_SUN( OP, I, J, N_SUN, DTAU, V  ) 
        
        Implicit none
        Integer,  Intent(In) :: I, J, N_SUN
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau, V
        Type(Operator), Intent(Out) :: OP

        Call OP_Make( Op,2 )

        Op%P(1)   = I
        Op%P(2)   = J
        Op%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
        Op%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
        Op%g      = SQRT(CMPLX(DTAU*V/real(N_SUN,kind(0.d0)), 0.D0, kind(0.D0))) 
        Op%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
        Op%type   = 2
        
        Call Op_set( Op )
        
      end Subroutine Predefined_Int_V_SUN

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
      Subroutine Predefined_Int_VJ_SUN( OP, I, J, N_SUN, DTAU, V  ) 
        
        Implicit none
        Integer,  Intent(In) :: I, J, N_SUN
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau, V
        Type(Operator), Intent(Out) :: OP

        Call OP_Make( Op,2 )

        Op%P(1)   = I
        Op%P(2)   = J
        Op%O(1,2) = cmplx(0.d0 ,  1.d0, kind(0.D0)) 
        Op%O(2,1) = cmplx(0.d0 , -1.d0, kind(0.D0))
        Op%g      = SQRT(CMPLX(DTAU*V/real(N_SUN,kind(0.d0)), 0.D0, kind(0.D0))) 
        Op%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
        Op%type   = 2
        
        Call Op_set( Op )
        
      end Subroutine Predefined_Int_VJ_SUN

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
      Subroutine Predefined_Int_Ising_SUN( OP, I, J, DTAU, Xi  ) 
        
        Implicit none
        Integer,  Intent(In) :: I, J
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau, Xi
        Type(Operator), Intent(Out) :: OP

        Call OP_Make( Op,2 )

        Op%P(1)   = I
        Op%P(2)   = J
        Op%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
        Op%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
        Op%g      = cmplx(-dtau*xi,0.D0,kind(0.D0))
        Op%alpha  = cmplx(0d0,0.d0, kind(0.D0)) 
        Op%type   = 1
        
        Call Op_set( Op )
        
      end Subroutine Predefined_Int_Ising_SUN

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
      Subroutine Predefined_Int_LRC( OP, I, DTAU  ) 
        
        Implicit none
        Integer,  Intent(In) :: I
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau
        Type(Operator), Intent(Out) :: OP

        Call OP_Make( Op,1 )


        Op%P(1)   = I
        Op%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
        Op%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
        Op%g      = cmplx(0.d0  ,Dtau, kind(0.D0)) 
        Op%type   = 3
        
        Call Op_set( Op )
        
      end Subroutine Predefined_Int_LRC

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
      Subroutine Predefined_Int_Jz (  OP_up,  Op_do, I, J,  DTAU, Jz  ) 

        Implicit none
        Integer,  Intent(In) :: I,J
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau, Jz
        Type(Operator), Intent(Out) :: OP_up, Op_do

        Call OP_Make( Op_up,2 )
        Call OP_Make( Op_do,2 )
        
        Op_up%P(1)   =  I
        Op_up%P(2)   =  J
        Op_up%O(1,1) =  cmplx(        1.d0  ,0.d0, kind(0.D0))
        Op_up%O(2,2) =  cmplx(- Jz/Abs(Jz)  ,0.d0, kind(0.D0))
        Op_up%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_up%g      =  SQRT(CMPLX(DTAU*Abs(Jz)/8.d0, 0.D0, kind(0.D0))) 
        Op_up%type   =  2

        Op_do%P(1)   =  I
        Op_do%P(2)   =  J
        Op_do%O(1,1) =  cmplx(        1.d0  ,0.d0, kind(0.D0))
        Op_do%O(2,2) =  cmplx(- Jz/Abs(Jz)  ,0.d0, kind(0.D0))
        Op_do%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_do%g      = -SQRT(CMPLX(DTAU*Abs(Jz)/8.d0, 0.D0, kind(0.D0))) 
        Op_do%type   =  2

        Call Op_set( Op_up ) 
        Call Op_set( Op_do )
        
      end Subroutine Predefined_Int_Jz

!!$       Still has to be done when implementing Hamiltoninas.
!!$      Subroutine Predefined_Int_J_SUN( )
!!$        Implicit none
!!$      end Subroutine Predefined_Int_J_SUN


!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Sets a spin-spin interaction term of the form \f$ J S_a S_b \f$  in the operator format.
!> Here, a and b can be x,y or z.
!> This formulation assumes that all charges fluctuations are frozen out, e.g. by a large Hubbard U.
!> Thereby, \f$ c^\dagger_{i,\uparrow} c_{i,\uparrow} + c^\dagger_{i,\downarrow} c_{i,\downarrow} = 1 \f$.
!> The spin operators are represented in terms of fermionic operators as follows:
!> \f$ S^{\alpha}_i = \frac{1}{2}(c^\dagger_i \sigma^\alpha c_i)\f$
!> Here, \f$ c^\dagger_i = (c^\dagger_{i,\uparrow}, c^\dagger_{i,\downarrow}) \f$ is a spinor of fermionic creation operators and \f$ \sigma^\alpha \f$ are the Pauli matrices.
!> The exact term implemented here is \f$g \abs{J} / 8 (S^a_i + conjg(g) J/\abs{J} S^b_j)^2 = S^a_i S^b_j + const.\f$.
!> The term includes a gauge factor g = gauge, which can possibly be used to alleviate the sign problem.
!> The gauge factor is a pure phase factor, defined as \f$ g = e^{i \phi} \f$ with \f$ \phi = \text{gauge\_phi} \pi \f$.
!--------------------------------------------------------------------
        subroutine Predefined_Int_SaSb(Op, a, b, I1_up, I2_up, I1_do, I2_do, J, dtau, gauge_phi)
          implicit none
          type(Operator), intent(out) :: Op
          Character (len=1), intent(in) :: a, b
          integer, intent(in) :: I1_up, I2_up, I1_do, I2_do
          real(Kind=Kind(0.d0)), intent(in) :: J, dtau
          real(Kind=Kind(0.d0)), intent(in) :: gauge_phi

          real(Kind=Kind(0.d0)), parameter :: pi = acos(-1.d0)
          real(Kind=Kind(0.d0)) :: J_sgn
          complex(Kind=Kind(0.d0)) :: gauge

            gauge = cmplx(cos(gauge_phi*pi), sin(gauge_phi*pi), kind(0.D0))
            if (abs(abs(gauge)-1.d0) > 1.d-8) then
                write(error_unit, *) 'Error in Predefined_Int_SaSb: gauge must be a pure phase factor with absolute value 1.'
                call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
            endif
            if (abs(J).le.1D-12)then
              J_sgn = 0.d0
            else
              J_sgn = J/abs(J)
            end if

            call Op_make(Op,4)
            Op%P(1) = I1_up
            Op%P(2) = I2_up
            Op%P(3) = I1_do
            Op%P(4) = I2_do
            select case (a)
            case('x')
               Op%O(1,3) = cmplx(1.d0, 0.d0, kind(0.D0))
               Op%O(3,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
            case('y')
               Op%O(1,3) = cmplx(0.d0, -1.d0, kind(0.D0))
               Op%O(3,1) = cmplx(0.d0 ,1.d0, kind(0.D0))
            case('z')
               Op%O(1,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
               Op%O(3,3) = cmplx(-1.d0 ,0.d0, kind(0.D0))
            end select
            select case (b)
            case('x')
               Op%O(2,4) = conjg(gauge)*cmplx(J_sgn*1.d0 ,0.d0, kind(0.D0))
               Op%O(4,2) = conjg(gauge)*cmplx(J_sgn*1.d0 ,0.d0, kind(0.D0))
            case('y')
               Op%O(2,4) = conjg(gauge)*cmplx(0.d0, -J_sgn*1.d0, kind(0.D0))
               Op%O(4,2) = conjg(gauge)*cmplx(0.d0 ,J_sgn*1.d0, kind(0.D0))
            case('z')
               Op%O(2,2) = conjg(gauge)*cmplx(J_sgn*1.d0 ,0.d0, kind(0.D0))
               Op%O(4,4) = conjg(gauge)*cmplx(-J_sgn*1.d0 ,0.d0, kind(0.D0))
            end select
            Op%g     = SQRT(gauge*CMPLX(-dtau*abs(J)/8d0, 0.D0, kind(0.D0)))
            Op%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
            Op%type  = 2
            call Op_set(Op)
        end subroutine Predefined_Int_SaSb


     end Module Predefined_Int
