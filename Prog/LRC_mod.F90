!  Copyright (C) 2016 - 2018 The ALF project
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
!> This module handles the  long range Coulomb repulsion.
!> @details
!>  The intercation reads:  \f$ H_V = \frac{1}{4} \sum_{i,j} (n_i - 1) V_{i,j} (n_j -1) \f$  and the action
!> \f$ S = \sum_{i,j,\tau} \Delta \tau \phi_{i,\tau} V^{-1}_{i,j} \phi_{j,\tau} + \sum_{i} i \Delta \tau \Phi_i(n_i -1) \f$
!>
!
!--------------------------------------------------------------------

    Module LRC_mod

      Use Lattices_v3
      Use MyMats
      Use Random_wrap
      use iso_fortran_env, only: output_unit, error_unit

      Implicit none

      !> Space for the interaction matrix, orthogonal transformation, and spectrum.
      Real (Kind=Kind(0.d0)), allocatable,  Private :: V_int(:,:),  U_int(:,:), E_int(:), A_tmp(:), V_int_m1(:,:)

    contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets the Coulomb repulsion.
!>
!> @details
!> V(r)/Uhub  =        1 for r=0
!> V(r)/Uhub  = alpha*d1/r for r \= 0
!> Subroutine returns V(r)
!>
!-------------------------------------------------------------------

      Real ( Kind=Kind(0.d0) ) function  LRC_V_func(X_p, Uhub, alpha, d1)

        Implicit none

        !> Point
        Real (Kind=Kind(0.d0)), Intent(IN) :: X_p(2)
        !> Parameters
        Real (Kind=Kind(0.d0)), Intent(IN) :: Uhub,alpha, d1

        ! Local
        Real (Kind=Kind(0.d0)) :: X

        LRC_V_func = 0.d0
        X  = sqrt( X_p(1)**2   + X_p(2)**2   )
        if (  abs(X) < 1.D-10 ) then
           LRC_V_func = Uhub
        else
           LRC_V_func = Uhub*alpha*d1/X
        endif

      end function LRC_V_func

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Returns the value  private array V_int
!>
!> @details
!>
!-------------------------------------------------------------------
      Real ( Kind=Kind(0.d0) ) function  LRC_V_int(I,J)

        Implicit none

        Integer,   Intent(IN) :: I,J

        LRC_V_int = V_int(I,J)

      end function LRC_V_int



!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Prints the Coulomb repulsion as well as eigenvectors in the file Coulomb_Rep
!>
!> @details

!> @param [in] Lattice
!> \verbatim
!>  Type (Lattice)
!> \endverbatim
!> @param [in] Latt_unit
!> \verbatim
!>  Type (Unit_cell)
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>  Type  Integer(:,:)
!>  List(I=1.. Ndim,1)    =   Unit cell of site I
!>  List(I=1.. Ndim,2)    =   Orbital index  of site I
!>  Invlist(Unit_cell,Orbital) = site I
!> \endverbatim
!-------------------------------------------------------------------
      Subroutine LRC_Print(Latt, Latt_unit, list, invlist)

        Use Lattices_v3
        Implicit none

        !  Lattice
        Type (Lattice)  , intent(in) :: Latt
        !  Unit cell
        Type (Unit_cell), intent(in) :: Latt_unit
        !  List(I=1.. Ndim,1)    =   Unit cell of site I
        !  List(I=1.. Ndim,2)    =   Orbital index  of site I
        !  Invlist(Unit_cell,Orbital) = site I
        Integer, intent(in), Dimension(:,:) :: List, Invlist

        ! Local
        Real (Kind=Kind(0.d0)) :: X_p(2),  X0_p(2)
        Integer :: I,J, no_J, Ju, no_I, Iu, I0, imj

        Open (Unit = 25,file="Coulomb_Rep",status="unknown")


        I0=1
        Iu  = I0
        no_I= 1
        I = Invlist(Iu,no_I)
        Do Ju = 1, Latt%N
           do no_J = 1,Latt_unit%Norb
              J    = invlist(Ju,no_J)
              ImJ  = Latt%imj(Iu,Ju)
              X_p(:) = dble(Latt%list(ImJ,1))*latt%a1_p(:)  +  dble(Latt%list(ImJ,2))*latt%a2_p(:) + &
                   &   Latt_unit%Orb_pos_p(no_i,:)  -  Latt_unit%Orb_pos_p(no_j,:)
              Write(25,"(F16.8,2x,F16.8)") sqrt(X_p(1)**2 + X_p(2)**2), V_int(I,J)
           enddo
        Enddo
        Do J = 1,Latt%N * Latt_unit%Norb
           Write(25,*) E_int(J)
        Enddo
        Close(25)

      end Subroutine LRC_Print
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets the Coulomb repulsion.
!>
!> @details
!> \verbatim
!> Allocates and sets V_int, U_int and E_int.
!> If Norb > 2 then d1, the minimal distance, is set to the norm of Latt_unit%Orb_pos_p(2,:).
!> If Norb = 1 then d1, the minimal distance, is set to the norm of Latt%a1_p(:)
!> The intercation reads:  1/ N  \sum_{i,j} (n_i - 1) V_{i,j} (n_j -1)
!> \endverbatim
!> @param [in] Lattice
!> \verbatim
!>  Type (Lattice)
!> \endverbatim
!> @param [in] Latt_unit
!> \verbatim
!>  Type (Unit_cell)
!> \endverbatim
!> @param [in] V0,V1
!> \verbatim
!>  Type Real
!>  V(r) = Uhub  (r=0), V(r) =  Uhub*alpha*d1/r  (r \= 0)
!> \endverbatim
!> @param [in]  List, Invlist
!> \verbatim
!>  Type  Integer(:,:)
!>  List(I=1.. Ndim,1)    =   Unit cell of site I
!>  List(I=1.. Ndim,2)    =   Orbital index  of site I
!>  Invlist(Unit_cell,Orbital) = site I
!> \endverbatim
!-------------------------------------------------------------------
      Subroutine LRC_Set_VIJ(Latt, Latt_unit, Uhub, alpha, list, invlist)

        Use Lattices_v3
        Implicit none

        !  Lattice
        Type (Lattice)  , intent(in) :: Latt
        !  Unit cell
        Type (Unit_cell), intent(in) :: Latt_unit
        Real (Kind=Kind(0.d0)), intent(in) :: Uhub, alpha
        !  List(I=1.. Ndim,1)    =   Unit cell of site I
        !  List(I=1.. Ndim,2)    =   Orbital index  of site I
        !  Invlist(Unit_cell,Orbital) = site I
        Integer, intent(in), Dimension(:,:) :: List, Invlist

        !Local
        Integer ::   I,J,no_i,no_j, n, m, no, imj
        Real (Kind=Kind(0.d0)) :: X_p(2), X0_p(2), X1_p(2), d1, X, X_min, Xmean,Xmax
        Real (Kind=Kind(0.d0)), allocatable :: M_Tmp(:,:), M_Tmp1(:,:)
        Logical :: L_test=.false.


        ! Set d1, the minimal distance.
        If (Latt_unit%Norb > 1 ) then
           no = 2
           X_min = sqrt( Latt_unit%Orb_pos_p(no,1)**2  + Latt_unit%Orb_pos_p(no,2)**2 )
           d1    =  X_min
           do no = 3, Latt_unit%Norb
              X_min = sqrt( Latt_unit%Orb_pos_p(no,1)**2  + Latt_unit%Orb_pos_p(no,2)**2 )
              if (X_min <=  d1) d1 = X_min
           enddo
        else
           X_min  = sqrt( Latt%a1_p(1)**2  + Latt%a1_p(2)**2 )
           d1     =  X_min
           X_min  = sqrt( Latt%a2_p(1)**2  + Latt%a2_p(2)**2 )
           if (X_min <=  d1) d1 = X_min
        endif

        ! Allocate space
        Allocate (V_int   (Latt%N*Latt_unit%Norb, Latt%N*Latt_unit%Norb), &
             &    U_int   (Latt%N*Latt_unit%Norb, Latt%N*Latt_unit%Norb), &
             &    V_int_m1(Latt%N*Latt_unit%Norb, Latt%N*Latt_unit%Norb), &
             &    E_int   (Latt%N*Latt_unit%Norb), A_tmp(Latt%N*Latt_unit%Norb) )
        V_int = 0.d0; U_int =0.d0; E_int = 0.d0

        ! Set Potential
        Do i = 1, Latt%N
           do j = 1, Latt%N
              !Write(6,*) I,J
              imj = Latt%imj(i,j)
              X0_p = dble(Latt%list(imj,1))*Latt%a1_p + dble(Latt%list(imj,2))*Latt%a2_p
              do no_i = 1,Latt_unit%Norb
                 do no_j = 1,Latt_unit%Norb
                    n = invlist(i,no_i)
                    m = invlist(j,no_j)
                    X_p(:) = X0_p(:) +  Latt_unit%Orb_pos_p(no_i,:) - Latt_unit%Orb_pos_p(no_j,:)
                    V_int(n,m) = V_int(n,m) + LRC_V_func(X_p,Uhub,alpha,d1)
                    !Write(25,*) sqrt(X_p(1)**2 + X_p(2)**2), LRC_V_func(X_p,UHub,alpha,d1)
                 enddo
              enddo
           enddo
        enddo
        Call Diag(V_int,U_int,E_int)

        Do I = 1,size(E_int,1)
           !Write(25,*) E_int(I)
           if ( E_int(i) < 1.D-10 ) then
              Write(error_unit,*) 'LRC_Set_VIJ: V_int(i,j) is not positive definite '
              error stop 1
           endif
        enddo

        V_int_m1 = 0.d0
        Do M = 1,size(E_int,1)
           DO J = 1,size(E_int,1)
              X = U_int(j,m) / E_int(m)
              DO I = 1,size(E_int,1)
                 V_int_m1(i,j) = V_int_m1(i,j)   +  U_int(i,m) * X
              Enddo
           Enddo
        Enddo

        ! Test
        If (L_Test) then
           Allocate (M_Tmp (Latt%N*Latt_unit%Norb, Latt%N*Latt_unit%Norb))
           Allocate (M_Tmp1(Latt%N*Latt_unit%Norb, Latt%N*Latt_unit%Norb))
           M_Tmp = 0.d0; M_Tmp1 = 0.d0
           Do I = 1, Size( M_TMP,1)
              M_TMP1(I,I) = 1.d0
           Enddo
           Call MMULT(M_TMP, V_int, V_int_m1)
           Call Compare  (M_Tmp, M_Tmp1, Xmean,Xmax)
           X_min = 1.d0
           do  I = 1, Size( M_TMP,1)
              Do J = 1, Size( M_TMP,2)
                 X = Abs(V_int(I,J) - V_int(J,I)) + Abs(V_int_m1(I,J) - V_int_m1(J,I))
                 if (X <= X_min) X_min = X
              Enddo
           Enddo
           Write(6,*) 'Test LRC: ', Xmean, Xmax, X
           Deallocate (M_Tmp, M_tmp1)
        Endif


      end Subroutine LRC_Set_VIJ

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Single spin flip S0 ratio
!> @details
!> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
!> a spin flip of Operator n on time slice nt. The new value of this field is: A_n_new
!> @details
!>  @param [in]  n
!> \verbatim
!>  Integer.  Operator index of field to be fliped
!> \endverbatim
!>  @param [in]  dtau
!> \verbatim
!>  Real.  Imaginary time step
!> \endverbatim
!>  @param [in]  A_old(:)
!> \verbatim
!>  Real, dimension(:). Old configuration,  just the spacial part on the considered time slice
!> \endverbatim
!>  @param [in]  A_n_new
!> \verbatim
!>  Real.  Value of the new field for operation n.
!> \endverbatim
!>  @param [in]  N_SUN
!> \verbatim
!>  Integer.  Number of colors
!> \endverbatim
!--------------------------------------------------------------------
      Real (Kind=kind(0.d0)) function LRC_S0(n,dtau,A_old,A_n_new,N_SUN)

        Integer, intent(IN) :: n
        Integer, intent(IN) ::N_SUN
        Real (Kind=kind(0.d0)), intent(in)  :: A_n_new
        Real (Kind=kind(0.d0)), intent(in)  :: dtau
        Real (Kind=kind(0.d0)), dimension(:), intent(in)  :: A_old

        Real (Kind=Kind(0.d0)) :: X, Delta
        Integer                :: J

        Delta = A_n_new - A_old(n)
        X = 0.d0
        !Write(6,*) 'In LRC_S0:', Size(A_old,1)
        Do J = 1,Size(A_old,1)
           X = X + A_old(J)*V_int_m1(J,n) *Delta
        Enddo
        X = 2.d0 * X
        X = X + V_int_m1(n,n)*(Delta**2)

        LRC_S0 = exp( -Dtau*real(N_SUN,kind(0.d0)) * X /4.d0 )

      end function LRC_S0
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Deallocates private arrays
!--------------------------------------------------------------------
      Subroutine LRC_Clear

        Deallocate (V_int, U_int, E_int, V_int_m1, A_tmp )

      End Subroutine LRC_Clear
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Draws new fields for the long range Coulomb repulsion.
!>
!> @details
!> @param [in] A_old
!> \verbatim
!>  Real, Dimension(:)
!>  Old fields
!> \endverbatim
!> @param [out] A_new
!> \verbatim
!>  Real, Dimension(:)
!>  New fields
!> \endverbatim
!> @param [in] Percent_change
!> \verbatim
!>  Real
!>  Precent of "diagonal" fields that will be changed
!> \endverbatim
!> @param [in] Dtau
!> \verbatim
!>  Real
!>  Imaginary time step
!> \endverbatim
!--------------------------------------------------------------------
      Subroutine LRC_draw_field(Percent_change, Dtau, A_old, A_new,N_SUN)


        Implicit none

        Real (Kind=Kind(0.d0)), Intent(IN)  :: Percent_change, Dtau
        Integer               , Intent(IN)  :: N_SUN
        Real (Kind=Kind(0.d0)), Intent(IN) , dimension(:) :: A_old
        Real (Kind=Kind(0.d0)), Intent(OUT), dimension(:) :: A_new

        Integer :: n, n1,  i, m
        Real (Kind=Kind(0.d0)) :: X, Alpha, Beta

        M    =  Size(E_int,1)
        do n = 1,M
           if ( ranf_wrap() <= Percent_change )  then
              ! X = Dtau/E_int(n)
              ! X = rang(iseed)/sqrt(2.d0*gk)
              ! Distribution of X is P(X) = sqrt(gk/3.141)* exp(-gk*x**2)
              A_tmp(n) = rang_wrap()*2.d0*sqrt(E_int(n)/(2.d0*Dtau*real(N_SUN,kind(0.d0)) ))
           else
              A_tmp(n) = 0.d0
              do n1 = 1,Size(E_int,1)
                 A_tmp(n) = A_tmp(n) + U_int(n1,n)*A_old(n1)
              enddo
           endif
        enddo

        Alpha = 1.d0
        Beta  = 0.d0
        A_new = 0.d0
        Call dgemv('N', M, M, alpha, U_int, M, A_tmp, 1, beta, A_new, 1)

!!$        Real (Kind=Kind(0.d0)), allocatable :: A_test_new(:)
!!$        Logical :: Test
!!$        if (Test) then
!!$           Allocate ( A_test_new(M))
!!$           A_test_new = 0.d0
!!$           do n = 1,m
!!$              do i = 1, m
!!$                 A_test_new(i) = A_test_new(i)  + U_int(i,n)*A_tmp(n)
!!$              enddo
!!$           enddo
!!$           X = 0.d0
!!$           do i = 1,m
!!$              X = X +  ABS(A_test_new(i) -   A_new(i))
!!$           enddo
!!$           If (X >= 1.D-12 ) then
!!$              Write(6,*) X
!!$              error stop 1
!!$           Endif
!!$           Deallocate( A_test_new)
!!$        Endif

      end Subroutine LRC_Draw_Field

    end Module LRC_mod
