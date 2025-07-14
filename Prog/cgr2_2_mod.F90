!  Copyright (C) 2016 - 2022  The ALF project
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

module cgr2_2_mod
   implicit none
   private
   public :: CGR2_2, solve_extended_system, get_blocks
   contains

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief
!> This function separates out function blocks from the matrix inp
!> The assumed layout for the result is:
!>  (A, B)
!>  (C, D)
!>
!> @param[in] inp The matrix that is separated into four equal sized blocks
!> @param[inout] A
!> @param[inout] B
!> @param[inout] C
!> @param[inout] D
!> @param[in] LQ The dimension of the matrices A, B, C, D
!--------------------------------------------------------------------
      Subroutine get_blocks(A, B, C, D, INP, LQ)
        Implicit none
        Integer, intent(in) :: LQ
        Complex (Kind=Kind(0.D0)), intent(inout) :: A(LQ, LQ), B(LQ, LQ), C(LQ, LQ), D(LQ, LQ)
        Complex (Kind=Kind(0.D0)), intent(in) :: INP(2*LQ, 2*LQ)
        Integer :: I, J, I1, J1

        DO I = 1,LQ
          I1 = I+LQ
          DO J = 1,LQ
              J1 = J + LQ
              A(I,J) = INP(I ,J )
              D(I,J) = INP(I1,J1)
              C(I,J) = INP(I1,J )
              B(I,J) = INP(I,J1 )
          ENDDO
        ENDDO

      end Subroutine

#if (defined(STAB2) || defined(STAB1)) && !defined(STABLOG)
!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief
!> This functions constructs an extended matrix H1 from UCT and VINV for use as the rhs of
!> an equation. Then the solution X of the System V3 * X = H1 is sought.
!> Then a scaling by D3 and a multiplication by U3 is performed.
!> Usually the matrices U3, D3, V3 stem from a UDV decomposition.
!
!> The rationale for constructing this extended matrix is that Fakher says it's more stable.
!>
!> @param[inout] HLP The result matrix
!> @param[in] UCT Matrix, dimension(LQ, LQ)
!> @param[in] V1INV Matrix, dimension(LQ, LQ)
!> @param[in] U3 Matrix, dimension(2*LQ, 2*LQ)
!> @param[in] D3 Matrix, dimension(2*LQ, 2*LQ)
!> @param[in] V3 Matrix, dimension(2*LQ, 2*LQ)
!> @param[in] LQ The dimension of the matrices UCT and V1INV
!--------------------------------------------------------------------
           Subroutine solve_extended_System(HLP, UCT, VINV, U3, D3, V3, LQ)
           Use MyMats
           Implicit none
           Integer, intent(in) :: LQ
           Complex (Kind=Kind(0.D0)), intent(in) :: U3(2*LQ, 2*LQ), D3(2*LQ), V3(2*LQ,2*LQ), UCT(LQ,LQ), VINV(LQ,LQ)
           Complex (Kind=Kind(0.D0)), intent(inout) ::HLP(2*LQ, 2*LQ)
           Complex (Kind=Kind(0.D0)), Allocatable, Dimension(:) :: TMPVEC
           Complex (Kind=Kind(0.D0)), Allocatable, Dimension(:, :) :: HLPB1
           INTEGER, Dimension(:), Allocatable :: IPVT
           Integer :: LQ2, info, I, j
           Complex (Kind=Kind(0.D0)) :: zero

           zero = 0.D0
           LQ2 = 2*LQ
           ALLOCATE(TMPVEC(LQ2), HLPB1(LQ2, LQ2), IPVT(LQ2))
           TMPVEC = conjg(1.D0/D3)
           ! set HLPB1 equal to zero
           call zlaset('A', LQ2, LQ2, zero, zero, HLPB1, LQ2)
           DO J = 1, LQ
              DO I = 1, LQ
                 HLPB1(I   , J    ) =  UCT(I, J)
                 HLPB1(I+LQ, J+LQ ) =  VINV(I,J)
              ENDDO
           ENDDO
           CALL ZGETRF(LQ2, LQ2, V3, LQ2, IPVT, info)
           CALL ZGETRS('C', LQ2, LQ2, V3, LQ2, IPVT, HLPB1, LQ2, info)! Block structure of HLPB1 is not exploited
           DO J = 1,LQ2
              DO I = 1,LQ2
                 HLPB1(I,J)  = TMPVEC(I)*HLPB1(I,J)
              ENDDO
           ENDDO
           deallocate (TMPVEC)
           CALL MMULT(HLP, U3, HLPB1)
           deallocate(HLPB1)
         end subroutine solve_extended_System
#else
!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief
!> This functions constructs an extended matrix H1 from UCT and VINV for use as the rhs of
!> an equation. Then the solution X of the System P^* V^* * X = H1 is sought.
!> Then a scaling by D and a multiplication by a unitary matrix is performed.
!
!> The rationale for constructing this extended matrix is that Fakher says it's more stable.
!>
!> @param[inout] HLP The result matrix
!> @param[in] UCT Matrix, dimension(LQ, LQ)
!> @param[in] V1INV Matrix, dimension(LQ, LQ)
!> @param[in] A The input matrix that contains the Householder reflectors on its subdiagonal
!>              and the upper triangular matrix R in its upper part as returned by ?GETRF
!> @param[in] D A vector containing the scales.
!> @param[in] TAU A vector containing the scalar factors of the Householder decomposition
!> @param[in] IPVT A permutation vector usable by ?LAPMR/?LAPMT
!> @param[in] LQ The dimension of the matrices UCT and V1INV
!> @param[in] WORK work space
!> @param[in] LWORK The size of the work space.
!--------------------------------------------------------------------
         Subroutine solve_extended_System(HLP, UCT, VINV, A, D, TAU, PIVT, LQ, WORK, LWORK)
           Implicit none
           Integer, intent(in) :: LQ, LWORK
           Complex (Kind=Kind(0.D0)), intent(in) :: A(2*LQ, 2*LQ), D(2*LQ), TAU(2*LQ), UCT(LQ,LQ), VINV(LQ,LQ), WORK(LWORK)
           Integer, intent(in) :: PIVT(2*LQ)
           Complex (Kind=Kind(0.D0)), intent(inout) :: HLP(2*LQ, 2*LQ)
           Complex (Kind=Kind(0.D0)), Allocatable, Dimension(:) :: TMPVEC
           Integer :: LQ2, info, I, j
           Complex (Kind=Kind(0.D0)) :: z
           LOGICAL :: FORWRD

           z = 0.D0
           LQ2 = 2*LQ
           ALLOCATE(TMPVEC(LQ2))
           TMPVEC = conjg(1.D0/D)
           ! set HLP equal to zero
           call zlaset('A', LQ2, LQ2, z, z, HLP, LQ2)
           DO I = 1, LQ
              DO J = 1, LQ
                 HLP(I   , J    ) =  UCT(I, J)
                 HLP(I+LQ, J+LQ ) =  VINV(I,J)
              ENDDO
           ENDDO
           z = 1.D0
           !P * V3^* X = H1
           FORWRD = .true.
           CALL ZLAPMR(FORWRD, LQ2, LQ2, HLP(1, 1), LQ2, PIVT(1))
           CALL ZTRSM('L', 'U', 'C', 'N', LQ2, LQ2, z, A(1, 1), LQ2, HLP(1, 1), LQ2)! Block structure of HLPB1 is not exploited
           DO J = 1,LQ2
              DO I = 1,LQ2
                 HLP(I, J)  = TMPVEC(I)*HLP(I, J)
              ENDDO
           ENDDO
           deallocate (TMPVEC)
           ! Res = U * HLP
           CALL ZUNMQR('L', 'N', LQ2, LQ2, LQ2, A(1, 1), LQ2, TAU(1), HLP(1, 1), LQ2, WORK(1), LWORK, INFO)
         end subroutine solve_extended_System
#endif

!--------------------------------------------------------------------

      SUBROUTINE CGR2_2(GRT0, GR00, GRTT, GR0T, udv2, udv1, LQ)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given      B2 = U2*D2*V2 is right (i.e. from time slice 0 to tau) propagation to time tau 
!>            B1 = V1*D1*U1 is left  (i.e. from time slice Ltrot to tau) propagation to time tau
!> Calc:      (  1   B1 )^-1       ( G00    G0T )   
!>            (-B2   1  )     ==   ( GT0    GTT )
!>            
!>       G00 = (1 + B1*B2)^-1   G0T = -(1 -  G00)*B2^-1
!>       GT0  =   B2 * G00      GTT = (1 + B2*B1)^-1
!>
!>(  1   V1*D1*U1 )^-1      ( ( V1   0 )   ( V1^-1     D1*U1 ) )^-1
!>(-U2*D2*V2   1  )      == ( ( 0   U2 ) * (-D2*V2     U2^-1 ) )       == I
!>  Transpose before carrying out the singular value decomposition 
!>
!>     
!>       ( ( V1   0 )   ( V1^-1     D1*U1 )^*^* )^-1                      ( V1^-1    0 ) 
!> I ==  ( ( 0   U2 ) * (-D2*V2     U2^-1 )     )      =   (UDV^*)^(-1) * ( 0     U2^-1) =
!>
!>                                      ( V1^-1    0 )
!>  ==   U * D^(*,-1) * V^(*,-1)     *  ( 0     U2^-1)
!
!--------------------------------------------------------------------

        Use MyMats
        Use UDV_State_mod
#if (defined(STAB2) || defined(STAB1)) && !defined(STABLOG)   

        Use UDV_WRAP_mod
        Implicit none

        !  Arguments
        Integer,  intent(in) :: LQ
        CLASS(UDV_State), intent(in) :: udv1, udv2
        Complex (Kind=Kind(0.d0)), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)

        ! Local::
        Complex  (Kind=Kind(0.d0)), allocatable :: V1INV(:,:)
        Complex  (Kind=Kind(0.d0)), allocatable :: D3B(:)
        Complex  (Kind=Kind(0.d0)) :: Z, alpha, beta
        Complex(Kind = Kind(0.D0)), allocatable, Dimension(:, :) :: MYU2, HLPB1, HLPB2, U3B, V3B
        Integer :: LQ2, I,J, NCON

        allocate( V1INV(LQ,LQ), D3B(2*LQ) )
        
        if(udv1%side .ne. "L" .and. udv1%side .ne. "l" ) then
          write(*,*) "calling wrong decompose"
        endif
        if(udv2%side .ne. "R" .and. udv2%side .ne. "r" ) then
          write(*,*) "calling wrong decompose"
        endif
        
        LQ2 = LQ*2
        NCON = 0
        alpha = 1.D0
        beta = 0.D0
        ALLOCATE(MYU2(LQ, LQ), HLPB1(LQ2, LQ2), HLPB2(LQ2, LQ2), U3B(LQ2, LQ2), V3B(LQ2, LQ2))
        MYU2 = CONJG(TRANSPOSE(udv2%U))
        CALL INV(udv1%V,V1INV,Z)
        If (dble(udv1%D(1)) >  dble(udv2%D(1)) ) Then 
           !Write(6,*) "D1(1) >  D2(1)", dble(D1(1)), dble(D2(1))
           call zlacpy('A',LQ,LQ,V1INV(1,1) ,LQ,HLPB2(1   ,1   ),LQ2)
           call zlacpy('A',LQ,LQ,MYU2(1,1)  ,LQ,HLPB2(1+LQ,1+LQ),LQ2)
           DO J = 1,LQ
              DO I = 1,LQ
                 HLPB2(I   , J+LQ ) =  udv1%D(I)*conjg(udv1%U(J,I))!udv1%U(I,J)
                 HLPB2(I+LQ, J    ) = -udv2%D(I)*udv2%V(I,J)
              ENDDO
           ENDDO
           HLPB1 = CT(HLPB2)

           !CALL UDV_wrap(HLPB1,U3B,D3B,V3B,NCON)
           CALL UDV_wrap_Pivot(HLPB1,U3B,D3B,V3B,NCON,LQ2,LQ2)
! !!$!!!!!!!!!!!!!  Tests
! !!$        Xmax = 0.d0
! !!$        DO I = 1,LQ2
! !!$           DO J = 1,LQ2
! !!$              Z = cmplx(0.d0,0.d0)
! !!$              DO N = 1,LQ2
! !!$                 Z = Z + U3B(I,N) *conjg(U3B(J,N))
! !!$              ENDDO
! !!$              if (I == J)  Z = Z - cmplx(1.d0,0.d0)
! !!$              X = real(SQRT( Z* conjg(Z)),Kind=Kind(0.d0))
! !!$              if (X > Xmax) Xmax = X
! !!$           ENDDO
! !!$        ENDDO
! !!$        !Write(6,*) 'Cgr2_2, ortho: ', Xmax
! !!$        DO I = 1,LQ2
! !!$           Z =  D3B(I)
! !!$           if (I == 1)  Xmax = real(SQRT( Z* conjg(Z)),Kind=Kind(0.d0)) 
! !!$           if ( real(SQRT( Z* conjg(Z)),Kind=Kind(0.d0))  < Xmax ) Xmax = &
! !!$                & real(SQRT( Z* conjg(Z)),Kind=Kind(0.d0))
! !!$        ENDDO
! !!$        !Write(6,*) 'Cgr2_2, Cutoff: ', Xmax
! !!$!!!!!!!!!!!!! End Tests
           call solve_extended_System(HLPB1, V1INV, MYU2, U3B, D3B, V3B, LQ)
           call get_blocks(GR00, GR0T, GRT0, GRTT, HLPB1, LQ)
        Else
           !Write(6,*) "D1(1) <  D2(1)", dble(D1(1)), dble(D2(1))
           call zlacpy('A',LQ,LQ,MYU2(1,1) ,LQ,HLPB2(1   ,1   ),LQ2)
           call zlacpy('A',LQ,LQ,V1INV(1,1),LQ,HLPB2(1+LQ,1+LQ),LQ2)
           DO J = 1,LQ
              DO I = 1,LQ
                 HLPB2(I   , J+LQ ) = -udv2%D(I)*udv2%V(I,J)
                 HLPB2(I+LQ, J    ) =  udv1%D(I)*conjg(udv1%U(J,I))!udv1%U(I,J)
              ENDDO
           ENDDO
           ! HLPB1 = CT(HLPB2)
           HLPB1 = conjg(transpose(HLPB2))
           
           !CALL UDV_wrap(HLPB1,U3B,D3B,V3B,NCON)
           CALL UDV_wrap_Pivot(HLPB1,U3B,D3B,V3B,NCON,LQ2,LQ2)
           call solve_extended_System(HLPB1, MYU2, V1INV, U3B, D3B, V3B, LQ)
           call get_blocks(GRTT, GRT0, GR0T, GR00, HLPB1, LQ)
        Endif
        DEALLOCATE(MYU2, HLPB1, HLPB2, U3B, V3B)
#else
        Use QDRP_mod
        Implicit none

        !  Arguments
        Integer,  intent(in) :: LQ
        CLASS(UDV_State), intent(in) :: udv1, udv2
        Complex (Kind=Kind(0.d0)), intent(inout) :: GRT0(LQ,LQ), GR0T(LQ,LQ), GR00(LQ,LQ), GRTT(LQ,LQ)


        ! Local::
        Complex  (Kind=Kind(0.d0)), allocatable, Dimension(:) :: D3, D1m, D2m
        Complex  (Kind=Kind(0.d0)) :: Z
        Complex(Kind = Kind(0.D0)), allocatable, Dimension(:, :) :: MYU2, HLPB1, HLPB2, V1INV
        Integer :: LQ2, I, J, LWORK
        
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK
        INTEGER, Dimension(:), Allocatable :: IPVT
        
        if(udv1%side .ne. "L" .and. udv1%side .ne. "l" ) then
          write(*,*) "calling wrong decompose"
        endif
        if(udv2%side .ne. "R" .and. udv2%side .ne. "r" ) then
          write(*,*) "calling wrong decompose"
        endif
        
        LQ2 = LQ*2
        ALLOCATE(MYU2(LQ, LQ), V1INV(LQ,LQ), HLPB1(LQ2, LQ2), HLPB2(LQ2, LQ2), D3(LQ2), D1m(LQ), D2m(LQ))
        Allocate(IPVT(LQ2), TAU(LQ2))
        IPVT = 0
        MYU2 = CONJG(TRANSPOSE(udv2%U))
        CALL INV(udv1%V, V1INV,Z)
#if defined(STAB3) || defined(STABLOG)
#if defined(STABLOG)
        DO J=1,LQ
          !keep scales smaller than 1.0 in D1*U1 and D2*V2
          !bring scales larger that 1.0 with V1^-1 and U2^-1
          If(udv1%L(J) <= 0.d0) then
            D1m(J)=cmplx(exp(udv1%L(J)),0.d0,kind(0.d0))
          else
            D1m(J)=cmplx(1.d0,0.d0,kind(0.d0))
            call zscal(LQ,cmplx(exp(-udv1%L(J)),0.d0,kind(0.d0)),V1INV(J,1),LQ)
          endif
          If(udv2%L(J) <= 0.d0) then
            D2m(J)=cmplx(exp(udv2%L(J)),0.d0,kind(0.d0))
          else
            D2m(J)=cmplx(1.d0,0.d0,kind(0.d0))
            call zscal(LQ,cmplx(exp(-udv2%L(J)),0.d0,kind(0.d0)),MYU2(J,1),LQ)
          endif
        ENDDO
#else
        DO J=1,LQ
          !keep scales smaller than 1.0 in D1*U1 and D2*V2
          !bring scales larger that 1.0 with V1^-1 and U2^-1
          If(dble(udv1%D(J)) <= 1.d0) then
            D1m(J)=udv1%D(J)
          else
            D1m(J)=cmplx(1.d0,0.d0,kind(0.d0))
            call zscal(LQ,1.d0/udv1%D(J),V1INV(J,1),LQ)
          endif
          If(dble(udv2%D(J)) <=1.d0) then
            D2m(J)=udv2%D(J)
          else
            D2m(J)=cmplx(1.d0,0.d0,kind(0.d0))
            call zscal(LQ,1.d0/udv2%D(J),MYU2(J,1),LQ)
          endif
        ENDDO
#endif
#else
        D1m=udv1%D
        D2m=udv2%D
#endif
#if defined(STABLOG)
        If (udv1%L(1) >  udv2%L(1) ) Then
#else 
        If (dble(udv1%D(1)) >  dble(udv2%D(1)) ) Then
#endif
           !Write(6,*) "D1(1) >  D2(1)", dble(D1(1)), dble(D2(1))
           call zlacpy('A',LQ,LQ,V1INV(1,1) ,LQ,HLPB2(1   ,1   ),LQ2)
           call zlacpy('A',LQ,LQ,MYU2(1,1)  ,LQ,HLPB2(1+LQ,1+LQ),LQ2)
           DO J = 1,LQ
              DO I = 1,LQ
                 HLPB2(I   , J+LQ ) =  D1m(I)*conjg(udv1%U(J,I))!udv1%U(I,J)
                 HLPB2(I+LQ, J    ) = -D2m(I)*udv2%V(I,J)
              ENDDO
           ENDDO
           ! HLPB1 = CT(HLPB2)
           HLPB1 = conjg(transpose(HLPB2))
           call QDRP_decompose(LQ2, LQ2, HLPB1, D3, IPVT, TAU, WORK, LWORK)
           call solve_extended_System(HLPB2, V1INV, MYU2, HLPB1, D3, TAU, IPVT, LQ, WORK, LWORK)
           call get_blocks(GR00, GR0T, GRT0, GRTT, HLPB2, LQ)
        Else
           !Write(6,*) "D1(1) <  D2(1)", dble(D1(1)), dble(D2(1))
           call zlacpy('A',LQ,LQ,MYU2(1,1) ,LQ,HLPB2(1   ,1   ),LQ2)
           call zlacpy('A',LQ,LQ,V1INV(1,1),LQ,HLPB2(1+LQ,1+LQ),LQ2)
           DO J = 1,LQ
              DO I = 1,LQ
                 HLPB2(I   , J+LQ ) = -D2m(I)*udv2%V(I,J)
                 HLPB2(I+LQ, J    ) =  D1m(I)*conjg(udv1%U(J,I))!udv1%U(I,J)
              ENDDO
           ENDDO
           ! HLPB1 = CT(HLPB2)
           HLPB1 = conjg(transpose(HLPB2))
           call QDRP_decompose(LQ2, LQ2, HLPB1, D3, IPVT, TAU, WORK, LWORK)
           call solve_extended_System(HLPB2, MYU2, V1INV, HLPB1, D3, TAU, IPVT, LQ, WORK, LWORK)
           call get_blocks(GRTT, GRT0, GR0T, GR00, HLPB2, LQ)
        Endif
        DEALLOCATE(MYU2, V1INV, HLPB1, HLPB2, WORK, IPVT, TAU, D3, D1m, D2m)
#endif
      END SUBROUTINE CGR2_2

end module cgr2_2_mod
