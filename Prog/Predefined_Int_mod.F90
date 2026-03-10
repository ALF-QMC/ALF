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
!> Predefined interaction operators for auxiliary-field QMC.
!>
!> Each routine in this module configures a single \c Operator instance
!> that represents the Hubbard-Stratonovich (HS) decoupling of a quartic
!> fermionic interaction.  The common pattern is:
!>
!> \verbatim
!>   Op%P(k)    = global site index of the k-th active degree of freedom
!>   Op%O(i,j)  = single-particle matrix A of the bilinear  c† A c
!>   Op%g       = HS coupling, absorbs Δτ and the physical coupling constant
!>   Op%alpha   = constant shift inside the bilinear: c† A c + alpha
!>   Op%type    = HS type (1=Ising ±1; 2=discrete 4-field; 3=continuous scalar)
!> \endverbatim
!>
!> After filling these fields the routine calls \c Op_set, which
!> pre-computes the exponentials used during the QMC sweep.
!>
!> @see Operator_mod for the definition of the \c Operator type.
!
!--------------------------------------------------------------------

    Module Predefined_Int
      
      Use Operator_mod
      
      Implicit none
      
      
    contains
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Hubbard U interaction decoupled with a discrete 4-field HS transformation, SU(N) flavour symmetry.
!> Represents \f$ \frac{U}{N_{\mathrm{SUN}}} \left[\sum_{\sigma=1}^{N_{\mathrm{SUN}}}\left(n_{i,\sigma}-\frac{1}{2}\right)\right]^2 \f$.
!> One \c Operator per site is sufficient because all N_SUN flavours couple to the same HS field.
!>
!> @param[out] OP     Operator instance configured for site \p I.
!> @param[in]  I      Global site index of the interaction site.
!> @param[in]  N_SUN  Number of SU(N) flavours (= number of fermion flavours N_FL).
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  U      Hubbard coupling constant.  \f$U > 0\f$ is repulsive.
!--------------------------------------------------------------------
      Subroutine Predefined_Int_U_SUN( OP, I, N_SUN, DTAU, U  )
        
        Implicit none
        Integer,  Intent(In) :: I, N_SUN 
        Real (Kind=Kind(0.d0)), Intent(IN) :: Dtau, U
        Type(Operator), Intent(Out) :: OP

        Call OP_Make( Op,1 )

        ! One active site; O = 1 acts on (n - 1/2); alpha = -1/2 centres the density.
        ! g = sqrt(-Dtau*U/N_SUN): from (n-1/2)^2 summed over N_SUN flavours the HS
        ! coupling picks up a factor 1/N_SUN.  The imaginary unit from sqrt(-x) is
        ! absorbed into the complex coupling g so that exp(g*phi*(n-1/2)) is real
        ! for phi in {±1, ±2} (discrete 4-field HS, type=2).
        Op%P(1)   = I
        Op%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
        Op%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
        Op%g      = SQRT(CMPLX(-DTAU*U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
        Op%type   = 2

        Call Op_set( Op )

      end Subroutine Predefined_Int_U_SUN

!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Hubbard U interaction decoupled with a continuouos HS transformation, SU(N) flavour symmetry.
!> Same interaction as Predefined_Int_U_SUN but uses the continuous Gaussian identity
!> \f$ e^{A^2} = \int_{-\infty}^{\infty} \frac{dx}{\sqrt{2\pi}} e^{-x^2/2 + \sqrt{2}A\,x} \f$.
!> The extra factor of 2 in \c g relative to the discrete variant comes from \f$\sqrt{2}\f$ in the
!> Hubbard-Stratonovich identity.
!>
!> @param[out] OP     Operator instance configured for site \p I.
!> @param[in]  I      Global site index of the interaction site.
!> @param[in]  N_SUN  Number of SU(N) flavours.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  U      Hubbard coupling constant.
!> @see Predefined_Int_U_SUN for the discrete-HS variant.
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
!>
!> @brief
!> Hubbard U in the \f$M_z\f$ channel using a discrete 4-field HS transformation.
!> Decouples \f$ -\frac{U}{2}\left[(n_{i,\uparrow}-\tfrac{1}{2})-(n_{i,\downarrow}-\tfrac{1}{2})\right]^2 \f$.
!> Requires N_FL = 2.  Sets both the spin-up and spin-down operators with
!> opposite signs of \c g so they couple to the same HS field with opposite weights.
!>
!> @param[out] OP_up  Operator instance for spin-up electrons at site \p I.
!> @param[out] Op_do  Operator instance for spin-down electrons at site \p I.
!> @param[in]  I      Global site index.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  U      Hubbard coupling constant.
!> @see Predefined_Int_U_MZ_continuous_HS for the continuous-HS variant.
!> @see Predefined_Int_U_SUN for the SU(N)-symmetric single-operator variant.
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
        Op_up%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))   ! placeholder — overwritten below
        Op_up%g      =  SQRT(CMPLX(DTAU*U/2.d0, 0.D0, kind(0.D0)))
        Op_up%alpha  =  cmplx(-0.5d0,0.d0, kind(0.D0))  ! final value: centres density at 1/2
        Op_up%type   =  2

        ! Spin-down couples to the same HS field with g_do = -g_up, which ensures
        ! that the product exp(g_up*phi*(n_up-1/2)) * exp(g_do*phi*(n_do-1/2)) =
        ! exp(g_up*phi*[(n_up-1/2)-(n_do-1/2)]) decouples the desired Mz^2 channel.
        Op_do%P(1)   =  I
        Op_do%O(1,1) =  cmplx(1.d0  ,0.d0, kind(0.D0))
        Op_do%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))   ! placeholder — overwritten below
        Op_do%g      = -SQRT(CMPLX(DTAU*U/2.d0, 0.D0, kind(0.D0)))
        Op_do%alpha  =  cmplx(-0.5d0,0.d0, kind(0.D0))  ! final value
        Op_do%type   =  2

        Call Op_set( Op_up )
        Call Op_set( Op_do )
        
      end Subroutine Predefined_Int_U_MZ
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Hubbard U in the \f$M_z\f$ channel using a continuous HS transformation.
!> Same physical interaction as Predefined_Int_U_MZ.
!> Uses the Gaussian identity \f$ e^{A^2} = \int dx\, e^{-x^2/2+\sqrt{2}Ax}/\sqrt{2\pi} \f$
!> with \f$A^2 = \Delta\tau (U/2)(n_{\uparrow}-n_{\downarrow})^2 \f$.
!>
!> @param[out] OP_up  Operator instance for spin-up electrons at site \p I.
!> @param[out] Op_do  Operator instance for spin-down electrons at site \p I.
!> @param[in]  I      Global site index.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  U      Hubbard coupling constant.
!> @see Predefined_Int_U_MZ for the discrete-HS variant.
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
!>
!> @brief
!> Density-density (bond-kinetic) interaction decoupled via a discrete 4-field HS, SU(N) symmetric.
!> Represents \f$ -\frac{V}{N_{\mathrm{SUN}}}\left[\sum_{s=1}^{N_{\mathrm{SUN}}}\left(c^\dagger_{i,s}c_{j,s}+c^\dagger_{j,s}c_{i,s}\right)\right]^2 \f$.
!> The two-site operator matrix has off-diagonal entries 1 to encode the hopping bilinear.
!>
!> @param[out] OP     Operator instance for the (i,j) bond.
!> @param[in]  I      Global site index of the first site.
!> @param[in]  J      Global site index of the second site.
!> @param[in]  N_SUN  Number of SU(N) flavours.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  V      Density-density coupling constant.
!> @see Predefined_Int_VJ_SUN for the imaginary (current-current) variant.
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
!>
!> @brief
!> Current-current (imaginary hopping) interaction, SU(N) symmetric.
!> Represents \f$ -\frac{V}{N_{\mathrm{SUN}}}\left[\sum_{s=1}^{N_{\mathrm{SUN}}}\left(ic^\dagger_{i,s}c_{j,s}-ic^\dagger_{j,s}c_{i,s}\right)\right]^2 \f$.
!> The operator matrix is anti-Hermitian: \f$O_{12}=+i,\; O_{21}=-i\f$, encoding the
!> current (imaginary hopping) on the bond.
!>
!> @param[out] OP     Operator instance for the (i,j) bond.
!> @param[in]  I      Global site index of the first site.
!> @param[in]  J      Global site index of the second site.
!> @param[in]  N_SUN  Number of SU(N) flavours.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  V      Current-current coupling constant.
!> @see Predefined_Int_V_SUN for the real (density-density) variant.
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
!>
!> @brief
!> Ising-bond interaction with an Ising (type=1) HS field, SU(N) symmetric.
!> Represents the bond operator
!> \f$ \hat{Z}_{i,j}\,\xi\,\sum_{s=1}^{N_{\mathrm{SUN}}}\left(c^\dagger_{i,s}c_{j,s}+c^\dagger_{j,s}c_{i,s}\right) \f$
!> where \f$\hat{Z}_{i,j}\in\{\pm 1\}\f$ is the Ising auxiliary field.
!> The coupling \c g = -Dtau*Xi is real; the Ising field type=1 restricts phi to \f$\{+1,-1\}\f$.
!>
!> @param[out] OP     Operator instance for the (i,j) Ising bond.
!> @param[in]  I      Global site index of the first site.
!> @param[in]  J      Global site index of the second site.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  Xi     Ising coupling strength \f$\xi\f$.
!> @see Predefined_Int_V_SUN for a quartic density-density decoupling.
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
!>
!> @brief
!> Long-range Coulomb (LRC) vertex.
!> Couples each flavour to a continuous bosonic field \f$\Phi\f$ via
!> \f$ \sum_{s=1}^{N_{\mathrm{SUN}}} i\Phi\,\left(c^\dagger_{i,s}c_{i,s}-\tfrac{1}{2}\right) \f$.
!> The purely imaginary coupling \f$g = i\,\Delta\tau\f$ means that the product
!> of all LRC vertices in a configuration gives a real Boltzmann weight after
!> integrating out the fermions.
!>
!> @param[out] OP     Operator instance for site \p I.
!> @param[in]  I      Global site index.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
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
!>
!> @brief
!> Ising (Jz) spin–spin interaction via a discrete 4-field HS.
!> Decouples
!> \f$ -\frac{|J_z|}{2}\left(S^z_i - \operatorname{sgn}(J_z)S^z_j\right)^2
!>   = J_z S^z_i S^z_j - \frac{|J_z|}{2}\left[(S^z_i)^2+(S^z_j)^2\right] \f$.
!> Requires N_FL = 2 (spin-up and spin-down flavours).  The sign of Jz controls
!> ferro (Jz > 0) vs antiferro (Jz < 0) coupling via the \f$O_{22}\f$ entry.
!>
!> @param[out] OP_up  Operator for spin-up at sites I and J.
!> @param[out] Op_do  Operator for spin-down at sites I and J.
!> @param[in]  I      Global site index of the first site.
!> @param[in]  J      Global site index of the second site.
!> @param[in]  DTAU   Imaginary-time step \f$\Delta\tau\f$.
!> @param[in]  Jz     Ising coupling.  \f$J_z > 0\f$ → ferromagnetic; \f$J_z < 0\f$ → antiferromagnetic.
!> @see Predefined_Int_U_MZ for charging-energy decoupling in the Mz channel.
!--------------------------------------------------------------------
      Subroutine Predefined_Int_Jz (  OP_up,  Op_do, I, J,  DTAU, Jz  ) 

        Implicit none
        Integer,  Intent(In) :: I,J
        Real (Kind=Kind(0.d0)), Intent(IN) ::  Dtau, Jz
        Type(Operator), Intent(Out) :: OP_up, Op_do

        Call OP_Make( Op_up,2 )
        Call OP_Make( Op_do,2 )
        
        ! O = diag(1, -sgn(Jz)): the sign of Jz encodes ferro (Jz>0) vs antiferro (Jz<0)
        ! coupling.  alpha = 0 since the decomposition uses (n_up - n_do) = 2*S^z.
        ! g = sqrt(Dtau*|Jz|/8): the factor 1/8 arises from expanding
        !   (S^z_i - sgn(Jz)*S^z_j)^2 = ((n_up-n_do)_i/2 - sgn(Jz)*(n_up-n_do)_j/2)^2
        ! and collecting the 4 cross terms, each contributing a factor of 1/2 (S^z = (n_up-n_do)/2).
        ! Spin-up and spin-down operators share the same g magnitude but opposite sign
        ! so they couple to the same HS field and reproduce the desired Sz*Sz product.
        Op_up%P(1)   =  I
        Op_up%P(2)   =  J
        Op_up%O(1,1) =  cmplx(        1.d0  ,0.d0, kind(0.D0))
        Op_up%O(2,2) =  cmplx(- Jz/Abs(Jz)  ,0.d0, kind(0.D0))  ! -sgn(Jz): ferro vs antiferro
        Op_up%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_up%g      =  SQRT(CMPLX(DTAU*Abs(Jz)/8.d0, 0.D0, kind(0.D0)))
        Op_up%type   =  2

        Op_do%P(1)   =  I
        Op_do%P(2)   =  J
        Op_do%O(1,1) =  cmplx(        1.d0  ,0.d0, kind(0.D0))
        Op_do%O(2,2) =  cmplx(- Jz/Abs(Jz)  ,0.d0, kind(0.D0))  ! same sign as Op_up
        Op_do%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
        Op_do%g      = -SQRT(CMPLX(DTAU*Abs(Jz)/8.d0, 0.D0, kind(0.D0)))  ! opposite sign
        Op_do%type   =  2

        Call Op_set( Op_up ) 
        Call Op_set( Op_do )
        
      end Subroutine Predefined_Int_Jz

!!$       Still has to be done when implementing Hamiltoninas.
!!$      Subroutine Predefined_Int_J_SUN( )
!!$        Implicit none
!!$      end Subroutine Predefined_Int_J_SUN


     end Module Predefined_Int
