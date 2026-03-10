!  Copyright (C) 2016 - 2023 The ALF project
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
!> Predefined observable measurement routines for auxiliary-field QMC.
!>
!> Provides equal-time and imaginary-time-displaced correlation functions
!> (Green, density-density, spin-spin) as well as second Rényi entropy and
!> mutual information estimators.  Each routine accumulates a symmetry-resolved
!> structure factor into an \c Obser_Latt accumulator.
!>
!> **Indexing convention** (common to all measurement routines):
!> \verbatim
!>   I1, J1   = composite site index running over Size(List,1)
!>   I = List(I1,1) = unit-cell index    (1..Latt%%N)
!>   no_I = List(I1,2) = orbital index   (1..Latt_Unit%%Norb)
!>   imj  = latt%%imj(I,J) = translation-class index used for the
!>          momentum-space structure factor accumulation
!>   GRC(I1,J1,nf) = <c_{I1} c^dag_{J1}> for fermion flavour nf
!>   GR (I1,J1,nf) = <c^dag_{I1} c_{J1}> = 1 - GRC (at equal time)
!>   ZP  = Phase / Re(Phase)  -- pure phase factor (unit modulus) derived from the
!>          complex configuration weight; used to reweight observables under the
!>          assumption that the partition function is real.
!>   ZS  = sign( Re(Phase) ) * Mc_step_weight
!>          -- sign of the real part of the weight, multiplied by an additional
!>          MC step weight.  In standard Metropolis sampling Mc_step_weight = 1;
!>          for Langevin-type updates it encodes a step-size-dependent correction
!>          factor.  ZS is accumulated into the average sign.
!>   Physical estimator:  <O> = <O * ZP * ZS>_{|w|} / <ZS>_{|w|}
!>   Obs%%Obs_Latt(imj,it,no_I,no_J) += quantity * ZP * ZS
!> \endverbatim
!>
!> For time-displaced routines the four propagators are:
!> \verbatim
!>   GT0(I1,J1,nf) = G(tau, 0) = <c_{I1}(tau) c^dag_{J1}(0)>
!>   G0T(I1,J1,nf) = G(0, tau) = <c_{I1}(0)   c^dag_{J1}(tau)>
!>   G00(I1,J1,nf) = G(0, 0)
!>   GTT(I1,J1,nf) = G(tau, tau)
!> \endverbatim
!>   NT is the imaginary-time index (0..Ltau-1); counts and sign are
!>   accumulated only when NT == 0 to avoid double-counting.
!>
!> @see Operator_mod for the Operator type used in Predefined_Obs_V_Int.
!
!--------------------------------------------------------------------


    Module Predefined_Obs

      use runtime_error_mod
      use Operator_mod
      Use Observables
      Use Lattices_v3
      Use entanglement_mod
      use iso_fortran_env, only: output_unit, error_unit

      Implicit none
      
      INTERFACE Predefined_Obs_scal_Renyi_Ent
        MODULE PROCEDURE Predefined_Obs_scal_Renyi_Ent_gen_all, Predefined_Obs_scal_Renyi_Ent_indep, &
        & Predefined_Obs_scal_Renyi_Ent_gen_fl
      END INTERFACE
      INTERFACE Predefined_Obs_scal_Mutual_Inf
        MODULE PROCEDURE Predefined_Obs_scal_Mutual_Inf_indep, Predefined_Obs_scal_Mutual_Inf_gen_fl, &
        & Predefined_Obs_scal_Mutual_Inf_gen_all
      END INTERFACE

    contains
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Measure SU(N) spin-spin correlations (equal time).
!> When N_FL = 1 returns:
!> \f$ \frac{2N}{N^2-1}\sum_{a=1}^{N^2-1}\langle c^\dagger_i T^a c_i\; c^\dagger_j T^a c_j\rangle \f$
!> where \f$T^a\f$ are the SU(N) generators normalised as \f$\mathrm{Tr}[T^a T^b]=\delta_{ab}/2\f$.
!> Using SU(N) symmetry this reduces to \f$N_{\mathrm{SUN}}\cdot G_{\rm RC}(i,j)G_{\rm R}(i,j)\f$.
!>
!> @param[in]  Latt       Lattice geometry (used for imj table).
!> @param[in]  Latt_unit  Unit-cell descriptor (Norb etc.).
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  GR         Equal-time propagator \f$\langle c^\dagger c\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  GRC        Equal-time propagator \f$\langle c c^\dagger\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  N_SUN      Number of SU(N) flavours.
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs      Accumulator (File_Latt must equal "SpinZ").
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
        Complex (Kind=Kind(0.d0)) :: ZZ


        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2  )   then
           Write(error_unit,*) 'List in  Predefined_Obs_eq_SpinSUN_measure has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        If ( Obs%File_Latt .ne. "SpinZ" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_SpinSUN_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))

        ! Measure
        N_FL = Size(GR,3)
        If (N_FL == 1)   then
           Do I1 = 1,Size(List,1)
              I    = List(I1,1)
              If  ( I > 0 ) then 
                 no_I = List(I1,2)
                 Do J1 = 1,Size(List,1)
                    J    = List(J1,1)
                    if  (  J  > 0 )  then 
                       no_J = List(J1,2)
                       imj  = latt%imj(I,J)
                       ZZ   = GRC(I1,J1,1) * GR(I1,J1,1) * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                       Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + ZZ*ZP*ZS
                    endif
                 enddo
              endif
           enddo
        endif


      end Subroutine Predefined_Obs_eq_SpinSUN_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Measure SU(2) spin-spin correlations in the Mz channel (equal time).
!> For N_FL = 2, N_SUN = 1 returns three observables:
!> \f$ S_z       = 4\langle S^z_i S^z_j\rangle \f$,
!> \f$ S_{xy}    = 2(\langle S^x_i S^x_j\rangle + \langle S^y_i S^y_j\rangle) \f$,
!> \f$ S_{\rm T} = (2S_{xy}+S_z)/3 \f$ (isotropic total spin correlator).
!>
!> @param[in]  Latt       Lattice geometry.
!> @param[in]  Latt_unit  Unit-cell descriptor.
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  GR         Equal-time propagator \f$\langle c^\dagger c\rangle\f$, shape (Ndim, Ndim, 2).
!> @param[in]  GRC        Equal-time propagator \f$\langle c c^\dagger\rangle\f$, shape (Ndim, Ndim, 2).
!> @param[in]  N_SUN      Number of SU(N) colour flavours (must be 1 for this routine).
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] ObsZ     Accumulator for \f$S_z\f$ (File_Latt = "SpinZ").
!> @param[inout] ObsXY    Accumulator for \f$S_{xy}\f$ (File_Latt = "SpinXY").
!> @param[inout] ObsXYZ   Accumulator for \f$S_{\rm T}\f$ (File_Latt = "SpinT").
!> @see Predefined_Obs_eq_SpinSUN_measure for the SU(N) variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_eq_SpinMz_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, ObsZ, ObsXY, ObsXYZ )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: ObsZ, ObsXY, ObsXYZ

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
        Complex (Kind=Kind(0.d0)) :: ZXY, ZZ

        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in Predefined_Obs_eq_SpinMz_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        
        If ( ObsZ%File_Latt .ne. "SpinZ" .and. ObsXY%File_Latt .ne. "SpinXY" .and.  &
           & ObsXYZ%File_Latt .ne. "SpinT"  )   then
           Write(error_unit,*) 'Predefined_Obs_eq_SpinMz_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        ObsZ%N          = ObsZ%N + 1
        ObsZ%Ave_sign   = ObsZ%Ave_sign + real(ZS,kind(0.d0))
        ObsXY%N         = ObsXY%N + 1
        ObsXY%Ave_sign  = ObsXY%Ave_sign + real(ZS,kind(0.d0))
        ObsXYZ%N        = ObsXYZ%N + 1
        ObsXYZ%Ave_sign = ObsXYZ%Ave_sign + real(ZS,kind(0.d0))


        ! Measure
        N_FL = Size(GR,3)
        If (N_FL == 2 .and. N_SUN == 1 ) Then
           Do I1 = 1,Size(List,1) 
              I    = List(I1,1)
              If (I  > 0 )  then 
                 no_I = List(I1,2)
                 Do J1 = 1,Size(List,1) 
                    J    = List(J1,1)
                    If (J > 0 )  then 
                       no_J = List(J1,2)
                       imj  = latt%imj(I,J)
                       ZXY  = GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1)
                       ZZ   = GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                            (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))
                       
                       ObsZ  %Obs_Latt(imj,1,no_I,no_J) =  ObsZ  %Obs_Latt(imj,1,no_I,no_J) +  ZZ  *ZP*ZS
                       ObsXY %Obs_Latt(imj,1,no_I,no_J) =  ObsXY %Obs_Latt(imj,1,no_I,no_J) +  ZXY *ZP*ZS
                       ObsXYZ%Obs_Latt(imj,1,no_I,no_J) =  ObsXYZ%Obs_Latt(imj,1,no_I,no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                    endif
                 enddo
                 ObsZ  %Obs_Latt0(no_I) =  ObsZ  %Obs_Latt0(no_I) +   (GRC(I1,I1,2) - GRC(I1,I1,1)) * ZP*ZS
                 ObsXYZ%Obs_Latt0(no_I) =  ObsXYZ%Obs_Latt0(no_I) +   (GRC(I1,I1,2) - GRC(I1,I1,1)) * ZP*ZS/sqrt(3.d0)
              endif
           enddo
        endif
      End Subroutine Predefined_Obs_eq_SpinMz_measure


!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Measure equal-time single-particle Green function, summed over flavours:
!> \f$ \sum_{\sigma=1}^{N_{\rm SUN}}\sum_{nf=1}^{N_{\rm FL}}\langle c^\dagger_{i,\sigma,nf} c_{j,\sigma,nf}\rangle \f$.
!>
!> @param[in]  Latt       Lattice geometry.
!> @param[in]  Latt_unit  Unit-cell descriptor.
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  GR         Equal-time propagator \f$\langle c^\dagger c\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  GRC        Equal-time propagator \f$\langle c c^\dagger\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  N_SUN      Number of SU(N) flavours (overall prefactor).
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs      Accumulator (File_Latt must equal "Green").
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_eq_Green_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z

        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in  Predefined_Obs_eq_Green_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        If ( Obs%File_Latt .ne. "Green" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_Green_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))

        ! Measure
        N_FL = Size(GR,3)
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If ( I  > 0 ) then 
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  ( J > 0 ) then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    Z = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1, N_FL
                       Z = Z +  GRC(I1,J1,nf)
                    enddo
                    Z   = Z * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + Z*ZP*ZS
                 endif
              enddo
           endif
        enddo

      end Subroutine Predefined_Obs_eq_Green_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Measure equal-time density-density connected correlator:
!> \f$ \langle N_i N_j\rangle - \langle N_i\rangle\langle N_j\rangle \f$
!> where \f$N_i = \sum_{\sigma,nf} c^\dagger_{i,\sigma,nf} c_{i,\sigma,nf}\f$.
!>
!> @param[in]  Latt       Lattice geometry.
!> @param[in]  Latt_unit  Unit-cell descriptor.
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  GR         Equal-time propagator \f$\langle c^\dagger c\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  GRC        Equal-time propagator \f$\langle c c^\dagger\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  N_SUN      Number of SU(N) colour flavours.
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs      Accumulator (File_Latt must equal "Den").
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_eq_Den_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: ZI, ZJ, Z


        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in  Predefined_Obs_eq_Den_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        If ( Obs%File_Latt .ne. "Den" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_Den_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))

        ! Measure

        N_FL = Size(GR,3)
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If (I > 0 )  then 
              no_I = List(I1,2)
              ZI = cmplx(0.d0,0.d0,Kind(0.d0))
              Do nf = 1,N_FL
                 ZI  = ZI + GRC(I1,I1,nf)
              enddo
              ZI  = ZI * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  (J > 0 ) then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    ZJ = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       ZJ  = ZJ + GRC(J1,J1,nf)
                    enddo
                    ZJ  = ZJ * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Z   = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       Z = Z + GRC(I1,J1,nf) * GR(I1,J1,nf)
                    Enddo
                    Z   =  Z * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + (ZI*ZJ + Z)*ZP*ZS
                 endif
              enddo
              Obs%Obs_Latt0(no_I) =  Obs%Obs_Latt0(no_I) +  ZI * ZP*ZS
           endif
        enddo
      end Subroutine Predefined_Obs_eq_Den_measure


!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Measure time-displaced single-particle Green function:
!> \f$ \sum_{\sigma,nf}\langle c^\dagger_{i,\sigma,nf}(\tau)\, c_{j,\sigma,nf}(0)\rangle \f$.
!> Accumulates into \c Obs%%Obs_Latt(:, NT+1, :, :).
!> Counts and sign are updated only at NT = 0.
!>
!> @param[in]  Latt       Lattice geometry.
!> @param[in]  Latt_unit  Unit-cell descriptor.
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  NT         Imaginary-time slice index (0..Ltau-1).
!> @param[in]  GT0        \f$G(\tau,0)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  G0T        \f$G(0,\tau)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  G00        \f$G(0,0)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  GTT        \f$G(\tau,\tau)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  N_SUN      Number of SU(N) flavours.
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs      Accumulator (File_Latt = "Green").
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_tau_Green_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z


        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in  Predefined_Obs_tau_Green_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        If ( Obs%File_Latt .ne. "Green" )   then
           Write(error_unit,*) 'Predefined_Obs_tau_Green_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        ! Count and average sign
        If (NT == 0 ) then
           Obs%N        = Obs%N + 1
           Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        endif

        ! Measure
        N_FL = Size(GT0,3)
        Do I1 = 1,Size(List,1) 
           I    = List(I1,1)
           If  ( I > 0 ) then 
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1) 
                 J    = List(J1,1)
                 If  (J > 0 )  then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    Z = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1, N_FL
                       Z = Z +  GT0(I1,J1,nf)
                    enddo
                    Z   = Z * cmplx(1.d0/dble(N_FL), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs%Obs_Latt(imj,NT+1,no_I,no_J) + Z*ZP*ZS
                 endif
              enddo
           endif
        enddo

      end Subroutine Predefined_Obs_tau_Green_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Measure time-displaced SU(N) spin-spin correlations (N_FL = 1).
!> Returns
!> \f$ \frac{2N}{N^2-1}\sum_a \langle c^\dagger_i(\tau)T^a c_i(\tau)\; c^\dagger_j(0)T^a c_j(0)\rangle \f$.
!>
!> @param[in]  Latt       Lattice geometry.
!> @param[in]  Latt_unit  Unit-cell descriptor.
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  NT         Imaginary-time slice index (0..Ltau-1).
!> @param[in]  GT0        \f$G(\tau,0)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  G0T        \f$G(0,\tau)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  G00        \f$G(0,0)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  GTT        \f$G(\tau,\tau)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  N_SUN      Number of SU(N) flavours (N_FL must equal 1).
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs      Accumulator (File_Latt = "SpinZ").
!> @see Predefined_Obs_tau_SpinMz_measure for the N_FL = 2 spin-1/2 variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_tau_SpinSUN_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z

        If ( Obs%File_Latt .ne. "SpinZ" )   then
           Write(error_unit,*) 'Predefined_Obs_tau_SpinSUN_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2  )   then
           Write(error_unit,*) 'List in  Predefined_Obs_tau_SpinSUN_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        ! Count and average sign
        If (NT == 0 ) then
           Obs%N        = Obs%N + 1
           Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        endif

        ! Measure
        N_FL = Size(GT0,3)
        Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
        Do I1 = 1,Size(List,1)
           I   = List(I1,1)
           If (I > 0 ) then
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  (J > 0)  then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    Obs%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs%Obs_Latt(imj,NT+1,no_I,no_J) &
                         &  -  Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                 endif
              enddo
           endif
        enddo

      end Subroutine Predefined_Obs_tau_SpinSUN_measure


!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Measure time-displaced SU(2) spin correlations for Mz models (N_FL = 2, N_SUN = 1).
!> Returns three time-displaced structure factors:
!> \f$ S_z(\tau)    = 4\langle S^z_i(\tau)S^z_j(0)\rangle \f$,
!> \f$ S_{xy}(\tau) = 2(\langle S^x_i(\tau)S^x_j\rangle+\langle S^y_i(\tau)S^y_j\rangle) \f$,
!> \f$ S_{\rm T}(\tau) = (2S_{xy}+S_z)/3 \f$.
!>
!> @param[in]  Latt       Lattice geometry.
!> @param[in]  Latt_unit  Unit-cell descriptor.
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  NT         Imaginary-time slice index (0..Ltau-1).
!> @param[in]  GT0        \f$G(\tau,0)\f$, shape (Ndim, Ndim, 2).
!> @param[in]  G0T        \f$G(0,\tau)\f$, shape (Ndim, Ndim, 2).
!> @param[in]  G00        \f$G(0,0)\f$, shape (Ndim, Ndim, 2).
!> @param[in]  GTT        \f$G(\tau,\tau)\f$, shape (Ndim, Ndim, 2).
!> @param[in]  N_SUN      Number of SU(N) colour flavours (must be 1).
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] ObsZ     Accumulator for \f$S_z(\tau)\f$.
!> @param[inout] ObsXY    Accumulator for \f$S_{xy}(\tau)\f$.
!> @param[inout] ObsXYZ   Accumulator for \f$S_{\rm T}(\tau)\f$.
!> @see Predefined_Obs_tau_SpinSUN_measure for the SU(N) variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_tau_SpinMz_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, ObsZ, ObsXY, ObsXYZ )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: ObsZ, ObsXY, ObsXYZ

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: ZZ, ZXY

        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in Predefined_Obs_tau_SpinMz_measure has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif


        ! Count and average sign
        If ( ObsZ%File_Latt .ne. "SpinZ" .and. ObsXY%File_Latt .ne. "SpinXY" .and.  &
           & ObsXYZ%File_Latt .ne. "SpinT"  )   then
           Write(error_unit,*) 'Predefined_Obs_tau_SpinMz_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        If (NT == 0 ) then
           ! Count and average sign
           ObsZ%N          = ObsZ%N + 1
           ObsZ%Ave_sign   = ObsZ%Ave_sign + real(ZS,kind(0.d0))
           ObsXY%N         = ObsXY%N + 1
           ObsXY%Ave_sign  = ObsXY%Ave_sign + real(ZS,kind(0.d0))
           ObsXYZ%N        = ObsXYZ%N + 1
           ObsXYZ%Ave_sign = ObsXYZ%Ave_sign + real(ZS,kind(0.d0))
        endif

        ! Measure
        N_FL = Size(GT0,3)
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If  (I > 0 ) then
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 if (J > 0 ) then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    ZZ   =      (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                         &    -  G0T(J1,I1,1) * GT0(I1,J1,1)  -  G0T(J1,I1,2) * GT0(I1,J1,2) 
                    ZXY  =    -  G0T(J1,I1,1) * GT0(I1,J1,2)  -  G0T(J1,I1,2) * GT0(I1,J1,1) 
                    ObsZ  %Obs_Latt(imj,NT+1,no_I,no_J) =  ObsZ %Obs_Latt(imj,NT+1,no_I,no_J)  +  ZZ*ZP*ZS
                    ObsXY %Obs_Latt(imj,NT+1,no_I,no_J) =  ObsXY%Obs_Latt(imj,NT+1,no_I,no_J)  +  ZXY*ZP*ZS
                    ObsXYZ%Obs_Latt(imj,NT+1,no_I,no_J) =  ObsXYZ%Obs_Latt(imj,NT+1,no_I,no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                 endif
              enddo
              ObsZ  %Obs_Latt0(no_I) =  ObsZ  %Obs_Latt0(no_I) +   (GTT(I1,I1,2) - GTT(I1,I1,1)) * ZP*ZS
              ObsXYZ%Obs_Latt0(no_I) =  ObsXYZ%Obs_Latt0(no_I) +   (GTT(I1,I1,2) - GTT(I1,I1,1)) * ZP*ZS/sqrt(3.d0)
           endif
        enddo

      end Subroutine Predefined_Obs_tau_SpinMz_measure


!-------------------------------------------------------------------
!> @Author
!> ALF-project
!>
!> @brief
!> Measure time-displaced density-density connected correlator:
!> \f$ \langle N_i(\tau) N_j(0)\rangle - \langle N_i\rangle\langle N_j\rangle \f$.
!>
!> @param[in]  Latt       Lattice geometry.
!> @param[in]  Latt_unit  Unit-cell descriptor.
!> @param[in]  List       Site-to-(unit-cell, orbital) table, shape (Ndim, 2).
!> @param[in]  NT         Imaginary-time slice index (0..Ltau-1).
!> @param[in]  GT0        \f$G(\tau,0)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  G0T        \f$G(0,\tau)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  G00        \f$G(0,0)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  GTT        \f$G(\tau,\tau)\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  N_SUN      Number of SU(N) colour flavours.
!> @param[in]  ZS         Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP         Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs      Accumulator (File_Latt = "Den").
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_tau_Den_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z, ZI, ZJ
        
        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in Predefined_Obs_tau_Den_measure has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        
        If ( Obs%File_Latt .ne. "Den" )   then
           Write(error_unit,*) 'Predefined_Obs_tau_Den_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        ! Count and average sign
        If (NT == 0 ) then
           Obs%N        = Obs%N + 1
           Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        endif
        N_FL = Size(GT0,3)
        !Measure
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If ( I > 0 )  then 
              no_I = List(I1,2)
              ZI = cmplx(0.d0,0.d0,Kind(0.d0))
              Do nf = 1,N_FL
                 ZI  = ZI  + cmplx(1.d0,0.d0,kind(0.d0)) -  GTT(I1,I1,nf)
              enddo
              ZI  = ZI * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  ( J > 0 ) then
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    ZJ = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       ZJ  = ZJ + cmplx(1.d0,0.d0,kind(0.d0)) -  G00(J1,J1,nf)
                    enddo
                    ZJ  = ZJ * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Z   = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       Z = Z - G0T(J1,I1,nf) * GT0(I1,J1,nf)
                    Enddo
                    Z   =  Z * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs%Obs_Latt(imj,NT+1,no_I,no_J) + (ZI*ZJ + Z)*ZP*ZS
                 endif
              enddo
              Obs%Obs_Latt0(no_I) =  Obs%Obs_Latt0(no_I) +  ZI * ZP*ZS
           endif
        enddo

      end Subroutine Predefined_Obs_tau_Den_measure
      
#include  "Cotunneling_dimer_obs.F90"
!-------------------------------------------------------------------
!> @Author
!> ALF-project
!>
!> @brief
!> Accumulate a type-2 bilinear vertex expectation value:
!> \f$ \left\langle\!\left\langle \left(\sum_{nf,\sigma}\sum_{x,y}
!>   c^\dagger_{x,nf,\sigma} V^{nf}_{x,y} c_{y,nf,\sigma} + \alpha_{nf}\right)^2
!> \right\rangle\!\right\rangle \f$
!> where \f$V^{nf}\f$ is the operator matrix of \c OP_Vint(nf) and the double-bracket
!> denotes a connected expectation value.
!>
!> @param[in]  OP_Vint  Array of type-2 Operator instances (one per fermion flavour).
!>                      OP_Vint(nf)%%O gives \f$V^{nf}\f$; OP_Vint(nf)%%P gives site indices.
!> @param[in]  GR       Equal-time propagator \f$\langle c^\dagger c\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  GRC      Equal-time propagator \f$\langle c c^\dagger\rangle\f$, shape (Ndim, Ndim, N_FL).
!> @param[in]  N_SUN    Number of SU(N) colour flavours (overall prefactor).
!> @return     Complex scalar: the expectation value of the squared bilinear.
!--------------------------------------------------------------------
      Complex  (Kind=Kind(0.d0))  function  Predefined_Obs_V_Int(OP_Vint, GR, GRC, N_SUN )

        Implicit none
        type (Operator)          , Intent(In) :: OP_Vint(:) 
        Integer                  , Intent(In) :: N_SUN
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:)

        ! Local
        Complex  (Kind=Kind(0.d0))  ::  Z,  Z_tmp,  ZC  
        Integer  ::  N_FL, N,  I, J,  I1,  J1, nf
        Real     (Kind=Kind(0.d0))  ::   Zero=1.0D-16

        If ( OP_Vint(1)%type  .ne. 2 )   then
           Write(error_unit,*) 'Predefined_Obs_V_Int  routine is   defined  for  type  2  vertices.  '
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        
        N_FL  =   Size(GR,3)
        N     =   Size(OP_Vint(1)%O,1)

        Z    =  cmplx(0.d0,0.d0,kind(0.d0))
        do  nf  =  1,N_FL 
           Z_tmp     =   cmplx(0.d0,0.d0,kind(0.d0))
           do  J =  1,N 
              do  I  =  1,N
                 Z_tmp  =  Z_tmp  +   OP_Vint(nf)%O(i,j)  * GRC(OP_Vint(nf)%P(I), OP_Vint(nf)%P(J), nf)
              enddo
           enddo
           Z =   Z   +  Z_tmp   +   OP_Vint(nf)%alpha
        enddo
        Z  =  Z*dble(N_SUN)

        ZC  =  cmplx(0.d0,0.d0,kind(0.d0))
        Do nf = 1, N_FL
           Do  J =  1,N
              Do I =  1,N
                 Z_tmp  =   OP_Vint(nf)%O(I,J)
                 If  (  real(Z_tmp*conjg(Z_tmp),  kind(0.d0)) >  Zero )  then
                    Do J1 = 1,N
                       Do I1  = 1,N
                          ZC  =  ZC  +  Z_tmp*OP_Vint(nf)%O(I1,J1)* &
                               &        GRC(OP_Vint(nf)%P(I),OP_Vint(nf)%P(J1),nf)*&
                               &         GR(OP_Vint(nf)%P(J),OP_Vint(nf)%P(I1),nf)
                       Enddo
                    Enddo
                 Endif
              Enddo
           Enddo
        Enddo
        ZC =  ZC * dble(N_SUN)
        Predefined_Obs_V_Int =   Z*Z   + ZC 
        
      end function Predefined_Obs_V_Int

      
      
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Accumulate the second Rényi entropy \f$S_2(A)\f$ of subsystem A
!> assuming independent fermion flavours (N_SUN copies all with the
!> same reduced density matrix).
!> Calls \c Calc_Renyi_Ent_indep from \c entanglement_mod.
!>
!> @param[in]  GRC     Equal-time propagator \f$\langle c c^\dagger\rangle\f$, shape (Ndim_full, Ndim_full, N_FL).
!> @param[in]  List    1-D array of composite site indices defining subsystem A.
!> @param[in]  Nsites  Number of sites in subsystem A.
!> @param[in]  N_SUN   Number of independent SU(N) flavour copies.
!> @param[in]  ZS      Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP      Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs   Scalar accumulator; Obs%%Obs_vec(1) accumulates \f$S_2\f$.
!> @see Predefined_Obs_scal_Renyi_Ent_gen_fl for the per-flavour-general variant.
!> @see Predefined_Obs_scal_Renyi_Ent_gen_all for the fully general variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_scal_Renyi_Ent_indep(GRC, List, Nsites, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:), INTENT(IN)     :: List
        Integer, INTENT(IN)                   :: Nsites, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi

        Renyi = Calc_Renyi_Ent_indep(GRC, List, Nsites, N_SUN)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Renyi_Ent_indep
      
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> General per-flavour Rényi entropy estimator.
!> Allows each fermion flavour to have a different subsystem specification
!> (\c List is 2-D, \c Nsites and \c N_SUN are arrays).
!> Calls \c Calc_Renyi_Ent_gen_fl from \c entanglement_mod.
!>
!> @param[in]  GRC     Equal-time propagator, shape (Ndim, Ndim, N_FL).
!> @param[in]  List    2-D array (Nsites_max, N_FL) of site indices per flavour.
!> @param[in]  Nsites  Array of length N_FL: number of active sites per flavour.
!> @param[in]  N_SUN   Array of length N_FL: degeneracy weight per flavour.
!> @param[in]  ZS      Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP      Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs   Scalar accumulator; Obs%%Obs_vec(1) accumulates \f$S_2\f$.
!> @see Predefined_Obs_scal_Renyi_Ent_indep for the all-flavours-equal variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_scal_Renyi_Ent_gen_fl(GRC, List, Nsites, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:), INTENT(IN)   :: List
        Integer, Dimension(:), INTENT(IN)     :: Nsites, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi

        Renyi = Calc_Renyi_Ent_gen_fl(GRC, List, Nsites, N_SUN)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Renyi_Ent_gen_fl
      
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Fully general Rényi entropy estimator.
!> Handles the most general case: \c List is 3-D (site, index, flavour),
!> \c Nsites is 2-D (index, flavour).  Calls \c Calc_Renyi_Ent_gen_all.
!>
!> @param[in]  GRC     Equal-time propagator, shape (Ndim, Ndim, N_FL).
!> @param[in]  List    3-D array (Nsites_max, ..., N_FL) of site indices.
!> @param[in]  Nsites  2-D array of active site counts, shape (..., N_FL).
!> @param[in]  ZS      Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP      Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs   Scalar accumulator; Obs%%Obs_vec(1) accumulates \f$S_2\f$.
!> @see Predefined_Obs_scal_Renyi_Ent_gen_fl for the per-flavour variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_scal_Renyi_Ent_gen_all(GRC, List, Nsites, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:,:), INTENT(IN)   :: List
        Integer, Dimension(:,:), INTENT(IN)     :: Nsites
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi

        Renyi = Calc_Renyi_Ent_gen_all(GRC, List, Nsites)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Renyi_Ent_gen_all
      
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Accumulate second Rényi entropy for two subsystems A and B, and their
!> union A∪B, to compute the mutual information
!> \f$I_2(A:B) = S_2(A)+S_2(B)-S_2(A\cup B)\f$.
!> Assumes independent flavours (scalar \c N_SUN and 1-D site lists).
!> Calls \c Calc_Mutual_Inf_indep from \c entanglement_mod.
!>
!> @param[in]  GRC       Equal-time propagator, shape (Ndim, Ndim, N_FL).
!> @param[in]  List_A    1-D array of site indices for subsystem A.
!> @param[in]  Nsites_A  Number of sites in A.
!> @param[in]  List_B    1-D array of site indices for subsystem B.
!> @param[in]  Nsites_B  Number of sites in B.
!> @param[in]  N_SUN     Number of independent SU(N) flavour copies.
!> @param[in]  ZS        Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP        Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs     Scalar accumulator (length ≥ 3):
!>                       Obs%%Obs_vec(1) = \f$S_2(A)\f$,
!>                       Obs%%Obs_vec(2) = \f$S_2(B)\f$,
!>                       Obs%%Obs_vec(3) = \f$S_2(A\cup B)\f$.
!> @see Predefined_Obs_scal_Renyi_Ent_indep for single-subsystem variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_scal_Mutual_Inf_indep(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:), INTENT(IN)     :: List_A, List_B
        Integer, INTENT(IN)                   :: Nsites_A ,Nsites_B, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

        call Calc_Mutual_Inf_indep(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, Renyi_A, Renyi_B, Renyi_AB)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi_A* ZP*ZS
        Obs%Obs_vec(2)   = Obs%Obs_vec(2) + Renyi_B* ZP*ZS
        Obs%Obs_vec(3)   = Obs%Obs_vec(3) + Renyi_AB* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Mutual_Inf_indep
      
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Per-flavour-general mutual information estimator.
!> \c List_A and \c List_B are 2-D; \c Nsites_A, \c Nsites_B, and \c N_SUN are arrays.
!> Calls \c Calc_Mutual_Inf_gen_fl from \c entanglement_mod.
!>
!> @param[in]  GRC       Equal-time propagator, shape (Ndim, Ndim, N_FL).
!> @param[in]  List_A    2-D site-index array for subsystem A.
!> @param[in]  Nsites_A  Array of active site counts for A, length N_FL.
!> @param[in]  List_B    2-D site-index array for subsystem B.
!> @param[in]  Nsites_B  Array of active site counts for B, length N_FL.
!> @param[in]  N_SUN     Array of flavour degeneracy weights, length N_FL.
!> @param[in]  ZS        Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP        Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs     Scalar accumulator (length ≥ 3): \f$S_2(A)\f$, \f$S_2(B)\f$, \f$S_2(A\cup B)\f$.
!> @see Predefined_Obs_scal_Mutual_Inf_indep for the scalar-parameter variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_scal_Mutual_Inf_gen_fl(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:), INTENT(IN)   :: List_A, List_B
        Integer, Dimension(:), INTENT(IN)     :: Nsites_A ,Nsites_B, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

        call Calc_Mutual_Inf_gen_fl(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, Renyi_A, Renyi_B, Renyi_AB)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi_A* ZP*ZS
        Obs%Obs_vec(2)   = Obs%Obs_vec(2) + Renyi_B* ZP*ZS
        Obs%Obs_vec(3)   = Obs%Obs_vec(3) + Renyi_AB* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Mutual_Inf_gen_fl
      
!-------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Fully general mutual information estimator.
!> Handles the most general case: \c List_A and \c List_B are 3-D;
!> \c Nsites_A and \c Nsites_B are 2-D.  Calls \c Calc_Mutual_Inf_gen_all.
!>
!> @param[in]  GRC       Equal-time propagator, shape (Ndim, Ndim, N_FL).
!> @param[in]  List_A    3-D site-index array for subsystem A.
!> @param[in]  Nsites_A  2-D array of active site counts for A.
!> @param[in]  List_B    3-D site-index array for subsystem B.
!> @param[in]  Nsites_B  2-D array of active site counts for B.
!> @param[in]  ZS        Sign of Re(Phase) times the MC step weight; accumulated into the average sign.
!> @param[in]  ZP        Phase reweighting factor Phase/Re(Phase); applied to each observable contribution.
!> @param[inout] Obs     Scalar accumulator (length ≥ 3): \f$S_2(A)\f$, \f$S_2(B)\f$, \f$S_2(A\cup B)\f$.
!> @see Predefined_Obs_scal_Mutual_Inf_gen_fl for the per-flavour variant.
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_scal_Mutual_Inf_gen_all(GRC, List_A, Nsites_A, List_B, Nsites_B, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:,:), INTENT(IN)   :: List_A, List_B
        Integer, Dimension(:,:), INTENT(IN)     :: Nsites_A ,Nsites_B
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

        call Calc_Mutual_Inf_gen_all(GRC, List_A, Nsites_A, List_B, Nsites_B, Renyi_A, Renyi_B, Renyi_AB)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi_A* ZP*ZS
        Obs%Obs_vec(2)   = Obs%Obs_vec(2) + Renyi_B* ZP*ZS
        Obs%Obs_vec(3)   = Obs%Obs_vec(3) + Renyi_AB* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Mutual_Inf_gen_all


      
    End Module Predefined_Obs
