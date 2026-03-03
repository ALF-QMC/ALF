  !  Copyright (C) 2018-2026 The ALF project
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
  !     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

    !> @file natural_constants_mod.F90
    !> @author ALF Collaboration
    !> @brief Physical constants and unit conversion factors for ALF simulations.
    !> @details
    !> This module defines fundamental physical constants (π, ℏ) and provides
    !> conversion factors for common units used in quantum materials simulations:
    !> - Energy: electron volt (eV) ↔ Joules
    !> - Mass: atomic mass unit (amu) ↔ kilograms  
    !> - Length: Angstrom (Å) ↔ meters
    !> - Action: reduced Planck constant (ℏ) in SI units (J·s)
    !>
    !> All constants are initialized to SI base units via Set_NC() subroutine.
    !> @see Set_NC()
    Module Natural_Constants

      !> π (pi): Mathematical constant, computed to machine precision
      !> Used in conversions (e.g., ℏ = h/(2π))
      Real (Kind=Kind(0.d0))  :: pi

      !> eV: Electron volt to Joules conversion factor
      !> 1 eV = 1.602176565e-19 J (CODATA 2010), computed as 1/6.24150974 × 10^-18 J
      !> Used to convert energies: E[Joules] = E[eV] × eV
      Real (Kind=Kind(0.d0))  :: eV

      !> amu: Atomic mass unit to kilograms conversion factor
      !> 1 amu = 1.66053886e-27 kg (historical CODATA value)
      !> Used to convert masses: m[kg] = m[amu] × amu
      Real (Kind=Kind(0.d0))  :: amu

      !> Ang: Angstrom to meters conversion factor  
      !> 1 Å = 1e-10 m
      !> Used to convert lengths: L[meters] = L[Å] × Ang
      Real (Kind=Kind(0.d0))  :: Ang

      !> hbar: Reduced Planck constant (ℏ = h/(2π)) in SI units (J·s)
      !> h = 6.62607555e-34 J·s; ℏ = h/(2π) ≈ 1.054571817e-34 J·s
      !> Used in quantum mechanical calculations and action integrals
      Real (Kind=Kind(0.d0))  :: hbar

    contains
     
      !> @brief Initialize natural constants to SI units.
      !> @details Computes fundamental physical constants (π, ℏ) and conversion
      !> factors (eV, amu, Å) to SI base units. These constants are used throughout
      !> ALF for energy scales, mass conversions, and length scales in Hamiltonian
      !> definitions and observable calculations.
      !>
      !> Constants initialized:
      !> - π via acos(-1) to obtain machine-precision value
      !> - eV: electron volt → Joules (inverse of 6.24150974 × 10^18 eV/J)
      !> - amu: atomic mass unit → kilograms
      !> - Å: Angstrom → meters
      !> - ℏ: reduced Planck constant h/(2π) in J·s
      !> @see Natural_Constants module constants
      subroutine Set_NC
        
        ! Compute π to machine precision from acos(-1.0), ensuring accurate
        ! downstream calculations of ℏ and other periodic constants
        pi = acos(-1.d0)
        
        ! Electron volt conversion: 1 eV ≈ 1.602176565e-19 J
        ! Computed as (1.0 / 6.24150974) × 10^-18 J, where 6.24150974 is eV/J
        eV  = (1.0/6.24150974) *( 10.0**(-18) )
        
        ! Atomic mass unit: 1 amu ≈ 1.66053886e-27 kg (CODATA reference)
        ! Used to convert atomic/molecular masses to SI kilograms
        amu =  1.66053886 * (10.0**(-27))
        
        ! Angstrom: 1 Å = 1e-10 m (exact by definition)
        ! Common length scale in condensed matter for lattice spacings
        Ang =   10.0**(-10)
        
        ! Reduced Planck constant: ℏ = h/(2π) with h = 6.6260755e-34 J·s
        ! After computing π above, divide (h/(2π)) to obtain ℏ in J·s
        ! Used in kinetic energy terms and commutation relations
        hbar = 6.6260755*(10.0**(-34))/(2.0*pi)

      end subroutine Set_NC
    end Module Natural_Constants
