    Module Natural_Constants

      ! Mathematical constants
      Real (Kind=Kind(0.d0)), parameter :: pi    = acos(-1.d0)
      Real (Kind=Kind(0.d0)), parameter :: twopi = 2.d0 * pi

      ! Tolerance constants
      Real (Kind=Kind(0.d0)), parameter :: Eps_machine     = 1.D-12  ! Near machine precision: eigenvalue zero-detection, floating-point identity
      Real (Kind=Kind(0.d0)), parameter :: Eps_small        = 1.D-10  ! Negligibility: parameter on/off gating, singularity protection
      Real (Kind=Kind(0.d0)), parameter :: Eps_convergence = 1.D-6   ! Iterative convergence, coarse k-space detection

      ! Physical constants
      Real (Kind=Kind(0.d0)), parameter :: eV = (1.0/6.24150974) *( 10.0**(-18) )
      Real (Kind=Kind(0.d0)), parameter :: amu =  1.66053886 * (10.0**(-27))
      Real (Kind=Kind(0.d0)), parameter :: Ang =   10.0**(-10)
      Real (Kind=Kind(0.d0)), parameter :: hbar = 6.6260755*(10.0**(-34))/(2.0*pi)

    end Module Natural_Constants
