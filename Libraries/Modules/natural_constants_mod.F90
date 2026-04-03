    Module Natural_Constants

      Real (Kind=Kind(0.d0)), parameter :: pi = acos(-1.d0)
      Real (Kind=Kind(0.d0)), parameter :: eV = (1.0/6.24150974) *( 10.0**(-18) )
      Real (Kind=Kind(0.d0)), parameter :: amu =  1.66053886 * (10.0**(-27))
      Real (Kind=Kind(0.d0)), parameter :: Ang =   10.0**(-10)
      Real (Kind=Kind(0.d0)), parameter :: hbar = 6.6260755*(10.0**(-34))/(2.0*pi)

    contains
     
      subroutine Set_NC
         ! This subroutine is now redundant since all constants are parameters
      end subroutine Set_NC

    end Module Natural_Constants
