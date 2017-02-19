    Module Natural_Constants

      Real (Kind=Kind(0.d0))  :: eV, amu, Ang, hbar, pi

    contains
     
      subroutine Set_NC
        
        pi = acos(-1.d0)
        eV  = (1.0/6.24150974) *( 10.0**(-18) )
        amu =  1.66053886 * (10.0**(-27))
        Ang =   10.0**(-10)
        hbar = 6.6260755*(10.0**(-34))/(2.0*pi)

      end subroutine Set_NC
    end Module Natural_Constants
