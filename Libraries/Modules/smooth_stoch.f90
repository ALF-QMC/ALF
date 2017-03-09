    Program Trans
      
      Implicit Real (Kind=Kind(0.d0)) (A-G,O-Z)
      Implicit Integer (H-N)

      parameter (Ndis=650)
      Real (Kind=Kind(0.d0)) :: Xn_m(Ndis), om(Ndis), Xn_m_new(Ndis)

      open  (Unit=10, File="Aom_ps_20",status="unknown")
      do i = 1,Ndis
         read(10,*) om(i), Xn_m(i), X, Y, Z
      enddo
      close(10)
      
      pi = acos(-1.0)
      Xn_m_new = 0.d0
      Del = om(2) - om(1)
      do nd = 1,Ndis
         weight = Xn_m(nd)
         x_0    = om(nd)
         do i = 1,Ndis
            x = om(i)
            Xn_m_new(i) = Xn_m_new(i) + weight*del*g(10.0*del,x_0,x,pi)
         enddo
      enddo

      do i = 1,Ndis
         write(20,*) Om(i), Xn_m_new(i)
      enddo
    end Program Trans

 
    real (Kind=Kind(0.d0)) function g(del,a,om,pi)
      
      implicit none 
      real (Kind=Kind(0.d0)) :: del, a, om, pi

      g = exp( -((om - a)/del)**2)/(sqrt(pi)*del)

    end function g
