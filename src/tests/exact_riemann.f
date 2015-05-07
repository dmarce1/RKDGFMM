  
      subroutine riemann(rhol, pl, ul, rhor, pr,ur, gamma,t, x, rho,u,p)
      implicit none

c..solves the exact riemann problem
c..original code from bruce fryxell

c..declare
      integer    itmax,iter
      real*8      x, rho, u, p
      real*8       rhol,pl,ul,rhor,pr,ur,gamma,xi,t
      real*8          rho1,p1,u1,rho5,p5,u5,p40,p41,f0,f,eps
      real*8            f1,p4,error,z,c5,gm1,gp1,gmfac1,gmfac2,fact
      real*8           u4,rho4,w,p3,u3,rho3,c1,c3,xsh,xcd,xft,xhd


         if (ul .ne. 0. .or. ur .ne. 0.) then
          write (6,*) 'must have ul = ur = 0.'
          stop
         endif

c..location of discontinuity at t = 0
         xi = 0.0



c..begin solution
      if (pl .gt. pr) then
         rho1 = rhol
         p1   = pl
         u1   = ul
         rho5 = rhor
         p5   = pr
         u5   = ur
      else
         rho1 = rhor
         p1   = pr
         u1   = ur
         rho5 = rhol
         p5   = pl
         u5   = ul
      endif


c..solve for post-shock pressure by secant method
c..initial guesses

         p40 = p1
         p41 = p5
         f0  = f(p40, p1, p5, rho1, rho5, gamma)

c..maximum number of iterations and maxium allowable relative error
         itmax = 20
         eps   = 1.e-5

         do iter = 1, itmax
          f1 = f(p41, p1, p5, rho1, rho5, gamma)
          if (f1 .eq. f0) go to 10

          p4 = p41 - (p41 - p40) * f1 / (f1 - f0)

          error = abs (p4 - p41) / p41
          if (error .lt. eps) go to 10

          p40 = p41
          p41 = p4
          f0  = f1
         enddo
         write (6,*) 'iteration failed to converge'
         stop 'abnormal termination'

10       continue


c..compute post-shock density and velocity
         z  = (p4 / p5 - 1.)
         c5 = sqrt (gamma * p5 / rho5)

         gm1 = gamma - 1.
         gp1 = gamma + 1.
         gmfac1 = 0.5 * gm1 / gamma
         gmfac2 = 0.5 * gp1 / gamma

         fact = sqrt (1. + gmfac2 * z)

         u4 = c5 * z / (gamma * fact)
         rho4 = rho5 * (1. + gmfac2 * z)
     1               / (1. + gmfac1 * z)

c..shock speed
         w = c5 * fact


c..compute values at foot of rarefaction
         p3 = p4
         u3 = u4
         rho3 = rho1 * (p3 / p1)**(1. /gamma)

c..compute positions of waves
      if (pl .gt. pr) then
       c1 = sqrt (gamma * p1 / rho1)
       c3 = sqrt (gamma * p3 / rho3)

       xsh = xi + w * t
       xcd = xi + u3 * t
       xft = xi + (u3 - c3) * t
       xhd = xi - c1 * t




c..compute solution as a function of position
      
        if (x .lt. xhd) then
         rho = rho1
         p   = p1
         u   = u1
        else if (x .lt. xft) then
         u   = 2. / gp1 * (c1 + (x - xi) / t)
         fact   = 1. - 0.5 * gm1 * u / c1
         rho = rho1 * fact ** (2. / gm1)
         p   = p1 * fact ** (2. * gamma / gm1)
        else if (x .lt. xcd) then
         rho = rho3
         p   = p3
         u   = u3
        else if (x .lt. xsh) then
         rho = rho4
         p   = p4
         u   = u4
        else
         rho = rho5
         p   = p5
         u   = u5
        endif
      endif    
 

c..if pr > pl, reverse solution
      if (pr .gt. pl) then
       c1 = sqrt (gamma * p1 / rho1)
       c3 = sqrt (gamma * p3 / rho3)

       xsh = xi - w * t
       xcd = xi - u3 * t
       xft = xi - (u3 - c3) * t
       xhd = xi + c1 * t

      
        if (x .lt. xsh) then
           rho = rho5
           p   = p5
           u   = -u5
        else if (x .lt. xcd) then
           rho = rho4
           p   = p4
           u   = -u4
        else if (x .lt. xft) then
           rho = rho3
           p   = p3
           u   = -u3
        else if (x .lt. xhd) then
           u   = -2. / gp1 * (c1 + (xi - x) / t)
           fact   = 1. + 0.5 * gm1 * u / c1
           rho = rho1 * fact ** (2. / gm1)
           p   = p1 * fact ** (2. * gamma / gm1)
        else
           rho = rho1
           p   = p1
           u   = -u1
        endif
      endif

      end




      real*8 function f(p4, p1, p5, rho1, rho5, gamma)
      implicit none

c..shock tube equation

c..declare the pass
      real*8 p4,p1,p5,rho1,rho5,gamma

c..local variables
      real*8 z,c1,c5,gm1,gp1,g2,fact

      z = (p4 / p5 - 1.)
      c1 = sqrt (gamma * p1 / rho1)
      c5 = sqrt (gamma * p5 / rho5)

      gm1 = gamma - 1.
      gp1 = gamma + 1.
      g2  = 2. * gamma

      fact = gm1 / g2 * (c5 / c1) * z
     1     / sqrt (1. + gp1 / g2 * z)
      fact = (1. - fact) ** (g2 / gm1)

      f = p1 * fact - p4

      return
      end

