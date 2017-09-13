ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This file contains the functions and subroutines for the calculation
c   of heavy quark energy loss due to single gluon radiation
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function alphas(kT0,temp0)

c This is a function to calculate the strong coupling costant, temp is the 
c medium temperature, and alphas=alphas(pi*temp) if kT<pi*temp due to 
c screening effect. It's a first order calculation.

      implicit none

      include 'rad_coms.in'
      double precision alphas,kT0,kT1,temp0,nflavor,lambdas,error_para

      error_para=1.d0

      if(kT0.lt.PI*temp0*error_para) then
         kT1=PI*temp0*error_para
      else
         kT1=kT0
      endif

      alphas=4.d0*PI/(11.d0-2.d0*nflavor(kT1)/3.d0)/2.d0/
     &       Log(kT1/lambdas(kT1))

c      write(*,*) alphas

c      alphas=fixAlphas

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function nflavor(kT0)

      implicit none

      include 'rad_coms.in'
      double precision kT0,nflavor

      if (kT0.lt.cMass) then
         nflavor=3.d0
      elseif (kT0.lt.bMass) then
         nflavor=4.d0
      else
         nflavor=5.d0
      endif

c      write(*,*) nflavor

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function lambdas(kT0)

      implicit none

      include 'rad_coms.in'
      double precision kT0,lambdas

      if (kT0.lt.cMass) then
         lambdas=0.2d0
      elseif (kT0.lt.bMass) then
         lambdas=0.172508d0
      else
         lambdas=0.130719d0
      endif

c      write(*,*) lambdas

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function splittingP(z0)

      implicit none

      include 'rad_coms.in'
      
      double precision splittingP,z0

      splittingP = (2.d0-2.d0*z0+z0*z0)*C_F/z0
c      splittingP = (1d0-z0)*(2.d0-2.d0*z0+z0*z0)*C_F/z0

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function tau_f(x0,y0)

      implicit none
      
      include 'rad_coms.in'

      double precision tau_f,x0,y0
      tau_f = 2.d0*HQenergy*x0*(1.d0-x0)/((x0*y0*HQenergy)**2+x0*x0*
     &        HQmass*HQmass)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function dNg_over_dxdydt(argument)

c single gluon radiation formula from Wang's paper

      implicit none

      include 'rad_coms.in'

      double precision dNg_over_dxdydt,alphas,splittingP,tau_f
      double precision x0,y0,argument(2),qhat

      x0=argument(1)
      y0=argument(2)

      qhat=4.d0*alpha*temp**3/C_F

      if(x0*HQenergy.lt.PI*temp) then
c      if(x0*HQenergy.lt.sqrt(6d0*PI*fixAlphas)*temp) then  ! mu_D
         dNg_over_dxdydt = 1.0E-12
         goto 111
      endif

      if(tau_f(x0,y0).lt.1.d0/PI/temp) then
         dNg_over_dxdydt = 1.0E-12
         goto 111
      endif

c      if(x0*y0*HQenergy.lt.sqrt(2d0/3d0)*PI*temp) then
c         dNg_over_dxdydt = 1.0E-12
c         goto 111
c      endif

      dNg_over_dxdydt = 4.d0/PI*N_c*alphas(x0*y0*HQenergy,temp)*
     &                 splittingP(x0)*qhat*y0**5*(sin((time-time_init)
     &                 /2.d0/tau_f(x0,y0)/hbarc))**2*
     &                 (HQenergy*HQenergy/(y0*y0*HQenergy*HQenergy+
     &                 HQmass*HQmass))**4/x0/x0/HQenergy/HQenergy/hbarc

c      write(6,*) dNg_over_dxdydt

 111  continue

      if(dNg_over_dxdydt.gt.max_dNg_over_dxdydt) 
     &     max_dNg_over_dxdydt=dNg_over_dxdydt

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
