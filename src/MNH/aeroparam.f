!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 aerosol 2006/05/18 13:07:25
!-----------------------------------------------------------------
      subroutine aeroparam(tempE, iaero, Ustar, Kn, Sc, St, Vsett)
C     This subroutine calculates the Knudson, Schmidt, and Stokes numbers
C     as a function of the inpput temperature and particle size.
C     Written by Zhaoyue Meng in May 1996 at Caltech.
C
C     Inputs:
C        tempE  : Temperature (Centigrade)
C        iaero  : Aerosol size section index
C        Ustar  : Friction velocity
C
C     Outputs:
C        Kn     : Knudson number
C        Sc     : Particle Schmidt number
C        St     : Stokes number
C        Vsett  : Particle settling velocity
C
C     Variables used:
C        p       : Atmospheric pressure
C        Rgas    : Gas constant
C        viscos  : Viscosity of air
C        delmu   : Interval of particle size sections on a logrithmic scale
C        Dp      : Particle diameter (m)
C        mfp     : Mean free path of air molecules
C        DiffuP  : Diffusivity of particles
C        densP   : Density of particles
C        densA   : Density of air
C        Vsett   : Particle settling velocity
C        Ustar   : Friction velocity

      INCLUDE 'parameter.h'
      real Kn, Mwair, KBolzm, mfp
      common / mfp / mfp
      external fundiff
      data p,Rgas,Mwair/1.013e5,8.314,29.e-3/
      data pi/3.141593/
      data viscos/1.8e-5/ !in unit of kg/m/s
      data KBolzm,grav/1.38e-23, 9.80/
      data alfa,beta,gama/1.257,0.40,-1.10/
         tempK = tempE + 273.15
         delmu=(alog(Dpup)-alog(Dplow))/float(nasect)
         Dp1=alog(Dplow)+(float(iaero)-1.0)*delmu
         Dp2=Dp1+delmu
         Dp1=exp(Dp1)*1.e-6
         Dp2=exp(Dp2)*1.e-6
         mfp = viscos/0.499/p/sqrt(8.*Mwair/pi/Rgas/tempK)
         aa = DensP*1000.*grav/18./viscos
         Vsett = aa/2./delmu*((Dp2*Dp2-Dp1*Dp1) + 4.*mfp*alfa*(Dp2-Dp1)
     &                  +8.*beta*mfp*mfp/gama*(exp(gama*Dp2/2./mfp)
     &                                        -exp(gama*Dp1/2./mfp)))
         Dp = sqrt(Dp1*Dp2)
         Kn = 2.*mfp/Dp
c        CcKn = Cc(Kn)
c        fff = 3.*pi*Dp*Viscos/CcKn
c        DiffuP1 = KBolzm*TempK/fff
c        Vsett1 = Dp*Dp*DensP*1000.*grav*CcKn/18./viscos
         densA = p/Rgas/tempK*Mwair
         St = Vsett/grav*Ustar*Ustar/Viscos*densA
         bb = KBolzm*TempK/(3.*pi*Viscos)
         call qgaus(fundiff,Dp1,Dp2,ss)
         DiffuP = bb/delmu*((1./Dp1-1./Dp2) + 
     &        alfa*mfp*(1./Dp1/Dp1-1./Dp2/Dp2) + 2.*mfp*beta*ss)
         SC = Viscos/densA/DiffuP
      return
      end
C**********************
c     function Cc(Kn)
c     real Kn
c     data alfa,beta,gama/1.257,0.40,1.10/
c     Cc = 1. + Kn*(alfa + beta*exp(-gama/Kn))
c     return
c     end
C**********************
      function fundiff(Dp)
      real mfp
      common / mfp / mfp
      data alfa,beta,gama/1.257,0.40,-1.10/
      fundiff = exp(gama*Dp/2./mfp)/Dp/Dp/Dp
      return
      end
c Gaussian quadrature integration solver
      subroutine qgaus(func,a,b,ss)
      real a, b, ss, func
      external func
      real dx,xm,xr,w(5),x(5)
      data w/.2955242247,.2692667193,.2190863625,.1494513491,.066671344/
      data x/.1488743389,.4333953941,.6794095682,.8650633666,.973906529/
      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0.
      do 11 j = 1, 5
            dx = xr*x(j)
            ss = ss + w(j)*(func(xm+dx)+func(xm-dx))
 11   continue
      ss = xr*ss
      return
      end
