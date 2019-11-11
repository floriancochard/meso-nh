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
      subroutine aeroeq(gas,caero,tempk,rh,ICOLUMN,err)
C***********************************************************************
C   main program for gas/particle equilibrium
C   Using the total bulk aerosol conc. and the gybrid method for dist.
C   calculation including crustal mateials written  by Zhaoyue Meng.   *
C***********************************************************************
C
C...PROGRAM VARIABLES...
C
c             Dp: particle diameter in unit of meter
c             OM: organic mass (ug/m3)
C
C     INPUTS
C     ------
C
c         gas(naspec): Bulk gas-phasec concentrations in units of micro g/m3
c         caero(nasect,naspec): aerosol concentrations in units of micro g/m3
c
c             gas(1) :  Not used in this module
c             gas(2) :  H2SO4
c             gas(3) :  NH3
c             gas(4) :  HNO3
c             gas(5) :  HCl
c             gas(21):  SB1
c             gas(23):  SB2
c             gas(25):  SB3
c             gas(26):  SB4
c             gas(27):  SB5
c             gas(28):  SA1
c             gas(31):  SA2
c             gas(34):  SA3
c             gas(35):  SA4
c             gas(37):  SA5
c             caero(1) : Na
c             caero(2) : SO4
c             caero(3) : NH4
c             caero(4) : NO3
c             caero(5) : Cl 
c             caero(6) : K   
c             caero(7) : Ca  
c             caero(8) : Mg  
c             caero(9) : H2CO3  
c             caero(10): H2O 
c             caero(11): other inorganic species
c             caero(12): EC  
c             caero(13): POA1
c             caero(14): POA2
c             caero(15): POA3
c             caero(16): POA4
c             caero(17): POA5
c             caero(18): POA6
c             caero(19): POA7
c             caero(20): POA8
c             caero(21): SB1
c             caero(22): SB1-
c             caero(23): SB2
c             caero(24): SB2-
c             caero(25): SB3
c             caero(26): SB4
c             caero(27): SB5
c             caero(28): SA1
c             caero(29): SA1-
c             caero(30): SA1--
c             caero(31): SA2
c             caero(32): SA2-
c             caero(33): SA2--
c             caero(34): SA3
c             caero(35): SA4
c             caero(36): SA5
c             caero(37): SA5-
c             tempk    : temperature in Kelvin
c             rh       : relative humidity
c             ICOLUMN  : Column index
C
C     OUTPUTS
C     -------
C             The updated caero and gas based on equilibrium calculations.
c
      include 'parameter.h'
      dimension gas(naspec),caero(nasect,naspec),caero0(nasect,naspec)
      real nh3,nitrate,magn,w(9),newvol,negconc,temp,hconc
      real mass(nasect),NN(nasect),Kn,csoa(nasect,17)
      real dmass(nasect),factor(nasect),lwc,zsrconc(10),dellwc
      real tollwc,denom(8),frac(8),dpindex(nasect)
      real cpart(10),cpt(8),sumsoa,ctsoa(17)
      real org(10),aq(10),am(10),ad(10),vp(10),mwpart(10)
      logical done
      common/cond/temp,lwc,hconc
      common/concsoa/cpart,cpt
      common/props/vp,mwpart
      common/fraccalc/NN,DPINDEX
      data Diffus,alfa/1.e-5, 0.01/
      data pi/3.1415926/
c Calculate the bulk aerosol composition
      do j=1, 8
        w(j) = gas(j)
        do inasect = 1, nasect
           w(j) = w(j) + caero(inasect,j)
        enddo
      enddo
      do j=9, 10 !Use gas(9) and (10) hold the old total carbonate and wat.
        gas(j) = 0.
        do inasect = 1, nasect
           gas(j) = gas(j) + caero(inasect,j)
        enddo
      enddo
c
        w(9) = 9.12e5
        call eq(iapart,iaopt,iaapp,iasec,iatx,aeps1,aeps2,tempk,rh,
     1  w,ammon,sodium,sulfate,nitrate,chloride,bisulf,pott,calcium,
     2  magn,nh3,hno3,hcl,rh2so4,hydro,roh,rna2so4,rnh42s4,rnh4cl,
     3  rnacl,rnano3,rnh4no3,rclc,rnh4hso4,rnahso4,rkcl,rk2so4,rkhso4,
     4  rkno3,rcacl2,rcaso4,rcano32,rmgcl2,rmgso4,rmgno32,rwater,
     5  rna2co3,rnahco3,rk2co3,rkhco3,rcaco3,rmgco3,rcarbamate,rnh4hco3,
     6  co2,rco2,carbonate,bicarbon,carbamation,rnh3,rhno3,rhcl,ierr,xh)
c
        if(ierr .eq. 1) then
          write(6,*)w,dtdummy,deltat
          write(6,*)caero
          write(6,*)nn(inasect),mass(inasect),inasect
          return
        endif
c Reserve the old aerosol concentrations
      do i = 1, nasect
        do j = 1, naspec
           caero0(i,j) = caero(i,j)
        enddo
      enddo
c Calculate mass transport factor
      delmu=(alog(Dpup)-alog(Dplow))/float(nasect)
      fact = 0.
      do inasect = 1, nasect
        Dp=alog(Dplow)+(float(inasect)-0.5)*delmu
        Dp=exp(Dp)*1.e-6
        DPINDEX(INASECT) = dp
        totmass=0.
        do j=1,naspec
          totmass=totmass+caero(inasect,j)
        enddo
        mass(inasect)=pi/6.*Dp**3*densp*1.e12
        NN(inasect)=totmass/mass(inasect)
        NN(inasect) = max(NN(inasect),1.)
        call aeroparam(tempK-273.15,inasect,1.,Kn,Sc,St,Vsett)
        beta=Kn/alfa
c Using the hybrid method to distribute sulfate
        factor(inasect) = 2.*pi*Dp*Diffus*NN(inasect)/(1.+beta)
        fact = fact + factor(inasect)
      enddo
      do i = 1, nasect
        factor(i) = factor(i)/fact
        caero(i,2)=caero(i,2) + gas(2)*factor(i)
        caero(i,3)=caero(i,3) + (gas(3)-nh3)*factor(i)
        caero(i,4)=caero(i,4) + (gas(4)-hno3)*factor(i)
        caero(i,5)=caero(i,5) + (gas(5)-hcl)*factor(i)
c w(9) has been changed in 'eq' to for aerosol phase only.
        caero(i,9)=caero(i,9) + (w(9)-gas(9))*factor(i)
        caero(i,10)=caero(i,10) + (rwater-gas(10))*factor(i)
      enddo
      gas(2) = 0.
      gas(3) = nh3
      gas(4) = hno3
      gas(5) = hcl
c Now if any particles have negative concentrations redistribute the
c materials over other sized particles and set the minimum conc to 0.
      do j = 1, naspec
         negconc = 0.
         totconc = 0.
         do i = 1, nasect
            if(caero(i,j) .lt. 0.) then
              negconc = negconc + caero(i,j)
              caero(i,j) = 0.
            endif
            totconc = totconc + caero(i,j)
            totconc = max(1.e-6,totconc)
         enddo
         do i = 1, nasect
            caero(i,j) = caero(i,j)*(1.+negconc/totconc)
         enddo
      enddo
c    
c     if(icolumn .eq. 34) stop
c     if(icolumn .eq. 33) then
c     write(6,*) 'gas(2),gas(3),(4), (5) = ',(gas(i),i=2,5)
c     write(6,*) 'nh3,hno3,hcl = ',nh3,hno3,hcl
c     write(6,*) 'Inorganic aerosol conc = '
c     write(6,'(11e8.3)') (w(i),i=1,9),rh,rwater
c     endif
c
c
 45   continue
c  
      do i = 1, nasect
         dmass(i) = 0.
         do j = 1, naspec
            dmass(i) = dmass(i) + caero(i,j) - caero0(i,j)
         enddo
      enddo
c move the aerosol concentrations to caero0 and then redistribute the
c moving sections to the fixed caero(i,j)
      do i = 1, nasect
         do j = 1, naspec
            caero0(i,j) = caero(i,j)
            caero(i,j) = 0.
         enddo
      enddo
c        
      do inasect = 1, nasect
        Dp=alog(Dplow)+(float(inasect)-0.5)*delmu
        Dp=exp(Dp)*1.e-6
         newvol = Dp**3 + 6./pi*dmass(inasect)/NN(inasect)*1.e-12/densp
         if(newvol .lt. 0.) newvol = 0.
         Dp1 = newvol**(1./3.)
         Dp1 = Dp1*1.e6
c        write(6,*)icolumn,inasect,Dp*1.e6,Dp1
         Dp1 = max(Dp1, Dplow)
         Dpmove = (alog(Dp1) - alog(Dplow))/delmu + 0.5
         if(Dpmove .lt. 1.) Dpmove = 1.000001
         if(Dpmove .gt. float(nasect)) Dpmove = float(nasect)+0.000001
         imove = int(Dpmove)
         distr = Dpmove - float(imove)
         do j = 1, naspec
           caero(imove,j)=caero(imove,j) + (1.-distr)*caero0(inasect,j)
           if(imove .ne. nasect) caero(imove+1,j) = caero(imove+1,j) 
     &                                    + distr*caero0(inasect,j)
         enddo
      enddo
C
      do k = 13,20
         cpt(k-12) = 0.0
         do j = 1,nasect
            cpt(k-12) = cpt(k-12) + Caero(j,k)
         enddo   
      enddo   
C
      CALL MAINPART(ORG,AQ,AM,AD,prevsol)
C     
      DO K = 1,10
         ZSRCONC(K) = AQ(K) + AM(K)*MWPART(K)/(MWPART(K)-1.) +
     +       AD(K)*MWPART(K)/(MWPART(K)-2.)
      ENDDO   
C
c
       CALL ZSR(ZSRCONC,RH,DELLWC)       
C
       CTSOA(1) = ORG(1) + AQ(1)
       CTSOA(2) = AM(1)*MWPART(1)/(MWPART(1)-1.)
       CTSOA(3) = ORG(2) + AQ(2)
       CTSOA(4) = AM(2)*MWPART(2)/(MWPART(2)-1.)
       CTSOA(5) = ORG(3) + AQ(3)
       CTSOA(6) = ORG(4) + AQ(4)
       CTSOA(7) = ORG(5) + AQ(5)
       CTSOA(8) = ORG(6) + AQ(6)
       CTSOA(9) = AM(6)*MWPART(6)/(MWPART(6)-1.)
       CTSOA(10) = AD(6)*MWPART(6)/(MWPART(6)-2.)
       CTSOA(11) = ORG(7) + AQ(7)
       CTSOA(12) = AM(7)*MWPART(7)/(MWPART(7)-1.)
       CTSOA(13) = AD(7)*MWPART(7)/(MWPART(7)-2.)
       CTSOA(14) = ORG(8) + AQ(8)
       CTSOA(15) = ORG(9) + AQ(9)
       CTSOA(16) = AM(9)*MWPART(9)/(MWPART(9)-1.)
       CTSOA(17) = ORG(10) + AQ(10)
C
      DO J = 1,NASECT
         DENOM(J) = 0.0
      enddo   
C
      DO I = 1,NASECT
         DENOM(I) = PI*DPINDEX(I)*DPINDEX(I)*NN(I)
      ENDDO   
C
      BOT = 0.0
      DO I = 1,NASECT
         BOT = BOT + DENOM(I)
      ENDDO
C
      DO I = 1,NASECT
         FRAC(I) = DENOM(I)/BOT
      ENDDO   
C
      DO I = 1,17
      DO J = 1,NASECT
         CSOA(J,I) = CTSOA(I)*FRAC(J)
      ENDDO   
      ENDDO
C
C
       DO J = 1,NASECT
          CAERO(J,10) = CAERO(J,10) + DELLWC*FRAC(J)
       ENDDO   
C
      return
      end
