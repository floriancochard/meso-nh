!MNH_LIC Copyright 2000-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     algorithme initial créé par Michael Mishchenko (2000)
!
!     algorithme modifié par Corinne Burlaud (2000) puis Olivier Brunau (2002)
!     calcul des paramètres radar en fonction de la forme des gouttes, du taux
!     de pluies, de l'élévation du radar, de l'oscillation des gouttes, 
!     de la longueur d'onde émise, pour une distribution de taille
!     de goutte réelle ou théorique.
!     l'essentiel des paramètres modifiables sont demandés à l'utilisateur
!     pendant l'execution du programme.
!
!     Modif par Olivier Caumont (04/2008) pour interfaçage avec diagnostic 
!     radar de Méso-NH.
!     P. Wautelet 22/01/2019: replace double precision declarations by real(kind(0.0d0)) (to allow compilation by NAG compiler)
!
!****************************************************************************

!     This test result was calculated by the code in its
!     curent setting.
!
!     ICHOICE=1  NCHECK=1
!     RAT= .963711
!     PROLATE SPHEROIDS, A/B=   .5000000
!     LAM=  6.283185   MRR= .1500D+01   MRI= .2000D-01
!     ACCURACY OF COMPUTATIONS DDELT =  .10D-02
!     EQUAL-SURFACE-AREA-SPHERE RADIUS= 10.0000
!     thet0= 56.00  thet= 65.00  phi0=114.00  phi=128.00  alpha=145.00
!     beta= 52.00
!     AMPLITUDE MATRIX
!     S11=-.50941D+01 + i* .24402D+02
!     S12=-.19425D+01 + i* .19971D+01
!     S21=-.11521D+01 + i*-.30977D+01
!     S22=-.69323D+01 + i* .24748D+02
!     PHASE MATRIX
!     650  .3172  -17.9846   10.0498  -12.7580
!     -21.1462  631.6322 -127.3059   87.2144
!     6    .8322  132.6131  635.2767  -34.7730
!     -9.6629  -78.1229   51.4094  643.1738
!     time =     .03 min
!
!   CALCULATION OF THE AMPLITUDE AND PHASE MATRICES FOR                 
!   A PARTICLE WITH AN AXIALLY SYMMETRIC SHAPE                   
!                                                                      
!   This version of the code uses REAL variables,          
!   is applicable to spheroids, finite circular cylinders,            
!   Chebyshev particles, and generalized Chebyshev particles
!   (distorted water drops), and must be used along with the 
!   accompanying file ampld.par.f.                
!                                                                      
!   The code has been developed by Michael Mishchenko at the NASA      
!   Goddard Institute for Space Studies, New York.  The development    
!   of the code was supported by the NASA Radiation Science Program
!   and the NASA FIRE III program.                                                  
!   The code may be used without limitations in any not-for-      
!   profit scientific research.  The only request is that in any       
!   publication using the code the source of the code be acknowledged  
!   and relevant references be made.                                   
!                                                                       
!   The computational method is based on the T-matrix approach         
!   [P. C. Waterman, Phys. Rev. D 3, 825 (1971)], also known as        
!   the extended boundary condition method.  The method was            
!   improved in the following papers:                         
!                                                                       
!  1.  M. I. Mishchenko and L. D. Travis, T-matrix computations       
!      of light scattering by large spheroidal particles,             
!      Opt. Commun., vol. 109, 16-21 (1994).                          
!                                                                      
!   2.  M. I. Mishchenko, L. D. Travis, and A. Macke, Scattering       
!       of light by polydisperse, randomly oriented, finite            
!       circular cylinders, Appl. Opt., vol. 35, 4927-4940 (1996).     
!                                                                      
!   3.  D. J. Wielaard, M. I. Mishchenko, A. Macke, and B. E. Carlson, 
!       Improved T-matrix computations for large, nonabsorbing and     
!       weakly absorbing nonspherical particles and comparison         
!       with geometrical optics approximation, Appl. Opt., vol. 36,    
!       4305-4313 (1997).                                             
!                                                                      
!   A general review of the T-matrix approach can be found in          
!                                                                      
!   4.  M. I. Mishchenko, L. D. Travis, and D. W. Mackowski,           
!       T-matrix computations of light scattering by nonspherical      
!       particles:  a review, J. Quant. Spectrosc. Radiat.             
!       Transfer, vol. 55, 535-575 (1996).                             
!                                                                      
!   Additional useful information is contained in the paper
!                                                                      
!   5.  M. I. Mishchenko and L. D. Travis, Capabilities and            
!       limitations of a current FORTRAN implementation of the         
!       T-matrix method for randomly oriented, rotationally            
!       symmetric scatterers, J. Quant. Spectrosc. Radiat. Transfer,   
!       vol. 60, 309-324 (1998).                                       
!                                                                      
!   The definitions and notation used can also be found in
!
!   6.  M. I. Mishchenko, J. W. Hovenier, and L. D. Travis, Concepts,
!       terms, notation.  In "Light Scattering by Nonspherical
!       Particles:  Theory, Measurements, and Applications," 
!       edited by M. I. Mishchenko, J. W. Hovenier, and L. D. Travis,
!      Acedemic Press, San Diego, 1999, pp. 3-27
!
!   and especially
!
!   7.  M. I. Mishchenko, Calculation of the amplitude matrix
!       for a nonspherical particle in a fixed orientation,
!       Appl. Opt. vol. 39, 1026-1031 (2000).
!
!   Copies of these papers are available upon request from Michael     
!   Mishchenko.  Please send your request to crmim@giss.nasa.gov.      
!   They are also available in the PDF format at 
!   http://www/giss/nasa.gov/~crmim (button "Publications On-Line")
!                                                                      
!   One of the main advantages of the T-matrix method is that the      
!   T-matrix for a given nonspherical particle needs to be computed    
!   only once and then can be used any number of times for computing   
!   the amplitude and phase matrices for any directions of the incident 
!   and scattered beams.  This makes the T-matrix method           
!   extremely convenient and efficient in averaging over particle      
!   orientations and/or directions of incidence (or scattering).       
!                                                                       
!   The use of extended precision variables (Ref. 1) can significantly
!   increase the maximum convergent equivalent-sphere size parameter 
!   and make it well larger than 100 (depending on refractive index    
!   and aspect ratio).  The extended-precision code is also available. 
!   However, the use of extended precision variables results in a      
!   greater consumption of CPU time.                                   
!   On IBM RISC workstations, that code is approximately               
!   five times slower than this double-precision code.  The            
!   CPU time difference between the double-precision and extended-     
!   precision codes can be larger on supercomputers.                   
!
!                                                              
!   As described in Ref. 3 above, these changes will affect the        
!   performance of the code only for nonabsorbing or weakly            
!   absorbing particles (imaginary part of the refractive              
!   index smaller than about 0.001).                                   
! 
!   
!      
!   INPUT PARAMETERS:                                                  
!                                                                     
!      AXI - equivalent-sphere radius                                  
!      RAT = 1 - particle size is specified in terms of the            
!                equal-volume-sphere radius                             
!      RAT.ne.1 - particle size is specified in terms of the           
!                equal-surface-area-sphere radius                      
!      LAM - WAVELENGTH OF INCIDENT LIGHT                                      
!      MRR and MRI - real and imaginary parts of the refractive        
!                  index                                               
!      EPS and NP - specify the shape of the particles.                
!             For spheroids NP=-1 and EPS is the ratio of the          
!                 horizontal to rotational axes.  EPS is larger than   
!                 1 for oblate spheroids and smaller than 1 for        
!                 prolate spheroids.                                   
!             For cylinders NP=-2 and EPS is the ratio of the          
!                 diameter to the length.                              
!             For Chebyshev particles NP must be even and positive and 
!                 is the degree of the Chebyshev polynomial, while     
!                 EPS is the deformation parameter (Ref. 5).                   
!             For generalized Chebyshev particles (describing the shape
!                 of distorted water drops) NP=-3.  The coefficients
!                 of the Chebyshev polynomial expansion of the particle
!                 shape (Ref. 7) are specified in subroutine DROP.
!      DDELT - accuracy of the computations                            
!      NDGS - parameter controlling the number of division points      
!             in computing integrals over the particle surface (Ref. 5).       
!             For compact particles, the recommended value is 2.       
!             For highly aspherical particles larger values (3, 4,...) 
!             may be necessary to obtain convergence.                  
!             The code does not check convergence over this parameter. 
!             Therefore, control comparisons of results obtained with  
!             different NDGS-values are recommended.                   
!      ALPHA and BETA - Euler angles (in degrees) specifying the orientation 
!      of the scattering particle relative to the laboratory reference
!         frame (Refs. 6 and 7).
!      THET0 - zenith angle of the incident beam in degrees
!      THET - zenith angle of the scattered beam in degrees    
!      PHI0 - azimuth angle of the incident beam in degrees    
!      PHI - azimuth angle of the scattered beam in degrees   
!            (Refs. 6 and 7)
!      ALPHA, BETA, THET0, THET, PHI0, and PHI are specified at the end of
!      the main program before the line                                    
!                                                                      
!       "CALL AMPL (NMAX,...)"                     
!                                                                      
!      The part of the main program following the line 
!                                                                      
!       "COMPUTATION OF THE AMPLITUDE AND PHASE MATRICES"               
!                                                                      
!      can be repeated any number of times.  At this point the T-matrix 
!      for the given scattering particle has already     
!      been fully computed and can be repeatedly used in computations  
!      for any directions of illumination and scattering and any particle
!      orientations.              
!                                                                       
!   OUTPUT PARAMETERS:                                                 
!                                                                      
!      Elements of the 2x2 amplitude matrix       
!      Elements of the 4x4 phase matrix
!                                                                      
!   Note that LAM and AXI must be given in the same units of length        
!   (e.g., microns).                                                          
!                                                                      
!   The convergence of the T-matrix method for particles with          
!   different sizes, refractive indices, and aspect ratios can be      
!   dramatically different.  Usually, large sizes and large aspect     
!   ratios cause problems.  The user of this code                      
!   should "play" a little with different input parameters in          
!   order to get an idea of the range of applicability of this         
!   technique.  Sometimes decreasing the aspect ratio                  
!   from 3 to 2 can increase the maximum convergent equivalent-        
!   sphere size parameter by a factor of several.                      
!                                                                      
!   The CPU time required rapidly increases with increasing ratio      
!   radius/wavelength and/or with increasing particle asphericity.     
!   This should be taken into account in planning massive computations.
!                                                                      
!   Execution can be automatically terminated if dimensions of certain 
!   arrays are not big enough or if the convergence procedure decides  
!   that the accuracy of double-precision variables is insufficient  
!   to obtain a converged T-matrix solution for given particles.       
!   In all cases, a message appears explaining                         
!   the cause of termination.                                          
!                                                                      
!   The message                                                        
!        "WARNING:  W IS GREATER THAN 1"                               
!   means that the single-scattering albedo exceeds the maximum        
!   possible value 1.  If W is greater than 1 by more than             
!   DDELT, this message can be an indication of numerical              
!   instability caused by extreme values of particle parameters.       
!                                                                      
!   The message "WARNING: NGAUSS=NPNG1" means that convergence over    
!   the parameter NG (see Ref. 2) cannot be obtained for the NPNG1     
!   value specified in the PARAMETER statement in the file ampld.par.f.
!   Often this is not a serious problem, especially for compact
!   particles.
!                                                                      
!   Larger and/or more aspherical particles may require larger
!   values of the parameters NPN1, NPN4, and NPNG1 in the file
!   ampld.par.f.  It is recommended to keep NPN1=NPN4+25 and
!   NPNG1=3*NPN1.  Note that the memory requirement increases
!   as the third power of NPN4. If the memory of
!   a computer is too small to accomodate the code in its current
!   setting, the parameters NPN1, NPN4, and NPNG1 should be
!   decreased. However, this will decrease the maximum size parameter
!   that can be handled by the code.
!
!   In some cases any increases of NPN1 will not make the T-matrix     
!   computations convergent.  This means that the particle is just     
!   too "bad" (extreme size parameter and/or extreme aspect ratio      
!   and/or extreme refractive index).                                  
!   The main program contains several PRINT statements which are       
!   currently commentd out.  If uncommented, these statements will     
!   produce numbers which show the convergence rate and can be         
!   used to determine whether T-matrix computations for given particle 
!   parameters will converge at all, whatever the parameter NPN1 is.   
!                                                                       
!   Some of the common blocks are used to save memory rather than      
!   to transfer data.  Therefore, if a compiler produces a warning     
!   message that the lengths of a common block are different in        
!   different subroutines, this is not a real problem.                 
!                                                                       
!   The recommended value of DDELT is 0.001.  For bigger values,       
!   false convergence can be obtained.                                 
!                                                                       
!   In computations for spheres, use EPS=1.000001 instead of EPS=1.    
!   EPS=1 can cause overflows in some rare cases.                      
!                                                                       
!   For some compilers, DACOS must be raplaced by DARCOS and DASIN     
!   by DARSIN.                                                         
!                                                                       
!   I would highly appreciate informing me of any problems encountered 
!   with this code.  Please send your message to the following         
!   e-mail address:  CRMIM@GISS.NASA.GOV.                              
!
!   WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
!   IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
!   INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
!   VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
!   THE AUTHORS AND THEIR ORGANIZATION DISCLAIM ALL LIABILITY FOR
!   ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM. 


      SUBROUTINE TMD(Deq,LAM,MRR,MRI,NDGS,OSCIL,ELEV,qbackHH,&
          qbackVV,kdp,QQEXT,EPS)


      USE MODD_TMAT,ONLY : XRT11,XRT12,XRT21,XRT22,&
                           XIT11,XIT12,XIT21,XIT22,&
                           XTR1,XTI1,NPN1,NPNG1,NPNG2,NPN2,NPN4,NPN6

      IMPLICIT REAL*8 (A-H,O-Z)

!!      Parameter (NPN1=200, NPNG1=600, NPNG2=2*NPNG1, NPN2=2*NPN1,&
!!          NPN4=NPN1, NPN6=NPN4+1)
      


      INTEGER param,oscil
      REAL*8 SIGBETA,Fbeta
      REAL*8 PDtotal,PDalpha,PDbeta,Poids

      REAL*8  LAM,MRR,MRI,Deq,X(NPNG2),W(NPNG2),S(NPNG2),SS(NPNG2),&
             AN(NPN1),R(NPNG2),DR(NPNG2),&
             DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)

      REAL*8 NUMZ,DENZ,NUMTZ,DENTZ,ZDRT,NTotal
      REAL*8 NUMD,DEND,DELTA,NUMTD,DENTD,DELTAT
      REAL*8 PAS,KDP

      REAL*8 THET0,THET,Elev,AXI

      REAL*8 MYS11cr,MYS22cr


!!      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
!!     REAL*4 RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),&
!!          RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),&
! !         IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),&
!!          IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
      COMPLEX*16 S11,S12,S21,S22
      COMPLEX*16 S11u,S12u,S21u,S22u

      REAL*8 S11carre,S22carre
      COMPLEX*16 NUMrhoAB,NUMTrhoAB
      
!!      COMMON /CT/ TR1,TI1
!!      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
      LOGICAL SORTIE1
                
!****************************************************************

! initialisation des variables

!c      IF (freq.EQ.1) then
!
!c     I-Variables depending on the wavelength:
!C     Complex refractive index
!c         write(8,9887)
!c 9887    format(20x,'bande S')         
!C     1)  BANDE S
!c         MRR=9.006 D0
!c         MRI=0.9531 D0
!c         LAM=10.7D-2
!c         K2=0.9313
!C     imag(-K)
!c         IMK=0.00688         
!c      elseif(freq.EQ.2) then         
!c         write(8,9888)
!c 9888    format(20x,'bande X')          
!C     2) BANDE X
!c         MRR=7.813D0
!c         MRI=2.427D0
!c         LAM=3.20D-2
!c         K2=0.9282 
!c         IMK=0.0247         
!c      endif
!      
!c      K2=ABS(((MRR+(0,1)*MRI)**2-1.)/((MRR+(0,1)*MRI)**2+2.))**2
!C      IMK=AIMAG(((MRR+(0,1)*MRI)**2-1.)/((MRR+(0,1)*MRI)**2+2.))
!c      PRINT*,K2,IMK,'>0'

      P=ACOS(-1D0)             !calcul de pi!
      
!****  Lecture du fichier de DSD**************************************       
!*********************spectre théorique**********************
!*     constante en m-3.mm-1
!c      N0=8000D0
!      
!*     taux de precipitation en mm/h
      
!!c      R_lu(1)=0.1D0
!c      R_lu(2)=0.5D0
!c      R_lu(3)=1D0
!c      R_lu(4)=5D0
!c      DO I=1,5   
!c         R_lu(I+4)=R_lu(I+3)+5D0
!c      ENDDO
!      
!c      DO I=1,5
!c         R_lu(I+9)=R_lu(I+8)+10D0
!c      ENDDO
!      
!c      R_lu(15)=100D0
!c      R_lu(16)=125D0
!c      R_lu(17)=150D0
!      
!      
!c      DO echant=1,17
!c         Rmmh(echant)=R_lu(echant)
!c         LAMDAf(echant)=4.1D0*(Rmmh(echant)**(-0.21D0))
!c      enddo
!      
!c      AB=1
!      
!c      do echant=1,AB
!c     Rmmh(echant)=0
!*     initial diameter:Deq=0.1 millimètre
!*     initial radius:AXI=Deq/2
!         
!c         AXI=0.05D0
!         
!c         NDeq(echant,1)=1
!c         Deq_lu(echant,1)=Deq
!c     print*,Deq_lu(echant,1)
!c     calcul de zthéorique
!c         Deq_lu2(echant,1)=Deq_lu(echant,1)**6
!c         somm(echant)=somm(echant)+
!c     &        (NDeq(echant,1)*Deq_lu2(echant,1))
!                           
!c     recalcul de Rmmh ********************
!               
!               
!c     Vit(echant,I)=3.78*Deq_lu(echant,I)**0.67
!c     Volgtte(echant,I)= (4/3)*P*(AXI*1D-3)**3
!               
!c     Rgtte(echant,I)=Volgtte(echant,I)*Vit(echant,I)*
!c     &             (NDeq(echant,I)*1D3)*(deltaD*1D-3)
!c     Rmmh(echant)=Rmmh(echant)+Rgtte(echant,I)*1D3*3600
!               
!               
!******************************************
!c         AXI=AXI+0.05D0
!            
!c         sample(echant,1)=echant
!c         Ztheorique(echant)=10*lOG10(somm(echant))
!c      enddo
!      
!
!*****************************************************************
!*****************début du fichier de sortie**********************
          
      THET0=90D0-Elev
      do param=1,2

         IF(param.EQ.1) then
            THET=THET0
         else 
!c param=2
            THET=180-THET0
         endif
         
!********Initialisation ***************
         ALPHA=0D0
         BETA=0D0
         
         Ntotal=0
         NUMTZ=0
         DENTZ=0
         NUMTD=0
         DENTD=0
         MYS11cr=0D0
         MYS22cr=0D0
         NUMTrhoAB=0D0
!c     ZHHT=0D0
!c     ZVVT=0D0
!         
!c     NDeqtot=0
!c???????
!         
!****************************************
!*********debut de la boucle pour chaque diamètre de goutte mesuré
!*********32 diamètres par echantillon réel, 64 pour echantillon théorique
!         
!C     Diametre equivolumetrique en mètre           
!c         Deq=Deq_lu(echant,1)*1D-3
!            
!         AXI=Deq/2
!         
!c         IF(PARF.EQ."PB70") THEN
!c     détermination de la forme de la goutte (Pruppacher & Beard, 1970)
!c            if (gtte.EQ.2.AND.Deq.GE.5D-4) then               
!c     oblate case
!c               EPS=1./(1.03-62.*Deq)       
!c            ELSE 
!c*     spherical case
!c               EPS=1.
!c            ENDIF               
!c         ELSE
!c     détermination de la forme de la goutte (Andsager et al., 1999)
!c            if (gtte.EQ.2) then               
!c     oblate case
!c               if(Deq.ge.1.1E-3.and.Deq.le.4.4E-3) THEN
!c                  EPS=1./(1.012-14.4*Deq-1.03E4*Deq**2)
!c                  IF(EPS.LT.1.) EPS=1.
!c               ELSE
!c                  EPS=1.0048+.57*Deq-2.628E4*Deq**2+3.682E6*Deq**3
!c     &                 -1.677E8*Deq**4
!c                  IF(EPS.LT.0.2) THEN
!c                     EPS=5.
!c                  ELSE
!c                     EPS=1./EPS
!c                  ENDIF
!c               ENDIF
!c           ELSE 
!c*     spherical case
!c               EPS=1.
!c            ENDIF               
!c         ENDIF
!C         print*,eps
!c     WRITE(*,854) 1
!c     WRITE(*,855) Deq
!c     854  FORMAT('Resultat pour la',I3,' eme goutte')
!c     855  FORMAT('Deq=',F11.7)
!      
!************algorithme NASA*************
      
         DDELT=0.001D0 
!c     NDGS=2
!         
!c     NCHECK=1
!         
!C     IF (NP.EQ.-1) NCHECK=1
!C     IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1
!         
!c     WRITE (10,5454) NCHECK
!c     5454 FORMAT ('  NCHECK=',I1)
!         
!!c     IF(ABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA(EPS,RAT)
!c     IF(ABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
!C     IF(NP.EQ.-3) CALL DROP(RAT)
!         
!c     8000 FORMAT ('RAT=',F8.6)
!         
!C     IF(NP.EQ.-1.AND.EPS.GE.1D0) WRITE(10,7000)EPS
!C     IF(NP.EQ.-1.AND.EPS.LT.1D0) WRITE(10,7001)EPS
!C     IF(NP.GE.0) WRITE(10,7100)NP,EPS
!C     IF(NP.EQ.-3)WRITE(10,7160) 
!C     WRITE(10,7400) LAM,MRR,MRI
!C     WRITE(10,7200) DDELT
         
! 7000    FORMAT('OBLATE SPHEROIDS, A/B=',F11.7)
! 7001    FORMAT('PROLATE SPHEROIDS, A/B=',F11.7)
! 7100    FORMAT('CHEBYSHEV PARTICLES, T',I1,'(',F5.2,')')
! 7160    FORMAT('GENERALIZED CHEBYSHEV PARTICLES')
! 7200    FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
! 7400    FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
         
         DDELT=0.1D0*DDELT
         
!c     IF (ABS(RAT-1D0).LE.1D-6)WRITE(*,8003)AXI 
!c     IF (ABS(RAT-1D0).GT.1D-6) WRITE(*,8004) AXI
!         
!c     8003 FORMAT('EQUAL-VOLUME-SPHERE RADIUS=',F9.7)
!c     8004 FORMAT('EQUAL-SURFACE-AREA-SPHERE RADIUS=',F9.7)
         
      A=AXI
      XEV=2D0*P*A/LAM
      IXXX=XEV+4.05D0*XEV**0.333333D0
      INM1=MAX0(4,IXXX)
      
!C     IF (INM1.GE.NPN1)WRITE(10,7333) NPN1
      IF (INM1.GE.NPN1) STOP
      
! 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  &
!          '.  EXECUTION TERMINATED')
      
      QEXT1=0D0
      QSCA1=0D0
      
!****  Initialisation      
      NMA=INM1
      SORTIE1=.FALSE.
      
      DO WHILE ((NMA.LE.NPN1).AND.(SORTIE1.EQV..FALSE.))
         NMAX=NMA
!c         MMAX=1
         NGAUSS=NMAX*NDGS
         
!C     IF (NGAUSS.GT.NPNG1) WRITE(10,7340) NGAUSS
         IF (NGAUSS.GT.NPNG1) STOP
         
!c 7340    FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',
!c     &        '  EXECUTION TERMINATED')
!c 7334    FORMAT(' NMAX =', I3,'  DSCA=',D8.2,'   DEXT=',D8.2)
         
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS)
         CALL VARY(LAM,MRR,MRI,A,EPS,NGAUSS,X,P,PPI,PIR,PII,R,&
             DR,DDR,DRR,DRI,NMAX)
         CALL TMATR0(NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,&
             DDR,DRR,DRI,NMAX)
         QEXT=0D0
         QSCA=0D0
         
         DO N=1,NMAX
            N1=N+NMAX
            TR1NN=XTR1(N,N)
            TI1NN=XTI1(N,N)
            TR1NN1=XTR1(N1,N1)
            TI1NN1=XTI1(N1,N1)
            DN1=FLOAT(2*N+1)
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+&
                TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
!            TR1NN=TR1(N,N)
!            TI1NN=TI1(N,N)
!            TR1NN1=TR1(N1,N1)
!            TI1NN1=TI1(N1,N1)
!            DN1=FLOAT(2*N+1)
!            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+&
!                TR1NN1*TR1NN1+TI1NN1*TI1NN1)
!            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
         ENDDO
            
         DSCA=ABS((QSCA1-QSCA)/QSCA)
         DEXT=ABS((QEXT1-QEXT)/QEXT)
         QEXT1=QEXT
         QSCA1=QSCA
         
!c     PRINT 7334, NMAX,DSCA,DEXT
         
         IF(.not.(DSCA.LE.DDELT.AND.DEXT.LE.DDELT)) THEN
!C     IF (NMA.EQ.NPN1) WRITE(10,7333) NPN1
            IF (NMA.EQ.NPN1) STOP      
         ELSE
            SORTIE1=.TRUE.
         ENDIF
         
         NMA=NMA+1
         
      ENDDO
      
      NNNGGG=NGAUSS+1
!C      MMAX=NMAX
      
!C     IF (NGAUSS.EQ.NPNG1) WRITE(10,7336) 
      
!****  Initialisation
      SORTIE1=.FALSE.
      NGAUS=NNNGGG
      
      IF (NGAUSS.NE.NPNG1) THEN
         DO WHILE ((SORTIE1.EQV..FALSE.).AND.(NGAUS.LE.NPNG1))
            NGAUSS=NGAUS
!c            NGGG=2*NGAUSS
 7336       FORMAT('WARNING: NGAUSS=NPNG1')
 7337       FORMAT(' NG=',I3,'  DSCA=',D8.2,'   DEXT=',D8.2)
            CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS)
            CALL VARY(LAM,MRR,MRI,A,EPS,NGAUSS,X,P,PPI,PIR,PII,R,&
                DR,DDR,DRR,DRI,NMAX)
            CALL TMATR0(NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,&
                DDR,DRR,DRI,NMAX)
            QEXT=0D0
            QSCA=0D0
            
            DO N=1,NMAX
               N1=N+NMAX
               TR1NN=XTR1(N,N)
               TI1NN=XTI1(N,N)
               TR1NN1=XTR1(N1,N1)
               TI1NN1=XTI1(N1,N1)
               DN1=FLOAT(2*N+1)
               QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+&
                    TR1NN1*TR1NN1+TI1NN1*TI1NN1)
               QEXT=QEXT+(TR1NN+TR1NN1)*DN1
            ENDDO
            
            DSCA=ABS((QSCA1-QSCA)/QSCA)
            DEXT=ABS((QEXT1-QEXT)/QEXT)
            
!c     PRINT 7337, NGGG,DSCA,DEXT
            
            QEXT1=QEXT
            QSCA1=QSCA
            
            IF(.NOT.(DSCA.LE.DDELT.AND.DEXT.LE.DDELT)) THEN
!C     IF (NGAUS.EQ.NPNG1) WRITE(10,7336)
            ELSE
               SORTIE1=.TRUE.
            ENDIF
            
            NGAUS=NGAUS+1
            
         ENDDO
         
      ENDIF
      
      QSCA=0D0
      QEXT=0D0
      NNM=NMAX*2
      
      DO N=1,NNM
         QEXT=QEXT+XTR1(N,N)
      ENDDO
      
      DO N2=1,NMAX
         NN2=N2+NMAX
         DO N1=1,NMAX
            NN1=N1+NMAX
            ZZ1=XTR1(N1,N2)
            XRT11(1,N1,N2)=ZZ1
            ZZ2=XTI1(N1,N2)
            XIT11(1,N1,N2)=ZZ2
            ZZ3=XTR1(N1,NN2)
            XRT12(1,N1,N2)=ZZ3
            ZZ4=XTI1(N1,NN2)
            XIT12(1,N1,N2)=ZZ4
            ZZ5=XTR1(NN1,N2)
            XRT21(1,N1,N2)=ZZ5
            ZZ6=XTI1(NN1,N2)
            XIT21(1,N1,N2)=ZZ6
            ZZ7=XTR1(NN1,NN2)
            XRT22(1,N1,N2)=ZZ7
            ZZ8=XTI1(NN1,NN2)
            XIT22(1,N1,N2)=ZZ8
            QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4+&
                 ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
         ENDDO
      ENDDO
      
!c     PRINT 7800,0,ABS(QEXT),QSCA,NMAX
      
      DO  M=1,NMAX
         CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,&
               DDR,DRR,DRI,NMAX)
         NM=NMAX-M+1
         M1=M+1
         QSC=0D0
         DO N2=1,NM
            NN2=N2+M-1
            N22=N2+NM
            DO N1=1,NM
               NN1=N1+M-1
               N11=N1+NM
               ZZ1=XTR1(N1,N2)
               XRT11(M1,NN1,NN2)=ZZ1
               ZZ2=XTI1(N1,N2)
               XIT11(M1,NN1,NN2)=ZZ2
               ZZ3=XTR1(N1,N22)
               XRT12(M1,NN1,NN2)=ZZ3
               ZZ4=XTI1(N1,N22)
               XIT12(M1,NN1,NN2)=ZZ4
               ZZ5=XTR1(N11,N2)
               XRT21(M1,NN1,NN2)=ZZ5
               ZZ6=XTI1(N11,N2)
               XIT21(M1,NN1,NN2)=ZZ6
               ZZ7=XTR1(N11,N22)
               XRT22(M1,NN1,NN2)=ZZ7
               ZZ8=XTI1(N11,N22)
               XIT22(M1,NN1,NN2)=ZZ8
               QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4+&
                  ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2D0
            ENDDO
         ENDDO
         
         NNM=2*NM
         QXT=0D0
         
         DO N=1,NNM
            QXT=QXT+XTR1(N,N)*2D0
         ENDDO
         
         QSCA=QSCA+QSC
         QEXT=QEXT+QXT
         
!c     PRINT 7800,M,ABS(QXT),QSC,NMAX
! 7800    FORMAT(' m=',I3,'  qxt=',D12.6,'  qsc=',D12.6,&
!            '  nmax=',I3)
         
      ENDDO
      
!c      WALB=-QSCA/QEXT
      
!C     IF (WALB.GT.1D0+DDELT)WRITE(10,9111) 
! 9111 FORMAT ('WARNING: W IS GREATER THAN 1')
!*****************************************************
!C     COMPUTATION OF THE AMPLITUDE AND PHASE MATRICES
!***************************************************
      
!*     Initialisation
      ALPHA=0D0
      PDalpha=1D0/9D0
      Poids=0D0
      PDtotal=0
      S11=0D0
      S12=0D0
      S21=0D0
      S22=0D0
      
      S11carre=0D0
      S22carre=0D0
      NUMrhoAB=0D0
      
      
!c     NDeq(echant,I)=NDeq(echant,I)*1D3
      
      if (param.EQ.1) then
!c     calcul paramètres de propagation
         
!C     PHI0 restera toujours fixe
         PHI0=0D0
         PHI=PHI0
         
!*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         if(oscil.EQ.1) then
!c     algorithme de l'oscillation des gouttes
            
!*     ecart type de beta
            IF (Deq.LE.2D-3) THEN
               SIGBETA=90D0*EXP(-0.95D0*((Deq*1D3)**2))
            ELSE
               SIGBETA=2D0
            ENDIF
            
            PAS=SIGBETA/10D0
!******************************************
            BETA=0D0
            DO WHILE (BETA.LE.(2*SIGBETA))
               Fbeta=EXP(-(BETA**2)/(2*SIGBETA**2))/&
                   (SQRT(2*P*SIGBETA**2))
               Poids=Poids+Fbeta
               BETA=BETA+PAS
            ENDDO
!*****************************************
            
!***********boucle  pour la variation de alpha
            
            DO WHILE (ALPHA.LT.90)
               BETA=0D0
               
!***********boucle    pour la variation de beta
               
               DO WHILE (BETA.LE.(2*SIGBETA))
                  
                  Fbeta=EXP(-(BETA**2)/(2*SIGBETA**2))/&
                  (SQRT(2*P*SIGBETA**2))
                  PDbeta=Fbeta/Poids 
                  
                  
!C     AMPLITUDE MATRIX [Eqs. (2)-(4) of Ref. 6]
                  CALL AMPL(NMAX,LAM,THET0,THET,PHI0,PHI,ALPHA,BETA,&
                      S11u,S12u,S21u,S22u)     
                  
                  S11=S11+S11u*PDalpha*PDbeta
                  S12=S12+S12u*PDalpha*PDbeta
                  S21=S21+S21u*PDalpha*PDbeta
                  S22=S22+S22u*PDalpha*PDbeta
                  
                  
                  BETA=BETA+PAS
                  
                  PDtotal=PDbeta*PDalpha+PDtotal
                  
               ENDDO
               
!*     end de la boucle de variation de beta
               
               ALPHA=ALPHA+10D0
               
            ENDDO
            
!*     end de la boucle de variation de alpha
            
            S11=S11/PDtotal
            S12=S12/PDtotal
            S21=S21/PDtotal
            S22=S22/PDtotal
            
            
!*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c     fin oscillation des gouttes
!c            print*,'youve been here'

         elseif(oscil.EQ.2) then
!C     AMPLITUDE MATRIX [Eqs. (2)-(4) of Ref. 6]
            CALL AMPL(NMAX,LAM,THET0,THET,PHI0,PHI,ALPHA,BETA,&
                S11,S12,S21,S22)    
         endif      
         
!c         DIFFER=DBLE(S22-S11)*LAM**2*1D8          
         
         
!C     Attenuation (H) en dB/km
!C     Ah=4343*LAM*AIMAG(S22)
         
!C         Ah=8686*LAM*AIMAG(S22)
!C         print*,aimag(s22),Deq
         QQEXT=8.*LAM*AIMAG(S22)/P/Deq**2

!C     A reference en db/km
!C     Aref=(0.4343*P**2*Deq**3/LAM)*AIMAG(-K)
!C     LAM et Deq sont en cm
            
!C         Aref=(0.8686*P**2*(Deq*1D2)**3/(LAM*1D2))*IMK
!         
!C     Ah-v=4343*LAM*AIMAG(S22-S11)
            
!C         Ahv=8686*LAM*AIMAG(S22-S11)
            
!C     Changement de phase KDP en degres/km
            
!C         KDP=(180/P)*1D3*LAM*DBLE(S22-S11) 
!C         KDP=Deq*DBLE(S22-S11)
         KDP=LAM**2*DBLE(S22-S11)/2./P            
         
!*************************************************

      elseif(param.EQ.2) then
!c     calcul paramètres de retrodiffusion
         
!C     PHI0 restera toujours fixe
         PHI0=0D0
         PHI=180-PHI0
         
!*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         if(oscil.EQ.1) then
            
!c     algorithme de l'oscillation des gouttes
            
!*     ecart type de beta
            IF (Deq.LE.2D-3) THEN
               SIGBETA=90D0*EXP(-0.95D0*((Deq*1D3)**2))
            ELSE
               SIGBETA=2D0
            ENDIF
            
            PAS=SIGBETA/10D0
!******************************************
            BETA=0D0
            DO WHILE (BETA.LE.(2*SIGBETA))
               Fbeta=EXP(-(BETA**2)/(2*SIGBETA**2))/&
                   (SQRT(2*P*SIGBETA**2))
               Poids=Poids+Fbeta
               BETA=BETA+PAS
            ENDDO
!*****************************************
            
!***********boucle  pour la variation de alpha
               
            DO WHILE (ALPHA.LT.90)
               BETA=0D0
               
!***********boucle    pour la variation de beta
               
               DO WHILE (BETA.LE.(2*SIGBETA))
                  
                  Fbeta=EXP(-(BETA**2)/(2*SIGBETA**2))/&
                      (SQRT(2*P*SIGBETA**2))
                  PDbeta=Fbeta/Poids 
                  
                  
!C     AMPLITUDE MATRIX [Eqs. (2)-(4) of Ref. 6]
                  CALL AMPL(NMAX,LAM,THET0,THET,PHI0,PHI,ALPHA,BETA,&
                      S11u,S12u,S21u,S22u)     
                  
                  S11=S11+S11u*PDalpha*PDbeta
                  S12=S12+S12u*PDalpha*PDbeta
                  S21=S21+S21u*PDalpha*PDbeta
                  S22=S22+S22u*PDalpha*PDbeta
                  
!c                  print*,'NUMrhoA',NUMrhoAB,pdbeta,pdalpha

                  NUMrhoAB=NUMrhoAB+S11u*CONJG(S22u)*PDalpha*PDbeta
                  
                  S11carre=S11carre+((DBLE(S11u))**2+(AIMAG(S11u))**2)*&
                      PDalpha*PDbeta
                  
                  S22carre=S22carre+((DBLE(S22u))**2+(AIMAG(S22u))**2)*&
                      PDalpha*PDbeta
                  
                  BETA=BETA+PAS
!c                  print*,'NUMrhoA',NUMrhoAB,pdbeta,pdalpha
                  
                  PDtotal=PDbeta*PDalpha+PDtotal
                  
               ENDDO
               
!*     end de la boucle de variation de beta
               
               ALPHA=ALPHA+10D0
               
            ENDDO
            
!*     end de la boucle de variation de alpha
            
!c            print*,'pdtotal',pdtotal
            
            S11=S11/PDtotal
            S12=S12/PDtotal
            S21=S21/PDtotal
            S22=S22/PDtotal
            
            S11carre=S11carre/PDtotal
            S22carre=S22carre/PDtotal
            NUMrhoAB=NUMrhoAB/PDtotal
            
!*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c     fin oscillation des gouttes
!c            print*,'youve been there'

         elseif(oscil.EQ.2) then
               
!C     AMPLITUDE MATRIX [Eqs. (2)-(4) of Ref. 6]
            CALL AMPL(NMAX,LAM,THET0,THET,PHI0,PHI,ALPHA,BETA,&
                S11,S12,S21,S22)    
            
            NUMrhoAB=NUMrhoAB+S11*CONJG(S22)
            S11carre=S11carre+((DBLE(S11))**2+(AIMAG(S11))**2)
            S22carre=S22carre+((DBLE(S22))**2+(AIMAG(S22))**2)
         endif
            
!C     ZDR=|Shh|2/|Svv|2=|HH|2/|VV|2
!C     ZDRtotal=Somme(|Shh|2)/Somme(|Svv|2)
         NUMZ=(DBLE(S22))**2+(AIMAG(S22))**2
         DENZ=(DBLE(S11))**2+(AIMAG(S11))**2
!c     ZDR=NUMZ/DENZ
!*************
!c     ZDR=10*LOG10(ZDR)
!*************
!c     PRINT*,'ZDR (DB) =',ZDR
!               
!C     RHO
               
               
!c         DENrho=S11carre*S22carre
         
!c         RHO=NUMrhoAB/(DENrho**0.5)
         
!c         norme=((DBLE(NUMrhoAB))**2+(AIMAG(NUMrhoAB))**2)**0.5
         
!C************
!c         RHOhv=((DBLE(RHO))**2+(AIMAG(RHO))**2)**0.5
         
!C***********


!C     Calcul de ZH et ZV a partir de la matrice d'amplitude
               
!C         ZHH=(LAM**4)*(((DBLE(S22))**2+(AIMAG(S22))**2)*4*P)/(P**5*K2)
         
!C         ZVV=(LAM**4)*(((DBLE(S11))**2+(AIMAG(S11))**2)*4*P)/(P**5*K2)
          
         QBACKHH=16.*ABS(S22)**2/Deq**2
         QBACKVV=16.*ABS(S11)**2/Deq**2
!c         print*,Deq,qback

!c     print *, 10*LOG10(ZHH)+180
          
!C     DELTA EN RADIANS
!****************
         NUMD=AIMAG(S22*CONJG(S11))
         DEND=DBLE(S22*CONJG(S11))  
         DELTA=NUMD/DEND
         
!**************
         DELTA=ATAN(DELTA)
!C     DELTA EN DEGRES
         DELTA=DELTA*180/P
         
         
         NTotal=Ntotal+1
         
         NUMTZ=NUMTZ+NUMZ
         DENTZ=DENTZ+DENZ
         
         
         NUMTD=NUMTD+NUMD
         DENTD=DENTD+DEND
         
!c         ZHHT=ZHHT+ZHH
!c         ZVVT=ZVVT+ZVV
         
         NUMTrhoAB=NUMTrhoAB+NUMrhoAB
         
         MYS11cr=MYS11cr+S11carre
         
         MYS22cr=MYS22cr+S22carre
!***********************************************
            
            
!C     PHASE MATRIX [Eqs. (13)-(29) of Ref. 6]
!            
!c      Z11=0.5D0*(S11*CONJG(S11)+S12*CONJG(S12)
!c     &     +S21*CONJG(S21)+S22*CONJG(S22))
!c      Z12=0.5D0*(S11*CONJG(S11)-S12*CONJG(S12)
!c     &     +S21*CONJG(S21)-S22*CONJG(S22))
!c      Z13=-S11*CONJG(S12)-S22*CONJG(S21)
!c      Z14=(0D0,1D0)*(S11*CONJG(S12)-S22*CONJG(S21))
!c      Z21=0.5D0*(S11*CONJG(S11)+S12*CONJG(S12)
!c     &     -S21*CONJG(S21)-S22*CONJG(S22))
!c      Z22=0.5D0*(S11*CONJG(S11)-S12*CONJG(S12)
!c     &     -S21*CONJG(S21)+S22*CONJG(S22))
!c      Z23=-S11*CONJG(S12)+S22*CONJG(S21)
!c      Z24=(0D0,1D0)*(S11*CONJG(S12)+S22*CONJG(S21))
!c      Z31=-S11*CONJG(S21)-S22*CONJG(S12)
!c      Z32=-S11*CONJG(S21)+S22*CONJG(S12)
!c      Z33=S11*CONJG(S22)+S12*CONJG(S21)
!c      Z34=(0D0,-1D0)*(S11*CONJG(S22)+S21*CONJG(S12))
!c      Z41=(0D0,1D0)*(S21*CONJG(S11)+S22*CONJG(S12))
!c      Z42=(0D0,1D0)*(S21*CONJG(S11)-S22*CONJG(S12))
!c      Z43=(0D0,-1D0)*(S22*CONJG(S11)-S12*CONJG(S21))
!c      Z44=S22*CONJG(S11)-S12*CONJG(S21)
      
!C     WRITE (10,5000) 
!C     WRITE (10,5001) Z11,Z12,Z13,Z14
!C     WRITE (10,5001) Z21,Z22,Z23,Z24
!C     WRITE (10,5001) Z31,Z32,Z33,Z34
!C     WRITE (10,5001) Z41,Z42,Z43,Z44
            
! 5000 FORMAT ('PHASE MATRIX')
! 5001 FORMAT (4F11.7)
            
!****mclock gives the time of execution (CPU)in centieme of second
!c            ITIME=MCLOCK()
!c            TIME=FLOAT(ITIME)/6000D0
!C     WRITE(10,1001)TIME
!c 1001       FORMAT (' time =',F8.2,' min')
            
!*************************************************************
            
!***   enddo de la boucle des gouttes
!                  
!c     ZHHT=ZHHT/NTotal
!c     ZVVT=ZVVT/NTotal
!         
!C     ZH  et ZV en mm6.m-3
!         
!c      ZHHT=ZHHT*1D18
!c      ZVVT=ZVVT*1D18
!         
!C     ZH et ZV en dBZ
!         
!c      ZHHT=10*LOG10(ZHHT)
!c      ZVVT=10*LOG10(ZVVT)
      
      ZDRT= NUMTZ/DENTZ
      ZDRT=10*LOG10(ZDRT)
      
      DELTAT= NUMTD/DENTD
      DELTAT=ATAN(DELTAT)
      DELTAT=DELTAT*180/P
      
      NUMTrhoAB=NUMTrhoAB/NTOTAL
      
      MYS11cr=MYS11cr/NTOTAL
      MYS22cr=MYS22cr/NTOTAL
!c      DENTrhoAB=MYS11cr*MYS22cr
!c      normTOT=((DBLE(NUMTrhoAB))**2+(AIMAG(NUMTrhoAB))**2)**0.5
      
      endif
      
!c      RHOT=NUMTrhoAB/((DENTrhoAB)**0.5)
      
!c      RHOhvT=((DBLE(RHOT))**2+(AIMAG(RHOT))**2)**0.5
      
!c     ecriture du fichier
               
!c     ENDDO
      
!c     CLOSE(8)
      
      ENDDO
      
      RETURN
      END

!************fin du main program***********************************
! 
!C********************************************************************
!C   CALCULATION OF THE AMPLITUDE MATRIX                             *
!C******************************************************************** 
      SUBROUTINE AMPL(NMAX,DLAM,TL,TL1,PL,PL1,ALPHA,BETA,&
                      VV,VH,HV,HH)  

      USE MODD_TMAT,ONLY :XRT11,XRT12,XRT21,XRT22,XIT11,XIT12,&
                          XIT21,XIT22,NPN1,NPN4,NPN6

!c      INCLUDE 'ampld.par.f'
!!      Parameter (NPN1=200,NPN4=NPN1, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-B,D-H,O-Z), COMPLEX*16 (C)
      REAL*8 AL(3,2),AL1(3,2),AP(2,3),AP1(2,3),B(3,3),&
            R(2,2),R1(2,2),C(3,2),CA,CB,CT,CP,CTP,CPP,CT1,CP1,&
            CTP1,CPP1
!C       REAL*8 ZDR,NUM,DEN
      REAL*8 DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
!!      REAL*4 TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),&
!!          TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),&
!!          TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),&
!!          TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CAL(NPN4,NPN4),VV,VH,HV,HH
!!      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22

      IF (ALPHA.LT.0D0.OR.ALPHA.GT.360D0.OR.&
          BETA.LT.0D0.OR.BETA.GT.180D0.OR.&
          TL.LT.0D0.OR.TL.GT.180D0.OR.&
          TL1.LT.0D0.OR.TL1.GT.180D0.OR.&
          PL.LT.0D0.OR.PL.GT.360D0.OR.&
          PL1.LT.0D0.OR.PL1.GT.360D0) THEN 
!C      WRITE (10,2000)
         STOP
      ENDIF  

! 2000 FORMAT ('AN ANGULAR PARAMETER IS OUTSIDE ITS',&
!             ' ALLOWABLE RANGE')
      PIN=ACOS(-1D0)

      PIN2=PIN*0.5D0
      PI=PIN/180D0
      ALPH=ALPHA*PI
      BET=BETA*PI
      THETL=TL*PI
      PHIL=PL*PI
      THETL1=TL1*PI
      PHIL1=PL1*PI

      EPSILON=1D-7
      IF (THETL.LT.PIN2) THETL=THETL+EPSILON
      IF (THETL.GT.PIN2) THETL=THETL-EPSILON
      IF (THETL1.LT.PIN2) THETL1=THETL1+EPSILON
      IF (THETL1.GT.PIN2) THETL1=THETL1-EPSILON
      IF (PHIL.LT.PIN) PHIL=PHIL+EPSILON
      IF (PHIL.GT.PIN) PHIL=PHIL-EPSILON
      IF (PHIL1.LT.PIN) PHIL1=PHIL1+EPSILON
      IF (PHIL1.GT.PIN) PHIL1=PHIL1-EPSILON
      
!C_____________COMPUTE THETP, PHIP, THETP1, AND PHIP1, EQS. (8), (19), AND (20)

      CB=COS(BET)
      SB=SIN(BET)
      CT=COS(THETL)
      ST=SIN(THETL)
      CP=COS(PHIL-ALPH)
      SP=SIN(PHIL-ALPH)
      CTP=CT*CB+ST*SB*CP
      THETP=ACOS(CTP)
      CPP=CB*ST*CP-SB*CT
      SPP=ST*SP
      PHIP=ATAN(SPP/CPP)

      IF (PHIP.GT.0D0.AND.SP.LT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0.AND.SP.GT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0) PHIP=PHIP+2D0*PIN

      CT1=COS(THETL1)
      ST1=SIN(THETL1)
      CP1=COS(PHIL1-ALPH)
      SP1=SIN(PHIL1-ALPH)
      CTP1=CT1*CB+ST1*SB*CP1
      THETP1=ACOS(CTP1)
      CPP1=CB*ST1*CP1-SB*CT1
      SPP1=ST1*SP1
      PHIP1=ATAN(SPP1/CPP1)
      IF (PHIP1.GT.0D0.AND.SP1.LT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0.AND.SP1.GT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0) PHIP1=PHIP1+2D0*PIN

!C____________COMPUTE MATRIX BETA, EQ. (21)

      CA=COS(ALPH)
      SA=SIN(ALPH)
      B(1,1)=CA*CB
      B(1,2)=SA*CB
      B(1,3)=-SB
      B(2,1)=-SA
      B(2,2)=CA
      B(2,3)=0D0
      B(3,1)=CA*SB
      B(3,2)=SA*SB
      B(3,3)=CB

!C____________COMPUTE MATRICES AL AND AL1, EQ. (14) 

      CP=COS(PHIL)
      SP=SIN(PHIL)
      CP1=COS(PHIL1)
      SP1=SIN(PHIL1)
      AL(1,1)=CT*CP
      AL(1,2)=-SP
      AL(2,1)=CT*SP
      AL(2,2)=CP
      AL(3,1)=-ST
      AL(3,2)=0D0
      AL1(1,1)=CT1*CP1
      AL1(1,2)=-SP1
      AL1(2,1)=CT1*SP1
      AL1(2,2)=CP1
      AL1(3,1)=-ST1
      AL1(3,2)=0D0

!C____________COMPUTE MATRICES AP^(-1) AND AP1^(-1), EQ. (15) 

      CT=CTP
      ST=SIN(THETP) 
      CP=COS(PHIP)
      SP=SIN(PHIP)
      CT1=CTP1
      ST1=SIN(THETP1)
      CP1=COS(PHIP1)
      SP1=SIN(PHIP1)
      AP(1,1)=CT*CP
      AP(1,2)=CT*SP
      AP(1,3)=-ST  
      AP(2,1)=-SP
      AP(2,2)=CP 
      AP(2,3)=0D0
      AP1(1,1)=CT1*CP1
      AP1(1,2)=CT1*SP1
      AP1(1,3)=-ST1   
      AP1(2,1)=-SP1
      AP1(2,2)=CP1 
      AP1(2,3)=0D0

!C____________COMPUTE MATRICES R AND R^(-1), EQ. (13)
      DO I=1,3
          DO J=1,2
          X=0D0
              DO K=1,3
              X=X+B(I,K)*AL(K,J)
              ENDDO
          C(I,J)=X
          ENDDO
      ENDDO

      DO I=1,2
          DO J=1,2
          X=0D0
              DO K=1,3
              X=X+AP(I,K)*C(K,J)
              ENDDO
          R(I,J)=X
          ENDDO
      ENDDO

      DO I=1,3
          DO J=1,2
          X=0D0
              DO K=1,3
              X=X+B(I,K)*AL1(K,J)
              ENDDO
          C(I,J)=X
          ENDDO
      ENDDO

      DO I=1,2
          DO J=1,2
          X=0D0
              DO K=1,3
              X=X+AP1(I,K)*C(K,J)
              ENDDO
          R1(I,J)=X
          ENDDO
      ENDDO

      D=1D0/(R1(1,1)*R1(2,2)-R1(1,2)*R1(2,1))
      X=R1(1,1)
      R1(1,1)=R1(2,2)*D
      R1(1,2)=-R1(1,2)*D
      R1(2,1)=-R1(2,1)*D
      R1(2,2)=X*D

      CI=(0D0,1D0)

      DO NN=1,NMAX
          DO N=1,NMAX
          CN=CI**(NN-N-1)
          DNN=FLOAT((2*N+1)*(2*NN+1)) 
          DNN=DNN/FLOAT( N*NN*(N+1)*(NN+1) ) 
          RN=SQRT(DNN)
          CAL(N,NN)=CN*RN
          ENDDO
      ENDDO
    
      DCTH0=CTP
      DCTH=CTP1 
      PH=PHIP1-PHIP
      VV=(0D0,0D0)
      VH=(0D0,0D0)
      HV=(0D0,0D0)
      HH=(0D0,0D0)

      DO M=0,NMAX
         M1=M+1
         NMIN=MAX(M,1)
         
         CALL VIGAMPL(DCTH, NMAX, M, DV1, DV2)
         CALL VIGAMPL(DCTH0, NMAX, M, DV01, DV02)
    
         FC=2D0*COS(M*PH)
         FS=2D0*SIN(M*PH)
         
         DO NN=NMIN,NMAX
            DV1NN=M*DV01(NN)
            DV2NN=DV02(NN)
            DO N=NMIN,NMAX
               DV1N=M*DV1(N)
               DV2N=DV2(N)
               
               CT11=CMPLX(XRT11(M1,N,NN),XIT11(M1,N,NN))
               CT22=CMPLX(XRT22(M1,N,NN),XIT22(M1,N,NN))
               
               IF (M.EQ.0) THEN
                  
                  CN=CAL(N,NN)*DV2N*DV2NN
                  
                  VV=VV+CN*CT22  
                  HH=HH+CN*CT11
                  
               ELSE
                  
                  CT12=CMPLX(XRT12(M1,N,NN),XIT12(M1,N,NN))
                  CT21=CMPLX(XRT21(M1,N,NN),XIT21(M1,N,NN))
                  
                  CN1=CAL(N,NN)*FC
                  CN2=CAL(N,NN)*FS
                  
                  D11=DV1N*DV1NN
                  D12=DV1N*DV2NN
                  D21=DV2N*DV1NN
                  D22=DV2N*DV2NN
                  
                  VV=VV+(CT11*D11+CT21*D21+CT12*D12+CT22*D22)*CN1   
                  
                  VH=VH+(CT11*D12+CT21*D22+CT12*D11+CT22*D21)*CN2
                  
                  HV=HV-(CT11*D21+CT21*D11+CT12*D22+CT22*D12)*CN2   
                  
                  HH=HH+(CT11*D22+CT21*D12+CT12*D21+CT22*D11)*CN1      
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DK=2D0*PIN/DLAM
      VV=VV/DK
      VH=VH/DK
      HV=HV/DK
      HH=HH/DK
      CVV=VV*R(1,1)+VH*R(2,1)
      CVH=VV*R(1,2)+VH*R(2,2)
      CHV=HV*R(1,1)+HH*R(2,1)
      CHH=HV*R(1,2)+HH*R(2,2)
      VV=R1(1,1)*CVV+R1(1,2)*CHV
      VH=R1(1,1)*CVH+R1(1,2)*CHH
      HV=R1(2,1)*CVV+R1(2,2)*CHV
      HH=R1(2,1)*CVH+R1(2,2)*CHH
      
!C      WRITE (10,1005) TL,TL1,PL,PL1,ALPHA,BETA 
!C      WRITE (10,1006)
!C      WRITE (10,1101) VV
!C      WRITE (10,1102) VH
!C      WRITE (10,1103) HV
!C      WRITE (10,1104) HH
! 1101 FORMAT ('S11=',D11.5,' + i*',D11.5)
! 1102 FORMAT ('S12=',D11.5,' + i*',D11.5)
! 1103 FORMAT ('S21=',D11.5,' + i*',D11.5)
! 1104 FORMAT ('S22=',D11.5,' + i*',D11.5)
!
!C     ZDR=|Shh|2/|Svv|2=|HH|2/|VV|2
!C     ZDRtotal=Somme(|Shh|2)/Somme(|Svv|2)
!C        NUM=(REAL(HH))**2+(AIMAG(HH))**2
!C        DEN=(REAL(VV))**2+(AIMAG(VV))**2
!C        ZDR=NUM/DEN
!C       PRINT*,'ZDR=',ZDR
!C        NUMT=NUMT+NUM
!C        DENT=DENT+DEN
!C        ZDRT=NUMT/DENT
        



! 1005 FORMAT ('thet0=',F6.2,'  thet=',F6.2,'  phi0=',F6.2,&
!             '  phi=',F6.2,'  alpha=',F6.2,'  beta=',F6.2)
! 1006 FORMAT ('AMPLITUDE MATRIX')
      RETURN
      END
      
!C*****************************************************************
!C     Calculation of the functions                               *
!C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)                  *
!C     and                                                        *
!C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)                *
!C     1.LE.N.LE.NMAX                                             *
!C     0.LE.X.LE.1                                                *
!C*****************************************************************
      SUBROUTINE VIGAMPL(X, NMAX, M, DV1, DV2)

      USE MODD_TMAT,   ONLY: NPN1,NPN4,NPN6
!!      Parameter (NPN1=200,NPN4=NPN1, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN6), DV2(NPN6)

      DO N=1,NMAX
      DV1(N)=0D0
      DV2(N)=0D0
      ENDDO

      DX=ABS(X)

      IF (ABS(1D0-DX).GT.1D-10) THEN
      A=1D0
      QS=SQRT(1D0-X*X)
      QS1=1D0/QS
      DSI=QS1
           IF (M.EQ.0) THEN
           D1=1D0
           D2=X  
               DO N=1,NMAX  
               QN=FLOAT(N)
               QN1=FLOAT(N+1)
               QN2=FLOAT(2*N+1)
               D3=(QN2*X*D2-QN*D1)/QN1 
               DER=QS1*(QN1*QN/QN2)*(-D1+D3)
               DV1(N)=D2*DSI
               DV2(N)=DER
               D1=D2
               D2=D3
               ENDDO
           ELSE
           QMM=FLOAT(M*M)
               DO I=1,M
               I2=I*2
               A=A*SQRT(FLOAT(I2-1)/FLOAT(I2))*QS
               ENDDO
           D1=0D0
           D2=A 
               DO N=M,NMAX
               QN=FLOAT(N)
               QN2=FLOAT(2*N+1)
               QN1=FLOAT(N+1)
               QNM=SQRT(QN*QN-QMM)
               QNM1=SQRT(QN1*QN1-QMM)
               D3=(QN2*X*D2-QNM*D1)/QNM1
               DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
               DV1(N)=D2*DSI
               DV2(N)=DER
               D1=D2
               D2=D3
               ENDDO
           ENDIF
      ELSE
          IF (M.EQ.1) THEN

              DO N=1,NMAX
              DN=FLOAT(N*(N+1))
              DN=0.5D0*SQRT(DN)
                  IF (X.LT.0D0) DN=DN*(-1)**(N+1)
              DV1(N)=DN
                  IF (X.LT.0D0) DN=-DN
              DV2(N)=DN
              ENDDO
           ENDIF
      ENDIF

      RETURN
      END 

!C**********************************************************************

      SUBROUTINE CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS)

      USE MODD_TMAT,   ONLY: NPN1,NPNG1,NPNG2
      IMPLICIT REAL*8 (A-H,O-Z)
!!      Parameter (NPN1=200, NPNG1=600, NPNG2=2*NPNG1)
      REAL*8 X(NPNG2),W(NPNG2),&
             S(NPNG2),SS(NPNG2),&
             AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)
 
      DO N=1,NMAX
         NN=N*(N+1)
         AN(N)=FLOAT(NN)
         D=SQRT(FLOAT(2*N+1)/FLOAT(NN))
         DD(N)=D
         DO N1=1,N
            DDD=D*DD(N1)*0.5D0
            ANN(N,N1)=DDD
            ANN(N1,N)=DDD
         ENDDO
      ENDDO
      
      NG=2*NGAUSS
      
      CALL GAUSS(NG,0,0,X,W)
      
      DO I=1,NGAUSS
         Y=X(I)
         Y=1D0/(1D0-Y*Y)
         SS(I)=Y
         SS(NG-I+1)=Y
         Y=SQRT(Y)
         S(I)=Y
         S(NG-I+1)=Y
      ENDDO
      
      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE VARY(LAM,MRR,MRI,A,EPS,NGAUSS,X,P,PPI,PIR,PII,&
                     R,DR,DDR,DRR,DRI,NMAX)
       USE MODD_TMAT,   ONLY: NPN1,NPNG1,NPNG2
              
!!      Parameter (NPN1=200, NPNG1=600, NPNG2=2*NPNG1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),R(NPNG2),DR(NPNG2),MRR,MRI,LAM,&
             Z(NPNG2),ZR(NPNG2),ZI(NPNG2),&
             DDR(NPNG2),DRR(NPNG2),DRI(NPNG2)
             
!!      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI

      NG=NGAUSS*2

!c      IF (NP.GT.0) CALL RSP2(X,NG,A,EPS,NP,R,DR)
!C      IF (NP.EQ.-1) 
      CALL RSP1(X,NG,NGAUSS,A,EPS,R,DR)
!c      IF (NP.EQ.-3) CALL RSP4(X,NG,A,R,DR)
!*****PI=2*PI/lambda
      PI=P*2D0/LAM 
!*****PPI=(2*PI/lamda)**2  
      PPI=PI*PI
      PIR=PPI*MRR
      PII=PPI*MRI
      V=1D0/(MRR*MRR+MRI*MRI)
      PRR=MRR*V
      PRI=-MRI*V
      TA=0D0


      DO I=1,NG
         VV=SQRT(R(I))
         V=VV*PI
         TA=MAX(TA,V)
         VV=1D0/V
         DDR(I)=VV
         DRR(I)=PRR*VV
         DRI(I)=PRI*VV
         V1=V*MRR
         V2=V*MRI
         Z(I)=V
         ZR(I)=V1
         ZI(I)=V2
      ENDDO


!C      IF (NMAX.GT.NPN1) WRITE (10,9000) NMAX,NPN1 
      IF (NMAX.GT.NPN1) STOP

 9000 FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
      TB=TA*SQRT(MRR*MRR+MRI*MRI)
      TB=MAX(TB,FLOAT(NMAX))
      NNMAX1=1.2D0*SQRT(MAX(TA,FLOAT(NMAX)))+3D0
      NNMAX2=(TB+4D0*(TB**0.33333D0)+1.2D0*SQRT(TB))
      NNMAX2=NNMAX2-NMAX+5
      CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)

      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE RSP1(X,NG,NGAUSS,REV,EPS,R,DR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)

      A=REV*EPS**(1D0/3D0)
      AA=A*A
      EE=EPS*EPS
      EE1=EE-1D0

      DO I=1,NGAUSS
         C=X(I)
         CC=C*C
         SS=1D0-CC
         S=SQRT(SS)
         RR=1D0/(SS+EE*CC)
         R(I)=AA*RR
         R(NG-I+1)=R(I)
         DR(I)=RR*C*S*EE1
         DR(NG-I+1)=-DR(I)
      ENDDO

      RETURN
      END
 

!C*********************************************************************
 
      SUBROUTINE BESS (X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)

      USE MODD_TMAT,ONLY : XJ,XY,XJR,XJI,XDJ,XDY,XDJR,XDJI,NPN1,NPNG1,NPNG2
!!      Parameter (NPN1=200, NPNG1=600, NPNG2=2*NPNG1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),XR(NG),XI(NG),&
!!             J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),&
!!             JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),&
!!             DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),&
             AJ(NPN1),AY(NPN1),AJR(NPN1),AJI(NPN1),&
             ADJ(NPN1),ADY(NPN1),ADJR(NPN1),&
             ADJI(NPN1)
!!      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
 
      DO I=1,NG
      XX=X(I)

      CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
      CALL RYB(XX,AY,ADY,NMAX)

      YR=XR(I)
      YI=XI(I)

      CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,NNMAX2)

          DO N=1,NMAX
          XJ(I,N)=AJ(N)
          XY(I,N)=AY(N)
          XJR(I,N)=AJR(N)
          XJI(I,N)=AJI(N)
          XDJ(I,N)=ADJ(N)
          XDY(I,N)=ADY(N)
          XDJR(I,N)=ADJR(N)
          XDJI(I,N)=ADJI(N)
          ENDDO
      ENDDO

      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE RJB(X,Y,U,NMAX,NNMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),U(NMAX),Z(800)

      L=NMAX+NNMAX
      XX=1D0/X
      Z(L)=1D0/(FLOAT(2*L+1)*XX)
      L1=L-1

      DO I=1,L1
      I1=L-I
      Z(I1)=1D0/(FLOAT(2*I1+1)*XX-Z(I1+1))
      ENDDO

      Z0=1D0/(XX-Z(1))
      Y0=Z0*COS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1

      DO I=2,NMAX
      YI1=Y(I-1)
      YI=YI1*Z(I)
      U(I)=YI1-FLOAT(I)*YI*XX
      Y(I)=YI
      ENDDO

      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE RYB(X,Y,V,NMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),V(NMAX)

      C=COS(X)
      S=SIN(X)
      X1=1D0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3D0*X3+X1)*C-3D0*X2*S
      NMAX1=NMAX-1

      DO I=2,NMAX1
      Y(I+1)=FLOAT(2*I+1)*X1*Y(I)-Y(I-1)
      ENDDO

      V(1)=-X1*(C+Y1)
      DO  I=2,NMAX
      V(I)=Y(I-1)-FLOAT(I)*X1*Y(I)
      ENDDO

      RETURN
      END
 
!C**********************************************************************
!C                                                                     *
!C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
!C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
!C   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
!C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
!C                                                                     *
!C**********************************************************************
 
      SUBROUTINE CJB (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)

      USE MODD_TMAT,ONLY:NPN1
!!      Parameter (NPN1=200)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 YR(NMAX),YI(NMAX),UR(NMAX),UI(NMAX)
      REAL*8 CYR(NPN1),CYI(NPN1),CZR(1200),CZI(1200)

      L=NMAX+NNMAX
      XRXI=1D0/(XR*XR+XI*XI)
      CXXR=XR*XRXI
      CXXI=-XI*XRXI 
      QF=1D0/FLOAT(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1

      DO I=1,L1
      I1=L-I
      QF=FLOAT(2*I1+1)
      AR=QF*CXXR-CZR(I1+1)
      AI=QF*CXXI-CZI(I1+1)
      ARI=1D0/(AR*AR+AI*AI)
      CZR(I1)=AR*ARI
      CZI(I1)=-AI*ARI
      ENDDO   

      AR=CXXR-CZR(1)
      AI=CXXI-CZI(1)
      ARI=1D0/(AR*AR+AI*AI)
      CZ0R=AR*ARI
      CZ0I=-AI*ARI
      CR=COS(XR)*COSH(XI)
      CI=-SIN(XR)*SINH(XI)
      AR=CZ0R*CR-CZ0I*CI
      AI=CZ0I*CR+CZ0R*CI
      CY0R=AR*CXXR-AI*CXXI
      CY0I=AI*CXXR+AR*CXXI
      CY1R=CY0R*CZR(1)-CY0I*CZI(1)
      CY1I=CY0I*CZR(1)+CY0R*CZI(1)
      AR=CY1R*CXXR-CY1I*CXXI
      AI=CY1I*CXXR+CY1R*CXXI
      CU1R=CY0R-AR
      CU1I=CY0I-AI
      CYR(1)=CY1R
      CYI(1)=CY1I
!c      CUR(1)=CU1R
!c      CUI(1)=CU1I
      YR(1)=CY1R
      YI(1)=CY1I
      UR(1)=CU1R
      UI(1)=CU1I

      DO I=2,NMAX
      QI=FLOAT(I)
      CYI1R=CYR(I-1)
      CYI1I=CYI(I-1)
      CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
      CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
      AR=CYIR*CXXR-CYII*CXXI
      AI=CYII*CXXR+CYIR*CXXI
      CUIR=CYI1R-QI*AR
      CUII=CYI1I-QI*AI
      CYR(I)=CYIR
      CYI(I)=CYII
!c      CUR(I)=CUIR
!c      CUI(I)=CUII
      YR(I)=CYIR
      YI(I)=CYII
      UR(I)=CUIR
      UI(I)=CUII
      ENDDO   

      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE TMATR0(NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,&
                       DRR,DRI,NMAX)

      USE MODD_TMAT,ONLY : XJ,XY,XJR,XJI,XDJ,XDY,XDJR,XDJI,&
                           XQR,XQI,XRGQR,XRGQI,&
                           NPN1,NPNG1,NPNG2,NPN2,NPN4,NPN6
!!      Parameter (NPN1=200, NPNG1=600, NPNG2=2*NPNG1, NPN2=2*NPN1,&   
!!          NPN4=NPN1,NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),&
             R(NPNG2),DR(NPNG2),SIG(NPN2),&
!!             J(NPNG2,NPN1),Y(NPNG2,NPN1),&
!!             JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),&
!!             DY(NPNG2,NPN1),DJR(NPNG2,NPN1),&
!!             DJI(NPNG2,NPN1),
             DDR(NPNG2),DRR(NPNG2),&
             D1(NPNG2,NPN1),D2(NPNG2,NPN1),&
             DRI(NPNG2),RR(NPNG2),&
             DV1(NPN1),DV2(NPN1)
        REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: R11,R12,R21,R22
        REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: I11,I12,I21,I22
        REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: RG11,RG12,RG21,RG22
        REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: IG11,IG12,IG21,IG22

     REAL*8  ANN(NPN1,NPN1),&
!!             QR(NPN2,NPN2),QI(NPN2,NPN2),&
!!             RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),&
             TQR(NPN2,NPN2),TQI(NPN2,NPN2),&
             TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
             
!!      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
!!      COMMON /CTT/ QR,QI,RGQR,RGQI

      MM1=1
      NNMAX=NMAX+NMAX
!c      NG=2*NGAUSS
!c      NGSS=NG
!c      FACTOR=1D0

!c      IF (NCHECK.EQ.1) THEN
      NGSS=NGAUSS
      FACTOR=2D0
!c      ENDIF

      SI=1D0



! Allocation des tableaux :
ALLOCATE(R11(NPN1,NPN1))
ALLOCATE(R12(NPN1,NPN1))
ALLOCATE(R21(NPN1,NPN1))
ALLOCATE(R22(NPN1,NPN1))
ALLOCATE(I11(NPN1,NPN1))
ALLOCATE(I12(NPN1,NPN1))
ALLOCATE(I21(NPN1,NPN1))
ALLOCATE(I22(NPN1,NPN1))
ALLOCATE(RG11(NPN1,NPN1))
ALLOCATE(RG12(NPN1,NPN1))
ALLOCATE(RG21(NPN1,NPN1))
ALLOCATE(RG22(NPN1,NPN1))
ALLOCATE(IG11(NPN1,NPN1))
ALLOCATE(IG12(NPN1,NPN1))
ALLOCATE(IG21(NPN1,NPN1))
ALLOCATE(IG22(NPN1,NPN1))


      DO N=1,NNMAX
      SI=-SI
      SIG(N)=SI
      ENDDO

   20 DO I=1,NGAUSS
      I1=NGAUSS+I
      I2=NGAUSS-I+1
      CALL VIG ( X(I1), NMAX, 0, DV1, DV2)
          DO N=1,NMAX
          SI=SIG(N)
          DD1=DV1(N)
          DD2=DV2(N)
          D1(I1,N)=DD1
          D2(I1,N)=DD2
          D1(I2,N)=DD1*SI
          D2(I2,N)=-DD2*SI
          ENDDO
      ENDDO

   30 DO I=1,NGSS
           RR(I)=W(I)*R(I)
      ENDDO
 
      DO N1=MM1,NMAX
      AN1=AN(N1)
          DO N2=MM1,NMAX
          AN2=AN(N2)
          AR12=0D0
          AR21=0D0
          AI12=0D0
          AI21=0D0
          GR12=0D0
          GR21=0D0
          GI12=0D0
          GI21=0D0
              IF (.NOT.(SIG(N1+N2).LT.0D0)) THEN
                  DO I=1,NGSS
                  D1N1=D1(I,N1)
                  D2N1=D2(I,N1)
                  D1N2=D1(I,N2)
                  D2N2=D2(I,N2)
                  A12=D1N1*D2N2
                  A21=D2N1*D1N2
                  A22=D2N1*D2N2
                  AA1=A12+A21

                  QJ1=XJ(I,N1)
                  QY1=XY(I,N1)
                  QJR2=XJR(I,N2)
                  QJI2=XJI(I,N2)
                  QDJR2=XDJR(I,N2)
                  QDJI2=XDJI(I,N2)
                  QDJ1=XDJ(I,N1)
                  QDY1=XDY(I,N1)
 
                  C1R=QJR2*QJ1
                  C1I=QJI2*QJ1
                  B1R=C1R-QJI2*QY1
                  B1I=C1I+QJR2*QY1
 
                  C2R=QJR2*QDJ1
                  C2I=QJI2*QDJ1
                  B2R=C2R-QJI2*QDY1
                  B2I=C2I+QJR2*QDY1

                  DDRI=DDR(I)
                  C3R=DDRI*C1R
                  C3I=DDRI*C1I
                  B3R=DDRI*B1R
                  B3I=DDRI*B1I

                  C4R=QDJR2*QJ1
                  C4I=QDJI2*QJ1
                  B4R=C4R-QDJI2*QY1
                  B4I=C4I+QDJR2*QY1
  
                  DRRI=DRR(I)
                  DRII=DRI(I)
                  C5R=C1R*DRRI-C1I*DRII
                  C5I=C1I*DRRI+C1R*DRII
                  B5R=B1R*DRRI-B1I*DRII
                  B5I=B1I*DRRI+B1R*DRII
 
                  URI=DR(I)
                  RRI=RR(I)
 
                  F1=RRI*A22
                  F2=RRI*URI*AN1*A12
                  AR12=AR12+F1*B2R+F2*B3R
                  AI12=AI12+F1*B2I+F2*B3I
                  GR12=GR12+F1*C2R+F2*C3R
                  GI12=GI12+F1*C2I+F2*C3I
                    
                  F2=RRI*URI*AN2*A21
                  AR21=AR21+F1*B4R+F2*B5R
                  AI21=AI21+F1*B4I+F2*B5I
                  GR21=GR21+F1*C4R+F2*C5R
                  GI21=GI21+F1*C4I+F2*C5I
                  ENDDO
              ENDIF
          AN12=ANN(N1,N2)*FACTOR
          R12(N1,N2)=AR12*AN12
          R21(N1,N2)=AR21*AN12
          I12(N1,N2)=AI12*AN12
          I21(N1,N2)=AI21*AN12
          RG12(N1,N2)=GR12*AN12
          RG21(N1,N2)=GR21*AN12
          IG12(N1,N2)=GI12*AN12
          IG21(N1,N2)=GI21*AN12
          ENDDO
      ENDDO

      TPIR=PIR
      TPII=PII
      TPPI=PPI
 
      NM=NMAX

      DO N1=MM1,NMAX
      K1=N1-MM1+1
      KK1=K1+NM
          DO N2=MM1,NMAX
          K2=N2-MM1+1
          KK2=K2+NM
 
          TAR12= I12(N1,N2)
          TAI12=-R12(N1,N2)
          TGR12= IG12(N1,N2)
          TGI12=-RG12(N1,N2)
 
          TAR21=-I21(N1,N2)
          TAI21= R21(N1,N2)
          TGR21=-IG21(N1,N2)
          TGI21= RG21(N1,N2)
 
          TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
          TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
          TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
          TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
          TQR(K1,KK2)=0D0
          TQI(K1,KK2)=0D0
          TRGQR(K1,KK2)=0D0
          TRGQI(K1,KK2)=0D0
 
          TQR(KK1,K2)=0D0
          TQI(KK1,K2)=0D0
          TRGQR(KK1,K2)=0D0
          TRGQI(KK1,K2)=0D0

          TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
          TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
          TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
          TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
          ENDDO
      ENDDO
 
      NNMAX=2*NM

      DO N1=1,NNMAX
          DO N2=1,NNMAX
          XQR(N1,N2)=TQR(N1,N2)
          XQI(N1,N2)=TQI(N1,N2)
          XRGQR(N1,N2)=TRGQR(N1,N2)
          XRGQI(N1,N2)=TRGQI(N1,N2)
          ENDDO
      ENDDO
DEALLOCATE(R12)
DEALLOCATE(R21)
DEALLOCATE(R22)
DEALLOCATE(I11)
DEALLOCATE(I12)
DEALLOCATE(I21)
DEALLOCATE(I22)
DEALLOCATE(RG11)
DEALLOCATE(RG12)
DEALLOCATE(RG21)
DEALLOCATE(RG22)
DEALLOCATE(IG11)
DEALLOCATE(IG12)
DEALLOCATE(IG21)
DEALLOCATE(IG22)

      CALL TT(NMAX)

      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,&
                      DRR,DRI,NMAX)

      USE MODD_TMAT,ONLY : XJ,XY,XJR,XJI,XDJ,XDY,XDJR,XDJI,&
                           XQR,XQI,XRGQR,XRGQI,&
                           NPN1,NPNG1,NPNG2,NPN2,NPN4,NPN6

!c      INCLUDE 'ampld.par.f'
!!      Parameter (NPN1=200, NPNG1=600, NPNG2=2*NPNG1, NPN2=2*NPN1, & 
!!           NPN4=NPN1, NPN6=NPN4+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),&
             R(NPNG2),DR(NPNG2),SIG(NPN2),&
!!             J(NPNG2,NPN1),Y(NPNG2,NPN1),&
!!             JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),&
!!             DY(NPNG2,NPN1),DJR(NPNG2,NPN1),&
!!             DJI(NPNG2,NPN1),DDR(NPNG2),
             DRR(NPNG2),&
             D1(NPNG2,NPN1),D2(NPNG2,NPN1),&
             DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),&
             DV1(NPN1),DV2(NPN1)
      REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: R11,R12,R21,R22
      REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: I11,I12,I21,I22
      REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: RG11,RG12,RG21,RG22
      REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: IG11,IG12,IG21,IG22

      REAL*8  ANN(NPN1,NPN1),&
!!             QR(NPN2,NPN2),QI(NPN2,NPN2),&
!!             RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),&
             TQR(NPN2,NPN2),TQI(NPN2,NPN2),&
             TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
!!      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
!!      COMMON /CTT/ QR,QI,RGQR,RGQI

      MM1=M
      QM=FLOAT(M)
      QMM=QM*QM
!c      NG=2*NGAUSS
!c      NGSS=NG
!c      FACTOR=1D0

!c      IF (NCHECK.EQ.1) THEN
         NGSS=NGAUSS
         FACTOR=2D0
!c      ENDIF

      SI=1D0
      NM=NMAX+NMAX

! Allocation des tableaux :
ALLOCATE(R11(NPN1,NPN1))
ALLOCATE(R12(NPN1,NPN1))
ALLOCATE(R21(NPN1,NPN1))
ALLOCATE(R22(NPN1,NPN1))
ALLOCATE(I11(NPN1,NPN1))
ALLOCATE(I12(NPN1,NPN1))
ALLOCATE(I21(NPN1,NPN1))
ALLOCATE(I22(NPN1,NPN1))
ALLOCATE(RG11(NPN1,NPN1))
ALLOCATE(RG12(NPN1,NPN1))
ALLOCATE(RG21(NPN1,NPN1))
ALLOCATE(RG22(NPN1,NPN1))
ALLOCATE(IG11(NPN1,NPN1))
ALLOCATE(IG12(NPN1,NPN1))
ALLOCATE(IG21(NPN1,NPN1))
ALLOCATE(IG22(NPN1,NPN1))



      DO N=1,NM
         SI=-SI
         SIG(N)=SI
      ENDDO

   20 DO I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG(X(I1),NMAX,M,DV1,DV2)
         DO N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
         ENDDO
      ENDDO

 30   DO I=1,NGSS
         WR=W(I)*R(I)
         DS(I)=S(I)*QM*WR
         DSS(I)=SS(I)*QMM
         RR(I)=WR
      ENDDO
 
      DO N1=MM1,NMAX
         AN1=AN(N1)
         DO N2=MM1,NMAX
            AN2=AN(N2)
            AR11=0D0
            AR12=0D0
            AR21=0D0
            AR22=0D0
            AI11=0D0
            AI12=0D0
            AI21=0D0
            AI22=0D0
            GR11=0D0
            GR12=0D0
            GR21=0D0
            GR22=0D0
            GI11=0D0
            GI12=0D0
            GI21=0D0
            GI22=0D0
            SI=SIG(N1+N2)
            DO  I=1,NGSS
               D1N1=D1(I,N1)
               D2N1=D2(I,N1)
               D1N2=D1(I,N2)
               D2N2=D2(I,N2)
               A11=D1N1*D1N2
               A12=D1N1*D2N2
               A21=D2N1*D1N2
               A22=D2N1*D2N2
               AA1=A12+A21
               AA2=A11*DSS(I)+A22
               QJ1=XJ(I,N1)
               QY1=XY(I,N1)
               QJR2=XJR(I,N2)
               QJI2=XJI(I,N2)
               QDJR2=XDJR(I,N2)
               QDJI2=XDJI(I,N2)
               QDJ1=XDJ(I,N1)
               QDY1=XDY(I,N1)

               C1R=QJR2*QJ1
               C1I=QJI2*QJ1
               B1R=C1R-QJI2*QY1
               B1I=C1I+QJR2*QY1
               
               C2R=QJR2*QDJ1
               C2I=QJI2*QDJ1
               B2R=C2R-QJI2*QDY1
               B2I=C2I+QJR2*QDY1
               
               DDRI=DDR(I)
               C3R=DDRI*C1R
               C3I=DDRI*C1I
               B3R=DDRI*B1R
               B3I=DDRI*B1I
 
               C4R=QDJR2*QJ1
               C4I=QDJI2*QJ1
               B4R=C4R-QDJI2*QY1
               B4I=C4I+QDJR2*QY1
               
               DRRI=DRR(I)
               DRII=DRI(I)
               C5R=C1R*DRRI-C1I*DRII
               C5I=C1I*DRRI+C1R*DRII
               B5R=B1R*DRRI-B1I*DRII
               B5I=B1I*DRRI+B1R*DRII
               
               C6R=QDJR2*QDJ1
               C6I=QDJI2*QDJ1
               B6R=C6R-QDJI2*QDY1
               B6I=C6I+QDJR2*QDY1
               
               C7R=C4R*DDRI
               C7I=C4I*DDRI
               B7R=B4R*DDRI
               B7I=B4I*DDRI
               
               C8R=C2R*DRRI-C2I*DRII
               C8I=C2I*DRRI+C2R*DRII
               B8R=B2R*DRRI-B2I*DRII
               B8I=B2I*DRRI+B2R*DRII
               
               URI=DR(I)
               DSI=DS(I)
!c               DSSI=DSS(I)
               RRI=RR(I)
               IF (SI.GT.0D0) THEN    
                  F1=RRI*AA2
                  F2=RRI*URI*AN1*A12
                  AR12=AR12+F1*B2R+F2*B3R
                  AI12=AI12+F1*B2I+F2*B3I
                  GR12=GR12+F1*C2R+F2*C3R
                  GI12=GI12+F1*C2I+F2*C3I
                  
                  F2=RRI*URI*AN2*A21
                  AR21=AR21+F1*B4R+F2*B5R
                  AI21=AI21+F1*B4I+F2*B5I
                  GR21=GR21+F1*C4R+F2*C5R
                  GI21=GI21+F1*C4I+F2*C5I
!c                  IF (NCHECK.NE.1) THEN
!c                     E2=DSI*URI*A11
!c                     E3=E2*AN2
!c                     E2=E2*AN1
!c                     AR22=AR22+E1*B6R+E2*B7R+E3*B8R
!c                     AI22=AI22+E1*B6I+E2*B7I+E3*B8I
!c                     GR22=GR22+E1*C6R+E2*C7R+E3*C8R
!c                     GI22=GI22+E1*C6I+E2*C7I+E3*C8I
!c                  ENDIF
               ELSE
!c                  print*,'aa1',aa1
                  E1=DSI*AA1
                  AR11=AR11+E1*B1R
                  AI11=AI11+E1*B1I
                  GR11=GR11+E1*C1R
                  GI11=GI11+E1*C1I
!c                  IF (NCHECK.EQ.1) THEN
                     E2=DSI*URI*A11
                     E3=E2*AN2
                     E2=E2*AN1
                     AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                     AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                     GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                     GI22=GI22+E1*C6I+E2*C7I+E3*C8I
!c                  ENDIF
               ENDIF
            ENDDO
            
            AN12=ANN(N1,N2)*FACTOR
            R11(N1,N2)=AR11*AN12
            R12(N1,N2)=AR12*AN12
            R21(N1,N2)=AR21*AN12
            R22(N1,N2)=AR22*AN12
            I11(N1,N2)=AI11*AN12
            I12(N1,N2)=AI12*AN12
            I21(N1,N2)=AI21*AN12
            I22(N1,N2)=AI22*AN12
            RG11(N1,N2)=GR11*AN12
            RG12(N1,N2)=GR12*AN12
            RG21(N1,N2)=GR21*AN12
            RG22(N1,N2)=GR22*AN12
            IG11(N1,N2)=GI11*AN12
            IG12(N1,N2)=GI12*AN12
            IG21(N1,N2)=GI21*AN12
            IG22(N1,N2)=GI22*AN12
            
         ENDDO
      ENDDO
      
      TPIR=PIR
      TPII=PII
      TPPI=PPI
      NM=NMAX-MM1+1

      DO N1=MM1,NMAX
         K1=N1-MM1+1
         KK1=K1+NM
         DO N2=MM1,NMAX
            K2=N2-MM1+1
            KK2=K2+NM
            
            TAR11=-R11(N1,N2)
            TAI11=-I11(N1,N2)
            TGR11=-RG11(N1,N2)
            TGI11=-IG11(N1,N2)
            
            TAR12= I12(N1,N2)
            TAI12=-R12(N1,N2)
            TGR12= IG12(N1,N2)
            TGI12=-RG12(N1,N2)
            
            TAR21=-I21(N1,N2)
            TAI21= R21(N1,N2)
            TGR21=-IG21(N1,N2)
            TGI21= RG21(N1,N2)
            
            TAR22=-R22(N1,N2)
            TAI22=-I22(N1,N2)
            TGR22=-RG22(N1,N2)
            TGI22=-IG22(N1,N2)
            
            TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
            TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
            TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
            TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
            
            TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
            TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
            TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
            TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22
            
            TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
            TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
            TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
            TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11
            
            TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
            TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
            TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
            TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
         ENDDO
      ENDDO
      
      NNMAX=2*NM
      
      DO N1=1,NNMAX
         DO N2=1,NNMAX
            XQR(N1,N2)=TQR(N1,N2)
            XQI(N1,N2)=TQI(N1,N2)
            XRGQR(N1,N2)=TRGQR(N1,N2)
            XRGQI(N1,N2)=TRGQI(N1,N2)
         ENDDO
      ENDDO
DEALLOCATE(R12)
DEALLOCATE(R21)
DEALLOCATE(R22)
DEALLOCATE(I11)
DEALLOCATE(I12)
DEALLOCATE(I21)
DEALLOCATE(I22)
DEALLOCATE(RG11)
DEALLOCATE(RG12)
DEALLOCATE(RG21)
DEALLOCATE(RG22)
DEALLOCATE(IG11)
DEALLOCATE(IG12)
DEALLOCATE(IG21)
DEALLOCATE(IG22)

      
      CALL TT(NM)

 
      RETURN
      END
 
!C*****************************************************************
 
      SUBROUTINE VIG(X,NMAX,M,DV1,DV2)

      USE MODD_TMAT, ONLY:NPN1
!!      Parameter (NPN1=200)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN1),DV2(NPN1)
 
      A=1D0
      QS=SQRT(1D0-X*X)
      QS1=1D0/QS
      DO N=1,NMAX
      DV1(N)=0D0
      DV2(N)=0D0
      ENDDO   

      IF (M.NE.0) THEN
      QMM=FLOAT(M*M)
          DO I=1,M
          I2=I*2
          A=A*SQRT(FLOAT(I2-1)/FLOAT(I2))*QS
          ENDDO   
      D1=0D0
      D2=A
          DO N=M,NMAX
          QN=FLOAT(N)
          QN2=FLOAT(2*N+1)
          QN1=FLOAT(N+1)
          QNM=SQRT(QN*QN-QMM)
          QNM1=SQRT(QN1*QN1-QMM)
          D3=(QN2*X*D2-QNM*D1)/QNM1
          DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
          DV1(N)=D2
          DV2(N)=DER
          D1=D2
          D2=D3
          ENDDO   
      ELSE
      D1=1D0
      D2=X  
          DO N=1,NMAX  
          QN=FLOAT(N)
          QN1=FLOAT(N+1)
          QN2=FLOAT(2*N+1)
          D3=(QN2*X*D2-QN*D1)/QN1 
          DER=QS1*(QN1*QN/QN2)*(-D1+D3)
          DV1(N)=D2
          DV2(N)=DER
          D1=D2
          D2=D3
          ENDDO
      ENDIF

      RETURN
      END 
 
!C**********************************************************************
!C                                                                     *
!C   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              *
!C                                                                     *
!C   INPUT INFORTMATION IS IN COMMON /CTT/                             *
!C   OUTPUT INFORMATION IS IN COMMON /CT/                              *
!C                                                                     *
!C**********************************************************************
 
      SUBROUTINE TT(NMAX)
      USE MODD_TMAT,ONLY : XTR1,XTI1,XQR,XQI,XRGQR,XRGQI,NPN1,NPN2


!!      Parameter (NPN1=200, NPN2=2*NPN1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(NPN2,NPN2),&
!!            QR(NPN2,NPN2),QI(NPN2,NPN2),&
!!            RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),&
            A(NPN2,NPN2),C(NPN2,NPN2),D(NPN2,NPN2),E(NPN2,NPN2)
!!      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
!c      INTEGER IPVT(NPN2)
!!      COMMON /CT/ TR1,TI1
!!      COMMON /CTT/ QR,QI,RGQR,RGQI

       NDIM=NPN2
       NNMAX=2*NMAX
 
!C        Gaussian elimination

          DO N1=1,NNMAX
              DO N2=1,NNMAX
              F(N1,N2)=XQI(N1,N2)
              ENDDO
          ENDDO

!c          IF (NCHECK.EQ.1) THEN
          CALL INV1(NMAX,F,A)
!c          ELSE
!c          CALL INVERT(NDIM,NNMAX,F,A,COND,IPVT,WORK,B) 
!c          ENDIF
 
       CALL PROD(XQR,A,C,NDIM,NNMAX)
       CALL PROD(C,XQR,D,NDIM,NNMAX)

          DO N1=1,NNMAX
              DO N2=1,NNMAX
              C(N1,N2)=D(N1,N2)+XQI(N1,N2)
              ENDDO
          ENDDO

!c          IF (NCHECK.EQ.1) THEN
          CALL INV1(NMAX,C,XQI)
!c          ELSE
!c          CALL INVERT(NDIM,NNMAX,C,QI,COND,IPVT,WORK,B) 
!c          ENDIF

       CALL PROD(A,XQR,D,NDIM,NNMAX)
       CALL PROD(D,XQI,XQR,NDIM,NNMAX)
 
       CALL PROD(XRGQR,XQR,A,NDIM,NNMAX)
       CALL PROD(XRGQI,XQI,C,NDIM,NNMAX)
       CALL PROD(XRGQR,XQI,D,NDIM,NNMAX)
       CALL PROD(XRGQI,XQR,E,NDIM,NNMAX)

          DO N1=1,NNMAX
              DO N2=1,NNMAX
              XTR1(N1,N2)=-A(N1,N2)-C(N1,N2)
              XTI1(N1,N2)= D(N1,N2)-E(N1,N2)
              ENDDO
          ENDDO
           RETURN
           END
 
!C********************************************************************
!C     Calcul of produit de 2 matrices de dimension NDIM*N           *
!C********************************************************************
      SUBROUTINE PROD(A,B,C,NDIM,N)

      REAL*8 A(NDIM,N),B(NDIM,N),C(NDIM,N),cij

      DO I=1,N
          DO J=1,N
          CIJ=0d0
              DO K=1,N
              CIJ=CIJ+A(I,K)*B(K,J)
              ENDDO
          C(I,J)=CIJ
          ENDDO
      ENDDO

      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE INV1(NMAX,F,A)

      USE MODD_TMAT,ONLY : NPN1,NPN2
      IMPLICIT REAL*8 (A-H,O-Z)
!!      Parameter (NPN1=200, NPN2=2*NPN1)
      REAL*8  A(NPN2,NPN2),F(NPN2,NPN2),B(NPN1),&
             WORK(NPN1),Q1(NPN1,NPN1),Q2(NPN1,NPN1),&
             P1(NPN1,NPN1),P2(NPN1,NPN1)
      INTEGER IPVT(NPN1),IND1(NPN1),IND2(NPN1)

      NDIM=NPN1
      NN1=(FLOAT(NMAX)-0.1D0)*0.5D0+1D0 
      NN2=NMAX-NN1
      DO I=1,NMAX
      IND1(I)=2*I-1
          IF(I.GT.NN1) IND1(I)=NMAX+2*(I-NN1)
      IND2(I)=2*I
          IF(I.GT.NN2) IND2(I)=NMAX+2*(I-NN2)-1
      ENDDO

      NNMAX=2*NMAX

      DO I=1,NMAX
      I1=IND1(I)
      I2=IND2(I)
          DO J=1,NMAX
          J1=IND1(J)
          J2=IND2(J)
          Q1(J,I)=F(J1,I1)
          Q2(J,I)=F(J2,I2)
          ENDDO
      ENDDO

      CALL INVERT(NDIM,NMAX,Q1,P1,COND,IPVT,WORK,B)
      CALL INVERT(NDIM,NMAX,Q2,P2,COND,IPVT,WORK,B)

      DO I=1,NNMAX
          DO J=1,NNMAX
          A(J,I)=0D0
          ENDDO
      ENDDO

      DO I=1,NMAX
      I1=IND1(I)
      I2=IND2(I)
          DO J=1,NMAX
          J1=IND1(J)
          J2=IND2(J)
          A(J1,I1)=P1(J,I)
          A(J2,I2)=P2(J,I)
          ENDDO
      ENDDO

      RETURN
      END
 
!C*********************************************************************
 
      SUBROUTINE INVERT(NDIM,N,A,X,COND,IPVT,WORK,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),X(NDIM,N),WORK(N),B(N)
      INTEGER IPVT(N)

      CALL DECOMP (NDIM,N,A,COND,IPVT,WORK)

!C     IF (COND+1D0.EQ.COND) WRITE(10,5) COND

!C     IF (COND+1D0.EQ.COND) STOP

!C    5 FORMAT('THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY',&
!C           'COND = ',D12.6)

      DO I=1,N
          DO J=1,N
          B(J)=0D0
              IF (J.EQ.I) B(J)=1D0
          ENDDO
      CALL SOLVE (NDIM,N,A,B,IPVT)
           DO J=1,N
           X(J,I)=B(J)
           ENDDO
      ENDDO
   
      RETURN
      END
 
!C********************************************************************
 
      SUBROUTINE DECOMP (NDIM,N,A,COND,IPVT,WORK)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)

      IPVT(N)=1
!*****IF 1
      IF(N.NE.1) THEN
      NM1=N-1
      ANORM=0D0
          DO J=1,N
          T=0D0
               DO I=1,N
               T=T+ABS(A(I,J))
               ENDDO
               IF (T.GT.ANORM) ANORM=T
          ENDDO

          DO K=1,NM1
          KP1=K+1
          M=K
              DO I=KP1,N
                   IF (ABS(A(I,K)).GT.ABS(A(M,K))) M=I
              ENDDO
          IPVT(K)=M
              IF (M.NE.K) IPVT(N)=-IPVT(N)
          T=A(M,K)
          A(M,K)=A(K,K)
          A(K,K)=T
              IF (T.NE.0d0) THEN
                  DO I=KP1,N
                  A(I,K)=-A(I,K)/T
                  ENDDO
                  DO J=KP1,N
                  T=A(M,J)
                  A(M,J)=A(K,J)
                  A(K,J)=T
                      IF (T.NE.0D0) THEN
                          DO I=KP1,N
                          A(I,J)=A(I,J)+A(I,K)*T
                          ENDDO
                      ENDIF
                 ENDDO
             ENDIF
         ENDDO

      DO K=1,N
      T=0D0
          IF (K.NE.1) THEN
          KM1=K-1
              DO I=1,KM1
              T=T+A(I,K)*WORK(I)
              ENDDO
          ELSE
          EK=1D0
              IF (T.LT.0D0) EK=-1D0
              IF (A(K,K).EQ.0D0) THEN 
!               COND=1D52
              RETURN
              ELSE
              WORK(K)=-(EK+T)/A(K,K)
              ENDIF
          ENDIF
      ENDDO

      DO KB=1,NM1
      K=N-KB
      T=0D0
      KP1=K+1
          DO I=KP1,N
          T=T+A(I,K)*WORK(K)
          ENDDO
      WORK(K)=T
      M=IPVT(K)
          IF (M.EQ.K) THEN
          T=WORK(M)
          WORK(M)=WORK(K)
          WORK(K)=T
          ENDIF
      ENDDO
      YNORM=0D0

      DO I=1,N
      YNORM=YNORM+ABS(WORK(I))
      ENDDO
  
      CALL SOLVE (NDIM,N,A,WORK,IPVT)
      ZNORM=0D0

      DO I=1,N
      ZNORM=ZNORM+ABS(WORK(I))
      ENDDO

      COND=ANORM*ZNORM/YNORM
      IF (COND.LT.1d0) COND=1D0
      RETURN

      ELSE
!****ELSE DU IF 1
      COND=1D0
          IF (A(1,1).EQ.0D0) THEN
          COND=1D52
          ENDIF
      ENDIF
!****END DU IF 1   

      RETURN
      END
 
!C**********************************************************************
 
      SUBROUTINE SOLVE (NDIM,N,A,B,IPVT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),B(N)
      INTEGER IPVT(N)

      IF (N.NE.1) THEN
      NM1=N-1
          DO K=1,NM1
          KP1=K+1
          M=IPVT(K)
          T=B(M)
          B(M)=B(K)
          B(K)=T
              DO I=KP1,N
              B(I)=B(I)+A(I,K)*T
              ENDDO
          ENDDO
          DO KB=1,NM1
          KM1=N-KB
          K=KM1+1
          B(K)=B(K)/A(K,K)
          T=-B(K)
              DO I=1,KM1
              B(I)=B(I)+A(I,K)*T
              ENDDO
          ENDDO
      ENDIF

      B(1)= B(1)/A(1,1)

      RETURN
      END
 
!C**********************************************************************
!C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
!C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
!C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
!C    N - NUMBER OF POINTS                                             *
!C    Z - DIVISION POINTS                                              *
!C    W - WEIGHTS                                                      *
!C**********************************************************************
 
      SUBROUTINE GAUSS(N,IND1,IND2,Z,W)
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      DATA A,B,C /1D0,2D0,3D0/
    
      IND=MOD(N,2)
      K=N/2+IND
      F=FLOAT(N)

!*****DO 1
      DO I=1,K
      M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
      NITER=0
      CHECK=1D-16
      PB=1D17

          DO WHILE(ABS(PB).GT.CHECK*ABS(X))
          PB=1D0
          NITER=NITER+1
              IF (NITER.GT.100) CHECK=CHECK*10D0
              PC=X
              DJ=A
              DO J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
              PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
              ENDDO
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          ENDDO  

      Z(M)=X
      W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M) 
          IF(.NOT.(I.EQ.K.AND.IND.EQ.1))THEN
          Z(I)=-Z(M)
          W(I)=W(M)
          ENDIF

      ENDDO
!****END DU DO 1
!!      IF(IND2.EQ.1) THEN
!!      
!! 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',&
!!      ' OF ',I4,'-TH ORDER')
!!          DO I=1,K
!!c          ZZ=-Z(I)
!!C          WRITE(10,1200) I,ZZ,I,W(I)
!!          ENDDO
!! 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
!!      ELSE
!!C     PRINT 1300,N
!! 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
!!      ENDIF

      IF(IND1.NE.0) THEN
          DO  I=1,N
          Z(I)=(A+Z(I))/B
          ENDDO 
      ENDIF

      RETURN
      END








