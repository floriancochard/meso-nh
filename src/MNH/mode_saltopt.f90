!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!        ###################
         MODULE MODE_SALTOPT
!        ###################
!
!!
!!    PURPOSE
!!    -------
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!
  IMPLICIT NONE
  PUBLIC
  PRIVATE :: SALTOPT_LKT

CONTAINS
  
  !****************************************************
  SUBROUTINE SALTOPT_GET(     &
       PSVT                   & !I [moments/molec_{air}] Transported moments of sea salts
       ,PZZ                   & !I [m] height of layers
       ,PRHODREF              & !I [kg/m3] density of air
       ,PTHT                  &
       ,PPABST                &
       ,PRT                   &
       ,PPIZA_WVL             & !O [-] single scattering albedo of sea salt layer for all SW wavelengths
       ,PCGA_WVL              & !O [-] assymetry factor for sea salt layer for all SW wavelengths
       ,PTAUREL_WVL           & !O [-] opt.depth/opt.depth(550) for sea salt layer for all SW wvl 
       ,PTAU550               & !O [-] opt.depth at 550nm for all sea salt layer
       ,KSWB                  & !I [nbr] number of shortwave bands
       )
    

    USE MODE_SALT_PSD_WET   !Conversion procedures from moments to radius, ,number, mass and sigma
    USE MODE_SALT_PSD
    USE MODD_SALT, ONLY : NMODE_SLT

    IMPLICIT NONE
    
    !INPUT
    REAL, DIMENSION(:,:,:,:),INTENT(IN)      :: PSVT       !I [moments/molec_{air}] transported moments of sea salt
    REAL, DIMENSION(:,:,:),INTENT(IN)        :: PZZ        !I [m] height of layers
    REAL, DIMENSION(:,:,:),INTENT(IN)        :: PRHODREF   !I [kg/m3] density of air
    REAL, DIMENSION(:,:,:),INTENT(IN)        :: PTHT, PPABST   !I
    REAL, DIMENSION(:,:,:,:),INTENT(IN)      :: PRT
    INTEGER, INTENT(IN)                      :: KSWB       !I [nbr] number of shortwave wavelengths
    REAL, PARAMETER                          :: EPSILON=1.e-8 !a very low number for optical depth in a layer

    !OUTPUT
    REAL, DIMENSION(:,:,:,:),INTENT(OUT)     :: PPIZA_WVL   !O [-] single scattering albedo of sea salt layer for all SW wavelengths
    REAL, DIMENSION(:,:,:,:),INTENT(OUT)     :: PCGA_WVL    !O [-] assymetry factor for sea salt layer for all SW wavelengths
    REAL, DIMENSION(:,:,:,:),INTENT(OUT)     :: PTAUREL_WVL !O [-] opt.depth/opt.depth(550) for sea salt layer for all SW wvl 
    REAL, DIMENSION(:,:,:), INTENT(OUT)      :: PTAU550     !O [-] opt.depth at 550nm for all sea salt layer

    !LOCAL VARIABLES
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),SIZE(PSVT,4))  ::    ZSVT
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NMODE_SLT) :: ZMASS         ![kg/m3] mass of one sea salt mode
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NMODE_SLT) :: ZRADIUS       ![um] number median radius of one sea salt mode
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NMODE_SLT) :: ZSIGMA        ![-] dispersion coefficient one sea salt mode
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NMODE_SLT) :: ZDENSITY        ![-] [g/m2] density of wet aerosol (water + sea salt)
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)                   :: ZTAU550_MDE   ![-] opt.depth 550nm one mode
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                 :: ZTAU_WVL_MDE  ![-] opt.depth @ wvl, one mode
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                 :: ZPIZA_WVL_MDE ![-] single scattering albedo @ wvl, one mode
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                 :: ZCGA_WVL_MDE  ![-] assymetry factor @ wvl, one mode
    INTEGER                                                 :: JMDE        ![idx] counter for modes
    INTEGER                                                 :: JWVL        ![idx] counter for wavelengths
    !Allocate arrays which size depend on number of modes
    ALLOCATE(ZTAU550_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),NMODE_SLT)) 
    ALLOCATE(ZTAU_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_SLT))
    ALLOCATE(ZPIZA_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_SLT))
    ALLOCATE(ZCGA_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_SLT))
    
    ZSVT(:,:,:,:)=PSVT(:,:,:,:)
    CALL PPP2SALT_WET(     &
         ZSVT                                   & !I [moments/molec_{air}] moments of sea salt for all modes
         ,PRHODREF                              & !I [kg/m3] air density
         ,PPABST                                & !I Pression
         ,PTHT                                  & !I Potential temperature
         ,PRT                                   & !I Large scale vapor mixing ratio
         ,PSIG3D=ZSIGMA                         & !O [-] dispersion coefficient
         ,PRG3D=ZRADIUS                         & !O [um] number median radius
         ,PMASS3D=ZMASS                         & !O [kg/m3] mass of sea salt
         ,PDENSITY_WET=ZDENSITY                 & !0 [g/m2] density of wet aerosol (water + salt)
         )

!    CALL PPP2SALT(     &
!         ZSVT                                   & !I [moments/molec_{air}] moments of sea salt for all modes
!         ,PRHODREF                              & !I [kg/m3] air density
!         ,PSIG3D=ZSIGMA                         & !O [-] dispersion coefficient
!         ,PRG3D=ZRADIUS                         & !O [um] number median radius
!         ,PMASS3D=ZMASS                         & !O [kg/m3] mass of sea salt
!         )

    DO JMDE=1,NMODE_SLT
       !Get sea salt optical properties from look up tables
       CALL SALTOPT_LKT(                     &
            ZRADIUS(:,:,:,JMDE)                   &  !I [um] number median radius for current mode
            ,ZSIGMA(:,:,:,JMDE)                   &  !I [none] dispersion coefficient for current mode
            ,ZMASS(:,:,:, JMDE)                    &  !I [kg/m3] Mass of sea salt for current mode
            ,ZTAU550_MDE(:,:,:,JMDE)         &  !O [-] optical depth at 550 nm wavelength
            ,ZTAU_WVL_MDE(:,:,:,:,JMDE)      &  !O [-] opt.depth(lambda)/opt.depth(550nm)
            ,ZPIZA_WVL_MDE(:,:,:,:,JMDE)     &  !O [-] single scattering coefficient at any wavelength
            ,ZCGA_WVL_MDE(:,:,:,:,JMDE)      &  !O [-] assymetry factor at any wavelength
            ,PZZ(:,:,:)                      &  !I [m] height of layers
            ,KSWB                            &  !I [nbr] number of shortwave bands
            )
    ENDDO  !Loop on modes

    !Erase earlier value of optical depth at 550 nm
    PTAU550(:,:,:)=0.d0
    
    !Get total at 550 nm from all modes 
    DO JMDE=1,NMODE_SLT
       PTAU550(:,:,:) =        &        !sea salt optical depth at 550 nm for all sea salt
            PTAU550(:,:,:)     &        !sea salt optical depth at 550 nm for all sea salt
            + ZTAU550_MDE(:,:,:,JMDE)   !Optical depth for one mode at 550 nm      
    ENDDO
    
    !Initialize output variables
    PTAUREL_WVL(:,:,:,:)=0.d0         !Initialize opt.depth at wvl=lambda
    PCGA_WVL(:,:,:,:)=0.d0           !Initialize assym.factor at wvl=lambda
    PPIZA_WVL(:,:,:,:)=0.d0          !Initialize single scattering albedo at wvl=lambda


    !Find the numerator in the expression for the average of the optical properties   
    DO JMDE=1,NMODE_SLT             !Number of modes
       DO JWVL=1,KSWB                    !Number of SW wavelengths

          !Get sum of optical depth from all modes at wvl
          PTAUREL_WVL(:,:,:,JWVL)  =                &  !new opt.depth(lambda) / opt.depth(550)
               PTAUREL_WVL(:,:,:,JWVL)              &  !old sum for all modes at wvl=lambda
               +ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)       !optical depth for one mode at wvl=lambda
          
          !Get sum of all assymmetry factors  from all modes at wvl=lambda 
          PCGA_WVL(:,:,:,JWVL) =                     &  !New sum of assymetry factors
               PCGA_WVL(:,:,:,JWVL)                  &  !old sum of assymetry factors
               +ZCGA_WVL_MDE(:,:,:,JWVL,JMDE)     &  !Assymetry factor for one mode and one wavelength
               *ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)     &  !Optical depth of this wavelength and mode
               *ZPIZA_WVL_MDE(:,:,:,JWVL,JMDE)       !Fraction of radiation scattered

          !Get sum of single scattering albdedo at wvl=lambda
          PPIZA_WVL(:,:,:,JWVL)  =                  & !New sum of single scattering albedo
               PPIZA_WVL(:,:,:,JWVL)             & !Old sum of single scattering albedo
               +ZPIZA_WVL_MDE(:,:,:,JWVL,JMDE)   & !SSA for onen mode and one wavelength
               *ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)        !Optical depth for this wavelength and mode
          
       ENDDO
    ENDDO
  
    !Compute the output values for sea salt optical properties
    DO JWVL=1,KSWB
       
       !Divide total single scattering albdeo by total optical depth at this wavelength
       !This is needed since we weight all single scattering alebdos by wavelengths just above
       PPIZA_WVL(:,:,:,JWVL) =               &     !The value we want is ....
            PPIZA_WVL(:,:,:,JWVL)            &     !..value weighted by optical depths of all wvl and modes
            /max(epsilon,PTAUREL_WVL(:,:,:,JWVL))            !..divided by the optical depth for all wvl 
       
       
       !Divide total assymetry factor by total optical depth at this wavelength
       !This is needed since we weight all assymetry factors by wavelengths just above
       PCGA_WVL(:,:,:,JWVL) =              &  !The value we want is ....
            PCGA_WVL(:,:,:,JWVL)           &  !..value weighted by optical depths of all wvl and modes
            /                              &
            (max(epsilon,                  &
            (PTAUREL_WVL(:,:,:,JWVL)       &  !..divided scattered fraction of by the optical depth
            *PPIZA_WVL(:,:,:,JWVL))))

       !Finally convert PTAUREL_WVL which was until now an optical depth to a fraction of optical depth
       PTAUREL_WVL(:,:,:,JWVL) =                    &
            PTAUREL_WVL(:,:,:,JWVL)                 &  !Opt.depth at lambda with contr. from all modes
            /max(epsilon,PTAU550(:,:,:))               !Optical depth at 550 contr. from all modes
       
    ENDDO !Loop on wavelenghts
       

    !DEALLOCATE local  arrays which size depend on number of modes
    DEALLOCATE(ZTAU550_MDE)
    DEALLOCATE(ZTAU_WVL_MDE)
    DEALLOCATE(ZPIZA_WVL_MDE)
    DEALLOCATE(ZCGA_WVL_MDE)
    
  END SUBROUTINE SALTOPT_GET

  !*****************************************************************

  SUBROUTINE SALTOPT_LKT(              &
        PRG                            & !I [um] number median radius of aerosol mode
       ,PSIGMA                         & !I [-] lognormal dispersion coefficient
       ,PMASS                          & !I [kg/m3] Mass concentration of sea salt
       ,PTAU550                        & !O [optical depth at 550 nm
       ,PTAU_WVL                       & !O [-] opt.depth(lambda)/opt.depth(550nm)
       ,PPIZA_WVL                      & !O [-] single scattering coefficient at any wavelength
       ,PCGA_WVL                       & !O [-] assymetry factor at any wavelength
       ,PZZ                            & !I [m] height of layers
       ,KSWB                           & !I [nbr] number of short wave bands
       )

    !Purpose: Get optical properties of one sea salt mode from the mass concentration, 
    !dispersion coefficient and number median radius.

    !Use the module with the sea salt optical properties look up tables
    USE MODD_SALT_OPT_LKT  

    IMPLICIT NONE
    !INPUT
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PRG          !I [um] number median radius for one mode
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PSIGMA       !I [-] dispersion coefficient for one mode
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PMASS        !I [kg/m3] mass of sea salt
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ          !I [m] height of layers
    INTEGER, INTENT(IN)                       :: KSWB         !I [nbr] number of shortwave bands

    !OUTPUT
    REAL, DIMENSION(:,:,:), INTENT(OUT)       :: PTAU550      !O [-] optical depth at 550 nm
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PTAU_WVL     !O [-] optical depth at wvl
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PPIZA_WVL    !O [-] single scattering albedo @ wvl
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PCGA_WVL     !O [-] assymetry factor @ wvl

    !LOCALS
    REAL, DIMENSION(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB) :: ZEXT_COEFF_WVL    ![m2/kg] Extinction coefficient at wvl
    REAL, DIMENSION(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3)     ) :: ZEXT_COEFF_550    ![m2/kg] Extinction coefficient at 550nm
    REAL                      :: FACT_SIGMA  ![-] factor needed to get right index in look up table for sigma
    REAL                      :: FACT_RADIUS ![-] factor needed to get right index in look up table for radius
    INTEGER                   :: WVL_IDX     ![idx] counter for wavelengths
    INTEGER                   :: JI, JJ, JK  ![idx] counters for lon, lat and lev
    INTEGER                   :: RG_IDX      ![idx] index for radius to get in look up table
    INTEGER                   :: SG_IDX      ![idx] index for sigma to get in look up table
    REAL,DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZRG     ![um] bounded value for number median radius
    REAL,DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZSIGMA  ![um] bounded value for sigma
    REAL, PARAMETER           :: EPSILON=1.d-8                     ![um] a small number used to avoid zero
    REAL                      :: ZRADIUS_LKT_MAX, ZRADIUS_LKT_MIN  ![um] values limited at midpoint values of bin
    REAL                      :: ZSIGMA_LKT_MAX, ZSIGMA_LKT_MIN    ![-] values limited at midpoint of bin
    INTEGER                   :: JKRAD                             !Index valid for radiation code

    !Limit max and min values to be midpont of bin to avoid 0 or NMAX+1 values
    ZRADIUS_LKT_MAX=exp(log(XRADIUS_LKT_MAX)   &
         - 0.5d0/DBLE(NMAX_RADIUS_LKT)*log(XRADIUS_LKT_MAX/XRADIUS_LKT_MIN))
    ZRADIUS_LKT_MIN=exp(log(XRADIUS_LKT_MIN)   &
         + 0.5d0/DBLE(NMAX_RADIUS_LKT)*log(XRADIUS_LKT_MAX/XRADIUS_LKT_MIN))
    ZSIGMA_LKT_MAX=XSIGMA_LKT_MAX - 0.5d0/DBLE(NMAX_SIGMA_LKT)*(XSIGMA_LKT_MAX-XSIGMA_LKT_MIN)
    ZSIGMA_LKT_MIN=XSIGMA_LKT_MIN + 0.5d0/DBLE(NMAX_SIGMA_LKT)*(XSIGMA_LKT_MAX-XSIGMA_LKT_MIN)

    !Begin code
    FACT_SIGMA = DBLE(NMAX_SIGMA_LKT)/(XSIGMA_LKT_MAX-XSIGMA_LKT_MIN)
    FACT_RADIUS = DBLE(NMAX_RADIUS_LKT)/(LOG(XRADIUS_LKT_MAX/XRADIUS_LKT_MIN))

    !Remove unphysical values for rg
    ZRG(:,:,:) = min( max(ZRADIUS_LKT_MIN,PRG(:,:,:)), ZRADIUS_LKT_MAX)
    ZSIGMA(:,:,:) = min( max(ZSIGMA_LKT_MIN,PSIGMA(:,:,:)), ZSIGMA_LKT_MAX)

    !Initilalize arrays to make sure, they are intent(OUT), 
    !and may be initialized strangely by the computer
    PTAU550(:,:,:)=EPSILON 
    PTAU_WVL(:,:,:,:)=EPSILON     
    PPIZA_WVL(:,:,:,:)=EPSILON
    PCGA_WVL(:,:,:,:)=EPSILON

    DO WVL_IDX = 1,KSWB
       DO JK=2,SIZE(PMASS,3)-1
          JKRAD = JK - 1  !Index in radiation code
          DO JJ=1,SIZE(PMASS,2)
             DO JI=1,SIZE(PMASS,1)

                !Get the correct indexes for the look up tables
                RG_IDX =  nint(                                  &
                     log(ZRG(JI,JJ,JK)/XRADIUS_LKT_MIN)          &
                     *FACT_RADIUS                                &
                     +0.5)
                   
                SG_IDX = nint((ZSIGMA(JI,JJ,JK)-XSIGMA_LKT_MIN)*FACT_SIGMA + 0.5) 

                !Open the look up tables and get the right values out of them
                !The extinction coefficient for 
                ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)         = XEXT_COEFF_WVL_LKT(RG_IDX,SG_IDX,WVL_IDX)
                ZEXT_COEFF_550(JI,JJ,JKRAD)                 = XEXT_COEFF_550_LKT(RG_IDX,SG_IDX)
                !Switch to radiation code indexes for the output values
                PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)              = XPIZA_LKT(RG_IDX,SG_IDX,WVL_IDX)
                PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)               = XCGA_LKT(RG_IDX,SG_IDX,WVL_IDX)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !Get the optical depth of this mode using the looked up extinction coeffient
    DO JK=2,SIZE(PZZ,3)-1
       JKRAD = JK - 1           !Index in radiation code
       PTAU550(:,:,JKRAD) = ZEXT_COEFF_550(:,:,JKRAD) &
            * PMASS(:,:,JK)                    &
            * (PZZ(:,:,JK+1) - PZZ(:,:,JK))
    ENDDO

    !Get the optical depth of whatever wavelength using the looked up tables
    DO WVL_IDX=1,KSWB
       DO JK=2,SIZE(PZZ,3)-1
          JKRAD = JK -1
          PTAU_WVL(:,:,JKRAD,WVL_IDX) =              &
               PMASS(:,:,JK)                      & ![kg/m3] Mass in this mode
               *ZEXT_COEFF_WVL(:,:,JKRAD,WVL_IDX)    & ![m2/kg] mass exinction coefficient
               *(PZZ(:,:,JK+1) - PZZ(:,:,JK))       ![m] Height of layer
       ENDDO !Loop on levels
    ENDDO    !Loop on wavelengths

    !Avoid unphysical values (which might occur on grid edges) on grid edges
    PTAU550(:,:,:)=max(PTAU550(:,:,:),EPSILON)
    PTAU_WVL(:,:,:,:)=max(PTAU_WVL(:,:,:,:),EPSILON)
    PPIZA_WVL(:,:,:,:)=max(PPIZA_WVL(:,:,:,:),EPSILON)
    PCGA_WVL(:,:,:,:)=max(PCGA_WVL(:,:,:,:),EPSILON)

  END SUBROUTINE SALTOPT_LKT

  !***************************************************************************
  SUBROUTINE SALT_OPT_LKT_SET1()

!Purpose: Read the look up tables for sea salt optical properties
!All variables are defined in the module MODD_SALT_OPT_LKT
! Based upon refractive indexes at wavelength intervals as:
!
! New tables from Mallet (LA) and Tulet (CNRM) :
! 0.185-0.25  RI[1]="(1.448,-0.00292)"
! 0.25-0.44   RI[2]="(1.448,-0.00292)"
! 0.44-0.69   RI[3]="(1.448,-0.00292)"
! 0.69-1.19   RI[4]="(1.44023,-0.00116)"
! 1.19-2.38   RI[5]="(1.41163,-0.00106)"
! 2.38-4.0    RI[6]="(1.41163,-0.00106)"
!             RI550="(1.44412,-0.00204)"

    USE MODD_SALT_OPT_LKT

    IMPLICIT NONE
    
    !Here are the output values from the mie program:
XEXT_COEFF_WVL_LKT(1,1,1:6)=(/ 44.560000,11.638000,0.691110,0.260970,0.532810,132.670000 /)
XPIZA_LKT(1,1,1:6)=(/ 0.236252,0.196169,0.900,0.900,0.900,0.000010 /)
XCGA_LKT(1,1,1:6)=(/ 0.008550,0.003990,0.900,0.900,0.900,0.000077 /)
XEXT_COEFF_550_LKT(1,1)=1.860100 !rg=0.0104084 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,2,1:6)=(/ 47.412000,12.247000,0.799480,0.276270,0.534440,132.670000 /)
XPIZA_LKT(1,2,1:6)=(/ 0.279956,0.235019,0.900,0.900,0.900,0.000012 /)
XCGA_LKT(1,2,1:6)=(/ 0.010863,0.005070,0.900,0.900,0.900,0.000097 /)
XEXT_COEFF_550_LKT(1,2)=2.159700 !rg=0.0104084 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,3,1:6)=(/ 54.339000,13.730000,1.064000,0.313590,0.538410,132.680000 /)
XPIZA_LKT(1,3,1:6)=(/ 0.367602,0.315487,0.740010,0.900,0.900,0.000018 /)
XCGA_LKT(1,3,1:6)=(/ 0.016550,0.007733,0.003303,0.900,0.900,0.000150 /)
XEXT_COEFF_550_LKT(1,3)=2.891000 !rg=0.0104084 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,4,1:6)=(/ 68.897000,16.861000,1.624200,0.392640,0.546790,132.700000 /)
XPIZA_LKT(1,4,1:6)=(/ 0.495600,0.439542,0.829287,0.900,0.900,0.000031 /)
XCGA_LKT(1,4,1:6)=(/ 0.028200,0.013210,0.005650,0.900,0.900,0.000257 /)
XEXT_COEFF_550_LKT(1,4)=4.439900 !rg=0.0104084 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,5,1:6)=(/ 97.382000,23.043000,2.734500,0.549220,0.563330,132.740000 /)
XPIZA_LKT(1,5,1:6)=(/ 0.636881,0.586319,0.898213,0.900,0.900,0.000057 /)
XCGA_LKT(1,5,1:6)=(/ 0.047230,0.022183,0.009507,0.900,0.900,0.000430 /)
XEXT_COEFF_550_LKT(1,5)=7.508500 !rg=0.0104084 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,6,1:6)=(/ 148.690000,34.380000,4.780800,0.837730,0.593720,132.790000 /)
XPIZA_LKT(1,6,1:6)=(/ 0.756241,0.719164,0.941444,0.900,0.900,0.000104 /)
XCGA_LKT(1,6,1:6)=(/ 0.073803,0.034727,0.014910,0.900,0.900,0.000677 /)
XEXT_COEFF_550_LKT(1,6)=13.159000 !rg=0.0104084 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,7,1:6)=(/ 237.320000,54.712000,8.480800,1.359600,0.648600,132.860000 /)
XPIZA_LKT(1,7,1:6)=(/ 0.842106,0.820255,0.966711,0.849812,0.900,0.000190 /)
XCGA_LKT(1,7,1:6)=(/ 0.110600,0.052033,0.022380,0.008447,0.900,0.001020 /)
XEXT_COEFF_550_LKT(1,7)=23.353000 !rg=0.0104084 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,8,1:6)=(/ 380.690000,90.131000,15.037000,2.285800,0.745810,132.980000 /)
XPIZA_LKT(1,8,1:6)=(/ 0.897341,0.888105,0.980997,0.910242,0.900,0.000342 /)
XCGA_LKT(1,8,1:6)=(/ 0.162557,0.076233,0.032820,0.012407,0.900,0.001500 /)
XEXT_COEFF_550_LKT(1,8)=41.317000 !rg=0.0104084 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,9,1:6)=(/ 598.500000,150.890000,26.698000,3.940900,0.919320,133.160000 /)
XPIZA_LKT(1,9,1:6)=(/ 0.931199,0.930874,0.989112,0.947576,0.900,0.000612 /)
XCGA_LKT(1,9,1:6)=(/ 0.237567,0.110900,0.047703,0.018063,0.900,0.002187 /)
XEXT_COEFF_550_LKT(1,9)=72.858000 !rg=0.0104084 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,10,1:6)=(/ 911.710000,249.590000,47.113000,6.872900,1.226600,133.430000 /)
XPIZA_LKT(1,10,1:6)=(/ 0.951701,0.956398,0.993682,0.969639,0.566268,0.001091 /)
XCGA_LKT(1,10,1:6)=(/ 0.339350,0.161273,0.069060,0.026187,0.008507,0.003180 /)
XEXT_COEFF_550_LKT(1,10)=126.490000 !rg=0.0104084 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,11,1:6)=(/ 1366.300000,400.050000,82.608000,12.115000,1.776800,133.850000 /)
XPIZA_LKT(1,11,1:6)=(/ 0.965026,0.971332,0.996279,0.982528,0.699111,0.001947 /)
XCGA_LKT(1,11,1:6)=(/ 0.441663,0.235600,0.100097,0.037943,0.012350,0.004617 /)
XEXT_COEFF_550_LKT(1,11)=214.270000 !rg=0.0104084 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,12,1:6)=(/ 1929.500000,616.170000,141.770000,21.419000,2.759200,134.500000 /)
XPIZA_LKT(1,12,1:6)=(/ 0.973261,0.980077,0.997741,0.989916,0.804842,0.003467 /)
XCGA_LKT(1,12,1:6)=(/ 0.517987,0.339687,0.145627,0.054980,0.017927,0.006713 /)
XEXT_COEFF_550_LKT(1,12)=346.630000 !rg=0.0104084 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,13,1:6)=(/ 2585.400000,931.540000,233.530000,37.746000,4.511400,135.520000 /)
XPIZA_LKT(1,13,1:6)=(/ 0.978385,0.985649,0.998557,0.994117,0.879398,0.006157 /)
XCGA_LKT(1,13,1:6)=(/ 0.588277,0.447203,0.213210,0.079743,0.026027,0.009760 /)
XEXT_COEFF_550_LKT(1,13)=535.890000 !rg=0.0104084 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,14,1:6)=(/ 3246.400000,1322.500000,364.340000,65.547000,7.614900,137.130000 /)
XPIZA_LKT(1,14,1:6)=(/ 0.981540,0.989084,0.999015,0.996487,0.927491,0.010868 /)
XCGA_LKT(1,14,1:6)=(/ 0.640103,0.522337,0.311457,0.115893,0.037767,0.014190 /)
XEXT_COEFF_550_LKT(1,14)=810.600000 !rg=0.0104084 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,15,1:6)=(/ 3829.700000,1770.300000,553.770000,110.530000,13.093000,139.680000 /)
XPIZA_LKT(1,15,1:6)=(/ 0.983356,0.991165,0.999295,0.997820,0.956950,0.019046 /)
XCGA_LKT(1,15,1:6)=(/ 0.679310,0.591493,0.427100,0.169283,0.054790,0.020630 /)
XEXT_COEFF_550_LKT(1,15)=1145.900000 !rg=0.0104084 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,16,1:6)=(/ 4221.900000,2225.400000,807.350000,177.600000,22.687000,143.760000 /)
XPIZA_LKT(1,16,1:6)=(/ 0.984073,0.992462,0.999475,0.998567,0.974442,0.033063 /)
XCGA_LKT(1,16,1:6)=(/ 0.704653,0.642720,0.510100,0.248870,0.079537,0.030003 /)
XEXT_COEFF_550_LKT(1,16)=1530.500000 !rg=0.0104084 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,17,1:6)=(/ 4318.100000,2623.700000,1089.700000,271.120000,38.929000,150.290000 /)
XPIZA_LKT(1,17,1:6)=(/ 0.983733,0.993197,0.999580,0.998991,0.984546,0.056260 /)
XCGA_LKT(1,17,1:6)=(/ 0.716767,0.681200,0.578217,0.360207,0.115597,0.043607 /)
XEXT_COEFF_550_LKT(1,17)=1910.100000 !rg=0.0104084 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,18,1:6)=(/ 4110.000000,2886.600000,1389.700000,406.750000,65.036000,160.730000 /)
XPIZA_LKT(1,18,1:6)=(/ 0.982309,0.993474,0.999646,0.999265,0.990318,0.093070 /)
XCGA_LKT(1,18,1:6)=(/ 0.717300,0.705453,0.633967,0.469667,0.168723,0.063427 /)
XEXT_COEFF_550_LKT(1,18)=2237.600000 !rg=0.0104084 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,19,1:6)=(/ 3723.600000,2951.500000,1665.100000,570.310000,103.620000,177.330000 /)
XPIZA_LKT(1,19,1:6)=(/ 0.979843,0.993314,0.999685,0.999435,0.993581,0.147741 /)
XCGA_LKT(1,19,1:6)=(/ 0.711833,0.716653,0.675073,0.540017,0.247783,0.092527 /)
XEXT_COEFF_550_LKT(1,19)=2447.200000 !rg=0.0104084 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,20,1:6)=(/ 3309.600000,2805.500000,1861.900000,752.610000,156.950000,202.790000 /)
XPIZA_LKT(1,20,1:6)=(/ 0.977031,0.992694,0.999702,0.999537,0.995451,0.220490 /)
XCGA_LKT(1,20,1:6)=(/ 0.711727,0.715933,0.702040,0.606347,0.359080,0.135837 /)
XEXT_COEFF_550_LKT(1,20)=2483.500000 !rg=0.0104084 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,1,1:6)=(/ 47.546000,12.274000,0.804030,0.276910,0.534520,132.670000 /)
XPIZA_LKT(2,1,1:6)=(/ 0.281570,0.236515,0.900,0.900,0.900,0.000012 /)
XCGA_LKT(2,1,1:6)=(/ 0.010030,0.004680,0.900,0.900,0.900,0.000090 /)
XEXT_COEFF_550_LKT(2,1)=2.172300 !rg=0.011276 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,2,1:6)=(/ 51.162000,13.048000,0.941810,0.296360,0.536590,132.680000 /)
XPIZA_LKT(2,2,1:6)=(/ 0.329919,0.280516,0.900,0.900,0.900,0.000015 /)
XCGA_LKT(2,2,1:6)=(/ 0.012740,0.005950,0.900,0.900,0.900,0.000113 /)
XEXT_COEFF_550_LKT(2,2)=2.553200 !rg=0.011276 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,3,1:6)=(/ 59.945000,14.931000,1.278100,0.343820,0.541620,132.690000 /)
XPIZA_LKT(2,3,1:6)=(/ 0.423716,0.368928,0.783330,0.900,0.900,0.000023 /)
XCGA_LKT(2,3,1:6)=(/ 0.019403,0.009073,0.003877,0.900,0.900,0.000177 /)
XEXT_COEFF_550_LKT(2,3)=3.483000 !rg=0.011276 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,4,1:6)=(/ 78.388000,18.907000,1.990500,0.444310,0.552260,132.720000 /)
XPIZA_LKT(2,4,1:6)=(/ 0.553549,0.498453,0.860491,0.900,0.900,0.000040 /)
XCGA_LKT(2,4,1:6)=(/ 0.033043,0.015490,0.006630,0.900,0.900,0.000300 /)
XEXT_COEFF_550_LKT(2,4)=5.452500 !rg=0.011276 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,5,1:6)=(/ 114.350000,26.754000,3.402400,0.643400,0.573270,132.750000 /)
XPIZA_LKT(2,5,1:6)=(/ 0.687801,0.641972,0.918016,0.900,0.900,0.000072 /)
XCGA_LKT(2,5,1:6)=(/ 0.055300,0.026000,0.011150,0.900,0.900,0.000507 /)
XEXT_COEFF_550_LKT(2,5)=9.353600 !rg=0.011276 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,6,1:6)=(/ 178.600000,41.128000,6.004000,1.010200,0.611880,132.810000 /)
XPIZA_LKT(2,6,1:6)=(/ 0.794457,0.763624,0.953228,0.798291,0.900,0.000133 /)
XCGA_LKT(2,6,1:6)=(/ 0.086390,0.040687,0.017483,0.006593,0.900,0.000797 /)
XEXT_COEFF_550_LKT(2,6)=16.533000 !rg=0.011276 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,7,1:6)=(/ 287.750000,66.819000,10.706000,1.673800,0.681600,132.910000 /)
XPIZA_LKT(2,7,1:6)=(/ 0.867618,0.851415,0.973513,0.877792,0.900,0.000241 /)
XCGA_LKT(2,7,1:6)=(/ 0.129530,0.060947,0.026230,0.009907,0.900,0.001197 /)
XEXT_COEFF_550_LKT(2,7)=29.466000 !rg=0.011276 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,8,1:6)=(/ 459.410000,111.230000,19.030000,2.851600,0.805140,133.040000 /)
XPIZA_LKT(2,8,1:6)=(/ 0.913183,0.908162,0.984888,0.927866,0.900,0.000434 /)
XCGA_LKT(2,8,1:6)=(/ 0.190560,0.089307,0.038457,0.014550,0.900,0.001760 /)
XEXT_COEFF_550_LKT(2,8)=52.174000 !rg=0.011276 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,9,1:6)=(/ 713.510000,186.180000,33.797000,4.955600,1.025700,133.260000 /)
XPIZA_LKT(2,9,1:6)=(/ 0.940778,0.943045,0.991321,0.958156,0.482397,0.000778 /)
XCGA_LKT(2,9,1:6)=(/ 0.277407,0.130047,0.055887,0.021180,0.006873,0.002567 /)
XEXT_COEFF_550_LKT(2,9)=91.720000 !rg=0.011276 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,10,1:6)=(/ 1080.500000,304.420000,59.495000,8.681300,1.416200,133.580000 /)
XPIZA_LKT(2,10,1:6)=(/ 0.957868,0.963519,0.994936,0.975836,0.623672,0.001386 /)
XCGA_LKT(2,10,1:6)=(/ 0.384113,0.189440,0.080920,0.030693,0.009980,0.003730 /)
XEXT_COEFF_550_LKT(2,10)=157.800000 !rg=0.011276 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,11,1:6)=(/ 1587.100000,479.220000,103.580000,15.333000,2.115700,134.080000 /)
XPIZA_LKT(2,11,1:6)=(/ 0.968922,0.975446,0.996985,0.986091,0.746638,0.002472 /)
XCGA_LKT(2,11,1:6)=(/ 0.474010,0.276253,0.117400,0.044460,0.014483,0.005420 /)
XEXT_COEFF_550_LKT(2,11)=262.690000 !rg=0.011276 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,12,1:6)=(/ 2186.000000,732.200000,175.150000,27.097000,3.364300,134.870000 /)
XPIZA_LKT(2,12,1:6)=(/ 0.975601,0.982649,0.998135,0.991945,0.839333,0.004399 /)
XCGA_LKT(2,12,1:6)=(/ 0.547573,0.387420,0.171150,0.064423,0.021020,0.007873 /)
XEXT_COEFF_550_LKT(2,12)=416.010000 !rg=0.011276 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,13,1:6)=(/ 2858.000000,1087.300000,282.170000,47.573000,5.590400,136.100000 /)
XPIZA_LKT(2,13,1:6)=(/ 0.979881,0.987295,0.998777,0.995267,0.902141,0.007803 /)
XCGA_LKT(2,13,1:6)=(/ 0.610830,0.480540,0.250883,0.093483,0.030507,0.011450 /)
XEXT_COEFF_550_LKT(2,13)=637.860000 !rg=0.011276 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,14,1:6)=(/ 3505.100000,1497.900000,433.730000,81.819000,9.528000,138.050000 /)
XPIZA_LKT(2,14,1:6)=(/ 0.982447,0.990043,0.999145,0.997135,0.941599,0.013743 /)
XCGA_LKT(2,14,1:6)=(/ 0.657780,0.551220,0.360907,0.136073,0.044257,0.016643 /)
XEXT_COEFF_550_LKT(2,14)=945.150000 !rg=0.011276 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,15,1:6)=(/ 4021.800000,1959.600000,653.790000,135.440000,16.455000,141.140000 /)
XPIZA_LKT(2,15,1:6)=(/ 0.983782,0.991786,0.999381,0.998182,0.965377,0.024000 /)
XCGA_LKT(2,15,1:6)=(/ 0.691000,0.614440,0.466437,0.199300,0.064207,0.024197 /)
XEXT_COEFF_550_LKT(2,15)=1296.800000 !rg=0.011276 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,16,1:6)=(/ 4299.100000,2404.400000,918.780000,212.640000,28.489000,146.120000 /)
XPIZA_LKT(2,16,1:6)=(/ 0.984069,0.992830,0.999525,0.998770,0.979355,0.041417 /)
XCGA_LKT(2,16,1:6)=(/ 0.710930,0.660337,0.537647,0.292707,0.093253,0.035190 /)
XEXT_COEFF_550_LKT(2,16)=1689.700000 !rg=0.011276 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,17,1:6)=(/ 4264.200000,2755.600000,1217.900000,322.070000,48.499000,154.080000 /)
XPIZA_LKT(2,17,1:6)=(/ 0.983272,0.993363,0.999612,0.999119,0.987370,0.069804 /)
XCGA_LKT(2,17,1:6)=(/ 0.717813,0.692710,0.604603,0.410823,0.135760,0.051170 /)
XEXT_COEFF_550_LKT(2,17)=2058.600000 !rg=0.011276 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,18,1:6)=(/ 3964.400000,2939.000000,1509.800000,473.950000,79.610000,166.800000 /)
XPIZA_LKT(2,18,1:6)=(/ 0.981385,0.993462,0.999665,0.999348,0.991918,0.113745 /)
XCGA_LKT(2,18,1:6)=(/ 0.714930,0.711577,0.652200,0.501890,0.198747,0.074517 /)
XEXT_COEFF_550_LKT(2,18)=2343.700000 !rg=0.011276 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,19,1:6)=(/ 3547.400000,2911.400000,1758.300000,643.240000,123.930000,186.830000 /)
XPIZA_LKT(2,19,1:6)=(/ 0.978757,0.993113,0.999694,0.999482,0.994486,0.176498 /)
XCGA_LKT(2,19,1:6)=(/ 0.711247,0.717380,0.687730,0.568360,0.291830,0.109003 /)
XEXT_COEFF_550_LKT(2,19)=2483.200000 !rg=0.011276 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,20,1:6)=(/ 3157.600000,2697.100000,1909.300000,829.100000,186.150000,216.860000 /)
XPIZA_LKT(2,20,1:6)=(/ 0.975754,0.992315,0.999702,0.999568,0.996019,0.255288 /)
XCGA_LKT(2,20,1:6)=(/ 0.714530,0.713313,0.709220,0.628327,0.411073,0.160780 /)
XEXT_COEFF_550_LKT(2,20)=2441.500000 !rg=0.011276 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,1,1:6)=(/ 51.331000,13.082000,0.947590,0.297180,0.536680,132.680000 /)
XPIZA_LKT(3,1,1:6)=(/ 0.331667,0.282184,0.900,0.900,0.900,0.000016 /)
XCGA_LKT(3,1,1:6)=(/ 0.011767,0.005493,0.900,0.900,0.900,0.000107 /)
XEXT_COEFF_550_LKT(3,1)=2.569200 !rg=0.0122159 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,2,1:6)=(/ 55.917000,14.064000,1.122800,0.321910,0.539310,132.690000 /)
XPIZA_LKT(3,2,1:6)=(/ 0.383882,0.330940,0.753502,0.900,0.900,0.000020 /)
XCGA_LKT(3,2,1:6)=(/ 0.014943,0.006980,0.002980,0.900,0.900,0.000133 /)
XEXT_COEFF_550_LKT(3,2)=3.053600 !rg=0.0122159 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,3,1:6)=(/ 67.053000,16.455000,1.550400,0.382240,0.545700,132.700000 /)
XPIZA_LKT(3,3,1:6)=(/ 0.481639,0.425689,0.821151,0.900,0.900,0.000030 /)
XCGA_LKT(3,3,1:6)=(/ 0.022750,0.010640,0.004550,0.900,0.900,0.000207 /)
XEXT_COEFF_550_LKT(3,3)=4.235700 !rg=0.0122159 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,4,1:6)=(/ 90.401000,21.504000,2.456300,0.510000,0.559210,132.730000 /)
XPIZA_LKT(3,4,1:6)=(/ 0.609727,0.557253,0.886744,0.900,0.900,0.000050 /)
XCGA_LKT(3,4,1:6)=(/ 0.038707,0.018163,0.007777,0.900,0.900,0.000353 /)
XEXT_COEFF_550_LKT(3,4)=6.739800 !rg=0.0122159 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,5,1:6)=(/ 135.700000,31.463000,4.251500,0.763130,0.585890,132.780000 /)
XPIZA_LKT(3,5,1:6)=(/ 0.734036,0.693831,0.934225,0.900,0.900,0.000092 /)
XCGA_LKT(3,5,1:6)=(/ 0.064727,0.030470,0.013077,0.900,0.900,0.000593 /)
XEXT_COEFF_550_LKT(3,5)=11.699000 !rg=0.0122159 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,6,1:6)=(/ 215.720000,49.667000,7.558900,1.229600,0.634950,132.850000 /)
XPIZA_LKT(3,6,1:6)=(/ 0.827365,0.802714,0.962715,0.834041,0.900,0.000169 /)
XCGA_LKT(3,6,1:6)=(/ 0.101103,0.047657,0.020497,0.007733,0.900,0.000933 /)
XEXT_COEFF_550_LKT(3,6)=20.815000 !rg=0.0122159 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,7,1:6)=(/ 348.740000,82.048000,13.534000,2.073400,0.723540,132.960000 /)
XPIZA_LKT(3,7,1:6)=(/ 0.888778,0.877675,0.978936,0.901140,0.900,0.000307 /)
XCGA_LKT(3,7,1:6)=(/ 0.151690,0.071380,0.030743,0.011620,0.900,0.001403 /)
XEXT_COEFF_550_LKT(3,7)=37.204000 !rg=0.0122159 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,8,1:6)=(/ 551.750000,137.420000,24.093000,3.570800,0.880550,133.120000 /)
XPIZA_LKT(3,8,1:6)=(/ 0.926077,0.924591,0.987976,0.942224,0.900,0.000552 /)
XCGA_LKT(3,8,1:6)=(/ 0.223140,0.104640,0.045060,0.017063,0.900,0.002067 /)
XEXT_COEFF_550_LKT(3,8)=65.830000 !rg=0.0122159 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,9,1:6)=(/ 848.230000,228.870000,42.755000,6.245400,1.160800,133.370000 /)
XPIZA_LKT(3,9,1:6)=(/ 0.948710,0.952825,0.993069,0.966654,0.542007,0.000989 /)
XCGA_LKT(3,9,1:6)=(/ 0.320877,0.152570,0.065473,0.024827,0.008063,0.003013 /)
XEXT_COEFF_550_LKT(3,9)=115.110000 !rg=0.0122159 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,10,1:6)=(/ 1274.500000,368.480000,74.954000,10.977000,1.657300,133.760000 /)
XPIZA_LKT(3,10,1:6)=(/ 0.963054,0.969185,0.995925,0.980772,0.677734,0.001761 /)
XCGA_LKT(3,10,1:6)=(/ 0.423877,0.222503,0.094840,0.035970,0.011703,0.004377 /)
XEXT_COEFF_550_LKT(3,10)=195.590000 !rg=0.0122159 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,11,1:6)=(/ 1819.000000,571.210000,129.190000,19.408000,2.546400,134.370000 /)
XPIZA_LKT(3,11,1:6)=(/ 0.972036,0.978791,0.997541,0.988916,0.788847,0.003138 /)
XCGA_LKT(3,11,1:6)=(/ 0.503520,0.321590,0.137783,0.052097,0.016983,0.006357 /)
XEXT_COEFF_550_LKT(3,11)=319.060000 !rg=0.0122159 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,12,1:6)=(/ 2460.100000,867.510000,214.460000,34.241000,4.133200,135.310000 /)
XPIZA_LKT(3,12,1:6)=(/ 0.977590,0.984825,0.998444,0.993550,0.868641,0.005579 /)
XCGA_LKT(3,12,1:6)=(/ 0.575833,0.430313,0.201303,0.075493,0.024640,0.009240 /)
XEXT_COEFF_550_LKT(3,12)=496.630000 !rg=0.0122159 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,13,1:6)=(/ 3128.200000,1249.100000,338.090000,59.765000,6.960200,136.800000 /)
XPIZA_LKT(3,13,1:6)=(/ 0.981071,0.988599,0.998951,0.996174,0.920901,0.009881 /)
XCGA_LKT(3,13,1:6)=(/ 0.630783,0.509427,0.294347,0.109640,0.035753,0.013430 /)
XEXT_COEFF_550_LKT(3,13)=756.470000 !rg=0.0122159 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,14,1:6)=(/ 3735.600000,1687.300000,515.930000,101.440000,11.950000,139.160000 /)
XPIZA_LKT(3,14,1:6)=(/ 0.983113,0.990855,0.999255,0.997645,0.953023,0.017357 /)
XCGA_LKT(3,14,1:6)=(/ 0.672643,0.579873,0.409377,0.159930,0.051860,0.019523 /)
XEXT_COEFF_550_LKT(3,14)=1084.000000 !rg=0.0122159 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,15,1:6)=(/ 4167.500000,2143.500000,760.230000,164.280000,20.686000,142.930000 /)
XPIZA_LKT(3,15,1:6)=(/ 0.984011,0.992272,0.999450,0.998466,0.972123,0.030178 /)
XCGA_LKT(3,15,1:6)=(/ 0.700290,0.633983,0.497710,0.234790,0.075247,0.028380 /)
XEXT_COEFF_550_LKT(3,15)=1458.900000 !rg=0.0122159 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,16,1:6)=(/ 4324.400000,2561.600000,1037.900000,253.090000,35.689000,149.000000 /)
XPIZA_LKT(3,16,1:6)=(/ 0.983854,0.993104,0.999565,0.998934,0.983257,0.051702 /)
XCGA_LKT(3,16,1:6)=(/ 0.714887,0.675123,0.566530,0.341793,0.109393,0.041280 /)
XEXT_COEFF_550_LKT(3,16)=1844.900000 !rg=0.0122159 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,17,1:6)=(/ 4167.300000,2853.200000,1340.600000,381.680000,60.043000,158.710000 /)
XPIZA_LKT(3,17,1:6)=(/ 0.982617,0.993456,0.999637,0.999228,0.989597,0.086129 /)
XCGA_LKT(3,17,1:6)=(/ 0.717400,0.701700,0.625880,0.455097,0.159620,0.060070 /)
XEXT_COEFF_550_LKT(3,17)=2187.700000 !rg=0.0122159 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,18,1:6)=(/ 3795.400000,2952.800000,1622.500000,542.310000,96.471000,174.160000 /)
XPIZA_LKT(3,18,1:6)=(/ 0.980426,0.993366,0.999680,0.999413,0.993172,0.137827 /)
XCGA_LKT(3,18,1:6)=(/ 0.713457,0.715037,0.668790,0.529333,0.234337,0.087627 /)
XEXT_COEFF_550_LKT(3,18)=2420.000000 !rg=0.0122159 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,19,1:6)=(/ 3388.800000,2839.400000,1836.100000,721.810000,147.260000,198.130000 /)
XPIZA_LKT(3,19,1:6)=(/ 0.977562,0.992832,0.999700,0.999523,0.995213,0.208269 /)
XCGA_LKT(3,19,1:6)=(/ 0.712053,0.716420,0.698140,0.596630,0.341537,0.128640 /)
XEXT_COEFF_550_LKT(3,19)=2487.500000 !rg=0.0122159 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,20,1:6)=(/ 3004.000000,2582.400000,1931.300000,902.270000,220.520000,233.080000 /)
XPIZA_LKT(3,20,1:6)=(/ 0.974244,0.991844,0.999701,0.999591,0.996506,0.291052 /)
XCGA_LKT(3,20,1:6)=(/ 0.716293,0.710487,0.714023,0.646823,0.457140,0.190710 /)
XEXT_COEFF_550_LKT(3,20)=2374.700000 !rg=0.0122159 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,1,1:6)=(/ 56.129000,14.108000,1.130100,0.322950,0.539420,132.690000 /)
XPIZA_LKT(4,1,1:6)=(/ 0.385729,0.332758,0.755065,0.900,0.900,0.000020 /)
XCGA_LKT(4,1,1:6)=(/ 0.013800,0.006443,0.002753,0.900,0.900,0.000123 /)
XEXT_COEFF_550_LKT(4,1)=3.073900 !rg=0.0132342 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,2,1:6)=(/ 61.944000,15.354000,1.352900,0.354380,0.542760,132.700000 /)
XPIZA_LKT(4,2,1:6)=(/ 0.440661,0.385476,0.795181,0.900,0.900,0.000025 /)
XCGA_LKT(4,2,1:6)=(/ 0.017523,0.008187,0.003497,0.900,0.900,0.000157 /)
XEXT_COEFF_550_LKT(4,2)=3.689700 !rg=0.0132342 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,3,1:6)=(/ 76.063000,18.391000,1.896600,0.431080,0.550880,132.720000 /)
XPIZA_LKT(4,3,1:6)=(/ 0.539793,0.484357,0.853578,0.900,0.900,0.000038 /)
XCGA_LKT(4,3,1:6)=(/ 0.026667,0.012480,0.005337,0.900,0.900,0.000243 /)
XEXT_COEFF_550_LKT(4,3)=5.192900 !rg=0.0132342 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,4,1:6)=(/ 105.590000,24.802000,3.048600,0.593530,0.568030,132.750000 /)
XPIZA_LKT(4,4,1:6)=(/ 0.662747,0.614327,0.908559,0.900,0.900,0.000064 /)
XCGA_LKT(4,4,1:6)=(/ 0.045330,0.021290,0.009123,0.900,0.900,0.000413 /)
XEXT_COEFF_550_LKT(4,4)=8.376700 !rg=0.0132342 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,5,1:6)=(/ 162.470000,37.433000,5.331200,0.915380,0.601930,132.800000 /)
XPIZA_LKT(4,5,1:6)=(/ 0.775085,0.740980,0.947391,0.900,0.900,0.000117 /)
XCGA_LKT(4,5,1:6)=(/ 0.075733,0.035697,0.015333,0.900,0.900,0.000697 /)
XEXT_COEFF_550_LKT(4,5)=14.678000 !rg=0.0132342 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,6,1:6)=(/ 261.390000,60.454000,9.535200,1.508500,0.664260,132.890000 /)
XPIZA_LKT(4,6,1:6)=(/ 0.855240,0.836446,0.970318,0.864506,0.900,0.000214 /)
XCGA_LKT(4,6,1:6)=(/ 0.118290,0.055813,0.024027,0.009073,0.900,0.001097 /)
XEXT_COEFF_550_LKT(4,6)=26.249000 !rg=0.0132342 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,7,1:6)=(/ 421.600000,101.120000,17.122000,2.581500,0.776840,133.020000 /)
XPIZA_LKT(4,7,1:6)=(/ 0.906158,0.899521,0.983251,0.920406,0.900,0.000390 /)
XCGA_LKT(4,7,1:6)=(/ 0.177573,0.083593,0.036027,0.013627,0.900,0.001650 /)
XEXT_COEFF_550_LKT(4,7)=46.979000 !rg=0.0132342 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,8,1:6)=(/ 659.610000,169.620000,30.500000,4.485300,0.976390,133.210000 /)
XPIZA_LKT(4,8,1:6)=(/ 0.936605,0.937931,0.990421,0.953843,0.900,0.000701 /)
XCGA_LKT(4,8,1:6)=(/ 0.260397,0.122630,0.052787,0.020003,0.900,0.002423 /)
XEXT_COEFF_550_LKT(4,8)=82.923000 !rg=0.0132342 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,9,1:6)=(/ 1006.900000,279.770000,54.017000,7.884000,1.332700,133.510000 /)
XPIZA_LKT(4,9,1:6)=(/ 0.955383,0.960640,0.994451,0.973453,0.600381,0.001256 /)
XCGA_LKT(4,9,1:6)=(/ 0.364720,0.179060,0.076707,0.029103,0.009460,0.003537 /)
XEXT_COEFF_550_LKT(4,9)=143.840000 !rg=0.0132342 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,10,1:6)=(/ 1486.800000,442.710000,94.105000,13.889000,1.963800,133.980000 /)
XPIZA_LKT(4,10,1:6)=(/ 0.967296,0.973712,0.996705,0.984696,0.727355,0.002237 /)
XCGA_LKT(4,10,1:6)=(/ 0.457603,0.260857,0.111193,0.042150,0.013727,0.005137 /)
XEXT_COEFF_550_LKT(4,10)=240.520000 !rg=0.0132342 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,11,1:6)=(/ 2067.800000,679.630000,160.040000,24.557000,3.093900,134.710000 /)
XPIZA_LKT(4,11,1:6)=(/ 0.974578,0.981578,0.997977,0.991153,0.825590,0.003983 /)
XCGA_LKT(4,11,1:6)=(/ 0.533303,0.368647,0.161830,0.061040,0.019913,0.007460 /)
XEXT_COEFF_550_LKT(4,11)=384.150000 !rg=0.0132342 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,12,1:6)=(/ 2734.200000,1017.400000,260.070000,43.189000,5.110100,135.850000 /)
XPIZA_LKT(4,12,1:6)=(/ 0.979232,0.986620,0.998687,0.994819,0.893204,0.007072 /)
XCGA_LKT(4,12,1:6)=(/ 0.599820,0.465633,0.236790,0.088487,0.028883,0.010837 /)
XEXT_COEFF_550_LKT(4,12)=591.820000 !rg=0.0132342 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,13,1:6)=(/ 3395.000000,1419.500000,403.090000,74.739000,8.697800,137.650000 /)
XPIZA_LKT(4,13,1:6)=(/ 0.982063,0.989641,0.999093,0.996888,0.936239,0.012500 /)
XCGA_LKT(4,13,1:6)=(/ 0.649387,0.538207,0.342343,0.128680,0.041900,0.015753 /)
XEXT_COEFF_550_LKT(4,13)=886.720000 !rg=0.0132342 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,14,1:6)=(/ 3944.100000,1878.100000,611.120000,124.700000,15.012000,140.520000 /)
XPIZA_LKT(4,14,1:6)=(/ 0.983598,0.991533,0.999347,0.998044,0.962222,0.021887 /)
XCGA_LKT(4,14,1:6)=(/ 0.685303,0.604440,0.451193,0.188177,0.060767,0.022897 /)
XEXT_COEFF_550_LKT(4,14)=1230.200000 !rg=0.0132342 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,15,1:6)=(/ 4270.200000,2326.800000,869.450000,197.330000,25.987000,145.100000 /)
XPIZA_LKT(4,15,1:6)=(/ 0.984070,0.992672,0.999505,0.998690,0.977505,0.037848 /)
XCGA_LKT(4,15,1:6)=(/ 0.707690,0.652347,0.525540,0.276307,0.088207,0.033287 /)
XEXT_COEFF_550_LKT(4,15)=1620.000000 !rg=0.0132342 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,16,1:6)=(/ 4295.000000,2703.100000,1165.100000,300.780000,44.536000,152.510000 /)
XPIZA_LKT(4,16,1:6)=(/ 0.983480,0.993292,0.999599,0.999070,0.986347,0.064267 /)
XCGA_LKT(4,16,1:6)=(/ 0.717037,0.687500,0.594247,0.392637,0.128423,0.048437 /)
XEXT_COEFF_550_LKT(4,16)=1998.300000 !rg=0.0132342 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,17,1:6)=(/ 4032.200000,2921.900000,1460.800000,447.240000,73.730000,164.350000 /)
XPIZA_LKT(4,17,1:6)=(/ 0.981830,0.993465,0.999657,0.999318,0.991348,0.105559 /)
XCGA_LKT(4,17,1:6)=(/ 0.716193,0.708803,0.644547,0.489940,0.187917,0.070553 /)
XEXT_COEFF_550_LKT(4,17)=2303.000000 !rg=0.0132342 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,18,1:6)=(/ 3633.200000,2929.200000,1719.800000,613.230000,115.740000,183.020000 /)
XPIZA_LKT(4,18,1:6)=(/ 0.979244,0.993197,0.999690,0.999464,0.994159,0.165321 /)
XCGA_LKT(4,18,1:6)=(/ 0.711910,0.716743,0.682363,0.557200,0.276137,0.103173 /)
XEXT_COEFF_550_LKT(4,18)=2470.400000 !rg=0.0132342 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,19,1:6)=(/ 3230.900000,2744.800000,1890.800000,799.620000,174.700000,211.420000 /)
XPIZA_LKT(4,19,1:6)=(/ 0.976332,0.992462,0.999702,0.999556,0.995817,0.242377 /)
XCGA_LKT(4,19,1:6)=(/ 0.713503,0.714147,0.706050,0.620307,0.393647,0.152137 /)
XEXT_COEFF_550_LKT(4,19)=2459.200000 !rg=0.0132342 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,20,1:6)=(/ 2836.700000,2462.100000,1929.400000,974.610000,258.680000,251.510000 /)
XPIZA_LKT(4,20,1:6)=(/ 0.972756,0.991328,0.999694,0.999611,0.996913,0.326614 /)
XCGA_LKT(4,20,1:6)=(/ 0.719310,0.707503,0.716467,0.664303,0.493227,0.226470 /)
XEXT_COEFF_550_LKT(4,20)=2286.900000 !rg=0.0132342 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,1,1:6)=(/ 62.213000,15.409000,1.362200,0.355710,0.542910,132.700000 /)
XPIZA_LKT(5,1,1:6)=(/ 0.442562,0.387409,0.796549,0.900,0.900,0.000025 /)
XCGA_LKT(5,1,1:6)=(/ 0.016187,0.007560,0.003230,0.900,0.900,0.000147 /)
XEXT_COEFF_550_LKT(5,1)=3.715600 !rg=0.0143374 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,2,1:6)=(/ 69.588000,16.993000,1.645400,0.395670,0.547140,132.710000 /)
XPIZA_LKT(5,2,1:6)=(/ 0.498802,0.442941,0.831364,0.900,0.900,0.000032 /)
XCGA_LKT(5,2,1:6)=(/ 0.020547,0.009603,0.004103,0.900,0.900,0.000187 /)
XEXT_COEFF_550_LKT(5,2)=4.498600 !rg=0.0143374 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,3,1:6)=(/ 87.482000,20.849000,2.336800,0.493180,0.557460,132.730000 /)
XPIZA_LKT(5,3,1:6)=(/ 0.596573,0.543314,0.880952,0.900,0.900,0.000048 /)
XCGA_LKT(5,3,1:6)=(/ 0.031253,0.014637,0.006260,0.900,0.900,0.000283 /)
XEXT_COEFF_550_LKT(5,3)=6.409900 !rg=0.0143374 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,4,1:6)=(/ 124.770000,28.988000,3.801700,0.699730,0.579240,132.770000 /)
XPIZA_LKT(5,4,1:6)=(/ 0.711521,0.668226,0.926497,0.900,0.900,0.000082 /)
XCGA_LKT(5,4,1:6)=(/ 0.053070,0.024957,0.010700,0.900,0.900,0.000487 /)
XEXT_COEFF_550_LKT(5,4)=10.458000 !rg=0.0143374 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,5,1:6)=(/ 195.850000,44.995000,6.703900,1.109000,0.622300,132.830000 /)
XPIZA_LKT(5,5,1:6)=(/ 0.810804,0.782892,0.958020,0.816098,0.900,0.000149 /)
XCGA_LKT(5,5,1:6)=(/ 0.088580,0.041813,0.017980,0.006780,0.900,0.000820 /)
XEXT_COEFF_550_LKT(5,5)=18.462000 !rg=0.0143374 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,6,1:6)=(/ 317.020000,74.042000,12.046000,1.863200,0.701510,132.930000 /)
XPIZA_LKT(5,6,1:6)=(/ 0.878524,0.865078,0.976389,0.890087,0.900,0.000272 /)
XCGA_LKT(5,6,1:6)=(/ 0.138347,0.065353,0.028160,0.010640,0.900,0.001287 /)
XEXT_COEFF_550_LKT(5,6)=33.131000 !rg=0.0143374 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,7,1:6)=(/ 507.730000,124.850000,21.674000,3.227500,0.844580,133.090000 /)
XPIZA_LKT(5,7,1:6)=(/ 0.920360,0.917495,0.986676,0.936160,0.900,0.000496 /)
XCGA_LKT(5,7,1:6)=(/ 0.207627,0.097900,0.042210,0.015980,0.900,0.001933 /)
XEXT_COEFF_550_LKT(5,7)=59.291000 !rg=0.0143374 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,8,1:6)=(/ 786.130000,208.770000,38.592000,5.647500,1.098200,133.320000 /)
XPIZA_LKT(5,8,1:6)=(/ 0.945300,0.948681,0.992356,0.963194,0.516198,0.000891 /)
XCGA_LKT(5,8,1:6)=(/ 0.301403,0.143753,0.061837,0.023453,0.007617,0.002847 /)
XEXT_COEFF_550_LKT(5,8)=104.180000 !rg=0.0143374 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,9,1:6)=(/ 1190.400000,339.590000,68.105000,9.964900,1.551100,133.680000 /)
XPIZA_LKT(5,9,1:6)=(/ 0.960987,0.966869,0.995541,0.978875,0.655976,0.001596 /)
XCGA_LKT(5,9,1:6)=(/ 0.404883,0.210123,0.089883,0.034107,0.011097,0.004150 /)
XEXT_COEFF_550_LKT(5,9)=178.680000 !rg=0.0143374 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,10,1:6)=(/ 1712.000000,528.970000,117.590000,17.579000,2.353300,134.240000 /)
XPIZA_LKT(5,10,1:6)=(/ 0.970700,0.977378,0.997318,0.987810,0.771832,0.002840 /)
XCGA_LKT(5,10,1:6)=(/ 0.487977,0.304043,0.130437,0.049387,0.016097,0.006027 /)
XEXT_COEFF_550_LKT(5,10)=293.130000 !rg=0.0143374 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,11,1:6)=(/ 2336.400000,806.850000,196.580000,31.044000,3.789700,135.110000 /)
XPIZA_LKT(5,11,1:6)=(/ 0.976729,0.983931,0.998319,0.992923,0.857018,0.005053 /)
XCGA_LKT(5,11,1:6)=(/ 0.562440,0.412483,0.190210,0.071523,0.023347,0.008750 /)
XEXT_COEFF_550_LKT(5,11)=459.700000 !rg=0.0143374 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,12,1:6)=(/ 3005.000000,1175.100000,312.580000,54.321000,6.350500,136.490000 /)
XPIZA_LKT(5,12,1:6)=(/ 0.980540,0.988051,0.998879,0.995819,0.913551,0.008959 /)
XCGA_LKT(5,12,1:6)=(/ 0.620603,0.495603,0.277963,0.103753,0.033853,0.012713 /)
XEXT_COEFF_550_LKT(5,12)=703.470000 !rg=0.0143374 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,13,1:6)=(/ 3636.700000,1604.100000,479.940000,92.890000,10.899000,138.680000 /)
XPIZA_LKT(5,13,1:6)=(/ 0.982831,0.990514,0.999211,0.997449,0.948689,0.015796 /)
XCGA_LKT(5,13,1:6)=(/ 0.665330,0.567447,0.391047,0.151163,0.049097,0.018480 /)
XEXT_COEFF_550_LKT(5,13)=1022.600000 !rg=0.0143374 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,14,1:6)=(/ 4108.700000,2062.900000,714.670000,151.810000,18.870000,142.160000 /)
XPIZA_LKT(5,14,1:6)=(/ 0.983923,0.992068,0.999423,0.998357,0.969598,0.027546 /)
XCGA_LKT(5,14,1:6)=(/ 0.695550,0.624887,0.484620,0.221587,0.071210,0.026857 /)
XEXT_COEFF_550_LKT(5,14)=1388.800000 !rg=0.0143374 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,15,1:6)=(/ 4318.300000,2492.500000,984.830000,235.340000,32.583000,147.760000 /)
XPIZA_LKT(5,15,1:6)=(/ 0.983951,0.992988,0.999548,0.998868,0.981786,0.047316 /)
XCGA_LKT(5,15,1:6)=(/ 0.712523,0.668250,0.554157,0.323483,0.103447,0.039047 /)
XEXT_COEFF_550_LKT(5,15)=1775.500000 !rg=0.0143374 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,16,1:6)=(/ 4217.400000,2814.300000,1289.800000,357.080000,55.264000,156.790000 /)
XPIZA_LKT(5,16,1:6)=(/ 0.982931,0.993423,0.999627,0.999187,0.988788,0.079476 /)
XCGA_LKT(5,16,1:6)=(/ 0.717540,0.697483,0.617117,0.439187,0.150920,0.056850 /)
XEXT_COEFF_550_LKT(5,16)=2134.300000 !rg=0.0143374 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,17,1:6)=(/ 3880.900000,2950.400000,1578.100000,514.780000,89.661000,171.190000 /)
XPIZA_LKT(5,17,1:6)=(/ 0.980810,0.993408,0.999674,0.999389,0.992722,0.128333 /)
XCGA_LKT(5,17,1:6)=(/ 0.713893,0.713147,0.662007,0.518400,0.221473,0.082940 /)
XEXT_COEFF_550_LKT(5,17)=2390.000000 !rg=0.0143374 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,18,1:6)=(/ 3465.200000,2872.500000,1805.600000,690.010000,137.800000,193.600000 /)
XPIZA_LKT(5,18,1:6)=(/ 0.978108,0.992948,0.999697,0.999507,0.994947,0.195987 /)
XCGA_LKT(5,18,1:6)=(/ 0.711847,0.716557,0.693663,0.586033,0.323967,0.121683 /)
XEXT_COEFF_550_LKT(5,18)=2486.800000 !rg=0.0143374 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,19,1:6)=(/ 3075.600000,2630.000000,1923.700000,873.210000,207.210000,226.800000 /)
XPIZA_LKT(5,19,1:6)=(/ 0.974982,0.992046,0.999701,0.999583,0.996334,0.277839 /)
XCGA_LKT(5,19,1:6)=(/ 0.716870,0.711287,0.711940,0.639553,0.441960,0.180320 /)
XEXT_COEFF_550_LKT(5,19)=2402.600000 !rg=0.0143374 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,20,1:6)=(/ 2667.900000,2345.100000,1901.300000,1037.700000,297.980000,272.230000 /)
XPIZA_LKT(5,20,1:6)=(/ 0.970653,0.990839,0.999684,0.999626,0.997238,0.360990 /)
XCGA_LKT(5,20,1:6)=(/ 0.720617,0.707787,0.716880,0.679190,0.521823,0.268437 /)
XEXT_COEFF_550_LKT(5,20)=2190.800000 !rg=0.0143374 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,1,1:6)=(/ 69.928000,17.062000,1.657300,0.397350,0.547330,132.710000 /)
XPIZA_LKT(6,1,1:6)=(/ 0.500697,0.444930,0.832534,0.900,0.900,0.000032 /)
XCGA_LKT(6,1,1:6)=(/ 0.018983,0.008870,0.003790,0.900,0.900,0.000170 /)
XEXT_COEFF_550_LKT(6,1)=4.531400 !rg=0.0155325 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,2,1:6)=(/ 79.281000,19.073000,2.017400,0.448160,0.552710,132.720000 /)
XPIZA_LKT(6,2,1:6)=(/ 0.556699,0.501841,0.862234,0.900,0.900,0.000040 /)
XCGA_LKT(6,2,1:6)=(/ 0.024090,0.011267,0.004817,0.900,0.900,0.000217 /)
XEXT_COEFF_550_LKT(6,2)=5.527100 !rg=0.0155325 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,3,1:6)=(/ 101.950000,23.971000,2.896600,0.572130,0.565800,132.750000 /)
XPIZA_LKT(6,3,1:6)=(/ 0.650530,0.600916,0.903760,0.900,0.900,0.000061 /)
XCGA_LKT(6,3,1:6)=(/ 0.036623,0.017163,0.007347,0.900,0.900,0.000333 /)
XEXT_COEFF_550_LKT(6,3)=7.957300 !rg=0.0155325 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,4,1:6)=(/ 148.910000,34.299000,4.759300,0.834760,0.593480,132.790000 /)
XPIZA_LKT(6,4,1:6)=(/ 0.755347,0.717836,0.941121,0.900,0.900,0.000104 /)
XCGA_LKT(6,4,1:6)=(/ 0.062110,0.029243,0.012550,0.900,0.900,0.000570 /)
XEXT_COEFF_550_LKT(6,4)=13.102000 !rg=0.0155325 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,5,1:6)=(/ 237.210000,54.561000,8.448500,1.355100,0.648180,132.870000 /)
XPIZA_LKT(6,5,1:6)=(/ 0.841339,0.819401,0.966556,0.849270,0.900,0.000189 /)
XCGA_LKT(6,5,1:6)=(/ 0.103550,0.048967,0.021077,0.007957,0.900,0.000960 /)
XEXT_COEFF_550_LKT(6,5)=23.265000 !rg=0.0155325 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,6,1:6)=(/ 384.070000,91.094000,15.234000,2.314200,0.748830,132.990000 /)
XPIZA_LKT(6,6,1:6)=(/ 0.897769,0.889044,0.981224,0.911307,0.900,0.000346 /)
XCGA_LKT(6,6,1:6)=(/ 0.161693,0.076507,0.032997,0.012480,0.900,0.001510 /)
XEXT_COEFF_550_LKT(6,6)=41.834000 !rg=0.0155325 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,7,1:6)=(/ 608.980000,154.140000,27.437000,4.048700,0.930670,133.170000 /)
XPIZA_LKT(6,7,1:6)=(/ 0.931978,0.932144,0.989390,0.948943,0.900,0.000630 /)
XCGA_LKT(6,7,1:6)=(/ 0.242000,0.114653,0.049447,0.018737,0.900,0.002270 /)
XEXT_COEFF_550_LKT(6,7)=74.727000 !rg=0.0155325 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,8,1:6)=(/ 935.120000,255.710000,48.776000,7.124400,1.253100,133.450000 /)
XPIZA_LKT(6,8,1:6)=(/ 0.952575,0.957296,0.993886,0.970687,0.575311,0.001132 /)
XCGA_LKT(6,8,1:6)=(/ 0.343567,0.168537,0.072437,0.027490,0.008933,0.003340 /)
XEXT_COEFF_550_LKT(6,8)=130.370000 !rg=0.0155325 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,9,1:6)=(/ 1393.800000,409.240000,85.599000,12.605000,1.828700,133.890000 /)
XPIZA_LKT(6,9,1:6)=(/ 0.965586,0.971846,0.996401,0.983187,0.707530,0.002027 /)
XCGA_LKT(6,9,1:6)=(/ 0.439733,0.246190,0.105343,0.039967,0.013013,0.004867 /)
XEXT_COEFF_550_LKT(6,9)=220.330000 !rg=0.0155325 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,10,1:6)=(/ 1953.600000,630.410000,146.000000,22.245000,2.848500,134.550000 /)
XPIZA_LKT(6,10,1:6)=(/ 0.973469,0.980412,0.997801,0.990276,0.810862,0.003605 /)
XCGA_LKT(6,10,1:6)=(/ 0.518137,0.349797,0.153100,0.057863,0.018877,0.007070 /)
XEXT_COEFF_550_LKT(6,10)=354.130000 !rg=0.0155325 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,11,1:6)=(/ 2610.600000,949.990000,239.210000,39.180000,4.673700,135.610000 /)
XPIZA_LKT(6,11,1:6)=(/ 0.978522,0.985882,0.998587,0.994322,0.883499,0.006407 /)
XCGA_LKT(6,11,1:6)=(/ 0.587827,0.449537,0.223607,0.083817,0.027367,0.010267 /)
XEXT_COEFF_550_LKT(6,11)=548.610000 !rg=0.0155325 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,12,1:6)=(/ 3277.100000,1340.500000,373.460000,68.041000,7.924300,137.270000 /)
XPIZA_LKT(6,12,1:6)=(/ 0.981624,0.989190,0.999034,0.996607,0.930242,0.011339 /)
XCGA_LKT(6,12,1:6)=(/ 0.640057,0.524500,0.324103,0.121723,0.039673,0.014913 /)
XEXT_COEFF_550_LKT(6,12)=828.380000 !rg=0.0155325 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,13,1:6)=(/ 3858.400000,1795.100000,569.890000,114.540000,13.683000,139.930000 /)
XPIZA_LKT(6,13,1:6)=(/ 0.983385,0.991251,0.999311,0.997889,0.958735,0.019932 /)
XCGA_LKT(6,13,1:6)=(/ 0.678940,0.593533,0.434817,0.177750,0.057530,0.021673 /)
XEXT_COEFF_550_LKT(6,13)=1164.800000 !rg=0.0155325 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,14,1:6)=(/ 4234.600000,2248.500000,821.820000,182.980000,23.712000,144.170000 /)
XPIZA_LKT(6,14,1:6)=(/ 0.984045,0.992502,0.999483,0.998602,0.975490,0.034585 /)
XCGA_LKT(6,14,1:6)=(/ 0.703933,0.643877,0.513193,0.260813,0.083467,0.031500 /)
XEXT_COEFF_550_LKT(6,14)=1550.800000 !rg=0.0155325 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,15,1:6)=(/ 4314.400000,2641.600000,1109.600000,279.900000,40.719000,150.990000 /)
XPIZA_LKT(6,15,1:6)=(/ 0.983662,0.993205,0.999585,0.999016,0.985179,0.058922 /)
XCGA_LKT(6,15,1:6)=(/ 0.715710,0.681513,0.582793,0.373803,0.121400,0.045813 /)
XEXT_COEFF_550_LKT(6,15)=1931.600000 !rg=0.0155325 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,16,1:6)=(/ 4102.800000,2896.900000,1410.200000,420.440000,68.063000,162.010000 /)
XPIZA_LKT(6,16,1:6)=(/ 0.982173,0.993458,0.999649,0.999284,0.990709,0.097667 /)
XCGA_LKT(6,16,1:6)=(/ 0.716457,0.705397,0.636403,0.476850,0.177570,0.066760 /)
XEXT_COEFF_550_LKT(6,16)=2256.500000 !rg=0.0155325 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,17,1:6)=(/ 3709.900000,2943.100000,1680.500000,584.100000,107.930000,179.440000 /)
XPIZA_LKT(6,17,1:6)=(/ 0.979819,0.993274,0.999686,0.999445,0.993802,0.154532 /)
XCGA_LKT(6,17,1:6)=(/ 0.712757,0.715903,0.676620,0.545923,0.261027,0.097610 /)
XEXT_COEFF_550_LKT(6,17)=2453.200000 !rg=0.0155325 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,18,1:6)=(/ 3309.400000,2785.700000,1868.200000,768.570000,163.550000,206.090000 /)
XPIZA_LKT(6,18,1:6)=(/ 0.976914,0.992623,0.999701,0.999544,0.995595,0.229253 /)
XCGA_LKT(6,18,1:6)=(/ 0.713947,0.714807,0.702417,0.611423,0.375517,0.143793 /)
XEXT_COEFF_550_LKT(6,18)=2471.800000 !rg=0.0155325 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,19,1:6)=(/ 2917.900000,2512.100000,1930.700000,946.540000,244.110000,244.360000 /)
XPIZA_LKT(6,19,1:6)=(/ 0.973321,0.991594,0.999697,0.999603,0.996772,0.313499 /)
XCGA_LKT(6,19,1:6)=(/ 0.718473,0.709800,0.715070,0.657493,0.481050,0.214053 /)
XEXT_COEFF_550_LKT(6,19)=2324.900000 !rg=0.0155325 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,20,1:6)=(/ 2488.100000,2242.200000,1851.900000,1093.000000,337.800000,295.440000 /)
XPIZA_LKT(6,20,1:6)=(/ 0.968717,0.990251,0.999671,0.999635,0.997490,0.393709 /)
XCGA_LKT(6,20,1:6)=(/ 0.721313,0.708163,0.715567,0.691297,0.548780,0.315560 /)
XEXT_COEFF_550_LKT(6,20)=2084.600000 !rg=0.0155325 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,1,1:6)=(/ 79.713000,19.161000,2.032500,0.450290,0.552950,132.720000 /)
XPIZA_LKT(7,1,1:6)=(/ 0.558539,0.503833,0.863219,0.900,0.900,0.000041 /)
XCGA_LKT(7,1,1:6)=(/ 0.022257,0.010403,0.004447,0.900,0.900,0.000200 /)
XEXT_COEFF_550_LKT(7,1)=5.568700 !rg=0.0168272 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,2,1:6)=(/ 91.571000,21.715000,2.490400,0.514880,0.559770,132.740000 /)
XPIZA_LKT(7,2,1:6)=(/ 0.612788,0.560538,0.888187,0.900,0.900,0.000051 /)
XCGA_LKT(7,2,1:6)=(/ 0.028243,0.013213,0.005650,0.900,0.900,0.000257 /)
XEXT_COEFF_550_LKT(7,2)=6.834800 !rg=0.0168272 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,3,1:6)=(/ 120.250000,27.935000,3.608400,0.672520,0.576410,132.770000 /)
XPIZA_LKT(7,3,1:6)=(/ 0.700510,0.655675,0.922558,0.900,0.900,0.000077 /)
XCGA_LKT(7,3,1:6)=(/ 0.042907,0.020123,0.008617,0.900,0.900,0.000390 /)
XEXT_COEFF_550_LKT(7,3)=9.924900 !rg=0.0168272 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,4,1:6)=(/ 179.220000,41.035000,5.976800,1.006400,0.611560,132.820000 /)
XPIZA_LKT(7,4,1:6)=(/ 0.793908,0.762433,0.952960,0.797441,0.900,0.000132 /)
XCGA_LKT(7,4,1:6)=(/ 0.072660,0.034263,0.014717,0.005547,0.900,0.000670 /)
XEXT_COEFF_550_LKT(7,4)=16.463000 !rg=0.0168272 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,5,1:6)=(/ 288.050000,66.634000,10.665000,1.668200,0.681070,132.910000 /)
XPIZA_LKT(7,5,1:6)=(/ 0.867064,0.850657,0.973385,0.877330,0.900,0.000240 /)
XCGA_LKT(7,5,1:6)=(/ 0.120980,0.057327,0.024703,0.009330,0.900,0.001127 /)
XEXT_COEFF_550_LKT(7,5)=29.355000 !rg=0.0168272 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,6,1:6)=(/ 464.090000,112.380000,19.278000,2.887600,0.808980,133.050000 /)
XPIZA_LKT(7,6,1:6)=(/ 0.913575,0.908869,0.985065,0.928733,0.900,0.000440 /)
XCGA_LKT(7,6,1:6)=(/ 0.188730,0.089550,0.038660,0.014637,0.900,0.001770 /)
XEXT_COEFF_550_LKT(7,6)=52.810000 !rg=0.0168272 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,7,1:6)=(/ 728.130000,189.910000,34.720000,5.092600,1.040100,133.270000 /)
XPIZA_LKT(7,7,1:6)=(/ 0.941556,0.943989,0.991539,0.959254,0.489460,0.000801 /)
XCGA_LKT(7,7,1:6)=(/ 0.280070,0.134273,0.057917,0.021967,0.007133,0.002663 /)
XEXT_COEFF_550_LKT(7,7)=93.969000 !rg=0.0168272 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,8,1:6)=(/ 1108.000000,311.260000,61.536000,8.999900,1.449900,133.610000 /)
XPIZA_LKT(7,8,1:6)=(/ 0.958660,0.964177,0.995093,0.976669,0.632281,0.001439 /)
XCGA_LKT(7,8,1:6)=(/ 0.383357,0.197540,0.084857,0.032217,0.010480,0.003917 /)
XEXT_COEFF_550_LKT(7,8)=162.280000 !rg=0.0168272 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,9,1:6)=(/ 1612.200000,490.340000,107.120000,15.951000,2.181700,134.130000 /)
XPIZA_LKT(7,9,1:6)=(/ 0.969295,0.975866,0.997078,0.986612,0.754180,0.002574 /)
XCGA_LKT(7,9,1:6)=(/ 0.471153,0.287050,0.123510,0.046827,0.015263,0.005713 /)
XEXT_COEFF_550_LKT(7,9)=269.400000 !rg=0.0168272 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,10,1:6)=(/ 2216.000000,749.780000,179.860000,28.128000,3.477700,134.930000 /)
XPIZA_LKT(7,10,1:6)=(/ 0.975800,0.982961,0.998179,0.992228,0.844478,0.004574 /)
XCGA_LKT(7,10,1:6)=(/ 0.547993,0.393780,0.179810,0.067793,0.022130,0.008293 /)
XEXT_COEFF_550_LKT(7,10)=424.960000 !rg=0.0168272 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,11,1:6)=(/ 2882.400000,1102.900000,288.480000,49.326000,5.796200,136.200000 /)
XPIZA_LKT(7,11,1:6)=(/ 0.979962,0.987450,0.998799,0.995426,0.905535,0.008118 /)
XCGA_LKT(7,11,1:6)=(/ 0.609637,0.480797,0.262493,0.098253,0.032077,0.012043 /)
XEXT_COEFF_550_LKT(7,11)=653.360000 !rg=0.0168272 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,12,1:6)=(/ 3529.900000,1519.700000,445.170000,84.751000,9.919000,138.220000 /)
XPIZA_LKT(7,12,1:6)=(/ 0.982496,0.990134,0.999162,0.997227,0.943828,0.014336 /)
XCGA_LKT(7,12,1:6)=(/ 0.657133,0.554073,0.372273,0.142913,0.046487,0.017493 /)
XEXT_COEFF_550_LKT(7,12)=960.570000 !rg=0.0168272 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,13,1:6)=(/ 4043.000000,1981.300000,669.790000,139.930000,17.195000,141.450000 /)
XPIZA_LKT(7,13,1:6)=(/ 0.983798,0.991842,0.999393,0.998234,0.966803,0.025107 /)
XCGA_LKT(7,13,1:6)=(/ 0.690257,0.615150,0.470487,0.209190,0.067410,0.025420 /)
XEXT_COEFF_550_LKT(7,13)=1318.900000 !rg=0.0168272 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,14,1:6)=(/ 4306.200000,2422.400000,933.900000,218.780000,29.751000,146.620000 /)
XPIZA_LKT(7,14,1:6)=(/ 0.984019,0.992858,0.999530,0.998798,0.980180,0.043295 /)
XCGA_LKT(7,14,1:6)=(/ 0.709787,0.660880,0.541597,0.305867,0.097863,0.036947 /)
XEXT_COEFF_550_LKT(7,14)=1707.400000 !rg=0.0168272 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,15,1:6)=(/ 4262.400000,2767.100000,1235.400000,332.680000,50.633000,154.940000 /)
XPIZA_LKT(7,15,1:6)=(/ 0.983182,0.993368,0.999616,0.999141,0.987862,0.073023 /)
XCGA_LKT(7,15,1:6)=(/ 0.716933,0.692550,0.607343,0.421913,0.142597,0.053763 /)
XEXT_COEFF_550_LKT(7,15)=2074.800000 !rg=0.0168272 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,16,1:6)=(/ 3954.200000,2941.500000,1530.300000,486.930000,83.056000,168.350000 /)
XPIZA_LKT(7,16,1:6)=(/ 0.981342,0.993439,0.999667,0.999362,0.992217,0.119116 /)
XCGA_LKT(7,16,1:6)=(/ 0.715340,0.710863,0.654550,0.506863,0.209167,0.078453 /)
XEXT_COEFF_550_LKT(7,16)=2355.600000 !rg=0.0168272 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,17,1:6)=(/ 3549.300000,2900.900000,1772.300000,658.690000,128.820000,189.330000 /)
XPIZA_LKT(7,17,1:6)=(/ 0.978697,0.993054,0.999694,0.999490,0.994659,0.184010 /)
XCGA_LKT(7,17,1:6)=(/ 0.712910,0.716237,0.688757,0.574870,0.306750,0.115043 /)
XEXT_COEFF_550_LKT(7,17)=2482.300000 !rg=0.0168272 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,18,1:6)=(/ 3157.700000,2682.400000,1911.800000,843.060000,194.110000,220.630000 /)
XPIZA_LKT(7,18,1:6)=(/ 0.975583,0.992230,0.999702,0.999573,0.996145,0.264243 /)
XCGA_LKT(7,18,1:6)=(/ 0.715553,0.712760,0.709257,0.631787,0.425430,0.170287 /)
XEXT_COEFF_550_LKT(7,18)=2428.900000 !rg=0.0168272 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,19,1:6)=(/ 2738.400000,2400.500000,1913.000000,1013.300000,282.920000,264.180000 /)
XPIZA_LKT(7,19,1:6)=(/ 0.971721,0.991051,0.999688,0.999620,0.997125,0.348294 /)
XCGA_LKT(7,19,1:6)=(/ 0.720637,0.707823,0.716433,0.673427,0.511433,0.253917 /)
XEXT_COEFF_550_LKT(7,19)=2229.000000 !rg=0.0168272 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,20,1:6)=(/ 2319.400000,2131.800000,1788.000000,1136.900000,380.250000,321.150000 /)
XPIZA_LKT(7,20,1:6)=(/ 0.967574,0.989747,0.999653,0.999642,0.997695,0.424862 /)
XCGA_LKT(7,20,1:6)=(/ 0.725443,0.711243,0.712723,0.701087,0.577450,0.364270 /)
XEXT_COEFF_550_LKT(7,20)=1990.600000 !rg=0.0168272 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,1,1:6)=(/ 92.122000,21.826000,2.509600,0.517600,0.560070,132.740000 /)
XPIZA_LKT(8,1,1:6)=(/ 0.614530,0.562476,0.889005,0.900,0.900,0.000052 /)
XCGA_LKT(8,1,1:6)=(/ 0.026097,0.012203,0.005217,0.900,0.900,0.000237 /)
XEXT_COEFF_550_LKT(8,1)=6.887700 !rg=0.0182298 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,2,1:6)=(/ 107.150000,25.070000,3.091900,0.599720,0.568750,132.760000 /)
XPIZA_LKT(8,2,1:6)=(/ 0.665700,0.617419,0.909736,0.900,0.900,0.000065 /)
XCGA_LKT(8,2,1:6)=(/ 0.033103,0.015497,0.006630,0.900,0.900,0.000300 /)
XEXT_COEFF_550_LKT(8,2)=8.497600 !rg=0.0182298 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,3,1:6)=(/ 143.390000,32.968000,4.513500,0.800140,0.589870,132.790000 /)
XPIZA_LKT(8,3,1:6)=(/ 0.745717,0.706387,0.937913,0.900,0.900,0.000098 /)
XCGA_LKT(8,3,1:6)=(/ 0.050263,0.023590,0.010107,0.900,0.900,0.000460 /)
XEXT_COEFF_550_LKT(8,3)=12.426000 !rg=0.0182298 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,4,1:6)=(/ 217.110000,49.570000,7.524600,1.224700,0.634530,132.860000 /)
XPIZA_LKT(8,4,1:6)=(/ 0.827214,0.801684,0.962493,0.833298,0.900,0.000168 /)
XCGA_LKT(8,4,1:6)=(/ 0.084967,0.040130,0.017253,0.006507,0.900,0.000787 /)
XEXT_COEFF_550_LKT(8,4)=20.732000 !rg=0.0182298 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,5,1:6)=(/ 350.010000,81.829000,13.481000,2.066100,0.722860,132.960000 /)
XPIZA_LKT(8,5,1:6)=(/ 0.888483,0.877012,0.978830,0.900749,0.900,0.000306 /)
XCGA_LKT(8,5,1:6)=(/ 0.141223,0.067093,0.028950,0.010943,0.900,0.001323 /)
XEXT_COEFF_550_LKT(8,5)=37.066000 !rg=0.0182298 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,6,1:6)=(/ 558.980000,138.770000,24.402000,3.616500,0.885410,133.130000 /)
XPIZA_LKT(8,6,1:6)=(/ 0.926538,0.925103,0.988112,0.942923,0.900,0.000559 /)
XCGA_LKT(8,6,1:6)=(/ 0.219660,0.104790,0.045287,0.017160,0.900,0.002080 /)
XEXT_COEFF_550_LKT(8,6)=66.602000 !rg=0.0182298 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,7,1:6)=(/ 868.460000,233.090000,43.899000,6.419000,1.179200,133.390000 /)
XPIZA_LKT(8,7,1:6)=(/ 0.949527,0.953512,0.993238,0.967531,0.549006,0.001017 /)
XCGA_LKT(8,7,1:6)=(/ 0.319873,0.157230,0.067830,0.025750,0.008367,0.003127 /)
XEXT_COEFF_550_LKT(8,7)=117.760000 !rg=0.0182298 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,8,1:6)=(/ 1301.700000,376.330000,77.422000,11.380000,1.700100,133.790000 /)
XPIZA_LKT(8,8,1:6)=(/ 0.963660,0.969682,0.996046,0.981433,0.685715,0.001828 /)
XCGA_LKT(8,8,1:6)=(/ 0.418963,0.231197,0.099417,0.037753,0.012290,0.004597 /)
XEXT_COEFF_550_LKT(8,8)=200.640000 !rg=0.0182298 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,9,1:6)=(/ 1847.600000,585.620000,133.290000,20.185000,2.630300,134.420000 /)
XPIZA_LKT(8,9,1:6)=(/ 0.972312,0.979175,0.997609,0.989326,0.795461,0.003268 /)
XCGA_LKT(8,9,1:6)=(/ 0.501963,0.330997,0.144877,0.054863,0.017897,0.006703 /)
XEXT_COEFF_550_LKT(8,9)=326.600000 !rg=0.0182298 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,10,1:6)=(/ 2488.700000,885.720000,219.610000,35.519000,4.277300,135.390000 /)
XPIZA_LKT(8,10,1:6)=(/ 0.977754,0.985080,0.998476,0.993771,0.872974,0.005801 /)
XCGA_LKT(8,10,1:6)=(/ 0.574687,0.432093,0.211223,0.079437,0.025943,0.009730 /)
XEXT_COEFF_550_LKT(8,10)=508.150000 !rg=0.0182298 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,11,1:6)=(/ 3157.700000,1263.700000,345.540000,61.872000,7.220900,136.930000 /)
XPIZA_LKT(8,11,1:6)=(/ 0.981147,0.988697,0.998969,0.996296,0.923681,0.010279 /)
XCGA_LKT(8,11,1:6)=(/ 0.629940,0.510113,0.306527,0.115230,0.037590,0.014127 /)
XEXT_COEFF_550_LKT(8,11)=772.300000 !rg=0.0182298 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,12,1:6)=(/ 3763.000000,1709.000000,529.640000,104.800000,12.443000,139.370000 /)
XPIZA_LKT(8,12,1:6)=(/ 0.983133,0.990934,0.999270,0.997714,0.954816,0.018101 /)
XCGA_LKT(8,12,1:6)=(/ 0.671727,0.581473,0.417280,0.167943,0.054470,0.020517 /)
XEXT_COEFF_550_LKT(8,12)=1098.900000 !rg=0.0182298 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,13,1:6)=(/ 4189.000000,2167.600000,774.580000,169.260000,21.609000,143.300000 /)
XPIZA_LKT(8,13,1:6)=(/ 0.983994,0.992315,0.999459,0.998505,0.973255,0.031555 /)
XCGA_LKT(8,13,1:6)=(/ 0.699540,0.634757,0.500220,0.246180,0.079003,0.029813 /)
XEXT_COEFF_550_LKT(8,13)=1480.300000 !rg=0.0182298 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,14,1:6)=(/ 4328.700000,2578.300000,1055.400000,260.490000,37.225000,149.610000 /)
XPIZA_LKT(8,14,1:6)=(/ 0.983810,0.993111,0.999570,0.998957,0.983903,0.054003 /)
XCGA_LKT(8,14,1:6)=(/ 0.714053,0.675093,0.570737,0.355067,0.114810,0.043347 /)
XEXT_COEFF_550_LKT(8,14)=1864.500000 !rg=0.0182298 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,15,1:6)=(/ 4163.400000,2862.300000,1356.400000,393.240000,62.530000,159.760000 /)
XPIZA_LKT(8,15,1:6)=(/ 0.982556,0.993439,0.999640,0.999246,0.989975,0.089970 /)
XCGA_LKT(8,15,1:6)=(/ 0.716933,0.701290,0.627557,0.462317,0.167677,0.063120 /)
XEXT_COEFF_550_LKT(8,15)=2202.900000 !rg=0.0182298 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,16,1:6)=(/ 3800.500000,2952.400000,1638.400000,554.910000,100.330000,176.020000 /)
XPIZA_LKT(8,16,1:6)=(/ 0.980301,0.993330,0.999681,0.999423,0.993401,0.143973 /)
XCGA_LKT(8,16,1:6)=(/ 0.713963,0.714413,0.670277,0.534377,0.246497,0.092287 /)
XEXT_COEFF_550_LKT(8,16)=2429.400000 !rg=0.0182298 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,17,1:6)=(/ 3391.600000,2827.700000,1843.900000,737.240000,153.030000,201.050000 /)
XPIZA_LKT(8,17,1:6)=(/ 0.977460,0.992765,0.999700,0.999530,0.995355,0.216319 /)
XCGA_LKT(8,17,1:6)=(/ 0.713307,0.715503,0.698493,0.601793,0.357180,0.135840 /)
XEXT_COEFF_550_LKT(8,17)=2482.300000 !rg=0.0182298 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,18,1:6)=(/ 2990.600000,2569.700000,1928.400000,916.620000,229.450000,237.310000 /)
XPIZA_LKT(8,18,1:6)=(/ 0.974279,0.991758,0.999698,0.999595,0.996614,0.299851 /)
XCGA_LKT(8,18,1:6)=(/ 0.719287,0.709727,0.713270,0.650117,0.467503,0.202030 /)
XEXT_COEFF_550_LKT(8,18)=2358.300000 !rg=0.0182298 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,19,1:6)=(/ 2563.100000,2287.000000,1873.500000,1071.100000,322.200000,286.410000 /)
XPIZA_LKT(8,19,1:6)=(/ 0.969838,0.990569,0.999676,0.999631,0.997400,0.381568 /)
XCGA_LKT(8,19,1:6)=(/ 0.722783,0.709130,0.715803,0.686370,0.538447,0.299423 /)
XEXT_COEFF_550_LKT(8,19)=2134.300000 !rg=0.0182298 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,20,1:6)=(/ 2171.300000,2025.900000,1710.100000,1167.200000,425.410000,348.650000 /)
XPIZA_LKT(8,20,1:6)=(/ 0.965219,0.989179,0.999634,0.999643,0.997870,0.454588 /)
XCGA_LKT(8,20,1:6)=(/ 0.726970,0.714187,0.709613,0.708480,0.604927,0.409123 /)
XEXT_COEFF_550_LKT(8,20)=1897.400000 !rg=0.0182298 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,1,1:6)=(/ 107.860000,25.212000,3.116300,0.603170,0.569130,132.760000 /)
XPIZA_LKT(9,1,1:6)=(/ 0.667316,0.619256,0.910410,0.900,0.900,0.000066 /)
XCGA_LKT(9,1,1:6)=(/ 0.030593,0.014313,0.006120,0.900,0.900,0.000277 /)
XEXT_COEFF_550_LKT(9,1)=8.565000 !rg=0.0197494 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,2,1:6)=(/ 126.890000,29.332000,3.856800,0.707580,0.580140,132.780000 /)
XPIZA_LKT(9,2,1:6)=(/ 0.714390,0.671073,0.927448,0.900,0.900,0.000083 /)
XCGA_LKT(9,2,1:6)=(/ 0.038800,0.018173,0.007777,0.900,0.900,0.000353 /)
XEXT_COEFF_550_LKT(9,2)=10.612000 !rg=0.0197494 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,3,1:6)=(/ 172.600000,39.357000,5.664500,0.962420,0.606980,132.820000 /)
XPIZA_LKT(9,3,1:6)=(/ 0.785747,0.752243,0.950365,0.900,0.900,0.000125 /)
XCGA_LKT(9,3,1:6)=(/ 0.058870,0.027653,0.011857,0.900,0.900,0.000540 /)
XEXT_COEFF_550_LKT(9,3)=15.606000 !rg=0.0197494 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,4,1:6)=(/ 264.180000,60.372000,9.492300,1.502300,0.663720,132.900000 /)
XPIZA_LKT(9,4,1:6)=(/ 0.855528,0.835592,0.970135,0.863864,0.900,0.000213 /)
XCGA_LKT(9,4,1:6)=(/ 0.099307,0.046993,0.020227,0.007633,0.900,0.000923 /)
XEXT_COEFF_550_LKT(9,4)=26.151000 !rg=0.0197494 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,5,1:6)=(/ 424.840000,100.880000,17.056000,2.572200,0.775960,133.030000 /)
XPIZA_LKT(9,5,1:6)=(/ 0.906170,0.898956,0.983163,0.920079,0.900,0.000388 /)
XCGA_LKT(9,5,1:6)=(/ 0.164653,0.078493,0.033920,0.012837,0.900,0.001553 /)
XEXT_COEFF_550_LKT(9,5)=46.810000 !rg=0.0197494 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,6,1:6)=(/ 671.210000,171.190000,30.886000,4.543200,0.982560,133.230000 /)
XPIZA_LKT(9,6,1:6)=(/ 0.937210,0.938288,0.990526,0.954404,0.900,0.000711 /)
XCGA_LKT(9,6,1:6)=(/ 0.254133,0.122593,0.053040,0.020120,0.900,0.002440 /)
XEXT_COEFF_550_LKT(9,6)=83.846000 !rg=0.0197494 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,7,1:6)=(/ 1031.600000,284.580000,55.418000,8.104000,1.356000,133.540000 /)
XPIZA_LKT(9,7,1:6)=(/ 0.956156,0.961140,0.994580,0.974150,0.607124,0.001293 /)
XCGA_LKT(9,7,1:6)=(/ 0.358603,0.184023,0.079437,0.030180,0.009813,0.003670 /)
XEXT_COEFF_550_LKT(9,7)=146.900000 !rg=0.0197494 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,8,1:6)=(/ 1512.600000,452.420000,97.033000,14.398000,2.018100,134.020000 /)
XPIZA_LKT(9,8,1:6)=(/ 0.967711,0.974124,0.996796,0.985218,0.734574,0.002321 /)
XCGA_LKT(9,8,1:6)=(/ 0.451550,0.269453,0.116497,0.044233,0.014413,0.005393 /)
XEXT_COEFF_550_LKT(9,8)=246.150000 !rg=0.0197494 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,9,1:6)=(/ 2104.200000,697.840000,164.650000,25.528000,3.200400,134.770000 /)
XPIZA_LKT(9,9,1:6)=(/ 0.974839,0.981940,0.998027,0.991475,0.831282,0.004148 /)
XCGA_LKT(9,9,1:6)=(/ 0.532550,0.374373,0.170010,0.064273,0.020983,0.007863 /)
XEXT_COEFF_550_LKT(9,9)=393.180000 !rg=0.0197494 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,10,1:6)=(/ 2761.700000,1033.100000,265.760000,44.753000,5.292800,135.940000 /)
XPIZA_LKT(9,10,1:6)=(/ 0.979336,0.986792,0.998710,0.994990,0.896805,0.007353 /)
XCGA_LKT(9,10,1:6)=(/ 0.597693,0.464683,0.247853,0.093097,0.030407,0.011417 /)
XEXT_COEFF_550_LKT(9,10)=606.300000 !rg=0.0197494 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,11,1:6)=(/ 3420.800000,1437.500000,412.550000,77.217000,9.027300,137.800000 /)
XPIZA_LKT(9,11,1:6)=(/ 0.982121,0.989722,0.999108,0.996981,0.938497,0.013002 /)
XCGA_LKT(9,11,1:6)=(/ 0.648160,0.539950,0.353537,0.135227,0.044050,0.016573 /)
XEXT_COEFF_550_LKT(9,11)=900.230000 !rg=0.0197494 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,12,1:6)=(/ 3967.200000,1896.400000,625.120000,128.460000,15.631000,140.770000 /)
XPIZA_LKT(9,12,1:6)=(/ 0.983627,0.991584,0.999359,0.998095,0.963656,0.022818 /)
XCGA_LKT(9,12,1:6)=(/ 0.684167,0.604430,0.455017,0.197523,0.063827,0.024067 /)
XEXT_COEFF_550_LKT(9,12)=1248.200000 !rg=0.0197494 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,13,1:6)=(/ 4285.600000,2348.400000,883.740000,202.960000,27.128000,145.560000 /)
XPIZA_LKT(9,13,1:6)=(/ 0.984056,0.992709,0.999510,0.998720,0.978399,0.039552 /)
XCGA_LKT(9,13,1:6)=(/ 0.706593,0.652790,0.528647,0.289000,0.092613,0.034970 /)
XEXT_COEFF_550_LKT(9,13)=1638.400000 !rg=0.0197494 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,14,1:6)=(/ 4299.000000,2716.700000,1181.400000,309.880000,46.372000,153.250000 /)
XPIZA_LKT(9,14,1:6)=(/ 0.983422,0.993300,0.999603,0.999092,0.986849,0.067060 /)
XCGA_LKT(9,14,1:6)=(/ 0.716083,0.687147,0.596867,0.403967,0.134797,0.050863 /)
XEXT_COEFF_550_LKT(9,14)=2014.300000 !rg=0.0197494 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,15,1:6)=(/ 4037.700000,2924.800000,1477.700000,458.190000,76.560000,165.620000 /)
XPIZA_LKT(9,15,1:6)=(/ 0.981725,0.993450,0.999659,0.999332,0.991636,0.110072 /)
XCGA_LKT(9,15,1:6)=(/ 0.715843,0.707960,0.646267,0.494343,0.197400,0.074157 /)
XEXT_COEFF_550_LKT(9,15)=2313.800000 !rg=0.0197494 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,16,1:6)=(/ 3635.800000,2924.100000,1734.900000,627.220000,120.090000,185.240000 /)
XPIZA_LKT(9,16,1:6)=(/ 0.979221,0.993154,0.999691,0.999472,0.994337,0.172184 /)
XCGA_LKT(9,16,1:6)=(/ 0.713097,0.715797,0.683270,0.563123,0.289980,0.108707 /)
XEXT_COEFF_550_LKT(9,16)=2472.900000 !rg=0.0197494 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,17,1:6)=(/ 3234.300000,2731.700000,1896.200000,812.960000,181.690000,214.770000 /)
XPIZA_LKT(9,17,1:6)=(/ 0.976229,0.992402,0.999701,0.999562,0.995942,0.250698 /)
XCGA_LKT(9,17,1:6)=(/ 0.716197,0.713300,0.706097,0.623607,0.407970,0.160740 /)
XEXT_COEFF_550_LKT(9,17)=2450.000000 !rg=0.0197494 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,18,1:6)=(/ 2823.300000,2451.700000,1921.400000,986.750000,267.540000,256.210000 /)
XPIZA_LKT(9,18,1:6)=(/ 0.972488,0.991309,0.999692,0.999614,0.996998,0.334965 /)
XCGA_LKT(9,18,1:6)=(/ 0.720897,0.709257,0.715693,0.667010,0.500200,0.239733 /)
XEXT_COEFF_550_LKT(9,18)=2274.700000 !rg=0.0197494 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,19,1:6)=(/ 2389.400000,2184.400000,1813.800000,1120.400000,363.420000,311.170000 /)
XPIZA_LKT(9,19,1:6)=(/ 0.968167,0.990031,0.999661,0.999640,0.997619,0.413255 /)
XCGA_LKT(9,19,1:6)=(/ 0.723710,0.711427,0.713653,0.697067,0.566693,0.347867 /)
XEXT_COEFF_550_LKT(9,19)=2033.000000 !rg=0.0197494 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,20,1:6)=(/ 2021.100000,1910.500000,1632.400000,1182.800000,469.010000,376.320000 /)
XPIZA_LKT(9,20,1:6)=(/ 0.962927,0.988427,0.999611,0.999641,0.998013,0.482388 /)
XCGA_LKT(9,20,1:6)=(/ 0.728670,0.715513,0.707560,0.713653,0.627227,0.446527 /)
XEXT_COEFF_550_LKT(9,20)=1810.600000 !rg=0.0197494 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,1,1:6)=(/ 127.800000,29.512000,3.887700,0.711960,0.580620,132.780000 /)
XPIZA_LKT(10,1,1:6)=(/ 0.715855,0.672760,0.927994,0.900,0.900,0.000084 /)
XCGA_LKT(10,1,1:6)=(/ 0.035863,0.016787,0.007180,0.900,0.900,0.000327 /)
XEXT_COEFF_550_LKT(10,1)=10.698000 !rg=0.0213956 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,2,1:6)=(/ 151.880000,34.745000,4.829400,0.844720,0.594610,132.810000 /)
XPIZA_LKT(10,2,1:6)=(/ 0.758177,0.720407,0.941883,0.900,0.900,0.000105 /)
XCGA_LKT(10,2,1:6)=(/ 0.045473,0.021310,0.009123,0.900,0.900,0.000413 /)
XEXT_COEFF_550_LKT(10,2)=13.301000 !rg=0.0213956 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,3,1:6)=(/ 209.350000,47.465000,7.128100,1.168800,0.628700,132.860000 /)
XPIZA_LKT(10,3,1:6)=(/ 0.820530,0.792817,0.960404,0.825311,0.900,0.000158 /)
XCGA_LKT(10,3,1:6)=(/ 0.068947,0.032407,0.013903,0.005237,0.900,0.000633 /)
XEXT_COEFF_550_LKT(10,3)=19.648000 !rg=0.0213956 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,4,1:6)=(/ 322.280000,74.016000,11.993000,1.855300,0.700800,132.950000 /)
XPIZA_LKT(10,4,1:6)=(/ 0.879278,0.864416,0.976240,0.889539,0.900,0.000271 /)
XCGA_LKT(10,4,1:6)=(/ 0.116003,0.055013,0.023710,0.008953,0.900,0.001083 /)
XEXT_COEFF_550_LKT(10,4)=33.024000 !rg=0.0213956 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,5,1:6)=(/ 514.500000,124.620000,21.590000,3.215600,0.843450,133.100000 /)
XPIZA_LKT(10,5,1:6)=(/ 0.920711,0.917029,0.986603,0.935885,0.900,0.000494 /)
XCGA_LKT(10,5,1:6)=(/ 0.191530,0.091797,0.039733,0.015050,0.900,0.001823 /)
XEXT_COEFF_550_LKT(10,5)=59.088000 !rg=0.0213956 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,6,1:6)=(/ 803.510000,210.630000,39.066000,5.720900,1.106000,133.340000 /)
XPIZA_LKT(10,6,1:6)=(/ 0.946035,0.948929,0.992435,0.963640,0.519508,0.000903 /)
XCGA_LKT(10,6,1:6)=(/ 0.290857,0.143353,0.062107,0.023583,0.007660,0.002863 /)
XEXT_COEFF_550_LKT(10,6)=105.260000 !rg=0.0213956 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,7,1:6)=(/ 1216.100000,345.350000,69.794000,10.243000,1.580700,133.710000 /)
XPIZA_LKT(10,7,1:6)=(/ 0.961599,0.967254,0.995639,0.979427,0.662292,0.001642 /)
XCGA_LKT(10,7,1:6)=(/ 0.394610,0.215063,0.093030,0.035363,0.011513,0.004307 /)
XEXT_COEFF_550_LKT(10,7)=182.140000 !rg=0.0213956 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,8,1:6)=(/ 1741.800000,541.860000,120.980000,18.219000,2.422400,134.290000 /)
XPIZA_LKT(10,8,1:6)=(/ 0.971014,0.977762,0.997386,0.988220,0.778215,0.002947 /)
XCGA_LKT(10,8,1:6)=(/ 0.483357,0.311047,0.136547,0.051820,0.016903,0.006330 /)
XEXT_COEFF_550_LKT(10,8)=299.540000 !rg=0.0213956 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,9,1:6)=(/ 2374.600000,826.700000,201.690000,32.249000,3.924800,135.190000 /)
XPIZA_LKT(10,9,1:6)=(/ 0.976960,0.984237,0.998355,0.993175,0.861836,0.005261 /)
XCGA_LKT(10,9,1:6)=(/ 0.560460,0.413337,0.199527,0.075303,0.024597,0.009223 /)
XEXT_COEFF_550_LKT(10,9)=471.320000 !rg=0.0213956 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,10,1:6)=(/ 3039.400000,1189.400000,319.270000,56.205000,6.581800,136.600000 /)
XPIZA_LKT(10,10,1:6)=(/ 0.980633,0.988159,0.998897,0.995950,0.916510,0.009313 /)
XCGA_LKT(10,10,1:6)=(/ 0.618953,0.494717,0.289623,0.109147,0.035637,0.013390 /)
XEXT_COEFF_550_LKT(10,10)=718.970000 !rg=0.0213956 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,11,1:6)=(/ 3664.700000,1623.800000,491.670000,95.726000,11.315000,138.860000 /)
XPIZA_LKT(10,11,1:6)=(/ 0.982852,0.990589,0.999225,0.997518,0.950508,0.016427 /)
XCGA_LKT(10,11,1:6)=(/ 0.663820,0.568417,0.398993,0.158810,0.051613,0.019437 /)
XEXT_COEFF_550_LKT(10,11)=1034.700000 !rg=0.0213956 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,12,1:6)=(/ 4132.700000,2083.100000,726.940000,155.940000,19.642000,142.480000 /)
XPIZA_LKT(10,12,1:6)=(/ 0.983912,0.992103,0.999431,0.998394,0.970737,0.028706 /)
XCGA_LKT(10,12,1:6)=(/ 0.694420,0.624770,0.486200,0.232353,0.074793,0.028227 /)
XEXT_COEFF_550_LKT(10,12)=1407.600000 !rg=0.0213956 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,13,1:6)=(/ 4334.900000,2511.400000,1001.600000,242.060000,33.976000,148.310000 /)
XPIZA_LKT(10,13,1:6)=(/ 0.983919,0.993002,0.999553,0.998893,0.982486,0.049409 /)
XCGA_LKT(10,13,1:6)=(/ 0.711717,0.668037,0.557990,0.336603,0.108620,0.041023 /)
XEXT_COEFF_550_LKT(10,13)=1795.600000 !rg=0.0213956 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,14,1:6)=(/ 4224.700000,2824.400000,1303.600000,367.380000,57.413000,157.700000 /)
XPIZA_LKT(10,14,1:6)=(/ 0.982865,0.993411,0.999630,0.999205,0.989173,0.082824 /)
XCGA_LKT(10,14,1:6)=(/ 0.716897,0.696793,0.618300,0.446763,0.158420,0.059707 /)
XEXT_COEFF_550_LKT(10,14)=2148.200000 !rg=0.0213956 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,15,1:6)=(/ 3884.500000,2951.500000,1591.400000,524.900000,92.821000,172.710000 /)
XPIZA_LKT(10,15,1:6)=(/ 0.980828,0.993377,0.999675,0.999399,0.992940,0.133530 /)
XCGA_LKT(10,15,1:6)=(/ 0.714683,0.712227,0.663087,0.522300,0.232553,0.087200 /)
XEXT_COEFF_550_LKT(10,15)=2398.000000 !rg=0.0213956 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,16,1:6)=(/ 3475.000000,2865.200000,1815.800000,704.860000,142.880000,196.200000 /)
XPIZA_LKT(10,16,1:6)=(/ 0.978150,0.992898,0.999698,0.999515,0.995091,0.203415 /)
XCGA_LKT(10,16,1:6)=(/ 0.714583,0.715607,0.694033,0.591230,0.338823,0.128267 /)
XEXT_COEFF_550_LKT(10,16)=2487.000000 !rg=0.0213956 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,17,1:6)=(/ 3077.400000,2621.700000,1923.200000,886.540000,215.290000,230.580000 /)
XPIZA_LKT(10,17,1:6)=(/ 0.974904,0.992002,0.999700,0.999586,0.996443,0.286106 /)
XCGA_LKT(10,17,1:6)=(/ 0.718407,0.711863,0.711247,0.642410,0.452753,0.190570 /)
XEXT_COEFF_550_LKT(10,17)=2393.000000 !rg=0.0213956 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,18,1:6)=(/ 2639.600000,2345.900000,1891.100000,1047.200000,306.340000,277.440000 /)
XPIZA_LKT(10,18,1:6)=(/ 0.970767,0.990774,0.999681,0.999627,0.997299,0.368771 /)
XCGA_LKT(10,18,1:6)=(/ 0.721950,0.709037,0.715570,0.680913,0.527773,0.283347 /)
XEXT_COEFF_550_LKT(10,18)=2175.500000 !rg=0.0213956 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,19,1:6)=(/ 2230.900000,2076.000000,1744.400000,1155.300000,407.910000,338.050000 /)
XPIZA_LKT(10,19,1:6)=(/ 0.966192,0.989392,0.999641,0.999643,0.997806,0.443533 /)
XCGA_LKT(10,19,1:6)=(/ 0.726027,0.713207,0.710993,0.705217,0.595110,0.394160 /)
XEXT_COEFF_550_LKT(10,19)=1940.300000 !rg=0.0213956 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,20,1:6)=(/ 1882.100000,1789.700000,1557.700000,1183.100000,510.560000,402.780000 /)
XPIZA_LKT(10,20,1:6)=(/ 0.959695,0.987687,0.999585,0.999634,0.998124,0.507321 /)
XCGA_LKT(10,20,1:6)=(/ 0.731343,0.716537,0.705047,0.716240,0.645837,0.477590 /)
XEXT_COEFF_550_LKT(10,20)=1723.800000 !rg=0.0213956 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,1,1:6)=(/ 153.080000,34.974000,4.868700,0.850290,0.595220,132.810000 /)
XPIZA_LKT(11,1,1:6)=(/ 0.759496,0.721926,0.942324,0.900,0.900,0.000106 /)
XCGA_LKT(11,1,1:6)=(/ 0.042037,0.019687,0.008423,0.900,0.900,0.000383 /)
XEXT_COEFF_550_LKT(11,1)=13.410000 !rg=0.0231791 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,2,1:6)=(/ 183.480000,41.618000,6.066200,1.019100,0.612990,132.840000 /)
XPIZA_LKT(11,2,1:6)=(/ 0.796762,0.764728,0.953567,0.799812,0.900,0.000134 /)
XCGA_LKT(11,2,1:6)=(/ 0.053293,0.024983,0.010700,0.004027,0.900,0.000487 /)
XEXT_COEFF_550_LKT(11,2)=16.719000 !rg=0.0231791 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,3,1:6)=(/ 255.430000,57.747000,8.989000,1.431100,0.656310,132.900000 /)
XPIZA_LKT(11,3,1:6)=(/ 0.850259,0.828038,0.968460,0.857088,0.900,0.000201 /)
XCGA_LKT(11,3,1:6)=(/ 0.080740,0.037973,0.016303,0.006147,0.900,0.000740 /)
XEXT_COEFF_550_LKT(11,3)=24.784000 !rg=0.0231791 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,4,1:6)=(/ 393.410000,91.208000,15.170000,2.304100,0.747920,133.010000 /)
XPIZA_LKT(11,4,1:6)=(/ 0.898993,0.888585,0.981102,0.910845,0.900,0.000344 /)
XCGA_LKT(11,4,1:6)=(/ 0.135413,0.064380,0.027787,0.010503,0.900,0.001270 /)
XEXT_COEFF_550_LKT(11,4)=41.732000 !rg=0.0231791 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,5,1:6)=(/ 621.260000,154.030000,27.333000,4.033600,0.929220,133.190000 /)
XPIZA_LKT(11,5,1:6)=(/ 0.932651,0.931784,0.989330,0.948715,0.900,0.000627 /)
XCGA_LKT(11,5,1:6)=(/ 0.221830,0.107297,0.046533,0.017647,0.900,0.002140 /)
XEXT_COEFF_550_LKT(11,5)=74.500000 !rg=0.0231791 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,6,1:6)=(/ 957.600000,258.090000,49.355000,7.217200,1.263000,133.470000 /)
XPIZA_LKT(11,6,1:6)=(/ 0.953323,0.957483,0.993944,0.971039,0.578536,0.001148 /)
XCGA_LKT(11,6,1:6)=(/ 0.327920,0.167507,0.072713,0.027643,0.008987,0.003360 /)
XEXT_COEFF_550_LKT(11,6)=131.630000 !rg=0.0231791 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,7,1:6)=(/ 1419.700000,416.830000,87.600000,12.956000,1.866300,133.920000 /)
XPIZA_LKT(11,7,1:6)=(/ 0.966025,0.972185,0.996473,0.983623,0.713300,0.002086 /)
XCGA_LKT(11,7,1:6)=(/ 0.428487,0.250387,0.108947,0.041433,0.013500,0.005053 /)
XEXT_COEFF_550_LKT(11,7)=224.260000 !rg=0.0231791 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,8,1:6)=(/ 1992.800000,647.210000,149.840000,23.045000,2.936200,134.610000 /)
XPIZA_LKT(11,8,1:6)=(/ 0.973772,0.980785,0.997850,0.990598,0.816395,0.003740 /)
XCGA_LKT(11,8,1:6)=(/ 0.514853,0.353013,0.160080,0.060707,0.019820,0.007427 /)
XEXT_COEFF_550_LKT(11,8)=361.970000 !rg=0.0231791 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,9,1:6)=(/ 2648.900000,968.300000,244.930000,40.662000,4.845100,135.700000 /)
XPIZA_LKT(11,9,1:6)=(/ 0.978691,0.986100,0.998614,0.994517,0.887523,0.006670 /)
XCGA_LKT(11,9,1:6)=(/ 0.584813,0.447113,0.233953,0.088233,0.028830,0.010823 /)
XEXT_COEFF_550_LKT(11,9)=563.520000 !rg=0.0231791 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,10,1:6)=(/ 3311.200000,1358.100000,381.990000,70.264000,8.216900,137.410000 /)
XPIZA_LKT(11,10,1:6)=(/ 0.981713,0.989276,0.999050,0.996707,0.932653,0.011786 /)
XCGA_LKT(11,10,1:6)=(/ 0.638337,0.524920,0.334963,0.128023,0.041760,0.015707 /)
XEXT_COEFF_550_LKT(11,10)=842.030000 !rg=0.0231791 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,11,1:6)=(/ 3886.700000,1811.700000,582.320000,117.700000,14.205000,140.150000 /)
XPIZA_LKT(11,11,1:6)=(/ 0.983424,0.991302,0.999323,0.997940,0.960192,0.020723 /)
XCGA_LKT(11,11,1:6)=(/ 0.677397,0.592767,0.438367,0.186650,0.060473,0.022800 /)
XEXT_COEFF_550_LKT(11,11)=1179.400000 !rg=0.0231791 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,12,1:6)=(/ 4255.200000,2268.800000,833.220000,187.590000,24.666000,144.550000 /)
XPIZA_LKT(11,12,1:6)=(/ 0.984052,0.992536,0.999488,0.998631,0.976389,0.036025 /)
XCGA_LKT(11,12,1:6)=(/ 0.702700,0.643750,0.514980,0.272873,0.087663,0.033107 /)
XEXT_COEFF_550_LKT(11,12)=1566.500000 !rg=0.0231791 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,13,1:6)=(/ 4329.500000,2660.500000,1126.600000,288.220000,42.392000,151.670000 /)
XPIZA_LKT(11,13,1:6)=(/ 0.983631,0.993218,0.999589,0.999038,0.985724,0.061471 /)
XCGA_LKT(11,13,1:6)=(/ 0.714917,0.681077,0.585463,0.385507,0.127477,0.048133 /)
XEXT_COEFF_550_LKT(11,13)=1950.600000 !rg=0.0231791 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,14,1:6)=(/ 4110.500000,2904.000000,1425.200000,430.430000,70.514000,163.110000 /)
XPIZA_LKT(11,14,1:6)=(/ 0.982173,0.993449,0.999651,0.999298,0.991001,0.101622 /)
XCGA_LKT(11,14,1:6)=(/ 0.716590,0.704577,0.637560,0.481067,0.186387,0.070127 /)
XEXT_COEFF_550_LKT(11,14)=2269.000000 !rg=0.0231791 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,15,1:6)=(/ 3727.600000,2941.500000,1692.300000,594.970000,111.460000,181.260000 /)
XPIZA_LKT(11,15,1:6)=(/ 0.979826,0.993236,0.999687,0.999452,0.993968,0.160376 /)
XCGA_LKT(11,15,1:6)=(/ 0.714623,0.714740,0.677030,0.550747,0.273717,0.102660 /)
XEXT_COEFF_550_LKT(11,15)=2456.500000 !rg=0.0231791 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,16,1:6)=(/ 3323.700000,2780.400000,1876.000000,781.690000,169.710000,209.100000 /)
XPIZA_LKT(11,16,1:6)=(/ 0.976835,0.992578,0.999701,0.999549,0.995719,0.237019 /)
XCGA_LKT(11,16,1:6)=(/ 0.715633,0.714530,0.702433,0.614700,0.389710,0.151650 /)
XEXT_COEFF_550_LKT(11,16)=2468.300000 !rg=0.0231791 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,17,1:6)=(/ 2900.400000,2512.900000,1927.900000,959.100000,252.400000,248.570000 /)
XPIZA_LKT(11,17,1:6)=(/ 0.973409,0.991499,0.999694,0.999607,0.996859,0.321416 /)
XCGA_LKT(11,17,1:6)=(/ 0.720967,0.709577,0.714377,0.660093,0.488130,0.226123 /)
XEXT_COEFF_550_LKT(11,17)=2312.800000 !rg=0.0231791 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,18,1:6)=(/ 2462.500000,2234.500000,1840.800000,1101.000000,346.470000,301.190000 /)
XPIZA_LKT(11,18,1:6)=(/ 0.969297,0.990306,0.999667,0.999636,0.997536,0.401010 /)
XCGA_LKT(11,18,1:6)=(/ 0.724100,0.711007,0.714350,0.692470,0.555547,0.330973 /)
XEXT_COEFF_550_LKT(11,18)=2081.000000 !rg=0.0231791 sigma=2.75 wvl=0.55
 
END SUBROUTINE SALT_OPT_LKT_SET1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET2()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE

 
XEXT_COEFF_WVL_LKT(11,19,1:6)=(/ 2089.700000,1963.500000,1667.800000,1177.500000,452.220000,365.670000 /)
XPIZA_LKT(11,19,1:6)=(/ 0.963428,0.988751,0.999619,0.999642,0.997962,0.472113 /)
XCGA_LKT(11,19,1:6)=(/ 0.728703,0.715453,0.707507,0.711427,0.619120,0.433767 /)
XEXT_COEFF_550_LKT(11,19)=1854.800000 !rg=0.0231791 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,20,1:6)=(/ 1736.900000,1678.300000,1483.400000,1167.800000,551.750000,427.980000 /)
XPIZA_LKT(11,20,1:6)=(/ 0.957422,0.986387,0.999562,0.999623,0.998214,0.528858 /)
XCGA_LKT(11,20,1:6)=(/ 0.732047,0.717120,0.706287,0.716947,0.663407,0.506160 /)
XEXT_COEFF_550_LKT(11,20)=1623.500000 !rg=0.0231791 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,1,1:6)=(/ 185.070000,41.911000,6.116200,1.026200,0.613760,132.840000 /)
XPIZA_LKT(12,1,1:6)=(/ 0.797943,0.766062,0.953920,0.801149,0.900,0.000135 /)
XCGA_LKT(12,1,1:6)=(/ 0.049273,0.023083,0.009883,0.003720,0.900,0.000447 /)
XEXT_COEFF_550_LKT(12,1)=16.858000 !rg=0.0251112 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,2,1:6)=(/ 223.370000,50.345000,7.639000,1.240800,0.636340,132.870000 /)
XPIZA_LKT(12,2,1:6)=(/ 0.830162,0.803723,0.962973,0.835318,0.900,0.000170 /)
XCGA_LKT(12,2,1:6)=(/ 0.062460,0.029287,0.012553,0.004727,0.900,0.000570 /)
XEXT_COEFF_550_LKT(12,2)=21.065000 !rg=0.0251112 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,3,1:6)=(/ 312.910000,70.773000,11.355000,1.764700,0.691370,132.950000 /)
XPIZA_LKT(12,3,1:6)=(/ 0.875314,0.858110,0.974901,0.883869,0.900,0.000256 /)
XCGA_LKT(12,3,1:6)=(/ 0.094550,0.044487,0.019117,0.007210,0.900,0.000870 /)
XEXT_COEFF_550_LKT(12,3)=31.308000 !rg=0.0251112 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,4,1:6)=(/ 479.690000,112.790000,19.203000,2.874800,0.807800,133.070000 /)
XPIZA_LKT(12,4,1:6)=(/ 0.915224,0.908618,0.984968,0.928345,0.900,0.000438 /)
XCGA_LKT(12,4,1:6)=(/ 0.157900,0.075307,0.032557,0.012317,0.900,0.001490 /)
XEXT_COEFF_550_LKT(12,4)=52.742000 !rg=0.0251112 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,5,1:6)=(/ 747.430000,190.130000,34.592000,5.073300,1.038200,133.290000 /)
XPIZA_LKT(12,5,1:6)=(/ 0.942456,0.943747,0.991489,0.959065,0.488411,0.000797 /)
XCGA_LKT(12,5,1:6)=(/ 0.254997,0.125333,0.054483,0.020690,0.006717,0.002510 /)
XEXT_COEFF_550_LKT(12,5)=93.743000 !rg=0.0251112 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,6,1:6)=(/ 1133.300000,314.650000,62.230000,9.117100,1.462500,133.630000 /)
XPIZA_LKT(12,6,1:6)=(/ 0.959293,0.964355,0.995136,0.976946,0.635336,0.001458 /)
XCGA_LKT(12,6,1:6)=(/ 0.364107,0.195427,0.085113,0.032393,0.010543,0.003943 /)
XEXT_COEFF_550_LKT(12,6)=163.770000 !rg=0.0251112 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,7,1:6)=(/ 1643.800000,501.060000,109.440000,16.392000,2.229400,134.170000 /)
XPIZA_LKT(12,7,1:6)=(/ 0.969649,0.976210,0.997131,0.986954,0.759328,0.002648 /)
XCGA_LKT(12,7,1:6)=(/ 0.461730,0.289107,0.127587,0.048540,0.015833,0.005927 /)
XEXT_COEFF_550_LKT(12,7)=274.060000 !rg=0.0251112 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,8,1:6)=(/ 2260.300000,768.790000,184.150000,29.122000,3.589100,135.000000 /)
XPIZA_LKT(12,8,1:6)=(/ 0.976083,0.983285,0.998215,0.992480,0.849192,0.004746 /)
XCGA_LKT(12,8,1:6)=(/ 0.544023,0.391917,0.187657,0.071110,0.023233,0.008713 /)
XEXT_COEFF_550_LKT(12,8)=435.310000 !rg=0.0251112 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,9,1:6)=(/ 2929.500000,1120.000000,295.220000,51.120000,6.013500,136.320000 /)
XPIZA_LKT(12,9,1:6)=(/ 0.980111,0.987593,0.998820,0.995576,0.908856,0.008452 /)
XCGA_LKT(12,9,1:6)=(/ 0.607223,0.478110,0.273363,0.103407,0.033790,0.012697 /)
XEXT_COEFF_550_LKT(12,9)=670.100000 !rg=0.0251112 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,10,1:6)=(/ 3566.200000,1540.600000,456.060000,87.306000,10.288000,138.390000 /)
XPIZA_LKT(12,10,1:6)=(/ 0.982545,0.990216,0.999177,0.997301,0.945776,0.014898 /)
XCGA_LKT(12,10,1:6)=(/ 0.655170,0.554273,0.380083,0.150257,0.048927,0.018423 /)
XEXT_COEFF_550_LKT(12,10)=972.570000 !rg=0.0251112 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,11,1:6)=(/ 4071.600000,1998.800000,680.630000,143.380000,17.846000,141.720000 /)
XPIZA_LKT(12,11,1:6)=(/ 0.983804,0.991873,0.999402,0.998271,0.967961,0.026094 /)
XCGA_LKT(12,11,1:6)=(/ 0.688723,0.614023,0.471100,0.219430,0.070860,0.026740 /)
XEXT_COEFF_550_LKT(12,11)=1335.700000 !rg=0.0251112 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,12,1:6)=(/ 4329.200000,2439.200000,947.350000,224.230000,30.918000,147.080000 /)
XPIZA_LKT(12,12,1:6)=(/ 0.984005,0.992873,0.999534,0.998821,0.980885,0.045069 /)
XCGA_LKT(12,12,1:6)=(/ 0.708760,0.660117,0.544420,0.318530,0.102787,0.038837 /)
XEXT_COEFF_550_LKT(12,12)=1724.000000 !rg=0.0251112 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,13,1:6)=(/ 4277.900000,2781.300000,1249.900000,342.450000,52.602000,155.770000 /)
XPIZA_LKT(12,13,1:6)=(/ 0.983171,0.993369,0.999619,0.999160,0.988280,0.076094 /)
XCGA_LKT(12,13,1:6)=(/ 0.716573,0.691743,0.608313,0.430110,0.149740,0.056497 /)
XEXT_COEFF_550_LKT(12,13)=2090.600000 !rg=0.0251112 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,14,1:6)=(/ 3975.100000,2946.600000,1543.500000,495.870000,85.793000,169.670000 /)
XPIZA_LKT(12,14,1:6)=(/ 0.981314,0.993414,0.999669,0.999372,0.992437,0.123703 /)
XCGA_LKT(12,14,1:6)=(/ 0.716073,0.709787,0.655373,0.509940,0.219477,0.082430 /)
XEXT_COEFF_550_LKT(12,14)=2364.000000 !rg=0.0251112 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,15,1:6)=(/ 3572.200000,2897.900000,1781.700000,670.860000,132.880000,191.470000 /)
XPIZA_LKT(12,15,1:6)=(/ 0.978653,0.993023,0.999695,0.999497,0.994791,0.190388 /)
XCGA_LKT(12,15,1:6)=(/ 0.714293,0.715563,0.688787,0.579597,0.320600,0.121047 /)
XEXT_COEFF_550_LKT(12,15)=2483.700000 !rg=0.0251112 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,16,1:6)=(/ 3156.700000,2680.000000,1914.500000,855.520000,201.420000,224.040000 /)
XPIZA_LKT(12,16,1:6)=(/ 0.975642,0.992170,0.999701,0.999576,0.996254,0.272053 /)
XCGA_LKT(12,16,1:6)=(/ 0.718870,0.712137,0.708763,0.634223,0.436677,0.179653 /)
XEXT_COEFF_550_LKT(12,16)=2422.000000 !rg=0.0251112 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,17,1:6)=(/ 2721.100000,2397.200000,1906.400000,1022.900000,290.730000,268.840000 /)
XPIZA_LKT(12,17,1:6)=(/ 0.971827,0.991045,0.999686,0.999622,0.997189,0.355691 /)
XCGA_LKT(12,17,1:6)=(/ 0.722750,0.709857,0.715340,0.675067,0.516820,0.267670 /)
XEXT_COEFF_550_LKT(12,17)=2225.500000 !rg=0.0251112 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,18,1:6)=(/ 2307.600000,2129.800000,1775.000000,1140.900000,389.870000,327.280000 /)
XPIZA_LKT(12,18,1:6)=(/ 0.966671,0.989729,0.999650,0.999641,0.997734,0.431830 /)
XCGA_LKT(12,18,1:6)=(/ 0.725177,0.713793,0.711710,0.701513,0.584430,0.378180 /)
XEXT_COEFF_550_LKT(12,18)=1987.600000 !rg=0.0251112 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,19,1:6)=(/ 1936.400000,1847.200000,1589.800000,1183.400000,494.070000,392.470000 /)
XPIZA_LKT(12,19,1:6)=(/ 0.960896,0.987884,0.999597,0.999637,0.998083,0.498108 /)
XCGA_LKT(12,19,1:6)=(/ 0.728990,0.716090,0.706760,0.714787,0.638540,0.466367 /)
XEXT_COEFF_550_LKT(12,19)=1761.800000 !rg=0.0251112 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,20,1:6)=(/ 1604.100000,1564.300000,1417.500000,1139.000000,587.900000,452.510000 /)
XPIZA_LKT(12,20,1:6)=(/ 0.955590,0.985484,0.999536,0.999607,0.998284,0.547369 /)
XCGA_LKT(12,20,1:6)=(/ 0.737770,0.717270,0.708173,0.715890,0.678493,0.534780 /)
XEXT_COEFF_550_LKT(12,20)=1524.100000 !rg=0.0251112 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,1,1:6)=(/ 225.520000,50.720000,7.702600,1.249800,0.637320,132.880000 /)
XPIZA_LKT(13,1,1:6)=(/ 0.831226,0.804876,0.963254,0.836462,0.900,0.000172 /)
XCGA_LKT(13,1,1:6)=(/ 0.057760,0.027063,0.011593,0.004363,0.900,0.000527 /)
XEXT_COEFF_550_LKT(13,1)=21.242000 !rg=0.0272044 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,2,1:6)=(/ 273.590000,61.421000,9.639000,1.522700,0.666000,132.920000 /)
XPIZA_LKT(13,2,1:6)=(/ 0.858633,0.837415,0.970514,0.865555,0.900,0.000216 /)
XCGA_LKT(13,2,1:6)=(/ 0.073217,0.034330,0.014720,0.005547,0.900,0.000670 /)
XEXT_COEFF_550_LKT(13,2)=26.589000 !rg=0.0272044 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,3,1:6)=(/ 384.120000,87.253000,14.362000,2.188900,0.735940,133.000000 /)
XPIZA_LKT(13,3,1:6)=(/ 0.896183,0.883431,0.980036,0.906152,0.900,0.000326 /)
XCGA_LKT(13,3,1:6)=(/ 0.110733,0.052113,0.022413,0.008457,0.900,0.001020 /)
XEXT_COEFF_550_LKT(13,3)=39.588000 !rg=0.0272044 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,4,1:6)=(/ 583.330000,139.760000,24.320000,3.600300,0.883900,133.160000 /)
XPIZA_LKT(13,4,1:6)=(/ 0.928514,0.925067,0.988035,0.942601,0.900,0.000556 /)
XCGA_LKT(13,4,1:6)=(/ 0.183780,0.088050,0.038137,0.014443,0.900,0.001750 /)
XEXT_COEFF_550_LKT(13,4)=66.630000 !rg=0.0272044 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,5,1:6)=(/ 894.830000,234.070000,43.746000,6.394600,1.176800,133.410000 /)
XPIZA_LKT(13,5,1:6)=(/ 0.950490,0.953402,0.993196,0.967375,0.547956,0.001013 /)
XCGA_LKT(13,5,1:6)=(/ 0.290133,0.146273,0.063773,0.024250,0.007880,0.002947 /)
XEXT_COEFF_550_LKT(13,5)=117.610000 !rg=0.0272044 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,6,1:6)=(/ 1330.400000,381.710000,78.242000,11.528000,1.716100,133.820000 /)
XPIZA_LKT(13,6,1:6)=(/ 0.964166,0.969898,0.996076,0.981649,0.688533,0.001852 /)
XCGA_LKT(13,6,1:6)=(/ 0.399533,0.227233,0.099607,0.037953,0.012363,0.004627 /)
XEXT_COEFF_550_LKT(13,6)=202.520000 !rg=0.0272044 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,7,1:6)=(/ 1890.400000,600.240000,135.920000,20.736000,2.690900,134.470000 /)
XPIZA_LKT(13,7,1:6)=(/ 0.972668,0.979529,0.997647,0.989594,0.799958,0.003362 /)
XCGA_LKT(13,7,1:6)=(/ 0.494523,0.328950,0.149410,0.056857,0.018563,0.006957 /)
XEXT_COEFF_550_LKT(13,7)=332.650000 !rg=0.0272044 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,8,1:6)=(/ 2535.700000,904.000000,224.460000,36.741000,4.418500,135.470000 /)
XPIZA_LKT(13,8,1:6)=(/ 0.977979,0.985318,0.998502,0.993966,0.876932,0.006018 /)
XCGA_LKT(13,8,1:6)=(/ 0.569940,0.426627,0.219780,0.083303,0.027233,0.010220 /)
XEXT_COEFF_550_LKT(13,8)=521.820000 !rg=0.0272044 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,9,1:6)=(/ 3208.900000,1284.400000,354.150000,64.003000,7.495800,137.060000 /)
XPIZA_LKT(13,9,1:6)=(/ 0.981297,0.988810,0.998987,0.996410,0.926395,0.010699 /)
XCGA_LKT(13,9,1:6)=(/ 0.627850,0.508920,0.316647,0.121237,0.039597,0.014893 /)
XEXT_COEFF_550_LKT(13,9)=788.040000 !rg=0.0272044 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,10,1:6)=(/ 3804.100000,1728.000000,541.720000,107.650000,12.907000,139.570000 /)
XPIZA_LKT(13,10,1:6)=(/ 0.983198,0.990996,0.999283,0.997768,0.956378,0.018807 /)
XCGA_LKT(13,10,1:6)=(/ 0.669937,0.579977,0.420470,0.176463,0.057327,0.021610 /)
XEXT_COEFF_550_LKT(13,10)=1113.000000 !rg=0.0272044 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,11,1:6)=(/ 4219.100000,2187.800000,784.000000,173.060000,22.415000,143.630000 /)
XPIZA_LKT(13,11,1:6)=(/ 0.984020,0.992347,0.999464,0.998533,0.974170,0.032782 /)
XCGA_LKT(13,11,1:6)=(/ 0.698210,0.633920,0.500563,0.257670,0.083040,0.031363 /)
XEXT_COEFF_550_LKT(13,11)=1494.700000 !rg=0.0272044 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,12,1:6)=(/ 4352.200000,2597.200000,1070.200000,267.310000,38.628000,150.170000 /)
XPIZA_LKT(13,12,1:6)=(/ 0.983805,0.993120,0.999574,0.998978,0.984450,0.056171 /)
XCGA_LKT(13,12,1:6)=(/ 0.713147,0.674153,0.572933,0.366703,0.120587,0.045563 /)
XEXT_COEFF_550_LKT(13,12)=1882.100000 !rg=0.0272044 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,13,1:6)=(/ 4186.500000,2876.600000,1371.600000,403.130000,64.792000,160.760000 /)
XPIZA_LKT(13,13,1:6)=(/ 0.982548,0.993433,0.999642,0.999261,0.990294,0.093619 /)
XCGA_LKT(13,13,1:6)=(/ 0.717137,0.700517,0.628250,0.466710,0.176067,0.066340 /)
XEXT_COEFF_550_LKT(13,13)=2219.500000 !rg=0.0272044 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,14,1:6)=(/ 3824.200000,2953.900000,1649.000000,564.040000,103.370000,177.590000 /)
XPIZA_LKT(13,14,1:6)=(/ 0.980311,0.993310,0.999682,0.999430,0.993567,0.149168 /)
XCGA_LKT(13,14,1:6)=(/ 0.714960,0.713373,0.670353,0.538207,0.258353,0.096997 /)
XEXT_COEFF_550_LKT(13,14)=2436.400000 !rg=0.0272044 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,15,1:6)=(/ 3410.400000,2826.800000,1850.100000,748.240000,157.950000,203.520000 /)
XPIZA_LKT(13,15,1:6)=(/ 0.977587,0.992727,0.999700,0.999535,0.995468,0.223041 /)
XCGA_LKT(13,15,1:6)=(/ 0.716607,0.714687,0.698080,0.604750,0.370850,0.142993 /)
XEXT_COEFF_550_LKT(13,15)=2480.200000 !rg=0.0272044 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,16,1:6)=(/ 2989.300000,2565.900000,1928.800000,929.390000,237.230000,241.110000 /)
XPIZA_LKT(13,16,1:6)=(/ 0.974262,0.991777,0.999697,0.999598,0.996704,0.307406 /)
XCGA_LKT(13,16,1:6)=(/ 0.721143,0.711410,0.712570,0.652523,0.474890,0.213097 /)
XEXT_COEFF_550_LKT(13,16)=2355.100000 !rg=0.0272044 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,17,1:6)=(/ 2542.800000,2292.200000,1864.700000,1079.900000,329.990000,291.580000 /)
XPIZA_LKT(13,17,1:6)=(/ 0.969855,0.990542,0.999673,0.999632,0.997447,0.388483 /)
XCGA_LKT(13,17,1:6)=(/ 0.723103,0.711437,0.714657,0.687450,0.544293,0.314030 /)
XEXT_COEFF_550_LKT(13,17)=2128.100000 !rg=0.0272044 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,18,1:6)=(/ 2149.300000,2021.100000,1701.600000,1169.300000,434.550000,354.640000 /)
XPIZA_LKT(13,18,1:6)=(/ 0.964439,0.988951,0.999630,0.999642,0.997903,0.461125 /)
XCGA_LKT(13,18,1:6)=(/ 0.726783,0.714407,0.709777,0.708577,0.610133,0.419883 /)
XEXT_COEFF_550_LKT(13,18)=1896.900000 !rg=0.0272044 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,19,1:6)=(/ 1791.200000,1723.000000,1520.500000,1174.200000,535.730000,417.990000 /)
XPIZA_LKT(13,19,1:6)=(/ 0.958456,0.986974,0.999570,0.999628,0.998179,0.520806 /)
XCGA_LKT(13,19,1:6)=(/ 0.733497,0.716420,0.705853,0.716457,0.656547,0.495327 /)
XEXT_COEFF_550_LKT(13,19)=1670.000000 !rg=0.0272044 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,20,1:6)=(/ 1484.400000,1457.500000,1346.500000,1100.500000,619.460000,475.420000 /)
XPIZA_LKT(13,20,1:6)=(/ 0.952222,0.985070,0.999504,0.999586,0.998328,0.563615 /)
XCGA_LKT(13,20,1:6)=(/ 0.739847,0.723867,0.710077,0.712933,0.690747,0.561990 /)
XEXT_COEFF_550_LKT(13,20)=1420.500000 !rg=0.0272044 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,1,1:6)=(/ 276.540000,61.905000,9.719900,1.534100,0.667240,132.920000 /)
XPIZA_LKT(14,1,1:6)=(/ 0.859608,0.838397,0.970736,0.866516,0.900,0.000218 /)
XCGA_LKT(14,1,1:6)=(/ 0.067717,0.031727,0.013597,0.005120,0.900,0.000617 /)
XEXT_COEFF_550_LKT(14,1)=26.817000 !rg=0.0294721 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,2,1:6)=(/ 336.540000,75.469000,12.182000,1.881200,0.703690,132.970000 /)
XPIZA_LKT(14,2,1:6)=(/ 0.882586,0.866071,0.976537,0.890935,0.900,0.000275 /)
XCGA_LKT(14,2,1:6)=(/ 0.085853,0.040233,0.017263,0.006507,0.900,0.000783 /)
XEXT_COEFF_550_LKT(14,2)=33.609000 !rg=0.0294721 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,3,1:6)=(/ 471.560000,108.060000,18.184000,2.728400,0.792560,133.070000 /)
XPIZA_LKT(14,3,1:6)=(/ 0.913391,0.904501,0.984121,0.924498,0.900,0.000414 /)
XCGA_LKT(14,3,1:6)=(/ 0.129713,0.061037,0.026273,0.009920,0.900,0.001200 /)
XEXT_COEFF_550_LKT(14,3)=50.089000 !rg=0.0294721 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,4,1:6)=(/ 706.530000,173.250000,30.804000,4.522800,0.980620,133.250000 /)
XPIZA_LKT(14,4,1:6)=(/ 0.939355,0.938471,0.990466,0.954136,0.900,0.000707 /)
XCGA_LKT(14,4,1:6)=(/ 0.213230,0.102897,0.044663,0.016937,0.900,0.002053 /)
XEXT_COEFF_550_LKT(14,4)=84.089000 !rg=0.0294721 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,5,1:6)=(/ 1064.600000,287.080000,55.244000,8.072900,1.352900,133.560000 /)
XPIZA_LKT(14,5,1:6)=(/ 0.957055,0.961176,0.994546,0.974021,0.606097,0.001287 /)
XCGA_LKT(14,5,1:6)=(/ 0.326547,0.170470,0.074620,0.028420,0.009243,0.003457 /)
XEXT_COEFF_550_LKT(14,5)=146.980000 !rg=0.0294721 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,6,1:6)=(/ 1550.200000,461.060000,97.989000,14.583000,2.038500,134.050000 /)
XPIZA_LKT(14,6,1:6)=(/ 0.968170,0.974403,0.996817,0.985386,0.737106,0.002352 /)
XCGA_LKT(14,6,1:6)=(/ 0.434803,0.262397,0.116537,0.044460,0.014500,0.005427 /)
XEXT_COEFF_550_LKT(14,6)=248.790000 !rg=0.0294721 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,7,1:6)=(/ 2155.700000,715.050000,167.630000,26.212000,3.277400,134.830000 /)
XPIZA_LKT(14,7,1:6)=(/ 0.975187,0.982260,0.998055,0.991683,0.835136,0.004266 /)
XCGA_LKT(14,7,1:6)=(/ 0.525240,0.367143,0.174907,0.066590,0.021763,0.008160 /)
XEXT_COEFF_550_LKT(14,7)=401.690000 !rg=0.0294721 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,8,1:6)=(/ 2819.800000,1050.800000,271.590000,46.234000,5.471800,136.040000 /)
XPIZA_LKT(14,8,1:6)=(/ 0.979541,0.986953,0.998731,0.995139,0.900085,0.007628 /)
XCGA_LKT(14,8,1:6)=(/ 0.593847,0.458797,0.256617,0.097593,0.031917,0.011990 /)
XEXT_COEFF_550_LKT(14,8)=622.220000 !rg=0.0294721 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,9,1:6)=(/ 3475.300000,1463.200000,423.700000,79.689000,9.374300,137.960000 /)
XPIZA_LKT(14,9,1:6)=(/ 0.982234,0.989830,0.999126,0.997066,0.940694,0.013532 /)
XCGA_LKT(14,9,1:6)=(/ 0.645983,0.539103,0.360720,0.142200,0.046393,0.017470 /)
XEXT_COEFF_550_LKT(14,9)=914.640000 !rg=0.0294721 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,10,1:6)=(/ 4008.700000,1915.900000,636.050000,131.580000,16.210000,141.020000 /)
XPIZA_LKT(14,10,1:6)=(/ 0.983673,0.991624,0.999369,0.998135,0.964899,0.023701 /)
XCGA_LKT(14,10,1:6)=(/ 0.682457,0.602347,0.454627,0.207293,0.067167,0.025347 /)
XEXT_COEFF_550_LKT(14,10)=1265.700000 !rg=0.0294721 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,11,1:6)=(/ 4318.500000,2365.600000,894.630000,207.420000,28.113000,145.950000 /)
XPIZA_LKT(14,11,1:6)=(/ 0.984068,0.992728,0.999514,0.998741,0.979115,0.041067 /)
XCGA_LKT(14,11,1:6)=(/ 0.705400,0.651447,0.530173,0.301170,0.097343,0.036790 /)
XEXT_COEFF_550_LKT(14,11)=1652.600000 !rg=0.0294721 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,12,1:6)=(/ 4326.300000,2731.700000,1194.000000,318.160000,48.025000,153.940000 /)
XPIZA_LKT(14,12,1:6)=(/ 0.983444,0.993308,0.999606,0.999110,0.987269,0.069682 /)
XCGA_LKT(14,12,1:6)=(/ 0.715897,0.685940,0.597260,0.412360,0.141573,0.053473 /)
XEXT_COEFF_550_LKT(14,12)=2028.600000 !rg=0.0294721 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,13,1:6)=(/ 4062.700000,2935.700000,1493.000000,467.100000,79.098000,166.820000 /)
XPIZA_LKT(14,13,1:6)=(/ 0.981774,0.993438,0.999661,0.999342,0.991876,0.114331 /)
XCGA_LKT(14,13,1:6)=(/ 0.716643,0.706823,0.646930,0.496903,0.207203,0.077957 /)
XEXT_COEFF_550_LKT(14,13)=2325.900000 !rg=0.0294721 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,14,1:6)=(/ 3663.900000,2928.100000,1745.300000,637.810000,123.550000,187.090000 /)
XPIZA_LKT(14,14,1:6)=(/ 0.979336,0.993125,0.999691,0.999478,0.994466,0.177900 /)
XCGA_LKT(14,14,1:6)=(/ 0.715913,0.714897,0.683053,0.567407,0.303090,0.114300 /)
XEXT_COEFF_550_LKT(14,14)=2477.100000 !rg=0.0294721 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,15,1:6)=(/ 3254.100000,2732.400000,1899.900000,822.610000,187.660000,217.570000 /)
XPIZA_LKT(14,15,1:6)=(/ 0.976346,0.992395,0.999701,0.999565,0.996042,0.257500 /)
XCGA_LKT(14,15,1:6)=(/ 0.718837,0.713870,0.705500,0.625320,0.419263,0.169257 /)
XEXT_COEFF_550_LKT(14,15)=2448.700000 !rg=0.0294721 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,16,1:6)=(/ 2806.900000,2459.600000,1918.900000,996.760000,274.940000,260.420000 /)
XPIZA_LKT(14,16,1:6)=(/ 0.972554,0.991287,0.999690,0.999616,0.997065,0.342058 /)
XCGA_LKT(14,16,1:6)=(/ 0.722320,0.710753,0.714610,0.668607,0.505237,0.252473 /)
XEXT_COEFF_550_LKT(14,16)=2268.000000 !rg=0.0294721 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,17,1:6)=(/ 2371.500000,2184.700000,1806.700000,1125.500000,372.140000,316.770000 /)
XPIZA_LKT(14,17,1:6)=(/ 0.967862,0.989990,0.999658,0.999640,0.997658,0.419841 /)
XCGA_LKT(14,17,1:6)=(/ 0.724570,0.712790,0.713020,0.697500,0.573217,0.361567 /)
XEXT_COEFF_550_LKT(14,17)=2033.400000 !rg=0.0294721 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,18,1:6)=(/ 2003.400000,1898.900000,1628.800000,1181.100000,476.960000,381.680000 /)
XPIZA_LKT(14,18,1:6)=(/ 0.961639,0.988326,0.999606,0.999639,0.998036,0.488134 /)
XCGA_LKT(14,18,1:6)=(/ 0.730277,0.716603,0.706817,0.712903,0.630737,0.454327 /)
XEXT_COEFF_550_LKT(14,18)=1811.400000 !rg=0.0294721 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,19,1:6)=(/ 1654.600000,1612.700000,1447.200000,1151.700000,573.980000,442.770000 /)
XPIZA_LKT(14,19,1:6)=(/ 0.956115,0.985838,0.999549,0.999613,0.998258,0.540304 /)
XCGA_LKT(14,19,1:6)=(/ 0.735683,0.719647,0.708130,0.715853,0.672663,0.523907 /)
XEXT_COEFF_550_LKT(14,19)=1566.800000 !rg=0.0294721 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,20,1:6)=(/ 1366.500000,1356.900000,1273.000000,1053.300000,644.930000,494.670000 /)
XPIZA_LKT(14,20,1:6)=(/ 0.948517,0.984161,0.999477,0.999564,0.998360,0.577538 /)
XCGA_LKT(14,20,1:6)=(/ 0.741197,0.724960,0.712587,0.710427,0.700703,0.585160 /)
XEXT_COEFF_550_LKT(14,20)=1326.500000 !rg=0.0294721 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,1,1:6)=(/ 340.720000,76.100000,12.285000,1.895700,0.705260,132.980000 /)
XPIZA_LKT(15,1,1:6)=(/ 0.883503,0.866903,0.976713,0.891733,0.900,0.000277 /)
XCGA_LKT(15,1,1:6)=(/ 0.079420,0.037190,0.015947,0.006007,0.900,0.000723 /)
XEXT_COEFF_550_LKT(15,1)=33.902000 !rg=0.0319288 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,2,1:6)=(/ 415.040000,93.271000,15.415000,2.337000,0.751570,133.030000 /)
XPIZA_LKT(15,2,1:6)=(/ 0.902517,0.890129,0.981337,0.911982,0.900,0.000350 /)
XCGA_LKT(15,2,1:6)=(/ 0.100720,0.047153,0.020243,0.007633,0.900,0.000920 /)
XEXT_COEFF_550_LKT(15,2)=42.527000 !rg=0.0319288 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,3,1:6)=(/ 577.760000,134.250000,23.039000,3.414300,0.864520,133.150000 /)
XPIZA_LKT(15,3,1:6)=(/ 0.927465,0.921864,0.987366,0.939470,0.900,0.000526 /)
XCGA_LKT(15,3,1:6)=(/ 0.151993,0.071477,0.030790,0.011637,0.900,0.001407 /)
XEXT_COEFF_550_LKT(15,3)=63.386000 !rg=0.0319288 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,4,1:6)=(/ 851.290000,214.560000,39.006000,5.695400,1.103500,133.370000 /)
XPIZA_LKT(15,4,1:6)=(/ 0.948175,0.949327,0.992391,0.963419,0.518173,0.000898 /)
XCGA_LKT(15,4,1:6)=(/ 0.246260,0.120177,0.052290,0.019857,0.006447,0.002407 /)
XEXT_COEFF_550_LKT(15,4)=105.930000 !rg=0.0319288 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,5,1:6)=(/ 1257.800000,350.540000,69.616000,10.204000,1.576800,133.740000 /)
XPIZA_LKT(15,5,1:6)=(/ 0.962421,0.967442,0.995612,0.979320,0.661314,0.001635 /)
XCGA_LKT(15,5,1:6)=(/ 0.363960,0.198153,0.087277,0.033300,0.010843,0.004057 /)
XEXT_COEFF_550_LKT(15,5)=182.780000 !rg=0.0319288 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,6,1:6)=(/ 1793.800000,554.530000,122.100000,18.449000,2.448200,134.330000 /)
XPIZA_LKT(15,6,1:6)=(/ 0.971494,0.978088,0.997400,0.988349,0.780440,0.002986 /)
XCGA_LKT(15,6,1:6)=(/ 0.469603,0.299397,0.136293,0.052070,0.017000,0.006370 /)
XEXT_COEFF_550_LKT(15,6)=303.690000 !rg=0.0319288 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,7,1:6)=(/ 2433.200000,844.100000,205.180000,33.089000,4.022500,135.260000 /)
XPIZA_LKT(15,7,1:6)=(/ 0.977263,0.984479,0.998376,0.993334,0.865087,0.005411 /)
XCGA_LKT(15,7,1:6)=(/ 0.553130,0.402533,0.204540,0.077987,0.025510,0.009573 /)
XEXT_COEFF_550_LKT(15,7)=483.120000 !rg=0.0319288 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,8,1:6)=(/ 3106.600000,1210.900000,326.900000,57.966000,6.808300,136.730000 /)
XPIZA_LKT(15,8,1:6)=(/ 0.980847,0.988288,0.998916,0.996063,0.919197,0.009660 /)
XCGA_LKT(15,8,1:6)=(/ 0.615947,0.490557,0.297393,0.114360,0.037403,0.014067 /)
XEXT_COEFF_550_LKT(15,8)=734.560000 !rg=0.0319288 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,9,1:6)=(/ 3728.200000,1649.700000,504.580000,98.521000,11.751000,139.060000 /)
XPIZA_LKT(15,9,1:6)=(/ 0.982975,0.990679,0.999241,0.997581,0.952275,0.017092 /)
XCGA_LKT(15,9,1:6)=(/ 0.662050,0.566110,0.401407,0.166863,0.054357,0.020490 /)
XEXT_COEFF_550_LKT(15,9)=1051.400000 !rg=0.0319288 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,10,1:6)=(/ 4180.000000,2107.300000,736.310000,159.380000,20.360000,142.770000 /)
XPIZA_LKT(15,10,1:6)=(/ 0.983969,0.992143,0.999437,0.998423,0.971719,0.029805 /)
XCGA_LKT(15,10,1:6)=(/ 0.693170,0.623220,0.485033,0.243303,0.078700,0.029730 /)
XEXT_COEFF_550_LKT(15,10)=1423.900000 !rg=0.0319288 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,11,1:6)=(/ 4369.300000,2531.400000,1014.700000,247.690000,35.162000,148.790000 /)
XPIZA_LKT(15,11,1:6)=(/ 0.983958,0.993011,0.999557,0.998912,0.983040,0.051267 /)
XCGA_LKT(15,11,1:6)=(/ 0.710970,0.666537,0.559470,0.348027,0.114157,0.043163 /)
XEXT_COEFF_550_LKT(15,11)=1812.800000 !rg=0.0319288 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,12,1:6)=(/ 4257.300000,2841.000000,1316.000000,376.020000,59.309000,158.530000 /)
XPIZA_LKT(15,12,1:6)=(/ 0.982903,0.993406,0.999632,0.999220,0.989490,0.085952 /)
XCGA_LKT(15,12,1:6)=(/ 0.717087,0.695683,0.618080,0.451010,0.166363,0.062780 /)
XEXT_COEFF_550_LKT(15,12)=2164.100000 !rg=0.0319288 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,13,1:6)=(/ 3916.900000,2961.800000,1603.700000,533.550000,95.634000,174.150000 /)
XPIZA_LKT(15,13,1:6)=(/ 0.980955,0.993367,0.999676,0.999406,0.993120,0.138391 /)
XCGA_LKT(15,13,1:6)=(/ 0.717027,0.711510,0.663027,0.525277,0.243860,0.091693 /)
XEXT_COEFF_550_LKT(15,13)=2410.400000 !rg=0.0319288 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,14,1:6)=(/ 3511.500000,2869.000000,1822.300000,715.070000,147.030000,198.340000 /)
XPIZA_LKT(15,14,1:6)=(/ 0.978161,0.992883,0.999698,0.999520,0.995199,0.209492 /)
XCGA_LKT(15,14,1:6)=(/ 0.716890,0.715287,0.693360,0.594080,0.352137,0.134917 /)
XEXT_COEFF_550_LKT(15,14)=2488.700000 !rg=0.0319288 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,15,1:6)=(/ 3080.400000,2632.400000,1924.700000,896.920000,221.850000,233.700000 /)
XPIZA_LKT(15,15,1:6)=(/ 0.974940,0.991958,0.999699,0.999589,0.996527,0.292709 /)
XCGA_LKT(15,15,1:6)=(/ 0.720947,0.712053,0.710247,0.644143,0.460200,0.200647 /)
XEXT_COEFF_550_LKT(15,15)=2390.200000 !rg=0.0319288 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,16,1:6)=(/ 2623.100000,2348.500000,1886.800000,1056.500000,313.480000,282.120000 /)
XPIZA_LKT(15,16,1:6)=(/ 0.970881,0.990805,0.999679,0.999628,0.997347,0.375382 /)
XCGA_LKT(15,16,1:6)=(/ 0.723433,0.711337,0.714897,0.681867,0.532773,0.297177 /)
XEXT_COEFF_550_LKT(15,16)=2177.400000 !rg=0.0319288 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,17,1:6)=(/ 2222.000000,2074.300000,1737.800000,1158.900000,416.720000,343.680000 /)
XPIZA_LKT(15,17,1:6)=(/ 0.965145,0.989362,0.999638,0.999642,0.997839,0.449785 /)
XCGA_LKT(15,17,1:6)=(/ 0.727330,0.715550,0.710080,0.705307,0.600393,0.405100 /)
XEXT_COEFF_550_LKT(15,17)=1946.900000 !rg=0.0319288 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,18,1:6)=(/ 1847.100000,1782.400000,1552.400000,1179.000000,518.700000,407.520000 /)
XPIZA_LKT(15,18,1:6)=(/ 0.959590,0.987111,0.999584,0.999631,0.998141,0.512016 /)
XCGA_LKT(15,18,1:6)=(/ 0.731483,0.717130,0.707030,0.715563,0.649130,0.484007 /)
XEXT_COEFF_550_LKT(15,18)=1713.200000 !rg=0.0319288 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,19,1:6)=(/ 1527.500000,1501.600000,1379.700000,1116.100000,606.950000,466.500000 /)
XPIZA_LKT(15,19,1:6)=(/ 0.953466,0.985211,0.999522,0.999595,0.998310,0.557356 /)
XCGA_LKT(15,19,1:6)=(/ 0.738257,0.721560,0.710863,0.714130,0.685767,0.551860 /)
XEXT_COEFF_550_LKT(15,19)=1463.300000 !rg=0.0319288 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,20,1:6)=(/ 1264.900000,1256.100000,1197.700000,1007.300000,662.470000,510.100000 /)
XPIZA_LKT(15,20,1:6)=(/ 0.944595,0.982891,0.999432,0.999535,0.998367,0.588662 /)
XCGA_LKT(15,20,1:6)=(/ 0.745833,0.725190,0.712500,0.707160,0.708227,0.604897 /)
XEXT_COEFF_550_LKT(15,20)=1246.100000 !rg=0.0319288 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,1,1:6)=(/ 421.080000,94.103000,15.547000,2.355500,0.753560,133.040000 /)
XPIZA_LKT(16,1,1:6)=(/ 0.903408,0.890831,0.981476,0.912636,0.900,0.000353 /)
XCGA_LKT(16,1,1:6)=(/ 0.093200,0.043590,0.018700,0.007050,0.900,0.000850 /)
XEXT_COEFF_550_LKT(16,1)=42.906000 !rg=0.0345903 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,2,1:6)=(/ 512.130000,115.800000,19.526000,2.916600,0.812420,133.110000 /)
XPIZA_LKT(16,2,1:6)=(/ 0.918944,0.910106,0.985154,0.929263,0.900,0.000444 /)
XCGA_LKT(16,2,1:6)=(/ 0.118263,0.055260,0.023733,0.008953,0.900,0.001080 /)
XEXT_COEFF_550_LKT(16,2)=53.848000 !rg=0.0345903 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,3,1:6)=(/ 705.080000,167.090000,29.202000,4.286500,0.955960,133.250000 /)
XPIZA_LKT(16,3,1:6)=(/ 0.938899,0.936058,0.989939,0.951604,0.900,0.000668 /)
XCGA_LKT(16,3,1:6)=(/ 0.178140,0.083700,0.036080,0.013647,0.900,0.001650 /)
XEXT_COEFF_550_LKT(16,3)=80.193000 !rg=0.0345903 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,4,1:6)=(/ 1019.500000,265.060000,49.358000,7.185700,1.259800,133.510000 /)
XPIZA_LKT(16,4,1:6)=(/ 0.955345,0.958084,0.993914,0.970858,0.577206,0.001142 /)
XCGA_LKT(16,4,1:6)=(/ 0.282750,0.140257,0.061197,0.023273,0.007563,0.002827 /)
XEXT_COEFF_550_LKT(16,4)=133.110000 !rg=0.0345903 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,5,1:6)=(/ 1476.000000,426.050000,87.465000,12.906000,1.861400,133.950000 /)
XPIZA_LKT(16,5,1:6)=(/ 0.966829,0.972507,0.996453,0.983535,0.712389,0.002077 /)
XCGA_LKT(16,5,1:6)=(/ 0.402100,0.229190,0.102033,0.039007,0.012717,0.004760 /)
XEXT_COEFF_550_LKT(16,5)=226.040000 !rg=0.0345903 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,6,1:6)=(/ 2058.300000,663.010000,151.230000,23.329000,2.968900,134.660000 /)
XPIZA_LKT(16,6,1:6)=(/ 0.974254,0.981098,0.997860,0.990696,0.818308,0.003790 /)
XCGA_LKT(16,6,1:6)=(/ 0.502630,0.336340,0.159303,0.060977,0.019933,0.007473 /)
XEXT_COEFF_550_LKT(16,6)=368.700000 !rg=0.0345903 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,7,1:6)=(/ 2722.700000,986.290000,249.380000,41.677000,4.968800,135.780000 /)
XPIZA_LKT(16,7,1:6)=(/ 0.978981,0.986273,0.998631,0.994638,0.890226,0.006860 /)
XCGA_LKT(16,7,1:6)=(/ 0.579047,0.436103,0.238523,0.091330,0.029900,0.011230 /)
XEXT_COEFF_550_LKT(16,7)=577.810000 !rg=0.0345903 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,8,1:6)=(/ 3385.200000,1386.000000,392.160000,72.312000,8.502700,137.560000 /)
XPIZA_LKT(16,8,1:6)=(/ 0.981897,0.989402,0.999068,0.996791,0.934833,0.012223 /)
XCGA_LKT(16,8,1:6)=(/ 0.635627,0.521697,0.339677,0.134037,0.043823,0.016500 /)
XEXT_COEFF_550_LKT(16,8)=856.800000 !rg=0.0345903 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,9,1:6)=(/ 3952.600000,1838.600000,594.770000,120.800000,14.751000,140.390000 /)
XPIZA_LKT(16,9,1:6)=(/ 0.983542,0.991368,0.999335,0.997986,0.961600,0.021555 /)
XCGA_LKT(16,9,1:6)=(/ 0.675900,0.589783,0.436687,0.195840,0.063680,0.024033 /)
XEXT_COEFF_550_LKT(16,9)=1200.600000 !rg=0.0345903 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,10,1:6)=(/ 4305.600000,2291.600000,843.690000,191.590000,25.545000,144.920000 /)
XPIZA_LKT(16,10,1:6)=(/ 0.984110,0.992569,0.999492,0.998653,0.977158,0.037384 /)
XCGA_LKT(16,10,1:6)=(/ 0.701613,0.641943,0.515030,0.284513,0.092233,0.034873 /)
XEXT_COEFF_550_LKT(16,10)=1582.400000 !rg=0.0345903 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,11,1:6)=(/ 4371.700000,2679.000000,1138.200000,295.280000,43.790000,152.260000 /)
XPIZA_LKT(16,11,1:6)=(/ 0.983679,0.993233,0.999592,0.999055,0.986147,0.063726 /)
XCGA_LKT(16,11,1:6)=(/ 0.714733,0.679493,0.585207,0.393943,0.133963,0.050650 /)
XEXT_COEFF_550_LKT(16,11)=1965.200000 !rg=0.0345903 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,12,1:6)=(/ 4150.000000,2917.800000,1439.200000,438.090000,72.635000,164.120000 /)
XPIZA_LKT(16,12,1:6)=(/ 0.982266,0.993448,0.999653,0.999309,0.991237,0.105291 /)
XCGA_LKT(16,12,1:6)=(/ 0.717957,0.703247,0.637560,0.482773,0.195653,0.073750 /)
XEXT_COEFF_550_LKT(16,12)=2282.000000 !rg=0.0345903 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,13,1:6)=(/ 3769.300000,2951.700000,1705.300000,605.060000,114.630000,182.960000 /)
XPIZA_LKT(16,13,1:6)=(/ 0.979915,0.993224,0.999687,0.999458,0.994106,0.165776 /)
XCGA_LKT(16,13,1:6)=(/ 0.716880,0.713973,0.676663,0.554580,0.286343,0.107987 /)
XEXT_COEFF_550_LKT(16,13)=2465.400000 !rg=0.0345903 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,14,1:6)=(/ 3344.500000,2790.100000,1882.100000,790.220000,174.830000,211.520000 /)
XPIZA_LKT(16,14,1:6)=(/ 0.977035,0.992550,0.999700,0.999552,0.995816,0.243226 /)
XCGA_LKT(16,14,1:6)=(/ 0.719110,0.714243,0.701707,0.615973,0.401217,0.159560 /)
XEXT_COEFF_550_LKT(16,14)=2470.800000 !rg=0.0345903 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,15,1:6)=(/ 2900.500000,2520.000000,1926.000000,967.560000,258.660000,252.010000 /)
XPIZA_LKT(16,15,1:6)=(/ 0.973482,0.991538,0.999693,0.999608,0.996922,0.327602 /)
XCGA_LKT(16,15,1:6)=(/ 0.722720,0.711553,0.713380,0.661283,0.492637,0.237797 /)
XEXT_COEFF_550_LKT(16,15)=2316.300000 !rg=0.0345903 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,16,1:6)=(/ 2455.600000,2242.800000,1836.300000,1107.700000,354.350000,306.320000 /)
XPIZA_LKT(16,16,1:6)=(/ 0.968458,0.990253,0.999664,0.999637,0.997574,0.407283 /)
XCGA_LKT(16,16,1:6)=(/ 0.724337,0.713520,0.713420,0.692920,0.561447,0.344427 /)
XEXT_COEFF_550_LKT(16,16)=2085.700000 !rg=0.0345903 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,17,1:6)=(/ 2061.500000,1960.300000,1661.600000,1177.400000,459.850000,370.820000 /)
XPIZA_LKT(16,17,1:6)=(/ 0.962996,0.988523,0.999619,0.999641,0.997986,0.477727 /)
XCGA_LKT(16,17,1:6)=(/ 0.728267,0.716173,0.709090,0.710820,0.622500,0.441567 /)
XEXT_COEFF_550_LKT(16,17)=1852.900000 !rg=0.0345903 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,18,1:6)=(/ 1704.000000,1658.300000,1484.400000,1161.500000,558.770000,432.530000 /)
XPIZA_LKT(16,18,1:6)=(/ 0.957645,0.986496,0.999560,0.999619,0.998228,0.532603 /)
XCGA_LKT(16,18,1:6)=(/ 0.735083,0.717907,0.708317,0.715693,0.666177,0.512560 /)
XEXT_COEFF_550_LKT(16,18)=1613.600000 !rg=0.0345903 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,19,1:6)=(/ 1413.500000,1398.000000,1308.100000,1074.700000,635.370000,487.180000 /)
XPIZA_LKT(16,19,1:6)=(/ 0.949486,0.984392,0.999486,0.999572,0.998347,0.572204 /)
XCGA_LKT(16,19,1:6)=(/ 0.741253,0.723857,0.711507,0.710890,0.696620,0.576357 /)
XEXT_COEFF_550_LKT(16,19)=1371.800000 !rg=0.0345903 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,20,1:6)=(/ 1159.700000,1169.000000,1116.300000,958.600000,671.880000,522.180000 /)
XPIZA_LKT(16,20,1:6)=(/ 0.941181,0.981109,0.999397,0.999507,0.998359,0.597429 /)
XCGA_LKT(16,20,1:6)=(/ 0.747183,0.727270,0.713353,0.705410,0.713597,0.622953 /)
XEXT_COEFF_550_LKT(16,20)=1159.100000 !rg=0.0345903 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,1,1:6)=(/ 521.050000,116.910000,19.695000,2.940100,0.814940,133.120000 /)
XPIZA_LKT(17,1,1:6)=(/ 0.919839,0.910703,0.985263,0.929793,0.900,0.000448 /)
XCGA_LKT(17,1,1:6)=(/ 0.109473,0.051093,0.021930,0.008270,0.900,0.000997 /)
XEXT_COEFF_550_LKT(17,1)=54.342000 !rg=0.0374736 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,2,1:6)=(/ 630.930000,144.250000,24.750000,3.653700,0.889740,133.190000 /)
XPIZA_LKT(17,2,1:6)=(/ 0.932369,0.926546,0.988185,0.943335,0.900,0.000565 /)
XCGA_LKT(17,2,1:6)=(/ 0.139027,0.064767,0.027823,0.010503,0.900,0.001270 /)
XEXT_COEFF_550_LKT(17,2)=68.207000 !rg=0.0374736 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,3,1:6)=(/ 855.530000,208.040000,37.020000,5.395400,1.072200,133.360000 /)
XPIZA_LKT(17,3,1:6)=(/ 0.948140,0.947586,0.991978,0.961382,0.504073,0.000849 /)
XCGA_LKT(17,3,1:6)=(/ 0.208767,0.098017,0.042273,0.016003,0.005190,0.001937 /)
XEXT_COEFF_550_LKT(17,3)=101.370000 !rg=0.0374736 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,4,1:6)=(/ 1212.800000,326.230000,62.379000,9.078900,1.458400,133.670000 /)
XPIZA_LKT(17,4,1:6)=(/ 0.961180,0.965129,0.995118,0.976798,0.634047,0.001450 /)
XCGA_LKT(17,4,1:6)=(/ 0.322327,0.163493,0.071597,0.027277,0.008870,0.003317 /)
XEXT_COEFF_550_LKT(17,4)=166.680000 !rg=0.0374736 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,5,1:6)=(/ 1719.700000,515.160000,109.460000,16.330000,2.223100,134.200000 /)
XPIZA_LKT(17,5,1:6)=(/ 0.970468,0.976612,0.997116,0.986880,0.758496,0.002637 /)
XCGA_LKT(17,5,1:6)=(/ 0.440147,0.262907,0.119210,0.045687,0.014910,0.005583 /)
XEXT_COEFF_550_LKT(17,5)=277.920000 !rg=0.0374736 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,6,1:6)=(/ 2339.900000,786.260000,186.060000,29.469000,3.630500,135.060000 /)
XPIZA_LKT(17,6,1:6)=(/ 0.976534,0.983541,0.998223,0.992553,0.850809,0.004808 /)
XCGA_LKT(17,6,1:6)=(/ 0.533350,0.372300,0.185973,0.071390,0.023367,0.008767 /)
XEXT_COEFF_550_LKT(17,6)=445.410000 !rg=0.0374736 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,7,1:6)=(/ 3018.200000,1143.000000,301.430000,52.325000,6.169900,136.420000 /)
XPIZA_LKT(17,7,1:6)=(/ 0.980416,0.987739,0.998836,0.995667,0.911075,0.008691 /)
XCGA_LKT(17,7,1:6)=(/ 0.603030,0.469270,0.276343,0.106957,0.035037,0.013177 /)
XEXT_COEFF_550_LKT(17,7)=684.760000 !rg=0.0374736 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,8,1:6)=(/ 3653.700000,1571.000000,468.250000,89.629000,10.647000,138.570000 /)
XPIZA_LKT(17,8,1:6)=(/ 0.982738,0.990329,0.999193,0.997363,0.947530,0.015447 /)
XCGA_LKT(17,8,1:6)=(/ 0.653227,0.550057,0.379873,0.157137,0.051343,0.019353 /)
XEXT_COEFF_550_LKT(17,8)=989.910000 !rg=0.0374736 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,9,1:6)=(/ 4146.800000,2032.300000,691.850000,146.830000,18.524000,142.010000 /)
XPIZA_LKT(17,9,1:6)=(/ 0.983922,0.991936,0.999409,0.998305,0.969075,0.027131 /)
XCGA_LKT(17,9,1:6)=(/ 0.687910,0.611800,0.468193,0.229673,0.074607,0.028190 /)
XEXT_COEFF_550_LKT(17,9)=1357.600000 !rg=0.0374736 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,10,1:6)=(/ 4384.600000,2465.300000,960.690000,229.300000,31.979000,147.530000 /)
XPIZA_LKT(17,10,1:6)=(/ 0.984091,0.992892,0.999538,0.998840,0.981480,0.046739 /)
XCGA_LKT(17,10,1:6)=(/ 0.708460,0.658190,0.544980,0.329590,0.108130,0.040910 /)
XEXT_COEFF_550_LKT(17,10)=1744.000000 !rg=0.0374736 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,11,1:6)=(/ 4325.300000,2801.800000,1260.600000,350.110000,54.204000,156.480000 /)
XPIZA_LKT(17,11,1:6)=(/ 0.983259,0.993369,0.999620,0.999175,0.988599,0.078794 /)
XCGA_LKT(17,11,1:6)=(/ 0.717183,0.690290,0.607113,0.434153,0.157323,0.059453 /)
XEXT_COEFF_550_LKT(17,11)=2107.100000 !rg=0.0374736 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,12,1:6)=(/ 4024.100000,2962.700000,1554.800000,502.840000,88.126000,170.880000 /)
XPIZA_LKT(17,12,1:6)=(/ 0.981413,0.993409,0.999670,0.999379,0.992611,0.127916 /)
XCGA_LKT(17,12,1:6)=(/ 0.717840,0.708857,0.654800,0.511613,0.230163,0.086710 /)
XEXT_COEFF_550_LKT(17,12)=2377.900000 !rg=0.0374736 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,13,1:6)=(/ 3606.900000,2910.700000,1791.400000,681.450000,136.660000,193.450000 /)
XPIZA_LKT(17,13,1:6)=(/ 0.978914,0.993009,0.999695,0.999502,0.994904,0.196192 /)
XCGA_LKT(17,13,1:6)=(/ 0.718160,0.715133,0.688067,0.582487,0.333737,0.127373 /)
XEXT_COEFF_550_LKT(17,13)=2493.100000 !rg=0.0374736 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,14,1:6)=(/ 3177.500000,2689.600000,1917.500000,864.630000,207.250000,226.740000 /)
XPIZA_LKT(17,14,1:6)=(/ 0.975748,0.992201,0.999700,0.999579,0.996338,0.278121 /)
XCGA_LKT(17,14,1:6)=(/ 0.721510,0.713513,0.707560,0.635360,0.444490,0.189017 /)
XEXT_COEFF_550_LKT(17,14)=2426.200000 !rg=0.0374736 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,15,1:6)=(/ 2719.300000,2413.100000,1905.400000,1030.000000,296.490000,272.650000 /)
XPIZA_LKT(17,15,1:6)=(/ 0.971493,0.991084,0.999684,0.999623,0.997233,0.361406 /)
XCGA_LKT(17,15,1:6)=(/ 0.722970,0.712623,0.714283,0.675523,0.520700,0.280547 /)
XEXT_COEFF_550_LKT(17,15)=2229.200000 !rg=0.0374736 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,16,1:6)=(/ 2286.300000,2135.800000,1771.700000,1145.800000,398.330000,332.580000 /)
XPIZA_LKT(17,16,1:6)=(/ 0.966134,0.989624,0.999648,0.999641,0.997768,0.437833 /)
XCGA_LKT(17,16,1:6)=(/ 0.725333,0.714600,0.712073,0.701543,0.589727,0.389323 /)
XEXT_COEFF_550_LKT(17,16)=1991.600000 !rg=0.0374736 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,17,1:6)=(/ 1903.500000,1832.600000,1591.700000,1181.800000,501.610000,396.990000 /)
XPIZA_LKT(17,17,1:6)=(/ 0.961326,0.987737,0.999595,0.999634,0.998100,0.502772 /)
XCGA_LKT(17,17,1:6)=(/ 0.732633,0.716990,0.707737,0.714050,0.641397,0.472327 /)
XEXT_COEFF_550_LKT(17,17)=1762.400000 !rg=0.0374736 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,18,1:6)=(/ 1581.300000,1545.100000,1413.900000,1132.300000,593.350000,456.840000 /)
XPIZA_LKT(17,18,1:6)=(/ 0.954163,0.986043,0.999534,0.999603,0.998290,0.550523 /)
XCGA_LKT(17,18,1:6)=(/ 0.737320,0.722467,0.709800,0.714450,0.680260,0.540960 /)
XEXT_COEFF_550_LKT(17,18)=1511.100000 !rg=0.0374736 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,19,1:6)=(/ 1299.500000,1301.300000,1230.300000,1026.200000,655.520000,503.960000 /)
XPIZA_LKT(17,19,1:6)=(/ 0.945968,0.983056,0.999454,0.999548,0.998364,0.584349 /)
XCGA_LKT(17,19,1:6)=(/ 0.742250,0.725117,0.713153,0.708387,0.704910,0.596957 /)
XEXT_COEFF_550_LKT(17,19)=1278.500000 !rg=0.0374736 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,20,1:6)=(/ 1062.400000,1080.100000,1049.100000,914.230000,672.450000,529.500000 /)
XPIZA_LKT(17,20,1:6)=(/ 0.937880,0.980020,0.999292,0.999479,0.998327,0.604219 /)
XCGA_LKT(17,20,1:6)=(/ 0.753787,0.727717,0.713730,0.706490,0.716307,0.638713 /)
XEXT_COEFF_550_LKT(17,20)=1074.600000 !rg=0.0374736 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,1,1:6)=(/ 644.280000,145.760000,24.967000,3.683500,0.892930,133.210000 /)
XPIZA_LKT(18,1,1:6)=(/ 0.933297,0.927063,0.988270,0.943762,0.900,0.000570 /)
XCGA_LKT(18,1,1:6)=(/ 0.128767,0.059893,0.025710,0.009700,0.900,0.001170 /)
XEXT_COEFF_550_LKT(18,1)=68.859000 !rg=0.0405973 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,2,1:6)=(/ 774.220000,180.060000,31.387000,4.591000,0.988000,133.300000 /)
XPIZA_LKT(18,2,1:6)=(/ 0.943257,0.939979,0.990589,0.954719,0.900,0.000718 /)
XCGA_LKT(18,2,1:6)=(/ 0.163707,0.075923,0.032617,0.012320,0.900,0.001490 /)
XEXT_COEFF_550_LKT(18,2)=86.395000 !rg=0.0405973 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,3,1:6)=(/ 1030.700000,258.720000,46.923000,6.805400,1.219900,133.500000 /)
XPIZA_LKT(18,3,1:6)=(/ 0.955585,0.956898,0.993593,0.969225,0.563375,0.001079 /)
XCGA_LKT(18,3,1:6)=(/ 0.244417,0.114793,0.049520,0.018767,0.006087,0.002273 /)
XEXT_COEFF_550_LKT(18,3)=127.970000 !rg=0.0405973 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,4,1:6)=(/ 1433.100000,399.610000,78.688000,11.483000,1.710800,133.870000 /)
XPIZA_LKT(18,4,1:6)=(/ 0.965944,0.970792,0.996070,0.981530,0.687314,0.001842 /)
XCGA_LKT(18,4,1:6)=(/ 0.364270,0.190180,0.083723,0.031960,0.010403,0.003890 /)
XEXT_COEFF_550_LKT(18,4)=207.800000 !rg=0.0405973 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,5,1:6)=(/ 1987.600000,619.020000,136.310000,20.660000,2.682800,134.510000 /)
XPIZA_LKT(18,5,1:6)=(/ 0.973476,0.979938,0.997640,0.989533,0.799212,0.003347 /)
XCGA_LKT(18,5,1:6)=(/ 0.476987,0.298393,0.139167,0.053493,0.017483,0.006553 /)
XEXT_COEFF_550_LKT(18,5)=339.720000 !rg=0.0405973 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,6,1:6)=(/ 2637.200000,924.260000,227.420000,37.158000,4.471000,135.540000 /)
XPIZA_LKT(18,6,1:6)=(/ 0.978430,0.985522,0.998512,0.994019,0.878274,0.006097 /)
XCGA_LKT(18,6,1:6)=(/ 0.562130,0.407620,0.216547,0.083563,0.027387,0.010287 /)
XEXT_COEFF_550_LKT(18,6)=534.760000 !rg=0.0405973 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,7,1:6)=(/ 3310.600000,1315.200000,362.860000,65.402000,7.693100,137.180000 /)
XPIZA_LKT(18,7,1:6)=(/ 0.981585,0.988959,0.999004,0.996476,0.928197,0.011001 /)
XCGA_LKT(18,7,1:6)=(/ 0.624640,0.501723,0.316170,0.125257,0.041050,0.015457 /)
XEXT_COEFF_550_LKT(18,7)=802.860000 !rg=0.0405973 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,8,1:6)=(/ 3898.700000,1761.100000,553.970000,110.250000,13.357000,139.790000 /)
XPIZA_LKT(18,8,1:6)=(/ 0.983400,0.991085,0.999295,0.997812,0.957777,0.019495 /)
XCGA_LKT(18,8,1:6)=(/ 0.668637,0.575310,0.415850,0.184213,0.060143,0.022703 /)
XEXT_COEFF_550_LKT(18,8)=1135.800000 !rg=0.0405973 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,9,1:6)=(/ 4298.500000,2222.700000,796.380000,177.100000,23.249000,143.990000 /)
XPIZA_LKT(18,9,1:6)=(/ 0.984153,0.992406,0.999470,0.998558,0.975043,0.034068 /)
XCGA_LKT(18,9,1:6)=(/ 0.697733,0.631803,0.498877,0.268530,0.087417,0.033067 /)
XEXT_COEFF_550_LKT(18,9)=1516.800000 !rg=0.0405973 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,10,1:6)=(/ 4415.100000,2624.800000,1083.200000,273.830000,39.882000,150.710000 /)
XPIZA_LKT(18,10,1:6)=(/ 0.983905,0.993148,0.999577,0.998996,0.984904,0.058204 /)
XCGA_LKT(18,10,1:6)=(/ 0.713367,0.672367,0.571997,0.374987,0.126830,0.048003 /)
XEXT_COEFF_550_LKT(18,10)=1901.600000 !rg=0.0405973 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,11,1:6)=(/ 4243.400000,2896.300000,1384.800000,409.930000,66.578000,161.620000 /)
XPIZA_LKT(18,11,1:6)=(/ 0.982662,0.993445,0.999643,0.999272,0.990528,0.096800 /)
XCGA_LKT(18,11,1:6)=(/ 0.718517,0.699113,0.627417,0.467530,0.184890,0.069827 /)
XEXT_COEFF_550_LKT(18,11)=2235.500000 !rg=0.0405973 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,12,1:6)=(/ 3873.600000,2970.700000,1661.200000,572.080000,105.970000,179.030000 /)
XPIZA_LKT(18,12,1:6)=(/ 0.980549,0.993308,0.999682,0.999435,0.993698,0.153879 /)
XCGA_LKT(18,12,1:6)=(/ 0.718367,0.712577,0.669437,0.540960,0.270350,0.102063 /)
XEXT_COEFF_550_LKT(18,12)=2448.000000 !rg=0.0405973 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,13,1:6)=(/ 3447.500000,2840.600000,1859.800000,757.260000,162.630000,205.780000 /)
XPIZA_LKT(18,13,1:6)=(/ 0.977765,0.992742,0.999699,0.999539,0.995568,0.229039 /)
XCGA_LKT(18,13,1:6)=(/ 0.720237,0.715243,0.697310,0.605863,0.382703,0.150510 /)
XEXT_COEFF_550_LKT(18,13)=2488.700000 !rg=0.0405973 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,14,1:6)=(/ 2999.000000,2587.100000,1931.300000,937.780000,242.940000,244.090000 /)
XPIZA_LKT(18,14,1:6)=(/ 0.974091,0.991778,0.999696,0.999600,0.996768,0.313107 /)
XCGA_LKT(18,14,1:6)=(/ 0.722313,0.713113,0.711713,0.653433,0.479257,0.223997 /)
XEXT_COEFF_550_LKT(18,14)=2359.800000 !rg=0.0405973 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,15,1:6)=(/ 2533.400000,2307.300000,1863.300000,1086.100000,336.110000,295.770000 /)
XPIZA_LKT(18,15,1:6)=(/ 0.969508,0.990542,0.999672,0.999633,0.997480,0.393850 /)
XCGA_LKT(18,15,1:6)=(/ 0.723877,0.713090,0.714130,0.687543,0.549067,0.326913 /)
XEXT_COEFF_550_LKT(18,15)=2136.300000 !rg=0.0405973 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,16,1:6)=(/ 2129.200000,2014.000000,1703.700000,1171.400000,442.070000,359.610000 /)
XPIZA_LKT(18,16,1:6)=(/ 0.963957,0.989029,0.999627,0.999641,0.997929,0.466631 /)
XCGA_LKT(18,16,1:6)=(/ 0.729690,0.717423,0.709420,0.708197,0.613513,0.427837 /)
XEXT_COEFF_550_LKT(18,16)=1905.200000 !rg=0.0405973 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,17,1:6)=(/ 1760.000000,1710.100000,1518.100000,1170.500000,542.950000,422.240000 /)
XPIZA_LKT(18,17,1:6)=(/ 0.958415,0.987072,0.999572,0.999624,0.998194,0.524512 /)
XCGA_LKT(18,17,1:6)=(/ 0.733683,0.720117,0.708427,0.715330,0.659190,0.500990 /)
XEXT_COEFF_550_LKT(18,17)=1661.100000 !rg=0.0405973 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,18,1:6)=(/ 1454.700000,1445.500000,1343.400000,1092.300000,624.160000,478.840000 /)
XPIZA_LKT(18,18,1:6)=(/ 0.950595,0.984743,0.999503,0.999583,0.998331,0.566274 /)
XCGA_LKT(18,18,1:6)=(/ 0.738360,0.723150,0.712233,0.712350,0.691957,0.566740 /)
XEXT_COEFF_550_LKT(18,18)=1409.700000 !rg=0.0405973 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,19,1:6)=(/ 1195.500000,1201.900000,1154.600000,980.300000,668.680000,517.490000 /)
XPIZA_LKT(18,19,1:6)=(/ 0.942079,0.981520,0.999398,0.999520,0.998361,0.593925 /)
XCGA_LKT(18,19,1:6)=(/ 0.748330,0.724397,0.713287,0.706940,0.711287,0.615563 /)
XEXT_COEFF_550_LKT(18,19)=1198.500000 !rg=0.0405973 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,20,1:6)=(/ 977.110000,993.880000,974.590000,873.790000,664.270000,531.950000 /)
XPIZA_LKT(18,20,1:6)=(/ 0.933573,0.978969,0.999300,0.999444,0.998278,0.608548 /)
XCGA_LKT(18,20,1:6)=(/ 0.756980,0.734473,0.716830,0.706690,0.717177,0.652037 /)
XEXT_COEFF_550_LKT(18,20)=990.130000 !rg=0.0405973 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,1,1:6)=(/ 794.180000,182.170000,31.667000,4.628900,0.992050,133.310000 /)
XPIZA_LKT(19,1,1:6)=(/ 0.944238,0.940438,0.990656,0.955060,0.900,0.000724 /)
XCGA_LKT(19,1,1:6)=(/ 0.151773,0.070223,0.030140,0.011380,0.900,0.001373 /)
XEXT_COEFF_550_LKT(19,1)=87.266000 !rg=0.0439813 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,2,1:6)=(/ 943.910000,224.970000,39.813000,5.782800,1.112900,133.420000 /)
XPIZA_LKT(19,2,1:6)=(/ 0.952020,0.950890,0.992495,0.963880,0.521824,0.000912 /)
XCGA_LKT(19,2,1:6)=(/ 0.193187,0.089030,0.038227,0.014450,0.004683,0.001747 /)
XEXT_COEFF_550_LKT(19,2)=109.390000 !rg=0.0439813 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,3,1:6)=(/ 1232.100000,320.850000,59.445000,8.597600,1.407700,133.660000 /)
XPIZA_LKT(19,3,1:6)=(/ 0.961573,0.964386,0.994872,0.975496,0.620852,0.001370 /)
XCGA_LKT(19,3,1:6)=(/ 0.285383,0.134470,0.058000,0.022000,0.007143,0.002667 /)
XEXT_COEFF_550_LKT(19,3)=161.170000 !rg=0.0439813 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,4,1:6)=(/ 1681.200000,486.670000,98.999000,14.533000,2.031700,134.110000 /)
XPIZA_LKT(19,4,1:6)=(/ 0.969846,0.975341,0.996822,0.985291,0.735976,0.002339 /)
XCGA_LKT(19,4,1:6)=(/ 0.407423,0.220450,0.097853,0.037440,0.012203,0.004567 /)
XEXT_COEFF_550_LKT(19,4)=257.710000 !rg=0.0439813 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,5,1:6)=(/ 2277.600000,738.350000,168.820000,26.120000,3.267000,134.880000 /)
XPIZA_LKT(19,5,1:6)=(/ 0.975965,0.982633,0.998054,0.991632,0.834476,0.004247 /)
XCGA_LKT(19,5,1:6)=(/ 0.512010,0.335050,0.162277,0.062613,0.020497,0.007687 /)
XEXT_COEFF_550_LKT(19,5)=412.820000 !rg=0.0439813 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,6,1:6)=(/ 2944.200000,1078.200000,276.410000,46.727000,5.537900,136.120000 /)
XPIZA_LKT(19,6,1:6)=(/ 0.980011,0.987146,0.998743,0.995177,0.901182,0.007727 /)
XCGA_LKT(19,6,1:6)=(/ 0.588750,0.442837,0.250753,0.097797,0.032093,0.012067 /)
XEXT_COEFF_550_LKT(19,6)=636.640000 !rg=0.0439813 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,7,1:6)=(/ 3596.300000,1499.100000,434.570000,81.280000,9.622200,138.110000 /)
XPIZA_LKT(19,7,1:6)=(/ 0.982535,0.989971,0.999142,0.997113,0.942143,0.013910 /)
XCGA_LKT(19,7,1:6)=(/ 0.644117,0.531737,0.355160,0.146683,0.048090,0.018130 /)
XEXT_COEFF_550_LKT(19,7)=932.910000 !rg=0.0439813 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,8,1:6)=(/ 4116.900000,1957.500000,647.500000,134.500000,16.769000,141.280000 /)
XPIZA_LKT(19,8,1:6)=(/ 0.983874,0.991710,0.999377,0.998167,0.966005,0.024559 /)
XCGA_LKT(19,8,1:6)=(/ 0.682153,0.598797,0.448533,0.215787,0.070453,0.026627 /)
XEXT_COEFF_550_LKT(19,8)=1291.100000 !rg=0.0439813 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,9,1:6)=(/ 4406.100000,2404.400000,910.580000,212.530000,29.123000,146.390000 /)
XPIZA_LKT(19,9,1:6)=(/ 0.984221,0.992772,0.999519,0.998764,0.979792,0.042651 /)
XCGA_LKT(19,9,1:6)=(/ 0.705907,0.649340,0.529507,0.311477,0.102450,0.038787 /)
XEXT_COEFF_550_LKT(19,9)=1680.100000 !rg=0.0439813 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,10,1:6)=(/ 4396.300000,2761.300000,1206.100000,325.560000,49.468000,154.590000 /)
XPIZA_LKT(19,10,1:6)=(/ 0.983578,0.993324,0.999608,0.999126,0.987609,0.072127 /)
XCGA_LKT(19,10,1:6)=(/ 0.717030,0.684340,0.595153,0.416090,0.148853,0.056340 /)
XEXT_COEFF_550_LKT(19,10)=2049.800000 !rg=0.0439813 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,11,1:6)=(/ 4125.900000,2958.800000,1504.800000,472.900000,81.048000,167.850000 /)
XPIZA_LKT(19,11,1:6)=(/ 0.981963,0.993443,0.999662,0.999349,0.992047,0.118002 /)
XCGA_LKT(19,11,1:6)=(/ 0.719230,0.705763,0.645787,0.497160,0.217367,0.082070 /)
XEXT_COEFF_550_LKT(19,11)=2342.600000 !rg=0.0439813 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,12,1:6)=(/ 3720.000000,2946.400000,1755.900000,646.920000,126.610000,188.750000 /)
XPIZA_LKT(19,12,1:6)=(/ 0.979575,0.993136,0.999692,0.999483,0.994572,0.182995 /)
XCGA_LKT(19,12,1:6)=(/ 0.719753,0.714723,0.681950,0.569787,0.315750,0.120303 /)
XEXT_COEFF_550_LKT(19,12)=2490.400000 !rg=0.0439813 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,13,1:6)=(/ 3279.400000,2755.200000,1906.800000,831.770000,193.180000,220.090000 /)
XPIZA_LKT(19,13,1:6)=(/ 0.976390,0.992391,0.999700,0.999567,0.996131,0.263434 /)
XCGA_LKT(19,13,1:6)=(/ 0.721280,0.714747,0.704407,0.625980,0.427697,0.178150 /)
XEXT_COEFF_550_LKT(19,13)=2458.000000 !rg=0.0439813 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,14,1:6)=(/ 2807.600000,2479.800000,1920.300000,1003.300000,280.050000,263.710000 /)
XPIZA_LKT(19,14,1:6)=(/ 0.972487,0.991307,0.999688,0.999617,0.997110,0.347298 /)
XCGA_LKT(19,14,1:6)=(/ 0.723613,0.712637,0.713600,0.668743,0.508320,0.264653 /)
XEXT_COEFF_550_LKT(19,14)=2279.000000 !rg=0.0439813 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,15,1:6)=(/ 2369.900000,2194.700000,1807.900000,1129.300000,379.040000,321.170000 /)
XPIZA_LKT(19,15,1:6)=(/ 0.966890,0.990014,0.999655,0.999639,0.997687,0.424988 /)
XCGA_LKT(19,15,1:6)=(/ 0.726427,0.716127,0.712330,0.697097,0.578007,0.372577 /)
XEXT_COEFF_550_LKT(19,15)=2049.300000 !rg=0.0439813 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,16,1:6)=(/ 1967.100000,1894.200000,1628.000000,1181.600000,484.010000,386.050000 /)
XPIZA_LKT(19,16,1:6)=(/ 0.961854,0.988018,0.999607,0.999637,0.998054,0.492800 /)
XCGA_LKT(19,16,1:6)=(/ 0.730920,0.718127,0.708833,0.712217,0.633183,0.459993 /)
XEXT_COEFF_550_LKT(19,16)=1810.200000 !rg=0.0439813 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,17,1:6)=(/ 1624.500000,1594.600000,1450.000000,1146.100000,579.450000,446.960000 /)
XPIZA_LKT(19,17,1:6)=(/ 0.955470,0.986279,0.999547,0.999610,0.998266,0.543370 /)
XCGA_LKT(19,17,1:6)=(/ 0.735497,0.720717,0.710540,0.714927,0.674357,0.529603 /)
XEXT_COEFF_550_LKT(19,17)=1557.400000 !rg=0.0439813 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,18,1:6)=(/ 1342.600000,1339.700000,1269.600000,1049.400000,647.260000,497.050000 /)
XPIZA_LKT(19,18,1:6)=(/ 0.946717,0.983533,0.999463,0.999558,0.998358,0.579474 /)
XCGA_LKT(19,18,1:6)=(/ 0.743023,0.723487,0.712570,0.709577,0.701167,0.588410 /)
XEXT_COEFF_550_LKT(19,18)=1324.500000 !rg=0.0439813 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,19,1:6)=(/ 1094.700000,1115.600000,1075.800000,935.800000,672.440000,526.790000 /)
XPIZA_LKT(19,19,1:6)=(/ 0.938662,0.979827,0.999350,0.999489,0.998339,0.601546 /)
XCGA_LKT(19,19,1:6)=(/ 0.751413,0.729027,0.712670,0.705517,0.714793,0.632213 /)
XEXT_COEFF_550_LKT(19,19)=1105.700000 !rg=0.0439813 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,20,1:6)=(/ 897.040000,914.430000,907.260000,829.710000,648.420000,530.050000 /)
XPIZA_LKT(19,20,1:6)=(/ 0.928517,0.977544,0.999196,0.999414,0.998207,0.610970 /)
XCGA_LKT(19,20,1:6)=(/ 0.759373,0.736263,0.720420,0.710657,0.716233,0.663627 /)
XEXT_COEFF_550_LKT(19,20)=916.850000 !rg=0.0439813 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,1,1:6)=(/ 973.270000,227.980000,40.179000,5.831100,1.118000,133.440000 /)
XPIZA_LKT(20,1,1:6)=(/ 0.953062,0.951311,0.992548,0.964152,0.523902,0.000920 /)
XCGA_LKT(20,1,1:6)=(/ 0.179393,0.082363,0.035330,0.013347,0.004323,0.001613 /)
XEXT_COEFF_550_LKT(20,1)=110.570000 !rg=0.0476474 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,2,1:6)=(/ 1140.600000,280.910000,50.503000,7.298300,1.271600,133.570000 /)
XPIZA_LKT(20,2,1:6)=(/ 0.959016,0.959708,0.994004,0.971222,0.580732,0.001158 /)
XCGA_LKT(20,2,1:6)=(/ 0.228527,0.104460,0.044800,0.016947,0.005493,0.002050 /)
XEXT_COEFF_550_LKT(20,2)=138.370000 !rg=0.0476474 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,3,1:6)=(/ 1461.200000,396.180000,75.235000,10.875000,1.646300,133.860000 /)
XPIZA_LKT(20,3,1:6)=(/ 0.966396,0.970389,0.995884,0.980493,0.675053,0.001740 /)
XCGA_LKT(20,3,1:6)=(/ 0.331403,0.157570,0.067927,0.025790,0.008380,0.003130 /)
XEXT_COEFF_550_LKT(20,3)=202.350000 !rg=0.0476474 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,4,1:6)=(/ 1957.600000,588.770000,124.120000,18.399000,2.439600,134.400000 /)
XPIZA_LKT(20,4,1:6)=(/ 0.973057,0.979000,0.997416,0.988276,0.779415,0.002969 /)
XCGA_LKT(20,4,1:6)=(/ 0.450413,0.254290,0.114310,0.043847,0.014310,0.005357 /)
XEXT_COEFF_550_LKT(20,4)=317.700000 !rg=0.0476474 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,5,1:6)=(/ 2587.400000,874.000000,207.840000,32.983000,4.009200,135.320000 /)
XPIZA_LKT(20,5,1:6)=(/ 0.978037,0.984821,0.998382,0.993293,0.864510,0.005386 /)
XCGA_LKT(20,5,1:6)=(/ 0.544930,0.372657,0.188827,0.073267,0.024023,0.009020 /)
XEXT_COEFF_550_LKT(20,5)=498.230000 !rg=0.0476474 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,6,1:6)=(/ 3253.700000,1248.400000,334.270000,58.542000,6.891600,136.820000 /)
XPIZA_LKT(20,6,1:6)=(/ 0.981310,0.988489,0.998931,0.996089,0.920083,0.009785 /)
XCGA_LKT(20,6,1:6)=(/ 0.612953,0.477320,0.287397,0.114423,0.037600,0.014157 /)
XEXT_COEFF_550_LKT(20,6)=750.930000 !rg=0.0476474 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,7,1:6)=(/ 3863.600000,1691.200000,515.990000,100.320000,12.061000,139.240000 /)
XPIZA_LKT(20,7,1:6)=(/ 0.983294,0.990801,0.999254,0.997615,0.953428,0.017566 /)
XCGA_LKT(20,7,1:6)=(/ 0.661383,0.559013,0.391427,0.171727,0.056330,0.021267 /)
XEXT_COEFF_550_LKT(20,7)=1076.300000 !rg=0.0476474 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,8,1:6)=(/ 4296.400000,2153.700000,749.160000,162.860000,21.047000,143.090000 /)
XPIZA_LKT(20,8,1:6)=(/ 0.984197,0.992229,0.999443,0.998448,0.972585,0.030870 /)
XCGA_LKT(20,8,1:6)=(/ 0.693540,0.620290,0.480253,0.252087,0.082533,0.031233 /)
XEXT_COEFF_550_LKT(20,8)=1451.000000 !rg=0.0476474 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,9,1:6)=(/ 4464.400000,2575.100000,1031.800000,254.350000,36.363000,149.330000 /)
XPIZA_LKT(20,9,1:6)=(/ 0.984131,0.993064,0.999561,0.998933,0.983558,0.053200 /)
XCGA_LKT(20,9,1:6)=(/ 0.712133,0.664830,0.557707,0.355673,0.120110,0.045510 /)
XEXT_COEFF_550_LKT(20,9)=1842.200000 !rg=0.0476474 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,10,1:6)=(/ 4335.600000,2873.000000,1331.100000,382.860000,60.923000,159.320000 /)
XPIZA_LKT(20,10,1:6)=(/ 0.983086,0.993433,0.999633,0.999232,0.989741,0.088847 /)
XCGA_LKT(20,10,1:6)=(/ 0.719227,0.694447,0.616403,0.450903,0.174810,0.066157 /)
XEXT_COEFF_550_LKT(20,10)=2187.700000 !rg=0.0476474 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,11,1:6)=(/ 3989.100000,2986.600000,1615.900000,540.040000,97.784000,175.370000 /)
XPIZA_LKT(20,11,1:6)=(/ 0.981155,0.993380,0.999677,0.999411,0.993246,0.142520 /)
XCGA_LKT(20,11,1:6)=(/ 0.720393,0.710730,0.661510,0.526680,0.255283,0.096553 /)
XEXT_COEFF_550_LKT(20,11)=2427.400000 !rg=0.0476474 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,12,1:6)=(/ 3559.600000,2893.900000,1832.900000,722.930000,150.860000,200.250000 /)
XPIZA_LKT(20,12,1:6)=(/ 0.978419,0.992896,0.999698,0.999523,0.995292,0.214788 /)
XCGA_LKT(20,12,1:6)=(/ 0.720540,0.715697,0.692190,0.594653,0.363890,0.142037 /)
XEXT_COEFF_550_LKT(20,12)=2502.500000 !rg=0.0476474 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,13,1:6)=(/ 3096.900000,2654.800000,1932.000000,906.500000,227.490000,236.490000 /)
XPIZA_LKT(20,13,1:6)=(/ 0.974971,0.991995,0.999698,0.999591,0.996598,0.298339 /)
XCGA_LKT(20,13,1:6)=(/ 0.723387,0.713767,0.709350,0.644860,0.464780,0.211033 /)
XEXT_COEFF_550_LKT(20,13)=2403.200000 !rg=0.0476474 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,14,1:6)=(/ 2628.300000,2370.900000,1889.500000,1063.100000,318.630000,285.750000 /)
XPIZA_LKT(20,14,1:6)=(/ 0.970021,0.990835,0.999678,0.999628,0.997379,0.380258 /)
XCGA_LKT(20,14,1:6)=(/ 0.724400,0.714530,0.714053,0.681683,0.536537,0.309683 /)
XEXT_COEFF_550_LKT(20,14)=2194.300000 !rg=0.0476474 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,15,1:6)=(/ 2198.200000,2081.500000,1739.200000,1161.600000,423.030000,347.840000 /)
XPIZA_LKT(20,15,1:6)=(/ 0.965229,0.989343,0.999639,0.999641,0.997862,0.454574 /)
XCGA_LKT(20,15,1:6)=(/ 0.728060,0.717387,0.711323,0.704767,0.603463,0.412940 /)
XEXT_COEFF_550_LKT(20,15)=1954.200000 !rg=0.0476474 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,16,1:6)=(/ 1813.600000,1763.900000,1558.400000,1177.200000,525.990000,411.550000 /)
XPIZA_LKT(20,16,1:6)=(/ 0.959595,0.987447,0.999584,0.999629,0.998155,0.515740 /)
XCGA_LKT(20,16,1:6)=(/ 0.732630,0.718903,0.709450,0.714513,0.651560,0.489010 /)
XEXT_COEFF_550_LKT(20,16)=1712.000000 !rg=0.0476474 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,17,1:6)=(/ 1503.800000,1486.300000,1381.000000,1112.200000,612.030000,470.080000 /)
XPIZA_LKT(20,17,1:6)=(/ 0.951633,0.985352,0.999514,0.999591,0.998314,0.559997 /)
XCGA_LKT(20,17,1:6)=(/ 0.739023,0.722380,0.711123,0.713030,0.686877,0.556487 /)
XEXT_COEFF_550_LKT(20,17)=1459.100000 !rg=0.0476474 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,18,1:6)=(/ 1230.300000,1244.200000,1186.700000,1002.300000,663.750000,511.950000 /)
XPIZA_LKT(20,18,1:6)=(/ 0.943585,0.981951,0.999427,0.999531,0.998361,0.589966 /)
XCGA_LKT(20,18,1:6)=(/ 0.745593,0.725817,0.713433,0.706960,0.708350,0.607580 /)
XEXT_COEFF_550_LKT(20,18)=1232.000000 !rg=0.0476474 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,19,1:6)=(/ 1004.400000,1026.000000,1001.400000,892.540000,667.700000,531.070000 /)
XPIZA_LKT(20,19,1:6)=(/ 0.935000,0.978630,0.999343,0.999460,0.998298,0.606823 /)
XCGA_LKT(20,19,1:6)=(/ 0.754897,0.731723,0.719357,0.707677,0.716653,0.646280 /)
XEXT_COEFF_550_LKT(20,19)=1020.700000 !rg=0.0476474 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,20,1:6)=(/ 826.040000,842.780000,845.480000,786.280000,626.530000,523.410000 /)
XPIZA_LKT(20,20,1:6)=(/ 0.923327,0.975913,0.999160,0.999377,0.998114,0.611596 /)
XCGA_LKT(20,20,1:6)=(/ 0.764700,0.737943,0.719927,0.712687,0.713423,0.673390 /)
XEXT_COEFF_550_LKT(20,20)=851.930000 !rg=0.0476474 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,1,1:6)=(/ 1182.100000,285.310000,50.984000,7.359800,1.278100,133.590000 /)
XPIZA_LKT(21,1,1:6)=(/ 0.960110,0.960111,0.994047,0.971437,0.582748,0.001168 /)
XCGA_LKT(21,1,1:6)=(/ 0.212833,0.096667,0.041413,0.015653,0.005073,0.001893 /)
XEXT_COEFF_550_LKT(21,1)=140.020000 !rg=0.0516191 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,2,1:6)=(/ 1363.500000,350.030000,64.048000,9.225100,1.473400,133.750000 /)
XPIZA_LKT(21,2,1:6)=(/ 0.964555,0.966805,0.995201,0.977086,0.637357,0.001471 /)
XCGA_LKT(21,2,1:6)=(/ 0.270910,0.122673,0.052503,0.019870,0.006447,0.002407 /)
XEXT_COEFF_550_LKT(21,2)=174.750000 !rg=0.0516191 sigma=1.15 wvl=0.55


XEXT_COEFF_WVL_LKT(21,3,1:6)=(/ 1719.900000,486.250000,95.074000,13.769000,1.949700,134.090000 /)
XPIZA_LKT(21,3,1:6)=(/ 0.970298,0.975185,0.996684,0.984470,0.724860,0.002210 /)
XCGA_LKT(21,3,1:6)=(/ 0.381353,0.184673,0.079540,0.030227,0.009830,0.003673 /)
XEXT_COEFF_550_LKT(21,3)=252.940000 !rg=0.0516191 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,4,1:6)=(/ 2261.400000,707.190000,154.940000,23.290000,2.958000,134.740000 /)
XPIZA_LKT(21,4,1:6)=(/ 0.975711,0.981946,0.997886,0.990641,0.817396,0.003767 /)
XCGA_LKT(21,4,1:6)=(/ 0.491923,0.291527,0.133440,0.051337,0.016780,0.006287 /)
XEXT_COEFF_550_LKT(21,4)=389.020000 !rg=0.0516191 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,5,1:6)=(/ 2911.300000,1027.000000,254.380000,41.567000,4.951900,135.850000 /)
XPIZA_LKT(21,5,1:6)=(/ 0.979760,0.986610,0.998645,0.994605,0.889726,0.006828 /)
XCGA_LKT(21,5,1:6)=(/ 0.575297,0.410853,0.218837,0.085697,0.028157,0.010583 /)
XEXT_COEFF_550_LKT(21,5)=596.620000 !rg=0.0516191 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,6,1:6)=(/ 3560.100000,1432.300000,401.890000,72.986000,8.607000,137.670000 /)
XPIZA_LKT(21,6,1:6)=(/ 0.982380,0.989600,0.999083,0.996808,0.935538,0.012378 /)
XCGA_LKT(21,6,1:6)=(/ 0.634863,0.509767,0.324547,0.133830,0.044047,0.016607 /)
XEXT_COEFF_550_LKT(21,6)=878.480000 !rg=0.0516191 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,7,1:6)=(/ 4107.200000,1891.500000,606.120000,122.890000,15.135000,140.600000 /)
XPIZA_LKT(21,7,1:6)=(/ 0.983867,0.991489,0.999343,0.998011,0.962508,0.022147 /)
XCGA_LKT(21,7,1:6)=(/ 0.676673,0.584470,0.425377,0.200857,0.065977,0.024943 /)
XEXT_COEFF_550_LKT(21,7)=1230.300000 !rg=0.0516191 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,8,1:6)=(/ 4433.900000,2344.100000,860.710000,196.110000,26.378000,145.300000 /)
XPIZA_LKT(21,8,1:6)=(/ 0.984358,0.992642,0.999498,0.998676,0.977827,0.038697 /)
XCGA_LKT(21,8,1:6)=(/ 0.703203,0.639377,0.511780,0.292490,0.096690,0.036640 /)
XEXT_COEFF_550_LKT(21,8)=1616.200000 !rg=0.0516191 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,9,1:6)=(/ 4473.800000,2725.400000,1155.000000,303.130000,45.184000,152.900000 /)
XPIZA_LKT(21,9,1:6)=(/ 0.983892,0.993278,0.999595,0.999074,0.986537,0.066058 /)
XCGA_LKT(21,9,1:6)=(/ 0.716980,0.678117,0.582227,0.396933,0.140877,0.053407 /)
XEXT_COEFF_550_LKT(21,9)=1996.900000 !rg=0.0516191 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,10,1:6)=(/ 4235.900000,2953.400000,1454.600000,443.880000,74.400000,165.060000 /)
XPIZA_LKT(21,10,1:6)=(/ 0.982499,0.993471,0.999654,0.999316,0.991419,0.108650 /)
XCGA_LKT(21,10,1:6)=(/ 0.721167,0.702327,0.635917,0.481557,0.205350,0.077730 /)
XEXT_COEFF_550_LKT(21,10)=2306.200000 !rg=0.0516191 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,11,1:6)=(/ 3841.400000,2980.700000,1718.200000,613.010000,117.140000,184.380000 /)
XPIZA_LKT(21,11,1:6)=(/ 0.980149,0.993250,0.999688,0.999462,0.994205,0.170265 /)
XCGA_LKT(21,11,1:6)=(/ 0.720780,0.714053,0.675147,0.556180,0.298500,0.113737 /)
XEXT_COEFF_550_LKT(21,11)=2484.400000 !rg=0.0516191 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,12,1:6)=(/ 3385.800000,2816.900000,1891.400000,797.660000,179.490000,213.650000 /)
XPIZA_LKT(21,12,1:6)=(/ 0.977242,0.992595,0.999700,0.999555,0.995900,0.248484 /)
XCGA_LKT(21,12,1:6)=(/ 0.723043,0.715557,0.700503,0.615737,0.409843,0.167973 /)
XEXT_COEFF_550_LKT(21,12)=2486.200000 !rg=0.0516191 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,13,1:6)=(/ 2911.600000,2548.000000,1932.900000,975.180000,263.750000,255.090000 /)
XPIZA_LKT(21,13,1:6)=(/ 0.973069,0.991592,0.999693,0.999610,0.996972,0.332796 /)
XCGA_LKT(21,13,1:6)=(/ 0.723967,0.714650,0.712470,0.661300,0.495237,0.249523 /)
XEXT_COEFF_550_LKT(21,13)=2332.800000 !rg=0.0516191 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,14,1:6)=(/ 2442.000000,2264.200000,1839.800000,1111.800000,360.320000,310.210000 /)
XPIZA_LKT(21,14,1:6)=(/ 0.968210,0.990291,0.999665,0.999637,0.997601,0.411954 /)
XCGA_LKT(21,14,1:6)=(/ 0.725163,0.715827,0.713643,0.692287,0.565753,0.355500 /)
XEXT_COEFF_550_LKT(21,14)=2100.400000 !rg=0.0516191 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,15,1:6)=(/ 2032.700000,1953.200000,1671.500000,1178.300000,465.330000,374.380000 /)
XPIZA_LKT(21,15,1:6)=(/ 0.963446,0.988628,0.999617,0.999639,0.998001,0.481819 /)
XCGA_LKT(21,15,1:6)=(/ 0.731090,0.718557,0.710047,0.709770,0.624240,0.446710 /)
XEXT_COEFF_550_LKT(21,15)=1863.900000 !rg=0.0516191 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,16,1:6)=(/ 1681.100000,1644.200000,1488.100000,1159.400000,564.530000,436.560000 /)
XPIZA_LKT(21,16,1:6)=(/ 0.956211,0.986726,0.999560,0.999616,0.998238,0.535629 /)
XCGA_LKT(21,16,1:6)=(/ 0.735123,0.720963,0.709970,0.714887,0.667830,0.517683 /)
XEXT_COEFF_550_LKT(21,16)=1610.700000 !rg=0.0516191 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,17,1:6)=(/ 1378.900000,1387.600000,1303.500000,1069.000000,638.390000,489.770000 /)
XPIZA_LKT(21,17,1:6)=(/ 0.948417,0.983810,0.999486,0.999570,0.998349,0.574238 /)
XCGA_LKT(21,17,1:6)=(/ 0.739877,0.723987,0.713873,0.710790,0.697097,0.579407 /)
XEXT_COEFF_550_LKT(21,17)=1359.400000 !rg=0.0516191 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,18,1:6)=(/ 1125.400000,1146.300000,1110.300000,956.700000,670.950000,523.180000 /)
XPIZA_LKT(21,18,1:6)=(/ 0.940756,0.981168,0.999368,0.999505,0.998349,0.598398 /)
XCGA_LKT(21,18,1:6)=(/ 0.750733,0.727077,0.715510,0.707183,0.712873,0.625057 /)
XEXT_COEFF_550_LKT(21,18)=1139.300000 !rg=0.0516191 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,19,1:6)=(/ 925.700000,941.940000,935.450000,852.090000,655.310000,531.070000 /)
XPIZA_LKT(21,19,1:6)=(/ 0.929796,0.977953,0.999294,0.999428,0.998234,0.609955 /)
XCGA_LKT(21,19,1:6)=(/ 0.759213,0.734520,0.720007,0.709583,0.716113,0.658537 /)
XEXT_COEFF_550_LKT(21,19)=949.460000 !rg=0.0516191 sigma=2.85 wvl=0.55
 
END SUBROUTINE SALT_OPT_LKT_SET2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET3()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE
 
XEXT_COEFF_WVL_LKT(21,20,1:6)=(/ 754.760000,779.690000,780.170000,737.560000,600.310000,512.970000 /)
XPIZA_LKT(21,20,1:6)=(/ 0.918452,0.973211,0.999124,0.999331,0.998011,0.610514 /)
XCGA_LKT(21,20,1:6)=(/ 0.766737,0.741147,0.719747,0.713717,0.710980,0.681940 /)
XEXT_COEFF_550_LKT(21,20)=782.640000 !rg=0.0516191 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,1,1:6)=(/ 1418.300000,356.570000,64.690000,9.303600,1.481600,133.770000 /)
XPIZA_LKT(22,1,1:6)=(/ 0.965664,0.967206,0.995235,0.977256,0.639257,0.001484 /)
XCGA_LKT(22,1,1:6)=(/ 0.253660,0.113570,0.048540,0.018357,0.005950,0.002220 /)
XEXT_COEFF_550_LKT(22,1)=177.100000 !rg=0.0559219 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,2,1:6)=(/ 1611.900000,434.430000,81.178000,11.675000,1.729800,133.960000 /)
XPIZA_LKT(22,2,1:6)=(/ 0.968910,0.972493,0.996148,0.981759,0.690334,0.001868 /)
XCGA_LKT(22,2,1:6)=(/ 0.321287,0.144250,0.061530,0.023297,0.007563,0.002823 /)
XEXT_COEFF_550_LKT(22,2)=220.150000 !rg=0.0559219 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,3,1:6)=(/ 2010.700000,592.430000,119.880000,17.442000,2.335400,134.380000 /)
XPIZA_LKT(22,3,1:6)=(/ 0.973485,0.979012,0.997318,0.987627,0.769559,0.002805 /)
XCGA_LKT(22,3,1:6)=(/ 0.433110,0.216390,0.093140,0.035420,0.011530,0.004310 /)
XEXT_COEFF_550_LKT(22,3)=314.410000 !rg=0.0559219 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,4,1:6)=(/ 2590.100000,843.150000,192.410000,29.465000,3.616900,135.150000 /)
XPIZA_LKT(22,4,1:6)=(/ 0.977918,0.984327,0.998258,0.992514,0.850010,0.004779 /)
XCGA_LKT(22,4,1:6)=(/ 0.530803,0.331723,0.155620,0.060087,0.019670,0.007377 /)
XEXT_COEFF_550_LKT(22,4)=472.870000 !rg=0.0559219 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,5,1:6)=(/ 3243.300000,1197.500000,309.500000,52.236000,6.148400,136.490000 /)
XPIZA_LKT(22,5,1:6)=(/ 0.981189,0.988081,0.998854,0.995641,0.910647,0.008648 /)
XCGA_LKT(22,5,1:6)=(/ 0.602927,0.448687,0.251800,0.100190,0.032990,0.012413 /)
XEXT_COEFF_550_LKT(22,5)=708.670000 !rg=0.0559219 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,6,1:6)=(/ 3853.100000,1627.900000,479.280000,90.456000,10.777000,138.690000 /)
XPIZA_LKT(22,6,1:6)=(/ 0.983244,0.990514,0.999206,0.997375,0.948083,0.015641 /)
XCGA_LKT(22,6,1:6)=(/ 0.654433,0.539947,0.360837,0.156437,0.051590,0.019480 /)
XEXT_COEFF_550_LKT(22,6)=1020.000000 !rg=0.0559219 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,7,1:6)=(/ 4315.900000,2094.400000,705.320000,149.470000,18.997000,142.270000 /)
XPIZA_LKT(22,7,1:6)=(/ 0.984282,0.992061,0.999416,0.998325,0.969782,0.027866 /)
XCGA_LKT(22,7,1:6)=(/ 0.689823,0.607850,0.458560,0.234340,0.077267,0.029260 /)
XEXT_COEFF_550_LKT(22,7)=1391.400000 !rg=0.0559219 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,8,1:6)=(/ 4522.700000,2526.300000,980.290000,235.340000,32.971000,147.990000 /)
XPIZA_LKT(22,8,1:6)=(/ 0.984365,0.992975,0.999544,0.998862,0.981990,0.048346 /)
XCGA_LKT(22,8,1:6)=(/ 0.710933,0.656403,0.541197,0.334777,0.113300,0.042983 /)
XEXT_COEFF_550_LKT(22,8)=1782.800000 !rg=0.0559219 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,9,1:6)=(/ 4435.500000,2853.700000,1281.100000,357.820000,55.782000,157.260000 /)
XPIZA_LKT(22,9,1:6)=(/ 0.983511,0.993421,0.999623,0.999189,0.988887,0.081569 /)
XCGA_LKT(22,9,1:6)=(/ 0.720437,0.689563,0.604637,0.432820,0.165307,0.062700 /)
XEXT_COEFF_550_LKT(22,9)=2143.700000 !rg=0.0559219 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,10,1:6)=(/ 4114.500000,3001.000000,1570.500000,509.070000,90.068000,172.000000 /)
XPIZA_LKT(22,10,1:6)=(/ 0.981703,0.993445,0.999671,0.999384,0.992744,0.131718 /)
XCGA_LKT(22,10,1:6)=(/ 0.722137,0.708557,0.652830,0.511497,0.241050,0.091407 /)
XEXT_COEFF_550_LKT(22,10)=2405.100000 !rg=0.0559219 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,11,1:6)=(/ 3674.900000,2944.700000,1803.900000,688.710000,139.810000,195.070000 /)
XPIZA_LKT(22,11,1:6)=(/ 0.979185,0.993045,0.999695,0.999506,0.994991,0.200882 /)
XCGA_LKT(22,11,1:6)=(/ 0.722747,0.715773,0.686480,0.582437,0.345253,0.134177 /)
XEXT_COEFF_550_LKT(22,11)=2513.000000 !rg=0.0559219 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,12,1:6)=(/ 3208.200000,2723.200000,1927.700000,873.230000,212.170000,229.080000 /)
XPIZA_LKT(22,12,1:6)=(/ 0.975768,0.992269,0.999699,0.999581,0.996408,0.283105 /)
XCGA_LKT(22,12,1:6)=(/ 0.724077,0.716153,0.706460,0.635373,0.448953,0.198847 /)
XEXT_COEFF_550_LKT(22,12)=2445.900000 !rg=0.0559219 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,13,1:6)=(/ 2712.800000,2444.500000,1912.000000,1038.100000,301.380000,276.070000 /)
XPIZA_LKT(22,13,1:6)=(/ 0.971172,0.991082,0.999684,0.999623,0.997267,0.366218 /)
XCGA_LKT(22,13,1:6)=(/ 0.724343,0.714937,0.713960,0.675180,0.523617,0.292860 /)
XEXT_COEFF_550_LKT(22,13)=2246.300000 !rg=0.0559219 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,14,1:6)=(/ 2267.600000,2142.000000,1781.400000,1149.700000,404.170000,336.350000 /)
XPIZA_LKT(22,14,1:6)=(/ 0.966795,0.989728,0.999646,0.999640,0.997791,0.442251 /)
XCGA_LKT(22,14,1:6)=(/ 0.729097,0.718303,0.712093,0.700830,0.592690,0.397357 /)
XEXT_COEFF_550_LKT(22,14)=2012.000000 !rg=0.0559219 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,15,1:6)=(/ 1881.900000,1824.800000,1599.000000,1181.000000,507.500000,400.110000 /)
XPIZA_LKT(22,15,1:6)=(/ 0.960298,0.988004,0.999596,0.999633,0.998111,0.505982 /)
XCGA_LKT(22,15,1:6)=(/ 0.731717,0.720227,0.709560,0.713193,0.643130,0.476360 /)
XEXT_COEFF_550_LKT(22,15)=1766.900000 !rg=0.0559219 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,16,1:6)=(/ 1546.900000,1539.600000,1417.100000,1128.900000,598.620000,460.560000 /)
XPIZA_LKT(22,16,1:6)=(/ 0.952758,0.985473,0.999532,0.999600,0.998293,0.553149 /)
XCGA_LKT(22,16,1:6)=(/ 0.735933,0.721827,0.712677,0.713797,0.681237,0.545413 /)
XEXT_COEFF_550_LKT(22,16)=1501.900000 !rg=0.0559219 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,17,1:6)=(/ 1264.300000,1280.600000,1226.800000,1024.600000,657.620000,506.010000 /)
XPIZA_LKT(22,17,1:6)=(/ 0.945512,0.982679,0.999432,0.999545,0.998359,0.585731 /)
XCGA_LKT(22,17,1:6)=(/ 0.746460,0.723117,0.713947,0.709317,0.705027,0.599247 /)
XEXT_COEFF_550_LKT(22,17)=1270.700000 !rg=0.0559219 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,18,1:6)=(/ 1038.000000,1055.200000,1032.900000,915.740000,670.320000,529.440000 /)
XPIZA_LKT(22,18,1:6)=(/ 0.935978,0.980375,0.999350,0.999471,0.998314,0.604665 /)
XCGA_LKT(22,18,1:6)=(/ 0.753817,0.732887,0.716927,0.706863,0.715637,0.639983 /)
XEXT_COEFF_550_LKT(22,18)=1055.200000 !rg=0.0559219 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,19,1:6)=(/ 847.270000,872.440000,868.160000,805.860000,635.580000,526.320000 /)
XPIZA_LKT(22,19,1:6)=(/ 0.924883,0.976042,0.999200,0.999392,0.998153,0.611403 /)
XCGA_LKT(22,19,1:6)=(/ 0.760870,0.737247,0.719860,0.712370,0.714553,0.668957 /)
XEXT_COEFF_550_LKT(22,19)=874.990000 !rg=0.0559219 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,20,1:6)=(/ 689.790000,714.230000,726.610000,688.780000,574.080000,499.070000 /)
XPIZA_LKT(22,20,1:6)=(/ 0.913334,0.971804,0.998854,0.999295,0.997877,0.608252 /)
XCGA_LKT(22,20,1:6)=(/ 0.774433,0.741513,0.722460,0.715067,0.707310,0.689423 /)
XEXT_COEFF_550_LKT(22,20)=722.120000 !rg=0.0559219 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,1,1:6)=(/ 1676.600000,444.260000,82.049000,11.775000,1.740300,133.990000 /)
XPIZA_LKT(23,1,1:6)=(/ 0.969959,0.972904,0.996177,0.981892,0.692078,0.001884 /)
XCGA_LKT(23,1,1:6)=(/ 0.303757,0.133630,0.056897,0.021523,0.006983,0.002607 /)
XEXT_COEFF_550_LKT(23,1)=223.590000 !rg=0.0605833 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,2,1:6)=(/ 1887.100000,535.920000,102.780000,14.788000,2.055800,134.210000 /)
XPIZA_LKT(23,2,1:6)=(/ 0.972334,0.977033,0.996899,0.985475,0.738665,0.002372 /)
XCGA_LKT(23,2,1:6)=(/ 0.379540,0.169923,0.072123,0.027313,0.008870,0.003313 /)
XEXT_COEFF_550_LKT(23,2)=276.360000 !rg=0.0605833 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,3,1:6)=(/ 2335.300000,715.770000,150.670000,22.102000,2.825600,134.710000 /)
XPIZA_LKT(23,3,1:6)=(/ 0.976119,0.982063,0.997819,0.990133,0.808828,0.003559 /)
XCGA_LKT(23,3,1:6)=(/ 0.483793,0.253230,0.109077,0.041500,0.013520,0.005057 /)
XEXT_COEFF_550_LKT(23,3)=388.130000 !rg=0.0605833 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,4,1:6)=(/ 2938.700000,997.750000,237.530000,37.235000,4.454000,135.650000 /)
XPIZA_LKT(23,4,1:6)=(/ 0.979759,0.986260,0.998553,0.993996,0.877585,0.006058 /)
XCGA_LKT(23,4,1:6)=(/ 0.566243,0.374050,0.181167,0.070297,0.023057,0.008653 /)
XEXT_COEFF_550_LKT(23,4)=570.300000 !rg=0.0605833 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,5,1:6)=(/ 3575.600000,1384.400000,374.070000,65.394000,7.666000,137.270000 /)
XPIZA_LKT(23,5,1:6)=(/ 0.982373,0.989293,0.999024,0.996458,0.927831,0.010945 /)
XCGA_LKT(23,5,1:6)=(/ 0.627817,0.485097,0.286813,0.117067,0.038647,0.014563 /)
XEXT_COEFF_550_LKT(23,5)=835.260000 !rg=0.0605833 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,6,1:6)=(/ 4125.400000,1833.900000,566.260000,111.370000,13.517000,139.940000 /)
XPIZA_LKT(23,6,1:6)=(/ 0.983922,0.991274,0.999306,0.997823,0.958204,0.019735 /)
XCGA_LKT(23,6,1:6)=(/ 0.671840,0.568207,0.396270,0.182660,0.060413,0.022850 /)
XEXT_COEFF_550_LKT(23,6)=1173.600000 !rg=0.0605833 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,7,1:6)=(/ 4484.100000,2294.800000,814.790000,180.750000,23.819000,144.290000 /)
XPIZA_LKT(23,7,1:6)=(/ 0.984537,0.992523,0.999476,0.998578,0.975585,0.034973 /)
XCGA_LKT(23,7,1:6)=(/ 0.701147,0.628860,0.491397,0.271773,0.090487,0.034323 /)
XEXT_COEFF_550_LKT(23,7)=1559.700000 !rg=0.0605833 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,8,1:6)=(/ 4561.800000,2691.300000,1103.800000,281.210000,41.036000,151.280000 /)
XPIZA_LKT(23,8,1:6)=(/ 0.984222,0.993229,0.999581,0.999016,0.985287,0.060147 /)
XCGA_LKT(23,8,1:6)=(/ 0.717157,0.671243,0.567280,0.375387,0.132793,0.050440 /)
XEXT_COEFF_550_LKT(23,8)=1944.500000 !rg=0.0605833 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,9,1:6)=(/ 4357.500000,2952.300000,1407.700000,416.850000,68.328000,162.550000 /)
XPIZA_LKT(23,9,1:6)=(/ 0.982982,0.993498,0.999646,0.999282,0.990738,0.100039 /)
XCGA_LKT(23,9,1:6)=(/ 0.723117,0.698807,0.625373,0.464607,0.194013,0.073650 /)
XEXT_COEFF_550_LKT(23,9)=2273.700000 !rg=0.0605833 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,10,1:6)=(/ 3964.700000,3015.100000,1679.600000,580.110000,108.220000,180.330000 /)
XPIZA_LKT(23,10,1:6)=(/ 0.980871,0.993352,0.999683,0.999440,0.993801,0.158044 /)
XCGA_LKT(23,10,1:6)=(/ 0.723403,0.713037,0.667643,0.541577,0.281973,0.107610 /)
XEXT_COEFF_550_LKT(23,10)=2477.300000 !rg=0.0605833 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,11,1:6)=(/ 3505.700000,2880.100000,1873.400000,763.720000,166.600000,207.590000 /)
XPIZA_LKT(23,11,1:6)=(/ 0.978006,0.992802,0.999699,0.999541,0.995649,0.233711 /)
XCGA_LKT(23,11,1:6)=(/ 0.724330,0.717083,0.695983,0.604673,0.391363,0.158537 /)
XEXT_COEFF_550_LKT(23,11)=2513.200000 !rg=0.0605833 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,12,1:6)=(/ 3008.700000,2625.500000,1941.200000,944.870000,247.320000,246.660000 /)
XPIZA_LKT(23,12,1:6)=(/ 0.974191,0.991843,0.999696,0.999601,0.996818,0.317668 /)
XCGA_LKT(23,12,1:6)=(/ 0.725127,0.715870,0.710740,0.652940,0.481040,0.235153 /)
XEXT_COEFF_550_LKT(23,12)=2382.000000 !rg=0.0605833 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,13,1:6)=(/ 2522.900000,2328.300000,1873.600000,1092.300000,341.790000,299.490000 /)
XPIZA_LKT(23,13,1:6)=(/ 0.969773,0.990625,0.999671,0.999633,0.997508,0.398451 /)
XCGA_LKT(23,13,1:6)=(/ 0.727133,0.717197,0.713983,0.686873,0.552887,0.338210 /)
XEXT_COEFF_550_LKT(23,13)=2161.300000 !rg=0.0605833 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,14,1:6)=(/ 2104.600000,2018.100000,1711.600000,1173.200000,446.930000,362.850000 /)
XPIZA_LKT(23,14,1:6)=(/ 0.964200,0.989108,0.999629,0.999640,0.997944,0.470481 /)
XCGA_LKT(23,14,1:6)=(/ 0.729803,0.719803,0.711003,0.707020,0.614840,0.432807 /)
XEXT_COEFF_550_LKT(23,14)=1918.100000 !rg=0.0605833 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,15,1:6)=(/ 1732.700000,1704.900000,1529.200000,1169.600000,547.850000,425.330000 /)
XPIZA_LKT(23,15,1:6)=(/ 0.957199,0.986964,0.999572,0.999622,0.998203,0.527002 /)
XCGA_LKT(23,15,1:6)=(/ 0.732590,0.719730,0.711157,0.714240,0.660433,0.505070 /)
XEXT_COEFF_550_LKT(23,15)=1662.400000 !rg=0.0605833 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,16,1:6)=(/ 1423.000000,1428.200000,1345.100000,1091.500000,628.110000,481.690000 /)
XPIZA_LKT(23,16,1:6)=(/ 0.949078,0.984234,0.999497,0.999580,0.998335,0.568410 /)
XCGA_LKT(23,16,1:6)=(/ 0.740967,0.722307,0.713313,0.712073,0.692447,0.569650 /)
XEXT_COEFF_550_LKT(23,16)=1408.600000 !rg=0.0605833 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,17,1:6)=(/ 1160.400000,1181.300000,1142.400000,980.480000,668.700000,518.990000 /)
XPIZA_LKT(23,17,1:6)=(/ 0.941695,0.981623,0.999403,0.999516,0.998355,0.594976 /)
XCGA_LKT(23,17,1:6)=(/ 0.749113,0.728767,0.714483,0.707273,0.710737,0.617450 /)
XEXT_COEFF_550_LKT(23,17)=1176.900000 !rg=0.0605833 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,18,1:6)=(/ 950.880000,974.560000,964.140000,871.200000,660.720000,531.320000 /)
XPIZA_LKT(23,18,1:6)=(/ 0.930930,0.978486,0.999295,0.999444,0.998261,0.608600 /)
XCGA_LKT(23,18,1:6)=(/ 0.755660,0.733547,0.718457,0.710037,0.715957,0.652900 /)
XEXT_COEFF_550_LKT(23,18)=974.490000 !rg=0.0605833 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,19,1:6)=(/ 776.740000,801.010000,808.320000,758.710000,612.010000,517.500000 /)
XPIZA_LKT(23,19,1:6)=(/ 0.919452,0.973813,0.999101,0.999358,0.998048,0.610979 /)
XCGA_LKT(23,19,1:6)=(/ 0.767797,0.737470,0.720320,0.714410,0.711280,0.678000 /)
XEXT_COEFF_550_LKT(23,19)=810.470000 !rg=0.0605833 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,20,1:6)=(/ 632.970000,654.460000,666.990000,645.810000,546.040000,482.690000 /)
XPIZA_LKT(23,20,1:6)=(/ 0.907609,0.969453,0.999004,0.999216,0.997756,0.605024 /)
XCGA_LKT(23,20,1:6)=(/ 0.778833,0.748563,0.725947,0.714493,0.706167,0.696457 /)
XEXT_COEFF_550_LKT(23,20)=660.680000 !rg=0.0605833 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,1,1:6)=(/ 1952.400000,550.660000,103.990000,14.917000,2.069200,134.240000 /)
XPIZA_LKT(24,1,1:6)=(/ 0.973202,0.977463,0.996923,0.985580,0.740226,0.002391 /)
XCGA_LKT(24,1,1:6)=(/ 0.364860,0.157580,0.066703,0.025237,0.008193,0.003060 /)
XEXT_COEFF_550_LKT(24,1)=281.480000 !rg=0.0656333 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,2,1:6)=(/ 2196.200000,655.650000,129.930000,18.743000,2.470300,134.520000 /)
XPIZA_LKT(24,2,1:6)=(/ 0.975069,0.980639,0.997494,0.988426,0.781754,0.003010 /)
XCGA_LKT(24,2,1:6)=(/ 0.443177,0.200613,0.084560,0.032017,0.010407,0.003887 /)
XEXT_COEFF_550_LKT(24,2)=345.170000 !rg=0.0656333 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,3,1:6)=(/ 2692.600000,857.250000,188.580000,28.007000,3.448800,135.120000 /)
XPIZA_LKT(24,3,1:6)=(/ 0.978326,0.984498,0.998215,0.992118,0.842687,0.004514 /)
XCGA_LKT(24,3,1:6)=(/ 0.530520,0.295390,0.127760,0.048613,0.015857,0.005933 /)
XEXT_COEFF_550_LKT(24,3)=475.240000 !rg=0.0656333 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,4,1:6)=(/ 3300.100000,1171.700000,291.270000,46.973000,5.517400,136.240000 /)
XPIZA_LKT(24,4,1:6)=(/ 0.981296,0.987837,0.998788,0.995168,0.900598,0.007676 /)
XCGA_LKT(24,4,1:6)=(/ 0.597907,0.417270,0.210277,0.082210,0.027023,0.010153 /)
XEXT_COEFF_550_LKT(24,4)=682.390000 !rg=0.0656333 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,5,1:6)=(/ 3898.800000,1586.300000,448.640000,81.481000,9.588200,138.220000 /)
XPIZA_LKT(24,5,1:6)=(/ 0.983339,0.990292,0.999160,0.997103,0.941833,0.013838 /)
XCGA_LKT(24,5,1:6)=(/ 0.650023,0.519610,0.323133,0.136683,0.045263,0.017083 /)
XEXT_COEFF_550_LKT(24,5)=976.770000 !rg=0.0656333 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,6,1:6)=(/ 4366.300000,2045.300000,663.370000,136.220000,16.964000,141.460000 /)
XPIZA_LKT(24,6,1:6)=(/ 0.984432,0.991904,0.999386,0.998178,0.966329,0.024854 /)
XCGA_LKT(24,6,1:6)=(/ 0.687027,0.594197,0.431493,0.212777,0.070730,0.026800 /)
XEXT_COEFF_550_LKT(24,6)=1337.100000 !rg=0.0656333 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,7,1:6)=(/ 4604.700000,2489.400000,933.110000,217.680000,29.802000,146.760000 /)
XPIZA_LKT(24,7,1:6)=(/ 0.984637,0.992899,0.999526,0.998784,0.980201,0.043757 /)
XCGA_LKT(24,7,1:6)=(/ 0.710517,0.647753,0.522313,0.311507,0.105967,0.040263 /)
XEXT_COEFF_550_LKT(24,7)=1731.300000 !rg=0.0656333 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,8,1:6)=(/ 4549.900000,2836.600000,1231.200000,333.110000,50.780000,155.280000 /)
XPIZA_LKT(24,8,1:6)=(/ 0.983939,0.993410,0.999612,0.999142,0.987890,0.074444 /)
XCGA_LKT(24,8,1:6)=(/ 0.721927,0.684203,0.591200,0.411860,0.155680,0.059207 /)
XEXT_COEFF_550_LKT(24,8)=2100.400000 !rg=0.0656333 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,9,1:6)=(/ 4246.000000,3019.300000,1528.600000,480.280000,82.996000,168.960000 /)
XPIZA_LKT(24,9,1:6)=(/ 0.982279,0.993510,0.999664,0.999355,0.992201,0.121695 /)
XCGA_LKT(24,9,1:6)=(/ 0.724587,0.706383,0.643597,0.495273,0.227553,0.086573 /)
XEXT_COEFF_550_LKT(24,9)=2386.100000 !rg=0.0656333 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,10,1:6)=(/ 3804.400000,2995.400000,1774.000000,655.030000,129.450000,190.250000 /)
XPIZA_LKT(24,10,1:6)=(/ 0.979892,0.993197,0.999693,0.999487,0.994661,0.187381 /)
XCGA_LKT(24,10,1:6)=(/ 0.725047,0.716050,0.680193,0.569063,0.326893,0.126850 /)
XEXT_COEFF_550_LKT(24,10)=2522.200000 !rg=0.0656333 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,11,1:6)=(/ 3316.200000,2800.800000,1921.200000,839.750000,197.530000,222.080000 /)
XPIZA_LKT(24,11,1:6)=(/ 0.976646,0.992465,0.999701,0.999569,0.996200,0.267857 /)
XCGA_LKT(24,11,1:6)=(/ 0.725417,0.717317,0.703103,0.625117,0.431987,0.187520 /)
XEXT_COEFF_550_LKT(24,11)=2485.500000 !rg=0.0656333 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,12,1:6)=(/ 2809.100000,2514.900000,1932.700000,1010.600000,284.000000,266.550000 /)
XPIZA_LKT(24,12,1:6)=(/ 0.972617,0.991433,0.999689,0.999617,0.997141,0.351466 /)
XCGA_LKT(24,12,1:6)=(/ 0.726323,0.717043,0.713273,0.667837,0.509940,0.276530 /)
XEXT_COEFF_550_LKT(24,12)=2308.800000 !rg=0.0656333 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,13,1:6)=(/ 2346.900000,2213.000000,1817.900000,1135.300000,385.080000,324.910000 /)
XPIZA_LKT(24,13,1:6)=(/ 0.967478,0.990083,0.999657,0.999639,0.997712,0.429407 /)
XCGA_LKT(24,13,1:6)=(/ 0.727633,0.719077,0.713470,0.696327,0.580997,0.381040 /)
XEXT_COEFF_550_LKT(24,13)=2069.100000 !rg=0.0656333 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,14,1:6)=(/ 1939.700000,1888.700000,1642.100000,1183.400000,489.140000,388.770000 /)
XPIZA_LKT(24,14,1:6)=(/ 0.961318,0.988288,0.999609,0.999635,0.998064,0.495831 /)
XCGA_LKT(24,14,1:6)=(/ 0.730350,0.719613,0.711417,0.711297,0.634307,0.463363 /)
XEXT_COEFF_550_LKT(24,14)=1821.700000 !rg=0.0656333 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,15,1:6)=(/ 1602.100000,1586.900000,1460.900000,1145.700000,583.510000,449.970000 /)
XPIZA_LKT(24,15,1:6)=(/ 0.953575,0.986006,0.999546,0.999608,0.998269,0.545476 /)
XCGA_LKT(24,15,1:6)=(/ 0.736697,0.720830,0.711850,0.714167,0.674837,0.533360 /)
XEXT_COEFF_550_LKT(24,15)=1559.800000 !rg=0.0656333 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,16,1:6)=(/ 1302.600000,1325.700000,1261.700000,1048.300000,650.060000,499.320000 /)
XPIZA_LKT(24,16,1:6)=(/ 0.945992,0.982716,0.999464,0.999555,0.998355,0.580978 /)
XCGA_LKT(24,16,1:6)=(/ 0.743807,0.724853,0.714770,0.709513,0.701203,0.590333 /)
XEXT_COEFF_550_LKT(24,16)=1307.500000 !rg=0.0656333 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,17,1:6)=(/ 1063.900000,1089.300000,1063.900000,935.580000,671.590000,527.320000 /)
XPIZA_LKT(24,17,1:6)=(/ 0.937356,0.980589,0.999382,0.999490,0.998328,0.602208 /)
XCGA_LKT(24,17,1:6)=(/ 0.751443,0.730840,0.718127,0.708507,0.714057,0.633320 /)
XEXT_COEFF_550_LKT(24,17)=1087.800000 !rg=0.0656333 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,18,1:6)=(/ 873.720000,896.860000,900.600000,828.540000,644.460000,528.640000 /)
XPIZA_LKT(24,18,1:6)=(/ 0.926011,0.976205,0.999215,0.999413,0.998186,0.610854 /)
XCGA_LKT(24,18,1:6)=(/ 0.761773,0.734120,0.719537,0.712500,0.714783,0.664057 /)
XEXT_COEFF_550_LKT(24,18)=906.780000 !rg=0.0656333 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,19,1:6)=(/ 710.120000,736.920000,745.580000,710.690000,584.580000,505.060000 /)
XPIZA_LKT(24,19,1:6)=(/ 0.914324,0.971636,0.999096,0.999306,0.997941,0.609320 /)
XCGA_LKT(24,19,1:6)=(/ 0.772040,0.742400,0.720700,0.714017,0.709190,0.685937 /)
XEXT_COEFF_550_LKT(24,19)=740.970000 !rg=0.0656333 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,20,1:6)=(/ 580.250000,599.410000,611.980000,602.130000,521.260000,464.500000 /)
XPIZA_LKT(24,20,1:6)=(/ 0.900792,0.967813,0.998935,0.999153,0.997624,0.601328 /)
XCGA_LKT(24,20,1:6)=(/ 0.781930,0.752450,0.730700,0.715257,0.706813,0.703327 /)
XEXT_COEFF_550_LKT(24,20)=606.940000 !rg=0.0656333 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,1,1:6)=(/ 2249.300000,677.240000,131.640000,18.910000,2.487200,134.550000 /)
XPIZA_LKT(25,1,1:6)=(/ 0.975617,0.981091,0.997515,0.988508,0.783118,0.003035 /)
XCGA_LKT(25,1,1:6)=(/ 0.437000,0.186380,0.078223,0.029587,0.009610,0.003590 /)
XEXT_COEFF_550_LKT(25,1)=352.880000 !rg=0.0711042 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,2,1:6)=(/ 2550.700000,793.820000,163.840000,23.765000,2.997100,134.880000 /)
XPIZA_LKT(25,2,1:6)=(/ 0.977351,0.983485,0.997966,0.990767,0.819390,0.003819 /)
XCGA_LKT(25,2,1:6)=(/ 0.506353,0.237423,0.099190,0.037523,0.012207,0.004563 /)
XEXT_COEFF_550_LKT(25,2)=428.180000 !rg=0.0711042 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,3,1:6)=(/ 3077.200000,1018.000000,234.730000,35.474000,4.240900,135.610000 /)
XPIZA_LKT(25,3,1:6)=(/ 0.980195,0.986449,0.998529,0.993691,0.871419,0.005723 /)
XCGA_LKT(25,3,1:6)=(/ 0.571273,0.342430,0.149687,0.056940,0.018590,0.006963 /)
XEXT_COEFF_550_LKT(25,3)=576.700000 !rg=0.0711042 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,4,1:6)=(/ 3665.100000,1365.000000,354.580000,59.109000,6.867300,136.970000 /)
XPIZA_LKT(25,4,1:6)=(/ 0.982575,0.989130,0.998975,0.996095,0.919596,0.009717 /)
XCGA_LKT(25,4,1:6)=(/ 0.625837,0.459997,0.242980,0.096090,0.031663,0.011913 /)
XEXT_COEFF_550_LKT(25,4)=810.110000 !rg=0.0711042 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,5,1:6)=(/ 4203.100000,1801.300000,533.680000,100.970000,12.019000,139.360000 /)
XPIZA_LKT(25,5,1:6)=(/ 0.984114,0.991123,0.999270,0.997613,0.953165,0.017471 /)
XCGA_LKT(25,5,1:6)=(/ 0.669700,0.551910,0.360447,0.159403,0.052997,0.020040 /)
XEXT_COEFF_550_LKT(25,5)=1132.300000 !rg=0.0711042 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,6,1:6)=(/ 4567.700000,2258.100000,771.330000,165.650000,21.282000,143.310000 /)
XPIZA_LKT(25,6,1:6)=(/ 0.984782,0.992421,0.999452,0.998463,0.972824,0.031231 /)
XCGA_LKT(25,6,1:6)=(/ 0.700203,0.617777,0.466343,0.246603,0.082797,0.031440 /)
XEXT_COEFF_550_LKT(25,6)=1509.900000 !rg=0.0711042 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,7,1:6)=(/ 4673.500000,2670.300000,1057.200000,260.890000,37.152000,149.770000 /)
XPIZA_LKT(25,7,1:6)=(/ 0.984590,0.993194,0.999567,0.998953,0.983860,0.054537 /)
XCGA_LKT(25,7,1:6)=(/ 0.718253,0.664440,0.550340,0.350730,0.124100,0.047243 /)
XEXT_COEFF_550_LKT(25,7)=1901.200000 !rg=0.0711042 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,8,1:6)=(/ 4493.500000,2954.600000,1360.800000,389.890000,62.387000,160.140000 /)
XPIZA_LKT(25,8,1:6)=(/ 0.983493,0.993526,0.999637,0.999242,0.989945,0.091558 /)
XCGA_LKT(25,8,1:6)=(/ 0.725553,0.695010,0.613397,0.444810,0.182507,0.069523 /)
XEXT_COEFF_550_LKT(25,8)=2242.500000 !rg=0.0711042 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,9,1:6)=(/ 4102.400000,3053.200000,1644.100000,549.570000,100.050000,176.670000 /)
XPIZA_LKT(25,9,1:6)=(/ 0.981521,0.993454,0.999679,0.999417,0.993365,0.146602 /)
XCGA_LKT(25,9,1:6)=(/ 0.726490,0.712117,0.659700,0.526000,0.266120,0.101860 /)
XEXT_COEFF_550_LKT(25,9)=2473.600000 !rg=0.0711042 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,10,1:6)=(/ 3630.000000,2948.800000,1854.000000,730.270000,154.510000,201.920000 /)
XPIZA_LKT(25,10,1:6)=(/ 0.978719,0.992976,0.999698,0.999526,0.995376,0.219192 /)
XCGA_LKT(25,10,1:6)=(/ 0.725683,0.718097,0.690897,0.592590,0.372397,0.149740 /)
XEXT_COEFF_550_LKT(25,10)=2539.200000 !rg=0.0711042 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,11,1:6)=(/ 3118.000000,2702.000000,1947.100000,913.880000,231.390000,238.660000 /)
XPIZA_LKT(25,11,1:6)=(/ 0.975169,0.992115,0.999699,0.999593,0.996648,0.302354 /)
XCGA_LKT(25,11,1:6)=(/ 0.726443,0.717820,0.708607,0.643790,0.465717,0.221690 /)
XEXT_COEFF_550_LKT(25,11)=2436.700000 !rg=0.0711042 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,12,1:6)=(/ 2616.300000,2405.200000,1903.300000,1069.900000,323.140000,288.880000 /)
XPIZA_LKT(25,12,1:6)=(/ 0.970450,0.990951,0.999678,0.999629,0.997403,0.384207 /)
XCGA_LKT(25,12,1:6)=(/ 0.725970,0.718730,0.714253,0.680623,0.539240,0.320840 /)
XEXT_COEFF_550_LKT(25,12)=2223.200000 !rg=0.0711042 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,13,1:6)=(/ 2172.000000,2086.400000,1755.700000,1165.900000,428.200000,351.180000 /)
XPIZA_LKT(25,13,1:6)=(/ 0.965133,0.989397,0.999640,0.999641,0.997880,0.458541 /)
XCGA_LKT(25,13,1:6)=(/ 0.728830,0.719680,0.713180,0.703737,0.604653,0.418053 /)
XEXT_COEFF_550_LKT(25,13)=1977.200000 !rg=0.0711042 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,14,1:6)=(/ 1793.900000,1760.100000,1574.000000,1178.200000,530.820000,414.190000 /)
XPIZA_LKT(25,14,1:6)=(/ 0.957879,0.987485,0.999584,0.999628,0.998165,0.518023 /)
XCGA_LKT(25,14,1:6)=(/ 0.732927,0.720173,0.711023,0.713490,0.652507,0.492220 /)
XEXT_COEFF_550_LKT(25,14)=1721.300000 !rg=0.0711042 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,15,1:6)=(/ 1465.400000,1482.100000,1382.400000,1111.500000,615.630000,472.420000 /)
XPIZA_LKT(25,15,1:6)=(/ 0.950842,0.984609,0.999520,0.999590,0.998316,0.561749 /)
XCGA_LKT(25,15,1:6)=(/ 0.738517,0.723147,0.715110,0.712687,0.687000,0.558853 /)
XEXT_COEFF_550_LKT(25,15)=1452.100000 !rg=0.0711042 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,16,1:6)=(/ 1192.900000,1218.400000,1179.200000,1003.000000,665.070000,513.900000 /)
XPIZA_LKT(25,16,1:6)=(/ 0.942928,0.981661,0.999432,0.999531,0.998357,0.591098 /)
XCGA_LKT(25,16,1:6)=(/ 0.747313,0.726320,0.716877,0.709180,0.708030,0.609183 /)
XEXT_COEFF_550_LKT(25,16)=1213.000000 !rg=0.0711042 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,17,1:6)=(/ 980.050000,1001.800000,997.140000,894.700000,665.690000,531.060000 /)
XPIZA_LKT(25,17,1:6)=(/ 0.932174,0.979208,0.999299,0.999459,0.998284,0.607042 /)
XCGA_LKT(25,17,1:6)=(/ 0.756430,0.732180,0.717517,0.710083,0.715507,0.646933 /)
XEXT_COEFF_550_LKT(25,17)=1009.600000 !rg=0.0711042 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,18,1:6)=(/ 797.010000,827.330000,831.440000,781.010000,622.140000,521.550000 /)
XPIZA_LKT(25,18,1:6)=(/ 0.921363,0.974622,0.999162,0.999372,0.998098,0.611208 /)
XCGA_LKT(25,18,1:6)=(/ 0.764760,0.738977,0.718363,0.713883,0.712877,0.673637 /)
XEXT_COEFF_550_LKT(25,18)=832.250000 !rg=0.0711042 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,19,1:6)=(/ 649.250000,674.060000,686.130000,662.860000,558.910000,489.790000 /)
XPIZA_LKT(25,19,1:6)=(/ 0.909352,0.969462,0.999067,0.999255,0.997808,0.606537 /)
XCGA_LKT(25,19,1:6)=(/ 0.776500,0.745580,0.728687,0.714540,0.707147,0.693233 /)
XEXT_COEFF_550_LKT(25,19)=678.670000 !rg=0.0711042 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,20,1:6)=(/ 533.520000,550.020000,568.070000,564.710000,497.180000,445.250000 /)
XPIZA_LKT(25,20,1:6)=(/ 0.893520,0.965398,0.998816,0.999047,0.997479,0.597506 /)
XCGA_LKT(25,20,1:6)=(/ 0.787943,0.754383,0.730350,0.717893,0.707490,0.710457 /)
XEXT_COEFF_550_LKT(25,20)=560.890000 !rg=0.0711042 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,1,1:6)=(/ 2591.500000,824.030000,166.340000,23.983000,3.018600,134.920000 /)
XPIZA_LKT(26,1,1:6)=(/ 0.977509,0.983954,0.997985,0.990833,0.820560,0.003850 /)
XCGA_LKT(26,1,1:6)=(/ 0.515010,0.221317,0.091783,0.034683,0.011273,0.004213 /)
XEXT_COEFF_550_LKT(26,1)=439.770000 !rg=0.0770312 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,2,1:6)=(/ 2958.300000,949.810000,205.900000,30.138000,3.666900,135.320000 /)
XPIZA_LKT(26,2,1:6)=(/ 0.979364,0.985720,0.998339,0.992624,0.851683,0.004842 /)
XCGA_LKT(26,2,1:6)=(/ 0.561283,0.281530,0.116443,0.043980,0.014317,0.005353 /)
XEXT_COEFF_550_LKT(26,2)=526.460000 !rg=0.0770312 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,3,1:6)=(/ 3478.600000,1199.500000,290.200000,44.894000,5.247600,136.190000 /)
XPIZA_LKT(26,3,1:6)=(/ 0.981780,0.988023,0.998777,0.994937,0.895472,0.007251 /)
XCGA_LKT(26,3,1:6)=(/ 0.605487,0.393003,0.175423,0.066683,0.021797,0.008170 /)
XEXT_COEFF_550_LKT(26,3)=693.350000 !rg=0.0770312 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,4,1:6)=(/ 4022.900000,1577.100000,428.290000,74.133000,8.579900,137.840000 /)
XPIZA_LKT(26,4,1:6)=(/ 0.983629,0.990197,0.999124,0.996827,0.935140,0.012289 /)
XCGA_LKT(26,4,1:6)=(/ 0.650347,0.500940,0.279153,0.112257,0.037093,0.013973 /)
XEXT_COEFF_550_LKT(26,4)=954.130000 !rg=0.0770312 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,5,1:6)=(/ 4477.900000,2025.300000,629.900000,124.380000,15.085000,140.750000 /)
XPIZA_LKT(26,5,1:6)=(/ 0.984712,0.991811,0.999359,0.998018,0.962287,0.022022 /)
XCGA_LKT(26,5,1:6)=(/ 0.686967,0.581567,0.398517,0.185543,0.062037,0.023507 /)
XEXT_COEFF_550_LKT(26,5)=1300.800000 !rg=0.0770312 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,6,1:6)=(/ 4721.700000,2467.400000,889.030000,200.440000,26.657000,145.560000 /)
XPIZA_LKT(26,6,1:6)=(/ 0.984975,0.992845,0.999508,0.998693,0.977999,0.039133 /)
XCGA_LKT(26,6,1:6)=(/ 0.711373,0.639057,0.499520,0.283053,0.096897,0.036880 /)
XEXT_COEFF_550_LKT(26,6)=1688.200000 !rg=0.0770312 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,7,1:6)=(/ 4688.600000,2833.600000,1186.900000,310.130000,46.083000,153.440000 /)
XPIZA_LKT(26,7,1:6)=(/ 0.984392,0.993416,0.999601,0.999090,0.986755,0.067651 /)
XCGA_LKT(26,7,1:6)=(/ 0.724333,0.679153,0.576290,0.387337,0.145333,0.055443 /)
XEXT_COEFF_550_LKT(26,7)=2067.100000 !rg=0.0770312 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,8,1:6)=(/ 4394.600000,3042.500000,1486.900000,451.550000,76.050000,166.030000 /)
XPIZA_LKT(26,8,1:6)=(/ 0.982892,0.993576,0.999658,0.999324,0.991568,0.111752 /)
XCGA_LKT(26,8,1:6)=(/ 0.727917,0.704093,0.633167,0.476547,0.213807,0.081690 /)
XEXT_COEFF_550_LKT(26,8)=2368.900000 !rg=0.0770312 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,9,1:6)=(/ 3944.100000,3051.900000,1747.200000,623.550000,120.000000,185.880000 /)
XPIZA_LKT(26,9,1:6)=(/ 0.980545,0.993341,0.999690,0.999468,0.994308,0.174615 /)
XCGA_LKT(26,9,1:6)=(/ 0.727660,0.716467,0.673593,0.554613,0.308887,0.119983 /)
XEXT_COEFF_550_LKT(26,9)=2535.100000 !rg=0.0770312 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,10,1:6)=(/ 3436.100000,2877.700000,1913.400000,806.700000,183.680000,215.490000 /)
XPIZA_LKT(26,10,1:6)=(/ 0.977467,0.992691,0.999701,0.999557,0.995974,0.252682 /)
XCGA_LKT(26,10,1:6)=(/ 0.727473,0.718987,0.699307,0.613987,0.413847,0.176943 /)
XEXT_COEFF_550_LKT(26,10)=2526.900000 !rg=0.0770312 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,11,1:6)=(/ 2917.200000,2597.600000,1951.300000,982.470000,267.050000,257.480000 /)
XPIZA_LKT(26,11,1:6)=(/ 0.973229,0.991720,0.999693,0.999610,0.997004,0.336409 /)
XCGA_LKT(26,11,1:6)=(/ 0.726003,0.719167,0.712127,0.659790,0.495460,0.260963 /)
XEXT_COEFF_550_LKT(26,11)=2369.800000 !rg=0.0770312 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,12,1:6)=(/ 2426.400000,2287.500000,1858.200000,1118.100000,365.500000,313.400000 /)
XPIZA_LKT(26,12,1:6)=(/ 0.968316,0.990363,0.999666,0.999637,0.997624,0.415786 /)
XCGA_LKT(26,12,1:6)=(/ 0.727007,0.719237,0.714710,0.691083,0.568210,0.364010 /)
XEXT_COEFF_550_LKT(26,12)=2132.900000 !rg=0.0770312 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,13,1:6)=(/ 2011.400000,1954.400000,1688.900000,1182.700000,470.490000,377.210000 /)
XPIZA_LKT(26,13,1:6)=(/ 0.961808,0.988682,0.999619,0.999638,0.998013,0.485026 /)
XCGA_LKT(26,13,1:6)=(/ 0.730407,0.720290,0.712090,0.708840,0.624890,0.449700 /)
XEXT_COEFF_550_LKT(26,13)=1882.600000 !rg=0.0770312 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,14,1:6)=(/ 1646.100000,1647.200000,1500.400000,1160.500000,568.230000,439.280000 /)
XPIZA_LKT(26,14,1:6)=(/ 0.955236,0.986081,0.999561,0.999615,0.998242,0.537501 /)
XCGA_LKT(26,14,1:6)=(/ 0.734223,0.720940,0.713597,0.714100,0.668000,0.520843 /)
XEXT_COEFF_550_LKT(26,14)=1608.500000 !rg=0.0770312 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,15,1:6)=(/ 1339.900000,1366.800000,1304.400000,1071.300000,640.540000,491.510000 /)
XPIZA_LKT(26,15,1:6)=(/ 0.948245,0.983812,0.999485,0.999569,0.998346,0.575446 /)
XCGA_LKT(26,15,1:6)=(/ 0.743403,0.723487,0.716103,0.711660,0.696707,0.580557 /)
XEXT_COEFF_550_LKT(26,15)=1349.900000 !rg=0.0770312 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,16,1:6)=(/ 1098.300000,1122.000000,1099.300000,961.290000,671.330000,524.330000 /)
XPIZA_LKT(26,16,1:6)=(/ 0.938233,0.980996,0.999394,0.999502,0.998340,0.599263 /)
XCGA_LKT(26,16,1:6)=(/ 0.751010,0.729777,0.716583,0.708423,0.712183,0.626010 /)
XEXT_COEFF_550_LKT(26,16)=1127.300000 !rg=0.0770312 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,17,1:6)=(/ 894.590000,926.580000,925.790000,850.140000,652.210000,530.450000 /)
XPIZA_LKT(26,17,1:6)=(/ 0.927747,0.976845,0.999260,0.999423,0.998221,0.610065 /)
XCGA_LKT(26,17,1:6)=(/ 0.757980,0.734663,0.718257,0.711767,0.715267,0.658857 /)
XEXT_COEFF_550_LKT(26,17)=931.350000 !rg=0.0770312 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,18,1:6)=(/ 727.030000,756.200000,771.090000,732.260000,597.910000,510.760000 /)
XPIZA_LKT(26,18,1:6)=(/ 0.917053,0.973023,0.999011,0.999330,0.997981,0.610172 /)
XCGA_LKT(26,18,1:6)=(/ 0.771687,0.739780,0.721750,0.714407,0.709750,0.682110 /)
XEXT_COEFF_550_LKT(26,18)=762.200000 !rg=0.0770312 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,19,1:6)=(/ 597.010000,615.530000,633.060000,622.790000,532.700000,472.420000 /)
XPIZA_LKT(26,19,1:6)=(/ 0.902565,0.968527,0.998849,0.999158,0.997678,0.603170 /)
XCGA_LKT(26,19,1:6)=(/ 0.781920,0.750390,0.728530,0.715933,0.706303,0.700270 /)
XEXT_COEFF_550_LKT(26,19)=628.060000 !rg=0.0770312 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,20,1:6)=(/ 487.560000,507.170000,522.120000,521.290000,473.180000,425.420000 /)
XPIZA_LKT(26,20,1:6)=(/ 0.887425,0.961953,0.998661,0.998996,0.997324,0.593823 /)
XCGA_LKT(26,20,1:6)=(/ 0.790753,0.758507,0.730670,0.719880,0.710813,0.717983 /)
XEXT_COEFF_550_LKT(26,20)=513.340000 !rg=0.0770312 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,1,1:6)=(/ 3024.300000,988.930000,209.600000,30.423000,3.694200,135.370000 /)
XPIZA_LKT(27,1,1:6)=(/ 0.979280,0.986186,0.998358,0.992676,0.852669,0.004881 /)
XCGA_LKT(27,1,1:6)=(/ 0.584903,0.264030,0.107787,0.040653,0.013223,0.004943 /)
XEXT_COEFF_550_LKT(27,1)=543.560000 !rg=0.0834522 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,2,1:6)=(/ 3409.200000,1123.200000,257.500000,38.213000,4.518300,135.840000 /)
XPIZA_LKT(27,2,1:6)=(/ 0.981196,0.987464,0.998635,0.994096,0.878971,0.006137 /)
XCGA_LKT(27,2,1:6)=(/ 0.602403,0.333740,0.136857,0.051540,0.016787,0.006280 /)
XEXT_COEFF_550_LKT(27,2)=640.330000 !rg=0.0834522 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,3,1:6)=(/ 3882.100000,1403.500000,355.940000,56.735000,6.526700,136.900000 /)
XPIZA_LKT(27,3,1:6)=(/ 0.983115,0.989305,0.998973,0.995923,0.915383,0.009179 /)
XCGA_LKT(27,3,1:6)=(/ 0.633960,0.444760,0.205580,0.078090,0.025550,0.009583 /)
XEXT_COEFF_550_LKT(27,3)=826.090000 !rg=0.0834522 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,4,1:6)=(/ 4361.400000,1805.800000,513.230000,92.581000,10.750000,138.900000 /)
XPIZA_LKT(27,4,1:6)=(/ 0.984483,0.991083,0.999244,0.997405,0.947766,0.015523 /)
XCGA_LKT(27,4,1:6)=(/ 0.671853,0.539013,0.318460,0.131057,0.043443,0.016393 /)
XEXT_COEFF_550_LKT(27,4)=1114.600000 !rg=0.0834522 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,5,1:6)=(/ 4712.800000,2254.000000,737.730000,152.330000,18.938000,142.440000 /)
XPIZA_LKT(27,5,1:6)=(/ 0.985146,0.992381,0.999432,0.998340,0.969598,0.027702 /)
XCGA_LKT(27,5,1:6)=(/ 0.702000,0.608480,0.436583,0.215167,0.072593,0.027573 /)
XEXT_COEFF_550_LKT(27,5)=1480.800000 !rg=0.0834522 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,6,1:6)=(/ 4821.200000,2666.600000,1014.600000,241.170000,33.297000,148.310000 /)
XPIZA_LKT(27,6,1:6)=(/ 0.985020,0.993184,0.999553,0.998881,0.982110,0.048863 /)
XCGA_LKT(27,6,1:6)=(/ 0.720703,0.658017,0.530340,0.320240,0.113373,0.043267 /)
XEXT_COEFF_550_LKT(27,6)=1868.100000 !rg=0.0834522 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,7,1:6)=(/ 4651.700000,2972.000000,1320.000000,364.750000,56.797000,157.900000 /)
XPIZA_LKT(27,7,1:6)=(/ 0.984042,0.993570,0.999629,0.999201,0.989043,0.083428 /)
XCGA_LKT(27,7,1:6)=(/ 0.729130,0.691727,0.600370,0.421483,0.170150,0.065087 /)
XEXT_COEFF_550_LKT(27,7)=2222.200000 !rg=0.0834522 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,8,1:6)=(/ 4258.500000,3097.300000,1609.100000,519.190000,92.025000,173.140000 /)
XPIZA_LKT(27,8,1:6)=(/ 0.982191,0.993561,0.999674,0.999391,0.992860,0.135157 /)
XCGA_LKT(27,8,1:6)=(/ 0.730200,0.711303,0.650803,0.508163,0.249833,0.096060 /)
XEXT_COEFF_550_LKT(27,8)=2472.600000 !rg=0.0834522 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,9,1:6)=(/ 3758.400000,3021.600000,1837.200000,698.950000,143.540000,196.760000 /)
XPIZA_LKT(27,9,1:6)=(/ 0.979490,0.993152,0.999697,0.999510,0.995087,0.205308 /)
XCGA_LKT(27,9,1:6)=(/ 0.728753,0.719360,0.685587,0.579513,0.353117,0.141493 /)
XEXT_COEFF_550_LKT(27,9)=2568.400000 !rg=0.0834522 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,10,1:6)=(/ 3236.800000,2787.000000,1952.400000,882.790000,216.080000,231.100000 /)
XPIZA_LKT(27,10,1:6)=(/ 0.975807,0.992379,0.999701,0.999583,0.996463,0.286929 /)
XCGA_LKT(27,10,1:6)=(/ 0.727540,0.720450,0.706097,0.633777,0.449007,0.209043 /)
XEXT_COEFF_550_LKT(27,10)=2491.800000 !rg=0.0834522 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,11,1:6)=(/ 2707.800000,2488.000000,1933.100000,1046.100000,305.010000,278.710000 /)
XPIZA_LKT(27,11,1:6)=(/ 0.971220,0.991221,0.999685,0.999624,0.997289,0.369604 /)
XCGA_LKT(27,11,1:6)=(/ 0.726293,0.719427,0.714450,0.673680,0.524943,0.303787 /)
XEXT_COEFF_550_LKT(27,11)=2288.600000 !rg=0.0834522 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,12,1:6)=(/ 2254.800000,2160.500000,1802.000000,1155.700000,408.690000,339.210000 /)
XPIZA_LKT(27,12,1:6)=(/ 0.965442,0.989742,0.999649,0.999640,0.997808,0.445760 /)
XCGA_LKT(27,12,1:6)=(/ 0.728723,0.720763,0.713840,0.699680,0.593360,0.402263 /)
XEXT_COEFF_550_LKT(27,12)=2043.800000 !rg=0.0834522 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,13,1:6)=(/ 1844.800000,1828.500000,1617.000000,1185.100000,512.980000,402.790000 /)
XPIZA_LKT(27,13,1:6)=(/ 0.959180,0.987587,0.999599,0.999632,0.998122,0.508417 /)
XCGA_LKT(27,13,1:6)=(/ 0.730640,0.720073,0.713357,0.712213,0.643870,0.478880 /)
XEXT_COEFF_550_LKT(27,13)=1773.500000 !rg=0.0834522 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,14,1:6)=(/ 1506.300000,1527.500000,1427.600000,1131.800000,602.390000,462.860000 /)
XPIZA_LKT(27,14,1:6)=(/ 0.952900,0.985222,0.999536,0.999599,0.998295,0.554758 /)
XCGA_LKT(27,14,1:6)=(/ 0.739583,0.721263,0.715203,0.713770,0.681073,0.547437 /)
XEXT_COEFF_550_LKT(27,14)=1501.000000 !rg=0.0834522 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,15,1:6)=(/ 1234.200000,1256.700000,1219.100000,1029.800000,659.270000,507.600000 /)
XPIZA_LKT(27,15,1:6)=(/ 0.943790,0.983119,0.999451,0.999542,0.998355,0.586549 /)
XCGA_LKT(27,15,1:6)=(/ 0.745907,0.727733,0.716317,0.709527,0.704507,0.600040 /)
XEXT_COEFF_550_LKT(27,15)=1256.400000 !rg=0.0834522 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,16,1:6)=(/ 1004.500000,1038.100000,1024.900000,916.130000,669.250000,530.030000 /)
XPIZA_LKT(27,16,1:6)=(/ 0.933472,0.979269,0.999349,0.999476,0.998305,0.605094 /)
XCGA_LKT(27,16,1:6)=(/ 0.752513,0.731200,0.717633,0.710490,0.714647,0.640410 /)
XEXT_COEFF_550_LKT(27,16)=1038.800000 !rg=0.0834522 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,17,1:6)=(/ 817.080000,849.520000,860.780000,804.020000,633.200000,525.210000 /)
XPIZA_LKT(27,17,1:6)=(/ 0.923019,0.975477,0.999176,0.999392,0.998134,0.611273 /)
XCGA_LKT(27,17,1:6)=(/ 0.765817,0.734967,0.720043,0.714023,0.713330,0.669050 /)
XEXT_COEFF_550_LKT(27,17)=860.960000 !rg=0.0834522 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,18,1:6)=(/ 668.840000,690.760000,710.070000,686.200000,570.580000,496.760000 /)
XPIZA_LKT(27,18,1:6)=(/ 0.910742,0.971738,0.998878,0.999268,0.997869,0.607942 /)
XCGA_LKT(27,18,1:6)=(/ 0.775777,0.747317,0.724453,0.714447,0.707803,0.689723 /)
XEXT_COEFF_550_LKT(27,18)=701.500000 !rg=0.0834522 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,19,1:6)=(/ 546.430000,567.670000,582.140000,579.540000,508.560000,453.710000 /)
XPIZA_LKT(27,19,1:6)=(/ 0.895368,0.965876,0.998871,0.999125,0.997539,0.599514 /)
XCGA_LKT(27,19,1:6)=(/ 0.784433,0.753867,0.729683,0.716670,0.708133,0.707440 /)
XEXT_COEFF_550_LKT(27,19)=576.150000 !rg=0.0834522 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,20,1:6)=(/ 445.830000,463.200000,480.520000,481.140000,449.070000,405.310000 /)
XPIZA_LKT(27,20,1:6)=(/ 0.881427,0.959654,0.998619,0.998915,0.997123,0.590387 /)
XCGA_LKT(27,20,1:6)=(/ 0.798490,0.759647,0.735690,0.722893,0.711933,0.725900 /)
XEXT_COEFF_550_LKT(27,20)=471.950000 !rg=0.0834522 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,1,1:6)=(/ 3571.000000,1168.200000,263.070000,38.593000,4.553200,135.910000 /)
XPIZA_LKT(28,1,1:6)=(/ 0.981242,0.987894,0.998654,0.994138,0.879789,0.006185 /)
XCGA_LKT(28,1,1:6)=(/ 0.628990,0.316443,0.126750,0.047650,0.015507,0.005797 /)
XEXT_COEFF_550_LKT(28,1)=664.500000 !rg=0.0904084 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,2,1:6)=(/ 3870.000000,1315.800000,319.920000,48.428000,5.600600,136.480000 /)
XPIZA_LKT(28,2,1:6)=(/ 0.982832,0.988832,0.998869,0.995262,0.901736,0.007772 /)
XCGA_LKT(28,2,1:6)=(/ 0.630117,0.393533,0.161103,0.060400,0.019687,0.007370 /)
XEXT_COEFF_550_LKT(28,2)=769.330000 !rg=0.0904084 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,3,1:6)=(/ 4272.500000,1630.800000,432.680000,71.546000,8.151400,137.760000 /)
XPIZA_LKT(28,3,1:6)=(/ 0.984217,0.990365,0.999128,0.996703,0.931710,0.011609 /)
XCGA_LKT(28,3,1:6)=(/ 0.658213,0.494733,0.240713,0.091440,0.029947,0.011243 /)
XEXT_COEFF_550_LKT(28,3)=976.190000 !rg=0.0904084 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,4,1:6)=(/ 4668.700000,2047.500000,610.180000,115.040000,13.496000,140.190000 /)
XPIZA_LKT(28,4,1:6)=(/ 0.985151,0.991819,0.999340,0.997864,0.957962,0.019580 /)
XCGA_LKT(28,4,1:6)=(/ 0.690743,0.573517,0.360213,0.152860,0.050863,0.019230 /)
XEXT_COEFF_550_LKT(28,4)=1291.100000 !rg=0.0904084 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,5,1:6)=(/ 4897.900000,2481.700000,856.600000,185.470000,23.759000,144.490000 /)
XPIZA_LKT(28,5,1:6)=(/ 0.985419,0.992851,0.999492,0.998599,0.975436,0.034759 /)
XCGA_LKT(28,5,1:6)=(/ 0.714873,0.632643,0.473513,0.247830,0.084913,0.032343 /)
XEXT_COEFF_550_LKT(28,5)=1669.100000 !rg=0.0904084 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,6,1:6)=(/ 4863.400000,2850.100000,1147.600000,287.900000,41.417000,151.650000 /)
XPIZA_LKT(28,6,1:6)=(/ 0.984900,0.993448,0.999590,0.999032,0.985367,0.060748 /)
XCGA_LKT(28,6,1:6)=(/ 0.728177,0.674807,0.559177,0.356653,0.132607,0.050763 /)
XEXT_COEFF_550_LKT(28,6)=2046.100000 !rg=0.0904084 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,7,1:6)=(/ 4566.600000,3081.400000,1452.100000,424.840000,69.506000,163.320000 /)
XPIZA_LKT(28,7,1:6)=(/ 0.983509,0.993660,0.999652,0.999290,0.990853,0.102157 /)
XCGA_LKT(28,7,1:6)=(/ 0.732330,0.702463,0.622053,0.454690,0.199030,0.076443 /)
XEXT_COEFF_550_LKT(28,7)=2363.400000 !rg=0.0904084 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,8,1:6)=(/ 4098.700000,3116.400000,1721.300000,592.040000,110.750000,181.650000 /)
XPIZA_LKT(28,8,1:6)=(/ 0.981250,0.993488,0.999687,0.999447,0.993903,0.161716 /)
XCGA_LKT(28,8,1:6)=(/ 0.731260,0.717037,0.666257,0.537940,0.290040,0.113057 /)
XEXT_COEFF_550_LKT(28,8)=2551.700000 !rg=0.0904084 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,9,1:6)=(/ 3560.600000,2960.600000,1908.400000,775.960000,171.030000,209.480000 /)
XPIZA_LKT(28,9,1:6)=(/ 0.978229,0.992924,0.999702,0.999544,0.995736,0.238002 /)
XCGA_LKT(28,9,1:6)=(/ 0.729790,0.721533,0.695400,0.602090,0.394643,0.167020 /)
XEXT_COEFF_550_LKT(28,9)=2572.600000 !rg=0.0904084 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,10,1:6)=(/ 3017.000000,2688.400000,1968.900000,954.310000,250.640000,248.890000 /)
XPIZA_LKT(28,10,1:6)=(/ 0.974118,0.991968,0.999698,0.999603,0.996852,0.321103 /)
XCGA_LKT(28,10,1:6)=(/ 0.727567,0.720940,0.710827,0.650973,0.479813,0.246137 /)
XEXT_COEFF_550_LKT(28,10)=2432.400000 !rg=0.0904084 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,11,1:6)=(/ 2516.900000,2365.000000,1898.400000,1099.700000,346.280000,302.260000 /)
XPIZA_LKT(28,11,1:6)=(/ 0.968749,0.990734,0.999673,0.999633,0.997528,0.401768 /)
XCGA_LKT(28,11,1:6)=(/ 0.727773,0.721590,0.715210,0.685257,0.554537,0.346653 /)
XEXT_COEFF_550_LKT(28,11)=2204.800000 !rg=0.0904084 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,12,1:6)=(/ 2072.000000,2030.800000,1735.200000,1179.400000,451.130000,365.180000 /)
XPIZA_LKT(28,12,1:6)=(/ 0.963111,0.988892,0.999632,0.999640,0.997954,0.473329 /)
XCGA_LKT(28,12,1:6)=(/ 0.728990,0.720617,0.714427,0.705853,0.614600,0.435050 /)
XEXT_COEFF_550_LKT(28,12)=1940.800000 !rg=0.0904084 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,13,1:6)=(/ 1691.900000,1698.100000,1548.000000,1173.900000,552.210000,428.210000 /)
XPIZA_LKT(28,13,1:6)=(/ 0.957011,0.986617,0.999575,0.999622,0.998209,0.528967 /)
XCGA_LKT(28,13,1:6)=(/ 0.735733,0.719847,0.714417,0.713920,0.660490,0.507703 /)
XEXT_COEFF_550_LKT(28,13)=1665.600000 !rg=0.0904084 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,14,1:6)=(/ 1382.200000,1410.100000,1345.300000,1095.700000,630.460000,483.390000 /)
XPIZA_LKT(28,14,1:6)=(/ 0.948977,0.984606,0.999503,0.999578,0.998335,0.569572 /)
XCGA_LKT(28,14,1:6)=(/ 0.741297,0.726083,0.716393,0.712040,0.691840,0.570290 /)
XEXT_COEFF_550_LKT(28,14)=1395.200000 !rg=0.0904084 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,15,1:6)=(/ 1128.500000,1161.900000,1136.900000,984.550000,669.190000,520.100000 /)
XPIZA_LKT(28,15,1:6)=(/ 0.939140,0.981549,0.999417,0.999518,0.998348,0.595640 /)
XCGA_LKT(28,15,1:6)=(/ 0.747443,0.728240,0.717150,0.709967,0.709677,0.617777 /)
XEXT_COEFF_550_LKT(28,15)=1162.700000 !rg=0.0904084 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,16,1:6)=(/ 920.300000,951.720000,959.610000,873.630000,659.630000,531.510000 /)
XPIZA_LKT(28,16,1:6)=(/ 0.928335,0.977363,0.999257,0.999445,0.998248,0.608908 /)
XCGA_LKT(28,16,1:6)=(/ 0.759043,0.731173,0.718517,0.712850,0.715073,0.653093 /)
XEXT_COEFF_550_LKT(28,16)=963.040000 !rg=0.0904084 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,17,1:6)=(/ 747.790000,778.390000,791.350000,756.610000,608.780000,516.140000 /)
XPIZA_LKT(28,17,1:6)=(/ 0.918104,0.973231,0.999153,0.999340,0.998040,0.610875 /)
XCGA_LKT(28,17,1:6)=(/ 0.770167,0.741207,0.721053,0.713737,0.711457,0.678093 /)
XEXT_COEFF_550_LKT(28,17)=786.740000 !rg=0.0904084 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,18,1:6)=(/ 612.580000,635.260000,652.830000,639.780000,545.210000,480.400000 /)
XPIZA_LKT(28,18,1:6)=(/ 0.903848,0.969278,0.998921,0.999197,0.997745,0.604962 /)
XCGA_LKT(28,18,1:6)=(/ 0.778487,0.749243,0.726840,0.713617,0.707720,0.697010 /)
XEXT_COEFF_550_LKT(28,18)=643.700000 !rg=0.0904084 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,19,1:6)=(/ 500.790000,519.700000,539.760000,536.080000,485.600000,434.190000 /)
XPIZA_LKT(28,19,1:6)=(/ 0.888597,0.962947,0.998735,0.999108,0.997388,0.595900 /)
XCGA_LKT(28,19,1:6)=(/ 0.791707,0.754853,0.732023,0.722927,0.709307,0.714947 /)
XEXT_COEFF_550_LKT(28,19)=529.250000 !rg=0.0904084 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,20,1:6)=(/ 409.000000,423.840000,439.310000,447.890000,421.060000,385.020000 /)
XPIZA_LKT(28,20,1:6)=(/ 0.875785,0.956175,0.998544,0.998866,0.996936,0.587129 /)
XCGA_LKT(28,20,1:6)=(/ 0.803397,0.767100,0.737967,0.723493,0.713617,0.734023 /)
XEXT_COEFF_550_LKT(28,20)=431.210000 !rg=0.0904084 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,1,1:6)=(/ 4145.100000,1359.100000,328.300000,48.942000,5.645100,136.550000 /)
XPIZA_LKT(29,1,1:6)=(/ 0.983304,0.989170,0.998889,0.995297,0.902408,0.007832 /)
XCGA_LKT(29,1,1:6)=(/ 0.643397,0.380110,0.149337,0.055853,0.018187,0.006803 /)
XEXT_COEFF_550_LKT(29,1)=801.130000 !rg=0.0979445 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,2,1:6)=(/ 4300.700000,1533.300000,394.090000,61.316000,6.976200,137.240000 /)
XPIZA_LKT(29,2,1:6)=(/ 0.984200,0.989929,0.999053,0.996186,0.920529,0.009834 /)
XCGA_LKT(29,2,1:6)=(/ 0.650240,0.457680,0.190043,0.070797,0.023080,0.008643 /)
XEXT_COEFF_550_LKT(29,2)=913.000000 !rg=0.0979445 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,3,1:6)=(/ 4635.500000,1880.300000,521.090000,89.954000,10.214000,138.790000 /)
XPIZA_LKT(29,3,1:6)=(/ 0.985101,0.991253,0.999251,0.997321,0.945002,0.014665 /)
XCGA_LKT(29,3,1:6)=(/ 0.679707,0.540150,0.281157,0.107080,0.035090,0.013190 /)
XEXT_COEFF_550_LKT(29,3)=1145.100000 !rg=0.0979445 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,4,1:6)=(/ 4932.700000,2297.000000,719.710000,142.100000,16.962000,141.750000 /)
XPIZA_LKT(29,4,1:6)=(/ 0.985644,0.992433,0.999419,0.998227,0.966157,0.024652 /)
XCGA_LKT(29,4,1:6)=(/ 0.707283,0.604207,0.403283,0.178003,0.059533,0.022553 /)
XEXT_COEFF_550_LKT(29,4)=1482.200000 !rg=0.0979445 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,5,1:6)=(/ 5024.700000,2702.000000,985.640000,224.340000,29.752000,146.990000 /)
XPIZA_LKT(29,5,1:6)=(/ 0.985531,0.993232,0.999541,0.998807,0.980084,0.043477 /)
XCGA_LKT(29,5,1:6)=(/ 0.725670,0.654157,0.508620,0.282650,0.099277,0.037940 /)
XEXT_COEFF_550_LKT(29,5)=1862.000000 !rg=0.0979445 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,6,1:6)=(/ 4844.400000,3011.000000,1285.500000,340.510000,51.244000,155.710000 /)
XPIZA_LKT(29,6,1:6)=(/ 0.984634,0.993641,0.999621,0.999154,0.987947,0.075124 /)
XCGA_LKT(29,6,1:6)=(/ 0.734060,0.689397,0.585900,0.392153,0.155013,0.059573 /)
XEXT_COEFF_550_LKT(29,6)=2216.200000 !rg=0.0979445 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,7,1:6)=(/ 4433.500000,3157.800000,1581.800000,491.190000,84.473000,169.860000 /)
XPIZA_LKT(29,7,1:6)=(/ 0.982864,0.993684,0.999670,0.999364,0.992293,0.124028 /)
XCGA_LKT(29,7,1:6)=(/ 0.734940,0.711310,0.641580,0.487627,0.232257,0.089830 /)
XEXT_COEFF_550_LKT(29,7)=2484.300000 !rg=0.0979445 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,8,1:6)=(/ 3904.700000,3103.300000,1821.800000,667.480000,132.830000,191.740000 /)
XPIZA_LKT(29,8,1:6)=(/ 0.980246,0.993342,0.999697,0.999493,0.994759,0.191128 /)
XCGA_LKT(29,8,1:6)=(/ 0.732357,0.721167,0.679760,0.564377,0.332303,0.133187 /)
XEXT_COEFF_550_LKT(29,8)=2602.700000 !rg=0.0979445 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,9,1:6)=(/ 3348.700000,2880.200000,1960.200000,853.760000,201.940000,224.170000 /)
XPIZA_LKT(29,9,1:6)=(/ 0.976613,0.992626,0.999704,0.999573,0.996267,0.271842 /)
XCGA_LKT(29,9,1:6)=(/ 0.729097,0.723257,0.703503,0.623080,0.430847,0.197120 /)
XEXT_COEFF_550_LKT(29,9)=2550.000000 !rg=0.0979445 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,10,1:6)=(/ 2800.900000,2571.800000,1963.700000,1021.900000,287.490000,269.050000 /)
XPIZA_LKT(29,10,1:6)=(/ 0.972477,0.991571,0.999691,0.999618,0.997165,0.354679 /)
XCGA_LKT(29,10,1:6)=(/ 0.728317,0.722413,0.714250,0.666030,0.509737,0.287117 /)
XEXT_COEFF_550_LKT(29,10)=2361.700000 !rg=0.0979445 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,11,1:6)=(/ 2320.200000,2239.400000,1846.400000,1143.800000,389.250000,327.430000 /)
XPIZA_LKT(29,11,1:6)=(/ 0.966806,0.990117,0.999660,0.999639,0.997729,0.432517 /)
XCGA_LKT(29,11,1:6)=(/ 0.728373,0.722220,0.716003,0.694997,0.581050,0.385677 /)
XEXT_COEFF_550_LKT(29,11)=2107.500000 !rg=0.0979445 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,12,1:6)=(/ 1898.800000,1890.600000,1669.100000,1189.200000,494.020000,390.880000 /)
XPIZA_LKT(29,12,1:6)=(/ 0.961089,0.988028,0.999611,0.999636,0.998073,0.497924 /)
XCGA_LKT(29,12,1:6)=(/ 0.732397,0.720397,0.714667,0.710377,0.634317,0.464730 /)
XEXT_COEFF_550_LKT(29,12)=1836.800000 !rg=0.0979445 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,13,1:6)=(/ 1552.800000,1575.700000,1470.800000,1151.300000,588.050000,452.690000 /)
XPIZA_LKT(29,13,1:6)=(/ 0.953531,0.985946,0.999548,0.999607,0.998272,0.547220 /)
XCGA_LKT(29,13,1:6)=(/ 0.737267,0.724220,0.715657,0.713930,0.674507,0.535173 /)
XEXT_COEFF_550_LKT(29,13)=1550.400000 !rg=0.0979445 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,14,1:6)=(/ 1266.100000,1300.000000,1260.900000,1053.100000,652.280000,500.930000 /)
XPIZA_LKT(29,14,1:6)=(/ 0.944691,0.983292,0.999473,0.999559,0.998351,0.581727 /)
XCGA_LKT(29,14,1:6)=(/ 0.743193,0.725930,0.717287,0.711770,0.700500,0.590477 /)
XEXT_COEFF_550_LKT(29,14)=1297.400000 !rg=0.0979445 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,15,1:6)=(/ 1035.900000,1068.900000,1064.400000,942.770000,671.290000,527.900000 /)
XPIZA_LKT(29,15,1:6)=(/ 0.934576,0.979649,0.999359,0.999491,0.998321,0.602540 /)
XCGA_LKT(29,15,1:6)=(/ 0.753320,0.728900,0.717230,0.711133,0.713247,0.633070 /)
XEXT_COEFF_550_LKT(29,15)=1079.400000 !rg=0.0979445 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,16,1:6)=(/ 839.660000,877.350000,886.780000,828.160000,642.570000,528.310000 /)
XPIZA_LKT(29,16,1:6)=(/ 0.923640,0.975181,0.999221,0.999404,0.998176,0.611016 /)
XCGA_LKT(29,16,1:6)=(/ 0.763053,0.735743,0.717673,0.713357,0.714227,0.664007 /)
XEXT_COEFF_550_LKT(29,16)=882.860000 !rg=0.0979445 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,17,1:6)=(/ 684.560000,711.370000,728.680000,705.410000,584.120000,503.540000 /)
XPIZA_LKT(29,17,1:6)=(/ 0.912275,0.971901,0.999117,0.999298,0.997923,0.609270 /)
XCGA_LKT(29,17,1:6)=(/ 0.773353,0.744867,0.726563,0.714597,0.709500,0.686123 /)
XEXT_COEFF_550_LKT(29,17)=721.860000 !rg=0.0979445 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,18,1:6)=(/ 562.090000,583.540000,603.460000,597.690000,521.220000,462.370000 /)
XPIZA_LKT(29,18,1:6)=(/ 0.897239,0.965805,0.998836,0.999118,0.997598,0.601552 /)
XCGA_LKT(29,18,1:6)=(/ 0.785410,0.750487,0.729080,0.718847,0.707303,0.704267 /)
XEXT_COEFF_550_LKT(29,18)=595.350000 !rg=0.0979445 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,19,1:6)=(/ 457.740000,477.570000,493.550000,496.990000,459.310000,414.240000 /)
XPIZA_LKT(29,19,1:6)=(/ 0.882108,0.959652,0.998635,0.999021,0.997246,0.592468 /)
XCGA_LKT(29,19,1:6)=(/ 0.796317,0.760393,0.731827,0.722297,0.713243,0.722847 /)
XEXT_COEFF_550_LKT(29,19)=482.430000 !rg=0.0979445 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,20,1:6)=(/ 375.040000,387.870000,400.670000,413.080000,393.310000,364.600000 /)
XPIZA_LKT(29,20,1:6)=(/ 0.870015,0.953212,0.998496,0.998802,0.996771,0.583896 /)
XCGA_LKT(29,20,1:6)=(/ 0.806617,0.771900,0.745047,0.723883,0.716197,0.742110 /)
XEXT_COEFF_550_LKT(29,20)=394.000000 !rg=0.0979445 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,1,1:6)=(/ 4573.400000,1566.600000,406.520000,62.027000,7.033100,137.330000 /)
XPIZA_LKT(30,1,1:6)=(/ 0.985016,0.990120,0.999074,0.996216,0.921075,0.009909 /)
XCGA_LKT(30,1,1:6)=(/ 0.645513,0.454283,0.176433,0.065477,0.021323,0.007980 /)
XEXT_COEFF_550_LKT(30,1)=950.390000 !rg=0.106109 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,2,1:6)=(/ 4682.200000,1783.900000,480.350000,77.517000,8.724300,138.160000 /)
XPIZA_LKT(30,2,1:6)=(/ 0.985251,0.990853,0.999198,0.996919,0.935909,0.012431 /)
XCGA_LKT(30,2,1:6)=(/ 0.669777,0.519643,0.224730,0.083000,0.027060,0.010140 /)
XEXT_COEFF_550_LKT(30,2)=1072.500000 !rg=0.106109 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,3,1:6)=(/ 4959.500000,2147.700000,621.860000,112.640000,12.830000,140.050000 /)
XPIZA_LKT(30,3,1:6)=(/ 0.985779,0.992004,0.999349,0.997809,0.955756,0.018500 /)
XCGA_LKT(30,3,1:6)=(/ 0.699287,0.579307,0.326717,0.125423,0.041117,0.015470 /)
XEXT_COEFF_550_LKT(30,3)=1333.700000 !rg=0.106109 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,4,1:6)=(/ 5142.100000,2547.600000,842.040000,174.400000,21.325000,143.650000 /)
XPIZA_LKT(30,4,1:6)=(/ 0.985967,0.992942,0.999483,0.998516,0.972720,0.030968 /)
XCGA_LKT(30,4,1:6)=(/ 0.721600,0.631213,0.446303,0.206693,0.069653,0.026453 /)
XEXT_COEFF_550_LKT(30,4)=1684.900000 !rg=0.106109 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,5,1:6)=(/ 5086.200000,2907.900000,1123.900000,269.310000,37.147000,150.040000 /)
XPIZA_LKT(30,5,1:6)=(/ 0.985482,0.993535,0.999582,0.998975,0.983777,0.054173 /)
XCGA_LKT(30,5,1:6)=(/ 0.734477,0.673157,0.541640,0.318840,0.116003,0.044503 /)
XEXT_COEFF_550_LKT(30,5)=2055.100000 !rg=0.106109 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,6,1:6)=(/ 4768.700000,3143.600000,1424.900000,399.260000,63.014000,160.650000 /)
XPIZA_LKT(30,6,1:6)=(/ 0.984178,0.993770,0.999646,0.999253,0.989991,0.092303 /)
XCGA_LKT(30,6,1:6)=(/ 0.738300,0.701987,0.610177,0.427383,0.181017,0.069930 /)
XEXT_COEFF_550_LKT(30,6)=2374.200000 !rg=0.106109 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,7,1:6)=(/ 4267.200000,3197.400000,1703.700000,563.170000,102.090000,177.710000 /)
XPIZA_LKT(30,7,1:6)=(/ 0.981986,0.993649,0.999685,0.999425,0.993453,0.149073 /)
XCGA_LKT(30,7,1:6)=(/ 0.736253,0.718493,0.658890,0.518860,0.269493,0.105630 /)
XEXT_COEFF_550_LKT(30,7)=2581.300000 !rg=0.106109 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,8,1:6)=(/ 3695.800000,3056.100000,1905.500000,745.270000,158.680000,203.600000 /)
XPIZA_LKT(30,8,1:6)=(/ 0.978999,0.993153,0.999703,0.999530,0.995468,0.222832 /)
XCGA_LKT(30,8,1:6)=(/ 0.732820,0.724433,0.691150,0.588483,0.373097,0.157010 /)
XEXT_COEFF_550_LKT(30,8)=2624.300000 !rg=0.106109 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,9,1:6)=(/ 3118.600000,2781.600000,1989.300000,928.250000,235.340000,241.010000 /)
XPIZA_LKT(30,9,1:6)=(/ 0.974992,0.992253,0.999702,0.999595,0.996693,0.306003 /)
XCGA_LKT(30,9,1:6)=(/ 0.729250,0.723900,0.709613,0.641593,0.462800,0.232000 /)
XEXT_COEFF_550_LKT(30,9)=2502.900000 !rg=0.106109 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,10,1:6)=(/ 2595.500000,2451.200000,1938.100000,1080.800000,327.650000,291.580000 /)
XPIZA_LKT(30,10,1:6)=(/ 0.970060,0.991071,0.999682,0.999630,0.997424,0.387395 /)
XCGA_LKT(30,10,1:6)=(/ 0.727317,0.723720,0.716287,0.678840,0.539873,0.329067 /)
XEXT_COEFF_550_LKT(30,10)=2276.400000 !rg=0.106109 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,11,1:6)=(/ 2133.000000,2098.500000,1788.100000,1174.500000,431.860000,353.190000 /)
XPIZA_LKT(30,11,1:6)=(/ 0.964637,0.989417,0.999643,0.999641,0.997890,0.461101 /)
XCGA_LKT(30,11,1:6)=(/ 0.729887,0.722510,0.716283,0.702353,0.603480,0.419487 /)
XEXT_COEFF_550_LKT(30,11)=2009.300000 !rg=0.106109 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,12,1:6)=(/ 1745.200000,1755.900000,1594.600000,1185.800000,534.910000,416.530000 /)
XPIZA_LKT(30,12,1:6)=(/ 0.957709,0.987473,0.999591,0.999628,0.998171,0.519621 /)
XCGA_LKT(30,12,1:6)=(/ 0.733200,0.723170,0.715793,0.712997,0.652057,0.493727 /)
XEXT_COEFF_550_LKT(30,12)=1723.000000 !rg=0.106109 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,13,1:6)=(/ 1419.300000,1456.800000,1390.700000,1118.200000,619.230000,474.610000 /)
XPIZA_LKT(30,13,1:6)=(/ 0.949764,0.985049,0.999522,0.999591,0.998319,0.563140 /)
XCGA_LKT(30,13,1:6)=(/ 0.738673,0.724557,0.717330,0.713857,0.686357,0.559223 /)
XEXT_COEFF_550_LKT(30,13)=1443.000000 !rg=0.106109 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,14,1:6)=(/ 1162.800000,1195.800000,1179.700000,1012.200000,666.100000,515.320000 /)
XPIZA_LKT(30,14,1:6)=(/ 0.940109,0.982057,0.999425,0.999533,0.998353,0.591738 /)
XCGA_LKT(30,14,1:6)=(/ 0.747920,0.726573,0.716140,0.711163,0.706877,0.609047 /)
XEXT_COEFF_550_LKT(30,14)=1207.600000 !rg=0.106109 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,15,1:6)=(/ 943.560000,985.420000,989.200000,898.480000,665.120000,531.510000 /)
XPIZA_LKT(30,15,1:6)=(/ 0.930279,0.977834,0.999279,0.999461,0.998275,0.607237 /)
XCGA_LKT(30,15,1:6)=(/ 0.756023,0.732170,0.716860,0.712303,0.714430,0.646510 /)
XEXT_COEFF_550_LKT(30,15)=992.380000 !rg=0.106109 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,16,1:6)=(/ 765.630000,800.260000,815.940000,779.170000,621.740000,521.050000 /)
XPIZA_LKT(30,16,1:6)=(/ 0.919527,0.973587,0.999185,0.999366,0.998080,0.611337 /)
XCGA_LKT(30,16,1:6)=(/ 0.768093,0.738307,0.723573,0.714707,0.712270,0.673643 /)
XEXT_COEFF_550_LKT(30,16)=807.020000 !rg=0.106109 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,17,1:6)=(/ 629.160000,651.690000,674.960000,662.170000,558.290000,488.300000 /)
XPIZA_LKT(30,17,1:6)=(/ 0.905475,0.969907,0.998932,0.999184,0.997798,0.606717 /)
XCGA_LKT(30,17,1:6)=(/ 0.779300,0.747150,0.726170,0.714757,0.707737,0.693717 /)
XEXT_COEFF_550_LKT(30,17)=665.720000 !rg=0.106109 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,18,1:6)=(/ 512.530000,535.650000,554.700000,553.890000,496.700000,443.270000 /)
XPIZA_LKT(30,18,1:6)=(/ 0.890833,0.963696,0.998750,0.999110,0.997464,0.598064 /)
XCGA_LKT(30,18,1:6)=(/ 0.789000,0.756073,0.728307,0.721013,0.710300,0.711797 /)
XEXT_COEFF_550_LKT(30,18)=543.050000 !rg=0.106109 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,19,1:6)=(/ 418.630000,435.530000,451.500000,459.870000,433.310000,394.010000 /)
XPIZA_LKT(30,19,1:6)=(/ 0.877057,0.956630,0.998583,0.998765,0.997045,0.589203 /)
XCGA_LKT(30,19,1:6)=(/ 0.801497,0.763923,0.739927,0.721657,0.714687,0.731000 /)
XEXT_COEFF_550_LKT(30,19)=439.390000 !rg=0.106109 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,20,1:6)=(/ 344.720000,355.480000,370.360000,382.930000,367.160000,344.130000 /)
XPIZA_LKT(30,20,1:6)=(/ 0.863642,0.950081,0.998304,0.998380,0.996546,0.580553 /)
XCGA_LKT(30,20,1:6)=(/ 0.812207,0.774667,0.745410,0.725830,0.716133,0.749990 /)
XEXT_COEFF_550_LKT(30,20)=362.070000 !rg=0.106109 sigma=2.95 wvl=0.55

END SUBROUTINE SALT_OPT_LKT_SET3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET4()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE

XEXT_COEFF_WVL_LKT(31,1,1:6)=(/ 4821.300000,1811.400000,498.090000,78.525000,8.797300,138.270000 /)
XPIZA_LKT(31,1,1:6)=(/ 0.985983,0.990884,0.999220,0.996946,0.936350,0.012524 /)
XCGA_LKT(31,1,1:6)=(/ 0.659013,0.532013,0.209203,0.076780,0.025003,0.009363 /)
XEXT_COEFF_550_LKT(31,1)=1109.600000 !rg=0.114954 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,2,1:6)=(/ 5022.900000,2071.300000,578.440000,97.777000,10.945000,139.280000 /)
XPIZA_LKT(31,2,1:6)=(/ 0.986000,0.991675,0.999311,0.997500,0.948408,0.015695 /)
XCGA_LKT(31,2,1:6)=(/ 0.692153,0.571677,0.266337,0.097357,0.031720,0.011897 /)
XEXT_COEFF_550_LKT(31,2)=1252.000000 !rg=0.114954 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,3,1:6)=(/ 5233.500000,2424.800000,735.920000,140.300000,16.144000,141.580000 /)
XPIZA_LKT(31,3,1:6)=(/ 0.986266,0.992638,0.999428,0.998196,0.964416,0.023299 /)
XCGA_LKT(31,3,1:6)=(/ 0.717070,0.611983,0.376357,0.146943,0.048167,0.018147 /)
XEXT_COEFF_550_LKT(31,3)=1541.700000 !rg=0.114954 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,4,1:6)=(/ 5286.300000,2791.700000,976.910000,212.480000,26.793000,145.960000 /)
XPIZA_LKT(31,4,1:6)=(/ 0.986119,0.993360,0.999536,0.998746,0.977958,0.038796 /)
XCGA_LKT(31,4,1:6)=(/ 0.733700,0.654890,0.487920,0.238980,0.081457,0.031027 /)
XEXT_COEFF_550_LKT(31,4)=1895.100000 !rg=0.114954 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,5,1:6)=(/ 5079.400000,3092.200000,1269.100000,320.650000,46.191000,153.760000 /)
XPIZA_LKT(31,5,1:6)=(/ 0.985251,0.993766,0.999616,0.999110,0.986707,0.067181 /)
XCGA_LKT(31,5,1:6)=(/ 0.741230,0.689790,0.572130,0.356047,0.135450,0.052207 /)
XEXT_COEFF_550_LKT(31,5)=2242.900000 !rg=0.114954 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,6,1:6)=(/ 4636.200000,3242.400000,1563.700000,464.650000,77.003000,166.640000 /)
XPIZA_LKT(31,6,1:6)=(/ 0.983543,0.993835,0.999667,0.999335,0.991618,0.112526 /)
XCGA_LKT(31,6,1:6)=(/ 0.740883,0.712593,0.632150,0.462337,0.210903,0.082117 /)
XEXT_COEFF_550_LKT(31,6)=2513.900000 !rg=0.114954 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,7,1:6)=(/ 4061.600000,3201.500000,1815.500000,638.860000,122.880000,187.070000 /)
XPIZA_LKT(31,7,1:6)=(/ 0.980958,0.993544,0.999696,0.999475,0.994399,0.177101 /)
XCGA_LKT(31,7,1:6)=(/ 0.736670,0.724020,0.674180,0.547200,0.309157,0.124290 /)
XEXT_COEFF_550_LKT(31,7)=2650.800000 !rg=0.114954 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,8,1:6)=(/ 3466.500000,2984.800000,1970.600000,824.720000,187.990000,217.380000 /)
XPIZA_LKT(31,8,1:6)=(/ 0.977444,0.992882,0.999706,0.999561,0.996047,0.256074 /)
XCGA_LKT(31,8,1:6)=(/ 0.731770,0.726683,0.700743,0.610930,0.409827,0.185053 /)
XEXT_COEFF_550_LKT(31,8)=2616.800000 !rg=0.114954 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,9,1:6)=(/ 2897.600000,2666.000000,1996.900000,999.530000,271.200000,260.170000 /)
XPIZA_LKT(31,9,1:6)=(/ 0.972642,0.991856,0.999697,0.999613,0.997035,0.339881 /)
XCGA_LKT(31,9,1:6)=(/ 0.728427,0.725450,0.714230,0.657930,0.493477,0.270877 /)
XEXT_COEFF_550_LKT(31,9)=2439.300000 !rg=0.114954 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,10,1:6)=(/ 2390.700000,2319.200000,1896.600000,1131.000000,370.140000,315.990000 /)
XPIZA_LKT(31,10,1:6)=(/ 0.967737,0.990434,0.999670,0.999638,0.997643,0.418875 /)
XCGA_LKT(31,10,1:6)=(/ 0.727733,0.723430,0.717820,0.689750,0.567593,0.368267 /)
XEXT_COEFF_550_LKT(31,10)=2182.600000 !rg=0.114954 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,11,1:6)=(/ 1962.400000,1956.800000,1720.400000,1192.100000,474.980000,378.940000 /)
XPIZA_LKT(31,11,1:6)=(/ 0.961327,0.988744,0.999625,0.999639,0.998020,0.486872 /)
XCGA_LKT(31,11,1:6)=(/ 0.729947,0.723067,0.716367,0.708117,0.624003,0.449793 /)
XEXT_COEFF_550_LKT(31,11)=1903.400000 !rg=0.114954 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,12,1:6)=(/ 1597.000000,1629.900000,1520.000000,1168.800000,572.360000,441.700000 /)
XPIZA_LKT(31,12,1:6)=(/ 0.954121,0.986489,0.999566,0.999616,0.998243,0.538925 /)
XCGA_LKT(31,12,1:6)=(/ 0.734483,0.722873,0.717460,0.714327,0.667100,0.521853 /)
XEXT_COEFF_550_LKT(31,12)=1608.000000 !rg=0.114954 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,13,1:6)=(/ 1305.400000,1341.200000,1307.600000,1081.600000,643.950000,493.590000 /)
XPIZA_LKT(31,13,1:6)=(/ 0.945195,0.983844,0.999481,0.999569,0.998345,0.576414 /)
XCGA_LKT(31,13,1:6)=(/ 0.742820,0.724687,0.716553,0.712750,0.695933,0.580230 /)
XEXT_COEFF_550_LKT(31,13)=1346.900000 !rg=0.114954 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,14,1:6)=(/ 1059.700000,1106.500000,1096.900000,967.800000,672.360000,525.290000 /)
XPIZA_LKT(31,14,1:6)=(/ 0.936001,0.980089,0.999373,0.999507,0.998334,0.599691 /)
XCGA_LKT(31,14,1:6)=(/ 0.749607,0.728840,0.716360,0.711670,0.711230,0.625310 /)
XEXT_COEFF_550_LKT(31,14)=1111.300000 !rg=0.114954 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,15,1:6)=(/ 859.170000,900.540000,919.470000,853.610000,652.010000,530.560000 /)
XPIZA_LKT(31,15,1:6)=(/ 0.926492,0.976629,0.999186,0.999426,0.998209,0.610265 /)
XCGA_LKT(31,15,1:6)=(/ 0.763103,0.733113,0.719707,0.714380,0.714387,0.658237 /)
XEXT_COEFF_550_LKT(31,15)=910.240000 !rg=0.114954 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,16,1:6)=(/ 703.830000,730.210000,752.630000,731.120000,596.760000,510.070000 /)
XPIZA_LKT(31,16,1:6)=(/ 0.913430,0.972479,0.999127,0.999305,0.997977,0.610417 /)
XCGA_LKT(31,16,1:6)=(/ 0.773263,0.743600,0.724527,0.714637,0.710177,0.682200 /)
XEXT_COEFF_550_LKT(31,16)=745.560000 !rg=0.114954 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,17,1:6)=(/ 574.350000,600.650000,620.140000,614.640000,533.080000,471.040000 /)
XPIZA_LKT(31,17,1:6)=(/ 0.898959,0.966718,0.998767,0.999185,0.997678,0.603619 /)
XCGA_LKT(31,17,1:6)=(/ 0.781933,0.751057,0.725183,0.717263,0.708993,0.701133 /)
XEXT_COEFF_550_LKT(31,17)=610.060000 !rg=0.114954 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,18,1:6)=(/ 467.670000,488.330000,510.130000,512.530000,472.670000,423.550000 /)
XPIZA_LKT(31,18,1:6)=(/ 0.884803,0.961429,0.998607,0.999046,0.997320,0.594672 /)
XCGA_LKT(31,18,1:6)=(/ 0.796700,0.757683,0.732347,0.721143,0.712467,0.719693 /)
XEXT_COEFF_550_LKT(31,18)=495.930000 !rg=0.114954 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,19,1:6)=(/ 384.770000,397.180000,413.680000,427.870000,404.920000,373.570000 /)
XPIZA_LKT(31,19,1:6)=(/ 0.871402,0.954344,0.998348,0.998703,0.996852,0.585973 /)
XCGA_LKT(31,19,1:6)=(/ 0.806830,0.770337,0.742810,0.724257,0.715213,0.739180 /)
XEXT_COEFF_550_LKT(31,19)=405.810000 !rg=0.114954 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,20,1:6)=(/ 315.530000,327.180000,340.200000,350.300000,342.870000,323.780000 /)
XPIZA_LKT(31,20,1:6)=(/ 0.858497,0.945883,0.997955,0.998576,0.996303,0.577079 /)
XCGA_LKT(31,20,1:6)=(/ 0.815180,0.779217,0.745133,0.729503,0.717550,0.757650 /)
XEXT_COEFF_550_LKT(31,20)=332.080000 !rg=0.114954 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,1,1:6)=(/ 5091.300000,2125.800000,602.130000,99.237000,11.039000,139.410000 /)
XPIZA_LKT(32,1,1:6)=(/ 0.986283,0.991634,0.999333,0.997524,0.948762,0.015810 /)
XCGA_LKT(32,1,1:6)=(/ 0.693987,0.597433,0.249170,0.090083,0.029310,0.010980 /)
XEXT_COEFF_550_LKT(32,1)=1281.400000 !rg=0.124536 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,2,1:6)=(/ 5331.000000,2385.100000,687.930000,122.920000,13.764000,140.640000 /)
XPIZA_LKT(32,2,1:6)=(/ 0.986526,0.992422,0.999399,0.997960,0.958509,0.019785 /)
XCGA_LKT(32,2,1:6)=(/ 0.715527,0.609343,0.315870,0.114277,0.037177,0.013953 /)
XEXT_COEFF_550_LKT(32,2)=1458.600000 !rg=0.124536 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,3,1:6)=(/ 5446.300000,2701.300000,864.350000,173.600000,20.336000,143.430000 /)
XPIZA_LKT(32,3,1:6)=(/ 0.986573,0.993170,0.999491,0.998501,0.971362,0.029279 /)
XCGA_LKT(32,3,1:6)=(/ 0.732717,0.639223,0.428023,0.172203,0.056417,0.021283 /)
XEXT_COEFF_550_LKT(32,3)=1765.500000 !rg=0.124536 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,4,1:6)=(/ 5356.300000,3021.100000,1123.300000,256.880000,33.609000,148.770000 /)
XPIZA_LKT(32,4,1:6)=(/ 0.986095,0.993695,0.999579,0.998930,0.982129,0.048438 /)
XCGA_LKT(32,4,1:6)=(/ 0.743513,0.675660,0.526973,0.274760,0.095217,0.036390 /)
XEXT_COEFF_550_LKT(32,4)=2107.300000 !rg=0.124536 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,5,1:6)=(/ 5001.000000,3247.600000,1418.400000,378.800000,57.151000,158.270000 /)
XPIZA_LKT(32,5,1:6)=(/ 0.984841,0.993929,0.999644,0.999219,0.989032,0.082831 /)
XCGA_LKT(32,5,1:6)=(/ 0.746013,0.704193,0.599863,0.394053,0.157980,0.061250 /)
XEXT_COEFF_550_LKT(32,5)=2419.600000 !rg=0.124536 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,6,1:6)=(/ 4453.700000,3303.300000,1696.900000,536.150000,93.569000,173.840000 /)
XPIZA_LKT(32,6,1:6)=(/ 0.982724,0.993834,0.999684,0.999402,0.992924,0.135909 /)
XCGA_LKT(32,6,1:6)=(/ 0.742437,0.721300,0.651767,0.495767,0.244527,0.096457 /)
XEXT_COEFF_550_LKT(32,6)=2630.500000 !rg=0.124536 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,7,1:6)=(/ 3832.300000,3168.100000,1912.200000,717.910000,147.220000,198.110000 /)
XPIZA_LKT(32,7,1:6)=(/ 0.979692,0.993384,0.999704,0.999517,0.995176,0.207688 /)
XCGA_LKT(32,7,1:6)=(/ 0.736597,0.728213,0.687330,0.573343,0.348483,0.146297 /)
XEXT_COEFF_550_LKT(32,7)=2689.800000 !rg=0.124536 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,8,1:6)=(/ 3222.800000,2887.900000,2013.700000,902.290000,220.100000,233.260000 /)
XPIZA_LKT(32,8,1:6)=(/ 0.975759,0.992550,0.999706,0.999587,0.996513,0.290068 /)
XCGA_LKT(32,8,1:6)=(/ 0.731177,0.727890,0.708387,0.630980,0.442930,0.217553 /)
XEXT_COEFF_550_LKT(32,8)=2582.600000 !rg=0.124536 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,9,1:6)=(/ 2666.400000,2542.900000,1982.100000,1063.700000,310.360000,281.710000 /)
XPIZA_LKT(32,9,1:6)=(/ 0.970724,0.991355,0.999690,0.999627,0.997317,0.373112 /)
XCGA_LKT(32,9,1:6)=(/ 0.727693,0.725867,0.717587,0.672090,0.524243,0.311367 /)
XEXT_COEFF_550_LKT(32,9)=2355.100000 !rg=0.124536 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,10,1:6)=(/ 2204.000000,2176.100000,1843.100000,1168.900000,412.860000,341.370000 /)
XPIZA_LKT(32,10,1:6)=(/ 0.964557,0.989790,0.999654,0.999642,0.997820,0.448407 /)
XCGA_LKT(32,10,1:6)=(/ 0.728470,0.724010,0.718023,0.698427,0.591327,0.402823 /)
XEXT_COEFF_550_LKT(32,10)=2084.300000 !rg=0.124536 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,11,1:6)=(/ 1793.200000,1822.100000,1648.700000,1196.000000,517.220000,404.780000 /)
XPIZA_LKT(32,11,1:6)=(/ 0.957963,0.987648,0.999604,0.999632,0.998129,0.509753 /)
XCGA_LKT(32,11,1:6)=(/ 0.730300,0.721647,0.717800,0.711777,0.642830,0.479040 /)
XEXT_COEFF_550_LKT(32,11)=1786.900000 !rg=0.124536 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,12,1:6)=(/ 1465.800000,1506.200000,1441.800000,1143.000000,606.320000,464.900000 /)
XPIZA_LKT(32,12,1:6)=(/ 0.949864,0.985425,0.999531,0.999600,0.998299,0.555957 /)
XCGA_LKT(32,12,1:6)=(/ 0.737880,0.723397,0.716907,0.714437,0.680033,0.547037 /)
XEXT_COEFF_550_LKT(32,12)=1501.100000 !rg=0.124536 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,13,1:6)=(/ 1188.800000,1239.500000,1218.300000,1038.800000,661.780000,509.750000 /)
XPIZA_LKT(32,13,1:6)=(/ 0.941338,0.982080,0.999447,0.999547,0.998355,0.587392 /)
XCGA_LKT(32,13,1:6)=(/ 0.744013,0.725800,0.717180,0.712377,0.703517,0.599573 /)
XEXT_COEFF_550_LKT(32,13)=1244.600000 !rg=0.124536 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,14,1:6)=(/ 964.600000,1012.200000,1025.700000,924.890000,669.870000,530.970000 /)
XPIZA_LKT(32,14,1:6)=(/ 0.932234,0.978710,0.999288,0.999477,0.998298,0.605349 /)
XCGA_LKT(32,14,1:6)=(/ 0.757343,0.728563,0.717690,0.713430,0.713603,0.639537 /)
XEXT_COEFF_550_LKT(32,14)=1023.900000 !rg=0.124536 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,15,1:6)=(/ 788.940000,821.340000,844.730000,806.680000,632.690000,525.290000 /)
XPIZA_LKT(32,15,1:6)=(/ 0.920622,0.975509,0.999117,0.999380,0.998130,0.611435 /)
XCGA_LKT(32,15,1:6)=(/ 0.766977,0.740253,0.720820,0.714427,0.713197,0.668513 /)
XEXT_COEFF_550_LKT(32,15)=835.380000 !rg=0.124536 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,16,1:6)=(/ 643.330000,673.150000,693.190000,680.900000,571.560000,496.110000 /)
XPIZA_LKT(32,16,1:6)=(/ 0.906920,0.970010,0.999064,0.999245,0.997864,0.608375 /)
XCGA_LKT(32,16,1:6)=(/ 0.775620,0.745990,0.725480,0.713887,0.709660,0.690163 /)
XEXT_COEFF_550_LKT(32,16)=683.850000 !rg=0.124536 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,17,1:6)=(/ 524.720000,548.030000,571.710000,571.000000,510.170000,452.450000 /)
XPIZA_LKT(32,17,1:6)=(/ 0.892157,0.964651,0.998830,0.999107,0.997527,0.600299 /)
XCGA_LKT(32,17,1:6)=(/ 0.790277,0.751723,0.730493,0.720983,0.709650,0.708740 /)
XEXT_COEFF_550_LKT(32,17)=559.920000 !rg=0.124536 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,18,1:6)=(/ 430.300000,444.880000,466.030000,476.440000,445.410000,403.430000 /)
XPIZA_LKT(32,18,1:6)=(/ 0.878543,0.958924,0.998453,0.998915,0.997136,0.591412 /)
XCGA_LKT(32,18,1:6)=(/ 0.801330,0.766143,0.735433,0.722023,0.714387,0.727870 /)
XEXT_COEFF_550_LKT(32,18)=453.590000 !rg=0.124536 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,19,1:6)=(/ 352.680000,365.770000,380.020000,392.640000,377.840000,353.000000 /)
XPIZA_LKT(32,19,1:6)=(/ 0.864947,0.950629,0.998315,0.998711,0.996673,0.582642 /)
XCGA_LKT(32,19,1:6)=(/ 0.809630,0.774207,0.743663,0.724873,0.716377,0.747197 /)
XEXT_COEFF_550_LKT(32,19)=372.160000 !rg=0.124536 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,20,1:6)=(/ 288.970000,299.430000,311.870000,321.530000,322.010000,303.770000 /)
XPIZA_LKT(32,20,1:6)=(/ 0.853031,0.942933,0.998043,0.998182,0.995768,0.573559 /)
XCGA_LKT(32,20,1:6)=(/ 0.821697,0.781017,0.751843,0.734313,0.719153,0.765183 /)
XEXT_COEFF_550_LKT(32,20)=304.810000 !rg=0.124536 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,1,1:6)=(/ 5524.700000,2514.800000,716.420000,125.090000,13.886000,140.790000 /)
XPIZA_LKT(33,1,1:6)=(/ 0.986602,0.992472,0.999420,0.997984,0.958794,0.019926 /)
XCGA_LKT(33,1,1:6)=(/ 0.733147,0.634190,0.298200,0.105777,0.034360,0.012880 /)
XEXT_COEFF_550_LKT(33,1)=1481.000000 !rg=0.134917 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,2,1:6)=(/ 5591.600000,2699.800000,809.370000,153.820000,17.341000,142.280000 /)
XPIZA_LKT(33,2,1:6)=(/ 0.986892,0.993079,0.999468,0.998325,0.966639,0.024897 /)
XCGA_LKT(33,2,1:6)=(/ 0.735837,0.634527,0.373393,0.134290,0.043570,0.016367 /)
XEXT_COEFF_550_LKT(33,2)=1696.700000 !rg=0.134917 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,3,1:6)=(/ 5585.600000,2966.900000,1008.100000,213.120000,25.623000,145.680000 /)
XPIZA_LKT(33,3,1:6)=(/ 0.986702,0.993606,0.999543,0.998743,0.976918,0.036702 /)
XCGA_LKT(33,3,1:6)=(/ 0.745790,0.662617,0.478880,0.201813,0.066070,0.024963 /)
XEXT_COEFF_550_LKT(33,3)=1999.000000 !rg=0.134917 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,4,1:6)=(/ 5345.600000,3227.400000,1279.200000,308.100000,42.050000,152.200000 /)
XPIZA_LKT(33,4,1:6)=(/ 0.985879,0.993954,0.999615,0.999077,0.985446,0.060228 /)
XCGA_LKT(33,4,1:6)=(/ 0.750927,0.693887,0.562617,0.313730,0.111240,0.042680 /)
XEXT_COEFF_550_LKT(33,4)=2315.200000 !rg=0.134917 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,5,1:6)=(/ 4855.300000,3367.600000,1568.500000,444.040000,70.325000,163.760000 /)
XPIZA_LKT(33,5,1:6)=(/ 0.984219,0.994025,0.999666,0.999309,0.990882,0.101417 /)
XCGA_LKT(33,5,1:6)=(/ 0.748907,0.716460,0.624847,0.432183,0.183917,0.071863 /)
XEXT_COEFF_550_LKT(33,5)=2579.100000 !rg=0.134917 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,6,1:6)=(/ 4232.500000,3323.800000,1821.600000,612.590000,113.160000,182.470000 /)
XPIZA_LKT(33,6,1:6)=(/ 0.981614,0.993771,0.999697,0.999457,0.993982,0.162394 /)
XCGA_LKT(33,6,1:6)=(/ 0.741980,0.728313,0.669200,0.526853,0.280873,0.113330 /)
XEXT_COEFF_550_LKT(33,6)=2719.300000 !rg=0.134917 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,7,1:6)=(/ 3584.300000,3103.100000,1991.100000,799.400000,175.010000,211.030000 /)
XPIZA_LKT(33,7,1:6)=(/ 0.978058,0.993152,0.999709,0.999551,0.995809,0.240201 /)
XCGA_LKT(33,7,1:6)=(/ 0.734537,0.731320,0.698613,0.597683,0.385260,0.172123 /)
XEXT_COEFF_550_LKT(33,7)=2698.100000 !rg=0.134917 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,8,1:6)=(/ 2984.400000,2771.400000,2034.700000,977.540000,254.940000,251.410000 /)
XPIZA_LKT(33,8,1:6)=(/ 0.973482,0.992165,0.999704,0.999607,0.996888,0.324163 /)
XCGA_LKT(33,8,1:6)=(/ 0.729223,0.729307,0.714427,0.648853,0.474690,0.253967 /)
XEXT_COEFF_550_LKT(33,8)=2526.100000 !rg=0.134917 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,9,1:6)=(/ 2449.300000,2400.000000,1950.500000,1119.900000,352.270000,305.320000 /)
XPIZA_LKT(33,9,1:6)=(/ 0.968683,0.990801,0.999679,0.999636,0.997554,0.405288 /)
XCGA_LKT(33,9,1:6)=(/ 0.728530,0.726413,0.719660,0.684283,0.553060,0.350117 /)
XEXT_COEFF_550_LKT(33,9)=2262.700000 !rg=0.134917 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,10,1:6)=(/ 2009.600000,2033.800000,1776.900000,1194.500000,456.190000,367.100000 /)
XPIZA_LKT(33,10,1:6)=(/ 0.961935,0.988869,0.999638,0.999641,0.997963,0.475317 /)
XCGA_LKT(33,10,1:6)=(/ 0.728247,0.723030,0.719330,0.705417,0.612813,0.433793 /)
XEXT_COEFF_550_LKT(33,10)=1969.200000 !rg=0.134917 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,11,1:6)=(/ 1645.100000,1685.200000,1574.900000,1186.700000,556.300000,430.500000 /)
XPIZA_LKT(33,11,1:6)=(/ 0.954049,0.986680,0.999579,0.999623,0.998212,0.530159 /)
XCGA_LKT(33,11,1:6)=(/ 0.733973,0.721957,0.718043,0.714230,0.658983,0.507667 /)
XEXT_COEFF_550_LKT(33,11)=1672.700000 !rg=0.134917 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,12,1:6)=(/ 1334.200000,1393.600000,1352.700000,1107.600000,634.020000,485.340000 /)
XPIZA_LKT(33,12,1:6)=(/ 0.946365,0.983667,0.999503,0.999582,0.998334,0.570405 /)
XCGA_LKT(33,12,1:6)=(/ 0.739220,0.723737,0.718367,0.714123,0.690623,0.568987 /)
XEXT_COEFF_550_LKT(33,12)=1390.800000 !rg=0.134917 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,13,1:6)=(/ 1082.200000,1136.700000,1139.800000,996.240000,671.720000,521.890000 /)
XPIZA_LKT(33,13,1:6)=(/ 0.937609,0.980682,0.999367,0.999523,0.998345,0.596385 /)
XCGA_LKT(33,13,1:6)=(/ 0.751490,0.725483,0.716813,0.713487,0.708737,0.616840 /)
XEXT_COEFF_550_LKT(33,13)=1149.400000 !rg=0.134917 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,14,1:6)=(/ 882.690000,924.420000,946.000000,880.910000,660.120000,532.300000 /)
XPIZA_LKT(33,14,1:6)=(/ 0.927188,0.977336,0.999267,0.999441,0.998243,0.609267 /)
XCGA_LKT(33,14,1:6)=(/ 0.761187,0.735437,0.719530,0.713837,0.714367,0.652127 /)
XEXT_COEFF_550_LKT(33,14)=938.490000 !rg=0.134917 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,15,1:6)=(/ 721.750000,753.920000,777.450000,754.000000,610.500000,516.190000 /)
XPIZA_LKT(33,15,1:6)=(/ 0.914271,0.973276,0.999132,0.999337,0.998031,0.611256 /)
XCGA_LKT(33,15,1:6)=(/ 0.769650,0.741823,0.722717,0.715170,0.711903,0.677707 /)
XEXT_COEFF_550_LKT(33,15)=765.820000 !rg=0.134917 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,16,1:6)=(/ 589.540000,615.110000,641.040000,633.210000,547.510000,479.810000 /)
XPIZA_LKT(33,16,1:6)=(/ 0.899930,0.967431,0.998908,0.999231,0.997727,0.605677 /)
XCGA_LKT(33,16,1:6)=(/ 0.783137,0.746950,0.726727,0.719810,0.708520,0.697830 /)
XEXT_COEFF_550_LKT(33,16)=627.620000 !rg=0.134917 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,17,1:6)=(/ 480.470000,501.670000,522.400000,531.370000,484.300000,433.010000 /)
XPIZA_LKT(33,17,1:6)=(/ 0.885892,0.961526,0.998754,0.999024,0.997387,0.596970 /)
XCGA_LKT(33,17,1:6)=(/ 0.795607,0.759177,0.731127,0.719977,0.712603,0.716660 /)
XEXT_COEFF_550_LKT(33,17)=511.530000 !rg=0.134917 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,18,1:6)=(/ 394.620000,409.020000,425.970000,439.710000,417.470000,383.000000 /)
XPIZA_LKT(33,18,1:6)=(/ 0.871878,0.955441,0.998496,0.998863,0.996974,0.588193 /)
XCGA_LKT(33,18,1:6)=(/ 0.804270,0.769660,0.739810,0.721637,0.715937,0.736140 /)
XEXT_COEFF_550_LKT(33,18)=416.250000 !rg=0.134917 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,19,1:6)=(/ 323.610000,334.650000,350.500000,359.600000,354.610000,332.450000 /)
XPIZA_LKT(33,19,1:6)=(/ 0.859273,0.947092,0.998098,0.998702,0.996347,0.579152 /)
XCGA_LKT(33,19,1:6)=(/ 0.816007,0.776067,0.746077,0.732503,0.717200,0.754990 /)
XEXT_COEFF_550_LKT(33,19)=341.390000 !rg=0.134917 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,20,1:6)=(/ 265.250000,273.910000,284.410000,297.580000,298.130000,284.320000 /)
XPIZA_LKT(33,20,1:6)=(/ 0.848053,0.939900,0.997802,0.998304,0.995625,0.570086 /)
XCGA_LKT(33,20,1:6)=(/ 0.826337,0.787700,0.753277,0.734487,0.718733,0.772690 /)
XEXT_COEFF_550_LKT(33,20)=278.460000 !rg=0.134917 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,1,1:6)=(/ 5861.500000,2898.600000,838.630000,157.070000,17.501000,142.460000 /)
XPIZA_LKT(34,1,1:6)=(/ 0.987237,0.993317,0.999486,0.998348,0.966869,0.025068 /)
XCGA_LKT(34,1,1:6)=(/ 0.752633,0.643480,0.358097,0.124363,0.040273,0.015107 /)
XEXT_COEFF_550_LKT(34,1)=1735.700000 !rg=0.146163 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,2,1:6)=(/ 5777.200000,2989.300000,945.620000,191.260000,21.873000,144.280000 /)
XPIZA_LKT(34,2,1:6)=(/ 0.987109,0.993619,0.999522,0.998613,0.973159,0.031259 /)
XCGA_LKT(34,2,1:6)=(/ 0.751070,0.653800,0.436687,0.158050,0.051060,0.019193 /)
XEXT_COEFF_550_LKT(34,2)=1959.700000 !rg=0.146163 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,3,1:6)=(/ 5639.900000,3212.100000,1166.900000,259.340000,32.271000,148.430000 /)
XPIZA_LKT(34,3,1:6)=(/ 0.986645,0.993952,0.999587,0.998934,0.981352,0.045864 /)
XCGA_LKT(34,3,1:6)=(/ 0.755987,0.683517,0.525980,0.236343,0.077370,0.029277 /)
XEXT_COEFF_550_LKT(34,3)=2233.700000 !rg=0.146163 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,4,1:6)=(/ 5250.100000,3402.400000,1441.600000,366.630000,52.419000,156.390000 /)
XPIZA_LKT(34,4,1:6)=(/ 0.985458,0.994141,0.999645,0.999196,0.988081,0.074513 /)
XCGA_LKT(34,4,1:6)=(/ 0.755903,0.709800,0.594473,0.355247,0.129873,0.050053 /)
XEXT_COEFF_550_LKT(34,4)=2511.900000 !rg=0.146163 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,5,1:6)=(/ 4647.600000,3446.400000,1715.200000,516.120000,86.053000,170.400000 /)
XPIZA_LKT(34,5,1:6)=(/ 0.983335,0.994055,0.999685,0.999383,0.992361,0.123141 /)
XCGA_LKT(34,5,1:6)=(/ 0.749497,0.726660,0.647123,0.469287,0.213333,0.084320 /)
XEXT_COEFF_550_LKT(34,5)=2715.200000 !rg=0.146163 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,6,1:6)=(/ 3968.700000,3304.000000,1932.600000,693.550000,136.100000,192.720000 /)
XPIZA_LKT(34,6,1:6)=(/ 0.980348,0.993631,0.999707,0.999503,0.994843,0.191698 /)
XCGA_LKT(34,6,1:6)=(/ 0.740983,0.733437,0.684390,0.555910,0.318060,0.133157 /)
XEXT_COEFF_550_LKT(34,6)=2776.400000 !rg=0.146163 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,7,1:6)=(/ 3314.200000,3008.800000,2048.600000,880.450000,205.870000,226.010000 /)
XPIZA_LKT(34,7,1:6)=(/ 0.976338,0.992835,0.999711,0.999579,0.996319,0.273921 /)
XCGA_LKT(34,7,1:6)=(/ 0.733107,0.732810,0.707907,0.619643,0.419513,0.202010 /)
XEXT_COEFF_550_LKT(34,7)=2675.100000 !rg=0.146163 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,8,1:6)=(/ 2736.100000,2642.000000,2032.400000,1047.200000,293.160000,271.960000 /)
XPIZA_LKT(34,8,1:6)=(/ 0.971241,0.991667,0.999698,0.999623,0.997197,0.357898 /)
XCGA_LKT(34,8,1:6)=(/ 0.727720,0.729170,0.719083,0.664570,0.506343,0.292400 /)
XEXT_COEFF_550_LKT(34,8)=2445.700000 !rg=0.146163 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,9,1:6)=(/ 2247.600000,2251.800000,1900.600000,1164.900000,395.060000,330.260000 /)
XPIZA_LKT(34,9,1:6)=(/ 0.965596,0.990179,0.999666,0.999642,0.997748,0.435719 /)
XCGA_LKT(34,9,1:6)=(/ 0.727187,0.726460,0.720987,0.694373,0.578167,0.385027 /)
XEXT_COEFF_550_LKT(34,9)=2157.400000 !rg=0.146163 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,10,1:6)=(/ 1830.800000,1882.600000,1707.800000,1205.900000,499.490000,393.100000 /)
XPIZA_LKT(34,10,1:6)=(/ 0.959583,0.987963,0.999618,0.999638,0.998083,0.499407 /)
XCGA_LKT(34,10,1:6)=(/ 0.731470,0.722207,0.720090,0.710373,0.632740,0.463403 /)
XEXT_COEFF_550_LKT(34,10)=1852.400000 !rg=0.146163 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,11,1:6)=(/ 1494.300000,1562.700000,1488.900000,1165.900000,592.660000,454.820000 /)
XPIZA_LKT(34,11,1:6)=(/ 0.950948,0.985252,0.999555,0.999610,0.998275,0.548314 /)
XCGA_LKT(34,11,1:6)=(/ 0.735797,0.722950,0.720020,0.715210,0.673010,0.533867 /)
XEXT_COEFF_550_LKT(34,11)=1550.400000 !rg=0.146163 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,12,1:6)=(/ 1213.100000,1275.600000,1267.800000,1068.600000,655.770000,503.180000 /)
XPIZA_LKT(34,12,1:6)=(/ 0.943291,0.982568,0.999438,0.999562,0.998352,0.582438 /)
XCGA_LKT(34,12,1:6)=(/ 0.746083,0.723047,0.717470,0.714540,0.699373,0.589100 /)
XEXT_COEFF_550_LKT(34,12)=1287.300000 !rg=0.146163 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,13,1:6)=(/ 988.270000,1039.900000,1054.900000,953.740000,673.530000,529.700000 /)
XPIZA_LKT(34,13,1:6)=(/ 0.933079,0.979321,0.999352,0.999492,0.998317,0.603091 /)
XCGA_LKT(34,13,1:6)=(/ 0.755133,0.731560,0.718043,0.713467,0.712253,0.631930 /)
XEXT_COEFF_550_LKT(34,13)=1052.400000 !rg=0.146163 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,14,1:6)=(/ 806.870000,845.670000,871.290000,830.640000,644.360000,529.110000 /)
XPIZA_LKT(34,14,1:6)=(/ 0.921349,0.975784,0.999173,0.999408,0.998170,0.611363 /)
XCGA_LKT(34,14,1:6)=(/ 0.764010,0.737800,0.721727,0.716293,0.713990,0.663123 /)
XEXT_COEFF_550_LKT(34,14)=860.460000 !rg=0.146163 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,15,1:6)=(/ 661.380000,690.780000,719.840000,704.730000,586.370000,503.750000 /)
XPIZA_LKT(34,15,1:6)=(/ 0.908387,0.970202,0.999057,0.999275,0.997915,0.609855 /)
XCGA_LKT(34,15,1:6)=(/ 0.776787,0.742913,0.724403,0.716973,0.709803,0.686110 /)
XEXT_COEFF_550_LKT(34,15)=707.690000 !rg=0.146163 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,16,1:6)=(/ 537.730000,564.870000,587.330000,589.950000,522.200000,461.880000 /)
XPIZA_LKT(34,16,1:6)=(/ 0.893130,0.964266,0.998879,0.999175,0.997611,0.602592 /)
XCGA_LKT(34,16,1:6)=(/ 0.787823,0.752653,0.726153,0.720210,0.710777,0.705577 /)
XEXT_COEFF_550_LKT(34,16)=572.700000 !rg=0.146163 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,17,1:6)=(/ 440.180000,457.830000,476.160000,489.620000,458.450000,413.010000 /)
XPIZA_LKT(34,17,1:6)=(/ 0.879836,0.959035,0.998713,0.998969,0.997244,0.593727 /)
XCGA_LKT(34,17,1:6)=(/ 0.799357,0.763777,0.738923,0.719730,0.714910,0.724873 /)
XEXT_COEFF_550_LKT(34,17)=465.230000 !rg=0.146163 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,18,1:6)=(/ 362.290000,375.250000,392.290000,405.980000,390.770000,362.360000 /)
XPIZA_LKT(34,18,1:6)=(/ 0.866038,0.950780,0.998397,0.998712,0.996713,0.584874 /)
XCGA_LKT(34,18,1:6)=(/ 0.810957,0.771213,0.743737,0.726387,0.715197,0.744293 /)
XEXT_COEFF_550_LKT(34,18)=383.360000 !rg=0.146163 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,19,1:6)=(/ 296.390000,307.620000,319.490000,331.550000,330.570000,312.160000 /)
XPIZA_LKT(34,19,1:6)=(/ 0.853587,0.943203,0.998033,0.998559,0.995933,0.575562 /)
XCGA_LKT(34,19,1:6)=(/ 0.820167,0.781640,0.747377,0.733873,0.717053,0.762633 /)
XEXT_COEFF_550_LKT(34,19)=310.730000 !rg=0.146163 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,20,1:6)=(/ 243.810000,251.140000,258.930000,272.930000,274.820000,265.520000 /)
XPIZA_LKT(34,20,1:6)=(/ 0.843219,0.936917,0.997743,0.998030,0.995470,0.566689 /)
XCGA_LKT(34,20,1:6)=(/ 0.829280,0.792453,0.761703,0.735360,0.725010,0.780190 /)
XEXT_COEFF_550_LKT(34,20)=254.140000 !rg=0.146163 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,1,1:6)=(/ 5909.400000,3161.900000,969.750000,196.160000,22.087000,144.490000 /)
XPIZA_LKT(35,1,1:6)=(/ 0.987517,0.993968,0.999534,0.998637,0.973347,0.031467 /)
XCGA_LKT(35,1,1:6)=(/ 0.756620,0.646007,0.429213,0.146490,0.047207,0.017720 /)
XEXT_COEFF_550_LKT(35,1)=2057.400000 !rg=0.158347 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,2,1:6)=(/ 5870.200000,3245.300000,1101.900000,235.850000,27.606000,146.720000 /)
XPIZA_LKT(35,2,1:6)=(/ 0.987147,0.994024,0.999568,0.998841,0.978375,0.039143 /)
XCGA_LKT(35,2,1:6)=(/ 0.762037,0.673753,0.500220,0.186397,0.059843,0.022510 /)
XEXT_COEFF_550_LKT(35,2)=2227.500000 !rg=0.158347 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,3,1:6)=(/ 5599.200000,3429.000000,1338.800000,312.660000,40.588000,151.780000 /)
XPIZA_LKT(35,3,1:6)=(/ 0.986382,0.994213,0.999623,0.999086,0.984883,0.057100 /)
XCGA_LKT(35,3,1:6)=(/ 0.763153,0.702600,0.567190,0.276163,0.090597,0.034337 /)
XEXT_COEFF_550_LKT(35,3)=2460.500000 !rg=0.158347 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,4,1:6)=(/ 5072.000000,3538.400000,1606.300000,432.840000,65.046000,161.480000 /)
XPIZA_LKT(35,4,1:6)=(/ 0.984782,0.994257,0.999670,0.999293,0.990174,0.091631 /)
XCGA_LKT(35,4,1:6)=(/ 0.758217,0.723490,0.622567,0.398237,0.151493,0.058703 /)
XEXT_COEFF_550_LKT(35,4)=2690.400000 !rg=0.158347 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,5,1:6)=(/ 4382.700000,3479.800000,1854.200000,594.500000,104.710000,178.390000 /)
XPIZA_LKT(35,5,1:6)=(/ 0.982230,0.994014,0.999700,0.999443,0.993549,0.148077 /)
XCGA_LKT(35,5,1:6)=(/ 0.748660,0.734790,0.666860,0.504627,0.245840,0.098933 /)
XEXT_COEFF_550_LKT(35,5)=2822.300000 !rg=0.158347 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,6,1:6)=(/ 3688.100000,3242.300000,2026.400000,777.830000,162.470000,204.800000 /)
XPIZA_LKT(35,6,1:6)=(/ 0.978716,0.993433,0.999714,0.999541,0.995541,0.223343 /)
XCGA_LKT(35,6,1:6)=(/ 0.738433,0.737247,0.697557,0.582910,0.354537,0.156330 /)
XEXT_COEFF_550_LKT(35,6)=2800.100000 !rg=0.158347 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,7,1:6)=(/ 3051.700000,2886.000000,2083.300000,960.210000,239.820000,243.250000 /)
XPIZA_LKT(35,7,1:6)=(/ 0.973993,0.992471,0.999710,0.999602,0.996732,0.308187 /)
XCGA_LKT(35,7,1:6)=(/ 0.730410,0.734143,0.715500,0.639413,0.452747,0.235583 /)
XEXT_COEFF_550_LKT(35,7)=2625.900000 !rg=0.158347 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,8,1:6)=(/ 2502.800000,2490.000000,2010.400000,1109.700000,334.410000,294.700000 /)
XPIZA_LKT(35,8,1:6)=(/ 0.969086,0.991143,0.999689,0.999634,0.997455,0.390788 /)
XCGA_LKT(35,8,1:6)=(/ 0.727917,0.729587,0.722203,0.678283,0.536300,0.330000 /)
XEXT_COEFF_550_LKT(35,8)=2352.700000 !rg=0.158347 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,9,1:6)=(/ 2045.500000,2097.900000,1839.100000,1198.300000,438.710000,355.900000 /)
XPIZA_LKT(35,9,1:6)=(/ 0.962502,0.989339,0.999651,0.999644,0.997905,0.463724 /)
XCGA_LKT(35,9,1:6)=(/ 0.726847,0.724567,0.722557,0.702637,0.600843,0.416610 /)
XEXT_COEFF_550_LKT(35,9)=2041.800000 !rg=0.158347 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,10,1:6)=(/ 1673.400000,1738.900000,1627.900000,1204.400000,540.230000,419.280000 /)
XPIZA_LKT(35,10,1:6)=(/ 0.955752,0.987320,0.999596,0.999630,0.998177,0.520980 /)
XCGA_LKT(35,10,1:6)=(/ 0.732053,0.724263,0.720843,0.713880,0.650097,0.492473 /)
XEXT_COEFF_550_LKT(35,10)=1727.200000 !rg=0.158347 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,11,1:6)=(/ 1358.000000,1430.800000,1402.200000,1135.700000,623.390000,476.710000 /)
XPIZA_LKT(35,11,1:6)=(/ 0.947981,0.984423,0.999524,0.999595,0.998321,0.563957 /)
XCGA_LKT(35,11,1:6)=(/ 0.739847,0.723140,0.720457,0.716147,0.684730,0.556853 /)
XEXT_COEFF_550_LKT(35,11)=1434.600000 !rg=0.158347 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,12,1:6)=(/ 1108.200000,1167.000000,1174.800000,1027.800000,669.600000,517.490000 /)
XPIZA_LKT(35,12,1:6)=(/ 0.938577,0.981635,0.999423,0.999536,0.998351,0.592485 /)
XCGA_LKT(35,12,1:6)=(/ 0.748903,0.728720,0.718123,0.713987,0.705677,0.607377 /)
XEXT_COEFF_550_LKT(35,12)=1182.100000 !rg=0.158347 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,13,1:6)=(/ 901.920000,948.770000,975.530000,906.560000,667.630000,533.340000 /)
XPIZA_LKT(35,13,1:6)=(/ 0.927866,0.978010,0.999314,0.999464,0.998273,0.607919 /)
XCGA_LKT(35,13,1:6)=(/ 0.758070,0.733317,0.721400,0.716147,0.714067,0.645407 /)
XEXT_COEFF_550_LKT(35,13)=966.750000 !rg=0.158347 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,14,1:6)=(/ 740.280000,772.920000,806.340000,781.450000,623.600000,521.970000 /)
XPIZA_LKT(35,14,1:6)=(/ 0.915385,0.973762,0.999124,0.999343,0.998081,0.611939 /)
XCGA_LKT(35,14,1:6)=(/ 0.770540,0.739223,0.720660,0.716120,0.712533,0.673017 /)
XEXT_COEFF_550_LKT(35,14)=793.670000 !rg=0.158347 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,15,1:6)=(/ 602.340000,634.100000,660.670000,654.920000,560.990000,488.670000 /)
XPIZA_LKT(35,15,1:6)=(/ 0.902201,0.967946,0.998896,0.999246,0.997805,0.607655 /)
XCGA_LKT(35,15,1:6)=(/ 0.780323,0.748230,0.722503,0.717823,0.710367,0.694140 /)
XEXT_COEFF_550_LKT(35,15)=644.980000 !rg=0.158347 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,16,1:6)=(/ 490.600000,513.630000,537.520000,546.410000,498.270000,442.820000 /)
XPIZA_LKT(35,16,1:6)=(/ 0.887439,0.962381,0.998758,0.999098,0.997472,0.599390 /)
XCGA_LKT(35,16,1:6)=(/ 0.794190,0.756180,0.733117,0.720210,0.712943,0.713560 /)
XEXT_COEFF_550_LKT(35,16)=520.450000 !rg=0.158347 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,17,1:6)=(/ 404.630000,418.720000,438.860000,454.680000,430.520000,392.610000 /)
XPIZA_LKT(35,17,1:6)=(/ 0.872739,0.956394,0.998565,0.998790,0.997040,0.590508 /)
XCGA_LKT(35,17,1:6)=(/ 0.805493,0.767053,0.739440,0.722580,0.715420,0.733230 /)
XEXT_COEFF_550_LKT(35,17)=428.390000 !rg=0.158347 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,18,1:6)=(/ 330.780000,344.550000,359.370000,371.440000,364.430000,341.640000 /)
XPIZA_LKT(35,18,1:6)=(/ 0.860543,0.947645,0.998116,0.998635,0.996415,0.581376 /)
XCGA_LKT(35,18,1:6)=(/ 0.814247,0.777353,0.742177,0.729397,0.715657,0.752237 /)
XEXT_COEFF_550_LKT(35,18)=349.690000 !rg=0.158347 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,19,1:6)=(/ 271.210000,280.380000,292.000000,304.720000,307.990000,292.340000 /)
XPIZA_LKT(35,19,1:6)=(/ 0.848714,0.940551,0.997783,0.998312,0.995495,0.571980 /)
XCGA_LKT(35,19,1:6)=(/ 0.825297,0.785247,0.754580,0.733757,0.720950,0.770233 /)
XEXT_COEFF_550_LKT(35,19)=282.960000 !rg=0.158347 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,20,1:6)=(/ 223.920000,230.070000,238.940000,251.120000,254.570000,247.420000 /)
XPIZA_LKT(35,20,1:6)=(/ 0.837387,0.934758,0.997449,0.997725,0.995030,0.563332 /)
XCGA_LKT(35,20,1:6)=(/ 0.834060,0.795297,0.763103,0.739967,0.726090,0.787643 /)
XEXT_COEFF_550_LKT(35,20)=233.410000 !rg=0.158347 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,1,1:6)=(/ 5987.800000,3319.900000,1119.500000,243.140000,27.895000,146.970000 /)
XPIZA_LKT(36,1,1:6)=(/ 0.987272,0.994293,0.999571,0.998866,0.978532,0.039393 /)
XCGA_LKT(36,1,1:6)=(/ 0.764740,0.664183,0.507110,0.173007,0.055333,0.020780 /)
XEXT_COEFF_550_LKT(36,1)=2390.500000 !rg=0.171546 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,2,1:6)=(/ 5862.000000,3475.600000,1282.500000,287.850000,34.840000,149.680000 /)
XPIZA_LKT(36,2,1:6)=(/ 0.986975,0.994309,0.999607,0.999019,0.982540,0.048858 /)
XCGA_LKT(36,2,1:6)=(/ 0.770067,0.696593,0.556167,0.220347,0.070140,0.026400 /)
XEXT_COEFF_550_LKT(36,2)=2476.900000 !rg=0.171546 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,3,1:6)=(/ 5457.900000,3610.200000,1519.100000,373.510000,50.927000,155.880000 /)
XPIZA_LKT(36,3,1:6)=(/ 0.985878,0.994397,0.999655,0.999206,0.987691,0.070767 /)
XCGA_LKT(36,3,1:6)=(/ 0.767207,0.719837,0.601817,0.321160,0.106097,0.040273 /)
XEXT_COEFF_550_LKT(36,3)=2671.400000 !rg=0.171546 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,4,1:6)=(/ 4814.500000,3628.200000,1768.600000,506.890000,80.274000,167.680000 /)
XPIZA_LKT(36,4,1:6)=(/ 0.983836,0.994302,0.999691,0.999372,0.991839,0.111871 /)
XCGA_LKT(36,4,1:6)=(/ 0.757887,0.734940,0.647197,0.441330,0.176433,0.068843 /)
XEXT_COEFF_550_LKT(36,4)=2843.600000 !rg=0.171546 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,5,1:6)=(/ 4081.700000,3464.900000,1980.500000,678.620000,126.610000,187.950000 /)
XPIZA_LKT(36,5,1:6)=(/ 0.980736,0.993904,0.999712,0.999493,0.994506,0.176103 /)
XCGA_LKT(36,5,1:6)=(/ 0.745397,0.741000,0.684150,0.537910,0.280557,0.116040 /)
XEXT_COEFF_550_LKT(36,5)=2895.300000 !rg=0.171546 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,6,1:6)=(/ 3387.200000,3146.100000,2099.300000,863.230000,192.160000,218.910000 /)
XPIZA_LKT(36,6,1:6)=(/ 0.976745,0.993139,0.999717,0.999572,0.996104,0.256705 /)
XCGA_LKT(36,6,1:6)=(/ 0.734327,0.739350,0.708663,0.607467,0.390077,0.183090 /)
XEXT_COEFF_550_LKT(36,6)=2788.300000 !rg=0.171546 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,7,1:6)=(/ 2780.300000,2746.100000,2093.400000,1035.600000,277.300000,262.870000 /)
XPIZA_LKT(36,7,1:6)=(/ 0.971753,0.991988,0.999707,0.999620,0.997072,0.342440 /)
XCGA_LKT(36,7,1:6)=(/ 0.727543,0.733897,0.721490,0.656983,0.485727,0.271397 /)
XEXT_COEFF_550_LKT(36,7)=2547.600000 !rg=0.171546 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,8,1:6)=(/ 2280.600000,2332.900000,1967.100000,1162.200000,377.180000,319.070000 /)
XPIZA_LKT(36,8,1:6)=(/ 0.966161,0.990490,0.999678,0.999643,0.997668,0.422139 /)
XCGA_LKT(36,8,1:6)=(/ 0.725850,0.728823,0.724557,0.689940,0.562933,0.364767 /)
XEXT_COEFF_550_LKT(36,8)=2241.300000 !rg=0.171546 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,9,1:6)=(/ 1868.600000,1939.700000,1769.800000,1217.600000,482.930000,382.100000 /)
XPIZA_LKT(36,9,1:6)=(/ 0.958588,0.988486,0.999631,0.999643,0.998036,0.489027 /)
XCGA_LKT(36,9,1:6)=(/ 0.728177,0.723913,0.722660,0.709003,0.621963,0.446720 /)
XEXT_COEFF_550_LKT(36,9)=1921.400000 !rg=0.171546 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,10,1:6)=(/ 1521.100000,1603.100000,1544.200000,1189.900000,578.730000,444.590000 /)
XPIZA_LKT(36,10,1:6)=(/ 0.951695,0.986216,0.999570,0.999620,0.998249,0.540291 /)
XCGA_LKT(36,10,1:6)=(/ 0.733083,0.723103,0.721807,0.716250,0.665283,0.519577 /)
XEXT_COEFF_550_LKT(36,10)=1603.000000 !rg=0.171546 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,11,1:6)=(/ 1244.200000,1308.000000,1308.800000,1100.700000,648.780000,496.180000 /)
XPIZA_LKT(36,11,1:6)=(/ 0.943133,0.983477,0.999488,0.999574,0.998346,0.577119 /)
XCGA_LKT(36,11,1:6)=(/ 0.742733,0.725657,0.719703,0.715787,0.694600,0.577803 /)
XEXT_COEFF_550_LKT(36,11)=1327.800000 !rg=0.171546 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,12,1:6)=(/ 1010.700000,1068.100000,1089.700000,982.680000,675.650000,527.540000 /)
XPIZA_LKT(36,12,1:6)=(/ 0.933357,0.980055,0.999389,0.999512,0.998334,0.600314 /)
XCGA_LKT(36,12,1:6)=(/ 0.751453,0.729567,0.720097,0.715817,0.710353,0.623427 /)
XEXT_COEFF_550_LKT(36,12)=1084.800000 !rg=0.171546 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,13,1:6)=(/ 827.750000,867.560000,903.370000,859.170000,654.950000,532.380000 /)
XPIZA_LKT(36,13,1:6)=(/ 0.921859,0.976292,0.999234,0.999428,0.998208,0.610998 /)
XCGA_LKT(36,13,1:6)=(/ 0.764203,0.735010,0.719817,0.717107,0.714217,0.657233 /)
XEXT_COEFF_550_LKT(36,13)=890.780000 !rg=0.171546 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,14,1:6)=(/ 674.040000,711.790000,740.050000,726.960000,599.930000,511.170000 /)
XPIZA_LKT(36,14,1:6)=(/ 0.909774,0.970758,0.999079,0.999318,0.997987,0.611257 /)
XCGA_LKT(36,14,1:6)=(/ 0.773120,0.743113,0.721170,0.716803,0.712303,0.681960 /)
XEXT_COEFF_550_LKT(36,14)=724.830000 !rg=0.171546 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,15,1:6)=(/ 548.620000,577.310000,607.650000,609.890000,537.580000,471.610000 /)
XPIZA_LKT(36,15,1:6)=(/ 0.896126,0.966070,0.998827,0.999183,0.997671,0.604914 /)
XCGA_LKT(36,15,1:6)=(/ 0.788747,0.749510,0.726693,0.718970,0.710927,0.702097 /)
XEXT_COEFF_550_LKT(36,15)=588.590000 !rg=0.171546 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,16,1:6)=(/ 451.000000,467.450000,490.660000,507.750000,471.690000,423.020000 /)
XPIZA_LKT(36,16,1:6)=(/ 0.880925,0.960168,0.998719,0.998963,0.997303,0.596193 /)
XCGA_LKT(36,16,1:6)=(/ 0.799757,0.763260,0.736480,0.720237,0.714497,0.721833 /)
XEXT_COEFF_550_LKT(36,16)=479.010000 !rg=0.171546 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,17,1:6)=(/ 370.090000,385.020000,402.410000,416.830000,402.320000,371.890000 /)
XPIZA_LKT(36,17,1:6)=(/ 0.866970,0.952252,0.998341,0.998814,0.996843,0.587191 /)
XCGA_LKT(36,17,1:6)=(/ 0.808427,0.772127,0.739250,0.725233,0.715787,0.741523 /)
XEXT_COEFF_550_LKT(36,17)=392.130000 !rg=0.171546 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,18,1:6)=(/ 302.670000,314.200000,329.870000,341.700000,341.790000,321.070000 /)
XPIZA_LKT(36,18,1:6)=(/ 0.854796,0.945088,0.998025,0.998559,0.996049,0.577732 /)
XCGA_LKT(36,18,1:6)=(/ 0.821387,0.779610,0.747720,0.731250,0.718463,0.760007 /)
XEXT_COEFF_550_LKT(36,18)=319.190000 !rg=0.171546 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,19,1:6)=(/ 249.570000,255.980000,266.920000,281.720000,283.840000,273.140000 /)
XPIZA_LKT(36,19,1:6)=(/ 0.843859,0.938045,0.997778,0.998164,0.995103,0.568452 /)
XCGA_LKT(36,19,1:6)=(/ 0.829993,0.792020,0.760210,0.736923,0.722323,0.777833 /)
XEXT_COEFF_550_LKT(36,19)=261.320000 !rg=0.171546 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,20,1:6)=(/ 205.480000,211.320000,219.180000,228.710000,235.160000,230.070000 /)
XPIZA_LKT(36,20,1:6)=(/ 0.832738,0.931539,0.997233,0.998002,0.994761,0.560009 /)
XCGA_LKT(36,20,1:6)=(/ 0.836923,0.799837,0.763437,0.743673,0.727353,0.795017 /)
XEXT_COEFF_550_LKT(36,20)=214.190000 !rg=0.171546 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,1,1:6)=(/ 6044.700000,3530.700000,1308.300000,298.310000,35.238000,149.980000 /)
XPIZA_LKT(37,1,1:6)=(/ 0.987245,0.994386,0.999605,0.999045,0.982674,0.049156 /)
XCGA_LKT(37,1,1:6)=(/ 0.778383,0.702590,0.578643,0.205040,0.064867,0.024367 /)
XEXT_COEFF_550_LKT(37,1)=2632.800000 !rg=0.185845 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,2,1:6)=(/ 5742.800000,3683.200000,1484.000000,347.130000,43.934000,153.300000 /)
XPIZA_LKT(37,2,1:6)=(/ 0.986560,0.994509,0.999643,0.999159,0.985861,0.060748 /)
XCGA_LKT(37,2,1:6)=(/ 0.775273,0.719553,0.598520,0.261087,0.082230,0.030960 /)
XEXT_COEFF_550_LKT(37,2)=2697.900000 !rg=0.185845 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,3,1:6)=(/ 5215.400000,3747.800000,1701.500000,442.400000,63.673000,160.870000 /)
XPIZA_LKT(37,3,1:6)=(/ 0.985083,0.994506,0.999680,0.999303,0.989921,0.087232 /)
XCGA_LKT(37,3,1:6)=(/ 0.767990,0.734830,0.630560,0.370407,0.124263,0.047240 /)
XEXT_COEFF_550_LKT(37,3)=2859.600000 !rg=0.185845 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,4,1:6)=(/ 4491.500000,3665.900000,1923.100000,588.650000,98.446000,175.190000 /)
XPIZA_LKT(37,4,1:6)=(/ 0.982544,0.994273,0.999707,0.999436,0.993164,0.135428 /)
XCGA_LKT(37,4,1:6)=(/ 0.754953,0.744070,0.668787,0.483170,0.204910,0.080727 /)
XEXT_COEFF_550_LKT(37,4)=2964.800000 !rg=0.185845 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,5,1:6)=(/ 3743.600000,3402.900000,2089.600000,767.160000,151.950000,199.320000 /)
XPIZA_LKT(37,5,1:6)=(/ 0.978949,0.993703,0.999721,0.999534,0.995277,0.206903 /)
XCGA_LKT(37,5,1:6)=(/ 0.740663,0.745030,0.699190,0.568713,0.316677,0.135983 /)
XEXT_COEFF_550_LKT(37,5)=2930.200000 !rg=0.185845 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,6,1:6)=(/ 3084.900000,3013.900000,2148.100000,948.420000,225.330000,235.280000 /)
XPIZA_LKT(37,6,1:6)=(/ 0.974490,0.992768,0.999718,0.999598,0.996561,0.291134 /)
XCGA_LKT(37,6,1:6)=(/ 0.730573,0.739927,0.717917,0.629710,0.425313,0.213227 /)
XEXT_COEFF_550_LKT(37,6)=2742.900000 !rg=0.185845 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,7,1:6)=(/ 2526.000000,2581.600000,2080.800000,1105.000000,318.030000,284.790000 /)
XPIZA_LKT(37,7,1:6)=(/ 0.969230,0.991430,0.999700,0.999634,0.997355,0.376094 /)
XCGA_LKT(37,7,1:6)=(/ 0.725857,0.733043,0.725920,0.672497,0.517113,0.307223 /)
XEXT_COEFF_550_LKT(37,7)=2450.100000 !rg=0.185845 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,8,1:6)=(/ 2068.200000,2165.000000,1909.500000,1203.700000,421.260000,344.560000 /)
XPIZA_LKT(37,8,1:6)=(/ 0.962963,0.989695,0.999664,0.999647,0.997841,0.451274 /)
XCGA_LKT(37,8,1:6)=(/ 0.725107,0.726640,0.726407,0.699703,0.587143,0.396830 /)
XEXT_COEFF_550_LKT(37,8)=2119.400000 !rg=0.185845 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,9,1:6)=(/ 1691.700000,1794.700000,1687.300000,1223.700000,525.330000,408.740000 /)
XPIZA_LKT(37,9,1:6)=(/ 0.955524,0.987216,0.999611,0.999637,0.998142,0.511819 /)
XCGA_LKT(37,9,1:6)=(/ 0.728967,0.722670,0.724083,0.713777,0.640633,0.476297 /)
XEXT_COEFF_550_LKT(37,9)=1786.500000 !rg=0.185845 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,10,1:6)=(/ 1390.000000,1469.600000,1455.500000,1166.300000,612.500000,467.890000 /)
XPIZA_LKT(37,10,1:6)=(/ 0.947038,0.984932,0.999534,0.999605,0.998305,0.557151 /)
XCGA_LKT(37,10,1:6)=(/ 0.737020,0.722687,0.720590,0.717397,0.678240,0.543647 /)
XEXT_COEFF_550_LKT(37,10)=1488.600000 !rg=0.185845 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,11,1:6)=(/ 1130.800000,1202.000000,1216.900000,1058.600000,666.580000,512.580000 /)
XPIZA_LKT(37,11,1:6)=(/ 0.938168,0.981692,0.999450,0.999553,0.998356,0.588241 /)
XCGA_LKT(37,11,1:6)=(/ 0.744260,0.725140,0.719473,0.716727,0.702103,0.597090 /)
XEXT_COEFF_550_LKT(37,11)=1220.800000 !rg=0.185845 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,12,1:6)=(/ 925.890000,974.370000,1013.800000,937.960000,674.130000,533.550000 /)
XPIZA_LKT(37,12,1:6)=(/ 0.927824,0.978401,0.999301,0.999483,0.998298,0.606128 /)
XCGA_LKT(37,12,1:6)=(/ 0.757713,0.730323,0.718633,0.717153,0.713043,0.637833 /)
XEXT_COEFF_550_LKT(37,12)=1000.700000 !rg=0.185845 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,13,1:6)=(/ 753.320000,797.870000,829.220000,805.850000,636.690000,527.330000 /)
XPIZA_LKT(37,13,1:6)=(/ 0.916556,0.973408,0.999153,0.999387,0.998133,0.612398 /)
XCGA_LKT(37,13,1:6)=(/ 0.766400,0.738173,0.719130,0.717547,0.714203,0.667883 /)
XEXT_COEFF_550_LKT(37,13)=813.040000 !rg=0.185845 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,14,1:6)=(/ 613.960000,647.810000,683.580000,676.700000,576.690000,497.450000 /)
XPIZA_LKT(37,14,1:6)=(/ 0.903870,0.969055,0.998762,0.999278,0.997867,0.609602 /)
XCGA_LKT(37,14,1:6)=(/ 0.782297,0.743470,0.723040,0.718607,0.711253,0.690447 /)
XEXT_COEFF_550_LKT(37,14)=663.200000 !rg=0.185845 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,15,1:6)=(/ 504.270000,524.970000,555.240000,568.360000,511.880000,453.130000 /)
XPIZA_LKT(37,15,1:6)=(/ 0.888910,0.963908,0.998397,0.999080,0.997542,0.601928 /)
XCGA_LKT(37,15,1:6)=(/ 0.793803,0.758407,0.728487,0.719463,0.712880,0.710217 /)
XEXT_COEFF_550_LKT(37,15)=536.580000 !rg=0.185845 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,16,1:6)=(/ 413.020000,430.690000,449.920000,467.340000,444.090000,402.670000 /)
XPIZA_LKT(37,16,1:6)=(/ 0.873985,0.956611,0.998522,0.998932,0.997137,0.592990 /)
XCGA_LKT(37,16,1:6)=(/ 0.802720,0.766417,0.737633,0.720143,0.715790,0.730280 /)
XEXT_COEFF_550_LKT(37,16)=439.190000 !rg=0.185845 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,17,1:6)=(/ 338.750000,352.080000,369.290000,381.550000,377.330000,351.000000 /)
XPIZA_LKT(37,17,1:6)=(/ 0.861009,0.948608,0.998272,0.998763,0.996437,0.583681 /)
XCGA_LKT(37,17,1:6)=(/ 0.815633,0.773600,0.745263,0.730273,0.715863,0.749630 /)
XEXT_COEFF_550_LKT(37,17)=360.020000 !rg=0.185845 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,18,1:6)=(/ 278.600000,286.700000,301.160000,315.670000,316.740000,300.880000 /)
XPIZA_LKT(37,18,1:6)=(/ 0.849686,0.941664,0.997649,0.998406,0.995920,0.574038 /)
XCGA_LKT(37,18,1:6)=(/ 0.825613,0.787670,0.749620,0.733717,0.719800,0.767710 /)
XEXT_COEFF_550_LKT(37,18)=292.020000 !rg=0.185845 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,19,1:6)=(/ 229.080000,236.060000,244.990000,257.290000,262.100000,254.620000 /)
XPIZA_LKT(37,19,1:6)=(/ 0.838241,0.934784,0.997477,0.998042,0.995121,0.564963 /)
XCGA_LKT(37,19,1:6)=(/ 0.832657,0.795443,0.760867,0.737087,0.724387,0.785400 /)
XEXT_COEFF_550_LKT(37,19)=239.600000 !rg=0.185845 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,20,1:6)=(/ 188.520000,194.130000,201.120000,208.510000,218.140000,213.540000 /)
XPIZA_LKT(37,20,1:6)=(/ 0.827877,0.928355,0.997156,0.997812,0.994069,0.556766 /)
XCGA_LKT(37,20,1:6)=(/ 0.841953,0.801657,0.769767,0.750683,0.729017,0.802330 /)
XEXT_COEFF_550_LKT(37,20)=196.800000 !rg=0.185845 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,1,1:6)=(/ 5803.400000,3841.000000,1550.200000,361.230000,44.496000,153.660000 /)
XPIZA_LKT(38,1,1:6)=(/ 0.986762,0.994552,0.999643,0.999186,0.985980,0.061104 /)
XCGA_LKT(38,1,1:6)=(/ 0.783553,0.738717,0.625667,0.244077,0.076067,0.028573 /)
XEXT_COEFF_550_LKT(38,1)=2772.800000 !rg=0.201336 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,2,1:6)=(/ 5505.500000,3854.400000,1692.600000,413.410000,55.308000,157.720000 /)
XPIZA_LKT(38,2,1:6)=(/ 0.985844,0.994646,0.999674,0.999268,0.988506,0.075185 /)
XCGA_LKT(38,2,1:6)=(/ 0.776913,0.738680,0.627057,0.309663,0.096447,0.036310 /)
XEXT_COEFF_550_LKT(38,2)=2896.000000 !rg=0.201336 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,3,1:6)=(/ 4879.200000,3833.100000,1879.000000,520.040000,79.223000,166.970000 /)
XPIZA_LKT(38,3,1:6)=(/ 0.983918,0.994541,0.999702,0.999382,0.991691,0.106841 /)
XCGA_LKT(38,3,1:6)=(/ 0.765190,0.747143,0.654943,0.421973,0.145587,0.055420 /)
XEXT_COEFF_550_LKT(38,3)=3018.700000 !rg=0.201336 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,4,1:6)=(/ 4118.900000,3647.300000,2064.400000,677.590000,119.890000,184.250000 /)
XPIZA_LKT(38,4,1:6)=(/ 0.980802,0.994162,0.999720,0.999489,0.994223,0.162352 /)
XCGA_LKT(38,4,1:6)=(/ 0.748700,0.750773,0.687727,0.522577,0.236983,0.094650 /)
XEXT_COEFF_550_LKT(38,4)=3047.900000 !rg=0.201336 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,5,1:6)=(/ 3398.200000,3292.600000,2176.900000,858.430000,180.920000,212.740000 /)
XPIZA_LKT(38,5,1:6)=(/ 0.976748,0.993424,0.999726,0.999569,0.995899,0.239980 /)
XCGA_LKT(38,5,1:6)=(/ 0.734870,0.747180,0.712047,0.596767,0.353830,0.159027 /)
XEXT_COEFF_550_LKT(38,5)=2925.100000 !rg=0.201336 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,6,1:6)=(/ 2798.100000,2854.700000,2170.900000,1030.600000,262.260000,254.050000 /)
XPIZA_LKT(38,6,1:6)=(/ 0.971574,0.992295,0.999716,0.999618,0.996937,0.325983 /)
XCGA_LKT(38,6,1:6)=(/ 0.724980,0.739567,0.725263,0.649593,0.460310,0.245780 /)
XEXT_COEFF_550_LKT(38,6)=2666.400000 !rg=0.201336 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,7,1:6)=(/ 2292.200000,2407.800000,2044.500000,1165.400000,360.900000,308.640000 /)
XPIZA_LKT(38,7,1:6)=(/ 0.965781,0.990773,0.999690,0.999644,0.997587,0.408428 /)
XCGA_LKT(38,7,1:6)=(/ 0.722377,0.731570,0.728950,0.685910,0.545623,0.341410 /)
XEXT_COEFF_550_LKT(38,7)=2332.600000 !rg=0.201336 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,8,1:6)=(/ 1882.400000,1996.500000,1838.700000,1231.600000,466.370000,370.950000 /)
XPIZA_LKT(38,8,1:6)=(/ 0.958815,0.988831,0.999645,0.999647,0.997986,0.477874 /)
XCGA_LKT(38,8,1:6)=(/ 0.725073,0.725137,0.726747,0.707600,0.609717,0.427587 /)
XEXT_COEFF_550_LKT(38,8)=1989.700000 !rg=0.201336 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,9,1:6)=(/ 1531.100000,1643.100000,1601.100000,1216.500000,565.920000,434.970000 /)
XPIZA_LKT(38,9,1:6)=(/ 0.952587,0.986253,0.999586,0.999629,0.998223,0.532328 /)
XCGA_LKT(38,9,1:6)=(/ 0.733153,0.721657,0.724317,0.717300,0.657100,0.504247 /)
XEXT_COEFF_550_LKT(38,9)=1654.700000 !rg=0.201336 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,10,1:6)=(/ 1259.500000,1350.300000,1356.900000,1132.700000,641.300000,489.000000 /)
XPIZA_LKT(38,10,1:6)=(/ 0.943163,0.983075,0.999503,0.999588,0.998339,0.571511 /)
XCGA_LKT(38,10,1:6)=(/ 0.738787,0.722437,0.720883,0.718107,0.689270,0.565563 /)
XEXT_COEFF_550_LKT(38,10)=1369.500000 !rg=0.201336 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,11,1:6)=(/ 1033.400000,1096.300000,1132.800000,1016.000000,677.090000,524.920000 /)
XPIZA_LKT(38,11,1:6)=(/ 0.932867,0.979975,0.999393,0.999529,0.998348,0.597224 /)
XCGA_LKT(38,11,1:6)=(/ 0.751113,0.725430,0.718577,0.717830,0.708007,0.614203 /)
XEXT_COEFF_550_LKT(38,11)=1124.800000 !rg=0.201336 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,12,1:6)=(/ 841.650000,896.120000,931.910000,887.450000,664.880000,535.000000 /)
XPIZA_LKT(38,12,1:6)=(/ 0.922902,0.975832,0.999261,0.999448,0.998247,0.610238 /)
XCGA_LKT(38,12,1:6)=(/ 0.760087,0.733280,0.718550,0.718027,0.714550,0.650600 /)
XEXT_COEFF_550_LKT(38,12)=914.650000 !rg=0.201336 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,13,1:6)=(/ 686.070000,726.390000,765.720000,752.250000,615.940000,518.380000 /)
XPIZA_LKT(38,13,1:6)=(/ 0.910937,0.971710,0.999076,0.999351,0.998032,0.612509 /)
XCGA_LKT(38,13,1:6)=(/ 0.775783,0.738260,0.721550,0.718450,0.712867,0.677497 /)
XEXT_COEFF_550_LKT(38,13)=746.370000 !rg=0.201336 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,14,1:6)=(/ 562.320000,590.480000,621.580000,632.640000,551.260000,481.370000 /)
XPIZA_LKT(38,14,1:6)=(/ 0.897215,0.966423,0.998923,0.999192,0.997755,0.607276 /)
XCGA_LKT(38,14,1:6)=(/ 0.787743,0.752073,0.725957,0.718157,0.712280,0.698707 /)
XEXT_COEFF_550_LKT(38,14)=603.020000 !rg=0.201336 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,15,1:6)=(/ 462.160000,481.310000,505.530000,523.410000,486.470000,433.660000 /)
XPIZA_LKT(38,15,1:6)=(/ 0.881559,0.961037,0.998724,0.999025,0.997399,0.598844 /)
XCGA_LKT(38,15,1:6)=(/ 0.796967,0.762003,0.733673,0.718367,0.715123,0.718583 /)
XEXT_COEFF_550_LKT(38,15)=491.360000 !rg=0.201336 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,16,1:6)=(/ 379.020000,393.830000,414.940000,429.100000,416.600000,381.900000 /)
XPIZA_LKT(38,16,1:6)=(/ 0.867567,0.952963,0.998300,0.998890,0.996893,0.589678 /)
XCGA_LKT(38,16,1:6)=(/ 0.809713,0.768297,0.739717,0.727853,0.715407,0.738710 /)
XEXT_COEFF_550_LKT(38,16)=402.310000 !rg=0.201336 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,17,1:6)=(/ 310.640000,322.150000,336.570000,353.140000,351.410000,330.140000 /)
XPIZA_LKT(38,17,1:6)=(/ 0.855283,0.945138,0.998168,0.998611,0.996218,0.579976 /)
XCGA_LKT(38,17,1:6)=(/ 0.820803,0.780870,0.745867,0.729900,0.715817,0.757550 /)
XEXT_COEFF_550_LKT(38,17)=328.080000 !rg=0.201336 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,18,1:6)=(/ 256.100000,263.660000,274.010000,289.050000,292.240000,281.260000 /)
XPIZA_LKT(38,18,1:6)=(/ 0.843823,0.938672,0.997686,0.998306,0.995767,0.570373 /)
XCGA_LKT(38,18,1:6)=(/ 0.828503,0.791863,0.757190,0.733177,0.724080,0.775410 /)
XEXT_COEFF_550_LKT(38,18)=267.990000 !rg=0.201336 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,19,1:6)=(/ 210.560000,216.150000,225.870000,234.780000,243.630000,236.830000 /)
XPIZA_LKT(38,19,1:6)=(/ 0.833233,0.932079,0.997024,0.998019,0.994718,0.561495 /)
XCGA_LKT(38,19,1:6)=(/ 0.837647,0.797637,0.763423,0.745637,0.726693,0.792903 /)
XEXT_COEFF_550_LKT(38,19)=219.880000 !rg=0.201336 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,20,1:6)=(/ 173.330000,177.630000,183.500000,192.780000,200.420000,197.860000 /)
XPIZA_LKT(38,20,1:6)=(/ 0.823175,0.926077,0.996833,0.997626,0.993843,0.553635 /)
XCGA_LKT(38,20,1:6)=(/ 0.845870,0.807353,0.770613,0.750987,0.729093,0.809583 /)
XEXT_COEFF_550_LKT(38,20)=180.160000 !rg=0.201336 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,1,1:6)=(/ 5572.700000,4028.400000,1812.800000,430.620000,56.123000,158.160000 /)
XPIZA_LKT(39,1,1:6)=(/ 0.985760,0.994814,0.999682,0.999294,0.988615,0.075609 /)
XCGA_LKT(39,1,1:6)=(/ 0.781477,0.753207,0.642120,0.291950,0.089240,0.033510 /)
XEXT_COEFF_550_LKT(39,1)=2935.900000 !rg=0.218119 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,2,1:6)=(/ 5150.600000,3970.700000,1889.500000,486.910000,69.431000,163.130000 /)
XPIZA_LKT(39,2,1:6)=(/ 0.984747,0.994719,0.999701,0.999353,0.990609,0.092547 /)
XCGA_LKT(39,2,1:6)=(/ 0.774317,0.752697,0.647300,0.366310,0.113207,0.042590 /)
XEXT_COEFF_550_LKT(39,2)=3075.800000 !rg=0.218119 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,3,1:6)=(/ 4462.000000,3858.300000,2045.200000,607.000000,97.957000,174.400000 /)
XPIZA_LKT(39,3,1:6)=(/ 0.982312,0.994498,0.999719,0.999445,0.993094,0.129871 /)
XCGA_LKT(39,3,1:6)=(/ 0.758720,0.756517,0.676467,0.473080,0.170610,0.065030 /)
XEXT_COEFF_550_LKT(39,3)=3142.100000 !rg=0.218119 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,4,1:6)=(/ 3710.200000,3570.500000,2187.000000,772.530000,144.910000,195.140000 /)
XPIZA_LKT(39,4,1:6)=(/ 0.978665,0.993960,0.999730,0.999534,0.995069,0.192496 /)
XCGA_LKT(39,4,1:6)=(/ 0.740807,0.754877,0.704297,0.558650,0.272563,0.110930 /)
XEXT_COEFF_550_LKT(39,4)=3087.700000 !rg=0.218119 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,5,1:6)=(/ 3055.700000,3142.000000,2238.300000,950.510000,213.720000,228.440000 /)
XPIZA_LKT(39,5,1:6)=(/ 0.973933,0.993034,0.999728,0.999597,0.996405,0.274702 /)
XCGA_LKT(39,5,1:6)=(/ 0.726463,0.747370,0.722803,0.622077,0.391803,0.185177 /)
XEXT_COEFF_550_LKT(39,5)=2878.900000 !rg=0.218119 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,6,1:6)=(/ 2511.100000,2675.800000,2166.400000,1107.700000,302.690000,275.240000 /)
XPIZA_LKT(39,6,1:6)=(/ 0.968512,0.991665,0.999711,0.999635,0.997248,0.360550 /)
XCGA_LKT(39,6,1:6)=(/ 0.721073,0.736807,0.731043,0.667280,0.493877,0.279270 /)
XEXT_COEFF_550_LKT(39,6)=2558.600000 !rg=0.218119 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,7,1:6)=(/ 2062.900000,2227.300000,1988.200000,1215.200000,405.680000,334.030000 /)
XPIZA_LKT(39,7,1:6)=(/ 0.962314,0.989898,0.999677,0.999651,0.997778,0.438780 /)
XCGA_LKT(39,7,1:6)=(/ 0.720967,0.727873,0.731320,0.697393,0.571883,0.373877 /)
XEXT_COEFF_550_LKT(39,7)=2197.900000 !rg=0.218119 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,8,1:6)=(/ 1698.700000,1838.700000,1754.000000,1246.000000,510.490000,398.100000 /)
XPIZA_LKT(39,8,1:6)=(/ 0.955106,0.987604,0.999625,0.999644,0.998103,0.502021 /)
XCGA_LKT(39,8,1:6)=(/ 0.724977,0.722413,0.727797,0.713820,0.629923,0.457823 /)
XEXT_COEFF_550_LKT(39,8)=1845.600000 !rg=0.218119 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,9,1:6)=(/ 1394.000000,1499.700000,1503.700000,1198.000000,602.650000,459.650000 /)
XPIZA_LKT(39,9,1:6)=(/ 0.947807,0.985404,0.999555,0.999616,0.998288,0.550422 /)
XCGA_LKT(39,9,1:6)=(/ 0.734647,0.723450,0.723773,0.719283,0.671407,0.529460 /)
XEXT_COEFF_550_LKT(39,9)=1525.400000 !rg=0.218119 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,10,1:6)=(/ 1140.000000,1228.500000,1265.100000,1094.000000,663.180000,507.430000 /)
XPIZA_LKT(39,10,1:6)=(/ 0.939816,0.981872,0.999434,0.999569,0.998358,0.583763 /)
XCGA_LKT(39,10,1:6)=(/ 0.746450,0.721797,0.719400,0.719380,0.698107,0.585907 /)
XEXT_COEFF_550_LKT(39,10)=1256.700000 !rg=0.218119 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,11,1:6)=(/ 936.980000,1004.200000,1045.500000,969.080000,679.660000,533.300000 /)
XPIZA_LKT(39,11,1:6)=(/ 0.928381,0.977832,0.999332,0.999500,0.998321,0.604098 /)
XCGA_LKT(39,11,1:6)=(/ 0.755103,0.728557,0.717483,0.718283,0.711787,0.629610 /)
XEXT_COEFF_550_LKT(39,11)=1025.400000 !rg=0.218119 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,12,1:6)=(/ 764.690000,814.660000,860.690000,834.950000,650.460000,532.210000 /)
XPIZA_LKT(39,12,1:6)=(/ 0.918110,0.974242,0.999071,0.999412,0.998174,0.612570 /)
XCGA_LKT(39,12,1:6)=(/ 0.769550,0.733337,0.718697,0.718813,0.714583,0.662093 /)
XEXT_COEFF_550_LKT(39,12)=836.010000 !rg=0.218119 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,13,1:6)=(/ 626.960000,662.240000,698.850000,702.810000,591.520000,506.180000 /)
XPIZA_LKT(39,13,1:6)=(/ 0.905118,0.969071,0.999045,0.999272,0.997936,0.611482 /)
XCGA_LKT(39,13,1:6)=(/ 0.781293,0.746123,0.723033,0.716850,0.712967,0.686533 /)
XEXT_COEFF_550_LKT(39,13)=677.630000 !rg=0.218119 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,14,1:6)=(/ 514.880000,538.250000,566.940000,585.310000,526.980000,463.550000 /)
XPIZA_LKT(39,14,1:6)=(/ 0.889888,0.964474,0.998850,0.999126,0.997626,0.604561 /)
XCGA_LKT(39,14,1:6)=(/ 0.791403,0.756587,0.731757,0.718363,0.714133,0.707040 /)
XEXT_COEFF_550_LKT(39,14)=550.890000 !rg=0.218119 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,15,1:6)=(/ 423.900000,441.170000,465.700000,483.760000,459.440000,413.460000 /)
XPIZA_LKT(39,15,1:6)=(/ 0.875230,0.956911,0.998570,0.998882,0.997193,0.595693 /)
XCGA_LKT(39,15,1:6)=(/ 0.804350,0.763620,0.736690,0.721827,0.715137,0.727127 /)
XEXT_COEFF_550_LKT(39,15)=452.520000 !rg=0.218119 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,16,1:6)=(/ 346.430000,361.460000,378.110000,394.260000,387.930000,360.840000 /)
XPIZA_LKT(39,16,1:6)=(/ 0.861443,0.948962,0.998342,0.998765,0.996638,0.586155 /)
XCGA_LKT(39,16,1:6)=(/ 0.814123,0.774880,0.740133,0.728413,0.715473,0.746987 /)
XEXT_COEFF_550_LKT(39,16)=366.470000 !rg=0.218119 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,17,1:6)=(/ 285.030000,294.580000,305.960000,323.970000,326.400000,309.570000 /)
XPIZA_LKT(39,17,1:6)=(/ 0.850394,0.941725,0.998066,0.998426,0.996066,0.576169 /)
XCGA_LKT(39,17,1:6)=(/ 0.824233,0.785440,0.754990,0.729970,0.722173,0.765373 /)
XEXT_COEFF_550_LKT(39,17)=298.400000 !rg=0.218119 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,18,1:6)=(/ 235.280000,242.070000,252.230000,265.580000,271.230000,262.290000 /)
XPIZA_LKT(39,18,1:6)=(/ 0.838660,0.935048,0.997517,0.998051,0.995325,0.566738 /)
XCGA_LKT(39,18,1:6)=(/ 0.834073,0.793673,0.760840,0.737817,0.723510,0.783100 /)
XEXT_COEFF_550_LKT(39,18)=247.000000 !rg=0.218119 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,19,1:6)=(/ 193.150000,199.030000,206.220000,215.420000,224.350000,219.850000 /)
XPIZA_LKT(39,19,1:6)=(/ 0.828098,0.928578,0.997113,0.997706,0.994095,0.558081 /)
XCGA_LKT(39,19,1:6)=(/ 0.840960,0.802697,0.765557,0.748277,0.725437,0.800343 /)
XEXT_COEFF_550_LKT(39,19)=200.570000 !rg=0.218119 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,20,1:6)=(/ 159.360000,163.290000,167.380000,176.700000,183.950000,183.040000 /)
XPIZA_LKT(39,20,1:6)=(/ 0.819193,0.922724,0.996789,0.997217,0.993261,0.550628 /)
XCGA_LKT(39,20,1:6)=(/ 0.848037,0.811317,0.778590,0.751777,0.736143,0.816753 /)
XEXT_COEFF_550_LKT(39,20)=164.380000 !rg=0.218119 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,1,1:6)=(/ 5135.100000,4043.400000,2017.000000,504.970000,70.636000,163.660000 /)
XPIZA_LKT(40,1,1:6)=(/ 0.984841,0.994863,0.999716,0.999375,0.990715,0.093051 /)
XCGA_LKT(40,1,1:6)=(/ 0.777103,0.756280,0.643860,0.350510,0.104783,0.039300 /)
XEXT_COEFF_550_LKT(40,1)=3195.300000 !rg=0.236301 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,2,1:6)=(/ 4691.000000,4021.900000,2064.700000,569.200000,86.794000,169.740000 /)
XPIZA_LKT(40,2,1:6)=(/ 0.983129,0.994715,0.999722,0.999420,0.992281,0.113188 /)
XCGA_LKT(40,2,1:6)=(/ 0.766850,0.762780,0.666427,0.429150,0.133023,0.049963 /)
XEXT_COEFF_550_LKT(40,2)=3227.500000 !rg=0.236301 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,3,1:6)=(/ 3992.500000,3817.100000,2194.600000,703.310000,120.210000,183.440000 /)
XPIZA_LKT(40,3,1:6)=(/ 0.980070,0.994365,0.999732,0.999499,0.994206,0.156485 /)
XCGA_LKT(40,3,1:6)=(/ 0.747750,0.762810,0.696043,0.520757,0.199950,0.076333 /)
XEXT_COEFF_550_LKT(40,3)=3222.500000 !rg=0.236301 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,4,1:6)=(/ 3303.900000,3435.900000,2285.800000,871.670000,173.790000,208.130000 /)
XPIZA_LKT(40,4,1:6)=(/ 0.975815,0.993658,0.999736,0.999570,0.995749,0.225491 /)
XCGA_LKT(40,4,1:6)=(/ 0.729597,0.756467,0.718627,0.590957,0.311363,0.129900 /)
XEXT_COEFF_550_LKT(40,4)=3080.400000 !rg=0.236301 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,5,1:6)=(/ 2719.800000,2955.500000,2270.700000,1040.800000,250.560000,246.620000 /)
XPIZA_LKT(40,5,1:6)=(/ 0.970863,0.992512,0.999727,0.999620,0.996818,0.310340 /)
XCGA_LKT(40,5,1:6)=(/ 0.719933,0.745067,0.731540,0.644660,0.429960,0.214047 /)
XEXT_COEFF_550_LKT(40,5)=2792.500000 !rg=0.236301 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,6,1:6)=(/ 2255.500000,2474.100000,2136.400000,1176.800000,345.950000,298.650000 /)
XPIZA_LKT(40,6,1:6)=(/ 0.965112,0.990981,0.999702,0.999647,0.997504,0.394068 /)
XCGA_LKT(40,6,1:6)=(/ 0.718417,0.734243,0.734890,0.682727,0.525103,0.312570 /)
XEXT_COEFF_550_LKT(40,6)=2430.100000 !rg=0.236301 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,7,1:6)=(/ 1866.200000,2040.000000,1916.200000,1251.900000,451.920000,360.740000 /)
XPIZA_LKT(40,7,1:6)=(/ 0.957952,0.988986,0.999659,0.999654,0.997937,0.466794 /)
XCGA_LKT(40,7,1:6)=(/ 0.721230,0.725533,0.731827,0.706927,0.596357,0.405500 /)
XEXT_COEFF_550_LKT(40,7)=2055.900000 !rg=0.236301 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,8,1:6)=(/ 1533.300000,1676.200000,1662.900000,1246.500000,553.300000,425.250000 /)
XPIZA_LKT(40,8,1:6)=(/ 0.951936,0.986423,0.999600,0.999639,0.998196,0.523866 /)
XCGA_LKT(40,8,1:6)=(/ 0.730337,0.720513,0.727307,0.718617,0.647927,0.486693 /)
XEXT_COEFF_550_LKT(40,8)=1707.300000 !rg=0.236301 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,9,1:6)=(/ 1263.800000,1368.500000,1403.100000,1168.000000,634.790000,482.440000 /)
XPIZA_LKT(40,9,1:6)=(/ 0.942917,0.983859,0.999522,0.999602,0.998332,0.566025 /)
XCGA_LKT(40,9,1:6)=(/ 0.736587,0.721400,0.723077,0.721243,0.683717,0.552527 /)
XEXT_COEFF_550_LKT(40,9)=1401.000000 !rg=0.236301 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,10,1:6)=(/ 1040.000000,1116.300000,1166.000000,1052.000000,678.180000,522.110000 /)
XPIZA_LKT(40,10,1:6)=(/ 0.934628,0.980843,0.999414,0.999543,0.998359,0.593924 /)
XCGA_LKT(40,10,1:6)=(/ 0.750137,0.727470,0.719403,0.719357,0.705197,0.604200 /)
XEXT_COEFF_550_LKT(40,10)=1146.000000 !rg=0.236301 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,11,1:6)=(/ 850.310000,910.470000,963.090000,919.340000,674.860000,537.280000 /)
XPIZA_LKT(40,11,1:6)=(/ 0.924472,0.976556,0.999265,0.999467,0.998279,0.609261 /)
XCGA_LKT(40,11,1:6)=(/ 0.762370,0.730270,0.720280,0.719643,0.714300,0.643433 /)
XEXT_COEFF_550_LKT(40,11)=933.160000 !rg=0.236301 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,12,1:6)=(/ 699.340000,740.450000,784.680000,782.150000,630.380000,525.350000 /)
XPIZA_LKT(40,12,1:6)=(/ 0.912302,0.972209,0.999124,0.999361,0.998096,0.613550 /)
XCGA_LKT(40,12,1:6)=(/ 0.774743,0.741620,0.720810,0.717603,0.714670,0.672517 /)
XEXT_COEFF_550_LKT(40,12)=759.040000 !rg=0.236301 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,13,1:6)=(/ 573.320000,602.460000,636.360000,652.010000,567.820000,491.270000 /)
XPIZA_LKT(40,13,1:6)=(/ 0.898234,0.967288,0.998996,0.999222,0.997823,0.609667 /)
XCGA_LKT(40,13,1:6)=(/ 0.785250,0.750673,0.728560,0.717403,0.713717,0.695213 /)
XEXT_COEFF_550_LKT(40,13)=617.380000 !rg=0.236301 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,14,1:6)=(/ 472.960000,492.040000,522.670000,542.970000,501.410000,444.470000 /)
XPIZA_LKT(40,14,1:6)=(/ 0.882471,0.961675,0.998707,0.999012,0.997466,0.601625 /)
XCGA_LKT(40,14,1:6)=(/ 0.798520,0.759053,0.731037,0.719300,0.714690,0.715557 /)
XEXT_COEFF_550_LKT(40,14)=506.380000 !rg=0.236301 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,15,1:6)=(/ 386.880000,404.860000,425.540000,444.120000,429.920000,392.690000 /)
XPIZA_LKT(40,15,1:6)=(/ 0.869113,0.953682,0.998136,0.998727,0.997022,0.592406 /)
XCGA_LKT(40,15,1:6)=(/ 0.807913,0.769787,0.735227,0.724123,0.716300,0.735700 /)
XEXT_COEFF_550_LKT(40,15)=411.790000 !rg=0.236301 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,16,1:6)=(/ 316.700000,328.640000,345.380000,361.840000,363.370000,339.700000 /)
XPIZA_LKT(40,16,1:6)=(/ 0.856480,0.946389,0.998128,0.998693,0.996221,0.582400 /)
XCGA_LKT(40,16,1:6)=(/ 0.820367,0.778620,0.747473,0.729087,0.717747,0.755070 /)
XEXT_COEFF_550_LKT(40,16)=333.020000 !rg=0.236301 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,17,1:6)=(/ 261.950000,269.660000,281.710000,297.620000,302.060000,289.490000 /)
XPIZA_LKT(40,17,1:6)=(/ 0.844277,0.939372,0.997828,0.998345,0.995845,0.572348 /)
XCGA_LKT(40,17,1:6)=(/ 0.829597,0.789410,0.756893,0.735100,0.722630,0.773190 /)
XEXT_COEFF_550_LKT(40,17)=274.650000 !rg=0.236301 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,18,1:6)=(/ 215.290000,222.440000,230.660000,242.400000,250.600000,244.030000 /)
XPIZA_LKT(40,18,1:6)=(/ 0.833762,0.932056,0.997335,0.997836,0.994644,0.563110 /)
XCGA_LKT(40,18,1:6)=(/ 0.836910,0.798920,0.761373,0.741503,0.723717,0.790737 /)
XEXT_COEFF_550_LKT(40,18)=225.520000 !rg=0.236301 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,19,1:6)=(/ 177.220000,181.620000,188.790000,197.590000,207.070000,203.710000 /)
XPIZA_LKT(40,19,1:6)=(/ 0.824098,0.926506,0.996666,0.997592,0.993566,0.554774 /)
XCGA_LKT(40,19,1:6)=(/ 0.845543,0.805793,0.771220,0.749597,0.729760,0.807727 /)
XEXT_COEFF_550_LKT(40,19)=183.180000 !rg=0.236301 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,20,1:6)=(/ 146.480000,149.630000,154.340000,162.270000,168.640000,169.090000 /)
XPIZA_LKT(40,20,1:6)=(/ 0.814279,0.921024,0.996599,0.997207,0.992969,0.547754 /)
XCGA_LKT(40,20,1:6)=(/ 0.851847,0.814060,0.781060,0.757583,0.739023,0.823823 /)
XEXT_COEFF_550_LKT(40,20)=151.320000 !rg=0.236301 sigma=2.95 wvl=0.55
 
END SUBROUTINE SALT_OPT_LKT_SET4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET5()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE
 
XEXT_COEFF_WVL_LKT(41,1,1:6)=(/ 4618.900000,4120.800000,2134.300000,584.510000,88.602000,170.390000 /)
XPIZA_LKT(41,1,1:6)=(/ 0.982625,0.994759,0.999736,0.999435,0.992386,0.113793 /)
XCGA_LKT(41,1,1:6)=(/ 0.768413,0.766463,0.655080,0.420410,0.123183,0.046097 /)
XEXT_COEFF_550_LKT(41,1)=3384.000000 !rg=0.255998 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,2,1:6)=(/ 4146.800000,4002.100000,2221.100000,663.370000,107.850000,177.820000 /)
XPIZA_LKT(41,2,1:6)=(/ 0.980808,0.994621,0.999736,0.999475,0.993606,0.137398 /)
XCGA_LKT(41,2,1:6)=(/ 0.753643,0.770050,0.688383,0.493040,0.156543,0.058637 /)
XEXT_COEFF_550_LKT(41,2)=3334.900000 !rg=0.255998 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,3,1:6)=(/ 3489.200000,3705.800000,2322.400000,807.790000,146.250000,194.380000 /)
XPIZA_LKT(41,3,1:6)=(/ 0.977162,0.994128,0.999742,0.999543,0.995088,0.186660 /)
XCGA_LKT(41,3,1:6)=(/ 0.732707,0.765907,0.713843,0.562720,0.234180,0.089640 /)
XEXT_COEFF_550_LKT(41,3)=3253.200000 !rg=0.255998 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,4,1:6)=(/ 2896.200000,3250.700000,2355.700000,972.580000,206.800000,223.490000 /)
XPIZA_LKT(41,4,1:6)=(/ 0.972495,0.993218,0.999739,0.999601,0.996298,0.260746 /)
XCGA_LKT(41,4,1:6)=(/ 0.717447,0.755070,0.730720,0.619473,0.352757,0.151843 /)
XEXT_COEFF_550_LKT(41,4)=3023.600000 !rg=0.255998 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,5,1:6)=(/ 2423.300000,2739.600000,2272.300000,1126.600000,291.290000,267.380000 /)
XPIZA_LKT(41,5,1:6)=(/ 0.966961,0.991878,0.999723,0.999639,0.997157,0.346084 /)
XCGA_LKT(41,5,1:6)=(/ 0.712357,0.741750,0.738143,0.664687,0.467160,0.244983 /)
XEXT_COEFF_550_LKT(41,5)=2671.800000 !rg=0.255998 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,6,1:6)=(/ 2014.700000,2272.400000,2079.400000,1235.700000,391.790000,324.070000 /)
XPIZA_LKT(41,6,1:6)=(/ 0.961419,0.990055,0.999690,0.999655,0.997714,0.425876 /)
XCGA_LKT(41,6,1:6)=(/ 0.714997,0.729893,0.737467,0.696137,0.554280,0.345497 /)
XEXT_COEFF_550_LKT(41,6)=2276.900000 !rg=0.255998 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,7,1:6)=(/ 1674.200000,1866.600000,1824.900000,1274.900000,497.970000,388.550000 /)
XPIZA_LKT(41,7,1:6)=(/ 0.954406,0.987684,0.999640,0.999653,0.998067,0.492441 /)
XCGA_LKT(41,7,1:6)=(/ 0.721667,0.721970,0.732450,0.714707,0.618470,0.436727 /)
XEXT_COEFF_550_LKT(41,7)=1898.100000 !rg=0.255998 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,8,1:6)=(/ 1387.100000,1525.500000,1558.300000,1233.900000,593.070000,451.370000 /)
XPIZA_LKT(41,8,1:6)=(/ 0.947556,0.985283,0.999570,0.999628,0.998271,0.543303 /)
XCGA_LKT(41,8,1:6)=(/ 0.732177,0.721320,0.726513,0.721827,0.663790,0.513197 /)
XEXT_COEFF_550_LKT(41,8)=1563.900000 !rg=0.255998 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,9,1:6)=(/ 1152.400000,1243.200000,1304.000000,1132.400000,660.740000,502.880000 /)
XPIZA_LKT(41,9,1:6)=(/ 0.937688,0.982103,0.999476,0.999583,0.998360,0.579461 /)
XCGA_LKT(41,9,1:6)=(/ 0.742967,0.720900,0.720540,0.722320,0.693963,0.574070 /)
XEXT_COEFF_550_LKT(41,9)=1288.500000 !rg=0.255998 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,10,1:6)=(/ 946.360000,1015.700000,1074.800000,1003.300000,685.060000,532.920000 /)
XPIZA_LKT(41,10,1:6)=(/ 0.928851,0.978964,0.999369,0.999519,0.998343,0.601937 /)
XCGA_LKT(41,10,1:6)=(/ 0.753083,0.727880,0.720307,0.721093,0.710307,0.620740 /)
XEXT_COEFF_550_LKT(41,10)=1045.600000 !rg=0.255998 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,11,1:6)=(/ 780.150000,825.460000,879.830000,867.490000,663.170000,536.860000 /)
XPIZA_LKT(41,11,1:6)=(/ 0.918157,0.975304,0.999235,0.999423,0.998222,0.612608 /)
XCGA_LKT(41,11,1:6)=(/ 0.767583,0.736887,0.720430,0.718950,0.715507,0.655887 /)
XEXT_COEFF_550_LKT(41,11)=852.390000 !rg=0.255998 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,12,1:6)=(/ 639.410000,674.110000,717.390000,726.450000,608.760000,514.920000 /)
XPIZA_LKT(41,12,1:6)=(/ 0.905482,0.970421,0.998957,0.999309,0.997998,0.613267 /)
XCGA_LKT(41,12,1:6)=(/ 0.778317,0.745267,0.724283,0.717580,0.714780,0.682220 /)
XEXT_COEFF_550_LKT(41,12)=692.640000 !rg=0.255998 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,13,1:6)=(/ 526.750000,549.720000,585.510000,607.530000,542.800000,474.270000 /)
XPIZA_LKT(41,13,1:6)=(/ 0.890384,0.964992,0.998784,0.999121,0.997697,0.607300 /)
XCGA_LKT(41,13,1:6)=(/ 0.792420,0.753473,0.728057,0.718530,0.714187,0.703863 /)
XEXT_COEFF_550_LKT(41,13)=566.750000 !rg=0.255998 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,14,1:6)=(/ 431.820000,452.450000,477.060000,497.310000,473.920000,424.430000 /)
XPIZA_LKT(41,14,1:6)=(/ 0.876314,0.957717,0.998505,0.998901,0.997301,0.598542 /)
XCGA_LKT(41,14,1:6)=(/ 0.801727,0.764183,0.731363,0.720793,0.716407,0.724233 /)
XEXT_COEFF_550_LKT(41,14)=461.460000 !rg=0.255998 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,15,1:6)=(/ 353.530000,368.670000,390.490000,407.570000,402.960000,371.500000 /)
XPIZA_LKT(41,15,1:6)=(/ 0.862980,0.950723,0.998144,0.998537,0.996646,0.588888 /)
XCGA_LKT(41,15,1:6)=(/ 0.815883,0.772000,0.740167,0.725683,0.716143,0.744150 /)
XEXT_COEFF_550_LKT(41,15)=375.760000 !rg=0.255998 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,16,1:6)=(/ 291.550000,299.540000,315.050000,334.740000,337.070000,318.740000 /)
XPIZA_LKT(41,16,1:6)=(/ 0.850662,0.943381,0.998003,0.998496,0.995964,0.578478 /)
XCGA_LKT(41,16,1:6)=(/ 0.825320,0.786463,0.752590,0.731140,0.719663,0.763040 /)
XEXT_COEFF_550_LKT(41,16)=306.600000 !rg=0.255998 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,17,1:6)=(/ 240.430000,247.740000,259.010000,271.570000,278.800000,270.050000 /)
XPIZA_LKT(41,17,1:6)=(/ 0.838862,0.936231,0.997240,0.998263,0.995533,0.568543 /)
XCGA_LKT(41,17,1:6)=(/ 0.832567,0.794703,0.755030,0.737010,0.722777,0.781007 /)
XEXT_COEFF_550_LKT(41,17)=251.150000 !rg=0.255998 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,18,1:6)=(/ 197.520000,203.310000,212.240000,221.490000,231.990000,226.560000 /)
XPIZA_LKT(41,18,1:6)=(/ 0.828365,0.929728,0.997105,0.997823,0.994416,0.559518 /)
XCGA_LKT(41,18,1:6)=(/ 0.842567,0.801443,0.765133,0.745027,0.727370,0.798310 /)
XEXT_COEFF_550_LKT(41,18)=206.450000 !rg=0.255998 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,19,1:6)=(/ 163.120000,166.180000,172.610000,182.090000,189.820000,188.460000 /)
XPIZA_LKT(41,19,1:6)=(/ 0.819413,0.923911,0.996699,0.997333,0.992760,0.551583 /)
XCGA_LKT(41,19,1:6)=(/ 0.849160,0.811927,0.777400,0.752417,0.732523,0.815043 /)
XEXT_COEFF_550_LKT(41,19)=168.980000 !rg=0.255998 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,20,1:6)=(/ 134.730000,137.540000,142.180000,147.560000,155.090000,156.010000 /)
XPIZA_LKT(41,20,1:6)=(/ 0.809725,0.918497,0.996110,0.997009,0.992379,0.545036 /)
XCGA_LKT(41,20,1:6)=(/ 0.854293,0.818230,0.779603,0.760463,0.740400,0.830780 /)
XEXT_COEFF_550_LKT(41,20)=138.720000 !rg=0.255998 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,1,1:6)=(/ 3985.900000,4116.500000,2249.800000,674.420000,110.580000,178.620000 /)
XPIZA_LKT(42,1,1:6)=(/ 0.980128,0.994754,0.999742,0.999481,0.993714,0.138138 /)
XCGA_LKT(42,1,1:6)=(/ 0.749597,0.779167,0.688210,0.498003,0.145083,0.054087 /)
XEXT_COEFF_550_LKT(42,1)=3406.600000 !rg=0.277337 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,2,1:6)=(/ 3556.500000,3904.300000,2363.300000,772.240000,132.960000,187.670000 /)
XPIZA_LKT(42,2,1:6)=(/ 0.977501,0.994422,0.999746,0.999523,0.994653,0.165350 /)
XCGA_LKT(42,2,1:6)=(/ 0.733330,0.774330,0.711740,0.550263,0.184593,0.068853 /)
XEXT_COEFF_550_LKT(42,2)=3388.200000 !rg=0.277337 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,3,1:6)=(/ 3000.900000,3525.500000,2423.300000,917.860000,176.320000,207.560000 /)
XPIZA_LKT(42,3,1:6)=(/ 0.973342,0.993764,0.999748,0.999581,0.995787,0.220142 /)
XCGA_LKT(42,3,1:6)=(/ 0.714437,0.765643,0.729543,0.598100,0.273687,0.105330 /)
XEXT_COEFF_550_LKT(42,3)=3228.100000 !rg=0.277337 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,4,1:6)=(/ 2531.200000,3018.400000,2392.700000,1072.400000,244.170000,241.490000 /)
XPIZA_LKT(42,4,1:6)=(/ 0.968555,0.992646,0.999739,0.999626,0.996743,0.297487 /)
XCGA_LKT(42,4,1:6)=(/ 0.706017,0.751083,0.740510,0.644490,0.395697,0.176960 /)
XEXT_COEFF_550_LKT(42,4)=2918.400000 !rg=0.277337 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,5,1:6)=(/ 2139.500000,2512.000000,2241.400000,1205.100000,335.600000,290.720000 /)
XPIZA_LKT(42,5,1:6)=(/ 0.963026,0.991026,0.999715,0.999653,0.997437,0.381113 /)
XCGA_LKT(42,5,1:6)=(/ 0.706750,0.735577,0.742870,0.682257,0.502617,0.277443 /)
XEXT_COEFF_550_LKT(42,5)=2516.100000 !rg=0.277337 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,6,1:6)=(/ 1803.200000,2062.100000,2001.400000,1281.900000,439.570000,351.290000 /)
XPIZA_LKT(42,6,1:6)=(/ 0.957449,0.989083,0.999674,0.999661,0.997890,0.455567 /)
XCGA_LKT(42,6,1:6)=(/ 0.714507,0.724660,0.738570,0.707470,0.581423,0.378370 /)
XEXT_COEFF_550_LKT(42,6)=2112.500000 !rg=0.277337 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,7,1:6)=(/ 1505.900000,1688.500000,1725.000000,1282.900000,543.310000,416.780000 /)
XPIZA_LKT(42,7,1:6)=(/ 0.950937,0.986594,0.999615,0.999649,0.998171,0.515766 /)
XCGA_LKT(42,7,1:6)=(/ 0.725597,0.719020,0.731597,0.720907,0.638373,0.466823 /)
XEXT_COEFF_550_LKT(42,7)=1743.200000 !rg=0.277337 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,8,1:6)=(/ 1255.400000,1381.000000,1450.800000,1208.700000,628.770000,475.980000 /)
XPIZA_LKT(42,8,1:6)=(/ 0.942943,0.983992,0.999538,0.999616,0.998324,0.560271 /)
XCGA_LKT(42,8,1:6)=(/ 0.735017,0.719370,0.725393,0.724623,0.677623,0.537707 /)
XEXT_COEFF_550_LKT(42,8)=1430.300000 !rg=0.277337 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,9,1:6)=(/ 1042.500000,1136.100000,1200.500000,1088.600000,680.070000,519.900000 /)
XPIZA_LKT(42,9,1:6)=(/ 0.933335,0.979992,0.999428,0.999561,0.998370,0.590814 /)
XCGA_LKT(42,9,1:6)=(/ 0.746370,0.722350,0.719367,0.722960,0.702323,0.593690 /)
XEXT_COEFF_550_LKT(42,9)=1170.400000 !rg=0.277337 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,10,1:6)=(/ 866.360000,922.810000,991.640000,953.220000,684.610000,539.480000 /)
XPIZA_LKT(42,10,1:6)=(/ 0.923043,0.977056,0.999294,0.999489,0.998311,0.608188 /)
XCGA_LKT(42,10,1:6)=(/ 0.760833,0.728830,0.718140,0.721710,0.713907,0.635747 /)
XEXT_COEFF_550_LKT(42,10)=958.920000 !rg=0.277337 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,11,1:6)=(/ 711.780000,755.840000,805.120000,809.350000,646.690000,532.250000 /)
XPIZA_LKT(42,11,1:6)=(/ 0.911660,0.972628,0.999170,0.999382,0.998148,0.614506 /)
XCGA_LKT(42,11,1:6)=(/ 0.770467,0.738360,0.720370,0.719033,0.716320,0.667253 /)
XEXT_COEFF_550_LKT(42,11)=777.470000 !rg=0.277337 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,12,1:6)=(/ 586.760000,614.720000,659.330000,677.800000,585.030000,501.420000 /)
XPIZA_LKT(42,12,1:6)=(/ 0.898386,0.967953,0.998939,0.999213,0.997886,0.612069 /)
XCGA_LKT(42,12,1:6)=(/ 0.786023,0.747547,0.723210,0.717970,0.714457,0.691467 /)
XEXT_COEFF_550_LKT(42,12)=636.480000 !rg=0.277337 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,13,1:6)=(/ 480.880000,505.430000,535.080000,557.720000,516.900000,455.680000 /)
XPIZA_LKT(42,13,1:6)=(/ 0.883726,0.961270,0.998601,0.999096,0.997555,0.604579 /)
XCGA_LKT(42,13,1:6)=(/ 0.795537,0.758593,0.727500,0.719773,0.716380,0.712597 /)
XEXT_COEFF_550_LKT(42,13)=516.970000 !rg=0.277337 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,14,1:6)=(/ 394.660000,412.390000,437.040000,456.780000,445.890000,403.660000 /)
XPIZA_LKT(42,14,1:6)=(/ 0.869892,0.954937,0.998472,0.998869,0.997061,0.595278 /)
XCGA_LKT(42,14,1:6)=(/ 0.810260,0.765660,0.737110,0.724293,0.716323,0.732963 /)
XEXT_COEFF_550_LKT(42,14)=422.680000 !rg=0.277337 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,15,1:6)=(/ 325.300000,335.960000,355.330000,375.760000,374.270000,350.100000 /)
XPIZA_LKT(42,15,1:6)=(/ 0.857213,0.947319,0.998058,0.998402,0.996511,0.585092 /)
XCGA_LKT(42,15,1:6)=(/ 0.820723,0.780990,0.742727,0.727653,0.716860,0.752413 /)
XEXT_COEFF_550_LKT(42,15)=342.690000 !rg=0.277337 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,16,1:6)=(/ 267.480000,276.160000,288.220000,305.360000,311.780000,298.200000 /)
XPIZA_LKT(42,16,1:6)=(/ 0.844699,0.939776,0.997910,0.998440,0.995911,0.574495 /)
XCGA_LKT(42,16,1:6)=(/ 0.828227,0.789840,0.754547,0.730843,0.722107,0.770980 /)
XEXT_COEFF_550_LKT(42,16)=282.180000 !rg=0.277337 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,17,1:6)=(/ 220.360000,227.250000,237.060000,247.110000,259.070000,251.300000 /)
XPIZA_LKT(42,17,1:6)=(/ 0.834017,0.932655,0.997333,0.998206,0.995021,0.564738 /)
XCGA_LKT(42,17,1:6)=(/ 0.837963,0.796300,0.762830,0.744630,0.724697,0.788787 /)
XEXT_COEFF_550_LKT(42,17)=231.370000 !rg=0.277337 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,18,1:6)=(/ 181.830000,186.150000,193.340000,203.620000,212.490000,209.950000 /)
XPIZA_LKT(42,18,1:6)=(/ 0.824078,0.926717,0.996883,0.997744,0.994163,0.556010 /)
XCGA_LKT(42,18,1:6)=(/ 0.845920,0.807827,0.767827,0.750140,0.728277,0.805840 /)
XEXT_COEFF_550_LKT(42,18)=189.000000 !rg=0.277337 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,19,1:6)=(/ 149.910000,153.360000,158.600000,165.980000,173.690000,174.080000 /)
XPIZA_LKT(42,19,1:6)=(/ 0.814479,0.921027,0.996569,0.997217,0.993012,0.548525 /)
XCGA_LKT(42,19,1:6)=(/ 0.851410,0.814853,0.778163,0.754113,0.736483,0.822260 /)
XEXT_COEFF_550_LKT(42,19)=155.670000 !rg=0.277337 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,20,1:6)=(/ 123.690000,126.570000,130.330000,134.470000,143.340000,143.770000 /)
XPIZA_LKT(42,20,1:6)=(/ 0.805664,0.915322,0.996319,0.996829,0.991268,0.542493 /)
XCGA_LKT(42,20,1:6)=(/ 0.857580,0.819770,0.787023,0.768283,0.742830,0.837607 /)
XEXT_COEFF_550_LKT(42,20)=127.910000 !rg=0.277337 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,1,1:6)=(/ 3339.600000,3937.200000,2442.600000,786.600000,137.010000,188.690000 /)
XPIZA_LKT(43,1,1:6)=(/ 0.975775,0.994461,0.999748,0.999521,0.994767,0.166283 /)
XCGA_LKT(43,1,1:6)=(/ 0.722097,0.781673,0.728420,0.571240,0.171317,0.063493 /)
XEXT_COEFF_550_LKT(43,1)=3461.400000 !rg=0.300455 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,2,1:6)=(/ 2960.200000,3724.600000,2485.700000,894.590000,162.270000,199.670000 /)
XPIZA_LKT(43,2,1:6)=(/ 0.972948,0.994091,0.999754,0.999566,0.995478,0.197041 /)
XCGA_LKT(43,2,1:6)=(/ 0.706190,0.774830,0.732370,0.594347,0.218190,0.080913 /)
XEXT_COEFF_550_LKT(43,2)=3382.300000 !rg=0.300455 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,3,1:6)=(/ 2552.200000,3281.700000,2491.900000,1029.700000,210.640000,223.320000 /)
XPIZA_LKT(43,3,1:6)=(/ 0.968427,0.993238,0.999751,0.999613,0.996344,0.256409 /)
XCGA_LKT(43,3,1:6)=(/ 0.691417,0.761737,0.742683,0.627467,0.318397,0.123847 /)
XEXT_COEFF_550_LKT(43,3)=3144.300000 !rg=0.300455 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,4,1:6)=(/ 2212.000000,2757.000000,2393.400000,1167.800000,285.970000,262.310000 /)
XPIZA_LKT(43,4,1:6)=(/ 0.963743,0.991882,0.999735,0.999647,0.997106,0.334811 /)
XCGA_LKT(43,4,1:6)=(/ 0.693043,0.744387,0.747883,0.666413,0.438827,0.205333 /)
XEXT_COEFF_550_LKT(43,4)=2767.700000 !rg=0.300455 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,5,1:6)=(/ 1899.300000,2268.100000,2180.600000,1273.400000,383.210000,316.540000 /)
XPIZA_LKT(43,5,1:6)=(/ 0.959072,0.990059,0.999704,0.999664,0.997667,0.414719 /)
XCGA_LKT(43,5,1:6)=(/ 0.705150,0.728887,0.745407,0.697553,0.536040,0.311170 /)
XEXT_COEFF_550_LKT(43,5)=2340.200000 !rg=0.300455 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,6,1:6)=(/ 1627.100000,1862.700000,1903.000000,1313.700000,488.040000,380.020000 /)
XPIZA_LKT(43,6,1:6)=(/ 0.952433,0.987911,0.999653,0.999662,0.998034,0.482982 /)
XCGA_LKT(43,6,1:6)=(/ 0.715027,0.720243,0.737677,0.716930,0.606137,0.411170 /)
XEXT_COEFF_550_LKT(43,6)=1942.100000 !rg=0.300455 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,7,1:6)=(/ 1366.600000,1524.900000,1611.700000,1276.700000,586.310000,444.510000 /)
XPIZA_LKT(43,7,1:6)=(/ 0.945662,0.985464,0.999585,0.999641,0.998257,0.536665 /)
XCGA_LKT(43,7,1:6)=(/ 0.727987,0.718193,0.729457,0.725363,0.656090,0.494990 /)
XEXT_COEFF_550_LKT(43,7)=1591.100000 !rg=0.300455 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,8,1:6)=(/ 1146.200000,1249.700000,1342.600000,1175.300000,658.960000,498.570000 /)
XPIZA_LKT(43,8,1:6)=(/ 0.937086,0.982465,0.999492,0.999598,0.998362,0.575027 /)
XCGA_LKT(43,8,1:6)=(/ 0.741140,0.718807,0.721860,0.726037,0.689430,0.560737 /)
XEXT_COEFF_550_LKT(43,8)=1307.500000 !rg=0.300455 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,9,1:6)=(/ 943.280000,1026.300000,1107.400000,1040.400000,691.480000,533.190000 /)
XPIZA_LKT(43,9,1:6)=(/ 0.929635,0.978622,0.999355,0.999535,0.998365,0.600010 /)
XCGA_LKT(43,9,1:6)=(/ 0.755720,0.722393,0.718890,0.723993,0.708843,0.611560 /)
XEXT_COEFF_550_LKT(43,9)=1061.900000 !rg=0.300455 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,10,1:6)=(/ 786.820000,845.450000,904.560000,896.570000,676.800000,541.560000 /)
XPIZA_LKT(43,10,1:6)=(/ 0.917825,0.974290,0.999235,0.999452,0.998264,0.612618 /)
XCGA_LKT(43,10,1:6)=(/ 0.763887,0.732380,0.717210,0.721717,0.716403,0.649323 /)
XEXT_COEFF_550_LKT(43,10)=870.660000 !rg=0.300455 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,11,1:6)=(/ 651.680000,688.280000,740.670000,753.110000,626.590000,523.740000 /)
XPIZA_LKT(43,11,1:6)=(/ 0.905111,0.969704,0.999046,0.999339,0.998055,0.615055 /)
XCGA_LKT(43,11,1:6)=(/ 0.779267,0.739593,0.720797,0.720310,0.715960,0.677760 /)
XEXT_COEFF_550_LKT(43,11)=712.990000 !rg=0.300455 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,12,1:6)=(/ 534.860000,564.900000,600.780000,625.190000,559.980000,485.460000 /)
XPIZA_LKT(43,12,1:6)=(/ 0.891695,0.964391,0.998768,0.999166,0.997771,0.610155 /)
XCGA_LKT(43,12,1:6)=(/ 0.789290,0.752610,0.723230,0.719133,0.716107,0.700543 /)
XEXT_COEFF_550_LKT(43,12)=578.530000 !rg=0.300455 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,13,1:6)=(/ 439.440000,460.330000,489.620000,511.770000,490.360000,435.880000 /)
XPIZA_LKT(43,13,1:6)=(/ 0.876841,0.958695,0.998649,0.999018,0.997374,0.601607 /)
XCGA_LKT(43,13,1:6)=(/ 0.804557,0.759663,0.733173,0.722500,0.716650,0.721460 /)
XEXT_COEFF_550_LKT(43,13)=473.020000 !rg=0.300455 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,14,1:6)=(/ 362.060000,376.410000,396.770000,421.960000,415.560000,382.330000 /)
XPIZA_LKT(43,14,1:6)=(/ 0.863874,0.950954,0.998395,0.998793,0.996823,0.591746 /)
XCGA_LKT(43,14,1:6)=(/ 0.815843,0.774450,0.739853,0.724220,0.715933,0.741597 /)
XEXT_COEFF_550_LKT(43,14)=384.720000 !rg=0.300455 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,15,1:6)=(/ 298.910000,308.530000,322.420000,343.940000,347.620000,328.750000 /)
XPIZA_LKT(43,15,1:6)=(/ 0.850674,0.943839,0.998111,0.998535,0.996394,0.581073 /)
XCGA_LKT(43,15,1:6)=(/ 0.823793,0.785813,0.750840,0.727213,0.721560,0.760537 /)
XEXT_COEFF_550_LKT(43,15)=314.050000 !rg=0.300455 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,16,1:6)=(/ 245.920000,253.210000,265.300000,278.970000,289.180000,278.240000 /)
XPIZA_LKT(43,16,1:6)=(/ 0.839071,0.936174,0.997534,0.998213,0.995537,0.570501 /)
XCGA_LKT(43,16,1:6)=(/ 0.833867,0.791923,0.757457,0.738593,0.722167,0.778930 /)
XEXT_COEFF_550_LKT(43,16)=258.220000 !rg=0.300455 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,17,1:6)=(/ 202.720000,208.100000,215.840000,228.050000,238.640000,233.320000 /)
XPIZA_LKT(43,17,1:6)=(/ 0.828573,0.929823,0.997279,0.997717,0.994656,0.560942 /)
XCGA_LKT(43,17,1:6)=(/ 0.842393,0.802367,0.764773,0.744970,0.723790,0.796513 /)
XEXT_COEFF_550_LKT(43,17)=210.960000 !rg=0.300455 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,18,1:6)=(/ 167.290000,171.100000,176.190000,186.450000,194.510000,194.210000 /)
XPIZA_LKT(43,18,1:6)=(/ 0.819260,0.924120,0.996864,0.997477,0.993927,0.552620 /)
XCGA_LKT(43,18,1:6)=(/ 0.848307,0.812013,0.775810,0.748967,0.735413,0.813303 /)
XEXT_COEFF_550_LKT(43,18)=173.700000 !rg=0.300455 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,19,1:6)=(/ 138.000000,140.650000,145.820000,151.360000,160.470000,160.580000 /)
XPIZA_LKT(43,19,1:6)=(/ 0.810074,0.918241,0.996289,0.997057,0.992484,0.545625 /)
XCGA_LKT(43,19,1:6)=(/ 0.855027,0.816927,0.780943,0.762083,0.739750,0.829373 /)
XEXT_COEFF_550_LKT(43,19)=142.760000 !rg=0.300455 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,20,1:6)=(/ 114.000000,115.980000,118.860000,124.130000,130.930000,132.370000 /)
XPIZA_LKT(43,20,1:6)=(/ 0.800851,0.913199,0.996209,0.996561,0.991101,0.540134 /)
XCGA_LKT(43,20,1:6)=(/ 0.860903,0.824160,0.788407,0.769733,0.743427,0.844287 /)
XEXT_COEFF_550_LKT(43,20)=117.180000 !rg=0.300455 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,1,1:6)=(/ 2660.700000,3773.100000,2611.000000,931.260000,168.100000,200.980000 /)
XPIZA_LKT(44,1,1:6)=(/ 0.969851,0.994077,0.999760,0.999565,0.995597,0.198258 /)
XCGA_LKT(44,1,1:6)=(/ 0.682650,0.778853,0.750110,0.621850,0.202997,0.074590 /)
XEXT_COEFF_550_LKT(44,1)=3488.500000 !rg=0.3255 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,2,1:6)=(/ 2424.400000,3464.800000,2575.700000,1022.500000,195.720000,214.200000 /)
XPIZA_LKT(44,2,1:6)=(/ 0.966588,0.993588,0.999758,0.999605,0.996125,0.232226 /)
XCGA_LKT(44,2,1:6)=(/ 0.669410,0.770830,0.747913,0.624300,0.258507,0.095190 /)
XEXT_COEFF_550_LKT(44,2)=3311.200000 !rg=0.3255 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,3,1:6)=(/ 2160.000000,2986.300000,2522.700000,1139.100000,249.520000,241.980000 /)
XPIZA_LKT(44,3,1:6)=(/ 0.962994,0.992502,0.999750,0.999640,0.996790,0.294667 /)
XCGA_LKT(44,3,1:6)=(/ 0.673873,0.753590,0.752933,0.652273,0.367440,0.145707 /)
XEXT_COEFF_550_LKT(44,3)=3001.200000 !rg=0.3255 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,4,1:6)=(/ 1933.700000,2477.000000,2356.000000,1255.500000,332.170000,286.100000 /)
XPIZA_LKT(44,4,1:6)=(/ 0.958937,0.990878,0.999727,0.999663,0.997404,0.371783 /)
XCGA_LKT(44,4,1:6)=(/ 0.687960,0.734063,0.752767,0.685650,0.480787,0.236957 /)
XEXT_COEFF_550_LKT(44,4)=2577.600000 !rg=0.3255 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,5,1:6)=(/ 1698.500000,2034.800000,2089.700000,1328.500000,433.370000,344.700000 /)
XPIZA_LKT(44,5,1:6)=(/ 0.953984,0.988910,0.999687,0.999671,0.997859,0.446388 /)
XCGA_LKT(44,5,1:6)=(/ 0.703003,0.721613,0.745890,0.710660,0.567013,0.345963 /)
XEXT_COEFF_550_LKT(44,5)=2145.500000 !rg=0.3255 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,6,1:6)=(/ 1459.600000,1679.700000,1787.600000,1329.400000,536.450000,409.630000 /)
XPIZA_LKT(44,6,1:6)=(/ 0.947471,0.986382,0.999628,0.999660,0.998152,0.508035 /)
XCGA_LKT(44,6,1:6)=(/ 0.716860,0.713827,0.736380,0.724540,0.628527,0.443193 /)
XEXT_COEFF_550_LKT(44,6)=1764.000000 !rg=0.3255 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,7,1:6)=(/ 1234.300000,1378.300000,1493.500000,1255.600000,625.890000,471.160000 /)
XPIZA_LKT(44,7,1:6)=(/ 0.940505,0.983670,0.999549,0.999630,0.998321,0.555109 /)
XCGA_LKT(44,7,1:6)=(/ 0.731017,0.714280,0.727187,0.728887,0.671730,0.521440 /)
XEXT_COEFF_550_LKT(44,7)=1443.700000 !rg=0.3255 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,8,1:6)=(/ 1038.100000,1137.700000,1232.200000,1131.200000,682.850000,518.090000 /)
XPIZA_LKT(44,8,1:6)=(/ 0.932131,0.980296,0.999444,0.999579,0.998382,0.587672 /)
XCGA_LKT(44,8,1:6)=(/ 0.743607,0.719110,0.720163,0.727093,0.699297,0.581967 /)
XEXT_COEFF_550_LKT(44,8)=1182.600000 !rg=0.3255 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,9,1:6)=(/ 863.310000,926.340000,1010.400000,988.340000,695.330000,542.340000 /)
XPIZA_LKT(44,9,1:6)=(/ 0.923585,0.977420,0.999290,0.999503,0.998343,0.607380 /)
XCGA_LKT(44,9,1:6)=(/ 0.761033,0.729663,0.717923,0.723440,0.713740,0.627923 /)
XEXT_COEFF_550_LKT(44,9)=964.070000 !rg=0.3255 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,10,1:6)=(/ 714.390000,765.710000,830.540000,838.270000,663.650000,539.310000 /)
XPIZA_LKT(44,10,1:6)=(/ 0.912961,0.972500,0.999099,0.999413,0.998196,0.615492 /)
XCGA_LKT(44,10,1:6)=(/ 0.774773,0.732657,0.717220,0.721533,0.717467,0.661783 /)
XEXT_COEFF_550_LKT(44,10)=791.460000 !rg=0.3255 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,11,1:6)=(/ 592.870000,629.140000,673.840000,697.560000,603.000000,511.820000 /)
XPIZA_LKT(44,11,1:6)=(/ 0.898993,0.967395,0.998986,0.999293,0.997962,0.614550 /)
XCGA_LKT(44,11,1:6)=(/ 0.784247,0.746480,0.719503,0.719447,0.717090,0.687723 /)
XEXT_COEFF_550_LKT(44,11)=645.710000 !rg=0.3255 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,12,1:6)=(/ 487.850000,513.600000,550.880000,575.520000,534.880000,467.530000 /)
XPIZA_LKT(44,12,1:6)=(/ 0.884980,0.962237,0.998408,0.999129,0.997635,0.607741 /)
XCGA_LKT(44,12,1:6)=(/ 0.799053,0.753560,0.727247,0.721043,0.716893,0.709603 /)
XEXT_COEFF_550_LKT(44,12)=527.930000 !rg=0.3255 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,13,1:6)=(/ 402.650000,420.360000,445.160000,472.960000,459.820000,415.130000 /)
XPIZA_LKT(44,13,1:6)=(/ 0.870762,0.954985,0.998556,0.998913,0.997201,0.598377 /)
XCGA_LKT(44,13,1:6)=(/ 0.810523,0.768423,0.734420,0.721390,0.717810,0.730373 /)
XEXT_COEFF_550_LKT(44,13)=430.850000 !rg=0.3255 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,14,1:6)=(/ 332.360000,343.980000,360.060000,385.560000,386.410000,360.640000 /)
XPIZA_LKT(44,14,1:6)=(/ 0.857620,0.947785,0.998332,0.998691,0.996722,0.587897 /)
XCGA_LKT(44,14,1:6)=(/ 0.819387,0.780213,0.748250,0.723870,0.719970,0.750057 /)
XEXT_COEFF_550_LKT(44,14)=350.810000 !rg=0.3255 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,15,1:6)=(/ 274.580000,283.000000,296.990000,315.690000,323.460000,307.730000 /)
XPIZA_LKT(44,15,1:6)=(/ 0.845039,0.940052,0.997915,0.998285,0.996041,0.576926 /)
XCGA_LKT(44,15,1:6)=(/ 0.830033,0.787910,0.753853,0.731077,0.721113,0.768617 /)
XEXT_COEFF_550_LKT(44,15)=289.350000 !rg=0.3255 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,16,1:6)=(/ 225.230000,232.670000,242.340000,254.530000,266.520000,258.960000 /)
XPIZA_LKT(44,16,1:6)=(/ 0.834109,0.933226,0.997427,0.998208,0.995020,0.566492 /)
XCGA_LKT(44,16,1:6)=(/ 0.837297,0.798160,0.758433,0.742573,0.721437,0.786860 /)
XEXT_COEFF_550_LKT(44,16)=235.140000 !rg=0.3255 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,17,1:6)=(/ 186.110000,190.920000,196.620000,208.960000,218.930000,216.190000 /)
XPIZA_LKT(44,17,1:6)=(/ 0.824594,0.926350,0.997109,0.997620,0.994199,0.557217 /)
XCGA_LKT(44,17,1:6)=(/ 0.845210,0.806150,0.773177,0.745113,0.730213,0.804190 /)
XEXT_COEFF_550_LKT(44,17)=192.460000 !rg=0.3255 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,18,1:6)=(/ 153.860000,157.230000,162.470000,171.060000,179.480000,179.370000 /)
XPIZA_LKT(44,18,1:6)=(/ 0.814585,0.921217,0.996607,0.997229,0.993287,0.549358 /)
XCGA_LKT(44,18,1:6)=(/ 0.852463,0.813903,0.778453,0.754657,0.735863,0.820683 /)
XEXT_COEFF_550_LKT(44,18)=159.550000 !rg=0.3255 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,19,1:6)=(/ 126.670000,129.590000,133.670000,138.800000,147.050000,147.960000 /)
XPIZA_LKT(44,19,1:6)=(/ 0.805799,0.915645,0.996261,0.996646,0.991396,0.542899 /)
XCGA_LKT(44,19,1:6)=(/ 0.857250,0.820947,0.783390,0.766840,0.738893,0.836360 /)
XEXT_COEFF_550_LKT(44,19)=130.300000 !rg=0.3255 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,20,1:6)=(/ 104.790000,106.780000,108.890000,114.030000,119.850000,121.760000 /)
XPIZA_LKT(44,20,1:6)=(/ 0.797111,0.910168,0.995972,0.996041,0.990053,0.537978 /)
XCGA_LKT(44,20,1:6)=(/ 0.862500,0.826847,0.794347,0.769437,0.750587,0.850793 /)
XEXT_COEFF_550_LKT(44,20)=107.270000 !rg=0.3255 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,1,1:6)=(/ 2069.500000,3436.100000,2642.700000,1093.000000,203.620000,215.920000 /)
XPIZA_LKT(45,1,1:6)=(/ 0.961197,0.993591,0.999766,0.999613,0.996246,0.233866 /)
XCGA_LKT(45,1,1:6)=(/ 0.627913,0.773137,0.754397,0.641423,0.241580,0.087727 /)
XEXT_COEFF_550_LKT(45,1)=3344.100000 !rg=0.352632 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,2,1:6)=(/ 1956.900000,3134.300000,2624.600000,1144.500000,233.160000,231.660000 /)
XPIZA_LKT(45,2,1:6)=(/ 0.959190,0.992854,0.999759,0.999638,0.996629,0.270372 /)
XCGA_LKT(45,2,1:6)=(/ 0.635543,0.761757,0.759007,0.645127,0.306610,0.112163 /)
XEXT_COEFF_550_LKT(45,2)=3171.000000 !rg=0.352632 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,3,1:6)=(/ 1871.100000,2653.800000,2511.000000,1241.900000,293.330000,263.790000 /)
XPIZA_LKT(45,3,1:6)=(/ 0.956807,0.991524,0.999745,0.999660,0.997151,0.333909 /)
XCGA_LKT(45,3,1:6)=(/ 0.658953,0.741710,0.760120,0.674070,0.418943,0.171470 /)
XEXT_COEFF_550_LKT(45,3)=2803.700000 !rg=0.352632 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,4,1:6)=(/ 1722.900000,2191.000000,2280.900000,1332.100000,382.470000,312.870000 /)
XPIZA_LKT(45,4,1:6)=(/ 0.953437,0.989679,0.999714,0.999675,0.997650,0.407551 /)
XCGA_LKT(45,4,1:6)=(/ 0.686170,0.723317,0.754877,0.702483,0.520373,0.271677 /)
XEXT_COEFF_550_LKT(45,4)=2361.700000 !rg=0.352632 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,5,1:6)=(/ 1520.000000,1812.100000,1974.600000,1368.100000,485.130000,374.840000 /)
XPIZA_LKT(45,5,1:6)=(/ 0.949189,0.987434,0.999665,0.999674,0.998017,0.475802 /)
XCGA_LKT(45,5,1:6)=(/ 0.706487,0.711873,0.744670,0.721670,0.595240,0.381373 /)
XEXT_COEFF_550_LKT(45,5)=1943.200000 !rg=0.352632 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,6,1:6)=(/ 1319.500000,1500.000000,1664.000000,1328.900000,583.230000,439.320000 /)
XPIZA_LKT(45,6,1:6)=(/ 0.942702,0.984788,0.999594,0.999654,0.998248,0.530616 /)
XCGA_LKT(45,6,1:6)=(/ 0.726047,0.710337,0.732500,0.730433,0.648563,0.473830 /)
XEXT_COEFF_550_LKT(45,6)=1600.000000 !rg=0.352632 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,7,1:6)=(/ 1125.300000,1239.100000,1376.100000,1223.400000,660.490000,496.100000 /)
XPIZA_LKT(45,7,1:6)=(/ 0.934673,0.981859,0.999502,0.999615,0.998368,0.571296 /)
XCGA_LKT(45,7,1:6)=(/ 0.739863,0.713273,0.722737,0.731057,0.685290,0.546453 /)
XEXT_COEFF_550_LKT(45,7)=1312.600000 !rg=0.352632 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,8,1:6)=(/ 940.500000,1026.200000,1132.000000,1081.400000,699.130000,534.060000 /)
XPIZA_LKT(45,8,1:6)=(/ 0.927796,0.978311,0.999374,0.999554,0.998387,0.598159 /)
XCGA_LKT(45,8,1:6)=(/ 0.755603,0.718247,0.717827,0.727883,0.707340,0.601483 /)
XEXT_COEFF_550_LKT(45,8)=1075.000000 !rg=0.352632 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,9,1:6)=(/ 787.340000,842.870000,920.060000,928.420000,691.730000,546.960000 /)
XPIZA_LKT(45,9,1:6)=(/ 0.917008,0.975148,0.999263,0.999471,0.998304,0.612918 /)
XCGA_LKT(45,9,1:6)=(/ 0.764490,0.731330,0.717893,0.724137,0.717313,0.642827 /)
XEXT_COEFF_550_LKT(45,9)=873.930000 !rg=0.352632 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,10,1:6)=(/ 654.160000,693.610000,752.670000,781.030000,644.540000,532.880000 /)
XPIZA_LKT(45,10,1:6)=(/ 0.906723,0.970556,0.999044,0.999358,0.998124,0.616948 /)
XCGA_LKT(45,10,1:6)=(/ 0.780803,0.742810,0.719040,0.719490,0.718603,0.673277 /)
XEXT_COEFF_550_LKT(45,10)=715.860000 !rg=0.352632 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,11,1:6)=(/ 539.680000,570.660000,614.870000,645.160000,579.170000,497.030000 /)
XPIZA_LKT(45,11,1:6)=(/ 0.893144,0.965302,0.998841,0.999234,0.997845,0.613172 /)
XCGA_LKT(45,11,1:6)=(/ 0.793067,0.749190,0.724553,0.719527,0.717797,0.697370 /)
XEXT_COEFF_550_LKT(45,11)=584.480000 !rg=0.352632 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,12,1:6)=(/ 447.460000,467.780000,498.340000,531.160000,505.790000,448.070000 /)
XPIZA_LKT(45,12,1:6)=(/ 0.878310,0.958952,0.998679,0.999041,0.997492,0.604942 /)
XCGA_LKT(45,12,1:6)=(/ 0.805190,0.763413,0.730360,0.719703,0.718747,0.718723 /)
XEXT_COEFF_550_LKT(45,12)=479.370000 !rg=0.352632 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,13,1:6)=(/ 369.250000,383.390000,403.280000,432.880000,429.630000,393.650000 /)
XPIZA_LKT(45,13,1:6)=(/ 0.864584,0.951839,0.998494,0.998831,0.997020,0.594832 /)
XCGA_LKT(45,13,1:6)=(/ 0.814443,0.774117,0.743513,0.721680,0.719240,0.739207 /)
XEXT_COEFF_550_LKT(45,13)=390.440000 !rg=0.352632 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,14,1:6)=(/ 305.510000,314.950000,331.820000,355.280000,359.820000,338.870000 /)
XPIZA_LKT(45,14,1:6)=(/ 0.850859,0.944789,0.998117,0.998339,0.996469,0.583760 /)
XCGA_LKT(45,14,1:6)=(/ 0.825460,0.783330,0.748927,0.726787,0.719900,0.758363 /)
XEXT_COEFF_550_LKT(45,14)=321.780000 !rg=0.352632 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,15,1:6)=(/ 251.200000,260.220000,271.200000,287.280000,298.360000,287.220000 /)
XPIZA_LKT(45,15,1:6)=(/ 0.839932,0.936408,0.997726,0.998189,0.995586,0.572728 /)
XCGA_LKT(45,15,1:6)=(/ 0.833137,0.793597,0.754540,0.734930,0.720463,0.776707 /)
XEXT_COEFF_550_LKT(45,15)=264.110000 !rg=0.352632 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,16,1:6)=(/ 206.540000,212.010000,221.480000,233.300000,246.910000,240.430000 /)
XPIZA_LKT(45,16,1:6)=(/ 0.829173,0.930734,0.997015,0.998026,0.994519,0.562478 /)
XCGA_LKT(45,16,1:6)=(/ 0.842713,0.801207,0.765257,0.743537,0.724707,0.794743 /)
XEXT_COEFF_550_LKT(45,16)=214.450000 !rg=0.352632 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,17,1:6)=(/ 170.930000,174.820000,181.220000,191.380000,200.580000,199.950000 /)
XPIZA_LKT(45,17,1:6)=(/ 0.819656,0.924337,0.996792,0.997655,0.993995,0.553600 /)
XCGA_LKT(45,17,1:6)=(/ 0.849183,0.810213,0.775777,0.751967,0.733100,0.811817 /)
XEXT_COEFF_550_LKT(45,17)=177.500000 !rg=0.352632 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,18,1:6)=(/ 141.190000,144.900000,148.610000,156.210000,164.890000,165.440000 /)
XPIZA_LKT(45,18,1:6)=(/ 0.810517,0.918148,0.996535,0.996797,0.992329,0.546251 /)
XCGA_LKT(45,18,1:6)=(/ 0.854833,0.817867,0.780537,0.757030,0.735690,0.827960 /)
XEXT_COEFF_550_LKT(45,18)=146.580000 !rg=0.352632 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,19,1:6)=(/ 116.490000,118.440000,122.220000,127.430000,135.250000,136.190000 /)
XPIZA_LKT(45,19,1:6)=(/ 0.801576,0.913790,0.996117,0.996553,0.990800,0.540370 /)
XCGA_LKT(45,19,1:6)=(/ 0.860927,0.823520,0.788527,0.767453,0.743857,0.843200 /)
XEXT_COEFF_550_LKT(45,19)=119.500000 !rg=0.352632 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,20,1:6)=(/ 96.419000,97.918000,100.540000,104.600000,109.610000,111.920000 /)
XPIZA_LKT(45,20,1:6)=(/ 0.792627,0.908262,0.995552,0.996169,0.989708,0.536032 /)
XCGA_LKT(45,20,1:6)=(/ 0.865390,0.829597,0.795813,0.776033,0.755420,0.857113 /)
XEXT_COEFF_550_LKT(45,20)=98.825000 !rg=0.352632 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,1,1:6)=(/ 1637.800000,3084.300000,2673.200000,1224.500000,242.860000,233.980000 /)
XPIZA_LKT(46,1,1:6)=(/ 0.950450,0.992617,0.999762,0.999656,0.996745,0.272643 /)
XCGA_LKT(46,1,1:6)=(/ 0.575677,0.763000,0.761143,0.643603,0.288887,0.103343 /)
XEXT_COEFF_550_LKT(46,1)=3213.300000 !rg=0.382026 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,2,1:6)=(/ 1655.400000,2750.800000,2628.200000,1253.400000,274.670000,252.460000 /)
XPIZA_LKT(46,2,1:6)=(/ 0.951104,0.991799,0.999756,0.999663,0.997022,0.310644 /)
XCGA_LKT(46,2,1:6)=(/ 0.610620,0.746487,0.767047,0.664080,0.362807,0.132433 /)
XEXT_COEFF_550_LKT(46,2)=2961.700000 !rg=0.382026 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,3,1:6)=(/ 1635.900000,2314.400000,2453.900000,1334.800000,342.420000,288.960000 /)
XPIZA_LKT(46,3,1:6)=(/ 0.951234,0.990161,0.999735,0.999677,0.997447,0.373016 /)
XCGA_LKT(46,3,1:6)=(/ 0.657710,0.724273,0.764140,0.693850,0.470167,0.201690 /)
XEXT_COEFF_550_LKT(46,3)=2558.900000 !rg=0.382026 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,4,1:6)=(/ 1535.600000,1927.200000,2169.100000,1394.500000,436.210000,342.520000 /)
XPIZA_LKT(46,4,1:6)=(/ 0.948683,0.988090,0.999697,0.999683,0.997855,0.441425 /)
XCGA_LKT(46,4,1:6)=(/ 0.690950,0.709060,0.754463,0.717073,0.556667,0.309097 /)
XEXT_COEFF_550_LKT(46,4)=2122.200000 !rg=0.382026 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,5,1:6)=(/ 1381.300000,1604.600000,1841.100000,1390.200000,537.420000,406.400000 /)
XPIZA_LKT(46,5,1:6)=(/ 0.943437,0.985846,0.999636,0.999673,0.998148,0.502756 /)
XCGA_LKT(46,5,1:6)=(/ 0.714263,0.705567,0.740523,0.730613,0.620723,0.416657 /)
XEXT_COEFF_550_LKT(46,5)=1748.600000 !rg=0.382026 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,6,1:6)=(/ 1190.900000,1347.400000,1527.000000,1312.000000,627.200000,468.420000 /)
XPIZA_LKT(46,6,1:6)=(/ 0.938602,0.982769,0.999556,0.999644,0.998324,0.550726 /)
XCGA_LKT(46,6,1:6)=(/ 0.732570,0.708363,0.728687,0.734470,0.666390,0.503090 /)
XEXT_COEFF_550_LKT(46,6)=1433.100000 !rg=0.382026 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,7,1:6)=(/ 1017.100000,1123.100000,1253.700000,1179.000000,689.140000,518.310000 /)
XPIZA_LKT(46,7,1:6)=(/ 0.930516,0.979564,0.999451,0.999595,0.998398,0.585312 /)
XCGA_LKT(46,7,1:6)=(/ 0.745747,0.714390,0.719687,0.731950,0.696887,0.569753 /)
XEXT_COEFF_550_LKT(46,7)=1178.400000 !rg=0.382026 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,8,1:6)=(/ 856.460000,927.820000,1027.400000,1025.800000,707.710000,546.000000 /)
XPIZA_LKT(46,8,1:6)=(/ 0.922777,0.976285,0.999329,0.999522,0.998375,0.606756 /)
XCGA_LKT(46,8,1:6)=(/ 0.762150,0.725530,0.716007,0.726810,0.713690,0.619500 /)
XEXT_COEFF_550_LKT(46,8)=967.770000 !rg=0.382026 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,9,1:6)=(/ 721.270000,765.220000,840.810000,868.420000,681.320000,547.120000 /)
XPIZA_LKT(46,9,1:6)=(/ 0.911053,0.972172,0.999181,0.999424,0.998250,0.616805 /)
XCGA_LKT(46,9,1:6)=(/ 0.774250,0.732707,0.716303,0.723167,0.719460,0.656560 /)
XEXT_COEFF_550_LKT(46,9)=798.830000 !rg=0.382026 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,10,1:6)=(/ 598.730000,630.740000,683.490000,721.450000,623.130000,522.690000 /)
XPIZA_LKT(46,10,1:6)=(/ 0.899394,0.968372,0.998921,0.999301,0.998033,0.617203 /)
XCGA_LKT(46,10,1:6)=(/ 0.784667,0.746990,0.721973,0.718897,0.719650,0.684117 /)
XEXT_COEFF_550_LKT(46,10)=650.470000 !rg=0.382026 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,11,1:6)=(/ 496.800000,516.960000,556.960000,597.100000,552.000000,479.870000 /)
XPIZA_LKT(46,11,1:6)=(/ 0.885593,0.963400,0.998784,0.999141,0.997723,0.611134 /)
XCGA_LKT(46,11,1:6)=(/ 0.799153,0.759100,0.728203,0.718983,0.718913,0.706877 /)
XEXT_COEFF_550_LKT(46,11)=533.330000 !rg=0.382026 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,12,1:6)=(/ 410.400000,426.500000,452.090000,486.120000,476.240000,427.390000 /)
XPIZA_LKT(46,12,1:6)=(/ 0.871401,0.956566,0.998603,0.998958,0.997335,0.601784 /)
XCGA_LKT(46,12,1:6)=(/ 0.808953,0.769293,0.737957,0.719150,0.719760,0.727860 /)
XEXT_COEFF_550_LKT(46,12)=436.690000 !rg=0.382026 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,13,1:6)=(/ 339.430000,350.490000,370.470000,398.100000,399.780000,371.650000 /)
XPIZA_LKT(46,13,1:6)=(/ 0.857225,0.948822,0.998325,0.998627,0.996810,0.590916 /)
XCGA_LKT(46,13,1:6)=(/ 0.820920,0.778047,0.744827,0.724210,0.718650,0.747877 /)
XEXT_COEFF_550_LKT(46,13)=358.950000 !rg=0.382026 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,14,1:6)=(/ 279.830000,289.840000,303.720000,322.410000,333.600000,317.320000 /)
XPIZA_LKT(46,14,1:6)=(/ 0.845477,0.940384,0.997747,0.998456,0.996188,0.579428 /)
XCGA_LKT(46,14,1:6)=(/ 0.828773,0.788743,0.748800,0.730313,0.720460,0.766603 /)
XEXT_COEFF_550_LKT(46,14)=294.590000 !rg=0.382026 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,15,1:6)=(/ 230.360000,237.430000,249.560000,262.180000,276.370000,267.360000 /)
XPIZA_LKT(46,15,1:6)=(/ 0.833819,0.934041,0.997342,0.998158,0.995194,0.568498 /)
XCGA_LKT(46,15,1:6)=(/ 0.839400,0.796110,0.758463,0.738597,0.722613,0.784790 /)
XEXT_COEFF_550_LKT(46,15)=241.390000 !rg=0.382026 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,16,1:6)=(/ 190.340000,193.890000,202.510000,214.940000,225.760000,222.750000 /)
XPIZA_LKT(46,16,1:6)=(/ 0.824129,0.927698,0.996996,0.997623,0.993989,0.558506 /)
XCGA_LKT(46,16,1:6)=(/ 0.846670,0.808233,0.770403,0.746103,0.726883,0.802583 /)
XEXT_COEFF_550_LKT(46,16)=197.680000 !rg=0.382026 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,17,1:6)=(/ 157.430000,160.610000,166.690000,174.340000,184.250000,184.630000 /)
XPIZA_LKT(46,17,1:6)=(/ 0.814857,0.922072,0.996390,0.997445,0.993520,0.550111 /)
XCGA_LKT(46,17,1:6)=(/ 0.851767,0.814843,0.774530,0.753160,0.734323,0.819367 /)
XEXT_COEFF_550_LKT(46,17)=162.350000 !rg=0.382026 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,18,1:6)=(/ 129.850000,132.440000,137.090000,142.380000,151.600000,152.390000 /)
XPIZA_LKT(46,18,1:6)=(/ 0.805552,0.916111,0.996180,0.996845,0.991970,0.543328 /)
XCGA_LKT(46,18,1:6)=(/ 0.858817,0.820383,0.783223,0.763343,0.741293,0.835117 /)
XEXT_COEFF_550_LKT(46,18)=134.020000 !rg=0.382026 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,19,1:6)=(/ 107.430000,108.720000,112.290000,117.390000,123.500000,125.240000 /)
XPIZA_LKT(46,19,1:6)=(/ 0.797399,0.910922,0.995863,0.996259,0.989832,0.538048 /)
XCGA_LKT(46,19,1:6)=(/ 0.863730,0.828077,0.792463,0.770533,0.746277,0.849873 /)
XEXT_COEFF_550_LKT(46,19)=110.370000 !rg=0.382026 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,20,1:6)=(/ 88.786000,90.091000,92.540000,95.302000,100.580000,102.810000 /)
XPIZA_LKT(46,20,1:6)=(/ 0.788022,0.906306,0.995486,0.995926,0.988726,0.534305 /)
XCGA_LKT(46,20,1:6)=(/ 0.867137,0.832883,0.795637,0.777890,0.756647,0.863223 /)
XEXT_COEFF_550_LKT(46,20)=90.670000 !rg=0.382026 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,1,1:6)=(/ 1357.800000,2626.900000,2714.400000,1301.000000,284.960000,255.620000 /)
XPIZA_LKT(47,1,1:6)=(/ 0.939930,0.991440,0.999761,0.999681,0.997122,0.313816 /)
XCGA_LKT(47,1,1:6)=(/ 0.537567,0.739633,0.774913,0.652677,0.346790,0.122017 /)
XEXT_COEFF_550_LKT(47,1)=2948.700000 !rg=0.41387 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,2,1:6)=(/ 1454.200000,2341.600000,2582.200000,1350.300000,321.100000,276.890000 /)
XPIZA_LKT(47,2,1:6)=(/ 0.944783,0.990291,0.999749,0.999682,0.997332,0.351920 /)
XCGA_LKT(47,2,1:6)=(/ 0.610123,0.723857,0.772263,0.685680,0.425393,0.156777 /)
XEXT_COEFF_550_LKT(47,2)=2691.600000 !rg=0.41387 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,3,1:6)=(/ 1485.000000,1973.700000,2350.700000,1414.800000,396.820000,317.560000 /)
XPIZA_LKT(47,3,1:6)=(/ 0.946943,0.988490,0.999721,0.999689,0.997692,0.410892 /)
XCGA_LKT(47,3,1:6)=(/ 0.671177,0.704330,0.764807,0.711870,0.518123,0.236783 /)
XEXT_COEFF_550_LKT(47,3)=2284.000000 !rg=0.41387 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,4,1:6)=(/ 1390.200000,1672.600000,2028.000000,1439.500000,492.400000,374.740000 /)
XPIZA_LKT(47,4,1:6)=(/ 0.945071,0.986358,0.999672,0.999687,0.998025,0.472915 /)
XCGA_LKT(47,4,1:6)=(/ 0.704043,0.696770,0.750973,0.729433,0.589210,0.348527 /)
XEXT_COEFF_550_LKT(47,4)=1885.700000 !rg=0.41387 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,5,1:6)=(/ 1247.300000,1429.700000,1689.700000,1393.200000,588.750000,438.650000 /)
XPIZA_LKT(47,5,1:6)=(/ 0.939200,0.983742,0.999602,0.999668,0.998255,0.527141 /)
XCGA_LKT(47,5,1:6)=(/ 0.721873,0.698607,0.735907,0.737503,0.643477,0.451207 /)
XEXT_COEFF_550_LKT(47,5)=1551.700000 !rg=0.41387 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,6,1:6)=(/ 1080.500000,1201.100000,1392.700000,1279.000000,666.710000,496.120000 /)
XPIZA_LKT(47,6,1:6)=(/ 0.934332,0.981097,0.999514,0.999630,0.998381,0.568508 /)
XCGA_LKT(47,6,1:6)=(/ 0.740193,0.707113,0.724327,0.737360,0.681983,0.530937 /)
XEXT_COEFF_550_LKT(47,6)=1283.800000 !rg=0.41387 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,7,1:6)=(/ 921.170000,1007.300000,1141.900000,1125.500000,710.410000,537.150000 /)
XPIZA_LKT(47,7,1:6)=(/ 0.926917,0.977983,0.999385,0.999571,0.998413,0.597157 /)
XCGA_LKT(47,7,1:6)=(/ 0.756717,0.715300,0.717203,0.732513,0.706567,0.591387 /)
XEXT_COEFF_550_LKT(47,7)=1058.600000 !rg=0.41387 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,8,1:6)=(/ 781.070000,837.310000,930.720000,962.670000,708.500000,553.380000 /)
XPIZA_LKT(47,8,1:6)=(/ 0.917112,0.974765,0.999282,0.999489,0.998347,0.613485 /)
XCGA_LKT(47,8,1:6)=(/ 0.767157,0.729070,0.716933,0.726853,0.718540,0.636033 /)
XEXT_COEFF_550_LKT(47,8)=873.730000 !rg=0.41387 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,9,1:6)=(/ 656.320000,700.360000,762.730000,803.980000,664.780000,542.790000 /)
XPIZA_LKT(47,9,1:6)=(/ 0.905361,0.969367,0.999026,0.999383,0.998186,0.619195 /)
XCGA_LKT(47,9,1:6)=(/ 0.778547,0.738923,0.714550,0.721640,0.721567,0.669237 /)
XEXT_COEFF_550_LKT(47,9)=722.900000 !rg=0.41387 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,10,1:6)=(/ 549.800000,574.860000,625.380000,668.530000,598.780000,509.190000 /)
XPIZA_LKT(47,10,1:6)=(/ 0.892257,0.965685,0.998879,0.999207,0.997922,0.616435 /)
XCGA_LKT(47,10,1:6)=(/ 0.793477,0.749763,0.721333,0.718803,0.719860,0.694490 /)
XEXT_COEFF_550_LKT(47,10)=597.220000 !rg=0.41387 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,11,1:6)=(/ 455.260000,474.890000,506.860000,546.530000,524.020000,460.800000 /)
XPIZA_LKT(47,11,1:6)=(/ 0.877814,0.959981,0.998664,0.999066,0.997582,0.608553 /)
XCGA_LKT(47,11,1:6)=(/ 0.802500,0.762787,0.729827,0.717040,0.720257,0.716347 /)
XEXT_COEFF_550_LKT(47,11)=486.840000 !rg=0.41387 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,12,1:6)=(/ 377.290000,390.060000,415.700000,448.110000,445.300000,405.770000 /)
XPIZA_LKT(47,12,1:6)=(/ 0.864093,0.953165,0.998430,0.998789,0.997123,0.598231 /)
XCGA_LKT(47,12,1:6)=(/ 0.816080,0.772420,0.737953,0.720873,0.718917,0.736917 /)
XEXT_COEFF_550_LKT(47,12)=400.850000 !rg=0.41387 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,13,1:6)=(/ 310.950000,322.300000,338.950000,361.890000,371.620000,349.440000 /)
XPIZA_LKT(47,13,1:6)=(/ 0.851281,0.944611,0.998054,0.998607,0.996567,0.586648 /)
XCGA_LKT(47,13,1:6)=(/ 0.824210,0.783963,0.744647,0.726597,0.719317,0.756383 /)
XEXT_COEFF_550_LKT(47,13)=328.460000 !rg=0.41387 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,14,1:6)=(/ 256.460000,265.100000,277.850000,293.870000,309.850000,296.220000 /)
XPIZA_LKT(47,14,1:6)=(/ 0.839470,0.937364,0.997801,0.997988,0.995620,0.574995 /)
XCGA_LKT(47,14,1:6)=(/ 0.835093,0.790840,0.756707,0.735403,0.720507,0.774843 /)
XEXT_COEFF_550_LKT(47,14)=269.990000 !rg=0.41387 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,15,1:6)=(/ 212.010000,217.260000,226.890000,240.860000,253.740000,248.240000 /)
XPIZA_LKT(47,15,1:6)=(/ 0.829334,0.930775,0.997045,0.998060,0.995007,0.564236 /)
XCGA_LKT(47,15,1:6)=(/ 0.843340,0.803310,0.761697,0.743093,0.723677,0.792840 /)
XEXT_COEFF_550_LKT(47,15)=221.080000 !rg=0.41387 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,16,1:6)=(/ 174.770000,179.140000,185.430000,195.700000,206.430000,205.980000 /)
XPIZA_LKT(47,16,1:6)=(/ 0.819264,0.924797,0.996921,0.997688,0.993991,0.554635 /)
XCGA_LKT(47,16,1:6)=(/ 0.849030,0.811377,0.773333,0.747047,0.731407,0.810380 /)
XEXT_COEFF_550_LKT(47,16)=181.960000 !rg=0.41387 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,17,1:6)=(/ 144.420000,147.770000,152.840000,158.490000,170.220000,170.230000 /)
XPIZA_LKT(47,17,1:6)=(/ 0.810504,0.918585,0.996459,0.997240,0.992460,0.546780 /)
XCGA_LKT(47,17,1:6)=(/ 0.855490,0.816410,0.781350,0.762177,0.736377,0.826820 /)
XEXT_COEFF_550_LKT(47,17)=149.890000 !rg=0.41387 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,18,1:6)=(/ 119.540000,121.670000,125.200000,130.730000,138.390000,140.230000 /)
XPIZA_LKT(47,18,1:6)=(/ 0.801877,0.913478,0.995824,0.996688,0.991443,0.540602 /)
XCGA_LKT(47,18,1:6)=(/ 0.861370,0.824783,0.785497,0.768343,0.741773,0.842130 /)
XEXT_COEFF_550_LKT(47,18)=123.410000 !rg=0.41387 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,19,1:6)=(/ 98.653000,100.490000,103.100000,107.130000,112.430000,115.090000 /)
XPIZA_LKT(47,19,1:6)=(/ 0.792681,0.908544,0.995891,0.996050,0.989903,0.535949 /)
XCGA_LKT(47,19,1:6)=(/ 0.865293,0.830623,0.794520,0.772940,0.753013,0.856350 /)
XEXT_COEFF_550_LKT(47,19)=101.630000 !rg=0.41387 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,20,1:6)=(/ 81.578000,82.941000,84.950000,87.077000,92.644000,94.394000 /)
XPIZA_LKT(47,20,1:6)=(/ 0.783626,0.903012,0.995397,0.995757,0.987345,0.532806 /)
XCGA_LKT(47,20,1:6)=(/ 0.869530,0.834053,0.801573,0.785337,0.759783,0.869103 /)
XEXT_COEFF_550_LKT(47,20)=83.782000 !rg=0.41387 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,1,1:6)=(/ 1251.400000,2178.400000,2615.100000,1368.100000,329.950000,281.210000 /)
XPIZA_LKT(48,1,1:6)=(/ 0.936513,0.989502,0.999753,0.999691,0.997401,0.356310 /)
XCGA_LKT(48,1,1:6)=(/ 0.576413,0.712030,0.780967,0.683433,0.416070,0.144520 /)
XEXT_COEFF_550_LKT(48,1)=2651.500000 !rg=0.448369 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,2,1:6)=(/ 1378.000000,1933.900000,2482.600000,1438.500000,374.150000,305.100000 /)
XPIZA_LKT(48,2,1:6)=(/ 0.941955,0.988202,0.999735,0.999695,0.997585,0.392875 /)
XCGA_LKT(48,2,1:6)=(/ 0.645393,0.692600,0.773903,0.709083,0.489423,0.186167 /)
XEXT_COEFF_550_LKT(48,2)=2372.200000 !rg=0.448369 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,3,1:6)=(/ 1381.500000,1677.800000,2204.100000,1478.700000,455.920000,349.530000 /)
XPIZA_LKT(48,3,1:6)=(/ 0.942755,0.986362,0.999699,0.999697,0.997899,0.446600 /)
XCGA_LKT(48,3,1:6)=(/ 0.689237,0.682990,0.761867,0.727840,0.560477,0.276833 /)
XEXT_COEFF_550_LKT(48,3)=1992.600000 !rg=0.448369 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,4,1:6)=(/ 1278.200000,1461.700000,1860.600000,1464.700000,549.670000,409.000000 /)
XPIZA_LKT(48,4,1:6)=(/ 0.940059,0.984440,0.999639,0.999687,0.998166,0.501736 /)
XCGA_LKT(48,4,1:6)=(/ 0.716140,0.687347,0.744573,0.739493,0.617963,0.389040 /)
XEXT_COEFF_550_LKT(48,4)=1652.900000 !rg=0.448369 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,5,1:6)=(/ 1130.400000,1263.300000,1537.100000,1376.800000,637.680000,470.760000 /)
XPIZA_LKT(48,5,1:6)=(/ 0.936330,0.981858,0.999555,0.999660,0.998341,0.548976 /)
XCGA_LKT(48,5,1:6)=(/ 0.737327,0.695447,0.728670,0.742460,0.663670,0.484660 /)
XEXT_COEFF_550_LKT(48,5)=1378.400000 !rg=0.448369 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,6,1:6)=(/ 993.160000,1076.700000,1260.000000,1233.200000,700.520000,521.390000 /)
XPIZA_LKT(48,6,1:6)=(/ 0.928329,0.979504,0.999457,0.999611,0.998422,0.584025 /)
XCGA_LKT(48,6,1:6)=(/ 0.750497,0.709073,0.717553,0.738237,0.695527,0.557103 /)
XEXT_COEFF_550_LKT(48,6)=1153.600000 !rg=0.448369 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,7,1:6)=(/ 846.580000,905.180000,1030.100000,1065.300000,723.870000,551.980000 /)
XPIZA_LKT(48,7,1:6)=(/ 0.920651,0.976663,0.999331,0.999540,0.998411,0.607040 /)
XCGA_LKT(48,7,1:6)=(/ 0.764240,0.723227,0.713877,0.731080,0.714470,0.611470 /)
XEXT_COEFF_550_LKT(48,7)=954.340000 !rg=0.448369 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,8,1:6)=(/ 718.650000,758.720000,845.800000,898.830000,701.570000,556.110000 /)
XPIZA_LKT(48,8,1:6)=(/ 0.910170,0.972675,0.999188,0.999440,0.998303,0.618481 /)
XCGA_LKT(48,8,1:6)=(/ 0.776620,0.732257,0.713470,0.724727,0.721953,0.651317 /)
XEXT_COEFF_550_LKT(48,8)=795.210000 !rg=0.448369 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,9,1:6)=(/ 596.760000,634.320000,695.480000,742.320000,644.840000,534.330000 /)
XPIZA_LKT(48,9,1:6)=(/ 0.900083,0.967597,0.998950,0.999330,0.998097,0.620241 /)
XCGA_LKT(48,9,1:6)=(/ 0.789640,0.740630,0.716923,0.720693,0.722400,0.681137 /)
XEXT_COEFF_550_LKT(48,9)=653.430000 !rg=0.448369 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,10,1:6)=(/ 501.650000,528.230000,567.260000,612.150000,571.880000,492.860000 /)
XPIZA_LKT(48,10,1:6)=(/ 0.885556,0.961931,0.998780,0.999127,0.997802,0.614830 /)
XCGA_LKT(48,10,1:6)=(/ 0.797020,0.756033,0.721700,0.718593,0.721373,0.704593 /)
XEXT_COEFF_550_LKT(48,10)=541.520000 !rg=0.448369 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,11,1:6)=(/ 417.750000,434.660000,465.210000,500.040000,494.350000,440.210000 /)
XPIZA_LKT(48,11,1:6)=(/ 0.871229,0.955691,0.998515,0.998977,0.997386,0.605482 /)
XCGA_LKT(48,11,1:6)=(/ 0.810833,0.764607,0.733033,0.721280,0.719430,0.725770 /)
XEXT_COEFF_550_LKT(48,11)=447.290000 !rg=0.448369 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,12,1:6)=(/ 344.980000,359.010000,379.100000,406.920000,414.330000,383.450000 /)
XPIZA_LKT(48,12,1:6)=(/ 0.858015,0.948368,0.998116,0.998605,0.996912,0.594243 /)
XCGA_LKT(48,12,1:6)=(/ 0.819457,0.778380,0.738547,0.723110,0.718963,0.745813 /)
XEXT_COEFF_550_LKT(48,12)=365.820000 !rg=0.448369 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,13,1:6)=(/ 284.920000,294.890000,310.120000,329.080000,346.460000,327.310000 /)
XPIZA_LKT(48,13,1:6)=(/ 0.845401,0.941069,0.997981,0.998581,0.996135,0.582111 /)
XCGA_LKT(48,13,1:6)=(/ 0.830863,0.785710,0.752107,0.733053,0.719753,0.764803 /)
XEXT_COEFF_550_LKT(48,13)=301.370000 !rg=0.448369 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,14,1:6)=(/ 235.700000,242.430000,252.590000,270.730000,284.470000,275.730000 /)
XPIZA_LKT(48,14,1:6)=(/ 0.834296,0.933961,0.997542,0.998144,0.995433,0.570496 /)
XCGA_LKT(48,14,1:6)=(/ 0.839883,0.797903,0.758773,0.736123,0.719817,0.783093 /)
XEXT_COEFF_550_LKT(48,14)=246.580000 !rg=0.448369 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,15,1:6)=(/ 195.050000,199.490000,206.480000,220.230000,232.110000,229.950000 /)
XPIZA_LKT(48,15,1:6)=(/ 0.823958,0.927899,0.996863,0.997739,0.994797,0.559994 /)
XCGA_LKT(48,15,1:6)=(/ 0.845843,0.808083,0.769767,0.741290,0.729830,0.800850 /)
XEXT_COEFF_550_LKT(48,15)=202.680000 !rg=0.448369 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,16,1:6)=(/ 160.900000,164.400000,170.400000,178.670000,190.870000,190.130000 /)
XPIZA_LKT(48,16,1:6)=(/ 0.814796,0.921461,0.996683,0.997362,0.993411,0.550893 /)
XCGA_LKT(48,16,1:6)=(/ 0.852910,0.813277,0.777060,0.754643,0.732760,0.818113 /)
XEXT_COEFF_550_LKT(48,16)=166.900000 !rg=0.448369 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,17,1:6)=(/ 133.090000,135.720000,139.380000,146.080000,155.500000,156.760000 /)
XPIZA_LKT(48,17,1:6)=(/ 0.805532,0.916247,0.996381,0.996936,0.992290,0.543634 /)
XCGA_LKT(48,17,1:6)=(/ 0.858807,0.820997,0.784090,0.764033,0.736627,0.834157 /)
XEXT_COEFF_550_LKT(48,17)=136.770000 !rg=0.448369 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,18,1:6)=(/ 110.000000,111.910000,114.420000,119.830000,126.000000,128.910000 /)
XPIZA_LKT(48,18,1:6)=(/ 0.797259,0.911295,0.995881,0.996293,0.990993,0.538099 /)
XCGA_LKT(48,18,1:6)=(/ 0.863090,0.828390,0.792423,0.767813,0.750310,0.848973 /)
XEXT_COEFF_550_LKT(48,18)=113.220000 !rg=0.448369 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,19,1:6)=(/ 90.939000,92.151000,94.612000,98.034000,103.650000,105.690000 /)
XPIZA_LKT(48,19,1:6)=(/ 0.788519,0.905781,0.995715,0.995793,0.988978,0.534087 /)
XCGA_LKT(48,19,1:6)=(/ 0.867843,0.832190,0.797810,0.778393,0.756167,0.862617 /)
XEXT_COEFF_550_LKT(48,19)=93.253000 !rg=0.448369 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,20,1:6)=(/ 75.306000,76.353000,77.708000,80.376000,84.693000,86.632000 /)
XPIZA_LKT(48,20,1:6)=(/ 0.778836,0.901092,0.995319,0.995547,0.986566,0.531533 /)
XCGA_LKT(48,20,1:6)=(/ 0.871947,0.837473,0.804023,0.787380,0.761137,0.874740 /)
XEXT_COEFF_550_LKT(48,20)=76.708000 !rg=0.448369 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,1,1:6)=(/ 1305.800000,1711.100000,2513.000000,1481.500000,380.540000,310.960000 /)
XPIZA_LKT(49,1,1:6)=(/ 0.939539,0.986709,0.999734,0.999696,0.997613,0.398770 /)
XCGA_LKT(49,1,1:6)=(/ 0.657107,0.663477,0.778513,0.724357,0.493450,0.171887 /)
XEXT_COEFF_550_LKT(49,1)=2275.100000 !rg=0.485743 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,2,1:6)=(/ 1368.600000,1571.500000,2329.500000,1515.300000,435.490000,337.050000 /)
XPIZA_LKT(49,2,1:6)=(/ 0.941147,0.985445,0.999715,0.999704,0.997805,0.432109 /)
XCGA_LKT(49,2,1:6)=(/ 0.689213,0.657663,0.771263,0.730217,0.547260,0.221730 /)
XEXT_COEFF_550_LKT(49,2)=2028.000000 !rg=0.485743 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,3,1:6)=(/ 1293.800000,1424.800000,2020.100000,1523.200000,518.280000,384.620000 /)
XPIZA_LKT(49,3,1:6)=(/ 0.939714,0.983817,0.999669,0.999700,0.998074,0.479482 /)
XCGA_LKT(49,3,1:6)=(/ 0.713650,0.660193,0.755097,0.741303,0.596257,0.321370 /)
XEXT_COEFF_550_LKT(49,3)=1705.700000 !rg=0.485743 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,4,1:6)=(/ 1173.600000,1284.100000,1679.500000,1467.800000,606.400000,444.550000 /)
XPIZA_LKT(49,4,1:6)=(/ 0.935643,0.982141,0.999598,0.999683,0.998282,0.527788 /)
XCGA_LKT(49,4,1:6)=(/ 0.729167,0.678027,0.736023,0.747177,0.643187,0.429613 /)
XEXT_COEFF_550_LKT(49,4)=1439.900000 !rg=0.485743 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,5,1:6)=(/ 1039.400000,1125.900000,1378.200000,1341.600000,682.490000,501.740000 /)
XPIZA_LKT(49,5,1:6)=(/ 0.930890,0.980262,0.999505,0.999646,0.998408,0.568355 /)
XCGA_LKT(49,5,1:6)=(/ 0.747570,0.698990,0.720377,0.745107,0.681390,0.516627 /)
XEXT_COEFF_550_LKT(49,5)=1219.000000 !rg=0.485743 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,6,1:6)=(/ 906.380000,975.850000,1133.000000,1173.100000,727.150000,543.440000 /)
XPIZA_LKT(49,6,1:6)=(/ 0.922472,0.977052,0.999396,0.999588,0.998446,0.597332 /)
XCGA_LKT(49,6,1:6)=(/ 0.754223,0.710910,0.712923,0.738157,0.707000,0.581563 /)
XEXT_COEFF_550_LKT(49,6)=1029.000000 !rg=0.485743 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,7,1:6)=(/ 773.480000,823.540000,929.100000,995.570000,728.960000,562.200000 /)
XPIZA_LKT(49,7,1:6)=(/ 0.914159,0.974042,0.999263,0.999506,0.998393,0.614989 /)
XCGA_LKT(49,7,1:6)=(/ 0.768030,0.725820,0.711430,0.730137,0.720763,0.630013 /)
XEXT_COEFF_550_LKT(49,7)=859.460000 !rg=0.485743 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,8,1:6)=(/ 655.560000,695.450000,764.270000,829.250000,687.840000,554.070000 /)
XPIZA_LKT(49,8,1:6)=(/ 0.904370,0.969331,0.999047,0.999397,0.998247,0.621878 /)
XCGA_LKT(49,8,1:6)=(/ 0.779457,0.738237,0.712387,0.722753,0.724813,0.665443 /)
XEXT_COEFF_550_LKT(49,8)=719.490000 !rg=0.485743 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,9,1:6)=(/ 548.970000,573.800000,627.010000,685.330000,619.820000,522.090000 /)
XPIZA_LKT(49,9,1:6)=(/ 0.892856,0.965716,0.998796,0.999256,0.998003,0.620106 /)
XCGA_LKT(49,9,1:6)=(/ 0.796147,0.753050,0.718630,0.718473,0.723503,0.692410 /)
XEXT_COEFF_550_LKT(49,9)=591.860000 !rg=0.485743 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,10,1:6)=(/ 457.860000,480.290000,518.990000,559.080000,543.680000,474.210000 /)
XPIZA_LKT(49,10,1:6)=(/ 0.878813,0.959729,0.998332,0.999112,0.997663,0.612504 /)
XCGA_LKT(49,10,1:6)=(/ 0.807150,0.757423,0.725060,0.719713,0.721633,0.714527 /)
XEXT_COEFF_550_LKT(49,10)=493.560000 !rg=0.485743 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,11,1:6)=(/ 381.430000,398.090000,422.490000,455.900000,461.560000,418.390000 /)
XPIZA_LKT(49,11,1:6)=(/ 0.864678,0.952943,0.998423,0.998932,0.997217,0.601911 /)
XCGA_LKT(49,11,1:6)=(/ 0.815263,0.773133,0.732533,0.722547,0.720093,0.735090 /)
XEXT_COEFF_550_LKT(49,11)=406.100000 !rg=0.485743 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,12,1:6)=(/ 316.080000,327.780000,346.560000,370.290000,386.830000,360.740000 /)
XPIZA_LKT(49,12,1:6)=(/ 0.851491,0.945315,0.998134,0.998577,0.996493,0.589829 /)
XCGA_LKT(49,12,1:6)=(/ 0.827110,0.780387,0.745760,0.727247,0.718523,0.754540 /)
XEXT_COEFF_550_LKT(49,12)=335.230000 !rg=0.485743 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,13,1:6)=(/ 261.740000,269.870000,282.180000,302.680000,319.100000,305.560000 /)
XPIZA_LKT(49,13,1:6)=(/ 0.839233,0.937633,0.997806,0.998396,0.995875,0.577406 /)
XCGA_LKT(49,13,1:6)=(/ 0.836117,0.793290,0.753023,0.732663,0.718683,0.773213 /)
XEXT_COEFF_550_LKT(49,13)=274.490000 !rg=0.485743 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,14,1:6)=(/ 216.890000,222.470000,229.420000,246.740000,260.980000,255.960000 /)
XPIZA_LKT(49,14,1:6)=(/ 0.829342,0.930705,0.997461,0.997964,0.995236,0.565950 /)
XCGA_LKT(49,14,1:6)=(/ 0.842667,0.803217,0.768423,0.736977,0.726600,0.791323 /)
XEXT_COEFF_550_LKT(49,14)=225.030000 !rg=0.485743 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,15,1:6)=(/ 179.370000,183.500000,190.310000,201.720000,213.580000,212.570000 /)
XPIZA_LKT(49,15,1:6)=(/ 0.819121,0.924735,0.996830,0.997602,0.994253,0.555834 /)
XCGA_LKT(49,15,1:6)=(/ 0.850413,0.810133,0.772630,0.747580,0.729863,0.808827 /)
XEXT_COEFF_550_LKT(49,15)=186.300000 !rg=0.485743 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,16,1:6)=(/ 147.710000,151.250000,156.470000,163.030000,174.430000,175.240000 /)
XPIZA_LKT(49,16,1:6)=(/ 0.810471,0.918984,0.996251,0.997284,0.992807,0.547304 /)
XCGA_LKT(49,16,1:6)=(/ 0.855550,0.818323,0.777510,0.760560,0.732617,0.825753 /)
XEXT_COEFF_550_LKT(49,16)=152.360000 !rg=0.485743 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,17,1:6)=(/ 122.200000,124.620000,127.650000,134.050000,142.280000,144.180000 /)
XPIZA_LKT(49,17,1:6)=(/ 0.802021,0.913009,0.996116,0.996606,0.991347,0.540701 /)
XCGA_LKT(49,17,1:6)=(/ 0.860853,0.823417,0.789807,0.763910,0.743080,0.841350 /)
XEXT_COEFF_550_LKT(49,17)=125.240000 !rg=0.485743 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,18,1:6)=(/ 101.380000,103.120000,105.670000,110.060000,115.910000,118.420000 /)
XPIZA_LKT(49,18,1:6)=(/ 0.793126,0.908523,0.995738,0.996152,0.990159,0.535834 /)
XCGA_LKT(49,18,1:6)=(/ 0.866150,0.829733,0.794060,0.773060,0.752383,0.855627 /)
XEXT_COEFF_550_LKT(49,18)=103.740000 !rg=0.485743 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,19,1:6)=(/ 83.580000,84.977000,87.248000,89.733000,94.831000,97.004000 /)
XPIZA_LKT(49,19,1:6)=(/ 0.784301,0.903502,0.995469,0.995699,0.987966,0.532462 /)
XCGA_LKT(49,19,1:6)=(/ 0.869467,0.835233,0.798633,0.784747,0.756893,0.868650 /)
XEXT_COEFF_550_LKT(49,19)=85.451000 !rg=0.485743 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,20,1:6)=(/ 69.208000,70.238000,71.432000,73.897000,77.478000,79.481000 /)
XPIZA_LKT(49,20,1:6)=(/ 0.774563,0.897918,0.995227,0.995329,0.985427,0.530492 /)
XCGA_LKT(49,20,1:6)=(/ 0.873240,0.839063,0.807360,0.787920,0.767783,0.880110 /)
XEXT_COEFF_550_LKT(49,20)=70.436000 !rg=0.485743 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,1,1:6)=(/ 1400.700000,1345.000000,2330.300000,1592.300000,443.320000,344.710000 /)
XPIZA_LKT(50,1,1:6)=(/ 0.944352,0.982705,0.999717,0.999710,0.997799,0.439643 /)
XCGA_LKT(50,1,1:6)=(/ 0.733380,0.607647,0.774237,0.748783,0.567433,0.205513 /)
XEXT_COEFF_550_LKT(50,1)=1900.100000 !rg=0.526233 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,2,1:6)=(/ 1343.700000,1286.400000,2127.200000,1573.100000,504.600000,372.530000 /)
XPIZA_LKT(50,2,1:6)=(/ 0.941119,0.981916,0.999685,0.999709,0.998004,0.468349 /)
XCGA_LKT(50,2,1:6)=(/ 0.734197,0.616183,0.763613,0.746367,0.592220,0.264600 /)
XEXT_COEFF_550_LKT(50,2)=1678.100000 !rg=0.526233 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,3,1:6)=(/ 1218.700000,1227.300000,1810.300000,1544.900000,581.790000,422.350000 /)
XPIZA_LKT(50,3,1:6)=(/ 0.935034,0.981282,0.999626,0.999700,0.998221,0.509225 /)
XCGA_LKT(50,3,1:6)=(/ 0.738703,0.651203,0.743447,0.751900,0.625957,0.369203 /)
XEXT_COEFF_550_LKT(50,3)=1452.500000 !rg=0.526233 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,4,1:6)=(/ 1089.800000,1138.300000,1495.700000,1447.700000,660.750000,480.370000 /)
XPIZA_LKT(50,4,1:6)=(/ 0.929719,0.979942,0.999542,0.999674,0.998377,0.551119 /)
XCGA_LKT(50,4,1:6)=(/ 0.747587,0.678793,0.723167,0.752317,0.665297,0.469243 /)
XEXT_COEFF_550_LKT(50,4)=1262.000000 !rg=0.526233 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,5,1:6)=(/ 951.490000,1012.400000,1229.200000,1287.700000,721.580000,530.480000 /)
XPIZA_LKT(50,5,1:6)=(/ 0.925288,0.977971,0.999443,0.999627,0.998457,0.585353 /)
XCGA_LKT(50,5,1:6)=(/ 0.754300,0.700273,0.712680,0.746107,0.696837,0.546753 /)
XEXT_COEFF_550_LKT(50,5)=1081.900000 !rg=0.526233 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,6,1:6)=(/ 829.660000,879.340000,1021.100000,1104.500000,745.650000,561.430000 /)
XPIZA_LKT(50,6,1:6)=(/ 0.916746,0.974398,0.999299,0.999558,0.998453,0.608577 /)
XCGA_LKT(50,6,1:6)=(/ 0.769503,0.713427,0.706613,0.736693,0.716583,0.604320 /)
XEXT_COEFF_550_LKT(50,6)=930.210000 !rg=0.526233 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,7,1:6)=(/ 710.220000,746.660000,841.160000,924.620000,725.760000,567.510000 /)
XPIZA_LKT(50,7,1:6)=(/ 0.908030,0.971046,0.999161,0.999463,0.998358,0.621116 /)
XCGA_LKT(50,7,1:6)=(/ 0.780023,0.728463,0.708410,0.727673,0.725360,0.647177 /)
XEXT_COEFF_550_LKT(50,7)=783.320000 !rg=0.526233 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,8,1:6)=(/ 598.090000,631.270000,695.220000,762.310000,669.070000,547.450000 /)
XPIZA_LKT(50,8,1:6)=(/ 0.898005,0.967048,0.998980,0.999352,0.998167,0.623787 /)
XCGA_LKT(50,8,1:6)=(/ 0.791880,0.739063,0.714137,0.721767,0.726117,0.678633 /)
XEXT_COEFF_550_LKT(50,8)=655.050000 !rg=0.526233 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,9,1:6)=(/ 503.870000,524.560000,566.340000,626.220000,592.810000,506.540000 /)
XPIZA_LKT(50,9,1:6)=(/ 0.884701,0.962802,0.998835,0.999176,0.997888,0.618943 /)
XCGA_LKT(50,9,1:6)=(/ 0.799693,0.758110,0.723877,0.716943,0.724543,0.703253 /)
XEXT_COEFF_550_LKT(50,9)=538.380000 !rg=0.526233 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,10,1:6)=(/ 420.660000,437.140000,467.550000,512.910000,511.250000,453.640000 /)
XPIZA_LKT(50,10,1:6)=(/ 0.871979,0.956352,0.998585,0.999012,0.997505,0.609527 /)
XCGA_LKT(50,10,1:6)=(/ 0.813520,0.768827,0.730237,0.718193,0.722347,0.724317 /)
XEXT_COEFF_550_LKT(50,10)=447.290000 !rg=0.526233 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,11,1:6)=(/ 348.670000,362.370000,385.980000,415.670000,431.190000,395.660000 /)
XPIZA_LKT(50,11,1:6)=(/ 0.858328,0.949660,0.998178,0.998840,0.996839,0.597818 /)
XCGA_LKT(50,11,1:6)=(/ 0.823230,0.775957,0.739483,0.722943,0.718997,0.744237 /)
XEXT_COEFF_550_LKT(50,11)=368.420000 !rg=0.526233 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,12,1:6)=(/ 290.370000,299.670000,314.380000,340.670000,357.410000,337.980000 /)
XPIZA_LKT(50,12,1:6)=(/ 0.845612,0.941388,0.997858,0.998490,0.996304,0.585061 /)
XCGA_LKT(50,12,1:6)=(/ 0.832510,0.788907,0.748967,0.727497,0.718360,0.763160 /)
XEXT_COEFF_550_LKT(50,12)=305.260000 !rg=0.526233 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,13,1:6)=(/ 240.540000,247.190000,256.010000,276.170000,292.850000,284.360000 /)
XPIZA_LKT(50,13,1:6)=(/ 0.834378,0.933995,0.997637,0.998163,0.995637,0.572598 /)
XCGA_LKT(50,13,1:6)=(/ 0.839297,0.798187,0.763390,0.732840,0.723640,0.781637 /)
XEXT_COEFF_550_LKT(50,13)=249.590000 !rg=0.526233 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,14,1:6)=(/ 199.290000,204.010000,211.740000,226.400000,239.770000,237.020000 /)
XPIZA_LKT(50,14,1:6)=(/ 0.823407,0.928378,0.997142,0.997509,0.994733,0.561394 /)
XCGA_LKT(50,14,1:6)=(/ 0.847323,0.806170,0.770277,0.742057,0.727297,0.799520 /)
XEXT_COEFF_550_LKT(50,14)=206.630000 !rg=0.526233 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,15,1:6)=(/ 164.570000,168.930000,173.990000,184.010000,196.200000,196.150000 /)
XPIZA_LKT(50,15,1:6)=(/ 0.814727,0.921196,0.996814,0.996994,0.993281,0.551798 /)
XCGA_LKT(50,15,1:6)=(/ 0.853057,0.814413,0.774990,0.749570,0.729567,0.816750 /)
XEXT_COEFF_550_LKT(50,15)=171.090000 !rg=0.526233 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,16,1:6)=(/ 135.590000,138.200000,142.660000,149.470000,160.470000,161.290000 /)
XPIZA_LKT(50,16,1:6)=(/ 0.805768,0.916551,0.996412,0.996882,0.991912,0.543913 /)
XCGA_LKT(50,16,1:6)=(/ 0.859340,0.820630,0.784937,0.761663,0.736707,0.833283 /)
XEXT_COEFF_550_LKT(50,16)=139.800000 !rg=0.526233 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,17,1:6)=(/ 112.570000,114.230000,117.440000,122.850000,129.900000,132.490000 /)
XPIZA_LKT(50,17,1:6)=(/ 0.798029,0.911186,0.995962,0.996547,0.990770,0.538002 /)
XCGA_LKT(50,17,1:6)=(/ 0.863917,0.827373,0.793483,0.770987,0.747697,0.848380 /)
XEXT_COEFF_550_LKT(50,17)=115.770000 !rg=0.526233 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,18,1:6)=(/ 93.112000,94.858000,96.654000,100.610000,106.620000,108.710000 /)
XPIZA_LKT(50,18,1:6)=(/ 0.788735,0.905714,0.995793,0.995808,0.988282,0.533815 /)
XCGA_LKT(50,18,1:6)=(/ 0.867947,0.832660,0.797357,0.774813,0.751613,0.862060 /)
XEXT_COEFF_550_LKT(50,18)=95.953000 !rg=0.526233 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,19,1:6)=(/ 76.894000,77.878000,79.662000,82.652000,87.304000,88.997000 /)
XPIZA_LKT(50,19,1:6)=(/ 0.779648,0.901438,0.995117,0.995529,0.986537,0.531086 /)
XCGA_LKT(50,19,1:6)=(/ 0.871940,0.837120,0.803600,0.784950,0.761157,0.874427 /)
XEXT_COEFF_550_LKT(50,19)=78.598000 !rg=0.526233 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,20,1:6)=(/ 63.818000,64.448000,65.772000,67.879000,70.966000,72.905000 /)
XPIZA_LKT(50,20,1:6)=(/ 0.770286,0.895613,0.994979,0.995288,0.984435,0.529677 /)
XCGA_LKT(50,20,1:6)=(/ 0.875463,0.841433,0.810013,0.793440,0.772613,0.885210 /)
XEXT_COEFF_550_LKT(50,20)=65.141000 !rg=0.526233 sigma=2.95 wvl=0.55
 
END SUBROUTINE SALT_OPT_LKT_SET5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET6()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE

 
XEXT_COEFF_WVL_LKT(51,1,1:6)=(/ 1445.900000,1035.500000,2096.300000,1619.500000,524.410000,381.990000 /)
XPIZA_LKT(51,1,1:6)=(/ 0.944923,0.978048,0.999676,0.999719,0.997999,0.477284 /)
XCGA_LKT(51,1,1:6)=(/ 0.776920,0.553503,0.765090,0.754163,0.619797,0.247207 /)
XEXT_COEFF_550_LKT(51,1)=1498.900000 !rg=0.570098 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,2,1:6)=(/ 1279.000000,1074.300000,1885.200000,1606.200000,577.200000,411.240000 /)
XPIZA_LKT(51,2,1:6)=(/ 0.938619,0.978691,0.999641,0.999711,0.998182,0.500681 /)
XCGA_LKT(51,2,1:6)=(/ 0.765830,0.600297,0.750133,0.757897,0.622930,0.315377 /)
XEXT_COEFF_550_LKT(51,2)=1372.300000 !rg=0.570098 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,3,1:6)=(/ 1126.000000,1097.900000,1583.100000,1540.900000,644.000000,461.910000 /)
XPIZA_LKT(51,3,1:6)=(/ 0.931543,0.978692,0.999571,0.999695,0.998343,0.535837 /)
XCGA_LKT(51,3,1:6)=(/ 0.753570,0.649080,0.728363,0.759457,0.651003,0.418417 /)
XEXT_COEFF_550_LKT(51,3)=1222.500000 !rg=0.570098 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,4,1:6)=(/ 997.590000,1031.900000,1310.400000,1404.100000,710.820000,515.230000 /)
XPIZA_LKT(51,4,1:6)=(/ 0.926246,0.977243,0.999478,0.999660,0.998452,0.571847 /)
XCGA_LKT(51,4,1:6)=(/ 0.757107,0.683467,0.711303,0.754807,0.684693,0.507033 /)
XEXT_COEFF_550_LKT(51,4)=1098.800000 !rg=0.570098 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,5,1:6)=(/ 878.960000,914.220000,1094.400000,1219.500000,753.350000,555.940000 /)
XPIZA_LKT(51,5,1:6)=(/ 0.919176,0.975658,0.999362,0.999602,0.998490,0.600067 /)
XCGA_LKT(51,5,1:6)=(/ 0.769167,0.705623,0.702120,0.744790,0.710083,0.574890 /)
XEXT_COEFF_550_LKT(51,5)=971.660000 !rg=0.570098 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,6,1:6)=(/ 756.320000,802.250000,909.510000,1027.900000,755.090000,574.620000 /)
XPIZA_LKT(51,6,1:6)=(/ 0.912299,0.971370,0.999230,0.999519,0.998445,0.617793 /)
XCGA_LKT(51,6,1:6)=(/ 0.778213,0.723200,0.702820,0.732927,0.724343,0.625373 /)
XEXT_COEFF_550_LKT(51,6)=830.200000 !rg=0.570098 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,7,1:6)=(/ 647.290000,683.710000,755.110000,849.560000,714.300000,567.670000 /)
XPIZA_LKT(51,7,1:6)=(/ 0.902783,0.968287,0.999067,0.999415,0.998310,0.625513 /)
XCGA_LKT(51,7,1:6)=(/ 0.785637,0.737813,0.706790,0.723927,0.729057,0.663030 /)
XEXT_COEFF_550_LKT(51,7)=705.680000 !rg=0.570098 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,8,1:6)=(/ 547.510000,574.970000,625.310000,700.060000,644.200000,536.520000 /)
XPIZA_LKT(51,8,1:6)=(/ 0.892034,0.963860,0.998916,0.999262,0.998080,0.624339 /)
XCGA_LKT(51,8,1:6)=(/ 0.799267,0.751237,0.714860,0.717583,0.727533,0.691027 /)
XEXT_COEFF_550_LKT(51,8)=591.700000 !rg=0.570098 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,9,1:6)=(/ 462.520000,479.790000,517.970000,573.020000,562.970000,488.170000 /)
XPIZA_LKT(51,9,1:6)=(/ 0.877953,0.959075,0.998709,0.999041,0.997741,0.616861 /)
XCGA_LKT(51,9,1:6)=(/ 0.808877,0.760683,0.726060,0.717153,0.723943,0.713767 /)
XEXT_COEFF_550_LKT(51,9)=495.140000 !rg=0.570098 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,10,1:6)=(/ 386.430000,399.120000,423.150000,466.640000,478.350000,431.530000 /)
XPIZA_LKT(51,10,1:6)=(/ 0.864607,0.953825,0.998435,0.998899,0.997335,0.605913 /)
XCGA_LKT(51,10,1:6)=(/ 0.817177,0.775440,0.738183,0.717340,0.722503,0.733943 /)
XEXT_COEFF_550_LKT(51,10)=407.830000 !rg=0.570098 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,11,1:6)=(/ 321.470000,329.690000,350.360000,381.510000,398.330000,372.370000 /)
XPIZA_LKT(51,11,1:6)=(/ 0.851609,0.946437,0.997869,0.998579,0.996739,0.593208 /)
XCGA_LKT(51,11,1:6)=(/ 0.828463,0.786237,0.744167,0.723923,0.719413,0.753200 /)
XEXT_COEFF_550_LKT(51,11)=337.600000 !rg=0.570098 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,12,1:6)=(/ 267.010000,274.420000,284.950000,309.600000,328.800000,315.490000 /)
XPIZA_LKT(51,12,1:6)=(/ 0.839616,0.938198,0.997941,0.998386,0.996156,0.580050 /)
XCGA_LKT(51,12,1:6)=(/ 0.835773,0.795017,0.758990,0.727803,0.722900,0.771750 /)
XEXT_COEFF_550_LKT(51,12)=278.710000 !rg=0.570098 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,13,1:6)=(/ 221.050000,226.420000,235.650000,252.640000,268.850000,263.890000 /)
XPIZA_LKT(51,13,1:6)=(/ 0.828314,0.931573,0.997428,0.998062,0.995344,0.567713 /)
XCGA_LKT(51,13,1:6)=(/ 0.844270,0.802270,0.766153,0.739157,0.724310,0.790053 /)
XEXT_COEFF_550_LKT(51,13)=229.960000 !rg=0.570098 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,14,1:6)=(/ 183.000000,187.370000,194.150000,205.000000,219.400000,219.010000 /)
XPIZA_LKT(51,14,1:6)=(/ 0.818792,0.924891,0.996932,0.997778,0.994384,0.556902 /)
XCGA_LKT(51,14,1:6)=(/ 0.850247,0.811010,0.770723,0.746750,0.727723,0.807687 /)
XEXT_COEFF_550_LKT(51,14)=189.740000 !rg=0.570098 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,15,1:6)=(/ 151.300000,154.570000,160.440000,167.080000,180.200000,180.700000 /)
XPIZA_LKT(51,15,1:6)=(/ 0.809463,0.919090,0.996298,0.997344,0.992908,0.547922 /)
XCGA_LKT(51,15,1:6)=(/ 0.857303,0.817350,0.778160,0.757560,0.734147,0.824590 /)
XEXT_COEFF_550_LKT(51,15)=156.320000 !rg=0.570098 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,16,1:6)=(/ 125.130000,126.860000,131.280000,137.740000,146.010000,148.290000 /)
XPIZA_LKT(51,16,1:6)=(/ 0.801819,0.913556,0.995879,0.996329,0.991198,0.540744 /)
XCGA_LKT(51,16,1:6)=(/ 0.862147,0.825940,0.787833,0.764360,0.739970,0.840673 /)
XEXT_COEFF_550_LKT(51,16)=128.590000 !rg=0.570098 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,17,1:6)=(/ 103.600000,105.110000,107.970000,112.000000,118.810000,121.660000 /)
XPIZA_LKT(51,17,1:6)=(/ 0.793228,0.909032,0.995762,0.996310,0.990401,0.535556 /)
XCGA_LKT(51,17,1:6)=(/ 0.865683,0.830763,0.794257,0.772040,0.749660,0.855210 /)
XEXT_COEFF_550_LKT(51,17)=105.920000 !rg=0.570098 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,18,1:6)=(/ 85.768000,87.062000,89.305000,91.703000,97.737000,99.741000 /)
XPIZA_LKT(51,18,1:6)=(/ 0.784693,0.903815,0.995268,0.995805,0.987996,0.532060 /)
XCGA_LKT(51,18,1:6)=(/ 0.870370,0.835053,0.798937,0.782460,0.758127,0.868253 /)
XEXT_COEFF_550_LKT(51,18)=87.576000 !rg=0.570098 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,19,1:6)=(/ 71.009000,71.629000,73.523000,75.959000,79.493000,81.627000 /)
XPIZA_LKT(51,19,1:6)=(/ 0.775315,0.898341,0.995215,0.995142,0.985354,0.529955 /)
XCGA_LKT(51,19,1:6)=(/ 0.874227,0.840260,0.806437,0.787343,0.764270,0.879933 /)
XEXT_COEFF_550_LKT(51,19)=72.388000 !rg=0.570098 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,20,1:6)=(/ 58.741000,59.386000,60.480000,61.982000,65.038000,66.861000 /)
XPIZA_LKT(51,20,1:6)=(/ 0.764904,0.893271,0.994936,0.994909,0.983995,0.529083 /)
XCGA_LKT(51,20,1:6)=(/ 0.876797,0.843937,0.810587,0.794313,0.774830,0.890020 /)
XEXT_COEFF_550_LKT(51,20)=59.639000 !rg=0.570098 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,1,1:6)=(/ 1358.900000,896.310000,1818.700000,1633.400000,616.290000,422.100000 /)
XPIZA_LKT(52,1,1:6)=(/ 0.940417,0.974112,0.999630,0.999715,0.998218,0.510181 /)
XCGA_LKT(52,1,1:6)=(/ 0.795743,0.539953,0.747083,0.759790,0.641040,0.299013 /)
XEXT_COEFF_550_LKT(52,1)=1165.200000 !rg=0.617619 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,2,1:6)=(/ 1162.400000,978.030000,1618.700000,1611.800000,646.790000,452.870000 /)
XPIZA_LKT(52,2,1:6)=(/ 0.934846,0.976281,0.999579,0.999708,0.998334,0.528824 /)
XCGA_LKT(52,2,1:6)=(/ 0.782017,0.608747,0.729513,0.766257,0.644113,0.373237 /)
XEXT_COEFF_550_LKT(52,2)=1101.900000 !rg=0.617619 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,3,1:6)=(/ 1019.900000,992.810000,1364.300000,1509.200000,702.620000,502.060000 /)
XPIZA_LKT(52,3,1:6)=(/ 0.930161,0.976721,0.999495,0.999684,0.998441,0.559565 /)
XCGA_LKT(52,3,1:6)=(/ 0.775793,0.664067,0.707890,0.763847,0.672957,0.466693 /)
XEXT_COEFF_550_LKT(52,3)=1057.700000 !rg=0.617619 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,4,1:6)=(/ 905.300000,932.810000,1147.800000,1338.300000,754.710000,547.770000 /)
XPIZA_LKT(52,4,1:6)=(/ 0.924315,0.975278,0.999393,0.999639,0.998508,0.590106 /)
XCGA_LKT(52,4,1:6)=(/ 0.776300,0.693230,0.697143,0.754707,0.701680,0.542330 /)
XEXT_COEFF_550_LKT(52,4)=975.090000 !rg=0.617619 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,5,1:6)=(/ 802.690000,839.250000,965.500000,1137.300000,776.320000,577.050000 /)
XPIZA_LKT(52,5,1:6)=(/ 0.915008,0.972432,0.999276,0.999569,0.998506,0.612591 /)
XCGA_LKT(52,5,1:6)=(/ 0.774233,0.715250,0.696670,0.741360,0.721240,0.600930 /)
XEXT_COEFF_550_LKT(52,5)=867.790000 !rg=0.617619 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,6,1:6)=(/ 690.450000,724.490000,812.580000,945.650000,755.390000,582.510000 /)
XPIZA_LKT(52,6,1:6)=(/ 0.907953,0.969406,0.999167,0.999475,0.998419,0.625070 /)
XCGA_LKT(52,6,1:6)=(/ 0.785870,0.731303,0.704743,0.729513,0.730297,0.644807 /)
XEXT_COEFF_550_LKT(52,6)=746.050000 !rg=0.617619 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,7,1:6)=(/ 589.240000,619.930000,683.420000,776.490000,696.460000,562.720000 /)
XPIZA_LKT(52,7,1:6)=(/ 0.897960,0.966489,0.998895,0.999358,0.998239,0.628263 /)
XCGA_LKT(52,7,1:6)=(/ 0.796790,0.741660,0.709923,0.720900,0.731123,0.677743 /)
XEXT_COEFF_550_LKT(52,7)=634.990000 !rg=0.617619 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,8,1:6)=(/ 502.050000,522.140000,562.720000,636.670000,616.350000,521.720000 /)
XPIZA_LKT(52,8,1:6)=(/ 0.885122,0.961941,0.998862,0.999203,0.997969,0.623651 /)
XCGA_LKT(52,8,1:6)=(/ 0.803977,0.758740,0.723653,0.716183,0.728257,0.702783 /)
XEXT_COEFF_550_LKT(52,8)=534.660000 !rg=0.617619 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,9,1:6)=(/ 422.580000,441.140000,469.590000,518.930000,529.730000,467.440000 /)
XPIZA_LKT(52,9,1:6)=(/ 0.871115,0.955664,0.998524,0.998954,0.997589,0.613934 /)
XCGA_LKT(52,9,1:6)=(/ 0.812680,0.768810,0.725687,0.716930,0.724503,0.724007 /)
XEXT_COEFF_550_LKT(52,9)=450.010000 !rg=0.617619 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,10,1:6)=(/ 355.330000,365.610000,389.420000,427.340000,444.820000,408.270000 /)
XPIZA_LKT(52,10,1:6)=(/ 0.857286,0.950139,0.998280,0.998745,0.997109,0.601655 /)
XCGA_LKT(52,10,1:6)=(/ 0.824460,0.778653,0.738880,0.718913,0.720710,0.743363 /)
XEXT_COEFF_550_LKT(52,10)=375.280000 !rg=0.617619 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,11,1:6)=(/ 295.460000,303.860000,318.730000,346.730000,368.080000,348.870000 /)
XPIZA_LKT(52,11,1:6)=(/ 0.844703,0.942664,0.997837,0.998582,0.996513,0.588155 /)
XCGA_LKT(52,11,1:6)=(/ 0.831583,0.790767,0.749053,0.723173,0.721450,0.762037 /)
XEXT_COEFF_550_LKT(52,11)=309.020000 !rg=0.617619 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,12,1:6)=(/ 245.510000,251.640000,262.880000,284.020000,302.290000,293.520000 /)
XPIZA_LKT(52,12,1:6)=(/ 0.833058,0.935397,0.997600,0.997861,0.995648,0.574877 /)
XCGA_LKT(52,12,1:6)=(/ 0.841113,0.798127,0.760590,0.731673,0.721760,0.780353 /)
XEXT_COEFF_550_LKT(52,12)=255.570000 !rg=0.617619 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,13,1:6)=(/ 203.250000,208.100000,217.010000,229.320000,246.590000,244.240000 /)
XPIZA_LKT(52,13,1:6)=(/ 0.823056,0.928313,0.996698,0.997981,0.994973,0.562794 /)
XCGA_LKT(52,13,1:6)=(/ 0.847317,0.807823,0.763927,0.741603,0.725093,0.798443 /)
XEXT_COEFF_550_LKT(52,13)=210.590000 !rg=0.617619 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,14,1:6)=(/ 168.040000,172.360000,177.980000,186.180000,202.590000,201.970000 /)
XPIZA_LKT(52,14,1:6)=(/ 0.813950,0.921532,0.996865,0.997469,0.993587,0.552534 /)
XCGA_LKT(52,14,1:6)=(/ 0.854317,0.813147,0.778103,0.754927,0.729827,0.815817 /)
XEXT_COEFF_550_LKT(52,14)=174.340000 !rg=0.617619 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,15,1:6)=(/ 139.260000,141.890000,145.900000,153.640000,164.300000,166.240000 /)
XPIZA_LKT(52,15,1:6)=(/ 0.805840,0.916071,0.996331,0.997027,0.992567,0.544243 /)
XCGA_LKT(52,15,1:6)=(/ 0.860117,0.822027,0.782513,0.762093,0.734617,0.832327 /)
XEXT_COEFF_550_LKT(52,15)=143.850000 !rg=0.617619 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,16,1:6)=(/ 114.920000,117.040000,120.230000,125.730000,132.930000,136.190000 /)
XPIZA_LKT(52,16,1:6)=(/ 0.797147,0.911405,0.995943,0.996458,0.990907,0.537826 /)
XCGA_LKT(52,16,1:6)=(/ 0.864007,0.828323,0.790940,0.766370,0.746470,0.847893 /)
XEXT_COEFF_550_LKT(52,16)=118.350000 !rg=0.617619 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,17,1:6)=(/ 95.218000,96.684000,99.528000,102.160000,109.450000,111.630000 /)
XPIZA_LKT(52,17,1:6)=(/ 0.789176,0.906299,0.995327,0.995975,0.989096,0.533382 /)
XCGA_LKT(52,17,1:6)=(/ 0.868233,0.832430,0.797303,0.780123,0.753037,0.861817 /)
XEXT_COEFF_550_LKT(52,17)=97.854000 !rg=0.617619 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,18,1:6)=(/ 78.951000,80.047000,81.402000,84.451000,88.977000,91.477000 /)
XPIZA_LKT(52,18,1:6)=(/ 0.780295,0.901057,0.995408,0.995594,0.987396,0.530566 /)
XCGA_LKT(52,18,1:6)=(/ 0.872440,0.837700,0.802750,0.785877,0.759800,0.874183 /)
XEXT_COEFF_550_LKT(52,18)=80.930000 !rg=0.617619 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,19,1:6)=(/ 65.211000,66.250000,67.440000,69.540000,72.424000,74.850000 /)
XPIZA_LKT(52,19,1:6)=(/ 0.770310,0.896271,0.995050,0.995188,0.984887,0.529068 /)
XCGA_LKT(52,19,1:6)=(/ 0.875373,0.842280,0.807977,0.790720,0.771533,0.885150 /)
XEXT_COEFF_550_LKT(52,19)=66.605000 !rg=0.617619 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,20,1:6)=(/ 54.049000,54.634000,55.809000,56.808000,59.815000,61.314000 /)
XPIZA_LKT(52,20,1:6)=(/ 0.759925,0.890019,0.994789,0.994806,0.982752,0.528705 /)
XCGA_LKT(52,20,1:6)=(/ 0.878627,0.845000,0.813807,0.800150,0.778367,0.894540 /)
XEXT_COEFF_550_LKT(52,20)=55.143000 !rg=0.617619 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,1,1:6)=(/ 1144.100000,871.970000,1523.000000,1664.600000,692.590000,464.660000 /)
XPIZA_LKT(53,1,1:6)=(/ 0.932093,0.972751,0.999548,0.999713,0.998413,0.537374 /)
XCGA_LKT(53,1,1:6)=(/ 0.784187,0.574117,0.717463,0.773463,0.643613,0.362410 /)
XEXT_COEFF_550_LKT(53,1)=907.260000 !rg=0.669101 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,2,1:6)=(/ 1017.900000,930.990000,1350.400000,1587.300000,709.070000,496.910000 /)
XPIZA_LKT(53,2,1:6)=(/ 0.930097,0.975231,0.999491,0.999700,0.998453,0.553280 /)
XCGA_LKT(53,2,1:6)=(/ 0.786003,0.645160,0.701113,0.771830,0.663023,0.434903 /)
XEXT_COEFF_550_LKT(53,2)=928.720000 !rg=0.669101 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,3,1:6)=(/ 934.820000,923.130000,1153.300000,1449.300000,755.720000,541.110000 /)
XPIZA_LKT(53,3,1:6)=(/ 0.924524,0.975598,0.999403,0.999667,0.998517,0.580739 /)
XCGA_LKT(53,3,1:6)=(/ 0.785007,0.689763,0.684913,0.764913,0.692857,0.511850 /)
XEXT_COEFF_550_LKT(53,3)=929.340000 !rg=0.669101 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,4,1:6)=(/ 835.330000,849.750000,995.490000,1253.800000,790.590000,576.590000 /)
XPIZA_LKT(53,4,1:6)=(/ 0.918458,0.974286,0.999303,0.999610,0.998546,0.606009 /)
XCGA_LKT(53,4,1:6)=(/ 0.785220,0.713887,0.684393,0.751343,0.716417,0.574770 /)
XEXT_COEFF_550_LKT(53,4)=872.590000 !rg=0.669101 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,5,1:6)=(/ 728.800000,762.750000,861.390000,1046.900000,789.430000,592.930000 /)
XPIZA_LKT(53,5,1:6)=(/ 0.911689,0.970727,0.999159,0.999528,0.998504,0.622965 /)
XCGA_LKT(53,5,1:6)=(/ 0.789780,0.721370,0.691813,0.736707,0.730313,0.624877 /)
XEXT_COEFF_550_LKT(53,5)=781.410000 !rg=0.669101 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,6,1:6)=(/ 637.770000,657.030000,728.730000,865.000000,746.100000,584.720000 /)
XPIZA_LKT(53,6,1:6)=(/ 0.901200,0.968448,0.999041,0.999413,0.998376,0.630451 /)
XCGA_LKT(53,6,1:6)=(/ 0.795607,0.741730,0.702637,0.723577,0.734553,0.662687 /)
XEXT_COEFF_550_LKT(53,6)=682.700000 !rg=0.669101 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,7,1:6)=(/ 543.880000,560.700000,612.780000,708.620000,671.290000,552.860000 /)
XPIZA_LKT(53,7,1:6)=(/ 0.890615,0.965037,0.998762,0.999279,0.998154,0.629450 /)
XCGA_LKT(53,7,1:6)=(/ 0.803640,0.756590,0.712857,0.716197,0.732413,0.691433 /)
XEXT_COEFF_550_LKT(53,7)=577.910000 !rg=0.669101 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,8,1:6)=(/ 462.650000,476.480000,512.990000,580.910000,584.280000,503.530000 /)
XPIZA_LKT(53,8,1:6)=(/ 0.876986,0.959625,0.998718,0.999086,0.997831,0.621813 /)
XCGA_LKT(53,8,1:6)=(/ 0.812203,0.763737,0.724630,0.714290,0.727593,0.714013 /)
XEXT_COEFF_550_LKT(53,8)=491.220000 !rg=0.669101 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,9,1:6)=(/ 386.260000,401.760000,429.650000,471.220000,496.020000,444.820000 /)
XPIZA_LKT(53,9,1:6)=(/ 0.864077,0.953191,0.998359,0.998931,0.997367,0.610195 /)
XCGA_LKT(53,9,1:6)=(/ 0.821683,0.771023,0.731433,0.717607,0.723017,0.733987 /)
XEXT_COEFF_550_LKT(53,9)=409.760000 !rg=0.669101 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,10,1:6)=(/ 325.120000,336.960000,354.420000,386.020000,411.630000,384.220000 /)
XPIZA_LKT(53,10,1:6)=(/ 0.851069,0.945264,0.998009,0.998595,0.996871,0.596774 /)
XCGA_LKT(53,10,1:6)=(/ 0.827917,0.785197,0.740107,0.720983,0.720273,0.752580 /)
XEXT_COEFF_550_LKT(53,10)=342.530000 !rg=0.669101 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,11,1:6)=(/ 271.610000,279.490000,292.430000,315.640000,340.560000,325.540000 /)
XPIZA_LKT(53,11,1:6)=(/ 0.838884,0.937611,0.997730,0.998382,0.996084,0.582767 /)
XCGA_LKT(53,11,1:6)=(/ 0.837873,0.792400,0.755870,0.729637,0.719617,0.770823 /)
XEXT_COEFF_550_LKT(53,11)=284.970000 !rg=0.669101 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,12,1:6)=(/ 225.270000,231.450000,240.980000,256.560000,277.010000,272.240000 /)
XPIZA_LKT(53,12,1:6)=(/ 0.827918,0.931135,0.997180,0.998136,0.995426,0.569591 /)
XCGA_LKT(53,12,1:6)=(/ 0.844470,0.803533,0.760160,0.736703,0.721860,0.788963 /)
XEXT_COEFF_550_LKT(53,12)=234.670000 !rg=0.669101 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,13,1:6)=(/ 186.440000,191.210000,198.430000,207.850000,227.110000,225.520000 /)
XPIZA_LKT(53,13,1:6)=(/ 0.818453,0.924572,0.996921,0.997895,0.994269,0.557919 /)
XCGA_LKT(53,13,1:6)=(/ 0.851580,0.809580,0.773170,0.750837,0.726430,0.806820 /)
XEXT_COEFF_550_LKT(53,13)=194.140000 !rg=0.669101 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,14,1:6)=(/ 154.760000,157.760000,162.410000,171.810000,184.620000,185.950000 /)
XPIZA_LKT(53,14,1:6)=(/ 0.809621,0.919099,0.996393,0.997330,0.993305,0.548327 /)
XCGA_LKT(53,14,1:6)=(/ 0.857893,0.818483,0.779087,0.755883,0.729657,0.823873 /)
XEXT_COEFF_550_LKT(53,14)=159.800000 !rg=0.669101 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,15,1:6)=(/ 128.170000,130.450000,133.420000,140.550000,149.320000,152.750000 /)
XPIZA_LKT(53,15,1:6)=(/ 0.801320,0.913932,0.995904,0.996619,0.992233,0.540806 /)
XCGA_LKT(53,15,1:6)=(/ 0.861940,0.826193,0.788907,0.761953,0.743417,0.839930 /)
XEXT_COEFF_550_LKT(53,15)=131.930000 !rg=0.669101 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,16,1:6)=(/ 105.970000,107.800000,110.370000,114.910000,122.620000,124.990000 /)
XPIZA_LKT(53,16,1:6)=(/ 0.793459,0.908437,0.995851,0.996245,0.990098,0.535185 /)
XCGA_LKT(53,16,1:6)=(/ 0.866657,0.830087,0.795493,0.772903,0.747837,0.854917 /)
XEXT_COEFF_550_LKT(53,16)=109.130000 !rg=0.669101 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,17,1:6)=(/ 87.722000,89.072000,91.039000,94.208000,99.904000,102.380000 /)
XPIZA_LKT(53,17,1:6)=(/ 0.784543,0.903880,0.995368,0.995838,0.988467,0.531486 /)
XCGA_LKT(53,17,1:6)=(/ 0.870443,0.835350,0.799760,0.782647,0.753793,0.868173 /)
XEXT_COEFF_550_LKT(53,17)=89.453000 !rg=0.669101 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,18,1:6)=(/ 72.697000,73.721000,74.609000,77.339000,81.031000,83.872000 /)
XPIZA_LKT(53,18,1:6)=(/ 0.775615,0.899029,0.995308,0.995514,0.986438,0.529340 /)
XCGA_LKT(53,18,1:6)=(/ 0.873697,0.840697,0.807830,0.787537,0.768387,0.879827 /)
XEXT_COEFF_550_LKT(53,18)=74.181000 !rg=0.669101 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,19,1:6)=(/ 60.177000,60.809000,61.934000,63.765000,66.733000,68.627000 /)
XPIZA_LKT(53,19,1:6)=(/ 0.766095,0.893293,0.994822,0.994736,0.984013,0.528422 /)
XCGA_LKT(53,19,1:6)=(/ 0.877163,0.843503,0.811180,0.794310,0.774363,0.890070 /)
XEXT_COEFF_550_LKT(53,19)=61.446000 !rg=0.669101 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,20,1:6)=(/ 49.860000,50.418000,51.161000,52.430000,54.857000,56.226000 /)
XPIZA_LKT(53,20,1:6)=(/ 0.754736,0.887401,0.994509,0.994696,0.981867,0.528527 /)
XCGA_LKT(53,20,1:6)=(/ 0.880477,0.847300,0.815460,0.802027,0.779663,0.898760 /)
XEXT_COEFF_550_LKT(53,20)=50.558000 !rg=0.669101 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,1,1:6)=(/ 932.320000,923.090000,1218.400000,1611.600000,737.510000,510.240000 /)
XPIZA_LKT(54,1,1:6)=(/ 0.907758,0.974112,0.999434,0.999706,0.998538,0.559206 /)
XCGA_LKT(54,1,1:6)=(/ 0.768927,0.652100,0.678250,0.780897,0.651727,0.435953 /)
XEXT_COEFF_550_LKT(54,1)=770.530000 !rg=0.724875 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,2,1:6)=(/ 911.090000,924.170000,1098.700000,1530.100000,764.400000,541.980000 /)
XPIZA_LKT(54,2,1:6)=(/ 0.920816,0.975124,0.999370,0.999685,0.998539,0.575177 /)
XCGA_LKT(54,2,1:6)=(/ 0.791623,0.692380,0.662903,0.773940,0.684457,0.494640 /)
XEXT_COEFF_550_LKT(54,2)=829.220000 !rg=0.724875 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,3,1:6)=(/ 848.810000,867.200000,978.230000,1362.300000,801.550000,577.080000 /)
XPIZA_LKT(54,3,1:6)=(/ 0.917960,0.974043,0.999293,0.999642,0.998573,0.599616 /)
XCGA_LKT(54,3,1:6)=(/ 0.783697,0.712083,0.665123,0.762463,0.710997,0.552410 /)
XEXT_COEFF_550_LKT(54,3)=846.180000 !rg=0.724875 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,4,1:6)=(/ 764.670000,784.060000,871.960000,1152.300000,816.730000,600.380000 /)
XPIZA_LKT(54,4,1:6)=(/ 0.911523,0.972170,0.999207,0.999573,0.998566,0.619617 /)
XCGA_LKT(54,4,1:6)=(/ 0.785893,0.726957,0.677313,0.745787,0.728923,0.604257 /)
XEXT_COEFF_550_LKT(54,4)=793.270000 !rg=0.724875 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,5,1:6)=(/ 671.890000,690.870000,762.460000,953.300000,791.720000,602.830000 /)
XPIZA_LKT(54,5,1:6)=(/ 0.905365,0.969793,0.999065,0.999473,0.998485,0.631236 /)
XCGA_LKT(54,5,1:6)=(/ 0.796967,0.742293,0.690433,0.728780,0.737380,0.646783 /)
XEXT_COEFF_550_LKT(54,5)=707.090000 !rg=0.724875 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,6,1:6)=(/ 584.780000,606.260000,654.770000,780.820000,728.030000,581.150000 /)
XPIZA_LKT(54,6,1:6)=(/ 0.894306,0.965555,0.998965,0.999352,0.998314,0.633988 /)
XCGA_LKT(54,6,1:6)=(/ 0.796910,0.750457,0.703677,0.718477,0.737530,0.679137 /)
XEXT_COEFF_550_LKT(54,6)=619.030000 !rg=0.724875 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,7,1:6)=(/ 499.800000,515.550000,553.300000,639.630000,641.410000,538.440000 /)
XPIZA_LKT(54,7,1:6)=(/ 0.882501,0.961952,0.998753,0.999199,0.998046,0.629155 /)
XCGA_LKT(54,7,1:6)=(/ 0.806450,0.762497,0.716763,0.712980,0.732967,0.704237 /)
XEXT_COEFF_550_LKT(54,7)=525.980000 !rg=0.724875 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,8,1:6)=(/ 423.460000,439.480000,465.800000,523.330000,549.240000,482.450000 /)
XPIZA_LKT(54,8,1:6)=(/ 0.870229,0.955565,0.998509,0.999036,0.997674,0.618897 /)
XCGA_LKT(54,8,1:6)=(/ 0.815257,0.771710,0.725090,0.713477,0.727277,0.724793 /)
XEXT_COEFF_550_LKT(54,8)=448.550000 !rg=0.724875 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,9,1:6)=(/ 355.940000,365.760000,388.570000,430.000000,459.220000,420.750000 /)
XPIZA_LKT(54,9,1:6)=(/ 0.857216,0.949756,0.998148,0.998816,0.997161,0.605659 /)
XCGA_LKT(54,9,1:6)=(/ 0.827410,0.782767,0.736337,0.717973,0.721767,0.743700 /)
XEXT_COEFF_550_LKT(54,9)=372.740000 !rg=0.724875 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,10,1:6)=(/ 298.140000,307.780000,324.390000,349.510000,381.800000,359.830000 /)
XPIZA_LKT(54,10,1:6)=(/ 0.844226,0.942457,0.997899,0.998579,0.996465,0.591343 /)
XCGA_LKT(54,10,1:6)=(/ 0.834987,0.787423,0.748097,0.726260,0.719317,0.761643 /)
XEXT_COEFF_550_LKT(54,10)=313.560000 !rg=0.724875 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,11,1:6)=(/ 248.470000,256.270000,267.300000,285.720000,311.650000,302.670000 /)
XPIZA_LKT(54,11,1:6)=(/ 0.832998,0.935043,0.997608,0.998329,0.995804,0.577147 /)
XCGA_LKT(54,11,1:6)=(/ 0.841243,0.800113,0.754920,0.734630,0.719090,0.779620 /)
XEXT_COEFF_550_LKT(54,11)=259.620000 !rg=0.724875 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,12,1:6)=(/ 206.740000,212.360000,220.400000,232.800000,255.610000,251.790000 /)
XPIZA_LKT(54,12,1:6)=(/ 0.822350,0.927872,0.997293,0.997877,0.994795,0.564244 /)
XCGA_LKT(54,12,1:6)=(/ 0.849253,0.805910,0.769343,0.744100,0.722843,0.797563 /)
XEXT_COEFF_550_LKT(54,12)=215.270000 !rg=0.724875 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,13,1:6)=(/ 171.780000,175.270000,180.730000,191.540000,207.280000,207.830000 /)
XPIZA_LKT(54,13,1:6)=(/ 0.813149,0.921780,0.996872,0.997256,0.993924,0.553161 /)
XCGA_LKT(54,13,1:6)=(/ 0.855457,0.815187,0.775347,0.751717,0.725577,0.815167 /)
XEXT_COEFF_550_LKT(54,13)=177.340000 !rg=0.724875 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,14,1:6)=(/ 142.380000,145.300000,148.160000,157.100000,168.190000,170.950000 /)
XPIZA_LKT(54,14,1:6)=(/ 0.805727,0.915557,0.996474,0.996926,0.992616,0.544333 /)
XCGA_LKT(54,14,1:6)=(/ 0.859613,0.822200,0.787470,0.756740,0.737567,0.831827 /)
XEXT_COEFF_550_LKT(54,14)=145.890000 !rg=0.724875 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,15,1:6)=(/ 118.100000,120.190000,123.300000,129.140000,137.010000,140.220000 /)
XPIZA_LKT(54,15,1:6)=(/ 0.797291,0.911036,0.995941,0.996376,0.991491,0.537638 /)
XCGA_LKT(54,15,1:6)=(/ 0.865017,0.827603,0.790647,0.766760,0.745353,0.847363 /)
XEXT_COEFF_550_LKT(54,15)=121.020000 !rg=0.724875 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,16,1:6)=(/ 97.309000,99.199000,101.650000,105.240000,111.990000,114.630000 /)
XPIZA_LKT(54,16,1:6)=(/ 0.789584,0.906262,0.995681,0.995913,0.989461,0.532834 /)
XCGA_LKT(54,16,1:6)=(/ 0.868423,0.833800,0.795430,0.778217,0.749207,0.861710 /)
XEXT_COEFF_550_LKT(54,16)=99.856000 !rg=0.724875 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,17,1:6)=(/ 80.660000,81.886000,83.498000,86.875000,91.320000,93.856000 /)
XPIZA_LKT(54,17,1:6)=(/ 0.780584,0.901137,0.995440,0.995169,0.987213,0.529881 /)
XCGA_LKT(54,17,1:6)=(/ 0.872133,0.836870,0.804537,0.780663,0.761047,0.874253 /)
XEXT_COEFF_550_LKT(54,17)=82.018000 !rg=0.724875 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,18,1:6)=(/ 67.078000,67.921000,69.141000,71.338000,74.611000,76.884000 /)
XPIZA_LKT(54,18,1:6)=(/ 0.771375,0.896292,0.994822,0.995203,0.985236,0.528380 /)
XCGA_LKT(54,18,1:6)=(/ 0.875920,0.841253,0.808880,0.790113,0.771463,0.885173 /)
XEXT_COEFF_550_LKT(54,18)=68.294000 !rg=0.724875 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,19,1:6)=(/ 55.347000,56.075000,57.189000,58.508000,61.058000,62.916000 /)
XPIZA_LKT(54,19,1:6)=(/ 0.761098,0.890605,0.994633,0.994709,0.983344,0.528003 /)
XCGA_LKT(54,19,1:6)=(/ 0.878563,0.845763,0.810780,0.799683,0.776320,0.894683 /)
XEXT_COEFF_550_LKT(54,19)=56.336000 !rg=0.724875 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,20,1:6)=(/ 45.867000,46.394000,47.104000,48.377000,50.097000,51.561000 /)
XPIZA_LKT(54,20,1:6)=(/ 0.749538,0.884116,0.994512,0.994225,0.980974,0.528543 /)
XCGA_LKT(54,20,1:6)=(/ 0.881737,0.848347,0.817993,0.801537,0.785667,0.902683 /)
XEXT_COEFF_550_LKT(54,20)=46.407000 !rg=0.724875 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,1,1:6)=(/ 784.070000,991.520000,939.700000,1547.200000,774.780000,559.720000 /)
XPIZA_LKT(55,1,1:6)=(/ 0.901745,0.975831,0.999270,0.999683,0.998582,0.578110 /)
XCGA_LKT(55,1,1:6)=(/ 0.742550,0.727443,0.624113,0.778897,0.681240,0.511457 /)
XEXT_COEFF_550_LKT(55,1)=712.010000 !rg=0.785298 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,2,1:6)=(/ 828.240000,914.860000,887.330000,1440.100000,814.710000,585.160000 /)
XPIZA_LKT(55,2,1:6)=(/ 0.911407,0.974403,0.999223,0.999661,0.998599,0.595655 /)
XCGA_LKT(55,2,1:6)=(/ 0.774963,0.729610,0.628200,0.771823,0.707867,0.546267 /)
XEXT_COEFF_550_LKT(55,2)=787.520000 !rg=0.785298 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,3,1:6)=(/ 783.640000,811.320000,842.520000,1252.300000,838.370000,608.010000 /)
XPIZA_LKT(55,3,1:6)=(/ 0.910522,0.972547,0.999161,0.999607,0.998610,0.616279 /)
XCGA_LKT(55,3,1:6)=(/ 0.795757,0.729850,0.644587,0.756013,0.727113,0.587930 /)
XEXT_COEFF_550_LKT(55,3)=793.610000 !rg=0.785298 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,4,1:6)=(/ 704.420000,721.520000,775.540000,1043.400000,831.640000,618.020000 /)
XPIZA_LKT(55,4,1:6)=(/ 0.905326,0.969874,0.999086,0.999523,0.998568,0.630937 /)
XCGA_LKT(55,4,1:6)=(/ 0.798163,0.736443,0.670357,0.736867,0.739143,0.630903 /)
XEXT_COEFF_550_LKT(55,4)=733.650000 !rg=0.785298 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,5,1:6)=(/ 617.090000,634.910000,682.090000,855.420000,783.070000,606.280000 /)
XPIZA_LKT(55,5,1:6)=(/ 0.897624,0.967264,0.998997,0.999414,0.998446,0.637403 /)
XCGA_LKT(55,5,1:6)=(/ 0.798633,0.751390,0.694267,0.721860,0.742437,0.666750 /)
XEXT_COEFF_550_LKT(55,5)=644.890000 !rg=0.785298 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,6,1:6)=(/ 536.120000,554.580000,596.960000,703.200000,702.330000,571.890000 /)
XPIZA_LKT(55,6,1:6)=(/ 0.887070,0.962199,0.998781,0.999290,0.998225,0.635714 /)
XCGA_LKT(55,6,1:6)=(/ 0.807613,0.751963,0.706477,0.714697,0.738417,0.694257 /)
XEXT_COEFF_550_LKT(55,6)=568.730000 !rg=0.785298 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,7,1:6)=(/ 458.550000,473.180000,505.400000,578.060000,607.100000,519.960000 /)
XPIZA_LKT(55,7,1:6)=(/ 0.875860,0.957355,0.998628,0.999100,0.997904,0.627441 /)
XCGA_LKT(55,7,1:6)=(/ 0.815907,0.764493,0.722023,0.712243,0.731523,0.716270 /)
XEXT_COEFF_550_LKT(55,7)=485.750000 !rg=0.785298 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,8,1:6)=(/ 388.050000,401.710000,426.060000,472.380000,512.980000,459.030000 /)
XPIZA_LKT(55,8,1:6)=(/ 0.862514,0.952623,0.998403,0.998963,0.997466,0.614947 /)
XCGA_LKT(55,8,1:6)=(/ 0.823653,0.772737,0.733473,0.716373,0.724953,0.735170 /)
XEXT_COEFF_550_LKT(55,8)=411.340000 !rg=0.785298 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,9,1:6)=(/ 327.570000,335.850000,351.400000,389.240000,423.760000,395.680000 /)
XPIZA_LKT(55,9,1:6)=(/ 0.849472,0.946206,0.998179,0.998631,0.996990,0.600365 /)
XCGA_LKT(55,9,1:6)=(/ 0.830550,0.789043,0.746293,0.716157,0.722273,0.753170 /)
XEXT_COEFF_550_LKT(55,9)=341.130000 !rg=0.785298 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,10,1:6)=(/ 274.240000,281.680000,293.960000,320.320000,349.940000,335.490000 /)
XPIZA_LKT(55,10,1:6)=(/ 0.838105,0.938151,0.997582,0.998425,0.996231,0.585473 /)
XCGA_LKT(55,10,1:6)=(/ 0.840247,0.796010,0.752990,0.727083,0.718353,0.770637 /)
XEXT_COEFF_550_LKT(55,10)=286.400000 !rg=0.785298 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,11,1:6)=(/ 227.990000,234.010000,245.020000,260.280000,286.810000,280.470000 /)
XPIZA_LKT(55,11,1:6)=(/ 0.826759,0.932056,0.997275,0.998211,0.995174,0.571369 /)
XCGA_LKT(55,11,1:6)=(/ 0.847160,0.802667,0.762003,0.736627,0.720027,0.788430 /)
XEXT_COEFF_550_LKT(55,11)=237.000000 !rg=0.785298 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,12,1:6)=(/ 190.260000,194.360000,200.760000,214.540000,232.900000,232.300000 /)
XPIZA_LKT(55,12,1:6)=(/ 0.817542,0.924783,0.997046,0.997726,0.994502,0.558916 /)
XCGA_LKT(55,12,1:6)=(/ 0.853337,0.811907,0.771090,0.744643,0.722193,0.806157 /)
XEXT_COEFF_550_LKT(55,12)=197.390000 !rg=0.785298 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,13,1:6)=(/ 157.860000,161.100000,164.810000,175.210000,189.070000,191.180000 /)
XPIZA_LKT(55,13,1:6)=(/ 0.809469,0.918135,0.996674,0.997206,0.993280,0.548577 /)
XCGA_LKT(55,13,1:6)=(/ 0.857687,0.818517,0.783610,0.752043,0.732743,0.823450 /)
XEXT_COEFF_550_LKT(55,13)=161.990000 !rg=0.785298 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,14,1:6)=(/ 130.870000,133.280000,136.780000,143.930000,153.390000,156.960000 /)
XPIZA_LKT(55,14,1:6)=(/ 0.801120,0.913815,0.996282,0.996854,0.992284,0.540597 /)
XCGA_LKT(55,14,1:6)=(/ 0.862850,0.824777,0.790197,0.763853,0.741133,0.839653 /)
XEXT_COEFF_550_LKT(55,14)=134.350000 !rg=0.785298 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,15,1:6)=(/ 108.440000,110.520000,112.840000,117.660000,125.770000,128.610000 /)
XPIZA_LKT(55,15,1:6)=(/ 0.793230,0.907998,0.995845,0.996220,0.989759,0.534768 /)
XCGA_LKT(55,15,1:6)=(/ 0.867003,0.830630,0.793560,0.769467,0.744593,0.854593 /)
XEXT_COEFF_550_LKT(55,15)=111.840000 !rg=0.785298 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,16,1:6)=(/ 89.534000,90.806000,93.049000,96.489000,102.950000,105.080000 /)
XPIZA_LKT(55,16,1:6)=(/ 0.784711,0.903956,0.995535,0.995880,0.987597,0.530794 /)
XCGA_LKT(55,16,1:6)=(/ 0.870953,0.835447,0.801020,0.780433,0.753983,0.868240 /)
XEXT_COEFF_550_LKT(55,16)=91.612000 !rg=0.785298 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,17,1:6)=(/ 74.384000,75.223000,76.773000,79.485000,83.433000,86.019000 /)
XPIZA_LKT(55,17,1:6)=(/ 0.776907,0.899171,0.995142,0.995319,0.986249,0.528565 /)
XCGA_LKT(55,17,1:6)=(/ 0.874220,0.840187,0.807627,0.788700,0.765777,0.880040 /)
XEXT_COEFF_550_LKT(55,17)=75.883000 !rg=0.785298 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,18,1:6)=(/ 61.597000,62.452000,63.255000,65.262000,68.566000,70.470000 /)
XPIZA_LKT(55,18,1:6)=(/ 0.766511,0.893428,0.994937,0.994936,0.983401,0.527676 /)
XCGA_LKT(55,18,1:6)=(/ 0.877320,0.843710,0.811127,0.792083,0.771330,0.890207 /)
XEXT_COEFF_550_LKT(55,18)=63.004000 !rg=0.785298 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,19,1:6)=(/ 50.971000,51.497000,52.311000,53.886000,56.279000,57.681000 /)
XPIZA_LKT(55,19,1:6)=(/ 0.755781,0.887888,0.994756,0.994653,0.981756,0.527807 /)
XCGA_LKT(55,19,1:6)=(/ 0.880300,0.847163,0.816057,0.799863,0.780067,0.898987 /)
XEXT_COEFF_550_LKT(55,19)=51.915000 !rg=0.785298 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,20,1:6)=(/ 42.348000,42.632000,43.333000,44.372000,46.023000,47.286000 /)
XPIZA_LKT(55,20,1:6)=(/ 0.744632,0.881033,0.994345,0.994392,0.980206,0.528736 /)
XCGA_LKT(55,20,1:6)=(/ 0.883620,0.850600,0.820290,0.807340,0.789973,0.906313 /)
XEXT_COEFF_550_LKT(55,20)=42.999000 !rg=0.785298 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,1,1:6)=(/ 760.010000,992.830000,744.230000,1445.200000,837.590000,610.670000 /)
XPIZA_LKT(56,1,1:6)=(/ 0.914722,0.975923,0.999045,0.999664,0.998606,0.598016 /)
XCGA_LKT(56,1,1:6)=(/ 0.800940,0.771713,0.560803,0.775000,0.722347,0.573463 /)
XEXT_COEFF_550_LKT(56,1)=745.800000 !rg=0.850757 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,2,1:6)=(/ 772.020000,861.790000,745.610000,1319.200000,858.730000,622.350000 /)
XPIZA_LKT(56,2,1:6)=(/ 0.907892,0.973394,0.999049,0.999627,0.998641,0.615108 /)
XCGA_LKT(56,2,1:6)=(/ 0.797353,0.754943,0.593503,0.764800,0.729243,0.586530 /)
XEXT_COEFF_550_LKT(56,2)=784.140000 !rg=0.850757 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,3,1:6)=(/ 714.860000,758.820000,734.690000,1124.300000,864.280000,632.260000 /)
XPIZA_LKT(56,3,1:6)=(/ 0.906834,0.968962,0.999049,0.999559,0.998628,0.630626 /)
XCGA_LKT(56,3,1:6)=(/ 0.793217,0.748763,0.644160,0.745197,0.740750,0.618870 /)
XEXT_COEFF_550_LKT(56,3)=741.220000 !rg=0.850757 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,4,1:6)=(/ 642.870000,672.050000,690.370000,928.000000,834.090000,628.650000 /)
XPIZA_LKT(56,4,1:6)=(/ 0.900642,0.965959,0.998976,0.999461,0.998550,0.639932 /)
XCGA_LKT(56,4,1:6)=(/ 0.799237,0.751347,0.675287,0.725247,0.746997,0.654920 /)
XEXT_COEFF_550_LKT(56,4)=674.050000 !rg=0.850757 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,5,1:6)=(/ 567.130000,581.640000,620.300000,765.510000,763.550000,603.030000 /)
XPIZA_LKT(56,5,1:6)=(/ 0.891603,0.964470,0.998861,0.999333,0.998383,0.641466 /)
XCGA_LKT(56,5,1:6)=(/ 0.809047,0.755920,0.696593,0.713167,0.745270,0.684893 /)
XEXT_COEFF_550_LKT(56,5)=596.530000 !rg=0.850757 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,6,1:6)=(/ 490.720000,510.470000,538.620000,630.960000,668.550000,557.270000 /)
XPIZA_LKT(56,6,1:6)=(/ 0.880123,0.958624,0.998709,0.999201,0.998123,0.635665 /)
XCGA_LKT(56,6,1:6)=(/ 0.813597,0.764167,0.708520,0.707980,0.738670,0.708183 /)
XEXT_COEFF_550_LKT(56,6)=515.250000 !rg=0.850757 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,7,1:6)=(/ 419.270000,435.040000,458.900000,518.620000,568.120000,498.000000 /)
XPIZA_LKT(56,7,1:6)=(/ 0.868657,0.954742,0.998427,0.999032,0.997752,0.624367 /)
XCGA_LKT(56,7,1:6)=(/ 0.819833,0.775347,0.721453,0.710293,0.730703,0.727627 /)
XEXT_COEFF_550_LKT(56,7)=442.250000 !rg=0.850757 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,8,1:6)=(/ 356.330000,367.900000,386.340000,429.640000,473.600000,433.820000 /)
XPIZA_LKT(56,8,1:6)=(/ 0.855425,0.948434,0.998320,0.998813,0.997242,0.610000 /)
XCGA_LKT(56,8,1:6)=(/ 0.829860,0.783113,0.735543,0.714157,0.723183,0.745183 /)
XEXT_COEFF_550_LKT(56,8)=374.390000 !rg=0.850757 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,9,1:6)=(/ 300.850000,308.640000,324.130000,354.640000,390.520000,370.130000 /)
XPIZA_LKT(56,9,1:6)=(/ 0.842858,0.942150,0.998002,0.998480,0.996674,0.594393 /)
XCGA_LKT(56,9,1:6)=(/ 0.837233,0.791457,0.750143,0.720917,0.719607,0.762453 /)
XEXT_COEFF_550_LKT(56,9)=315.070000 !rg=0.850757 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,10,1:6)=(/ 252.370000,258.200000,266.520000,290.500000,319.490000,311.560000 /)
XPIZA_LKT(56,10,1:6)=(/ 0.831848,0.935148,0.997782,0.998294,0.996042,0.579290 /)
XCGA_LKT(56,10,1:6)=(/ 0.843287,0.802400,0.763723,0.727363,0.721777,0.779627 /)
XEXT_COEFF_550_LKT(56,10)=262.030000 !rg=0.850757 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,11,1:6)=(/ 210.240000,214.060000,223.440000,239.250000,260.860000,259.130000 /)
XPIZA_LKT(56,11,1:6)=(/ 0.821873,0.928177,0.997018,0.997863,0.995081,0.565497 /)
XCGA_LKT(56,11,1:6)=(/ 0.850880,0.810327,0.766087,0.740050,0.721593,0.797247 /)
XEXT_COEFF_550_LKT(56,11)=217.440000 !rg=0.850757 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,12,1:6)=(/ 175.140000,178.820000,182.760000,195.620000,211.760000,213.860000 /)
XPIZA_LKT(56,12,1:6)=(/ 0.813261,0.921222,0.996774,0.997426,0.994151,0.553703 /)
XCGA_LKT(56,12,1:6)=(/ 0.855440,0.816517,0.780690,0.746073,0.729763,0.814730 /)
XEXT_COEFF_550_LKT(56,12)=179.950000 !rg=0.850757 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,13,1:6)=(/ 145.040000,147.700000,152.070000,160.300000,171.940000,175.600000 /)
XPIZA_LKT(56,13,1:6)=(/ 0.805011,0.916054,0.996353,0.997227,0.993002,0.544220 /)
XCGA_LKT(56,13,1:6)=(/ 0.861043,0.822320,0.786557,0.760277,0.736127,0.831643 /)
XEXT_COEFF_550_LKT(56,13)=149.400000 !rg=0.850757 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,14,1:6)=(/ 120.510000,122.430000,126.170000,130.700000,140.210000,143.980000 /)
XPIZA_LKT(56,14,1:6)=(/ 0.796930,0.911263,0.995693,0.996696,0.991606,0.537157 /)
XCGA_LKT(56,14,1:6)=(/ 0.865133,0.828680,0.788850,0.767397,0.742460,0.847310 /)
XEXT_COEFF_550_LKT(56,14)=123.290000 !rg=0.850757 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,15,1:6)=(/ 99.825000,101.520000,104.190000,107.130000,115.310000,117.880000 /)
XPIZA_LKT(56,15,1:6)=(/ 0.789350,0.906171,0.995358,0.996228,0.989309,0.532223 /)
XCGA_LKT(56,15,1:6)=(/ 0.869400,0.833563,0.796213,0.778040,0.750683,0.861587 /)
XEXT_COEFF_550_LKT(56,15)=102.010000 !rg=0.850757 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,16,1:6)=(/ 82.627000,83.477000,85.666000,88.686000,93.329000,96.286000 /)
XPIZA_LKT(56,16,1:6)=(/ 0.781105,0.901071,0.995155,0.995442,0.987746,0.529067 /)
XCGA_LKT(56,16,1:6)=(/ 0.872910,0.838847,0.802243,0.783813,0.758480,0.874480 /)
XEXT_COEFF_550_LKT(56,16)=84.375000 !rg=0.850757 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,17,1:6)=(/ 68.477000,69.310000,70.700000,72.537000,76.358000,78.823000 /)
XPIZA_LKT(56,17,1:6)=(/ 0.771676,0.896720,0.994972,0.995313,0.985619,0.527537 /)
XCGA_LKT(56,17,1:6)=(/ 0.875560,0.842740,0.808743,0.791020,0.768737,0.885507 /)
XEXT_COEFF_550_LKT(56,17)=69.532000 !rg=0.850757 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,18,1:6)=(/ 56.788000,57.427000,58.451000,59.669000,62.927000,64.588000 /)
XPIZA_LKT(56,18,1:6)=(/ 0.762095,0.890937,0.994745,0.994925,0.983091,0.527225 /)
XCGA_LKT(56,18,1:6)=(/ 0.878960,0.845537,0.812757,0.799073,0.777413,0.894917 /)
XEXT_COEFF_550_LKT(56,18)=57.686000 !rg=0.850757 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,19,1:6)=(/ 47.048000,47.414000,48.339000,49.446000,51.278000,52.884000 /)
XPIZA_LKT(56,19,1:6)=(/ 0.751015,0.884420,0.994434,0.994420,0.981129,0.527815 /)
XCGA_LKT(56,19,1:6)=(/ 0.882040,0.849320,0.816893,0.802390,0.783427,0.902977 /)
XEXT_COEFF_550_LKT(56,19)=47.761000 !rg=0.850757 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,20,1:6)=(/ 38.974000,39.325000,39.833000,40.607000,42.183000,43.371000 /)
XPIZA_LKT(56,20,1:6)=(/ 0.738617,0.877764,0.994109,0.994263,0.979787,0.529090 /)
XCGA_LKT(56,20,1:6)=(/ 0.884703,0.852447,0.821070,0.808967,0.792240,0.909650 /)
XEXT_COEFF_550_LKT(56,20)=39.388000 !rg=0.850757 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,1,1:6)=(/ 800.800000,898.850000,594.050000,1300.700000,902.230000,653.000000 /)
XPIZA_LKT(57,1,1:6)=(/ 0.916237,0.975080,0.998832,0.999618,0.998667,0.620994 /)
XCGA_LKT(57,1,1:6)=(/ 0.820203,0.785740,0.523073,0.766413,0.748137,0.609933 /)
XEXT_COEFF_550_LKT(57,1)=809.060000 !rg=0.921673 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,2,1:6)=(/ 715.940000,786.960000,645.560000,1173.600000,892.220000,650.070000 /)
XPIZA_LKT(57,2,1:6)=(/ 0.904574,0.968024,0.998921,0.999577,0.998668,0.632832 /)
XCGA_LKT(57,2,1:6)=(/ 0.800893,0.765810,0.595377,0.752070,0.745717,0.616993 /)
XEXT_COEFF_550_LKT(57,2)=771.680000 !rg=0.921673 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,3,1:6)=(/ 649.460000,695.070000,669.830000,987.040000,877.340000,648.710000 /)
XPIZA_LKT(57,3,1:6)=(/ 0.902990,0.966836,0.998928,0.999494,0.998627,0.642441 /)
XCGA_LKT(57,3,1:6)=(/ 0.805877,0.749693,0.652950,0.730533,0.751540,0.646183 /)
XEXT_COEFF_550_LKT(57,3)=690.480000 !rg=0.921673 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,4,1:6)=(/ 584.870000,615.300000,631.040000,816.470000,823.440000,631.710000 /)
XPIZA_LKT(57,4,1:6)=(/ 0.895619,0.964628,0.998814,0.999384,0.998511,0.646538 /)
XCGA_LKT(57,4,1:6)=(/ 0.809963,0.752903,0.684190,0.713133,0.752297,0.676547 /)
XEXT_COEFF_550_LKT(57,4)=615.840000 !rg=0.921673 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,5,1:6)=(/ 518.340000,539.440000,560.400000,677.330000,733.570000,593.150000 /)
XPIZA_LKT(57,5,1:6)=(/ 0.885391,0.960617,0.998739,0.999250,0.998300,0.643405 /)
XCGA_LKT(57,5,1:6)=(/ 0.811400,0.767340,0.702073,0.704447,0.746443,0.701347 /)
XEXT_COEFF_550_LKT(57,5)=547.090000 !rg=0.921673 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,6,1:6)=(/ 448.850000,464.830000,485.420000,564.140000,630.170000,537.800000 /)
XPIZA_LKT(57,6,1:6)=(/ 0.874365,0.955951,0.998640,0.999076,0.997985,0.633872 /)
XCGA_LKT(57,6,1:6)=(/ 0.819647,0.771453,0.725907,0.704383,0.737057,0.721033 /)
XEXT_COEFF_550_LKT(57,6)=466.490000 !rg=0.921673 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,7,1:6)=(/ 383.360000,396.890000,419.750000,467.130000,528.300000,473.200000 /)
XPIZA_LKT(57,7,1:6)=(/ 0.861493,0.952439,0.998246,0.998931,0.997536,0.619983 /)
XCGA_LKT(57,7,1:6)=(/ 0.827937,0.777707,0.731483,0.709470,0.727687,0.738393 /)
XEXT_COEFF_550_LKT(57,7)=401.470000 !rg=0.921673 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,8,1:6)=(/ 327.480000,336.340000,348.840000,387.660000,435.820000,407.380000 /)
XPIZA_LKT(57,8,1:6)=(/ 0.848802,0.944970,0.998254,0.998658,0.997034,0.604116 /)
XCGA_LKT(57,8,1:6)=(/ 0.833740,0.789767,0.750290,0.713967,0.722557,0.754883 /)
XEXT_COEFF_550_LKT(57,8)=339.450000 !rg=0.921673 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,9,1:6)=(/ 275.480000,284.280000,295.660000,320.020000,357.390000,344.540000 /)
XPIZA_LKT(57,9,1:6)=(/ 0.836496,0.937808,0.997836,0.998177,0.996347,0.587868 /)
XCGA_LKT(57,9,1:6)=(/ 0.840700,0.798160,0.750947,0.724043,0.718667,0.771647 /)
XEXT_COEFF_550_LKT(57,9)=287.790000 !rg=0.921673 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,10,1:6)=(/ 232.080000,237.240000,246.430000,266.000000,291.960000,288.310000 /)
XPIZA_LKT(57,10,1:6)=(/ 0.825257,0.931981,0.997524,0.997968,0.995555,0.572889 /)
XCGA_LKT(57,10,1:6)=(/ 0.848270,0.805080,0.765550,0.732083,0.720210,0.788640 /)
XEXT_COEFF_550_LKT(57,10)=240.610000 !rg=0.921673 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,11,1:6)=(/ 193.530000,197.490000,203.710000,217.400000,237.190000,238.770000 /)
XPIZA_LKT(57,11,1:6)=(/ 0.816486,0.924914,0.996989,0.997682,0.994759,0.559628 /)
XCGA_LKT(57,11,1:6)=(/ 0.853540,0.814420,0.772907,0.740363,0.726053,0.806063 /)
XEXT_COEFF_550_LKT(57,11)=200.140000 !rg=0.921673 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,12,1:6)=(/ 160.890000,164.130000,168.890000,179.250000,193.170000,196.530000 /)
XPIZA_LKT(57,12,1:6)=(/ 0.807702,0.919087,0.996572,0.997070,0.993662,0.548674 /)
XCGA_LKT(57,12,1:6)=(/ 0.859193,0.819207,0.783140,0.752757,0.731647,0.823263 /)
XEXT_COEFF_550_LKT(57,12)=165.600000 !rg=0.921673 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,13,1:6)=(/ 133.650000,135.720000,140.130000,145.830000,157.050000,161.100000 /)
XPIZA_LKT(57,13,1:6)=(/ 0.800555,0.913853,0.995934,0.996959,0.992405,0.540147 /)
XCGA_LKT(57,13,1:6)=(/ 0.863293,0.826557,0.785260,0.761850,0.737820,0.839707 /)
XEXT_COEFF_550_LKT(57,13)=136.880000 !rg=0.921673 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,14,1:6)=(/ 110.740000,112.890000,115.580000,118.970000,129.090000,131.960000 /)
XPIZA_LKT(57,14,1:6)=(/ 0.793458,0.908046,0.995971,0.996485,0.990398,0.534052 /)
XCGA_LKT(57,14,1:6)=(/ 0.867390,0.830330,0.796487,0.776073,0.745763,0.854757 /)
XEXT_COEFF_550_LKT(57,14)=113.810000 !rg=0.921673 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,15,1:6)=(/ 91.917000,93.192000,94.886000,98.754000,104.810000,108.000000 /)
XPIZA_LKT(57,15,1:6)=(/ 0.785313,0.903494,0.995629,0.995897,0.988948,0.530013 /)
XCGA_LKT(57,15,1:6)=(/ 0.871613,0.835997,0.800540,0.781310,0.752257,0.868307 /)
XEXT_COEFF_550_LKT(57,15)=94.368000 !rg=0.921673 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,16,1:6)=(/ 75.954000,77.116000,78.555000,81.436000,85.009000,88.208000 /)
XPIZA_LKT(57,16,1:6)=(/ 0.776496,0.899060,0.995139,0.995391,0.986552,0.527659 /)
XCGA_LKT(57,16,1:6)=(/ 0.874323,0.840700,0.805487,0.785203,0.766000,0.880407 /)
XEXT_COEFF_550_LKT(57,16)=77.656000 !rg=0.921673 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,17,1:6)=(/ 63.123000,63.686000,65.277000,66.507000,70.212000,72.222000 /)
XPIZA_LKT(57,17,1:6)=(/ 0.767975,0.894035,0.994835,0.995085,0.984227,0.526793 /)
XCGA_LKT(57,17,1:6)=(/ 0.877500,0.843903,0.810910,0.796467,0.772990,0.890643 /)
XEXT_COEFF_550_LKT(57,17)=64.293000 !rg=0.921673 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,18,1:6)=(/ 52.292000,52.843000,53.413000,55.023000,57.386000,59.198000 /)
XPIZA_LKT(57,18,1:6)=(/ 0.757013,0.888076,0.994766,0.994694,0.982626,0.527007 /)
XCGA_LKT(57,18,1:6)=(/ 0.880643,0.847270,0.815827,0.801393,0.779687,0.899300 /)
XEXT_COEFF_550_LKT(57,18)=53.210000 !rg=0.921673 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,19,1:6)=(/ 43.279000,43.811000,44.356000,45.526000,46.848000,48.490000 /)
XPIZA_LKT(57,19,1:6)=(/ 0.745431,0.881874,0.994252,0.994282,0.980349,0.528012 /)
XCGA_LKT(57,19,1:6)=(/ 0.883267,0.850810,0.818390,0.804953,0.790100,0.906657 /)
XEXT_COEFF_550_LKT(57,19)=43.987000 !rg=0.921673 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,20,1:6)=(/ 35.939000,36.167000,36.834000,37.303000,38.714000,39.786000 /)
XPIZA_LKT(57,20,1:6)=(/ 0.733499,0.873934,0.994039,0.993927,0.978897,0.529591 /)
XCGA_LKT(57,20,1:6)=(/ 0.886410,0.853310,0.822997,0.812167,0.796120,0.912707 /)
XEXT_COEFF_550_LKT(57,20)=36.437000 !rg=0.921673 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,1,1:6)=(/ 729.120000,729.770000,549.020000,1137.200000,920.090000,674.940000 /)
XPIZA_LKT(58,1,1:6)=(/ 0.910557,0.969863,0.998746,0.999564,0.998715,0.643931 /)
XCGA_LKT(58,1,1:6)=(/ 0.812550,0.790200,0.557633,0.750267,0.754203,0.626097 /)
XEXT_COEFF_550_LKT(58,1)=850.510000 !rg=0.9985 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,2,1:6)=(/ 652.730000,682.130000,611.890000,1011.200000,911.790000,667.500000 /)
XPIZA_LKT(58,2,1:6)=(/ 0.901841,0.967196,0.998861,0.999506,0.998677,0.647474 /)
XCGA_LKT(58,2,1:6)=(/ 0.803660,0.769413,0.627837,0.732507,0.757490,0.642603 /)
XEXT_COEFF_550_LKT(58,2)=727.800000 !rg=0.9985 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,3,1:6)=(/ 598.920000,621.380000,615.970000,851.340000,875.910000,656.790000 /)
XPIZA_LKT(58,3,1:6)=(/ 0.896229,0.966365,0.998855,0.999405,0.998605,0.651502 /)
XCGA_LKT(58,3,1:6)=(/ 0.810613,0.776840,0.672410,0.709577,0.759300,0.670783 /)
XEXT_COEFF_550_LKT(58,3)=636.130000 !rg=0.9985 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,4,1:6)=(/ 539.020000,555.110000,572.090000,715.010000,799.370000,626.980000 /)
XPIZA_LKT(58,4,1:6)=(/ 0.888548,0.963795,0.998738,0.999283,0.998448,0.650678 /)
XCGA_LKT(58,4,1:6)=(/ 0.815230,0.776040,0.697820,0.697057,0.755023,0.695997 /)
XEXT_COEFF_550_LKT(58,4)=563.660000 !rg=0.9985 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,5,1:6)=(/ 473.490000,493.230000,516.390000,600.150000,695.400000,576.960000 /)
XPIZA_LKT(58,5,1:6)=(/ 0.878440,0.959062,0.998519,0.999156,0.998182,0.643201 /)
XCGA_LKT(58,5,1:6)=(/ 0.820050,0.767447,0.712067,0.698320,0.745087,0.716240 /)
XEXT_COEFF_550_LKT(58,5)=498.860000 !rg=0.9985 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,6,1:6)=(/ 413.710000,422.680000,442.230000,508.340000,586.200000,514.150000 /)
XPIZA_LKT(58,6,1:6)=(/ 0.866648,0.954789,0.998510,0.998937,0.997817,0.630370 /)
XCGA_LKT(58,6,1:6)=(/ 0.826490,0.781590,0.732770,0.700773,0.733950,0.732927 /)
XEXT_COEFF_550_LKT(58,6)=432.350000 !rg=0.9985 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,7,1:6)=(/ 353.780000,361.140000,380.410000,423.940000,485.140000,446.230000 /)
XPIZA_LKT(58,7,1:6)=(/ 0.853829,0.949519,0.997965,0.998720,0.997324,0.614345 /)
XCGA_LKT(58,7,1:6)=(/ 0.833317,0.790080,0.738720,0.709187,0.724660,0.748653 /)
XEXT_COEFF_550_LKT(58,7)=367.150000 !rg=0.9985 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,8,1:6)=(/ 301.190000,307.970000,320.670000,353.130000,398.720000,380.300000 /)
XPIZA_LKT(58,8,1:6)=(/ 0.840847,0.942307,0.998080,0.998518,0.996764,0.597395 /)
XCGA_LKT(58,8,1:6)=(/ 0.839710,0.794580,0.754207,0.718487,0.719643,0.764357 /)
XEXT_COEFF_550_LKT(58,8)=313.810000 !rg=0.9985 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,9,1:6)=(/ 252.910000,259.590000,271.940000,289.270000,327.680000,319.330000 /)
XPIZA_LKT(58,9,1:6)=(/ 0.829031,0.935055,0.997497,0.998337,0.995862,0.580934 /)
XCGA_LKT(58,9,1:6)=(/ 0.846627,0.800847,0.757330,0.729977,0.717600,0.780823 /)
XEXT_COEFF_550_LKT(58,9)=263.410000 !rg=0.9985 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,10,1:6)=(/ 212.980000,218.360000,225.940000,240.120000,266.060000,265.940000 /)
XPIZA_LKT(58,10,1:6)=(/ 0.820216,0.927399,0.997100,0.997910,0.995133,0.566364 /)
XCGA_LKT(58,10,1:6)=(/ 0.851653,0.810240,0.765833,0.737620,0.720110,0.797677 /)
XEXT_COEFF_550_LKT(58,10)=220.870000 !rg=0.9985 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,11,1:6)=(/ 178.030000,181.890000,187.070000,198.430000,217.380000,219.530000 /)
XPIZA_LKT(58,11,1:6)=(/ 0.811367,0.920698,0.996931,0.997574,0.994099,0.553872 /)
XCGA_LKT(58,11,1:6)=(/ 0.857670,0.815913,0.779883,0.748830,0.725190,0.814887 /)
XEXT_COEFF_550_LKT(58,11)=184.320000 !rg=0.9985 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,12,1:6)=(/ 148.020000,150.660000,155.250000,162.330000,176.090000,180.330000 /)
XPIZA_LKT(58,12,1:6)=(/ 0.803634,0.915926,0.996368,0.997224,0.993063,0.543896 /)
XCGA_LKT(58,12,1:6)=(/ 0.861917,0.823583,0.783300,0.757727,0.732457,0.831710 /)
XEXT_COEFF_550_LKT(58,12)=152.170000 !rg=0.9985 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,13,1:6)=(/ 122.720000,124.940000,128.490000,132.510000,144.500000,147.640000 /)
XPIZA_LKT(58,13,1:6)=(/ 0.796767,0.910594,0.996042,0.996732,0.991244,0.536410 /)
XCGA_LKT(58,13,1:6)=(/ 0.865913,0.828297,0.792527,0.771843,0.740747,0.847603 /)
XEXT_COEFF_550_LKT(58,13)=126.450000 !rg=0.9985 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,14,1:6)=(/ 102.080000,103.430000,105.520000,109.870000,117.420000,120.880000 /)
XPIZA_LKT(58,14,1:6)=(/ 0.788778,0.906115,0.995836,0.996213,0.990058,0.531301 /)
XCGA_LKT(58,14,1:6)=(/ 0.870130,0.833970,0.797880,0.777970,0.746323,0.861957 /)
XEXT_COEFF_550_LKT(58,14)=104.330000 !rg=0.9985 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,15,1:6)=(/ 84.605000,85.865000,86.967000,90.381000,95.262000,98.913000 /)
XPIZA_LKT(58,15,1:6)=(/ 0.781371,0.901530,0.995158,0.995749,0.988221,0.528153 /)
XCGA_LKT(58,15,1:6)=(/ 0.872767,0.839227,0.805913,0.783177,0.761713,0.874723 /)
XEXT_COEFF_550_LKT(58,15)=86.480000 !rg=0.9985 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,16,1:6)=(/ 70.127000,71.015000,72.070000,74.395000,78.445000,80.796000 /)
XPIZA_LKT(58,16,1:6)=(/ 0.772605,0.896258,0.995190,0.994982,0.985399,0.526569 /)
XCGA_LKT(58,16,1:6)=(/ 0.876223,0.842117,0.810077,0.790997,0.767837,0.886000 /)
XEXT_COEFF_550_LKT(58,16)=71.787000 !rg=0.9985 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,17,1:6)=(/ 58.120000,58.745000,59.659000,61.289000,64.306000,66.174000 /)
XPIZA_LKT(58,17,1:6)=(/ 0.762671,0.891374,0.994741,0.994765,0.983306,0.526316 /)
XCGA_LKT(58,17,1:6)=(/ 0.878857,0.845730,0.813310,0.798820,0.773847,0.895440 /)
XEXT_COEFF_550_LKT(58,17)=58.870000 !rg=0.9985 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,18,1:6)=(/ 48.186000,48.764000,49.149000,50.428000,52.456000,54.262000 /)
XPIZA_LKT(58,18,1:6)=(/ 0.751924,0.885484,0.994577,0.994617,0.981415,0.527014 /)
XCGA_LKT(58,18,1:6)=(/ 0.881737,0.849527,0.818667,0.803507,0.786587,0.903353 /)
XEXT_COEFF_550_LKT(58,18)=48.899000 !rg=0.9985 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,19,1:6)=(/ 39.951000,40.286000,40.831000,41.676000,43.263000,44.467000 /)
XPIZA_LKT(58,19,1:6)=(/ 0.740227,0.878190,0.994312,0.994149,0.979481,0.528385 /)
XCGA_LKT(58,19,1:6)=(/ 0.884923,0.851927,0.821780,0.808480,0.792070,0.910033 /)
XEXT_COEFF_550_LKT(58,19)=40.667000 !rg=0.9985 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,20,1:6)=(/ 33.132000,33.373000,33.735000,34.458000,35.673000,36.503000 /)
XPIZA_LKT(58,20,1:6)=(/ 0.727445,0.870180,0.993821,0.993658,0.977637,0.530219 /)
XCGA_LKT(58,20,1:6)=(/ 0.887843,0.854743,0.824263,0.813393,0.796813,0.915493 /)
XEXT_COEFF_550_LKT(58,20)=33.404000 !rg=0.9985 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,1,1:6)=(/ 599.050000,565.510000,578.290000,952.260000,926.710000,678.720000 /)
XPIZA_LKT(59,1,1:6)=(/ 0.893762,0.964320,0.998784,0.999471,0.998696,0.660223 /)
XCGA_LKT(59,1,1:6)=(/ 0.793713,0.754217,0.629560,0.720770,0.759350,0.640933 /)
XEXT_COEFF_550_LKT(59,1)=773.180000 !rg=1.08173 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,2,1:6)=(/ 596.060000,592.650000,603.800000,847.150000,915.840000,676.150000 /)
XPIZA_LKT(59,2,1:6)=(/ 0.895268,0.965657,0.998837,0.999403,0.998664,0.658047 /)
XCGA_LKT(59,2,1:6)=(/ 0.809667,0.781523,0.670690,0.704150,0.766030,0.667973 /)
XEXT_COEFF_550_LKT(59,2)=665.510000 !rg=1.08173 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,3,1:6)=(/ 551.580000,566.550000,576.410000,720.170000,858.760000,656.220000 /)
XPIZA_LKT(59,3,1:6)=(/ 0.887899,0.964106,0.998803,0.999300,0.998558,0.657652 /)
XCGA_LKT(59,3,1:6)=(/ 0.811963,0.781423,0.697663,0.688960,0.763903,0.693187 /)
XEXT_COEFF_550_LKT(59,3)=581.910000 !rg=1.08173 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,4,1:6)=(/ 496.840000,510.110000,523.430000,619.770000,762.810000,614.540000 /)
XPIZA_LKT(59,4,1:6)=(/ 0.879808,0.960783,0.998701,0.999179,0.998354,0.652274 /)
XCGA_LKT(59,4,1:6)=(/ 0.817220,0.781737,0.714667,0.686597,0.755063,0.713423 /)
XEXT_COEFF_550_LKT(59,4)=516.260000 !rg=1.08173 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,5,1:6)=(/ 436.120000,447.790000,467.040000,534.920000,648.760000,555.060000 /)
XPIZA_LKT(59,5,1:6)=(/ 0.870982,0.956747,0.998365,0.999033,0.998039,0.640843 /)
XCGA_LKT(59,5,1:6)=(/ 0.825760,0.784110,0.724047,0.690800,0.742187,0.729710 /)
XEXT_COEFF_550_LKT(59,5)=454.120000 !rg=1.08173 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,6,1:6)=(/ 379.830000,391.310000,403.650000,453.370000,540.080000,487.100000 /)
XPIZA_LKT(59,6,1:6)=(/ 0.858501,0.951318,0.998406,0.998861,0.997617,0.625206 /)
XCGA_LKT(59,6,1:6)=(/ 0.829620,0.787913,0.736980,0.699640,0.730403,0.743993 /)
XEXT_COEFF_550_LKT(59,6)=397.290000 !rg=1.08173 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,7,1:6)=(/ 325.820000,332.950000,345.210000,381.200000,443.540000,417.820000 /)
XPIZA_LKT(59,7,1:6)=(/ 0.845514,0.945735,0.998142,0.998657,0.997099,0.607541 /)
XCGA_LKT(59,7,1:6)=(/ 0.836510,0.795443,0.748153,0.708597,0.722367,0.758510 /)
XEXT_COEFF_550_LKT(59,7)=337.070000 !rg=1.08173 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,8,1:6)=(/ 276.480000,283.690000,295.200000,317.980000,363.510000,353.130000 /)
XPIZA_LKT(59,8,1:6)=(/ 0.833908,0.938268,0.997349,0.998451,0.996451,0.589981 /)
XCGA_LKT(59,8,1:6)=(/ 0.843447,0.801443,0.751903,0.721260,0.717910,0.773710 /)
XEXT_COEFF_550_LKT(59,8)=287.610000 !rg=1.08173 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,9,1:6)=(/ 233.000000,237.850000,246.990000,265.020000,296.930000,294.850000 /)
XPIZA_LKT(59,9,1:6)=(/ 0.823642,0.930832,0.997221,0.998173,0.995600,0.573721 /)
XCGA_LKT(59,9,1:6)=(/ 0.850960,0.808120,0.762860,0.734833,0.716377,0.790033 /)
XEXT_COEFF_550_LKT(59,9)=241.270000 !rg=1.08173 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,10,1:6)=(/ 195.660000,200.320000,206.980000,217.810000,243.980000,244.640000 /)
XPIZA_LKT(59,10,1:6)=(/ 0.814707,0.924268,0.997076,0.997590,0.994481,0.559832 /)
XCGA_LKT(59,10,1:6)=(/ 0.855663,0.813080,0.774953,0.745280,0.720757,0.806737 /)
XEXT_COEFF_550_LKT(59,10)=202.750000 !rg=1.08173 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,11,1:6)=(/ 163.180000,167.020000,171.910000,180.360000,197.910000,201.470000 /)
XPIZA_LKT(59,11,1:6)=(/ 0.807029,0.918211,0.996509,0.997257,0.993380,0.548322 /)
XCGA_LKT(59,11,1:6)=(/ 0.860050,0.821463,0.778773,0.754827,0.725370,0.823677 /)
XEXT_COEFF_550_LKT(59,11)=168.590000 !rg=1.08173 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,12,1:6)=(/ 136.040000,138.930000,142.460000,147.440000,161.950000,165.250000 /)
XPIZA_LKT(59,12,1:6)=(/ 0.799673,0.912514,0.996297,0.997058,0.992074,0.539442 /)
XCGA_LKT(59,12,1:6)=(/ 0.864560,0.825630,0.790640,0.767193,0.735557,0.840033 /)
XEXT_COEFF_550_LKT(59,12)=140.130000 !rg=1.08173 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,13,1:6)=(/ 113.250000,114.890000,117.300000,122.270000,131.430000,135.210000 /)
XPIZA_LKT(59,13,1:6)=(/ 0.792228,0.908280,0.995939,0.996457,0.990940,0.533042 /)
XCGA_LKT(59,13,1:6)=(/ 0.868590,0.831913,0.795380,0.774373,0.740913,0.855280 /)
XEXT_COEFF_550_LKT(59,13)=115.570000 !rg=1.08173 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,14,1:6)=(/ 93.919000,95.303000,96.789000,100.990000,107.040000,110.670000 /)
XPIZA_LKT(59,14,1:6)=(/ 0.785596,0.903060,0.995580,0.995655,0.988937,0.528931 /)
XCGA_LKT(59,14,1:6)=(/ 0.871153,0.836220,0.803023,0.777727,0.754987,0.868867 /)
XEXT_COEFF_550_LKT(59,14)=95.649000 !rg=1.08173 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,15,1:6)=(/ 78.105000,79.069000,80.443000,83.365000,87.488000,90.572000 /)
XPIZA_LKT(59,15,1:6)=(/ 0.777398,0.899130,0.995193,0.995406,0.986991,0.526646 /)
XCGA_LKT(59,15,1:6)=(/ 0.874953,0.839670,0.807230,0.786013,0.765303,0.880807 /)
XEXT_COEFF_550_LKT(59,15)=79.460000 !rg=1.08173 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,16,1:6)=(/ 64.451000,65.409000,66.536000,68.443000,71.722000,74.004000 /)
XPIZA_LKT(59,16,1:6)=(/ 0.768326,0.894130,0.994875,0.994997,0.984429,0.525784 /)
XCGA_LKT(59,16,1:6)=(/ 0.877507,0.844807,0.810073,0.794957,0.769600,0.891240 /)
XEXT_COEFF_550_LKT(59,16)=65.641000 !rg=1.08173 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,17,1:6)=(/ 53.446000,54.042000,54.919000,56.453000,58.781000,60.635000 /)
XPIZA_LKT(59,17,1:6)=(/ 0.757684,0.888405,0.994743,0.994442,0.982645,0.526099 /)
XCGA_LKT(59,17,1:6)=(/ 0.880343,0.846867,0.816000,0.798677,0.780607,0.899890 /)
XEXT_COEFF_550_LKT(59,17)=54.167000 !rg=1.08173 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,18,1:6)=(/ 44.493000,44.850000,45.508000,46.578000,48.240000,49.744000 /)
XPIZA_LKT(59,18,1:6)=(/ 0.747088,0.882289,0.994456,0.994298,0.980628,0.527224 /)
XCGA_LKT(59,18,1:6)=(/ 0.883483,0.849793,0.820583,0.805083,0.790330,0.907083 /)
XEXT_COEFF_550_LKT(59,18)=45.058000 !rg=1.08173 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,19,1:6)=(/ 36.779000,37.172000,37.629000,38.413000,39.592000,40.785000 /)
XPIZA_LKT(59,19,1:6)=(/ 0.734739,0.874940,0.994012,0.993918,0.979128,0.528910 /)
XCGA_LKT(59,19,1:6)=(/ 0.886173,0.853797,0.821133,0.811217,0.794587,0.913117 /)
XEXT_COEFF_550_LKT(59,19)=37.243000 !rg=1.08173 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,20,1:6)=(/ 30.493000,30.737000,31.111000,31.665000,32.584000,33.497000 /)
XPIZA_LKT(59,20,1:6)=(/ 0.721264,0.865997,0.993755,0.993573,0.977569,0.530959 /)
XCGA_LKT(59,20,1:6)=(/ 0.889247,0.855693,0.826253,0.814017,0.802570,0.918020 /)
XEXT_COEFF_550_LKT(59,20)=30.780000 !rg=1.08173 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,1,1:6)=(/ 540.430000,533.570000,628.680000,767.610000,945.550000,680.900000 /)
XPIZA_LKT(60,1,1:6)=(/ 0.887653,0.953102,0.998878,0.999340,0.998684,0.666349 /)
XCGA_LKT(60,1,1:6)=(/ 0.800720,0.757233,0.706250,0.683500,0.772890,0.670043 /)
XEXT_COEFF_550_LKT(60,1)=645.350000 !rg=1.1719 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,2,1:6)=(/ 547.140000,553.730000,602.900000,687.940000,902.900000,677.260000 /)
XPIZA_LKT(60,2,1:6)=(/ 0.887515,0.960515,0.998813,0.999266,0.998628,0.664638 /)
XCGA_LKT(60,2,1:6)=(/ 0.812347,0.782247,0.709617,0.669163,0.771803,0.694207 /)
XEXT_COEFF_550_LKT(60,2)=582.680000 !rg=1.1719 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,3,1:6)=(/ 505.170000,514.920000,546.170000,613.530000,825.620000,646.960000 /)
XPIZA_LKT(60,3,1:6)=(/ 0.882346,0.961262,0.998668,0.999164,0.998483,0.660809 /)
XCGA_LKT(60,3,1:6)=(/ 0.820480,0.779747,0.716367,0.666087,0.765210,0.713483 /)
XEXT_COEFF_550_LKT(60,3)=537.330000 !rg=1.1719 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,4,1:6)=(/ 454.860000,466.450000,488.070000,544.920000,715.140000,594.820000 /)
XPIZA_LKT(60,4,1:6)=(/ 0.873855,0.957493,0.998546,0.999045,0.998224,0.651254 /)
XCGA_LKT(60,4,1:6)=(/ 0.825400,0.781600,0.726400,0.677420,0.751970,0.728943 /)
XEXT_COEFF_550_LKT(60,4)=479.390000 !rg=1.1719 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,5,1:6)=(/ 402.260000,411.130000,424.520000,474.090000,598.020000,528.250000 /)
XPIZA_LKT(60,5,1:6)=(/ 0.861989,0.953728,0.998342,0.998888,0.997855,0.636334 /)
XCGA_LKT(60,5,1:6)=(/ 0.828873,0.790213,0.738353,0.689427,0.737503,0.741907 /)
XEXT_COEFF_550_LKT(60,5)=415.750000 !rg=1.1719 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,6,1:6)=(/ 348.580000,359.210000,373.910000,405.680000,493.950000,457.530000 /)
XPIZA_LKT(60,6,1:6)=(/ 0.850247,0.947388,0.998154,0.998794,0.997351,0.618452 /)
XCGA_LKT(60,6,1:6)=(/ 0.835863,0.788037,0.745197,0.707603,0.724223,0.754367 /)
XEXT_COEFF_550_LKT(60,6)=364.050000 !rg=1.1719 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,7,1:6)=(/ 299.000000,307.120000,318.350000,345.570000,404.310000,388.680000 /)
XPIZA_LKT(60,7,1:6)=(/ 0.838427,0.940149,0.998004,0.998453,0.996735,0.599699 /)
XCGA_LKT(60,7,1:6)=(/ 0.842847,0.796497,0.756900,0.716347,0.717117,0.768080 /)
XEXT_COEFF_550_LKT(60,7)=311.900000 !rg=1.1719 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,8,1:6)=(/ 253.690000,260.360000,270.540000,286.500000,331.720000,326.370000 /)
XPIZA_LKT(60,8,1:6)=(/ 0.827226,0.934228,0.997529,0.998386,0.995982,0.582049 /)
XCGA_LKT(60,8,1:6)=(/ 0.848273,0.803070,0.763130,0.731963,0.715330,0.783037 /)
XEXT_COEFF_550_LKT(60,8)=264.820000 !rg=1.1719 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,9,1:6)=(/ 214.540000,218.630000,224.880000,241.290000,268.980000,271.330000 /)
XPIZA_LKT(60,9,1:6)=(/ 0.817546,0.927744,0.997030,0.997719,0.995361,0.566361 /)
XCGA_LKT(60,9,1:6)=(/ 0.853570,0.813960,0.771733,0.733203,0.720940,0.799287 /)
XEXT_COEFF_550_LKT(60,9)=221.300000 !rg=1.1719 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,10,1:6)=(/ 180.120000,183.530000,188.470000,200.820000,220.960000,224.530000 /)
XPIZA_LKT(60,10,1:6)=(/ 0.809837,0.920831,0.996936,0.997536,0.994191,0.553423 /)
XCGA_LKT(60,10,1:6)=(/ 0.859420,0.818140,0.777740,0.747533,0.720280,0.815817 /)
XEXT_COEFF_550_LKT(60,10)=185.900000 !rg=1.1719 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,11,1:6)=(/ 150.120000,152.890000,157.850000,164.680000,180.840000,184.610000 /)
XPIZA_LKT(60,11,1:6)=(/ 0.801813,0.915662,0.996356,0.997192,0.992881,0.543065 /)
XCGA_LKT(60,11,1:6)=(/ 0.863770,0.824063,0.784967,0.758430,0.730343,0.832397 /)
XEXT_COEFF_550_LKT(60,11)=154.700000 !rg=1.1719 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,12,1:6)=(/ 125.340000,127.220000,130.090000,136.240000,147.030000,151.300000 /)
XPIZA_LKT(60,12,1:6)=(/ 0.795433,0.910468,0.995982,0.996758,0.991776,0.535364 /)
XCGA_LKT(60,12,1:6)=(/ 0.867410,0.829803,0.791410,0.768750,0.735597,0.848183 /)
XEXT_COEFF_550_LKT(60,12)=128.510000 !rg=1.1719 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,13,1:6)=(/ 104.030000,105.710000,107.700000,112.330000,119.850000,123.750000 /)
XPIZA_LKT(60,13,1:6)=(/ 0.789110,0.905182,0.995582,0.996096,0.989882,0.530083 /)
XCGA_LKT(60,13,1:6)=(/ 0.869997,0.833753,0.799600,0.774097,0.749147,0.862697 /)
XEXT_COEFF_550_LKT(60,13)=106.000000 !rg=1.1719 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,14,1:6)=(/ 86.421000,87.544000,89.579000,92.642000,97.585000,101.300000 /)
XPIZA_LKT(60,14,1:6)=(/ 0.781628,0.901233,0.994915,0.995790,0.988512,0.526949 /)
XCGA_LKT(60,14,1:6)=(/ 0.873350,0.838453,0.804267,0.785060,0.760827,0.875453 /)
XEXT_COEFF_550_LKT(60,14)=88.165000 !rg=1.1719 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,15,1:6)=(/ 71.723000,72.686000,73.780000,76.200000,80.279000,82.926000 /)
XPIZA_LKT(60,15,1:6)=(/ 0.772824,0.896320,0.995051,0.994877,0.985257,0.525484 /)
XCGA_LKT(60,15,1:6)=(/ 0.876450,0.842410,0.808947,0.787933,0.765687,0.886530 /)
XEXT_COEFF_550_LKT(60,15)=73.316000 !rg=1.1719 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,16,1:6)=(/ 59.345000,59.983000,61.135000,62.765000,65.923000,67.783000 /)
XPIZA_LKT(60,16,1:6)=(/ 0.763047,0.891808,0.994747,0.994836,0.983081,0.525297 /)
XCGA_LKT(60,16,1:6)=(/ 0.879230,0.846020,0.813713,0.797107,0.775247,0.896117 /)
XEXT_COEFF_550_LKT(60,16)=60.503000 !rg=1.1719 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,17,1:6)=(/ 49.323000,49.734000,50.578000,51.855000,53.883000,55.566000 /)
XPIZA_LKT(60,17,1:6)=(/ 0.753352,0.885998,0.994506,0.994562,0.981311,0.526119 /)
XCGA_LKT(60,17,1:6)=(/ 0.882010,0.849487,0.819127,0.804220,0.784200,0.903993 /)
XEXT_COEFF_550_LKT(60,17)=50.072000 !rg=1.1719 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,18,1:6)=(/ 40.930000,41.299000,41.689000,42.678000,44.315000,45.609000 /)
XPIZA_LKT(60,18,1:6)=(/ 0.741257,0.879092,0.994066,0.993995,0.979501,0.527619 /)
XCGA_LKT(60,18,1:6)=(/ 0.884827,0.851870,0.821007,0.806673,0.791443,0.910497 /)
XEXT_COEFF_550_LKT(60,18)=41.532000 !rg=1.1719 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,19,1:6)=(/ 33.885000,34.137000,34.646000,35.295000,36.535000,37.415000 /)
XPIZA_LKT(60,19,1:6)=(/ 0.728445,0.871129,0.993879,0.993928,0.977396,0.529573 /)
XCGA_LKT(60,19,1:6)=(/ 0.887663,0.854720,0.824713,0.812577,0.797870,0.915917 /)
XEXT_COEFF_550_LKT(60,19)=34.339000 !rg=1.1719 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,20,1:6)=(/ 28.149000,28.307000,28.652000,29.195000,30.037000,30.744000 /)
XPIZA_LKT(60,20,1:6)=(/ 0.715902,0.862004,0.993402,0.993473,0.976070,0.531793 /)
XCGA_LKT(60,20,1:6)=(/ 0.890843,0.857577,0.828240,0.818277,0.805047,0.920303 /)
XEXT_COEFF_550_LKT(60,20)=28.450000 !rg=1.1719 sigma=2.95 wvl=0.55
 
END SUBROUTINE SALT_OPT_LKT_SET6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET7()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE
 
 
XEXT_COEFF_WVL_LKT(61,1,1:6)=(/ 543.860000,533.850000,645.780000,596.720000,917.960000,689.090000 /)
XPIZA_LKT(61,1,1:6)=(/ 0.889026,0.961030,0.998933,0.999152,0.998660,0.667581 /)
XCGA_LKT(61,1,1:6)=(/ 0.825880,0.771323,0.765880,0.632543,0.781003,0.710170 /)
XEXT_COEFF_550_LKT(61,1)=502.570000 !rg=1.26958 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,2,1:6)=(/ 501.110000,515.150000,579.990000,560.960000,871.480000,670.050000 /)
XPIZA_LKT(61,2,1:6)=(/ 0.880341,0.958319,0.998742,0.999087,0.998562,0.667992 /)
XCGA_LKT(61,2,1:6)=(/ 0.819377,0.772020,0.741373,0.631253,0.774173,0.719020 /)
XEXT_COEFF_550_LKT(61,2)=526.580000 !rg=1.26958 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,3,1:6)=(/ 461.730000,478.370000,504.770000,521.220000,777.060000,629.170000 /)
XPIZA_LKT(61,3,1:6)=(/ 0.875278,0.956862,0.998594,0.999019,0.998372,0.660919 /)
XCGA_LKT(61,3,1:6)=(/ 0.825073,0.787433,0.729970,0.648277,0.763027,0.731500 /)
XEXT_COEFF_550_LKT(61,3)=491.720000 !rg=1.26958 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,4,1:6)=(/ 416.060000,431.390000,447.490000,478.260000,658.170000,568.540000 /)
XPIZA_LKT(61,4,1:6)=(/ 0.865908,0.953595,0.998423,0.998922,0.998058,0.647557 /)
XCGA_LKT(61,4,1:6)=(/ 0.829227,0.789910,0.731253,0.672923,0.746667,0.742660 /)
XEXT_COEFF_550_LKT(61,4)=440.370000 !rg=1.26958 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,5,1:6)=(/ 368.660000,377.540000,394.460000,427.230000,544.250000,497.550000 /)
XPIZA_LKT(61,5,1:6)=(/ 0.854614,0.950162,0.998219,0.998728,0.997615,0.629710 /)
XCGA_LKT(61,5,1:6)=(/ 0.835787,0.791557,0.745463,0.691863,0.729917,0.752983 /)
XEXT_COEFF_550_LKT(61,5)=386.240000 !rg=1.26958 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,6,1:6)=(/ 320.050000,329.780000,341.790000,366.420000,446.260000,426.340000 /)
XPIZA_LKT(61,6,1:6)=(/ 0.841576,0.943073,0.998100,0.998628,0.997069,0.610229 /)
XCGA_LKT(61,6,1:6)=(/ 0.840133,0.796450,0.745097,0.709113,0.719743,0.764207 /)
XEXT_COEFF_550_LKT(61,6)=332.150000 !rg=1.26958 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,7,1:6)=(/ 273.680000,281.820000,292.330000,310.990000,365.410000,359.470000 /)
XPIZA_LKT(61,7,1:6)=(/ 0.830887,0.937257,0.997551,0.998323,0.996408,0.591002 /)
XCGA_LKT(61,7,1:6)=(/ 0.846430,0.804560,0.753650,0.722297,0.714643,0.777500 /)
XEXT_COEFF_550_LKT(61,7)=284.670000 !rg=1.26958 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,8,1:6)=(/ 233.730000,238.800000,246.410000,262.840000,299.650000,300.420000 /)
XPIZA_LKT(61,8,1:6)=(/ 0.820230,0.930215,0.997515,0.997870,0.995573,0.573770 /)
XCGA_LKT(61,8,1:6)=(/ 0.852810,0.809520,0.765307,0.733643,0.712880,0.792413 /)
XEXT_COEFF_550_LKT(61,8)=241.870000 !rg=1.26958 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,9,1:6)=(/ 197.250000,201.450000,208.180000,220.720000,244.790000,249.000000 /)
XPIZA_LKT(61,9,1:6)=(/ 0.812003,0.923777,0.996901,0.997704,0.994831,0.559003 /)
XCGA_LKT(61,9,1:6)=(/ 0.857813,0.815787,0.774637,0.741353,0.719940,0.808593 /)
XEXT_COEFF_550_LKT(61,9)=203.830000 !rg=1.26958 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,10,1:6)=(/ 165.910000,168.950000,171.790000,182.950000,199.920000,205.690000 /)
XPIZA_LKT(61,10,1:6)=(/ 0.805329,0.917225,0.996763,0.997404,0.993842,0.547264 /)
XCGA_LKT(61,10,1:6)=(/ 0.861533,0.822797,0.786967,0.749347,0.728997,0.824887 /)
XEXT_COEFF_550_LKT(61,10)=170.310000 !rg=1.26958 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,11,1:6)=(/ 138.420000,140.510000,144.730000,151.370000,163.750000,168.960000 /)
XPIZA_LKT(61,11,1:6)=(/ 0.798669,0.911845,0.995874,0.996939,0.992484,0.538180 /)
XCGA_LKT(61,11,1:6)=(/ 0.866057,0.828493,0.786070,0.764423,0.732650,0.840997 /)
XEXT_COEFF_550_LKT(61,11)=142.050000 !rg=1.26958 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,12,1:6)=(/ 115.370000,117.410000,119.070000,124.880000,133.630000,138.420000 /)
XPIZA_LKT(61,12,1:6)=(/ 0.792225,0.906950,0.995811,0.996183,0.990863,0.531718 /)
XCGA_LKT(61,12,1:6)=(/ 0.868513,0.832567,0.798093,0.769297,0.745260,0.856107 /)
XEXT_COEFF_550_LKT(61,12)=117.640000 !rg=1.26958 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,13,1:6)=(/ 95.836000,96.974000,99.163000,102.920000,109.080000,113.220000 /)
XPIZA_LKT(61,13,1:6)=(/ 0.785603,0.903319,0.995503,0.996050,0.989207,0.527551 /)
XCGA_LKT(61,13,1:6)=(/ 0.872280,0.837070,0.803423,0.781990,0.754797,0.869803 /)
XEXT_COEFF_550_LKT(61,13)=97.953000 !rg=1.26958 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,14,1:6)=(/ 79.662000,80.588000,82.464000,84.413000,89.343000,92.711000 /)
XPIZA_LKT(61,14,1:6)=(/ 0.777499,0.899536,0.995042,0.995486,0.987350,0.525358 /)
XCGA_LKT(61,14,1:6)=(/ 0.874827,0.841327,0.804113,0.786983,0.762560,0.881680 /)
XEXT_COEFF_550_LKT(61,14)=80.885000 !rg=1.26958 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,15,1:6)=(/ 66.120000,66.966000,68.089000,69.676000,73.756000,75.925000 /)
XPIZA_LKT(61,15,1:6)=(/ 0.768883,0.894156,0.994934,0.994745,0.984525,0.524665 /)
XCGA_LKT(61,15,1:6)=(/ 0.877993,0.844350,0.811257,0.794927,0.771457,0.891883 /)
XEXT_COEFF_550_LKT(61,15)=67.324000 !rg=1.26958 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,16,1:6)=(/ 54.812000,55.219000,56.239000,57.572000,59.863000,62.093000 /)
XPIZA_LKT(61,16,1:6)=(/ 0.759315,0.888463,0.994664,0.994703,0.983035,0.525087 /)
XCGA_LKT(61,16,1:6)=(/ 0.880673,0.848130,0.814257,0.800803,0.779313,0.900630 /)
XEXT_COEFF_550_LKT(61,16)=55.686000 !rg=1.26958 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,17,1:6)=(/ 45.407000,45.913000,46.650000,47.476000,49.356000,50.929000 /)
XPIZA_LKT(61,17,1:6)=(/ 0.747628,0.883030,0.994246,0.994418,0.980748,0.526358 /)
XCGA_LKT(61,17,1:6)=(/ 0.883190,0.851413,0.818637,0.806567,0.787670,0.907753 /)
XEXT_COEFF_550_LKT(61,17)=45.959000 !rg=1.26958 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,18,1:6)=(/ 37.717000,38.069000,38.512000,39.151000,40.804000,41.826000 /)
XPIZA_LKT(61,18,1:6)=(/ 0.735884,0.875566,0.994153,0.993982,0.978513,0.528180 /)
XCGA_LKT(61,18,1:6)=(/ 0.886193,0.853230,0.823170,0.812050,0.795337,0.913600 /)
XEXT_COEFF_550_LKT(61,18)=38.273000 !rg=1.26958 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,19,1:6)=(/ 31.281000,31.440000,31.942000,32.420000,33.305000,34.330000 /)
XPIZA_LKT(61,19,1:6)=(/ 0.723068,0.866818,0.993699,0.993728,0.977617,0.530352 /)
XCGA_LKT(61,19,1:6)=(/ 0.889300,0.856263,0.825400,0.815310,0.801297,0.918447 /)
XEXT_COEFF_550_LKT(61,19)=31.597000 !rg=1.26958 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,20,1:6)=(/ 25.926000,26.137000,26.376000,26.755000,27.522000,28.224000 /)
XPIZA_LKT(61,20,1:6)=(/ 0.709316,0.857571,0.993253,0.993210,0.975867,0.532704 /)
XCGA_LKT(61,20,1:6)=(/ 0.892153,0.859303,0.828500,0.819917,0.807493,0.922353 /)
XEXT_COEFF_550_LKT(61,20)=26.143000 !rg=1.26958 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,1,1:6)=(/ 492.460000,534.120000,620.920000,464.380000,880.680000,685.560000 /)
XPIZA_LKT(62,1,1:6)=(/ 0.878567,0.963010,0.998778,0.998882,0.998553,0.671598 /)
XCGA_LKT(62,1,1:6)=(/ 0.829823,0.801230,0.775657,0.563090,0.779343,0.741020 /)
XEXT_COEFF_550_LKT(62,1)=427.840000 !rg=1.37541 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,2,1:6)=(/ 459.680000,478.750000,531.440000,461.480000,821.410000,653.020000 /)
XPIZA_LKT(62,2,1:6)=(/ 0.871678,0.956376,0.998602,0.998883,0.998458,0.668456 /)
XCGA_LKT(62,2,1:6)=(/ 0.825560,0.788260,0.748777,0.595223,0.772347,0.739813 /)
XEXT_COEFF_550_LKT(62,2)=474.000000 !rg=1.37541 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,3,1:6)=(/ 424.120000,438.240000,471.110000,456.990000,715.340000,603.310000 /)
XPIZA_LKT(62,3,1:6)=(/ 0.865987,0.955633,0.998215,0.998879,0.998212,0.657895 /)
XCGA_LKT(62,3,1:6)=(/ 0.829303,0.783703,0.741710,0.644820,0.756800,0.747073 /)
XEXT_COEFF_550_LKT(62,3)=447.550000 !rg=1.37541 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,4,1:6)=(/ 381.850000,393.880000,416.600000,427.670000,596.670000,536.740000 /)
XPIZA_LKT(62,4,1:6)=(/ 0.856833,0.951994,0.998121,0.998810,0.997829,0.641143 /)
XCGA_LKT(62,4,1:6)=(/ 0.834497,0.789200,0.742850,0.676493,0.737597,0.754703 /)
XEXT_COEFF_550_LKT(62,4)=401.430000 !rg=1.37541 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,5,1:6)=(/ 337.260000,348.510000,360.710000,382.220000,489.310000,464.070000 /)
XPIZA_LKT(62,5,1:6)=(/ 0.846117,0.945772,0.998112,0.998558,0.997343,0.621057 /)
XCGA_LKT(62,5,1:6)=(/ 0.840007,0.798533,0.745867,0.696940,0.723060,0.763140 /)
XEXT_COEFF_550_LKT(62,5)=353.230000 !rg=1.37541 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,6,1:6)=(/ 293.490000,301.260000,310.750000,331.600000,403.100000,394.420000 /)
XPIZA_LKT(62,6,1:6)=(/ 0.834336,0.939599,0.997876,0.998432,0.996688,0.600712 /)
XCGA_LKT(62,6,1:6)=(/ 0.845657,0.800177,0.762167,0.712187,0.714687,0.773687 /)
XEXT_COEFF_550_LKT(62,6)=301.670000 !rg=1.37541 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,7,1:6)=(/ 251.320000,257.430000,268.180000,282.240000,331.360000,330.790000 /)
XPIZA_LKT(62,7,1:6)=(/ 0.823048,0.933965,0.997539,0.998219,0.995908,0.581672 /)
XCGA_LKT(62,7,1:6)=(/ 0.851877,0.806947,0.763770,0.726673,0.711730,0.786900 /)
XEXT_COEFF_550_LKT(62,7)=260.720000 !rg=1.37541 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,8,1:6)=(/ 214.850000,219.370000,223.980000,239.360000,270.990000,275.570000 /)
XPIZA_LKT(62,8,1:6)=(/ 0.815253,0.925739,0.997363,0.997790,0.995173,0.565336 /)
XCGA_LKT(62,8,1:6)=(/ 0.855737,0.813307,0.777470,0.734347,0.716920,0.801867 /)
XEXT_COEFF_550_LKT(62,8)=220.460000 !rg=1.37541 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,9,1:6)=(/ 181.050000,185.320000,190.530000,200.460000,221.910000,227.980000 /)
XPIZA_LKT(62,9,1:6)=(/ 0.807192,0.919405,0.996859,0.997376,0.993963,0.551812 /)
XCGA_LKT(62,9,1:6)=(/ 0.860823,0.820157,0.776947,0.745767,0.719270,0.817953 /)
XEXT_COEFF_550_LKT(62,9)=187.260000 !rg=1.37541 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,10,1:6)=(/ 152.500000,155.370000,159.290000,168.250000,181.830000,188.150000 /)
XPIZA_LKT(62,10,1:6)=(/ 0.800367,0.915006,0.996034,0.996484,0.993082,0.541458 /)
XCGA_LKT(62,10,1:6)=(/ 0.864587,0.824950,0.788423,0.755373,0.730637,0.833903 /)
XEXT_COEFF_550_LKT(62,10)=156.200000 !rg=1.37541 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,11,1:6)=(/ 127.470000,129.470000,132.040000,138.470000,148.170000,154.500000 /)
XPIZA_LKT(62,11,1:6)=(/ 0.794740,0.909599,0.995870,0.996718,0.992143,0.533745 /)
XCGA_LKT(62,11,1:6)=(/ 0.867967,0.831773,0.794323,0.764423,0.741527,0.849420 /)
XEXT_COEFF_550_LKT(62,11)=131.100000 !rg=1.37541 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,12,1:6)=(/ 106.090000,107.730000,110.070000,114.470000,121.580000,126.580000 /)
XPIZA_LKT(62,12,1:6)=(/ 0.788251,0.905424,0.995721,0.996255,0.990516,0.528536 /)
XCGA_LKT(62,12,1:6)=(/ 0.871010,0.834723,0.800577,0.777133,0.750440,0.863753 /)
XEXT_COEFF_550_LKT(62,12)=108.340000 !rg=1.37541 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,13,1:6)=(/ 88.237000,89.229000,91.305000,93.907000,99.621000,103.570000 /)
XPIZA_LKT(62,13,1:6)=(/ 0.781207,0.901350,0.995303,0.995818,0.988711,0.525457 /)
XCGA_LKT(62,13,1:6)=(/ 0.873657,0.839843,0.804163,0.783120,0.757413,0.876557 /)
XEXT_COEFF_550_LKT(62,13)=89.727000 !rg=1.37541 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,14,1:6)=(/ 73.233000,74.282000,75.620000,77.144000,82.280000,84.844000 /)
XPIZA_LKT(62,14,1:6)=(/ 0.773555,0.896375,0.994971,0.995372,0.985692,0.524157 /)
XCGA_LKT(62,14,1:6)=(/ 0.876440,0.842413,0.810047,0.794477,0.766360,0.887527 /)
XEXT_COEFF_550_LKT(62,14)=74.781000 !rg=1.37541 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,15,1:6)=(/ 60.921000,61.460000,62.237000,64.293000,67.205000,69.521000 /)
XPIZA_LKT(62,15,1:6)=(/ 0.764621,0.891976,0.994908,0.994823,0.983927,0.524166 /)
XCGA_LKT(62,15,1:6)=(/ 0.879590,0.846003,0.814307,0.797433,0.773650,0.896850 /)
XEXT_COEFF_550_LKT(62,15)=62.058000 !rg=1.37541 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,16,1:6)=(/ 50.415000,50.994000,51.592000,53.210000,54.639000,56.888000 /)
XPIZA_LKT(62,16,1:6)=(/ 0.753969,0.886289,0.994532,0.994090,0.982118,0.525135 /)
XCGA_LKT(62,16,1:6)=(/ 0.881947,0.849540,0.817307,0.800630,0.786933,0.904773 /)
XEXT_COEFF_550_LKT(62,16)=51.352000 !rg=1.37541 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,17,1:6)=(/ 41.913000,42.154000,42.922000,43.672000,45.301000,46.687000 /)
XPIZA_LKT(62,17,1:6)=(/ 0.742931,0.879759,0.994263,0.993928,0.980088,0.526797 /)
XCGA_LKT(62,17,1:6)=(/ 0.884803,0.852237,0.821087,0.809247,0.792723,0.911183 /)
XEXT_COEFF_550_LKT(62,17)=42.514000 !rg=1.37541 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,18,1:6)=(/ 34.784000,34.997000,35.297000,36.172000,37.293000,38.364000 /)
XPIZA_LKT(62,18,1:6)=(/ 0.730446,0.872263,0.994050,0.993739,0.978020,0.528884 /)
XCGA_LKT(62,18,1:6)=(/ 0.887807,0.854413,0.825077,0.813307,0.797760,0.916410 /)
XEXT_COEFF_550_LKT(62,18)=35.248000 !rg=1.37541 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,19,1:6)=(/ 28.796000,29.039000,29.304000,29.963000,30.499000,31.506000 /)
XPIZA_LKT(62,19,1:6)=(/ 0.716900,0.863114,0.993505,0.993424,0.976777,0.531230 /)
XCGA_LKT(62,19,1:6)=(/ 0.890507,0.857413,0.826560,0.815667,0.806847,0.920727 /)
XEXT_COEFF_550_LKT(62,19)=29.129000 !rg=1.37541 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,20,1:6)=(/ 23.931000,24.039000,24.377000,24.644000,25.289000,25.916000 /)
XPIZA_LKT(62,20,1:6)=(/ 0.703538,0.852596,0.992996,0.993035,0.975325,0.533677 /)
XCGA_LKT(62,20,1:6)=(/ 0.893863,0.860067,0.829897,0.821520,0.811270,0.924190 /)
XEXT_COEFF_550_LKT(62,20)=24.179000 !rg=1.37541 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,1,1:6)=(/ 433.700000,471.130000,513.880000,368.860000,825.550000,656.570000 /)
XPIZA_LKT(63,1,1:6)=(/ 0.862474,0.959796,0.998693,0.998637,0.998470,0.674541 /)
XCGA_LKT(63,1,1:6)=(/ 0.824233,0.798883,0.770900,0.533183,0.775630,0.753910 /)
XEXT_COEFF_550_LKT(63,1)=434.150000 !rg=1.49006 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,2,1:6)=(/ 420.990000,437.080000,463.600000,401.520000,753.660000,626.110000 /)
XPIZA_LKT(63,2,1:6)=(/ 0.865828,0.954590,0.998467,0.998721,0.998304,0.665664 /)
XCGA_LKT(63,2,1:6)=(/ 0.832270,0.790147,0.769037,0.595283,0.765703,0.756080 /)
XEXT_COEFF_550_LKT(63,2)=429.560000 !rg=1.49006 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,3,1:6)=(/ 389.980000,400.450000,421.670000,415.120000,642.930000,570.250000 /)
XPIZA_LKT(63,3,1:6)=(/ 0.858683,0.952047,0.998344,0.998735,0.997998,0.651600 /)
XCGA_LKT(63,3,1:6)=(/ 0.835820,0.792890,0.754890,0.647113,0.746520,0.760193 /)
XEXT_COEFF_550_LKT(63,3)=404.750000 !rg=1.49006 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,4,1:6)=(/ 351.730000,359.950000,376.610000,390.300000,530.740000,500.700000 /)
XPIZA_LKT(63,4,1:6)=(/ 0.848786,0.947976,0.997809,0.998664,0.997555,0.632021 /)
XCGA_LKT(63,4,1:6)=(/ 0.840487,0.797443,0.754177,0.681837,0.726877,0.765257 /)
XEXT_COEFF_550_LKT(63,4)=364.510000 !rg=1.49006 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,5,1:6)=(/ 309.840000,317.890000,334.060000,345.320000,438.350000,428.990000 /)
XPIZA_LKT(63,5,1:6)=(/ 0.836347,0.943479,0.997458,0.998559,0.996972,0.610539 /)
XCGA_LKT(63,5,1:6)=(/ 0.845160,0.800323,0.754853,0.706267,0.713473,0.772603 /)
XEXT_COEFF_550_LKT(63,5)=323.600000 !rg=1.49006 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,6,1:6)=(/ 269.980000,274.960000,283.890000,304.880000,360.380000,362.610000 /)
XPIZA_LKT(63,6,1:6)=(/ 0.826750,0.936629,0.997781,0.998229,0.996329,0.590149 /)
XCGA_LKT(63,6,1:6)=(/ 0.850717,0.807847,0.771003,0.718483,0.709573,0.783000 /)
XEXT_COEFF_550_LKT(63,6)=279.680000 !rg=1.49006 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,7,1:6)=(/ 231.830000,235.860000,245.210000,259.280000,296.980000,303.130000 /)
XPIZA_LKT(63,7,1:6)=(/ 0.816959,0.929229,0.997126,0.998020,0.995578,0.571954 /)
XCGA_LKT(63,7,1:6)=(/ 0.855987,0.813363,0.767097,0.733893,0.709853,0.796387 /)
XEXT_COEFF_550_LKT(63,7)=238.690000 !rg=1.49006 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,8,1:6)=(/ 197.240000,201.140000,206.650000,219.080000,244.420000,252.090000 /)
XPIZA_LKT(63,8,1:6)=(/ 0.809480,0.922902,0.997090,0.997811,0.994826,0.556944 /)
XCGA_LKT(63,8,1:6)=(/ 0.859690,0.817883,0.781533,0.744877,0.717250,0.811427 /)
XEXT_COEFF_550_LKT(63,8)=203.850000 !rg=1.49006 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,9,1:6)=(/ 166.580000,169.780000,175.680000,181.460000,202.780000,208.370000 /)
XPIZA_LKT(63,9,1:6)=(/ 0.802483,0.916827,0.996014,0.997510,0.993213,0.544950 /)
XCGA_LKT(63,9,1:6)=(/ 0.863977,0.823810,0.781420,0.756313,0.722697,0.827330 /)
XEXT_COEFF_550_LKT(63,9)=171.120000 !rg=1.49006 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,10,1:6)=(/ 140.220000,142.570000,146.330000,152.090000,165.260000,171.930000 /)
XPIZA_LKT(63,10,1:6)=(/ 0.796800,0.911738,0.996184,0.997056,0.992528,0.536110 /)
XCGA_LKT(63,10,1:6)=(/ 0.867130,0.828820,0.789307,0.762610,0.731977,0.842800 /)
XEXT_COEFF_550_LKT(63,10)=143.880000 !rg=1.49006 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,11,1:6)=(/ 117.410000,119.260000,121.580000,126.890000,136.080000,141.190000 /)
XPIZA_LKT(63,11,1:6)=(/ 0.790902,0.906686,0.995777,0.996337,0.991006,0.529814 /)
XCGA_LKT(63,11,1:6)=(/ 0.870463,0.833260,0.798157,0.771377,0.742853,0.857607 /)
XEXT_COEFF_550_LKT(63,11)=120.520000 !rg=1.49006 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,12,1:6)=(/ 97.817000,99.041000,101.690000,104.170000,111.120000,115.720000 /)
XPIZA_LKT(63,12,1:6)=(/ 0.784563,0.903266,0.995209,0.996052,0.989593,0.525844 /)
XCGA_LKT(63,12,1:6)=(/ 0.872877,0.838047,0.799063,0.780227,0.752100,0.871063 /)
XEXT_COEFF_550_LKT(63,12)=99.467000 !rg=1.49006 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,13,1:6)=(/ 81.150000,82.115000,84.213000,85.823000,91.654000,94.732000 /)
XPIZA_LKT(63,13,1:6)=(/ 0.777548,0.898832,0.994830,0.995432,0.987160,0.523807 /)
XCGA_LKT(63,13,1:6)=(/ 0.875403,0.841403,0.806757,0.790787,0.761687,0.882923 /)
XEXT_COEFF_550_LKT(63,13)=83.006000 !rg=1.49006 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,14,1:6)=(/ 67.639000,68.370000,69.320000,71.323000,75.093000,77.651000 /)
XPIZA_LKT(63,14,1:6)=(/ 0.769154,0.894768,0.994916,0.995116,0.985041,0.523327 /)
XCGA_LKT(63,14,1:6)=(/ 0.878333,0.845113,0.812583,0.796583,0.768063,0.892970 /)
XEXT_COEFF_550_LKT(63,14)=68.634000 !rg=1.49006 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,15,1:6)=(/ 56.107000,56.732000,57.242000,58.813000,61.304000,63.665000 /)
XPIZA_LKT(63,15,1:6)=(/ 0.759733,0.889171,0.994642,0.994849,0.982873,0.523973 /)
XCGA_LKT(63,15,1:6)=(/ 0.880540,0.848323,0.816807,0.800623,0.781500,0.901423 /)
XEXT_COEFF_550_LKT(63,15)=56.936000 !rg=1.49006 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,16,1:6)=(/ 46.510000,46.882000,47.505000,48.515000,50.485000,52.130000 /)
XPIZA_LKT(63,16,1:6)=(/ 0.749015,0.883238,0.994538,0.994461,0.981107,0.525420 /)
XCGA_LKT(63,16,1:6)=(/ 0.883527,0.850537,0.820530,0.806997,0.788900,0.908557 /)
XEXT_COEFF_550_LKT(63,16)=47.282000 !rg=1.49006 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,17,1:6)=(/ 38.594000,38.925000,39.421000,40.174000,41.656000,42.808000 /)
XPIZA_LKT(63,17,1:6)=(/ 0.737198,0.876517,0.994074,0.994034,0.979026,0.527408 /)
XCGA_LKT(63,17,1:6)=(/ 0.885943,0.853537,0.823310,0.812510,0.793107,0.914287 /)
XEXT_COEFF_550_LKT(63,17)=38.985000 !rg=1.49006 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,18,1:6)=(/ 32.032000,32.295000,32.546000,33.110000,34.144000,35.197000 /)
XPIZA_LKT(63,18,1:6)=(/ 0.724332,0.868087,0.993742,0.993809,0.977582,0.529711 /)
XCGA_LKT(63,18,1:6)=(/ 0.888827,0.856260,0.826213,0.815743,0.802863,0.918937 /)
XEXT_COEFF_550_LKT(63,18)=32.352000 !rg=1.49006 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,19,1:6)=(/ 26.580000,26.742000,27.028000,27.392000,28.201000,28.921000 /)
XPIZA_LKT(63,19,1:6)=(/ 0.710993,0.858547,0.993398,0.993243,0.975786,0.532188 /)
XCGA_LKT(63,19,1:6)=(/ 0.892137,0.858510,0.829027,0.819517,0.808393,0.922767 /)
XEXT_COEFF_550_LKT(63,19)=26.940000 !rg=1.49006 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,20,1:6)=(/ 22.056000,22.188000,22.372000,22.696000,23.320000,23.802000 /)
XPIZA_LKT(63,20,1:6)=(/ 0.697198,0.847759,0.992797,0.992810,0.974239,0.534695 /)
XCGA_LKT(63,20,1:6)=(/ 0.895233,0.861317,0.831163,0.823473,0.811623,0.925830 /)
XEXT_COEFF_550_LKT(63,20)=22.208000 !rg=1.49006 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,1,1:6)=(/ 413.180000,395.190000,402.520000,342.200000,743.320000,621.280000 /)
XPIZA_LKT(64,1,1:6)=(/ 0.866761,0.941772,0.998391,0.998496,0.998267,0.668220 /)
XCGA_LKT(64,1,1:6)=(/ 0.842647,0.767250,0.771937,0.551327,0.767370,0.763347 /)
XEXT_COEFF_550_LKT(64,1)=456.120000 !rg=1.61427 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,2,1:6)=(/ 389.720000,397.640000,403.830000,380.670000,671.750000,590.640000 /)
XPIZA_LKT(64,2,1:6)=(/ 0.857698,0.951358,0.998410,0.998619,0.998079,0.659073 /)
XCGA_LKT(64,2,1:6)=(/ 0.837290,0.791603,0.770650,0.616943,0.753403,0.768910 /)
XEXT_COEFF_550_LKT(64,2)=405.640000 !rg=1.61427 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,3,1:6)=(/ 360.830000,366.890000,378.580000,381.400000,565.630000,531.320000 /)
XPIZA_LKT(64,3,1:6)=(/ 0.849469,0.949268,0.997996,0.998638,0.997704,0.641887 /)
XCGA_LKT(64,3,1:6)=(/ 0.840340,0.799260,0.769023,0.668850,0.731893,0.771027 /)
XEXT_COEFF_550_LKT(64,3)=370.700000 !rg=1.61427 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,4,1:6)=(/ 324.870000,330.250000,341.500000,355.220000,467.840000,461.890000 /)
XPIZA_LKT(64,4,1:6)=(/ 0.839196,0.944770,0.997788,0.998484,0.997209,0.620284 /)
XCGA_LKT(64,4,1:6)=(/ 0.844677,0.803657,0.765833,0.692607,0.714560,0.774577 /)
XEXT_COEFF_550_LKT(64,4)=333.920000 !rg=1.61427 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,5,1:6)=(/ 285.590000,291.690000,302.170000,317.710000,387.620000,393.460000 /)
XPIZA_LKT(64,5,1:6)=(/ 0.828481,0.938245,0.997783,0.998397,0.996596,0.598416 /)
XCGA_LKT(64,5,1:6)=(/ 0.850700,0.805603,0.764527,0.714423,0.705687,0.781633 /)
XEXT_COEFF_550_LKT(64,5)=295.000000 !rg=1.61427 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,6,1:6)=(/ 248.250000,254.260000,260.610000,276.620000,322.480000,331.620000 /)
XPIZA_LKT(64,6,1:6)=(/ 0.818883,0.932131,0.997666,0.998150,0.995952,0.578855 /)
XCGA_LKT(64,6,1:6)=(/ 0.854950,0.811770,0.772570,0.721537,0.706397,0.792330 /)
XEXT_COEFF_550_LKT(64,6)=257.920000 !rg=1.61427 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,7,1:6)=(/ 213.540000,217.460000,222.660000,236.150000,266.740000,276.840000 /)
XPIZA_LKT(64,7,1:6)=(/ 0.810576,0.925489,0.997264,0.997881,0.995257,0.562119 /)
XCGA_LKT(64,7,1:6)=(/ 0.858933,0.818410,0.778273,0.733710,0.712797,0.806023 /)
XEXT_COEFF_550_LKT(64,7)=220.030000 !rg=1.61427 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,8,1:6)=(/ 181.760000,184.730000,191.170000,199.080000,220.880000,230.110000 /)
XPIZA_LKT(64,8,1:6)=(/ 0.803929,0.919690,0.996376,0.997676,0.994327,0.548813 /)
XCGA_LKT(64,8,1:6)=(/ 0.862873,0.823010,0.778250,0.747750,0.717850,0.821083 /)
XEXT_COEFF_550_LKT(64,8)=186.900000 !rg=1.61427 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,9,1:6)=(/ 153.370000,155.900000,159.410000,167.380000,183.070000,190.200000 /)
XPIZA_LKT(64,9,1:6)=(/ 0.798296,0.912991,0.996498,0.997292,0.993079,0.538551 /)
XCGA_LKT(64,9,1:6)=(/ 0.866813,0.827227,0.787483,0.761910,0.723550,0.836667 /)
XEXT_COEFF_550_LKT(64,9)=157.860000 !rg=1.61427 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,10,1:6)=(/ 128.920000,131.520000,134.350000,138.300000,151.790000,156.990000 /)
XPIZA_LKT(64,10,1:6)=(/ 0.793064,0.908243,0.996166,0.996749,0.991471,0.531311 /)
XCGA_LKT(64,10,1:6)=(/ 0.869233,0.831320,0.796200,0.772173,0.735930,0.851513 /)
XEXT_COEFF_550_LKT(64,10)=132.330000 !rg=1.61427 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,11,1:6)=(/ 107.760000,109.760000,111.590000,116.010000,124.190000,128.990000 /)
XPIZA_LKT(64,11,1:6)=(/ 0.787479,0.904281,0.995773,0.995924,0.989796,0.526427 /)
XCGA_LKT(64,11,1:6)=(/ 0.871987,0.836500,0.800093,0.775563,0.743383,0.865483 /)
XEXT_COEFF_550_LKT(64,11)=110.560000 !rg=1.61427 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,12,1:6)=(/ 89.853000,91.319000,93.121000,95.072000,102.260000,105.780000 /)
XPIZA_LKT(64,12,1:6)=(/ 0.781480,0.900215,0.995376,0.995850,0.988176,0.523657 /)
XCGA_LKT(64,12,1:6)=(/ 0.874230,0.839360,0.806390,0.788413,0.756460,0.877990 /)
XEXT_COEFF_550_LKT(64,12)=91.908000 !rg=1.61427 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,13,1:6)=(/ 74.824000,75.749000,77.083000,79.283000,83.668000,86.653000 /)
XPIZA_LKT(64,13,1:6)=(/ 0.773332,0.896679,0.994903,0.995278,0.986276,0.522586 /)
XCGA_LKT(64,13,1:6)=(/ 0.877117,0.843543,0.809300,0.793273,0.762797,0.888877 /)
XEXT_COEFF_550_LKT(64,13)=75.936000 !rg=1.61427 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,14,1:6)=(/ 62.174000,62.939000,63.824000,65.610000,68.625000,71.076000 /)
XPIZA_LKT(64,14,1:6)=(/ 0.765082,0.891761,0.994817,0.994882,0.983860,0.522855 /)
XCGA_LKT(64,14,1:6)=(/ 0.879153,0.846143,0.814683,0.797327,0.775387,0.897993 /)
XEXT_COEFF_550_LKT(64,14)=63.095000 !rg=1.61427 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,15,1:6)=(/ 51.726000,52.230000,52.941000,54.374000,56.378000,58.315000 /)
XPIZA_LKT(64,15,1:6)=(/ 0.755219,0.886690,0.994490,0.994270,0.981741,0.524059 /)
XCGA_LKT(64,15,1:6)=(/ 0.882103,0.848627,0.818943,0.802033,0.785923,0.905607 /)
XEXT_COEFF_550_LKT(64,15)=52.405000 !rg=1.61427 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,16,1:6)=(/ 42.809000,43.297000,43.709000,44.879000,46.325000,47.780000 /)
XPIZA_LKT(64,16,1:6)=(/ 0.743748,0.880485,0.994207,0.993850,0.979732,0.525912 /)
XCGA_LKT(64,16,1:6)=(/ 0.884653,0.852777,0.821123,0.808040,0.790173,0.911990 /)
XEXT_COEFF_550_LKT(64,16)=43.465000 !rg=1.61427 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,17,1:6)=(/ 35.541000,35.801000,36.296000,36.970000,38.174000,39.261000 /)
XPIZA_LKT(64,17,1:6)=(/ 0.731541,0.872983,0.993928,0.993780,0.978209,0.528172 /)
XCGA_LKT(64,17,1:6)=(/ 0.887447,0.854483,0.824083,0.813257,0.798360,0.917087 /)
XEXT_COEFF_550_LKT(64,17)=35.897000 !rg=1.61427 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,18,1:6)=(/ 29.554000,29.753000,30.149000,30.649000,31.424000,32.299000 /)
XPIZA_LKT(64,18,1:6)=(/ 0.718742,0.864224,0.993569,0.993335,0.976635,0.530641 /)
XCGA_LKT(64,18,1:6)=(/ 0.890413,0.856857,0.827783,0.816590,0.806763,0.921207 /)
XEXT_COEFF_550_LKT(64,18)=29.872000 !rg=1.61427 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,19,1:6)=(/ 24.487000,24.674000,24.867000,25.342000,25.851000,26.554000 /)
XPIZA_LKT(64,19,1:6)=(/ 0.704771,0.854031,0.992958,0.993067,0.975141,0.533208 /)
XCGA_LKT(64,19,1:6)=(/ 0.893510,0.860207,0.829063,0.821007,0.810270,0.924587 /)
XEXT_COEFF_550_LKT(64,19)=24.745000 !rg=1.61427 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,20,1:6)=(/ 20.318000,20.433000,20.629000,20.899000,21.385000,21.865000 /)
XPIZA_LKT(64,20,1:6)=(/ 0.690744,0.842344,0.992439,0.992716,0.973745,0.535745 /)
XCGA_LKT(64,20,1:6)=(/ 0.896813,0.862287,0.831573,0.824120,0.816120,0.927290 /)
XEXT_COEFF_550_LKT(64,20)=20.473000 !rg=1.61427 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,1,1:6)=(/ 372.210000,367.790000,339.390000,351.330000,652.280000,588.100000 /)
XPIZA_LKT(65,1,1:6)=(/ 0.852539,0.950221,0.997899,0.998527,0.998020,0.658380 /)
XCGA_LKT(65,1,1:6)=(/ 0.839070,0.788897,0.711410,0.619160,0.752150,0.778903 /)
XEXT_COEFF_550_LKT(65,1)=410.890000 !rg=1.74883 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,2,1:6)=(/ 357.160000,366.030000,361.370000,371.060000,579.970000,548.190000 /)
XPIZA_LKT(65,2,1:6)=(/ 0.849457,0.948470,0.998253,0.998600,0.997761,0.648308 /)
XCGA_LKT(65,2,1:6)=(/ 0.843210,0.796380,0.768040,0.665650,0.734377,0.779270 /)
XEXT_COEFF_550_LKT(65,2)=375.730000 !rg=1.74883 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,3,1:6)=(/ 330.520000,337.980000,352.000000,360.020000,487.590000,488.240000 /)
XPIZA_LKT(65,3,1:6)=(/ 0.840784,0.945757,0.997918,0.998471,0.997311,0.628661 /)
XCGA_LKT(65,3,1:6)=(/ 0.845743,0.801680,0.768503,0.692743,0.711603,0.779903 /)
XEXT_COEFF_550_LKT(65,3)=344.840000 !rg=1.74883 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,4,1:6)=(/ 297.710000,304.820000,316.910000,329.340000,408.720000,421.860000 /)
XPIZA_LKT(65,4,1:6)=(/ 0.830299,0.940541,0.997820,0.998289,0.996768,0.606164 /)
XCGA_LKT(65,4,1:6)=(/ 0.850057,0.805523,0.770060,0.708640,0.698943,0.783010 /)
XEXT_COEFF_550_LKT(65,4)=310.240000 !rg=1.74883 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,5,1:6)=(/ 263.260000,267.650000,275.820000,289.900000,343.560000,358.520000 /)
XPIZA_LKT(65,5,1:6)=(/ 0.819787,0.935122,0.996736,0.998124,0.996207,0.585070 /)
XCGA_LKT(65,5,1:6)=(/ 0.854263,0.812477,0.772467,0.717993,0.700390,0.790530 /)
XEXT_COEFF_550_LKT(65,5)=270.900000 !rg=1.74883 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,6,1:6)=(/ 228.490000,233.760000,241.570000,249.740000,290.070000,302.040000 /)
XPIZA_LKT(65,6,1:6)=(/ 0.811773,0.927918,0.997243,0.998114,0.995411,0.567198 /)
XCGA_LKT(65,6,1:6)=(/ 0.858723,0.814513,0.776943,0.741743,0.702207,0.801840 /)
XEXT_COEFF_550_LKT(65,6)=235.660000 !rg=1.74883 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,7,1:6)=(/ 196.360000,200.590000,205.720000,216.660000,241.460000,252.190000 /)
XPIZA_LKT(65,7,1:6)=(/ 0.804789,0.920707,0.997090,0.997453,0.994568,0.552456 /)
XCGA_LKT(65,7,1:6)=(/ 0.862753,0.820040,0.783980,0.743957,0.710973,0.815850 /)
XEXT_COEFF_550_LKT(65,7)=203.090000 !rg=1.74883 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,8,1:6)=(/ 166.790000,169.990000,175.100000,180.200000,201.630000,209.740000 /)
XPIZA_LKT(65,8,1:6)=(/ 0.799896,0.915482,0.996547,0.997567,0.993478,0.541154 /)
XCGA_LKT(65,8,1:6)=(/ 0.865523,0.825400,0.787633,0.761770,0.720360,0.830800 /)
XEXT_COEFF_550_LKT(65,8)=172.440000 !rg=1.74883 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,9,1:6)=(/ 141.170000,143.530000,145.670000,153.150000,165.040000,173.480000 /)
XPIZA_LKT(65,9,1:6)=(/ 0.794596,0.910492,0.996211,0.996996,0.992751,0.532747 /)
XCGA_LKT(65,9,1:6)=(/ 0.868543,0.831613,0.794973,0.761863,0.734657,0.845887 /)
XEXT_COEFF_550_LKT(65,9)=144.930000 !rg=1.74883 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,10,1:6)=(/ 118.850000,120.380000,122.740000,128.130000,137.580000,143.290000 /)
XPIZA_LKT(65,10,1:6)=(/ 0.789704,0.906337,0.995746,0.996352,0.991148,0.527119 /)
XCGA_LKT(65,10,1:6)=(/ 0.871563,0.834653,0.796837,0.773500,0.736537,0.859963 /)
XEXT_COEFF_550_LKT(65,10)=121.670000 !rg=1.74883 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,11,1:6)=(/ 99.222000,100.680000,102.980000,106.160000,113.510000,117.820000 /)
XPIZA_LKT(65,11,1:6)=(/ 0.783296,0.902172,0.995465,0.995850,0.989586,0.523617 /)
XCGA_LKT(65,11,1:6)=(/ 0.874120,0.838897,0.802037,0.780530,0.751687,0.872990 /)
XEXT_COEFF_550_LKT(65,11)=101.570000 !rg=1.74883 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,12,1:6)=(/ 83.003000,83.890000,85.135000,88.028000,93.103000,96.700000 /)
XPIZA_LKT(65,12,1:6)=(/ 0.777445,0.898693,0.995310,0.995377,0.987697,0.521967 /)
XCGA_LKT(65,12,1:6)=(/ 0.876313,0.842340,0.808330,0.790243,0.757657,0.884493 /)
XEXT_COEFF_550_LKT(65,12)=84.371000 !rg=1.74883 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,13,1:6)=(/ 68.829000,69.708000,70.837000,73.210000,76.448000,79.273000 /)
XPIZA_LKT(65,13,1:6)=(/ 0.769558,0.894175,0.994994,0.994619,0.984978,0.521783 /)
XCGA_LKT(65,13,1:6)=(/ 0.878293,0.844613,0.812703,0.791223,0.770730,0.894390 /)
XEXT_COEFF_550_LKT(65,13)=69.749000 !rg=1.74883 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,14,1:6)=(/ 57.323000,57.782000,58.788000,60.288000,62.827000,65.072000 /)
XPIZA_LKT(65,14,1:6)=(/ 0.760922,0.889464,0.994617,0.994837,0.982914,0.522713 /)
XCGA_LKT(65,14,1:6)=(/ 0.880900,0.847967,0.816933,0.802573,0.781317,0.902603 /)
XEXT_COEFF_550_LKT(65,14)=58.313000 !rg=1.74883 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,15,1:6)=(/ 47.620000,48.048000,48.599000,49.631000,51.731000,53.427000 /)
XPIZA_LKT(65,15,1:6)=(/ 0.750024,0.883902,0.994243,0.994349,0.980715,0.524397 /)
XCGA_LKT(65,15,1:6)=(/ 0.883423,0.850800,0.819437,0.805123,0.787610,0.909410 /)
XEXT_COEFF_550_LKT(65,15)=48.317000 !rg=1.74883 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,16,1:6)=(/ 39.444000,39.783000,40.336000,41.115000,42.592000,43.804000 /)
XPIZA_LKT(65,16,1:6)=(/ 0.737887,0.877188,0.994099,0.993948,0.979273,0.526591 /)
XCGA_LKT(65,16,1:6)=(/ 0.886127,0.853767,0.823110,0.810510,0.795490,0.915083 /)
XEXT_COEFF_550_LKT(65,16)=40.046000 !rg=1.74883 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,17,1:6)=(/ 32.775000,32.985000,33.442000,34.144000,35.061000,36.016000 /)
XPIZA_LKT(65,17,1:6)=(/ 0.725920,0.869243,0.993734,0.993369,0.977194,0.529062 /)
XCGA_LKT(65,17,1:6)=(/ 0.888973,0.856240,0.826990,0.815360,0.800830,0.919597 /)
XEXT_COEFF_550_LKT(65,17)=33.227000 !rg=1.74883 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,18,1:6)=(/ 27.216000,27.435000,27.602000,28.083000,28.919000,29.647000 /)
XPIZA_LKT(65,18,1:6)=(/ 0.712344,0.860163,0.993379,0.993292,0.975834,0.531652 /)
XCGA_LKT(65,18,1:6)=(/ 0.891767,0.858697,0.828867,0.819107,0.807913,0.923230 /)
XEXT_COEFF_550_LKT(65,18)=27.482000 !rg=1.74883 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,19,1:6)=(/ 22.570000,22.710000,22.960000,23.247000,23.845000,24.387000 /)
XPIZA_LKT(65,19,1:6)=(/ 0.698346,0.849078,0.992808,0.992888,0.974290,0.534275 /)
XCGA_LKT(65,19,1:6)=(/ 0.895013,0.861163,0.831013,0.821990,0.813020,0.926207 /)
XEXT_COEFF_550_LKT(65,19)=22.812000 !rg=1.74883 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,20,1:6)=(/ 18.747000,18.828000,19.020000,19.316000,19.680000,20.090000 /)
XPIZA_LKT(65,20,1:6)=(/ 0.684758,0.836891,0.992220,0.992377,0.972738,0.536814 /)
XCGA_LKT(65,20,1:6)=(/ 0.898457,0.863693,0.833827,0.826130,0.817383,0.928583 /)
XEXT_COEFF_550_LKT(65,20)=18.947000 !rg=1.74883 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,1,1:6)=(/ 340.300000,366.110000,347.750000,379.600000,546.460000,538.770000 /)
XPIZA_LKT(66,1,1:6)=(/ 0.842281,0.950247,0.997835,0.998638,0.997611,0.647877 /)
XCGA_LKT(66,1,1:6)=(/ 0.843810,0.815640,0.767237,0.697650,0.723130,0.789687 /)
XEXT_COEFF_550_LKT(66,1)=337.040000 !rg=1.89461 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,2,1:6)=(/ 328.120000,339.050000,344.070000,366.300000,486.530000,500.460000 /)
XPIZA_LKT(66,2,1:6)=(/ 0.839785,0.942741,0.997971,0.998612,0.997301,0.633056 /)
XCGA_LKT(66,2,1:6)=(/ 0.847457,0.800910,0.775770,0.713747,0.706470,0.787393 /)
XEXT_COEFF_550_LKT(66,2)=338.630000 !rg=1.89461 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,3,1:6)=(/ 302.230000,312.220000,320.800000,334.800000,413.560000,442.990000 /)
XPIZA_LKT(66,3,1:6)=(/ 0.830802,0.940149,0.997737,0.998431,0.996832,0.612004 /)
XCGA_LKT(66,3,1:6)=(/ 0.850930,0.806620,0.761200,0.716770,0.691283,0.787260 /)
XEXT_COEFF_550_LKT(66,3)=314.580000 !rg=1.89461 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,4,1:6)=(/ 272.630000,280.600000,289.870000,300.330000,355.310000,382.120000 /)
XPIZA_LKT(66,4,1:6)=(/ 0.820587,0.935612,0.997696,0.998132,0.996287,0.590088 /)
XCGA_LKT(66,4,1:6)=(/ 0.855140,0.810730,0.765570,0.724130,0.688250,0.790993 /)
XEXT_COEFF_550_LKT(66,4)=282.490000 !rg=1.89461 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,5,1:6)=(/ 241.680000,247.120000,255.750000,268.040000,304.900000,325.060000 /)
XPIZA_LKT(66,5,1:6)=(/ 0.811419,0.930341,0.997385,0.997941,0.995678,0.571001 /)
XCGA_LKT(66,5,1:6)=(/ 0.858850,0.814337,0.775750,0.730650,0.692927,0.799587 /)
XEXT_COEFF_550_LKT(66,5)=250.560000 !rg=1.89461 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,6,1:6)=(/ 209.860000,214.810000,220.440000,229.070000,259.120000,274.290000 /)
XPIZA_LKT(66,6,1:6)=(/ 0.804855,0.923184,0.997285,0.997762,0.994753,0.555585 /)
XCGA_LKT(66,6,1:6)=(/ 0.861697,0.819537,0.778367,0.747840,0.701087,0.811647 /)
XEXT_COEFF_550_LKT(66,6)=215.690000 !rg=1.89461 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,7,1:6)=(/ 179.940000,184.170000,188.900000,196.720000,217.470000,229.350000 /)
XPIZA_LKT(66,7,1:6)=(/ 0.799836,0.916900,0.996794,0.997234,0.993864,0.543259 /)
XCGA_LKT(66,7,1:6)=(/ 0.865537,0.824840,0.783723,0.752553,0.711397,0.825857 /)
XEXT_COEFF_550_LKT(66,7)=185.530000 !rg=1.89461 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,8,1:6)=(/ 153.950000,156.270000,159.700000,166.350000,182.310000,190.990000 /)
XPIZA_LKT(66,8,1:6)=(/ 0.794799,0.912210,0.996436,0.997294,0.992966,0.534141 /)
XCGA_LKT(66,8,1:6)=(/ 0.868367,0.829203,0.790380,0.765090,0.719693,0.840500 /)
XEXT_COEFF_550_LKT(66,8)=157.320000 !rg=1.89461 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,9,1:6)=(/ 130.040000,132.450000,135.080000,140.940000,150.560000,158.150000 /)
XPIZA_LKT(66,9,1:6)=(/ 0.790726,0.907604,0.995718,0.996643,0.992000,0.527637 /)
XCGA_LKT(66,9,1:6)=(/ 0.870877,0.832703,0.795717,0.767773,0.737420,0.854900 /)
XEXT_COEFF_550_LKT(66,9)=132.910000 !rg=1.89461 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,10,1:6)=(/ 109.410000,111.260000,112.430000,117.480000,124.710000,130.760000 /)
XPIZA_LKT(66,10,1:6)=(/ 0.786786,0.903024,0.995620,0.995709,0.990384,0.523588 /)
XCGA_LKT(66,10,1:6)=(/ 0.872283,0.837263,0.803247,0.774900,0.748390,0.868070 /)
XEXT_COEFF_550_LKT(66,10)=111.380000 !rg=1.89461 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,11,1:6)=(/ 91.476000,92.723000,94.434000,97.357000,103.110000,107.620000 /)
XPIZA_LKT(66,11,1:6)=(/ 0.780690,0.899798,0.995192,0.995817,0.988836,0.521384 /)
XCGA_LKT(66,11,1:6)=(/ 0.875320,0.841207,0.803100,0.787733,0.753913,0.880073 /)
XEXT_COEFF_550_LKT(66,11)=93.593000 !rg=1.89461 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,12,1:6)=(/ 76.296000,77.327000,78.272000,80.877000,84.933000,88.409000 /)
XPIZA_LKT(66,12,1:6)=(/ 0.774008,0.895831,0.995044,0.995169,0.986202,0.520766 /)
XCGA_LKT(66,12,1:6)=(/ 0.876990,0.843710,0.811317,0.790750,0.766570,0.890540 /)
XEXT_COEFF_550_LKT(66,12)=77.504000 !rg=1.89461 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,13,1:6)=(/ 63.502000,64.081000,65.200000,66.973000,69.870000,72.538000 /)
XPIZA_LKT(66,13,1:6)=(/ 0.766071,0.892219,0.994718,0.994757,0.984038,0.521371 /)
XCGA_LKT(66,13,1:6)=(/ 0.879843,0.847177,0.815660,0.799457,0.776187,0.899457 /)
XEXT_COEFF_550_LKT(66,13)=64.569000 !rg=1.89461 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,14,1:6)=(/ 52.803000,53.294000,54.094000,55.138000,57.593000,59.590000 /)
XPIZA_LKT(66,14,1:6)=(/ 0.755959,0.887358,0.994561,0.994401,0.982378,0.522871 /)
XCGA_LKT(66,14,1:6)=(/ 0.881870,0.849987,0.817637,0.803353,0.783553,0.906793 /)
XEXT_COEFF_550_LKT(66,14)=53.381000 !rg=1.89461 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,15,1:6)=(/ 43.887000,44.358000,44.926000,45.520000,47.639000,48.961000 /)
XPIZA_LKT(66,15,1:6)=(/ 0.744920,0.880831,0.994339,0.994019,0.979713,0.524961 /)
XCGA_LKT(66,15,1:6)=(/ 0.884750,0.852200,0.821693,0.810033,0.791370,0.912840 /)
XEXT_COEFF_550_LKT(66,15)=44.559000 !rg=1.89461 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,16,1:6)=(/ 36.419000,36.683000,37.120000,37.786000,38.811000,40.170000 /)
XPIZA_LKT(66,16,1:6)=(/ 0.733059,0.873574,0.993882,0.993879,0.978733,0.527426 /)
XCGA_LKT(66,16,1:6)=(/ 0.887463,0.855107,0.823463,0.814170,0.798277,0.917860 /)
XEXT_COEFF_550_LKT(66,16)=36.892000 !rg=1.89461 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,17,1:6)=(/ 30.199000,30.480000,30.849000,31.276000,32.134000,33.048000 /)
XPIZA_LKT(66,17,1:6)=(/ 0.719682,0.865436,0.993520,0.993369,0.976580,0.530057 /)
XCGA_LKT(66,17,1:6)=(/ 0.890183,0.857887,0.826263,0.817203,0.804453,0.921833 /)
XEXT_COEFF_550_LKT(66,17)=30.497000 !rg=1.89461 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,18,1:6)=(/ 25.092000,25.279000,25.515000,25.759000,26.679000,27.219000 /)
XPIZA_LKT(66,18,1:6)=(/ 0.706241,0.855382,0.993172,0.993167,0.974767,0.532727 /)
XCGA_LKT(66,18,1:6)=(/ 0.893187,0.859537,0.829947,0.822203,0.810460,0.925027 /)
XEXT_COEFF_550_LKT(66,18)=25.401000 !rg=1.89461 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,19,1:6)=(/ 20.840000,20.946000,21.135000,21.413000,21.789000,22.401000 /)
XPIZA_LKT(66,19,1:6)=(/ 0.692631,0.843806,0.992544,0.992721,0.974041,0.535372 /)
XCGA_LKT(66,19,1:6)=(/ 0.896563,0.862410,0.831467,0.824600,0.815853,0.927640 /)
XEXT_COEFF_550_LKT(66,19)=21.031000 !rg=1.89461 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,20,1:6)=(/ 17.282000,17.402000,17.526000,17.711000,18.056000,18.463000 /)
XPIZA_LKT(66,20,1:6)=(/ 0.678217,0.831383,0.991646,0.992099,0.972107,0.537889 /)
XCGA_LKT(66,20,1:6)=(/ 0.899953,0.865307,0.833143,0.827393,0.819640,0.929723 /)
XEXT_COEFF_550_LKT(66,20)=17.402000 !rg=1.89461 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,1,1:6)=(/ 317.030000,318.460000,371.700000,400.560000,441.820000,483.800000 /)
XPIZA_LKT(67,1,1:6)=(/ 0.835803,0.943159,0.996958,0.998698,0.997025,0.626757 /)
XCGA_LKT(67,1,1:6)=(/ 0.852807,0.809370,0.769277,0.751187,0.687033,0.794653 /)
XEXT_COEFF_550_LKT(67,1)=309.020000 !rg=2.05254 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,2,1:6)=(/ 301.670000,307.330000,322.060000,353.510000,395.790000,449.780000 /)
XPIZA_LKT(67,2,1:6)=(/ 0.830146,0.941634,0.997971,0.998554,0.996694,0.613133 /)
XCGA_LKT(67,2,1:6)=(/ 0.851520,0.810767,0.764270,0.745020,0.672910,0.793350 /)
XEXT_COEFF_550_LKT(67,2)=308.720000 !rg=2.05254 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,3,1:6)=(/ 278.020000,283.260000,295.690000,308.970000,352.390000,397.650000 /)
XPIZA_LKT(67,3,1:6)=(/ 0.820030,0.938595,0.996622,0.998395,0.996207,0.592307 /)
XCGA_LKT(67,3,1:6)=(/ 0.856063,0.813313,0.771967,0.737267,0.666757,0.793677 /)
XEXT_COEFF_550_LKT(67,3)=288.220000 !rg=2.05254 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,4,1:6)=(/ 250.690000,255.370000,267.150000,274.010000,312.700000,343.980000 /)
XPIZA_LKT(67,4,1:6)=(/ 0.810843,0.932730,0.996449,0.998213,0.995674,0.572731 /)
XCGA_LKT(67,4,1:6)=(/ 0.859837,0.815850,0.771973,0.739863,0.677107,0.799040 /)
XEXT_COEFF_550_LKT(67,4)=259.400000 !rg=2.05254 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,5,1:6)=(/ 221.740000,227.250000,234.590000,243.440000,271.060000,293.720000 /)
XPIZA_LKT(67,5,1:6)=(/ 0.804175,0.924746,0.997118,0.997683,0.995101,0.556816 /)
XCGA_LKT(67,5,1:6)=(/ 0.863270,0.818547,0.774897,0.740940,0.691413,0.809073 /)
XEXT_COEFF_550_LKT(67,5)=228.960000 !rg=2.05254 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,6,1:6)=(/ 192.910000,196.270000,202.340000,209.460000,234.180000,248.630000 /)
XPIZA_LKT(67,6,1:6)=(/ 0.799454,0.919457,0.996606,0.997730,0.994034,0.544452 /)
XCGA_LKT(67,6,1:6)=(/ 0.865807,0.822163,0.783280,0.751437,0.704460,0.821803 /)
XEXT_COEFF_550_LKT(67,6)=197.130000 !rg=2.05254 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,7,1:6)=(/ 165.670000,168.590000,173.870000,179.390000,198.000000,208.380000 /)
XPIZA_LKT(67,7,1:6)=(/ 0.794329,0.913515,0.996517,0.997290,0.993253,0.534807 /)
XCGA_LKT(67,7,1:6)=(/ 0.868787,0.828367,0.787697,0.759260,0.716710,0.835983 /)
XEXT_COEFF_550_LKT(67,7)=170.520000 !rg=2.05254 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,8,1:6)=(/ 141.430000,143.750000,146.430000,152.990000,165.240000,173.840000 /)
XPIZA_LKT(67,8,1:6)=(/ 0.792301,0.908170,0.996083,0.996861,0.992136,0.527934 /)
XCGA_LKT(67,8,1:6)=(/ 0.869920,0.831020,0.795167,0.764277,0.730193,0.850080 /)
XEXT_COEFF_550_LKT(67,8)=144.140000 !rg=2.05254 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,9,1:6)=(/ 119.380000,121.500000,123.950000,128.480000,137.400000,144.160000 /)
XPIZA_LKT(67,9,1:6)=(/ 0.787797,0.904040,0.995671,0.996267,0.990271,0.523292 /)
XCGA_LKT(67,9,1:6)=(/ 0.872753,0.835747,0.798140,0.771473,0.737850,0.863610 /)
XEXT_COEFF_550_LKT(67,9)=122.730000 !rg=2.05254 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,10,1:6)=(/ 100.600000,102.150000,103.950000,107.820000,113.620000,119.340000 /)
XPIZA_LKT(67,10,1:6)=(/ 0.782965,0.901532,0.995525,0.995906,0.989768,0.520735 /)
XCGA_LKT(67,10,1:6)=(/ 0.874413,0.838847,0.805713,0.782143,0.754143,0.875757 /)
XEXT_COEFF_550_LKT(67,10)=102.410000 !rg=2.05254 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,11,1:6)=(/ 84.178000,85.300000,86.486000,89.573000,93.479000,98.328000 /)
XPIZA_LKT(67,11,1:6)=(/ 0.777133,0.897805,0.995298,0.995490,0.988181,0.519725 /)
XCGA_LKT(67,11,1:6)=(/ 0.876510,0.843597,0.809727,0.787187,0.764980,0.886683 /)
XEXT_COEFF_550_LKT(67,11)=86.274000 !rg=2.05254 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,12,1:6)=(/ 70.258000,70.991000,72.336000,74.303000,77.552000,80.850000 /)
XPIZA_LKT(67,12,1:6)=(/ 0.770392,0.894030,0.994741,0.995205,0.985572,0.520027 /)
XCGA_LKT(67,12,1:6)=(/ 0.878707,0.845520,0.812750,0.797023,0.773387,0.896110 /)
XEXT_COEFF_550_LKT(67,12)=71.406000 !rg=2.05254 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,13,1:6)=(/ 58.456000,59.042000,60.021000,61.201000,63.985000,66.394000 /)
XPIZA_LKT(67,13,1:6)=(/ 0.760974,0.889865,0.994541,0.994747,0.983459,0.521316 /)
XCGA_LKT(67,13,1:6)=(/ 0.880823,0.849190,0.816437,0.801857,0.779423,0.904073 /)
XEXT_COEFF_550_LKT(67,13)=59.159000 !rg=2.05254 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,14,1:6)=(/ 48.618000,49.051000,49.913000,50.606000,52.957000,54.585000 /)
XPIZA_LKT(67,14,1:6)=(/ 0.751126,0.884256,0.994463,0.994317,0.980988,0.523301 /)
XCGA_LKT(67,14,1:6)=(/ 0.883247,0.850920,0.820637,0.808410,0.787493,0.910580 /)
XEXT_COEFF_550_LKT(67,14)=49.413000 !rg=2.05254 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,15,1:6)=(/ 40.467000,40.720000,41.125000,42.128000,43.497000,44.881000 /)
XPIZA_LKT(67,15,1:6)=(/ 0.739866,0.878106,0.994108,0.993967,0.979139,0.525715 /)
XCGA_LKT(67,15,1:6)=(/ 0.886163,0.853417,0.823513,0.811563,0.794040,0.915917 /)
XEXT_COEFF_550_LKT(67,15)=41.033000 !rg=2.05254 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,16,1:6)=(/ 33.530000,33.795000,34.067000,34.856000,35.569000,36.847000 /)
XPIZA_LKT(67,16,1:6)=(/ 0.727047,0.870164,0.993793,0.993626,0.977730,0.528392 /)
XCGA_LKT(67,16,1:6)=(/ 0.888713,0.856357,0.826720,0.813753,0.804260,0.920333 /)
XEXT_COEFF_550_LKT(67,16)=33.972000 !rg=2.05254 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,17,1:6)=(/ 27.857000,28.007000,28.393000,28.762000,29.484000,30.333000 /)
XPIZA_LKT(67,17,1:6)=(/ 0.713868,0.861135,0.993304,0.993440,0.976369,0.531136 /)
XCGA_LKT(67,17,1:6)=(/ 0.891637,0.858563,0.828013,0.820157,0.809017,0.923820 /)
XEXT_COEFF_550_LKT(67,17)=28.165000 !rg=2.05254 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,18,1:6)=(/ 23.158000,23.252000,23.444000,23.852000,24.392000,24.996000 /)
XPIZA_LKT(67,18,1:6)=(/ 0.700343,0.850551,0.992873,0.992809,0.974427,0.533847 /)
XCGA_LKT(67,18,1:6)=(/ 0.894817,0.860810,0.830883,0.823390,0.813163,0.926617 /)
XEXT_COEFF_550_LKT(67,18)=23.390000 !rg=2.05254 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,19,1:6)=(/ 19.199000,19.316000,19.424000,19.770000,20.022000,20.582000 /)
XPIZA_LKT(67,19,1:6)=(/ 0.686000,0.838605,0.992300,0.992393,0.973019,0.536486 /)
XCGA_LKT(67,19,1:6)=(/ 0.898080,0.863640,0.833013,0.824150,0.819283,0.928907 /)
XEXT_COEFF_550_LKT(67,19)=19.376000 !rg=2.05254 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,20,1:6)=(/ 15.944000,16.007000,16.173000,16.332000,16.625000,16.971000 /)
XPIZA_LKT(67,20,1:6)=(/ 0.672213,0.825113,0.991403,0.992012,0.971787,0.538963 /)
XCGA_LKT(67,20,1:6)=(/ 0.901580,0.866117,0.834640,0.829000,0.823090,0.930727 /)
XEXT_COEFF_550_LKT(67,20)=16.069000 !rg=2.05254 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,1,1:6)=(/ 290.610000,281.080000,332.900000,377.670000,345.710000,431.740000 /)
XPIZA_LKT(68,1,1:6)=(/ 0.823180,0.938619,0.997776,0.998644,0.996182,0.602739 /)
XCGA_LKT(68,1,1:6)=(/ 0.855533,0.808063,0.780657,0.787457,0.636460,0.798993 /)
XEXT_COEFF_550_LKT(68,1)=313.510000 !rg=2.22363 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,2,1:6)=(/ 275.930000,281.830000,295.740000,326.510000,323.490000,398.830000 /)
XPIZA_LKT(68,2,1:6)=(/ 0.819720,0.935937,0.997650,0.998453,0.995869,0.588796 /)
XCGA_LKT(68,2,1:6)=(/ 0.857723,0.810573,0.770780,0.754193,0.632797,0.797657 /)
XEXT_COEFF_550_LKT(68,2)=283.540000 !rg=2.22363 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,3,1:6)=(/ 256.660000,261.360000,268.530000,287.910000,298.540000,354.190000 /)
XPIZA_LKT(68,3,1:6)=(/ 0.810230,0.931611,0.997663,0.998175,0.995567,0.570434 /)
XCGA_LKT(68,3,1:6)=(/ 0.861693,0.812713,0.777920,0.749247,0.651113,0.799927 /)
XEXT_COEFF_550_LKT(68,3)=263.890000 !rg=2.22363 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,4,1:6)=(/ 231.410000,235.310000,242.010000,254.760000,273.930000,308.490000 /)
XPIZA_LKT(68,4,1:6)=(/ 0.802277,0.926331,0.997381,0.998033,0.995119,0.555012 /)
XCGA_LKT(68,4,1:6)=(/ 0.864763,0.817270,0.779340,0.750987,0.673657,0.807687 /)
XEXT_COEFF_550_LKT(68,4)=237.770000 !rg=2.22363 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,5,1:6)=(/ 203.940000,207.550000,214.770000,220.690000,244.720000,264.940000 /)
XPIZA_LKT(68,5,1:6)=(/ 0.797701,0.921081,0.996413,0.997828,0.994283,0.543183 /)
XCGA_LKT(68,5,1:6)=(/ 0.866763,0.823840,0.783213,0.755277,0.691827,0.819170 /)
XEXT_COEFF_550_LKT(68,5)=209.570000 !rg=2.22363 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,6,1:6)=(/ 177.390000,179.930000,185.070000,194.350000,209.770000,225.180000 /)
XPIZA_LKT(68,6,1:6)=(/ 0.794775,0.915160,0.996696,0.997382,0.993328,0.534209 /)
XCGA_LKT(68,6,1:6)=(/ 0.868800,0.827927,0.791983,0.758567,0.708543,0.832287 /)
XEXT_COEFF_550_LKT(68,6)=182.470000 !rg=2.22363 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,7,1:6)=(/ 152.740000,154.940000,159.180000,165.270000,178.400000,189.270000 /)
XPIZA_LKT(68,7,1:6)=(/ 0.791543,0.909520,0.996205,0.997200,0.992844,0.527319 /)
XCGA_LKT(68,7,1:6)=(/ 0.870827,0.831290,0.788927,0.768800,0.719593,0.846120 /)
XEXT_COEFF_550_LKT(68,7)=156.960000 !rg=2.22363 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,8,1:6)=(/ 130.070000,131.930000,135.030000,140.170000,149.380000,158.230000 /)
XPIZA_LKT(68,8,1:6)=(/ 0.789315,0.905825,0.995862,0.996815,0.991808,0.522640 /)
XCGA_LKT(68,8,1:6)=(/ 0.872167,0.834910,0.799857,0.774950,0.736853,0.859417 /)
XEXT_COEFF_550_LKT(68,8)=133.250000 !rg=2.22363 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,9,1:6)=(/ 110.040000,111.740000,114.270000,116.890000,126.050000,131.420000 /)
XPIZA_LKT(68,9,1:6)=(/ 0.785254,0.902029,0.995431,0.996228,0.989809,0.519758 /)
XCGA_LKT(68,9,1:6)=(/ 0.874273,0.838920,0.801543,0.781453,0.745147,0.871913 /)
XEXT_COEFF_550_LKT(68,9)=112.180000 !rg=2.22363 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,10,1:6)=(/ 92.749000,93.845000,95.904000,98.072000,103.890000,108.940000 /)
XPIZA_LKT(68,10,1:6)=(/ 0.780076,0.899650,0.995195,0.995758,0.988816,0.518558 /)
XCGA_LKT(68,10,1:6)=(/ 0.875940,0.841783,0.805540,0.786310,0.756003,0.882957 /)
XEXT_COEFF_550_LKT(68,10)=94.171000 !rg=2.22363 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,11,1:6)=(/ 77.619000,78.543000,79.710000,82.196000,86.156000,89.859000 /)
XPIZA_LKT(68,11,1:6)=(/ 0.773664,0.895690,0.994759,0.995208,0.986898,0.518616 /)
XCGA_LKT(68,11,1:6)=(/ 0.878160,0.844733,0.811600,0.792577,0.767933,0.892783 /)
XEXT_COEFF_550_LKT(68,11)=79.066000 !rg=2.22363 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,12,1:6)=(/ 64.765000,65.381000,66.579000,67.839000,71.108000,73.960000 /)
XPIZA_LKT(68,12,1:6)=(/ 0.766315,0.892443,0.994710,0.994837,0.984613,0.519715 /)
XCGA_LKT(68,12,1:6)=(/ 0.879710,0.847663,0.813363,0.798263,0.775460,0.901190 /)
XEXT_COEFF_550_LKT(68,12)=65.493000 !rg=2.22363 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,13,1:6)=(/ 53.916000,54.289000,55.470000,56.240000,58.814000,60.788000 /)
XPIZA_LKT(68,13,1:6)=(/ 0.757313,0.887309,0.994411,0.994519,0.982076,0.521587 /)
XCGA_LKT(68,13,1:6)=(/ 0.882313,0.850140,0.818207,0.805770,0.784280,0.908243 /)
XEXT_COEFF_550_LKT(68,13)=54.715000 !rg=2.22363 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,14,1:6)=(/ 44.867000,45.302000,45.811000,46.739000,48.595000,50.016000 /)
XPIZA_LKT(68,14,1:6)=(/ 0.746099,0.881804,0.994190,0.994190,0.980082,0.523964 /)
XCGA_LKT(68,14,1:6)=(/ 0.884797,0.852757,0.822217,0.810360,0.789287,0.913977 /)
XEXT_COEFF_550_LKT(68,14)=45.339000 !rg=2.22363 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,15,1:6)=(/ 37.269000,37.658000,37.856000,38.558000,39.843000,41.153000 /)
XPIZA_LKT(68,15,1:6)=(/ 0.734109,0.874588,0.994013,0.993850,0.978543,0.526635 /)
XCGA_LKT(68,15,1:6)=(/ 0.887040,0.855137,0.825100,0.813957,0.799613,0.918663 /)
XEXT_COEFF_550_LKT(68,15)=37.660000 !rg=2.22363 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,16,1:6)=(/ 30.943000,31.133000,31.464000,32.005000,32.903000,33.809000 /)
XPIZA_LKT(68,16,1:6)=(/ 0.721368,0.866220,0.993584,0.993555,0.976880,0.529466 /)
XCGA_LKT(68,16,1:6)=(/ 0.890280,0.857180,0.828127,0.818307,0.806317,0.922530 /)
XEXT_COEFF_550_LKT(68,16)=31.272000 !rg=2.22363 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,17,1:6)=(/ 25.678000,25.843000,26.103000,26.487000,27.157000,27.848000 /)
XPIZA_LKT(68,17,1:6)=(/ 0.707603,0.856673,0.993201,0.993125,0.975218,0.532275 /)
XCGA_LKT(68,17,1:6)=(/ 0.892937,0.859757,0.829893,0.822523,0.809660,0.925577 /)
XEXT_COEFF_550_LKT(68,17)=25.904000 !rg=2.22363 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,18,1:6)=(/ 21.330000,21.490000,21.600000,21.841000,22.431000,22.960000 /)
XPIZA_LKT(68,18,1:6)=(/ 0.693666,0.845415,0.992539,0.992796,0.973775,0.534997 /)
XCGA_LKT(68,18,1:6)=(/ 0.896057,0.862290,0.831867,0.824707,0.816197,0.928020 /)
XEXT_COEFF_550_LKT(68,18)=21.491000 !rg=2.22363 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,19,1:6)=(/ 17.712000,17.790000,17.931000,18.127000,18.516000,18.914000 /)
XPIZA_LKT(68,19,1:6)=(/ 0.679803,0.832721,0.991945,0.992324,0.972401,0.537606 /)
XCGA_LKT(68,19,1:6)=(/ 0.899757,0.864580,0.834253,0.827677,0.821267,0.930020 /)
XEXT_COEFF_550_LKT(68,19)=17.879000 !rg=2.22363 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,20,1:6)=(/ 14.699000,14.779000,14.862000,15.019000,15.312000,15.603000 /)
XPIZA_LKT(68,20,1:6)=(/ 0.665858,0.818952,0.991096,0.991664,0.970528,0.540023 /)
XCGA_LKT(68,20,1:6)=(/ 0.903140,0.867430,0.835747,0.830520,0.823280,0.931610 /)
XEXT_COEFF_550_LKT(68,20)=14.797000 !rg=2.22363 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,1,1:6)=(/ 268.600000,283.510000,269.010000,322.190000,265.720000,375.500000 /)
XPIZA_LKT(69,1,1:6)=(/ 0.815882,0.935225,0.997833,0.998392,0.994983,0.573726 /)
XCGA_LKT(69,1,1:6)=(/ 0.857760,0.818070,0.767987,0.782043,0.568407,0.800123 /)
XEXT_COEFF_550_LKT(69,1)=282.610000 !rg=2.40898 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,2,1:6)=(/ 253.890000,257.690000,267.310000,292.260000,263.730000,350.310000 /)
XPIZA_LKT(69,2,1:6)=(/ 0.807515,0.933504,0.997685,0.998105,0.994974,0.561086 /)
XCGA_LKT(69,2,1:6)=(/ 0.863480,0.817703,0.780423,0.761657,0.599553,0.801360 /)
XEXT_COEFF_550_LKT(69,2)=261.260000 !rg=2.40898 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,3,1:6)=(/ 236.170000,239.960000,245.060000,261.600000,261.950000,314.210000 /)
XPIZA_LKT(69,3,1:6)=(/ 0.799770,0.928355,0.997501,0.997997,0.994946,0.547818 /)
XCGA_LKT(69,3,1:6)=(/ 0.865213,0.820297,0.784240,0.746387,0.646690,0.806970 /)
XEXT_COEFF_550_LKT(69,3)=242.960000 !rg=2.40898 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,4,1:6)=(/ 212.850000,216.420000,221.050000,233.590000,244.510000,276.330000 /)
XPIZA_LKT(69,4,1:6)=(/ 0.794447,0.922717,0.996797,0.997721,0.994670,0.538035 /)
XCGA_LKT(69,4,1:6)=(/ 0.867803,0.824160,0.785780,0.746690,0.679427,0.817383 /)
XEXT_COEFF_550_LKT(69,4)=218.830000 !rg=2.40898 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,5,1:6)=(/ 188.100000,190.950000,195.410000,204.890000,219.170000,238.930000 /)
XPIZA_LKT(69,5,1:6)=(/ 0.791757,0.915478,0.996971,0.997642,0.993878,0.530752 /)
XCGA_LKT(69,5,1:6)=(/ 0.870193,0.826250,0.787583,0.762467,0.695117,0.829937 /)
XEXT_COEFF_550_LKT(69,5)=193.370000 !rg=2.40898 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,6,1:6)=(/ 163.190000,166.120000,170.500000,177.400000,189.180000,203.970000 /)
XPIZA_LKT(69,6,1:6)=(/ 0.790301,0.910761,0.996560,0.997380,0.993284,0.525212 /)
XCGA_LKT(69,6,1:6)=(/ 0.871417,0.831180,0.792733,0.759300,0.715640,0.842977 /)
XEXT_COEFF_550_LKT(69,6)=168.980000 !rg=2.40898 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,7,1:6)=(/ 140.530000,142.900000,145.200000,151.910000,160.650000,171.950000 /)
XPIZA_LKT(69,7,1:6)=(/ 0.788483,0.906440,0.996150,0.996822,0.992555,0.520971 /)
XCGA_LKT(69,7,1:6)=(/ 0.872437,0.835407,0.797633,0.765767,0.732820,0.856103 /)
XEXT_COEFF_550_LKT(69,7)=144.610000 !rg=2.40898 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,8,1:6)=(/ 119.920000,121.200000,124.450000,128.040000,135.970000,144.060000 /)
XPIZA_LKT(69,8,1:6)=(/ 0.785722,0.903516,0.995565,0.996581,0.991231,0.518326 /)
XCGA_LKT(69,8,1:6)=(/ 0.873773,0.838137,0.800113,0.775720,0.740403,0.868380 /)
XEXT_COEFF_550_LKT(69,8)=122.110000 !rg=2.40898 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,9,1:6)=(/ 101.300000,102.660000,103.990000,108.080000,114.400000,119.850000 /)
XPIZA_LKT(69,9,1:6)=(/ 0.782270,0.899711,0.995279,0.995983,0.989547,0.517033 /)
XCGA_LKT(69,9,1:6)=(/ 0.875790,0.840843,0.805240,0.784880,0.747050,0.879727 /)
XEXT_COEFF_550_LKT(69,9)=103.540000 !rg=2.40898 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,10,1:6)=(/ 85.219000,86.620000,88.058000,89.589000,95.714000,99.480000 /)
XPIZA_LKT(69,10,1:6)=(/ 0.777226,0.896638,0.995247,0.995619,0.987445,0.517041 /)
XCGA_LKT(69,10,1:6)=(/ 0.876957,0.843070,0.810703,0.794120,0.761320,0.889620 /)
XEXT_COEFF_550_LKT(69,10)=87.007000 !rg=2.40898 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,11,1:6)=(/ 71.430000,72.574000,73.166000,75.495000,79.072000,82.150000 /)
XPIZA_LKT(69,11,1:6)=(/ 0.770306,0.893609,0.994910,0.994831,0.985225,0.518017 /)
XCGA_LKT(69,11,1:6)=(/ 0.879217,0.846720,0.814190,0.793813,0.768303,0.898360 /)
XEXT_COEFF_550_LKT(69,11)=72.998000 !rg=2.40898 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,12,1:6)=(/ 59.576000,60.249000,61.258000,62.189000,65.549000,67.681000 /)
XPIZA_LKT(69,12,1:6)=(/ 0.762085,0.889441,0.994725,0.994563,0.982933,0.519796 /)
XCGA_LKT(69,12,1:6)=(/ 0.880943,0.848530,0.817770,0.804300,0.779447,0.905783 /)
XEXT_COEFF_550_LKT(69,12)=60.631000 !rg=2.40898 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,13,1:6)=(/ 49.665000,50.104000,50.743000,51.892000,53.978000,55.676000 /)
XPIZA_LKT(69,13,1:6)=(/ 0.752067,0.884767,0.994328,0.994181,0.981037,0.522139 /)
XCGA_LKT(69,13,1:6)=(/ 0.883377,0.851470,0.820503,0.807880,0.785270,0.911983 /)
XEXT_COEFF_550_LKT(69,13)=50.152000 !rg=2.40898 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,14,1:6)=(/ 41.281000,41.678000,42.231000,43.162000,44.415000,45.844000 /)
XPIZA_LKT(69,14,1:6)=(/ 0.740776,0.878633,0.994091,0.993692,0.979049,0.524829 /)
XCGA_LKT(69,14,1:6)=(/ 0.885770,0.853417,0.823543,0.809883,0.795150,0.917003 /)
XEXT_COEFF_550_LKT(69,14)=41.699000 !rg=2.40898 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,15,1:6)=(/ 34.373000,34.616000,35.045000,35.668000,36.636000,37.747000 /)
XPIZA_LKT(69,15,1:6)=(/ 0.728735,0.871071,0.993921,0.993447,0.977616,0.527689 /)
XCGA_LKT(69,15,1:6)=(/ 0.888610,0.855667,0.826803,0.815000,0.804293,0.921097 /)
XEXT_COEFF_550_LKT(69,15)=34.725000 !rg=2.40898 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,16,1:6)=(/ 28.486000,28.763000,28.895000,29.501000,30.197000,31.030000 /)
XPIZA_LKT(69,16,1:6)=(/ 0.715124,0.862399,0.993564,0.993277,0.975871,0.530619 /)
XCGA_LKT(69,16,1:6)=(/ 0.891467,0.858840,0.829293,0.818837,0.807580,0.924470 /)
XEXT_COEFF_550_LKT(69,16)=28.792000 !rg=2.40898 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,17,1:6)=(/ 23.670000,23.826000,24.025000,24.370000,24.987000,25.572000 /)
XPIZA_LKT(69,17,1:6)=(/ 0.701427,0.852016,0.992925,0.992943,0.974499,0.533458 /)
XCGA_LKT(69,17,1:6)=(/ 0.894500,0.860937,0.830603,0.822940,0.813290,0.927127 /)
XEXT_COEFF_550_LKT(69,17)=23.861000 !rg=2.40898 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,18,1:6)=(/ 19.682000,19.780000,19.979000,20.219000,20.626000,21.095000 /)
XPIZA_LKT(69,18,1:6)=(/ 0.687535,0.840012,0.992347,0.992511,0.973210,0.536161 /)
XCGA_LKT(69,18,1:6)=(/ 0.897803,0.863103,0.833227,0.825663,0.819773,0.929253 /)
XEXT_COEFF_550_LKT(69,18)=19.835000 !rg=2.40898 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,19,1:6)=(/ 16.326000,16.432000,16.507000,16.750000,16.993000,17.386000 /)
XPIZA_LKT(69,19,1:6)=(/ 0.673438,0.826969,0.991618,0.991973,0.971239,0.538719 /)
XCGA_LKT(69,19,1:6)=(/ 0.901263,0.866160,0.835043,0.828297,0.822207,0.930997 /)
XEXT_COEFF_550_LKT(69,19)=16.450000 !rg=2.40898 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,20,1:6)=(/ 13.555000,13.614000,13.699000,13.835000,14.102000,14.348000 /)
XPIZA_LKT(69,20,1:6)=(/ 0.659733,0.812198,0.990620,0.991366,0.969886,0.541063 /)
XCGA_LKT(69,20,1:6)=(/ 0.904863,0.868637,0.836183,0.830983,0.826387,0.932383 /)
XEXT_COEFF_550_LKT(69,20)=13.635000 !rg=2.40898 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,1,1:6)=(/ 244.370000,245.990000,242.770000,252.350000,213.160000,326.200000 /)
XPIZA_LKT(70,1,1:6)=(/ 0.803479,0.928890,0.997453,0.998099,0.993861,0.538466 /)
XCGA_LKT(70,1,1:6)=(/ 0.864483,0.817553,0.755683,0.768143,0.538060,0.802847 /)
XEXT_COEFF_550_LKT(70,1)=240.590000 !rg=2.60978 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,2,1:6)=(/ 233.730000,238.870000,244.280000,259.880000,230.120000,306.540000 /)
XPIZA_LKT(70,2,1:6)=(/ 0.796414,0.928389,0.997585,0.997809,0.994222,0.532207 /)
XCGA_LKT(70,2,1:6)=(/ 0.868687,0.820487,0.782803,0.733177,0.596923,0.806063 /)
XEXT_COEFF_550_LKT(70,2)=245.630000 !rg=2.60978 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,3,1:6)=(/ 216.860000,221.900000,227.650000,242.490000,235.640000,278.810000 /)
XPIZA_LKT(70,3,1:6)=(/ 0.788737,0.923279,0.997299,0.997555,0.994325,0.526400 /)
XCGA_LKT(70,3,1:6)=(/ 0.869853,0.822297,0.784420,0.753950,0.649440,0.815813 /)
XEXT_COEFF_550_LKT(70,3)=223.690000 !rg=2.60978 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,4,1:6)=(/ 195.730000,200.020000,205.330000,215.950000,221.560000,247.800000 /)
XPIZA_LKT(70,4,1:6)=(/ 0.786927,0.917120,0.996909,0.997496,0.994057,0.522925 /)
XCGA_LKT(70,4,1:6)=(/ 0.871323,0.825847,0.786567,0.757563,0.684367,0.828377 /)
XEXT_COEFF_550_LKT(70,4)=201.240000 !rg=2.60978 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,5,1:6)=(/ 172.980000,176.060000,178.440000,187.840000,197.260000,215.670000 /)
XPIZA_LKT(70,5,1:6)=(/ 0.788180,0.911502,0.996376,0.997382,0.993614,0.520092 /)
XCGA_LKT(70,5,1:6)=(/ 0.871973,0.831263,0.794860,0.759150,0.710487,0.841247 /)
XEXT_COEFF_550_LKT(70,5)=177.630000 !rg=2.60978 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,6,1:6)=(/ 150.470000,152.690000,157.360000,161.420000,173.360000,184.890000 /)
XPIZA_LKT(70,6,1:6)=(/ 0.787215,0.906830,0.996108,0.997146,0.992570,0.517714 /)
XCGA_LKT(70,6,1:6)=(/ 0.872987,0.834397,0.793687,0.773450,0.722090,0.853677 /)
XEXT_COEFF_550_LKT(70,6)=154.350000 !rg=2.60978 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,7,1:6)=(/ 129.380000,131.610000,134.090000,139.310000,147.250000,156.320000 /)
XPIZA_LKT(70,7,1:6)=(/ 0.785387,0.903454,0.995871,0.996587,0.991710,0.515860 /)
XCGA_LKT(70,7,1:6)=(/ 0.874317,0.836877,0.800133,0.774363,0.736697,0.865767 /)
XEXT_COEFF_550_LKT(70,7)=132.580000 !rg=2.60978 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,8,1:6)=(/ 110.210000,111.570000,114.630000,116.620000,125.300000,131.230000 /)
XPIZA_LKT(70,8,1:6)=(/ 0.783652,0.900906,0.995003,0.996303,0.989863,0.515014 /)
XCGA_LKT(70,8,1:6)=(/ 0.875210,0.840437,0.802920,0.786050,0.746600,0.876847 /)
XEXT_COEFF_550_LKT(70,8)=113.070000 !rg=2.60978 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,9,1:6)=(/ 93.292000,94.649000,95.468000,99.079000,103.700000,109.350000 /)
XPIZA_LKT(70,9,1:6)=(/ 0.779302,0.897845,0.995276,0.995842,0.988773,0.515103 /)
XCGA_LKT(70,9,1:6)=(/ 0.876737,0.843460,0.810137,0.786860,0.760333,0.886963 /)
XEXT_COEFF_550_LKT(70,9)=95.219000 !rg=2.60978 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,10,1:6)=(/ 78.647000,79.395000,80.496000,82.943000,87.182000,90.884000 /)
XPIZA_LKT(70,10,1:6)=(/ 0.773464,0.895494,0.995085,0.995384,0.986884,0.516136 /)
XCGA_LKT(70,10,1:6)=(/ 0.878747,0.845503,0.811900,0.796007,0.762917,0.895717 /)
XEXT_COEFF_550_LKT(70,10)=79.936000 !rg=2.60978 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,11,1:6)=(/ 65.879000,66.449000,67.661000,69.077000,72.426000,75.133000 /)
XPIZA_LKT(70,11,1:6)=(/ 0.766347,0.891550,0.994521,0.994812,0.984765,0.517889 /)
XCGA_LKT(70,11,1:6)=(/ 0.880627,0.848467,0.814757,0.799610,0.776953,0.903397 /)
XEXT_COEFF_550_LKT(70,11)=66.848000 !rg=2.60978 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,12,1:6)=(/ 55.005000,55.534000,56.133000,57.490000,59.919000,61.958000 /)
XPIZA_LKT(70,12,1:6)=(/ 0.757536,0.887795,0.994550,0.994402,0.982343,0.520216 /)
XCGA_LKT(70,12,1:6)=(/ 0.882560,0.850620,0.819823,0.806433,0.781647,0.909900 /)
XEXT_COEFF_550_LKT(70,12)=55.677000 !rg=2.60978 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,13,1:6)=(/ 45.688000,46.127000,46.780000,47.799000,49.332000,51.011000 /)
XPIZA_LKT(70,13,1:6)=(/ 0.746843,0.881827,0.994328,0.993922,0.980404,0.522938 /)
XCGA_LKT(70,13,1:6)=(/ 0.884567,0.852323,0.822227,0.808160,0.792013,0.915310 /)
XEXT_COEFF_550_LKT(70,13)=46.219000 !rg=2.60978 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,14,1:6)=(/ 38.132000,38.318000,38.839000,39.621000,40.854000,42.034000 /)
XPIZA_LKT(70,14,1:6)=(/ 0.735842,0.875467,0.994038,0.993933,0.978202,0.525859 /)
XCGA_LKT(70,14,1:6)=(/ 0.887367,0.855147,0.825637,0.815207,0.799443,0.919683 /)
XEXT_COEFF_550_LKT(70,14)=38.635000 !rg=2.60978 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,15,1:6)=(/ 31.662000,31.908000,32.179000,32.605000,33.730000,34.632000 /)
XPIZA_LKT(70,15,1:6)=(/ 0.722622,0.867601,0.993693,0.993592,0.976662,0.528847 /)
XCGA_LKT(70,15,1:6)=(/ 0.889870,0.857477,0.827847,0.817853,0.805087,0.923243 /)
XEXT_COEFF_550_LKT(70,15)=31.954000 !rg=2.60978 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,16,1:6)=(/ 26.248000,26.429000,26.699000,27.063000,27.804000,28.487000 /)
XPIZA_LKT(70,16,1:6)=(/ 0.708597,0.857726,0.993249,0.993132,0.975470,0.531833 /)
XCGA_LKT(70,16,1:6)=(/ 0.892920,0.859793,0.829913,0.820630,0.811470,0.926177 /)
XEXT_COEFF_550_LKT(70,16)=26.514000 !rg=2.60978 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,17,1:6)=(/ 21.845000,21.925000,22.181000,22.554000,22.985000,23.489000 /)
XPIZA_LKT(70,17,1:6)=(/ 0.695251,0.846655,0.992710,0.992655,0.973565,0.534667 /)
XCGA_LKT(70,17,1:6)=(/ 0.896233,0.862030,0.832490,0.824507,0.814373,0.928483 /)
XEXT_COEFF_550_LKT(70,17)=22.079000 !rg=2.60978 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,18,1:6)=(/ 18.144000,18.251000,18.370000,18.536000,19.023000,19.385000 /)
XPIZA_LKT(70,18,1:6)=(/ 0.681126,0.834636,0.992059,0.992369,0.971992,0.537327 /)
XCGA_LKT(70,18,1:6)=(/ 0.899310,0.864667,0.834017,0.827907,0.819753,0.930337 /)
XEXT_COEFF_550_LKT(70,18)=18.273000 !rg=2.60978 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,19,1:6)=(/ 15.050000,15.128000,15.255000,15.380000,15.703000,15.984000 /)
XPIZA_LKT(70,19,1:6)=(/ 0.666948,0.820431,0.991211,0.991651,0.970919,0.539817 /)
XCGA_LKT(70,19,1:6)=(/ 0.902893,0.867230,0.835840,0.829023,0.824607,0.931850 /)
XEXT_COEFF_550_LKT(70,19)=15.161000 !rg=2.60978 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,20,1:6)=(/ 12.508000,12.550000,12.652000,12.811000,12.977000,13.196000 /)
XPIZA_LKT(70,20,1:6)=(/ 0.653816,0.805296,0.990218,0.991073,0.968955,0.542077 /)
XCGA_LKT(70,20,1:6)=(/ 0.906687,0.869823,0.837537,0.832430,0.827053,0.933060 /)
XEXT_COEFF_550_LKT(70,20)=12.608000 !rg=2.60978 sigma=2.95 wvl=0.55
 
 END SUBROUTINE SALT_OPT_LKT_SET7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET8()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE

 
XEXT_COEFF_WVL_LKT(71,1,1:6)=(/ 224.530000,231.910000,245.620000,203.950000,196.790000,282.010000 /)
XPIZA_LKT(71,1,1:6)=(/ 0.790713,0.921400,0.997382,0.997912,0.993126,0.505764 /)
XCGA_LKT(71,1,1:6)=(/ 0.870790,0.819473,0.784873,0.762850,0.546440,0.804927 /)
XEXT_COEFF_550_LKT(71,1)=237.470000 !rg=2.82732 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,2,1:6)=(/ 216.230000,219.990000,227.320000,227.330000,216.400000,269.130000 /)
XPIZA_LKT(71,2,1:6)=(/ 0.782978,0.923579,0.997056,0.997794,0.993707,0.505644 /)
XCGA_LKT(71,2,1:6)=(/ 0.872930,0.825657,0.783663,0.773793,0.614090,0.813797 /)
XEXT_COEFF_550_LKT(71,2)=221.410000 !rg=2.82732 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,3,1:6)=(/ 199.660000,203.110000,210.230000,217.140000,217.830000,248.440000 /)
XPIZA_LKT(71,3,1:6)=(/ 0.781110,0.917288,0.996673,0.997214,0.993854,0.508292 /)
XCGA_LKT(71,3,1:6)=(/ 0.874513,0.826130,0.785223,0.762947,0.668460,0.827160 /)
XEXT_COEFF_550_LKT(71,3)=204.640000 !rg=2.82732 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,4,1:6)=(/ 179.880000,183.290000,189.160000,196.030000,202.420000,222.870000 /)
XPIZA_LKT(71,4,1:6)=(/ 0.783015,0.911115,0.996334,0.997089,0.993420,0.510623 /)
XCGA_LKT(71,4,1:6)=(/ 0.874430,0.829523,0.788157,0.763127,0.695187,0.840563 /)
XEXT_COEFF_550_LKT(71,4)=184.490000 !rg=2.82732 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,5,1:6)=(/ 159.360000,162.360000,165.530000,173.060000,180.170000,195.040000 /)
XPIZA_LKT(71,5,1:6)=(/ 0.783823,0.907227,0.996435,0.997068,0.992938,0.511576 /)
XCGA_LKT(71,5,1:6)=(/ 0.874403,0.832783,0.795443,0.767500,0.717610,0.852827 /)
XEXT_COEFF_550_LKT(71,5)=162.520000 !rg=2.82732 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,6,1:6)=(/ 138.320000,140.960000,143.800000,147.890000,157.800000,167.800000 /)
XPIZA_LKT(71,6,1:6)=(/ 0.785453,0.903415,0.996018,0.996895,0.991322,0.511850 /)
XCGA_LKT(71,6,1:6)=(/ 0.874063,0.836663,0.798860,0.780843,0.723490,0.864130 /)
XEXT_COEFF_550_LKT(71,6)=140.980000 !rg=2.82732 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,7,1:6)=(/ 119.000000,121.210000,122.990000,127.730000,134.460000,142.230000 /)
XPIZA_LKT(71,7,1:6)=(/ 0.783754,0.900068,0.995763,0.995853,0.990392,0.512007 /)
XCGA_LKT(71,7,1:6)=(/ 0.875660,0.839550,0.803513,0.776780,0.737997,0.874923 /)
XEXT_COEFF_550_LKT(71,7)=122.320000 !rg=2.82732 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,8,1:6)=(/ 101.560000,102.850000,104.840000,107.840000,114.100000,119.620000 /)
XPIZA_LKT(71,8,1:6)=(/ 0.780393,0.898715,0.995044,0.996018,0.989461,0.512669 /)
XCGA_LKT(71,8,1:6)=(/ 0.876643,0.842313,0.806197,0.789083,0.747037,0.884703 /)
XEXT_COEFF_550_LKT(71,8)=103.060000 !rg=2.82732 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,9,1:6)=(/ 86.058000,87.202000,88.285000,91.362000,95.148000,99.831000 /)
XPIZA_LKT(71,9,1:6)=(/ 0.776512,0.896120,0.995134,0.995324,0.987582,0.513917 /)
XCGA_LKT(71,9,1:6)=(/ 0.878030,0.843733,0.812050,0.790237,0.765423,0.893583 /)
XEXT_COEFF_550_LKT(71,9)=87.250000 !rg=2.82732 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,10,1:6)=(/ 72.407000,73.348000,74.094000,76.342000,79.556000,83.071000 /)
XPIZA_LKT(71,10,1:6)=(/ 0.770510,0.892889,0.994737,0.994773,0.985461,0.515794 /)
XCGA_LKT(71,10,1:6)=(/ 0.879067,0.847053,0.814733,0.796183,0.772633,0.901220 /)
XEXT_COEFF_550_LKT(71,10)=73.417000 !rg=2.82732 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,11,1:6)=(/ 60.666000,61.297000,62.100000,63.434000,66.052000,68.747000 /)
XPIZA_LKT(71,11,1:6)=(/ 0.762694,0.889673,0.994360,0.994755,0.983700,0.518171 /)
XCGA_LKT(71,11,1:6)=(/ 0.881573,0.849703,0.815137,0.804600,0.778510,0.907910 /)
XEXT_COEFF_550_LKT(71,11)=61.823000 !rg=2.82732 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,12,1:6)=(/ 50.583000,51.158000,51.801000,52.963000,54.816000,56.743000 /)
XPIZA_LKT(71,12,1:6)=(/ 0.752996,0.885007,0.994264,0.994251,0.981147,0.520935 /)
XCGA_LKT(71,12,1:6)=(/ 0.883273,0.851373,0.820783,0.807067,0.788133,0.913557 /)
XEXT_COEFF_550_LKT(71,12)=51.175000 !rg=2.82732 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,13,1:6)=(/ 42.187000,42.465000,43.075000,43.936000,45.335000,46.754000 /)
XPIZA_LKT(71,13,1:6)=(/ 0.742311,0.879270,0.994147,0.994041,0.978985,0.523940 /)
XCGA_LKT(71,13,1:6)=(/ 0.886000,0.854313,0.824880,0.812963,0.795360,0.918253 /)
XEXT_COEFF_550_LKT(71,13)=42.735000 !rg=2.82732 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,14,1:6)=(/ 35.088000,35.327000,35.771000,36.236000,37.493000,38.554000 /)
XPIZA_LKT(71,14,1:6)=(/ 0.729625,0.872107,0.993723,0.993800,0.977801,0.527022 /)
XCGA_LKT(71,14,1:6)=(/ 0.888260,0.856563,0.826100,0.816410,0.801750,0.922043 /)
XEXT_COEFF_550_LKT(71,14)=35.363000 !rg=2.82732 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,15,1:6)=(/ 29.181000,29.415000,29.722000,30.000000,31.096000,31.784000 /)
XPIZA_LKT(71,15,1:6)=(/ 0.716644,0.863361,0.993514,0.993262,0.975510,0.530087 /)
XCGA_LKT(71,15,1:6)=(/ 0.891183,0.858120,0.829187,0.820917,0.807713,0.925127 /)
XEXT_COEFF_550_LKT(71,15)=29.532000 !rg=2.82732 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,16,1:6)=(/ 24.228000,24.396000,24.571000,24.907000,25.431000,26.159000 /)
XPIZA_LKT(71,16,1:6)=(/ 0.702828,0.853228,0.992931,0.993116,0.974613,0.533085 /)
XCGA_LKT(71,16,1:6)=(/ 0.894317,0.860963,0.830273,0.823850,0.813247,0.927673 /)
XEXT_COEFF_550_LKT(71,16)=24.502000 !rg=2.82732 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,17,1:6)=(/ 20.126000,20.271000,20.401000,20.654000,21.034000,21.581000 /)
XPIZA_LKT(71,17,1:6)=(/ 0.688618,0.841755,0.992459,0.992529,0.973294,0.535886 /)
XCGA_LKT(71,17,1:6)=(/ 0.897557,0.863533,0.832040,0.826220,0.818503,0.929673 /)
XEXT_COEFF_550_LKT(71,17)=20.305000 !rg=2.82732 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,18,1:6)=(/ 16.724000,16.820000,16.976000,17.080000,17.521000,17.818000 /)
XPIZA_LKT(71,18,1:6)=(/ 0.674716,0.828477,0.991707,0.992064,0.971399,0.538484 /)
XCGA_LKT(71,18,1:6)=(/ 0.900847,0.865500,0.835220,0.828963,0.822163,0.931280 /)
XEXT_COEFF_550_LKT(71,18)=16.873000 !rg=2.82732 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,19,1:6)=(/ 13.893000,13.957000,14.039000,14.203000,14.368000,14.698000 /)
XPIZA_LKT(71,19,1:6)=(/ 0.661075,0.814179,0.990813,0.991495,0.970238,0.540892 /)
XCGA_LKT(71,19,1:6)=(/ 0.904617,0.868413,0.836210,0.831523,0.826413,0.932597 /)
XEXT_COEFF_550_LKT(71,19)=14.001000 !rg=2.82732 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,20,1:6)=(/ 11.533000,11.592000,11.640000,11.740000,11.908000,12.139000 /)
XPIZA_LKT(71,20,1:6)=(/ 0.647735,0.798596,0.989701,0.990625,0.968166,0.543057 /)
XCGA_LKT(71,20,1:6)=(/ 0.908320,0.871363,0.837340,0.833130,0.829373,0.933650 /)
XEXT_COEFF_550_LKT(71,20)=11.606000 !rg=2.82732 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,1,1:6)=(/ 209.070000,210.530000,221.430000,210.160000,199.720000,247.020000 /)
XPIZA_LKT(72,1,1:6)=(/ 0.773554,0.920293,0.997169,0.997635,0.993240,0.476788 /)
XCGA_LKT(72,1,1:6)=(/ 0.877087,0.827000,0.794170,0.767000,0.614703,0.814727 /)
XEXT_COEFF_550_LKT(72,1)=215.770000 !rg=3.06299 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,2,1:6)=(/ 197.530000,202.850000,206.220000,210.290000,211.500000,238.730000 /)
XPIZA_LKT(72,2,1:6)=(/ 0.771865,0.918028,0.997228,0.997837,0.993619,0.485426 /)
XCGA_LKT(72,2,1:6)=(/ 0.877493,0.827313,0.796670,0.776550,0.663023,0.826277 /)
XEXT_COEFF_550_LKT(72,2)=203.010000 !rg=3.06299 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,3,1:6)=(/ 183.650000,186.360000,190.280000,196.850000,205.940000,222.940000 /)
XPIZA_LKT(72,3,1:6)=(/ 0.776222,0.911275,0.996949,0.997042,0.993212,0.495181 /)
XCGA_LKT(72,3,1:6)=(/ 0.877410,0.831380,0.796520,0.771833,0.690957,0.840977 /)
XEXT_COEFF_550_LKT(72,3)=187.680000 !rg=3.06299 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,4,1:6)=(/ 165.620000,168.290000,172.490000,177.160000,188.220000,201.240000 /)
XPIZA_LKT(72,4,1:6)=(/ 0.780405,0.906017,0.995755,0.997303,0.992394,0.501680 /)
XCGA_LKT(72,4,1:6)=(/ 0.876437,0.834817,0.795487,0.775403,0.709393,0.853477 /)
XEXT_COEFF_550_LKT(72,4)=168.960000 !rg=3.06299 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,5,1:6)=(/ 146.320000,148.770000,152.590000,158.010000,164.670000,176.790000 /)
XPIZA_LKT(72,5,1:6)=(/ 0.782123,0.902145,0.995873,0.996060,0.991964,0.505358 /)
XCGA_LKT(72,5,1:6)=(/ 0.876153,0.836327,0.797337,0.770747,0.723477,0.864270 /)
XEXT_COEFF_550_LKT(72,5)=149.880000 !rg=3.06299 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,6,1:6)=(/ 127.430000,128.700000,131.860000,136.120000,145.020000,152.510000 /)
XPIZA_LKT(72,6,1:6)=(/ 0.782999,0.900944,0.995659,0.996680,0.990605,0.507639 /)
XCGA_LKT(72,6,1:6)=(/ 0.876273,0.839087,0.800910,0.781127,0.734983,0.874073 /)
XEXT_COEFF_550_LKT(72,6)=129.720000 !rg=3.06299 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,7,1:6)=(/ 109.720000,111.030000,113.630000,116.300000,123.500000,129.540000 /)
XPIZA_LKT(72,7,1:6)=(/ 0.780879,0.898261,0.995334,0.996250,0.989821,0.509374 /)
XCGA_LKT(72,7,1:6)=(/ 0.877173,0.842627,0.804297,0.785660,0.749110,0.883420 /)
XEXT_COEFF_550_LKT(72,7)=111.870000 !rg=3.06299 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,8,1:6)=(/ 93.417000,94.633000,96.181000,99.941000,104.290000,109.130000 /)
XPIZA_LKT(72,8,1:6)=(/ 0.778424,0.896501,0.995149,0.995137,0.988106,0.511235 /)
XCGA_LKT(72,8,1:6)=(/ 0.877503,0.843317,0.808860,0.785267,0.758300,0.891873 /)
XEXT_COEFF_550_LKT(72,8)=94.723000 !rg=3.06299 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,9,1:6)=(/ 79.171000,80.013000,81.194000,83.405000,87.276000,91.196000 /)
XPIZA_LKT(72,9,1:6)=(/ 0.773549,0.894113,0.994583,0.995013,0.986106,0.513402 /)
XCGA_LKT(72,9,1:6)=(/ 0.879197,0.846057,0.812723,0.792940,0.766887,0.899547 /)
XEXT_COEFF_550_LKT(72,9)=80.692000 !rg=3.06299 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,10,1:6)=(/ 66.668000,67.356000,68.570000,70.216000,72.754000,75.970000 /)
XPIZA_LKT(72,10,1:6)=(/ 0.766799,0.891317,0.994586,0.994876,0.984873,0.515950 /)
XCGA_LKT(72,10,1:6)=(/ 0.880550,0.848233,0.815703,0.801743,0.780093,0.906137 /)
XEXT_COEFF_550_LKT(72,10)=67.651000 !rg=3.06299 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,11,1:6)=(/ 55.844000,56.461000,57.055000,58.397000,60.199000,62.930000 /)
XPIZA_LKT(72,11,1:6)=(/ 0.758128,0.887706,0.994443,0.994395,0.982737,0.518812 /)
XCGA_LKT(72,11,1:6)=(/ 0.882450,0.851383,0.819847,0.805260,0.787870,0.911913 /)
XEXT_COEFF_550_LKT(72,11)=56.837000 !rg=3.06299 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,12,1:6)=(/ 46.704000,46.990000,47.692000,48.670000,50.311000,51.986000 /)
XPIZA_LKT(72,12,1:6)=(/ 0.748642,0.882423,0.994276,0.994133,0.980038,0.521902 /)
XCGA_LKT(72,12,1:6)=(/ 0.884760,0.852907,0.823073,0.811470,0.793517,0.916787 /)
XEXT_COEFF_550_LKT(72,12)=47.344000 !rg=3.06299 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,13,1:6)=(/ 38.837000,39.212000,39.734000,40.287000,41.576000,42.868000 /)
XPIZA_LKT(72,13,1:6)=(/ 0.736317,0.876203,0.993815,0.993889,0.978373,0.525105 /)
XCGA_LKT(72,13,1:6)=(/ 0.887033,0.855953,0.824147,0.815300,0.798773,0.920840 /)
XEXT_COEFF_550_LKT(72,13)=39.213000 !rg=3.06299 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,14,1:6)=(/ 32.380000,32.536000,33.069000,33.354000,34.383000,35.372000 /)
XPIZA_LKT(72,14,1:6)=(/ 0.724283,0.868266,0.993719,0.993468,0.976896,0.528288 /)
XCGA_LKT(72,14,1:6)=(/ 0.889823,0.857347,0.828000,0.818717,0.805610,0.924113 /)
XEXT_COEFF_550_LKT(72,14)=32.726000 !rg=3.06299 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,15,1:6)=(/ 26.934000,27.060000,27.315000,27.756000,28.467000,29.179000 /)
XPIZA_LKT(72,15,1:6)=(/ 0.710651,0.859143,0.993202,0.993169,0.975309,0.531380 /)
XCGA_LKT(72,15,1:6)=(/ 0.892773,0.859443,0.829847,0.822517,0.811137,0.926777 /)
XEXT_COEFF_550_LKT(72,15)=27.178000 !rg=3.06299 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,16,1:6)=(/ 22.329000,22.459000,22.607000,22.963000,23.368000,24.028000 /)
XPIZA_LKT(72,16,1:6)=(/ 0.696309,0.848186,0.992855,0.992735,0.973754,0.534358 /)
XCGA_LKT(72,16,1:6)=(/ 0.895760,0.862327,0.832597,0.823623,0.817450,0.928980 /)
XEXT_COEFF_550_LKT(72,16)=22.528000 !rg=3.06299 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,17,1:6)=(/ 18.563000,18.634000,18.811000,19.053000,19.363000,19.832000 /)
XPIZA_LKT(72,17,1:6)=(/ 0.682316,0.835885,0.992132,0.992359,0.972547,0.537103 /)
XCGA_LKT(72,17,1:6)=(/ 0.899110,0.864413,0.833777,0.827470,0.821273,0.930710 /)
XEXT_COEFF_550_LKT(72,17)=18.750000 !rg=3.06299 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,18,1:6)=(/ 15.432000,15.500000,15.599000,15.785000,16.079000,16.382000 /)
XPIZA_LKT(72,18,1:6)=(/ 0.668493,0.822311,0.991333,0.991823,0.970856,0.539621 /)
XCGA_LKT(72,18,1:6)=(/ 0.902590,0.866867,0.835480,0.830577,0.824823,0.932103 /)
XEXT_COEFF_550_LKT(72,18)=15.539000 !rg=3.06299 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,19,1:6)=(/ 12.808000,12.868000,12.929000,13.072000,13.225000,13.518000 /)
XPIZA_LKT(72,19,1:6)=(/ 0.654848,0.807296,0.990381,0.991094,0.969111,0.541936 /)
XCGA_LKT(72,19,1:6)=(/ 0.906260,0.869747,0.837313,0.830980,0.828423,0.933247 /)
XEXT_COEFF_550_LKT(72,19)=12.889000 !rg=3.06299 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,20,1:6)=(/ 10.640000,10.671000,10.748000,10.834000,10.978000,11.168000 /)
XPIZA_LKT(72,20,1:6)=(/ 0.641969,0.791176,0.989100,0.990351,0.967296,0.544003 /)
XCGA_LKT(72,20,1:6)=(/ 0.910060,0.872413,0.838517,0.834223,0.831517,0.934163 /)
XEXT_COEFF_550_LKT(72,20)=10.716000 !rg=3.06299 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,1,1:6)=(/ 189.330000,194.670000,186.000000,223.810000,213.940000,220.280000 /)
XPIZA_LKT(73,1,1:6)=(/ 0.753830,0.916140,0.997026,0.996396,0.993823,0.460817 /)
XCGA_LKT(73,1,1:6)=(/ 0.883670,0.835627,0.791023,0.762710,0.699080,0.830850 /)
XEXT_COEFF_550_LKT(73,1)=193.660000 !rg=3.31831 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,2,1:6)=(/ 182.080000,184.680000,189.220000,195.290000,209.840000,214.880000 /)
XPIZA_LKT(73,2,1:6)=(/ 0.762718,0.912836,0.995536,0.997678,0.993534,0.474507 /)
XCGA_LKT(73,2,1:6)=(/ 0.881890,0.831703,0.791820,0.775637,0.710820,0.843500 /)
XEXT_COEFF_550_LKT(73,2)=187.010000 !rg=3.31831 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,3,1:6)=(/ 169.320000,171.140000,174.540000,182.710000,192.320000,201.660000 /)
XPIZA_LKT(73,3,1:6)=(/ 0.772732,0.906008,0.996274,0.997253,0.992834,0.487724 /)
XCGA_LKT(73,3,1:6)=(/ 0.879103,0.836347,0.795723,0.775520,0.711147,0.856287 /)
XEXT_COEFF_550_LKT(73,3)=173.830000 !rg=3.31831 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,4,1:6)=(/ 152.350000,154.600000,157.310000,164.280000,172.600000,182.480000 /)
XPIZA_LKT(73,4,1:6)=(/ 0.777978,0.901282,0.996347,0.997150,0.992181,0.496138 /)
XCGA_LKT(73,4,1:6)=(/ 0.877597,0.838167,0.797973,0.779680,0.718617,0.866373 /)
XEXT_COEFF_550_LKT(73,4)=156.510000 !rg=3.31831 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,5,1:6)=(/ 134.850000,137.220000,139.800000,143.430000,152.810000,160.630000 /)
XPIZA_LKT(73,5,1:6)=(/ 0.781547,0.898536,0.995918,0.996578,0.990909,0.501368 /)
XCGA_LKT(73,5,1:6)=(/ 0.877120,0.840073,0.802353,0.781293,0.734537,0.875143 /)
XEXT_COEFF_550_LKT(73,5)=137.470000 !rg=3.31831 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,6,1:6)=(/ 117.430000,118.530000,122.010000,125.840000,131.650000,138.820000 /)
XPIZA_LKT(73,6,1:6)=(/ 0.781498,0.897494,0.995254,0.996297,0.989529,0.504973 /)
XCGA_LKT(73,6,1:6)=(/ 0.877447,0.842600,0.803023,0.784817,0.743490,0.883270 /)
XEXT_COEFF_550_LKT(73,6)=120.380000 !rg=3.31831 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,7,1:6)=(/ 100.950000,102.300000,103.990000,107.100000,112.450000,118.100000 /)
XPIZA_LKT(73,7,1:6)=(/ 0.779758,0.896356,0.994966,0.995945,0.989379,0.507852 /)
XCGA_LKT(73,7,1:6)=(/ 0.877787,0.844227,0.805080,0.792123,0.751657,0.891143 /)
XEXT_COEFF_550_LKT(73,7)=103.620000 !rg=3.31831 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,8,1:6)=(/ 86.105000,87.045000,88.744000,91.232000,94.910000,99.635000 /)
XPIZA_LKT(73,8,1:6)=(/ 0.776651,0.895053,0.994898,0.995314,0.987559,0.510615 /)
XCGA_LKT(73,8,1:6)=(/ 0.878560,0.846127,0.812963,0.795437,0.766650,0.898313 /)
XEXT_COEFF_550_LKT(73,8)=87.603000 !rg=3.31831 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,9,1:6)=(/ 72.907000,73.895000,74.987000,76.285000,80.413000,83.359000 /)
XPIZA_LKT(73,9,1:6)=(/ 0.770546,0.892148,0.994781,0.994833,0.985071,0.513485 /)
XCGA_LKT(73,9,1:6)=(/ 0.880040,0.847890,0.815013,0.799767,0.773407,0.904853 /)
XEXT_COEFF_550_LKT(73,9)=74.087000 !rg=3.31831 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,10,1:6)=(/ 61.451000,62.039000,63.127000,64.010000,66.854000,69.513000 /)
XPIZA_LKT(73,10,1:6)=(/ 0.763173,0.890035,0.994368,0.994686,0.983686,0.516531 /)
XCGA_LKT(73,10,1:6)=(/ 0.881490,0.850123,0.815740,0.803887,0.781717,0.910487 /)
XEXT_COEFF_550_LKT(73,10)=62.125000 !rg=3.31831 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,11,1:6)=(/ 51.564000,52.031000,52.696000,53.768000,55.602000,57.633000 /)
XPIZA_LKT(73,11,1:6)=(/ 0.753964,0.885055,0.994350,0.994245,0.981595,0.519753 /)
XCGA_LKT(73,11,1:6)=(/ 0.883933,0.851870,0.821120,0.808393,0.791243,0.915433 /)
XEXT_COEFF_550_LKT(73,11)=52.147000 !rg=3.31831 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,12,1:6)=(/ 43.000000,43.332000,43.896000,44.536000,46.193000,47.648000 /)
XPIZA_LKT(73,12,1:6)=(/ 0.743050,0.879952,0.994194,0.993996,0.979717,0.523069 /)
XCGA_LKT(73,12,1:6)=(/ 0.885657,0.854430,0.823863,0.812560,0.795863,0.919617 /)
XEXT_COEFF_550_LKT(73,12)=43.312000 !rg=3.31831 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,13,1:6)=(/ 35.860000,36.007000,36.571000,37.118000,38.164000,39.318000 /)
XPIZA_LKT(73,13,1:6)=(/ 0.731335,0.872680,0.993874,0.993391,0.977762,0.526401 /)
XCGA_LKT(73,13,1:6)=(/ 0.888493,0.856477,0.826453,0.816660,0.803883,0.923100 /)
XEXT_COEFF_550_LKT(73,13)=36.286000 !rg=3.31831 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,14,1:6)=(/ 29.849000,30.033000,30.292000,30.875000,31.755000,32.463000 /)
XPIZA_LKT(73,14,1:6)=(/ 0.718012,0.864449,0.993520,0.993126,0.975465,0.529626 /)
XCGA_LKT(73,14,1:6)=(/ 0.891190,0.858553,0.829113,0.819453,0.806440,0.925917 /)
XEXT_COEFF_550_LKT(73,14)=30.041000 !rg=3.31831 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,15,1:6)=(/ 24.799000,25.002000,25.133000,25.432000,26.148000,26.794000 /)
XPIZA_LKT(73,15,1:6)=(/ 0.704114,0.854531,0.992965,0.992948,0.974273,0.532708 /)
XCGA_LKT(73,15,1:6)=(/ 0.893807,0.860963,0.831140,0.823480,0.813527,0.928217 /)
XEXT_COEFF_550_LKT(73,15)=25.021000 !rg=3.31831 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,16,1:6)=(/ 20.608000,20.703000,20.893000,21.176000,21.591000,22.075000 /)
XPIZA_LKT(73,16,1:6)=(/ 0.690150,0.842982,0.992451,0.992614,0.973239,0.535637 /)
XCGA_LKT(73,16,1:6)=(/ 0.897393,0.863060,0.833890,0.826413,0.819360,0.930117 /)
XEXT_COEFF_550_LKT(73,16)=20.778000 !rg=3.31831 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,17,1:6)=(/ 17.118000,17.213000,17.347000,17.546000,17.816000,18.229000 /)
XPIZA_LKT(73,17,1:6)=(/ 0.675953,0.830247,0.991811,0.992193,0.971789,0.538305 /)
XCGA_LKT(73,17,1:6)=(/ 0.900643,0.865640,0.834707,0.829857,0.822500,0.931613 /)
XEXT_COEFF_550_LKT(73,17)=17.234000 !rg=3.31831 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,18,1:6)=(/ 14.230000,14.307000,14.358000,14.491000,14.779000,15.064000 /)
XPIZA_LKT(73,18,1:6)=(/ 0.662237,0.815954,0.990888,0.991385,0.969929,0.540732 /)
XCGA_LKT(73,18,1:6)=(/ 0.904093,0.868293,0.836493,0.830940,0.825853,0.932820 /)
XEXT_COEFF_550_LKT(73,18)=14.318000 !rg=3.31831 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,19,1:6)=(/ 11.820000,11.865000,11.937000,12.046000,12.236000,12.435000 /)
XPIZA_LKT(73,19,1:6)=(/ 0.648972,0.800317,0.989841,0.990814,0.968466,0.542946 /)
XCGA_LKT(73,19,1:6)=(/ 0.908063,0.870863,0.838417,0.833397,0.830103,0.933813 /)
XEXT_COEFF_550_LKT(73,19)=11.895000 !rg=3.31831 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,20,1:6)=(/ 9.814100,9.851700,9.905300,9.989200,10.103000,10.277000 /)
XPIZA_LKT(73,20,1:6)=(/ 0.636294,0.783948,0.988582,0.989868,0.966165,0.544909 /)
XCGA_LKT(73,20,1:6)=(/ 0.911780,0.873783,0.839237,0.835703,0.832497,0.934610 /)
XEXT_COEFF_550_LKT(73,20)=9.866700 !rg=3.31831 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,1,1:6)=(/ 175.340000,180.900000,187.770000,204.320000,225.350000,201.310000 /)
XPIZA_LKT(74,1,1:6)=(/ 0.744113,0.908911,0.995994,0.997618,0.994141,0.459177 /)
XCGA_LKT(74,1,1:6)=(/ 0.888803,0.830623,0.795833,0.786537,0.754303,0.854260 /)
XEXT_COEFF_550_LKT(74,1)=181.350000 !rg=3.59491 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,2,1:6)=(/ 167.210000,169.590000,173.190000,181.590000,202.490000,196.210000 /)
XPIZA_LKT(74,2,1:6)=(/ 0.763051,0.905212,0.996524,0.997403,0.993336,0.473076 /)
XCGA_LKT(74,2,1:6)=(/ 0.882647,0.837030,0.797797,0.780753,0.742923,0.863053 /)
XEXT_COEFF_550_LKT(74,2)=172.060000 !rg=3.59491 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,3,1:6)=(/ 155.710000,158.650000,159.460000,167.700000,176.490000,183.680000 /)
XPIZA_LKT(74,3,1:6)=(/ 0.773701,0.899177,0.996527,0.997104,0.992793,0.485374 /)
XCGA_LKT(74,3,1:6)=(/ 0.879487,0.838680,0.802360,0.770487,0.741983,0.871500 /)
XEXT_COEFF_550_LKT(74,3)=159.160000 !rg=3.59491 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,4,1:6)=(/ 140.320000,142.900000,143.860000,150.940000,156.410000,166.060000 /)
XPIZA_LKT(74,4,1:6)=(/ 0.778766,0.897120,0.996077,0.996871,0.992153,0.493606 /)
XCGA_LKT(74,4,1:6)=(/ 0.878170,0.841090,0.804410,0.775807,0.745860,0.878467 /)
XEXT_COEFF_550_LKT(74,4)=143.770000 !rg=3.59491 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,5,1:6)=(/ 124.100000,125.740000,127.400000,133.010000,139.560000,146.270000 /)
XPIZA_LKT(74,5,1:6)=(/ 0.780697,0.896521,0.995732,0.996486,0.990825,0.499340 /)
XCGA_LKT(74,5,1:6)=(/ 0.877800,0.842380,0.804700,0.784320,0.737697,0.885077 /)
XEXT_COEFF_550_LKT(74,5)=126.890000 !rg=3.59491 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,6,1:6)=(/ 107.910000,109.360000,112.050000,114.690000,119.380000,126.550000 /)
XPIZA_LKT(74,6,1:6)=(/ 0.780214,0.895630,0.995381,0.996145,0.989937,0.503677 /)
XCGA_LKT(74,6,1:6)=(/ 0.878007,0.844833,0.806327,0.788007,0.755657,0.891543 /)
XEXT_COEFF_550_LKT(74,6)=110.640000 !rg=3.59491 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,7,1:6)=(/ 92.934000,94.234000,95.298000,98.546000,101.840000,107.790000 /)
XPIZA_LKT(74,7,1:6)=(/ 0.777579,0.894895,0.994897,0.995361,0.988923,0.507308 /)
XCGA_LKT(74,7,1:6)=(/ 0.878480,0.846363,0.810997,0.791530,0.767197,0.898017 /)
XEXT_COEFF_550_LKT(74,7)=95.149000 !rg=3.59491 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,8,1:6)=(/ 79.333000,80.087000,81.735000,83.279000,86.890000,91.037000 /)
XPIZA_LKT(74,8,1:6)=(/ 0.773093,0.893387,0.994602,0.995208,0.986946,0.510698 /)
XCGA_LKT(74,8,1:6)=(/ 0.879297,0.848053,0.813623,0.797760,0.770427,0.904003 /)
XEXT_COEFF_550_LKT(74,8)=80.287000 !rg=3.59491 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,9,1:6)=(/ 67.200000,67.767000,68.402000,70.584000,73.324000,76.245000 /)
XPIZA_LKT(74,9,1:6)=(/ 0.767325,0.891113,0.994738,0.994528,0.984409,0.514068 /)
XCGA_LKT(74,9,1:6)=(/ 0.881263,0.849003,0.817643,0.801087,0.775567,0.909520 /)
XEXT_COEFF_550_LKT(74,9)=68.270000 !rg=3.59491 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,10,1:6)=(/ 56.549000,57.203000,57.962000,58.767000,61.686000,63.636000 /)
XPIZA_LKT(74,10,1:6)=(/ 0.759080,0.887106,0.994345,0.994461,0.981876,0.517472 /)
XCGA_LKT(74,10,1:6)=(/ 0.882447,0.850710,0.820787,0.808813,0.785713,0.914293 /)
XEXT_COEFF_550_LKT(74,10)=57.436000 !rg=3.59491 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,11,1:6)=(/ 47.438000,47.957000,48.236000,49.434000,51.327000,52.804000 /)
XPIZA_LKT(74,11,1:6)=(/ 0.748815,0.882784,0.994419,0.993900,0.979292,0.520933 /)
XCGA_LKT(74,11,1:6)=(/ 0.884847,0.853390,0.824043,0.808597,0.790237,0.918510 /)
XEXT_COEFF_550_LKT(74,11)=48.292000 !rg=3.59491 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,12,1:6)=(/ 39.618000,39.891000,40.567000,40.978000,42.426000,43.688000 /)
XPIZA_LKT(74,12,1:6)=(/ 0.737688,0.876680,0.993793,0.993606,0.978446,0.524395 /)
XCGA_LKT(74,12,1:6)=(/ 0.886993,0.855320,0.825430,0.815647,0.799560,0.922083 /)
XEXT_COEFF_550_LKT(74,12)=40.098000 !rg=3.59491 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,13,1:6)=(/ 33.021000,33.258000,33.634000,34.143000,35.152000,36.074000 /)
XPIZA_LKT(74,13,1:6)=(/ 0.725188,0.869295,0.993693,0.993540,0.976656,0.527786 /)
XCGA_LKT(74,13,1:6)=(/ 0.889570,0.857550,0.828403,0.819530,0.804353,0.925070 /)
XEXT_COEFF_550_LKT(74,13)=33.299000 !rg=3.59491 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,14,1:6)=(/ 27.477000,27.673000,27.952000,28.356000,28.999000,29.802000 /)
XPIZA_LKT(74,14,1:6)=(/ 0.711574,0.860117,0.993429,0.993152,0.975486,0.531012 /)
XCGA_LKT(74,14,1:6)=(/ 0.892380,0.859347,0.830523,0.820660,0.811817,0.927487 /)
XEXT_COEFF_550_LKT(74,14)=27.707000 !rg=3.59491 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,15,1:6)=(/ 22.885000,22.997000,23.256000,23.520000,24.044000,24.611000 /)
XPIZA_LKT(74,15,1:6)=(/ 0.697970,0.849564,0.992851,0.992878,0.974003,0.534051 /)
XCGA_LKT(74,15,1:6)=(/ 0.895550,0.861643,0.832340,0.824910,0.818083,0.929467 /)
XEXT_COEFF_550_LKT(74,15)=23.099000 !rg=3.59491 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,16,1:6)=(/ 18.991000,19.140000,19.202000,19.500000,19.861000,20.287000 /)
XPIZA_LKT(74,16,1:6)=(/ 0.683582,0.837822,0.992214,0.992421,0.972273,0.536906 /)
XCGA_LKT(74,16,1:6)=(/ 0.898837,0.864620,0.834170,0.826823,0.820283,0.931103 /)
XEXT_COEFF_550_LKT(74,16)=19.153000 !rg=3.59491 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,17,1:6)=(/ 15.786000,15.855000,15.942000,16.180000,16.449000,16.759000 /)
XPIZA_LKT(74,17,1:6)=(/ 0.669589,0.824106,0.991478,0.991670,0.970947,0.539483 /)
XCGA_LKT(74,17,1:6)=(/ 0.902270,0.866853,0.835523,0.829207,0.824493,0.932400 /)
XEXT_COEFF_550_LKT(74,17)=15.859000 !rg=3.59491 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,18,1:6)=(/ 13.129000,13.176000,13.265000,13.403000,13.618000,13.854000 /)
XPIZA_LKT(74,18,1:6)=(/ 0.656128,0.809103,0.990493,0.991301,0.969351,0.541808 /)
XCGA_LKT(74,18,1:6)=(/ 0.905960,0.869277,0.837077,0.832317,0.828980,0.933443 /)
XEXT_COEFF_550_LKT(74,18)=13.220000 !rg=3.59491 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,19,1:6)=(/ 10.897000,10.953000,10.994000,11.112000,11.241000,11.441000 /)
XPIZA_LKT(74,19,1:6)=(/ 0.642994,0.793320,0.989371,0.990443,0.967340,0.543916 /)
XCGA_LKT(74,19,1:6)=(/ 0.909693,0.872363,0.838753,0.834010,0.831183,0.934303 /)
XEXT_COEFF_550_LKT(74,19)=10.963000 !rg=3.59491 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,20,1:6)=(/ 9.052400,9.080800,9.117400,9.210400,9.330100,9.458100 /)
XPIZA_LKT(74,20,1:6)=(/ 0.630800,0.776241,0.987938,0.989380,0.965015,0.545774 /)
XCGA_LKT(74,20,1:6)=(/ 0.913513,0.875220,0.839813,0.835790,0.834057,0.934993 /)
XEXT_COEFF_550_LKT(74,20)=9.085100 !rg=3.59491 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,1,1:6)=(/ 162.320000,165.330000,167.660000,170.600000,217.010000,187.260000 /)
XPIZA_LKT(75,1,1:6)=(/ 0.763211,0.901168,0.996725,0.995435,0.993861,0.469403 /)
XCGA_LKT(75,1,1:6)=(/ 0.883553,0.836117,0.807567,0.757897,0.786523,0.877713 /)
XEXT_COEFF_550_LKT(75,1)=163.890000 !rg=3.89457 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,2,1:6)=(/ 153.820000,156.950000,160.530000,165.570000,185.380000,180.890000 /)
XPIZA_LKT(75,2,1:6)=(/ 0.766276,0.897673,0.996436,0.997383,0.993113,0.478388 /)
XCGA_LKT(75,2,1:6)=(/ 0.882423,0.842247,0.799820,0.780173,0.760793,0.881403 /)
XEXT_COEFF_550_LKT(75,2)=159.720000 !rg=3.89457 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,3,1:6)=(/ 143.280000,145.790000,147.970000,154.160000,162.620000,168.080000 /)
XPIZA_LKT(75,3,1:6)=(/ 0.774493,0.894795,0.996138,0.995260,0.992177,0.486660 /)
XCGA_LKT(75,3,1:6)=(/ 0.879707,0.841187,0.803007,0.776337,0.754883,0.885117 /)
XEXT_COEFF_550_LKT(75,3)=144.890000 !rg=3.89457 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,4,1:6)=(/ 129.270000,131.600000,133.240000,138.110000,144.300000,151.540000 /)
XPIZA_LKT(75,4,1:6)=(/ 0.778821,0.894117,0.995783,0.996481,0.991346,0.493390 /)
XCGA_LKT(75,4,1:6)=(/ 0.878417,0.842303,0.805670,0.782513,0.754843,0.889160 /)
XEXT_COEFF_550_LKT(75,4)=131.090000 !rg=3.89457 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,5,1:6)=(/ 114.400000,116.160000,116.880000,122.160000,126.140000,133.430000 /)
XPIZA_LKT(75,5,1:6)=(/ 0.780403,0.893268,0.995412,0.996292,0.990552,0.498914 /)
XCGA_LKT(75,5,1:6)=(/ 0.878087,0.844697,0.809357,0.784560,0.760520,0.893813 /)
XEXT_COEFF_550_LKT(75,5)=116.750000 !rg=3.89457 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,6,1:6)=(/ 99.512000,100.420000,102.460000,105.410000,110.580000,115.500000 /)
XPIZA_LKT(75,6,1:6)=(/ 0.778987,0.893797,0.995127,0.995667,0.989084,0.503522 /)
XCGA_LKT(75,6,1:6)=(/ 0.878510,0.846567,0.809273,0.790397,0.763113,0.898800 /)
XEXT_COEFF_550_LKT(75,6)=101.530000 !rg=3.89457 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,7,1:6)=(/ 85.744000,86.904000,88.038000,90.537000,93.955000,98.470000 /)
XPIZA_LKT(75,7,1:6)=(/ 0.775684,0.893128,0.994561,0.995275,0.987870,0.507582 /)
XCGA_LKT(75,7,1:6)=(/ 0.879580,0.847090,0.811950,0.795950,0.771710,0.904027 /)
XEXT_COEFF_550_LKT(75,7)=86.890000 !rg=3.89457 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,8,1:6)=(/ 73.141000,73.642000,75.413000,76.385000,80.238000,83.241000 /)
XPIZA_LKT(75,8,1:6)=(/ 0.771238,0.891922,0.994363,0.994980,0.985200,0.511370 /)
XCGA_LKT(75,8,1:6)=(/ 0.880400,0.849277,0.814960,0.802130,0.776663,0.908963 /)
XEXT_COEFF_550_LKT(75,8)=74.341000 !rg=3.89457 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,9,1:6)=(/ 61.857000,62.601000,63.063000,64.654000,66.892000,69.777000 /)
XPIZA_LKT(75,9,1:6)=(/ 0.763196,0.888864,0.994475,0.994612,0.983591,0.515067 /)
XCGA_LKT(75,9,1:6)=(/ 0.881750,0.851060,0.819633,0.805000,0.784990,0.913577 /)
XEXT_COEFF_550_LKT(75,9)=62.730000 !rg=3.89457 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,10,1:6)=(/ 52.234000,52.635000,53.184000,54.376000,56.370000,58.286000 /)
XPIZA_LKT(75,10,1:6)=(/ 0.754619,0.885663,0.994279,0.994169,0.981417,0.518695 /)
XCGA_LKT(75,10,1:6)=(/ 0.883910,0.852703,0.822670,0.810600,0.788573,0.917600 /)
XEXT_COEFF_550_LKT(75,10)=52.823000 !rg=3.89457 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,11,1:6)=(/ 43.770000,44.123000,44.651000,45.346000,47.008000,48.400000 /)
XPIZA_LKT(75,11,1:6)=(/ 0.744337,0.880239,0.994120,0.993757,0.979421,0.522307 /)
XCGA_LKT(75,11,1:6)=(/ 0.885983,0.854707,0.823857,0.814077,0.798237,0.921173 /)
XEXT_COEFF_550_LKT(75,11)=44.172000 !rg=3.89457 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,12,1:6)=(/ 36.541000,36.838000,37.188000,37.881000,39.087000,40.072000 /)
XPIZA_LKT(75,12,1:6)=(/ 0.732132,0.873714,0.993797,0.993592,0.977267,0.525836 /)
XCGA_LKT(75,12,1:6)=(/ 0.888397,0.856640,0.826810,0.817030,0.801363,0.924220 /)
XEXT_COEFF_550_LKT(75,12)=36.938000 !rg=3.89457 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,13,1:6)=(/ 30.419000,30.615000,30.991000,31.441000,32.214000,33.108000 /)
XPIZA_LKT(75,13,1:6)=(/ 0.719150,0.865421,0.993554,0.993303,0.975886,0.529237 /)
XCGA_LKT(75,13,1:6)=(/ 0.890927,0.858400,0.828927,0.820253,0.808960,0.926773 /)
XEXT_COEFF_550_LKT(75,13)=30.670000 !rg=3.89457 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,14,1:6)=(/ 25.372000,25.490000,25.722000,26.143000,26.792000,27.367000 /)
XPIZA_LKT(75,14,1:6)=(/ 0.705876,0.855804,0.993094,0.992984,0.973949,0.532423 /)
XCGA_LKT(75,14,1:6)=(/ 0.893903,0.860847,0.831940,0.824297,0.813710,0.928847 /)
XEXT_COEFF_550_LKT(75,14)=25.640000 !rg=3.89457 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,15,1:6)=(/ 21.099000,21.243000,21.384000,21.561000,22.186000,22.612000 /)
XPIZA_LKT(75,15,1:6)=(/ 0.691401,0.844810,0.992498,0.992734,0.972647,0.535390 /)
XCGA_LKT(75,15,1:6)=(/ 0.896980,0.863213,0.832767,0.826920,0.817997,0.930547 /)
XEXT_COEFF_550_LKT(75,15)=21.258000 !rg=3.89457 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,16,1:6)=(/ 17.518000,17.585000,17.737000,17.897000,18.308000,18.647000 /)
XPIZA_LKT(75,16,1:6)=(/ 0.677238,0.831521,0.991863,0.992220,0.971857,0.538156 /)
XCGA_LKT(75,16,1:6)=(/ 0.900427,0.865533,0.835103,0.828793,0.823557,0.931960 /)
XEXT_COEFF_550_LKT(75,16)=17.612000 !rg=3.89457 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,17,1:6)=(/ 14.572000,14.615000,14.730000,14.905000,15.122000,15.411000 /)
XPIZA_LKT(75,17,1:6)=(/ 0.663480,0.817481,0.991039,0.991595,0.970138,0.540629 /)
XCGA_LKT(75,17,1:6)=(/ 0.904103,0.867807,0.836583,0.831213,0.825183,0.933080 /)
XEXT_COEFF_550_LKT(75,17)=14.669000 !rg=3.89457 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,18,1:6)=(/ 12.107000,12.173000,12.226000,12.305000,12.548000,12.744000 /)
XPIZA_LKT(75,18,1:6)=(/ 0.649960,0.802478,0.989981,0.990980,0.968283,0.542846 /)
XCGA_LKT(75,18,1:6)=(/ 0.907610,0.870890,0.837563,0.833653,0.828533,0.933983 /)
XEXT_COEFF_550_LKT(75,18)=12.174000 !rg=3.89457 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,19,1:6)=(/ 10.048000,10.086000,10.143000,10.212000,10.369000,10.528000 /)
XPIZA_LKT(75,19,1:6)=(/ 0.637051,0.785657,0.988759,0.989970,0.966404,0.544844 /)
XCGA_LKT(75,19,1:6)=(/ 0.911443,0.873570,0.839283,0.834600,0.832840,0.934730 /)
XEXT_COEFF_550_LKT(75,19)=10.097000 !rg=3.89457 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,20,1:6)=(/ 8.355300,8.372600,8.418000,8.493200,8.580300,8.705900 /)
XPIZA_LKT(75,20,1:6)=(/ 0.625542,0.768440,0.987215,0.988929,0.963775,0.546596 /)
XCGA_LKT(75,20,1:6)=(/ 0.915413,0.876410,0.840503,0.837110,0.834537,0.935327 /)
XEXT_COEFF_550_LKT(75,20)=8.391800 !rg=3.89457 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,1,1:6)=(/ 149.090000,151.670000,151.040000,146.880000,188.820000,175.300000 /)
XPIZA_LKT(76,1,1:6)=(/ 0.774738,0.894691,0.996312,0.997076,0.992083,0.484236 /)
XCGA_LKT(76,1,1:6)=(/ 0.881120,0.843497,0.800200,0.764077,0.774950,0.896713 /)
XEXT_COEFF_550_LKT(76,1)=151.890000 !rg=4.21921 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,2,1:6)=(/ 142.150000,143.520000,146.650000,151.850000,166.440000,167.230000 /)
XPIZA_LKT(76,2,1:6)=(/ 0.776020,0.889856,0.996134,0.996758,0.991834,0.486260 /)
XCGA_LKT(76,2,1:6)=(/ 0.880200,0.844670,0.805003,0.777377,0.766417,0.895990 /)
XEXT_COEFF_550_LKT(76,2)=146.090000 !rg=4.21921 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,3,1:6)=(/ 132.010000,133.370000,136.590000,140.080000,147.860000,154.110000 /)
XPIZA_LKT(76,3,1:6)=(/ 0.778965,0.889791,0.995839,0.996690,0.991585,0.489871 /)
XCGA_LKT(76,3,1:6)=(/ 0.879403,0.844983,0.806730,0.782833,0.757000,0.896227 /)
XEXT_COEFF_550_LKT(76,3)=134.230000 !rg=4.21921 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,4,1:6)=(/ 118.950000,120.550000,122.940000,126.690000,132.100000,138.520000 /)
XPIZA_LKT(76,4,1:6)=(/ 0.780173,0.890416,0.994818,0.995744,0.990214,0.494735 /)
XCGA_LKT(76,4,1:6)=(/ 0.878930,0.845417,0.807787,0.782893,0.756253,0.898127 /)
XEXT_COEFF_550_LKT(76,4)=121.400000 !rg=4.21921 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,5,1:6)=(/ 105.230000,106.990000,108.080000,112.110000,116.200000,121.880000 /)
XPIZA_LKT(76,5,1:6)=(/ 0.779749,0.892606,0.995131,0.995252,0.989049,0.499686 /)
XCGA_LKT(76,5,1:6)=(/ 0.878447,0.845233,0.811733,0.788930,0.766320,0.901253 /)
XEXT_COEFF_550_LKT(76,5)=106.590000 !rg=4.21921 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,6,1:6)=(/ 91.598000,92.977000,94.957000,96.581000,101.390000,105.540000 /)
XPIZA_LKT(76,6,1:6)=(/ 0.778019,0.892658,0.994639,0.995587,0.987919,0.504258 /)
XCGA_LKT(76,6,1:6)=(/ 0.878807,0.847417,0.810043,0.799160,0.761653,0.905023 /)
XEXT_COEFF_550_LKT(76,6)=93.061000 !rg=4.21921 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,7,1:6)=(/ 78.896000,79.901000,80.526000,83.110000,86.752000,90.031000 /)
XPIZA_LKT(76,7,1:6)=(/ 0.773075,0.891582,0.994825,0.994978,0.985289,0.508504 /)
XCGA_LKT(76,7,1:6)=(/ 0.880263,0.848520,0.816977,0.795953,0.769757,0.909197 /)
XEXT_COEFF_550_LKT(76,7)=80.895000 !rg=4.21921 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,8,1:6)=(/ 67.379000,67.980000,68.980000,70.645000,73.447000,76.165000 /)
XPIZA_LKT(76,8,1:6)=(/ 0.767373,0.890399,0.994386,0.994567,0.984487,0.512507 /)
XCGA_LKT(76,8,1:6)=(/ 0.881167,0.850287,0.818367,0.804750,0.777430,0.913237 /)
XEXT_COEFF_550_LKT(76,8)=67.982000 !rg=4.21921 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,9,1:6)=(/ 57.052000,57.588000,58.394000,59.823000,61.523000,63.894000 /)
XPIZA_LKT(76,9,1:6)=(/ 0.759593,0.887372,0.994357,0.993864,0.982226,0.516393 /)
XCGA_LKT(76,9,1:6)=(/ 0.882920,0.851153,0.821430,0.805777,0.790837,0.917070 /)
XEXT_COEFF_550_LKT(76,9)=57.721000 !rg=4.21921 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,10,1:6)=(/ 48.056000,48.567000,49.136000,50.021000,51.684000,53.410000 /)
XPIZA_LKT(76,10,1:6)=(/ 0.749877,0.883054,0.994182,0.994021,0.979907,0.520142 /)
XCGA_LKT(76,10,1:6)=(/ 0.884533,0.853430,0.823203,0.812017,0.793960,0.920450 /)
XEXT_COEFF_550_LKT(76,10)=48.659000 !rg=4.21921 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,11,1:6)=(/ 40.325000,40.641000,40.973000,41.721000,42.904000,44.382000 /)
XPIZA_LKT(76,11,1:6)=(/ 0.738797,0.877271,0.993969,0.993615,0.978811,0.523817 /)
XCGA_LKT(76,11,1:6)=(/ 0.887203,0.855647,0.825373,0.816157,0.800360,0.923473 /)
XEXT_COEFF_550_LKT(76,11)=40.953000 !rg=4.21921 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,12,1:6)=(/ 33.642000,33.916000,34.318000,34.949000,35.715000,36.767000 /)
XPIZA_LKT(76,12,1:6)=(/ 0.726135,0.870098,0.993826,0.993138,0.976582,0.527357 /)
XCGA_LKT(76,12,1:6)=(/ 0.889440,0.857283,0.828590,0.816910,0.806690,0.926063 /)
XEXT_COEFF_550_LKT(76,12)=33.940000 !rg=4.21921 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,13,1:6)=(/ 28.062000,28.219000,28.553000,29.056000,29.669000,30.396000 /)
XPIZA_LKT(76,13,1:6)=(/ 0.713107,0.861287,0.993305,0.992915,0.974862,0.530722 /)
XCGA_LKT(76,13,1:6)=(/ 0.892423,0.859763,0.831067,0.821913,0.810670,0.928243 /)
XEXT_COEFF_550_LKT(76,13)=28.399000 !rg=4.21921 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,14,1:6)=(/ 23.371000,23.533000,23.701000,23.981000,24.569000,25.138000 /)
XPIZA_LKT(76,14,1:6)=(/ 0.698988,0.851056,0.992903,0.992826,0.973775,0.533838 /)
XCGA_LKT(76,14,1:6)=(/ 0.895210,0.862370,0.832293,0.825723,0.815847,0.930023 /)
XEXT_COEFF_550_LKT(76,14)=23.535000 !rg=4.21921 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,15,1:6)=(/ 19.444000,19.555000,19.746000,19.889000,20.406000,20.779000 /)
XPIZA_LKT(76,15,1:6)=(/ 0.684934,0.838944,0.992336,0.992399,0.971907,0.536715 /)
XCGA_LKT(76,15,1:6)=(/ 0.898373,0.863820,0.834530,0.827640,0.820137,0.931483 /)
XEXT_COEFF_550_LKT(76,15)=19.613000 !rg=4.21921 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,16,1:6)=(/ 16.158000,16.232000,16.334000,16.500000,16.773000,17.143000 /)
XPIZA_LKT(76,16,1:6)=(/ 0.670827,0.825572,0.991501,0.992012,0.971075,0.539375 /)
XCGA_LKT(76,16,1:6)=(/ 0.902077,0.866660,0.835223,0.830423,0.824603,0.932703 /)
XEXT_COEFF_550_LKT(76,16)=16.280000 !rg=4.21921 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,17,1:6)=(/ 13.429000,13.498000,13.570000,13.695000,13.888000,14.174000 /)
XPIZA_LKT(76,17,1:6)=(/ 0.657075,0.811072,0.990666,0.991280,0.969643,0.541737 /)
XCGA_LKT(76,17,1:6)=(/ 0.905613,0.869387,0.836957,0.832537,0.828547,0.933670 /)
XEXT_COEFF_550_LKT(76,17)=13.515000 !rg=4.21921 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,18,1:6)=(/ 11.166000,11.205000,11.284000,11.348000,11.546000,11.725000 /)
XPIZA_LKT(76,18,1:6)=(/ 0.644044,0.794902,0.989507,0.990538,0.967211,0.543841 /)
XCGA_LKT(76,18,1:6)=(/ 0.909310,0.871737,0.838927,0.833963,0.830513,0.934453 /)
XEXT_COEFF_550_LKT(76,18)=11.235000 !rg=4.21921 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,19,1:6)=(/ 9.275200,9.302500,9.348800,9.432500,9.529100,9.689000 /)
XPIZA_LKT(76,19,1:6)=(/ 0.631778,0.778149,0.988097,0.989616,0.965246,0.545728 /)
XCGA_LKT(76,19,1:6)=(/ 0.913203,0.874847,0.839610,0.836063,0.833863,0.935100 /)
XEXT_COEFF_550_LKT(76,19)=9.320600 !rg=4.21921 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,20,1:6)=(/ 7.703800,7.729700,7.761500,7.812700,7.893000,8.014600 /)
XPIZA_LKT(76,20,1:6)=(/ 0.620238,0.760734,0.986484,0.988387,0.962684,0.547375 /)
XCGA_LKT(76,20,1:6)=(/ 0.917067,0.878147,0.840790,0.837670,0.836547,0.935617 /)
XEXT_COEFF_550_LKT(76,20)=7.737700 !rg=4.21921 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,1,1:6)=(/ 136.700000,137.910000,143.290000,149.370000,151.110000,163.080000 /)
XPIZA_LKT(77,1,1:6)=(/ 0.782851,0.883533,0.996199,0.997048,0.990521,0.496339 /)
XCGA_LKT(77,1,1:6)=(/ 0.878547,0.847507,0.811413,0.790070,0.739320,0.908873 /)
XEXT_COEFF_550_LKT(77,1)=139.650000 !rg=4.57091 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,2,1:6)=(/ 130.530000,132.950000,134.710000,138.480000,145.390000,154.120000 /)
XPIZA_LKT(77,2,1:6)=(/ 0.782648,0.882875,0.995536,0.996867,0.990714,0.493176 /)
XCGA_LKT(77,2,1:6)=(/ 0.877970,0.847673,0.805543,0.788393,0.749703,0.906033 /)
XEXT_COEFF_550_LKT(77,2)=133.210000 !rg=4.57091 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,3,1:6)=(/ 121.330000,123.890000,125.580000,128.100000,137.300000,141.300000 /)
XPIZA_LKT(77,3,1:6)=(/ 0.782114,0.885448,0.995605,0.995993,0.989737,0.493623 /)
XCGA_LKT(77,3,1:6)=(/ 0.878140,0.846887,0.808210,0.788787,0.764043,0.904653 /)
XEXT_COEFF_550_LKT(77,3)=124.010000 !rg=4.57091 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,4,1:6)=(/ 109.580000,111.660000,113.280000,115.410000,122.820000,126.740000 /)
XPIZA_LKT(77,4,1:6)=(/ 0.781571,0.888921,0.995142,0.995988,0.988974,0.496982 /)
XCGA_LKT(77,4,1:6)=(/ 0.878697,0.847293,0.809490,0.791713,0.764947,0.905340 /)
XEXT_COEFF_550_LKT(77,4)=111.360000 !rg=4.57091 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,5,1:6)=(/ 97.019000,97.956000,99.502000,102.120000,106.380000,111.450000 /)
XPIZA_LKT(77,5,1:6)=(/ 0.779642,0.891331,0.994637,0.995757,0.988447,0.501283 /)
XCGA_LKT(77,5,1:6)=(/ 0.879243,0.847653,0.812620,0.792740,0.767777,0.907423 /)
XEXT_COEFF_550_LKT(77,5)=98.696000 !rg=4.57091 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,6,1:6)=(/ 84.335000,85.027000,86.524000,89.232000,93.338000,96.512000 /)
XPIZA_LKT(77,6,1:6)=(/ 0.776138,0.892115,0.994529,0.995252,0.986668,0.505661 /)
XCGA_LKT(77,6,1:6)=(/ 0.879457,0.848877,0.815710,0.798587,0.770377,0.910260 /)
XEXT_COEFF_550_LKT(77,6)=85.919000 !rg=4.57091 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,7,1:6)=(/ 72.779000,73.595000,74.762000,75.950000,79.541000,82.375000 /)
XPIZA_LKT(77,7,1:6)=(/ 0.770957,0.890937,0.994331,0.994922,0.985020,0.509927 /)
XCGA_LKT(77,7,1:6)=(/ 0.880907,0.850623,0.815663,0.803567,0.780450,0.913580 /)
XEXT_COEFF_550_LKT(77,7)=73.636000 !rg=4.57091 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,8,1:6)=(/ 61.941000,62.570000,63.538000,65.180000,67.295000,69.731000 /)
XPIZA_LKT(77,8,1:6)=(/ 0.763568,0.888629,0.994448,0.994077,0.983280,0.514004 /)
XCGA_LKT(77,8,1:6)=(/ 0.882030,0.850903,0.819863,0.804017,0.785167,0.916877 /)
XEXT_COEFF_550_LKT(77,8)=62.754000 !rg=4.57091 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,9,1:6)=(/ 52.494000,52.946000,53.544000,54.494000,56.618000,58.537000 /)
XPIZA_LKT(77,9,1:6)=(/ 0.755011,0.885600,0.994341,0.994145,0.981123,0.517962 /)
XCGA_LKT(77,9,1:6)=(/ 0.883800,0.852873,0.823147,0.809730,0.792010,0.920060 /)
XEXT_COEFF_550_LKT(77,9)=53.156000 !rg=4.57091 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,10,1:6)=(/ 44.320000,44.600000,45.263000,46.069000,47.450000,48.963000 /)
XPIZA_LKT(77,10,1:6)=(/ 0.745153,0.880443,0.993959,0.993964,0.979134,0.521750 /)
XCGA_LKT(77,10,1:6)=(/ 0.885920,0.854517,0.824517,0.814797,0.799997,0.922890 /)
XEXT_COEFF_550_LKT(77,10)=44.935000 !rg=4.57091 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,11,1:6)=(/ 37.139000,37.471000,37.664000,38.306000,39.360000,40.712000 /)
XPIZA_LKT(77,11,1:6)=(/ 0.733049,0.874469,0.993983,0.993733,0.977560,0.525425 /)
XCGA_LKT(77,11,1:6)=(/ 0.888067,0.857140,0.827847,0.818577,0.805733,0.925447 /)
XEXT_COEFF_550_LKT(77,11)=37.576000 !rg=4.57091 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,12,1:6)=(/ 31.076000,31.231000,31.572000,32.115000,32.956000,33.747000 /)
XPIZA_LKT(77,12,1:6)=(/ 0.720624,0.866527,0.993541,0.993329,0.975526,0.528926 /)
XCGA_LKT(77,12,1:6)=(/ 0.890993,0.858630,0.830067,0.821710,0.809640,0.927650 /)
XEXT_COEFF_550_LKT(77,12)=31.392000 !rg=4.57091 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,13,1:6)=(/ 25.862000,26.077000,26.317000,26.628000,27.219000,27.913000 /)
XPIZA_LKT(77,13,1:6)=(/ 0.706458,0.857097,0.993103,0.992910,0.974187,0.532219 /)
XCGA_LKT(77,13,1:6)=(/ 0.893657,0.861283,0.830320,0.823620,0.813977,0.929513 /)
XEXT_COEFF_550_LKT(77,13)=26.068000 !rg=4.57091 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,14,1:6)=(/ 21.576000,21.654000,21.938000,22.130000,22.576000,23.096000 /)
XPIZA_LKT(77,14,1:6)=(/ 0.692825,0.845753,0.992637,0.992647,0.973255,0.535244 /)
XCGA_LKT(77,14,1:6)=(/ 0.896903,0.863010,0.833573,0.826470,0.819443,0.931037 /)
XEXT_COEFF_550_LKT(77,14)=21.753000 !rg=4.57091 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,15,1:6)=(/ 17.942000,18.032000,18.138000,18.389000,18.718000,19.100000 /)
XPIZA_LKT(77,15,1:6)=(/ 0.678394,0.833352,0.992029,0.992125,0.971672,0.538012 /)
XCGA_LKT(77,15,1:6)=(/ 0.900150,0.865370,0.834760,0.829497,0.823087,0.932293 /)
XEXT_COEFF_550_LKT(77,15)=18.070000 !rg=4.57091 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,16,1:6)=(/ 14.896000,14.973000,15.040000,15.183000,15.422000,15.764000 /)
XPIZA_LKT(77,16,1:6)=(/ 0.664302,0.819421,0.991195,0.991696,0.970111,0.540556 /)
XCGA_LKT(77,16,1:6)=(/ 0.903650,0.868093,0.836807,0.831313,0.827207,0.933343 /)
XEXT_COEFF_550_LKT(77,16)=15.003000 !rg=4.57091 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,17,1:6)=(/ 12.390000,12.431000,12.524000,12.618000,12.806000,13.038000 /)
XPIZA_LKT(77,17,1:6)=(/ 0.651001,0.803893,0.990170,0.990984,0.968619,0.542801 /)
XCGA_LKT(77,17,1:6)=(/ 0.907363,0.870467,0.838197,0.832767,0.830533,0.934180 /)
XEXT_COEFF_550_LKT(77,17)=12.484000 !rg=4.57091 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,18,1:6)=(/ 10.303000,10.342000,10.390000,10.507000,10.627000,10.790000 /)
XPIZA_LKT(77,18,1:6)=(/ 0.638212,0.787813,0.988948,0.990159,0.966495,0.544791 /)
XCGA_LKT(77,18,1:6)=(/ 0.911127,0.873330,0.839400,0.835333,0.832583,0.934857 /)
XEXT_COEFF_550_LKT(77,18)=10.358000 !rg=4.57091 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,19,1:6)=(/ 8.554900,8.582000,8.616100,8.681500,8.771500,8.918400 /)
XPIZA_LKT(77,19,1:6)=(/ 0.626213,0.770517,0.987430,0.989001,0.963952,0.546566 /)
XCGA_LKT(77,19,1:6)=(/ 0.915010,0.876440,0.840530,0.836377,0.835380,0.935417 /)
XEXT_COEFF_550_LKT(77,19)=8.594200 !rg=4.57091 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,20,1:6)=(/ 7.108000,7.126200,7.161400,7.205500,7.283900,7.379200 /)
XPIZA_LKT(77,20,1:6)=(/ 0.615273,0.752664,0.985654,0.987759,0.961218,0.548110 /)
XCGA_LKT(77,20,1:6)=(/ 0.000000,0.879567,0.841593,0.838080,0.838217,0.935863 /)
XEXT_COEFF_550_LKT(77,20)=7.142400 !rg=4.57091 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,1,1:6)=(/ 126.160000,127.950000,125.830000,140.500000,119.450000,149.660000 /)
XPIZA_LKT(78,1,1:6)=(/ 0.787336,0.871082,0.995798,0.992080,0.990054,0.501783 /)
XCGA_LKT(78,1,1:6)=(/ 0.877357,0.853087,0.811857,0.793993,0.766273,0.914903 /)
XEXT_COEFF_550_LKT(78,1)=128.740000 !rg=4.95192 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,2,1:6)=(/ 120.560000,121.470000,124.170000,128.040000,134.740000,141.250000 /)
XPIZA_LKT(78,2,1:6)=(/ 0.786633,0.879881,0.995583,0.996679,0.988101,0.497536 /)
XCGA_LKT(78,2,1:6)=(/ 0.877643,0.850760,0.813587,0.792747,0.753553,0.912140 /)
XEXT_COEFF_550_LKT(78,2)=122.020000 !rg=4.95192 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,3,1:6)=(/ 111.570000,113.060000,114.410000,119.080000,125.130000,129.430000 /)
XPIZA_LKT(78,3,1:6)=(/ 0.784063,0.886141,0.995070,0.996117,0.989701,0.497128 /)
XCGA_LKT(78,3,1:6)=(/ 0.878203,0.848460,0.810113,0.789253,0.754780,0.910720 /)
XEXT_COEFF_550_LKT(78,3)=113.420000 !rg=4.95192 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,4,1:6)=(/ 100.790000,102.160000,102.980000,107.130000,112.030000,116.010000 /)
XPIZA_LKT(78,4,1:6)=(/ 0.782008,0.889780,0.994867,0.995712,0.989082,0.499643 /)
XCGA_LKT(78,4,1:6)=(/ 0.878713,0.848303,0.812700,0.792763,0.762770,0.910980 /)
XEXT_COEFF_550_LKT(78,4)=102.590000 !rg=4.95192 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,5,1:6)=(/ 89.231000,90.764000,91.952000,93.446000,98.582000,101.980000 /)
XPIZA_LKT(78,5,1:6)=(/ 0.778632,0.890125,0.994654,0.994939,0.987332,0.503409 /)
XCGA_LKT(78,5,1:6)=(/ 0.879257,0.848800,0.814010,0.799000,0.775023,0.912433 /)
XEXT_COEFF_550_LKT(78,5)=90.676000 !rg=4.95192 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,6,1:6)=(/ 77.958000,78.475000,80.338000,82.166000,85.194000,88.324000 /)
XPIZA_LKT(78,6,1:6)=(/ 0.774094,0.890355,0.994386,0.994689,0.985066,0.507524 /)
XCGA_LKT(78,6,1:6)=(/ 0.880733,0.850350,0.815060,0.799443,0.774923,0.914607 /)
XEXT_COEFF_550_LKT(78,6)=79.229000 !rg=4.95192 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,7,1:6)=(/ 66.958000,67.666000,68.270000,69.909000,72.552000,75.417000 /)
XPIZA_LKT(78,7,1:6)=(/ 0.767765,0.889732,0.994367,0.994593,0.984628,0.511705 /)
XCGA_LKT(78,7,1:6)=(/ 0.881697,0.851180,0.818597,0.806577,0.782207,0.917260 /)
XEXT_COEFF_550_LKT(78,7)=68.341000 !rg=4.95192 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,8,1:6)=(/ 57.199000,57.623000,58.527000,59.749000,61.716000,63.875000 /)
XPIZA_LKT(78,8,1:6)=(/ 0.760217,0.887512,0.994299,0.994267,0.981761,0.515760 /)
XCGA_LKT(78,8,1:6)=(/ 0.883150,0.852913,0.822607,0.809903,0.789830,0.919957 /)
XEXT_COEFF_550_LKT(78,8)=57.974000 !rg=4.95192 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,9,1:6)=(/ 48.370000,48.913000,49.480000,50.012000,52.215000,53.652000 /)
XPIZA_LKT(78,9,1:6)=(/ 0.750305,0.883112,0.994236,0.993880,0.979433,0.519711 /)
XCGA_LKT(78,9,1:6)=(/ 0.884737,0.853597,0.823990,0.814327,0.795920,0.922600 /)
XEXT_COEFF_550_LKT(78,9)=49.111000 !rg=4.95192 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,10,1:6)=(/ 40.837000,41.181000,41.587000,42.181000,43.658000,44.904000 /)
XPIZA_LKT(78,10,1:6)=(/ 0.739743,0.878216,0.993954,0.993457,0.978484,0.523466 /)
XCGA_LKT(78,10,1:6)=(/ 0.886743,0.855940,0.826077,0.815680,0.801350,0.924973 /)
XEXT_COEFF_550_LKT(78,10)=41.156000 !rg=4.95192 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,11,1:6)=(/ 34.316000,34.496000,34.911000,35.409000,36.339000,37.358000 /)
XPIZA_LKT(78,11,1:6)=(/ 0.727642,0.870669,0.993651,0.993453,0.976409,0.527091 /)
XCGA_LKT(78,11,1:6)=(/ 0.889687,0.857340,0.829087,0.819143,0.808643,0.927137 /)
XEXT_COEFF_550_LKT(78,11)=34.679000 !rg=4.95192 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,12,1:6)=(/ 28.604000,28.828000,29.060000,29.415000,30.248000,30.983000 /)
XPIZA_LKT(78,12,1:6)=(/ 0.713883,0.862560,0.993467,0.993289,0.975144,0.530514 /)
XCGA_LKT(78,12,1:6)=(/ 0.891983,0.860167,0.830693,0.823267,0.811720,0.929010 /)
XEXT_COEFF_550_LKT(78,12)=28.794000 !rg=4.95192 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,13,1:6)=(/ 23.859000,23.966000,24.235000,24.504000,24.981000,25.641000 /)
XPIZA_LKT(78,13,1:6)=(/ 0.700215,0.852211,0.992886,0.992977,0.974035,0.533711 /)
XCGA_LKT(78,13,1:6)=(/ 0.895137,0.861800,0.832023,0.825733,0.817987,0.930603 /)
XEXT_COEFF_550_LKT(78,13)=24.088000 !rg=4.95192 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,14,1:6)=(/ 19.886000,19.975000,20.122000,20.382000,20.832000,21.226000 /)
XPIZA_LKT(78,14,1:6)=(/ 0.686156,0.840420,0.992414,0.992422,0.972120,0.536622 /)
XCGA_LKT(78,14,1:6)=(/ 0.898313,0.864110,0.834373,0.828170,0.819597,0.931910 /)
XEXT_COEFF_550_LKT(78,14)=20.005000 !rg=4.95192 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,15,1:6)=(/ 16.545000,16.628000,16.692000,16.870000,17.218000,17.560000 /)
XPIZA_LKT(78,15,1:6)=(/ 0.671905,0.827246,0.991602,0.991921,0.970774,0.539272 /)
XCGA_LKT(78,15,1:6)=(/ 0.901650,0.866587,0.835460,0.830173,0.824343,0.932993 /)
XEXT_COEFF_550_LKT(78,15)=16.645000 !rg=4.95192 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,16,1:6)=(/ 13.747000,13.806000,13.914000,14.032000,14.252000,14.499000 /)
XPIZA_LKT(78,16,1:6)=(/ 0.658101,0.812657,0.990806,0.991385,0.969423,0.541694 /)
XCGA_LKT(78,16,1:6)=(/ 0.905460,0.869090,0.837903,0.832353,0.828360,0.933900 /)
XEXT_COEFF_550_LKT(78,16)=13.846000 !rg=4.95192 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,17,1:6)=(/ 11.425000,11.478000,11.535000,11.648000,11.765000,11.996000 /)
XPIZA_LKT(78,17,1:6)=(/ 0.644863,0.796932,0.989660,0.990609,0.967921,0.543819 /)
XCGA_LKT(78,17,1:6)=(/ 0.909083,0.871860,0.838277,0.834667,0.831483,0.934620 /)
XEXT_COEFF_550_LKT(78,17)=11.490000 !rg=4.95192 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,18,1:6)=(/ 9.503800,9.536100,9.579300,9.645700,9.779600,9.930000 /)
XPIZA_LKT(78,18,1:6)=(/ 0.632474,0.780111,0.988285,0.989706,0.965481,0.545693 /)
XCGA_LKT(78,18,1:6)=(/ 0.912893,0.874663,0.839587,0.836050,0.833773,0.935207 /)
XEXT_COEFF_550_LKT(78,18)=9.544400 !rg=4.95192 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,19,1:6)=(/ 7.894300,7.917200,7.962400,8.020900,8.110400,8.210300 /)
XPIZA_LKT(78,19,1:6)=(/ 0.621058,0.762711,0.986699,0.988476,0.962979,0.547359 /)
XCGA_LKT(78,19,1:6)=(/ 0.916803,0.877747,0.841377,0.837610,0.836697,0.935690 /)
XEXT_COEFF_550_LKT(78,19)=7.926300 !rg=4.95192 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,20,1:6)=(/ 6.557200,6.577700,6.594700,6.650900,6.702600,6.795100 /)
XPIZA_LKT(78,20,1:6)=(/ 0.610406,0.744823,0.984716,0.987131,0.959832,0.548802 /)
XCGA_LKT(78,20,1:6)=(/ 0.000000,0.881147,0.841843,0.839387,0.838827,0.936077 /)
XEXT_COEFF_550_LKT(78,20)=6.583400 !rg=4.95192 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,1,1:6)=(/ 115.280000,118.760000,120.020000,117.200000,117.860000,135.700000 /)
XPIZA_LKT(79,1,1:6)=(/ 0.790928,0.871632,0.994864,0.996432,0.989765,0.501156 /)
XCGA_LKT(79,1,1:6)=(/ 0.875433,0.853543,0.806667,0.790380,0.761807,0.916737 /)
XEXT_COEFF_550_LKT(79,1)=116.430000 !rg=5.36469 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,2,1:6)=(/ 111.580000,112.240000,115.500000,120.200000,121.130000,128.910000 /)
XPIZA_LKT(79,2,1:6)=(/ 0.788940,0.878718,0.995007,0.995591,0.989070,0.499692 /)
XCGA_LKT(79,2,1:6)=(/ 0.877030,0.852307,0.808483,0.786340,0.765570,0.915573 /)
XEXT_COEFF_550_LKT(79,2)=114.000000 !rg=5.36469 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,3,1:6)=(/ 103.000000,104.720000,105.160000,109.890000,112.310000,118.440000 /)
XPIZA_LKT(79,3,1:6)=(/ 0.785465,0.885532,0.994349,0.995329,0.989459,0.500179 /)
XCGA_LKT(79,3,1:6)=(/ 0.877080,0.849947,0.813603,0.792427,0.779767,0.915010 /)
XEXT_COEFF_550_LKT(79,3)=104.790000 !rg=5.36469 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,4,1:6)=(/ 93.000000,94.157000,94.809000,98.464000,101.250000,106.210000 /)
XPIZA_LKT(79,4,1:6)=(/ 0.781222,0.888728,0.994464,0.995565,0.988560,0.502440 /)
XCGA_LKT(79,4,1:6)=(/ 0.878687,0.850107,0.815683,0.796873,0.780667,0.915320 /)
XEXT_COEFF_550_LKT(79,4)=94.675000 !rg=5.36469 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,5,1:6)=(/ 82.206000,83.009000,83.815000,86.527000,89.971000,93.363000 /)
XPIZA_LKT(79,5,1:6)=(/ 0.777160,0.891127,0.994173,0.994955,0.986974,0.505828 /)
XCGA_LKT(79,5,1:6)=(/ 0.879980,0.849737,0.816193,0.800020,0.774920,0.916457 /)
XEXT_COEFF_550_LKT(79,5)=83.474000 !rg=5.36469 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,6,1:6)=(/ 71.531000,72.375000,73.750000,75.120000,77.446000,80.874000 /)
XPIZA_LKT(79,6,1:6)=(/ 0.771204,0.890615,0.994305,0.994526,0.985069,0.509686 /)
XCGA_LKT(79,6,1:6)=(/ 0.880837,0.851173,0.817210,0.805377,0.784927,0.918177 /)
XEXT_COEFF_550_LKT(79,6)=72.751000 !rg=5.36469 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,7,1:6)=(/ 61.647000,62.325000,62.737000,64.169000,66.179000,69.082000 /)
XPIZA_LKT(79,7,1:6)=(/ 0.764034,0.888740,0.994350,0.994485,0.983616,0.513725 /)
XCGA_LKT(79,7,1:6)=(/ 0.882190,0.852570,0.822060,0.809930,0.791453,0.920323 /)
XEXT_COEFF_550_LKT(79,7)=62.657000 !rg=5.36469 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,8,1:6)=(/ 52.633000,53.161000,54.077000,54.721000,56.610000,58.539000 /)
XPIZA_LKT(79,8,1:6)=(/ 0.755252,0.885542,0.993943,0.994100,0.980962,0.517690 /)
XCGA_LKT(79,8,1:6)=(/ 0.883800,0.854310,0.821747,0.812853,0.793763,0.922543 /)
XEXT_COEFF_550_LKT(79,8)=53.147000 !rg=5.36469 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,9,1:6)=(/ 44.665000,44.892000,45.369000,46.314000,47.637000,49.197000 /)
XPIZA_LKT(79,9,1:6)=(/ 0.745706,0.881021,0.993917,0.993438,0.979151,0.521570 /)
XCGA_LKT(79,9,1:6)=(/ 0.886130,0.854877,0.825080,0.815173,0.799440,0.924747 /)
XEXT_COEFF_550_LKT(79,9)=45.194000 !rg=5.36469 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,10,1:6)=(/ 37.640000,37.887000,38.405000,38.827000,40.135000,41.197000 /)
XPIZA_LKT(79,10,1:6)=(/ 0.734069,0.874721,0.993992,0.993473,0.977057,0.525249 /)
XCGA_LKT(79,10,1:6)=(/ 0.888043,0.856470,0.828293,0.818787,0.804577,0.926747 /)
XEXT_COEFF_550_LKT(79,10)=38.083000 !rg=5.36469 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,11,1:6)=(/ 31.559000,31.829000,31.928000,32.535000,33.493000,34.292000 /)
XPIZA_LKT(79,11,1:6)=(/ 0.721104,0.867329,0.993651,0.993214,0.975065,0.528781 /)
XCGA_LKT(79,11,1:6)=(/ 0.890757,0.858817,0.830293,0.819727,0.808800,0.928577 /)
XEXT_COEFF_550_LKT(79,11)=31.961000 !rg=5.36469 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,12,1:6)=(/ 26.411000,26.534000,26.909000,27.150000,27.719000,28.455000 /)
XPIZA_LKT(79,12,1:6)=(/ 0.708123,0.857971,0.993250,0.992997,0.974344,0.532100 /)
XCGA_LKT(79,12,1:6)=(/ 0.893603,0.860720,0.831920,0.824040,0.815273,0.930173 /)
XEXT_COEFF_550_LKT(79,12)=26.638000 !rg=5.36469 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,13,1:6)=(/ 21.997000,22.110000,22.294000,22.573000,23.031000,23.559000 /)
XPIZA_LKT(79,13,1:6)=(/ 0.693549,0.847264,0.992742,0.992680,0.972931,0.535178 /)
XCGA_LKT(79,13,1:6)=(/ 0.896533,0.862943,0.833380,0.827630,0.818837,0.931540 /)
XEXT_COEFF_550_LKT(79,13)=22.169000 !rg=5.36469 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,14,1:6)=(/ 18.329000,18.415000,18.567000,18.758000,19.115000,19.511000 /)
XPIZA_LKT(79,14,1:6)=(/ 0.679425,0.834724,0.992028,0.992341,0.971822,0.537966 /)
XCGA_LKT(79,14,1:6)=(/ 0.899870,0.865123,0.834780,0.828973,0.823777,0.932660 /)
XEXT_COEFF_550_LKT(79,14)=18.444000 !rg=5.36469 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,15,1:6)=(/ 15.268000,15.313000,15.451000,15.602000,15.879000,16.148000 /)
XPIZA_LKT(79,15,1:6)=(/ 0.665668,0.820831,0.991259,0.991757,0.970278,0.540489 /)
XCGA_LKT(79,15,1:6)=(/ 0.903440,0.867550,0.836237,0.831257,0.827623,0.933593 /)
XEXT_COEFF_550_LKT(79,15)=15.387000 !rg=5.36469 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,16,1:6)=(/ 12.671000,12.738000,12.778000,12.913000,13.126000,13.338000 /)
XPIZA_LKT(79,16,1:6)=(/ 0.651615,0.805992,0.990353,0.991049,0.968264,0.542784 /)
XCGA_LKT(79,16,1:6)=(/ 0.907110,0.870477,0.838263,0.832713,0.829073,0.934377 /)
XEXT_COEFF_550_LKT(79,16)=12.752000 !rg=5.36469 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,17,1:6)=(/ 10.538000,10.575000,10.618000,10.726000,10.869000,11.038000 /)
XPIZA_LKT(79,17,1:6)=(/ 0.638840,0.789478,0.989087,0.990166,0.966801,0.544787 /)
XCGA_LKT(79,17,1:6)=(/ 0.910833,0.873107,0.839103,0.834140,0.832583,0.935003 /)
XEXT_COEFF_550_LKT(79,17)=10.601000 !rg=5.36469 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,18,1:6)=(/ 8.769600,8.788900,8.846300,8.907700,9.027900,9.140300 /)
XPIZA_LKT(79,18,1:6)=(/ 0.627098,0.772310,0.987645,0.989224,0.964438,0.546547 /)
XCGA_LKT(79,18,1:6)=(/ 0.914727,0.875917,0.840467,0.836917,0.835670,0.935510 /)
XEXT_COEFF_550_LKT(79,18)=8.824200 !rg=5.36469 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,19,1:6)=(/ 7.279900,7.307800,7.327100,7.387100,7.459500,7.559400 /)
XPIZA_LKT(79,19,1:6)=(/ 0.615791,0.754875,0.985860,0.987964,0.961446,0.548106 /)
XCGA_LKT(79,19,1:6)=(/ 0.000000,0.879367,0.841483,0.838110,0.837420,0.935927 /)
XEXT_COEFF_550_LKT(79,19)=7.309300 !rg=5.36469 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,20,1:6)=(/ 6.049500,6.064800,6.081800,6.128200,6.184600,6.258000 /)
XPIZA_LKT(79,20,1:6)=(/ 0.605744,0.736734,0.983775,0.986304,0.958091,0.549452 /)
XCGA_LKT(79,20,1:6)=(/ 0.000000,0.882677,0.842343,0.839127,0.839710,0.936260 /)
XEXT_COEFF_550_LKT(79,20)=6.070100 !rg=5.36469 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,1,1:6)=(/ 106.800000,110.130000,108.700000,113.640000,123.170000,122.710000 /)
XPIZA_LKT(80,1,1:6)=(/ 0.790139,0.877568,0.995145,0.995780,0.989875,0.498742 /)
XCGA_LKT(80,1,1:6)=(/ 0.876840,0.852420,0.817833,0.785270,0.778897,0.916957 /)
XEXT_COEFF_550_LKT(80,1)=110.800000 !rg=5.81187 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,2,1:6)=(/ 102.180000,103.450000,105.620000,107.160000,112.390000,117.520000 /)
XPIZA_LKT(80,2,1:6)=(/ 0.788312,0.884551,0.994631,0.996173,0.986810,0.501095 /)
XCGA_LKT(80,2,1:6)=(/ 0.877490,0.852050,0.812810,0.803767,0.776463,0.917690 /)
XEXT_COEFF_550_LKT(80,2)=103.790000 !rg=5.81187 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,3,1:6)=(/ 94.536000,96.165000,97.496000,99.908000,103.310000,108.340000 /)
XPIZA_LKT(80,3,1:6)=(/ 0.784770,0.888138,0.993895,0.995072,0.988438,0.502903 /)
XCGA_LKT(80,3,1:6)=(/ 0.877580,0.849473,0.816420,0.799927,0.782333,0.918113 /)
XEXT_COEFF_550_LKT(80,3)=95.867000 !rg=5.81187 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,4,1:6)=(/ 85.585000,86.819000,87.660000,90.370000,93.263000,97.239000 /)
XPIZA_LKT(80,4,1:6)=(/ 0.780236,0.890103,0.994342,0.994812,0.986530,0.505238 /)
XCGA_LKT(80,4,1:6)=(/ 0.878903,0.849427,0.818847,0.799853,0.783657,0.918667 /)
XEXT_COEFF_550_LKT(80,4)=86.592000 !rg=5.81187 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,5,1:6)=(/ 75.766000,76.798000,77.245000,79.530000,82.013000,85.504000 /)
XPIZA_LKT(80,5,1:6)=(/ 0.774662,0.890265,0.994210,0.994379,0.986029,0.508396 /)
XCGA_LKT(80,5,1:6)=(/ 0.879930,0.851377,0.819103,0.804177,0.785257,0.919673 /)
XEXT_COEFF_550_LKT(80,5)=76.700000 !rg=5.81187 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,6,1:6)=(/ 66.023000,66.625000,67.360000,69.126000,71.749000,74.087000 /)
XPIZA_LKT(80,6,1:6)=(/ 0.768457,0.889725,0.994308,0.993955,0.984051,0.512027 /)
XCGA_LKT(80,6,1:6)=(/ 0.881523,0.852387,0.820370,0.805617,0.788650,0.921093 /)
XEXT_COEFF_550_LKT(80,6)=67.144000 !rg=5.81187 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,7,1:6)=(/ 56.993000,57.505000,58.051000,59.353000,61.090000,63.309000 /)
XPIZA_LKT(80,7,1:6)=(/ 0.760695,0.887125,0.993935,0.994043,0.982276,0.515895 /)
XCGA_LKT(80,7,1:6)=(/ 0.883403,0.852430,0.823150,0.809867,0.795680,0.922863 /)
XEXT_COEFF_550_LKT(80,7)=57.628000 !rg=5.81187 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,8,1:6)=(/ 48.621000,48.823000,49.677000,50.440000,51.998000,53.671000 /)
XPIZA_LKT(80,8,1:6)=(/ 0.751261,0.883408,0.994111,0.993489,0.980086,0.519730 /)
XCGA_LKT(80,8,1:6)=(/ 0.885127,0.854667,0.824113,0.813420,0.799787,0.924710 /)
XEXT_COEFF_550_LKT(80,8)=49.230000 !rg=5.81187 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,9,1:6)=(/ 41.067000,41.538000,41.823000,42.427000,43.751000,45.128000 /)
XPIZA_LKT(80,9,1:6)=(/ 0.739778,0.878167,0.993837,0.993453,0.977885,0.523494 /)
XCGA_LKT(80,9,1:6)=(/ 0.886620,0.856307,0.826810,0.817837,0.803460,0.926560 /)
XEXT_COEFF_550_LKT(80,9)=41.480000 !rg=5.81187 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,10,1:6)=(/ 34.732000,35.010000,35.281000,35.863000,36.915000,37.808000 /)
XPIZA_LKT(80,10,1:6)=(/ 0.728293,0.871876,0.993799,0.993449,0.976118,0.527058 /)
XCGA_LKT(80,10,1:6)=(/ 0.889527,0.858067,0.829487,0.820697,0.807233,0.928250 /)
XEXT_COEFF_550_LKT(80,10)=35.033000 !rg=5.81187 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,11,1:6)=(/ 29.117000,29.286000,29.573000,29.906000,30.760000,31.487000 /)
XPIZA_LKT(80,11,1:6)=(/ 0.715133,0.863296,0.993447,0.993257,0.974935,0.530473 /)
XCGA_LKT(80,11,1:6)=(/ 0.892023,0.859730,0.830630,0.823863,0.813937,0.929807 /)
XEXT_COEFF_550_LKT(80,11)=29.351000 !rg=5.81187 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,12,1:6)=(/ 24.351000,24.473000,24.650000,25.026000,25.654000,26.139000 /)
XPIZA_LKT(80,12,1:6)=(/ 0.701355,0.853416,0.992946,0.992740,0.973070,0.533664 /)
XCGA_LKT(80,12,1:6)=(/ 0.895030,0.861840,0.832710,0.825240,0.815797,0.931173 /)
XEXT_COEFF_550_LKT(80,12)=24.502000 !rg=5.81187 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,13,1:6)=(/ 20.281000,20.396000,20.527000,20.783000,21.208000,21.652000 /)
XPIZA_LKT(80,13,1:6)=(/ 0.686941,0.841957,0.992476,0.992502,0.972246,0.536611 /)
XCGA_LKT(80,13,1:6)=(/ 0.898103,0.864077,0.834113,0.828000,0.821753,0.932340 /)
XEXT_COEFF_550_LKT(80,13)=20.417000 !rg=5.81187 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,14,1:6)=(/ 16.911000,16.980000,17.123000,17.348000,17.625000,17.938000 /)
XPIZA_LKT(80,14,1:6)=(/ 0.673072,0.828840,0.991823,0.992069,0.970805,0.539265 /)
XCGA_LKT(80,14,1:6)=(/ 0.901553,0.866373,0.836627,0.830663,0.824453,0.933310 /)
XEXT_COEFF_550_LKT(80,14)=17.070000 !rg=5.81187 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,15,1:6)=(/ 14.072000,14.152000,14.211000,14.325000,14.615000,14.852000 /)
XPIZA_LKT(80,15,1:6)=(/ 0.659044,0.814559,0.990898,0.991504,0.969446,0.541656 /)
XCGA_LKT(80,15,1:6)=(/ 0.905000,0.869090,0.837173,0.832583,0.827113,0.934113 /)
XEXT_COEFF_550_LKT(80,15)=14.154000 !rg=5.81187 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,16,1:6)=(/ 11.691000,11.731000,11.808000,11.888000,12.097000,12.271000 /)
XPIZA_LKT(80,16,1:6)=(/ 0.645737,0.798538,0.989852,0.990743,0.967918,0.543822 /)
XCGA_LKT(80,16,1:6)=(/ 0.908773,0.871540,0.838770,0.834383,0.832020,0.934790 /)
XEXT_COEFF_550_LKT(80,16)=11.745000 !rg=5.81187 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,17,1:6)=(/ 9.726500,9.750600,9.818300,9.888000,10.002000,10.159000 /)
XPIZA_LKT(80,17,1:6)=(/ 0.633300,0.781767,0.988486,0.989747,0.965666,0.545705 /)
XCGA_LKT(80,17,1:6)=(/ 0.912673,0.874290,0.839857,0.835643,0.833197,0.935330 /)
XEXT_COEFF_550_LKT(80,17)=9.778600 !rg=5.81187 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,18,1:6)=(/ 8.087100,8.115300,8.148400,8.198700,8.298800,8.414500 /)
XPIZA_LKT(80,18,1:6)=(/ 0.621675,0.764762,0.986915,0.988675,0.963219,0.547352 /)
XCGA_LKT(80,18,1:6)=(/ 0.916417,0.877480,0.840803,0.837513,0.835440,0.935770 /)
XEXT_COEFF_550_LKT(80,18)=8.120800 !rg=5.81187 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,19,1:6)=(/ 6.716200,6.734500,6.761300,6.798400,6.888400,6.961000 /)
XPIZA_LKT(80,19,1:6)=(/ 0.610932,0.746575,0.985044,0.987339,0.960174,0.548807 /)
XCGA_LKT(80,19,1:6)=(/ 0.000000,0.880810,0.842047,0.838953,0.839127,0.936130 /)
XEXT_COEFF_550_LKT(80,19)=6.739600 !rg=5.81187 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,20,1:6)=(/ 5.582400,5.593100,5.615000,5.656300,5.697300,5.763900 /)
XPIZA_LKT(80,20,1:6)=(/ 0.601415,0.728525,0.982714,0.985643,0.956246,0.550061 /)
XCGA_LKT(80,20,1:6)=(/ 0.000000,0.884203,0.842857,0.840177,0.840350,0.936417 /)
XEXT_COEFF_550_LKT(80,20)=5.604000 !rg=5.81187 sigma=2.95 wvl=0.55
 
END SUBROUTINE SALT_OPT_LKT_SET8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET9()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE
 
XEXT_COEFF_WVL_LKT(81,1,1:6)=(/ 98.064000,100.610000,102.210000,105.730000,117.820000,111.720000 /)
XPIZA_LKT(81,1,1:6)=(/ 0.787577,0.891066,0.994135,0.995947,0.984825,0.499238 /)
XCGA_LKT(81,1,1:6)=(/ 0.877200,0.851233,0.812057,0.798647,0.789913,0.917877 /)
XEXT_COEFF_550_LKT(81,1)=101.220000 !rg=6.29632 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,2,1:6)=(/ 94.163000,95.226000,96.161000,98.847000,104.460000,107.330000 /)
XPIZA_LKT(81,2,1:6)=(/ 0.786335,0.888594,0.994155,0.995421,0.987942,0.503010 /)
XCGA_LKT(81,2,1:6)=(/ 0.878107,0.851250,0.816270,0.797973,0.780473,0.919490 /)
XEXT_COEFF_550_LKT(81,2)=96.531000 !rg=6.29632 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,3,1:6)=(/ 87.199000,87.945000,89.112000,91.223000,94.622000,99.128000 /)
XPIZA_LKT(81,3,1:6)=(/ 0.783017,0.890566,0.994074,0.995442,0.987834,0.505509 /)
XCGA_LKT(81,3,1:6)=(/ 0.878473,0.850407,0.817803,0.801660,0.781603,0.920510 /)
XEXT_COEFF_550_LKT(81,3)=88.089000 !rg=6.29632 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,4,1:6)=(/ 78.876000,79.501000,80.289000,82.336000,85.536000,89.046000 /)
XPIZA_LKT(81,4,1:6)=(/ 0.778169,0.890947,0.994108,0.995005,0.986123,0.507987 /)
XCGA_LKT(81,4,1:6)=(/ 0.879923,0.851120,0.819853,0.802720,0.784047,0.921293 /)
XEXT_COEFF_550_LKT(81,4)=79.945000 !rg=6.29632 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,5,1:6)=(/ 69.684000,70.562000,71.342000,73.402000,75.258000,78.334000 /)
XPIZA_LKT(81,5,1:6)=(/ 0.772029,0.890444,0.994165,0.993926,0.984796,0.511015 /)
XCGA_LKT(81,5,1:6)=(/ 0.880670,0.851087,0.821507,0.805093,0.791580,0.922253 /)
XEXT_COEFF_550_LKT(81,5)=70.557000 !rg=6.29632 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,6,1:6)=(/ 60.842000,61.461000,62.708000,63.762000,65.870000,67.897000 /)
XPIZA_LKT(81,6,1:6)=(/ 0.764754,0.888402,0.994008,0.993865,0.983050,0.514446 /)
XCGA_LKT(81,6,1:6)=(/ 0.882327,0.853080,0.818780,0.810037,0.789383,0.923480 /)
XEXT_COEFF_550_LKT(81,6)=61.793000 !rg=6.29632 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,7,1:6)=(/ 52.342000,52.856000,53.098000,54.380000,56.518000,58.042000 /)
XPIZA_LKT(81,7,1:6)=(/ 0.755471,0.885592,0.994221,0.993780,0.979717,0.518136 /)
XCGA_LKT(81,7,1:6)=(/ 0.884183,0.853937,0.825530,0.810723,0.794297,0.924970 /)
XEXT_COEFF_550_LKT(81,7)=53.331000 !rg=6.29632 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,8,1:6)=(/ 44.757000,45.113000,45.731000,46.456000,47.833000,49.226000 /)
XPIZA_LKT(81,8,1:6)=(/ 0.745853,0.881400,0.993990,0.993633,0.978886,0.521815 /)
XCGA_LKT(81,8,1:6)=(/ 0.885913,0.855590,0.826833,0.816847,0.800550,0.926527 /)
XEXT_COEFF_550_LKT(81,8)=45.138000 !rg=6.29632 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,9,1:6)=(/ 37.915000,38.142000,38.637000,39.269000,40.223000,41.410000 /)
XPIZA_LKT(81,9,1:6)=(/ 0.734477,0.875203,0.993896,0.993274,0.977168,0.525440 /)
XCGA_LKT(81,9,1:6)=(/ 0.888217,0.856563,0.828263,0.818113,0.809337,0.928087 /)
XEXT_COEFF_550_LKT(81,9)=38.265000 !rg=6.29632 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,10,1:6)=(/ 31.970000,32.244000,32.602000,33.049000,33.787000,34.709000 /)
XPIZA_LKT(81,10,1:6)=(/ 0.721893,0.868211,0.993327,0.993088,0.975111,0.528865 /)
XCGA_LKT(81,10,1:6)=(/ 0.890550,0.858647,0.829930,0.821080,0.810943,0.929523 /)
XEXT_COEFF_550_LKT(81,10)=32.234000 !rg=6.29632 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,11,1:6)=(/ 26.845000,26.990000,27.158000,27.575000,28.161000,28.919000 /)
XPIZA_LKT(81,11,1:6)=(/ 0.708705,0.858988,0.993346,0.992971,0.974314,0.532140 /)
XCGA_LKT(81,11,1:6)=(/ 0.893557,0.860653,0.832480,0.825370,0.816463,0.930853 /)
XEXT_COEFF_550_LKT(81,11)=27.134000 !rg=6.29632 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,12,1:6)=(/ 22.425000,22.552000,22.752000,22.984000,23.433000,24.018000 /)
XPIZA_LKT(81,12,1:6)=(/ 0.694358,0.848274,0.992840,0.992745,0.973132,0.535191 /)
XCGA_LKT(81,12,1:6)=(/ 0.896390,0.862697,0.833573,0.826537,0.820383,0.932023 /)
XEXT_COEFF_550_LKT(81,12)=22.595000 !rg=6.29632 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,13,1:6)=(/ 18.718000,18.775000,18.970000,19.245000,19.521000,19.903000 /)
XPIZA_LKT(81,13,1:6)=(/ 0.680412,0.835966,0.992213,0.992234,0.971493,0.537996 /)
XCGA_LKT(81,13,1:6)=(/ 0.899890,0.865003,0.835600,0.829410,0.822277,0.933030 /)
XEXT_COEFF_550_LKT(81,13)=18.894000 !rg=6.29632 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,14,1:6)=(/ 15.590000,15.688000,15.770000,15.922000,16.176000,16.496000 /)
XPIZA_LKT(81,14,1:6)=(/ 0.666172,0.822791,0.991250,0.991704,0.970184,0.540512 /)
XCGA_LKT(81,14,1:6)=(/ 0.903170,0.867920,0.835993,0.831630,0.826310,0.933867 /)
XEXT_COEFF_550_LKT(81,14)=15.671000 !rg=6.29632 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,15,1:6)=(/ 12.977000,13.028000,13.129000,13.195000,13.439000,13.662000 /)
XPIZA_LKT(81,15,1:6)=(/ 0.652636,0.807401,0.990430,0.991134,0.968430,0.542771 /)
XCGA_LKT(81,15,1:6)=(/ 0.906723,0.869893,0.838400,0.833070,0.829220,0.934560 /)
XEXT_COEFF_550_LKT(81,15)=13.066000 !rg=6.29632 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,16,1:6)=(/ 10.785000,10.820000,10.873000,10.974000,11.098000,11.292000 /)
XPIZA_LKT(81,16,1:6)=(/ 0.639637,0.791118,0.989203,0.990359,0.967097,0.544808 /)
XCGA_LKT(81,16,1:6)=(/ 0.910660,0.872860,0.839257,0.835720,0.833053,0.935147 /)
XEXT_COEFF_550_LKT(81,16)=10.853000 !rg=6.29632 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,17,1:6)=(/ 8.967400,9.004400,9.040100,9.106400,9.206400,9.351300 /)
XPIZA_LKT(81,17,1:6)=(/ 0.627514,0.774325,0.987818,0.989342,0.964755,0.546572 /)
XCGA_LKT(81,17,1:6)=(/ 0.914397,0.875887,0.840110,0.836647,0.835783,0.935613 /)
XEXT_COEFF_550_LKT(81,17)=9.014100 !rg=6.29632 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,18,1:6)=(/ 7.460500,7.481200,7.524400,7.552900,7.646400,7.747500 /)
XPIZA_LKT(81,18,1:6)=(/ 0.616470,0.756574,0.986149,0.988073,0.961912,0.548110 /)
XCGA_LKT(81,18,1:6)=(/ 0.918210,0.878863,0.841657,0.837943,0.837337,0.935993 /)
XEXT_COEFF_550_LKT(81,18)=7.494600 !rg=6.29632 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,19,1:6)=(/ 6.197200,6.211700,6.232200,6.274500,6.328400,6.410800 /)
XPIZA_LKT(81,19,1:6)=(/ 0.606265,0.738510,0.984058,0.986617,0.958550,0.549464 /)
XCGA_LKT(81,19,1:6)=(/ 0.000000,0.882270,0.842383,0.839820,0.839673,0.936307 /)
XEXT_COEFF_550_LKT(81,19)=6.220200 !rg=6.29632 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,20,1:6)=(/ 5.149500,5.163500,5.179700,5.205700,5.246800,5.309500 /)
XPIZA_LKT(81,20,1:6)=(/ 0.597109,0.720638,0.981619,0.984824,0.954505,0.550630 /)
XCGA_LKT(81,20,1:6)=(/ 0.000000,0.886017,0.843233,0.840697,0.842013,0.936550 /)
XEXT_COEFF_550_LKT(81,20)=5.168600 !rg=6.29632 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,1,1:6)=(/ 90.615000,91.915000,92.106000,92.566000,96.244000,102.590000 /)
XPIZA_LKT(82,1,1:6)=(/ 0.787269,0.893941,0.993798,0.995726,0.988708,0.503572 /)
XCGA_LKT(82,1,1:6)=(/ 0.877653,0.848833,0.817230,0.803240,0.778107,0.920173 /)
XEXT_COEFF_550_LKT(82,1)=92.517000 !rg=6.82116 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,2,1:6)=(/ 86.803000,88.160000,89.447000,91.896000,95.228000,98.243000 /)
XPIZA_LKT(82,2,1:6)=(/ 0.785213,0.892494,0.993744,0.994640,0.987289,0.505733 /)
XCGA_LKT(82,2,1:6)=(/ 0.878330,0.851650,0.816467,0.799487,0.775630,0.921387 /)
XEXT_COEFF_550_LKT(82,2)=88.577000 !rg=6.82116 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,3,1:6)=(/ 80.262000,81.782000,82.510000,83.427000,87.585000,90.742000 /)
XPIZA_LKT(82,3,1:6)=(/ 0.780416,0.891632,0.993913,0.995045,0.986399,0.508150 /)
XCGA_LKT(82,3,1:6)=(/ 0.879107,0.849960,0.819017,0.807933,0.786537,0.922510 /)
XEXT_COEFF_550_LKT(82,3)=81.246000 !rg=6.82116 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,4,1:6)=(/ 72.507000,73.751000,74.551000,75.373000,78.945000,81.566000 /)
XPIZA_LKT(82,4,1:6)=(/ 0.775057,0.891206,0.993997,0.994597,0.984974,0.510697 /)
XCGA_LKT(82,4,1:6)=(/ 0.880043,0.851193,0.819593,0.808890,0.789583,0.923407 /)
XEXT_COEFF_550_LKT(82,4)=73.572000 !rg=6.82116 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,5,1:6)=(/ 64.213000,64.708000,65.551000,66.649000,69.409000,71.787000 /)
XPIZA_LKT(82,5,1:6)=(/ 0.768239,0.890331,0.994110,0.994240,0.983609,0.513625 /)
XCGA_LKT(82,5,1:6)=(/ 0.881600,0.852340,0.822663,0.809270,0.791797,0.924347 /)
XEXT_COEFF_550_LKT(82,5)=64.967000 !rg=6.82116 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,6,1:6)=(/ 56.023000,56.486000,57.099000,58.776000,60.529000,62.245000 /)
XPIZA_LKT(82,6,1:6)=(/ 0.760089,0.887695,0.994253,0.993808,0.980896,0.516890 /)
XCGA_LKT(82,6,1:6)=(/ 0.882993,0.853753,0.824560,0.810253,0.794487,0.925433 /)
XEXT_COEFF_550_LKT(82,6)=57.005000 !rg=6.82116 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,7,1:6)=(/ 48.315000,48.645000,49.251000,49.947000,51.743000,53.231000 /)
XPIZA_LKT(82,7,1:6)=(/ 0.751079,0.883802,0.993996,0.993627,0.979631,0.520401 /)
XCGA_LKT(82,7,1:6)=(/ 0.885143,0.854967,0.824953,0.815690,0.801627,0.926717 /)
XEXT_COEFF_550_LKT(82,7)=48.776000 !rg=6.82116 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,8,1:6)=(/ 41.212000,41.506000,42.068000,42.764000,43.864000,45.165000 /)
XPIZA_LKT(82,8,1:6)=(/ 0.740291,0.878853,0.993943,0.993405,0.977676,0.523910 /)
XCGA_LKT(82,8,1:6)=(/ 0.887087,0.856210,0.826920,0.817993,0.805120,0.928043 /)
XEXT_COEFF_550_LKT(82,8)=41.621000 !rg=6.82116 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,9,1:6)=(/ 34.944000,35.214000,35.525000,35.847000,37.145000,38.010000 /)
XPIZA_LKT(82,9,1:6)=(/ 0.728355,0.872490,0.993759,0.993345,0.975695,0.527374 /)
XCGA_LKT(82,9,1:6)=(/ 0.889357,0.858160,0.829327,0.821303,0.809240,0.929377 /)
XEXT_COEFF_550_LKT(82,9)=35.245000 !rg=6.82116 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,10,1:6)=(/ 29.534000,29.632000,29.961000,30.444000,31.214000,31.873000 /)
XPIZA_LKT(82,10,1:6)=(/ 0.715908,0.864102,0.993529,0.993208,0.974235,0.530646 /)
XCGA_LKT(82,10,1:6)=(/ 0.892143,0.859793,0.831530,0.824123,0.813997,0.930607 /)
XEXT_COEFF_550_LKT(82,10)=29.860000 !rg=6.82116 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,11,1:6)=(/ 24.745000,24.925000,25.033000,25.299000,25.934000,26.568000 /)
XPIZA_LKT(82,11,1:6)=(/ 0.701919,0.854667,0.993070,0.992952,0.972969,0.533769 /)
XCGA_LKT(82,11,1:6)=(/ 0.894797,0.862093,0.832720,0.826737,0.818683,0.931747 /)
XEXT_COEFF_550_LKT(82,11)=24.932000 !rg=6.82116 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,12,1:6)=(/ 20.709000,20.779000,20.971000,21.248000,21.659000,22.075000 /)
XPIZA_LKT(82,12,1:6)=(/ 0.688170,0.843058,0.992574,0.992600,0.971840,0.536671 /)
XCGA_LKT(82,12,1:6)=(/ 0.898040,0.863987,0.835090,0.828927,0.821447,0.932757 /)
XEXT_COEFF_550_LKT(82,12)=20.890000 !rg=6.82116 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,13,1:6)=(/ 17.250000,17.357000,17.438000,17.630000,17.888000,18.299000 /)
XPIZA_LKT(82,13,1:6)=(/ 0.673522,0.830452,0.991892,0.992048,0.971107,0.539329 /)
XCGA_LKT(82,13,1:6)=(/ 0.901340,0.866493,0.835250,0.830733,0.825653,0.933623 /)
XEXT_COEFF_550_LKT(82,13)=17.371000 !rg=6.82116 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,14,1:6)=(/ 14.385000,14.439000,14.559000,14.687000,14.887000,15.172000 /)
XPIZA_LKT(82,14,1:6)=(/ 0.659894,0.816038,0.991002,0.991591,0.969916,0.541704 /)
XCGA_LKT(82,14,1:6)=(/ 0.904857,0.868753,0.837433,0.832680,0.829157,0.934347 /)
XEXT_COEFF_550_LKT(82,14)=14.483000 !rg=6.82116 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,15,1:6)=(/ 11.977000,12.030000,12.079000,12.217000,12.392000,12.570000 /)
XPIZA_LKT(82,15,1:6)=(/ 0.646578,0.800677,0.989950,0.990900,0.967861,0.543829 /)
XCGA_LKT(82,15,1:6)=(/ 0.908533,0.871480,0.838837,0.834573,0.831143,0.934947 /)
XEXT_COEFF_550_LKT(82,15)=12.041000 !rg=6.82116 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,16,1:6)=(/ 9.946800,9.985500,10.021000,10.099000,10.238000,10.393000 /)
XPIZA_LKT(82,16,1:6)=(/ 0.633775,0.783706,0.988664,0.990033,0.965796,0.545739 /)
XCGA_LKT(82,16,1:6)=(/ 0.912333,0.874407,0.839923,0.836210,0.834193,0.935453 /)
XEXT_COEFF_550_LKT(82,16)=9.994900 !rg=6.82116 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,17,1:6)=(/ 8.276400,8.297400,8.348900,8.389200,8.501900,8.608900 /)
XPIZA_LKT(82,17,1:6)=(/ 0.622281,0.766289,0.987177,0.988844,0.963584,0.547387 /)
XCGA_LKT(82,17,1:6)=(/ 0.916197,0.877190,0.841113,0.837367,0.837377,0.935860 /)
XEXT_COEFF_550_LKT(82,17)=8.321300 !rg=6.82116 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,18,1:6)=(/ 6.883800,6.905900,6.925900,6.985300,7.052100,7.134200 /)
XPIZA_LKT(82,18,1:6)=(/ 0.611537,0.748666,0.985277,0.987529,0.960513,0.548819 /)
XCGA_LKT(82,18,1:6)=(/ 0.000000,0.880523,0.841970,0.839010,0.838413,0.936187 /)
XEXT_COEFF_550_LKT(82,18)=6.908800 !rg=6.82116 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,19,1:6)=(/ 5.717500,5.732100,5.744200,5.780700,5.825900,5.904700 /)
XPIZA_LKT(82,19,1:6)=(/ 0.601787,0.730488,0.983071,0.985888,0.956682,0.550079 /)
XCGA_LKT(82,19,1:6)=(/ 0.000000,0.884033,0.842877,0.840240,0.840740,0.936457 /)
XEXT_COEFF_550_LKT(82,19)=5.738500 !rg=6.82116 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,20,1:6)=(/ 4.751900,4.761200,4.777200,4.800600,4.846700,4.891500 /)
XPIZA_LKT(82,20,1:6)=(/ 0.593163,0.712591,0.980441,0.983890,0.952459,0.551161 /)
XCGA_LKT(82,20,1:6)=(/ 0.000000,0.000000,0.843840,0.841190,0.843367,0.936663 /)
XEXT_COEFF_550_LKT(82,20)=4.768300 !rg=6.82116 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,1,1:6)=(/ 83.502000,83.199000,85.489000,88.307000,84.549000,94.407000 /)
XPIZA_LKT(83,1,1:6)=(/ 0.783549,0.896571,0.993250,0.995418,0.986828,0.508578 /)
XCGA_LKT(83,1,1:6)=(/ 0.878677,0.850373,0.820587,0.808783,0.760790,0.922860 /)
XEXT_COEFF_550_LKT(83,1)=85.589000 !rg=7.38975 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,2,1:6)=(/ 79.952000,80.505000,81.732000,84.301000,86.976000,90.039000 /)
XPIZA_LKT(83,2,1:6)=(/ 0.781038,0.895008,0.993519,0.995042,0.982283,0.508807 /)
XCGA_LKT(83,2,1:6)=(/ 0.878920,0.850823,0.821663,0.802413,0.781210,0.923320 /)
XEXT_COEFF_550_LKT(83,2)=81.753000 !rg=7.38975 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,3,1:6)=(/ 74.110000,74.687000,75.410000,77.171000,79.793000,83.116000 /)
XPIZA_LKT(83,3,1:6)=(/ 0.777472,0.893878,0.993277,0.994671,0.986349,0.510843 /)
XCGA_LKT(83,3,1:6)=(/ 0.880257,0.851137,0.821007,0.809627,0.790770,0.924267 /)
XEXT_COEFF_550_LKT(83,3)=75.153000 !rg=7.38975 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,4,1:6)=(/ 67.022000,67.461000,68.015000,69.711000,71.885000,74.740000 /)
XPIZA_LKT(83,4,1:6)=(/ 0.771881,0.892071,0.993694,0.993610,0.984493,0.513360 /)
XCGA_LKT(83,4,1:6)=(/ 0.881513,0.852020,0.821620,0.809820,0.792660,0.925167 /)
XEXT_COEFF_550_LKT(83,4)=67.945000 !rg=7.38975 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,5,1:6)=(/ 59.166000,59.940000,60.637000,61.228000,64.010000,65.808000 /)
XPIZA_LKT(83,5,1:6)=(/ 0.764070,0.889309,0.994122,0.993521,0.981822,0.516202 /)
XCGA_LKT(83,5,1:6)=(/ 0.882407,0.852573,0.823797,0.813963,0.795450,0.926067 /)
XEXT_COEFF_550_LKT(83,5)=59.929000 !rg=7.38975 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,6,1:6)=(/ 51.788000,52.089000,52.887000,53.764000,55.253000,57.083000 /)
XPIZA_LKT(83,6,1:6)=(/ 0.756106,0.885554,0.994015,0.993566,0.980300,0.519309 /)
XCGA_LKT(83,6,1:6)=(/ 0.884440,0.854600,0.824217,0.813057,0.797020,0.927053 /)
XEXT_COEFF_550_LKT(83,6)=52.313000 !rg=7.38975 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,7,1:6)=(/ 44.531000,44.788000,45.089000,45.945000,47.268000,48.835000 /)
XPIZA_LKT(83,7,1:6)=(/ 0.745823,0.881723,0.994081,0.993461,0.978731,0.522643 /)
XCGA_LKT(83,7,1:6)=(/ 0.886460,0.855667,0.827600,0.817677,0.804310,0.928173 /)
XEXT_COEFF_550_LKT(83,7)=45.042000 !rg=7.38975 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,8,1:6)=(/ 38.030000,38.272000,38.739000,39.475000,40.373000,41.452000 /)
XPIZA_LKT(83,8,1:6)=(/ 0.734646,0.876005,0.993740,0.992873,0.976362,0.525978 /)
XCGA_LKT(83,8,1:6)=(/ 0.888487,0.857483,0.829237,0.819070,0.806603,0.929320 /)
XEXT_COEFF_550_LKT(83,8)=38.531000 !rg=7.38975 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,9,1:6)=(/ 32.188000,32.435000,32.813000,33.063000,34.139000,34.900000 /)
XPIZA_LKT(83,9,1:6)=(/ 0.721984,0.868539,0.993680,0.993095,0.974629,0.529274 /)
XCGA_LKT(83,9,1:6)=(/ 0.890563,0.858527,0.830873,0.822920,0.811730,0.930467 /)
XEXT_COEFF_550_LKT(83,9)=32.567000 !rg=7.38975 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,10,1:6)=(/ 27.184000,27.359000,27.591000,27.881000,28.671000,29.277000 /)
XPIZA_LKT(83,10,1:6)=(/ 0.708968,0.860114,0.993318,0.993093,0.974018,0.532382 /)
XCGA_LKT(83,10,1:6)=(/ 0.893157,0.861110,0.831900,0.825153,0.815777,0.931527 /)
XEXT_COEFF_550_LKT(83,10)=27.330000 !rg=7.38975 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,11,1:6)=(/ 22.850000,22.918000,23.171000,23.406000,23.891000,24.414000 /)
XPIZA_LKT(83,11,1:6)=(/ 0.695562,0.849403,0.992893,0.992725,0.972596,0.535347 /)
XCGA_LKT(83,11,1:6)=(/ 0.896487,0.862403,0.834390,0.827313,0.821490,0.932510 /)
XEXT_COEFF_550_LKT(83,11)=23.058000 !rg=7.38975 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,12,1:6)=(/ 19.079000,19.202000,19.330000,19.509000,19.883000,20.293000 /)
XPIZA_LKT(83,12,1:6)=(/ 0.680964,0.837664,0.992189,0.992377,0.971392,0.538093 /)
XCGA_LKT(83,12,1:6)=(/ 0.899500,0.865557,0.834930,0.830113,0.823263,0.933383 /)
XEXT_COEFF_550_LKT(83,12)=19.201000 !rg=7.38975 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,13,1:6)=(/ 15.915000,15.969000,16.091000,16.266000,16.477000,16.828000 /)
XPIZA_LKT(83,13,1:6)=(/ 0.667093,0.823976,0.991520,0.991837,0.970382,0.540603 /)
XCGA_LKT(83,13,1:6)=(/ 0.902983,0.867377,0.836750,0.831487,0.827940,0.934133 /)
XEXT_COEFF_550_LKT(83,13)=16.048000 !rg=7.38975 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,14,1:6)=(/ 13.264000,13.326000,13.384000,13.514000,13.724000,13.957000 /)
XPIZA_LKT(83,14,1:6)=(/ 0.653350,0.809320,0.990579,0.991271,0.968691,0.542837 /)
XCGA_LKT(83,14,1:6)=(/ 0.906507,0.870017,0.838067,0.833910,0.829327,0.934760 /)
XEXT_COEFF_550_LKT(83,14)=13.343000 !rg=7.38975 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,15,1:6)=(/ 11.045000,11.087000,11.144000,11.226000,11.369000,11.567000 /)
XPIZA_LKT(83,15,1:6)=(/ 0.640419,0.793205,0.989431,0.990501,0.967227,0.544830 /)
XCGA_LKT(83,15,1:6)=(/ 0.910227,0.872617,0.839157,0.835357,0.832640,0.935280 /)
XEXT_COEFF_550_LKT(83,15)=11.095000 !rg=7.38975 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,16,1:6)=(/ 9.180600,9.203100,9.266300,9.316300,9.439500,9.566500 /)
XPIZA_LKT(83,16,1:6)=(/ 0.628222,0.775862,0.988050,0.989509,0.964872,0.546616 /)
XCGA_LKT(83,16,1:6)=(/ 0.914233,0.875490,0.840757,0.836680,0.835497,0.935720 /)
XEXT_COEFF_550_LKT(83,16)=9.229200 !rg=7.38975 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,17,1:6)=(/ 7.633400,7.655400,7.688900,7.750600,7.801500,7.926500 /)
XPIZA_LKT(83,17,1:6)=(/ 0.616908,0.758411,0.986351,0.988337,0.962371,0.548152 /)
XCGA_LKT(83,17,1:6)=(/ 0.918003,0.878707,0.841160,0.838487,0.837927,0.936070 /)
XEXT_COEFF_550_LKT(83,17)=7.663200 !rg=7.38975 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,18,1:6)=(/ 6.350300,6.367900,6.391800,6.429100,6.478900,6.570200 /)
XPIZA_LKT(83,18,1:6)=(/ 0.606750,0.740466,0.984329,0.986843,0.959018,0.549483 /)
XCGA_LKT(83,18,1:6)=(/ 0.000000,0.882040,0.842327,0.839577,0.839383,0.936353 /)
XEXT_COEFF_550_LKT(83,18)=6.371600 !rg=7.38975 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,19,1:6)=(/ 5.276200,5.288200,5.307500,5.339100,5.384700,5.439200 /)
XPIZA_LKT(83,19,1:6)=(/ 0.597598,0.722361,0.981997,0.985075,0.954974,0.550651 /)
XCGA_LKT(83,19,1:6)=(/ 0.000000,0.885637,0.843543,0.840760,0.841640,0.936583 /)
XEXT_COEFF_550_LKT(83,19)=5.295300 !rg=7.38975 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,20,1:6)=(/ 4.384500,4.393700,4.405100,4.431200,4.456700,4.506800 /)
XPIZA_LKT(83,20,1:6)=(/ 0.589329,0.704698,0.979126,0.982983,0.950222,0.551654 /)
XCGA_LKT(83,20,1:6)=(/ 0.000000,0.000000,0.844127,0.841980,0.843707,0.936760 /)
XEXT_COEFF_550_LKT(83,20)=4.395900 !rg=7.38975 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,1,1:6)=(/ 76.978000,77.979000,77.664000,81.696000,84.764000,86.541000 /)
XPIZA_LKT(84,1,1:6)=(/ 0.779845,0.897349,0.992418,0.993936,0.986700,0.511534 /)
XCGA_LKT(84,1,1:6)=(/ 0.879723,0.848690,0.822850,0.799500,0.792493,0.924897 /)
XEXT_COEFF_550_LKT(84,1)=78.061000 !rg=8.00573 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,2,1:6)=(/ 73.595000,74.477000,74.986000,76.160000,79.264000,82.535000 /)
XPIZA_LKT(84,2,1:6)=(/ 0.777536,0.894983,0.993309,0.994799,0.985422,0.511750 /)
XCGA_LKT(84,2,1:6)=(/ 0.879620,0.851780,0.821343,0.810367,0.787633,0.925117 /)
XEXT_COEFF_550_LKT(84,2)=74.585000 !rg=8.00573 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,3,1:6)=(/ 68.200000,69.148000,69.732000,71.216000,73.626000,76.169000 /)
XPIZA_LKT(84,3,1:6)=(/ 0.772738,0.893145,0.993389,0.993473,0.983893,0.513560 /)
XCGA_LKT(84,3,1:6)=(/ 0.880510,0.852783,0.822520,0.811700,0.788080,0.925833 /)
XEXT_COEFF_550_LKT(84,3)=68.970000 !rg=8.00573 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,4,1:6)=(/ 61.564000,62.390000,62.897000,63.946000,66.146000,68.511000 /)
XPIZA_LKT(84,4,1:6)=(/ 0.766727,0.891106,0.993632,0.993592,0.982845,0.515981 /)
XCGA_LKT(84,4,1:6)=(/ 0.881573,0.853597,0.823607,0.814293,0.793223,0.926657 /)
XEXT_COEFF_550_LKT(84,4)=62.173000 !rg=8.00573 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,5,1:6)=(/ 54.693000,54.970000,55.470000,56.513000,58.213000,60.346000 /)
XPIZA_LKT(84,5,1:6)=(/ 0.759688,0.888408,0.993931,0.993491,0.981584,0.518717 /)
XCGA_LKT(84,5,1:6)=(/ 0.883913,0.853710,0.824947,0.816220,0.800363,0.927507 /)
XEXT_COEFF_550_LKT(84,5)=55.358000 !rg=8.00573 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,6,1:6)=(/ 47.557000,48.107000,48.562000,49.764000,50.522000,52.365000 /)
XPIZA_LKT(84,6,1:6)=(/ 0.750141,0.884439,0.993985,0.993191,0.978899,0.521675 /)
XCGA_LKT(84,6,1:6)=(/ 0.885070,0.855157,0.825310,0.814927,0.803927,0.928403 /)
XEXT_COEFF_550_LKT(84,6)=48.166000 !rg=8.00573 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,7,1:6)=(/ 40.995000,41.393000,41.624000,42.126000,43.446000,44.816000 /)
XPIZA_LKT(84,7,1:6)=(/ 0.739626,0.879351,0.993917,0.993522,0.976875,0.524837 /)
XCGA_LKT(84,7,1:6)=(/ 0.887200,0.857150,0.828170,0.821233,0.807027,0.929397 /)
XEXT_COEFF_550_LKT(84,7)=41.414000 !rg=8.00573 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,8,1:6)=(/ 35.031000,35.371000,35.763000,36.199000,37.065000,38.055000 /)
XPIZA_LKT(84,8,1:6)=(/ 0.728022,0.872956,0.993699,0.993034,0.975352,0.527996 /)
XCGA_LKT(84,8,1:6)=(/ 0.889480,0.858903,0.828337,0.821583,0.810520,0.930397 /)
XEXT_COEFF_550_LKT(84,8)=35.343000 !rg=8.00573 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,9,1:6)=(/ 29.712000,29.868000,30.107000,30.552000,31.312000,32.053000 /)
XPIZA_LKT(84,9,1:6)=(/ 0.715607,0.864720,0.993556,0.993040,0.974030,0.531119 /)
XCGA_LKT(84,9,1:6)=(/ 0.892140,0.859830,0.831387,0.825090,0.816083,0.931393 /)
XEXT_COEFF_550_LKT(84,9)=29.983000 !rg=8.00573 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,10,1:6)=(/ 25.099000,25.210000,25.522000,25.694000,26.242000,26.899000 /)
XPIZA_LKT(84,10,1:6)=(/ 0.702504,0.855209,0.993115,0.992772,0.973316,0.534061 /)
XCGA_LKT(84,10,1:6)=(/ 0.894887,0.861797,0.833470,0.826490,0.818887,0.932310 /)
XEXT_COEFF_550_LKT(84,10)=25.305000 !rg=8.00573 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,11,1:6)=(/ 21.046000,21.158000,21.241000,21.526000,21.997000,22.440000 /)
XPIZA_LKT(84,11,1:6)=(/ 0.688306,0.844351,0.992525,0.992425,0.971849,0.536863 /)
XCGA_LKT(84,11,1:6)=(/ 0.897853,0.864003,0.834657,0.827400,0.822507,0.933167 /)
XEXT_COEFF_550_LKT(84,11)=21.204000 !rg=8.00573 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,12,1:6)=(/ 17.610000,17.666000,17.834000,18.014000,18.289000,18.659000 /)
XPIZA_LKT(84,12,1:6)=(/ 0.674510,0.831500,0.992002,0.991949,0.971119,0.539453 /)
XCGA_LKT(84,12,1:6)=(/ 0.901220,0.866197,0.836670,0.830633,0.826407,0.933920 /)
XEXT_COEFF_550_LKT(84,12)=17.740000 !rg=8.00573 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,13,1:6)=(/ 14.676000,14.748000,14.833000,14.981000,15.162000,15.479000 /)
XPIZA_LKT(84,13,1:6)=(/ 0.660315,0.817597,0.991148,0.991689,0.969752,0.541814 /)
XCGA_LKT(84,13,1:6)=(/ 0.904673,0.868660,0.837420,0.833443,0.829260,0.934573 /)
XEXT_COEFF_550_LKT(84,13)=14.762000 !rg=8.00573 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,14,1:6)=(/ 12.234000,12.278000,12.342000,12.446000,12.649000,12.842000 /)
XPIZA_LKT(84,14,1:6)=(/ 0.647053,0.802027,0.990089,0.990962,0.968149,0.543910 /)
XCGA_LKT(84,14,1:6)=(/ 0.908280,0.871253,0.838517,0.834263,0.832020,0.935113 /)
XEXT_COEFF_550_LKT(84,14)=12.295000 !rg=8.00573 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,15,1:6)=(/ 10.194000,10.220000,10.282000,10.367000,10.505000,10.646000 /)
XPIZA_LKT(84,15,1:6)=(/ 0.634623,0.785522,0.988865,0.990045,0.966136,0.545774 /)
XCGA_LKT(84,15,1:6)=(/ 0.912130,0.873860,0.840010,0.836187,0.834470,0.935567 /)
XEXT_COEFF_550_LKT(84,15)=10.256000 !rg=8.00573 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,16,1:6)=(/ 8.464900,8.494600,8.520700,8.584800,8.673800,8.807100 /)
XPIZA_LKT(84,16,1:6)=(/ 0.622556,0.768181,0.987374,0.988904,0.963795,0.547440 /)
XCGA_LKT(84,16,1:6)=(/ 0.915980,0.877063,0.841087,0.837133,0.836187,0.935950 /)
XEXT_COEFF_550_LKT(84,16)=8.504300 !rg=8.00573 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,17,1:6)=(/ 7.041400,7.061200,7.087200,7.129800,7.209700,7.299000 /)
XPIZA_LKT(84,17,1:6)=(/ 0.611876,0.750211,0.985532,0.987733,0.961005,0.548868 /)
XCGA_LKT(84,17,1:6)=(/ 0.000000,0.880210,0.841803,0.838397,0.838797,0.936253 /)
XEXT_COEFF_550_LKT(84,17)=7.069700 !rg=8.00573 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,18,1:6)=(/ 5.860200,5.872700,5.893500,5.938500,5.991800,6.051600 /)
XPIZA_LKT(84,18,1:6)=(/ 0.602259,0.732282,0.983389,0.986150,0.957320,0.550102 /)
XCGA_LKT(84,18,1:6)=(/ 0.000000,0.883587,0.842900,0.840367,0.840870,0.936497 /)
XEXT_COEFF_550_LKT(84,18)=5.882900 !rg=8.00573 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,19,1:6)=(/ 4.867300,4.879100,4.889900,4.917200,4.957600,5.010900 /)
XPIZA_LKT(84,19,1:6)=(/ 0.593471,0.714385,0.980805,0.984221,0.952971,0.551184 /)
XCGA_LKT(84,19,1:6)=(/ 0.000000,0.000000,0.843777,0.841067,0.842430,0.936693 /)
XEXT_COEFF_550_LKT(84,19)=4.880700 !rg=8.00573 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,20,1:6)=(/ 4.045300,4.053000,4.063900,4.079600,4.114900,4.152700 /)
XPIZA_LKT(84,20,1:6)=(/ 0.585720,0.696733,0.977722,0.981925,0.947976,0.552111 /)
XCGA_LKT(84,20,1:6)=(/ 0.000000,0.000000,0.844537,0.842013,0.844607,0.936843 /)
XEXT_COEFF_550_LKT(84,20)=4.055300 !rg=8.00573 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,1,1:6)=(/ 71.203000,71.628000,70.857000,75.422000,80.344000,79.135000 /)
XPIZA_LKT(85,1,1:6)=(/ 0.775950,0.897012,0.993041,0.993842,0.984716,0.513167 /)
XCGA_LKT(85,1,1:6)=(/ 0.880210,0.850457,0.826907,0.802630,0.801380,0.926340 /)
XEXT_COEFF_550_LKT(85,1)=73.056000 !rg=8.67306 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,2,1:6)=(/ 67.874000,68.535000,69.427000,70.828000,72.459000,75.649000 /)
XPIZA_LKT(85,2,1:6)=(/ 0.773319,0.895082,0.993119,0.993892,0.984646,0.514476 /)
XCGA_LKT(85,2,1:6)=(/ 0.880993,0.850883,0.820940,0.808287,0.795450,0.926670 /)
XEXT_COEFF_550_LKT(85,2)=68.898000 !rg=8.67306 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,3,1:6)=(/ 62.849000,63.331000,64.264000,65.906000,67.112000,69.831000 /)
XPIZA_LKT(85,3,1:6)=(/ 0.768712,0.893313,0.993874,0.993781,0.983784,0.516252 /)
XCGA_LKT(85,3,1:6)=(/ 0.881637,0.852470,0.824607,0.811227,0.798680,0.927227 /)
XEXT_COEFF_550_LKT(85,3)=63.732000 !rg=8.67306 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,4,1:6)=(/ 56.812000,57.229000,57.932000,59.457000,60.625000,62.822000 /)
XPIZA_LKT(85,4,1:6)=(/ 0.762501,0.890295,0.993912,0.993192,0.981696,0.518550 /)
XCGA_LKT(85,4,1:6)=(/ 0.883053,0.853270,0.825777,0.811817,0.801763,0.927947 /)
XEXT_COEFF_550_LKT(85,4)=57.408000 !rg=8.67306 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,5,1:6)=(/ 50.215000,50.848000,51.243000,51.918000,53.668000,55.354000 /)
XPIZA_LKT(85,5,1:6)=(/ 0.753519,0.886474,0.993955,0.993001,0.979236,0.521161 /)
XCGA_LKT(85,5,1:6)=(/ 0.884227,0.855277,0.826607,0.817707,0.800410,0.928727 /)
XEXT_COEFF_550_LKT(85,5)=50.918000 !rg=8.67306 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,6,1:6)=(/ 43.944000,44.244000,44.689000,45.448000,46.745000,48.051000 /)
XPIZA_LKT(85,6,1:6)=(/ 0.744938,0.882020,0.994027,0.992787,0.977573,0.523973 /)
XCGA_LKT(85,6,1:6)=(/ 0.886533,0.856260,0.827963,0.817857,0.805667,0.929543 /)
XEXT_COEFF_550_LKT(85,6)=44.611000 !rg=8.67306 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,7,1:6)=(/ 37.913000,38.020000,38.468000,39.048000,40.070000,41.139000 /)
XPIZA_LKT(85,7,1:6)=(/ 0.734417,0.875947,0.993859,0.993208,0.975963,0.526965 /)
XCGA_LKT(85,7,1:6)=(/ 0.888807,0.856857,0.829980,0.820600,0.811740,0.930427 /)
XEXT_COEFF_550_LKT(85,7)=38.265000 !rg=8.67306 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,8,1:6)=(/ 32.324000,32.472000,32.874000,33.273000,33.949000,34.946000 /)
XPIZA_LKT(85,8,1:6)=(/ 0.721903,0.869164,0.993562,0.993213,0.975127,0.529950 /)
XCGA_LKT(85,8,1:6)=(/ 0.890883,0.859047,0.830260,0.823727,0.815387,0.931310 /)
XEXT_COEFF_550_LKT(85,8)=32.663000 !rg=8.67306 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,9,1:6)=(/ 27.375000,27.583000,27.725000,28.033000,28.774000,29.445000 /)
XPIZA_LKT(85,9,1:6)=(/ 0.708747,0.860517,0.993288,0.992597,0.972996,0.532897 /)
XCGA_LKT(85,9,1:6)=(/ 0.893257,0.861207,0.832567,0.825833,0.817083,0.932183 /)
XEXT_COEFF_550_LKT(85,9)=27.599000 !rg=8.67306 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,10,1:6)=(/ 23.139000,23.273000,23.414000,23.753000,24.332000,24.720000 /)
XPIZA_LKT(85,10,1:6)=(/ 0.695584,0.850546,0.992912,0.992562,0.971740,0.535672 /)
XCGA_LKT(85,10,1:6)=(/ 0.896383,0.863000,0.833893,0.826797,0.819670,0.932983 /)
XEXT_COEFF_550_LKT(85,10)=23.294000 !rg=8.67306 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,11,1:6)=(/ 19.407000,19.502000,19.642000,19.824000,20.287000,20.629000 /)
XPIZA_LKT(85,11,1:6)=(/ 0.681441,0.838587,0.992416,0.992457,0.970934,0.538312 /)
XCGA_LKT(85,11,1:6)=(/ 0.899453,0.864980,0.835830,0.830543,0.824783,0.933727 /)
XEXT_COEFF_550_LKT(85,11)=19.589000 !rg=8.67306 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,12,1:6)=(/ 16.231000,16.307000,16.400000,16.572000,16.882000,17.159000 /)
XPIZA_LKT(85,12,1:6)=(/ 0.667493,0.825465,0.991652,0.991989,0.970021,0.540746 /)
XCGA_LKT(85,12,1:6)=(/ 0.902843,0.867360,0.837093,0.832170,0.826413,0.934383 /)
XEXT_COEFF_550_LKT(85,12)=16.326000 !rg=8.67306 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,13,1:6)=(/ 13.535000,13.589000,13.640000,13.814000,14.002000,14.240000 /)
XPIZA_LKT(85,13,1:6)=(/ 0.653806,0.810708,0.990736,0.991139,0.968919,0.542961 /)
XCGA_LKT(85,13,1:6)=(/ 0.906353,0.869877,0.837910,0.832827,0.830607,0.934950 /)
XEXT_COEFF_550_LKT(85,13)=13.585000 !rg=8.67306 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,14,1:6)=(/ 11.289000,11.322000,11.402000,11.523000,11.648000,11.818000 /)
XPIZA_LKT(85,14,1:6)=(/ 0.640939,0.794569,0.989628,0.990647,0.967256,0.544922 /)
XCGA_LKT(85,14,1:6)=(/ 0.910173,0.872403,0.839590,0.835630,0.832440,0.935423 /)
XEXT_COEFF_550_LKT(85,14)=11.366000 !rg=8.67306 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,15,1:6)=(/ 9.400600,9.434400,9.477700,9.542200,9.668500,9.799800 /)
XPIZA_LKT(85,15,1:6)=(/ 0.628688,0.777970,0.988268,0.989535,0.965258,0.546660 /)
XCGA_LKT(85,15,1:6)=(/ 0.913893,0.875370,0.840380,0.836587,0.834480,0.935817 /)
XEXT_COEFF_550_LKT(85,15)=9.440300 !rg=8.67306 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,16,1:6)=(/ 7.809000,7.830600,7.871700,7.907600,8.020100,8.109000 /)
XPIZA_LKT(85,16,1:6)=(/ 0.617223,0.759980,0.986599,0.988496,0.962734,0.548210 /)
XCGA_LKT(85,16,1:6)=(/ 0.917773,0.878410,0.841493,0.838403,0.838140,0.936147 /)
XEXT_COEFF_550_LKT(85,16)=7.840400 !rg=8.67306 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,17,1:6)=(/ 6.497700,6.510700,6.539900,6.581300,6.637700,6.722100 /)
XPIZA_LKT(85,17,1:6)=(/ 0.607035,0.741955,0.984644,0.987058,0.959438,0.549535 /)
XCGA_LKT(85,17,1:6)=(/ 0.000000,0.881660,0.842447,0.839530,0.839510,0.936410 /)
XEXT_COEFF_550_LKT(85,17)=6.519100 !rg=8.67306 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,18,1:6)=(/ 5.405700,5.418400,5.438100,5.466900,5.517400,5.574500 /)
XPIZA_LKT(85,18,1:6)=(/ 0.597833,0.724156,0.982305,0.985327,0.955647,0.550678 /)
XCGA_LKT(85,18,1:6)=(/ 0.000000,0.885300,0.843143,0.840500,0.841290,0.936617 /)
XEXT_COEFF_550_LKT(85,18)=5.420400 !rg=8.67306 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,19,1:6)=(/ 4.491000,4.498900,4.513600,4.533200,4.576400,4.616800 /)
XPIZA_LKT(85,19,1:6)=(/ 0.589570,0.706128,0.979529,0.983247,0.950996,0.551679 /)
XCGA_LKT(85,19,1:6)=(/ 0.000000,0.000000,0.844177,0.841723,0.843860,0.936783 /)
XEXT_COEFF_550_LKT(85,19)=4.502800 !rg=8.67306 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,20,1:6)=(/ 3.732900,3.738400,3.749500,3.766200,3.790800,3.826800 /)
XPIZA_LKT(85,20,1:6)=(/ 0.582357,0.688935,0.976211,0.980789,0.945285,0.552534 /)
XCGA_LKT(85,20,1:6)=(/ 0.000000,0.000000,0.845030,0.842850,0.845173,0.936910 /)
XEXT_COEFF_550_LKT(85,20)=3.741600 !rg=8.67306 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,1,1:6)=(/ 65.730000,65.769000,66.658000,68.902000,68.190000,72.499000 /)
XPIZA_LKT(86,1,1:6)=(/ 0.772204,0.896637,0.993178,0.990931,0.983472,0.515652 /)
XCGA_LKT(86,1,1:6)=(/ 0.883160,0.850767,0.822157,0.801177,0.797753,0.927697 /)
XEXT_COEFF_550_LKT(86,1)=65.499000 !rg=9.39601 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,2,1:6)=(/ 62.393000,63.093000,63.578000,64.714000,67.573000,69.353000 /)
XPIZA_LKT(86,2,1:6)=(/ 0.767954,0.895022,0.994065,0.993609,0.983197,0.517135 /)
XCGA_LKT(86,2,1:6)=(/ 0.881557,0.852727,0.824717,0.818440,0.795097,0.927993 /)
XEXT_COEFF_550_LKT(86,2)=63.276000 !rg=9.39601 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,3,1:6)=(/ 57.992000,58.395000,59.258000,59.764000,62.274000,64.041000 /)
XPIZA_LKT(86,3,1:6)=(/ 0.763760,0.892439,0.993838,0.993437,0.982248,0.518875 /)
XCGA_LKT(86,3,1:6)=(/ 0.882800,0.853447,0.824527,0.815383,0.799850,0.928447 /)
XEXT_COEFF_550_LKT(86,3)=58.595000 !rg=9.39601 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,4,1:6)=(/ 52.316000,52.746000,53.394000,53.854000,56.209000,57.626000 /)
XPIZA_LKT(86,4,1:6)=(/ 0.756708,0.889136,0.994024,0.993300,0.980136,0.521044 /)
XCGA_LKT(86,4,1:6)=(/ 0.884080,0.854490,0.826790,0.817040,0.801697,0.929067 /)
XEXT_COEFF_550_LKT(86,4)=52.916000 !rg=9.39601 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,5,1:6)=(/ 46.398000,46.637000,47.260000,48.133000,49.249000,50.791000 /)
XPIZA_LKT(86,5,1:6)=(/ 0.748602,0.884430,0.994045,0.993120,0.978862,0.523524 /)
XCGA_LKT(86,5,1:6)=(/ 0.885867,0.855403,0.827473,0.817750,0.808593,0.929773 /)
XEXT_COEFF_550_LKT(86,5)=46.941000 !rg=9.39601 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,6,1:6)=(/ 40.475000,40.856000,41.329000,42.049000,42.909000,44.105000 /)
XPIZA_LKT(86,6,1:6)=(/ 0.738848,0.879437,0.993890,0.993035,0.976803,0.526183 /)
XCGA_LKT(86,6,1:6)=(/ 0.887587,0.857470,0.826947,0.819147,0.809373,0.930517 /)
XEXT_COEFF_550_LKT(86,6)=40.842000 !rg=9.39601 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,7,1:6)=(/ 34.851000,35.079000,35.256000,35.850000,36.929000,37.774000 /)
XPIZA_LKT(86,7,1:6)=(/ 0.727133,0.873011,0.993632,0.992902,0.974555,0.529009 /)
XCGA_LKT(86,7,1:6)=(/ 0.889760,0.858433,0.830823,0.821090,0.812727,0.931307 /)
XEXT_COEFF_550_LKT(86,7)=35.232000 !rg=9.39601 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,8,1:6)=(/ 29.794000,29.947000,30.238000,30.659000,31.299000,32.099000 /)
XPIZA_LKT(86,8,1:6)=(/ 0.715076,0.865322,0.993562,0.992982,0.973760,0.531823 /)
XCGA_LKT(86,8,1:6)=(/ 0.892140,0.860123,0.832073,0.825653,0.816647,0.932090 /)
XEXT_COEFF_550_LKT(86,8)=30.070000 !rg=9.39601 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,9,1:6)=(/ 25.276000,25.360000,25.625000,25.962000,26.493000,27.056000 /)
XPIZA_LKT(86,9,1:6)=(/ 0.702302,0.855653,0.993123,0.992869,0.972557,0.534600 /)
XCGA_LKT(86,9,1:6)=(/ 0.895010,0.861797,0.832970,0.826713,0.821350,0.932857 /)
XEXT_COEFF_550_LKT(86,9)=25.503000 !rg=9.39601 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,10,1:6)=(/ 21.315000,21.433000,21.601000,21.820000,22.200000,22.722000 /)
XPIZA_LKT(86,10,1:6)=(/ 0.688417,0.845109,0.992763,0.992490,0.971949,0.537209 /)
XCGA_LKT(86,10,1:6)=(/ 0.897773,0.863780,0.835233,0.828080,0.823737,0.933560 /)
XEXT_COEFF_550_LKT(86,10)=21.465000 !rg=9.39601 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,11,1:6)=(/ 17.908000,17.963000,18.067000,18.328000,18.601000,18.969000 /)
XPIZA_LKT(86,11,1:6)=(/ 0.674770,0.832692,0.992121,0.992247,0.970462,0.539688 /)
XCGA_LKT(86,11,1:6)=(/ 0.901190,0.865937,0.836817,0.831623,0.826703,0.934210 /)
XEXT_COEFF_550_LKT(86,11)=18.049000 !rg=9.39601 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,12,1:6)=(/ 14.968000,15.026000,15.126000,15.265000,15.515000,15.783000 /)
XPIZA_LKT(86,12,1:6)=(/ 0.660896,0.818686,0.991266,0.991713,0.969663,0.541972 /)
XCGA_LKT(86,12,1:6)=(/ 0.904463,0.868480,0.837563,0.832617,0.829803,0.934783 /)
XEXT_COEFF_550_LKT(86,12)=15.041000 !rg=9.39601 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,13,1:6)=(/ 12.494000,12.526000,12.600000,12.734000,12.884000,13.102000 /)
XPIZA_LKT(86,13,1:6)=(/ 0.647623,0.803431,0.990218,0.991083,0.968158,0.544043 /)
XCGA_LKT(86,13,1:6)=(/ 0.908270,0.870840,0.838900,0.834693,0.831023,0.935280 /)
XEXT_COEFF_550_LKT(86,13)=12.567000 !rg=9.39601 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,14,1:6)=(/ 10.412000,10.455000,10.489000,10.568000,10.697000,10.877000 /)
XPIZA_LKT(86,14,1:6)=(/ 0.634797,0.787316,0.989042,0.990176,0.966421,0.545872 /)
XCGA_LKT(86,14,1:6)=(/ 0.911920,0.873920,0.839493,0.836120,0.834230,0.935690 /)
XEXT_COEFF_550_LKT(86,14)=10.465000 !rg=9.39601 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,15,1:6)=(/ 8.671800,8.694700,8.753600,8.780000,8.899900,9.021900 /)
XPIZA_LKT(86,15,1:6)=(/ 0.623049,0.769801,0.987612,0.989111,0.964212,0.547491 /)
XCGA_LKT(86,15,1:6)=(/ 0.915717,0.876647,0.841207,0.837470,0.836607,0.936030 /)
XEXT_COEFF_550_LKT(86,15)=8.716400 !rg=9.39601 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,16,1:6)=(/ 7.204900,7.224800,7.250500,7.308900,7.367700,7.467200 /)
XPIZA_LKT(86,16,1:6)=(/ 0.612175,0.751911,0.985808,0.987889,0.961272,0.548930 /)
XCGA_LKT(86,16,1:6)=(/ 0.000000,0.879917,0.842190,0.839483,0.839030,0.936317 /)
XEXT_COEFF_550_LKT(86,16)=7.234800 !rg=9.39601 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,17,1:6)=(/ 5.994500,6.011600,6.025900,6.066400,6.110500,6.191400 /)
XPIZA_LKT(86,17,1:6)=(/ 0.602396,0.733969,0.983649,0.986368,0.957927,0.550157 /)
XCGA_LKT(86,17,1:6)=(/ 0.000000,0.883390,0.842593,0.840087,0.840987,0.936543 /)
XEXT_COEFF_550_LKT(86,17)=6.014900 !rg=9.39601 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,18,1:6)=(/ 4.988300,4.997500,5.020200,5.033100,5.087200,5.135500 /)
XPIZA_LKT(86,18,1:6)=(/ 0.593723,0.715904,0.981164,0.984459,0.953768,0.551213 /)
XCGA_LKT(86,18,1:6)=(/ 0.000000,0.000000,0.843777,0.841177,0.842877,0.936720 /)
XEXT_COEFF_550_LKT(86,18)=5.005100 !rg=9.39601 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,19,1:6)=(/ 4.143700,4.151100,4.163500,4.183100,4.211600,4.254100 /)
XPIZA_LKT(86,19,1:6)=(/ 0.585916,0.698259,0.978137,0.982242,0.948628,0.552136 /)
XCGA_LKT(86,19,1:6)=(/ 0.000000,0.000000,0.844573,0.842553,0.844587,0.936863 /)
XEXT_COEFF_550_LKT(86,19)=4.154800 !rg=9.39601 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,20,1:6)=(/ 0.000000,3.450800,3.457600,3.472100,3.493300,3.526800 /)
XPIZA_LKT(86,20,1:6)=(/ 0.900,0.681465,0.974592,0.979579,0.942658,0.552922 /)
XCGA_LKT(86,20,1:6)=(/ 0.900,0.000000,0.845273,0.843267,0.846260,0.936967 /)
XEXT_COEFF_550_LKT(86,20)=3.453300 !rg=9.39601 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,1,1:6)=(/ 60.116000,60.639000,61.057000,62.434000,63.800000,66.550000 /)
XPIZA_LKT(87,1,1:6)=(/ 0.765444,0.894168,0.994422,0.993517,0.983387,0.518882 /)
XCGA_LKT(87,1,1:6)=(/ 0.882637,0.854343,0.827827,0.816713,0.808517,0.928957 /)
XEXT_COEFF_550_LKT(87,1)=60.867000 !rg=10.1792 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,2,1:6)=(/ 57.752000,58.038000,58.806000,60.106000,62.485000,63.607000 /)
XPIZA_LKT(87,2,1:6)=(/ 0.764000,0.892814,0.993668,0.993010,0.980868,0.519774 /)
XCGA_LKT(87,2,1:6)=(/ 0.883073,0.854217,0.825163,0.813623,0.799280,0.929110 /)
XEXT_COEFF_550_LKT(87,2)=58.590000 !rg=10.1792 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,3,1:6)=(/ 53.350000,53.929000,54.662000,55.089000,57.538000,58.748000 /)
XPIZA_LKT(87,3,1:6)=(/ 0.758151,0.890543,0.994197,0.992805,0.979554,0.521420 /)
XCGA_LKT(87,3,1:6)=(/ 0.883827,0.853803,0.828367,0.818007,0.800170,0.929503 /)
XEXT_COEFF_550_LKT(87,3)=54.034000 !rg=10.1792 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,4,1:6)=(/ 48.222000,48.674000,49.316000,49.655000,51.643000,52.876000 /)
XPIZA_LKT(87,4,1:6)=(/ 0.750978,0.886961,0.993815,0.992933,0.978322,0.523457 /)
XCGA_LKT(87,4,1:6)=(/ 0.885157,0.854837,0.828263,0.820073,0.803990,0.930047 /)
XEXT_COEFF_550_LKT(87,4)=48.824000 !rg=10.1792 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,5,1:6)=(/ 42.745000,43.138000,43.598000,43.796000,45.588000,46.618000 /)
XPIZA_LKT(87,5,1:6)=(/ 0.742175,0.882442,0.993970,0.993152,0.976739,0.525793 /)
XCGA_LKT(87,5,1:6)=(/ 0.886933,0.856717,0.828770,0.821327,0.808017,0.930683 /)
XEXT_COEFF_550_LKT(87,5)=43.160000 !rg=10.1792 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,6,1:6)=(/ 37.295000,37.553000,37.898000,38.692000,39.620000,40.494000 /)
XPIZA_LKT(87,6,1:6)=(/ 0.731831,0.876254,0.993941,0.993099,0.974644,0.528304 /)
XCGA_LKT(87,6,1:6)=(/ 0.888720,0.857943,0.830413,0.820567,0.811347,0.931353 /)
XEXT_COEFF_550_LKT(87,6)=37.748000 !rg=10.1792 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,7,1:6)=(/ 32.152000,32.339000,32.661000,32.984000,33.973000,34.693000 /)
XPIZA_LKT(87,7,1:6)=(/ 0.720579,0.869215,0.993766,0.993124,0.973198,0.530966 /)
XCGA_LKT(87,7,1:6)=(/ 0.891223,0.859340,0.831550,0.825260,0.815420,0.932063 /)
XEXT_COEFF_550_LKT(87,7)=32.545000 !rg=10.1792 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,8,1:6)=(/ 27.460000,27.654000,27.825000,28.232000,28.809000,29.491000 /)
XPIZA_LKT(87,8,1:6)=(/ 0.708162,0.861056,0.993430,0.992865,0.972890,0.533613 /)
XCGA_LKT(87,8,1:6)=(/ 0.893577,0.861260,0.832870,0.826320,0.819497,0.932760 /)
XEXT_COEFF_550_LKT(87,8)=27.666000 !rg=10.1792 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,9,1:6)=(/ 23.275000,23.474000,23.553000,23.760000,24.471000,24.866000 /)
XPIZA_LKT(87,9,1:6)=(/ 0.695029,0.851209,0.992925,0.992823,0.971808,0.536220 /)
XCGA_LKT(87,9,1:6)=(/ 0.896313,0.863407,0.833853,0.828937,0.820883,0.933440 /)
XEXT_COEFF_550_LKT(87,9)=23.446000 !rg=10.1792 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,10,1:6)=(/ 19.686000,19.749000,19.890000,20.128000,20.565000,20.891000 /)
XPIZA_LKT(87,10,1:6)=(/ 0.681906,0.839501,0.992423,0.992375,0.970674,0.538669 /)
XCGA_LKT(87,10,1:6)=(/ 0.899477,0.865013,0.836327,0.830797,0.824740,0.934057 /)
XEXT_COEFF_550_LKT(87,10)=19.844000 !rg=10.1792 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,11,1:6)=(/ 16.506000,16.579000,16.661000,16.794000,17.123000,17.445000 /)
XPIZA_LKT(87,11,1:6)=(/ 0.667745,0.826411,0.991705,0.992092,0.970185,0.540991 /)
XCGA_LKT(87,11,1:6)=(/ 0.902680,0.867347,0.836637,0.832747,0.828313,0.934627 /)
XEXT_COEFF_550_LKT(87,11)=16.592000 !rg=10.1792 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,12,1:6)=(/ 13.809000,13.853000,13.958000,14.116000,14.306000,14.520000 /)
XPIZA_LKT(87,12,1:6)=(/ 0.654287,0.811877,0.990843,0.991376,0.968860,0.543127 /)
XCGA_LKT(87,12,1:6)=(/ 0.906290,0.869623,0.838680,0.833983,0.830290,0.935130 /)
XEXT_COEFF_550_LKT(87,12)=13.926000 !rg=10.1792 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,13,1:6)=(/ 11.518000,11.566000,11.617000,11.708000,11.852000,12.057000 /)
XPIZA_LKT(87,13,1:6)=(/ 0.641204,0.796293,0.989775,0.990671,0.967602,0.545060 /)
XCGA_LKT(87,13,1:6)=(/ 0.909900,0.872433,0.839230,0.835683,0.833690,0.935560 /)
XEXT_COEFF_550_LKT(87,13)=11.582000 !rg=10.1792 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,14,1:6)=(/ 9.605400,9.628400,9.689000,9.755800,9.862200,10.012000 /)
XPIZA_LKT(87,14,1:6)=(/ 0.628962,0.779268,0.988390,0.989838,0.965571,0.546763 /)
XCGA_LKT(87,14,1:6)=(/ 0.913737,0.875067,0.840420,0.836980,0.835900,0.935920 /)
XEXT_COEFF_550_LKT(87,14)=9.664000 !rg=10.1792 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,15,1:6)=(/ 8.000900,8.028900,8.056600,8.122000,8.204200,8.306800 /)
XPIZA_LKT(87,15,1:6)=(/ 0.617701,0.762013,0.986820,0.988648,0.962896,0.548266 /)
XCGA_LKT(87,15,1:6)=(/ 0.917563,0.878263,0.841597,0.838240,0.837270,0.936217 /)
XEXT_COEFF_550_LKT(87,15)=8.031800 !rg=10.1792 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,16,1:6)=(/ 6.648100,6.666500,6.683000,6.725400,6.800000,6.877000 /)
XPIZA_LKT(87,16,1:6)=(/ 0.607272,0.743762,0.984905,0.987231,0.959721,0.549600 /)
XCGA_LKT(87,16,1:6)=(/ 0.000000,0.881543,0.842180,0.839527,0.839830,0.936463 /)
XEXT_COEFF_550_LKT(87,16)=6.668600 !rg=10.1792 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,17,1:6)=(/ 5.531700,5.543400,5.567300,5.584200,5.652500,5.703300 /)
XPIZA_LKT(87,17,1:6)=(/ 0.598048,0.725656,0.982628,0.985577,0.956217,0.550734 /)
XCGA_LKT(87,17,1:6)=(/ 0.000000,0.884997,0.843420,0.840580,0.842260,0.936657 /)
XEXT_COEFF_550_LKT(87,17)=5.556300 !rg=10.1792 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,18,1:6)=(/ 4.602600,4.612000,4.624900,4.649800,4.687400,4.731600 /)
XPIZA_LKT(87,18,1:6)=(/ 0.589807,0.707937,0.979923,0.983563,0.951510,0.551708 /)
XCGA_LKT(87,18,1:6)=(/ 0.000000,0.000000,0.844170,0.841753,0.843053,0.936810 /)
XEXT_COEFF_550_LKT(87,18)=4.615300 !rg=10.1792 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,19,1:6)=(/ 3.823800,3.830500,3.837800,3.854600,3.887100,3.920200 /)
XPIZA_LKT(87,19,1:6)=(/ 0.582529,0.690476,0.976652,0.981160,0.946077,0.552558 /)
XCGA_LKT(87,19,1:6)=(/ 0.000000,0.000000,0.844933,0.842817,0.845267,0.936930 /)
XEXT_COEFF_550_LKT(87,19)=3.832500 !rg=10.1792 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,20,1:6)=(/ 0.000000,3.183800,3.190600,3.200500,3.229400,3.250500 /)
XPIZA_LKT(87,20,1:6)=(/ 0.900,0.673954,0.972856,0.978258,0.939832,0.553278 /)
XCGA_LKT(87,20,1:6)=(/ 0.900,0.000000,0.000000,0.843697,0.847407,0.937013 /)
XEXT_COEFF_550_LKT(87,20)=3.187500 !rg=10.1792 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,1,1:6)=(/ 55.551000,56.081000,56.729000,56.570000,59.983000,61.059000 /)
XPIZA_LKT(88,1,1:6)=(/ 0.760717,0.892689,0.994255,0.993101,0.981055,0.521607 /)
XCGA_LKT(88,1,1:6)=(/ 0.882867,0.855393,0.827213,0.820560,0.803997,0.929883 /)
XEXT_COEFF_550_LKT(88,1)=55.980000 !rg=11.0277 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,2,1:6)=(/ 53.055000,53.702000,54.024000,55.569000,56.430000,58.358000 /)
XPIZA_LKT(88,2,1:6)=(/ 0.757252,0.890716,0.994438,0.992357,0.981331,0.522350 /)
XCGA_LKT(88,2,1:6)=(/ 0.883893,0.854780,0.827307,0.817713,0.808060,0.930060 /)
XEXT_COEFF_550_LKT(88,2)=53.606000 !rg=11.0277 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,3,1:6)=(/ 49.412000,49.836000,50.033000,50.903000,52.039000,53.909000 /)
XPIZA_LKT(88,3,1:6)=(/ 0.752928,0.888336,0.994323,0.992722,0.979731,0.523861 /)
XCGA_LKT(88,3,1:6)=(/ 0.885487,0.855657,0.828583,0.821493,0.808713,0.930427 /)
XEXT_COEFF_550_LKT(88,3)=49.901000 !rg=11.0277 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,4,1:6)=(/ 44.564000,44.771000,45.171000,45.774000,47.007000,48.532000 /)
XPIZA_LKT(88,4,1:6)=(/ 0.745197,0.884203,0.994198,0.992889,0.977387,0.525769 /)
XCGA_LKT(88,4,1:6)=(/ 0.886757,0.856003,0.828877,0.822950,0.810710,0.930907 /)
XEXT_COEFF_550_LKT(88,4)=45.061000 !rg=11.0277 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,5,1:6)=(/ 39.364000,39.661000,40.147000,40.530000,41.846000,42.799000 /)
XPIZA_LKT(88,5,1:6)=(/ 0.735743,0.879057,0.994090,0.992721,0.974908,0.527965 /)
XCGA_LKT(88,5,1:6)=(/ 0.888020,0.856917,0.830917,0.822030,0.810097,0.931473 /)
XEXT_COEFF_550_LKT(88,5)=39.855000 !rg=11.0277 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,6,1:6)=(/ 34.455000,34.555000,35.012000,35.321000,36.066000,37.188000 /)
XPIZA_LKT(88,6,1:6)=(/ 0.725958,0.872506,0.993902,0.993160,0.974602,0.530326 /)
XCGA_LKT(88,6,1:6)=(/ 0.890323,0.858840,0.830777,0.823890,0.814343,0.932080 /)
XEXT_COEFF_550_LKT(88,6)=34.725000 !rg=11.0277 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,7,1:6)=(/ 29.669000,29.751000,29.979000,30.470000,31.125000,31.871000 /)
XPIZA_LKT(88,7,1:6)=(/ 0.713998,0.865231,0.993686,0.992875,0.972644,0.532828 /)
XCGA_LKT(88,7,1:6)=(/ 0.892793,0.860033,0.833370,0.826103,0.818787,0.932720 /)
XEXT_COEFF_550_LKT(88,7)=29.970000 !rg=11.0277 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,8,1:6)=(/ 25.363000,25.435000,25.723000,26.134000,26.517000,27.101000 /)
XPIZA_LKT(88,8,1:6)=(/ 0.701501,0.855991,0.993275,0.992669,0.972050,0.535312 /)
XCGA_LKT(88,8,1:6)=(/ 0.895370,0.861947,0.834493,0.827503,0.819627,0.933340 /)
XEXT_COEFF_550_LKT(88,8)=25.628000 !rg=11.0277 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,9,1:6)=(/ 21.468000,21.567000,21.775000,21.925000,22.436000,22.859000 /)
XPIZA_LKT(88,9,1:6)=(/ 0.687986,0.845387,0.992766,0.992564,0.970651,0.537758 /)
XCGA_LKT(88,9,1:6)=(/ 0.897860,0.863777,0.835797,0.829273,0.822937,0.933940 /)
XEXT_COEFF_550_LKT(88,9)=21.632000 !rg=11.0277 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,10,1:6)=(/ 18.138000,18.238000,18.327000,18.512000,18.865000,19.210000 /)
XPIZA_LKT(88,10,1:6)=(/ 0.674587,0.833717,0.992149,0.992219,0.970481,0.540049 /)
XCGA_LKT(88,10,1:6)=(/ 0.900977,0.866523,0.836520,0.831767,0.826100,0.934487 /)
XEXT_COEFF_550_LKT(88,10)=18.236000 !rg=11.0277 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,11,1:6)=(/ 15.234000,15.287000,15.412000,15.537000,15.771000,16.047000 /)
XPIZA_LKT(88,11,1:6)=(/ 0.661248,0.819936,0.991367,0.991800,0.969467,0.542221 /)
XCGA_LKT(88,11,1:6)=(/ 0.904450,0.868313,0.837803,0.833107,0.830650,0.934990 /)
XEXT_COEFF_550_LKT(88,11)=15.353000 !rg=11.0277 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,12,1:6)=(/ 12.736000,12.799000,12.851000,12.954000,13.125000,13.361000 /)
XPIZA_LKT(88,12,1:6)=(/ 0.647679,0.805069,0.990302,0.991111,0.968097,0.544213 /)
XCGA_LKT(88,12,1:6)=(/ 0.907997,0.871210,0.838337,0.834567,0.831607,0.935430 /)
XEXT_COEFF_550_LKT(88,12)=12.802000 !rg=11.0277 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,13,1:6)=(/ 10.627000,10.658000,10.729000,10.791000,10.932000,11.097000 /)
XPIZA_LKT(88,13,1:6)=(/ 0.635118,0.788481,0.989178,0.990326,0.966560,0.546013 /)
XCGA_LKT(88,13,1:6)=(/ 0.911777,0.873607,0.840230,0.835910,0.835417,0.935807 /)
XEXT_COEFF_550_LKT(88,13)=10.697000 !rg=11.0277 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,14,1:6)=(/ 8.861200,8.891900,8.928800,8.996600,9.077700,9.217300 /)
XPIZA_LKT(88,14,1:6)=(/ 0.623292,0.771608,0.987785,0.989287,0.964405,0.547596 /)
XCGA_LKT(88,14,1:6)=(/ 0.915550,0.876557,0.841010,0.838110,0.836650,0.936120 /)
XEXT_COEFF_550_LKT(88,14)=8.902600 !rg=11.0277 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,15,1:6)=(/ 7.380200,7.400400,7.431000,7.476500,7.538900,7.649400 /)
XPIZA_LKT(88,15,1:6)=(/ 0.612449,0.753620,0.986042,0.988061,0.961669,0.548989 /)
XCGA_LKT(88,15,1:6)=(/ 0.000000,0.879690,0.841943,0.838813,0.838563,0.936377 /)
XEXT_COEFF_550_LKT(88,15)=7.404900 !rg=11.0277 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,16,1:6)=(/ 6.134000,6.147200,6.171600,6.206600,6.265100,6.334100 /)
XPIZA_LKT(88,16,1:6)=(/ 0.602609,0.735447,0.983963,0.986600,0.958334,0.550222 /)
XCGA_LKT(88,16,1:6)=(/ 0.000000,0.883000,0.842977,0.840353,0.841013,0.936590 /)
XEXT_COEFF_550_LKT(88,16)=6.158200 !rg=11.0277 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,17,1:6)=(/ 5.103400,5.115900,5.126100,5.159000,5.186100,5.254200 /)
XPIZA_LKT(88,17,1:6)=(/ 0.593849,0.717548,0.981459,0.984692,0.954310,0.551268 /)
XCGA_LKT(88,17,1:6)=(/ 0.000000,0.886713,0.843460,0.841303,0.842793,0.936757 /)
XEXT_COEFF_550_LKT(88,17)=5.117400 !rg=11.0277 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,18,1:6)=(/ 4.246400,4.253900,4.264300,4.286200,4.312600,4.359800 /)
XPIZA_LKT(88,18,1:6)=(/ 0.586098,0.699801,0.978581,0.982565,0.949317,0.552166 /)
XCGA_LKT(88,18,1:6)=(/ 0.000000,0.000000,0.844503,0.842147,0.844237,0.936887 /)
XEXT_COEFF_550_LKT(88,18)=4.255300 !rg=11.0277 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,19,1:6)=(/ 3.528500,3.534100,3.542900,3.556600,3.585900,3.612800 /)
XPIZA_LKT(88,19,1:6)=(/ 0.579307,0.682779,0.975080,0.979943,0.943488,0.552946 /)
XCGA_LKT(88,19,1:6)=(/ 0.000000,0.000000,0.845420,0.843243,0.846247,0.936983 /)
XEXT_COEFF_550_LKT(88,19)=3.537200 !rg=11.0277 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,20,1:6)=(/ 0.000000,2.937700,2.941600,2.953800,2.971100,2.996100 /)
XPIZA_LKT(88,20,1:6)=(/ 0.900,0.666673,0.970962,0.976853,0.936660,0.553603 /)
XCGA_LKT(88,20,1:6)=(/ 0.900,0.000000,0.000000,0.844247,0.847880,0.937053 /)
XEXT_COEFF_550_LKT(88,20)=2.939200 !rg=11.0277 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,1,1:6)=(/ 51.464000,51.443000,53.238000,53.539000,54.071000,56.011000 /)
XPIZA_LKT(89,1,1:6)=(/ 0.756816,0.889756,0.993673,0.991737,0.975244,0.523942 /)
XCGA_LKT(89,1,1:6)=(/ 0.886263,0.854527,0.825217,0.818143,0.796827,0.930640 /)
XEXT_COEFF_550_LKT(89,1)=51.673000 !rg=11.9469 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,2,1:6)=(/ 49.254000,49.483000,50.315000,50.811000,51.783000,53.557000 /)
XPIZA_LKT(89,2,1:6)=(/ 0.752687,0.888548,0.994474,0.992017,0.979936,0.524796 /)
XCGA_LKT(89,2,1:6)=(/ 0.885953,0.855690,0.827593,0.818887,0.808330,0.930883 /)
XEXT_COEFF_550_LKT(89,2)=49.545000 !rg=11.9469 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,3,1:6)=(/ 45.322000,45.790000,46.299000,46.544000,48.061000,49.483000 /)
XPIZA_LKT(89,3,1:6)=(/ 0.745351,0.885523,0.994058,0.992464,0.976874,0.526193 /)
XCGA_LKT(89,3,1:6)=(/ 0.886117,0.856133,0.828513,0.822047,0.808510,0.931227 /)
XEXT_COEFF_550_LKT(89,3)=46.016000 !rg=11.9469 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,4,1:6)=(/ 40.968000,41.388000,41.692000,42.040000,43.320000,44.557000 /)
XPIZA_LKT(89,4,1:6)=(/ 0.738005,0.881435,0.994106,0.992166,0.975315,0.527975 /)
XCGA_LKT(89,4,1:6)=(/ 0.887530,0.857490,0.830367,0.822817,0.810913,0.931663 /)
XEXT_COEFF_550_LKT(89,4)=41.508000 !rg=11.9469 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,5,1:6)=(/ 36.382000,36.706000,36.864000,37.387000,38.327000,39.303000 /)
XPIZA_LKT(89,5,1:6)=(/ 0.729441,0.876113,0.994043,0.992746,0.974021,0.530031 /)
XCGA_LKT(89,5,1:6)=(/ 0.889657,0.858823,0.831040,0.824847,0.815800,0.932173 /)
XEXT_COEFF_550_LKT(89,5)=36.706000 !rg=11.9469 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,6,1:6)=(/ 31.704000,31.939000,32.119000,32.945000,33.123000,34.160000 /)
XPIZA_LKT(89,6,1:6)=(/ 0.718346,0.869122,0.993748,0.992658,0.973167,0.532246 /)
XCGA_LKT(89,6,1:6)=(/ 0.891490,0.859733,0.831727,0.822060,0.819050,0.932717 /)
XEXT_COEFF_550_LKT(89,6)=31.974000 !rg=11.9469 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,7,1:6)=(/ 27.313000,27.523000,27.650000,27.894000,28.592000,29.286000 /)
XPIZA_LKT(89,7,1:6)=(/ 0.706486,0.860911,0.993436,0.992965,0.972272,0.534593 /)
XCGA_LKT(89,7,1:6)=(/ 0.893923,0.861603,0.833077,0.828540,0.820643,0.933290 /)
XEXT_COEFF_550_LKT(89,7)=27.488000 !rg=11.9469 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,8,1:6)=(/ 23.358000,23.516000,23.618000,23.911000,24.321000,24.910000 /)
XPIZA_LKT(89,8,1:6)=(/ 0.693933,0.851364,0.993039,0.992590,0.971668,0.536920 /)
XCGA_LKT(89,8,1:6)=(/ 0.896670,0.863433,0.833770,0.829150,0.823493,0.933843 /)
XEXT_COEFF_550_LKT(89,8)=23.525000 !rg=11.9469 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,9,1:6)=(/ 19.812000,19.925000,20.020000,20.285000,20.651000,21.018000 /)
XPIZA_LKT(89,9,1:6)=(/ 0.680889,0.840170,0.992477,0.992520,0.970464,0.539208 /)
XCGA_LKT(89,9,1:6)=(/ 0.899603,0.865333,0.836093,0.830917,0.826310,0.934380 /)
XEXT_COEFF_550_LKT(89,9)=19.939000 !rg=11.9469 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,10,1:6)=(/ 16.741000,16.796000,16.967000,17.106000,17.345000,17.668000 /)
XPIZA_LKT(89,10,1:6)=(/ 0.667778,0.827213,0.991807,0.991968,0.970163,0.541350 /)
XCGA_LKT(89,10,1:6)=(/ 0.902803,0.867317,0.837747,0.832063,0.828837,0.934863 /)
XEXT_COEFF_550_LKT(89,10)=16.856000 !rg=11.9469 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,11,1:6)=(/ 14.043000,14.120000,14.154000,14.287000,14.518000,14.764000 /)
XPIZA_LKT(89,11,1:6)=(/ 0.654142,0.813358,0.990929,0.991428,0.968788,0.543375 /)
XCGA_LKT(89,11,1:6)=(/ 0.906197,0.869893,0.838277,0.833787,0.831227,0.935303 /)
XEXT_COEFF_550_LKT(89,11)=14.123000 !rg=11.9469 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,12,1:6)=(/ 11.749000,11.787000,11.865000,11.950000,12.098000,12.295000 /)
XPIZA_LKT(89,12,1:6)=(/ 0.641328,0.797522,0.989882,0.990882,0.967723,0.545231 /)
XCGA_LKT(89,12,1:6)=(/ 0.909817,0.872137,0.839677,0.835783,0.833817,0.935690 /)
XEXT_COEFF_550_LKT(89,12)=11.820000 !rg=11.9469 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,13,1:6)=(/ 9.800900,9.838400,9.875900,9.958000,10.045000,10.215000 /)
XPIZA_LKT(89,13,1:6)=(/ 0.629060,0.780816,0.988603,0.989908,0.965789,0.546904 /)
XCGA_LKT(89,13,1:6)=(/ 0.913603,0.875080,0.840410,0.837317,0.836127,0.936020 /)
XEXT_COEFF_550_LKT(89,13)=9.846900 !rg=11.9469 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,14,1:6)=(/ 8.173800,8.197100,8.223300,8.286100,8.386400,8.486900 /)
XPIZA_LKT(89,14,1:6)=(/ 0.617776,0.763321,0.987061,0.988812,0.963271,0.548372 /)
XCGA_LKT(89,14,1:6)=(/ 0.917363,0.877997,0.841440,0.838387,0.838030,0.936293 /)
XEXT_COEFF_550_LKT(89,14)=8.196600 !rg=11.9469 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,15,1:6)=(/ 6.811000,6.826300,6.848400,6.905400,6.975400,7.044800 /)
XPIZA_LKT(89,15,1:6)=(/ 0.607527,0.745433,0.985170,0.987462,0.960230,0.549661 /)
XCGA_LKT(89,15,1:6)=(/ 0.000000,0.881180,0.842533,0.839823,0.839820,0.936517 /)
XEXT_COEFF_550_LKT(89,15)=6.839000 !rg=11.9469 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,16,1:6)=(/ 5.658400,5.672600,5.689700,5.719500,5.774200,5.834800 /)
XPIZA_LKT(89,16,1:6)=(/ 0.598053,0.727180,0.982908,0.985792,0.956627,0.550799 /)
XCGA_LKT(89,16,1:6)=(/ 0.000000,0.884753,0.843217,0.840577,0.841570,0.936700 /)
XEXT_COEFF_550_LKT(89,16)=5.674000 !rg=11.9469 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,17,1:6)=(/ 4.708400,4.718200,4.733800,4.757000,4.792300,4.841000 /)
XPIZA_LKT(89,17,1:6)=(/ 0.589806,0.709273,0.980299,0.983867,0.952172,0.551762 /)
XCGA_LKT(89,17,1:6)=(/ 0.000000,0.000000,0.844180,0.841673,0.843487,0.936840 /)
XEXT_COEFF_550_LKT(89,17)=4.720900 !rg=11.9469 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,18,1:6)=(/ 3.918700,3.924900,3.935800,3.956500,3.987400,4.017600 /)
XPIZA_LKT(89,18,1:6)=(/ 0.582649,0.691930,0.977131,0.981501,0.946857,0.552588 /)
XCGA_LKT(89,18,1:6)=(/ 0.000000,0.000000,0.845053,0.842927,0.845270,0.936950 /)
XEXT_COEFF_550_LKT(89,18)=3.928400 !rg=11.9469 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,19,1:6)=(/ 0.000000,3.261100,3.267500,3.280700,3.304000,3.329800 /)
XPIZA_LKT(89,19,1:6)=(/ 0.900,0.675268,0.973400,0.978654,0.940638,0.553302 /)
XCGA_LKT(89,19,1:6)=(/ 0.900,0.000000,0.000000,0.843613,0.846867,0.937030 /)
XEXT_COEFF_550_LKT(89,19)=3.262000 !rg=11.9469 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,20,1:6)=(/ 0.000000,2.710600,2.716300,2.725200,2.743800,2.761900 /)
XPIZA_LKT(89,20,1:6)=(/ 0.900,0.659573,0.968969,0.975349,0.933300,0.553900 /)
XCGA_LKT(89,20,1:6)=(/ 0.900,0.000000,0.000000,0.844600,0.848700,0.937087 /)
XEXT_COEFF_550_LKT(89,20)=2.711000 !rg=11.9469 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,1,1:6)=(/ 46.643000,47.866000,47.784000,48.401000,50.849000,51.418000 /)
XPIZA_LKT(90,1,1:6)=(/ 0.745894,0.888380,0.994650,0.990208,0.977180,0.526340 /)
XCGA_LKT(90,1,1:6)=(/ 0.885337,0.856127,0.830620,0.818910,0.813113,0.931403 /)
XEXT_COEFF_550_LKT(90,1)=47.660000 !rg=12.9427 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,2,1:6)=(/ 45.000000,45.490000,45.667000,47.154000,47.379000,49.165000 /)
XPIZA_LKT(90,2,1:6)=(/ 0.743811,0.885801,0.994573,0.991831,0.978479,0.527102 /)
XCGA_LKT(90,2,1:6)=(/ 0.886233,0.856327,0.831243,0.817860,0.813397,0.931610 /)
XEXT_COEFF_550_LKT(90,2)=45.736000 !rg=12.9427 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,3,1:6)=(/ 42.016000,42.042000,42.437000,43.266000,44.412000,45.432000 /)
XPIZA_LKT(90,3,1:6)=(/ 0.740030,0.882354,0.994188,0.992439,0.976021,0.528405 /)
XCGA_LKT(90,3,1:6)=(/ 0.888363,0.857037,0.829550,0.821867,0.813730,0.931933 /)
XEXT_COEFF_550_LKT(90,3)=42.386000 !rg=12.9427 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,4,1:6)=(/ 37.865000,37.971000,38.449000,39.001000,39.968000,40.918000 /)
XPIZA_LKT(90,4,1:6)=(/ 0.731704,0.877731,0.994129,0.992723,0.974531,0.530070 /)
XCGA_LKT(90,4,1:6)=(/ 0.889493,0.857923,0.830027,0.823317,0.816950,0.932330 /)
XEXT_COEFF_550_LKT(90,4)=38.277000 !rg=12.9427 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,5,1:6)=(/ 33.493000,33.762000,33.977000,34.284000,35.255000,36.102000 /)
XPIZA_LKT(90,5,1:6)=(/ 0.722027,0.872170,0.993784,0.992759,0.972416,0.531991 /)
XCGA_LKT(90,5,1:6)=(/ 0.890823,0.859577,0.831653,0.826030,0.816813,0.932787 /)
XEXT_COEFF_550_LKT(90,5)=33.796000 !rg=12.9427 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,6,1:6)=(/ 29.279000,29.367000,29.703000,29.979000,30.600000,31.387000 /)
XPIZA_LKT(90,6,1:6)=(/ 0.711614,0.864441,0.993674,0.992985,0.971942,0.534065 /)
XCGA_LKT(90,6,1:6)=(/ 0.893173,0.860540,0.833073,0.826837,0.821133,0.933277 /)
XEXT_COEFF_550_LKT(90,6)=29.649000 !rg=12.9427 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,7,1:6)=(/ 25.224000,25.317000,25.629000,25.805000,26.333000,26.916000 /)
XPIZA_LKT(90,7,1:6)=(/ 0.699811,0.855842,0.993307,0.992754,0.971442,0.536262 /)
XCGA_LKT(90,7,1:6)=(/ 0.895620,0.862173,0.834487,0.828030,0.823743,0.933790 /)
XEXT_COEFF_550_LKT(90,7)=25.472000 !rg=12.9427 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,8,1:6)=(/ 21.551000,21.621000,21.792000,22.079000,22.357000,22.901000 /)
XPIZA_LKT(90,8,1:6)=(/ 0.686967,0.845662,0.992787,0.992532,0.970941,0.538438 /)
XCGA_LKT(90,8,1:6)=(/ 0.898250,0.864107,0.835510,0.829987,0.825780,0.934283 /)
XEXT_COEFF_550_LKT(90,8)=21.755000 !rg=12.9427 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,9,1:6)=(/ 18.275000,18.368000,18.490000,18.618000,18.972000,19.328000 /)
XPIZA_LKT(90,9,1:6)=(/ 0.673656,0.834052,0.992127,0.992213,0.969998,0.540574 /)
XCGA_LKT(90,9,1:6)=(/ 0.901280,0.866393,0.836437,0.832097,0.827420,0.934760 /)
XEXT_COEFF_550_LKT(90,9)=18.386000 !rg=12.9427 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,10,1:6)=(/ 15.434000,15.494000,15.596000,15.741000,16.033000,16.253000 /)
XPIZA_LKT(90,10,1:6)=(/ 0.660726,0.820649,0.991390,0.991771,0.969015,0.542571 /)
XCGA_LKT(90,10,1:6)=(/ 0.904470,0.868490,0.837910,0.833330,0.828610,0.935187 /)
XEXT_COEFF_550_LKT(90,10)=15.513000 !rg=12.9427 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,11,1:6)=(/ 12.954000,12.996000,13.075000,13.148000,13.432000,13.585000 /)
XPIZA_LKT(90,11,1:6)=(/ 0.647606,0.805932,0.990498,0.991236,0.968060,0.544458 /)
XCGA_LKT(90,11,1:6)=(/ 0.907933,0.870800,0.839147,0.834957,0.832747,0.935577 /)
XEXT_COEFF_550_LKT(90,11)=13.055000 !rg=12.9427 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,12,1:6)=(/ 10.838000,10.879000,10.928000,11.016000,11.141000,11.317000 /)
XPIZA_LKT(90,12,1:6)=(/ 0.635096,0.789979,0.989346,0.990390,0.966589,0.546183 /)
XCGA_LKT(90,12,1:6)=(/ 0.911600,0.873600,0.839960,0.836883,0.834303,0.935917 /)
XEXT_COEFF_550_LKT(90,12)=10.892000 !rg=12.9427 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,13,1:6)=(/ 9.041900,9.068400,9.099000,9.165200,9.275600,9.404700 /)
XPIZA_LKT(90,13,1:6)=(/ 0.623200,0.772741,0.987920,0.989370,0.964672,0.547735 /)
XCGA_LKT(90,13,1:6)=(/ 0.915450,0.876420,0.840900,0.836967,0.836977,0.936207 /)
XEXT_COEFF_550_LKT(90,13)=9.079800 !rg=12.9427 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,14,1:6)=(/ 7.544100,7.558900,7.590500,7.657900,7.717700,7.815300 /)
XPIZA_LKT(90,14,1:6)=(/ 0.612621,0.755122,0.986272,0.988269,0.962029,0.549093 /)
XCGA_LKT(90,14,1:6)=(/ 0.919277,0.879277,0.842043,0.839430,0.838460,0.936443 /)
XEXT_COEFF_550_LKT(90,14)=7.572000 !rg=12.9427 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,15,1:6)=(/ 6.281900,6.298900,6.318200,6.353700,6.418700,6.488800 /)
XPIZA_LKT(90,15,1:6)=(/ 0.602697,0.737136,0.984266,0.986780,0.958731,0.550283 /)
XCGA_LKT(90,15,1:6)=(/ 0.000000,0.882857,0.842790,0.840033,0.840343,0.936637 /)
XEXT_COEFF_550_LKT(90,15)=6.302600 !rg=12.9427 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,16,1:6)=(/ 5.220600,5.231700,5.250700,5.269900,5.335900,5.375400 /)
XPIZA_LKT(90,16,1:6)=(/ 0.593822,0.718778,0.981841,0.985010,0.954865,0.551332 /)
XCGA_LKT(90,16,1:6)=(/ 0.000000,0.886387,0.843800,0.841317,0.842937,0.936793 /)
XEXT_COEFF_550_LKT(90,16)=5.239100 !rg=12.9427 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,17,1:6)=(/ 4.345200,4.350900,4.363600,4.389200,4.419500,4.460600 /)
XPIZA_LKT(90,17,1:6)=(/ 0.586124,0.701043,0.978967,0.982847,0.950017,0.552218 /)
XCGA_LKT(90,17,1:6)=(/ 0.000000,0.000000,0.844543,0.842433,0.844460,0.936913 /)
XEXT_COEFF_550_LKT(90,17)=4.356300 !rg=12.9427 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,18,1:6)=(/ 3.615500,3.622100,3.631400,3.645600,3.670000,3.702600 /)
XPIZA_LKT(90,18,1:6)=(/ 0.579312,0.684176,0.975598,0.980310,0.944207,0.552975 /)
XCGA_LKT(90,18,1:6)=(/ 0.000000,0.000000,0.845307,0.843290,0.845677,0.937003 /)
XEXT_COEFF_550_LKT(90,18)=3.624100 !rg=12.9427 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,19,1:6)=(/ 0.000000,3.008700,3.015200,3.024900,3.050700,3.069200 /)
XPIZA_LKT(90,19,1:6)=(/ 0.900,0.667844,0.971557,0.977278,0.937647,0.553627 /)
XCGA_LKT(90,19,1:6)=(/ 0.900,0.000000,0.000000,0.844167,0.847953,0.937067 /)
XEXT_COEFF_550_LKT(90,19)=3.011000 !rg=12.9427 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,20,1:6)=(/ 0.000000,2.500400,2.505000,2.515000,2.527200,2.546100 /)
XPIZA_LKT(90,20,1:6)=(/ 0.900,0.652527,0.966800,0.973696,0.929603,0.554169 /)
XCGA_LKT(90,20,1:6)=(/ 0.900,0.000000,0.000000,0.845157,0.849413,0.937110 /)
XEXT_COEFF_550_LKT(90,20)=2.502400 !rg=12.9427 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,1,1:6)=(/ 43.530000,43.751000,44.095000,43.688000,45.516000,47.211000 /)
XPIZA_LKT(91,1,1:6)=(/ 0.741748,0.885092,0.994426,0.991913,0.977295,0.528573 /)
XCGA_LKT(91,1,1:6)=(/ 0.887127,0.857527,0.828520,0.827320,0.809117,0.932090 /)
XEXT_COEFF_550_LKT(91,1)=43.994000 !rg=14.0216 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,2,1:6)=(/ 41.760000,41.879000,42.227000,42.842000,43.754000,45.145000 /)
XPIZA_LKT(91,2,1:6)=(/ 0.738994,0.882406,0.994436,0.992526,0.976573,0.529276 /)
XCGA_LKT(91,2,1:6)=(/ 0.888437,0.857240,0.831213,0.822953,0.815903,0.932263 /)
XEXT_COEFF_550_LKT(91,2)=42.155000 !rg=14.0216 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,3,1:6)=(/ 38.548000,38.999000,39.145000,39.510000,40.953000,41.723000 /)
XPIZA_LKT(91,3,1:6)=(/ 0.732462,0.879492,0.994371,0.992217,0.974160,0.530496 /)
XCGA_LKT(91,3,1:6)=(/ 0.888813,0.858353,0.831717,0.823783,0.813233,0.932557 /)
XEXT_COEFF_550_LKT(91,3)=39.007000 !rg=14.0216 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,4,1:6)=(/ 34.873000,35.230000,35.493000,35.575000,36.938000,37.585000 /)
XPIZA_LKT(91,4,1:6)=(/ 0.724332,0.874619,0.993988,0.992940,0.972440,0.532048 /)
XCGA_LKT(91,4,1:6)=(/ 0.890637,0.859623,0.832230,0.826420,0.815983,0.932917 /)
XEXT_COEFF_550_LKT(91,4)=35.177000 !rg=14.0216 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,5,1:6)=(/ 30.951000,31.005000,31.372000,31.831000,32.551000,33.170000 /)
XPIZA_LKT(91,5,1:6)=(/ 0.715376,0.867670,0.993716,0.992838,0.971462,0.533844 /)
XCGA_LKT(91,5,1:6)=(/ 0.892677,0.860110,0.831970,0.825720,0.820553,0.933330 /)
XEXT_COEFF_550_LKT(91,5)=31.261000 !rg=14.0216 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,6,1:6)=(/ 26.978000,27.177000,27.309000,27.813000,28.078000,28.844000 /)
XPIZA_LKT(91,6,1:6)=(/ 0.704038,0.860323,0.993392,0.992702,0.971181,0.535780 /)
XCGA_LKT(91,6,1:6)=(/ 0.894587,0.862083,0.833457,0.827530,0.823493,0.933770 /)
XEXT_COEFF_550_LKT(91,6)=27.174000 !rg=14.0216 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,7,1:6)=(/ 23.236000,23.406000,23.472000,23.772000,24.288000,24.743000 /)
XPIZA_LKT(91,7,1:6)=(/ 0.691965,0.851054,0.993047,0.992629,0.970658,0.537833 /)
XCGA_LKT(91,7,1:6)=(/ 0.897043,0.863843,0.835310,0.829543,0.824693,0.934230 /)
XEXT_COEFF_550_LKT(91,7)=23.391000 !rg=14.0216 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,8,1:6)=(/ 19.872000,19.978000,20.094000,20.324000,20.573000,21.058000 /)
XPIZA_LKT(91,8,1:6)=(/ 0.679476,0.840157,0.992552,0.992541,0.970644,0.539864 /)
XCGA_LKT(91,8,1:6)=(/ 0.899877,0.865400,0.836367,0.831910,0.827803,0.934673 /)
XEXT_COEFF_550_LKT(91,8)=19.996000 !rg=14.0216 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,9,1:6)=(/ 16.862000,16.903000,17.045000,17.189000,17.503000,17.778000 /)
XPIZA_LKT(91,9,1:6)=(/ 0.666764,0.827349,0.991842,0.992074,0.969527,0.541857 /)
XCGA_LKT(91,9,1:6)=(/ 0.903053,0.867290,0.837507,0.832707,0.829660,0.935093 /)
XEXT_COEFF_550_LKT(91,9)=17.004000 !rg=14.0216 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,10,1:6)=(/ 14.233000,14.289000,14.377000,14.475000,14.674000,14.953000 /)
XPIZA_LKT(91,10,1:6)=(/ 0.653809,0.813733,0.990952,0.991626,0.969007,0.543716 /)
XCGA_LKT(91,10,1:6)=(/ 0.906187,0.869637,0.838233,0.834147,0.832007,0.935473 /)
XEXT_COEFF_550_LKT(91,10)=14.305000 !rg=14.0216 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,11,1:6)=(/ 11.953000,11.984000,12.048000,12.176000,12.300000,12.502000 /)
XPIZA_LKT(91,11,1:6)=(/ 0.641259,0.798485,0.990013,0.990845,0.967378,0.545470 /)
XCGA_LKT(91,11,1:6)=(/ 0.909763,0.872033,0.839790,0.836187,0.834287,0.935817 /)
XEXT_COEFF_550_LKT(91,11)=12.024000 !rg=14.0216 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,12,1:6)=(/ 9.997500,10.029000,10.076000,10.145000,10.293000,10.418000 /)
XPIZA_LKT(91,12,1:6)=(/ 0.629029,0.782034,0.988761,0.990047,0.965823,0.547070 /)
XCGA_LKT(91,12,1:6)=(/ 0.913467,0.874903,0.840610,0.837090,0.836277,0.936113 /)
XEXT_COEFF_550_LKT(91,12)=10.037000 !rg=14.0216 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,13,1:6)=(/ 8.344000,8.363400,8.407700,8.459300,8.539100,8.659400 /)
XPIZA_LKT(91,13,1:6)=(/ 0.617821,0.764537,0.987237,0.988899,0.963455,0.548508 /)
XCGA_LKT(91,13,1:6)=(/ 0.917303,0.877713,0.841587,0.838297,0.837350,0.936367 /)
XEXT_COEFF_550_LKT(91,13)=8.380700 !rg=14.0216 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,14,1:6)=(/ 6.958000,6.976800,7.000700,7.040600,7.098600,7.197600 /)
XPIZA_LKT(91,14,1:6)=(/ 0.607531,0.746817,0.985453,0.987660,0.960794,0.549763 /)
XCGA_LKT(91,14,1:6)=(/ 0.000000,0.881070,0.842307,0.839667,0.840133,0.936570 /)
XEXT_COEFF_550_LKT(91,14)=6.985900 !rg=14.0216 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,15,1:6)=(/ 5.797000,5.809200,5.835200,5.851000,5.918100,5.977200 /)
XPIZA_LKT(91,15,1:6)=(/ 0.598292,0.728715,0.983281,0.985981,0.957288,0.550860 /)
XCGA_LKT(91,15,1:6)=(/ 0.000000,0.884373,0.843390,0.840383,0.841983,0.936737 /)
XEXT_COEFF_550_LKT(91,15)=5.818200 !rg=14.0216 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,16,1:6)=(/ 4.817600,4.825600,4.840400,4.873300,4.901900,4.952500 /)
XPIZA_LKT(91,16,1:6)=(/ 0.589824,0.710514,0.980635,0.984149,0.952763,0.551825 /)
XCGA_LKT(91,16,1:6)=(/ 0.000000,0.000000,0.844183,0.842170,0.843547,0.936873 /)
XEXT_COEFF_550_LKT(91,16)=4.830200 !rg=14.0216 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,17,1:6)=(/ 4.009200,4.016300,4.024400,4.042400,4.067100,4.110500 /)
XPIZA_LKT(91,17,1:6)=(/ 0.582607,0.693162,0.977527,0.981783,0.947558,0.552638 /)
XCGA_LKT(91,17,1:6)=(/ 0.000000,0.000000,0.844807,0.842573,0.845103,0.936973 /)
XEXT_COEFF_550_LKT(91,17)=4.020200 !rg=14.0216 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,18,1:6)=(/ 0.000000,3.341700,3.350900,3.361100,3.389800,3.412600 /)
XPIZA_LKT(91,18,1:6)=(/ 0.900,0.676511,0.973913,0.979024,0.941650,0.553330 /)
XCGA_LKT(91,18,1:6)=(/ 0.900,0.000000,0.845847,0.843590,0.846980,0.937047 /)
XEXT_COEFF_550_LKT(91,18)=3.344600 !rg=14.0216 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,19,1:6)=(/ 0.000000,2.776100,2.781300,2.795100,2.807700,2.829200 /)
XPIZA_LKT(91,19,1:6)=(/ 0.900,0.660620,0.969584,0.975818,0.934172,0.553923 /)
XCGA_LKT(91,19,1:6)=(/ 0.900,0.000000,0.000000,0.844897,0.848590,0.937097 /)
XEXT_COEFF_550_LKT(91,19)=2.777200 !rg=14.0216 sigma=2.85 wvl=0.55
 
END SUBROUTINE SALT_OPT_LKT_SET9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SALT_OPT_LKT_SET10()

  USE MODD_SALT_OPT_LKT
  
  IMPLICIT NONE 
 
XEXT_COEFF_WVL_LKT(91,20,1:6)=(/ 0.000000,0.000000,2.311000,2.317600,2.329300,2.347300 /)
XPIZA_LKT(91,20,1:6)=(/ 0.900,0.900,0.964503,0.971907,0.925774,0.554415 /)
XCGA_LKT(91,20,1:6)=(/ 0.900,0.900,0.000000,0.845430,0.850140,0.937130 /)
XEXT_COEFF_550_LKT(91,20)=2.309200 !rg=14.0216 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,1,1:6)=(/ 40.518000,39.817000,40.565000,40.265000,41.646000,43.355000 /)
XPIZA_LKT(92,1,1:6)=(/ 0.736740,0.879959,0.994549,0.992698,0.976214,0.530647 /)
XCGA_LKT(92,1,1:6)=(/ 0.889670,0.855327,0.832003,0.826907,0.820620,0.932693 /)
XEXT_COEFF_550_LKT(92,1)=40.939000 !rg=15.1904 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,2,1:6)=(/ 38.450000,38.846000,38.972000,39.604000,40.457000,41.464000 /)
XPIZA_LKT(92,2,1:6)=(/ 0.731370,0.879371,0.994422,0.992541,0.974694,0.531323 /)
XCGA_LKT(92,2,1:6)=(/ 0.889557,0.858903,0.832057,0.823273,0.818933,0.932850 /)
XEXT_COEFF_550_LKT(92,2)=38.667000 !rg=15.1904 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,3,1:6)=(/ 35.584000,35.720000,36.152000,36.578000,37.299000,38.327000 /)
XPIZA_LKT(92,3,1:6)=(/ 0.725364,0.874896,0.994132,0.992556,0.972356,0.532465 /)
XCGA_LKT(92,3,1:6)=(/ 0.890413,0.858650,0.833377,0.825567,0.817247,0.933113 /)
XEXT_COEFF_550_LKT(92,3)=36.017000 !rg=15.1904 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,4,1:6)=(/ 32.148000,32.278000,32.661000,32.991000,33.810000,34.532000 /)
XPIZA_LKT(92,4,1:6)=(/ 0.717318,0.869808,0.993987,0.992673,0.970641,0.533914 /)
XCGA_LKT(92,4,1:6)=(/ 0.891903,0.859623,0.833757,0.826137,0.818517,0.933440 /)
XEXT_COEFF_550_LKT(92,4)=32.505000 !rg=15.1904 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,5,1:6)=(/ 28.478000,28.732000,28.792000,29.090000,30.028000,30.482000 /)
XPIZA_LKT(92,5,1:6)=(/ 0.707498,0.863804,0.993715,0.992895,0.970694,0.535587 /)
XCGA_LKT(92,5,1:6)=(/ 0.893843,0.861670,0.833910,0.828163,0.820357,0.933810 /)
XEXT_COEFF_550_LKT(92,5)=28.682000 !rg=15.1904 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,6,1:6)=(/ 24.865000,24.995000,25.290000,25.482000,26.014000,26.514000 /)
XPIZA_LKT(92,6,1:6)=(/ 0.696514,0.855028,0.993227,0.992600,0.970340,0.537394 /)
XCGA_LKT(92,6,1:6)=(/ 0.896037,0.862780,0.834990,0.828097,0.823767,0.934207 /)
XEXT_COEFF_550_LKT(92,6)=25.095000 !rg=15.1904 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,7,1:6)=(/ 21.435000,21.515000,21.685000,21.841000,22.426000,22.749000 /)
XPIZA_LKT(92,7,1:6)=(/ 0.684620,0.845211,0.992759,0.992587,0.969697,0.539309 /)
XCGA_LKT(92,7,1:6)=(/ 0.898680,0.864423,0.836127,0.830940,0.826040,0.934617 /)
XEXT_COEFF_550_LKT(92,7)=21.666000 !rg=15.1904 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,8,1:6)=(/ 18.319000,18.399000,18.472000,18.757000,19.001000,19.367000 /)
XPIZA_LKT(92,8,1:6)=(/ 0.672181,0.833969,0.992263,0.991995,0.969959,0.541202 /)
XCGA_LKT(92,8,1:6)=(/ 0.901447,0.866490,0.836797,0.831053,0.828807,0.935013 /)
XEXT_COEFF_550_LKT(92,8)=18.385000 !rg=15.1904 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,9,1:6)=(/ 15.542000,15.627000,15.675000,15.824000,16.133000,16.355000 /)
XPIZA_LKT(92,9,1:6)=(/ 0.659564,0.821243,0.991484,0.991798,0.969116,0.543057 /)
XCGA_LKT(92,9,1:6)=(/ 0.904677,0.868737,0.837917,0.833513,0.829253,0.935390 /)
XEXT_COEFF_550_LKT(92,9)=15.628000 !rg=15.1904 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,10,1:6)=(/ 13.133000,13.170000,13.260000,13.393000,13.548000,13.760000 /)
XPIZA_LKT(92,10,1:6)=(/ 0.647172,0.806619,0.990631,0.991356,0.968034,0.544786 /)
XCGA_LKT(92,10,1:6)=(/ 0.908067,0.870793,0.839713,0.835550,0.832627,0.935723 /)
XEXT_COEFF_550_LKT(92,10)=13.225000 !rg=15.1904 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,11,1:6)=(/ 11.023000,11.063000,11.101000,11.164000,11.353000,11.507000 /)
XPIZA_LKT(92,11,1:6)=(/ 0.634723,0.790875,0.989448,0.990582,0.966794,0.546413 /)
XCGA_LKT(92,11,1:6)=(/ 0.911630,0.873517,0.839870,0.836677,0.835237,0.936027 /)
XEXT_COEFF_550_LKT(92,11)=11.075000 !rg=15.1904 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,12,1:6)=(/ 9.224900,9.245900,9.304500,9.389200,9.476700,9.591000 /)
XPIZA_LKT(92,12,1:6)=(/ 0.623217,0.773861,0.988154,0.989606,0.964973,0.547896 /)
XCGA_LKT(92,12,1:6)=(/ 0.915390,0.876183,0.841357,0.838280,0.836687,0.936287 /)
XEXT_COEFF_550_LKT(92,12)=9.278500 !rg=15.1904 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,13,1:6)=(/ 7.696000,7.722400,7.745600,7.792700,7.867500,7.974300 /)
XPIZA_LKT(92,13,1:6)=(/ 0.612388,0.756509,0.986450,0.988385,0.962432,0.549226 /)
XCGA_LKT(92,13,1:6)=(/ 0.919127,0.879380,0.841780,0.838923,0.839397,0.936507 /)
XEXT_COEFF_550_LKT(92,13)=7.729700 !rg=15.1904 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,14,1:6)=(/ 6.419600,6.435200,6.461700,6.495300,6.558800,6.629500 /)
XPIZA_LKT(92,14,1:6)=(/ 0.602680,0.738413,0.984533,0.986960,0.959289,0.550383 /)
XCGA_LKT(92,14,1:6)=(/ 0.000000,0.882577,0.842987,0.840073,0.841440,0.936683 /)
XEXT_COEFF_550_LKT(92,14)=6.446200 !rg=15.1904 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,15,1:6)=(/ 5.348400,5.359600,5.372600,5.411000,5.455800,5.506600 /)
XPIZA_LKT(92,15,1:6)=(/ 0.593886,0.720371,0.982147,0.985255,0.955258,0.551393 /)
XCGA_LKT(92,15,1:6)=(/ 0.000000,0.886057,0.843653,0.841237,0.842020,0.936827 /)
XEXT_COEFF_550_LKT(92,15)=5.365300 !rg=15.1904 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,16,1:6)=(/ 4.445000,4.452400,4.463100,4.483700,4.521000,4.563400 /)
XPIZA_LKT(92,16,1:6)=(/ 0.586058,0.702349,0.979323,0.983171,0.950739,0.552279 /)
XCGA_LKT(92,16,1:6)=(/ 0.000000,0.000000,0.844373,0.842370,0.844363,0.936940 /)
XEXT_COEFF_550_LKT(92,16)=4.454900 !rg=15.1904 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,17,1:6)=(/ 3.699500,3.706000,3.714900,3.726500,3.758100,3.788200 /)
XPIZA_LKT(92,17,1:6)=(/ 0.579288,0.685324,0.976035,0.980629,0.945064,0.553024 /)
XCGA_LKT(92,17,1:6)=(/ 0.000000,0.000000,0.845230,0.843130,0.846187,0.937023 /)
XEXT_COEFF_550_LKT(92,17)=3.710700 !rg=15.1904 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,18,1:6)=(/ 0.000000,3.083400,3.089100,3.103200,3.123300,3.145400 /)
XPIZA_LKT(92,18,1:6)=(/ 0.900,0.669045,0.972094,0.977706,0.938448,0.553655 /)
XCGA_LKT(92,18,1:6)=(/ 0.900,0.000000,0.000000,0.844173,0.847300,0.937083 /)
XEXT_COEFF_550_LKT(92,18)=3.085200 !rg=15.1904 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,19,1:6)=(/ 0.000000,2.561700,2.565100,2.574800,2.590500,2.608100 /)
XPIZA_LKT(92,19,1:6)=(/ 0.900,0.653625,0.967449,0.974168,0.930746,0.554193 /)
XCGA_LKT(92,19,1:6)=(/ 0.900,0.000000,0.000000,0.845083,0.849303,0.937123 /)
XEXT_COEFF_550_LKT(92,19)=2.562700 !rg=15.1904 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,20,1:6)=(/ 0.000000,0.000000,2.132900,2.138400,2.151500,2.164100 /)
XPIZA_LKT(92,20,1:6)=(/ 0.900,0.900,0.962057,0.970022,0.921853,0.554637 /)
XCGA_LKT(92,20,1:6)=(/ 0.900,0.900,0.000000,0.845900,0.851120,0.937143 /)
XEXT_COEFF_550_LKT(92,20)=2.130700 !rg=15.1904 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,1,1:6)=(/ 37.154000,37.990000,37.754000,37.761000,38.429000,39.827000 /)
XPIZA_LKT(93,1,1:6)=(/ 0.728518,0.878660,0.994112,0.993338,0.973181,0.532628 /)
XCGA_LKT(93,1,1:6)=(/ 0.891187,0.860690,0.832210,0.828720,0.820147,0.933233 /)
XEXT_COEFF_550_LKT(93,1)=37.434000 !rg=16.4566 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,2,1:6)=(/ 35.397000,35.550000,35.940000,36.373000,37.457000,38.092000 /)
XPIZA_LKT(93,2,1:6)=(/ 0.723477,0.874814,0.994222,0.993093,0.972659,0.533248 /)
XCGA_LKT(93,2,1:6)=(/ 0.891000,0.859043,0.832627,0.824943,0.817600,0.933373 /)
XEXT_COEFF_550_LKT(93,2)=35.574000 !rg=16.4566 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,3,1:6)=(/ 32.830000,33.149000,33.333000,33.628000,34.662000,35.215000 /)
XPIZA_LKT(93,3,1:6)=(/ 0.718071,0.871632,0.993855,0.992992,0.969491,0.534314 /)
XCGA_LKT(93,3,1:6)=(/ 0.892007,0.860787,0.833090,0.827593,0.819207,0.933613 /)
XEXT_COEFF_550_LKT(93,3)=33.112000 !rg=16.4566 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,4,1:6)=(/ 29.658000,29.926000,30.048000,30.520000,31.212000,31.734000 /)
XPIZA_LKT(93,4,1:6)=(/ 0.709530,0.865964,0.993825,0.992969,0.969413,0.535665 /)
XCGA_LKT(93,4,1:6)=(/ 0.893523,0.861537,0.833857,0.827943,0.822497,0.933907 /)
XEXT_COEFF_550_LKT(93,4)=29.894000 !rg=16.4566 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,5,1:6)=(/ 26.261000,26.364000,26.705000,26.906000,27.482000,28.018000 /)
XPIZA_LKT(93,5,1:6)=(/ 0.699882,0.858050,0.993333,0.992689,0.969389,0.537226 /)
XCGA_LKT(93,5,1:6)=(/ 0.895363,0.862067,0.835303,0.828957,0.822977,0.934240 /)
XEXT_COEFF_550_LKT(93,5)=26.528000 !rg=16.4566 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,6,1:6)=(/ 22.959000,23.073000,23.272000,23.485000,23.715000,24.376000 /)
XPIZA_LKT(93,6,1:6)=(/ 0.689377,0.849828,0.993105,0.992738,0.970394,0.538909 /)
XCGA_LKT(93,6,1:6)=(/ 0.897783,0.863967,0.835667,0.830157,0.826650,0.934593 /)
XEXT_COEFF_550_LKT(93,6)=23.155000 !rg=16.4566 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,7,1:6)=(/ 19.779000,19.830000,19.986000,20.281000,20.518000,20.921000 /)
XPIZA_LKT(93,7,1:6)=(/ 0.677456,0.839150,0.992564,0.992412,0.969615,0.540693 /)
XCGA_LKT(93,7,1:6)=(/ 0.900360,0.865507,0.836847,0.832107,0.828987,0.934960 /)
XEXT_COEFF_550_LKT(93,7)=19.908000 !rg=16.4566 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,8,1:6)=(/ 16.918000,16.962000,17.098000,17.274000,17.488000,17.815000 /)
XPIZA_LKT(93,8,1:6)=(/ 0.665258,0.827363,0.991857,0.992155,0.969493,0.542455 /)
XCGA_LKT(93,8,1:6)=(/ 0.903413,0.867327,0.837993,0.833257,0.829263,0.935313 /)
XEXT_COEFF_550_LKT(93,8)=17.030000 !rg=16.4566 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,9,1:6)=(/ 14.337000,14.387000,14.505000,14.576000,14.791000,15.048000 /)
XPIZA_LKT(93,9,1:6)=(/ 0.652576,0.813870,0.991126,0.991513,0.968685,0.544179 /)
XCGA_LKT(93,9,1:6)=(/ 0.906453,0.869690,0.839070,0.834217,0.831743,0.935647 /)
XEXT_COEFF_550_LKT(93,9)=14.433000 !rg=16.4566 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,10,1:6)=(/ 12.112000,12.175000,12.220000,12.303000,12.457000,12.664000 /)
XPIZA_LKT(93,10,1:6)=(/ 0.640473,0.799425,0.990054,0.990899,0.967376,0.545784 /)
XCGA_LKT(93,10,1:6)=(/ 0.909873,0.872497,0.839630,0.835990,0.833563,0.935943 /)
XEXT_COEFF_550_LKT(93,10)=12.176000 !rg=16.4566 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,11,1:6)=(/ 10.172000,10.204000,10.266000,10.330000,10.459000,10.593000 /)
XPIZA_LKT(93,11,1:6)=(/ 0.628683,0.783043,0.988903,0.990132,0.966098,0.547292 /)
XCGA_LKT(93,11,1:6)=(/ 0.913537,0.874773,0.840863,0.837247,0.837020,0.936210 /)
XEXT_COEFF_550_LKT(93,11)=10.224000 !rg=16.4566 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,12,1:6)=(/ 8.509500,8.537800,8.564100,8.623000,8.698400,8.831200 /)
XPIZA_LKT(93,12,1:6)=(/ 0.617551,0.765972,0.987452,0.989040,0.963925,0.548664 /)
XCGA_LKT(93,12,1:6)=(/ 0.917217,0.877767,0.841513,0.838480,0.838203,0.936437 /)
XEXT_COEFF_550_LKT(93,12)=8.547200 !rg=16.4566 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,13,1:6)=(/ 7.101700,7.119300,7.151900,7.182800,7.270700,7.344200 /)
XPIZA_LKT(93,13,1:6)=(/ 0.607373,0.748050,0.985680,0.987805,0.961134,0.549890 /)
XCGA_LKT(93,13,1:6)=(/ 0.000000,0.880833,0.842590,0.839633,0.840893,0.936627 /)
XEXT_COEFF_550_LKT(93,13)=7.134000 !rg=16.4566 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,14,1:6)=(/ 5.923100,5.938600,5.949500,5.996200,6.036200,6.107000 /)
XPIZA_LKT(93,14,1:6)=(/ 0.598083,0.730091,0.983506,0.986289,0.957765,0.550957 /)
XCGA_LKT(93,14,1:6)=(/ 0.000000,0.884223,0.843240,0.841190,0.841847,0.936780 /)
XEXT_COEFF_550_LKT(93,14)=5.943400 !rg=16.4566 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,15,1:6)=(/ 4.934800,4.943800,4.957500,4.984800,5.020600,5.073500 /)
XPIZA_LKT(93,15,1:6)=(/ 0.589813,0.712010,0.981008,0.984411,0.953474,0.551884 /)
XCGA_LKT(93,15,1:6)=(/ 0.000000,0.000000,0.844137,0.841877,0.843410,0.936900 /)
XEXT_COEFF_550_LKT(93,15)=4.945800 !rg=16.4566 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,16,1:6)=(/ 4.101700,4.108900,4.121900,4.140100,4.169200,4.205300 /)
XPIZA_LKT(93,16,1:6)=(/ 0.582529,0.694368,0.977969,0.982137,0.948352,0.552696 /)
XCGA_LKT(93,16,1:6)=(/ 0.000000,0.000000,0.844977,0.842833,0.845617,0.936997 /)
XEXT_COEFF_550_LKT(93,16)=4.112400 !rg=16.4566 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,17,1:6)=(/ 0.000000,3.419600,3.425400,3.440900,3.459900,3.491400 /)
XPIZA_LKT(93,17,1:6)=(/ 0.900,0.677605,0.974375,0.979427,0.942324,0.553376 /)
XCGA_LKT(93,17,1:6)=(/ 0.900,0.000000,0.845630,0.843687,0.846863,0.937063 /)
XEXT_COEFF_550_LKT(93,17)=3.420300 !rg=16.4566 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,18,1:6)=(/ 0.000000,2.845100,2.849100,2.862200,2.877600,2.899500 /)
XPIZA_LKT(93,18,1:6)=(/ 0.900,0.661742,0.970156,0.976248,0.935274,0.553950 /)
XCGA_LKT(93,18,1:6)=(/ 0.900,0.000000,0.000000,0.844500,0.848367,0.937113 /)
XEXT_COEFF_550_LKT(93,18)=2.845900 !rg=16.4566 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,19,1:6)=(/ 0.000000,2.363800,2.368500,2.375200,2.388700,2.404500 /)
XPIZA_LKT(93,19,1:6)=(/ 0.900,0.646895,0.965199,0.972479,0.926950,0.554438 /)
XCGA_LKT(93,19,1:6)=(/ 0.900,0.000000,0.000000,0.845497,0.850160,0.937140 /)
XEXT_COEFF_550_LKT(93,19)=2.364600 !rg=16.4566 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,20,1:6)=(/ 0.000000,0.000000,1.967600,1.973300,1.981500,1.995400 /)
XPIZA_LKT(93,20,1:6)=(/ 0.900,0.900,0.959395,0.967990,0.917499,0.554839 /)
XCGA_LKT(93,20,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.851747,0.937157 /)
XEXT_COEFF_550_LKT(93,20)=0.000000 !rg=16.4566 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,1,1:6)=(/ 33.908000,34.053000,34.410000,35.171000,35.899000,36.594000 /)
XPIZA_LKT(94,1,1:6)=(/ 0.719037,0.872160,0.994277,0.993246,0.972139,0.534472 /)
XCGA_LKT(94,1,1:6)=(/ 0.891543,0.859857,0.831927,0.827883,0.822620,0.933717 /)
XEXT_COEFF_550_LKT(94,1)=34.347000 !rg=17.8284 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,2,1:6)=(/ 32.637000,32.937000,33.070000,33.489000,34.051000,35.002000 /)
XPIZA_LKT(94,2,1:6)=(/ 0.716258,0.871427,0.994119,0.993368,0.969389,0.535054 /)
XCGA_LKT(94,2,1:6)=(/ 0.892250,0.860737,0.833267,0.827590,0.821300,0.933847 /)
XEXT_COEFF_550_LKT(94,2)=33.126000 !rg=17.8284 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,3,1:6)=(/ 30.261000,30.499000,30.725000,31.081000,31.658000,32.363000 /)
XPIZA_LKT(94,3,1:6)=(/ 0.710160,0.867035,0.993649,0.993002,0.968731,0.536049 /)
XCGA_LKT(94,3,1:6)=(/ 0.893407,0.861407,0.833680,0.828573,0.823877,0.934060 /)
XEXT_COEFF_550_LKT(94,3)=30.479000 !rg=17.8284 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,4,1:6)=(/ 27.362000,27.542000,27.710000,27.888000,28.647000,29.169000 /)
XPIZA_LKT(94,4,1:6)=(/ 0.701991,0.861026,0.993577,0.993010,0.968664,0.537308 /)
XCGA_LKT(94,4,1:6)=(/ 0.895083,0.862547,0.834123,0.830273,0.824233,0.934323 /)
XEXT_COEFF_550_LKT(94,4)=27.540000 !rg=17.8284 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,5,1:6)=(/ 24.237000,24.414000,24.513000,24.815000,25.445000,25.759000 /)
XPIZA_LKT(94,5,1:6)=(/ 0.692359,0.853547,0.993136,0.992857,0.968885,0.538761 /)
XCGA_LKT(94,5,1:6)=(/ 0.897060,0.863870,0.835527,0.830037,0.825393,0.934623 /)
XEXT_COEFF_550_LKT(94,5)=24.395000 !rg=17.8284 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,6,1:6)=(/ 21.163000,21.285000,21.340000,21.707000,21.900000,22.416000 /)
XPIZA_LKT(94,6,1:6)=(/ 0.681553,0.844132,0.992789,0.992499,0.969299,0.540328 /)
XCGA_LKT(94,6,1:6)=(/ 0.899397,0.865127,0.836540,0.829113,0.828473,0.934937 /)
XEXT_COEFF_550_LKT(94,6)=21.268000 !rg=17.8284 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,7,1:6)=(/ 18.227000,18.322000,18.398000,18.512000,18.928000,19.243000 /)
XPIZA_LKT(94,7,1:6)=(/ 0.669799,0.833068,0.992178,0.992357,0.969391,0.541987 /)
XCGA_LKT(94,7,1:6)=(/ 0.901957,0.866870,0.837043,0.833103,0.829857,0.935267 /)
XEXT_COEFF_550_LKT(94,7)=18.322000 !rg=17.8284 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,8,1:6)=(/ 15.592000,15.666000,15.730000,15.874000,16.089000,16.390000 /)
XPIZA_LKT(94,8,1:6)=(/ 0.657934,0.820823,0.991560,0.991917,0.969355,0.543625 /)
XCGA_LKT(94,8,1:6)=(/ 0.905003,0.868957,0.838177,0.834613,0.832187,0.935580 /)
XEXT_COEFF_550_LKT(94,8)=15.673000 !rg=17.8284 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,9,1:6)=(/ 13.225000,13.287000,13.333000,13.486000,13.687000,13.848000 /)
XPIZA_LKT(94,9,1:6)=(/ 0.645816,0.806793,0.990634,0.991380,0.967985,0.545226 /)
XCGA_LKT(94,9,1:6)=(/ 0.908313,0.871200,0.839397,0.835407,0.832973,0.935873 /)
XEXT_COEFF_550_LKT(94,9)=13.291000 !rg=17.8284 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,10,1:6)=(/ 11.176000,11.211000,11.281000,11.364000,11.488000,11.657000 /)
XPIZA_LKT(94,10,1:6)=(/ 0.634168,0.791563,0.989626,0.990556,0.967126,0.546712 /)
XCGA_LKT(94,10,1:6)=(/ 0.911743,0.873493,0.840733,0.836803,0.835713,0.936137 /)
XEXT_COEFF_550_LKT(94,10)=11.241000 !rg=17.8284 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,11,1:6)=(/ 9.382400,9.413800,9.450200,9.500900,9.640400,9.752900 /)
XPIZA_LKT(94,11,1:6)=(/ 0.622777,0.775051,0.988275,0.989753,0.964747,0.548108 /)
XCGA_LKT(94,11,1:6)=(/ 0.915330,0.876237,0.841417,0.838333,0.836933,0.936370 /)
XEXT_COEFF_550_LKT(94,11)=9.425100 !rg=17.8284 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,12,1:6)=(/ 7.850500,7.870700,7.908200,7.954000,8.022900,8.132500 /)
XPIZA_LKT(94,12,1:6)=(/ 0.612118,0.757568,0.986674,0.988538,0.962763,0.549375 /)
XCGA_LKT(94,12,1:6)=(/ 0.919087,0.879127,0.842233,0.838913,0.839487,0.936567 /)
XEXT_COEFF_550_LKT(94,12)=7.889500 !rg=17.8284 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,13,1:6)=(/ 6.551200,6.567600,6.589900,6.636300,6.674400,6.764600 /)
XPIZA_LKT(94,13,1:6)=(/ 0.602404,0.739626,0.984736,0.987200,0.959679,0.550505 /)
XCGA_LKT(94,13,1:6)=(/ 0.000000,0.882413,0.842773,0.840467,0.841173,0.936733 /)
XEXT_COEFF_550_LKT(94,13)=6.573100 !rg=17.8284 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,14,1:6)=(/ 5.464600,5.477000,5.490100,5.523000,5.569300,5.626100 /)
XPIZA_LKT(94,14,1:6)=(/ 0.593696,0.721606,0.982462,0.985381,0.955895,0.551486 /)
XCGA_LKT(94,14,1:6)=(/ 0.000000,0.885877,0.843637,0.841003,0.842673,0.936863 /)
XEXT_COEFF_550_LKT(94,14)=5.478700 !rg=17.8284 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,15,1:6)=(/ 4.553900,4.561300,4.575300,4.604000,4.637300,4.674800 /)
XPIZA_LKT(94,15,1:6)=(/ 0.586056,0.703739,0.979726,0.983467,0.951244,0.552336 /)
XCGA_LKT(94,15,1:6)=(/ 0.000000,0.000000,0.844627,0.842530,0.844360,0.936967 /)
XEXT_COEFF_550_LKT(94,15)=4.566700 !rg=17.8284 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,16,1:6)=(/ 3.784600,3.791700,3.799800,3.814900,3.842800,3.875500 /)
XPIZA_LKT(94,16,1:6)=(/ 0.579126,0.686412,0.976461,0.980960,0.945891,0.553079 /)
XCGA_LKT(94,16,1:6)=(/ 0.000000,0.000000,0.845187,0.843173,0.845987,0.937043 /)
XEXT_COEFF_550_LKT(94,16)=3.792200 !rg=17.8284 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,17,1:6)=(/ 0.000000,3.154800,3.161600,3.171900,3.196000,3.218100 /)
XPIZA_LKT(94,17,1:6)=(/ 0.900,0.669999,0.972615,0.978071,0.939441,0.553699 /)
XCGA_LKT(94,17,1:6)=(/ 0.900,0.000000,0.000000,0.843900,0.847687,0.937100 /)
XEXT_COEFF_550_LKT(94,17)=3.155300 !rg=17.8284 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,18,1:6)=(/ 0.000000,2.625000,2.630100,2.640600,2.655800,2.672900 /)
XPIZA_LKT(94,18,1:6)=(/ 0.900,0.654624,0.968112,0.974701,0.931797,0.554220 /)
XCGA_LKT(94,18,1:6)=(/ 0.900,0.000000,0.000000,0.845077,0.849283,0.937133 /)
XEXT_COEFF_550_LKT(94,18)=2.627100 !rg=17.8284 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,19,1:6)=(/ 0.000000,0.000000,2.185000,2.191000,2.202200,2.216900 /)
XPIZA_LKT(94,19,1:6)=(/ 0.900,0.900,0.962772,0.970598,0.923021,0.554660 /)
XCGA_LKT(94,19,1:6)=(/ 0.900,0.900,0.000000,0.845953,0.850667,0.937157 /)
XEXT_COEFF_550_LKT(94,19)=2.181500 !rg=17.8284 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,20,1:6)=(/ 0.000000,0.000000,1.816100,1.821100,1.829700,1.839900 /)
XPIZA_LKT(94,20,1:6)=(/ 0.900,0.900,0.956575,0.965822,0.912994,0.555022 /)
XCGA_LKT(94,20,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.852567,0.937163 /)
XEXT_COEFF_550_LKT(94,20)=0.000000 !rg=17.8284 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,1,1:6)=(/ 31.244000,31.628000,31.843000,32.448000,32.732000,33.630000 /)
XPIZA_LKT(95,1,1:6)=(/ 0.710380,0.869150,0.993932,0.993441,0.968865,0.536203 /)
XCGA_LKT(95,1,1:6)=(/ 0.892513,0.862113,0.834873,0.828093,0.826273,0.934153 /)
XEXT_COEFF_550_LKT(95,1)=31.943000 !rg=19.3145 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,2,1:6)=(/ 30.122000,30.396000,30.598000,31.111000,31.320000,32.170000 /)
XPIZA_LKT(95,2,1:6)=(/ 0.708508,0.866835,0.993937,0.993075,0.967164,0.536746 /)
XCGA_LKT(95,2,1:6)=(/ 0.893970,0.862037,0.834793,0.825003,0.825457,0.934267 /)
XEXT_COEFF_550_LKT(95,2)=30.139000 !rg=19.3145 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,3,1:6)=(/ 27.939000,27.990000,28.251000,28.682000,29.363000,29.748000 /)
XPIZA_LKT(95,3,1:6)=(/ 0.702650,0.861729,0.993792,0.992882,0.967219,0.537673 /)
XCGA_LKT(95,3,1:6)=(/ 0.895120,0.861730,0.835560,0.829937,0.823933,0.934460 /)
XEXT_COEFF_550_LKT(95,3)=28.268000 !rg=19.3145 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,4,1:6)=(/ 25.251000,25.300000,25.553000,25.855000,26.426000,26.817000 /)
XPIZA_LKT(95,4,1:6)=(/ 0.694452,0.855300,0.993441,0.992985,0.968140,0.538846 /)
XCGA_LKT(95,4,1:6)=(/ 0.896813,0.862980,0.835693,0.829860,0.826470,0.934697 /)
XEXT_COEFF_550_LKT(95,4)=25.537000 !rg=19.3145 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,5,1:6)=(/ 22.345000,22.467000,22.623000,22.856000,23.241000,23.686000 /)
XPIZA_LKT(95,5,1:6)=(/ 0.684628,0.847757,0.992993,0.992773,0.969199,0.540198 /)
XCGA_LKT(95,5,1:6)=(/ 0.898583,0.864637,0.836173,0.831597,0.827497,0.934963 /)
XEXT_COEFF_550_LKT(95,5)=22.465000 !rg=19.3145 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,6,1:6)=(/ 19.525000,19.597000,19.726000,19.884000,20.212000,20.616000 /)
XPIZA_LKT(95,6,1:6)=(/ 0.674122,0.837829,0.992559,0.992551,0.969317,0.541655 /)
XCGA_LKT(95,6,1:6)=(/ 0.901117,0.865977,0.837750,0.832330,0.830550,0.935243 /)
XEXT_COEFF_550_LKT(95,6)=19.679000 !rg=19.3145 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,7,1:6)=(/ 16.826000,16.876000,17.043000,17.130000,17.420000,17.702000 /)
XPIZA_LKT(95,7,1:6)=(/ 0.662589,0.826371,0.991911,0.992158,0.969307,0.543196 /)
XCGA_LKT(95,7,1:6)=(/ 0.903900,0.867753,0.838270,0.833643,0.831923,0.935537 /)
XEXT_COEFF_550_LKT(95,7)=16.925000 !rg=19.3145 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,8,1:6)=(/ 14.385000,14.433000,14.546000,14.620000,14.820000,15.082000 /)
XPIZA_LKT(95,8,1:6)=(/ 0.650908,0.813527,0.991103,0.991690,0.968556,0.544716 /)
XCGA_LKT(95,8,1:6)=(/ 0.906873,0.869920,0.839397,0.834603,0.833750,0.935813 /)
XEXT_COEFF_550_LKT(95,8)=14.491000 !rg=19.3145 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,9,1:6)=(/ 12.199000,12.243000,12.307000,12.386000,12.536000,12.745000 /)
XPIZA_LKT(95,9,1:6)=(/ 0.639111,0.799126,0.990148,0.990992,0.967562,0.546199 /)
XCGA_LKT(95,9,1:6)=(/ 0.910130,0.872397,0.839907,0.836120,0.834407,0.936077 /)
XEXT_COEFF_550_LKT(95,9)=12.256000 !rg=19.3145 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,10,1:6)=(/ 10.307000,10.341000,10.379000,10.466000,10.592000,10.731000 /)
XPIZA_LKT(95,10,1:6)=(/ 0.627759,0.783642,0.989002,0.990207,0.965934,0.547575 /)
XCGA_LKT(95,10,1:6)=(/ 0.913647,0.874870,0.840893,0.837627,0.835550,0.936307 /)
XEXT_COEFF_550_LKT(95,10)=10.355000 !rg=19.3145 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,11,1:6)=(/ 8.654900,8.678900,8.725600,8.761500,8.885800,8.980400 /)
XPIZA_LKT(95,11,1:6)=(/ 0.616930,0.766650,0.987560,0.989194,0.964070,0.548865 /)
XCGA_LKT(95,11,1:6)=(/ 0.917257,0.877587,0.841807,0.838197,0.838613,0.936510 /)
XEXT_COEFF_550_LKT(95,11)=8.697600 !rg=19.3145 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,12,1:6)=(/ 7.243500,7.264100,7.287200,7.336700,7.388700,7.490000 /)
XPIZA_LKT(95,12,1:6)=(/ 0.606944,0.749196,0.985868,0.988000,0.961480,0.550033 /)
XCGA_LKT(95,12,1:6)=(/ 0.000000,0.880697,0.842547,0.840223,0.840163,0.936683 /)
XEXT_COEFF_550_LKT(95,12)=7.270300 !rg=19.3145 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,13,1:6)=(/ 6.044000,6.057800,6.076700,6.108800,6.167200,6.231400 /)
XPIZA_LKT(95,13,1:6)=(/ 0.597769,0.730999,0.983769,0.986490,0.958214,0.551073 /)
XCGA_LKT(95,13,1:6)=(/ 0.000000,0.884020,0.843173,0.840570,0.842070,0.936823 /)
XEXT_COEFF_550_LKT(95,13)=6.062700 !rg=19.3145 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,14,1:6)=(/ 5.042800,5.051100,5.066300,5.099400,5.130900,5.183600 /)
XPIZA_LKT(95,14,1:6)=(/ 0.589646,0.713044,0.981264,0.984646,0.953925,0.551972 /)
XCGA_LKT(95,14,1:6)=(/ 0.000000,0.000000,0.844063,0.841980,0.843190,0.936933 /)
XEXT_COEFF_550_LKT(95,14)=5.058800 !rg=19.3145 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,15,1:6)=(/ 4.201600,4.209500,4.221400,4.237000,4.266100,4.307900 /)
XPIZA_LKT(95,15,1:6)=(/ 0.582396,0.695646,0.978388,0.982407,0.948974,0.552751 /)
XCGA_LKT(95,15,1:6)=(/ 0.000000,0.000000,0.844957,0.842843,0.844820,0.937020 /)
XEXT_COEFF_550_LKT(95,15)=4.213300 !rg=19.3145 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,16,1:6)=(/ 3.492200,3.497100,3.506000,3.517000,3.551300,3.571900 /)
XPIZA_LKT(95,16,1:6)=(/ 0.576009,0.678488,0.974856,0.979747,0.943249,0.553429 /)
XCGA_LKT(95,16,1:6)=(/ 0.000000,0.000000,0.845680,0.843627,0.847050,0.937083 /)
XEXT_COEFF_550_LKT(95,16)=3.500400 !rg=19.3145 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,17,1:6)=(/ 0.000000,2.910500,2.915900,2.928300,2.944500,2.966500 /)
XPIZA_LKT(95,17,1:6)=(/ 0.900,0.662561,0.970724,0.976693,0.936285,0.553993 /)
XCGA_LKT(95,17,1:6)=(/ 0.900,0.000000,0.000000,0.844593,0.848577,0.937127 /)
XEXT_COEFF_550_LKT(95,17)=2.912400 !rg=19.3145 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,18,1:6)=(/ 0.000000,2.422800,2.427200,2.434000,2.446000,2.464200 /)
XPIZA_LKT(95,18,1:6)=(/ 0.900,0.647837,0.965908,0.972979,0.928043,0.554464 /)
XCGA_LKT(95,18,1:6)=(/ 0.900,0.000000,0.000000,0.845413,0.849750,0.937153 /)
XEXT_COEFF_550_LKT(95,18)=2.424200 !rg=19.3145 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,19,1:6)=(/ 0.000000,0.000000,2.016100,2.020500,2.034400,2.044000 /)
XPIZA_LKT(95,19,1:6)=(/ 0.900,0.900,0.960197,0.968602,0.918858,0.554861 /)
XCGA_LKT(95,19,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.851760,0.937167 /)
XEXT_COEFF_550_LKT(95,19)=2.013600 !rg=19.3145 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,20,1:6)=(/ 0.000000,0.000000,1.675300,1.679400,1.687000,1.696700 /)
XPIZA_LKT(95,20,1:6)=(/ 0.900,0.900,0.953543,0.963492,0.908157,0.555186 /)
XCGA_LKT(95,20,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.853470,0.937167 /)
XEXT_COEFF_550_LKT(95,20)=0.000000 !rg=19.3145 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,1,1:6)=(/ 29.057000,29.134000,29.175000,29.562000,30.274000,30.913000 /)
XPIZA_LKT(96,1,1:6)=(/ 0.704649,0.863921,0.993612,0.993500,0.962865,0.537822 /)
XCGA_LKT(96,1,1:6)=(/ 0.894290,0.860303,0.833870,0.829613,0.825200,0.934543 /)
XEXT_COEFF_550_LKT(96,1)=29.582000 !rg=20.9245 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,2,1:6)=(/ 27.813000,27.870000,28.168000,28.403000,28.977000,29.573000 /)
XPIZA_LKT(96,2,1:6)=(/ 0.700939,0.861241,0.993776,0.993370,0.965693,0.538328 /)
XCGA_LKT(96,2,1:6)=(/ 0.895510,0.861683,0.836557,0.829287,0.828353,0.934647 /)
XEXT_COEFF_550_LKT(96,2)=28.112000 !rg=20.9245 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,3,1:6)=(/ 25.762000,25.903000,26.016000,26.279000,27.009000,27.350000 /)
XPIZA_LKT(96,3,1:6)=(/ 0.694541,0.857054,0.993555,0.993246,0.967197,0.539192 /)
XCGA_LKT(96,3,1:6)=(/ 0.896760,0.863453,0.835327,0.830273,0.825960,0.934817 /)
XEXT_COEFF_550_LKT(96,3)=25.847000 !rg=20.9245 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,4,1:6)=(/ 23.275000,23.426000,23.430000,23.731000,24.383000,24.659000 /)
XPIZA_LKT(96,4,1:6)=(/ 0.686554,0.850483,0.993235,0.992964,0.968137,0.540283 /)
XCGA_LKT(96,4,1:6)=(/ 0.898330,0.864490,0.836030,0.831387,0.826023,0.935027 /)
XEXT_COEFF_550_LKT(96,4)=23.370000 !rg=20.9245 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,5,1:6)=(/ 20.625000,20.687000,20.850000,21.064000,21.462000,21.784000 /)
XPIZA_LKT(96,5,1:6)=(/ 0.677088,0.841548,0.992763,0.992562,0.968710,0.541539 /)
XCGA_LKT(96,5,1:6)=(/ 0.900403,0.865480,0.837477,0.832253,0.829213,0.935267 /)
XEXT_COEFF_550_LKT(96,5)=20.809000 !rg=20.9245 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,6,1:6)=(/ 18.003000,18.102000,18.153000,18.406000,18.542000,18.965000 /)
XPIZA_LKT(96,6,1:6)=(/ 0.666515,0.831619,0.992235,0.992288,0.968675,0.542892 /)
XCGA_LKT(96,6,1:6)=(/ 0.902773,0.867457,0.837873,0.833253,0.831690,0.935517 /)
XEXT_COEFF_550_LKT(96,6)=18.094000 !rg=20.9245 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,7,1:6)=(/ 15.516000,15.595000,15.635000,15.764000,16.093000,16.288000 /)
XPIZA_LKT(96,7,1:6)=(/ 0.655240,0.819852,0.991533,0.992001,0.968194,0.544322 /)
XCGA_LKT(96,7,1:6)=(/ 0.905620,0.869303,0.838773,0.835170,0.831713,0.935777 /)
XEXT_COEFF_550_LKT(96,7)=15.603000 !rg=20.9245 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,8,1:6)=(/ 13.265000,13.326000,13.366000,13.495000,13.614000,13.880000 /)
XPIZA_LKT(96,8,1:6)=(/ 0.643853,0.806296,0.990706,0.991423,0.968347,0.545731 /)
XCGA_LKT(96,8,1:6)=(/ 0.908713,0.871397,0.839550,0.836223,0.834597,0.936023 /)
XEXT_COEFF_550_LKT(96,8)=13.332000 !rg=20.9245 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,9,1:6)=(/ 11.257000,11.288000,11.347000,11.453000,11.613000,11.732000 /)
XPIZA_LKT(96,9,1:6)=(/ 0.632662,0.791320,0.989618,0.990691,0.966942,0.547104 /)
XCGA_LKT(96,9,1:6)=(/ 0.912083,0.873617,0.840483,0.837277,0.836097,0.936250 /)
XEXT_COEFF_550_LKT(96,9)=11.325000 !rg=20.9245 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,10,1:6)=(/ 9.508000,9.536800,9.574800,9.638000,9.766500,9.880100 /)
XPIZA_LKT(96,10,1:6)=(/ 0.621853,0.775383,0.988401,0.989819,0.965344,0.548376 /)
XCGA_LKT(96,10,1:6)=(/ 0.915443,0.876277,0.841307,0.837700,0.837953,0.936453 /)
XEXT_COEFF_550_LKT(96,10)=9.543200 !rg=20.9245 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,11,1:6)=(/ 7.985600,8.006100,8.044900,8.101800,8.169900,8.270100 /)
XPIZA_LKT(96,11,1:6)=(/ 0.611528,0.758239,0.986825,0.988686,0.963009,0.549565 /)
XCGA_LKT(96,11,1:6)=(/ 0.919123,0.879027,0.841973,0.839467,0.839633,0.936630 /)
XEXT_COEFF_550_LKT(96,11)=8.014300 !rg=20.9245 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,12,1:6)=(/ 6.682100,6.699300,6.716200,6.764900,6.828100,6.899000 /)
XPIZA_LKT(96,12,1:6)=(/ 0.601991,0.740637,0.984991,0.987271,0.960059,0.550641 /)
XCGA_LKT(96,12,1:6)=(/ 0.000000,0.882270,0.842850,0.840233,0.841147,0.936780 /)
XEXT_COEFF_550_LKT(96,12)=6.697400 !rg=20.9245 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,13,1:6)=(/ 5.577200,5.586600,5.606600,5.639600,5.680000,5.740800 /)
XPIZA_LKT(96,13,1:6)=(/ 0.593362,0.722415,0.982736,0.985698,0.956394,0.551596 /)
XCGA_LKT(96,13,1:6)=(/ 0.000000,0.885610,0.843807,0.841463,0.842650,0.936900 /)
XEXT_COEFF_550_LKT(96,13)=5.591600 !rg=20.9245 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,14,1:6)=(/ 4.652700,4.662500,4.676100,4.694200,4.726700,4.776400 /)
XPIZA_LKT(96,14,1:6)=(/ 0.585716,0.704816,0.980083,0.983709,0.951972,0.552420 /)
XCGA_LKT(96,14,1:6)=(/ 0.000000,0.000000,0.844473,0.842243,0.844430,0.936993 /)
XEXT_COEFF_550_LKT(96,14)=4.666700 !rg=20.9245 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,15,1:6)=(/ 3.877100,3.883300,3.897000,3.906900,3.943000,3.970100 /)
XPIZA_LKT(96,15,1:6)=(/ 0.579019,0.687498,0.976927,0.981322,0.946709,0.553132 /)
XCGA_LKT(96,15,1:6)=(/ 0.000000,0.000000,0.845413,0.843120,0.846133,0.937067 /)
XEXT_COEFF_550_LKT(96,15)=3.886800 !rg=20.9245 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,16,1:6)=(/ 0.000000,3.227200,3.233700,3.248400,3.267200,3.292300 /)
XPIZA_LKT(96,16,1:6)=(/ 0.900,0.670828,0.973117,0.978489,0.940180,0.553749 /)
XCGA_LKT(96,16,1:6)=(/ 0.900,0.000000,0.000000,0.844227,0.847617,0.937117 /)
XEXT_COEFF_550_LKT(96,16)=3.228800 !rg=20.9245 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,17,1:6)=(/ 0.000000,2.686100,2.690100,2.700200,2.713300,2.734700 /)
XPIZA_LKT(96,17,1:6)=(/ 0.900,0.655427,0.968704,0.975127,0.932723,0.554260 /)
XCGA_LKT(96,17,1:6)=(/ 0.900,0.000000,0.000000,0.844850,0.848877,0.937150 /)
XEXT_COEFF_550_LKT(96,17)=2.687400 !rg=20.9245 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,18,1:6)=(/ 0.000000,0.000000,2.239400,2.245300,2.259000,2.271900 /)
XPIZA_LKT(96,18,1:6)=(/ 0.900,0.900,0.963527,0.971165,0.924248,0.554686 /)
XCGA_LKT(96,18,1:6)=(/ 0.900,0.900,0.000000,0.845773,0.850863,0.937167 /)
XEXT_COEFF_550_LKT(96,18)=2.236300 !rg=20.9245 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,19,1:6)=(/ 0.000000,0.000000,1.860000,1.865400,1.872800,1.884700 /)
XPIZA_LKT(96,19,1:6)=(/ 0.900,0.900,0.957415,0.966484,0.914258,0.555043 /)
XCGA_LKT(96,19,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.852430,0.937173 /)
XEXT_COEFF_550_LKT(96,19)=0.000000 !rg=20.9245 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,20,1:6)=(/ 0.000000,0.000000,1.546000,1.549600,1.555200,1.564600 /)
XPIZA_LKT(96,20,1:6)=(/ 0.900,0.900,0.950340,0.960987,0.902968,0.555332 /)
XCGA_LKT(96,20,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.853980,0.937167 /)
XEXT_COEFF_550_LKT(96,20)=0.000000 !rg=20.9245 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,1,1:6)=(/ 26.943000,27.058000,27.126000,27.416000,27.906000,28.421000 /)
XPIZA_LKT(97,1,1:6)=(/ 0.697682,0.859316,0.993796,0.993591,0.962135,0.539336 /)
XCGA_LKT(97,1,1:6)=(/ 0.896830,0.863100,0.836997,0.832063,0.830930,0.934893 /)
XEXT_COEFF_550_LKT(97,1)=26.504000 !rg=22.6687 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,2,1:6)=(/ 25.642000,25.834000,25.898000,26.267000,26.457000,27.191000 /)
XPIZA_LKT(97,2,1:6)=(/ 0.692728,0.856376,0.993547,0.993370,0.965455,0.539807 /)
XCGA_LKT(97,2,1:6)=(/ 0.897153,0.863830,0.835633,0.830157,0.829283,0.934987 /)
XEXT_COEFF_550_LKT(97,2)=25.722000 !rg=22.6687 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,3,1:6)=(/ 23.733000,23.826000,24.141000,24.239000,24.706000,25.150000 /)
XPIZA_LKT(97,3,1:6)=(/ 0.686297,0.850854,0.993348,0.993125,0.967601,0.540609 /)
XCGA_LKT(97,3,1:6)=(/ 0.898280,0.864120,0.837007,0.831763,0.828877,0.935137 /)
XEXT_COEFF_550_LKT(97,3)=23.905000 !rg=22.6687 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,4,1:6)=(/ 21.453000,21.532000,21.799000,21.913000,22.298000,22.679000 /)
XPIZA_LKT(97,4,1:6)=(/ 0.678407,0.843975,0.992889,0.992598,0.968168,0.541624 /)
XCGA_LKT(97,4,1:6)=(/ 0.899890,0.865173,0.837330,0.831977,0.828903,0.935327 /)
XEXT_COEFF_550_LKT(97,4)=21.617000 !rg=22.6687 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,5,1:6)=(/ 19.020000,19.115000,19.172000,19.392000,19.795000,20.039000 /)
XPIZA_LKT(97,5,1:6)=(/ 0.669231,0.835577,0.992493,0.992482,0.968797,0.542790 /)
XCGA_LKT(97,5,1:6)=(/ 0.902137,0.866957,0.837447,0.832863,0.828950,0.935537 /)
XEXT_COEFF_550_LKT(97,5)=19.102000 !rg=22.6687 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,6,1:6)=(/ 16.600000,16.657000,16.795000,16.878000,17.193000,17.448000 /)
XPIZA_LKT(97,6,1:6)=(/ 0.658790,0.824517,0.991861,0.991981,0.968846,0.544045 /)
XCGA_LKT(97,6,1:6)=(/ 0.904553,0.868407,0.838740,0.833213,0.831877,0.935760 /)
XEXT_COEFF_550_LKT(97,6)=16.697000 !rg=22.6687 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,7,1:6)=(/ 14.308000,14.365000,14.468000,14.531000,14.790000,14.989000 /)
XPIZA_LKT(97,7,1:6)=(/ 0.647985,0.812382,0.991125,0.991645,0.968282,0.545371 /)
XCGA_LKT(97,7,1:6)=(/ 0.907417,0.870347,0.839637,0.834800,0.833297,0.935990 /)
XEXT_COEFF_550_LKT(97,7)=14.402000 !rg=22.6687 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,8,1:6)=(/ 12.238000,12.279000,12.315000,12.409000,12.579000,12.775000 /)
XPIZA_LKT(97,8,1:6)=(/ 0.637023,0.798506,0.990173,0.990984,0.967696,0.546674 /)
XCGA_LKT(97,8,1:6)=(/ 0.910593,0.872580,0.839933,0.835633,0.835443,0.936207 /)
XEXT_COEFF_550_LKT(97,8)=12.295000 !rg=22.6687 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,9,1:6)=(/ 10.381000,10.417000,10.466000,10.528000,10.661000,10.801000 /)
XPIZA_LKT(97,9,1:6)=(/ 0.626307,0.783403,0.989038,0.990279,0.966234,0.547944 /)
XCGA_LKT(97,9,1:6)=(/ 0.913873,0.875073,0.840790,0.837287,0.836173,0.936407 /)
XEXT_COEFF_550_LKT(97,9)=10.424000 !rg=22.6687 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,10,1:6)=(/ 8.773700,8.794800,8.835200,8.922400,8.997500,9.097800 /)
XPIZA_LKT(97,10,1:6)=(/ 0.616116,0.766971,0.987733,0.989377,0.964409,0.549117 /)
XCGA_LKT(97,10,1:6)=(/ 0.917413,0.877593,0.841920,0.839190,0.838437,0.936583 /)
XEXT_COEFF_550_LKT(97,10)=8.814500 !rg=22.6687 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,11,1:6)=(/ 7.368300,7.387200,7.402200,7.449900,7.533100,7.616900 /)
XPIZA_LKT(97,11,1:6)=(/ 0.606368,0.749827,0.986023,0.988074,0.961746,0.550213 /)
XCGA_LKT(97,11,1:6)=(/ 0.000000,0.880647,0.842520,0.839697,0.840317,0.936737 /)
XEXT_COEFF_550_LKT(97,11)=7.391900 !rg=22.6687 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,12,1:6)=(/ 6.166500,6.177500,6.195100,6.242900,6.289600,6.355300 /)
XPIZA_LKT(97,12,1:6)=(/ 0.597353,0.731947,0.984016,0.986680,0.958548,0.551201 /)
XCGA_LKT(97,12,1:6)=(/ 0.000000,0.883730,0.843360,0.841357,0.841690,0.936867 /)
XEXT_COEFF_550_LKT(97,12)=6.185700 !rg=22.6687 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,13,1:6)=(/ 5.146100,5.157200,5.168800,5.196300,5.229300,5.289300 /)
XPIZA_LKT(97,13,1:6)=(/ 0.589219,0.714063,0.981597,0.984873,0.954589,0.552077 /)
XCGA_LKT(97,13,1:6)=(/ 0.000000,0.000000,0.843957,0.841873,0.843677,0.936970 /)
XEXT_COEFF_550_LKT(97,13)=5.160300 !rg=22.6687 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,14,1:6)=(/ 4.293400,4.301400,4.313900,4.330000,4.368500,4.401500 /)
XPIZA_LKT(97,14,1:6)=(/ 0.582124,0.696585,0.978764,0.982668,0.949770,0.552831 /)
XCGA_LKT(97,14,1:6)=(/ 0.000000,0.000000,0.844953,0.842660,0.845737,0.937047 /)
XEXT_COEFF_550_LKT(97,14)=4.305300 !rg=22.6687 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,15,1:6)=(/ 3.577600,3.583500,3.590000,3.606800,3.628800,3.659100 /)
XPIZA_LKT(97,15,1:6)=(/ 0.575856,0.679674,0.975324,0.980151,0.943866,0.553479 /)
XCGA_LKT(97,15,1:6)=(/ 0.000000,0.000000,0.845577,0.843697,0.846187,0.937103 /)
XEXT_COEFF_550_LKT(97,15)=3.585400 !rg=22.6687 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,16,1:6)=(/ 0.000000,2.977900,2.981500,2.993300,3.013100,3.034900 /)
XPIZA_LKT(97,16,1:6)=(/ 0.900,0.663348,0.971258,0.977095,0.937103,0.554041 /)
XCGA_LKT(97,16,1:6)=(/ 0.900,0.000000,0.000000,0.844563,0.848213,0.937143 /)
XEXT_COEFF_550_LKT(97,16)=2.979100 !rg=22.6687 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,17,1:6)=(/ 0.000000,2.478600,2.482200,2.490400,2.506100,2.521100 /)
XPIZA_LKT(97,17,1:6)=(/ 0.900,0.648505,0.966516,0.973482,0.929189,0.554503 /)
XCGA_LKT(97,17,1:6)=(/ 0.900,0.000000,0.000000,0.845330,0.850020,0.937167 /)
XEXT_COEFF_550_LKT(97,17)=2.480000 !rg=22.6687 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,18,1:6)=(/ 0.000000,0.000000,2.065400,2.071400,2.080300,2.094800 /)
XPIZA_LKT(97,18,1:6)=(/ 0.900,0.900,0.960970,0.969214,0.919945,0.554886 /)
XCGA_LKT(97,18,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.851290,0.937177 /)
XEXT_COEFF_550_LKT(97,18)=2.063600 !rg=22.6687 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,19,1:6)=(/ 0.000000,0.000000,1.715600,1.719800,1.727600,1.738000 /)
XPIZA_LKT(97,19,1:6)=(/ 0.900,0.900,0.954436,0.964175,0.909513,0.555206 /)
XCGA_LKT(97,19,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.853170,0.937177 /)
XEXT_COEFF_550_LKT(97,19)=0.000000 !rg=22.6687 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,20,1:6)=(/ 0.000000,0.000000,0.000000,1.429900,1.436000,1.442900 /)
XPIZA_LKT(97,20,1:6)=(/ 0.900,0.900,0.900,0.958307,0.897667,0.555464 /)
XCGA_LKT(97,20,1:6)=(/ 0.900,0.900,0.900,0.000000,0.855040,0.937167 /)
XEXT_COEFF_550_LKT(97,20)=0.000000 !rg=22.6687 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,1,1:6)=(/ 24.757000,24.745000,24.869000,25.077000,26.126000,26.135000 /)
XPIZA_LKT(98,1,1:6)=(/ 0.689265,0.853487,0.993474,0.993186,0.964157,0.540749 /)
XCGA_LKT(98,1,1:6)=(/ 0.898633,0.865067,0.834387,0.830470,0.830640,0.935207 /)
XEXT_COEFF_550_LKT(98,1)=24.604000 !rg=24.5583 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,2,1:6)=(/ 23.624000,23.709000,23.990000,24.055000,24.502000,25.006000 /)
XPIZA_LKT(98,2,1:6)=(/ 0.684488,0.850465,0.993246,0.993261,0.965915,0.541186 /)
XCGA_LKT(98,2,1:6)=(/ 0.898647,0.864580,0.836913,0.832327,0.830557,0.935290 /)
XEXT_COEFF_550_LKT(98,2)=23.736000 !rg=24.5583 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,3,1:6)=(/ 21.889000,22.030000,22.079000,22.364000,23.018000,23.131000 /)
XPIZA_LKT(98,3,1:6)=(/ 0.678266,0.845258,0.993017,0.993024,0.966806,0.541932 /)
XCGA_LKT(98,3,1:6)=(/ 0.899900,0.865637,0.836963,0.830957,0.827637,0.935423 /)
XEXT_COEFF_550_LKT(98,3)=22.078000 !rg=24.5583 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,4,1:6)=(/ 19.791000,19.903000,19.968000,20.247000,20.675000,20.861000 /)
XPIZA_LKT(98,4,1:6)=(/ 0.670585,0.838068,0.992570,0.992619,0.968215,0.542873 /)
XCGA_LKT(98,4,1:6)=(/ 0.901720,0.866743,0.837527,0.832280,0.830133,0.935590 /)
XEXT_COEFF_550_LKT(98,4)=19.920000 !rg=24.5583 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,5,1:6)=(/ 17.535000,17.600000,17.776000,17.816000,18.093000,18.436000 /)
XPIZA_LKT(98,5,1:6)=(/ 0.661311,0.828500,0.992189,0.992391,0.969104,0.543955 /)
XCGA_LKT(98,5,1:6)=(/ 0.903900,0.867960,0.838693,0.834293,0.831840,0.935777 /)
XEXT_COEFF_550_LKT(98,5)=17.651000 !rg=24.5583 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,6,1:6)=(/ 15.315000,15.371000,15.452000,15.624000,15.729000,16.056000 /)
XPIZA_LKT(98,6,1:6)=(/ 0.651521,0.817763,0.991508,0.992016,0.968904,0.545118 /)
XCGA_LKT(98,6,1:6)=(/ 0.906447,0.869503,0.839020,0.835410,0.833880,0.935973 /)
XEXT_COEFF_550_LKT(98,6)=15.422000 !rg=24.5583 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,7,1:6)=(/ 13.200000,13.245000,13.310000,13.446000,13.608000,13.795000 /)
XPIZA_LKT(98,7,1:6)=(/ 0.641038,0.804742,0.990632,0.991321,0.968128,0.546345 /)
XCGA_LKT(98,7,1:6)=(/ 0.909307,0.871557,0.839473,0.836277,0.835707,0.936177 /)
XEXT_COEFF_550_LKT(98,7)=13.267000 !rg=24.5583 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,8,1:6)=(/ 11.292000,11.324000,11.386000,11.459000,11.586000,11.760000 /)
XPIZA_LKT(98,8,1:6)=(/ 0.630691,0.790573,0.989666,0.990696,0.966966,0.547550 /)
XCGA_LKT(98,8,1:6)=(/ 0.912467,0.873740,0.840690,0.837237,0.835710,0.936367 /)
XEXT_COEFF_550_LKT(98,8)=11.354000 !rg=24.5583 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,9,1:6)=(/ 9.578300,9.603600,9.668200,9.685500,9.813400,9.945000 /)
XPIZA_LKT(98,9,1:6)=(/ 0.620241,0.774950,0.988465,0.989796,0.965598,0.548722 /)
XCGA_LKT(98,9,1:6)=(/ 0.915873,0.876383,0.841623,0.838077,0.838217,0.936543 /)
XEXT_COEFF_550_LKT(98,9)=9.630100 !rg=24.5583 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,10,1:6)=(/ 8.094400,8.117000,8.140700,8.187700,8.256800,8.378400 /)
XPIZA_LKT(98,10,1:6)=(/ 0.610555,0.758637,0.986956,0.988818,0.963235,0.549803 /)
XCGA_LKT(98,10,1:6)=(/ 0.919303,0.879170,0.842013,0.839387,0.839197,0.936697 /)
XEXT_COEFF_550_LKT(98,10)=8.125900 !rg=24.5583 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,11,1:6)=(/ 6.798100,6.813200,6.836500,6.881300,6.945400,7.016000 /)
XPIZA_LKT(98,11,1:6)=(/ 0.601328,0.741176,0.985172,0.987549,0.960519,0.550810 /)
XCGA_LKT(98,11,1:6)=(/ 0.000000,0.882107,0.843067,0.840757,0.841610,0.936830 /)
XEXT_COEFF_550_LKT(98,11)=6.821200 !rg=24.5583 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,12,1:6)=(/ 5.689200,5.702100,5.718500,5.745300,5.791000,5.855100 /)
XPIZA_LKT(98,12,1:6)=(/ 0.592938,0.723450,0.983010,0.985875,0.956938,0.551717 /)
XCGA_LKT(98,12,1:6)=(/ 0.000000,0.885537,0.843777,0.841403,0.842847,0.936940 /)
XEXT_COEFF_550_LKT(98,12)=5.708900 !rg=24.5583 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,13,1:6)=(/ 4.748600,4.757700,4.774000,4.785600,4.837600,4.873800 /)
XPIZA_LKT(98,13,1:6)=(/ 0.585340,0.705617,0.980388,0.983963,0.952649,0.552519 /)
XCGA_LKT(98,13,1:6)=(/ 0.000000,0.000000,0.844593,0.842260,0.844887,0.937027 /)
XEXT_COEFF_550_LKT(98,13)=4.763900 !rg=24.5583 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,14,1:6)=(/ 3.961500,3.968800,3.976300,3.996900,4.018700,4.056400 /)
XPIZA_LKT(98,14,1:6)=(/ 0.578690,0.688374,0.977313,0.981660,0.947242,0.553207 /)
XCGA_LKT(98,14,1:6)=(/ 0.000000,0.000000,0.845247,0.843333,0.845910,0.937090 /)
XEXT_COEFF_550_LKT(98,14)=3.969600 !rg=24.5583 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,15,1:6)=(/ 0.000000,3.306500,3.311700,3.326500,3.351000,3.372700 /)
XPIZA_LKT(98,15,1:6)=(/ 0.900,0.671829,0.973613,0.978898,0.941154,0.553797 /)
XCGA_LKT(98,15,1:6)=(/ 0.900,0.000000,0.845893,0.843960,0.847673,0.937133 /)
XEXT_COEFF_550_LKT(98,15)=3.306300 !rg=24.5583 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,16,1:6)=(/ 0.000000,2.748000,2.752900,2.763400,2.779700,2.797700 /)
XPIZA_LKT(98,16,1:6)=(/ 0.900,0.656120,0.969265,0.975579,0.933831,0.554306 /)
XCGA_LKT(98,16,1:6)=(/ 0.900,0.000000,0.000000,0.845003,0.849373,0.937163 /)
XEXT_COEFF_550_LKT(98,16)=2.748800 !rg=24.5583 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,17,1:6)=(/ 0.000000,0.000000,2.289900,2.297500,2.309100,2.324400 /)
XPIZA_LKT(98,17,1:6)=(/ 0.900,0.900,0.964182,0.971689,0.925232,0.554723 /)
XCGA_LKT(98,17,1:6)=(/ 0.900,0.900,0.000000,0.845843,0.850550,0.937180 /)
XEXT_COEFF_550_LKT(98,17)=2.287100 !rg=24.5583 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,18,1:6)=(/ 0.000000,0.000000,1.905700,1.911700,1.921000,1.931500 /)
XPIZA_LKT(98,18,1:6)=(/ 0.900,0.900,0.958257,0.967136,0.915702,0.555067 /)
XCGA_LKT(98,18,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.852313,0.937183 /)
XEXT_COEFF_550_LKT(98,18)=0.000000 !rg=24.5583 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,19,1:6)=(/ 0.000000,0.000000,1.584000,1.587600,1.594100,1.602700 /)
XPIZA_LKT(98,19,1:6)=(/ 0.900,0.900,0.951300,0.961747,0.904484,0.555353 /)
XCGA_LKT(98,19,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.854170,0.937177 /)
XEXT_COEFF_550_LKT(98,19)=0.000000 !rg=24.5583 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,20,1:6)=(/ 0.000000,0.000000,0.000000,1.319300,1.324400,1.330800 /)
XPIZA_LKT(98,20,1:6)=(/ 0.900,0.900,0.900,0.955468,0.891981,0.555581 /)
XCGA_LKT(98,20,1:6)=(/ 0.900,0.900,0.900,0.000000,0.855827,0.937163 /)
XEXT_COEFF_550_LKT(98,20)=0.000000 !rg=24.5583 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,1,1:6)=(/ 22.756000,22.611000,23.353000,23.174000,24.129000,24.038000 /)
XPIZA_LKT(99,1,1:6)=(/ 0.680044,0.846362,0.993304,0.993384,0.966069,0.542066 /)
XCGA_LKT(99,1,1:6)=(/ 0.900190,0.864757,0.840393,0.834987,0.829477,0.935490 /)
XEXT_COEFF_550_LKT(99,1)=23.267000 !rg=26.6054 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,2,1:6)=(/ 21.775000,21.851000,21.978000,22.188000,22.477000,23.000000 /)
XPIZA_LKT(99,2,1:6)=(/ 0.676166,0.844397,0.993040,0.993005,0.966730,0.542472 /)
XCGA_LKT(99,2,1:6)=(/ 0.900413,0.865363,0.837320,0.832047,0.830987,0.935560 /)
XEXT_COEFF_550_LKT(99,2)=21.978000 !rg=26.6054 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,3,1:6)=(/ 20.192000,20.268000,20.446000,20.650000,20.886000,21.278000 /)
XPIZA_LKT(99,3,1:6)=(/ 0.670367,0.838737,0.992813,0.992803,0.968631,0.543164 /)
XCGA_LKT(99,3,1:6)=(/ 0.901607,0.866480,0.838583,0.833373,0.830470,0.935680 /)
XEXT_COEFF_550_LKT(99,3)=20.290000 !rg=26.6054 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,4,1:6)=(/ 18.258000,18.319000,18.427000,18.599000,18.889000,19.193000 /)
XPIZA_LKT(99,4,1:6)=(/ 0.662911,0.831248,0.992278,0.992593,0.969009,0.544036 /)
XCGA_LKT(99,4,1:6)=(/ 0.903427,0.867690,0.838507,0.834303,0.831587,0.935827 /)
XEXT_COEFF_550_LKT(99,4)=18.351000 !rg=26.6054 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,5,1:6)=(/ 16.178000,16.254000,16.303000,16.495000,16.792000,16.964000 /)
XPIZA_LKT(99,5,1:6)=(/ 0.653890,0.821780,0.991728,0.992246,0.968592,0.545039 /)
XCGA_LKT(99,5,1:6)=(/ 0.905757,0.869247,0.839027,0.834493,0.832110,0.935993 /)
XEXT_COEFF_550_LKT(99,5)=16.265000 !rg=26.6054 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,6,1:6)=(/ 14.134000,14.182000,14.235000,14.353000,14.522000,14.777000 /)
XPIZA_LKT(99,6,1:6)=(/ 0.644511,0.810138,0.991019,0.991691,0.968229,0.546114 /)
XCGA_LKT(99,6,1:6)=(/ 0.908320,0.870947,0.839390,0.834873,0.834493,0.936167 /)
XEXT_COEFF_550_LKT(99,6)=14.171000 !rg=26.6054 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,7,1:6)=(/ 12.181000,12.221000,12.245000,12.327000,12.509000,12.698000 /)
XPIZA_LKT(99,7,1:6)=(/ 0.634380,0.797060,0.990174,0.990963,0.967606,0.547248 /)
XCGA_LKT(99,7,1:6)=(/ 0.911210,0.873013,0.840340,0.836610,0.835503,0.936343 /)
XEXT_COEFF_550_LKT(99,7)=12.230000 !rg=26.6054 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,8,1:6)=(/ 10.416000,10.456000,10.485000,10.554000,10.677000,10.828000 /)
XPIZA_LKT(99,8,1:6)=(/ 0.624308,0.782624,0.989080,0.990339,0.966524,0.548361 /)
XCGA_LKT(99,8,1:6)=(/ 0.914380,0.875350,0.840913,0.837843,0.838107,0.936510 /)
XEXT_COEFF_550_LKT(99,8)=10.471000 !rg=26.6054 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,9,1:6)=(/ 8.836200,8.862600,8.895300,8.955800,9.059800,9.158000 /)
XPIZA_LKT(99,9,1:6)=(/ 0.614445,0.766701,0.987793,0.989399,0.964309,0.549441 /)
XCGA_LKT(99,9,1:6)=(/ 0.917757,0.877760,0.842007,0.838923,0.838030,0.936663 /)
XEXT_COEFF_550_LKT(99,9)=8.872400 !rg=26.6054 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,10,1:6)=(/ 7.467200,7.484500,7.517500,7.556100,7.636700,7.716800 /)
XPIZA_LKT(99,10,1:6)=(/ 0.605194,0.749884,0.986155,0.988287,0.962199,0.550435 /)
XCGA_LKT(99,10,1:6)=(/ 0.921233,0.880613,0.842740,0.840043,0.840630,0.936797 /)
XEXT_COEFF_550_LKT(99,10)=7.499700 !rg=26.6054 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,11,1:6)=(/ 6.271600,6.288200,6.309200,6.336800,6.396100,6.463200 /)
XPIZA_LKT(99,11,1:6)=(/ 0.596552,0.732566,0.984250,0.986848,0.958731,0.551361 /)
XCGA_LKT(99,11,1:6)=(/ 0.000000,0.883860,0.843407,0.841227,0.841447,0.936910 /)
XEXT_COEFF_550_LKT(99,11)=6.292500 !rg=26.6054 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,12,1:6)=(/ 5.249300,5.260500,5.279300,5.302500,5.347100,5.394700 /)
XPIZA_LKT(99,12,1:6)=(/ 0.588720,0.714825,0.981887,0.985093,0.955109,0.552192 /)
XCGA_LKT(99,12,1:6)=(/ 0.000000,0.887227,0.844257,0.842013,0.844207,0.937003 /)
XEXT_COEFF_550_LKT(99,12)=5.266300 !rg=26.6054 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,13,1:6)=(/ 4.381400,4.390000,4.397200,4.419300,4.440800,4.491300 /)
XPIZA_LKT(99,13,1:6)=(/ 0.581687,0.697292,0.979054,0.982917,0.950375,0.552924 /)
XCGA_LKT(99,13,1:6)=(/ 0.000000,0.000000,0.844790,0.842817,0.845207,0.937073 /)
XEXT_COEFF_550_LKT(99,13)=4.390700 !rg=26.6054 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,14,1:6)=(/ 3.655500,3.661400,3.669300,3.681400,3.709700,3.738600 /)
XPIZA_LKT(99,14,1:6)=(/ 0.575495,0.680285,0.975756,0.980486,0.944783,0.553551 /)
XCGA_LKT(99,14,1:6)=(/ 0.000000,0.000000,0.845570,0.843490,0.846830,0.937123 /)
XEXT_COEFF_550_LKT(99,14)=3.662200 !rg=26.6054 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,15,1:6)=(/ 0.000000,3.050200,3.056300,3.071100,3.088700,3.108900 /)
XPIZA_LKT(99,15,1:6)=(/ 0.900,0.664211,0.971819,0.977548,0.938058,0.554086 /)
XCGA_LKT(99,15,1:6)=(/ 0.900,0.000000,0.000000,0.844607,0.848283,0.937160 /)
XEXT_COEFF_550_LKT(99,15)=3.053000 !rg=26.6054 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,16,1:6)=(/ 0.000000,2.535500,2.540200,2.547900,2.563400,2.579300 /)
XPIZA_LKT(99,16,1:6)=(/ 0.900,0.649043,0.967138,0.973966,0.930194,0.554546 /)
XCGA_LKT(99,16,1:6)=(/ 0.900,0.000000,0.000000,0.845497,0.849837,0.937180 /)
XEXT_COEFF_550_LKT(99,16)=2.536200 !rg=26.6054 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,17,1:6)=(/ 0.000000,0.000000,2.113600,2.118500,2.130700,2.143200 /)
XPIZA_LKT(99,17,1:6)=(/ 0.900,0.900,0.961696,0.969764,0.921233,0.554921 /)
XCGA_LKT(99,17,1:6)=(/ 0.900,0.900,0.000000,0.846073,0.851463,0.937187 /)
XEXT_COEFF_550_LKT(99,17)=2.111000 !rg=26.6054 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,18,1:6)=(/ 0.000000,0.000000,1.758900,1.764500,1.771200,1.781100 /)
XPIZA_LKT(99,18,1:6)=(/ 0.900,0.900,0.955357,0.964897,0.910963,0.555229 /)
XCGA_LKT(99,18,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.853117,0.937183 /)
XEXT_COEFF_550_LKT(99,18)=0.000000 !rg=26.6054 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,19,1:6)=(/ 0.000000,0.000000,0.000000,1.464800,1.471200,1.478100 /)
XPIZA_LKT(99,19,1:6)=(/ 0.900,0.900,0.900,0.959111,0.899265,0.555483 /)
XCGA_LKT(99,19,1:6)=(/ 0.900,0.900,0.900,0.000000,0.854907,0.937173 /)
XEXT_COEFF_550_LKT(99,19)=0.000000 !rg=26.6054 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,20,1:6)=(/ 0.000000,0.000000,0.000000,1.217300,1.221900,1.227400 /)
XPIZA_LKT(99,20,1:6)=(/ 0.900,0.900,0.900,0.952409,0.886113,0.555684 /)
XCGA_LKT(99,20,1:6)=(/ 0.900,0.900,0.900,0.000000,0.856750,0.937157 /)
XEXT_COEFF_550_LKT(99,20)=0.000000 !rg=26.6054 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,1,1:6)=(/ 20.907000,20.975000,20.774000,21.286000,21.678000,22.112000 /)
XPIZA_LKT(100,1,1:6)=(/ 0.670784,0.841000,0.992887,0.993151,0.969149,0.543292 /)
XCGA_LKT(100,1,1:6)=(/ 0.901090,0.865393,0.836740,0.834730,0.832437,0.935740 /)
XEXT_COEFF_550_LKT(100,1)=21.268000 !rg=28.8231 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,2,1:6)=(/ 20.091000,20.184000,20.218000,20.507000,20.808000,21.159000 /)
XPIZA_LKT(100,2,1:6)=(/ 0.668161,0.838205,0.992724,0.992989,0.968803,0.543669 /)
XCGA_LKT(100,2,1:6)=(/ 0.902077,0.866713,0.837950,0.833600,0.832410,0.935803 /)
XEXT_COEFF_550_LKT(100,2)=20.188000 !rg=28.8231 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,3,1:6)=(/ 18.633000,18.699000,18.858000,18.983000,19.366000,19.577000 /)
XPIZA_LKT(100,3,1:6)=(/ 0.662333,0.832393,0.992440,0.992743,0.968932,0.544310 /)
XCGA_LKT(100,3,1:6)=(/ 0.903617,0.867477,0.839427,0.834537,0.832737,0.935907 /)
XEXT_COEFF_550_LKT(100,3)=18.758000 !rg=28.8231 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,4,1:6)=(/ 16.846000,16.904000,17.014000,17.187000,17.380000,17.660000 /)
XPIZA_LKT(100,4,1:6)=(/ 0.655179,0.824359,0.991973,0.992313,0.969161,0.545118 /)
XCGA_LKT(100,4,1:6)=(/ 0.905367,0.868740,0.839207,0.835077,0.833923,0.936037 /)
XEXT_COEFF_550_LKT(100,4)=16.967000 !rg=28.8231 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,5,1:6)=(/ 14.923000,14.966000,15.061000,15.158000,15.352000,15.612000 /)
XPIZA_LKT(100,5,1:6)=(/ 0.646562,0.814102,0.991406,0.991924,0.968959,0.546045 /)
XCGA_LKT(100,5,1:6)=(/ 0.907543,0.870387,0.839830,0.835503,0.834047,0.936183 /)
XEXT_COEFF_550_LKT(100,5)=14.979000 !rg=28.8231 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,6,1:6)=(/ 13.038000,13.077000,13.155000,13.283000,13.417000,13.601000 /)
XPIZA_LKT(100,6,1:6)=(/ 0.637424,0.802618,0.990629,0.991408,0.968427,0.547039 /)
XCGA_LKT(100,6,1:6)=(/ 0.910243,0.872020,0.840880,0.836960,0.835670,0.936337 /)
XEXT_COEFF_550_LKT(100,6)=13.113000 !rg=28.8231 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,7,1:6)=(/ 11.237000,11.268000,11.321000,11.408000,11.540000,11.690000 /)
XPIZA_LKT(100,7,1:6)=(/ 0.627806,0.788995,0.989651,0.990820,0.967222,0.548086 /)
XCGA_LKT(100,7,1:6)=(/ 0.913163,0.874167,0.841097,0.838070,0.837483,0.936490 /)
XEXT_COEFF_550_LKT(100,7)=11.293000 !rg=28.8231 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,8,1:6)=(/ 9.610900,9.638600,9.691500,9.728600,9.858300,9.970000 /)
XPIZA_LKT(100,8,1:6)=(/ 0.618324,0.774241,0.988529,0.989929,0.965670,0.549111 /)
XCGA_LKT(100,8,1:6)=(/ 0.916297,0.876640,0.841913,0.838673,0.839543,0.936637 /)
XEXT_COEFF_550_LKT(100,8)=9.663500 !rg=28.8231 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,9,1:6)=(/ 8.152000,8.171300,8.201000,8.245200,8.319100,8.434200 /)
XPIZA_LKT(100,9,1:6)=(/ 0.608892,0.757955,0.987071,0.988867,0.963494,0.550105 /)
XCGA_LKT(100,9,1:6)=(/ 0.919663,0.879217,0.842373,0.839410,0.839623,0.936770 /)
XEXT_COEFF_550_LKT(100,9)=8.177500 !rg=28.8231 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,10,1:6)=(/ 6.889500,6.907600,6.927700,6.968200,7.021700,7.108300 /)
XPIZA_LKT(100,10,1:6)=(/ 0.600253,0.741331,0.985335,0.987609,0.960692,0.551019 /)
XCGA_LKT(100,10,1:6)=(/ 0.000000,0.882213,0.843090,0.840927,0.840980,0.936883 /)
XEXT_COEFF_550_LKT(100,10)=6.914900 !rg=28.8231 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,11,1:6)=(/ 5.787000,5.797200,5.821600,5.841800,5.902500,5.954500 /)
XPIZA_LKT(100,11,1:6)=(/ 0.592166,0.723695,0.983228,0.986065,0.957293,0.551867 /)
XCGA_LKT(100,11,1:6)=(/ 0.000000,0.885380,0.843903,0.841203,0.842993,0.936980 /)
XEXT_COEFF_550_LKT(100,11)=5.805300 !rg=28.8231 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,12,1:6)=(/ 4.843400,4.853400,4.862600,4.892400,4.919500,4.970900 /)
XPIZA_LKT(100,12,1:6)=(/ 0.584805,0.706241,0.980686,0.984189,0.953133,0.552627 /)
XCGA_LKT(100,12,1:6)=(/ 0.000000,0.000000,0.844523,0.842593,0.844443,0.937057 /)
XEXT_COEFF_550_LKT(100,12)=4.855800 !rg=28.8231 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,13,1:6)=(/ 4.042700,4.049800,4.060400,4.077200,4.103100,4.139300 /)
XPIZA_LKT(100,13,1:6)=(/ 0.578211,0.688948,0.977687,0.981945,0.947965,0.553295 /)
XCGA_LKT(100,13,1:6)=(/ 0.000000,0.000000,0.845273,0.843227,0.845950,0.937113 /)
XEXT_COEFF_550_LKT(100,13)=4.050700 !rg=28.8231 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,14,1:6)=(/ 0.000000,3.377700,3.384700,3.399300,3.419500,3.446000 /)
XPIZA_LKT(100,14,1:6)=(/ 0.900,0.672384,0.974086,0.979243,0.941845,0.553864 /)
XCGA_LKT(100,14,1:6)=(/ 0.900,0.000000,0.845977,0.844197,0.847347,0.937153 /)
XEXT_COEFF_550_LKT(100,14)=3.379900 !rg=28.8231 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,15,1:6)=(/ 0.000000,2.815400,2.819800,2.829400,2.843400,2.866000 /)
XPIZA_LKT(100,15,1:6)=(/ 0.900,0.656921,0.969859,0.976039,0.934714,0.554349 /)
XCGA_LKT(100,15,1:6)=(/ 0.900,0.000000,0.000000,0.844973,0.848857,0.937177 /)
XEXT_COEFF_550_LKT(100,15)=2.816200 !rg=28.8231 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,16,1:6)=(/ 0.000000,0.000000,2.343600,2.349000,2.364800,2.378000 /)
XPIZA_LKT(100,16,1:6)=(/ 0.900,0.900,0.964847,0.972180,0.926502,0.554763 /)
XCGA_LKT(100,16,1:6)=(/ 0.900,0.900,0.000000,0.845590,0.850817,0.937190 /)
XEXT_COEFF_550_LKT(100,16)=2.340400 !rg=28.8231 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,17,1:6)=(/ 0.000000,0.000000,1.949600,1.955700,1.964000,1.976200 /)
XPIZA_LKT(100,17,1:6)=(/ 0.900,0.900,0.959012,0.967716,0.916828,0.555100 /)
XCGA_LKT(100,17,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.852293,0.937193 /)
XEXT_COEFF_550_LKT(100,17)=0.000000 !rg=28.8231 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,18,1:6)=(/ 0.000000,0.000000,1.622800,1.626800,1.632700,1.642500 /)
XPIZA_LKT(100,18,1:6)=(/ 0.900,0.900,0.952246,0.962491,0.905939,0.555375 /)
XCGA_LKT(100,18,1:6)=(/ 0.900,0.900,0.000000,0.000000,0.853743,0.937183 /)
XEXT_COEFF_550_LKT(100,18)=0.000000 !rg=28.8231 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,19,1:6)=(/ 0.000000,0.000000,0.000000,1.351200,1.357600,1.363200 /)
XPIZA_LKT(100,19,1:6)=(/ 0.900,0.900,0.900,0.956304,0.893704,0.555599 /)
XCGA_LKT(100,19,1:6)=(/ 0.900,0.900,0.900,0.000000,0.855813,0.937170 /)
XEXT_COEFF_550_LKT(100,19)=0.000000 !rg=28.8231 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,20,1:6)=(/ 0.000000,0.000000,0.000000,1.123200,1.126800,1.132100 /)
XPIZA_LKT(100,20,1:6)=(/ 0.900,0.900,0.900,0.949148,0.879877,0.555776 /)
XCGA_LKT(100,20,1:6)=(/ 0.900,0.900,0.900,0.000000,0.000000,0.937150 /)
XEXT_COEFF_550_LKT(100,20)=0.000000 !rg=28.8231 sigma=2.95 wvl=0.55
 

END SUBROUTINE SALT_OPT_LKT_SET10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MODE_SALTOPT
