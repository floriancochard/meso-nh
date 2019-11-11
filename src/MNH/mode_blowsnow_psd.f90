!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!   ########################
     MODULE MODE_BLOWSNOW_PSD
!!   ########################
!!
!!    PURPOSE
!!    -------
!! MODULE BLOWSNOW PSD (Particle Size Distribution)
!! Purpose: Contains subroutines to convert from transported variables (#/kg_{air},kg_{snow}/kg_{air})
!! to understandable BLOWSNOW variables, e.g. #/m3, kg/m3, sigma, R_{n}
!!
!!    AUTHOR
!!    ------
!!      Vincent Vionnet (GMME) based on Alf Grini (CNRM/GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODD_CSTS_BLOWSNOW         !Constants which are important for drifting snow calculations
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
!
IMPLICIT NONE
!
CONTAINS
!
!!   ############################################################
  SUBROUTINE PPP2SNOW(             &
       PSVT                         & !I [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PBET3D                     & !O [m] scale parameter of snow distribution
       , PRG3D                      & !O [um] mean radius of snow  distribution
       , PN3D                       & !O [#/m3] number concentration of snow particles
       , PMASS3D                    & !O [kg/m3] mass concentration of snow particles
       , PMOB3D                     & !O [-] mobility index
       ,PM3D                        & !O   BLOWSNOW moments
         )
!!   ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the two moments M0 and, M3 
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Vincent Vionnet based on routine from Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!    None
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)  :: PSVT      !I [ppp] first moment
REAL,       DIMENSION(:,:,:),    INTENT(IN)     :: PRHODREF !I [kg/m3] density of air

REAL,       DIMENSION(:,:,:),  OPTIONAL, INTENT(OUT)     :: PBET3D   !O [-] scale parameter
REAL,       DIMENSION(:,:,:),  OPTIONAL, INTENT(OUT)     :: PRG3D    !O [um] mean radius
REAL,       DIMENSION(:,:,:),  OPTIONAL, INTENT(OUT)     :: PN3D     !O [#/m3] number concentration
REAL,       DIMENSION(:,:,:),  OPTIONAL, INTENT(OUT)     :: PMASS3D  !O [kg_{snw}/m3] mass concentration
REAL,       DIMENSION(:,:,:),  OPTIONAL, INTENT(OUT)     :: PMOB3D   !O [-] mobility index
REAL,       DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)   :: PM3D     !O   BLOWSNOW moments
!
!
!*      0.2    declarations local variables
!
  REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSV                ! [snow particles moment concentration]
  REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM                 ! [aerosol units] local array which goes to output later
  REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZBETA              ! [-] standard deviation
  REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZRG                ! [m] number median diameter
!  REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZMOB               ! [-] mobility index
  REAL,DIMENSION(2)                     :: ZMMIN              ! [BLOWSNOW units] minimum values for N and M
  REAL                                  :: ZRGMIN             ! [m] minimum radius accepted
!
!-------------------------------------------------------------------------------
!
!        1.1    initialisation
!
ALLOCATE (ZBETA(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)))
ALLOCATE (ZRG(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)))
!ALLOCATE (ZMOB(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)))
ALLOCATE (ZSV(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), SIZE(PSVT,4)))
ALLOCATE (ZM(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), SIZE(PSVT,4)))

ZSV(:,:,:,:) = MAX(PSVT(:,:,:,:), 1.E-80)
!Get minimum values possible
  ZMMIN(1) = XN0MIN_SNW
  ZRGMIN   = XINIRADIUS_SNW
  ZMMIN(2) = 4*XPI*XRHOLI/3*(ZRGMIN/XALPHA_SNOW)**3.*(XALPHA_SNOW+2)*(XALPHA_SNOW+1)*XALPHA_SNOW*XN0MIN_SNW

!Get number concentration (#/m3) from M0 (#/kg_{air})
ZM(:,:,:,1)=ZSV(:,:,:,1)*PRHODREF(:,:,:)
!Get mass concentration (kg_{snow}/m3) from M3 (kg_{snow})/kg_{air})
ZM(:,:,:,2)=ZSV(:,:,:,2)*PRHODREF(:,:,:)
!Get mobility index in term of mass concentration (kg_{snow}/m3) from M3 (kg_{snow})/kg_{air})
!ZM(:,:,:,3)=ZSV(:,:,:,3)*PRHODREF(:,:,:)


! Limit concentration to minimum values
    WHERE ((ZM(:,:,:,1) < ZMMIN(1) ).OR. &
           (ZM(:,:,:,2) < ZMMIN(2)))
       ZM(:,:,:,1) = ZMMIN(1)
       ZM(:,:,:,2) = ZMMIN(2)
!       ZM(:,:,:,3) = 1.2*ZMMIN(2)
       PSVT(:,:,:,1) = ZM(:,:,:,1)/ PRHODREF(:,:,:) 
       PSVT(:,:,:,2) = ZM(:,:,:,2)/ PRHODREF(:,:,:)  
!       PSVT(:,:,:,3) = ZM(:,:,:,3)/ PRHODREF(:,:,:) 
    ENDWHERE

ZBETA(:,:,:)=((3*ZM(:,:,:,2))/(4*XPI*XRHOLI*(XALPHA_SNOW+2)*(XALPHA_SNOW+1)*XALPHA_SNOW*ZM(:,:,:,1)))**(1./3.)

ZRG(:,:,:)  = ZBETA(:,:,:)*XALPHA_SNOW

!ZMOB(:,:,:) = ZM(:,:,:,3)/ZM(:,:,:,2)

!Give the beta-values to the passed array
IF(PRESENT(PBET3D))THEN
        PBET3D(:,:,:) = ZBETA(:,:,:)
ENDIF

!Get the mean radius
IF(PRESENT(PRG3D))THEN
     PRG3D(:,:,:)= ZRG(:,:,:)
ENDIF

!Get the mobility index
!IF(PRESENT(PMOB3D))THEN
!     PMOB3D(:,:,:)= ZMOB(:,:,:)
!ENDIF

!Get the number concentration
IF(PRESENT(PN3D))THEN
     PN3D(:,:,:)= ZM(:,:,:,1)
ENDIF

!Get the mass concentration
IF(PRESENT(PMASS3D))THEN
     PMASS3D(:,:,:)= ZM(:,:,:,2)
ENDIF

!Get number, mass concentration and mobility index
IF(PRESENT(PM3D))THEN
     PM3D(:,:,:,:)= ZM(:,:,:,:)
ENDIF


    DEALLOCATE(ZSV)
    DEALLOCATE(ZRG)
    DEALLOCATE(ZBETA)
    DEALLOCATE(ZM)
!    DEALLOCATE(ZMOB)
!
!
END SUBROUTINE PPP2SNOW


END MODULE MODE_BLOWSNOW_PSD
