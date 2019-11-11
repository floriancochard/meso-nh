!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_POLAR_MEAN
!     ######################
INTERFACE POLAR_MEAN
!
      SUBROUTINE POLAR_MEAN_P(PVARIN,PR0,PICEN,PJCEN,PHDR0MOY)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARIN !3-D array of input field
REAL, DIMENSION(:,:),  INTENT(IN)  :: PR0  ! 2-D array of the radius of the 
                                           ! filtering domain based on the
                                           ! tangential wind criterion
INTEGER, DIMENSION(:), INTENT(IN)  :: PICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)  :: PJCEN !at level pk
REAL, DIMENSION(:), INTENT(OUT)    :: PHDR0MOY ! 1-D array of the average along 
                                               ! the periphery of the filter
                                               ! at the radius ro(p)
!
END SUBROUTINE POLAR_MEAN_P
!
      SUBROUTINE POLAR_MEAN_RP(PVARIN,PR0,PICEN,PJCEN,PHDR0MOY)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARIN !3D array of input field
REAL, DIMENSION(:,:,:), INTENT(IN) :: PR0    !3D array of the radius ro(r,phi,p)
INTEGER, DIMENSION(:), INTENT(IN)  :: PICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)  :: PJCEN !at level pk
REAL, DIMENSION(:,:),  INTENT(OUT) :: PHDR0MOY ! 2-D array of the average along 
                                               !peripheries of ro(r,phi,p)
!
END SUBROUTINE POLAR_MEAN_RP
!
END INTERFACE POLAR_MEAN
!
END MODULE MODI_POLAR_MEAN
!
!
!
!     ######################################################
      SUBROUTINE POLAR_MEAN_P(PVARIN,PR0,PICEN,PJCEN,PHDR0MOY)
!     ######################################################
!
!!****  *POLAR_MEAN* - This subroutine calculates Hdbar(r0) following
!!                          the Kurihara et al. (1993) formulation 
!!
!!    PURPOSE
!!    -------
!         This subroutine calculates the average of field
!           varin around the center along the periphery of the
!          filter r0(p)
!
!!**  METHOD
!!    ------
!!      
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      R. Rogers             * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY: XPI
USE MODD_GRID_n, ONLY: XXHAT,XYHAT
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARIN !3-D array of input field to convert
REAL, DIMENSION(:,:),  INTENT(IN)  :: PR0  ! 2-D array of the radius of the 
                                           ! filtering domain based on the
                                           ! tangential wind criterion
INTEGER, DIMENSION(:), INTENT(IN)  :: PICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)  :: PJCEN ! at level pk
REAL, DIMENSION(:), INTENT(OUT)    :: PHDR0MOY ! 1-D array of the average along 
                                               ! the periphery of the filter
                                               ! at the radius ro(p)
!
!*       0.2   Declarations of local variables
!
REAL                                          :: ZPHI,ZDPHI,ZDELTAX,ZDELTAY
REAL                                          :: ZXI0,ZYJ0,ZX00,ZY00
REAL                                          :: ZXK,ZYK     
REAL                                          :: ZHDR0
INTEGER                                       :: IIX,IJY 
INTEGER                                       :: IP,IPHI
INTEGER                                       :: ICOUNT
INTEGER                                       :: IIU,IJU,IIB,IJB,IIE,IJE
INTEGER                                       :: JPHI,JP
!
!-------------------------------------------------------------------------------
!
!*	 1. INITIALIZATIONS
!           ---------------
!
!
IP= SIZE(PVARIN,3)
IPHI= SIZE(PR0,1)
!
IIU=SIZE(PVARIN,1)
IJU=SIZE(PVARIN,2)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
ZDELTAX = XXHAT(3) - XXHAT(2)
ZDELTAY = XYHAT(3) - XYHAT(2)
ZDPHI = 2.*XPI /IPHI
!
PHDR0MOY(:) = 0.
!
!-------------------------------------------------------------------------------
!
!*	 2. 
!
!
DO JP = 1, IP
  PHDR0MOY(JP) = 0.
  ZXI0 = XXHAT(PICEN(JP)) + (ZDELTAX / 2.)
  ZYJ0 = XYHAT(PJCEN(JP)) + (ZDELTAY / 2.)
  ZX00 = XXHAT(1) + (ZDELTAX / 2.) 
  ZY00 = XYHAT(1) + (ZDELTAY / 2.)
  ICOUNT = 0
  DO JPHI = 1, IPHI
    ZPHI = (JPHI - 1) * ZDPHI
    ZXK = PR0(JPHI,JP) * COS(ZPHI) + ZXI0
    ZYK = PR0(JPHI,JP) * SIN(ZPHI) + ZYJ0
    IIX = (ZXK - ZX00) / ZDELTAX + 1
    IJY = (ZYK - ZY00) / ZDELTAY + 1
    !
    IF (IIX >= IIB .AND. (IIX+1) <= IIE      &
        .AND. IJY >= IJB .AND. (IJY+1) <= IJE) THEN 
      ZXK = (ZXK - XXHAT(IIX)) / ZDELTAX - 0.5
      ZYK = (ZYK - XYHAT(IJY)) / ZDELTAY - 0.5
      ZHDR0 = PVARIN(IIX,IJY,JP)*(1-ZXK)*(1-ZYK) &
             +PVARIN(IIX+1,IJY,JP)*ZXK*(1-ZYK)   &
             +PVARIN(IIX,IJY+1,JP)*(1-ZXK)*ZYK   &
             +PVARIN(IIX+1,IJY+1,JP)*ZXK*ZYK   
      !
      PHDR0MOY(JP) = PHDR0MOY(JP) + ZHDR0
      ICOUNT = ICOUNT+1
    END IF
  END DO
  IF (ICOUNT /= 0) PHDR0MOY(JP) = PHDR0MOY(JP) / ICOUNT
END DO
!
END SUBROUTINE POLAR_MEAN_P
!-------------------------------------------------------------------------------
!
!     #########################################################
      SUBROUTINE POLAR_MEAN_RP(PVARIN,PR0,PICEN,PJCEN,PHDR0MOY)
!     #########################################################
!
!!****  *POLAR_MEAN3D* - This subroutine interpolates into cyclindrical 
!!     coordinates (r,phi,k)  and does azimuthal average of field
!!
!!    PURPOSE
!!    -------
!         This subroutine calculates the average of field
!           varin around the center along the periphery of the
!          filter r0(p)
!
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      R. Rogers             * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY: XPI
USE MODD_GRID_n, ONLY: XXHAT,XYHAT
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARIN !3D array of input field
REAL, DIMENSION(:,:,:), INTENT(IN) :: PR0    !3D array of the radius ro(r,phi,p)
INTEGER, DIMENSION(:), INTENT(IN)  :: PICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)  :: PJCEN !at level pk
REAL, DIMENSION(:,:),  INTENT(OUT) :: PHDR0MOY ! 2-D array of the average along 
                                               !peripheries of ro(r,phi,p)
!
!*       0.2   Declarations of local variables
!
REAL                                          :: ZPHI,ZDPHI,ZDELTAX,ZDELTAY
REAL                                          :: ZXI0,ZYJ0,ZX00,ZY00
REAL                                          :: ZXK,ZYK     
INTEGER                                       :: IIX,IJY
INTEGER                                       :: JR,JPHI,JP
INTEGER                                       :: INR,IPHI,IP
INTEGER                                       :: ICOUNT
INTEGER                                       :: IIU,IJU,IIB,IJB,IIE,IJE
REAL                                          :: ZHDR0     
!
!-------------------------------------------------------------------------------
!
!*	 1. INITIALIZATIONS
!           ---------------
!
!
IP= SIZE(PVARIN,3)
INR = SIZE(PR0,1)
IPHI= SIZE(PR0,2)
!
IIU=SIZE(PVARIN,1)
IJU=SIZE(PVARIN,2)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
ZDELTAX = XXHAT(3) - XXHAT(2)
ZDELTAY = XYHAT(3) - XYHAT(2)
ZDPHI = 2.*XPI /IPHI
!
PHDR0MOY(:,:) = 0.
!
!-------------------------------------------------------------------------------
!
!*	 2. 
!
DO JP = 1, IP
  ZXI0 = XXHAT(PICEN(JP)) + (ZDELTAX / 2.)
  ZYJ0 = XYHAT(PJCEN(JP)) + (ZDELTAY / 2.)
  ZX00 = XXHAT(1) + (ZDELTAX / 2.) 
  ZY00 = XYHAT(1) + (ZDELTAY / 2.)
  DO JR = 1, INR
    ICOUNT = 0
    DO JPHI = 1, IPHI
      ZPHI = (JPHI - 1) * ZDPHI
      ZXK = PR0(JR,JPHI,JP) * COS(ZPHI) + ZXI0
      ZYK = PR0(JR,JPHI,JP) * SIN(ZPHI) + ZYJ0
      IIX = (ZXK - ZX00) / ZDELTAX + 1
      IJY = (ZYK - ZY00) / ZDELTAY + 1
      !
      IF (IIX >= IIB .AND. (IIX+1) <= IIE      &
          .AND. IJY >= IJB .AND. (IJY+1) <= IJE) THEN 
        ZXK = (ZXK - XXHAT(IIX)) / ZDELTAX - 0.5
        ZYK = (ZYK - XYHAT(IJY)) / ZDELTAY - 0.5
        ZHDR0 = PVARIN(IIX,IJY,JP)*(1-ZXK)*(1-ZYK) &
               +PVARIN(IIX+1,IJY,JP)*ZXK*(1-ZYK)   &
               +PVARIN(IIX,IJY+1,JP)*(1-ZXK)*ZYK   &
               +PVARIN(IIX+1,IJY+1,JP)*ZXK*ZYK   
        PHDR0MOY(JR,JP) = PHDR0MOY(JR,JP) + ZHDR0
        ICOUNT = ICOUNT + 1
      END IF
      !
    END DO
    IF (ICOUNT /= 0) PHDR0MOY(JR,JP) = PHDR0MOY(JR,JP) / ICOUNT    
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE POLAR_MEAN_RP
